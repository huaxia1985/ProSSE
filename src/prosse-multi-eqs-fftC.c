#include "config.h"

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> /* why? */
#include <Rmath.h> /* for dnorm() */

#ifdef HAVE_FFTW3_H
#include <complex.h>
#include <fftw3.h>
#include "rfftw.h"
#include "quasse-eqs-fftC.h"
#include "prosse-multi-eqs-fftC.h"

static void prosse_multi_fft_finalize(SEXP extPtr);

SEXP r_make_prosse_multi_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  prosse_multi_fft *obj;
  SEXP extPtr;
  int nx = INTEGER(r_nx)[0];
  double dx = REAL(r_dx)[0];
  int n_fft = LENGTH(r_nd);
  int i;
  int flags;
  int *nd = (int*)calloc(n_fft, sizeof(int));
  for ( i = 0; i < n_fft; i++ )
    nd[i] = INTEGER(r_nd)[i];
  
  if ( INTEGER(r_flags)[0] == -1 )
    flags = FFTW_ESTIMATE;
  else if ( INTEGER(r_flags)[0] == 1 )
    flags = FFTW_PATIENT;
  else if ( INTEGER(r_flags)[0] == 2 )
    flags = FFTW_EXHAUSTIVE;
  else
    flags = FFTW_MEASURE;

  obj = make_prosse_multi_fft(n_fft, nx, dx, nd, flags);

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, prosse_multi_fft_finalize);

  return extPtr;
}

SEXP r_do_integrate_prosse_multi(SEXP extPtr, SEXP vars, SEXP Et, SEXP star, SEXP b, SEXP mu, SEXP lambda, 
		    SEXP Q, SEXP Qr, SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  prosse_multi_fft *obj = (prosse_multi_fft*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int ndat = LENGTH(lambda);
  int c_star = INTEGER(star)[0], c_nt = INTEGER(nt)[0];
  double c_dt = REAL(dt)[0];
  double *c_lambda=REAL(lambda);
    double c_Et=REAL(Et)[0], c_b=REAL(b)[0], c_mu=REAL(mu)[0], c_drift=REAL(drift)[0], c_diffusion=REAL(diffusion)[0];
    int i, idx, nd;
  if ( obj == NULL )
    error("Corrupt ProSSE integrator: ptr is NULL (are you using multicore?)");
  nd = LENGTH(vars) / obj->nx;
  
  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x_prosse_multi(obj, REAL(vars), nd, 1);

  obj->lambda = REAL(lambda);
  obj->Q = REAL(Q);
  obj->Qr = REAL(Qr);

  for ( i = 0; i < ndat; i++ )
    obj->z[i] = exp(- c_dt * c_lambda[i]);

  qf_setup_kern_prosse_multi(obj, c_drift, c_diffusion, c_dt, nkl, nkr);

  do_integrate_prosse_multi(obj, c_star, c_Et, c_b, c_mu, c_nt, c_dt, idx);

  obj->lambda = NULL;
  obj->Q = NULL;
  obj->Qr = NULL;
    
  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x_prosse_multi(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  
  return ret;
}

/* This does the memory allocation and plans the FFT transforms */
prosse_multi_fft* make_prosse_multi_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags) {
  prosse_multi_fft *obj = calloc(1, sizeof(prosse_multi_fft));
  int ny = (((int)floor(nx/2)) + 1);
  int i, max_nd=1;
  for ( i = 0; i < n_fft; i++ )
    if ( nd[i] > max_nd )
      max_nd = nd[i];

  obj->n_fft = n_fft;
  obj->nx = nx;
  obj->ny = ny;
  obj->dx = dx;
  obj->nd = nd;

  obj->x   = fftw_malloc(max_nd *  nx    * sizeof(double));
  obj->y   = fftw_malloc(max_nd * (ny+1) * sizeof(fftw_complex));
  
  obj->z   = (double*)calloc(nx, sizeof(double));
  obj->wrk = (double*)calloc(nx, sizeof(double));
  obj->wrkd = (double*)calloc(max_nd *  nx, sizeof(double));

  obj->fft = (rfftw_plan_real**)calloc(n_fft, sizeof(rfftw_plan_real*));

  for ( i = 0; i < n_fft; i++ ) {
    obj->fft[i] = make_rfftw_plan_real(nd[i], nx, DIR_COLS,
				       obj->x, obj->y, flags);
  }
  
  /* Brownian kernel */
  obj->kern_x = fftw_malloc(nx     * sizeof(double));
  obj->kern_y = fftw_malloc((ny+1) * sizeof(fftw_complex));
  obj->kernel = make_rfftw_plan_real(1, nx, DIR_COLS, 
				     obj->kern_x, obj->kern_y, flags);
  
  return obj;
}

static void prosse_multi_fft_finalize(SEXP extPtr) {
  prosse_multi_fft *obj = (prosse_multi_fft*)R_ExternalPtrAddr(extPtr);
  int i;
  /* Rprintf("Cleaning up\n"); */

  for ( i = 0; i < obj->n_fft; i++ ) {
    fftw_destroy_plan(obj->fft[i]->plan_f);
    fftw_destroy_plan(obj->fft[i]->plan_b);
  }
  free(obj->fft);
  free(obj->nd);

  fftw_free(obj->x);
  fftw_free(obj->y);
  
  free(obj->z);
  free(obj->wrk);
  free(obj->wrkd);

  fftw_destroy_plan(obj->kernel->plan_f);
  fftw_destroy_plan(obj->kernel->plan_b);

  fftw_free(obj->kern_x);
  fftw_free(obj->kern_y);

  free(obj);
}


void qf_copy_x_prosse_multi(prosse_multi_fft *obj, double *x, int nd, int copy_in) {
  int i, n = obj->nx * nd;
  double *fft_x = obj->x;
  if ( copy_in )
    for ( i = 0; i < n; i++ )
      fft_x[i] = x[i];
  else {
    for ( i = 0; i < n; i++ ) {
      x[i] = fft_x[i];
    }
  }
}

void qf_setup_kern_prosse_multi(prosse_multi_fft *obj, double drift, double diffusion,
		   double dt, int nkl, int nkr) {
  const int nx = obj->nx;  
  int i;
  double x, *kern_x=obj->kern_x, tot=0, dx=obj->dx;
  double mean, sd;
  
  obj->nkl  = nkl;
  obj->nkr  = nkr;
  obj->npad = nkl + 1 + nkr;
  obj->ndat = nx - obj->npad;
  obj->drift = drift;
  obj->diffusion = diffusion;

  tot   = 0;
  mean  = - dt * drift;
  sd = sqrt(dt * diffusion);

  for ( i = 0, x = 0; i <= nkr; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;
  for ( i = nkr + 1; i < nx - nkl; i++ )
    kern_x[i] = 0;
  for ( i = nx - nkl, x = -nkl*dx; i < nx; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;

  for ( i = 0;        i <= nkr; i++ ) kern_x[i] /= tot;
  for ( i = nx - nkl; i < nx;   i++ ) kern_x[i] /= tot;

  fftw_execute(obj->kernel->plan_f);
}

void do_integrate_prosse_multi(prosse_multi_fft *obj, int star, double Et, double b, double mu, int nt, double dt, int idx) {
  int i, nkl=obj->nkl;
  for ( i = 0; i < nt; i++ ) {
    propagate_t_prosse_multi(obj, i, star, Et, b, mu, dt, idx);
    propagate_x_prosse_multi(obj, idx);
    if ( ISNAN(obj->x[nkl]) )
      error("Integration failure at step %d\n", i);
  }
}

/* Lower level functions */
void propagate_t_prosse_multi(prosse_multi_fft *obj, int i, int star, double Et, double b, double mu, double dt, int idx) {
    int ix, id, ik, nx=obj->nx, ndat=obj->ndat, nd=obj->nd[idx], nd2=nd/2;
    double *d = obj->wrkd, *dd = obj->wrk, *vars=obj->x;
  double tmp1, tmp2, z_x, Q_x, d_x, z1, e;
  
  z1 = exp(i * dt * (b - mu));
  e = (mu + z1*(Et - 1)*mu - b*Et) / (mu + z1*(Et - 1)*b - b*Et);
  
  z1 = exp(dt * (b - mu));
  tmp1 = z1 * b - mu + b * e * (1 - z1);
  tmp2 = z1 * (b - mu) * (b - mu) / tmp1 / tmp1;
  	
  for ( ix = 0; ix < ndat; ix++ ) {
    z_x = obj->z[ix];
    /* Here is the D scaling factor */
    dd[ix] = tmp2 * z_x;
  }

  /* Update the D values */
  if (star == 1) { 
  
  for ( id = 0; id < nd; id++ ) {
    for ( ix = 0; ix < ndat; ix++ ) {
        if (vars[ nx * id + ix ] < 0) {
            vars[ nx * id + ix ] = 0;
        } else {
            vars[ nx * id + ix ] *= dd[ix];
        }
    }
  }
      
  for ( id = 0; id < nd2; id++ ) {
  	for ( ix = 0; ix < ndat; ix++ ) {
  		d[ nx * id + ix ] = 0;
  		for (ik = 0; ik < nd2; ik++ ) {
            Q_x = obj->Q[ ik * nd2 + id ];
            d_x = vars[ nx * ik + ix ];
  			d[ nx * id + ix ] += Q_x * d_x;
  		}
  	}
  }
      
  for ( id = 0; id < nd2; id++ ) {
      for ( ix = 0; ix < ndat; ix++ ) {
          d[ nx * (id+nd2) + ix ] = 0;
          for (ik = 0; ik < nd2; ik++ ) {
            Q_x = obj->Q[ ik * nd2 + id ];
            d_x = vars[ nx * (ik+nd2) + ix ];
            d[ nx * (id+nd2) + ix ] += Q_x * d_x;
          }
      }
  }

  for ( id = 0; id < nd; id++ ) {
      for ( ix = 0; ix < ndat; ix++ ) {
          vars[ nx * id + ix ] = d[ nx * id + ix ];
      }
  }
  
  } else if (star == 2) {
  	
  /* Update the DR values */
      for ( id = nd2; id < nd; id++ ) {
        for ( ix = 0; ix < ndat; ix++ ) {
            if (vars[ nx * id + ix ] < 0) {
                vars[ nx * id + ix ] = 0;
            } else {
                vars[ nx * id + ix ] *= tmp2;
            }
        }
      }
      
  for ( id = 0; id < nd2; id++ ) {
  	for ( ix = 0; ix < ndat; ix++ ) {
  		d[ nx * (id+nd2) + ix ] = 0;
  		for (ik = 0; ik < nd2; ik++ ) {
  			Q_x = obj->Qr[ ik * nd2 + id ];
  			d_x = vars[ nx * (ik + nd2) + ix ];
  			d[ nx * (id+nd2) + ix ] += d_x * Q_x;
  		}
  	}
  }
  
  /* Update the DI values */
  for ( id = 0; id < nd2; id++ ) {
    for ( ix = 0; ix < ndat; ix++ ) {
        if (vars[ nx * id + ix ] < 0) {
            vars[ nx * id + ix ] = 0;
        } else {
            vars[ nx * id + ix ] *= dd[ix];
        }
	    d_x = vars[ nx * (id + nd2) + ix ];
	    z_x = obj->z[ix];
	    vars[ nx * id + ix ] -= d_x * z_x;
    }
  }
  
  for ( id = 0; id < nd2; id++ ) {
  	for ( ix = 0; ix < ndat; ix++ ) {
  		d[ nx * id + ix ] = 0;
  		for (ik = 0; ik < nd2; ik++ ) {
  			Q_x = obj->Q[ ik * nd2 + id ];
  			d_x = vars[ nx * ik + ix ];
  			d[ nx * id + ix ] += d_x * Q_x;
  		}
  	}
  }

  /* copy the DR values */
  for ( id = nd2; id < nd; id++ ) {
      for ( ix = 0; ix < ndat; ix++ ) {
          vars[ nx * id + ix ] = d[ nx * id + ix ];
      }
  }
  /* copy the DI values */
  for ( id = 0; id < nd2; id++ ) {
      for ( ix = 0; ix < ndat; ix++ ) {
      	d_x = d[ nx * id + ix ];
        vars[ nx * id + ix ] = vars[ nx * (id + nd2) + ix ] + d_x;
      }
  }
  }
      
}
  
void propagate_x_prosse_multi(prosse_multi_fft *obj, int idx) {
  double *x = obj->x, *d = obj->wrkd;
  int i, id, nx = obj->nx, nd=obj->nd[idx];
  int nkl = obj->nkl, nkr = obj->nkr, npad = obj->npad;

  /* TODO: I am not sure if these are flipped nkl/nkr */
    
  for ( id = 0; id < nd; id++ ) {
  for ( i = 0; i < nkl; i++ )
    d[i + nx * id] = x[i + nx * id];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    d[i + nx * id] = x[i + nx * id];
  }

  convolve(obj->fft[idx], obj->kern_y);

  for ( id = 0; id < nd; id++ ) {
  for ( i = 0; i < nkl; i++ )
    x[i + nx * id] = d[i + nx * id];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    x[i + nx * id] = d[i + nx * id];
  }

    /*for ( id = 0; id < (nx*nd); id++ ) {
      if ( x[id] < 0 )
        x[id] = 0;
    }
    
  /* Zeroing takes a little more work, now.  We might be able to get
     away with just zeroing the E though */
  for ( id = 0; id < nd; id++ ) {
    x = obj->x + (obj->nx)*(id+1) - npad;
    for ( i = 0; i < npad; i++ ) 
      x[i] = 0;
  }
    
}

#else
SEXP r_do_integrate_prosse_multi(SEXP extPtr, SEXP vars, SEXP star, SEXP b, SEXP mu, SEXP lambda,  SEXP Q, SEXP Qr,
		    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  Rf_error("FFTW support not included");
}
SEXP r_make_prosse_multi_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  Rf_error("FFTW support not included");
}
#endif
