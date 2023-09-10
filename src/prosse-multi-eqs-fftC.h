/* 
   In this more simple version of the integrator, there is just one
   structure: 
     "prosse_multi_fft"
   This holds everything needed to perform an integration.  It will
   probably be the case that we will want to make many of these to do
   things efficiently, but the C code is agnostic about this, and will
   just provide a few functions to get data in and out of these, as
   well as perform the integrations.

   This version allows there to be multiple nested plans.
*/
typedef struct {
  int n_fft; /* number of different plans */

  int nx;    /* x-extent */
  double dx; /* distance between x's */

  int *nd;    /* number of dimensions for each plan */

  /* Data */
  double *x;
  fftw_complex *y;

  /* Diversification information */
  double *lambda;
  double *Q;
  double *Qr;

  /* Drift and diffusion parameters */
  double drift;
  double diffusion;

  /* Scratch space */
  double *z;

  double *wrk;
  double *wrkd;

  /* Transform for x propagation */
  rfftw_plan_real **fft;

  /* Kernel information (space propagation) */
  double ny;      /* Fourier space extent */
  fftw_complex *fkern; /* Gaussian kernel transformation */
  int nkl;        /* Kernel width to the left */
  int nkr;        /* Kernel width to the right */
  int npad;       /* nkl + nkr + 1 */
  int ndat;       /* nx - npad */

  /* The kernel itself */
  double  *kern_x;
  fftw_complex *kern_y;
  rfftw_plan_real *kernel;
} prosse_multi_fft;

prosse_multi_fft* make_prosse_multi_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags);
void qf_copy_x_prosse_multi(prosse_multi_fft *obj, double *x, int nd, int copy_in);
void qf_setup_kern_prosse_multi(prosse_multi_fft *obj, double drift, double diffusion,
		   double dt, int nkl, int nkr);
void do_integrate_prosse_multi(prosse_multi_fft *obj, int star, double Et, double b, double mu, int nt, double dt, int idx);
void propagate_t_prosse_multi(prosse_multi_fft *obj, int i, int star, double Et, double b, double mu, double dt, int idx);
void propagate_x_prosse_multi(prosse_multi_fft *obj, int idx);
