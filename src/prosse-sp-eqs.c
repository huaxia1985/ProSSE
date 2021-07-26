/*
 * These are the protracted_representative equations, implemented in c
 */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_prosse_sp(double *pars, const double *y, double *ydot) {
  double Ei = y[0], Er = y[1], D = y[2], Dr = y[3], Di = y[4], Dir = y[5];
  double b = pars[0], mu = pars[1], lambda = pars[2];

  ydot[0] = -(b+mu) * Ei + b * Ei * Ei + mu;
  ydot[1] = -(b+mu+lambda) * Er + b * Er * Er + mu + lambda * Ei;
  ydot[2] = -(b+mu+lambda) * D + b * D * Er;
  ydot[3] = -(b+mu) * Dr + 2 * b * Dr * Ei + lambda * D;
  ydot[4] = -(b+mu+lambda) * Di + 2* b * Di * Er + lambda * (D+Dr);
  ydot[5] = -(b+mu+lambda) * Dir + 2 * b * Dir * Ei + lambda * (D+Dr);
}

/* 
   Wrap these up for GslOde, which allows time and neq to be here too.
   
   TODO: Eventually move the contents of do_derivs_* into these
   functions, as no other backend possible now.
 */
void derivs_prosse_sp_gslode(int neqs, double t, double *pars, 
			 const double *y, double *dydt) {
  do_derivs_prosse_sp(pars, y, dydt);
}
