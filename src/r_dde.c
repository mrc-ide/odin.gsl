#include "dde.h"
#include "sir.h"
#include <omp.h>

SEXP r_dde_example(SEXP r_n_replicates, SEXP r_steps, SEXP r_y_initial,
                   SEXP r_n_threads, SEXP user) {
  
  size_t n = 3;     // number of equations
  size_t n_out = 2; // number of "output" variables
  
  size_t n_replicates = INTEGER(r_n_replicates)[0];
  size_t n_steps = length(r_steps);
  // This is required to avoid the issue that int* and size_t* need
  // not be the same.  I should probably just relax this within run
  // instead
  size_t *steps = (size_t*) R_alloc(n_steps, sizeof(size_t));
  int* tmp = INTEGER(r_steps);
  for (size_t i = 0; i < n_steps; ++i) {
    steps[i] = (size_t) tmp[i];
  }

  size_t nt = n_steps;
  SEXP r_y = PROTECT(alloc3DArray(REALSXP, n, nt, n_replicates));
  SEXP r_out = PROTECT(alloc3DArray(REALSXP, n_out, nt, n_replicates));
  setAttrib(r_y, install("output"), r_out);

  double *y_initial = REAL(r_y_initial);
  double *y = REAL(r_y);
  double *out = REAL(r_out);

  size_t n_threads = INTEGER(r_n_threads);
  omp_set_num_threads(n_threads);
  #pragma omp parallel
  {
    SEXP sir = PROTECT(sir_create(user));
    sir_set_user(sir, user);

    sir_internal* data = sir_get_internal(sir, true);
    difeq_target* target = &sir_rhs_dde;

    // NOTE: this will leak on R error as the destructor will not be
    // called! That's fine for now, but usually I use a smart pointer
    // technique to prevent this.
    difeq_data* obj = difeq_data_alloc(target, n, n_out, data);

    //GetRNGstate();
    size_t i;
    #pragma omp for private(i) schedule(static)
    for (i = 0; i < n_replicates; ++i) {
      difeq_run(obj, y_initial, steps, n_steps, y, out);
      y += n * nt;
      out += n_out * nt;
    }
    //PutRNGstate();

    difeq_data_free(obj);

    UNPROTECT(3);
  }

  return r_y;
}
