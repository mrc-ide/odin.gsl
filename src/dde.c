#include "dde.h"
#include <R.h>

void fill_na(double *x, size_t n);

difeq_data* difeq_data_alloc(difeq_target* target,
                             size_t n, size_t n_out, const void *data) {
  difeq_data *ret = (difeq_data*) R_Calloc(1, difeq_data);

  ret->target = target;
  ret->data = data;

  ret->n = n;
  ret->n_out = n_out;

  ret->n_steps = 0;
  ret->steps = NULL;

  // State vectors
  ret->y0 = R_Calloc(n, double); // initial
  ret->y1 = R_Calloc(n, double); // final (only used on restart?)

  return ret;
}


void difeq_data_free(difeq_data *obj) {
  R_Free(obj->y0);
  R_Free(obj->y1);
  R_Free(obj);
}


void difeq_data_reset(difeq_data *obj, const double *y,
                      const size_t *steps, size_t n_steps) {
  memcpy(obj->y0, y, obj->n * sizeof(double));

  obj->n_steps = n_steps;
  obj->steps = steps;
  obj->steps_idx = 1; // skipping the first step!

  // TODO: I don't check that there is at least one time anywhere in
  // *this* routine, but it is checked in r_difeq which calls
  // this.
  if (steps[n_steps - 1] == steps[0]) {
    Rf_error("Initialisation failure: Beginning and end steps are the same");
  }

  for (size_t i = 0; i < n_steps - 1; ++i) {
    if (steps[i + 1] <= steps[i]) { // NOTE: Disallows ties
      Rf_error("Initialisation failure: Steps not strictly increasing");
    }
  }

  obj->step0 = steps[0];
  obj->step = steps[0];
  obj->step1 = steps[n_steps - 1];
}


void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out) {
  difeq_data_reset(obj, y, steps, n_steps);

  double *y_next = NULL, *out_next = NULL;
  double *write_y = y_out, *write_out = out;

  bool has_output = obj->n_out > 0;
  bool store_next_output = false;

  memcpy(write_y, y, obj->n * sizeof(double));
  store_next_output = true;
  write_y += obj->n;

  y_next = write_y;
  out_next = write_out;

  // This is used for scratch space in the case where we have output
  // (see the end of loop body).
  double *ytmp = has_output ? (double*) R_alloc(obj->n, sizeof(double)) : NULL;

  while (true) {
    obj->target(obj->n, obj->step, y, y_next, obj->n_out, out_next, obj->data);
    obj->step++;
    y = y_next;

    if (has_output && store_next_output) {
      out_next += obj->n_out;
      write_out += obj->n_out;
      store_next_output = false;
    }

    if (obj->step == obj->steps[obj->steps_idx]) {
      y_next += obj->n;
      write_y += obj->n;
      store_next_output = true;
      obj->steps_idx++;
    }

    if (obj->step == obj->step1) {
      // The final 'y' state is used on restart
      memcpy(obj->y1, y, obj->n * sizeof(double));
      break;
    }
  }

  if (has_output && store_next_output) {
    obj->target(obj->n, obj->step, y, ytmp, obj->n_out, write_out,
                obj->data);
  }
}


void fill_na(double *x, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    x[i] = NA_REAL;
  }
}
