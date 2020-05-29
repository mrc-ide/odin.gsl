#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// This is *very* similar to the format used by deSolve and dopri, but
// with a couple of differences:
//
//   double time has been replaced by size_t step
//
//   const void *data is passed through (vs deSolve)
//
//   output is always passed through (vs dopri)
//
typedef void difeq_target(size_t n_eq, size_t step,
                          const double *y, double *ynext,
                          size_t n_out, double *output, const void *data,
                          gsl_rng* rng);

typedef struct {
  difeq_target * target;
  const void *data;

  size_t n;     // number of equations
  size_t n_out; // number of output variables (possibly zero)

  size_t step0; // initial step number; used in models with delays
  size_t step;  // current step number
  size_t step1; // final step number

  // Steps for integration to report at
  const size_t *steps; // Set of steps to stop at
  size_t n_steps;
  size_t steps_idx;

  double * y0; // initial state
  double * y1; // final state
} difeq_data;

typedef struct {
  gsl_rng** generators;
  unsigned long* seeds;
  const gsl_rng_type* gen_type;
  const size_t n_generators;
} rng_array;

rng_array rng_init(const size_t n_generators,
                   const gsl_rng_type* gen_type,
                   unsigned long* seeds);
void rng_free(rng_array* rngs);

difeq_data* difeq_data_alloc(difeq_target * target,
                             size_t n, size_t n_out, const void *data);
void difeq_data_free(difeq_data *obj);
void difeq_data_reset(difeq_data *obj, const double *y,
                      const size_t *steps, size_t n_steps);
void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out, gsl_rng *rng);
