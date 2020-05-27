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

difeq_data* difeq_data_alloc(difeq_target * target,
                             size_t n, size_t n_out, const void *data,
                             size_t n_history, bool grow_history);
void difeq_data_free(difeq_data *obj);
void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out,
               bool return_initial);
