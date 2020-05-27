#include "dde.h"
void fill_na(double *x, size_t n);

difeq_data* difeq_data_alloc(difeq_target* target,
                             size_t n, size_t n_out, const void *data,
                             size_t n_history, bool grow_history) {
  difeq_data *ret = (difeq_data*) R_Calloc(1, difeq_data);
  overflow_action on_overflow =
    grow_history ? OVERFLOW_GROW : OVERFLOW_OVERWRITE;

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


void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out,
               bool return_initial) {
  difeq_data_reset(obj, y, steps, n_steps);

  double *y_next = NULL, *out_next = NULL;
  double *write_y = y_out, *write_out = out;

  bool has_output = obj->n_out > 0;
  bool store_next_output = false;
  bool store_y_in_history = obj->history_len > 0;

  if (store_y_in_history) {
    difeq_global_obj = obj;
    bool first_entry = ring_buffer_is_empty(obj->history);
    double *h = (double*) obj->history->head;
    memcpy(h + obj->history_idx_y, y, obj->n * sizeof(double));
    // Not 100% sure if we should do this, but it's probably not that unsafe.
    fill_na(h + obj->history_idx_out, obj->n_out);

    if (first_entry) {
      // No need to store the time as that should be checked ahead of
      // time
      difeq_store_time(obj);
      // Don't think that we should error data out here.
      // And we're off:
      h = ring_buffer_head_advance(obj->history);
    }

    y_next = h + obj->history_idx_y;
    out_next = y_next + obj->n;
  }

  // This is used to detect when the buffer has grown - because this
  // triggers a reallocation we can end up with pointers that are
  // invalidated and we need to update that.  This is only ever
  // accessed when 'store_y_in_history' is true.  Set this *after*
  // running advance above as otherwise it could be immediately
  // invalidated.
  const data_t * buffer_data = store_y_in_history ? obj->history->data : NULL;

  // If requested, copy initial conditions into the output space
  if (return_initial) {
    memcpy(write_y, y, obj->n * sizeof(double));
    store_next_output = true;
    write_y += obj->n;
  }

  if (!store_y_in_history) {
    y_next = write_y;
    out_next = write_out;
  }

  // This is used for scratch space in the case where we have output
  // (see the end of loop body).
  double *ytmp = has_output ? (double*) R_alloc(obj->n, sizeof(double)) : NULL;

  while (true) {
    // Try to be clever about where we store things to avoid copies
    // that will always happen.
    //
    // If we are using the ring buffer then storing 'y' anywhere but
    // the ring buffer will suffer a copy, so we might as well store
    // it there.
    //
    // Otherwise we'll use the final output array as a place to store
    // things.
    //
    // There is a *huge* bit of fencepost/off by one drama with
    // output; does output computed against variables y(i) go with
    // y(i) or with y(i + 1).  Hopefully this does not need to be
    // configurable because it's a real headache.
    obj->target(obj->n, obj->step, y, y_next, obj->n_out, out_next, obj->data);
    obj->step++;
    y = y_next;

    if (has_output && store_next_output) {
      if (store_y_in_history) {
        memcpy(write_out, out_next, obj->n_out * sizeof(double));
        out_next = y_next + obj->n;
      } else {
        out_next += obj->n_out;
      }
      write_out += obj->n_out;
      store_next_output = false;
    }

    if (obj->step == obj->steps[obj->steps_idx]) {
      if (store_y_in_history) {
        memcpy(write_y, y_next, obj->n * sizeof(double));
        y_next = ((double*) obj->history->head) + obj->history_idx_y;
      } else {
        y_next += obj->n;
      }
      write_y += obj->n;
      store_next_output = true;
      obj->steps_idx++;
    }

    if (store_y_in_history) {
      difeq_store_time(obj);
      double* h = (double*) ring_buffer_head_advance(obj->history);
      y_next = h + obj->history_idx_y;
      if (buffer_data != obj->history->data) {
        buffer_data = obj->history->data;
        // I think that this is always allowed because advance will
        // always move the head forward on grow, which is the only
        // case where this is triggered.
        y = y_next - obj->history_len;
        // This bit is less controversial
        out_next = y_next + obj->n_out;
      }
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
    if (store_y_in_history) { // NOTE: opposite direction to usual.
      memcpy(out_next, write_out, obj->n_out * sizeof(double));
    }
  }

  difeq_global_obj = NULL;
}


void fill_na(double *x, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    x[i] = NA_REAL;
  }
}
