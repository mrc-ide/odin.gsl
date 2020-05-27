#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <R_ext/Rdynload.h>
typedef struct sir_internal {
  double beta;
  double dt;
  double gamma;
  double I0;
  double initial_I;
  double initial_R;
  double initial_S;
  double p_IR;
  double S0;
  double steps_per_day;
} sir_internal;
sir_internal* sir_get_internal(SEXP internal_p, int closed_error);
static void sir_finalise(SEXP internal_p);
SEXP sir_create(SEXP user);
void sir_initmod_desolve(void(* odeparms) (int *, double *));
SEXP sir_contents(SEXP internal_p);
SEXP sir_set_user(SEXP internal_p, SEXP user);
SEXP sir_metadata(SEXP internal_p);
SEXP sir_initial_conditions(SEXP internal_p, SEXP step_ptr);
void sir_rhs(sir_internal* internal, size_t step, double * state, double * state_next, double * output);
void sir_rhs_dde(size_t n_eq, size_t step, double * state, double * state_next, size_t n_out, double * output, void * internal);
SEXP sir_rhs_r(SEXP internal_p, SEXP step, SEXP state);
double user_get_scalar_double(SEXP user, const char *name,
                              double default_value, double min, double max);
int user_get_scalar_int(SEXP user, const char *name,
                        int default_value, double min, double max);
void user_check_values_double(double * value, size_t len,
                                  double min, double max, const char *name);
void user_check_values_int(int * value, size_t len,
                               double min, double max, const char *name);
void user_check_values(SEXP value, double min, double max,
                           const char *name);
SEXP user_list_element(SEXP list, const char *name);
sir_internal* sir_get_internal(SEXP internal_p, int closed_error) {
  sir_internal *internal = NULL;
  if (TYPEOF(internal_p) != EXTPTRSXP) {
    Rf_error("Expected an external pointer");
  }
  internal = (sir_internal*) R_ExternalPtrAddr(internal_p);
  if (!internal && closed_error) {
    Rf_error("Pointer has been invalidated");
  }
  return internal;
}
void sir_finalise(SEXP internal_p) {
  sir_internal *internal = sir_get_internal(internal_p, 0);
  if (internal_p) {
    Free(internal);
    R_ClearExternalPtr(internal_p);
  }
}
SEXP sir_create(SEXP user) {
  sir_internal *internal = (sir_internal*) Calloc(1, sir_internal);
  internal->initial_R = 0;
  internal->beta = 0.20000000000000001;
  internal->gamma = 0.10000000000000001;
  internal->I0 = 10;
  internal->S0 = 1000;
  internal->steps_per_day = 4;
  SEXP ptr = PROTECT(R_MakeExternalPtr(internal, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(ptr, sir_finalise);
  UNPROTECT(1);
  return ptr;
}
static sir_internal *sir_internal_ds;
void sir_initmod_desolve(void(* odeparms) (int *, double *)) {
  static DL_FUNC get_desolve_gparms = NULL;
  if (get_desolve_gparms == NULL) {
    get_desolve_gparms =
      R_GetCCallable("deSolve", "get_deSolve_gparms");
  }
  sir_internal_ds = sir_get_internal(get_desolve_gparms(), 1);
}
SEXP sir_contents(SEXP internal_p) {
  sir_internal *internal = sir_get_internal(internal_p, 1);
  SEXP contents = PROTECT(allocVector(VECSXP, 10));
  SET_VECTOR_ELT(contents, 0, ScalarReal(internal->beta));
  SET_VECTOR_ELT(contents, 1, ScalarReal(internal->dt));
  SET_VECTOR_ELT(contents, 2, ScalarReal(internal->gamma));
  SET_VECTOR_ELT(contents, 3, ScalarReal(internal->I0));
  SET_VECTOR_ELT(contents, 4, ScalarReal(internal->initial_I));
  SET_VECTOR_ELT(contents, 5, ScalarReal(internal->initial_R));
  SET_VECTOR_ELT(contents, 6, ScalarReal(internal->initial_S));
  SET_VECTOR_ELT(contents, 7, ScalarReal(internal->p_IR));
  SET_VECTOR_ELT(contents, 8, ScalarReal(internal->S0));
  SET_VECTOR_ELT(contents, 9, ScalarReal(internal->steps_per_day));
  SEXP nms = PROTECT(allocVector(STRSXP, 10));
  SET_STRING_ELT(nms, 0, mkChar("beta"));
  SET_STRING_ELT(nms, 1, mkChar("dt"));
  SET_STRING_ELT(nms, 2, mkChar("gamma"));
  SET_STRING_ELT(nms, 3, mkChar("I0"));
  SET_STRING_ELT(nms, 4, mkChar("initial_I"));
  SET_STRING_ELT(nms, 5, mkChar("initial_R"));
  SET_STRING_ELT(nms, 6, mkChar("initial_S"));
  SET_STRING_ELT(nms, 7, mkChar("p_IR"));
  SET_STRING_ELT(nms, 8, mkChar("S0"));
  SET_STRING_ELT(nms, 9, mkChar("steps_per_day"));
  setAttrib(contents, R_NamesSymbol, nms);
  UNPROTECT(2);
  return contents;
}
SEXP sir_set_user(SEXP internal_p, SEXP user) {
  sir_internal *internal = sir_get_internal(internal_p, 1);
  internal->beta = user_get_scalar_double(user, "beta", internal->beta, NA_REAL, NA_REAL);
  internal->gamma = user_get_scalar_double(user, "gamma", internal->gamma, NA_REAL, NA_REAL);
  internal->I0 = user_get_scalar_double(user, "I0", internal->I0, NA_REAL, NA_REAL);
  internal->S0 = user_get_scalar_double(user, "S0", internal->S0, NA_REAL, NA_REAL);
  internal->steps_per_day = user_get_scalar_double(user, "steps_per_day", internal->steps_per_day, NA_REAL, NA_REAL);
  internal->dt = 1 / (double) internal->steps_per_day;
  internal->initial_I = internal->I0;
  internal->initial_S = internal->S0;
  internal->p_IR = 1 - exp(-(internal->gamma));
  return R_NilValue;
}
SEXP sir_metadata(SEXP internal_p) {
  sir_internal *internal = sir_get_internal(internal_p, 1);
  SEXP ret = PROTECT(allocVector(VECSXP, 4));
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, mkChar("variable_order"));
  SET_STRING_ELT(nms, 1, mkChar("output_order"));
  SET_STRING_ELT(nms, 2, mkChar("n_out"));
  SET_STRING_ELT(nms, 3, mkChar("interpolate_t"));
  setAttrib(ret, R_NamesSymbol, nms);
  SEXP variable_length = PROTECT(allocVector(VECSXP, 3));
  SEXP variable_names = PROTECT(allocVector(STRSXP, 3));
  setAttrib(variable_length, R_NamesSymbol, variable_names);
  SET_VECTOR_ELT(variable_length, 0, R_NilValue);
  SET_VECTOR_ELT(variable_length, 1, R_NilValue);
  SET_VECTOR_ELT(variable_length, 2, R_NilValue);
  SET_STRING_ELT(variable_names, 0, mkChar("S"));
  SET_STRING_ELT(variable_names, 1, mkChar("I"));
  SET_STRING_ELT(variable_names, 2, mkChar("R"));
  SET_VECTOR_ELT(ret, 0, variable_length);
  UNPROTECT(2);
  SEXP output_length = PROTECT(allocVector(VECSXP, 2));
  SEXP output_names = PROTECT(allocVector(STRSXP, 2));
  setAttrib(output_length, R_NamesSymbol, output_names);
  SET_VECTOR_ELT(output_length, 0, R_NilValue);
  SET_VECTOR_ELT(output_length, 1, R_NilValue);
  SET_STRING_ELT(output_names, 0, mkChar("incid"));
  SET_STRING_ELT(output_names, 1, mkChar("day"));
  SET_VECTOR_ELT(ret, 1, output_length);
  UNPROTECT(2);
  SET_VECTOR_ELT(ret, 2, ScalarInteger(2));
  UNPROTECT(2);
  return ret;
}
SEXP sir_initial_conditions(SEXP internal_p, SEXP step_ptr) {
  sir_internal *internal = sir_get_internal(internal_p, 1);
  SEXP r_state = PROTECT(allocVector(REALSXP, 3));
  double * state = REAL(r_state);
  state[0] = internal->initial_S;
  state[1] = internal->initial_I;
  state[2] = internal->initial_R;
  UNPROTECT(1);
  return r_state;
}
void sir_rhs(sir_internal* internal, size_t step, double * state, double * state_next, double * output) {
  double S = state[0];
  double I = state[1];
  double R = state[2];
  double N = S + I + R;
  double n_IR = Rf_rbinom(round(I), internal->p_IR * internal->dt);
  double p_SI = 1 - exp(-(internal->beta) * I / (double) N);
  double n_SI = Rf_rbinom(round(S), p_SI * internal->dt);
  state_next[2] = R + n_IR;
  state_next[1] = I + n_SI - n_IR;
  state_next[0] = S - n_SI;
  output[1] = step * internal->dt;
  output[0] = n_SI;
}
void sir_rhs_dde(size_t n_eq, size_t step, double * state, double * state_next, size_t n_out, double * output, void * internal) {
  sir_rhs((sir_internal*)internal, step, state, state_next, output);
}
SEXP sir_rhs_r(SEXP internal_p, SEXP step, SEXP state) {
  SEXP state_next = PROTECT(allocVector(REALSXP, LENGTH(state)));
  sir_internal *internal = sir_get_internal(internal_p, 1);
  SEXP output_ptr = PROTECT(allocVector(REALSXP, 2));
  setAttrib(state_next, install("output"), output_ptr);
  UNPROTECT(1);
  double *output = REAL(output_ptr);
  GetRNGstate();
  sir_rhs(internal, INTEGER(step)[0], REAL(state), REAL(state_next), output);
  PutRNGstate();
  UNPROTECT(1);
  return state_next;
}
double user_get_scalar_double(SEXP user, const char *name,
                              double default_value, double min, double max) {
  double ret = default_value;
  SEXP el = user_list_element(user, name);
  if (el != R_NilValue) {
    if (length(el) != 1) {
      Rf_error("Expected a scalar numeric for '%s'", name);
    }
    if (TYPEOF(el) == REALSXP) {
      ret = REAL(el)[0];
    } else if (TYPEOF(el) == INTSXP) {
      ret = INTEGER(el)[0];
    } else {
      Rf_error("Expected a numeric value for %s", name);
    }
  }
  if (ISNA(ret)) {
    Rf_error("Expected a value for '%s'", name);
  }
  user_check_values_double(&ret, 1, min, max, name);
  return ret;
}
int user_get_scalar_int(SEXP user, const char *name,
                        int default_value, double min, double max) {
  int ret = default_value;
  SEXP el = user_list_element(user, name);
  if (el != R_NilValue) {
    if (length(el) != 1) {
      Rf_error("Expected scalar integer for %d", name);
    }
    if (TYPEOF(el) == REALSXP) {
      double tmp = REAL(el)[0];
      if (fabs(tmp - round(tmp)) > 2e-8) {
        Rf_error("Expected '%s' to be integer-like", name);
      }
    }
    ret = INTEGER(coerceVector(el, INTSXP))[0];
  }
  if (ret == NA_INTEGER) {
    Rf_error("Expected a value for '%s'", name);
  }
  user_check_values_int(&ret, 1, min, max, name);
  return ret;
}
void user_check_values_double(double * value, size_t len,
                                  double min, double max, const char *name) {
  for (size_t i = 0; i < len; ++i) {
    if (ISNA(value[i])) {
      Rf_error("'%s' must not contain any NA values", name);
    }
  }
  if (min != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] < min) {
        Rf_error("Expected '%s' to be at least %g", name, min);
      }
    }
  }
  if (max != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] > max) {
        Rf_error("Expected '%s' to be at most %g", name, max);
      }
    }
  }
}
void user_check_values_int(int * value, size_t len,
                               double min, double max, const char *name) {
  for (size_t i = 0; i < len; ++i) {
    if (ISNA(value[i])) {
      Rf_error("'%s' must not contain any NA values", name);
    }
  }
  if (min != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] < min) {
        Rf_error("Expected '%s' to be at least %g", name, min);
      }
    }
  }
  if (max != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] > max) {
        Rf_error("Expected '%s' to be at most %g", name, max);
      }
    }
  }
}
void user_check_values(SEXP value, double min, double max,
                           const char *name) {
  size_t len = (size_t)length(value);
  if (TYPEOF(value) == INTSXP) {
    user_check_values_int(INTEGER(value), len, min, max, name);
  } else {
    user_check_values_double(REAL(value), len, min, max, name);
  }
}
SEXP user_list_element(SEXP list, const char *name) {
  SEXP ret = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); ++i) {
    if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      ret = VECTOR_ELT(list, i);
      break;
    }
  }
  return ret;
}
