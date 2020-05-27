dde_example <- function(replicates, steps, y, user = NULL) {
  replicates <- as.integer(replicates)
  steps <- as.integer(steps)
  y <- as.double(y)
  .Call("r_dde_example", replicates, steps, y, user, PACKAGE = "odin.gsl")
}
