dde_example <- function(replicates, steps, y, threads, user = NULL) {
  replicates <- as.integer(replicates)
  steps <- as.integer(steps)
  y <- as.double(y)
  .Call("r_dde_example", replicates, steps, y, as.integer(threads), user, PACKAGE = "odin.gsl")
}
