# Parallel map with a text progress bar.
# Uses a FORK cluster on macOS/Linux and a PSOCK cluster on Windows.
#
# x         list (or vector) of tasks
# fun       function of one task; must be defined at top level
# n_cores   number of worker processes
# export    names of global objects that `fun` needs (used with PSOCK workers)
#
# Returns a list with one element per task; a task that fails returns its error.

.safe_call <- function(el, .fn) tryCatch(.fn(el), error = function(e) e)

parallel_map <- function(x, fun, n_cores, export = character(0),
                         chunk_size = 2 * n_cores,
                         type = if (.Platform$OS.type == "unix") "FORK" else "PSOCK") {
  cl <- parallel::makeCluster(n_cores, type = type)
  on.exit(parallel::stopCluster(cl))
  if (type == "PSOCK") {
    parallel::clusterEvalQ(cl, library(DCCA))
    parallel::clusterExport(cl, c(export, ".safe_call"), envir = globalenv())
  }
  pb <- txtProgressBar(max = length(x), style = 3)
  out <- vector("list", length(x))
  done <- 0
  for (idx in split(seq_along(x), ceiling(seq_along(x) / chunk_size))) {
    out[idx] <- parallel::parLapplyLB(cl, x[idx], .safe_call, .fn = fun)
    done <- done + length(idx)
    setTxtProgressBar(pb, done)
  }
  close(pb)
  out
}
