# obtain object name for printing in figure titles
string <- function(o) {
  newname <- substitute(o)
  # Travel up the frame stack until we hit the top.
  for(i in seq_len(sys.nframe())) {
    oldname <- do.call("substitute", list(as.name(newname), parent.frame(i)))
    newname <- oldname
  }
  deparse(newname)
}