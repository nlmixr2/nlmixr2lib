PK_1cmt <- function() {
  ini({
    ka <- 0.45 # Log Ka
    cl <- 1    # Log Cl
    v  <- 3.45  # log V
    prop.err <- 0.5
  })
  model({
    ka <- exp(ka)
    cl <- exp(cl)
    v  <- exp(v)
    linCmt() ~ prop(prop.err)
  })
}
