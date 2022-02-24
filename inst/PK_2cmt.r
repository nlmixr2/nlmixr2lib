PK_2cmt <- function() {
  ini({
    ka <- 0.45  # Log Ka
    cl <- 1     # Log Cl
    v2 <- 3     # Log V2
    v3 <- 5     # Log V3
    q  <- 0.1   # Log Q
    prop.err <- 0.5
  })
  model({
    ka <- exp(ka)
    cl <- exp(cl)
    v2 <- exp(v2)
    v3 <- exp(v3)
    q  <- exp(q)
    linCmt() ~ prop(prop.err)
  })
}
