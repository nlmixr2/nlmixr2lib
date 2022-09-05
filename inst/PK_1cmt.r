# Description: One compartment PK model with linear clearance
PK_1cmt <- function() {
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lv  <- 3.45 ; label("Central volume of distribution (V)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    v  <- exp(lv)
    linCmt() ~ prop(prop.err)
  })
}
