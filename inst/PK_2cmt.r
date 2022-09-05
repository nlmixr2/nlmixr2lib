# Description: Two compartment PK model with linear clearance
PK_2cmt <- function() {
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lv  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    v <- exp(lv)
    vp <- exp(lvp)
    q  <- exp(lq)
    linCmt() ~ prop(prop.err)
  })
}
