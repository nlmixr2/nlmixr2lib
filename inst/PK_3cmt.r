# Description: Three compartment PK model with linear clearance
PK_3cmt <- function() {
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lv  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lvp2  <- 8 ; label("Second peripheral volume of distribution (Vp2)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    lq2  <- 0.5 ; label("Second intercompartmental clearance (Q2)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    v <- exp(lv)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q  <- exp(lq)
    q2  <- exp(lq2)
    linCmt() ~ prop(prop.err)
  })
}
