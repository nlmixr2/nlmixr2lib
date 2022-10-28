PK_1cmt <- function() {
  description <- "One compartment PK model with linear clearance"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)

    cp <- linCmt()
    cp ~ prop(propSd)
  })
}
