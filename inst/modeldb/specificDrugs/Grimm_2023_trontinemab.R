Grimm_2023_trontinemab <- function() {
  description <- "Trontinemab PK model (Grimm 2017)"
  reference <- "Grimm HP, Schumacher ,Vanessa, Schäfer ,Martin, et al. Delivery of the Brainshuttle™ amyloid-beta antibody fusion trontinemab to non-human primate brain and projected efficacious dose regimens in humans. mAbs. 2023;15(1):2261509. doi:10.1080/19420862.2023.2261509"
  units <- list(time = "hour", dosing = "mg/kg", dv = "ng/mL")
  covariateData <-
    list(
      WT = "Body weight (kg)"
    )
  ini({
    lcl <- log(1.01); label("Clearance (mL/h/kg)") # from table 1
    lvc <- log(45.4); label("Central volume of distribution (mL/kg)") # from table 1
    lvp <- log(63.2); label("Peripheral volume of distribution (mL/kg)") # from table 1
    lq <- log(0.746); label("Intercompartmental clearance (mL/h/kg)") # from table 1
    lvm <- log(78.8); label("Michaelis-Menten maximum clearance (mL/h/kg)") # CLmax from supplement
    lkm <- log(78.6); label("Michaelis-Menten half-maximal concentration (ng/mL)") # from supplement
    allo_cl <- fixed(0.85); label("Allometric exponent for clearance") # from text on page 11
    allo_v <- fixed(1); label("Allometric exponent for clearance") # from text on page 11
    
    CcpropSd <- 0; label("Proportional residual error (fraction)")
    CcaddSd <- 0; label("Additive residual error (ug/mL)")
  })
  model({
    cl <- exp(lcl + log(WT/70)*allo_cl)
    vc <- exp(lvc + log(WT/70)*allo_v)
    vp <- exp(lvp + log(WT/70)*allo_v)
    q <- exp(lq + log(WT/70)*allo_cl)
    vm <- exp(lvm)
    km <- exp(lkm)
    
    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(central) <- -kel*central - (vm * central/vc*1e6)/(km + central/vc*1e6) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # *1e6 for unit conversion from mg/mL to ng/mL
    Cc <- central / vc * 1e6
    
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
