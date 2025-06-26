Grimm_2023_gantenerumab <- function() {
  description <- "Gantenerumab PK model (Grimm 2017)"
  reference <- "Grimm HP, Schumacher ,Vanessa, Schäfer ,Martin, et al. Delivery of the Brainshuttle™ amyloid-beta antibody fusion trontinemab to non-human primate brain and projected efficacious dose regimens in humans. mAbs. 2023;15(1):2261509. doi:10.1080/19420862.2023.2261509"
  units <-
    list(
      time = "hour",
      dosing = "mg/kg",
      Cc = "ng/mL",
      Ccerebellum = "ng/g",
      Cbrain_cerebellum = "ng/g",
      Chippocampus = "ng/g",
      Cbrain_hippocampus = "ng/g",
      Cstriatum = "ng/g",
      Cbrain_striatum = "ng/g",
      Ccortex = "ng/g",
      Cbrain_cortex = "ng/g",
      Cchoroid_plexus = "ng/g",
      Cbrain_choroid_plexus = "ng/g",
      Ccsf = "ng/g",
      Cbrain_csf = "ng/g"
    )
  covariateData <-
    list(
      WT = "Body weight (kg)"
    )
  ini({
    lvc <- log(50.8); label("Central volume of distribution (mL/kg)") # from table 1
    lcl <- log(0.537); label("Clearance (mL/h/kg)") # from table 1
    lvp <- log(97.1); label("Peripheral volume of distribution (mL/kg)") # from table 1
    lq <- log(2.86); label("Intercompartmental clearance (mL/h/kg)") # from table 1
    allo_cl <- fixed(0.85); label("Allometric exponent for clearance") # from text on page 11
    allo_v <- fixed(1); label("Allometric exponent for clearance") # from text on page 11
    
    # From supplementary Table 1
    kp_cerebellum <- 0.726e-3; label("Brain distribution coefficient, cerebellum (unitless)")
    kout_cerebellum <- 0.0624; label("Brain outflow rate, cerebellum (1/h)")
    fpla_cerebellum <- 1.26e-3; label("Residual plasma fraction, cerebellum")
    bsv_fpla_cerebellum ~ 1.93

    kp_hippocampus <- 0.298e-3; label("Brain distribution coefficient, hippocampus (unitless)")
    kout_hippocampus <- 0.0433; label("Brain outflow rate, hippocampus (1/h)")
    fpla_hippocampus <- 0.621e-3; label("Residual plasma fraction, hippocampus")
    bsv_fpla_hippocampus ~ 1.04
    
    kp_striatum <- 0.235e-3; label("Brain distribution coefficient, striatum (unitless)")
    kout_striatum <- 0.0371; label("Brain outflow rate, striatum (1/h)")
    fpla_striatum <- 0.298e-3; label("Residual plasma fraction, striatum")
    bsv_fpla_striatum ~ 1.96
    
    kp_cortex <- 0.406e-3; label("Brain distribution coefficient, cortex (unitless)")
    kout_cortex <- 0.0344; label("Brain outflow rate, cortex (1/h)")
    fpla_cortex <- 0.782e-3; label("Residual plasma fraction, cortex")
    bsv_fpla_cortex ~ 1.57
    
    kp_choroid_plexus <- 1.66e-3; label("Brain distribution coefficient, choroid plexus (unitless)")
    kout_choroid_plexus <- 0.0318; label("Brain outflow rate, choroid plexus (1/h)")
    fpla_choroid_plexus <- 18.4e-3; label("Residual plasma fraction, choroid plexus")
    bsv_fpla_choroid_plexus ~ 1.09
    
    kp_csf <- 3.00e-3; label("Brain distribution coefficient, cerebrospinal fluid (unitless)")
    kout_csf <- 0.0134; label("Brain outflow rate, cerebrospinal fluid (1/h)")
    fpla_csf <- fixed(0); label("Residual plasma fraction, cerebrospinal fluid")
    
    CcpropSd <- 0; label("Proportional residual error (fraction)")
    CcaddSd <- 0; label("Additive residual error (ug/mL)")
  })
  model({
    cl <- exp(lcl + log(WT/5)*allo_cl)
    vc <- exp(lvc + log(WT/5)*allo_v)
    vp <- exp(lvp + log(WT/5)*allo_v)
    q <- exp(lq + log(WT/5)*allo_cl)

    fpla_cerebellum_i <- fpla_cerebellum * exp(bsv_fpla_cerebellum)
    fpla_hippocampus_i <- fpla_hippocampus * exp(bsv_fpla_hippocampus)
    fpla_striatum_i <- fpla_striatum * exp(bsv_fpla_striatum)
    fpla_cortex_i <- fpla_cortex * exp(bsv_fpla_cortex)
    fpla_choroid_plexus_i <- fpla_choroid_plexus * exp(bsv_fpla_choroid_plexus)
    
    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(central) <- -kel*central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Unit conversion from mg/kg to ng/mL
    Cc <- central / vc * 1e6
    
    # Equation 1 (page 11): dCext/dt = kout*(Kp*Cpla - Cext)
    d/dt(Ccerebellum) <- kout_cerebellum*(kp_cerebellum*Cc - Ccerebellum)
    # Equation 2 (page 11): Cbrn = fpla*Cpla + Cext
    Cbrain_cerebellum <- fpla_cerebellum_i*Cc + Ccerebellum
    
    d/dt(Chippocampus) <- kout_hippocampus*(kp_hippocampus*Cc - Chippocampus)
    Cbrain_hippocampus <- fpla_hippocampus_i*Cc + Chippocampus
    
    d/dt(Cstriatum) <- kout_striatum*(kp_striatum*Cc - Cstriatum)
    Cbrain_striatum <- fpla_striatum_i*Cc + Cstriatum
    
    d/dt(Ccortex) <- kout_cortex*(kp_cortex*Cc - Ccortex)
    Cbrain_cortex <- fpla_cortex_i*Cc + Ccortex
    
    d/dt(Cchoroid_plexus) <- kout_choroid_plexus*(kp_choroid_plexus*Cc - Cchoroid_plexus)
    Cbrain_choroid_plexus <- fpla_choroid_plexus_i*Cc + Cchoroid_plexus
    
    d/dt(Ccsf) <- kout_csf*(kp_csf*Cc - Ccsf)
    Cbrain_csf <- fpla_csf*Cc + Ccsf
    
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
