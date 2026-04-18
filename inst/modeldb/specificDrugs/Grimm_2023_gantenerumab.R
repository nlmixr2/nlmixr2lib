Grimm_2023_gantenerumab <- function() {
  description <- "Gantenerumab PK model in cynomolgus monkeys (Grimm 2023): two-compartment plasma PK with brain extracellular distribution across six brain regions (cerebellum, hippocampus, striatum, cortex, choroid plexus, CSF)."
  reference <- "Grimm HP, Schumacher V, Schafer M, et al. Delivery of the Brainshuttle(TM) amyloid-beta antibody fusion trontinemab to non-human primate brain and projected efficacious dose regimens in humans. mAbs. 2023;15(1):2261509. doi:10.1080/19420862.2023.2261509"
  units <- list(time = "hour", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of CL, Vc, Vp, and Q; reference weight 5 kg (cynomolgus monkey).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    species        = "cynomolgus monkey",
    age_range      = "TODO: from source paper",
    weight_range   = "TODO: from source paper (reference weight 5 kg)",
    sex_female_pct = "TODO: from source paper",
    disease_state  = "Healthy non-human primates (cynomolgus monkey) used for gantenerumab brain-distribution PK characterization.",
    dose_range     = "TODO: from source paper",
    regions        = "Preclinical (non-human primate study)",
    notes          = "Model is the reference gantenerumab (non-Brainshuttle) arm of Grimm 2023, used for comparison against trontinemab. Structural parameters from Table 1; brain-region parameters from Supplementary Table 1."
  )

  ini({
    # Structural plasma PK parameters — Table 1 (reference WT = 5 kg cynomolgus monkey)
    lvc <- log(50.8); label("Central volume of distribution (mL/kg)")                          # Table 1
    lcl <- log(0.537); label("Clearance (mL/h/kg)")                                            # Table 1
    lvp <- log(97.1); label("Peripheral volume of distribution (mL/kg)")                       # Table 1
    lq <- log(2.86); label("Intercompartmental clearance (mL/h/kg)")                           # Table 1
    allo_cl <- fixed(0.85); label("Allometric exponent for clearance")                         # text on page 11
    allo_v <- fixed(1); label("Allometric exponent for volume")                                # text on page 11

    # Brain-region parameters — Supplementary Table 1
    # fpla_* are log-transformed so IIV follows the eta + transformed-name convention:
    # fpla_X = exp(lfpla_X + etalfpla_X) is mathematically equivalent to fpla_X * exp(eta).
    kp_cerebellum <- 0.726e-3; label("Brain distribution coefficient, cerebellum (unitless)")  # Supp. Table 1
    kout_cerebellum <- 0.0624; label("Brain outflow rate, cerebellum (1/h)")                   # Supp. Table 1
    lfpla_cerebellum <- log(1.26e-3); label("log residual plasma fraction, cerebellum")        # Supp. Table 1
    etalfpla_cerebellum ~ 1.93

    kp_hippocampus <- 0.298e-3; label("Brain distribution coefficient, hippocampus (unitless)") # Supp. Table 1
    kout_hippocampus <- 0.0433; label("Brain outflow rate, hippocampus (1/h)")                 # Supp. Table 1
    lfpla_hippocampus <- log(0.621e-3); label("log residual plasma fraction, hippocampus")     # Supp. Table 1
    etalfpla_hippocampus ~ 1.04

    kp_striatum <- 0.235e-3; label("Brain distribution coefficient, striatum (unitless)")      # Supp. Table 1
    kout_striatum <- 0.0371; label("Brain outflow rate, striatum (1/h)")                       # Supp. Table 1
    lfpla_striatum <- log(0.298e-3); label("log residual plasma fraction, striatum")           # Supp. Table 1
    etalfpla_striatum ~ 1.96

    kp_cortex <- 0.406e-3; label("Brain distribution coefficient, cortex (unitless)")          # Supp. Table 1
    kout_cortex <- 0.0344; label("Brain outflow rate, cortex (1/h)")                           # Supp. Table 1
    lfpla_cortex <- log(0.782e-3); label("log residual plasma fraction, cortex")               # Supp. Table 1
    etalfpla_cortex ~ 1.57

    kp_choroid_plexus <- 1.66e-3; label("Brain distribution coefficient, choroid plexus (unitless)") # Supp. Table 1
    kout_choroid_plexus <- 0.0318; label("Brain outflow rate, choroid plexus (1/h)")           # Supp. Table 1
    lfpla_choroid_plexus <- log(18.4e-3); label("log residual plasma fraction, choroid plexus") # Supp. Table 1
    etalfpla_choroid_plexus ~ 1.09

    kp_csf <- 3.00e-3; label("Brain distribution coefficient, cerebrospinal fluid (unitless)") # Supp. Table 1
    kout_csf <- 0.0134; label("Brain outflow rate, cerebrospinal fluid (1/h)")                 # Supp. Table 1
    fpla_csf <- fixed(0); label("Residual plasma fraction, cerebrospinal fluid")               # Supp. Table 1

    CcpropSd <- 0; label("Proportional residual error (fraction)")
    CcaddSd <- 0; label("Additive residual error (ng/mL)")
  })
  model({
    cl <- exp(lcl + log(WT/5)*allo_cl)
    vc <- exp(lvc + log(WT/5)*allo_v)
    vp <- exp(lvp + log(WT/5)*allo_v)
    q <- exp(lq + log(WT/5)*allo_cl)

    # Individual residual plasma fractions — exp(lfpla + etalfpla) keeps eta on log-scale
    # and matches the naming convention; equivalent to the original fpla * exp(bsv_fpla).
    fpla_cerebellum_i <- exp(lfpla_cerebellum + etalfpla_cerebellum)
    fpla_hippocampus_i <- exp(lfpla_hippocampus + etalfpla_hippocampus)
    fpla_striatum_i <- exp(lfpla_striatum + etalfpla_striatum)
    fpla_cortex_i <- exp(lfpla_cortex + etalfpla_cortex)
    fpla_choroid_plexus_i <- exp(lfpla_choroid_plexus + etalfpla_choroid_plexus)

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
