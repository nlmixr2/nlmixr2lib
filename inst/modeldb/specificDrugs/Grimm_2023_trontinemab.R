Grimm_2023_trontinemab <- function() {
  description <- "Trontinemab PK model in non-human primates (Grimm 2023): two-compartment plasma PK with Michaelis-Menten elimination and brain-region effect-compartment distribution (cerebellum, hippocampus, striatum, cortex, choroid plexus, CSF)."
  reference <- "Grimm HP, Schumacher V, Schafer M, et al. Delivery of the Brainshuttle amyloid-beta antibody fusion trontinemab to non-human primate brain and projected efficacious dose regimens in humans. mAbs. 2023;15(1):2261509. doi:10.1080/19420862.2023.2261509"
  vignette <- "Grimm_2023"
  units <- list(
    time = "hour",
    dosing = "mg/kg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling; reference weight 5 kg (cynomolgus monkey). Allometric exponent 0.85 on clearances, 1 on volumes.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    age_range      = "TODO: from source paper",
    weight_range   = "TODO: from source paper",
    sex_female_pct = "TODO: from source paper",
    species        = "Cynomolgus monkey (Macaca fascicularis)",
    disease_state  = "Healthy non-human primates (no amyloid pathology)",
    dose_range     = "TODO: from source paper (intravenous dosing, mg/kg)",
    regions        = "Preclinical (non-human primate study)",
    notes          = "Non-human primate (NHP) dataset used to characterize plasma PK and brain-region distribution of trontinemab. Reference body weight 5 kg. Demographic details (exact N, age, weight, sex distribution) TODO: extract from Grimm 2023 methods / supplement."
  )

  ini({
    # Structural plasma PK parameters — Grimm 2023 Table 1 (reference WT = 5 kg)
    lvc <- log(45.4);  label("Central volume of distribution (mL/kg)")                   # Table 1
    lcl <- log(1.01);  label("Clearance (mL/h/kg)")                                      # Table 1
    lvp <- log(63.2);  label("Peripheral volume of distribution (mL/kg)")                # Table 1
    lq  <- log(0.746); label("Intercompartmental clearance (mL/h/kg)")                   # Table 1
    lvm <- log(78.8);  label("Michaelis-Menten maximum clearance (mL/h/kg)")             # CLmax from supplement
    lkm <- log(78.6);  label("Michaelis-Menten half-maximal concentration (ng/mL)")      # from supplement

    # Allometric exponents — fixed per Grimm 2023 text on page 11
    allo_cl <- fixed(0.85); label("Allometric exponent for clearance (unitless)")        # page 11
    allo_v  <- fixed(1);    label("Allometric exponent for volume (unitless)")           # page 11

    # Brain-region effect-compartment parameters — Grimm 2023 supplementary Table 1
    # Cerebellum
    kp_cerebellum   <- 5.17e-3; label("Brain distribution coefficient, cerebellum (unitless)")  # supp. Table 1
    kout_cerebellum <- 0.0624;  label("Brain outflow rate, cerebellum (1/h)")                    # supp. Table 1
    lfpla_cerebellum <- log(1.26e-3); label("Residual plasma fraction, cerebellum (log)")        # supp. Table 1: fpla = 1.26e-3
    etalfpla_cerebellum ~ 1.93                                                                   # supp. Table 1: IIV variance on fpla_cerebellum

    # Hippocampus
    kp_hippocampus   <- 5.59e-3;  label("Brain distribution coefficient, hippocampus (unitless)") # supp. Table 1
    kout_hippocampus <- 0.0433;   label("Brain outflow rate, hippocampus (1/h)")                  # supp. Table 1
    lfpla_hippocampus <- log(0.621e-3); label("Residual plasma fraction, hippocampus (log)")      # supp. Table 1: fpla = 0.621e-3
    etalfpla_hippocampus ~ 1.04                                                                   # supp. Table 1

    # Striatum
    kp_striatum   <- 7.76e-3; label("Brain distribution coefficient, striatum (unitless)")        # supp. Table 1
    kout_striatum <- 0.0371;  label("Brain outflow rate, striatum (1/h)")                         # supp. Table 1
    lfpla_striatum <- log(0.298e-3); label("Residual plasma fraction, striatum (log)")            # supp. Table 1: fpla = 0.298e-3
    etalfpla_striatum ~ 1.96                                                                      # supp. Table 1

    # Cortex
    kp_cortex   <- 4.62e-3; label("Brain distribution coefficient, cortex (unitless)")            # supp. Table 1
    kout_cortex <- 0.0344;  label("Brain outflow rate, cortex (1/h)")                             # supp. Table 1
    lfpla_cortex <- log(0.782e-3); label("Residual plasma fraction, cortex (log)")                # supp. Table 1: fpla = 0.782e-3
    etalfpla_cortex ~ 1.57                                                                        # supp. Table 1

    # Choroid plexus
    kp_choroid_plexus   <- 4.31e-3; label("Brain distribution coefficient, choroid plexus (unitless)") # supp. Table 1
    kout_choroid_plexus <- 0.0318;  label("Brain outflow rate, choroid plexus (1/h)")                  # supp. Table 1
    lfpla_choroid_plexus <- log(18.4e-3); label("Residual plasma fraction, choroid plexus (log)")      # supp. Table 1: fpla = 18.4e-3
    etalfpla_choroid_plexus ~ 1.09                                                                     # supp. Table 1

    # Cerebrospinal fluid (no residual plasma fraction)
    kp_csf   <- 2.34e-3; label("Brain distribution coefficient, cerebrospinal fluid (unitless)")  # supp. Table 1
    kout_csf <- 0.0373;  label("Brain outflow rate, cerebrospinal fluid (1/h)")                   # supp. Table 1
    fpla_csf <- fixed(0); label("Residual plasma fraction, cerebrospinal fluid")                  # supp. Table 1: fpla_csf fixed at 0

    # Residual error — source paper does not report residual error magnitudes; left at 0 as placeholders
    CcpropSd <- 0; label("Proportional residual error, plasma Cc (fraction)")
    CcaddSd  <- 0; label("Additive residual error, plasma Cc (ng/mL)")
  })

  model({
    # Individual plasma PK parameters with allometric scaling to 5 kg
    cl <- exp(lcl + log(WT / 5) * allo_cl)
    vc <- exp(lvc + log(WT / 5) * allo_v)
    vp <- exp(lvp + log(WT / 5) * allo_v)
    q  <- exp(lq  + log(WT / 5) * allo_cl)
    vm <- exp(lvm)
    km <- exp(lkm)

    # Individual residual plasma fractions per brain region (log-normal IIV)
    fpla_cerebellum     <- exp(lfpla_cerebellum     + etalfpla_cerebellum)
    fpla_hippocampus    <- exp(lfpla_hippocampus    + etalfpla_hippocampus)
    fpla_striatum       <- exp(lfpla_striatum       + etalfpla_striatum)
    fpla_cortex         <- exp(lfpla_cortex         + etalfpla_cortex)
    fpla_choroid_plexus <- exp(lfpla_choroid_plexus + etalfpla_choroid_plexus)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Plasma two-compartment model with parallel linear + Michaelis-Menten elimination
    d/dt(central)     <- -kel * central - central * (vm / vc) / (1 + (central / vc * 1e6) / km) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration: unit conversion from mg/kg compartment amount to ng/mL
    Cc <- central / vc * 1e6

    # Brain-region effect compartments — Grimm 2023 page 11
    # Equation 1: dCext/dt = kout * (Kp * Cpla - Cext)
    # Equation 2: Cbrn     = fpla * Cpla + Cext
    d/dt(Ccerebellum)     <- kout_cerebellum     * (kp_cerebellum     * Cc - Ccerebellum)
    Cbrain_cerebellum     <- fpla_cerebellum     * Cc + Ccerebellum

    d/dt(Chippocampus)    <- kout_hippocampus    * (kp_hippocampus    * Cc - Chippocampus)
    Cbrain_hippocampus    <- fpla_hippocampus    * Cc + Chippocampus

    d/dt(Cstriatum)       <- kout_striatum       * (kp_striatum       * Cc - Cstriatum)
    Cbrain_striatum       <- fpla_striatum       * Cc + Cstriatum

    d/dt(Ccortex)         <- kout_cortex         * (kp_cortex         * Cc - Ccortex)
    Cbrain_cortex         <- fpla_cortex         * Cc + Ccortex

    d/dt(Cchoroid_plexus) <- kout_choroid_plexus * (kp_choroid_plexus * Cc - Cchoroid_plexus)
    Cbrain_choroid_plexus <- fpla_choroid_plexus * Cc + Cchoroid_plexus

    d/dt(Ccsf)            <- kout_csf            * (kp_csf            * Cc - Ccsf)
    Cbrain_csf            <- fpla_csf            * Cc + Ccsf

    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
