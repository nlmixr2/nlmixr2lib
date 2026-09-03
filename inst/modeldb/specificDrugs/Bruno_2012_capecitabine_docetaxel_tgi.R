Bruno_2012_capecitabine_docetaxel_tgi <- function() {
  description <- paste(
    "Tumor growth inhibition (TGI) model for the capecitabine plus docetaxel",
    "combination in second-line locally advanced / metastatic breast cancer",
    "(Bruno 2012, pooled studies SO14999 and NO16853). Tumor size (sum of the",
    "longest diameters of measurable / target lesions, mm) grows exponentially",
    "at rate KL and is reduced by two additive, independent cell-kill terms,",
    "one per drug. Each drug's dose enters its own K-PD virtual biophase",
    "compartment with first-order elimination rate KP,X; the elimination flux",
    "KP,X * depot_kpd_X is the 'virtual infusion rate' that drives the kill",
    "term, and each drug's cell-kill coefficient decays exponentially from",
    "KD,0,X with the resistance-appearance rate lambda_X. There is no PK model:",
    "dose is the exposure metric. Doses are administered amounts in mg; the",
    "capecitabine kill coefficient is expressed per gram, so its biophase flux",
    "is divided by 1000. The additive residual error is study-specific.",
    "Time is in weeks. Companion models: modellib('Bruno_2012_capecitabine_docetaxel_os')",
    "and modellib('Bruno_2012_capecitabine_docetaxel_pfs').",
    sep = " "
  )
  reference <- paste(
    "Bruno R, Lindbom L, Schaedeli Stark F, Chanu P, Gilberg F, Frey N, Claret L.",
    "Simulations to assess phase II noninferiority trials of different doses of",
    "capecitabine in combination with docetaxel for metastatic breast cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2012;1(12):e19.",
    "doi:10.1038/psp.2012.20.",
    "Structural equations, parameter definitions and the NONMEM control stream",
    "are in the Supplementary Information; final parameter estimates are in",
    "Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Bruno_2012_capecitabine_docetaxel_mbc"

  units <- list(
    time          = "week",
    dosing        = "mg",
    concentration = "mm (the observed variable is tumor size, the sum of the longest diameters, not a drug concentration)"
  )

  compartmentData <- list(
    depot_kpd_docetaxel = list(
      analyte  = "docetaxel",
      units    = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    depot_kpd_capecitabine = list(
      analyte  = "capecitabine",
      units    = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_size = list(
      analyte  = "tumor size (sum of the longest diameters of measurable / target lesions)",
      units    = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUM_SLD = list(
      description        = "Observed baseline tumor size: the sum of the longest diameters of all measurable lesions (WHO criteria, study SO14999) or of the target lesions (RECIST, study NO16853).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Sets the per-subject initial condition of the tumor-size ODE (`tumor_size(0) <- TUM_SLD`), reproducing the source NONMEM `A_0(3) = BASE`. Bruno 2012 explicitly kept the observed baseline out of the dataset as an observation and supplied it as a covariate instead (Supplementary Information). Observed range in the pooled 888-patient dataset: 10 to 520 mm (Bruno 2012 Results, TGI model).",
      source_name        = "BASE"
    ),
    STUDY_NO16853 = list(
      description        = "Study indicator for the randomized phase II noninferiority study NO16853: 1 = NO16853, 0 = the pivotal phase III study SO14999.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study SO14999)",
      notes              = "Time-fixed per subject. Selects the study-specific additive residual error variance: 332 mm^2 in SO14999 versus 112 mm^2 in NO16853 (Supplementary Table S1). The source NONMEM control stream codes this as `IF(STUD.EQ.2) REE = SQRT(THETA(9))`, i.e. STUD = 2 is NO16853. Bruno 2012 attributes the higher residual variance in SO14999 to that study's older WHO response criteria and measurement practice. Distinct from the `1 / 2`-coded study effect on survival and PFS time, which is carried by the same indicator in the companion models `Bruno_2012_capecitabine_docetaxel_os` and `Bruno_2012_capecitabine_docetaxel_pfs`.",
      source_name        = "STUD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 888L,
    n_studies      = 2L,
    n_observations = 2988L,
    disease_state  = "locally advanced / metastatic breast cancer, second line after anthracycline pretreatment",
    dose_range     = "capecitabine 825 or 1,250 mg/m^2 twice daily on days 1-14 of each 3-week cycle plus docetaxel 75 mg/m^2 on day 1 of each cycle; study SO14999 additionally contributed a single-agent docetaxel every-3-weeks arm (Bruno 2012 does not state that arm's docetaxel dose)",
    notes          = "Pooled longitudinal tumor-size database of 888 patients (463 from the pivotal phase III study SO14999, 425 from the randomized phase II noninferiority study NO16853) contributing 2,988 tumor-size observations, a mean of 3.4 measurements per patient (Bruno 2012 Results, TGI model). Baseline tumor size ranged from 10 to 520 mm. Tumor response was assessed with WHO criteria in SO14999 and with RECIST in NO16853; the sum of the longest diameters of all measurable lesions (WHO) or of the target lesions (RECIST) is the modelled variable. Bruno 2012 does not tabulate age, weight, sex or race for this cohort."
  )

  ini({
    # Structural parameters. Supplementary Table S1 ("Tumor Growth Inhibition
    # Final Model Parameter Estimates") reports each value with its RSE; the
    # source NONMEM $THETA block (Supplementary Information) names the same
    # parameters KL, KD1, KD2, DM1, DM2, KP1, KP2. All rate constants are in
    # week^-1 or week^-1 per unit virtual infusion rate; time is in weeks.
    lkgrow <- log(0.00437)
    label("Tumor growth rate KL (1/week)")                                                          # Bruno 2012 Supplementary Table S1: KL = 0.00437 /week (RSE 15%); source $THETA(1)

    lkkill_docetaxel <- log(0.00128)
    label("Docetaxel cell-kill rate coefficient KD,0,Doc at time zero (per mg/week of virtual infusion rate)") # Bruno 2012 Supplementary Table S1: KD,0,Doc = 0.00128 (RSE 31%); source $THETA(2) = KD1

    lkkill_capecitabine <- log(0.00470)
    label("Capecitabine cell-kill rate coefficient KD,0,Cap at time zero (per g/week of virtual infusion rate)") # Bruno 2012 Supplementary Table S1: KD,0,Cap = 0.00470 (RSE 20%); source $THETA(3) = KD2

    lres_docetaxel <- log(0.0450)
    label("Docetaxel resistance-appearance rate lambda_Doc (1/week)")                               # Bruno 2012 Supplementary Table S1: lambda_Doc = 0.0450 /week (RSE 39%, the least well estimated parameter, cited in Results); source $THETA(4) = DM1

    lres_capecitabine <- log(0.240)
    label("Capecitabine resistance-appearance rate lambda_Cap (1/week)")                            # Bruno 2012 Supplementary Table S1: lambda_Cap = 0.240 /week (RSE 23%); source $THETA(5) = DM2

    # K-PD biophase elimination rate constants. Supplementary Table S1 reports
    # both as "Fixed" in place of an RSE, and the source $THETA block carries
    # the NONMEM FIX flag on both: "0.1 FIXED ; KP DOCE" and "3 FIXED ; KP CAPE".
    # The Table S1 footnote states they "were not identifiable and were fixed to
    # optimal values based on likelihood profiling".
    lkel_docetaxel <- fixed(log(0.1))
    label("Docetaxel K-PD biophase elimination rate constant KP,Doc (1/week)")                      # Bruno 2012 Supplementary Table S1: KP,Doc = 0.1 /week, Fixed; source $THETA(7) = 0.1 FIXED

    lkel_capecitabine <- fixed(log(3))
    label("Capecitabine K-PD biophase elimination rate constant KP,Cap (1/week)")                   # Bruno 2012 Supplementary Table S1: KP,Cap = 3 /week, Fixed; source $THETA(8) = 3 FIXED

    # IIV. Supplementary Information: "Inter-patient variability of parameters
    # was modeled using exponential random effects", omega^2 reported directly
    # (not as CV%) in Supplementary Table S1. The source $OMEGA block is
    # "0.5 0.5 0.5 0 FIXED 0 FIXED 1 FIXED": the first three are the estimated
    # variances on KL, KD1 and KD2 (initial values 0.5), and the two 0-FIXED
    # entries are ETA(4) / ETA(5) on DM1 / DM2, i.e. the model carries NO
    # inter-patient variability on either resistance-appearance rate. Bruno 2012
    # states "No covariance between the random effects was considered (diagonal
    # covariance matrix)", so no block is used.
    etalkgrow              ~ 1.60                                                                   # Bruno 2012 Supplementary Table S1: omega^2_KL = 1.60 (RSE 5%)
    etalkkill_docetaxel    ~ 1.31                                                                   # Bruno 2012 Supplementary Table S1: omega^2_KD,0,Doc = 1.31 (RSE 15%)
    etalkkill_capecitabine ~ 1.30                                                                   # Bruno 2012 Supplementary Table S1: omega^2_KD,0,Cap = 1.30 (RSE 19%)

    # Residual error. Supplementary Information: "Residual variability was
    # modeled using an additive error model". Supplementary Table S1 reports the
    # study-specific residual VARIANCES sigma^2 in mm^2; the source $ERROR block
    # takes their square root ("REE = SQRT(THETA(6))" / "REE = SQRT(THETA(9))")
    # with $SIGMA 1 FIXED, so the values below are the corresponding SDs.
    addSdStudySO14999 <- sqrt(332)
    label("Additive residual SD on tumor size in study SO14999 (mm)")                               # Bruno 2012 Supplementary Table S1: sigma^2_Study SO14999 = 332 mm^2 (RSE 11%) -> SD = sqrt(332) = 18.22 mm; source $THETA(6)

    addSdStudyNO16853 <- sqrt(112)
    label("Additive residual SD on tumor size in study NO16853 (mm)")                               # Bruno 2012 Supplementary Table S1: sigma^2_Study NO16853 = 112 mm^2 (RSE 18%) -> SD = sqrt(112) = 10.58 mm; source $THETA(9)
  })
  model({
    # 1. Individual parameters (exponential inter-patient random effects;
    #    source $PK: KL = THETA(1)*EXP(ETA(1)), KD1 = THETA(2)*EXP(ETA(2)),
    #    KD2 = THETA(3)*EXP(ETA(3)), DM1 = THETA(4)*EXP(ETA(4)) with
    #    ETA(4) / ETA(5) fixed to zero variance).
    kgrow              <- exp(lkgrow + etalkgrow)
    kkill_docetaxel    <- exp(lkkill_docetaxel + etalkkill_docetaxel)
    kkill_capecitabine <- exp(lkkill_capecitabine + etalkkill_capecitabine)
    res_docetaxel      <- exp(lres_docetaxel)
    res_capecitabine   <- exp(lres_capecitabine)
    kel_docetaxel      <- exp(lkel_docetaxel)
    kel_capecitabine   <- exp(lkel_capecitabine)

    # 2. Time since the first dose, floored at zero. Source $DES:
    #    "TTD = TFD; IF (TFD.LT.0) TTD = 0", where TFD is time from first dose.
    #    In a forward simulation started at the first dose, ttd equals t.
    ttd <- ifelse(t < 0, 0, t)

    # 3. Virtual infusion rates. The K-PD biophase elimination flux
    #    KP,X * depot_kpd_X is the "virtual infusion rate" IR_X that drives the
    #    kill term ("In this implementation, the kill rate is proportional to
    #    the dose present in the biophase compartment", Supplementary
    #    Information). The capecitabine flux is divided by 1000 to convert mg to
    #    g, matching the g^-1 unit of KD,0,Cap in Supplementary Table S1; source
    #    $DES: "DR2 = A(2)/1000 * KP2 * EXP(-DM2*TTD)".
    ir_docetaxel    <- kel_docetaxel * depot_kpd_docetaxel
    ir_capecitabine <- kel_capecitabine * depot_kpd_capecitabine / 1000

    # 4. ODE system. Source $DES:
    #      DADT(1) = -KP1 * A(1)
    #      DADT(2) = -KP2 * A(2)
    #      DADT(3) =  KL*A(3) - KD1*DR1*A(3) - KD2*DR2*A(3)
    #    with the two drug effects additive, as Bruno 2012 Methods states
    #    ("assuming additive effect of capecitabine and docetaxel").
    d/dt(depot_kpd_docetaxel)    <- -kel_docetaxel * depot_kpd_docetaxel
    d/dt(depot_kpd_capecitabine) <- -kel_capecitabine * depot_kpd_capecitabine
    d/dt(tumor_size) <- kgrow * tumor_size -
      kkill_docetaxel * ir_docetaxel * exp(-res_docetaxel * ttd) * tumor_size -
      kkill_capecitabine * ir_capecitabine * exp(-res_capecitabine * ttd) * tumor_size

    # 5. Initial conditions. Source $PK: "A_0(3) = BASE; A_0(1) = 0; A_0(2) = 0".
    #    The Supplementary Information also describes an ATTEMPT to add a
    #    baseline random effect eta_i0 with the same variance as the residual
    #    error; that component is not referenced in the printed $PK block and no
    #    corresponding variance appears in Supplementary Table S1, so the final
    #    model encoded here initialises the tumor-size state at the observed
    #    baseline without a random effect.
    tumor_size(0) <- TUM_SLD

    # 6. Study-specific additive residual error (source $ERROR:
    #    "Y = F + REE * ERR(1)" with REE selected by the STUD flag).
    addSd <- (1 - STUDY_NO16853) * addSdStudySO14999 +
      STUDY_NO16853 * addSdStudyNO16853
    tumor_size ~ add(addSd)
  })
}
