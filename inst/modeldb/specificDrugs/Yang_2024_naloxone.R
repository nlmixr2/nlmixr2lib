Yang_2024_naloxone <- function() {
  description <- paste(
    "Two-compartment population PK of intramuscular / subcutaneous",
    "naloxone delivered by auto-injector (NAI), with a depot plus three",
    "transit absorption compartments in series (all sharing the single",
    "transit rate constant KTR) and linear first-order elimination from",
    "the central compartment. Apparent clearance is allometrically",
    "scaled on total body weight (reference 70 kg); body weight was the",
    "only covariate retained after a full stepwise covariate search",
    "(age, sex, race and naloxone dose were all screened and dropped).",
    "Developed by Yang 2024 from 2063 naloxone concentrations in 48",
    "healthy adults across two crossover auto-injector studies covering",
    "0.4, 0.8, 2 and 10 mg. This is the naloxone PK layer that Yang 2024",
    "embedded, held fixed, inside each of its four mechanistic",
    "opioid-induced respiratory depression models",
    "(Yang_2024_naloxone_buprenorphine, Yang_2024_naloxone_morphine,",
    "Yang_2024_naloxone_fentanyl, Yang_2024_naloxone_carfentanil)."
  )
  reference <- paste(
    "Yang TE, Del Bene F, Lavezzi SM, Iavarone L, Zhang J, Kim J,",
    "Gjurich B, Kessler C. Mechanistic pharmacokinetic-pharmacodynamic",
    "modeling and simulations of naloxone auto-injector 10 mg reversal",
    "of opioid-induced respiratory depression.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1722-1733.",
    "doi:10.1002/psp4.13215. PMCID PMC11494827.",
    "Parameter source: Appendix S1 Table S4 (Parameter Estimates for",
    "Final Naloxone Population PK Model using NAI, run 45).",
    "Structure confirmed against the Appendix S1 example NONMEM control",
    "streams ($MODEL DEPOTN / TRAN1 / TRAN2 / TRAN3 / CENTRALN / PERIN",
    "and the matching $DES block), which encode the identical naloxone",
    "PK layer with every THETA held FIX at the Table S4 estimate."
  )
  vignette <- "Yang_2024_naloxone_opioid_reversal"
  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Appendix S1 control-stream
  # $MODEL block, which names each compartment explicitly, and against
  # the "; AMT in ug ... ; PK volumes = L" unit header of the same file.
  compartmentData <- list(
    depot       = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Only covariate retained in the final model (Appendix S1",
        "Table S3, run 45). Enters as an allometric power on apparent",
        "clearance only: CL/F = 3.26 * (WT/70)^0.538. The reference",
        "weight of 70 kg is explicit in every Appendix S1 control",
        "stream ('WT= 70 ; kg' followed by",
        "'NTVCL = THETA(2)*(WT/70)**THETA(3)'), so it is source-traced",
        "rather than assumed. Observed range in the analysis population",
        "was 57.2-100.2 kg (Table S1). A weight effect on central volume",
        "was significant in the univariate step (run 26, dOFV -5.39) and",
        "in the full model, but was removed in the backward step",
        "(run 45) and is therefore NOT in this model."
      ),
      source_name        = "Body Weight"
    )
  )

  # Screened during covariate model building (Appendix S1 Table S3) but
  # not retained in the final model, so they are documented rather than
  # encoded. Listing them here preserves the provenance of the paper's
  # covariate search without raising a "declared but not referenced"
  # convention warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on CL (run 20, dOFV -0.60) and on V2 (run 25, dOFV -0.22); neither reached the p<0.05 inclusion threshold of 3.84. Range 23-54 years (Table S1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL (run 24, dOFV -0.04), V2 (run 27, dOFV -2.96) and KTR (run 31, dOFV -5.90). Sex on KTR entered the full model but was removed in the first backward step (run 43, dOFV +6.04 < 10.828). 47.9 percent female (Table S2)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL (run 23, dOFV -2.21) and KTR (run 30, dOFV -3.07); neither reached the inclusion threshold. 75 percent of the analysis population was Black or African American (Table S2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48L,
    n_studies      = 2L,
    age_range      = "23-54 years (mean 37.9)",
    weight_range   = "57.2-100.2 kg (mean 78.2)",
    sex_female_pct = 47.9,
    race_ethnicity = c(White = 20.8, Black = 75.0, Other = 4.2),
    disease_state  = "Healthy adult volunteers",
    dose_range     = paste(
      "Naloxone auto-injector 0.4 mg, 0.8 mg (two 0.4 mg injections),",
      "2 mg and 10 mg intramuscular / subcutaneous. 24 subjects in a",
      "six-sequence three-period crossover received 0.4 / 0.8 / 2 mg;",
      "24 subjects in a two-sequence two-period crossover received",
      "2 mg and 10 mg."
    ),
    regions        = NA_character_,
    n_observations = 2063L,
    bmi_range      = "18.8-31.6 kg/m^2 (mean 26.5)",
    iov_structure  = paste(
      "Table S4 additionally reports inter-occasion variability on KTR",
      "(omega^2 = 0.127, 35.6 percent CV, shrinkage 13.9 / 27.9 / 45.2",
      "percent across the three occasions of the crossover design).",
      "IOV is NOT encoded structurally in this model file: the source",
      "does not define an operational occasion column for the",
      "model-library use case, and the nlmixr2lib convention (Brooks",
      "2021 / Andrews_2017_tacrolimus precedent) is to omit IOV when no",
      "occasion mapping is defined. Yang 2024's own simulation control",
      "streams likewise carry only the three IIV terms and no IOV.",
      "See the vignette Assumptions and deviations section."
    ),
    notes          = paste(
      "Model selection (Table S3): run 14 base model; weight on CL and",
      "on V2 plus sex on KTR formed the full model (run 40); backward",
      "elimination removed sex on KTR (run 43) then weight on V2",
      "(run 45, FINAL). A CL-V2 omega block (run 46) lowered OFV by",
      "17.7 but shrank the weight-on-CL exponent to 0.328 with ~40",
      "percent RSE, so run 45 was retained as final and the IIV terms",
      "here are diagonal. The Discussion notes CL/F 3.26 L/min versus a",
      "literature IV naloxone CL of 3.45 L/min, and apparent total",
      "Vd/F ~486 L (404 + 81.8) versus 114 L after IV, consistent with",
      "an IM bioavailability of roughly 35-36 percent; because F is not",
      "identifiable from IM data alone, all disposition parameters here",
      "are apparent (CL/F, V/F, Q/F)."
    )
  )

  ini({
    # All values are the final population estimates of Appendix S1
    # Table S4 (run 45). They are point estimates with reported RSEs, not
    # fixed constants, so they are NOT wrapped in fixed() here. (The
    # Appendix S1 simulation control streams re-declare the identical
    # numbers with a NONMEM FIX flag purely because those runs are
    # $SIM-only; that is a simulation device, not an estimation
    # constraint. The four coupled Yang_2024_naloxone_<opioid> models do
    # encode them as fixed(), because there they are inherited and
    # re-used without re-fitting.)
    lktr <- log(0.696)
    label("Transit absorption rate constant (KTR, 1/min)")           # Table S4: KTR 0.696 1/min, RSE 6.0%
    lcl <- log(3.26)
    label("Apparent clearance at 70 kg (CL/F, L/min)")               # Table S4: CL/F 3.26 L/min, RSE 2.5%
    lvc <- log(404)
    label("Apparent central volume of distribution (V2/F, L)")       # Table S4: V2/F 404 L, RSE 4.2%
    lq <- log(0.0847)
    label("Apparent intercompartmental clearance (Q/F, L/min)")      # Table S4: Q/F 0.0847 L/min, RSE 7.7%
    lvp <- log(81.8)
    label("Apparent peripheral volume of distribution (V3/F, L)")    # Table S4: V3/F 81.8 L, RSE 10.5%
    e_wt_cl <- 0.538
    label("Allometric exponent of body weight on CL/F (unitless)")   # Table S4: Body weight on CL/F 0.538, RSE 14.4%

    # IIV variances. Table S4 prints the omega^2 variance first and the
    # approximate CV% in parentheses, e.g. "0.0129 (11.4)" where
    # sqrt(0.0129) = 0.114. The same three variances appear verbatim in
    # the $OMEGA block of every Appendix S1 control stream, confirming
    # that the tabulated number is the variance and not an SD.
    etalktr ~ 0.111
    etalcl ~ 0.0129
    etalvc ~ 0.0658

    # Table S4: proportional residual error 0.157 (39.6 CV%). The
    # tabulated 0.157 is the variance; nlmixr2's propSd is an SD, so
    # propSd = sqrt(0.157) = 0.3962, matching the printed 39.6%.
    propSd <- 0.3962
    label("Proportional residual error SD (fraction)")               # Table S4: 0.157 variance = 39.6 CV%, RSE 10.3%
  })

  model({
    ktr <- exp(lktr + etalktr)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # Depot plus three transit compartments in series, every transfer
    # governed by the single rate constant KTR (Appendix S1 control
    # stream $DES: DADT(1) = -KTR*A(1); DADT(2..4) = KTR*A(n-1)-KTR*A(n);
    # DADT(5) receives KTR*A(4)). Mean absorption time is therefore
    # 4/KTR = 5.75 min at the typical KTR.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # Amounts are in ug and volumes in L, so central/vc is ug/L = ng/mL,
    # matching the control-stream scaling S5 = NV and its
    # "; DV: conc = ng/mL for PK" header.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
