Ayyoub_2016_pyronaridine <- function() {
  description <- paste(
    "Pooled population PK model of oral pyronaridine in 349 pediatric malaria",
    "patients (0.51-15 years, 6.8-56.2 kg) from one phase II and five phase III",
    "studies of the pyronaridine-artesunate fixed-dose combination (Pyramax).",
    "Two-compartment disposition with first-order absorption and first-order",
    "elimination from the central compartment. Body weight enters as fixed",
    "allometric scaling (exponent 0.75 on CL/F and Q/F, 1.00 on V2/F and V3/F,",
    "centred on a 20 kg reference). Age enters as a power covariate on the",
    "peripheral volume V3/F (exponent 0.624, centred on a 7 yr reference).",
    "Formulation (1 = pediatric granule sachet, 0 = tablet) increases the",
    "absorption rate Ka by 1.63-fold over the tablet baseline. Residual error",
    "is additive on the natural-log concentration scale (equivalent to",
    "proportional in linear space). Dose is encoded as pyronaridine base in",
    "mg (paper Methods: pyronaridine tetraphosphate doses are multiplied by",
    "0.57 prior to modeling)."
  )
  reference <- paste(
    "Ayyoub A, Methaneethorn J, Ramharter M, Djimde AA, Tekete M, Duparc S,",
    "Borghini-Fuhrer I, Shin JS, Fleckenstein L. Population Pharmacokinetics of",
    "Pyronaridine in Pediatric Malaria Patients.",
    "Antimicrob Agents Chemother. 2016 Mar;60(3):1450-1458.",
    "doi:10.1128/AAC.02004-15"
  )
  vignette <- "Ayyoub_2016_pyronaridine"
  units    <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in the source analysis (paper does not state",
        "time-varying weight). Used for fixed-exponent allometric scaling",
        "centred on the 20 kg cohort-median reference: cl_typ = exp(lcl) *",
        "(WT/20)^0.75; vc_typ = exp(lvc) * (WT/20); q_typ = exp(lq) * (WT/20)^0.75;",
        "vp_typ = exp(lvp) * (WT/20) * (AGE/7)^e_age_vp. Allometric exponents",
        "(0.75 on CL/F and Q/F; 1.00 on V2/F and V3/F) were FIXED in the source",
        "analysis (Ayyoub 2016 Methods 'Base model development' p. 1453,",
        "Results 'Population pharmacokinetic model' p. 1454)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Power covariate on the peripheral volume",
        "V3/F centred on the 7 yr cohort-median reference: vp_typ = exp(lvp)",
        "* (WT/20) * (AGE/7)^0.624 (Ayyoub 2016 Results p. 1454 'V3/F (liters)",
        "= [3,230 * (weight/20) * (age/7)^0.624] * exp(eta)')."
      ),
      source_name        = "AGE"
    ),
    FORM_GRANULE = list(
      description        = "Pediatric granule-for-oral-suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = paste(
        "1 = subject received the pediatric granule (sachet) formulation of",
        "pyronaridine-artesunate, 0 = tablet (the reference). 67.0% of subjects",
        "received the tablet (Ayyoub 2016 Results p. 1454). Per-subject",
        "categorical indicator; in the source dataset every record carries the",
        "single per-subject formulation assignment (Pyramax granule sachet for",
        "<20 kg patients in studies SP-C-007-07 and SP-C-013-11; tablet for the",
        "20+ kg patients in studies SP-C-004-06 / SP-C-005-06 / SP-C-006-06,",
        "and either form in the phase II study SP-C-003-05). Linear",
        "(additive-multiplier) effect on Ka centred on the tablet reference:",
        "ka_typ = exp(lka) * (1 + e_form_granule_ka * FORM_GRANULE), with",
        "e_form_granule_ka = 1.63 (Ayyoub 2016 Results p. 1454 'Ka (hours) =",
        "[17.9 * (1 + formulation * 1.63)] * exp(eta)'; note the published",
        "headline 'Ka (hours)' is a typesetting error -- the structural value",
        "is 17.9 day^-1, matching Table 2 'K_a (day^-1) 17.9')."
      ),
      source_name        = "FORM"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 349L,
    n_observations  = 1085L,
    n_studies       = 6L,
    age_range       = "0.51-15 years (median 7)",
    age_median      = "7 years",
    weight_range    = "6.8-56.2 kg (median 20)",
    weight_median   = "20 kg",
    sex_female_pct  = 52.4,
    race_ethnicity  = "Predominantly African (Gabon, Mali, and other African sites in the phase II / III trials); race not used as a model covariate.",
    disease_state   = "Acute uncomplicated Plasmodium falciparum or Plasmodium vivax malaria.",
    dose_range      = paste(
      "Pyronaridine-artesunate (PA) 6:2 mg/kg to 12:4 mg/kg once daily for 3",
      "days. Tablet PA 180:60 mg dosed by weight band (1-4 tablets across",
      "20-90 kg); pediatric granule PA 60:20 mg per sachet dosed by weight",
      "band (1-3 sachets across 5-25 kg). Pyronaridine tetraphosphate doses",
      "were converted to pyronaridine base by multiplying by 0.57 prior to",
      "modeling (Ayyoub 2016 Methods 'Population pharmacokinetic and",
      "statistical analyses' p. 1452)."
    ),
    regions         = "Sub-Saharan Africa (Gabon and Mali phase II / III sites) and Southeast Asia (phase III studies).",
    n_studies_detail = "Phase II SP-C-003-05 (tablet + granule, n=57); phase III SP-C-004-06 (tablet, n=40); SP-C-005-06 (tablet, n=143); SP-C-006-06 (tablet, n=9); SP-C-007-07 (granule, n=83); SP-C-013-11 (granule, n=17). Two SP-C-013-11 subjects vomited and were excluded.",
    formulation_split_pct = "67.0% tablet, 33.0% granule (Ayyoub 2016 Results p. 1454).",
    sampling        = "Phase II: rich sampling, predose and 0.5, 1, 1.5, 2.5, 4, 8, 12 h plus 3, 7, 14, 21 days after the first dose. Phase III: sparse sampling, one to two samples per subject across the day-0-to-day-3 and day-4-to-day-42 windows.",
    bloq_handling   = "165/1252 raw observations (13.2%) below the 5.7 ng/mL assay LOQ were excluded from the analysis; 3 observations (0.24%) were excluded as outliers (Ayyoub 2016 Results p. 1454, Table 1).",
    notes           = "Demographics summarised from Ayyoub 2016 Table 1 (pooled row across 6 studies)."
  )

  ini({
    # Structural parameters (Ayyoub 2016 Table 2, 'Estimate' column; reference
    # individual: 20 kg body weight, 7 yr age, tablet formulation).
    lka <- log(17.9)
    label("Apparent first-order absorption rate constant Ka at FORM_GRANULE = 0 (1/day)")  # Ayyoub 2016 Table 2: Ka = 17.9 day^-1 (%RSE 11.7)
    lcl <- log(377)
    label("Apparent oral clearance CL/F at 20 kg reference (L/day)")  # Ayyoub 2016 Table 2: CL/F = 377 L/day (%RSE 6.58)
    lvc <- log(2230)
    label("Apparent central volume of distribution V2/F at 20 kg reference (L)")  # Ayyoub 2016 Table 2: V2/F = 2230 L (%RSE 6.59)
    lvp <- log(3230)
    label("Apparent peripheral volume of distribution V3/F at 20 kg and 7 yr reference (L)")  # Ayyoub 2016 Table 2: V3/F = 3230 L (%RSE 15.0)
    lq  <- log(804)
    label("Apparent intercompartmental clearance Q/F at 20 kg reference (L/day)")  # Ayyoub 2016 Table 2: Q/F = 804 L/day (%RSE 11.2)

    # Fixed-exponent allometric scaling on body weight (Ayyoub 2016 Methods
    # 'Base model development' p. 1453: exponents fixed at 0.75 on clearance
    # parameters and 1.0 on volume parameters; reference weight 20 kg).
    e_wt_cl_q  <- fixed(0.75)
    label("Fixed allometric WT exponent shared across CL/F and Q/F (unitless)")  # Ayyoub 2016 Methods p. 1453: fixed at the theoretical 0.75
    e_wt_vc_vp <- fixed(1.00)
    label("Fixed allometric WT exponent shared across V2/F and V3/F (unitless)")  # Ayyoub 2016 Methods p. 1453: fixed at the theoretical 1.00

    # Covariate effects retained after backward elimination (P < 0.001).
    # Age on V3/F: power covariate centred on the 7 yr cohort median.
    e_age_vp <- 0.624
    label("Power-covariate exponent of AGE on V3/F (unitless; reference age 7 yr)")  # Ayyoub 2016 Table 2: theta_6 = 0.624 (%RSE 38.6); Results p. 1454
    # Formulation on Ka: linear additive-multiplier effect centred on tablet.
    e_form_granule_ka <- 1.63
    label("Additive-multiplier effect of pediatric granule formulation on Ka (unitless)")  # Ayyoub 2016 Table 2: theta_7 = 1.63 (%RSE 37.8); Results p. 1454

    # Inter-individual variability. Ayyoub 2016 Methods (p. 1453) writes
    # Pi = P_pop * exp(eta_i) with eta ~ N(0, omega^2); Table 2 reports
    # omega^2 in the 'Estimate' column and the corresponding %CV computed
    # as 100 * sqrt(omega^2). A full variance-covariance matrix could not
    # be fit successfully (Methods p. 1453), so a diagonal omega is used.
    # IIV on Q/F was fixed to zero in the source (Results p. 1454: 'The
    # IIV on Q/F was fixed to zero, as the estimates could not be obtained
    # with good precision'), so etalq is absent.
    etalcl ~ 0.166  # Ayyoub 2016 Table 2: omega^2 on CL/F = 0.166 (%RSE 26.7); CV 40.7%
    etalvc ~ 0.993  # Ayyoub 2016 Table 2: omega^2 on V2/F = 0.993 (%RSE 8.76); CV 99.6%
    etalvp ~ 0.256  # Ayyoub 2016 Table 2: omega^2 on V3/F = 0.256 (%RSE 48.4); CV 50.6%
    etalka ~ 0.433  # Ayyoub 2016 Table 2: omega^2 on Ka  = 0.433 (%RSE 26.6); CV 65.8%

    # Residual variability. The source paper modelled natural-log
    # concentrations with additive residual on the log scale (Methods
    # p. 1453: 'ln Cij = ln Cpred,ij + eps_ij' with eps ~ N(0, sigma^2)).
    # By the standing convention in references/parameter-names.md, NONMEM
    # additive-on-log-scale residual maps to nlmixr2 proportional residual
    # in linear space, with propSd = sqrt(sigma^2) on the log scale
    # corresponding to ~CV in linear space to first order. Table 2 reports
    # sigma^2 = 0.195 (%RSE 12.1); propSd = sqrt(0.195) = 0.4416.
    propSd <- 0.4416
    label("Proportional residual SD for pyronaridine blood concentration")  # Ayyoub 2016 Table 2: RV (additive error) variance = 0.195 (%RSE 12.1); propSd = sqrt(0.195)
  })

  model({
    # Individual PK parameters with allometric and covariate effects.
    # Reference individual: WT = 20 kg, AGE = 7 yr, FORM_GRANULE = 0 (tablet).
    # Ka is shifted by a tablet-vs-granule additive-multiplier; V3/F carries
    # an age-on-volume power effect; CL, V2, Q, V3 all carry the fixed
    # allometric (WT/20) effect (exponent 0.75 on clearances, 1.0 on volumes).
    ka <- exp(lka + etalka) * (1 + e_form_granule_ka * FORM_GRANULE)
    cl <- exp(lcl + etalcl) * (WT / 20)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 20)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 20)^e_wt_vc_vp * (AGE / 7)^e_age_vp
    q  <- exp(lq)           * (WT / 20)^e_wt_cl_q

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system with first-order oral absorption from
    # depot into central, first-order distribution central <-> peripheral1,
    # and first-order elimination from central. Dose is administered to
    # depot as pyronaridine base in mg (the user/vignette is responsible
    # for the tetraphosphate -> base 0.57 conversion when applicable).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Blood pyronaridine concentration. Dose mg, volume L -> mg/L = ug/mL;
    # multiply by 1000 in the vignette to compare against the paper's
    # reported ng/mL values (assay LOQ 5.7 ng/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
