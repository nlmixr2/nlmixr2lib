Alabdulkarim_2026_amikacin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous amikacin in",
    "paediatric patients from birth to 18 years (Alabdulkarim 2026), developed",
    "in neonates and refitted to the combined neonatal and older paediatric",
    "cohort. All four disposition parameters carry fixed allometric weight",
    "scaling to 70 kg (exponent 0.75 on CL and Q, 1 on V1 and V2). Clearance",
    "additionally carries a fixed Rhodin postmenstrual-age sigmoid (TM50 47.7",
    "weeks, Hill 3.4), a power effect of serum creatinine relative to the",
    "Ceriotti age-expected typical value, an exponential postnatal-age",
    "maturation from 44% of the mature value at birth, and an exponential",
    "reduction with concomitant glycopeptide use. Full 4x4 IIV block on CL,",
    "V1, Q and V2; combined proportional and additive residual error."
  )
  reference <- paste(
    "Alabdulkarim N, Barker C, Indrakumar B, Hatcher J, McGarrity O, Best T,",
    "Cheng IL, Standing JF. A unified population pharmacokinetic model of",
    "amikacin in paediatric patients from birth to 18 years for model-informed",
    "precision dosing. J Antimicrob Chemother. 2026;81(11):dkag354.",
    "doi:10.1093/jac/dkag354. Parameter values are the final-model estimates",
    "printed in the authors' nlmixr2 code (Supplementary Code S1), which agree",
    "with Table 2 to the reported precision. Age-expected serum creatinine per",
    "Ceriotti F, et al. Clin Chem. 2008;54(3):559-566",
    "(doi:10.1373/clinchem.2007.099648); postmenstrual-age maturation per",
    "Rhodin MM, et al. Pediatr Nephrol. 2009;24(1):67-76",
    "(doi:10.1007/s00467-008-0997-5)."
  )
  vignette <- "Alabdulkarim_2026_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL and Q as (WT / 70)^0.75 and V1 and V2 as (WT / 70)^1,",
        "both exponents fixed (Methods Eq. 1; Table 2 footnote; Code S1",
        "beta.wt.cl = fix(0.75), beta.wt.v = fix(1); the code applies",
        "beta.wt.cl to Q as well). Observed range 0.6-143 kg (Table 1). The",
        "70 kg reference is an adult anchor; the pooled median weight was",
        "4.73 kg."
      ),
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age plus postnatal age).",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the fixed Rhodin sigmoid PMA^3.4 / (47.7^3.4 + PMA^3.4) on",
        "CL (Methods Eq. 2). Supplied in WEEKS, as in the paper and the",
        "Rhodin constants; the register default for PAGE is months, so read",
        "these units before supplying the column. Observed range 24-992",
        "weeks (Table 1). Gestational age was missing for 51.8% of the",
        "cohort (17.2% of neonates); the authors then assumed 40 weeks",
        "gestation when computing PMA (Methods 'Data collection')."
      ),
      source_name = "PMA"
    ),
    PNA = list(
      description = "Postnatal (chronological) age since birth.",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper uses postnatal age in DAYS in the maturation function",
        "(Eq. 5) and in YEARS in the Ceriotti typical-creatinine polynomial",
        "(Eq. 3). The canonical PNA column is in months; model() recovers",
        "days as PNA * 30.4375 and years as PNA / 12 (1 year = 365.25",
        "days). PNA must be strictly positive: the Ceriotti polynomial",
        "contains ln(PNA_years), so at PNA = 0 the typical creatinine is",
        "infinite and clearance diverges. The youngest patient in the",
        "dataset was about 1 day old (0.003 years, Table 1); use PNA of at",
        "least 1 day (0.0329 months) when simulating."
      ),
      source_name = "PNA"
    ),
    CREAT = list(
      description = "Measured serum creatinine (enzymatic assay).",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL as (CREAT / TSCR)^e_creat_cl, where TSCR is the",
        "age-expected typical serum creatinine computed INSIDE the model",
        "from postnatal age with the Ceriotti 2008 second-order fractional",
        "polynomial printed as Methods Eq. 3:",
        "TSCR = -2.37330 - 12.91367 * ln(PNA_years) + 23.93581 * PNA_years^0.5",
        "(umol/L). No separate CREAT_REF column is needed. Observed range",
        "7-770 umol/L (Table 1). A patient whose creatinine equals the",
        "age-expected value has a renal factor of 1."
      ),
      source_name = "Scr"
    ),
    CONMED_GLYCOPEPTIDE = list(
      description = paste(
        "Concomitant glycopeptide antibiotic use (1 = yes, 0 = no)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant glycopeptide)",
      notes = paste(
        "Enters CL as exp(e_conmed_glycopeptide_cl * CONMED_GLYCOPEPTIDE)",
        "(Code S1: exp(... + beta.glycopeptide * glycopeptide)), i.e. an",
        "18.6% lower clearance with a glycopeptide. Table 2 prints the",
        "effect as (glycopeptides)^beta, which is undefined for a 0/1",
        "indicator; the code is what was fitted. The paper does not name",
        "the agents pooled into the indicator (in UK paediatric practice",
        "vancomycin and teicoplanin). 58 of 409 patients (14.2%) received",
        "a glycopeptide (Table 1)."
      ),
      source_name = "glycopeptide"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 409,
    n_studies = 1,
    n_observations = 1091,
    age_range = "birth to 18 years (0.003-18 years; median 0.12 years)",
    age_median = "0.12 years",
    weight_range = "0.6-143 kg",
    weight_median = "4.73 kg",
    pma_range = "24-992 weeks (median 44 weeks)",
    sex_female_pct = 44.3,
    race_ethnicity = NULL,
    renal_function = "serum creatinine 7-770 umol/L (median 41 umol/L)",
    co_medication = "concomitant glycopeptide in 58/409 (14.2%)",
    disease_state = paste(
      "hospitalised neonates, infants, children and adolescents receiving",
      "intravenous amikacin for at least 3 days, monitored by routine",
      "therapeutic drug monitoring"
    ),
    dose_range = paste(
      "5-30 mg/kg per dose (median 15 mg/kg), almost all once daily",
      "(97.8% q24h); doses >= 20 mg/kg as 30-min infusions, lower doses as",
      "~3-min slow bolus injections"
    ),
    regions = "United Kingdom (Great Ormond Street Hospital for Children, London)",
    notes = paste(
      "Retrospective TDM data, September 2020 to December 2022 (Methods;",
      "Table 1). Model built on 204 neonates (birth to < 1 month), evaluated",
      "without re-estimation on 205 patients aged 1 month to 18 years (72",
      "aged 1 month-2 years, 41 aged 3-6, 52 aged 7-12, 40 aged 13-18), then",
      "refitted to all 409 patients; this file holds the refitted (final)",
      "model. Median 2 samples per patient (range 2-11); 290 of 1091",
      "concentrations (26.6%) were below the 0.8 mg/L quantification limit",
      "and were handled with the M3 method. Male 228/409 (55.7%). Total",
      "serum amikacin by particle-enhanced turbidimetric immunoassay",
      "(assay range 0.8-100 mg/L). Estimation with nlmixr2 7.0.1 (FOCEI)."
    )
  )

  ini({
    # Structural parameters: typical values for a 70 kg, fully mature subject
    # (PMA factor 1, creatinine at the age-expected value, postnatal factor 1,
    # no glycopeptide). Values are the Code S1 final estimates.
    lcl <- 1.72851276237038
    label("Clearance CL at 70 kg, fully mature (L/h)")
    # Code S1 tcl = 1.72851276237038 (log scale; exp = 5.633 L/h); Table 2 'CL (L/h/70 kg)' = 5.634

    lvc <- 3.3789418101558
    label("Central volume V1 at 70 kg (L)")
    # Code S1 tv = 3.3789418101558 (log scale; exp = 29.34 L); Table 2 'V1 (L/70 kg)' = 29.34

    lq <- 0.0630953315121936
    label("Intercompartmental clearance Q at 70 kg (L/h)")
    # Code S1 tq = 0.0630953315121936 (log scale; exp = 1.065 L/h); Table 2 'Q (L/h/70 kg)' = 1.07

    lvp <- 3.40976623849171
    label("Peripheral volume V2 at 70 kg (L)")
    # Code S1 tvp = 3.40976623849171 (log scale; exp = 30.26 L); Table 2 'V2 (L/70 kg)' = 30.26

    # Allometric exponents, fixed a priori
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)")
    # Methods Eq. 1 and Table 2 footnote: fixed 0.75 for CL and Q; Code S1 beta.wt.cl = fix(0.75)

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of body weight on V1 and V2 (unitless)")
    # Methods Eq. 1 and Table 2 footnote: fixed 1.0 for V1 and V2; Code S1 beta.wt.v = fix(1)

    # Postmenstrual-age renal maturation (Rhodin 2009), fixed
    tmat50 <- fixed(47.7)
    label("Postmenstrual age at 50% of mature renal function, PMA50 (weeks)")
    # Methods 'Covariate model' and Eq. 2: PMA50 fixed to 47.7 weeks

    hill_mat <- fixed(3.4)
    label("Hill coefficient of the postmenstrual-age maturation (unitless)")
    # Methods 'Covariate model' and Eq. 2: Hill coefficient fixed to 3.4

    # Serum creatinine relative to the Ceriotti age-expected value
    e_creat_cl <- -0.72449
    label("Exponent of (CREAT / age-expected typical creatinine) on CL (unitless)")
    # Code S1 beta.cr.cl = -0.72449; Table 2 'beta Scr (CL)' = -0.73

    # Postnatal-age maturation, Eq. 5: hM + (1 - hM) * (1 - exp(-PNA_days * hN)),
    # rewritten exactly as 1 + cl_pna_famp * exp(-cl_pna_kdes * PNA_days)
    # with cl_pna_famp = hM - 1 and cl_pna_kdes = hN.
    cl_pna_famp <- -0.559126789740313
    label("Fractional deficit of CL at birth relative to the postnatally matured value, hM - 1 (unitless)")
    # Code S1 hM = 0.440873210259687, so hM - 1 = -0.559126789740313; Table 2 'hM (CL)' = 0.44

    lcl_pna_kdes <- log(0.2839214638671)
    label("Postnatal maturation rate constant of CL, hN (1/day of postnatal age)")
    # Code S1 hN = 0.2839214638671; Table 2 'hN (CL)' = 0.28

    e_conmed_glycopeptide_cl <- -0.20576
    label("Exponential effect of concomitant glycopeptide use on CL (unitless)")
    # Code S1 beta.glycopeptide = -0.20576; Table 2 'beta glycopeptides (CL)' = -0.21

    # IIV: full 4x4 block, Code S1 order eta.cl, eta.v, eta.q, eta.vp.
    # Diagonal CVs sqrt(exp(omega^2) - 1) = 49.93, 22.80, 91.40, 32.70 % match Table S5 'Final model' BSV rows;
    # Table 2 prints 49.94, 22.79, 91.50, 32.70 (Q differs by 0.1 CV points).
    etalcl + etalvc + etalq + etalvp ~ c(
      0.222597520852288,
      0.0186169355158569, 0.0506970970546489,
      0.229709539534124, 0.0211160875664682, 0.607283061694239,
      0.0334585338716289, 0.011002182666226, 0.0218988956792476, 0.101563602367693
    ) # Code S1 'eta.cl + eta.v + eta.q + eta.vp ~ c(...)' lower triangle

    # Residual error, combined
    propSd <- 0.24890419860156
    label("Proportional residual error (fraction)")
    # Code S1 prop.err = 0.24890419860156; Table 2 'Proportional error (%)' = 24.89

    addSd <- 1.01995057448704
    label("Additive residual error (mg/L)")
    # Code S1 add.err = 1.01995057448704; Table 2 'Additive error (mg/L)' = 1.020
  })

  model({
    # 1. Derived covariate terms
    pna_days <- PNA * 30.4375
    pna_years <- PNA / 12

    # Methods Eq. 2 (Rhodin 2009 sigmoid on postmenstrual age in weeks)
    f_pma <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # Methods Eq. 3 (Ceriotti 2008) and Eq. 4
    creat_typ <- -2.37330 - 12.91367 * log(pna_years) + 23.93581 * pna_years^0.5
    f_creat <- (CREAT / creat_typ)^e_creat_cl

    # Methods Eq. 5
    cl_pna_kdes <- exp(lcl_pna_kdes)
    f_pna <- 1 + cl_pna_famp * exp(-cl_pna_kdes * pna_days)

    f_glyco <- exp(e_conmed_glycopeptide_cl * CONMED_GLYCOPEPTIDE)

    # 2. Individual parameters (Table 2 equations; Code S1 model block)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * f_pma * f_creat * f_pna * f_glyco
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (intravenous dosing into central)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
