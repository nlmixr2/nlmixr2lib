Pettersen_2009_pantoprazole <- function() {
  description <- "Two-compartment population PK model for intravenous pantoprazole in 20 paediatric intensive-care patients aged 10 days to 16.4 years (Pettersen 2009). Pantoprazole is given as a zero-order infusion (15-30 min) into the central compartment with first-order elimination. Body-weight allometric scaling is fixed (0.75 on CL/Q, 1 on Vc/V2, reference 20 kg). Clearance is further modified by age (power on AGE/5 years), and three binary clinical covariates retained at the forward-selection / backward-elimination step: systemic inflammatory response syndrome (DIS_SIRS), concomitant CYP2C19-inhibitor coadministration (CONMED_CYP2C19_INH, pooling fluconazole, voriconazole, and isoniazid), and clinically defined hepatic dysfunction (HEPIMP, paediatric criterion TBILI >= 4 mg/dL OR ALT > 2x ULN for age). Each of the three indicators reduces pantoprazole CL by 62.3%, 65.8%, and 50.5% respectively when present alone. The reference subject is a 20 kg / 5-year-old paediatric ICU patient without SIRS, hepatic dysfunction, or CYP2C19 inhibitor coadministration."
  reference   <- paste(
    "Pettersen G, Mouksassi M-S, Theoret Y, Labbe L, Faure C,",
    "Nguyen B, Litalien C.",
    "Population pharmacokinetics of intravenous pantoprazole in paediatric",
    "intensive care patients.",
    "Br J Clin Pharmacol. 2009;67(2):216-227.",
    "doi:10.1111/j.1365-2125.2008.03328.x.",
    sep = " "
  )
  vignette    <- "Pettersen_2009_pantoprazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling fixed a priori at 0.75 on clearances (CL, Q)",
        "and 1 on volumes (Vc, V2), standardised to a body weight of 20 kg",
        "(Pettersen 2009 Methods, citing common paediatric practice).",
        "Cohort range 2.7-84.5 kg, median 12.7 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pettersen 2009 final-model CL equation includes a power effect of",
        "(AGE / 5 years) with exponent 0.316; the AGE / 5 ratio is",
        "dimensionless and unaffected by the units of AGE provided AGE and",
        "the reference 5 carry the same unit. Cohort range 10 days to 16.4",
        "years, median 2.1 years (Table 1). The age term is structurally",
        "complementary to the body-weight allometric term and captures the",
        "paediatric CYP2C19 ontogeny that is not fully explained by size."
      ),
      source_name        = "AGE"
    ),
    DIS_SIRS = list(
      description        = "Systemic inflammatory response syndrome (SIRS) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no SIRS)",
      notes              = paste(
        "Pettersen 2009 follows the consensus paediatric SIRS criteria of",
        "the ACCP / SCCM (Bone 1992; reference 39 in the paper): >= 2 of",
        "(1) temperature > 38 C or < 36 C, (2) heart rate above the",
        "age-specific threshold, (3) respiratory rate or PaCO2 above",
        "age-specific threshold, (4) WBC > 12000 or < 4000 cells/mm3 or",
        ">= 10% bands. 7 of 20 patients met SIRS criteria at study entry.",
        "Treated as time-fixed at the start of the PK observation window;",
        "future use of this covariate over a time-varying ICU course should",
        "set DIS_SIRS on each day SIRS criteria are met."
      ),
      source_name        = "SIRS"
    ),
    CONMED_CYP2C19_INH = list(
      description        = "Concomitant CYP2C19-inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP2C19 inhibitor)",
      notes              = paste(
        "Pettersen 2009 pools fluconazole, voriconazole, and isoniazid into",
        "the inhibitor-positive group. 4 of 20 patients positive (all from",
        "cohort I); patient #4 received both fluconazole and the inducer",
        "rifampicin -- the inducer covariate was not retained in the final",
        "model (only 1 of 20 patients positive, no significant OFV change).",
        "Time-fixed indicator over the PK observation window; future models",
        "may treat it as time-varying when concomitant medication start /",
        "stop dates are resolved at the per-observation level."
      ),
      source_name        = "INH"
    ),
    HEPIMP = list(
      description        = "Clinical hepatic-dysfunction indicator (paediatric criterion)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hepatic dysfunction)",
      notes              = paste(
        "Pettersen 2009 defines hepatic dysfunction as total bilirubin",
        ">= 4 mg/dL OR ALT > 2x upper limit of normal for age (Methods,",
        "citing reference 39). This is the paediatric-clinical cut-point",
        "scheme of the HEPIMP canonical, distinct from the adult NCI ODWG",
        "and Child-Pugh schemes; per-model documentation of the criterion",
        "is required and recorded here. 5 of 20 patients positive (4 from",
        "cohort I, 1 from cohort II)."
      ),
      source_name        = "HEP"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 20L,
    n_studies        = 1L,
    n_observations   = 156L,
    age_range        = "10 days to 16.4 years",
    age_median       = "2.1 years (cohort I 9.4, cohort II 0.7)",
    weight_range     = "2.7-84.5 kg",
    weight_median    = "12.7 kg (cohort I 30.6, cohort II 6.8)",
    sex_female_pct   = 35,
    race_ethnicity   = "Not reported (single-centre Montreal ICU)",
    disease_state    = paste(
      "Paediatric intensive-care patients at risk for or with upper",
      "gastrointestinal bleeding. Cohort I (n = 8): retrospective",
      "convenience sample of patients prescribed intravenous pantoprazole",
      "in 2002 after dose-finding requests from attending physicians;",
      "indications were refractory epigastric pain (5) and active upper",
      "GI bleeding (3); underlying diseases mostly hepatic disease or",
      "transplantation. Cohort II (n = 12): single-centre open-label phase",
      "I / II study (started February 2004) of stress-ulcer prophylaxis;",
      "indications were respiratory failure, coagulopathy, or paediatric",
      "Risk of Mortality score >= 10. SIRS present in 7 / 20; hepatic",
      "dysfunction (TBILI >= 4 mg/dL OR ALT > 2x ULN for age) in 5 / 20;",
      "concomitant CYP2C19 inhibitor in 4 / 20 (all cohort I); concomitant",
      "CYP2C19 inducer in 1 / 20 (not retained in the final model)."
    ),
    dose_range       = "19.9-140.6 mg/1.73 m^2/day intravenously over 15-30 min, once daily (one patient twice daily). Initial regimen 20 mg/1.73 m^2/day in neonates and 40 mg/1.73 m^2/day above 1 month.",
    regions          = "Canada (Centre Hospitalier Universitaire Sainte-Justine, Montreal)",
    nonmem_method    = "NONMEM VI (Level 1.2), FOCE-I",
    sampling_schema  = paste(
      "Cohort I: predose then 0, 0.25, 0.75, 1, 2, 4, 6, 12 h after the end",
      "of pantoprazole infusion. Cohort II: predose then 0, 0.25, 0.5, 1,",
      "2, 4, 8, 12, 24 h after the end of infusion. Plasma assay HPLC-UV at",
      "290 nm; LLOQ 0.1 mg/L, ULOQ 25 mg/L, within- and between-run CV < 10%."
    ),
    notes            = paste(
      "Demographics from Pettersen 2009 Table 1. Underlying conditions",
      "summarised: cohort I (n = 8) included 4 hepatic disease /",
      "transplantation, 2 haematological disorders, 1 respiratory failure,",
      "1 polytrauma; cohort II (n = 12) included 10 open heart surgery, 1",
      "respiratory failure, 1 shock. Sex-female fraction (35%) derived from",
      "Results paragraph 1: 'Twenty patients (13 boys, 7 girls)'.",
      "Predictive-check P-value 0.52 (Figure 2); bootstrap 927 / 1000",
      "successful runs (Table 3 footnote)."
    )
  )

  ini({
    # ===== Structural PK (Pettersen 2009 Table 3 'Final model' column) =====
    # Reference subject: WT = 20 kg, AGE = 5 years, no SIRS, no hepatic
    # dysfunction, no CYP2C19 inhibitor (Pettersen 2009 Table 3 footnote).
    lcl <- log(5.28); label("Typical CL at reference (L/h)")              # Pettersen 2009 Table 3: CL = 5.28 L/h, RSE 10.9% (bootstrap median 5.08, 95% CI 3.88-6.90)
    lvc <- log(2.22); label("Typical central volume Vc at reference (L)") # Pettersen 2009 Table 3: Vc = 2.22 L,  RSE 12.3% (bootstrap median 2.20, 95% CI 1.54-2.83)
    lq  <- log(1.10); label("Typical inter-compartmental clearance Q at reference (L/h)") # Pettersen 2009 Table 3: Q = 1.1 L/h,  RSE 19.0% (bootstrap median 1.1, 95% CI 0.7-1.6)
    lvp <- log(2.73); label("Typical peripheral volume V2 at reference (L)")              # Pettersen 2009 Table 3: V2 = 2.73 L, RSE 25.3% (bootstrap median 2.69, 95% CI 1.76-6.04)

    # ===== Allometric exponents (Pettersen 2009 Methods, fixed) =====
    # "The power factor for body weight was fixed at 0.75 for clearances
    # and 1 for volumes, as is common practice in paediatric studies."
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/20) for CL (unitless, fixed)") # Pettersen 2009 Methods, Population PK analysis paragraph 4
    e_wt_vc <- fixed(1.00); label("Allometric exponent on (WT/20) for Vc (unitless, fixed)") # Pettersen 2009 Methods, Population PK analysis paragraph 4
    e_wt_q  <- fixed(0.75); label("Allometric exponent on (WT/20) for Q (unitless, fixed)")  # Pettersen 2009 Methods, Population PK analysis paragraph 4
    e_wt_vp <- fixed(1.00); label("Allometric exponent on (WT/20) for V2 (unitless, fixed)") # Pettersen 2009 Methods, Population PK analysis paragraph 4

    # ===== Covariate effects on CL (Pettersen 2009 Table 3, final-model equation) =====
    # CL = 5.28 * (WT/20)^0.75 * 0.377^SIRS * (AGE/5)^0.316 *
    #              0.342^INH * 0.495^HEP
    # The dichotomous covariates enter as a multiplicative factor raised to
    # the binary indicator (= the factor when indicator = 1, = 1 when 0).
    # AGE enters as a power on (AGE/5 years).
    e_age_cl                <- 0.316; label("Power exponent on (AGE / 5 years) for CL (unitless)")                                # Pettersen 2009 Table 3: Age covariate effect = 0.316, RSE 12.4%
    e_dis_sirs_cl           <- 0.377; label("Multiplicative factor on CL for SIRS positive (unitless)")                            # Pettersen 2009 Table 3: SIRS covariate effect = 0.377, RSE 28.5%; SIRS decreases CL by 62.3%
    e_conmed_cyp2c19_inh_cl <- 0.342; label("Multiplicative factor on CL for concomitant CYP2C19 inhibitor (unitless)")            # Pettersen 2009 Table 3: CYP2C19 inhibitor covariate effect = 0.342, RSE 37.1%; INH decreases CL by 65.8%
    e_hepimp_cl             <- 0.495; label("Multiplicative factor on CL for clinical hepatic dysfunction (unitless)")             # Pettersen 2009 Table 3: Hepatic dysfunction covariate effect = 0.495, RSE 20.9%; HEP decreases CL by 50.5%

    # ===== Inter-individual variability (Pettersen 2009 Table 3 final model) =====
    # Exponential IIV reported as approximate CV% (square root of variance,
    # Table 3 footnote). Conversion to log-normal variance for nlmixr2:
    # omega^2 = log(1 + CV^2). Independent etas (no correlation block in
    # the final model; full-block correlations were tested per Methods but
    # not retained).
    etalcl ~ 0.10043  # Pettersen 2009 Table 3: IIV CL = 32.5%; log(1 + 0.325^2) = 0.10043
    etalvc ~ 0.15267  # Pettersen 2009 Table 3: IIV Vc = 40.6%; log(1 + 0.406^2) = 0.15267
    etalq  ~ 0.05924  # Pettersen 2009 Table 3: IIV Q  = 24.7%; log(1 + 0.247^2) = 0.05924
    etalvp ~ 0.67413  # Pettersen 2009 Table 3: IIV V2 = 98.1%; log(1 + 0.981^2) = 0.67413

    # ===== Residual error (Pettersen 2009 Table 3 final model) =====
    # Combined additive + proportional. The additive SD was fixed in the
    # source NONMEM run at 1e-5 mg/L (Table 3 footnote: "The additive error
    # was fixed in the model"); the proportional SD was estimated at 19.5%
    # CV.
    addSd  <- fixed(0.00001); label("Additive residual SD (mg/L, fixed)")     # Pettersen 2009 Table 3: additive SD = 0.00001 mg/L, FIXED
    propSd <- 0.195;          label("Proportional residual SD (fraction)")    # Pettersen 2009 Table 3: proportional SD = 19.5%, RSE 22.5%
  })

  model({
    # ----- Individual PK parameters -----
    # CL: allometric on WT + power on AGE + multiplicative factors for
    # SIRS, CYP2C19 inhibitor, hepatic dysfunction.
    cl <- exp(lcl + etalcl) * (WT / 20)^e_wt_cl *
          (AGE / 5)^e_age_cl *
          e_dis_sirs_cl^DIS_SIRS *
          e_conmed_cyp2c19_inh_cl^CONMED_CYP2C19_INH *
          e_hepimp_cl^HEPIMP
    vc <- exp(lvc + etalvc) * (WT / 20)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 20)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 20)^e_wt_vp

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # Intravenous pantoprazole zero-order infusion (15-30 min per Methods);
    # dose enters the central compartment directly.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    # Plasma pantoprazole concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
