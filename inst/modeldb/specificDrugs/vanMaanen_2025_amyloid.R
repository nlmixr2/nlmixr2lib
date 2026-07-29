vanMaanen_2025_amyloid <- function() {
  description <- paste(
    "Population exposure-response model of brain amyloid plaque burden",
    "(Centiloid, CL) in patients with mild cognitive impairment or mild",
    "dementia due to Alzheimer's disease. A single indirect-response",
    "(turnover) ODE describes the natural plaque time course (zero-order",
    "formation Kin, first-order elimination Kout). The BACE1 inhibitor",
    "verubecestat inhibits plaque formation via a fixed-Imax AUC-driven",
    "sigmoid on Kin (Inh_verub, AUC50 in uM*h of verubecestat AUC over",
    "24 h at steady state; Imax fixed to 1 per Table 1). Each of four",
    "anti-A-beta monoclonal antibodies (aducanumab, donanemab,",
    "gantenerumab, lecanemab) stimulates plaque elimination via a linear",
    "AUC/molecular-weight term on Kout (Stim_mAb = slope_mAb * (AUC_mAb",
    "/ MW_mAb), with distinct slopes per mAb). AUC of each drug is a",
    "time-varying model input (covariate); the paper's simulations set",
    "each AUC to its steady-state value during dosing periods and to 0",
    "during washout. Baseline plaque burden is a per-subject covariate",
    "used as the plaque initial condition. Fit by NONMEM FOCE-I to 370",
    "individual verubecestat amyloid PET measurements from 188 aMCI due",
    "to AD subjects (phase 3 APECS, NCT01953601) pooled with 120",
    "summary-level PET measurements from four anti-A-beta mAb trial",
    "programmes; external validation on 1506 amyloid-PET measurements",
    "from 521 amyloid-positive ADNI subjects with up to 10 years of",
    "follow-up. No IIV is encoded (paper: 'no additional",
    "interindividual variability was included'; summary-level mAb data",
    "constrained the model to central-tendency predictions).",
    sep = " "
  )
  reference <- paste(
    "van Maanen E, Robey S, Bennacef I, Duffull S, Egan MF, Kennedy ME,",
    "Stone JA (2025).",
    "Modeling amyloid plaque turnover dynamics improves characterization",
    "of drug effects.",
    "Alzheimer's & Dementia: Translational Research & Clinical",
    "Interventions 11(3):e70169.",
    "doi:10.1002/trc2.70169.",
    sep = " "
  )
  vignette <- "vanMaanen_2025_amyloid"

  paper_specific_compartments <- c("plaque")

  units <- list(
    time          = "day",
    dosing        = paste(
      "(this model has no dosing events; drug exposure enters via the",
      "time-varying covariates AUC_VERUB, AUC_ADU, AUC_DON, AUC_GAN,",
      "AUC_LEC. See covariateData for units and reference values.)",
      sep = " "
    ),
    concentration = paste(
      "Amyloid plaque burden 'plaque' in Centiloid units (CL);",
      "additive residual error on plaque.",
      sep = " "
    )
  )

  covariateData <- list(
    PLAQUE_BL = list(
      description        = paste(
        "Per-subject baseline amyloid plaque burden on the Centiloid",
        "scale used as the initial condition for the plaque state at",
        "t = 0 (before any drug exposure). van Maanen 2025 Supplement",
        "Section 'Further details on the exposure-response model':",
        "'Individual participant \"observed\" baseline amyloid plaque",
        "burden (before the first dose) was an initial condition at",
        "time 0. No additional interindividual variability was",
        "included.' Reference baseline values in the modelling cohort",
        "ranged 70-104 CL across the pooled trials (Table S2) with 71.2",
        "CL in the APECS analysis dataset.",
        sep = " "
      ),
      units              = "CL (Centiloid; 0 = no plaque, 100 = typical AD plaque burden)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Required per-subject covariate. Set to the observed individual",
        "baseline value; for typical-subject simulations at the natural-",
        "progression starting condition, use e.g. 10 CL (paper Section 2.3",
        "starting condition for simulations); for treatment-onset",
        "simulations, use 85 or 100 CL (paper Section 2.3). Follows the",
        "HGB_BL / FERRITIN_BL pattern for baseline-as-initial-condition",
        "covariates.",
        sep = " "
      ),
      source_name        = "(paper does not report a specific NONMEM column name; individual APECS baselines came from the trial PET dataset)"
    ),
    AUC_VERUB = list(
      description        = paste(
        "Verubecestat plasma exposure (AUC over 24 h at steady state).",
        "Enters the model as the driver of the inhibitory sigmoid Emax",
        "on plaque formation Kin (Inh_verub = Imax * AUC_VERUB /",
        "(AUC_VERUB + AUC50); Eq 4 of the paper). Set to 0 when",
        "verubecestat is not on board; set to the observed / projected",
        "steady-state AUC (dose-adjusted from a upstream popPK model such",
        "as Dockendorf 2022, cited in Table S1) during verubecestat",
        "dosing periods. Typical AUC values from Table S5 (BACE1",
        "inhibition levels): 40 mg daily -> approximately 4.4 uM*h (91.8%",
        "inhibition); 11.9 mg -> approximately 1.44 uM*h (78.7%); 2.9 mg",
        "-> approximately 0.35 uM*h (47.2%); 1.7 mg -> approximately 0.21",
        "uM*h (34.5%). (Reverse-computed from the reported inhibition",
        "levels and AUC50 = 0.392 uM*h using Eq 4 rearranged.)",
        sep = " "
      ),
      units              = "uM*h (micromolar * hour; verubecestat AUC over 24 h at steady state)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying step-wise (held constant during a verubecestat",
        "treatment period, resets to 0 during washout). Must be in the",
        "SAME units as AUC50 (uM*h) so the sigmoid AUC_VERUB / (AUC_VERUB",
        "+ AUC50) is dimensionless. Member of the AUC_<DRUG> canonical",
        "family (compare AUC_CARBO, AUC_GEM, AUC_GCV, AUC_LCM, AUC_CBZ,",
        "AUC_PAZO, AUC_RTV, AUC_EMPA in inst/references/covariate-",
        "columns.md).",
        sep = " "
      ),
      source_name        = "AUC_verub (paper Eq 4)"
    ),
    AUC_ADU = list(
      description        = paste(
        "Aducanumab serum exposure (AUC over 4 weeks at steady state).",
        "Enters the model as the driver of the linear stimulation term",
        "on plaque elimination Kout for aducanumab (Stim_ADU =",
        "slope_adu * (AUC_ADU / MW_ADU); Eq 6). Set to 0 when",
        "aducanumab is not on board; set to the observed / projected",
        "steady-state AUC during aducanumab dosing periods. Reference",
        "regimen (Table S3): EMERGE/ENGAGE high-dose titration -> target",
        "10 mg/kg IV Q4W (aducanumab), giving a typical 4-week AUC of",
        "approximately 1944 mg*day/L for a 70 kg subject at CL",
        "approximately 0.36 L/day.",
        sep = " "
      ),
      units              = "mg*day/L (milligrams * day per litre; aducanumab AUC over 4 weeks at steady state). Divide by MW_ADU (145912 g/mol; hard-coded in model()) to obtain the paper's dimensionless AUC/MW term in mmol*day/L (= mM*day).",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying step-wise. Units MUST be mg*day/L (not mg*h/mL,",
        "not ug*h/mL) so that (AUC_ADU / MW_ADU_g_per_mol) evaluates in",
        "mmol*day/L (= mM*day) and slope_adu (mM^-1 * day^-1) *",
        "AUC/MW is dimensionless. Conversion from mg*h/L: divide by 24.",
        "Member of the AUC_<DRUG> canonical family.",
        sep = " "
      ),
      source_name        = "AUC_mAb for aducanumab (paper Eq 6)"
    ),
    AUC_DON = list(
      description        = paste(
        "Donanemab serum exposure (AUC over 4 weeks at steady state).",
        "Enters the model as the driver of the linear stimulation term",
        "on plaque elimination Kout for donanemab (Stim_DON =",
        "slope_don * (AUC_DON / MW_DON); Eq 6). Reference regimen",
        "(Table S3): TRAILBLAZER-ALZ 1/2/4 titration to 1400 mg IV Q4W",
        "after 3 loading doses of 10 mg/kg IV Q4W.",
        sep = " "
      ),
      units              = "mg*day/L (see AUC_ADU units notes)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying step-wise; same unit convention as AUC_ADU (mg*day/L).",
        "Divided in model() by MW_DON = 145087 g/mol.",
        sep = " "
      ),
      source_name        = "AUC_mAb for donanemab (paper Eq 6)"
    ),
    AUC_GAN = list(
      description        = paste(
        "Gantenerumab serum exposure (AUC over 4 weeks at steady state).",
        "Enters the model as the driver of the linear stimulation term",
        "on plaque elimination Kout for gantenerumab (Stim_GAN =",
        "slope_gan * (AUC_GAN / MW_GAN); Eq 6). Reference regimen",
        "(Table S3): GRADUATE 1/2 titration to 1020 mg every 4 weeks",
        "(delivered as 510 mg Q2W) after a 9-month step-up.",
        sep = " "
      ),
      units              = "mg*day/L (see AUC_ADU units notes)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying step-wise; same unit convention as AUC_ADU.",
        "Divided in model() by MW_GAN = 146300 g/mol.",
        sep = " "
      ),
      source_name        = "AUC_mAb for gantenerumab (paper Eq 6)"
    ),
    AUC_LEC = list(
      description        = paste(
        "Lecanemab serum exposure (AUC over 4 weeks at steady state).",
        "Enters the model as the driver of the linear stimulation term",
        "on plaque elimination Kout for lecanemab (Stim_LEC =",
        "slope_lec * (AUC_LEC / MW_LEC); Eq 6). Reference regimen",
        "(Table S3): Clarity AD 10 mg/kg IV Q2W (approximately",
        "20 mg/kg per 4 weeks).",
        sep = " "
      ),
      units              = "mg*day/L (see AUC_ADU units notes)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying step-wise; same unit convention as AUC_ADU.",
        "Divided in model() by MW_LEC = 150000 g/mol.",
        sep = " "
      ),
      source_name        = "AUC_mAb for lecanemab (paper Eq 6)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 188L,
    n_studies      = 8L,
    age_range      = "median age 70-75 years across contributing trials (Table S2); APECS individual dataset median 71.5 years",
    weight_range   = "(paper does not tabulate body-weight distribution; the modelled dataset covers roughly the age 55-85 body-mass range typical of late-life mild AD)",
    sex_female_pct = 51.5,
    race_ethnicity = c(White = 85.0, Black = 1.5, Asian = 8.0, Hispanic = 6.0),
    disease_state  = paste(
      "Mild cognitive impairment (aMCI, i.e., prodromal AD; APECS",
      "study of verubecestat) or mild dementia due to Alzheimer's",
      "disease (mAb studies). Baseline amyloid PET burden across the",
      "contributing cohorts spanned 70-104 CL (Table S2), consistent",
      "with amyloid-positive symptomatic AD. Approximately 60-70% APOE4",
      "carriers in the pooled trials; 53-72% APOE4 carriers depending",
      "on cohort.",
      sep = " "
    ),
    dose_range     = paste(
      "Model INPUTS are AUC covariates, not dose amounts, but the",
      "reference dosing regimens used to derive the input AUCs are:",
      "verubecestat 12-40 mg PO once daily (APECS phase 3);",
      "aducanumab 1-10 mg/kg IV Q4W with dose titration (EMERGE,",
      "ENGAGE, TRAILBLAZER-ALZ 4, PRIME); donanemab 700-1400 mg IV",
      "Q4W with 10 mg/kg loading (TRAILBLAZER-ALZ 1/2/4);",
      "gantenerumab 120-1020 mg SC / IV every 4 weeks with titration",
      "(SCarlet RoAD, Marguerite RoAD, GRADUATE 1/2); lecanemab",
      "10 mg/kg IV Q2W (Clarity AD, study 201). See Table S3 for the",
      "phase 3 simulation regimens.",
      sep = " "
    ),
    regions        = "predominantly North America and Europe; contributing trials also enrolled Asia-Pacific and Latin American sites (per Table S2 Hispanic and Asian race percentages up to 24% and 17% in some cohorts)",
    notes          = paste(
      "The n_subjects = 188 count reflects only the individual-level",
      "APECS verubecestat cohort; the mAb data (aducanumab, donanemab,",
      "gantenerumab, lecanemab) were pooled at the summary level from",
      "the 8 phase-1b / phase-2 / phase-3 trials listed in Table S1",
      "(PRIME, EMERGE, ENGAGE, TRAILBLAZER-ALZ 1/4, SCarlet RoAD +",
      "Marguerite RoAD OLEs, GRADUATE 1/2, Clarity AD / CTAD, Study",
      "201) contributing 120 summary-level plaque-time points to the",
      "fit. The n_studies = 8 count is the sum: 1 individual APECS +",
      "at least 8 summary-level trials across the four mAbs (aducanumab",
      "3 sources, donanemab 3 sources, gantenerumab 2 sources,",
      "lecanemab 2 sources per Table S1). External validation used a",
      "further 1506 individual PET measurements from 521 amyloid-",
      "positive ADNI subjects (up to 10 years of follow-up).",
      "Individual baseline plaque burden was carried as a per-subject",
      "initial condition (no IIV was estimated on any parameter).",
      "Baseline demographics were similar across the trials (Table S2),",
      "supporting the pooled meta-analysis; covariate effects on Kin,",
      "Kout, AUC50 or slope_mAb could not be estimated from the summary-",
      "level meta-data (Supplement 'Statistical methods' + 'Centiloid",
      "conversions to allow for cross-study pooling of data').",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural turnover parameters -- van Maanen 2025 Table 1 (right
    # column: 'Population estimate (%RSE)'). Kin and Kout were estimated
    # on the log scale in NONMEM (Table 1 rows: 'Log-transformed zero-
    # order formation rate ... K_in' and 'Log-transformed first-order
    # elimination rate ... K_out'). Direct back-transformation gives the
    # linear-scale values reported in the paper text (Section 3.1: 'K_in
    # was estimated to be 0.028 CL/day or 10 CL/year. K_out was
    # estimated to be 0.00029/day or 0.11/year, corresponding to a half-
    # life of amyloid plaque of ~6.4 years').
    # =====================================================================
    lkin  <- -3.58 ; label("Log zero-order amyloid plaque formation rate Kin (log-CL/day; exp(-3.58) = 0.0279 CL/day)")   # Table 1 Kin log-scale estimate -3.58 (RSE 2.4%)
    lkout <- -8.13 ; label("Log first-order amyloid plaque elimination rate Kout (log-1/day; exp(-8.13) = 0.000294/day; t1/2 = ln(2)/Kout = 6.46 years)") # Table 1 Kout log-scale estimate -8.13 (RSE 1.4%)

    # =====================================================================
    # Verubecestat Imax-sigmoid parameters on plaque formation Kin.
    # Imax is FIXED at 1 per Table 1 ('Maximal inhibition by verubecestat
    # (Imax) 1^a' with footnote 'a Fixed.'). AUC50 was estimated on the
    # linear scale in NONMEM at 0.392 uM*h (Table 1, RSE 59.7%); the
    # narrative gives AUC50 = 0.402 uM*h in Section 3.1 -- this appears to
    # be a rounding discrepancy between the Table 1 point estimate 0.392
    # and the derived value used in the 'reduce plaque formation by
    # 91.8%' calculation (91.8% = 4.4 / (4.4 + 0.402) with the target 40
    # mg AUC of approximately 4.4 uM*h). The Table 1 point estimate 0.392
    # is used here as the primary fitted value.
    # =====================================================================
    limax  <- fixed(log(1)) ; label("Log maximum verubecestat inhibition of plaque formation Imax (unitless; FIXED at 1 per Table 1)")  # Table 1 Imax = 1 (Fixed)
    lauc50 <- log(0.392)    ; label("Log verubecestat AUC producing 50% of maximal plaque-formation inhibition AUC50 (log-uM*h)")      # Table 1 AUC50 = 0.392 uM*h (RSE 59.7%)

    # =====================================================================
    # Anti-A-beta mAb linear stimulation slopes on plaque elimination
    # Kout -- van Maanen 2025 Table 1 (linear-scale point estimates). The
    # paper estimated these slopes on the linear scale in NONMEM (Table 1
    # RSEs 11.8-24.3% are quoted on linear scale); log-transform is
    # applied here for positivity constraint. Bootstrap 95% CIs (Table 1,
    # right-most column) are: aducanumab (573, 918), donanemab (761,
    # 1460), gantenerumab (249, 751), lecanemab (474, 780).
    # =====================================================================
    lslope_adu <- log(719)  ; label("Log linear slope on Kout for aducanumab (log-(mM^-1 * day^-1))")     # Table 1 slope_aducanumab = 719 mM^-1 day^-1 (RSE 11.9%)
    lslope_don <- log(1120) ; label("Log linear slope on Kout for donanemab (log-(mM^-1 * day^-1))")     # Table 1 slope_donanemab = 1120 mM^-1 day^-1 (RSE 11.8%)
    lslope_gan <- log(397)  ; label("Log linear slope on Kout for gantenerumab (log-(mM^-1 * day^-1))")  # Table 1 slope_gantenerumab = 397 mM^-1 day^-1 (RSE 24.3%)
    lslope_lec <- log(606)  ; label("Log linear slope on Kout for lecanemab (log-(mM^-1 * day^-1))")     # Table 1 slope_lecanemab = 606 mM^-1 day^-1 (RSE 12.9%)

    # =====================================================================
    # Additive residual error on plaque (CL). van Maanen 2025 Table 1
    # first row: 'Additive error, CL 6.39 (11.3%)'. Supplement Eq 3
    # (Statistical methods): y_ij = yhat_ij + eps_ij with eps ~ N(0,
    # sigma^2). Here sigma = 6.39 CL. No IIV was estimated on any
    # parameter (Supplement 'Further details on the exposure-response
    # model': 'No additional interindividual variability was included').
    # =====================================================================
    addSd <- 6.39 ; label("Additive residual SD on amyloid plaque (CL)")  # Table 1 additive error 6.39 CL (RSE 11.3%)
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Molecular weights (g/mol) -- hard-coded structural constants per
    #    Supplement Section 'Further details on the exposure-response
    #    model': aducanumab 145912; donanemab 145087; gantenerumab 146300;
    #    lecanemab 150000. These divide each mAb AUC input to convert
    #    mg*day/L to mmol*day/L (= mM*day), giving the AUC/MW term
    #    slope_mAb multiplies against.
    # ---------------------------------------------------------------------
    MW_ADU <- 145912
    MW_DON <- 145087
    MW_GAN <- 146300
    MW_LEC <- 150000

    # ---------------------------------------------------------------------
    # 2. Back-transformed structural parameters.
    # ---------------------------------------------------------------------
    kin        <- exp(lkin)
    kout       <- exp(lkout)
    imax       <- exp(limax)
    auc50      <- exp(lauc50)
    slope_adu  <- exp(lslope_adu)
    slope_don  <- exp(lslope_don)
    slope_gan  <- exp(lslope_gan)
    slope_lec  <- exp(lslope_lec)

    # ---------------------------------------------------------------------
    # 3. Drug-effect terms (Eq 4 and Eq 6). AUC_VERUB enters the sigmoid
    #    Imax; each mAb AUC enters a linear (AUC/MW)-driven stimulation.
    #    Each slope * (AUC / MW) evaluates in mmol*day/L * mM^-1 * day^-1
    #    = dimensionless as intended (units in covariateData notes).
    # ---------------------------------------------------------------------
    Inh_verub <- imax * AUC_VERUB / (AUC_VERUB + auc50)

    Stim_adu <- slope_adu * (AUC_ADU / MW_ADU)
    Stim_don <- slope_don * (AUC_DON / MW_DON)
    Stim_gan <- slope_gan * (AUC_GAN / MW_GAN)
    Stim_lec <- slope_lec * (AUC_LEC / MW_LEC)

    # The paper (Section 3.3, Figure 4B) simulates a single mAb at a
    # time -- the linear stimulation terms are additive on Kout so a
    # combined-mAb simulation is straightforward; treat overlap as an
    # extrapolation beyond the fit.
    Stim_mAb <- Stim_adu + Stim_don + Stim_gan + Stim_lec

    # ---------------------------------------------------------------------
    # 4. Plaque turnover ODE -- van Maanen 2025 Eq 5:
    #      d Plaque/dt = Kin * (1 - Inh_verub)
    #                    - Plaque * Kout * (1 + Stim_mAb)
    #    with Plaque(t=0) = Baseline_Plaque (per-subject covariate).
    # ---------------------------------------------------------------------
    d/dt(plaque) <- kin * (1 - Inh_verub) - plaque * kout * (1 + Stim_mAb)

    plaque(0) <- PLAQUE_BL

    # ---------------------------------------------------------------------
    # 5. Observation and residual error (Supplement Eq 3: additive
    #    error on plaque in CL units).
    # ---------------------------------------------------------------------
    plaque ~ add(addSd)
  })
}
