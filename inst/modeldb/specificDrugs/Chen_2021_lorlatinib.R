Chen_2021_lorlatinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral lorlatinib (a third-",
    "generation ALK/ROS1 tyrosine kinase inhibitor) in adult patients",
    "with advanced ALK-positive or ROS1-positive non-small cell lung",
    "cancer and healthy participants (Chen 2021; N = 425 across seven",
    "studies). Disposition is a two-compartment model with sequential",
    "zero-order and first-order oral absorption (dose enters the depot",
    "via a zero-order window of duration D1 = 1.15 h followed by",
    "first-order absorption at ka = 3.11 h^-1) and time-varying metabolic",
    "auto-induction of clearance: CL(t) = CLI + (CLMX - CLI) *",
    "(1 - exp(-cl_exp_kdes * t)), rising from a single-dose CLI = 9.04 L/h to",
    "a steady-state CLMX = 14.5 L/h with induction rate constant",
    "cl_exp_kdes = 0.020 h^-1 (~7.25 d to functional steady state; Chen 2021",
    "abstract, Table 4). CLI and CLMX share a fixed allometric",
    "exponent 0.75 on body weight (reference 70 kg) and both are",
    "modulated by a shared multiplicative covariate block:",
    "1 + e_alb_cl * (ALB - 40 g/L) with e_alb_cl = 0.00670 per g/L,",
    "1 + e_dose_lor_cl * (DOSE_LOR_MGD - 100 mg/day) with",
    "e_dose_lor_cl = 0.00100 per mg/day, and a power effect",
    "(CRCL / 100)^e_crcl_cl with e_crcl_cl = 0.235. V2 = 121 L carries",
    "a fixed allometric exponent 1.0 on body weight; V3 = 155 L and",
    "Q = 22.0 L/h have no allometric scaling. ka is modulated by",
    "proton-pump inhibitor co-administration: ka x (1 - 0.675 * CONMED_PPI),",
    "i.e. a 67.5% ka reduction on PPI (which reduces Cmax by ~30% with",
    "no effect on AUCinf per the Chen 2021 Discussion). Bioavailability",
    "F1 = 0.759. Inter-individual variability is a correlated CL/F",
    "block (log-scale variances 0.030 and 0.022; covariance -0.006),",
    "a correlated V2/V3 block (0.086, -0.017, 0.101), and an",
    "independent ka block (2.33). Residual error is a route-specific",
    "log-scale additive term (approximately proportional in linear",
    "space): propSd 11.5% CV for IV data (Study B7461007 absolute",
    "bioavailability arm) and 43.8% CV for oral data.",
    sep = " "
  )
  reference <- paste(
    "Chen J, Houk B, Pithavala YK, Ruiz-Garcia A. (2021).",
    "Population pharmacokinetic model with time-varying clearance for",
    "lorlatinib using pooled data from patients with non-small cell lung",
    "cancer and healthy participants.",
    "CPT Pharmacometrics Syst Pharmacol 10(2):148-160.",
    "doi:10.1002/psp4.12585",
    sep = " "
  )
  vignette <- "Chen_2021_lorlatinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Chen 2021. Cohort mean 70.53 kg",
        "(SD 16.89) across 425 subjects (Chen 2021 Table 3). Reference",
        "70 kg with fixed allometric exponents 0.75 on CL (both CLI and",
        "CLMX arms) and 1.0 on V2 (Chen 2021 Methods 'Model development'",
        "and NONMEM control stream ListS1: TVCLI = THETA(1)*(BWT/70)**0.75,",
        "TVV2 = THETA(2)*(BWT/70), TVCLMX = THETA(9)*(BWT/70)**0.75).",
        "V3, Q, ka, D1, F1 have no allometric scaling."
      ),
      source_name        = "BWT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Chen 2021. Retained as a significant",
        "covariate on CL in the final model. Enters as a linear centered",
        "multiplier on both CLI and CLMX: 1 + e_alb_cl * (ALB - 40 g/L)",
        "with e_alb_cl = 0.00670 per g/L (= 0.0670 per g/dL scaled to",
        "canonical SI units by a factor of 10; Chen 2021 Table 4",
        "theta_BALB_on_CL = 0.067 per g/dL centered at BALB = 4 g/dL).",
        "Cohort mean 3.92 g/dL (39.2 g/L) with SD 0.58 g/dL. Chen 2021",
        "reports a 5.3% CL reduction at the 10th-percentile ALB",
        "(3.2 g/dL = 32 g/L) and a 4% CL increase at the 90th-percentile",
        "ALB (4.6 g/dL = 46 g/L); the paper's Discussion attributes the",
        "trend to ALB as a marker of overall health status rather than a",
        "direct binding-driven effect."
      ),
      source_name        = "BALB (in g/dL; convert to canonical g/L by multiplying by 10)"
    ),
    DOSE_LOR_MGD = list(
      description        = "Current total daily lorlatinib dose",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record covariate; constant within an inter-dose",
        "interval and updated when the prescriber alters the daily dose.",
        "For q.d. regimens equals the single-dose amount (100 mg q.d.",
        "-> DOSE_LOR_MGD = 100); for b.i.d. regimens equals the sum",
        "across the day (75 mg b.i.d. -> DOSE_LOR_MGD = 150). Retained",
        "as a significant covariate on CL in the final model: enters as",
        "a linear centered multiplier 1 + e_dose_lor_cl * (DOSE_LOR_MGD",
        "- 100) shared between CLI and CLMX with e_dose_lor_cl =",
        "0.00100 per mg/day (Chen 2021 Table 4 theta_TDOSE_on_CL =",
        "0.001). Captures the paper's dose-nonlinearity finding that",
        "10 mg q.d. gives 12.4% lower CL and 200 mg q.d. gives 13.8%",
        "higher CL relative to the 100 mg q.d. reference (Chen 2021",
        "Results 'Covariate effects on lorlatinib steady-state",
        "exposure'), attributed by the Discussion to concentration-",
        "related increases in auto-induction potency at higher doses."
      ),
      source_name        = "TDOSE (total daily lorlatinib dose in mg)"
    ),
    CRCL = list(
      description        = "Baseline weight-normalized creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Chen 2021. Retained as a significant",
        "covariate on CL in the final model. Derived from the Cockcroft-",
        "Gault estimate BCCL by weight-normalization to remove the",
        "confounding of the a priori (BWT/70)^0.75 allometric CL scaling",
        "(Chen 2021 Methods 'Derivation of standardized creatinine",
        "clearance'; the source paper cites Rhodin et al. 2009). Enters",
        "as a power effect (CRCL / 100)^e_crcl_cl shared between CLI",
        "and CLMX with e_crcl_cl = 0.235 (Chen 2021 Table 4",
        "theta_WNCL_on_CL). Cohort mean 98.31 mL/min (SD 32.13) across",
        "the analysis dataset (Chen 2021 Table 3). For downstream users:",
        "in the reference 70-kg adult with BSA ~1.73 m^2, weight-",
        "normalized WNCL is numerically nearly identical to conventional",
        "BSA-normalized CRCL because BWT and BSA are highly correlated;",
        "the canonical CRCL column accepts either normalization form",
        "with the derivation documented per-model here."
      ),
      source_name        = "WNCL (weight-normalized creatinine clearance in mL/min)"
    ),
    CONMED_PPI = list(
      description        = "Concomitant proton-pump inhibitor use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no PPI co-administration)",
      notes              = paste(
        "Time-varying per dose record: 1 = subject received a PPI",
        "(rabeprazole in the founding study B7461008) concurrent with",
        "lorlatinib, 0 otherwise. Retained as a significant covariate",
        "on ka in the final model; PPI use had no significant effect on",
        "F or on downstream systemic exposure (AUCinf unchanged per the",
        "Chen 2021 Discussion). Enters as a linear multiplier on ka:",
        "ka * (1 + e_ppi_ka * CONMED_PPI) with e_ppi_ka = -0.675",
        "(Chen 2021 Table 4 theta_PPI_on_ka), i.e. a 67.5% ka reduction",
        "with PPI co-administration. Cohort prevalence 5% (23 of 425",
        "subjects; Chen 2021 Table 3 'No PPI use = 402 / 95%')."
      ),
      source_name        = "PPI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Baseline subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a candidate covariate on CL and V2 (Chen 2021",
        "Table 2 evaluated-covariates list) but not retained in the",
        "final model. Cohort mean 49.86 years (SD 13.20). Not required",
        "for simulation."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened as a candidate covariate on CL and V2 but not",
        "retained in the final model. Cohort 191 / 425 female (45%)",
        "per Chen 2021 Table 3."
      ),
      source_name        = "SEX (paper coded as SEX; 1-female / 0-male convention preserved)"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "Screened as a candidate covariate on CL and V2 but not",
        "retained in the final model. Cohort 113 / 425 Asian (27%) per",
        "Chen 2021 Table 3 (RACE: White 220 / 52%, Black 32 / 8%,",
        "Asian 113 / 27%, Other 29 / 7%)."
      ),
      source_name        = "RACE (Asian sub-level of the paper's RACE column)"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black)",
      notes              = paste(
        "Screened as a candidate covariate but not retained. Cohort",
        "32 / 425 Black (8%) per Chen 2021 Table 3."
      ),
      source_name        = "RACE (Black sub-level)"
    ),
    CYP3A5_EM = list(
      description        = "CYP3A5 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (poor / intermediate / ultra-rapid CYP3A5 metabolizer)",
      notes              = paste(
        "Screened but not retained in the final model. Cohort CYP3A5",
        "phenotype distribution per Chen 2021 Table 3: Poor 195 / 46%,",
        "Intermediate 66 / 16%, Extensive 17 / 4%, Ultra-rapid 0 / 0%.",
        "Because lorlatinib is metabolised primarily by CYP3A4 and",
        "UGT1A4 with lesser contribution from CYP3A5 (Chen 2021 Intro),",
        "a null CYP3A5 phenotype covariate is consistent with the",
        "mechanism."
      ),
      source_name        = "CYP3A5 (phenotype code encoded as ordinal 1-4)"
    ),
    CYP2C19_EM = list(
      description        = "CYP2C19 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (poor / intermediate / ultra-rapid CYP2C19 metabolizer)",
      notes              = paste(
        "Screened but not retained. Cohort CYP2C19 phenotype:",
        "Poor 18 / 4%, Intermediate 100 / 24%, Extensive 153 / 36%,",
        "Ultra-rapid 7 / 2% (Chen 2021 Table 3)."
      ),
      source_name        = "CYP2C19 (phenotype code encoded as ordinal 1-4)"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal food-effect indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Screened as a candidate covariate on ka and F but not",
        "retained. Study B7461008 also showed that administration of",
        "lorlatinib with a high-fat meal had no meaningful effect on",
        "lorlatinib exposure (Chen 2021 Intro paragraph 2)."
      ),
      source_name        = "FOOD"
    ),
    HEPIMP_MILD = list(
      description        = "Mild-hepatic-impairment indicator (NCI B1 or B2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = paste(
        "Screened but not retained in the final model. Cohort baseline",
        "hepatic impairment per Chen 2021 Table 3: Normal (A) 365 / 86%,",
        "Mild (B1) 50 / 12%, Mild (B2) 10 / 2%, Moderate-Severe (C-D)",
        "0 / 0%. The paper's Table S1 confirms no monotonic trend of",
        "CL with hepatic-impairment stage."
      ),
      source_name        = "BHGRADE (paper's ordinal 1-4 NCI classification)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 425L,
    n_studies       = 7L,
    n_observations  = 5806L,
    age_range       = "18-83 years (approximate; Chen 2021 Table 3 cohort mean 49.86 years, SD 13.20)",
    age_median      = "50 years",
    weight_range    = "35-136 kg (approximate; Chen 2021 Table 3 cohort mean 70.53 kg, SD 16.89)",
    weight_median   = "70 kg",
    sex_female_pct  = 45.0,
    race_ethnicity  = c(White = 52.0, Black = 8.0, Asian = 27.0, Other = 7.0),
    disease_state   = paste(
      "Pooled cohort of 330 patients with advanced anaplastic",
      "lymphoma kinase (ALK)-positive or c-ROS oncogene 1",
      "(ROS1)-positive non-small cell lung cancer (NSCLC) enrolled",
      "in the phase I/II dose-escalation and expansion study B7461001",
      "(NCT01970865), plus 95 healthy participants enrolled in six",
      "phase I clinical-pharmacology studies (B7461004 mass balance,",
      "B7461005 relative bioavailability, B7461007 absolute",
      "bioavailability, B7461008 rabeprazole/food effect, B7461011",
      "rifampin CYP3A4 induction, B7461016 bioequivalence)."
    ),
    dose_range      = paste(
      "Phase I patients: 10, 25, 50, 75, 100, 150, 200 mg orally q.d.,",
      "or 35, 75, or 100 mg orally b.i.d. Phase II patients: 100 mg",
      "orally q.d. (the labelled dose). Healthy volunteers: single",
      "oral 100 mg tablet (B7461005, B7461008, B7461011, B7461016) or",
      "single 50 mg intravenous vs 100 mg oral crossover (B7461007)."
    ),
    regions         = "Not summarised in Chen 2021 (multicentre phase I/II NSCLC study plus multi-site healthy-volunteer studies).",
    baseline_creatinine_clearance = "Cohort mean 98.31 mL/min (SD 32.13) per Chen 2021 Table 3 (Cockcroft-Gault estimate; renal impairment stages A normal 57%, B mild 31%, C moderate 12%, D severe 0.2%).",
    baseline_hepatic_function     = "Cohort: NCI A normal 86%, B1 mild 12%, B2 mild 2%, C-D moderate-severe 0% per Chen 2021 Table 3.",
    cyp3a5_phenotype              = "Poor 46%, Intermediate 16%, Extensive 4%, Ultra-rapid 0% per Chen 2021 Table 3.",
    cyp2c19_phenotype             = "Poor 4%, Intermediate 24%, Extensive 36%, Ultra-rapid 2% per Chen 2021 Table 3.",
    notes           = paste(
      "Chen 2021 Results 'Analysis dataset' and Table 1 summarise the",
      "pooled 5806-observation, 425-subject analysis dataset drawn from",
      "seven Pfizer B7461-series studies. The parameter estimates in",
      "Table 4 (and the sequential zero-first order absorption plus",
      "auto-induction structure this file implements) are from the",
      "final NONMEM 7.3 FOCEI covariate model after PsN stepwise",
      "covariate selection (SCM alpha 0.05 forward inclusion, 0.01",
      "backward elimination)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural PK -- Chen 2021 Table 4 final-model 'Value' column, in
    # hour / L / L/h units matching the paper. Concentration is reported
    # in ng/mL after the L->mL and mg->ug conversions inside model()
    # (dose mg / vc L = mg/L; multiply by 1000 to get ng/mL, matching
    # the assay units in Chen 2021 Figure 1 and Figure 2).
    # ---------------------------------------------------------------------
    lcl    <- log(9.035)   ; label("Initial (single-dose) clearance CLI at 70 kg (L/h)")                                                            # Chen 2021 Table 4 theta_CLI = 9.035 L/h (bootstrap 95% CI 8.01, 10.06)
    lcl_exp_inf <- log(14.472)  ; label("Maximum induced (steady-state) clearance CLMX at 70 kg (L/h)")                                                  # Chen 2021 Table 4 theta_CLMX = 14.472 L/h (bootstrap 95% CI 12.73, 16.22)
    lcl_exp_kdes  <- log(0.020)   ; label("Auto-induction rate constant cl_exp_kdes (1/h)")                                                                      # Chen 2021 Table 4 theta_IND = 0.020 1/h; abstract '~7.25 days to functional steady state (~5 induction half-lives)'
    lvc    <- log(120.511) ; label("Central volume of distribution V2 at 70 kg (L)")                                                                # Chen 2021 Table 4 theta_V2 = 120.511 L (bootstrap 95% CI 103.4, 137.7)
    lvp    <- log(154.905) ; label("Peripheral volume of distribution V3 (L)")                                                                       # Chen 2021 Table 4 theta_V3 = 154.905 L (bootstrap 95% CI 134.2, 175.6)
    lq     <- log(22.002)  ; label("Inter-compartmental clearance Q (L/h)")                                                                          # Chen 2021 Table 4 theta_Q = 22.002 L/h (bootstrap 95% CI 17.65, 26.36)
    lka    <- log(3.113)   ; label("First-order absorption rate constant k_a (1/h)")                                                                 # Chen 2021 Table 4 theta_ka = 3.113 1/h (bootstrap 95% CI 2.31, 3.91)
    ld1    <- log(1.148)   ; label("Zero-order absorption input duration D1 (h)")                                                                    # Chen 2021 Table 4 theta_D1 = 1.148 h (bootstrap 95% CI 1.03, 1.26)
    lfdepot <- log(0.759)  ; label("Absolute oral bioavailability F1 (fraction, unitless)")                                                          # Chen 2021 Table 4 theta_F = 0.759 (bootstrap 95% CI 0.67, 0.85; anchored by IV vs oral crossover B7461007)

    # ---------------------------------------------------------------------
    # Allometric size scaling on CL (both CLI and CLMX arms) and V2.
    # Chen 2021 Methods 'Model development': "Allometric BWT correction
    # was included a priori in the base model on CL and V2 by using a
    # scaling factor exponent of 0.75 and 1, respectively, to remove a
    # confounding effect observed between BWT and sex." The NONMEM
    # control stream ListS1 encodes these as fixed literal exponents
    # (no THETA slot). Wrapped in fixed() per fixed-parameter
    # conventions (parameter-names.md 'Fixed parameters').
    # ---------------------------------------------------------------------
    e_wt_cl <- fixed(0.75) ; label("Body-weight allometric exponent on CL (unitless)")                                                        # Chen 2021 NONMEM ListS1: TVCLI = THETA(1)*(BWT/70)**0.75; TVCLMX = THETA(9)*(BWT/70)**0.75
    e_wt_vc <- fixed(1.0)  ; label("Body-weight allometric exponent on V2 (unitless)")                                                        # Chen 2021 NONMEM ListS1: TVV2 = THETA(2)*(BWT/70)

    # ---------------------------------------------------------------------
    # Covariate effects on CL (linear centered ALB and DOSE_LOR_MGD; power
    # on CRCL). Chen 2021 NONMEM ListS1 encodes:
    #   CLBALB  = 1 + THETA(12) * (BALB - 4)      (BALB in g/dL)
    #   CLTDOSE = 1 + THETA(13) * (TDOSE - 100)   (mg/day)
    #   CLWNCL  = (WNCL / 100) ** THETA(14)       (mL/min)
    #   CLCOV   = CLBALB * CLTDOSE * CLWNCL       (shared between CLI/CLMX)
    # ALB centering value is expressed in g/L (canonical SI) so the
    # paper's 0.067 per g/dL becomes 0.00670 per g/L (divide by 10).
    # ---------------------------------------------------------------------
    e_alb_cl      <- 0.00670 ; label("Linear effect of baseline serum albumin on CL (per g/L, centered at 40 g/L)")                                  # Chen 2021 Table 4 theta_BALB_on_CL = 0.067 per g/dL centered at BALB = 4 g/dL; converted to per g/L
    e_dose_lor_cl <- 0.00100 ; label("Linear effect of total daily lorlatinib dose on CL (per mg/day, centered at 100 mg/day)")                       # Chen 2021 Table 4 theta_TDOSE_on_CL = 0.001 per mg/day centered at 100 mg/day
    e_crcl_cl     <- 0.235   ; label("Power effect of weight-normalized creatinine clearance on CL (unitless, reference 100 mL/min)")                # Chen 2021 Table 4 theta_WNCL_on_CL = 0.235 (bootstrap 95% CI 0.146, 0.324)

    # ---------------------------------------------------------------------
    # Covariate effect on ka (linear indicator on PPI co-administration).
    # Chen 2021 NONMEM ListS1: KAPPI = 1 if PPI = 0; KAPPI = 1 + THETA(15)
    # if PPI = 1. THETA(15) = -0.675 -> ka x 0.325 with PPI co-medication.
    # ---------------------------------------------------------------------
    e_ppi_ka <- -0.675 ; label("Linear effect of proton-pump inhibitor use on ka (unitless indicator)")                                              # Chen 2021 Table 4 theta_PPI_on_ka = -0.675 (bootstrap 95% CI -0.851, -0.499)

    # ---------------------------------------------------------------------
    # Inter-individual variability (Chen 2021 Table 4 IIV block, on the
    # log scale). CL and F share a 2x2 block; V2 and V3 share a 2x2
    # block; ka is independent. Bare NONMEM $OMEGA variances / covariances
    # from ListS1 line 102-108:
    #   BLOCK(2): var(CL) 0.0296, cov(CL,F) -0.00557, var(F) 0.0223
    #   BLOCK(2): var(V2) 0.0867, cov(V2,V3) -0.0167, var(V3) 0.1005
    #   diagonal: var(ka) 2.279
    # Note that Chen 2021 Table 4 shows the block off-diagonals as
    # 'omega_F omega_CL' and 'omega_V2 omega_V3' rather than direct
    # covariance labels; the ListS1 covariance values above are the
    # authoritative numbers.
    # ---------------------------------------------------------------------
    etalcl + etalfdepot ~ c(0.030,
                            -0.006, 0.022)                                                                                                          # Chen 2021 Table 4 / ListS1 CL-F block: var(CL) 0.030, cov(CL,F) -0.006, var(F) 0.022 (correlation ~ -0.234)
    etalvc + etalvp     ~ c(0.086,
                            -0.017, 0.101)                                                                                                          # Chen 2021 Table 4 / ListS1 V2-V3 block: var(V2) 0.086, cov(V2,V3) -0.017, var(V3) 0.101 (correlation ~ -0.183)
    etalka              ~ 2.329                                                                                                                     # Chen 2021 Table 4 / ListS1 ka: variance 2.329 (corresponds to 152.6% CV per Table 4 reporting)

    # ---------------------------------------------------------------------
    # Residual error (Chen 2021 ListS1 $ERROR: LOG(IPRED) is fitted with
    # ROUT-switched proportional SD; ROUT = 28 uses THETA(10) = 0.115
    # for IV data, all other ROUT codes use THETA(11) = 0.438 for oral
    # data). NONMEM 'additive-on-log-scale' is equivalent to nlmixr2's
    # proportional error on the linear scale for small SDs; encoded here
    # as a single propSd = the oral-route residual because the labelled
    # lorlatinib route of administration is oral (100 mg q.d. tablet).
    # The lower IV-route SD is documented in the vignette Assumptions
    # and deviations section for reviewers examining the B7461007
    # absolute-bioavailability arm.
    # ---------------------------------------------------------------------
    propSd <- 0.438 ; label("Proportional residual SD, oral route (fraction; log-scale additive in Chen 2021)")                                       # Chen 2021 Table 4 theta_Res_Error_for_PO = 0.438 (bootstrap 95% CI 0.409, 0.467); ListS1 $THETA line 11
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual PK parameters. CLI and CLMX share:
    #      allometric size term (BWT/70)^0.75
    #      covariate multiplier CLCOV = CLBALB * CLTDOSE * CLWNCL
    #      individual random deviation exp(etalcl) (same eta for both)
    #    per Chen 2021 NONMEM ListS1 lines 43, 53, 63-68. CL is not
    #    modified by PPI (that effect is on ka only, per Table 4).
    # -------------------------------------------------------------------
    cl_balb   <- 1 + e_alb_cl      * (ALB          - 40)
    cl_tdose  <- 1 + e_dose_lor_cl * (DOSE_LOR_MGD - 100)
    cl_wncl   <- (CRCL / 100) ^ e_crcl_cl
    cl_cov    <- cl_balb * cl_tdose * cl_wncl

    cli  <- exp(lcl    + etalcl) * (WT / 70) ^ e_wt_cl * cl_cov
    cl_exp_inf <- exp(lcl_exp_inf + etalcl) * (WT / 70) ^ e_wt_cl * cl_cov

    # 2. Time-varying total CL via the paper's auto-induction model
    #    CL(t) = CLI + (CLMX - CLI) * (1 - exp(-cl_exp_kdes * t))
    #    rising from CLI at t = 0 (single-dose regime) to CLMX at
    #    t -> infinity (steady-state / auto-induced regime). Chen 2021
    #    Methods 'Lorlatinib clearance estimation' and Results
    #    paragraph 3 (auto-induction to ~7.25 days = ~5 induction
    #    half-lives). See vignette Assumptions and deviations for the
    #    piecewise-SD-flag NONMEM implementation vs this continuous
    #    canonical form (mathematically equivalent at the SD anchor
    #    and at steady state; the continuous form supports arbitrary
    #    dosing histories in downstream simulation).
    #
    #    `time` in rxode2 corresponds to elapsed simulation time from
    #    t = 0 at the first modeled dose event.
    cl_exp_kdes <- exp(lcl_exp_kdes)
    cl    <- cli + (cl_exp_inf - cli) * (1 - exp(-cl_exp_kdes * time))

    # 3. Central and peripheral distribution (V2, V3, Q). V2 has fixed
    #    allometric exponent 1.0 on BWT; V3 and Q are constant.
    vc <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # 4. Absorption. ka carries a linear PPI multiplier (67.5% ka
    #    reduction with PPI co-administration). D1 sets the zero-order
    #    input duration; F1 = 0.759 is the oral bioavailability anchor
    #    from the IV vs oral crossover study B7461007.
    ka <- exp(lka + etalka) * (1 + e_ppi_ka * CONMED_PPI)
    d1 <- exp(ld1)
    fdepot <- exp(lfdepot + etalfdepot)

    # 5. Two-compartment ODE system with sequential zero- and first-
    #    order absorption. Dose enters `depot` over the zero-order
    #    window of duration D1, then is absorbed first-order at ka to
    #    `central`. Elimination is linear via CL(t) / V2 with time-
    #    varying CL(t).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot) <- d1
    f(depot)   <- fdepot

    # 6. Observation. dose in mg, vc in L -> central / vc is mg/L;
    #    multiply by 1000 to report in ng/mL (matching the assay range
    #    reported in Chen 2021 Figures 1-3 and the typical-patient
    #    steady-state Cmax = 606 ng/mL / AUCtau = 5180 ng.h/mL
    #    reference stated in the sensitivity-analysis Methods section).
    Cc <- central / vc * 1000

    # 7. Residual error. See ini() comment: encoded as a single
    #    proportional SD for the oral route (the labelled route).
    Cc ~ prop(propSd)
  })
}
