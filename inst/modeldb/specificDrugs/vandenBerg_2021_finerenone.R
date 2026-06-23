vandenBerg_2021_finerenone <- function() {
  description <- paste(
    "Two-compartment population PK model for oral finerenone (Bayer BAY",
    "94-8862, a non-steroidal selective mineralocorticoid receptor",
    "antagonist) in adults with chronic kidney disease and type 2 diabetes,",
    "developed on n=2284 subjects / 5057 sparse PK observations from the",
    "FIDELIO-DKD Phase III trial (NCT02540993). Absorption is modelled via",
    "a chain of four sequential first-order steps (depot + three transit",
    "buffers, all at common rate Ka = 22.5 1/h, mean transit time MTT =",
    "n_steps / Ka = 0.178 h with the depot counted as the first compartment",
    "in the chain) preceded by a fixed 0.215 h absorption lag time; the",
    "central-peripheral disposition is two-compartment with the peripheral",
    "volume fixed equal to the central volume Vp/F = Vc/F (ratio fixed at",
    "1). Covariates retained in the final model are body weight and Korean",
    "ethnicity on Vc/F; time-varying eGFR-CKD-EPI, body height, serum",
    "creatinine, smoking status (current or former vs never), long-term",
    "(>=50% of treatment period) SGLT2 inhibitor use, gamma",
    "glutamyl-transferase, and a two-tier CYP3A4-inhibitor coadministration",
    "categorisation (strong/moderate/weak inhibitor >=50% of treatment",
    "period vs any other inhibitor exposure) on CL/F; with each of the CL/F",
    "covariates (except GGT) ALSO applied inversely to the relative",
    "bioavailability F1 in the paper's NONMEM control stream (so the",
    "covariate appears on both CL/F and F simultaneously, the net effect on",
    "steady-state AUC scales as 1 / covariate-factor^2 and the net effect",
    "on Cmax scales as 1 / covariate-factor). Inter-individual variability",
    "is a 2x2 block on CL/F and Vc/F (omega^2 0.0961 / 0.104, covariance",
    "0.0442, correlation ~0.44); no IIV on Ka or absorption. Residual error",
    "is proportional (sigma^2 = 0.313, propSd = sqrt(0.313) ~= 0.5595)."
  )
  reference <- paste(
    "van den Berg P, Ruppert M, Mesic E, Snelder N, Seelmann A, Heinig R,",
    "Joseph A, Garmann D, Lippert J, Eissing T (2022).",
    "Finerenone Dose-Exposure-Response for the Primary Kidney Outcome in",
    "FIDELIO-DKD Phase III: Population Pharmacokinetic and Time-to-Event",
    "Analysis.",
    "Clin Pharmacokinet 61(7):943-955.",
    "doi:10.1007/s40262-021-01082-2",
    sep = " "
  )
  vignette <- "vandenBerg_2021_finerenone"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight used for power-form effect on Vc/F:",
        "(WT / 85)^0.501 with reference 85 kg per FIDELIO-DKD popPK NONMEM",
        "control stream in ESM ('CV1 = (BW0/85)**THETA(9)'). Effect on Vc/F",
        "is independent of the CL/F-and-F covariate cluster."
      ),
      source_name        = "BW0"
    ),
    HT = list(
      description        = "Body height (baseline)",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body height used for power-form effect on CL/F AND F:",
        "(HT / 167)^0.720 with reference 167 cm per FIDELIO-DKD popPK",
        "NONMEM control stream in ESM ('CV3 = (HGHT/167)**THETA(11)'). The",
        "same factor is applied to CL/F (multiplicative) and to F1",
        "(inversely), so the net effect on AUC at steady state scales as",
        "(HT / 167)^(-2 * 0.720)."
      ),
      source_name        = "HGHT"
    ),
    CRCL = list(
      description        = "Time-varying eGFR-CKD-EPI (BSA-normalized estimated glomerular filtration rate)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying estimated glomerular filtration rate computed via",
        "the CKD-EPI equation (BSA-normalized to mL/min/1.73 m^2) and used",
        "for power-form effect on CL/F AND F: (CRCL / 39.1)^0.155 with",
        "reference 39.1 mL/min/1.73 m^2 (FIDELIO-DKD cohort median) per",
        "ESM NONMEM control stream ('CV2 = (EGFREP/39.1)**THETA(10)').",
        "Renal function declines over the trial follow-up (median 2.6",
        "years); the time-varying form was retained over baseline-only",
        "eGFR. Same factor applied to CL/F (multiplicative) and F1",
        "(inversely)."
      ),
      source_name        = "EGFREP"
    ),
    CREAT = list(
      description        = "Serum creatinine (baseline)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline serum creatinine in mg/dL used for power-form effect on",
        "CL/F AND F: (CREAT / 1.51)^0.118 with reference 1.51 mg/dL per",
        "ESM NONMEM control stream ('CV4 = (CREA/1.51)**THETA(12)'). The",
        "paper interprets this as a residual marker of impaired CYP3A4",
        "metabolism beyond what eGFR alone captures (Discussion paragraph",
        "on 'kidney markers as surrogates representing a more complex",
        "underlying pathophysiology where likely accumulating uremic",
        "toxins ultimately impair CYP3A4 metabolism')."
      ),
      source_name        = "CREA"
    ),
    GGT = list(
      description        = "Serum gamma-glutamyl-transferase (baseline)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline serum gamma-glutamyl-transferase used for power-form",
        "effect on CL/F ONLY: (GGT / 25)^(-0.0694) with reference 25 U/L",
        "per ESM NONMEM control stream ('CV5 = (GGT/25)**THETA(16)').",
        "Unlike the other CL/F covariates, GGT is NOT applied to F1",
        "(consistent with the NONMEM code's CL-only effect): see model",
        "body 'TVCL = ... * CV5 * ...' but 'F1 = ... / (CV2*CV3*CV4*",
        "ESMOK*ESGLT*ECYPINHR)' with CV5 absent from F1."
      ),
      source_name        = "GGT"
    ),
    SMOKE_NEVER = list(
      description        = "Never-smoker indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (never smoker is the reference for the paper's ever-vs-never effect)",
      notes              = paste(
        "FIDELIO-DKD paper's 'effect of smoking (current or former",
        "smokers)' is encoded as a single (1 - SMOKE_NEVER) ever-smoker",
        "indicator: the paper's 1.04 effect applies to subjects with",
        "SMOKE_NEVER = 0 (i.e. current OR former smokers, NONMEM SMOK >=",
        "2) and the reference (no effect) is SMOKE_NEVER = 1 (never",
        "smokers, NONMEM SMOK = 1). The paper does not separately",
        "estimate a current-vs-former effect, so the paired canonical",
        "SMOKE_CURRENT (which would split former from current within the",
        "ever-smoker group) is not needed in this model. Effect on both",
        "CL/F (multiplicative) and F1 (inversely)."
      ),
      source_name        = "SMOK"
    ),
    CONMED_SGLT2I = list(
      description        = "Concomitant SGLT2 inhibitor coadministration indicator (long-term, >=50% of on-treatment period)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no long-term SGLT2 inhibitor coadministration)",
      notes              = paste(
        "1 = subject coadministered a sodium-glucose cotransporter 2",
        "inhibitor (canagliflozin, dapagliflozin, empagliflozin, or",
        "another systemic SGLT2i) for at least 50% of the on-finerenone",
        "treatment period; 0 otherwise. The >=50%-of-treatment-period",
        "threshold is the FIDELIO-DKD analysis-specific definition (see",
        "Table 2 'Effect of SGLT2 inhibitor use (>=50% of on treatment",
        "period) on CL/F and F') and the dataset column GLYCSI_CATN = 2",
        "in the ESM NONMEM control stream encodes it. Effect on CL/F",
        "(multiplicative, 1.10) and F1 (inversely)."
      ),
      source_name        = "SGLT"
    ),
    CONMED_CYP3A4_INH_HI = list(
      description        = "Concomitant CYP3A4 inhibitor (strong/moderate/weak) coadministration at >=50% of on-treatment period",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor coadministration in this strength-and-exposure category; CONMED_CYP3A4_INH_LO may still be 1)",
      notes              = paste(
        "1 = subject coadministered a CYP3A4 inhibitor classified as",
        "strong, moderate, or weak (according to FDA/EMA classifications)",
        "for at least 50% of the on-finerenone treatment period; 0",
        "otherwise. Pooled across strong/moderate/weak strengths because",
        "the FIDELIO-DKD analysis did not separate them. ESM NONMEM",
        "control stream encodes via CYPINH IN (4, 6, 8) corresponding to",
        "strong/moderate/weak inhibitor >=50% respectively. Effect on",
        "CL/F (multiplicative, 0.951) and F1 (inversely). Mutually",
        "exclusive with CONMED_CYP3A4_INH_LO (a subject is in at most",
        "one of the two categories)."
      ),
      source_name        = "CYPINH"
    ),
    CONMED_CYP3A4_INH_LO = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration in 'other' category (unclassified at any duration, or strong/moderate/weak below the 50% threshold)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor coadministration in this category; CONMED_CYP3A4_INH_HI may still be 1)",
      notes              = paste(
        "1 = subject coadministered a CYP3A4 inhibitor in any of the",
        "following sub-categories: unclassified inhibitor at any duration",
        "(CYPINH = 1 or 2 in ESM NONMEM code), or strong/moderate/weak",
        "inhibitor present for LESS THAN 50% of the on-finerenone",
        "treatment period (CYPINH = 3, 5, 7 in ESM NONMEM code); 0",
        "otherwise. Effect on CL/F (multiplicative, 0.996) and F1",
        "(inversely). Mutually exclusive with CONMED_CYP3A4_INH_HI."
      ),
      source_name        = "CYPINH"
    ),
    RACE_KOREAN = list(
      description        = "Korean-heritage race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Korean race, including other Asian, White, Black, etc.)",
      notes              = paste(
        "1 = Korean-heritage subject (RACA = 3.3 in ESM NONMEM control",
        "stream), 0 otherwise. The FIDELIO-DKD popPK analysis found a",
        "single-race effect on Vc/F for Korean subjects only (no other",
        "race or ethnicity reached significance). The paper Discussion",
        "notes 'this may be a spurious finding based on limited data",
        "(only 2.4% of subjects were Korean)'. Effect on Vc/F only (NOT",
        "applied to CL/F or F)."
      ),
      source_name        = "RACA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 2284L,
    n_studies       = 1L,
    n_observations  = 5057L,
    disease_state   = "chronic kidney disease and type 2 diabetes mellitus (FIDELIO-DKD eligibility: eGFR 25 to <75 mL/min/1.73 m^2 with persistent moderately or severely elevated albuminuria, on maximally tolerated renin-angiotensin system inhibitor)",
    dose_range      = "10 or 20 mg finerenone QD oral (starting dose 10 mg if eGFR <60, 20 mg if eGFR >=60; up- or down-titrated by potassium / eGFR), average 15.1 mg/day across follow-up",
    follow_up       = "median 2.6 years",
    egfr_range      = "median (5th-95th) 43.0 (26.7-66.9) mL/min/1.73 m^2 at baseline (FIDELIO-DKD population)",
    uacr_range      = "median (5th-95th) 852 (140-3366) mg/g at baseline (FIDELIO-DKD population)",
    reference_weight     = "85 kg (cohort median, used as ref for WT power-form on Vc/F)",
    reference_height     = "167 cm (cohort median, used as ref for HT power-form on CL/F and F)",
    reference_creatinine = "1.51 mg/dL (cohort median, used as ref for CREAT power-form on CL/F and F)",
    reference_egfr       = "39.1 mL/min/1.73 m^2 (cohort median time-varying value, used as ref for CRCL power-form on CL/F and F)",
    reference_ggt        = "25 U/L (cohort median, used as ref for GGT power-form on CL/F)",
    sampling_design = "sparse: trough at month 4, post-dose at any time on yearly visit days",
    regions         = "international (FIDELIO-DKD enrolled in North America, EU, Latin America, and Asia-Pacific)",
    notes           = paste(
      "Population is the FIDELIO-DKD per-protocol Phase III analysis",
      "set (5734 randomized total, 5674 valid for analysis, 2833 on",
      "finerenone of whom 2284 had at least one valid PK sample after",
      "outlier exclusions; see Results 'Clinical Study' and 'Population",
      "PK Modeling and Simulation' paragraphs). Demographic ranges",
      "(age, sex, race percentages) not separately reported in the popPK",
      "section; see Bakris et al. 2020 NEJM 383:2219-2229 (the main",
      "FIDELIO-DKD efficacy publication) for full Table 1 baseline",
      "demographics."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters - paper Table 2 final estimates
    # =========================================================================
    lka     <- log(22.5);  label("Common first-order transit-and-absorption rate constant Ka (1/h)")     # van den Berg 2021 Table 2 'Ka (1/h) = 22.5 (RSE 16.2%)'; also ESM NONMEM '$THETA TH2 22.6' (final-vs-draft difference 22.5 vs 22.6)
    lcl     <- log(29.9);  label("Apparent oral clearance CL/F at reference covariates (L/h)")            # van den Berg 2021 Table 2 'CL/F (L/h) = 29.9 (RSE 3.62%)'
    lvc     <- log(113);   label("Apparent central volume of distribution Vc/F at reference covariates (L)") # van den Berg 2021 Table 2 'Vc/F (L) = 113 (RSE 2.79%)'
    lq      <- log(0.335); label("Apparent inter-compartmental clearance Q/F (L/h)")                      # van den Berg 2021 Table 2 'Q/F (L/h) = 0.335 (RSE 9.28%)'

    # =========================================================================
    # Fixed structural parameters (FIDELIO-DKD popPK analysis fixed these from
    # prior Phase 1/2 information or for identifiability; see paper Results
    # 'Population PK Modeling and Simulation' paragraph and ESM 'Volume of
    # Distribution' subsection)
    # =========================================================================
    lvp_ratio <- fixed(log(1));     label("Log-ratio Vp/F to Vc/F (FIXED at 1, peripheral volume equals central volume)") # van den Berg 2021 Table 2 'Ratio Vp/F and Vc/F = 1 (fixed)' and ESM 'Volume of Distribution': 'two volumes of distribution were assumed equal as the data did not allow the estimation of two separate volumes ... supported by prior Phase 1 and 2 data analyses'
    lalag1    <- fixed(log(0.215)); label("Absorption lag time ALAG1 (h, FIXED at 0.215)")               # van den Berg 2021 Table 2 'Absorption lag time (h) = 0.215 (fixed)'; ESM NONMEM '$THETA TH7 0.215 FIX'
    lfdepot   <- fixed(log(1));     label("Reference relative bioavailability F1 (FIXED at 1; oral-only data)") # van den Berg 2021 Table 2 'Relative bioavailability = 1 (fixed)'; ESM NONMEM '$THETA TH8 1 FIX'

    # =========================================================================
    # Covariate effects - power-form exponents (continuous) and multiplicative
    # log-factors (categorical); all values are paper Table 2 final estimates.
    # NOTE on paper-specific parameterisation: every CL/F covariate EXCEPT GGT
    # is applied in IDENTICAL multiplicative form to BOTH CL/F (positive sense)
    # and to F1 (inverse sense). See ESM NONMEM control stream:
    #   TVCL = THETA(3) * CV1 * CV2 * CV3 * CV4 * ESMOK * CV5 * ESGLT * ECYPINHR
    #   F1   = THETA(8) / (CV2 * CV3 * CV4 * ESMOK * ESGLT * ECYPINHR)
    # The net effect on steady-state AUC = F1 * Dose / CL is 1 / cov_factor^2
    # for each shared covariate (matches the paper's reported 17.1% AUC drop
    # for SGLT2i+: 1/1.10^2 = 0.826 = 17.4% drop, see Results paragraph 4).
    # =========================================================================
    e_wt_vc      <- 0.501;   label("Power exponent of WT on Vc/F (unitless)")                            # van den Berg 2021 Table 2 'Effect of body weight on Vc/F = 0.501 (RSE 9.87%)'; centered at WT 85 kg
    e_egfr_clf   <- 0.155;   label("Power exponent of time-varying eGFR-CKD-EPI on CL/F and F (unitless)") # van den Berg 2021 Table 2 'Effect of eGFR-EPI (time-varying) on CL/F and F = 0.155 (RSE 20.1%)'; centered at CRCL 39.1 mL/min/1.73 m^2
    e_ht_clf     <- 0.720;   label("Power exponent of HT on CL/F and F (unitless)")                       # van den Berg 2021 Table 2 'Effect of body height on CL/F and F = 0.720 (RSE 16.1%)'; centered at HT 167 cm
    e_creat_clf  <- 0.118;   label("Power exponent of CREAT on CL/F and F (unitless)")                    # van den Berg 2021 Table 2 'Effect of creatinine on CL/F and F = 0.118 (RSE 38.4%)'; centered at CREAT 1.51 mg/dL
    e_korean_vc  <- log(1.29);  label("Log-multiplicative effect of Korean ethnicity on Vc/F (unitless)") # van den Berg 2021 Table 2 'Effect for Korean subjects on Vc/F = 1.29 (RSE 6.79%)'; encoded as log() because the NONMEM form is TH13^RACE_KOREAN i.e. multiplied by 1.29 only when RACE_KOREAN = 1
    e_sglt_clf   <- log(1.10);  label("Log-multiplicative effect of long-term SGLT2 inhibitor use (CONMED_SGLT2I = 1) on CL/F and F (unitless)") # van den Berg 2021 Table 2 'Effect of SGLT2 inhibitor use (>=50% of on treatment period) on CL/F and F = 1.10 (RSE 2.23%)'
    e_smoke_clf  <- log(1.04);  label("Log-multiplicative effect of ever-smoker status (current or former) on CL/F and F (unitless)") # van den Berg 2021 Table 2 'Effect of smoking (current or former smokers) on CL/F and F = 1.04 (RSE 1.11%)'
    e_ggt_clf    <- -0.0694; label("Power exponent of GGT on CL/F (unitless; CL/F only, NOT applied to F)") # van den Berg 2021 Table 2 'Effect of GGT on CL/F = -0.0694 (RSE 16.2%)'; centered at GGT 25 U/L
    e_cypinhi_clf <- log(0.951); label("Log-multiplicative effect of strong/moderate/weak CYP3A4-inhibitor use >=50% of treatment (CONMED_CYP3A4_INH_HI = 1) on CL/F and F (unitless)") # van den Berg 2021 Table 2 'Effect of CYP3A4 inhibitor use (weak, moderate, or strong >=50% of on-treatment period) on CL/F and F = 0.951 (RSE 1.66%)'
    e_cypinlo_clf <- log(0.996); label("Log-multiplicative effect of 'other' CYP3A4-inhibitor categories (CONMED_CYP3A4_INH_LO = 1) on CL/F and F (unitless)") # van den Berg 2021 Table 2 'Effect of CYP3A4 inhibitor use (other categories) on CL/F and F = 0.996 (RSE 2.17%)'

    # =========================================================================
    # Inter-individual variability - 2x2 block on (log)CL/F and (log)Vc/F.
    # The ARTS-DN model's IIV on Ka was DROPPED from the FIDELIO-DKD popPK
    # model (Results paragraph 2: 'the inter-individual variability parameter
    # on the absorption rate was dropped as data were not sufficiently
    # informative').
    # =========================================================================
    etalcl + etalvc ~ c(0.0961,
                        0.0442, 0.104)  # van den Berg 2021 Table 2 'omega^2 CL/F = 0.0961 (CV 31.8%, shrinkage 26.6%, RSE 7.30%)' / 'Covariance CL/F x Vc/F = 0.0442 (RSE 19.4%)' / 'omega^2 Vc/F = 0.104 (CV 33.0%, shrinkage 47.8%, RSE 18.4%)'; correlation ~= 0.442

    # =========================================================================
    # Residual error - proportional only. NONMEM SIGMA is FIXED at 1 with the
    # variance carried by a THETA: $ERROR sets W = SQRT(THETA(1)) * IPRED and
    # Y = IPRED + W * ERR(1), so the effective proportional SD is sqrt(sigma^2).
    # =========================================================================
    propSd <- sqrt(0.313);  label("Proportional residual error (fraction)")  # van den Berg 2021 Table 2 'sigma^2 = 0.313 (RSE 2.83%)' with NONMEM SIGMA fixed at 1 (ESM '$SIGMA 1 FIX') and W = sqrt(THETA(1)) * IPRED; propSd = sqrt(0.313) = 0.5595
  })

  model({
    # ------------------------------------------------------------------------
    # 1. Derived covariate terms
    #    Continuous covariates: power-form centered at cohort median.
    #    Categorical covariates: indicator gates the multiplicative factor.
    #    SMOKE_EVER derived from canonical SMOKE_NEVER (paper's "current or
    #    former" = (1 - SMOKE_NEVER), reference = never smoker).
    # ------------------------------------------------------------------------
    cv_wt    <- (WT    / 85  )^e_wt_vc
    cv_egfr  <- (CRCL  / 39.1)^e_egfr_clf
    cv_ht    <- (HT    / 167 )^e_ht_clf
    cv_creat <- (CREAT / 1.51)^e_creat_clf
    cv_ggt   <- (GGT   / 25  )^e_ggt_clf

    smoke_ever <- 1 - SMOKE_NEVER
    e_korean_eff  <- exp(e_korean_vc  * RACE_KOREAN)
    e_sglt_eff    <- exp(e_sglt_clf   * CONMED_SGLT2I)
    e_smoke_eff   <- exp(e_smoke_clf  * smoke_ever)
    e_cypinhi_eff <- exp(e_cypinhi_clf * CONMED_CYP3A4_INH_HI)
    e_cypinlo_eff <- exp(e_cypinlo_clf * CONMED_CYP3A4_INH_LO)

    cov_cl_and_f <- cv_egfr * cv_ht * cv_creat * e_smoke_eff * e_sglt_eff * e_cypinhi_eff * e_cypinlo_eff

    # ------------------------------------------------------------------------
    # 2. Individual PK parameters
    # ------------------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * cov_cl_and_f * cv_ggt
    vc <- exp(lvc + etalvc) * cv_wt * e_korean_eff
    vp <- vc * exp(lvp_ratio)
    q  <- exp(lq)

    fdepot <- exp(lfdepot) / cov_cl_and_f
    alag1  <- exp(lalag1)

    # ------------------------------------------------------------------------
    # 3. Micro-constants
    # ------------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------------
    # 4. ODE system - depot + three transit buffers (all rate ka) into central
    #    plus two-compartment disposition (central + peripheral1). NONMEM ESM:
    #      K14 = K45 = K56 = K62 = Ka
    #    with COMP = (DEPOT, CENTRAL, PERI, BUFFER, BUFFER2, BUFFER3)
    #    DEPOT (1) -> BUFFER (4) -> BUFFER2 (5) -> BUFFER3 (6) -> CENTRAL (2)
    #    plus CENTRAL <-> PERI (3) via Q.
    # ------------------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(transit1)     <-  ka * depot    - ka * transit1
    d/dt(transit2)     <-  ka * transit1 - ka * transit2
    d/dt(transit3)     <-  ka * transit2 - ka * transit3
    d/dt(central)      <-  ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------------
    # 5. Bioavailability and absorption lag time on the dose-receiving depot.
    # ------------------------------------------------------------------------
    f(depot)    <- fdepot
    alag(depot) <- alag1

    # ------------------------------------------------------------------------
    # 6. Observation and residual error
    #    Dose in mg, Vc/F in L gives central / vc in mg/L. Paper reports
    #    plasma finerenone in ug/L, so multiply by 1000.
    # ------------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
