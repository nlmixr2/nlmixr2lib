`Llanos-Paez_2020_gentamicin` <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in pediatric oncology and nononcology patients (Llanos-Paez 2020); body composition is described by normal fat mass (NFM = FFM + Ffat * (TBW - FFM)) with separate Ffat estimates for CL (0.48) and V1 (0.10) and Ffat fixed to 0 for Q and V2; CL is driven by Holford 2017 GFR-maturation (PMA-based Hill function) and a power ratio of age/sex-matched physiological mean serum creatinine (Ceriotti 2008) over individual SCR; oncology cohort has 15.4% lower V1 and 32.1% lower Q than nononcology."
  reference <- "Llanos-Paez CC, Staatz CE, Lawson R, Hennig S. Differences in the Pharmacokinetics of Gentamicin between Oncology and Nononcology Pediatric Patients. Antimicrob Agents Chemother. 2020;64(2):e01730-19. doi:10.1128/AAC.01730-19. Structural model carried over from Llanos-Paez et al. 2017 Antimicrob Agents Chemother 61:e00205-17 (doi:10.1128/AAC.00205-17); GFR maturation function from Holford NHG. 2017. Systems pharmacology learning from GAVamycin (PAGANZ TM50 = 46.5 weeks PMA, Hill = 3.43, adult plateau = 119 mL/min); serum-creatinine reference from Ceriotti F, et al. 2008. Clin Chem 54:559-566 (doi:10.1373/clinchem.2007.099648)."
  vignette <- "Llanos-Paez_2020_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (TBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; combined with FFM via NFM = FFM + Ffat * (WT - FFM) for size scaling on CL, V1, Q, V2.",
      source_name        = "TBW"
    ),
    FFM = list(
      description        = "Fat-free mass (Janmahasatian 2005 prediction equation in the source paper; equivalent estimators acceptable)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Reference adult value 56.1 kg. Combined with WT via NFM = FFM + Ffat * (WT - FFM).",
      source_name        = "FFM"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the GFR maturation Hill function on CL via PMA in weeks (PMA = PAGE * 4.35). Holford 2017 plateau parameters: TM50 = 46.5 weeks PMA, Hill = 3.43, adult GFR = 119 mL/min.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Individual patient serum creatinine concentration (SCR_i)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Numerator-vs-denominator role is paired with CREAT_REF: the model uses (CREAT_REF / CREAT)^0.58 as a multiplicative renal-function factor on CL. Source paper replaces values below the 30 umol/L laboratory limit with the Ceriotti 2008 age/sex-matched physiological mean (i.e., CREAT_REF) before fitting.",
      source_name        = "SCR"
    ),
    CREAT_REF = list(
      description        = "Age- and sex-matched physiological mean serum creatinine (Ceriotti 2008 reference)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Population-typical SCR for the patient's age and sex per Ceriotti et al. 2008 Clin Chem 54:559-566. The vignette derives this column from a Ceriotti age-band lookup; for typical-value simulation of a virtual patient with normal renal function, set CREAT_REF = CREAT so the (CREAT_REF / CREAT)^0.58 ratio collapses to 1.",
      source_name        = "SCR_mean"
    ),
    DIS_CANCER_PED = list(
      description        = "Pediatric oncology cohort indicator (1 = oncology cohort, 0 = nononcology)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nononcology pediatric admission; appendicitis and renal/urinary infection were most common in Llanos-Paez 2020)",
      notes              = "Time-fixed per subject. Multiplicative cohort shifts on V1 (factor 1 - 0.154 when 1) and Q (factor 1 - 0.321 when 1); CL has no oncology effect. Renamed from source column ONCOLOGY to the canonical DIS_CANCER_PED per covariate-columns.md.",
      source_name        = "ONCOLOGY"
    )
  )

  population <- list(
    n_subjects     = 538,
    n_studies      = 1,
    age_range      = "PNA 0.45 to 18.4 years (oncology) and 0.62 to 16.8 years (nononcology); PMA mean 361.9 weeks (SD 241.8 weeks)",
    age_median     = "PNA 6.19 years pooled; oncology 6.36 years, nononcology 5.53 years",
    weight_range   = "TBW 4.8 to 102.8 kg (oncology) and 3.4 to 121.0 kg (nononcology); pooled mean 24.6 kg (SD 17.5 kg)",
    weight_median  = "TBW 25.2 kg (oncology) and 22.3 kg (nononcology)",
    sex_female_pct = 45.9,
    race_ethnicity = "Not reported",
    disease_state  = "Pediatric oncology (n = 423; predominantly leukemia 45% and blastomas 13%) pooled with pediatric nononcology admissions (n = 115; appendicitis 12.2%, kidney disease / urinary tract infection 10.4%, multi-factorial others) at the Children's Hospital of Queensland, Brisbane, Australia.",
    dose_range     = "30-min IV infusion of 7.5 mg/kg q24h (patients < 10 years old) or 6 mg/kg q24h (patients >= 10 years old) per local clinical guidelines.",
    regions        = "Australia (single-center retrospective TDM dataset 2008-2013).",
    notes          = "Pooled 423-patient oncology cohort (Llanos-Paez 2017; 2,422 gentamicin concentrations) plus 115-patient nononcology cohort (487 concentrations) collected by routine TDM. 372 (15.4%) oncology and 60 (12.3%) nononcology samples were below the assay LLOQ and replaced by LLOQ/2 prior to model fitting. Demographics from Llanos-Paez 2020 Table 1.",
    bsv_caveat     = "Paper Table 2 reports cohort-stratified BSV CV% on V1 (23.8% oncology vs 26.0% nononcology) and Q (29.4% oncology vs 59.8% nononcology). nlmixr2's eta-variance is population-level and cannot be stratified by a covariate without splitting the model file. This implementation uses single eta variances calibrated to the larger oncology cohort (n = 423); see vignette Assumptions and deviations for the impact on stochastic simulation."
  )

  ini({
    # Structural population estimates -- nononcology baseline so the oncology
    # multiplier (1 + e_cancer_ped_*) acts as a reduction. Values in Table 2.
    lcl <- log(4.58); label("Population CL at NFM_std_CL = 62.8 kg, adult-mature GFR, normal SCR (L/h)")  # Llanos-Paez 2020 Table 2
    lvc <- log(21.4); label("Population V1 (nononcology) at NFM_std_V1 = 57.5 kg (L)")                    # Llanos-Paez 2020 Table 2
    lq  <- log(0.84); label("Population Q (nononcology) at NFM_std_Q = 56.1 kg (L/h)")                    # Llanos-Paez 2020 Table 2
    lvp <- log(18.2); label("Population V2 at NFM_std_V2 = 56.1 kg (L)")                                  # Llanos-Paez 2020 Table 2

    # Fat-mass contribution to NFM for each PK parameter (Holford NFM scheme)
    ffat_cl <- 0.48; label("Fat fraction Ffat for CL (unitless)")  # Llanos-Paez 2020 Table 2
    ffat_v1 <- 0.10; label("Fat fraction Ffat for V1 (unitless)")  # Llanos-Paez 2020 Table 2
    # Ffat for Q and V2 fixed to 0 in the published model; expressed inline in model().

    # Pediatric oncology cohort effects on V1 and Q (multiplicative-from-1 form)
    e_cancer_ped_vc <- -0.154; label("Pediatric oncology effect on V1: V1_oncology / V1_nononcology - 1 (fraction)")  # Llanos-Paez 2020 Table 2 ratio 18.1 / 21.4
    e_cancer_ped_q  <- -0.321; label("Pediatric oncology effect on Q: Q_oncology / Q_nononcology - 1 (fraction)")     # Llanos-Paez 2020 Table 2 ratio 0.57 / 0.84

    # Renal function on CL (power ratio of Ceriotti reference SCR over individual SCR)
    e_creat_cl <- 0.58; label("Exponent on (CREAT_REF / CREAT) ratio for CL (unitless)")  # Llanos-Paez 2020 Table 2 theta_serum_creatinine

    # Inter-individual variability (omega^2 = log(CV^2 + 1); diagonal omega used --
    # the paper extended to a full variance-covariance matrix but only diagonal
    # CV% values are tabulated in the publication; off-diagonals not reproducible).
    etalcl ~ 0.02686  # CV 16.5%   # Llanos-Paez 2020 Table 2 BSV
    etalvc ~ 0.05509  # CV 23.8%   # Llanos-Paez 2020 Table 2 BSV (oncology cohort; nononcology BSV 26.0% not encoded)
    etalq  ~ 0.08288  # CV 29.4%   # Llanos-Paez 2020 Table 2 BSV (oncology cohort; nononcology BSV 59.8% not encoded)
    etalvp ~ 0.37652  # CV 67.6%   # Llanos-Paez 2020 Table 2 BSV

    # Combined residual error
    propSd <- 0.293; label("Proportional residual error (fraction)")  # Llanos-Paez 2020 Table 2
    addSd  <- 0.05;  label("Additive residual error (mg/L)")          # Llanos-Paez 2020 Table 2
  })

  model({
    # Body composition: normal fat mass per Holford / Anderson scheme
    nfm_cl <- FFM + ffat_cl * (WT - FFM)
    nfm_v1 <- FFM + ffat_v1 * (WT - FFM)
    nfm_q  <- FFM
    nfm_v2 <- FFM
    nfm_std_cl <- 56.1 + ffat_cl * (70 - 56.1)
    nfm_std_v1 <- 56.1 + ffat_v1 * (70 - 56.1)
    nfm_std_q  <- 56.1
    nfm_std_v2 <- 56.1

    # Postmenstrual age in weeks from PAGE (months); Holford 2017 GFR maturation Hill
    pma_wks <- PAGE * 4.35
    gfr_mat <- (nfm_cl / nfm_std_cl)^0.75 * pma_wks^3.43 / (46.5^3.43 + pma_wks^3.43) * 119

    # Multiplicative pediatric-oncology cohort effects on V1 and Q (no effect on CL or V2)
    cohort_vc <- 1 + e_cancer_ped_vc * DIS_CANCER_PED
    cohort_q  <- 1 + e_cancer_ped_q  * DIS_CANCER_PED

    # Renal function: ratio of age/sex-matched physiological mean SCR to individual SCR
    scr_factor <- (CREAT_REF / CREAT)^e_creat_cl

    # Individual PK parameters
    cl <- exp(lcl + etalcl) * (gfr_mat / 100) * scr_factor
    vc <- exp(lvc + etalvc) * (nfm_v1 / nfm_std_v1)         * cohort_vc
    q  <- exp(lq  + etalq)  * (nfm_q  / nfm_std_q )^0.75    * cohort_q
    vp <- exp(lvp + etalvp) * (nfm_v2 / nfm_std_v2)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV; gentamicin is administered as a 30-min infusion into central
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
