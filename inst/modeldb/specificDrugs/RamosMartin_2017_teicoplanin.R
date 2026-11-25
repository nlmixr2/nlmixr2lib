RamosMartin_2017_teicoplanin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for teicoplanin in neonates, infants and",
    "children (0.01-16 years) with suspected or proven staphylococcal sepsis, fitted",
    "non-parametrically with NPAG in Pmetrics to serve as the Bayesian prior of a BestDose",
    "dose-individualisation controller. The elimination rate constant is allometrically",
    "scaled, kel = Ke0 * (WT/70)^-0.25 * Renal, where Renal is postnatal age (days) /",
    "serum creatinine (umol/L) below 3 months of age and eGFR^pw (bedside Schwartz-Haycock)",
    "otherwise; central volume scales linearly with weight. Ke0 and V0 take separate values",
    "for samples drawn before and from 96 h after the first dose (the paper's two",
    "'occasions'). Typical values are the NPAG medians; inter-individual variability is a",
    "log-normal approximation to the reported non-parametric SDs. Residual error is",
    "fixed(0) because the Pmetrics assay-error polynomial was not published.",
    sep = " "
  )
  reference <- paste(
    "Ramos-Martin V, Neely MN, Padmore K, Peak M, Beresford MW, Turner MA, Paulus S,",
    "Hope WW. Tools for the Individualized Therapy of Teicoplanin for Neonates and Children.",
    "Antimicrob Agents Chemother. 2017;61(10):e00707-17. doi:10.1128/AAC.00707-17.",
    "PMID: 28760897; PMCID: PMC5610524.",
    sep = " "
  )
  vignette <- "RamosMartin_2017_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size: kel scales as (WT/70)^-0.25 and the central volume as (WT/70)^1,",
        "both exponents fixed (Methods, equations 2a/2b and the Table 3 footnote). Cohort",
        "median 7.5 kg (range 0.69-62.2; Table 1).",
        sep = " "
      ),
      source_name = "wt"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the renal-function descriptor (PNA < 3 months, i.e. age < 0.25 year, uses",
        "PNA/SCr; otherwise eGFR^pw) and enters the neonatal/young-infant descriptor as",
        "PNA in DAYS divided by serum creatinine in umol/L. model() converts the canonical",
        "months to days with PNA * 30.4375. The Methods text says 'PNA (in years)'; the",
        "maintainers read days (Table 1 reports PNA in days), because PNA in years would put",
        "the young-infant Ke0 support points roughly 365-fold above the older children's,",
        "which the reported Ke0 distribution (mean 0.038, SD 0.04 /h; Table 3) excludes.",
        "See the vignette. Cohort PNA 4-5774 days; neonates median 17 days (Table 1).",
        sep = " "
      ),
      source_name = "PNA"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only below 3 months postnatal age, as the denominator of the PNA/SCr renal",
        "descriptor (no centring value). Cohort median 41 umol/L (range 21-308); neonates",
        "median 44.5 (21-265) (Table 1).",
        sep = " "
      ),
      source_name = "SCr"
    ),
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate by the Haycock-Schwartz bedside formula,",
        "k * height (cm) / serum creatinine, BSA-normalised"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only at postnatal age >= 3 months, raised UNCENTRED to the estimated power pw",
        "(Renal = eGFR^pw; Methods and Table 3 footnote). k = 0.33 (preterm neonates),",
        "0.45 (term neonates) and 0.41 (older children); missing heights were replaced by",
        "the UK growth-chart median for age and sex (Methods; Table 1 footnote b). Cohort",
        "median 62.06 (5.4-160.3); children > 1 month median 78.94 (6.43-160.3) (Table 1).",
        sep = " "
      ),
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened against the base-model Bayesian posterior clearance and volume; no",
        "relationship, so not included in the final model (Results, 'Population PK of",
        "teicoplanin in neonates and children').",
        sep = " "
      )
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      notes = paste(
        "A sigmoidal PMA maturation factor on Ke was explored and not incorporated in the",
        "final model (Methods; Fig. 1J). Neonates median 37 weeks (26-44; Table 1).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 57L,
    n_studies = 2L,
    n_observations = 394L,
    age_range = "0.01-15.82 years (median 0.88); 18 neonates, 16 infants/toddlers 1-23 months, 20 children 2-11 years, 3 children 11-16 years",
    weight_range = "0.69-62.2 kg (median 7.5); neonates median 2.04 kg (0.69-5.08)",
    ga_range = "Neonates: postmenstrual age 26-44 weeks (median 37)",
    sex_female_pct = 42.1,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Proven or suspected methicillin-resistant staphylococcal sepsis (coagulase-negative",
      "staphylococci or MRSA) and/or central-line-associated bloodstream infection; 40% from",
      "the intensive-care unit (mostly after cardiac surgery), 30% oncology patients with",
      "febrile neutropenia, the rest general medical or cardiac.",
      sep = " "
    ),
    renal_function = "eGFR (Haycock-Schwartz) median 62.06 mL/min/1.73 m^2 (5.4-160.3); serum creatinine median 41 umol/L (21-308)",
    dose_range = paste(
      "Neonates (PMA <= 44 weeks): 16 mg/kg loading then 8 mg/kg once daily as a 30-min",
      "infusion. Children > 1 month: 10 mg/kg every 12 h for three doses then 10 mg/kg once",
      "daily as a 2-5 min intravenous bolus.",
      sep = " "
    ),
    regions = "United Kingdom (Alder Hey Children's NHS Foundation Trust and Liverpool Women's Hospital, Liverpool)",
    notes = paste(
      "Pooled from two prospective studies (EudraCT 2012-005738-12), April 2013-January 2015,",
      "previously analysed as separate neonatal and paediatric models. Samples during the first",
      "dosing interval and at steady state (1, 3, 6, 24 h post-dose) on days 3-7; FPIA assay,",
      "LOQ 3 mg/L, range 3-100 mg/L.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural values are the MEDIANS of the NPAG non-parametric marginal
    # distributions of the final model, Ramos-Martin 2017 Table 3 (the Abstract
    # quotes the medians as 'the median population PK parameters'). The suffix
    # _lt96h / _ge96h is the paper's occasion 1 / occasion 2: samples drawn
    # < 96 h and >= 96 h after treatment initiation (Methods; Table 3 footnote).
    # ------------------------------------------------------------------------
    lkel_lt96h <- log(0.03); label("Weight-standardised elimination rate constant Ke0, < 96 h after the first dose (1/h)") # Table 3, Ke01 median 0.03 (mean 0.038, SD 0.04)
    lkel_ge96h <- log(0.025); label("Weight-standardised elimination rate constant Ke0, >= 96 h after the first dose (1/h)") # Table 3, Ke02 median 0.025 (mean 0.036, SD 0.03)
    lvc_lt96h <- log(19.52); label("Weight-standardised central volume V0 at 70 kg, < 96 h after the first dose (L)") # Table 3, V01 median 19.52 (mean 22.636, SD 14.23)
    lvc_ge96h <- log(22.7); label("Weight-standardised central volume V0 at 70 kg, >= 96 h after the first dose (L)") # Table 3, V02 median 22.7 (mean 22.472, SD 7.95)
    lk12 <- log(0.23); label("Central-to-peripheral rate constant Kcp (1/h)") # Table 3, Kcp median 0.23 (mean 0.490, SD 0.48)
    lk21 <- log(0.12); label("Peripheral-to-central rate constant Kpc (1/h)") # Table 3, Kpc median 0.12 (mean 0.214, SD 0.27)
    e_crcl_kel <- 0.07; label("Power exponent pw on uncentred eGFR in the renal descriptor of kel, age >= 3 months (unitless)") # Table 3, pw median 0.07 (mean 0.125, SD 0.15)

    # Allometric exponents fixed, not estimated (Methods: 'we chose to fix these
    # values'; equations 2a/2b; Table 3 footnote).
    e_wt_kel <- fixed(-0.25); label("Allometric exponent of (WT/70) on kel (unitless)") # Eq. 2a and Table 3 footnote, fixed
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on the central volume (unitless)") # Table 3 footnote 'Volume = V0 * (wt/70)', fixed

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete joint density, not a
    # parametric omega; Table 3 prints the mean and SD of each marginal. Carried as a
    # LOG-NORMAL approximation with CV = SD/mean and omega^2 = log(CV^2 + 1). No
    # covariance between parameters or between the two occasions is reported, so all
    # etas are independent.
    # ------------------------------------------------------------------------
    etalkel_lt96h ~ 0.745755 # Table 3 Ke01: CV = 0.04/0.038 = 1.0526 -> log(1.0526^2 + 1)
    etalkel_ge96h ~ 0.527355 # Table 3 Ke02: CV = 0.03/0.036 = 0.8333 -> log(0.8333^2 + 1)
    etalvc_lt96h ~ 0.333034 # Table 3 V01: CV = 14.23/22.636 = 0.6286 -> log(0.6286^2 + 1)
    etalvc_ge96h ~ 0.117921 # Table 3 V02: CV = 7.95/22.472 = 0.3538 -> log(0.3538^2 + 1)
    etalk12 ~ 0.672740 # Table 3 Kcp: CV = 0.48/0.490 = 0.9796 -> log(0.9796^2 + 1)
    etalk21 ~ 0.952369 # Table 3 Kpc: CV = 0.27/0.214 = 1.2617 -> log(1.2617^2 + 1)
    etae_crcl_kel ~ 0.891998 # Table 3 pw: CV = 0.15/0.125 = 1.2 -> log(1.2^2 + 1); multiplicative on pw

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported: observations were weighted by
    # the inverse of the estimated assay variance (Methods), but neither the Pmetrics
    # assay-error polynomial coefficients nor a gamma/lambda term are printed.
    # Carried as fixed(0) rather than invented.
    # ------------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    # 1. Occasion: the paper's IOV = 1 for samples < 96 h and IOV = 2 for samples
    #    >= 96 h after treatment initiation (Methods). Parameters switch as a step.
    occ_ge96h <- 0
    if (tafd() >= 96) occ_ge96h <- 1

    # 2. Renal-function descriptor (Methods, text after Eq. 2b): PNA (days) / SCr
    #    (umol/L) below 3 months (age < 0.25 year), otherwise eGFR^pw.
    pna_days <- PNA * 30.4375
    e_crcl_kel_i <- e_crcl_kel * exp(etae_crcl_kel)
    renal <- CRCL^e_crcl_kel_i
    if (PNA < 3) renal <- pna_days / CREAT

    # 3. Individual parameters; Ke0 and V0 are occasion-specific (Table 3).
    kel0_lt96h <- exp(lkel_lt96h + etalkel_lt96h)
    kel0_ge96h <- exp(lkel_ge96h + etalkel_ge96h)
    v0_lt96h <- exp(lvc_lt96h + etalvc_lt96h)
    v0_ge96h <- exp(lvc_ge96h + etalvc_ge96h)
    kel0 <- kel0_lt96h * (1 - occ_ge96h) + kel0_ge96h * occ_ge96h
    v0 <- v0_lt96h * (1 - occ_ge96h) + v0_ge96h * occ_ge96h
    kel <- kel0 * (WT / 70)^e_wt_kel * renal
    vc <- v0 * (WT / 70)^e_wt_vc
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 4. Two-compartment ODEs (Eqs. 2a and 2b); infusions enter central.
    d/dt(central) <- -(k12 + kel) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation: Y(1) = X(1)/V (text after Eq. 2b).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
