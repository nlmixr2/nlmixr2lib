Luo_2021_daratumumab <- function() {
  description <- "Two-compartment population PK model for subcutaneous (with rHuPH20) and intravenous daratumumab (anti-CD38 IgG1k) in adults with multiple myeloma, with first-order SC absorption and bioavailability, and parallel linear and Michaelis-Menten eliminations from the central compartment whose maximum velocity decays mono-exponentially over time (target depletion). Fit jointly to PAVO Part 2, MMY1008, COLUMBA (SC and IV arms) and PLEIADES (Luo 2021)."
  reference <- "Luo MM, Usmani SZ, Mateos MV, Nahi H, Chari A, San-Miguel J, Touzeau C, Suzuki K, Kaiser M, Carson R, Heuck C, Qi M, Zhou H, Sun YN, Parasrampuria DA. Exposure-Response and Population Pharmacokinetic Analyses of a Novel Subcutaneous Formulation of Daratumumab Administered to Multiple Myeloma Patients. J Clin Pharmacol. 2021;61(5):614-627. doi:10.1002/jcph.1771"
  vignette <- "Luo_2021_daratumumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "daratumumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "daratumumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "daratumumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on linear CL (exponent 1.24) and on V1 (exponent 0.91), normalised to 78.6 kg (Luo 2021 Table S1 footnote equations for TVCL and TVV1).",
      source_name = "BW / WT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on linear CL (exponent -1.85), normalised to 37 g/L (Luo 2021 Table S1 footnote). Albumin is reported in g/L in Luo 2021 Table 2 (median 39, range 19-53).",
      source_name = "ALB"
    ),
    MM_NIGG = list(
      description = "Multiple-myeloma immunoglobulin type indicator: 1 = non-IgG MM, 0 = IgG MM",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (non-IgG MM) -- the typical-value CL is anchored to non-IgG MM; an IgG-MM patient receives the additive shift TPMMCL = 1 + 1.26 (Luo 2021 Table S1 footnote). Same model-level orientation as Xu_2020_daratumumab.",
      notes = "Additive shift on linear CL: TPMMCL = 1 + e_igg_cl * (1 - MM_NIGG), e_igg_cl = 1.26 (Luo 2021 Table S1 'TPMM on CL'). IgG MM patients have 2.26-fold the linear CL of non-IgG MM patients.",
      source_name = "TPMM (type of myeloma, IgG vs non-IgG)"
    ),
    SEXF = list(
      description = "Biological sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) -- Luo 2021 Table S1 footnote: 'SEXV1 is a shift factor of 1 for female and 1 [-] 0.105 for male'. The male indicator of the source equation is encoded as (1 - SEXF).",
      notes = "Additive fractional shift on V1: SEXV1 = 1 + e_male_vc * (1 - SEXF), e_male_vc = -0.105 (Luo 2021 Table S1 'Sex on V1', RSE 45.4%). Males have 10.5% lower V1 than females as printed; this direction is opposite to the sex effect of the intravenous daratumumab fit in Xu_2020_daratumumab (females 20.5% lower) -- see the vignette Assumptions section.",
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 742L,
    n_studies = 4L,
    n_observations = 5159L,
    age_range = "33-92 years",
    age_median = "67 years",
    weight_range = "28.6-147.6 kg",
    sex_female_pct = 43.3,
    race_ethnicity = c(White = 71, Black = 3, Asian = 11, Other = 15),
    disease_state = "Multiple myeloma: relapsed or refractory MM (PAVO Part 2, MMY1008 [Japanese patients], COLUMBA, PLEIADES D-Rd arm) and newly diagnosed MM (PLEIADES D-VRd and D-VMP arms). IgG myeloma 58%; ECOG 0/1/2/3: 36/52/12/<1%.",
    dose_range = "Daratumumab SC 1800 mg flat dose (co-formulated with rHuPH20, 3-5 min abdominal injection; n = 487: monotherapy 288, combination therapy 199) or IV 16 mg/kg (n = 255, COLUMBA). Monotherapy schedule: weekly in cycles 1-2, every 2 weeks in cycles 3-6, every 4 weeks thereafter (28-day cycles).",
    regions = "Multinational (Europe, North and South America, Asia-Pacific, Israel, Russia, Ukraine; Luo 2021 Table 1).",
    albumin_range = "19-53 g/L (median 39)",
    notes = "Baseline demographics from Luo 2021 Table 2 (Total column, N = 742). Reference covariates for the typical-value equations (Table S1 footnote): WT 78.6 kg, ALB 37 g/L, non-IgG MM, female. NONMEM 7.2 FOCE; objective function value -1759.836; condition number 42.8."
  )

  ini({
    # Structural typical-value PK parameters for the reference patient
    # (Luo 2021 Table S1: WT 78.6 kg, ALB 37 g/L, non-IgG MM, female).
    lcl <- log(0.00496); label("Linear (nonspecific) clearance CL at reference covariates (L/h)") # Luo 2021 Table S1: CL 0.00496 L/h, RSE 8.4%
    lvc <- log(5.25); label("Central volume of distribution V1 at reference covariates (L)") # Luo 2021 Table S1: V1 5.25 L, RSE 4.1%
    lvp <- log(3.78); label("Peripheral volume of distribution V2 (L)") # Luo 2021 Table S1: V2 3.78 L, RSE 7.0%
    lq <- log(0.00955); label("Intercompartmental clearance Q (L/h)") # Luo 2021 Table S1: Q 0.00955 L/h, RSE 8.0%

    # Saturable (target-mediated) elimination from the central compartment;
    # Vmax decays mono-exponentially from its baseline at first-order rate KDES.
    lvmax <- log(1.15); label("Baseline maximum velocity of the saturable clearance Vmax (mg/h)") # Luo 2021 Table S1: Vmax 1.15 mg/h, RSE 9.7%
    lkdes <- log(0.0000783); label("First-order rate of decrease of Vmax, KDES (1/h)") # Luo 2021 Table S1: KDES 0.0000783 1/h, RSE 33.3%
    lkm <- log(2.56); label("Michaelis-Menten constant Km (ug/mL)") # Luo 2021 Table S1: Km 2.56 ug/mL, RSE 16.0%

    # Subcutaneous absorption
    lka <- log(0.0117); label("First-order SC absorption rate constant Ka (1/h)") # Luo 2021 Table S1: Ka 0.0117 1/h, RSE 5.8%
    lfdepot <- log(0.689); label("SC bioavailability F1 (fraction)") # Luo 2021 Table S1: F1 0.689, RSE 2.7%

    # Covariate effects (Luo 2021 Table S1 footnote equations for TVCL and TVV1)
    e_wt_cl <- 1.24; label("Power exponent of WT/78.6 on linear CL (unitless)") # Luo 2021 Table S1: WT on CL 1.24, RSE 10.6%
    e_alb_cl <- -1.85; label("Power exponent of ALB/37 on linear CL (unitless)") # Luo 2021 Table S1: ALB on CL -1.85, RSE 12.1%
    e_igg_cl <- 1.26; label("Additive shift on linear CL for IgG MM vs the non-IgG MM reference (fraction)") # Luo 2021 Table S1: TPMM on CL 1.26, RSE 14.4%
    e_wt_vc <- 0.91; label("Power exponent of WT/78.6 on V1 (unitless)") # Luo 2021 Table S1: WT on V1 0.91, RSE 11.5%
    e_male_vc <- -0.105; label("Additive fractional shift on V1 for male sex vs the female reference (fraction)") # Luo 2021 Table S1: Sex on V1 -0.105, RSE 45.4%

    # Inter-individual variability. Table S1 reports IIV as %CV; converted to
    # log-normal variance via omega^2 = log(CV^2 + 1). No IIV on V2, Q, Km or F1.
    etalcl ~ 0.29607 # Luo 2021 Table S1: IIV CL 58.7% CV -> log(0.587^2 + 1)
    etalvc ~ 0.12766 # Luo 2021 Table S1: IIV V1 36.9% CV -> log(0.369^2 + 1)
    etalvmax ~ 0.37451 # Luo 2021 Table S1: IIV Vmax 67.4% CV -> log(0.674^2 + 1)
    etalkdes ~ 1.1406 # Luo 2021 Table S1: IIV KDES 145.9% CV -> log(1.459^2 + 1)
    etalka ~ 0.1225 # Luo 2021 Table S1: IIV Ka 36.1% CV -> log(0.361^2 + 1)

    # Residual error: additive on the log scale (Table S1 'ADD ERR (%CV)' 34.8,
    # 'additive error term on the log-scale') -> log-normal error in nlmixr2.
    expSd <- 0.348; label("Additive residual error on the log scale (SD)") # Luo 2021 Table S1: ADD ERR 34.8 %CV, RSE 0.3%
  })
  model({
    # Individual PK parameters (Luo 2021 Table S1 footnote):
    #   TVCL = 0.00496 * (WT/78.6)^1.24 * (ALB/37)^-1.85 * TPMMCL
    #     TPMMCL = 1 for non-IgG MM (MM_NIGG = 1), 1 + 1.26 for IgG MM (MM_NIGG = 0)
    #   TVV1 = 5.25 * (WT/78.6)^0.91 * SEXV1
    #     SEXV1 = 1 for female (SEXF = 1), 1 - 0.105 for male (SEXF = 0)
    cl <- exp(lcl + etalcl) *
      (WT / 78.6)^e_wt_cl *
      (ALB / 37)^e_alb_cl *
      (1 + e_igg_cl * (1 - MM_NIGG))
    vc <- exp(lvc + etalvc) *
      (WT / 78.6)^e_wt_vc *
      (1 + e_male_vc * (1 - SEXF))
    vp <- exp(lvp)
    q <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    kdes <- exp(lkdes + etalkdes)
    km <- exp(lkm)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot)

    # Time-varying Vmax: first-order decay from the baseline at rate KDES
    # (time since the first dose), representing CD38 target depletion.
    vmax_t <- vmax * exp(-kdes * t)

    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (cl / vc) * central -
      vmax_t * Cc / (km + Cc) -
      (q / vc) * central +
      (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # SC doses enter the depot with bioavailability F1; IV doses go to central.
    f(depot) <- fdepot

    Cc ~ lnorm(expSd)
  })
}
