Pandey_2023_hydroxyurea <- function() {
  description <- paste(
    "Semi-mechanistic PK-PD model of oral hydroxyurea in children with",
    "sickle cell disease (HUSTLE trial, NCT00305175), used by Pandey 2023",
    "to study treatment nonadherence. PK: Savic transit-compartment",
    "absorption (Nt + 1 compartments, rate ktr) into a plasma compartment",
    "with first-order elimination ke; plasma volume from the Nadler",
    "blood-volume equation (or 70 mL/kg below 25 kg) times (1 - HCT/100).",
    "The pharmacodynamic driver is the daily average plasma concentration",
    "(piecewise constant within each dosing day, zero on a missed day).",
    "Erythropoiesis: stem-cell pool with Hill feedback from circulating",
    "red cells, three maturation compartments, and circulating red cells;",
    "hydroxyurea kills stem cells and prolongs red-cell lifespan",
    "(saturable). Mean cell volume (MCV) rises linearly with the daily",
    "average concentration. Mean cell fetal hemoglobin (HbF, pg) is",
    "induced through a Michaelis-Menten-formed intermediate acting via a",
    "Hill function. The model was fitted subject by subject (no population",
    "variability or residual error); typical values are the HUSTLE cohort",
    "central values of the parent paper (Pandey 2022)."
  )
  reference <- paste(
    "Pandey A, Raja R, Estepp JH, Ramkrishna D (2023). Leveraging",
    "mathematical modeling to analyze nonadherence for hydroxyurea therapy",
    "in sickle cell disease. CPT Pharmacometrics Syst Pharmacol",
    "12(6):748-758.",
    "doi:10.1002/psp4.12945.",
    "Model equations and parameter values from the parent paper:",
    "Pandey A, Estepp JH, Raja R, Kang G, Ramkrishna D (2022).",
    "Mathematical modeling of hydroxyurea therapy in individuals with",
    "sickle cell disease. Pharmaceutics 14(5):1065.",
    "doi:10.3390/pharmaceutics14051065.",
    sep = " "
  )
  vignette <- "Pandey_2023_hydroxyurea"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Mean cell volume, the hydroxyurea-derived HbF-inducing intermediate and
  # mean cell fetal hemoglobin are states specific to this model.
  paper_specific_compartments <- c("mcv", "hbf_intermediate", "hbf_mc")

  compartmentData <- list(
    depot = list(analyte = "hydroxyurea", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxyurea", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(
      analyte = "erythroid stem cells and early progenitors",
      units = "10^12 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor1 = list(
      analyte = "erythroid precursor cells",
      units = "10^12 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2 = list(
      analyte = "erythroid precursor cells",
      units = "10^12 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor3 = list(
      analyte = "erythroid precursor cells",
      units = "10^12 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    erythrocytes = list(
      analyte = "red blood cells",
      units = "10^12 cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    mcv = list(analyte = "mean cell volume", units = "fL", specimen = "blood cell", verified = TRUE),
    hbf_intermediate = list(
      analyte = "hydroxyurea metabolites that activate HbF production",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    hbf_mc = list(analyte = "fetal hemoglobin per red cell", units = "pg", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the blood-volume equation: Nadler formula for WT >= 25 kg,",
        "70 mL/kg for WT < 25 kg (Pandey 2023 Text S2 Eqs. S.5-S.6;",
        "Pandey 2022 Eqs. 3-4). Recomputed daily in the source analysis,",
        "so it may be supplied time-varying."
      ),
      source_name = "weight"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the Nadler blood-volume equation as (HT / 100)^3, i.e. height",
        "in metres cubed. Both source papers print the term as",
        "'height (cm)^3', which gives blood volumes of order 10^5 L; Nadler's",
        "equation uses metres, and the metre form reproduces the plasma",
        "volume reported in Pandey 2022 Table 2 (mean 1.77 L)."
      ),
      source_name = "height"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the male or female Nadler blood-volume coefficients (Pandey 2023 Text S2 Eq. S.5).",
      source_name = "sex"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Plasma volume Vp = Vb * (1 - HCT / 100) (Pandey 2023 Text S2 Eq. S.4).",
        "Observed hematocrit is used; it is not linked to the modelled red",
        "cell count or MCV. Recomputed daily in the source analysis, so it",
        "may be supplied time-varying. HUSTLE baseline mean 23.46%",
        "(Pandey 2022 Table 1)."
      ),
      source_name = "HCT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 1L,
    age_range = "children and adolescents; HUSTLE cohort starting age 10.12 (SD 4.72) years (Pandey 2022 Table 1)",
    sex_female_pct = 100 * 31 / 85,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Sickle cell disease on daily oral hydroxyurea (HUSTLE, NCT00305175,",
      "St. Jude Children's Research Hospital). The 40 participants analysed",
      "in Pandey 2023 showed declines in MCV or HbF suggestive of",
      "nonadherence; participants with transfusions, sparse or gapped data,",
      "or adherence adequately explained by pharmacy refill records were",
      "excluded."
    ),
    dose_range = "Daily oral dose computed per visit interval from pharmacy refill records (Pandey 2023 Text S1)",
    regions = "United States (Memphis, TN)",
    notes = paste(
      "The parent analysis (Pandey 2022) fitted the same equations to 85",
      "HUSTLE participants (54 male / 31 female; Table 1); the sex percentage",
      "above is that cohort's. Pandey 2023 refitted 40 of them with dosing",
      "profiles corrected for inferred nonadherence (Tables S1-S4 give the",
      "original and corrected per-participant estimates). Every parameter was",
      "estimated subject by subject by weighted least squares in MATLAB, so",
      "there is no population variability or residual error model. Baseline",
      "demographics and laboratory values: Pandey 2022 Table 1."
    )
  )

  ini({
    # PK (Pandey 2022 Table 2, 'Mean (SD)' over HUSTLE participants; the
    # per-participant PK estimates used in Pandey 2023 are not tabulated
    # there). Hourly rates converted to 1/day (x 24).
    lfdepot <- log(0.12); label("Oral bioavailability F (fraction)") # Pandey 2022 Table 2: F = 0.12 (0.04)
    lktr <- log(5.02 * 24); label("Transit-absorption rate constant ktr (1/day)") # Pandey 2022 Table 2: ktr = 5.02 (2.61) 1/h
    lnn <- log(1.14); label("Transit exponent Nt; the chain has Nt + 1 compartments (unitless)") # Pandey 2022 Table 2: Nt = 1.14 (1.08)
    lkel <- log(0.54 * 24); label("Plasma elimination rate constant ke (1/day)") # Pandey 2022 Table 2: ke = 0.54 (0.26) 1/h

    # Erythropoiesis and MCV (Pandey 2022 Table 3, 'Median' column;
    # per-participant values in Pandey 2023 Tables S1-S2).
    lpsi_rbc <- log(0.22); label("Red-cell feedback parameter Psi_e (10^12 cells/L)") # Pandey 2022 Table 3: Psi_e = 2.20 x 10^11 cells/L
    lgamma_rbc <- log(1.6); label("Steepness of the red-cell feedback gamma_e (unitless)") # Pandey 2022 Table 3: gamma_e = 1.6
    lkdsemax <- log(0.17); label("Maximum hydroxyurea-induced stem-cell death rate kdse,max (1/day)") # Pandey 2022 Table 3: kdse,max = 0.17 1/day
    lkdse50 <- log(0.43); label("Daily average concentration at half-maximal stem-cell death Kdse,50 (umol/L)") # Pandey 2022 Table 3: Kdse,50 = 0.43 uM
    lkte <- log(0.2); label("Erythroid maturation transit rate constant kte (1/day)") # Pandey 2022 Table 3: kte = 0.2 1/day
    lkdemax <- log(0.03); label("Red-cell death rate constant without drug kde,max (1/day)") # Pandey 2022 Table 3: kde,max = 0.03 1/day
    lkde50 <- log(310); label("Daily average concentration halving the red-cell death rate Kde,50 (umol/L)") # Pandey 2022 Table 3: Kde,50 = 310 uM
    lalpha_mcv <- log(0.37); label("Slope of new-cell MCV on daily average concentration alpha (fL per umol/L)") # Pandey 2022 Table 3: alpha = 0.37 fL/uM

    # Fetal hemoglobin (Pandey 2022 Table 5, 'Median' column;
    # per-participant values in Pandey 2023 Tables S3-S4).
    lkmet <- log(1.2); label("Maximum formation rate of the HbF-inducing intermediate kmet (umol/L/day)") # Pandey 2022 Table 5: kmet = 1.2 uM/day
    lkmmet <- log(27); label("Michaelis constant for intermediate formation Kmet (umol/L)") # Pandey 2022 Table 5: Kmet = 27 uM
    lkdi <- log(0.03); label("Degradation rate constant of the intermediate kdi (1/day)") # Pandey 2022 Table 5: kdi = 0.03 1/day
    lkaf <- log(1.5); label("Maximum hydroxyurea-activated HbF production rate kaf (pg/day)") # Pandey 2022 Table 5: kaf = 1.5 pg/day
    lkaf50 <- log(13); label("Intermediate concentration at half-maximal HbF activation Kaf (umol/L)") # Pandey 2022 Table 5: Kaf = 13 uM
    lhill_hbf <- log(3.2); label("Hill coefficient of HbF activation n (unitless)") # Pandey 2022 Table 5: n = 3.2
    lkdf <- log(0.07); label("HbF degradation rate constant kdf (1/day)") # Pandey 2022 Table 5: kdf = 0.07 1/day

    # Baselines. The source sets the initial conditions from each
    # participant's own pre-treatment data; the typical values are the
    # HUSTLE cohort baseline means.
    lbl_rbc <- log(2.777); label("Baseline red blood cell count Ne(0) (10^12 cells/L)") # Pandey 2022 Table 1: RBC at baseline 2.777 x 10^6 cells/uL
    lbl_mcv <- log(85.34); label("Baseline mean cell volume Vm0 (fL)") # Pandey 2022 Table 1: MCV at baseline 85.34 fL
    lbl_fm <- log(8.004 / 100 * 29.76); label("Baseline mean cell fetal hemoglobin Fm(0) (pg)") # Pandey 2022 Table 1 baseline HbF 8.004% x MCH 29.76 pg via Pandey 2023 Eq. S.1
  })

  model({
    # Individual parameters
    fdepot <- exp(lfdepot)
    ktr <- exp(lktr)
    nn <- exp(lnn)
    kel <- exp(lkel)

    psi_rbc <- exp(lpsi_rbc)
    gamma_rbc <- exp(lgamma_rbc)
    kdsemax <- exp(lkdsemax)
    kdse50 <- exp(lkdse50)
    kte <- exp(lkte)
    kdemax <- exp(lkdemax)
    kde50 <- exp(lkde50)
    alpha_mcv <- exp(lalpha_mcv)

    kmet <- exp(lkmet)
    kmmet <- exp(lkmmet)
    kdi <- exp(lkdi)
    kaf <- exp(lkaf)
    kaf50 <- exp(lkaf50)
    hill_hbf <- exp(lhill_hbf)
    kdf <- exp(lkdf)

    bl_rbc <- exp(lbl_rbc)
    bl_mcv <- exp(lbl_mcv)
    bl_fm <- exp(lbl_fm)

    # Plasma volume (Text S2 Eqs. S.4-S.6); Nadler height term in metres.
    vb <- (1 - SEXF) * (0.3669 * (HT / 100)^3 + 0.03219 * WT + 0.6041) +
      SEXF * (0.3561 * (HT / 100)^3 + 0.03308 * WT + 0.1833)
    if (WT < 25) {
      vb <- 70 / 1000 * WT
    }
    vp <- vb * (1 - HCT / 100)

    # Hydroxyurea molar mass (CH4N2O2, 76.05 g/mol) converts mg/L to umol/L;
    # the PD parameters are reported per umol/L.
    mw_hu <- 76.05

    # Transit absorption (Text S2 Eqs. S.2-S.3): the last of Nt + 1 transit
    # compartments holds F * D * (ktr t)^Nt exp(-ktr t) / Nt! after a dose and
    # empties into plasma at rate ktr. The dose is recorded on `depot`, which
    # only anchors podo()/tad(); daily doses are cleared well within a day, so
    # the last dose alone carries the input. Long daily-dosing solves need a
    # step cap below one day (e.g. rxSolve(..., hmax = 0.1)) so LSODA does not
    # step across a day's absorption peak and exhaust its step budget.
    tad_dose <- tad(depot)
    dose_amt <- podo(depot)
    transit_in <- 0
    Cavg <- 0
    if (!is.na(tad_dose) && !is.na(dose_amt) && dose_amt > 0) {
      if (tad_dose > 0) {
        transit_in <- exp(log(fdepot * dose_amt) + log(ktr) + nn * log(ktr * tad_dose) -
          ktr * tad_dose - lgamma(nn + 1))
      }
      # Daily average plasma concentration (Text S2: the PK model is
      # simulated each day with that day's dose and averaged). The day's
      # exposure F * D / (ke * Vp) spread over the 1-day dosing interval,
      # held constant for the day and zero on a day without a dose.
      if (tad_dose < 1) {
        Cavg <- fdepot * dose_amt / (kel * vp) / 1
      }
    }
    Cavg_uM <- Cavg * 1000 / mw_hu

    d/dt(depot) <- 0
    f(depot) <- 0
    d/dt(central) <- transit_in - kel * central
    Cc <- central / vp

    # Erythropoiesis (Text S2 Eqs. S.7-S.10). Pre-treatment steady state:
    # Nse(0) = Ne1(0) = Ne2(0) = Ne3(0) = kde,max * Ne(0) / kte, and kpe,max
    # follows from the Tables S1-S2 footnote formula.
    prol0 <- kdemax * bl_rbc / kte
    kpemax <- kte * prol0 * (psi_rbc^gamma_rbc + bl_rbc^gamma_rbc) / psi_rbc^gamma_rbc
    kpe <- kpemax * psi_rbc^gamma_rbc / (psi_rbc^gamma_rbc + erythrocytes^gamma_rbc)
    kdse <- kdsemax * Cavg_uM / (kdse50 + Cavg_uM)
    kde <- kdemax * (1 - Cavg_uM / (kde50 + Cavg_uM))

    d/dt(prol) <- kpe - kte * prol - kdse * prol
    d/dt(precursor1) <- kte * prol - kte * precursor1
    d/dt(precursor2) <- kte * precursor1 - kte * precursor2
    d/dt(precursor3) <- kte * precursor2 - kte * precursor3
    d/dt(erythrocytes) <- kte * precursor3 - kde * erythrocytes
    prol(0) <- prol0
    precursor1(0) <- prol0
    precursor2(0) <- prol0
    precursor3(0) <- prol0
    erythrocytes(0) <- bl_rbc

    # Mean cell volume (Text S2 Eq. S.12)
    d/dt(mcv) <- (alpha_mcv * Cavg_uM + bl_mcv - mcv) * kte * precursor3 / erythrocytes
    mcv(0) <- bl_mcv

    # Fetal hemoglobin (Text S2 Eqs. S.13-S.14); kbf = kdf * Fm(0) per the
    # Tables S3-S4 footnote.
    kbf <- kdf * bl_fm
    d/dt(hbf_intermediate) <- kmet * Cavg_uM / (kmmet + Cavg_uM) - kdi * hbf_intermediate
    d/dt(hbf_mc) <- kbf + kaf * hbf_intermediate^hill_hbf / (kaf50^hill_hbf + hbf_intermediate^hill_hbf) -
      kdf * hbf_mc
    hbf_mc(0) <- bl_fm

    RBC <- erythrocytes
    MCV <- mcv
  })
}
