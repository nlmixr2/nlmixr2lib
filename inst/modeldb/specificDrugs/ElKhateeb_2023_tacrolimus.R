ElKhateeb_2023_tacrolimus <- function() {
  description <- "Middle-out (PBPK-informed) one-compartment population PK model for oral tacrolimus whole-blood trough concentrations in adult renal transplant recipients with chronic kidney disease (base model, no covariate on intrinsic clearance). The only estimated structural parameter is the hepatic unbound intrinsic clearance CLint; hepatic blood flow, the blood:plasma ratio, the unbound fraction, renal clearance (filtration only), bioavailability (fa * fg * fh, well-stirred liver) and the steady-state blood volume of distribution are computed bottom-up from body weight, height, sex, hematocrit, serum albumin and eGFR at every occasion. Absorption is taken as complete within a dosing interval (ka >> kel), so the dose enters the central compartment directly. Between-subject and inter-occasion variability on CLint, where every trough measurement is its own occasion. Sibling covariate models: ElKhateeb_2023_tacrolimus_egfr and ElKhateeb_2023_tacrolimus_ckdstage."
  reference <- paste(
    "El-Khateeb E, Chinnadurai R, Al Qassabi J, Scotcher D, Darwich AS,",
    "Kalra PA, Rostami-Hodjegan A. Using Prior Knowledge on Systems Through",
    "PBPK to Gain Further Insight into Routine Clinical Data on Trough",
    "Concentrations: The Case of Tacrolimus in Chronic Kidney Disease.",
    "Ther Drug Monit. 2023;45(6):743-753. doi:10.1097/FTD.0000000000001108.",
    "Structural equations are Eq. 3-17; the bottom-up formulas and system",
    "constants are those of the authors' calculation spreadsheet",
    "(Supplementary File 1) and Monolix script (Supplementary File 2);",
    "estimates are the Table 1 'Base' row.",
    sep = " "
  )
  vignette <- "ElKhateeb_2023_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (recorded at every trough measurement). Enters the cardiac output CO = 187 * WT^0.81 mL/min (hepatic blood flow Qh = 0.227 * CO), the sex-specific blood volume (Eq. 14-15) and the tissue volumes, which are ICRP reference organ masses scaled linearly by WT / 60 kg (female) or WT / 73 kg (male).",
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only the total blood volume, Eq. 14 (male) and Eq. 15 (female), which take height in metres; the model divides HT by 100.",
      source_name = "HEIGHT"
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the female or male blood-volume equation (Eq. 14 vs 15) and the female (60 kg) or male (73 kg) ICRP reference organ masses. Supplementary File 1 codes Gender = 0 for female and 1 for male, so SEXF = 1 - Gender.",
      source_name = "Gender"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. The source records hematocrit as a volume fraction (median 0.36); the canonical column is in percent, so the model uses HCT / 100. Sets the blood:plasma ratio B/P = E:P * Hct + (1 - Hct) (Eq. 16 rearranged, E:P = 77.27) and splits blood volume into plasma (Eq. 12) and erythrocyte (Eq. 13) volumes.",
      source_name = "HCT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Scales the plasma unbound fraction from its healthy-volunteer value (0.013 at 45 g/L), Eq. 10.",
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI creatinine equation)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Computed by the authors with the 2009 CKD-EPI creatinine equation without a race coefficient (Supplementary File 1). The model converts it to L/h/1.73 m^2 (x 60/1000) and uses it, without de-normalising for body surface area, as the glomerular filtration rate in the renal blood clearance CLr,b = fub * eGFR (Eq. 7), as the Monolix script does.",
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on CLint",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- nesting level of the inter-occasion random effect",
      notes = "El-Khateeb 2023: 'Each Ctrough measurement was considered a separate occasion for each patient', so OCC increments with every trough (1167 occasions over 40 patients). Consumed by the ini() nesting syntax 'etaiov_lclint ~ ... | OCC', not by a model() expression. An event table needs at least two distinct OCC values.",
      source_name = "Occasion for each row"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 1L,
    n_observations = 1167L,
    age_range = "IQR 33.6-50.9 years at baseline",
    age_median = "41.6 years",
    weight_range = "40-141 kg over all measurements (IQR 65-85 kg; Supplementary File 1)",
    weight_median = "73 kg over all measurements (Supplementary File 1)",
    sex_female_pct = 42.5,
    race_ethnicity = c(White = 77.5, Other = 22.5),
    disease_state = "Adult renal transplant recipients more than 3 months post-transplant with chronic kidney disease of the graft: 20 with stable transplant function and 20 with failing transplant function (annual eGFR fall > 5 mL/min/1.73 m^2). eGFR 4.97-119.2 mL/min/1.73 m^2 over the follow-up (CKD stages 1-5); baseline eGFR median 45 (IQR 34.5-55.5).",
    dose_range = "Immediate-release tacrolimus (Prograf) capsules twice daily, 2-14 mg/day; troughs normalised to 5 mg every 12 h for the analysis",
    regions = "Salford Royal Hospital, United Kingdom (single centre)",
    notes = "Retrospective TDM cohort, median 29 steady-state whole-blood troughs per patient (1167 in total) over a median 3.3-year follow-up. Baseline demographics from El-Khateeb 2023 Results and Supplementary File 3 Table 2. Estimated in Monolix 2020 (SAEM)."
  )

  ini({
    lclint <- log(5886.7); label("Hepatic unbound intrinsic clearance CLint (L/h)") # El-Khateeb 2023 Table 1, Base row: CLintpop = 5886.7 L/h (RSE 5.22%); Supplementary File 3 Table 3 proportional row

    fa <- fixed(1); label("Fraction of the oral dose absorbed (unitless)") # El-Khateeb 2023 Methods below Eq. 5: fa = 1, held constant; Monolix script fa = 1
    egut <- fixed(0.86); label("Gut-wall extraction ratio, 1 - fg (unitless)") # El-Khateeb 2023 Methods below Eq. 5: fg = 0.14, held constant; Monolix script fg = 0.14; egut = 1 - 0.14

    etalclint ~ 0.1024 # El-Khateeb 2023 Table 1, Base row: BSV omega = 0.32 (Monolix SD); 0.32^2
    etaiov_lclint ~ 0.04 | OCC # El-Khateeb 2023 Table 1, Base row: IOV gamma = 0.2 (Monolix SD); 0.2^2; one occasion per trough

    # NOT FROM THE PAPER. The base model's proportional residual-error
    # magnitude is not reported anywhere in the article or supplements. The
    # value is the profile-likelihood estimate on the authors' own
    # per-measurement dataset (Supplementary File 1, 1167 troughs) with
    # CLintpop, omega and gamma held at the Table 1 Base values; the error is
    # proportional on the log-transformed observation, Monolix's proportional
    # error for log-normally distributed observations. At that value the
    # marginal -2LL is 7266 against the published 7289.1.
    propSd <- 0.101; label("Proportional residual error on log-transformed concentration (fraction)")
  })

  model({
    # Bottom-up system and drug parameters, recomputed at every occasion from
    # the formulas of the authors' spreadsheet (Supplementary File 1).
    sexm <- 1 - SEXF
    hct <- HCT / 100

    # Blood:plasma ratio, Eq. 16 rearranged; E:P = 77.27 is the value used in
    # Supplementary File 1 (the text rounds the fit to 77.4).
    ep_ratio <- 77.27
    bpr <- ep_ratio * hct + (1 - hct)

    # Plasma unbound fraction, Eq. 10: fup,healthy = 0.013, albumin 45 g/L.
    fup <- 1 / (1 + (ALB / 45) * (1 - 0.013) / 0.013)
    fub <- fup / bpr # Eq. 9

    # Total blood volume (L), Eq. 14 (male) and Eq. 15 (female), height in m.
    vb <- sexm * (0.3669 * (HT / 100)^3 + 0.03219 * WT + 0.6041) +
      SEXF * (0.3561 * (HT / 100)^3 + 0.03308 * WT + 0.1833)
    vplasma <- vb * (1 - hct) # Eq. 12
    verythro <- vb * hct # Eq. 13

    # Sum over tissues of Vt * Pt:p per kg body weight (L/kg). Organ mass (g)
    # at the ICRP reference weight (female 60 kg, male 73 kg) / tissue density
    # (g/L) * Rodgers Pt:p, in the order adipose, bone, brain, gut, heart,
    # kidney, liver, lung, pancreas, muscle, skin, spleen, other organs.
    # Masses and densities: Supplementary File 1; Pt:p: Supplementary File 3
    # Table 1 (Simcyp V18, Rodgers method).
    vtkp_female <- (22500 / 950 * 1.90 + 7800 / 1920 * 1.89 + 1300 / 1050 * 1.91 +
      1965 / 1040 * 1.56 + 250 / 1050 * 0.74 + 275 / 1050 * 0.96 +
      1400 / 1070 * 1.33 + 420 / 384 * 0.52 + 120 / 1050 * 1.37 +
      17500 / 1050 * 0.96 + 2300 / 1090 * 1.07 + 130 / 1040 * 0.97 +
      4040 / 1050 * 1.90) / 60
    vtkp_male <- (18200 / 950 * 1.90 + 10500 / 1920 * 1.89 + 1450 / 1050 * 1.91 +
      2110 / 1040 * 1.56 + 330 / 1050 * 0.74 + 310 / 1050 * 0.96 +
      1800 / 1070 * 1.33 + 500 / 384 * 0.52 + 140 / 1050 * 1.37 +
      29000 / 1050 * 0.96 + 3300 / 1090 * 1.07 + 150 / 1040 * 0.97 +
      5210 / 1050 * 1.90) / 73

    # Steady-state volume, Eq. 11, with every Pt:p multiplied by 10 to recover
    # the literature blood Vss (Methods, after Eq. 16).
    vss_plasma <- vplasma + verythro * ep_ratio +
      10 * WT * (SEXF * vtkp_female + sexm * vtkp_male)
    vss_b <- vss_plasma / bpr

    # Hepatic blood flow (L/h): CO = 187 * WT^0.81 mL/min, Qh = 0.227 * CO.
    co <- 187 * WT^0.81 * 0.001 * 60
    qh <- 0.227 * co

    clint <- exp(lclint + etalclint + etaiov_lclint)

    clh <- qh * fub * clint / (qh + fub * clint) # Eq. 8, hepatic blood clearance
    clr <- fub * CRCL * 60 / 1000 # Eq. 7, renal blood clearance
    eh <- clh / qh # Eq. 5, fh = 1 - eh
    fbio <- fa * (1 - egut) * (1 - eh) # Eq. 4
    kel <- (clh + clr) / vss_b # Eq. 6

    d / dt(central) <- -kel * central
    f(central) <- fbio

    # mg / L * 1000 = ng/mL whole blood
    Cc <- 1000 * central / vss_b
    Cc ~ lnorm(NA) + propT(propSd)
  })
}
