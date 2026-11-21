Krzyzanski_2026_dexamethasone_3cohort <- function() {
  description <- paste(
    "One-compartment population PK model with first-order intramuscular",
    "absorption for dexamethasone in parturient women given two 8 mg IM",
    "doses of dexamethasone phosphate 22 h apart before caesarean",
    "delivery, with the umbilical-cord venous (fetal) concentration",
    "proportional to the maternal concentration. Gauss-Hermite 3-cohort",
    "fit: apparent clearance takes one of three values (23.1, 38.8 or",
    "65.2 L/h) drawn with fixed probabilities 1/6, 2/3, 1/6 through a",
    "mixture, the between-cohort SD on log CL/F standing in for",
    "between-subject variability; separate log-scale residual errors for",
    "maternal and cord samples. Clearance and volume are apparent",
    "(divided by F) and relative to the dexamethasone phosphate dose",
    "amount (Krzyzanski 2026)."
  )
  reference <- paste(
    "Krzyzanski W. Cohort population analysis of sparse data:",
    "Dexamethasone pharmacokinetics in mother and fetus based on blood",
    "sampling at birth. J Pharmacokinet Pharmacodyn. 2026;53(5):36.",
    "doi:10.1007/s10928-026-10048-5.",
    "Data from Tsuei SE, Petersen MC, Ashley JJ, McBride WG, Moore RG.",
    "Disposition of synthetic glucocorticoids II. Dexamethasone in",
    "parturient women. Clin Pharmacol Ther. 1980;28(1):88-98.",
    "doi:10.1038/clpt.1980.136."
  )
  vignette <- "Krzyzanski_2026_dexamethasone"
  # The source codes Dose = 8000 ug of dexamethasone PHOSPHATE per injection
  # (Methods 'Pharmacokinetic model'; Appendix 1 DOSE = 8000.0), not the
  # 6.67 mg dexamethasone base it is equivalent to (Discussion). ug / L gives
  # ng/mL, the unit of the reported concentrations (Table S2).
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  dosing <- c("depot")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "dexamethasone",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "dexamethasone",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    n_observations = 28L,
    sex_female_pct = 100,
    age_range = "not reported in the source paper",
    weight_range = "not reported in the source paper",
    race_ethnicity = "not reported in the source paper",
    disease_state = paste(
      "Near-term pregnant women delivering by caesarean section and given",
      "antenatal dexamethasone; their fetuses were sampled through the",
      "umbilical cord vein at birth."
    ),
    dose_range = paste(
      "8 mg dexamethasone phosphate (6.67 mg dexamethasone) IM at noon the",
      "day before delivery and 8 mg IM at 10 am on the day of delivery",
      "(22 h apart)."
    ),
    regions = "not reported in the source paper (data from Tsuei 1980)",
    notes = paste(
      "One maternal peripheral-venous and one umbilical-cord venous plasma",
      "sample per woman, both drawn at birth, 0.78-4.42 h after the second",
      "dose (22.78-26.42 h after the first; Table S2 and Appendix 1",
      "DEX.csv). Data re-analysed from Tsuei 1980. Individual cohort",
      "assignment (highest posterior probability): 2 women in the low,",
      "10 in the middle and 2 in the high clearance cohort (Results,",
      "Fig. 5). The same author's naive-pooled fit of these data is",
      "Krzyzanski_2026_dexamethasone_1cohort."
    )
  )

  ini({
    # Table 3, column 'Estimate (%RSE) N = 3'. Natural-scale THETAs in the
    # control stream (Appendix 1 $THETA 1-4).
    lcl <- log(38.8)
    label("Typical (middle-cohort) apparent clearance CL/F (L/h)") # Table 3, CL/F = 38.8 L/h (RSE 21%); THETA(1)
    lvc <- log(79.6)
    label("Apparent volume of distribution V/F (L)") # Table 3, V/F = 79.6 L (RSE 59%); THETA(2)
    lka <- log(1.63)
    label("First-order IM absorption rate constant ka (1/h)") # Table 3, ka = 1.63 1/h (RSE 147%); THETA(3)
    ratio_cord_central <- 0.329
    label("Fetal (umbilical-cord venous) to maternal plasma concentration ratio FMR (unitless)") # Table 3, FMR = 0.329 (RSE 7%); THETA(4); Eq 25

    # Between-cohort SD of log(CL/F) (Eq 26). Estimated as a fixed effect,
    # THETA(5) = OM1 in the control stream, not as an OMEGA.
    omega_cohort_lcl <- 0.30
    label("Between-cohort SD of log CL/F (unitless)") # Table 3, omega = 0.30 (RSE 11%); Results: omega^2 = 0.09

    # Mixing probabilities = Gauss-Hermite weights w_i / sqrt(pi) for N = 3
    # (Table S1b: w1 = w3 = 0.2954089751509, w2 = 1.181635900604), fixed in
    # the $MIX block of Appendix 1. They evaluate to 1/6 and 2/3; the third
    # cohort takes the remainder 1/6 inside mix().
    p_cohort1 <- fixed(0.2954089751509 / sqrt(pi))
    label("Probability of the low-clearance cohort (fraction)") # Appendix 1 $MIX P(1) = W1; Table S1b
    p_cohort2 <- fixed(1.181635900604 / sqrt(pi))
    label("Probability of the middle (typical) clearance cohort (fraction)") # Appendix 1 $MIX P(2) = W2; Table S1b

    # Appendix 1 $OMEGA 0.0 FIX: the cohort values carry a zero-variance eta.
    # Kept because rxode2's mix() requires the model to carry an eta.
    etalcl ~ fixed(0)

    # Residual error on log-transformed concentrations (Eqs 28-29; Appendix 1
    # Y = IPRED + EPS). Table 3 reports the variances sigma_m^2 and
    # sigma_f^2, so the log-scale SDs are their square roots.
    expSd <- sqrt(0.0611)
    label("Log-scale residual SD, maternal plasma (unitless)") # Table 3, sigma_m^2 = 0.0611 (RSE 54%)
    expSd_Ccord <- sqrt(0.0695)
    label("Log-scale residual SD, umbilical-cord venous plasma (unitless)") # Table 3, sigma_f^2 = 0.0695 (RSE 38%)
  })

  model({
    # Eq 26 / Appendix 1: CL = THETA(1) * exp(sqrt(2) * xi_i * OM1), with the
    # N = 3 Gauss-Hermite nodes xi = -1.224744871391589, 0, +1.224744871391589
    # (Table S1a; control stream ksi1-ksi3). sqrt(2) * xi = sqrt(3), so the
    # cohort values are 38.8 * exp(-/+ sqrt(3) * 0.30) = 23.1 / 65.2 L/h
    # (Results). mix() draws one cohort per subject with the fixed
    # probabilities; a 'mixest' data column (1, 2 or 3) pins the cohort.
    xi_cohort <- sqrt(2) * 1.224744871391589 * omega_cohort_lcl
    cl <- mix(
      exp(lcl - xi_cohort + etalcl), p_cohort1,
      exp(lcl + etalcl), p_cohort2,
      exp(lcl + xi_cohort + etalcl)
    )
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- cl / vc

    # Eqs 22-23. The two IM doses (t = 0 and t2 = 22 h) are dose records on
    # depot rather than the delta functions of Eq 22.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Eq 24 (maternal) and Eq 25 (fetal = FMR x maternal).
    Cc <- central / vc
    Ccord <- ratio_cord_central * Cc

    Cc ~ lnorm(expSd)
    Ccord ~ lnorm(expSd_Ccord)
  })
}
