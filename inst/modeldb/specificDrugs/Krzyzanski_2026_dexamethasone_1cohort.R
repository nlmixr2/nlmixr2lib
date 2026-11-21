Krzyzanski_2026_dexamethasone_1cohort <- function() {
  description <- paste(
    "One-compartment population PK model with first-order intramuscular",
    "absorption for dexamethasone in parturient women given two 8 mg IM",
    "doses of dexamethasone phosphate 22 h apart before caesarean",
    "delivery, with the umbilical-cord venous (fetal) concentration",
    "proportional to the maternal concentration. Naive-pooled (1-cohort)",
    "fit: no between-subject variability; separate log-scale residual",
    "errors for maternal and cord samples. Clearance and volume are",
    "apparent (divided by F) and relative to the dexamethasone phosphate",
    "dose amount (Krzyzanski 2026)."
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
      "DEX.csv). Data re-analysed from Tsuei 1980. The same author's",
      "3-cohort fit of these data is Krzyzanski_2026_dexamethasone_3cohort."
    )
  )

  ini({
    # Table 3, column 'Estimate (%RSE) N = 1'. All four are natural-scale
    # THETAs in the control stream (Appendix 1 $THETA); no IIV ('For the
    # 1-cohort model IIV of log(CL/F) were fixed at 0', Data analysis).
    lcl <- log(34.7)
    label("Apparent clearance CL/F (L/h)") # Table 3, CL/F = 34.7 L/h (RSE 18%)
    lvc <- log(103)
    label("Apparent volume of distribution V/F (L)") # Table 3, V/F = 103 L (RSE 23%)
    lka <- log(2.92)
    label("First-order IM absorption rate constant ka (1/h)") # Table 3, ka = 2.92 1/h (RSE 127%)
    ratio_cord_central <- 0.329
    label("Fetal (umbilical-cord venous) to maternal plasma concentration ratio FMR (unitless)") # Table 3, FMR = 0.329 (RSE 7%); Eq 25

    # Residual error on log-transformed concentrations (Eqs 28-29; Appendix 1
    # Y = IPRED + EPS). Table 3 reports the variances sigma_m^2 and
    # sigma_f^2, so the log-scale SDs are their square roots.
    expSd <- sqrt(0.145)
    label("Log-scale residual SD, maternal plasma (unitless)") # Table 3, sigma_m^2 = 0.145 (RSE 30%)
    expSd_Ccord <- sqrt(0.161)
    label("Log-scale residual SD, umbilical-cord venous plasma (unitless)") # Table 3, sigma_f^2 = 0.161 (RSE 28%)
  })

  model({
    cl <- exp(lcl)
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
