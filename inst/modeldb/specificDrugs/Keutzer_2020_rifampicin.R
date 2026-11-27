Keutzer_2020_rifampicin <- function() {
  description <- "One-compartment population PK model for high-dose oral rifampicin in adults with pulmonary tuberculosis (Keutzer 2020 re-estimation of the Svensson 2018 HIGHRIF1 model with all inter-occasion variability removed): Erlang transit-compartment absorption, Michaelis-Menten clearance scaled by an auto-induced enzyme turnover compartment, fat-free-mass allometric scaling on Vmax (0.75) and volume (1.0), and saturable dose-dependent bioavailability above a 450 mg reference dose. This is the deliberately IOV-free comparator of the paper; the IOV-including original is Svensson_2018_rifampicin."
  reference <- paste(
    "Keutzer L, Simonsson USH. (2020). Individualized Dosing With High",
    "Inter-Occasion Variability Is Correctly Handled With Model-Informed",
    "Precision Dosing - Using Rifampicin as an Example.",
    "Front Pharmacol 11:794. doi:10.3389/fphar.2020.00794.",
    "Structural model re-estimated from Svensson RJ et al. (2018)",
    "Clin Pharmacol Ther 103(4):674-683 doi:10.1002/cpt.778;",
    "see modellib('Svensson_2018_rifampicin').",
    sep = " "
  )
  vignette <- "Keutzer_2020_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "Rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "Rifampicin", units = "mg", specimen = "plasma", verified = TRUE),
    enzyme = list(
      analyte = "Induced metabolising enzyme (relative amount, baseline 1)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on a 70 kg fat-free-mass reference with fixed exponents: Vmax scales (FFM/70)^0.75 and V scales (FFM/70)^1 (Keutzer 2020 Supplementary Appendix S1, $PK ALLMCL / ALLMV). In the Keutzer 2020 simulations FFM was derived from body weight with sex-specific power functions fitted to the HIGHRIF1 data (Eq. 1 male: 2.541 * WT^0.728; Eq. 2 female: 2.496 * WT^0.669).",
      source_name = "FFM"
    ),
    DOSE = list(
      description = "Administered rifampicin dose (mg) entering the saturable dose-dependent bioavailability function.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Bioavailability is 1 at the 450 mg reference dose and rises as 1 + femax * (DOSE - 450) / (fed50 + (DOSE - 450)) above it (Eq. 3; Supplementary Appendix S1 TVBIO). Keep DOSE equal to the dose of the current dosing record. Below 450 mg the published expression is undefined by the data (F450 was fixed to 1 for lack of lower-dose data); the model clamps the excess dose at zero so bioavailability stays 1 there.",
      source_name = "DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 1L,
    age_range = "adults",
    weight_range = "40.2-84.2 kg (mean 55.1 kg)",
    ffm_range = "28.5-57.9 kg (mean 40.0 kg)",
    sex_female_pct = 28.9,
    disease_state = "Drug-susceptible pulmonary tuberculosis.",
    dose_range = "Oral rifampicin 10, 20, 25, 30, and 35 mg/kg once daily (HIGHRIF1 dose-rising cohorts) with standard isoniazid, pyrazinamide and ethambutol from day 8.",
    regions = "Cape Town, South Africa (PanACEA HIGHRIF1, NCT01392911).",
    notes = "Keutzer 2020 re-estimated the Svensson 2018 population PK model on the original HIGHRIF1 data set after fixing all four inter-occasion variances to zero (Methods 'Evaluating the Predictive Performance of MIPD Incorporating or Ignoring IOV'); the OFV rose by 1,727.4 points. Demographics from Keutzer 2020 Table 1 (59 male, 24 female)."
  )

  ini({
    # Final estimates of the re-estimated POPPK model without IOV:
    # Keutzer 2020 Table 4, column 'Estimates model without IOV (RSE [%])'.
    # Supplementary Appendix S1 $THETA carries the same estimates to more
    # digits; those digits are used where they round to the Table 4 value.
    lvmax <- log(309.8)
    label("Maximal elimination rate Vmax at FFM 70 kg (mg/h)") # Table 4 Vmax 309.8 mg/h/70 kg (S1 THETA(1) 310)
    lkm <- log(15.8)
    label("Michaelis-Menten constant km (mg/L)") # Table 4 km 15.8 mg/L; S1 THETA(2) 15.8
    lvc <- log(93.6)
    label("Apparent volume of distribution at FFM 70 kg (L)") # Table 4 V 93.6 L/70 kg; S1 THETA(3) 93.6
    lka <- log(2.4)
    label("Absorption rate constant ka (1/h)") # Table 4 ka 2.4 1/h; S1 THETA(4) 2.4
    lmtt <- log(0.812)
    label("Mean transit time of the absorption chain MTT (h)") # Table 4 MTT 0.81 h; S1 THETA(8) 0.812
    lnn <- log(7.568)
    label("Number of transit compartments NN (unitless)") # Table 4 NN 7.6; S1 THETA(9) 7.568
    emax <- 1.219
    label("Maximal increase in enzyme production rate Emax (unitless)") # Table 4 Emax 1.2; S1 THETA(5) 1.219
    ec50 <- 0.05256
    label("Rifampicin concentration at half-maximal induction EC50 (mg/L)") # Table 4 EC50 0.053 mg/L; S1 THETA(6) 0.05256
    kenz <- 0.005328
    label("First-order enzyme formation and degradation rate constant kENZ (1/h)") # Table 4 kENZ 0.0053 1/h; S1 THETA(7) 0.005328
    femax <- 0.401
    label("Maximal increase in bioavailability above 450 mg Fmax (fraction)") # Table 4 Fmax 0.40; S1 THETA(10) 0.401
    fed50 <- 17.35
    label("Dose above 450 mg giving half of Fmax ED50 (mg)") # Table 4 ED50 17.4 mg; S1 THETA(11) 17.35

    # IIV: Table 4 lists IIV as % (SD of eta x 100), the same convention as its
    # 'original model' column, whose values equal sqrt(OMEGA) of Svensson 2018.
    # Supplementary Appendix S1 $OMEGA prints these SDs (and the km-Vmax
    # correlation 0.5957) rather than variances; variances below are SD^2 and
    # the covariance is corr * SD_km * SD_vmax = 0.5957 * 0.8062 * 0.5524.
    etalkm + etalvmax ~ c(0.64996, 0.26529, 0.30515) # Table 4 IIV km 80.6%, IIV Vmax 55.2%, correlation Vmax-km 59.6% (S1 0.8062 / 0.5957 / 0.5524)
    etalvc ~ 0.0074218 # Table 4 IIV V 8.6% (S1 0.08615)
    etalmtt ~ 0.34516 # Table 4 IIV MTT 58.8% (S1 0.5875)
    etalnn ~ 0.81866 # Table 4 IIV NN 90.5% (S1 0.9048)
    etalka ~ 0.082254 # Table 4 IIV ka 28.7% (S1 0.2868)

    # Residual error: additive on log-transformed concentrations, SD 0.569
    # (Table 4 'Additive error on log scale' 56.9%; S1 $SIGMA 0.57).
    expSd <- 0.569
    label("Residual SD of log-transformed concentrations (log-normal)") # Table 4 additive error on log scale 56.9%
  })

  model({
    # Individual parameters (Supplementary Appendix S1 $PK); fixed allometric
    # exponents 0.75 on Vmax and 1 on V with a 70 kg FFM reference.
    vmax <- exp(lvmax + etalvmax) * (FFM / 70)^0.75
    km <- exp(lkm + etalkm)
    vc <- exp(lvc + etalvc) * (FFM / 70)
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn + etalnn)

    # Dose-dependent bioavailability (Eq. 3, S1 TVBIO); excess dose clamped
    # at zero below the 450 mg reference.
    dose_xs <- max(DOSE - 450, 0)
    bio <- 1 * (1 + femax * dose_xs / (fed50 + dose_xs))

    Cc <- central / vc
    eff <- emax * Cc / (ec50 + Cc)

    # Erlang transit input (S1 $DES DADT(1)); the NONMEM code sets F1 = 0 so
    # the dose enters the depot only through the transit function.
    d / dt(depot) <- transit(nn, mtt, bio) - ka * depot
    d / dt(central) <- ka * depot - vmax / (km + Cc) / vc * central * enzyme
    d / dt(enzyme) <- kenz * (1 + eff) - kenz * enzyme

    # S1 $PK: A_0(2) = 0.0001, A_0(3) = 1
    central(0) <- 0.0001
    enzyme(0) <- 1
    f(depot) <- 0

    Cc ~ lnorm(expSd)
  })
}
