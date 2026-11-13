Saravanakumar_2026_plozalizumab_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model of plozalizumab (TAK-202, MLN1202), the humanized anti-CCR2 monoclonal antibody that is the parental antibody of the antibody-drug conjugate TAK-500, with parallel linear (kel) and Michaelis-Menten (CCR2 target-mediated) elimination from the central compartment (Saravanakumar 2026, Supplementary Table S3). Volume and Vmax are per kg body weight, so the model scales them linearly with WT and takes the dose in nmol per animal (1 mg = 6.725 nmol, MW 148,700 g/mol); concentrations in nM are then independent of weight. The authors used this model to show that the TAK-500 total-antibody profile overlays plozalizumab in monkeys (Figure 2C-E)."
  reference <- "Saravanakumar A, Deng R, Dong L, Haridas S, Xia CQ, Appleman VA, Abu-Yousif AO, Piatkov K, Wang H. Population-Based Modeling to Predict Human PK/PD of TAK-500, an Anti-CCR2 Antibody-Drug Conjugate for First-in-Human Study in Cancer Patients. Clin Transl Sci. 2026;19(4):e70522. doi:10.1111/cts.70522"
  vignette <- "Saravanakumar_2026_tak_500"
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Saravanakumar 2026 Equations 15-17 and Figure 2B: Ac/At = plozalizumab
  # amount in the central / peripheral compartment. Serum measured by ELISA.
  compartmentData <- list(
    central = list(analyte = "plozalizumab", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "plozalizumab", units = "nmol", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a fitted covariate effect: Supplementary Table S3 reports Vc in L/kg and Vmax in nmol/h/kg, so both are multiplied by WT (exponent 1). Because the dose was also given per kg, concentrations do not depend on the WT supplied. Individual monkey weights are not reported.",
      source_name = "per-kg parameterisation (Table S3 units)"
    )
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 18,
    n_studies = 1,
    sex_female_pct = 50,
    disease_state = "Healthy cynomolgus monkeys.",
    dose_range = "Single 30-min IV infusion of plozalizumab (TAK-202); Methods state 10, 30 and 100 mg/kg (3 male + 3 female per dose), while Results and Supplementary Table S2 report 1, 3, 10, 30 and 100 mg/kg.",
    notes = "Saravanakumar 2026 Methods 'Pharmacokinetic/Pharmacodynamic Study of TAK-202 (Anti-CCR2 mAb) in Cynomolgus Monkeys' and 'Two-Compartment PK Modeling of TAK-202 in Cynomolgus Monkeys'; serum sampled to day 36 and assayed by ELISA (LLOQ 0.4 ug/mL). The route is modelled as an IV bolus into the central compartment. n_subjects counts the three dose groups described in Methods; the Results text describes n = 6 per dose over five doses. Fitted in Phoenix NLME (FOCE-ELS)."
  )

  ini({
    # Saravanakumar 2026 Supplementary Table S3 'NLME estimates'.
    lvc <- log(0.051); label("Central volume of distribution per kg body weight (L/kg)") # Table S3 'TV (Vc), L/kg' 0.051 (%RSE 5)
    lkel <- log(0.0129); label("Linear elimination rate constant from the central compartment (1/h)") # Table S3 'TV (kel), 1/h' 0.0129 (%RSE 10)
    lvmax <- log(0.0397); label("Maximum Michaelis-Menten (CCR2-mediated) elimination rate per kg body weight (nmol/h/kg)") # Table S3 'TV (Vmax), nmole/h/kg' 0.0397 (%RSE 55)
    lkm <- log(13.9); label("Antibody concentration at half-maximal Michaelis-Menten elimination (nM)") # Table S3 'TV (Km), nM' 13.9 (%RSE 63)
    lk12 <- log(0.00801); label("Distribution rate constant from central to peripheral compartment, kpt (1/h)") # Table S3 'TV (kpt), 1/h' 0.00801 (%RSE 11)
    lk21 <- log(0.00811); label("Distribution rate constant from peripheral to central compartment, ktp (1/h)") # Table S3 'TV (ktp),1/h' 0.00811 (%RSE 11)

    # IIV -- Table S3 'Inter-subject variability, eta' column, read as the
    # Phoenix Omega variance x 100 (see the vignette for the evidence).
    etalkel ~ 0.88 # Table S3 'TV (kel)' eta 88 (shrinkage 10%)
    etalkm ~ 0.24 # Table S3 'TV (Km)' eta 24 (shrinkage 8%)

    propSd <- 0.30; label("Proportional residual error (fraction)") # Table S3 'Proportional error model (CV %)' 30
  })

  model({
    vc <- exp(lvc) * WT
    kel <- exp(lkel + etalkel)
    vmax <- exp(lvmax) * WT
    km <- exp(lkm + etalkm)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Equations 15-16: Michaelis-Menten elimination written on amounts.
    d/dt(central) <- k21 * peripheral1 - kel * central - vmax * central / (km * vc + central) - k12 * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equations 17-18
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
