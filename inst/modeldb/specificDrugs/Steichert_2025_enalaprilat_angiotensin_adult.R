Steichert_2025_enalaprilat_angiotensin_adult <- function() {
  description <- "Population PK/PD model of enalaprilat and the plasma angiotensin II / angiotensin I ratio in healthy adults after a single oral 20 mg dose of enalapril maleate (Steichert 2025, Pharmaceutics). Enalaprilat is modelled alone: Savic transit-compartment input (Monolix Mtt/Ktr parameterisation, number of transit stages n = Mtt * Ktr - 1) into a first-order absorption compartment whose rate constant lumps enalapril absorption and its hydrolysis to enalaprilat, then a two-compartment disposition model with linear elimination, reported as apparent (/F) parameters. An effect compartment (ke0) drives a full-inhibition sigmoid Imax model on the angiotensin II / angiotensin I ratio, a measure of in vivo ACE activity. IIV on ktr, Mtt, CL/F and V1/F (correlated), E0 and IC50; no covariates. Doses must be supplied as enalapril FREE BASE in ug: multiply a mass of enalapril maleate by 376.45/492.52 = 0.76433 (20 mg maleate = 15287 ug base). The PD data were restricted to 0-24 h after the dose."
  reference <- "Steichert M, Cawello W, Burckhardt BB, Suessenbach FK, Laeer S; LENA Consortium. Angiotensin II/Angiotensin I Ratio as a New Pharmacodynamic Parameter for Population Modelling in Healthy Adults and Children with Heart Failure Treated with Enalapril. Pharmaceutics. 2025;17(10):1345. doi:10.3390/pharmaceutics17101345"
  vignette <- "Steichert_2025_enalaprilat_angiotensin"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    depot = list(
      analyte = "enalapril (absorption compartment; absorption and conversion to enalaprilat lumped)",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "enalaprilat", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "enalaprilat", units = "ug", specimen = "tissue", verified = TRUE),
    effect = list(
      analyte = "enalaprilat (effect-site concentration)",
      units = "ug/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested with an exponential model on every parameter with a random effect in the COSSAC covariate search (Section 2.4.1); not retained because the p = 0.01 threshold was not reached (Section 3.3.1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as a power model scaled by the weighted mean in the COSSAC covariate search (Section 2.4.1); not retained (Section 3.3.1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested as a power model scaled by the weighted mean in the COSSAC covariate search (Section 2.4.1); not retained (Section 3.3.1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested as a power model scaled by the weighted mean in the COSSAC covariate search (Section 2.4.1); not retained (Section 3.3.1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 9,
    n_studies = 1,
    age_range = "19-30 years",
    age_median = "21 years",
    weight_range = "47-88 kg",
    weight_median = "74 kg",
    bmi_range = "18-25 kg/m^2",
    sex_female_pct = 33.3,
    disease_state = "Healthy volunteers, free of known organ disease",
    dose_range = "Single oral dose of 20 mg enalapril maleate (EnaHEXAL tablet) after a fast of at least 10 h; median 0.27 mg/kg (range 0.23-0.43 mg/kg)",
    regions = "Germany (Heinrich Heine University Duesseldorf, 2013; DRKS00037310)",
    n_observations = "288 serum enalaprilat concentrations (47 below the 0.70 ug/L LLOQ, handled as left-censored at 0, M4) and 99 angiotensin II / angiotensin I ratios (0-24 h; the 48 and 72 h ratios were excluded)",
    notes = "Demographics from Table 1 and Section 2.1. Angiotensin I by 125I radioimmunoassay and angiotensin II by solid-phase extraction plus ELISA; all samples drawn after at least 30 min supine."
  )

  ini({
    # Absorption: Savic transit chain (Monolix Mtt / Ktr parameterisation)
    # followed by a first-order absorption compartment (Table 2).
    lktr <- log(5.31); label("Transit rate constant ktr (1/h)") # Table 2: ktr = 5.31 1/h (RSE 25.7%)
    lmtt <- log(1.46); label("Mean transit time Mtt (h)") # Table 2: Mtt = 1.46 h (RSE 11.7%)
    lka <- log(1.19); label("Absorption rate constant ka, lumping enalapril absorption and conversion to enalaprilat (1/h)") # Table 2: ka = 1.19 1/h (RSE 8.6%)

    # Enalaprilat disposition, apparent parameters (Table 2).
    lcl <- log(36.39); label("Apparent clearance of enalaprilat CL/F (L/h)") # Table 2: CL/F = 36.39 L/h (RSE 10.2%)
    lvc <- log(223.71); label("Apparent central volume of enalaprilat V1/F (L)") # Table 2: V1/F = 223.71 L (RSE 15.3%)
    lq <- log(6.38); label("Apparent intercompartmental clearance of enalaprilat Q/F (L/h)") # Table 2: Q/F = 6.38 L/h (RSE 13.7%)
    lvp <- log(108.26); label("Apparent peripheral volume of enalaprilat V2/F (L)") # Table 2: V2/F = 108.26 L (RSE 27.2%)

    # Effect compartment and full-inhibition sigmoid Imax model (Table 2).
    lke0 <- log(0.48); label("Effect-compartment equilibration rate constant ke0 (1/h)") # Table 2: ke0 = 0.48 1/h (RSE 21.9%)
    lhill <- log(2.02); label("Sigmoidicity factor gamma (unitless)") # Table 2: gamma = 2.02 (RSE 14.9%)
    le0 <- log(0.043); label("Baseline angiotensin II / angiotensin I ratio E0 (unitless)") # Table 2: E0 = 0.043 (RSE 39.1%)
    lic50 <- log(30.01); label("Effect-site enalaprilat concentration giving half-maximal inhibition IC50 (ug/L)") # Table 2: IC50 = 30.01 ug/L (RSE 27.8%)

    # IIV. Table 2 reports each Monolix log-normal random effect as a CV%,
    # converted here with omega^2 = log(1 + CV^2); the vignette shows the
    # baseline VPC bands of Figure 4b favour this over omega = CV/100.
    etalktr ~ 0.42611 # Table 2: IIV ktr = 72.89 CV%
    etalmtt ~ 0.10757 # Table 2: IIV Mtt = 33.7 CV%
    # CL/F and V1/F block: variances from IIV CL/F = 30.17 CV% and
    # IIV V1/F = 46.61 CV%; covariance = 0.84 * sqrt(0.087116 * 0.19659),
    # using Table 2 correlation V1/F, CL/F = 0.84.
    etalcl + etalvc ~ c(0.087116, 0.10993, 0.19659)
    etale0 ~ 1.1031 # Table 2: IIV E0 = 141.9 CV%
    etalic50 ~ 0.48836 # Table 2: IIV IC50 = 79.35 CV%

    # Residual error (Table 2).
    propSd <- 0.072; label("Proportional residual error on enalaprilat (fraction)") # Table 2: PK proportional error = 0.072
    addSd <- 0.42; label("Additive residual error on enalaprilat (ug/L)") # Table 2: PK additive error = 0.42 ug/L
    propSd_ang2_ang1_ratio <- 0.41; label("Proportional residual error on the angiotensin II / angiotensin I ratio (fraction)") # Table 2: PD proportional error = 0.41
  })

  model({
    ktr <- exp(lktr + etalktr)
    mtt <- exp(lmtt + etalmtt)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0)
    hill <- exp(lhill)
    e0 <- exp(le0 + etale0)
    ic50 <- exp(lic50 + etalic50)

    # Monolix transit absorption: Mtt = (n + 1) / Ktr, so n = Mtt * Ktr - 1
    # (non-integer). rxode2's transit(n, mtt) uses (n + 1) / mtt internally,
    # which is Ktr again.
    ntr <- mtt * ktr - 1

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The dose enters depot only through the analytic transit input, so the
    # bolus itself is suppressed with f(depot) = 0.
    d/dt(depot) <- transit(ntr, mtt, 1) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    Cc <- central / vc

    # Effect compartment (Section 2.4.1 / 3.3.1), in concentration units.
    d/dt(effect) <- ke0 * (Cc - effect)
    Ce <- effect

    # Full-inhibition sigmoid Imax model (Imax fixed at 1; Section 3.3.1),
    # in the Monolix library form E = E0 * (1 - Ce^gamma / (Ce^gamma + IC50^gamma)).
    ang2_ang1_ratio <- e0 * (1 - Ce^hill / (Ce^hill + ic50^hill))

    # PK: Monolix combined error; the default 'combined1' form sums the SDs.
    Cc ~ add(addSd) + prop(propSd) + combined1()
    ang2_ang1_ratio ~ prop(propSd_ang2_ang1_ratio)
  })
}
