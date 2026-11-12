Dong_2024_adefovir <- function() {
  description <- paste(
    "One-compartment population PK model for adefovir after a single oral",
    "10 mg dose of the prodrug adefovir dipivoxil in healthy adults, given",
    "alone (reference period) and together with the four other components",
    "of a transporter phenotyping cocktail (metformin, sitagliptin,",
    "pitavastatin and digoxin; test period) (Dong 2024). First-order",
    "absorption with a lag time; elimination is the sum of a linear",
    "nonrenal clearance and a Michaelis-Menten renal arm that feeds a",
    "urinary-excretion compartment, with Vmax proportional to the absolute",
    "estimated GFR. The cocktail raises the apparent bioavailability and",
    "lowers the absorption rate constant. Inter-individual variability on V",
    "and Vmax; inter-occasion variability (two periods) on F, ka, lag time",
    "and V; proportional residual error on plasma concentration and on",
    "urinary amount."
  )
  reference <- paste(
    "Dong Q, Chen C, Taubert M, Bilal M, Kinzig M, Sorgel F,",
    "Scherf-Clavel O, Fuhr U, Dokos C. Understanding adefovir",
    "pharmacokinetics as a component of a transporter phenotyping cocktail.",
    "Eur J Clin Pharmacol. 2024;80(7):1069-1078.",
    "doi:10.1007/s00228-024-03673-x"
  )
  vignette <- "Dong_2024_adefovir"

  # Amounts are adefovir molar equivalents. The paper reports Vmax in
  # umol/h, Km in nmol/L, plasma concentration in nmol/L and urinary
  # amounts in umol (Table 3; Figure 1 axis), so the model carries amounts
  # in umol and volumes in L. A 10 mg adefovir dipivoxil dose (molecular
  # weight 501.47 g/mol) is 19.94 umol of adefovir. The reference-period F
  # of 0.590 is the published oral bioavailability of adefovir from
  # adefovir dipivoxil, so the dose enters as the prodrug's molar amount.
  units <- list(time = "h", dosing = "umol", concentration = "nmol/L")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Absolute estimated glomerular filtration rate (AGFR): the CKD-EPI",
        "2012 creatinine-cystatin C eGFR de-normalised to the subject's own",
        "Mosteller body surface area, in absolute mL/min (NOT indexed to",
        "1.73 m^2)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Proportional effect on Vmax, TVVmax = theta_Vmax * (AGFR / 105),",
        "where 105 mL/min is the cohort median (Results 'Final model with",
        "demographic and physiological covariates'; Supplemental Table 1s",
        "median 105, range 77.5-136 mL/min). A power form did not improve",
        "the fit."
      ),
      source_name = "AGFR"
    ),
    CONMED_DIGOXIN_METFORMIN_PITAVASTATIN_SITAGLIPTIN = list(
      description = paste(
        "Adefovir dipivoxil co-administered with metformin, sitagliptin,",
        "pitavastatin and digoxin (the test period of the transporter",
        "phenotyping cocktail study): 1 = test period, 0 = adefovir",
        "dipivoxil alone (reference period)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adefovir dipivoxil alone)",
      notes = paste(
        "The paper estimates period-specific F (FR = 0.590 fixed, FT =",
        "0.736) and ka (KaR = 5.18, KaT = 2.29 1/h) (Table 3, 'Impact of",
        "co-treatment'); encoded here as log-scale shifts from the",
        "reference period, so the typical values reproduce Table 3",
        "exactly. The paper attributes the effect to the cocktail as a",
        "whole, not to any one co-administered drug. Set per dosing",
        "occasion."
      ),
      source_name = "period (reference / test)"
    ),
    OCC = list(
      description = "Occasion (study period) index for inter-occasion variability: 1 or 2",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each subject received adefovir dipivoxil in two periods (reference",
        "and test); IOV on F, ka, lag time and V has one eta per period",
        "(Results 'Adefovir base model'). Records with OCC outside 1-2",
        "carry no IOV."
      ),
      source_name = "occasion"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on PK parameters and not retained (Results 'Final model with demographic and physiological covariates')."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    BSA = list(
      description = "Body surface area (Mosteller)",
      units = "m^2",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened and not retained; AGFR on Vmax was the only retained covariate."
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      notes = "Screened and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "adefovir dipivoxil", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "adefovir", units = "umol", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "adefovir", units = "umol", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "20-68 years (median 35)",
    weight_range = "55.5-94.9 kg (median 71.3)",
    sex_female_pct = 58.3,
    disease_state = "Healthy adults",
    dose_range = paste(
      "Single oral 10 mg adefovir dipivoxil, once alone (reference period)",
      "and once with metformin, sitagliptin, pitavastatin and digoxin",
      "(test period)"
    ),
    regions = "Germany",
    renal_function = "Absolute eGFR (CKD-EPI 2012, BSA-de-normalised) median 105 mL/min, range 77.5-136",
    n_observations = "1101 adefovir plasma concentrations and urinary amounts; 9 samples below the LLOQ removed",
    notes = paste(
      "Reanalysis of the transporter cocktail phenotyping trial NCT02743260",
      "(Cologne). Demographics from Supplemental Table 1s (10 male / 14",
      "female; BMI 24.5 (19.5-29.7) kg/m^2). Plasma and urine sampled up",
      "to 24 h; urine collected over 0-4, 4-8, 8-12, 12-16 and 16-24 h.",
      "Estimation in NONMEM 7.5."
    )
  )

  ini({
    # Structural parameters: Dong 2024 Table 3, 'Fixed effects'.
    ltlag <- log(0.122); label("Absorption lag time ALAG (h)") # Table 3, 'ALAG (h)' = 0.122 (RSE 17.2%)
    lvc <- log(235); label("Volume of distribution V (L)") # Table 3, 'V (L)' = 235 (RSE 14.1%)
    lcl_nonren <- log(11.8); label("Linear nonrenal clearance CLNR (L/h)") # Table 3, 'CLNR (L/h)' = 11.8 (RSE 25.7%)
    lkm <- log(0.170); label("Michaelis-Menten constant of renal elimination Km (umol/L)") # Table 3, 'Km (nmol/L)' = 170 = 0.170 umol/L (RSE 26.3%)
    lvmax <- log(2.40); label("Maximum renal elimination rate Vmax at AGFR 105 mL/min (umol/h)") # Table 3, 'Vmax (umol/h)' = 2.40 (RSE 23.5%)

    # Period-specific absorption: Dong 2024 Table 3, 'Impact of co-treatment'.
    lfdepot <- fixed(log(0.590)); label("Apparent bioavailability F, reference period (fraction)") # Table 3, 'FR' = 0.590 fixed (published value, ref 24)
    lka <- log(5.18); label("Absorption rate constant ka, reference period (1/h)") # Table 3, 'KaR (h-1)' = 5.18 (RSE 20.5%)
    e_conmed_digoxin_metformin_pitavastatin_sitagliptin_fdepot <- log(0.736 / 0.590); label("Cocktail co-administration shift on log F (unitless)") # Table 3, 'FT' = 0.736 (RSE 13.5%) vs FR 0.590: log(0.736/0.590)
    e_conmed_digoxin_metformin_pitavastatin_sitagliptin_ka <- log(2.29 / 5.18); label("Cocktail co-administration shift on log ka (unitless)") # Table 3, 'KaT (h-1)' = 2.29 (RSE 17.6%) vs KaR 5.18: log(2.29/5.18)

    # Inter-individual variability: Table 3 CV% converted with
    # omega^2 = log(CV^2 + 1) (Table 3 footnote).
    etalvc ~ 0.028822 # Table 3, IIV V = 17.1 CV%
    etalvmax ~ 0.035831 # Table 3, IIV Vmax = 19.1 CV%

    # Inter-occasion variability (two periods), same conversion.
    etaiov_fdepot_1 ~ 0.078095 # Table 3, IOV F = 28.5 CV%
    etaiov_fdepot_2 ~ fixed(0.078095) # same variance as occasion 1
    etaiov_ka_1 ~ 0.73314 # Table 3, IOV Ka = 104 CV%
    etaiov_ka_2 ~ fixed(0.73314) # same variance as occasion 1
    etaiov_tlag_1 ~ 0.46459 # Table 3, IOV ALAG = 76.9 CV%
    etaiov_tlag_2 ~ fixed(0.46459) # same variance as occasion 1
    etaiov_vc_1 ~ 0.0068654 # Table 3, IOV V = 8.30 CV%
    etaiov_vc_2 ~ fixed(0.0068654) # same variance as occasion 1

    # Residual error: Table 3 CV% converted with sigma = sqrt(log(CV^2 + 1))
    # (Table 3 footnote).
    propSd <- 0.094787; label("Proportional residual error, plasma (fraction)") # Table 3, 'Proportional error, plasma' = 9.50 CV%
    propSd_Aurine <- 0.27568; label("Proportional residual error, urinary amount (fraction)") # Table 3, 'Proportional error, urine' = 28.1 CV%
  })
  model({
    # Occasion indicators and inter-occasion etas.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # Individual parameters. Period effects on F and ka (Results 'Covariate
    # model for effects of co-treatment'); Vmax proportional to AGFR / 105
    # (Results 'Final model with demographic and physiological covariates').
    fdepot <- exp(lfdepot + e_conmed_digoxin_metformin_pitavastatin_sitagliptin_fdepot *
      CONMED_DIGOXIN_METFORMIN_PITAVASTATIN_SITAGLIPTIN + iov_fdepot)
    ka <- exp(lka + e_conmed_digoxin_metformin_pitavastatin_sitagliptin_ka *
      CONMED_DIGOXIN_METFORMIN_PITAVASTATIN_SITAGLIPTIN + iov_ka)
    tlag <- exp(ltlag + iov_tlag)
    vc <- exp(lvc + etalvc + iov_vc)
    cl_nonren <- exp(lcl_nonren)
    km <- exp(lkm)
    vmax <- exp(lvmax + etalvmax) * (CRCL / 105)

    # Table 1 'Model 2' (A1 = depot, A2 = central, A3 = urine); fu = 1
    # (Methods), so the Michaelis-Menten term acts on total concentration.
    conc <- central / vc
    renal_rate <- vmax * conc / (km + conc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl_nonren * conc - renal_rate
    d/dt(urine) <- renal_rate

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Plasma adefovir in nmol/L (umol/L x 1000). Aurine is the CUMULATIVE
    # amount excreted in urine (umol); the paper fits the amount excreted
    # within each collection interval, recovered by differencing Aurine
    # across the interval boundaries (or resetting urine at each boundary).
    Cc <- 1000 * conc
    Aurine <- urine

    Cc ~ prop(propSd)
    Aurine ~ prop(propSd_Aurine)
  })
}
