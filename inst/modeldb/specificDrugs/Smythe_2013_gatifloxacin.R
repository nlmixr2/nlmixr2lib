Smythe_2013_gatifloxacin <- function() {
  description <- paste(
    "One-compartment population PK model for oral gatifloxacin in adult",
    "African pulmonary tuberculosis patients co-administered rifampin,",
    "isoniazid, and pyrazinamide (Smythe 2013). Savic transit-compartment",
    "absorption (analytical form, N = 12.6, MTT = 0.65 h) feeds first-order",
    "absorption into a one-compartment disposition model. Apparent oral",
    "clearance is split into a GFR-mediated component scaled linearly with",
    "Cockcroft-Gault creatinine clearance and a non-GFR (other) component",
    "scaled allometrically with fat-free mass (FFM, Janmahasatian formula);",
    "apparent volume is scaled linearly with FFM. Age, sex, and HIV status",
    "modify the absorption rate constant. Relative bioavailability is fixed",
    "at 1 on the first dose and 11.7% lower at steady state.",
    sep = " "
  )
  reference <- paste(
    "Smythe W, Merle CS, Rustomjee R, Gninafon M, Bocar Lo M, Bah-Sow O,",
    "Olliaro PL, Lienhardt C, Horton J, Smith P, McIlleron H, Simonsson USH.",
    "Evaluation of Initial and Steady-State Gatifloxacin Pharmacokinetics",
    "and Dose in Pulmonary Tuberculosis Patients by Using Monte Carlo",
    "Simulations. Antimicrob Agents Chemother. 2013 Sep;57(9):4164-4171.",
    "doi:10.1128/AAC.00479-13.",
    sep = " "
  )
  vignette <- "Smythe_2013_gatifloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Raw Cockcroft-Gault creatinine clearance (not BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Smythe 2013 equation 2 reports CLCR (mL/min) computed by the",
        "Cockcroft-Gault formula with K = 1.23 for men and K = 1.04 for",
        "women. Cohort median CLCR = 94 mL/min (Table 1). Linear effect on",
        "the GFR-mediated component of apparent oral clearance:",
        "CL/F_GFR = (CL/F_GFR)_STD * CRCL / 94 (Smythe 2013 equation 4).",
        "Stored under the canonical CRCL with raw Cockcroft-Gault mL/min",
        "per inst/references/covariate-columns.md (the CLCR alias is the",
        "raw, non-BSA-normalised Cockcroft-Gault form).",
        sep = " "
      ),
      source_name        = "CLCR"
    ),
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Smythe 2013 equation 9 reports FFM derived from total body weight,",
        "height, and sex via the Janmahasatian semi-mechanistic formula",
        "(WHSMAX = 42.92, WHS50 = 30.93 in men; WHSMAX = 37.99,",
        "WHS50 = 35.98 in women, units kg/m^2). Cohort median FFM = 45 kg",
        "(Table 1). The reference patient is a typical 70-kg male with",
        "FFM = 55 kg (Table 2 footnote), so the normalisation denominator",
        "in the model is 55 kg. Allometric scaling on CL/F_other",
        "(exponent 0.75) and linear scaling on V/F (exponent 1).",
        sep = " "
      ),
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort range 18-58 years, median 29 years (Table 1). Fractional",
        "linear effect on ka centred at 29 years (Smythe 2013 Table 2",
        "footnote: 'AGEka, % increase in ka for every year change from",
        "the median AGE of 29 years').",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the cohort majority, 116/169)",
      notes              = paste(
        "Smythe 2013 Table 2 footnote: 'SEXka, % decrease in ka for female",
        "patients relative to male patients'. Encoded as a fractional",
        "decrease applied when SEXF = 1; reference category is male per",
        "the canonical SEXF orientation.",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    HIV_POS = list(
      description        = "HIV-1 antibody-positive status indicator: 1 = HIV+, 0 = HIV-",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative; 115/169 in the cohort)",
      notes              = paste(
        "Smythe 2013 Table 2 footnote: 'HIV+ -ka (%), % increase in ka for",
        "patients with HIV relative to patients without HIV'. 54/169",
        "subjects HIV-positive (51 of 99 in South Africa, 3 of 25 in",
        "Benin, none in Senegal or Guinea per Table 1). All HIV-positive",
        "subjects were antiretroviral-naive at enrolment.",
        sep = " "
      ),
      source_name        = "HIV"
    ),
    OCC = list(
      description        = "Sampling-occasion indicator: 1 = first-dose occasion, 2 = steady-state (day-28) occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Smythe 2013 sampled three plasma concentrations after the first",
        "dose (occasion 1) and three at steady state on approximately day",
        "28 (occasion 2). Decomposed inside model() into binary indicators",
        "oc1 and oc2 to multiplex (a) a single bioavailability shift",
        "(F is 11.7% lower at steady state, Table 2) and (b) the per-",
        "occasion IOV etas on log-CL, log-V, and log-MTT (each IOV",
        "variance is fixed equal across occasions, matching the source's",
        "single-variance IOV reporting per parameter).",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 169L,
    n_studies       = 1L,
    n_observations  = 954L,
    age_range       = "18-58 years (Table 1: median 29 years, IQR 24-35 years)",
    weight_range    = "35-80 kg (Table 1: median 55 kg, IQR 51-60 kg)",
    ffm_range       = "Table 1: median 45 kg, IQR 39-49 kg",
    crcl_range      = "Table 1: median 94 mL/min, IQR 81-110 mL/min (Cockcroft-Gault)",
    sex_female_pct  = 31.4,
    n_hiv_positive  = 54L,
    disease_state   = paste(
      "Newly diagnosed drug-sensitive pulmonary tuberculosis;",
      "antiretroviral-naive at enrolment in HIV-positive subjects."
    ),
    dose_range      = paste(
      "400 mg gatifloxacin (Lupin Pharmaceuticals) administered orally",
      "once daily for the first 2 months of treatment, irrespective of",
      "body weight, under directly observed therapy (DOT)."
    ),
    regions         = "Africa: South Africa (n=99), Senegal (n=26), Benin (n=25), Guinea (n=19)",
    co_medication   = paste(
      "Fixed-dose-combination rifampin 150 mg + isoniazid 75 mg +",
      "pyrazinamide 400 mg per tablet; 3 tablets if WT < 50 kg, 4 tablets",
      "if WT >= 50 kg. All four drugs given together once daily."
    ),
    notes           = paste(
      "OFLOTUB phase 3 randomised controlled trial",
      "(ClinicalTrials.gov NCT00216385); subset randomised to the",
      "4-month gatifloxacin-containing regimen. 12 of 954 observations",
      "were below the LLOQ of 0.1 ug/mL and replaced with LLOQ/2 in the",
      "source analysis; no more than one BLQ observation per subject in",
      "any individual absorption or elimination phase."
    )
  )

  ini({
    # Structural PK parameters. All values are FINAL ESTIMATES from
    # Smythe 2013 Table 2 (gatifloxacin final population PK model). The
    # parameter symbols (CL/F_GFR)_STD, (CL/F_Other)_STD, (V/F)_STD denote
    # the typical-value oral clearance / volume for the reference patient
    # (a typical 70-kg male with median CLCR 94 mL/min and FFM 55 kg).

    lclgfr  <- log(6.17);   label("Typical GFR-mediated oral clearance (CL/F_GFR)_STD at the reference patient (L/h)")    # Smythe 2013 Table 2: (CL/F_GFR)_STD = 6.17 L/h (RSE 9.7%)
    lclother<- log(5.11);   label("Typical non-GFR oral clearance (CL/F_Other)_STD at the reference patient (L/h)")        # Smythe 2013 Table 2: (CL/F_Other)_STD = 5.11 L/h (RSE 15.4%)
    lvc     <- log(141);    label("Typical apparent volume of distribution (V/F)_STD at the reference patient (L)")        # Smythe 2013 Table 2: (V/F)_STD = 141 L (RSE 2.7%)
    lka     <- log(4.13);   label("Absorption rate constant ka from the final transit compartment to central (1/h)")       # Smythe 2013 Table 2: ka = 4.13 1/h (RSE 13.5%)
    lmtt    <- log(0.65);   label("Mean transit time of the Savic transit-absorption chain (h)")                            # Smythe 2013 Table 2: MTT = 0.65 h (RSE 8.1%)
    lnn     <- log(12.6);   label("Number of transit compartments N (Savic chain, non-integer allowed) (unitless)")        # Smythe 2013 Table 2: N = 12.6 (RSE 19.7%)

    # Bioavailability: first-dose F is fixed at 1 (Table 2: 'F first dose
    # = 1 FIX'); steady-state F is reduced by 11.7% relative to first dose.
    lfdose1   <- fixed(log(1));    label("Bioavailability on the first dose, fixed at 1 (unitless)")                       # Smythe 2013 Table 2: F first dose = 1 FIX (not estimated)
    e_ss_fbio <- -0.117;            label("Fractional change in bioavailability at steady state vs first dose (unitless)") # Smythe 2013 Table 2: F steady-state = -11.7% change from F_first (RSE 17.4%)

    # Covariate effects on ka. Centred / fractional forms per Table 2 footnote.
    e_age_ka <-  0.032;     label("Fractional change in ka per +1 year of age above the cohort median 29 years (1/year)") # Smythe 2013 Table 2: AGE-ka = +3.2% per year from median 29 years (RSE 15.2%)
    e_sex_ka <- -0.548;     label("Fractional change in ka for females (SEXF = 1) relative to males (SEXF = 0) (unitless)")# Smythe 2013 Table 2: SEX-ka = -54.8% female vs male (RSE 10.7%)
    e_hiv_ka <-  0.619;     label("Fractional change in ka for HIV-positive (HIV_POS = 1) relative to HIV-negative (unitless)") # Smythe 2013 Table 2: HIV+ -ka = +61.9% (RSE 38.4%)

    # IIV. Source reports IIV as CV%; converted to log-normal variance
    # omega^2 = log(CV^2 + 1). Smythe 2013 fits IIV exponentially for all
    # parameters; the final model retains IIV on CL/F (combined GFR + other)
    # and on V/F (Table 2 'IIV' block). The two etas are kept diagonal
    # here -- the paper notes 'Covariance between parameters was also
    # tested' (Methods) but the final model in Table 2 reports no
    # covariance between CL and V.
    etalcl  ~ 0.10336       # Smythe 2013 Table 2: IIV CL/F = 33.0% CV (RSE 7.7%) -> log(0.33^2 + 1) = 0.10336
    etalvc  ~ 0.04769       # Smythe 2013 Table 2: IIV V/F  = 22.1% CV (RSE 10.9%) -> log(0.221^2 + 1) = 0.04769

    # IOV. The source reports a single IOV variance per parameter
    # (Table 2 'IOV' block). With only two sampling occasions (first-dose
    # and steady-state) in the design, we encode IOV via two per-occasion
    # etas with the same variance (the second fixed equal to the first),
    # matching the NONMEM '$OMEGA BLOCK(1) + SAME' idiom used in the
    # sister Wilkins 2008 rifampicin model (DDMODEL00000280). The OCC
    # column (1 or 2) is decomposed into binary indicators inside model().
    etaiov_cl_1  ~ 0.10336              # Smythe 2013 Table 2: IOV CL/F = 33.0% CV (RSE 5.7%) -> log(0.33^2 + 1) = 0.10336 (occasion 1)
    etaiov_cl_2  ~ fix(0.10336)         # IOV on log-CL, occasion 2; variance fixed equal to occasion 1 per source's single-variance IOV reporting
    etaiov_vc_1  ~ 0.01727              # Smythe 2013 Table 2: IOV V/F  = 13.2% CV (RSE 13.9%) -> log(0.132^2 + 1) = 0.01727 (occasion 1)
    etaiov_vc_2  ~ fix(0.01727)         # IOV on log-V,  occasion 2; variance fixed equal to occasion 1
    etaiov_mtt_1 ~ 0.18365              # Smythe 2013 Table 2: IOV MTT  = 44.9% CV (RSE 12.3%) -> log(0.449^2 + 1) = 0.18365 (occasion 1)
    etaiov_mtt_2 ~ fix(0.18365)         # IOV on log-MTT, occasion 2; variance fixed equal to occasion 1

    # Residual error. Combined additive + proportional on the linear
    # concentration scale (mg/L). The source paper also reports a
    # separately estimated 'predose additive error' of 0.0418 ug/mL
    # (Table 2) that is applied only to predose-occasion observations
    # following an unobserved dose; this third error component is omitted
    # from the packaged model and documented in the vignette Errata.
    addSd  <- 0.341;        label("Additive residual error (mg/L)")                                  # Smythe 2013 Table 2: additive error = 0.341 ug/mL = 0.341 mg/L (RSE 5.1%)
    propSd <- 0.0735;       label("Proportional residual error (fraction)")                          # Smythe 2013 Table 2: proportional error = 7.35% (RSE 12.5%)
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # 2. Per-occasion IOV etas. Each parameter carries one of two etas
    # depending on the current occasion; the IOV variances are fixed equal
    # across occasions (see ini()).
    iov_cl  <- oc1 * etaiov_cl_1  + oc2 * etaiov_cl_2
    iov_vc  <- oc1 * etaiov_vc_1  + oc2 * etaiov_vc_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # 3. Individual structural PK parameters.
    # GFR-mediated CL/F: linear in CLCR (raw Cockcroft-Gault mL/min), no
    # IIV by itself in the source (IIV is reported on total CL/F, applied
    # to the 'other' arm here so total CL/F variability matches the
    # paper's etalcl on CL/F overall).
    cl_gfr   <- exp(lclgfr) * (CRCL / 94)
    # Non-GFR CL/F: allometric on FFM with the paper's reference FFM 55 kg
    # (the FFM of the reference 70-kg male per Table 2 footnote). The IIV
    # eta etalcl and the IOV eta iov_cl ride on this arm so total CL/F
    # variability inherits the published 33% CV IIV and 33% CV IOV.
    cl_other <- exp(lclother + etalcl + iov_cl) * (FFM / 55)^0.75
    cl       <- cl_gfr + cl_other
    # Apparent volume: linear in FFM, reference 55 kg. IIV + IOV.
    vc       <- exp(lvc + etalvc + iov_vc) * (FFM / 55)
    # Absorption rate constant ka with three covariate effects.
    ka       <- exp(lka) *
                  (1 + e_age_ka * (AGE - 29)) *
                  (1 + e_sex_ka * SEXF) *
                  (1 + e_hiv_ka * HIV_POS)
    # Mean transit time and chain length. ktr derived as (N + 1)/MTT per
    # Smythe 2013 equation 1.
    mtt      <- exp(lmtt + iov_mtt)
    nn       <- exp(lnn)
    # ktr is the per-transit-compartment first-order rate; consumed
    # internally by rxode2's analytical transit() function.
    ktr      <- (nn + 1) / mtt

    # 4. Bioavailability. First-dose occasion gets the fixed F = 1; the
    # steady-state occasion gets F = 1 * (1 - 0.117) = 0.883. Encoded as
    # a fractional shift on the typical-value F driven by the oc2 indicator.
    fbio     <- exp(lfdose1) * (1 + e_ss_fbio * oc2)

    # 5. ODE system. Analytical Savic transit-compartment chain (rxode2
    # transit(n, mtt, bio) emits the gamma-PDF input rate accommodating
    # non-integer n via lgamma(n+1)); the depot then drains to central via
    # first-order ka. The dose lands on depot but the bolus is suppressed
    # via f(depot) <- 0 so the entire dose enters via the transit-chain
    # input rate (the same pattern as Wilkins 2008 / Vinnard 2017 /
    # vanderWalt 2013 transit-absorption models in nlmixr2lib).
    kel <- cl / vc
    d/dt(depot)   <- transit(nn, mtt, fbio) - ka * depot
    d/dt(central) <-                          ka * depot - kel * central
    f(depot) <- 0

    # 6. Observation and error.
    # Dose units mg, V units L -> Cc units mg/L (= ug/mL); the source
    # paper reports concentrations in ug/mL throughout (Table 1, Figures 1
    # and 2), which are numerically identical to mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
