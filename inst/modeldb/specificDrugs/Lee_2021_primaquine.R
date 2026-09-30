Lee_2021_primaquine <- function() {
  description <- "Minimal physiologically based (semi-physiological liver) population PK model for oral primaquine and its carboxyprimaquine metabolite in healthy Korean adult men of normal weight and with obesity. First-order absorption with lag time into a well-stirred liver compartment (volume from body weight and height) exchanging with a one-compartment primaquine plasma pool at the fixed hepatic plasma flow; hepatic metabolism is split into a monoamine-oxidase intrinsic clearance that forms carboxyprimaquine (one-compartment disposition) and a CYP2D6 intrinsic clearance whose value rises exponentially with CYP2D6 activity score and body weight."
  reference <- paste(
    "Lee WY, Chae DW, Kim CO, Lee SE, Kwak YG, Yeom JS, Park KS.",
    "Population Pharmacokinetics of Primaquine in the Korean Population.",
    "Pharmaceutics. 2021;13(5):652. doi:10.3390/pharmaceutics13050652.",
    "The well-stirred liver structure and the definition of the hepatic",
    "extraction ratio, which Lee 2021 uses but does not print, are from",
    "the model Lee 2021 cites as its basis: Goncalves BP et al. Age, weight,",
    "and CYP2D6 genotype are major determinants of primaquine",
    "pharmacokinetics in African children. Antimicrob Agents Chemother.",
    "2017;61(5):e02590-16. doi:10.1128/AAC.02590-16.",
    sep = " "
  )
  vignette <- "Lee_2021_primaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Enters twice: (1) exponentially on the CYP2D6 intrinsic clearance, centred at 77.45 kg (Section 3.2.2 equation; Table 2 COVBW = 0.041 per kg), and (2) in the liver volume equation of Yu 2004, V_liver (mL) = 21.585 * WT^0.732 * HT^0.225 (Section 2.5.2). Cohort means 67.6 +/- 7.6 kg (normal-weight group) and 83.3 +/- 6.7 kg (obese group), Table 1.",
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Enters only the liver volume equation V_liver (mL) = 21.585 * WT^0.732 * HT^0.225 (Section 2.5.2, citing Yu 2004 for the Korean standard liver volume). Cohort means 175.4 +/- 4.6 cm (normal-weight) and 173.1 +/- 5.1 cm (obese), Table 1.",
      source_name = "height"
    ),
    CYP2D6 = list(
      description = "CYP2D6 activity score (activity score model A, sum of the two allele activity values)",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = "Genotype-derived CYP2D6 activity score assigned from both alleles of 17 genotyped variants (Section 2.3). Enters exponentially on the CYP2D6 intrinsic clearance, centred at 1.5, the modal value in the study (Section 3.2.2 equation; Table 2 COVAS = 1.254). Observed values 0.5, 1.0, 1.5 and 2.0 (Table 1); no subject had an activity score of 0 or above 2.0, so the exponential form is an extrapolation outside 0.5-2.0. No value transformation: the activity score is the column value.",
      source_name = "AS"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "primaquine", units = "mg", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "primaquine", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "primaquine", units = "mg", specimen = "plasma", verified = TRUE),
    central_cpq = list(analyte = "carboxyprimaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    age_range = "19-50 years (inclusion criterion)",
    age_median = "mean 29.1 years (normal-weight group) and 26.4 years (obese group)",
    weight_range = "mean 67.6 +/- 7.6 kg (normal-weight group) and 83.3 +/- 6.7 kg (obese group)",
    weight_median = "77.45 kg (covariate centring value)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy Korean adult men, G6PD-normal, randomised 12:12 into a normal body-weight group (BMI 18.5-24.9 kg/m^2) and an obese group (BMI >= 25.0 kg/m^2, Asian WHO criteria); overall inclusion BMI 18.6-31.2 kg/m^2.",
    dose_range = "Primaquine 15 mg orally once daily for 4 days, co-administered with hydroxychloroquine (800 mg, then 400 mg 10 h later on day 1, 400 mg on days 2 and 3); primaquine alone on day 4.",
    regions = "Republic of Korea (Severance Hospital, Seoul)",
    cyp2d6_activity_score = "Activity score 0.5: n = 1; 1.0: n = 6; 1.5: n = 12; 2.0: n = 5 (Table 1).",
    notes = "Baseline demographics from Lee 2021 Table 1. Plasma primaquine and carboxyprimaquine sampled on day 4 pre-dose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12 and 24 h post-dose (274 samples). LC-MS/MS, LLOQ 1 ng/mL for both analytes. Concentrations were modelled on the mass scale (ng/mL) without molar conversion (Section 2.5.2)."
  )

  ini({
    # All values are the final estimates of Lee 2021 Table 2 unless noted.
    # Clearances and volumes are for the typical subject (activity score
    # 1.5, 77.45 kg); no allometric scaling was retained (Section 3.2.1:
    # 'The theoretical allometric approach using body weight did not
    # improve the model').
    lka <- log(1.7); label("First-order absorption rate constant from depot into the liver ka (1/h)") # Table 2 KA = 1.7 1/h (RSE 12.7%)
    ltlag <- log(0.45); label("Absorption lag time (h)") # Table 2 ALAG1 = 0.45 h (RSE 1.7%)
    lvc <- log(142.2); label("Primaquine volume of distribution V3 (L)") # Table 2 V3 = 142.2 L (RSE 3.9%)
    lclint_mao <- log(19.1); label("Monoamine-oxidase intrinsic hepatic clearance of primaquine forming carboxyprimaquine CL_MAO (L/h)") # Table 2 CLMAO = 19.1 L/h (RSE 7.5%)
    lclint_cyp2d6 <- log(7.5); label("CYP2D6 intrinsic hepatic clearance of primaquine CL_CYP at activity score 1.5 and 77.45 kg (L/h)") # Table 2 CLCYP = 7.5 L/h (RSE 27.6%); Section 3.2.2 equation
    lcl_cpq <- log(1.3); label("Carboxyprimaquine clearance CLM (L/h)") # Table 2 CLM = 1.3 L/h (RSE 8.7%)
    lvc_cpq <- log(30.1); label("Carboxyprimaquine volume of distribution V4 (L)") # Table 2 V4 = 30.1 L (RSE 8.4%)

    # Exponential covariate effects on the CYP2D6 intrinsic clearance:
    # CLCYP = 7.5 * exp(COVAS * (AS - 1.5) + COVBW * (BW - 77.45)).
    # Section 3.2.2 prints the equation with the coefficients rounded to
    # 1.25 and 0.04; Table 2 carries the unrounded 1.254 and 0.041.
    e_cyp2d6_clint_cyp2d6 <- 1.254; label("Exponential effect of CYP2D6 activity score on CL_CYP (per unit activity score)") # Table 2 COVAS = 1.254 (RSE 11%)
    e_wt_clint_cyp2d6 <- 0.041; label("Exponential effect of body weight on CL_CYP (per kg)") # Table 2 COVBW = 0.041 (RSE 9.1%)

    # Liver plasma flow, fixed from 90 L/h liver blood flow and a
    # haematocrit of 45% (90 * 0.55).
    lqh <- fixed(log(49.5)); label("Liver plasma flow QH (L/h)") # Section 2.5.2 'The liver plasma flow rate was fixed to 49.5 L/h'

    # Inter-individual variability. Table 2 reports each omega as 'CV%';
    # the %RSE column is on the SD scale (the V3 and KA RSEs of 13.4% and
    # 16.1% are below the sqrt(2/24) = 28.9% floor a variance estimate
    # from 24 subjects cannot beat), so the CV% is read as 100 * omega
    # and the variance is (CV/100)^2. See the vignette Assumptions.
    etalvc ~ 0.027556 # Table 2 'BSV on V3' 16.6 CV% -> 0.166^2
    etalclint_mao ~ 0.051529 # Table 2 'BSV on CLMAO' 22.7 CV% -> 0.227^2
    etalclint_cyp2d6 ~ 0.304704 # Table 2 'BSV on CLCYP' 55.2 CV% -> 0.552^2
    etalka ~ 0.687241 # Table 2 'BSV on KA' 82.9 CV% -> 0.829^2
    etalcl_cpq ~ 0.04 # Table 2 'BSV on CLM' 20 CV% -> 0.20^2
    etaltlag ~ 0.004624 # Table 2 'BSV on ALAG' 6.8 CV% -> 0.068^2

    # Residual error, Eq. (2) Y = PRED * (1 + eps_pro) + eps_add. Table 2
    # lists only a proportional term for primaquine and a combined
    # proportional + additive term for carboxyprimaquine.
    propSd <- 0.179; label("Proportional residual error for primaquine (fraction)") # Table 2 sigma pro1 = 17.9 CV% (RSE 13.9%)
    propSd_cpq <- 0.157; label("Proportional residual error for carboxyprimaquine (fraction)") # Table 2 sigma pro2 = 15.7 CV% (RSE 14.2%)
    addSd_cpq <- 22.3; label("Additive residual error for carboxyprimaquine (ng/mL)") # Table 2 sigma add = 22.3 SD (RSE 31.7%)
  })

  model({
    # 1. Liver volume (Section 2.5.2, Yu 2004 Korean standard liver
    #    volume): V_liver (mL) = 21.585 * BW^0.732 * height^0.225,
    #    divided by 1000 to give litres.
    vliver <- 21.585 * WT^0.732 * HT^0.225 / 1000

    # 2. Individual parameters (Eq. 1, P_i = theta * exp(eta_i)).
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    vc <- exp(lvc + etalvc)
    clint_mao <- exp(lclint_mao + etalclint_mao)
    clint_cyp2d6 <- exp(lclint_cyp2d6 + etalclint_cyp2d6 +
      e_cyp2d6_clint_cyp2d6 * (CYP2D6 - 1.5) + e_wt_clint_cyp2d6 * (WT - 77.45))
    cl_cpq <- exp(lcl_cpq + etalcl_cpq)
    vc_cpq <- exp(lvc_cpq)
    qh <- exp(lqh)

    # 3. Well-stirred liver (Goncalves 2017, the model Lee 2021 builds
    #    on): E_H = CL_int / (QH + CL_int) with CL_int the sum of the two
    #    pathway intrinsic clearances, and each pathway's hepatic
    #    clearance is its share of E_H * QH. Liver outflow is QH in
    #    total: QH * (1 - E_H) returns to plasma (Figure 2B K23) and
    #    QH * E_H is metabolised. Figure 2B's K24 = CL_MAO / V2 and
    #    K20 = CL_CYP / V2 are these hepatic clearances. See the vignette
    #    for why the paper's own Table 3 excludes the reading in which the
    #    intrinsic clearances drain the liver in addition to
    #    QH * (1 - E_H).
    clint <- clint_mao + clint_cyp2d6
    eh <- clint / (qh + clint)
    clh_mao <- qh * clint_mao / (qh + clint)
    clh_cyp2d6 <- qh * clint_cyp2d6 / (qh + clint)

    # 4. ODEs (Figure 2B). K12 = KA; K23 = QH * (1 - EH) / V2;
    #    K32 = QH / V3; K24 = CL_MAO / V2; K20 = CL_CYP / V2;
    #    K40 = CLM / V4. Amounts in mg; carboxyprimaquine formed
    #    mass-for-mass because the authors did not convert to molar
    #    units (Section 2.5.2).
    d/dt(depot) <- -ka * depot
    d/dt(liver) <- ka * depot + qh * central / vc -
      qh * (1 - eh) * liver / vliver -
      (clh_mao + clh_cyp2d6) * liver / vliver
    d/dt(central) <- qh * (1 - eh) * liver / vliver - qh * central / vc
    d/dt(central_cpq) <- clh_mao * liver / vliver - cl_cpq * central_cpq / vc_cpq

    alag(depot) <- tlag

    # 5. Observations: mg/L * 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc_cpq <- 1000 * central_cpq / vc_cpq

    Cc ~ prop(propSd)
    Cc_cpq ~ add(addSd_cpq) + prop(propSd_cpq)
  })
}
