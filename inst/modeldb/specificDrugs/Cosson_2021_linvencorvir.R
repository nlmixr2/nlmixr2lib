Cosson_2021_linvencorvir <- function() {
  description <- "Semiphysiological joint parent + metabolite population PK model for oral linvencorvir (RO7049389, RG7907), an HBV core protein allosteric modulator, and its active metabolite M5 in healthy volunteers and adults with chronic hepatitis B (Cosson 2021). The dose passes through a Savic transit-compartment absorption chain (non-integer number of transits) into a gut absorption compartment, from which drug enters the liver both passively (first order, ka) and by saturable OATP1B-mediated active uptake (Michaelis-Menten on the AMOUNT, vmax_uptake / km_uptake); the same saturable uptake also carries drug from the plasma compartment into the liver. Liver and plasma exchange passively (first-order k_liver_central / k_central_liver). Parent is eliminated from the liver by first-order formation of M5 (k_m5_form) and from plasma by a dose-dependent clearance representing the non-M5 metabolic pathway and possible biliary secretion. M5 has a one-compartment plasma disposition with volume equal to the parent plasma volume. Food raises relative bioavailability (fasted F = 0.439 relative to fed) and slows absorption (additive food effect on MTT). Asian ethnicity lowers vmax_uptake, cl and vc and raises k_m5_form; female sex lowers cl. All ODE states are in MILLIMOLES, as in the authors' NONMEM control stream; doses are entered in mg and converted inside model(). Because the amounts are small in mmol, solve with a tight absolute tolerance (for example rxSolve(..., atol = 1e-12)); the default leaves solver noise of order 1e-3 ng/mL at low troughs."
  reference <- paste(
    "Cosson V, Feng S, Jaminion F, Lemenuel-Diot A, Parrott N, Paehler A, Bo Q, Jin Y.",
    "How Semiphysiological Population Pharmacokinetic Modeling Incorporating Active",
    "Hepatic Uptake Supports Phase II Dose Selection of RO7049389, A Novel",
    "Anti-Hepatitis B Virus Drug.",
    "Clin Pharmacol Ther. 2021;109(4):1081-1091.",
    "doi:10.1002/cpt.2184.",
    "Parameter values from Table 1; model structure from the NONMEM control",
    "stream (Run71) in Supplementary Material S5.",
    sep = " "
  )
  vignette <- "Cosson_2021_linvencorvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_LINVENCORVIR_MG = list(
      description = "Nominal linvencorvir dose per administration in the treatment arm (mg)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on plasma clearance normalised to 200 mg: cl = theta_CL * (DOSE_LINVENCORVIR_MG / 200)^e_dose_cl. Source column TRT; the control stream computes NDOS = DOSE / (200 * BIO) with DOSE = TRT * BIO, so the bioavailability cancels and only the nominal dose enters. Constant within a treatment arm. In a data frame passed to rxSolve() place this column AFTER the event columns (id, time, evid, amt, cmt), because rxode2 can drop a dose-named covariate that precedes amt.",
      source_name = "TRT"
    ),
    FED = list(
      description = "Fed-vs-fasted dosing indicator, 1 = dosed with a standard meal, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Multiplies bioavailability by e_fasted_fdepot when fasted (BIO = 1*FOOD + THETA(8)*(1 - FOOD)) and adds e_fed_mtt to the mean transit time when fed (MTT = THETA(11) + FOOD*THETA(13)).",
      source_name = "FOOD"
    ),
    RACE_ASIAN = list(
      description = "Asian ethnicity indicator, 1 = Asian, 0 = non-Asian",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Power-style multiplicative effects on vmax_uptake, cl, vc (= vc_m5) and k_m5_form. The control stream sets ASIA = 1 when ETHN == 1; Supplementary Material S2 defines ETN as 0 for non-Asian and 1 for Asian.",
      source_name = "ETHN"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The control stream codes SEX with 1 = male and applies THETA(19)^(1 - SEX) to CL; 1 - SEX is SEXF. Direction confirmed by the text: 'Female subjects have a 63% lower CL compared with males' (THETA(19) = 0.369).",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "linvencorvir", units = "mmol", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "linvencorvir", units = "mmol", specimen = "tissue", verified = TRUE),
    central = list(analyte = "linvencorvir", units = "mmol", specimen = "plasma", verified = TRUE),
    central_m5 = list(analyte = "linvencorvir metabolite M5", units = "mmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 135L,
    n_studies = 3L,
    age_range = "18-60 years (median 29)",
    weight_range = "48.2-99.4 kg (median 71.7)",
    sex_female_pct = 13.3,
    race_ethnicity = c(Asian = 45.2, `non-Asian` = 54.8),
    disease_state = "105 healthy volunteers and 30 adults with chronic hepatitis B.",
    dose_range = "Single oral doses of 150-2,500 mg and multiple oral doses of 200-1,000 mg once daily or 200-800 mg twice daily, fasted or with a standard meal.",
    regions = "Global phase I/II (NCT02952924), China (NCT03570658) and a pitavastatin interaction study (NCT03717064).",
    notes = "Cosson 2021 Methods and Table S1. 2,994 linvencorvir and 1,983 M5 plasma concentrations. 7 of 105 healthy volunteers (all non-Asian) and 11 of 30 patients (all Asian) were female; 35 healthy volunteers and 26 patients were Asian. LLOQ 1.0 ng/mL for both analytes."
  )

  ini({
    # Saturable active hepatic uptake (OATP1B). Both constants are AMOUNTS in
    # mmol (Table 1 footnote b: molecular weight 598.69 g/mol); the control
    # stream divides VM by (KM + A) with A the compartment amount.
    lvmax_uptake <- log(0.123); label("Maximum liver uptake rate Vm (mmol/h)") # Cosson 2021 Table 1: theta1 Vm = 0.123 mMol/h, RSE 12.1%
    lkm_uptake <- log(0.00107); label("Michaelis-Menten constant of liver uptake Km (mmol)") # Cosson 2021 Table 1: theta2 Km = 0.00107 mMol, RSE 13.1%

    # Parent disposition.
    lcl <- log(71.7); label("Plasma clearance of parent at a 200 mg dose, non-Asian male (L/h)") # Cosson 2021 Table 1: theta3 CL = 71.7 L/h, RSE 13.9%
    lvc <- log(10.2); label("Plasma volume of parent V3 (= M5 volume V4), non-Asian (L)") # Cosson 2021 Table 1: theta4 V3 (=V4) = 10.2 L, RSE 9.19%
    lk_liver_central <- log(1.18); label("Liver-to-plasma rate constant K23 (1/h)") # Cosson 2021 Table 1: theta5 K23 = 1.18 1/h, RSE 5.68%
    lk_central_liver <- log(0.377); label("Plasma-to-liver passive rate constant K32 (1/h)") # Cosson 2021 Table 1: theta18 K32 = 0.377 1/h, RSE 18.1%

    # M5 formation and disposition.
    lk_m5_form <- log(0.0318); label("M5 formation rate constant from the liver KFM, non-Asian (1/h)") # Cosson 2021 Table 1: theta6 KFM = 0.0318 1/h, RSE 9.94%
    lcl_m5 <- log(1.27); label("M5 clearance CLM (L/h)") # Cosson 2021 Table 1: theta7 CLM = 1.27 L/h, RSE 7.60%

    # Absorption: Savic transit input into the absorption compartment.
    lka <- log(1.03); label("Passive first-order absorption rate constant into the liver Ka (1/h)") # Cosson 2021 Table 1: theta10 Ka = 1.03 1/h, RSE 11.7%
    lmtt <- log(0.360); label("Mean transit time, fasted (h)") # Cosson 2021 Table 1: theta11 MTT = 0.360 h, RSE 8.08%
    lntr <- log(7.24); label("Number of transit compartments N (unitless, non-integer)") # Cosson 2021 Table 1: theta12 N = 7.24, RSE 15.2%

    # Dose and food effects.
    e_dose_cl <- -0.684; label("Power exponent of (dose/200 mg) on plasma clearance (unitless)") # Cosson 2021 Table 1: theta9 Dose (mg) on CL = -0.684, RSE 8.35%
    e_fasted_fdepot <- 0.439; label("Relative bioavailability fasted vs fed (fraction)") # Cosson 2021 Table 1: theta8 rel BA in fasted state = 0.439, RSE 2.42%
    e_fed_mtt <- 0.805; label("Additive increase in MTT when fed (h)") # Cosson 2021 Table 1: theta13 Food on MTT = 0.805 h, RSE 7.45%

    # Ethnicity and sex effects (multiplicative, theta^indicator).
    e_race_asian_vmax_uptake <- 0.699; label("Asian vs non-Asian multiplier on Vm (ratio)") # Cosson 2021 Table 1: theta14 Asian on Vm = 0.699, RSE 20.0%
    e_race_asian_cl <- 0.460; label("Asian vs non-Asian multiplier on CL (ratio)") # Cosson 2021 Table 1: theta15 Asian on CL = 0.460, RSE 20.0%
    e_race_asian_vc <- 0.619; label("Asian vs non-Asian multiplier on V3 = V4 (ratio)") # Cosson 2021 Table 1: theta16 Asian on V3 = 0.619, RSE 13.7%
    e_race_asian_k_m5_form <- 1.62; label("Asian vs non-Asian multiplier on KFM (ratio)") # Cosson 2021 Table 1: theta17 Asian on KFM = 1.62, RSE 16.5%
    e_sexf_cl <- 0.369; label("Female vs male multiplier on CL (ratio)") # Cosson 2021 Table 1: theta19 Gender on CL = 0.369, RSE 24.1%

    # Between-subject variability: exponential, diagonal. The control stream's
    # ETA(11) on K32 is 0 FIX, so K32 carries no eta.
    etalvmax_uptake ~ 0.402 # Cosson 2021 Table 1: omega1^2 on Vm = 0.402 (63.4% CV)
    etalkm_uptake ~ 0.670 # Cosson 2021 Table 1: omega2^2 on Km = 0.670 (81.9% CV)
    etalcl ~ 0.542 # Cosson 2021 Table 1: omega3^2 on CL = 0.542 (73.6% CV)
    etalvc ~ 0.0853 # Cosson 2021 Table 1: omega4^2 on V3 = 0.0853 (29.2% CV)
    etalk_m5_form ~ 0.180 # Cosson 2021 Table 1: omega5^2 on KFM = 0.180 (42.4% CV)
    etalcl_m5 ~ 0.113 # Cosson 2021 Table 1: omega6^2 on CLM = 0.113 (33.6% CV)
    etalka ~ 0.277 # Cosson 2021 Table 1: omega7^2 on Ka = 0.277 (52.6% CV)
    etalmtt ~ 0.277 # Cosson 2021 Table 1: omega8^2 on MTT = 0.277 (52.6% CV)
    etalntr ~ 0.862 # Cosson 2021 Table 1: omega9^2 on N = 0.862 (92.8% CV)
    etalk_liver_central ~ 0.143 # Cosson 2021 Table 1: omega10^2 on K23 = 0.143 (37.8% CV)

    # Residual error: combined additive + proportional, Y = IPRED + IPRED*ERR1
    # + ERR3 (variances add, nlmixr2 combined2). Table 1 gives variances; the
    # SDs below are their square roots. The published proportional covariance
    # (0.0418, correlation 0.404) and additive covariance (34.8, correlation
    # 0.843) between the two analytes cannot be expressed in nlmixr2 and are
    # not carried.
    propSd <- 0.45277; label("Proportional residual SD, linvencorvir (fraction)") # Cosson 2021 Table 1: sigma1^2 = 0.205 (45.3%); sqrt(0.205) = 0.45277
    addSd <- 1.93391; label("Additive residual SD, linvencorvir (ng/mL)") # Cosson 2021 Table 1: sigma3^2 = 3.74 (SD 1.93 ng/mL); sqrt(3.74) = 1.93391
    propSd_m5 <- 0.22869; label("Proportional residual SD, M5 (fraction)") # Cosson 2021 Table 1: sigma2^2 = 0.0523 (22.8%); sqrt(0.0523) = 0.22869
    addSd_m5 <- 21.3776; label("Additive residual SD, M5 (ng/mL)") # Cosson 2021 Table 1: sigma4^2 = 457 (SD 21.4 ng/mL); sqrt(457) = 21.3776
  })

  model({
    # Molecular weights from the control stream $PK (MWRO, MWM5), g/mol.
    mw_parent <- 598.69
    mw_m5 <- 498.58

    # Individual parameters (control stream $PK, Run71).
    fdepot <- FED + e_fasted_fdepot * (1 - FED)
    vmax_uptake <- exp(lvmax_uptake + etalvmax_uptake) * e_race_asian_vmax_uptake^RACE_ASIAN
    km_uptake <- exp(lkm_uptake + etalkm_uptake)
    cl <- exp(lcl + etalcl) * (DOSE_LINVENCORVIR_MG / 200)^e_dose_cl * e_race_asian_cl^RACE_ASIAN * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * e_race_asian_vc^RACE_ASIAN
    k_liver_central <- exp(lk_liver_central + etalk_liver_central)
    k_central_liver <- exp(lk_central_liver)
    k_m5_form <- exp(lk_m5_form + etalk_m5_form) * e_race_asian_k_m5_form^RACE_ASIAN
    cl_m5 <- exp(lcl_m5 + etalcl_m5)
    vc_m5 <- vc
    ka <- exp(lka + etalka)
    mtt <- (exp(lmtt) + e_fed_mtt * FED) * exp(etalmtt)
    ntr <- exp(lntr + etalntr)

    kel <- cl / vc
    kel_m5 <- cl_m5 / vc_m5

    # Savic transit input, written exactly as the control stream $DES: the
    # most recent dose (mg, converted to mmol and scaled by the relative
    # bioavailability) enters the absorption compartment as a gamma-shaped
    # rate. LNFAC is the Stirling approximation of log(N!) used by the
    # authors; the 1e-5 terms are the stream's X = .00001 guards.
    ktr <- (ntr + 1) / mtt
    lnfac <- log(2.5066) + (ntr + 0.5) * log(ntr) - ntr + log(1 + 1 / (12 * ntr))
    tdose <- tad(depot)
    dose_mmol <- fdepot / mw_parent * podo(depot)
    ratein <- exp(log(dose_mmol + 1e-5) + log(ktr + 1e-5) + ntr * log(ktr * tdose + 1e-5) - ktr * tdose - lnfac)

    # Saturable active uptake into the liver, from the absorption compartment
    # (UPT1) and from plasma (UPT2); Michaelis-Menten on amounts.
    upt_depot <- vmax_uptake / (km_uptake + depot)
    upt_central <- vmax_uptake / (km_uptake + central)

    d/dt(depot) <- ratein - ka * depot - upt_depot * depot
    d/dt(liver) <- ka * depot + upt_depot * depot + upt_central * central + k_central_liver * central - k_liver_central * liver - k_m5_form * liver
    d/dt(central) <- k_liver_central * liver - upt_central * central - kel * central - k_central_liver * central
    d/dt(central_m5) <- k_m5_form * liver - kel_m5 * central_m5

    # The dose record only drives the transit input (control stream F1 = 0).
    f(depot) <- 0

    # Concentrations in ng/mL: mmol * g/mol = mg; mg/L * 1000 = ng/mL.
    Cc <- 1000 * mw_parent * central / vc
    Cc_m5 <- 1000 * mw_m5 * central_m5 / vc_m5

    Cc ~ add(addSd) + prop(propSd) + combined2()
    Cc_m5 ~ add(addSd_m5) + prop(propSd_m5) + combined2()
  })
}
