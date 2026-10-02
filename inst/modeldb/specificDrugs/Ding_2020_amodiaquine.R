# Joint parent-metabolite population PK model for oral amodiaquine and its
# active metabolite desethylamodiaquine in children aged 3-59 months
# receiving sulfadoxine-pyrimethamine + amodiaquine seasonal malaria
# chemoprevention in Magaria District, Niger (Ding 2020, Clin Pharmacol
# Ther 107(5):1179-1188; doi:10.1002/cpt.1707). The final NONMEM control
# stream is printed in Supplementary Material S1 ('NONMEM code for final
# population PK model'), which is the primary source for every value here.

Ding_2020_amodiaquine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral amodiaquine and",
    "its active CYP2C8-derived metabolite desethylamodiaquine in children",
    "aged 3-59 months receiving seasonal malaria chemoprevention",
    "(sulfadoxine-pyrimethamine + amodiaquine) in Niger (Ding 2020, n = 136).",
    "First-order absorption into a two-compartment amodiaquine disposition",
    "model with complete (molar) conversion to a three-compartment",
    "desethylamodiaquine disposition model. Allometric scaling on all",
    "apparent clearances (exponent 0.75) and volumes (exponent 1) at a",
    "reference weight of 10 kg, and a hyperbolic age-driven maturation",
    "function on both the amodiaquine and desethylamodiaquine clearances.",
    "Relative bioavailability is fixed to 1 with inter-individual",
    "variability. Predictions are capillary whole-blood (dried blood spot)",
    "concentrations in nmol/L.",
    sep = " "
  )
  reference <- paste(
    "Ding J, Coldiron ME, Assao B, Guindo O, Blessborn D, Winterberg M,",
    "Grais RF, Koscalova A, Langendorf C, Tarning J (2020). Adherence and",
    "population pharmacokinetic properties of amodiaquine when used for",
    "seasonal malaria chemoprevention in African children. Clinical",
    "Pharmacology & Therapeutics 107(5):1179-1188. doi:10.1002/cpt.1707.",
    sep = " "
  )
  vignette <- "Ding_2020_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  # What each ODE state holds. Doses are amodiaquine base in mg; the
  # desethylamodiaquine states hold mg of desethylamodiaquine base (the
  # molar conversion is applied on the formation flux in model()). Samples
  # were capillary blood dried on filter paper (Supplementary Material S1,
  # 'Blood samples').
  compartmentData <- list(
    depot = list(analyte = "amodiaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amodiaquine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "amodiaquine", units = "mg", specimen = "tissue", verified = TRUE),
    central_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Reference WT = 10 kg = predicted median body weight of",
        "the PK cohort (Supplementary Material S1, equations 5-6; Table 1",
        "footnote). Body weight was NOT measured in the study: it was",
        "PREDICTED from age with the pooled boys-and-girls regression of",
        "rural Niger children, WT (kg) = 5.46 + 0.162 * age (months)",
        "(Supplementary Material S1, equation 4). Sex-specific forms: boys",
        "5.75 + 0.160 * age, girls 5.13 + 0.162 * age (equations 2-3)."
      ),
      source_name = "WT"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the hyperbolic maturation function PNA / (pna50 + PNA) on",
        "both the amodiaquine and desethylamodiaquine clearances",
        "(Supplementary Material S1, equation 7 and $PK). The cohort spans",
        "3-59 months; the typical values in ini() are for full maturation."
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    WAZ = list(
      description = "Weight-for-age Z-score (WHO 2006 growth standard) from the predicted body weight",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a linear covariate on relative bioavailability (38.3%",
        "decrease per 1-unit WAZ decrease, dOFV = -7.5) but NOT retained",
        "because of poor precision (RSE 56%) (Results, 'Population PK of AQ",
        "and DEAQ in the PK cohort')."
      ),
      source_name = "WAZ"
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on all parameters by stepwise selection; not significant (Results).",
      source_name = "sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 136L,
    n_studies = 1L,
    n_observations = paste(
      "404 capillary blood samples; 367 amodiaquine samples after excluding",
      "37 collected >= 300 h post-dose (42 of them, 11.4%, below the LLOQ",
      "and omitted, M1) and 404 desethylamodiaquine samples (none below the",
      "LLOQ) (Supplementary Material S1, 'Population PK analysis')"
    ),
    age_range = "3-59 months (Methods); median 36 months (Methods, 'Assessment the predictive performance of the adherence method')",
    weight_range = "not measured; predicted from age, median 10 kg (Supplementary Material S1, equations 4-6)",
    sex_female_pct = 42.6,
    disease_state = "Healthy children receiving seasonal malaria chemoprevention (no acute malaria).",
    dose_range = paste(
      "Oral amodiaquine once daily for 3 days (57.5 mg for ages 3-11",
      "months, 153 mg for ages 12-59 months) with a single dose of",
      "sulfadoxine-pyrimethamine (250/12.5 or 500/25 mg) given with the",
      "first amodiaquine dose, all doses directly observed (Methods, 'Study",
      "design and drug regimen'). The 4 preceding monthly SMC rounds were",
      "assumed fully taken and included in the fit (Supplementary Material",
      "S1)."
    ),
    regions = "Magaria District, Niger",
    notes = paste(
      "PK cohort of the 2016 SMC case-control effectiveness study; 165",
      "enrolled, 12 with <= 1 blood draw and 17 who vomited within 60 min",
      "excluded. 78 of 136 (57.4%) male. Each child was sampled in 3 of 6",
      "pre-defined windows between 0 h and 35 days after the first dose.",
      "LLOQ 1.87 ng/mL (amodiaquine) and 2.95 ng/mL (desethylamodiaquine)."
    )
  )

  ini({
    # All structural values are Table 1 'NONMEM population estimates' and
    # the $THETA block of the Supplementary Material S1 control stream
    # (identical values). Typical values are for a 10 kg child with fully
    # mature metabolising enzymes (Table 1 footnote).
    lka <- log(2.85)
    label("First-order absorption rate constant ka (1/h)") # Table 1 ka = 2.85; S1 THETA(3)

    lcl <- log(101)
    label("Apparent amodiaquine clearance CL/F at WT = 10 kg, full maturation (L/h)") # Table 1 CL/F AQ = 101; S1 THETA(1)

    lvc <- log(314)
    label("Apparent amodiaquine central volume Vc/F at WT = 10 kg (L)") # Table 1 VC/F AQ = 314; S1 THETA(2)

    lq <- log(119)
    label("Apparent amodiaquine inter-compartmental clearance Q/F at WT = 10 kg (L/h)") # Table 1 Q/F AQ = 119; S1 THETA(5)

    lvp <- log(1820)
    label("Apparent amodiaquine peripheral volume Vp/F at WT = 10 kg (L)") # Table 1 VP/F AQ = 1,820; S1 THETA(6)

    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of amodiaquine (unitless)") # Table 1 F AQ = 100% fix; S1 THETA(4) = 1 FIX

    lcl_deaq <- log(2.33)
    label("Apparent desethylamodiaquine clearance CL/F at WT = 10 kg, full maturation (L/h)") # Table 1 CL/F DEAQ = 2.33; S1 THETA(7)

    lvc_deaq <- log(49.1)
    label("Apparent desethylamodiaquine central volume Vc/F at WT = 10 kg (L)") # Table 1 VC/F DEAQ = 49.1; S1 THETA(8)

    lq_deaq <- log(2.31)
    label("Apparent desethylamodiaquine inter-compartmental clearance Q1/F at WT = 10 kg (L/h)") # Table 1 Q1/F DEAQ = 2.31; S1 THETA(9) = 2.3

    lvp_deaq <- log(363)
    label("Apparent desethylamodiaquine first peripheral volume Vp1/F at WT = 10 kg (L)") # Table 1 VP1/F DEAQ = 363; S1 THETA(10)

    lq2_deaq <- log(4.34)
    label("Apparent desethylamodiaquine inter-compartmental clearance Q2/F at WT = 10 kg (L/h)") # Table 1 Q2/F DEAQ = 4.34; S1 THETA(11)

    lvp2_deaq <- log(98.1)
    label("Apparent desethylamodiaquine second peripheral volume Vp2/F at WT = 10 kg (L)") # Table 1 VP2/F DEAQ = 98.1; S1 THETA(12)

    lpna50 <- log(4.66)
    label("Postnatal age at 50% of amodiaquine clearance maturation (months)") # Table 1 Age50 on CL/F AQ = 4.66; S1 THETA(13)

    lpna50_deaq <- log(2.42)
    label("Postnatal age at 50% of desethylamodiaquine clearance maturation (months)") # Table 1 Age50 on CL/F DEAQ = 2.42; S1 THETA(14)

    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearances (unitless)") # S1 equation 5 and $PK (WT/10)**0.75

    e_wt_vc <- fixed(1)
    label("Allometric WT exponent on all apparent volumes (unitless)") # S1 equation 6 and $PK (WT/10)

    # IIV variances from the S1 $OMEGA block. The block is headed 'Initial
    # estimates' but its values are the final estimates: sqrt(omega)
    # reproduces every Table 1 'CV for IIV' entry (sqrt(0.0493) = 22.2%,
    # sqrt(0.646) = 80.4%, sqrt(3.0) = 173%, sqrt(0.141) = 37.5%,
    # sqrt(0.0232) = 15.2%, sqrt(0.466) = 68.3%). OMEGA(5), (6), (8), (9),
    # (11) and (12) are 0 FIX and are omitted.
    etalcl ~ 0.0493 # S1 OMEGA(1); Table 1 CV 22.2%
    etalvc ~ 0.646 # S1 OMEGA(2); Table 1 CV 80.4%
    etalka ~ 3.0 # S1 OMEGA(3); Table 1 CV 173%
    etalfdepot ~ 0.141 # S1 OMEGA(4); Table 1 CV 37.5%
    etalcl_deaq ~ 0.0232 # S1 OMEGA(7); Table 1 CV 15.2%
    etalvp_deaq ~ 0.466 # S1 OMEGA(10); Table 1 CV 68.3%

    # Additive residual error on natural-log molar concentrations
    # (S1 $ERROR, Y = IPRED + EPS); SD = sqrt(SIGMA).
    expSd <- 0.829
    label("Amodiaquine additive residual SD on the log scale (unitless)") # Table 1 sigma AQ = 0.829; S1 SIGMA(1,1) = 0.688

    expSd_deaq <- 0.204
    label("Desethylamodiaquine additive residual SD on the log scale (unitless)") # Table 1 sigma DEAQ = 0.204; S1 SIGMA(2,2) = 0.0417
  })

  model({
    # Molecular weights (g/mol) of the free bases, used to express the
    # 100% in vivo conversion of amodiaquine to desethylamodiaquine on a
    # molar basis and the molar output units. The source fitted molar
    # amounts (AMT in umol; S1 $INPUT) and scaled S = V/1000 to report
    # nmol/L; the same values are used by the other amodiaquine models in
    # this library.
    mwAQ <- 355.85
    mwDEAQ <- 327.81
    molarFactor <- mwDEAQ / mwAQ

    # Hyperbolic maturation of clearance on postnatal age in months
    # (S1 equation 7; $PK MF_AQ, MF_DEAQ).
    matAQ <- PNA / (exp(lpna50) + PNA)
    matDEAQ <- PNA / (exp(lpna50_deaq) + PNA)

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 10)^e_wt_cl * matAQ
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc
    q <- exp(lq) * (WT / 10)^e_wt_cl
    vp <- exp(lvp) * (WT / 10)^e_wt_vc

    cl_deaq <- exp(lcl_deaq + etalcl_deaq) * (WT / 10)^e_wt_cl * matDEAQ
    vc_deaq <- exp(lvc_deaq) * (WT / 10)^e_wt_vc
    q_deaq <- exp(lq_deaq) * (WT / 10)^e_wt_cl
    vp_deaq <- exp(lvp_deaq + etalvp_deaq) * (WT / 10)^e_wt_vc
    q2_deaq <- exp(lq2_deaq) * (WT / 10)^e_wt_cl
    vp2_deaq <- exp(lvp2_deaq) * (WT / 10)^e_wt_vc

    kel_aq <- cl / vc
    k12_aq <- q / vc
    k21_aq <- q / vp
    kel_deaq <- cl_deaq / vc_deaq
    k12_deaq <- q_deaq / vc_deaq
    k21_deaq <- q_deaq / vp_deaq
    k13_deaq <- q2_deaq / vc_deaq
    k31_deaq <- q2_deaq / vp2_deaq

    # ADVAN5 system of S1 $PK: CMT 1 dose, 2 AQ central, 4 AQ peripheral,
    # 3 DEAQ central, 5 and 6 DEAQ peripherals; K23 = CL/V2 routes all
    # amodiaquine elimination into desethylamodiaquine.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel_aq * central - k12_aq * central + k21_aq * peripheral1
    d/dt(peripheral1) <- k12_aq * central - k21_aq * peripheral1
    d/dt(central_deaq) <- molarFactor * kel_aq * central - kel_deaq * central_deaq -
      k12_deaq * central_deaq + k21_deaq * peripheral1_deaq -
      k13_deaq * central_deaq + k31_deaq * peripheral2_deaq
    d/dt(peripheral1_deaq) <- k12_deaq * central_deaq - k21_deaq * peripheral1_deaq
    d/dt(peripheral2_deaq) <- k13_deaq * central_deaq - k31_deaq * peripheral2_deaq

    # S1 $PK: F1 = THETA(4) * EXP(ETA(4)) on the dose compartment.
    fdepot <- exp(lfdepot + etalfdepot)
    f(depot) <- fdepot

    # mg / (g/mol) = mmol; mmol / L * 1e6 = nmol/L.
    Cc <- 1e6 * central / (mwAQ * vc)
    Cc_deaq <- 1e6 * central_deaq / (mwDEAQ * vc_deaq)

    Cc ~ lnorm(expSd)
    Cc_deaq ~ lnorm(expSd_deaq)
  })
}
