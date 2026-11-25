Tham_2017_insulin <- function() {
  description <- paste(
    "Joint one-compartment population PK model of subcutaneous insulins in",
    "healthy adults (Tham 2017): three independent disposition systems, each",
    "with its own apparent clearance and volume -- insulin glargine; human",
    "insulin (regular U-100, regular U-500 and NPH share CL/F and V/F); and",
    "insulin lispro (lispro and insulin lispro protamine suspension, ILPS,",
    "share CL/F and V/F). Each formulation has its own first-order absorption",
    "depot; NPH and lispro carry an absorption lag; glargine and lispro add a",
    "parallel zero-order input of part of the dose directly into the central",
    "compartment. Body weight on V/F (linear for glargine, power for human",
    "insulin and lispro) and on the endogenous insulin baseline; per-injection",
    "dose on glargine bioavailability, regular U-100 ka and human-insulin CL/F.",
    "Premixes (Mix 70/30, Mix 50/50, Mix 75/25) are dosed by component."
  )
  reference <- paste(
    "Tham LS, Schneck K, Ertekin A, Reviriego J. Modeling Pharmacokinetic",
    "Profiles of Insulin Regimens to Enhance Understanding of Subcutaneous",
    "Insulin Regimens. J Clin Pharmacol. 2017;57(9):1126-1137.",
    "doi:10.1002/jcph.899. Equations from the Supplement Methods",
    "('Pharmacokinetic Model Equations')."
  )
  vignette <- "Tham_2017_insulin"
  units <- list(
    time = "h",
    dosing = "U (international units of insulin)",
    concentration = "pmol/L"
  )

  compartmentData <- list(
    depot_glar = list(
      analyte = "insulin glargine",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    central_glar = list(analyte = "insulin glargine", units = "U", specimen = "serum", verified = TRUE),
    depot_regu100 = list(
      analyte = "human regular insulin U-100",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_regu500 = list(
      analyte = "human regular insulin U-500",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_nph = list(
      analyte = "human insulin isophane suspension (NPH)",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    central_hum = list(
      analyte = "exogenous human insulin (regular and NPH)",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    depot_lis = list(
      analyte = "insulin lispro",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_ilps = list(
      analyte = "insulin lispro protamine suspension (ILPS)",
      units = "U",
      specimen = "administration site",
      verified = TRUE
    ),
    central_lis = list(
      analyte = "insulin lispro (from lispro and ILPS)",
      units = "U",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NA_character_,
      notes = paste(
        "Centred at 71.4 kg (linear) for glargine V/F and at 80 kg (power)",
        "for human-insulin V/F, lispro V/F and the endogenous insulin",
        "baseline (Supplement Methods equations)."
      ),
      source_name = "BW"
    ),
    DOSE_GLARGINE_U = list(
      description = "Insulin glargine dose of the current injection",
      units = "U",
      type = "continuous",
      reference_category = NA_character_,
      notes = paste(
        "Drives glargine relative bioavailability, (Dose/35)^theta1.",
        "Time-varying across injections; set it to the units of each",
        "glargine injection from that dose onwards. A value of 0 (no",
        "glargine given) switches the effect off (reference 35 U used)."
      ),
      source_name = "Dose_Glar"
    ),
    DOSE_REGULARU100_U = list(
      description = "Human regular insulin U-100 dose of the current injection",
      units = "U",
      type = "continuous",
      reference_category = NA_character_,
      notes = paste(
        "Drives the regular U-100 absorption rate constant, (Dose/10)^theta2.",
        "For Mix 70/30 it is the regular (30%) component only. A value of 0",
        "switches the effect off (reference 10 U used)."
      ),
      source_name = "Dose_RegU100"
    ),
    DOSE_HUMANINSULIN_U = list(
      description = "Total human insulin (regular U-100 + regular U-500 + NPH) units of the current injection(s)",
      units = "U",
      type = "continuous",
      reference_category = NA_character_,
      notes = paste(
        "Drives the shared human-insulin CL/F, (Dose/10)^theta6. Because",
        "regular and NPH insulin share one central compartment with one CL/F,",
        "the dose term is the total human-insulin units given at that time",
        "(for Mix 70/30, the whole premix dose). A value of 0 switches the",
        "effect off (reference 10 U used)."
      ),
      source_name = "Dose_Reg/NPH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 317L,
    n_studies = 16L,
    age_range = "18-62 years",
    weight_range = "42.8-124 kg",
    bmi_range = "16.6-38.6 kg/m^2",
    disease_state = "Healthy volunteers without diabetes (PK and euglycemic clamp studies)",
    dose_range = "Single SC doses: regular U-100 3-100 U, regular U-500 50-100 U, NPH 21-42 U, Mix 70/30 19-31 U, lispro 19-24 U, Mix 50/50 16-25 U, Mix 75/25 17-24 U, glargine 14-58 U",
    regions = "Not reported (Eli Lilly and Company clinical pharmacology studies)",
    notes = paste(
      "Tham 2017 Table 1. Subjects by formulation: regular U-100 30, regular",
      "U-500 12, NPH 24, Mix 70/30 18, lispro 8, Mix 50/50 8, Mix 75/25 8,",
      "glargine (Basaglar) 209; n_subjects is the sum of the Table 1 rows and",
      "may double-count crossover subjects. Sex distribution not reported.",
      "Externally validated (Figure 3) in obese non-diabetic, T1DM and T2DM",
      "subjects."
    )
  )

  ini({
    # ---------------- Insulin glargine (Table 2, 'Insulin Glargine' column)
    lka_glar <- log(0.0830)
    label("Glargine first-order absorption rate constant (1/h)") # Table 2 Ka glargine 0.0830 (RSE 4.25%)
    lvc_glar <- log(768)
    label("Glargine apparent volume of distribution at 71.4 kg (L)") # Table 2 V glargine 768 L (RSE 4.41%)
    lcl_glar <- log(74.5)
    label("Glargine apparent clearance (L/h)") # Table 2 CL glargine 74.5 L/h (RSE 2.47%)
    lfdepot_glar <- fixed(log(1))
    label("Glargine relative bioavailability at 35 U (unitless)") # Table 2 FGlar 1 (Fixed)
    e_dose_fdepot_glar <- -0.300
    label("Power exponent of glargine dose (U/35) on bioavailability (unitless)") # Table 2 theta1 -0.300 (RSE 21.4%)
    logitffo_glar <- qlogis(0.822)
    label("Logit of the fraction of the glargine dose absorbed first-order (unitless)") # Table 2 fraction first-order glargine 0.822
    ld1_glar <- log(0.612)
    label("Duration of glargine zero-order input (h)") # Table 2 duration of 0-order absorption glargine 0.612 h (RSE 4.53%)
    e_wt_vc_glar <- 0.00728
    label("Linear slope of body weight (kg - 71.4) on glargine V/F (1/kg)") # Table 2 BW on V glargine 0.00728 (RSE 32.4%)

    # ---------------- Human insulin: regular U-100, regular U-500, NPH
    lka_regu100 <- log(0.67)
    label("Regular U-100 first-order absorption rate constant at 10 U (1/h)") # Table 2 Ka regular U100 0.67 (RSE 11.1%)
    e_dose_ka_regu100 <- -0.275
    label("Power exponent of regular U-100 dose (U/10) on its ka (unitless)") # Table 2 dose on ka regular U100 -0.275 (RSE 24.9%)
    lka_regu500 <- log(0.185)
    label("Regular U-500 first-order absorption rate constant (1/h)") # Table 2 Ka regular U500 0.185 (RSE 9.08%)
    lka_nph <- log(0.129)
    label("NPH first-order absorption rate constant (1/h)") # Table 2 Ka NPH 0.129 (RSE 14.3%)
    ltlag_nph <- log(0.378)
    label("NPH absorption lag time (h)") # Table 2 ALag NPH 0.378 h (RSE 28.0%)
    lvc_hum <- log(178)
    label("Human insulin (regular/NPH) apparent volume of distribution at 80 kg (L)") # Table 2 V regular/NPH 178 L (RSE 15.9%)
    e_wt_vc_hum <- 1.62
    label("Power exponent of body weight (WT/80) on human insulin V/F (unitless)") # Table 2 BW on V regular U100 1.62 (RSE 27.0%)
    lcl_hum <- log(127)
    label("Human insulin (regular/NPH) apparent clearance at 10 U (L/h)") # Table 2 CL regular/NPH 127 L/h (RSE 6.36%)
    e_dose_cl_hum <- -0.282
    label("Power exponent of human insulin dose (U/10) on CL/F (unitless)") # Table 2 dose on CL regular/NPH -0.282 (RSE 13.6%)

    # ---------------- Insulin lispro and ILPS
    lka_lis <- log(0.989)
    label("Lispro first-order absorption rate constant (1/h)") # Table 2 Ka lispro 0.989 (RSE 14.9%)
    ltlag_lis <- log(0.265)
    label("Lispro first-order absorption lag time (h)") # Table 2 ALag lispro 0.265 h (RSE 15.4%)
    logitffo_lis <- qlogis(0.729)
    label("Logit of the fraction of the lispro dose absorbed first-order (unitless)") # Table 2 fraction first-order lispro 0.729 (RSE 8.49%)
    ld1_lis <- log(1.06)
    label("Duration of lispro zero-order input (h)") # Table 2 duration of 0-order absorption lispro 1.06 h (RSE 0.206%)
    lka_ilps <- log(0.0365)
    label("ILPS first-order absorption rate constant (1/h)") # Table 2 Ka lispro protamine suspension 0.0365 (RSE 16.5%)
    lvc_lis <- log(43.0)
    label("Lispro/ILPS apparent volume of distribution at 80 kg (L)") # Table 2 V lispro/ILPS 43.0 L (RSE 14.8%)
    e_wt_vc_lis <- 2.48
    label("Power exponent of body weight (WT/80) on lispro V/F (unitless)") # Table 2 BW on V lispro/ILPS 2.48 (RSE 19.8%)
    lcl_lis <- log(30.5)
    label("Lispro/ILPS apparent clearance (L/h)") # Table 2 CL lispro/ILPS 30.5 L/h (RSE 4.49%)

    # ---------------- Endogenous insulin baseline (regular/NPH observations)
    lrbase <- log(79.7)
    label("Baseline endogenous insulin concentration at 80 kg (pmol/L)") # Table 2 baseline endogenous insulin 79.7 pmol/L (RSE 4.78%)
    e_wt_rbase <- 0.739
    label("Power exponent of body weight (WT/80) on endogenous insulin baseline (unitless)") # Table 2 BW on baseline endogenous insulin 0.739 (RSE 35.3%)

    # ---------------- IIV: only the two magnitudes printed in the Results text
    etalvc_hum ~ 0.5268 # Results text: IIV 82.4% on V regular insulin; log(1 + 0.824^2)
    etalka_lis ~ 0.4448 # Results text: IIV 74.1% on Ka lispro; log(1 + 0.741^2)

    # ---------------- Residual error (proportional)
    propSd_glar <- 0.314
    label("Proportional residual error, glargine (fraction)") # Table 2 proportional glargine 31.4% (RSE 5.37%)
    propSd_hum <- 0.298
    label("Proportional residual error, human insulin (fraction)") # Table 2 proportional regular/NPH 29.8% (RSE 5.51%)
    propSd_lis <- 0.298
    label("Proportional residual error, lispro/ILPS (fraction)") # Table 2 proportional lispro/ILPS 29.8% (RSE 5.51%), same estimate as human insulin
  })
  model({
    # 1 U of insulin = 6 nmol = 6000 pmol; amounts are in U, volumes in L.
    u_to_pmol <- 6000

    # A dose covariate of 0 (insulin not given) falls back to its reference
    # dose so the power term stays finite; the depot it acts on is then empty.
    dose_glar <- DOSE_GLARGINE_U
    if (dose_glar <= 0) dose_glar <- 35
    dose_regu100 <- DOSE_REGULARU100_U
    if (dose_regu100 <= 0) dose_regu100 <- 10
    dose_hum <- DOSE_HUMANINSULIN_U
    if (dose_hum <= 0) dose_hum <- 10

    # Glargine (Supplement Methods 'Equations: Glargine insulin')
    ka_glar <- exp(lka_glar)
    vc_glar <- exp(lvc_glar) * (1 + e_wt_vc_glar * (WT - 71.4))
    cl_glar <- exp(lcl_glar)
    fdepot_glar <- exp(lfdepot_glar) * (dose_glar / 35)^e_dose_fdepot_glar
    ffo_glar <- expit(logitffo_glar)
    d1_glar <- exp(ld1_glar)

    # Regular U-100, U-500 and NPH ('Equations: Regular Insulins ... and NPH')
    ka_regu100 <- exp(lka_regu100) * (dose_regu100 / 10)^e_dose_ka_regu100
    ka_regu500 <- exp(lka_regu500)
    ka_nph <- exp(lka_nph)
    tlag_nph <- exp(ltlag_nph)
    vc_hum <- exp(lvc_hum + etalvc_hum) * (WT / 80)^e_wt_vc_hum
    cl_hum <- exp(lcl_hum) * (dose_hum / 10)^e_dose_cl_hum

    # Lispro and ILPS ('Equations: Insulin Lispro')
    ka_lis <- exp(lka_lis + etalka_lis)
    tlag_lis <- exp(ltlag_lis)
    ffo_lis <- expit(logitffo_lis)
    d1_lis <- exp(ld1_lis)
    ka_ilps <- exp(lka_ilps)
    vc_lis <- exp(lvc_lis) * (WT / 80)^e_wt_vc_lis
    cl_lis <- exp(lcl_lis)

    # Endogenous baseline ('Equation: Common Endogenous Baseline')
    rbase <- exp(lrbase) * (WT / 80)^e_wt_rbase

    d/dt(depot_glar) <- -ka_glar * depot_glar
    d/dt(central_glar) <- ka_glar * depot_glar - cl_glar / vc_glar * central_glar

    d/dt(depot_regu100) <- -ka_regu100 * depot_regu100
    d/dt(depot_regu500) <- -ka_regu500 * depot_regu500
    d/dt(depot_nph) <- -ka_nph * depot_nph
    d/dt(central_hum) <- ka_regu100 * depot_regu100 + ka_regu500 * depot_regu500 +
      ka_nph * depot_nph - cl_hum / vc_hum * central_hum

    d/dt(depot_lis) <- -ka_lis * depot_lis
    d/dt(depot_ilps) <- -ka_ilps * depot_ilps
    d/dt(central_lis) <- ka_lis * depot_lis + ka_ilps * depot_ilps -
      cl_lis / vc_lis * central_lis

    # Glargine: FRAC of the bioavailable dose first-order from depot_glar, the
    # rest zero-order into central_glar over d1 (dose central_glar with rate = -2).
    f(depot_glar) <- fdepot_glar * ffo_glar
    f(central_glar) <- fdepot_glar * (1 - ffo_glar)
    dur(central_glar) <- d1_glar

    alag(depot_nph) <- tlag_nph

    # Lispro: lag applies to the first-order part only; the zero-order part
    # starts at the injection.
    alag(depot_lis) <- tlag_lis
    f(depot_lis) <- ffo_lis
    f(central_lis) <- 1 - ffo_lis
    dur(central_lis) <- d1_lis

    Cc_glar <- u_to_pmol * central_glar / vc_glar
    Cc_hum <- u_to_pmol * central_hum / vc_hum + rbase
    Cc_lis <- u_to_pmol * central_lis / vc_lis

    Cc_glar ~ prop(propSd_glar)
    Cc_hum ~ prop(propSd_hum)
    Cc_lis ~ prop(propSd_lis)
  })
}
