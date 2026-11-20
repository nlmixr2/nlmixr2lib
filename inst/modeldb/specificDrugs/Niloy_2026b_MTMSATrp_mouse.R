Niloy_2026b_MTMSATrp_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Tumor growth inhibition and time-to-event",
    "(TGI-TTE) model for MTMSA-Trp, a mithramycin analogue, in female",
    "athymic nude mice bearing subcutaneous TC-32 Ewing sarcoma",
    "xenografts (Niloy 2026, Cancer Chemother Pharmacol). Simeoni",
    "unperturbed growth (exponential rate lambda0 switching to linear",
    "rate lambda1, psi = 20) with a direct linear kill term k * Cp, driven",
    "by the typical-value one-compartment mouse PK model of the same",
    "group (dose-power effect on clearance; Niloy_2026_MTMSATrp_mouse).",
    "The second stage is a log-logistic time-to-event model (shape 4) for",
    "the model-derived time to a 1500 mm^3 tumor volume, whose scale",
    "depends on the individual lambda0 and lambda1 and on the regimen",
    "average plasma concentration (CAV)."
  )
  reference <- paste(
    "Niloy KK, Horn J, Bhuiyan NH, Bhosale SS, Shaaban KA, Prisinzano TE,",
    "Thorson JS, Rohr J, Leggas M. Bridging tumor growth dynamics and",
    "survival in preclinical oncology drug development - a Tumor Growth",
    "Inhibition and Time-to-Event (TGI-TTE) modeling framework. Cancer",
    "Chemother Pharmacol. 2026;96(1):100. doi:10.1007/s00280-026-04947-0.",
    "PK layer: Niloy KK et al. Nonlinear Mixed-Effects Modeling to",
    "Characterize the Pharmacokinetics of a Novel Mithramycin Analogue for",
    "Ewing Sarcoma in Mice. Res Sq. 2026. doi:10.21203/rs.3.rs-9035594/v1."
  )
  vignette <- "Niloy_2026b_MTMSATrp_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL", tumor = "mm^3")

  compartmentData <- list(
    central = list(analyte = "MTMSATrp", units = "mg/kg", specimen = "plasma", verified = TRUE),
    tumor_size = list(analyte = "TC-32 xenograft tumor volume", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list(
    TUM_VOL = list(
      description = "Observed tumor volume of the mouse at time zero (start of treatment), used as the initial condition of the tumor state",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Niloy 2026 Results, TGI modeling: TV0 'was included in the model as a",
        "regressor rather than estimated as a population parameter because",
        "initial tumor volume measurements were available for each mouse'.",
        "Initial volumes ranged from 37 to 251 mm^3 (Methods, Animal experiments)."
      ),
      source_name = "TV0"
    ),
    COHDOSE = list(
      description = "Per-dose amount of the MTMSA-Trp regimen (mg/kg); selects the dose-specific typical clearance of the PK layer",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Niloy 2026 Methods, TGI modeling: 'the typical population clearance",
        "was calculated for each dose level used in the xenograft studies",
        "based on this dose-clearance relationship'. Power effect",
        "(COHDOSE / 3)^-0.30 from the upstream PK model. Vehicle animals",
        "receive no dose; a COHDOSE of 0 is replaced by the 3 mg/kg",
        "reference inside model() so that clearance stays finite."
      ),
      source_name = "Dose"
    ),
    CAV = list(
      description = "Regimen-specific average MTMSA-Trp plasma concentration over the dosing interval (AUC_tau / tau), 0 for vehicle",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Niloy 2026 Methods, TTE modeling: exposure 'was summarized as",
        "regimen-specific average concentration' computed by NCA of the",
        "concentration-time profile simulated from the PK model for each",
        "regimen. Regimen-level (not individual) exposure. Enters the TTE",
        "scale as exp(beta_Cavg * CAV / median CAV) (Eq. 3); the median",
        "(393.7 ng/mL) is not printed and is reconstructed from the PK",
        "model and the Table S1 group sizes (see model() and the vignette)."
      ),
      source_name = "Cavg"
    )
  )

  population <- list(
    species = "mouse (female athymic nude, NU(NCr)-Foxn1nu; TC-32 Ewing sarcoma xenograft)",
    n_subjects = 104L,
    n_observations = 1377L,
    n_studies = 3L,
    age_range = "(not reported)",
    weight_range = "(not reported; doses in mg/kg)",
    sex_female_pct = 100,
    disease_state = "Subcutaneous TC-32 Ewing sarcoma xenograft (rear left flank)",
    dose_range = paste(
      "Vehicle; MTMSA-Trp IV bolus 0.3, 0.6 or 0.9 mg/kg once daily x 5 (experiment 1),",
      "2.40 mg/kg every 3 days x 6 (experiments 2 and 3), 2.85 mg/kg twice weekly x 6 (experiment 3)"
    ),
    regions = "United States (St. Jude Children's Research Hospital, Memphis TN)",
    tumor_volume_range = "37-251 mm^3 at the start of treatment",
    notes = paste(
      "TGI dataset: 1377 tumor volumes from 104 mice (three experiments",
      "pooled with historical TC-32 vehicle controls), followed up to 100",
      "days. TTE dataset (Table S1): 62 mice of experiments 1-3, 54 events",
      "and 8 right-censored. The event is the TGI-model-predicted time for",
      "the tumor to reach 1500 mm^3, so the TTE model is a second-stage",
      "description of a model-derived endpoint, not an independent",
      "survival analysis (Niloy 2026 Discussion)."
    )
  )

  ini({
    # ------------------------------------------------------------
    # PK layer: typical-value one-compartment model of the same group
    # (Niloy 2026 reference [19]; nlmixr2lib Niloy_2026_MTMSATrp_mouse,
    # Table 1 of that preprint). The TGI paper simulated the typical
    # profile of each regimen and did not re-estimate the PK.
    # ------------------------------------------------------------
    lcl <- fixed(log(39.18)); label("Apparent clearance at the 3 mg/kg reference dose, taken from the Niloy 2026 PK model (mL/h/kg)") # Niloy 2026 Res Sq Table 1: CL = 39.18
    lvc <- fixed(log(53.06)); label("Volume of distribution, taken from the Niloy 2026 PK model (mL/kg)") # Niloy 2026 Res Sq Table 1: V = 53.06
    e_cohdose_cl <- fixed(-0.30); label("Dose (mg/kg) power exponent on clearance, taken from the Niloy 2026 PK model (unitless)") # Niloy 2026 Res Sq Table 1: beta = -0.30

    # ------------------------------------------------------------
    # Simeoni TGI model -- Table 1 (estimates in 1/h, mm^3/h, mL/ng/h)
    # ------------------------------------------------------------
    ltumorExpGrowth <- log(0.0063); label("Exponential tumor growth rate lambda0 (1/h)") # Table 1: lambda0 = 0.0063 (RSE 8.56)
    ltumorLinGrowth <- log(5.93); label("Linear tumor growth rate lambda1 (mm^3/h)") # Table 1: lambda1 = 5.93 (RSE 11.8)
    ldrugSlope <- log(0.000024); label("Linear tumor kill rate constant k (mL/ng/h)") # Table 1: k = 0.000024 (RSE 13.4)
    psi <- fixed(20); label("Exponent of the exponential-to-linear growth switch psi (unitless)") # Eq. 1 and Fig. 1: exponent 20 and 1/20 are constants

    # ------------------------------------------------------------
    # Log-logistic TTE model -- Table 2 and Eq. 2-3
    # ------------------------------------------------------------
    llambda <- log(703.17); label("Typical log-logistic scale Te,pop (median event time at median covariates and zero exposure) (h)") # Table 2: Te = 703.17 (RSE 4.72)
    lgamma <- fixed(log(4)); label("Log-logistic shape s (unitless)") # Table 2: s = 4, not estimated; Table S6 note: s not identifiable
    e_tumorExpGrowth_lambda <- -0.75; label("Coefficient of median-normalized lambda0 on the TTE scale (unitless)") # Table 2: beta_lambda0 = -0.75 (RSE 8.3)
    e_tumorLinGrowth_lambda <- -0.46; label("Coefficient of median-normalized lambda1 on the TTE scale (unitless)") # Table 2: beta_lambda1 = -0.46 (RSE 17.9)
    e_cav_lambda <- 0.63; label("Coefficient of median-normalized regimen average concentration on the TTE scale (unitless)") # Table 2: beta_Cavg = 0.63 (RSE 10.4)

    # ------------------------------------------------------------
    # IIV. Table 1 and Table 2 report omega (SD of the log-normal
    # random effect) followed by the %CV, e.g. 0.85 (102.82):
    # sqrt(exp(0.85^2) - 1) = 1.029. Variance = omega^2.
    # ------------------------------------------------------------
    etaltumorExpGrowth ~ 0.7225 # Table 1: omega lambda0 = 0.85 (CV 102.82); 0.85^2
    etaltumorLinGrowth ~ 0.5329 # Table 1: omega lambda1 = 0.73 (CV 83.39); 0.73^2
    etaldrugSlope ~ 0.7569 # Table 1: omega k = 0.87 (CV 106.68); 0.87^2
    etallambda ~ 0.001849 # Table 2: omega Te = 0.043 (CV 4.26); 0.043^2

    # ------------------------------------------------------------
    # Residual error on tumor volume -- Table 1 (Monolix a, b)
    # ------------------------------------------------------------
    addSd_tumor_size <- 11.22; label("Additive residual error on tumor volume (mm^3)") # Table 1: a = 11.22
    propSd_tumor_size <- 0.22; label("Proportional residual error on tumor volume (fraction)") # Table 1: b = 0.22
  })

  model({
    # ------------------------------------------------------------
    # 1. PK layer (typical value; no IIV). Dose-power effect on CL
    #    with reference dose 3 mg/kg. Vehicle (COHDOSE = 0) has no
    #    dose, so the reference dose is used to keep CL finite.
    # ------------------------------------------------------------
    cohdose_cl <- COHDOSE
    if (COHDOSE <= 0) cohdose_cl <- 3
    cl <- exp(lcl) * (cohdose_cl / 3)^e_cohdose_cl
    vc <- exp(lvc)
    kel <- cl / vc

    d/dt(central) <- -kel * central
    # mg/kg / (mL/kg) = mg/mL; x 1e6 -> ng/mL
    Cc <- central / vc * 1e6

    # ------------------------------------------------------------
    # 2. Simeoni TGI (Eq. 1) with a direct linear kill term
    # ------------------------------------------------------------
    tumorExpGrowth <- exp(ltumorExpGrowth + etaltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth + etaltumorLinGrowth)
    drugSlope <- exp(ldrugSlope + etaldrugSlope)

    growth <- tumorExpGrowth * tumor_size / (1 + (tumorExpGrowth / tumorLinGrowth * tumor_size)^psi)^(1 / psi)
    d/dt(tumor_size) <- growth - drugSlope * Cc * tumor_size
    tumor_size(0) <- TUM_VOL

    # ------------------------------------------------------------
    # 3. Log-logistic TTE (Eq. 2) with covariates on the scale (Eq. 3):
    #    Te_i = Te_pop * exp(b_l0 * (lambda0_i / med_l0 - 1))
    #                  * exp(b_l1 * (lambda1_i / med_l1 - 1))
    #                  * exp(b_cavg * Cavg / med_Cavg) * exp(eta_i)
    #    The Cavg term has no "- 1" as printed in Eq. 3, so Te_pop is
    #    the scale at zero exposure (vehicle).
    #    The medians are not printed. med_l0 and med_l1 are taken as the
    #    Table 1 typical values; med_Cavg = 393.7 ng/mL is the Cavg of
    #    the 0.6 mg/kg once-daily regimen, which holds the 31st and 32nd
    #    of the 62 Table S1 mice when ranked by Cavg:
    #    0.6 / (39.18 * (0.6 / 3)^-0.30) / 24 * 1e6 = 393.7 ng/mL.
    # ------------------------------------------------------------
    med_tumorExpGrowth <- 0.0063
    med_tumorLinGrowth <- 5.93
    med_cav <- 393.7

    lambda <- exp(llambda + etallambda) *
      exp(e_tumorExpGrowth_lambda * (tumorExpGrowth / med_tumorExpGrowth - 1)) *
      exp(e_tumorLinGrowth_lambda * (tumorLinGrowth / med_tumorLinGrowth - 1)) *
      exp(e_cav_lambda * CAV / med_cav)
    gamma <- exp(lgamma)

    # Closed-form log-logistic survival; shape 4 > 1 so the hazard is 0 at t = 0
    tscaled <- t / lambda
    hazard <- (gamma / lambda) * tscaled^(gamma - 1) / (1 + tscaled^gamma)
    cumhaz <- log(1 + tscaled^gamma)
    sur <- 1 / (1 + tscaled^gamma)

    tumor_size ~ add(addSd_tumor_size) + prop(propSd_tumor_size) + combined1()
  })
}
