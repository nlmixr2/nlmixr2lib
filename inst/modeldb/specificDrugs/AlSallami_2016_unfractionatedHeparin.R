AlSallami_2016_unfractionatedHeparin <- function() {
  description <- paste(
    "One-compartment population PK + linear pharmacodynamic model for",
    "unfractionated heparin (UFH) in paediatric patients receiving a single",
    "high intravenous bolus dose during cardiac angiography (Al-Sallami 2016).",
    "Fat-free mass (FFM) scales heparin clearance linearly and total body",
    "weight (WT) scales the central volume of distribution linearly. The PD",
    "layer is a linear concentration-effect model relating activated partial",
    "thromboplastin time (aPTT) to plasma heparin concentration (E0 + slope x",
    "Cc). The IV bolus was modelled in the source paper as a 0.1 h zero-order",
    "input (theta_D1 = 0.1 h, fixed); reproduce this in simulation by dosing",
    "the central compartment with rate = -2 to engage the model-defined",
    "duration. PD parameters were estimated sequentially via PPP&D with the",
    "PK parameters fixed at the values reported in Table 2."
  )
  reference <- paste(
    "Al-Sallami HS, Newall F, Monagle P, Ignjatovic V, Cranswick N, Duffull S.",
    "Development of a population pharmacokinetic-pharmacodynamic model of a",
    "single bolus dose of unfractionated heparin in paediatric patients.",
    "Br J Clin Pharmacol. 2016 Jul;82(1):178-184. doi:10.1111/bcp.12930.",
    sep = " "
  )
  vignette <- "AlSallami_2016_unfractionatedHeparin"
  units <- list(time = "h", dosing = "IU", concentration = "IU/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (exponent fixed at 1) scaling on central volume of",
        "distribution V with reference WT = 20 kg (Al-Sallami 2016 Table 2",
        "covariate model V = theta_V * WT / 20). 20 kg is the cohort median",
        "weight reported in the paper's covariate-model footnote."
      ),
      source_name        = "WT"
    ),
    FFM = list(
      description        = "Fat-free mass at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (exponent fixed at 1) scaling on heparin clearance CL with",
        "reference FFM = 15 kg (Al-Sallami 2016 Table 2 covariate model",
        "CL = theta_CL * FFM / 15). 15 kg is the cohort median FFM reported",
        "in the paper's covariate-model footnote. The paediatric FFM",
        "derivation uses the Al-Sallami et al. 2015 (Clin Pharmacokinet",
        "54:1169-1178) children-adapted variant of the Janmahasatian 2005",
        "semi-mechanistic adult formula."
      ),
      source_name        = "FFM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 1L,
    age_range      = "0.5-15.5 years (mean 6.7)",
    age_median     = "6.7 years (mean)",
    weight_range   = "6.7-68.6 kg (mean 23.6)",
    weight_median  = "23.6 kg (mean)",
    height_range   = "65-176 cm (mean 115.7)",
    bmi_range      = "11.5-24.7 kg/m^2 (mean 16.1)",
    sex_female_pct = 53.1,
    race_ethnicity = "Not reported (single-centre Royal Children's Hospital, Melbourne, Australia cohort)",
    disease_state  = paste(
      "Paediatric patients undergoing diagnostic cardiac angiography",
      "(catheterization). None of the children were being treated for",
      "thromboses at the time of the study. 8 subjects were < 2 years old."
    ),
    dose_range     = paste(
      "Single intravenous bolus of UFH at 75-100 IU/kg (mean 91 IU/kg,",
      "range 47.9-105.4 IU/kg). Total dose range 600-5000 IU (mean 2020 IU)."
    ),
    regions        = "Single-centre study at the Royal Children's Hospital, Melbourne, Victoria, Australia",
    n_observations = "231 UFH concentration measurements; 290 aPTT measurements (43% above the upper limit of quantification of 999 s)",
    sampling       = "Blood samples drawn at baseline and at 15, 30, 45, and 120 min post-dose",
    notes          = paste(
      "Baseline demographics from Al-Sallami 2016 Table 1. UFH concentration",
      "was quantified by a modified protamine titration assay; aPTT was",
      "measured with the PTT-A kit (Diagnostica Stago) on the STA-R analyser",
      "with the upper limit of clot detection modified to measure up to 999 s",
      "to accommodate the high single dose. Censored aPTT data above the",
      "upper limit were accounted for using Beal's M3 method modified for",
      "right censoring. Two patients had missing height values that were",
      "imputed via multivariate linear regression on age, sex, and weight",
      "before deriving FFM. Estimation: NONMEM v7.2 with Wings for NONMEM",
      "v720, first-order conditional estimation with interaction, combined",
      "residual error model. Bootstrap 95% CIs from 1000 resamples. VPC",
      "validation from 1000 simulated datasets (paper Figures 1 and 2)."
    )
  )

  ini({
    # --- Pharmacokinetics (Al-Sallami 2016 Table 2 final model) ---
    # The IV bolus is modelled as a short fixed-duration zero-order input:
    # theta_D1 = 0.1 h was fixed during estimation (Table 2 footnote;
    # "the infusion rate parameter D1 was fixed"). CL and V scale linearly
    # (exponent 1) on FFM and WT, respectively, with reference FFM = 15 kg
    # and reference WT = 20 kg (Al-Sallami 2016 Table 2 covariate model).
    lcl  <- log(0.603)
    label("Clearance at reference FFM = 15 kg (CL, L/h)")                              # Al-Sallami 2016 Table 2 final model: theta_CL = 0.603 (L/h per 15 kg FFM)
    lvc  <- log(0.751)
    label("Central volume of distribution at reference WT = 20 kg (V, L)")             # Al-Sallami 2016 Table 2 final model: theta_V = 0.751 (L per 20 kg WT)
    tdur <- fixed(0.1)
    label("IV bolus zero-order input duration (D1, h); fixed")                          # Al-Sallami 2016 Table 2 footnote: theta_D1 = 0.1 h (fixed)

    # --- PK inter-individual variability (Al-Sallami 2016 Table 2) ---
    # Exponential ETA model. omega^2 = log(CV^2 + 1).
    #   CV(CL) = 50%  -> omega^2 = log(1 + 0.50^2) = 0.223144
    #   CV(V)  = 40%  -> omega^2 = log(1 + 0.40^2) = 0.148420
    # Block covariance with corr(CL, V) = 0.75:
    #   cov(CL, V) = 0.75 * sqrt(0.223144 * 0.148420) = 0.136490
    etalcl + etalvc ~ c(
      log(1 + 0.5^2),
      0.75 * sqrt(log(1 + 0.5^2) * log(1 + 0.4^2)), log(1 + 0.4^2)
    )  # Al-Sallami 2016 Table 2 final model: omega_CL = 50% CV, omega_V = 40% CV, Corr(CL, V) = 0.75

    # --- PK residual error (Al-Sallami 2016 Table 2 final model) ---
    propSd <- 0.17
    label("Heparin proportional residual error (fraction)")                            # Al-Sallami 2016 Table 2 final model: sigma_prop = 17%
    addSd  <- 90
    label("Heparin additive residual error (IU/L)")                                    # Al-Sallami 2016 Table 2 final model: sigma_add = 90 U/L

    # --- Pharmacodynamics (Al-Sallami 2016 Table 3 final model) ---
    # Linear concentration-effect: aPTT(t) = E0 + slope * Cc(t). PD parameters
    # were estimated sequentially (PPP&D method, Zhang et al. 2003) with the
    # PK parameters fixed at the Table 2 final-model estimates above. The
    # parameters marked with * in Table 3 are fixed from the PK model.
    lrbase <- log(35.6)
    label("Baseline aPTT (E0, s)")                                                     # Al-Sallami 2016 Table 3 final model: theta_E0 = 35.6 s
    lslope <- log(0.67)
    label("Linear slope of aPTT vs heparin concentration (slope, s per IU/L)")         # Al-Sallami 2016 Table 3 final model: theta_SLP = 0.67

    # --- PD inter-individual variability (Al-Sallami 2016 Table 3) ---
    # Table 3's omega column is labelled "(%)" with the legend "presented as
    # coefficient of variation percentage, CV%". The values are reported
    # inconsistently within the same column: omega_SLP = 64 is clearly 64%
    # CV (integer-percent form, matched by bootstrap 63 (50, 76)), while
    # omega_E0 = 0.43 is a decimal whose literal CV% reading (0.43% CV)
    # implies essentially no between-subject variability in baseline aPTT,
    # which is implausible for a paediatric cohort. The bootstrap CI for
    # omega_E0 is also reported on the decimal scale (0.44 (0.38, 0.51)),
    # consistent with the value being on a decimal-CV scale. We interpret
    # omega_E0 = 0.43 as decimal-form CV (= 43% CV) here; the alternative
    # interpretations are recorded in the vignette Errata section.
    #   CV(E0)  = 43% (interpreted) -> omega^2 = log(1 + 0.43^2) = 0.169658
    #   CV(SLP) = 64%               -> omega^2 = log(1 + 0.64^2) = 0.343306
    etalrbase ~ log(1 + 0.43^2)  # Al-Sallami 2016 Table 3 final model: omega_E0  = 0.43 (interpreted as 43% CV; see vignette Errata)
    etalslope ~ log(1 + 0.64^2)  # Al-Sallami 2016 Table 3 final model: omega_SLP = 64% CV

    # --- PD residual error (Al-Sallami 2016 Table 3 final model) ---
    # The Table 3 sigma_add unit label "(U l-1)" appears to be a copy-paste
    # carry-over from Table 2 (the PK additive RUV is in IU/L). The PD
    # response is aPTT (units of seconds), so the unit label is wrong for
    # this row; the numerical value 0.005 is interpreted here as the
    # additive SD on aPTT in seconds. It is essentially negligible relative
    # to the 30% proportional component and only contributes near a near-
    # zero predicted aPTT (an extrapolation regime not observed in the
    # study). The unit-label discrepancy is recorded in vignette Errata.
    propSd_aPTT <- 0.30
    label("aPTT proportional residual error (fraction)")                               # Al-Sallami 2016 Table 3 final model: sigma_prop = 30%
    addSd_aPTT  <- 0.005
    label("aPTT additive residual error (s; unit-label discrepancy: see Errata)")     # Al-Sallami 2016 Table 3 final model: sigma_add = 0.005 (table unit reads U/L; interpreted as seconds for the aPTT response)
  })

  model({
    # --- Individual PK parameters (FFM linear on CL, WT linear on V) ---
    cl <- exp(lcl + etalcl) * (FFM / 15)
    vc <- exp(lvc + etalvc) * (WT  / 20)
    kel <- cl / vc

    # --- ODE: single central compartment receiving zero-order IV input ---
    d/dt(central) <- -kel * central

    # --- Model-defined zero-order input duration ---
    # Dose records must specify rate = -2 to engage the model duration.
    # Doses in IU; vc in L; so central / vc has units of IU/L (= U/L)
    # matching the paper's reported residual error units.
    dur(central) <- tdur

    # --- Outputs (Cc concentration; linear PD on aPTT) ---
    Cc    <- central / vc
    rbase <- exp(lrbase + etalrbase)
    slope <- exp(lslope + etalslope)
    aPTT  <- rbase + slope * Cc

    # Multi-output residual error: simulation event rows specify which
    # endpoint is observed via cmt = "Cc" or cmt = "aPTT".
    Cc   ~ add(addSd)      + prop(propSd)
    aPTT ~ add(addSd_aPTT) + prop(propSd_aPTT)
  })
}
