Nanga_2019_tacrolimus_metaanalysis <- function() {
  description <- "MBMA. Two-compartment population PK meta-model for oral tacrolimus in solid organ transplantation (Nanga 2019), built from pooled individual-patient data across 7 historical NONMEM datasets (n = 281 paediatric + adult liver and kidney transplant recipients). Structural model: first-order absorption with fixed lag time, time-varying first-order elimination, allometric (WT/50 kg) scaling on apparent clearance and apparent central volume, multiplicative reduction of CL/F in hepatic-graft recipients, sigmoidal post-transplant-day recovery of CL/F, and reduced relative bioavailability for the oral syrup formulation. The literature-review summary table (Nanga 2019 Table 2: 76 published popPK models) is not used for parameter fitting and is not reproduced here."
  reference <- "Nanga TM, Doan TTP, Marquet P, Musuamba FT. Toward a robust tool for pharmacokinetic-based personalization of treatment with tacrolimus in solid organ transplantation: A model-based meta-analysis approach. Br J Clin Pharmacol. 2019;85(12):2793-2823. doi:10.1111/bcp.14110"
  vignette <- "Nanga_2019_tacrolimus_metaanalysis"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Nanga 2019. Allometric power scaling on CL/F (exponent 0.61, Table 3 WT_CL) and on the apparent volumes of distribution V2/F and V3/F (exponent 0.53, Table 3 WT_V) with reference WT = 50 kg (the pooled-cohort median per Table 5; range 5 - 128 kg covering paediatric to adult patients). The single WT_V coefficient is applied to both Vc and Vp because the paper's text describes the body-weight effect on 'the apparent volume of distribution' (V/F = Vc + Vp).",
      source_name        = "WT"
    ),
    POD = list(
      description        = "Days post-transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. CL/F follows a sigmoidal recovery from the immediate post-transplant low to twice that level with increasing POD: time_factor_CL = 1 + POD^8.88 / (POD^8.88 + 6.12^8.88), with Hill steepness 8.88 (Table 3 'Sigmoidity coefficient for time_CL') and half-recovery at POD = 6.12 days (Table 3 'Time 50% recovery'). The factor is 1 at POD = 0 and asymptotes to 2 as POD -> infinity.",
      source_name        = "POD"
    ),
    TX_LIVER = list(
      description        = "Liver (hepatic) graft indicator: 1 if the patient received a liver transplant, 0 if a non-liver solid-organ transplant (kidney, heart, or lung) -- the pooled model-building cohort includes only liver and kidney recipients (Nanga 2019 Table 5).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-liver solid-organ graft)",
      notes              = "Time-fixed at transplantation date per subject. Nanga 2019 model-building cohort (Table 5): 201 liver-transplant patients (71.5%) and 80 non-liver (kidney) patients (28.5%). Encoded as a multiplicative power coefficient theta^TX_LIVER on CL/F per the paper's general categorical-covariate form (Eq. 4): liver-transplant patients have CL/F reduced by factor 0.38 relative to non-liver recipients (Table 3 'Hepatic trans_CL' = 0.38; printed in the explicit CL equation on p.2815 as 0.39, a rounding/typesetting artifact -- the Table 3 value 0.38 is authoritative).",
      source_name        = "Hepatic trans_CL"
    ),
    FORM_SYRUP = list(
      description        = "Oral syrup-formulation indicator: 1 if the patient received the oral suspension / syrup formulation of tacrolimus on a given dose, 0 if the patient received the capsule (immediate-release) formulation. Per-dose-occasion in principle; in paediatric cohorts typically time-fixed per subject by clinical convention.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule / standard solid oral formulation)",
      notes              = "Time-fixed per-dose-occasion. Nanga 2019 Table 3 estimates the syrup-formulation relative bioavailability F_syrup = 0.53 (95% bootstrap CI 0.31 - 0.75). Encoded inside model() as f(depot) <- F_syrup^FORM_SYRUP per the paper's categorical-covariate form (Eq. 4), so capsule users have F = 1 and syrup users have F = 0.53. The paediatric subset of databases 1 - 3 (Table 1) contributed the syrup-formulation observations.",
      source_name        = "syrup formulation"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 281L,
    n_studies            = 7L,
    age_range            = "0.3 - 68 years",
    age_median           = "2.3 years",
    weight_range         = "5 - 128 kg",
    weight_median        = "50 kg",
    sex_female_pct       = 42.2,
    sex_distribution     = "104 male / 76 female of subjects with non-missing sex; 38% of records missing sex.",
    race_ethnicity       = "Not reported in the source paper.",
    disease_state        = "Paediatric and adult recipients of liver (201 patients, 71.5%) or kidney (80 patients, 28.5%) solid-organ transplants on oral tacrolimus immunosuppression. Observations span D1 - D394 post-transplant.",
    dose_range           = "Oral tacrolimus dosing varied by centre and age (capsule or syrup formulation); per-database dosing schedules are summarised in Nanga 2019 Table 1.",
    organ_breakdown      = "Model-building cohort (Table 5): 201 liver-transplant patients (71.5%); 80 kidney-transplant patients (28.5%). External-validation cohorts (113 additional patients, NOT used for parameter estimation): 61 adult lung (Stimmugrep trial), 20 adult heart (Pigrec trial), 32 adult kidney (PCCP trial).",
    formulations         = "Mixed capsule and oral-suspension (syrup) formulations across paediatric and adult datasets. Relative bioavailability for syrup vs capsule estimated at 0.53 (Table 3).",
    regions              = "Multi-centre European and international cohorts (databases 1 - 7 model-building from in-house and previously published data; databases 8 - 10 external evaluation).",
    notes                = "Pooled patient-level dataset from 7 historical NONMEM-format individual-patient-data sources (Nanga 2019 Table 1, databases 1 - 7). The paper also reviewed 76 published tacrolimus popPK models (Table 2) and used those for narrative comparison only -- they were NOT used for parameter fitting. Median body weight 50 kg, age 2.3 years (0.3 - 68). Haematocrit, plasma albumin, and CYP3A5 status had > 50% missing data in the pooled dataset, so CYP3A5 polymorphism was not included as a covariate."
  )

  ini({
    # Structural PK -- Nanga 2019 Table 3 final-model estimates. Time in hours;
    # apparent clearances (CL/F, Q/F) in L/h; apparent volumes (V2/F = central,
    # V3/F = peripheral) in L. The Table 3 reference subject is WT = 50 kg
    # (pooled-cohort median), TX_LIVER = 0 (non-liver graft), FORM_SYRUP = 0
    # (capsule), and POD = 0 (day of transplantation; the time-varying factor
    # below ramps CL/F up by a factor of 2 to its stable-period asymptote).
    lcl   <- log(22.5)        ; label("Apparent oral clearance CL/F at WT = 50 kg, non-liver graft, capsule, POD = 0 (L/h)")  # Nanga 2019 Table 3 final CL/F = 22.5 L/h
    lvc   <- log(246.2)       ; label("Apparent central volume V2/F at WT = 50 kg (L)")                                       # Nanga 2019 Table 3 final V2/F = 246.2 L
    lq    <- log(24.2)        ; label("Apparent inter-compartmental clearance Q/F (L/h)")                                     # Nanga 2019 Table 3 final Q/F = 24.2 L/h
    lvp   <- log(109.9)       ; label("Apparent peripheral volume V3/F at WT = 50 kg (L)")                                    # Nanga 2019 Table 3 final V3/F = 109.9 L
    lka   <- fixed(log(3.37)) ; label("Absorption rate constant ka (1/h; fixed)")                                              # Nanga 2019 Table 3 final KA = 3.37 1/h (fixed per Results section text, p.2815)
    ltlag <- fixed(log(0.32)) ; label("Absorption lag time ALAG1 (h; fixed)")                                                  # Nanga 2019 Table 3 final ALAG1 = 0.32 h (fixed per Results section text, p.2815)

    # Allometric exponents per Nanga 2019 Eq. (3): P_i = TVp * (WT_i / WT_med)^theta_WT.
    e_wt_cl <- 0.61  ; label("Allometric exponent of (WT / 50 kg) on CL/F (unitless)")            # Nanga 2019 Table 3 WT_CL = 0.61
    e_wt_vc_vp <- 0.53  ; label("Shared allometric exponent of (WT / 50 kg) on V2/F and V3/F (unitless)")  # Nanga 2019 Table 3 WT_V  = 0.53 (single shared exponent applied to both volumes)

    # Categorical covariate effects per Nanga 2019 Eq. (4): P_i = TVp * theta^COV.
    e_tx_liver_cl  <- 0.38   ; label("Hepatic-transplant multiplicative factor on CL/F (theta_HEPA_CL)")     # Nanga 2019 Table 3 Hepatic_trans_CL = 0.38 (printed as 0.39 in the explicit CL equation on p.2815; the Table 3 value is authoritative)
    e_form_syrup_f <- 0.53   ; label("Relative bioavailability of oral syrup vs. capsule (theta_F_syrup)") # Nanga 2019 Table 3 Bioavailability for syrup formulation = 0.53

    # Sigmoidal post-transplant-day recovery on CL/F (Nanga 2019 explicit CL
    # equation, p.2815): time_factor_CL = 1 + POD^sigm / (POD^sigm + t50^sigm).
    # Ranges from 1 at POD = 0 (day of surgery) to 2 as POD -> infinity, with
    # 50% recovery at POD = 6.12 days and Hill steepness 8.88.
    sigm_cl <- 8.88  ; label("Sigmoid Hill coefficient for the POD effect on CL/F (unitless)")   # Nanga 2019 Table 3 Sigmoidity coefficient for time_CL = 8.88
    t50_cl  <- 6.12  ; label("Half-recovery POD for the time-varying effect on CL/F (days)")   # Nanga 2019 Table 3 Time 50% recovery = 6.12 days

    # Inter-individual variability -- Nanga 2019 Table 3 reports IIV as %CV;
    # the paper's exponential / log-normal IIV model (Eq. 1) maps CV% to
    # log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F  CV  59.4% -> log(0.594^2 + 1) = 0.302
    #   V2/F  CV 133.2% -> log(1.332^2 + 1) = 1.020
    # IIV on Q/F, V3/F, ka, and ALAG1 was not estimated in the final model.
    etalcl ~ 0.302   # Nanga 2019 Table 3 IIV CL/F  = 59.4% CV
    etalvc ~ 1.020   # Nanga 2019 Table 3 IIV V2/F = 133.2% CV

    # Residual unexplained variability -- combined proportional + additive
    # (Nanga 2019 Eq. 2: Y = F * (1 + eps1) + eps2). Concentrations in ng/mL.
    propSd <- 0.2432 ; label("Proportional residual error (fraction)")  # Nanga 2019 Table 3 eps1 = 24.32%
    addSd  <- 3.22   ; label("Additive residual error (ng/mL)")           # Nanga 2019 Table 3 eps2 = 3.22 ng/mL
  })

  model({
    # Sigmoidal post-transplant-day factor on CL/F per the published explicit
    # CL equation (Nanga 2019 p.2815): time_factor_CL = 1 + POD^sigm /
    # (POD^sigm + t50^sigm). Ranges 1 -> 2 with 50% recovery at POD = 6.12 d.
    time_factor_cl <- 1 + POD ^ sigm_cl / (POD ^ sigm_cl + t50_cl ^ sigm_cl)

    # Individual PK parameters at reference WT = 50 kg, TX_LIVER = 0,
    # FORM_SYRUP = 0. Categorical theta^COV per Eq. (4): liver -> 0.38, non-liver -> 1.
    cl <- exp(lcl + etalcl) *
          (WT / 50) ^ e_wt_cl *
          e_tx_liver_cl ^ TX_LIVER *
          time_factor_cl
    vc <- exp(lvc + etalvc) * (WT / 50) ^ e_wt_vc_vp
    vp <- exp(lvp)          * (WT / 50) ^ e_wt_vc_vp
    q  <- exp(lq)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption + lag time. Relative
    # bioavailability of the syrup formulation is applied multiplicatively to
    # the depot per Eq. (4): capsule -> 1, syrup -> 0.53.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- e_form_syrup_f ^ FORM_SYRUP
    alag(depot) <- tlag

    # Tacrolimus whole-blood concentrations in ng/mL. Convert internal mg / L
    # (dose mg, vc L) to ng/mL by * 1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
