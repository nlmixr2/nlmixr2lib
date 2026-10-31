Chang_2024_mavacamten <- function() {
  description <- "Two-compartment population PK model for oral mavacamten in healthy adults and adults with obstructive or nonobstructive hypertrophic cardiomyopathy, with first-order absorption after a lag time, CYP2C19-phenotype, sex, eGFR, health-status and omeprazole/esomeprazole effects on clearance, and dose-, formulation-, food- and sex-dependent absorption (Chang 2024)"
  reference <- "Chang P, Perera V, Salinger DH, Merali S, Thanneer N, Back H, Seroogy JD, Gretler DD, Sehnert AJ, Palmisano M, Roy A. Characterization of mavacamten pharmacokinetics in patients with hypertrophic cardiomyopathy to inform dose titration. CPT Pharmacometrics Syst Pharmacol. 2024;13(9):1462-1475. doi:10.1002/psp4.13197"
  vignette <- "Chang_2024_mavacamten"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects with reference 84 kg (the reference participant of Chang 2024 Methods 'PopPK model application' and Figure 1): exponent fixed to 0.75 on CL/F and Q/F, and one estimated exponent of 0.813 shared by V2/F and V3/F (Table S5 'Body weight-Vx/F').",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F with exponent 0.205 and reference 95 mL/min/1.73 m^2 (reference participant, Figure 1). The paper does not name the estimating equation. Cohort range 29.5-148 (Results). Baseline value.",
      source_name = "eGFR"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Proportional effects (1 + theta)^SEXF on CL/F (+27.3%) and on F (+9.45%). The net change in apparent clearance is 1.273 / 1.0945 = 1.16, the 'approximately 16% higher' CL/F the paper reports for women.",
      source_name = "Female sex"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with obstructive or nonobstructive hypertrophic cardiomyopathy)",
      notes = "Proportional effect (1 + 0.272)^DIS_HEALTHY on CL/F. Also gates the food-status effects on KA and F, which were estimated only in healthy participants (food status was unknown in all patient records): see FED.",
      source_name = "Healthy individual / participant type"
    ),
    FED = list(
      description = "Fed state at dosing (healthy participants only)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted). Read only when DIS_HEALTHY = 1; patient records (food status unknown) take the reference absorption whatever FED holds.",
      notes = "Healthy participants dosed fasted (DIS_HEALTHY * (1 - FED) = 1) have KA x 1.873 and F x 0.686; healthy participants dosed fed (DIS_HEALTHY * FED = 1) have F x 0.567 and the reference KA (the fed-healthy KA effect was dropped from the final model, Results 'Final PopPK model'). Per dose record.",
      source_name = "Food status"
    ),
    FORM_SOLUTION = list(
      description = "Oral-solution formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet or capsule)",
      notes = "Proportional effects on KA (1 + 0.681)^FORM_SOLUTION and on F (1 - 0.194)^FORM_SOLUTION, relative to the pooled tablet/capsule reference. The solution was given only in the single-ascending-dose study MYK-461-002 (Table 2: 36 participants).",
      source_name = "Formulation (solution)"
    ),
    DOSE_MAVACAMTEN_MG = list(
      description = "Mavacamten amount of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects with reference 5 mg (the reference participant's 5 mg once-daily dose, Figure 1 'Dose - 5 (1-25) mg'): exponent -0.135 on KA, -0.087 on ALAG1, 0.0777 on F. Time-varying: the paper carries the value of each dose record forward to the following observation records, so with dose titration set it on every dose record to that record's amount. Place the column after the event columns so rxode2 does not read it as the dose amount.",
      source_name = "Dose"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with all CYP2C19 indicators 0 the reference is the normal metabolizer (*1/*1)",
      notes = "Proportional effect (1 - 0.722)^CYP2C19_PM on CL/F (72% lower). Genotypes *2/*2, *2/*3, *2A/*2A, *2A/*2B, *2A/*9 (Table S2).",
      source_name = "CYP2C19:PM"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with all CYP2C19 indicators 0 the reference is the normal metabolizer (*1/*1)",
      notes = "Proportional effect (1 - 0.331)^CYP2C19_IM on CL/F (33% lower). One loss-of-function allele with a *1 allele (*1/*2, *1/*3, *1/*4, *1/*6, *1/*8; Table S2, where this group is printed as 'EM'; Figure 1 n = 89). The Results report no difference in CL/F between *2/*17 and *1/*2 carriers, so *2/*17 diplotypes are also coded IM here.",
      source_name = "CYP2C19:IM"
    ),
    CYP2C19_RM = list(
      description = "CYP2C19 rapid-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with all CYP2C19 indicators 0 the reference is the normal metabolizer (*1/*1)",
      notes = "Proportional effect (1 - 0.00218)^CYP2C19_RM on CL/F (-0.2%, RSE 2047%); kept in the final model by the authors although the SCM did not select it (Results 'Final PopPK model'). Genotype *1/*17 only; *17/*17 is CYP2C19_UM.",
      source_name = "CYP2C19:RM"
    ),
    CYP2C19_UM = list(
      description = "CYP2C19 ultrarapid-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with all CYP2C19 indicators 0 the reference is the normal metabolizer (*1/*1)",
      notes = "Proportional effect (1 + 0.253)^CYP2C19_UM on CL/F (25% higher). Genotype *17/*17 only (Table S2); *1/*17 is CYP2C19_RM in this model.",
      source_name = "CYP2C19:UM"
    ),
    CYP2C19_MISSING = list(
      description = "CYP2C19 phenotype not determined",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phenotype known)",
      notes = "Proportional effect (1 - 0.14)^CYP2C19_MISSING on CL/F. 59 of 497 participants (11.9%) had no or conflicting genotype assignments (Table 2, Table S2 footnote c); the coefficient describes that cohort's phenotype mix and does not transfer to other populations.",
      source_name = "CYP2C19:missing"
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no omeprazole)",
      notes = "Proportional effect (1 - 0.332)^CONMED_OMEPRAZOLE on CL/F (33% lower; CYP2C19 inhibition). Time-varying: set on the records during which omeprazole was taken (Methods 'Covariate analysis').",
      source_name = "Omeprazole"
    ),
    CONMED_ESOMEPRAZOLE = list(
      description = "Concomitant esomeprazole",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no esomeprazole)",
      notes = "Proportional effect (1 - 0.421)^CONMED_ESOMEPRAZOLE on CL/F (42% lower; CYP2C19 inhibition). Time-varying like CONMED_OMEPRAZOLE.",
      source_name = "Esomeprazole"
    ),
    STUDY_MYK461014 = list(
      description = "Record from study MYK-461-014",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any of the other 11 studies)",
      notes = "Selects the separately estimated residual error of the food-effect / relative-bioavailability study MYK-461-014 (Table S5 'EP14'). Does not affect any PK parameter; set 0 for simulation of new subjects.",
      source_name = "Study"
    )
  )

  covariatesDataExcluded <- list(
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Added to CL/F in the SCM forward step and removed in backward elimination (Results 'Covariate analysis')."
    ),
    FORM_TABLET = list(
      description = "Tablet formulation",
      units = "(binary)",
      type = "binary",
      notes = "Tablet effect on F was added in the SCM forward step and removed in backward elimination; tablets and capsules share the reference in the final model."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mavacamten", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mavacamten", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mavacamten", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Exponential residual error with a separate magnitude for study MYK-461-014;
  # both are combined into expSdCc inside model().
  paper_specific_residual_sds <- c("expSdMyk461014")

  # Inter-individual variability on the residual error magnitude.
  paper_specific_etas <- c("etaruv")

  population <- list(
    species = "human",
    n_subjects = 497L,
    n_studies = 12L,
    n_observations = 9244L,
    age_range = "18-82 years (mean 49.4)",
    weight_range = "5th-95th percentile 58-108 kg (mean 81.1, SD 16.9)",
    sex_female_pct = 40.2,
    race_ethnicity = c(White = 79.5, Black = 5.0, Asian = 12.5, Other = 0.8, Missing = 2.2),
    disease_state = "192 healthy adults, 54 adults with nonobstructive hypertrophic cardiomyopathy and 251 adults with obstructive hypertrophic cardiomyopathy",
    dose_range = "Single oral doses 1-48 mg; multiple oral doses 2-25 mg/day (once or twice daily, titrated in the patient studies)",
    renal_function = "eGFR 29.5-148 mL/min/1.73 m^2 (mean 94.5)",
    cyp2c19 = "NM 39.8%, IM (*1/*2-type) 17.9%, RM (*1/*17) 19.1%, *2/*17 4.6%, PM 3.4%, UM (*17/*17) 3.2%, missing 11.9%",
    regions = "Not reported",
    notes = "Pooled seven phase 1 studies, two phase 2 studies, a phase 2 open-label extension, the phase 3 EXPLORER-HCM study and the MAVA-LTE extension (Chang 2024 Table 1). Baseline demographics in Table 2."
  )

  ini({
    # Structural parameters for the reference participant: male patient with
    # obstructive HCM, 84 kg, CYP2C19 NM, eGFR 95 mL/min/1.73 m^2, 5 mg
    # capsule/tablet, no PPI (Chang 2024 Table 3 / Table S5).
    lcl <- log(0.914); label("Apparent clearance CL/F (L/h)") # Table S5 'CL/F, L/h' 0.914
    lvc <- log(6.63); label("Apparent central volume V2/F (L)") # Table S5 'V2/F, L' 6.63
    lvp <- log(252); label("Apparent peripheral volume V3/F (L)") # Table S5 'V3/F, L' 252
    lq <- log(16.3); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S5 'Q/F, L/h' 16.3
    lka <- log(0.301); label("Absorption rate constant KA (1/h)") # Table S5 'KA, h-1' 0.301
    ltlag <- log(0.192); label("Absorption lag time ALAG1 (h)") # Table S5 'ALAG1, h' 0.192
    lfdepot <- fixed(log(1)); label("Relative bioavailability F (reference)") # Table S5 'F' 1 (fixed)

    # Body-weight power exponents (reference 84 kg)
    e_wt_cl <- fixed(0.75); label("Weight exponent on CL/F and Q/F (unitless)") # Table S5 'Body weight-CLx/F' 0.75 (fixed)
    e_wt_vc_vp <- 0.813; label("Weight exponent on V2/F and V3/F (unitless)") # Table S5 'Body weight-Vx/F' 0.813

    # Covariate effects on CL/F; categorical effects enter as (1 + theta)^X (Figure S2)
    e_cyp2c19_pm_cl <- -0.722; label("Fractional change in CL/F, CYP2C19 PM vs NM (unitless)") # Table S5 'CYP2C19:PM-CL/F' -0.722
    e_cyp2c19_im_cl <- -0.331; label("Fractional change in CL/F, CYP2C19 IM vs NM (unitless)") # Table S5 'CYP2C19:IM-CL/F' -0.331
    e_cyp2c19_rm_cl <- -0.00218; label("Fractional change in CL/F, CYP2C19 RM vs NM (unitless)") # Table S5 'CYP2C19:RM-CL/F' -0.00218
    e_cyp2c19_um_cl <- 0.253; label("Fractional change in CL/F, CYP2C19 UM vs NM (unitless)") # Table S5 'CYP2C19:UM-CL/F' 0.253
    e_cyp2c19_missing_cl <- -0.14; label("Fractional change in CL/F, CYP2C19 phenotype missing vs NM (unitless)") # Table S5 'CYP2C19:missing' -0.14
    e_healthy_cl <- 0.272; label("Fractional change in CL/F, healthy participant vs HCM patient (unitless)") # Table S5 'Healthy individual-CL/F' 0.272
    e_crcl_cl <- 0.205; label("Power exponent of eGFR/95 on CL/F (unitless)") # Table S5 'eGFR-CL/F' 0.205
    e_esomeprazole_cl <- -0.421; label("Fractional change in CL/F with esomeprazole (unitless)") # Table S5 'Esomeprazole-CL/F' -0.421
    e_omeprazole_cl <- -0.332; label("Fractional change in CL/F with omeprazole (unitless)") # Table S5 'Omeprazole-CL/F' -0.332
    e_sexf_cl <- 0.273; label("Fractional change in CL/F, female vs male (unitless)") # Table S5 'Female sex-CL/F' 0.273

    # Covariate effects on KA
    e_form_solution_ka <- 0.681; label("Fractional change in KA, oral solution vs tablet/capsule (unitless)") # Table S5 'Solution formulation-KA' 0.681
    e_healthy_fast_ka <- 0.873; label("Fractional change in KA, healthy fasted vs reference (unitless)") # Table S5 'Healthy individual, fasted-KA' 0.873
    e_dose_ka <- -0.135; label("Power exponent of dose/5 mg on KA (unitless)") # Table S5 'Dose-KA' -0.135

    # Covariate effect on ALAG1
    e_dose_tlag <- -0.087; label("Power exponent of dose/5 mg on ALAG1 (unitless)") # Table S5 'Dose-ALAG1' -0.087

    # Covariate effects on F
    e_form_solution_f <- -0.194; label("Fractional change in F, oral solution vs tablet/capsule (unitless)") # Table S5 'Solution formulation-F' -0.194
    e_dose_f <- 0.0777; label("Power exponent of dose/5 mg on F (unitless)") # Table S5 'Dose-F' 0.0777
    e_healthy_fast_f <- -0.314; label("Fractional change in F, healthy fasted vs reference (unitless)") # Table S5 'Healthy individual, fasted-F' -0.314
    e_healthy_fed_f <- -0.433; label("Fractional change in F, healthy fed vs reference (unitless)") # Table S5 'Healthy individual, fed-F' -0.433
    e_sexf_f <- 0.0945; label("Fractional change in F, female vs male (unitless)") # Table S5 'Female sex-F' 0.0945

    # IIV: Table S5 reports CV% = 100 * sqrt(exp(omega^2) - 1), so
    # omega^2 = log(1 + CV^2); covariances = r * omega_i * omega_j.
    # CL 57.2%, V2 217%, V3 22.9%; r(CL,V2) = -0.021, r(CL,V3) = 0.576,
    # r(V2,V3) = 0.0307 (Table S5 random effects).
    etalcl + etalvc + etalvp ~ c(
      0.283059,
      -0.0147464, 1.74203,
      0.0692825, 0.00916069, 0.0511122
    )
    etalka ~ 0.168933 # Table S5 'IIV on KA' 42.9%: log(1 + 0.429^2)
    etalq ~ 0.0528666 # Table S5 'IIV on Q/F' 23.3%: log(1 + 0.233^2)
    etaruv ~ 0.217567 # Table S5 'IIV on RUV' 49.3%: log(1 + 0.493^2)

    # Residual error on log-transformed concentrations; Table S5 EP / EP14 are
    # NONMEM SIGMA variances, so the SD is their square root.
    expSd <- 0.170880; label("Exponential residual error SD, all studies except MYK-461-014") # Table S5 'EP' 0.0292: sqrt(0.0292)
    expSdMyk461014 <- 0.691375; label("Exponential residual error SD, study MYK-461-014") # Table S5 'EP14' 0.478: sqrt(0.478)
  })

  model({
    # Reference participant: 84 kg, eGFR 95 mL/min/1.73 m^2, 5 mg dose
    wt_cl <- (WT / 84)^e_wt_cl
    wt_v <- (WT / 84)^e_wt_vc_vp
    dose_rel <- DOSE_MAVACAMTEN_MG / 5
    healthy_fasted <- DIS_HEALTHY * (1 - FED)
    healthy_fed <- DIS_HEALTHY * FED

    cl <- exp(lcl + etalcl) * wt_cl *
      (CRCL / 95)^e_crcl_cl *
      (1 + e_cyp2c19_pm_cl)^CYP2C19_PM *
      (1 + e_cyp2c19_im_cl)^CYP2C19_IM *
      (1 + e_cyp2c19_rm_cl)^CYP2C19_RM *
      (1 + e_cyp2c19_um_cl)^CYP2C19_UM *
      (1 + e_cyp2c19_missing_cl)^CYP2C19_MISSING *
      (1 + e_healthy_cl)^DIS_HEALTHY *
      (1 + e_omeprazole_cl)^CONMED_OMEPRAZOLE *
      (1 + e_esomeprazole_cl)^CONMED_ESOMEPRAZOLE *
      (1 + e_sexf_cl)^SEXF
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp + etalvp) * wt_v
    q <- exp(lq + etalq) * wt_cl
    ka <- exp(lka + etalka) *
      (1 + e_form_solution_ka)^FORM_SOLUTION *
      (1 + e_healthy_fast_ka)^healthy_fasted *
      dose_rel^e_dose_ka
    tlag <- exp(ltlag) * dose_rel^e_dose_tlag
    fdepot <- exp(lfdepot) *
      (1 + e_form_solution_f)^FORM_SOLUTION *
      dose_rel^e_dose_f *
      (1 + e_healthy_fast_f)^healthy_fasted *
      (1 + e_healthy_fed_f)^healthy_fed *
      (1 + e_sexf_f)^SEXF

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # mg / L -> ng/mL
    Cc <- central / vc * 1000
    expSdCc <- (expSd * (1 - STUDY_MYK461014) + expSdMyk461014 * STUDY_MYK461014) * exp(etaruv)
    Cc ~ lnorm(expSdCc)
  })
}
