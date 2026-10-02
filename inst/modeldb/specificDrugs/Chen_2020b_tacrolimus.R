Chen_2020b_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral tacrolimus whole-blood concentrations ",
    "in Chinese children after liver transplantation (Chen 2020). The ",
    "absorption rate constant ka is fixed at 4.48 1/h from earlier paediatric ",
    "liver-transplant tacrolimus models. Apparent oral clearance CL/F is ",
    "allometrically scaled by body weight (fixed exponent 0.75, reference ",
    "70 kg), multiplied by 1.61 in recipients carrying a CYP3A5*1 allele ",
    "(expressers) and reduced by 10.8% on concomitant Wuzhi capsule ",
    "(Schisandra sphenanthera extract). Apparent volume V/F scales linearly ",
    "with body weight (fixed exponent 1). Exponential IIV on V/F only; ",
    "combined proportional-plus-additive residual error."
  )
  reference <- paste0(
    "Chen X, Wang DD, Xu H, Li ZP. Population pharmacokinetics and ",
    "pharmacogenomics of tacrolimus in Chinese children receiving a liver ",
    "transplant: initial dose recommendation. Transl Pediatr. ",
    "2020;9(5):576-586. doi:10.21037/tp-20-84"
  )
  vignette <- "Chen_2020b_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to a standard weight of 70 kg (Chen 2020 ",
        "Methods Eq. 3; final-model Eqs. 7-8). Table 1: 13.79 +/- 5.12 kg, ",
        "median 13.00 (range 6.40-28.00) kg, so the 70 kg reference lies ",
        "well outside the observed range."
      ),
      source_name = "WT"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (1 = at least one CYP3A5*1 allele)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = paste0(
        "Chen 2020 text below Table 1: 'For patients who were CYP3A5*3/*3, ",
        "the CYP3A5 value = 0; for patients with a CYP3A5*1 allele, the ",
        "CYP3A5 value = 1', so the source column maps directly onto ",
        "CYP3A5_EXPR. Recipient (not donor) genotype, from next-generation ",
        "sequencing of residual TDM blood (PGxOne 160). Table 2: *1/*1 = 2, ",
        "*1/*3 = 5, *3/*3 = 5 of 12. Enters CL/F as 1.61^CYP3A5_EXPR ",
        "(Methods Eq. 4; Eq. 7)."
      ),
      source_name = "CYP3A5"
    ),
    CONMED_WUZHI = list(
      description = "Concomitant Wuzhi capsule indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no Wuzhi capsule)",
      notes = paste0(
        "Chen 2020 text below Table 1: 'for patients co-administered WZ, ",
        "the WZ value = 1, otherwise WZ value = 0'. Enters CL/F as the ",
        "fractional form (1 + theta_WZ * WZ) with theta_WZ = -0.108 ",
        "(Methods Eq. 6; Eq. 7 printed as (1 - 0.108 x WZ)). Table 1: 6 of ",
        "12 patients received Wuzhi capsule. The paper does not state ",
        "whether the indicator was time-varying within a patient."
      ),
      source_name = "WZ"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected (Table 1: 3.11 +/- 2.28 years, median 2.42, range 0.47-7.96) but not retained in the final model (Eqs. 7-8)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Collected (Table 1: median 39.80, range 22.20-49.00 g/L) but not retained."
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = "Collected (Table 1: median 33.71, range 15.84-43.30 %) but not retained, although haematocrit is a retained CL/F covariate in several other tacrolimus models."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Collected (Table 1: median 22.00, range 7.00-82.00 umol/L) but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "0.47-7.96 years",
    age_median = "2.42 years",
    weight_range = "6.40-28.00 kg",
    weight_median = "13.00 kg",
    sex_female_pct = 33.3,
    race_ethnicity = "Chinese (single centre, Children's Hospital of Fudan University, Shanghai).",
    disease_state = paste0(
      "Paediatric liver-transplant recipients treated with oral tacrolimus ",
      "between September 2014 and October 2019, analysed retrospectively. ",
      "Post-transplantation day median 121 (range 4-1877) days. 8 boys, ",
      "4 girls."
    ),
    dose_range = "Oral tacrolimus twice daily; the individual doses are not tabulated.",
    regions = "China (Shanghai).",
    co_medication = "Glucocorticoid 11/12, aspirin 9/12, Wuzhi capsule 6/12, fluconazole 2/12 (Table 1).",
    genotypes = "Recipient CYP3A5 *1/*1 2, *1/*3 5, *3/*3 5 (Table 2); ten other pharmacogenes were sequenced and none was retained.",
    assay = "Whole-blood tacrolimus by Emit 2000 Tacrolimus Assay (Siemens), range 2.0-30 ng/mL.",
    notes = paste0(
      "Routine therapeutic-drug-monitoring data, partly shared with two ",
      "earlier reports by the same group (Chen 2020 refs 23-24). The number ",
      "of concentrations and the sampling times are not reported. NONMEM 7 ",
      "FOCE-I; 1000-replicate bootstrap."
    )
  )

  ini({
    # Chen 2020 Table 3 'Parameter estimates of final model and bootstrap
    # validation'; final-model equations printed in Results as Eqs. 7-8:
    #   CL/F = 6.57 x (WT/70)^0.75 x 1.61^CYP3A5 x (1 - 0.108 x WZ)
    #   V/F  = 77.6 x (WT/70)
    lcl <- log(6.57); label("Apparent oral clearance CL/F at 70 kg, CYP3A5*3/*3, no Wuzhi capsule (L/h)") # Table 3 CL/F = 6.57 L/h (SE 0.414); Eq. 7
    lvc <- log(77.6); label("Apparent volume of distribution V/F at 70 kg (L)") # Table 3 V/F = 77.6 L (SE 1.198); Eq. 8
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)") # Table 3 Ka = 4.48 (fixed); Methods 'The value of Ka was fixed to 4.48/h (23,25)'

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Methods Eq. 3, POWER = 0.75 for CL/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Methods Eq. 3, POWER = 1 for V/F

    e_cyp3a5_expr_cl <- 1.61; label("CYP3A5 expresser multiplicative factor on CL/F (unitless)") # Table 3 theta CYP3A5 = 1.61 (SE 0.341); Eq. 4 and Eq. 7
    e_conmed_wuzhi_cl <- -0.108; label("Fractional change in CL/F with Wuzhi capsule (unitless)") # Table 3 theta WZ = -0.108 (SE 0.615); Eq. 6 and Eq. 7

    # Table 3 omega V/F = 0.396 is read as the SD of eta (variance
    # 0.396^2 = 0.156816). Deterministic replication of the Figure 3 Monte
    # Carlo target-attainment curves favours the SD reading over the
    # variance reading (see the vignette Assumptions section).
    etalvc ~ 0.156816 # Table 3 omega V/F = 0.396 (SD), SE 0.818

    # Table 3 sigma 1 (proportional) and sigma 2 (additive) are read as SDs,
    # the same scale as omega in the same table; Methods Eq. 2
    # C = Cpre x (1 + eps1) + eps2.
    propSd <- 0.288; label("Proportional residual error (fraction)") # Table 3 sigma 1 = 0.288 (SE 0.137)
    addSd <- 0.720; label("Additive residual error (ng/mL)") # Table 3 sigma 2 = 0.720 (SE 0.697)
  })

  model({
    ka <- exp(lka)
    # Eq. 7: CYP3A5 as theta^CYP3A5 (Eq. 4); Wuzhi capsule as (1 + theta * WZ)
    # (Eq. 6). No IIV on CL/F is reported (Table 3 lists omega V/F only).
    cl <- exp(lcl) * (WT / 70)^e_wt_cl *
      e_cyp3a5_expr_cl^CYP3A5_EXPR *
      (1 + e_conmed_wuzhi_cl * CONMED_WUZHI)
    # Eq. 8
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, vc in L -> mg/L; x 1000 gives ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
