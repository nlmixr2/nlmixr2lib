Chen_2021_sirolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral sirolimus whole-blood trough concentrations in Chinese pediatric ",
    "patients with tuberous sclerosis complex (TSC)-related epilepsy (Chen ",
    "2021). The absorption rate constant ka is fixed at 0.485 1/h from the ",
    "literature because only trough concentrations were available. Apparent ",
    "oral clearance CL/F is allometrically scaled by body weight (fixed ",
    "exponent 0.75, reference 70 kg) and multiplied by 1.16 with concomitant ",
    "oxcarbazepine (power form 1.16^CONMED_OXC). Apparent volume V/F is ",
    "scaled linearly by body weight (fixed exponent 1). Exponential IIV on ",
    "CL/F only; additive residual error."
  )
  reference <- paste0(
    "Chen X, Wang D, Zhu L, Lu J, Huang Y, Wang G, Zhu Y, Ye Q, Wang Y, Xu H, ",
    "Li Z. Population Pharmacokinetics and Initial Dose Optimization of ",
    "Sirolimus Improving Drug Blood Level for Seizure Control in Pediatric ",
    "Patients With Tuberous Sclerosis Complex. Front Pharmacol. ",
    "2021;12:647232. doi:10.3389/fphar.2021.647232."
  )
  vignette <- "Chen_2021_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to 70 kg (Chen 2021 Methods Eq. 5 and ",
        "Results Eqs. 8-9). Table 1: 23.50 +/- 11.71 kg, median 20.50 ",
        "(range 8.00-68.00) kg."
      ),
      source_name = "weight"
    ),
    CONMED_OXC = list(
      description = "Concomitant oxcarbazepine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant oxcarbazepine)",
      notes = paste0(
        "1 = co-administered oxcarbazepine, 0 = not (Chen 2021 Table 2 ",
        "footnote: 'when patients received oxcarbazepine, OXC was 1, ",
        "otherwise OXC was 0'). Power-form effect on CL/F, 1.16^OXC (Results ",
        "Eq. 8), i.e. 16% higher clearance on oxcarbazepine, attributed to ",
        "CYP3A4 induction. Table 1: 23 of 80 patients received ",
        "oxcarbazepine. The paper does not state whether the indicator was ",
        "time-varying within a patient."
      ),
      source_name = "OXC"
    )
  )

  covariatesDataExcluded <- list(
    DOSE_FORM = list(
      description = "Sirolimus dosage form (tablet vs oral solution)",
      units = "(categorical)",
      type = "categorical",
      notes = paste0(
        "Screened and not retained (Chen 2021 Results, Modeling; Discussion: ",
        "'sirolimus dosage forms, tablet or solution, had no significant ",
        "effect on sirolimus clearance rate'). 41 person-times each."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    n_concentrations = 188L,
    age_range = "0.61-16.61 years (median 5.76; mean +/- SD 6.35 +/- 3.83)",
    weight_range = "8.00-68.00 kg (median 20.50; mean +/- SD 23.50 +/- 11.71)",
    sex_female_pct = 56.25,
    race_ethnicity = c(Asian = 100),
    disease_state = "Tuberous sclerosis complex (TSC)-related epilepsy, pediatric",
    dose_range = "Oral sirolimus (tablet or oral solution), once daily, titrated by TDM; individual doses not tabulated",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    sampling_design = paste0(
      "Retrospective routine TDM, May 2016 - October 2020; 188 whole-blood ",
      "trough concentrations (2.35 per patient) measured by Emit 2000 ",
      "immunoassay (linear range 3.5-30 ng/mL)."
    ),
    notes = paste0(
      "35 boys and 45 girls. Estimated in NONMEM 7 by FOCE-I; bootstrap ",
      "(1,000 replicates) and pcVPC. Co-medications (Table 1): valproic acid ",
      "40, oxcarbazepine 23, vigabatrin 12, levetiracetam 10, topiramate 8, ",
      "lamotrigine 4, carbamazepine 2. Only oxcarbazepine was retained."
    )
  )

  ini({
    # Chen 2021 Table 2 (final model), Results Eqs. 8-9
    lka <- fixed(log(0.485)); label("Absorption rate constant ka (1/h)")  # Table 2 Ka = 0.485 (fixed), from Wang 2020 (Methods, Population Pharmacokinetic Model)
    lcl <- log(8.59); label("Apparent oral clearance CL/F at 70 kg without oxcarbazepine (L/h)")  # Table 2 CL/F = 8.59 L/h, SE 0.251; Eq. 8
    lvc <- log(294); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 2 V/F = 294 L, SE 0.806; Eq. 9

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods Eq. 5, power = 0.75 for CL/F; Eq. 8
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods Eq. 5, power = 1 for V/F; Eq. 9

    e_oxc_cl <- 1.16; label("Multiplicative factor on CL/F with concomitant oxcarbazepine (unitless)")  # Table 2 theta OXC = 1.16, SE 0.062; Eq. 8 1.16^OXC

    # Table 2 reports omega CL/F = 0.175 as the SD of eta (variance =
    # 0.175^2 = 0.030625). The SD reading reproduces the probability-of-target
    # percentages printed in Figures 2-5; the variance reading flattens them
    # (see the vignette Assumptions section). IIV on V/F was dropped (Results,
    # Evaluation).
    etalcl ~ 0.030625   # Table 2 omega CL/F = 0.175 (SD), SE 0.232

    addSd <- 1.913; label("Additive residual error (ng/mL)")  # Table 2 sigma 1 = 1.913, additive error, SE 0.065
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * e_oxc_cl^CONMED_OXC
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
