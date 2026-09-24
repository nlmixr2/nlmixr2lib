Wang_2019_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral tacrolimus whole-blood trough concentrations in Chinese pediatric ",
    "patients with systemic-onset juvenile idiopathic arthritis (Wang 2019). ",
    "The absorption rate constant ka is fixed at 4.48 1/h from the literature ",
    "because only trough concentrations were available. Apparent oral ",
    "clearance CL/F is allometrically scaled by body weight (fixed exponent ",
    "0.75, reference 70 kg) and reduced by three fractional co-medication ",
    "effects (omeprazole -36.2%, loratadine -32.2%, diltiazem -30.7%). ",
    "Apparent volume V/F is scaled linearly by body weight (fixed exponent ",
    "1). Exponential IIV on CL/F only; additive residual error."
  )
  reference <- paste0(
    "Wang D, Chen X, Xu H, Li Z. Population pharmacokinetics of tacrolimus ",
    "in pediatric patients with systemic-onset juvenile idiopathic ",
    "arthritis: Initial dosage recommendations. Exp Ther Med. ",
    "2019;18(6):4653-4660. doi:10.3892/etm.2019.8129."
  )
  vignette <- "Wang_2019_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to 70 kg (Wang 2019 Methods equation iii and ",
        "Results equations vi-vii). Table I: 29.83 +/- 10.66 kg, median ",
        "33.60 (range 13.50-46.00) kg."
      ),
      source_name = "Weight"
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no omeprazole)",
      notes = paste0(
        "1 = co-administered omeprazole, 0 = not (Wang 2019 Results, text ",
        "below equation vii). Fractional effect (1 + theta * CONMED) on CL/F ",
        "with theta = -0.362. Table II: 6 of 17 patients received omeprazole. ",
        "The paper does not state whether the indicator was time-varying ",
        "within a patient."
      ),
      source_name = "omeprazole"
    ),
    CONMED_LORATADINE = list(
      description = "Concomitant loratadine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no loratadine)",
      notes = paste0(
        "1 = co-administered loratadine, 0 = not. Fractional effect on CL/F ",
        "with theta = -0.322. Table II: 5 of 17 patients received loratadine."
      ),
      source_name = "loratadine"
    ),
    CONMED_DILTIAZEM = list(
      description = "Concomitant diltiazem indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diltiazem)",
      notes = paste0(
        "1 = co-administered diltiazem, 0 = not. Fractional effect on CL/F ",
        "with theta = -0.307. Table II: 3 of 17 patients received diltiazem. ",
        "Drug-specific indicator; distinct from the class-level CONMED_CCB."
      ),
      source_name = "diltiazem"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 17L,
    n_studies = 1L,
    n_concentrations = 86L,
    age_range = "3.20-14.60 years (median 9.50; mean +/- SD 8.23 +/- 3.30)",
    weight_range = "13.50-46.00 kg (median 33.60; mean +/- SD 29.83 +/- 10.66)",
    sex_female_pct = 52.9,
    race_ethnicity = c(Asian = 100),
    disease_state = "Systemic-onset juvenile idiopathic arthritis (SOJIA), pediatric (<16 years)",
    dose_range = "Oral tacrolimus; daily dose 1.00-4.00 mg (median 1.50; mean 1.66 +/- 0.71), titrated by TDM",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    sampling_design = paste0(
      "Retrospective routine TDM; 86 whole-blood trough concentrations ",
      "(1.3-9.2 ng/mL) measured by Emit 2000 immunoassay."
    ),
    notes = paste0(
      "Data January 2014 - December 2017; 8 males and 9 females. Estimated in ",
      "NONMEM 7 by FOCE-I; bootstrap (988/1000 successful) and pcVPC. ",
      "Common co-medications (Table II): prednisone 15/17, methotrexate 8/17, ",
      "cefdinir 6/17, methylprednisolone 5/17, loratadine 5/17, omeprazole ",
      "6/17, cefprozil 4/17, cefixime 3/17, diltiazem 3/17."
    )
  )

  ini({
    # Wang 2019 Table III (final model), equations vi-vii
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)")  # Table III Ka = 4.480 (fixed), from Yang 2015 (ref 15)
    lcl <- log(29.7); label("Apparent oral clearance CL/F at 70 kg without co-medication (L/h)")  # Table III CL/F = 29.700 L/h, SE 9.3%
    lvc <- log(1120); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table III V/F = 1120.000 L, SE 27.9%

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods equation iii, COE = 0.75 for CL/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods equation iii, COE = 1 for V/F

    e_conmed_omeprazole_cl <- -0.362; label("Fractional change in CL/F with omeprazole (unitless)")  # Table III theta omeprazole = -0.362, SE 16.8%
    e_conmed_loratadine_cl <- -0.322; label("Fractional change in CL/F with loratadine (unitless)")  # Table III theta loratadine = -0.322, SE 23.8%
    e_conmed_diltiazem_cl <- -0.307; label("Fractional change in CL/F with diltiazem (unitless)")  # Table III theta diltiazem = -0.307, SE 34.2%

    # Table III reports omega CL/F = 0.265 as the SD of eta (variance =
    # 0.265^2 = 0.070225). The SD reading reproduces the Table IV simulated
    # 95% intervals; the variance reading gives intervals ~3x too wide (see
    # the vignette Assumptions section).
    etalcl ~ 0.070225   # Table III omega CL/F = 0.265 (SD), SE 18.4%

    addSd <- 1.229; label("Additive residual error (ng/mL)")  # Table III sigma1 = 1.229, additive error, SE 5.1%
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_conmed_omeprazole_cl * CONMED_OMEPRAZOLE) *
      (1 + e_conmed_loratadine_cl * CONMED_LORATADINE) *
      (1 + e_conmed_diltiazem_cl * CONMED_DILTIAZEM)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
