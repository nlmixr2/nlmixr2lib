Chen_2020_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral tacrolimus whole-blood concentrations in Chinese children and ",
    "adolescents with lupus nephritis (Chen 2020). The absorption rate ",
    "constant ka is fixed at 4.48 1/h from the literature. Apparent oral ",
    "clearance CL/F is allometrically scaled by body weight (fixed exponent ",
    "0.75, reference 70 kg) and reduced by 29% with concomitant Wuzhi capsule; ",
    "apparent volume V/F is scaled linearly by body weight. Exponential IIV on ",
    "CL/F only; proportional residual error."
  )
  reference <- paste0(
    "Chen X, Wang DD, Xu H, Li ZP. Population pharmacokinetics model and ",
    "initial dose optimization of tacrolimus in children and adolescents with ",
    "lupus nephritis based on real-world data. Exp Ther Med. ",
    "2020;20(2):1423-1430. doi:10.3892/etm.2020.8821."
  )
  vignette <- "Chen_2020_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to 70 kg (Chen 2020 Methods equation iii and ",
        "Results equations vi-vii). Table I: 45.89 +/- 10.55 kg, median ",
        "47.00 (range 17.00-66.50) kg."
      ),
      source_name = "WT"
    ),
    CONMED_WUZHI = list(
      description = "Concomitant Wuzhi capsule indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no Wuzhi capsule)",
      notes = paste0(
        "1 = co-administered Wuzhi capsule, 0 = not (Chen 2020 Results, text ",
        "below equation vii). Fractional effect (1 + theta * WZ) on CL/F with ",
        "theta = -0.290, i.e. CL/F ratio 1:0.71. Table I: 12 of 32 patients ",
        "received Wuzhi capsule. The paper does not state whether the ",
        "indicator was time-varying within a patient."
      ),
      source_name = "WZ"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "2.86-17.99 years (median 13.87; mean +/- SD 13.44 +/- 2.86)",
    weight_range = "17.00-66.50 kg (median 47.00; mean +/- SD 45.89 +/- 10.55)",
    sex_female_pct = 84.4,
    race_ethnicity = c(Asian = 100),
    disease_state = "Lupus nephritis, pediatric and adolescent patients",
    dose_range = "Oral tacrolimus, twice daily, titrated by TDM (daily dose not tabulated)",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    sampling_design = paste0(
      "Retrospective routine TDM of tacrolimus whole-blood concentrations ",
      "(about 2-13 ng/mL in Figure 1); number of samples not reported."
    ),
    notes = paste0(
      "Data August 2014 - September 2019; 5 males and 27 females. All 32 ",
      "patients received a glucocorticoid and 12 received Wuzhi capsule ",
      "(Table I). Estimated in NONMEM 7 by FOCE-I; bootstrap (n = 1000) and ",
      "pcVPC. Some patients overlap with earlier cohorts from the same centre."
    )
  )

  ini({
    # Chen 2020 Table II (final model), equations vi-vii
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)")  # Table II Ka = 4.48 (fixed); Methods, refs 13, 16-18
    lcl <- log(15.5); label("Apparent oral clearance CL/F at 70 kg without Wuzhi capsule (L/h)")  # Table II CL/F = 15.5 L/h, SE 0.729
    lvc <- log(174); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table II V/F = 174 L, SE 1.908

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods equation iii, power = 0.75 for CL/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods equation iii, power = 1 for V/F

    e_conmed_wuzhi_cl <- -0.290; label("Fractional change in CL/F with Wuzhi capsule (unitless)")  # Table II theta WZ = -0.290, SE 0.390

    # Table II reports omega CL/F = 0.172 without saying whether it is a
    # variance or an SD. Read as the SD of eta (variance 0.172^2 = 0.029584):
    # that reading reproduces the Figure 4 probabilities of target attainment,
    # the variance reading does not (see the vignette Assumptions section).
    etalcl ~ 0.029584 # Table II omega CL/F = 0.172 (SD)

    propSd <- 0.281; label("Proportional residual error (fraction)")  # Table II sigma1 = 0.281, proportional error, SE 0.049
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_conmed_wuzhi_cl * CONMED_WUZHI)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
