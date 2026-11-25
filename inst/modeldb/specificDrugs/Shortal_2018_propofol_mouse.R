Shortal_2018_propofol_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Three-compartment pharmacokinetic model for",
    "intravenous propofol in adult male C57BL/6 mice (Shortal 2018), with",
    "compartment 1 = blood, compartment 2 = whole brain tissue (serial",
    "cortical biopsies) and compartment 3 = all other tissues; elimination",
    "from blood only. Parameterised by first-order rate constants and a",
    "blood volume scaling constant (sinf, L/kg). Fitted by unweighted least",
    "squares to blood and brain concentrations from 21 mice and used to",
    "drive brain target-controlled infusion (TCI). No between-animal or",
    "residual variability was estimated; the model is deterministic."
  )
  reference <- paste(
    "Shortal BP, Reitz SL, Aggarwal A, Meng QC, McKinstry-Wu AR, Kelz MB,",
    "Proekt A (2018). Development and validation of brain target controlled",
    "infusion of propofol in mice. PLoS One 13(4):e0194949.",
    "doi:10.1371/journal.pone.0194949.",
    "Authors' MATLAB TCI code: doi:10.5281/zenodo.1205136.",
    sep = " "
  )
  vignette <- "Shortal_2018_propofol_mouse"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the blood (compartment 1) volume, vc = sinf * WT. The paper",
        "dosed in mg/kg/min and divided the infusion rate by sinf (Eq 1), so",
        "for a dose given in mg (rate in mg/min = mg/kg/min * WT) the",
        "predicted concentrations do not depend on WT. Study mice weighed",
        "27.4 +/- 1.5 g (0.0274 kg) at surgery."
      ),
      source_name = "body weight"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "whole blood", verified = TRUE),
    brain = list(analyte = "propofol", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6, male)",
    n_subjects = 21L,
    n_studies = 1L,
    age_range = "adult; 2.90 +/- 0.45 months (mean +/- SD)",
    weight_range = "27.4 +/- 1.5 g (mean +/- SD)",
    sex_female_pct = 0,
    disease_state = paste(
      "healthy; jugular venous catheter, bilateral craniotomy and durotomy",
      "under isoflurane, reduced to a sub-anaesthetic 0.2 percent during",
      "propofol infusion"
    ),
    dose_range = paste(
      "IV propofol: 15 mg/kg bolus (150 mg/kg/min for 6 s), 2 mg/kg/min for",
      "1 h, and two generations of brain TCI targeting 10 ug/g"
    ),
    regions = "USA (University of Pennsylvania)",
    notes = paste(
      "Model fitted to 241 samples (193 brain, 48 blood) from 21 mice in",
      "four infusion paradigms (Results; Figs 3-4). Validated prospectively",
      "in 10 further mice (68 brain, 12 blood samples) under TCI targeting",
      "10 and 15 ug/g brain (Fig 5)."
    )
  )

  ini({
    # Table 1 'Constants of final model' (min^-1). The authors' TCI code
    # (TCIMaster2.m, doi:10.5281/zenodo.1205136) carries the same final
    # constants in 1/s: k = [0.0258 0.0452 0.0037 0.0007 0.0012 2.8373] for
    # [k12 k21 k13 k31 k10 kInf], i.e. 1.548, 2.712, 0.222, 0.042, 0.072
    # 1/min and kInf = 1/sinf = 1/0.3524 kg/L; Table 1 rounds these.
    lk_central_brain <- log(1.55); label("Blood-to-brain rate constant k12 (1/min)") # Table 1 'Blood to brain, k12' = 1.55
    lk_brain_central <- log(2.71); label("Brain-to-blood rate constant k21 (1/min)") # Table 1 'Brain to blood, k21' = 2.71
    lk12 <- log(0.22); label("Blood-to-other-tissues rate constant k13 (1/min)") # Table 1 'Blood to other tissues, k13' = 0.22
    lk21 <- log(0.04); label("Other-tissues-to-blood rate constant k31 (1/min)") # Table 1 'Other tissues to blood, k31' = 0.04
    lkel <- log(0.07); label("Elimination rate constant from blood k10 (1/min)") # Table 1 'Elimination, k10' = 0.07
    lvc <- log(0.35); label("Blood (compartment 1) volume per kg body weight sinf (L/kg)") # Table 1 'Scaling constant for infusion (L kg-1), sinf' = 0.35

    # Unweighted least-squares fit (Methods 'Fitting procedure'); no residual
    # variance was estimated or reported.
    addSd <- fixed(0); label("Additive residual error, blood (mg/L); not estimated")
    addSd_Cbrain <- fixed(0); label("Additive residual error, brain (mg/L, i.e. ug/g); not estimated")
  })

  model({
    k_central_brain <- exp(lk_central_brain)
    k_brain_central <- exp(lk_brain_central)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    vc <- exp(lvc) * WT

    # Eq 1 is written for concentrations x1, x2, x3 with input I(t)/sinf.
    # Multiplying through by vc = sinf * WT gives these amount-scale states,
    # with x_i = state_i / vc for every compartment, so the brain and
    # other-tissue concentrations are also obtained by dividing by vc.
    d/dt(central) <- -(k_central_brain + k12 + kel) * central +
      k_brain_central * brain + k21 * peripheral1
    d/dt(brain) <- k_central_brain * central - k_brain_central * brain
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cbrain <- brain / vc

    Cc ~ add(addSd)
    Cbrain ~ add(addSd_Cbrain)
  })
}
