Nel_2026_ceftriaxone <- function() {
  description <- paste(
    "Two-compartment population PK model for ceftriaxone given as a 2 g",
    "intravenous dose and as a 2 g subcutaneous infusion (50 mL normal saline",
    "into the abdominal wall over about 30 minutes) in non-critically ill adult",
    "inpatients (Nel 2026; 20 patients, self-controlled IV-then-SC cross-over,",
    "212 capillary dried-blood-spot concentrations). Subcutaneous drug enters a",
    "depot with bioavailability 0.958 and first-order absorption (ka 0.764",
    "1/h). Clearance is proportional to Cockcroft-Gault creatinine clearance",
    "(1.57 L/h per 100 mL/min, no weight term); central volume scales linearly",
    "with body weight and exponentially with age; peripheral volume and",
    "intercompartmental clearance carry fixed allometric weight scaling.",
    "Concentrations are WHOLE BLOOD (dried blood spot), not plasma.",
    sep = " "
  )
  reference <- paste(
    "Nel H, Murray F, Yoo O, Rawlins M, Raby E, Page-Sharp M, Moore BR,",
    "Salman S, Manning L. Safety, tolerability, and pharmacokinetics of a 2 g",
    "subcutaneous dose of ceftriaxone as an alternative to intravenous",
    "delivery. Antimicrob Agents Chemother. 2026;70(1):e01081-25.",
    "doi:10.1128/aac.01081-25. PMID 41277849; PMC12777559.",
    sep = " "
  )
  vignette <- "Nel_2026_ceftriaxone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Ceftriaxone was measured in capillary dried blood spots (LC-MS/MS, LOQ
  # 0.1 mg/L) and the model was fitted to those whole-blood concentrations
  # (Methods, 'The final model of DBS (whole blood) concentrations was used
  # for simulation'). Cc is therefore a whole-blood concentration.
  compartmentData <- list(
    depot = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Computed by the authors with the Cockcroft-Gault equation from the",
        "serum creatinine measured on the day of dosing (Methods, 'PK",
        "sampling'). Raw mL/min, NOT normalized to 1.73 m^2 -- the register's",
        "CRCL entry documents this Cockcroft-Gault alias. Enters clearance as",
        "the through-origin linear ratio CRCL / 100 (Results: 'CL = 1.57 x",
        "CrCl/100'; Table 2 unit 'l.h-1.(100 mL/min)-1'). The reference 100",
        "mL/min is not the cohort median (108 mL/min, IQR 71-147). Body weight",
        "does NOT additionally scale clearance in the final model: the",
        "Cockcroft-Gault estimate already carries weight.",
        sep = " "
      ),
      source_name = "CrCl"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling was applied a priori with fixed exponents of 0.75",
        "on clearance terms and 1 on volume terms (Methods). In the final",
        "model weight enters Vc and Vp (exponent 1, Table 2 units",
        "'l.70kg-1') and Q (exponent 0.75, Table 2 unit 'l.h-1.70kg-1') with",
        "reference 70 kg; it does not enter CL, which is scaled by CRCL alone.",
        "Cohort median 84.5 kg (IQR 74.3-97.8).",
        sep = " "
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the central volume as exp(0.0094 x (AGE - 60)) (Results;",
        "Table 2 row 'thetaAGE where Vc x EXP(thetaAGE*(AGE-60))'). The",
        "reference 60 years is the centring value printed by the authors;",
        "the cohort median was 63.5 years (IQR 43.5-76.5).",
        sep = " "
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not a covariate of the PK model (tested on volume, not retained).",
        "Used only outside the model, at the cohort median of 0.37, to",
        "convert simulated whole-blood concentrations to unbound plasma",
        "concentrations: Cu = fu / (H x rho + (1 - H)) x C_DBS (Methods,",
        "'Population pharmacokinetics modeling'). See",
        "population$protein_binding.",
        sep = " "
      ),
      source_name = "H"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and V in the forward-inclusion covariate search; not retained (Methods; Discussion notes few hypoalbuminemic patients). Cohort median 35 g/L (IQR 31-36).",
      source_name = "albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = 212L,
    age_range = "Median 63.5 years (IQR 43.5-76.5); adults >= 18 years (Results; Table 1).",
    weight_range = "Median 84.5 kg (IQR 74.3-97.8); BMI median 26.9 kg/m^2 (IQR 23.4-32.0) (Results).",
    sex_female_pct = 30,
    disease_state = paste(
      "Clinically stable (not in ICU/HDU and no medical emergency team call in",
      "the previous 24 h) adult inpatients already receiving IV ceftriaxone",
      "2 g for a variety of infections: genito-urinary 30%, intra-abdominal",
      "30%, primary bacteremia 10%, and one each of skin/soft tissue, bone and",
      "joint, infective endocarditis, respiratory, gynecological and CNS",
      "(Table 1). No participant had severe liver disease and few were",
      "markedly hypoalbuminemic.",
      sep = " "
    ),
    renal_function = "Cockcroft-Gault creatinine clearance median 108 mL/min (IQR 71-147); serum creatinine median 78 umol/L (IQR 62.3-88) (Results).",
    dose_range = paste(
      "Ceftriaxone 2 g IV (19 of 20 once daily, 1 twice daily), followed by a",
      "single 2 g SC dose diluted in 50 mL normal saline, gravity-fed over",
      "about 30 minutes into the lower anterior abdominal wall through a 24-gauge",
      "SC catheter, then a 10-20 mL saline flush over 5-10 minutes.",
      sep = " "
    ),
    regions = "Australia (Fiona Stanley Hospital, Western Australia)",
    hematocrit = "Median 37.5% (IQR 29.3-40) (Results); the simulation conversion used 0.37 (Methods).",
    protein_binding = paste(
      "Not part of the model. Cc is the WHOLE-BLOOD (dried blood spot)",
      "concentration. For Figures 2 and 3 the authors converted simulated",
      "whole-blood concentrations to unbound plasma concentrations with",
      "Cu = fu / (H x rho + (1 - H)) x C_DBS, H = 0.37 (cohort median",
      "hematocrit), and a red-cell partitioning ratio rho and plasma unbound",
      "fraction fu taken from Page-Sharp 2016 (Antimicrob Agents Chemother",
      "60:14-23, doi:10.1128/AAC.01740-15). Nel 2026 does not print the values",
      "of fu and rho; Page-Sharp 2016 Table 1 reports red-cell partitioning",
      "0.55-0.61 and plasma protein binding 94.5%, 95.3% and 90.1% at 10, 50",
      "and 200 mg/L. The validation vignette documents the values used to",
      "reproduce Figure 2.",
      sep = " "
    ),
    sampling = "Capillary dried blood spots before dosing and 1, 2, 4, 8 and 24 h after the start of each IV and SC dose (plus 12 h for the twice-daily patient); LC-MS/MS, LOQ 0.1 mg/L.",
    notes = "Prospective self-controlled cross-over trial (ACTRN12624000692538). Fitted in NONMEM 7.5.1 (FOCE-I, ADVAN6) with PsN/Pirana; 1,000-replicate nonparametric bootstrap. Covariate search by forward inclusion (P < 0.05) and backward elimination (P < 0.01)."
  )

  ini({
    lcl <- log(1.57)
    label("Clearance at CRCL = 100 mL/min (L/h)")
    # Table 2 'CL, l.h-1.(100 mL/min)-1' = 1.57 (RSE 10%); Results 'CL = 1.57 x CrCl/100'.

    lvc <- log(17.6)
    label("Central volume at WT = 70 kg and AGE = 60 years (L)")
    # Table 2 'Vc, l.70kg-1 for a 60 years old' = 17.6 (RSE 6%).

    lvp <- log(12.9)
    label("Peripheral volume at WT = 70 kg (L)")
    # Table 2 'VP, l.70kg-1' = 12.9 (RSE 21%).

    lq <- log(0.917)
    label("Intercompartmental clearance at WT = 70 kg (L/h)")
    # Table 2 'Q, l.h-1.70kg-1' = 0.917 (RSE 12%).

    lka <- log(0.764)
    label("Subcutaneous absorption rate constant (1/h)")
    # Table 2 'KASC, h-1' = 0.764 (RSE 13%).

    lfdepot <- log(0.958)
    label("Absolute subcutaneous bioavailability (fraction)")
    # Table 2 'FSC' = 0.958 (RSE 2%); the abstract quotes the bootstrap median 95.7%.

    e_age_vc <- 0.0094
    label("Exponential age coefficient on Vc, per year from 60 years (1/year)")
    # Table 2 'thetaAGE where Vc x EXP(thetaAGE*(AGE-60))' = 0.0094 (RSE 33%).

    e_crcl_cl <- fixed(1)
    label("Exponent on (CRCL / 100) for CL (unitless)")
    # Results 'CL = 1.57 x CrCl/100': linear, no exponent estimated.

    e_wt_vc <- fixed(1)
    label("Allometric exponent on (WT / 70) for Vc (unitless)")
    # Methods: a priori allometric exponent 1 on volumes; Results 'Vc = 17.6 x (WT/70) x ...'.

    e_wt_vp <- fixed(1)
    label("Allometric exponent on (WT / 70) for Vp (unitless)")
    # Methods: a priori allometric exponent 1 on volumes; Table 2 unit 'l.70kg-1'.

    e_wt_q <- fixed(0.75)
    label("Allometric exponent on (WT / 70) for Q (unitless)")
    # Methods: a priori allometric exponent 0.75 on clearances; Table 2 unit 'l.h-1.70kg-1'.

    # Table 2 footnote: IIV and residual error are presented as
    # 100 x sqrt(variance), so omega^2 = (value / 100)^2.
    etalcl ~ 0.205209 # Table 2 'IIV CL' = 45.3 (shrinkage 3%) -> 0.453^2
    etalvc ~ 0.0576 # Table 2 'IIV VC' = 24 (shrinkage 11%) -> 0.24^2

    propSd <- 0.157
    label("Proportional residual error (fraction)")
    # Table 2 'Proportional error (%)' = 15.7; Results 'overall mean residual variability of the model was 15.7%'.
  })

  model({
    # Reference values printed in Results / Table 2.
    crcl_ref <- 100 # mL/min
    wt_ref <- 70 # kg
    age_ref <- 60 # years

    # CL = 1.57 x CrCl/100 (no weight term)
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    # Vc = 17.6 x (WT/70) x exp(0.0094 x (AGE - 60))
    vc <- exp(lvc + etalvc) *
      (WT / wt_ref)^e_wt_vc *
      exp(e_age_vc * (AGE - age_ref))
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vp
    q <- exp(lq) * (WT / wt_ref)^e_wt_q
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # SC infusion enters the depot; IV doses go directly to central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Whole-blood (dried blood spot) ceftriaxone concentration.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
