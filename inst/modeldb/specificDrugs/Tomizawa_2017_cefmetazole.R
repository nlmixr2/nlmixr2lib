Tomizawa_2017_cefmetazole <- function() {
  description <- paste(
    "One-compartment IV population PK model for cefmetazole given as",
    "intraoperative antimicrobial prophylaxis to adults undergoing surgery",
    "for colorectal cancer (Tomizawa 2017). Clearance is directly",
    "proportional to Cockcroft-Gault creatinine clearance (raw mL/min) and",
    "volume is directly proportional to total body weight, both through the",
    "origin with no intercept; exponential IIV on CL and V and a",
    "proportional residual error."
  )
  reference <- paste(
    "Tomizawa A, Nakamura T, Komatsu T, Inano H, Kondo R, Watanabe M,",
    "Atsuda K.",
    "Optimal dosage of cefmetazole for intraoperative antimicrobial",
    "prophylaxis in patients undergoing surgery for colorectal cancer.",
    "J Pharm Health Care Sci. 2017;3:1.",
    "doi:10.1186/s40780-016-0071-6."
  )
  vignette <- "Tomizawa_2017_cefmetazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefmetazole", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance by the Cockcroft-Gault equation, raw mL/min",
        "(NOT BSA-normalized to mL/min/1.73 m^2)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source symbol Ccr ('Ccr was calculated using the Cockcroft-Gault",
        "equation', Methods, Covariate analysis). Enters CL as the",
        "through-origin product CL (L/h) = theta1 * Ccr (Table 3), with no",
        "reference value or intercept. Cohort mean 73.9 mL/min (SD 21.7,",
        "range 47.2-126.3; Table 1). The weight used inside Cockcroft-Gault",
        "is not stated (total body weight assumed)."
      ),
      source_name = "Ccr"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source symbol Bw. Enters V as the through-origin product",
        "Vd (L) = theta2 * BW (Table 3). Cohort mean 63.7 kg (SD 9.9, range",
        "47.5-89.0; Table 1)."
      ),
      source_name = "Bw"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL (Age >= 65 vs < 65, dOBJ -0.014) and Vd (dOBJ +1.116); not retained (Table 2).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source coded Gender male = 1, female = 0 (SEXF = 1 - Gender).",
        "Tested on CL (dOBJ -4.371) and Vd (dOBJ +0.010); not retained",
        "(Table 2)."
      ),
      source_name = "Gender"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL as theta1 + theta2 / Scr (dOBJ -4.347); not retained (Table 2).",
      source_name = "Scr"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Source reports g/dL. Tested on CL dichotomised at 3.8 g/dL (dOBJ 0); not retained (Table 2).",
      source_name = "Alb"
    ),
    TUMOR_STAGE = list(
      description = "Clinical pathological stage of the colorectal cancer (I-IV)",
      units = "(ordinal)",
      type = "categorical",
      reference_category = NULL,
      notes = "Tested on CL as theta1 + theta2 * (1 + (4 - stage)) (dOBJ -1.251); not retained (Table 2).",
      source_name = "Stage"
    ),
    SURG_LAPAROSCOPIC = list(
      description = "Laparoscopic (vs open) surgery indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Tested on CL as theta1 * theta2^Procedure (dOBJ -0.231); not retained (Table 2).",
      source_name = "Procedure"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23,
    n_studies = 1,
    n_observations = 86,
    age_range = "41-84 years (mean 69, SD 10)",
    weight_range = "47.5-89.0 kg (mean 63.7, SD 9.9)",
    sex_female_pct = 100 * 5 / 23,
    race_ethnicity = "Japanese (single-centre study in Japan; not tabulated)",
    disease_state = paste(
      "Colorectal cancer (colon 14, rectum 9; stage I/II/III/IV",
      "10/6/6/1) undergoing elective surgery (laparoscopic 13, open 10);",
      "emergency surgery, ileus and dialysis excluded."
    ),
    dose_range = paste(
      "Cefmetazole 1 g IV over 5-10 min after induction of anaesthesia and",
      "within 60 min before incision, then 1 g every 3 h intraoperatively."
    ),
    regions = "Japan (Kitasato University Hospital)",
    renal_function = "Cockcroft-Gault CrCl 47.2-126.3 mL/min (mean 73.9, SD 21.7); serum creatinine 0.57-1.3 mg/dL",
    notes = paste(
      "Table 1. Operation time 140-430 min (mean 238). Samples at the start",
      "of surgery, on completion of the anastomosis, immediately before",
      "additional doses and after abdominal closure; HPLC assay, LOD",
      "0.5 ug/mL. NONMEM VI, ADVAN1 TRANS2. November 2008 to December 2010."
    )
  )

  ini({
    lcl <- log(0.0704); label("CL per unit Cockcroft-Gault CrCl (L/h per mL/min)") # Table 3: CLCMZ = theta1 x Ccr, theta1 = 0.0704 (Table 4 SE 0.0029)
    lvc <- log(0.163); label("V per kg total body weight (L/kg)") # Table 3: VdCMZ = theta2 x BW, theta2 = 0.163 (Table 4 SE 0.0054)

    # Table 3 prints omega as CV% (21.0, 8.4) and Table 4 prints the same
    # omegas as 0.210 and 0.084, i.e. the SD of eta; variance = omega^2.
    etalcl ~ 0.0441 # Table 3 'omegaCL (%)' 21.0; Table 4 omegaCL 0.210 -> 0.210^2
    etalvc ~ 0.007056 # Table 3 'omegaVd (%)' 8.4; Table 4 omegaVd 0.084 -> 0.084^2

    propSd <- 0.135; label("Proportional residual error (fraction)") # Table 3 sigma 13.5%; Table 4 sigma 0.135 (Eq. 2 C = Cpred x (1 + eps))
  })
  model({
    cl <- exp(lcl + etalcl) * CRCL
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    d / dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
