Preijers_2021_factor_viii <- function() {
  description <- "Two-compartment population PK model for factor VIII (FVIII) concentrates in severe and moderate hemophilia A patients (children and adults) undergoing surgery (Preijers 2021). The model was re-estimated on the pooled Hazendonk 2016 perioperative dataset (119 patients, Netherlands) plus 87 children undergoing central-venous-access-device surgery at Great Ormond Street Hospital, London. PK parameters are allometrically scaled to a 68 kg reference body weight with fixed exponents of 0.75 on clearances and 1.0 on volumes; typical CL, V1, Q, V2 at 68 kg are 171 mL/h, 2930 mL, 172 mL/h, 1810 mL. CL carries a power age effect (exponent -0.12, centered at 40 years) and +14% for blood group O; V1 carries a power age effect (exponent -0.09). The measured endogenous FVIII level is added to the model prediction, and samples from patients on the B-domain-deleted product Refacto AF read 30% lower. IIV on CL (39.6%) and V1 (27.5%) is correlated (r = 0.566). Residual error is combined additive + proportional and differs by treatment-centre group (centres 1-3, 4-5 and 6)."
  reference <- "Preijers T, Liesner R, Hazendonk HCAM, Chowdary P, Driessens MHE, Hart DP, Laros-van Gorkom BAP, van der Meer FJM, Meijer K, Fijnvandraat K, Leebeek FWG, Mathot RAA, Cnossen MH; OPTI-CLOT study group. Validation of a perioperative population factor VIII pharmacokinetic model with a large cohort of pediatric hemophilia A patients. Br J Clin Pharmacol. 2021;87(11):4408-4420. doi:10.1111/bcp.14864. PMID:33884664."
  vignette <- "Preijers_2021_factor_viii"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  compartmentData <- list(
    central = list(analyte = "factor viii", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "factor viii", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling at reference 68 kg with fixed exponents 0.75 on CL and Q and 1.0 on V1 and V2 (Preijers 2021 Methods Eq. 2 and Table 2 footnote). Total-cohort median 30 kg (range 4-111 kg; Table 1). Body weight missing for 10 GOSH children was imputed from age in the source analysis (Supplemental Table S1). Treated as time-fixed per surgical procedure.",
      source_name = "BW"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect centered at 40 years: exponent -0.12 on CL and -0.09 on V1 (Preijers 2021 Methods Eq. 5, Table 2 final model and footnote). Total-cohort median 7.79 years (range 0.03-77.6; Table 1). Must be strictly positive (the power form diverges at AGE = 0).",
      source_name = "AGE"
    ),
    BLOOD_GROUP_O = list(
      description = "ABO blood group O indicator (1 = O, 0 = non-O: A, B, or AB)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-O ABO blood group: A, B, or AB)",
      notes = "Preijers 2021 Table 2 footnote: CL multiplied by 1.14^BG, i.e. 14% higher CL in blood group O. Blood group available for 175 of 206 patients; 39% blood group O (Table 1). Time-fixed per subject.",
      source_name = "BG"
    ),
    FORM_FVIII_BDD = list(
      description = "B-domain-deleted recombinant FVIII product indicator (1 = Refacto AF / moroctocog alfa; 0 = full-length recombinant or plasma-derived FVIII)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (full-length recombinant or plasma-derived FVIII)",
      notes = "Preijers 2021 Methods Eq. 1: measured FVIII = (predicted + endogenous baseline) reduced by the fraction theta_prod = 0.30 (Table 2 'B-domain deleted recombinant factor VIII') when the dichotomous covariate theta_Refacto AF = 1 (muroctocog alfa, i.e. Refacto AF). Per-observation indicator; a subject may have received different products across procedures.",
      source_name = "Refacto AF"
    ),
    FVIII_BL = list(
      description = "Measured endogenous (untreated) baseline plasma FVIII activity of the patient",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Preijers 2021 Methods Eq. 1: C_base,i, 'the measured endogenous FVIII level', is added to the model-predicted FVIII level before the BDD-product correction. Severe hemophilia A patients have < 0.01 IU/mL (83% of the total cohort) and moderate patients 0.01 to < 0.05 IU/mL (Table 1). Set to 0 to simulate the exogenous (drug-attributable) FVIII only. Time-fixed per subject.",
      source_name = "Cbase"
    ),
    STUDY_OPTICLOT_CTR45 = list(
      description = "Treatment-centre-group indicator: 1 = sample from OPTI-CLOT Netherlands centre 4 or 5; 0 otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Netherlands centres 1, 2, 3 when STUDY_GOSH is also 0)",
      notes = "Selects the centre-4/5 residual error (additive 0.06 IU/mL, proportional 21%; Preijers 2021 Table 2 final model). Affects residual error only, not the typical-value prediction. Mutually exclusive with STUDY_GOSH.",
      source_name = "Centres 4,5"
    ),
    STUDY_GOSH = list(
      description = "Treatment-centre indicator: 1 = sample from Great Ormond Street Hospital, London (centre 6, the pediatric validation cohort); 0 otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Netherlands centres 1, 2, 3 when STUDY_OPTICLOT_CTR45 is also 0)",
      notes = "Selects the centre-6 residual error (additive 0.17 IU/mL, proportional 22%; Preijers 2021 Table 2 final model). Affects residual error only. Mutually exclusive with STUDY_OPTICLOT_CTR45.",
      source_name = "Centre 6"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 206L,
    n_studies = 2L,
    n_surgeries = 342L,
    age_range = "0.03-77.6 years (new GOSH cohort 0.03-15.2; original Hazendonk 2016 cohort 0.24-77.6)",
    age_median = "7.79 years total cohort (new cohort 2.57; original cohort 39.6)",
    weight_range = "4-111 kg (new cohort 4-57; original cohort 5-111)",
    weight_median = "30 kg total cohort (new cohort 14; original cohort 75)",
    sex_female_pct = 0,
    race_ethnicity = "not reported; Netherlands and United Kingdom cohorts",
    disease_state = "Severe (83%, FVIII < 0.01 IU/mL) or moderate hemophilia A undergoing elective minor (72%) or major (28%) surgery. All 87 new-cohort children had severe hemophilia A and minor surgery to insert, replace or remove a central venous access device.",
    dose_range = "Perioperative FVIII replacement by bolus (66% of occasions) or continuous infusion (34%); median 10 doses per occasion (range 2-50).",
    regions = "Netherlands (5 hemophilia treatment centres, OPTI-CLOT) and United Kingdom (Great Ormond Street Hospital, London)",
    products = "Recombinant FVIII (88% of procedures: Advate, Recombinate, Kogenate FS, Refacto AF, Helixate FS, Octanate, Nuwiq, Innovate) or plasma-derived FVIII (12%)",
    notes = "Pooled dataset of the 119 patients / 197 procedures used to build Hazendonk 2016 plus 87 GOSH children / 145 procedures (Preijers 2021 Table 1). 2092 FVIII level measurements in total (one-stage clotting assay); no samples below the 0.01 IU/mL quantification limit. Hemophilia A is X-linked recessive; the cohort is taken as all-male (sex not tabulated)."
  )

  ini({
    # Structural PK parameters at the 68 kg reference (Preijers 2021 Table 2,
    # 'Final model' column, printed in mL and mL/h per 68 kg; stored in L and
    # L/h as in the Hazendonk 2016 sibling model).
    lcl <- log(0.171); label("Clearance CL at 68 kg, 40 years, non-O (L/h)")       # Table 2 final model: CL = 171 mL/h/68 kg (RSE 7%)
    lvc <- log(2.930); label("Central volume V1 at 68 kg, 40 years (L)")           # Table 2 final model: V1 = 2930 mL/68 kg (RSE 4%)
    lq  <- log(0.172); label("Inter-compartmental clearance Q2 at 68 kg (L/h)")    # Table 2 final model: Q2 = 172 mL/h/68 kg (RSE 19%)
    lvp <- log(1.810); label("Peripheral volume V2 at 68 kg (L)")                  # Table 2 final model: V2 = 1810 mL/68 kg (RSE 10%)

    # Fractional reduction of measured FVIII for B-domain-deleted Refacto AF
    # (theta_prod in Methods Eq. 1).
    theta_bdp <- 0.30; label("Fractional under-read of measured FVIII for B-domain-deleted Refacto AF (unitless)") # Table 2 final model: 'B-domain deleted recombinant factor VIII' = 0.30 (RSE 14%)

    # Covariate effects (Table 2 final model; equation form in Table 2 footnote)
    e_age_cl   <- -0.12; label("Power exponent of (AGE/40) on CL (unitless)")                    # Table 2 final model: CL - Age = -0.12 (RSE 26%)
    e_age_vc   <- -0.09; label("Power exponent of (AGE/40) on V1 (unitless)")                    # Table 2 final model: V1 - Age = -0.09 (RSE 24%)
    e_blood_cl <-  0.14; label("Fractional increase in CL for blood group O (unitless)")          # Table 2 final model: CL - Blood group O = 14% (RSE 6%); footnote 1.14^BG

    # Allometric exponents fixed a priori (Methods, after Eq. 2)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/68) on CL and Q (unitless)")  # Methods: 'fixed ... to 0.75 for all clearance parameters (CL, Q2)'
    e_wt_vc <- fixed(1.00); label("Allometric exponent of (WT/68) on V1 and V2 (unitless)") # Methods: 'fixed to 1 in case of a volume parameter (V1, V2)'

    # IIV (exponential). Table 2 reports %CV; omega^2 = log(1 + CV^2).
    # Covariance = r * sqrt(omega2_CL * omega2_V1) with r = 0.566.
    etalcl + etalvc ~ c(0.14568,
                        0.05833, 0.07290) # Table 2 final model: IIV CL 39.6 %CV, IIV V1 27.5 %CV, correlation CL-V1 56.6 %

    # Residual error, combined additive + proportional, by centre group
    # (Table 2 final model; Centres 1-5 Netherlands, Centre 6 GOSH London).
    addSd_center123  <- 0.12; label("Additive residual SD, centres 1, 2, 3 (IU/mL)")        # Table 2 final model: additive SD Centres 1,2,3 = 0.12 IU/mL (RSE 13%)
    addSd_center45   <- 0.06; label("Additive residual SD, centres 4, 5 (IU/mL)")           # Table 2 final model: additive SD Centres 4,5 = 0.06 IU/mL (RSE 24%)
    addSd_center6    <- 0.17; label("Additive residual SD, centre 6 (GOSH) (IU/mL)")        # Table 2 final model: additive SD Centre 6 = 0.17 IU/mL (RSE 24%)
    propSd_center123 <- 0.197; label("Proportional residual SD, centres 1, 2, 3 (fraction)") # Table 2 final model: proportional Centres 1,2,3 = 19.7 %CV (RSE 11%)
    propSd_center45  <- 0.21; label("Proportional residual SD, centres 4, 5 (fraction)")    # Table 2 final model: proportional Centres 4,5 = 0.21 (RSE 8%)
    propSd_center6   <- 0.22; label("Proportional residual SD, centre 6 (GOSH) (fraction)") # Table 2 final model: proportional Centre 6 = 0.22 (RSE 12%)
  })
  model({
    # Table 2 footnote:
    #   CL (mL/h) = 171 * (BW/68)^0.75 * (AGE/40)^-0.12 * 1.14^BG
    #   V1 (mL)   = 2930 * (BW/68)^1.0 * (AGE/40)^-0.09
    # 1.14^BG equals (1 + 0.14 * BG) for a binary BG.
    ws <- WT / 68
    cl <- exp(lcl + etalcl) * ws^e_wt_cl * (AGE / 40)^e_age_cl * (1 + e_blood_cl * BLOOD_GROUP_O)
    vc <- exp(lvc + etalvc) * ws^e_wt_vc * (AGE / 40)^e_age_vc
    q  <- exp(lq) * ws^e_wt_cl
    vp <- exp(lvp) * ws^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment IV model; bolus or continuous-infusion doses enter central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Methods Eq. 1 (read as a multiplicative BDD correction; see vignette):
    #   C_FVIII = (C_PRED + C_base) * (1 - theta_prod * RefactoAF)
    # Dose in IU and vc in L give IU/L; divide by 1000 for IU/mL.
    Cc <- (central / vc / 1000 + FVIII_BL) * (1 - theta_bdp * FORM_FVIII_BDD)

    # Centre-group residual error; the two indicators are mutually exclusive
    # and both 0 selects centres 1, 2, 3.
    ctr123 <- 1 - STUDY_OPTICLOT_CTR45 - STUDY_GOSH
    addSd_sel  <- addSd_center123 * ctr123 + addSd_center45 * STUDY_OPTICLOT_CTR45 + addSd_center6 * STUDY_GOSH
    propSd_sel <- propSd_center123 * ctr123 + propSd_center45 * STUDY_OPTICLOT_CTR45 + propSd_center6 * STUDY_GOSH
    Cc ~ add(addSd_sel) + prop(propSd_sel)
  })
}
