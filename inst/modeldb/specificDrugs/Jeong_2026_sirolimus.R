Jeong_2026_sirolimus <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for whole-blood sirolimus in paediatric patients (1 month to 18 years) with vascular anomalies (Jeong 2026). Clearance, intercompartmental clearance and both volumes are allometrically scaled to a 70 kg subject (exponent 0.75 on CL and Q, 1 on Vc and Vp); clearance also carries a fixed sigmoid postmenstrual-age maturation function (TM50 62.9 weeks, Hill 2.94) and a linear haemoglobin effect centred on 12.6 g/dL (higher haemoglobin gives lower CL/F). Between-subject variability on CL/F and Vc/F; residual error additive on the log scale."
  reference <- paste(
    "Jeong S, Kim J, Jung W, Jeon J, Jo H, Ryu KB, Hong JY, Kim KS, Yang S,",
    "Ju HY, Yun HY, Chae JW, Yoo KH, Lee S. Optimal dosing for vascular",
    "anomalies paediatric patients with population pharmacokinetic model of",
    "sirolimus. Br J Clin Pharmacol. 2026;92(9):3115.",
    "doi:10.1002/bcp.70595.",
    sep = " "
  )
  vignette <- "Jeong_2026_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Sirolimus was measured in WHOLE BLOOD (Methods: 'Whole blood samples for
  # PK analysis were collected using EDTA tubes to determine blood sirolimus
  # concentration'; LC-MS/MS, LLOQ 1.0 ng/mL).
  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "sirolimus", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to a 70 kg reference: exponent 0.75 on CL/F and Q/F, exponent 1 on Vc/F and Vp/F (Results final-model equation; Figure S1). Cohort mean 37 +/- 24.7 kg (range 6.3-83.8 kg; Table 1).",
      source_name = "Weight"
    ),
    PAGE = list(
      description = "Postmenstrual age.",
      units = "weeks (this model; canonical PAGE is months)",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the sigmoid maturation of CL/F, MF = PAGE^Hill / (TM50^Hill + PAGE^Hill), with TM50 = 62.9 weeks and Hill = 2.94 (both fixed). The paper had no gestational-age data and derived PMA from chronological age assuming a gestational age of 40 weeks (Results), i.e. PAGE = 40 + postnatal age in weeks. Supply PAGE in WEEKS.",
      source_name = "PMA"
    ),
    HGB = list(
      description = "Blood haemoglobin concentration.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL/F, CL multiplied by (1 - 0.155 * (HGB - 12.6)) (main-text final-model equation; Table 2 COVHGB = -0.155). The centring value 12.6 g/dL is the cohort mean (Table 1: 12.6 +/- 1.2 g/dL, range 10.2-14.6). Figure S1 prints the same equation centred on 12.9 g/dL instead; see the vignette. The linear factor reaches zero at HGB = 19.05 g/dL, so the model must not be used above the observed range. Treated as a per-subject baseline value (the paper does not say whether it was time-varying).",
      source_name = "HGB"
    )
  )

  # Covariates screened in the stepwise / bootstrap covariate search but not
  # retained (Results: 'No significant effects were observed for covariates
  # such as BIL, CYP genotype, FORM and DTT'). Documentation only.
  covariatesDataExcluded <- list(
    TBILI = list(
      description = "Total bilirubin.",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested on CL/F and Vc/F; not retained. Cohort 0.4 +/- 0.3 (0.2-1.4) (Table 1)."
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 genotype (expresser status).",
      units = "unitless",
      type = "categorical",
      notes = "CYP3A5 genotype and a CYP3A5 allele categorisation (CAT3A5), as well as CYP3A4 genotype (*1/*1 n = 25, *1/*18 n = 3), were tested on CL/F and not retained. CYP3A5 *1/*1 n = 1, *1/*3 n = 11, *3/*3 n = 16 (Table 1)."
    ),
    FORM_SIROLIMUS_SYRUP = list(
      description = "Sirolimus oral-syrup formulation indicator (vs tablet).",
      units = "unitless",
      type = "categorical",
      notes = "Formulation (FORM) and duration of treatment (DTT; cohort 593.9 +/- 771.6 days) were tested on Vc/F and not retained (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 2L,
    n_observations = 107L,
    age_range = "0.3-18.3 years (mean 8.7 +/- 5.6); 4 infants, 15 children, 9 adolescents",
    weight_range = "6.3-83.8 kg (mean 37 +/- 24.7)",
    sex_female_pct = 60.7,
    race_ethnicity = "Not reported; single-centre Korean cohort (Samsung Medical Center, Seoul).",
    disease_state = "Paediatric vascular anomalies (vascular tumours, including kaposiform haemangioendothelioma, and vascular malformations such as lymphatic and venous malformations) treated with off-label oral sirolimus.",
    dose_range = "0.27-1.04 mg/m2 or 0.01-0.04 mg/kg per dose, repeated dosing (Discussion); syrup or tablet.",
    regions = "Republic of Korea.",
    haemoglobin = "12.6 +/- 1.2 g/dL (10.2-14.6) (Table 1).",
    co_medication = "1 of 28 patients received a concomitant drug (Table 1).",
    assay = "LC-MS/MS in whole blood; calibration 1.0-26.3 ng/mL; LLOQ 1.0 ng/mL.",
    therapeutic_window = "Whole-blood steady-state trough 5-15 ng/mL (target used for the dosing simulations).",
    notes = "Pooled prospective (7 patients, 74 samples with intensive early post-dose sampling) and retrospective TDM (21 patients, trough samples on days 3, 7 and about 1 month) data; 107 samples in total (Results; Table 1). The structural parameters were estimated with a NONMEM $PRIOR NWPRI informed by the paediatric/adult two-compartment model of Miyazaki 2024 (doi:10.1016/j.dmpk.2024.101024), from which TM50 and Hill were also fixed. NONMEM 7.5, FOCE-I."
  )

  ini({
    # Final-model estimates: Jeong 2026 Table 2 ('Estimate' column). All
    # structural values are apparent (/F) values for a 70 kg subject.
    lka <- log(0.175); label("Absorption rate constant, Ka (1/h)") # Table 2 Ka = 0.175 1/h (RSE 56.6%)
    lcl <- log(14.7); label("Apparent clearance CL/F for a 70 kg mature subject at HGB 12.6 g/dL (L/h)") # Table 2 CL/F = 14.7 L/h/70 kg (RSE 9%)
    lvc <- log(102); label("Apparent central volume Vc/F for a 70 kg subject (L)") # Table 2 Vc/F = 102 L/70 kg (RSE 92%)
    lvp <- log(725); label("Apparent peripheral volume Vp/F for a 70 kg subject (L)") # Table 2 Vp/F = 725 L/70 kg (RSE 28.8%)
    lq <- log(24.7); label("Apparent intercompartmental clearance Q/F for a 70 kg subject (L/h)") # Table 2 Q/F = 24.7 L/h/70 kg (RSE 37.1%)

    # Allometric exponents (theory-based, not estimated): main-text final-model
    # equation and Figure S1 -- (Weight/70)^0.75 on CL and Q, (Weight/70)^1 on
    # Vc and Vp.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F and Q/F (unitless)") # Results equation: (Weight/70)^0.75 on CL and Q
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on Vc/F and Vp/F (unitless)") # Results equation and Figure S1: (Weight/70) on Vc and Vp

    # Sigmoid postmenstrual-age maturation of CL (fixed to literature values).
    ltm50_cl <- fixed(log(62.9)); label("Postmenstrual age at 50% CL maturation, TM50 (weeks)") # Table 2 TM50 = 62.9 weeks, fixed
    lhill_cl <- fixed(log(2.94)); label("Hill coefficient of the CL maturation function (unitless)") # Table 2 Hill = 2.94, fixed

    # Linear haemoglobin effect on CL: CL * (1 - 0.155 * (Hgb - 12.6)).
    e_hgb_cl <- -0.155; label("Fractional change in CL/F per g/dL of haemoglobin above 12.6 g/dL (1/(g/dL))") # Table 2 COVHGB = -0.155 (RSE 28.1%); Results equation

    # Between-subject variability: Table 2 OMEGA rows, read as NONMEM
    # variances (exponential IIV, Results equation e^eta).
    etalcl ~ 0.124 # Table 2 OMEGA 'IIV_CL/F' = 0.124 (RSE 15.5%)
    etalvc ~ 1.9 # Table 2 OMEGA 'IIV_Vc/F' = 1.9 (RSE 9.7%)

    # Residual error. Concentrations were natural-log transformed (Methods),
    # so the 'proportional' residual is additive on the log scale. Table 2
    # SIGMA 'RV (proportional)' = 0.13 is read as the NONMEM variance, so the
    # log-scale SD is sqrt(0.13) = 0.3606.
    expSd <- sqrt(0.13); label("Residual error SD on the natural-log concentration scale") # Table 2 SIGMA 'RV (proportional)' = 0.13 (variance; RSE 19.5%)
  })

  model({
    # Reference haemoglobin: main-text final-model equation (Hgb - 12.6),
    # equal to the cohort mean in Table 1.
    ref_hgb <- 12.6

    # Maturation function (Results equation): MF = PMA^Hill / (TM50^Hill + PMA^Hill)
    tm50_cl <- exp(ltm50_cl)
    hill_cl <- exp(lhill_cl)
    mf <- PAGE^hill_cl / (tm50_cl^hill_cl + PAGE^hill_cl)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * mf * (1 + e_hgb_cl * (HGB - ref_hgb))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
