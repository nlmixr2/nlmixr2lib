Friedrich_2014_evacetrapib <- function() {
  description <- paste0(
    "Two-compartment population PK model with first-order absorption for ",
    "oral evacetrapib, a cholesteryl ester transfer protein (CETP) ",
    "inhibitor, given once daily as monotherapy or with atorvastatin, ",
    "simvastatin or rosuvastatin to adults with hypercholesterolemia or low ",
    "HDL-C (Friedrich 2014). Apparent clearance rises linearly with the ",
    "evacetrapib dose (less than dose-proportional exposure) and with ",
    "Cockcroft-Gault creatinine clearance; the absorption rate constant is ",
    "fixed. Correlated IIV on CL/F and Vc/F, independent IIV on Q/F, and a ",
    "proportional residual error. The steady-state AUC this model gives ",
    "(1000 x DOSE / CL/F) drives the paper's two exposure-response models, ",
    "Friedrich_2014_evacetrapib_hdl.R and Friedrich_2014_evacetrapib_ldl.R."
  )
  reference <- paste0(
    "Friedrich S, Kastelein JJP, James D, Waterhouse T, Nissen SE, ",
    "Nicholls SJ, Krueger KA. The pharmacokinetics and ",
    "pharmacokinetic/pharmacodynamic relationships of evacetrapib ",
    "administered as monotherapy or in combination with statins. CPT ",
    "Pharmacometrics Syst Pharmacol. 2014;3(1):e94. doi:10.1038/psp.2013.70. ",
    "Final estimates: Table 1 and Supplementary Data ",
    "'Friedrich_PK_model_code_and_output' (NONMEM 7.2 $PK and final ",
    "estimates)."
  )
  vignette <- "Friedrich_2014_evacetrapib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "evacetrapib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "evacetrapib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "evacetrapib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    DOSE = list(
      description = "Assigned once-daily evacetrapib dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Use case (a) of the DOSE register entry: the subject's fixed ",
        "assigned daily dose (30, 100 or 500 mg in the source study). ",
        "Enters CL/F linearly as (1 + e_dose_cl * DOSE) with no centring ",
        "(Table 1 footnote c). The supplementary NONMEM $PK names the data ",
        "item DDI; Table 1 footnote c defines it as 'the evacetrapib dose ",
        "in mg'. The less-than-proportional rise in exposure is attributed ",
        "by the authors to a falling extent of absorption, which in this ",
        "parameterisation appears as a rising apparent clearance."
      ),
      source_name = "DOSE (DDI in the NONMEM $PK)"
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance, raw (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Raw Cockcroft-Gault estimate in mL/min, NOT normalised to ",
        "1.73 m^2 (the paper reports CGCL in ml/minute). Enters CL/F as ",
        "the linear centred term (1 + e_crcl_cl * (CRCL - 91.28)), 91.28 ",
        "mL/min being the population median CGCL (Table 1 footnote c). The ",
        "body-weight convention inside the Cockcroft-Gault equation is not ",
        "stated. Model-predicted CL/F is about 10% lower at 50 than at ",
        "100 mL/min; few patients had CGCL below 50 mL/min (Supplementary ",
        "Figure S4), so the model does not describe moderate or severe ",
        "renal impairment."
      ),
      source_name = "CGCL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the PK parameters (Methods) and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the PK parameters (Methods) and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on the PK parameters (Methods) and not retained."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the PK parameters (Methods) as gender and not retained."
    ),
    CONMED_STATIN = list(
      description = "Coadministered statin (atorvastatin 20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily)",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Each statin was screened on the PK parameters (Methods) and none ",
        "was retained; the Discussion states that statin coadministration ",
        "had no significant effect on evacetrapib clearance."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 226L,
    n_observations = 1629L,
    n_studies = 1L,
    age_range = "28-83 years (whole randomised study)",
    age_mean = "58.3 years (whole randomised study)",
    weight_range = "44-164 kg (whole randomised study)",
    weight_mean = "83.6 kg (whole randomised study)",
    sex_female_pct = 56.0,
    race_ethnicity = c(White = 92.9),
    disease_state = paste0(
      "Adults with hypercholesterolemia or low HDL-C (triglycerides below ",
      "400 mg/dL after lipid washout and dietary lead-in); mean baseline ",
      "LDL-C 144.3, HDL-C 55.1 and triglycerides 139.8 mg/dL."
    ),
    dose_range = paste0(
      "Evacetrapib 30, 100 or 500 mg orally once daily for 12 weeks as ",
      "monotherapy, or 100 mg once daily with atorvastatin 20 mg, ",
      "simvastatin 40 mg or rosuvastatin 10 mg once daily."
    ),
    renal_function = "Median Cockcroft-Gault creatinine clearance 91.28 mL/min (Table 1 footnote c).",
    co_medication = "Atorvastatin 20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily in the combination arms.",
    notes = paste0(
      "Phase II, randomised, double-blind, placebo- and active-controlled ",
      "parallel-group study (398 randomised, about 40 per arm). The PK ",
      "dataset held 1,629 concentrations from 226 patients (1-9 samples ",
      "each); the supplementary NONMEM output reports 233 patients in the ",
      "data file. Sampling: predose and 1-2 h postdose at week 2; predose ",
      "and 3-18 h postdose at weeks 4, 8 and 12; one follow-up sample 4-6 ",
      "weeks after the last dose. LLOQ 1 ng/mL. Demographics are given for ",
      "the whole randomised population in the Results text; no ",
      "PK-subset demographic table is published."
    )
  )

  ini({
    # Supplementary Data 'Friedrich_PK_model_code_and_output' (NONMEM 7.2,
    # ADVAN4 TRANS4, FOCE-I) gives the final estimates to more digits than
    # Table 1 and is the source used here; every value agrees with Table 1.
    lka <- fixed(log(0.3)); label("Absorption rate constant (1/h)") # Table 1 'Ka' = 0.3 hour-1 (fixed); NONMEM THETA(1) 0.300 FIXED. Results: Ka could not be estimated from the sparse early data, 0.3 gave the best objective function of the values tried
    lcl <- log(15.3); label("Apparent clearance CL/F at DOSE = 0 and CRCL = 91.28 mL/min (L/h)") # Table 1 'CL/F' = 15.3 l/hour (%SEE 4.57); NONMEM THETA(2) 15.3
    lvc <- log(228); label("Apparent central volume Vc/F (L)") # Table 1 'V2/F' = 228 l (%SEE 6.36); NONMEM THETA(3) 228
    lq <- log(9.79); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Q/F' = 9.79 l/hour (%SEE 24.9); NONMEM THETA(4) 9.79
    lvp <- log(827); label("Apparent peripheral volume Vp/F (L)") # Table 1 'V3/F' = 827 l (%SEE 18); NONMEM THETA(5) 827

    e_dose_cl <- 0.00105; label("Linear effect of evacetrapib dose on CL/F (1/mg)") # Table 1 'Effect of dose on CL/F' = 0.00105 (%SEE 21.6); NONMEM THETA(6); footnote c: CL/F = 15.3 * (1 + 0.00105 * DOSE) * (...)
    e_crcl_cl <- 0.00198; label("Linear effect of CGCL on CL/F, centred at 91.28 mL/min (1/(mL/min))") # Table 1 'Effect of CGCL on CL/F' = 0.00198 (%SEE 34.1); NONMEM THETA(7); footnote c: (1 + 0.00198 * (CGCL - 91.28))

    # NONMEM $OMEGA BLOCK(2): 0.144, 0.225, 0.426. Table 1 reports these as
    # %CV via 100*sqrt(exp(OMEGA)-1): 39.4% (CL/F), 72.9% (V2/F); the same
    # formula applied to the covariance 0.225 gives the printed '50.2%'
    # (correlation 0.225 / sqrt(0.144 * 0.426) = 0.91).
    etalcl + etalvc ~ c(0.144, 0.225, 0.426) # NONMEM OMEGA(2,2) 0.144, OMEGA(3,2) 0.225, OMEGA(3,3) 0.426; Table 1 IIV 39.4%, covariance 50.2%, 72.9%
    etalq ~ 1.04 # NONMEM OMEGA(4,4) 1.04; Table 1 Q/F IIV 135.2% CV. No IIV on Ka or V3/F (OMEGA held at 0)

    propSd <- 0.3302; label("Proportional residual error (fraction)") # NONMEM SIGMA 0.109 (variance) -> sqrt(0.109) = 0.3302; Table 1 'Residual error (proportional, %)' = 33 (%SEE 6.63)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 + e_dose_cl * DOSE) * (1 + e_crcl_cl * (CRCL - 91.28))
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; the NONMEM scale S2 = V2/1000
    # reports ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
