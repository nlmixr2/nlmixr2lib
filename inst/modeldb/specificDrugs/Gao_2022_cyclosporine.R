Gao_2022_cyclosporine <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral cyclosporine A in Chinese children with acquired aplastic anemia (Gao 2022), with fixed allometric body-weight scaling (reference 28 kg) on all clearances and volumes, a linear total-bilirubin effect on CL/F, peripheral volume and intercompartmental clearance fixed from the literature, and log-normal residual error on whole-blood concentrations."
  reference <- "Gao X, Bian ZL, Qiao XH, Qian XW, Li J, Shen GM, Miao H, Yu Y, Meng JH, Zhu XH, Jiang JY, Le J, Yu L, Wang HS, Zhai XW. Population Pharmacokinetics of Cyclosporine in Chinese Pediatric Patients With Acquired Aplastic Anemia. Front Pharmacol. 2022;13:933739. doi:10.3389/fphar.2022.933739"
  vignette <- "Gao_2022_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "cyclosporine", units = "mg", specimen = "administration site", verified = TRUE),
    # Gao 2022 Methods 'Pharmacokinetic Sampling and Measurement': whole-blood
    # concentrations measured with the Emit 2000 cyclosporine-specific assay.
    central = list(analyte = "cyclosporine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "cyclosporine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling (WT / 28)^0.75 on CL/F and Q/F and",
        "(WT / 28)^1 on Vc/F and Vp/F (Gao 2022 Eqs 2-3 and Table 2 footnote).",
        "28 kg is the stated median body weight used for centring (Methods);",
        "Table 1 prints the cohort median as 27.5 kg (range 12.0-91.0)."
      ),
      source_name = "BW"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on CL/F, centred at 10.65 umol/L: CL = CLtypical *",
        "(1 + (TBIL - 10.65) * -0.0107) (Table 2 footnote, which drops the",
        "'1 +'; see the vignette). Table 1 median 10.90 umol/L (range",
        "4.00-70.00). The Abstract and Discussion say 'per 1 nmol/L', a",
        "misprint for umol/L (Table 1 units)."
      ),
      source_name = "TBIL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 157L,
    n_studies = 1L,
    n_observations = 681L,
    age_range = "1.5-18.8 years",
    age_median = "7.8 years",
    weight_range = "12.0-91.0 kg",
    weight_median = "27.5 kg",
    race_ethnicity = "Chinese",
    disease_state = "Acquired aplastic anemia (59.9% non-severe, 25.4% severe, 13.4% very severe); stem-cell transplant recipients excluded",
    dose_range = "Oral cyclosporine 5 mg/kg/day twice daily initially, titrated to 5-8 (up to 10) mg/kg/day to trough 100-200 ng/mL and peak 300-400 ng/mL",
    regions = "China (Children's Hospital of Fudan University and Tongji Hospital of Tongji University, Shanghai)",
    co_medication = "Testosterone undecanoate 70.7%, methylprednisolone 29.9%, prednisone 15.3%, rabbit ATG 29.9%, G-CSF 17.2% (Table 1); none retained as covariates",
    notes = "Retrospective TDM data, 2014-2021: sparse whole-blood peak (2-4 h post-dose) and pre-dose troughs. Baseline demographics in Gao 2022 Table 1; sex distribution not reported."
  )

  ini({
    # Structural parameters at the 28-kg reference child, TBIL = 10.65 umol/L
    lka <- log(1.26)
    label("Absorption rate constant (Ka, 1/h)") # Table 2 'Ka (/h)' 1.26 (RSE 20.9%)
    lcl <- log(29.1)
    label("Apparent clearance at WT = 28 kg, TBIL = 10.65 umol/L (CL/F, L/h)") # Table 2 'CL/F (L/h)' 29.1 (RSE 3.8%)
    lvc <- log(325)
    label("Apparent central volume at WT = 28 kg (Vc/F, L)") # Table 2 'VC/F (L)' 325 (RSE 15.7%)
    lq <- fixed(log(3.1))
    label("Apparent intercompartmental clearance at WT = 28 kg (Q/F, L/h)") # Table 2 'Q/F (L/h)' 3.1 FIX; Results: fixed from Eljebari 2012
    lvp <- fixed(log(262))
    label("Apparent peripheral volume at WT = 28 kg (Vp/F, L)") # Table 2 'VP/F (L)' 262 FIX; Results: fixed from Eljebari 2012

    # Allometric exponents (fixed)
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of WT on CL/F and Q/F (unitless)") # Eq 2; Table 2 footnote
    e_wt_vc <- fixed(1)
    label("Allometric exponent of WT on Vc/F and Vp/F (unitless)") # Eq 3; Table 2 footnote

    # Covariate effect on CL/F
    e_tbili_cl <- -0.0107
    label("Linear slope of (TBIL - 10.65 umol/L) on CL/F (1/(umol/L))") # Table 2 'TBIL on CL (%)' -1.07 (RSE 20.4%); footnote slope -0.0107

    # IIV: omega^2 = log(CV^2 + 1), Table 2 footnote CV formula
    etalcl ~ 0.0754785 # Table 2 CL/F 'CV for IIV' 28.0%; log(0.280^2 + 1)
    etalvc ~ 0.3261628 # Table 2 VC/F 'CV for IIV' 62.1%; log(0.621^2 + 1)

    # Residual error: additive on log-transformed concentrations
    expSd <- 0.348
    label("Additive residual error on the log scale (SD)") # Table 2 'sigma' 0.348 (RSE 10.1%)
  })

  model({
    # TBILI centring value 10.65 umol/L from the Table 2 footnote equation
    cl <- exp(lcl + etalcl) * (WT / 28)^e_wt_cl * (1 + e_tbili_cl * (TBILI - 10.65))
    vc <- exp(lvc + etalvc) * (WT / 28)^e_wt_vc
    q <- exp(lq) * (WT / 28)^e_wt_cl
    vp <- exp(lvp) * (WT / 28)^e_wt_vc
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 -> ng/mL (paper units)
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
