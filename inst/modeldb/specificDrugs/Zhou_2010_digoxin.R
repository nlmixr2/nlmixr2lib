Zhou_2010_digoxin <- function() {
  description <- "One-compartment first-order oral absorption population PK model of digoxin in older Chinese patients (Zhou 2010, Acta Pharmacol Sin); concomitant spironolactone, body weight, and serum creatinine modify Cl/F via multiplicative linear-deviation terms."
  reference <- "Zhou XD, Gao Y, Guan Z, Li ZD, Li J. Population pharmacokinetic model of digoxin in older Chinese patients and its application in clinical practice. Acta Pharmacologica Sinica. 2010;31(6):753-758. doi:10.1038/aps.2010.51."
  vignette <- "Zhou_2010_digoxin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 62.9 kg (population mean, Zhou 2010 Table 1). Applied as a multiplicative linear-deviation effect on Cl/F: cl *= (1 - 0.0101 * (WT - 62.9)). Population observed range 34-91 kg.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine (baseline)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 126.8 umol/L (population mean, Zhou 2010 Table 1). Applied as a multiplicative linear-deviation effect on Cl/F: cl *= (1 - 0.00120 * (CREAT - 126.8)). Source paper column 'Cr' in umol/L; population observed range 36-686 umol/L. Patients with serious hepatic or renal dysfunction were excluded.",
      source_name        = "Cr"
    ),
    CONMED_SPIRON = list(
      description        = "Concomitant spironolactone coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant spironolactone)",
      notes              = "32 of 119 subjects (27%) in Zhou 2010 were coadministered spironolactone (Table 1). Applied as a multiplicative linear-deviation effect on Cl/F: cl *= (1 - 0.412 * CONMED_SPIRON), i.e. spironolactone reduces digoxin Cl/F by ~41%, consistent with spironolactone's inhibition of renal-tubular digoxin secretion (Zhou 2010 Discussion).",
      source_name        = "SPI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 119L,
    n_studies       = 1L,
    age_range       = "60-88 years (mean 71.0)",
    age_median      = "71 years (mean; Zhou 2010 Table 1)",
    weight_range    = "34-91 kg (mean 62.9)",
    weight_median   = "62.9 kg (mean; Zhou 2010 Table 1)",
    sex_female_pct  = 42,
    race_ethnicity  = "Han Chinese (PLA General Hospital of Air Force, Beijing)",
    disease_state   = "Older patients (>=60 years) on chronic oral digoxin (>=7 days); 113/119 (95%) with congestive heart failure of varying severity. Patients with serious hepatic or renal dysfunction were excluded.",
    dose_range      = "Oral digoxin tablets 0.25 mg/pellet (Sine Pharmaceutical, Shanghai). 173 trough-leaning samples; median 22.9 h (range 6-192 h) between last dose and phlebotomy.",
    regions         = "China (Beijing)",
    concomitant     = "SPI 32, nifedipine 53, diltiazem 1, nitrate 86, propafenone 27 (Zhou 2010 Table 1). Only SPI retained as covariate on Cl/F in the final model.",
    renal_function  = "Serum creatinine mean 126.8 umol/L (range 36-686); BUN mean 9.0 mmol/L (range 2.5-28.9). Patients with serious renal dysfunction excluded.",
    hepatic_function = "ALT mean 29.3 U/L (range 3-383); AST mean 33.9 U/L (range 8-402); ALB mean 63.3 g/L (range 31-89). Patients with serious hepatic dysfunction excluded.",
    assay           = "Fluorescence polarization immunoassay (TDx-FLx, Abbott); LLOQ 0.2 ug/L.",
    notes           = "Baseline demographics per Zhou 2010 Table 1. Retrospective single-center cohort. Digoxin serum concentration mean 1.11 ug/L (range 0.07-4.45). Validation: bootstrap intra-validation (10 groups, dOFV < 3.84) and an external 8-patient cohort with -4.3% to +25% percent difference between observed and individual-predicted concentrations."
  )

  ini({
    # Structural parameters (Zhou 2010 Table 7, final model; OFV = -11.354).
    lka <- fixed(log(1.63));      label("First-order oral absorption rate (1/h)")   # Table 7: Ka = 1.63 (fixed; absorption phase complete, fixed to avoid Ka-fluctuation instability per Results section)
    lcl <- log(5.90);             label("Apparent oral clearance Cl/F (L/h)")        # Table 7: Cl/F = 5.90 L/h (RSE 6.97%, 95% CI 5.09-6.71)
    lvc <- log(550);              label("Apparent oral volume of distribution Vd/F (L)") # Table 7: Vd/F = 550 L (RSE 19.6%, 95% CI 338-762)

    # Covariate-effect coefficients (Zhou 2010 Table 7). Sign convention follows the
    # paper's reported (1 - theta * cov) multiplicative linear-deviation form, so the
    # coefficients are stored with their paper-reported positive sign and the minus
    # sign appears explicitly in model() below.
    e_conmed_spiron_cl <- 0.412;  label("Multiplicative linear-deviation effect of spironolactone on Cl/F (unitless; applied as (1 - e * CONMED_SPIRON))") # Table 7: theta_SPI-Cl = 0.412 (RSE 26.5%, 95% CI 0.198-0.626)
    e_wt_cl <- 0.0101;            label("Multiplicative linear-deviation slope of WT on Cl/F (1/kg; applied as (1 - e * (WT - 62.9)))")                # Table 7: theta_WT-Cl = 0.0101 (RSE 120%, 95% CI -0.0136 to 0.0338)
    e_creat_cl <- 0.00120;        label("Multiplicative linear-deviation slope of CREAT on Cl/F (L/umol; applied as (1 - e * (CREAT - 126.8)))")       # Table 7: theta_Cr-Cl = 0.00120 (RSE 45.8%, 95% CI 0.000122-0.00228)

    # Inter-individual variability (Zhou 2010 Table 7 final model). Reported as CV%
    # on the linear scale; converted to log-normal variance via omega^2 = log(CV^2 + 1).
    # Ka IIV was held fixed at 0 (Zhou 2010 Results: "Ka and inter-individual variation
    # were fixed to 1.63 h-1 and 0, respectively, in the next calculation"); etalka is
    # therefore omitted from the IIV block.
    etalcl ~ 0.21528                                                                      # log(0.490^2 + 1) = 0.21528 (Table 7: 49.0% CV on Cl/F)
    etalvc ~ 0.63595                                                                      # log(0.943^2 + 1) = 0.63595 (Table 7: 94.3% CV on Vd/F)

    # Residual error: additive only in the final model (Zhou 2010 Table 7: proportional
    # CV reported as "-", indicating the proportional term was dropped during covariate
    # modelling; only the additive component was retained from the basic-model error
    # equation Y = YPRED + YPRED * ERR(1) + ERR(2)). See vignette Assumptions and
    # deviations for the basic-model error structure.
    addSd <- 0.365;               label("Additive residual standard deviation (ug/L)")     # Table 7: additive SD = 0.365 ug/L (residual variability between observed and predicted)
  })

  model({
    # Individual PK parameters. Cl/F applies the paper's multiplicative
    # linear-deviation covariate structure with reference covariate values 62.9 kg WT,
    # 126.8 umol/L CREAT, and 0 SPI. The final-model Cl/F formula reported in the
    # Discussion is:
    #   Cl/F = 5.9 * (1 - 0.412 * SPI) * (1 - 0.0101 * (WT - 62.9)) * (1 - 0.00120 * (Cr - 126.8))
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 - e_conmed_spiron_cl * CONMED_SPIRON) *
      (1 - e_wt_cl * (WT - 62.9)) *
      (1 - e_creat_cl * (CREAT - 126.8))
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, central in mg, vc in L -> central / vc has units mg/L = ug/mL;
    # to express the observation in ug/L (the source paper's unit) multiply by 1000.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
