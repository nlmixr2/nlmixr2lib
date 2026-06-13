Bogason_2011_daunorubicin <- function() {
  description <- "Two-compartment population PK model for daunorubicin (DNR) in adults with acute myeloid leukaemia, with baseline white blood cell count as a covariate on central volume of distribution (Bogason 2011)"
  reference <- "Bogason A, Quartino AL, Lafolie P, Masquelier M, Karlsson MO, Paul C, Gruber A, Vitols S. Inverse relationship between leukaemic cell burden and plasma concentrations of daunorubicin in patients with acute myeloid leukaemia. Br J Clin Pharmacol. 2011;71(4):514-521. doi:10.1111/j.1365-2125.2010.03894.x"
  vignette <- "Bogason_2011_daunorubicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WBC = list(
      description        = "Baseline white blood cell count at AML diagnosis (peripheral blood). Used here as a biomarker of leukaemic cell burden and entered as an additive-fractional centred deviation on central volume of distribution.",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline-only (single measurement at diagnosis, time-fixed per subject). The paper reports WBC as 10^6 cells/mL blood, which is numerically equivalent to 10^9 cells/L (the canonical unit in inst/references/covariate-columns.md). Reference value 39 (10^9 cells/L) is the cohort mean / median used in the published covariate equation (Bogason 2011 Table 2 footnote).",
      source_name        = "WBC"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 40,
    n_studies        = 1,
    age_range        = "33-83 years",
    age_mean         = "61.3 years",
    weight_range     = "Not reported individually (dose 60 mg/m^2 normalised to body surface area)",
    weight_median    = "Not reported",
    sex_female_pct   = 52.5,
    race_ethnicity   = "Not reported (single-centre Swedish cohort, Karolinska University Hospital)",
    disease_state    = "Acute myeloid leukaemia (33 de novo AML, 7 secondary AML following myelodysplastic syndrome or prior chemo/radiotherapy); 26 FAB-classified subtypes (M1-M6) and 14 unclassified",
    dose_range       = "60 mg/m^2 daunorubicin as 1-h IV infusion on days 1-3 (induction); 3 patients received reduced doses (one at 30, one at 45, one at 54 mg/m^2). Cytarabine 200 mg/m^2 2-h IV infusion days 1-7 was co-administered immediately before each DNR infusion.",
    regions          = "Sweden (Karolinska University Hospital, Huddinge and Solna)",
    wbc_range        = "1-219 x 10^9 cells/L (mean 39 x 10^9 cells/L) at diagnosis",
    notes            = "Baseline demographics per Bogason 2011 Patients section and Table 1. Plasma DNR concentrations were measured at end of 1-h infusion (1 h), 5 h, and 24 h after the first DNR infusion in all 40 patients (sparse sampling). The PK model was fitted to log-transformed DNR plasma concentrations only (DOL metabolite was reported separately and not included in the final structural model). Baseline WBC was the only covariate retained after stepwise covariate analysis on CL, Vc, Q, and Vp against WBC, age, gender, body weight, and body surface area."
  )

  ini({
    # Structural PK parameters - final estimates from Bogason 2011 Table 2,
    # "Population PK model with baseline WBC" column (the final published model).
    lcl  <- log(114);  label("Clearance CL (L/h)")                          # Bogason 2011 Table 2
    lvc  <- log(412);  label("Central volume of distribution Vc (L)")       # Bogason 2011 Table 2
    lq   <- log(254);  label("Intercompartmental clearance Q (L/h)")        # Bogason 2011 Table 2
    lvp  <- log(1120); label("Peripheral volume of distribution Vp (L)")    # Bogason 2011 Table 2

    # Covariate effect of baseline WBC on Vc (additive-fractional, centred at WBC = 39).
    # Vc,i = Vc * [1 + e_wbc_vc * (WBC - 39)]   -- Bogason 2011 Table 2 footnote.
    e_wbc_vc <- 0.0138; label("Additive-fractional effect of (WBC - 39) on Vc (per 10^9 cells/L)")  # Bogason 2011 Table 2

    # Inter-individual variability on CL and Vc (final model only - paper did not
    # retain IIV on Q or Vp). Paper reports IIV as %CV on linear scale, log-normal
    # variance model. omega^2 = log(CV^2 + 1):
    #   omega^2(CL) = log(0.51^2 + 1) = 0.23116
    #   omega^2(Vc) = log(1.20^2 + 1) = 0.89200
    # Correlation between CL and Vc random effects on log scale = 0.65
    # (Bogason 2011 Table 2). Covariance = 0.65 * sqrt(0.23116) * sqrt(0.89200) = 0.29508.
    etalcl + etalvc ~ c(0.23116,
                        0.29508, 0.89200)  # Bogason 2011 Table 2 (51% CV CL, 120% CV Vc, 65% correlation)

    # Residual error - proportional only, applied on log-transformed concentrations
    # (Methods: "plasma concentrations were log-transformed ... proportional error model").
    propSd <- 0.517; label("Proportional residual error (fraction)")        # Bogason 2011 Table 2
  })
  model({
    # Individual structural parameters.
    # WBC effect on Vc is centred at the cohort median/mean WBC of 39 x 10^9/L
    # and enters as an additive-fractional multiplier (not power-law).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (1 + e_wbc_vc * (WBC - 39))
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-rate constants for the two-compartment linear model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model. Dose enters the central compartment as an
    # IV infusion (RATE column in the event dataset; dose 60 mg/m^2 over 1 h).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Observation: plasma DNR concentration in mg/L. For comparison to the
    # paper's published nM concentrations, multiply by 1 / MW(DNR) where
    # MW(DNR) = 563.99 g/mol (i.e., 1 mg/L = 1773 nM).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
