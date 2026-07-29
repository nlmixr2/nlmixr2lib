Royer_2010_HuHMFG1 <- function() {
  description <- "Two-compartment population PK model with linear elimination for HuHMFG1 (AS1402), a humanised anti-MUC1 monoclonal antibody, in patients with metastatic breast cancer; serum AST enters the typical clearance equation additively (Royer 2010)"
  reference <- "Royer B, Yin W, Pegram M, Ibrahim N, Villanueva C, Mir D, Erlandsson F, Pivot X. Population pharmacokinetics of the humanised monoclonal antibody, HuHMFG1 (AS1402), derived from a phase I study on breast cancer. Br J Cancer. 2010;102(5):827-832. doi:10.1038/sj.bjc.6605560"
  vignette <- "Royer_2010_HuHMFG1"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the typical-CL equation additively as 0.0036 * (AST / ASTmed), with ASTmed = 29 U/L (study population median from Royer 2010 Table 2). The Royer 2010 covariate parameterisation is additive (intercept 0.016 L/h + 0.0036 L/h * (AST / 29)) rather than the more common power form (AST / ref)^exponent. Source column 'AST' maps to the canonical AST covariate.",
      source_name        = "AST"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    n_observations = 435L,
    age_range      = "32-72 years (median 55)",
    weight_range   = "50-108 kg (median 73.1)",
    height_range   = "155-179 cm (median 165)",
    bmi_range      = "18.00-40.65 kg/m^2 (median 25.75)",
    sex_female_pct = NA_real_,
    disease_state  = "Locally advanced or metastatic breast cancer, previously treated with up to three chemotherapeutic regimens (including neoadjuvant and adjuvant therapy)",
    dose_range     = "1, 3, 9, or 16 mg/kg IV infusion; infusion duration depended on the dose level (60 min at 1-3 mg/kg, 120 min at 9 mg/kg, 180 min at 16 mg/kg). Multiple administrations were given (up to 10 in one patient) per the study schedule.",
    regions        = "Multi-centre (France and United States; Antisoma-sponsored phase I study reported by Pegram 2009)",
    baseline_labs  = "Median (range) per Royer 2010 Table 2: total protein 68 g/L (54-82), albumin 41 g/L (27-48), creatinine clearance 85.3 mL/min (37.2-648.9), alkaline phosphatase 86 U/L (47-1127), ALT 27 U/L (6-252), AST 29 U/L (14-1392), GGT 37 U/L (9-1386), CA15-3 45 (9-16400), CA27.29 31.3 (3.5-1723).",
    notes          = "Patient population described in Royer 2010 Table 2 (n = 26; 3, 9, 6, and 8 patients in the 1, 3, 9, and 16 mg/kg dose cohorts respectively). Sex distribution is not explicitly tabulated in Royer 2010; the parent phase I study (Pegram 2009) enrolled women with breast cancer. Human anti-human antibody (HAHA) testing was negative in every assessable patient (HAHA detected in 0 of 23) and was therefore not investigated further. Data were analysed with NONMEM VI level 1.0 (ADVAN3 TRANS4, FOCE with INTERACTION)."
  )

  ini({
    # Structural parameters (Royer 2010 Table 3). The covariate equation for CL
    # is additive (Royer 2010 Results, page 829):
    #   CL (L/h) = 0.016 + 0.0036 * (AST / ASTmed),  ASTmed = 29 U/L (Table 2).
    # The 0.016 L/h intercept and the 0.0036 L/h AST coefficient are reported
    # separately in Table 3 (as theta_CL and COV_CL). Volumes and Q are reported
    # directly with no covariate adjustment.
    lcl <- log(0.016); label("Intercept component of the typical clearance (L/h); AST-independent term in the additive CL equation")  # Royer 2010 Table 3 (theta_CL = 0.016 L/h)
    lvc <- log(3.31);  label("Central volume of distribution V1 (L)")                                                                 # Royer 2010 Table 3 (theta_V1 = 3.31 L)
    lq  <- log(0.017); label("Inter-compartmental clearance Q (L/h)")                                                                 # Royer 2010 Table 3 (theta_Q = 0.017 L/h)
    lvp <- log(2.33);  label("Peripheral volume of distribution V2 (L)")                                                              # Royer 2010 Table 3 (theta_V2 = 2.33 L)

    # AST additive coefficient on CL (Royer 2010 Table 3 COV_CL row). Linear-scale,
    # may in principle be negative -- left untransformed, no log() prefix.
    e_ast_cl <- 0.0036; label("Additive AST coefficient on CL (L/h per unit AST/ASTmed)")  # Royer 2010 Table 3 (COV_CL = 0.0036)

    # Inter-individual variability (Royer 2010 Table 3). The paper reports
    # omega^2 (variance) on a log-normal exponential random effect for V1 and CL
    # only; random effects on Q and V2 could not be obtained and inter-occasion
    # variability was insignificant (no IOV stored). omega^2 values transcribed
    # directly as the variances of etalvc and etalcl.
    etalvc ~ 0.042  # Royer 2010 Table 3 (omega^2_V1 = 0.042)
    etalcl ~ 0.062  # Royer 2010 Table 3 (omega^2_CL = 0.062)

    # Residual error (Royer 2010 Table 3). The paper reports sigma^2 (variances)
    # for a combined proportional + additive error model; the additive sigma^2
    # was FIXED at 2.26. nlmixr2 stores residual error as SDs, so propSd and
    # addSd are sqrt() of the reported variances.
    propSd <- sqrt(0.034);        label("Proportional residual error (SD, fraction)")  # Royer 2010 Table 3 (sigma^2_prop = 0.034)
    addSd  <- fixed(sqrt(2.26));  label("Additive residual error (SD, mg/L)")          # Royer 2010 Table 3 (sigma^2_add = 2.26, FIXED)
  })
  model({
    # Typical clearance follows the Royer 2010 Results page-829 additive equation
    #     CL = 0.016 + 0.0036 * (AST / 29)
    # with the inter-individual variability applied multiplicatively as an
    # exponential random effect (paper Methods, page 828, "exponential random
    # effect"). The constant term is recovered as exp(lcl) so the parameter is
    # stored on the canonical log-positive scale; e_ast_cl is the additive
    # AST coefficient on the linear scale.
    tvcl <- exp(lcl) + e_ast_cl * (AST / 29)
    cl   <- tvcl * exp(etalcl)

    # Volumes and inter-compartmental clearance: no covariate effects in the
    # final Royer 2010 model. V1 carries IIV (exp(lvc + etalvc)); Q and V2 are
    # reported without estimable random effects.
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for the ADVAN3 TRANS4 (linear two-compartment) structure.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Dosing is intravenous (zero-order infusion via the event
    # table's rate / amt columns); no depot or absorption compartment is needed.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration in central compartment: dose in mg, volume in L -> mg/L.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
