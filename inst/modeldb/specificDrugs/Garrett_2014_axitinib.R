Garrett_2014_axitinib <- function() {
  description <- "Two-compartment population PK model for axitinib in healthy volunteers (Garrett 2014). First-order absorption with fixed lag time, allometric power-form effect of body weight on the central volume of distribution (reference 75 kg), linear-proportional fasting effects on the first-order absorption rate constant ka and on bioavailability F, and a linear-proportional reduction in F for the marketed crystal polymorph Form XLI relative to the earlier Form IV reference. Pooled data from 337 healthy subjects across ten Pfizer Phase I studies."
  reference   <- "Garrett M, Poland B, Brennan M, Hee B, Pithavala YK, Amantea MA. Population pharmacokinetic analysis of axitinib in healthy volunteers. Br J Clin Pharmacol. 2014;77(3):480-492. doi:10.1111/bcp.12206"
  vignette    <- "Garrett_2014_axitinib"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form allometric effect on Vc with reference 75 kg (population median); the exponent 0.758 was freely estimated rather than fixed at 1. Garrett 2014 Results 'Full model and final model' and Table 3.",
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-vs-fasted indicator at the dose record.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (fed; the typical-value reference for ka and F in this model).",
      notes              = "The fasting effect is applied internally as (1 - FED), so FED = 1 fed leaves ka and F at their typical-value Form IV / fed estimates and FED = 0 fasted activates the linear-proportional fasting increases (Garrett 2014 Results 'Full model and final model' and Table 3). Per-record covariate.",
      source_name        = "FED"
    ),
    FORM_AXI_XLI = list(
      description        = "Crystal polymorph form indicator for axitinib (Form XLI marketed vs Form IV earlier).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Form IV; the typical-value F reference in Garrett 2014 Table 3).",
      notes              = "1 = Form XLI (marketed commercial crystal polymorph); 0 = Form IV (earlier Phase I crystal polymorph and the typical-value F reference). Linear-proportional effect on F only (no effect on ka or any other parameter retained in the final model). Garrett 2014 Results 'Full model and final model' and Table 3.",
      source_name        = "FORM"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 337L,
    n_studies       = 10L,
    age_range       = "18+ years (healthy volunteers)",
    age_median      = "31.0 years (mean 34.1 +/- 11.6)",
    weight_range    = "approximately 50-110 kg (median 75.0; mean 76.7 +/- 11.6)",
    weight_median   = "75.0 kg",
    sex_female_pct  = 7,
    race_ethnicity  = c(White = 62, Black = 8, Asian = 18, Japanese = 6, Hispanic = 3, Other = 3),
    disease_state   = "Healthy volunteers (no axitinib indication; pooled Phase I clinical pharmacology studies).",
    dose_range      = "Single 5 mg oral dose (Form IV or Form XLI) in fed or fasted state; one study additionally administered a 1 mg intravenous dose of Form IV (n = 16) to estimate absolute bioavailability.",
    regions         = "United States (Austin TX, Fargo ND, La Mesa CA, Plantation FL), Belgium (Bruxelles), and Singapore.",
    n_observations  = "Pooled from ten Phase I clinical studies (Table 1). Study designs included single-dose absolute-bioavailability, food-effect, formulation-comparison (Form IV vs Form XLI), and drug-drug-interaction arms; analyses included only the single-dose-axitinib alone records, excluding doses administered with rifampicin or ketoconazole.",
    smoking_status  = c(`Non-smoker` = 88, `Ex-smoker` = 12),
    notes           = "Demographic counts reproduced from Garrett 2014 Table 2 (n = 337 unique subjects across ten studies). Genetic covariates UGT1A1*28 and CYP2C19 inferred phenotype were tested but not retained as significant in the final model; ALT, AST, bilirubin, creatinine clearance, age, sex, race, and smoking status were likewise tested and dropped (Garrett 2014 Results 'Full model and final model')."
  )

  ini({
    # Structural PK parameters -- Garrett 2014 Table 3 final-model estimates.
    # Reference subject: a 75 kg healthy volunteer receiving 5 mg axitinib
    # oral (Form IV) in the fed state. Apparent oral clearance CL/F and
    # apparent oral volumes are anchored by the 1 mg intravenous arm
    # (study 3) so the canonical CL and Vc are reported without the /F
    # abuse-of-notation.
    lka     <- log(0.523)        ; label("Absorption rate constant ka, Form IV fed reference (1/h)")        # Garrett 2014 Table 3: ka (h-1) Form IV, Fed = 0.523
    lcl     <- log(17.0)         ; label("Systemic clearance CL (L/h)")                                     # Garrett 2014 Table 3: CL = 17.0 L/h
    lvc     <- log(45.3)         ; label("Central volume of distribution Vc at WT = 75 kg (L)")             # Garrett 2014 Table 3: Vc = 45.3 L (at WT median 75 kg)
    lq      <- log(1.74)         ; label("Inter-compartmental clearance Q (L/h)")                           # Garrett 2014 Table 3: Q = 1.74 L/h
    lvp     <- log(45.9)         ; label("Peripheral volume of distribution Vp (L)")                        # Garrett 2014 Table 3: Vp = 45.9 L
    lfdepot <- log(0.465)        ; label("Absolute oral bioavailability F, Form IV fed reference (unitless)")  # Garrett 2014 Table 3: F, Form IV, Fed = 0.465
    ltlag   <- fixed(log(0.457)) ; label("Absorption lag time (h; FIXED)")                                  # Garrett 2014 Table 3 footnote: tlag fixed at 0.457 in the final run

    # Covariate effects -- Garrett 2014 Table 3 final-model estimates.
    # Weight on Vc is an estimated allometric exponent (75 kg reference).
    # Fasting and Form XLI effects are linear-proportional shifts in the
    # paper's notation theta1 * (1 + theta2 * indicator).
    e_wt_vc      <- 0.758    ; label("Allometric exponent of (WT / 75 kg) on Vc (unitless)")                            # Garrett 2014 Table 3: Weight effect on Vc = 0.758
    e_fast_ka    <- 2.07     ; label("Linear-proportional fasting effect on ka (relative to Form IV fed reference)")    # Garrett 2014 Table 3: ka Fasting = 2.07 (i.e. ka_fasted = ka_fed * (1 + 2.07))
    e_fast_f     <- 0.338    ; label("Linear-proportional fasting effect on F (relative to Form IV fed reference)")     # Garrett 2014 Table 3: F Fasting = 0.338 (i.e. F_fasted = F_fed * (1 + 0.338))
    e_xli_f      <- -0.150   ; label("Linear-proportional Form XLI effect on F (relative to Form IV reference)")        # Garrett 2014 Table 3: F Form XLI = -0.150 (i.e. F_XLI = F_IV * (1 - 0.150))

    # IIV -- Garrett 2014 Table 3 final-model omega^2 values. The paper
    # parameterises exponential IIV on each structural PK parameter:
    # %IIV = sqrt(omega^2) * 100. Block correlations between CL and Vc
    # and between Q and Vp were retained in the final model.
    etalcl + etalvc ~ c(0.272,
                        0.141, 0.0949)   # Garrett 2014 Table 3: omega^2 CL = 0.272, cov(CL,Vc) = 0.141, omega^2 Vc = 0.0949
    etalq  + etalvp ~ c(0.406,
                        0.619, 1.07)     # Garrett 2014 Table 3: omega^2 Q = 0.406, cov(Q,Vp) = 0.619, omega^2 Vp = 1.07
    etalka          ~ 0.506              # Garrett 2014 Table 3: omega^2 ka = 0.506

    # Residual error -- Garrett 2014 Methods 'Model development' fit the
    # residual additively in a logarithmic scale (NONMEM Y = LOG(IPRED) +
    # EPS(1)), which corresponds to a log-normal residual on the linear
    # concentration scale. Table 3 reports the residual SD on the log
    # scale as a percentage (oral 50.9%, intravenous 34.2%); only the
    # oral residual is encoded here because the model is dominantly used
    # for simulating oral dosing.
    expSd <- 0.509 ; label("Log-normal residual error SD (oral; log-scale fraction)")  # Garrett 2014 Table 3 final: oral residual = 50.9%
  })

  model({
    # Reference covariate values (Garrett 2014 Methods and Table 3 footnote).
    ref_wt <- 75

    # Individual PK parameters. Weight enters allometrically on Vc only;
    # CL, Q, and Vp have no retained covariates in the final model. The
    # fasting effect on ka is the linear-proportional theta1 * (1 + theta2 * (1-FED))
    # form from the Methods 'Statistical analysis' section.
    ka <- exp(lka + etalka) * (1 + e_fast_ka * (1 - FED))
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Absolute oral bioavailability with the two retained covariate
    # effects: fasting and crystal polymorph Form XLI. F1 = 1 is anchored
    # by the intravenous arm (study 3).
    fdepot <- exp(lfdepot) * (1 + e_fast_f * (1 - FED)) * (1 + e_xli_f * FORM_AXI_XLI)

    tlag <- exp(ltlag)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption and lag time.
    # NONMEM ADVAN4/TRANS4 equivalent. Concentrations are expressed in
    # ng/mL by scaling the central compartment amount-per-litre output by
    # 1000 (axitinib dose entered in mg, internal mass/volume in mg/L).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot)    <- fdepot

    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
