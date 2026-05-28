Golubovic_2019_sirolimus <- function() {
  description <- "Two-compartment population PK model for sirolimus in adult kidney transplant recipients on triple immunosuppressive therapy (sirolimus + mycophenolate mofetil + corticosteroids) developed from routine therapeutic-drug-monitoring trough data with the NONMEM informative-prior functionality (Golubovic 2019). Covariate effects on CL/F: aspartate aminotransferase greater than 37 IU/L as a binary indicator of elevated liver enzymes (-37 percent multiplicative effect via power form 0.63^AST_HIGH) and age as a linear-deviation effect on CL/F with reference age 44 years (coefficient -0.388 on AGE/44, reproducing the 49 percent CL/F decrease from age 16 to age 64 reported in the Discussion)."
  reference <- paste(
    "Golubovic B, Vucicevic K, Radivojevic D, Vezmar Kovacevic S,",
    "Prostran M, Miljkovic B. Exploring sirolimus pharmacokinetic",
    "variability using data available from the routine clinical care",
    "of renal transplant patients - population pharmacokinetic approach.",
    "J Med Biochem. 2019;38(3):323-331. doi:10.2478/jomb-2018-0030.",
    sep = " "
  )
  vignette <- "Golubovic_2019_sirolimus"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age (baseline; years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on CL/F with reference age 44 years (close to the cohort mean of 43.22 years and present in the paper's final equation): cl_eff = (1 + e_age_cl * AGE / 44). Coefficient e_age_cl = -0.388 reproduces the 49 percent CL/F decrease from age 16 to age 64 reported in the Discussion.",
      source_name        = "AGE"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase activity (baseline or time-varying).",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as ast_high <- (AST > 37) per Golubovic 2019 Results: AST entered as a 0/1 categorical indicator with 1 = AST greater than 37 IU/L (the laboratory upper limit of normal used in the source paper). Power-form effect on CL/F: 0.630^ast_high (-37 percent when ast_high = 1).",
      source_name        = "AST"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 25L,
    n_studies        = 1L,
    n_observations   = 250L,
    age_range        = "16-64 years",
    age_median       = "43.22 years (mean)",
    weight_range     = "44-128 kg",
    weight_median    = "77.07 kg (mean)",
    sex_female_pct   = 28.0,
    race_ethnicity   = "Not reported (single-country Serbian cohort).",
    disease_state    = "Adult kidney transplant recipients converted to sirolimus from a calcineurin inhibitor (tacrolimus or cyclosporine) as the second-line immunosuppressive treatment; all on triple immunosuppression with mycophenolate mofetil and corticosteroids.",
    dose_range       = "0.5-15 mg/day oral sirolimus titrated to maintain trough blood concentrations of 8-20 ng/mL.",
    regions          = "Serbia (single centre: Nephrology Clinic, Clinical Center of Serbia, University of Belgrade).",
    co_medication    = "Mycophenolate mofetil (mean 1104 mg/day, range 0-2000) and corticosteroids (mean 10.74 mg/day, range 0-50).",
    graft_origin     = "Living donor n = 23 (92 percent); cadaveric n = 2 (8 percent).",
    dialysis_pre     = "Pre-transplant dialysis n = 21 (84 percent); no pre-transplant dialysis n = 4 (16 percent).",
    assay            = "Architect Sirolimus chemiluminescent microparticle immunoassay (Abbott Laboratories), nominal measurement range 2-30 ng/mL; samples above range diluted per manufacturer protocol.",
    notes            = "Single-centre retrospective TDM cohort. Demographics summarised here are Table I model-development column (n = 25); an independent external-validation cohort of 13 newly converted patients is described in the source but is not encoded in this model. All samples are end-of-dosing-interval trough concentrations drawn before the morning dose. Sex: 18 male (72 percent) / 7 female (28 percent)."
  )

  ini({
    # Final-model THETAs from Golubovic 2019 Table IV (estimates column).

    # Apparent clearance for the reference patient (AST <= 37 IU/L, AGE -> 0).
    # Note: this is the THETA in the equation
    #   CL/F = 12.2 * 0.63^AST_HIGH * (1 - 0.388 * AGE / 44)
    # so the linear-deviation age term is centred such that THETA equals
    # CL/F at AGE = 0 rather than at the cohort mean. See vignette
    # "Assumptions and deviations" for the typical-value calculation at
    # cohort-mean age.
    lcl <- log(12.2); label("Apparent clearance, CL/F, at AGE = 0 and AST <= 37 IU/L (L/h)")  # Table IV: CL/F = 12.2 L/h
    lvc <- log(118);  label("Apparent central volume of distribution, Vc/F (L)")              # Table IV: Vc/F = 118 L
    lvp <- log(609);  label("Apparent peripheral volume of distribution, Vp/F (L)")           # Table IV: Vp/F = 609 L
    lq  <- log(5.07); label("Apparent inter-compartmental clearance, Q/F (L/h)")              # Table IV: Q/F = 5.07 L/h

    # ka and its IIV were held at the literature prior values (Jiao 2009 and
    # Dansirikul 2005 prior bundle): the Methods state "we used the same
    # value for this parameter and its interindividual variability as in
    # 2-COMP" and the Discussion notes "the change in value of ka was not
    # observed." Table IV reports ka with an SE of 4.79e-5 (relative SE
    # approximately 2e-5, i.e., numerical-noise level) and bootstrap
    # confidence interval 2.19-2.19, consistent with fix semantics rather
    # than true posterior uncertainty.
    lka <- fixed(log(2.19)); label("First-order absorption rate constant, ka (1/h; prior-fixed)")  # Table IV: ka = 2.19 1/h (prior-fixed)

    # Covariate effects on CL/F. The paper's final equation is
    #   CL/F = 12.2 * 0.63^AST_HIGH * (1 - 0.388 * AGE / 44)
    # with AST_HIGH = 1 when AST > 37 IU/L (0 otherwise) and AGE in years.
    # We store the AST coefficient on the log scale so the model() body can
    # write exp(e_ast_cl * ast_high), reproducing 0.63^ast_high; the AGE
    # coefficient is a bare linear factor used directly in (1 + e_age_cl *
    # AGE / 44).
    e_ast_cl <- log(0.630); label("Log-effect of elevated AST (> 37 IU/L) on CL/F (unitless)")           # Table IV: theta_AST = 0.630 (-37 percent on CL/F)
    e_age_cl <- -0.388;     label("Linear-deviation coefficient on AGE / 44 for CL/F (unitless)")        # Table IV: theta_AGE = -0.388 (-49 percent across 16 -> 64 years)

    # Inter-individual variability. The Methods state "Interindividual
    # variability was evaluated by an exponential model", so each reported
    # omega^2 is the variance on the natural-log scale and translates
    # directly to nlmixr2 syntax `~ omega2`. ka and its IIV were prior-fixed
    # (see lka comment); the remaining etas are estimated with informative
    # priors from Dansirikul 2005 / Jiao 2009.
    etalcl ~ 0.0547           # Table IV: omega^2 CL/F = 0.0547
    etalvc ~ 0.306            # Table IV: omega^2 Vc/F = 0.306
    etalvp ~ 0.0657           # Table IV: omega^2 Vp/F = 0.0657
    etalq  ~ 0.103            # Table IV: omega^2 Q/F  = 0.103
    etalka ~ fixed(0.145)     # Table IV: omega^2 ka  = 0.145 (prior-fixed; SE 1.28e-6, bootstrap CI 0.144-0.144)

    # Residual unexplained variability: "slope-intercept" (combined additive
    # + proportional) error model (Results, "the residual variability was
    # best described with slope-intercept error model"; Table IV Wa, Wp).
    addSd  <- 1.93;  label("Additive residual error, Wa (ng/mL)")                # Table IV: Wa = 1.93 ng/mL
    propSd <- 0.249; label("Proportional residual error, Wp (fraction)")         # Table IV: Wp = 0.249
  })

  model({
    # Reference values from Golubovic 2019 final equation (Results, p. 327).
    ref_age <- 44   # reference age in years (cohort mean 43.22 rounded; reproduces the published equation)
    ast_uln <- 37   # upper limit of normal for AST (IU/L) per Results

    # Binary AST indicator derived from continuous AST: 0 if AST <= 37 IU/L,
    # 1 if AST > 37 IU/L (Golubovic 2019 Results).
    ast_high <- (AST > ast_uln)

    # Individual structural parameters. CL/F carries the AST power-form
    # multiplier and the AGE linear-deviation factor; the other parameters
    # have no covariate effects in the final model.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * exp(e_ast_cl * ast_high) * (1 + e_age_cl * AGE / ref_age)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Micro-rate constants for the 2-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: oral 2-compartment model with first-order absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Observation: whole-blood sirolimus concentration in ng/mL. Dose is in
    # mg and volumes in L, so central / vc gives mg/L; multiply by 1000 to
    # convert to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
