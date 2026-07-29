Blair_2004_raltitrexed <- function() {
  description <- "Three-compartment population PK model for intravenous raltitrexed (Tomudex) in adult patients with advanced solid tumours, with linear-additive covariate effects of Cockcroft-Gault creatinine clearance on CL and of body weight and serum albumin on central volume (Blair 2004)"
  reference <- "Blair EYL, Rivory LP, Clarke SJ, McLachlan AJ. Population pharmacokinetics of raltitrexed in patients with advanced solid tumours. Br J Clin Pharmacol. 2004;57(4):416-426. doi:10.1111/j.1365-2125.2003.02050.x"
  vignette <- "Blair_2004_raltitrexed"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault formula in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md, following the raw-Cockcroft-Gault pattern of Delattre 2010 amikacin and Aoyama 2012 sepantronium. Population range 23.4-193.0 mL/min, median ~86.5 (Blair 2004 Table 1). Effect form: linear-additive on CL via CL = 0.54 + 0.02 * CRCL (Blair 2004 Table 3); missing CRCL values for 10 patients imputed at the cohort median.",
      source_name        = "CLCR"
    ),
    WT = list(
      description        = "Body weight (baseline, at the start of each treatment course)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Population range 39.0-145.0 kg, median 72.8 (Blair 2004 Table 1). Treatment courses considered as separate individuals per the paper's analysis (Methods: 'Clinical and pharmacokinetic data') with weights recorded per course; treat as baseline within each course. Effect form: linear-additive on V via V = 6.64 + 0.08 * WT - 0.16 * ALB (Blair 2004 Table 3).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin concentration (baseline, at the start of each treatment course)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SI g/L unit (Blair 2004 Table 1 reports albumin in g/L; population range 20.0-47.0, median 36.9). Missing albumin values for 27 patients imputed at the cohort median. Effect form: linear-additive on V via V = 6.64 + 0.08 * WT - 0.16 * ALB (Blair 2004 Table 3); negative coefficient implies higher albumin reduces V, consistent with raltitrexed being >90% protein-bound (paper Discussion).",
      source_name        = "ALB"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 112L,
    n_courses       = 135L,
    n_studies       = 4L,
    age_range       = "21-74 years",
    age_median      = "55.5 years",
    weight_range    = "39.0-145.0 kg",
    weight_median   = "72.8 kg",
    bsa_range       = "1.3-2.4 m^2",
    bsa_median      = "1.8 m^2",
    sex_female_pct  = 45.2,
    race_ethnicity  = "Not reported in source",
    disease_state   = "Advanced solid tumours, predominantly colorectal cancer; also breast and ovarian cancers",
    dose_range      = "0.1-4.5 mg/m^2 raltitrexed as a 15-30 min IV infusion every 3 weeks; absolute doses 0.2-9.2 mg (Blair 2004 Table 1)",
    regions         = "Pooled across one European Phase I dose-finding study, one US Phase I dose-finding study, a radiolabel mass-balance disposition study, and an open-label renal-function study",
    renal_function  = "Cockcroft-Gault creatinine clearance 23.4-193.0 mL/min (median ~86.5); raw mL/min, NOT BSA-normalized",
    hepatic_function = "Baseline ALT 1-50 U/L (median ~14.7); AST 5-86 U/L (median ~20.5); total bilirubin 2-23 umol/L (median ~10.1); no significant correlation between liver function tests and raltitrexed clearance found in the analysis",
    notes           = "Total cohort: 112 adult patients, 135 treatment courses, 2105 raltitrexed plasma concentration observations (3-23 per patient). 23 of 112 patients received a second course; each course was treated as a separate individual in the P-Pharm popPK analysis (Blair 2004 Methods). Data-splitting validation used a 2:1 random allocation into 90-patient (1374 observations) model-development and 45-patient (731 observations) model-validation datasets, then the final model was refit on the pooled total cohort (135 courses) for the parameter estimates reported here. Software: P-Pharm 1.5.1 (InnaPhase). Bioanalysis: radioimmunoassay with sheep antiserum, intra-assay CV 10.7%, inter-assay CV 12.9%, LLOQ 0.2 ng/mL. The clinical studies are referenced as Clarke 2000 (European Phase I, ref 6), Grem 1999 (US Phase I, ref 7), Beale 1998 (radiolabel mass-balance, ref 8) and Judson 1998 (renal-function study, ref 9)."
  )

  ini({
    # Structural PK parameters - reference values from Blair 2004 Table 3
    # (final model, total cohort with covariates) for the CL and V regression
    # intercepts, and from Table 2 (covariate model, model development dataset)
    # for the inter-compartmental first-order rate constants k12, k21, k13, k31.
    # Parameterisation follows Blair 2004's reported model: CL and V depend
    # linearly on covariates, the four rate constants are estimated directly
    # (analogous to the K-rate parameterisation in Ferron_2013_cabazitaxel).
    lcl  <- log(0.54);  label("Intercept of CL regression at CRCL = 0 (L/h)") # Blair 2004 Table 3 q1 = 0.54 +/- 0.12; CL = q1 + q2 * CRCL
    lvc  <- log(6.64);  label("Intercept of Vc regression at WT = 0, ALB = 0 (L)") # Blair 2004 Table 3 q3 = 6.64 +/- 1.26; V = q3 + q4 * WT + q5 * ALB
    lk12 <- log(0.99);  label("Rate constant central -> peripheral1 (k12, 1/h)") # Blair 2004 Table 2 covariate model
    lk21 <- log(0.97);  label("Rate constant peripheral1 -> central (k21, 1/h)") # Blair 2004 Table 2 covariate model
    lk13 <- log(0.96);  label("Rate constant central -> peripheral2 (k13, 1/h)") # Blair 2004 Table 2 covariate model
    lk31 <- log(0.01);  label("Rate constant peripheral2 -> central (k31, 1/h)") # Blair 2004 Table 2 covariate model

    # Linear covariate effects (paper's regression-coefficient parameterisation).
    # Coefficients have units of (parameter unit) per (covariate unit), and are
    # added in linear (un-transformed) space inside model() to reproduce the
    # paper's regression equations CL = q1 + q2 * CRCL and
    # V = q3 + q4 * WT + q5 * ALB.
    e_crcl_cl <- 0.02;  label("CRCL slope on CL (L/h per mL/min)")             # Blair 2004 Table 3 q2 = 0.02 +/- 0.003
    e_wt_vc   <- 0.08;  label("WT slope on Vc (L per kg)")                     # Blair 2004 Table 3 q4 = 0.08 +/- 0.02
    e_alb_vc  <- -0.16; label("ALB slope on Vc (L per g/L)")                   # Blair 2004 Table 3 q5 = -0.16 +/- 0.03

    # Inter-individual variability. Blair 2004 reports IIV as %CV on a
    # normal-distribution scale (P-Pharm software; paper Methods: "normal
    # distribution of interpatient variability in the clearance, volume of
    # distribution and the distributional first-order rate constants").
    # nlmixr2lib encodes IIV as log-normal (canonical convention); convert via
    # omega^2 = log(1 + CV^2), the small-CV approximation common across the
    # registry (e.g., Ferron_2013_cabazitaxel). The translation is exact for
    # small-to-moderate CV (<30%) and the resulting individual-parameter
    # distributions are very close to the paper's normal-IIV simulation.
    etalcl  ~ 0.0755 # Blair 2004 Table 3 CL %CV = 28; log(1 + 0.28^2) = 0.0755
    etalvc  ~ 0.0606 # Blair 2004 Table 3 V %CV = 25; log(1 + 0.25^2) = 0.0606
    etalk12 ~ 0.0384 # Blair 2004 Table 2 covariate model k12 %CV = 19.78
    etalk21 ~ 0.0188 # Blair 2004 Table 2 covariate model k21 %CV = 13.77
    etalk13 ~ 0.0786 # Blair 2004 Table 2 covariate model k13 %CV = 28.6
    etalk31 ~ 0.0750 # Blair 2004 Table 2 covariate model k31 %CV = 27.9

    # Residual error. The paper fits "heteroscedastic residual error (with a
    # weighing factor of 1/concentration^2)" with sigma = 0.062; the 1/Cp^2
    # weighting and the small numerical value identify this as a proportional
    # error structure with SD ~6.2% on the linear concentration scale (Blair
    # 2004 Methods and Table 2; Tables 2 and 3 both report Sigma = 0.062).
    propSd <- 0.062; label("Proportional residual error (fraction)") # Blair 2004 Tables 2-3 sigma = 0.062
  })

  model({
    # Individual PK parameters. CL and Vc apply the paper's linear-additive
    # regression equations on the linear scale; the lognormal-IIV multiplicative
    # term exp(etalcl) / exp(etalvc) represents the IIV translated to the
    # nlmixr2lib canonical form (the paper used additive normal IIV - see
    # ini() comments).
    cl  <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)
    vc  <- (exp(lvc) + e_wt_vc * WT + e_alb_vc * ALB) * exp(etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    k13 <- exp(lk13 + etalk13)
    k31 <- exp(lk31 + etalk31)

    kel <- cl / vc

    # Three-compartment IV model. No depot (raltitrexed was administered as a
    # 15-30 min IV infusion into central via the dose rate column).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Concentration. Dose in mg, vc in L -> mg/L = ug/mL = 1000 ng/mL. Multiply
    # by 1000 to express the central concentration in ng/mL (= ug/L), the unit
    # used by the paper (LLOQ 0.2 ng/mL; prediction errors reported in ug/L).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
