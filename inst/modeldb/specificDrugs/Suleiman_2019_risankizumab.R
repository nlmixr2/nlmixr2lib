Suleiman_2019_risankizumab <- function() {
  description <- "Two-compartment population PK model of risankizumab (anti-IL-23 mAb) with first-order SC absorption in healthy subjects and patients with moderate-to-severe plaque psoriasis (Suleiman 2019)"
  reference <- "Suleiman AA, Khatri A, Minocha M, Othman AA. Population Pharmacokinetics of Risankizumab in Healthy Volunteers and Subjects with Moderate to Severe Plaque Psoriasis: Integrated Analyses of Phase I-III Clinical Trials. Clin Pharmacokinet. 2019;58(10):1309-1321. doi:10.1007/s40262-019-00759-z"
  vignette <- "Suleiman_2019_risankizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, Vc, and Vp; normalized as WT/70 per Suleiman 2019 Sect. 2.3 (reference 70 kg was explicitly stated, not median).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALB/44 per Suleiman 2019 Sect. 2.3 (median of all subjects, Table 2). Source uses SI units (g/L); convert g/dL to g/L by x10 if needed.",
      source_name        = "ALB"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as CREAT/76 per Suleiman 2019 Sect. 2.3 (median of all subjects, Table 2).",
      source_name        = "CREAT"
    ),
    CRP = list(
      description        = "Baseline high-sensitivity C-reactive protein (hs-CRP assay; baseline, time-fixed per subject)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as CRP/2.8 per Suleiman 2019 Sect. 2.3 (median of all subjects, Table 2; hs-CRP assay). Source column 'hs-CRP' maps to the canonical general-scope CRP covariate; the assay type (hs-CRP) is documented here rather than via a separate hsCRP canonical.",
      source_name        = "hs-CRP"
    ),
    ADA_TITER = list(
      description        = "Antidrug-antibody reciprocal-dilution titer at the matched PK sample time. Linear-titer convention: 0 = ADA negative.",
      units              = "(reciprocal dilution; 0 = negative)",
      type               = "continuous",
      reference_category = "0 (ADA negative).",
      notes              = "Suleiman 2019 modelled ADA titer as a threshold effect on CL: CL increases by 43% once ADA_TITER >= 128; values below 128 (including 0 / negative) are reference (Sect. 3.2 and Eq. 8 with Titer_threshold = 128). The paper imputes ADA titers reported as < 1 with 0.5 for testing the continuous-power parameterization (Eq. 7), which was not retained; in this library model users may supply 0 for negative samples because only the >=128 threshold matters.",
      source_name        = "ADA titer"
    )
  )

  population <- list(
    n_subjects     = 1899L,
    n_studies      = 7L,
    n_healthy      = 67L,
    n_psoriasis    = 1832L,
    age_range      = "18-85 years (median 47)",
    age_median     = "47 years",
    weight_range   = "42.6-193 kg (median 87)",
    weight_median  = "87 kg",
    sex_female_pct = 29.4,
    race_ethnicity = "White and other 83%, Asian 17% (Table 2). Regions: USA / USA+Canada 62%, Europe 20%, Korea 6%, Japan 5%, Rest of world 4%, Taiwan 2%, China 1%.",
    disease_state  = "Moderate-to-severe plaque psoriasis (six phase I-III studies, n=1832 patients) pooled with healthy male volunteers (one phase I study, n=67).",
    dose_range     = "0.01-5 mg/kg IV, 200-1200 mg IV, 0.25-1 mg/kg SC, 18-300 mg SC (phase I-II); 150 mg SC at weeks 0 and 4 and every 12 weeks thereafter (phase III clinical regimen).",
    regions        = "Global (North America, Europe, East Asia).",
    notes          = "Baseline demographics from Suleiman 2019 Table 2. Data set: 13,123 plasma concentration measurements from 1899 subjects (after excluding BLQ samples and 12 subjects with no post-dose data). Reference covariate values (normalizers in the power-covariate terms): WT = 70 kg (explicitly stated, not median), ALB = 44 g/L (median of all subjects), CREAT = 76 umol/L (median of all subjects), CRP = 2.8 mg/L (median of all subjects, hs-CRP assay). Reference bioavailability is the phase III drug-supply formulation (F = 0.890, logit = 2.09); the phase I-II drug-supply formulation had a distinct F = 0.710 (logit = 0.896) that is not the default here because the approved clinical regimen uses the phase III drug supply. ADA-positive subjects with titer >= 128 were ~1.5% (28/1807) of all phase III ADA-evaluable subjects."
  )

  ini({
    # Structural parameters from Suleiman 2019 Table 3 (final population PK model).
    # Typical values are for the reference subject (70 kg, ALB 44 g/L, CREAT 76 umol/L,
    # CRP 2.8 mg/L, ADA titer < 128, phase III drug supply).
    lcl     <- log(0.243);  label("Clearance (CL, L/day)")                              # Suleiman 2019 Table 3
    lvc     <- log(4.86);   label("Central volume of distribution (Vc, L)")             # Suleiman 2019 Table 3
    lka     <- log(0.229);  label("First-order SC absorption rate (ka, 1/day)")         # Suleiman 2019 Table 3
    lq      <- log(0.656);  label("Intercompartmental clearance (Q, L/day)")            # Suleiman 2019 Table 3
    lvp     <- log(4.25);   label("Peripheral volume of distribution (Vp, L)")          # Suleiman 2019 Table 3

    # Absolute SC bioavailability is parameterized on the logit scale with additive
    # IIV in logit space (Suleiman 2019 Eq. 2 / Sect. 2.3). The phase III
    # formulation estimate (logit = 2.09 -> F = 0.890) is used as the default
    # because the approved clinical regimen uses the phase III drug supply.
    # The phase I-II drug supply had a distinct logit estimate of 0.896
    # (-> F = 0.710, Table 3 footnote c), not carried as a separate parameter here;
    # users wishing to simulate the early-phase formulation can set logitfdepot
    # to 0.896 when constructing the model object.
    logitfdepot <- 2.09; label("Logit of SC bioavailability, phase III drug supply (unitless)")  # Suleiman 2019 Table 3 footnote d

    # Covariate effects on CL (Suleiman 2019 Table 3; power exponents on
    # covariates normalized to their reference values, per Eq. 6).
    e_wt_cl     <-  0.933;   label("Power exponent of body weight on CL (unitless)")        # Suleiman 2019 Table 3
    e_alb_cl    <- -0.715;   label("Power exponent of serum albumin on CL (unitless)")      # Suleiman 2019 Table 3
    e_creat_cl  <- -0.253;   label("Power exponent of serum creatinine on CL (unitless)")   # Suleiman 2019 Table 3
    e_crp_cl    <-  0.044;   label("Power exponent of hs-CRP on CL (unitless)")             # Suleiman 2019 Table 3

    # Covariate effects on volumes (Suleiman 2019 Table 3; power exponents on WT).
    e_wt_vc     <-  1.17;    label("Power exponent of body weight on Vc (unitless)")        # Suleiman 2019 Table 3
    e_wt_vp     <-  0.377;   label("Power exponent of body weight on Vp (unitless)")        # Suleiman 2019 Table 3

    # ADA titer threshold effect on CL (Suleiman 2019 Eq. 8 with Titer_threshold = 128):
    #   ADA_eff = 1              if ADA_TITER <  128
    #   ADA_eff = 1 + e_ada_cl   if ADA_TITER >= 128
    e_ada_cl    <-  0.428;   label("Proportional increase in CL for ADA titer >= 128 (fraction)") # Suleiman 2019 Table 3

    # IIV. Table 3 reports %CV for the log-normal IIV entries with
    #   %CV = SQRT[exp(omega^2) - 1] * 100 (Suleiman 2019 Table 3 footnote e),
    # so we back-transform: omega^2 = log(CV^2 + 1).
    #   CL   24%  -> omega^2 = log(1 + 0.24^2) = 0.05600
    #   Vc   34%  -> omega^2 = log(1 + 0.34^2) = 0.10942
    #   ka   63%  -> omega^2 = log(1 + 0.63^2) = 0.33416
    # The correlation between IIV CL and IIV Vc is 39% (Suleiman 2019 Table 3);
    # covariance = 0.39 * sqrt(0.05600 * 0.10942) = 0.03053.
    etalcl + etalvc ~ c(0.05600,
                        0.03053, 0.10942)                                     # Suleiman 2019 Table 3 (IIV CL 24%, IIV Vc 34%, correlation 39%)
    etalka ~ 0.33416                                                          # Suleiman 2019 Table 3 (IIV ka 63%)

    # IIV on F: additive in logit domain with variance 0.492 (Suleiman 2019
    # Table 3 footnote f, "additive error model in logit domain").
    etalogitfdepot ~ 0.492                                                    # Suleiman 2019 Table 3

    # Residual error: proportional, 19% CV (Suleiman 2019 Table 3).
    # For the prop(propSd) model C_obs = C_pred * (1 + eps) with var(eps) = propSd^2,
    # so propSd = 0.19 corresponds to the 19% CV reported.
    propSd <- 0.19; label("Proportional residual error (SD, fraction)")        # Suleiman 2019 Table 3
  })
  model({
    # Individual PK parameters. Reference subject: 70 kg, ALB 44 g/L,
    # CREAT 76 umol/L, CRP 2.8 mg/L (hs-CRP), ADA titer < 128, phase III drug
    # supply. Covariate forms per Suleiman 2019 Eq. 6 (power models normalized
    # to each covariate's reference value) and Eq. 8 (ADA threshold effect).
    ada_eff <- 1 + e_ada_cl * (ADA_TITER >= 128)

    cl <- exp(lcl + etalcl) *
      (WT    / 70)^e_wt_cl *
      (ALB   / 44)^e_alb_cl *
      (CREAT / 76)^e_creat_cl *
      (CRP   / 2.8)^e_crp_cl *
      ada_eff

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp)          * (WT / 70)^e_wt_vp
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    # Two-compartment model with first-order SC absorption.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Absolute SC bioavailability on the logit scale with additive IIV
    # (Suleiman 2019 Eq. 2). The default logitfdepot uses the phase III
    # drug-supply value (logit 2.09, F = 0.890).
    logit_f  <- logitfdepot + etalogitfdepot
    fdepot   <- exp(logit_f) / (1 + exp(logit_f))
    f(depot) <- fdepot

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
