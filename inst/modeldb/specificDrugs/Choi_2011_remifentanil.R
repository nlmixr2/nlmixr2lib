Choi_2011_remifentanil <- function() {
  description <- paste(
    "Combined effect-and-tolerance pharmacodynamic model for the central-nervous-system effect of",
    "remifentanil on EEG-derived temporal linear mode complexity (TLMC) in healthy volunteers",
    "(Choi 2011). The PD model captures both depression of CNS activity during infusion and",
    "rebound during recovery via a sigmoid Emax driver from an effect compartment plus an",
    "opposing sigmoid term driven by a slow tolerance compartment. The TLMC observation is",
    "baseline-normalised so the readout is dimensionless with baseline E0 = 1 by construction.",
    "The 3-compartment IV remifentanil PK underneath is fixed from the upstream Kang 2007 BJCP",
    "popPK (Kang 2007 Table 2; rate-constant parameterisation; AGE and BSA additive covariate",
    "effects on the elimination rate constant k10). Choi 2011 itself did not re-estimate PK -- it",
    "read per-subject individual PK estimates from the Kang 2007 fit as input data columns.",
    "Reported model-selection metrics across the three competing PD models in Choi 2011 (combined",
    "effect-and-tolerance, feedback, sigmoid Emax) identified the combined model implemented here",
    "as the best by AIC (-6966) and positive predictive value of rebound (100%); the feedback and",
    "sigmoid Emax models did not capture rebound (Table 3).",
    sep = " "
  )
  reference <- paste(
    "Choi BM, Shin DH, Noh MH, Kim YH, Jeong YB, Lee SH, Lee EK, Noh GJ.",
    "Temporal linear mode complexity as a surrogate measure of the effect of remifentanil",
    "on the central nervous system in healthy volunteers.",
    "Br J Clin Pharmacol 2011; 71(6):879-888.",
    "doi:10.1111/j.1365-2125.2011.03904.x.",
    "PK structure and parameter values adapted from",
    "Kang SH, Poynton MR, Kim KM, Lee H, Kim DH, Lee SH, Bae KS, Linares O, Kern SE, Noh GJ.",
    "Population pharmacokinetic and pharmacodynamic models of remifentanil in healthy",
    "volunteers using artificial neural network analysis.",
    "Br J Clin Pharmacol 2007; 64(1):3-13.",
    "doi:10.1111/j.1365-2125.2007.02845.x.",
    sep = " "
  )
  vignette <- "Choi_2011_remifentanil"
  paper_specific_compartments <- c("tolerance", "tlmc")
  units <- list(
    time = "min",
    dosing = "ug",
    concentration = "ng/mL"
  )
  # k10/k12/k21/k13/k31 are in 1/min per Kang 2007 Table 2; the natural
  # time unit is the minute. Dose in micrograms (remifentanil is typically
  # infused at 1-8 ug/kg/min, so a 60 kg subject at 3 ug/kg/min = 180 ug/min)
  # divided by V1 in L gives ug/L = ng/mL, matching Kang 2007 Table 2's
  # concentration units (assay LLOQ 0.05 ng/mL, see Kang 2007 Methods).

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Enters the structural k10 elimination rate from Kang 2007 Table 2 via the additive form k10 = theta2 - theta9 * (AGE / 100) + theta10 * BSA (Kang 2007 Methods: 'Pharmacokinetic analysis', Table 2). Choi 2011 cohort spanned ages 20-79 years (Choi 2011 Table 1; n = 28 healthy volunteers).",
      source_name        = "AGE"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Enters the structural k10 elimination rate from Kang 2007 Table 2 via the additive form k10 = theta2 - theta9 * (AGE / 100) + theta10 * BSA. Choi 2011 Table 1 does not report BSA directly; downstream users compute BSA from height and weight (Du Bois, Mosteller, or Gehan-George; Kang 2007 Methods does not state the BSA formula either). Approximate cohort BSA derived from Choi 2011 Table 1 demographics (overall mean WT approx 60.6 kg, mean HT approx 162.4 cm) is BSA approx 1.65 m^2 via the Du Bois formula.",
      source_name        = "BSA"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 28L,
    n_studies        = 1L,
    age_range        = "20-79 years",
    age_subgroups    = "Young (<= 40 years, n = 9, M:F = 5:4, mean 27.8 +/- 5.7 y); middle-aged or elderly (> 40 years, n = 19, M:F = 9:10, mean 57.7 +/- 13.9 y).",
    weight_range     = "Young mean 63.8 +/- 11.3 kg; middle/elderly mean 59.1 +/- 10.6 kg (Choi 2011 Table 1)",
    height_range     = "Young mean 167.9 +/- 7.5 cm; middle/elderly mean 160.8 +/- 9.6 cm",
    sex_female_pct   = 50,
    disease_state    = "Healthy volunteers with no medical problems or abnormal laboratory test results.",
    dose_range       = "Zero-order remifentanil infusion. Young volunteers were randomized to a fixed rate of 1, 2, 3, 4, 5, 6, 7, or 8 ug/kg/min for 15-20 min (total dose 5.4 +/- 2.8 mg). Middle/elderly volunteers received 3 ug/kg/min until 95% spectral edge frequency stopped changing (mean infusion duration 15.3 +/- 3.9 min; total dose 2.7 +/- 0.8 mg).",
    regions          = "South Korea (Asan Medical Center, Seoul; healthy volunteers under IRB approval).",
    notes            = "Choi 2011 reused the blood-concentration and raw-EEG data from Noh 2006 Anesthesiology 104:921-32 [Choi 2011 ref 2]. The PK was not re-estimated in Choi 2011; per-subject individual PK estimates from the upstream Kang 2007 BJCP 64:3-13 popPK fit [Choi 2011 ref 3] were used as input data columns IV1, IK10, IK12, IK13, IK21, IK31 in the NONMEM data files (Choi 2011 Appendix 3 control streams). EEG was recorded on parietal montage P4 (right-handed cohort); TLMC and ApEn were derived from 10 s epochs and each normalised to per-subject baseline (Choi 2011 Methods 'Electroencephalographic analysis'), so the PD effect is dimensionless with baseline E0 = 1 by construction. The total observation window was up to 240 min after discontinuation of the infusion."
  )

  ini({
    # ====================================================================
    # PK PARAMETERS (FIXED FROM KANG 2007 BJCP TABLE 2 'Estimate' column).
    # Kang 2007 reports the 3-compartment IV remifentanil popPK fit on
    # healthy volunteers using a rate-constant parameterisation; Choi 2011
    # neither re-estimates nor reports these typical values inline -- the
    # numeric values therefore come from Kang 2007 and are fixed here.
    # IIV on the PK parameters is set to zero (Choi 2011 used per-subject
    # individual Bayesian PK estimates as fixed inputs to the PD fit).
    # ====================================================================

    lvc <- fixed(log(9.95))
    label("Central volume of distribution V1 (L; reference V1 of Kang 2007 cohort)")
    # Kang 2007 Table 2: V1 = 9.95 L (theta1), bootstrap median 10.00 (95% CI 7.95-17.00)

    lk10_int <- fixed(log(0.3))
    label("Elimination rate constant k10 intercept (1/min) at AGE = 0 and BSA = 0")
    # Kang 2007 Table 2: theta2 = 0.3 1/min; intercept of additive form
    # k10 = theta2 - theta9 * (AGE / 100) + theta10 * BSA

    e_age_k10 <- fixed(-0.0939)
    label("Additive AGE-on-k10 coefficient (1/min per unit of AGE / 100)")
    # Kang 2007 Table 2: theta9 = 0.0939 with negative sign in additive form
    # k10 = theta2 - theta9 * (AGE / 100) + theta10 * BSA

    e_bsa_k10 <- fixed(0.0491)
    label("Additive BSA-on-k10 coefficient (1/min per m^2)")
    # Kang 2007 Table 2: theta10 = 0.0491 1/min per m^2 of BSA

    lk12 <- fixed(log(0.159))
    label("Central -> peripheral1 distribution rate k12 (1/min)")
    # Kang 2007 Table 2: k12 = 0.159 1/min (theta3)

    lk21 <- fixed(log(0.136))
    label("Peripheral1 -> central return rate k21 (1/min)")
    # Kang 2007 Table 2: k21 = 0.136 1/min (theta4)

    lk13 <- fixed(log(0.0185))
    label("Central -> peripheral2 distribution rate k13 (1/min)")
    # Kang 2007 Table 2: k13 = 0.0185 1/min (theta5)

    lk31 <- fixed(log(0.00204))
    label("Peripheral2 -> central return rate k31 (1/min)")
    # Kang 2007 Table 2: k31 = 0.00204 1/min (theta6)

    # ====================================================================
    # PD PARAMETERS (CHOI 2011 TABLE 4, 'Combined effect and tolerance model'
    # column). TLMC is baseline-normalised so baseline E0 = 1 by construction
    # and is not estimated. The %CV column in Table 4 reports IIV on the
    # log scale; we convert each via omega^2 = log(1 + CV^2).
    # ====================================================================

    lemax <- log(0.48)
    label("Maximum fractional Emax of the sigmoid effect / tolerance terms (unitless, in [0,1])")
    # Choi 2011 Table 4 effect-tolerance col: Emax = 0.48, RSE 12.6%, IIV %CV 27.3

    lec50 <- log(6.02)
    label("Effect-compartment half-maximal concentration CE50 (ng/mL)")
    # Choi 2011 Table 4 effect-tolerance col: CE50 = 6.02 ng/mL, RSE 43.7%, IIV %CV 102

    lct50 <- log(1.96)
    label("Tolerance-compartment half-maximal concentration CT50 (ng/mL)")
    # Choi 2011 Table 4 effect-tolerance col: CT50 = 1.96 ng/mL, RSE 30.6%, IIV %CV 136

    lhill <- log(3.72)
    label("Hill (gamma) coefficient of the sigmoid effect / tolerance terms (unitless)")
    # Choi 2011 Table 4 effect-tolerance col: gamma = 3.72, RSE 31.5%, IIV %CV 74.0;
    # the same gamma is applied to both the effect and the tolerance sigmoid in
    # Choi 2011 Methods / Appendix 3 control stream.

    lke0 <- log(0.62)
    label("Effect-compartment equilibration rate ke0 (1/min)")
    # Choi 2011 Table 4 effect-tolerance col: ke0 = 0.62 1/min, RSE 35.0%, IIV %CV 119;
    # corresponding effect-site equilibration half-life log(2)/ke0 = 1.1 min
    # (Choi 2011 Results 'Pharmacodynamic modelling' confirms 'equilibration half-life
    # for the delay of effect ... in a typical person was 1.1 min').

    lkt0 <- log(0.0033)
    label("Tolerance-compartment equilibration rate kt0 (1/min)")
    # Choi 2011 Table 4 effect-tolerance col: kt0 = 0.0033 1/min, RSE 45.2%, IIV %CV 133;
    # corresponding tolerance-development half-life log(2)/kt0 = 210 min = 3.5 h
    # (Choi 2011 Results 'Pharmacodynamic modelling' confirms 'tolerance development
    # in a typical person was ... 3.5 h').

    # IIV (Choi 2011 Table 4 %CV column -> omega^2 = log(1 + CV^2)).
    # Inter-individual random variability and residual random variability were
    # modelled using a log-normal model and an additive error model respectively
    # (Choi 2011 Table 4 footnote and Methods 'Pharmacodynamic modelling').
    etalemax ~ log(1 + 0.273^2)   # Choi 2011 Table 4: Emax IIV %CV 27.3
    etalec50 ~ log(1 + 1.02^2)    # Choi 2011 Table 4: CE50 IIV %CV 102
    etalct50 ~ log(1 + 1.36^2)    # Choi 2011 Table 4: CT50 IIV %CV 136
    etalhill ~ log(1 + 0.740^2)   # Choi 2011 Table 4: gamma IIV %CV 74.0
    etalke0  ~ log(1 + 1.19^2)    # Choi 2011 Table 4: ke0  IIV %CV 119
    etalkt0  ~ log(1 + 1.33^2)    # Choi 2011 Table 4: kt0  IIV %CV 133

    # Residual error (additive on baseline-normalised TLMC scale).
    addSd <- sqrt(0.002)
    label("Additive residual SD on baseline-normalised TLMC (unitless)")
    # Choi 2011 Table 4 effect-tolerance col: sigma^2 = 0.002 (additive); SD = sqrt(0.002) = 0.04472
  })

  model({
    # -------------------------------------------------------------------
    # PK individual parameters (no IIV: Choi 2011 used per-subject
    # Bayesian PK estimates from Kang 2007 as fixed inputs to the PD fit,
    # so the typical population values are the only available structural
    # PK driving the PD here).
    # -------------------------------------------------------------------
    vc <- exp(lvc)

    # Kang 2007 Table 2 additive form for the elimination rate constant:
    #   k10 = theta2 - theta9 * (AGE / 100) + theta10 * BSA
    # The covariate effects e_age_k10 / e_bsa_k10 are applied additively
    # (NOT in exp(...) form) because Kang 2007's published formulation is
    # additive on the linear k10 scale.
    k10 <- exp(lk10_int) + e_age_k10 * (AGE / 100) + e_bsa_k10 * BSA
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)

    # Canonical CL-form expressions of the above rate constants
    # (downstream interpretation; not used in the ODEs below).
    cl  <- k10 * vc
    q   <- k12 * vc
    q2  <- k13 * vc
    vp  <- q  / k21
    vp2 <- q2 / k31

    # -------------------------------------------------------------------
    # PD individual parameters
    # -------------------------------------------------------------------
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    ct50 <- exp(lct50 + etalct50)
    hill <- exp(lhill + etalhill)
    ke0  <- exp(lke0  + etalke0)
    kt0  <- exp(lkt0  + etalkt0)

    # -------------------------------------------------------------------
    # ODE system. The 3-cmt IV remifentanil PK is in rate-constant form
    # per Kang 2007 (paper-faithful). The effect and tolerance compartments
    # track concentrations directly via the standard rxode2 effect-
    # compartment idiom with nominal unit volume (equivalent to Choi 2011's
    # NONMEM V4 = V5 = 0.0001 trick which makes mass loss into the
    # effect / tolerance compartments negligible).
    # -------------------------------------------------------------------
    d/dt(central) <- -k10 * central -
                      k12 * central - k13 * central +
                      k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    d/dt(effect)    <- ke0 * (Cc - effect)
    d/dt(tolerance) <- kt0 * (Cc - tolerance)

    # -------------------------------------------------------------------
    # Observation (Choi 2011 Methods / Figure 4 caption / Appendix 3
    # control stream IPRED line). Baseline-normalised TLMC. The negative
    # term is the effect-compartment-driven depression of TLMC; the
    # positive term is the tolerance-compartment-driven rebound that
    # appears as central drug clears faster than tolerance develops.
    #   net_effect = 1 - Emax * CE^gamma / (CE^gamma + CE50^gamma)
    #                + Emax * CT^gamma / (CT^gamma + CT50^gamma)
    # -------------------------------------------------------------------
    ce_gam <- effect^hill
    ct_gam <- tolerance^hill
    tlmc <- 1 -
      emax * ce_gam / (ce_gam + ec50^hill) +
      emax * ct_gam / (ct_gam + ct50^hill)

    tlmc ~ add(addSd)
  })
}
