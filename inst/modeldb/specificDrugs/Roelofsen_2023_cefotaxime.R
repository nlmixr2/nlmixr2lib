Roelofsen_2023_cefotaxime <- function() {
  description <- "Two-compartment intravenous population PK model for cefotaxime given as pre-emptive treatment (selective digestive decontamination) in critically ill adult ICU patients (Roelofsen 2023). Clearance scales as a power function of CKD-EPI estimated glomerular filtration rate (exponent 0.477, reference 57 mL/min/1.73 m^2) and of serum albumin (exponent 0.640, reference 26 g/L); together the two covariates explain 48% of the between-subject variability in clearance. Between-subject variability is carried on clearance, central volume and intercompartmental clearance, and residual variability is a combined additive-plus-proportional error model."
  reference <- "Roelofsen EE, Abdulla A, Muller AE, Endeman H, Gommers D, Dijkstra A, Hunfeld NGM, de Winter BCM, Koch BCP. Dose optimization of cefotaxime as pre-emptive treatment in critically ill adult patients: A population pharmacokinetic study. Br J Clin Pharmacol. 2023;89(2):705-713. doi:10.1111/bcp.15487"
  vignette <- "Roelofsen_2023_cefotaxime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefotaxime is given as a short IV infusion, so the dose
  # enters `central` directly and there is no depot state. `central` is
  # verified: Roelofsen 2023 Sect. 2.2 states that TOTAL CEFOTAXIME PLASMA
  # concentrations were measured by a validated UPLC-MS/MS assay, and the
  # supplementary control stream (S1) declares `COMP=(CENTRAL, DEFOBS)` with
  # `S1 = V1`, so the observed compartment is the central one. `peripheral1`
  # is left unverified because the paper never states what matrix the
  # peripheral distribution compartment represents; "plasma" follows the
  # repository default for a mathematical peripheral compartment and is not a
  # paper-sourced claim.
  compartmentData <- list(
    central     = list(analyte = "cefotaxime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefotaxime", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the Chronic Kidney Disease Epidemiology Collaboration (CKD-EPI) equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roelofsen 2023 Sect. 2.3 names the CKD-EPI equation as the estimating formula. Table 1 reports the cohort median 57 mL/min/1.73 m^2 (range 4-347), and 57 is exactly the normalizing constant printed in the Sect. 3.2 clearance equation, consistent with Sect. 2.6 ('Continuous covariates were normalized to the population median'). The cohort spans severe renal impairment through augmented renal clearance. One patient had an eGFR above 300 mL/min/1.73 m^2; the Discussion reports that capping that subject at 141 (the second-highest value) did not markedly change the estimates. The Sect. 2.7 Monte Carlo simulations evaluated only 10, 30, 50, 80 and 100 mL/min/1.73 m^2, explicitly to avoid extrapolating beyond the range in which the covariate predominantly occurred. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; the assay form here is the creatinine-based CKD-EPI estimate. Note that the paper's abstract glosses eGFR as '(creatinine clearance)', but Sect. 2.3 is the authority and specifies CKD-EPI.",
      source_name        = "eGFR"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roelofsen 2023 Table 1 reports the cohort median 26 g/L (range 11-47), and 26 is exactly the normalizing constant printed in the Sect. 3.2 clearance equation. The paper reports albumin in SI g/L, which is already the canonical unit, so no conversion is applied in model(). The effect is POSITIVE (higher albumin gives higher clearance), which the Discussion notes is the opposite of the direction expected from protein-binding displacement for a drug with only ~30% protein binding; the authors' preferred interpretation is that higher albumin marks less severe illness and therefore fewer physiological changes affecting PK, while noting that SOFA and APACHE II scores did not themselves reach significance on clearance. The Sect. 2.7 simulations evaluated only 20, 30 and 40 g/L, to avoid extrapolating beyond the range in which the covariate predominantly occurred.",
      source_name        = "albumin"
    )
  )

  # Screened in the Roelofsen 2023 covariate analysis (Sect. 2.6 lists the
  # candidate set; Sect. 3.2 reports "No other significant covariates were
  # found") but NOT retained in the final model, so these are documentation
  # only and are not referenced in model(). Forward inclusion used p < 0.05
  # (dOFV 3.84) and backward elimination p < 0.001 (dOFV 10.83).
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. Table 1 median 64 years (range 23-85).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened per Sect. 2.3/2.6 (supplementary control stream column `GEN`); not retained. Table 1 reports 57 male / 35 female.",
      source_name        = "GEN"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6 (supplementary control stream column `WGT`); not retained. Table 1 median 76 kg (range 45-150). The final model therefore carries NO body-size term, so clearance and volumes are absolute values for a typical ICU adult rather than weight-normalized.",
      source_name        = "WGT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.6 (supplementary control stream column `BMI`); not retained. Table 1 median 26 kg/m^2 (range 17.8-46.3).",
      source_name        = "BMI"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6 (supplementary control stream column `CRE`); not retained as a standalone covariate, the CKD-EPI eGFR built from it having entered instead. Table 1 median 98 umol/L (range 5-913).",
      source_name        = "CRE"
    ),
    UREA = list(
      description        = "Serum urea",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6 (supplementary control stream column `URE`); not retained. No summary value is tabulated.",
      source_name        = "URE"
    ),
    CRP = list(
      description        = "C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. Table 1 median 127 mg/L (range 0-488).",
      source_name        = "CRP"
    ),
    WBC = list(
      description        = "White blood cell (leucocyte) count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. Table 1 median 13 x 10^9 cells/L (range 0.9-100).",
      source_name        = "WBC"
    ),
    TEMP = list(
      description        = "Body temperature",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. No summary value is tabulated.",
      source_name        = "TEMP"
    ),
    SOFA = list(
      description        = "Sequential Organ Failure Assessment score",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. The Discussion states explicitly that SOFA showed NO change in objective function value during forward covariate analysis. Table 1 median 13 (range 1-21).",
      source_name        = "SOF"
    ),
    APACHE2 = list(
      description        = "Acute Physiology and Chronic Health Evaluation II score",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6; not retained. The Discussion states explicitly that APACHE II produced an INCREASE in objective function value during forward covariate analysis. Table 1 median 23 (range 7-71).",
      source_name        = "APA"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on CRRT)",
      notes              = "Screened per Sect. 2.3/2.6 as a binary covariate; not retained. Table 1 reports 5 of 92 patients (5.4%) on CRRT. The Discussion records two important caveats: CRRT was registered only at BASELINE because duration and continuation during sampling were not captured, and excluding the 5 CRRT patients did not markedly influence the PK estimates. The authors state the effect of CRRT on cefotaxime PK needs further investigation, so this model should not be used to describe patients on renal replacement therapy.",
      source_name        = "EPI"
    ),
    FLUIDBAL = list(
      description        = "Cumulative fluid balance",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Sect. 2.3/2.6 (supplementary control stream column `VBL`); not retained. No summary value is tabulated.",
      source_name        = "VBL"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 92L,
    n_studies        = 1L,
    n_sites          = 2L,
    n_concentrations = 437L,
    age_range        = "23-85 years (Table 1 median 64)",
    age_median       = "64 years (range 23-85)",
    weight_range     = "45-150 kg",
    weight_median    = "76 kg (range 45-150)",
    sex_female_pct   = 38.0,
    disease_state    = "Critically ill adults admitted to the intensive care unit with an expected stay of more than 72 hours, receiving intravenous cefotaxime as the systemic component of selective digestive decontamination (SDD) with or without additional treatment. The cohort is deliberately heterogeneous (general ICU plus trauma). Excluded: patients under 18 years, admission for burn wounds, cefotaxime discontinued before sampling, and absence of written informed consent. Table 1 severity: APACHE II median 23 (range 7-71); SOFA median 13 (range 1-21); C-reactive protein median 127 mg/L (range 0-488); leucocytes median 13 x 10^9 cells/L (range 0.9-100). Five patients (5.4%) were on CRRT at baseline.",
    dose_range       = "Cefotaxime 1 g intravenously every 6 h (80 patients) or every 4 h (12 patients), per the SDD protocol and at the discretion of the attending physician. Infusion durations in the study ranged from 1 minute to 1 hour.",
    regions          = "Netherlands (Erasmus University Medical Center and Maasstad Hospital, Rotterdam)",
    renal_function   = "CKD-EPI eGFR median 57 mL/min/1.73 m^2 (range 4-347); serum creatinine median 98 umol/L (range 5-913). The cohort spans severe renal impairment through augmented renal clearance.",
    notes            = "Prospective observational PK/PD sub-study of the EXPAT trial (Netherlands Trial Registry NTR 5632), enrolling January 2016 to June 2017. 93 patients were enrolled and 1 was excluded for a physiologically impossible concentration-time profile, leaving 92. Five samples per patient were drawn within a single dosing interval on day 2 of therapy (15-30 min pre-dose, 15-30 min post-administration, 1 h and 3 h after end of infusion, and immediately pre-next-dose); of 453 analysed samples, 16 were excluded and 7 were not drawn, giving 437 observations. Total (not free) plasma cefotaxime was assayed by validated UPLC-MS/MS with a calibration range of 0.25-12.5 mg/L. Two samples (0.5%) were below the limit of quantification and were dropped rather than imputed. Model fit in NONMEM 7.4.2 with FOCE-INTERACTION on untransformed data (supplementary control stream S1 uses ADVAN5). Covariate selection by forward inclusion at p < 0.05 (dOFV 3.84) then backward elimination at p < 0.001 (dOFV 10.83). Reported eta shrinkage in the final model: 2.2% on CL, 13.3% on Vc, 17.3% on Q. The desacetylcefotaxime metabolite was NOT measured; the authors argue its contribution is about 5% of cefotaxime's antimicrobial activity. A protein binding of 30% was assumed for the target-attainment simulations (taken from Aardema et al.), so free concentrations are 0.70 times the total concentrations this model predicts."
  )

  ini({
    # Structural parameters, final model including covariates (Roelofsen 2023
    # Table 2, "Final model including covariates" column). The reference
    # subject has an eGFR of 57 mL/min/1.73 m^2 and an albumin of 26 g/L,
    # which are the Table 1 cohort medians and the normalizing constants of
    # the Sect. 3.2 clearance equation.
    lcl <- log(7.08);  label("Clearance at CRCL=57 mL/min/1.73 m^2 and ALB=26 g/L (L/h)")  # Roelofsen 2023 Table 2 final model: CL 7.08 L/h (RSE 5.4%)
    lvc <- log(15.70); label("Central volume of distribution (L)")                          # Roelofsen 2023 Table 2 final model: V1 15.70 L (RSE 6.2%)
    lvp <- log(25.00); label("Peripheral volume of distribution (L)")                       # Roelofsen 2023 Table 2 final model: V2 25.00 L (RSE 37.0%)
    lq  <- log(4.81);  label("Intercompartmental clearance (L/h)")                          # Roelofsen 2023 Table 2 final model: Q 4.81 L/h (RSE 15.2%)

    # Covariate effects on clearance. Roelofsen 2023 Sect. 3.2 prints the
    # equation directly:
    #   CL = 7.08 * (eGFR/57)^0.477 * (albumin concentration/26)^0.64
    # Sect. 2.6 confirms the form ("Continuous covariates were normalized to
    # the population median and implemented by use of an exponential model").
    # Both are estimated, not fixed, so neither is wrapped in fixed().
    e_crcl_cl <- 0.477; label("Power exponent on (CRCL/57) for CL (unitless)")  # Roelofsen 2023 Table 2 final model: "Covariate eGFR on CL" 0.477 (RSE 15.7%); Sect. 3.2 equation
    e_alb_cl  <- 0.640; label("Power exponent on (ALB/26) for CL (unitless)")   # Roelofsen 2023 Table 2 final model: "Covariate albumin concentration on CL" 0.640 (RSE 24.8%); Sect. 3.2 equation prints 0.64

    # Between-subject variability, final model (Roelofsen 2023 Table 2,
    # reported as a percentage). The supplementary control stream (S1)
    # specifies exponential BSV -- CL = THETA(1)*EXP(ETA(1)), V1 =
    # THETA(2)*EXP(ETA(2)), Q = THETA(4)*EXP(ETA(3)) -- on exactly these three
    # parameters and none on V2, matching Sect. 3.2.
    #
    # The paper does not state which %CV convention Table 2 uses, but Sect. 3.2
    # does state that "the covariates could explain 48% of the IIV on
    # clearance", which pins it. Reading the percentages as the NONMEM
    # approximation omega = CV/100 gives a variance reduction of
    # (0.696^2 - 0.503^2)/0.696^2 = 47.8%, which rounds to the stated 48%.
    # Reading them as CV = sqrt(exp(omega^2)-1) instead gives 42.9%, which does
    # not. So omega = CV/100 and omega^2 = (CV/100)^2 is used below.
    #
    # Sect. 3.2 and the abstract both state that correlations between the IIV
    # were included via an OMEGA BLOCK, and the supplementary control stream
    # confirms `$OMEGA BLOCK(3)`. The final off-diagonal estimates are NOT
    # published anywhere in the paper or its supplement -- only the initial
    # values (0.001) appear in the control stream -- so the etas are declared
    # DIAGONAL here (zero covariance) rather than with invented correlations.
    # This is the one structural deviation from the published model; see the
    # vignette "Assumptions and deviations" section.
    etalcl ~ 0.253009  # 0.503^2; Roelofsen 2023 Table 2 final model: variability on CL 50.3% (RSE 11.8%, shrinkage 2.2%)
    etalvc ~ 0.121801  # 0.349^2; Roelofsen 2023 Table 2 final model: variability on V1 34.9% (RSE 16.1%, shrinkage 13.3%)
    etalq  ~ 0.848241  # 0.921^2; Roelofsen 2023 Table 2 final model: variability on Q  92.1% (RSE 13.8%, shrinkage 17.3%)

    # Residual variability. Sect. 2.5 specifies "a combined proportional and
    # additive model". The supplementary control stream writes it as
    #   Y = F + THETA(5)*EPS(1) + F*THETA(6)*EPS(2)
    # with both SIGMAs FIXED to 1, so THETA(5) and THETA(6) are the additive
    # and proportional residual STANDARD DEVIATIONS directly, in the same
    # scale as Table 2's "Additive error (mg L-1)" and "Proportional error"
    # rows. No variance-to-SD conversion is needed.
    addSd  <- 0.617; label("Additive residual SD (mg/L)")        # Roelofsen 2023 Table 2 final model: additive error 0.617 mg/L (RSE 25.9%)
    propSd <- 0.191; label("Proportional residual SD (fraction)")  # Roelofsen 2023 Table 2 final model: proportional error 0.191 (RSE 8.6%)
  })
  model({
    # Individual PK parameters. Roelofsen 2023 Sect. 3.2:
    #   CL = 7.08 * (eGFR/57)^0.477 * (albumin/26)^0.64
    # V1, V2 and Q carry no covariate effects (Sect. 3.2: "No other
    # significant covariates were found").
    cl <- exp(lcl + etalcl) * (CRCL / 57)^e_crcl_cl * (ALB / 26)^e_alb_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L. This is the TOTAL
    # (protein-bound plus free) plasma cefotaxime concentration that the
    # Sect. 2.2 assay measured; the paper's target-attainment work multiplies
    # it by 0.70 to obtain the free concentration.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
