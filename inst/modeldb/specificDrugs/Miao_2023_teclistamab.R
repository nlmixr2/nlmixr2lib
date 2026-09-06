Miao_2023_teclistamab <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order",
    "subcutaneous absorption and parallel time-independent plus",
    "time-dependent elimination for teclistamab, a B-cell maturation antigen",
    "(BCMA) x CD3 bispecific antibody, in patients with relapsed/refractory",
    "multiple myeloma (Miao 2023, MajesTEC-1 phase I/II; 4840 serum",
    "concentrations from 338 patients dosed intravenously or subcutaneously).",
    "Total clearance is CL(t) = CL1 + CL2 * exp(-KDES * t), where t is time",
    "since the start of teclistamab treatment; the time-dependent component",
    "CL2 contributes about 43% of total clearance at treatment start and",
    "decays with a 0.0292/day first-order rate constant. Miao 2023 describes",
    "this as an approximation of target-mediated drug disposition combined",
    "with a feedback loop from improving disease status onto clearance",
    "(quasi-steady-state TMDD and time-varying soluble BCMA parameterisations",
    "were tested and rejected). Retained covariates are body weight on CL1,",
    "Vc and Vp; ISS stage (II and III vs I) on CL1; and myeloma type",
    "(non-IgG vs IgG) on both CL1 and CL2. The paper's exposure-response",
    "analyses for overall response rate, duration of response,",
    "progression-free survival, overall survival and grade >= 3 adverse",
    "events are graphical / logistic-regression explorations whose",
    "coefficients are not reported; only the population PK model is",
    "reproducible and only it is encoded here."
  )
  reference <- paste(
    "Miao X, Wu LS, Wang Lin SX, Xu Y, Chen Y, Iwaki Y, Kobos R, Stephenson T,",
    "Kemmerer K, Uhlar CM, Banerjee A, Goldberg JD, Trancucci D, Apte A,",
    "Verona R, Pei L, Desai R, Hickey K, Su Y, Ouellet D, Samtani MN, Guo Y,",
    "Garfall AL, Krishnan A, Usmani SZ, Zhou H, Girgis S.",
    "Population pharmacokinetics and exposure-response with teclistamab in",
    "patients with relapsed/refractory multiple myeloma: results from",
    "MajesTEC-1. Target Oncol. 2023;18(5):667-684.",
    "doi:10.1007/s11523-023-00989-z"
  )
  vignette <- "Miao_2023_teclistamab"

  # Miao 2023 reports clearances in L/day, volumes in L, Ka and KDES in
  # 1/day, and serum concentrations in ug/mL (Fig. 2 and Fig. 3 axes).
  # Doses are in mg (weight-based mg/kg regimens are converted to mg per
  # patient). With dose in mg and Vc in L, central/vc is mg/L = ug/mL, so
  # no conversion factor is needed in the observation equation.
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "teclistamab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "teclistamab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "teclistamab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric-style) effect normalised to the cohort median of 74 kg, with separately estimated exponents on the time-independent clearance component (0.704), the central volume (0.358) and the peripheral volume (1.40). Miao 2023 Table 2 footnotes a, c and d. Body weight was carried into the structural model a priori rather than by stepwise selection, because it is an established covariate for IgG-based monoclonal antibodies (Methods 2.3). No exponent was applied to Q, which has no covariates. Cohort median 74.3 kg, range 41.0-139 kg (Table 1). Note the reference value is the rounded 74 kg that appears in the printed equations, not the tabulated 74.3 kg median.",
      source_name        = "BWT"
    ),
    ISS_II = list(
      description        = "Baseline International Staging System stage II multiple myeloma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ISS stage I, the model reference level; ISS_II and ISS_III are mutually exclusive so an ISS stage III patient has ISS_II = 0)",
      notes              = "Multiplicative factor 1.31 on the time-independent clearance component, encoded exactly as printed in Miao 2023 Table 2 footnote a: CL1 = 0.449 * (BWT/74)^0.704 * 1.31^(ISS = II) * 1.67^(ISS = III) * 0.689^(TPMM = Non-IgG). Higher clearance, hence lower exposure, for ISS stage II than stage I; the paper reports a corresponding ~17% lower geometric-mean Cave,1stdose for stage II vs stage I (Results 3.5). Baseline ISS was derived from serum beta-2-microglobulin and albumin (Table 1 footnote b). Cohort split: stage I 165 (48.8%), II 109 (32.2%), III 59 (17.5%), not reported 5 (1.5%).",
      source_name        = "ISS"
    ),
    ISS_III = list(
      description        = "Baseline International Staging System stage III multiple myeloma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ISS stage I, the model reference level; ISS_II and ISS_III are mutually exclusive so an ISS stage II patient has ISS_III = 0)",
      notes              = "Multiplicative factor 1.67 on the time-independent clearance component (Miao 2023 Table 2 footnote a). The paper reports a corresponding ~29% lower geometric-mean Cave,1stdose for stage III vs stage I (Results 3.5); that value is reproduced by this model because CL1 carries only part of the total clearance at the time Cave,1stdose is measured. See ISS_II for the shared staging definition.",
      source_name        = "ISS"
    ),
    MM_NIGG = list(
      description        = "Non-IgG-secreting multiple myeloma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IgG multiple myeloma) -- Miao 2023 anchors the typical values of CL1 and CL2 to IgG myeloma and applies a reduction factor to non-IgG patients, the same reference orientation as Fau_2020_isatuximab and the opposite of Xu_2020_daratumumab. Canonical column semantics (1 = non-IgG) are preserved.",
      notes              = "The only covariate retained on BOTH clearance components: factor 0.689 on the time-independent component CL1 (Table 2 footnote a) and 0.295 on the time-dependent component CL2 (Table 2 footnote b). Non-IgG myeloma therefore lowers total clearance and raises exposure; the paper reports geometric-mean Cave,1stdose ~32% lower in IgG than in non-IgG patients (Results 3.5). Mechanistic rationale (Discussion): endogenous monoclonal IgG in IgG-myeloma patients saturates FcRn-mediated recycling and accelerates catabolism of the therapeutic IgG. Cohort split: IgG 173 (51.2%), non-IgG 165 (48.8%) (Table 1).",
      source_name        = "TPMM"
    )
  )

  # Screened in the stepwise covariate search (Methods 2.3) but NOT retained
  # in the final model, so they are documented rather than encoded. Miao 2023
  # reports no point estimate for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened as a demographic covariate; not retained. Cohort median 64 years (range 24-84)."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "Screened as a demographic covariate; not retained. Cohort 44.4% female."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance", units = "mL/min", type = "continuous",
      notes = "Screened as a clinical covariate; not retained. Renal function distribution in Table 1; only 1 of 338 patients had eGFR < 30 mL/min/1.73 m2."
    ),
    ALB = list(
      description = "Baseline serum albumin", units = "g/dL", type = "continuous",
      notes = "Screened as a clinical covariate; not retained in its own right. Albumin does enter the model indirectly, because baseline ISS stage is derived from serum beta-2-microglobulin and albumin (Table 1 footnote b)."
    ),
    SBCMA = list(
      description = "Soluble B-cell maturation antigen, baseline and time-varying", units = "ng/mL", type = "continuous",
      notes = "Screened both as a baseline covariate and as a time-varying driver of the time-dependent clearance; neither was retained (models relating time-dependent clearance to time-varying sBCMA had significant increases in objective function value). Baseline sBCMA WAS a significant prognostic factor for overall response rate in the multivariate logistic regression (odds ratio 0.99 per ng/mL), which is an exposure-response finding and not part of this PK model. Discussion and Results 3.6."
    ),
    ADA_POS = list(
      description = "Anti-teclistamab antibody positivity", units = "(binary)", type = "binary",
      notes = "Assessed both as a subject-level and as a time-varying covariate in a sensitivity analysis on the final model; no significant effect on PK. Only 2 of 321 ADA-evaluable patients developed anti-drug antibodies, both at a titre of 20 (Discussion)."
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator", units = "(binary)", type = "binary",
      notes = "Screened as a disease characteristic; not retained. Cohort: ECOG 0 in 121 (35.8%), 1 in 216 (63.9%), 3 in 1 (0.3%)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 338L,
    n_studies      = 1L,
    age_range      = "24-84 years",
    age_median     = "64 years",
    weight_range   = "41.0-139 kg",
    weight_median  = "74.3 kg",
    sex_female_pct = 44.4,
    race_ethnicity = c(White = 83.4, Black = 9.2, Asian = 1.8, Other = 5.6),
    disease_state  = "Relapsed/refractory multiple myeloma diagnosed per International Myeloma Working Group criteria. Phase I patients were relapsed, refractory or intolerant to established therapies; phase II patients had received at least three prior lines including an immunomodulatory agent, a proteasome inhibitor and an anti-CD38 antibody (cohort A) or additionally prior anti-BCMA therapy (cohort C). Baseline ISS stage I 48.8%, II 32.2%, III 17.5%; IgG myeloma 51.2%; triple-refractory 79.0%; >3 prior lines of therapy 79.9%.",
    dose_range     = "Intravenous 0.0003-0.0192 mg/kg Q2W and up to 0.01/0.06/0.24 then 0.72 mg/kg QW; subcutaneous 0.02 mg/kg up to 0.06/0.3/1.5 then 6 mg/kg QW, plus flat 2/6/30 then 150/300 mg QW/Q2W. Recommended phase II dose 1.5 mg/kg SC weekly preceded by step-up doses of 0.06 and 0.3 mg/kg.",
    regions        = "Multicenter international (not tabulated by region in Miao 2023).",
    renal_function = "eGFR >= 90 in 29.9%, 60 to < 90 in 45.3%, 30 to < 60 in 24.6%, < 30 in 0.3% (mL/min/1.73 m2).",
    hepatic_function = "Normal in 88.2%; mild impairment in 11.8% (NCI Organ Dysfunction Working Group criteria). No moderate or severe hepatic impairment.",
    notes          = "Baseline demographics from Miao 2023 Table 1, 'Combined IV and SC dosing' column (n = 338). The analysis dataset held 4840 measurable serum concentrations: 83 patients / 1976 observations intravenous and 255 patients / 2864 observations subcutaneous. MajesTEC-1 is a single phase I/II first-in-human trial carrying two ClinicalTrials.gov registrations (NCT03145181 for phase I, NCT04557098 for phase II), which is why n_studies is 1. PK data cutoff 1 December 2021."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. All point estimates are Miao 2023 Table 2
    # ("Parameter estimates of teclistamab for the final population PK
    # model"); the covariate equations are Table 2 footnotes a-d, read from
    # the MathML in the EuropePMC JATS full text of PMC10518021 (the
    # pdftotext / docling renderings of footnote c drop the division bar
    # and print "V1 = 4.13 x (BWT)^0.358 / 74" on two lines).
    #
    # Total clearance decomposition (Table 2 footnote b):
    #
    #     CL(t) = CL1 + CL2 * exp(-KDES * t),  t in days
    #
    # mapped onto the registered cl_exp_ family:
    #
    #     cl_exp_inf       = CL1  = 0.449 L/day  (time-independent arm)
    #     cl_exp_component = CL2  = 0.547 L/day  (decaying arm, value at t = 0)
    #     cl_exp_kdes      = KDES = 0.0292 1/day
    #
    # Every parameter in Table 2 was estimated; none carries a FIX flag or a
    # 0% relative standard error, so nothing here is wrapped in fixed().
    # ------------------------------------------------------------------
    lcl_exp_inf       <- log(0.449);  label("Time-independent clearance component CL1 at 74 kg, ISS I, IgG myeloma (L/day)")       # Table 2: CL1 = 0.449 L/day, RSE 8.87%
    lcl_exp_component <- log(0.547);  label("Time-dependent clearance component CL2 at t = 0, IgG myeloma (L/day)")               # Table 2: CL2 = 0.547 L/day, RSE 15.6%
    lcl_exp_kdes      <- log(0.0292); label("First-order rate constant for the decay of CL2 over time (1/day)")                   # Table 2: KDES = 0.0292 1/day, RSE 13.0%
    lvc               <- log(4.13);   label("Central volume of distribution V1 at 74 kg (L)")                                     # Table 2: V1 = 4.13 L, RSE 4.40%
    lvp               <- log(1.34);   label("Peripheral volume of distribution V2 at 74 kg (L)")                                  # Table 2: V2 = 1.34 L, RSE 26.1%
    lq                <- log(0.0390); label("Intercompartmental clearance Q (L/day)")                                             # Table 2: Q = 0.0390 L/day, RSE 55.5%
    lka               <- log(0.133);  label("First-order subcutaneous absorption rate constant Ka (1/day)")                       # Table 2: Ka = 0.133 1/day, RSE 7.73%
    lfdepot           <- log(0.718);  label("Subcutaneous bioavailability F (fraction)")                                          # Table 2: Bioavailability = 0.718, RSE 7.38%

    # ------------------------------------------------------------------
    # Covariate effects. Body weight enters as a power model normalised to
    # 74 kg; ISS stage and myeloma type enter as multiplicative factors
    # raised to their 0/1 indicator, exactly as printed in the footnotes.
    # ------------------------------------------------------------------
    e_wt_cl_exp_inf         <- 0.704;  label("Body-weight exponent on CL1, normalised to 74 kg (unitless)")                       # Table 2 footnote a, RSE 21.8%
    e_wt_vc                 <- 0.358;  label("Body-weight exponent on V1, normalised to 74 kg (unitless)")                        # Table 2 footnote c, RSE 60.9%
    e_wt_vp                 <- 1.40;   label("Body-weight exponent on V2, normalised to 74 kg (unitless)")                        # Table 2 footnote d, RSE 25.5%
    e_iss_ii_cl_exp_inf     <- 1.31;   label("Multiplicative factor on CL1 for ISS stage II vs stage I (unitless)")               # Table 2 footnote a, RSE 7.83%
    e_iss_iii_cl_exp_inf    <- 1.67;   label("Multiplicative factor on CL1 for ISS stage III vs stage I (unitless)")              # Table 2 footnote a, RSE 11.1%
    e_nigg_cl_exp_inf       <- 0.689;  label("Multiplicative factor on CL1 for non-IgG vs IgG myeloma (unitless)")                # Table 2 footnote a, RSE 7.76%
    e_nigg_cl_exp_component <- 0.295;  label("Multiplicative factor on CL2 for non-IgG vs IgG myeloma (unitless)")                # Table 2 footnote b, RSE 21.6%

    # ------------------------------------------------------------------
    # Interindividual variability. Table 2 reports IIV only as a percent
    # coefficient of variation, one value per block of rows; per Methods
    # 3.2 the four IIV terms sit on CL1, V1, CL2 and Ka. IIV on Q, KDES and
    # F was attempted and dropped (non-convergent or worse OFV), so those
    # three carry no eta. Variances are the log-normal back-transform
    # omega^2 = log(CV^2 + 1):
    #
    #     CL1 53.6% -> log(0.536^2 + 1) = 0.2526
    #     CL2 107%  -> log(1.07^2  + 1) = 0.7631
    #     V1  48.8% -> log(0.488^2 + 1) = 0.2136
    #     Ka  45.2% -> log(0.452^2 + 1) = 0.1859
    #
    # Miao 2023 does not state whether its CV% column is the exact
    # log-normal CV or the sqrt(omega^2) approximation; see the vignette's
    # Assumptions section. No off-diagonal covariances are reported, so the
    # four etas are encoded as independent.
    # ------------------------------------------------------------------
    etalcl_exp_inf       ~ 0.2526  # Table 2: IIV on CL1, CV 53.6%, RSE 14.3%, shrinkage 14.4%
    etalcl_exp_component ~ 0.7631  # Table 2: IIV on CL2, CV 107%,  RSE 20.5%, shrinkage 33.8%
    etalvc               ~ 0.2136  # Table 2: IIV on V1,  CV 48.8%, RSE 50.6%, shrinkage 29.5%
    etalka               ~ 0.1859  # Table 2: IIV on Ka,  CV 45.2%, RSE 32.1%, shrinkage 44.3%

    # ------------------------------------------------------------------
    # Residual error. Miao 2023 fitted "an additive model in the log-scale"
    # (Results 3.2) and reports it in Table 2 as an additive log-scale term
    # expressed as a percent coefficient of variation, 41.7%. An additive
    # error on the log scale is a proportional error in nlmixr2's linear
    # space with the same coefficient of variation, so propSd = 0.417.
    # ------------------------------------------------------------------
    propSd <- 0.417; label("Proportional residual error (fraction)")  # Table 2: additive error term on the log-scale, CV 41.7%, RSE 4.35%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # The 0/1 covariate indicators are used as EXPONENTS of their
    # multiplicative factors, which is the printed form of Miao 2023
    # Table 2 footnotes a and b and is exact for binary indicators
    # (1.31^0 = 1, 1.31^1 = 1.31). ISS_II and ISS_III are mutually
    # exclusive level indicators, so at most one stage factor is ever
    # active; both zero means ISS stage I, the reference level.
    # ------------------------------------------------------------------
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) *
      (WT / 74)^e_wt_cl_exp_inf *
      e_iss_ii_cl_exp_inf^ISS_II *
      e_iss_iii_cl_exp_inf^ISS_III *
      e_nigg_cl_exp_inf^MM_NIGG

    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component) *
      e_nigg_cl_exp_component^MM_NIGG

    cl_exp_kdes <- exp(lcl_exp_kdes)

    vc <- exp(lvc + etalvc) * (WT / 74)^e_wt_vc
    vp <- exp(lvp) * (WT / 74)^e_wt_vp
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    # ------------------------------------------------------------------
    # 2. Total clearance, decaying towards the time-independent arm.
    # `time` is the rxode2 solver time; with the first teclistamab dose at
    # t = 0 it equals the paper's "Time (in days)" since treatment start.
    # ------------------------------------------------------------------
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * time)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 4. ODE system. Subcutaneous doses go into `depot`; intravenous doses
    # (bolus or infusion) go directly into `central`, which is how the
    # single Miao 2023 model describes both routes.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # 5. Bioavailability applies to the subcutaneous depot only.
    f(depot) <- exp(lfdepot)

    # 6. Observation. Dose in mg and vc in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
