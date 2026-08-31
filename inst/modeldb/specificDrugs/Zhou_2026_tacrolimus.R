Zhou_2026_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for oral immediate-release tacrolimus in adult Chinese patients with nephrotic syndrome, built from steady-state trough (Cmin) therapeutic drug monitoring data. Because only troughs were available the absorption rate constant is fixed at a literature value. Apparent clearance CL/F carries two multiplicative covariate effects: a reduction with concomitant Wuzhi capsule and a reduction in CYP3A5*3/*3 non-expressers; no covariate was retained on apparent volume of distribution. Exponential inter-individual variability on CL/F only, with a proportional residual error. This is the population PK half of a paper whose second half feeds the individual CL/F estimate into a machine-learning ensemble; only the pharmacokinetic model is represented here."
  reference <- paste(
    "Zhou Y, Zhou Z, Chen S, Zhu L, Yun Y, Yuan Y, Chen C, Zou J, Zhao J.",
    "An Integrated Population Pharmacokinetic and Machine Learning Model for",
    "Predicting Tacrolimus Exposure in Adult Patients with Nephrotic Syndrome.",
    "Drug Des Devel Ther. 2026;20.",
    "doi:10.2147/DDDT.S576598.",
    "Parameter values are the 'Final model' Estimate column of Supplementary",
    "Table S1; the structural covariate equation is the displayed equation in",
    "Results, 'Population Pharmacokinetic Model'.",
    sep = " "
  )
  vignette <- "Zhou_2026_tacrolimus"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. MATRIX IS INTERNALLY INCONSISTENT IN THE SOURCE.
  # Zhou 2026 calls the measurement "TAC plasma concentration" throughout
  # (title of Table 1's first row, the abbreviation list, and every Results
  # heading), but Methods "Measurement of TAC Plasma Concentration" states the
  # assay is the Abbott ARCHITECT i1000sr with the "ARCHITECT Tacrolimus
  # Reagent Kit" run as a chemiluminescent microparticle immunoassay -- that
  # kit is a WHOLE-BLOOD assay (tacrolimus partitions heavily into
  # erythrocytes and is universally monitored in whole blood, which is also
  # what the paper's own 5-10 ng/mL therapeutic target refers to). The
  # specimen is therefore recorded as whole blood and `verified` is FALSE
  # because the source contradicts itself; V/F and CL/F are whole-blood
  # apparent parameters either way, so the numerical model is unaffected.
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    CONMED_WUZHI = list(
      description        = "Concomitant Wuzhi capsule (Schisandra sphenanthera extract) indicator (1 = receiving Wuzhi capsule, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant Wuzhi capsule)",
      notes              = "Source label WZC. Zhou 2026 Results states 'If combined with WZC, theta_WZC = 0.731; if not combined with WZC, theta_WZC = 1', so the effect enters as the multiplier 0.731^CONMED_WUZHI on CL/F -- a 26.9% reduction in apparent clearance. 27 of 410 samples (6.6%) were taken under concomitant Wuzhi capsule (Table 1). The paper attributes the direction to inhibition of CYP3A and P-glycoprotein by the Schisandra extract (Discussion). Note that Zhou 2026 spells out WZC as 'compound Salvia miltiorrhiza polyphenolic acid capsule' in Results; this is an error in the source -- Wuzhi capsule is a Schisandra sphenanthera preparation, and the paper's own Discussion and reference 27 both discuss it as such.",
      source_name        = "WZC"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status (1 = carries at least one CYP3A5*1 allele, genotype *1/*1 or *1/*3; 0 = CYP3A5*3/*3 non-expresser)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expresser)",
      notes              = "VALUE INVERSION relative to the source. Zhou 2026 Results codes the opposite orientation: 'if CYP3A5 genotype is *3*3, theta_CYP3A5 = 0.768; if CYP3A5 genotype is *1*1/*1*3, theta_CYP3A5 = 1'. The canonical register mandates the expresser-equals-1 orientation and explicitly instructs papers using a *3/*3 indicator to record values under CYP3A5_EXPR with the inversion documented, so the published multiplier is applied here to (1 - CYP3A5_EXPR). The direction of the effect is unchanged from the paper: non-expressers have 0.768 times the apparent clearance of expressers (23.2% lower CL/F), equivalently expressers clear 1/0.768 = 1.30 times faster. Zhou 2026 pooled *1/*1 (4.4% of samples) with *1/*3 (35.6%) into the single expresser reference group; 60.0% of samples are *3/*3 non-expressers. Genotype was determined by PCR and Sanger sequencing of rs776746 and the distribution satisfies Hardy-Weinberg equilibrium (Table 3: chi-squared 0.395, p = 0.821).",
      source_name        = "CYP3A5"
    )
  )

  # Covariates that Zhou 2026 collected and screened but did NOT retain in the
  # population PK model. Most were carried forward only into the downstream
  # machine-learning feature set (Results, "Feature Selection"), which is not
  # represented in this model file; the stepwise PK covariate search retained
  # only CYP3A5 genotype and Wuzhi capsule on CL/F.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 71 kg [IQR 61-79] (Zhou 2026 Table 1). Screened but not retained on CL/F or V/F by the stepwise PK covariate search, and subsequently dropped from the machine-learning feature set for multicollinearity with BMI (Table S3: VIF 4.945)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Median 24.40 kg/m^2 [IQR 22.45-27.90] (Zhou 2026 Table 1). Retained in the machine-learning feature set in preference to body weight, but not retained in the population PK model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Median 35.05 g/L [IQR 29.00-40.20] (Zhou 2026 Table 1). A machine-learning feature only; not retained on CL/F or V/F. Hypoalbuminaemia is prominent in nephrotic syndrome and the Discussion argues it should raise the free fraction, but the effect did not reach the OFV criterion in the PK covariate search."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Source column CCR; median 109.61 mL/min [IQR 82.43-143.35] (Zhou 2026 Table 1). A machine-learning feature only; not retained in the population PK model."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Median 6.40 mmol/L [IQR 5.00-8.57] (Zhou 2026 Table 1). A machine-learning feature only; not retained in the population PK model."
    ),
    TBIL = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Median 10.00 umol/L [IQR 7.29-13.59] (Zhou 2026 Table 1). A machine-learning feature only; not retained in the population PK model."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 44 years [IQR 31-53] (Zhou 2026 Table 1). Screened but not retained in the population PK model and not selected into the machine-learning feature set."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "41.5% female (Zhou 2026 Table 1). Screened but not retained in the population PK model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 141L,
    n_studies      = 1L,
    age_range      = "18 years and older; median 44 years (IQR 31-53)",
    age_median     = "44 years",
    weight_range   = "median 71 kg (IQR 61-79)",
    weight_median  = "71 kg",
    sex_female_pct = 41.5,
    race_ethnicity = c(Han = 67.1, `Kazakh/Uyghur` = 27.1, Hui = 4.6, Other = 1.2),
    disease_state  = "Adult patients with nephrotic syndrome receiving oral immediate-release tacrolimus as immunosuppressive therapy, sampled at steady state during routine therapeutic drug monitoring",
    dose_range     = "Oral immediate-release tacrolimus capsules, initial total daily dose 0.05-0.1 mg/kg divided into a morning and an evening dose (q12h) on an empty stomach, then adjusted to the therapeutic drug monitoring result; observed total daily dose median 3.00 mg (IQR 3.00-4.00)",
    regions        = "Single centre, First Affiliated Hospital of Xinjiang Medical University, Urumqi, Xinjiang, People's Republic of China",
    notes          = "Retrospective cohort of 182 consecutive patients screened; 141 patients contributing 410 steady-state trough samples formed the internal dataset (enrolled January 2018 to December 2019) from which this model was estimated. A temporal external validation set (enrolled January 2020 to December 2023) contributed a further 12 patients and 41 samples. Table 1 percentages are per-sample (N = 410), not per-patient. Samples were drawn after at least 3 days of continuous dosing; trough target for adult nephrotic syndrome is 5-10 ng/mL, and 46.6% of samples were below, 49.3% within and 4.1% above that window. Median time on tacrolimus at sampling 59 days (IQR 19-139). CYP3A5 rs776746 genotype distribution *1/*1 4.4%, *1/*3 35.6%, *3/*3 60.0%. Estimated in NONMEM 7.3 by FOCE; 500-sample non-parametric bootstrap and 1000-dataset NPDE used for evaluation, with all bootstrap-vs-estimate deviations within 5%. Ethics approval K201912-07."
  )

  ini({
    # -------------------------------------------------------------------------
    # Structural parameters. Zhou 2026 Supplementary Table S1, 'Final model'
    # Estimate column. The structural equation is the displayed equation in
    # Results, "Population Pharmacokinetic Model":
    #   CL/F (L/h) = theta_CL/F * theta_WZC * theta_CYP3A5 * e^0.285
    #   V/F (L) = 519 ; ka (1/h) = 4.5
    #
    # Only steady-state trough concentrations were available, so the absorption
    # phase carries no information and ka was fixed. Zhou 2026 Results:
    # "The peak time for TAC is 1-2 hours. As all data in this study
    # represented trough concentrations following steady-state dosing, accurate
    # fitting of absorption parameters was unfeasible. Therefore, based on
    # reference, 18 Ka was fixed at 4.5 h-1." Table S1 lists "4.5(fixed)" in
    # the base, final and bootstrap columns alike.
    # -------------------------------------------------------------------------
    lka <- fixed(log(4.5)); label("Absorption rate constant ka (1/h)")                    # Zhou 2026 Table S1: Ka = 4.5 h^-1 (fixed); literature value, no absorption-phase data
    lvc <- log(519);        label("Apparent volume of distribution V/F (L)")              # Zhou 2026 Table S1 final model: V/F = 519 L (RSE 14.1%; bootstrap median 522.32, 95% CI 351.0-671.5)
    lcl <- log(20.9);       label("Apparent clearance CL/F in a CYP3A5 expresser not taking Wuzhi capsule (L/h)")  # Zhou 2026 Table S1 final model: theta_CL/F = 20.9 L/h (RSE 4.8%; bootstrap median 20.73, 95% CI 18.58-23.01)

    # -------------------------------------------------------------------------
    # Covariate effects on CL/F. Zhou 2026 reports these as fold-change
    # MULTIPLIERS, not as exponential coefficients, so they are applied in the
    # ratio^covariate form (the repository convention for a tabulated
    # fold-change; cf. Goel_2016_Sonidegib.R and Shu_2024_posaconazole.R).
    # Each factor collapses to exactly 1 at its reference category.
    #
    # The CYP3A5 multiplier is applied to (1 - CYP3A5_EXPR) because the paper
    # codes the *3/*3 non-expresser as the affected group while the canonical
    # covariate carries the expresser-equals-1 orientation; see
    # covariateData$CYP3A5_EXPR.
    # -------------------------------------------------------------------------
    e_conmed_wuzhi_cl <- 0.731; label("Multiplicative Wuzhi-capsule CL/F ratio (applied as ratio^CONMED_WUZHI)")            # Zhou 2026 Table S1 final model: theta_Wuzhi = 0.731 (RSE 7.5%; bootstrap median 0.760, 95% CI 0.69-0.85)
    e_cyp3a5_expr_cl  <- 0.768; label("Multiplicative CYP3A5*3/*3 non-expresser CL/F ratio (applied as ratio^(1 - CYP3A5_EXPR))")  # Zhou 2026 Table S1 final model: theta_CYP3A5 = 0.768 (RSE 5.2%; bootstrap median 0.770, 95% CI 0.69-0.86)

    # -------------------------------------------------------------------------
    # Inter-individual variability on CL/F only; no IIV was reported on V/F.
    #
    # SCALE DETERMINATION (this is the one number the source does not label
    # unambiguously). Table S1 reports the row "omega_CL/F" as a bare 0.285
    # (RSE 17.1%; bootstrap median 0.282, 95% CI 0.23-0.33) and the displayed
    # Results equation substitutes it straight into the exponential as
    # "e^0.285" where the eta belongs. Nothing states whether 0.285 is the
    # log-scale SD (omega) or the NONMEM variance (omega^2). The paper's own
    # observed trough distribution decides it, and it is not close:
    #
    #   Simulating this model at the published median dose (3 mg/day q12h) with
    #   the published genotype mix (40.0% expressers), Wuzhi prevalence (6.6%)
    #   and 20.4% proportional residual error gives a trough CV of
    #     45.4% reading 0.285 as the SD      (IIV = 29.1% CV)
    #     74.6% reading 0.285 as the variance (IIV = 57.4% CV)
    #   against an observed 5.39 +- 2.31 ng/mL, i.e. CV 42.9% (Table 1) -- and
    #   the observed figure additionally contains real dose variability (total
    #   daily dose IQR 3-4 mg) that the simulation holds fixed. The SD reading
    #   lands within 6% of the observed spread; the variance reading overshoots
    #   it by 74% and would require the TDM loop to have removed two-thirds of
    #   the variance, in a cohort where TDM control was demonstrably loose
    #   (46.6% of samples below the 5-10 ng/mL target, only 49.3% within it).
    #   The same test run against the published TDM-subgroup percentages puts
    #   the SD reading near the observed 4.1% above 10 ng/mL and the variance
    #   reading at several times that. Corroboration: the symbol omega (not
    #   omega^2) conventionally denotes the SD, and the closest sibling model
    #   in this library -- Xiang_2025_tacrolimus, same journal, same drug, same
    #   two covariates, same one-compartment trough-only structure -- reports
    #   IIV CL/F = 32.6% CV, against 29.1% here under the SD reading and 57.4%
    #   under the variance reading.
    #
    #   The bootstrap-CI-width test is inconclusive for this row and was not
    #   used: the relative width (0.33 - 0.23)/0.282 = 0.355 sits between the
    #   SD prediction 3.92*sqrt(1/(2*141)) = 0.233 and the variance prediction
    #   3.92*sqrt(2/141) = 0.467. The reported RSE of 17.1% is likewise
    #   uninformative here (nearer the variance asymptote, but omega RSEs from
    #   sparse trough-only data routinely run 2-3x above asymptotic).
    #
    #   Encoded as the variance 0.285^2 = 0.081225. See the vignette's
    #   "Assumptions and deviations" section, which reproduces this test.
    # -------------------------------------------------------------------------
    etalcl ~ 0.081225  # Zhou 2026 Table S1 final model: omega_CL/F = 0.285 read as the log-scale SD (RSE 17.1%; bootstrap median 0.282, 95% CI 0.23-0.33) -> variance 0.285^2

    # -------------------------------------------------------------------------
    # Residual error. Zhou 2026 Results: the base model was selected
    # "employing a proportional residual model". Table S1 reports the row as
    # "sigma_prop err(%)", i.e. already expressed as a percentage, so 20.4
    # corresponds to a proportional SD of 0.204.
    # -------------------------------------------------------------------------
    propSd <- 0.204; label("Proportional residual error for tacrolimus trough concentration")  # Zhou 2026 Table S1 final model: sigma_prop err = 20.4% (RSE 10.9%; bootstrap median 20.372, 95% CI 18.16-22.65)
  })

  model({
    # -----------------------------------------------------------------------
    # Individual PK parameters. Written in the multiplicative form the paper
    # prints, so each factor maps one-to-one onto a published term:
    #   CL/F = 20.9 * 0.731^WZC * 0.768^(*3/*3) * exp(eta_CL)
    # No covariate and no random effect were retained on V/F.
    # -----------------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl) *
      e_conmed_wuzhi_cl^CONMED_WUZHI *
      e_cyp3a5_expr_cl^(1 - CYP3A5_EXPR)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # -----------------------------------------------------------------------
    # Observation. The dose is in mg and V/F is in L, so central / vc is in
    # mg/L; tacrolimus troughs are reported in ng/mL and 1 mg/L = 1000 ng/mL.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
