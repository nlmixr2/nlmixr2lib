Siccardi_2012_efavirenz <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for",
    "efavirenz in European HIV-1-positive adults and healthy volunteers,",
    "with CYP2B6 c.516G>T (rs3745274) genotype as the only retained",
    "covariate on apparent oral clearance (CL/F = 13.3 / 9.7 / 2.6 L/h for",
    "516GG / 516GT / 516TT). Two window-average exposure accumulator states",
    "reconstruct the paper's pharmacodynamic driver C8-16h (the mean plasma",
    "concentration between 8 and 16 h after a dose), which feeds two binary",
    "logistic exposure-response models: the probability of viral suppression",
    "and the probability of central-nervous-system side effects. The paper's",
    "companion Simcyp whole-body IVIVE/PBPK arm is not reproduced here (the",
    "platform system parameters are not published); only the NONMEM",
    "population PK model and the logistic PK/PD relationships are encoded.",
    sep = " "
  )
  reference <- paste(
    "Siccardi M, Almond L, Schipani A, Csajka C, Marzolini C, Wyen C,",
    "Brockmeyer NH, Boffito M, Owen A, Back D.",
    "Pharmacokinetic and pharmacodynamic analysis of efavirenz dose",
    "reduction using an in vitro-in vivo extrapolation model.",
    "Clin Pharmacol Ther. 2012;92(4):494-502. doi:10.1038/clpt.2012.61.",
    sep = " "
  )
  vignette <- "Siccardi_2012_efavirenz"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Doses are carried in mg and volumes in L, so `central / vc` is in mg/L;
  # `Cc` multiplies by 1000 to reach the ng/mL used for every concentration
  # in the paper. The unit is load-bearing rather than cosmetic: the two
  # logistic exposure-response models take log10(C8-16h) with C8-16h in
  # ng/mL, so a 1000-fold unit slip would move the linear predictor of the
  # viral-suppression model by 3 * 3.12 = 9.4 logits.

  covariateData <- list(
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = "0 (516GG homozygous wild-type)",
      notes              = paste(
        "Time-invariant germline genotype; the only covariate retained in the",
        "final population PK model (Siccardi 2012 Results, 'Population PK",
        "model': CL/F was significantly (P < 0.001) correlated with CYP2B6",
        "516G>T genotype and no demographic covariate was statistically",
        "significant). The paper's Methods writes the effect as a linear",
        "additive decomposition of the typical value,",
        "CL = CL0 + theta1 * GT + theta2 * TT, with GT and TT mutually",
        "exclusive indicators over the three genotypes and 516GG as the",
        "reference; the canonical count column reconstructs them via",
        "(SNP_CYP2B6_RS3745274_T_COUNT == 1) and (... == 2). The paper",
        "reports the three resulting typical values rather than theta1 and",
        "theta2, so the two shifts in ini() are the exact differences",
        "9.7 - 13.3 and 2.6 - 13.3 L/h. Cohort distribution (Siccardi 2012",
        "Results, 'Population PK model', n = 157): 41 subjects 516GT, 10",
        "subjects 516TT, therefore 106 subjects 516GG.",
        sep = " "
      ),
      source_name        = "CYP2B6 516 G>T genotype"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Siccardi 2012 Results, 'Population PK",
        "model': weight was one of five covariates entered into the stepwise",
        "backward elimination and 'no demographic covariates were",
        "statistically significant'. Median 70.5 kg (IQR 49-98 kg).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained (same backward-elimination step as WT).",
        "Median 41 years (IQR 22-67 years).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained (reported as 'gender'). 110 of 157",
        "subjects (70%) were men, so SEXF = 1 in 30% of the cohort.",
        sep = " "
      )
    ),
    RACE_BLACK = list(
      description = "Black-race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained (reported as 'ethnicity'). Siccardi 2012",
        "does not tabulate the ethnic composition of the population PK",
        "cohort; the IVIVE arm was run in a North European Caucasian virtual",
        "population, so no reference-category coding can be recovered and no",
        "coefficient is available.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    auc_8_16 = list(analyte = "efavirenz", units = NA_character_, specimen = "not applicable", verified = TRUE),
    t_8_16   = list(analyte = "efavirenz", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  # Bookkeeping accumulators that reconstruct the PD driver C8-16h -- the
  # numerator (concentration integrated over the 8-16 h post-dose window)
  # and the denominator (elapsed time inside that window). Same
  # construction as the Assmus_2025_benznidazole_mouse.R auc_total /
  # t_above_mic pair, except that the window is relative to the most recent
  # dose via tad() rather than to absolute time.
  paper_specific_compartments <- c("auc_8_16", "t_8_16")

  population <- list(
    species        = "human",
    n_subjects     = 157L,
    n_studies      = 2L,
    age_range      = "IQR 22-67 years (population PK cohort)",
    age_median     = "41 years (population PK cohort)",
    weight_range   = "IQR 49-98 kg (population PK cohort)",
    weight_median  = "70.5 kg (population PK cohort)",
    sex_female_pct = 30,
    disease_state  = paste(
      "HIV-1 infection on an efavirenz-based regimen combined with two",
      "nucleoside reverse-transcriptase inhibitors, at steady state, plus",
      "nine healthy volunteers. Patients on a boosted protease inhibitor or",
      "on any other drug known to interact with efavirenz were excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "600 mg once daily (the observed standard regimen); 400 mg and 200 mg",
      "once daily were simulated dose reductions, not observed.",
      sep = " "
    ),
    regions        = "United Kingdom (nine healthy volunteers, Royal Free NHS Trust, London) and Germany (148 HIV-positive patients, KompNet cohort)",
    genotype_distribution = "516GG 106 (67.5%), 516GT 41 (26.1%), 516TT 10 (6.4%)",
    n_observations = "202 plasma efavirenz samples spanning 530-26,020 ng/mL; six samples per healthy volunteer and one random sample per HIV-positive patient",
    pd_cohorts     = paste(
      "The two logistic exposure-response models were fitted to separate",
      "cohorts drawn from Csajka 2003 and Marzolini 2001, not to the",
      "population PK cohort above. Viral suppression: 93 patients (14 with",
      "therapeutic failure), 65 men (68%), median age 38 years (IQR 33-44),",
      "median weight 66 kg (IQR 60-75). CNS side effects: 121 patients (30",
      "with CNS side effects), 88 men (73%), median age 39 years",
      "(IQR 35-45), median weight 66 kg (IQR 57-75). Both cohorts were",
      "co-treated with the older nucleoside reverse-transcriptase",
      "inhibitors zidovudine, stavudine and didanosine.",
      sep = " "
    ),
    notes          = paste(
      "Baseline demographics are in Siccardi 2012 Results, 'Population PK",
      "model' (PK cohort) and 'PD of dose reduction' (the two PD cohorts).",
      "The paper's companion IVIVE arm simulated 500 virtual North European",
      "Caucasian subjects per genotype (20-50 years old, 0.5 proportion",
      "female) in Simcyp 10.1; those virtual-population characteristics",
      "belong to the PBPK arm, which this model does not reproduce.",
      sep = " "
    )
  )

  ini({
    # ---- Structural population PK parameters (final model) -------------
    # Siccardi 2012 Table 2 reports the final population PK model beside
    # the IVIVE predictions. The base model (Results, 'Population PK
    # model') had CL/F 11.6 L/h, V/F 317 L and ka 0.54 1/h; the final
    # model adds the CYP2B6 516G>T effect on CL/F only.
    lka <- log(0.36)
    label("First-order absorption rate constant (ka, 1/h)")  # Siccardi 2012 Table 2, Population PK column: ka = 0.36 1/h (RSE 76%). Very imprecise (one random sample per HIV-positive patient); the bootstrap 90% CI is 0.4-1.1 1/h and the base model gave 0.54 1/h.

    lcl <- log(13.3)
    label("Typical apparent oral clearance CL/F for CYP2B6 516GG homozygotes (L/h)")  # Siccardi 2012 Table 2, Population PK column: CL/F 516GG = 13.3 L/h (RSE 7%); same value in Results, 'Population PK model'. Bootstrap 90% CI 10.6-14.1 L/h.

    lvc <- log(4.3 * 70.5)
    label("Typical apparent volume of distribution V/F (L)")  # Siccardi 2012 Table 2, Population PK column: Vd = 4.3 (RSE 10%). The column header reads '(l)' but the value is per kg: 4.3 L/kg * 70.5 kg (the cohort median weight, Results 'Population PK model') = 303 L, which sits inside the paper's own bootstrap 90% CI for the volume (260.5-339.5 L) and beside the base-model 317 L, whereas 4.3 L absolute would give a 0.22 h half-life that cannot produce the reported once-daily troughs. The matching IVIVE entry (8.1) is a Simcyp Vss, natively reported in L/kg, and the paper's 'less than twofold difference' claim (8.1 / 4.3 = 1.9) only holds if both columns share one unit.

    lfdepot <- fixed(log(1))
    label("Bioavailability of the depot compartment (fraction)")  # Structural anchor: Siccardi 2012 reports apparent parameters CL/F and V/F, so F is not separately identifiable and is held at unity.

    # ---- CYP2B6 516G>T effect on CL/F ---------------------------------
    # Siccardi 2012 Methods, 'Population PK model at standard regimens':
    #   CL = CL0 + theta1 * GT + theta2 * TT
    # with GT / TT mutually exclusive indicators and 516GG the reference.
    # The paper prints the three resulting typical values instead of
    # theta1 and theta2, so each shift below is an exact difference of
    # printed values.
    e_516gt_cl <- 9.7 - 13.3
    label("Linear-additive shift on CL/F for CYP2B6 516GT heterozygotes (L/h)")  # Siccardi 2012 Table 2, Population PK column: CL/F 516GT = 9.7 L/h (RSE 12%) vs 13.3 L/h for 516GG, i.e. -3.6 L/h (27.1% lower); the Discussion reports the corresponding IVIVE decrement as 26%.

    e_516tt_cl <- 2.6 - 13.3
    label("Linear-additive shift on CL/F for CYP2B6 516TT homozygotes (L/h)")  # Siccardi 2012 Table 2, Population PK column: CL/F 516TT = 2.6 L/h (RSE 41%) vs 13.3 L/h for 516GG, i.e. -10.7 L/h (80.5% lower). Only 10 subjects were 516TT; the Discussion states the IVIVE arm underestimated this effect (IVIVE decrement 54%) and flags the disagreement.

    # ---- Pharmacodynamics: two binary logistic exposure-response models
    # Siccardi 2012 Equation 1 (page 497):
    #   p = 1 / (1 + e^-(A + B * log10(C8-16h)))
    # i.e. a natural-logit link with log10(C8-16h) in ng/mL as the single
    # predictor. The natural-log link is corroborated by the paper's own
    # reported odds ratios per log10 unit: exp(3.12) = 22.65 against the
    # reported OR 22.6 for viral suppression, and exp(1.68) = 5.37 against
    # the reported OR 5.4 for CNS side effects. The sentence following
    # Equation 1 ("output being log10(p/(1-p)) = A + B x") contradicts
    # Equation 1 and the odds ratios, and is not used -- see the vignette
    # Assumptions and deviations section.
    logite0_supp <- -8.38
    label("Logit of the probability of viral suppression at C8-16h = 1 ng/mL (unitless)")  # Siccardi 2012 Results, 'PD of dose reduction' (below Equation 1): 'For viral suppression, A = -8.38 ... with a SE of 100% for A'.

    e_c816_supp <- 3.12
    label("Log-odds of viral suppression per 10-fold increase in C8-16h (unitless)")  # Siccardi 2012 Results, 'PD of dose reduction': 'B = 3.12 ... 42% for B'; corroborated by the reported odds ratio for log10 C8-16h, OR = 22.6 (95% CI 1.37-372, P = 0.029), since exp(3.12) = 22.65.

    logite0_cns <- -6.65
    label("Logit of the probability of CNS side effects at C8-16h = 1 ng/mL (unitless)")  # Siccardi 2012 Results, 'PD of dose reduction': 'For CNS side effects, A = -6.65 ... with an SE of 43% for A'.

    e_c816_cns <- 1.68
    label("Log-odds of CNS side effects per 10-fold increase in C8-16h (unitless)")  # Siccardi 2012 Results, 'PD of dose reduction': 'B = 1.68 ... 51% for B'; corroborated by the reported odds ratio for log10 C8-16h, OR = 5.4 (95% CI 1.01-28.7, P = 0.05), since exp(1.68) = 5.37.

    # ---- Inter-individual variability ---------------------------------
    # Siccardi 2012 Results, 'Population PK model': interindividual random
    # effects were exponential, theta_i = theta * exp(eta_i), and the base
    # model carried 'an interindividual variability (expressed as
    # percentage of the coefficient of variation) of 57.4% on CL/F'. For a
    # lognormal parameter the stated coefficient of variation gives
    #   omega^2 = log(CV^2 + 1) = log(0.574^2 + 1) = 0.2847849
    # The conversion is written out rather than transcribed as a decimal so
    # the published 57.4% CV stays visible at the assignment site.
    # No IIV is reported on V/F or ka, and none is added here.
    etalcl ~ log(0.574^2 + 1)  # Siccardi 2012 Results, 'Population PK model': 57.4% CV on CL/F in the base model. The final model's residual IIV is not printed -- the paper states only that the genotype covariate 'explained 12.6% of the interindividual variability of CL/F' without saying whether that share is of the variance or of the CV -- so the printed base-model value is carried forward unreduced. See the vignette Assumptions and deviations section.

    # ---- Residual error ------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; 0 -- magnitude not reported in the source)")  # Siccardi 2012 Results, 'Population PK model': 'residual variability was best described by a proportional structure; the inclusion of an additive structure did not improve the model'. The structure is stated but no magnitude is reported anywhere in the paper, so the proportional SD is held at zero rather than invented.
  })

  model({
    # 1. Individual parameters. The genotype effect is a linear-additive
    #    shift on the typical value, and the exponential IIV multiplies the
    #    covariate-adjusted typical value (Siccardi 2012 Methods:
    #    theta_i = theta * exp(eta_i) with theta the population value).
    het_516 <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    hom_516 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)

    tvcl <- exp(lcl) +
      e_516gt_cl * het_516 +
      e_516tt_cl * hom_516

    ka <- exp(lka)
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. One-compartment disposition with first-order absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Bioavailability
    f(depot) <- exp(lfdepot)

    # 5. Observation. Amounts are in mg and vc in L, so the 1000 converts
    #    mg/L to the ng/mL used throughout the paper.
    Cc <- 1000 * central / vc

    # 6. Reconstruction of the PD driver C8-16h -- "mean EFV concentrations
    #    between 8 and 16 h after the dose", Siccardi 2012 Methods,
    #    'PK/PD analysis'. `auc_8_16` integrates concentration only inside
    #    that post-dose window and `t_8_16` accumulates the elapsed time
    #    inside it, so their ratio is the window mean. Because both states
    #    run from the start of the simulation, `c816` is the running mean
    #    over every elapsed window and converges to the steady-state value;
    #    the exact value for one steady-state interval is obtained either
    #    by restarting the two states at the top of that interval or by
    #    differencing them across it (both shown in the vignette). The
    #    max() guards keep the ratio and its log10 finite before the first
    #    window opens, where t_8_16 is still zero.
    inwin_8_16 <- (tad() >= 8) * (tad() <= 16)
    d/dt(auc_8_16) <- Cc * inwin_8_16
    d/dt(t_8_16)   <- inwin_8_16
    c816 <- auc_8_16 / max(t_8_16, 1e-12)

    # 7. Exposure-response (Siccardi 2012 Equation 1). Two independent
    #    binary logistic regressions on log10(C8-16h).
    lc816 <- log10(max(c816, 1e-12))
    psupp <- 1 / (1 + exp(-(logite0_supp + e_c816_supp * lc816)))
    pcns  <- 1 / (1 + exp(-(logite0_cns  + e_c816_cns  * lc816)))

    # 8. Residual error: proportional structure, magnitude unpublished.
    Cc ~ prop(propSd)
  })
}
