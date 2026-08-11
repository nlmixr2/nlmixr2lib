Yang_2024_dabigatran <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for total dabigatran after a single oral 150 mg dose of dabigatran etexilate in healthy Chinese adults, with a high-fat-meal (postprandial) effect on the absorption rate constant, the absorption lag time and apparent clearance, and an ABCB1 rs4148738 heterozygote (CT) effect on the apparent central volume of distribution (Yang 2024)."
  reference <- "Yang Z, Tan WR, Li Q, Wang Y, Liu S, Chen L, Zhou Y, Zeng C, Zeng Y, Xiong Y, Zhang Q, Li N, Du P, Liu L, Chen J, He Y. Population pharmacokinetic study of the effect of polymorphisms in the ABCB1 and CES1 genes on the pharmacokinetics of dabigatran. Front Pharmacol. 2024;15:1454612. doi:10.3389/fphar.2024.1454612"
  vignette <- "Yang_2024_dabigatran"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Yang 2024 Supplementary Material
  # ($MODEL block: COMP1 = DEFDOSE, COMP2 = CENTRAL, COMP3 = PERIPHER) and
  # Methods 2.3 (analyte = total dabigatran in plasma by HPLC-MS/MS).
  compartmentData <- list(
    depot       = list(analyte = "dabigatran", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dabigatran", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dabigatran", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FED_HIGHFAT = list(
      description        = "High-fat-meal (postprandial) dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted administration)",
      notes              = "Time-fixed per subject: Yang 2024 used a parallel design in which each subject was dosed once, either fasted (n = 61) or within 30 min after a standard high-fat meal (n = 62); the paper does not report the caloric or fat content of the meal beyond 'standard high-fat meal' (Methods 2.2). The source NONMEM data column was named `YS` (Supplementary Material $INPUT and the `KAYS` / `CLYS` / `ALAG1YS` covariate blocks); renamed to the canonical `FED_HIGHFAT` per inst/references/covariate-columns.md. Enters ka, CL/F and the absorption lag time as a linear proportional deviation, (1 + theta * FED_HIGHFAT).",
      source_name        = "YS"
    ),
    SNP_ABCB1_RS4148738_HET = list(
      description        = "ABCB1 rs4148738 heterozygote (CT) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-CT: the union of CC homozygotes, n = 20, and TT homozygotes, n = 31)",
      notes              = "Time-fixed per subject (germline genotype). 99 of the 123 subjects were genotyped (Results 3.2, Table 2): CC 20, CT 48, TT 31; minor (C) allele frequency 44.44%. The source NONMEM dataset carried three dummy indicators, `ABCB1C` / `ABCB1CT` / `ABCB1T`, and only `ABCB1CT` survived stepwise covariate modelling (Supplementary Material $INPUT and the `V2ABCB1CT` block), so the reference group is the union of the two homozygous strata rather than the wild-type homozygote alone. Enters V2/F as a linear proportional deviation, (1 + theta * SNP_ABCB1_RS4148738_HET). For subjects who were not genotyped, set to 0 to obtain the non-CT typical value.",
      source_name        = "ABCB1CT"
    )
  )

  # Covariates screened by Yang 2024 but not retained in the final model
  # (Methods 2.6 lists the tested set; Results 3.5 and Discussion report that
  # only food intake and ABCB1 rs4148738 survived forward selection at
  # dOFV 3.84 and backward elimination at dOFV 10.83). Documented here for
  # provenance; no point estimates are reported for them anywhere in the paper.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Methods 2.6) but not retained. Range 18-43 years; the Discussion attributes the null result to the narrow healthy-volunteer age range."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (Methods 2.6) but not retained. Range 45.2-82.0 kg."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened (Methods 2.6) but not retained. Range 145.5-182.5 cm."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened (Methods 2.6) but not retained. Range 19.2-25.9 kg/m^2; the Discussion notes the range was too narrow to resolve an effect."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Methods 2.6) but not retained. 29 of 123 subjects (23.6%) were female; the Discussion states 'the smaller number of female subjects precluded a thorough analysis of gender effects'. Non-compartmental analysis (Table 3) showed mean AUC about 8% higher in women."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 rs1045642 (c.3435C>T) variant indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Genotyped (Table 2: AA 15, AG 47, GG 36; minor A allele frequency 39.29%) and screened, but not retained. Results 3.3: no significant effect on Cmax or AUC."
    ),
    SNP_CES1_RS2244613 = list(
      description = "CES1 rs2244613 variant indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Genotyped (Table 2: GG 40, GT 45, TT 14; minor T allele frequency 36.87%) and screened, but not retained. TT carriers had numerically higher mean Cmax and AUC in the non-compartmental analysis (Table 3), but the effect did not survive stepwise covariate modelling."
    ),
    SNP_CES1_RS8192935 = list(
      description = "CES1 rs8192935 variant indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Genotyped (Table 2: AA 49, AG 47, GG 3; minor G allele frequency 26.77%) and screened, but not retained. The three GG carriers had markedly higher mean Cmax and AUC (Table 3), but the Discussion states the sample size was too small to support a modelled effect."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 123L,
    n_studies      = 1L,
    n_observations = 1926L,
    age_range      = "18-43 years",
    age_median     = "25 years in both arms (mean 25.20 fasting, 24.94 postprandial)",
    weight_range   = "45.2-82.0 kg",
    weight_median  = "58.8 kg fasting arm, 59.7 kg postprandial arm",
    height_range   = "145.5-182.5 cm",
    bmi_range      = "19.2-25.9 kg/m^2",
    sex_female_pct = 23.6,
    race_ethnicity = "All subjects Chinese; the paper reports no sub-ethnicity breakdown (Methods 2.1, single centre in Guiyang, China).",
    disease_state  = "Healthy volunteers. Data came from the reference-formulation arm of a dabigatran etexilate bioequivalence study; renal impairment and co-medication were excluded by the inclusion criteria, which the Discussion notes removes the two covariates that dominate dabigatran PK in atrial-fibrillation patients.",
    dose_range     = "Single oral 150 mg dabigatran etexilate (Pradaxa) with 240 mL warm water, either fasted (n = 61) or within 30 min after a standard high-fat meal (n = 62). Sampling at 0 h predose and 0.25, 0.5, 0.75, 1, 1.33, 1.67, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 12, 24, 36 and 48 h (19 points per subject).",
    regions        = "China (Clinical Trials Center, Affiliated Hospital of Guizhou Medical University, Guiyang).",
    genotyped_n    = 99L,
    notes          = "Baseline demographics from Yang 2024 Table 1; genotype distributions from Table 2; non-compartmental exposure summaries from Table 3. Registered as NCT06387407; ethics approval 2024037K. 99 of the 123 subjects were genotyped for the four candidate SNPs, so simulations that exercise SNP_ABCB1_RS4148738_HET should reproduce the genotyped subset (CC 20 / CT 48 / TT 31) rather than the full 123."
  )

  ini({
    # ---- Structural parameters ------------------------------------------
    # Yang 2024 Table 4, "Final model / Estimate (RSE%)" column. All apparent
    # (CL/F, V/F): the source model carries no bioavailability parameter
    # (Supplementary Material $PK has no F1 statement), so F is implicitly 1.
    # The typical values below are the fasting, ABCB1 rs4148738 non-CT values;
    # the covariate terms in model() shift them to the other strata.
    lka   <- log(0.38);  label("First-order absorption rate constant, fasting (ka, 1/h)")                              # Yang 2024 Table 4 KA = 0.38 (RSE 16.2%, 95% CI 0.26-0.50)
    ltlag <- log(0.48);  label("Absorption lag time, fasting (ALAG, h)")                                               # Yang 2024 Table 4 ALAG = 0.48 (RSE 5%, 95% CI 0.44-0.53)
    lcl   <- log(122);   label("Apparent clearance, fasting (CL/F, L/h)")                                              # Yang 2024 Table 4 CL = 122 (RSE 5.6%, 95% CI 108.61-135.39)
    lvc   <- log(245);   label("Apparent central volume of distribution, ABCB1 rs4148738 non-CT (V2/F, L)")            # Yang 2024 Table 4 V2 = 245 (RSE 14.1%, 95% CI 177.18-312.82)
    lq    <- log(70.80); label("Apparent inter-compartmental clearance (Q/F, L/h)")                                    # Yang 2024 Table 4 Q = 70.80 (RSE 7.9%, 95% CI 59.80-81.80)
    lvp   <- log(756);   label("Apparent peripheral volume of distribution (V3/F, L)")                                 # Yang 2024 Table 4 V3 = 756 (RSE 10.6%, 95% CI 599.20-912.80)

    # ---- Covariate effects ----------------------------------------------
    # Linear proportional deviations, theta_i = theta_tv * (1 + theta_cov * COV).
    # The form is fixed by the Yang 2024 Supplementary Material control stream,
    # which writes each block as e.g.
    #   IF(YS.EQ.0) CLYS = 1
    #   IF(YS.EQ.1) CLYS = ( 1 + THETA(10) )
    #   TVCL = CLCOV*TVCL
    # so the tabulated coefficients are FRACTIONS, not percentages. The
    # Abstract and Results phrase them as "increased ALAG and CL by 2.65% and
    # 0.51%", which is an author slip: a 2.65% lag-time increase (0.48 ->
    # 0.49 h) cannot produce the observed 2.05 h postprandial Tmax delay
    # (Table 3), whereas 1 + 2.65 = 3.65-fold (0.48 -> 1.75 h) reproduces it.
    e_fed_highfat_ka   <- -0.24; label("High-fat meal effect on ka (fractional deviation)")                            # Yang 2024 Table 4 "Postprandial on KA" = -0.24 (RSE 30%, 95% CI -0.39 to -0.10)
    e_fed_highfat_tlag <-  2.65; label("High-fat meal effect on absorption lag time (fractional deviation)")           # Yang 2024 Table 4 "Postprandial on ALAG" = 2.65 (RSE 8.7%, 95% CI 2.20-3.10)
    e_fed_highfat_cl   <-  0.51; label("High-fat meal effect on CL/F (fractional deviation)")                          # Yang 2024 Table 4 "Postprandial on CL" = 0.51 (RSE 25.5%, 95% CI 0.26-0.77)
    e_abcb1_het_vc     <-  0.38; label("ABCB1 rs4148738 CT-heterozygote effect on V2/F (fractional deviation)")        # Yang 2024 Table 4 "rs4148738 on V2" = 0.38 (RSE 45.5%, 95% CI 0.04-0.71)

    # ---- Inter-individual variability ------------------------------------
    # Yang 2024 Table 4 "Inter-individual variability (IIV)" block reports
    # 0.31 (CL), 0.27 (V2) and 0.19 (ALAG) with no stated scale. These are read
    # here as omega on the STANDARD-DEVIATION scale (equivalently, apparent
    # CV), so the internal variances are omega^2. Rationale, from the paper's
    # own data: for this model AUCinf = Dose / (CL/F) exactly, so the observed
    # between-subject CV of AUCinf is an upper bound on the CL/F IIV CV.
    # Table 3 reports AUC0-inf CV of 33.5% (fasting) and 29.3% (postprandial),
    # and Cmax CV of 35.9% / 32.2%. Reading 0.31 as an SD gives 31.8% CV for
    # CL/F -- a near-exact match; reading it as a variance gives 60.3% CV,
    # roughly double the observed spread and therefore arithmetically
    # incompatible with Table 3. The same comparison rules out the variance
    # reading for V2 (27.5% vs 55.7% against an observed Cmax CV of 35.9%).
    # See the vignette "Assumptions and deviations" section, which also records
    # that the Supplementary Material $OMEGA block holds larger values
    # (0.2584 / 0.4350 / 0.1392) -- those are the run's INITIAL estimates from
    # an earlier fit (the structural $THETA initials likewise differ from the
    # Table 4 finals: KA 0.273 vs 0.38, CL 160.3 vs 122, ALAG1 1.056 vs 0.48),
    # not the final estimates, so they do not fix the scale of Table 4.
    etalcl   ~ 0.0961  # Yang 2024 Table 4 IIV CL   = 0.31 -> omega^2 = 0.31^2 (shrinkage 7.2%)
    etalvc   ~ 0.0729  # Yang 2024 Table 4 IIV V2   = 0.27 -> omega^2 = 0.27^2 (shrinkage 0.1%)
    etaltlag ~ 0.0361  # Yang 2024 Table 4 IIV ALAG = 0.19 -> omega^2 = 0.19^2 (shrinkage 14.6%)

    # ---- Residual error ---------------------------------------------------
    # Combined additive + proportional. The Supplementary Material $ERROR block
    #   W = SQRT(THETA(6)**2*IPRED**2 + THETA(7)**2); Y = IPRED + W*ERR(1)
    # with $SIGMA 1 FIX makes THETA(6) and THETA(7) residual standard
    # deviations, matching nlmixr2's add() + prop() parameterisation directly.
    # Table 4 labels THETA(6) "Exponential error"; Methods Equation 2 and the
    # $ERROR block both show it multiplying IPRED, i.e. it is proportional.
    propSd <- 0.08; label("Proportional residual error (fraction)")                                                    # Yang 2024 Table 4 "Exponential error" = 0.08 (RSE 24.9%, 95% CI 0.04-0.11)
    addSd  <- 6.31; label("Additive residual error (ng/mL)")                                                           # Yang 2024 Table 4 "Additive error" = 6.31 (RSE 14.7%, 95% CI 4.50-8.12)
  })
  model({
    # Individual parameters. The covariate term multiplies the typical value
    # and the exponential random effect is applied on top, exactly as in the
    # Yang 2024 control stream (TVCL = CLCOV*TVCL; CL = TVCL*EXP(ETA(1))).
    # ka, Q and V3 carry no IIV in the source model.
    ka   <- exp(lka)                * (1 + e_fed_highfat_ka   * FED_HIGHFAT)
    tlag <- exp(ltlag + etaltlag)   * (1 + e_fed_highfat_tlag * FED_HIGHFAT)
    cl   <- exp(lcl   + etalcl)     * (1 + e_fed_highfat_cl   * FED_HIGHFAT)
    vc   <- exp(lvc   + etalvc)     * (1 + e_abcb1_het_vc     * SNP_ABCB1_RS4148738_HET)
    q    <- exp(lq)
    vp   <- exp(lvp)

    # Micro-constants; Yang 2024 Supplementary Material $DES writes
    # K = CL/V2, K23 = Q/V2, K32 = Q/V3.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag on the dosing compartment (NONMEM ALAG1 on COMP1, the
    # DEFDOSE compartment).
    alag(depot) <- tlag

    # Dose in mg and vc in L give central/vc in mg/L = ug/mL; Yang 2024 reports
    # total dabigatran in ng/mL, so rescale by 1000. This reproduces the
    # source model's S2 = V2/1000 scaling factor.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
