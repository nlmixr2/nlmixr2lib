Knights_2015_aripiprazole <- function() {
  description <- "Two-compartment population PK model for oral aripiprazole in adult psychiatric patients (Knights 2015), with first-order absorption, linear-deviation weight (gated by WT < 115 kg) and age effects on apparent oral clearance, multiplicative CYP2D6 poor-metabolizer effect on CL/F, linear weight (gated by WT < 115 kg) and age effects with multiplicative female-sex effect on the peripheral volume, linear weight (gated by WT < 115 kg) effect with multiplicative female-sex effect on apparent inter-compartmental clearance, correlated inter-individual variability across Vc/F, Q/F, and Vp/F, independent IIV on ka and CL/F, and a proportional residual error."
  reference <- "Knights J, Rohatagi S. Development and application of an aggregate adherence metric derived from population pharmacokinetics to inform clinical trial enrichment. J Pharmacokinet Pharmacodyn. 2015;42(3):263-273. doi:10.1007/s10928-015-9414-4"
  vignette <- "Knights_2015_aripiprazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Knights 2015. Centred at 74.19 kg (the supplement NONMEM control stream centring value matches the population mean of the 448-subject 24-study pooled aripiprazole dataset). The weight effect on CL/F, Q/F, and Vp/F is gated by an indicator that is 1 only when WT < 115 kg and 0 otherwise (supplement variables WTA, WTPA, WTOKQ -- all three are identical 'WT-below-115' thresholds, plus a defensive non-missing guard WT > 0 that is not needed when WT is non-missing in modern datasets). The paper's Eq. 2 prose phrases the threshold as 'WT <= 115 kg' but the supplement control stream uses strict WT < 115; the strict-inequality form is used here as the more authoritative source.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Knights 2015. Centred at 32 years (Knights 2015 Eq. 2 and supplement NONMEM control stream). Linear-deviation effect on CL/F (per year) and on Vp/F (per year). Studied population age range 18-55 years (Knights 2015 Clinical data section, 47-patient validation cohort; the 24-study 448-subject popPK building cohort is described as having a wider but unspecified age distribution).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex, female indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "1 = female, 0 = male. The Knights 2015 source dataset's SEX column uses the opposite convention (SEX = 0 -> female, SEX = 1 -> male) confirmed by the supplement cross-tabulation of SEX vs. MMAS04 / MMAS06 in the clinical-trial 47-subject cohort: the supplement reports 31 SEX=1 subjects and 16 SEX=0 subjects, matching the paper's text 'enrolled 47 patients (31 male)'. The supplement control stream's intermediate variable SEXO = (SEX == 0) therefore equals 1 for females, which is identical to the nlmixr2lib canonical SEXF. To use a dataset that encodes SEX with 0 = male / 1 = female, pass the value unchanged as SEXF; to use a dataset that follows the Knights 2015 source convention (0 = female / 1 = male), pass SEXF = 1 - SEX.",
      source_name        = "SEX (0 = female, 1 = male in the source) recoded to SEXF = 1 - SEX"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, intermediate, or ultrarapid metabolizer)",
      notes              = "1 = subject is a CYP2D6 poor metabolizer (genotype encoding no functional CYP2D6 activity), 0 otherwise. Knights 2015 Eq. 2 reports the binary as `2D6PM` (PM = 1, non-PM = 0). The supplement NONMEM control stream stores the dataset column as CYP2D6EM (extensive-metabolizer indicator) and derives PM internally as PM = (CYP2D6EM == 0); the supplement explicitly treats CYP2D6EM = -99 (missing genotype) as not-extensive (i.e. PM = 1 if EM = 0 OR genotype is missing). In modern datasets without the -99 missing-value sentinel, set CYP2D6_PM = 1 only for genotypes encoding no functional CYP2D6 activity and CYP2D6_PM = 0 otherwise (including unknowns where defaulting to non-PM is the safer choice). CYP2D6 PMs have 47.8% lower apparent oral aripiprazole clearance than non-PMs in the Knights 2015 final popPK model.",
      source_name        = "2D6PM (paper text and Figure 1B); derived in the supplement control stream from CYP2D6EM as `PM = (CYP2D6EM == 0)`"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 448L,
    n_studies        = 24L,
    n_observations   = 13500L,
    age_range        = "Studied range not reported numerically for the 24-study popPK cohort. The 47-subject MMAS8 validation cohort spanned 18-55 years (Knights 2015 Clinical data).",
    weight_range     = "Studied range not reported numerically; the supplement covariate-effect equations and the WT < 115 kg gating indicator both imply that the cohort included subjects with WT both below and at-or-above 115 kg (otherwise the WTPA / WTOKQ indicators would never have differentiated subjects).",
    sex_female_pct   = NA_real_,
    sex_balance      = "Not reported for the 24-study popPK cohort. The 47-subject MMAS8 validation cohort was 16 female / 31 male (34.0% female).",
    race_ethnicity   = "Not reported in the Knights 2015 abstract, methods, or supplement.",
    disease_state    = "Pooled across the 24 clinical studies (referred to as 'the entire family of PK studies that were conducted during the approval process'); spans healthy volunteers and patients with psychiatric disorders treated with oral aripiprazole. The 47-subject MMAS8 validation cohort had a current diagnosis of bipolar 1 disorder (n = 15) or schizophrenia (n = 32) per DSM-IV-TR criteria.",
    dose_range       = "Oral aripiprazole, dose levels and dosing regimens spanning the family of 24 clinical studies. The 47-subject MMAS8 validation cohort had been on stable oral doses of 10, 15, 20, or 30 mg once daily for at least 2 weeks before sampling.",
    regions          = "Not specified (multi-study pooled dataset).",
    cyp2d6_pm_pct    = NA_real_,
    notes            = "Knights 2015 was designed to apply, rather than to develop, the underlying popPK model: the paper focuses on a novel 'reverse' application of the popPK model to compute an aggregate adherence metric (ADHMET) from steady-state sparse-sample plasma concentrations. The popPK model itself was built from 24 clinical studies, 448 individuals, and over 13,500 plasma aripiprazole observations submitted as part of the original aripiprazole new-drug-application package; the per-study composition is not separately tabulated in the paper. Parameter point estimates were transcribed from Figure 1B and the explicit Eq. 2 (CL/F equation); covariate-effect functional forms came from the supplement NONMEM control stream ($PK block)."
  )

  ini({
    # Structural PK -- Knights 2015 Figure 1B 'Abilify Model Parameter Estimates'
    # (verified line by line against the supplement NONMEM control stream
    # $THETA initial-value block). Apparent oral clearances (CL/F, Q/F) in
    # L/h; apparent volumes (Vc/F, Vp/F) in L; absorption rate constant in
    # 1/h. The reference subject is a 74.19-kg, 32-year-old, non-poor-CYP2D6-
    # metabolizer male: WT_indicator = 1, AGE-32 = 0, CYP2D6_PM = 0, SEXF = 0;
    # so all covariate-deviation factors reduce to 1.0 and the typical values
    # equal the structural thetas.
    lka <- log(0.515) ; label("Absorption rate constant ka (1/h)")                                                  # Knights 2015 Figure 1B KA = 0.515 1/h; supplement $THETA(1)
    lvc <- log(192)   ; label("Apparent central volume Vc/F at reference covariates (L)")                            # Knights 2015 Figure 1B VC = 192 L; supplement $THETA(2)
    lq  <- log(12.2)  ; label("Apparent inter-compartmental clearance Q/F at reference covariates (L/h)")            # Knights 2015 Figure 1B Q = 12.20 L/hr; supplement $THETA(3)
    lvp <- log(151)   ; label("Apparent peripheral volume Vp/F at reference covariates (L)")                         # Knights 2015 Figure 1B VP = 151 L; supplement $THETA(4)
    lcl <- log(3.88)  ; label("Apparent oral clearance CL/F at reference covariates (L/h)")                          # Knights 2015 Figure 1B CL = 3.88 L/hr; supplement $THETA(5)

    # Covariate effects -- the WT covariate enters as a linear-deviation
    # additive shift to the typical CL/F, Q/F, and Vp/F (in linear-scale
    # L/h or L), gated by an indicator that is 1 only when WT < 115 kg
    # (supplement WTPA / WTOKQ / WTA indicators are all identical
    # 'WT-below-115' thresholds). The AGE covariate enters as a linear-
    # deviation additive shift to the typical CL/F and Vp/F, always applied.
    # The CYP2D6_PM and SEXF covariates enter as proportional shifts on
    # the linear-scale typical value: TVCL *= (1 + e_2d6pm_cl * CYP2D6_PM);
    # TVQ *= (1 + e_sexf_q * SEXF); TVVP *= (1 + e_sexf_vp * SEXF).
    e_wt_vp     <- 4.07     ; label("WT linear-deviation slope on Vp/F (L/kg; gated by WT < 115 kg, centred at 74.19 kg)")    # Knights 2015 Figure 1B WT_VP = 4.07; supplement $THETA(6)
    e_age_vp    <- 0.915    ; label("AGE linear-deviation slope on Vp/F (L/year; centred at 32 years; always applied)")      # Knights 2015 Figure 1B AGE_VP = 0.915; supplement $THETA(7)
    e_wt_cl     <- 0.0251   ; label("WT linear-deviation slope on CL/F (L/h/kg; gated by WT < 115 kg, centred at 74.19 kg)") # Knights 2015 Eq. 2 WT coefficient 0.0251; Figure 1B WT_CL = 0.025 (Eq. 2 carries the more precise value); supplement $THETA(8)
    e_wt_q      <- 0.425    ; label("WT linear-deviation slope on Q/F (L/h/kg; gated by WT < 115 kg, centred at 74.19 kg)")  # Knights 2015 Figure 1B WT_Q = 0.425; supplement $THETA(9)
    e_age_cl    <- -0.0167  ; label("AGE linear-deviation slope on CL/F (L/h/year; centred at 32 years; always applied)")    # Knights 2015 Eq. 2 AGE coefficient -0.0167; Figure 1B AGE_CL = -0.017 (Eq. 2 carries the more precise value); supplement $THETA(10)
    e_2d6pm_cl  <- -0.478   ; label("CYP2D6 poor-metabolizer proportional shift on CL/F (fraction; -0.478 = 47.8% lower CL/F in PMs)") # Knights 2015 Eq. 2 and Figure 1B 2D6PM_CL = -0.478 (fold-decrease); supplement $THETA(11)
    e_sexf_q    <- 0.543    ; label("Female-sex proportional shift on Q/F (fraction; +0.543 = 54.3% higher Q/F in females)")  # Knights 2015 Figure 1B SEXF_Q = 0.543; supplement $THETA(12)
    e_sexf_vp   <- 0.341    ; label("Female-sex proportional shift on Vp/F (fraction; +0.341 = 34.1% higher Vp/F in females)") # Knights 2015 Figure 1B SEXF_VP = 0.341; supplement $THETA(13)

    # Inter-individual variability -- Knights 2015 Figure 1B 'ETA VARIANCE'
    # row reports omega^2 on the internal log scale directly; the 'THETA
    # IIV(%)' row is sqrt(omega^2) * 100 (the IIV expressed as a log-scale
    # standard deviation in percent). The supplement $OMEGA block carries
    # the more precise covariances; values below are transcribed verbatim
    # from the supplement.
    #
    # ETA structure (supplement $OMEGA blocks):
    #   ETA1 (KA): diagonal, omega^2 = 0.398
    #   ETA2-4 (VC, Q, VP): BLOCK(3) with the lower-triangle entries
    #     omega^2_VC                       = 4.93e-2
    #     cov(VC, Q)        omega^2_Q      = 3.33e-2  5.53e-2
    #     cov(VC, VP)  cov(Q, VP)  omega^2_VP = 4.59e-2  1.21e-2  7.67e-2
    #   ETA5 (CL): diagonal, omega^2 = 0.153
    etalka ~ 0.398                                              # Knights 2015 supplement $OMEGA ETA1 KA = 0.398 (Figure 1B ETA_KA = 0.398; IIV% = 63.09%)
    etalvc + etalq + etalvp ~ c(4.93e-2,
                                3.33e-2,  5.53e-2,
                                4.59e-2,  1.21e-2,  7.67e-2)    # Knights 2015 supplement $OMEGA BLOCK(3) ETA2-4 VC/Q/VP (Figure 1B ETA_VC = 0.049, ETA_Q = 0.055, ETA_VP = 0.077; cov(VC,Q) = 0.033, cov(VC,VP) = 0.046, cov(Q,VP) = 0.012)
    etalcl ~ 0.153                                              # Knights 2015 supplement $OMEGA ETA5 CL = 0.153 (Figure 1B ETA_CL = 0.153; IIV% = 39.12%)

    # Residual unexplained variability -- Knights 2015 Figure 1B reports
    # SIGMA = 0.0538 (variance on the proportional-error scale). The
    # supplement $SIGMA confirms (0.0538) and the $ERROR block writes
    # Y = IPRED * (1 + ERR(1)) -- i.e., a pure-proportional linear-scale
    # residual error. nlmixr2's ~ prop(propSd) parameter is the SD (not
    # the variance), so propSd = sqrt(0.0538) ~ 0.232.
    propSd <- sqrt(0.0538) ; label("Proportional residual error (fraction)")  # Knights 2015 Figure 1B SIGMA = 0.0538 (variance); supplement $SIGMA = 0.0538; propSd = sqrt(0.0538) = 0.232 (EPSshrink 5.91%)
  })

  model({
    # Indicator: WT effects on CL/F, Q/F, and Vp/F are gated by the
    # supplement's strict WT < 115 kg threshold (WTPA / WTOKQ / WTA;
    # paper text phrases the same threshold as 'WT <= 115 kg' but the
    # supplement uses the strict-inequality form, which is taken as
    # authoritative). At WT >= 115 kg the WT covariate effect is set
    # to zero on all three parameters -- producing a small step
    # discontinuity at the threshold that is faithful to the source
    # model.
    ind_wt_lt_115 <- (WT < 115)

    # Typical-value structural parameters (linear scale) with the full
    # Knights 2015 covariate equations:
    #   TVCL  = (CL  + WTPA  * e_wt_cl  * (WT - 74.19) + e_age_cl * (AGE - 32)) * (1 + e_2d6pm_cl * CYP2D6_PM)
    #   TVQ   = (Q   + WTOKQ * e_wt_q   * (WT - 74.19))                          * (1 + e_sexf_q  * SEXF)
    #   TVVP  = (VP  + WTPA  * e_wt_vp  * (WT - 74.19) + e_age_vp * (AGE - 32)) * (1 + e_sexf_vp * SEXF)
    #   TVVC  =  VC                                                                                                  (no covariates)
    #   TVKA  =  KA                                                                                                  (no covariates)
    tvcl <- (exp(lcl) + ind_wt_lt_115 * e_wt_cl * (WT - 74.19) + e_age_cl * (AGE - 32)) *
            (1 + e_2d6pm_cl * CYP2D6_PM)
    tvq  <- (exp(lq)  + ind_wt_lt_115 * e_wt_q  * (WT - 74.19)) *
            (1 + e_sexf_q * SEXF)
    tvvp <- (exp(lvp) + ind_wt_lt_115 * e_wt_vp * (WT - 74.19) + e_age_vp * (AGE - 32)) *
            (1 + e_sexf_vp * SEXF)
    tvvc <- exp(lvc)
    tvka <- exp(lka)

    # Individual PK parameters -- the supplement's MU_i = LOG(TVi);
    # PARAM_i = EXP(MU_i + ETA_i) is mathematically PARAM_i = TVi *
    # EXP(ETA_i). IIV thus applies multiplicatively to the linear-scale
    # covariate-adjusted typical value.
    ka <- tvka * exp(etalka)
    cl <- tvcl * exp(etalcl)
    vc <- tvvc * exp(etalvc)
    q  <- tvq  * exp(etalq)
    vp <- tvvp * exp(etalvp)

    # Micro-constants for the two-compartment oral-absorption ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption. Dose lands in
    # `depot`. Bioavailability is implicit (CL, Vc, Q, Vp are all 'apparent'
    # /F values -- the popPK model was developed from oral-only data and
    # cannot separate F from CL or volumes; see Knights 2015 Methods note
    # that the diagnostic VPC underpredicted the highest observed
    # concentrations 'which is not uncommon for oral medications without
    # intravenous data').
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Aripiprazole plasma concentration in ng/mL. Supplement $PK scaling:
    # S2 = VC / 1000 with units comment 'C_UNITS IN NG/ML, DOSE IN MG'.
    # That is equivalent to Cc (ng/mL) = (central amount, mg) / (Vc, L)
    # * 1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
