vanHasselt_2014_cefazolin_semiphysiological <- function() {
  description <- "Two-compartment population PK model for free and total cefazolin in pregnant women (semiphysiological gestational covariate model). Clearance is the sum of a non-renal arm and a renal arm scaled by the gestational rise in creatinine clearance, which is generated inline by a fixed hyperbolic CrCL trajectory in gestational age rather than fitted to the PK data; disposition parameters are referenced to the UNBOUND concentration and total cefazolin is reconstructed algebraically as Cunbound/fu."
  reference <- paste(
    "van Hasselt JGC, Allegaert K, van Calsteren K, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Semiphysiological versus empirical modelling of the population pharmacokinetics of free and total cefazolin during pregnancy.",
    "Biomed Res Int. 2014;2014:897216. doi:10.1155/2014/897216.",
    "Corrigendum: Biomed Res Int. 2015;2015:124035. doi:10.1155/2015/124035",
    "(corrects the Table 2 covariate-equation footnotes, which are switched in the original,",
    "and clarifies that CL, Vc, Vp and Q were fitted on free cefazolin so CLtotal = CLfree * fu).",
    sep = " "
  )
  vignette <- "vanHasselt_2014_cefazolin"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    EGA = list(
      description = "Maternal estimated gestational age at the time of the observation",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives clearance INDIRECTLY, through the inline creatinine-clearance trajectory",
        "of Methods eq 4, whose value is then normalised by its own prepregnancy anchor",
        "to give the CrCL(t)/CrCL(0) ratio of Methods eq 6.",
        "Unlike the empirical sibling model, this form collapses exactly to the",
        "non-pregnant clearance cl_nonren + cl_renal at EGA = 0, because the ratio is 1",
        "there by construction. That is precisely the extrapolation property the paper",
        "argues for in its Discussion.",
        "Observed range 17-40 weeks, median 33 weeks (Table 1).",
        "A fixed EGA of 40 weeks was assigned to the term-pregnancy caesarean cohort (Methods 2.1).",
        "NO CrCL DATA COLUMN IS REQUIRED. The model consumes only the RATIO of the",
        "trajectory to its own EGA = 0 value, so the mL/min scale cancels and the",
        "trajectory is generated from EGA alone."
      ),
      source_name = "GA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 3L,
    n_observations = 187L,
    age_range = "20-42 years",
    age_median = "31 years",
    weight_range = "54-99 kg",
    weight_median = "72 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    disease_state = "pregnant women undergoing in utero surgical intervention, elective caesarean delivery, or fetal intervention; cefazolin given as surgical prophylaxis",
    dose_range = "1 g or 2 g intravenously; 2 g every 8 h for 2 days in the prospective cohort, single 1 g or 2 g bolus in the two literature cohorts",
    ga_range = "17-40 weeks (median 33)",
    regions = "Belgium (University Hospitals Leuven) plus two previously published cohorts",
    renal_function = "serum creatinine median 0.64 mg/dL (range 0.33-0.88); creatinine clearance computed by Cockcroft-Gault using body weight",
    notes = paste(
      "Pooled from one prospective study and two published studies with individual-level data (Methods 2.1).",
      "Prospective cohort: 41 pregnant women, 153 cefazolin observations, median GA 25 weeks (range 17-34),",
      "2 g every 8 h for 2 days during in utero surgery, free cefazolin available for 84% of observations.",
      "Fiore Mitchell et al.: 24 plasma samples at term caesarean delivery after 1 g i.v. bolus, mean sampling time 1.85 h,",
      "GA fixed at 40 weeks for all patients in that study.",
      "Brown et al.: 10 fetal interventions in 7 women, single 2 g i.v. bolus, mean GA 27 weeks, mean sampling time 0.5 h.",
      "Demographics in Table 1. Substantial missing data imputed by pooled medians (Methods 2.5);",
      "in this semiphysiological arm specifically, missing CrCL values were imputed by the typical",
      "change in CrCL predicted by the embedded trajectory rather than by a pooled median.",
      "TWO-STAGE FIT. The CrCL trajectory was developed separately (Methods 2.4.1) and its individual",
      "empirical Bayes estimates were then carried into the PK model as a known individual covariate",
      "time course (Methods 2.4.2), so only the PK parameters below were estimated against the",
      "cefazolin concentrations."
    )
  )

  ini({
    # Structural parameters. Source: van Hasselt 2014 Table 2, column
    # 'Semiphysiological Model CL ~ CrCL'. CL and Q are reported in L/min and
    # the volumes in L, so the model time unit is minutes and the values are
    # carried unconverted.
    #
    # CLEARANCE DECOMPOSITION. Methods eq 6, verbatim from the JATS MathML:
    #     CL_i = theta_CL0 + theta_CLpreg * (CrCL_i(t) / CrCL_i0)
    # and Methods 2.4.2 names the pair directly: 'theta_CL0 represents
    # non-CrCL-related clearance, and theta_CLpreg represents GFR-related
    # clearance, which is multiplied by the normalized change in CrCL'. That
    # sentence is the direct warrant for the registered renal / non-renal
    # clearance-arm canonicals used here. Neither arm is the total clearance:
    # the total is their sum.
    lcl_nonren <- log(0.142)
    label("Non-CrCL-related (non-renal) clearance arm, referenced to unbound cefazolin (L/min)") # Table 2, row 'Clearance' theta_CL0, semiphysiological column = 0.142 L/min (RSE 44%)
    lcl_renal <- log(0.212)
    label("GFR-related (renal) clearance arm at the prepregnancy CrCL anchor, referenced to unbound cefazolin (L/min)") # Table 2, row 'Gestation effect on clearance' theta_CLPreg, semiphysiological column = 0.212 (RSE 38%)
    lvc <- log(14.1)
    label("Central volume, referenced to unbound cefazolin (L)") # Table 2, row 'Central volume' V_C, semiphysiological column = 14.1 L (RSE 25%)
    lvp <- log(17.1)
    label("Peripheral volume, referenced to unbound cefazolin (L)") # Table 2, row 'Peripheral volume' V_P, semiphysiological column = 17.1 L (RSE 7%)
    lq <- log(0.436)
    label("Intercompartmental clearance, referenced to unbound cefazolin (L/min)") # Table 2, row 'Intercompartmental clearance' Q, semiphysiological column = 0.436 L/min (RSE 10%)

    # Protein binding. Results 3.1 eq 7: C_total = C_free / fu, i.e. a constant
    # (linear) binding model; nonlinear binding could not be identified.
    lfu <- log(0.291)
    label("Fraction of cefazolin unbound in plasma (unitless)") # Table 2, row 'Free fraction' F_U, semiphysiological column = 0.291 (RSE 9%)

    # ---- Embedded gestational creatinine-clearance trajectory ----
    # Methods eq 4: CrCL(t) = CrCL_0 + (CrCL_MAX * t) / (CrCL_50 + t), t in
    # gestational weeks. Carried here from the separately developed mixed
    # effect model (Table 3) so this PK model is self-contained; it also ships
    # standalone as inst/modeldb/endogenous/vanHasselt_2014_crcl_pregnancy.R.
    #
    # All three structural values are FIXED (Table 3 asterisk footnote: 'These
    # values were fixed during estimation of the mixed effect model; that is,
    # only random effects were estimated'), having come from the upstream
    # meta-analysis of gestational CrCL dynamics cited as reference 13.
    #
    # Names follow the operator ruling of 2026-09-21: this is the additive
    # Anderson-Holford shape (hill = 1) on a GESTATIONAL-AGE axis, which the
    # registered postnatal-age family does not cover, so every token stays
    # inside the crcl_ namespace. crcl_matspan is the SPAN above baseline, not
    # the plateau.
    lcrcl_ega0 <- fixed(log(97.83))
    label("Baseline nonpregnant creatinine clearance at EGA 0 (mL/min)") # Table 3, row 'Baseline CrCL' CrCL0 = 97.83 mL/min (upstream RSE 3.91%), FIXED
    lcrcl_matspan <- fixed(log(83.83))
    label("Maximum gestational increase in creatinine clearance, the span above baseline (mL/min)") # Table 3, row 'Maximum CrCL' CrCLMAX = 83.83 (upstream RSE 12.48%), FIXED
    lcrcl_ega50 <- fixed(log(13.3))
    label("Gestational age at half of the creatinine-clearance span (weeks)") # Table 3, row 'Time of half-maximum CrCL' CrCL50 = 13.3 weeks (upstream RSE 37.59%), FIXED

    # IIV. Methods eq 1 is exponential: P_i = P * exp(eta_i). Tables 2 and 3
    # report the between-subject variability as CV%, so the internal variance
    # is taken as (CV/100)^2. See the vignette Errata for the alternative
    # exact-lognormal reading log(1 + CV^2), which is material at the 101.5%
    # and 111.8% CVs carried by this model.
    etalcl ~ 0.010816 # Table 2, row 'Clearance' omega_CL, semiphysiological column = 10.4 CV% (RSE 70%) -> 0.104^2
    etalvc ~ 1.030225 # Table 2, row 'Central volume' omega_V1, semiphysiological column = 101.5 CV% (RSE 21%) -> 1.015^2
    etalvp ~ 0.462400 # Table 2, row 'Peripheral volume' omega_V2, semiphysiological column = 68 CV% (RSE 26%) -> 0.68^2
    etalfu ~ 0.036481 # Table 2, row 'Free fraction' omega_FU, semiphysiological column = 19.1 CV% (RSE 38%) -> 0.191^2

    # IIV on the embedded CrCL trajectory (Table 3). These are retained because
    # the paper's two-stage procedure carries INDIVIDUAL predicted CrCL curves
    # into the PK model (Methods 2.4.1 last paragraph: 'we estimated individual
    # empirical Bayes estimates describing the individual predicted changes in
    # CrCL during gestation'). They are also why omega_CL falls from 19.9 CV%
    # in the empirical sibling to 10.4 CV% here: part of the clearance
    # variability has moved into the CrCL trajectory.
    etalcrcl_ega0 ~ 0.098596 # Table 3, row 'Baseline CrCL' omega_CrCL0 = 31.4 CV% (RSE 16%) -> 0.314^2
    etalcrcl_matspan ~ 0.123201 # Table 3, row 'Maximum CrCL' omega_CrCLMAX = 35.1 CV% (RSE 17%) -> 0.351^2
    etalcrcl_ega50 ~ 1.249924 # Table 3, row 'Time of half-maximum CrCL' omega_CrCL50 = 111.8 CV% (RSE 121%) -> 1.118^2

    # Residual error. Methods eq 2 is a combined proportional + additive model on
    # the linear concentration scale, fitted separately for the free and total
    # cefazolin outputs. Table 2 heads this block 'Residual unexplained
    # variability variances', so each tabulated value is a VARIANCE and the SD
    # entered here is its square root.
    propSd <- 0.181108
    label("Proportional residual error, total cefazolin (fraction)") # Table 2, row 'Proportional, total concentration' sigma_TP, semiphysiological column = 0.0328 (variance) -> sqrt = 0.181108
    addSd <- 0.825227
    label("Additive residual error, total cefazolin (mg/L)") # Table 2, row 'Additive, total concentration' sigma_TA, semiphysiological column = 0.681 (variance) -> sqrt = 0.825227
    propSd_Cunbound <- 0.123693
    label("Proportional residual error, free cefazolin (fraction)") # Table 2, row 'Proportional, free concentration' sigma_FP, semiphysiological column = 0.0153 (variance) -> sqrt = 0.123693
    addSd_Cunbound <- 0.465833
    label("Additive residual error, free cefazolin (mg/L)") # Table 2, row 'Additive, free concentration' sigma_FA, semiphysiological column = 0.217 (variance) -> sqrt = 0.465833
  })

  model({
    # 1. Embedded creatinine-clearance trajectory, Methods eq 4.
    #    crcl_ega0 cancels out of the ratio formed in step 2, so only the
    #    SHAPE of this curve reaches the PK model; its mL/min scale does not.
    crcl_ega0 <- exp(lcrcl_ega0 + etalcrcl_ega0)
    crcl_matspan <- exp(lcrcl_matspan + etalcrcl_matspan)
    crcl_ega50 <- exp(lcrcl_ega50 + etalcrcl_ega50)
    CrCL <- crcl_ega0 + crcl_matspan * EGA / (crcl_ega50 + EGA)

    # 2. Normalised gestational change in renal function.
    #    Methods eq 6 divides by 'the baseline CrCL_i0 prior to start of
    #    pregnancy', which is exactly the EGA = 0 value of eq 4, i.e.
    #    crcl_ega0. The ratio is therefore 1 at EGA = 0 by construction.
    preg_cl <- CrCL / crcl_ega0

    # 3. Individual parameters.
    #    Corrigendum (Biomed Res Int 2015;2015:124035) Table 2 footnote b:
    #      CL = theta_CL0 + theta_CLPreg * (CrCL_ij(t) / CRCL_0)
    #    The ORIGINAL paper prints this equation against the EMPIRICAL column;
    #    the corrigendum states footnotes a and b were switched. The decisive
    #    evidence for the corrected assignment is that this equation IS Methods
    #    eq 6, which appears in Methods 2.4.2, the section that builds the
    #    semiphysiological model. It is also the reading that reproduces the
    #    base model's CL of 0.49 L/min at the cohort median GA of 33 weeks
    #    (0.142 + 0.212 * 1.6107 = 0.483 L/min, vs 0.529 under the switched
    #    reading) -- though note that clearance at the median GA does not on
    #    its own discriminate the two readings for BOTH models; see the
    #    vignette's corrigendum section for the comparison that does.
    #
    #    Methods eq 1 places the exponential IIV on the covariate-adjusted
    #    typical value, so the single reported omega_CL multiplies the SUM of
    #    the two clearance arms rather than either arm.
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal)
    cl <- (cl_nonren + cl_renal * preg_cl) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    fu <- exp(lfu + etalfu)

    # 4. Micro-constants. cl, q, vc and vp are all referenced to the unbound
    #    concentration, so the compartment amounts are TOTAL cefazolin and the
    #    micro-constants are formed in the usual way.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. Two-compartment disposition (Results 3.1: 'A two-compartmental model
    #    best described the data').
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 6. Observations. The model was fitted on the FREE cefazolin concentration
    #    (corrigendum clarification), so central/vc is the unbound concentration
    #    and the total concentration follows from the constant binding model of
    #    Results eq 7, C_total = C_free / fu. Equivalently CLtotal = CLfree * fu,
    #    as the corrigendum states.
    Cunbound <- central / vc
    Cc <- Cunbound / fu

    Cc ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
  })
}
