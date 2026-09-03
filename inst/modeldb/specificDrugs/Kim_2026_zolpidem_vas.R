Kim_2026_zolpidem_vas <- function() {
  description <- paste(
    "Sequential population PK-PD model linking oral zolpidem exposure to",
    "self-rated sedation on a 100 mm visual analog scale (VAS) in 30",
    "healthy Korean volunteers given a single 10 mg dose (Kim 2026). The PK",
    "layer is the companion `Kim_2026_zolpidem` model held entirely fixed;",
    "only the VAS parameters are estimated. VAS is a direct (no effect",
    "compartment) sigmoid Emax response bounded by the top of the scale,",
    "VAS = VASbaseline + (100 - VASbaseline) * Cp^HILL /",
    "(EC50^HILL + Cp^HILL). Baseline sedation showed no interpretable time",
    "trend, so it is described empirically by a piecewise-linear spline",
    "through four knots at nominal times 1, 2, 4 and 12 h -- measured on",
    "Day -1 and reapplied to time after dose on Day 1. Each knot is",
    "estimated as a fraction of the 100 mm scale held on the logit scale,",
    "which both keeps it inside (0, 100) and is the scale on which its",
    "additive inter-individual variability is placed. VAS is the most",
    "potent of the three endpoints (EC50 146 ug/L vs 205 for DSST and 282",
    "for CRT), which the authors read as subjective sedation preceding any",
    "measurable objective impairment. No covariate was retained. The",
    "`effect` compartment declared in the deposited control stream is",
    "vestigial: it carries no differential equation and is not referenced,",
    "so it is not reproduced here."
  )
  reference <- paste(
    "Kim HC, Sunwoo J, Yoon S, Jang I-J, Chung J-Y. Population",
    "Pharmacokinetic-Pharmacodynamic Analysis of Zolpidem in Healthy",
    "Volunteers: Covariate Contributions to Variability in CNS Responses.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70208.",
    "doi:10.1002/psp4.70208.",
    "PK layer fixed from the companion final PK model; see",
    "modellib('Kim_2026_zolpidem')."
  )
  vignette <- "Kim_2026_zolpidem"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot   = list(analyte = "zolpidem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "zolpidem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened in the stepwise covariate analysis (Kim 2026 Methods",
        "2.4, Tables S3 and S4) but not retained. Results 3.2.4: 'no",
        "significant covariates were identified that explained the IIV in",
        "VAS parameters'. The Discussion notes that observed baseline VAS",
        "scores were nonetheless higher in females (Figure S3). No point",
        "estimate reported. Source column SEX coded 1 = M, 2 = F."
      ),
      source_name = "SEX"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any VAS parameter; no point estimate reported.",
      source_name = "WEI"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any VAS parameter; no point estimate reported.",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL as reported by Kim 2026 (canonical register unit is g/L)",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any VAS parameter; no point estimate reported.",
      source_name = "ALB"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL as reported by Kim 2026 (canonical register unit is umol/L)",
      type = "continuous",
      notes = paste(
        "Kim 2026 Results 3.2.4: 'Total bilirubin was identified as a",
        "statistically significant covariate on H1 but was excluded from",
        "the final PK-VAS model because the relationship was not",
        "considered biologically plausible (Tables S3 and S4).' No point",
        "estimate is reported for the excluded relationship, so it cannot",
        "be encoded. Source column TBIL."
      ),
      source_name = "TBIL"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "CREA"
    ),
    DSST_BL = list(
      description = "Baseline digit symbol substitution test score at 1 h on Day -1",
      units = "score",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any VAS parameter; no point estimate reported.",
      source_name = "DSSTB1"
    ),
    CRT_BL = list(
      description = "Baseline choice reaction time at 1 h on Day -1",
      units = "msec",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any VAS parameter; no point estimate reported.",
      source_name = "CRTB1"
    ),
    VAS_SEDATION_BL = list(
      description = "Baseline sedation visual analog scale at 1 h on Day -1",
      units = "mm",
      type = "continuous",
      notes = paste(
        "Screened as a candidate covariate (Kim 2026 Methods 2.4, Table",
        "S1; control-stream $INPUT column VASB1, 'VAS baseline 1h') but",
        "not retained -- in this model the Day -1 baseline trajectory is",
        "instead estimated structurally as the four spline knots",
        "logitrbase_t1..t4, of which the 1 h knot is the model's own",
        "counterpart to this observed column. Cohort median 30.5 mm (range",
        "0-97). No point estimate reported."
      ),
      source_name = "VASB1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    age_range = "20-44 years",
    age_median = "29.5 years",
    weight_range = "50.0-83.0 kg",
    weight_median = "60.05 kg",
    sex_female_pct = 50,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral 10 mg immediate-release zolpidem",
    regions = "Korea",
    n_observations = "325 plasma zolpidem concentrations and 388 VAS observations",
    laboratory = paste(
      "Median (range) at baseline (Kim 2026 Table 1): sedation VAS at 1 h",
      "on Day -1 30.5 (0-97) mm."
    ),
    notes = paste(
      "Sedation VAS was measured at -23, -22, -20 and -12 h on Day -1",
      "(baseline) and at 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h on Day 1 (Kim",
      "2026 Methods 2.1). For the analysis the -23 h Day -1 time point was",
      "redefined as t = 1 h, placing the dose at t = 24 h; the four Day -1",
      "assessments therefore fall at nominal 1, 2, 4 and 12 h, which are",
      "exactly the spline knots. Results 3.2.4 notes VAS 'demonstrated the",
      "least consistent time-dependent trend and the highest variability'",
      "of the three endpoints."
    )
  )

  ini({
    # ====================================================================
    # PHARMACOKINETIC LAYER -- FIXED FROM THE FINAL PK MODEL
    # ====================================================================
    # Kim 2026 Results 3.2: each PD model was linked sequentially to the
    # final PK model. Data S4 $THETA(1)-$THETA(8) and $OMEGA(1,1)-
    # $OMEGA(6,6) all carry the FIX flag at the Data S1 / Table 2 final
    # values. See modellib("Kim_2026_zolpidem") for the estimating run.

    lka <- fixed(log(11.7))
    label("Absorption rate constant, depot -> central (Ka, 1/h)")                    # Kim 2026 Table 2: Ka = 11.7 1/h. Data S4 $THETA(3) = "(11.7) FIX".

    lcl <- fixed(log(18.0))
    label("Apparent clearance (CL/F, L/h)")                                          # Kim 2026 Table 2: CL/F = 18.0 L/h. Data S4 $THETA(1) = "(18) FIX".

    lvc <- fixed(log(64.0))
    label("Apparent central volume of distribution (Vc/F, L)")                       # Kim 2026 Table 2: Vc/F = 64.0 L. Data S4 $THETA(2) = "(64) FIX".

    lmtt <- fixed(log(0.25))
    label("Mean transit time of the absorption chain (MTT, h)")                      # Kim 2026 Table 2: MTT = 0.25 h. Data S4 $THETA(4) = "(0.25) FIX".

    lnn <- fixed(log(19.4))
    label("Number of transit compartments (NN, unitless)")                           # Kim 2026 Table 2: NN = 19.4. Data S4 $THETA(5) = "(19.4) FIX".

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (BIO, fraction)")                                # Kim 2026 Table 2: BIO = "1 Fixed". Data S4 $THETA(6) = "(1) FIX".

    etalka ~ fixed(2.48)
    # Kim 2026 Table 2: IIV Ka = 330.8 %CV. Data S4 $OMEGA(3,3) = "2.48 FIX".

    etalcl ~ fixed(0.0565)
    # Kim 2026 Table 2: IIV CL/F = 24.1 %CV. Data S4 $OMEGA(1,1) = "0.0565 FIX".

    etalmtt ~ fixed(0.133)
    # Kim 2026 Table 2: IIV MTT = 37.7 %CV. Data S4 $OMEGA(4,4) = "0.133 FIX".

    etalfdepot ~ fixed(0.0647)
    # Kim 2026 Table 2: IIV BIO = 25.9 %CV. Data S4 $OMEGA(6,6) = "0.0647 FIX".

    propSd <- fixed(0.189)
    label("Proportional residual error on plasma zolpidem (fraction)")               # Kim 2026 Table 2: proportional residual error = 18.9%. Data S4 $THETA(7) = "(0.189) FIX".

    addSd <- fixed(0.00001)
    label("Additive residual error on plasma zolpidem (ug/L)")                       # Data S4 $THETA(8) = "(0.00001) FIX"; a numerical guard, not reported in Kim 2026 Table 2.

    # ====================================================================
    # VAS DRUG-EFFECT PARAMETERS
    # ====================================================================
    # Typical values are the Final model column of Kim 2026 Table 3,
    # Sedation visual analog scale (VAS) block, reproduced in the Data S4
    # $THETA block. Results 3.2.4 confirms EC50 146 ug/L and HILL 1.92.

    lec50 <- log(146)
    label("Half-maximal effective concentration for sedation VAS (EC50, ug/L)")      # Kim 2026 Table 3: EC50 = 146 ug/L (RSE 12.7%; SIR median 147, 95% CI 118-184). Data S4 $THETA(9). The lowest of the three endpoint potencies -- Discussion: "the rank order of potency (VAS > DSST > CRT) suggests that zolpidem produces measurable sedative effects at lower concentrations".

    lhill <- log(1.92)
    label("Hill coefficient of the VAS concentration-response (HILL, unitless)")     # Kim 2026 Table 3: HILL = 1.92 (RSE 16.3%; SIR median 1.89, 95% CI 1.47-2.49). Data S4 $THETA(10).

    # ====================================================================
    # PIECEWISE-LINEAR BASELINE VAS SPLINE (LOGIT SCALE)
    # ====================================================================
    # Kim 2026 Results 3.2.4 / Methods 2.3: baseline VAS had no clear time
    # trend, so it is described by linear splines through four knots at
    # nominal times 1, 2, 4 and 12 h. Each typical knot TVHi is a FRACTION
    # of the 100 mm scale and is logit-transformed,
    #   TVLGTHi = LOG(TVHi / (1 - TVHi)),
    # which constrains the back-transformed value to (0, 100) mm and is
    # the scale carrying the additive IIV. Table 3 reports the knots on
    # the untransformed fraction scale as H1-H4 (mm/100); the ini() values
    # below are those fractions passed through the logit, so that the
    # back-transform expit(logitrbase_ti) recovers the printed Hi exactly.
    # Results 3.2.4 states the corresponding baselines as 30.3, 36.2, 19.8
    # and 12.7 mm.

    logitrbase_t1 <- log(0.303 / (1 - 0.303))
    label("Logit of the typical baseline VAS fraction at the 1 h knot (H1)")         # Kim 2026 Table 3: H1 = 0.303 mm/100 (RSE 17%; SIR median 0.302, 95% CI 0.229-0.382). Data S4 $THETA(11). Figure 1 legend: "H1, VAS measurement at 1 h on Day -1".

    logitrbase_t2 <- log(0.362 / (1 - 0.362))
    label("Logit of the typical baseline VAS fraction at the 2 h knot (H2)")         # Kim 2026 Table 3: H2 = 0.362 mm/100 (RSE 13.5%; SIR median 0.363, 95% CI 0.288-0.431). Data S4 $THETA(12). Figure 1 legend: "H2, VAS measurement at 2 h on Day -1".

    logitrbase_t3 <- log(0.198 / (1 - 0.198))
    label("Logit of the typical baseline VAS fraction at the 4 h knot (H3)")         # Kim 2026 Table 3: H3 = 0.198 mm/100 (RSE 26.6%; SIR median 0.198, 95% CI 0.124-0.279), corroborated by Results 3.2.4 ("19.8" mm). Data S4 $THETA(13) prints 0.199; the paper's value is used here because it is stated three times independently (final estimate, SIR median, Results text) against the control stream's one. Figure 1 legend: "H3, VAS measurement at 4 h on Day -1".

    logitrbase_t4 <- log(0.127 / (1 - 0.127))
    label("Logit of the typical baseline VAS fraction at the 12 h knot (H4)")        # Kim 2026 Table 3: H4 = 0.127 mm/100 (RSE 20.9%; SIR median 0.127, 95% CI 0.083-0.171). Data S4 $THETA(14). Figure 1 legend: "H4, VAS measurement at 12 h on Day -1".

    # ====================================================================
    # VAS INTER-INDIVIDUAL VARIABILITY AND RESIDUAL ERROR
    # ====================================================================
    # Data S4 $OMEGA(8,8) (HILL) is "0 FIX", so no etalhill appears.
    # Kim 2026 Table 3 footnote f: for H1-H4 "The IIV was modeled as an
    # additive effect on the logit-transformed baseline VAS (Day -1) and
    # expressed as omega^2 instead of %CV" -- so the four knot variances
    # below are read straight from Table 3, whereas the EC50 row is a %CV
    # that must be converted.

    etalec50 ~ 0.255
    # Kim 2026 Table 3: IIV EC50 = 53.9 %CV (RSE 44.7%, shrinkage 20%). Data S4 $OMEGA(7,7) = 0.255; 100*sqrt(exp(0.255)-1) = 53.9%.

    etalogitrbase_t1 ~ 1.03
    # Kim 2026 Table 3: IIV H1 = 1.03 (RSE 44.6%, shrinkage 20.6%; SIR median 1.13, 95% CI 0.60-1.81), reported as omega^2 per footnote f. Data S4 $OMEGA(9,9) = 1.03.

    etalogitrbase_t2 ~ 0.653
    # Kim 2026 Table 3: IIV H2 = 0.653 (RSE 49.0%, shrinkage 22.1%; SIR median 0.698, 95% CI 0.307-1.259), omega^2 per footnote f. Data S4 $OMEGA(10,10) = 0.653.

    etalogitrbase_t3 ~ 1.69
    # Kim 2026 Table 3: IIV H3 = 1.7 (RSE 37.9%, shrinkage 19.9%; SIR median 1.8, 95% CI 0.9-3.2), omega^2 per footnote f. Data S4 $OMEGA(11,11) = 1.69, the unrounded value behind Table 3's 1.7.

    etalogitrbase_t4 ~ 0.699
    # Kim 2026 Table 3: IIV H4 = 0.699 (RSE 40.9%, shrinkage 30.6%; SIR median 0.765, 95% CI 0.278-1.402), omega^2 per footnote f. Data S4 $OMEGA(12,12) = 0.699.

    propSd_VAS <- fixed(0.00001)
    label("Proportional residual error on sedation VAS (fraction)")                  # Data S4 $THETA(15) "V_Prop.RE" = "(0, 0.00001) FIX". Results 3.2.4 describes the VAS residual error as additive; the fixed 1e-5 proportional term is a numerical guard and is not reported in Kim 2026 Table 3.

    addSd_VAS <- 14.4
    label("Additive residual error on sedation VAS (mm)")                            # Kim 2026 Table 3: additive residual error = 14.4 mm (RSE 5.5%; SIR median 14.4, 95% CI 13.4-15.6). Data S4 $THETA(16) "V_Add.RE".
  })

  model({
    # --- 1. Individual parameters -----------------------------------------
    # No covariate was retained (Kim 2026 Results 3.2.4).
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    # Data S4 $PK holds each knot on the logit scale, adds its eta there,
    # then back-transforms:
    #   LGTHi = TVLGTHi + ETA(8+i);  Hi = 1 / (1 + EXP(-LGTHi))
    h1 <- expit(logitrbase_t1 + etalogitrbase_t1)
    h2 <- expit(logitrbase_t2 + etalogitrbase_t2)
    h3 <- expit(logitrbase_t3 + etalogitrbase_t3)
    h4 <- expit(logitrbase_t4 + etalogitrbase_t4)

    # --- 2. Micro-constants ------------------------------------------------
    kel <- cl / vc

    # --- 3. ODE system ------------------------------------------------------
    # Identical to the companion PK model; see modellib("Kim_2026_zolpidem").
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    # --- 4. Bioavailability -------------------------------------------------
    f(depot) <- 0

    # --- 5. Baseline VAS spline ---------------------------------------------
    # Data S4 $PK folds the Day -1 baseline period and the Day 1 post-dose
    # period onto one within-day clock:
    #   IF (NTIME .LT. 24) THEN  NCKT = NTIME        (Day -1 baseline)
    #   ELSE                     NCKT = NTIME - TDOS (time after dose)
    # with TDOS = 24 h by the Methods 2.1 time redefinition, so both
    # periods are indexed by the same knots at 1, 2, 4 and 12 h. Event
    # tables driving this model must therefore place the dose at t = 24 h.
    if (time < 24) {
      ckt <- time
    } else {
      ckt <- time - 24
    }

    # Data S4 $ERROR interpolates linearly between adjacent knots. The
    # control stream divides by the ACTUAL times of the knot observations
    # (CKT1..CKT4); those are the nominal 1, 2, 4 and 12 h whenever
    # sampling is on schedule, and they are written as such here. Below
    # the first knot the H1-H2 segment is extrapolated backwards, exactly
    # as the control stream's "NCKT .GE. 0 .AND. NCKT .LT. 2" branch does.
    # At and beyond the last knot the baseline is held flat at 100 * H4;
    # the control stream tests "NCKT .EQ. 12" only because no observation
    # falls after 12 h.
    if (ckt < 2) {
      vasbase <- 100 * (h1 + (h2 - h1) / (2 - 1) * (ckt - 1))
    } else if (ckt < 4) {
      vasbase <- 100 * (h2 + (h3 - h2) / (4 - 2) * (ckt - 2))
    } else if (ckt < 12) {
      vasbase <- 100 * (h3 + (h4 - h3) / (12 - 4) * (ckt - 4))
    } else {
      vasbase <- 100 * h4
    }

    # --- 6. Observations and residual error ---------------------------------
    Cc <- 1000 * central / vc

    # Data S4 $ERROR:
    #   EDRUG = (100 - VAS_BASE) * ((CP**HILL) / (EC50**HILL + CP**HILL))
    #   VAS   = VAS_BASE + EDRUG
    # Results 3.2.4: "the drug effect was also bounded by 100 mm by
    # multiplying the concentration-effect relationship by
    # (100 - VASbaseline)".
    VAS <- vasbase + (100 - vasbase) * Cc^hill / (ec50^hill + Cc^hill)

    Cc  ~ prop(propSd) + add(addSd)
    VAS ~ prop(propSd_VAS) + add(addSd_VAS)
  })
}
