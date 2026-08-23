Jian_2025_peginterferon_alfa_2b <- function() {
  description <- "One-compartment quasi-equilibrium TMDD population PK model with first-order subcutaneous absorption and time-dependent interferon-receptor downregulation for peginterferon alfa-2b (Pegbing) in healthy subjects and chronic hepatitis B patients (Jian 2025 final model)"
  reference <- "Jian W, Yin Y, Chen R, Luo P, Wang T, Gu J, Du Z, Cai L, Bao T, Xue J, He R, Zhou T. Population Pharmacokinetic Model of Pegbing in Healthy Subjects and Chronic Hepatitis B Patients. CPT Pharmacometrics Syst Pharmacol. 2025;14:2014-2025. doi:10.1002/psp4.70104 (final-model NONMEM control stream: Supporting Information Data S2, file PSP4-14-2014-s002.docx; covariate-screening Table S1 and Figures S1-S9: Data S1, file PSP4-14-2014-s001.docx)"
  vignette <- "Jian_2025_peginterferon_alfa_2b"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Pegbing is a 40 kDa Y-shaped branched-PEG peginterferon alfa-2b (Xiamen
  # Amoytop Biotech). It shares the INN "peginterferon alfa-2b" with the
  # 12 kDa linear-PEG product (PegIntron) modelled in
  # Gupta_2006_peginterferon_alfa_2b.R, but is a structurally distinct
  # molecule with a much longer half-life (~7 days).

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance by the Cockcroft-Gault equation, raw (NOT BSA-normalized). Jian 2025 Table 1 footnotes (d) and (e) give the exact form used, with serum creatinine in umol/L: male CrCL = (140 - Age) * Weight / (0.818 * Creatinine); female CrCL = 0.85 * (140 - Age) * Weight / (0.818 * Creatinine). The 0.818 denominator is the standard Cockcroft-Gault 72 divided by the 88.4 umol/L-per-mg/dL creatinine conversion (72 / 88.4 = 0.814), so this is ordinary Cockcroft-Gault expressed for umol/L creatinine.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "LINEAR effect on apparent linear clearance, centered at 110.38 mL/min: CL = CL_TV * (1 + 0.00504 * (CRCL - 110.38)). This is the form in the executed final-model NONMEM control stream (Data S2): 'CLB_CLCR=(1+THETA(11)*(B_CLCR-110.38))' followed by 'TVCL = CLCOV*TVCL' -- i.e. Equation 1 (linear) of Methods Section 2.3, the PsN stepwise-covariate-model linear relation. FORM DISCREPANCY: the Jian 2025 Table 2 row LABEL prints the EXPONENTIAL form 'CL x exp(theta x (CrCL-110.38))' (Equation 2), which the control stream contradicts; the executed code is used. Baseline (not time-varying) CrCL is the model input -- the control stream reads the B_CLCR column, not the time-varying CLCR. Higher CrCL raises clearance. The effect reached statistical significance in the stepwise covariate model (Table S1 round F2, dOFV -7.43; backward elimination round B1 retained it at +10.01) but was judged NOT clinically relevant (Jian 2025 Section 3.5: the 10th / 90th CrCL percentile ratios were 0.904 and 1.21, inside the 80-120% relevance band).",
      source_name        = "CrCL"
    ),
    DIS_HEALTHY = list(
      description        = "Health-status indicator: 1 = healthy volunteer (Phase I SAD trial), 0 = chronic hepatitis B (CHB) patient (Phase II MAD trial NCT01143662). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (healthy volunteer)",
      notes              = paste(
        "POLARITY NOTE: the canonical column is DIS_HEALTHY (1 = healthy), but Jian 2025 uses HEALTHY as the",
        "reference category and reports the effect on the CHB arm. Jian 2025 Table 2 prints the row label",
        "'theta CHB on ka [(1 + theta) x ka for CHB patients]' with theta = -0.405, and Figure 4 confirms",
        "'Ratios were normalized by the typical values with health status of healthy'. The model therefore",
        "evaluates the CHB indicator as (1 - DIS_HEALTHY), so ka = ka_TV * (1 + e_chb_ka * (1 - DIS_HEALTHY))",
        "with ka_TV = 0.0101 /h being the HEALTHY typical value and CHB patients getting",
        "0.0101 * (1 - 0.405) = 0.00601 /h, i.e. the 40.5% lower ka reported in Section 3.4.",
        "This was the only covariate effect judged clinically relevant (Section 3.5)."
      ),
      source_name        = "Health status (healthy v.s. CHB)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on V/F and initially included during forward selection, but dropped in backward elimination and absent from the Jian 2025 Table 2 final model (Discussion, 'Based on the predefined inclusion and exclusion criteria, body weight was initially included as a covariate on V, but was later excluded during the stepwise selection process'). The paper attributes the non-significance to the narrow weight range (49-90 kg) and to CrCL - which is itself weight-derived - already carrying part of the effect. NOTE: the Abstract's claim that 'body weight affected the volume of distribution' contradicts both the Discussion and Table 2; the final model has no weight effect."
    ),
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes       = "Screened in the stepwise covariate model (Section 2.3) but not retained in the Jian 2025 Table 2 final model. Enters the model only indirectly through the Cockcroft-Gault CRCL calculation."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes       = "Screened (Section 2.3) but not retained in the final model. Enters only indirectly through the Cockcroft-Gault CRCL calculation (0.85 multiplier for females)."
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous",
      notes       = "Screened as one of the biochemical indicators in Section 2.3; not retained in the final model."
    ),
    ALT = list(
      description = "Alanine transaminase", units = "U/L", type = "continuous",
      notes       = "Screened as one of the biochemical indicators in Section 2.3; not retained in the final model. Strongly separates the two cohorts (healthy median 20 U/L vs CHB median 130 U/L, Table 1), so its effect is largely confounded with the retained health-status covariate."
    ),
    PLT = list(
      description = "Platelet count", units = "10^9/L", type = "continuous",
      notes       = "Screened as one of the hematological indicators in Section 2.3; not retained in the final model."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "peginterferon alfa-2b", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "peginterferon alfa-2b", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    total_target = list(
      # NOTE: unlike depot / central (amounts), this state holds a
      # drug-equivalent CONCENTRATION in ug/L. Jian 2025 Eq. 9 is written on
      # the concentration scale (its Ctol - Cfree term is a concentration) and
      # Table 2 reports R0 in ug/L, even though the paper's symbol glossary
      # loosely calls Rtol an "amount".
      analyte = "type I interferon receptor (IFNAR)", units = "ug/L",
      specimen = "not applicable", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 67L,
    n_studies      = 2L,
    age_range      = "18-46 years (Phase I healthy 20-40, median 30; Phase II CHB dense-PK 18-46, median 25)",
    age_median     = "Phase I healthy 30 years; Phase II CHB 25 years",
    weight_range   = "49-90 kg (Phase I healthy 49.0-76.0, median 59.0; Phase II CHB dense-PK 50.0-90.0, median 60.0)",
    weight_median  = "Phase I healthy 59.0 kg; Phase II CHB 60.0 kg",
    sex_female_pct = 43.3,
    disease_state  = "Pooled modelling dataset of 28 healthy volunteers (Phase I single-ascending-dose trial) and 39 chronic hepatitis B (CHB) patients with dense PK sampling (Phase II multiple-ascending-dose trial NCT01143662).",
    dose_range     = "Phase I: single SC doses of 45, 90, 180, or 270 ug. Phase II: 90, 135, or 180 ug SC once weekly for 48 weeks.",
    renal_function = "Normal; creatinine clearance (Cockcroft-Gault) 73.21-156.02 mL/min in healthy subjects (median 108.95) and 77.57-187.04 mL/min in CHB patients (median 120.49). Jian 2025 Table 1.",
    regions        = "Not stated in the paper; sponsor and investigators are China-based (Xiamen Amoytop Biotech Co. Ltd.; Peking University, Beijing).",
    external_validation = "An independent Phase II sparse-PK dataset of 115 CHB patients with 464 observations (90 / 135 / 180 ug weekly) was used for external pcVPC validation (Jian 2025 Figure 3b) and did NOT contribute to parameter estimation.",
    notes          = paste(
      "Baseline demographics are in Jian 2025 Table 1. The modelling dataset comprises 67 individuals and",
      "1013 observations (428 of them from the 28 healthy subjects). Plasma Pegbing was assayed by ELISA",
      "(BMS216/BMS216TEN, eBiosciences) with LLOQ 300 pg/mL; 16.5% of Phase I observations were below the",
      "LLOQ and none in Phase II, giving <5% BLQ overall, and BLQ records were excluded from the modelling",
      "dataset. Estimation was by FOCEI in NONMEM 7.5. The sex_female_pct field is the pooled modelling-set",
      "value (18 of 28 healthy + 11 of 39 CHB = 29 of 67).",
      "TABLE 1 ERRATUM: the validation-dataset sex percentages are transposed in the published Table 1",
      "(printed as Male 84 (27.0%) / Female 31 (73.0%), but 84/115 = 73.0% and 31/115 = 27.0%). This",
      "affects only the external-validation cohort description, not the modelling dataset or any parameter."
    )
  )

  ini({
    # ---- Structural PK (Jian 2025 Table 2, "Final model" block) -------------
    # CL/F and V/F are APPARENT: Jian 2025 Eq. 5 sets Xa(0) = Dose with no
    # separate bioavailability term, so F is absorbed into CL/F and V/F and
    # the model is run with F = 1. The Data S2 control stream confirms this --
    # it carries an FR term but pins it: THETA(8) is "(0, 1) FIX", its ETA(8)
    # is "0 FIX", and Phase I is hard-set to TVFR=1 -- so FR == 1 for every
    # subject in the final model.
    lcl  <- log(0.136);   label("Apparent linear clearance of free drug CL/F (L/h), at the reference CrCL of 110.38 mL/min")  # Jian 2025 Table 2 final model: CL/F = 0.136 L/h (RSE 5.3%; bootstrap median 0.136, 95% CI 0.124-0.148)
    lvc  <- log(2.95);    label("Apparent central volume of distribution V/F (L)")                                            # Jian 2025 Table 2 final model: V/F = 2.95 L (RSE 1.2%; bootstrap median 2.94, 95% CI 2.90-2.96)
    lka  <- log(0.0101);  label("First-order SC absorption rate constant ka (1/h), healthy-subject reference")                # Jian 2025 Table 2 final model: ka = 0.0101 /h (RSE 10.7%; bootstrap median 0.0100, 95% CI 0.0082-0.0119)

    # ---- Quasi-equilibrium TMDD receptor block ------------------------------
    # R0 and KD are drug-equivalent concentrations (ug/L), the same scale as
    # Ctol = central / vc, so the QE algebra runs on a single ug/L scale.
    lrbase <- log(0.328);          label("Baseline total IFN-receptor level R0 (ug/L, drug-equivalent)")                       # Jian 2025 Table 2 final model: R0 = 0.328 ug/L (RSE 25.3%; bootstrap 95% CI 0.255-0.474). The Table 2 bootstrap-median cell prints "0.0328", an evident decimal-point slip: it sits outside its own 95% CI (0.255-0.474) and is 10x below the estimate. The point estimate 0.328 is used.
    lkint  <- log(0.0827);         label("Internalization rate constant of the drug-receptor complex kint (1/h)")              # Jian 2025 Table 2 final model: kint = 0.0827 /h (RSE 10.5%; bootstrap median 0.0822, 95% CI 0.0706-0.0999)
    lkdeg  <- fixed(log(0.544));   label("First-order IFN-receptor degradation rate constant kdeg (1/h)")                      # Jian 2025 Table 2 final model: kdeg = 0.544 /h FIXED. Table 2 footnote (a): "fixed according to the PopPK modeling result of another similar PegIFNalpha" (ropeginterferon alfa-2b / Besremi, reference [1]); Discussion confirms the TMDD parameters were fixed to Besremi values.
    lkd    <- fixed(log(0.0493));  label("Quasi-equilibrium dissociation constant KD (ug/L)")                                  # Jian 2025 Table 2 final model: KD = 0.0493 ug/L FIXED. Table 2 footnote (b): "fixed according to the PopPK model in healthy subjects". Discussion: "The lower KD of Pegbing (0.0493 v.s. 0.142 ug/L) indicates a higher binding affinity to the IFN receptor compared to Besremi."
    lkdes  <- log(0.0068);         label("Exponential decay rate of IFN-receptor synthesis kR (1/day)")                        # Jian 2025 Table 2 final model: kR = 0.0068 /day (RSE 19.9%; bootstrap median 0.0066, 95% CI 0.0054-0.0095). NOTE the published units are 1/DAY while the model time unit is hours; the /24 conversion is applied in model().

    # ---- Covariate effects (Jian 2025 Table 2, final model) -----------------
    e_crcl_cl <-  0.00504;  label("Linear coefficient of CrCL on CL/F (per mL/min), centered at 110.38 mL/min")       # Jian 2025 Table 2 = 0.00504 (RSE 34.0%; bootstrap median 0.00490, 95% CI 0.00228-0.00794); THETA(11) of the Data S2 control stream, applied there as the LINEAR relation CLB_CLCR=(1+THETA(11)*(B_CLCR-110.38)). The Table 2 row LABEL prints the exponential form instead; the executed control stream governs. See the vignette Errata.
    e_chb_ka  <- -0.405;    label("Linear fractional effect of CHB disease status on ka (unitless)")                  # Jian 2025 Table 2 row "theta CHB on ka [(1 + theta) x ka for CHB patients]" = -0.405 (RSE 20.9%; bootstrap median -0.389, 95% CI -0.506 to -0.236); Section 3.4: "CHB patients had a 40.5% lower ka compared to healthy subjects"

    # ---- Inter-individual variability ---------------------------------------
    # SCALE OF THE REPORTED PERCENTAGES, settled by the Data S2 control stream.
    # Jian 2025 Table 2 reports IIV as a percentage and its abbreviations
    # footnote gives the back-transform as "SQRT(omega^2 - 1)", which is
    # impossible as printed (it would need omega^2 > 1). The final-model
    # $OMEGA block in Data S2 settles it: the omega VARIANCES are
    #   0.118 (CL), 0.232 (V), 0.197 (ka), 0.272 (kint),
    # and 100 * SQRT(variance) reproduces Table 2's 34.4 / 48.2 / 44.4 / 52.2
    # EXACTLY, to every printed digit. The log-normal CV reading
    # 100 * SQRT(exp(omega^2) - 1) would instead give 35.4 / 51.1 / 46.7 / 55.9
    # and matches none of them. The percentages are therefore the log-scale SD
    # (an approximate CV), and the variances below are taken from $OMEGA
    # directly rather than back-transformed.
    #
    # $OMEGA also fixes IIV on R0, kdeg, KD and F to 0, so only these four
    # parameters carry a random effect.
    etalcl   ~ 0.118  # Data S2 $OMEGA 1; Jian 2025 Table 2 final model IIV CL/F  = 34.4% (RSE 14.3%, shrinkage 14%) = 100*sqrt(0.118)
    etalvc   ~ 0.232  # Data S2 $OMEGA 2; Jian 2025 Table 2 final model IIV V/F   = 48.2% (RSE 16.0%, shrinkage 17%) = 100*sqrt(0.232)
    etalka   ~ 0.197  # Data S2 $OMEGA 3; Jian 2025 Table 2 final model IIV ka    = 44.4% (RSE 11.9%, shrinkage  9%) = 100*sqrt(0.197)
    etalkint ~ 0.272  # Data S2 $OMEGA 5; Jian 2025 Table 2 final model IIV kint  = 52.2% (RSE 18.8%, shrinkage 19%) = 100*sqrt(0.272)

    # ---- Residual unexplained variability -----------------------------------
    # Methods Section 2.2: "RUV was described by a mixed model". Data S2
    # $ERROR: Y = IPRED*(1 + EPS(2)) + EPS(1) with IPRED = A(2)/V*1000, i.e.
    # NONMEM works in pg/mL while this model works in ug/L (1 ug/L = 1000
    # pg/mL). $SIGMA gives the VARIANCES 15100 (additive, pg^2/mL^2) and 0.037
    # (proportional), whose square roots 122.88 pg/mL and 0.19235 reproduce
    # Table 2's final-model "123 pg/mL" and "19.2%" -- the same sqrt(variance)
    # scale as the IIVs above.
    propSd <- 0.1923538; label("Proportional residual error (fraction)")  # Data S2 $SIGMA 2 = 0.037 -> sqrt = 0.1923538; Jian 2025 Table 2 final model RUVprop = 19.2% (RSE 5.7%, shrinkage 10%; bootstrap median 19.1, 95% CI 17.2-20.9)
    addSd  <- 0.1228821; label("Additive residual error (ug/L)")          # Data S2 $SIGMA 1 = 15100 pg^2/mL^2 -> sqrt = 122.8821 pg/mL = 0.1228821 ug/L; Jian 2025 Table 2 final model RUVadd = 123 pg/mL (RSE 35.4%, shrinkage 10%; bootstrap median 119, 95% CI 56.0-202)
  })

  model({
    # ---- 1. Individual parameters -------------------------------------------
    # CrCL enters CL/F LINEARLY, centered at 110.38 mL/min. Data S2 control
    # stream: CLB_CLCR=(1+THETA(11)*(B_CLCR-110.38)); CLCOV=CLB_CLCR;
    # TVCL = CLCOV*TVCL. Figure 4 names 110.38 mL/min as the reference. The
    # Jian 2025 Table 2 row LABEL prints an exp() form instead -- see the
    # e_crcl_cl source-trace comment and the vignette Errata.
    cl   <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 110.38))
    vc   <- exp(lvc + etalvc)

    # Health status on ka. The canonical register column is DIS_HEALTHY
    # (1 = healthy), while Jian 2025 parameterises the effect on the CHB arm
    # with healthy as the reference; chb is therefore the complement.
    chb  <- 1 - DIS_HEALTHY
    ka   <- exp(lka + etalka) * (1 + e_chb_ka * chb)

    rbase <- exp(lrbase)
    kint  <- exp(lkint + etalkint)
    kdeg  <- exp(lkdeg)
    kd    <- exp(lkd)
    kdes  <- exp(lkdes)

    # ---- 2. IFN-receptor synthesis with exponential downregulation ----------
    # Jian 2025 Eq. 10 (ksyn = kdeg * R0) sets the baseline receptor level at
    # steady state; Eq. 11 makes the synthesis rate decay from the start of
    # treatment. kdes is published in 1/day and t is in hours, hence t / 24.
    #
    # SIGN OF Eq. 11. As printed, Eq. 11 reads ksyn = kdeg * R0 * exp(kR * t)
    # -- with NO minus sign (verified against the typeset PDF, and minus signs
    # do render correctly in Eqs. 5, 6 and 9 of the same paper). That printed
    # sign is a typographical error. The Data S2 final-model control stream
    # settles it directly and unambiguously:
    #     KSYN = REC0*KDEG*EXP(-THETA(9)*TIME)
    # with THETA(9) = 0.0068 -- the executed model decays. (Control-stream
    # TIME is in days: every rate THETA is multiplied by 24, and the example
    # dataset Data S3 runs TIME 0-335, i.e. 48 weeks.) The decay reading is
    # independently falsified four more ways by the paper's own content:
    #   (i)   the sentence introducing Eq. 11: "it is assumed that the
    #         synthesis rate of receptors DECREASES exponentially from the
    #         beginning of treatment";
    #   (ii)  the Figure 2 caption: "kR, exponential DECREASE rate of IFN
    #         receptor";
    #   (iii) Discussion: "IFN receptor levels would decrease by half after
    #         14.5 weeks" -- log(2) / 0.0068 per day = 101.9 d = 14.56 weeks,
    #         which reproduces 14.5 weeks only for a decay;
    #   (iv)  Discussion: "After 48 weeks of repeated dosing, receptor levels
    #         were reduced to approximately 10% of baseline" --
    #         exp(-0.0068 * 336 d) = 0.102, again only for a decay.
    # With the printed "+" sign the receptor would instead grow to ~9.8x
    # baseline by week 48, inverting the very observation the mechanism was
    # added to explain (higher exposure at week 48 than week 8, Figure 1c).
    # The decay form is therefore used. See the vignette's "Assumptions and
    # deviations" section.
    ksyn <- kdeg * rbase * exp(-kdes * t / 24)

    # ---- 3. Quasi-equilibrium free-drug algebra (Jian 2025 Eq. 7) -----------
    # central holds a drug AMOUNT (ug); total_target holds a drug-equivalent
    # CONCENTRATION (ug/L), the same scale as ctot.
    ctot   <- central / vc
    qecore <- ctot - total_target - kd
    cfree  <- 0.5 * (qecore + sqrt(qecore * qecore + 4 * kd * ctot))

    # Bound (complex) concentration. Jian 2025 writes the nonlinear
    # elimination term as Rtol * kint * Drugfree / (KD + Cfree); under the
    # quasi-equilibrium assumption Rtol * Cfree / (KD + Cfree) is identically
    # the complex concentration Ctol - Cfree, which is what Eq. 7 solves for
    # and what Eq. 9 uses. The mass-balance form below is algebraically
    # identical and better conditioned numerically.
    bound  <- ctot - cfree

    # ---- 4. ODE system (Jian 2025 Eqs. 5, 6, 9) -----------------------------
    # Eq. 5: dXa/dt = -ka * Xa, Xa(0) = Dose
    # Eq. 6: dDrugtol/dt = ka*Xa - (CL/V)*Drugfree - Rtol*kint*Drugfree/(KD+Cfree)
    #        with Drugfree = Cfree * V (Eq. 8), so (CL/V)*Drugfree = CL*Cfree
    #        and the TMDD term = kint * bound * V.
    # Eq. 9: dRtol/dt = ksyn - kdeg*Rtol - (kint - kdeg)*(Ctol - Cfree)
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - cl * cfree - kint * bound * vc
    d/dt(total_target) <-  ksyn - kdeg * total_target - (kint - kdeg) * bound

    # Jian 2025 Eq. 9 initial condition: Rtol(0) = R0.
    total_target(0) <- rbase

    # ---- 5. Observation and residual error ----------------------------------
    # The ELISA measures TOTAL drug, matching the paper's integrated state
    # Drugtol / Ctol (Eq. 6). Free drug Cfree is available as a derived
    # variable for users who need it.
    Cc <- ctot
    Cc ~ add(addSd) + prop(propSd)
  })
}
