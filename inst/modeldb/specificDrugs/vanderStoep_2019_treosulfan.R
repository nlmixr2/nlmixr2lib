vanderStoep_2019_treosulfan <- function() {
  description <- "Two-compartment IV-infusion population PK model for treosulfan in paediatric patients undergoing haematopoietic stem cell transplantation (van der Stoep 2019). Allometric body-weight scaling normalised to a 70 kg adult, exponents fixed at 0.75 on CL and Q and at 1 on V1 and V2, combined with a sigmoid Emax maturation function of postmenstrual age applied to BOTH CL and Q (TM50 = 38 weeks, Hill = 1.22; clearance reaches 90 percent of adult values at about 4 postnatal years). Full 4x4 correlated IIV block on CL, V1, V2 and Q, plus interoccasion variability on CL across the three conditioning days. Proportional residual error."
  reference <- "van der Stoep MYEC, Zwaveling J, Bertaina A, Locatelli F, Guchelaar HJ, Lankester AC, Moes DJAR. Population pharmacokinetics of treosulfan in paediatric patients undergoing hematopoietic stem cell transplantation. Br J Clin Pharmacol. 2019;85(9):2033-2044. doi:10.1111/bcp.13995"
  vignette <- "vanderStoep_2019_treosulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. van der Stoep 2019 Methods Section 2.2 collected blood
  # in serum tubes and assayed SERUM treosulfan by reversed-phase HPLC-UV, so
  # the assayed matrix is serum rather than plasma.
  compartmentData <- list(
    central = list(analyte = "treosulfan", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "treosulfan", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling normalised to a 70 kg adult (van der Stoep 2019 Equation 4, F_size = (BW/70 kg)^alpha), with alpha fixed at 0.75 for CL and Q and at 1 for V1 and V2 per Anderson and Holford. Applied to ALL four disposition parameters -- unlike the sibling Danielak 2017 treosulfan model, which leaves Q unscaled. Cohort median 15.6 kg (range 3.8-75.0; Table 1). The Discussion notes that standardising to 70 kg, well outside the observed paediatric weight range, inflates the RSE on the V2 and Q IIV terms (201 and 153 percent) relative to standardising to the cohort median 15.6 kg (53 and 46 percent), but 70 kg was retained in the published final model for comparability with other publications, so 70 kg is what this file encodes.",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "WEEKS, not the register-default months: van der Stoep 2019 Equation 5 and Table 2 state TM50 = 38 weeks, and the sigmoid maturation function is only meaningful on the week scale (the PAGE register entry in inst/references/covariate-columns.md explicitly provides for this; same treatment as Chen_2023_vancomycin.R and DAgate_2024_aciclovir.R). Methods Section 2.5: 'PMA was estimated by adding a gestational age of 40 weeks to postnatal age' -- i.e. the cohort carries an ASSUMED 40-week gestation rather than a recorded one, so PAGE = 40 + postnatal age in weeks. Drives F_mat = 1/(1 + (PMA/TM50)^-Hill) on BOTH CL and Q. Cohort postnatal age median 4.3 years (range 0.1-18.2), i.e. PAGE roughly 45-990 weeks; 33 of 91 subjects (36 percent) were aged 2 years or less.",
      source_name = "PMA"
    ),
    OCC = list(
      description = "Conditioning day (interoccasion-variability occasion index)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Values 1, 2, 3 identify the conditioning day on which the dose was given and the samples drawn, matching the DAY column and the IOV = DAY1*ETA(5) + DAY2*ETA(6) + DAY3*ETA(7) construction of the Supporting Information 4 control stream. Methods Section 2.4: 'each dose and subsequent sampling defined as a separate occasion'; interoccasion variability could be evaluated in 24 of the 91 patients, who consented to day-3 sampling in addition to the routine day-1 sampling. The three occasion variances are equal ($OMEGA BLOCK(1) 0.0194 followed by two SAME blocks), so occasions 2 and 3 are encoded with fixed() at the occasion-1 value. Decomposed inside model() into binary indicators oc1..oc3. For a single-occasion simulation pass OCC = 1 so the first IOV eta applies.",
      source_name = "DAY"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 91L,
    n_studies = 1L,
    age_range = "0.1-18.2 years (median 4.3); 33 of 91 (36%) were infants aged 2 years or less",
    weight_range = "3.8-75.0 kg (median 15.6)",
    bsa_range = "0.3-1.9 m^2 (median 0.7)",
    sex_female_pct = 36.3,
    disease_state = "Paediatric patients receiving treosulfan-based conditioning prior to allogeneic haematopoietic stem cell transplantation. Underlying disease: haemoglobinopathy 35 (38.5%), primary immune deficiency 26 (28.6%), haematological malignancy 17 (18.7%), bone marrow failure 11 (12.1%), other 2 (2.2%). Conditioning regimen: treosulfan + fludarabine + thiotepa 59 (64.8%), treosulfan + fludarabine 29 (31.9%), treosulfan + other (e.g. melphalan) 3 (3.3%).",
    dose_range = "Intravenous treosulfan over a 3 h infusion on 3 consecutive days. Patients older than 1 year received 14 g/m^2 per day (total 42 g/m^2); patients younger than 1 year received 10 or 12 g/m^2 per day (total 30 or 36 g/m^2). Administered doses: 14 g/m^2 in 73 (80.2%), 10 g/m^2 in 16 (17.6%), 12 g/m^2 in 2 (2.2%). Body surface area by the Mosteller formula (Supporting Information 2).",
    regions = "The Netherlands (Leiden University Medical Center, n = 63) and Italy (IRCCS Bambino Gesu Children's Hospital, Rome, n = 28)",
    renal_function = "eGFR by the revised Schwartz formula, capped at 120 mL/min/1.73 m^2: median 111 (range 16-120). Only 5 patients had an eGFR below 60. eGFR was a statistically significant covariate on CL in the stepwise screen (dOFV = -16.72) but was REJECTED from the final model because the prediction-corrected VPC worsened and the IIV of the PK parameters increased.",
    n_observations = "410 serum treosulfan concentrations across 91 patients (a subset of 35 full profiles from 28 patients was used for the limited-sampling analysis). 4 of 410 (1%) samples were below the 6.8 mg/L lower limit of quantification and were retained at their actual measured values.",
    notes = "Prospective, observational, multicentre study of patients treated between June 2011 and March 2017. Seven patients underwent a second transplantation with treosulfan conditioning a median of 8.5 months later and were entered twice, as distinct individuals, giving 91 analysis records from 84 unique patients. Patients without permanent central venous access were excluded. Baseline laboratory values (Table 1): creatinine 26 umol/L (8-166), albumin 38 g/L (20-52), haematocrit 0.291 L/L (0.199-0.384), haemoglobin 6.6 mmol/L (4.6-10.5). Observed treosulfan AUC(0-inf) 1658 mg*h/L (range 643-3371). Sex, underlying disease, conditioning regimen, haemoglobin, haematocrit and serum albumin were screened and none was retained. Parameter estimates from Table 2 'Final model' and the Supporting Information 4 NONMEM control stream; NONMEM 7.3.0, FOCE with INTERACTION, ADVAN6."
  )

  ini({
    # Structural parameters, all normalised to a 70 kg adult (Table 2 footnote
    # a). The Supporting Information 4 control stream carries these same values
    # in $THETA, so the two sources agree to the precision each prints.
    lcl <- log(18.8); label("Clearance at 70 kg reference, fully mature (L/h)") # Table 2 final model: Cl 18.8 L/h/70 kg (RSE 7%; bootstrap median 19.4, 95% CI 16.6-26.2); Data S4 $THETA (0, 18.8) TH_CL
    lvc <- log(20.2); label("Central volume at 70 kg reference (L)") # Table 2 final model: V1 20.2 L/70 kg (RSE 18%; bootstrap median 19.8, 95% CI 5.1-29.6); Data S4 $THETA (0, 20.2) TH_V
    lq <- log(21.3); label("Intercompartmental clearance at 70 kg reference, fully mature (L/h)") # Table 2 final model: Q 21.3 L/h/70 kg (RSE 31%; bootstrap median 22.0, 95% CI 9.7-68.9); Data S4 $THETA (0, 21.3) TH_Q
    lvp <- log(16.8); label("Peripheral volume at 70 kg reference (L)") # Table 2 final model: V2 16.8 L/70 kg (RSE 16%; bootstrap median 16.8, 95% CI 10.9-29.6); Data S4 $THETA (0, 16.8) TH_V2

    # Allometric exponents on body weight, Equation 4. The Methods fix them
    # rather than estimating them: 'When scaling clearance (Cl) and
    # intercompartmental clearance (Q) alpha is fixed to 0.75 and for volume of
    # distribution of the central (V1) and peripheral compartment (V2) alpha is
    # fixed to 1.' Neither Table 2 nor the control stream reports an exponent
    # estimate; Data S4 hardcodes **0.75 and (WT/70).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL (unitless)") # Methods Section 2.5 under Equation 4: alpha fixed to 0.75 for Cl; Data S4 FSIZE=(WT/70)**0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent on (WT/70) for Q (unitless)") # Methods Section 2.5 under Equation 4: alpha fixed to 0.75 for Q; Data S4 Q=THETA(3)*EXP(ETA(4))*FSIZE*FMAT re-uses the same FSIZE
    e_wt_vc <- fixed(1); label("Allometric exponent on (WT/70) for V1 (unitless)") # Methods Section 2.5 under Equation 4: alpha fixed to 1 for V1; Data S4 V=TVV*EXP(ETA(2))*(WT/70)
    e_wt_vp <- fixed(1); label("Allometric exponent on (WT/70) for V2 (unitless)") # Methods Section 2.5 under Equation 4: alpha fixed to 1 for V2; Data S4 V2=THETA(4)*EXP(ETA(3))*(WT/70)

    # Sigmoid Emax maturation of clearance on postmenstrual age, Equation 5:
    #   F_mat = 1 / (1 + (PMA / TM50)^-Hill)
    # Both terms are estimated and carry an RSE in Table 2. Table 2 prints the
    # Hill coefficient rounded to 1.2; the Data S4 control stream carries the
    # unrounded 1.22, which is used here (it reproduces the Table 3 dosing
    # recommendations slightly better -- see the vignette source-trace section).
    pma_tm50 <- 38; label("Postmenstrual age at 50% of adult CL (weeks)") # Table 2 final model: TM50 = 38 weeks (RSE 19%; bootstrap median 43, 95% CI 22.2-74.4); Data S4 $THETA (0, 38) TM50
    pma_hill <- 1.22; label("Hill coefficient of the postmenstrual-age CL maturation function (unitless)") # Data S4 $THETA (0, 1.22) Hill; Table 2 final model prints the same estimate rounded to 1.2 (RSE 34%; bootstrap median 1.1, 95% CI 0.3-3.2)

    # Correlated interindividual variability on all four disposition
    # parameters, from the Data S4 $OMEGA BLOCK(4) lower triangle in the
    # control stream's eta order ET_CL, ET_Vc, ET_Vp, ET_Q:
    #     0.101
    #     0.131   0.211
    #     0.0349  0.026   0.0299
    #    -0.0307 -0.0354 -0.05    0.171
    # Table 2 reports only the diagonals, as percentages, and does so on the
    # sqrt(omega^2) scale rather than the log-normal sqrt(exp(omega^2)-1)
    # scale: sqrt(0.101) = 31.8%, sqrt(0.211) = 45.9%, sqrt(0.0299) = 17.3%,
    # sqrt(0.171) = 41.4%, matching Table 2 exactly to the printed precision.
    # That agreement is what identifies this control stream as carrying the
    # FINAL estimates rather than initial values. Implied correlations:
    # CL-V1 0.897, CL-V2 0.635, CL-Q -0.234, V1-V2 0.327, V1-Q -0.186,
    # V2-Q -0.699. The block is positive definite but ill-conditioned
    # (smallest eigenvalue 8.4e-05), a direct consequence of the 70 kg
    # standardisation discussed in covariateData$WT$notes.
    etalcl + etalvc + etalvp + etalq ~ c(
      0.101,
      0.131, 0.211,
      0.0349, 0.026, 0.0299,
      -0.0307, -0.0354, -0.05, 0.171
    ) # Data S4 $OMEGA BLOCK(4); diagonals match Table 2 IIV CL 31.8%, V1 45.9%, V2 17.3%, Q 41.4% as sqrt(variance)

    # Interoccasion variability on CL across the three conditioning days,
    # Equation 2 and the Data S4 IOV = DAY1*ETA(5) + DAY2*ETA(6) + DAY3*ETA(7)
    # construction. $OMEGA BLOCK(1) 0.0194 followed by two BLOCK(1) SAME
    # blocks, so all three occasions share one variance; occasions 2 and 3 are
    # therefore fixed() at the occasion-1 value.
    etaiov_cl_1 ~ 0.0194 # Data S4 $OMEGA BLOCK(1) 0.0194; Table 2 final model IOV Cl 13.9% = sqrt(0.0194) (RSE 23%; bootstrap median 13.0, 95% CI 9.8-17.1)
    etaiov_cl_2 ~ fixed(0.0194) # Data S4 $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance
    etaiov_cl_3 ~ fixed(0.0194) # Data S4 $OMEGA BLOCK(1) SAME: occasion 3 shares the occasion-1 variance

    # Proportional residual error, Equation 3 (Y = YPRED * (1 + eps)). The
    # paper compared a proportional and a combined proportional-plus-additive
    # model and retained the proportional-only form.
    propSd <- 0.1233; label("Proportional residual error (fraction)") # Data S4 $SIGMA 0.0152 ER_Prop -> SD = sqrt(0.0152) = 0.1233; Table 2 final model sigma (proportional error) 12.3% (RSE 4%; bootstrap median 12.0, 95% CI 9.4-14.6)
  })

  model({
    # Reference body weight for the allometric terms (Equation 4 standardises
    # to a 70 kg individual).
    ref_wt <- 70

    # Size scaling, Equation 4: F_size = (BW / 70 kg)^alpha.
    size_cl <- (WT / ref_wt)^e_wt_cl
    size_q <- (WT / ref_wt)^e_wt_q
    size_vc <- (WT / ref_wt)^e_wt_vc
    size_vp <- (WT / ref_wt)^e_wt_vp

    # Sigmoid Emax maturation, Equation 5, written exactly as the Data S4
    # control stream writes it (FMAT = 1/(1+(PMA/TM50)**(-HILL))). With a
    # positive Hill this rises from 0 towards 1 as PMA increases, passing 0.5
    # at PMA = TM50. Applied to BOTH CL and Q: Equation 6 gives
    # Cl_tot = Cl_pop * F_size * F_mat and is followed by 'A similar model was
    # used for Q', and Results Section 3.3 reports 'The addition of maturation
    # of treosulfan Cl based on PMA on Cl and Q improved the model even
    # further'.
    fmat <- 1 / (1 + (PAGE / pma_tm50)^(-pma_hill))

    # Decompose the integer conditioning-day column into binary indicators to
    # multiplex the three interoccasion-variability etas on log-CL. For a
    # single-occasion simulation pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual disposition parameters, Data S4 $PK block.
    cl <- exp(lcl + etalcl + iov_cl) * size_cl * fmat
    vc <- exp(lvc + etalvc) * size_vc
    q <- exp(lq + etalq) * size_q * fmat
    vp <- exp(lvp + etalvp) * size_vp

    # Micro-constants for the linear two-compartment ODE (Data S4: K = CL/V,
    # K12 = Q/V, K21 = Q/V2).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment IV-infusion ODE system, Data S4 $DES. Dose enters
    # `central` directly through the rxode2 event-table infusion (rate or
    # duration); the paper administers every dose as a 3 h infusion.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Observation. Dose in mg / V in L -> Cc in mg/L, matching the serum
    # treosulfan concentrations quantified by HPLC-UV over a 6.8-500 mg/L
    # calibrated range (Methods Section 2.2). Data S4 sets S1 = V.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
