Chandasana_2023_dolutegravir <- function() {
  description <- "One-compartment oral population PK model with first-order absorption and elimination, estimated allometric weight scaling on CL/F and V/F, a fixed postmenstrual-age enzyme-maturation function on CL/F, formulation-specific absorption rate constant and formulation- and food-specific relative bioavailability for dolutegravir in infants, children and adolescents with HIV-1 pooled from IMPAACT P1093 and ODYSSEY (Chandasana 2023)"
  reference <- paste(
    "Chandasana H, Thapar M, Hayes S, Baker M, Gibb DM, Turkova A, Ford D,",
    "Ruel T, Wiznia A, Fairlie L, Bwakura-Dangarembizi M, Mujuru H, Alvero C,",
    "Farhad M, Hazra R, Townley E, Buchanan A, Bollen P, Waalewijn H,",
    "Colbers A, Burger D, Acosta EP, Singh R; for the IMPAACT P1093, ODYSSEY",
    "(PENTA 20) Study Teams. Population pharmacokinetic modeling of",
    "dolutegravir to optimize pediatric dosing in HIV-1-infected infants,",
    "children, and adolescents. Clin Pharmacokinet. 2023;62(10):1445-1459.",
    "doi:10.1007/s40262-023-01289-5.",
    "Structural equations, parameter estimates and the full random-effects",
    "covariance are taken from Chandasana 2023 Table 2 and from the final",
    "NONMEM control stream reproduced verbatim in Supplementary Text 1 of the",
    "Electronic Supplementary Material.",
    "The same fit was later applied without re-estimation to the",
    "abacavir/dolutegravir/lamivudine fixed-dose combination in IMPAACT 2019;",
    "see modellib('Chandasana_2024b_dolutegravir').",
    sep = " "
  )
  vignette <- "Chandasana_2023_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated allometric power effect on CL/F (exponent 0.455, 95% CI 0.418-0.492) and on V/F (exponent 0.556, 95% CI 0.514-0.598) with reference weight 70 kg (Chandasana 2023 Table 2 and the covariate-relationship equations printed beneath it; NONMEM control stream CLWT = (WT/70)**THETA(5), VWT = (WT/70)**THETA(6)). Neither exponent was fixed at the canonical 0.75 / 1.0; the paper's Discussion notes both came out lower than the canonical values while preserving their relative ratio. Evaluated as a time-varying covariate. Baseline weight 3.9-91.0 kg; across all analysis records weight spanned 4-96 kg (Chandasana 2023 Sect. 3.2).",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the sigmoidal enzyme-maturation function on CL/F, FMAT = PMA^Hill / (PMA^Hill + TM50^Hill), with TM50 = 52.2 postmenstrual weeks and Hill = 3.43, both FIXED from Anderson and Holford (2009) paracetamol values because paracetamol shares dolutegravir's UGT-dominated clearance pathways (Chandasana 2023 Table 2 footnote a and Sect. 3.2). The source model is written in postmenstrual WEEKS and derives them from postnatal age assuming term birth: PMA (weeks) = PNA (years) * 52 + 40 (Chandasana 2023 covariate-relationship block; control stream PMA = AGE*52+40). The canonical PAGE column is in months, so model() converts months to weeks with 1 month = 30.4375 / 7 = 4.348125 weeks before evaluating FMAT. FMAT saturates towards 1 in older children, so the maturation term is only material in infants; it is what separates the 10 mg (<6 months) and 15 mg (>=6 months) dispersible-tablet doses in the 6 to <10 kg weight band.",
      source_name        = "PMA (weeks)"
    ),
    FORM_DTG_DT = list(
      description        = "Dolutegravir dispersible-tablet / granule formulation indicator (1 = dispersible tablet or granules for oral suspension, 0 = film-coated tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (film-coated tablet, FCT)",
      notes              = "The dispersible tablet and the granules share a single set of absorption and bioavailability estimates because a healthy-adult study showed the two to be bioequivalent at the same dose (Chandasana 2023 Sect. 2.3, citing Buchanan 2017), so one indicator covers both. Two simultaneous effects, both multiplicative and both taken from the control stream: the absorption rate constant is multiplied by 2.04, giving Ka = 0.854 * 2.04 = 1.74 1/h for the dispersible tablet / granules against 0.854 1/h for the film-coated tablet (control stream KAFORM = THETA(17)**FORMK with FORMK = 1 when FORM >= 2; Chandasana 2023 covariate-relationship block and Table S1, which prints the product 1.74 with 95% CI 1.20-2.28); and relative bioavailability is multiplied by 1.53 (control stream FFFLAG = THETA(9)**FFLAG). Note that 2.04 and 1.53 are RATIOS to the film-coated-tablet reference, not absolute values.",
      source_name        = "FORM"
    ),
    FED = list(
      description        = "Dose administered without regard to food (1 = dosed without regard to food, 0 = dosed fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Relative bioavailability is multiplied by 1.10 (95% CI 1.03-1.17) when the dose was given without regard to food rather than fasted (Chandasana 2023 Table 2). The effect is NOT formulation-specific: the control stream forms F1 = 1 * THETA(4)**SFLAG * THETA(9)**FFLAG, whose own comment block enumerates fasted FCT = 1, fasted DT = THETA(9), fed FCT = THETA(4) and fed DT = THETA(9) * THETA(4), and Table S1 accordingly prints F for the dispersible tablet without regard to food as 1.68 = 1.10 * 1.53. The source category is 'without regard to food' (the protocol simply did not require fasting) rather than a controlled fed challenge with a defined meal, so the general FED indicator applies rather than FED_HIGHFAT; the paper attributes the smaller pediatric food effect relative to adults to exactly this lack of a defined meal. Operationally SFLAG = 1 for every sparse sample and for P1093 intensive samples from week 4 onward.",
      source_name        = "SFLAG"
    ),
    STUDY_ODYSSEY = list(
      description        = "ODYSSEY study indicator (1 = record from the ODYSSEY trial, 0 = record from IMPAACT P1093)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IMPAACT P1093)",
      notes              = "Selects between the two study-specific residual-error magnitudes of this pooled two-study analysis. Separate residual errors were estimated because the two trials used different bioanalytical methods, with lower limits of quantification of 5 ng/mL (P1093) and 9.7 ng/mL (ODYSSEY) (Chandasana 2023 Sects. 2.2 and 3.2). The control stream writes the combined error as Y = F + F*ERR(1)*(1-STUD) + ERR(2)*(1-STUD) + F*ERR(3)*STUD + ERR(4)*STUD, i.e. a proportional plus additive pair per study. The between-study difference is confined to residual error; no structural or covariate parameter differs between the studies.",
      source_name        = "STDY"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion variability on CL/F and Ka",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "The control stream carries FOUR inter-occasion CL/F slots (ETA(4)-ETA(7), $OMEGA BLOCK(1) plus SAME(3)) and TWO inter-occasion Ka slots (ETA(8)-ETA(9), $OMEGA BLOCK(1) plus SAME(1)), matching the four omega^2 IOV,CL rows and two omega^2 IOV,KA rows of Chandasana 2023 Table 2. All slots within a parameter share one magnitude, so occasions after the first are fixed to the first occasion's variance. In the source the Ka occasions were additionally gated on SERIAL = 1, i.e. only intensive-PK occasions carried inter-occasion variability on absorption; that gating is not reproduced here because it is a study-design flag rather than a model parameter. Set OCC = 0 on every record to switch inter-occasion variability off entirely, or OCC = 1 to reproduce the single steady-state occasion simulated for Chandasana 2023 Table 3.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 239L,
    n_studies      = 2L,
    age_range      = "0.170-17.5 years",
    age_median     = "6.00 years",
    weight_range   = "3.9-91.0 kg at baseline (4-96 kg across all analysis records)",
    weight_median  = "21.1 kg (mean)",
    sex_female_pct = 49.8,
    race_ethnicity = c(White = 6.3, Black = 79.5, Asian = 7.5, Other = 3.8, Unknown = 2.9),
    disease_state  = "Infants, children and adolescents living with HIV-1, treatment naive or treatment experienced, receiving oral dolutegravir once daily with optimized background antiretroviral therapy",
    dose_range     = "Oral dolutegravir once daily under fasting conditions: granules 4.8-32 mg, dispersible tablet 5-30 mg and film-coated tablet 20-50 mg across the WHO weight bands tested in P1093 and ODYSSEY (Chandasana 2023 Table S4)",
    regions        = "North America, Europe, Africa (South Africa, Zimbabwe, Uganda), Asia (Thailand)",
    notes          = "Pooled analysis of two trials: IMPAACT P1093 (N = 151), a phase I/II multicentre open-label single-arm trial in children and adolescents aged 4 weeks to <18 years, and ODYSSEY / PENTA 20 (N = 88), a phase II/III multicentre open-label randomised non-inferiority trial in children aged at least 28 days and <18 years, including nested PK sub-studies. 240 participants provided 2714 plasma concentrations; 64 were excluded (21 below the limit of quantification, 7 outlier or very low, 4 sample mix-up, 13 haemolysed, 6 non-adherent, 13 from a participant with multiple co-morbidities), leaving 2650 samples (1909 intensive, 741 sparse) from 239 participants (Chandasana 2023 Sect. 3.1 and Table 1). Baseline demographics are Chandasana 2023 Table 1; the ODYSSEY cohort was 100% Black and contributed no albumin, race or ethnicity detail. Estimation used NONMEM 7.3.0 with Monte Carlo importance sampling EM (IMPMAP). Predefined adult-matching exposure targets were a geometric-mean C24 of 0.995 ug/mL (target range 0.697-2.260 ug/mL, individual lower limit 0.500 ug/mL) and a geometric-mean AUC0-24 of 46 ug*h/mL (target range 37-134 ug*h/mL), with C24 the primary target (Chandasana 2023 Sect. 2.3)."
  )

  ini({
    # Structural parameters. Reference subject: 70 kg, fully mature (FMAT = 1),
    # fasted film-coated tablet. Chandasana 2023 Table 2; the control stream
    # (ESM Supplementary Text 1) estimates THETA(1)-THETA(3) on the log scale
    # (MU_1/MU_2/MU_3), which is why they are log-transformed here.
    lcl <- log(1.03);  label("Apparent oral clearance CL/F at WT = 70 kg and full maturation (L/h)")  # Chandasana 2023 Table 2: CL/F = 1.03 L/h (95% CI 0.980-1.07, %RSE 2.31)
    lvc <- log(13.6);  label("Apparent volume of distribution V/F at WT = 70 kg (L)")                 # Chandasana 2023 Table 2: V/F = 13.6 L (95% CI 13.0-14.3, %RSE 2.42)
    lka <- log(0.854); label("Absorption rate constant Ka, film-coated tablet (1/h)")                 # Chandasana 2023 Table 2: KA, FCT = 0.854 1/h (95% CI 0.686-1.06, %RSE 11.2)

    # Allometric weight exponents; both estimated, neither fixed at 0.75 / 1.0.
    e_wt_cl <- 0.455; label("Power exponent of body weight on CL/F, reference 70 kg (unitless)")  # Chandasana 2023 Table 2: CL/F~WT = 0.455 (95% CI 0.418-0.492, %RSE 4.15)
    e_wt_vc <- 0.556; label("Power exponent of body weight on V/F, reference 70 kg (unitless)")   # Chandasana 2023 Table 2: V/F~WT = 0.556 (95% CI 0.514-0.598, %RSE 3.87). Sect. 3.2 prints 0.566 in prose, a typographical slip: the tabulated point estimate and its 95% CI both centre on 0.556.

    # Sigmoidal enzyme-maturation function on CL/F, in postmenstrual WEEKS:
    #   FMAT = PMA^hill_mat / (PMA^hill_mat + tmat50^hill_mat)
    # Both parameters were held FIX at the paracetamol values of Anderson and
    # Holford (2009); neither was estimated from the dolutegravir data.
    tmat50   <- fixed(52.2); label("Postmenstrual age at half-maximal CL/F maturation TM50 (weeks)")     # Chandasana 2023 Table 2: TM50 = 52.2 week FIX (control stream THETA(8) = 52.2 FIX)
    hill_mat <- fixed(3.43); label("Hill coefficient of the CL/F enzyme-maturation function (unitless)") # Chandasana 2023 Table 2: Hill = 3.43 FIX (control stream THETA(7) = 3.43 FIX)

    # Formulation effect on the absorption rate constant. THETA(17) is a RATIO
    # to the film-coated tablet, not an absolute Ka: the control stream forms
    # KAFORM = THETA(17)**FORMK and KA = EXP(MU_3 + ETA(3)) * IOVKA * KAFORM,
    # so Ka for the dispersible tablet / granules is 0.854 * 2.04 = 1.74 1/h.
    e_form_dt_ka <- log(2.04); label("log ratio of the dispersible-tablet / granule absorption rate constant to the film-coated-tablet value (unitless; exp = 2.04)")  # Chandasana 2023 Table 2: KA~DT and granules = 2.04 (95% CI 1.41-2.67, %RSE 15.7); covariate-relationship block: KA (DT and granules) = 1.74 (95% CI 1.20-2.28) = 0.854 x 2.04. Table S1 tabulates 1.74 (1.20-2.28).

    # Formulation and food effects on relative bioavailability, both
    # multiplicative ratios to the fasted film-coated tablet (F = 1 FIX), and
    # both applied independently of each other: the control stream forms
    # F1 = 1 * THETA(4)**SFLAG * THETA(9)**FFLAG.
    e_food_fdepot    <- log(1.10); label("log relative bioavailability when dosed without regard to food, any formulation (unitless; exp = 1.10)")  # Chandasana 2023 Table 2: F, without regard to food, FCT = 1.10 (95% CI 1.03-1.17, %RSE 3.03); control stream THETA(4), not gated on formulation
    e_form_dt_fdepot <- log(1.53); label("log relative bioavailability of the dispersible tablet / granules versus the film-coated tablet (unitless; exp = 1.53)") # Chandasana 2023 Table 2: F, fasted DT/granules = 1.53 (95% CI 1.43-1.63, %RSE 3.26); control stream THETA(9)

    # IIV - a full 3x3 NONMEM $OMEGA BLOCK(3) on log-CL/F, log-V/F and log-Ka.
    # Chandasana 2023 Table 2 reports the omega^2 variances and the three
    # covariances directly, so no back-transformation from CV% is needed. The
    # tabulated CV% column follows Table 2's own footnote (CV% = 100*sqrt(w^2)
    # for w^2 <= 0.15, and 100*sqrt(exp(w^2)-1) above it); every value below
    # reproduces its printed CV% and its printed correlation exactly.
    etalcl + etalvc + etalka ~ c(0.0863,
                                 0.0499, 0.0698,
                                 0.0953, 0.138,  0.762)
    # Chandasana 2023 Table 2: omega^2 CL = 0.0863 (CV 29.4%, Shr 21.5%);
    #   Covar eta_CL,eta_V  = 0.0499 (R = 0.643); omega^2 V  = 0.0698 (CV 26.4%, Shr 22.2%);
    #   Covar eta_CL,eta_KA = 0.0953 (R = 0.372); Covar eta_V,eta_KA = 0.138 (R = 0.598);
    #   omega^2 KA = 0.762 (CV 107%, Shr 33.2%). Control stream: $OMEGA BLOCK(3).

    # IOV on CL/F over four occasions and on Ka over two, one shared magnitude
    # per parameter (NONMEM $OMEGA BLOCK(1) followed by SAME); occasions after
    # the first are fixed to the first occasion's variance.
    etaiov_lcl_1 ~ 0.115         # Chandasana 2023 Table 2: omega^2 IOV,CL = 0.115 (CV 33.9%), first of four occasions; control stream ETA(4)
    etaiov_lcl_2 ~ fix(0.115)    # Chandasana 2023 Table 2: second omega^2 IOV,CL row, same magnitude; control stream ETA(5) via $OMEGA BLOCK(1) SAME(3)
    etaiov_lcl_3 ~ fix(0.115)    # Chandasana 2023 Table 2: third omega^2 IOV,CL row, same magnitude; control stream ETA(6)
    etaiov_lcl_4 ~ fix(0.115)    # Chandasana 2023 Table 2: fourth omega^2 IOV,CL row, same magnitude; control stream ETA(7)
    etaiov_lka_1 ~ 0.610         # Chandasana 2023 Table 2: omega^2 IOV,KA = 0.610 (CV 91.7%), first of two occasions; control stream ETA(8)
    etaiov_lka_2 ~ fix(0.610)    # Chandasana 2023 Table 2: second omega^2 IOV,KA row, same magnitude; control stream ETA(9) via $OMEGA BLOCK(1) SAME(1)

    # Residual error - combined additive plus proportional, estimated
    # separately for each pooled study. Table 2's "Point estimate" column holds
    # the NONMEM $SIGMA VARIANCES; the SD / CV% column holds their square roots,
    # which are the standard-deviation-scale values nlmixr2 expects.
    propSd         <- 0.286;  label("Proportional residual error, IMPAACT P1093 (fraction)")  # Chandasana 2023 Table 2: proportional error P1093, variance 0.0818 (95% CI 0.0695-0.0941, %RSE 7.67), CV% = 28.6
    addSd          <- 0.0405; label("Additive residual error, IMPAACT P1093 (ug/mL)")         # Chandasana 2023 Table 2: additive error P1093, variance 0.00164 (95% CI -0.00142 to 0.00470, %RSE 95.1), SD = 0.0405
    propSd_odyssey <- 0.111;  label("Proportional residual error, ODYSSEY (fraction)")        # Chandasana 2023 Table 2: proportional error ODYSSEY, variance 0.0123 (95% CI 0.00787-0.0167, %RSE 18.4), CV% = 11.1
    addSd_odyssey  <- 0.300;  label("Additive residual error, ODYSSEY (ug/mL)")               # Chandasana 2023 Table 2: additive error ODYSSEY, variance 0.0900 (95% CI 0.0677-0.112, %RSE 12.7), SD = 0.300
  })

  model({
    # Inter-occasion variability, multiplexed by the occasion indicator. OCC = 0
    # zeroes every indicator and switches inter-occasion variability off.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 +
              oc3 * etaiov_lcl_3 + oc4 * etaiov_lcl_4
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2

    # Enzyme maturation on CL/F. The source model is written in postmenstrual
    # WEEKS; the canonical PAGE column carries postmenstrual age in MONTHS, so
    # convert with 1 month = 30.4375 / 7 = 4.348125 weeks.
    pma_weeks <- PAGE * 4.348125
    fmat <- pma_weeks^hill_mat / (pma_weeks^hill_mat + tmat50^hill_mat)

    # Formulation- and food-dependent relative bioavailability. The reference is
    # the fasted film-coated tablet (F = 1 FIX). The two ratios multiply:
    # fasted FCT 1.00, FCT without regard to food 1.10, fasted DT/granules 1.53,
    # DT/granules without regard to food 1.68.
    frel <- exp(e_food_fdepot * FED + e_form_dt_fdepot * FORM_DTG_DT)

    # Individual PK parameters at the 70 kg reference.
    ka <- exp(lka + e_form_dt_ka * FORM_DTG_DT + etalka + iov_ka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- frel

    # Study-specific residual error magnitudes.
    propSd_i <- propSd * (1 - STUDY_ODYSSEY) + propSd_odyssey * STUDY_ODYSSEY
    addSd_i  <- addSd  * (1 - STUDY_ODYSSEY) + addSd_odyssey  * STUDY_ODYSSEY

    Cc <- central / vc
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
