Yoon_2023_vancomycin <- function() {
  description <- "Joint population PK/PD model for intravenous vancomycin in Korean hospitalised adults and children (Yoon 2023). Two-compartment PK with allometric scaling on CL/Q (0.75) and V/Vp (1) referenced to 70 kg; clearance is further scaled by a renal-function factor built from a creatinine-production-rate model that declines with age above 30 years and rises with age below it, by a sigmoid postmenstrual-age maturation factor, and by an exponential blood-urea-nitrogen effect together with proportional shifts for female sex, diabetes and renal disease. Vancomycin-induced nephrotoxicity is carried as a first-order exponential decline of creatinine clearance with time on therapy. The PD layer describes C-reactive protein with a proliferation compartment feeding two transit compartments and a circulating CRP compartment, whose production is stimulated linearly by a latent disease-severity state; severity grows first-order and is suppressed by cumulative vancomycin AUC, and the transit rate constant increases additively in patients with pneumonia. Vancomycin concentrations carry separate combined proportional-plus-additive residual errors for peak and trough samples."
  reference <- "Yoon S, Guk J, Lee S-G, Chae D, Kim J-H, Park K. Model-informed precision dosing in vancomycin treatment. Front Pharmacol. 2023;14:1252757. doi:10.3389/fphar.2023.1252757"
  vignette <- "Yoon_2023_vancomycin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # `severity` and `auc_total` are paper-mechanistic bookkeeping states with no
  # canonical role in inst/references/compartment-names.md. Per the standing
  # convention a new canonical compartment needs a second independent paper, so
  # they are declared here rather than registered on the strength of this one
  # extraction.
  paper_specific_compartments <- c("severity", "auc_total")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling referenced to 70 kg, exponent 0.75 on CL and Q and 1 on V and Vp (Yoon 2023 Eq. 1-2 and Methods section 2.3: 'The same allometry scaling was also applied to the inter-compartmental clearance (Q) and peripheral volume of distribution (Vp)'). Cohort median 59 kg (range 2.6-106 kg), Table 1.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Chronological age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model twice. (1) The creatinine production rate RCr = 64.2 * exp(kCr * (AGE - 30)) (Yoon 2023 Eq. 7), where kCr takes one of two estimated values according to whether AGE is at or above 30 years (-0.0127/yr) or below it (0.0193/yr) - Table 3 reports the two values as separate rows. (2) A centred exponential effect on the central volume, COVV = exp(kV * (AGE - 40)) (Eq. 16). Cohort median 60 years (range 0-93), Table 1.",
      source_name        = "age"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the sigmoid clearance-maturation factor Fmat = PMA^gamma / (PMA50^gamma + PMA^gamma) (Yoon 2023 Eq. 8, after Holford 2013), with PMA50 = 43.9 and gamma = 2.08 from Table 3. UNITS ARE WEEKS, not the register default of months: Yoon 2023 Methods section 2.4 states the maturation factor is 'associated with postmenstrual age (PMA) of up to 48 weeks (Anderson et al., 2007)', PMA50 = 43.9 is only meaningful on the weeks scale (43.9 months would place 50% maturation at 3.7 years and leave a 4-year-old only 54% mature, contradicting the paper's own restriction of Fmat to children under 4), and Table 1's PMA row - median 70, range 39-232 - is the under-4 paediatric subgroup in weeks (39 weeks is term gestation; 232 weeks is 4.45 years). Table 1 labels that row '(month)', which is a units typo; see the vignette Errata. For an adult, PAGE (weeks) is approximately AGE * 52.18 + 40, at which point Fmat is within 0.02% of 1.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum (plasma) creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as the denominator of the creatinine-clearance model CLCr = RCr / Cr * exp(-ktox * t) (Yoon 2023 Eq. 6). Measured by a rate-blanked compensated kinetic Jaffe method (Methods section 2.2). Cohort median 0.7 mg/dL (range 0.2-12.9), Table 1. The mg/dL unit is fixed by the paper's own arithmetic: RCr is a creatinine production rate in mg/h, so RCr/Cr has units of dL/h, and 64.2/1.0 = 64.2 dL/h = 107 mL/min is the expected creatinine clearance of a healthy 30-year-old.",
      source_name        = "Cr"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential effect on clearance centred at 15 mg/dL, exp(kBUN * (BUN - 15)) with kBUN = -0.00874 (Yoon 2023 Eq. 15, Table 3). The centring constant 15 is close to the cohort median of 15.25 mg/dL (Table 1). Note that the paper's Results describe a biphasic BUN relationship ('Vancomycin CL exhibited a gradual increase with BUN levels up to 15 mg/dL, followed by a subsequent decrease'); the printed final-model equation is monotonically decreasing and is what is encoded here - see the vignette Errata.",
      source_name        = "BUN"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional shift on clearance, (1 + thetaFEM * FEM) with thetaFEM = -0.199, i.e. 19.9% lower clearance in women (Yoon 2023 Eq. 15, Table 3). The paper defines 'FEM = 1 for female and 0 for male' in Results section 3.2, so the canonical SEXF orientation matches the source column with no value transformation. 41.1% of the PK cohort were female (Table 1). The corresponding dose reduction is quoted as 20% in the Table 5 footnote.",
      source_name        = "FEM"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Proportional shift on clearance, (1 + thetaDM * DM) with thetaDM = -0.151 (Yoon 2023 Eq. 15, Table 3). The paper defines 'DM = 1 for diabetes and 0 for no diabetes' in Results section 3.2. 28.4% of the PK cohort were diabetic (Table 1). The corresponding dose reduction is quoted as 15% in the Table 5 footnote.",
      source_name        = "DM"
    ),
    DIS_RENAL = list(
      description        = "History of renal disease indicator (any of acute kidney disease, chronic kidney disease, or other renal disease)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no renal disease)",
      notes              = "Proportional shift on clearance, (1 + thetaREN * REN) with thetaREN = -0.237 (Yoon 2023 Eq. 15, Table 3). The paper defines 'REN = 1 for renal disease and 0 for no renal disease' in Results section 3.2. Table 1 resolves the source column into four levels - none (72.9%), acute kidney disease (9.04%), chronic kidney disease (8.86%) and others (9.23%) - which the final model pools into a single binary indicator; the 27.1% with any renal disease are the REN = 1 group. The corresponding dose reduction is quoted as 23% in the Table 5 footnote. This effect is multiplicative WITH, not a replacement for, the continuous creatinine-driven renal-function factor Fren.",
      source_name        = "REN"
    ),
    DIS_PNEUMONIA = list(
      description        = "Pneumonia indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no pneumonia)",
      notes              = "Additive - not multiplicative - shift on the CRP transit rate constant, ktr = theta_ktr + thetaPNE * PNE with thetaPNE = 0.0058/h (Yoon 2023 Results section 3.3 equation, Table 4). The paper defines 'PNE = 1 indicates pneumonia and 0 indicates no pneumonia'. 32.0% of the 128-patient PD cohort had pneumonia (Table 2). Recorded as a general pneumonia indicator, not the hospital- or ventilator-acquired infection-type indicators DIS_HABP / DIS_VABP: Yoon 2023 uses pneumonia as a comorbidity / secondary-infection flag ('The hospitalization of up to 113 days was due to a secondary infection caused by pneumonia'), in the same medical-history family as DIS_DIAB and DIS_HYPERT, and never classifies it as hospital- or ventilator-associated.",
      source_name        = "PNE"
    )
  )

  covariatesDataExcluded <- list(
    DIS_HYPERT = list(
      description = "Hypertension comorbidity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a candidate covariate on CL and V (Yoon 2023 Methods section 2.4) but not retained in the final model. 52.4% of the PK cohort (Table 1)."
    ),
    DIS_SEPSIS = list(
      description = "Sepsis indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Yoon 2023 Methods section 2.4). 17.5% of the PK cohort (Table 1)."
    ),
    DIS_EDEMA = list(
      description = "Pleural effusion / oedema indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Yoon 2023 Methods section 2.4). Yoon 2023 Table 1 pools pleural effusion and oedema into a single row; 9.0% of the PK cohort."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Collected and tabulated (Yoon 2023 Table 1, median 2.9 g/dL, range 1-4.4) but not retained in the final model."
    )
    # Cardiovascular disease (47.2%), haematological malignancy (22.5%),
    # neutropenia (3.5%) and total serum protein (median 5.7 g/dL) were also
    # tabulated in Yoon 2023 Table 1 and screened as candidate covariates
    # (Methods section 2.4) without being retained. They are recorded in
    # `population$notes` rather than here because none of them has a canonical
    # column in inst/references/covariate-columns.md, and a covariate the final
    # model never uses is not a reason to mint one.
  )

  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    auc_total   = list(analyte = "vancomycin", units = "mg*h/L", specimen = "not applicable", verified = TRUE),
    prol        = list(analyte = "C-reactive protein", units = "mg/L", specimen = "not applicable", verified = TRUE),
    transit1    = list(analyte = "C-reactive protein", units = "mg/L", specimen = "not applicable", verified = TRUE),
    transit2    = list(analyte = "C-reactive protein", units = "mg/L", specimen = "not applicable", verified = TRUE),
    crp         = list(analyte = "C-reactive protein", units = "mg/L", specimen = "serum", verified = TRUE),
    severity    = list(analyte = "disease severity", units = "unitless", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 542L,
    n_subjects_pd  = 128L,
    n_studies      = 1L,
    age_range      = "0-93 years",
    age_median     = "60 years",
    weight_range   = "2.6-106 kg",
    weight_median  = "59 kg",
    sex_female_pct = 41.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Hospitalised patients treated with intravenous vancomycin and enrolled in therapeutic drug monitoring; comorbidities include hypertension (52.4%), cardiovascular disease (47.2%), diabetes (28.4%), haematological malignancy (22.5%), sepsis (17.5%), renal disease (27.1%) and pneumonia (32.0% of the PD subset)",
    dose_range     = "500-1500 mg per dose by intravenous infusion, dosing intervals 6-24 h",
    regions        = "Republic of Korea (Severance Hospital, Seoul)",
    notes          = "Retrospective analysis of electronic medical records; 1,526 vancomycin concentrations from 542 patients for PK (22 aged under 4 years, 18 aged 4-19 years, 502 adults) and 845 CRP measurements from 128 of those patients for PD. Demographics are Yoon 2023 Tables 1 (PK) and 2 (PD). The PD subset has median age 63 years, median weight 57.15 kg, 38.3% female and median CRP 73 mg/L. Mean length of stay 20 days (range 2-113 days). Vancomycin assayed by a KIMS immunoassay on a Roche Cobas c702 with a lower limit of quantitation of 4.0 mg/L; the cohort is Korean, so race is recorded as Asian although the paper does not tabulate it."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Reference subject is 70 kg, 30 years old
    # with a serum creatinine of 1.0 mg/dL, a BUN of 15 mg/dL, male,
    # non-diabetic, no renal disease and fully mature.
    # All values from Yoon 2023 Table 3.
    # ------------------------------------------------------------------
    lcl <- log(4.32);  label("Systemic clearance CL (L/h)")                             # Yoon 2023 Table 3 thetaCL = 4.32 (RSE 4.33%)
    lvc <- log(38.6);  label("Central volume of distribution V (L)")                    # Yoon 2023 Table 3 thetaV = 38.6 (RSE 2.69%)
    lq  <- log(3.93);  label("Inter-compartmental clearance Q (L/h)")                   # Yoon 2023 Table 3 thetaQ = 3.93 (RSE 9.31%)
    lvp <- log(66.8);  label("Peripheral volume of distribution Vp (L)")                # Yoon 2023 Table 3 thetaV2 = 66.8 (RSE 9.78%)

    # Allometric exponents. Yoon 2023 Eq. 1-2 write the exponents into the
    # model as the canonical 0.75 / 1 without reporting an estimate or an
    # RSE for either, so both are fixed.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")      # Yoon 2023 Eq. 2 (WT/70)^0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on V and Vp (unitless)")      # Yoon 2023 Eq. 1 (WT/70)

    # ------------------------------------------------------------------
    # Renal-function submodel (Yoon 2023 Eq. 5-7).
    #   RCr  = 64.2 * exp(kCr * (age - 30))            [mg/h]
    #   CLCr = RCr / Cr * exp(-ktox * t)               [dL/h]
    #   Fren = (CLCr / CLCr_TV)^lambda
    # kCr takes two separately estimated values either side of age 30;
    # Table 3 prints them as two rows, so both carry a stratum suffix.
    # ------------------------------------------------------------------
    kcr_agege30 <- -0.0127; label("Age coefficient of creatinine production rate for age >= 30 years (1/year)")  # Yoon 2023 Table 3 kCr (if age >= 30) = -0.0127 (RSE 14.3%)
    kcr_agelt30 <-  0.0193; label("Age coefficient of creatinine production rate for age < 30 years (1/year)")   # Yoon 2023 Table 3 kCr (if age < 30) = 0.0193 (RSE 37.9%)
    e_crcl_cl   <-  0.655;  label("Exponent on the creatinine-clearance ratio for CL (unitless)")                # Yoon 2023 Table 3 lambda = 0.655 (RSE 3.66%)
    ktox        <-  0.00598; label("First-order rate of nephrotoxicity-driven decline in creatinine clearance (1/day)")  # Yoon 2023 Table 3 ktox = 0.00598 (RSE 19.9%)

    # Maturation of clearance with postmenstrual age (Yoon 2023 Eq. 8).
    pma50_cl <- 43.9; label("Postmenstrual age at 50% mature CL (weeks)")               # Yoon 2023 Table 3 PMA50 = 43.9 (RSE 16.3%)
    h_pma_cl <- 2.08; label("Steepness factor of the CL maturation sigmoid (unitless)") # Yoon 2023 Table 3 gamma = 2.08 (RSE 63.5%)

    # ------------------------------------------------------------------
    # Covariate effects on CL and V (Yoon 2023 Eq. 15-16, Table 3).
    # ------------------------------------------------------------------
    e_age_vc       <-  0.00957; label("Exponential coefficient of age on V, centred at 40 years (1/year)")  # Yoon 2023 Table 3 kV = 0.00957 (RSE 8.88%)
    e_bun_cl       <- -0.00874; label("Exponential coefficient of BUN on CL, centred at 15 mg/dL (dL/mg)")  # Yoon 2023 Table 3 kBUN = -0.00874 (RSE 13.2%)
    e_dis_renal_cl <- -0.237;   label("Proportional change in CL with a history of renal disease (fraction)")  # Yoon 2023 Table 3 thetaREN = -0.237 (RSE 13.0%)
    e_sexf_cl      <- -0.199;   label("Proportional change in CL for female sex (fraction)")                   # Yoon 2023 Table 3 thetaFEM = -0.199 (RSE 12.7%)
    e_dis_diab_cl  <- -0.151;   label("Proportional change in CL with diabetes (fraction)")                    # Yoon 2023 Table 3 thetaDM = -0.151 (RSE 22.3%)

    # ------------------------------------------------------------------
    # PD: CRP proliferation-transit cascade (Yoon 2023 Eq. 9-14, Table 4).
    # ------------------------------------------------------------------
    lktr   <- log(0.0129);  label("Transit rate constant of the CRP cascade in the absence of pneumonia (1/h)")  # Yoon 2023 Table 4 theta_ktr = 0.0129 (RSE 4.90%)
    lkout  <- fixed(log(0.0365)); label("First-order elimination rate constant of circulating CRP (1/h)")        # Yoon 2023 Results 3.3: "kCRP was fixed at 0.0365 h-1 based on the prior knowledge that CRP's half-life was 19 h"; log(2)/19 = 0.0365. Table 4 prints "0.365 FIX", a decimal-point slip - see vignette Errata.
    lrbase <- log(110);     label("Initial circulating CRP concentration (mg/L)")                                # Yoon 2023 Table 4 CRP0 = 110 (RSE 8.4%)
    lkdis  <- log(0.00192); label("First-order disease-progression rate constant (1/h)")                         # Yoon 2023 Table 4 kD = 0.00192 (RSE 29.2%)
    scrp   <- 102;          label("Scaling factor of disease severity on CRP production (unitless)")             # Yoon 2023 Table 4 SCRP = 102 (RSE 6.56%)

    e_dis_pneumonia_ktr <- 0.0058;   label("Additive change in the CRP transit rate constant with pneumonia (1/h)")  # Yoon 2023 Table 4 thetaPNE = 0.0058 (RSE 17.2%)
    e_auc_kdis          <- 0.000239; label("Coefficient of cumulative vancomycin AUC on the drug effect (L/(mg*h))") # Yoon 2023 Table 4 alpha = 0.000239 (RSE 8.12%)

    # ------------------------------------------------------------------
    # Between-subject variability. Table 3 and Table 4 report these on the
    # CV(%) scale despite the omega-squared row labels, so the log-scale
    # variance is omega^2 = log(1 + CV^2).
    # ------------------------------------------------------------------
    etalcl    ~ 0.081277  # log(1 + 0.291^2); Yoon 2023 Table 3 omega2 CL  CV 29.1% (RSE 4.18%)
    etalvp    ~ 0.703147  # log(1 + 1.010^2); Yoon 2023 Table 3 omega2 V2  CV 101%  (RSE 8.04%)
    etalkdis  ~ 1.156430  # log(1 + 1.476^2); Yoon 2023 Table 4 omega2 KD  CV 147.6% (RSE 19.5%)
    etalrbase ~ 0.765089  # log(1 + 1.072^2); Yoon 2023 Table 4 omega2 CRP0 CV 107.2% (RSE 7.65%)

    # ------------------------------------------------------------------
    # Residual error. Yoon 2023 Table 3 estimates a separate combined
    # proportional-plus-additive error for trough samples (taken at the
    # start of an infusion) and peak samples (taken at the end of one) -
    # see Methods section 2.3. The two are encoded as two endpoints of the
    # same underlying concentration, selected by dvid in the data.
    # The additive rows carry the unit "mg/L" rather than "(mg/L)^2", and
    # the variability rows throughout both tables are reported on the
    # CV / SD scale, so these are standard deviations.
    # ------------------------------------------------------------------
    propSd_Cctrough <- 0.178; label("Proportional residual error on trough concentrations (fraction)")  # Yoon 2023 Table 3 sigma2 proportional trough CV 17.8% (RSE 5.57%)
    addSd_Cctrough  <- 0.956; label("Additive residual error on trough concentrations (mg/L)")          # Yoon 2023 Table 3 sigma2 additive trough = 0.956 mg/L (RSE 23.0%)
    propSd_Ccpeak   <- 0.110; label("Proportional residual error on peak concentrations (fraction)")    # Yoon 2023 Table 3 sigma2 proportional peak CV 11.0% (RSE 11.6%)
    addSd_Ccpeak    <- 4.47;  label("Additive residual error on peak concentrations (mg/L)")            # Yoon 2023 Table 3 sigma2 additive peak = 4.47 mg/L (RSE 8.59%)
    propSd_crp      <- 0.549; label("Proportional residual error on CRP (fraction)")                    # Yoon 2023 Table 4 sigma2 proportional CV 54.9% (RSE 2.71%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Renal-function factor Fren (Yoon 2023 Eq. 5-7).
    #
    # RCr is the creatinine production rate in mg/h; kCr is positive below
    # age 30 (production still rising with growth) and negative at or
    # above it (production falling with age), exactly as the paper states:
    # "kCr > 0 for age < 30 and kCr <= 0 for age >= 30".
    #
    # CLCr_TV, "the typical value of CLCr corresponding to the subject
    # with age of 30 years old", is not given a number anywhere in the
    # paper. It is taken here as 64.2 dL/h - the value of RCr at age 30
    # divided by a serum creatinine of 1.0 mg/dL - which is the only
    # reading that leaves no unreported constant, since the 64.2 then
    # cancels out of Fren entirely. It is also the value that best
    # reproduces the paper's own optimal-dose simulations: across all 480
    # cells of Table 5 it puts 98% of the published doses within one
    # 0.2 g/day grid step of the dose this model needs to hit the stated
    # 15 mg/L adult trough target. See the vignette Errata.
    #
    # ktox is reported per DAY (Table 3) while the model time unit is
    # hours, hence t/24.
    # ------------------------------------------------------------------
    kcr  <- kcr_agege30 * (AGE >= 30) + kcr_agelt30 * (AGE < 30)
    rcr  <- 64.2 * exp(kcr * (AGE - 30))
    crcl <- rcr / CREAT * exp(-ktox * t / 24)
    fren <- (crcl / 64.2)^e_crcl_cl

    # 2. Maturation factor Fmat (Yoon 2023 Eq. 8); PAGE in weeks.
    #    Applied continuously rather than gated at age 4: the paper's prose
    #    restricts Fmat to children under 4 but its printed equation has no
    #    gate, and a gate would introduce a 2.7% discontinuity at exactly
    #    4 years. Fmat is 0.973 at 4 years and within 0.02% of 1 in adults.
    fmat <- PAGE^h_pma_cl / (pma50_cl^h_pma_cl + PAGE^h_pma_cl)

    # 3. Remaining covariate functions (Yoon 2023 Eq. 15-16).
    covcl <- exp(e_bun_cl * (BUN - 15)) *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_dis_diab_cl * DIS_DIAB) *
      (1 + e_dis_renal_cl * DIS_RENAL)
    covv <- exp(e_age_vc * (AGE - 40))

    # 4. Individual PK parameters (Yoon 2023 Eq. 1-2).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * fren * fmat * covcl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp * covv
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # 5. Individual PD parameters (Yoon 2023 Table 4 and Results 3.3).
    #    ktr rises ADDITIVELY with pneumonia: ktr = theta_ktr + thetaPNE * PNE.
    ktr   <- exp(lktr) + e_dis_pneumonia_ktr * DIS_PNEUMONIA
    kcrp  <- exp(lkout)
    crp0  <- exp(lrbase + etalrbase)
    kdis  <- exp(lkdis + etalkdis)

    # 6. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 7. PK ODEs.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc

    # 8. Cumulative vancomycin AUC drives the drug effect (Yoon 2023
    #    Eq. 14: EDrug = alpha * AUC, with "AUC obtained from the
    #    developed PK model"). The AUC is cumulative from the start of
    #    treatment: alpha = 0.000239 makes EDrug reach 1 at an AUC of
    #    about 4,200 mg*h/L, i.e. roughly ten days of typical dosing,
    #    which is exactly where Figure 6 shows the dose arms separating.
    d/dt(auc_total) <- Cc

    edrug <- e_auc_kdis * auc_total

    # 9. CRP cascade (Yoon 2023 Eq. 9-13). Kin = Kout = ktr, "assumed for
    #    numerical simplicity" (Methods section 2.5).
    d/dt(prol)     <- ktr * (1 + scrp * severity) - ktr * prol
    d/dt(transit1) <- ktr * (prol - transit1)
    d/dt(transit2) <- ktr * (transit1 - transit2)
    d/dt(crp)      <- ktr * transit2 - kcrp * crp
    d/dt(severity) <- kdis * severity * (1 - edrug)

    # 10. Initial conditions. Disease severity starts at 1 and the
    #     proliferation and transit chain starts at the steady state that
    #     severity implies, 1 + SCRP; circulating CRP starts at the
    #     separately estimated CRP0. The chain is therefore NOT in
    #     equilibrium with CRP0, which is what produces the steep early
    #     fall from 110 mg/L to a plateau near ktr * (1 + SCRP) / kCRP
    #     that Figure 6 shows over the first 200 h. Initialising the chain
    #     to be consistent with CRP0 instead removes that fall and is
    #     falsified by Figure 6.
    severity(0) <- 1
    prol(0)     <- 1 + scrp
    transit1(0) <- 1 + scrp
    transit2(0) <- 1 + scrp
    crp(0)      <- crp0

    # 11. Observations. Trough and peak vancomycin samples share the
    #     structural model but carry separate residual errors.
    Cctrough <- Cc
    Ccpeak   <- Cc

    Cctrough ~ prop(propSd_Cctrough) + add(addSd_Cctrough)
    Ccpeak   ~ prop(propSd_Ccpeak) + add(addSd_Ccpeak)
    crp      ~ prop(propSd_crp)
  })
}
