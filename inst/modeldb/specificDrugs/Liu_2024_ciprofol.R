Liu_2024_ciprofol <- function() {
  description <- "Three-compartment intravenous population PK model for ciprofol (HSK3486) pooled across 12 clinical trials in healthy adults, elective surgical patients and mechanically ventilated adults sedated in the intensive care unit (Liu 2024; 462 subjects, 3918 plasma concentrations). Clearance depends on body weight (power, 62 kg reference), sex, blood sampling site (venous samples 17.7% lower than arterial), patient population (non-ICU 21.1% lower than ICU) and study design (short-term bolus / short-infusion studies vs long-term maintenance-infusion studies); the central volume depends on age (power, 33 years reference); the study-design flag also shifts the deep peripheral volume, the shallow inter-compartmental clearance and the shallow peripheral volume. Diagonal log-normal IIV on all six disposition parameters; additive residual error on log-transformed concentrations with a separate magnitude for venous and arterial samples."
  reference <- paste(
    "Liu L, Wang K, Sun Z, Yan P, Hu M, Liu X, Chen M, Wu N, Xiang X.",
    "Pharmacokinetics and exposure-safety relationship of ciprofol for",
    "sedation in mechanically ventilated patients in the intensive care unit.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(5):823-836.",
    "doi:10.1002/psp4.13121.",
    "PMCID: PMC11098162.",
    sep = " "
  )
  vignette <- "Liu_2024_ciprofol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: the paper measures ciprofol in plasma
  # by a validated LC-MS/MS assay (Methods, 'Study design and patient
  # population'), and the disposition compartments hold unchanged ciprofol.
  compartmentData <- list(
    central = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL only, normalised to 62 kg: CL ~ (WT / 62)^0.578 (Equation 3; supplementary control stream 'THETA(7)*LOG(WT/62)'). Allometric scaling of weight on every parameter was tested on the base model and not retained (Table S3 run 3, p = 0.037 against the 0.01 entry threshold). Analysis range 41.6-97.9 kg, median 61.7 kg (Table 1).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the central volume V1, normalised to 33 years: V1 ~ (AGE / 33)^0.496 (Equation 3; control stream 'THETA(8)*LOG(AGE/33)'). Analysis range 18-77 years, median 38 years (Table 1).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on CL, exp(0.089 * SEXF), i.e. female CL 9.3% HIGHER than male (Equation 3 '+0.089 * sex = female'; Table 2 'Impact of gender (female) on CL' = 1.09). The Results text says 'Women exhibited 9.31% lower CL compared to men'; the equation, Table 2 and the Figure 1 sex panel (female empirical-Bayes CL above male) all show the opposite direction, and the covariate sensitivity statement that sex changes exposure by < 7.1% reproduces only with the equation's direction, so the equation is encoded. The control stream column is GEND with a +THETA(14) coefficient.",
      source_name = "GEND"
    ),
    SAMPLE_ARTERIAL = list(
      description = "Per-observation blood sampling-site indicator, 1 = arterial sample, 0 = venous sample.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (venous)",
      notes = "The paper's flag is the reverse: the control stream sets SITE0 = 1 for a venous sample (SITE != 1) and 0 for an arterial one, so venous = 1 - SAMPLE_ARTERIAL. Two roles: (1) an exponential effect on CL, exp(-0.195 * (1 - SAMPLE_ARTERIAL)), venous CL 17.7% lower than arterial (Equation 3; Table 2 'Impact of sampling site (vein) on CL' = 0.823); (2) the switch between the venous and arterial log-scale residual SDs (control stream 'Y = IPRED + EPS(1)*(1-SITE) + EPS(2)*SITE', SITE = 1 arterial). In NONMEM $PK the CL shift is evaluated per record, so a subject who contributed both arterial and venous samples carries the CL of the site of each record; supply the column per observation row to reproduce that, or hold it constant per subject for a simulation that assumes one sampling site. Of the 462 subjects, 230 were sampled venously only, 146 arterially only and 86 at both sites (Table 1).",
      source_name = "SITE"
    ),
    DIS_CRITILL = list(
      description = "Critical-illness indicator, 1 = mechanically ventilated patient sedated in the intensive care unit, 0 = healthy volunteer or elective surgical patient.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-ICU)",
      notes = "The paper's flag is the reverse: the control stream variable ICU equals 1 for the non-ICU studies and 0 for studies HSK3486-205 and HSK3486-307 (the two ICU studies), so the paper's ICU = 1 - DIS_CRITILL. Effect on CL: exp(-0.237 * (1 - DIS_CRITILL)), i.e. non-ICU CL 21.1% lower than ICU CL (Equation 3; Table 2 'Impact of ICU on CL' = 0.789; Results 'The model estimated CL values were 21.1% lower ... in long-term infusion in non-ICU'). The typical-value parameters of the model are those of an ICU patient. Set per study, not per subject characteristic: all 112 subjects of the two ICU studies are 1, all 350 other subjects are 0 (Table 1).",
      source_name = "ICU"
    ),
    REGIMEN_CONTINF = list(
      description = "Study-design indicator, 1 = long-term maintenance-infusion study, 0 = bolus or short-term-infusion study.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (bolus / short-term infusion)",
      notes = "The paper's flag FLG is the reverse: FLG = 0 for the long-term infusion studies HSK3486-102, -204, -205, -306 and -307, FLG = 1 for the remaining seven studies (control stream $PK; Table 2 footnotes), so FLG = 1 - REGIMEN_CONTINF. The flag is assigned per study, not derived from the infusion duration of each dose: HSK3486 SAD_03 gave a 1 min bolus plus a 30 min infusion and is FLG = 1 (short-term). Effects: CL exp(0.229 * FLG), V2 exp(-2.08 * FLG), Q3 exp(0.802 * FLG), V3 exp(-1.62 * FLG) (Equation 3; Table 2 multiplicative factors 1.26, 0.125, 2.23 and 0.198), so the typical-value parameters of the model are those of a long-term infusion study. The paper attributes the much larger peripheral volumes in long-term infusions to progressive distribution of the lipophilic drug into adipose tissue (Discussion).",
      source_name = "FLG"
    )
  )

  # Covariates screened in the stepwise search (Methods 'PopPK model
  # development'; Table S2) but not retained in the final model.
  # Documentation only -- not referenced in model().
  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Race indicator, 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V1 (Table S2); not retained. 394 Asian / 68 non-Asian subjects, the non-Asian subjects all from the three Australian phase 1 studies (Table 1, Table S1)."
    ),
    DIS_HEALTHY = list(
      description = "Health-status indicator, 1 = healthy volunteer, 0 = patient.",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V1 (Table S2 'Health status (healthy subject vs patient)'); not retained. 125 healthy volunteers (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault).",
      units = "mL/min",
      type = "continuous",
      notes = "Tested on CL (Table S2); not retained. Computed as (140 - age) x weight / (72 x creatinine [umol/L] / 88.4), x 0.85 for women (Table 1 footnote). Range 44.4-212 mL/min, median 110."
    ),
    TPRO = list(
      description = "Total serum protein.",
      units = "g/L",
      type = "continuous",
      notes = "Tested on CL (Table S2); not retained. Range 31.5-87.9 g/L, median 70.6 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL (Table S2); not retained. Range 7-408 U/L, median 21 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL (Table S2); not retained. Range 3-416 U/L, median 19 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = "Tested on CL (Table S2); not retained. Range 3-197 umol/L, median 11 (Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase.",
      units = "U/L",
      type = "continuous",
      notes = "Listed among the covariates examined in Methods 'PopPK model development' (absent from Table S2 and Table 1); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 462L,
    n_studies = 12L,
    age_range = "18-77 years",
    age_median = "38 years",
    weight_range = "41.6-97.9 kg",
    weight_median = "61.7 kg",
    sex_female_pct = 51.08,
    race_ethnicity = c(Asian = 85.28, `Non-Asian` = 14.72),
    disease_state = "Pooled: 125 healthy volunteers (six phase 1 studies), 225 adults undergoing elective surgery (four phase 2/3 studies) and 112 mechanically ventilated adults sedated in the intensive care unit (studies HSK3486-205 and HSK3486-307).",
    dose_range = "Single intravenous injections of 0.128-0.9 mg/kg; induction boluses of 0.1-1.0 mg/kg followed by maintenance infusions of 0.06-2.4 mg/kg/h for 0.5 h to 24 h; ICU sedation with a 0.1-0.2 mg/kg loading dose and a 0.3 mg/kg/h initial maintenance rate titrated within 0.06-0.8 mg/kg/h for 6-24 h (Table S1).",
    regions = "China (nine studies) and Australia (three phase 1 studies).",
    notes = "Table 1 baseline characteristics: 3918 plasma concentrations (3457 non-ICU, 461 ICU); BLQ samples (5.59%, LLOQ 5 ng/mL) excluded. ICU subjects were older (median 57 vs 32 years) and had lower total protein (55.8 vs 72.1 g/L) than non-ICU subjects. Sampling site: 230 venous only, 146 arterial only, 86 both. Fitted in NONMEM 7.4.3 (FOCE-I) with PsN 4.2; evaluated by GOF, pcVPC, NPC and bootstrap (Table 2)."
  )

  ini({
    # Typical values for the reference subject of Equation 3: a 62 kg,
    # 33-year-old male ICU patient in a long-term infusion study, sampled
    # arterially. The paper's numbering maps to the canonical names as
    #   V1 -> vc, Q2 -> q, V2 -> vp, Q3 -> q2, V3 -> vp2.
    # Equation 3 prints the same values on the log scale (4.39, 1.26, 4.30,
    # 3.95, 3.35, 6.04); Table 2 carries one more significant figure.
    lcl <- log(80.9); label("Clearance CL (L/h)") # Table 2: CL = 80.9 L/h (RSE 3.33%); Eq. 3 intercept 4.39
    lvc <- log(3.53); label("Central volume V1 (L)") # Table 2: V1 = 3.53 L (RSE 5.71%); Eq. 3 intercept 1.26
    lq <- log(73.9); label("Inter-compartmental clearance Q2, central to deep peripheral (L/h)") # Table 2: Q2 = 73.9 L/h (RSE 5.45%); Eq. 3 intercept 4.30
    lvp <- log(51.9); label("Peripheral volume V2, deep peripheral (L)") # Table 2: V2 = 51.9 L (RSE 7.49%); Eq. 3 intercept 3.95
    lq2 <- log(28.6); label("Inter-compartmental clearance Q3, central to shallow peripheral (L/h)") # Table 2: Q3 = 28.6 L/h (RSE 5.33%); Eq. 3 intercept 3.35
    lvp2 <- log(419); label("Peripheral volume V3, shallow peripheral (L)") # Table 2: V3 = 419 L (RSE 8.70%); Eq. 3 intercept 6.04

    # Covariate effects, all on the log scale as printed in Equation 3. The
    # categorical effects keep the paper's sign and are applied to the
    # complement of the canonical indicator (see covariateData notes).
    e_wt_cl <- 0.578; label("Power exponent of (WT / 62 kg) on CL (unitless)") # Table 2: 'Impact of body weight on CL' = 0.578 (RSE 11.7%); Eq. 3
    e_age_vc <- 0.496; label("Power exponent of (AGE / 33 years) on V1 (unitless)") # Table 2: 'Impact of age on V1' = 0.496 (RSE 25.7%); Eq. 3
    e_sexf_cl <- 0.089; label("Log-scale effect of female sex on CL (unitless)") # Eq. 3: +0.089 * (sex = female); Table 2 factor 1.09 (RSE 2.52%)
    e_venous_cl <- -0.195; label("Log-scale effect of a venous sample on CL (unitless)") # Eq. 3: -0.195 * (blood sampling site = vein); Table 2 factor 0.823 (RSE 2.44%)
    e_noncritill_cl <- -0.237; label("Log-scale effect of a non-ICU population on CL (unitless)") # Eq. 3: -0.237 * ICU; Table 2 'Impact of ICU on CL' factor 0.789 (RSE 4.36%)
    e_shortinf_cl <- 0.229; label("Log-scale effect of a short-term-infusion study on CL (unitless)") # Eq. 3: +0.229 * FLG; Table 2 factor 1.26 (RSE 3.18%)
    e_shortinf_vp <- -2.08; label("Log-scale effect of a short-term-infusion study on V2 (unitless)") # Eq. 3: -2.08 * FLG; Table 2 factor 0.125 (RSE 8.59%)
    e_shortinf_q2 <- 0.802; label("Log-scale effect of a short-term-infusion study on Q3 (unitless)") # Eq. 3: +0.802 * FLG; Table 2 factor 2.23 (RSE 6.35%)
    e_shortinf_vp2 <- -1.62; label("Log-scale effect of a short-term-infusion study on V3 (unitless)") # Eq. 3: -1.62 * FLG; Table 2 factor 0.198 (RSE 9.62%)

    # Diagonal IIV (control stream $OMEGA has six diagonal elements). Table 2
    # prints IIV in percent; the control stream's $OMEGA values (0.0361,
    # 0.574, 0.249, 0.198, 0.0673, 0.0847) are the square of the printed
    # percentages to within the last digit (sqrt = 19.0, 75.8, 49.9, 44.5,
    # 25.9, 29.1%), so the percentages are omega, the log-scale SD, and the
    # variances below are (percent / 100)^2.
    etalcl ~ 0.191^2 # Table 2: IIV CL = 19.1%, log-scale SD
    etalvc ~ 0.750^2 # Table 2: IIV V1 = 75.0%, log-scale SD
    etalq ~ 0.509^2 # Table 2: IIV Q2 = 50.9%, log-scale SD
    etalvp ~ 0.449^2 # Table 2: IIV V2 = 44.9%, log-scale SD
    etalq2 ~ 0.258^2 # Table 2: IIV Q3 = 25.8%, log-scale SD
    etalvp2 ~ 0.290^2 # Table 2: IIV V3 = 29.0%, log-scale SD

    # Residual error: additive on log-transformed concentrations (Methods
    # Eq. 2; control stream 'Y = IPRED + EPS(1)*(1-SITE) + EPS(2)*SITE' with
    # IPRED = LOG(F)), i.e. log-normal error in nlmixr2. Table 2 labels the
    # two magnitudes 'proportional residual' in percent; on the log scale
    # they are the SD of epsilon, as the $SIGMA values 0.102 and 0.0474
    # (sqrt = 0.319 and 0.218) confirm.
    expSdVenous <- 0.319; label("Log-scale residual SD, venous samples (unitless)") # Table 2: 'Proportional residual of sampling from vein' = 31.9% (RSE 5.33%)
    expSdArterial <- 0.218; label("Log-scale residual SD, arterial samples (unitless)") # Table 2: 'Proportional residual of sampling from arterial' = 21.8% (RSE 7.77%)
  })

  model({
    # Indicators in the paper's own coding (Equation 3 and the control
    # stream): venous sample, non-ICU population, short-term-infusion study.
    venous <- 1 - SAMPLE_ARTERIAL
    noncritill <- 1 - DIS_CRITILL
    shortinf <- 1 - REGIMEN_CONTINF

    # Individual PK parameters (Equation 3).
    cl <- exp(
      lcl +
        etalcl +
        e_wt_cl * log(WT / 62) +
        e_venous_cl * venous +
        e_sexf_cl * SEXF +
        e_shortinf_cl * shortinf +
        e_noncritill_cl * noncritill
    )
    vc <- exp(lvc + etalvc + e_age_vc * log(AGE / 33))
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp + e_shortinf_vp * shortinf)
    q2 <- exp(lq2 + etalq2 + e_shortinf_q2 * shortinf)
    vp2 <- exp(lvp2 + etalvp2 + e_shortinf_vp2 * shortinf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment linear disposition with first-order elimination
    # (NONMEM ADVAN11 TRAN4). Doses are intravenous into central.
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg and V1 in L give mg/L; x 1000 gives ng/mL (control stream
    # scaling S1 = V1 / 1000).
    Cc <- 1000 * central / vc

    expSdSite <- expSdVenous * venous + expSdArterial * SAMPLE_ARTERIAL
    Cc ~ lnorm(expSdSite)
  })
}
