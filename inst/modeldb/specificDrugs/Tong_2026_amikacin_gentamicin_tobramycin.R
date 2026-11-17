Tong_2026_amikacin_gentamicin_tobramycin <- function() {
  description <- paste(
    "Unified two-compartment population PK model for intravenous",
    "amikacin, gentamicin and tobramycin across the human lifespan",
    "(neonates, children and adults; Tong 2026). One jointly fitted",
    "model: clearance, central volume, inter-compartmental clearance",
    "and peripheral volume are shared by the three aminoglycosides, with",
    "drug factors on clearance (tobramycin 1.17, gentamicin 1.11 relative",
    "to amikacin) and on the amikacin inter-compartmental clearance (0.84)",
    "and peripheral volume (1.41). Allometric scaling (exponent 0.75 on",
    "clearances, 1 on volumes, 70 kg reference) uses a size descriptor",
    "that is total body weight below 1 year of age, Al-Sallami fat-free",
    "mass from 3 years, and a linear blend in between. Clearance also",
    "carries sigmoidal postmenstrual-age maturation and decline functions",
    "normalised to a 35-year-old, an exponential effect of serum",
    "creatinine relative to an age-typical value, and a cystic-fibrosis",
    "factor (1.06). IIV on clearance and central volume; combined",
    "proportional and additive residual error with drug-specific",
    "magnitudes."
  )
  reference <- paste(
    "Tong DMH, Faldasz JD, Keizer RJ, Hughes JH. One Model to Dose Them",
    "All: A Population Pharmacokinetic Model for Multiple Aminoglycosides",
    "Across the Human Lifespan. J Clin Pharmacol. 2026;66(9):e70231.",
    "doi:10.1002/jcph.70231. Final estimates in Table 2; NONMEM control",
    "stream in the Supporting Information ('Supplemental Code')."
  )
  vignette <- "Tong_2026_amikacin_gentamicin_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The same states carry whichever aminoglycoside the subject received
  # (TRT_AMIKACIN / TRT_GENTAMICIN / TRT_TOBRAMYCIN). Serum concentrations
  # (Methods, 'Patient Data Collection').
  compartmentData <- list(
    central = list(analyte = "amikacin, gentamicin or tobramycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "amikacin, gentamicin or tobramycin",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  dosing <- "central"

  # One residual magnitude per aminoglycoside (Table 2; control stream
  # SIGMA(1)-SIGMA(6)). The residual-SD suffix names the drug, not a
  # model output, so the canonical propSd_<output> matcher does not
  # recognise it.
  paper_specific_residual_sds <- c(
    "propSd_tobramycin",
    "addSd_tobramycin",
    "propSd_gentamicin",
    "addSd_gentamicin",
    "propSd_amikacin",
    "addSd_amikacin"
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source analysis (Methods, 'Model",
        "Development'). Enters the allometric size descriptor WTallo",
        "directly below 1 year of age and through the Janmahasatian",
        "fat-free mass at older ages (control stream 'WT_USE'). Range",
        "0.3-200 kg, median 55 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    HT = list(
      description = "Body height (length in infants)",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only to form BMI = WT / (HT/100)^2 for the Janmahasatian",
        "fat-free mass that underlies the Al-Sallami (M2) FFM. The source",
        "dataset supplied BMI and the Janmahasatian FFM as precomputed",
        "columns ($INPUT BMI, FFM); this model recomputes both from WT,",
        "HT and SEXF with the equation printed in Methods. Only affects",
        "subjects aged 1 year or older. Range 20-210 cm (Table 1)."
      ),
      source_name = "HT"
    ),
    AGE = list(
      description = "Postnatal (chronological) age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the allometric size descriptor (WT below 1 year, a",
        "linear WT-to-FFM blend from 1 to 3 years, Al-Sallami FFM from 3",
        "years), drives the Al-Sallami age-dependent FFM factor, and gives",
        "postmenstrual age in years as AGE + 40/52 once PAGE >= 60 weeks",
        "(control stream 'PMA_YRS = AGE + (40/52)'). Range 0-100 years",
        "(Table 1)."
      ),
      source_name = "AGE"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only while PAGE < 60 weeks, where postmenstrual age in",
        "years is PAGE / 52 (control stream 'IF(PMA.LT.60) PMA_YRS =",
        "PMA/52'); from 60 weeks the model uses AGE + 40/52 instead.",
        "Gestational age was collected for patients younger than 2 years",
        "and assumed to be 40 weeks for older patients (Methods,",
        "'Patient Data Collection'). Units are weeks because the control",
        "stream's 60-week switch and /52 conversion are written in weeks."
      ),
      source_name = "PMA"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Source column SEX with 1 = male (control stream 'IF (SEX.EQ.1)",
        "THEN ; male'), so SEXF = 1 - SEX. Selects the sex-specific",
        "Janmahasatian and Al-Sallami fat-free-mass equations. The",
        "age-typical serum creatinine in the control stream does not use",
        "sex. Female 3169 / 5972 treatment courses (Table 1)."
      ),
      source_name = "SEX"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Exponential effect on clearance",
        "exp(-0.631 * (CREAT - SCRstd)), where SCRstd is an age-typical",
        "serum creatinine computed from postmenstrual age (Table 2;",
        "control stream 'Fscr'). Median 0.70 mg/dL, range 0.06-13 mg/dL",
        "(Table 1)."
      ),
      source_name = "CR"
    ),
    DIS_CF = list(
      description = "Cystic fibrosis, 1 = present, 0 = absent",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cystic fibrosis)",
      notes = paste(
        "Multiplies clearance by 1.06 (Table 2; control stream",
        "'THETA(12) ** CF'). All 74 patients with cystic fibrosis received",
        "tobramycin (Results, 'Population Pharmacokinetic Model';",
        "Table 1)."
      ),
      source_name = "CF"
    ),
    TRT_TOBRAMYCIN = list(
      description = "Tobramycin indicator, 1 = the subject received tobramycin",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (amikacin or gentamicin)",
      notes = paste(
        "Multiplies clearance by 1.17 relative to amikacin and selects the",
        "tobramycin residual-error magnitudes (Table 2; control stream",
        "'IF (TOBRA.EQ.1) Fdrug = THETA(10)'). Exactly one of",
        "TRT_TOBRAMYCIN, TRT_GENTAMICIN and TRT_AMIKACIN must be 1."
      ),
      source_name = "TOBRA"
    ),
    TRT_GENTAMICIN = list(
      description = "Gentamicin indicator, 1 = the subject received gentamicin",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (amikacin or tobramycin)",
      notes = paste(
        "Multiplies clearance by 1.11 relative to amikacin and selects the",
        "gentamicin residual-error magnitudes (Table 2; control stream",
        "'IF (GENTA.EQ.1) Fdrug = THETA(11)')."
      ),
      source_name = "GENTA"
    ),
    TRT_AMIKACIN = list(
      description = "Amikacin indicator, 1 = the subject received amikacin",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (gentamicin or tobramycin)",
      notes = paste(
        "Amikacin is the reference drug for clearance. The indicator",
        "multiplies the inter-compartmental clearance by 0.84 and the",
        "peripheral volume by 1.41, and selects the amikacin",
        "residual-error magnitudes (Table 2; control stream",
        "'THETA(13) ** AMI', 'THETA(14) ** AMI')."
      ),
      source_name = "AMI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 5658,
    n_treatment_courses = 5972,
    n_studies = 1,
    n_training = 4179,
    age_range = "0-100 years (median 22)",
    weight_range = "0.3-200 kg (median 55)",
    sex_female_pct = 53.1,
    disease_state = paste(
      "Hospitalised patients receiving amikacin, gentamicin or tobramycin",
      "with therapeutic drug monitoring in routine care; 74 treatment",
      "courses in patients with cystic fibrosis."
    ),
    dose_range = "Routine clinical intravenous dosing (doses not tabulated).",
    regions = "United States (154 hospitals)",
    age_groups = "3235 adult, 510 paediatric and 2227 neonatal treatment courses",
    drugs = "552 amikacin, 3972 gentamicin and 1448 tobramycin treatment courses",
    n_observations = 13822,
    notes = paste(
      "Real-world de-identified data from the InsightRX Nova",
      "precision-dosing platform, January 2013 to December 2025 (Results,",
      "'Patients and Data Collection'; Table 1). The model was fitted to a",
      "random 70% training split (4179 treatment courses, 9681",
      "concentrations); 647 of 13,822 concentrations (4.7%) were below the",
      "limit of quantification and were handled with the M3 method."
    )
  )

  ini({
    # Structural parameters at the reference patient: a 70 kg (WTallo)
    # 35-year-old receiving amikacin, without cystic fibrosis and with
    # serum creatinine equal to the age-typical value.
    lcl <- log(5.1); label("Clearance at the reference patient, amikacin (L/h)") # Table 2 theta_CL = 5.1 L/h (RSE 2.3%)
    lvc <- log(28.3); label("Central volume of distribution at 70 kg WTallo (L)") # Table 2 theta_Vc = 28.3 L (RSE 1.2%)
    lq <- log(1.45); label("Inter-compartmental clearance at 70 kg WTallo, gentamicin and tobramycin (L/h)") # Table 2 theta_Q = 1.45 L/h (RSE 4.2%)
    lvp <- log(24.7); label("Peripheral volume of distribution at 70 kg WTallo, gentamicin and tobramycin (L)") # Table 2 theta_Vp = 24.7 L (RSE 4.5%)

    # Allometric exponents, fixed (Results, 'Population Pharmacokinetic
    # Model': 'exponents fixed at 0.75 for clearance and inter-compartmental
    # clearance and 1 for central and peripheral volumes').
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WTallo on CL and Q (unitless)") # Table 2 '(WTallo/70)^0.75'; control stream '** 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of WTallo on Vc and Vp (unitless)") # Table 2 '(WTallo/70)'; control stream 'WT_USE/70'

    # Drug and disease factors on the PK parameters (multiplicative).
    e_trt_tobramycin_cl <- 1.17; label("Tobramycin factor on clearance relative to amikacin (unitless)") # Table 2 'Tobramycin factor on CL' = 1.17 (RSE 2.5%); control stream THETA(10)
    e_trt_gentamicin_cl <- 1.11; label("Gentamicin factor on clearance relative to amikacin (unitless)") # Table 2 'Gentamicin factor on CL' = 1.11 (RSE 2.3%); control stream THETA(11)
    e_dis_cf_cl <- 1.06; label("Cystic-fibrosis factor on clearance (unitless)") # Table 2 'Cystic fibrosis factor on CL' = 1.06 (RSE 3.8%); control stream THETA(12)
    e_trt_amikacin_q <- 0.84; label("Amikacin factor on inter-compartmental clearance (unitless)") # Table 2 'Amikacin factor on Q' = 0.84 (RSE 8.4%); control stream THETA(13)
    e_trt_amikacin_vp <- 1.41; label("Amikacin factor on peripheral volume (unitless)") # Table 2 'Amikacin factor on Vp' = 1.41 (RSE 12.6%); control stream THETA(14)

    # Serum creatinine effect: exp(-e_creat_cl * (CREAT - SCRstd)).
    e_creat_cl <- 0.631; label("Exponential serum-creatinine coefficient on clearance (dL/mg)") # Table 2 equation 'FCR = exp(-0.631*(CR-CRstd))'; parameter row 'Serum creatinine exponent' = 0.63 (RSE 11.2%); control stream EXPCR = EXP(THETA(5))

    # Postmenstrual-age maturation and decline of clearance, both
    # normalised to a 35-year-old (PMA 35 + 40/52 years).
    ltmat50 <- log(1.03); label("Postmenstrual age at 50% clearance maturation (log years)") # Table 2 'Maturation T50' = 1.03 (RSE 59.3%); control stream MAT50 = EXP(THETA(8)); Figure S19 MAT_50 about 1.03
    lhill_mat <- log(2.61); label("Hill exponent of clearance maturation (log unitless)") # Table 2 'Maturation Hill exponent' = 2.61 (RSE 2.9%); control stream GAMM = EXP(THETA(9))
    ltdec50 <- log(66.1); label("Postmenstrual age at 50% age-related decline of clearance (log years)") # Table 2 'Age decay T50' = 66.1 (RSE 0.5%); control stream AGE50 = EXP(THETA(6))
    lhill_dec <- log(2.34); label("Hill exponent of the age-related decline of clearance, magnitude (log unitless)") # Table 2 'Age decay Hill exponent' = -2.34 (RSE 4.7%); control stream GAMMA2 = -EXP(THETA(7)), sign carried by the decline form in model()

    # IIV: Table 2 reports %CV only; variances from omega^2 = log(1 + CV^2).
    # The control stream estimates an OMEGA BLOCK(2), but the covariance
    # is not reported, so it is set to zero here.
    etalcl ~ 0.092352 # Table 2 'IIV on CL' = 31.1 %CV (RSE 4.8%); log(1 + 0.311^2)
    etalvc ~ 0.043964 # Table 2 'IIV on Vc' = 21.2 %CV (RSE 13.8%); log(1 + 0.212^2)

    # Drug-specific combined residual error, Y = IPRED*(1 + EPSprop) +
    # EPSadd (control stream $ERROR).
    propSd_tobramycin <- 0.275; label("Proportional residual SD, tobramycin (fraction)") # Table 2 'PROP TOBRA' = 27.5 %CV (RSE 5.8%)
    addSd_tobramycin <- 0.0648; label("Additive residual SD, tobramycin (mg/L)") # Table 2 'ADD TOBRA' = 0.0648 mg/L (RSE 75.3%)
    propSd_gentamicin <- 0.243; label("Proportional residual SD, gentamicin (fraction)") # Table 2 'PROP GENTA' = 24.3 %CV (RSE 5.6%)
    addSd_gentamicin <- 0.137; label("Additive residual SD, gentamicin (mg/L)") # Table 2 'ADD GENTA' = 0.137 mg/L (RSE 10.7%)
    propSd_amikacin <- 0.197; label("Proportional residual SD, amikacin (fraction)") # Table 2 'PROP AMI' = 19.7 %CV (RSE 14.5%)
    addSd_amikacin <- 0.827; label("Additive residual SD, amikacin (mg/L)") # Table 2 'ADD AMI' = 0.827 mg/L (RSE 21.4%)
  })

  model({
    # Postmenstrual age in years (control stream PMA_YRS): PAGE (weeks)
    # below 60 weeks, otherwise postnatal age plus 40 weeks.
    if (PAGE < 60) {
      pma_yrs <- PAGE / 52
    } else {
      pma_yrs <- AGE + 40 / 52
    }
    pma_ref <- 35 + 40 / 52

    # Maturation and decline, each normalised to the reference
    # postmenstrual age (Table 2; control stream Fmat/Fmat_R and
    # Fage/Fage_R). The control stream writes the decline as
    # pma^g2 / (pma^g2 + tdec50^g2) with g2 = -2.34, which equals
    # 1 - pma^2.34 / (pma^2.34 + tdec50^2.34).
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    tdec50 <- exp(ltdec50)
    hill_dec <- exp(lhill_dec)
    fmat <- pma_yrs^hill_mat / (pma_yrs^hill_mat + tmat50^hill_mat)
    fmat_ref <- pma_ref^hill_mat / (pma_ref^hill_mat + tmat50^hill_mat)
    fdec <- 1 - pma_yrs^hill_dec / (pma_yrs^hill_dec + tdec50^hill_dec)
    fdec_ref <- 1 - pma_ref^hill_dec / (pma_ref^hill_dec + tdec50^hill_dec)

    # Age-typical serum creatinine (mg/dL) from postmenstrual age in years
    # (control stream LnSCRstd; constants are fixed in the published code,
    # with a base-10 logarithm of age).
    creat_std <- exp(-1.22839 + 0.67190 * log10(pma_yrs) + 6.27017 * exp(-3.10940 * pma_yrs))
    f_creat <- exp(-e_creat_cl * (CREAT - creat_std))

    # Allometric size descriptor WTallo (Methods equation; control stream
    # WT_USE). The Al-Sallami (M2) fat-free mass is the Janmahasatian FFM
    # times an age- and sex-dependent factor; the control stream uses 13.2
    # years in the male factor (the Methods equation prints 13.4).
    bmi <- WT / (HT / 100)^2
    ffm_janma <- (1 - SEXF) * 9270 * WT / (6680 + 216 * bmi) +
      SEXF * 9270 * WT / (8780 + 244 * bmi)
    if (AGE < 1) {
      fac_ffm <- 1
    } else {
      fac_ffm <- (1 - SEXF) * (0.88 + (1 - 0.88) / (1 + (AGE / 13.2)^(-12.7))) +
        SEXF * (1.11 + (1 - 1.11) / (1 + (AGE / 7.1)^(-1.1)))
    }
    ffm <- fac_ffm * ffm_janma
    if (AGE < 1) {
      wt_allo <- WT
    } else if (AGE < 3) {
      wt_allo <- (-0.5 * AGE + 1.5) * WT + (0.5 * AGE - 0.5) * ffm
    } else {
      wt_allo <- ffm
    }

    f_drug_cl <- e_trt_tobramycin_cl^TRT_TOBRAMYCIN * e_trt_gentamicin_cl^TRT_GENTAMICIN

    cl <- exp(lcl + etalcl) * (wt_allo / 70)^e_wt_cl * f_creat *
      (fdec / fdec_ref) * (fmat / fmat_ref) * f_drug_cl * e_dis_cf_cl^DIS_CF
    vc <- exp(lvc + etalvc) * (wt_allo / 70)^e_wt_vc
    q <- exp(lq) * (wt_allo / 70)^e_wt_cl * e_trt_amikacin_q^TRT_AMIKACIN
    vp <- exp(lvp) * (wt_allo / 70)^e_wt_vc * e_trt_amikacin_vp^TRT_AMIKACIN

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Drug-specific residual magnitudes (control stream VARPROP / VARADD).
    w_prop <- TRT_TOBRAMYCIN * propSd_tobramycin + TRT_GENTAMICIN * propSd_gentamicin +
      TRT_AMIKACIN * propSd_amikacin
    w_add <- TRT_TOBRAMYCIN * addSd_tobramycin + TRT_GENTAMICIN * addSd_gentamicin +
      TRT_AMIKACIN * addSd_amikacin
    Cc ~ add(w_add) + prop(w_prop) + combined2()
  })
}
