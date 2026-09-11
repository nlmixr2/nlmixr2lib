Jayanti_2026_pyrazinamide <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral pyrazinamide in Korean and Indonesian adults with drug-susceptible tuberculosis (Jayanti 2026); lean body weight is an allometric covariate on CL/F and Vd/F (fixed exponents 0.75 and 1), apparent clearance carries a separate typical value and a separate interindividual variance for each ethnicity, and diabetes mellitus raises CL/F by 23% in Indonesian patients and by 26% in Korean patients aged 60 years or older"
  reference <- paste(
    "Jayanti RP, Cho Y-S, Soedarsono S, Kim H-J, Kang J, Kim J, Oh JY, Kang BH,",
    "Ha JH, Kim J-W, Mertaniasih NM, Kusmiati T, Permatasari A, Yuliwulandari R,",
    "Kim R, Seong H-J, Ghim J-L, Kim D-H, Shin J-G; on behalf of the cPMTb.",
    "Population pharmacokinetics model of pyrazinamide to optimize tuberculosis",
    "treatment: An interethnic cohort study of diabetes mellitus effect on drug",
    "exposure. PLoS One. 2026;21(1):e0340133. doi:10.1371/journal.pone.0340133.",
    "Correction: PLoS One. 2026;21(4):e0347490. doi:10.1371/journal.pone.0347490",
    "(corrects the funding statement only; no model parameter is affected)."
  )
  vignette <- "Jayanti_2026_pyrazinamide"

  # Table 2 and S2 Table report concentrations in mg/L and exposures in
  # mg.h/L; the bioanalytical section quotes the same calibration range as
  # 2.0-80.0 mg/L and the LLOQ as 2.0 ug/mL. mg/L and ug/mL are the same unit,
  # so no conversion is involved anywhere in this file.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Apparent clearance carries TWO separate interindividual variances, one per
  # ethnicity (Table 2 rows 'omega^2; CL/F Indonesian (%)' and 'omega^2; CL/F
  # Korean (%)'), selected by RACE_KOREAN rather than a single eta modulated by
  # a covariate. Neither name matches the bare canonical `etalcl`, so both are
  # declared here. Same mechanism and same rationale as the cohort-stratified
  # clearance etas of Du_2025_repotrectinib.R.
  paper_specific_etas <- c("etalclIndonesian", "etalclKorean")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrazinamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor on both CL/F and Vd/F with fixed exponents 0.75 and 1 respectively, normalised to a 45 kg reference subject. The 45 kg reference is stated in the Table 2 footnote d: 'Allometric scaling was applied to the CL/F and Vd/F data, and typical values reported here refer to the typical patient, with lean body weight of 45 kg.' (That footnote is present in the publisher's rendered Table 2 but is dropped by PDF-to-markdown conversion, so it must be read off the table image.) It is the rounded population median lean body weight, reported as 45.5 kg in Table 1 and computed as 45.56 kg from the deposited dataset (S2 File). Lean body weight rather than total body weight was chosen because it 'described the size effect sufficiently and significantly reduced the OFV' (Results, 'Population pharmacokinetics model of PZA'). Derived per subject from the Boer equation given in Methods: LBW(male) = 0.407 * WT + 0.267 * HT - 19.2 and LBW(female) = 0.252 * WT + 0.473 * HT - 48.3, with WT in kg and HT in cm; reproducing that formula against the WEIGHT / HEIGHT / SEX columns of the deposited dataset recovers the LBW column to within 0.63 kg, which also establishes that SEX = 1 codes male in that file.",
      source_name        = "LBW"
    ),
    RACE_KOREAN = list(
      description        = "Korean-heritage indicator (1 = Korean, 0 = Indonesian).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Indonesian)",
      notes              = "The cohort contains exactly two ethnicities, matched 1:1 on body weight and age (160 Korean and 160 Indonesian patients; Table 1), so this single indicator fully partitions the population and the reference category is specifically Indonesian rather than the generic 'non-Korean' of the register entry. It is load-bearing in three distinct places in the model, which is why it is a structural stratifier rather than an ordinary covariate effect. (1) It selects which of two separately estimated typical apparent clearances applies (Table 2 'CL/F Indonesian; = theta1' = 3.18 L/h versus 'CL/F Korean; = theta2' = 3.5 L/h). (2) It selects which of two separately estimated clearance interindividual variances applies (Table 2 'omega^2; CL/F Indonesian' versus 'omega^2; CL/F Korean'). (3) It selects which diabetes covariate form applies: plain diabetes in Indonesians versus diabetes crossed with age >= 60 years in Koreans. The authors stress that ethnicity itself is NOT a covariate effect here: 'Ethnicity did not significantly affect the CL/F and Vd/F of PZA' and 'we did not find inter-ethnicity differences in the PK of PZA'. The per-ethnicity split exists because the age distribution of diabetes differs between the two countries ('The Indonesian TB patients with DM tended to be < 60 years old. Most Korean TB-DM patients were older adults'), so a single pooled diabetes term could not describe both. In the deposited dataset (S2 File) the source column Ethnic codes 1 = Indonesian and 2 = Korean, so RACE_KOREAN = Ethnic - 1.",
      source_name        = "Ethnic"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator (1 = diabetes mellitus, 0 = no diabetes mellitus).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes mellitus)",
      notes              = "Diabetes mellitus recorded as a comorbidity at study entry; 77 of the 320 patients (24.1%) were diabetic, 55 of 160 Indonesians and 22 of 160 Koreans (Table 1). Type 1 versus type 2 is not distinguished. In Indonesian patients this column enters the model as a main effect on CL/F; in Korean patients it enters ONLY through the product AGE_GE60 * DIS_DIAB, so a Korean diabetic patient younger than 60 takes the Korean reference clearance unmodified. Jayanti 2026 attributes the higher apparent clearance to diabetes-induced elevation of xanthine oxidase, the enzyme that converts pyrazinamide to 5-hydroxypyrazinoic acid (Discussion), the same mechanism cited by the group's earlier Korean-only model (Kim_2023_pyrazinamide.R).",
      source_name        = "DM"
    ),
    AGE_GE60 = list(
      description        = "Advanced-age indicator (1 = aged 60 years or older at baseline, 0 = younger than 60 years).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (younger than 60 years)",
      notes              = "Enters the model ONLY through the product AGE_GE60 * DIS_DIAB ('old DM'), never as a main effect; 25 of the 320 patients (7.8%) were both >= 60 years old and diabetic (Table 1). The paper is internally inconsistent about whether the cut point is strict: Table 1's row header and the S3 / S4 Table footnotes print 'Age > 60 years old', while Methods ('geriatric (>= 60 years) with DM'), Results ('older patients (>= 60 years old) with DM') and the Fig 3 caption ('old DM: >= 60 years old patients with DM, Other patients: patients who aged < 60 years old') print the inclusive form. The deposited dataset (S2 File) settles it: the source column OLD spans ages 17-59 where it is 0 and 60-80 where it is 1, so the threshold is inclusive and AGE_GE60 (not AGE_GT60) is the correct canonical. The dataset also confirms OLDDM = OLD * DM exactly for all 300 deposited subjects.",
      source_name        = "OLD"
    )
  )

  # Covariates that Jayanti 2026 screened but did not retain in the final
  # model, plus the three demographics that are inputs to the lean-body-weight
  # derivation rather than model terms. Documentation only -- none of these is
  # referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as an allometric size descriptor but rejected in favour of lean body weight, which gave the better fit and stability (Methods; Results). Still required upstream as an input to the Boer lean-body-weight equation. Cohort median 55 kg (IQR 50-59.8); Korean median 56 kg, Indonesian median 50 kg (Table 1). Body weight was one of the two variables matched between ethnicities, to a maximum difference of 5 kg."
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Not screened as a covariate in its own right; required as the second input to the Boer lean-body-weight equation. Cohort median 165 cm (IQR 160-170) and identical in both ethnicities (Table 1)."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate and not retained. Sex was deliberately NOT matched between ethnicities: the authors argue that 'most of the effect of sex on PK is a result of the different body fat compositions of males and females' and that carrying lean body weight in the model already absorbs it, 'thereby addressing sex-based variability in PZA PK' (Discussion). Sex is nonetheless required as an input to the Boer equation, which uses different coefficients for males and females. Cohort was 38.4% female (Table 1). The deposited dataset codes SEX = 1 for male, so SEXF = 1 - SEX."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate and not retained; kept in the model only as the >= 60 year threshold crossed with diabetes mellitus (see AGE_GE60 in covariateData). Cohort median 46 years (IQR 31-57); range 17-80 years in the deposited dataset. Age was the second of the two variables matched between ethnicities, to a maximum difference of 5 years."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL (as reported in Table 1; the canonical SI unit is g/L)",
      type        = "continuous",
      notes       = "Significant in univariate analysis (p < 0.05) but 'excluded due to large standard errors' (Results), so it is absent from the final model. Cohort median 3.6 g/dL (IQR 3.0-4.1), with a marked between-country difference (Korean 4.1, Indonesian 3.0) that the paper does not model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Significant in univariate analysis (p < 0.05) but 'excluded due to large standard errors' (Results), so it is absent from the final model. Cohort median 27 U/L (IQR 23-36.75)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not significant and not retained. Cohort median 24 U/L (IQR 18-33)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL (as reported in Table 1; the canonical SI unit is umol/L)",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not significant and not retained. Cohort median 0.52 mg/dL (IQR 0.4-0.69)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not significant and not retained. Cohort median 10.4 mg/dL (IQR 8.4-13.3)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, computed with the CKD-EPI 2009 equation.",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not significant and not retained. The cohort had essentially normal renal function: median 105.5 mL/min/1.73 m^2 (IQR 97.7-121.8), so the data carry little information about renal impairment. Serum creatinine was screened separately (median 0.7 mg/dL) and likewise not retained. Liver disease was screened as a categorical covariate and not retained; it affected only 7 of 320 patients (2.2%) and has no canonical covariate column in this library."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 320L,
    n_observations = 407L,
    n_studies      = 1L,
    n_sites        = 23L,
    age_range      = "17-80 years (deposited dataset); Table 1 reports median 46 years, IQR 31-57",
    age_median     = "46 years",
    weight_range   = "Table 1 reports median 55 kg, IQR 50-59.8 kg",
    weight_median  = "55 kg",
    lbm_median     = "45.5 kg (IQR 40.8-49.3)",
    height_median  = "165 cm (IQR 160-170)",
    sex_female_pct = 38.4,
    race_ethnicity = c(Korean = 50, Indonesian = 50),
    disease_state  = "Adults (>= 18 years) with drug-susceptible tuberculosis receiving a pyrazinamide-based regimen for at least 2 weeks. Comorbidities: 77 of 320 (24.1%) diabetes mellitus, of whom 25 (7.8% of the cohort) were also >= 60 years old; 7 of 320 (2.2%) liver impairment. Non-adherent, pregnant and multi-drug-resistant patients were excluded.",
    dose_range     = "Oral pyrazinamide once daily as prescribed under WHO weight-band guidance: 500 mg (0.6%), 750 mg (2.2%), 800 mg (0.3%), 1,000 mg (25.9%), 1,200 mg (6.9%), 1,250 mg (12.2%), 1,275 mg (0.6%), 1,500 mg (42.8%), 1,600 mg (7.5%) and 2,000 mg (0.9%). Sampling was at steady state; the deposited dataset codes every dose record as SS = 1 with II = 24 h.",
    regions        = "Republic of Korea (22 hospitals) and Indonesia (one hospital, Surabaya); the cPMTb multinational prospective tuberculosis cohort, patients recruited 2018-08-20 to 2021-11-03.",
    renal_function = "Essentially normal: median eGFR 105.5 mL/min/1.73 m^2 (IQR 97.7-121.8).",
    co_medication  = "Companion anti-tuberculosis regimens, chiefly RHZE (88.1%), plus RHZES, RHZ, RZEM and HZEL. Interactions among the anti-tuberculosis drugs were not modelled. 77 patients received antidiabetic co-medication (insulin, biguanide, sulfonylurea, DPP-IV inhibitor or SGLT2 inhibitor).",
    notes          = "The two ethnic groups were matched 1:1 on body weight and age (maximum difference 5 kg and 5 years), using the smaller Indonesian sample as the reference, to remove demographic confounding from the interethnic comparison; 160 patients per ethnicity contributed 407 concentrations (240 Indonesian, 167 Korean). Blood was drawn at random times 0-24 h after the last dose: one sample per outpatient, at least two per inpatient. Pyrazinamide was quantified by validated HPLC-ESI-MS/MS over 2.0-80.0 mg/L (LLOQ 2.0 mg/L); 117 below-LLOQ samples and 28 outlier patients were excluded, the authors having verified by refitting that including the below-LLOQ records worsened the goodness-of-fit diagnostics and parameter precision. Baseline demographics from Table 1; final parameter estimates and the 45 kg allometric reference from Table 2 including its footnote d; post-hoc exposure summaries from S2 Table; probability-of-target-attainment results from S3 and S4 Tables; covariate-model selection from S1 Table; subject-level covariate coding verified against the deposited NONMEM dataset (S2 File), which contains 300 of the 320 modelled subjects and 404 of the 407 concentrations."
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters. Jayanti 2026 fits pyrazinamide with a
    # one-compartment model with first-order absorption and first-order
    # elimination: "A one-compartment model with first-order
    # absorption-elimination with additive residual error adequately
    # described the PK of PZA" (Results). Lag-time, sequential zero- then
    # first-order, and transit-compartment absorption models were all
    # evaluated and rejected -- "None of the absorption models evaluated had
    # improved performance" (Results) -- so plain first-order absorption is
    # the final structure.
    #
    # ALL typical values below refer to a patient with a lean body weight of
    # 45 kg. That reference is given only in the Table 2 footnote d:
    # "Allometric scaling was applied to the CL/F and Vd/F data, and typical
    # values reported here refer to the typical patient, with lean body
    # weight of 45 kg." It is the rounded population median (45.5 kg in
    # Table 1; 45.56 kg recomputed from the deposited dataset).
    #
    # Apparent clearance has TWO typical values, one per ethnicity, because
    # "Model performance was significantly improved by estimating the
    # population parameters of CL/F for each ethnicity separately" (Results).
    # This is a structural split, not an ethnicity covariate effect: the
    # paper is explicit that "Ethnicity did not significantly affect the
    # CL/F and Vd/F of PZA".
    # ---------------------------------------------------------------
    lclIndonesian <- log(3.18); label("Apparent oral clearance for a 45 kg lean-body-weight Indonesian patient without diabetes mellitus (L/h)")  # Jayanti 2026 Table 2 'CL/F Indonesian; = theta1' = 3.18 (RSE 4.7%; bootstrap median 3.18, 95% CI 2.9-3.5)
    lclKorean     <- log(3.5);  label("Apparent oral clearance for a 45 kg lean-body-weight Korean patient who is not both diabetic and aged 60 years or older (L/h)")  # Jayanti 2026 Table 2 'CL/F Korean; = theta2' = 3.5 (RSE 4.5%; bootstrap median 3.5, 95% CI 3.2-3.8)
    lvc           <- log(52.8); label("Apparent volume of distribution for a 45 kg lean-body-weight patient (L)")  # Jayanti 2026 Table 2 'Vd/F (L) = theta5' = 52.8 (RSE 6.6%; bootstrap median 52.7, 95% CI 46.2-60.6)
    lka           <- log(2.0);  label("Absorption rate constant (1/h)")  # Jayanti 2026 Table 2 'Ka (h-1) = theta6' = 2.0 (RSE 15%; bootstrap median 2.05, 95% CI 1.3-2.6)

    # Allometric exponents on lean body weight. Both are FIXED to the
    # canonical Anderson and Holford (2008) values, not estimated:
    # "LBW-based scaling using fixed exponents of 0.75 for CL/F and 1 for
    # Vd/F was incorporated due to significant improvement in model fit and
    # stability" (Methods). Consistent with that, neither appears as a theta
    # row in Table 2 and neither carries an %RSE.
    e_lbm_cl <- fixed(0.75); label("Allometric exponent of lean body weight on CL/F (unitless)")  # Jayanti 2026 Methods, 'fixed exponents of 0.75 for CL/F'
    e_lbm_vc <- fixed(1.0);  label("Allometric exponent of lean body weight on Vd/F (unitless)")  # Jayanti 2026 Methods, 'and 1 for Vd/F'

    # ---------------------------------------------------------------
    # Diabetes effects on CL/F. Both are encoded as a fractional increase
    # from 1, matching the categorical covariate form given in Methods,
    # P = theta_i * (1 + theta_(i+1) * Cat_Cov), and written out per row in
    # Table 2 as "CL/F; DM = theta1 x (1 + theta3)" and
    # "CL/F; Old DM = theta2 x (1 + theta4)". Note which base theta each row
    # multiplies: theta3 acts on the INDONESIAN typical value and theta4 on
    # the KOREAN one. The Results state the retained combination explicitly:
    # "we included DM (in Indonesians) - OldDM (in Koreans) as covariates of
    # CL/F in the final model."
    #
    # Arithmetic check against the Results narrative: 3.18 * 1.23 = 3.91 and
    # 3.5 * 1.26 = 4.41, against the quoted "CL/F estimates for Indonesian
    # and older Korean patients with DM were 3.88 and 4.38 L/h" -- agreeing
    # to 0.8% and 0.7%, the rounding of the two-significant-figure table.
    # ---------------------------------------------------------------
    e_diab_cl          <- 0.23; label("Fractional increase in CL/F for Indonesian patients with diabetes mellitus (unitless)")  # Jayanti 2026 Table 2 'CL/F; DM = theta1 x (1 + theta3)', theta3 = 0.23 (RSE 41.6%; bootstrap 95% CI 0.07-0.48)
    e_age_ge60_diab_cl <- 0.26; label("Fractional increase in CL/F for Korean patients aged 60 years or older with diabetes mellitus (unitless)")  # Jayanti 2026 Table 2 'CL/F; Old DM = theta2 x (1 + theta4)', theta4 = 0.26 (RSE 42.3%; bootstrap 95% CI 0.06-0.5)

    # ---------------------------------------------------------------
    # Interindividual variability. Table 2 reports these rows as
    # 'omega^2; <parameter> (%)', and footnote a defines omega^2 as the
    # "variance of interindividual variability", so the tabulated number is
    # the variance multiplied by 100:
    #   CL/F Indonesian: 20.5 -> 0.205 (RSE 9.7%,  shrinkage 35.9%)
    #   CL/F Korean:      5   -> 0.05  (RSE 19.6%, shrinkage 56.9%)
    #   Vd/F:            23.6 -> 0.236 (RSE 12.7%, shrinkage 24.6%)
    #   Ka:               0   -> fixed to zero, see below
    # The same cPMTb group uses the identical table convention in
    # Kim_2023_pyrazinamide.R, where it is confirmed against the literal
    # '(0.03 FIX)' of the deposited NONMEM control stream. It is also the
    # only reading consistent with this paper's own post-hoc spread: S2
    # Table gives the Indonesian CL/F interquartile range as 2.61-4.44 L/h
    # about a 3.23 L/h median, i.e. an empirical-Bayes log-scale spread of
    # about 0.41, which a variance of 0.205 (omega = 0.45) can shrink down to
    # but a variance of 0.0413 (omega = 0.203, the reading in which 20.5 is
    # a CV%) cannot reach at all.
    #
    # Clearance variability is nearly four times larger in the Indonesian
    # arm; the Discussion notes the matching observation that "the CL/F and
    # AUC 0-24 ranges were wider in Indonesian than in Korean patients". The
    # Korean shrinkage of 56.9% reflects the one-sample-per-outpatient
    # design (167 concentrations from 160 Korean patients, against 240 from
    # 160 Indonesian patients).
    #
    # IIV on Ka is reported as "0 (FIX)" and is therefore NOT represented by
    # an eta here. Declaring `etalka ~ fixed(0)` would put a structural zero
    # on the OMEGA diagonal and make the matrix singular, which breaks
    # simulation at the Cholesky factorisation; omitting the eta is the
    # faithful encoding of a variance fixed to zero.
    # ---------------------------------------------------------------
    etalclIndonesian ~ 0.205; label("IIV on apparent clearance, Indonesian patients (variance, log scale)")  # Jayanti 2026 Table 2 'omega^2; CL/F Indonesian (%)' = 20.5 -> 0.205
    etalclKorean     ~ 0.05;  label("IIV on apparent clearance, Korean patients (variance, log scale)")      # Jayanti 2026 Table 2 'omega^2; CL/F Korean (%)' = 5 -> 0.05
    etalvc           ~ 0.236; label("IIV on apparent volume of distribution (variance, log scale)")          # Jayanti 2026 Table 2 'omega^2; Vd/F (%)' = 23.6 -> 0.236

    # Residual error. The final model is additive only, on the mg/L scale of
    # the assay: Table 2 lists a single "Additive" residual row and the
    # Results describe the model as having "additive residual error". No
    # proportional component is reported. Additive, proportional and combined
    # error models were all tested (Methods). The value is a standard
    # deviation, following the same group's Kim_2023_pyrazinamide.R, where
    # the deposited control stream writes W = SQRT(THETA(4)**2 + ...) with
    # $SIGMA 1 FIX, i.e. the tabulated residual theta is an SD in
    # concentration units.
    addSd <- 1.29; label("Additive residual standard deviation (mg/L)")  # Jayanti 2026 Table 2 'Residual variability / Additive' = 1.29 (RSE 14.6%; bootstrap median 1.28, 95% CI 0.92-1.7)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Ethnicity and composite covariate indicators.
    #
    # RACE_KOREAN partitions the cohort completely (160 Korean, 160
    # Indonesian), so exactly one of isKorean / isIndonesian is 1 for any
    # subject. The deposited dataset (S2 File) carries the composite
    # geriatric-diabetes flag as its own column OLDDM; it is reconstructed
    # here as the product of the two canonical columns, which reproduces
    # that column exactly for all 300 deposited subjects.
    # ---------------------------------------------------------------
    isKorean     <- RACE_KOREAN
    isIndonesian <- 1 - RACE_KOREAN
    oldDiab      <- AGE_GE60 * DIS_DIAB

    # ---------------------------------------------------------------
    # 2. Ethnicity-stratified clearance terms.
    #
    # Typical value, interindividual variability and the diabetes covariate
    # form are each selected by ethnicity (Table 2; Results). Written out:
    #
    #   Indonesian: CL/F = 3.18 * (LBW/45)^0.75 * (1 + 0.23 * DM)
    #   Korean:     CL/F = 3.5  * (LBW/45)^0.75 * (1 + 0.26 * OLD * DM)
    #
    # A Korean diabetic patient younger than 60 therefore takes the Korean
    # reference clearance unmodified, which is what the Results describe as
    # "3.5 L/h for younger Korean patients with DM"; an Indonesian diabetic
    # patient takes the 23% increase at any age.
    # ---------------------------------------------------------------
    lclSel    <- lclIndonesian * isIndonesian + lclKorean * isKorean
    etalclSel <- etalclIndonesian * isIndonesian + etalclKorean * isKorean
    diabCl    <- (1 + e_diab_cl * DIS_DIAB) * isIndonesian +
      (1 + e_age_ge60_diab_cl * oldDiab) * isKorean

    # ---------------------------------------------------------------
    # 3. Individual parameters. Allometric scaling on lean body weight is
    #    applied to both disposition parameters, normalised to the 45 kg
    #    typical patient of the Table 2 footnote d.
    # ---------------------------------------------------------------
    cl <- exp(lclSel + etalclSel) * (LBM / 45)^e_lbm_cl * diabCl
    vc <- exp(lvc + etalvc) * (LBM / 45)^e_lbm_vc
    ka <- exp(lka)

    kel <- cl / vc

    # ODEs: one-compartment disposition with first-order absorption,
    # equivalent to NONMEM ADVAN2 TRANS2 with S2 = V and K = CL/V.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. CL/F and Vd/F are apparent (oral) parameters, so
    # bioavailability is absorbed into them and no f(depot) term is needed.
    # Dose in mg divided by volume in L gives mg/L, the assay unit.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
