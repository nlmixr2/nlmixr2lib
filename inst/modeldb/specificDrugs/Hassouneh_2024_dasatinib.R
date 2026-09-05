Hassouneh_2024_dasatinib <- function() {
  description <- paste(
    "Two-compartment population PK model with Savic transit-compartment",
    "absorption for oral dasatinib (SPRYCEL) 140 mg in healthy Middle",
    "Eastern male volunteers (Hassouneh 2024). An analytical transit chain",
    "(Ktr = 18.8 1/h, Mtt = 0.48 h, derived N = 8.02) feeds first-order",
    "absorption (Ka) into a two-compartment disposition model. Body mass",
    "index is the only retained covariate and scales Ka with a power",
    "exponent of -0.85, so absorption slows as BMI rises. Interindividual",
    "variability on Ktr, Mtt, Ka, and CL; interoccasion variability across",
    "the two SPRYCEL dosing occasions on the same four parameters.",
    "Combined additive-plus-proportional residual error. Data come from the",
    "two reference-product periods of a four-period full-replicate",
    "bioequivalence study; CL, V1, Q, and V2 are apparent (/F) quantities.",
    sep = " "
  )
  reference <- paste(
    "Hassouneh WB, Al-Ghazawi MA, Saleh MI, Najib N. Population",
    "Pharmacokinetics of Dasatinib in Healthy Subjects. Pharmaceuticals",
    "(Basel). 2024 May 23;17(6):671. doi:10.3390/ph17060671.",
    sep = " "
  )
  vignette <- "Hassouneh_2024_dasatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the Hassouneh 2024 final model",
        "(Results 2.2: 'Body mass index (BMI) index was the only identified",
        "covariate associated with the absorption rate constant'). Table 3",
        "reports the coefficient as beta_BMI = -0.85 (RSE 36.9%), i.e.",
        "absorption slows as BMI rises, which the Discussion frames as",
        "'the higher value of BMI, the lower the Ka value'.",
        "",
        "CENTERING ERRATUM. The equation printed immediately above Table 3",
        "is log(Ka) = log(0.37) - 0.85*log(BMI), with NO centering value",
        "inside the logarithm. Taken literally that makes Ka = 0.37 the",
        "absorption rate at BMI = 1 kg/m^2 and the typical-subject Ka only",
        "0.37 * 22.9^-0.85 = 0.026 1/h. That reading is falsified by the",
        "paper's own visual predictive check (Figure 5): it predicts a",
        "typical 140 mg Cmax of 10.7 ng/mL, whereas the observed and",
        "simulated median profiles in Figure 5 both peak near 105-130 ng/mL",
        "at about 1-1.5 h. The centered reading, Ka = 0.37 * (BMI/22.9)^-0.85,",
        "reproduces Figure 5 closely (typical Cmax 131 ng/mL at Tmax 0.90 h),",
        "so the printed equation is read as having dropped the centering",
        "term, which is Monolix's standard log-transformed-covariate form",
        "log(BMI/BMI_ref). The reference is taken as the cohort MEDIAN BMI",
        "of 22.9 kg/m^2 (Table 2). The cohort MEAN of 23.6 kg/m^2 is the",
        "other plausible Monolix default; because the covariate enters as a",
        "power, that choice only rescales the typical Ka by",
        "(23.6/22.9)^0.85 = 2.6%, which is immaterial next to the 14-fold",
        "difference between the centered and uncentered readings. The",
        "vignette quantifies both comparisons.",
        "",
        "A coefficient on log(BMI/BMI_ref) added to a log-scale parameter is",
        "algebraically identical to the power form (BMI/BMI_ref)^beta used in",
        "model(), which is the form standard across nlmixr2lib.",
        sep = " "
      ),
      source_name        = "BMI"
    ),
    OCC = list(
      description        = "Dosing-occasion indicator: 1 = first SPRYCEL period, 2 = second SPRYCEL period",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "The source study is a four-period full-replicate bioequivalence",
        "study in which 'in two periods of the study, the subjects were",
        "administered an oral dose of SPRYCEL 140 mg film coated tablet, the",
        "reference product' (Methods 4.1). Only those two reference periods",
        "enter the popPK dataset -- the abstract's '4180 plasma observations",
        "from 110 subjects who were administered SPRYCEL on two separate",
        "occasions' is exactly 110 x 2 x 19 samples. Each subject therefore",
        "contributes two occasions, and Results 2.1 notes 'this separate dual",
        "administration raises the ability to investigate the IOV'.",
        "Decomposed inside model() into binary indicators oc1 and oc2 that",
        "multiplex the per-occasion IOV etas on log-Ktr, log-Mtt, log-Ka, and",
        "log-CL. Table 3 reports a single IOV standard deviation per",
        "parameter (not one per occasion), so the occasion-2 variances are",
        "fixed equal to the occasion-1 variances, mirroring the NONMEM",
        "'$OMEGA BLOCK(1) SAME' idiom already used by",
        "Jiang_2024_empagliflozin.R and Chen_2023_nemonoxacin.R. Both",
        "occasions administer the identical product at the identical dose,",
        "so occasion enters only through the IOV random effects.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates screened by Hassouneh 2024 but NOT retained in the final model.
  # Methods 4.3.2 and Discussion list the full screen: "subject demographics
  # (age; body weight; body mass index; smoking status), laboratory tests
  # (total bilirubin, ALP, AST, ALT, Cr and BUN, glucose, white blood cell
  # count, red blood cells count, platelets count, neutrophils, lymphocytes
  # and hemoglobin level); and concurrent medications (such as paracetamol
  # and diclofenac)". Only BMI survived COSSAC forward inclusion and backward
  # elimination, so the rest are recorded here as documentation and are
  # deliberately never referenced in model(). Hassouneh 2024 publishes no
  # point estimate for any of them, so none can be encoded even optionally.
  # The paper reports no numeric summary for the laboratory covariates; the
  # inclusion criteria (Methods 4.1) required every laboratory value to lie
  # "within laboratory reference ranges", which the Conclusions call out as a
  # limitation ("the data used for analysis were taken from healthy
  # volunteers, who rendered the data of many covariates within the reference
  # ranges").
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. Table 2: median 33 years (range 18-49), mean 32 +/- 8.36. Study eligibility 18-55 years."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. Table 2: median 70 kg (range 51-100), mean 72 +/- 12.71. Note that BMI -- which combines weight and height -- WAS retained, on Ka."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Not listed among the screened covariates, but reported in Table 2 (median 175 cm, range 160-192, mean 175 +/- 5.86) and load-bearing because BMI = WT / (HT/100)^2. Recorded so a data assembler can reconstruct the retained BMI covariate."
    ),
    SMOKE = list(
      description        = "Current-smoker indicator (1 = smoker, 0 = non-smoker)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker)",
      notes              = paste(
        "Screened (Methods 4.3.2), not retained. Table 2: 85.34% smokers,",
        "14.66% non-smokers. Methods 4.3.2 defines the positive level",
        "narrowly -- 'smoker subject was a subject who smoked more than 10",
        "cigarettes per day' -- yet exclusion criterion 4 excluded any",
        "'heavy smoker (more than 10 cigarettes per day)'. The two",
        "statements cannot both hold, so the operational definition behind",
        "the 85.34% figure is ambiguous; see the vignette Errata. Immaterial",
        "to the packaged model because the covariate was not retained.",
        sep = " "
      )
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as 'total bilirubin' (Methods 4.3.2), not retained. No numeric summary published; inclusion required values within laboratory reference ranges. Units not stated by the source; mg/dL recorded as the register's default assay unit."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published. Inclusion criterion 7 accepted ALP below the reference range."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required values within laboratory reference ranges."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required values within laboratory reference ranges."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as 'Cr' (Methods 4.3.2), not retained. No numeric summary published. Inclusion criterion 7 accepted creatinine below the reference range. Units not stated by the source."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required values within laboratory reference ranges."
    ),
    FPG = list(
      description = "Fasting plasma glucose (screening clinical chemistry)",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as 'glucose' (Methods 4.3.2), not retained. No numeric summary published. The source says only 'glucose'; screening chemistry in a fasting bioequivalence study is a fasting sample, so it is filed under the FPG canonical."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^12/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    NEUT = list(
      description = "Absolute neutrophil count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    LYMPH_ABS = list(
      description = "Absolute peripheral-blood lymphocyte count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened as 'lymphocytes' (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened as 'hemoglobin level' (Methods 4.3.2), not retained. No numeric summary published; inclusion required hematology within 5% of the reference limits."
    ),
    CONMED_PARACETAMOL = list(
      description        = "Concomitant paracetamol (acetaminophen) administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant paracetamol)",
      notes              = "Screened (Methods 4.3.2: 'concurrent medications (such as paracetamol and diclofenac)'), not retained. No numeric summary or point estimate published. Concomitant medications in this study were rescue treatments for adverse events, not scheduled therapy."
    ),
    CONMED_DICLOFENAC = list(
      description        = "Concomitant diclofenac administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant diclofenac)",
      notes              = "Screened (Methods 4.3.2), not retained. No numeric summary or point estimate published. Concomitant medications in this study were rescue treatments for adverse events, not scheduled therapy."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "dasatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dasatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dasatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species         = "human",
    n_subjects      = 110L,
    n_studies       = 1L,
    n_observations  = 4180L,
    age_range       = "18-49 years (Table 2: median 33, mean 32 +/- 8.36); study eligibility 18-55 years",
    weight_range    = "51-100 kg (Table 2: median 70 kg, mean 72 +/- 12.71)",
    height_range    = "160-192 cm (Table 2: median 175 cm, mean 175 +/- 5.86)",
    bmi_range       = "18.6-29.8 kg/m^2 (Table 2: median 22.9, mean 23.6 +/- 3.50); study eligibility 18.5-30.0 kg/m^2. About 35% of subjects were overweight or obese (BMI > 25 kg/m^2) per the Discussion",
    sex_female_pct  = 0,
    race_ethnicity  = c(`Middle Eastern` = 100),
    disease_state   = "Healthy volunteers (medical history, physical examination, 12-lead ECG and laboratory investigations within acceptable limits; no hepatic, renal, cardiovascular, gastrointestinal or hematological disease)",
    dose_range      = "Single oral dose of SPRYCEL 140 mg film-coated tablet under fasting conditions with 240 mL of water, given on each of two separate occasions (one tablet per occasion); no food until 4 h post-dose",
    smoking_status  = "85.34% smokers, 14.66% non-smokers (Table 2); exclusion criterion 4 barred heavy smokers (more than 10 cigarettes per day)",
    regions         = "Jordan (single centre, International Pharmaceutical Research Center, Amman)",
    co_medication   = "None scheduled; subjects agreed to take no prescription or non-prescription drugs with systemic absorption for at least two weeks before the first dose, and no CYP3A4-affecting medication or food, grapefruit, alcohol or methylxanthines during the study. Paracetamol and diclofenac appear only as rescue medication and were screened as covariates without being retained.",
    notes           = paste(
      "116 healthy subjects enrolled, 110 completed and contributed the",
      "4180 plasma dasatinib observations (110 subjects x 2 SPRYCEL",
      "occasions x 19 samples). The dataset was split before analysis: 88",
      "subjects (3344 observations) were used to estimate the model and 22",
      "subjects (20%) were held out for internal validation. Samples were",
      "drawn pre-dose and at 0.167, 0.333, 0.50, 0.667, 1.00, 1.33, 1.67,",
      "2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 12.00, 16.00 and",
      "24.00 h. Plasma dasatinib was measured by LC-MS/MS with a",
      "dasatinib-d8 internal standard, validated over 0.50-500.00 ng/mL",
      "(r > 0.99). The model was built in Monolix 2020R1 by SAEM, with the",
      "base structure (one vs two vs three compartments; no-delay vs lag",
      "time vs transit chain; with vs without IOV) chosen on the corrected",
      "Bayesian information criterion, and covariates screened by COSSAC",
      "(conditional sampling for stepwise approach based on correlation",
      "tests) with backward elimination at p = 0.01. Evaluation was by",
      "visual predictive check (1000 simulated datasets) against three",
      "datasets: the estimation data (Figure 5), the 22-subject internal",
      "hold-out (Figure 6), and an external dataset of 90 subjects from a",
      "separate bioequivalence study (Figure 7). The source parent study is",
      "a four-period full-replicate bioequivalence trial; only the two",
      "reference-product (SPRYCEL) periods enter this analysis, so the",
      "packaged model describes the originator product.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are FINAL ESTIMATES from Hassouneh
    # 2024 Table 3 ("Summary of final model parameters"), column "Population
    # Parameter (%RSE)". Every %RSE is below 10% except the BMI coefficient
    # (36.9%); Results 2.2 states the final model gave "accurate estimates of
    # population estimated parameters, and interindividual and inter-occasion
    # variability with a RSE% of less than 50%".
    #
    # The study is oral-only with no intravenous reference arm, so CL, V1, Q
    # and V2 are all APPARENT (/F) quantities; Hassouneh 2024 estimates no
    # bioavailability term, and F is implicitly 1 throughout.
    # ---------------------------------------------------------------------

    # Absorption. Hassouneh 2024 parameterizes the transit chain by BOTH the
    # transit rate constant Ktr and the mean transit time Mtt, which is the
    # Monolix `depot(target = ..., Ktr, Mtt)` macro. Monolix derives the
    # number of transit compartments from those two as N = Mtt*Ktr - 1
    # (equivalently Ktr = (N+1)/Mtt, the Savic 2007 relation), giving
    # N = 18.8*0.48 - 1 = 8.02 here. N is therefore a derived quantity and is
    # computed in model() rather than stored here. Hassouneh 2024 does not
    # print N; the same Monolix Ktr/Mtt pairing is spelled out explicitly by
    # the sibling extraction Jiang_2024_empagliflozin.R, whose paper states
    # the relation as Formula (1) and quotes the derived N.
    lktr <- log(18.8); label("Transit rate constant Ktr of the absorption chain (1/h)")                    # Hassouneh 2024 Table 3: Ktr = 18.8 (RSE 9.34%), described as "Transit time rate constant"
    lmtt <- log(0.48); label("Mean transit time Mtt for absorption (h)")                                   # Hassouneh 2024 Table 3: Mtt = 0.48 h (RSE 4.63%), described as "Meant transit time" [sic]
    lka  <- log(0.37); label("First-order absorption rate constant Ka into the central compartment at the reference BMI of 22.9 kg/m2 (1/h)") # Hassouneh 2024 Table 3: Ka = 0.37 (RSE 4.8%), and the covariate equation above Table 3, log(Ka) = log(0.37) - 0.85*log(BMI)

    # Disposition. Hassouneh 2024 uses the Monolix {Cl, V1, Q, V2}
    # parameterization; Table 3 labels V1 "Volume of the central compartment
    # (L)" and V2 "Volume of the peripheral compartment (L)" explicitly, so
    # V1 maps to the canonical lvc and V2 to lvp.
    lcl <- log(273.14); label("Apparent clearance CL/F (L/h)")                                             # Hassouneh 2024 Table 3: Cl = 273.14 (RSE 7.37%)
    lvc <- log(18.98);  label("Apparent central volume of distribution V1/F (L)")                          # Hassouneh 2024 Table 3: V1 = 18.98 L (RSE 8.94%)
    lq  <- log(64.62);  label("Apparent intercompartmental clearance Q/F (L/h)")                           # Hassouneh 2024 Table 3: Q = 64.62 L/h (RSE 3.92%)
    lvp <- log(487.9);  label("Apparent peripheral volume of distribution V2/F (L)")                       # Hassouneh 2024 Table 3: V2 = 487.9 L (RSE 3.18%)

    # Covariate effect. Power exponent of BMI on Ka; see covariateData$BMI for
    # the centering erratum (the published equation drops the reference BMI)
    # and for the numeric falsification of the uncentered reading.
    e_bmi_ka <- -0.85; label("Exponent of body mass index on the absorption rate constant (unitless)")     # Hassouneh 2024 Table 3: beta_BMI = -0.85 (RSE 36.9%), "Regression coefficient for the effect of body mass index (BMI) on Absorption rate constant"

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 3's column header is explicit:
    # "Between-Subject Variability (Standard Deviation (%RSE))", so the
    # published numbers are the STANDARD DEVIATIONS of the log-normal random
    # effects and the nlmixr2 omega entries below are those values SQUARED.
    # They are NOT %CV, so log(CV^2 + 1) must not be applied. Methods 4.3.1:
    # "A lognormal distribution was assumed for individual parameters", with
    # eta_i ~ N(0, omega) where "omega ... [is the] standard deviation of the
    # interindividual ... variability".
    # There is no IIV on V1, Q or V2 in the final model (those Table 3 cells
    # are blank).
    # ---------------------------------------------------------------------
    etalktr ~ 0.2401; # Hassouneh 2024 Table 3 BSV Ktr = 0.49 SD (RSE 26%)   -> variance 0.49^2 = 0.2401 (52.1% CV)
    etalmtt ~ 0.09;   # Hassouneh 2024 Table 3 BSV Mtt = 0.3 SD  (RSE 16.1%) -> variance 0.3^2  = 0.09   (30.7% CV)
    etalka  ~ 0.1296; # Hassouneh 2024 Table 3 BSV Ka  = 0.36 SD (RSE 12.7%) -> variance 0.36^2 = 0.1296 (37.1% CV)
    etalcl  ~ 0.3844; # Hassouneh 2024 Table 3 BSV Cl  = 0.62 SD (RSE 9.82%) -> variance 0.62^2 = 0.3844 (68.5% CV)

    # ---------------------------------------------------------------------
    # Interoccasion variability across the two SPRYCEL dosing occasions, on
    # Ktr, Mtt, Ka and CL. Table 3's column header is again explicit:
    # "Inter-Occasion Variability (Standard Deviation)", and Methods 4.3.1
    # gives eta_ki ~ N(0, gamma) with gamma the IOV standard deviation, so
    # the encoded omegas are the published SDs squared. Table 3 reports ONE
    # IOV SD per parameter, shared across occasions, so occasion 2 is
    # fixed() equal to occasion 1 (the NONMEM "$OMEGA BLOCK(1) SAME" idiom;
    # nlmixr2 has no SAME shortcut).
    # ---------------------------------------------------------------------
    etaiov_ktr_1 ~ 0.7056;          # Hassouneh 2024 Table 3 IOV Ktr = 0.84 SD (RSE 9.54%) -> variance 0.84^2 = 0.7056 (occasion 1; 100.5% CV)
    etaiov_ktr_2 ~ fixed(0.7056);   # IOV on log-Ktr, occasion 2; variance held equal to occasion 1 per the source's single-SD IOV reporting
    etaiov_mtt_1 ~ 0.16;            # Hassouneh 2024 Table 3 IOV Mtt = 0.4 SD  (RSE 8.04%) -> variance 0.4^2  = 0.16   (occasion 1; 41.6% CV)
    etaiov_mtt_2 ~ fixed(0.16);     # IOV on log-Mtt, occasion 2; variance held equal to occasion 1
    etaiov_ka_1  ~ 0.0961;          # Hassouneh 2024 Table 3 IOV Ka  = 0.31 SD (RSE 10.1%) -> variance 0.31^2 = 0.0961 (occasion 1; 31.7% CV)
    etaiov_ka_2  ~ fixed(0.0961);   # IOV on log-Ka, occasion 2; variance held equal to occasion 1
    etaiov_cl_1  ~ 0.1764;          # Hassouneh 2024 Table 3 IOV Cl  = 0.42 SD (RSE 7.94%) -> variance 0.42^2 = 0.1764 (occasion 1; 43.9% CV)
    etaiov_cl_2  ~ fixed(0.1764);   # IOV on log-CL, occasion 2; variance held equal to occasion 1

    # ---------------------------------------------------------------------
    # Residual variability. Results 2.1 prints the error model as
    #   Y = F + sqrt(a^2 + b^2 * F^2) * eps,   eps ~ N(0, 1)
    # which is Monolix's "combined2" model and is exactly nlmixr2's
    # add(addSd) + prop(propSd). The additive term carries the units of F;
    # every concentration axis in the paper (Figures 3-5) is ng/mL and the
    # LC-MS/MS calibration range is 0.50-500.00 ng/mL, so a = 0.78 is in
    # ng/mL. The model() block therefore reports Cc in ng/mL.
    # ---------------------------------------------------------------------
    addSd  <- 0.78; label("Additive residual error (ng/mL)")                                              # Hassouneh 2024 Table 3: a (constant) = 0.78 (RSE 5%)
    propSd <- 0.22; label("Proportional residual error (fraction)")                                       # Hassouneh 2024 Table 3: b (proportional) = 0.22 (RSE 1.72%)
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # 2. Per-occasion IOV etas. Each affected parameter picks up one of two
    # etas depending on the current occasion; the two variances are equal
    # (see ini()), so this reproduces a single shared IOV magnitude.
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_ka  <- oc1 * etaiov_ka_1  + oc2 * etaiov_ka_2
    iov_cl  <- oc1 * etaiov_cl_1  + oc2 * etaiov_cl_2

    # 3. Individual structural parameters. The BMI term is the power form of
    # Hassouneh 2024's log-transformed-BMI covariate model:
    # exp(beta * log(BMI/BMI_ref)) == (BMI/BMI_ref)^beta. The reference BMI of
    # 22.9 kg/m^2 is the cohort median (Table 2); the published equation omits
    # the centering entirely, which is falsified by the paper's own Figure 5
    # VPC -- see covariateData$BMI and the vignette Errata.
    ktr <- exp(lktr + etalktr + iov_ktr)
    mtt <- exp(lmtt + etalmtt + iov_mtt)
    ka  <- exp(lka  + etalka  + iov_ka) * (BMI / 22.9)^e_bmi_ka
    cl  <- exp(lcl  + etalcl  + iov_cl)
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)

    # Absorption-chain shape, derived from Ktr and Mtt exactly as Monolix
    # does: N = Mtt*Ktr - 1 (Savic 2007). At the typical values this is
    # 18.8*0.48 - 1 = 8.02. Because Ktr and Mtt each carry both IIV and IOV,
    # N varies by subject and occasion, which is the behaviour of the source
    # Monolix model. The gamma density below stays well defined for any
    # N > -1 (its shape parameter is N + 1 = Ktr*Mtt > 0), but for N < 0 the
    # input is singular at t -> 0 and the solver cannot integrate it; the
    # vignette counts and reports those occasions rather than dropping them
    # silently.
    nn <- ktr * mtt - 1

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. ODE system. The analytical Savic transit-compartment chain feeds the
    # depot, which then drains into central by first-order ka. The dose record
    # targets depot and the bolus is suppressed with f(depot) <- 0, so the
    # entire dose enters through the gamma-density input rate below.
    #
    # The gamma density is written out explicitly rather than via rxode2's
    # transit() macro. The two are algebraically identical -- the line below
    # is exactly what transit(nn, mtt, 1) expands to, and the vignette asserts
    # the two agree on Cmax / Tmax / AUC -- but the transit() macro silently
    # returns ZERO input for a model in nlmixr2 UI (ini()/model()) form when
    # combined with the conventional f(depot) <- 0 bolus suppression, which
    # would make the whole dose vanish without any error. podo(depot) is not
    # bioavailability-adjusted, so it still returns the full dose amount under
    # f(depot) <- 0. Hassouneh 2024 estimates no bioavailability term (CL and
    # the volumes are apparent /F quantities), so no F multiplies podo.
    d/dt(depot)       <- exp(log(podo(depot)) + log(ktr) + nn * log(ktr * tad(depot)) -
                               ktr * tad(depot) - lgamma(nn + 1)) - ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    # 6. Observation and error. Dose in mg and volumes in L give mg/L; the
    # source reports plasma dasatinib in ng/mL, which is 1000x that value
    # (1 mg/L = 1000 ng/mL). Reporting Cc directly in ng/mL keeps the
    # published additive residual term a = 0.78 ng/mL in its native units.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
