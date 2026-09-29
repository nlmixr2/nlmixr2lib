Daskapan_2019_darunavir <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for ritonavir-boosted darunavir in HIV-1-infected adult outpatients, developed by an iterative two-stage Bayesian procedure in MWPharm for therapeutic drug monitoring. Clearance is scaled by total body weight and volume by an MWPharm fat-corrected lean body mass; bioavailability and absorption rate are fixed at literature values and the renal clearance fraction is fixed at zero."
  reference <- paste(
    "Daskapan A, Tran QTD, Cattaneo D, Gervasoni C, Resnati C, Stienstra Y,",
    "Bierman WFW, Kosterink JGW, van der Werf TS, Proost JH, Alffenaar JWC,",
    "Touw DJ. Darunavir Population Pharmacokinetic Model Based on HIV",
    "Outpatient Data. Ther Drug Monit. 2019;41(1):59-65.",
    "doi:10.1097/FTD.0000000000000576. PMCID: PMC6358182.",
    sep = " "
  )
  vignette <- "Daskapan_2019_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters twice. (1) Linearly on the metabolic clearance arm as",
        "cl_nonren = CLm * (WT / 70), the MWPharm per-70-kg clearance",
        "convention (Methods, CL equation); the reference weight 70 kg is the",
        "unit denominator of the published CLm, not a cohort median (the",
        "hospital A median was 72.0 kg). (2) Through the fat-corrected lean",
        "body mass LBMc = LBM + (WT - LBM) * fd that scales the volume.",
        "Recorded at the outpatient visit of the drug-level measurement",
        "(Methods, Data Collection). Hospital A median 72.0 kg, range 40-123 kg",
        "(Table 1). Missing values were imputed at national-statistics averages",
        "(80 kg male / 70 kg female CBS, 75 kg male / 65 kg female ISTAT), but",
        "no weight was missing in the hospital A development set.",
        "DOMAIN BOUND: see the LBMc note under HT -- the fd = 5 setting makes",
        "LBMc, and hence the volume, negative for sufficiently lean subjects."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Input to the sex-specific lean-body-mass formula",
        "LBM = 50.0 + 0.9 * (HT - 152) for men and",
        "LBM = 45.5 + 0.9 * (HT - 152) for women (Methods, V equation), which",
        "is then fat-corrected as LBMc = LBM + (WT - LBM) * fd with fd = 5 and",
        "multiplies the per-kg volume Vd.",
        "DOMAIN BOUND: with fd = 5 the correction EXTRAPOLATES rather than",
        "interpolates (fd = 1 would give LBMc = WT, fd = 0 would give",
        "LBMc = LBM), so LBMc = 5 * WT - 4 * LBM is positive only when",
        "WT > 0.8 * LBM. For a man that is WT > 40 + 0.72 * (HT - 152) kg and",
        "for a woman WT > 36.4 + 0.72 * (HT - 152) kg; a 173 cm man needs",
        "WT > 55.1 kg and a 193 cm man WT > 69.5 kg. Below the bound the",
        "published model returns a negative volume. The bound is reachable",
        "inside the source cohort's own reported ranges (Table 1 gives",
        "height up to 193 cm, weight down to 40 kg and BMI down to",
        "16.9 kg/m2), so a virtual cohort must be restricted to",
        "physiologically paired weight-height combinations and checked for",
        "LBMc > 0. This is a property of the published parameterisation, not",
        "of the transcription.",
        "Hospital A median 173.0 cm, range 150-193 cm (Table 1). Missing values",
        "were imputed at national-statistics averages (180 cm male / 170 cm",
        "female CBS, 175 cm male / 165 cm female ISTAT), but no height was",
        "missing in the hospital A development set."
      ),
      source_name = "Height"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the intercept of the lean-body-mass formula only",
        "(50.0 for men, 45.5 for women; Methods, V equation). Sex has no",
        "separate effect on clearance or on the volume coefficient. Hospital A",
        "was 29% female (57 of 198; Table 1)."
      ),
      source_name = "sex"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the CKD-EPI equation.",
        "Multiplied by the fixed renal fraction e_crcl_cl_renal, which the",
        "source set to zero, so this column is INERT in the final model."
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The published clearance equation is",
        "CL = CLm * (WT / 70) + fr * CLcr, with CLcr from the CKD-EPI formula",
        "converted to L/h (Methods, CL equation). The model file converts the",
        "canonical mL/min/1.73 m^2 column to L/h inside model() by the factor",
        "0.06 and carries fr as e_crcl_cl_renal.",
        "fr was FIXED AT ZERO in the final model: Supplement 4 grid-searched",
        "fr = 0.12 (the darunavir SPC renal fraction) against fr = 0 and the",
        "latter gave the lower AIC (1584.89 vs 1611.08), which the Discussion",
        "attributes to darunavir being ~80% hepatically eliminated. The arm is",
        "retained in the model file so the tested-and-rejected renal pathway",
        "stays visible and a user can restore fr = 0.12 by overriding",
        "e_crcl_cl_renal; at the published value of 0 the CRCL column has no",
        "effect on any prediction.",
        "Serum creatinine (the CKD-EPI input) was recorded at the drug-level",
        "visit or within 15 days of it; hospital A median 83.5 umol/L,",
        "range 44.2-230.7 umol/L (Table 1). CKD-EPI additionally requires age",
        "and sex, so a user computing this column needs both; the source",
        "does not state which CKD-EPI race coefficient (if any) was applied."
      ),
      source_name = "CLcr"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 198,
    n_studies = 1,
    n_observations = 198,
    age_range = "24-74 years",
    age_median = "54 years",
    weight_range = "40-123 kg",
    weight_median = "72.0 kg",
    height_range = "150-193 cm",
    height_median = "173.0 cm",
    bmi_range = "16.9-35.3 kg/m^2",
    bmi_median = "24.6 kg/m^2",
    sex_female_pct = 29,
    race_ethnicity = "Not reported.",
    disease_state = "HIV-1 infection, adult outpatients on ritonavir-boosted darunavir",
    renal_function = "Serum creatinine median 83.5 umol/L, range 44.2-230.7 umol/L; renal function not used as a covariate in the final model (fr = 0)",
    dose_range = "darunavir/ritonavir 800/100 mg once daily (162 of 198, 82%) or 600/100 mg twice daily (36 of 198, 18%)",
    regions = "Italy (development set, ASST Fatebenefratelli Sacco University Hospital, Milano)",
    notes = paste(
      "Development set = 'hospital A' (ASST Fatebenefratelli Sacco, Milano),",
      "a retrospective therapeutic-drug-monitoring record review over",
      "April 2015 - August 2017: 198 unique adult patients contributing",
      "198 darunavir plasma samples, i.e. ONE sample per patient (Results,",
      "Data set; Table 1). Samples with unknown intake or sampling time, or",
      "below the 0.2 mg/L lower limit of quantification, were excluded, and",
      "few samples fell in the 0-4 h absorption phase, which is why the",
      "absorption rate constant could not be estimated and was fixed.",
      "A separate external validation set of 170 patients / 170 samples from",
      "University Medical Center Groningen ('hospital B', January 2010 - May",
      "2017; Netherlands) was used for Passing-Bablok and Bland-Altman",
      "agreement analysis and did NOT contribute to parameter estimation.",
      "Baseline demographics for both hospitals are in Table 1."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Table 2, 'Model 2 AIC = 1584.89', marked by footnote * as the chosen
    # final population pharmacokinetic model. Model 1 (AIC 945.31, only CLm
    # Bayesian) is a rejected development submodel -- its volume was held at
    # the literature value so it 'implies that the volume of distribution
    # is the same for each patient, which does not seem logical' (Results)
    # -- and is therefore not packaged (base-vs-final policy).
    #
    # Reported as 'Mean (95% CI)' with a separate interindividual
    # 'SD (95% CI)' column, both on the NATURAL scale of the parameter
    # (Supplement 1 and Supplement 3 print the same 'value (SD)' pairs in
    # the same units). Parameters were assumed log-normally distributed
    # (Methods, Population Pharmacokinetic Model Development), so the
    # reported mean is carried as the typical value and the reported SD is
    # converted to a log-scale variance by omega^2 = log(1 + CV^2) with
    # CV = SD / mean. See the vignette Errata for the alternative reading
    # in which the printed mean is the arithmetic mean of the log-normal.
    # ---------------------------------------------------------------------

    # Absorption. Not estimable from these data: 'the low number of
    # darunavir samples drawn in the absorption phase; 0-4 hours after drug
    # intake' (Discussion), so it was held at the Arab-Alameddine 2014
    # literature value carried in as the ITSB prior.
    lka <- fixed(log(1.04))
    label("Absorption rate constant (1/h; literature value, Arab-Alameddine 2014)") # Table 2 Ka = 1.04, footnotes + and ++ 'literature value' / 'set on fixed value'; Supplement 1 one-compartment prior Ka 1.04 (0.35)

    # Metabolic (non-renal) clearance arm, per 70 kg total body weight.
    # Named cl_nonren because the published equation adds a renal arm
    # fr * CLcr on top of it; total cl is their sum inside model().
    lcl_nonren <- log(9.47)
    label("Metabolic clearance CLm at 70 kg total body weight (L/h/70 kg)") # Table 2, Model 2, CLm = 9.47 (95% CI 8.24-10.65) L/h/70kgBW

    # Volume coefficient, per kg of fat-corrected lean body mass. Kept as
    # the paper-named canonical 'vd' rather than 'lvc' because the estimate
    # is a per-kg coefficient (L/kgLBMc), not an absolute central volume;
    # the absolute volume vc = vd * LBMc is derived inside model().
    lvd <- log(2.13)
    label("Volume of distribution per kg fat-corrected lean body mass (L/kg LBMc)") # Table 2, Model 2, Vd = 2.13 (95% CI 1.39-3.26) L/kgLBMc

    # Bioavailability. 'Because of the absence of data on drug
    # concentrations after parenteral darunavir administration as a
    # comparison for oral administration to measure bioavailability,
    # bioavailability was fixed in all parameterizations at the literature
    # value of 0.82' (Results, Population Pharmacokinetic Model).
    lfdepot <- fixed(log(0.82))
    label("Oral bioavailability (fraction; literature value, darunavir SPC)") # Table 2 F = 0.82, footnote section-sign 'literature value from SPC'

    # fd, the MWPharm 'fat distribution' factor weighting body weight in
    # excess of lean body mass inside LBMc = LBM + (WT - LBM) * fd. This is
    # NOT an allometric exponent -- it is a linear interpolation /
    # extrapolation coefficient between LBM (fd = 0) and WT (fd = 1), and
    # the selected value of 5 extrapolates well beyond WT. Grid-searched,
    # not estimated: Supplement 4 tried 1, 2, 5, 10 and 20 and 5 gave the
    # best AIC that also stayed inside the paper's own 2-point selection
    # threshold. See the vignette Errata for the negative-LBMc domain bound
    # this value creates.
    e_wt_vd <- fixed(5)
    label("Fat-distribution factor fd weighting (WT - LBM) inside LBMc (unitless interpolation coefficient, not an exponent)") # Table 2 'Fat distribution' = 5; Supplement 4 AIC grid 1:1 3743.29, 1:2 3676.84, 1:5 3603.29, 1:10 3630.27, 1:20 3602.20; Discussion 'a ratio of fat distribution (fd) of 5 ... provided better AIC scores'

    # fr, the ratio of darunavir renal clearance to creatinine clearance,
    # i.e. the slope of the renal clearance arm on CLcr. Fixed at zero in
    # the final model, so the CRCL covariate is inert; see
    # covariateData$CRCL$notes.
    e_crcl_cl_renal <- fixed(0)
    label("Ratio fr of darunavir renal clearance to creatinine clearance (unitless slope of the renal CL arm on CLcr in L/h)") # Table 2, Model 2, fr = 0; Supplement 4 renal-fraction grid fr = 0.12 AIC 1611.08 vs fr = 0 AIC 1584.89

    # IIV. Log-scale variances from the Table 2 natural-scale
    # interindividual SD column via omega^2 = log(1 + (SD / mean)^2):
    #   CLm  SD 6.19 / mean 9.47 -> CV 65.4 % -> log(1 + 0.653643^2) = 0.355749
    #   Vd   SD 2.60 / mean 2.13 -> CV 122.1 % -> log(1 + 1.220657^2) = 0.912284
    # No covariance between them is reported, so the block is diagonal.
    etalcl_nonren ~ 0.355749 # Table 2, Model 2, CLm 'SD (95% CI)' = 6.19 (4.85-7.76) L/h/70kgBW -> 65.4 % CV
    etalvd ~ 0.912284 # Table 2, Model 2, Vd 'SD (95% CI)' = 2.60 (1.43-4.66) L/kgLBMc -> 122.1 % CV

    # Residual error. 'the residual error was assumed to be normally
    # distributed and equal to the SD of the assay, which was estimated as
    # 0.2 + 0.05 * C, where C is the observed darunavir plasma
    # concentration' (Methods, Population Pharmacokinetic Model
    # Development). The LINEAR sum of an additive and a proportional term
    # is nlmixr2's combined1() form, not the default quadrature
    # combination. Both terms are assay characteristics rather than fitted
    # quantities, so both are fixed. The additive term equals the reported
    # 0.2 mg/L lower limit of quantification.
    addSd <- fixed(0.2)
    label("Additive residual SD (mg/L; assay SD intercept, not estimated)") # Methods 'the SD of the assay, which was estimated as 0.2 + 0.05 * C'
    propSd <- fixed(0.05)
    label("Proportional residual SD (fraction; assay SD slope, not estimated)") # Methods 'the SD of the assay, which was estimated as 0.2 + 0.05 * C'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariate terms
    # ------------------------------------------------------------------
    # Sex-specific lean body mass (kg) from body height in cm.
    # Methods, V equation: 'LBM is calculated from 50.0 + 0.9 * (Height -
    # 152) for male patients and 45.5 + 0.9 * (height - 152) for female
    # patients'.
    lbm <- (50.0 - 4.5 * SEXF) + 0.9 * (HT - 152)

    # Fat-corrected lean body mass (kg), the size descriptor for the
    # volume. Methods, V equation: 'LBMc = LBM + (BW - LBM) * fd'.
    lbmc <- lbm + (WT - lbm) * e_wt_vd

    # Creatinine clearance on the paper's scale. The canonical CRCL column
    # is mL/min/1.73 m^2; the published equation wants L/h, so convert by
    # 60 min/h / 1000 mL/L = 0.06 (Methods, 'converted to unit L/h').
    clcr <- 0.06 * CRCL

    # ------------------------------------------------------------------
    # 2. Individual parameters
    # ------------------------------------------------------------------
    ka <- exp(lka)

    # Methods, CL equation: CL = CLm * (BW / 70) + fr * CLcr. The two arms
    # are summed into `cl` so that the total -- not one arm -- is the
    # quantity driving elimination.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (WT / 70)
    cl_renal <- e_crcl_cl_renal * clcr
    cl <- cl_nonren + cl_renal

    # Methods, V equation: V = V1 * LBMc, with V1 (Table 2 'Vd') in
    # L/kgLBMc.
    vd <- exp(lvd + etalvd)
    vc <- vd * lbmc

    # ------------------------------------------------------------------
    # 3. Micro-constants
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. ODE system -- one compartment with first-order absorption and
    #    first-order elimination (Results: 'A 1-compartment model with a
    #    first-order absorption and elimination ... resulted in the best
    #    model'). A second compartment was explored and rejected: its
    #    peripheral volume estimated to 0.051 L/kg, 'which is negligible
    #    as a significant peripheral compartment'.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ------------------------------------------------------------------
    # 5. Bioavailability
    # ------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 6. Observation and error
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
