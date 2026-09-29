Song_2017_vancomycin <- function() {
  description <- "Two-compartment intravenous population pharmacokinetic model for vancomycin in 316 Chinese neonates and young infants (postnatal age under 60 days at admission) treated at the Children's Hospital of Chongqing Medical University (Song 2017). Clearance carries two power covariates referenced to the cohort medians -- birth body weight (exponent 0.888, reference 3.22 kg) and postnatal age (exponent 0.449, reference 29 days) -- and is the only parameter retaining interindividual variability; the shrinkage on V1, V2 and Q exceeded 0.5 so their etas were dropped during model building. Residual error is purely additive (2.187 ug/mL). The authors report median clearance 0.106 L/h/kg and median volume of distribution 0.935 L/kg, both larger than previously published Caucasian neonatal values. The disposition is written with explicit k12/k21 micro-constant ODEs; the validation vignette solves with useLinCmt = FALSE defensively, which on rxode2 5.1.8 reproduces the default solve to eight significant figures."
  reference <- "Song L, He CY, Yin NG, Liu F, Jia YT, Liu Y. A population pharmacokinetic model for individualised dosage regimens of vancomycin in Chinese neonates and young infants. Oncotarget. 2017;8(62):105211-105221. doi:10.18632/oncotarget.22114"
  vignette <- "Song_2017_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given by intravenous infusion (Materials
  # and Methods, 'Patients and data collection': 'suspected or confirmed
  # bacterial infection that necessitated intravenously infusion of
  # vancomycin'), so the dose enters `central` directly and there is no depot
  # state. The specimen is verified as SERUM: the same section states 'Serum
  # vancomycin concentrations were measured by our in-house Clinical
  # Pharmacokinetic Service, using a fluorescence polarization immunoassay
  # method', with a lower limit of quantification of 1 ug/mL.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT_BIRTH = list(
      description = "Birth body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at birth. Enters clearance as the power term (WT_BIRTH / 3.22)^0.888, where 3.22 kg is the cohort MEDIAN birth body weight stated in the Results 'PPK model' paragraph ('where 3.22 (kg) is the median BBW') and confirmed by Table 1 (birth body weight 3.22 kg, range 1.25-5.38). Song 2017 calls the column BBW. Note that the retained size descriptor is birth weight, NOT the current body weight (Table 1 median 3.95 kg, range 1.25-7.62), which was screened separately and not retained -- the Discussion argues birth body weight combined with postnatal age was preferable 'as weight is partly correlated to PNA'. Supplying current weight in this column would misstate every clearance. The Discussion notes only one of twelve previous neonatal vancomycin studies identified birth body weight as a significant covariate on clearance.",
      source_name = "BBW"
    ),
    PNA = list(
      description = "Postnatal age at the time the vancomycin concentration was determined",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Song 2017 reports postnatal age in DAYS (Table 1 row 'Postnatal age at vancomycin determined, d': median 29, range 2-77), but the canonical nlmixr2lib PNA column carries MONTHS, so the paper's clearance term (PNA_days / 29)^0.449 is reparameterised inside model() as (PNA_months / 0.952772)^0.449 using 1 month = 30.4375 days, giving the reference 29 / 30.4375 = 0.952772 months. Numerator and denominator carry the same units factor so it cancels exactly and the exponent is unchanged. This follows the same reparameterisation precedent as Zhao_2018_omeprazole.R (days) and Bardhi_2026_ampicillin_foal.R (hours). Users must supply PNA in MONTHS. Note this is postnatal age at sampling, which the paper distinguishes from postnatal age at admission (Table 1 median 24 days, range 0-60) used for the study inclusion criterion of under 60 days.",
      source_name = "PNA"
    )
  )

  # Screened in the Song 2017 stepwise forward-inclusion / backward-elimination
  # covariate search but NOT retained in the final model, so they are
  # documentation only and are not referenced in model(). The full candidate
  # list is given in Materials and Methods, 'Model development': 'the candidate
  # covariates were sex, GA, PNA (at vancomycin concentrations determined),
  # PMA, BBW, BW, HT, BSA, SCr, GFR, BUN, ALT, ALB and concomitant drug
  # therapy'. Only BBW and PNA survived (Table 2 and the Results 'PPK model'
  # final equation).
  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      notes = "The ONLY covariate that entered the full model and was then removed by backward-elimination. Table 2 records it being added at forward-inclusion step 3 (model 3, OFV 2502.56, delta-OFV -7.26 versus model 2) and then removed from model 4 with a delta-OFV of only 4.71, which failed the pre-specified backward-elimination threshold of an OFV increase exceeding 10.828 (p < 0.001, chi-squared with 1 df); the Table 2 P column marks this row '> 0.001'. The final model is therefore model 6 = basic + PNA + BBW. Song 2017 calls the column PMA. The paper does not tabulate postmenstrual age in Table 1, so no cohort summary is available for it.",
      source_name = "PMA"
    ),
    WT = list(
      description = "Current body weight at the time of sampling",
      units = "kg",
      type = "continuous",
      notes = "Screened but not retained; birth body weight was the retained size descriptor (see covariateData$WT_BIRTH$notes). Cohort median 3.95 kg, range 1.25-7.62 (Table 1). This value is nonetheless load-bearing for reading the paper's summary statistics: the reported median clearance of 0.106 L/h/kg is 0.42 L/h divided by 3.95 kg, and the reported median volume of distribution of 0.935 L/kg is (1.27 + 2.422) L divided by 3.95 kg. Song 2017 calls the column BW.",
      source_name = "BW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 37 weeks, range 28-41 (Table 1). 102 of the 316 patients (32%) were premature.",
      source_name = "GA"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "",
      type = "binary",
      notes = "Screened but not retained. Cohort 115 female of 316 (36.4%); Table 1 reports counts as male/female pairs (201/115).",
      source_name = "sex"
    ),
    HT = list(
      description = "Body height / length",
      units = "cm",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 49 cm, range 35-62.3 (Table 1). Also an input to the Schwartz glomerular filtration rate formula used to derive the screened GFR covariate.",
      source_name = "HT"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened but not retained. Derived by Song 2017 rather than measured, using the formula quoted in Materials and Methods: BSA (m^2) = body weight(kg)^0.5378 * height(cm)^0.3964 * 0.024265. Not summarised in Table 1.",
      source_name = "BSA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "umol/L as reported by Song 2017 Table 1 (median 28.6, range 12-151); divide by 88.4 to obtain mg/dL. Screened but not retained, and the Discussion explains the exclusion on physiological grounds rather than statistical ones: serum creatinine 'is known to be influenced by age, sex, muscle mass and diet - limiting its utility as a marker of the glomerular filtration rate, especially for neonates - it was excluded from our modelling process; this practice is common in the related studies'. The Discussion also notes that 4 of 12 previous neonatal vancomycin studies did retain serum creatinine on clearance, so this is a genuine point of divergence between published models rather than a universal finding.",
      source_name = "SCr"
    ),
    CRCL = list(
      description = "Creatinine-based estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened but not retained. Derived by Song 2017 using the Schwartz formula quoted in Materials and Methods: GFR (mL/min/1.73 m^2) = k * HT(cm) / SCr(mg/dL), with k = 0.45 for full-term infants and 0.33 for preterm infants. Note the formula consumes serum creatinine in mg/dL while Table 1 reports it in umol/L. Not summarised in Table 1.",
      source_name = "GFR"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "mmol/L as reported by Song 2017 Table 1 (median 2.35, range 0.48-7.54), NOT the US-convention mg/dL; multiply by 2.80 to obtain mg/dL. Screened but not retained.",
      source_name = "BUN"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 25.7 U/L, range 5.4-130.9 (Table 1).",
      source_name = "ALT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "g/L, matching the register-canonical SI unit, as reported by Song 2017 Table 1 (median 33.1, range 13.5-52.6). Screened but not retained.",
      source_name = "ALB"
    ),
    CONMED_PANIPENEM = list(
      description = "Concomitant panipenem therapy indicator",
      units = "",
      type = "binary",
      notes = "Screened but not retained. One of the two drugs by which Song 2017 operationalised the 'concomitant drug therapy' candidate covariate (Results 'PPK model': 'the concomitant drug therapy covariate represented by panipenem and furosemide'). Entered as a categorical effect of the form P_i = P_pop * exp(dPdCov * Cov) with Cov a 0/1 dummy (Materials and Methods, 'Model development'). 168 of 316 patients (53.2%) received panipenem (Table 1). The Discussion confirms the outcome: 'Our modelling did, however, consider concomitant therapies (i.e. panipenem and furosemide) as categorical covariates, but the factors produced no significant impact on the final model.' Concomitant drug use overall involved only about 20% of the studied patients.",
      source_name = "panipenem treatment"
    ),
    CONMED_FUROSEMIDE = list(
      description = "Concomitant furosemide therapy indicator",
      units = "",
      type = "binary",
      notes = "Screened but not retained. The second of the two drugs representing the 'concomitant drug therapy' candidate covariate; see the CONMED_PANIPENEM note for the shared encoding and the Discussion's statement that neither had a significant impact. 74 of 316 patients (23.4%) received furosemide (Table 1).",
      source_name = "furosemide treatment"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 316L,
    n_studies = 1L,
    age_range = "Postnatal age at admission 0-60 days (median 24); postnatal age at vancomycin sampling 2-77 days (median 29); gestational age at birth 28-41 weeks (median 37). Study inclusion required age under 60 days at admission.",
    age_median = "Postnatal age 29 days at vancomycin determination (Table 1)",
    weight_range = "Current body weight 1.25-7.62 kg (median 3.95); birth body weight 1.25-5.38 kg (median 3.22)",
    weight_median = "3.95 kg current body weight; 3.22 kg birth body weight (Table 1)",
    sex_female_pct = 36.4,
    race_ethnicity = "Chinese (single-centre cohort in Chongqing, China). The paper's central comparative claim is that clearance and volume of distribution in this population exceed previously published Caucasian and Malaysian neonatal values.",
    disease_state = "Neonates and young infants receiving intravenous vancomycin for suspected or documented Gram-positive bacterial infection. Pathogenic culture was positive in 137 of 316 cases (59 from blood, the remainder from sputum, fester and puncture fluid); 179 were culture-negative or, for a few admitted before 2012, had no culture. 102 of 316 (32%) were premature.",
    dose_range = "13.7-73.5 mg/kg/day (median 36.7), divided into 2 to 4 intravenous infusions per day",
    regions = "China (Children's Hospital of Chongqing Medical University, Chongqing; external evaluation at the same centre and a preliminary two-patient validation at Southwest Hospital, Third Military Medical University)",
    renal_function = "Serum creatinine 12-151 umol/L (median 28.6) and blood urea nitrogen 0.48-7.54 mmol/L (median 2.35) (Table 1). Patients receiving renal replacement therapy were excluded, so the model's domain of applicability is the observed range; no renal-impairment stratum is defined.",
    notes = "Retrospective single-centre analysis of routine clinical therapeutic-drug-monitoring data collected between November 2011 and December 2016, comprising 421 serum vancomycin concentrations from 316 patients (1-6 samples per patient; median observed concentration 9.33 ug/mL, range 1.34-38.65). Estimated in Phoenix NLME version 1.3 (Certara), with initial estimates from Phoenix WinNonlin 6.4. Exclusion criteria were renal replacement therapy, vancomycin treatment for under 24 h, and missing demographic data. Assay was fluorescence polarization immunoassay with a lower limit of quantification of 1 ug/mL. A one-compartment model was rejected in favour of two compartments on objective function value (2764.31 versus 2623.54). Internal evaluation used diagnostic scatter plots, a 2000-replicate nonparametric bootstrap (Table 3) and a visual predictive check from 1000 simulations (Figure 4). External evaluation used 27 samples from 19 further patients admitted January-June 2017, giving mean prediction error -0.29 +/- 0.99, mean absolute error 1.388 +/- 0.71 and mean squared prediction error 1.928 +/- 1.665 (units printed as ng/mL in the Results, which is inconsistent with the ug/mL concentration scale used throughout; see the vignette Errata). 233 of 347 observed trough concentrations were below the 10 ug/mL target."
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters -- Song 2017 Table 3 (final model column) and the
    # Results 'PPK model' final equation, which prints the same four values:
    #
    #   V1 (l) = 1.27
    #   V2 (l) = 2.422
    #   CL (l/h) = 0.42 * (BBW/3.22)^0.888 * (PNA/29)^0.449 * exp(eta_CL)
    #   Q  (l/h) = 1.161
    #
    # The reference subject is therefore a neonate at the cohort median birth
    # body weight of 3.22 kg and median postnatal age of 29 days, at which
    # both covariate terms equal exactly 1 and clearance equals 0.42 L/h.
    # Volumes and intercompartmental clearance carry no covariates.
    #
    # Consistency check on the transcription, using the cohort median CURRENT
    # weight of 3.95 kg: 0.42 / 3.95 = 0.1063 L/h/kg and (1.27 + 2.422) / 3.95
    # = 0.9347 L/kg, reproducing the Results sentence 'about 0.106 l/h/kg ...
    # and 0.935 l/kg'.
    # ------------------------------------------------------------------------
    lvc <- log(1.27); label("Central volume of distribution V1 (L)")  # Song 2017 Table 3 (tvV1 = 1.27 L, SE 0.191, 95% CI 0.895-1.644; bootstrap median 1.255)
    lvp <- log(2.422); label("Peripheral volume of distribution V2 (L)")  # Song 2017 Table 3 (tvV2 = 2.422 L, SE 0.425, 95% CI 1.586-3.258; bootstrap median 2.386)
    lcl <- log(0.42); label("Clearance at 3.22 kg birth body weight and 29 days postnatal age (L/h)")  # Song 2017 Table 3 (tvCL = 0.42 L/h, SE 0.0124, 95% CI 0.395-0.444; bootstrap median 0.416)
    lq <- log(1.161); label("Intercompartmental clearance Q (L/h)")  # Song 2017 Table 3 (tvQ = 1.161 L/h, SE 0.177, 95% CI 0.814-1.509; bootstrap median 1.159)

    # ------------------------------------------------------------------------
    # Covariate power exponents on clearance. Both were ESTIMATED, not held
    # constant: Table 3 reports a standard error, a 95% confidence interval and
    # a bootstrap distribution for each. The Table 3 legend calls them 'fixed
    # parameter coefficient of birth body weight / postnatal age', where
    # 'fixed' is the fixed-effect (THETA) sense standard in Phoenix NLME, not
    # the held-constant sense -- so neither is wrapped in fixed() here.
    #
    # Continuous covariates were implemented as P_i = P_pop *
    # (Cov / Cov_median)^dPdCov (Materials and Methods, 'Model development').
    # ------------------------------------------------------------------------
    e_wtbirth_cl <- 0.888; label("Power exponent on (WT_BIRTH / 3.22 kg) for clearance (unitless)")  # Song 2017 Table 3 (dCldBBW = 0.888, SE 0.12, 95% CI 0.652-1.124; bootstrap median 0.885)
    e_pna_cl <- 0.449; label("Power exponent on (PNA / 29 days) for clearance (unitless)")  # Song 2017 Table 3 (dCldPNA = 0.449, SE 0.058, 95% CI 0.336-0.563; bootstrap median 0.453)

    # ------------------------------------------------------------------------
    # Interindividual variability, exponential: P_i = P_pop * exp(eta_i), with
    # eta normally distributed with mean 0 and variance omega^2 (Materials and
    # Methods, 'Model development').
    #
    # Clearance is the ONLY parameter retaining an eta. Results 'PPK model':
    # the shrinkage factors of V1, V2 and Q were all above 0.5, 'indicating
    # minor inter-individual variability of the parameters that could be
    # eliminated without significantly altering the parameter values or OFV;
    # as such, each were excluded in the model building process'.
    #
    # 0.317 is a VARIANCE, not a standard deviation. The Table 3 legend
    # distinguishes the two scales explicitly within a single sentence --
    # 'stdev0 = standard deviation; omega-CL = variance of the
    # inter-individual variability of CL' -- so the authors were labelling the
    # residual row as an SD and this row as a variance deliberately. A
    # variance of 0.317 corresponds to an SD of 0.563 on the log scale and an
    # apparent coefficient of variation of sqrt(exp(0.317) - 1) = 61.1%.
    # ------------------------------------------------------------------------
    etalcl ~ 0.317  # Song 2017 Table 3, row 'omega-CL' = 0.317 variance (SE 0.015, 95% CI 0.288-0.346; bootstrap median 0.316)

    # ------------------------------------------------------------------------
    # Residual error, purely additive. Results 'PPK model': 'the additive model
    # best described the residual variability', chosen over the proportional
    # and combined alternatives that were also evaluated (Materials and
    # Methods, 'Model development'). Table 3's legend confirms stdev0 is a
    # standard deviation, so it is carried directly as addSd in the paper's
    # ug/mL concentration units.
    # ------------------------------------------------------------------------
    addSd <- 2.187; label("Additive residual error (ug/mL)")  # Song 2017 Table 3 (stdev0 = 2.187 ug/mL, SE 0.194, 95% CI 1.807-2.568; bootstrap median 2.158)
  })
  model({
    # ----------------------------------------------------------------------
    # 1. Derived covariate terms -- Song 2017 Results 'PPK model':
    #
    #      CL (l/h) = 0.42 * (BBW/3.22)^0.888 * (PNA/29)^0.449 * exp(eta_CL)
    #
    #    Reference constants are the cohort medians: birth body weight
    #    3.22 kg and postnatal age 29 DAYS. The canonical PNA column is in
    #    MONTHS, so the age reference is rescaled once here:
    #    29 days / 30.4375 days per month = 0.952772 months. Both numerator
    #    and denominator of the ratio carry the same units factor, so it
    #    cancels exactly and the exponent 0.449 is unchanged.
    # ----------------------------------------------------------------------
    wtbirth_ref <- 3.22
    pna_ref_months <- 29 / 30.4375

    f_wtbirth <- (WT_BIRTH / wtbirth_ref)^e_wtbirth_cl
    f_pna <- (PNA / pna_ref_months)^e_pna_cl

    # 2. Individual parameters. Only clearance carries an eta and only
    #    clearance carries covariates.
    cl <- exp(lcl + etalcl) * f_wtbirth * f_pna
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 4. ODE system -- two compartments with first-order elimination from the
    #    central compartment and intravenous dosing only (Materials and
    #    Methods, 'Model development'; vancomycin was given as an intravenous
    #    infusion, so there is no absorption step).
    #
    #    The disposition uses k12/k21 micro-constants, a shape for which
    #    rxSolve()'s default useLinCmt = TRUE auto-converts the ODEs to
    #    linCmt(). That conversion has been observed elsewhere in this library
    #    to drop the peripheral state silently; it does NOT do so here -- on
    #    rxode2 5.1.8, with the etas zeroed, useLinCmt = TRUE and FALSE agree
    #    on Cc to eight significant figures and both retain peripheral1. The
    #    vignette nonetheless passes useLinCmt = FALSE defensively and pairs
    #    its exposure gate with an independent terminal-half-life gate, which
    #    is the check that would go red if a future version did collapse the
    #    model (the collapse preserves Dose/CL, so an AUC gate cannot see it).
    # ----------------------------------------------------------------------
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error. Dose in mg, vc in L -> mg/L = ug/mL, which is
    #    the unit Song 2017 uses for every reported concentration.
    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
