Tseng_2026_vancomycin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for vancomycin in adult",
    "neurosurgical patients with intracranial hemorrhage managed with an",
    "external ventricular drain (Tseng 2026). Fit in Phoenix NLME by FOCE-ELS",
    "to plasma concentrations from nine prospectively sampled patients in",
    "northern Taiwan. Twelve candidate covariates were screened by stepwise",
    "search and none reached significance on clearance or volume, so the model",
    "carries no covariate effects. Interindividual variability is retained on",
    "clearance and intercompartmental clearance only; the central and",
    "peripheral volume random effects were dropped from the final model for",
    "shrinkage above 0.9. Residual error is multiplicative (proportional).",
    "The paper's headline cerebrospinal-fluid penetration results are",
    "noncompartmental and linear-regression analyses that sit OUTSIDE this",
    "population PK model, which is plasma-only and has no CSF compartment."
  )
  reference <- paste(
    "Tseng YJ, Juan L, Wu CC, Lan YX, Chen GY, Huang APH, Chen KW, Wang KC,",
    "Luh HT, Lin SW. Population pharmacokinetics and cerebrospinal fluid",
    "penetration of intravenous vancomycin in intracranial hemorrhage patients",
    "with external ventricular drains: implications for dosing and therapeutic",
    "drug monitoring. Drug Des Devel Ther. 2026;20:1-14.",
    "doi:10.2147/DDDT.S574548. PMCID: PMC12965103."
  )
  vignette <- "Tseng_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482. Vancomycin was given intravenously, so the dose enters
  # `central` directly and there is no depot state. Tseng 2026 Sampling
  # Protocol: "3 mL of blood ... were collected ... Immediately after
  # collection, all samples were centrifuged"; "Plasma vancomycin
  # concentrations were quantified using a validated LC-MS/MS method with a
  # linear calibration range of 0.78-100 ug/mL". The assayed matrix is
  # therefore explicitly plasma, so `central` is verified. `peripheral1` is a
  # mathematical distribution compartment that was never sampled, so it is
  # recorded as unverified.
  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  # No covariate is used by this model. Tseng 2026 Results, PopPK Modeling:
  # "patient-specific factors including age, sex, disease state, CLcr, eGFR,
  # culture results, CSF WBC, CSF glucose level, CSF total protein,
  # cell-index, drainage output, and urine output levels were incorporated
  # into the model. However, none of these variables demonstrated a
  # statistically significant impact on Vd or CL (p <= 0.05)."
  covariateData <- list()

  # The twelve screened-but-not-retained covariates, for those concepts that
  # already have a canonical column in inst/references/covariate-columns.md.
  # These are documentation only -- none is referenced in model(). The
  # remaining screened concepts (CSF white blood cell count, CSF glucose, the
  # CSF cell index, external-ventricular-drain output volume, CSF culture
  # positivity, and the intracerebral-versus-subarachnoid hemorrhage
  # indication) have no canonical register entry; because they are
  # documentation only and no new canonical is being ratified here, they are
  # recorded narratively in population$notes rather than invented as columns.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 2 median 67 years (IQR 18). Inclusion required age 20 years or older. Screened on CL and Vd in the stepwise covariate search (Results, PopPK Modeling) and not retained.",
      source_name        = "age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Table 2 reports 8 of 9 patients (88.9%) male, so 11.1% female. Screened and not retained. With one female patient in the cohort a sex effect was not estimable in practice.",
      source_name        = "sex"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 2 median 82.4 mL/min (IQR 23.9). The Table 1 and Table 2 footnotes define it as (140 - age) x BW / (72 x SCr), with a serum creatinine at or below 0.8 mg/dL floored to 0.8 mg/dL. Screened and not retained. Tseng 2026 Discussion attributes the absence of a renal-function effect to the NTUH dosing nomogram (Table 1) having already adjusted doses for creatinine clearance, so the covariate's influence was largely designed out of the observed data, and to the narrow range of laboratory values in a nine-patient cohort. Note that the paper screened both this Cockcroft-Gault CLcr and a separate eGFR; it does not state which eGFR equation was used, and neither was retained.",
      source_name        = "CLcr"
    ),
    CSF_TPRO = list(
      description        = "Cerebrospinal-fluid total protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed among the screened covariates in Results, PopPK Modeling as 'CSF total protein'; not retained. Tseng 2026 does not tabulate its distribution. The paper's ventriculitis case definition uses a CSF protein threshold of 50 mg/dL (= 0.5 g/L in the canonical SI units of this column).",
      source_name        = "CSF total protein"
    ),
    URINE_VOL_24H = list(
      description        = "Daily urine output",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Recorded per Methods, Data Collection ('Physiological parameters, such as daily urine output and EVD drainage were also documented') and screened on CL and Vd; not retained in the popPK model. It IS the strongest single predictor in the paper's separate linear-regression analysis of CSF exposure, where AUC_CSF = -5.084 x (urine output) + 29.089 with R^2 = 0.510 and p = 0.0136, and it also enters the multivariate regression. That regression is a statistical model of an NCA-derived AUC, not a pharmacokinetic covariate relationship, and it is not encoded here; see the vignette. Tseng 2026 does not state the units in which urine output entered the regression, so the regression coefficient is not directly interpretable.",
      source_name        = "urine output"
    ),
    WBC = list(
      description        = "Blood white blood cell count",
      units              = "10^3/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 2 median 8.29 K/uL (IQR 1.31). Reported as a baseline characteristic; the separately screened covariate in the popPK search is the CSF white blood cell count, not this blood count. Not retained.",
      source_name        = "WBC"
    ),
    CRP = list(
      description        = "C-reactive protein",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 2 median 2.31 mg/dL (IQR 2.69), measured as high-sensitivity CRP per Methods, Data Collection. Reported as a baseline characteristic and not among the covariates listed as entered into the popPK search; not in the final model.",
      source_name        = "CRP"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 9L,
    n_studies        = 1L,
    n_sites          = 2L,
    age_median       = "67 years (IQR 18)",
    age_range        = "adults 20 years or older",
    weight_median    = "71.2 kg (IQR 11.2)",
    height_median    = "163 cm (IQR 15.375)",
    sex_female_pct   = 11.1,
    race_ethnicity   = "Not reported. The cohort was enrolled at two medical centers in northern Taiwan and is described in the Conclusion as a Taiwanese neurosurgical population.",
    disease_state    = "Adults with intracranial hemorrhage requiring placement of an external ventricular drain and receiving intravenous vancomycin. EVD indication was intracerebral hemorrhage in 6 of 9 (66.7%) and subarachnoid hemorrhage in 3 of 9 (33.3%); 1 of 9 (11.1%) also had a lumbar drain. Meningitis/ventriculitis was microbiologically confirmed in 2 of 9 (22.2%, Enterococcus faecium and Serratia marcescens) and suspected in the remaining 7 (77.8%). Exclusions: history of metastatic or intracranial tumor, prior cranial radiation, prior craniotomy, or ongoing immunosuppressive therapy.",
    renal_function   = "Cockcroft-Gault creatinine clearance median 82.4 mL/min (IQR 23.9).",
    dose_range       = "Total daily doses of 1000 to 3000 mg intravenously, given per the National Taiwan University Hospital AUC-based nomogram (Table 1), which selects 250-1250 mg at q8h, q12h, once-daily or every-other-day intervals from creatinine clearance and body weight. All patients received a loading dose of 20-25 mg/kg (25 mg/kg for critically ill or end-stage-renal-disease patients), capped at 3000 mg per dose. 8 of 9 patients (88.9%) adhered to the protocol; one was underdosed relative to it. The infusion duration is NOT reported anywhere in the paper.",
    regions          = "Taiwan (two medical centers in the north of the country)",
    n_concentrations = "Not reported as a total. The design specifies 5 plasma and 5 CSF samples per patient during the fifth dosing interval (pre-dose, end of infusion, and 4, 6 and 8 hours post-dose), so the plasma dataset the popPK model was fit to is on the order of 45 observations.",
    notes            = "Prospective observational study, enrollment January 2023 to December 2024, NTUH REC approval 202302126RINA. Sampling was performed after the fourth dose to approximate steady state. Plasma vancomycin was assayed by validated LC-MS/MS over 0.78-100 ug/mL (R^2 > 0.99) with intra- and inter-day precision and accuracy within +/-15% RSD. The popPK model was fit in Phoenix NLME by FOCE-ELS; one-, two- and three-compartment infusion models and additive, multiplicative and combined residual-error models were compared (Supplementary Table S2, not on disk), and a two-compartment model with multiplicative residual error won. Final-model fit statistics: LogLik = -137.185, -2LL = 274.370, AIC = 288.370, BIC = 301.017, nParm = 7, EPS shrinkage = 0.130, condition number = 23.800. Evaluation was by bootstrap (Table 5) and visual predictive check plus goodness-of-fit plots (Supplementary Figures S1-S2, not on disk). TWELVE COVARIATES WERE SCREENED and none was significant on CL or Vd: age, sex, disease state, CLcr, eGFR, culture results, CSF WBC, CSF glucose level, CSF total protein, cell index, EVD drainage output, and urine output. Six of those twelve concepts (CSF WBC, CSF glucose, cell index, EVD drainage output, culture results, and the hemorrhage-type disease state) have no canonical covariate column in this library and, being documentation-only here, are recorded in this note rather than as invented columns. THE PAPER'S CSF RESULTS ARE NOT PART OF THIS MODEL: CSF penetration was characterized by noncompartmental analysis (AUC_CSF/plasma 0.84-14.22%; median plasma AUC0-24 478.10 mg*hr/L versus median CSF AUC 10.29 mg*hr/L) and by Spearman correlation and linear regression, with the end-of-infusion CSF/plasma concentration ratio the best single-time-point surrogate (Spearman r = 0.791, p = 0.004). No CSF compartment was fit, so none is encoded. Tseng 2026 Study Limitations state that the nine-patient sample size limits power, raises overfitting risk, and makes the popPK and regression analyses exploratory and hypothesis-generating rather than confirmatory."
  )

  ini({
    # STRUCTURAL PARAMETERS -- Tseng 2026 Table 5, "Final Model" block. Table 4
    # reports the same two-compartment model BEFORE the covariate search and
    # before the zero-shrinkage random effects were dropped (tvV 6.085,
    # tvV2 45.117, tvCL 4.288, tvCL2 20.157); Table 5 supersedes it and is the
    # final model, so Table 5 is used throughout.
    #
    # Phoenix NLME names the two-compartment clearance parameterization
    # V / V2 / CL / CL2, where V is the central volume, V2 the peripheral
    # volume, CL the elimination clearance and CL2 the intercompartmental
    # clearance. These map onto the library canonicals vc / vp / cl / q.
    #
    # INDEPENDENT CROSS-CHECK of that mapping: tvV + tvV2 = 6.065 + 45.117 =
    # 51.18 L, against the Table 3 median individual volume of distribution of
    # 50.46 L from the separate WinNonlin noncompartmental/individual analysis
    # of the same nine patients -- so the two volumes are central and
    # peripheral and their sum is Vss, not two alternative whole-body volumes.
    # Likewise tvCL 4.289 L/hr against the Table 3 median individual CL of
    # 4.63 L/hr.
    lvc <- log(6.065);  label("Central volume of distribution (L)")               # Tseng 2026 Table 5, Final Model, tvV = 6.065 L (Stderr 1.409, RSE 23.234%, 95% CI 3.212-8.918); bootstrap median 6.080, 95% CI 0.009-19.579
    lvp <- log(45.117); label("Peripheral volume of distribution (L)")            # Tseng 2026 Table 5, Final Model, tvV2 = 45.117 L (Stderr 8.232, RSE 18.222%, 95% CI 28.512-61.843); bootstrap median 46.929, 95% CI 32.462-265.621
    lcl <- log(4.289);  label("Clearance (L/h)")                                  # Tseng 2026 Table 5, Final Model, tvCL = 4.289 L/hr (Stderr 0.656, RSE 15.297%, 95% CI 2.960-5.617); bootstrap median 4.239, 95% CI 0.362-6.011
    lq  <- log(20.208); label("Intercompartmental clearance (L/h)")               # Tseng 2026 Table 5, Final Model, tvCL2 = 20.208 L/hr (Stderr 7.274, RSE 35.995%, 95% CI 5.483-34.932); bootstrap median 18.031, 95% CI 8.673-46.945

    # NO COVARIATE EFFECTS. Twelve covariates were screened and none reached
    # the forward-addition criterion of a 3.84-point OFV drop (Results, PopPK
    # Modeling). Accordingly there is no e_<cov>_<param> term in this model.

    # INTERINDIVIDUAL VARIABILITY -- Tseng 2026 Table 5, Final Model, the two
    # 'omega' rows. Only CL and CL2 carry a random effect: Results, PopPK
    # Modeling states "Since both omega-V and omega-V2 exhibited shrinkage
    # values greater than 0.9, these random effects were excluded from the
    # final model" (Table 4 gives omega-V shrinkage 0.937 and omega-V2
    # shrinkage 0.980). So there is deliberately no etalvc and no etalvp here.
    #
    # SCALE: these are VARIANCES of a log-normally distributed random effect,
    # not standard deviations and not CV fractions, so they are used as
    # written on the nlmixr2 eta variance scale. Three independent lines of
    # evidence, none of which relies on the Phoenix default alone:
    #  (1) Phoenix NLME reports the Omega matrix itself in the Random Effect
    #      block, whose diagonal entries are variances, and it names the
    #      residual term separately and explicitly as 'stdev()' -- the output
    #      distinguishes the two scales by name.
    #  (2) The printed relative standard errors match a VARIANCE. For n = 9
    #      subjects the expected RSE of a variance is about
    #      sqrt(2/(n-1)) = 50%, versus about 1/sqrt(2*(n-1)) = 25% for a
    #      standard deviation. Table 5 prints 62.649% for omega-CL and
    #      49.436% for omega-CL2 -- both consistent with a variance and about
    #      twice what an SD would give.
    #  (3) The magnitude matches the observed spread. The nine individual
    #      clearances in Table 3 (4.63, 6.17, 2.87, 4.65, 3.46, 13.81, 4.18,
    #      2.91, 4.83 L/hr) have a log-scale standard deviation of 0.477, i.e.
    #      a variance of 0.227 -- close to omega-CL = 0.188 read as a
    #      variance, and irreconcilable with 0.188 read as an SD (which would
    #      imply a variance of 0.035, roughly six-fold too small).
    # Read as variances, omega-CL = 0.188 is an apparent CV of about 45% and
    # omega-CL2 = 0.796 an apparent CV of about 103%. See the vignette Errata.
    etalcl ~ 0.188  # Tseng 2026 Table 5, Final Model, 'omega CL' = 0.188 variance (Stderr 0.118, RSE 62.649%, shrinkage 0.020); bootstrap median 0.171, 95% CI 0.012-2.043
    etalq  ~ 0.796  # Tseng 2026 Table 5, Final Model, 'omega CL2' = 0.796 variance (Stderr 0.394, RSE 49.436%, shrinkage 0.183); bootstrap median 0.624, 95% CI 0.000-2.759

    # RESIDUAL ERROR. Tseng 2026 selected "a multiplicative residual error
    # structure" (Results, PopPK Modeling; also the title of Tables 4 and 5).
    # Phoenix NLME's multiplicative error is Cobs = Cpred * (1 + eps) with
    # eps ~ N(0, stdev^2), which is exactly nlmixr2's prop(propSd). The row is
    # labelled 'stdev()' in the source table, so this one IS on the standard
    # deviation scale -- unlike the two omega rows above.
    propSd <- 0.164; label("Proportional residual error (fraction)")  # Tseng 2026 Table 5, Final Model, 'stdev()' = 0.164 (Stderr 0.014, RSE 8.594%, 95% CI 0.136-0.193); bootstrap median 0.159, 95% CI 0.118-0.198; EPS shrinkage 0.130
  })
  model({
    # Individual PK parameters. No covariates enter the model (none was
    # significant), and only CL and Q carry interindividual variability.
    cl <- exp(lcl + etalcl)
    q  <- exp(lq + etalq)
    vc <- exp(lvc)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L, so central/vc is mg/L == ug/mL, the units
    # Tseng 2026 reports vancomycin concentrations in (LC-MS/MS calibration
    # range 0.78-100 ug/mL) and consistent with the reported AUC units of
    # mg*hr/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
