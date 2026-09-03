Ma_2026_colistinSulfate <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous colistin sulfate in",
    "critically ill adults with carbapenem-resistant organism infections",
    "(Ma 2026; n = 51 Chinese ICU patients; 123 sparse therapeutic-drug-",
    "monitoring serum samples spanning 0.12-4.40 mg/L). Linear elimination",
    "from a single central compartment with a 1-hour intravenous infusion",
    "input. Cockcroft-Gault creatinine clearance enters clearance as a power",
    "function centred on 94.76 mL/min (exponent 0.232); sex was screened but",
    "eliminated in backward elimination. Inter-individual variability was",
    "estimated on CL only; residual error is proportional. Colistin sulfate",
    "is administered as the active drug and must not be confused with",
    "colistimethate sodium (CMS), the inactive prodrug modelled in",
    "Plachouras 2009, Mohamed 2012, Jacobs 2016 and Karaiskos 2015."
  )
  reference <- paste(
    "Ma Y, Wang Y, Wu X, Wang J, Pang Y, Jia Y, Yang X, Gu J (2026).",
    "Optimizing colistin sulfate dosing in severe infections: a population",
    "pharmacokinetic model-guided approach.",
    "Drug Des Devel Ther 20:1-16.",
    "doi:10.2147/DDDT.S600942.",
    sep = " "
  )
  vignette <- "Ma_2026_colistinSulfate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Methods "Dosage Administration and Blood Sample Collection": 200 uL of
    # SERUM was extracted and "the serum mass concentration was analyzed"
    # by HPLC-MS/MS (linear 0.1-10 ug/mL, LLOQ 0.1 ug/mL). The assayed
    # analyte is colistin sulfate itself, not a colistimethate metabolite.
    central = list(analyte = "colistin sulfate", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation,",
        "reported as raw mL/min and NOT normalised to 1.73 m^2 body surface",
        "area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on CL:",
        "CL = 1.66 * (CrCl / 94.76)^0.232 * exp(eta) per the final-model",
        "equation printed on Ma 2026 p. 7, with the exponent 0.232 also",
        "tabulated as the 'CL_CCR' row of Table 3 (RSE 11.1%; bootstrap",
        "median 0.23, 95% CI 0.10-0.35). The centring constant 94.76 mL/min",
        "appears ONLY inside that printed equation -- it is not tabulated,",
        "and it is not the cohort mean (Table 1 gives CrCL 100.42 +/- 68.87",
        "mL/min); it is most consistent with the cohort median, which the",
        "paper does not report. The Cockcroft-Gault form itself is not",
        "written out; Clinical Data Collection states only that CrCL was",
        "'calculated using the Cockcroft-Gault formula'. Fitted over an",
        "observed range of 12.14-246.91 mL/min (Discussion); the paper's own",
        "Monte Carlo simulations extrapolate marginally below that, to 10",
        "mL/min, and flag the extrapolation explicitly. Patients on CRRT or",
        "ECMO were excluded, so the model carries no information about renal",
        "replacement therapy."
      ),
      source_name        = "CrCL"
    )
  )

  # Screened during stepwise covariate modelling but not retained in the final
  # model. Sex was carried into the full regression model and then removed in
  # backward elimination ("Gender was eliminated during the backward
  # elimination process of the full regression model due to its lack of
  # significant impact"); the paper reports no coefficient for it. The
  # remaining laboratory covariates were screened by SCM but the paper reports
  # neither their forward-inclusion dOFV nor any coefficient.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "unitless",
      type               = "categorical",
      reference_category = "male",
      notes              = paste(
        "Screened as an exponential categorical effect in the stepwise",
        "covariate model and entered the full regression model, but removed",
        "during backward elimination (Results, 'The Population",
        "Pharmacokinetic Model of Colistin Sulfate'). No point estimate is",
        "reported, so nothing is encoded."
      ),
      source_name        = "Gender"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    n_observations = 123L,
    age_range      = "18 years or older (inclusion criterion); interquartile range 42.00-70.00 years",
    age_median     = "65.00 years",
    weight_range   = "not reported; BMI interquartile range 21.43-26.10 kg/m^2",
    weight_median  = "not reported; BMI median 24.05 kg/m^2",
    sex_female_pct = 100 * 19 / 51,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill intensive-care patients with confirmed",
      "carbapenem-resistant organism (CRO) infections -- predominantly",
      "pulmonary, abdominal, central-nervous-system, bloodstream and urinary",
      "tract. Isolates were Acinetobacter baumannii (58.00%), Klebsiella",
      "pneumoniae (28.00%), Pseudomonas aeruginosa (6.00%), Escherichia coli",
      "(4.00%) and Enterobacter cloacae (4.00%). Mean APACHE II score",
      "21.22 +/- 7.80. All patients received one or two concomitant",
      "antibacterials, most often tigecycline (58.82%), a carbapenem",
      "(17.65%) or cefoperazone-sulbactam (13.73%). Patients on continuous",
      "renal replacement therapy or ECMO, pregnant or lactating patients,",
      "and patients who died within 24 h of the first dose were excluded."
    ),
    dose_range     = paste(
      "Intravenous colistin sulfate at the label-recommended 1 to 1.5",
      "million IU per day divided into 2 or 3 doses, with the exact regimen",
      "chosen by the attending physician. Table 1 gives a median daily dose",
      "of 100.00 (100.00, 150.00) x 10^4 IU; the Discussion restates this as",
      "a median 0.68 mg/kg/day. Treatment lasted 3 to 31 days (median 11.50",
      "days). The paper's Monte Carlo simulations standardise the infusion",
      "duration at 1 hour for every regimen."
    ),
    sampling       = paste(
      "Sparse therapeutic drug monitoring after steady state was reached",
      "(five doses), at up to three time points per patient: peak 0.5 h",
      "after the end of the infusion, an intermediate sample 6 h after the",
      "end of the infusion, and trough within 0.5 h before the next dose.",
      "Of the 123 samples, 41.46% were troughs, 34.96% peaks and 23.58%",
      "6-hour samples."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 100.42 +/- 68.87 mL/min (Table",
      "1), observed range 12.14-246.91 mL/min (Discussion). Serum creatinine",
      "median 71.00 (46.00, 118.00) umol/L. Patients on CRRT were excluded."
    ),
    regions        = "People's Republic of China (single centre; Second Hospital of Hebei Medical University, Shijiazhuang, Hebei).",
    notes          = paste(
      "Baseline demographics from Ma 2026 Table 1 and Results 'Basic",
      "Information of the Patient'. Retrospective single-centre cohort",
      "collected June 2021 to June 2023 (ethics approval 2020-R551). The",
      "final model was fit in NONMEM 7.2 with FOCE-I via Pirana 2.8.0 and",
      "PsN 5.16.3, and evaluated with a 1000-sample bootstrap and a 1000-",
      "replicate prediction- and variability-corrected VPC. Serum, not",
      "plasma, was the assayed matrix. The unbound fraction of 0.5 used for",
      "the paper's fAUC/MIC target-attainment simulations is a literature",
      "assumption, not a fitted parameter, and is therefore not encoded",
      "here; the paper notes as a limitation that free concentrations were",
      "calculated rather than measured."
    )
  )

  ini({
    # Structural parameters. Typical values refer to the covariate reference
    # subject, CrCl 94.76 mL/min (Ma 2026 final-model equation, p. 7).

    lcl <- log(1.66);  label("Clearance (CL, L/h) at CrCl 94.76 mL/min")  # Table 3, TVCL 1.660 (RSE 7.7%; bootstrap median 1.65, 95% CI 1.44-1.92); also the leading coefficient of the p. 7 equation
    lvc <- log(10.10); label("Central volume of distribution (V, L)")     # Table 3, TVV 10.10 (RSE 7.1%; bootstrap median 10.05, 95% CI 8.32-12.87); also the p. 7 equation "V (L) = 10.10"

    # Covariate effect: creatinine clearance on CL as a power function of the
    # covariate divided by 94.76 mL/min (p. 7 equation
    # CL (L/h) = 1.66 * (CrCl/94.76)^0.232 * e^ETA1).
    e_crcl_cl <- 0.232; label("Power exponent of CrCl on CL (unitless)")  # Table 3, "CL_CCR" 0.232 (RSE 11.1%; bootstrap median 0.23, 95% CI 0.10-0.35)

    # Inter-individual variability, exponential on CL only (Methods: "Inter-
    # individual variability was assessed employing the exponential model";
    # the p. 7 equation carries e^ETA1 on CL and none on V).
    #
    # SCALE. Table 3 heads this row "Between-subject variability omega^2 CL"
    # and reports 0.282, which read literally is a VARIANCE (omega = 0.531,
    # 57.4% CV). Two independent published outputs refute that reading and
    # show 0.282 is the standard deviation (28.7% CV); it is encoded here as
    # 0.282^2 = 0.079524.
    #   1. Table 4. Css,avg = daily dose / (CL * 24), so within a CrCL stratum
    #      the z-scores of P(Css,avg > 2 mg/L) across dose levels differ by
    #      log(dose ratio) / omega -- a dose-unit-free determination. Every
    #      pair in the CCR = 10 and CCR = 50 columns returns omega = 0.27-0.30
    #      (e.g. CCR = 50: 5.60% and 82.50% for 1 and 2 MU/day give
    #      log(2) / (qnorm(0.825) - qnorm(0.056)) = 0.275).
    #   2. Figure 3, the prediction- and variability-corrected VPC. At the
    #      ~11.5 h trough the observed 5th and 95th percentiles are ~0.22 and
    #      ~1.92 mg/L, a spread of 0.658 on the log scale. At steady state
    #      dlog(C)/dlog(CL) = -1.57 there, so the predicted spread is
    #      sqrt((1.57*omega)^2 + propSd^2): 0.591 for omega = 0.282 versus
    #      0.922 for omega = 0.531.
    # The vignette runs both checks ("iiv-scale-check") and ships the rejected
    # omega = 0.531 counterfactual alongside them.
    etalcl ~ 0.079524  # omega = 0.282 (Table 3; RSE 24.5%, shrinkage 10.4%; bootstrap median 0.272, 95% CI 0.196-0.361)

    # Residual error: proportional only. Table 3 heads the row "Residual
    # variability sigma^2", but the Results text defines the model as
    # DV = IPRED * (1 + theta_prop * eps) with "the variance of EPS (1) ...
    # fixed to define the scale of the random error term, whereas the
    # proportional residual error coefficient theta_prop was estimated".
    # Under that parameterisation sigma^2 is fixed at 1 and is not an
    # estimated quantity, so the tabulated 0.391 can only be theta_prop --
    # i.e. a proportional SD of 39.1%, which is what is encoded.
    propSd <- 0.391; label("Proportional residual error (fraction)")  # Table 3, theta_prop 0.391 (RSE 7.7%, shrinkage 12.4%; bootstrap median 0.391, 95% CI 0.321-0.450)
  })

  model({
    # Individual PK parameters. The covariate enters as a power function
    # centred on 94.76 mL/min (Ma 2026 p. 7 final-model equation). V carries
    # no covariate and no inter-individual variability.
    cl <- exp(lcl + etalcl) * (CRCL / 94.76)^e_crcl_cl
    vc <- exp(lvc)

    # Colistin sulfate is given intravenously; dose records target `central`
    # directly (a 1-hour infusion in the paper's simulations), so there is no
    # depot compartment.
    d/dt(central) <- -(cl / vc) * central

    # Concentration in mg/L (= ug/mL, the assay's reporting unit).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
