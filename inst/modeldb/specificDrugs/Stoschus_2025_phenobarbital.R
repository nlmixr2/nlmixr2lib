Stoschus_2025_phenobarbital <- function() {
  description <- paste(
    "One-compartment population PK model for oral and intravenous",
    "phenobarbital in critically ill adults with refractory and",
    "superrefractory status epilepticus (Stoschus 2025), with first-order",
    "absorption (ka fixed at 1.9 /h from Viswanathan 1979 because the",
    "therapeutic-drug-monitoring samples were taken after the absorption",
    "phase), estimated oral bioavailability F = 0.96, and ideal body weight",
    "entering as allometric scaling on both volume of distribution",
    "(exponent 1) and clearance (exponent 0.75), centered on the cohort",
    "median IBW of 68.8 kg. Interindividual variability is carried on V and",
    "CL, interoccasion variability on CL across consecutive TDM sampling",
    "intervals, with a proportional residual-error model."
  )
  reference <- paste(
    "Stoschus M, Schmidbauer ML, Starp J, Kunst S, Gakis G, Paal M,",
    "Vogeser M, Scharf-Janssen C, Liebchen U, Dimitriadis K (2025).",
    "Optimizing phenobarbital dosing in critically ill patients with",
    "refractory and superrefractory status epilepticus using a population",
    "pharmacokinetic model.",
    "Epilepsia 66(11):3757-3768.",
    "doi:10.1111/epi.18517.",
    sep = " "
  )
  vignette <- "Stoschus_2025_phenobarbital"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "phenobarbital", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "phenobarbital", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    IBW = list(
      description        = "Ideal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. The only covariate retained in the final",
        "model (Stoschus 2025 Results 3.2): it enters as allometric scaling",
        "on V (exponent 1, Equation 1) and on CL (exponent 0.75, Equation",
        "2), each centered on the cohort median IBW of 68.8 kg (Table 1",
        "row 'Ideal body weight, kg = 68.8 (53.3, 81.9)'). Adding IBW",
        "lowered the residual unexplained variability by 10.8% and the OFV",
        "by 71.9 versus the base model. Stoschus 2025 Methods 2.1 states",
        "IBW was computed with the formula of Brower et al. (paper",
        "reference 25); the formula itself is not reproduced in the",
        "article, so a user reproducing the covariate must consult that",
        "reference rather than assume a Devine-family variant. Simulations",
        "and the published dosing table were restricted to IBW 55-80 kg,",
        "within the observed cohort range of 52.4-82.4 kg; the paper warns",
        "that extrapolation beyond this range could mislead."
      ),
      source_name        = "IBW"
    ),
    OCC = list(
      description        = "Integer-valued sampling-occasion indicator for interoccasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Stoschus 2025 Methods 2.2 defines an occasion as the interval",
        "between consecutive TDM samples, with a cohort median interval of",
        "48 h. The paper reports a single IOV magnitude on CL (Table 2,",
        "36.85% CV) and does not state a fixed number of occasions, so this",
        "file encodes eight occasions of 48 h each, spanning 0-384 h. That",
        "covers the 336 h at which the paper assessed steady-state Cmin",
        "for the maintenance-dose simulations (336 h falls in occasion 8).",
        "The occasion index is decomposed inside model() into",
        "mutually-exclusive binary indicators oc1..oc8 multiplexing the",
        "eight per-occasion CL etas; occasions 2-8 are fixed equal to",
        "occasion 1 because a single shared IOV variance was estimated.",
        "Records with OCC outside 1..8 carry no IOV. For a single-occasion",
        "simulation pass OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Stoschus 2025 Methods 2.2 'patient characteristics') but not retained; Results 3.2 states that apart from IBW 'other covariates did not significantly reduce the OFV'."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a patient characteristic (Stoschus 2025 Methods 2.2) but not retained in the final model."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened alongside height and IBW (Stoschus 2025 Methods 2.2). Ideal body weight was the size descriptor retained; the paper's Discussion argues IBW 'might better reflect physiological factors like organ function and distribution, which do not always correlate with changes in total body weight'."
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a patient characteristic (Stoschus 2025 Methods 2.2) but not retained; it enters the model only indirectly through the derivation of IBW."
    ),
    RRT = list(
      description = "Renal replacement therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a clinical parameter (Stoschus 2025 Methods 2.2); 3 of 37 patients (8%) received renal replacement therapy (Table 1). Not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as laboratory data (Stoschus 2025 Methods 2.2); cohort median 0.7 mg/dL (Table 1). Not retained in the final model."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as laboratory data (Stoschus 2025 Methods 2.2) but not retained in the final model; the paper does not tabulate its cohort distribution."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened as laboratory data (Stoschus 2025 Methods 2.2) but not retained in the final model; the paper does not tabulate its cohort distribution."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as laboratory data (Stoschus 2025 Methods 2.2); cohort median 27 U/L (Table 1). Not retained in the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as laboratory data (Stoschus 2025 Methods 2.2); cohort median 28 U/L (Table 1). Not retained in the final model."
    ),
    CONMED_VALPROATE = list(
      description = "Concomitant valproate indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as time-varying comedication recorded per day of phenobarbital treatment (Stoschus 2025 Methods 2.2) but not retained. The Discussion notes that earlier studies reported reduced CL with concomitant valproic acid and that, given the common coadministration of other antiepileptics in this cohort, such an effect 'cannot be ruled out and might be included in the estimated CL parameter'."
    ),
    CONMED_PHENYTOIN = list(
      description = "Concomitant phenytoin indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as time-varying comedication (Stoschus 2025 Methods 2.2) but not retained; see the CONMED_VALPROATE note on the possibility that such effects are absorbed into the estimated CL."
    ),
    CONMED_METAMIZOLE = list(
      description = "Concomitant metamizole indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as time-varying comedication (Stoschus 2025 Methods 2.2) but not retained in the final model."
    ),
    CONMED_CENOBAMATE = list(
      description = "Concomitant cenobamate indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as time-varying comedication (Stoschus 2025 Methods 2.2) but not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 37L,
    n_studies      = 1L,
    age_range      = "24.0-83.0 years (5th-95th percentile); all patients >= 18 years by inclusion criterion",
    age_median     = "63.0 years",
    weight_range   = "51.7-93.3 kg total body weight (5th-95th percentile)",
    weight_median  = "78.5 kg total body weight; 68.8 kg ideal body weight (the allometric-scaling reference)",
    sex_female_pct = 35.1,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Refractory (2 patients, 5%) and superrefractory (35 patients, 95%)",
      "status epilepticus in a neurointensive care unit. Median Status",
      "Epilepticus Severity Score 5. Patients with a history of",
      "hypoxic-ischemic encephalopathy were excluded. Renal replacement",
      "therapy in 3 patients (8%); median estimated glomerular filtration",
      "rate 102 mL/min by CKD-EPI, median serum creatinine 0.7 mg/dL."
    ),
    dose_range     = paste(
      "Guideline-directed phenobarbital: loading doses 10-20 mg/kg and",
      "maintenance doses 1-4 mg/kg, median administered dose 200 mg (5th-95th",
      "percentile 100-600 mg) across 1242 doses. 824 doses (66%) were oral",
      "(tablets dispersed in 5 mL water and given via gastric tube) and 418",
      "(34%) intravenous as 5-min infusions. Median treatment duration 16",
      "days."
    ),
    regions        = "Germany (single center: LMU University Hospital, Munich).",
    notes          = paste(
      "Retrospective single-center cohort of patients treated 2015-2024,",
      "contributing 301 therapeutic-drug-monitoring samples (median 6 per",
      "patient, 5th-95th percentile 1-26) with a median phenobarbital",
      "concentration of 32.8 mg/L. Trough sampling was routine on Mondays,",
      "Wednesdays, and Fridays; serum concentrations were measured by",
      "HPLC-UV over a quantification range of 0.6-150 mg/L. Height was",
      "missing for 4 patients and weight for 2; these were imputed with the",
      "population median. The model was fit in MONOLIX 2024R1 (SAEM).",
      "Developed for adults only: Stoschus 2025 states the model 'is not",
      "intended for pediatric use' and that extrapolation beyond the adult",
      "population is not appropriate."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters (Stoschus 2025 Table 2, 'Fixed-effects
    # parameters'). V and CL are the typical values at the cohort median
    # ideal body weight of 68.8 kg, because the allometric terms in
    # Equations 1 and 2 are centered on medIBW.
    # =========================================================================
    lfdepot <- log(0.96)
    label("Oral bioavailability (unitless)")
    # Table 2 row 'F = .96 (3.95; .85-.99)'. Estimable because the cohort
    # received both oral (66%) and intravenous (34%) doses; Discussion
    # paragraph 2 notes this was previously unreported in adult ICU patients.

    lka <- fixed(log(1.9))
    label("First-order absorption rate constant from depot to central (1/h), literature value from Viswanathan 1979")
    # Table 2 row 'k a, h-1 = 1.9' with footnote 'a k a was fixed according
    # to literature values.' Methods 2.2: 'As the majority of the TDMs were
    # taken after the absorption phase of oral medication, estimation of the
    # absorption rate constant (k a) was not possible. Therefore, k a was
    # fixed at 1.9 h-1 according to the literature.' (paper reference 29,
    # Viswanathan et al., cited in the Discussion as having 'specifically
    # studied the bioavailability and absorption rate of oral phenobarbital,
    # capturing samples within the absorption phase').

    lvc <- log(34.3)
    label("Volume of distribution at the reference ideal body weight of 68.8 kg (L)")
    # Table 2 row 'V, L = 34.3 (12.7; 26.81-43.89)'

    lcl <- log(0.38)
    label("Total body clearance at the reference ideal body weight of 68.8 kg (L/h)")
    # Table 2 row 'CL, L/h = .38 (7.81; .32-.44)'

    # =========================================================================
    # Allometric exponents on ideal body weight. Both are printed as literal
    # exponents in the covariate equations rather than as estimated values
    # with uncertainty, so both are encoded as fixed.
    # =========================================================================
    e_ibw_vc <- fixed(1)
    label("Allometric exponent on V with ideal body weight (unitless)")
    # Results 3.2 Equation 1: V_i = V * (IBW_i / medIBW)^1

    e_ibw_cl <- fixed(0.75)
    label("Allometric exponent on CL with ideal body weight (unitless)")
    # Results 3.2 Equation 2: CL_i = CL * (IBW_i / medIBW)^.75

    # =========================================================================
    # Interindividual variability (Stoschus 2025 Table 2, 'IIV parameters').
    # A lognormal distribution was assumed for IIV and IOV (Methods 2.2), and
    # the CV% column is derived from the SD column via the standard lognormal
    # relation (Appendix S1: Equation 1). nlmixr2lib stores variances on the
    # log scale, so omega^2 = log(1 + CV^2). The CV% is printed to four
    # significant figures whereas the SD column is rounded to two, so the
    # variances below are back-solved from the CV; each reproduces the printed
    # SD to its stated precision (sqrt(0.50799) = 0.713 vs printed .71,
    # sqrt(0.157914) = 0.397 vs printed .4, sqrt(0.127331) = 0.357 vs
    # printed .36).
    # =========================================================================
    etalvc ~ 0.50799
    # Table 2 IIV row 'V, L: SD .71 (12.1; .56-.9), CV% 81.36';
    # omega^2 = log(1 + 0.8136^2) = 0.50799

    etalcl ~ 0.157914
    # Table 2 IIV row 'CL, L/h: SD .4 (16.2; .29-.54), CV% 41.36';
    # omega^2 = log(1 + 0.4136^2) = 0.157914

    # =========================================================================
    # Interoccasion variability on CL (Stoschus 2025 Table 2, 'IOV
    # parameters'). One IOV magnitude is reported and occasions are defined
    # as the intervals between consecutive TDM samples (median 48 h, Methods
    # 2.2). Eight 48-h occasions are encoded, spanning 0-384 h so that the
    # 336-h steady-state assessment time of the maintenance-dose simulations
    # is covered. Occasions 2-8 share occasion 1's variance because a single
    # value was estimated.
    # =========================================================================
    etaiov_cl_1 ~ 0.127331
    # Table 2 IOV row 'CL, L/h: SD .36 (6.64; .31-.41), CV% 36.85';
    # omega^2 = log(1 + 0.3685^2) = 0.127331
    etaiov_cl_2 ~ fixed(0.127331)
    etaiov_cl_3 ~ fixed(0.127331)
    etaiov_cl_4 ~ fixed(0.127331)
    etaiov_cl_5 ~ fixed(0.127331)
    etaiov_cl_6 ~ fixed(0.127331)
    etaiov_cl_7 ~ fixed(0.127331)
    etaiov_cl_8 ~ fixed(0.127331)

    # =========================================================================
    # Residual unexplained variability (Stoschus 2025 Table 2, 'RUV
    # parameters'). Additive, proportional, and combined models were assessed
    # (Methods 2.2); the proportional model was retained.
    # =========================================================================
    propSd <- 0.082
    label("Proportional residual error (fraction)")
    # Table 2 row 'Proportional = .082 (7.76; .071-.096)'
  })

  model({
    # --- 1. Occasion indicators and interoccasion variability on CL --------
    # Mutually-exclusive indicators over the eight encoded occasions; a
    # record with OCC outside 1..8 zeroes every term and so carries no IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 +
      oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4 +
      oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8

    # --- 2. Individual parameters -----------------------------------------
    # Allometric scaling on ideal body weight centered at the cohort median
    # 68.8 kg (Equations 1 and 2). IIV is lognormal on both V and CL; the
    # single IOV term is additive with the IIV on the log-CL scale.
    vc <- exp(lvc + etalvc) * (IBW / 68.8)^e_ibw_vc
    cl <- exp(lcl + etalcl + iov_cl) * (IBW / 68.8)^e_ibw_cl
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    # --- 3. Micro-constants -----------------------------------------------
    kel <- cl / vc

    # --- 4. ODE system ----------------------------------------------------
    # One compartment with first-order absorption and first-order
    # elimination (Methods 2.2, Results 3.2). Oral doses go to depot,
    # intravenous doses (5-min infusions in the source study) go directly to
    # central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # --- 5. Bioavailability -----------------------------------------------
    # F applies only to the depot (oral) route; intravenous doses entering
    # central are complete by construction, which is what makes F estimable
    # in this cohort.
    f(depot) <- fdepot

    # --- 6. Observation and residual error --------------------------------
    # Dose in mg / volume in L gives mg/L, the unit the paper reports.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
