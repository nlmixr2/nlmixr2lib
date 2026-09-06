Setiawan_2023_sulbactam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for sulbactam in hospitalised adult patients",
    "on non-ICU general wards in Surabaya, Indonesia, covering a wide range of renal function",
    "(median eGFR 42.2, range 5.9-108.4 mL/min/1.73 m2). Fitted non-parametrically with the NPAG",
    "algorithm in Pmetrics 1.9.7. Clearance, central volume and the two intercompartmental rate",
    "constants k12 (KCP) and k21 (KPC) are primary parameters, each carrying its own",
    "inter-individual variability. Serum creatinine is the only retained covariate and enters",
    "clearance as the inverse ratio CL = CLpop * (1.4 / CREAT), with 1.4 mg/dL the cohort median",
    "(paper Eq. 2); the authors tested eGFR CKD-EPI with and without allometric scaling and both",
    "were outperformed by serum creatinine. The unbound concentration Cu = 0.62 * Cc is exposed",
    "using the fixed 38% protein binding the authors applied in their Monte Carlo %fT>MIC",
    "simulations against Acinetobacter baumannii. Ampicillin and sulbactam were fitted in two",
    "separate NPAG runs and are supplied as two separate model files; see",
    "modellib('Setiawan_2023_ampicillin') for the partner component of the fixed 2:1",
    "ampicillin-sulbactam combination. Residual unexplained variability is carried as fixed(0)",
    "because the Pmetrics assay-error polynomial and the selected lambda/gamma term were never",
    "published.",
    sep = " "
  )
  reference <- paste(
    "Setiawan E, Cotta MO, Abdul-Aziz MH, Widjanarko D, Sosilya H, Lukas DL, Wallis SC, Parker S,",
    "Roberts JA. Population pharmacokinetics and dosing simulations of ampicillin and sulbactam in",
    "hospitalised adult patients. Clin Pharmacokinet. 2023;62(4):573-586.",
    "doi:10.1007/s40262-023-01219-5. PMCID: PMC10085897.",
    sep = " "
  )
  vignette <- "Setiawan_2023_ampicillin_sulbactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate retained in the final model. Enters clearance as the simple inverse",
        "ratio CL = CLpop * (1.4 / CREAT) (Setiawan 2023 Eq. 2); the exponent is structurally -1,",
        "not an estimated power. The reference 1.4 mg/dL is the cohort median serum creatinine",
        "(Table 1, median 1.4, range 0.6-6.4 mg/dL), so CLpop = 4.79 L/h is the clearance of a",
        "patient at the cohort median. Because the effect is a raw reciprocal with no exponent and",
        "no lower bound, clearance is unbounded as CREAT approaches zero -- do not simulate outside",
        "the observed 0.6-6.4 mg/dL range. The authors also tested eGFR CKD-EPI on CL both with and",
        "without allometric (0.75) scaling; both were worse than serum creatinine by -2LL and AIC",
        "(Table 2, model-selection block: 249/260 for SeCr vs 256/268 and 252/263 for eGFR).",
        "For the dosing simulations the authors mapped serum creatinine onto creatinine clearance",
        "as CREAT 6, 2, 1.5, 1, 0.7 mg/dL <-> CLcr 10, 20, 30, 70, 100 mL/min/1.73 m2 (Results 3.3).",
        sep = " "
      ),
      source_name        = "SeCr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate on sulbactam PK (Methods 2.4.2) but not retained in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as a candidate covariate (Methods 2.4.2, 'gender') but not retained in the final model."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Methods 2.4.2) but not retained; no parameter in this",
        "model is scaled by body size. Table 1 footnote b records that weight was directly measured",
        "in only 8 of the 16 patients.",
        sep = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.4.2) but not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.4.2) but not retained in the final model; no value is reported in Table 1."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI, BSA-normalized)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Tested on CL both with and without allometric 0.75 scaling and rejected in favour of the",
        "raw serum-creatinine reciprocal (Table 2 model-selection block). Cohort median 42.2,",
        "range 5.9-108.4 mL/min/1.73 m2 (Table 1).",
        sep = " "
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 16L,
    n_studies        = 1L,
    n_concentrations = 60L,
    age_range        = "40-82 years",
    age_median       = "68 years",
    weight_range     = "40-82 kg",
    weight_median    = "62 kg",
    sex_female_pct   = 62.5,
    race_ethnicity   = "Not reported (single-centre Indonesian cohort)",
    disease_state    = paste(
      "Hospitalised adults (>= 18 years) receiving intravenous ampicillin-sulbactam on general",
      "wards (NOT the intensive care unit) of a referral hospital in Surabaya, Indonesia. Patients",
      "on or planned for renal replacement therapy at the time of sampling, and pregnant women,",
      "were excluded. 11 of 16 patients had eGFR CKD-EPI < 60 mL/min/1.73 m2.",
      sep = " "
    ),
    renal_function   = "Serum creatinine median 1.4 mg/dL (range 0.6-6.4); eGFR CKD-EPI median 42.2 mL/min/1.73 m2 (range 5.9-108.4)",
    dose_range       = paste(
      "1000 mg ampicillin + 500 mg sulbactam (1.5 g total) as a ~3-minute intravenous bolus",
      "injection: q8h in 12 patients (75%) and q6h in 4 patients (25%). Only one",
      "ampicillin-sulbactam product was available at the site.",
      sep = " "
    ),
    regions          = "Indonesia (single referral hospital, Surabaya, East Java)",
    notes            = paste(
      "Baseline demographics from Setiawan 2023 Table 1. Serial sampling within one dosing",
      "interval at 5, 20, 120 and 240 minutes after injection plus a pre-dose trough; 60 sulbactam",
      "concentrations from 16 patients entered the sulbactam analysis (59 entered the ampicillin",
      "analysis). Total (not unbound) plasma sulbactam was measured by UHPLC-MS/MS, linear from",
      "0.5 to 100 mg/L with an LLOQ of 0.5 mg/L. Weight was directly measured in only 8 of 16",
      "patients (Table 1 footnote b). The reported weight range, 40-82 kg, is numerically identical",
      "to the reported age range, 40-82 years, in Table 1; both are transcribed here as printed.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # All structural values are the MEAN of the NPAG non-parametric parameter
    # distribution reported in Setiawan 2023 Table 2 (sulbactam columns).
    # Pmetrics reports the support-point distribution of the PRIMARY parameters,
    # so CL here is the clearance at the covariate reference (CREAT = 1.4 mg/dL),
    # not a covariate-free typical value.
    # ------------------------------------------------------------------------
    lcl <- log(4.79);  label("Clearance at the reference serum creatinine of 1.4 mg/dL (L/h)")
    # Table 2, sulbactam CL mean = 4.79 L/h (SD 2.05, CV 42.8%, median 4.47); also Abstract and Discussion 4.2
    lvc <- log(15.36); label("Central volume of distribution (L)")
    # Abstract: "the volumes of distribution were 12.6 L and 15.36 L". Table 2 rounds this to 15.4;
    # the extra digit is confirmed by the CV column, since 4.73/15.36 = 30.8% (printed) while
    # 4.73/15.4 = 30.7%. SD 4.73, CV 30.8%, median 14.6.

    lk12 <- log(0.42); label("Transfer rate constant central -> peripheral1 (KCP, 1/h)")
    # Table 2, sulbactam KCP mean = 0.42 (SD 0.41, CV 96.8%, median 0.32). Table 2 labels KCP/KPC
    # "(L/h)"; the same table's own footnote calls them "the rate constant from the central
    # compartment to the peripheral compartment", so the unit is 1/h, not L/h -- see vignette Errata.
    lk21 <- log(0.69); label("Transfer rate constant peripheral1 -> central (KPC, 1/h)")
    # Table 2, sulbactam KPC mean = 0.69 (SD 0.56, CV 81.7%, median 0.79)

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution rather than a parametric omega; Setiawan 2023 Methods 2.4.2
    # states "The %CV represented the inter-individual variability", so the
    # Table 2 CV% column is carried here as a LOG-NORMAL approximation using
    # omega^2 = log(CV^2 + 1). This is an approximation imposed on a
    # non-parametric distribution -- see vignette Assumptions and deviations.
    # Off-diagonal covariances are not reported and are therefore absent.
    # ------------------------------------------------------------------------
    # Table 2 sulbactam CL CV% = 42.8 -> log(0.428^2 + 1)
    etalcl  ~ 0.168209
    # Table 2 sulbactam V CV% = 30.8 -> log(0.308^2 + 1)
    etalvc  ~ 0.090630
    # Table 2 sulbactam KCP CV% = 96.8 -> log(0.968^2 + 1)
    etalk12 ~ 0.661153
    # Table 2 sulbactam KPC CV% = 81.7 -> log(0.817^2 + 1)
    etalk21 ~ 0.511319

    # ------------------------------------------------------------------------
    # Protein binding. Not estimated: a literature value applied by the authors
    # to convert total to unbound concentrations for the %fT>MIC simulations.
    # ------------------------------------------------------------------------
    fu <- fixed(0.62); label("Fraction of sulbactam unbound in plasma (unitless)")
    # Methods 2.4.4: "Fixed protein binding values of ... 38% for ... sulbactam ... were used in the dosing simulations"; 1 - 0.38 = 0.62

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported. Methods 2.4.1 states only
    # that "Both lambda and gamma error models were tested"; neither the selected
    # error model nor the Pmetrics assay-error polynomial coefficients (C0-C3)
    # appear in the paper, and there is no supplement on disk carrying them.
    # Carried as fixed(0) rather than invented -- see vignette Errata.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })

  model({
    # 1. Individual parameters. Serum creatinine acts on clearance only, as the
    #    raw inverse ratio of Setiawan 2023 Eq. 2:
    #      CL_SUL = CL_SUL x (1.4 / SeCr)
    #    with 1.4 mg/dL the cohort median serum creatinine (Results 3.2).
    cl <- exp(lcl + etalcl) * (1.4 / CREAT)
    vc <- exp(lvc + etalvc)

    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # 2. Intercompartmental clearance and peripheral volume, derived from the
    #    primary rate constants. The ODEs below MUST be driven by q / vp rather
    #    than by kcp / kpc directly: rxSolve() defaults to useLinCmt = TRUE and,
    #    when the peripheral transfer is written straight from stored
    #    micro-constants, that rewrite can silently drop peripheral1 and solve a
    #    one-compartment model. Routing through q and vp keeps the closed-form
    #    and ODE solvers in agreement (the vignette asserts this identity).
    q  <- kcp * vc
    vp <- q / kpc

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment ODE system. Dosing is intravenous into central: either a
    #    ~3-minute bolus injection (as administered in the study) or a 4-h
    #    prolonged infusion (as simulated in Methods 2.4.4).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation. The assay measured TOTAL plasma sulbactam (Methods 2.3),
    #    so Cc is the total concentration and Cu is the unbound concentration
    #    that drives the %fT>MIC targets of Methods 2.4.4.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd)
  })
}
