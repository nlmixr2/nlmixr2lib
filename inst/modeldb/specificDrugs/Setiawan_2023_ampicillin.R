Setiawan_2023_ampicillin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for ampicillin in hospitalised adult patients",
    "on non-ICU general wards in Surabaya, Indonesia, covering a wide range of renal function",
    "(median eGFR 42.2, range 5.9-108.4 mL/min/1.73 m2). Fitted non-parametrically with the NPAG",
    "algorithm in Pmetrics 1.9.7. Clearance, central volume and the two intercompartmental rate",
    "constants k12 (KCP) and k21 (KPC) are primary parameters, each carrying its own",
    "inter-individual variability. Serum creatinine is the only retained covariate and enters",
    "clearance as the inverse ratio CL = CLpop * (1.4 / CREAT), with 1.4 mg/dL the cohort median",
    "(paper Eq. 1); the authors tested eGFR CKD-EPI with and without allometric scaling and both",
    "were outperformed by serum creatinine. The unbound concentration Cu = 0.72 * Cc is exposed",
    "using the fixed 28% protein binding the authors applied in their Monte Carlo %fT>MIC",
    "simulations. Ampicillin and sulbactam were fitted in two separate NPAG runs and are supplied",
    "as two separate model files; see modellib('Setiawan_2023_sulbactam') for the partner",
    "component of the fixed 2:1 ampicillin-sulbactam combination. NOTE: the ampicillin k21 (KPC)",
    "mean printed in Table 2 as 0.17 1/h is a typographical error; 1.17 1/h is used here and the",
    "correction is derived from the paper's own numbers -- see the ini() comment and the vignette",
    "Errata. Residual unexplained variability is carried as fixed(0) because the Pmetrics",
    "assay-error polynomial and the selected lambda/gamma term were never published.",
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
    central     = list(analyte = "ampicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ampicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate retained in the final model. Enters clearance as the simple inverse",
        "ratio CL = CLpop * (1.4 / CREAT) (Setiawan 2023 Eq. 1); the exponent is structurally -1,",
        "not an estimated power. The reference 1.4 mg/dL is the cohort median serum creatinine",
        "(Table 1, median 1.4, range 0.6-6.4 mg/dL), so CLpop = 5.58 L/h is the clearance of a",
        "patient at the cohort median. Because the effect is a raw reciprocal with no exponent and",
        "no lower bound, clearance is unbounded as CREAT approaches zero -- do not simulate outside",
        "the observed 0.6-6.4 mg/dL range. The authors also tested eGFR CKD-EPI on CL both with and",
        "without allometric (0.75) scaling; both were worse than serum creatinine by -2LL and AIC",
        "(Table 2, model-selection block: 306/317 for SeCr vs 312/323 and 314/325 for eGFR).",
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
      notes       = "Screened as a candidate covariate on ampicillin PK (Methods 2.4.2) but not retained in the final model."
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
    n_concentrations = 59L,
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
      "interval at 5, 20, 120 and 240 minutes after injection plus a pre-dose trough; 59 ampicillin",
      "concentrations from 16 patients entered the ampicillin analysis (60 entered the sulbactam",
      "analysis). Total (not unbound) plasma ampicillin was measured by UHPLC-MS/MS, linear from",
      "1 to 200 mg/L with an LLOQ of 1 mg/L. Weight was directly measured in only 8 of 16 patients",
      "(Table 1 footnote b). The reported weight range, 40-82 kg, is numerically identical to the",
      "reported age range, 40-82 years, in Table 1; both are transcribed here as printed.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # All structural values are the MEAN of the NPAG non-parametric parameter
    # distribution reported in Setiawan 2023 Table 2 (ampicillin columns).
    # Pmetrics reports the support-point distribution of the PRIMARY parameters,
    # so CL here is the clearance at the covariate reference (CREAT = 1.4 mg/dL),
    # not a covariate-free typical value.
    # ------------------------------------------------------------------------
    lcl <- log(5.58);  label("Clearance at the reference serum creatinine of 1.4 mg/dL (L/h)")
    # Table 2, ampicillin CL mean = 5.58 L/h (SD 2.57, CV 46.1%, median 5.07); also Abstract and Discussion 4.2
    lvc <- log(12.6);  label("Central volume of distribution (L)")
    # Table 2, ampicillin V mean = 12.6 L (SD 2.16, CV 17.2%, median 12.7); also Abstract "volumes of distribution were 12.6 L"

    lk12 <- log(0.90); label("Transfer rate constant central -> peripheral1 (KCP, 1/h)")
    # Table 2, ampicillin KCP mean = 0.90 (SD 0.78, CV 86.6%, median 0.44). Table 2 labels KCP/KPC
    # "(L/h)"; the same table's own footnote calls them "the rate constant from the central
    # compartment to the peripheral compartment", so the unit is 1/h, not L/h -- see vignette Errata.

    lk21 <- log(1.17); label("Transfer rate constant peripheral1 -> central (KPC, 1/h)")
    # Table 2, ampicillin KPC. CORRECTED from the printed mean of 0.17 to 1.17. The printed row is
    # internally impossible: (i) for a non-negative parameter Markov's inequality forces
    # mean >= median/2 = 1.11/2 = 0.555, and 0.17 < 0.555; (ii) Pmetrics reports CV% = SD/mean, and
    # every other row of Table 2 satisfies it to the printed rounding, but 0.79/0.17 = 465%, not the
    # printed 67.3% -- the CV column back-solves the mean as 0.79/0.673 = 1.174; (iii) with k21 =
    # 0.17 the peripheral volume becomes Vp = k12*Vc/k21 = 66.7 L and Vss = 79 L (1.28 L/kg at the
    # 62 kg cohort median), roughly four-fold above the ~0.3 L/kg known for a small hydrophilic
    # beta-lactam, while k21 = 1.17 gives Vss = 22.3 L (0.36 L/kg). A single lost leading digit
    # reconciles all three; 1.17 also reproduces the printed CV (0.787/1.17 = 67.3%). See the
    # vignette Errata section.

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution rather than a parametric omega; Setiawan 2023 Methods 2.4.2
    # states "The %CV represented the inter-individual variability", so the
    # Table 2 CV% column is carried here as a LOG-NORMAL approximation using
    # omega^2 = log(CV^2 + 1). This is an approximation imposed on a
    # non-parametric distribution -- see vignette Assumptions and deviations.
    # Off-diagonal covariances are not reported and are therefore absent.
    # ------------------------------------------------------------------------
    # Table 2 ampicillin CL CV% = 46.1 -> log(0.461^2 + 1)
    etalcl  ~ 0.192702
    # Table 2 ampicillin V CV% = 17.2 -> log(0.172^2 + 1)
    etalvc  ~ 0.029155
    # Table 2 ampicillin KCP CV% = 86.6 -> log(0.866^2 + 1)
    etalk12 ~ 0.559591
    # Table 2 ampicillin KPC CV% = 67.3 -> log(0.673^2 + 1)
    etalk21 ~ 0.373582

    # ------------------------------------------------------------------------
    # Protein binding. Not estimated: a literature value applied by the authors
    # to convert total to unbound concentrations for the %fT>MIC simulations.
    # ------------------------------------------------------------------------
    fu <- fixed(0.72); label("Fraction of ampicillin unbound in plasma (unitless)")
    # Methods 2.4.4: "Fixed protein binding values of 28% ... for ampicillin ... were used in the dosing simulations"; 1 - 0.28 = 0.72

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
    #    raw inverse ratio of Setiawan 2023 Eq. 1:
    #      CL_AMP = CL_AMP x (1.4 / SeCr)
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

    # 5. Observation. The assay measured TOTAL plasma ampicillin (Methods 2.3),
    #    so Cc is the total concentration and Cu is the unbound concentration
    #    that drives the %fT>MIC targets of Methods 2.4.4.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd)
  })
}
