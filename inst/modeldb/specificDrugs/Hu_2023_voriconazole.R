Hu_2023_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and elimination for intravenous and oral voriconazole in Chinese paediatric haematology patients with invasive fungal infection (Hu 2023); CYP2C19 metabolizer phenotype is the only retained covariate and enters multiplicatively on clearance, with the normal-metabolizer phenotype as the implicit reference."
  reference <- "Hu L, Huang S, Huang Q, Huang J, Feng Z, He G. Population pharmacokinetics of voriconazole and the role of CYP2C19 genotype on treatment optimization in pediatric patients. PLoS ONE. 2023;18(9):e0288794. doi:10.1371/journal.pone.0288794"
  vignette <- "Hu_2023_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer, CYP2C19*1/*1; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = paste(
        "Hu 2023 Table 2 reports the covariate as a multiplicative fraction on clearance",
        "('IM on CL' = 0.582), i.e. CL_IM = 7.35 * 0.582 = 4.28 L/h, with the normal",
        "metabolizer (NM) phenotype as the paper's own reference category -- no",
        "reparameterization was needed. IM genotypes pooled by Hu 2023 (Methods,",
        "'Measurement of VRC trough plasma concentrations and CYP2C19 phenotype') were",
        "CYP2C19*1/*2 and *1/*3. 43 of 91 subjects (47.3%) were IM. Note that Hu 2023",
        "pooled only *2 and *3 into the reduced-function set; no CYP2C19*17 allele was",
        "observed in the cohort, so no ultrarapid-metabolizer stratum exists and the",
        "NM reference is unambiguous."
      ),
      source_name        = "CYP2C19 phenotype (IM)"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer, CYP2C19*1/*1; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = paste(
        "Companion to CYP2C19_IM. Hu 2023 Table 2 reports 'PM on CL' = 0.381, i.e.",
        "CL_PM = 7.35 * 0.381 = 2.80 L/h. PM genotypes pooled by Hu 2023 were",
        "CYP2C19*2/*2, *2/*3 and *3/*3. 11 of 91 subjects (12.1%) were PM.",
        "The monotone ordering NM (1.0) > IM (0.582) > PM (0.381) matches the",
        "monotone decrease in observed dose-normalized trough concentration across",
        "the three phenotypes (Hu 2023 Fig 1A and Results, 'CYP2C19 phenotypes')."
      ),
      source_name        = "CYP2C19 phenotype (PM)"
    )
  )

  # Covariates that Hu 2023 screened in the stepwise covariate search but did
  # NOT retain in the final model. Documented here for provenance; they are
  # deliberately absent from model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Hu 2023 Methods, 'Covariate model') but not",
        "retained in the final model, so CL and Vc are absolute (not weight-scaled)",
        "values. Cohort median 31.0 kg (range 9.5-85.0), Hu 2023 Table 1. Weight still",
        "matters operationally because the paper's dosing regimens are expressed in",
        "mg/kg; it enters through the dose amount, not through the PK parameters."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained (Hu 2023 Methods, 'Covariate model'). Cohort median 10 years, range 2-14 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Hu 2023 Methods, 'Covariate model'). 33 of 91 subjects (36.3%) were female (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a liver-function indicator but not retained. Cohort median 35.20 g/L, range 20.30-49.00 (Hu 2023 Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a liver-function indicator but not retained. Cohort median 8.25 umol/L, range 2.10-105.20 (Hu 2023 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a liver-function indicator but not retained. Cohort median 31.70 U/L, range 4.60-346.40 (Hu 2023 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a liver-function indicator but not retained. Cohort median 31.05 U/L, range 3.70-255.90 (Hu 2023 Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a kidney-function indicator but not retained. Cohort median 48.00 umol/L, range 27.00-252.00 (Hu 2023 Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as a kidney-function indicator but not retained. Cohort median 3.36 mmol/L, range 0.98-11.72 (Hu 2023 Table 1)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as 'combination therapy' but not retained (Hu 2023 Methods,",
        "'Covariate model'). 55 of 91 subjects (60.4%) received a PPI: omeprazole",
        "(n = 15), pantoprazole (n = 34), lansoprazole (n = 6) (Hu 2023 Table 1",
        "footnote a). Hu 2023 pooled the three PPIs into a single indicator rather",
        "than modelling them separately, and reports no minimum-duration threshold",
        "for the flag, so it is an ever-versus-never subject-level indicator."
      )
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic glucocorticoid therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as 'combination therapy' but not retained. 44 of 91 subjects (48.4%)",
        "received a glucocorticoid: methylprednisolone (n = 17), dexamethasone",
        "(n = 17), prednisone (n = 10) (Hu 2023 Table 1 footnote b). Hu 2023 pooled",
        "the three glucocorticoids into a single ever-versus-never indicator."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 91L,
    n_studies      = 1L,
    n_observations = 210L,
    age_range      = "2-14 years",
    age_median     = "10 years",
    age_strata     = c(UnderSix_pct = 23.1, SixToTwelve_pct = 35.2, OverTwelve_pct = 41.8),
    weight_range   = "9.5-85.0 kg",
    weight_median  = "31.0 kg",
    sex_female_pct = 36.3,
    race_ethnicity = c(Chinese = 100),
    cyp2c19_phenotype = c(NM_pct = 40.7, IM_pct = 47.3, PM_pct = 12.1, UM_pct = 0.0),
    disease_state  = paste(
      "Paediatric patients (14 years of age or younger) with malignant haematological",
      "disease and invasive fungal infection treated with voriconazole. Acute",
      "lymphoblastic leukaemia 57.1%, acute myeloid leukaemia 14.3%, lymphoma 12.1%,",
      "thalassaemia 7.7%, aplastic anaemia 3.3%, other 5.5%. Invasive fungal infection",
      "was proven in 6.6%, probable in 18.7% and possible in 74.7%; lung was the most",
      "common site of infection (60.4%). Treatment indication was therapeutic in 17.6%,",
      "empirical in 53.8% and prophylactic in 28.6%."
    ),
    dose_range     = paste(
      "Voriconazole dosed per the manufacturer's paediatric labelling: intravenous",
      "loading 9 mg/kg and maintenance 8 mg/kg twice daily; oral maintenance 9 mg/kg",
      "twice daily with no oral loading dose. Most patients received no loading dose",
      "because they were dosed orally. Maintenance dose was subsequently adjusted by",
      "the treating physician on the basis of therapeutic drug monitoring, efficacy or",
      "adverse drug reactions. Route: oral only 81.3%, intravenous only 8.8%,",
      "intravenous switched to oral 6.6%, oral switched to intravenous 3.3%. Median",
      "duration of voriconazole use 15 days (range 3-148)."
    ),
    regions        = "Single centre: department of paediatric haematology, Xiangya Hospital, Central South University, Changsha, Hunan, China.",
    notes          = paste(
      "Retrospective observational study of medical records collected 1 January 2018 to",
      "31 December 2021. 210 steady-state trough concentrations from 91 children; median",
      "1 measurement per patient (range 1-21 per Table 1; the Results text states range",
      "1-19). All samples were troughs drawn 30 minutes before the next dose, so the",
      "dataset carries essentially no information about the absorption or distribution",
      "phase -- ka was therefore fixed to a literature value and Vc is estimated with",
      "very large interindividual variability. Median observed trough 1.23 mg/L (range",
      "0.02-8.58); 52.9% of troughs were within the 1.0-5.5 mg/L target range, 40.9%",
      "subtherapeutic and 6.2% supratherapeutic. Voriconazole was assayed by HPLC over",
      "0.02-19.60 mg/L. CYP2C19 phenotype was assigned by DNA microarray. Model built in",
      "NONMEM 7.5 with FOCE-I; final-model parameter estimates and 1000-replicate",
      "bootstrap per Hu 2023 Table 2."
    )
  )

  ini({
    # Reference subject for the typical-value equations is a CYP2C19 normal
    # metabolizer (CYP2C19_IM = CYP2C19_PM = 0). Hu 2023 retained no
    # continuous covariate, so CL and Vc are absolute values that are NOT
    # weight-scaled; see covariatesDataExcluded$WT.

    # Absorption. Hu 2023 Methods, 'Population pharmacokinetic modeling':
    # "Therefore, Ka was fixed at 1.19 h-1 as previously determined by Friberg
    # et al. [18]." Table 2 records it as "1.19 (fixed)" in both the base and
    # the final model, so it is encoded with fixed().
    lka <- fixed(log(1.19)); label("Absorption rate constant (1/h)")  # Hu 2023 Table 2 final model: ka = 1.19 (fixed)

    # Clearance for the CYP2C19 normal-metabolizer reference.
    lcl <- log(7.35); label("Clearance for the CYP2C19 normal-metabolizer reference (L/h)")  # Hu 2023 Table 2 final model: CL = 7.35 L/h (RSE 15%; bootstrap median 7.31, 95% CI 3.5-11.7)

    # Central volume of distribution. No covariate was retained on Vc.
    lvc <- log(376); label("Central volume of distribution (L)")  # Hu 2023 Table 2 final model: Vc = 376 L (RSE 21%; bootstrap median 393, 95% CI 206-854)

    # Oral bioavailability. Hu 2023 Table 2 reports F as a percentage (52.2%);
    # the Discussion restates it as "the typical value of VRC bioavailability
    # was 52.2%".
    lfdepot <- log(0.522); label("Oral bioavailability (fraction)")  # Hu 2023 Table 2 final model: F = 52.2% (RSE 15%; bootstrap median 51.4, 95% CI 28.2-82.8)

    # CYP2C19 phenotype effects on clearance. Hu 2023 Table 2 reports these as
    # multiplicative fractions relative to the NM reference ("IM on CL" 0.582,
    # "PM on CL" 0.381), so CL_IM = 7.35 * 0.582 = 4.28 L/h and
    # CL_PM = 7.35 * 0.381 = 2.80 L/h. They are converted to log-scale additive
    # shifts here so that they compose with the exponential IIV in the usual
    # way; the paper's own fractions remain visible inside log().
    #
    # The multiplicative-fraction reading (rather than an exponent, which would
    # give exp(0.582) = 1.79 and make intermediate metabolizers CLEAR FASTER
    # than normal metabolizers) is confirmed three ways:
    #   1. Direction: reduced CYP2C19 function must lower, not raise, CL.
    #   2. Observed data (Hu 2023 Results, 'CYP2C19 phenotypes'): median trough
    #      0.77 / 1.45 / 2.19 mg/L for NM / IM / PM, i.e. IM:NM = 1.88 and
    #      PM:NM = 2.84, bracketing 1/0.582 = 1.72 and 1/0.381 = 2.62 from
    #      above as steady-state accumulation makes trough supra-proportional
    #      to 1/CL.
    #   3. The paper's own simulations (Table 3, 9 mg/kg oral): median predicted
    #      trough 0.88 / 1.71 / 2.48 mg/L, ratios 1.94 and 2.82, same pattern.
    e_im_cl <- log(0.582); label("Log-scale CL shift for CYP2C19_IM vs normal metabolizer (unitless)")  # Hu 2023 Table 2 final model: IM on CL = 0.582 (RSE 12%; bootstrap median 0.611, 95% CI 0.343-0.887)
    e_pm_cl <- log(0.381); label("Log-scale CL shift for CYP2C19_PM vs normal metabolizer (unitless)")  # Hu 2023 Table 2 final model: PM on CL = 0.381 (RSE 14%; bootstrap median 0.387, 95% CI 0.200-0.826)

    # Interindividual variability. Hu 2023 Methods states the exponential form
    # explicitly: "The inter-individual variability (IIV) of PK parameters was
    # described by an exponential error model: Pi = Ptv * exp(eta_i) ... eta_i
    # is a normal distribution with a mean of 0 and a variance of omega^2."
    # Table 2's IIV block reports omega_CL = 24.7 and omega_Vc = 233.2 for the
    # final model, without stating whether those are omega itself expressed as
    # a percent or a CV% requiring back-transformation. The standard NONMEM
    # reporting convention CV% ~= sqrt(omega^2) * 100 is used here, giving
    # internal variances 0.247^2 for CL and 2.332^2 for Vc.
    #
    # That reading is confirmed against the paper's OWN simulation output
    # rather than assumed: re-running Hu 2023 Table 3 (28 days of q12h dosing,
    # the published cohort's body weights from S1 Raw data) reproduces the
    # twelve 9 mg/kg mean/median cells with a mean absolute error of 4.8%
    # (max 10.4%). The two competing readings are decisively worse -- the
    # exact-lognormal back-transformation omega^2 = log(1 + CV^2)
    # (omega_Vc = 1.365) gives a mean error of 32.1%, and reading the tabulated
    # numbers as omega^2 * 100 (omega_Vc = 1.527) gives 29.6%. Both overshoot
    # by roughly a factor of six. See the vignette for the full comparison.
    #
    # The omega_Vc value is genuinely extreme (bootstrap 95% CI 148-263, so it
    # is not a typographic slip). It reflects a trough-only dataset in which Vc
    # is very poorly identified; it is transcribed as published, not tempered.
    etalcl ~ 0.061009  # Hu 2023 Table 2 final model: omega_CL = 24.7 (bootstrap median 19.8, 95% CI 0.242-62.2); var = 0.247^2
    etalvc ~ 5.438224  # Hu 2023 Table 2 final model: omega_Vc = 233.2 (bootstrap median 223, 95% CI 148-263); var = 2.332^2

    # Residual error. Hu 2023 Methods gives the residual model as an additive
    # error on the natural-log scale ("Cobs,ij = ln(Cpred,ij) + epsilon"; the
    # left-hand ln is dropped in the typeset equation but the intended form is
    # ln(Cobs) = ln(Cpred) + epsilon), which is the proportional error model in
    # nlmixr2's linear concentration space. Table 2 labels the corresponding row
    # "Proportion residual error (%)" and reports 94.7% for the final model.
    propSd <- 0.947; label("Proportional residual error (fraction)")  # Hu 2023 Table 2 final model: proportional residual error = 94.7% (RSE 9%; bootstrap median 93.9, 95% CI 74.7-117.0)
  })

  model({
    # Individual parameters. CYP2C19 normal metabolizer is the implicit
    # reference (CYP2C19_IM = CYP2C19_PM = 0).
    ka     <- exp(lka)
    cl     <- exp(lcl + e_im_cl * CYP2C19_IM + e_pm_cl * CYP2C19_PM + etalcl)
    vc     <- exp(lvc + etalvc)
    fdepot <- exp(lfdepot)

    # One-compartment disposition with first-order absorption and first-order
    # elimination. Intravenous doses are given directly into central and
    # bypass the depot compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Oral bioavailability applies only to doses entering via depot.
    f(depot) <- fdepot

    # Observation. Vc is in L and doses are in mg, so Cc is in mg/L, the unit
    # used for every voriconazole trough concentration reported by Hu 2023.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
