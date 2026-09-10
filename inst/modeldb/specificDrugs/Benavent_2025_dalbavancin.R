Benavent_2025_dalbavancin <- function() {
  description <- paste(
    "One-compartment population PK model with intravenous infusion and first-order linear",
    "elimination for TOTAL plasma dalbavancin after a SINGLE 1500 mg dose given as sequencing",
    "therapy to elderly adults with chronic Gram-positive prosthetic joint infection managed by",
    "two-stage exchange with vancomycin/gentamicin-loaded cement spacers. Interindividual",
    "variability was estimated on both clearance and volume; residual error is proportional. No",
    "covariate improved the fit, so the structural model IS the final model: age, sex, height,",
    "body weight, body mass index, glomerular filtration rate, baseline creatinine clearance and",
    "same-day serum albumin were all screened and rejected. The model returns TOTAL dalbavancin",
    "Cc; the paper's PK/PD target attainment analysis derives unbound exposure by scaling Cc",
    "with a free fraction swept over four theoretical protein-binding scenarios (93, 95, 97 and",
    "99 percent), so no single free fraction is packaged here. Clearance is about 30 percent",
    "lower than reported for healthy volunteers and younger patients with acute infection",
    "(0.036 vs 0.050 L/h), which the authors attribute to the cohort's older age, physiological",
    "albumin and lower creatinine clearance.",
    sep = " "
  )
  reference <- paste(
    "Benavent E, Lora-Tamayo J, Ulldemolins M, Pons-Oltra P, Gregoire M,",
    "Mancheno-Losa M, Hernandez-Jimenez P, Melendez-Carmona MA, Casals V,",
    "Roberts JA, Rigo-Bonnin R, Murillo O. Efficacy, safety, and population",
    "pharmacokinetics of a single 1500 mg dose of dalbavancin for short-term",
    "therapy in patients with chronic prosthetic joint infections.",
    "Antimicrob Agents Chemother. 2025;69(12):e00773-25.",
    "doi:10.1128/aac.00773-25",
    sep = " "
  )
  vignette <- "Benavent_2025_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Benavent 2025 Results ("Population pharmacokinetic
  # analysis"): "a one-compartment model with intravenous infusion and
  # first-order linear elimination" fitted to TOTAL dalbavancin plasma
  # concentrations measured by UHPLC-MS/MS (Methods, "Sample handling and
  # bioanalysis").
  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # The final model carries NO covariates, so covariateData is empty. Every
  # covariate the paper screened is documented in covariatesDataExcluded below.
  covariateData <- list()

  # Covariates screened by Benavent 2025 but NOT retained. Documentation only
  # -- none of these is referenced in model(). Benavent 2025 Methods
  # ("Structural model building, covariate analysis, and model evaluation"):
  # "the effects of the following covariates on dalbavancin PK parameters were
  # evaluated on the structural model: age, gender, height, body weight, body
  # mass index, glomerular filtration rate (19), baseline creatinine clearance
  # (CrCL) (20), and albumin serum concentrations on the day of sampling."
  # Results ("Population pharmacokinetic analysis and Monte Carlo dosing
  # simulations"): "The covariate analysis did not result in model improvements
  # for which the structural model is the final model, summarized in Table 2."
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on the structural model but not retained. Median 75.5 years, interquartile interval 69-79 years (Benavent 2025 Results, first paragraph)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Printed as 'gender' in Benavent 2025 Methods; screened but not retained. 70% of the 20 enrolled patients were female (n = 14; Results, first paragraph). Recorded here on the canonical female-indicator convention."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Benavent 2025 Methods). The paper does not print a cohort height summary in the main text."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained (Benavent 2025 Methods). The paper does not print a cohort weight summary in the main text, and no allometric scaling is applied; the packaged V and CL are absolute values for this elderly cohort."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Benavent 2025 Methods). The paper does not print a cohort body-mass-index summary in the main text."
    ),
    CRCL = list(
      description = "Renal function: BSA-normalized glomerular filtration rate and, separately, baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Benavent 2025 screened TWO renal-function covariates as separate",
        "candidates -- glomerular filtration rate (their reference 19) and",
        "baseline creatinine clearance by the Cockcroft-Gault-style equation of",
        "their reference 20 -- and retained neither. Both fold onto the single",
        "canonical CRCL column, whose register entry explicitly spans",
        "creatinine-based estimated GFR and measured creatinine clearance. The",
        "cohort had uniformly preserved renal function (median glomerular",
        "filtration rate 90 mL/min, interquartile interval 75.8-96.3 mL/min;",
        "Results, first paragraph), which is the likely reason no renal effect",
        "was identifiable. The paper reports the GFR summary without stating",
        "whether it is BSA-normalized, so the units above are given as the",
        "printed mL/min. The Discussion nevertheless attributes this cohort's",
        "30% lower clearance versus healthy volunteers partly to its lower",
        "creatinine clearance, i.e. the effect is believed real but was not",
        "estimable in 18 patients spanning a narrow renal range."
      )
    ),
    ALB = list(
      description = "Serum albumin on the day of sampling",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Benavent 2025 Methods, which specifies 'albumin serum concentrations on the day of sampling', i.e. a time-varying candidate). The cohort had normal albumin throughout: median 45 g/L, interquartile interval 38-47 g/L (Results, first paragraph). The Discussion attributes part of this cohort's lower clearance to its physiological albumin relative to the comparator populations."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_median     = "75.5 years",
    age_range      = "interquartile interval 69-79 years",
    sex_female_pct = 70,
    race_ethnicity = "Not reported; two-centre Spanish cohort.",
    disease_state  = paste(
      "Chronic prosthetic joint infection caused by low-virulence Gram-positive",
      "bacteria susceptible to dalbavancin, managed by two-stage prosthetic",
      "exchange with a vancomycin- plus gentamicin-loaded cement spacer.",
      "Affected joints: hip 55% (n = 11), knee 30% (n = 6), shoulder 10%",
      "(n = 2), ankle 5% (n = 1). Isolates (Table 1, 24 isolates): coagulase-",
      "negative staphylococci 70.8% (Staphylococcus epidermidis 50%,",
      "S. lugdunensis 4.2%, other CoNS 16.7%), Cutibacterium acnes 25%,",
      "Enterococcus faecalis 4.2%; four polymicrobial infections."
    ),
    dose_range     = paste(
      "A single 1500 mg intravenous dose of dalbavancin, given as a 30 min",
      "short infusion at Hospital 12 de Octubre or a 2 h extended infusion at",
      "Hospital Universitari de Bellvitge, after a median 11.5 days",
      "(interquartile interval 10-16) of prior intravenous therapy",
      "(vancomycin 75%, daptomycin 25%; switched to an oxazolidinone before",
      "dalbavancin in 30% of cases)."
    ),
    regions        = "Spain (Hospital Universitari de Bellvitge, Barcelona; Hospital Universitario 12 de Octubre, Madrid)",
    renal_function = "Uniformly preserved: median glomerular filtration rate 90 mL/min, interquartile interval 75.8-96.3 mL/min (Results, first paragraph). No patient had renal impairment, and no change in renal function occurred during follow-up.",
    notes          = paste(
      "Retrospective, observational, two-centre clinical and PK study run",
      "1 January 2022 to 31 May 2023 (ethics reference EOM017/23). Twenty",
      "patients were enrolled and reported for the efficacy and safety",
      "endpoints; the population PK model was fitted to TOTAL plasma",
      "dalbavancin from 18 of them (Results, 'Population pharmacokinetic",
      "analysis'), each contributing 1-3 concentrations, so the analysis",
      "dataset holds between 18 and 54 observations -- the exact count is given",
      "only in supplementary Table S2, which was not retrievable (see the",
      "vignette Errata). Sampling was opportunistic at weekly-to-biweekly",
      "outpatient visits through week 4 post-dose rather than on a fixed",
      "schedule. Total dalbavancin was measured by UHPLC-MS/MS with a lower",
      "limit of quantification of 1.0 mg/L over a 1.0-250 mg/L measuring",
      "interval, imprecision <= 8.6% and absolute relative bias <= 7.3%.",
      "Estimation was SAEM in Monolix 2024R1; the final model was checked by",
      "prediction-corrected VPC (500 simulations) and by nonparametric",
      "bootstrap (n = 1000). Serum albumin was normal throughout (median 45",
      "g/L, interquartile interval 38-47 g/L)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Benavent 2025 Table 2 ("Population
    # pharmacokinetic estimates for dalbavancin"), 'Fixed effects' block,
    # 'Estimate (% RSE) [Shrinkage %]' column. Both are absolute values
    # for this cohort: the final model carries no covariate and no
    # allometric size term, so there is no reference weight.
    # ------------------------------------------------------------------
    # Table 2, 'V (L)' = 17.9 (6.2% RSE) [18.1% shrinkage]; bootstrap median 18.0 (95% CI 16.0-20.3)
    lvc <- log(17.9); label("Central volume of distribution (L)")
    # Table 2, 'CL (L/h)' = 0.036 (8.2% RSE) [12.6% shrinkage]; bootstrap median 0.037 (95% CI 0.031-0.043)
    lcl <- log(0.036); label("Clearance (L/h)")

    # ------------------------------------------------------------------
    # Interindividual variability -- Benavent 2025 Table 2, 'Random
    # effects' block. Methods: "All individual parameters were assumed to
    # be log-normally distributed. Inter-individual variability (IIV) was
    # described using an exponential model." Monolix reports omega on the
    # STANDARD DEVIATION scale and Table 2 says so explicitly in the row
    # labels and in footnote a ('SD, standard deviation'), so the
    # variances below are the printed SDs squared.
    # ------------------------------------------------------------------
    etalvc ~ 0.200^2  # Table 2, 'IIV V (SD)' = 0.200 (27.4% RSE); bootstrap median 0.200 (95% CI 0.061-0.280). Squared to a variance.
    etalcl ~ 0.290^2  # Table 2, 'IIV CL (SD)' = 0.290 (20.8% RSE); bootstrap median 0.280 (95% CI 0.150-0.380). Squared to a variance.

    # ------------------------------------------------------------------
    # Residual error -- Benavent 2025 Table 2, 'Residual variability'
    # block. Methods tested constant, proportional and combined error
    # models; Results reports "the residual error was modeled as
    # proportional", and Table 2 prints the single Monolix proportional
    # coefficient b, which maps directly onto propSd.
    # ------------------------------------------------------------------
    # Table 2, 'b (proportional)' = 0.120 (28.6% RSE); bootstrap median 0.110 (95% CI 0.053-0.170)
    propSd <- 0.120; label("Proportional residual error (fraction)")
  })

  model({
    # Individual parameters -- log-normal, per Methods.
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # Micro-constant
    kel <- cl / vc

    # One-compartment disposition. Dalbavancin is given intravenously, so
    # the dose enters `central` directly; the infusion duration is an
    # event-table property (30 min or 2 h in this study) rather than a
    # model parameter.
    d/dt(central) <- -kel * central

    # Observation: TOTAL plasma dalbavancin. The paper's PK/PD analysis
    # multiplies this by a free fraction of 1 - protein binding, swept
    # over 93/95/97/99% protein binding; that scaling belongs to the
    # simulation scenario, not to the fitted model.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
