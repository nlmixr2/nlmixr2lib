Assmus_2025_benznidazole <- function() {
  description <- paste(
    "One-compartment population PK model with transit-compartment",
    "absorption (one transit compartment) and linear elimination for",
    "oral benznidazole in adults with chronic indeterminate Chagas",
    "disease (Assmus 2025; BENDITA trial NCT03378661, n = 175 across",
    "six active treatment arms). Concentrations were measured in dried",
    "blood spots, so CL/F = 1.30 L/h and V/F = 31.6 L are apparent",
    "whole-blood values for a 65 kg woman. Mean transit time MTT = 0.75",
    "h. Body weight is scaled allometrically on CL/F (exponent 0.75)",
    "and V/F (exponent 1) about a 65 kg reference. Co-administration of",
    "fosravuconazole increases CL/F by 17.7%, and relative oral",
    "bioavailability is 12.9% lower in men than in women; the authors",
    "judged both effects not clinically relevant. Inter-individual",
    "variability is on MTT (61.6% CV), CL/F (18.9% CV) and relative",
    "bioavailability (10.2% CV); none was retained on V/F. Residual",
    "error is additive on the log-transformed concentration scale",
    "(variance 0.076), i.e. log-normal on the arithmetic scale."
  )
  reference <- paste(
    "Assmus F, Cruz C, Watson JA, White NJ, Adehin A, Hoglund RM,",
    "Blum de Oliveira B, Barreira F, Scandale I, Tarning J.",
    "Population pharmacokinetic-pharmacodynamic analysis of",
    "benznidazole monotherapy and combination therapy with",
    "fosravuconazole in chronic Chagas disease (BENDITA).",
    "PLoS Neglected Tropical Diseases. 2025;19(9):e0013522.",
    "doi:10.1371/journal.pntd.0013522.",
    sep = " "
  )
  vignette <- "Assmus_2025_benznidazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling applied a priori (not estimated) on CL/F",
        "(exponent 0.75) and V/F (exponent 1.0), standardized to 65 kg",
        "-- the cohort median (Assmus 2025 Methods, Population",
        "pharmacokinetic analysis (i); Table 1). Inclusion of the",
        "allometric function gave dOFV = -34.7 (Assmus 2025 Results,",
        "Pharmacokinetics). Enrolment was restricted to 50-80 kg, so",
        "extrapolation outside that band is unsupported by these data."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (male) is the register-canonical reference, but Assmus 2025",
        "reports typical values for a FEMALE subject; the male effect is",
        "therefore applied via (1 - SEXF) to preserve the published",
        "structural estimates verbatim."
      ),
      notes              = paste(
        "Retained on relative oral bioavailability F only (dOFV = -15.5);",
        "men had 12.9% lower F than women (Assmus 2025 Table 2 'Sex",
        "effect on F (reference: female)'). The PK analysis dataset was",
        "68% female (56/175 male). The authors judged the effect not",
        "clinically relevant but retained it in the final model",
        "(Assmus 2025 Discussion, Pharmacokinetic properties)."
      ),
      source_name        = "SEX"
    ),
    CONMED_FOSRAVUCONAZOLE = list(
      description        = "Concomitant fosravuconazole (E1224) administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (benznidazole monotherapy)",
      notes              = paste(
        "1 = subject randomised to a benznidazole + fosravuconazole",
        "combination arm (150 mg benznidazole once daily for 4 weeks +",
        "E1224, or 300 mg benznidazole once weekly for 8 weeks + E1224).",
        "Fosravuconazole was dosed 300 mg once daily for 3 days then 300",
        "mg once weekly for 8 weeks (Assmus 2025 Methods, Treatment",
        "regimens). Encoded as a time-fixed arm-level indicator because",
        "the paper models it as a single categorical covariate on CL/F",
        "(+17.7%, dOFV = -15.7; Table 2), not as a time-varying",
        "co-medication flag."
      ),
      source_name        = "E1224"
    )
  )

  # Covariates screened in the step-wise covariate model (scm, PsN) that
  # were not retained after backward elimination (p = 0.001). Listed for
  # provenance only; not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Preselected covariate, tested by scm; not retained (Assmus 2025 Methods, Population pharmacokinetic analysis (i); Results, Pharmacokinetics). PK dataset median 34 years (range 18-50).",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Marker of kidney function; estimated from blood creatinine by the Cockcroft and Gault equation. Tested by scm; not retained. PK dataset median 102 mL/min (range 54.8-180).",
      source_name        = "CLCR"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested by scm; not retained. Of specific interest here because concentrations were measured in dried blood spots, where hematocrit can influence the blood-to-DBS relationship.",
      source_name        = "HCT"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Marker of liver function; tested by scm, not retained. PK dataset median 22 IU/L (range 9-46).",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Marker of liver function; tested by scm, not retained. PK dataset median 23 IU/L (range 12-43.7).",
      source_name        = "AST"
    ),
    ALKP = list(
      description        = "Alkaline phosphatase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Marker of liver function; tested by scm, not retained. PK dataset median 195 IU/L (range 50-339).",
      source_name        = "ALP"
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Marker of liver function; tested by scm, not retained. PK dataset median 13.7 umol/L (range 6.84-20.5).",
      source_name        = "TBILI"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "benznidazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "benznidazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "benznidazole", units = "mg",
      specimen = "whole blood", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 175L,
    n_studies      = 1L,
    n_observations = 986L,
    age_range      = "18-50 years",
    age_median     = "34 years",
    weight_range   = "50-80 kg",
    weight_median  = "65 kg",
    sex_female_pct = 68,
    race_ethnicity = "All participants were Bolivian.",
    disease_state  = paste(
      "Adults (18-50 years, 50-80 kg) with chronic indeterminate Chagas",
      "disease, confirmed by serological testing and a positive",
      "qualitative PCR result. Subjects with chronic health conditions",
      "and with signs or symptoms of the chronic cardiac or digestive",
      "form of Chagas disease were excluded."
    ),
    dose_range     = paste(
      "Six active oral benznidazole arms (Abarax, Laboratorios ELEA):",
      "150 mg twice daily for 8, 4 or 2 weeks; 150 mg once daily for 4",
      "weeks alone or with fosravuconazole; and 300 mg once weekly",
      "(split into two doses) for 8 weeks with fosravuconazole.",
      "Fosravuconazole 300 mg once daily for 3 days then 300 mg once",
      "weekly for 8 weeks. A seventh placebo arm contributed no PK data."
    ),
    regions        = "Bolivia (Cochabamba, Tarija and Sucre).",
    sampling       = paste(
      "Sparse dried blood spot (DBS) sampling: predose and at follow-up",
      "visits on days 1-3 and weeks 2, 3, 4, 6 and 10 (9 protocol time",
      "points per patient). HPLC-MS/MS assay, LLOQ 50 ng/mL. 986",
      "concentrations from 175 subjects, 969 above LLOQ (1.7% BQL).",
      "Data were censored at 120 h after dose (about 10 half-lives) and",
      "modelled as natural logarithms."
    ),
    notes          = paste(
      "BENDITA trial, ClinicalTrials.gov NCT03378661, conducted",
      "2016-2018. 210 adults were randomised (30 per arm); the PK",
      "analysis retained 175 of the 180 actively treated subjects after",
      "excluding 5 with severe PK outliers possibly reflecting treatment",
      "misallocation. Estimation in NONMEM 7.4.3 with FOCE-I;",
      "uncertainty from a 1000-sample nonparametric bootstrap.",
      "Demographics from Assmus 2025 Table 1 (PK analysis dataset",
      "column); final parameter estimates from Assmus 2025 Table 2."
    )
  )

  ini({
    # Structural parameters. Assmus 2025 Table 2; population estimates
    # are reported "for a female adult weighting 65 kg". Because all
    # data were oral and no intravenous reference was available,
    # absolute bioavailability could not be estimated and CL and V are
    # apparent values (CL/F, V/F). Concentrations are dried blood spot
    # (whole blood) concentrations, so V/F is a blood-referenced
    # apparent volume.
    lmtt <- log(0.75);  label("Mean absorption transit time MTT (h)")            # Assmus 2025 Table 2: MTT 0.75 h (RSE 9.1%; bootstrap 95% CI 0.62-0.89)
    lcl  <- log(1.30);  label("Apparent oral clearance CL/F (L/h)")              # Assmus 2025 Table 2: CL/F 1.30 L/h (RSE 2.7%; bootstrap 95% CI 1.24-1.38)
    lvc  <- log(31.6);  label("Apparent central volume of distribution V/F (L)") # Assmus 2025 Table 2: V/F 31.6 L (RSE 2.9%; bootstrap 95% CI 29.9-33.5)

    # Relative oral bioavailability was FIXED to unity for the
    # population because no intravenous data were available (Assmus
    # 2025 Methods, Population pharmacokinetic analysis (i)). Only its
    # inter-individual variability and the sex effect were estimated.
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (fraction)")  # Assmus 2025 Table 2: "Relative oral bioavailability, F: 1 fixed"

    # Allometric scaling on body weight, applied a priori with the
    # standard theoretical exponents rather than estimated (Assmus 2025
    # Methods, Population pharmacokinetic analysis (i)).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/65) for CL/F (unitless)")  # Assmus 2025 Methods: "allometric function on all clearance (exponent 0.75) ... standardized to a body weight of 65 kg"
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on (WT/65) for V/F (unitless)")   # Assmus 2025 Methods: "... and volume (exponent 1) parameters"

    # Covariate effects retained after backward elimination (p = 0.001).
    # Both are expressed by the authors as a percentage change relative
    # to the reference subject, so they enter as fractional multipliers.
    e_conmed_fosravuconazole_cl <- 0.177;  label("Fractional increase in CL/F with concomitant fosravuconazole (unitless)")  # Assmus 2025 Table 2: "Co-administration of E1224 on CL/F (%) 17.7" (RSE 28.9%; bootstrap 95% CI 8.18-27.3)
    e_sexf_fdepot               <- -0.129; label("Fractional change in F in men (SEXF = 0) relative to the female reference (unitless)")  # Assmus 2025 Table 2: "Sex effect on F (reference: female) (%) -12.9" (RSE 23.0%; bootstrap 95% CI -18.6 to -6.9)

    # Inter-individual variability, implemented as an exponential
    # (log-normal) model. Assmus 2025 Table 2 reports IIV as %CV
    # computed from the variance as 100 * sqrt(exp(omega^2) - 1), so
    # omega^2 = log(1 + CV^2). IIV estimates below 10% were fixed to
    # zero in the final model; none was retained on V/F.
    etalmtt    ~ log(1 + 0.616^2)  # Assmus 2025 Table 2: IIV MTT 61.6% CV (RSE 15.6%; bootstrap 95% CI 35.0-75.3)
    etalcl     ~ log(1 + 0.189^2)  # Assmus 2025 Table 2: IIV CL/F 18.9% CV (RSE 16.8%; bootstrap 95% CI 12.6-25.3)
    etalfdepot ~ log(1 + 0.102^2)  # Assmus 2025 Table 2: IIV F 10.2% CV (RSE 31.9%; bootstrap 95% CI 2.97-15.5)

    # Residual unexplained variability was "an additive error on the
    # log-transformed observed concentrations (equivalent to an
    # exponential residual error on an arithmetic scale)" (Assmus 2025
    # Methods, Population pharmacokinetic analysis (i)). Table 2 reports
    # the VARIANCE, so the log-scale SD is sqrt(0.076) = 0.2757, i.e.
    # about 28% CV on the arithmetic scale.
    expSd <- sqrt(0.076); label("Residual error SD on the log-transformed concentration scale (log-normal)")  # Assmus 2025 Table 2: "Variance of residual error, sigma 0.076" (RSE 9.3%; bootstrap 95% CI 0.051-0.106)
  })

  model({
    # Individual parameters. Allometric weight scaling is centred on the
    # 65 kg cohort median; the fosravuconazole effect is a fractional
    # increase on CL/F. Assmus 2025 reports typical values for a FEMALE
    # subject, so the published male effect on F is applied through
    # (1 - SEXF) to keep the structural estimates verbatim while
    # retaining the canonical SEXF = 1 (female) orientation.
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl *
      (1 + e_conmed_fosravuconazole_cl * CONMED_FOSRAVUCONAZOLE)
    vc  <- exp(lvc) * (WT / 65)^e_wt_vc

    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_sexf_fdepot * (1 - SEXF))

    # Transit-absorption rate constant. Assmus 2025 Fig 2: ktr =
    # (1 + n) / MTT with n = 1 transit compartment, i.e. the chain is
    # gut -> transit1 -> central with two sequential transfers at rate
    # ktr, so the mean absorption transit time is 2 / ktr = MTT.
    ktr <- (1 + 1) / mtt
    kel <- cl / vc

    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1
    d/dt(central)  <-  ktr * transit1 - kel * central

    f(depot) <- fdepot

    # Dried blood spot concentration: mg / L, matching the HPLC-MS/MS
    # assay (LLOQ 50 ng/mL = 0.05 mg/L).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
