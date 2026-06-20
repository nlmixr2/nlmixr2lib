Agoram_2006_darbepoetin_alfa <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order subcutaneous",
    "absorption for darbepoetin alfa in healthy adult subjects (Agoram 2006).",
    "Both IV and SC routes are supported. SC bioavailability is a linear",
    "function of the SC dose amount (in ug). Body weight modifies clearance",
    "and central volume via a normalized power model (reference 70 kg);",
    "subject age modifies the absorption rate constant via a normalized power",
    "model (reference 47 years, the development-cohort mean). Total measured",
    "serum concentration is the sum of the simulated darbepoetin alfa and an",
    "individual-specific endogenous-erythropoietin (eEPO) constant that the",
    "ELISA assay cross-detects. Exponential (log-normal) residual error."
  )
  reference <- paste(
    "Agoram B, Sutjandra L, Sullivan JT. Population pharmacokinetics of",
    "darbepoetin alfa in healthy subjects.",
    "British Journal of Clinical Pharmacology. 2006;63(1):41-52.",
    "doi:10.1111/j.1365-2125.2006.02752.x"
  )
  vignette <- "Agoram_2006_darbepoetin_alfa"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on clearance and central volume; reference",
        "weight 70 kg per Agoram 2006 Discussion ('For an average 70-kg human,",
        "the estimated mean s.c. relative bioavailability ...'). Development",
        "cohort mean 68.6 +/- 10.3 kg (Agoram 2006 Table 2)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on the first-order SC absorption rate Ka.",
        "Reference age 47 years is the development-cohort mean (Agoram 2006",
        "Table 2: 47 +/- 17 years); the paper does not explicitly state the",
        "reference age, so the cohort mean is used as a defensible default,",
        "consistent with the unstated-reference precedent in",
        "Bi_2017_peginterferon_alfa_2a and Cirincione_2017_exenatide.",
        "Choosing the cohort mean as reference makes the reported typical Ka",
        "= 0.0212 1/h equal the typical-subject value (and matches the",
        "paper's reported absorption half-life ln(2)/Ka = 33 h, derivable",
        "from Ka = 0.0212 only). See vignette Assumptions and deviations."
      ),
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 140L,
    n_studies      = 6L,
    n_observations = 1664L,
    age_range      = paste(
      "Development cohort 47 +/- 17 years (mean +/- SD); evaluation cohort",
      "55 +/- 18 years (Agoram 2006 Table 2)."
    ),
    age_median     = "Mean 47 years in the development cohort (median not separately tabulated)",
    weight_range   = paste(
      "Development cohort 68.6 +/- 10.3 kg; evaluation cohort 71.0 +/- 10.8",
      "kg (Agoram 2006 Table 2)."
    ),
    weight_median  = "Mean 68.6 kg in the development cohort (median not separately tabulated)",
    sex_female_pct = 55.7,
    race_ethnicity = "Not explicitly tabulated in Agoram 2006.",
    disease_state  = paste(
      "Healthy adult volunteers (age >= 18 years; normal physical exam and",
      "12-lead ECG; transferrin saturation >= 15%; normal serum vitamin B12",
      "and folate; screening haemoglobin <= 15.0 g/dL; no infection with HIV,",
      "hepatitis B or C virus; no clinically significant cardiovascular,",
      "hepatic, or renal impairment; no primary haematological disorder; no",
      "recent ESA / blood donation / transfusion exposure; not pregnant)."
    ),
    dose_range     = paste(
      "IV 0.75 ug/kg single dose; SC 0.75-8.0 ug/kg single or multiple doses",
      "(0.75, 2.0, 3.0, 5.0, 6.5, 8.0 ug/kg single; 1.0 ug/kg Q6W x 2; 3.0",
      "ug/kg Q3W x 2; 6.5 ug/kg Q3W x 2); SC 80 ug Q4W x 2; SC 2.0 ug/kg QW",
      "x 4; SC 500 ug Q3W x 2 (Agoram 2006 Table 1)."
    ),
    regions        = paste(
      "Six Amgen-sponsored clinical studies (study numbers 20010262, 990134,",
      "20010198, 20010174, 20030163, 20000250); geographic sites not stated",
      "in the publication."
    ),
    notes          = paste(
      "Total N = 140 healthy subjects randomly split 50:50 into model",
      "development (N = 70, 1664 plasma samples) and model evaluation",
      "(N = 70) sets. Serum darbepoetin alfa quantified by the validated",
      "Quantikine IVD human erythropoietin ELISA (R&D Systems, Minneapolis",
      "MN); standard curve 0.125-5.0 ng/mL, LLOQ 0.14 ng/mL. Endogenous EPO",
      "cross-reacted with the assay and was modelled as an individual-",
      "specific additive constant (typical value 0.0867 ng/mL). NONMEM V",
      "FOCE-I (PREDPP IV) was used for estimation. Demographics from",
      "Agoram 2006 Table 2; trial summary from Table 1."
    )
  )

  ini({
    # Structural PK parameters (Agoram 2006 Table 3 'Final population covariate
    # model parameter estimates'). TVP values evaluated at the reference
    # covariates: WT = 70 kg, AGE = 47 years.
    lcl   <- log(0.164);  label("Typical clearance at WT 70 kg (L/h)")                                 # Agoram 2006 Table 3: CL = 0.164 L/h (%SEP 11.6)
    lvc   <- log(5.98);   label("Typical central volume at WT 70 kg (L)")                              # Agoram 2006 Table 3: V1 = 5.98 L (%SEP 13.7)
    lvp   <- log(1.21);   label("Peripheral volume (L)")                                               # Agoram 2006 Table 3: V2 = 1.21 L (%SEP 40.4)
    lq    <- log(0.0153); label("Inter-compartmental clearance (L/h)")                                 # Agoram 2006 Table 3: Q  = 0.0153 L/h (%SEP 38.1)
    lka   <- log(0.0212); label("Typical first-order SC absorption rate Ka at AGE 47 years (1/h)")     # Agoram 2006 Table 3: Ka = 0.0212 1/h (%SEP 3.54)
    leepo <- log(0.0867); label("Typical endogenous-EPO contribution to serum concentration (ng/mL)")  # Agoram 2006 Table 3: eEPO = 0.0867 ng/mL (%SEP 4.19)

    # SC bioavailability as a linear function of SC dose amount (Agoram 2006
    # Results equation 7: F = F0 + p1 * Dose; Dose in ug). IIV on F0 and p1
    # was fixed at 0 in the final model (Table 3 footnote: 'Interindividual
    # random variance was fixed at 0 in the PK model').
    f0depot <- 0.448;    label("Intercept of SC bioavailability vs. dose (F0, fraction)")               # Agoram 2006 Table 3: F0 = 0.448 (%SEP 10.2; IIV fixed at 0)
    p1depot <- 0.000586; label("Linear slope of SC bioavailability vs. SC dose amount (per ug)")        # Agoram 2006 Table 3: p1 = 5.86e-4 per ug (%SEP 37.2; IIV fixed at 0)

    # Covariate effects (normalized power model; Agoram 2006 Results equations
    # 8-10 with reference values WT 70 kg, AGE 47 years).
    e_wt_cl  <-  1.19;   label("Power exponent of (WT/70 kg) on clearance (unitless)")                 # Agoram 2006 Table 3: r1 (BWT on CL) = 1.19 (%SEP 21.9)
    e_wt_vc  <-  0.983;  label("Power exponent of (WT/70 kg) on central volume (unitless)")            # Agoram 2006 Table 3: r2 (BWT on V1) = 0.983 (%SEP 44.9)
    e_age_ka <- -0.951;  label("Power exponent of (AGE/47 years) on Ka (unitless)")                    # Agoram 2006 Table 3: r3 (Age on Ka) = -0.951 (%SEP 7.48)

    # IIV - exponential / log-normal, diagonal Omega. NONMEM omega^2 reported
    # directly in Table 3 as variances on log-transformed parameters. V2, Q,
    # F0, p1 had IIV fixed at 0 in the source (Table 3 footnote 'Inter-
    # individual random variance was fixed at 0 in the PK model'); they are
    # therefore not given etas here.
    etalcl   ~ 0.075    # Agoram 2006 Table 3: omega^2 CL    = 0.075   (~CV 27.4%)
    etalvc   ~ 0.227    # Agoram 2006 Table 3: omega^2 V1    = 0.227   (~CV 47.6%)
    etalka   ~ 0.0832   # Agoram 2006 Table 3: omega^2 Ka    = 0.0832  (~CV 28.8%)
    etaleepo ~ 0.132    # Agoram 2006 Table 3: omega^2 eEPO  = 0.132   (~CV 36.3%)

    # Residual error - exponential (NONMEM Y = F * exp(eps), eps ~ N(0,
    # sigma^2)); encoded here as Cc ~ lnorm(expSd) with SD on the log scale.
    # sigma^2 = 0.324 -> SD = sqrt(0.324) = 0.5692.
    expSd <- sqrt(0.324); label("Log-scale residual SD (~CV 56.9%)")                                   # Agoram 2006 Table 3: sigma1^2 = 0.324 (CV ~56.9%)
  })

  model({
    # Individual PK parameters with normalized power covariate effects.
    cl   <- exp(lcl + etalcl) * (WT  / 70)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT  / 70)^e_wt_vc
    vp   <- exp(lvp)
    q    <- exp(lq)
    ka   <- exp(lka + etalka) * (AGE / 47)^e_age_ka
    eepo <- exp(leepo + etaleepo)

    # Micro-constants for the two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV doses enter central directly (dataset cmt = central; default
    # f(central) = 1). SC doses enter depot (dataset cmt = depot) with the
    # dose-dependent bioavailability defined below.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # SC bioavailability is linear in the SC dose amount (Agoram 2006
    # Results equation 7). podo(depot) returns the current dose in dose
    # units (ug) at the moment the dose is administered.
    f(depot) <- f0depot + p1depot * podo(depot)

    # Total measured serum concentration = simulated darbepoetin alfa +
    # individual eEPO. Dose in ug, central in ug, vc in L -> central/vc in
    # ug/L = ng/mL; eepo also in ng/mL.
    Cc <- central / vc + eepo
    Cc ~ lnorm(expSd)
  })
}
