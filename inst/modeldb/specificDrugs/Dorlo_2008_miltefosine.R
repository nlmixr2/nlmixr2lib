Dorlo_2008_miltefosine <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption",
    "and linear elimination for miltefosine in 31 Dutch military personnel",
    "with Old World (Leishmania major) cutaneous leishmaniasis contracted",
    "in Afghanistan (Dorlo 2008), treated with oral miltefosine 50 mg three",
    "times daily (150 mg/day, median 1.76 mg/kg/day) for 28 days with",
    "post-treatment follow-up to a maximum of 202 days. CL/F, Vc/F, Q/F,",
    "and Vp/F are estimated apparent parameters; relative bioavailability F",
    "is unidentifiable from oral-only data and is structurally fixed at 1.",
    "Inter-individual variability is log-normal on ka, CL/F, and Vc/F",
    "(diagonal in this implementation; see Assumptions in the vignette for",
    "the unreported CL/Vc correlation noted by the authors). IIV on Q/F and",
    "Vp/F was not estimable from the data. Residual error is proportional",
    "(31.5% CV). No covariate effects were retained in the final model.",
    "This is the structural model later re-used as the base PK structure",
    "in Dorlo 2017 visceral-leishmaniasis miltefosine work."
  )
  reference <- paste(
    "Dorlo TPC, van Thiel PPAM, Huitema ADR, Keizer RJ, de Vries HJC,",
    "Beijnen JH, de Vries PJ. Pharmacokinetics of miltefosine in Old",
    "World cutaneous leishmaniasis patients. Antimicrob Agents Chemother.",
    "2008;52(8):2855-2860. doi:10.1128/AAC.00014-08.",
    sep = " "
  )
  vignette <- "Dorlo_2008_miltefosine"
  units    <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported in Dorlo 2008 Table 1 (median 24 years, IQR 23-29) but",
        "no age effect was retained in the final pharmacokinetic model.",
        "The cohort was a narrow young-adult military population so the",
        "study had limited power to detect an age effect."
      ),
      source_name        = "Age (yr)"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported in Dorlo 2008 Table 1 (median 85 kg, IQR 78-89). No",
        "weight effect was retained in the final pharmacokinetic model.",
        "The cohort spans a narrow adult weight range so the study had",
        "limited power to detect an allometric effect; the dose was a",
        "flat 50 mg three times daily independent of weight (median",
        "1.76 mg/kg/day, IQR 1.69-1.92)."
      ),
      source_name        = "Weight (kg)"
    ),
    HEIGHT = list(
      description        = "Body height at baseline",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported in Dorlo 2008 Table 1 (median 184 cm, IQR 180-188).",
        "Not retained as a covariate in the final pharmacokinetic model."
      ),
      source_name        = "Height (cm)"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Reported in Dorlo 2008 Table 1 as 'No. of males/no. of females:",
        "30/1'. Not assessed as a covariate in the final model: with only",
        "one female subject the data carried no information for a sex",
        "effect, and the cohort is effectively male-only."
      ),
      source_name        = "No. of males/no. of females"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    age_range      = "median 24 years (IQR 23-29)",
    age_median     = "24 years",
    weight_range   = "median 85 kg (IQR 78-89)",
    weight_median  = "85 kg",
    height_range   = "median 184 cm (IQR 180-188)",
    height_median  = "184 cm",
    sex_female_pct = 100 / 31,
    race_ethnicity = "Dutch military personnel and embedded civilians (no further breakdown reported)",
    disease_state  = paste(
      "Old World cutaneous leishmaniasis (Leishmania major), parasite",
      "infection contracted during deployment in northern Afghanistan and",
      "confirmed by microscopy with PCR/NASBA genotyping in 27 of 31",
      "patients. 30 of 31 patients had received prior intralesional",
      "pentavalent antimony (SbV) before initiating miltefosine; one",
      "patient was treatment-naive. Patients were treated with miltefosine",
      "as second-line therapy for insufficient response to intralesional",
      "SbV or as primary systemic therapy for extensive disease."
    ),
    dose_range     = paste(
      "Oral miltefosine (Impavido, Zentaris GmbH) 50 mg three times daily",
      "for 28 days (total 150 mg/day; median 1.76 mg/kg/day, IQR 1.69-",
      "1.92). Administered with a meal or snack to mitigate gastrointestinal",
      "side effects."
    ),
    regions        = "The Netherlands (Academic Medical Center, Amsterdam); patients were Dutch ISAF military personnel deployed to Mazar-e-Sharif, Afghanistan",
    co_medication  = "Most patients had received prior intralesional pentavalent antimony (SbV); none received concurrent systemic anti-leishmaniasis treatment during the miltefosine course.",
    samples        = paste(
      "382 miltefosine plasma concentrations from 31 patients (median 13",
      "samples per patient, range 9-20; 8 during treatment and 5 after",
      "treatment per patient at the median). Sampling on day 1 captured",
      "2, 4, and 6 h post-first-dose; further during-treatment samples",
      "were collected on an outpatient basis; post-treatment sampling",
      "continued up to day 202 (maximum follow-up). All samples assayed",
      "by validated LC-MS/MS (LLOQ 4 ng/mL); all post-baseline samples",
      "were above LLOQ."
    ),
    notes          = paste(
      "Demographics from Dorlo 2008 Table 1. 3 additional patients were",
      "screened but excluded from the PK dataset due to inconsistent",
      "labelling at the sampling site (final analysis n = 31). NONMEM",
      "VI was used with FOCE-INTERACTION; visual predictive check based",
      "on 2,000 simulated subjects (Dorlo 2008 Figure 1). All 30/31",
      "patients except one had prior intralesional antimony; none of",
      "these baseline characteristics was retained as a covariate in the",
      "final pharmacokinetic model (see covariatesDataExcluded)."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Dorlo 2008 Table 2 'Final parameter
    # estimates for pharmacokinetic model for miltefosine'. NONMEM
    # names V2/V3 map to nlmixr2 conventions Vc/Vp (linCmt() macro
    # V1 = depot, V2 = central, V3 = peripheral). F is unidentifiable
    # from oral-only data and is structurally fixed at 1 (the paper
    # reports parameters relative to F: CL/F, V2/F, Q/F, V3/F). The
    # absorption rate is reported in 1/h in the paper; here it is
    # rescaled to 1/day for consistency with CL/F and Q/F (also in
    # day units) so the model time axis is uniformly in days.
    # ============================================================
    lka <- log(0.36 * 24)
    label("First-order oral absorption rate constant ka (1/day)")          # Table 2: ka = 0.36 /h (RSE 10.1%) -> 0.36 x 24 = 8.64 /day
    lcl <- log(3.87)
    label("Apparent oral clearance CL/F (L/day)")                          # Table 2: CL/F = 3.87 L/day (RSE 5.3%)
    lvc <- log(39.6)
    label("Apparent central volume Vc/F (L)")                              # Table 2: V2/F = 39.6 L (RSE 4.0%); V2 in NONMEM notation is the central compartment for two-compartment oral models
    lq  <- log(0.0375)
    label("Apparent inter-compartmental clearance Q/F (L/day)")            # Table 2: Q/F = 0.0375 L/day (RSE 22.0%)
    lvp <- log(1.65)
    label("Apparent peripheral volume Vp/F (L)")                           # Table 2: V3/F = 1.65 L (RSE 12.4%)

    # Bioavailability anchor (F is unidentifiable from oral-only data,
    # so it is structurally fixed at 1; the structural parameters above
    # are all apparent CL/F, V/F, Q/F).
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless, FIXED at 1)")              # Methods 'Pharmacokinetic data analysis' paragraph 3: 'Bioavailability (F) was unknown, and therefore, parameters relative to the bioavailability were estimated (CL/F, V/F, etc.)'

    # ============================================================
    # Inter-individual variability -- Dorlo 2008 Table 2 column
    # '% Interindividual variability (relative SE [%])'. The paper
    # uses an exponential (log-normal) IIV model and reports CV%; the
    # internal log-normal variance is omega^2 = log(1 + CV^2).
    # IIV on Q/F and Vp/F is reported as 'NE' (not estimable) and is
    # therefore not declared here. The 'a' footnote in Table 2 marks
    # CL/F and V2/F IIVs as highly correlated (Results 'Pharmacokinetic
    # data analysis' paragraph 2: 'The high correlation between the
    # interindividual variability for CL and V_2 can probably be
    # attributed to variability in bioavailability or in the unbound
    # drug fraction'); the paper does NOT report the correlation
    # value, so the IIV is encoded here as diagonal -- see the
    # vignette's Assumptions and deviations section.
    # ============================================================
    etalka ~ log(1 + 0.242^2)
    # Table 2 IIV row 'Absorption rate (k_a) (h^-1)': 24.2% CV (RSE 63.3%) -> omega^2 = log(1 + 0.242^2)
    etalcl ~ log(1 + 0.232^2)
    # Table 2 IIV row 'Clearance (CL/F) (liters/day)': 23.2% CV (RSE 15.4%) -> omega^2 = log(1 + 0.232^2). Correlated with etalvc per footnote a; correlation value not published.
    etalvc ~ log(1 + 0.183^2)
    # Table 2 IIV row 'Volume of central compartment (V_2/F) (liters)': 18.3% CV (RSE 25.0%) -> omega^2 = log(1 + 0.183^2). Correlated with etalcl per footnote a; correlation value not published.

    # ============================================================
    # Residual error -- Dorlo 2008 Table 2: proportional residual
    # error 31.5% (Methods 'Pharmacokinetic data analysis': 'Residual
    # variability was modeled with a proportional error model').
    # ============================================================
    propSd <- 0.315
    label("Proportional residual error on Cc (fraction)")                  # Table 2: Residual variability 31.5% (RSE 6.4%)
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters. No covariates retained in the final
    # model (Dorlo 2008 does not report any covariate effects; the
    # cohort is a narrow young-adult military population).
    # ------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------
    # Two-compartment oral PK (depot -> central <-> peripheral1)
    # with first-order absorption and linear elimination from the
    # central compartment (Dorlo 2008 Results 'Pharmacokinetic data
    # analysis' paragraph 3: 'An open two-compartmental model with
    # first-order absorption and linear elimination from the central
    # compartment best fitted the data').
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------
    # Observation. Dose units mg; Vc units L -> central/vc in mg/L
    # which equals ug/mL (the units declared in `units` above). The
    # paper reports concentrations in ng/mL; the vignette converts
    # for the side-by-side comparison.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
