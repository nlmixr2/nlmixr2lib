Yu_2026_tenofovir_pregnancy_tfv <- function() {
  description <- paste(
    "Plasma tenofovir (TFV) two-compartment first-order-absorption popPK model in",
    "non-pregnant, pregnant and postpartum women receiving oral tenofovir disoproxil",
    "fumarate (TDF) 300 mg once daily. Non-pregnant (CONRAD 137, DOT-DBS) and",
    "pregnant / postpartum (IMPAACT P1026s TDF arm) plasma TFV data were fitted",
    "jointly, with apparent clearance CL/F estimated separately for each of four",
    "pregnancy states (non-pregnant, second trimester, third trimester, 6-12 weeks",
    "postpartum) and all disposition parameters shared across states. Tenofovir",
    "disoproxil fumarate itself is not carried as a state: its plasma half-life is",
    "about 24 seconds, so a 300 mg TDF dose is given directly into the depot as its",
    "136 mg TFV-equivalent (473 umol). Doses and states are in umol and",
    "concentrations in umol/L. This is the covariate-identification sub-model that",
    "supplied the pregnancy effect on plasma TFV clearance used by the",
    "semi-mechanistic model Yu_2026_tenofovir; note that the absolute clearances",
    "here imply smaller trimester effects than the percentages the paper carried",
    "into its simulations (see the vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Yu Y, Brooks KM, Doncel GF, Best BM, Marzinke MA, Mirochnick M, Anderson P,",
    "Myer L, Celum C, Heffron R, Coleman J, Joseph Davey D, Hendrix CW, Momper JD,",
    "Bies R, Scott RK. Development of a Semi-Mechanistic Population Pharmacokinetic",
    "Model for Predicting Tenofovir Disoproxil Fumarate and Tenofovir Alafenamide",
    "Exposure in Plasma and Cellular Matrices During Pregnancy and Postpartum.",
    "Clin Pharmacokinet. 2026;65(1):133-148. doi:10.1007/s40262-025-01589-y.",
    "Parameters from Electronic Supplementary Material Table S1.",
    sep = " "
  )
  vignette <- "Yu_2026_tenofovir"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    EGA = list(
      description        = "Estimated gestational age of the mother",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's second- and third-trimester indicator",
        "variables TRI2 and TRI3; the model has no continuous EGA effect. Per the",
        "EGA register entry, a source paper that reports trimesters rather than",
        "weeks records the week values used per trimester here rather than",
        "introducing a trimester-indicator canonical. This model uses the standard",
        "obstetric boundaries TRI2 = 14 <= EGA < 28 weeks and TRI3 = EGA >= 28",
        "weeks. Yu 2026 does not print its own trimester cut-offs; the contributing",
        "IMPAACT P1026s protocol sampled the second trimester at 20-26 weeks and",
        "the third trimester at 30-38 weeks, both of which fall unambiguously",
        "inside these boundaries. Set EGA = 0 for non-pregnant and for postpartum",
        "subjects (the postpartum state is carried by TPP, not by EGA).",
        "No first-trimester data were available, so the model makes no prediction",
        "for 0 < EGA < 14 weeks; such a record is scored as non-pregnant.",
        sep = " "
      ),
      source_name        = "TRI2 / TRI3 (trimester indicator variables)"
    ),
    TPP = list(
      description        = "Time postpartum (time elapsed since delivery)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's postpartum indicator PP; the model has no",
        "continuous TPP effect and applies a single step change for any TPP > 0.",
        "Yu 2026 sampled the postpartum state at 6-12 weeks after delivery, so the",
        "estimated postpartum clearance describes that window and should not be",
        "read as an immediately-post-delivery value. Set TPP = 0 during pregnancy",
        "and for non-pregnant subjects.",
        sep = " "
      ),
      source_name        = "PP (postpartum indicator variable)"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "tenofovir", units = "umol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 46L,
    n_studies     = 3L,
    disease_state = paste(
      "Pregnant and postpartum women with HIV taking TDF 300 mg once daily",
      "(IMPAACT P1026s TDF arm, 46 participants sampled in the second trimester,",
      "the third trimester and 6-12 weeks postpartum), jointly fitted with plasma",
      "TFV data from healthy HIV-negative non-pregnant women in CONRAD 137 (TDF",
      "arm) and the DOT-DBS study.",
      sep = " "
    ),
    dose_range = "TDF 300 mg orally once daily (modelled as 136 mg = 473 umol TFV)",
    regions    = "USA and international IMPAACT network sites; CONRAD 137 and DOT-DBS conducted in the USA",
    notes      = paste(
      "n_subjects counts the 46 pregnant IMPAACT P1026s TDF-arm participants",
      "(Table 1). The number of non-pregnant participants contributing plasma TFV",
      "to this particular fit is not stated separately; Table 1 reports 24 and 73",
      "non-pregnant women in the CONRAD 137 single- and multiple-dose phases",
      "(which may overlap) and 28 in DOT-DBS, across TDF and TAF arms combined.",
      "Age, weight, sex and race distributions are not reported in Yu 2026; the",
      "cohort is entirely female by design. Populations were deliberately not",
      "matched on age or race because covariate data were limited in the pooled",
      "datasets (Methods 2.1).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Apparent clearance, estimated separately in each pregnancy state.
    # ESM Table S1 prints one absolute CL/F per state rather than a
    # reference value plus offsets, so each state carries an explicit
    # stratum suffix (parameter-names.md, "Stratum-suffixed parameters")
    # and none keeps the bare `lcl`. The paper's covariate equation
    # (Methods 2.4.4) is additive, TVP = th1 + TRI2*th2 + TRI3*th3 +
    # PP*th4, which is exactly a per-state absolute value.
    # ------------------------------------------------------------------
    lcl_nonpreg <- log(55.1); label("Apparent clearance CL/F of plasma TFV, non-pregnant (L/h)")            # ESM Table S1 (CL/F Nonpregnant = 55.1 L/h, RSE 5%)
    lcl_tri2    <- log(65.3); label("Apparent clearance CL/F of plasma TFV, second trimester (L/h)")        # ESM Table S1 (CL/F 2nd Trimester = 65.3 L/h, RSE 18%)
    lcl_tri3    <- log(59.4); label("Apparent clearance CL/F of plasma TFV, third trimester (L/h)")         # ESM Table S1 (CL/F 3rd Trimester = 59.4 L/h, RSE 5%)
    lcl_pp      <- log(47.6); label("Apparent clearance CL/F of plasma TFV, 6-12 weeks postpartum (L/h)")   # ESM Table S1 (CL/F Postpartum = 47.6 L/h, RSE 5%)

    # Disposition parameters shared across all four pregnancy states.
    lka <- log(1.02); label("First-order absorption rate constant of TFV after oral TDF (1/h)")             # ESM Table S1 (Ka = 1.02 /h, RSE 25%)
    lvc <- log(427);  label("Apparent central volume of distribution V2/F of plasma TFV (L)")               # ESM Table S1 (V2/F = 427 L, RSE 13%)
    lvp <- log(711);  label("Apparent peripheral volume of distribution V3/F of plasma TFV (L)")            # ESM Table S1 (V3/F = 711 L, RSE 6%)
    lq  <- log(135);  label("Apparent intercompartmental clearance Q/F of plasma TFV (L/h)")                # ESM Table S1 (Q/F = 135 L/h, RSE 16%)

    # ------------------------------------------------------------------
    # Between-subject variability. ESM Table S1 reports BSV as %CV for
    # the exponential model P = TVP * exp(eta); omega^2 = log(CV^2 + 1).
    # The single CL/F BSV of 29% is printed against all four CL/F rows,
    # i.e. one shared eta on clearance rather than one per state.
    # ------------------------------------------------------------------
    etalcl ~ 0.0808   # CV 29%:    log(0.290^2 + 1)   (ESM Table S1, BSV on CL/F, RSE 8%)
    etalka ~ 1.0433   # CV 135.6%: log(1.356^2 + 1)   (ESM Table S1, BSV on Ka, RSE 14%)
    etalvc ~ 0.2088   # CV 48.2%:  log(0.482^2 + 1)   (ESM Table S1, BSV on V2/F, RSE 23%)

    # Residual error, plasma TFV.
    propSd <- 0.327;  label("Proportional residual error, plasma TFV (fraction)")                           # ESM Table S1 (PROP (TFV) = 32.7%, RSE 2%)
    addSd  <- 0.0269; label("Additive residual error, plasma TFV (umol/L)")                                 # ESM Table S1 (ADD (TFV) = 0.0269 umol/L, RSE 1%)
  })

  model({
    # 1. Pregnancy state. The paper's TRI2 / TRI3 / PP indicators are
    #    derived here from the canonical EGA and TPP covariate columns.
    #    Order matters: postpartum is tested first so that a record with
    #    TPP > 0 is scored postpartum regardless of any residual EGA.
    lclState <-
      ifelse(TPP > 0, lcl_pp,
        ifelse(EGA >= 28, lcl_tri3,
          ifelse(EGA >= 14, lcl_tri2, lcl_nonpreg)))

    # 2. Individual parameters
    cl <- exp(lclState + etalcl)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system -- two-compartment with first-order oral absorption.
    #    The depot receives the TFV-equivalent of the TDF dose.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation and error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
