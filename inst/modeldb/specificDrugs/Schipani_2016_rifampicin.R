Schipani_2016_rifampicin <- function() {
  description <- "Simultaneous population pharmacokinetic model for oral rifampicin in a mixed Malawian cohort of adults (n=115) and children (n=50) with tuberculosis. One-compartment disposition with first-order absorption (depot to central). Allometric scaling of CL/F and V/F to a 70 kg reference body weight with canonical Anderson and Holford (2008) exponents (0.75 on CL, 1.0 on V; both fixed). Estimated power-form effect of age on CL/F (exponent 0.517) centered at AGE_median. Children (defined as body weight 5-29 kg, age < 15 y per Schipani 2016 Results) carry a relative bioavailability factor F = 0.517 vs adults (F fixed at 1). Inter-individual variability is carried on CL/F (46.6% approx CV) and V/F (87.4% approx CV); ka and F have no estimated IIV. Proportional residual error 48%."
  reference <- "Schipani A, Pertinez H, Mlota R, Molyneux E, Lopez N, Dzinjalamala FK, van Oosterhout JJ, Ward SA, Khoo S, Davies G. (2016). A simultaneous population pharmacokinetic analysis of rifampicin in Malawian adults and children. British Journal of Clinical Pharmacology 81(4):679-687. doi:10.1111/bcp.12848"
  vignette <- "Schipani_2016_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the Schipani 2016 cohort. Used for allometric scaling on CL/F (exponent 0.75 fixed) and V/F (exponent 1.0 fixed) with reference weight 70 kg per Schipani 2016 Methods 'For weight as a covariate an allometric model was applied to standardize the CL and V pharmacokinetic parameters using a standard weight (WTstd) of 70 kg ... and fixing the exponent to 0.75 for CL and 1 for V.' Cohort range 4.8 - 87 kg (Table 2).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power-form normalized covariate effect on CL/F per Schipani 2016 equation 4: CL_i = (theta_CL * (WT_i/70)^0.75 * (AGE_i/AGE_median)^theta_age) * exp(eta_CL). The paper does not report AGE_median explicitly; we use 33 years (the adult median per Table 2) so that the equation evaluated at the adult typical demographics (WT = 70 kg, AGE = 33 y, F = 1) gives the paper's stated typical CL/F = 23.9 L/h ('mean population estimate for CL/F in the adult subpopulation (F fixed to 1)'). Cohort age range 0.58 (7 months) - 65 years (Table 2). See vignette Errata for the centering assumption.",
      source_name        = "AGE"
    ),
    CHILD = list(
      description        = "Binary indicator for child subpopulation (1 = child, 0 = adult)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Schipani 2016 Results: 'Individuals with a weight between 5 and 29 kg were considered children, which corresponded with age < 15 years.' Used to apply the relative bioavailability factor F = 0.517^CHILD on the depot. Adults: CHILD = 0 implies F = 1; children: CHILD = 1 implies F = 0.517 (children have 48.3% lower bioavailability than adults per Discussion). The paper used the weight cutoff < 30 kg; in the implementation either WT < 30 kg or AGE < 15 y can derive CHILD.",
      source_name        = "CHILD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 165L,
    n_studies      = 1L,
    age_range      = "0.58 - 65 years (overall); children 0.58 - 14 y (median 6.125); adults 14 - 65 y (median 33)",
    weight_range   = "4.8 - 87 kg (overall); children 4.8 - 29 kg (median 15); adults 30 - 87 kg (median 49)",
    sex_female_pct = round(68 / 165 * 100, 1),
    race_ethnicity = "Malawian (not stratified further)",
    disease_state  = "Tuberculosis patients receiving rifampicin in fixed-dose-combination anti-TB therapy. HIV co-infection: 62% of children, 70% of adults; HIV status not a significant covariate.",
    dose_range     = "Oral rifampicin in fixed-dose-combination tablets per Malawian weight-banded guidelines (Schipani 2016 Table 1): children 60-300 mg/day (RHZ R60/H30/Z150 or RH R60/H60 tablets, 1-5 tablets per weight band 0-29 kg); adults 300-750 mg/day (RHZE R150/H75/Z400/E275 tablets, 2-5 tablets per weight band 30 to >=75 kg).",
    regions        = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    sampling_design = "Mixed rich + sparse design. Rich PK: 40 adults and 22 children (5-6 samples per patient at 0 pre-dose, 0.5/1/2/3/4/6/8 h, and 24 h post-dose). Sparse PK: 75 adults and 28 children (1-2 samples per patient in a window up to 8 h post-dose). Total 608 plasma concentrations.",
    assay          = "Plasma rifampicin by validated HPLC; LLOQ 0.5 mg/L.",
    notes          = "Mixed adult+pediatric simultaneous popPK fit. Covariates evaluated by stepwise forward-backward selection: weight, age, gender, HIV status. Only weight (allometric, fixed) and age (estimated power form) retained on CL/F; a fixed child-vs-adult relative bioavailability factor was retained on F. Patients enrolled at least 2 weeks after starting intensive-phase TB treatment."
  )

  ini({
    # ============================================================
    # Structural PK at 70 kg reference body weight and AGE_median
    # (Schipani 2016 Table 4 'Final parameter estimates'). The
    # paper's equation 4 expresses CL_i in terms of theta_CL with
    # the per-individual (WT/70)^0.75 and (AGE/AGE_median)^0.517
    # multipliers. theta_CL = 23.9 L/h is described as the typical
    # adult CL/F (F fixed to 1), so AGE_median is taken as the
    # adult median age 33 y (see covariateData[[AGE]]$notes and
    # vignette Errata).
    # ============================================================
    lcl <- log(23.9)
    label("Apparent oral clearance CL/F at 70 kg, AGE_median 33 y (L/h)")
    # Schipani 2016 Table 4 'CL/F (l h-1) = 23.9 (RSE 6%)'.

    lvc <- log(44.6)
    label("Apparent central volume of distribution V/F at 70 kg (L)")
    # Schipani 2016 Table 4 'V/F (l) = 44.6 (RSE 11%)'.

    lka <- log(0.236)
    label("First-order absorption rate constant ka (1/h)")
    # Schipani 2016 Table 4 'ka (h-1) = 0.236 (RSE 3%)'. The paper
    # reports no IIV on ka ('IIV was not estimated for absorption
    # constant (ka), probably due to the limited data').

    lfdepot <- fixed(log(1))
    label("Bioavailability anchor (adult; fixed at unity)")
    # Schipani 2016 Results: 'for the adults this bioavailability
    # factor was fixed at unity.' The reported child relative-F
    # value (0.517) is encoded as a multiplicative factor below.

    # ============================================================
    # Allometric exponents fixed at canonical values (Schipani 2016
    # Methods: 'fixing the exponent to 0.75 for CL and 1 for V').
    # ============================================================
    allo_cl <- fixed(0.75)
    label("Allometric exponent on CL (unitless; fixed)")
    # Schipani 2016 Methods, citing reference 18 (Anderson and
    # Holford 2008): 0.75 for CL.

    allo_v <- fixed(1.0)
    label("Allometric exponent on V (unitless; fixed)")
    # Schipani 2016 Methods, citing reference 18 (Anderson and
    # Holford 2008): 1.0 for V.

    # ============================================================
    # Estimated covariate effects (Schipani 2016 Table 4).
    # ============================================================
    e_age_cl <- 0.517
    label("Power exponent of age on CL/F (unitless)")
    # Schipani 2016 Table 4 'Factor associated with age on RIF CL/F
    # = 0.517 (RSE 18%)'. Applied per equation 4 as
    # (AGE/AGE_median)^0.517 with AGE_median = 33 y (see ini-block
    # header note).

    e_child_fdepot <- 0.517
    label("Child relative bioavailability ratio (F_child / F_adult)")
    # Schipani 2016 Table 4 'F (% relative bioavailability,
    # children) = 51.7 (RSE 18%)'. Encoded as 0.517 (fractional)
    # and applied as F_i = 0.517^CHILD; CHILD = 0 (adult) gives F =
    # 1, CHILD = 1 (child) gives F = 0.517. The Discussion's '48.3%'
    # is the complementary 1 - 0.517 reduction.

    # ============================================================
    # Inter-individual variability (Schipani 2016 Table 4).
    # The paper reports IIV as approximate CV% (= 100 * sqrt(omega^2)
    # rather than the rigorous 100 * sqrt(exp(omega^2) - 1)): the
    # Table 3 final-row CL/F IIV omega^2 = 0.215 and Table 4 IIV
    # CL/F (%) = 46.6 are self-consistent only if 46.6% is treated
    # as 100 * sqrt(0.215) ~= 46.4% (the approximate convention).
    # The same convention back-converts Table 4 IIV V/F = 87.4% to
    # omega^2 = 0.874^2 ~= 0.764, which is the value used here
    # (Table 3 final-row 0.81 appears to be a rounded display of
    # the same estimate). No IIV on ka per Results body text.
    # ============================================================
    etalcl ~ 0.217156
    # Schipani 2016 Table 4 IIV CL/F (%) = 46.6 ->
    # omega^2 = 0.466^2 = 0.217156 (approximate CV convention).

    etalvc ~ 0.763876
    # Schipani 2016 Table 4 IIV V/F (%) = 87.4 ->
    # omega^2 = 0.874^2 = 0.763876 (approximate CV convention).

    # ============================================================
    # Residual error (Schipani 2016 Table 4).
    # ============================================================
    propSd <- 0.48
    label("Proportional residual error SD on rifampicin Cc (fraction)")
    # Schipani 2016 Table 4 'Residual error - Proportional (%) = 48
    # (RSE 8%)'. Table 3 final-row epsilon_ccv = 0.223 is the
    # sigma^2 form (sqrt(0.223) ~= 0.472 ~= 48%).
  })

  model({
    # ------------------------------------------------------------
    # 1. Body-weight allometric factors (Schipani 2016 Methods,
    #    Anderson and Holford 2008); reference 70 kg.
    # ------------------------------------------------------------
    bw_cl <- (WT / 70) ^ allo_cl
    bw_v  <- (WT / 70) ^ allo_v

    # ------------------------------------------------------------
    # 2. Age power-form factor on CL/F (Schipani 2016 equation 4);
    #    centered at AGE_median = 33 y (see ini-block header note).
    # ------------------------------------------------------------
    age_cl <- (AGE / 33) ^ e_age_cl

    # ------------------------------------------------------------
    # 3. Individual PK parameters with allometric weight scaling
    #    and the age power-form covariate on CL/F (Schipani 2016
    #    equations 4-5). No IIV on ka per Results body text.
    # ------------------------------------------------------------
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl) * bw_cl * age_cl
    vc     <- exp(lvc + etalvc) * bw_v
    fdepot <- exp(lfdepot) * (e_child_fdepot ^ CHILD)

    # ------------------------------------------------------------
    # 4. Micro-constant.
    # ------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------
    # 5. ODE system. One-compartment first-order absorption from
    #    depot to central, first-order elimination.
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Relative bioavailability on the dose entering the depot:
    # F = 1 for adults; F = 0.517 for children.
    f(depot) <- fdepot

    # ------------------------------------------------------------
    # 6. Observation: plasma rifampicin concentration (mg/L) with
    #    proportional residual error.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
