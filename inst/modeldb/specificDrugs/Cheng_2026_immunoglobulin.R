Cheng_2026_immunoglobulin <- function() {
  description <- "Two-compartment population PK model for intravenous immunoglobulin (IVIG) replacement therapy in pediatric primary-immunodeficiency and secondary-antibody-deficiency patients (Cheng 2026)"
  reference <- "Cheng IL, Huang ZH, Worth A, Booth C, Standing JF. Pharmacokinetic modelling of intravenous immunoglobulin in children with primary immunodeficiencies and secondary antibody deficiencies. Br J Clin Pharmacol. 2025;1-11. doi:10.1002/bcp.70420"
  vignette <- "Cheng_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying or baseline; used for theory-based allometric scaling on CL, Q (exponent 0.75) and on V1, V2 (exponent 1) with reference weight 70 kg. Cheng 2026 fixed both exponents in the final model after a base-model evaluation found the empirically estimated exponents (0.788 for CL, 0.743 for V) included the theory-based values within their 95% CIs but caused parameter collinearity.",
      source_name        = "WT"
    ),
    DIS_SAD = list(
      description        = "Secondary antibody deficiency indicator (1 = SAD, 0 = primary immunodeficiency)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (primary immunodeficiency, PID)",
      notes              = "Time-fixed per subject. In the Cheng 2026 cohort the SAD subgroup (n = 20 / 64) is 75% post-rituximab and 25% post-CAR-T cell therapy. Multiplicative theta^DIS_SAD effects on both CL and on the baseline IgG (CBAS).",
      source_name        = "DIS_SAD"
    ),
    IGM = list(
      description        = "Baseline serum immunoglobulin M (IgM) concentration (proxy for B-cell humoral capacity)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline-only. Power scaling (IGM / 0.21)^0.11 on the baseline IgG (CBAS); 0.21 g/L is the pooled PID + SAD pediatric-cohort median in Cheng 2026 Table 2.",
      source_name        = "IGM"
    )
  )

  population <- list(
    n_subjects     = 64L,
    n_studies      = 1L,
    age_range      = "0.06-16.8 years (3 weeks to 16.8 years)",
    age_median     = "4.08 years",
    weight_range   = "3.15-95.3 kg",
    weight_median  = "18.6 kg",
    sex_female_pct = round(100 * (64 - 41) / 64, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Pediatric patients with primary immunodeficiency (PID, n = 44) or secondary antibody deficiency (SAD, n = 20; 15 post-rituximab, 5 post-CAR-T cell therapy) on intravenous immunoglobulin replacement therapy",
    dose_range     = "0.24-1.38 g/kg IV every 21 or 28 days (median 0.56 g/kg per dose)",
    regions        = "United Kingdom (single tertiary pediatric centre)",
    notes          = "Retrospective electronic-health-record analysis of children treated with intravenous Ig at a tertiary paediatric hospital between April 2019 and April 2024 (Cheng 2026 Table 2). 444 plasma IgG samples, predominantly trough; lower limit of quantification 0.07 g/L. PID patients received 0.3 g/kg every 3 weeks; SAD patients received 0.5 g/kg every 4 weeks. The IVIG products were Privigen (n = 53), Octagam (n = 9), and Gamunex (n = 2) — Cheng 2026 Table 3."
  )

  ini({
    # Structural parameters — typical values for a 70 kg PID patient at the
    # population median IgM = 0.21 g/L. Cheng 2026 Table 4 (page 7).
    lcl   <- log(0.308); label("CL for a 70 kg PID patient (L/day)")                                 # Cheng 2026 Table 4
    lvc   <- log(3.59);  label("Central volume of distribution for a 70 kg patient (V1, L)")        # Cheng 2026 Table 4
    lq    <- log(1.08);  label("Intercompartmental clearance for a 70 kg patient (Q, L/day)")        # Cheng 2026 Table 4
    lvp   <- log(7.37);  label("Peripheral volume of distribution for a 70 kg patient (V2, L)")     # Cheng 2026 Table 4
    lcbas <- log(5.67);  label("Baseline endogenous IgG for a typical PID patient at IGM = 0.21 g/L (CBAS, g/L)")  # Cheng 2026 Table 4

    # Allometric exponents — fixed to theory-based values in the final model
    # after a base-model evaluation found the empirically estimated exponents
    # included 0.75 / 1.0 within their 95% CIs but caused parameter collinearity
    # (Cheng 2026 Methods, page 7). Q follows the CL exponent.
    allo_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Cheng 2026 Methods, fixed in final model
    allo_v  <- fixed(1.0);  label("Allometric exponent on V1 and V2 (unitless)") # Cheng 2026 Methods, fixed in final model

    # Covariate effects — multiplicative theta^DIS_SAD form for the binary
    # disease-type covariate (PID = 0, SAD = 1) and a power form (IGM/0.21)^e
    # for the continuous IgM covariate. Cheng 2026 Methods (eqs 2-3) and
    # Table 4. The SAD/PID ratio on CL (0.542) implies SAD CL is ~46% lower
    # than PID CL on the same allometric weight basis; the discussion text
    # describes this as "SAD clearance being half that of PID."
    e_sad_cl    <- 0.542; label("SAD-vs-PID multiplicative ratio on CL (CL_SAD / CL_PID)")          # Cheng 2026 Table 4
    e_sad_cbas  <- 0.541; label("SAD-vs-PID multiplicative ratio on baseline IgG (CBAS_SAD / CBAS_PID)")  # Cheng 2026 Table 4
    e_igm_cbas  <- 0.11;  label("Power exponent for IgM on baseline IgG: (IGM/0.21)^e_igm_cbas")     # Cheng 2026 Table 4

    # Inter-individual variability (omega^2 = log(CV^2 + 1)).
    # Cheng 2026 reports a "variance-covariance matrix" parameterisation but
    # only the diagonal CV% values appear in Table 4; off-diagonal covariances
    # are unavailable, so the etas are modelled as independent (no block).
    etalcl   ~ 0.16749  # 42.7% CV; Cheng 2026 Table 4 IIV CL
    etalvp   ~ 1.07192  # 138.6% CV; Cheng 2026 Table 4 IIV V2 (IIV is on V2 / vp, not V1)
    etalcbas ~ 0.21757  # 49.3% CV; Cheng 2026 Table 4 IIV CBAS

    # Combined residual error on observed total IgG (Cheng 2026 Table 4).
    addSd  <- 0.812;  label("Additive residual error on total IgG (g/L)")    # Cheng 2026 Table 4
    propSd <- 0.117;  label("Proportional residual error on total IgG (fraction)")  # Cheng 2026 Table 4
  })
  model({
    # Disease-type effect on CL and on baseline IgG (CBAS). PID = reference.
    cl_dis   <- e_sad_cl^DIS_SAD
    cbas_dis <- e_sad_cbas^DIS_SAD

    # IgM effect on baseline IgG (power scaling, reference IGM = 0.21 g/L).
    cbas_igm <- (IGM / 0.21)^e_igm_cbas

    # Individual PK and baseline parameters with theory-based allometric
    # weight scaling (reference 70 kg). IIV is on CL, V2 (peripheral volume),
    # and CBAS only — V1 (vc) and Q have typical-only values per Cheng 2026
    # Table 4.
    cl   <- exp(lcl   + etalcl)   * (WT / 70)^allo_cl * cl_dis
    vc   <- exp(lvc)              * (WT / 70)^allo_v
    q    <- exp(lq)               * (WT / 70)^allo_cl
    vp   <- exp(lvp   + etalvp)   * (WT / 70)^allo_v
    cbas <- exp(lcbas + etalcbas) * cbas_dis * cbas_igm

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system. The `central` state holds the *exogenous*
    # (therapeutic) IgG amount; endogenous IgG is added at the observation
    # step via `cbas`.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    # Dose in g, volume in L -> g/L; cbas in g/L. Cheng 2026 Methods page 5:
    # "measured IgG was assumed to be the sum of endogenous IgG, the baseline
    # IgG (CBAS) level prior to treatment and exogenous therapeutic Ig."
    Cc <- central / vc + cbas
    Cc ~ add(addSd) + prop(propSd)
  })
}
