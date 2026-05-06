Jonsson_2011_ethambutol <- function() {
  description <- "Two-compartment population PK model for oral ethambutol in adult South African pulmonary tuberculosis patients (Jonsson 2011), with one transit compartment preceding first-order absorption, allometric scaling on clearance and volume terms (theory-based exponents on a 50 kg reference), an HIV-status effect on bioavailability, and 4-occasion inter-occasion variability on apparent oral clearance."
  reference <- paste(
    "Jonsson S, Davidse A, Wilkins J, Van der Walt JS, Simonsson US,",
    "Karlsson MO, Smith P, McIlleron H. (2011). Population pharmacokinetics",
    "of ethambutol in South African tuberculosis patients.",
    "Antimicrob Agents Chemother 55(9):4230-7.",
    "doi:10.1128/AAC.00274-11.",
    "DDMORE Foundation Model Repository: DDMODEL00000220.",
    sep = " "
  )
  vignette <- "Jonsson_2011_ethambutol"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000220"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with theory-based exponents (0.75 on CL and Q, 1.0 on Vc and Vp) and a 50 kg reference weight. Mean WT 47 kg, range 29-86 kg per the DDMODEL00000220 RDF model-has-description-long abstract.",
      source_name        = "WT"
    ),
    HIV_POS = list(
      description        = "HIV-1 antibody-positive comorbidity indicator (1 = HIV-positive, 0 = HIV-negative).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = "12% HIV-positive in the Jonsson 2011 cohort per the RDF abstract. Multiplicative shift on bioavailability `f(transit1) <- 1 + e_hiv_pos_f * HIV_POS`; the source THETA(9) = -0.155 corresponds to a 15.5% reduction in ethambutol bioavailability for HIV-positive subjects.",
      source_name        = "HIV"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4 identify the dosing / sampling occasion within subject. Decomposed inside `model()` into binary indicators `oc1` .. `oc4` that multiplex the 4 IOV etas on log-CL.",
      source_name        = "OCC"
    )
  )

  population <- list(
    n_subjects     = 189L,
    n_studies      = 2L,
    age_range      = "16-72 years (mean 36)",
    age_median     = "36 years (mean reported)",
    weight_range   = "29-86 kg (mean 47)",
    weight_median  = "47 kg (mean reported)",
    sex_female_pct = 46,
    disease_state  = "Adults with pulmonary tuberculosis. 12% HIV-positive (HIV is a within-cohort comorbidity rather than the primary indication).",
    dose_range     = "Oral ethambutol 800-1500 mg daily, multiple-dose at steady state, combined with a standard antitubercular backbone.",
    regions        = "South Africa (two centers).",
    notes          = "Population descriptors are reproduced from the DDMODEL00000220 RDF `model-has-description-long` abstract, which mirrors the Jonsson 2011 paper's Methods. Estimated baseline creatinine clearance 79 mL/min (range 23-150 mL/min); renal function was not retained as a PK covariate. The Jonsson 2011 publication itself is not on disk in this worktree, so the demographics here come from the RDF abstract rather than the paper's Table 1; see the validation vignette's Errata section for the full caveat list."
  )

  ini({
    # Structural PK parameters - DDMODEL00000220 Output_real_run32150.lst FINAL PARAMETER ESTIMATE
    # block (THETA vector), captured after the `MINIMIZATION SUCCESSFUL` marker. Reference weight
    # is 50 kg (Jonsson 2011 .mod line 41-46: `(WT/50)**0.75` on CL/Q, `(WT/50)**1` on V2/V3).
    lcl    <- log(39.9)  ; label("Apparent oral clearance CL/F at WT = 50 kg (L/h)")               # THETA(1) FINAL = 3.99E+01
    lvc    <- log(82.4)  ; label("Apparent central volume of distribution V2/F at WT = 50 kg (L)") # THETA(2) FINAL = 8.24E+01
    lka    <- log(0.474) ; label("Absorption rate constant from depot to central, ka (1/h)")       # THETA(3) FINAL = 4.74E-01
    lvp    <- log(623)   ; label("Apparent peripheral volume of distribution V3/F at WT = 50 kg (L)") # THETA(6) FINAL = 6.23E+02
    lq     <- log(34.3)  ; label("Apparent inter-compartmental clearance Q/F at WT = 50 kg (L/h)") # THETA(7) FINAL = 3.43E+01
    lmtt   <- log(0.789) ; label("Mean transit time through the absorption transit compartment, MTT (h)") # THETA(8) FINAL = 7.89E-01

    # HIV covariate effect on bioavailability. f(transit1) = 1 + e_hiv_pos_f * HIV_POS.
    e_hiv_pos_f <- -0.154 ; label("HIV-positive multiplicative effect on bioavailability (fractional change)") # THETA(9) FINAL = -1.54E-01

    # Inter-individual variability (IIV). NONMEM `$OMEGA` diagonals (with V2, V3, Q
    # all `0 FIX` in the .mod, so no IIV is estimated on those parameters).
    etalcl  ~ 0.0381  # OMEGA(1,1) FINAL = 3.81E-02 — IIV CL (log-normal variance)
    etalka  ~ 0.153   # OMEGA(3,3) FINAL = 1.53E-01 — IIV Ka  (log-normal variance)
    etalmtt ~ 0.862   # OMEGA(6,6) FINAL = 8.62E-01 — IIV MTT (log-normal variance)

    # Inter-occasion variability (IOV) on log-CL. Source `$OMEGA BLOCK(1)` followed by
    # three `BLOCK(1) SAME` re-uses the same single-element variance across the four
    # occasions; the FINAL estimate of that variance is OMEGA(7,7) = 1.27E-01. nlmixr2
    # has no `SAME` shortcut, so each occasion gets its own eta with the variance fixed
    # to the shared value after the first (matching the Xie_2019_agomelatine pattern).
    etaiov_cl_1 ~ 0.127            # OMEGA(7,7)  FINAL = 1.27E-01 — estimated occasion-1 IOV variance
    etaiov_cl_2 ~ fix(0.127)       # OMEGA(8,8)  fixed equal to OMEGA(7,7) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_3 ~ fix(0.127)       # OMEGA(9,9)  fixed equal to OMEGA(7,7) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_4 ~ fix(0.127)       # OMEGA(10,10) fixed equal to OMEGA(7,7) per `$OMEGA BLOCK(1) SAME`

    # Combined residual error on the linear (mg/L) scale. The Jonsson .mod uses log-
    # transformed observations (`Y = LOG(F) + W*EPS(1)` with `W = SQRT(THETA(4)**2 +
    # (THETA(5)/F)**2)` and `$SIGMA 1 FIX`); on the back-transformed linear scale this
    # is the combined-error variance `var(eps_lin) = (THETA(4)*F)^2 + THETA(5)^2`,
    # i.e. nlmixr2's default Pythagorean-SD `combined2` form with proportional SD =
    # THETA(4) and additive SD = THETA(5).
    propSd <- 0.318  ; label("Proportional residual error (fraction)") # THETA(4) FINAL = 3.18E-01
    addSd  <- 0.107  ; label("Additive residual error (mg/L)")         # THETA(5) FINAL = 1.07E-01
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for IOV
    # multiplexing on log-CL (matches the .mod `OC1..OC4` derivation in $PK).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Individual PK parameters with allometric scaling on a 50 kg reference weight
    # (theory-based exponents 0.75 on CL/Q and 1.0 on V2/V3 per the Jonsson 2011
    # .mod $PK block lines 41-46).
    cl  <- exp(lcl  + etalcl  + iov_cl) * (WT / 50)^0.75
    vc  <- exp(lvc)                     * (WT / 50)
    ka  <- exp(lka  + etalka)
    vp  <- exp(lvp)                     * (WT / 50)
    q   <- exp(lq)                      * (WT / 50)^0.75
    mtt <- exp(lmtt + etalmtt)
    ktr <- 1 / mtt

    # Two-compartment oral PK with one transit compartment preceding first-order
    # absorption. Dose lands in `transit1` (the .mod's COMP=(TRANSIT) is compartment
    # 1, the NONMEM dosing compartment, with F1 = 1 * FCOV applied there); the
    # transit compartment empties into `depot` at rate `ktr`, which then absorbs
    # into `central` at rate `ka`. The 2-cmt central / peripheral block reproduces
    # the .mod K34 = Q/V2 / K43 = Q/V3 micro-constants and K30 = CL/V2 elimination.
    d/dt(transit1)    <- -ktr * transit1
    d/dt(depot)       <-  ktr * transit1 - ka * depot
    d/dt(central)     <-  ka  * depot   - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q   / vc      * central - q / vp * peripheral1

    # Bioavailability on the dosing (transit1) compartment. The .mod sets
    # F1 = 1 * FCOV with FCOV = 1 + THETA(9) on HIV-positive subjects (THETA(9) = -0.154
    # → 15.4% lower bioavailability on HIV-positive); HIV-negative reference keeps F1 = 1.
    f(transit1) <- 1 + e_hiv_pos_f * HIV_POS

    # Concentration in plasma. Dose units mg, Vc units L → Cc units mg/L (= ug/mL).
    Cc <- central / vc

    # Combined proportional + additive residual error (default Pythagorean / combined2
    # form, matching the linearized form of the .mod's log-transform-both-sides $ERROR).
    Cc ~ add(addSd) + prop(propSd)
  })
}
