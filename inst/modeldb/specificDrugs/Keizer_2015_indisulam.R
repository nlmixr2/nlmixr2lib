Keizer_2015_indisulam <- function() {
  description <- "Two-compartment linear IV population PK model for indisulam (E7070) in 34 adult solid-tumour patients receiving 250-525 mg/m2 as a 2-hour infusion in combination with irinotecan (Keizer 2015, Table 3, 'All data' column). This is the real-data model from a methodological study comparing four ways of handling concentrations below the limit of quantification; the 'All data' column is the paper's advocated method, in which extrapolated concentrations between the limit of detection and the LLOQ are used as continuous observations. Interindividual variability is on clearance only; no covariates were investigated. Indisulam is known to have nonlinear (saturable) disposition, but nonlinearity was not supported by this limited data set and the fitted model is linear."
  reference <- paste(
    "Keizer RJ, Jansen RS, Rosing H, Thijssen B, Beijnen JH,",
    "Schellens JHM, Huitema ADR.",
    "Incorporation of concentration data below the limit of quantification",
    "in population pharmacokinetic analyses.",
    "Pharmacol Res Perspect. 2015;3(2):e00131.",
    "doi:10.1002/prp2.131.",
    sep = " "
  )
  vignette <- "Keizer_2015_indisulam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "indisulam", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "indisulam", units = "mg", specimen = "plasma", verified = FALSE)
  )

  # No covariates. Keizer 2015 Methods, 'Real PopPK data set': "we only
  # performed basic compartmental PK modeling on this limited data set, and,
  # for example, the influence of covariates was not investigated."
  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 34L,
    n_studies        = 1L,
    disease_state    = paste(
      "Adult patients with advanced solid tumours enrolled in a phase I",
      "dose-escalation trial of indisulam (E7070) in combination with",
      "irinotecan (Ryan et al. 2005); the trial supplying the PK data was",
      "sponsored by Eisai."
    ),
    dose_range       = "250-525 mg/m^2 indisulam administered as a 2-hour intravenous infusion",
    n_concentrations = 231L,
    notes            = paste(
      "One PK curve per patient, sampled over 120 h (Keizer 2015,",
      "'Real PopPK data set'). Of 231 PK samples (excluding pre-first-dose",
      "samples), 17 (7.4%) were below the LLOQ. Plasma indisulam was assayed",
      "by a validated LC-MS/MS method over 0.1-20 microg/mL (Beumer et al.",
      "2004); the limit of detection was taken as 30% of the LLOQ. A",
      "two-compartment linear model fitted better than a one-compartment",
      "model; additional peripheral compartments and nonlinear elimination",
      "were not supported by this data set, and visual predictive checks",
      "revealed no relevant model misspecification. Estimation used NONMEM",
      "VI level 2.0 with the Laplacian method. Age, sex, weight, race and",
      "other demographics were not reported in this paper (the trial",
      "demographics are in Ryan et al. 2005) and no covariates were tested.",
      "A more extensive semiphysiological model for indisulam, including its",
      "nonlinear disposition, was published separately (Zandvliet et al.",
      "2006) and is NOT the model encoded here."
    )
  )

  ini({
    # ===== Structural PK (Keizer 2015 Table 3, 'All data' column) =====
    # The four Table 3 columns are the SAME structural model re-fitted under
    # four different BLQ-handling methods ('Discard', 'LLOQ/2', 'LIKE',
    # 'All data'). The 'All data' column is the paper's advocated method and
    # the one with by far the best estimation stability (91.5% successful
    # minimization, 85.7% successful covariance step vs 66.7%/11.4% for
    # 'LIKE'), so it is taken as the final model. The other three columns are
    # reproduced in the validation vignette as a sensitivity comparison; all
    # four agree to within ~1% on the structural parameters.
    lcl <- log(0.828); label("Typical clearance CL (L/h)")                          # Keizer 2015 Table 3, row 'CL (L/h)', 'All data' column = 0.828 (Discard 0.823, LLOQ/2 0.819, LIKE 0.822)
    lvc <- log(5.63);  label("Typical central volume V (L)")                        # Keizer 2015 Table 3, row 'V (L)', 'All data' column = 5.63 (Discard 5.61, LLOQ/2 5.65, LIKE 5.61)
    lq  <- log(1.88);  label("Typical intercompartmental clearance Q (L/h)")        # Keizer 2015 Table 3, row 'Q (L/h)', 'All data' column = 1.88 (Discard 1.90, LLOQ/2 1.89, LIKE 1.9)
    lvp <- log(12.6);  label("Typical peripheral volume Vper (L)")                  # Keizer 2015 Table 3, row 'Vper (L)', 'All data' column = 12.6 (Discard 12.4, LLOQ/2 12.5, LIKE 12.4)

    # ===== Interindividual variability (Keizer 2015 Table 3, 'All data') =====
    # Table 3 labels this row 'eta_CL' and gives it as a percentage, in the
    # same percentage style as the 'sigma_prop' row immediately below it.
    # sigma_prop is unambiguously the proportional residual standard
    # deviation itself (0.267 -> prop(0.267)), so by parallel construction
    # eta_CL = 60.6% is the standard deviation of the log-scale random
    # effect, i.e. omega = 0.606 and omega^2 = 0.606^2 = 0.367236. The
    # paper's own simulation section uses the same convention, defining
    # simulated BSV as '25%' (i.e. omega = 0.25, $OMEGA 0.0625). Encoding it
    # instead as an exact log-normal CV, log(1 + 0.606^2) = 0.3127, would
    # give omega = 0.559 and is inconsistent with the sigma_prop row.
    # IIV was estimated on clearance only; no IIV on V, Q or Vper is
    # reported in Table 3.
    etalcl ~ 0.367236  # Keizer 2015 Table 3, row 'eta_CL', 'All data' column = 60.6%; 0.606^2 = 0.367236

    # ===== Residual error (Keizer 2015 Table 3, 'All data' column) =====
    # Combined additive + proportional. Note the additive term differs
    # noticeably between BLQ methods (0.0364 for 'All data' vs 0.071 for
    # 'Discard') -- that spread is the paper's point, since the additive
    # term is the one the low concentrations inform.
    propSd <- 0.267;  label("Proportional residual error (fraction)")  # Keizer 2015 Table 3, row 'sigma_prop', 'All data' column = 26.7%
    addSd  <- 0.0364; label("Additive residual error (mg/L)")          # Keizer 2015 Table 3, row 'sigma_add (mg/L)', 'All data' column = 0.0364
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # Indisulam is given as a 2-hour intravenous infusion directly into the
    # central compartment (NONMEM ADVAN3 in the source).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # ----- Output -----
    # Plasma indisulam concentration: dose in mg, vc in L -> mg/L
    # (equivalently microg/mL, the unit of the LC-MS/MS assay range).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
