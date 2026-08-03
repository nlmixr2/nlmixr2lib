Gebhard_2023_methotrexate <- function() {
  description <- "Low-dose oral methotrexate (MTX) population PK during maintenance therapy for childhood acute lymphoblastic leukemia, predicting the erythrocyte methotrexate concentration E-MTX (Gebhard 2023, model PK^MTX_{fix,bio,lin,peri}). The plasma layer is a FIXED two-compartment model with first-order absorption whose parameters were taken from the literature (Supplementary Table S1) and a saturable dose-dependent bioavailability F = 1 - 0.77 * D / (15.01 + D) with D the MTX dose in mg/m^2. A fourth compartment carries methotrexate (predominantly methotrexate polyglutamates) inside red blood cells, filled by LINEAR influx driven by the PERIPHERAL plasma concentration and drained by a first-order efflux rate constant; both red-cell rate constants were estimated on 452 children with 4192 E-MTX observations (Table 2) and carry IIV. Because every patient's record begins mid-maintenance-therapy, the red-cell compartment is initialised at the patient's first observed E-MTX value (covariate BL_MTX_RBC) rather than at zero. Volumes are per m^2 body-surface area and doses are supplied in umol/m^2, so BSA cancels and no BSA covariate is needed."
  reference <- paste(
    "Gebhard A, Lilienthal P, Metzler M, Rauh M, Sager S, Schmiegelow K, Toksvang LN, Zierk J. (2023).",
    "Pharmacokinetic-pharmacodynamic modeling of maintenance therapy for childhood acute lymphoblastic leukemia.",
    "Sci Rep 13:11749. doi:10.1038/s41598-023-38414-0.",
    "Fixed plasma-PK parameters are reproduced in Gebhard 2023 Supplementary Table S1 from",
    "Panetta JC, Sparreboom A, Pui C-H, Relling MV, Evans WE (2010) and the sources cited there.",
    "The 6-mercaptopurine companion models are modellib('Gebhard_2023_mercaptopurine') and",
    "modellib('Gebhard_2023_mercaptopurine_anc').",
    sep = " "
  )
  vignette <- "Gebhard_2023_leukemia_maintenance_therapy"
  units <- list(
    time          = "day",
    dosing        = "umol/m^2",
    concentration = "umol/L"
  )
  # Unit note. Gebhard 2023 states that the E-MTX OBSERVATIONS were converted to
  # umol/L (assuming a hemoglobin molecular weight of 64458 g/mol and 330 g Hb/L
  # erythrocytes) but never states the units of the PLASMA concentrations. They
  # are fixed by dimensional analysis of equation MTX4c,
  #   dX_E/dt = K_in * (X_P / V_P) - K_eff * X_E,
  # in which K_in and K_eff are both reported in 1/day (Table 2). For
  # K_in * C_plasma to carry the units of dX_E/dt (umol/L/day) when K_in is
  # 1/day, C_plasma must be in umol/L -- the same units as X_E. Volumes are
  # L/m^2 (Supplementary Table S1), so doses must be supplied in umol/m^2 and
  # BSA cancels out of the model entirely. Convert a clinical mg/m^2 dose with
  # the methotrexate molecular weight 454.44 g/mol:
  #   amt [umol/m^2] = dose [mg/m^2] * 1000 / 454.44.
  # That conversion factor is a derived quantity, not printed in the paper.

  covariateData <- list(
    BL_MTX_RBC = list(
      description        = "Baseline (first observed) erythrocyte methotrexate concentration, used as the initial condition of the rbc_mtx compartment.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Gebhard 2023 Methods: 'X_E^MTX(0) = INIMTX with INIMTX being the first observation in the data set at time point 0'. Every patient's record begins after weeks of maintenance therapy, 28 days after the last high-dose intravenous MTX course, so the red-cell pool has already accumulated and cannot be initialised at zero. Cohort median 0.026 umol/L, range 0-0.10 umol/L (Table 1). The paper's own Discussion names an alternative remedy -- restructure the timeline so estimation can start where the red-cell compartment is genuinely zero -- but that is not the model as published. Originally measured in nmol/mmol hemoglobin and converted to umol/L assuming Hb MW 64458 g/mol and 330 g Hb/L erythrocytes (Gebhard 2023 Methods).",
      source_name        = "INIMTX"
    ),
    DOSE_MTX_MGM2 = list(
      description        = "Methotrexate dose administered on the current dose record, normalised to body-surface area.",
      units              = "mg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the saturable dose-dependent bioavailability F^MTX = 1 - 0.77 * DOSE_MTX_MGM2 / (15.01 + DOSE_MTX_MGM2) (Gebhard 2023 Supplementary Table S1, model PK^MTX_bio). A dedicated mg/m^2 column is required because this model's amt is supplied in umol/m^2; convert with the methotrexate molecular weight 454.44 g/mol via amt [umol/m^2] = DOSE_MTX_MGM2 * 1000 / 454.44. Cohort median weekly dose 15.0 mg/m^2, range 1.3-45.0 mg/m^2 (Table 1). At the median dose F = 1 - 0.77 * 15 / 30.01 = 0.615.",
      source_name        = "MTX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 452L,
    n_studies      = 1L,
    age_range      = "2.4-16.9 years (median 5.9)",
    weight_range   = "10.3-105.5 kg (median 21.5)",
    height_range   = "81.5-180.0 cm (median 114.0)",
    disease_state  = "Precursor-B-cell acute lymphoblastic leukemia in children, during oral maintenance therapy with 6-mercaptopurine plus low-dose oral methotrexate. Treated on the NOPHO ALL-92 protocol (Nordic Society for Paediatric Haematology and Oncology); Copenhagen ethics approval V.200.2080/91.",
    dose_range     = "Methotrexate 1.3-45.0 mg/m^2 weekly by mouth (median 15.0 mg/m^2); the concurrent 6-mercaptopurine dose was 5.4-175.0 mg/m^2 daily (median 57.1 mg/m^2). Doses were titrated to a target white-blood-cell count of 1.5-3.5 G/L.",
    regions        = "Nordic countries (NOPHO ALL-92)",
    observations   = "4192 E-MTX observations across 452 patients; the full data set also holds 4624 E-TGN and 9808 ANC observations (Table 1).",
    notes          = "Each patient's record begins with the first paired E-TGN / E-MTX measurement taken 28 days after the last high-dose intravenous methotrexate course, so the high-dose-MTX period is excluded. Patients were dropped when parameter estimation would have been impossible: fewer than two E-TGN or E-MTX observations, no 6MP or MTX dose, or no height, weight or ANC observation. Estimation used NONMEM 7.5 with FOCEi for this PK submodel."
  )

  ini({
    # =========================================================================
    # Plasma layer -- FIXED two-compartment model, Gebhard 2023 Supplementary
    # Table S1 (model PK^MTX_fix). The paper reports micro-constants
    # (ka, ke, k_cp, k_pc, V_C, V_P); they are reparameterised here to the
    # canonical lka / lcl / lvc / lq / lvp. The mapping is exact and is
    # verified by the paper's own internal consistency check: Table S1 states
    # V_P was computed as V_C * k12 / k21 = 11.606 * 0.9888 / 2.64 = 4.3466,
    # so q = k_cp * V_C reproduces k_pc = q / V_P = 2.64000 1/day exactly.
    # The products are written out inline so the arithmetic is auditable.
    # =========================================================================
    lka <- fixed(log(26.64));                label("ka: first-order absorption rate constant (1/day)")                  # Gebhard 2023 Suppl Table S1: ka^MTX = 26.64 1/day
    lvc <- fixed(log(11.606));               label("V_C: central volume per body-surface area (L/m^2)")                 # Gebhard 2023 Suppl Table S1: V_C^MTX = 11.606 L/m^2
    lvp <- fixed(log(4.347));                label("V_P: peripheral volume per body-surface area (L/m^2)")               # Gebhard 2023 Suppl Table S1: V_P^MTX = 4.347 L/m^2
    lcl <- fixed(log(14.2944 * 11.606));     label("CL: total clearance per body-surface area (L/day/m^2) = ke * V_C")   # Gebhard 2023 Suppl Table S1: ke^MTX = 14.2944 1/day, V_C^MTX = 11.606 L/m^2
    lq  <- fixed(log(0.9888 * 11.606));      label("Q: intercompartmental clearance per body-surface area (L/day/m^2) = k_cp * V_C") # Gebhard 2023 Suppl Table S1: k_cp = 0.9888 1/day, V_C^MTX = 11.606 L/m^2

    # -------------------------------------------------------------------------
    # Dose-dependent bioavailability (model variant PK^MTX_bio),
    # Gebhard 2023 Suppl Table S1: F^MTX = 1 - 0.77 * MTX / (15.01 + MTX)
    # with MTX the dose in mg/m^2. Encoded in the library's established
    # Imax / D50 form (cf. VasquezBahena_2009_lumiracoxib_rat.R): lfdepot is
    # the F0 anchor fixed at 1 and the saturable reduction is applied on top.
    # All three constants are fixed literature values, not estimated here.
    # -------------------------------------------------------------------------
    lfdepot      <- fixed(log(1));     label("F0: bioavailability anchor before the dose-dependent reduction (fraction)") # Gebhard 2023 Suppl Table S1: the printed formula is 1 - 0.77 * MTX / (15.01 + MTX), i.e. an anchor of 1
    limax_fdepot <- fixed(log(0.77));  label("Imax: maximum fractional reduction in bioavailability (dimensionless)")     # Gebhard 2023 Suppl Table S1: numerator coefficient 0.77
    ld50_fdepot  <- fixed(log(15.01)); label("D50: MTX dose giving half the maximal bioavailability reduction (mg/m^2)")  # Gebhard 2023 Suppl Table S1: denominator constant 15.01 mg/m^2

    # =========================================================================
    # Red-cell layer -- ESTIMATED on the whole data set, Gebhard 2023 Table 2
    # (model PK^MTX_{fix,bio,lin,peri}). Values are the whole-data-set column,
    # not the cross-validation column and not initial estimates.
    # =========================================================================
    lkinf_rbc <- log(0.031); label("K_in: first-order influx rate constant, peripheral plasma -> red cells (1/day)")   # Gebhard 2023 Table 2: K_in^MTX  = 0.031 1/day (RSE 5%)
    lkeff_rbc <- log(0.018); label("K_eff: first-order efflux / elimination rate constant out of red cells (1/day)")   # Gebhard 2023 Table 2: K_eff^MTX = 0.018 1/day (RSE 5%); the Discussion cross-checks this against literature red-cell MTX half-lives of 30-40 days (0.017-0.023 1/day)

    # IIV. Table 2 reports the coefficient of variation, and its footnote gives
    # the definition CV = sqrt(exp(omega^2) - 1); inverting, omega^2 =
    # log(CV^2 + 1). The inversion is written out inline so it is auditable.
    etalkinf_rbc ~ log(0.34^2 + 1)  # Gebhard 2023 Table 2: CV(K_in^MTX)  = 34% [29, 38], eta-shrinkage 29% -> omega^2 = 0.10938
    etalkeff_rbc ~ log(0.31^2 + 1)  # Gebhard 2023 Table 2: CV(K_eff^MTX) = 31% [26, 35], eta-shrinkage 32% -> omega^2 = 0.09176

    # -------------------------------------------------------------------------
    # Residual unexplained variability on E-MTX. Gebhard 2023 does not state
    # whether the tabulated numbers are variances or standard deviations; they
    # are VARIANCES. Evidence (also recorded in the vignette Errata): read as
    # variances the combined residual SD at the cohort median E-MTX of
    # 0.026 umol/L is sqrt(0.000018 + (0.1549 * 0.026)^2) = 0.0058 umol/L
    # against an observed median individual-prediction RMSE of 0.0042 (Table 2)
    # -- the right order of magnitude. Read as standard deviations it would be
    # sqrt(0.000018^2 + (0.024 * 0.026)^2) = 0.00062, seven times too small.
    # The 6-mercaptopurine arm decides the question even more sharply; see
    # Gebhard_2023_mercaptopurine.R. The square roots are written out inline.
    # -------------------------------------------------------------------------
    addSd_Crbc_mtx  <- sqrt(0.000018); label("Additive residual SD on E-MTX (umol/L)")        # Gebhard 2023 Table 2: additive residual error 0.000018 (RSE 7%), read as a variance
    propSd_Crbc_mtx <- sqrt(0.024);    label("Proportional residual SD on E-MTX (fraction)")  # Gebhard 2023 Table 2: proportional residual error 0.024 (RSE 8%), read as a variance
  })

  model({
    # ---- 1. Dose-dependent bioavailability (Suppl Table S1) ----
    imax_fdepot <- exp(limax_fdepot)
    d50_fdepot  <- exp(ld50_fdepot)
    fdepot <- exp(lfdepot) *
      (1 - imax_fdepot * DOSE_MTX_MGM2 / (d50_fdepot + DOSE_MTX_MGM2))

    # ---- 2. Individual parameters ----
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kinf_rbc <- exp(lkinf_rbc + etalkinf_rbc)
    keff_rbc <- exp(lkeff_rbc + etalkeff_rbc)

    # ---- 3. Micro-constants ----
    kel <- cl / vc   # = ke^MTX  = 14.2944 1/day
    k12 <- q  / vc   # = k_cp    =  0.9888 1/day
    k21 <- q  / vp   # = k_pc    =  2.6400 1/day

    # ---- 4. ODE system: equations MTX1, MTX2, MTX3 and MTX4c ----
    # MTX4c is the LINEAR-influx-from-the-PERIPHERAL-compartment variant, which
    # is the one the final model PK^MTX_{fix,bio,lin,peri} uses (the subscripts
    # lin,peri name exactly that choice among the four MTX4a-d variants).
    d/dt(depot)       <- -ka * depot                                                   # MTX1
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1  # MTX2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1                            # MTX3
    d/dt(rbc_mtx)     <-  kinf_rbc * (peripheral1 / vp) - keff_rbc * rbc_mtx            # MTX4c

    # ---- 5. Bioavailability ----
    f(depot) <- fdepot

    # ---- 6. Initial conditions ----
    # X_GI(0) = X_C(0) = X_P(0) = 0 (rxode2 default). The red-cell pool starts
    # at the patient's first observed E-MTX because the record begins weeks
    # into maintenance therapy (Gebhard 2023 Methods, MTX section).
    rbc_mtx(0) <- BL_MTX_RBC

    # ---- 7. Observations ----
    # Cc is the plasma methotrexate concentration (umol/L). It is NOT observed
    # in this study -- Gebhard 2023 fits only the red-cell concentration -- but
    # it is exposed so downstream users and the validation vignette can inspect
    # and NCA-check the fixed plasma layer that drives the red-cell influx.
    Cc   <- central / vc
    # E-MTX, the observed endpoint: methotrexate inside red blood cells.
    Crbc_mtx <- rbc_mtx
    Crbc_mtx ~ add(addSd_Crbc_mtx) + prop(propSd_Crbc_mtx)
  })
}
