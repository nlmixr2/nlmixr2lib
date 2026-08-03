Gebhard_2023_mercaptopurine <- function() {
  description <- "Oral 6-mercaptopurine (6MP) population PK during maintenance therapy for childhood acute lymphoblastic leukemia, predicting the erythrocyte 6-thioguanine-nucleotide concentration E-TGN (Gebhard 2023, model PK^6MP_{mm,pop}). The plasma layer is a FIXED one-compartment model with first-order absorption and a fixed bioavailability of 0.12, all taken from the literature (Supplementary Table S2). A third compartment carries 6-thioguanine nucleotides -- the active metabolites of 6MP -- inside red blood cells, filled by SATURABLE Michaelis-Menten influx driven by the plasma 6MP concentration and drained by a first-order efflux rate constant. V_mm and K_eff were estimated with IIV on 452 children with 4624 E-TGN observations (Table 3); K_mm was estimated as a population-only parameter with no IIV, which is what the 'pop' subscript of the model name denotes. Because every patient's record begins mid-maintenance-therapy, the red-cell compartment is initialised at the patient's first observed E-TGN value (covariate BL_TGN_RBC) rather than at zero. Volumes are per m^2 body-surface area and doses are supplied in umol/m^2, so BSA cancels and no BSA covariate is needed."
  reference <- paste(
    "Gebhard A, Lilienthal P, Metzler M, Rauh M, Sager S, Schmiegelow K, Toksvang LN, Zierk J. (2023).",
    "Pharmacokinetic-pharmacodynamic modeling of maintenance therapy for childhood acute lymphoblastic leukemia.",
    "Sci Rep 13:11749. doi:10.1038/s41598-023-38414-0.",
    "Fixed plasma-PK parameters are reproduced in Gebhard 2023 Supplementary Table S2 from",
    "Lennard L, Keen D, Lilleyman JS (1986) and Brunton LL, Lazo JS, Parker K (2005).",
    "The joint PKPD model that adds the Friberg myelosuppression chain to this 6MP arm is",
    "modellib('Gebhard_2023_mercaptopurine_anc'); the methotrexate companion is",
    "modellib('Gebhard_2023_methotrexate').",
    sep = " "
  )
  vignette <- "Gebhard_2023_leukemia_maintenance_therapy"
  units <- list(
    time          = "day",
    dosing        = "umol/m^2",
    concentration = "umol/L"
  )
  # Unit note. Gebhard 2023 states that the E-TGN OBSERVATIONS were converted to
  # umol/L (assuming a hemoglobin molecular weight of 64458 g/mol and 330 g Hb/L
  # erythrocytes) but never states the units of the PLASMA concentrations. They
  # are fixed by dimensional analysis of equation 6MP3b: V_mm is reported in
  # umol/L/day and K_mm in umol/L (Table 3), and K_mm is added directly to the
  # plasma concentration X_C / V_C in the denominator, so the plasma
  # concentration must also be in umol/L. Volumes are L/m^2 (Supplementary
  # Table S2), so doses must be supplied in umol/m^2 and BSA cancels out of the
  # model entirely. Convert a clinical mg/m^2 dose with the 6-mercaptopurine
  # molecular weight 152.18 g/mol:
  #   amt [umol/m^2] = dose [mg/m^2] * 1000 / 152.18.
  # That conversion factor is a derived quantity, not printed in the paper.

  covariateData <- list(
    BL_TGN_RBC = list(
      description        = "Baseline (first observed) erythrocyte 6-thioguanine-nucleotide concentration, used as the initial condition of the rbc_tgn compartment.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Gebhard 2023 Methods: 'we again assume the influx and efflux of 6MP into and out of the red blood cells to not influence plasma kinetics and initialize the compartments similar to the MTX PK model with the first observation INITGN', giving X_E^6MP(0) = INITGN. Every patient's record begins after weeks of maintenance therapy, 28 days after the last high-dose intravenous MTX course, so the red-cell pool has already accumulated and cannot be initialised at zero. Cohort median 0.83 umol/L, range 0-7.6 umol/L (Table 1). Originally measured in nmol/mmol hemoglobin and converted to umol/L assuming Hb MW 64458 g/mol and 330 g Hb/L erythrocytes (Gebhard 2023 Methods).",
      source_name        = "INITGN"
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
    dose_range     = "6-Mercaptopurine 5.4-175.0 mg/m^2 daily by mouth (median 57.1 mg/m^2); the concurrent methotrexate dose was 1.3-45.0 mg/m^2 weekly (median 15.0 mg/m^2). Doses were titrated to a target white-blood-cell count of 1.5-3.5 G/L.",
    regions        = "Nordic countries (NOPHO ALL-92)",
    observations   = "4624 E-TGN observations across 452 patients; the full data set also holds 4192 E-MTX and 9808 ANC observations (Table 1).",
    notes          = "Each patient's record begins with the first paired E-TGN / E-MTX measurement taken 28 days after the last high-dose intravenous methotrexate course, so the high-dose-MTX period is excluded. Patients were dropped when parameter estimation would have been impossible: fewer than two E-TGN or E-MTX observations, no 6MP or MTX dose, or no height, weight or ANC observation. Estimation used NONMEM 7.5 with FOCEi for this PK submodel."
  )

  ini({
    # =========================================================================
    # Plasma layer -- FIXED one-compartment model, Gebhard 2023 Supplementary
    # Table S2. The paper reports micro-constants (ka, ke, V_C); they are
    # reparameterised here to the canonical lka / lcl / lvc. The product is
    # written out inline so the arithmetic is auditable.
    # =========================================================================
    lka     <- fixed(log(21.07));           label("ka: first-order absorption rate constant (1/day)")                # Gebhard 2023 Suppl Table S2: ka^6MP = 21.07 1/day
    lvc     <- fixed(log(20.1911));         label("V_C: central volume per body-surface area (L/m^2)")               # Gebhard 2023 Suppl Table S2: V_C^6MP = 20.1911 L/m^2
    lcl     <- fixed(log(15.4 * 20.1911));  label("CL: total clearance per body-surface area (L/day/m^2) = ke * V_C") # Gebhard 2023 Suppl Table S2: ke^6MP = 15.4 1/day, V_C^6MP = 20.1911 L/m^2
    lfdepot <- fixed(log(0.12));            label("F: oral bioavailability of 6-mercaptopurine (fraction)")          # Gebhard 2023 Suppl Table S2: F^6MP = 0.12

    # =========================================================================
    # Red-cell layer -- ESTIMATED on the whole data set, Gebhard 2023 Table 3
    # (model PK^6MP_{mm,pop}). Values are the whole-data-set column, not the
    # cross-validation column and not initial estimates.
    # =========================================================================
    lvmax_rbc <- log(0.096); label("V_mm: maximum saturable influx rate, plasma -> red cells (umol/L/day)")          # Gebhard 2023 Table 3: V_mm^6MP  = 0.096 umol/L/day (RSE 8%)
    lkm_rbc   <- log(0.016); label("K_mm: plasma 6MP concentration at half-maximal red-cell influx (umol/L)")        # Gebhard 2023 Table 3: K_mm^6MP  = 0.016 umol/L (RSE 10%)
    lkeff_rbc <- log(0.041); label("K_eff: first-order efflux / elimination rate constant out of red cells (1/day)") # Gebhard 2023 Table 3: K_eff^6MP = 0.041 1/day (RSE 8%)

    # IIV. Table 3 reports the coefficient of variation, and its footnote gives
    # the definition CV = sqrt(exp(omega^2) - 1); inverting, omega^2 =
    # log(CV^2 + 1). The inversion is written out inline so it is auditable.
    # There is deliberately NO eta on K_mm: the 'pop' subscript of the model
    # name PK^6MP_{mm,pop} denotes exactly that, and Table 3 leaves the
    # K_mm^6MP row of the CV block empty.
    etalvmax_rbc ~ log(0.52^2 + 1)  # Gebhard 2023 Table 3: CV(V_mm^6MP)  = 52% [44, 59], eta-shrinkage 26% -> omega^2 = 0.23933
    etalkeff_rbc ~ log(0.50^2 + 1)  # Gebhard 2023 Table 3: CV(K_eff^6MP) = 50% [42, 57], eta-shrinkage 30% -> omega^2 = 0.22314

    # -------------------------------------------------------------------------
    # Residual unexplained variability on E-TGN. Gebhard 2023 does not state
    # whether the tabulated numbers are variances or standard deviations; they
    # are VARIANCES, and this arm is what decides the question. Read as
    # variances the combined residual SD at the cohort median E-TGN of
    # 0.83 umol/L is sqrt(0.023 + (0.25298 * 0.83)^2) = 0.259 umol/L against an
    # observed median / mean individual-prediction RMSE of 0.21 / 0.24
    # (Table 3) -- a direct match. Read as standard deviations it would be
    # sqrt(0.023^2 + (0.064 * 0.83)^2) = 0.058, about 3.6 times too small.
    # The implied CVs are also the deciding evidence: as variances they are
    # 25.3% here and 15.5% for E-MTX, plausible for a monthly red-cell
    # metabolite assay; as standard deviations they would be 6.4% and 2.4%.
    # The square roots are written out inline so the reading is auditable.
    # -------------------------------------------------------------------------
    addSd_Crbc_tgn  <- sqrt(0.023); label("Additive residual SD on E-TGN (umol/L)")        # Gebhard 2023 Table 3: additive residual error 0.023 (RSE 9%), read as a variance
    propSd_Crbc_tgn <- sqrt(0.064); label("Proportional residual SD on E-TGN (fraction)")  # Gebhard 2023 Table 3: proportional residual error 0.064 (RSE 5%), read as a variance
  })

  model({
    # ---- 1. Individual parameters ----
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    vmax_rbc <- exp(lvmax_rbc + etalvmax_rbc)
    km_rbc   <- exp(lkm_rbc)                    # population-only: no eta (the 'pop' variant)
    keff_rbc <- exp(lkeff_rbc + etalkeff_rbc)

    # ---- 2. Micro-constants ----
    kel <- cl / vc   # = ke^6MP = 15.4 1/day

    # ---- 3. ODE system: equations 6MP1, 6MP2 and 6MP3b ----
    # Cc must be defined before the red-cell ODE that consumes it.
    # (Note: the main text's typeset equation 6MP2 prints its left-hand side as
    # dX_GI/dt, a superscript / subscript garble. Supplementary Information
    # section I and the mass balance both give dX_C/dt, used here.)
    d/dt(depot)   <- -ka * depot                # 6MP1
    d/dt(central) <-  ka * depot - kel * central # 6MP2
    Cc <- central / vc
    # 6MP3b -- SATURABLE Michaelis-Menten influx, which is the variant the
    # final model PK^6MP_{mm,pop} uses (the 'mm' subscript names that choice
    # over the linear alternative 6MP3a).
    d/dt(rbc_tgn) <- vmax_rbc * Cc / (km_rbc + Cc) - keff_rbc * rbc_tgn

    # ---- 4. Bioavailability ----
    f(depot) <- exp(lfdepot)

    # ---- 5. Initial conditions ----
    # X_GI(0) = X_C(0) = 0 (rxode2 default). The red-cell pool starts at the
    # patient's first observed E-TGN because the record begins weeks into
    # maintenance therapy (Gebhard 2023 Methods, 6MP section).
    rbc_tgn(0) <- BL_TGN_RBC

    # ---- 6. Observations ----
    # Cc (defined above) is the plasma 6-mercaptopurine concentration
    # (umol/L). It is NOT observed in this study -- Gebhard 2023 fits only the
    # red-cell concentration -- but it is exposed so downstream users and the
    # validation vignette can inspect and NCA-check the fixed plasma layer that
    # drives the red-cell influx.
    # E-TGN, the observed endpoint: 6-thioguanine nucleotides inside red cells.
    Crbc_tgn <- rbc_tgn
    Crbc_tgn ~ add(addSd_Crbc_tgn) + prop(propSd_Crbc_tgn)
  })
}
