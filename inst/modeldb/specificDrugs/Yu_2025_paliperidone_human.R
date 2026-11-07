Yu_2025_paliperidone_human <- function() {
  description <- paste(
    "Semi-mechanistic intramuscular muscle-compartment model for paliperidone",
    "palmitate long-acting injectable nanocrystal suspensions in humans,",
    "extrapolated from the beagle-dog model (Yu 2025). The injected solid is split by mass into twelve",
    "particle-size bins; after a bin-specific lag time each bin dissolves into",
    "the injection-site fluid by the diffusion-layer model (particle number fixed,",
    "radius shrinking with the remaining mass, driving force = solubility minus",
    "the local dissolved concentration). Dissolved drug is carried to a",
    "one-compartment central disposition model by local blood flow. Covers the",
    "Invega Sustenna reference product and two Livzon test products (T1, T2)",
    "through formulation indicators that switch the bin particle sizes, bin lag",
    "times and bioavailability. Bin lag times, particle sizes, bioavailability",
    "and physicochemistry are carried over from the dog model; the diffusion",
    "coefficient and the injection-site clearance are scaled by the human-to-dog",
    "muscle blood-flow ratio (73/350) and disposition is from Samtani 2009.",
    "Typical-value (deterministic) model."
  )
  reference <- paste(
    "Yu Y, et al. Prediction of drug concentrations in humans for long-acting",
    "injectable suspensions by a semi-mechanical muscle compartment model: a",
    "case study of paliperidone palmitate. Front Pharmacol. 2025;16:1507828.",
    "doi:10.3389/fphar.2025.1507828"
  )
  vignette <- "Yu_2025_paliperidone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Twelve per-bin cumulative-dissolved bookkeeping states and the
  # injection-site dissolved pool are specific to this particle-binned
  # dissolution model.
  paper_specific_compartments <- c(
    paste0("diss_bin", 1:12),
    "injsite"
  )

  covariateData <- list(
    FORM_PALIPERIDONE_T1 = list(
      description = "1 = Livzon paliperidone palmitate test suspension T1, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Invega Sustenna reference suspension R, when FORM_PALIPERIDONE_T2 is also 0)",
      notes = "Switches the twelve bin particle sizes (Table 7), the twelve bin lag times (Table 2) and the bioavailability (Table 6: 83 percent) to the T1 values. Mutually exclusive with FORM_PALIPERIDONE_T2.",
      source_name = "T1"
    ),
    FORM_PALIPERIDONE_T2 = list(
      description = "1 = Livzon paliperidone palmitate test suspension T2, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Invega Sustenna reference suspension R, when FORM_PALIPERIDONE_T1 is also 0)",
      notes = "Switches the twelve bin particle sizes (Table 7), the twelve bin lag times (Table 2) and the bioavailability (Table 6: 72 percent) to the T2 values. Mutually exclusive with FORM_PALIPERIDONE_T1.",
      source_name = "T2"
    )
  )

  compartmentData <- c(
    list(
      depot = list(analyte = "paliperidone palmitate", units = "mg", specimen = "administration site", verified = TRUE)
    ),
    stats::setNames(
      rep(list(list(analyte = "paliperidone", units = "mg", specimen = "not applicable", verified = TRUE)), 12),
      paste0("diss_bin", 1:12)
    ),
    list(
      injsite = list(analyte = "paliperidone", units = "mg", specimen = "administration site", verified = TRUE),
      central = list(analyte = "paliperidone", units = "mg", specimen = "plasma", verified = TRUE)
    )
  )

  population <- list(
    species = "human",
    n_subjects = 300,
    disease_state = "virtual adults (prediction; disposition from patients with schizophrenia in Samtani 2009)",
    dose_range = "50 mg paliperidone palmitate, single intramuscular injection",
    formulations = "Invega Sustenna (reference R); Livzon test suspensions T1 and T2",
    notes = paste(
      "No human data were fitted. Section 2.6.2: 300 virtual individuals",
      "received a single 50 mg dose; the observed comparator for product R is",
      "the mean profile of Samtani 2009 (Figure 5, Table 5). The",
      "between-subject variability the authors simulated is not reported, so",
      "this file carries the typical-value model only."
    )
  )

  ini({
    # ---- Disposition (Section 3.3, Table 4) ----
    # ke = CL / Vd = 4.95 L/h / 391 L from Samtani 2009 (Section 3.3).
    lkel <- fixed(log(0.0127)); label("Elimination rate constant from central (1/h)")  # Table 4 Ke,human
    lvc <- fixed(log(391)); label("Volume of distribution (L)")  # Table 4 Vd,human
    # (73 / 350) x 35 = 7.3; the dog value is in mL/h, see the dog model.
    lcltrans <- fixed(log(7.3)); label("Injection-site clearance by local blood flow (mL/h)")  # Table 4 CLtrans '7.3 L/h' (read as mL/h)

    # ---- Injection site and physicochemistry ----
    # Carried over from the dog model (Table 1, Section 3.1.1); the
    # diffusion coefficient is scaled by the muscle blood-flow ratio.
    lvregion <- fixed(log(0.5)); label("Injection-site fluid volume (mL)")  # Section 3.1.1 depot volume 0.5 mL
    ldcoef <- fixed(log(9.78e-13)); label("Diffusion coefficient of paliperidone at the injection site (m^2/s)")  # Table 4 D; Section 3.3 (73/350) x 4.69e-12
    lhdl <- fixed(log(3)); label("Diffusion-layer thickness (um)")  # Table 1 h = 3 um (dog value carried over)
    rho <- fixed(1.2); label("Solid particle density (g/cm^3)")  # Table 1 rho
    sol <- fixed(0.00704); label("Solubility of paliperidone at the injection site (mg/mL)")  # Table 1 Cs

    # ---- Bioavailability by product (Table 6, dog) ----
    # Dog values carried over; the simulated human AUC0-t of T1 and T2 in
    # Table 5 are 0.83 and 0.72 times that of R.
    lfdepot <- fixed(log(1)); label("Bioavailability, reference product R (fraction)")  # Table 6 R 100%
    lfdepot_t1 <- fixed(log(0.83)); label("Bioavailability, test product T1 (fraction)")  # Table 6 T1 83%
    lfdepot_t2 <- fixed(log(0.72)); label("Bioavailability, test product T2 (fraction)")  # Table 6 T2 72%

    # ---- Mass fraction of the dose in each particle-size bin (Section 2.6.1) ----
    # Bins are cumulative-mass percentiles of the particle-size distribution:
    # 0-2, 2-5, 5-10, 10-20, ..., 90-100 percent of the total mass.
    frac_bin1 <- fixed(0.02); label("Fraction of the dose in particle-size bin 1 (0-2 percent of mass)")  # Section 2.6.1 bin 0-2%
    frac_bin2 <- fixed(0.03); label("Fraction of the dose in particle-size bin 2 (2-5 percent of mass)")  # Section 2.6.1 bin 2-5%
    frac_bin3 <- fixed(0.05); label("Fraction of the dose in particle-size bin 3 (5-10 percent of mass)")  # Section 2.6.1 bin 5-10%
    frac_bin4 <- fixed(0.1); label("Fraction of the dose in particle-size bin 4 (10-20 percent of mass)")  # Section 2.6.1 bin 10-20%
    frac_bin5 <- fixed(0.1); label("Fraction of the dose in particle-size bin 5 (20-30 percent of mass)")  # Section 2.6.1 bin 20-30%
    frac_bin6 <- fixed(0.1); label("Fraction of the dose in particle-size bin 6 (30-40 percent of mass)")  # Section 2.6.1 bin 30-40%
    frac_bin7 <- fixed(0.1); label("Fraction of the dose in particle-size bin 7 (40-50 percent of mass)")  # Section 2.6.1 bin 40-50%
    frac_bin8 <- fixed(0.1); label("Fraction of the dose in particle-size bin 8 (50-60 percent of mass)")  # Section 2.6.1 bin 50-60%
    frac_bin9 <- fixed(0.1); label("Fraction of the dose in particle-size bin 9 (60-70 percent of mass)")  # Section 2.6.1 bin 60-70%
    frac_bin10 <- fixed(0.1); label("Fraction of the dose in particle-size bin 10 (70-80 percent of mass)")  # Section 2.6.1 bin 70-80%
    frac_bin11 <- fixed(0.1); label("Fraction of the dose in particle-size bin 11 (80-90 percent of mass)")  # Section 2.6.1 bin 80-90%
    frac_bin12 <- fixed(0.1); label("Fraction of the dose in particle-size bin 12 (90-100 percent of mass)")  # Section 2.6.1 bin 90-100%

    # ---- Particle size per bin, reference product R (Table 7) ----
    rad_bin1 <- fixed(0.584); label("Particle radius in bin 1, reference product (um)")  # Table 7 row R, Bin 1
    rad_bin2 <- fixed(1.24); label("Particle radius in bin 2, reference product (um)")  # Table 7 row R, Bin 2
    rad_bin3 <- fixed(1.66); label("Particle radius in bin 3, reference product (um)")  # Table 7 row R, Bin 3
    rad_bin4 <- fixed(2.25); label("Particle radius in bin 4, reference product (um)")  # Table 7 row R, Bin 4
    rad_bin5 <- fixed(3.03); label("Particle radius in bin 5, reference product (um)")  # Table 7 row R, Bin 5
    rad_bin6 <- fixed(3.81); label("Particle radius in bin 6, reference product (um)")  # Table 7 row R, Bin 6
    rad_bin7 <- fixed(4.67); label("Particle radius in bin 7, reference product (um)")  # Table 7 row R, Bin 7
    rad_bin8 <- fixed(5.68); label("Particle radius in bin 8, reference product (um)")  # Table 7 row R, Bin 8
    rad_bin9 <- fixed(6.94); label("Particle radius in bin 9, reference product (um)")  # Table 7 row R, Bin 9
    rad_bin10 <- fixed(8.66); label("Particle radius in bin 10, reference product (um)")  # Table 7 row R, Bin 10
    rad_bin11 <- fixed(11.4); label("Particle radius in bin 11, reference product (um)")  # Table 7 row R, Bin 11
    rad_bin12 <- fixed(17.1); label("Particle radius in bin 12, reference product (um)")  # Table 7 row R, Bin 12

    # ---- Particle size per bin, test product T1 (Table 7) ----
    rad_t1_bin1 <- fixed(0.592); label("Particle radius in bin 1, test product T1 (um)")  # Table 7 row T1, Bin 1
    rad_t1_bin2 <- fixed(1.18); label("Particle radius in bin 2, test product T1 (um)")  # Table 7 row T1, Bin 2
    rad_t1_bin3 <- fixed(1.56); label("Particle radius in bin 3, test product T1 (um)")  # Table 7 row T1, Bin 3
    rad_t1_bin4 <- fixed(2.08); label("Particle radius in bin 4, test product T1 (um)")  # Table 7 row T1, Bin 4
    rad_t1_bin5 <- fixed(2.75); label("Particle radius in bin 5, test product T1 (um)")  # Table 7 row T1, Bin 5
    rad_t1_bin6 <- fixed(3.42); label("Particle radius in bin 6, test product T1 (um)")  # Table 7 row T1, Bin 6
    rad_t1_bin7 <- fixed(4.14); label("Particle radius in bin 7, test product T1 (um)")  # Table 7 row T1, Bin 7
    rad_t1_bin8 <- fixed(4.98); label("Particle radius in bin 8, test product T1 (um)")  # Table 7 row T1, Bin 8
    rad_t1_bin9 <- fixed(6.03); label("Particle radius in bin 9, test product T1 (um)")  # Table 7 row T1, Bin 9
    rad_t1_bin10 <- fixed(7.45); label("Particle radius in bin 10, test product T1 (um)")  # Table 7 row T1, Bin 10
    rad_t1_bin11 <- fixed(9.65); label("Particle radius in bin 11, test product T1 (um)")  # Table 7 row T1, Bin 11
    rad_t1_bin12 <- fixed(14.4); label("Particle radius in bin 12, test product T1 (um)")  # Table 7 row T1, Bin 12

    # ---- Particle size per bin, test product T2 (Table 7) ----
    rad_t2_bin1 <- fixed(0.547); label("Particle radius in bin 1, test product T2 (um)")  # Table 7 row T2, Bin 1
    rad_t2_bin2 <- fixed(1.10); label("Particle radius in bin 2, test product T2 (um)")  # Table 7 row T2, Bin 2
    rad_t2_bin3 <- fixed(1.46); label("Particle radius in bin 3, test product T2 (um)")  # Table 7 row T2, Bin 3
    rad_t2_bin4 <- fixed(1.96); label("Particle radius in bin 4, test product T2 (um)")  # Table 7 row T2, Bin 4
    rad_t2_bin5 <- fixed(2.60); label("Particle radius in bin 5, test product T2 (um)")  # Table 7 row T2, Bin 5
    rad_t2_bin6 <- fixed(3.24); label("Particle radius in bin 6, test product T2 (um)")  # Table 7 row T2, Bin 6
    rad_t2_bin7 <- fixed(3.94); label("Particle radius in bin 7, test product T2 (um)")  # Table 7 row T2, Bin 7
    rad_t2_bin8 <- fixed(4.76); label("Particle radius in bin 8, test product T2 (um)")  # Table 7 row T2, Bin 8
    rad_t2_bin9 <- fixed(5.77); label("Particle radius in bin 9, test product T2 (um)")  # Table 7 row T2, Bin 9
    rad_t2_bin10 <- fixed(7.14); label("Particle radius in bin 10, test product T2 (um)")  # Table 7 row T2, Bin 10
    rad_t2_bin11 <- fixed(9.29); label("Particle radius in bin 11, test product T2 (um)")  # Table 7 row T2, Bin 11
    rad_t2_bin12 <- fixed(13.9); label("Particle radius in bin 12, test product T2 (um)")  # Table 7 row T2, Bin 12

    # ---- Dissolution lag time per bin, reference product (Table 2, dog values carried over) ----
    tlag_bin1 <- fixed(4.8); label("Dissolution lag time of bin 1, reference product (h)")  # Table 2 row R, Bin 1
    tlag_bin2 <- fixed(21.6); label("Dissolution lag time of bin 2, reference product (h)")  # Table 2 row R, Bin 2
    tlag_bin3 <- fixed(33.6); label("Dissolution lag time of bin 3, reference product (h)")  # Table 2 row R, Bin 3
    tlag_bin4 <- fixed(52.8); label("Dissolution lag time of bin 4, reference product (h)")  # Table 2 row R, Bin 4
    tlag_bin5 <- fixed(79.2); label("Dissolution lag time of bin 5, reference product (h)")  # Table 2 row R, Bin 5
    tlag_bin6 <- fixed(84); label("Dissolution lag time of bin 6, reference product (h)")  # Table 2 row R, Bin 6
    tlag_bin7 <- fixed(84); label("Dissolution lag time of bin 7, reference product (h)")  # Table 2 row R, Bin 7
    tlag_bin8 <- fixed(72); label("Dissolution lag time of bin 8, reference product (h)")  # Table 2 row R, Bin 8
    tlag_bin9 <- fixed(72); label("Dissolution lag time of bin 9, reference product (h)")  # Table 2 row R, Bin 9
    tlag_bin10 <- fixed(72); label("Dissolution lag time of bin 10, reference product (h)")  # Table 2 row R, Bin 10
    tlag_bin11 <- fixed(72); label("Dissolution lag time of bin 11, reference product (h)")  # Table 2 row R, Bin 11
    tlag_bin12 <- fixed(72); label("Dissolution lag time of bin 12, reference product (h)")  # Table 2 row R, Bin 12

    # ---- Dissolution lag time per bin, test product T1 (Table 2, dog values carried over) ----
    tlag_t1_bin1 <- fixed(1.2); label("Dissolution lag time of bin 1, test product T1 (h)")  # Table 2 row T1, Bin 1
    tlag_t1_bin2 <- fixed(4.8); label("Dissolution lag time of bin 2, test product T1 (h)")  # Table 2 row T1, Bin 2
    tlag_t1_bin3 <- fixed(12); label("Dissolution lag time of bin 3, test product T1 (h)")  # Table 2 row T1, Bin 3
    tlag_t1_bin4 <- fixed(19.2); label("Dissolution lag time of bin 4, test product T1 (h)")  # Table 2 row T1, Bin 4
    tlag_t1_bin5 <- fixed(28.8); label("Dissolution lag time of bin 5, test product T1 (h)")  # Table 2 row T1, Bin 5
    tlag_t1_bin6 <- fixed(60); label("Dissolution lag time of bin 6, test product T1 (h)")  # Table 2 row T1, Bin 6
    tlag_t1_bin7 <- fixed(72); label("Dissolution lag time of bin 7, test product T1 (h)")  # Table 2 row T1, Bin 7
    tlag_t1_bin8 <- fixed(72); label("Dissolution lag time of bin 8, test product T1 (h)")  # Table 2 row T1, Bin 8
    tlag_t1_bin9 <- fixed(72); label("Dissolution lag time of bin 9, test product T1 (h)")  # Table 2 row T1, Bin 9
    tlag_t1_bin10 <- fixed(72); label("Dissolution lag time of bin 10, test product T1 (h)")  # Table 2 row T1, Bin 10
    tlag_t1_bin11 <- fixed(144); label("Dissolution lag time of bin 11, test product T1 (h)")  # Table 2 row T1, Bin 11
    tlag_t1_bin12 <- fixed(144); label("Dissolution lag time of bin 12, test product T1 (h)")  # Table 2 row T1, Bin 12

    # ---- Dissolution lag time per bin, test product T2 (Table 2, dog values carried over) ----
    tlag_t2_bin1 <- fixed(4.8); label("Dissolution lag time of bin 1, test product T2 (h)")  # Table 2 row T2, Bin 1
    tlag_t2_bin2 <- fixed(19.2); label("Dissolution lag time of bin 2, test product T2 (h)")  # Table 2 row T2, Bin 2
    tlag_t2_bin3 <- fixed(28.8); label("Dissolution lag time of bin 3, test product T2 (h)")  # Table 2 row T2, Bin 3
    tlag_t2_bin4 <- fixed(48); label("Dissolution lag time of bin 4, test product T2 (h)")  # Table 2 row T2, Bin 4
    tlag_t2_bin5 <- fixed(64.8); label("Dissolution lag time of bin 5, test product T2 (h)")  # Table 2 row T2, Bin 5
    tlag_t2_bin6 <- fixed(84); label("Dissolution lag time of bin 6, test product T2 (h)")  # Table 2 row T2, Bin 6
    tlag_t2_bin7 <- fixed(96); label("Dissolution lag time of bin 7, test product T2 (h)")  # Table 2 row T2, Bin 7
    tlag_t2_bin8 <- fixed(96); label("Dissolution lag time of bin 8, test product T2 (h)")  # Table 2 row T2, Bin 8
    tlag_t2_bin9 <- fixed(120); label("Dissolution lag time of bin 9, test product T2 (h)")  # Table 2 row T2, Bin 9
    tlag_t2_bin10 <- fixed(120); label("Dissolution lag time of bin 10, test product T2 (h)")  # Table 2 row T2, Bin 10
    tlag_t2_bin11 <- fixed(120); label("Dissolution lag time of bin 11, test product T2 (h)")  # Table 2 row T2, Bin 11
    tlag_t2_bin12 <- fixed(120); label("Dissolution lag time of bin 12, test product T2 (h)")  # Table 2 row T2, Bin 12
  })
  model({
    kel <- exp(lkel)
    vc <- exp(lvc)
    cltrans <- exp(lcltrans)
    vregion <- exp(lvregion)
    hdl <- exp(lhdl)

    # Unit conversions to a um / mg / h frame: D from m^2/s to um^2/h,
    # density from g/cm^3 to mg/um^3, solubility and the local dissolved
    # concentration from mg/mL to mg/um^3.
    dcoef_um <- exp(ldcoef) * 3.6e15
    rho_um <- rho * 1e-9
    sol_um <- sol * 1e-12
    csite <- injsite / vregion
    csite_um <- csite * 1e-12

    # Formulation selection: reference product R unless a test-product
    # indicator is set.
    form_r <- 1 - FORM_PALIPERIDONE_T1 - FORM_PALIPERIDONE_T2
    r0_1 <- form_r * rad_bin1 + FORM_PALIPERIDONE_T1 * rad_t1_bin1 + FORM_PALIPERIDONE_T2 * rad_t2_bin1
    r0_2 <- form_r * rad_bin2 + FORM_PALIPERIDONE_T1 * rad_t1_bin2 + FORM_PALIPERIDONE_T2 * rad_t2_bin2
    r0_3 <- form_r * rad_bin3 + FORM_PALIPERIDONE_T1 * rad_t1_bin3 + FORM_PALIPERIDONE_T2 * rad_t2_bin3
    r0_4 <- form_r * rad_bin4 + FORM_PALIPERIDONE_T1 * rad_t1_bin4 + FORM_PALIPERIDONE_T2 * rad_t2_bin4
    r0_5 <- form_r * rad_bin5 + FORM_PALIPERIDONE_T1 * rad_t1_bin5 + FORM_PALIPERIDONE_T2 * rad_t2_bin5
    r0_6 <- form_r * rad_bin6 + FORM_PALIPERIDONE_T1 * rad_t1_bin6 + FORM_PALIPERIDONE_T2 * rad_t2_bin6
    r0_7 <- form_r * rad_bin7 + FORM_PALIPERIDONE_T1 * rad_t1_bin7 + FORM_PALIPERIDONE_T2 * rad_t2_bin7
    r0_8 <- form_r * rad_bin8 + FORM_PALIPERIDONE_T1 * rad_t1_bin8 + FORM_PALIPERIDONE_T2 * rad_t2_bin8
    r0_9 <- form_r * rad_bin9 + FORM_PALIPERIDONE_T1 * rad_t1_bin9 + FORM_PALIPERIDONE_T2 * rad_t2_bin9
    r0_10 <- form_r * rad_bin10 + FORM_PALIPERIDONE_T1 * rad_t1_bin10 + FORM_PALIPERIDONE_T2 * rad_t2_bin10
    r0_11 <- form_r * rad_bin11 + FORM_PALIPERIDONE_T1 * rad_t1_bin11 + FORM_PALIPERIDONE_T2 * rad_t2_bin11
    r0_12 <- form_r * rad_bin12 + FORM_PALIPERIDONE_T1 * rad_t1_bin12 + FORM_PALIPERIDONE_T2 * rad_t2_bin12
    tlag1 <- form_r * tlag_bin1 + FORM_PALIPERIDONE_T1 * tlag_t1_bin1 + FORM_PALIPERIDONE_T2 * tlag_t2_bin1
    tlag2 <- form_r * tlag_bin2 + FORM_PALIPERIDONE_T1 * tlag_t1_bin2 + FORM_PALIPERIDONE_T2 * tlag_t2_bin2
    tlag3 <- form_r * tlag_bin3 + FORM_PALIPERIDONE_T1 * tlag_t1_bin3 + FORM_PALIPERIDONE_T2 * tlag_t2_bin3
    tlag4 <- form_r * tlag_bin4 + FORM_PALIPERIDONE_T1 * tlag_t1_bin4 + FORM_PALIPERIDONE_T2 * tlag_t2_bin4
    tlag5 <- form_r * tlag_bin5 + FORM_PALIPERIDONE_T1 * tlag_t1_bin5 + FORM_PALIPERIDONE_T2 * tlag_t2_bin5
    tlag6 <- form_r * tlag_bin6 + FORM_PALIPERIDONE_T1 * tlag_t1_bin6 + FORM_PALIPERIDONE_T2 * tlag_t2_bin6
    tlag7 <- form_r * tlag_bin7 + FORM_PALIPERIDONE_T1 * tlag_t1_bin7 + FORM_PALIPERIDONE_T2 * tlag_t2_bin7
    tlag8 <- form_r * tlag_bin8 + FORM_PALIPERIDONE_T1 * tlag_t1_bin8 + FORM_PALIPERIDONE_T2 * tlag_t2_bin8
    tlag9 <- form_r * tlag_bin9 + FORM_PALIPERIDONE_T1 * tlag_t1_bin9 + FORM_PALIPERIDONE_T2 * tlag_t2_bin9
    tlag10 <- form_r * tlag_bin10 + FORM_PALIPERIDONE_T1 * tlag_t1_bin10 + FORM_PALIPERIDONE_T2 * tlag_t2_bin10
    tlag11 <- form_r * tlag_bin11 + FORM_PALIPERIDONE_T1 * tlag_t1_bin11 + FORM_PALIPERIDONE_T2 * tlag_t2_bin11
    tlag12 <- form_r * tlag_bin12 + FORM_PALIPERIDONE_T1 * tlag_t1_bin12 + FORM_PALIPERIDONE_T2 * tlag_t2_bin12

    # Total solid mass injected so far (mg): what is still undissolved plus
    # what every bin has already released. Each bin starts with frac_bin<n>
    # of it, as particles of the bin's radius; the particle count is set by
    # that starting mass and stays constant while the particles shrink.
    minj <- depot + diss_bin1 + diss_bin2 + diss_bin3 + diss_bin4 + diss_bin5 + diss_bin6 +
      diss_bin7 + diss_bin8 + diss_bin9 + diss_bin10 + diss_bin11 + diss_bin12

    # Time since the injection; NA before the first dose.
    tinj <- tad(depot)
    if (is.na(tinj)) tinj <- -1

    m0_1 <- frac_bin1 * minj
    xrem1 <- 0
    if (m0_1 > 0) xrem1 <- (m0_1 - diss_bin1) / m0_1
    if (xrem1 < 0) xrem1 <- 0
    npart1 <- m0_1 / (rho_um * 4 / 3 * pi * r0_1^3)
    rad1 <- r0_1 * xrem1^(1 / 3)
    rdiss1 <- 0
    if (tinj > tlag1) rdiss1 <- npart1 * dcoef_um / hdl * 4 * pi * rad1 * (rad1 + hdl) * (sol_um - csite_um)
    m0_2 <- frac_bin2 * minj
    xrem2 <- 0
    if (m0_2 > 0) xrem2 <- (m0_2 - diss_bin2) / m0_2
    if (xrem2 < 0) xrem2 <- 0
    npart2 <- m0_2 / (rho_um * 4 / 3 * pi * r0_2^3)
    rad2 <- r0_2 * xrem2^(1 / 3)
    rdiss2 <- 0
    if (tinj > tlag2) rdiss2 <- npart2 * dcoef_um / hdl * 4 * pi * rad2 * (rad2 + hdl) * (sol_um - csite_um)
    m0_3 <- frac_bin3 * minj
    xrem3 <- 0
    if (m0_3 > 0) xrem3 <- (m0_3 - diss_bin3) / m0_3
    if (xrem3 < 0) xrem3 <- 0
    npart3 <- m0_3 / (rho_um * 4 / 3 * pi * r0_3^3)
    rad3 <- r0_3 * xrem3^(1 / 3)
    rdiss3 <- 0
    if (tinj > tlag3) rdiss3 <- npart3 * dcoef_um / hdl * 4 * pi * rad3 * (rad3 + hdl) * (sol_um - csite_um)
    m0_4 <- frac_bin4 * minj
    xrem4 <- 0
    if (m0_4 > 0) xrem4 <- (m0_4 - diss_bin4) / m0_4
    if (xrem4 < 0) xrem4 <- 0
    npart4 <- m0_4 / (rho_um * 4 / 3 * pi * r0_4^3)
    rad4 <- r0_4 * xrem4^(1 / 3)
    rdiss4 <- 0
    if (tinj > tlag4) rdiss4 <- npart4 * dcoef_um / hdl * 4 * pi * rad4 * (rad4 + hdl) * (sol_um - csite_um)
    m0_5 <- frac_bin5 * minj
    xrem5 <- 0
    if (m0_5 > 0) xrem5 <- (m0_5 - diss_bin5) / m0_5
    if (xrem5 < 0) xrem5 <- 0
    npart5 <- m0_5 / (rho_um * 4 / 3 * pi * r0_5^3)
    rad5 <- r0_5 * xrem5^(1 / 3)
    rdiss5 <- 0
    if (tinj > tlag5) rdiss5 <- npart5 * dcoef_um / hdl * 4 * pi * rad5 * (rad5 + hdl) * (sol_um - csite_um)
    m0_6 <- frac_bin6 * minj
    xrem6 <- 0
    if (m0_6 > 0) xrem6 <- (m0_6 - diss_bin6) / m0_6
    if (xrem6 < 0) xrem6 <- 0
    npart6 <- m0_6 / (rho_um * 4 / 3 * pi * r0_6^3)
    rad6 <- r0_6 * xrem6^(1 / 3)
    rdiss6 <- 0
    if (tinj > tlag6) rdiss6 <- npart6 * dcoef_um / hdl * 4 * pi * rad6 * (rad6 + hdl) * (sol_um - csite_um)
    m0_7 <- frac_bin7 * minj
    xrem7 <- 0
    if (m0_7 > 0) xrem7 <- (m0_7 - diss_bin7) / m0_7
    if (xrem7 < 0) xrem7 <- 0
    npart7 <- m0_7 / (rho_um * 4 / 3 * pi * r0_7^3)
    rad7 <- r0_7 * xrem7^(1 / 3)
    rdiss7 <- 0
    if (tinj > tlag7) rdiss7 <- npart7 * dcoef_um / hdl * 4 * pi * rad7 * (rad7 + hdl) * (sol_um - csite_um)
    m0_8 <- frac_bin8 * minj
    xrem8 <- 0
    if (m0_8 > 0) xrem8 <- (m0_8 - diss_bin8) / m0_8
    if (xrem8 < 0) xrem8 <- 0
    npart8 <- m0_8 / (rho_um * 4 / 3 * pi * r0_8^3)
    rad8 <- r0_8 * xrem8^(1 / 3)
    rdiss8 <- 0
    if (tinj > tlag8) rdiss8 <- npart8 * dcoef_um / hdl * 4 * pi * rad8 * (rad8 + hdl) * (sol_um - csite_um)
    m0_9 <- frac_bin9 * minj
    xrem9 <- 0
    if (m0_9 > 0) xrem9 <- (m0_9 - diss_bin9) / m0_9
    if (xrem9 < 0) xrem9 <- 0
    npart9 <- m0_9 / (rho_um * 4 / 3 * pi * r0_9^3)
    rad9 <- r0_9 * xrem9^(1 / 3)
    rdiss9 <- 0
    if (tinj > tlag9) rdiss9 <- npart9 * dcoef_um / hdl * 4 * pi * rad9 * (rad9 + hdl) * (sol_um - csite_um)
    m0_10 <- frac_bin10 * minj
    xrem10 <- 0
    if (m0_10 > 0) xrem10 <- (m0_10 - diss_bin10) / m0_10
    if (xrem10 < 0) xrem10 <- 0
    npart10 <- m0_10 / (rho_um * 4 / 3 * pi * r0_10^3)
    rad10 <- r0_10 * xrem10^(1 / 3)
    rdiss10 <- 0
    if (tinj > tlag10) rdiss10 <- npart10 * dcoef_um / hdl * 4 * pi * rad10 * (rad10 + hdl) * (sol_um - csite_um)
    m0_11 <- frac_bin11 * minj
    xrem11 <- 0
    if (m0_11 > 0) xrem11 <- (m0_11 - diss_bin11) / m0_11
    if (xrem11 < 0) xrem11 <- 0
    npart11 <- m0_11 / (rho_um * 4 / 3 * pi * r0_11^3)
    rad11 <- r0_11 * xrem11^(1 / 3)
    rdiss11 <- 0
    if (tinj > tlag11) rdiss11 <- npart11 * dcoef_um / hdl * 4 * pi * rad11 * (rad11 + hdl) * (sol_um - csite_um)
    m0_12 <- frac_bin12 * minj
    xrem12 <- 0
    if (m0_12 > 0) xrem12 <- (m0_12 - diss_bin12) / m0_12
    if (xrem12 < 0) xrem12 <- 0
    npart12 <- m0_12 / (rho_um * 4 / 3 * pi * r0_12^3)
    rad12 <- r0_12 * xrem12^(1 / 3)
    rdiss12 <- 0
    if (tinj > tlag12) rdiss12 <- npart12 * dcoef_um / hdl * 4 * pi * rad12 * (rad12 + hdl) * (sol_um - csite_um)

    # Equation 1 (diffusion-layer dissolution of each bin), Equation 2
    # (injection-site dissolved drug) and Equation 3 (central compartment).
    d/dt(depot) <- -(rdiss1 + rdiss2 + rdiss3 + rdiss4 + rdiss5 + rdiss6 +
      rdiss7 + rdiss8 + rdiss9 + rdiss10 + rdiss11 + rdiss12)
    d/dt(diss_bin1) <- rdiss1
    d/dt(diss_bin2) <- rdiss2
    d/dt(diss_bin3) <- rdiss3
    d/dt(diss_bin4) <- rdiss4
    d/dt(diss_bin5) <- rdiss5
    d/dt(diss_bin6) <- rdiss6
    d/dt(diss_bin7) <- rdiss7
    d/dt(diss_bin8) <- rdiss8
    d/dt(diss_bin9) <- rdiss9
    d/dt(diss_bin10) <- rdiss10
    d/dt(diss_bin11) <- rdiss11
    d/dt(diss_bin12) <- rdiss12
    d/dt(injsite) <- (rdiss1 + rdiss2 + rdiss3 + rdiss4 + rdiss5 + rdiss6 +
      rdiss7 + rdiss8 + rdiss9 + rdiss10 + rdiss11 + rdiss12) - csite * cltrans
    d/dt(central) <- csite * cltrans - kel * central

    f(depot) <- exp(form_r * lfdepot + FORM_PALIPERIDONE_T1 * lfdepot_t1 + FORM_PALIPERIDONE_T2 * lfdepot_t2)

    # Dissolved paliperidone concentration at the injection site (ng/mL),
    # the human counterpart of Figure 4B.
    Cinj <- csite * 1e6
    # Plasma paliperidone (ng/mL): mg / L = ug/mL, times 1000.
    Cc <- central / vc * 1000
  })
}
