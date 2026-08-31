Ward_2026_magnesium_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 4-compartment linear ODE).",
    "Long-term magnesium (Mg(II)) homeostasis in an average healthy 70 kg adult",
    "perturbed by a biodegradable Mg (alloy) bone implant. Exchangeable ('exposed')",
    "Mg is tracked in blood serum, bone, bulk tissue (muscle + soft tissue + liver),",
    "and a local peri-implant tissue zone; Mg enters from diet (after intestinal",
    "absorption) and from constant-rate implant corrosion, and leaves by",
    "concentration-proportional urinary excretion. Deterministic: no IIV and no",
    "residual error are reported. Supports multiple/larger implants (n_implant)",
    "and dietary intake control (f_diet), and is used to predict the number of",
    "implants and the time scale required to reach hypermagnesemia.",
    sep = " "
  )
  reference <- paste(
    "Ward JP, Ahmed SK, Liu Y. Physiologically Based Pharmacokinetic Model of",
    "Magnesium Implant Absorption and Distribution in Tissue and Organs.",
    "ACS Omega. 2026 Jan 22;11(4):5144-5153. doi:10.1021/acsomega.5c06910.",
    "PMCID: PMC12878740. Parameter derivations and the full (8-ODE) pathway model",
    "are in Supporting Information A-C (ao5c06910_si_001.pdf).",
    sep = " "
  )
  vignette <- "Ward_2026_magnesium_pbpk"
  units <- list(
    time          = "day",
    dosing        = "(none -- Mg enters as two zero-order rates: diet and implant corrosion)",
    concentration = "mmol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. NOTE: unlike most nlmixr2lib models, these states are
  # CONCENTRATIONS (mmol/L), not amounts -- the published system (eqs 1-4) is
  # written directly in concentration, with the compartment volume (and its
  # exposed/unexposed capacitance factor) appearing as the divisor on the LHS.
  compartmentData <- list(
    serum = list(
      analyte = "Magnesium(II) ion, exposed (exchangeable) fraction",
      units = "mmol/L", specimen = "serum", verified = TRUE
    ),
    bone = list(
      analyte = "Magnesium(II) ion, exposed (exchangeable) fraction",
      units = "mmol/L", specimen = "tissue", verified = TRUE
    ),
    # The paper's bulk soft-tissue store T (muscle + soft tissue + liver, lumped
    # because they carry similar Mg concentrations) maps onto the canonical
    # lumped-remainder compartment `other`, not onto a compartment named
    # `tissue`: the role is defined by the lumping, and the paper never resolves
    # muscle from liver from generic soft tissue. The parameter names below keep
    # the paper's own `tissue` wording (v_tissue, q_serum_tissue, ...).
    other = list(
      analyte = "Magnesium(II) ion, exposed (exchangeable) fraction",
      units = "mmol/L", specimen = "tissue", verified = TRUE
    ),
    implant_zone = list(
      analyte = "Magnesium(II) ion, exposed (exchangeable) fraction",
      units = "mmol/L", specimen = "tissue", verified = TRUE
    )
  )

  population <- list(
    species     = "human",
    n_subjects  = NA_integer_,
    n_studies   = NA_integer_,
    age_range   = "adult",
    weight_median = "70 kg",
    disease_state = paste(
      "Healthy adult reference physiology (Table 1). Reduced renal function is",
      "represented by lowering the urinary excretion clearance `cl_renal` (with a matched",
      "reduction in `rate_diet` so pre-implant homeostasis is preserved); chronic",
      "renal failure is taken as roughly one third of the healthy value",
      "(Section 2.2.3).",
      sep = " "
    ),
    dose_range = paste(
      "Implant Mg release rate 0.05 mmol/day per 3.2 x 32 mm screw;",
      "simulations span n_implant = 1-40 screws and a scaled-up plate at",
      "10 mmol/day (Figures 2-6).",
      sep = " "
    ),
    notes = paste(
      "NOT fitted to an individual-level data set. Parameters are assembled from",
      "the physiology literature for an 'average' healthy 70 kg adult carrying",
      "1000.7 mmol of total body Mg (Supporting Information B, Tables S2-S3), then",
      "the exchange clearances were tuned so the model reproduces the stable-isotope",
      "tracer kinetics of Sojka 1997 and Sabatier 2003 (Supporting Information B,",
      "Figure S2). The authors state the values are 'representative rather than",
      "definitive'. Model validation against long-term human data is not possible:",
      "the paper notes no human study has monitored Mg over 2-12 months post-implant.",
      sep = " "
    )
  )

  ini({
    # ---- Compartment volumes (Table 1) ------------------------------------
    # Serum / RBC / bone volumes are the wet weights of Supporting Information
    # Table S2 read as litres; tissue bundles muscle + soft tissue + liver.
    v_serum <- fixed(3)
    label("Volume of the serum compartment, Vs (L)")  # Table 1, row Vs
    v_rbc <- fixed(2)
    label("Volume of the red blood cells, Vr (L)")  # Table 1, row Vr
    v_bone <- fixed(12.3)
    label("Volume of the bone compartment, VN (L)")  # Table 1, row VN
    v_tissue <- fixed(52.7)
    label("Volume of the tissue compartment, VTtot (L)")  # Table 1, row VTtot
    # VI is not measurable; chosen so VI/VTtot = 1e-4, about 20x the volume of a
    # 3.2 x 32 mm screw (Supporting Information B, S15). CI scales as 1/VI but
    # serum/bone/tissue predictions are insensitive to it.
    v_implant_zone <- fixed(0.00527)
    label("Volume of the implant zone per implant, VI (L)")  # Table 1, row VI

    # ---- Clearances (Table 1) ---------------------------------------------
    # All are volumetric exchange/excretion clearances: (L/day) x (mmol/L) =
    # mmol/day. The paper calls them "exchange rate constants" (ERC), but their
    # units are L/day, i.e. clearances, not first-order rate constants.
    lcl_renal <- fixed(log(10.9))
    label("Log urinary excretion clearance of Mg(II) from serum, gamma (L/day)")  # Table 1, row gamma; = rate_diet / Cse (Supporting Information B, S11)
    q_serum_bone <- fixed(6.05)
    label("Serum-to-bone exchange clearance, mu1 (L/day)")  # Table 1, row mu1; = 0.02 x Brown 1997 bone blood flow (Supporting Information B, S12)
    q_bone_serum <- fixed(0.0775)
    label("Bone-to-serum exchange clearance, mu-1 (L/day)")  # Table 1, row mu-1; fixed by mu-1/mu1 = Cse/CNe (Supporting Information B, S12)
    q_serum_tissue <- fixed(138)
    label("Serum-to-tissue exchange clearance, k1 (L/day)")  # Table 1, row k1; = 0.02 x Brown 1997 non-bone blood flow (Supporting Information B, S12)
    q_tissue_serum <- fixed(34.7)
    label("Tissue-to-serum exchange clearance, k-1 (L/day)")  # Table 1, row k-1; fixed by k-1/k1 = Cse/CTe = 0.252 (Supporting Information B, S12)

    # ---- Zero-order Mg inputs (Table 1) -----------------------------------
    rate_diet <- fixed(6)
    label("Dietary Mg(II) intake rate after intestinal absorption, phiD (mmol/day)")  # Table 1, row phiD; chosen to match the Sojka 1997 tracer urine jump (Supporting Information B, S11)
    rate_implant <- fixed(0.05)
    label("Mg(II) release rate per implant, sigma (mmol/day)")  # Table 1, row sigma; 3.2 x 32 mm screw, from in vivo/in vitro corrosion data (Supporting Information B, S13-S14)

    # ---- Exposed/unexposed partitioning (Table 1) -------------------------
    # "Exposed" Mg is free to exchange between compartments; "unexposed" is not.
    # The fast exposed<->unexposed equilibration is assumed instantaneous, so
    # these ratios enter only as capacitance multipliers (Supporting
    # Information A.1, eq A.9).
    f_exposed_bone <- fixed(0.3)
    label("Volume fraction of exposed (exchangeable) Mg in bone, phi (unitless)")  # Table 1, row phi; = bone exchangeable fraction eN (Supporting Information Table S3)
    kp_unexposed_serum <- fixed(0.538)
    label("Unexposed:exposed serum Mg partition, xi1 (unitless)")  # Table 1, row xi1; = (1 - es)/es with es = 0.65 (Supporting Information Table S3)
    kp_rbc_serum <- fixed(0.452)
    label("Exposed RBC:serum Mg partition, xi2 (unitless)")  # Table 1, row xi2; = er*Mr*Vs/(es*Ms*Vr) (Supporting Information B, S10)
    kp_unexposed_tissue <- fixed(3)
    label("Unexposed:exposed tissue Mg partition, xi3 (unitless)")  # Table 1, row xi3; = (1 - eT)/eT with eT = 0.25 (Supporting Information Table S3)

    # ---- Scenario controls ------------------------------------------------
    # Both are simulation levers the paper varies; they are not estimated.
    n_implant <- fixed(1)
    label("Number of implants (or surface-area-equivalent multiple of one screw), n (unitless)")  # eq 9, Section 2.2.2
    f_diet <- fixed(1)
    label("External (dietary) Mg intake control factor, rho: 1 = normal intake, 0 = none (unitless)")  # eq 10, Section 2.2.3
  })

  model({
    # ---- Multiple / larger implants (eq 9) --------------------------------
    # n non-overlapping implants scale the implant zone volume and the release
    # rate together, so xi = VI/VTtot scales with n as well.
    v_iz  <- n_implant * v_implant_zone
    xi    <- v_iz / v_tissue
    sigma <- n_implant * rate_implant

    # ---- Effective serum capacitance (LHS of eq 1) ------------------------
    # Serum exposed Mg is buffered by the unexposed serum pool (xi1) and by the
    # exposed RBC pool (xi2), both in instantaneous equilibrium with it.
    v_serum_eff <- v_serum * (1 + kp_unexposed_serum) + v_rbc * kp_rbc_serum

    cl_renal <- exp(lcl_renal)

    # ---- Pre-implant homeostasis (eq 5) -----------------------------------
    # The steady state of eqs 1-4 with sigma = 0. Note this is evaluated with
    # the UNCONTROLLED diet (rate_diet, not f_diet * rate_diet): the patient is
    # at homeostasis before the implant goes in at t = 0, and any dietary
    # control starts at t = 0 (Section 2.2.3).
    cse <- rate_diet / cl_renal
    cne <- q_serum_bone * cse / q_bone_serum
    cte <- q_serum_tissue * cse / q_tissue_serum
    serum(0)        <- cse
    bone(0)         <- cne
    other(0)        <- cte
    implant_zone(0) <- cte

    # ---- Mg(II) exchange system (eqs 1-4; eq 10 for the dietary-control
    #      form of the serum equation) --------------------------------------
    d/dt(serum) <-
      (f_diet * rate_diet -
         (cl_renal + q_serum_bone + q_serum_tissue) * serum +
         q_bone_serum * bone +
         q_tissue_serum * ((1 - xi) * other + xi * implant_zone)) / v_serum_eff
    d/dt(bone) <-
      (q_serum_bone * serum - q_bone_serum * bone) / (f_exposed_bone * v_bone)
    d/dt(other) <-
      (1 - xi) * (q_serum_tissue * serum - q_tissue_serum * other) /
        (v_tissue * (1 + kp_unexposed_tissue))
    d/dt(implant_zone) <-
      (sigma + xi * (q_serum_tissue * serum - q_tissue_serum * implant_zone)) /
        (v_iz * (1 + kp_unexposed_tissue))

    # ---- Observations ------------------------------------------------------
    # Cc is TOTAL serum Mg (exposed + unexposed) in mmol/L -- the clinically
    # measured quantity, and the one the 0.65-1.05 mmol/L reference range and
    # the 1.05 mmol/L hypermagnesemia threshold refer to (Table 1, Chom/Chyp).
    # The states themselves carry only the exposed fraction.
    Cc <- serum * (1 + kp_unexposed_serum)

    # Normalised concentrations Cj* = Cj / Cje, as plotted throughout the paper
    # (Section 2.1.1); all four equal 1 at t = 0 by construction.
    Cs_norm <- serum / cse
    CN_norm <- bone / cne
    CT_norm <- other / cte
    CI_norm <- implant_zone / cte
    # No residual-error or IIV model: the publication is a deterministic
    # mechanistic simulation and reports neither.
  })
}
