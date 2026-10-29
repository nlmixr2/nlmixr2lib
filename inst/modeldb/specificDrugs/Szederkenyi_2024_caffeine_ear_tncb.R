Szederkenyi_2024_caffeine_ear_tncb <- function() {
  description <- paste(
    "Ex vivo (mouse ear skin with TNCB-induced allergic contact dermatitis,",
    "skin-on-a-chip diffusion chamber). Three-compartment linear",
    "diffusion-convection model of transdermal caffeine penetration from a 2% O/W",
    "cream through allergic (TNCB-sensitised and -challenged) mouse ear tissue, over",
    "a 5-h collection. A constant zero-order input kin enters the outermost of three",
    "skin layers obtained by spatially discretising a diffusion-convection equation;",
    "neighbouring layers exchange by diffusion (kdiff, both directions) and by",
    "inward convection (kconv), and the deepest layer loses drug by first-order",
    "outflow kout. The deepest-layer state is the cumulative caffeine mass permeated",
    "per cm^2 of diffusion surface (Szederkenyi 2024 Eq 1; S1 Table row 'contact",
    "dermatitis C57BL-6J-WT-TNCB', Figure 10B). One of 13 independent least-squares",
    "fits of the same structure, one per experiment; no between-sample variability",
    "or residual error was estimated."
  )
  reference <- paste(
    "Szederkenyi G, Kocsis D, Vaghy MA, Czaran D, Sasvari P, Lengyel M, Naszlady MB,",
    "Kreis F, Antal I, Csepanyi-Komi R, Erdo F.",
    "Mathematical modeling of transdermal delivery of topical drug formulations in a",
    "dynamic microfluidic diffusion chamber in health and disease.",
    "PLoS One. 2024;19(4):e0299501. doi:10.1371/journal.pone.0299501. PMCID: PMC11008853.",
    "Model: Section 2.7, Eq 1; estimates: S1 Table (row 'contact dermatitis C57BL-6J-WT-TNCB'); fit: Figure 10B.",
    sep = " "
  )
  vignette <- "Szederkenyi_2024_transdermal_skinonachip"
  # The cream in the donor chamber is encoded as the constant zero-order input
  # kin, so the model takes no dose records.
  units <- list(time = "h", dosing = "none", amount = "ug/cm^2", concentration = "ug/cm^2")

  covariateData <- list()

  # Each layer state holds drug mass per cm^2 of diffusion surface; the paper
  # calls the states concentrations but fits x3 to the cumulative permeated
  # mass per cm^2 (Figure axes, Section 3.2).
  compartmentData <- list(
    skin_slab1 = list(analyte = "caffeine", units = "ug/cm^2", specimen = "tissue", verified = TRUE),
    skin_slab2 = list(analyte = "caffeine", units = "ug/cm^2", specimen = "tissue", verified = TRUE),
    skin_slab3 = list(analyte = "caffeine", units = "ug/cm^2", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "ex vivo (mouse ear skin, male C57BL/6J)",
    n_subjects = 3L,
    n_studies = 1L,
    system = paste(
      "Skin-on-a-chip microfluidic diffusion chamber (polydimethylsiloxane in an",
      "aluminium frame): a 0.5 cm^2 diffusion surface separates the cream-filled",
      "donor chamber from a receptor chamber perfused continuously with peripheral",
      "perfusion fluid (147 mM NaCl, 4 mM KCl, 2.3 mM CaCl2, pH 7.4) at 4 uL/min by",
      "a syringe pump (Sections 2.1 and 2.4)."
    ),
    tissue = paste(
      "Full-thickness ear (auricular) tissue of age-matched male C57BL/6J mice",
      "weighing 24.2-31.2 g (Sections 2.5 and 2.6)."
    ),
    dose_range = paste(
      "O/W cream (46% petrolatum basis ointment, 4% paraffin oil, 38% purified",
      "water, 10% propylene glycol) containing 2% active ingredient (Section 2.2)",
      "applied to the donor chamber; the amount applied is not stated. Modelled as a",
      "constant zero-order input kin."
    ),
    temperature = "not stated for the disease-model experiments",
    sampling = paste(
      "Receptor fluid collected every 30 or 15 min for 5 h; caffeine in the",
      "collected receptor fluid by NanoDrop spectrophotometry at 273 nm (Section",
      "2.4)."
    ),
    disease_state = paste(
      "TNCB-induced allergic contact dermatitis: abdominal skin sensitised with 100",
      "uL of 3% (m/V) 2,4,6-trinitro-1-chlorobenzene in acetone, both ears",
      "challenged with 20 uL of 1% (m/V) TNCB on day 5, ears excised 24 h later",
      "(Section 2.5)"
    ),
    regions = "Hungary (Budapest)",
    notes = paste(
      "Mean of n = 3 replicate skin samples (Figure 10B caption); the model was",
      "fitted to the mean curve by minimising the relative 2-norm residual V = ||ym",
      "- yc|| / ||ym|| (Eq 7) with Matlab fmincon; V_min = 0.0406 (S1 Table). S1",
      "Table labels this row Figure '10.C', but Figure 10 has only panels A and B;",
      "the TNCB curve is Figure 10B. Section 3.6: allergic contact dermatitis",
      "increased ear permeability by 30% over the acetone control",
      "(Szederkenyi_2024_caffeine_ear_acetone)."
    )
  )

  ini({
    # Estimated on the linear scale by fmincon and reported without
    # uncertainty (S1 Table); log-transformed here for positivity.
    lkin <- log(1.1352); label("Zero-order input flux into the outermost skin layer (ug/cm^2/h)") # S1 Table row 'contact dermatitis C57BL-6J-WT-TNCB' (Figure 10.C): kin = 1.1352
    lkdiff <- log(0.0347); label("Diffusion rate constant between neighbouring skin layers (1/h)") # S1 Table row 'contact dermatitis C57BL-6J-WT-TNCB': kd = 0.0347
    lkconv <- log(9.0611); label("Convection rate constant from each skin layer to the next deeper one (1/h)") # S1 Table row 'contact dermatitis C57BL-6J-WT-TNCB': kc = 9.0611
    lkout <- log(0.1236); label("First-order outflow rate constant from the deepest skin layer (1/h)") # S1 Table row 'contact dermatitis C57BL-6J-WT-TNCB': kout = 0.1236
  })

  model({
    kin <- exp(lkin)
    kdiff <- exp(lkdiff)
    kconv <- exp(lkconv)
    kout <- exp(lkout)

    # Eq 1, with x1..x3 the outermost-to-deepest layers of the discretised
    # diffusion-convection equation (paper symbols kd = kdiff, kc = kconv).
    # Diffusion runs both ways between neighbours, convection only inward;
    # all layers start empty.
    d/dt(skin_slab1) <- kin + kdiff * skin_slab2 - (kdiff + kconv) * skin_slab1
    d/dt(skin_slab2) <- (kdiff + kconv) * skin_slab1 + kdiff * skin_slab3 - (2 * kdiff + kconv) * skin_slab2
    d/dt(skin_slab3) <- (kdiff + kconv) * skin_slab2 - (kdiff + kout) * skin_slab3

    # Model output y = x3 (Section 2.7): cumulative mass permeated per cm^2
    # of diffusion surface (ug/cm^2).
    cum_mass <- skin_slab3
  })
}
