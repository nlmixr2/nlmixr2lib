Hyland_2008_maraviroc_rcyp3a4 <- function() {
  description <- "In vitro (recombinant human CYP3A4 Supersomes). Michaelis-Menten enzyme-kinetic model of the CYP3A4-mediated N-dealkylation of maraviroc to its secondary-amine metabolite UK-408,027 in heterologously expressed CYP3A4, plus the parallel CYP3A4-mediated oxidative routes that make up the remainder of maraviroc's intrinsic clearance in the same system. The UK-408,027 route accounts for Vmax/Km = 0.23 uL/min/pmol CYP3A4 of the 1.7 uL/min/pmol CYP3A4 total depletion intrinsic clearance; the balance is assigned to the other oxidative pathways as a first-order route, because substrate depletion was measured only at 1 uM, far below Km, and so characterises a linear clearance. The model also carries the intersystem extrapolation factor of 0.2 with which the paper scales the recombinant intrinsic clearance to a human liver microsome basis, reproducing the 0.34 uL/min/pmol and 40.8 uL/min/mg values that were the CLint input to the maraviroc Simcyp model. No impurity term is carried, because in the recombinant system the substrate contaminant was subtracted by comparison against control Supersomes rather than fitted. Sibling model: Hyland_2008_maraviroc_hlm, the same reaction characterised in pooled human liver microsomes."
  reference <- "Hyland R, Dickins M, Collins C, Jones H, Jones B. Maraviroc: in vitro assessment of drug-drug interaction potential. Br J Clin Pharmacol. 2008 Oct;66(4):498-507. doi:10.1111/j.1365-2125.2008.03198.x. PMID: 18492127. PMCID: PMC2561101. Km and Vmax estimates: Results, 'Kinetics of maraviroc N-dealkylation by rCYP3A4', and Figure 4. Substrate-depletion intrinsic clearance, the intersystem extrapolation factor, the CYP3A4 content of the human liver microsome batch and the fraction unbound in microsomes: Results, 'CLint estimates from HLM and rCYP', and Table 1. Incubation design and the control-Supersome contaminant subtraction: Materials and methods, 'Metabolism of maraviroc by expressed recombinant CYPs', and Results, 'Metabolism of maraviroc in expressed recombinant CYPs'. The intersystem extrapolation factor is taken from Proctor NJ, Tucker GT, Rostami-Hodjegan A. Xenobiotica 2004;34:151-78, reference 13 of the source paper."
  vignette <- "Hyland_2008_maraviroc_invitro"
  units <- list(time = "min", dosing = "uM (incubation concentration)", concentration = "uM")

  # The "dose" of this static in vitro system is the maraviroc concentration
  # spiked into the incubation at time zero, so the dosing target is neither
  # `depot` nor `central`. Same convention as
  # HernandezLozano_2025_apramycin_invitro.R.
  dosing <- c("maraviroc")

  # See the sibling file Hyland_2008_maraviroc_hlm.R for the rationale: a
  # Supersome incubation has no body compartments, so each chemical species is
  # a state named after the species itself.
  paper_specific_compartments <- c("maraviroc", "uk408027")

  compartmentData <- list(
    maraviroc = list(
      analyte = "maraviroc",
      units = "uM",
      specimen = "administration site",
      verified = TRUE
    ),
    uk408027 = list(
      analyte = "UK-408,027 (the secondary amine formed by N-dealkylation of maraviroc)",
      units = "uM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (recombinant human CYP3A4 Supersomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Supersomes (BD Biosciences) expressing recombinant human CYP3A4, incubated in 50 mM phosphate buffer pH 7.4 with 1 mM MgCl2 and an NADPH-regenerating isocitric acid / isocitric acid dehydrogenase system; final incubation volume 100 uL",
    temperature = "37 C",
    kinetic_incubation = "10 pmol CYP3A4/mL for 15 min, within the ranges over which UK-408,027 formation was linear with time (up to 20 min) and with CYP content (up to 50 pmol/mL)",
    depletion_incubation = "1 uM maraviroc at 100 pmol recombinant CYP for up to 60 min; rates more stable than the 0.06 uL/pmol/min assay cut-off were reported as below that limit",
    concentration_range = "1 to 1000 uM maraviroc for the kinetic characterisation; 1 uM for substrate depletion and 50 uM for metabolite formation in the CYP-panel experiments",
    replication = "Km and Vmax are the mean of four determinations (Figure 4); the recombinant CYP panel results are the mean of quadruplicate determinations (Figure 3)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "At the 10 pmol/mL CYP concentration used for the kinetic study the authors calculate maraviroc to be 98% free, so the apparent Km is close to the real Km; the separately calculated fraction unbound in microsomes, 0.86, applies to the 100 pmol/mL (0.63 mg/mL protein) intrinsic-clearance incubation and is the value the paper divides by to place the intrinsic clearance on an unbound basis.",
      "Across a panel of Supersomes expressing CYP1A2, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP3A4 and CYP3A5, only rCYP3A4 depleted maraviroc at a rate significantly above the 0.06 uL/pmol/min assay cut-off. CYP2B6 formed a small but statistically significant amount of UK-408,027 at 0.1 pmol/pmol/min, 17-fold below the 1.7 pmol/pmol/min of rCYP3A4, and CYP3A5's intrinsic clearance was 25-fold below CYP3A4's; neither is carried in this model, because the paper concludes CYP3A4 is the only clinically relevant enzyme.",
      "The Km of 13 uM in recombinant CYP3A4 is consistent with the 21 uM measured in pooled human liver microsomes, and the two systems' unbound intrinsic clearances agree to 3% (47.4 versus 48.6 uL/min/mg microsomal protein) after the intersystem extrapolation factor and microsomal-binding corrections."
    )
  )

  ini({
    # --- Michaelis-Menten characterisation of UK-408,027 formation ---
    # Estimated by Grafit (version 4) from four determinations across 1-1000 uM.
    km_cyp3a4 <- 13
    label("Michaelis constant for maraviroc N-dealkylation to UK-408,027 in recombinant CYP3A4 (uM)") # Results, 'Kinetics of maraviroc N-dealkylation by rCYP3A4'; Abstract; Figure 4
    vmax_cyp3a4 <- 3
    label("Maximum velocity of UK-408,027 formation (pmol/pmol CYP3A4/min)") # Results, 'Kinetics of maraviroc N-dealkylation by rCYP3A4'; Abstract; Figure 4

    # --- Substrate-depletion intrinsic clearance ---
    clint <- 1.7
    label("Maraviroc depletion intrinsic clearance in recombinant CYP3A4, before the intersystem extrapolation factor (uL/min/pmol CYP3A4)") # Results, 'CLint estimates from HLM and rCYP'; Figure 3; Table 1 reports the same value as the Simcyp CLint input

    # --- Scaling constants used to bring the recombinant clearance onto a
    #     human liver microsome basis ---
    isef <- fixed(0.2)
    label("Intersystem extrapolation factor for CYP3A4 (unitless)") # Results, 'CLint estimates from HLM and rCYP'; value from Proctor 2004, reference 13
    cyp3a4_content <- fixed(120)
    label("CYP3A4 content of the human liver microsome batch (pmol/mg protein)") # Results, 'CLint estimates from HLM and rCYP'
    fumic <- fixed(0.86)
    label("Fraction of maraviroc unbound in the recombinant CYP3A4 incubation at 0.63 mg/mL protein (unitless)") # Results, 'CLint estimates from HLM and rCYP'; Table 1 'fu (mic)'

    # --- Assay design ---
    cyp_inc <- fixed(10)
    label("CYP3A4 concentration of the kinetic incubation (pmol/mL)") # Results, 'Kinetics of maraviroc N-dealkylation by rCYP3A4'

    # --- Residual error ---
    # The paper reports no fitted residual-error model; the value below is the
    # bioanalytical precision of the UK-408,027 LC-MS/MS assay, which is the
    # only variability the source quantifies for the fitted endpoint. See the
    # vignette's Assumptions and deviations section.
    propSd <- fixed(0.051)
    label("Proportional residual SD of the UK-408,027 formation velocity, from the assay CV at the 50 ng standard (fraction)") # Materials and methods, 'Metabolite formation assay' (CV 5.7, 5.1 and 4.2% at 5, 50 and 150 ng)
  })

  model({
    # 1. The paper's intrinsic-clearance scaling chain, reproduced so the
    #    derived values are visible to a downstream user: 1.7 * 0.2 = 0.34
    #    uL/min/pmol CYP3A4, 0.34 * 120 = 40.8 uL/min/mg microsomal protein,
    #    and 40.8 / 0.86 = 47.4 uL/min/mg on an unbound basis (Results,
    #    'CLint estimates from HLM and rCYP').
    clint_isef <- clint * isef
    clint_mg <- clint_isef * cyp3a4_content
    clint_u_mg <- clint_mg / fumic

    # 2. Published formation velocity. This is the quantity plotted against
    #    maraviroc concentration in Figure 4 and its Eadie-Hofstee inset, in
    #    pmol UK-408,027 formed per pmol CYP3A4 per min. No impurity term is
    #    carried: in the recombinant system the substrate contaminant was
    #    removed by comparison against control Supersomes rather than fitted,
    #    unlike the human liver microsome fit in the sibling model.
    vUK408027 <- vmax_cyp3a4 * maraviroc / (km_cyp3a4 + maraviroc)

    # 3. Route split, within the recombinant system and so before the
    #    intersystem extrapolation factor. The low-substrate limit of the
    #    Michaelis-Menten route is Vmax/Km = 3/13 = 0.231 uL/min/pmol CYP3A4,
    #    13.6% of the 1.7 uL/min/pmol depletion clearance; the remainder is
    #    carried as a first-order route, because substrate depletion was
    #    measured only at 1 uM, well below Km.
    clint_uk408027 <- vmax_cyp3a4 / km_cyp3a4
    clint_other <- clint - clint_uk408027

    # 4. Volumetric rates, in uM/min. vmax_cyp3a4 * cyp_inc has units
    #    pmol/(mL*min) = nmol/L/min, and clint * cyp_inc has units
    #    uL/(mL*min); both need the factor 1/1000 to reach uM/min and 1/min
    #    respectively.
    rate_uk408027 <- vmax_cyp3a4 * maraviroc / (km_cyp3a4 + maraviroc) *
      cyp_inc / 1000
    rate_other <- clint_other * maraviroc * cyp_inc / 1000

    d/dt(maraviroc) <- -(rate_uk408027 + rate_other)
    d/dt(uk408027) <- rate_uk408027

    # 5. Observation. The fitted endpoint is the formation velocity of
    #    Figure 4; the paper publishes no maraviroc concentration-time data, so
    #    the substrate concentration is returned as a derived quantity without
    #    a residual-error model. Its assay CV was 14% at 1 uM (Materials and
    #    methods, substrate-depletion assay).
    Cmaraviroc <- maraviroc

    vUK408027 ~ prop(propSd)
  })
}
