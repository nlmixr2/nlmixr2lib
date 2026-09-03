Hyland_2008_maraviroc_hlm <- function() {
  description <- "In vitro (pooled human liver microsomes, 60 donors). Michaelis-Menten enzyme-kinetic model of the CYP3A4-mediated N-dealkylation of maraviroc to its secondary-amine metabolite UK-408,027, plus the parallel oxidative routes that make up the remainder of maraviroc's microsomal intrinsic clearance. The published fit is a standard Michaelis-Menten velocity with an additive impurity term, v = Vmax * [S] / (Km + [S]) + C * [S], where C absorbs a concentration-dependent UK-408,027 contaminant present in the maraviroc substrate; the authors report that fitting impurity-corrected data to a plain Michaelis-Menten equation gave similar Km and Vmax, and C itself is not published, so it is carried here fixed at zero. The UK-408,027 route accounts for Vmax/Km = 0.0214 uL/min/pmol CYP, which is 20% of the 0.106 uL/min/pmol total depletion intrinsic clearance measured by substrate loss; the balance is assigned to the other CYP3A4-mediated oxidative pathways as a first-order route, because substrate depletion was measured only at 1 uM, far below Km, and so characterises a linear clearance. Sibling model: Hyland_2008_maraviroc_rcyp3a4, the same reaction characterised in recombinant CYP3A4 Supersomes."
  reference <- "Hyland R, Dickins M, Collins C, Jones H, Jones B. Maraviroc: in vitro assessment of drug-drug interaction potential. Br J Clin Pharmacol. 2008 Oct;66(4):498-507. doi:10.1111/j.1365-2125.2008.03198.x. PMID: 18492127. PMCID: PMC2561101. Michaelis-Menten equation with the impurity term and the Km / Vmax estimates: Results, 'Kinetics of maraviroc N-dealkylation in human liver microsomes', and Figure 2. Substrate-depletion intrinsic clearance, microsomal CYP content and fraction unbound in microsomes: Results, 'CLint estimates from HLM and rCYP'. Incubation design and bioanalytical precision: Materials and methods, 'Assays for maraviroc metabolism'. The statement that UK-408,027 formation is approximately 20% of the depletion intrinsic clearance is in the Discussion, second paragraph."
  vignette <- "Hyland_2008_maraviroc_invitro"
  units <- list(time = "min", dosing = "uM (incubation concentration)", concentration = "uM")

  # The "dose" of this static in vitro system is the maraviroc concentration
  # spiked into the incubation at time zero, so the dosing target is neither
  # `depot` nor `central`. Same convention as
  # HernandezLozano_2025_apramycin_invitro.R.
  dosing <- c("maraviroc")

  # A microsomal incubation has no body compartments. Following the in vitro
  # convention of HernandezLozano_2025_apramycin_invitro.R (and, upstream of
  # it, Nielsen_2007_semimechanistic_antibiotic_pd.R), each chemical species
  # in the incubation is a state named after the species itself. `uk408027`
  # is the secondary-amine N-dealkylation product the paper calls UK-408,027;
  # it is not a registered metabolite suffix, so the parent + metabolite
  # `<canonical>_<metab>` scheme in compartment-names.md does not apply.
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
    species = "in vitro (pooled human liver microsomes, 60 donors)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Human liver microsomes prepared from a pool of 60 donors (BD Biosciences), incubated in 50 mM phosphate buffer pH 7.4 with 1 mM MgCl2 and an NADPH-regenerating isocitric acid / isocitric acid dehydrogenase system; final incubation volume 100 uL",
    temperature = "37 C",
    kinetic_incubation = "0.1 mg/mL microsomal protein for 15 min, within the ranges over which UK-408,027 formation was linear with time (up to 20 min) and with protein (up to 0.5 mg/mL)",
    depletion_incubation = "1 uM maraviroc at 0.5 uM total CYP for up to 60 min; the first-order rate constant was taken as the gradient of ln(maraviroc / midazolam peak-area ratio) against incubation time",
    concentration_range = "1 to 1000 uM maraviroc for the kinetic characterisation; 1 uM for substrate depletion and 50 uM for metabolite formation in the chemical-inhibition experiments",
    replication = "Km and Vmax are the mean of seven determinations (Figure 2); the chemical-inhibition results are the mean of quadruplicate determinations (Table 2)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "At the 0.1 mg/mL protein concentration used for the kinetic study the authors calculate maraviroc to be 98% free, so the apparent Km is close to the real Km; the separately measured fraction unbound in microsomes, 0.72, was determined at the 1.5 mg/mL protein concentration of the intrinsic-clearance assay and is the value the paper divides by to place the intrinsic clearance on an unbound basis.",
      "Chemical inhibition (Table 2) localised the reaction to CYP3A4: 1 uM ketoconazole left 16.1 +/- 5.7% of control UK-408,027 formation and 17.3 +/- 3.4% of control maraviroc depletion (both P < 0.00001), while furafylline, sulphaphenazole, benzylnirvanol and quinidine had no significant effect. Uninhibited control activities were 75.9 pmol/mg/min for UK-408,027 formation at 50 uM maraviroc and 0.106 uL/pmol/min for maraviroc depletion at 1 uM.",
      "Maraviroc itself was a weak CYP inhibitor, with IC50 > 30 uM against probe substrates for CYP1A2, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6 and CYP3A4; no IC50 or Ki point estimate is reported, so no inhibition term is carried in this model."
    )
  )

  ini({
    # --- Michaelis-Menten characterisation of UK-408,027 formation ---
    # Estimated by Grafit (version 4) from seven determinations across 1-1000 uM.
    km_cyp3a4 <- 21
    label("Michaelis constant for maraviroc N-dealkylation to UK-408,027 (uM)") # Results, 'Kinetics of maraviroc N-dealkylation in human liver microsomes'; Abstract; Figure 2
    vmax_cyp3a4 <- 0.45
    label("Maximum velocity of UK-408,027 formation (pmol/pmol total CYP/min)") # Results, 'Kinetics of maraviroc N-dealkylation in human liver microsomes'; Abstract; Figure 2

    # The paper's fitted equation carries an additive impurity term C * [S] that
    # absorbs a concentration-dependent UK-408,027 contaminant in the substrate.
    # C is named in the equation but its value is not published, and the authors
    # state that fitting impurity-corrected data to a plain Michaelis-Menten
    # equation gave a similar Km and Vmax. It is therefore carried here at zero
    # so the published equation's structure is preserved and a user who obtains
    # C can set it; see the vignette's Assumptions and deviations section.
    c_impurity <- fixed(0)
    label("Impurity proportionality constant of the substrate contaminant, value not published (pmol/pmol total CYP/min per uM)") # Results, 'Kinetics of maraviroc N-dealkylation in human liver microsomes' (equation only)

    # --- Substrate-depletion intrinsic clearance ---
    clint <- 0.106
    label("Maraviroc depletion intrinsic clearance, not corrected for microsomal binding (uL/min/pmol total CYP)") # Results, 'CLint estimates from HLM and rCYP'; Table 1 ('CLint (rhCYP3A4)' row reports the recombinant value; this is the HLM value quoted in the same paragraph) and the Table 2 footnote

    # --- Measured properties of the microsomal batch and of the incubation ---
    cyp_total <- fixed(330)
    label("Total CYP content of the human liver microsome pool (pmol/mg protein)") # Results, 'CLint estimates from HLM and rCYP'
    fumic <- fixed(0.72)
    label("Fraction of maraviroc unbound in human liver microsomes, measured at 1.5 mg/mL protein (unitless)") # Results, 'CLint estimates from HLM and rCYP'
    prot_inc <- fixed(0.1)
    label("Microsomal protein concentration of the kinetic incubation (mg/mL)") # Results, 'Kinetics of maraviroc N-dealkylation in human liver microsomes'

    # --- Residual error ---
    # The paper reports no fitted residual-error model; the value below is the
    # bioanalytical precision of the UK-408,027 LC-MS/MS assay, which is the
    # only variability the source quantifies for the fitted endpoint. See the
    # vignette's Assumptions and deviations section.
    propSd <- fixed(0.051)
    label("Proportional residual SD of the UK-408,027 formation velocity, from the assay CV at the 50 ng standard (fraction)") # Materials and methods, 'Metabolite formation assay' (CV 5.7, 5.1 and 4.2% at 5, 50 and 150 ng)
  })

  model({
    # 1. Assay-design constants and unit bridges.
    #    Total CYP concentration of the kinetic incubation, in pmol CYP/mL.
    cyp_inc <- prot_inc * cyp_total

    #    The paper's intrinsic-clearance arithmetic, reproduced so the derived
    #    values are visible to a downstream user: 0.106 * 330 = 35.0 uL/min/mg
    #    microsomal protein, and 35.0 / 0.72 = 48.6 uL/min/mg on an unbound
    #    basis (Results, 'CLint estimates from HLM and rCYP').
    clint_mg <- clint * cyp_total
    clint_u_mg <- clint_mg / fumic

    # 2. Published formation velocity. This is the quantity plotted against
    #    maraviroc concentration in Figure 2 and its Eadie-Hofstee inset, in
    #    pmol UK-408,027 formed per pmol total CYP per min.
    vUK408027 <- vmax_cyp3a4 * maraviroc / (km_cyp3a4 + maraviroc) +
      c_impurity * maraviroc

    # 3. Route split. The low-substrate limit of the Michaelis-Menten route is
    #    Vmax/Km = 0.45/21 = 0.0214 uL/min/pmol CYP, which is 20.2% of the
    #    0.106 uL/min/pmol total depletion clearance -- the "approximately 20%"
    #    of the Discussion. The remainder is carried as a first-order route,
    #    because substrate depletion was measured only at 1 uM, well below Km,
    #    and so characterises a linear clearance rather than a saturable one.
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
    #    Figure 2; the paper publishes no maraviroc concentration-time data, so
    #    the substrate concentration is returned as a derived quantity without
    #    a residual-error model. Its assay CV was 14% at 1 uM (Materials and
    #    methods, substrate-depletion assay).
    Cmaraviroc <- maraviroc

    vUK408027 ~ prop(propSd)
  })
}
