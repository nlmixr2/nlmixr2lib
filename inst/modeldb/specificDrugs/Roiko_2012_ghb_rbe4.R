Roiko_2012_ghb_rbe4 <- function() {
  description <- "In vitro (RBE4 rat brain capillary endothelial cell line). Michaelis-Menten model of the carrier-mediated uptake of the drug of abuse gamma-hydroxybutyric acid (GHB) across an in vitro model of the rat blood-brain barrier. The published fit is the single-transporter Michaelis-Menten velocity v = Vmax * C / (Km + C) (equation 2 of the source), where C is the GHB concentration in the uptake buffer and v the uptake rate normalised to cell protein; a Michaelis-Menten route plus a parallel diffusional clearance (equation 3) and a two-transporter model (equation 4) were both fitted and rejected on Akaike information criterion, coefficient of variation and residual plots, so no passive-diffusion term is carried here. Uptake was linear through 30 s, and the concentration-response was measured at a 15 s incubation, so the buffer concentration is held static and the model describes initial-rate conditions. GHB is a substrate of the proton-dependent monocarboxylate transporters MCT1, MCT2 and MCT4, of which MCT1 predominates in RBE4 cells; the fitted Km of 23.3 mM lies above the peak plasma GHB concentrations of 8.4 to 15.6 mM reached at the 400 to 800 mg/kg intravenous doses studied in the same paper, which is the quantitative basis for the paper's conclusion that GHB brain uptake is not capacity-limited over that dose range. Sibling model: Roiko_2012_ghb_hcmecd3, the same uptake reaction characterised in the human brain capillary endothelial cell line hCMEC/D3."
  reference <- "Roiko SA, Felmlee MA, Morris ME. Brain uptake of the drug of abuse gamma-hydroxybutyric acid in rats. Drug Metab Dispos. 2012 Jan;40(1):212-218. doi:10.1124/dmd.111.041749. PMID: 22031624. PMCID: PMC3250048. Michaelis-Menten equation and the two rejected alternatives: Materials and Methods, 'Data and Statistical Analysis', equations 2, 3 and 4. Km and Vmax estimates for RBE4 cells: Results, 'GHB Uptake in Brain Endothelial Cells', and the Abstract; the fitted curve is Figure 6B. Uptake time course and incubation design: Materials and Methods, 'GHB Cell Uptake Studies', and Figure 6A. Plasma and brain extracellular-fluid concentrations used as the in vivo context: Tables 1 and 2."
  vignette <- "Roiko_2012_ghb_brain_uptake"
  units <- list(time = "min", dosing = "mM (incubation concentration)", concentration = "mM")

  # The "dose" of this static in vitro system is the GHB concentration placed
  # in the uptake buffer at time zero, so the dosing target is neither `depot`
  # nor `central`. Same convention as HernandezLozano_2025_apramycin_invitro.R
  # and Hyland_2008_maraviroc_hlm.R.
  dosing <- c("ghb_buffer")

  # A cell-uptake incubation has no body compartments. Following the in vitro
  # convention of Hyland_2008_maraviroc_hlm.R (and, upstream of it,
  # HernandezLozano_2025_apramycin_invitro.R), the states are named after the
  # compartments of the incubation rather than the canonical body
  # compartments. `ghb_buffer` holds the GHB concentration of the uptake
  # buffer in mM; `ghb_cell` holds the GHB accumulated inside the cell
  # monolayer, normalised to cell protein, which is the quantity the
  # scintillation-counting assay returns.
  paper_specific_compartments <- c("ghb_buffer", "ghb_cell")

  compartmentData <- list(
    ghb_buffer = list(
      analyte = "gamma-hydroxybutyric acid (GHB)",
      units = "mM",
      specimen = "administration site",
      verified = TRUE
    ),
    ghb_cell = list(
      analyte = "gamma-hydroxybutyric acid (GHB)",
      units = "pmol/mg cell protein",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (RBE4 rat brain capillary endothelial cell line)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Immortalised rat brain capillary endothelial cell line RBE4, passages 39 to 44, provided by Prof. P. Couraud (University Rene Descartes, Paris)",
    system = "Confluent monolayers plated on individual type-I rat-tail-collagen-coated 35 mm wells, grown in 1:1 alpha-minimum essential medium / Ham's F10 with L-glutamine 2.0 mM, geneticin 300 ug/mL, human recombinant fibroblast growth factor 1 ng/mL, gentamicin 50 ug/mL and 10 percent qualified fetal bovine serum",
    medium = "Uptake buffer containing NaCl 138 mM, CaCl2 1.8 mM, KCl 5.4 mM, MgSO4 0.8 mM, Na2HPO4 1.0 mM, D-glucose 5.5 mM and HEPES 20 mM, pH 7.4",
    temperature = "Cells equilibrated 30 min at 37 C, then equilibrated to room temperature for 5 min; all uptake incubations were run at room temperature",
    incubation_time = "15 s for the concentration-dependence experiment, chosen to minimise loss due to metabolism and loss of the radiolabel; the separate time-course experiment used 0.25, 0.5, 1, 2, 5 and 10 min with 58 nM [3H]GHB and showed rapid linear uptake through 30 s",
    concentration_range = "Materials and Methods lists 0.01, 0.1, 1, 3, 5, 10, 30 and 50 mM GHB and the Results describe the fitted range as 100 uM to 50 mM; Figure 6B additionally plots a point near 20 mM that the Methods list omits, so the RBE4 fit spans roughly 0 to 50 mM",
    replication = "Mean +/- S.D. of three experiments, each performed in triplicate (Figure 6)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "The companion in vivo arm of the same paper dosed male Sprague-Dawley rats (280 to 320 g) with GHB 400, 600 or 800 mg/kg intravenously, n = 4 per dose, and sampled plasma plus frontal-cortex extracellular fluid by microdialysis for 6 h. That arm was analysed noncompartmentally in Phoenix WinNonlin 6.0 and yields no compartmental model, so it is not extracted; its numbers are reproduced here because they are the in vivo context in which this uptake model's Km is interpreted.",
      "Peak plasma concentrations were 874 +/- 114, 1220 +/- 142 and 1620 +/- 113 ug/mL (8.4, 11.7 and 15.6 mM) at 5 min, and peak brain extracellular-fluid concentrations 41.3 +/- 14, 57.8 +/- 30 and 69.2 +/- 19 ug/mL (0.39, 0.55 and 0.66 mM) at 50 min, for the 400, 600 and 800 mg/kg doses respectively (Table 1).",
      "Unbound brain extracellular-fluid to plasma partition coefficients were 0.079 +/- 0.01, 0.070 +/- 0.02 and 0.070 +/- 0.02 across the same doses and did not differ significantly (Table 2); GHB has negligible plasma protein binding, so plasma concentrations are unbound concentrations.",
      "MCT involvement was confirmed but not quantified as a model term: uptake of 10 mM GHB at pH 6.5 was 126 percent of uptake at pH 7.4 in RBE4 cells, and the MCT inhibitor alpha-cyano-4-hydroxycinnamate at 2.5 mM reduced uptake of 10 mM GHB to 60 percent of control in RBE4 cells and 66 percent in hCMEC/D3 cells. Each of these is a single-concentration observation, which cannot distinguish an effect on Vmax from an effect on Km, so no pH or inhibitor term is carried in the model; see the vignette's Assumptions and deviations section.",
      "Kinetic parameters were estimated by weighted nonlinear regression in Phoenix WinNonlin 6.0. The weighting scheme is not stated and no residual-error model is reported."
    )
  )

  ini({
    # --- Michaelis-Menten characterisation of GHB uptake (equation 2) -----
    # Both parameters were estimated by weighted nonlinear regression; the
    # reported +/- values are the standard deviations of the estimates, not
    # inter-experiment variability, so they are recorded in the labels rather
    # than encoded as random effects.
    lkm_uptake <- log(23.3)
    label("Log Michaelis constant for GHB uptake into RBE4 cells (mM)")                 # Results, 'GHB Uptake in Brain Endothelial Cells': Km = 23.3 +/- 5 mM for RBE4 cells; also the Abstract
    lvmax_uptake <- log(258)
    label("Log maximum velocity of GHB uptake into RBE4 cells (pmol/mg protein/min)")   # Results, 'GHB Uptake in Brain Endothelial Cells': Vmax = 258 +/- 41 pmol/mg/min for RBE4 cells; also the Abstract

    # --- Residual error ---
    # The source reports no residual-error model and no assay coefficient of
    # variation for the [3H]GHB scintillation-counting endpoint; the only
    # variability it quantifies for the fitted endpoint is the S.D. of three
    # experiments plotted as error bars in Figure 6B, whose values are not
    # printed. Held at zero rather than invented, so the file remains a
    # typical-value model; see the vignette's Assumptions and deviations
    # section.
    propSd <- fixed(0)
    label("Proportional residual SD of the GHB uptake velocity (fraction); magnitude not reported")  # Materials and Methods, 'Data and Statistical Analysis': "weighted nonlinear regression analysis" with no weighting scheme or residual model stated
  })

  model({
    # 1. Back-transform. Km is in mM, matching the units of the buffer state;
    #    Vmax is in pmol per mg cell protein per min, matching the units of
    #    the cell state divided by the time unit.
    km_uptake <- exp(lkm_uptake)
    vmax_uptake <- exp(lvmax_uptake)

    # 2. Equation 2 of the source, the retained single-transporter
    #    Michaelis-Menten velocity. This is the quantity plotted against GHB
    #    concentration in Figure 6B.
    #
    #    Equation 3, v = Vmax * C / (Km + C) + P * C, adds a parallel
    #    diffusional clearance P, and equation 4,
    #    v = Vmax1 * C / (Km1 + C) + Vmax2 * C / (Km2 + C), adds a second
    #    transporter. Both were fitted and rejected: "Equation 2 yielded the
    #    smallest coefficient of variation percentage and AIC value and was
    #    subsequently used to estimate kinetic uptake parameters". Neither P
    #    nor the second transporter's constants are reported, so no term for
    #    them is carried.
    #
    #    (Note on the source PDF: its symbol font substitutes the glyphs for
    #    "times" and "plus" in the displayed equations, so equation 2 extracts
    #    as "v = (Vmax + C) / (Km x C)". The intended form is the standard
    #    Michaelis-Menten expression below, which is what the Figure 6B and
    #    6D fits and the text description of "a single-transporter model"
    #    require.)
    vGhbUptake <- vmax_uptake * ghb_buffer / (km_uptake + ghb_buffer)

    # 3. Static exposure. The concentration-dependence experiment used a 15 s
    #    incubation within the window over which uptake was linear in time,
    #    i.e. initial-rate conditions in which buffer depletion is negligible.
    #    Neither the buffer volume nor the cell protein per well is reported,
    #    so a depletion term could not be computed even if it were wanted.
    d/dt(ghb_buffer) <- 0

    # 4. Intracellular accumulation, in pmol per mg cell protein. This is what
    #    the scintillation-counting assay measures; the paper converts it to
    #    the velocity above by dividing by the incubation time.
    d/dt(ghb_cell) <- vGhbUptake

    # 5. Observation. The fitted endpoint is the uptake velocity of Figure 6B.
    #    It is not named `Cc` because it is a reaction velocity in
    #    pmol/mg/min, not a drug concentration; see the vignette's Assumptions
    #    and deviations section. Same decision as Hyland_2008_maraviroc_hlm.R.
    vGhbUptake ~ prop(propSd)
  })
}
