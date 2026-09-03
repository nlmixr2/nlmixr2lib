Roiko_2012_ghb_hcmecd3 <- function() {
  description <- "In vitro (hCMEC/D3 human brain capillary endothelial cell line). Michaelis-Menten model of the carrier-mediated uptake of the drug of abuse gamma-hydroxybutyric acid (GHB) across an in vitro model of the human blood-brain barrier. The published fit is the single-transporter Michaelis-Menten velocity v = Vmax * C / (Km + C) (equation 2 of the source), where C is the GHB concentration in the uptake buffer and v the uptake rate normalised to cell protein; a Michaelis-Menten route plus a parallel diffusional clearance (equation 3) and a two-transporter model (equation 4) were both fitted and rejected on Akaike information criterion, coefficient of variation and residual plots, so no passive-diffusion term is carried here. Uptake was linear through 30 s, and the concentration-response was measured at a 15 s incubation, so the buffer concentration is held static and the model describes initial-rate conditions. GHB is a substrate of the proton-dependent monocarboxylate transporters MCT1, MCT2 and MCT4, of which MCT1 and MCT4 messenger RNA have been detected in hCMEC/D3 cells; the human Km of 18.1 mM is slightly lower and the human Vmax slightly lower than the rat estimates, so the two cell lines describe near-identical uptake kinetics. Sibling model: Roiko_2012_ghb_rbe4, the same uptake reaction characterised in the rat brain capillary endothelial cell line RBE4."
  reference <- "Roiko SA, Felmlee MA, Morris ME. Brain uptake of the drug of abuse gamma-hydroxybutyric acid in rats. Drug Metab Dispos. 2012 Jan;40(1):212-218. doi:10.1124/dmd.111.041749. PMID: 22031624. PMCID: PMC3250048. Michaelis-Menten equation and the two rejected alternatives: Materials and Methods, 'Data and Statistical Analysis', equations 2, 3 and 4. Km and Vmax estimates for hCMEC/D3 cells: Results, 'GHB Uptake in Brain Endothelial Cells', and the Abstract; the fitted curve is Figure 6D. Uptake time course and incubation design: Materials and Methods, 'GHB Cell Uptake Studies', and Figure 6C."
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
    species = "in vitro (hCMEC/D3 human brain capillary endothelial cell line)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Immortalised human brain capillary endothelial cell line hCMEC/D3, passages 28 to 33, provided by Prof. P. Couraud (University Rene Descartes, Paris)",
    system = "Confluent monolayers plated on individual type-I rat-tail-collagen-coated 35 mm wells, grown in EBM-2 medium supplemented with 2 percent fetal bovine serum and the EGM-2 bullet-kit growth factors",
    medium = "Uptake buffer containing NaCl 138 mM, CaCl2 1.8 mM, KCl 5.4 mM, MgSO4 0.8 mM, Na2HPO4 1.0 mM, D-glucose 5.5 mM and HEPES 20 mM, pH 7.4",
    temperature = "Cells equilibrated 30 min at 37 C, then equilibrated to room temperature for 5 min; all uptake incubations were run at room temperature",
    incubation_time = "15 s for the concentration-dependence experiment, chosen to minimise loss due to metabolism and loss of the radiolabel; the separate time-course experiment used 0.25, 0.5, 1, 2, 5 and 10 min with 58 nM [3H]GHB and showed rapid linear uptake through 30 s",
    concentration_range = "Materials and Methods lists 0.01, 0.1, 1, 3, 5, 10, 30 and 50 mM GHB for both cell lines and the Results describe the fitted range as 100 uM to 50 mM, but the hCMEC/D3 panel (Figure 6D) has an x axis that stops at 30 mM and plots no point above it, so this fit is supported only to about 30 mM; extrapolation above that concentration rests on the Michaelis-Menten form rather than on data",
    replication = "Mean +/- S.D. of three experiments, each performed in triplicate (Figure 6)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "This human cell line was studied alongside the rat RBE4 line in a paper whose in vivo arm is entirely rat; the human arm exists to test whether the rat blood-brain-barrier model transfers, and the authors conclude that the similarity of the in vitro uptake characteristics indicates both serve as useful systems to study GHB brain uptake.",
      "The estimated Km values of 18.1 mM (human) and 23.3 mM (rat) are higher than previously reported values for GHB uptake, which the authors attribute to running these experiments at pH 7.4 to mimic physiological conditions at the blood-brain barrier rather than at the lower pH of 6 or 6.5 used previously.",
      "MCT involvement was confirmed but not quantified as a model term: the MCT inhibitor alpha-cyano-4-hydroxycinnamate at 2.5 mM reduced uptake of 10 mM GHB to 66 percent of control in hCMEC/D3 cells and to 60 percent in RBE4 cells, and uptake at pH 6.5 was 126 percent of that at pH 7.4 in RBE4 cells. Each of these is a single-concentration observation, which cannot distinguish an effect on Vmax from an effect on Km, so no pH or inhibitor term is carried in the model; see the vignette's Assumptions and deviations section.",
      "MCT1 and MCT4 messenger RNA have both been detected in hCMEC/D3 cells; MCT1 is the predominant MCT in RBE4 cells, with MCT2 and MCT4 detected to a lesser degree.",
      "Kinetic parameters were estimated by weighted nonlinear regression in Phoenix WinNonlin 6.0. The weighting scheme is not stated and no residual-error model is reported."
    )
  )

  ini({
    # --- Michaelis-Menten characterisation of GHB uptake (equation 2) -----
    # Both parameters were estimated by weighted nonlinear regression; the
    # reported +/- values are the standard deviations of the estimates, not
    # inter-experiment variability, so they are recorded in the labels rather
    # than encoded as random effects.
    lkm_uptake <- log(18.1)
    label("Log Michaelis constant for GHB uptake into hCMEC/D3 cells (mM)")                 # Results, 'GHB Uptake in Brain Endothelial Cells': Km = 18.1 +/- 3 mM for hCMEC/D3 cells; also the Abstract
    lvmax_uptake <- log(248)
    label("Log maximum velocity of GHB uptake into hCMEC/D3 cells (pmol/mg protein/min)")   # Results, 'GHB Uptake in Brain Endothelial Cells': Vmax = 248 +/- 34 pmol/mg/min for hCMEC/D3 cells; also the Abstract

    # --- Residual error ---
    # The source reports no residual-error model and no assay coefficient of
    # variation for the [3H]GHB scintillation-counting endpoint; the only
    # variability it quantifies for the fitted endpoint is the S.D. of three
    # experiments plotted as error bars in Figure 6D, whose values are not
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
    #    concentration in Figure 6D.
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

    # 5. Observation. The fitted endpoint is the uptake velocity of Figure 6D.
    #    It is not named `Cc` because it is a reaction velocity in
    #    pmol/mg/min, not a drug concentration; see the vignette's Assumptions
    #    and deviations section. Same decision as Hyland_2008_maraviroc_hlm.R.
    vGhbUptake ~ prop(propSd)
  })
}
