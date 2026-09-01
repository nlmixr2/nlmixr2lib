Mukker_2026_tuvusertib_mouse_tgi <- function() {
  description <- paste(
    "Preclinical (mouse bearing a CTG-3021 ARID1A-mutated gastric-cancer",
    "patient-derived xenograft). Semi-mechanistic PK-efficacy model of the ATR",
    "inhibitor tuvusertib: a two-compartment linear plasma PK model with",
    "first-order absorption drives a standard Simeoni 2004 tumor-growth-inhibition",
    "model in which unperturbed growth is exponential then linear and drug-damaged",
    "cells pass through a three-stage transit chain before dying. Fitted to",
    "25-50 mg/kg continuous daily and 2 weeks on / 1 week off dosing."
  )
  reference <- paste(
    "Mukker JK, Diderichsen PM, Hellmann F, Yap TA, Plummer R, Tolcher AW,",
    "de Bono JS, Gounaris I, Szucs Z, Zimmermann A, Kareva I, Bolleddula J,",
    "Seithel-Keuth A, Locatelli G, Enderlin M, Hicking C, Zutshi A, Gao W,",
    "Strotmann R, Benincosa L, Venkatakrishnan K.",
    "Integrated Population Pharmacokinetic, Pharmacodynamic, and Safety Analyses",
    "to Inform Dosage Selection in the Clinical Development of the ATR Inhibitor",
    "Tuvusertib. Clin Pharmacol Ther. 2026;119(3):618-628. doi:10.1002/cpt.70029.",
    "Parameter values from Supplementary Table S1. Structural TGI formulation from",
    "Simeoni M, Magni P, Cammia C, De Nicolao G, Croci V, Pesenti E, Germani M,",
    "Poggesi I, Rocchetti M. Predictive pharmacokinetic-pharmacodynamic modeling of",
    "tumor growth kinetics in xenograft models after administration of anticancer",
    "agents. Cancer Res. 2004;64(3):1094-1101. doi:10.1158/0008-5472.CAN-03-2524.",
    sep = " "
  )
  vignette <- "Mukker_2026_tuvusertib"

  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot          = list(analyte = "tuvusertib", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "tuvusertib", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1    = list(analyte = "tuvusertib", units = "mg/kg", specimen = "tissue", verified = TRUE),
    cycling_cells  = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "mouse (CTG-3021 ARID1A-mutated gastric-cancer patient-derived xenograft)",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    disease_state = "subcutaneous CTG-3021 ARID1A-mutated gastric-cancer patient-derived xenograft; ARID1A loss causes topoisomerase 2A and cell-cycle defects that make the cells more reliant on ATR activity",
    dose_range    = "25-50 mg/kg tuvusertib, continuous daily and 2 weeks on / 1 week off intermittent regimens",
    regions       = "preclinical in-vivo xenograft study",
    notes         = paste(
      "Mukker 2026 does not report the number of animals per arm, the",
      "randomisation tumor volume, or any interindividual / residual variability",
      "for this fit; Supplementary Table S1 lists point estimates only. The same",
      "fitted TGI system was subsequently driven by HUMAN tuvusertib plasma",
      "concentrations to produce the translational multicycle simulations in",
      "Figure 4c."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Mouse plasma PK -- Mukker 2026 Supplementary Table S1.
    #
    # UNIT NOTE (see the vignette's Assumptions and deviations section):
    # Table S1 labels the volumes "mL/kg" and the clearances "mL/h/kg". Taken
    # literally, V1 = 9.834 mL/kg is below mouse plasma volume (~40-50 mL/kg),
    # which is impossible for an APPARENT (oral, V/F >= V) volume. The
    # micro-constants Table S1 also prints are self-consistent with the ratios
    # (Cl1/V1 = 0.734 vs the tabulated k10 = 0.74; Cl2/V1 = 0.104 vs k12 = 0.11;
    # Cl2/V2 = 0.100 vs k21 = 0.1), so the NUMBERS are right and only the
    # metric prefix is wrong. Reading them as L/kg and L/h/kg puts the
    # concentration on the ng/mL scale that the Simeoni kill constant k2
    # requires: the Simeoni threshold concentration lambda0/k2 = 0.0022/0.00002
    # = 110 ng/mL then sits just below the human average steady-state
    # concentrations at the active clinical doses, which is the only reading
    # under which the paper's own Figure 4c translational simulation (human
    # ng/mL concentrations driving this TGI system) produces any tumor effect
    # at all. Encoded here as L/kg accordingly.
    lvc  <- fixed(log(9.834));  label("Apparent central volume of distribution (V1, L/kg; Table S1 prints mL/kg -- see file note)")     # Table S1: V1 = 9.834
    lvp  <- fixed(log(10.202)); label("Apparent peripheral volume of distribution (V2, L/kg; Table S1 prints mL/kg -- see file note)")  # Table S1: V2 = 10.202
    lcl  <- fixed(log(7.216));  label("Apparent clearance from the central compartment (Cl1, L/h/kg; Table S1 prints mL/h/kg)")         # Table S1: Cl1 = 7.216
    lq   <- fixed(log(1.02));   label("Apparent intercompartmental clearance (Cl2, L/h/kg; Table S1 prints mL/h/kg)")                   # Table S1: Cl2 = 1.02
    lka  <- fixed(log(1.3));    label("Absorption rate constant (k01, 1/h)")                                                            # Table S1: k01 = 1.3 1/h

    # ---------------------------------------------------------------------
    # Simeoni tumor-growth-inhibition parameters -- Table S1.
    # ---------------------------------------------------------------------
    ltumorExpGrowth  <- fixed(log(0.0022));  label("Tumor growth rate during the initial exponential phase (lambda0, 1/h)")            # Table S1: lambda0 = 0.0022 1/h
    ltumorLinGrowth  <- fixed(log(6.4));     label("Tumor growth rate during the later linear phase (lambda1, mm^3/h)")                # Table S1: lambda1 = 6.4 Vol/h
    ldamageTransit   <- fixed(log(0.0003));  label("Rate at which damaged cells progress through the death chain (k1, 1/h)")           # Table S1: k1 = 0.0003 1/h
    ldrugSlope       <- fixed(log(0.00002)); label("Rate of drug-induced irreversible cell damage (k2, 1/(h*ng/mL))")                  # Table S1: k2 = 0.00002 1/h/conc.
    psi              <- fixed(20);           label("Switch sharpness between the exponential and linear growth phases (unitless)")     # Table S1: psi = 20

    # Randomisation tumor volume. NOT reported anywhere in Mukker 2026 or its
    # supplement; assumed from the Figure 4c y-axis, whose curves all originate
    # just under 300 mm^3. Override via the data or ini() when a measured
    # starting volume is available.
    lrbase_tumor <- fixed(log(300)); label("Tumor volume at the start of treatment (mm^3; assumed, not reported in the source)")

    # Residual error is NOT reported for this fit -- Table S1 gives point
    # estimates only, with no residual-error or interindividual-variability
    # terms. Both SDs are therefore carried as fixed(0) rather than filled in
    # with a borrowed magnitude: the model simulates the published typical-value
    # tumor-growth trajectory exactly, and a user who wants a stochastic
    # xenograft cohort must supply their own residual magnitude. Flagged in the
    # vignette's Assumptions and deviations section.
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on tumor volume (fraction; ZERO -- not reported in the source)")
    addSd_tumor_vol  <- fixed(0); label("Additive residual error on tumor volume (mm^3; ZERO -- not reported in the source)")
  })

  model({
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    cl  <- exp(lcl)
    q   <- exp(lq)
    ka  <- exp(lka)

    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    damageTransit  <- exp(ldamageTransit)
    drugSlope      <- exp(ldrugSlope)
    rbase_tumor    <- exp(lrbase_tumor)

    # --- Mouse plasma PK ------------------------------------------------------
    # Amounts are mg/kg and volumes L/kg, so central/vc is mg/L; x1000 gives the
    # ng/mL scale that drugSlope (k2) is expressed on.
    Cc <- 1000 * central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # --- Simeoni 2004 tumor growth inhibition ---------------------------------
    # Unperturbed growth is first-order at lambda0 while the tumor is small and
    # zero-order at lambda1 once it is large; psi = 20 makes the switch abrupt.
    # The cytotoxic effect moves cells out of the cycling pool at rate k2*Cc and
    # they then transit three damage stages at k1 before leaving the system.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    drugEffectCyclingCells <- drugSlope * Cc

    d/dt(cycling_cells)  <- tumorExpGrowth * cycling_cells /
                             (1 + (tumorExpGrowth / tumorLinGrowth * tumor_vol)^psi)^(1 / psi) -
                            drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    cycling_cells(0)  <- rbase_tumor
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    tumor_vol ~ prop(propSd_tumor_vol) + add(addSd_tumor_vol)
  })
}
