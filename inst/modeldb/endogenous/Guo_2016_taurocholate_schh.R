Guo_2016_taurocholate_schh <- function() {
  description <- "In vitro (sandwich-cultured human hepatocytes, SCHH). Mechanistic hepatobiliary disposition model of the model bile acid taurocholate (as its deuterated tracer d8-TCA) in SCHH incubated with a physiologic 4 percent BSA medium, together with the transporter-inhibition overlay the paper uses to predict drug-bile acid interactions. Three states hold total TCA amount: buffer (incubation medium), cells (hepatocytes) and bile (canalicular networks). Uptake from buffer into cells is a single linear clearance CLUptake; efflux out of cells splits into a basolateral clearance CLBL back to the buffer and a biliary clearance CLBile into the canalicular networks, from which TCA returns to the buffer at first-order rate KFlux as the canaliculi contract. The binary covariate BUFFER_CA_FREE switches between the paper's two experimental configurations: with standard Ca2+-containing HBSS the tight junctions are intact and bile is a separate accumulating state (observed as Cells + Bile), whereas with Ca2+/Mg2+-free HBSS the junctions are disrupted so the biliary flux discharges straight into the buffer and only Cells is observed. Inhibition of each transport pathway enters as a competitive term CL / (1 + [I] / IC50): uptake is driven by the unbound inhibitor concentration in the medium and is split 70 percent NTCP / 30 percent OATP, whereas biliary (BSEP) and basolateral (MRP3) efflux are driven by the intracellular inhibitor concentration. Inhibitory potencies for the paper's two model inhibitors, telmisartan and bosentan, are carried as fixed ini() entries and selected by CONMED_TELMISARTAN / CONMED_BOSENTAN, so a user can re-point the model at a different inhibitor by overriding them. The model is deterministic: the paper propagated parameter uncertainty by Monte Carlo sampling of the Table 2 clearances rather than by a hierarchical random-effects structure."
  reference <- "Guo C, Yang K, Brouwer KR, St Claire RL 3rd, Brouwer KLR. Prediction of Altered Bile Acid Disposition Due to Inhibition of Multiple Transporters: An Integrated Approach Using Sandwich-Cultured Hepatocytes, Mechanistic Modeling, and Simulation. J Pharmacol Exp Ther. 2016 Aug;358(2):324-333. doi:10.1124/jpet.116.231928. PMID: 27233294. PMCID: PMC4959093. Model structure and differential equations: Materials and Methods, 'Determination of Kinetic Parameters for d8-TCA Using Mechanistic Modeling', eqs. 1-5, plus the model schemes in Figure 1A. Clearance and KFlux estimates: Table 2. Inhibition equations: Materials and Methods, 'Simulation of Inhibitor Effects on TCA Disposition and Comparison with Experimental Results', eqs. 7, 8 and 10. Inhibitory potencies: Table 1. Measured inhibitor concentrations: Table 3. Predicted-versus-observed fold changes: Table 4."
  vignette <- "Guo_2016_taurocholate_schh"
  units <- list(time = "min", dosing = "pmol", concentration = "umol/L")

  # The dose is the d8-TCA placed into the incubation medium at time 0, so the
  # dosing target is the buffer state rather than depot or central.
  dosing <- c("buffer")

  # Paper-mechanistic in vitro states: the three matrices of the SCHH assay.
  # These are not instances of any canonical PK compartment role.
  paper_specific_compartments <- c("buffer", "cells", "bile")

  compartmentData <- list(
    buffer = list(analyte = "taurocholate (d8-TCA)", units = "pmol", specimen = "administration site", verified = TRUE),
    cells  = list(analyte = "taurocholate (d8-TCA)", units = "pmol", specimen = "tissue", verified = TRUE),
    bile   = list(analyte = "taurocholate (d8-TCA)", units = "pmol", specimen = "bile", verified = TRUE)
  )

  covariateData <- list(
    PROTEIN_MG = list(
      description = "Cellular protein content of the incubated well, used to scale the per-mg-protein system constants (cellular volume and the three clearances) to the absolute per-well amounts the ODE states carry.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods state that VCells 'was calculated and fixed using the protein content of each preparation and a value of 7.4 uL/mg protein', so protein content is an explicit per-preparation experimental input rather than a model parameter with one true value.",
        "The source paper does not tabulate the protein content per well for the SCHH preparations used here; the seeding density is reported (0.4e6 cells/well in 24-well plates) but the protein-per-cell factor is not.",
        "It affects only the degree to which the medium is depleted during the incubation. The vignette uses 0.4 mg/well, recovered by inverting the paper's own two quantitative Figure 2 claims (Ct,Cells falls to 0.01-fold of baseline at 99 percent uptake inhibition and rises to approximately 15-fold at 99 percent efflux inhibition), which are reproduced only for roughly 0.2-0.4 mg/well; that value is independently consistent with the reported 0.4e6 cells/well seeding density. This is not a fitted parameter - no model parameter is adjusted - but it is not a printed value either, so it is recorded as non-paper-derived provenance in the vignette. The Table 4 fold changes are insensitive to it over at least 0.1-2 mg/well.",
        sep = " "
      ),
      source_name = "protein content of each preparation"
    ),
    BUFFER_CA_FREE = list(
      description = "Buffer-composition flag for the SCHH incubation: 1 = Ca2+/Mg2+-free HBSS containing EGTA, which disrupts the tight junctions so the contents of the bile canaliculi wash into the medium and the measured cell lysate is Cells alone; 0 = standard Ca2+-containing HBSS, which maintains the tight junctions so bile is a separate accumulating compartment and the measured cell lysate is Cells + Bile.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = standard Ca2+-containing HBSS",
      notes = paste(
        "Structural rather than a coefficient multiplier: the flag routes the biliary clearance term. With BUFFER_CA_FREE = 1 the CLBile flux out of the cells discharges directly into the buffer (eq. 2) and the bile state stays at zero; with BUFFER_CA_FREE = 0 it accumulates in the bile state, which drains back to the buffer at KFlux (eqs. 1 and 4).",
        "The cells equation (eq. 3) is identical under both configurations; only the destination of the biliary flux and the buffer return path differ.",
        sep = " "
      ),
      source_name = "'+' and '-' superscripts on the buffer and cell states in eqs. 1-5"
    ),
    WASH_ACTIVE = list(
      description = "Buffer-wash indicator: 1 while the incubation buffer is being washed off the cell monolayer, 0 otherwise. Time-varying within a subject.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no wash in progress",
      notes = paste(
        "Methods: 'To represent the 1-minute wash step, Kwash was activated for 1 minute at the end of the 20-minute uptake phase using an if-then statement.' The covariate is the data-side form of that if-then, so the wash window is a property of the experimental design supplied by the user rather than being hard-coded to the 20-minute uptake this paper used.",
        "The paired rate constant kwash is fixed at 1e4 per minute, which empties the buffer state essentially instantaneously over the 1-minute window.",
        sep = " "
      ),
      source_name = "KWash if-then statement"
    ),
    CONC_INHIBITOR_MEDIUM_UNBOUND_UM = list(
      description = "Unbound concentration of the transporter inhibitor in the incubation medium, [I]u,med. Drives competitive inhibition of the uptake clearance (eq. 7). Set to 0 for the no-inhibitor condition, which collapses the uptake term to its uninhibited form.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 3 values in the presence of 4 percent BSA: telmisartan 0.012 uM at the 1 uM dose level and 0.20 uM at 10 uM; bosentan 0.031 uM at 0.8 uM and 0.45 uM at 8 uM.",
        "The unbound rather than total medium concentration is used because, per the free-drug hypothesis quoted in the Introduction, uptake transporters at the sinusoidal membrane see the medium unbound concentration.",
        sep = " "
      ),
      source_name = "[I]u,med"
    ),
    CONC_INHIBITOR_CELL_UM = list(
      description = "Intracellular concentration of the transporter inhibitor at the efflux-transporter site, [I]cell. Drives competitive inhibition of the biliary (BSEP) and basolateral (MRP3) efflux clearances (eqs. 8 and 10). Set to 0 for the no-inhibitor condition.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper deliberately compares four different measured quantities in this one slot: cellular total [I]t,cell, cellular unbound [I]u,cell, cytosolic total [I]t,cyt and cytosolic unbound [I]u,cyt (Table 3). Which one a given simulation supplies is the paper's central question, so each model run must record it.",
        "Table 3 values, in the order [I]t,cell / [I]u,cell / [I]t,cyt / [I]u,cyt: telmisartan 1 uM = 16 / 2.1 / 16 / 0.85; telmisartan 10 uM = 40 / 3.7 / 35 / 2.8; bosentan 0.8 uM = 1.9 / 0.79 / 1.7 / 0.21; bosentan 8 uM = 17 / 3.8 / 14 / not available.",
        sep = " "
      ),
      source_name = "[I]cell"
    ),
    CONMED_TELMISARTAN = list(
      description = "1 = telmisartan is the co-incubated transporter inhibitor, 0 = it is not. Selects the telmisartan column of the Table 1 inhibitory potencies.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = telmisartan absent",
      notes = "Carries no magnitude information; the inhibitor exposure enters through CONC_INHIBITOR_MEDIUM_UNBOUND_UM and CONC_INHIBITOR_CELL_UM. CONMED_TELMISARTAN and CONMED_BOSENTAN are mutually exclusive; setting both to 0 gives the no-inhibitor condition.",
      source_name = "telmisartan"
    ),
    CONMED_BOSENTAN = list(
      description = "1 = bosentan is the co-incubated transporter inhibitor, 0 = it is not. Selects the bosentan column of the Table 1 inhibitory potencies.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = bosentan absent",
      notes = "Carries no magnitude information; the inhibitor exposure enters through CONC_INHIBITOR_MEDIUM_UNBOUND_UM and CONC_INHIBITOR_CELL_UM. CONMED_TELMISARTAN and CONMED_BOSENTAN are mutually exclusive; setting both to 0 gives the no-inhibitor condition.",
      source_name = "bosentan"
    )
  )

  population <- list(
    species = "in vitro (sandwich-cultured human hepatocytes)",
    n_subjects = 3L,
    n_studies = 1L,
    donors = "Three cryopreserved human hepatocyte lots (HUM4045, HUM4061B, HUM4059; Triangle Research Laboratories): one Caucasian male, one Caucasian female and one Hispanic female, aged 2 to 44 years, body mass index 18.3 to 30",
    system = "B-CLEAR-HU Transporter Certified cryopreserved human hepatocytes seeded at 0.4e6 cells/well in 24-well BioCoat plates (1.75e6 cells/well in 6-well plates for the inhibitor-concentration work), overlaid with Matrigel in a sandwich configuration and studied on day 6 of culture",
    medium = "Hanks' balanced salt solution with 4 percent bovine serum albumin, either standard (Ca2+-containing, tight junctions intact) or Ca2+/Mg2+-free with EGTA (tight junctions disrupted); 0.3 mL per well",
    temperature = "37 C",
    dose_range = "1 umol/L d8-TCA in 0.3 mL per well (0.3 nmol/well)",
    duration = "up to 20 minutes uptake, a 1-minute wash, then 15 minutes efflux; simulations of the inhibitor effect used a 10-minute uptake and the sensitivity analyses a 120-minute steady state",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "The three SCHH preparations were fitted as three independent data sets and Table 2 reports the mean and SD of the three sets of best-fit parameters; the profiles in Figure 1B were generated from the mean parameter estimates.",
      "Fitting was performed in Phoenix WinNonlin v6.3 with the stiff estimation method and a power residual-error model whose parameters are not reported.",
      "The inhibitor-concentration measurements in Table 3 come from a single SCHH preparation (n = 1) in triplicate, and the observed fold changes in Table 4 from a single preparation in duplicate.",
      "Linear (first-order) transport was assumed because the unbound medium TCA concentration (0.15 uM in 4 percent BSA) and the 5.6 uM cellular total concentration after 20 minutes both sit well below the reported Km values of the uptake transporters NTCP (5-20 uM) and OATPs (5.8-71.8 uM) and of the efflux transporters BSEP (6.2 uM), MRP4 (7.7 uM) and MRP3 (30 uM). Passive diffusion was omitted because active uptake dominates hepatocellular accumulation of TCA.",
      sep = " "
    )
  )

  ini({
    # ---- Estimated system parameters -------------------------------------
    # Table 2: mean of the best-fit estimates from n = 3 independent SCHH
    # preparations fitted separately, in the presence of 4 percent BSA. The
    # accompanying SD (used for the Monte Carlo in the vignette) is recorded in
    # the trailing comment; the paper reports no hierarchical random effects.
    lcluptake <- log(0.63);  label("Total uptake clearance of taurocholate (mL/min/g liver)")                       # Table 2: CLUptake 0.63, SD 0.12, CV 20%
    lclbl <- log(0.034);     label("Total basolateral efflux clearance of taurocholate (mL/min/g liver)")           # Table 2: CLBL 0.034, SD 0.011, CV 32%
    lclbile <- log(0.074);   label("Total biliary clearance of taurocholate (mL/min/g liver)")                      # Table 2: CLBile 0.074, SD 0.030, CV 36%
    lkflux <- log(0.018);    label("First-order rate constant of flux from the bile networks into the medium (1/min)")  # Table 2: KFlux 0.018, SD 0.0015, CV 8%

    # ---- Fixed system constants ------------------------------------------
    vcellsPerProtein <- fixed(7.4);  label("Cellular volume per unit cellular protein (uL/mg protein)")             # Methods: VCells calculated from protein content and 7.4 uL/mg protein
    vbuffer <- fixed(300);           label("Incubation buffer volume per well (uL)")                                # Methods: VBuffer set as a constant, 0.3 mL
    proteinPerLiver <- fixed(90);    label("Liver cellular protein content used for the CL unit conversion (mg protein/g liver)")  # Methods: 90 mg protein/g liver, Sohlenius-Sternbeck 2006
    kwash <- fixed(10000);           label("First-order rate constant of the buffer wash step (1/min)")             # Methods: KWash fixed at 1e4 /min
    f_ntcp <- fixed(0.7);            label("Fraction of uptake clearance mediated by NTCP (unitless)")              # Methods: NTCP contributes 70% and OATPs 30% to CLUptake

    # ---- Inhibitory potencies, telmisartan (Table 1) ----------------------
    # Methods state that the mean IC50 values for each transporter in Table 1
    # were used, so a row reporting a range is entered as the mean of its
    # endpoints.
    ic50NtcpTel <- fixed(60);      label("Telmisartan inhibitory potency against NTCP (umol/L)")                    # Table 1: 60 (Ki), Dong 2014
    ic50Oatp1b1Tel <- fixed(0.44); label("Telmisartan inhibitory potency against OATP1B1 (umol/L)")                 # Table 1: 0.44 (Ki), Hirano 2006
    ic50BsepTel <- fixed(16.1);    label("Telmisartan inhibitory potency against BSEP (umol/L)")                    # Table 1: 16-16.2 (IC50), mean of the range
    ic50Mrp3Tel <- fixed(60);      label("Telmisartan inhibitory potency against MRP3 (umol/L)")                    # Table 1: 60 (IC50), Morgan 2013

    # ---- Inhibitory potencies, bosentan (Table 1) -------------------------
    ic50NtcpBos <- fixed(27);      label("Bosentan inhibitory potency against NTCP (umol/L)")                       # Table 1: 18 (Ki) and 36 (IC50), mean of the two reported values
    ic50Oatp1b1Bos <- fixed(18);   label("Bosentan inhibitory potency against OATP1B1 (umol/L)")                    # Table 1: 18, footnote a - not available, assumed the same as NTCP
    ic50BsepBos <- fixed(32.5);    label("Bosentan inhibitory potency against BSEP (umol/L)")                       # Table 1: 23-42 (IC50), mean of the range
    ic50Mrp3Bos <- fixed(42);      label("Bosentan inhibitory potency against MRP3 (umol/L)")                       # Table 1: 42 (IC50), Morgan 2013

    # ---- Residual error ---------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD on the cellular total concentration (fraction); magnitude not reported")  # Methods: a power residual-error model was used in Phoenix WinNonlin, but neither its coefficient nor its exponent is reported
  })

  model({
    # ---- 1. Unit conversion to the per-well basis --------------------------
    # Table 2 reports clearances per gram of liver. The Methods conversion runs
    # the other way: uL/min/mg protein * (90 mg protein/g liver) / 1000 =
    # mL/min/g liver. Inverting it and multiplying by the protein content of
    # the well gives the per-well clearances the ODE states need, in uL/min.
    cluptake <- exp(lcluptake) * 1000 / proteinPerLiver * PROTEIN_MG
    clbl <- exp(lclbl) * 1000 / proteinPerLiver * PROTEIN_MG
    clbile <- exp(lclbile) * 1000 / proteinPerLiver * PROTEIN_MG
    kflux <- exp(lkflux)
    vcells <- vcellsPerProtein * PROTEIN_MG

    # ---- 2. Inhibitor selection -------------------------------------------
    # Mutually exclusive selection of the Table 1 potency column. With neither
    # inhibitor present the placeholder potency of 1 is harmless because both
    # inhibitor concentrations are then 0, so every inhibition ratio is 0.
    noInhibitor <- (1 - CONMED_TELMISARTAN) * (1 - CONMED_BOSENTAN)
    ic50Ntcp <- CONMED_TELMISARTAN * ic50NtcpTel + CONMED_BOSENTAN * ic50NtcpBos + noInhibitor
    ic50Oatp1b1 <- CONMED_TELMISARTAN * ic50Oatp1b1Tel + CONMED_BOSENTAN * ic50Oatp1b1Bos + noInhibitor
    ic50Bsep <- CONMED_TELMISARTAN * ic50BsepTel + CONMED_BOSENTAN * ic50BsepBos + noInhibitor
    ic50Mrp3 <- CONMED_TELMISARTAN * ic50Mrp3Tel + CONMED_BOSENTAN * ic50Mrp3Bos + noInhibitor

    # ---- 3. Competitive inhibition of each transport pathway --------------
    # eq. 7: uptake is split 70% NTCP / 30% OATP and is driven by the medium
    # unbound inhibitor concentration.
    cluptakeInhib <-
      cluptake * f_ntcp / (1 + CONC_INHIBITOR_MEDIUM_UNBOUND_UM / ic50Ntcp) +
      cluptake * (1 - f_ntcp) / (1 + CONC_INHIBITOR_MEDIUM_UNBOUND_UM / ic50Oatp1b1)
    # eq. 8: basolateral efflux is attributed entirely to MRP3, driven by the
    # intracellular inhibitor concentration.
    clblInhib <- clbl / (1 + CONC_INHIBITOR_CELL_UM / ic50Mrp3)
    # eq. 10: biliary efflux is attributed entirely to BSEP.
    clbileInhib <- clbile / (1 + CONC_INHIBITOR_CELL_UM / ic50Bsep)

    # ---- 4. Concentrations ------------------------------------------------
    ctBuffer <- buffer / vbuffer
    ctCells <- cells / vcells

    # ---- 5. ODE system (eqs. 1-4) -----------------------------------------
    # Standard HBSS (BUFFER_CA_FREE = 0): only the basolateral flux and the
    # canalicular KFlux return feed the buffer (eq. 1). Ca2+-free HBSS
    # (BUFFER_CA_FREE = 1): the tight junctions are open, so the biliary flux
    # discharges into the buffer as well and no bile state accumulates (eq. 2).
    d/dt(buffer) <-
      clblInhib * ctCells +
      BUFFER_CA_FREE * clbileInhib * ctCells +
      (1 - BUFFER_CA_FREE) * kflux * bile -
      cluptakeInhib * ctBuffer -
      WASH_ACTIVE * kwash * buffer
    # eq. 3, identical under both buffer configurations.
    d/dt(cells) <- cluptakeInhib * ctBuffer - (clblInhib + clbileInhib) * ctCells
    # eq. 4.
    d/dt(bile) <- (1 - BUFFER_CA_FREE) * (clbileInhib * ctCells - kflux * bile)

    # ---- 6. Observables ---------------------------------------------------
    # eq. 5: what the standard-HBSS cell lysate assay actually measures.
    xCellsBile <- cells + bile
    # Results: the TCA Ct,Cells was chosen as the model output in subsequent
    # simulations to reflect the altered hepatobiliary disposition of TCA in
    # the presence of inhibitors.
    Cc <- ctCells
    Cc ~ prop(propSd)
  })
}
