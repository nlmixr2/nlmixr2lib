Yamada_2025_oncolyticVirus_mouse_qsp <- function() {
  description <- paste(
    "Preclinical (mouse). QSP. Semi-mechanistic quantitative systems",
    "pharmacology model of CT26.WT-bearing mouse tumour dynamics under an",
    "intratumorally injected oncolytic vaccinia virus carrying human IL-7 and",
    "murine IL-12 genes (hIL-7/mIL-12-VV), alone and in combination with an",
    "anti-PD-1 antibody (DX400/RMP1-14, the murine surrogate for pembrolizumab).",
    "Ten published ODEs (the total-tumour equation is encoded algebraically as",
    "the sum of its parts): the injected tumour (noninfected and infected tumour",
    "volume, virus amount, nondifferentiated and differentiated effector T cells,",
    "a PD-L1 suppression variable, an anti-PD-1 depot/central PK pair) and a",
    "distant (noninjected) tumour driven by the systemic differentiated",
    "effector-T-cell pool. The cancer-immunity cycle (tumour-cell death ->",
    "dendritic-cell maturation -> systemic antigen -> T-cell influx,",
    "proliferation and differentiation, with IL-7 and IL-12 stimulating",
    "proliferation and differentiation and PD-L1/immune-suppressive-cell",
    "down-regulation of the immune activation rate) follows the Kosinsky 2018",
    "radiation/anti-PD-(L)1 CT26 model (Ref 3), extended with the oncolytic-virus",
    "infection and immunostimulation terms of Kim 2015 (Ref 2). The antagonistic",
    "effect of the anti-PD-1 antibody is encoded hypothetically as enhancement of",
    "virus clearance.",
    "Deterministic mechanism model: parameters were fitted in NONMEM to pooled",
    "preclinical tumour-growth data (minimisation terminated before convergence;",
    "the authors used the terminated estimates), so no etas and no residual-error",
    "model are encoded.",
    sep = " "
  )
  reference <- paste(
    "Yamada A, Tanaka S, Nakamura S, Asano N, Nomura T, Takahashi M, Caunt MP,",
    "Boyd PL, Brennan FA, Ortega F, Markham H, Orrell D, Chakrabarti C, Kato K.",
    "A Multiple-Model-Informed Drug-Development Approach for Optimal Regimen",
    "Selection of an Oncolytic Virus in Combination With Pembrolizumab.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):572-583.",
    "doi:10.1002/psp4.13297.",
    sep = " "
  )
  vignette <- "Yamada_2025_oncolytic_virus_pembrolizumab"

  units <- list(
    time = "day",
    dosing = "pfu (virus); nmol (anti-PD-1)",
    volume = "mm^3",
    concentration = "nM (anti-PD-1 serum concentration)"
  )

  covariateData <- list()

  compartmentData <- list(
    virus = list(
      analyte = "Oncolytic-virus amount at the injected tumour",
      units = "pfu",
      specimen = "tumor",
      verified = FALSE
    ),
    tumNi = list(
      analyte = "Volume of noninfected tumour cells in the virus-injected tumour",
      units = "mm^3",
      specimen = "tumor",
      verified = FALSE
    ),
    tumI = list(
      analyte = "Volume of infected tumour cells in the virus-injected tumour",
      units = "mm^3",
      specimen = "tumor",
      verified = FALSE
    ),
    nTeff = list(
      analyte = "Nondifferentiated precursors of cytotoxic effector T cells",
      units = "cells",
      specimen = "tumor",
      verified = FALSE
    ),
    dTeff = list(
      analyte = "Differentiated cytotoxic effector T cells",
      units = "cells",
      specimen = "tumor",
      verified = FALSE
    ),
    pdl1 = list(
      analyte = "PD-L1-mediated immune-suppression variable (relative, 0-1)",
      units = "fraction",
      specimen = "tumor",
      verified = FALSE
    ),
    mAbDep = list(
      analyte = "Anti-PD-1 antibody in the injection-site (depot) compartment",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    mAbC = list(
      analyte = "Anti-PD-1 antibody in the central (serum) compartment",
      units = "nmol",
      specimen = "serum",
      verified = FALSE
    ),
    tumDis = list(
      analyte = "Volume of the distant (noninjected) tumour",
      units = "mm^3",
      specimen = "tumor",
      verified = FALSE
    )
  )

  population <- list(
    species = "mouse (BALB/c, CT26.WT syngeneic colon carcinoma)",
    disease_state = "CT26.WT subcutaneous tumour (injected and bilateral distant tumour models)",
    dose_range = paste(
      "hIL-7/mIL-12-VV 2e4 to 2e7 pfu intratumorally on fixed days (see Table S1);",
      "anti-PD-1 antibody (DX400/RMP1-14) 100 ug intraperitoneally.",
      sep = " "
    ),
    notes = paste(
      "Calibrated to four independent CT26.WT tumour-growth-inhibition",
      "experiments (Exp1-Exp4, 19 experimental conditions including 6 vehicle",
      "arms). Tumour volume was converted to tumour-cell counts assuming 1e6",
      "cells/mm^3. Mice were randomised to treatment when tumour volume exceeded",
      "100 mm^3 (Day 0).",
      sep = " "
    )
  )

  # Mechanistic QSP states that are not canonical nlmixr2lib compartments.
  paper_specific_compartments <- c(
    "virus",
    "tumNi",
    "tumI",
    "nTeff",
    "dTeff",
    "pdl1",
    "mAbDep",
    "mAbC",
    "tumDis"
  )

  ini({
    # ---- tumour growth and spontaneous death ----
    # r was fitted per experiment (Exp1-Exp4: 0.36, 0.284, 0.348, 0.394 1/day);
    # the Exp1 value is the default. Override r to reproduce Exp2-Exp4.
    r <- fixed(0.36); label("Tumour-growth rate constant [1/day]")  # Table S2 (r = 0.36 0.284 0.348 0.394 1/day, Fitted; Exp1 value used)
    d0 <- fixed(0.01); label("Spontaneous death rate of tumour cells [1/day]")  # Table S2 (d0 = 0.01 1/day, Ref3)

    # ---- T-cell influx, proliferation, differentiation, elimination ----
    kLN <- fixed(279); label("Maximum influx rate of nTeff cells [cells/day]")  # Table S2 (kLN = 279 cells/day, Ref3)
    SL <- fixed(2.76); label("Sensitivity of nTeff influx to systemic antigen [cells]")  # Table S2 (SL = 2.76 cells, Ref3)
    kpro <- fixed(3.0); label("nTeff proliferation rate constant [1/day]")  # Table S2 (kpro = 3.0 1/day, Ref3)
    kdif <- fixed(3.2); label("nTeff differentiation rate constant [1/day]")  # Table S2 (kdif = 3.2 1/day, Ref3)
    kel_nTeff <- fixed(0.2); label("nTeff elimination rate constant [1/day]")  # Table S2 (kel = 0.2 1/day, Ref3)
    kapo <- fixed(2.0); label("dTeff apoptosis rate constant [1/day]")  # Table S2 (kapo = 2.0 1/day, Ref3)
    e <- fixed(0.0013); label("Rate of tumour-cell kill by dTeff [1/day/cell]")  # Table S2 (e = 0.0013 1/day/cell, Ref3)

    # ---- PD-L1 immune-suppression module ----
    K_pdl <- fixed(478); label("Sensitivity of PD-L1 up-regulation to dTeff count [cells]")  # Table S2 (K_pdl = 478 cells, Ref3)
    kpdl <- fixed(1.0); label("PD-L1 up-regulation rate constant [1/day]")  # Table S2 (kpdl = 1.0 1/day, Ref3)
    K_tcd <- fixed(0.2); label("Sensitivity of dendritic-cell maturation to tumour-cell death [1/day]")  # Table S2 (K_tcd = 0.2 1/day, Ref3)
    SR <- fixed(20.3); label("Sensitivity of immune-suppressive-cell accumulation to systemic antigen [cells]")  # Table S2 (SR = 20.3 cells, Ref3)

    # ---- anti-PD-1 antibody PK and binding ----
    kab <- fixed(8.0); label("Anti-PD-1 absorption rate constant from the depot compartment [1/day]")  # Table S2 (kab = 8.0 1/day, Ref3)
    Vd <- fixed(0.003); label("Anti-PD-1 distribution volume [L]")  # Table S2 (Vd = 0.003 L, Ref3)
    kel_mAb <- fixed(0.15); label("Anti-PD-1 elimination rate constant from the central compartment [1/day]")  # Table S2 (kel_mAb = 0.15 1/day, Ref3)
    KD <- fixed(30); label("Anti-PD-1 binding affinity to PD-1 [nM]")  # Table S2 (KD = 30 nM, Ref3)

    # ---- IL-7 / IL-12 immunostimulation (oncolytic-virus specific) ----
    kIL7 <- fixed(2.09); label("Maximum activation rate constant for nTeff proliferation by IL-7 [1/day]")  # Table S2 (kIL7 = 2.09 1/day, Fitted)
    kIL12 <- fixed(5.01); label("Maximum activation rate constant for nTeff differentiation by IL-12 [1/day]")  # Table S2 (kIL12 = 5.01 1/day, Fitted)
    K_IL7 <- fixed(0.01); label("Michaelis-Menten constant for nTeff proliferation activation by IL-7 [cells]")  # Table S2 (K_IL7 = 0.01 cells, Fixed)
    K_IL12 <- fixed(0.01); label("Michaelis-Menten constant for nTeff differentiation by IL-12 [cells]")  # Table S2 (K_IL12 = 0.01 cells, Fixed)

    # ---- oncolytic-virus kinetics ----
    V_alpha <- fixed(8.6); label("Viral production size [dimensionless]")  # Table S2 (V_alpha = 8.6, Fitted)
    V_beta <- fixed(0.271); label("Viral infection rate constant [1/day/cell]")  # Table S2 (V_beta = 0.271 1/day/cell, Fitted)
    delta_i <- fixed(57.4); label("Infected-cell lysis rate constant [1/day]")  # Table S2 (delta_i = 57.4 1/day, Fitted)
    delta_v <- fixed(2.39); label("Virus elimination rate constant [1/day]")  # Table S2 (delta_v = 2.39 1/day, Fitted)
    IC50 <- fixed(30); label("Half-maximal anti-PD-1 serum concentration enhancing virus clearance [nM]")  # Table S2 (IC50 = 30 nM, Fitted)
    Imax <- fixed(0.139); label("Maximum fold enhancement of virus clearance by anti-PD-1 [dimensionless]")  # Table S2 (Imax = 0.139, Fitted)
  })
  model({
    # Units: time in days; tumour volumes in mm^3; virus in pfu; anti-PD-1 in
    # nmol; nTeff/dTeff in cells. Tumour volume is converted to cell counts with
    # 1e6 cells/mm^3 where the virus-infection rate laws require counts.

    # ---- algebraic quantities ----
    # injected-tumour total volume (Supplement eq 9: Tum_tot = Tum_ni + Tum_i)
    tumTot <- tumNi + tumI
    # Eradication threshold of one tumour cell (1e-6 mm^3 at 1e6 cells/mm^3).
    # The infection flux V_beta * Virus * Tum_ni / Tum_tot does not shrink with
    # tumour size while Tum_ni dominates, so the injected tumour reaches zero in
    # finite time and the Tum_ni/Tum_tot and Tum_i/Tum_tot ratios become
    # infinitely stiff. Below one cell the injected tumour is treated as
    # eradicated: its two volumes are frozen and the ratio denominators are
    # floored at one cell.
    alive <- 1 * (tumTot >= 1e-6)
    tumTotD <- max(tumTot, 1e-6)
    # anti-PD-1 serum concentration (amount / distribution volume) [nM]
    mAbConc <- mAbC / Vd
    # tumour-cell-death rate, dendritic-cell maturation, systemic antigen level
    TCD <- delta_i * tumI / tumTotD + e * dTeff + d0   # Table S2 rule (TCD)
    DCm <- TCD^2 / (TCD^2 + K_tcd^2)   # Table S2 rule (DCm)
    Agsys <- DCm * tumTot   # Table S2 rule (Agsys)
    # free PD-L1 after anti-PD-1 blockade and immune-suppressive-cell fraction
    PD_L1free <- pdl1 / (1 + mAbConc / KD)   # Table S2 rule (PD_L1free)
    ISC <- Agsys / (Agsys + SR)   # Table S2 rule (ISC)
    # immune activation rate (0-1)
    IAR <- (1 - PD_L1free) * (1 - ISC)   # Table S2 rule (IAR)

    # ---- ODEs (Supplement eqs 1-10) ----
    # virus: produced by lysing infected cells, cleared with anti-PD-1-enhanced
    # elimination (eq 1)
    d/dt(virus) <- V_alpha * delta_i * tumI * 1e6 -
      delta_v * (1 + Imax * mAbConc / (IC50 + mAbConc)) * virus
    # noninfected injected-tumour volume: growth, immune + spontaneous kill,
    # loss to infection (eq 2; the kill term uses dTeff as in the Kosinsky source)
    d/dt(tumNi) <- alive * (r * tumNi - (e * dTeff + d0) * tumNi -
      V_beta * virus * tumNi / tumTotD / 1e6)
    # infected injected-tumour volume: gained from infection, lost to lysis and
    # immune/spontaneous kill (eq 3)
    d/dt(tumI) <- alive * (V_beta * virus * tumNi / tumTotD / 1e6 -
      delta_i * tumI - (e * dTeff + d0) * tumI)
    # nondifferentiated effector T cells: antigen-driven influx, IAR-gated net
    # proliferation minus differentiation (IL-7 and IL-12 modulated), elimination
    # (eq 4)
    d/dt(nTeff) <- kLN * Agsys / (Agsys + SL) +
      IAR * nTeff * ((kpro + kIL7 * tumI / (K_IL7 + tumI)) -
        (kdif + kIL12 * tumI / (K_IL12 + tumI))) -
      kel_nTeff * nTeff
    # differentiated cytotoxic effector T cells: IAR-gated differentiation,
    # apoptosis (eq 5)
    d/dt(dTeff) <- IAR * nTeff * (kdif + kIL12 * tumI / (K_IL12 + tumI)) -
      kapo * dTeff
    # PD-L1 suppression variable driven by dTeff count (eq 6; the Kosinsky 2018
    # source eq 16 uses dTeff in the numerator -- the Yamada supplement rendering
    # of eq 6 printed PDL1 there, which is a transcription artifact of dTeff)
    d/dt(pdl1) <- kpdl * (dTeff / (dTeff + K_pdl) - pdl1)
    # anti-PD-1 antibody depot -> central PK (eqs 7-8)
    d/dt(mAbDep) <- -kab * mAbDep
    d/dt(mAbC) <- kab * mAbDep - kel_mAb * mAbC
    # distant (noninjected) tumour: growth and kill by the systemic dTeff pool,
    # no local virus (eq 10)
    d/dt(tumDis) <- r * tumDis - (e * dTeff + d0) * tumDis

    # ---- initial conditions ----
    # Tumours were randomised to treatment when the volume exceeded 100 mm^3
    # (Day 0); 100 mm^3 is used as the representative baseline for the injected
    # and distant tumours. The immune and virus states start from rest.
    tumNi(0) <- 100
    tumDis(0) <- 100

    # ---- outputs ----
    tumour_volume <- tumTot   # injected-tumour volume [mm^3]
    tumour_volume_distant <- tumDis   # distant-tumour volume [mm^3]
  })
}
