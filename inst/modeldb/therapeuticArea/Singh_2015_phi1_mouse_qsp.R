Singh_2015_phi1_mouse_qsp <- function() {
  description <- paste(
    "QSP (preclinical, mouse). Systems-pharmacology model of erythropoiesis",
    "induced by an oral small-molecule prolyl-hydroxylase-2 (PHD2) inhibitor",
    "(the undisclosed Pfizer compound the paper calls PHI-1) in male BALB/c",
    "mice (Singh 2015). A two-compartment oral PK model (first-order ka, with",
    "a dose-dependent apparent volume vc: 9.15 L at 10 mg/kg, 2.4 L at 30",
    "mg/kg) drives a direct-effect inhibition of active PHD2; the resulting",
    "rise in kidney HIF1alpha (turnover with PHD2-mediated degradation only)",
    "stimulates kidney EPO mRNA through a power function, which is translated",
    "to plasma EPO through a single transit compartment. Plasma EPO linearly",
    "stimulates a cell-lifespan model of reticulocytes (mean circulating",
    "lifespan Tr = 81 h), which mature into red blood cells (assumed lifespan",
    "TRBC = 24 days); hemoglobin is the mean corpuscular hemoglobin times the",
    "red-cell count. The active-PHD2 fraction (PHD2/PHD2tot0) and the RBC and",
    "hemoglobin time courses are simulated outputs. Naive-pooled fit, no",
    "between-subject variability; the residual-error SDs were estimated but",
    "not reported, so this is packaged as a typical-value simulation model."
  )
  reference <- paste(
    "Singh AP, Zheng J, Lin J, Wu C, Harrold JW, Yan X, Mager DE, Barrett JS,",
    "Mendel D, Davis CB, Krzyzanski W. A Systems Pharmacology Model of",
    "Erythropoiesis in Mice Induced by Small Molecule Inhibitor of Prolyl",
    "Hydroxylase Enzymes. CPT Pharmacometrics Syst Pharmacol. 2015;4(1):e12.",
    "doi:10.1002/psp4.12"
  )
  vignette <- "Singh_2015_erythropoiesis"
  units <- list(time = "h", dosing = "umol", concentration = "uM")

  # Mechanistic erythropoiesis states that are not canonical PK compartments;
  # each is documented in compartmentData above/below.
  paper_specific_compartments <- c(
    "hif1",
    "epo_mrna",
    "epo_transit",
    "epo_plasma",
    "retic",
    "rbc_m"
  )

  compartmentData <- list(
    depot = list(
      analyte = "PHI-1 (PHD2 inhibitor)",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "PHI-1 (PHD2 inhibitor)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "PHI-1 (PHD2 inhibitor)", units = "umol", specimen = "tissue", verified = TRUE),
    hif1 = list(
      analyte = "hypoxia-inducible factor 1-alpha protein (kidney)",
      units = "AU",
      specimen = "tissue",
      verified = TRUE
    ),
    epo_mrna = list(
      analyte = "erythropoietin mRNA (kidney, relative to baseline)",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    epo_transit = list(
      analyte = "erythropoietin (translation/secretion transit)",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    epo_plasma = list(analyte = "endogenous erythropoietin", units = "ng/mL", specimen = "plasma", verified = TRUE),
    retic = list(
      analyte = "reticulocytes",
      units = "% of red blood cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    rbc_m = list(analyte = "mature red blood cells", units = "cells/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "mouse (male BALB/c)",
    n_subjects = "3-6 per timepoint (naive-pooled; terminal sampling)",
    age_range = "8-10 weeks at study start",
    disease_state = "healthy (non-anemic) mice; erythropoiesis stimulated pharmacologically",
    dose_range = paste(
      "Single oral doses of 10 and 30 mg/kg PHI-1 (coded as 25.8 and 77.4",
      "umol amounts in the authors' ADAPT model, giving plasma concentrations",
      "in uM when divided by the apparent volume)."
    ),
    regions = "preclinical (Pfizer, Cambridge, MA)",
    notes = paste(
      "Data: plasma PHI-1 (0.25-96 h), kidney HIF1alpha, kidney EPO mRNA (mean",
      "only), plasma EPO, reticulocytes (% RBC), and RBC and hemoglobin on day",
      "0 and day 10. Fitted sequentially (PK, then HIF1alpha, then EPO",
      "mRNA/plasma, then reticulocytes) by naive-pooled maximum likelihood in",
      "ADAPT 5; RBC, hemoglobin and the PHD2 activity ratio are predictive",
      "simulations, not fitted endpoints."
    )
  )

  ini({
    # --- Pharmacokinetics (Singh 2015 Table 1, two-compartment oral) --------
    lka <- log(8.5)
    label("Oral absorption rate constant (1/h)") # Table 1: ka 8.5 hr-1 (CV 54%)
    lkel <- log(0.66)
    label("Elimination rate constant from central (1/h)") # Table 1: kel 0.66 hr-1 (CV 9%)
    lkcp <- log(0.014)
    label("Central-to-peripheral rate constant (1/h)") # Table 1: kcp 0.014 hr-1 (CV 31%)
    lkpc <- log(0.005)
    label("Peripheral-to-central rate constant (1/h)") # Table 1: kpc 0.005 hr-1 (CV 31%)
    # Dose-dependent apparent volume of distribution (central). Default is the
    # 30 mg/kg value V2/F; set vc = 9.15 to simulate the 10 mg/kg dose.
    lvc <- log(2.4)
    label("Apparent central volume vc = V/F (L)") # Table 1: V2/F 2.4 L at 30 mg/kg (CV 20%); V1/F 9.15 L at 10 mg/kg (CV 35%)

    # --- PHD2 / HIF1alpha (Singh 2015 Table 1, Eqs. 1-4, S15-S18) -----------
    kePHD2 <- 9.9
    label("Hybrid PHD2-mediated HIF1alpha degradation rate KE x PHD2tot0 (1/h)") # Table 1: KE x PHD2tot0 9.9 (CV 14%)
    IC50 <- 1.7
    label("PHI-1 concentration for 50% PHD2 inhibition (uM)") # Table 1: IC50 1.7 uM (CV 6%)
    g1 <- 2.3
    label("Hill coefficient of PHD2 inhibition (unitless)") # Table 1: gamma1 2.3 (CV 11%)
    kdeg_HIF <- fixed(0)
    label("Drug-independent HIF1alpha degradation rate (1/h)") # Table 1: kdeg,HIF1alpha 0 (fixed); PHD2-mediated degradation dominates
    HIF1o <- 351
    label("Baseline kidney HIF1alpha (AU)") # Table 1: HIF1o 351 AU (CV 7%)

    # --- EPO mRNA and plasma EPO (Singh 2015 Table 1, Eqs. 5-13) -----------
    kdeg_mRNA <- 0.43
    label("EPO mRNA degradation rate (1/h)") # Table 1: kdeg,EPO mRNA 0.43 hr-1 (CV 15%)
    tau <- 2.7
    label("Mean transit time EPO mRNA to plasma EPO (h)") # Table 1: tau 2.7 hr (CV 7%); kdeg,EPO plasma = 1/tau = 0.37 hr-1
    EPO_mRNA0 <- 0.8
    label("Baseline kidney EPO mRNA (fold)") # Table 1: EPO_mRNA0 0.8 (CV 13%)
    EPO_plasma0 <- 0.047
    label("Baseline plasma EPO (ng/mL)") # Table 1: EPO_Plasma0 0.047 ng/mL (CV 4%)
    g2 <- 2.8
    label("Power exponent of HIF1alpha on EPO mRNA synthesis (unitless)") # Table 1: gamma2 2.8 (CV 5%)

    # --- Reticulocytes (Singh 2015 Table 1, Eqs. 14-16, S1) ----------------
    Retic0 <- 4.6
    label("Baseline reticulocytes (% of RBC)") # Table 1: Retic0 4.6 % of RBC (CV 4%)
    Tr <- 81
    label("Mean circulating reticulocyte lifespan (h)") # Table 1: Tr 81 hr (CV 4%)
    SL <- 0.2
    label("Linear EPO slope on reticulocyte production (mL/ng)") # Table 1: SL 0.2 (CV 18%)

    # --- Red blood cells and hemoglobin (simulation; Eqs. 17-22, Suppl. S1)
    RBC0 <- fixed(10.76e12)
    label("Baseline red blood cell count (cells/L)") # Supplement S1: RBC(0) fixed at observed mean 10.76 x 10^12 cells/L
    MCH_average <- fixed(16.04e-12)
    label("Mean corpuscular hemoglobin (g/cell)") # Table 1: MCHaverage 16 x 10^-12 g/cell (fixed); ADAPT value 16.04e-12
    Nkin <- 0.5226e10
    label("Mature-RBC production rate (cells/L/h)") # Supplement S1: Nkin = kin x RBC0, ADAPT-coded 0.5226e10 (kin 0.05 %cells/hr)
  })

  model({
    # Secondary (baseline-balance) parameters. kdeg_HIF is fixed to 0, so
    # ksyn_hif reduces to kePHD2 * HIF1o = 3473 AU/hr (Table 1 secondary).
    ksyn_hif <- (kdeg_HIF + kePHD2) * HIF1o
    ksyn_mrna <- kdeg_mRNA * EPO_mRNA0 # Table 1 secondary: ksyn,EPO mRNA 0.35
    ksyn_epo <- (1 / tau) * EPO_plasma0 / EPO_mRNA0 # Table 1 secondary: ksyn,EPO 0.02
    kdeg_epopl <- 1 / tau # Eq. 11: kdeg,EPO plasma = 1/tau (Table 1 secondary 0.37)
    kin_ret <- Retic0 / (Tr * (1 + SL * EPO_plasma0)) # Eq. 18: reticulocyte production (Table 1 secondary kin 0.05)

    ka <- exp(lka)
    kel <- exp(lkel)
    kcp <- exp(lkcp)
    kpc <- exp(lkpc)
    vc <- exp(lvc)

    # Plasma drug concentration (uM) and direct-effect active-PHD2 fraction
    # PHD2tot/PHD2tot0 (Eq. 1 / S15; simulated as Fig. 4c).
    Cc <- central / vc
    fPHD2 <- 1 - Cc^g1 / (IC50^g1 + Cc^g1)

    # Two-compartment oral PK
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + kpc * peripheral1 - (kel + kcp) * central
    d/dt(peripheral1) <- kcp * central - kpc * peripheral1

    # Kidney HIF1alpha turnover (Eq. 3 / S17). Natural degradation kdeg_HIF is
    # fixed to 0; the drug raises hif1 by lowering PHD2-mediated degradation.
    d/dt(hif1) <- ksyn_hif - (kdeg_HIF + kePHD2 * fPHD2) * hif1
    hif1(0) <- HIF1o

    # Kidney EPO mRNA: power-function stimulation by HIF1alpha (Eq. 5)
    d/dt(epo_mrna) <- ksyn_mrna * (hif1 / HIF1o)^g2 - kdeg_mRNA * epo_mrna
    epo_mrna(0) <- EPO_mRNA0

    # Plasma EPO: one signal-transduction transit compartment then plasma
    # (Eqs. 6-7). The transit compartment holds EPO_plasma0 at baseline (Eq. 9).
    d/dt(epo_transit) <- ksyn_epo * epo_mrna - (1 / tau) * epo_transit
    epo_transit(0) <- EPO_plasma0
    d/dt(epo_plasma) <- (1 / tau) * epo_transit - kdeg_epopl * epo_plasma
    epo_plasma(0) <- EPO_plasma0

    # Reticulocyte cell-lifespan model (Eqs. 14-16). S(t) = 1 + SL*epo_plasma;
    # delay() returns the epo_plasma value Tr hours earlier (= EPO_plasma0
    # before dosing), so the baseline is a steady state.
    S_now <- 1 + SL * epo_plasma
    S_tr <- 1 + SL * delay(epo_plasma, Tr)
    d/dt(retic) <- kin_ret * S_now - kin_ret * S_tr
    retic(0) <- Retic0

    # Mature RBC. The authors' ADAPT simulation code that produced Fig. 4a/b
    # (subroutine DIFFEQ, XP(42)) drives mature-RBC production from the CURRENT
    # plasma EPO with no reticulocyte-maturation delay and no RBC-death term,
    # rather than the delayed Eq. 19 form; the packaged model follows the code
    # so that the RBC and hemoglobin time courses reproduce Fig. 4. Written in
    # a baseline-stable form (the ADAPT `-Nkin` term leaves a negligible
    # baseline drift of ~Nkin*SL*EPO_plasma0). The assumed 24-day RBC lifespan
    # (TRBC) therefore does not remove cells within the published 10-day
    # window; see the vignette Errata.
    d/dt(rbc_m) <- Nkin * (1 + SL * epo_plasma) - Nkin * (1 + SL * EPO_plasma0)
    rbc_m(0) <- RBC0 * (1 - Retic0 / 100)

    # Total RBC (Eq. 20) and hemoglobin (Eq. 22). The reticulocyte % is
    # converted to a count with the baseline RBC denominator.
    RBC <- rbc_m + (retic / 100) * RBC0
    HGB <- MCH_average * RBC / 10

    # Active-PHD2 ratio reported in Fig. 4c
    PHD2ratio <- fPHD2
  })
}
