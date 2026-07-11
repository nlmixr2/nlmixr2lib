Betts_2019_pf_06671008_qsp <- function() {
  description <- "QSP. Preclinical (mouse). Translational quantitative systems pharmacology model for the CD3 x P-cadherin LP DART bispecific PF-06671008, fit to HCT-116 xenografts in NSG mice with human PBMC engraftment. Couples bispecific PK (2-cpt linear), permeability-diffusion tumor drug disposition, trimer (drug-CD3-Pcad) formation in the TME, and tumor growth via a 4-compartment transduction (Simeoni-based) model (Betts 2019 AAPS J)."
  reference <- "Betts A, Haddish-Berhane N, Shah DK, van der Graaf PH, Barletta F, King L, Clark T, Kamperschroer C, Root A, Hooper A, Chen X. A translational quantitative systems pharmacology model for CD3 bispecific molecules: application to quantify T cell-mediated tumor cell killing by P-Cadherin LP DART. AAPS J. 2019;21(4):66. doi:10.1208/s12248-019-0332-z"
  vignette <- "Betts_2019_pf_06671008"
  # Dose is entered as the *concentration equivalent* in the central compartment
  # (nM = dose_nmol / V1_L), not the raw nmol amount. The model states central,
  # peripheral1 and tumor carry concentrations (nM) to match Betts 2019 Eqs 1-4;
  # entering an evid=1 amt corresponds to the C1(0) = dose/V1 step of a bolus.
  # The vignette shows the mg/kg -> nM conversion for a 0.025 kg mouse.
  units <- list(time = "day", dosing = "nM (central-compartment concentration equivalent for IV bolus)", concentration = "nM")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. V1, V2, CL and CLd (reported in Table II mL/kg / mL/h/kg) are scaled by WT to convert to per-mouse L and L/day.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "mouse (female NSG NOD-scid IL-2rg-null; HCT-116 colorectal xenograft with human PBMC engraftment)",
    n_subjects     = "PK study n = 3/timepoint/dose; TGI n = 10/dose",
    n_studies      = 1,
    age_range      = "6-8 weeks at implantation",
    weight_median  = "0.025 kg (typical NSG mouse; used as ref for volume/CL scaling)",
    sex_female_pct = 100,
    disease_state  = "Established subcutaneous HCT-116 colorectal xenografts (~0.5 g at randomization for PK; ~7 days established for TGI). Human PBMCs inoculated 7 days prior to randomization (5e6 or 2.5e6 cells IP).",
    dose_range     = "IV bolus 0.05 and 0.5 mg/kg (PK study); 0.01-0.5 mg/kg q7d x 2 (T cell engrafted TGI). MW PF-06671008 = 105 kDa.",
    notes          = "Primary fit parameters (kmax, kc50, tau, kg0, kg, Mmax; Table II p6-7) are from HCT-116 T cell engrafted model (footnote a in Table II). Companion cell-line fits reported by the paper: HCT-116 adoptive transfer (kmax 1.32 /day, kc50 6.9e-5 nM, tau 3.99 d, TSC 0.011 pM) and SUM-149 adoptive transfer (kmax 0.74 /day, kc50 1.0e-4 nM, tau 4.78 d, TSC 0.0092 pM); see vignette narrative."
  )

  ini({
    # Structural PK, per-kg values (Betts 2019 Table II p6, mouse HCT-116 engrafted PK fit; Monolix v4.3.3)
    lv1_perkg  <- log(49.6);  label("Central volume of distribution per kg body weight (mL/kg)")           # Table II p6: V1 = 49.6 (9% CV)
    lv2_perkg  <- log(60.7);  label("Peripheral volume of distribution per kg body weight (mL/kg)")        # Table II p6: V2 = 60.7 (16% CV)
    lcl_perkg  <- log(0.45);  label("Clearance per kg body weight (mL/h/kg)")                              # Table II p6: CL = 0.45 (12% CV)
    lcld_perkg <- log(4.95);  label("Inter-compartmental clearance per kg body weight (mL/h/kg)")          # Table II p6: CLd = 4.95 (28% CV)

    # Binding kinetics (Betts 2019 Table II p6 and Table III p10; source ref 7 in-house SPR)
    kon_cd3    <- fixed(1.72);   label("Association rate for PF-06671008 to CD3 (1/nM/h)")                  # Table II p6
    koff_cd3   <- fixed(19.66);  label("Dissociation rate for PF-06671008 to CD3 (1/h)")                   # Table II p6
    kon_pcad   <- fixed(1.57);   label("Association rate for PF-06671008 to P-cadherin (1/nM/h)")          # Table II p6
    koff_pcad  <- fixed(0.74);   label("Dissociation rate for PF-06671008 to P-cadherin (1/h)")            # Table II p6

    # Tumor-penetration system parameters (Betts 2019 Table II p6; source ref 19)
    p_perm     <- fixed(334);    label("Drug permeability into tumor (um/day)")                              # Table II p6
    d_diff     <- fixed(0.022);  label("Drug diffusivity into tumor (cm^2/day)")                             # Table II p6
    void_frac  <- fixed(0.24);   label("Tumor void fraction for drug (unitless)")                            # Table II p6
    r_cap      <- fixed(8);      label("Capillary radius (um)")                                              # Table II p6
    r_krogh    <- fixed(75);     label("Krogh cylinder radius (average distance between 2 capillaries, um)") # Table II p6

    # Receptor densities and internalization (Betts 2019 Table II p6; source refs 17,18 CD3; 15 mPcad)
    cd3_receptors <- fixed(100000); label("CD3 receptors per T cell")                                        # Table II p6
    mpcad         <- fixed(28706);  label("P-cadherin receptors per HCT-116 tumor cell")                     # Table II p6 (SUM-149 = 17500)
    tumor_cells_g <- fixed(1e8);    label("Tumor cells per gram of tumor (density)")                         # Table II p6
    kint          <- fixed(0.1728); label("P-cadherin internalization rate (1/day; 96 h half-life)")         # Table II p6

    # Tumor T cell density: fixed system parameter (steady-state approximation).
    # Paper (Betts 2019 Eqs 7-11) uses a 5-day lag on plasma migration then piecewise
    # proliferation for 7 days at Prate = 0.014/(4+dose_ugkg) + 1.5e-5*dose_ugkg
    # then contraction at k_exhaust = 0.0412/h. Setting tcellst_tumor as a scalar
    # keeps the shipped model simulable without dose-conditional logic; users can
    # override via rxode2::rxSetPipelineData() or by editing this constant to
    # reproduce time-varying TIL dynamics.
    tcellst_tumor <- fixed(5e8);    label("Steady-state T cell density in tumor (cells/L)")                  # Betts 2019 Fig 3 baseline

    # Tumor growth (HCT-116 T-cell engrafted; Betts 2019 Table II p6 col a)
    lkg0        <- log(0.30);        label("Exponential tumor growth rate kg0 (1/day)")                       # Table II p6 col a: kg0 = 0.30 (fixed from unperturbed fit)
    lkg         <- log(105);         label("Linear tumor growth rate kg (mm^3/day)")                          # Table II p6 col a: kg = 105 (4% CV)
    lmmax       <- fixed(log(3800)); label("Maximum tumor volume Mmax (mm^3)")                                # Table II p6 col a: Mmax = 3.8e3 (fixed from unperturbed fit)
    psi_switch  <- fixed(20);        label("Switch exponent between exponential and linear tumor growth")     # Table II p6: psi = 20 (fixed based on Simeoni ref 22)

    # Tumor kill (HCT-116 T-cell engrafted; Betts 2019 Table II p7 col a)
    lkmax  <- log(2.71);           label("Maximum tumor kill rate kmax (1/day)")                              # Table II p7 col a: kmax = 2.71 (14% CV)
    lkc50  <- log(2.0e-4);         label("Trimer concentration at half-maximum kill rate kc50 (nM)")          # Table II p7 col a: kC50 = 2.0e-4 (15% CV)
    ltau   <- log(2.25);           label("Transduction time between tumor compartments tau (day)")            # Table II p7 col a: tau = 2.25 (3% CV)

    # IIV - HCT-116 T-cell engrafted (Betts 2019 Table II p7 col a)
    etalkg0 ~ 0.12   # Omega kg0 = 0.12 (25% CV)                                                              # Table II p7 col a
    etalkg  ~ 0.16   # Omega kg  = 0.16 (28% CV)                                                              # Table II p7 col a

    # Residual error - mouse serum PK (Betts 2019 Table II p6). Tumor-volume
    # residuals from Table II p7 col a (a = 60 mm^3, b = 0.01) are documented
    # in the vignette Source Trace but not carried as ini() parameters here
    # because the model exposes tumorVol only as a derived output, not an
    # explicit observation (single-output error model, per pattern 5b of the
    # skill's known-vignette-failure-patterns).
    addSd  <- fixed(0.067);         label("Additive residual error on serum drug (nM)")                       # Table II p6: a = 0.067 (32%)
    propSd <- fixed(0.207);         label("Proportional residual error on serum drug (fraction)")             # Table II p6: b = 0.207 (15%)
  })

  model({
    # ============================================================
    # Individual PK parameters - mL/kg and mL/h/kg -> L and L/day
    # 24 h/day and 1000 mL/L scale factors baked in
    # ============================================================
    v1  <- exp(lv1_perkg)  * WT / 1000              # L
    v2  <- exp(lv2_perkg)  * WT / 1000              # L
    cl  <- exp(lcl_perkg)  * WT * 24 / 1000         # L/day
    cld <- exp(lcld_perkg) * WT * 24 / 1000         # L/day

    kel <- cl  / v1                                  # 1/day
    k12 <- cld / v1                                  # 1/day
    k21 <- cld / v2                                  # 1/day

    # ============================================================
    # Individual PD parameters
    # ============================================================
    kg0_i  <- exp(lkg0 + etalkg0)                    # 1/day
    kg_i   <- exp(lkg  + etalkg)                     # mm^3/day
    mmax   <- exp(lmmax)                             # mm^3
    kmax_i <- exp(lkmax)                             # 1/day
    kc50_i <- exp(lkc50)                             # nM
    tau_i  <- exp(ltau)                              # day

    # ============================================================
    # Tumor geometry - sphere radius from total tumor volume w (mm^3)
    # w_safe floors at 0.001 mm^3 to guard the (1/r^2) diffusion term
    # ============================================================
    w      <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3   # mm^3
    w_safe <- w + 0.001
    r_tumor_cm <- (3.0 * w_safe / (4.0 * 3.141592653589793))^(1.0 / 3.0) / 10.0     # cm
    tv_l   <- w_safe / 1e6                            # L (1 mm^3 = 1e-6 L)

    # ============================================================
    # Tumor drug disposition rate (Betts 2019 Eq 3, p5).
    # (2 P Rcap / Rkrogh^2)  units: (um/day * um) / um^2 = 1/day
    # (6 D / Rtumor^2)        units: cm^2/day / cm^2   = 1/day
    # Multiplied by (C1 - C3/eps) [nM] -> nM/day. C3 == state 'tumor'.
    # ============================================================
    perm_term <- 2.0 * p_perm * r_cap / (r_krogh * r_krogh)              # 1/day
    diff_term <- 6.0 * d_diff / (r_tumor_cm * r_tumor_cm)                # 1/day
    kdisp     <- perm_term + diff_term                                   # 1/day
    tumor_disposition <- kdisp * (central - tumor / void_frac)           # nM/day

    # ============================================================
    # Tumor total Pcad and CD3 (Betts 2019 Eqs 12-13).
    # 1e9 converts nmol/L to nM after dividing by Avogadro.
    # TotPcadt uses fixed density (tumor_cells_g cells/g x mpcad receptors/cell).
    # TotCD3t uses fixed steady-state tcellst_tumor (cells/L).
    # ============================================================
    n_avogadro <- 6.023e23                                                # 1/mol
    totPcadt <- (tumor_cells_g * mpcad) / n_avogadro * 1e9                # nM
    totCD3t  <- (tcellst_tumor * cd3_receptors) / n_avogadro * 1e9        # nM

    # ============================================================
    # Binding rate constants (1/nM/h -> 1/nM/day; 1/h -> 1/day).
    # ============================================================
    kon_cd3_d   <- kon_cd3   * 24.0
    koff_cd3_d  <- koff_cd3  * 24.0
    kon_pcad_d  <- kon_pcad  * 24.0
    koff_pcad_d <- koff_pcad * 24.0

    # ============================================================
    # Free receptor/receptor complex concentrations in the tumor void.
    # Paper divides (TotX - dimer - trimer) by epsilon (void fraction) to get accessible concentration.
    # ============================================================
    free_pcad <- (totPcadt - drug_pcad_tumor - trimer) / void_frac
    free_cd3  <- (totCD3t  - drug_cd3_tumor  - trimer) / void_frac
    c3        <- tumor / void_frac                                        # accessible drug in tumor void (nM)

    # ============================================================
    # ODE system.
    # 'central' holds C1 (nM); 'peripheral1' holds C2 (nM); 'tumor' holds C3 (nM).
    # sPcad binding (Eqs 25-26) is excluded per Betts 2019 assumption "PF-06671008 does not bind to sPcad in mouse".
    # CD3 dimer in plasma (dcd3p, Eq 6) is excluded because plasma CD3 saturation is negligible (<1%).
    # T cell trafficking ODEs (Eqs 7-8) are collapsed into the tcellst_tumor scalar.
    # ============================================================

    # C1 in central (Betts 2019 Eq 1, mouse form without sPcad; dcd3p contribution dropped)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 * v2 / v1 - tumor_disposition * (tv_l / v1)

    # C2 in peripheral (Betts 2019 Eq 2)
    d/dt(peripheral1) <-  k12 * central * v1 / v2 - k21 * peripheral1

    # C3 in tumor (Betts 2019 Eq 4)
    d/dt(tumor)       <-  tumor_disposition - kon_pcad_d * c3 * free_pcad + koff_pcad_d * drug_pcad_tumor - kon_cd3_d * c3 * free_cd3 + koff_cd3_d * drug_cd3_tumor

    # drug_cd3_tumor (Betts 2019 Eq 14)
    d/dt(drug_cd3_tumor)       <-  kon_cd3_d * c3 * free_cd3 - koff_cd3_d * drug_cd3_tumor - kon_pcad_d * drug_cd3_tumor * free_pcad + koff_pcad_d * trimer

    # drug_pcad_tumor (Betts 2019 Eq 15) - includes internalization of drug-Pcad dimer at rate kint
    d/dt(drug_pcad_tumor)      <-  kon_pcad_d * c3 * free_pcad - koff_pcad_d * drug_pcad_tumor - kon_cd3_d * drug_pcad_tumor * free_cd3 + koff_cd3_d * trimer - kint * drug_pcad_tumor

    # trimer (Betts 2019 Eq 16)
    d/dt(trimer)      <-  kon_cd3_d * drug_pcad_tumor * free_cd3 + kon_pcad_d * drug_cd3_tumor * free_pcad - koff_cd3_d * trimer - koff_pcad_d * trimer

    # ============================================================
    # Tumor kill and 4-compartment transduction (Betts 2019 Eqs 17-22).
    # kkill (Eq 17) drives loss of the growth cpt cycling_cells; the killed mass propagates through
    # damaged_cells1..3 with residence tau (paper's M1..M4 = cycling_cells + damaged_cells1..3).
    # ============================================================
    kkill <- kmax_i * trimer / (kc50_i + trimer)                                       # 1/day (Eq 17)

    # Modified Simeoni logistic switch (Eq 18): exponential-then-linear growth minus kill
    grow_num <- kg0_i * (1 - w_safe / mmax) * cycling_cells
    grow_den <- (1 + (kg0_i / kg_i)^psi_switch * w_safe^psi_switch)^(1.0 / psi_switch)
    d/dt(cycling_cells)  <- grow_num / grow_den - kkill * cycling_cells                 # (Eq 18)
    d/dt(damaged_cells1) <- kkill * cycling_cells - damaged_cells1 / tau_i               # (Eq 19)
    d/dt(damaged_cells2) <- (damaged_cells1 - damaged_cells2) / tau_i                    # (Eq 20)
    d/dt(damaged_cells3) <- (damaged_cells2 - damaged_cells3) / tau_i                    # (Eq 21)

    # ============================================================
    # Observations. Cc is serum drug (nM); tumorVol is total tumor volume w (mm^3).
    # Only Cc carries an explicit residual-error model to keep the endpoint
    # single-output (simplifies dvid handling in simulations); tumorVol is
    # returned as a derived column with its Table II residuals available in
    # ini() as addSd_tumorVol / propSd_tumorVol for downstream error simulation.
    # ============================================================
    Cc       <- central
    tumorVol <- w

    Cc       ~ add(addSd) + prop(propSd)
  })
}
