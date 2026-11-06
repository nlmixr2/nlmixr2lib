Feldmuller_2025_zabedosertib <- function() {
  description <- "One-compartment population PK model for oral zabedosertib (BAY 1834845, an IRAK4 inhibitor) in healthy adult men (Feldmuller 2025). First-order absorption through one transit compartment (depot -> transit1 -> central, both steps at rate ka), a dose-dependent relative bioavailability that falls exponentially with the administered dose towards a fixed floor of 1%, linear elimination of UNBOUND drug, and capacity-limited binding of unbound drug to a high-affinity plasma binding partner ('AGP') in the central compartment with explicit association / dissociation kinetics (Kd fixed to the in vitro AGP affinity, 0.4 uM; dissociation rate and binding capacity estimated). The observed total plasma concentration is unbound plus bound. A target-occupancy sub-model links the model-predicted unbound plasma concentration to IRAK4 binding (Kd taken as the in vitro whole-blood IL-6 IC50, 86 nM; target residence time 5 min). IIV on unbound CL, unbound V, ka and binding capacity; combined proportional plus study-specific additive residual error."
  reference <- "Feldmuller M, Jodl SJ, Ploeger B, Wagenfeld A, Wiesinger H, Zollmann FS, Klein S, Zhang R, Rohde B, Hochel J. Zabedosertib, a novel interleukin-1 receptor-associated kinase-4 inhibitor, shows a favorable pharmacokinetic and safety profile across multiple phase 1 studies. Front Pharmacol. 2025;16:1521505. doi:10.3389/fphar.2025.1521505. PMCID PMC12163019. Model structure from Supplementary Section 1.5 (Figure S1 and text); parameter estimates from Supplementary Table S4; target-occupancy inputs from Section 2.11, Section 3.5 and Supplementary Section 1.6."
  vignette <- "Feldmuller_2025_zabedosertib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The AGP-bound pool and the IRAK4 occupancy fraction are specific to this
  # paper; neither is a registered canonical state.
  paper_specific_compartments <- c("complex_agp", "occ_irak4")

  compartmentData <- list(
    depot = list(
      analyte = "zabedosertib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral tablet dose. Bioavailability f(depot) is applied here (Supplementary Section 1.5, F1 equation)."
    ),
    transit1 = list(
      analyte = "zabedosertib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "The single 'Transit' box between Depot and the central compartment in Supplementary Figure S1. Table S4 lists only one absorption rate constant, so both depot -> transit1 and transit1 -> central run at ka."
    ),
    central = list(
      analyte = "zabedosertib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "UNBOUND zabedosertib in the central compartment (Figure S1 'Unbound' box). central / vc is the unbound plasma concentration; clearance acts on this pool only."
    ),
    complex_agp = list(
      analyte = "zabedosertib bound to the plasma binding partner ('AGP')",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Figure S1 'Bound' box, inside the same central compartment of volume Vunb. Carried as zabedosertib mass; the molar binding-partner concentration is computed in model(). Not eliminated directly: bound drug returns to the unbound pool at koff."
    ),
    occ_irak4 = list(
      analyte = "IRAK4 target occupancy",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE,
      notes = "Fraction (0-1) of IRAK4 bound by zabedosertib, driven by the unbound plasma concentration with second-order association and first-order dissociation, target in excess-drug conditions (Section 2.11). Not observed; a projection only."
    )
  )

  covariateData <- list(
    STUDY_NCT03493269 = list(
      description = "1 = participant from the multiple-ascending-dose study NCT03493269, 0 = single-ascending-dose study NCT03054402",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SAD study NCT03054402)",
      notes = "Selects the additive residual error only: Supplementary Table S4 fixes separate additive residual variances for the SAD (0.25) and MAD (2.5) studies, which used bioanalytical methods with different lower limits of quantitation (1.0 and 10.0 ug/L, Supplementary Table S2). Acts on no structural parameter. For simulation, either value is negligible against the observed concentrations.",
      source_name = "study"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 2L,
    age_range = "18-50 years (inclusion criterion; SAD 21-50, MAD 20-50)",
    age_mean = "38 years (SAD), 35.8 years (MAD)",
    bmi_range = "19.4-30.1 kg/m^2",
    sex_female_pct = 0,
    race_ethnicity = "Mainly White (SAD: one Black/African American and one Native Hawaiian/Pacific Islander participant); MAD all White",
    disease_state = "Healthy adult men",
    dose_range = "SAD: single oral tablet doses 15-480 mg fasted; MAD: 45 or 120 mg QD and 60, 120 or 200 mg BID for 10 days after a light meal, plus a single-dose period",
    regions = "Germany (Berlin)",
    notes = "PopPK data set (Supplementary Section 1.5, Table S3): 45 SAD participants (1065 samples) and 32 MAD participants (1435 samples), 2500 concentrations in total. SAD data after the liquid formulation or in the fed state were excluded. Demographics from main Table 1 (no body weight reported; BMI only). No covariate analysis was performed."
  )

  ini({
    # Structural PK parameters (Supplementary Table S4)
    lka <- log(1.2); label("Absorption rate constant, depot -> transit1 and transit1 -> central (1/h)") # Table S4 'ka' = 1.2 1/h (RSE 9.17%)
    lcl <- log(15.4); label("Apparent unbound clearance Cl_unb (L/h)") # Table S4 'Cl unb' = 15.4 L/h (RSE 12.3%)
    lvc <- log(7.84); label("Apparent unbound volume of distribution V_unb (L)") # Table S4 'V unb' = 7.84 L (RSE 8.36%)

    # Capacity-limited binding to the plasma binding partner (Supplementary Table S4)
    lkd_pb <- fixed(log(0.4)); label("Equilibrium dissociation constant for the binding partner Kd (uM)") # Table S4 'K d' = 0.4 uM, fixed to the in vitro AGP-binding affinity
    lkoff_pb <- log(0.058); label("Dissociation rate constant koff of the drug-binding-partner complex (1/h)") # Table S4 'k off' = 0.058 1/h (RSE 8.71%)
    lbmax_pb <- log(14.9); label("Binding-partner ('AGP') concentration in plasma, i.e. binding capacity (uM)") # Table S4 'AGP' = 14.9 uM (RSE 8.00%)

    # Dose-dependent relative bioavailability (Supplementary Section 1.5 F1 equation)
    ld50_fdepot <- log(204); label("Dose at which the relative bioavailability is reduced by 50% (mg)") # Table S4 'd50' = 204 mg (RSE 11.7%)
    fmin_fdepot <- fixed(0.01); label("Lowest relative bioavailability at an infinitely high dose (fraction)") # Table S4 'f min' = 0.01, fixed

    # IRAK4 target-occupancy projection (Sections 2.11 and 3.5)
    lkd_irak4 <- fixed(log(0.086)); label("IRAK4 binding Kd, set equal to the whole-blood IL-6 release IC50 (uM)") # Section 3.5: IC50 86 (60-122) nM, assumed to give 50% IRAK4 occupancy
    lkoff_irak4 <- fixed(log(60 / 5)); label("IRAK4 dissociation rate constant, 1 / target residence time (1/h)") # Section 2.11: target residence time ~5 min; koff = 1/(5/60 h) = 12 1/h

    # IIV (Supplementary Table S4, omega^2)
    etalcl ~ 0.293 # Table S4 'Cl unb (omega2)' = 0.293
    etalvc ~ 0.071 # Table S4 'V unb (omega2)' = 0.071
    etalka ~ 0.39 # Table S4 'ka (omega2)' = 0.39
    etalbmax_pb ~ 0.124 # Table S4 'AGP (omega2)' = 0.124

    # Residual error (Supplementary Table S4, sigma^2)
    propSd <- 0.19235; label("Proportional residual error (fraction)") # Table S4 sigma2 prop = 0.037 -> sqrt(0.037) = 0.19235
    addSd <- fixed(0.0005); label("Additive residual error, SAD study NCT03054402 (mg/L)") # Table S4 sigma2 add,SAD = 0.25 fixed -> SD 0.5 ug/L = 0.0005 mg/L
    addSd_mad <- fixed(0.0015811); label("Additive residual error, MAD study NCT03493269 (mg/L)") # Table S4 sigma2 add,MAD = 2.5 fixed -> SD 1.5811 ug/L = 0.0015811 mg/L
  })

  model({
    # Molecular weight of zabedosertib (C20H21F3N4O4S), used to convert the
    # mg/L concentrations to the uM scale of Kd and the binding capacity.
    mw <- 470.5

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    kd_pb <- exp(lkd_pb)
    koff_pb <- exp(lkoff_pb)
    kon_pb <- koff_pb / kd_pb
    bmax_pb <- exp(lbmax_pb + etalbmax_pb)

    d50_fdepot <- exp(ld50_fdepot)

    kd_irak4 <- exp(lkd_irak4)
    koff_irak4 <- exp(lkoff_irak4)
    kon_irak4 <- koff_irak4 / kd_irak4

    # Unbound and bound concentrations in uM (mg/L / (g/mol) * 1000)
    cu_um <- central / vc / mw * 1000
    cb_um <- complex_agp / vc / mw * 1000

    # Net binding flux in uM/h, converted to mg/h for the amount states
    bind_um <- kon_pb * cu_um * (bmax_pb - cb_um) - koff_pb * cb_um
    bind_mg <- bind_um * vc * mw / 1000

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - bind_mg
    d/dt(complex_agp) <- bind_mg
    d/dt(occ_irak4) <- kon_irak4 * cu_um * (1 - occ_irak4) - koff_irak4 * occ_irak4

    # Relative bioavailability falls with the per-administration dose amount;
    # podo(depot) is the unscaled amount of the current dose.
    f(depot) <- fmin_fdepot + (1 - fmin_fdepot) * exp(-log(2) / d50_fdepot * podo(depot))

    Cu <- central / vc
    Cc <- (central + complex_agp) / vc

    addSd_cc <- addSd * (1 - STUDY_NCT03493269) + addSd_mad * STUDY_NCT03493269
    Cc ~ add(addSd_cc) + prop(propSd)
  })
}
