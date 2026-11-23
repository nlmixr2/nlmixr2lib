Li_2014_oseltamivir_amantadine_noimmune_qsp <- function() {
  description <- "QSP. Reduced mechanistic model (the '-IM' model) of influenza A virus infection in an average adult without immune responses -- a target-cell-limited model -- driven by oral oseltamivir phosphate (OP) and oral amantadine (AMT). The PK modules are those of the full model: each drug has a one-compartment PK module with a first-order dissolution delay and first-order absorption, and OP is converted first-order to oseltamivir carboxylate (OC). The viral life cycle runs uninfected lung epithelial cells -> eclipse-phase cells -> virus-producing cells -> free virus, and is inhibited at the replication step by AMT (1/(1 + C_AMT/IC50_AMT)) and at the virus-release step by OC (1/(1 + C_OC/IC50_OC)). The infection ends only when target cells are depleted. Li 2014 fitted this model separately from the full model, so its viral-dynamic parameters and drug potencies differ from those of modellib('Li_2014_oseltamivir_amantadine_qsp'). Deterministic typical-value ODE model: no between-subject variability and no residual error. Time 0 is the initial virus exposure; the paper starts treatment 96 h later."
  reference <- paste(
    "Li Z, Zhou H, Lu Y, Colatsky T. (2014). A Critical Role for Immune System",
    "Response in Mediating Anti-influenza Drug Synergies Assessed by Mechanistic",
    "Modeling. CPT Pharmacometrics Syst Pharmacol 3(9):e135.",
    "doi:10.1038/psp.2014.32. Equations from Supplementary Table S1 (model",
    "without immune systems), parameter values from Supplementary Table S2",
    "(-IM value, the second of each '+IM; -IM' pair) and initial values from",
    "Supplementary Table S3."
  )
  vignette <- "Li_2014_oseltamivir_amantadine_influenza"

  # Oral oseltamivir phosphate is dosed into `depot` and oral amantadine into
  # `depot_amantadine`.
  dosing <- c("depot", "depot_amantadine")

  # Viral life-cycle state with no canonical analogue in
  # inst/references/compartment-names.md (`target`, `infected` and `virus` are
  # canonical); the amantadine PK states carry a sibling-drug suffix that is
  # not registered as a metabolite suffix.
  paper_specific_compartments <- c(
    "depot_amantadine",
    "transit1_amantadine",
    "central_amantadine",
    "eclipse"
  )

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL for the plasma drug concentrations Cc, Cc_oselcarb and Cc_amantadine; TCID50/mL for virus; cells for the epithelial-cell states"
  )

  # verified = TRUE -- each state checked against Li 2014 Supplementary
  # Table S1 (equations), Table S3 (initial values and units) and Figure 1.
  compartmentData <- list(
    depot = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "plasma", verified = TRUE),
    central_oselcarb = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "plasma", verified = TRUE),
    depot_amantadine = list(analyte = "amantadine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1_amantadine = list(analyte = "amantadine", units = "mg", specimen = "administration site", verified = TRUE),
    central_amantadine = list(analyte = "amantadine", units = "mg", specimen = "plasma", verified = TRUE),
    target = list(
      analyte = "uninfected susceptible lung epithelial cells (E1)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    eclipse = list(
      analyte = "lung epithelial cells carrying influenza A virus particles, eclipse phase (E3)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    infected = list(
      analyte = "lung epithelial cells with replicated influenza A virus, virus-producing (E4)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(analyte = "free influenza A virus", units = "TCID50/mL", specimen = "tissue", verified = TRUE)
  )

  # A typical-value model for an "average" adult; no subject-level covariates.
  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = "Typical ('average') adult calibrated against published summary data; no individual-level data",
    n_studies = 6L,
    disease_state = "Acute uncomplicated influenza A virus infection (natural); healthy volunteers for the PK modules",
    dose_range = "Oseltamivir phosphate 75-150 mg twice daily and amantadine 100 mg twice daily in the calibration data; 0-400 mg twice daily of each in the combination simulations",
    notes = paste(
      "Li 2014 calibrated the model to published summary data rather than fitting individual data.",
      "PK modules (shared with the full model): oseltamivir loading-dose (150 + 75 mg) multidose study in",
      "Thai volunteers (Wattanagoon 2009), 150 mg single dose in Japanese and Caucasian volunteers",
      "(Jhee 2008), amantadine 100 mg once daily (Reuman 1989) and 200 mg single dose averaged over young",
      "and elderly adults (Hayden 1985) -- Figure 2. Treatment calibration: oseltamivir 75 and 150 mg twice",
      "daily versus placebo (Treanor 2000) and amantadine versus placebo (Knight 1970) -- Figure 4.",
      "Parameters were estimated with a genetic algorithm (GAPSM) within literature-based ranges",
      "(Supplementary Table S2)."
    )
  )

  ini({
    # ==================================================================
    # PK modules -- Li 2014 Supplementary Table S2 (shared by the +IM and
    # -IM models). Clearances and volumes are tabulated in mL/h and mL for
    # a 70-kg adult and are kept in those units; model() converts the
    # mg/mL ratio of a dose in mg to ng/mL.
    # ==================================================================
    lktr <- log(6.412)
    label("OP dissolution (delay) rate constant from depot to transit1, K1a (1/h)") # Table S2: K1a = 6.412
    lka <- log(0.3652)
    label("OP absorption rate constant from transit1 to central, K1b (1/h)") # Table S2: K1b = 0.3652
    lkm <- log(1.768)
    label("OP to OC first-order conversion rate constant, K2 (1/h)") # Table S2: K2 = 1.768
    lcl <- log(78060)
    label("Apparent OP clearance, CLOP (mL/h)") # Table S2: CLOP = 78060
    lvc <- log(212600)
    label("Apparent OP volume of distribution, VDOP (mL)") # Table S2: VDOP = 212600
    lcl_oselcarb <- log(26460)
    label("Apparent OC clearance, CLOC (mL/h)") # Table S2: CLOC = 26460
    lvc_oselcarb <- log(115000)
    label("Apparent OC volume of distribution, VDOC (mL)") # Table S2: VDOC = 115000
    lktr_amantadine <- log(62.69)
    label("AMT dissolution (delay) rate constant from depot to transit1, K19 (1/h)") # Table S2: K19 = 62.69
    lka_amantadine <- log(3.356)
    label("AMT absorption rate constant from transit1 to central, K20 (1/h)") # Table S2: K20 = 3.356
    lcl_amantadine <- log(12850)
    label("Apparent AMT clearance, CLAMT (mL/h)") # Table S2: CLAMT = 12850
    lvc_amantadine <- log(325900)
    label("Apparent AMT volume of distribution, VDAMT (mL)") # Table S2: VDAMT = 325900

    # ==================================================================
    # Drug effects -- Li 2014 Supplementary Table S2 (-IM value, the second
    # of each '+IM; -IM' pair). 1/IC50 is K17 for OC and K22 for AMT.
    # ==================================================================
    lic50_oselcarb <- log(1 / 0.0237)
    label("In vivo plasma OC IC50 for inhibition of virus release, 1/K17 (ng/mL)") # Table S2: K17 = 0.0237 (-IM)
    lic50_amantadine <- log(1 / 0.25710)
    label("In vivo plasma AMT IC50 for inhibition of viral replication, 1/K22 (ng/mL)") # Table S2: K22 = 0.25710 (-IM)

    # ==================================================================
    # Virus life-cycle module -- Li 2014 Supplementary Table S2, -IM value.
    # ==================================================================
    lbeta <- log(4.74e-6)
    label("Epithelial cell infection rate constant, K11 (1/(TCID50/mL)/h)") # Table S2: K11 = 4.74*10-6 (-IM)
    lk_eclipse <- log(0.11661)
    label("Eclipse-phase transition rate constant to virus-producing cells, K13 (1/h)") # Table S2: K13 = 0.11661 (-IM)
    ldelta <- log(0.35438)
    label("Death rate constant of virus-producing cells, K14 (1/h)") # Table S2: K14 = 0.35438 (-IM)
    lp_virus <- log(0.000933)
    label("Virus release rate per virus-producing cell, K16 ((TCID50/mL)/cell/h)") # Table S2: K16 = 0.000933 (-IM)
    lc_virus <- log(0.31658)
    label("Free-virus elimination rate constant, K18 (1/h)") # Table S2: K18 = 0.31658 (-IM)

    target0 <- fixed(4e8)
    label("Initial uninfected lung epithelial cells (cells)") # Table S3: EpitheliumLung = 4*108
    # Initial free-virus titre. Li 2014 fitted it separately for each
    # simulated data set (Table S3, -IM values: 0.03996, 0.02006, 0.1098
    # for the Figure 4a placebo, 75-mg and 150-mg groups; 5.686*10-4 for the
    # Figure 4b placebo group and Figure 4c; 0.001157 for the Figure 4b
    # amantadine group) and fixed it at 10 for the drug-combination
    # simulations.
    virus0 <- fixed(10)
    label("Initial free influenza A virus titre at exposure (TCID50/mL)") # Methods 'Calculation of synergy': initial titre of 10
  })

  model({
    ktr <- exp(lktr)
    ka <- exp(lka)
    km <- exp(lkm)
    cl <- exp(lcl)
    vc <- exp(lvc)
    cl_oselcarb <- exp(lcl_oselcarb)
    vc_oselcarb <- exp(lvc_oselcarb)
    ktr_amantadine <- exp(lktr_amantadine)
    ka_amantadine <- exp(lka_amantadine)
    cl_amantadine <- exp(lcl_amantadine)
    vc_amantadine <- exp(lvc_amantadine)
    ic50_oselcarb <- exp(lic50_oselcarb)
    ic50_amantadine <- exp(lic50_amantadine)
    beta <- exp(lbeta)
    k_eclipse <- exp(lk_eclipse)
    delta <- exp(ldelta)
    p_virus <- exp(lp_virus)
    c_virus <- exp(lc_virus)

    # ------------------------------------------------------------------
    # PK modules -- Table S1. The OP to OC conversion carries no
    # molecular-weight correction, as printed.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - km * central - cl / vc * central
    d/dt(central_oselcarb) <- km * central - cl_oselcarb / vc_oselcarb * central_oselcarb
    d/dt(depot_amantadine) <- -ktr_amantadine * depot_amantadine
    d/dt(transit1_amantadine) <- ktr_amantadine * depot_amantadine - ka_amantadine * transit1_amantadine
    d/dt(central_amantadine) <- ka_amantadine * transit1_amantadine - cl_amantadine / vc_amantadine * central_amantadine

    # Amounts in mg over volumes in mL give mg/mL; x 1e6 gives ng/mL, the
    # units of the 1/IC50 parameters K17 and K22.
    Cc <- central / vc * 1e6
    Cc_oselcarb <- central_oselcarb / vc_oselcarb * 1e6
    Cc_amantadine <- central_amantadine / vc_amantadine * 1e6

    inh_replication <- 1 / (1 + Cc_amantadine / ic50_amantadine)
    inh_release <- 1 / (1 + Cc_oselcarb / ic50_oselcarb)

    # ------------------------------------------------------------------
    # Virus life-cycle module -- Table S1 (model without immune systems).
    # ------------------------------------------------------------------
    d/dt(target) <- -beta * virus * target
    d/dt(eclipse) <- beta * virus * target - inh_replication * k_eclipse * eclipse
    d/dt(infected) <- inh_replication * k_eclipse * eclipse - delta * infected
    d/dt(virus) <- inh_release * p_virus * infected - c_virus * virus

    target(0) <- target0
    virus(0) <- virus0

    # Li 2014 plots titres as log10(titre + 1) (Table S3 note).
    log10_viral_load <- log10(virus + 1)
    # Fraction of epithelial cells still uninfected, the dashed curve of
    # Figure 4c.
    frac_uninfected <- target / target0
  })
}
