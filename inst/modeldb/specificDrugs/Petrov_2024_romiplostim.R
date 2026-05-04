Petrov_2024_romiplostim <- function() {
  description <- "Population PK/PD model for romiplostim in adults with chronic immune thrombocytopenia (ITP). One-compartment first-order subcutaneous PK plus an Emax stimulation of platelet precursor production into a 4-transit-compartment Friberg-style chain feeding circulating platelets, with first-order platelet degradation. PK/PD backbone is the healthy-volunteer population PK/PD model (Makarenko 2024); ITP-specific platelet production (kin) and degradation (kdeg) constants and IIV(kdeg) come from Petrov 2024 supplement Table S1. Default parameters are non-splenectomized ITP patients with mechanism 1 (increased platelet degradation, normal precursor production); see vignette for the other 3 subpopulation variants (non-splenectomized mechanism 2; splenectomized mechanism 1; splenectomized mechanism 2)."
  reference <- paste(
    "Petrov A, Makarenko I, Sokolov V, Drai R, Bondareva I, Sigaev V, Stuchkov M, Galustyan A, Stepanenko I, Lebedev V, Mishchenko A. Optimization of Romiplostim Biosimilar Efficacy Trial Using In Silico Clinical Trial Approach for Patients With Immune Thrombocytopenia. Clin Pharmacol Drug Dev. 2025 Feb;14(2):116-126. doi:10.1002/cpdd.1494 (PMID 39702972).",
    "PK/PD backbone (healthy volunteers):",
    "Makarenko I, Petrov A, Sokolov V, Drai R, Mishchenko A, Bondareva I, Galustyan A, Sigaev V. Population Pharmacokinetic and Pharmacodynamic Modeling of Romiplostim Biosimilar GP40141 and Reference Product in Healthy Volunteers to Evaluate Biosimilarity. Clin Pharmacol Drug Dev. 2024. doi:10.1002/cpdd.1367 (PMID 38168134; reference 20 in Petrov 2024)."
  )
  vignette <- "Petrov_2024_romiplostim"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL", platelet = "10^9/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power covariate on apparent volume of distribution V/F: V = TVV * (WT/77)^1.04 (Petrov 2024 supplement Table S1; Monolix 2021R1 default form for a continuous covariate on a log-normal parameter). Reference 77 kg = population mean weight per Petrov 2024 Methods.",
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug (neutralizing) antibody status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NAB-negative)",
      notes              = "Log-additive covariate on apparent first-order elimination constant: kel = TVKEL * exp(0.25 * ADA_POS) (Petrov 2024 supplement Table S1; Monolix 2021R1 default form for a binary covariate on a log-normal parameter). NAB+ subjects therefore have ~28% higher kel (faster elimination via neutralizing-antibody-mediated clearance). The source paper uses 'NAB' (neutralizing antibody, as opposed to total ADA); the canonical column name in nlmixr2lib is ADA_POS. For datasets that record only total ADA, the effect strictly applies to the NAB-positive subset.",
      source_name        = "NAB"
    )
  )

  population <- list(
    n_subjects     = 83L,
    n_studies      = 2L,
    age_range      = ">= 18 years",
    weight_mean    = "77 kg (CV 20%, normally distributed in the simulated cohort per Petrov 2024 Methods)",
    sex_female_pct = NA,
    race_ethnicity = NULL,
    disease_state  = "Adults with chronic immune thrombocytopenia (ITP). Default parameter set: non-splenectomized, ITP mechanism 1 (increased platelet degradation, normal precursor production). Three additional subpopulations exist (non-splenectomized mechanism 2 with reduced precursor production; splenectomized mechanism 1; splenectomized mechanism 2) with different (kin, kdeg, IIV(kdeg)) tuples - documented in the vignette.",
    dose_range     = "Subcutaneous romiplostim with weekly dose-titration algorithm targeting platelets 50-200 x 10^9/L (validation simulation, max 15 ug/kg) or 50-150 x 10^9/L (biosimilar efficacy ISCT, max 10 ug/kg). Initial dose 1 ug/kg.",
    nab_pos_pct    = 5,
    regions        = "Russian Federation (clinical development of GP40141 biosimilar, Makarenko 2024 healthy-volunteer source data). Validation against international romiplostim-ref pivotal trials NCT00102323 (splenectomized adults) and NCT00102336 (non-splenectomized adults).",
    notes          = "n_subjects = 83 reflects the validation cohort (41 non-splenectomized + 42 splenectomized adults from Petrov 2024 Methods / supplement Covariate distribution). PK/PD backbone parameter values are inherited from Makarenko 2024 (healthy-volunteer popPK/PD analysis cited as reference 20 in Petrov 2024); the ITP modifications scale kin and kdeg using published platelet kinetic data (Ballem 1987 and Stoll 1985, references 31 and 32 in Petrov 2024) and were tuned to match the baseline platelet distribution of the romiplostim-ref pivotal trial (Kuter 2008, reference 13 in Petrov 2024). Steady-state baseline platelet count (simulated, 1000 virtual subjects without drug): mean ~18 x 10^9/L (CV 35%) for non-splenectomized; mean ~14 x 10^9/L (CV 50%) for splenectomized."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK parameters - Petrov 2024 supplement Table S1.
    # Footnote 'a' indicates 'parameters identical to healthy subjects',
    # i.e. inherited unchanged from Makarenko 2024 (reference 20 in
    # Petrov 2024). Reference subject is 77 kg, NAB-negative.
    # ---------------------------------------------------------------------
    lka  <- log(0.02);   label("First-order absorption rate ka (1/h)")                                            # Petrov 2024 Table S1
    lvc  <- log(2565);   label("Apparent central volume of distribution V/F (L) at reference WT 77 kg")           # Petrov 2024 Table S1
    lkel <- log(0.03);   label("First-order apparent elimination rate kel (1/h) for NAB-negative subject")        # Petrov 2024 Table S1

    # ---------------------------------------------------------------------
    # PD - drug stimulation - Petrov 2024 supplement Table S1.
    # EC50 originally 42 pg/mL; here converted to ng/mL (= 0.042 ng/mL)
    # so the same Cc = central(ug)/V(L) = ng/mL feeds the Emax expression
    # without intermediate unit conversion.
    # ---------------------------------------------------------------------
    lec50 <- log(0.042); label("Half-maximal effective concentration EC50 (ng/mL); Table S1 reports 42 pg/mL")    # Petrov 2024 Table S1
    lemax <- log(9);     label("Maximum stimulatory effect Emax on platelet precursor production (unitless)")     # Petrov 2024 Table S1

    # ---------------------------------------------------------------------
    # PD - platelet kinetics - Petrov 2024 supplement Table S1.
    # Default subpopulation = non-splenectomized, ITP mechanism 1
    # (increased platelet degradation only; normal precursor production).
    # The other 3 (kin, kdeg, IIV_kdeg) tuples are listed in the vignette.
    # ---------------------------------------------------------------------
    lktr  <- log(0.02);  label("Platelet precursor first-order transit constant ktr (1/h)")                       # Petrov 2024 Table S1
    lkin  <- log(4.3);   label("Platelet precursor production rate constant kin (10^9 cells/L/h) - non-splen mech 1")  # Petrov 2024 Table S1
    lkdeg <- log(0.20);  label("Platelet first-order degradation constant kdeg (1/h) - non-splen mech 1")         # Petrov 2024 Table S1

    # ---------------------------------------------------------------------
    # Covariate effects - Petrov 2024 supplement Table S1.
    # No explicit functional form is printed in Table S1; per Petrov 2024
    # Methods the model was fit in Monolix 2021R1, whose default
    # parameterization for log-normal parameters is allometric power on a
    # mean-centered continuous covariate and log-additive on a categorical
    # covariate. These forms reproduce the magnitude of the reported
    # effects: V/F scales near-linearly with body weight and NAB+ subjects
    # have ~28% higher kel.
    # ---------------------------------------------------------------------
    e_wt_vc   <- 1.04;  label("Allometric power exponent of (WT/77 kg) on Vc/F (unitless)")                       # Petrov 2024 Table S1
    e_nab_kel <- 0.25;  label("Log-additive effect of NAB+ on kel: kel_NABpos = kel * exp(0.25)")                 # Petrov 2024 Table S1

    # ---------------------------------------------------------------------
    # Inter-individual variability - Petrov 2024 supplement Table S1.
    # IIV reported as %CV (coefficient of variation, log-normal parameters);
    # internal omega^2 = log(CV^2 + 1).
    # No IIV on ka (Table S1: '-').
    # Diagonal Omega: no inter-parameter correlations are reported.
    # ---------------------------------------------------------------------
    etalvc   ~ 0.08618; label("IIV variance on log-V/F (Petrov 2024 Table S1: 30% CV)")                          # log(0.30^2 + 1) = 0.08618
    etalkel  ~ 0.11556; label("IIV variance on log-kel (Petrov 2024 Table S1: 35% CV)")                          # log(0.35^2 + 1) = 0.11556
    etalec50 ~ 0.42716; label("IIV variance on log-EC50 (Petrov 2024 Table S1: 73% CV)")                         # log(0.73^2 + 1) = 0.42716
    etalemax ~ 0.02225; label("IIV variance on log-Emax (Petrov 2024 Table S1: 15% CV)")                         # log(0.15^2 + 1) = 0.02225
    etalktr  ~ 0.01203; label("IIV variance on log-ktr (Petrov 2024 Table S1: 11% CV)")                          # log(0.11^2 + 1) = 0.01203
    etalkin  ~ 0.01941; label("IIV variance on log-kin (Petrov 2024 Table S1: 14% CV)")                          # log(0.14^2 + 1) = 0.01941
    etalkdeg ~ 0.22314; label("IIV variance on log-kdeg (Petrov 2024 Table S1: 50% CV) - non-splen")             # log(0.50^2 + 1) = 0.22314

    # ---------------------------------------------------------------------
    # Residual error - Petrov 2024 supplement Table S1 (b = 0.093).
    # Table S1 lists a single 'Proportional error model parameter' without
    # naming the output. The only modeled observation in the validation
    # and ISCT is platelet count (no concentration data are simulated as
    # observations), so b is interpreted here as the proportional residual
    # error on platelet count. See vignette Errata for the ambiguity.
    # ---------------------------------------------------------------------
    propSd_PLT <- 0.093; label("Proportional residual error on platelet count (fraction)")                        # Petrov 2024 Table S1
  })

  model({
    # Individual PK/PD parameters
    ka   <- exp(lka)
    vc   <- exp(lvc + etalvc)   * (WT / 77)^e_wt_vc
    kel  <- exp(lkel + etalkel) * exp(e_nab_kel * ADA_POS)
    ec50 <- exp(lec50 + etalec50)
    emax <- exp(lemax + etalemax)
    ktr  <- exp(lktr  + etalktr)
    kin  <- exp(lkin  + etalkin)
    kdeg <- exp(lkdeg + etalkdeg)

    # Concentration in ng/mL (= ug/L). Doses enter `depot` in micrograms;
    # vc is in litres so central(ug)/vc(L) = ug/L = ng/mL.
    Cc <- central / vc

    # Drug-driven stimulation of platelet precursor production
    stim <- emax * Cc / (ec50 + Cc)

    # 1-compartment SC PK + Friberg-style 4-transit precursor chain feeding
    # circulating platelets, with first-order platelet degradation.
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  ka * depot - kel * central
    d/dt(precursor1) <-  kin * (1 + stim) - ktr * precursor1
    d/dt(precursor2) <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <-  ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <-  ktr * precursor4 - kdeg * circ

    # Steady-state baseline (drug-free): each precursor pool at kin/ktr,
    # circulating platelets at kin/kdeg. Each subject's individual baseline
    # is set by their (kin, kdeg) draws.
    precursor1(0) <- kin / ktr
    precursor2(0) <- kin / ktr
    precursor3(0) <- kin / ktr
    precursor4(0) <- kin / ktr
    circ(0)       <- kin / kdeg

    # Observation: circulating platelet count in 10^9 cells/L
    PLT <- circ
    PLT ~ prop(propSd_PLT)
  })
}
