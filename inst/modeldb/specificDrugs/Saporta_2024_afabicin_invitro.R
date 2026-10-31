Saporta_2024_afabicin_invitro <- function() {
  description <- "In vitro (Staphylococcus aureus, 21 strains). Semi-mechanistic time-kill pharmacodynamic model of afabicin desphosphono, the active moiety of the FabI-inhibitor prodrug afabicin, fitted to 162 static time-kill curves. Bacteria are either growing and drug-susceptible (S) or dormant and drug-insusceptible (D); S enters D at kSD = (S + D) * (kgrowth - kdeath) / Bmax. Afabicin desphosphono adds a sigmoid Emax killing rate kdrug to S, with EC50 scaled by the strain MIC as EC50 = beta * MIC^gamma (beta fixed to 1). A reversible adaptation pair (aroff -> aron at kon * C, back at koff) reduces Emax as Emax * (1 - aron / (aron + AS50)). No inter-experiment variability; additive residual error on log10 CFU/mL. Sibling model: Saporta_2024_afabicin_invitro_mouse (the same structure re-estimated jointly on these in vitro data and neutropenic mouse thigh infection data)."
  reference <- paste(
    "Saporta R, Nielsen EI, Menetrey A, Cameron DR, Nicolas-Metral V,",
    "Friberg LE. Model-based translation of results from in vitro to in vivo",
    "experiments for afabicin activity against Staphylococcus aureus. J",
    "Antimicrob Chemother. 2024;79(12):3150-3159. doi:10.1093/jac/dkae334.",
    "PMID: 39315768. PMCID: PMC11638087. Model structure: Materials and",
    "methods, Eqs 1-6, and Figure 1. Parameter estimates: Table 1, column",
    "'PK/PD model based on in vitro data'. Strain MICs: Supplementary Table",
    "S1. kdeath fixed to 0.179 per hour from Nielsen EI et al., Antimicrob",
    "Agents Chemother 2007;51:128-136; koff fixed to 0.0139 per hour from",
    "Mohamed AF et al., Antimicrob Agents Chemother 2012;56:179-188.",
    sep = " "
  )
  vignette <- "Saporta_2024_afabicin"
  units <- list(time = "h", dosing = "mg/L (bath concentration)", concentration = "mg/L")

  # Declared explicitly because the dosing target is neither `depot` nor
  # `central`: this is a static in vitro system, so the "dose" is the bath
  # concentration placed directly into the `afabicin_desphosphono` state at
  # time 0.
  dosing <- c("afabicin_desphosphono")

  # S and D (Figure 1) follow the bact_<state> naming of the sibling
  # Nielsen-lineage files (HernandezLozano_2025_apramycin_invitro.R);
  # `afabicin_desphosphono` holds the static bath concentration in mg/L.
  # `aroff` / `aron` are the registered adaptive-resistance canonicals for
  # the paper's ASoff / ASon.
  paper_specific_compartments <- c(
    "afabicin_desphosphono",
    "bact_s",
    "bact_d"
  )

  compartmentData <- list(
    afabicin_desphosphono = list(
      analyte = "afabicin desphosphono",
      units = "mg/L",
      specimen = "administration site",
      verified = TRUE
    ),
    bact_s = list(
      analyte = "Staphylococcus aureus, growing drug-susceptible (S) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_d = list(
      analyte = "Staphylococcus aureus, dormant drug-insusceptible (D) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff = list(
      analyte = "fraction of bacteria without adaptation (ASoff)",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron = list(
      analyte = "fraction of bacteria with adaptive reduction of susceptibility (ASon)",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Staphylococcus aureus, 21 strains)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Staphylococcus aureus, 21 strains (12 MRSA, plus MSSA, VISA and linezolid-resistant isolates; Supplementary Table S1), afabicin desphosphono MIC 0.004 to 0.03 mg/L",
    system = "Static time-kill experiments in Mueller-Hinton Broth II at 35 C in 24-well plates (2 mL, 24 curves), 15 mL tubes (10 mL, 44 curves) or 125 mL Erlenmeyer flasks (20 mL, 94 curves), sampled up to 24 or 48 h",
    n_curves = 162L,
    n_observations = 722L,
    starting_inoculum = "target 1.6e5 to 2.0e7 CFU/mL",
    concentration_range = "afabicin desphosphono 0.004 to 1 mg/L, plus drug-free growth controls",
    limit_of_detection = "200 CFU/mL (8 of 722 counts below it)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "MICs were determined by broth microdilution according to CLSI M07; where an MIC was measured during a time-kill experiment that value was used, so several strains carry a range of MICs across experiments (Supplementary Table S1).",
      "The growth-control experiments were fitted first to estimate kgrowth and Bmax, which were then fixed while the drug-effect parameters were estimated on all time-kill experiments.",
      "Estimation in NONMEM 7.5.0 (FOCE with interaction) on log10-transformed counts with an additive error; inter-individual variability was not explored."
    )
  )

  ini({
    # ---- Simulation design inputs ---------------------------------------
    # The strain MIC scales EC50 (Eq 6). It is a fixed ini() entry, as in
    # HernandezLozano_2025_apramycin_mouse.R, so a user re-points the model
    # at another isolate with rxSolve(params = c(mic = ...)) or a per-id
    # params data frame. Default: S. aureus ATCC 29213, the strain with the
    # richest time-kill data.
    mic <- fixed(0.008); label("Afabicin desphosphono MIC of the simulated strain (mg/L)")         # Supplementary Table S1: ATCC 29213 MIC 0.004 - 0.008 mg/L; Results quote EC50 at MIC 0.008 mg/L
    # Starting inoculum is set by the experiment, not estimated. Default is
    # the round value nearest the median starting inoculum drawn in Figure 2
    # (upper broken line, about 6 log10 CFU/mL).
    linoc <- fixed(log(1e6)); label("Log starting bacterial inoculum (CFU/mL; design input)")      # Methods: target inocula 1.6e5 to 2.0e7 CFU/mL; Figure 2 median starting inoculum about 6 log10 CFU/mL

    # ---- Bacterial growth (Table 1, in vitro column) --------------------
    lkgrowth <- fixed(log(1.13)); label("Log bacterial growth rate constant kgrowth (1/h)")        # Table 1: kgrowth = 1.13 1/h (RSE 2.0%), footnote a: fixed from the growth-control fit
    lkdeath <- fixed(log(0.179)); label("Log bacterial natural death rate constant kdeath (1/h)")  # Table 1: kdeath = 0.179 1/h (Fixed); Results: literature value (Nielsen 2007, ref 40)
    lbmax <- fixed(log(10^9.85)); label("Log maximum system capacity Bmax (CFU/mL)")               # Table 1: Bmax vitro = 9.85 log10 CFU/mL (RSE 1.0%), footnote a: fixed from the growth-control fit

    # ---- Drug effect (Table 1, in vitro column) -------------------------
    lemax <- log(4.10); label("Log maximum afabicin desphosphono killing rate constant Emax (1/h)") # Table 1: Emax = 4.10 1/h (RSE 7.0%)
    # Eq 6: EC50 = beta * MIC^gamma. beta is therefore the EC50 at an MIC
    # of 1 mg/L, carried as the canonical lec50; gamma is the power of the
    # MIC effect on EC50.
    lec50 <- fixed(log(1)); label("Log EC50 scaling slope beta, the EC50 at an MIC of 1 mg/L (mg/L)") # Table 1: beta vitro = 1 (Fixed)
    e_mic_ec50 <- 0.160; label("Power of the MIC effect on EC50, gamma (unitless)")                # Table 1: gamma = 0.160 (RSE 36%)
    lhill <- log(0.360); label("Log sigmoidicity (Hill) factor of the Emax model (unitless)")       # Table 1: Hill = 0.360 (RSE 6.4%)

    # ---- Adaptation (Table 1, in vitro column; Eqs 3-5) -----------------
    lkon <- log(0.406); label("Log adaptation onset rate constant kon (L/(mg*h))")                 # Table 1: kon = 0.406 mL/(ug*h) (RSE 24%); mL/(ug*h) = L/(mg*h)
    lkoff <- fixed(log(0.0139)); label("Log adaptation offset rate constant koff (1/h)")           # Table 1: koff = 0.0139 1/h (Fixed); Methods: half-life of 50 h, from Mohamed 2012 (ref 41)
    lar50 <- log(0.861); label("Log fraction in ASon giving a 50% reduction in Emax, AS50 (unitless)") # Table 1: AS50 = 0.861 (RSE 13%)

    # ---- Residual error --------------------------------------------------
    addSd <- 0.869; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")       # Table 1: RES vitro = 0.869 log10 CFU/mL (RSE 5.0%); footnote b: SD scale
  })

  model({
    # ---- 1. Individual parameters ---------------------------------------
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    emax <- exp(lemax)
    hill <- exp(lhill)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ar50 <- exp(lar50)

    # Eq 6, one-parameter MIC implementation (beta fixed to 1).
    ec50 <- exp(lec50) * mic^e_mic_ec50

    # ---- 2. Drug effect ---------------------------------------------------
    # Static bath concentration in protein-free broth: the state is the
    # unbound concentration C of Figure 1, in mg/L.
    conc <- afabicin_desphosphono
    # Eq 5: adaptation reduces Emax.
    emax_t <- emax * (1 - aron / (aron + ar50))
    kdrug <- emax_t * conc^hill / (ec50^hill + conc^hill)

    # ---- 3. Density-dependent transfer into the dormant state -------------
    btot <- bact_s + bact_d
    ksd <- btot * (kgrowth - kdeath) / bmax

    # ---- 4. Drug in the growth medium ---------------------------------------
    d/dt(afabicin_desphosphono) <- 0

    # ---- 5. Bacterial system (Eqs 1-2) ------------------------------------
    d/dt(bact_s) <- kgrowth * bact_s - (kdeath + kdrug) * bact_s - ksd * bact_s
    d/dt(bact_d) <- ksd * bact_s - kdeath * bact_d

    # ---- 6. Adaptation (Eqs 3-4) ------------------------------------------
    d/dt(aroff) <- koff * aron - kon * conc * aroff
    d/dt(aron) <- kon * conc * aroff - koff * aron

    # ---- 7. Initial conditions --------------------------------------------
    # Methods: all bacteria start in S; ASoff = 1 and ASon = 0 at start.
    bact_s(0) <- exp(linoc)
    aroff(0) <- 1

    # ---- 8. Observation ----------------------------------------------------
    # Log10-transform-both-sides with an additive error (Data analysis).
    log_cfu <- log10(btot)
    log_cfu ~ add(addSd)
  })
}
