Saporta_2024_afabicin_invitro_mouse <- function() {
  description <- "In vitro (Staphylococcus aureus) and preclinical (mouse, neutropenic thigh infection). Semi-mechanistic pharmacodynamic model of afabicin desphosphono, the active moiety of the FabI-inhibitor prodrug afabicin, re-estimated jointly on in vitro static time-kill data and neutropenic mouse thigh infection data. Two parallel bacterial systems share one parameter set except for Bmax, the EC50 scaling slope beta and the residual error: an in vitro system driven by a static bath concentration, and a thigh system driven by the unbound plasma concentration fu * CP_AFABICIN_DESPHOSPHONO_MGL with fu = 0.02. Each system has growing drug-susceptible bacteria (S), dormant insusceptible bacteria (D) entered at kSD = (S + D) * (kgrowth - kdeath) / Bmax, a sigmoid Emax kill on S with EC50 = beta * MIC^gamma, and a reversible adaptation pair (aroff -> aron at kon * C, back at koff) that reduces Emax as Emax * (1 - aron / (aron + AS50)). The source mouse PK model (Bader 2017, a conference poster) reports no parameter values, so the plasma concentration is a time-varying input. No inter-experiment variability; additive residual error on log10 CFU. Sibling model: Saporta_2024_afabicin_invitro (the fit to the in vitro data alone)."
  reference <- paste(
    "Saporta R, Nielsen EI, Menetrey A, Cameron DR, Nicolas-Metral V,",
    "Friberg LE. Model-based translation of results from in vitro to in vivo",
    "experiments for afabicin activity against Staphylococcus aureus. J",
    "Antimicrob Chemother. 2024;79(12):3150-3159. doi:10.1093/jac/dkae334.",
    "PMID: 39315768. PMCID: PMC11638087. Model structure: Materials and",
    "methods, Eqs 1-6, section 'PK/PD model for in vivo data', and Figure 1.",
    "Parameter estimates: Table 1, column 'PK/PD model based on in vitro and",
    "in vivo data'. Strain MICs: Supplementary Table S1. The mouse PK model",
    "that supplied the plasma concentration is Bader JC, Lakota EA, Bravo J",
    "et al., ASM Microbe 2017, New Orleans, abstract 2176, whose parameter",
    "values are unpublished.",
    sep = " "
  )
  vignette <- "Saporta_2024_afabicin"
  units <- list(time = "h", dosing = "mg/L (in vitro bath concentration)", concentration = "mg/L")

  # The only dosed state is the static in vitro bath concentration; the thigh
  # system is driven by the CP_AFABICIN_DESPHOSPHONO_MGL covariate.
  dosing <- c("afabicin_desphosphono")

  # The in vitro states keep the names used in Saporta_2024_afabicin_invitro.R;
  # the thigh states carry a `_thigh` suffix. `aroff` / `aron` are the
  # registered adaptive-resistance canonicals for the paper's ASoff / ASon;
  # their thigh copies and the two log-count outputs are whitelisted here.
  paper_specific_compartments <- c(
    "afabicin_desphosphono",
    "bact_s",
    "bact_d",
    "bact_s_thigh",
    "bact_d_thigh",
    "aroff_thigh",
    "aron_thigh",
    "log_cfu_invitro",
    "log_cfu_thigh"
  )

  compartmentData <- list(
    afabicin_desphosphono = list(
      analyte = "afabicin desphosphono",
      units = "mg/L",
      specimen = "administration site",
      verified = TRUE
    ),
    bact_s = list(
      analyte = "Staphylococcus aureus in vitro, growing drug-susceptible (S) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_d = list(
      analyte = "Staphylococcus aureus in vitro, dormant drug-insusceptible (D) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff = list(
      analyte = "in vitro fraction of bacteria without adaptation (ASoff)",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron = list(
      analyte = "in vitro fraction of bacteria with adaptive reduction of susceptibility (ASon)",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_s_thigh = list(
      analyte = "Staphylococcus aureus in the mouse thigh, growing drug-susceptible (S) state",
      units = "CFU/g",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_d_thigh = list(
      analyte = "Staphylococcus aureus in the mouse thigh, dormant drug-insusceptible (D) state",
      units = "CFU/g",
      specimen = "tissue",
      verified = TRUE
    ),
    aroff_thigh = list(
      analyte = "thigh fraction of bacteria without adaptation (ASoff)",
      units = "(fraction)",
      specimen = "tissue",
      verified = TRUE
    ),
    aron_thigh = list(
      analyte = "thigh fraction of bacteria with adaptive reduction of susceptibility (ASon)",
      units = "(fraction)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_AFABICIN_DESPHOSPHONO_MGL = list(
      description = "Total afabicin desphosphono plasma concentration in the mouse, supplied per record as the time-varying driver of the thigh bacterial system.",
      units = "mg/L",
      type = "continuous",
      source_name = "not named in the source (derived from the Bader 2017 mouse PK model)",
      reference_category = NULL,
      notes = paste(
        "Methods, 'PK/PD model for in vivo data': 'A mouse PK model for afabicin was used to derive the afabicin desphosphono concentration-time profiles for the mouse efficacy studies'; all mice were assumed to have the typical PK parameters, the same after IV or IP afabicin.",
        "That model (Bader 2017, ASM Microbe abstract 2176) is cited only as a conference abstract; Saporta 2026 (doi:10.1128/aac.00959-25) Supplementary text 2 prints its structure (IP depot -> absorption transit -> two-compartment afabicin desphosphono with an afabicin-to-desphosphono conversion factor) but no source reports its parameter values, so the concentration enters as a covariate.",
        "The model multiplies this total concentration by fu = 0.02 to obtain the unbound concentration that drives killing and adaptation. Set to 0 for growth-control animals, and for records that simulate only the in vitro system."
      )
    )
  )

  population <- list(
    species = "in vitro (Staphylococcus aureus) + mouse (female ICR or Swiss Webster, neutropenic)",
    n_subjects = 952L,
    n_studies = NA_integer_,
    organism = "Staphylococcus aureus: 21 strains in vitro (MIC 0.004 to 0.03 mg/L) and 9 strains in vivo (6 MRSA; MIC 0.004 to 0.06 mg/L), with ATCC 29213 studied in both (Supplementary Table S1)",
    system = "In vitro: 162 static time-kill curves in Mueller-Hinton Broth II at 35 C (722 counts). In vivo: neutropenic mouse thigh infection model; cyclophosphamide 150 mg/kg IP on Day -4 and 100 mg/kg on Day -1, intramuscular inoculum of 8.0e4 to 1.3e7 CFU into the left thigh on Day 0, thighs harvested at 2, 26, 50 or 74 h after infection.",
    dose_range = "Afabicin (prodrug) started 2 h after infection: a single 129 mg/kg dose, or 0.011 to 190 mg/kg every 6 h, given IV (n = 82) or IP (n = 430)",
    sex_female_pct = 100,
    disease_state = "Neutropenic thigh infection with Staphylococcus aureus (in vivo arm)",
    notes = paste(
      "n_subjects counts the mice (568 ICR, 384 Swiss Webster); the in vitro arm adds 162 time-kill curves.",
      "The bacteria-related parameters were re-estimated first on growth-control data, with kgrowth and kdeath shared between in vitro and in vivo and Bmax setting-specific, and then fixed. Drug-effect parameters were then re-estimated by stepwise unfixing of in vivo / in vitro ratios; only the ratio on beta was retained.",
      "Estimation in NONMEM 7.5.0 (FOCE with interaction) on log10-transformed counts with setting-specific additive errors; inter-individual variability was not explored."
    )
  )

  ini({
    # ---- Simulation design inputs ---------------------------------------
    # Strain MICs scale EC50 (Eq 6); fixed ini() entries, as in
    # HernandezLozano_2025_apramycin_mouse.R, so a user re-points either
    # system at another isolate with rxSolve(params = ...). Default: S. aureus
    # ATCC 29213, the one strain studied both in vitro and in vivo.
    mic <- fixed(0.008); label("Afabicin desphosphono MIC of the strain in the in vitro system (mg/L)")   # Supplementary Table S1: ATCC 29213 MIC 0.004 - 0.008 mg/L; Results quote EC50 at MIC 0.008 mg/L
    mic_thigh <- fixed(0.008); label("Afabicin desphosphono MIC of the strain in the thigh system (mg/L)") # Supplementary Table S1: ATCC 29213 MIC 0.004 - 0.008 mg/L
    # Starting inocula are set by the experiment, not estimated.
    linoc <- fixed(log(1e6)); label("Log starting in vitro bacterial inoculum (CFU/mL; design input)")        # Methods: target inocula 1.6e5 to 2.0e7 CFU/mL; Figure 2 median starting inoculum about 6 log10 CFU/mL
    linoc_thigh <- fixed(log(1e6)); label("Log thigh bacterial count at infection (CFU/g; design input)")      # Methods: inoculum 8.0e4 to 1.3e7 CFU per thigh, chosen per strain
    fu <- fixed(0.02); label("Unbound fraction of afabicin desphosphono in mouse plasma (fraction)")           # Methods, 'PK/PD model for in vivo data': 'The unbound (fu = 0.02) afabicin desphosphono concentration'

    # ---- Bacterial growth (Table 1, in vitro and in vivo column) ---------
    lkgrowth <- fixed(log(1.10)); label("Log bacterial growth rate constant kgrowth, shared (1/h)")             # Table 1: kgrowth = 1.10 1/h (RSE 2.0%), footnote a: fixed from the growth-control fit
    lkdeath <- fixed(log(0.179)); label("Log bacterial natural death rate constant kdeath, shared (1/h)")       # Table 1: kdeath = 0.179 1/h (Fixed)
    lbmax <- fixed(log(10^9.86)); label("Log maximum system capacity in vitro, Bmax vitro (CFU/mL)")           # Table 1: Bmax vitro = 9.86 log10 CFU/mL (RSE 1.0%), footnote a
    lbmax_thigh <- fixed(log(10^9.08)); label("Log maximum system capacity in the thigh, Bmax vivo (CFU/g)")   # Table 1: Bmax vivo = 9.08 log10 CFU/g (RSE 1.0%), footnote a

    # ---- Drug effect (Table 1, in vitro and in vivo column) -------------
    lemax <- log(3.37); label("Log maximum afabicin desphosphono killing rate constant Emax, shared (1/h)")     # Table 1: Emax = 3.37 1/h (RSE 7.2%)
    # Eq 6: EC50 = beta * MIC^gamma; beta is the EC50 at an MIC of 1 mg/L,
    # carried as the canonical lec50. Only beta differs between settings.
    lec50 <- fixed(log(1)); label("Log EC50 scaling slope beta in vitro, the EC50 at an MIC of 1 mg/L (mg/L)")  # Table 1: beta vitro = 1 (Fixed)
    lec50_thigh <- log(0.707); label("Log EC50 scaling slope beta in vivo, the EC50 at an MIC of 1 mg/L (mg/L)") # Table 1: beta vivo = 0.707 (RSE 6.0%)
    e_mic_ec50 <- 0.204; label("Power of the MIC effect on EC50, gamma, shared (unitless)")                    # Table 1: gamma = 0.204 (RSE 21%)
    lhill <- log(0.324); label("Log sigmoidicity (Hill) factor of the Emax model, shared (unitless)")          # Table 1: Hill = 0.324 (RSE 6.3%)

    # ---- Adaptation (Table 1, in vitro and in vivo column) --------------
    lkon <- log(0.386); label("Log adaptation onset rate constant kon, shared (L/(mg*h))")                      # Table 1: kon = 0.386 mL/(ug*h) (RSE 19%); mL/(ug*h) = L/(mg*h)
    lkoff <- fixed(log(0.0139)); label("Log adaptation offset rate constant koff, shared (1/h)")                # Table 1: koff = 0.0139 1/h (Fixed)
    lar50 <- log(1.18); label("Log fraction in ASon giving a 50% reduction in Emax, AS50, shared (unitless)")   # Table 1: AS50 = 1.18 (RSE 12%)

    # ---- Residual error --------------------------------------------------
    addSd_log_cfu_invitro <- 0.877; label("Additive residual SD on the in vitro log10 count (log10 CFU/mL)")    # Table 1: RES vitro = 0.877 log10 CFU/mL (RSE 5.1%); footnote b: SD scale
    addSd_log_cfu_thigh <- 0.621; label("Additive residual SD on the thigh log10 count (log10 CFU/g)")         # Table 1: RES vivo = 0.621 log10 CFU/g (RSE 4.1%); footnote b: SD scale
  })

  model({
    # ---- 1. Individual parameters ---------------------------------------
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    bmax_thigh <- exp(lbmax_thigh)
    emax <- exp(lemax)
    hill <- exp(lhill)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ar50 <- exp(lar50)

    # Eq 6 with a setting-specific beta (Paramvivo = Paramvitro * Ratio,
    # retained for beta only).
    ec50 <- exp(lec50) * mic^e_mic_ec50
    ec50_thigh <- exp(lec50_thigh) * mic_thigh^e_mic_ec50

    # ---- 2. Drug-effect drivers -----------------------------------------
    # In vitro: static bath concentration in protein-free broth (mg/L).
    conc <- afabicin_desphosphono
    # In vivo: unbound plasma concentration of the active moiety. The
    # prodrug afabicin is assumed to have no killing effect (Methods).
    conc_thigh <- fu * CP_AFABICIN_DESPHOSPHONO_MGL

    # ---- 3. Drug effect, Eq 5 and sigmoid Emax ---------------------------
    emax_t <- emax * (1 - aron / (aron + ar50))
    kdrug <- emax_t * conc^hill / (ec50^hill + conc^hill)
    emax_t_thigh <- emax * (1 - aron_thigh / (aron_thigh + ar50))
    kdrug_thigh <- emax_t_thigh * conc_thigh^hill / (ec50_thigh^hill + conc_thigh^hill)

    # ---- 4. Density-dependent transfer into the dormant state -------------
    btot <- bact_s + bact_d
    ksd <- btot * (kgrowth - kdeath) / bmax
    btot_thigh <- bact_s_thigh + bact_d_thigh
    ksd_thigh <- btot_thigh * (kgrowth - kdeath) / bmax_thigh

    # ---- 5. In vitro system (Eqs 1-4) --------------------------------------
    d/dt(afabicin_desphosphono) <- 0
    d/dt(bact_s) <- kgrowth * bact_s - (kdeath + kdrug) * bact_s - ksd * bact_s
    d/dt(bact_d) <- ksd * bact_s - kdeath * bact_d
    d/dt(aroff) <- koff * aron - kon * conc * aroff
    d/dt(aron) <- kon * conc * aroff - koff * aron

    # ---- 6. Thigh system (Eqs 1-4) -----------------------------------------
    d/dt(bact_s_thigh) <- kgrowth * bact_s_thigh -
      (kdeath + kdrug_thigh) * bact_s_thigh - ksd_thigh * bact_s_thigh
    d/dt(bact_d_thigh) <- ksd_thigh * bact_s_thigh - kdeath * bact_d_thigh
    d/dt(aroff_thigh) <- koff * aron_thigh - kon * conc_thigh * aroff_thigh
    d/dt(aron_thigh) <- kon * conc_thigh * aroff_thigh - koff * aron_thigh

    # ---- 7. Initial conditions --------------------------------------------
    # Methods: all bacteria start in S; ASoff = 1 and ASon = 0 at start.
    # Thigh time zero is the time of infection; treatment began 2 h later.
    bact_s(0) <- exp(linoc)
    aroff(0) <- 1
    bact_s_thigh(0) <- exp(linoc_thigh)
    aroff_thigh(0) <- 1

    # ---- 8. Observations ---------------------------------------------------
    # Log10-transform-both-sides with setting-specific additive errors.
    log_cfu_invitro <- log10(btot)
    log_cfu_thigh <- log10(btot_thigh)
    log_cfu_invitro ~ add(addSd_log_cfu_invitro)
    log_cfu_thigh ~ add(addSd_log_cfu_thigh)
  })
}
