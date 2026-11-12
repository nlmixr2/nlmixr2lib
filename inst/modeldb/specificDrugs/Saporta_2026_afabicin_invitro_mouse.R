Saporta_2026_afabicin_invitro_mouse <- function() {
  description <- "In vitro (Staphylococcus aureus) and preclinical (mouse, immunocompetent thigh infection). Semi-mechanistic pharmacodynamic model of afabicin desphosphono, the active moiety of the FabI-inhibitor prodrug afabicin, fitted jointly to static in vitro time-kill data and to immunocompetent mouse thigh infection data, with a neutrophil-mediated immune response in the thigh. Two parallel bacterial systems share kgrowth, kdeath, gamma, Hill, kon, koff and AS50: an in vitro system driven by a static bath concentration, and a thigh system driven by the unbound plasma concentration fu * CP_AFABICIN_DESPHOSPHONO_MGL with fu = 0.02. Each system has growing drug-susceptible bacteria (S) and dormant drug-insusceptible bacteria (D) entered at kSD = (S + D) * (kgrowth - kdeath) / Bmax, a sigmoid Emax kill on S with EC50 = beta * MIC^gamma, and a reversible adaptation pair (aroff -> aron at kon * C, back at koff) that reduces Emax as Emax * (1 - aron / (aron + AS50)). In the thigh, Emax and beta are lower than in vitro (Emax vivo = Emax vitro * RatioEmax; separate beta vivo), and neutrophils phagocytose S and D bacteria at kphag = kphag,max * (1 - (P/N)^H / ((P/N)^H + (P/N50)^H)) into a phagocytosed state (P) that is digested at kdig = kphag,max; non-digested P bacteria are counted as CFU. Thigh neutrophils come from a fixed Friberg maturation chain (prol, transit1-3, circ; MTT 37.7 h, feedback 0.149) with a post-infection recruitment surge (amplitude 44.2, peak 50.5 h, width 28.2 h). The ATCC 33591 starting inoculum carries an estimated bioavailability of 0.437. The mouse PK model (Bader 2017, a conference poster) reports no parameter values, so the plasma concentration is a time-varying input. No inter-experiment variability; additive residual error on log10 CFU. Sibling models: Saporta_2024_afabicin_invitro_mouse (the joint fit with neutropenic mice) and Saporta_2024_afabicin_invitro."
  reference <- paste(
    "Saporta R, Nielsen EI, Menetrey A, Cameron DR, Nicolas-Metral V,",
    "Friberg LE. Model-based quantification of immune response and",
    "anti-staphylococcal activity of afabicin in immunocompetent mouse thigh",
    "infections to enable predictions of clinical efficacy. Antimicrob Agents",
    "Chemother. 2026;70(4):e00959-25. doi:10.1128/aac.00959-25. PMID:",
    "41778916. PMCID: PMC13041350. Model structure: Materials and Methods,",
    "Eqs 1-6, and Figure 1; neutrophil dynamics: Supplementary text 1.",
    "Parameter estimates: Table 1. Strain MICs: Supplementary Table S1.",
    "Unbound fraction fu = 0.02 and the shared effect structure: Saporta R et",
    "al. J Antimicrob Chemother. 2024;79(12):3150-3159.",
    "doi:10.1093/jac/dkae334. The mouse PK model that supplied the plasma",
    "concentration (structure in Supplementary text 2) is Bader JC, Lakota",
    "EA, Bravo J et al., ASM Microbe 2017, whose parameter values are",
    "unpublished.",
    sep = " "
  )
  vignette <- "Saporta_2026_afabicin"
  units <- list(time = "h", dosing = "mg/L (in vitro bath concentration)", concentration = "mg/L")

  # The only dosed state is the static in vitro bath concentration; the thigh
  # system is driven by the CP_AFABICIN_DESPHOSPHONO_MGL covariate.
  dosing <- c("afabicin_desphosphono")

  # State names follow Saporta_2024_afabicin_invitro_mouse.R; the thigh
  # phagocytosed state and the thigh neutrophil pool are added downstream of
  # the canonical Friberg chain (prol, transit1-3, circ), as neut_lung is in
  # Thorsted_2023_neutrophil_mouse.R.
  paper_specific_compartments <- c(
    "afabicin_desphosphono",
    "bact_s",
    "bact_d",
    "bact_s_thigh",
    "bact_d_thigh",
    "bact_phagocytosed_thigh",
    "aroff_thigh",
    "aron_thigh",
    "neut_thigh",
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
    prol = list(
      analyte = "proliferating neutrophil progenitor cells (Nprol)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "maturing neutrophils, transit 1 (Ntr1)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing neutrophils, transit 2 (Ntr2)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing neutrophils, transit 3 (Ntr3)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(
      analyte = "circulating neutrophils (Ncirc)",
      units = "cells/mL",
      specimen = "whole blood",
      verified = TRUE
    ),
    neut_thigh = list(
      analyte = "neutrophils at the thigh infection site (Nthigh)",
      units = "cells/mL",
      specimen = "tissue",
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
    bact_phagocytosed_thigh = list(
      analyte = "Staphylococcus aureus in the mouse thigh phagocytosed by neutrophils and not yet digested (P)",
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
        "Materials and Methods, 'PKPD modeling of afabicin desphosphono activity': 'Unbound afabicin desphosphono (active metabolite) concentrations were derived from a PK model for afabicin in mice (Text S2) (31)'.",
        "Supplementary text 2 prints the structure of that model (IP depot -> absorption transit -> two-compartment afabicin desphosphono with an afabicin-to-desphosphono conversion factor) but no source reports its parameter values (Bader 2017, ASM Microbe abstract), so the concentration enters as a covariate.",
        "The model multiplies this total concentration by fu = 0.02 (Saporta 2024, the effect model this paper re-uses) to obtain the unbound concentration that drives killing and adaptation. Set to 0 for growth-control animals, before the first dose, and for records that simulate only the in vitro system."
      )
    )
  )

  population <- list(
    species = "in vitro (Staphylococcus aureus) + mouse (female ICR or Swiss Webster, immunocompetent)",
    n_subjects = 819L,
    n_studies = NA_integer_,
    organism = "Staphylococcus aureus: 5 strains in vivo (570-493, ATCC 29213, ATCC 33591, IHMA 1073118 controls only, IHMA 1074670; MIC 0.004 to 0.06 mg/L; Supplementary Table S1), analysed jointly with the in vitro time-kill data of Saporta 2024",
    system = "Immunocompetent mouse thigh infection model: intramuscular inoculum of 6.76e5 to 1.26e8 CFU per thigh (controls) or 6.92e6 to 1.15e7 CFU per thigh (treated); thighs harvested at 2, 4, 6, 10, 12, 20, 26, 50 or 74 h after infection; counts in CFU per gram of thigh. In vitro arm: static time-kill experiments in broth (Saporta 2024).",
    dose_range = "Afabicin (prodrug) 0.011 to 150 mg/kg intraperitoneally every 6 h, started 2 h after infection",
    sex_female_pct = 100,
    disease_state = "Thigh infection with Staphylococcus aureus in immunocompetent mice (in vivo arm)",
    notes = paste(
      "n_subjects counts the mice (Materials and Methods: n = 819).",
      "Bacterial growth parameters (kgrowth, kdeath, Bmax) and the neutrophil-dynamics parameters (Thorsted 2023) were fixed; immune-response and drug-effect parameters were estimated jointly on the control and afabicin-treated immunocompetent data and the in vitro data. The neutropenic mouse data of Saporta 2024 were not used.",
      "Estimation in NONMEM 7.5.0 (FOCE) on log10-transformed counts with setting-specific additive errors; inter-individual variability was not evaluated."
    )
  )

  ini({
    # ---- Simulation design inputs ---------------------------------------
    # Strain MICs scale EC50 (Eq 6); fixed ini() entries, as in
    # Saporta_2024_afabicin_invitro_mouse.R, so a user re-points either
    # system at another isolate with rxSolve(params = ...). Default: S. aureus
    # ATCC 29213, the strain shared with the in vitro and neutropenic studies.
    mic <- fixed(0.008); label("Afabicin desphosphono MIC of the strain in the in vitro system (mg/L)")   # Saporta 2024 Supplementary Table S1: ATCC 29213 MIC 0.004 - 0.008 mg/L
    mic_thigh <- fixed(0.008); label("Afabicin desphosphono MIC of the strain in the thigh system (mg/L)") # Supplementary Table S1: ATCC 29213 MIC 0.008 mg/L
    # Starting inocula are set by the experiment, not estimated.
    linoc <- fixed(log(1e6)); label("Log starting in vitro bacterial inoculum (CFU/mL; design input)")    # Saporta 2024 Methods: in vitro inocula about 6 log10 CFU/mL
    linoc_thigh <- fixed(log(10^7.2)); label("Log thigh bacterial count at infection (CFU/g; design input)") # Table 2: treatment-group inoculum sizes 7.14 to 7.36 log10 CFU/g
    strain_atcc33591 <- fixed(0); label("Indicator that the thigh strain is S. aureus ATCC 33591 (1) or another strain (0); design input") # Results: F33591 applies to ATCC 33591 only ('did not improve the model fit when evaluated for other strains')
    fu <- fixed(0.02); label("Unbound fraction of afabicin desphosphono in mouse plasma (fraction)")      # Saporta 2024 (reference 28) Methods: 'The unbound (fu = 0.02) afabicin desphosphono concentration'

    # ---- Bacterial dynamics (Table 1) --------------------------------------
    lkgrowth <- fixed(log(1.10)); label("Log bacterial growth rate constant kgrowth, shared (1/h)")             # Table 1: kgrowth = 1.10 1/h (Fixed)
    lkdeath <- fixed(log(0.179)); label("Log bacterial natural death rate constant kdeath, shared (1/h)")       # Table 1: kdeath = 0.179 1/h (Fixed)
    lbmax <- fixed(log(10^9.86)); label("Log maximum system capacity in vitro, Bmax vitro (CFU/mL)")           # Table 1: Bmax vitro = 9.86 log10 CFU/mL (Fixed)
    lbmax_thigh <- fixed(log(10^9.08)); label("Log maximum system capacity in the thigh, Bmax vivo (CFU/g)")   # Table 1: Bmax vivo = 9.08 log10 CFU/g (Fixed)

    # ---- Immune response (Table 1) -----------------------------------------
    finoc_atcc33591 <- 0.437; label("Inoculum bioavailability for S. aureus ATCC 33591 in immunocompetent mice (fraction)") # Table 1: F33591 = 0.437 (RSE 6.2%)
    lcirc0 <- fixed(log(10^6.23)); label("Log baseline circulating neutrophil count Ncirc,T0 (cells/mL)")       # Table 1: Ncirc,T0 = 6.23 log10 N/mL (Fixed)
    lmtt <- fixed(log(37.7)); label("Log neutrophil maturation mean transit time MTT (h)")                      # Table 1: MTT = 37.7 h (Fixed)
    lgamma <- fixed(log(0.149)); label("Log feedback exponent on (Ncirc,T0 / Ncirc) (unitless)")               # Supplementary text 1: dNprol/dt feedback exponent 0.149
    lnamp <- fixed(log(44.2)); label("Log amplitude of the post-infection neutrophil recruitment surge Namp (unitless)") # Table 1: Namp = 44.2 (Fixed)
    lnt0 <- fixed(log(50.5)); label("Log time after infection of the recruitment surge peak NT0 (h)")          # Table 1: NT0 = 50.5 h (Fixed)
    lnsw <- fixed(log(28.2)); label("Log width of the recruitment surge NSW (h)")                              # Table 1: NSW = 28.2 h (Fixed)
    # kphag,max and kdig share one estimate (Table 1 footnote b; Results
    # 'kdig could share the same estimate as kphag,max'), carried as the
    # Thorsted_2023_neutrophil_mouse.R canonical lkn_max.
    lkn_max <- log(0.963); label("Log maximum phagocytosis rate constant kphag,max, equal to the digestion rate constant kdig (1/h)") # Table 1: kphag,max = kdig = 0.963 1/h (RSE 2.0%)
    lpn50 <- log(114); label("Log phagocytosed-bacteria-to-neutrophil ratio halving the phagocytosis rate P/N50 (CFU/cell)") # Table 1: P/N50 = 114 CFU/N (RSE 16%)
    lhill_pn <- log(0.475); label("Log sigmoidicity factor H of the phagocytosis saturation (unitless)")     # Table 1: H = 0.475 (RSE 12%)

    # ---- Afabicin desphosphono effect (Table 1) -----------------------------
    lemax <- log(2.99); label("Log maximum afabicin desphosphono killing rate constant in vitro, Emax vitro (1/h)") # Table 1: Emax vitro = 2.99 1/h (RSE 5.2%)
    lratio_emax_thigh <- log(0.0796); label("Log ratio of the in vivo to the in vitro Emax, RatioEmax (unitless)") # Table 1: RatioEmax = 0.0796 (RSE 9.3%); Emax vivo = Emax vitro * RatioEmax
    # Eq 6: EC50 = beta * MIC^gamma; beta is the EC50 at an MIC of 1 mg/L,
    # carried as the canonical lec50.
    lec50 <- fixed(log(1)); label("Log EC50 scaling slope beta in vitro, the EC50 at an MIC of 1 mg/L (mg/L)")   # Table 1: beta vitro = 1 (Fixed)
    lec50_thigh <- log(0.182); label("Log EC50 scaling slope beta in vivo, the EC50 at an MIC of 1 mg/L (mg/L)")  # Table 1: beta vivo = 0.182 (RSE 16%)
    e_mic_ec50 <- 0.185; label("Power of the MIC effect on EC50, gamma, shared (unitless)")                     # Table 1: gamma = 0.185 (RSE 24%)
    lhill <- log(0.515); label("Log sigmoidicity (Hill) factor of the Emax model, shared (unitless)")           # Table 1: Hill = 0.515 (RSE 7.4%)

    # ---- Adaptation (Table 1) ----------------------------------------------
    lkon <- log(0.556); label("Log adaptation onset rate constant kon, shared (L/(mg*h))")                       # Table 1: kon = 0.556 mL/(ug*h) (RSE 25%); mL/(ug*h) = L/(mg*h)
    lkoff <- fixed(log(0.0139)); label("Log adaptation offset rate constant koff, shared (1/h)")                 # Table 1: koff = 0.0139 1/h (Fixed)
    lar50 <- log(0.784); label("Log fraction in ASon giving a 50% reduction in Emax, AS50, shared (unitless)")   # Table 1: AS50 = 0.784 (RSE 10%)

    # ---- Residual error ----------------------------------------------------
    addSd_log_cfu_invitro <- 0.877; label("Additive residual SD on the in vitro log10 count (log10 CFU/mL)")     # Table 1: RES vitro = 0.877 log10 CFU/mL (RSE 4.9%); footnote a: SD scale
    addSd_log_cfu_thigh <- 0.555; label("Additive residual SD on the thigh log10 count (log10 CFU/g)")          # Table 1: RES vivo = 0.555 log10 CFU/g (RSE 4.2%); footnote a: SD scale
  })

  model({
    # ---- 1. Individual parameters ---------------------------------------
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    bmax_thigh <- exp(lbmax_thigh)
    hill <- exp(lhill)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ar50 <- exp(lar50)

    circ0 <- exp(lcirc0)
    ktr <- 4 / exp(lmtt)
    gamma <- exp(lgamma)
    namp <- exp(lnamp)
    nt0 <- exp(lnt0)
    nsw <- exp(lnsw)
    kn_max <- exp(lkn_max)
    pn50 <- exp(lpn50)
    hill_pn <- exp(lhill_pn)

    # Emax vivo = Emax vitro * RatioEmax (Table 1).
    emax <- exp(lemax)
    emax_thigh <- emax * exp(lratio_emax_thigh)
    # Eq 6 with a setting-specific beta.
    ec50 <- exp(lec50) * mic^e_mic_ec50
    ec50_thigh <- exp(lec50_thigh) * mic_thigh^e_mic_ec50

    # ---- 2. Drug-effect drivers -----------------------------------------
    # In vitro: static bath concentration in protein-free broth (mg/L).
    conc <- afabicin_desphosphono
    # In vivo: unbound plasma concentration of the active moiety.
    conc_thigh <- fu * CP_AFABICIN_DESPHOSPHONO_MGL

    # ---- 3. Drug effect: sigmoid Emax reduced by adaptation ---------------
    emax_t <- emax * (1 - aron / (aron + ar50))
    kdrug <- emax_t * conc^hill / (ec50^hill + conc^hill)
    emax_t_thigh <- emax_thigh * (1 - aron_thigh / (aron_thigh + ar50))
    kdrug_thigh <- emax_t_thigh * conc_thigh^hill / (ec50_thigh^hill + conc_thigh^hill)

    # ---- 4. Density-dependent transfer into the dormant state -------------
    # Bmax limits the extracellular S + D bacteria; phagocytosed bacteria do
    # not count (Figure 1 draws kSD driven by 'S+D').
    btot <- bact_s + bact_d
    ksd <- btot * (kgrowth - kdeath) / bmax
    bext_thigh <- bact_s_thigh + bact_d_thigh
    ksd_thigh <- bext_thigh * (kgrowth - kdeath) / bmax_thigh

    # ---- 5. Neutrophil dynamics (Supplementary text 1, Eq 1) ---------------
    # Model time is the time since infection, so the surge starts at t = 0
    # (Figure S1: 'The surge begins at time of infection').
    d/dt(prol) <- -ktr * prol + ktr * prol * (circ0 / circ)^gamma
    d/dt(transit1) <- -ktr * transit1 + ktr * prol
    d/dt(transit2) <- -ktr * transit2 + ktr * transit1
    d/dt(transit3) <- -ktr * transit3 + ktr * transit2
    d/dt(circ) <- -ktr * circ + ktr * transit3
    surge <- namp / (((t - nt0) / nsw)^6 + 1)
    d/dt(neut_thigh) <- -ktr * neut_thigh + ktr * (1 + surge) * circ

    # ---- 6. Saturable phagocytosis (Eq 5); digestion is not saturable ------
    pn_ratio <- bact_phagocytosed_thigh / neut_thigh
    kphag <- kn_max * (1 - pn_ratio^hill_pn / (pn_ratio^hill_pn + pn50^hill_pn))
    kdig <- kn_max

    # ---- 7. In vitro system -----------------------------------------------
    d/dt(afabicin_desphosphono) <- 0
    d/dt(bact_s) <- kgrowth * bact_s - (kdeath + kdrug) * bact_s - ksd * bact_s
    d/dt(bact_d) <- ksd * bact_s - kdeath * bact_d
    d/dt(aroff) <- koff * aron - kon * conc * aroff
    d/dt(aron) <- kon * conc * aroff - koff * aron

    # ---- 8. Thigh system (Eqs 2-4 with the drug effect on S) ---------------
    d/dt(bact_s_thigh) <- kgrowth * bact_s_thigh -
      (kdeath + kdrug_thigh) * bact_s_thigh - ksd_thigh * bact_s_thigh -
      kphag * bact_s_thigh
    d/dt(bact_d_thigh) <- ksd_thigh * bact_s_thigh - kdeath * bact_d_thigh -
      kphag * bact_d_thigh
    d/dt(bact_phagocytosed_thigh) <- kphag * bact_s_thigh + kphag * bact_d_thigh -
      kdig * bact_phagocytosed_thigh
    d/dt(aroff_thigh) <- koff * aron_thigh - kon * conc_thigh * aroff_thigh
    d/dt(aron_thigh) <- kon * conc_thigh * aroff_thigh - koff * aron_thigh

    # ---- 9. Initial conditions --------------------------------------------
    # All bacteria start in S; ASoff = 1 and ASon = 0. Every neutrophil state
    # starts at Ncirc,T0 (Supplementary text 1). The ATCC 33591 inoculum
    # enters with bioavailability F33591.
    bact_s(0) <- exp(linoc)
    aroff(0) <- 1
    finoc <- 1 + strain_atcc33591 * (finoc_atcc33591 - 1)
    bact_s_thigh(0) <- finoc * exp(linoc_thigh)
    aroff_thigh(0) <- 1
    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0
    neut_thigh(0) <- circ0

    # ---- 10. Observations --------------------------------------------------
    # Phagocytic cells were lysed before plating, so non-digested P bacteria
    # are counted. Log10-transform-both-sides with setting-specific errors.
    log_cfu_invitro <- log10(btot)
    log_cfu_thigh <- log10(bext_thigh + bact_phagocytosed_thigh)
    log_cfu_invitro ~ add(addSd_log_cfu_invitro)
    log_cfu_thigh ~ add(addSd_log_cfu_thigh)
  })
}
