HernandezLozano_2025_apramycin_mouse <- function() {
  description <- "Preclinical (mouse, C3H/HeJ female). Translational pharmacokinetic-pharmacodynamic model of apramycin in a murine model of complicated urinary tract infection, describing the time course of Escherichia coli burden in kidneys and bladder tissue after twice-daily subcutaneous dosing. A one-compartment subcutaneous PK model with a dose-dependent absorption rate constant (Sou 2021) drives an unbound plasma concentration into two independent organ pharmacodynamic systems, one for kidneys and one for bladder. Each organ carries the in vitro time-kill structure re-estimated on the in vivo data: two bacterial subpopulations (1 = main, apramycin-susceptible; 2 = decreased susceptibility), each with a growing drug-susceptible state S and a dormant drug-insusceptible state D entered at rate ksr = (kg - kd) * Btot / Bmax. Apramycin adds to the death rate of the S states through a power model normalized to the in vitro MIC, kdrug = Slope * (Cu/MIC)^gamma with gamma fixed to 1, and drives transfer from subpopulation 1 to subpopulation 2 at rate kada * (Cu/MIC). The kidney is assumed to be at pH 7.4 and the bladder at pH 6, so each organ normalizes by the MIC measured at its own pH. Net bacterial growth is 76-99% lower and apramycin potency 3- to 145-fold higher in vivo than in vitro. Model time zero is 6 h after bacterial inoculation. Sibling models: HernandezLozano_2025_apramycin_invitro (the underlying time-kill fit) and HernandezLozano_2025_apramycin_human (the same in vivo PD component driven by the human population PK model of Zhao 2022)."
  reference <- "Hernandez-Lozano I, Aranzana-Climent V, Cao S, Matias C, Hansen JU, Liepinsh E, Hughes D, Hobbie SN, Vingsbo Lundberg C, Friberg LE. Model-informed drug development for antimicrobials: translational pharmacokinetic-pharmacodynamic modelling of apramycin to facilitate prediction of efficacious dose in complicated urinary tract infections. J Antimicrob Chemother. 2025 Feb 3;80(2):302-311. doi:10.1093/jac/dkae409. PMID: 39545353. PMCID: PMC11695905. PD structure: Materials and methods, 'PKPD modelling', plus the schematic in Figure 1. In vivo PD estimates and 95% CIs: Table 1, section 'In vivo PD parameters'. MIC values and the kidney-pH-7.4 / bladder-pH-6 assumption: Results, 'In vitro time-kill curves and PD modelling', and Discussion. Unbound fraction in mouse plasma (91.6%): Materials and methods, 'PKPD modelling'. Mouse PK parameters are NOT re-estimated here; they are taken from Sou T, Hansen J, Liepinsh E, et al. Model-informed drug development for antimicrobials: translational PK and PK/PD modeling to predict an efficacious human dose for apramycin. Clin Pharmacol Ther 2021;109:1063-73. doi:10.1002/cpt.2104, Table 1 (mouse column) and Eqs. 1-3."
  vignette <- "HernandezLozano_2025_apramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Two independent organ PD systems, each replicating the in vitro
  # structure. Suffix k = kidneys (cfu per two kidneys), b = bladder
  # (cfu per bladder). Within an organ, subpopulation 1 is the main
  # apramycin-susceptible population and 2 the less-susceptible one; S is the
  # growing drug-susceptible state and D the dormant drug-insusceptible one
  # (Figure 1). The names do not match the canonical
  # bact_<phenotype><lifecycle-digit> scheme because that trailing digit
  # indexes the Bulitta / Wicha replication life cycle, a different
  # distinction from this model's growing-versus-dormant states, and because
  # that scheme has no organ dimension.
  paper_specific_compartments <- c(
    "bact_s1k", "bact_d1k", "bact_s2k", "bact_d2k",
    "bact_s1b", "bact_d1b", "bact_s2b", "bact_d2b"
  )

  compartmentData <- list(
    depot = list(analyte = "apramycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apramycin", units = "mg", specimen = "plasma", verified = TRUE),
    bact_s1k = list(analyte = "Escherichia coli in the kidneys, main apramycin-susceptible subpopulation, growing drug-susceptible (S) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_d1k = list(analyte = "Escherichia coli in the kidneys, main apramycin-susceptible subpopulation, dormant drug-insusceptible (D) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_s2k = list(analyte = "Escherichia coli in the kidneys, subpopulation with decreased apramycin susceptibility, growing drug-susceptible (S) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_d2k = list(analyte = "Escherichia coli in the kidneys, subpopulation with decreased apramycin susceptibility, dormant drug-insusceptible (D) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_s1b = list(analyte = "Escherichia coli in the bladder, main apramycin-susceptible subpopulation, growing drug-susceptible (S) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_d1b = list(analyte = "Escherichia coli in the bladder, main apramycin-susceptible subpopulation, dormant drug-insusceptible (D) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_s2b = list(analyte = "Escherichia coli in the bladder, subpopulation with decreased apramycin susceptibility, growing drug-susceptible (S) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE),
    bact_d2b = list(analyte = "Escherichia coli in the bladder, subpopulation with decreased apramycin susceptibility, dormant drug-insusceptible (D) state", units = "CFU/organ", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial-strain indicator selecting the strain-specific in vivo growth and death rate constants in each organ, and the strain's in vitro MICs at pH 7.4 (used in the kidney) and pH 6 (used in the bladder).",
      units = "(categorical)",
      type = "categorical",
      source_name = "STRAIN",
      reference_category = NULL,
      notes = paste(
        "591 = Escherichia coli EN591, an MDR rmtB clinical isolate; MIC 8 mg/L at pH 7.4 and 32 mg/L at pH 6.",
        "700336 = Escherichia coli ATCC 700336 (EN1085), a trimethoprim/sulfamethoxazole-resistant urinary tract infection isolate; MIC 4 mg/L at pH 7.4 and 16 mg/L at pH 6.",
        "There is no reference level: each level selects its own kg, Redukg and kd per organ rather than contributing a contrast against a baseline. The drug-effect parameters SlopeS, SlopeR and kada are SHARED across strains in both organs (Table 1, Strain column 'Both'); strain differences in apramycin effect are carried entirely by the strain-specific MIC used to normalize the concentration.",
        "The sibling file HernandezLozano_2025_apramycin_invitro.R uses the same two strain codes, but additionally carries PH_MEDIUM, because there the growth-medium pH is an experimental factor that indexes its own drug-effect parameter set. This file omits PH_MEDIUM because in vivo pH is a property of the organ (kidney assumed 7.4, bladder assumed 6) rather than something the experiment varied, so it is baked into which per-site MIC applies."
      )
    ),
    WT = list(
      description = "Mouse body weight, used to allometrically scale apramycin clearance and volume of distribution from the 70 kg reference of the upstream preclinical PK model.",
      units = "kg",
      type = "continuous",
      source_name = "WT",
      reference_category = NULL,
      notes = "Sou 2021 Eqs. 1-2 parameterize every preclinical species as CL_i = CL * (WT/70)^0.75 and V_i = V * (WT/70)^1, where CL and V are the typical parameters of a 70 kg human adult. The C3H/HeJ mice of this study had a mean weight of 16.8 g (range 12.3-21.2 g; Supplementary methods, 'Mouse PD experiments'), so the default here is 0.0168 kg. The paper's separate PBPK cross-check used a 20 g average-weight mouse instead; that PBPK model is not part of this file."
    ),
    DOSE_APRAMYCIN_MGKG = list(
      description = "Administered apramycin dose level in mg/kg, driving the dose-dependent subcutaneous absorption rate constant.",
      units = "mg/kg",
      type = "continuous",
      source_name = "DOSE",
      reference_category = NULL,
      notes = "Sou 2021 Eq. 3: ka = ka30 * (Dose/30)^pow, where Dose is the mg/kg dose level and ka30 is the absorption rate constant at 30 mg/kg. The dose level, not the absolute amount, enters the relationship, so this must be carried as a data column separate from the event-table amt (which is in mg). Dose levels used in this study: 1.5, 3, 5, 10, 15 and 30 mg/kg twice daily (Supplementary methods), plus 0.03, 0.1, 0.3 and 1 mg/kg in the EN591 low-dose validation study and 0.05, 0.2, 0.8, 3.2, 12.8 and 51.2 mg/kg in the ATCC 700336 validation study."
    )
  )

  population <- list(
    species = "mouse (C3H/HeJ, female, 6 weeks old)",
    n_subjects = 286L,
    n_studies = 3L,
    weight_range = "12.3-21.2 g",
    weight_median = "16.8 g (mean)",
    sex_female_pct = 100,
    disease_state = "Complicated urinary tract infection (ascending cUTI) established by transurethral inoculation of 5x10^7 CFU of Escherichia coli EN591 or ATCC 700336 into the bladder; mice were immunocompetent and received 5% glucose in the drinking water to induce diuresis",
    dose_range = "Apramycin 1.5-30 mg/kg subcutaneously twice daily in the primary studies (3, 10 and 30 mg/kg for EN591; 1.5, 5 and 15 mg/kg for ATCC 700336, halved because of the lower MIC), given at 24, 30, 48, 54, 72 and 78 h after inoculation. Validation studies used 0.03-1 mg/kg (EN591) and 0.05-51.2 mg/kg (ATCC 700336)",
    regions = "Statens Serum Institut, Denmark (Studies 1 and 2); Pharmacology Discovery Services Taiwan (Study 3)",
    notes = paste(
      "Animal counts by strain, dose and sampling time are in Supplementary Table S1: 286 mice in total across three studies.",
      "Sampling times were 6, 10, 24, 30, 48, 72 and 96 h after bacterial inoculation (168 h in Study 3); the primary outcome was CFU per two kidneys or per bladder.",
      "Study 3 (ATCC 700336, Pharmacology Discovery Services Taiwan) infected animals 96 h rather than 24 h before the start of treatment and used a 10^9 CFU target inoculum.",
      "Plasma PK was measured in EN591-infected mice at 15 min, 30 min and 1 h after the first and last dose and was described by the previously published one-compartment model with dose-dependent absorption; the PK parameters were not re-estimated in this paper.",
      "MODEL TIME ZERO IS 6 h AFTER BACTERIAL INOCULATION, the first sampling time, at which the initial bacterial densities Inoc_k and Inoc_b apply. Treatment in the primary studies therefore starts at model time 18 h and doses fall at 18, 24, 42, 48, 66 and 72 h."
    )
  )

  ini({
    # ---- Mouse population PK (Sou 2021 Table 1, mouse column) ------------
    # Not re-estimated in Hernandez-Lozano 2025; carried in from the cited
    # upstream model, so every PK parameter is fixed here.
    lcl <- fixed(log(8.49)); label("Log apramycin clearance for a 70 kg reference body weight (L/h)")            # Sou 2021 Table 1, mouse: CL/F = 8.49 L/hour/70 kg
    lvc <- fixed(log(6.55)); label("Log apramycin volume of distribution for a 70 kg reference body weight (L)") # Sou 2021 Table 1, mouse: V/F = 6.55 L/70 kg
    lka30 <- fixed(log(2.17)); label("Log subcutaneous absorption rate constant at the 30 mg/kg dose level (1/h)") # Sou 2021 Table 1, mouse: ka30 = 2.17 /hour
    powka <- fixed(-0.160); label("Power of the dose level on the absorption rate constant (unitless)")          # Sou 2021 Table 1, mouse: pow = -0.160; Eq. 3: ka = ka30 * (Dose/30)^pow
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on clearance (unitless)")                  # Sou 2021 Eq. 1: CL_i = CL * (WT/70)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on volume of distribution (unitless)")        # Sou 2021 Eq. 2: V_i = V * (WT/70)^1
    lfdepot <- fixed(log(1)); label("Log subcutaneous bioavailability (unitless)")                               # Sou 2021: "Because only s.c. data were available for mice, rats, and guinea pigs, bioavailability (F) was set to 1"
    fu <- fixed(0.916); label("Fraction of apramycin unbound in mouse plasma (unitless)")                        # Methods, PKPD modelling: "using an unbound fraction of apramycin in plasma of 91.6%"

    # ---- Strain MICs used to normalize the drug effect --------------------
    # Kidney is assumed to be at pH 7.4 and bladder at pH 6 (Discussion), so
    # each organ normalizes by the MIC measured at its own pH.
    micEn591ph74 <- fixed(8); label("MIC of apramycin for EN591 at pH 7.4 (mg/L)")            # Results: MIC_EN591 = 8 mg/L at pH 7.4
    micEn591ph6 <- fixed(32); label("MIC of apramycin for EN591 at pH 6 (mg/L)")              # Results: MIC_EN591 = 32 mg/L at pH 6
    micAtccph74 <- fixed(4); label("MIC of apramycin for ATCC 700336 at pH 7.4 (mg/L)")       # Results: MIC_ATCC700336 = 4 mg/L at pH 7.4
    micAtccph6 <- fixed(16); label("MIC of apramycin for ATCC 700336 at pH 6 (mg/L)")         # Results: MIC_ATCC700336 = 16 mg/L at pH 6

    # ---- In vivo PD, kidney (Table 1, In vivo PD parameters) -------------
    lkg1kEn591 <- log(0.694); label("Log kidney growth rate constant of the susceptible subpopulation, EN591 (1/h)")      # Table 1: kg1k = 0.694 (95% CI 0.432-1.064), EN591
    lkg1kAtcc <- log(0.205); label("Log kidney growth rate constant of the susceptible subpopulation, ATCC 700336 (1/h)") # Table 1: kg1k = 0.205 (95% CI 0.187-0.228), ATCC700336
    redukgkEn591 <- 25.4; label("Kidney fractional reduction in growth rate for the resistant subpopulation, EN591 (%)")      # Table 1: Reduckgk = 25.4 % (95% CI 16.8-34.5), EN591
    redukgkAtcc <- 6.9; label("Kidney fractional reduction in growth rate for the resistant subpopulation, ATCC 700336 (%)")  # Table 1: Reduckgk = 6.9 % (95% CI 0.5-16.9), ATCC700336
    lkdkEn591 <- log(0.526); label("Log kidney bacterial death rate constant, EN591 (1/h)")                               # Table 1: kd,k = 0.526 (95% CI 0.348-0.745), EN591
    lkdkAtcc <- fixed(log(0.179)); label("Log kidney bacterial death rate constant, ATCC 700336 (1/h)")                   # Table 1: kd,k = 0.179, Fixed for ATCC700336 because "the scarcity and variability of the observed data did not allow estimation of physiologically reasonable values"
    slopeSk <- 9.05; label("Apramycin effect on the susceptible subpopulation in kidney (1/h)")                           # Table 1: SlopeSk = 9.05 (95% CI 5.41-13.11), both strains
    slopeRk <- 0.066; label("Apramycin effect on the resistant subpopulation in kidney (1/h)")                            # Table 1: SlopeRk = 0.066 (95% CI 0.026-0.109), both strains
    kadak <- 0.031; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation in kidney (1/h)") # Table 1: kada,k = 0.031 (95% CI 0.006-0.062), both strains
    lbmaxk <- fixed(log(10^6.49)); label("Log maximum bacterial density in kidneys (CFU/organ)")                          # Table 1: Bmax,k = 6.49 log10 CFU/organ, Fixed
    linock <- log(10^5.70); label("Log initial bacterial density in kidneys at 6 h post-inoculation (CFU/organ)")         # Table 1: Inock = 5.70 log10 CFU/organ (95% CI 5.38-6.04)

    # ---- In vivo PD, bladder (Table 1, In vivo PD parameters) ------------
    lkg1bEn591 <- log(1.51); label("Log bladder growth rate constant of the susceptible subpopulation, EN591 (1/h)")       # Table 1: kg1b = 1.51 (95% CI 1.09-1.84), EN591
    lkg1bAtcc <- log(0.228); label("Log bladder growth rate constant of the susceptible subpopulation, ATCC 700336 (1/h)") # Table 1: kg1b = 0.228 (95% CI 0.206-0.254), ATCC700336
    redukgbEn591 <- 25.2; label("Bladder fractional reduction in growth rate for the resistant subpopulation, EN591 (%)")      # Table 1: Reduckgb = 25.2 % (95% CI 19.2-30.3), EN591
    redukgbAtcc <- 6.9; label("Bladder fractional reduction in growth rate for the resistant subpopulation, ATCC 700336 (%)")  # Table 1: Reduckgb = 6.9 % (95% CI 0.5-16.9), ATCC700336
    lkdbEn591 <- log(1.16); label("Log bladder bacterial death rate constant, EN591 (1/h)")                                # Table 1: kd,b = 1.16 (95% CI 0.85-1.39), EN591
    lkdbAtcc <- fixed(log(0.179)); label("Log bladder bacterial death rate constant, ATCC 700336 (1/h)")                   # Table 1: kd,b = 0.179, Fixed for ATCC700336
    slopeSb <- 191; label("Apramycin effect on the susceptible subpopulation in bladder (1/h)")                            # Table 1: SlopeSb = 191 (95% CI 58-1147), both strains; %RSE 159%
    slopeRb <- 0.276; label("Apramycin effect on the resistant subpopulation in bladder (1/h)")                            # Table 1: SlopeRb = 0.276 (95% CI 0.153-0.403), both strains
    kadab <- 1.35; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation in bladder (1/h)") # Table 1: kada,b = 1.35 (95% CI 0.43-8.86), both strains; %RSE 178%
    lbmaxb <- fixed(log(10^7.07)); label("Log maximum bacterial density in bladder (CFU/organ)")                           # Table 1: Bmax,b = 7.07 log10 CFU/organ, Fixed
    linocb <- log(10^4.42); label("Log initial bacterial density in bladder at 6 h post-inoculation (CFU/organ)")          # Table 1: Inocb = 4.42 log10 CFU/organ (95% CI 3.77-4.96)

    gam <- fixed(1); label("Power on the MIC-normalized drug effect, kidneys and bladder (unitless)")                      # Table 1: gamma = 1, Fixed (kidneys and bladder)

    # ---- Residual error ---------------------------------------------------
    propSd <- fixed(0.49); label("Proportional residual error on the mouse plasma concentration (fraction)")               # Sou 2021 Table 1, mouse: ERR = 49%. Not re-estimated in Hernandez-Lozano 2025
    # "For residual unexplained variability, additive error models were used"
    # (Materials and methods, Data analysis and software), but no sigma value
    # for the cfu observations is reported in the paper or the supplement.
    # Held at 0 rather than invented; see the vignette Errata.
    addSd_cfuKidney <- fixed(0); label("Additive residual SD on the kidney log10 bacterial count (log10 CFU/organ); magnitude not reported")   # Methods: additive residual error model; no value reported
    addSd_cfuBladder <- fixed(0); label("Additive residual SD on the bladder log10 bacterial count (log10 CFU/organ); magnitude not reported") # Methods: additive residual error model; no value reported
  })

  model({
    # ---- 1. Strain-specific parameter selection --------------------------
    if (BACT == 591) {
      kg1k <- exp(lkg1kEn591)
      redukgk <- redukgkEn591
      kdk <- exp(lkdkEn591)
      kg1b <- exp(lkg1bEn591)
      redukgb <- redukgbEn591
      kdb <- exp(lkdbEn591)
      micKidney <- micEn591ph74
      micBladder <- micEn591ph6
    } else {
      kg1k <- exp(lkg1kAtcc)
      redukgk <- redukgkAtcc
      kdk <- exp(lkdkAtcc)
      kg1b <- exp(lkg1bAtcc)
      redukgb <- redukgbAtcc
      kdb <- exp(lkdbAtcc)
      micKidney <- micAtccph74
      micBladder <- micAtccph6
    }
    kg2k <- kg1k * (1 - redukgk / 100)
    kg2b <- kg1b * (1 - redukgb / 100)

    bmaxk <- exp(lbmaxk)
    bmaxb <- exp(lbmaxb)

    # ---- 2. Mouse PK (Sou 2021 Eqs. 1-3) ---------------------------------
    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka30) * (DOSE_APRAMYCIN_MGKG / 30)^powka
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    # The PKPD model assumed that the unbound apramycin concentration in the
    # central compartment corresponds to the site-of-action concentration in
    # both organs; a PBPK cross-check showed kidney interstitial and plasma
    # profiles to be comparable (Results, Figure S4a).
    cu <- fu * Cc

    # ---- 3. Drug effect, kidney (pH 7.4) ---------------------------------
    cmick <- cu / micKidney
    kdrug1k <- slopeSk * cmick^gam
    kdrug2k <- slopeRk * cmick^gam
    # kada carries units of 1/h (Table 1), so the multiplier must be the
    # dimensionless MIC-normalized concentration for the transfer to be
    # "drug-driven" as the Methods describe it.
    kada1to2k <- kadak * cmick

    btotk <- bact_s1k + bact_d1k + bact_s2k + bact_d2k
    ksr1k <- (kg1k - kdk) * btotk / bmaxk
    ksr2k <- (kg2k - kdk) * btotk / bmaxk

    # ---- 4. Drug effect, bladder (pH 6) ----------------------------------
    cmicb <- cu / micBladder
    kdrug1b <- slopeSb * cmicb^gam
    kdrug2b <- slopeRb * cmicb^gam
    kada1to2b <- kadab * cmicb

    btotb <- bact_s1b + bact_d1b + bact_s2b + bact_d2b
    ksr1b <- (kg1b - kdb) * btotb / bmaxb
    ksr2b <- (kg2b - kdb) * btotb / bmaxb

    # ---- 5. Bacterial system in the kidneys (Figure 1) -------------------
    d/dt(bact_s1k) <- kg1k * bact_s1k - (kdk + kdrug1k) * bact_s1k -
      ksr1k * bact_s1k - kada1to2k * bact_s1k
    d/dt(bact_d1k) <- ksr1k * bact_s1k - kdk * bact_d1k
    d/dt(bact_s2k) <- kg2k * bact_s2k - (kdk + kdrug2k) * bact_s2k -
      ksr2k * bact_s2k + kada1to2k * bact_s1k
    d/dt(bact_d2k) <- ksr2k * bact_s2k - kdk * bact_d2k

    # ---- 6. Bacterial system in the bladder (Figure 1) -------------------
    d/dt(bact_s1b) <- kg1b * bact_s1b - (kdb + kdrug1b) * bact_s1b -
      ksr1b * bact_s1b - kada1to2b * bact_s1b
    d/dt(bact_d1b) <- ksr1b * bact_s1b - kdb * bact_d1b
    d/dt(bact_s2b) <- kg2b * bact_s2b - (kdb + kdrug2b) * bact_s2b -
      ksr2b * bact_s2b + kada1to2b * bact_s1b
    d/dt(bact_d2b) <- ksr2b * bact_s2b - kdb * bact_d2b

    # ---- 7. Initial conditions -------------------------------------------
    # Model time zero is 6 h after bacterial inoculation, the first sampling
    # time, at which Inoc_k and Inoc_b were estimated. All bacteria start in
    # the growing state of the main subpopulation: kada was included "rather
    # than assuming initial pre-existing percentages of each subpopulation".
    bact_s1k(0) <- exp(linock)
    bact_s1b(0) <- exp(linocb)

    # ---- 8. Observations --------------------------------------------------
    cfuKidney <- log10(btotk + 1e-6)
    cfuBladder <- log10(btotb + 1e-6)

    Cc ~ prop(propSd)
    cfuKidney ~ add(addSd_cfuKidney)
    cfuBladder ~ add(addSd_cfuBladder)
  })
}
