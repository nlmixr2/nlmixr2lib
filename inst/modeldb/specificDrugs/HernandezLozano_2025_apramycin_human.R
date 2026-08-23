HernandezLozano_2025_apramycin_human <- function() {
  description <- "Translational prediction of apramycin efficacy against Escherichia coli in human complicated urinary tract infection. The published human population PK model of Zhao 2022 (mammillary four-compartment model with linear elimination, absolute eGFR proportional on clearance, allometric total body weight on all fixed effects, and 90% renal fractional excretion into a urine compartment) drives the unbound plasma concentration into the mouse-derived pharmacodynamic component of Hernandez-Lozano 2025. Kidneys and bladder are modelled as two independent organ systems, each with two bacterial subpopulations (1 = main, apramycin-susceptible; 2 = decreased susceptibility) and, within each, a growing drug-susceptible state S and a dormant drug-insusceptible state D entered at rate ksr = (kg - kd) * Btot / Bmax. Apramycin adds to the death rate of the S states through a power model normalized to the in vitro MIC, kdrug = Slope * (Cu/MIC)^gamma with gamma fixed to 1, and drives transfer from subpopulation 1 to subpopulation 2 at rate kada * (Cu/MIC). The kidney is assumed to be at pH 7.4 and the bladder at pH 6, so each organ normalizes by the MIC measured at its own pH. Initial bacterial densities default to 10^6 and 10^5 CFU per organ for kidneys and bladder as used in the paper's human simulations, in which a single 10.8 mg/kg 30 min intravenous infusion produced stasis in both organs for both strains. Sibling models: HernandezLozano_2025_apramycin_invitro (the underlying time-kill fit) and HernandezLozano_2025_apramycin_mouse (the same PD component driven by a mouse subcutaneous PK model)."
  reference <- "Hernandez-Lozano I, Aranzana-Climent V, Cao S, Matias C, Hansen JU, Liepinsh E, Hughes D, Hobbie SN, Vingsbo Lundberg C, Friberg LE. Model-informed drug development for antimicrobials: translational pharmacokinetic-pharmacodynamic modelling of apramycin to facilitate prediction of efficacious dose in complicated urinary tract infections. J Antimicrob Chemother. 2025 Feb 3;80(2):302-311. doi:10.1093/jac/dkae409. PMID: 39545353. PMCID: PMC11695905. Human simulation design (75 kg adult, creatinine clearance 120 mL/min, 30 min intravenous infusions of 0.3-30 mg/kg, initial densities 10^6 and 10^5 CFU per organ, pH 7.4 in kidney and pH 6 in bladder, unbound fraction 92.9%): Materials and methods, 'Prediction of human efficacy', and Figure 5. PD estimates: Table 1, section 'In vivo PD parameters'. Human PK parameters are NOT re-estimated here; they are taken from Zhao C, Chirkova A, Rosenborg S, et al. Population pharmacokinetics of apramycin from first-in-human plasma and urine data to support prediction of efficacious dose. J Antimicrob Chemother 2022;77:2718-28. doi:10.1093/jac/dkac225, Table 2 (the final 'Plasma+urine data' column) and its footnote a."
  vignette <- "HernandezLozano_2025_apramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Two independent organ PD systems, each replicating the in vitro
  # structure. Suffix k = kidneys, b = bladder. Within an organ,
  # subpopulation 1 is the main apramycin-susceptible population and 2 the
  # less-susceptible one; S is the growing drug-susceptible state and D the
  # dormant drug-insusceptible one (Figure 1). The names do not match the
  # canonical bact_<phenotype><lifecycle-digit> scheme because that trailing
  # digit indexes the Bulitta / Wicha replication life cycle, a different
  # distinction from this model's growing-versus-dormant states, and because
  # that scheme has no organ dimension.
  paper_specific_compartments <- c(
    "bact_s1k", "bact_d1k", "bact_s2k", "bact_d2k",
    "bact_s1b", "bact_d1b", "bact_s2b", "bact_d2b"
  )

  compartmentData <- list(
    central = list(analyte = "apramycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apramycin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "apramycin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral3 = list(analyte = "apramycin", units = "mg", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "apramycin", units = "mg", specimen = "urine", verified = TRUE),
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
        "591 = Escherichia coli EN591, an MDR rmtB clinical isolate; MIC 8 mg/L at pH 7.4 and 32 mg/L at pH 6. Panel (a) of Figure 5.",
        "700336 = Escherichia coli ATCC 700336 (EN1085), a trimethoprim/sulfamethoxazole-resistant urinary tract infection isolate; MIC 4 mg/L at pH 7.4 and 16 mg/L at pH 6. Panel (b) of Figure 5.",
        "There is no reference level. The drug-effect parameters SlopeS, SlopeR and kada are SHARED across strains in both organs (Table 1, Strain column 'Both'); strain differences in apramycin effect are carried entirely by the strain-specific MIC used to normalize the concentration. The in vivo growth and death rate constants, which do differ by strain, come from the mouse fit and are applied unchanged to the human prediction."
      )
    ),
    WT = list(
      description = "Total body weight, allometrically scaling every fixed-effect PK parameter of the upstream human population PK model.",
      units = "kg",
      type = "continuous",
      source_name = "TBW",
      reference_category = NULL,
      notes = "Zhao 2022 scaled all clearance parameters with exponent 0.75 and all volume parameters with exponent 1, both fixed, against a 70 kg reference. Hernandez-Lozano 2025 simulated a 75 kg adult ('Prediction of human efficacy'), which is the default here. The Zhao 2022 Phase I cohort was 30 healthy volunteers."
    ),
    CRCL = list(
      description = "Absolute (NOT body-surface-area-normalized) glomerular filtration rate estimated with the CKD-EPI equation and multiplied back by BSA/1.73, proportional on apramycin clearance.",
      units = "mL/min",
      type = "continuous",
      source_name = "absolute eGFR",
      reference_category = NULL,
      notes = paste(
        "RAW renal function in mL/min, not normalized to 1.73 m^2, matching the parameterization of the upstream model (precedent for the un-normalized form under this canonical: Delattre_2010_amikacin.R, Georges_2009_ceftazidime.R). Zhao 2022 Table 2 footnote a: CL = 5.54 x (absolute eGFR/124) x (TBW/70)^0.75; absolute eGFR is the CKD-EPI estimate corrected by body surface area. The reference value 124 mL/min is the typical value of the healthy-volunteer cohort (observed range 91.9-157 mL/min).",
        "Hernandez-Lozano 2025 describes the simulated subject as having a 'creatinine clearance of 120 mL/min'. The covariate the upstream model actually consumes is absolute eGFR, so 120 is entered here as absolute eGFR; the wording difference is recorded in the vignette Errata."
      )
    )
  )

  population <- list(
    species = "human (simulated)",
    n_subjects = NA_integer_,
    n_studies = 0L,
    weight_median = "75 kg (single simulated typical subject)",
    renal_function = "Absolute eGFR 120 mL/min (described in the source as a creatinine clearance of 120 mL/min)",
    disease_state = "Simulated complicated urinary tract infection with Escherichia coli EN591 or ATCC 700336, with kidneys assumed to be at pH 7.4 and bladder at pH 6",
    dose_range = "Single 30 min intravenous infusions of 0.3, 1.2, 3.6, 10.8 and 30 mg/kg, the dose levels of the apramycin Phase I trial; bacterial burden was followed for 48 h",
    initial_inoculum = "10^6 CFU per organ in kidneys and 10^5 CFU per organ in bladder, set from the mouse data",
    notes = paste(
      "The PD component is the mouse in vivo model of Table 1 applied unchanged; only the PK layer is replaced, by the human population PK model of Zhao 2022 developed from 30 healthy volunteers in the Phase I single-ascending-dose trial JUV18-01 (480 plasma and 179 urine observations).",
      "An unbound fraction of 92.9% in human plasma was applied. The assumption that unbound central-compartment concentration equals the site-of-action concentration was supported by a PBPK cross-check showing comparable plasma and kidney interstitial profiles in humans (Figure S4b); that PBPK model is not part of this file.",
      "A single dose of 10.8 mg/kg was predicted to be sufficient to reach stasis in both organs for both strains, supporting an 11 mg/kg daily dose. Figure 5 shows medians and 95% CIs of simulations that included inter-individual variability, residual error and parameter uncertainty (relative standard errors) in all model parameters; the parameter-uncertainty component is not carried in this file."
    )
  )

  ini({
    # ---- Human population PK (Zhao 2022 Table 2, plasma+urine column) ----
    # Not re-estimated in Hernandez-Lozano 2025; carried in from the cited
    # upstream model, so every structural PK parameter is fixed here.
    lcl <- fixed(log(5.54)); label("Log apramycin clearance at the reference renal function and body weight (L/h)") # Zhao 2022 Table 2, plasma+urine: CL = 5.54 L/h (RSE 2.38%)
    lvc <- fixed(log(8.61)); label("Log central volume of distribution (L)")                                        # Zhao 2022 Table 2, plasma+urine: Vc = 8.61 L (RSE 5.27%)
    lq <- fixed(log(0.127)); label("Log intercompartmental clearance to the first peripheral compartment (L/h)")    # Zhao 2022 Table 2, plasma+urine: Q2 = 0.127 L/h (RSE 7.29%)
    lvp <- fixed(log(2.29)); label("Log volume of the first peripheral compartment (L)")                            # Zhao 2022 Table 2, plasma+urine: V2 = 2.29 L (RSE 3.58%)
    lq2 <- fixed(log(13.6)); label("Log intercompartmental clearance to the second peripheral compartment (L/h)")   # Zhao 2022 Table 2, plasma+urine: Q3 = 13.6 L/h (RSE 7.34%)
    lvp2 <- fixed(log(2.81)); label("Log volume of the second peripheral compartment (L)")                          # Zhao 2022 Table 2, plasma+urine: V3 = 2.81 L (RSE 13.2%)
    lq3 <- fixed(log(1.01)); label("Log intercompartmental clearance to the third peripheral compartment (L/h)")    # Zhao 2022 Table 2, plasma+urine: Q4 = 1.01 L/h (RSE 13.4%)
    lvp3 <- fixed(log(2.38)); label("Log volume of the third peripheral compartment (L)")                           # Zhao 2022 Table 2, plasma+urine: V4 = 2.38 L (RSE 5.71%)
    fe <- fixed(0.900); label("Fraction of eliminated apramycin excreted into urine (unitless)")                    # Zhao 2022 Table 2, plasma+urine: Fe = 0.900 (RSE 2.65%)
    crclRef <- fixed(124); label("Reference absolute eGFR for apramycin clearance (mL/min)")                        # Zhao 2022 Table 2 footnote a: CL = 5.54 x (absolute eGFR/124) x (TBW/70)^0.75
    e_wt_cl <- fixed(0.75); label("Allometric exponent of total body weight on clearance parameters (unitless)")    # Zhao 2022 Methods: "exponents being fixed to 0.75 for clearance parameters, and 1 for volume parameters"
    e_wt_vc <- fixed(1); label("Allometric exponent of total body weight on volume parameters (unitless)")          # Zhao 2022 Methods: as above
    fu <- fixed(0.929); label("Fraction of apramycin unbound in human plasma (unitless)")                           # Methods, Prediction of human efficacy: "An unbound fraction of 92.9% in human plasma was applied"

    # ---- Strain MICs used to normalize the drug effect --------------------
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
    lkdkAtcc <- fixed(log(0.179)); label("Log kidney bacterial death rate constant, ATCC 700336 (1/h)")                   # Table 1: kd,k = 0.179, Fixed for ATCC700336
    slopeSk <- 9.05; label("Apramycin effect on the susceptible subpopulation in kidney (1/h)")                           # Table 1: SlopeSk = 9.05 (95% CI 5.41-13.11), both strains
    slopeRk <- 0.066; label("Apramycin effect on the resistant subpopulation in kidney (1/h)")                            # Table 1: SlopeRk = 0.066 (95% CI 0.026-0.109), both strains
    kadak <- 0.031; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation in kidney (1/h)") # Table 1: kada,k = 0.031 (95% CI 0.006-0.062), both strains
    lbmaxk <- fixed(log(10^6.49)); label("Log maximum bacterial density in kidneys (CFU/organ)")                          # Table 1: Bmax,k = 6.49 log10 CFU/organ, Fixed
    linock <- fixed(log(10^6)); label("Log initial bacterial density in kidneys (CFU/organ)")                             # Methods, Prediction of human efficacy: "initial bacterial densities were set to 10^6 and 10^5 per organ for kidneys and bladder, respectively"

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
    linocb <- fixed(log(10^5)); label("Log initial bacterial density in bladder (CFU/organ)")                              # Methods, Prediction of human efficacy: "initial bacterial densities were set to 10^6 and 10^5 per organ for kidneys and bladder, respectively"

    gam <- fixed(1); label("Power on the MIC-normalized drug effect, kidneys and bladder (unitless)")                      # Table 1: gamma = 1, Fixed (kidneys and bladder)

    # ---- Inter-individual variability (Zhao 2022 Table 2, plasma+urine) --
    # Reported as CV%; converted with omega^2 = log(CV^2 + 1) per the table
    # footnote c ("CV, coefficient of variation, calculated by sqrt(exp(OMEGA)
    # - 1)"). Off-diagonals are corr x sqrt(var_A x var_B) per footnote d.
    # CL~Vc 0.52, CL~V3 -0.40, Vc~V3 -0.92, so etalcl / etalvc / etalvp2 form
    # one 3x3 block; V4 is independent.
    etalcl + etalvc + etalvp2 ~ c(
      0.0205237,
      0.0240203, 0.1039655,
      -0.0323211, -0.1673124, 0.3181215
    ) # Zhao 2022 Table 2: IIV CL 14.4% CV -> 0.0205237; Vc 33.1% -> 0.1039655; V3 61.2% -> 0.3181215; correlations 0.52, -0.40, -0.92
    etalvp3 ~ 0.0191367 # Zhao 2022 Table 2: IIV V4 13.9% CV -> log(1 + 0.139^2) = 0.0191367; no correlation reported with the other three

    # ---- Residual error ---------------------------------------------------
    propSd <- fixed(0.0879); label("Proportional residual error on the human plasma concentration (fraction)") # Zhao 2022 Table 2, plasma+urine: Prop plasma = 8.79% (RSE 8.96%)
    # "For residual unexplained variability, additive error models were used"
    # (Hernandez-Lozano Materials and methods), but no sigma value for the cfu
    # observations is reported in the paper or the supplement. Held at 0
    # rather than invented; see the vignette Errata.
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

    # ---- 2. Human PK (Zhao 2022 Table 2 and footnote a) ------------------
    cl <- exp(lcl + etalcl) * (CRCL / crclRef) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    q2 <- exp(lq2) * (WT / 70)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc
    q3 <- exp(lq3) * (WT / 70)^e_wt_cl
    vp3 <- exp(lvp3 + etalvp3) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    d/dt(central) <- -(kel + k12 + k13 + k14) * central +
      k21 * peripheral1 + k31 * peripheral2 + k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3
    # Renal fractional excretion: fe of the eliminated drug reaches urine.
    d/dt(urine) <- fe * kel * central

    Cc <- central / vc
    cu <- fu * Cc

    # ---- 3. Drug effect, kidney (pH 7.4) ---------------------------------
    cmick <- cu / micKidney
    kdrug1k <- slopeSk * cmick^gam
    kdrug2k <- slopeRk * cmick^gam
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
    # Treatment starts at model time zero, with all bacteria in the growing
    # state of the main subpopulation.
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
