Li_2025_doxycycline_duck_ff40 <- function() {
  description <- "Preclinical (duck, Tadorna tadorna). One-compartment first-order-absorption plasma PK model for intramuscular doxycycline in Riemerella anatipestifer-infected ducks, coupled to the inhibitory sigmoid Emax exposure-response for the 24 h change in bacterial count when doxycycline is co-administered with florfenicol at 40 mg/kg. Li 2025 fitted plasma doxycycline with WinNonlin 6.1 using a one-compartment model with first-order absorption (Methods, 'Pharmacokinetics (PK) of DOX in RA-infected ducks'); Table 1 reports mean T1/2ka = 0.60 h, mean T1/2kel = 11.21 h and mean Cl/F = 0.40 L/h/kg, from which V/F = (Cl/F)/kel = 6.47 L/kg. Doses and all volume/clearance terms are per kg bodyweight, so the dosing amount is in mg/kg and concentrations are in ug/mL. The exposure-response is the inhibitory sigmoid Emax model of Methods, 'PK and PD analyses', indexed on the paper's primary and best-correlating PK/PD index AUC24h/MIC (R^2 = 0.912, versus 0.913 for Cmax/MIC and 0.566 for %T>MIC); Table 3 gives the florfenicol 40 mg/kg arm. Two sign conventions in the source must be read carefully and are handled here (see the vignette): the row Li 2025 labels 'Emax' is the change in the drug-free model group (bacterial GROWTH of 0.06 log10 CFU/mL), and the row labelled 'E0' is the maximum antibacterial effect (a REDUCTION of 4.76 log10 CFU/mL), which is the reverse of the usual naming and the reverse of the same group's Chen 2023 tilmicosin paper. The orientation is fixed numerically on the companion florfenicol 20 mg/kg arm, where substituting Table 2 returns exactly a 3.000 log10 CFU/mL reduction at each of that arm's three reported breakpoints. NOTE that Table 3's AUC24h/MIC column does NOT reproduce its own reported breakpoint: substituting Table 3 as printed returns a 2.474 log10 CFU/mL reduction at 19.98 h rather than 3.000, and an E0 of 5.761 would be required to return 3.000. The printed value 4.76 is packaged here unchanged because parameters are never tuned to hit a validation target; the 19.98 h breakpoint is corroborated three times in the paper (Abstract, Conclusions, and the 12.76 mg/kg dose prediction) whereas 4.76 appears once, so a digit slip in Table 3 is the most likely explanation. See the vignette Errata. The signed 24 h change in bacterial count is carried as the state dlog10cfu, which starts at 0 and reaches the predicted change at t = 24 h; no absolute starting bacterial density is packaged because Li 2025 reported only changes. Li 2025 fitted naive individual/mean profiles in WinNonlin and reported neither between-subject variability nor a residual error magnitude, so there are no eta parameters and propSd is FIXED at 0; the model is intended for typical-value simulation. The companion florfenicol 20 mg/kg parameterisation is Li_2025_doxycycline_duck_ff20, and the lung and liver matrices are Li_2025_doxycycline_duck_lung and Li_2025_doxycycline_duck_liver."
  reference <- paste(
    "Li FL, He CY, Chen HY, Cheng SM, Liu Y, Ding HZ, Zhang HL. (2025).",
    "In vivo pharmacokinetic/pharmacodynamic relationship of florfenicol in",
    "combination with doxycycline against Riemerella anatipestifer in ducks",
    "and the effect upon resistance development.",
    "Poultry Science 104:104922.",
    "doi:10.1016/j.psj.2025.104922.",
    sep = " "
  )
  vignette <- "Li_2025_doxycycline_florfenicol_ducks"
  units <- list(
    time = "h",
    dosing = "mg/kg (doxycycline; all disposition terms are per kg bodyweight)",
    concentration = "ug/mL (doxycycline, Cc); log10 CFU/mL (bacterial count change, dlog10cfu)"
  )

  # auc_dox    -- cumulative doxycycline plasma AUC (ug*h/mL). Diagnostic only:
  #   the exposure driving the PD is the closed form dose/(Cl/F) below, which
  #   this state converges to. Accumulator idiom follows Beredaki_2023_micafungin_clsi.
  # t_gt_mic   -- cumulative time (h) with plasma doxycycline above the MIC; the
  #   numerator of the paper's %T>MIC index, emitted for the vignette.
  # dlog10cfu  -- SIGNED change in bacterial count (log10 CFU/mL) accrued since
  #   dosing; negative = kill. Starts at 0 and reaches the Emax-model prediction
  #   at t = 24 h, the window Li 2025 actually counted bacteria over.
  paper_specific_compartments <- c("auc_dox", "t_gt_mic", "dlog10cfu")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot     = list(analyte = "doxycycline", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central   = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE),
    auc_dox   = list(analyte = "doxycycline", units = "ug*h/mL", specimen = "plasma", verified = TRUE),
    t_gt_mic  = list(analyte = "doxycycline", units = "h", specimen = "plasma", verified = TRUE),
    dlog10cfu = list(analyte = "Riemerella anatipestifer", units = "log10 CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "duck (common shelduck, Tadorna tadorna), 7 days old",
    n_subjects     = 360L,
    n_studies      = 1L,
    age_median     = "7 days at the start of the study",
    weight_range   = "130-150 g",
    sex_female_pct = NA_real_,
    organism       = paste(
      "Riemerella anatipestifer strain CVCC3857 (Chinese Veterinary Microorganism",
      "Culture Collection Center, Beijing). Li 2025 Results, 'MIC and MPC of FF and",
      "DOX against RA': doxycycline MIC 1 ug/mL and MPC 8 ug/mL; florfenicol MIC",
      "1 ug/mL and MPC 4 ug/mL. Table 4 independently lists the doxycycline MIC of",
      "the original 3857 strain as 1 ug/mL. The Abstract instead states a",
      "doxycycline MIC of 2 ug/mL; the value of 1 ug/mL is the one the paper's own",
      "dose predictions were computed with (see the vignette Errata). A second,",
      "less susceptible strain RA38 (florfenicol MIC 4 ug/mL, doxycycline MIC",
      "2 ug/mL) was used only for the twice-daily dosing experiment, for which no",
      "exposure-response model was fitted."
    ),
    disease_state  = paste(
      "Experimental systemic Riemerella anatipestifer infection established by",
      "intraperitoneal injection of 10^9 CFU/mL; the target bacterial load was",
      "reached 12 h after inoculation and drug was given at that point"
    ),
    dose_range     = paste(
      "Doxycycline 1, 2.5, 5, 10 and 20 mg/kg as a single intramuscular injection",
      "into the thigh (PK, five groups of 72 ducks). Pharmacodynamics: florfenicol",
      "20 or 40 mg/kg combined with doxycycline 1, 2.5, 5, 10 or 20 mg/kg, plus",
      "florfenicol 20 and 40 mg/kg and doxycycline 10 and 20 mg/kg monotherapy",
      "arms and an untreated model group (eight ducks per group)"
    ),
    sampling       = paste(
      "PK: blood, lung and liver at 0.5, 1, 2, 4, 6, 8, 12, 24 and 36 h after",
      "dosing, assayed by HPLC-MS/MS (LLOQ 0.005 ug/mL in plasma and lung,",
      "0.01 ug/g in liver). PD: bacterial load in blood, heart, liver, spleen,",
      "lung, kidney and brain counted 24 h after administration"
    ),
    regions        = "China (Fuyang Normal University, Anhui; South China Agricultural University, Guangzhou)",
    notes          = paste(
      "Ducks were obtained from a commercial farm in Guangxi, China. The 360",
      "subjects counted here are the PK cohort (five dose groups of 72); a further",
      "120 ducks (15 groups of 8) entered the single-dose PD study and 72 ducks",
      "(9 groups of 8) the twice-daily RA38 study. Doxycycline plasma protein",
      "binding was measured by equilibrium dialysis at 0.1, 1 and 10 ug/mL as",
      "37.84%, 29.92% and 44.33% (mean 37.36%), giving fu = 0.6264. Li 2025",
      "reports parameters at each of the five dose levels as well as the means",
      "used here; apparent clearance rose monotonically with dose (Cl/F 0.29 to",
      "0.54 L/h/kg from 1 to 20 mg/kg), which the paper does not comment on. See",
      "the vignette for the per-dose reproduction of Table 1 and for the finding",
      "that the tabulated 'AUC24h' values are AUC0-infinity, not a 0-24 h partial",
      "area."
    )
  )

  ini({
    # =================================================================
    # Plasma pharmacokinetics -- Li 2025 Table 1 (plasma, "Mean" column)
    # Methods, "Pharmacokinetics (PK) of DOX in RA-infected ducks":
    # "DOX concentration-time data were analyzed using a one-compartmental
    # model with a first-order absorption model by employing WinNonlin 6.1".
    # All disposition terms are per kg bodyweight, as reported.
    # =================================================================
    lka <- log(log(2) / 0.60)
    label("Log first-order absorption rate constant ka (1/h)")  # Li 2025 Table 1, plasma mean T1/2ka = 0.60 +/- 0.20 h; ka = ln(2)/T1/2ka

    lcl <- log(0.40)
    label("Log apparent clearance Cl/F (L/h/kg)")  # Li 2025 Table 1, plasma mean Cl/F = 0.40 +/- 0.08 L/h/kg (also quoted in Methods, "Dose calculations")

    # V/F is not tabulated. It is recovered from the two parameters that are,
    # using the one-compartment identity V/F = (Cl/F)/kel with kel = ln(2)/T1/2kel.
    # This reproduces Table 1's Tmax and Cmax at every dose level to two decimal
    # places (see the vignette source-trace table), which validates the recovery.
    lvc <- log(0.40 / (log(2) / 11.21))
    label("Log apparent central volume of distribution V/F (L/kg)")  # derived: Li 2025 Table 1 plasma mean Cl/F = 0.40 L/h/kg and mean T1/2kel = 11.21 +/- 0.99 h

    # =================================================================
    # Plasma protein binding
    # =================================================================
    # FIXED: a measured physicochemical property, not an estimated parameter.
    # Li 2025 used it to convert total exposure into the free exposure that
    # drives the Toutain dose equation in Methods, "Dose calculations".
    fu <- fixed(1 - 0.3736)
    label("Unbound fraction of doxycycline in duck plasma (measured)")  # Li 2025 Results, "PK of DOX in ducks": binding 37.84%, 29.92%, 44.33% at 0.1, 1, 10 ug/mL, mean 37.36%; Methods, "Dose calculations": "fu is the unbound fraction, which was determined to be 62.64%"

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # The PK/PD index driving the sigmoid is AUC24h/MIC, so the MIC is required
    # numerically. FIXED because it is a measured property of the challenge
    # strain, not an estimated parameter. Change it to apply the model to an
    # isolate with a different susceptibility.
    mic <- fixed(1)
    label("Doxycycline MIC against R. anatipestifer CVCC3857 (ug/mL;, measured)")  # Li 2025 Results, "MIC and MPC Of FF and DOX against RA": "The MIC and MPC of DOX against the RA strain CVCC3857 were 1 ug/mL and 8 ug/mL"; Table 4 row "3857(original)", doxycycline column = 1. The Abstract's "MIC of DOX = 2 ug/mL" is inconsistent with both and with the paper's own dose predictions -- see the vignette Errata

    # =================================================================
    # Inhibitory sigmoid Emax exposure-response
    # Li 2025 Methods, "PK and PD analyses" and Table 3 (florfenicol 40 mg/kg)
    # =================================================================
    # Li 2025 prints the equation as a formula image that the text glosses as
    #   "E is the antibacterial effect; Emax is the change in the model group
    #    (absence of drugs); E0 is the maximum antibacterial effect; Ce is the
    #    PK/PD index; EC50 is the corresponding PK/PD index that produces a 50%
    #    reduction in the maximum antibacterial effect; N is the Hill coefficient"
    # -- i.e. the labels "Emax" and "E0" carry the OPPOSITE roles to their usual
    # meaning and to the same group's Chen 2023 tilmicosin paper. Written as a
    # reduction magnitude the fitted model is therefore
    #   R = e0 + (emax - e0) * Ce^N / (EC50^N + Ce^N)
    # with e0 = -(Li 2025 "Emax" row) because the drug-free control GREW, and
    # emax = (Li 2025 "E0" row). This orientation is confirmed numerically on
    # the companion florfenicol 20 mg/kg arm, where substituting Table 2 returns
    # R = 3.000 at each of that table's three reported "Decrease in 3 Log10
    # CFU/mL" breakpoints. For THIS table the Cmax/MIC column also returns
    # R = 3.024 at its reported breakpoint of 2.20, but the AUC24h/MIC column as
    # printed returns only R = 2.474 at its reported breakpoint of 19.98 h.
    e0 <- -0.06
    label("Change in bacterial count in the untreated model group over 24 h, as a reduction (log10 CFU/mL; negative = growth)")  # Li 2025 Table 3 row "E max (Log 10 CFU/mL)" = 0.06, which the Methods define as the change in the drug-free model group; sign flipped because Results report the model group GREW by 0.04-0.93 log10 CFU/mL. Kept on the natural scale because it is a signed net change

    lemax <- log(4.76)
    label("Log maximum antibacterial effect over 24 h Emax (log10 CFU/mL reduction)")  # Li 2025 Table 3 row "E 0 (Log 10 CFU/mL)" = 4.76, which the Methods define as the maximum antibacterial effect. Packaged as printed; it does not reproduce the paper's own 19.98 h breakpoint (5.761 would) -- see the vignette Errata

    lec50 <- log(16.97)
    label("Log AUC24h/MIC producing 50% of the maximal effect EC50 (h)")  # Li 2025 Table 3 row "EC 50", AUC24h/MIC column = 16.97

    lhill <- log(0.63)
    label("Log Hill coefficient N of the exposure-response relationship (unitless)")  # Li 2025 Table 3 row "Hill's slope", AUC24h/MIC column = 0.63

    # =================================================================
    # Residual error
    # =================================================================
    # Li 2025 fitted mean profiles in WinNonlin and reported no between-subject
    # variability and no residual standard deviation for either the PK or the
    # exposure-response fit (only R^2 = 0.912 for the AUC24h/MIC sigmoid), so
    # the concentration residual SD is held at zero for deterministic
    # typical-value simulation. See the vignette Assumptions and deviations.
    propSd <- fixed(0)
    label("Proportional residual error on plasma doxycycline concentration (0; not reported in Li 2025)")  # Li 2025 reported no residual error model
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---------------------------------------------------------------
    # One-compartment model with first-order absorption
    # (Li 2025 Methods, "Pharmacokinetics (PK) of DOX in RA-infected ducks";
    # the printed closed form is C(t) = F*D/Vd * ka/(ka-kel) *
    # (exp(-kel*t) - exp(-ka*t)), whose ODE form is written out here so the
    # model also supports multiple-dose regimens such as the paper's
    # twice-in-24-h RA38 experiment). F is not separately identifiable from
    # intramuscular data alone, so cl and vc are the apparent terms Cl/F and
    # V/F exactly as Li 2025 reports them, and no f(depot) is applied.
    # ---------------------------------------------------------------
    kel <- cl / vc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Ccu <- fu * Cc

    # Diagnostic exposure accumulators. auc_dox converges to dose/(Cl/F);
    # t_gt_mic is the numerator of the paper's %T>MIC index.
    d/dt(auc_dox) <- Cc
    d/dt(t_gt_mic) <- (Cc > mic)

    # ---------------------------------------------------------------
    # PK/PD index. Li 2025 labels it AUC24h/MIC, but the tabulated values are
    # AUC0-infinity: Table 1's Cl/F equals dose divided by the tabulated "AUC24"
    # at all five dose levels to two decimal places, whereas the true 0-24 h
    # partial area is about 75% of that (see the vignette Errata). The index is
    # therefore formed in closed form as dose/(Cl/F)/MIC, which is exact for
    # this linear model and is available from t = 0 so the effect trajectory
    # below lands on the fitted 24 h endpoint. Idiom follows
    # Beredaki_2023_micafungin_clsi.
    #
    # The index is TOTAL, not unbound, plasma exposure. Two independent checks
    # fix this: the Discussion equates "AUC24 h/MIC was 37.11 h" at 20 mg/kg
    # with Table 1's total AUC of 37.11, and Table 2's Cmax/MIC breakpoint of
    # 1.72 exceeds the highest unbound Cmax reachable in the study
    # (fu * 1.86 = 1.17 at 20 mg/kg) so it can only be a total-drug index.
    # Note that the paper's own dose predictions nevertheless divide by fu (see
    # the vignette Errata), which is inconsistent with a total-drug breakpoint.
    # ---------------------------------------------------------------
    aucmic <- podo(depot) / cl / mic

    # Inhibitory sigmoid Emax model (Li 2025 Methods, "PK and PD analyses").
    # kill_log10 is the MAGNITUDE of the log10 CFU/mL reduction accrued over
    # the 24 h window (positive = bacterial kill); it equals e0 at zero
    # exposure (i.e. the control growth, as a negative reduction) and
    # approaches emax at saturating exposure.
    kill_log10 <- e0 + (emax - e0) * aucmic^hill / (ec50^hill + aucmic^hill)

    # Li 2025 counted bacteria only at 24 h after dosing, so the fitted effect
    # is a per-24-h quantity. Spreading it uniformly across the window makes
    # dlog10cfu equal exactly -kill_log10 at t = 24 h, matching the paper's
    # model at the time it was fitted, while remaining well defined for the
    # twice-daily regimen. No absolute starting density is packaged because
    # Li 2025 reported only changes.
    d/dt(dlog10cfu) <- -kill_log10 / 24
    dlog10cfu(0) <- 0

    Cc ~ prop(propSd)
  })
}
