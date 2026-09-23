Jiang_2016_harmaline_meodmt_mouse <- function() {
  description <- "Preclinical (mouse). Mechanism-based PK/PD model of core body temperature for the serotonergic drug-drug interaction between the MAO-A inhibitor harmaline and the 5-HT receptor agonist 5-MeO-DMT, in wild-type and CYP2D6-humanized transgenic mice. Harmaline two-compartment PK with linear murine plus CYP2D6 elimination; 5-MeO-DMT two-compartment PK with three parallel Michaelis-Menten routes (MAO-A, other murine, O-demethylation), two of them competitively inhibited by harmaline. Thermoregulation is an indirect-response model with adaptive feedback, a handling/injection stress signal, harmaline stimulation of heat loss via 5-HT1A, and delayed 5-MeO-DMT stimulation of thermogenesis via 5-HT2A through three transit compartments."
  reference <- "Jiang XL, Shen HW, Mager DE, Schmidt S, Yu AM. Development of a mechanism-based pharmacokinetic/pharmacodynamic model to characterize the thermoregulatory effects of serotonergic drugs in mice. Acta Pharm Sin B. 2016;6(5):492-503. doi:10.1016/j.apsb.2016.07.007"
  vignette <- "Jiang_2016_harmaline_meodmt_mouse"

  # Neither drug is a metabolite of the other, so both carry an explicit drug
  # suffix rather than one of them taking the bare canonical names. `stress` is
  # the paper's handling/injection signal (Fig. 1, lower panel); `temp` is the
  # core-body-temperature turnover state. `temp` is held paper-specific rather
  # than promoted to a canonical compartment: the standing operator ruling is
  # that a compartment canonical needs a second independent paper, and `temp`
  # is a collision-prone token besides. See the vignette Assumptions and
  # deviations for the residual convention warning this leaves.
  paper_specific_compartments <- c(
    "depot_harmaline",
    "central_harmaline",
    "peripheral1_harmaline",
    "depot_meodmt",
    "central_meodmt",
    "peripheral1_meodmt",
    "stress",
    "temp"
  )

  compartmentData <- list(
    depot_harmaline = list(analyte = "harmaline", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central_harmaline = list(analyte = "harmaline", units = "umol/kg", specimen = "serum", verified = TRUE),
    peripheral1_harmaline = list(analyte = "harmaline", units = "umol/kg", specimen = "tissue", verified = TRUE),
    depot_meodmt = list(analyte = "5-MeO-DMT", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central_meodmt = list(analyte = "5-MeO-DMT", units = "umol/kg", specimen = "serum", verified = TRUE),
    peripheral1_meodmt = list(analyte = "5-MeO-DMT", units = "umol/kg", specimen = "tissue", verified = TRUE),
    stress = list(
      analyte = "handling/injection stress signal",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "5-HT2A transduction signal",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "5-HT2A transduction signal",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "5-HT2A transduction signal",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    temp = list(
      analyte = "core body temperature (implanted intraperitoneal telemetry)",
      units = "degC",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # Amounts are per kg body weight because every published volume, clearance and
  # Vmax is weight-normalised (Methods p.494). Dose records are therefore
  # umol/kg; see the vignette for the mg/kg -> umol/kg conversion.
  units <- list(time = "min", dosing = "umol/kg", concentration = "umol/L")

  covariateData <- list(
    CYP2D6_TG = list(
      description = "CYP2D6-humanized transgenic animal indicator; 1 = Tg-CYP2D6 mouse expressing functional human CYP2D6, 0 = wild-type mouse lacking it",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type FVB/N mouse; murine elimination pathways only)",
      notes = "Switches on the CYP2D6-mediated harmaline clearance arm (CLCYP2D6-H, absent in wild-type) and scales the 5-MeO-DMT O-demethylation Vmax by (1 + fcyp2d6_odm_meodmt). Both genotypes share one set of PD parameters: Discussion p.500 states the two strains share a genetic background, so the PD difference between them is entirely PK-driven.",
      source_name = "Tg-CYP2D6 vs wild-type (Methods 2.2)"
    ),
    INJ_REPEAT = list(
      description = "Repeat-handling indicator; 1 = the animal received a second handling/injection event after the index dose, 0 = a single injection",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single i.p. injection)",
      notes = "Selects the harmaline heat-loss sensitivity constant: ks_harmaline_single when 0, ks_harmaline_ddi when 1. Results p.497 and Discussion p.500: the second handling/injection in the DDI studies reduced kS-H by more than 50% (0.0347 -> 0.0136 L/umol). The trigger is the handling event itself, not the identity of what was injected - the paper reports the same reduction whether saline or 5-MeO-DMT was given 15 min after harmaline. Distinct from CONMED_HARMALINE, which keys on whether harmaline was present at all.",
      source_name = "single-dose vs DDI dosing regimen (Methods 2.4, Table 1)"
    ),
    CONMED_HARMALINE = list(
      description = "Concomitant harmaline coadministration indicator; 1 = 5-MeO-DMT given to a harmaline-pretreated animal, 0 = 5-MeO-DMT given alone or after vehicle",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no harmaline pretreatment)",
      notes = "Selects the 5-MeO-DMT potency on thermogenesis: ec50_meodmt_nonddi when 0, ec50_meodmt_ddi when 1 (Table 1 definitions, 1.88 vs 0.496 umol/L - the ~4-fold potentiation that is the paper's headline result). This is a pharmacodynamic interaction term only; the pharmacokinetic arm of the harmaline interaction is carried mechanistically by the two competitive inhibition constants and needs no covariate.",
      source_name = "pretreated with harmaline vs with vehicle (Table 1)"
    ),
    DOSE_HARMALINE_MGKG = list(
      description = "Administered harmaline dose level per kg body weight",
      units = "mg/kg",
      type = "continuous",
      reference_category = "n/a - used as a dose-group selector for the dose-dependent bioavailability, not as a normalized continuous term",
      notes = "Harmaline i.p. bioavailability rises with dose (Methods p.494: 34.6 / 74.2 / 90.3 percent at 2 / 5 / 15 mg/kg in wild-type mice and 36.8 / 68.8 / 80.1 percent in Tg-CYP2D6 mice). The model reproduces those six published values exactly at the three published dose levels and interpolates linearly between them; outside 2-15 mg/kg the nearest published value is held. Values at doses other than 2, 5 and 15 mg/kg are an encoding choice, not a published quantity - see the vignette Errata.",
      source_name = "harmaline dose (Methods 2.4)"
    )
  )

  population <- list(
    species = "mouse (male FVB/N wild-type and Tg-CYP2D6 humanized)",
    n_per_group = "26 per genotype for the saline / baseline arms (Fig. 2); 14 per genotype for the single-drug arms (Fig. 3); 11 wild-type and 12 Tg-CYP2D6 for the combination arms (Figs. 4, 6); 4-14 per genotype for the external validation arms (Fig. 5)",
    n_studies = 1L,
    weight_range = "25-35 g",
    sex_female_pct = 0,
    disease_state = "healthy; drug-induced thermoregulatory perturbation",
    dose_range = "harmaline 2-15 mg/kg i.p. and 5-MeO-DMT 2-20 mg/kg i.p., alone and in combination",
    regions = "United States (University at Buffalo, SUNY)",
    notes = "Telemetric core body temperature (Physiotel TA10TA-F20) sampled 10 times per minute and averaged over 5 or 10 min for modelling; ambient temperature 20 +/- 2 degC, 12 h light/dark cycle, observation window 10:30 a.m. to 3:30 p.m. Fitted in ADAPT V by naive-pooled maximum likelihood, so the model carries typical values only and no inter-individual variability."
  )

  ini({
    # ---------------------------------------------------------------------
    # Harmaline PK. Every value is fixed: Methods p.494 states 'All PK
    # parameters were fixed and linked to PD model', the estimates being
    # carried over from the companion PK-interaction study (Jiang 2013,
    # Drug Metab Dispos 41:975-86, reference 27).
    # ---------------------------------------------------------------------
    lka_harmaline <- fixed(log(0.307))
    label("Harmaline first-order absorption rate constant (ka-H, 1/min)")
    lvc_harmaline <- fixed(log(2.43))
    label("Harmaline central volume (VC-H, L/kg)")
    lvp_harmaline <- fixed(log(2.86))
    label("Harmaline peripheral volume (VP-H, L/kg)")
    lq_harmaline <- fixed(log(0.879))
    label("Harmaline distribution clearance (CLD-H, L/min/kg)")
    lcl_other_harmaline <- fixed(log(0.0962))
    label("Harmaline intrinsic murine clearance, both genotypes (CLother-H, L/min/kg)")
    lcl_cyp2d6_harmaline <- fixed(log(0.0608))
    label("Harmaline CYP2D6-mediated clearance, Tg-CYP2D6 only (CLCYP2D6-H, L/min/kg)")

    # Dose-dependent i.p. bioavailability, Methods p.494. Six published values.
    fdepot_harmaline_wt_d2 <- fixed(0.346)
    label("Harmaline bioavailability at 2 mg/kg, wild-type (FH, fraction)")
    fdepot_harmaline_wt_d5 <- fixed(0.742)
    label("Harmaline bioavailability at 5 mg/kg, wild-type (FH, fraction)")
    fdepot_harmaline_wt_d15 <- fixed(0.903)
    label("Harmaline bioavailability at 15 mg/kg, wild-type (FH, fraction)")
    fdepot_harmaline_tg_d2 <- fixed(0.368)
    label("Harmaline bioavailability at 2 mg/kg, Tg-CYP2D6 (FH, fraction)")
    fdepot_harmaline_tg_d5 <- fixed(0.688)
    label("Harmaline bioavailability at 5 mg/kg, Tg-CYP2D6 (FH, fraction)")
    fdepot_harmaline_tg_d15 <- fixed(0.801)
    label("Harmaline bioavailability at 15 mg/kg, Tg-CYP2D6 (FH, fraction)")

    # ---------------------------------------------------------------------
    # 5-MeO-DMT PK, also fixed from the same upstream study (Methods p.494).
    # ---------------------------------------------------------------------
    lka_meodmt <- fixed(log(0.0748))
    label("5-MeO-DMT first-order absorption rate constant (ka-M, 1/min)")
    lvc_meodmt <- fixed(log(0.460))
    label("5-MeO-DMT central volume (VC-M, L/kg)")
    lvp_meodmt <- fixed(log(2.29))
    label("5-MeO-DMT peripheral volume (VP-M, L/kg)")
    lq_meodmt <- fixed(log(0.301))
    label("5-MeO-DMT distribution clearance (CLD-M, L/min/kg)")
    fdepot_meodmt <- fixed(0.748)
    label("5-MeO-DMT i.p. bioavailability (FM, fraction)")

    lvmax_mao_meodmt <- fixed(log(2.69))
    label("5-MeO-DMT maximum metabolic rate by MAO-A (Vmax(M)-M, umol/min/kg)")
    lkm_mao_meodmt <- fixed(log(16.7))
    label("5-MeO-DMT Michaelis constant for MAO-A (Km(M)-M, umol/L)")
    lvmax_other_meodmt <- fixed(log(0.0610))
    label("5-MeO-DMT maximum metabolic rate by other murine pathways (Vmax(O)-M, umol/min/kg)")
    lkm_other_meodmt <- fixed(log(0.446))
    label("5-MeO-DMT Michaelis constant for other murine pathways (Km(O)-M, umol/L)")
    lvmax_odm_meodmt <- fixed(log(0.0334))
    label("5-MeO-DMT maximum metabolic rate by O-demethylation (Vmax(D)-M, umol/min/kg)")
    lkm_odm_meodmt <- fixed(log(1.27))
    label("5-MeO-DMT Michaelis constant for O-demethylation (Km(D)-M, umol/L)")
    fcyp2d6_odm_meodmt <- fixed(0.301)
    label("Fractional increment in 5-MeO-DMT O-demethylation capacity contributed by human CYP2D6 in Tg-CYP2D6 mice (fmCYP2D6(D)-M, fraction)")

    # Harmaline inhibition of the two 5-MeO-DMT routes it perpetrates on.
    # Fig. 1 draws both inhibition bars terminating on the Km labels, i.e.
    # competitive inhibition.
    lki_mao_harmaline <- fixed(log(0.048))
    label("Harmaline inhibition constant against MAO-A, literature value (Ki(M)-H, umol/L)")
    lki_odm_harmaline <- fixed(log(7.13))
    label("Harmaline inhibition constant against 5-MeO-DMT O-demethylase, from the upstream PK model fitting (Ki(D)-H, umol/L)")

    # ---------------------------------------------------------------------
    # Thermoregulation PD. These are the parameters this paper estimated;
    # every value and its precision is Table 1.
    # ---------------------------------------------------------------------
    lrbase <- log(35.8)
    label("Baseline core body temperature at the start of the observation window (Tempbasal(0), degC)")
    drift_rbase <- 0.00322
    label("Linear drift of baseline core body temperature with time (lambda, degC/min)")
    lkout <- log(0.0291)
    label("First-order rate constant of heat loss (kout, 1/min)")
    s0_stress <- 0.265
    label("Magnitude of the handling/injection stress signal on thermogenesis (S0, unitless)")
    lkel_stress <- log(0.103)
    label("First-order rate constant for loss of the stress signal (kel-S, 1/min)")
    ltau <- log(6.20)
    label("Mean transit time of each of the three 5-HT2A signal-transduction compartments (tau, min)")
    lemax_meodmt <- log(0.134)
    label("Maximum 5-MeO-DMT stimulation of thermogenesis (Smax-M, unitless)")
    lec50_meodmt_nonddi <- log(1.88)
    label("5-MeO-DMT concentration giving half-maximal thermogenesis when dosed alone or after vehicle (SC50-M-non-DDI, umol/L)")
    lec50_meodmt_ddi <- log(0.496)
    label("5-MeO-DMT concentration giving half-maximal thermogenesis when dosed after harmaline (SC50-M-DDI, umol/L)")
    lks_harmaline_single <- log(0.0347)
    label("Harmaline sensitivity constant on heat loss after a single injection (kS-H-single, L/umol)")
    lks_harmaline_ddi <- log(0.0136)
    label("Harmaline sensitivity constant on heat loss after a repeat-handling regimen (kS-H-DDI, L/umol)")

    # Eq. (10) defines the variance model as VARi = (sigma1 + sigma2 * Y)^2,
    # i.e. a combined additive-plus-proportional SD, but Table 1 tabulates only
    # the PD structural parameters and no sigma is printed anywhere in the
    # paper. The structure is preserved and both magnitudes are held at zero
    # rather than invented; see the vignette Errata.
    addSd_temp <- fixed(0)
    label("Additive residual SD on core body temperature, magnitude not published (degC)")
    propSd_temp <- fixed(0)
    label("Proportional residual SD on core body temperature, magnitude not published (fraction)")
  })

  model({
    # ---- harmaline disposition -------------------------------------------
    ka_harmaline <- exp(lka_harmaline)
    vc_harmaline <- exp(lvc_harmaline)
    vp_harmaline <- exp(lvp_harmaline)
    q_harmaline <- exp(lq_harmaline)
    cl_harmaline <- exp(lcl_other_harmaline) +
      exp(lcl_cyp2d6_harmaline) * CYP2D6_TG

    # Dose-dependent bioavailability: exact at the three published dose levels,
    # linear between them, nearest published value held outside 2-15 mg/kg.
    fd2 <- fdepot_harmaline_wt_d2 * (1 - CYP2D6_TG) +
      fdepot_harmaline_tg_d2 * CYP2D6_TG
    fd5 <- fdepot_harmaline_wt_d5 * (1 - CYP2D6_TG) +
      fdepot_harmaline_tg_d5 * CYP2D6_TG
    fd15 <- fdepot_harmaline_wt_d15 * (1 - CYP2D6_TG) +
      fdepot_harmaline_tg_d15 * CYP2D6_TG
    doseh <- DOSE_HARMALINE_MGKG
    fdepot_harmaline <- (doseh <= 2) * fd2 +
      (doseh > 2) * (doseh <= 5) * (fd2 + (fd5 - fd2) * (doseh - 2) / 3) +
      (doseh > 5) * (doseh <= 15) * (fd5 + (fd15 - fd5) * (doseh - 5) / 10) +
      (doseh > 15) * fd15

    Cc_harmaline <- central_harmaline / vc_harmaline
    Cp_harmaline <- peripheral1_harmaline / vp_harmaline

    d/dt(depot_harmaline) <- -ka_harmaline * depot_harmaline
    d/dt(central_harmaline) <- ka_harmaline * depot_harmaline -
      cl_harmaline * Cc_harmaline -
      q_harmaline * (Cc_harmaline - Cp_harmaline)
    d/dt(peripheral1_harmaline) <- q_harmaline *
      (Cc_harmaline - Cp_harmaline)
    f(depot_harmaline) <- fdepot_harmaline

    # ---- 5-MeO-DMT disposition -------------------------------------------
    ka_meodmt <- exp(lka_meodmt)
    vc_meodmt <- exp(lvc_meodmt)
    vp_meodmt <- exp(lvp_meodmt)
    q_meodmt <- exp(lq_meodmt)

    Cc_meodmt <- central_meodmt / vc_meodmt
    Cp_meodmt <- peripheral1_meodmt / vp_meodmt

    # Competitive inhibition by harmaline raises the apparent Km of the MAO-A
    # and O-demethylation routes (Fig. 1: both Ki arrows terminate on Km). The
    # third route, other murine elimination, carries no inhibition arrow.
    km_mao_app <- exp(lkm_mao_meodmt) *
      (1 + Cc_harmaline / exp(lki_mao_harmaline))
    km_odm_app <- exp(lkm_odm_meodmt) *
      (1 + Cc_harmaline / exp(lki_odm_harmaline))
    vmax_odm <- exp(lvmax_odm_meodmt) * (1 + fcyp2d6_odm_meodmt * CYP2D6_TG)

    elim_meodmt <- exp(lvmax_mao_meodmt) * Cc_meodmt /
      (km_mao_app + Cc_meodmt) +
      exp(lvmax_other_meodmt) * Cc_meodmt /
        (exp(lkm_other_meodmt) + Cc_meodmt) +
      vmax_odm * Cc_meodmt / (km_odm_app + Cc_meodmt)

    d/dt(depot_meodmt) <- -ka_meodmt * depot_meodmt
    d/dt(central_meodmt) <- ka_meodmt * depot_meodmt - elim_meodmt -
      q_meodmt * (Cc_meodmt - Cp_meodmt)
    d/dt(peripheral1_meodmt) <- q_meodmt * (Cc_meodmt - Cp_meodmt)
    f(depot_meodmt) <- fdepot_meodmt

    # ---- thermoregulation -------------------------------------------------
    rbase <- exp(lrbase)
    kout <- exp(lkout)
    kin <- kout * rbase # Eq. (3), the physiological steady state at time zero
    tau <- exp(ltau)
    emax_meodmt <- exp(lemax_meodmt)
    ec50_meodmt <- exp(lec50_meodmt_nonddi) * (1 - CONMED_HARMALINE) +
      exp(lec50_meodmt_ddi) * CONMED_HARMALINE
    ks_harmaline <- exp(lks_harmaline_single) * (1 - INJ_REPEAT) +
      exp(lks_harmaline_ddi) * INJ_REPEAT

    # Eq. (4). Carried as a state so that repeated injections superpose, which
    # is what Fig. 2 contrasts (single vs double saline injection). Dose one
    # unit into `stress` at each handling/injection event; the bioavailability
    # below scales it to the published magnitude S0.
    d/dt(stress) <- -exp(lkel_stress) * stress
    f(stress) <- s0_stress

    tempBasal <- rbase + drift_rbase * t # Eq. (1)
    sH <- ks_harmaline * Cc_harmaline # Eq. (5), 5-HT1A heat loss

    # Eqs. (6)-(8), the 5-HT2A transduction delay.
    d/dt(transit1) <- (emax_meodmt * Cc_meodmt / (ec50_meodmt + Cc_meodmt) -
      transit1) / tau
    d/dt(transit2) <- (transit1 - transit2) / tau
    d/dt(transit3) <- (transit2 - transit3) / tau

    # Eq. (9). The (tempBasal / temp) ratio is the adaptive negative feedback.
    d/dt(temp) <- kin * (tempBasal / temp) * (1 + stress) * (1 + transit3) -
      kout * (1 + sH) * temp
    temp(0) <- rbase

    temp ~ add(addSd_temp) + prop(propSd_temp)
  })
}
