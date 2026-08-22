Willemin_2024_interleukin6_cyp_talquetamab <- function() {
  description <- "Reduced from a Simcyp Simulator V21 minimal-PBPK analysis. Interleukin-6 (IL-6) disposition driving concentration- and time-dependent modulation of five hepatic cytochrome P450 activities (CYP1A2, 2C9, 2C19, 3A4, 3A5), developed to assess the drug-interaction risk created by the transient IL-6 elevation of cytokine release syndrome after talquetamab (GPRC5D x CD3 bispecific antibody) step-up and first treatment dosing in MonumenTAL-1, at both recommended phase 2 doses (0.4 mg/kg weekly and 0.8 mg/kg every other week). Talquetamab itself is never modelled: because IL-6 is endogenous, its appearance is represented by a series of zero-order IV IL-6 infusions whose rates the authors adjusted to recover the observed MonumenTAL-1 IL-6 profile, so the model is an IL-6 exposure driver rather than a talquetamab PK model. IL-6 is described here as a two-compartment IV model recovered from the paper's own simulated IL-6 profiles (Figs 1 and 2); a one-compartment reduction cannot reproduce them. Each CYP activity follows the enzyme-turnover equation d(E)/dt = kdeg * (1 + emax * C / (ec50 + C) - E) with activity relative to an untreated baseline of 1, suppressing activity for CYP2C9, 2C19, 3A4 and 3A5 and inducing it for CYP1A2. The downstream victim-drug exposure ratios (caffeine, S-warfarin, omeprazole, midazolam, cyclosporine, simvastatin) are NOT part of this model: those used proprietary Simcyp V21 compound files whose in vivo dispositions cannot be reconstructed from the published inputs."
  reference   <- "Willemin ME, Gong J, Hilder BW, Masterson T, Tolbert J, Renaud T, Heuck C, Kane C, De Zwart L, Girgis S, Ma X, Ouellet D. Evaluation of drug-drug interaction potential of talquetamab, a T-cell-redirecting GPRC5D x CD3 bispecific antibody, as a result of cytokine release syndrome in patients with relapsed/refractory multiple myeloma in MonumenTAL-1, using a physiologically based pharmacokinetic model. Target Oncol. 2024;19(6):965-975. doi:10.1007/s11523-024-01093-6. IL-6 disposition recovered from the simulated profiles of Figs 1 and 2; hepatic CYP turnover rate constants recovered from the activity time courses of Figs 3 and 4 and gated against Table 3. The IL-6 model itself is stated by this paper to be the previously published one of Willemin ME et al. CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1117-1129 (doi:10.1002/psp4.13144), from which the interaction potencies (Indmax, IndC50) are carried; those in turn are attributed to Dickmann LJ et al. Drug Metab Dispos. 2011;39:1415-1422 and Jiang X et al. AAPS J. 2016;18:767-776, and the turnover equation form to Machavaram KK et al. Clin Pharmacol Ther. 2013;94:260-268."
  vignette    <- "Willemin_2024_interleukin6_cyp_talquetamab"
  units       <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The IL-6 volumes are expressed per kilogram, so both compartment volumes scale linearly with body weight; this is the scaling implied by the L/kg unit of the source IL-6 model, not a fitted allometric exponent. Clearances are absolute and are NOT weight-scaled. The vignette uses 70 kg, the weight at which the reduction reproduces the published IL-6 peaks.",
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "interleukin-6",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "interleukin-6",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    enzyme_1a2 = list(
      analyte  = "cytochrome P450 1A2",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_2c9 = list(
      analyte  = "cytochrome P450 2C9",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_2c19 = list(
      analyte  = "cytochrome P450 2C19",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_3a4 = list(
      analyte  = "cytochrome P450 3A4",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_3a5 = list(
      analyte  = "cytochrome P450 3A5",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 195,
    n_studies      = 1,
    age_range      = "20-50 years (Simcyp healthy-volunteer simulation population)",
    sex_female_pct = 50,
    disease_state  = "relapsed/refractory multiple myeloma with cytokine release syndrome (source of the observed IL-6 data); the simulations themselves were run in a healthy-volunteer population",
    dose_range     = "talquetamab 0.01 and 0.06 mg/kg step-up doses then 0.4 mg/kg subcutaneous weekly; or 0.01, 0.06 and 0.3 mg/kg step-up doses then 0.8 mg/kg subcutaneous every other week (the IL-6 source regimens). IL-6 itself is dosed as a series of zero-order IV infusions.",
    regions        = "MonumenTAL-1 was a multinational phase I/II study (NCT03399799 / NCT04634552)",
    notes          = "Observed IL-6 concentration-time data come from 100 patients in the 0.4 mg/kg weekly cohort and 95 patients in the 0.8 mg/kg every-other-week cohort who experienced cytokine release syndrome and who either received no tocilizumab in cycle 1 or whose IL-6 Cmax occurred before tocilizumab was given. Two IL-6 scenarios are modelled per dosing schedule: scenario 1 is the median IL-6 profile (Cmax 18.4 pg/mL weekly, 7.07 pg/mL every other week) and scenario 2 is the single patient with the highest observed IL-6 Cmax (213 and 3503 pg/mL respectively). Prospective simulations used five trials of 100 subjects aged 20-50 years, 50 percent female. Cycle 1 (the first full treatment dose) begins 96 h after the first step-up dose on the weekly schedule and 168 h after it on the every-other-week schedule; those are the time origins for the Table 3 enzyme-activity timings."
  )

  ini({
    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition.
    #
    # This paper reports no IL-6 disposition parameters at all; it states only
    # that "the PBPK model of IL-6 used in this analysis has been previously
    # published" (Sect. 2.2, citing the teclistamab analysis) and that
    # distribution was modelled by a minimal PBPK model with non-mechanistic
    # elimination. The three rate constants below were therefore recovered from
    # this paper's OWN published output: the simulated IL-6 concentration-time
    # profiles of Figs 1 and 2 (four scenarios), fitted jointly so that a single
    # disposition serves all four, with the peaks anchored to the Cmax values
    # stated in the Table 1 footnotes. A one-compartment reduction was tested
    # first and rejected: it misses the Fig 2b decay by up to 4.5-fold (see the
    # vignette Errata).
    #
    #   kel = cl/vc = 0.0411 /h    k12 = q/vc = 0.0236 /h    k21 = q/vp = 0.0437 /h
    #
    # The figures identify those three rate constants and the ratio
    # Vss/Vc = 1.540, but not the absolute volume scale: the IL-6 infusion
    # amounts are themselves unpublished, so any volume can be absorbed into
    # them. The scale is therefore set by the one disposition number the cited
    # upstream IL-6 model does publish, Vss = 0.43 L/kg, which splits by the
    # recovered ratio into the vc and vp below.
    #
    # Fixed rather than estimated: they are recovered constants of a published
    # deterministic simulation, not estimates from a fit to observed data.
    # -----------------------------------------------------------------------
    lcl <- fixed(log(0.803)); label("IL-6 clearance (L/h)")  # recovered from Figs 1-2 as kel = 0.0411 /h, times vc
    lvc <- fixed(log(0.279)); label("IL-6 central volume of distribution (L/kg)")  # Vss 0.43 L/kg divided by the Vss/Vc = 1.540 recovered from Figs 1-2; scaled by body weight
    lq  <- fixed(log(0.461)); label("IL-6 intercompartmental clearance (L/h)")  # recovered from Figs 1-2 as k12 = 0.0236 /h, times vc
    lvp <- fixed(log(0.151)); label("IL-6 peripheral volume of distribution (L/kg)")  # Vss 0.43 L/kg minus vc; gives k21 = q/vp = 0.0437 /h

    # -----------------------------------------------------------------------
    # Layer B -- IL-6 potency on hepatic CYP synthesis.
    #
    # Amplitude parameterisation
    #
    #   fold = 1 + emax * C / (ec50 + C)
    #
    # so `emax` is the Simcyp Indmax minus 1: POSITIVE for the CYP1A2 net
    # induction and NEGATIVE for the four suppressed isoenzymes, and `ec50` is
    # the Simcyp IndC50 in pg/mL.
    #
    # These ten constants are NOT tabulated in this paper. It states that the
    # IL-6 model it uses is the previously published one and that, "other than
    # CYP1A2, for which an induction effect was modeled per clinical data, a
    # suppression effect was modeled for all studied CYP enzymes", citing the
    # same in vitro sources. The values are therefore carried from that model
    # (Table 1 of the teclistamab analysis, doi:10.1002/psp4.13144). They are
    # independently confirmed against THIS paper's own output: with the IL-6
    # baseline of about 2.6 pg/mL that Fig 1a settles at, these potencies give
    # terminal activities of 108.3, 98.0, 97.3 and 97.4 percent for CYP1A2,
    # 2C9, 2C19 and 3A4, matching the 109, 97.5, 97.5 and 97.5 percent plateaus
    # of Fig 3a. See the vignette source-trace table.
    # -----------------------------------------------------------------------
    emax_1a2 <- fixed(0.34); label("Maximal fractional change in CYP1A2 synthesis rate (unitless; net INDUCTION)")  # Indmax 1.34, so emax = 1.34 - 1
    ec50_1a2 <- fixed(8); label("IL-6 concentration producing half-maximal CYP1A2 induction (pg/mL)")  # IndC50 3.81e-7 uM = 8 pg/mL at MW 21,000 g/mol
    kdeg_1a2 <- fixed(0.0197); label("Hepatic CYP1A2 degradation rate constant (1/h)")  # recovered from the Fig 3a/3b/4a CYP1A2 activity time courses; turnover half-life 35.2 h

    emax_2c9 <- fixed(-0.947); label("Maximal fractional change in CYP2C9 synthesis rate (unitless; SUPPRESSION)")  # Indmax 0.053, so emax = 0.053 - 1
    ec50_2c9 <- fixed(121); label("IL-6 concentration producing half-maximal CYP2C9 suppression (pg/mL)")  # IndC50 5.76e-6 uM = 121 pg/mL
    kdeg_2c9 <- fixed(0.00583); label("Hepatic CYP2C9 degradation rate constant (1/h)")  # recovered from the Fig 3a/3b/4a CYP2C9 activity time courses; turnover half-life 119.0 h

    emax_2c19 <- fixed(-0.786); label("Maximal fractional change in CYP2C19 synthesis rate (unitless; SUPPRESSION)")  # Indmax 0.214, so emax = 0.214 - 1
    ec50_2c19 <- fixed(71.3); label("IL-6 concentration producing half-maximal CYP2C19 suppression (pg/mL)")  # IndC50 3.40e-6 uM = 71.3 pg/mL
    kdeg_2c19 <- fixed(0.0240); label("Hepatic CYP2C19 degradation rate constant (1/h)")  # recovered from the Fig 3a/3b/4a CYP2C19 activity time courses; turnover half-life 28.9 h

    emax_3a4 <- fixed(-0.76); label("Maximal fractional change in CYP3A4 synthesis rate (unitless; SUPPRESSION)")  # Indmax 0.24, so emax = 0.24 - 1
    ec50_3a4 <- fixed(73.2); label("IL-6 concentration producing half-maximal CYP3A4 suppression (pg/mL)")  # IndC50 3.48e-6 uM = 73.2 pg/mL
    kdeg_3a4 <- fixed(0.0161); label("Hepatic CYP3A4 degradation rate constant (1/h)")  # recovered from the Fig 3a/3b/4a CYP3A4/CYP3A5 activity time courses; turnover half-life 43.0 h

    # CYP3A5 shares CYP3A4's potencies by the source IL-6 model's explicit
    # assumption, and Table 3 of this paper reports activity minima and timings
    # for the two that are identical apart from a single percentage point in one
    # of four scenarios, so the recovered kdeg is the same. Kept as separate
    # parameters and a separate compartment because the two enzymes drive
    # different victim drugs (midazolam/simvastatin vs cyclosporine) and a user
    # may want to break the assumption.
    emax_3a5 <- fixed(-0.76); label("Maximal fractional change in CYP3A5 synthesis rate (unitless; SUPPRESSION)")  # Indmax 0.24, same as CYP3A4
    ec50_3a5 <- fixed(73.2); label("IL-6 concentration producing half-maximal CYP3A5 suppression (pg/mL)")  # IndC50 3.48e-6 uM = 73.2 pg/mL, same as CYP3A4
    kdeg_3a5 <- fixed(0.0161); label("Hepatic CYP3A5 degradation rate constant (1/h)")  # recovered together with CYP3A4; turnover half-life 43.0 h

    # -----------------------------------------------------------------------
    # Between-subject variability on IL-6 disposition.
    #
    # This paper prints no variance estimate. It does state that the IL-6 model
    # it uses is the previously published one, whose Table 1 sets a population
    # CV of 50 percent on both IL-6 clearance and volume; those are carried
    # here, converted with the log-normal identity omega^2 = log(1 + CV^2). The
    # magnitude is corroborated by this paper's own Fig 1a, whose 90 percent
    # confidence band at the median peak spans roughly 9.4 to 31.5 pg/mL about
    # a median of 18.4, implying a CV of about 38 percent on peak concentration.
    # No variability is carried on the interaction potencies: the canonical
    # amplitude `emax` is negative for the four suppressed isoenzymes, so a
    # log-normal eta is not defined on it. See the vignette Errata.
    # -----------------------------------------------------------------------
    etalcl ~ 0.2231  # CV 50 percent on IL-6 clearance, so omega^2 = log(1 + 0.50^2) = 0.2231
    etalvc ~ 0.2231  # CV 50 percent on IL-6 volume, so omega^2 = log(1 + 0.50^2) = 0.2231

    # -----------------------------------------------------------------------
    # Residual error. The source is a PBPK simulation, not a fitted population
    # model, and reports no residual error model for IL-6. Fixed to zero rather
    # than invented; see the vignette Errata.
    # -----------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: the source performs simulation only
  })

  model({
    # -----------------------------------------------------------------------
    # Unit constant. Amounts are in mg and volumes in L, so central/vc is in
    # mg/L; 1 mg/L = 1e6 pg/mL, the unit in which IL-6 and every IndC50 in this
    # paper are reported.
    # -----------------------------------------------------------------------
    mgPerLToPgPerML <- 1e6

    # -----------------------------------------------------------------------
    # Individual parameters. Both volumes are per kilogram and so scale
    # linearly with body weight; both clearances are absolute.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq)
    vp <- exp(lvp) * WT

    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition. Dosed as zero-order IV infusions of
    # hypothetical IL-6 amounts (mg) into `central`; the vignette gives the
    # infusion schedules recovered for each of the four published scenarios.
    # -----------------------------------------------------------------------
    d/dt(central)     <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    Cc <- central / vc * mgPerLToPgPerML

    # -----------------------------------------------------------------------
    # Layer B -- hepatic CYP turnover under IL-6 modulation, with activity
    # expressed relative to the untreated baseline so that E(0) = 1:
    #
    #   dE/dt = kdeg * fold(C) - kdeg * E,  fold(C) = 1 + emax * C / (ec50 + C)
    #         = kdeg * (fold(C) - E)
    #
    # fold is the fractional synthesis rate: 1 with no IL-6 present, tending to
    # Indmax as IL-6 rises. Each pool therefore relaxes toward fold(C) with a
    # time constant of 1/kdeg, which is what makes the enzyme nadir lag the
    # IL-6 peak and makes a transient IL-6 pulse produce far less suppression
    # than a sustained one at the same concentration.
    #
    # The fold expression is written inline in each d/dt() rather than via a
    # shared named intermediate, so that no intermediate can be evaluated
    # outside the ODE context.
    #
    # IL-6 was given intravenously, so Simcyp propagated the modulation to
    # hepatic enzymes only; the paper handled the gut by editing intestinal and
    # colonic abundances offline between two runs. That offline step is a
    # static population edit, not a differential equation, and is therefore not
    # part of this model -- see the vignette Errata.
    # -----------------------------------------------------------------------
    d/dt(enzyme_1a2)  <- kdeg_1a2 * (1 + emax_1a2 * Cc / (ec50_1a2 + Cc) - enzyme_1a2)
    d/dt(enzyme_2c9)  <- kdeg_2c9 * (1 + emax_2c9 * Cc / (ec50_2c9 + Cc) - enzyme_2c9)
    d/dt(enzyme_2c19) <- kdeg_2c19 * (1 + emax_2c19 * Cc / (ec50_2c19 + Cc) - enzyme_2c19)
    d/dt(enzyme_3a4)  <- kdeg_3a4 * (1 + emax_3a4 * Cc / (ec50_3a4 + Cc) - enzyme_3a4)
    d/dt(enzyme_3a5)  <- kdeg_3a5 * (1 + emax_3a5 * Cc / (ec50_3a5 + Cc) - enzyme_3a5)

    enzyme_1a2(0)  <- 1
    enzyme_2c9(0)  <- 1
    enzyme_2c19(0) <- 1
    enzyme_3a4(0)  <- 1
    enzyme_3a5(0)  <- 1

    Cc ~ prop(propSd)
  })
}
