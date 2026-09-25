Geerts_2018_amyloidNeurotransmission_qsp <- function() {
  description <- paste(
    "QSP. Amyloid-beta modulation of glutamatergic and nicotinic",
    "neurotransmission in Alzheimer's disease and mild cognitive impairment.",
    "Four ODE states carry the two low-order A-beta aggregate pools (A-beta40",
    "and A-beta42 monomer + dimer + trimer, on the paper's own arbitrary",
    "0-16 'unit' load scale) and the progressive loss of cortical neuron and",
    "synapse density. A-beta deposition is linear (zero-order) at 1 unit per",
    "13 weeks for the APOE4+/- heterozygote, scaled +/-50% per epsilon-4",
    "allele; neuron and synapse density decline linearly at 0.35%/week and",
    "0.04%/week respectively and are assumed independent of A-beta. The two",
    "A-beta loads drive two algebraic neurophysiological readouts: the",
    "relative maximum excitatory-excitatory NMDA conductance (gNmdaRel,",
    "Equations 1A/1B -- BIPHASIC in A-beta40, rising linearly to 1 + delta at",
    "x0 then falling with slope alpha, and monotonically falling in A-beta42",
    "with slope alpha*), and the relative alpha-7 nicotinic acetylcholine",
    "receptor activation (a7Rel, Equation 2, falling linearly in the summed",
    "load). Amyloid-lowering treatment enters as a fractional reduction of",
    "each deposition rate (fracRed40 / fracRed42; 0 = placebo), which is",
    "exactly how the paper introduces the solanezumab, verubecestat (BACE",
    "inhibitor) and semagacestat (gamma-secretase inhibitor) target-engagement",
    "data. Disease state is a scenario switch (SW_MCI): mild-to-moderate AD",
    "by default (30% cholinergic deficit), MCI when set to 1 (30% cholinergic",
    "INCREASE and a 3% lower neuron and synapse density).",
    "SCOPE LIMIT -- this file encodes the paper's own calibrated A-beta",
    "coupling and pathology-progression layer ONLY. The downstream transfer",
    "function from (NMDA conductance, alpha-7 activation) to an ADAS-Cog",
    "score is the proprietary In Silico Biosciences biophysical cortical",
    "network (80 pyramidal cells + 40 GABAergic interneurons simulated in",
    "NEURON 7.2), which is described in the upstream reference [16] and is",
    "not reproducible from this paper; the model therefore outputs the",
    "neurophysiological coupling variables, NOT ADAS-Cog. Deterministic:",
    "no IIV, no residual error and no dosing events are reported.",
    sep = " "
  )
  reference <- paste(
    "Geerts H, Spiros A, Roberts P (2018).",
    "Impact of amyloid-beta changes on cognitive outcomes in Alzheimer's",
    "disease: analysis of clinical trials using a quantitative systems",
    "pharmacology model.",
    "Alzheimer's Research & Therapy 10:14.",
    "doi:10.1186/s13195-018-0343-5.",
    "The cortical network platform that converts the coupling variables",
    "encoded here into an ADAS-Cog readout is Roberts PD, Spiros A, Geerts H",
    "(2012) Alzheimers Res Ther 4(6):50 (reference [16] of the 2018 paper)",
    "and is NOT encoded in this file.",
    sep = " "
  )
  vignette <- "Geerts_2018_amyloidNeurotransmission_qsp"

  # Paper-mechanistic states. None maps onto a canonical nlmixr2lib
  # compartment role: abeta40 / abeta42 are cortical low-order aggregate
  # loads on the paper's arbitrary 0-16 discretized scale, and neuron /
  # synapse are relative cortical densities rather than drug amounts.
  paper_specific_compartments <- c("abeta40", "abeta42", "neuron", "synapse")

  units <- list(
    time = "week",
    dosing = "none",
    concentration = paste(
      "A-beta loads are in the paper's own arbitrary 'units' (Methods,",
      "'Alzheimer pathology and amyloid deposition': 1 unit is the amount of",
      "A-beta accumulated in 13 weeks by an APOE4+/- heterozygote), discretized",
      "0-16. neuron and synapse are dimensionless fractions of the",
      "trial-baseline density. gNmdaRel and a7Rel are dimensionless ratios to",
      "the A-beta-free reference (g0 and 'activation0' of Equations 1A/1B",
      "and 2).",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Geerts 2018 Methods ('Alzheimer
  # pathology and amyloid deposition') and the Equation 1A/1B definitions of
  # x and y.
  compartmentData <- list(
    abeta40 = list(
      analyte = "amyloid-beta 1-40 low-order aggregates (monomer, dimer, trimer)",
      units = "arbitrary A-beta load units (0-16 scale)",
      specimen = "tissue",
      verified = TRUE
    ),
    abeta42 = list(
      analyte = "amyloid-beta 1-42 low-order aggregates (monomer, dimer, trimer)",
      units = "arbitrary A-beta load units (0-16 scale)",
      specimen = "tissue",
      verified = TRUE
    ),
    neuron = list(
      analyte = "cortical neurons",
      units = "fraction of trial-baseline density",
      specimen = "tissue",
      verified = TRUE
    ),
    synapse = list(
      analyte = "cortical synapses",
      units = "fraction of trial-baseline density",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    APOE4_COUNT = list(
      description = paste(
        "APOE-epsilon4 allele count (0 = APOE4-/-, 1 = APOE4+/-,",
        "2 = APOE4+/+). Drives two distinct effects in Geerts 2018 Methods,",
        "'APOE genotype': the A-beta deposition rate (1.50, 1.00 and 0.50",
        "units/13 weeks for APOE4+/+, +/- and -/- respectively) and the",
        "baseline cortical SYNAPSE density (-20% for APOE4+/+ and +20% for",
        "APOE4-/- relative to the heterozygote). Neuron density is NOT",
        "affected by APOE in this paper.",
        sep = " "
      ),
      units = "(count, 0 / 1 / 2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CENTRING VALUE = 1 allele, i.e. the APOE4+/- HETEROZYGOTE, not the",
        "cohort mean and not the non-carrier. Geerts 2018 states both effects",
        "explicitly relative to the heterozygote ('a 50% difference compared",
        "to the deposition rates of the heterozygote APOE4+/- form at 1",
        "unit/13 weeks'; 'as compared to the wild-type APOE4+/- genotype').",
        "Both published effects are exactly linear in allele count over the",
        "three genotypes (deposition 0.5 / 1.0 / 1.5; synapse density",
        "1.20 / 1.00 / 0.80), which is why the additive APOE4_COUNT canonical",
        "is used here rather than the non-additive APOE4_HET + APOE4_HOM",
        "pair. Time-invariant (germline genotype).",
        sep = " "
      ),
      source_name = "APOE4"
    )
  )

  covariatesDataExcluded <- list()

  population <- list(
    species = "human (in silico virtual patients)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Two simulated populations. (1) Mild-to-moderate Alzheimer's disease,",
      "MMSE 18-24, followed for 78 weeks -- the calibration population of the",
      "upstream ADAS-Cog platform. (2) Minimal / mild cognitive impairment",
      "(MCI), simulated as a 3% decrease in synapse and neuron density plus a",
      "30% INCREASE in cholinergic tone.",
      sep = " "
    ),
    dose_range = paste(
      "No dosing events. Amyloid-lowering interventions enter as the",
      "fractional reduction in A-beta deposition rate reported as clinical",
      "target engagement (Results, 'Therapeutic amyloid-beta interventions'):",
      "verubecestat / BACE inhibitor 80-90% (A-beta40) and 60-80%",
      "(A-beta42); semagacestat / gamma-secretase inhibitor 30-50% and",
      "15-30%; solanezumab 5-10% and 30-50%.",
      sep = " "
    ),
    notes = paste(
      "No individual-level cohort: this is a deterministic typical-value QSP",
      "model. The upstream cortical network was calibrated against 28",
      "retrospective historical treatment outcomes (placebo from the",
      "flurbiprofen / tarenflurbil trials at 72 weeks; donepezil and",
      "rivastigmine at 2 doses x 3 time points; galantamine at 3 doses x 3",
      "time points; SB742457 at 2 doses x 2 time points) in mild-to-moderate",
      "AD, giving R-squared above 0.6 (Results / Methods 'Calibration of the",
      "network', citing reference [16]). The coupling parameters encoded here",
      "were then constrained against three independent clinical datasets:",
      "baseline ADAS-Cog in A-beta+ vs A-beta- MCI subjects (Doraiswamy 2014",
      "Mol Psychiatry, reference [17] Table 1: 10.8 vs 8.5, with normal",
      "elderly controls at 5.6 and 4.1); the scopolamine dose-response in MCI",
      "(Lim 2015 Neurobiol Aging, reference [18]); and the APOE effect on",
      "placebo cognitive trajectory in AD (Samtani 2015, reference [19]).",
      "Absolute baseline ADAS-Cog is 20-22 in the AD population depending on",
      "APOE genotype, and 4.5 for healthy cognitively normal controls in the",
      "MCI model. No datasets were generated or analysed by the paper itself",
      "('Availability of data and materials').",
      sep = " "
    )
  )

  ini({
    # =================================================================
    # Every value is from the MAIN TEXT of Geerts 2018 (the paper has no
    # supplement -- EuropePMC hasSuppl = N, and 'Availability of data and
    # materials' states no datasets were generated or analysed).
    #
    # These are mechanistic constants of a deterministic platform,
    # reported as point values with no standard errors, no IIV and no
    # residual error, so they are kept on the LINEAR scale under names
    # derived from the paper's own symbols (see parameter-names.md
    # 'Endogenous / mechanistic parameters') and every one is wrapped in
    # fixed(): the paper reports a single calibrated parameterisation and
    # estimates nothing from a user's data.
    # =================================================================

    # ----- A-beta / neurotransmission coupling (Equations 1A, 1B, 2) -----
    # Final calibrated set, Results 'Sensitivity analysis' closing sentence:
    # 'We use the values x0 = 2, delta = 0.025, beta = 0.03, alpha = 0.002
    # and alpha* = 0.002 for the analysis in the next section.'
    # NOTE this supersedes the earlier exploratory set quoted in Results
    # 'Constraining system parameters using clinical data' (delta = 0.015,
    # alpha = 0.0015, alpha* = 0.00035, beta = 0.025) and the intermediate
    # set of Figure 3 (delta = 0.025, alpha = 0.002, alpha* = 0.002).
    x0_ab40 <- fixed(2)
    label("A-beta40 load at the maximal positive NMDA effect, x0 (A-beta load units)")
    delta_g <- fixed(0.025)
    label("Maximal relative increase in NMDA conductance from A-beta40, delta (unitless)")
    alpha_ab40 <- fixed(0.002)
    label("Slope of NMDA conductance decline above x0 per A-beta40 unit, alpha (1/unit)")
    alpha_ab42 <- fixed(0.002)
    label("Slope of NMDA conductance decline per A-beta42 unit, alpha* (1/unit)")
    beta_a7 <- fixed(0.03)
    label("Coupling factor of summed A-beta load on alpha-7 nAChR activation, beta (1/unit)")

    # Results 'Constraining system parameters using clinical data': the
    # A-beta+/A-beta- partition of the 17x17 matrix is 'x > 2 and y > 2'
    # versus 'x <= 2 and y <= 2' (Figure 2 legend), and the corresponding
    # 'cutoff value between A-beta+ and A-beta- load' is quoted as 3 units.
    # On the integer matrix the two statements are the same partition, so a
    # >= 3 test is used below. 3 units corresponds to a florbetapir SUVR of
    # 1.34 (Results, citing reference [33]); that anchor is NOT implemented
    # as a units->SUVR conversion because the paper gives only the single
    # anchor plus cohort-average SUVRs (1.0 for A-beta-, 1.5 for A-beta+),
    # which do not determine the assumed linear slopes.
    abetaCutoff <- fixed(3)
    label("A-beta+ / A-beta- classification cutoff on each load axis (A-beta load units)")

    # ----- A-beta deposition -----
    # Methods, 'Alzheimer pathology and amyloid deposition': 'amyloid
    # deposition as a function of trial duration is arbitrarily set at 1
    # unit/13 weeks for both A-beta40 and A-beta42 ... this value assumes
    # linear growth'. Rate is for the APOE4+/- heterozygote.
    #
    # INTERNAL INCONSISTENCY IN THE PAPER, adjudicated here. The Discussion
    # restates the genotype-specific rates as '1.5 units/13 weeks, 1
    # unit/week and 0.5 units/13 weeks for the APOE4+/+, APOE4+/- and
    # APOE4-/- genotype' -- i.e. 'per week', not 'per 13 weeks', for the
    # heterozygote only. The Methods value (1 unit/13 weeks) is the correct
    # one and is used here: it is stated twice in Methods, the two flanking
    # Discussion values both use /13 weeks, and it is the only value
    # consistent with the paper's own arithmetic in Results -- both the
    # worked GSI example (a 40% reduction giving 'an increase of 0.6 units
    # along the A-beta40 axis ... /13 weeks') and the published 3x3
    # weighting factors (which imply mass-average loads of 4.21 and 4.39
    # after one 13-week step under an 80%/60% reduction from a load of 4).
    # A 'per week' rate would make both of those wrong by 13-fold.
    kdep40 <- fixed(1 / 13)
    label("A-beta40 deposition rate in the APOE4+/- heterozygote (load units/week)")
    kdep42 <- fixed(1 / 13)
    label("A-beta42 deposition rate in the APOE4+/- heterozygote (load units/week)")
    # Methods, 'APOE genotype': 1.50 units/13 weeks for APOE4+/+ and 0.50
    # for APOE4-/-, 'a 50% difference compared to the deposition rates of
    # the heterozygote APOE4+/- form at 1 unit/13 weeks'.
    e_apoe4_kdep <- fixed(0.50)
    label("Fractional change in A-beta deposition rate per APOE4 allele above the heterozygote (unitless)")

    # ----- Progressive neurodegeneration -----
    # Methods, 'Alzheimer pathology and amyloid deposition': 'Progressive
    # neurodegeneration is simulated as a linear loss of neurons (at
    # 0.35%/week) and synapses (0.04%/week), values that were constrained
    # from historical clinical trials. We assume this neuronal degeneration
    # is independent from A-beta deposition.'
    klossNeuron <- fixed(0.0035)
    label("Linear cortical neuron loss rate, fraction of baseline density per week (1/week)")
    klossSynapse <- fixed(0.0004)
    label("Linear cortical synapse loss rate, fraction of baseline density per week (1/week)")
    # Methods, 'APOE genotype': synaptic density '-20% for APOE4+/+ and
    # +20% for APOE4-/- relative to the APOE4+/- heterozygote genotype'.
    # Applies to SYNAPSES ONLY -- neuron density carries no APOE effect.
    e_apoe4_syn <- fixed(0.20)
    label("Fractional change in baseline synapse density per APOE4 allele above the heterozygote (unitless)")

    # ----- Disease state -----
    # Methods, 'Alzheimer pathology and amyloid deposition': 'AD pathology
    # is introduced as a cholinergic deficit of 30% ... except for the case
    # of MCI where a compensatory increase of 30% is used'. Results: 'We
    # simulate an MCI patient population in the cognitive model using a 3%
    # decrease in synapse and neuron density in addition to a 30% increase
    # in cholinergic tone'.
    SW_MCI <- fixed(0)
    label("Scenario switch: 0 = mild-to-moderate AD (default), 1 = MCI population")
    densLossMci <- fixed(0.03)
    label("Baseline neuron and synapse density decrease in the MCI population (unitless fraction)")
    cholDefAd <- fixed(0.30)
    label("Cholinergic tone deficit in the mild-to-moderate AD population (unitless fraction)")
    cholIncMci <- fixed(0.30)
    label("Compensatory cholinergic tone increase in the MCI population (unitless fraction)")

    # ----- Baseline A-beta load (initial conditions) -----
    # Default is the paper's own 'mild baseline A-beta+ load' of 4 units
    # (Results / Figure 4a: "patients with a 'mild baseline' A-beta+ load
    # (4 units which is above the cutoff level of 3)"). Figure 6b uses a
    # high load of 8 units and Figure 6a a low load below the 3-unit
    # cutoff; override these to reproduce those arms.
    abeta40Bl <- fixed(4)
    label("Baseline A-beta40 load at trial start (A-beta load units)")
    abeta42Bl <- fixed(4)
    label("Baseline A-beta42 load at trial start (A-beta load units)")

    # ----- Amyloid-lowering treatment -----
    # Results, 'Therapeutic amyloid-beta interventions on anticipated
    # clinical outcome': the reported clinical target engagement is applied
    # as a proportional reduction of the DEPOSITION RATE, not of the
    # standing load -- 'Patients on the active drugs have a proportionally
    # lower A-beta40 and A-beta42 increase according to their biomarker
    # change; for instance, in the case of a low dose of GSI, a 40%
    # reduction in A-beta40 and a 20% reduction in A-beta42 corresponds to
    # an increase of 0.6 units along the A-beta40 axis and 0.80 units along
    # the A-beta42 axis/13 weeks.'
    # Reported ranges (defaults are placebo = 0):
    #   verubecestat (BACE-I)  fracRed40 0.80-0.90, fracRed42 0.60-0.80
    #   semagacestat (GSI)     fracRed40 0.30-0.50, fracRed42 0.15-0.30
    #   solanezumab            fracRed40 0.05-0.10, fracRed42 0.30-0.50
    fracRed40 <- fixed(0)
    label("Fractional reduction of the A-beta40 deposition rate by treatment, 0 = placebo (unitless)")
    fracRed42 <- fixed(0)
    label("Fractional reduction of the A-beta42 deposition rate by treatment, 0 = placebo (unitless)")
  })

  model({
    # ---- 1. Derived covariate and scenario terms ----
    # APOE4 effects are centred on the HETEROZYGOTE (APOE4_COUNT = 1), which
    # is the genotype the paper states both effects relative to.
    # apoeDep: 0.50 / 1.00 / 1.50 for 0 / 1 / 2 epsilon-4 alleles.
    apoeDep <- 1 + e_apoe4_kdep * (APOE4_COUNT - 1)
    # apoeSyn: 1.20 / 1.00 / 0.80 for 0 / 1 / 2 epsilon-4 alleles.
    apoeSyn <- 1 - e_apoe4_syn * (APOE4_COUNT - 1)

    # Disease-state baselines. MCI lowers BOTH densities by 3%; APOE acts on
    # synapses only. The AD arm starts at the trial baseline (1), since the
    # paper's neurodegeneration rates are quoted per week of trial duration.
    neuron0 <- 1 - densLossMci * SW_MCI
    synapse0 <- apoeSyn * (1 - densLossMci * SW_MCI)

    # Cholinergic tone relative to a cognitively normal control: 0.70 in the
    # AD arm (30% deficit), 1.30 in the MCI arm (30% compensatory increase).
    # Reported for completeness -- it is an input to the upstream cortical
    # network, not to Equations 1A/1B or 2.
    cholTone <- (1 - cholDefAd) * (1 - SW_MCI) + (1 + cholIncMci) * SW_MCI

    # ---- 2. ODE system ----
    # Linear (zero-order) A-beta accumulation, reduced proportionally by
    # treatment target engagement.
    d/dt(abeta40) <- kdep40 * apoeDep * (1 - fracRed40)
    d/dt(abeta42) <- kdep42 * apoeDep * (1 - fracRed42)
    abeta40(0) <- abeta40Bl
    abeta42(0) <- abeta42Bl

    # Linear neurodegeneration, assumed independent of A-beta deposition.
    # 'Linear loss' means a CONSTANT absolute rate, referenced to that
    # subject's own baseline density -- not first-order decay.
    d/dt(neuron) <- -klossNeuron * neuron0
    d/dt(synapse) <- -klossSynapse * synapse0
    neuron(0) <- neuron0
    synapse(0) <- synapse0

    # ---- 3. Neurophysiological coupling readouts ----
    # Equation (1A), x <= x0:
    #   g(x, y) = g0 * [1 + delta * (x / x0) - y * alpha*]
    # Equation (1B), x > x0:
    #   g(x, y) = g0 * [1 + delta + (x0 - x) * alpha - y * alpha*]
    # Continuous at x = x0, where both give g0 * [1 + delta - y * alpha*].
    # Expressed as the ratio g / g0.
    gNmdaRel <- 1 + delta_g * (abeta40 / x0_ab40) - alpha_ab42 * abeta42
    if (abeta40 > x0_ab40) {
      gNmdaRel <- 1 + delta_g + alpha_ab40 * (x0_ab40 - abeta40) -
        alpha_ab42 * abeta42
    }

    # Equation (2): alpha7(x, y) activation
    #   = alpha7 activation0 * [1 - beta * (x + y)]
    # 'here the effects for the two A-beta forms are identical', so the
    # summed load drives it. Expressed as the ratio to activation0.
    a7Rel <- 1 - beta_a7 * (abeta40 + abeta42)

    # ---- 4. A-beta+ / A-beta- classification ----
    # Figure 2 legend: 'A-beta- subjects with x <= 2 and y <= 2, A-beta+
    # subjects with x > 2 and y > 2'; equivalently a cutoff of 3 units on
    # the integer load matrix.
    abetaPos <- (abeta40 >= abetaCutoff) * (abeta42 >= abetaCutoff)
  })
}
