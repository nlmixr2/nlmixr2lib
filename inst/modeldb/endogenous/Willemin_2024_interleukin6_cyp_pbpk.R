Willemin_2024_interleukin6_cyp_pbpk <- function() {
  description <- "PBPK (reduced from Simcyp Simulator V21). Interleukin-6 (IL-6) disposition driving concentration- and time-dependent modulation of five hepatic cytochrome P450 activities (CYP1A2, 2C9, 2C19, 3A4, 3A5), developed to assess the drug-interaction risk created by the transient IL-6 elevation of cytokine release syndrome after teclistamab (BCMAxCD3 bispecific antibody) step-up and first treatment dosing in MajesTEC-1. Teclistamab itself is never modelled: because IL-6 is endogenous, its appearance is represented by a hand-calibrated series of six zero-order IV IL-6 infusions whose rates were adjusted to recover the observed MajesTEC-1 IL-6 profile, so the model is an IL-6 exposure driver rather than a teclistamab PK model. IL-6 is described as a one-compartment IV model (Vss 0.43 L/kg, CLiv 1 L/h) reduced from the Simcyp minimal-PBPK topology. Each CYP activity follows the Machavaram enzyme-turnover equation d(E)/dt = kdeg * (1 + (Indmax - 1) * [IL-6] / (IndC50 + [IL-6]) - E) with activity relative to an untreated baseline of 1, suppressing activity for CYP2C9, 2C19, 3A4 and 3A5 and inducing it for CYP1A2. The downstream victim-drug exposure ratios (caffeine, s-warfarin, omeprazole, midazolam, cyclosporine, simvastatin) are NOT part of this model: those used proprietary Simcyp V21 compound files whose in vivo dispositions cannot be reconstructed from the published inputs."
  reference   <- "Willemin ME, Wang Lin SX, De Zwart L, Wu LS, Miao X, Verona R, Banerjee A, Liu B, Kobos R, Qi M, Ouellet D, Goldberg JD, Girgis S. Evaluating drug interaction potential from cytokine release syndrome using a physiologically based pharmacokinetic model: A case study of teclistamab. CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1117-1129. doi:10.1002/psp4.13144. IL-6 disposition and interaction potencies from Table 1; IL-6 dosing regimens from supplement Table S3; steady-state enzyme activities from supplement Table S2; enzyme time-course targets from Table 3 and Figure 2. Enzyme-turnover equation form attributed by the paper to Machavaram KK et al. Clin Pharmacol Ther. 2013;94:260-268 and Machavaram KK et al. AAPS J. 2019;21:42; in vitro potencies to Dickmann LJ et al. Drug Metab Dispos. 2011;39:1415-1422 and Jiang X et al. AAPS J. 2016;18:767-776."
  vignette    <- "Willemin_2024_interleukin6_cyp_pbpk"
  units       <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Vss is reported in L/kg (Table 1), so the one-compartment volume is the per-kg value multiplied by body weight; this is the linear (exponent 1) scaling implied by the reported unit, not a fitted allometric exponent. CLiv is reported as an absolute 1 L/h and is NOT weight-scaled. The vignette uses 70 kg, the weight at which the reduction reproduces the published IL-6 peaks.",
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
    n_subjects     = 112,
    n_studies      = 1,
    age_range      = "20-50 years (Simcyp healthy-volunteer simulation population)",
    sex_female_pct = 50,
    disease_state  = "triple-class-exposed relapsed/refractory multiple myeloma with cytokine release syndrome (source of the observed IL-6 data); simulations themselves were run in a healthy-volunteer population",
    dose_range     = "teclistamab 0.06 and 0.3 mg/kg step-up doses followed by the 1.5 mg/kg subcutaneous first treatment dose (the IL-6 source regimen); IL-6 itself is dosed as six zero-order IV infusions of 0.0001-0.0114 mg",
    regions        = "MajesTEC-1 was a multinational phase I/II study (NCT03145181 / NCT04557098)",
    notes          = "Observed IL-6 concentration-time data come from up to 112 of the 119 patients (of 165 treated at the recommended phase II dose) who experienced cytokine release syndrome in MajesTEC-1 and whose IL-6 Cmax occurred before any tocilizumab administration, or who received no tocilizumab in cycle 1. Two IL-6 scenarios are modelled: scenario 1 is the mean IL-6 profile (Cmax 21 pg/mL) and scenario 2 the single patient with the highest observed IL-6 Cmax (288 pg/mL). Prospective simulations used 10 trials of 75 subjects aged 20-50 years, 50 percent female; the CYP-potency verification runs used 10 trials of 12 subjects at a clamped IL-6 of 50 pg/mL. Cycle 1 (the first 1.5 mg/kg treatment dose) begins 168 h after the first step-up dose, which is the time origin for the Table 3 enzyme-activity timings."
  )

  ini({
    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition.
    #
    # Table 1 specifies a Simcyp minimal-PBPK distribution model but reports
    # only Vss and CLiv; no compartmental volumes, tissue-partition
    # coefficients or inter-compartmental rates are given anywhere in the
    # paper or the supplement. The faithful reduction of what IS reported is
    # therefore a one-compartment IV model using Vss directly as the volume.
    # No Simcyp default was substituted for the unreported internals; see the
    # vignette Errata for the 6-7 percent consequence.
    # -----------------------------------------------------------------------
    lcl <- log(1); label("IL-6 clearance (L/h)")  # Table 1: CLiv 1 L/h (in vivo CL elimination model), attributed to Machavaram 2019; CV set to 50 percent
    lvc <- log(0.43); label("IL-6 volume of distribution at steady state (L/kg)")  # Table 1: Vss 0.43 L/kg (minimal PBPK distribution model), attributed to Machavaram 2019; CV set to 50 percent

    # -----------------------------------------------------------------------
    # Layer B -- IL-6 potency on hepatic CYP synthesis.
    #
    # Table 1 reports, per isoenzyme, a maximal fold change in the enzyme
    # synthesis rate (Simcyp "Indmax") and the IL-6 concentration producing
    # half of that change (Simcyp "IndC50", tabulated in both uM and pg/mL).
    # This file uses the Heathman-style amplitude parameterisation
    #
    #   fold = 1 + emax * C / (ec50 + C)
    #
    # so `emax` is Indmax - 1: POSITIVE for the CYP1A2 net induction and
    # NEGATIVE for the four suppressed isoenzymes. Each printed Indmax is
    # carried in the trailing comment so the paper value stays traceable.
    # ec50 is entered in pg/mL, the paper's own conversion of its uM value at
    # MW 21,000 g/mol, because every IL-6 concentration in the paper and in
    # this model is expressed in pg/mL.
    #
    # All ten constants are fixed inputs, not estimates: they are in vitro
    # potencies carried over from Dickmann 2011 and Jiang 2016.
    # -----------------------------------------------------------------------
    emax_1a2 <- fixed(0.34); label("Maximal fractional change in CYP1A2 synthesis rate (unitless; net INDUCTION)")  # Table 1: Indmax 1.34, so emax = 1.34 - 1 = 0.34; Jiang 2016, to account for the net induction observed with caffeine; CV = 30 percent
    ec50_1a2 <- fixed(8); label("IL-6 concentration producing half-maximal CYP1A2 induction (pg/mL)")  # Table 1: IndC50 3.81e-7 uM, which the table itself states "Corresponds to 8 pg/mL"; Jiang 2016; CV = 30 percent
    kdeg_1a2 <- fixed(0.0151); label("Hepatic CYP1A2 degradation rate constant (1/h)")  # NOT REPORTED by Willemin; back-solved from the Table 3 CYP1A2 activity extrema (118 and 128 percent) -- see vignette Errata. Implies a turnover half-life of 45.9 h

    emax_2c9 <- fixed(-0.947); label("Maximal fractional change in CYP2C9 synthesis rate (unitless; SUPPRESSION)")  # Table 1: Indmax 0.053, so emax = 0.053 - 1 = -0.947; Dickmann 2011; CV = 30 percent
    ec50_2c9 <- fixed(121); label("IL-6 concentration producing half-maximal CYP2C9 suppression (pg/mL)")  # Table 1: IndC50 5.76e-6 uM, which the table itself states "Corresponds to 121 pg/mL"; Dickmann 2011; CV = 30 percent
    kdeg_2c9 <- fixed(0.0059); label("Hepatic CYP2C9 degradation rate constant (1/h)")  # NOT REPORTED by Willemin; back-solved from the Table 3 CYP2C9 activity minima (95 and 77 percent) -- see vignette Errata. Implies a turnover half-life of 116.5 h

    emax_2c19 <- fixed(-0.786); label("Maximal fractional change in CYP2C19 synthesis rate (unitless; SUPPRESSION)")  # Table 1: Indmax 0.214, so emax = 0.214 - 1 = -0.786; Dickmann 2011; CV set to 50 percent to account for the variability observed by Jiang 2016
    ec50_2c19 <- fixed(71.3); label("IL-6 concentration producing half-maximal CYP2C19 suppression (pg/mL)")  # Table 1: IndC50 3.40e-6 uM, which the table itself states "Corresponds to 71.3 pg/mL"; Dickmann 2011; CV set to 50 percent
    kdeg_2c19 <- fixed(0.0233); label("Hepatic CYP2C19 degradation rate constant (1/h)")  # NOT REPORTED by Willemin; back-solved from the Table 3 CYP2C19 activity minima (87 and 56 percent) -- see vignette Errata. Implies a turnover half-life of 29.7 h

    emax_3a4 <- fixed(-0.76); label("Maximal fractional change in CYP3A4 synthesis rate (unitless; SUPPRESSION)")  # Table 1: Indmax 0.24, so emax = 0.24 - 1 = -0.76; Dickmann 2011; CV set to 50 percent to account for the variability observed by Jiang 2016
    ec50_3a4 <- fixed(73.2); label("IL-6 concentration producing half-maximal CYP3A4 suppression (pg/mL)")  # Table 1: IndC50 3.48e-6 uM, which the table itself states "Corresponds to 73.2 pg/mL"; Dickmann 2011; CV set to 50 percent
    kdeg_3a4 <- fixed(0.0153); label("Hepatic CYP3A4 degradation rate constant (1/h)")  # NOT REPORTED by Willemin; back-solved from the Table 3 CYP3A4 activity minima (89 and 63 percent) -- see vignette Errata. Implies a turnover half-life of 45.4 h

    # CYP3A5 shares CYP3A4's potencies by explicit assumption, and Table 3
    # reports identical activity minima and timings for the two isoenzymes, so
    # its back-solved kdeg is identical as well. Kept as separate parameters
    # and a separate compartment because the two enzymes drive different
    # victim drugs (midazolam/simvastatin vs cyclosporine) and a user may want
    # to break the assumption.
    emax_3a5 <- fixed(-0.76); label("Maximal fractional change in CYP3A5 synthesis rate (unitless; SUPPRESSION)")  # Table 1: Indmax 0.24, "Same as CYP3A4 (Machavaram et al. 2019)"; CV set to 50 percent
    ec50_3a5 <- fixed(73.2); label("IL-6 concentration producing half-maximal CYP3A5 suppression (pg/mL)")  # Table 1: IndC50 3.48e-6 uM = 73.2 pg/mL, "Same as CYP3A4 (Machavaram et al. 2019)"; CV set to 50 percent
    kdeg_3a5 <- fixed(0.0153); label("Hepatic CYP3A5 degradation rate constant (1/h)")  # NOT REPORTED by Willemin; back-solved from the Table 3 CYP3A5 activity minima (89 and 63 percent), which are identical to CYP3A4 -- see vignette Errata. Implies a turnover half-life of 45.4 h

    # -----------------------------------------------------------------------
    # Between-subject variability on IL-6 disposition.
    #
    # Table 1 reports Simcyp population %CV, converted with the log-normal
    # identity omega^2 = log(1 + CV^2). Table 1 also states CVs of 30-50
    # percent on the Indmax / IndC50 interaction inputs; those are NOT carried
    # as etas here (see vignette Errata) because the canonical amplitude
    # `emax` is negative for the four suppressed isoenzymes, so a log-normal
    # eta is not defined on it, and the paper does not report how Simcyp
    # correlates Indmax with IndC50. Every published quantity the vignette
    # reproduces is a typical-value or mean quantity.
    # -----------------------------------------------------------------------
    etalcl ~ 0.2231  # Table 1: CLiv CV set to 50 percent, so omega^2 = log(1 + 0.50^2) = 0.2231
    etalvc ~ 0.2231  # Table 1: Vss CV set to 50 percent, so omega^2 = log(1 + 0.50^2) = 0.2231

    # -----------------------------------------------------------------------
    # Residual error. The source is a PBPK simulation, not a fitted
    # population model, and reports no residual error model for IL-6. Fixed
    # to zero rather than invented; see the vignette Errata.
    # -----------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: the source performs simulation only
  })

  model({
    # -----------------------------------------------------------------------
    # Unit constant. Amounts are in mg and volumes in L, so central/vc is in
    # mg/L; 1 mg/L = 1e6 pg/mL, the unit in which IL-6 and every IndC50 in
    # this paper are reported.
    # -----------------------------------------------------------------------
    mgPerLToPgPerML <- 1e6

    # -----------------------------------------------------------------------
    # Individual parameters. Vss is reported in L/kg, so the volume scales
    # linearly with body weight; CLiv is an absolute 1 L/h and is not scaled.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * WT

    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition. Dosed as zero-order IV infusions of
    # hypothetical IL-6 amounts (mg) into `central`, per supplement Table S3.
    # -----------------------------------------------------------------------
    d/dt(central) <- -cl / vc * central

    Cc <- central / vc * mgPerLToPgPerML

    # -----------------------------------------------------------------------
    # Layer B -- hepatic CYP turnover under IL-6 modulation. Machavaram
    # enzyme-turnover form, with activity expressed relative to the untreated
    # baseline so that E(0) = 1:
    #
    #   dE/dt = kdeg * fold(C) - kdeg * E,  fold(C) = 1 + emax * C / (ec50 + C)
    #         = kdeg * (fold(C) - E)
    #
    # fold is the fractional synthesis rate: 1 with no IL-6 present, tending
    # to Indmax as IL-6 rises. Each pool therefore relaxes toward fold(C)
    # with a time constant of 1/kdeg, which is what makes the enzyme nadir
    # lag the IL-6 peak and makes a transient IL-6 pulse produce far less
    # suppression than a sustained one at the same concentration.
    #
    # The fold expression is written inline in each d/dt() rather than via a
    # shared named intermediate, so that no intermediate can be evaluated
    # outside the ODE context.
    #
    # IL-6 was given intravenously, so Simcyp propagated the modulation to
    # hepatic enzymes only; the paper handled the gut by editing intestinal
    # and colonic abundances offline between runs (Table S1). That offline
    # step is a static population edit, not a differential equation, and is
    # therefore not part of this model -- see the vignette Errata.
    # -----------------------------------------------------------------------
    d/dt(enzyme_1a2) <- kdeg_1a2 * (1 + emax_1a2 * Cc / (ec50_1a2 + Cc) - enzyme_1a2)
    d/dt(enzyme_2c9) <- kdeg_2c9 * (1 + emax_2c9 * Cc / (ec50_2c9 + Cc) - enzyme_2c9)
    d/dt(enzyme_2c19) <- kdeg_2c19 * (1 + emax_2c19 * Cc / (ec50_2c19 + Cc) - enzyme_2c19)
    d/dt(enzyme_3a4) <- kdeg_3a4 * (1 + emax_3a4 * Cc / (ec50_3a4 + Cc) - enzyme_3a4)
    d/dt(enzyme_3a5) <- kdeg_3a5 * (1 + emax_3a5 * Cc / (ec50_3a5 + Cc) - enzyme_3a5)

    enzyme_1a2(0) <- 1
    enzyme_2c9(0) <- 1
    enzyme_2c19(0) <- 1
    enzyme_3a4(0) <- 1
    enzyme_3a5(0) <- 1

    Cc ~ prop(propSd)
  })
}
