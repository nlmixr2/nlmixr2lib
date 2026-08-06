Chen_2024_interleukin6_cyp3a_pbpk <- function() {
  description <- "PBPK (reduced from Simcyp Simulator V17). Interleukin-6 (IL-6) disposition driving concentration- and time-dependent suppression of CYP3A enzyme activity in liver and gut, developed to assess the CYP3A drug-drug interaction risk created by the transient IL-6 elevation that follows mosunetuzumab (CD3/CD20 bispecific antibody) dosing. IL-6 is described as a one-compartment intravenous model (CL 2.8 L/h determined top-down; Vss 0.43 L/kg) reduced from the Simcyp minimal-PBPK-plus-single-adjusting-compartment topology, and is driven by hypothetical zero-order IL-6 infusions rather than by mosunetuzumab pharmacokinetics (mosunetuzumab itself is never modelled). CYP3A activity in each tissue follows the Machavaram enzyme-turnover equation d(ENZact)/dt = kdeg * ENZ0 * [1 + (Emin - 1) * [IL-6] / (EC50 + [IL-6])] - kdeg * ENZact with in vitro suppression constants EC50 43.7 pg/mL and Emin 0.217 (geometric means over five human hepatocyte donors, Dickmann 2011) and Simcyp library degradation rate constants kdeg 0.0193/h (liver) and 0.03/h (gut). The tissue-specific kdeg values are what separate the liver and gut responses to a transient IL-6 pulse. The downstream midazolam and simvastatin exposure-ratio predictions are NOT part of this model: those used unmodified proprietary Simcyp V17 compound files whose in vivo clearances cannot be reconstructed from the published inputs."
  reference   <- "Chen Y, Ma F, Jones N, Deng R, Li C, Li C-C. Assessment of CYP3A-mediated drug interaction via cytokine (IL-6) elevation for mosunetuzumab using physiologically-based pharmacokinetic modeling. CPT Pharmacometrics Syst Pharmacol. 2024;13(2):234-246. doi:10.1002/psp4.13073. Enzyme-turnover equation and in vitro suppression constants from Appendix S1 (Supplemental Material); equation attributed there to Machavaram KK et al. Clin Pharmacol Ther. 2013;94:260-268, in vitro constants to Dickmann LJ et al. Drug Metab Dispos. 2011;39:1415-1422."
  vignette    <- "Chen_2024_interleukin6_cyp3a_pbpk"
  units       <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Vss is reported in L/kg (Table 1), so the central volume is the per-kg value multiplied by body weight; this is a linear (exponent 1) weight scaling implied by the reported unit, not a fitted allometric exponent. Clearance is reported as an absolute 2.8 L/h and is NOT weight-scaled.",
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
    enzyme_liver = list(
      analyte  = "cytochrome P450 3A4",
      units    = "pmol/mg microsomal protein",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_gut = list(
      analyte  = "cytochrome P450 3A4",
      units    = "nmol/small intestine",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 212,
    n_studies      = 1,
    age_range      = "19-75 years (Simcyp healthy-volunteer simulation population)",
    sex_female_pct = 35,
    disease_state  = "relapsed/refractory non-Hodgkin lymphoma (source of the observed IL-6 data); CYP3A suppression verified against rheumatoid arthritis and post-hip-surgery cohorts",
    dose_range     = "mosunetuzumab 1/2/60 mg intravenous on cycle 1 days 1/8/15 (the IL-6 source regimen); IL-6 itself is dosed as hypothetical zero-order infusions of 0.0017-0.0968 mg",
    notes          = "Observed plasma IL-6 concentration-time profiles come from 212 patients in dose-expansion cohort B11 of study GO29781 (NCT02500407), measured by Quantikine ELISA over a validated range of 3.13-300 pg/mL (Appendix S1). Simulations used Simcyp virtual populations: healthy volunteers 10 trials x 10 subjects, age 19-75 years, 35 percent female for the mosunetuzumab application; rheumatoid arthritis 10 x 12, age 28-72 years, 67 percent female, and healthy volunteers 10 x 10, age 19-75 years, 60 percent female for the verification runs (Methods). The enzyme layer is validated in the vignette against Figure 4c (clamped IL-6) and Figure 4d (transient IL-6 after hip surgery)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition.
    #
    # Reduced to one compartment from the Simcyp minimal PBPK model with a
    # single adjusting compartment (SAC). Table 1 reports Vsac 0.05 L/kg but
    # also kin/kout = 0 (Default, not used), so the SAC is kinetically
    # disconnected and contributes no exchange; it only reduces the effective
    # systemic volume through the footnote relation Vsys = Vss - Vliver - Vsac,
    # and neither Vliver nor the hepatic blood flow is reported. Vss is
    # therefore used directly as the one-compartment volume.
    # -----------------------------------------------------------------------
    lcl <- log(2.8); label("IL-6 clearance (L/h)")  # Table 1 and Methods: CL determined top-down as 2.8 L/h; CLint 0.3 uL/min/10^6 cells back-calculated from it (Table 1 footnote a)
    lvc <- log(0.43); label("IL-6 volume of distribution at steady state (L/kg)")  # Table 1: Vss 0.43 L/kg (minimal PBPK distribution model)

    # -----------------------------------------------------------------------
    # Layer B -- CYP3A enzyme turnover and IL-6 suppression.
    #
    # All four constants are fixed inputs, not estimates: EC50 and Emin are
    # in vitro geometric means over five human hepatocyte donors, and both
    # kdeg values are Simcyp library values.
    # -----------------------------------------------------------------------
    emin <- fixed(0.217); label("Minimum CYP3A activity at maximal IL-6 suppression, as a fraction of vehicle control (unitless)")  # Appendix S1 donor table: geometric mean 21.7 percent, entered in Simcyp as Indmax 0.217; note this is the activity FLOOR (suppressed TO 21.7 percent), not 21.7 percent suppression
    ec50 <- fixed(43.7); label("IL-6 concentration producing half-maximal CYP3A suppression (pg/mL)")  # Appendix S1 donor table: geometric mean 43.7 pg/mL, equal to 2.08e-6 uM at MW 21000 (Table 1 footnote)
    kdeg_liver <- fixed(0.0193); label("Hepatic CYP3A degradation rate constant (1/h)")  # Table 1 and Appendix S1: Simcyp library value 0.0193/h
    kdeg_gut <- fixed(0.03); label("Intestinal CYP3A degradation rate constant (1/h)")  # Table 1: CYP3A4 turnover rate in gut 0.03/h (Simcyp library)

    # Baseline enzyme abundances. Healthy-volunteer values are used because
    # the paper reset the rheumatoid-arthritis abundances to the healthy
    # levels whenever CYP3A suppression was predicted from the in vitro data
    # (Table 1 footnote f), and Figure 4c reports the reduction from the
    # healthy-volunteer baseline. Both cancel out of the relative-activity
    # outputs cyp3aLiver and cyp3aGut.
    bl_enzyme_liver <- fixed(137); label("Baseline hepatic CYP3A4 abundance in healthy volunteers (pmol/mg microsomal protein)")  # Table 1: CYP3A4 in liver in HVs 137 pmol/mg-protein; the rheumatoid-arthritis population value is 82.2
    bl_enzyme_gut <- fixed(66.2); label("Baseline intestinal CYP3A4 abundance in healthy volunteers (nmol/small intestine)")  # Table 1: CYP3A4 in gut in HVs 66.2 nmol/small intestine; the rheumatoid-arthritis population value is 40.0

    # -----------------------------------------------------------------------
    # Between-subject variability.
    #
    # Table 1 reports these as Simcyp %CV, converted here with the log-normal
    # identity omega^2 = log(1 + CV^2). The values below are those used for
    # the mosunetuzumab application, where the CVs were deliberately inflated
    # above the Simcyp defaults so that the simulated 5th-95th percentile band
    # would cover the observed IL-6 data (Table 1 footnote b). The model
    # verification runs instead used the default 30 percent CV for both
    # (omega^2 = 0.0862); see the vignette.
    #
    # The CV on clearance is reported against CLint rather than CL; CL is a
    # monotone function of CLint, so the CV carries across.
    # -----------------------------------------------------------------------
    etalcl ~ 1.6094  # Table 1: CLint 200 %CV, so omega^2 = log(1 + 2.00^2) = 1.6094
    etalvc ~ 0.6931  # Table 1: Vss 100 %CV, so omega^2 = log(1 + 1.00^2) = 0.6931

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
    # mg/L; 1 mg/L = 1e6 pg/mL, the unit in which IL-6, EC50 and every
    # published IL-6 concentration in this paper are reported.
    # -----------------------------------------------------------------------
    mgPerLToPgPerML <- 1e6

    # -----------------------------------------------------------------------
    # Individual parameters. Vss is reported in L/kg, so the central volume
    # scales linearly with body weight; clearance is reported as an absolute
    # 2.8 L/h and is not weight-scaled.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * WT

    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition. Dosed as zero-order intravenous infusions
    # of hypothetical IL-6 amounts (mg) into `central`.
    # -----------------------------------------------------------------------
    d/dt(central) <- -cl / vc * central

    Cc <- central / vc * mgPerLToPgPerML

    # -----------------------------------------------------------------------
    # Layer B -- CYP3A enzyme turnover with IL-6 suppression, Appendix S1
    # Equation 1 (attributed to Machavaram 2013), applied to liver and gut:
    #
    #   d(ENZact)/dt = kdeg * ENZ0 * [1 + (Emin - 1) * [IL-6] / (EC50 + [IL-6])]
    #                  - kdeg * ENZact
    #
    # The bracketed term is the fractional synthesis rate: it equals 1 with no
    # IL-6 present and falls toward Emin as IL-6 rises, so each enzyme pool
    # relaxes toward ENZ0 * fsupp with a time constant of 1/kdeg. The gut
    # kdeg is faster than the liver kdeg, which is why a transient IL-6 pulse
    # suppresses gut CYP3A more deeply than hepatic CYP3A.
    #
    # Suppression is driven by systemic IL-6 in both tissues: the paper states
    # that the degree of gut CYP3A suppression is assumed similar to that in
    # the liver when the gut effect is considered (Methods, Appendix S1).
    # -----------------------------------------------------------------------
    fsupp <- 1 + (emin - 1) * Cc / (ec50 + Cc)

    d/dt(enzyme_liver) <- kdeg_liver * bl_enzyme_liver * fsupp - kdeg_liver * enzyme_liver
    d/dt(enzyme_gut) <- kdeg_gut * bl_enzyme_gut * fsupp - kdeg_gut * enzyme_gut

    enzyme_liver(0) <- bl_enzyme_liver
    enzyme_gut(0) <- bl_enzyme_gut

    # -----------------------------------------------------------------------
    # Relative CYP3A activity, as a fraction of the untreated baseline. These
    # are the quantities the paper plots (Figures 3b, 4c, 4d) and the ones the
    # vignette validates; the absolute abundances cancel out of them.
    # -----------------------------------------------------------------------
    cyp3aLiver <- enzyme_liver / bl_enzyme_liver
    cyp3aGut <- enzyme_gut / bl_enzyme_gut

    Cc ~ prop(propSd)
  })
}
