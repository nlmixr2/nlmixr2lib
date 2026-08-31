Machavaram_2013_interleukin6_cyp3a_pbpk <- function() {
  description <- "PBPK (reduced from Simcyp Population-Based Simulator V11.0). Interleukin-6 (IL-6) disposition driving concentration- and time-dependent suppression of hepatic CYP3A4 activity, developed to predict disease-drug interactions in which elevated endogenous IL-6 suppresses CYP3A4 and a cytokine-modulating therapeutic protein reverses that suppression. This is the founding source of the enzyme-turnover equation d(Enzact)/dt = kdeg * Enz0 * [1 + (Emin - 1) * [IL-6] / (EC50 + [IL-6])] - kdeg * Enzact that later cytokine-CYP3A models cite. IL-6 is described as a one-compartment intravenous model (CLiv 1.0 L/h, Vss 0.43 L/kg, both estimated by fitting) driven by zero-order IL-6 infusions of 0.0009-1.5 ug/h that clamp systemic IL-6 to 1-1500 pg/mL, or by a finite-duration infusion representing a step change from baseline synthesis R0 (assumed zero) to a perturbed synthesis rate R1 for duration T followed by exponential decline. Hepatic CYP3A4 activity is carried relative to an untreated baseline of 1 with in vitro suppression constants EC50 73.2 pg/mL and Emin 0.24 (means over five human hepatocyte donors, Dickmann 2011) and a Simcyp library degradation rate constant kdeg 0.0193/h. Intestinal CYP3A4 is NOT a separate state: the paper assumed the magnitude of gut suppression equals the hepatic magnitude and implemented it by lowering intestinal abundance directly in the Simcyp population library rather than by a second turnover equation. The downstream simvastatin and cyclosporine exposure predictions are NOT part of this model: those used proprietary Simcyp V11.0 compound files whose in vivo clearances cannot be reconstructed from the published inputs, which report CLint without the MPPGL, liver weight, CYP3A4 abundance and hepatic blood flow needed to scale it."
  reference   <- "Machavaram KK, Almond LM, Rostami-Hodjegan A, Gardner I, Jamei M, Tay S, Wong S, Joshi A, Kenny JR. A physiologically based pharmacokinetic modeling approach to predict disease-drug interactions: suppression of CYP3A by IL-6. Clin Pharmacol Ther. 2013;94(2):260-268. doi:10.1038/clpt.2013.79. Enzyme-turnover equation from Eq. 1 (Methods) and Supplementary Methods Eq. 3; IL-6 disposition parameters and infusion-rate range from Methods 'Modeling of IL-6 profiles' and Supplementary Methods; IL-6 perturbation-input equations from Supplementary Methods Eq. 1 and Eq. 2; meta-analysis IL-6 population distributions from Results 'Meta-analysis and distribution of IL-6 in patients' and Figure 1; case-study designs from Methods and Supplementary Methods. In vitro suppression constants attributed to Dickmann LJ et al. Drug Metab Dispos. 2011;39:1415-1422; kdeg attributed to Rowland Yeo K et al. Eur J Pharm Sci. 2011;43:160-173."
  vignette    <- "Machavaram_2013_interleukin6_cyp3a_pbpk"
  units       <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Vss is reported in L/kg (Methods, 'Modeling of IL-6 profiles'), so the central volume is the per-kg value multiplied by body weight; this is a linear (exponent 1) weight scaling implied by the reported unit, not a fitted allometric exponent. Clearance is reported as an absolute 1.0 L/h and is NOT weight-scaled.",
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
    enzyme_3a4 = list(
      analyte  = "cytochrome P450 3A4 (hepatic pool)",
      units    = "fraction of untreated baseline (relative-to-baseline normalisation; the paper never prints an absolute hepatic CYP3A4 abundance, and every reported outcome is a percentage of baseline)",
      specimen = "tissue",
      verified = TRUE,
      notes    = "The bare isoenzyme form is used rather than the composed enzyme_3a4_liver because this model resolves CYP3A4 in a single tissue; the organ suffix is reserved for models that resolve the same isoenzyme across several organs (as in Chen_2024_interleukin6_cyp3a_pbpk, which carries both liver and gut pools). This paper assumed gut suppression equalled the hepatic magnitude and implemented it as a static abundance reduction in the Simcyp population library, so no gut state exists here."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 445,
    n_studies      = 4,
    age_range      = "18-83 years across the four case studies (28-72 rheumatoid arthritis; 18-38 and 20/34 bone marrow transplant; 48-83 postsurgical)",
    weight_range   = "not reported (Simcyp virtual populations generate weight from the built-in demographic distributions)",
    sex_female_pct = 45,
    disease_state  = "rheumatoid arthritis (case study 1), bone marrow transplant (case studies 2 and 3), postsurgical trauma (case study 4); IL-6 meta-analysis additionally covers healthy subjects",
    dose_range     = "IL-6 as zero-order intravenous infusion 0.0009-1.5 ug/h, clamping systemic IL-6 to 1, 5, 10, 50, 100, 200, 500, 1000 or 1500 pg/mL; simulation durations 10-19 days. Victim drugs (NOT modelled here) were simvastatin 40 mg orally and cyclosporine 1.5 mg/kg over 4 h or 1-3 mg/kg/day continuously",
    regions        = "not reported",
    notes          = "Simulations are Simcyp virtual populations, not fitted individuals. Case study 1: 10 trials x 12 rheumatoid arthritis subjects, aged 28-72 years, 66 percent female, 18-day simulation with simvastatin on day 15. Case study 2: 10 trials x 5 male bone marrow transplant subjects, aged 18-38 years, 11-day simulation with cyclosporine on day 10. Case study 3: 10 trials x 1 male subject matched to patients 1 (age 20) and 5 (age 34) of Chen 1994, 19-day continuous cyclosporine infusion; these two patients were selected because their observed IL-6 profiles could be fitted by zero-order input with first-order elimination. Case study 4: 10 trials x 16 postsurgical subjects, aged 48-83 years, 31 percent female. n_subjects is the sum over case studies (120 + 50 + 15 + 160 = 345 simulated subjects) plus the 100-subject notional meta-analysis contribution; sex_female_pct is the trial-size-weighted mean of the reported proportions. A literature meta-analysis gave weighted geometric mean systemic IL-6 of 4 pg/mL (95 percent CI 2-8) in healthy subjects, 54 pg/mL (43-69) in rheumatoid arthritis and 229 pg/mL (174-300) postsurgery, with lognormal population distributions parameterised as mu 3.99 / sigma 1.2 (rheumatoid arthritis) and mu 5.4 / sigma 1.0 (postsurgery) on the natural-log scale; those distributions are reproduced in the vignette and are not part of model()."
  )

  ini({
    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition.
    #
    # One-compartment intravenous model. Both values are described by the
    # paper as "estimated parameters by fitting", obtained by fitting the
    # individual IL-6 profiles of Chen 1994 with a zero-order input /
    # first-order elimination model. They are APPARENT parameters that
    # describe the slow rise and decline of the disease IL-6 profile, not the
    # true turnover of an IL-6 bolus (Vss 0.43 L/kg with CL 1.0 L/h implies a
    # ~21 h half-life at 70 kg, far longer than the minutes-to-hours half-life
    # of exogenous IL-6); see the vignette Errata.
    # -----------------------------------------------------------------------
    lcl <- fixed(log(1.0)); label("IL-6 clearance (L/h)")  # Methods, "Modeling of IL-6 profiles": CL i.v. = 1.0 L/h, estimated by fitting; same value in Supplementary Methods ("CLiv = 1.0 L/h")
    lvc <- fixed(log(0.43)); label("IL-6 volume of distribution at steady state (L/kg)")  # Methods, "Modeling of IL-6 profiles": Vss = 0.43 l/kg, estimated by fitting; same value in Supplementary Methods

    # -----------------------------------------------------------------------
    # Layer B -- hepatic CYP3A4 turnover and IL-6 suppression.
    #
    # All three constants are fixed inputs, not estimates: EC50 and Emin are
    # in vitro means over five human hepatocyte donors taken "directly from
    # the in vitro data reported by Dickmann et al.", and kdeg is a Simcyp
    # library value cited to Rowland Yeo 2011.
    # -----------------------------------------------------------------------
    emin <- fixed(0.24); label("Minimum CYP3A4 activity at maximal IL-6 suppression, as a fraction of vehicle control (unitless)")  # Methods "Modeling of enzyme dynamics" and Supplementary Methods: Emin 0.24. This is the activity FLOOR (suppressed TO 24 percent of control), not 24 percent suppression
    ec50 <- fixed(73.2); label("IL-6 concentration producing half-maximal CYP3A4 suppression (pg/mL)")  # Methods "Modeling of enzyme dynamics" and Supplementary Methods: EC50 73.2 pg/mL. Defined by the paper as "the concentration that supports half-Emin (i.e., half of the maximum suppressive effect)"
    kdeg <- fixed(0.0193); label("Hepatic CYP3A4 degradation rate constant (1/h)")  # Methods "Modeling of enzyme dynamics": "Mean degradation rate constant (kdeg) value in the liver used in the simulations was 0.0193/h", cited to Rowland Yeo 2011 (ref 38). Implies a turnover half-life of 35.9 h

    # -----------------------------------------------------------------------
    # Between-subject variability.
    #
    # The paper reports NO variance, %CV or confidence interval for any of the
    # five parameters above. IL-6 CL and Vss are single fitted point estimates;
    # Emin, EC50 and kdeg are single fixed inputs. Simcyp would have applied
    # its own population-library variability, but none of it is printed, so no
    # eta is declared rather than inventing one. The model is therefore
    # deterministic; the vignette generates its cohort by sampling the paper's
    # own published IL-6 population distributions instead. See vignette Errata.
    # -----------------------------------------------------------------------

    # -----------------------------------------------------------------------
    # Residual error. The source performs simulation only and reports no
    # residual error model for IL-6. Fixed to zero rather than invented.
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
    # 1.0 L/h and is not weight-scaled.
    # -----------------------------------------------------------------------
    cl <- exp(lcl)
    vc <- exp(lvc) * WT

    # -----------------------------------------------------------------------
    # Layer A -- IL-6 disposition. Dosed as zero-order intravenous infusions
    # of IL-6 amounts (mg) into `central`. Supplementary Methods Eq. 1 and
    # Eq. 2 describe the input as a step from baseline synthesis R0 (assumed
    # zero) to a perturbed synthesis rate R1 held for duration T, after which
    # concentrations decline exponentially:
    #
    #   during perturbation:  C(t) = R1/CL * (1 - exp(-CL/V * t))
    #   after perturbation:   C(t) = R1/CL * (1 - exp(-CL/V * T))
    #                                     * exp(-CL/V * (t - T))
    #
    # That is exactly a zero-order infusion of rate R1 and duration T into a
    # one-compartment model, so it needs no extra structure here: set T to the
    # study duration for the clamped-IL-6 case studies (1, 2 and 4) and to the
    # fitted perturbation duration for the individual profiles of case study 3.
    # The steady state reached during a maintained infusion is Css = R1/CL,
    # which is how the paper's 0.0009-1.5 ug/h rates map onto its 1-1500 pg/mL
    # target concentrations at CL = 1.0 L/h.
    # -----------------------------------------------------------------------
    d/dt(central) <- -cl / vc * central

    Cc <- central / vc * mgPerLToPgPerML

    # -----------------------------------------------------------------------
    # Layer B -- hepatic CYP3A4 turnover with IL-6 suppression. Eq. 1 of the
    # paper (Methods; restated as Eq. 3 of the Supplementary Methods):
    #
    #   d(Enzact,H-3A)/dt = kdeg,H-3A * Enz0,H-3A
    #                         * [1 + (Emin - 1) * [I]t / (EC50 + [I]t)]
    #                       - kdeg,H-3A * Enzact,H-3A
    #
    # IL-6 is assumed to act directly on the rate of enzyme SYNTHESIS; the
    # degradation rate constant is assumed unaffected by the perpetrator, and
    # in the absence of IL-6 the synthesis rate equals kdeg * Enz0 so the pool
    # sits at Enz0. The bracketed term is the fractional synthesis rate: it
    # equals 1 with no IL-6 present and falls toward Emin as IL-6 rises, so
    # the pool relaxes toward Enz0 * fsupp with a time constant of 1/kdeg.
    #
    # Enz0 is carried as 1 because the paper never prints an absolute hepatic
    # CYP3A4 abundance and reports every outcome as a percentage of baseline;
    # the state is therefore the fraction of untreated baseline activity and
    # Enz0 cancels out of the equation exactly.
    # -----------------------------------------------------------------------
    fsupp <- 1 + (emin - 1) * Cc / (ec50 + Cc)

    d/dt(enzyme_3a4) <- kdeg * fsupp - kdeg * enzyme_3a4

    enzyme_3a4(0) <- 1

    # -----------------------------------------------------------------------
    # Hepatic CYP3A4 activity as a percentage of the untreated baseline. This
    # is the quantity the paper plots in Figure 6 (case study 4) and the one
    # the vignette validates.
    #
    # The paper assumed that intestinal CYP3A4 suppression matches the hepatic
    # magnitude, but implemented it by lowering gut abundance directly in the
    # Simcyp population library rather than with a second turnover equation and
    # never reported an intestinal kdeg. No gut state is declared here: adding
    # one would require a kdeg this paper does not publish.
    # -----------------------------------------------------------------------
    cyp3a4Liver <- enzyme_3a4 * 100

    Cc ~ prop(propSd)
  })
}
