Munana_2018_levetiracetam_dog <- function() {
  description <- "Veterinary (dog). One-compartment population PK model with first-order absorption and first-order elimination for extended-release levetiracetam (LEV-XR) in client-owned dogs with idiopathic epilepsy, sampled at steady state on their established q12h maintenance regimen. Parameterised in rate-constant form (k01, k10, V/F) per kg body weight, as published. Absorption is slower than elimination, so the terminal slope is absorption-limited (flip-flop) and the apparent 5 h terminal half-life is the absorption half-life. Concomitant phenobarbital raises the apparent volume of distribution V/F about 2.6-fold, lowering Cmax and AUC; concomitant zonisamide has no effect (Munana 2018)"
  reference <- paste(
    "Munana KR, Otamendi AJ, Nettifee JA, Papich MG. Population pharmacokinetics",
    "of extended-release levetiracetam in epileptic dogs when administered alone,",
    "with phenobarbital or zonisamide. J Vet Intern Med. 2018;32(5):1677-1683.",
    "doi:10.1111/jvim.15298.",
    "Structural and random-effect forms taken from the Supporting Information",
    "(JVIM-32-1677-s001.pdf, Equations 2-4); parameter values from Table 1.",
    sep = " "
  )
  vignette <- "Munana_2018_levetiracetam_dog"

  # Every published parameter is per kg body weight (V/F in L/kg, CL/F in
  # L/kg/hr), and the paper normalised every observed concentration to the mean
  # study dose of 29.4 mg/kg (Results 'AED administration'). Doses are therefore
  # given in mg/kg and compartment amounts carried in mg/kg, so central / vc is
  # mg/L, which is the ug/mL of the HPLC assay and of Table 1.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(
      analyte = "levetiracetam",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "levetiracetam",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONMED_PB = list(
      description = "Indicator for concomitant phenobarbital maintenance therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extended-release levetiracetam alone, the L group)",
      notes = paste(
        "Time-fixed cohort indicator. CONMED_PB = 1 identifies the six dogs of the LP group,",
        "on stable phenobarbital plus LEV-XR (mean phenobarbital dose 2.86 mg/kg q12h, SD",
        "1.40; mean trough serum phenobarbital 28.48 ug/mL, SD 11.67; Results 'AED",
        "administration'). Primidone is not mentioned by the source and no dog received it,",
        "so the pooling question raised in the register entry does not arise here. Every dog",
        "was required to be at steady state on all of its AEDs, with no dose change for at",
        "least five half-lives, so the indicator does not vary within the observation window.",
        "The source models it as one level of a three-level 'treatment group' categorical",
        "covariate (L / LP / LZ) entered multiplicatively on the log scale of V/F",
        "(Supporting Information Equation 4); this column and CONMED_ZNS together carry that",
        "three-level factor, with the L group as the shared reference. Phenobarbital is a",
        "broad-spectrum CYP and UGT inducer and the source attributes the effect to increased",
        "presystemic metabolism of levetiracetam reducing the absorbed fraction F, which",
        "inflates both V/F and CL/F (Discussion).",
        sep = " "
      ),
      source_name = "treatment group (LP)"
    ),
    CONMED_ZNS = list(
      description = "Indicator for concomitant zonisamide maintenance therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extended-release levetiracetam alone, the L group)",
      notes = paste(
        "Time-fixed cohort indicator. CONMED_ZNS = 1 identifies the LZ group, on stable",
        "zonisamide plus LEV-XR (mean zonisamide dose 7.82 mg/kg q12h, SD 2.22; mean trough",
        "serum zonisamide 55.09 ug/mL, SD 33.63; Results 'AED administration'). Six LZ dogs",
        "were enrolled and five analysed: one was excluded before the population fit for",
        "receiving more than 200 mg/kg of LEV-XR (Methods 'Population pharmacokinetics').",
        "The second level of the source's three-level 'treatment group' categorical covariate;",
        "see the CONMED_PB notes. Its retained effect on V/F is null at the precision Table 1",
        "prints - the LEV-XR + zonisamide and LEV-XR alone columns both give theta V = 0.15",
        "L/kg - which is the paper's headline negative finding: 'coadministration of",
        "zonisamide was not shown to contribute to the variability' (Conclusions).",
        sep = " "
      ),
      source_name = "treatment group (LZ)"
    )
  )

  # Body weight, age and sex were screened as covariates on k01, k10 and V and
  # none reduced the -2LL enough to be retained (Results 'Pharmacokinetic
  # analysis using population model'). The source reports these screens only
  # graphically (eta box plots for sex, eta-versus-covariate scatter plots for
  # weight and age) and publishes no point estimate for any of them, so they are
  # documented rather than encoded. Body weight is nonetheless present in the
  # model implicitly: every structural parameter is published per kg, which is
  # exact linear weight scaling of V/F and CL/F. The paper's negative finding is
  # that weight explains no variability BEYOND that per-kg normalisation.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a continuous covariate on k01, k10 and V via eta-versus-covariate",
        "scatter plots; not retained. Median 25.7 kg, range 7.8-45.5 kg (Results 'Dog",
        "demographics'). Discussion: 'Body weight also was explored as a possible covariate",
        "while considering that dog size can affect drug absorption or clearance, but it was",
        "not shown to be significant in this model.' Note that the published parameters are",
        "already per kg, so this is a statement about residual weight dependence only.",
        sep = " "
      ),
      source_name = "body weight (kg)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a continuous covariate on k01, k10 and V; not retained. Range 3-12",
        "years, median 6 years (Results 'Dog demographics').",
        sep = " "
      ),
      source_name = "age (years)"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neutered male)",
      notes = paste(
        "Screened as a categorical covariate on k01, k10 and V via eta box plots; not",
        "retained. Every dog was neutered or spayed: 11 spayed females and 6 neutered males",
        "among the 17 analysed (Results 'Dog demographics'), so the source's two-level",
        "'gender' covariate maps onto SEXF without an intact/neutered stratum.",
        sep = " "
      ),
      source_name = "gender (neutered male or spayed female)"
    )
  )

  population <- list(
    species = "dog (client-owned; mixed breed n = 6, Labrador retriever n = 3, Australian shepherd n = 2, and one each of Basset hound, Golden retriever, Pembroke Welsh corgi, Vizsla, Curly-coated retriever, English springer spaniel)",
    n_subjects = 17L,
    n_enrolled = 18L,
    n_studies = 1L,
    n_observations = 85L,
    age_range = "3-12 years",
    age_median = "6 years",
    weight_range = "7.8-45.5 kg",
    weight_median = "25.7 kg",
    sex_female_pct = 64.7,
    disease_state = "idiopathic epilepsy, at least International Veterinary Epilepsy Task Force tier 1 confidence, median duration 1 year",
    dose_range = "500, 750, 1000 or 1500 mg LEV-XR PO q12h per dog (group means 31.86, 30.91 and 23.52 mg/kg in the L, LP and LZ groups; study mean 29.4 mg/kg)",
    co_medication = "6 dogs LEV-XR alone (L), 6 LEV-XR + phenobarbital (LP), 6 enrolled / 5 analysed LEV-XR + zonisamide (LZ); no drugs other than monthly parasite preventatives",
    regions = "United States (NC State Veterinary Hospital and one regional veterinary hospital)",
    notes = paste(
      "Demographics from Results 'Dog demographics' and 'AED administration'. 18 dogs were",
      "enrolled, 6 per group; one LZ dog given more than 200 mg/kg of LEV-XR was excluded",
      "before the population fit, leaving 17 (Methods 'Population pharmacokinetics'), and the",
      "reported demographics (11 spayed females + 6 neutered males) are for those 17. Five",
      "samples per dog at 0, 2, 4, 8 and 12 h after the morning dose, giving 85 observations;",
      "the 0 h sample is the pre-dose steady-state trough. Dogs were fasted overnight and fed",
      "at the time of the morning dose. All dogs received generic LEV-XR; four manufacturers",
      "were represented (Apotex in over half the dogs, Qualitest, Lupin, BluePoint) and",
      "manufacturer was investigated as a covariate but the unequal distribution across four",
      "levels precluded a valid analysis, so it was neither retained nor reported (Methods",
      "'Population pharmacokinetics'; Discussion). LEV dose per dog (500 / 750 / 1000 / 1500",
      "mg, in 8 / 7 / 1 / 1 dogs) was screened as a four-level categorical covariate and not",
      "retained: the eta box plots hinted at an effect on V for the 1000 and 1500 mg levels,",
      "but each contained a single dog (Results 'Pharmacokinetic analysis using population",
      "model'). Concentrations were normalised to the study mean dose of 29.4 mg/kg before",
      "fitting, which is why the model's dose unit is mg/kg.",
      sep = " "
    )
  )

  ini({
    # ---- Structural parameters (Table 1, 'Overall (all groups)' column) ------
    # The source parameterises the one-compartment oral model by rate constants
    # and an apparent volume (Methods Equation 1: C = (F*D*k01) / (V*(k01-k10))
    # * (exp(-k10*t) - exp(-k01*t))), and places an independent exponential eta
    # on each of k01, k10 and V (Supporting Information Equation 2). That
    # parameterisation is kept rather than being recast as cl + vc, because the
    # equivalent clearance form would need etalcl = etalkel + etalvc, a
    # correlated eta block the authors did not fit.

    lka <- log(0.138)
    label("Absorption rate constant k01 (1/h)")
    # Table 1, row 'theta k 01', 'Overall (all groups)' Value = 0.138 1/hr.
    # Cross-checks against the same column's 'k 01 half-life' row:
    # log(2) / 0.138 = 5.02 h versus the printed 5.01 h.

    lkel <- log(0.505)
    label("Elimination rate constant k10 (1/h)")
    # Table 1, row 'theta k 10', 'Overall (all groups)' Value = 0.505 1/hr.
    # Cross-check: log(2) / 0.505 = 1.372 h versus the printed 'k 10 half-life'
    # of 1.37 h. Because k01 < k10 the profile is flip-flop, so this is the
    # FASTER of the two half-lives and is not the terminal slope (Discussion).

    lvc <- log(0.151)
    label("Apparent volume of distribution V/F without concomitant phenobarbital (L/kg)")
    # Table 1, row 'theta V', 'Overall (all groups)' Value = 0.151 L/kg. This is
    # the reference level of the treatment covariate: Table 1 prints theta V =
    # 0.15 L/kg for BOTH the 'LEV-XR alone' and the 'LEV-XR + Zonisamide'
    # columns, which is 0.151 at the two significant figures those columns
    # carry, so the reference covers the L and LZ groups alike. Table 1's
    # 'Overall' secondary parameters are all reproduced by this triple of
    # primary values: CL/F = k10 * V = 0.505 * 0.151 = 0.076 L/kg/hr (printed
    # 0.08), AUC = 29.4 / 0.076 = 386 h*ug/mL (printed 388.72), Tmax =
    # log(k01/k10) / (k01-k10) = 3.53 h (printed 3.53) and Cmax = 32.7 ug/mL
    # (printed 32.99). It is V/F, not V: no intravenous dose was given, so any
    # change in F appears here (Results; Discussion).

    # ---- Covariate effects on V/F (Supporting Information Equation 4) --------
    # Equation 4 is V_i = theta_V * exp(dvd_treatment) * exp(eta_i,V), i.e. the
    # three-level treatment factor enters multiplicatively on the log scale with
    # one coefficient per level. Treatment group was the ONLY covariate retained,
    # and it acts on V alone.

    e_conmed_pb_vc <- log(0.39 / 0.151)
    label("Effect of concomitant phenobarbital on log apparent volume of distribution")
    # Table 1, row 'theta V', 'LEV-XR + phenobarbital' Value = 0.39 L/kg against
    # the 0.151 L/kg reference, so the coefficient is log(0.39 / 0.151) = 0.949
    # and V/F is 2.58-fold higher on phenobarbital. Corroborated independently by
    # two of the same table's secondary parameters, which under a V-only
    # covariate must scale as 1/V: the Cmax ratio is 33.01 / 13.38 = 2.47 and the
    # AUC ratio is 352.95 / 134.86 = 2.62, bracketing 2.58.

    e_conmed_zns_vc <- 0
    label("Effect of concomitant zonisamide on log apparent volume of distribution")
    # Table 1, row 'theta V', 'LEV-XR + Zonisamide' Value = 0.15 L/kg, the same
    # as the 'LEV-XR alone' column, so the zonisamide contrast is
    # log(0.15 / 0.15) = 0 at the printed precision. Retained at zero rather than
    # dropped because the source's covariate is a single three-level factor and
    # this is one of its levels; the null is the paper's stated finding, not an
    # omission.

    # ---- Inter-individual variability ---------------------------------------
    # Exponential IIV on every primary parameter (Supporting Information
    # Equation 2: P_i = P_pop * exp(eta_i,P), eta ~ N(0, omega^2)). Entries below
    # are the omega^2 values printed in Table 1's 'Overall (all groups)' Omega^2
    # sub-column, which is a log-scale variance: each reproduces its neighbouring
    # CV% as sqrt(exp(omega^2) - 1), e.g. sqrt(exp(0.221) - 1) = 49.7% against
    # the printed 49.66%.
    etalka ~ 0.058 # Table 1, row 'theta k 01', 'Overall (all groups)' Omega 2 = 0.058 (CV 24.39%)
    etalkel ~ 0.001 # Table 1, row 'theta k 10', 'Overall (all groups)' Omega 2 = 0.001 (CV 3.32%)
    etalvc ~ 0.221 # Table 1, row 'theta V', 'Overall (all groups)' Omega 2 = 0.221 (CV 49.66%); the largest IIV in the model and the one the treatment covariate was introduced to explain

    # ---- Residual error ------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual SD (fraction); magnitude not reported by the source")
    # Supporting Information Equation 3 specifies a multiplicative residual error
    # model, C_obs,ij = C_pred,ij * (1 + eps_ij) with eps ~ N(0, sigma^2), chosen
    # over additive, log-additive, power and mixed alternatives. The form is
    # published but sigma itself is not: Table 1 carries Omega^2 and CV% columns
    # for the random effects only, and neither the main text nor the Supporting
    # Information prints a residual variance. Set to zero so simulations are not
    # given an invented residual magnitude; supply a value via ini() before using
    # this model for anything that needs realistic scatter.
  })

  model({
    # ---- Individual parameters ----------------------------------------------
    ka <- exp(lka + etalka)
    kel <- exp(lkel + etalkel)
    # Equation 4: the treatment factor is multiplicative on V, i.e. additive on
    # log V, with the LEV-XR-alone group as reference.
    vc <- exp(lvc + e_conmed_pb_vc * CONMED_PB + e_conmed_zns_vc * CONMED_ZNS + etalvc)

    # ---- ODE system ----------------------------------------------------------
    # One compartment with first-order absorption and first-order elimination;
    # the analytic solution of this system is the source's Equation 1. F is not
    # identifiable without an intravenous reference dose, so no f(depot) term is
    # applied: the whole dose enters `depot` and the fitted volume is V/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- Observation ---------------------------------------------------------
    # central in mg/kg over vc in L/kg gives mg/L, i.e. the ug/mL of the HPLC
    # assay and of Table 1.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
