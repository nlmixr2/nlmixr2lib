Xu_2013_ranibizumab <- function() {
  description <- "One-compartment population PK model for total ranibizumab in serum following intravitreal injection in adults with neovascular age-related macular degeneration (Xu et al. 2013, IOVS). Vitreous humor acts as a slow-release depot with first-order absorption Ka into a systemic central compartment with first-order clearance CL/F; a small parallel fraction of each intravitreal dose reaches the central compartment via a rapid needle-track shunt bypassing the vitreous. Covariates: Cockcroft-Gault creatinine clearance on CL/F (power) and concomitant verteporfin PDT on Ka (multiplicative). Serum ranibizumab was measured; vitreous concentration is computed algebraically from vitreous amount over an assumed 4 mL vitreous humor volume for downstream simulation and is not observed. Data pooled from two Phase 1, two Phase 1/2, and one Phase 3 trial (MARINA / FOCUS)."
  reference   <- "Xu L, Lu T, Tuomi L, Jumbe N, Lu J, Eppler S, Kuebler P, Damico-Beyer LA, Joshi A. Pharmacokinetics of ranibizumab in patients with neovascular age-related macular degeneration: a population approach. Invest Ophthalmol Vis Sci. 2013;54(3):1616-1624. doi:10.1167/iovs.12-10260."
  vignette    <- "Xu_2013_ranibizumab"
  units       <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault). Xu 2013 median = 65.22 mL/min in the analysis population (adults with neovascular AMD).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline value). Reference value for the power model is 65.22 mL/min (Xu 2013 median). Cockcroft-Gault used (Xu 2013 Fig 1 legend); values NOT BSA-normalized, mirroring the raw-Cockcroft-Gault convention in Delattre_2010_amikacin.R.",
      source_name        = "CrCL"
    ),
    CONMED_VERTEPORFIN = list(
      description        = "Concomitant verteporfin photodynamic therapy indicator. 1 = subject received one or more concomitant PDT procedures during the ranibizumab treatment window; 0 = no concomitant PDT. 42% of the analysis population received concomitant PDT (Xu 2013 Results).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant verteporfin PDT)",
      notes              = "Time-fixed baseline covariate in Xu 2013 (no time-varying covariates in the analysis). See CONMED_VERTEPORFIN entry in inst/references/covariate-columns.md for provenance.",
      source_name        = "PDT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline (years). Xu 2013 median = 78 years.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate on CL/F but not retained in the final model: 'Age was not statistically significantly correlated with CL/F after correction of CrCL, as the age effect was already reflected in CrCL' (Xu 2013 Results)."
    ),
    CNV_TYPE = list(
      description = "Choroidal neovascularization subtype: predominantly classic / minimally classic / occult.",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened per Fig 1 but not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 229L,
    n_studies      = 5L,
    age_median     = "78 years",
    weight_median  = "76 kg",
    sex_female_pct = 56,
    race_ethnicity = c(Caucasian = 95),
    disease_state  = "Neovascular (wet) age-related macular degeneration (AMD).",
    dose_range     = "Single or repeated intravitreal (ITV) doses of 0.05-2.0 mg/eye. Phase 3 (MARINA) monthly 0.3 or 0.5 mg/eye for up to 24 months; FOCUS 0.5 mg/eye monthly with concomitant verteporfin PDT (administered 7 days before treatment); Phase 1/1-2 single-dose or intra-subject escalation q2w / q4w.",
    regions        = "Multi-center clinical trials sponsored by Genentech.",
    notes          = "Analysis included 696 evaluable ranibizumab serum concentration records from 229 subjects with at least one measurable serum concentration. Pooled from Rosenfeld 2006 MARINA (Phase 3), Heier 2006 FOCUS (Phase 1/2), Heier 2006 (Phase 1/2), Rosenfeld 2005 (Phase 1 intra-subject escalation), Rosenfeld 2005 (Phase 1 single-dose). See Xu 2013 Table 1 for study details and Table 2 for sample counts. 42% of subjects received concomitant PDT; 22% received local or systemic IOP-lowering medication. Median CrCL = 65.22 mL/min; 16% of CrCL values were imputed to the median."
  )

  ini({
    # Structural parameters (Xu 2013 Table 3). h1..h4 are the population
    # typical values. CL/F is the apparent systemic clearance, Vc/F is the
    # apparent central volume of distribution, Ka is the first-order rate of
    # vitreous elimination (equivalently the first-order rate of systemic
    # absorption from the vitreous depot), and frac is the fraction of each
    # intravitreal dose that reaches the systemic circulation via a rapid
    # needle-track shunt bypassing the vitreous. The typical fraction estimate
    # of 2.51% is reported in Table 3 (theta4 = 0.0251); the paper mentions
    # 2.35% as the pre-final estimate in the base-model description.
    lcl  <- log(24.1)  ; label("Apparent clearance CL/F (L/day)")                     # Table 3, theta1 = 24.1 (RSE 4.52%)
    lvc  <- log(3.01)  ; label("Apparent central volume Vc/F (L)")                    # Table 3, theta2 = 3.01 (RSE 13.3%)
    lka  <- log(0.0806); label("Vitreous elimination / systemic absorption rate Ka (1/day)")  # Table 3, theta3 = 0.0806 (RSE 7.33%)
    lfrac <- log(0.0251); label("Fraction of ITV dose bypassing vitreous depot via needle-track shunt (unitless, in (0,1))") # Table 3, theta4 = 0.0251 (RSE 23.2%)

    # Covariate effects (Xu 2013 Table 3, theta5 and theta6). CRCL enters
    # CL/F as a power on the median-centered ratio; verteporfin PDT enters
    # Ka as the fractional-reduction magnitude in the multiplier
    # (1 - theta6 * CONMED_VERTEPORFIN), giving Ka * (1 - 0.353) = 0.647 * Ka
    # for subjects with concomitant PDT (35.3% lower Ka; Xu 2013 Results
    # narrative). Table 3 stores theta6 = 0.353 as a positive number with
    # the header "Covariate multiplier for Ka"; the subtractive formulation
    # in model() preserves that stored value.
    e_crcl_cl                <- 0.303 ; label("Power exponent for CRCL on CL/F (unitless)")             # Table 3, theta5 = 0.303 (RSE 32.0%)
    e_conmed_verteporfin_ka  <- 0.353 ; label("Fractional reduction in Ka on CONMED_VERTEPORFIN (unitless)")  # Table 3, theta6 = 0.353 (RSE 20.7%)

    # Inter-individual variability -- Xu 2013 reports omega values as
    # inter-individual CV%. Convert to log-normal variance via
    # omega^2 = log(1 + CV^2).
    etalcl   ~ 0.0940   # log(1 + 0.314^2) = 0.0940; Table 3, omega_CL/F = 31.4% (RSE 22.0%)
    etalvc   ~ 0.386    # log(1 + 0.686^2) = 0.386;  Table 3, omega_Vc/F = 68.6% (RSE 27.2%)
    etalka   ~ 0.0464   # log(1 + 0.218^2) = 0.0464; Table 3, omega_Ka   = 21.8% (RSE 96.2%)
    etalfrac ~ 0.580    # log(1 + 0.887^2) = 0.580;  Table 3, omega_frac = 88.7% (RSE 34.8%)

    # Residual variability (Xu 2013 Table 3). Combined additive + proportional
    # on the linear concentration scale (`ng/mL`).
    addSd  <- 0.145 ; label("Additive residual SD (ng/mL)")               # Table 3, theta7 = 0.145 (RSE 42.8%)
    propSd <- 0.334 ; label("Proportional residual SD (fraction)")        # Table 3, theta8 = 33.4% (RSE 8.44%)
  })

  model({
    # 1. Individual parameters
    # CL/F carries the CRCL power effect on the median-centered ratio.
    # Ka carries the verteporfin-PDT linear multiplier (1 + theta6 * PDT).
    # frac is bounded in (0, 1) at the typical value; log-normal sampling on
    # the log scale may cross 1 in the extreme upper tail (see vignette
    # Assumptions and deviations); in that case the fdepot bioavailability
    # would fall below 0, so upstream users of very-large-eta sims should
    # cap frac at the (0, 1) boundary.
    cl   <- exp(lcl  + etalcl)  * (CRCL / 65.22) ^ e_crcl_cl
    vc   <- exp(lvc  + etalvc)
    ka   <- exp(lka  + etalka)  * (1 - e_conmed_verteporfin_ka * CONMED_VERTEPORFIN)
    frac <- exp(lfrac + etalfrac)

    # 2. Micro-constants
    kel  <- cl / vc

    # 3. ODE system
    # Compartments: vitreous (paper: A_vit, drug amount in vitreous humor),
    #               central  (paper: A_sys, drug amount in systemic circulation).
    d/dt(vitreous) <- -ka * vitreous
    d/dt(central)  <-  ka * vitreous - kel * central

    # 4. Bioavailability -- split the dose across the two compartments so
    #    that (1 - frac) enters the slow vitreous depot and frac enters the
    #    central compartment directly via the needle-track shunt. In the
    #    event table, users provide the SAME nominal dose amount on both
    #    the vitreous and central dose rows at each ITV injection; the
    #    F() split allocates the correct proportion to each compartment.
    f(vitreous) <- 1 - frac
    f(central)  <-     frac

    # 5. Derived output: vitreous humor concentration (ng/mL). Xu 2013
    #    assumes a homogeneous vitreous volume of 4 mL in the population
    #    simulations (`with a vitreous volume of 4 mL`); we reproduce
    #    that assumption here so downstream simulations can display the
    #    vitreous concentration alongside serum. Units: vitreous amount
    #    is in mg by dose-unit convention; vitreous volume is 4 mL;
    #    (mg / mL) * 1e6 = ng/mL.
    Cvit <- vitreous / 4 * 1e6

    # 6. Observation and residual error -- serum total ranibizumab.
    #    Units: central amount is in mg; vc is in L; (mg / L) * 1000 = ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
