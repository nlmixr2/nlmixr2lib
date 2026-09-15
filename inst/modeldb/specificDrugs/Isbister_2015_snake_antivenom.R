Isbister_2015_snake_antivenom <- function() {
  description <- "Two-compartment population PK model for Indian polyvalent F(ab')2 snake antivenom (VINS Bioproducts Ltd) in adults with Russell's viper (Daboia russelii) envenoming (Isbister 2015): zero-order intravenous input, linear elimination from the central compartment, and a power effect of body weight on central volume. Relative bioavailability is fixed to 1 with between-subject variability estimated; that random effect absorbs the per-patient uncertainty in the delivered antivenom dose caused by variable losses during reconstitution of the freeze-dried vials. Fit in MONOLIX 4.2 (SAEM, M3 handling of below-limit-of-quantification data) to 411 quantifiable antivenom concentrations from 75 patients. The authors selected a combined (additive plus proportional) residual-error model but publish no residual-error magnitudes, so both are encoded as zero."
  reference <- "Isbister GK, Maduwage K, Saiao A, Buckley NA, Jayamanne SF, Seyed S, et al. Population pharmacokinetics of an Indian F(ab')2 snake antivenom in patients with Russell's viper (Daboia russelii) bites. PLoS Negl Trop Dis. 2015;9(7):e0003873. doi:10.1371/journal.pntd.0003873"
  vignette <- "Isbister_2015_snake_antivenom"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the central volume only, centred on the cohort average weight: V = theta_V * (WT/WTav)^f_wt (Isbister 2015, Methods 'Pharmacokinetic analysis', unnumbered equation). The paper says the covariate was 'centred to the average weight' but never prints the average; Table 1 reports only the MEDIAN weight of 57 kg (range 40 to 70 kg), so 57 kg is used as WTav here. See the vignette 'Assumptions and deviations' section. Age, sex and pre-antivenom venom concentration were screened by visual inspection of individual parameter estimates and were NOT retained (Methods and Results); they are listed in covariatesDataExcluded.",
      source_name        = "wt"
    )
  )

  # Covariates the source screened but did not retain in the final model.
  # Documentation only -- these are deliberately not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = "Isbister 2015 Methods: 'Age, sex and pre-antivenom concentrations were not included in the final model evaluation due to the absence of an association visually.' Cohort median 38 years (range 16 to 64)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Isbister 2015 Methods: screened by visual inspection of individual parameter estimates, no association seen, not carried into the final model. Cohort was 64/75 (85%) male."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Antivenom was measured by sandwich EIA in SERUM (Methods, 'Data
    # collection': blood collected in serum tubes for the antivenom EIA).
    central     = list(analyte = "Indian polyvalent F(ab')2 snake antivenom", units = "mg", specimen = "serum", verified = TRUE),
    # Mathematical distribution compartment; the paper assigns it no anatomical
    # identity beyond noting that V + Vp is 'consistent with a large molecule
    # which does not have a large volume of distribution' (Discussion).
    peripheral1 = list(analyte = "Indian polyvalent F(ab')2 snake antivenom", units = "mg", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 75,
    n_studies      = 1,
    age_range      = "16 to 64 years",
    age_median     = "38 years",
    weight_range   = "40 to 70 kg",
    weight_median  = "57 kg",
    sex_female_pct = 14.7,
    race_ethnicity = "Not reported; single-centre Sri Lankan cohort.",
    disease_state  = "Snake envenoming with coagulopathy (abnormal 20-minute whole blood clotting test). 71 of 75 Russell's viper (Daboia russelii) envenoming, of whom 52 had detectable venom before antivenom; 4 hump-nosed viper (Hypnale spp.) envenoming. Local envenoming 97%, coagulopathy 100%, systemic bleeding 35%, neurotoxicity (ptosis) 43%.",
    dose_range     = "8 to 40 vials (median 18) of Indian polyvalent antivenom given intravenously; each 10-vial dose is reconstituted in 100 mL and administered in a total of 500 mL normal saline over 1 hour. 21 of 75 patients (28%) received a repeat dose. The paper does NOT report the antivenom mass per vial, so the model's dosing unit is the assay's calibrator mass -- see the vignette for the Figure 2 back-solve.",
    regions        = "Sri Lanka (Base Hospital Polonnaruwa, Central Eastern Province)",
    notes          = "Isbister 2015 Table 1. Patients >15 years old recruited October 2010 to March 2012 from a prospective snakebite cohort and enrolled in a fresh-frozen-plasma dose-finding randomised trial; the PK sub-study required serial serum sampling and complete demographics. 510 antivenom samples were drawn, of which 411 had quantifiable antivenom (limit of quantification 40 ug/mL); median 5 samples per patient for the 54 single-dose patients and 7 for the 21 multiple-dose patients. Five antivenom batches were used (1060, 1096, 1102, 01015/10-11, 01AS11112); S1 Fig shows no relationship between batch and the F random effect."
  )

  ini({
    # Structural parameters -- Isbister 2015 Table 2, column
    # 'Model 3 (Final) Including F and weight on V'. MONOLIX 4.2 estimates.
    lcl <- log(0.0779); label("Clearance (L/h)")                              # Table 2 Model 3: CL = 0.0779 L/h (rse 34%)
    lvc <- log(2.16);   label("Central volume of distribution (L)")           # Table 2 Model 3: V  = 2.16 L (rse 10%)
    lq  <- log(0.178);  label("Intercompartmental clearance (L/h)")           # Table 2 Model 3: Q  = 0.178 L/h (rse 31%)
    lvp <- log(8.33);   label("Peripheral volume of distribution (L)")        # Table 2 Model 3: Vp = 8.33 L (rse 52%)

    # Power exponent on body weight for the central volume, centred on the
    # cohort average weight. Methods equation: V = theta_V * (wt/wt_av)^f_wt.
    e_wt_vc <- 0.132; label("Power exponent on (WT/57 kg) for central volume (unitless)")   # Table 2 Model 3: f_wt = 0.132 (rse 84%)

    # Relative bioavailability of the intravenous input. The authors fixed F
    # to 1 and estimated only its between-subject variability, deliberately
    # using that random effect as a per-patient dose-uncertainty term
    # (Methods: 'F was fixed to 1 and the BSV was estimated for each patient
    # similar to including uncertainty on dose'). There is no depot
    # compartment in this model; lfdepot is the library's canonical name for
    # a bioavailability multiplier and is applied to the central compartment
    # here. See the vignette 'Assumptions and deviations' section.
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the intravenous input (unitless)")  # Table 2 Model 3: F = 1 (fixed, no rse reported)

    # Between-subject variability. Isbister 2015 Table 2 block headed
    # 'Between subject variance (omega)'. The values are MONOLIX 4.2
    # omega_<parameter> outputs, which are the STANDARD DEVIATIONS of the
    # log-normal random effects, so the nlmixr2 variance is omega^2. The
    # table's 'variance' wording is an author slip; three independent checks
    # settle the scale:
    #   1. MONOLIX 4.2 reports omega_X as an SD, not a variance.
    #   2. Reading them as SDs reproduces the paper's own reported half-life
    #      distribution (Results): simulated medians 4.6 h (distribution) and
    #      133 h (elimination) against the published 4.6 h and 140 h, with a
    #      spread wider than the published 10th-90th percentiles only by the
    #      amount expected from shrinkage of the empirical Bayes estimates
    #      the paper actually tabulated. Reading them as variances
    #      over-disperses both half-lives by roughly a further 50%.
    #   3. omega_F = 0.197 as an SD is a 20% coefficient of variation in the
    #      delivered dose, which matches the authors' stated mechanism
    #      ('variable losses occurring during reconstitution of the
    #      individual freeze dried vials'); as a variance it would be a 47%
    #      CV, far larger than reconstitution losses can plausibly be.
    # SINGLE quotes only in these comments -- rxode2 promotes an unlabelled
    # ini() trailing comment into label().
    etalcl     ~ 0.715^2  # Table 2 Model 3 row 'Cl' = 0.715 (rse 46%), MONOLIX omega (SD) -> variance 0.5112, 81.7% CV
    etalvc     ~ 0.188^2  # Table 2 Model 3 row 'V'  = 0.188 (rse 126%), MONOLIX omega (SD) -> variance 0.0353, 19.0% CV
    etalq      ~ 0.533^2  # Table 2 Model 3 row 'Q'  = 0.533 (rse 57%), MONOLIX omega (SD) -> variance 0.2841, 57.4% CV
    etalvp     ~ 0.836^2  # Table 2 Model 3 row 'Vp' = 0.836 (rse 125%), MONOLIX omega (SD) -> variance 0.6989, 99.0% CV
    etalfdepot ~ 0.197^2  # Table 2 Model 3 row 'F'  = 0.197 (rse 42%), MONOLIX omega (SD) -> variance 0.0388, 19.9% CV

    # Residual variability. Results: 'a combined error model best described
    # the data', i.e. MONOLIX's additive-plus-proportional form, which maps
    # to add(addSd) + prop(propSd) in nlmixr2. Table 2 tabulates the
    # structural thetas, the five omegas and the objective function but NO
    # residual-error row, and there is no supplement carrying one (S1 to S5
    # are goodness-of-fit and covariate-screening figures). Both magnitudes
    # are therefore encoded as zero rather than invented, so simulations from
    # this model are residual-error-free. See the vignette Errata.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Isbister 2015: combined error model selected, no magnitude published
    addSd  <- fixed(0); label("Additive residual SD (ug/mL; 0 -- not reported in the source)")         # Isbister 2015: combined error model selected, no magnitude published
  })

  model({
    # 1. Reference weight for the central-volume covariate. The paper centres
    #    on the cohort 'average weight' but prints only the median, 57 kg
    #    (Table 1). See covariateData$WT$notes.
    wtav <- 57  # kg, Isbister 2015 Table 1 median weight

    # 2. Individual parameters. Isbister 2015 Table 2 Model 3 plus the
    #    Methods covariate equation V = theta_V * (wt/wt_av)^f_wt.
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc) * (WT / wtav)^e_wt_vc
    q      <- exp(lq + etalq)
    vp     <- exp(lvp + etalvp)
    fdepot <- exp(lfdepot + etalfdepot)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Antivenom enters the central compartment directly by
    #    zero-order intravenous infusion; the infusion duration is a design
    #    variable supplied in the event table (rate/dur), not an estimated
    #    parameter.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # 5. Relative bioavailability of the intravenous input (dose-uncertainty
    #    random effect; typical value fixed to 1).
    f(central) <- fdepot

    # 6. Observation and error. amt in mg, vc in L, so central/vc is mg/L,
    #    which is the ug/mL unit of the antivenom sandwich EIA (Fig 1 and
    #    Fig 2 axis labels).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
