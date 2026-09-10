Yu_2025b_vancomycin <- function() {
  description <- "One-compartment intravenous population PK model for vancomycin in 42 neonates and infants treated in a US neonatal intensive care unit (Yu 2025). Clearance is the product of a fully-mature typical value (0.46 L/h at 2.8 kg dosing weight and 0.3 mg/dL serum creatinine), a Hill maturation fraction in postmenstrual age (T50 42.6 weeks, Hill coefficient 2.24), an allometric weight term with the exponent fixed at 0.75, and the reciprocal power term (0.3 / SCr)^0.543 so clearance falls as creatinine rises. Volume of distribution is 2.12 L at 2.8 kg and scales linearly with weight. Interoccasion variability was retained on clearance alongside interindividual variability on clearance and volume. The paper's parallel machine-learning and metabolomics analysis predicts the same individual clearance estimates from clinical covariates and contributes nothing to this structural model."
  reference <- "Yu H, Xiao J, Zhu HJ. Predicting vancomycin clearance in neonates and infants by integrating machine learning and metabolomics with population pharmacokinetics. Clin Transl Sci. 2025;18(7):e70293. doi:10.1111/cts.70293"
  vignette <- "Yu_2025b_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given as an intermittent 60-minute
  # intravenous infusion (Yu 2025 Results 3.1), so the dose enters `central`
  # directly and there is no depot state. The specimen is verified: Methods
  # 2.1.1 states "at least one recorded serum vancomycin concentration" and
  # "Serum vancomycin levels were quantified by Michigan Medicine Laboratories
  # using a kinetic immunoassay method", so the sampled matrix is SERUM.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current dosing weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying: Yu 2025 Methods 2.1.2 lists weight among the covariates that 'were tested by incorporating them as regressors in the Monolix software', which is Monolix's mechanism for a time-varying covariate. Referenced to the cohort MEDIAN dosing weight of 2.8 kg, not the 70 kg adult standard (Results 3.2: 'CL and Vd were standardized to a median weight of 2.8 kg using an allometric scaling model'), so the Table 2 thetas 0.46 L/h and 2.12 L are already neonate-sized. Both allometric exponents were FIXED, not estimated -- 0.75 for CL and 1 for Vd (Results 3.2, citing the Anderson & Holford allometry reference) -- and neither appears in Table 2, which lists only estimated parameters with an S.E. Weight is the only retained covariate on Vd. Cohort dosing weight median 2.80 kg, 5th-95th percentile 0.654-6.53 kg (Table 1). Table 1 calls the column 'Dose weight (kg)' and Methods 2.1.1 calls it 'current dosing weight (weight)'; the printed CL and Vd equations use the bare symbol 'weight'.",
      source_name        = "weight"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months. Yu 2025 writes the Hill maturation function on CL directly in weeks (PMA_T50 = 42.6 weeks, PMA_Hill = 2.24; Results 3.2 CL equation and Table 2), and the register explicitly permits a weeks-scaled PAGE for models whose source equations are written that way. Time-varying (Methods 2.1.2 lists PMA among the regressor-coded time-varying covariates). The maturation term is the BARE Hill fraction PMA^2.24 / (42.6^2.24 + PMA^2.24), NOT normalised to a reference PMA, so it equals 0.5 exactly at PMA 42.6 weeks and 0.468 at the cohort median PMA of 40.2 weeks -- the Table 2 Cl_pop of 0.46 L/h is therefore the fully-mature clearance, roughly twice the typical clearance of a median subject. Cohort PMA median 40.2 weeks, 5th-95th percentile 26.5-66.3 weeks (Table 1). PMA is strongly collinear with weight in this cohort (r = 0.93, Results 3.1 and Figure S1), which the Discussion invokes to explain why weight -- retained by the population PK covariate search -- did not rank in the top 10 features of any ensemble machine-learning model.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "mg/dL, NOT the SI umol/L (divide umol/L by 88.4 to convert). Table 1 reports 'Serum creatinine (mg/dL)' with median 0.305 and 5th-95th percentile 0.113-1.18. Time-varying (Methods 2.1.2 regressor list). Enters CL as the RECIPROCAL power ratio (0.3 / SCr)^0.543 as printed in the Results 3.2 CL equation, so clearance FALLS as creatinine rises; this is identical to (CREAT / 0.3)^(-0.543), and the model keeps the paper's printed orientation so that the Table 2 estimate 0.543 (row SCr_pop) appears verbatim with a positive sign. The reference constant is the 0.3 mg/dL printed inside the equation, which is the cohort median 0.305 rounded to one significant figure; the model uses the printed 0.3, not 0.305. The assay is a Michigan Medicine Laboratories clinical creatinine measurement; the paper does not state whether it is Jaffe or enzymatic, so no inter-assay conversion is applied. Together with PMA, SCr was one of the two covariates that both the population PK covariate search and every ensemble machine-learning feature-importance ranking selected (Results 3.4, Figure 4).",
      source_name        = "SCr"
    ),
    OCC = list(
      description        = "Occasion index for the interoccasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Decomposed inside model() into the indicators occ1..occ7 that select the per-occasion IOV etas on CL. Yu 2025 defines an occasion explicitly: Methods 2.1.1, 'A new occasion was defined as a discontinuation of dosing exceeding 4 days'. Results 3.1 gives the observed range directly -- 18 of 42 patients (43%) experienced multiple treatment occasions, range 2-7, with a mean of 2.6 vancomycin concentrations per occasion -- so SEVEN occasion slots are encoded here, which covers the whole observed range rather than a guess. A user with more occasions can extend the pattern by adding further etaiov_cl_<k> slots at the same fixed variance. Adding IOV on CL to the base model dropped the objective function from 1447.63 to 1347.04 (Results 3.2), which is why it is in the final model. Records outside 1..7 contribute no IOV, i.e. they behave as the typical occasion.",
      source_name        = "occasion"
    )
  )

  # Screened in the Yu 2025 stepwise forward-selection / backward-elimination
  # covariate search (Methods 2.1.2: "each covariate was tested individually";
  # the collected covariate list is in Methods 2.1.1) but NOT retained in the
  # final model, so they are documentation only and are not referenced in
  # model(). Only SCr, PMA and weight survived on CL, and only weight on Vd
  # (Results 3.2, with the selection steps in Table S2). Every entry below was
  # also one of the ten clinical covariates fed to the machine-learning models
  # in Results 3.4.
  covariatesDataExcluded <- list(
    WT_BIRTH = list(
      description = "Birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained; current dosing weight was the retained size descriptor. Cohort median 0.835 kg, 5th-95th percentile 0.495-3.37 kg (Table 1). One of the ten clinical machine-learning input features (Results 3.4).",
      source_name = "BW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened but not retained; it is a component of the retained PAGE. Cohort median 28.3 weeks, 5th-95th percentile 23.9-38.3 weeks (Table 1). One of the ten clinical machine-learning input features (Results 3.4).",
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "weeks",
      type        = "continuous",
      notes       = "WEEKS in this paper, not the register-default months -- Table 1 reports 'Postnatal Age (week)' with median 8.64 and 5th-95th percentile 0.430-37.9. Screened but not retained; it is a component of the retained PAGE. One of the ten clinical machine-learning input features (Results 3.4).",
      source_name = "PNA"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 20 of 42 female (47.6%), 22 male (52.4%) (Table 1). One of the ten clinical machine-learning input features (Results 3.4), where the paper calls the column 'gender'.",
      source_name = "gender"
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "g/dL as reported by Yu 2025 Table 1 ('Albumin (g/dL)', median 3.00, 5th-95th percentile 2.21-4.00), NOT the register-canonical g/L -- multiply by 10 to convert. Screened but not retained. One of the ten clinical machine-learning input features (Results 3.4).",
      source_name = "ALB"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "mg/dL. Cohort median 19.0, 5th-95th percentile 6.08-65.6 (Table 1, 'Urea nitrogen (mg/dL)'). Screened but not retained; it was highly correlated with the retained SCr (r = 0.71, Results 3.1 and Figure S1). One of the ten clinical machine-learning input features (Results 3.4).",
      source_name = "BUN"
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 31 of 42 (73.8%) (Table 1). Race was one of the ten clinical machine-learning input features (Results 3.4); the paper does not state how the three-level race column was numerically encoded for the machine-learning models.",
      source_name = "race"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 9 of 42 (21.4%) (Table 1).",
      source_name = "race"
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 2 of 42 (4.76%) (Table 1).",
      source_name = "race"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    age_range      = "Gestational age 23.9-38.3 weeks (5th-95th percentile; median 28.3); postnatal age 0.430-37.9 weeks (median 8.64); postmenstrual age 26.5-66.3 weeks (median 40.2). The cohort is 25 neonates and 17 infants.",
    weight_range   = "Current dosing weight 0.654-6.53 kg (5th-95th percentile; median 2.80); birth weight 0.495-3.37 kg (median 0.835)",
    sex_female_pct = 47.6,
    race_ethnicity = "White 31 (73.8%), Black 9 (21.4%), Other 2 (4.76%) (Table 1)",
    disease_state  = "Neonates and infants admitted to the University of Michigan Neonatal Intensive Care Unit between 2019 and 2022 and treated with intravenous vancomycin for suspected or confirmed Gram-positive bacterial infection. Inclusion required at least one recorded serum vancomycin concentration; no renal-function exclusion is stated.",
    dose_range     = "3.5-25 mg/kg per dose, given every 6, 8, 12, 18 or 24 h as a 60-minute intravenous infusion (Results 3.1)",
    regions        = "United States (single centre: University of Michigan Neonatal Intensive Care Unit, Ann Arbor, MI)",
    renal_function = "Serum creatinine 0.113-1.18 mg/dL (5th-95th percentile; median 0.305) and blood urea nitrogen 6.08-65.6 mg/dL (median 19.0) (Table 1). Renal replacement therapy is not mentioned as an exclusion criterion, and no renal-impairment stratum is defined, so the model's domain of applicability is the observed creatinine range.",
    notes          = "Retrospective single-centre electronic-medical-record study approved by the University of Michigan IRB. 214 serum vancomycin concentrations from 42 patients, a mean of 5 per patient (range 1-22), predominantly steady-state troughs with some peak and random levels; the assay LLOQ was 4.0 ug/mL and below-LLOQ records were excluded from the population PK analysis. 18 of 42 patients (43%) contributed multiple occasions (range 2-7), an occasion being a dosing gap of more than 4 days, at a mean of 2.6 concentrations per occasion. Estimated by SAEM in Monolix 2024R1. A two-compartment model was tested and rejected: it did not lower the objective function and the R.S.E.% for V1, Q and V2 were much higher (Results 3.2). Approximately 47.7% of the analysed trough concentrations lay in the 10-20 ug/mL target range (Results 3.1). Vd is poorly informed by these largely trough-only data -- its shrinkage is 52.7% versus 15.8% for CL (Table 2) -- so individual Vd predictions from this model are close to the population mean. The paper's second half compares eleven machine-learning regressors, trained on the model's own empirical-Bayes CL estimates, using clinical covariates and/or untargeted plasma metabolomics; the best (gradient boosting on the ten clinical covariates) reached R^2 0.830, metabolomics added nothing, and none of that analysis alters the structural model carried here."
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters -- Yu 2025 Table 2 (final model), in the
    # parameterisation of the two Results 3.2 equations. The reference subject
    # is a neonate at the cohort MEDIAN dosing weight of 2.8 kg with serum
    # creatinine 0.3 mg/dL.
    #
    # NOTE on lcl: 0.46 L/h is the theta multiplying the BARE Hill maturation
    # FRACTION, i.e. the fully-mature (PMA -> infinity) clearance at 2.8 kg and
    # SCr 0.3 mg/dL, not the clearance of any subject in the studied PMA range.
    # At the cohort median PMA of 40.2 weeks the Hill fraction is 0.468, so the
    # typical median subject clears 0.215 L/h (0.077 L/h/kg). This conflicts
    # with the Discussion, which states "the estimated median vancomycin CL was
    # 0.16 L/h/kg for this cohort" -- a figure that is exactly 0.46 / 2.8 and so
    # drops both the maturation and the creatinine factor. The printed equation
    # is the correct reading: it puts a median subject's steady-state trough on
    # 15 mg/kg q12h inside the paper's own 10-20 ug/mL target window, whereas
    # the Discussion reading halves every exposure and lands the same subject
    # near 5 ug/mL, which cannot be reconciled with Results 3.1 ("Approximately
    # 47.7% of the analyzed trough concentrations were within the target
    # concentration range of 10-20 ug/mL"). See the vignette Errata.
    # ------------------------------------------------------------------------
    lcl <- log(0.46); label("Fully mature clearance at 2.8 kg dosing weight and 0.3 mg/dL serum creatinine (L/h)")  # Yu 2025 Table 2 (Cl_pop = 0.460 L/h, S.E. 0.0423, R.S.E. 9.21%, 95% CI 0.384-0.550) and Results 3.2 CL equation
    lvc <- log(2.12); label("Volume of distribution at 2.8 kg dosing weight (L)")                                   # Yu 2025 Table 2 (V_pop = 2.12 L, S.E. 0.286, R.S.E. 13.5%, 95% CI 1.63-2.75) and Results 3.2 Vd equation

    # Allometric exponents on current dosing weight, referenced to the cohort
    # median 2.8 kg. Both were FIXED, not estimated -- Results 3.2: "with the
    # power exponent fixed at 1 for Vd and 0.75 for CL". Neither appears in
    # Table 2, which lists only estimated parameters with an S.E.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT / 2.8 kg) for CL (unitless)")  # Yu 2025 Results 3.2 (fixed)
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT / 2.8 kg) for Vd (unitless)")  # Yu 2025 Results 3.2 (fixed)

    # Hill maturation of clearance with postmenstrual age, in WEEKS. Results
    # 3.2: "PMA was incorporated using the Hill function in the covariate
    # model"; the printed CL equation gives the bare fraction
    # PMA^2.24 / (42.6^2.24 + PMA^2.24), with no reference-PMA normalisation.
    pma_tm50 <- 42.6; label("Postmenstrual age at 50% of mature CL (PMA_T50, weeks)")  # Yu 2025 Table 2 (PMA_T50 = 42.6 weeks, S.E. 3.57, R.S.E. 8.38%, 95% CI 36.2-50.2)
    pma_hill <- 2.24; label("Hill coefficient for CL maturation with postmenstrual age (unitless)")  # Yu 2025 Table 2 (PMA_Hill = 2.24, S.E. 0.0354, R.S.E. 1.58%, 95% CI 2.17-2.31)

    # Serum-creatinine effect on clearance. The printed Results 3.2 CL equation
    # uses the RECIPROCAL ratio (0.3 / SCr)^0.543, so a POSITIVE exponent means
    # clearance falls as creatinine rises; identical to (CREAT / 0.3)^(-0.543).
    e_creat_cl <- 0.543; label("Power exponent on (0.3 mg/dL / CREAT) for CL (unitless)")  # Yu 2025 Table 2 (SCr_pop = 0.543, S.E. 0.0746, R.S.E. 13.7%, 95% CI 0.417-0.708) and Results 3.2 CL equation

    # ------------------------------------------------------------------------
    # Interindividual variability, exponential (Methods 2.1.2: "All individual
    # parameters were assumed to be log-normally distributed. Inter-individual
    # variability (IIV) and inter-occasion variability (IOV) were modeled for
    # each PK parameter using an exponential model").
    #
    # Table 2's "Standard deviation of the random effects" block reports the
    # omegas on the LOG SCALE and prints the apparent CV alongside, which pins
    # the convention with no ambiguity: sqrt(exp(0.284^2) - 1) = 0.290 = the
    # printed C.V. 29.0% for V, and sqrt(exp(0.195^2) - 1) = 0.197 = the
    # printed 19.6% for CL. The variances below are therefore the squares of
    # the printed values, not log(1 + CV^2).
    #
    # No CL-V correlation is reported, so the etas are independent.
    # ------------------------------------------------------------------------
    etalcl ~ 0.038025  # Yu 2025 Table 2 (omega_Cl = 0.195, C.V. 19.6%, R.S.E. 21.8%, 95% CI 0.129-0.294); 0.195^2
    etalvc ~ 0.080656  # Yu 2025 Table 2 (omega_V  = 0.284, C.V. 29.0%, R.S.E. 35.4%, 95% CI 0.150-0.540); 0.284^2

    # ------------------------------------------------------------------------
    # Interoccasion variability on clearance, exponential (the kappa term of
    # the Results 3.2 CL equation, printed as e^(eta_CL,i + kappa_CL,ik)).
    # Table 2 reports gamma_Cl = 0.161 on the same log SD scale as the two IIV
    # rows above -- sqrt(exp(0.161^2) - 1) = 0.162 reproduces the printed C.V.
    # 16.2% -- so the variance is 0.161^2 = 0.025921.
    #
    # Implemented by the occasion-indicator expansion rather than rxode2's
    # `~ var | OCC` syntax, which parses but cannot be simulated from an rxUi.
    # Occasions after the first repeat the same variance, the analogue of
    # NONMEM's $OMEGA BLOCK(1) SAME, and are therefore fixed().
    #
    # Seven slots, because Results 3.1 states the observed occasion count
    # directly: "18 (43%) experienced multiple treatment occasions (range:
    # 2-7)". This is a transcribed bound, not an assumption.
    # ------------------------------------------------------------------------
    etaiov_cl_1 ~ 0.025921         # Yu 2025 Table 2 (gamma_Cl = 0.161, C.V. 16.2%, R.S.E. 22.5%, 95% CI 0.105-0.246); 0.161^2
    etaiov_cl_2 ~ fixed(0.025921)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.025921)  # SAME-equivalent
    etaiov_cl_4 ~ fixed(0.025921)  # SAME-equivalent
    etaiov_cl_5 ~ fixed(0.025921)  # SAME-equivalent
    etaiov_cl_6 ~ fixed(0.025921)  # SAME-equivalent
    etaiov_cl_7 ~ fixed(0.025921)  # SAME-equivalent

    # ------------------------------------------------------------------------
    # Residual error, combined additive plus proportional (Results 3.2: "The
    # proportional-plus-additive model was selected for the residual error
    # model"; Table 2's note defines "a: The additive error (ug/mL), b: The
    # proportional error (unitless)").
    #
    # The paper does not say which of Monolix's two combined error models it
    # used. `add() + prop()` below is nlmixr2's default combined2 form,
    # sd = sqrt(a^2 + (b*f)^2); Monolix's combined1 form, sd = a + b*f, is the
    # alternative reading and gives a larger residual SD (4.09 vs 2.89 ug/mL at
    # f = 15 ug/mL). Only the residual noise is affected -- no structural
    # prediction changes. See the vignette Errata.
    # ------------------------------------------------------------------------
    addSd  <- 2.15;  label("Additive residual error (ug/mL)")     # Yu 2025 Table 2 (a = 2.15 ug/mL, S.E. 0.244, R.S.E. 11.4%, 95% CI 1.72-2.68)
    propSd <- 0.129; label("Proportional residual error (fraction)")  # Yu 2025 Table 2 (b = 0.129, S.E. 0.0274, R.S.E. 21.4%, 95% CI 0.0862-0.193)
  })
  model({
    # ----------------------------------------------------------------------
    # 1. Derived covariate terms -- Yu 2025 Results 3.2, the printed CL and Vd
    #    equations:
    #
    #      CL = 0.46 * (PMA^2.24 / (42.6^2.24 + PMA^2.24))
    #                * (weight / 2.8)^0.75
    #                * (0.3 / SCr)^0.543
    #                * exp(eta_CL,i + kappa_CL,ik)
    #
    #      Vd = 2.12 * (weight / 2.8) * exp(eta_Vd,i)
    #
    #    Reference constants: median dosing weight 2.8 kg, reference serum
    #    creatinine 0.3 mg/dL. PAGE is in WEEKS for this model (see
    #    covariateData$PAGE$notes).
    # ----------------------------------------------------------------------
    maturation_cl <- PAGE^pma_hill / (pma_tm50^pma_hill + PAGE^pma_hill)
    creat_cl      <- (0.3 / CREAT)^e_creat_cl

    # Interoccasion variability on clearance, selected by the occasion index
    # OCC (1..7; see covariateData$OCC$notes).
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    occ7 <- (OCC == 7)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7

    # 2. Individual parameters
    cl <- exp(lcl + etalcl + iov_cl) * maturation_cl * (WT / 2.8)^e_wt_cl * creat_cl
    vc <- exp(lvc + etalvc) * (WT / 2.8)^e_wt_vc

    # 3. Micro-constant
    kel <- cl / vc

    # ----------------------------------------------------------------------
    # 4. ODE system -- one compartment with first-order elimination and
    #    intravenous dosing only (Results 3.2: "the one-compartment model with
    #    first-order elimination was selected"; vancomycin was given as an
    #    intermittent 60-minute IV infusion).
    # ----------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # 5. Observation and error. Dose in mg, vc in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
