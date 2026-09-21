Cagnardi_2018_cefazolin_dog <- function() {
  description <- paste(
    "Clinical veterinary (dog).",
    "Two-compartment population PK model for cefazolin given as a single 25 mg/kg",
    "intravenous bolus 30 min before surgery to 78 client-owned dogs, parameterised",
    "per kg body weight in terms of serum clearance, intercompartmental clearance and",
    "the two volumes of distribution. Between-subject variability is exponential on all",
    "four disposition parameters with a full 4x4 variance-covariance matrix; the",
    "residual is combined proportional plus additive. Seven covariates (sex, age, body",
    "weight, breed, health status, serum creatinine and surgery duration) were screened",
    "by a stepwise BIC search and none was retained in the final model, so the typical",
    "values apply unadjusted across the canine surgical population. The measured",
    "unbound fraction (0.64) is carried so that free serum concentrations can drive the",
    "fT>MIC target-attainment analysis that set the 2 mg/L PK/PD cut-off",
    "(Cagnardi 2018)",
    sep = " "
  )
  reference <- paste(
    "Cagnardi P, Di Cesare F, Toutain P-L, Bousquet-Melou A, Ravasio G, Villa R.",
    "Population pharmacokinetic study of cefazolin used prophylactically in canine",
    "surgery for susceptibility testing breakpoint determination.",
    "Front Pharmacol. 2018;9:1137. doi:10.3389/fphar.2018.01137.",
    sep = " "
  )
  vignette <- "Cagnardi_2018_cefazolin_dog"

  # Every structural parameter is published per kg body weight (CL and Q in
  # L/kg/min, V1 and V2 in L/kg; Table 2), so doses are given in mg/kg and the
  # compartment amounts are carried in mg/kg. central/vc is then mg/L, which is
  # numerically the ug/mL of the HPLC-UV assay and of the additive residual SD,
  # so no scaling factor appears in Cc. Time is kept in minutes, the unit in
  # which the paper reports every rate and every half-life.
  units <- list(time = "min", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "cefazolin",
      units = "mg/kg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "cefazolin",
      units = "mg/kg",
      specimen = "serum",
      verified = TRUE
    )
  )

  # The stepwise BIC search returned 27 statistically significant covariate
  # scenarios (Results, 'Population Pharmacokinetics and Monte Carlo
  # Simulation'), but the authors then quantified each effect as the
  # multiplicative factor obtained when the covariate moved by +/-50% and judged
  # every one of them clinically irrelevant. The final model therefore carries no
  # covariate: the abstract states that "none of the seven explored covariates
  # were able to reduce this variability by an amplitude clinically relevant" and
  # the Conclusion that "no adjustment of dose for special dog populations seems
  # necessary". Table 2 -- the parameter table transcribed into ini() below -- is
  # explicitly the covariate-free model. All seven screened covariates are
  # documented here rather than encoded.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on V1, V2, CL and Q as a power model normalised to a 20 kg scaling",
        "weight (Methods eq. 10; Results: 'the scaling values were BW = 20 kg'). The most",
        "frequently selected covariate in the stepwise search - 16 of the 27 significant",
        "two-covariate scenarios included body weight - but not retained. Results reports",
        "the exponent for body weight on clearance alone as -0.2368, i.e. a 10 kg dog has",
        "1.178 times and a 30 kg dog 0.908 times the typical clearance of a 20 kg dog,",
        "'not relevant from a clinical point of view'. NOTE this exponent is the RESIDUAL",
        "weight effect: the disposition parameters of Table 2 are already normalised per",
        "kg, so a per-kg dose already carries the dominant weight dependence. Table 1",
        "range 4.5-56 kg, median 27 kg.",
        sep = " "
      ),
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on V1, V2, CL and Q as a power model normalised to an 8-year scaling",
        "age (Methods eq. 10; Results: 'age = 8 years'). Not retained: Results and",
        "Discussion both state that the influence of age on cefazolin exposure was not",
        "clinically relevant. Table 1 range 0.66-14 years, median 8 years.",
        sep = " "
      ),
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on V1, V2, CL and Q as a power model normalised to a 0.9 mg/dL scaling",
        "value (Methods eq. 10; Results: 'creatinine level = 0.9 mg/dL'). Used as the",
        "marker of individual kidney condition because 'cefazolin is mainly excreted by",
        "the kidney (approximately 80%, Nishida et al., 1970)'. Not retained; the",
        "influence was not clinically relevant. Supplementary Figure S5 shows a trend",
        "between creatinine level and both age and body weight. Table 1 range",
        "0.3-1.88 mg/dL, median 0.9 mg/dL. Units are mg/dL, not umol/L.",
        sep = " "
      ),
      source_name = "Creatinine level"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The paper coded sex as a THREE-level categorical covariate (Table 1): male",
        "n = 32 (code 0), female n = 23 (code 1), female neutered n = 23 (code 2). The",
        "canonical binary SEXF corresponds to codes 1 and 2 pooled, i.e. SEXF = 1 for the",
        "46 female dogs; the neutering split is carried separately as NEUTERED below.",
        "Screened on V1, V2, CL and Q; not retained. Results: 'For breed and sex, the",
        "magnitude of the effects was also clinically irrelevant despite their statistical",
        "significance.' No point estimate for the sex effect is published, so the effect",
        "cannot be reconstructed even for a user who wanted it.",
        sep = " "
      ),
      source_name = "Sex"
    ),
    NEUTERED = list(
      description = "Surgical neutering / gonadectomy status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sexually intact)",
      notes = paste(
        "Derived from the third level of the paper's sex covariate (Table 1, 'Female",
        "neutered n = 23 (Code 2)'). Identifiable only for the female dogs: the paper does",
        "not report castration status for the 32 males, so NEUTERED = 1 marks the 23",
        "neutered females and NEUTERED = 0 conflates intact females with males of unknown",
        "status. Documented for provenance only; it was never fitted as a covariate in its",
        "own right, only as a level of the three-level sex factor described under SEXF.",
        sep = " "
      ),
      source_name = "Sex (code 2)"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dog with concomitant disease)",
      notes = paste(
        "POLARITY IS INVERTED relative to the source coding. Table 1 codes health status",
        "as 'Healthy n = 19 (Code 0)' and 'Diseased n = 59 (Code 1)', so the canonical",
        "DIS_HEALTHY = 1 - (paper code). The diseased condition was defined on clinical",
        "exam, anamnesis and haematological/biochemical blood tests; about 60% of the",
        "diseased dogs were oncological patients. Screened on V1, V2, CL and Q and the",
        "only covariate that statistically influenced ALL four disposition parameters,",
        "but not retained. Results gives the single published point estimate: the largest",
        "effect was on Q (-0.267), 'meaning that in diseased dogs, the intercompartmental",
        "CL that likely reflects tissular blood flow, was decreased by 26.7% compared with",
        "that in control dogs'; the effects on V1, V2 and CL are described as clinically",
        "irrelevant and no numbers for them are printed.",
        sep = " "
      ),
      source_name = "Health status"
    ),
    BREED_PUREBRED = list(
      description = "Purebred (non-mongrel) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (mongrel)",
      notes = paste(
        "Table 1 codes breed as 'Mongrel n = 27 (Code 0)' and 'Other breeds n = 51",
        "(Code 1)', so the indicator direction here matches the paper's coding. The 51",
        "non-mongrel dogs spanned 26 breeds, of which 12 were represented by more than one",
        "dog and 14 by a single dog. Screened on V1, V2, CL and Q; not retained (Results:",
        "clinically irrelevant despite statistical significance, no point estimate",
        "published). There is no canonical register entry for dog breed: the covariate was",
        "not retained by the source model, so no canonical name was proposed for it and",
        "this name is documentation only.",
        sep = " "
      ),
      source_name = "Breed"
    ),
    DUR_SURGERY = list(
      description = "Duration of the surgical procedure",
      units = "min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on V1, V2, CL and Q as a power model normalised to an 80 min scaling",
        "value (Methods eq. 10; Results: 'surgery time = 80 min'). Not retained; the",
        "influence was not clinically relevant. Table 1 range 20-260 min, median 80 min,",
        "mean 87.63 +/- 58.09 min. There is no canonical register entry for surgical",
        "procedure duration: the covariate was not retained by the source model, so no",
        "canonical name was proposed for it and this name is documentation only.",
        sep = " "
      ),
      source_name = "Surgery Time"
    )
  )

  population <- list(
    species = "dog (client-owned Canis lupus familiaris; 27 mongrels and 51 dogs across 26 named breeds)",
    n_subjects = 78L,
    n_studies = 1L,
    n_observations = 629L,
    age_range = "0.66-14 years",
    age_median = "8 years",
    weight_range = "4.5-56 kg",
    weight_median = "27 kg",
    sex_female_pct = 58.97,
    disease_state = paste(
      "19 healthy dogs presenting for gynaecological or andrological surgery and 59 dogs",
      "with concomitant disease (about 60% of them oncological patients) presenting for",
      "procedures ranging from oncological to ophthalmic surgery",
      sep = " "
    ),
    dose_range = "25 mg/kg single intravenous bolus, given 30 min before surgery",
    regions = "Italy (University Veterinary Hospital of Milan)",
    renal_function = "Serum creatinine 0.91 +/- 0.32 mg/dL, range 0.3-1.88 mg/dL, median 0.9 mg/dL (Table 1)",
    co_medication = paste(
      "Perioperative anaesthetic protocol: alpha-2 agonist plus opioid (e.g. methadone)",
      "premedication, propofol induction, isoflurane in 100% oxygen for maintenance;",
      "meloxicam 0.2 mg/kg subcutaneously 30 min after the start of surgery",
      "(Methods, 'Sample Collection and Analysis')",
      sep = " "
    ),
    notes = paste(
      "Demographics from Table 1 and Results 'Animals and Cefazolin Concentrations'. Two",
      "to 11 samples (median 9) per dog were drawn between 5 and 480 min after dosing",
      "across 14 nominal sampling times, for 629 samples in total; there were no censored",
      "data (HPLC-UV, LOQ 0.2 ug/mL, LOD 0.00024 ug/mL). Serum protein binding measured",
      "in vitro by ultrafiltration was 36.2 +/- 5.3%. Sex: 32 male, 23 female, 23 female",
      "neutered, so the female percentage is (23 + 23) / 78. Surgery duration 87.63 +/-",
      "58.09 min, range 20-260 min, median 80 min. CAUTION: Table 1 prints body weight as",
      "'26.13 0.88' kg; a standard deviation of 0.88 kg is incompatible with the stated",
      "4.5-56 kg range and with a median of 27 kg, so the second figure is a standard",
      "error or a typographical error and only the range and median should be relied on.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Disposition. Table 2, 'Population primary parameters' -- bootstrap means.
    # Parameterisation is the paper's own: "Parameterization was in terms of
    # serum clearance (CL), intercompartmental CL(s) (Q) and volume(s) of
    # distribution (V) with V1, V2, CL, and Q being the primary estimated
    # parameters" (Methods, 'Population Pharmacokinetics and Monte Carlo
    # Simulation'). Every value is per kg body weight.
    #
    # The bootstrap MEAN column is transcribed rather than the median column,
    # because the paper quotes the means in its own abstract ("Population
    # primary parameter estimates V1, V2, CL, and Q were (typical value +/- SE)
    # 0.116 +/- 0.013 L/kg, 0.177 +/- 0.011 L/kg, 0.0037 +/- 0.0002 L/kg/min,
    # and 0.0103 +/- 0.0013 L/kg/min") and because the means -- not the medians
    # -- reproduce the published secondary parameters (see the vignette source
    # trace: Vss, Vz, MRT, AUC and Beta all land within 0.5%).
    # ------------------------------------------------------------------------
    lvc <- log(0.116)
    label("Central volume of distribution V1 (L/kg)")
    # Table 2, row 'tvV1' = 0.116 L/kg (SE 0.013, CV 11.36%, 95% CI 0.084-0.137)

    lvp <- log(0.177)
    label("Peripheral volume of distribution V2 (L/kg)")
    # Table 2, row 'tvV2' = 0.177 L/kg (SE 0.011, CV 6.01%, 95% CI 0.158-0.194)

    lcl <- log(0.0037)
    label("Serum clearance CL (L/kg/min)")
    # Table 2, row 'tvCL' = 0.0037 L/kg/min (SE 0.0002, CV 4.26%, 95% CI 0.0034-0.0040)

    lq <- log(0.0103)
    label("Intercompartmental clearance Q (L/kg/min)")
    # Table 2, row 'tvQ' = 0.0103 L/kg/min (SE 0.0013, CV 12.82%, 95% CI 0.0073-0.0123)

    # ------------------------------------------------------------------------
    # Protein binding. Measured in vitro by ultrafiltration, NOT estimated by
    # the population model, so it is fixed.
    # ------------------------------------------------------------------------
    fu <- fixed(0.64)
    label("Fraction of cefazolin unbound in serum (unitless)")
    # Results: serum protein binding 36.2 +/- 5.3%, i.e. an unbound fraction of
    # 0.638. The Monte Carlo section states the value the authors actually
    # applied -- "considering the average percentage of unbound drug calculated
    # (i.e., 0.64 ...)" -- and Table 5 confirms it exactly: each MIC divided by
    # 0.64 reproduces the printed total serum concentration to two decimal
    # places (0.25 -> 0.39, 0.5 -> 0.78, 1 -> 1.56, 2 -> 3.12, 4 -> 6.25,
    # 8 -> 12.5). 0.64 is therefore used rather than 0.638.

    # ------------------------------------------------------------------------
    # Between-subject variability. Methods eq. 3: Cl_i = theta_median *
    # exp(eta_i), i.e. exponential (log-normal) BSV, with 'V1, V2, and Q ...
    # modeled using equations of the same form'. The FULL variance-covariance
    # omega matrix was selected: 'inclusion of covariance terms prevented the
    # risk of biased estimation of the variance terms'.
    #
    # Values are the lower triangle of Table 3, read in the table's own row
    # order nV1, nV2, nCL, nQ -- which is the order of the eta block below.
    # Variances are the diagonal (bold in Table 3); the off-diagonals are
    # covariances, not correlations.
    #
    # Two independent checks confirm the transcription. (i) Converting each
    # diagonal with Methods eq. 4, CV% = 100 * sqrt(exp(omega^2) - 1), gives
    # 31.15, 38.03, 37.34 and 46.61%, exactly the four percentages printed in
    # the Table 3 footnote. (ii) Dividing each covariance by the square root of
    # the product of its two variances gives 0.89, 0.28, 0.56, 0.04, 0.38 and
    # 0.43, exactly the correlation matrix printed in the lower half of Table 3.
    # The matrix is positive definite (smallest eigenvalue 0.0033).
    #
    # Table 4 reports slightly different bootstrap BSVs (31.42, 42.70, 36.83 and
    # 46.12% CV) obtained from the bootstrap replicates rather than the single
    # full-data run. Table 3 is used here because only Table 3 publishes the
    # covariances, and a matrix mixing Table 4 diagonals with Table 3
    # off-diagonals would not be guaranteed positive definite.
    etalvc + etalvp + etalcl + etalq ~ c(
      0.092598,
      0.099641, 0.135063,
      0.031062, 0.073881, 0.130511,
      0.00518, 0.062008, 0.068201, 0.196614
    )

    # ------------------------------------------------------------------------
    # Residual error. Methods eq. 6: C(t) = f(theta, Time) * (1 + eps1) + eps2,
    # i.e. combined proportional plus additive on the linear scale, which is
    # exactly nlmixr2's add() + prop(). Both Phoenix sigmas are reported as
    # standard deviations: 'the additive sigma is reported as its SD, noted
    # stdev, with the same units as serum concentration (ug/mL) and the
    # multiplicative sigma is called multStdev'.
    # ------------------------------------------------------------------------
    propSd <- 0.257
    label("Proportional residual SD (fraction)")
    # Table 2, row 'tvCMultStdev' = 0.257 (SE 0.016, 95% CI 0.226-0.285). The
    # Table 2 legend reads it back as 'a coefficient of variation of 25.7%'.

    addSd <- 0.564
    label("Additive residual SD (ug/mL)")
    # Table 2, row 'stdev (sigma)' = 0.564 ug/mL (SE 0.166, 95% CI 0.314-0.943).
    # Of the same order as the 0.2 ug/mL LOQ of the HPLC-UV assay.
  })

  model({
    # ---- Individual parameters (Methods eq. 3, exponential BSV) -------------
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)

    # ---- Micro-constants ----------------------------------------------------
    # Written as cl/vc, q/vc and q/vp rather than from stored rate constants:
    # rxSolve() defaults to useLinCmt = TRUE and can silently drop peripheral1
    # when the transfer terms come straight from micro-constant parameters,
    # solving a one-compartment model instead. Deriving them from the clearance
    # and volume parameters keeps the closed-form and ODE solvers in agreement.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Two-compartment disposition, intravenous bolus into central --------
    # Methods eq. 1 defines the terminal slope as
    #   Beta = 0.5 * [ (Q/V1 + Q/V2 + CL/V1)
    #                  - sqrt((Q/V1 + Q/V2 + CL/V1)^2 - 4*(Q/V2)*(CL/V1)) ],
    # and eq. 2 gives Alpha = (Q/V2)*(CL/V1)/Beta. Those are the eigenvalues of
    # exactly this mammillary system: their sum is k12 + k21 + kel and their
    # product is k21*kel, which pins the structure with no ambiguity.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Observation --------------------------------------------------------
    # central is in mg/kg and vc in L/kg, so central/vc is mg/L == ug/mL, the
    # unit of the HPLC-UV assay and of addSd. Cc is the TOTAL serum
    # concentration that was measured and fitted; Cu is the free concentration
    # that drives the fT>MIC target of the Monte Carlo analysis (Table 5).
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd)
  })
}
