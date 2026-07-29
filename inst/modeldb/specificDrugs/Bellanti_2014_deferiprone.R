Bellanti_2014_deferiprone <- function() {
  description <- "One-compartment population PK model for the oral iron chelator deferiprone in healthy adult subjects, with first-order absorption, absorption lag time, and a binary sex effect on the apparent volume of distribution (Bellanti 2014)."
  reference <- "Bellanti F, Danhof M, Della Pasqua O. Population pharmacokinetics of deferiprone in healthy subjects. Br J Clin Pharmacol. 2014 Dec;78(6):1397-1406. doi:10.1111/bcp.12473"
  vignette <- "Bellanti_2014_deferiprone"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Final-model covariate on V/F: females show a ~17% lower apparent volume of distribution relative to males (Bellanti 2014 Table 1: V/F males = 78.4 L, V/F females = 65.3 L). Encoded as multiplicative effect e_sexf_vc on V_typical (reference category = male).",
      source_name        = "Gender"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight (kg).",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F and V/F during univariate covariate selection (Bellanti 2014 Results). Significant on both parameters when entered separately, but excluded from the final model because inclusion destabilised the bootstrap, attributed by the authors to the narrow weight range (52-92 kg) of the healthy-adult cohort. Documented here so the provenance of the screen is preserved without triggering a `declared but not referenced` convention warning."
    ),
    AGE = list(
      description = "Subject age (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in the covariate screen on CL/F and V/F (Bellanti 2014 Results). Not retained in the final model."
    ),
    CRCL = list(
      description = "Creatinine-based renal function (mL/min/1.73 m^2).",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Tested on CL/F in the covariate screen but the adult healthy cohort lacked the renal-function spread needed to estimate the effect (Bellanti 2014 Discussion: 'the population available for the analysis was limited to healthy subjects, the impact of another important covariate could not be estimated'). The renal-impairment dosing recommendations in Bellanti 2014 Table 2 come from a simulation exercise that reduces CL/F to 80%, 50% and 25% of the healthy-population value rather than from an estimated covariate effect; not a parameter of the structural model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55,
    n_studies      = 2,
    age_range      = "19-55 years",
    age_median     = "39 years",
    weight_range   = "52-92 kg",
    weight_median  = "72 kg",
    sex_female_pct = 16 / 55 * 100,
    disease_state  = "Healthy adult subjects",
    dose_range     = "Single oral dose of 1500 mg deferiprone as a 100 mg/mL solution (LA20-BA, LA21-BE studies). The model is also applied via simulation to 25 mg/kg single doses and 75 mg/kg/day b.i.d. dosing scenarios in the paper.",
    regions        = "Not reported in the publication; data supplied by ApoPharma Inc. (Canada) and shared within the DEEP consortium (FP7).",
    notes          = "39 males and 16 females. Up to 17 post-dose plasma samples per subject (median 15) over 14 h; analytical method HPLC-UV with LLOQ 1 microM (0.14 microg/mL). NONMEM v7.2 FOCE-I; bootstrap 500 runs in PsN v3.5.3. See Bellanti 2014 Methods - Data."
  )

  ini({
    # Structural parameters (typical values; estimated unless otherwise marked)
    lcl   <- log(30.8);  label("Apparent clearance CL/F (L/h)")               # Bellanti 2014 Table 1
    lvc   <- log(78.4);  label("Apparent volume of distribution V/F in males, reference (L)")  # Bellanti 2014 Table 1 (V/F males)
    lka   <- log(8.2);   label("First-order absorption rate constant Ka (1/h)")  # Bellanti 2014 Table 1
    ltlag <- log(0.146); label("Absorption lag time (h)")                      # Bellanti 2014 Table 1 (Lag time)

    # Sex effect on V/F: V_female / V_male - 1 = 65.3 / 78.4 - 1 = -0.16709
    # Encoded as multiplicative deviation from the male reference; reference category is male (SEXF = 0).
    e_sexf_vc <- -0.16709
    label("Sex effect on V/F: fractional change in V/F for females vs males (reference = male)")  # Bellanti 2014 Table 1, derived from V/F males = 78.4 L and V/F females = 65.3 L

    # IIV - block correlation between etalcl and etalvc; etalka separate.
    # Variances (omega^2) and covariance reported in Bellanti 2014 Table 1.
    # IIV CL/F: 0.057 (CV 23.87%), IIV V/F: 0.0278 (CV 16.67%), covariance CL-V: 0.0345
    # (column header `Correlation CL-V` refers to the NONMEM $OMEGA off-diagonal,
    # which is a covariance; implied correlation = 0.0345 / sqrt(0.057 * 0.0278) = 0.866).
    etalcl + etalvc ~ c(0.057,
                        0.0345, 0.0278)
    etalka ~ 0.991  # IIV Ka variance 0.991 (CV 99.54%) -- very wide, consistent with the small number of early-time samples that inform Ka.

    # Residual error: Bellanti 2014 Eq. 2 (Y_ij = F_ij + epsilon_ij * W) with
    # epsilon ~ N(0, sigma^2 = 0.00566) and the proportional weighting factor
    # theta_W = 2.4 reported on a separate Table 1 row. The standard NONMEM
    # 'proportional with a weighting factor' construction (W = theta_W * F)
    # yields an effective proportional residual SD of theta_W * sqrt(sigma^2)
    # = 2.4 * 0.0752 = 0.1805 (~18.05% CV). Encoded here as a single propSd
    # because nlmixr2 expresses proportional residual error in that form; the
    # decomposition is documented in the validation vignette's Assumptions and
    # deviations section.
    propSd <- 0.1805
    label("Effective proportional residual SD (fraction); theta_W * sigma = 2.4 * 0.0752")  # Bellanti 2014 Table 1 (Error: weighting factor = 2.4) and (Residual error variance = 0.00566)
  })

  model({
    # Sex effect applied multiplicatively to V/F. Reference category is male (SEXF = 0),
    # so for males the sex term is 1 and V/F = exp(lvc + etalvc).
    sex_vc <- 1 + e_sexf_vc * SEXF

    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc) * sex_vc
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
