Kastrissios_2012_managlinatDialanetil_linked <- function() {
  description <- paste(
    "Linked four-moiety population pharmacokinetic model for the oral",
    "hypoglycemic prodrug CS-917 (INN managlinat dialanetil) and its three",
    "metabolites in 141 patients with type 2 diabetes mellitus pooled from",
    "six phase I / IIa studies (Kastrissios 2012). CS-917 is a",
    "fructose-1,6-bisphosphatase inhibitor prodrug: an esterase converts it",
    "to the inactive R-134450, hepatic phosphoramidase converts that to the",
    "active moiety R-125338, and R-125338 is cleared by renal excretion in",
    "parallel with N-acetylation to the slowly forming, long-half-life",
    "inactive R-143047. Structurally the cascade is CS-917 one-compartment",
    "with first-order absorption and an absorption lag, R-134450",
    "two-compartment, R-125338 two-compartment, R-143047 one-compartment,",
    "each metabolite formed by first-order transfer of the whole upstream",
    "elimination flux. Every clearance and volume is apparent in the",
    "compound sense CL/(F*fm): the source states the metabolite terms are",
    "relative to CS-917 oral bioavailability AND the unknown fraction of",
    "CS-917 metabolized to the moiety, so the unknown fraction is absorbed",
    "into the downstream apparent volume and the inter-moiety flux carries",
    "no fractional multiplier. Retained covariates: body weight on CS-917",
    "CL/F and Vc/F, food on CS-917 relative bioavailability, age and Black",
    "race on R-134450 CL/F, creatinine clearance and food on R-125338 CL/F,",
    "Black race on R-125338 Vc/F, body weight on R-143047 CL/F, and female",
    "sex on R-143047 Vc/F. Only creatinine clearance had a clinically",
    "significant effect on active-moiety exposure. The companion file",
    "Kastrissios_2012_managlinatDialanetil_activeMoiety.R holds the",
    "separate, independently fitted model for R-125338 data alone that the",
    "paper develops for comparison. See the vignette Assumptions and",
    "deviations for the variance-versus-SD reading of Table III, the",
    "unprinted age reference value, and two Table III typos.")
  reference <- "Kastrissios H, Walker JR, Carrothers TJ, Kshirsagar S, Khariton T, Habtemariam B, Mager DE, Rohatagi S. Population pharmacokinetic model for a novel oral hypoglycemic formed in vivo: comparing the use of active metabolite data alone versus using data of upstream and downstream metabolites. J Clin Pharmacol. 2012;52(3):404-415. doi:10.1177/0091270010396373"
  vignette <- "Kastrissios_2012_managlinatDialanetil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Kastrissios 2012
  # Figure 1 (metabolism pathways) and Results p. 408.
  compartmentData <- list(
    depot                 = list(analyte = "CS-917 (managlinat dialanetil)", units = "mg", specimen = "administration site", verified = TRUE),
    central               = list(analyte = "CS-917 (managlinat dialanetil)", units = "mg", specimen = "plasma", verified = TRUE),
    central_r134450       = list(analyte = "R-134450 (inactive metabolite)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_r134450   = list(analyte = "R-134450 (inactive metabolite)", units = "mg", specimen = "plasma", verified = TRUE),
    central_r125338       = list(analyte = "R-125338 (active moiety)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_r125338   = list(analyte = "R-125338 (active moiety)", units = "mg", specimen = "plasma", verified = TRUE),
    central_r143047       = list(analyte = "R-143047 (inactive N-acetyl metabolite)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Subject body weight recorded at the screening visit.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CS-917 CL/F and Vc/F and on R-143047 CL/F. Reference 86 kg, printed verbatim in the Kastrissios 2012 Table III footnote equations (`Cl/F = 85.1*(WT/86)^0.95`), consistent with the Table II cohort mean of 87 +/- 16 kg.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age recorded at the screening visit.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on R-134450 CL/F. The reference (median) age is NOT printed anywhere in Kastrissios 2012; 58 years is assumed from the Table II cohort mean of 58 +/- 9 years, consistent with the two medians that ARE printed sitting just below their means (WT 86 vs mean 87; CLCR 70 vs mean 73). See vignette Errata.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance computed with the Cockcroft-Gault equation, RAW and NOT BSA-normalized.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Kastrissios 2012 Methods p. 406: `Creatinine clearance (CLCR) was computed using the Cockcroft-Gault equation`. No BSA normalization is applied, so values are raw mL/min -- the same convention as the registered precedents Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R and Wada_2023_sparsentan.R, and NOT the mL/min/1.73 m^2 canonical default. Power scaling on R-125338 CL/F with reference 70 mL/min, printed verbatim in the Table III footnote (`Cl/F = 30.9*(ClCR/70)^0.36*exp(-0.22*KFood)`); Table II cohort mean 73 +/- 20 mL/min.",
      source_name        = "CLCR"
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on R-143047 Vc/F. Orientation matches the source exactly: Kastrissios 2012 Eq. (4) gives the worked coding `men coded as 0 and women as 1`, so the -0.28 coefficient transfers with no sign flip. Cohort 90 men / 51 women (Table II). NOTE: Table III labels this row `Female sex on CL/F`; the body text states Vc/F three times and exp(-0.28) = 0.756 reproduces its quoted `25% lower apparent volume of distribution`, so the effect is placed on Vc/F -- see vignette Errata.",
      source_name        = "sex"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; White or Hispanic in this cohort)",
      notes              = "Exponential effect on R-134450 CL/F and on R-125338 Vc/F. Source indicator is `KBlack = 1 for black race` per the Kastrissios 2012 Table III footnote. Cohort 87 White / 19 Black / 35 Hispanic (Table II), so Black race is only 13% of subjects; the paper flags the Vc/F effect as remaining to be confirmed in a larger data set.",
      source_name        = "KBlack"
    ),
    FED = list(
      description        = "Fed-state indicator for the dose record: 1 = tablet taken with food, 0 = capsule taken fasted.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted capsule)",
      notes              = "Exponential effect on CS-917 relative bioavailability and on R-125338 CL/F. Source indicator is `KFood = 1 for tablet with food` per the Kastrissios 2012 Table III footnote. Cohort 109 tablet-with-food / 32 fasted-capsule (Table II). The contrast is confounded with formulation, but the authors state that `food status but not formulation was evaluated since capsule and tablet were shown to be bioequivalent (data on file: CTR-917-101)`, so this is a food covariate rather than a formulation covariate. The paper also notes the fed-lower-bioavailability direction conflicts with two phase I crossover studies not included in the analysis.",
      source_name        = "KFood"
    )
  )

  # Screened during covariate selection but not retained in any final moiety
  # model, and reported only as a negative finding with no usable point
  # estimate. Documented here rather than in covariateData so they carry no
  # `declared but not referenced` convention warning.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index at screening (Table II: 31 +/- 5 kg/m^2).",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Significant for R-125338 CL/F and Vc/F on univariate screening (P < .05) and carried into the full model, but dropped in backward elimination; body weight was preferred over BMI for consistency across CL/F and Vc/F."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase at baseline (Table II: 27 +/- 12 IU/L).",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened as a marker of hepatic function. Kastrissios 2012 Results: `There were no effects of liver function tests detected on the pharmacokinetics of any of the 4 moieties`, with the caveat that all patients were within the normal liver-function range."
    ),
    AST = list(
      description = "Serum aspartate aminotransferase at baseline (Table II: 21 +/- 7 IU/L).",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened as a marker of hepatic function; not retained. See ALT."
    ),
    TBILI = list(
      description = "Total bilirubin at baseline (Table II: 0.69 +/- 0.26 mg/dL).",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Significant for CS-917 CL/F on univariate screening and carried into the full model, but dropped in backward elimination. See ALT for the liver-function caveat."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 141L,
    n_studies      = 6L,
    n_observations = 8961L,
    age_range      = "mean 58 +/- 9 years (range not reported)",
    weight_range   = "mean 87 +/- 16 kg (range not reported)",
    sex_female_pct = 36.2,
    race_ethnicity = c(White = 61.7, Black = 13.5, Hispanic = 24.8),
    disease_state  = "Type 2 diabetes mellitus; entry required hemoglobin A1c >= 6.5% and fasting plasma glucose 160-300 mg/dL at baseline. Renal and hepatic function were primarily normal with little variation within the cohort (creatinine clearance 73 +/- 20 mL/min).",
    dose_range     = "Multiple oral CS-917 50-400 mg, QD / BID / TID, over 10-28 days (or 14 days per period in the two crossover studies). All subjects received multiple doses; all but study C were dosed with food.",
    notes          = "Four phase I studies (A, B, E, F) and two phase IIa studies (C, D) per Table I. Only monotherapy observations were used for structural model development; combination-therapy records and records that were zero, below the limit of detection, or missing were excluded, as were 14 observations with absolute weighted residuals greater than 6. Assay range 5-2000 ng/mL for CS-917, R-134450 and R-125338 and 10-2000 ng/mL for R-143047; interbatch accuracy and precision were < 11.5% and < 12% for all analytes. Fitted in NONMEM 5.1.1 with the first-order (FO) method."
  )

  ini({
    # ------------------------------------------------------------------
    # CS-917 (prodrug, parent analyte -- bare canonical names)
    # All values: Kastrissios 2012 Table III, "Linked Model" / CS-917
    # column. The parenthesised number in that table is the RSE of the
    # estimate ("%CV, coefficient of variation of the estimate (absolute
    # value of 100*SE/estimate)"), NOT a between-subject CV.
    # ------------------------------------------------------------------
    lka <- log(5.40)
    label("CS-917 first-order absorption rate constant (ka, 1/h)")  # Table III CS-917 ka = 5.40 (RSE 21%); Table III prints the unit as "L/h", a typo -- see vignette Errata
    ltlag <- log(0.23)
    label("CS-917 absorption lag time (Tlag, h)")  # Table III CS-917 Tlag = 0.23 h (RSE 3%)
    lcl <- log(85.1)
    label("CS-917 apparent clearance (CL/F, L/h)")  # Table III CS-917 CL/F = 85.1 (RSE 8%); footnote equation Cl/F = 85.1*(WT/86)^0.95
    lvc <- log(68.0)
    label("CS-917 apparent central volume (Vc/F, L)")  # Table III CS-917 Vc/F = 68.0 (RSE 9%); footnote equation Vc/F = 68*(WT/86)^0.63
    lfdepot <- fixed(log(1))
    label("CS-917 relative oral bioavailability for the fasted-capsule reference (unitless)")  # Table III footnote anchors the reference at F = 1: "F = 1*exp(-0.27*KFood)"

    e_wt_cl <- 0.95
    label("Power exponent on (WT/86) for CS-917 CL/F (unitless)")  # Table III Effect, body weight on CL/F = 0.95 (RSE 24%)
    e_wt_vc <- 0.63
    label("Power exponent on (WT/86) for CS-917 Vc/F (unitless)")  # Table III Effect, body weight on Vc/F = 0.63 (RSE 45%)
    e_fed_fdepot <- -0.27
    label("Exponential coefficient on FED for CS-917 relative bioavailability (unitless)")  # Table III Effect, food on relative F = -0.27 (RSE 40%); exp(-0.27) = 0.763, matching the paper's "24% lower relative bioavailability"

    # IIV. Kastrissios 2012 Eq. (1) is theta_i = theta_typical*exp(eta_i)
    # with eta of "mean 0 and variance omega^2", and the source is NONMEM
    # 5.1.1 output, whose $OMEGA is reported as a variance. The Table III
    # "Interindividual variability" entries are therefore taken as
    # variances and used directly on nlmixr2's variance scale. See the
    # vignette Assumptions and deviations for the full four-way argument.
    etalka ~ 0.94    # Table III CS-917 IIV Ka = 0.94 (RSE 45%)
    etaltlag ~ 0.01  # Table III CS-917 IIV Tlag = 0.01 (RSE 95%)
    etalcl ~ 0.20    # Table III CS-917 IIV CL/F = 0.20 (RSE 21%)
    etalvc ~ 0.16    # Table III CS-917 IIV Vc/F = 0.16 (RSE 36%)

    # ------------------------------------------------------------------
    # R-134450 (first, inactive metabolite)
    # ------------------------------------------------------------------
    lcl_r134450 <- log(187)
    label("R-134450 apparent clearance (CL/F, L/h)")  # Table III R-134450 CL/F = 187 (RSE 15%)
    lvc_r134450 <- log(50.7)
    label("R-134450 apparent central volume (Vc/F, L)")  # Table III R-134450 Vc/F = 50.7 (RSE 30%)
    lq_r134450 <- log(197)
    label("R-134450 apparent intercompartmental clearance (Q/F, L/h)")  # Table III R-134450 Q/F = 197 (RSE 11%)
    lvp_r134450 <- log(249)
    label("R-134450 apparent peripheral volume (Vp/F, L)")  # Table III R-134450 Vp/F = 249 (RSE 9%)

    e_age_cl_r134450 <- -0.55
    label("Power exponent on (AGE/58) for R-134450 CL/F (unitless)")  # Table III Effect, age on CL/F = -0.55 (RSE 45%); reference age assumed, not printed -- see vignette Errata
    e_race_black_cl_r134450 <- -0.37
    label("Exponential coefficient on RACE_BLACK for R-134450 CL/F (unitless)")  # Table III Effect, black race on CL/F = -0.37 (RSE 46%); exp(-0.37) = 0.691, matching the paper's "31% lower apparent clearance"

    etalcl_r134450 ~ 0.95  # Table III R-134450 IIV CL/F = 0.95 (RSE 57%)
    etalvc_r134450 ~ 1.37  # Table III R-134450 IIV Vc/F = 1.37 (RSE 49%)
    etalq_r134450 ~ 0.39   # Table III R-134450 IIV Q/F = 0.39 (RSE 63%)
    etalvp_r134450 ~ 0.60  # Table III R-134450 IIV Vp/F = 0.60 (RSE 30%)

    # ------------------------------------------------------------------
    # R-125338 (active moiety)
    # ------------------------------------------------------------------
    lcl_r125338 <- log(30.9)
    label("R-125338 apparent clearance (CL/F, L/h)")  # Table III R-125338 CL/F = 30.9 (RSE 7%); footnote equation Cl/F = 30.9*(ClCR/70)^0.36*exp(-0.22*KFood)
    lvc_r125338 <- log(60.7)
    label("R-125338 apparent central volume (Vc/F, L)")  # Table III R-125338 Vc/F = 60.7 (RSE 10%); footnote equation Vc/F = 60.7*exp(-0.56*KBlack)
    lq_r125338 <- log(7.07)
    label("R-125338 apparent intercompartmental clearance (Q/F, L/h)")  # Table III R-125338 Q/F = 7.07 (RSE 14%)
    lvp_r125338 <- log(278)
    label("R-125338 apparent peripheral volume (Vp/F, L)")  # Table III R-125338 Vp/F = 278 (RSE 7%)

    e_crcl_cl_r125338 <- 0.36
    label("Power exponent on (CRCL/70) for R-125338 CL/F (unitless)")  # Table III Effect, CLCR on CL/F = 0.36 (RSE 27%); footnote equation
    e_fed_cl_r125338 <- -0.22
    label("Exponential coefficient on FED for R-125338 CL/F (unitless)")  # Table III Effect, food on CL/F = -0.22 (RSE 37%); exp(-0.22) = 0.803, matching the paper's "20% lower apparent clearance"
    e_race_black_vc_r125338 <- -0.56
    label("Exponential coefficient on RACE_BLACK for R-125338 Vc/F (unitless)")  # Table III Effect, black race on Vc/F = -0.56 (RSE 30%); exp(-0.56) = 0.571, matching the paper's "43% reduction in apparent volume of distribution"

    etalcl_r125338 ~ 0.14  # Table III R-125338 IIV CL/F = 0.14 (RSE 19%)
    etalvc_r125338 ~ 0.72  # Table III R-125338 IIV Vc/F = 0.72 (RSE 22%)
    etalq_r125338 ~ 0.44   # Table III R-125338 IIV Q/F = 0.44 (RSE 23%)
    etalvp_r125338 ~ 0.19  # Table III R-125338 IIV Vp/F = 0.19 (RSE 23%)

    # ------------------------------------------------------------------
    # R-143047 (terminal, inactive N-acetyl metabolite)
    # ------------------------------------------------------------------
    lcl_r143047 <- log(9.81)
    label("R-143047 apparent clearance (CL/F, L/h)")  # Table III R-143047 CL/F = 9.81 (RSE 7%)
    lvc_r143047 <- log(338)
    label("R-143047 apparent central volume (Vc/F, L)")  # Table III R-143047 Vc/F = 338 (RSE 8%)

    e_wt_cl_r143047 <- 0.61
    label("Power exponent on (WT/86) for R-143047 CL/F (unitless)")  # Table III Effect, body weight on CL/F = 0.61 (RSE 34%)
    e_sexf_vc_r143047 <- -0.28
    label("Exponential coefficient on SEXF for R-143047 Vc/F (unitless)")  # Table III Effect = -0.28 (RSE 45%), row labelled "Female sex on CL/F"; placed on Vc/F per the body text (stated 3x) and exp(-0.28) = 0.756 matching "25% lower apparent volume of distribution" -- see vignette Errata

    etalcl_r143047 ~ 0.36  # Table III R-143047 IIV CL/F = 0.36 (RSE 27%)
    etalvc_r143047 ~ 0.25  # Table III R-143047 IIV Vc/F = 0.25 (RSE 30%)

    # ------------------------------------------------------------------
    # Residual error. Kastrissios 2012 Eq. (2) is
    #   y_ij = yhat_ij*(1 + eps1_ij) + eps2_ij
    # i.e. combined proportional + additive in linear concentration
    # space, one pair per moiety ("Residual error was modeled with
    # proportional and additive terms for each of the 4 moieties").
    # The Table III entries are NONMEM $SIGMA variances, so each
    # nlmixr2 SD is sqrt(table value). Decisive check: R-143047's
    # proportional 0.02 read as an SD would be a 2% residual error,
    # below the assay's own reported < 12% interbatch precision, whereas
    # sqrt(0.02) = 14% sits right at it; and the additive terms are only
    # dimensionally sensible as variances (sqrt(1620) = 40.2 ng/mL for
    # CS-917 against a 5 ng/mL assay lower limit, versus an impossible
    # 1620 ng/mL). Additive SDs are in ng/mL to match `units`.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.28)
    label("Proportional residual SD for CS-917 (unitless)")  # Table III CS-917 proportional = 0.28 (RSE 14%), a variance -> SD = sqrt(0.28) = 0.529
    addSd <- sqrt(1620)
    label("Additive residual SD for CS-917 (ng/mL)")  # Table III CS-917 additive = 1620 (RSE 38%), a variance -> SD = sqrt(1620) = 40.2 ng/mL

    propSd_r134450 <- sqrt(0.15)
    label("Proportional residual SD for R-134450 (unitless)")  # Table III R-134450 proportional = 0.15 (RSE 28%) -> SD = 0.387
    addSd_r134450 <- sqrt(146)
    label("Additive residual SD for R-134450 (ng/mL)")  # Table III R-134450 additive = 146 (RSE 20%) -> SD = 12.1 ng/mL

    propSd_r125338 <- sqrt(0.07)
    label("Proportional residual SD for R-125338 (unitless)")  # Table III R-125338 proportional = 0.07 (RSE 12%) -> SD = 0.265
    addSd_r125338 <- sqrt(13.7)
    label("Additive residual SD for R-125338 (ng/mL)")  # Table III R-125338 additive = 13.7 (RSE 38%) -> SD = 3.70 ng/mL

    propSd_r143047 <- sqrt(0.02)
    label("Proportional residual SD for R-143047 (unitless)")  # Table III R-143047 proportional = 0.02 (RSE 20%) -> SD = 0.141
    addSd_r143047 <- sqrt(271)
    label("Additive residual SD for R-143047 (ng/mL)")  # Table III R-143047 additive = 271 (RSE 36%) -> SD = 16.5 ng/mL
  })

  model({
    # Covariate reference values. The two printed verbatim in the
    # Kastrissios 2012 Table III footnote equations are exact; the age
    # reference is assumed from the Table II cohort mean because the
    # paper prints no median age and no R-134450 covariate equation.
    wtRef   <- 86  # kg,     Table III footnote: (WT/86)
    crclRef <- 70  # mL/min, Table III footnote: (ClCR/70)
    ageRef  <- 58  # years,  assumed (Table II mean 58 +/- 9); not printed

    # Doses are in mg and volumes in L, so amount/volume is mg/L. The
    # source reports concentrations in ng/mL (assay range 5-2000 ng/mL,
    # Table IV Cmax in ng/mL), and 1 mg/L = 1000 ng/mL.
    mgPerLToNgPerMl <- 1000

    # ------------------------------------------------------------------
    # Individual parameters. Continuous covariates enter as a power of
    # the median-normalized value (Kastrissios 2012 Eq. 3) and
    # categorical covariates as exp(indicator * coefficient) (Eq. 4).
    # ------------------------------------------------------------------
    # CS-917 (prodrug)
    ka     <- exp(lka + etalka)
    tlag   <- exp(ltlag + etaltlag)
    cl     <- exp(lcl + etalcl) * (WT / wtRef)^e_wt_cl
    vc     <- exp(lvc + etalvc) * (WT / wtRef)^e_wt_vc
    fdepot <- exp(lfdepot + e_fed_fdepot * FED)

    # R-134450
    cl_r134450 <- exp(lcl_r134450 + etalcl_r134450) *
      (AGE / ageRef)^e_age_cl_r134450 *
      exp(e_race_black_cl_r134450 * RACE_BLACK)
    vc_r134450 <- exp(lvc_r134450 + etalvc_r134450)
    q_r134450  <- exp(lq_r134450 + etalq_r134450)
    vp_r134450 <- exp(lvp_r134450 + etalvp_r134450)

    # R-125338 (active moiety)
    cl_r125338 <- exp(lcl_r125338 + etalcl_r125338) *
      (CRCL / crclRef)^e_crcl_cl_r125338 *
      exp(e_fed_cl_r125338 * FED)
    vc_r125338 <- exp(lvc_r125338 + etalvc_r125338) *
      exp(e_race_black_vc_r125338 * RACE_BLACK)
    q_r125338  <- exp(lq_r125338 + etalq_r125338)
    vp_r125338 <- exp(lvp_r125338 + etalvp_r125338)

    # R-143047
    cl_r143047 <- exp(lcl_r143047 + etalcl_r143047) *
      (WT / wtRef)^e_wt_cl_r143047
    vc_r143047 <- exp(lvc_r143047 + etalvc_r143047) *
      exp(e_sexf_vc_r143047 * SEXF)

    # Micro-constants
    kel          <- cl / vc
    kel_r134450  <- cl_r134450 / vc_r134450
    k12_r134450  <- q_r134450 / vc_r134450
    k21_r134450  <- q_r134450 / vp_r134450
    kel_r125338  <- cl_r125338 / vc_r125338
    k12_r125338  <- q_r125338 / vc_r125338
    k21_r125338  <- q_r125338 / vp_r125338
    kel_r143047  <- cl_r143047 / vc_r143047

    # ------------------------------------------------------------------
    # ODE system (Kastrissios 2012 Figure 1). Each metabolite is formed
    # by first-order transfer of the WHOLE upstream central-compartment
    # elimination flux, with no fractional-conversion multiplier: the
    # source states the metabolite CL/F, Q/F, Vc/F and Vp/F terms "were
    # expressed relative to the oral bioavailability of CS-917 and the
    # fraction of CS-917 metabolized to moiety of interest, neither of
    # which were known", so the unknown fraction is already absorbed
    # into the downstream apparent volume and clearance. Writing the
    # apparent amounts A* = A/(F*fm) makes this exact: the fm cancels
    # and the transfer becomes kel_upstream * A*_upstream. Distribution
    # flux to a peripheral compartment is NOT part of the transfer.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    d/dt(central_r134450) <- kel * central -
      kel_r134450 * central_r134450 -
      k12_r134450 * central_r134450 + k21_r134450 * peripheral1_r134450
    d/dt(peripheral1_r134450) <- k12_r134450 * central_r134450 -
      k21_r134450 * peripheral1_r134450

    d/dt(central_r125338) <- kel_r134450 * central_r134450 -
      kel_r125338 * central_r125338 -
      k12_r125338 * central_r125338 + k21_r125338 * peripheral1_r125338
    d/dt(peripheral1_r125338) <- k12_r125338 * central_r125338 -
      k21_r125338 * peripheral1_r125338

    d/dt(central_r143047) <- kel_r125338 * central_r125338 -
      kel_r143047 * central_r143047

    alag(depot) <- tlag
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # Observations, one per measured moiety.
    # ------------------------------------------------------------------
    Cc         <- mgPerLToNgPerMl * central / vc
    Cc_r134450 <- mgPerLToNgPerMl * central_r134450 / vc_r134450
    Cc_r125338 <- mgPerLToNgPerMl * central_r125338 / vc_r125338
    Cc_r143047 <- mgPerLToNgPerMl * central_r143047 / vc_r143047

    Cc         ~ add(addSd) + prop(propSd)
    Cc_r134450 ~ add(addSd_r134450) + prop(propSd_r134450)
    Cc_r125338 ~ add(addSd_r125338) + prop(propSd_r125338)
    Cc_r143047 ~ add(addSd_r143047) + prop(propSd_r143047)
  })
}
