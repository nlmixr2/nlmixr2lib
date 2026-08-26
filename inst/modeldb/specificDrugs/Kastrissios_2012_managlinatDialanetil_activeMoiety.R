Kastrissios_2012_managlinatDialanetil_activeMoiety <- function() {
  description <- paste(
    "Reduced single-analyte population pharmacokinetic model for R-125338,",
    "the active moiety formed in vivo from the oral hypoglycemic",
    "fructose-1,6-bisphosphatase-inhibitor prodrug CS-917 (INN managlinat",
    "dialanetil), in 141 patients with type 2 diabetes mellitus pooled from",
    "six phase I / IIa studies (Kastrissios 2012). This is the paper's",
    "second, independently fitted model: only R-125338 observations were",
    "analyzed, so R-125338 is the parent analyte here and takes the bare",
    "canonical compartment and parameter names. Two-compartment disposition",
    "with first-order absorption and an absorption lag applied to the oral",
    "CS-917 dose; the apparent bioavailability folded into CL/F and V/F",
    "reflects formation of R-125338 from the bioavailable fraction of the",
    "administered prodrug, so no separate bioavailability parameter is",
    "estimated. Retained covariates: creatinine clearance on CL/F and food",
    "on the absorption rate constant. The paper develops this model to test",
    "whether fitting the active moiety alone reproduces the conclusions of",
    "the full linked cascade, and finds that it does: individual R-125338",
    "predictions from the two models correlate at 0.989, and halving",
    "creatinine clearance raises active-moiety exposure by 31% here versus",
    "27% under the linked model, while requiring substantially less",
    "analysis time. The companion file",
    "Kastrissios_2012_managlinatDialanetil_linked.R holds the full",
    "four-moiety cascade. See the vignette Assumptions and deviations for",
    "the variance-versus-SD reading of Table III and the covariate forms",
    "reconstructed from Eqs. (3) and (4).")
  reference <- "Kastrissios H, Walker JR, Carrothers TJ, Kshirsagar S, Khariton T, Habtemariam B, Mager DE, Rohatagi S. Population pharmacokinetic model for a novel oral hypoglycemic formed in vivo: comparing the use of active metabolite data alone versus using data of upstream and downstream metabolites. J Clin Pharmacol. 2012;52(3):404-415. doi:10.1177/0091270010396373"
  vignette <- "Kastrissios_2012_managlinatDialanetil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Kastrissios 2012
  # "Model for R-125338" (p. 407) and Table III "Single Moiety" column.
  # The dose administered is CS-917, but every amount in this model is an
  # APPARENT R-125338 amount -- the prodrug-to-active-moiety conversion is
  # absorbed into the apparent absorption and disposition parameters.
  compartmentData <- list(
    depot       = list(analyte = "R-125338 (apparent, dosed as CS-917 prodrug)", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "R-125338 (active moiety)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "R-125338 (active moiety)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance computed with the Cockcroft-Gault equation, RAW and NOT BSA-normalized.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Kastrissios 2012 Methods p. 406: `Creatinine clearance (CLCR) was computed using the Cockcroft-Gault equation`. No BSA normalization, so values are raw mL/min -- the same convention as the registered precedents Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R and Wada_2023_sparsentan.R, and NOT the mL/min/1.73 m^2 canonical default. Power scaling on CL/F with exponent 0.37 (Table III, Single Moiety column). The paper prints the covariate equation only for the LINKED model, whose R-125338 arm is centered at 70 mL/min; the same reference is used here per Eq. (3) (median normalization) since the cohort is identical. Table II cohort mean 73 +/- 20 mL/min.",
      source_name        = "CLCR"
    ),
    FED = list(
      description        = "Fed-state indicator for the dose record: 1 = tablet taken with food, 0 = capsule taken fasted.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted capsule)",
      notes              = "Exponential effect on the absorption rate constant, exp(-1.07 * FED), i.e. food slows apparent absorption of the active moiety about three-fold (2.54 -> 0.87 1/h). Source indicator is `KFood = 1 for tablet with food` per the Kastrissios 2012 Table III footnote. Cohort 109 tablet-with-food / 32 fasted-capsule (Table II). Note the food effect lands on a different parameter here than in the linked model, where it acts on CS-917 relative bioavailability and on R-125338 CL/F; the paper reports that food had a negligible impact on exposure under this reduced model. The contrast is confounded with formulation, but the authors state capsule and tablet were shown to be bioequivalent (data on file: CTR-917-101).",
      source_name        = "KFood"
    )
  )

  # Screened during covariate selection for this cohort but not retained in
  # the reduced single-moiety model. Documented here rather than in
  # covariateData so they carry no `declared but not referenced` warning.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at screening (Table II: 87 +/- 16 kg).",
      units       = "kg",
      type        = "continuous",
      notes       = "Retained on CS-917 and R-143047 clearance in the linked model, but not retained for R-125338 in either model. Body mass index was preferred over weight when the R-125338 full model was built, and neither survived backward elimination."
    ),
    RACE_BLACK = list(
      description = "Black race indicator (Table II: 19 of 141 patients).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Retained on R-125338 apparent central volume in the LINKED model (coefficient -0.56) but absent from the Table III Single Moiety column, so it is not part of this reduced model."
    ),
    BMI = list(
      description = "Body mass index at screening (Table II: 31 +/- 5 kg/m^2).",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Significant for R-125338 CL/F and Vc/F on univariate screening (P < .05) and carried into the full model, but dropped in backward elimination."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase at baseline (Table II: 27 +/- 12 IU/L).",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Kastrissios 2012 Results: `There were no effects of liver function tests detected on the pharmacokinetics of any of the 4 moieties`, with the caveat that all patients were within the normal liver-function range."
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
      notes       = "Screened as a marker of hepatic function; not retained. See ALT."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 141L,
    n_studies      = 6L,
    age_range      = "mean 58 +/- 9 years (range not reported)",
    weight_range   = "mean 87 +/- 16 kg (range not reported)",
    sex_female_pct = 36.2,
    race_ethnicity = c(White = 61.7, Black = 13.5, Hispanic = 24.8),
    disease_state  = "Type 2 diabetes mellitus; entry required hemoglobin A1c >= 6.5% and fasting plasma glucose 160-300 mg/dL at baseline. Renal and hepatic function were primarily normal with little variation within the cohort (creatinine clearance 73 +/- 20 mL/min).",
    dose_range     = "Multiple oral CS-917 50-400 mg, QD / BID / TID, over 10-28 days (or 14 days per period in the two crossover studies).",
    notes          = "Same six-study cohort as the linked model (four phase I: A, B, E, F; two phase IIa: C, D), but only the R-125338 observations were analyzed in this fit, so the observation count is a subset of the linked model's 8961. R-125338 assay range 5-2000 ng/mL. Fitted in NONMEM 5.1.1 with the first-order (FO) method."
  )

  ini({
    # ------------------------------------------------------------------
    # All values: Kastrissios 2012 Table III, "Single Moiety" / R-125338
    # column. The parenthesised number in that table is the RSE of the
    # estimate, NOT a between-subject CV.
    # ------------------------------------------------------------------
    lka <- log(2.54)
    label("R-125338 apparent first-order absorption rate constant (ka, 1/h)")  # Table III Single Moiety ka = 2.54 (RSE 17%); Table III prints the unit as "L/h", a typo -- see vignette Errata
    ltlag <- log(0.46)
    label("R-125338 apparent absorption lag time (Tlag, h)")  # Table III Single Moiety Tlag = 0.46 h (RSE 1%)
    lcl <- log(35.3)
    label("R-125338 apparent clearance (CL/F, L/h)")  # Table III Single Moiety CL/F = 35.3 (RSE 4%)
    lvc <- log(187)
    label("R-125338 apparent central volume (Vc/F, L)")  # Table III Single Moiety Vc/F = 187 (RSE 6%)
    lq <- log(5.61)
    label("R-125338 apparent intercompartmental clearance (Q/F, L/h)")  # Table III Single Moiety Q/F = 5.61 (RSE 11%)
    lvp <- log(313)
    label("R-125338 apparent peripheral volume (Vp/F, L)")  # Table III Single Moiety Vp/F = 313 (RSE 6%)

    e_crcl_cl <- 0.37
    label("Power exponent on (CRCL/70) for R-125338 CL/F (unitless)")  # Table III Effect, CLCR on CL/F = 0.37 (RSE 23%)
    e_fed_ka <- -1.07
    label("Exponential coefficient on FED for R-125338 ka (unitless)")  # Table III Effect, food on ka = -1.07 (RSE 25%); exp(-1.07) = 0.343

    # IIV. Kastrissios 2012 Eq. (1) is theta_i = theta_typical*exp(eta_i)
    # with eta of "mean 0 and variance omega^2", and the source is NONMEM
    # 5.1.1 output, whose $OMEGA is reported as a variance, so the
    # Table III "Interindividual variability" entries are used directly
    # on nlmixr2's variance scale. Table III reports no Tlag IIV for the
    # single-moiety model (the cell is a dash), so none is encoded.
    etalka ~ 0.52  # Table III Single Moiety IIV Ka = 0.52 (RSE 26%)
    etalcl ~ 0.12  # Table III Single Moiety IIV CL/F = 0.12 (RSE 14%)
    etalvc ~ 0.19  # Table III Single Moiety IIV Vc/F = 0.19 (RSE 18%)
    etalq ~ 0.34   # Table III Single Moiety IIV Q/F = 0.34 (RSE 20%)
    etalvp ~ 0.17  # Table III Single Moiety IIV Vp/F = 0.17 (RSE 25%)

    # Residual error. Kastrissios 2012 Eq. (2),
    # y_ij = yhat_ij*(1 + eps1_ij) + eps2_ij, i.e. combined proportional
    # + additive in linear concentration space. The Table III entries are
    # NONMEM $SIGMA variances, so each nlmixr2 SD is sqrt(table value);
    # sqrt(9.63) = 3.10 ng/mL sits just below the 5 ng/mL assay lower
    # limit, as an additive residual should, whereas 9.63 ng/mL read as
    # an SD would be twice it. Additive SD is in ng/mL to match `units`.
    propSd <- sqrt(0.08)
    label("Proportional residual SD for R-125338 (unitless)")  # Table III Single Moiety proportional = 0.08 (RSE 11%), a variance -> SD = 0.283
    addSd <- sqrt(9.63)
    label("Additive residual SD for R-125338 (ng/mL)")  # Table III Single Moiety additive = 9.63 (RSE 40%), a variance -> SD = 3.10 ng/mL
  })

  model({
    # Covariate reference value. The paper prints the median-normalized
    # covariate equation only for the linked model, whose R-125338 arm is
    # centered at 70 mL/min (Table III footnote); the same cohort median
    # applies here per Eq. (3).
    crclRef <- 70  # mL/min

    # Doses are in mg and volumes in L, so amount/volume is mg/L, while
    # the source reports concentrations in ng/mL. 1 mg/L = 1000 ng/mL.
    mgPerLToNgPerMl <- 1000

    # Individual parameters. Continuous covariates enter as a power of
    # the median-normalized value (Kastrissios 2012 Eq. 3) and
    # categorical covariates as exp(indicator * coefficient) (Eq. 4).
    ka   <- exp(lka + etalka) * exp(e_fed_ka * FED)
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl) * (CRCL / crclRef)^e_crcl_cl
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order absorption from the
    # oral CS-917 dose. No bioavailability parameter is estimated: per
    # Kastrissios 2012 p. 407, "the oral bioavailability parameter
    # reflected formation of R-125338 from the bioavailable fraction of
    # the orally administered dose of prodrug CS-917", and that
    # unidentifiable factor is absorbed into the apparent CL/F and V/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    Cc <- mgPerLToNgPerMl * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
