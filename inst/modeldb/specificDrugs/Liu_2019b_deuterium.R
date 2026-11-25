Liu_2019b_deuterium <- function() {
  description <- "Mother-infant dyad model for deuterium oxide (D2O) given orally to the lactating mother (the dose-to-the-mother technique for classifying exclusive breastfeeding), refitted to 565 mother-infant pairs pooled from nine countries. Two linked one-compartment models: the mother's body water receives the oral dose and loses D2O by a first-order total elimination rate constant, part of which is the lactation clearance into the breastfed infant; the infant's body water receives that flow and loses D2O by its own water clearance. The infant volume of distribution is not estimated but computed from infant body weight (total body water x 1.041) and grows with the infant. Covariates: maternal weight on the maternal volume, infant weight on the infant clearance. Fitted in a fully Bayesian framework (Stan) to maternal and infant saliva enrichment, with a country-specific combined residual error."
  reference <- paste(
    "Liu Z, Diana A, Slater C, Preston T, Gibson RS, Houghton L, Duffull SB.",
    "Development of a parsimonious design for optimal classification of exclusive breastfeeding.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(8):596-605. doi:10.1002/psp4.12428.",
    "Structure and Stan code from Liu et al., J Pharmacokinet Pharmacodyn 2019;46(1):1",
    "(doi:10.1007/s10928-018-9613-x); parameter values from Eqs. 6-7, Table 3 and",
    "Supplementary Material S3 (Table S2) of the 2019 CPT:PSP article.",
    sep = " "
  )
  vignette <- "Liu_2019b_deuterium"
  units <- list(time = "day", dosing = "g", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "deuterium oxide", units = "g", specimen = "saliva", verified = TRUE),
    infant_central = list(analyte = "deuterium oxide", units = "g", specimen = "saliva", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Maternal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Power effect on the maternal D2O volume of distribution, (WT / 70)^0.46 (Eq. 6; upstream Stan code 'pow(MWT[i]/70, theta[1])'). The paper calls 70 kg the median maternal weight, but Table 1 gives the model-building median as 56 kg (IQR 50-64); 70 kg is the value the equation divides by. Source column MWT.",
      source_name = "MWT"
    ),
    WT_INFANT = list(
      description = "Breastfed infant body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Enters twice. (1) The infant D2O volume of distribution Vb = 1.041 * exp(-0.427 + 0.963 * log(WT_INFANT)) (upstream Liu 2019 JPKPD Eqs. A1.9-A1.10 and Stan code, which interpolates infant weight linearly between the start-of-study and end-of-study weighings), so supply WT_INFANT at every record (or use linear covariate interpolation) to reproduce the growth of Vb. (2) A power effect on the infant water clearance, (WT_INFANT / 5)^0.55 (Eq. 7); the upstream Stan code evaluates this term at the infant's AVERAGE weight over the study (Stan 'BWT'), whereas here it follows the current weight. Table 1 (model-building data): start weight 6.2 kg (IQR 5.5-7.2), end weight 6.5 kg (IQR 5.8-7.4). Source columns Wt_s / Wt_e (volume) and BWT (clearance).",
      source_name = "BWT"
    )
  )

  covariatesDataExcluded <- list(
    COUNTRY = list(
      description = "Country of the field study (nine countries; one assay method per country)",
      units = "(categorical)",
      type = "categorical",
      notes = "Covariate on the combined residual error model only (Results, Model building and selection; Supplementary Material S4 Table S3 model 5), accounting for regional and assay (IRMS vs FTIR) differences. The country-specific sigma values are not reported, so the residual error is held at zero and country is not used."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 565,
    n_studies = 9,
    age_range = "Mothers IQR 21-31 years; infants IQR 2.9-5.4 months",
    age_median = "Mothers 26 years; infants 3.4 months",
    weight_range = "Mothers IQR 50-64 kg; infants IQR 5.5-7.2 kg at study start",
    weight_median = "Mothers 56 kg; infants 6.2 kg at study start and 6.5 kg at study end",
    sex_female_pct = NA_real_,
    disease_state = "Healthy lactating mothers and their breastfed infants: one calibration study of exclusively breastfeeding pairs (Indonesia, 113 pairs) and two-thirds of a composite of field studies whose breastfeeding status was not controlled.",
    dose_range = "Single oral dose of deuterium oxide to the mother, median 30 g (range 6-60 g); doses below 30 g in the IRMS studies and 30 g or more in the FTIR studies.",
    regions = "Indonesia, South Africa, Thailand, Kenya, Sri Lanka, India, Guatemala, Ghana, Chile.",
    infant_partner = "Infant sex 36% female, 36% male, 27% unknown (Table 1, model-building column).",
    notes = "Model-building data set: the whole calibration study (113 pairs, Indonesia) plus two-thirds of the composite field study, 565 pairs and 6485 observations, after excluding 2 pairs as outliers (infant water clearance > 40% of body weight). The remaining third of the field study (225 pairs, 2481 observations) was held out for design evaluation and was not used to fit the model. Saliva D2O was measured by FTIR in 71% and IRMS in 29% of pairs (Table 1). Mothers are female; sex_female_pct is left NA because the modelled dyad includes the infant."
  )

  ini({
    # Typical values are population means of the log-scale individual
    # parameters (Stan 'log_theta'), at maternal weight 70 kg and infant weight
    # 5 kg. Eqs. 6-7 print the final-model values for Vm and CLbo. The paper
    # does not print final-model values for kmm and CLmb (no covariates act on
    # them); their base-model posterior on the same data (Table S2) is used.
    # Supplementary Material S5 lists each pair's final-model posterior mean of
    # the four log parameters before the covariate effects; their averages are
    # 3.563, -2.168, -0.357 and -0.158, matching all four values below.
    lvc <- 3.54
    label("Maternal D2O volume of distribution Vm at maternal WT = 70 kg (L)") # Eq. 6 'ln(Vm,i) = N(3.54, 0.12) + N(0.46, 0.03) ln(MWTi/70kg)'
    lkel <- -2.17
    label("Maternal total D2O elimination rate constant kmm, including the lactation route (1/day)") # Supplementary Material S3 Table S2 'ln(kmm) mean -2.17' (base-model posterior); main-text Table 3 prints -2.18
    lcl_milk <- -0.36
    label("Lactation-route water clearance from mother to infant CLmb (L/day)") # Supplementary Material S3 Table S2 'ln(CLmb) mean -0.36' (base-model posterior)
    lcl_infant <- -0.16
    label("Infant water clearance to the outside CLbo at infant WT = 5 kg (L/day)") # Eq. 7 'ln(CLbo,i) = N(-0.16, 0.17) + N(0.55, 0.03) ln(BWTi/5kg)'

    e_wt_vc <- 0.46
    label("Power exponent of maternal weight on Vm (unitless)") # Eq. 6 'N(0.46, 0.03)'
    e_wt_infant_cl_infant <- 0.55
    label("Power exponent of infant weight on CLbo (unitless)") # Eq. 7 'N(0.55, 0.03)'

    # Between-subject variability: log-normal. The second argument of N() in
    # Eqs. 6-7 on the intercept is the between-subject SD on the log scale (the
    # convention of the upstream paper, whose Eq. 16 SD equals its tabulated
    # BSV). The kmm and CLmb SDs are the base-model 'BSV mean' of Table S2.
    # The fitted LKJ correlations are not reported and are set to zero.
    etalvc ~ 0.0144 # Eq. 6 intercept SD 0.12: 0.12^2
    etalkel ~ 0.0324 # Table S2 'ln(kmm) BSV mean 0.18': 0.18^2
    etalcl_milk ~ 0.0961 # Table S2 'ln(CLmb) BSV mean 0.31': 0.31^2
    etalcl_infant ~ 0.0289 # Eq. 7 intercept SD 0.17: 0.17^2

    # Combined residual error with country-specific sigma (Results, Model
    # building and selection). The sigma values are not reported anywhere in
    # the paper or its supplements, so both pairs are held at zero; rxode2
    # needs a separate pair per endpoint.
    propSd <- fixed(0)
    label("Proportional residual SD, maternal saliva (fraction; country-specific in the source; not reported, held at zero)") # Eq. 2; value not reported
    addSd <- fixed(0)
    label("Additive residual SD, maternal saliva (mg/L; country-specific in the source; not reported, held at zero)") # Eq. 2; value not reported
    propSd_Cinfant <- fixed(0)
    label("Proportional residual SD, infant saliva (fraction; country-specific in the source; not reported, held at zero)") # Eq. 2; value not reported
    addSd_Cinfant <- fixed(0)
    label("Additive residual SD, infant saliva (mg/L; country-specific in the source; not reported, held at zero)") # Eq. 2; value not reported
  })

  model({
    # Infant D2O volume of distribution from infant weight: total body water
    # ln(TBW) = a + b ln(WT) with a ~ N(-0.427, 0.012), b ~ N(0.963, 0.005)
    # held at their means (Supplementary Material S1 Eq. 17), enlarged by 4.1%
    # for non-aqueous isotope exchange (upstream Liu 2019 JPKPD Eq. A1.10 and
    # Stan 'Vb_b = 1.041*exp(a + b*log(Wt))'; the factor is not restated here).
    tbw_infant <- exp(-0.427 + 0.963 * log(WT_INFANT))
    vc_infant <- 1.041 * tbw_infant

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- exp(lkel + etalkel)
    cl_milk <- exp(lcl_milk + etalcl_milk)
    cl_infant <- exp(lcl_infant + etalcl_infant) * (WT_INFANT / 5)^e_wt_infant_cl_infant

    # Supplementary Material S1 Eqs. 1-2 and 5-6: kmb = CLmb / Vm,
    # kbo = CLbo / Vb;
    # kmm is the TOTAL maternal elimination rate constant, so the lactation
    # flow kmb * Am is part of kmm * Am, not added to it.
    kmilkinf <- cl_milk / vc
    kel_infant <- cl_infant / vc_infant

    d/dt(central) <- -kel * central
    d/dt(infant_central) <- kmilkinf * central - kel_infant * infant_central

    # Amounts are g of D2O; x 1000 gives mg/L, the saliva enrichment in the
    # paper's units (upstream Stan 'amt*1000/Vm').
    Cc <- 1000 * central / vc
    Cinfant <- 1000 * infant_central / vc_infant

    Cc ~ add(addSd) + prop(propSd)
    Cinfant ~ add(addSd_Cinfant) + prop(propSd_Cinfant)
  })
}
