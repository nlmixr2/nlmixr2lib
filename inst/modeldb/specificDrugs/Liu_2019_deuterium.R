Liu_2019_deuterium <- function() {
  description <- "Mother-infant dyad model for deuterium oxide (D2O) given orally to the lactating mother (the dose-to-the-mother technique for assessing exclusive breastfeeding). Two linked one-compartment models: the mother's body water receives the oral dose and loses D2O by a first-order total elimination rate constant, part of which is the lactation clearance into the breastfed infant; the infant's body water receives that flow and loses D2O by its own water clearance. The infant volume of distribution is not estimated but computed from infant body weight (total body water x 1.041) and grows with the infant. Covariates: maternal weight on the maternal volume, infant weight on the infant clearance. Fitted in a fully Bayesian framework (Stan) to maternal and infant saliva enrichment."
  reference <- paste(
    "Liu Z, Diana A, Slater C, Preston T, Gibson RS, Houghton L, Duffull SB.",
    "Development of a nonlinear hierarchical model to describe the disposition",
    "of deuterium in mother-infant pairs to assess exclusive breastfeeding practice.",
    "J Pharmacokinet Pharmacodyn. 2019;46(1):1. doi:10.1007/s10928-018-9613-x.",
    "Structure from Eqs. 3-4 and the Stan code in Supplementary material 4;",
    "parameter values from Eqs. 16-17 and Supplementary materials 6 and 7.",
    sep = " "
  )
  vignette <- "Liu_2019_deuterium"
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
      notes = "Baseline. Power effect on the maternal D2O volume of distribution, (WT / 70)^0.62 (Eq. 16; Stan code 'pow(MWT[i]/70, theta[1])'). The paper calls 70 kg the median maternal weight, but Table 1 gives the cohort median as 53.1 kg (34.5-93.1); 70 kg is the value the fitted code divides by. Source column MWT.",
      source_name = "MWT"
    ),
    WT_INFANT = list(
      description = "Breastfed infant body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Enters twice. (1) The infant D2O volume of distribution Vb = 1.041 * exp(-0.427 + 0.963 * log(WT_INFANT)) (Eqs. A1.9-A1.10), where the Stan code interpolates infant weight linearly between the start-of-study and end-of-study weighings ('Wt[j] = Wt_s[j] + (Wt_e[j]-Wt_s[j])*time[j]/time[Tstop[i]]'), so supply WT_INFANT at every record (or use linear covariate interpolation) to reproduce the growth of Vb. (2) A power effect on the infant water clearance, (WT_INFANT / 5)^0.17 (Eq. 17); the paper evaluates this term at the infant's AVERAGE weight over the study (Stan 'BWT', 'infant average weight during the study'), whereas here it follows the current weight. Over the study's 14 days that difference changes the infant clearance by well under 1%. Table 1: start weight 5.9 kg (3.9-8.4), end weight 6.2 kg (4.3-8.5). Source columns Wt_s / Wt_e (volume) and BWT (clearance).",
      source_name = "BWT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 113,
    n_studies = 1,
    age_range = "Mothers 16-42 years; infants 2.0-5.4 months",
    age_median = "Mothers 25 years; infants 3.3 months",
    weight_range = "Mothers 34.5-93.1 kg; infants 3.9-8.4 kg at study start",
    weight_median = "Mothers 53.1 kg; infants 5.9 kg at study start and 6.2 kg at study end",
    sex_female_pct = NA_real_,
    disease_state = "Healthy lactating mothers and their healthy, singleton, full-term (> 37 weeks gestation, birth weight > 2500 g) infants identified as exclusively breastfed and verified as such by home observation over six non-consecutive days during the 14-day protocol.",
    dose_range = "Single oral dose of about 30 g deuterium oxide to the mother (median 30.0 g, range 30.0-30.2), diluted in about 50 g drinking water.",
    regions = "Sumedang (Tanjunsari, Sukasari and Pamulihan subdistricts), West Java, Indonesia.",
    infant_partner = "113 infants, 54 male / 59 female (Table 1). Saliva sampled from mother and infant pre-dose and on days 1, 2, 5, 6, 13 and 14; duplicate samples within 1 h were averaged.",
    notes = "121 mother-infant pairs recruited; 7 removed as non-exclusively breastfed on field observation and 1 as biologically implausible (infant water clearance > 40% of body weight), leaving 113 pairs and 1500 paired observations (Results, Data; Table 1). All mothers are female; sex_female_pct is left NA because the modelled dyad includes the infant, of whom 52% were female."
  )

  ini({
    # Typical values are the population means of the log-scale individual
    # parameters (Stan 'log_theta_1'), at maternal weight 70 kg and infant
    # weight 5 kg. Supplementary Table S6.1 prints the arithmetic population
    # means exp(mu + omega^2 / 2); Eqs. 16-17 and Supplement 7 print mu.
    lvc <- 3.49
    label("Maternal D2O volume of distribution Vm at maternal WT = 70 kg (L)") # Eq. 16 'ln Vm,i = N(3.49, 0.07) + N(0.62, 0.04) ln(MWTi/70kg)'; exp(3.49 + 0.07^2/2) = 32.87 vs Table S6.1 'Vm (L) 32.86'
    lkel <- -2.2277
    label("Maternal total D2O elimination rate constant kmm, including the lactation route (1/day)") # Table S6.1 'Kmm (1/day) 0.109, BSV 15.0%': log(0.109) - 0.15^2/2 = -2.2277 (mu is not printed for kmm)
    lcl_milk <- -0.33
    label("Lactation-route water clearance from mother to infant CLmb (L/day)") # Supplement 7 'log(CLmb) ~ N(-0.33, 0.20)'; exp(-0.33 + 0.201^2/2) = 0.733 = Table S6.1 'CLmb (L/day) 0.733'
    lcl_infant <- -0.16
    label("Infant water clearance to the outside CLbo at infant WT = 5 kg (L/day)") # Eq. 17 'ln CLbo,i = N(-0.16, 0.17) + N(0.17, 0.03) ln(BWTi/5kg)'; exp(-0.16 + 0.171^2/2) = 0.864 = Table S6.1 'CLbo (L/day) 0.864'

    e_wt_vc <- 0.62
    label("Power exponent of maternal weight on Vm (unitless)") # Eq. 16 'N(0.62, 0.04)'
    e_wt_infant_cl_infant <- 0.17
    label("Power exponent of infant weight on CLbo (unitless)") # Eq. 17 'N(0.17, 0.03)'

    # Between-subject variability: log-normal (Eq. 8). Variances are the
    # squares of the log-scale SDs; Supplement 7 states that the SD in
    # 'N(mean, SD)' is the between-subject variability, and Table S6.1 prints
    # it as 'BSV %' = 100 x SD. Only the CLmb-CLbo correlation is reported
    # (Supplement 7, 0.99); every other correlation of the fitted LKJ matrix
    # is unreported and set to zero.
    etalvc ~ 0.0049 # Table S6.1 'Vm BSV 7.0%' and Eq. 16 SD 0.07: 0.070^2
    etalkel ~ 0.0225 # Table S6.1 'Kmm BSV 15.0%': 0.150^2
    etalcl_milk + etalcl_infant ~ c(0.040401, 0.034028, 0.029241) # Table S6.1 'CLmb BSV 20.1%' and 'CLbo BSV 17.1%': 0.201^2, 0.99 x 0.201 x 0.171, 0.171^2; correlation 0.99 from Supplement 7 note 2

    # Combined residual error. The paper uses ONE (sigma1, sigma2) pair for
    # both mother and infant (Stan 'sqrt(ipred^2 * sigma[1]^2 + sigma[2]^2)');
    # rxode2 needs a separate pair per endpoint, so the pair is written twice.
    # The posterior values of sigma are not reported anywhere in the paper or
    # supplements, so both pairs are held at zero.
    propSd <- fixed(0)
    label("Proportional residual SD, maternal saliva (fraction; not reported, held at zero)") # Stan code 'sigma[1]'; value not reported
    addSd <- fixed(0)
    label("Additive residual SD, maternal saliva (mg/L; not reported, held at zero)") # Stan code 'sigma[2]'; value not reported
    propSd_Cinfant <- fixed(0)
    label("Proportional residual SD, infant saliva (fraction; same sigma as the mother in the source; not reported, held at zero)") # Stan code 'sigma[1]'; value not reported
    addSd_Cinfant <- fixed(0)
    label("Additive residual SD, infant saliva (mg/L; same sigma as the mother in the source; not reported, held at zero)") # Stan code 'sigma[2]'; value not reported
  })

  model({
    # Infant D2O volume of distribution from infant weight: total body water
    # ln(TBW) = a + b ln(WT) (Eq. A1.9, a ~ N(-0.427, 0.012), b ~ N(0.963,
    # 0.005), held at their means), enlarged by 4.1% for non-aqueous isotope
    # exchange (Eq. A1.10). Stan: 'Vb_b[j] = 1.041*exp(a + b*log(Wt[j]))'.
    tbw_infant <- exp(-0.427 + 0.963 * log(WT_INFANT))
    vc_infant <- 1.041 * tbw_infant

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- exp(lkel + etalkel)
    cl_milk <- exp(lcl_milk + etalcl_milk)
    cl_infant <- exp(lcl_infant + etalcl_infant) * (WT_INFANT / 5)^e_wt_infant_cl_infant

    # Eqs. 1-4: kmb = CLmb / Vm, kbo = CLbo / Vb; kmm is the TOTAL maternal
    # elimination rate constant, so the lactation flow kmb * Am is part of
    # kmm * Am, not added to it.
    kmilkinf <- cl_milk / vc
    kel_infant <- cl_infant / vc_infant

    d/dt(central) <- -kel * central
    d/dt(infant_central) <- kmilkinf * central - kel_infant * infant_central

    # Amounts are g of D2O; x 1000 gives mg/L (Stan 'amt*1000/Vm'), the
    # saliva enrichment in the paper's units (density of dilute heavy water
    # taken as 1 kg/L).
    Cc <- 1000 * central / vc
    Cinfant <- 1000 * infant_central / vc_infant

    Cc ~ add(addSd) + prop(propSd)
    Cinfant ~ add(addSd_Cinfant) + prop(propSd_Cinfant)
  })
}
