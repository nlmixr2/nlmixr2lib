White_2022_alfaxalone_rat_pop1 <- function() {
  description <- paste(
    "Preclinical (rat). One-compartment population PK model for",
    "intravenous alfaxalone (Alfaxan, 2-hydroxypropyl-beta-cyclodextrin",
    "formulation) given as a loading infusion followed by constant-rate",
    "infusions to adult Lewis and Sprague-Dawley rats of both sexes.",
    "'Population 1' fit of White 2022 (28 rats: 16 Lewis plus 12",
    "Sprague-Dawley rats from White 2017), fitted in Phoenix NLME on total",
    "(not per-kg) dose. Clearance carries exponential sex and strain",
    "effects and volume an exponential strain effect; body weight is not a",
    "covariate. This is the model the authors used to design the adjusted",
    "female Sprague-Dawley infusion regimen (Figure 6). Between-animal",
    "variability and residual error magnitudes are not reported and are",
    "held at zero, so simulations are typical-value. See",
    "White_2022_alfaxalone_rat_pop2 for the 52-rat 'population 2' refit.",
    sep = " "
  )
  reference <- paste(
    "White K, Aldurdunji M, Harris J, Ortori C, Paine S.",
    "Alfaxalone population pharmacokinetics in the rat: Model application",
    "for pharmacokinetic and pharmacodynamic design in inbred and outbred",
    "strains and sexes.",
    "Pharmacol Res Perspect. 2022;10(6):e01031.",
    "doi:10.1002/prp2.1031.",
    sep = " "
  )
  vignette <- "White_2022_alfaxalone_rat"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Methods 2.8.1: 'Categorical covariates were implemented for sex",
        "(male = 0, female = 1)'. Same orientation as the canonical, so no",
        "value transformation."
      ),
      source_name = "sex covariate (gender)"
    ),
    STRAIN_SD = list(
      description = "Rat strain indicator, 1 = Sprague-Dawley, 0 = Lewis",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Methods 2.8.1: 'strain (Lewis = 0, SD = 1)'. The reference strain",
        "in this model is the inbred Lewis rat (not the Wistar-Kyoto",
        "lineage of the founding STRAIN_SD example)."
      ),
      source_name = "strain covariate (Strain)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "alfaxalone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Lewis and Sprague-Dawley)",
    n_subjects = 28L,
    n_studies = 2L,
    age_range = "8-12 weeks (Lewis); Sprague-Dawley adults from White 2017",
    weight_range = paste(
      "Lewis males 308 +/- 49 g, Lewis females 222 +/- 9 g (mean +/- SD);",
      "the White 2017 Sprague-Dawley weights are not reported in White 2022"
    ),
    sex_female_pct = 46.4,
    disease_state = paste(
      "Healthy, surgically instrumented (jugular venous and carotid arterial",
      "cannulae) under isoflurane, then maintained on alfaxalone anesthesia",
      "during a nociceptive withdrawal reflex study"
    ),
    dose_range = paste(
      "Loading infusion 1.67 mg/kg/min for 2.5 min, then constant-rate",
      "infusion: Lewis 0.75 mg/kg/min for 60 min then 0.52 mg/kg/min to",
      "the end of the experiment; White 2017 Sprague-Dawley regimen as",
      "described in that paper."
    ),
    regions = "United Kingdom (University of Nottingham)",
    notes = paste(
      "Population 1 = 16 Lewis rats (9 male, 7 female; Methods 2.2) plus 6",
      "male and 6 female Sprague-Dawley rats from White 2017 (Vet Anaesth",
      "Analg 44:865-875; Figure 1). Plasma alfaxalone by LC-MS/MS, LLOQ",
      "200 ng/mL (Methods 2.7). Female percentage = (7 + 6) / 28."
    )
  )

  ini({
    # Typical values are for a male Lewis rat (SEXF = 0, STRAIN_SD = 0) and
    # are per animal, not per kg: the NLME fit used total dose (Methods 2.8.1).
    lcl <- log(0.0252); label("Clearance, male Lewis rat (L/min)") # Section 3.4 text below Eq 2: CL_TV = 25.2 mL/min
    lvc <- log(0.57); label("Volume of distribution, Lewis rat (L)") # Section 3.4 text below Eq 2: Vd_TV = 0.57 L

    e_sexf_cl <- -0.841; label("Exponential effect of female sex on CL (unitless)") # Eq 1: exp(-0.841 * sex covariate)
    e_strain_sd_cl <- 0.478; label("Exponential effect of Sprague-Dawley strain on CL (unitless)") # Eq 1: exp(0.478 * strain covariate)
    e_strain_sd_vc <- -0.0237; label("Exponential effect of Sprague-Dawley strain on Vd (unitless)") # Eq 2: exp(-0.0237 * strain covariate)

    # Section 3.4: random effects on all parameters, diagonal omega. The
    # variances are not reported anywhere in the paper or its supplement.
    etalcl ~ fixed(0) # Eq 1 CL eta; variance not reported
    etalvc ~ fixed(0) # Eq 2 Vd eta; variance not reported

    # Methods 2.8.1: Phoenix 'mixed ratio' residual error,
    # C + CEps * (1 + C * CMixRatio): an additive SD plus a proportional SD
    # summed linearly (nlmixr2 combined1). Magnitudes are not reported.
    addSd <- fixed(0); label("Additive residual error (ug/mL)") # Methods 2.8.1 mixed ratio error; value not reported
    propSd <- fixed(0); label("Proportional residual error (fraction)") # Methods 2.8.1 mixed ratio error; value not reported
  })

  model({
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_strain_sd_cl * STRAIN_SD)
    vc <- exp(lvc + etalvc + e_strain_sd_vc * STRAIN_SD)
    kel <- cl / vc

    d / dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
