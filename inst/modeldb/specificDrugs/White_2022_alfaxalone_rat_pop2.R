White_2022_alfaxalone_rat_pop2 <- function() {
  description <- paste(
    "Preclinical (rat). One-compartment population PK model for",
    "intravenous alfaxalone (Alfaxan, 2-hydroxypropyl-beta-cyclodextrin",
    "formulation) given as a loading infusion followed by constant-rate",
    "infusions to adult Lewis and Sprague-Dawley rats of both sexes.",
    "'Population 2' (final) fit of White 2022 (52 rats: the 28 rats of",
    "population 1 plus 24 further Sprague-Dawley rats), fitted in Phoenix",
    "NLME on total (not per-kg) dose. Clearance carries a linear effect of",
    "log10-centred body weight (LCBW = log10(WT / 0.317 kg)) and an",
    "exponential sex effect; volume carries an exponential strain effect.",
    "Between-animal variability and residual error magnitudes are not",
    "reported and are held at zero, so simulations are typical-value. See",
    "White_2022_alfaxalone_rat_pop1 for the 28-rat 'population 1' fit used",
    "to design the female dosing regimen.",
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
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters clearance as LCBW = log10(WT / 0.317 kg), the paper's 'log",
        "of centralized body weight' (Methods 2.8.1, Eq 3). Neither the",
        "logarithm base nor the centring weight is printed. Both were",
        "recovered from the Figure S7 source data in the supplement, whose",
        "52 per-rat LCBW values all map to whole-gram weights only for",
        "log10 and a 317 g centre (weights 208-488 g). The same form",
        "reproduces the Table 4 male and female Lewis typical clearances",
        "(33.5 and 10.1 mL/min) at the Methods 2.2 group mean weights",
        "(308 and 222 g). Baseline, time-fixed per animal."
      ),
      source_name = "LCBW (log centralized body weight, in grams)"
    ),
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
    n_subjects = 52L,
    n_studies = 4L,
    age_range = "8-12 weeks (Lewis), 9-12 weeks (new Sprague-Dawley cohort)",
    weight_range = "208-488 g (Figure S7 source data)",
    sex_female_pct = 55.8,
    disease_state = paste(
      "Healthy, surgically instrumented (jugular venous and carotid arterial",
      "cannulae) under isoflurane, then maintained on alfaxalone anesthesia",
      "during nociceptive withdrawal reflex / diffuse noxious inhibitory",
      "control studies"
    ),
    dose_range = paste(
      "Loading infusion 1.67 mg/kg/min for 2.5 min, then constant-rate",
      "infusion: Lewis 0.75 mg/kg/min for 60 min then 0.52 mg/kg/min;",
      "new Sprague-Dawley males 0.75 mg/kg/min throughout; new",
      "Sprague-Dawley females 0.52 mg/kg/min for 60 min then",
      "0.42 mg/kg/min; White 2017 Sprague-Dawley regimen as described in",
      "that paper."
    ),
    regions = "United Kingdom (University of Nottingham)",
    notes = paste(
      "Population 2 = population 1 (16 Lewis: 9 male 308 +/- 49 g, 7",
      "female 222 +/- 9 g; 12 White 2017 Sprague-Dawley: 6 male, 6",
      "female) plus 24 Sprague-Dawley rats (8 male 422 +/- 41 g, 16 female",
      "304 +/- 15 g; Methods 2.2 and Figure 1). Female percentage =",
      "(7 + 6 + 16) / 52."
    )
  )

  ini({
    # Typical values are for a male Lewis rat at the 317 g centring weight
    # (LCBW = 0, SEXF = 0, STRAIN_SD = 0), per animal: the NLME fit used
    # total dose (Methods 2.8.1).
    lcl <- log(0.0352); label("Clearance, male Lewis rat at 317 g (L/min)") # Section 3.9 text below Eq 4: CL_TV = 35.2 mL/min
    lvc <- log(0.51); label("Volume of distribution, Lewis rat (L)") # Section 3.9 text below Eq 4: Vd_TV = 0.51 L

    e_wt_cl <- 3.64; label("Linear slope of CL on log10(WT / 0.317 kg) (unitless)") # Eq 3: (1 + LCBW * 3.64)
    e_sexf_cl <- -0.43; label("Exponential effect of female sex on CL (unitless)") # Eq 3: exp(-0.43 * sex covariate)
    e_strain_sd_vc <- -0.692; label("Exponential effect of Sprague-Dawley strain on Vd (unitless)") # Eq 4: exp(-0.692 * strain covariate); Table 4 SD Vd 0.26 L

    # Section 3.9: random effects on all parameters, diagonal omega. The
    # variances are not reported anywhere in the paper or its supplement.
    etalcl ~ fixed(0) # Eq 3 CL eta; variance not reported
    etalvc ~ fixed(0) # Eq 4 Vd eta; variance not reported

    # Methods 2.8.1: Phoenix 'mixed ratio' residual error,
    # C + CEps * (1 + C * CMixRatio): an additive SD plus a proportional SD
    # summed linearly (nlmixr2 combined1). Magnitudes are not reported.
    addSd <- fixed(0); label("Additive residual error (ug/mL)") # Methods 2.8.1 mixed ratio error; value not reported
    propSd <- fixed(0); label("Proportional residual error (fraction)") # Methods 2.8.1 mixed ratio error; value not reported
  })

  model({
    # LCBW: log10 of body weight centred at 317 g (0.317 kg); see the WT
    # covariate notes for how the base and centre were recovered.
    lcbw <- log10(WT / 0.317)
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF) * (1 + e_wt_cl * lcbw)
    vc <- exp(lvc + etalvc + e_strain_sd_vc * STRAIN_SD)
    kel <- cl / vc

    d / dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
