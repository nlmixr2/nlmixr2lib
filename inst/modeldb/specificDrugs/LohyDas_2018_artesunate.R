LohyDas_2018_artesunate <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model of oral artesunate (ARS) and",
    "its active metabolite dihydroartemisinin (DHA) in 50 adult patients with",
    "uncomplicated, artemisinin sensitive or resistant Plasmodium falciparum",
    "malaria in southern Myanmar (Lohy Das 2018, Malaria Journal). ARS",
    "absorption is described by a 3-transit-compartment chain (n = 3 fixed)",
    "followed by a one-compartment ARS disposition; complete in-vivo conversion",
    "of ARS to DHA is assumed (all ARS clearance is metabolic conversion).",
    "DHA disposition is one-compartment. Allometric body-weight scaling is",
    "applied to all CL (exponent 0.75) and V (exponent 1.0) parameters,",
    "centered on the population-median 50 kg. F is fixed at 1. The packaged",
    "model file omits the published time-varying parasite-density covariates",
    "on MTT and on F (Eqs. 3 and 4) and the entire PD layer (mixture-Emax",
    "parasite-killing model with effect compartment); both depend on the",
    "upstream Lohy Das 2017 AAPS J paper (ref [36]) which is not on disk.",
    "See the vignette's Assumptions and deviations section for the rationale."
  )
  reference <- paste(
    "Lohy Das JP, Kyaw MP, Nyunt MH, Chit K, Aye KH, Aye MM, Karlsson MO,",
    "Bergstrand M, Tarning J (2018).",
    "Population pharmacokinetic and pharmacodynamic properties of artesunate",
    "in patients with artemisinin sensitive and resistant infections in",
    "Southern Myanmar.",
    "Malaria Journal 17:126.",
    "doi:10.1186/s12936-018-2278-5.",
    sep = " "
  )
  vignette <- "LohyDas_2018_artesunate"
  units    <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at enrollment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at enrollment. Lohy Das 2018 Results p.5: 'Allometric",
        "scaling of all disposition parameters, centered by the median",
        "weight of 50 kg improved the model fit.' Methods p.4: 'Scaled body",
        "weight was raised to the power of 0.75 and 1 for clearance and",
        "volume parameters, respectively, and centered on the median weight",
        "of the population.' Applied to both ARS and DHA disposition",
        "parameters. Cohort weight range 40-60 kg (Table 1 IQR 46.0-53.5)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "18-55 years (Table 1 median 25.5, IQR 21.5-39.5)",
    weight_range   = "median 50.0 kg (Table 1 IQR 46.0-53.5)",
    weight_typical = "50 kg (population median; allometric reference)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Southeast Asian (Burmese)",
    disease_state  = paste(
      "Uncomplicated Plasmodium falciparum malaria; mono-infection with",
      "asexual parasite density 10,000-100,000/uL at enrollment. Inclusion",
      "criteria: 18-55 years, P. falciparum monoinfection, fever in last",
      "24 h, ability to tolerate oral artesunate, written informed consent.",
      "Exclusion criteria: severe malaria, severe malnutrition, pregnancy,",
      "lactation, mixed malaria infection, infection other than malaria,",
      "history of chronic medical illness, splenectomy, hypersensitivity",
      "to artesunate, prior anti-malarial drug use within 48 h."
    ),
    dose_range     = paste(
      "Oral artesunate monotherapy 4 mg/kg/day once daily for 7 days,",
      "administered with 8 oz of milk (directly observed); study drug",
      "Guilin Pharmaceutical Co. Ltd., lot AS091001. Plasma concentration",
      "measurements were collected after the FIRST dose only (pre-dose and",
      "at 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 3, 4, 6, and 8 h); the popPK fit",
      "thus describes single-dose disposition. A 50 kg subject receives",
      "200 mg artesunate per dose, equivalent to about 520,260 nmol using",
      "the ARS molar mass 384.42 g/mol."
    ),
    sampling       = paste(
      "Frequent first-dose sampling at 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 3,",
      "4, 6, and 8 h post-dose. ARS and DHA plasma concentrations were",
      "measured by LC-MS/MS (Lourens / Lindegardh assay); LLOQ 1.2 ng/mL",
      "(ARS) and 2.0 ng/mL (DHA), i.e. 3.12 nM (ARS) and 7.02 nM (DHA) in",
      "molar units. M3 method was used for samples below LLOQ."
    ),
    regions        = paste(
      "Myanmar (Kawthaung, southern Myanmar, Palm Tree plantation site",
      "hospital). Australian New Zealand Clinical Trials Registry",
      "ACTRN12610000896077; conducted in 2011."
    ),
    notes          = paste(
      "Of 53 patients recruited, 1 was excluded for not meeting",
      "inclusion/exclusion criteria and 2 were excluded from the PK",
      "analysis for missing covariates and PK/PD data, leaving 50 patients",
      "in the final PK analysis (Results p.5). 56.1% of patients were",
      "estimated to have artemisinin-resistant infections (paper's mixture",
      "model on Emax, Table 2 P_MIX,resistant = 56.1%); the PD mixture",
      "layer is not encoded in this PK-only extraction."
    )
  )

  ini({
    # Structural PK parameters (Lohy Das 2018 Table 2 PK section).
    # Concentrations were modelled on a natural-log molar scale (nmol/L)
    # with the ARS dose pre-converted to nmols (Methods p.4: 'Observed
    # concentrations of ARS and DHA (molar units) were transformed into
    # their natural logarithms and modelled simultaneously'). Typical
    # values correspond to apparent clearance / volume at the cohort
    # median WT = 50 kg (Results p.5).

    lfdepot <- fixed(log(1))
    label("Reference relative oral bioavailability of artesunate, F (unitless); fixed at 1 (the source paper fixes F at 100%, Table 2)")  # Lohy Das 2018 Table 2: F (%) = 100 fix

    lmtt    <- log(1.34)
    label("Mean transit time of the 3-compartment transit-absorption chain, MTT (h)")  # Lohy Das 2018 Table 2: MTT = 1.34 h (%RSE 18.8, 95% CI 1.04-1.96)

    lcl     <- log(1750)
    label("Apparent artesunate elimination clearance, CL_ARS/F at WT = 50 kg (L/h); equivalent to the ARS-to-DHA conversion clearance under the source paper's assumption of complete in-vivo conversion")  # Lohy Das 2018 Table 2: CL_ARS/F = 1750 L/h (%RSE 8.55, 95% CI 1570-2090)

    lvc     <- log(1300)
    label("Apparent artesunate central volume of distribution, V_ARS/F at WT = 50 kg (L)")  # Lohy Das 2018 Table 2: V_ARS/F = 1300 L (%RSE 12.6, 95% CI 1110-1660)

    lcl_dihydroart <- log(76.7)
    label("Apparent dihydroartemisinin elimination clearance, CL_DHA/F at WT = 50 kg (L/h)")  # Lohy Das 2018 Table 2: CL_DHA/F = 76.7 L/h (%RSE 6.99, 95% CI 69.9-87.8)

    lvc_dihydroart <- log(102)
    label("Apparent dihydroartemisinin central volume of distribution, V_DHA/F at WT = 50 kg (L)")  # Lohy Das 2018 Table 2: V_DHA/F = 102 L (%RSE 8.95, 95% CI 89.5-119.0)

    # Allometric exponents. Lohy Das 2018 Methods p.4: 'Clearance and
    # volume of distribution of both parent and metabolite were scaled
    # allometrically using body weight. Scaled body weight was raised to
    # the power of 0.75 and 1 for clearance and volume parameters,
    # respectively.' The paper does not explicitly state 'FIXED' but the
    # allometric exponents are conventional and treated here as fixed
    # structural choices following nlmixr2lib's encoding for the sibling
    # Birgersson_2019_artesunate and Hendriksen_2013_artesunate models.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on clearance (ARS and DHA), unitless")  # Lohy Das 2018 Methods p.4: allometric power 0.75 on clearance
    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on central volume (ARS and DHA), unitless")  # Lohy Das 2018 Methods p.4: allometric power 1 on volume

    # Inter-individual variability (log-normal). Lohy Das 2018 Table 2
    # reports BSV as %CV with the footnote 'Coefficient of variation
    # (%CV) of between subject variability (BSV) was calculated as
    # 100 * (variance-1)^(1/2)' (read as the standard log-normal
    # transformation %CV = 100 * sqrt(exp(variance) - 1), i.e.
    # variance = log(1 + (CV/100)^2) on the eta scale). The model file
    # uses the variance values derived from the published %CVs.
    etalfdepot ~ 0.0929   # Lohy Das 2018 Table 2: BSV F = 31.2% (%RSE 29.4); variance = log(1 + 0.312^2) = 0.0929
    etalmtt    ~ 0.5468   # Lohy Das 2018 Table 2: BSV MTT = 85.3% (%RSE 24.9); variance = log(1 + 0.853^2) = 0.5468
    etalcl     ~ 0.0693   # Lohy Das 2018 Table 2: BSV CL_ARS = 26.8% (%RSE 44.3); variance = log(1 + 0.268^2) = 0.0693
    etalvc     ~ 0.4434   # Lohy Das 2018 Table 2: BSV V_ARS = 74.7% (%RSE 27.3); variance = log(1 + 0.747^2) = 0.4434
    etalcl_dihydroart ~ 0.0444   # Lohy Das 2018 Table 2: BSV CL_DHA = 21.3% (%RSE 30.3); variance = log(1 + 0.213^2) = 0.0444
    etalvc_dihydroart ~ 0.0953   # Lohy Das 2018 Table 2: BSV V_DHA = 31.6% (%RSE 40.5); variance = log(1 + 0.316^2) = 0.0953

    # Residual error. Lohy Das 2018 Methods p.4: 'unexplained residual
    # variability (RUV) was estimated by separate additive error models
    # for log-transformed ARS and DHA concentrations (i.e. equal to
    # exponential error models on an arithmetic scale).' By the standing
    # nlmixr2lib convention (see references/parameter-names.md and the
    # sibling Tan_2009_artesunate, Birgersson_2019_artesunate models),
    # NONMEM additive-on-log-scale residual maps to nlmixr2 proportional
    # residual in linear space, with propSd = SD on the log scale ~= CV
    # in linear space to first order. Table 2 reports RUV as a CV%; the
    # value is used here directly as propSd on the decimal scale.
    propSd     <- 0.732
    label("Proportional residual SD for artesunate plasma concentration (CV in linear space, equivalent to SD on log scale)")  # Lohy Das 2018 Table 2: RUV ARS = 73.2% (%RSE 3.95)
    propSd_dihydroart <- 0.585
    label("Proportional residual SD for dihydroartemisinin plasma concentration (CV in linear space, equivalent to SD on log scale)")  # Lohy Das 2018 Table 2: RUV DHA = 58.5% (%RSE 3.34)
  })

  model({
    # Individual PK parameters. Allometric body-weight scaling on CL with
    # exponent 0.75 and on V with exponent 1, both centered on the
    # cohort median WT = 50 kg (Lohy Das 2018 Methods p.4, Results p.5).
    cl     <- exp(lcl     + etalcl)     * (WT / 50)^e_wt_cl
    vc     <- exp(lvc     + etalvc)     * (WT / 50)^e_wt_vc
    cl_dihydroart <- exp(lcl_dihydroart + etalcl_dihydroart) * (WT / 50)^e_wt_cl
    vc_dihydroart <- exp(lvc_dihydroart + etalvc_dihydroart) * (WT / 50)^e_wt_vc
    mtt    <- exp(lmtt    + etalmtt)

    # Transit-absorption chain rate (Savic 2007 parameterisation, also
    # used by the sibling Birgersson_2019_artesunate model). With
    # n_transit = 3 fixed (Lohy Das 2018 Results p.5: 'The absorption
    # was described by a transit compartment (n = 3) model'), the chain
    # has n_transit + 1 = 4 first-order steps and ktr = 4 / MTT.
    nn  <- 3
    ktr <- (nn + 1) / mtt

    # Conversion / elimination micro-constants. Under the source paper's
    # assumption of complete in-vivo conversion of ARS to DHA (Methods
    # p.4: 'Complete metabolic in vivo conversion of ARS into DHA was
    # assumed throughout modelling'), all ARS clearance is metabolic
    # conversion, so the rate from central -> central_dihydroart is cl/vc. DHA
    # is then eliminated linearly at cl_dihydroart/vc_dihydroart. Doses and
    # concentrations are tracked on a molar basis (nmol / nmol/L), so no
    # molecular-weight conversion is needed at the parent-to-metabolite
    # transfer step (the paper modelled on the same molar basis).
    kel_ars <- cl     / vc
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    # ODE system: oral depot -> 3 transit compartments -> ARS central ->
    # DHA central -> elimination. The same first-order rate ktr governs
    # every absorption step in the transit chain. Lohy Das 2018 Fig. 1
    # shows the same compartment layout (depot -> transit1 -> transit2 ->
    # transit3 -> ARS_central -> DHA_central, with the PD parasite
    # compartment and effect-compartment links omitted here per the
    # vignette Assumptions and deviations section).
    d/dt(depot)        <- -ktr * depot
    d/dt(transit1)     <-  ktr * depot    - ktr * transit1
    d/dt(transit2)     <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)     <-  ktr * transit2 - ktr * transit3
    d/dt(central)      <-  ktr * transit3 - kel_ars * central
    d/dt(central_dihydroart)  <-  kel_ars * central - kel_dihydroart * central_dihydroart

    # Bioavailability on the depot compartment. lfdepot is fixed at
    # log(1) (F = 100% fix in Table 2) so exp(lfdepot) = 1 at the
    # population mean; etalfdepot carries the 31.2% IIV on F.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Plasma concentrations in nmol/L (dose in nmol, V in L).
    Cc     <- central     / vc
    Cc_dihydroart <- central_dihydroart / vc_dihydroart

    # Residual error: NONMEM additive-on-log-scale maps to nlmixr2
    # proportional in linear space (see ini() comments).
    Cc     ~ prop(propSd)
    Cc_dihydroart ~ prop(propSd_dihydroart)
  })
}
