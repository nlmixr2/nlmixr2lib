Liu_2018_sorafenib <- function() {
  description <- "One-compartment population PK/PD model for oral sorafenib and its active N-oxide metabolite in adults with relapsed or refractory acute myeloid leukemia (Liu 2018). The parent structure is adopted from Jain 2011: a dosing compartment plus four catenary transit compartments sharing one rate constant ka, enterohepatic recirculation in which a fraction fENT of the central outflow CL * C enters a gallbladder that empties back into the last transit compartment at K64 only once the time after dose exceeds TENT (binary switch EHYN), and the remaining fraction 1 - fENT forming the N-oxide metabolite (one compartment, V7, CLm). Unbound sorafenib is the total concentration times an unbound fraction fu. FLT3-ITD and ERK activities (percent of each patient's baseline) follow inhibitory Emax models driven by a potency-adjusted concentration Cp + 14.59 * Cm * MW(sorafenib) / MW(N-oxide); Imax for FLT3 is fixed to 1. Body weight (reference 80 kg) scales both volumes linearly and both clearances with exponent 0.75."
  reference <- paste(
    "Liu T, Ivaturi V, Sabato P, Gobburu JVS, Greer JM, Wright JJ, Smith BD,",
    "Pratz KW, Rudek MA.",
    "Sorafenib Dose Recommendation in Acute Myeloid Leukemia Based on",
    "Exposure-FLT3 Relationship.",
    "Clin Transl Sci. 2018;11(4):435-443.",
    "doi:10.1111/cts.12555.",
    "Corrigendum: Clin Transl Sci. 2018;11(6):624. doi:10.1111/cts.12587",
    "(swaps the in-text references to Figures 2 and 3; no parameter or",
    "equation changes).",
    sep = " "
  )
  vignette <- "Liu_2018_sorafenib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Structural allometric covariate on all four disposition parameters, reference weight 80 kg: exponent 1 on V and V7, exponent 0.75 on CL and CLm (Methods, equations for CL, CLm, V, V7). The exponents were fixed by the authors, not estimated.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole (strong CYP3A4 inhibitor), 1 = yes",
      units = "(binary)",
      type = "binary",
      notes = "8 of 15 patients received voriconazole for fungal pneumonia. It was the only covariate with a clear effect (about 68.3% less N-oxide formed, Figure 3) but was not retained in the final model because of the small sample, unstable estimation and the known lack of effect of strong CYP3A4 inhibitors on parent exposure (Results, 'Population PK model')."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested, not retained (Results, 'Population PK model')."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested as 'gender', not retained (Results, 'Population PK model')."
    ),
    ECOG_GE1 = list(
      description = "Symptomatic baseline performance status (1 = symptomatic, 0 = asymptomatic)",
      units = "(binary)",
      type = "binary",
      notes = "The paper screened 'baseline disease status (asymptomatic and symptomatic)', 8:7 in Table 1; mapped here to the ECOG >= 1 indicator. Not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder = list(analyte = "sorafenib", units = "mg", specimen = "bile", verified = TRUE),
    central_noxsora = list(analyte = "sorafenib N-oxide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "37-74 years",
    age_median = "58 years",
    weight_range = "52.2-97.7 kg",
    weight_median = "74.6 kg",
    sex_female_pct = 46.7,
    disease_state = "Relapsed or refractory acute myeloid leukemia (FLT3-ITD activity measured ex vivo in TF-1 cells expressing FLT3-ITD).",
    dose_range = "Sorafenib 400 mg twice daily for 14 days (n = 3) or 21 days (n = 9), or 600 mg twice daily for 14 days (n = 3), in 28-day cycles.",
    regions = "USA (Johns Hopkins; NCT00131989).",
    co_medication = "Voriconazole in 8 of 15 patients (not retained as a covariate).",
    notes = "Table 1: BSA median 1.85 (1.56-2.26) m^2; BMI median 26.5 (18.5-36.3); asymptomatic:symptomatic 8:7; male:female 8:7. 141 total, 136 unbound sorafenib and 100 N-oxide concentrations from 15 patients; 35 FLT3 and 36 ERK activity measurements from 11 patients."
  )

  ini({
    # Table 2 FOCE-I estimates, typical values for an 80 kg patient.
    lka <- log(3.75); label("Absorption / transit rate constant Ka (1/h)") # Table 2 Ka = 3.75 1/h
    lcl <- log(14.1); label("Apparent total clearance of sorafenib from central CL (L/h)") # Table 2 CL = 14.10 L/h
    lvc <- log(171); label("Apparent volume of distribution of sorafenib V (L)") # Table 2 V = 171 L
    fent <- 0.502; label("Fraction of central outflow entering the gallbladder fENT (unitless)") # Table 2 fENT = 50.20%
    lkehc <- log(0.173); label("Gallbladder release rate constant back to the gut K64 (1/h)") # Table 2 K64 = 0.173 1/h
    ltprime <- log(5.48); label("Time after dose at which gallbladder release switches on TENT (h)") # Table 2 TENT = 5.48 h
    lfu <- log(0.0181); label("Unbound fraction of sorafenib in plasma fu (unitless)") # Table 2 fu = 1.81%
    lcl_noxsora <- log(116); label("Apparent clearance of sorafenib N-oxide CLm (L/h)") # Table 2 CLm = 116 L/h
    lvc_noxsora <- log(89.8); label("Apparent volume of distribution of sorafenib N-oxide V7 (L)") # Table 2 V7 = 89.80 L

    # Allometric exponents, fixed by the authors (Methods, CL / CLm / V / V7 equations).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL and CLm (unitless)") # Methods 'exponent of 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on V and V7 (unitless)") # Methods 'exponent of 1'

    # Potency of the N-oxide relative to sorafenib on FLT3-ITD, from in vitro data (Methods, Cadjusted equation).
    pot_noxsora <- fixed(14.59); label("FLT3-ITD potency of the N-oxide relative to sorafenib (unitless)") # Methods '14.59-fold higher potency'

    # Exposure-response (Table 2); no BSV could be estimated on the PD parameters.
    lic50_flt3 <- log(69.3); label("Cadjusted giving half-maximal FLT3 inhibition IC50FLT3 (ng/mL)") # Table 2 IC50FLT3 = 69.30 ng/mL
    imax_flt3 <- fixed(1); label("Maximal fractional inhibition of FLT3 activity (unitless)") # Results 'Imax of sorafenib on FLT3 was fixed to 100%'
    lic50_erk <- log(85.7); label("Cadjusted giving half-maximal ERK inhibition IC50ERK (ng/mL)") # Table 2 IC50ERK = 85.70 ng/mL
    imax_erk <- 0.84; label("Maximal fractional inhibition of ERK activity (unitless)") # Table 2 Imax,ERK = 0.84

    # BSV reported as CV%; omega^2 = log(1 + CV^2).
    etalka ~ 0.33787 # Table 2 Ka BSV 63.40% -> log(1 + 0.634^2)
    etalvc ~ 0.44053 # Table 2 V BSV 74.40% -> log(1 + 0.744^2)
    etalvc_noxsora ~ 0.91164 # Table 2 V7 BSV 122% -> log(1 + 1.22^2)
    etalcl ~ 0.26429 # Table 2 CL BSV 55.00% -> log(1 + 0.55^2)
    etalfu ~ 0.00015875 # Table 2 fu BSV 1.26% -> log(1 + 0.0126^2)
    etalcl_noxsora ~ 0.86243 # Table 2 CLm BSV 117% -> log(1 + 1.17^2)

    # Residual error (Table 2). The paper used one proportional error for total and unbound sorafenib (Methods).
    # nlmixr2 cannot share one error parameter between two endpoints, so the
    # unbound endpoint carries its own propSd_Cu at the same value.
    propSd <- 0.476; label("Proportional residual error, total sorafenib (fraction)") # Table 2 Proportional error sorafenib = 47.60%
    propSd_Cu <- 0.476; label("Proportional residual error, unbound sorafenib (fraction)") # Table 2 Proportional error sorafenib = 47.60% (shared with total sorafenib in the paper)
    propSd_noxsora <- 0.424; label("Proportional residual error, sorafenib N-oxide (fraction)") # Table 2 Proportional error N-oxide = 42.40%
    addSd_noxsora <- 1.06; label("Additive residual error, sorafenib N-oxide (ng/mL)") # Table 2 Additive error N-oxide = 1.06 ng/mL
    propSd_flt3_activity <- 0.879; label("Proportional residual error, FLT3 activity (fraction)") # Table 2 Proportional error FLT3 = 87.90%
    propSd_erk_activity <- 0.394; label("Proportional residual error, ERK activity (fraction)") # Table 2 Proportional error ERK = 39.40%
  })
  model({
    # Molecular weights used to convert the N-oxide mass concentration to
    # sorafenib-equivalent molar units (Methods, 'Assessment of
    # exposure-response relationship').
    mw_sorafenib <- 464.825
    mw_noxsora <- 480.825

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc
    cl_noxsora <- exp(lcl_noxsora + etalcl_noxsora) * (WT / 80)^e_wt_cl
    vc_noxsora <- exp(lvc_noxsora + etalvc_noxsora) * (WT / 80)^e_wt_vc
    kehc <- exp(lkehc)
    tprime <- exp(ltprime)
    fu <- exp(lfu + etalfu)

    kel <- cl / vc
    kel_noxsora <- cl_noxsora / vc_noxsora

    # EHYN: binary switch that opens gallbladder release once the time
    # after the most recent dose reaches TENT (Methods, Figure 1).
    ehyn <- 0
    if (tad() >= tprime) ehyn <- 1

    # Figure 1: Aa -> A1 -> A2 -> A3 -> A4 -> A5 all at Ka; A5 -> A6 at
    # CL * fENT; A6 -> A4 at K64 * EHYN; A5 -> A7 at CL * (1 - fENT); A7
    # eliminated at CLm. Metabolite formation is carried in parent mass
    # units, so V7 and CLm are apparent values (Methods, AUCm / AUC =
    # CL * (1 - fENT) / CLm).
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(transit4) <- ka * transit3 - ka * transit4 + ehyn * kehc * gallbladder
    d/dt(central) <- ka * transit4 - kel * central
    d/dt(gallbladder) <- fent * kel * central - ehyn * kehc * gallbladder
    d/dt(central_noxsora) <- (1 - fent) * kel * central - kel_noxsora * central_noxsora

    # Dose in mg, volumes in L: mg/L * 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cu <- Cc * fu
    Cc_noxsora <- 1000 * central_noxsora / vc_noxsora

    # Potency-adjusted sorafenib concentration (Methods, Cadjusted equation).
    cadjusted <- Cc + pot_noxsora * Cc_noxsora * mw_sorafenib / mw_noxsora
    ic50_flt3 <- exp(lic50_flt3)
    ic50_erk <- exp(lic50_erk)
    flt3_activity <- 100 * (1 - imax_flt3 * cadjusted / (ic50_flt3 + cadjusted))
    erk_activity <- 100 * (1 - imax_erk * cadjusted / (ic50_erk + cadjusted))

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
    Cc_noxsora ~ add(addSd_noxsora) + prop(propSd_noxsora)
    flt3_activity ~ prop(propSd_flt3_activity)
    erk_activity ~ prop(propSd_erk_activity)
  })
}
