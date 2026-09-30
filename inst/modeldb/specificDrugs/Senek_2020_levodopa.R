Senek_2020_levodopa <- function() {
  description <- "One-compartment population PK model for levodopa given as an intrajejunal gel infusion (levodopa-carbidopa intestinal gel, LCIG, or levodopa-entacapone-carbidopa intestinal gel, LECIG) in advanced Parkinson's disease, with fixed fast first-order absorption from the jejunum, a one-transit-compartment absorption branch for night-time oral levodopa-carbidopa tablets, allometric body-weight scaling of CL/F and V/F, and a fractional shift in CL/F (with its own IIV) during simultaneous entacapone infusion (Senek 2020)"
  reference <- paste(
    "Senek M, Nyholm D, Nielsen EI.",
    "Population pharmacokinetics of levodopa gel infusion in Parkinson's disease:",
    "effects of entacapone infusion and genetic polymorphism.",
    "Sci Rep. 2020;10:18057.",
    "doi:10.1038/s41598-020-75052-2"
  )
  vignette <- "Senek_2020_levodopa"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    depot_oral = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levodopa", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling with reference weight 70 kg and fixed exponents 0.75 on CL/F and 1 on V/F (Senek 2020 Methods, 'Model development'; Table 2 units L/h/70 kg and L/70 kg).",
      source_name = "Weight"
    ),
    CONMED_ENTACAPONE = list(
      description = "Simultaneous entacapone infusion indicator: 1 = the levodopa gel is the levodopa-entacapone-carbidopa intestinal gel (LECIG), 0 = levodopa-carbidopa intestinal gel (LCIG) without entacapone",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (LCIG, no entacapone)",
      notes = "Time-varying in the source crossover (one LCIG day and one LECIG day per patient). Enters as a fractional shift CL/F_i = TVCL/F_LCIG * exp(etalcl) * (WT/70)^0.75 * (1 + e_conmed_entacapone_cl * exp(etae_conmed_entacapone_cl) * CONMED_ENTACAPONE) (Senek 2020 Table 2 footnote a). Set it to the treatment of the most recent gel dose; the paper does not separate the entacapone effect from the gel formulation.",
      source_name = "(treatment arm: LECIG vs LCIG)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    age_range = "63-76 years",
    age_median = "70 years",
    weight_range = "51-99 kg",
    weight_median = "73 kg",
    sex_female_pct = 36.4,
    disease_state = "Advanced Parkinson's disease on established levodopa-carbidopa intestinal gel (LCIG) treatment (duration of PD 8-23 years, duration of LCIG 0.2-7.6 years)",
    dose_range = "LCIG morning bolus 41-217 mg and continuous maintenance 363-1367 mg levodopa over a 14 h infusion day; LECIG at 80% or 90% of the LCIG morning dose and 80% of the LCIG maintenance and extra-bolus doses; a ~3 mL (60 mg levodopa) end-of-day tube flush; night-time oral levodopa-carbidopa tablets allowed until 3 h before infusion start",
    regions = "Sweden (Uppsala)",
    notes = "Randomised, open-label, two-day LCIG/LECIG crossover (Senek 2017 Mov Disord trial). Demographics from Senek 2020 Table 1 (n = 11, 7 male / 4 female). Genotypes (COMT rs4680, DDC rs921451, DDC rs3837091) were explored graphically on empirical Bayes CL/F estimates only and are not model covariates."
  )

  ini({
    lka <- fixed(log(50)); label("Absorption rate constant from the jejunal gel depot (1/h)") # Table 2 'ka (h-1) 50 FIX'; Results: fixed to the lowest value not significantly increasing OFV
    lcl <- log(27.9); label("Apparent clearance CL/F for LCIG at 70 kg (L/h)") # Table 2 'CL/F LCIG (L/h/70 kg) 27.9 (7.31)'
    lvc <- log(74.5); label("Apparent volume of distribution V/F at 70 kg (L)") # Table 2 'VC/F (L/70 kg) 74.5 (7.60)' (abstract prints 74.4)
    lfdepot <- fixed(log(1)); label("Relative bioavailability of LCIG/LECIG intestinal gel (unitless)") # Table 2 'Frel,LCIG/LECIG 1 FIX'
    lktr <- fixed(log(2.4)); label("Transfer rate constant of the oral-tablet absorption chain (1/h)") # Table 2 'ktr oral (h-1) 2.4 FIX'; Methods: from Othman 2014, one transit compartment
    lfdepot_oral <- fixed(log(1.03)); label("Bioavailability of oral levodopa-carbidopa tablets relative to the gel (unitless)") # Table 2 'Frel,oral 1.03 FIX'
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Methods 'Model development': allometric exponent 0.75 for CL/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Methods 'Model development': allometric exponent 1 for V/F
    e_conmed_entacapone_cl <- -0.365; label("Fractional shift in CL/F during LECIG (entacapone) infusion (unitless)") # Table 2 'CL/F LECIG,Shift -0.365 (5.24)'

    etalcl ~ 0.07497 # Table 2 'IIV CL/F,LCIG 27.9 (19.8)' CV%; log(0.279^2 + 1)
    etae_conmed_entacapone_cl ~ 0.01291 # Table 2 'IIV CL/F,LECIG,Shift 11.4 (23.5)' CV%; log(0.114^2 + 1); exponential on the shift, footnote a
    etalvc ~ 0.1118 # Table 2 'IIV VC 34.4 (17.0)' CV%; log(0.344^2 + 1)

    propSd <- 0.110; label("Proportional residual error (fraction)") # Table 2 'Proportional error (%) 11.0 (27.4)'
    addSd <- 0.316; label("Additive residual error (ug/mL)") # Table 2 'Additive error (ug/mL) 0.316 (10.2)'
  })

  model({
    ka <- exp(lka)
    ktr <- exp(lktr)
    shift_cl <- e_conmed_entacapone_cl * exp(etae_conmed_entacapone_cl)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + shift_cl * CONMED_ENTACAPONE)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(depot_oral) <- -ktr * depot_oral
    d/dt(transit1) <- ktr * depot_oral - ktr * transit1
    d/dt(central) <- ka * depot + ktr * transit1 - kel * central

    f(depot) <- exp(lfdepot)
    f(depot_oral) <- exp(lfdepot_oral)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
