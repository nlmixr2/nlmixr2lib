Carmichael_2003_hydroxychloroquine <- function() {
  description <- "One-compartment first-order-absorption population PK model with an absorption lag time for oral hydroxychloroquine (HCQ) whole-blood concentrations in 123 adult rheumatoid arthritis patients (74 on HCQ alone plus 49 on HCQ + methotrexate) pooled from four Australian studies, with bioavailability fixed at the value 0.746 estimated from a nine-patient IV/oral crossover sub-study and a linear additive shift in central volume of distribution for concomitant methotrexate coadministration (V_MTX = 1070 L added to the base V of 605 L when MTX is present) (Carmichael 2003)."
  reference   <- "Carmichael SJ, Charles B, Tett SE. Population Pharmacokinetics of Hydroxychloroquine in Patients With Rheumatoid Arthritis. Ther Drug Monit. 2003;25(6):671-681. doi:10.1097/00007691-200312000-00005. PMID 14639053."
  vignette    <- "Carmichael_2003_hydroxychloroquine"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    CONMED_MTX = list(
      description        = "Concomitant methotrexate coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCQ single-agent therapy; no concomitant methotrexate)",
      notes              = "Time-fixed per subject within the analysis window (patients are stratified into the HCQ-alone cohort of 74 patients or the HCQ + MTX cohort of 49 patients per Materials and Methods 'Patients'). Enters as a linear additive shift on the typical central volume of distribution per equation 3: V_i = V + V_MTX * MTX where V_MTX = 1070 L (Results 'Population Pharmacokinetic Model Including Patients Taking MTX' paragraph 3 and Table 2 row d). The paper's Discussion attributes the 175% higher V in the HCQ + MTX cohort partly to sparser early samples in that group rather than to a true drug-drug interaction; clearance was unchanged.",
      source_name        = "MTX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 123L,
    n_subjects_hcq_alone = 74L,
    n_subjects_hcq_mtx   = 49L,
    n_studies      = 4L,
    n_observations = 780L,
    age_range      = "20-81 years",
    weight_range   = "44-89 kg",
    sex_female_pct = 71.5,
    race_ethnicity = "not reported",
    disease_state  = "Adults with rheumatoid arthritis (RA), diagnosed by 1958 American Rheumatoid Association criteria or 1987 American College of Rheumatology criteria. HCQ was used as second-line therapy. 49 patients received concomitant methotrexate.",
    dose_range     = "155 mg (n=44) or 310 mg (n=21) oral HCQ base every 24 hours in the single-agent cohort, and 155 mg (n=25) or 310 mg (n=24) oral HCQ base every 24 hours in the HCQ + MTX cohort. Nine bioavailability-study patients received a single 155-mg oral dose and, in a randomized crossover, a 30-minute IV infusion of 155 mg HCQ. Doses of 155 and 310 mg base correspond to 200 mg and 400 mg HCQ sulfate (Plaquenil tablets, Winthrop Laboratories).",
    regions        = "Australia (all studies coordinated from the Department of Clinical Pharmacology, St. Vincent's Hospital, Darlinghurst, Sydney, NSW).",
    notes          = "Demographics summarised from Table 1: bioavailability sub-study (n=9, 2 male / 7 female, age 50.2 +/- 16.1 yr, weight 72.6 +/- 13.3 kg, 590 samples), HCQ single-agent cohort (n=74, 23 male / 51 female, age 53.9 +/- 14.8 yr, weight 69.2 +/- 12.4 kg, 461 samples), HCQ + MTX cohort (n=49, 12 male / 37 female, age 56.8 +/- 12.6 yr, weight 71.0 +/- 16.3 kg, 319 samples). Sampling schedules per subgroup: (a) bioavailability sub-study patients had samples every 15 minutes for 8 hours plus 24-h and 32-h samples after single-dose HCQ; (b) 22 single-agent HCQ patients had 1-6 trough samples collected between 22 and 192 days after therapy start; (c) 43 single-agent HCQ patients had a single trough sample at 180 days; (d) HCQ + MTX patients had a single trough sample at each of 7 hospital visits over 6 months. All whole-blood HCQ concentrations were assayed by HPLC (Tett et al reference 20). Analyses used NONMEM v5.1 level 1.1 with first-order conditional estimation (FOCE) including eta-epsilon interaction. Race / ethnicity was not reported in the paper."
  )

  ini({
    # Structural parameters from Table 2 row (d) 'All 123 patients' final model.
    # Abstract: "Cl = 9.9 +/- 0.4 L/h, V = 605 +/- 91 L, ka = 0.77 +/- 0.22 hours -1, t tag = 0.44 +/- 0.02 hours."
    lcl <- log(9.89) ; label("Population apparent oral clearance (CL/F, L/h)")               # Table 2 row (d): Cl = 9.89 (SE 0.405) L/h
    lvc <- log(605)  ; label("Population apparent central volume without concomitant MTX (V/F, L)") # Table 2 row (d): V = 605 (SE 91.3) L
    lka <- log(0.765); label("First-order absorption rate constant (ka, 1/h)")               # Table 2 row (d): ka = 0.765 (SE 0.218) 1/h
    ltlag <- log(0.44) ; label("Absorption lag time (t_lag, h)")                              # Table 2 row (d): t_tag = 0.44 (SE 0.022) h
    # F was fixed to the 0.746 typical value estimated from the 9-patient bioavailability sub-study
    # and held fixed during all downstream population fits (Results 'Bioavailability' final paragraph:
    # "The typical value for F was 0.746" and Results 'Population Pharmacokinetic Model for Single-Agent
    # HCQ' paragraph 1: "bioavailability was fixed to 0.746, as estimated previously.").
    lfdepot <- fixed(log(0.746)) ; label("Oral bioavailability (F, bioavailability-sub-study estimate)") # Table 2 row (a): F = 0.746 (SE 0.068); fixed in row (d)

    # Additive shift in V per equation 3: V_i = V + V_MTX * MTX (linear, additive; L units).
    # See vignette Assumptions and deviations -- unusual coding relative to the more common
    # multiplicative "1 + e * X" convention; the paper's parameterization is preserved so that
    # the reported V_MTX standard error remains directly interpretable.
    e_conmed_mtx_vc <- 1070 ; label("Additional apparent volume when concomitant MTX is present (L, linear additive)") # Table 2 row (d): V_MTX = 1070 (SE 186) L

    # IIV variances (Table 2 row (d) omega-squared columns).
    # The Carmichael 2003 tables report omega-squared consistent with NONMEM's log-normal
    # ETA variance parameterisation (theta_i = theta * exp(eta_i), eta ~ N(0, omega^2)).
    etalcl ~ 0.127                                                                            # Table 2 row (d): omega^2_Cl = 0.127 (SE 0.023)
    etalvc ~ 0.25                                                                             # Table 2 row (d): omega^2_V  = 0.25  (SE 0.08)
    etalka ~ 0.94                                                                             # Table 2 row (d): omega^2_ka = 0.94  (SE 0.49)
    # IIV on t_lag, F, and the V_MTX shift were reported as "not estimated" ("ne") in Table 2 row (d)
    # and are omitted here per the paper's final model structure.

    # Combination proportional + additive residual error (Data Analysis paragraph and equation 2).
    # Table 2 row (d): sigma_1^2 = 0.044 (SE 0.007), sigma_2^2 = 365 (SE 140).
    # Results text confirms: "the proportional component was 21% (CV) and the standard deviation
    # of the additive component was 19 ug/L".
    propSd <- 0.21 ; label("Proportional residual error (fraction)")                          # Table 2 row (d): sqrt(0.044) = 0.21 (21% CV) matches Results text
    addSd  <- 19.1 ; label("Additive residual error standard deviation (ug/L)")               # Table 2 row (d): sqrt(365) = 19.1 ug/L matches Results text ("19 ug/L")
  })

  model({
    # Individual PK parameters.
    # V includes the linear additive shift for concomitant MTX (equation 3); the log-normal ETA
    # multiplies the typical V (matching NONMEM's V_i = (V + V_MTX * MTX) * exp(eta_V) form).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- (exp(lvc) + e_conmed_mtx_vc * CONMED_MTX) * exp(etalvc)

    # Elimination micro-constant (one-compartment).
    kel <- cl / vc

    # ODE system (one-compartment first-order oral absorption).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag time (applied to the depot). tlag units are hours; dosing time unit is hours.
    alag(depot) <- exp(ltlag)

    # Bioavailability (fixed anchor from the bioavailability sub-study).
    f(depot) <- exp(lfdepot)

    # Observation: whole-blood HCQ concentration in ug/L.
    # central holds mg (dose is mg); central/vc yields mg/L; multiply by 1000 to obtain ug/L
    # consistent with the paper's Figure 1-4 axes and Table 1 concentrations.
    Cc <- (central / vc) * 1000

    # Combination proportional + additive residual error (equation 2).
    Cc ~ add(addSd) + prop(propSd)
  })
}
