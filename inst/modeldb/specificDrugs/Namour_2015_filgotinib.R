Namour_2015_filgotinib <- function() {
  description <- "Joint population PK/PD model for the selective JAK1 inhibitor filgotinib (GLPG0634) and its active metabolite in healthy male volunteers and patients with rheumatoid arthritis: two-compartment filgotinib PK with first-order absorption feeding a one-compartment metabolite (with an empirical dose-dependent loss of metabolite formation at 300 and 450 mg once daily), and a direct-response competitive sigmoid inhibition of IL-6-induced STAT1 phosphorylation (pSTAT1) in CD4+ lymphocytes driven by both filgotinib and metabolite plasma concentrations."
  reference <- paste(
    "Namour F, Diderichsen PM, Cox E, Vayssiere B, Van der Aa A, Tasset C,",
    "Van't Klooster G. Pharmacokinetics and Pharmacokinetic/Pharmacodynamic",
    "Modeling of Filgotinib (GLPG0634), a Selective JAK1 Inhibitor, in Support",
    "of Phase IIB Dose Selection. Clin Pharmacokinet. 2015;54(8):859-874.",
    "doi:10.1007/s40262-015-0240-z",
    sep = " "
  )
  vignette <- "Namour_2015_filgotinib"

  # Dose in mg of filgotinib, volumes in L, so central / vc is mg/L; the
  # observation lines multiply by 1000 to give ng/mL, the unit of the paper's
  # assay, of Table 5's IC50 values and of the PD driver.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on apparent filgotinib clearance CL_P/F with exponent",
        "0.679 (Namour 2015 Table 4 'Effect of weight'). Namour 2015 states",
        "continuous covariates entered 'as centered power functions' (PsN scm",
        "centres on the dataset median) but does not report the centring weight.",
        "75 kg is used: it is the body weight of the 'typical subject (male",
        "healthy volunteer of 75 kg)' in which the paper ran all of its",
        "typical-value simulations (Fig. 8, Tables 6-8), so at WT = 75 kg the",
        "model returns the published typical CL_P/F exactly. See vignette",
        "Assumptions and deviations.",
        sep = " "
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Factor on apparent filgotinib central volume V_C/F:",
        "V_C/F * (1 + 2.95 * SEXF) (Namour 2015 Table 4 'Effect of sex'). The",
        "paper entered binary covariates 'as factors' with PsN scm, whose",
        "categorical form is (1 + theta) for the non-reference level -- the same",
        "form the paper states explicitly for the RA bioavailability effect",
        "('Relative increase in bioavailability'). Male is the reference: all 56",
        "phase I volunteers were male. Namour 2015 Discussion notes that sex is",
        "confounded with study and subject status (33/36 RA patients were",
        "female), so this effect may represent a study or subject effect.",
        sep = " "
      ),
      source_name = "SEX"
    ),
    DIS_RA = list(
      description = "Rheumatoid arthritis patient indicator (1 = RA patient, 0 = healthy volunteer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = paste(
        "Relative bioavailability of filgotinib F = 1 + 0.216 * DIS_RA (Namour",
        "2015 Table 4 'Relative increase in bioavailability in RA patients vs",
        "healthy volunteers'; the paper's 'patient status' covariate). Healthy",
        "volunteers are the reference with F fixed to 1, so all PK parameters",
        "are apparent (/F) values for healthy volunteers.",
        sep = " "
      ),
      source_name = "patient status"
    ),
    DOSE_FILGOTINIB_MGD = list(
      description = "Total daily filgotinib dose of the subject's regimen",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the fraction of filgotinib clearance forming the active",
        "metabolite. Namour 2015 included a 'secondary pathway of elimination'",
        "only 'at doses at and above 300 mg' and estimated it separately for",
        "the 300 mg q.d. and 450 mg q.d. cohorts (Table 4 FRAC300mg,",
        "FRAC450mg); every lower regimen, including 100 mg b.i.d. (200 mg/day)",
        "and the RA patients' 200 mg/day, has complete conversion. Encoded as",
        "fm = 1 below 300 mg/day, the 300 mg value for 300 <= dose < 450, and",
        "the 450 mg value from 450 mg/day. The paper studied only these dose",
        "levels; values between them are not supported by the source.",
        sep = " "
      ),
      source_name = "dose group"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Entered the full covariate model on the metabolite volume of distribution but was removed by backward elimination (p = 0.001); Namour 2015 Results, Population Pharmacokinetic Model."
    )
  )

  # The metabolite states hold filgotinib-MASS equivalents: Namour 2015 Eq. 1
  # forms the metabolite at FRAC * CL_P/F * C_P with no molecular-weight
  # conversion and fits it against metabolite concentrations in ng/mL, so the
  # metabolite-to-filgotinib molecular-weight ratio is absorbed into the
  # estimated CL_M/F and V_M/F. Cc_metab is therefore directly comparable with
  # the measured metabolite ng/mL. The paper does not name the metabolite.
  compartmentData <- list(
    depot = list(analyte = "filgotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "filgotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "filgotinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_metab = list(analyte = "filgotinib active metabolite", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 3L,
    age_range = "40-60 years (healthy volunteers); 18-70 years (RA patients, inclusion range)",
    weight_range = "not reported (healthy volunteers: BMI 18-30 kg/m^2)",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "56 healthy male volunteers from two phase I studies (NCT01179581,",
      "NCT01419990) and 24 patients with moderate-to-severe rheumatoid",
      "arthritis with inadequate response to methotrexate from a phase IIa",
      "proof-of-concept study (NCT01384422; 33/36 enrolled patients female).",
      "The pSTAT1 PK/PD component was fitted to the phase I healthy-volunteer",
      "data only.",
      sep = " "
    ),
    dose_range = paste(
      "Healthy volunteers: single oral doses 25-200 mg and 10-day multiple",
      "doses of 25, 50, 100 mg b.i.d. and 200, 300, 450 mg q.d. (capsules,",
      "fed). RA patients: 100 mg b.i.d. or 200 mg q.d. for 4 weeks. The 10 mg",
      "single dose was excluded from the population PK analysis.",
      sep = " "
    ),
    regions = "Belgium (phase I); Republic of Moldova (phase IIa)",
    notes = paste(
      "Study designs in Namour 2015 Table 1 and Methods 'Population",
      "Pharmacokinetic and Pharmacokinetic/Pharmacodynamic Modeling'. The",
      "paper does not tabulate the weight distribution or the sex split of",
      "the 24 RA patients who contributed PK data (3/36 of the full phase IIa",
      "cohort were male). Below-LLOQ samples (about 5% for filgotinib, <1% for",
      "the metabolite; LLOQ 1.00 ng/mL) were excluded.",
      sep = " "
    )
  )

  ini({
    # Namour 2015 Table 4 reports every footnote-'a' parameter on the natural-
    # log scale (ka = -0.733 is only meaningful as log(ka)); the back-
    # transformed values below are exp() of the tabulated numbers. The
    # metabolite IC50 cross-check holds: exp(5.68) = 293 ng/mL and
    # exp(7.43) = 1686 ng/mL are the IC50 values the paper quotes in its
    # Results text.
    lka <- -0.733; label("Log first-order absorption rate constant ka (log 1/h)") # Table 4 'ka' = -0.733 (log; RSE 5.64%) -> 0.480 1/h
    lcl <- 3.97; label("Log apparent filgotinib clearance CL_P/F at 75 kg (log L/h)") # Table 4 'CLP/F' = 3.97 (log; RSE 1.04%) -> 53.0 L/h
    lvc <- 3.08; label("Log apparent filgotinib central volume V_C/F in males (log L)") # Table 4 'VC/F' = 3.08 (log; RSE 9.77%) -> 21.8 L
    lq <- 2.02; label("Log apparent filgotinib inter-compartmental clearance Q/F (log L/h)") # Table 4 'Q/F' = 2.02 (log; RSE 7.16%) -> 7.54 L/h
    lvp <- 4.72; label("Log apparent filgotinib peripheral volume V_P/F (log L)") # Table 4 'VP/F' = 4.72 (log; RSE 3.29%) -> 112 L
    lcl_metab <- 1.04; label("Log apparent metabolite clearance CL_M/F (log L/h)") # Table 4 'CLM/F' = 1.04 (log; RSE 3.44%) -> 2.83 L/h
    lvc_metab <- 4.36; label("Log apparent metabolite volume V_M/F (log L)") # Table 4 'VM/F' = 4.36 (log; RSE 0.813%) -> 78.3 L

    # Fraction of filgotinib clearance forming the metabolite at the two
    # highest once-daily doses, held on the probit scale (Table 4 footnote b
    # 'Probit-transformed parameter'; the table reports the transformed value,
    # as it does for the log-scale rows). Eq. 1 multiplies metabolite formation
    # by FRAC, and the Table 4 legend defines FRAC as the 'fraction of
    # filgotinib metabolized to its active metabolite at the respective doses',
    # so fm = pnorm(probitfm). See vignette Assumptions and deviations for the
    # conflicting Fig. 3 / prose wording.
    probitfm_dose300 <- 0.375; label("Probit of the fraction of filgotinib clearance forming the metabolite, 300 mg/day (probit)") # Table 4 'FRAC300 mg' = 0.375 (probit; RSE 23.9%) -> fm = 0.646
    probitfm_dose450 <- 0.363; label("Probit of the fraction of filgotinib clearance forming the metabolite, 450 mg/day (probit)") # Table 4 'FRAC450 mg' = 0.363 (probit; RSE 38.5%) -> fm = 0.642

    # Covariate effects (final model after backward elimination)
    e_wt_cl <- 0.679; label("Power exponent of body weight on CL_P/F (unitless)") # Table 4 'Effect of weight' = 0.679 (RSE 32.5%)
    e_sexf_vc <- 2.95; label("Fractional increase in V_C/F for females vs males (unitless)") # Table 4 'Effect of sex' = 2.95 (RSE 47.7%)
    e_dis_ra_fdepot <- 0.216; label("Fractional increase in relative bioavailability in RA patients vs healthy volunteers (unitless)") # Table 4 'Relative increase in bioavailability in RA patients vs healthy volunteers' = 0.216 (RSE 39.0%)

    # pSTAT1 direct-response model (Table 5; Eq. 2)
    rbase <- 22.3; label("Typical baseline percentage of pSTAT1-positive CD4+ cells after IL-6 stimulation (%)") # Table 5 'pSTAT1BL' = 22.3 (RSE 9.283%)
    lic50 <- 5.68; label("Log filgotinib plasma concentration for half-maximal pSTAT1 inhibition IC50_P (log ng/mL)") # Table 5 'log IC50,P' = 5.68 (RSE 4.25%) -> 293 ng/mL
    lic50_metab <- 7.43; label("Log metabolite plasma concentration for half-maximal pSTAT1 inhibition IC50_M (log ng/mL)") # Table 5 'log IC50,M' = 7.43 (RSE 1.45%) -> 1686 ng/mL
    lhill <- 0.680; label("Log Hill coefficient shared by filgotinib and metabolite (log unitless)") # Table 5 'log H' = 0.680 (RSE 27.7%) -> 1.97

    # Between-subject variability. Table 4 BSV is log-normal (variances on the
    # log scale); no covariances were reported.
    etalcl ~ 0.102 # Table 4 'CLP/F' BSV variance 0.102 (RSE 19.9%)
    etalvc ~ 2.55 # Table 4 'VC/F' BSV variance 2.55 (RSE 12.5%)
    etalcl_metab ~ 0.0444 # Table 4 'CLM/F' BSV variance 0.0444 (RSE 9.50%)
    etalvc_metab ~ 0.0418 # Table 4 'VM/F' BSV variance 0.0418 (RSE 13.6%)
    # Table 5: baseline BSV is normally distributed (additive, %^2); IC50 BSV
    # is log-normal.
    etarbase ~ 117 # Table 5 'pSTAT1BL' BSV variance 117 (RSE 12.1%), footnote a 'Normally distributed'
    etalic50 ~ 1.87 # Table 5 'log IC50,P' BSV variance 1.87 (RSE 22.3%), footnote b 'Log-normally'
    etalic50_metab ~ 0.106 # Table 5 'log IC50,M' BSV variance 0.106 (RSE 49.5%), footnote b 'Log-normally'

    # Residual error. Eq. 1 and Table 4: proportional error on both analytes,
    # jointly multivariate normal; the covariance 0.0684 (correlation 0.437)
    # between the two proportional errors cannot be expressed across separate
    # nlmixr2 endpoints and is dropped (see vignette).
    propSd <- 0.5805; label("Proportional residual SD, filgotinib (fraction)") # Table 4 sigma2_P,prop = 0.337 (RSE 4.79%) -> sqrt = 0.5805
    propSd_metab <- 0.2694; label("Proportional residual SD, metabolite (fraction)") # Table 4 sigma2_M,prop = 0.0726 (RSE 8.55%) -> sqrt = 0.2694
    addSd_pstat1 <- 1.4526; label("Additive residual SD, pSTAT1 (% positive cells)") # Table 5 sigma2_pSTAT1,add = 2.11 (RSE 27.1%) -> sqrt = 1.4526
    propSd_pstat1 <- 0.4450; label("Proportional residual SD, pSTAT1 (fraction)") # Table 5 sigma2_pSTAT1,prop = 0.198 (RSE 13.1%) -> sqrt = 0.4450
  })

  model({
    # 1. Individual PK parameters. Weight centred on 75 kg (see covariateData).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc) * (1 + e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp)
    cl_metab <- exp(lcl_metab + etalcl_metab)
    vc_metab <- exp(lvc_metab + etalvc_metab)

    # 2. Fraction of filgotinib clearance forming the metabolite: complete
    #    below 300 mg/day; the probit-scale cohort estimates at 300 and 450
    #    mg/day (Table 4 FRAC300mg / FRAC450mg; Eq. 1). The remaining
    #    (1 - fm) of CL_P/F is the paper's empirical secondary pathway.
    fm <- 1
    if (DOSE_FILGOTINIB_MGD >= 300) fm <- pnorm(probitfm_dose300)
    if (DOSE_FILGOTINIB_MGD >= 450) fm <- pnorm(probitfm_dose450)

    # 3. ODEs (Namour 2015 Eq. 1 and Fig. 3), amounts in mg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - q * (central / vc - peripheral1 / vp) - cl * central / vc
    d/dt(peripheral1) <- q * (central / vc - peripheral1 / vp)
    d/dt(central_metab) <- fm * cl * central / vc - cl_metab * central_metab / vc_metab

    # 4. Relative bioavailability, 1 in healthy volunteers.
    f(depot) <- 1 + e_dis_ra_fdepot * DIS_RA

    # 5. Plasma concentrations (ng/mL).
    Cc <- 1000 * central / vc
    Cc_metab <- 1000 * central_metab / vc_metab

    # 6. pSTAT1 (Eq. 2): competitive sigmoid inhibition with complete maximal
    #    inhibition, each moiety's concentration normalised to its own IC50
    #    and raised to the shared Hill coefficient.
    rbase_i <- rbase + etarbase
    ic50 <- exp(lic50 + etalic50)
    ic50_metab <- exp(lic50_metab + etalic50_metab)
    hill <- exp(lhill)
    c0 <- (Cc / ic50)^hill
    c0_metab <- (Cc_metab / ic50_metab)^hill
    pstat1 <- rbase_i * (1 - (c0 + c0_metab) / (1 + c0 + c0_metab))

    Cc ~ prop(propSd)
    Cc_metab ~ prop(propSd_metab)
    pstat1 ~ add(addSd_pstat1) + prop(propSd_pstat1)
  })
}
