Sasaki_2022_delamanid <- function() {
  description <- "Joint population PK model for oral delamanid and its major metabolite DM-6705 in children and adolescents (0.67-17 years) with multidrug-resistant tuberculosis: two-compartment delamanid disposition with a three-transit-compartment absorption chain parameterised by the mean absorption time, first-order formation of a two-compartment DM-6705 in proportion to the fraction metabolised, allometric body-weight scaling (reference 33.5 kg), dispersible-tablet effects on the mean absorption time and bioavailability, a linear age effect on bioavailability below 2 years, a dose-level (<= 50 mg) effect on bioavailability, a linear age effect on the fraction metabolised below 6 years, and inter-occasion variability on bioavailability and mean absorption time."
  reference <- paste(
    "Sasaki T, Svensson EM, Wang X, Wang Y, Hafkin J, Karlsson MO, Mallikaarjun S.",
    "Population Pharmacokinetic and Concentration-QTc Analysis of Delamanid in",
    "Pediatric Participants with Multidrug-Resistant Tuberculosis.",
    "Antimicrob Agents Chemother. 2022;66(2):e01608-21.",
    "doi:10.1128/aac.01608-21"
  )
  vignette <- "Sasaki_2022_delamanid"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (Cc = delamanid; Cc_dm6705 = DM-6705)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of all clearances (exponent 0.75, fixed) and all volumes (exponent 1, fixed) of both delamanid and DM-6705, around a reference weight of 33.5 kg (Sasaki 2022 Table 2 footnote c). Cohort mean (SD) 19.2 (11.6) kg (Table 1).",
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Two piecewise-linear effects. Bioavailability: F1 multiplier 1 - 0.201 * (2 - AGE) when AGE < 2 years, 1 otherwise. Fraction metabolised to DM-6705: FM multiplier 1 - 0.0654 * (6 - AGE) when AGE < 6 years, 1 otherwise (Sasaki 2022 Table 2 and footnote c). The footnote prints the F1 condition as 'if age is > 2.0 years', which contradicts the Table 2 row 'Age on F1: linear slope below 2 yrs', the Results and the Discussion, and which would make F1 grow without bound with age; the below-2-years reading is used (see the vignette).",
      source_name = "age"
    ),
    FORM_DELAMANID_DT = list(
      description = "Delamanid formulation indicator: 1 = pediatric dispersible tablet (5 mg or 25 mg), 0 = adult film-coated tablet (50 mg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (film-coated tablet)",
      notes = "Multiplies the mean absorption time by (1 + 0.495) and bioavailability by (1 - 0.158); both coefficients were fixed to the estimates of a population PK model that included the adult bioequivalence trial 245 (Sasaki 2022 Table 2 footnote a). In the source trials, groups 1-2 (6-17 years) took the film-coated tablet and groups 3-4 (0-5 years) the dispersible tablet (Table 1).",
      source_name = "Formulation"
    ),
    DOSE_DELAMANID_MG = list(
      description = "Amount of delamanid, in mg, given on this dose record",
      units = "mg (per dose)",
      type = "continuous",
      reference_category = "doses above 50 mg (F1 multiplier 1)",
      notes = "Per dose record. Doses of 50 mg or less have their bioavailability multiplied by (1 + 0.580); the coefficient was fixed from the adult population PK model over a 50-400 mg dose range (Sasaki 2022 Table 2 row 'Dose on F1: increase dose = 50 mg and lower' and footnote b). Must be set on the dose records, where f(depot) is evaluated.",
      source_name = "Dose"
    ),
    OCC = list(
      description = "Occasion index for the inter-occasion variability on bioavailability and mean absorption time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Sasaki 2022 Methods define an occasion only as 'the i-th dosing occasion' and do not state how occasions were delimited or how many there were. Four occasions are encoded (OCC = 1..4), each with its own eta on F1 and on the mean absorption time, all sharing one variance. Any other OCC value (for example 0) switches the inter-occasion random effects off. The indicator must be set on the dose records for the F1 term.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Statistically significant on Vc/F in the stepwise search but removed because its significance was driven by one or a few individuals (Sasaki 2022 Results and Figure S1). No coefficient is reported."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Statistically significant on F1 in the stepwise search but removed because its significance was driven by one or a few individuals (Sasaki 2022 Results and Figure S1). No coefficient is reported."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Statistically significant on CL/F in the stepwise search but removed as collinear with body weight (Sasaki 2022 Results). No coefficient is reported."
    ),
    TPROT = list(
      description = "Serum total protein",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Statistically significant on CL/F in the stepwise search but removed as collinear with albumin (Sasaki 2022 Results). No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "delamanid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "delamanid", units = "mg", specimen = "plasma", verified = TRUE),
    central_dm6705 = list(
      analyte = "DM-6705",
      units = "mg delamanid-equivalents",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_dm6705 = list(
      analyte = "DM-6705",
      units = "mg delamanid-equivalents",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 2L,
    age_range = "0.67-17 years",
    age_median = "mean (SD) 6.36 (5.19) years; median not reported",
    weight_range = "not reported; mean (SD) by age group 39 (4.59), 24.9 (6.79), 14.2 (3.2) and 9.76 (1.83) kg",
    weight_median = "mean (SD) 19.2 (11.6) kg; median not reported",
    sex_female_pct = 51.4,
    race_ethnicity = c(Asian = 67.6, Black = 5.41, Other = 27),
    disease_state = "Multidrug-resistant tuberculosis, on an optimized background regimen",
    dose_range = "Oral delamanid by age group: 100 mg BID (12-17 years) and 50 mg BID (6-11 years) as the 50-mg film-coated tablet; 25 mg BID (3-5 years) and 10 mg BID, 5 mg BID or 5 mg QD by weight (0-2 years) as the 5- or 25-mg dispersible tablet; 10 days (trial 232) followed by up to 6 months (trial 233)",
    regions = "Philippines (67.6%), South Africa (32.4%)",
    notes = "Pooled phase 1 trial 232 and its phase 2 extension trial 233 (Sasaki 2022 Methods and Table S1); 634 delamanid and 706 DM-6705 plasma concentrations. Demographics from Table 1."
  )

  ini({
    # Structural parameters of delamanid, for a 33.5-kg child
    # (Sasaki 2022 Table 2 and footnote c).
    lcl <- log(17.2)
    label("Apparent clearance of delamanid CL/F (L/h)") # Table 2 'CL/F (L/h)' = 17.2 (RSE 3.36%)
    lvc <- log(346)
    label("Apparent central volume of delamanid Vc/F (L)") # Table 2 'Vc/F (L)' = 346 (RSE 8.15%)
    lq <- log(62.4)
    label("Apparent intercompartmental clearance of delamanid Q/F (L/h)") # Table 2 'Q/F (L/h)' = 62.4 (RSE 17.6%)
    lvp <- log(296)
    label("Apparent peripheral volume of delamanid Vp/F (L)") # Table 2 'Vp/F (L)' = 296 (RSE 13.8%)
    lmtt <- log(2.73)
    label("Mean absorption time through the depot and three transit compartments, film-coated tablet (h)") # Table 2 'MAT (h)' = 2.73 (RSE 7.41%)
    lfdepot <- fixed(log(1))
    label("Typical relative bioavailability F1 (reference: > 2 years, film-coated tablet, dose > 50 mg)") # Table 2 footnote c 'F1 = 1.0{...}'

    # Structural parameters of DM-6705, apparent with respect to the
    # delamanid bioavailability and fraction metabolised; the model was
    # fitted on molar concentrations (Methods).
    lcl_dm6705 <- log(54.2)
    label("Apparent clearance of DM-6705 CLM/F (L/h)") # Table 2 'CLM/F (L/h)' = 54.2 (RSE 6.61%)
    lvc_dm6705 <- log(77.0)
    label("Apparent central volume of DM-6705 VcM/F (L)") # Table 2 'VcM/F (L)' = 77.0 (RSE 46.0%)
    lq_dm6705 <- log(425)
    label("Apparent intercompartmental clearance of DM-6705 QM/F (L/h)") # Table 2 'QM/F (L/h)' = 425 (RSE 5.34%)
    lvp_dm6705 <- log(13150)
    label("Apparent peripheral volume of DM-6705 VpM/F (L)") # Table 2 'VpM/F (L)' = 13,150 (RSE 4.27%)
    lfm <- fixed(log(1))
    label("Typical fraction of delamanid clearance forming DM-6705 (reference >= 6 years)") # Results: 'FM was fixed to 1.0 (estimation would render the model structurally unidentifiable)'

    # Covariate effects (Sasaki 2022 Table 2 and footnote c).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F, Q/F, CLM/F and QM/F (unitless)") # Table 2 'Allometric exponent for wt on CL/F, Q/F, CLM/F, and QM/F' = 0.75 Fixed
    e_wt_vc <- fixed(1.00)
    label("Allometric exponent of body weight on Vc/F, Vp/F, VcM/F and VpM/F (unitless)") # Table 2 'Allometric exponent for wt on Vc/F, Vp/F, VcM/F, and VpM/F' = 1.00 Fixed
    e_form_dt_mtt <- fixed(0.495)
    label("Fractional change in mean absorption time for the dispersible tablet (unitless)") # Table 2 'Formulation (dispersible tablet) on MAT' = 0.495 Fixed (footnote a: from the model including adult trial 245)
    e_form_dt_fdepot <- fixed(-0.158)
    label("Fractional change in F1 for the dispersible tablet (unitless)") # Table 2 'Formulation (dispersible tablet) on F1' = -0.158 Fixed (footnote a)
    e_age_fdepot <- 0.201
    label("Linear decrease in F1 per year of age below 2 years (1/year)") # Table 2 'Age on F1: linear slope below 2 yrs (1/yrs)' = 0.201 (RSE 35.5%)
    e_dose_le50_fdepot <- fixed(0.580)
    label("Fractional increase in F1 for doses of 50 mg or less (unitless)") # Table 2 'Dose on F1: increase dose = 50 mg and lower' = 0.580 Fixed (footnote b: adult model, 50-400 mg)
    e_age_fm <- 0.0654
    label("Linear decrease in the fraction metabolised per year of age below 6 years (1/year)") # Table 2 'Age on FM: linear slope below 6 yrs (1/yrs)' = 0.0654 (RSE 17.9%)

    # Inter-individual variability. Table 2 reports CV%; variances are
    # omega^2 = log(CV^2 + 1). The CL/F - CLM/F entry of 71.0 is a
    # correlation coefficient: cov = 0.710 * sqrt(0.0262222 * 0.1087836).
    etalcl + etalcl_dm6705 ~ c(
      0.0262222,
      0.0379205, 0.1087836
    ) # Table 2 'IIV on CL/F' = 16.3%, 'IIV on CLM/F' = 33.9%, 'Correlation IIV term between CL/F and CLM/F' = 71.0
    etalvp ~ 0.2943287 # Table 2 'IIV on Vp/F' = 58.5% CV
    etalfm ~ 0.0180609 # Table 2 'IIV on FM' = 13.5% CV

    # Inter-occasion variability, one eta per encoded occasion with a
    # shared variance (occasions 2-4 fixed to the occasion-1 value).
    etaiov_fdepot_1 ~ 0.0693619 # Table 2 'IOV on F1' = 26.8% CV
    etaiov_fdepot_2 ~ fixed(0.0693619) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.0693619) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(0.0693619) # same variance as occasion 1
    etaiov_mtt_1 ~ 0.3136780 # Table 2 'IOV on MAT' = 60.7% CV
    etaiov_mtt_2 ~ fixed(0.3136780) # same variance as occasion 1
    etaiov_mtt_3 ~ fixed(0.3136780) # same variance as occasion 1
    etaiov_mtt_4 ~ fixed(0.3136780) # same variance as occasion 1

    # Residual error: proportional for each analyte (Methods,
    # Y = C * (1 + eps)). The 39.4% residual correlation between the two
    # analytes (Table 2) cannot be expressed and is omitted.
    propSd <- 0.308
    label("Proportional residual error, delamanid (fraction)") # Table 2 'Proportional error, delamanid' = 30.8% CV
    propSd_dm6705 <- 0.181
    label("Proportional residual error, DM-6705 (fraction)") # Table 2 'Proportional error, DM-6705' = 18.1% CV
  })

  model({
    # Occasion indicators for the inter-occasion variability.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4

    # Allometric scaling around 33.5 kg (Table 2 footnote c).
    wt_cl <- (WT / 33.5)^e_wt_cl
    wt_v <- (WT / 33.5)^e_wt_vc

    # Delamanid disposition.
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v

    # DM-6705 disposition.
    cl_dm6705 <- exp(lcl_dm6705 + etalcl_dm6705) * wt_cl
    vc_dm6705 <- exp(lvc_dm6705) * wt_v
    q_dm6705 <- exp(lq_dm6705) * wt_cl
    vp_dm6705 <- exp(lvp_dm6705) * wt_v

    # Absorption: dose into depot, then three transit compartments, all
    # with the same rate constant; MAT spans the four first-order steps.
    mtt <- exp(lmtt + iov_mtt) * (1 + e_form_dt_mtt * FORM_DELAMANID_DT)
    ktr <- 4 / mtt

    # Bioavailability (Table 2 footnote c): linear age effect below 2 years,
    # dispersible-tablet effect and low-dose (<= 50 mg) effect.
    f_age <- 1 - e_age_fdepot * (2 - AGE) * (AGE < 2)
    f_form <- 1 + e_form_dt_fdepot * FORM_DELAMANID_DT
    f_dose <- 1 + e_dose_le50_fdepot * (DOSE_DELAMANID_MG <= 50)
    fdepot <- exp(lfdepot + iov_fdepot) * f_age * f_form * f_dose

    # Fraction of delamanid clearance forming DM-6705 (linear age effect
    # below 6 years).
    age_fm <- 1 - e_age_fm * (6 - AGE) * (AGE < 6)
    fm <- exp(lfm + etalfm) * age_fm

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_dm6705 <- cl_dm6705 / vc_dm6705
    k12_dm6705 <- q_dm6705 / vc_dm6705
    k21_dm6705 <- q_dm6705 / vp_dm6705

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_dm6705) <- fm * kel * central - kel_dm6705 * central_dm6705 - k12_dm6705 * central_dm6705 + k21_dm6705 * peripheral1_dm6705
    d/dt(peripheral1_dm6705) <- k12_dm6705 * central_dm6705 - k21_dm6705 * peripheral1_dm6705

    f(depot) <- fdepot

    # The metabolite states hold delamanid-equivalent mass (the model was
    # fitted in molar units, so formation is mole-for-mole). Molecular
    # weights from PubChem (not printed in the paper): delamanid
    # C25H25F3N4O6 534.5 g/mol (CID 6480466), DM-6705 C23H26F3N3O4
    # 465.5 g/mol (CID 44511499).
    mw_ratio_dm6705 <- 465.5 / 534.5

    Cc <- central / vc * 1000
    Cc_dm6705 <- central_dm6705 / vc_dm6705 * 1000 * mw_ratio_dm6705

    Cc ~ prop(propSd)
    Cc_dm6705 ~ prop(propSd_dm6705)
  })
}
