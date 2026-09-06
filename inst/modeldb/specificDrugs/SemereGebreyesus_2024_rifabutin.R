SemereGebreyesus_2024_rifabutin <- function() {
  description <- "Joint parent-metabolite population pharmacokinetic model for rifabutin and its 25-O-desacetyl metabolite (des-rifabutin) in 28 HIV/TB co-infected Nigerian children aged 8 months to 15 years on lopinavir/ritonavir-based antiretroviral therapy (Semere Gebreyesus 2024). Two-compartment disposition for both parent and metabolite with first-order absorption and an absorption lag time. Rifabutin elimination is split into two parallel parent pathways: an inhibitable CYP3A4 pathway (fully switched off by lopinavir/ritonavir co-treatment) and a clearance conversion via arylacetamide deacetylase (AADAC) that generates des-rifabutin at 1:1 molar stoichiometry. Lopinavir/ritonavir co-treatment raises rifabutin bioavailability 2.58-fold and cuts des-rifabutin clearance by 76.6 percent. Severe underweight lowers bioavailability by 26.0 percent per weight-for-age z-score unit below -3, and children aged 3 years or younger absorb 72.3 percent more slowly. Body weight is allometrically scaled a priori on all clearances and volumes (exponents fixed at 0.75 and 1) normalized to the 10 kg cohort median."
  reference <- "Semere Gebreyesus M, Wasmann RE, McIlleron H, Oladokun R, Okonkwo P, Wiesner L, Denti P, Rawizza HE. Population pharmacokinetics of rifabutin among HIV/TB co-infected children on lopinavir/ritonavir-based antiretroviral therapy. Antimicrob Agents Chemother. 2024;68(8):e00354-24. doi:10.1128/aac.00354-24"
  vignette <- "SemereGebreyesus_2024_rifabutin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(
      analyte = "Rifabutin", units = "mg", specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "Rifabutin", units = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "Rifabutin", units = "mg", specimen = "plasma", verified = TRUE
    ),
    central_desacetylrbn = list(
      analyte = "25-O-desacetyl rifabutin", units = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral1_desacetylrbn = list(
      analyte = "25-O-desacetyl rifabutin", units = "mg", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling applied a priori to ALL clearance and volume",
        "parameters of BOTH parent and metabolite, with exponents fixed at",
        "0.75 (clearances) and 1.0 (volumes) and normalization to the 10 kg",
        "cohort median weight (Semere Gebreyesus 2024 Methods 'Population",
        "pharmacokinetic analysis' and Results 'Structural model'; Table 2",
        "footnote c 'Allometric scaling with total body weight. Values are",
        "reported for median weight of 10 kg.' -- footnote c is attached to",
        "the des-rifabutin rows as well as the rifabutin rows, so the",
        "metabolite parameters ARE weight-scaled here. This differs from",
        "Hennig_2015_rifabutin.R, where the metabolite parameters were NOT",
        "scaled.). Cohort median 11 kg (range 4.5-45.0 kg; Table 1).",
        "Baseline in the source analysis (Table 1 demographics are at first",
        "pharmacokinetic sampling, week 2)."
      ),
      source_name        = "WT"
    ),
    WAZ = list(
      description        = "Weight-for-age z-score",
      units              = "unitless (z-score; standard-deviation units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column name ZWFA. Calculated with WHO growth charts for",
        "children aged 10 years or younger and US CDC growth charts for those",
        "over 10 years (Semere Gebreyesus 2024 Methods 'Population",
        "pharmacokinetic analysis'). Drives a THRESHOLD effect on rifabutin",
        "bioavailability that is active only below the WHO severe-underweight",
        "cut-off of -3: bioavailability falls 26.0 percent for each z-score",
        "unit below -3, and is unaffected at WAZ >= -3 (Table 2 row 'ZWFA",
        "effect (each point below -3) on bioavailability' = -26.0 percent",
        "with footnote d 'ZWFA effect on bioavailability per unit decrease in",
        "children with a ZWFA of <-3'; Results 'Covariate model').",
        "Cohort median -3.33 (range -5.15 to -1.32); about 60 percent of the",
        "children were severely underweight (WAZ < -3), Results 'Participant",
        "characteristics'. Baseline (time-fixed) in the source analysis.",
        "NOTE: the paper reports the coefficient but NOT the functional form",
        "of the below-threshold effect; the compounding (power) form encoded",
        "here is an extraction-side reading -- see the model-file comment on",
        "e_waz_fdepot and the vignette 'Assumptions and deviations' section."
      ),
      source_name        = "ZWFA"
    ),
    CONMED_LPV = list(
      description        = "Concomitant lopinavir/ritonavir (LPV/r)-based antiretroviral therapy indicator (1 = on LPV/r, 0 = rifabutin-containing TB treatment alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant lopinavir/ritonavir)",
      notes              = paste(
        "TIME-VARYING within subject in the source design: the under-1-year",
        "and 1-to-3-year cohorts received 2 weeks of rifabutin-containing TB",
        "treatment alone and then started LPV/r-based ART, so the same child",
        "contributes both CONMED_LPV = 0 and CONMED_LPV = 1 occasions",
        "(Methods 'Drug administration and sampling'). The 3-to-15-year",
        "cohort was ART-experienced and contributes CONMED_LPV = 1 only.",
        "Carries THREE independent effects (Table 2 'Covariates' block):",
        "(1) rifabutin bioavailability +158 percent, i.e. a 2.58-fold",
        "increase (Results 'Covariate model': 'increased rifabutin",
        "bioavailability by 2.58-fold (1.93-3.46)'); (2) the inhibitable",
        "CYP3A4 rifabutin clearance pathway is FIXED to -100 percent, i.e.",
        "completely switched off; (3) des-rifabutin clearance -76.6 percent.",
        "Ritonavir is the CYP3A4-inhibiting booster; the register canonical",
        "CONMED_LPV is the general 'on lopinavir at all' flag and is the",
        "same column used by Hoglund_2015_lumefantrine.R,",
        "Hoglund_2015_artemether.R and Kay_2020_lumefantrine.R."
      ),
      source_name        = "LPV/r co-treatment"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used ONLY through the paper's own dichotomization at 3 years:",
        "children aged 3 years or younger absorb 72.3 percent more slowly",
        "(Table 2 row 'Age effect (<=3 years old) on Ka' = -72.3 percent;",
        "Results 'Covariate model'). Decomposed inside model() into the",
        "binary indicator `age_le3 <- (AGE <= 3)`, following the same",
        "derive-the-indicator-in-model() idiom the register prescribes for",
        "OCC. Cohort median age 10 years (range 0.67-15.0 years; Table 1).",
        "Maturation was tested on all clearance parameters and was NOT",
        "retained (dOFV = -1.52, df = 2, P > 0.05; complete maturation was",
        "estimated at 8 months, the youngest age in the data set), so AGE",
        "carries no clearance effect in this model."
      ),
      source_name        = "age group"
    ),
    OCC = list(
      description        = "Integer-valued pharmacokinetic sampling-occasion indicator (1..4)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Between-occasion variability (BOV) was estimated for all three",
        "absorption parameters -- bioavailability, absorption rate constant",
        "and absorption lag time (Semere Gebreyesus 2024 Results 'Structural",
        "model': 'Between-occasion variability was estimated for all",
        "absorption parameters'; Table 2 BOV column). The paper reports ONE",
        "BOV magnitude per parameter shared across occasions, which is the",
        "NONMEM $OMEGA BLOCK(1) SAME pattern; nlmixr2 has no SAME shortcut,",
        "so occasions 2-4 are fix()'d to the occasion-1 variance following",
        "the registered idiom (Jonsson_2011_ethambutol.R,",
        "Blackman_2026_methotrexate.R,",
        "Mascarenhas_2015_pentadecanoic_triheptadecanoic.R -- the last of",
        "which likewise carries BOV on both an absorption parameter and",
        "relative bioavailability). The paper does not state an occasion",
        "count; FOUR occasions are encoded to span the richest per-child",
        "sampling schedule described in Methods 'Drug administration and",
        "sampling' (intensive sampling at weeks 2 and 4 for every cohort,",
        "plus week 6 in the under-1-year cohort, week 8 in the 3-to-15-year",
        "cohort, and sparse sampling at weeks 6 and 12 in the 1-to-3-year",
        "cohort). Decomposed inside model() into binary indicators oc1..oc4.",
        "For single-occasion simulations pass OCC = 1 so the first BOV eta",
        "applies."
      ),
      source_name        = "occasion (visit week)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28,
    n_studies      = 3,
    age_range      = "0.67-15.0 years",
    age_median     = "10 years",
    weight_range   = "4.5-45.0 kg",
    weight_median  = "11 kg",
    sex_female_pct = NULL,
    race_ethnicity = c(Black_African = 100),
    disease_state  = "HIV/tuberculosis co-infected children requiring protease-inhibitor-based antiretroviral therapy. Approximately 60 percent were severely underweight (weight-for-age z-score < -3; overall median -3.33, range -5.15 to -1.32). In the 3-to-15-year cohort 4 of 15 were WHO HIV stage 4 and the remainder stage 3. Neutropenia occurred in 12 children on 26 occasions (10 grade 1, 6 grade 2, 8 grade 3, 2 grade 4).",
    dose_range     = "Oral rifabutin suspension (20 mg/mL, compounded from Mycobutin capsules). Under-1-year cohort: 20 mg/kg/day for 2 weeks of TB-only treatment, then 5 mg/kg/day once LPV/r-based ART started. 1-to-3-year cohort: 15-20 mg/kg/day for 2 weeks of TB-only treatment, then 2.5 mg/kg/day with LPV/r. 3-to-15-year cohort (ART-experienced): 2.5 mg/kg/day with LPV/r from study entry.",
    regions        = "Nigeria (APIN PEPFAR pediatric ART program); external validation cohort from South Africa",
    cohorts        = "Three prospective age cohorts: under 1 year (n = 3), 1-3 years (n = 10), 3-15 years (n = 15).",
    sampling_design = "Intensive sampling at 0, 2, 4, 8, 12 and 24 h post dose during weeks 2 and 4 for all age groups; additionally week 6 in the under-1-year cohort and week 8 in the 3-to-15-year cohort. Sparse samples at weeks 6 and 12 (0 h and either 3-5 h or 24-26 h post dose) in the 1-to-3-year cohort. Week-2 intensive samples in the two younger cohorts were taken during rifabutin-containing TB-only treatment; all other visits were during LPV/r co-treatment. 462 samples in total, of which 16 (1.7 percent) were below the limit of quantification.",
    assay          = "LC-MS/MS quantifying rifabutin and des-rifabutin concurrently at the University of Cape Town; calibration range 3.91-1000.0 ug/L for rifabutin and 0.780-200 ug/L for des-rifabutin.",
    external_validation = "Six South African children aged 2 (0.83-3.0) years, weight 11 (9.0-12.0) kg, weight-for-age z-score -1.06 (-1.85 to 0.951), on LPV/r-based ART with rifabutin 5 mg/kg/day thrice weekly (Moultrie et al.). 36 samples, none below the limit of quantification. After validation these data were added to the analysis and the parameters re-estimated, so the Table 2 estimates encoded here reflect the pooled 28 + 6 child data set.",
    notes          = "Baseline demographics are Table 1 of Semere Gebreyesus 2024, reported at first pharmacokinetic sampling (week 2). Sex was not reported in Table 1, so sex_female_pct is NULL. Creatinine clearance (modified Schwartz) and CD4 count were tested as covariates but not retained. Fifteen profiles had pre-dose concentrations below one third of the corresponding 24 h concentration and were handled by the Dansirikul B2 initialization method (dosing history discarded, model initialized to the observed concentration); that data-handling device is a fitting-time construct and is not part of the packaged structural model."
  )

  ini({
    # =========================================================================
    # Rifabutin (parent). Two-compartment disposition with first-order
    # absorption and an absorption lag time. All clearances and volumes are
    # allometrically scaled with body weight (exponents fixed a priori at 0.75
    # and 1.0) and are reported at the 10 kg cohort median weight
    # (Table 2 footnote c).
    #
    # Rifabutin elimination from `central` runs through TWO parallel arms:
    #   * `lcl`                   -- the inhibitable CYP3A4 pathway, which
    #                                lopinavir/ritonavir switches off entirely
    #   * `lcl_form_desacetylrbn` -- the AADAC "clearance conversion" pathway
    #                                that forms des-rifabutin
    # so total rifabutin elimination clearance is (cl + cl_form_desacetylrbn).
    # Splitting the two is only identifiable because LPV/r inhibits the first
    # but not the second (Results 'Structural model').
    # =========================================================================
    lcl <- log(13.6)
    label("Rifabutin inhibitable CYP3A4-pathway clearance (L/h per 10 kg)")
    # Table 2 'Clearance (Inhibitable CYP3A4 pathway) (L/h)' = 13.6 (8.77-18.8)

    lcl_form_desacetylrbn <- log(16.2)
    label("Rifabutin-to-des-rifabutin AADAC clearance conversion (L/h per 10 kg)")
    # Table 2 'Clearance conversion (AADAC pathway) (L/h)' = 16.2 (12.8-20.7)

    lvc <- log(185)
    label("Rifabutin central volume of distribution V_C,P (L per 10 kg)")
    # Table 2 'Central volume of distribution, V C,P (L)' = 185 (135-251)

    lvp <- log(232)
    label("Rifabutin peripheral volume of distribution V_P,P (L per 10 kg)")
    # Table 2 'Peripheral volume of distribution, V P,P (L)' = 232 (171-313)

    lq <- log(25.1)
    label("Rifabutin intercompartmental clearance Q_P (L/h per 10 kg)")
    # Table 2 'Intercompartmental clearance, Q P (L/h)' = 25.1 (17.3-33.1)

    lka <- log(1.27)
    label("Rifabutin first-order absorption rate constant Ka (1/h)")
    # Table 2 'Absorption rate constant, Ka (1/h)' = 1.27 (0.810-2.14)

    ltlag <- log(0.544)
    label("Rifabutin absorption lag time (h)")
    # Table 2 'Absorption lag time, Lag (h)' = 0.544 (0.338-0.805)

    lfdepot <- fixed(log(1))
    label("Rifabutin oral bioavailability F (typical value, reference conditions)")
    # Table 2 'Bioavailability, F' = 1 fixed. Absolute F is not identifiable
    # from oral-only data; the covariate effects below are relative to it.

    # =========================================================================
    # 25-O-desacetyl rifabutin (metabolite, suffix `_desacetylrbn`).
    # Two-compartment disposition, formed only from the parent central
    # compartment via the AADAC conversion clearance. Also allometrically
    # weight-scaled (Table 2 footnote c is attached to these rows too).
    # =========================================================================
    lcl_desacetylrbn <- log(106)
    label("Des-rifabutin elimination clearance (L/h per 10 kg)")
    # Table 2 'Clearance metabolite (L/h)' = 106 (82.1-142.0)

    lvc_desacetylrbn <- log(43.0)
    label("Des-rifabutin central volume of distribution V_C,M (L per 10 kg)")
    # Table 2 'Central volume of distribution, V C,M (L)' = 43.0 (25.9-64.4)

    lvp_desacetylrbn <- log(241)
    label("Des-rifabutin peripheral volume of distribution V_P,M (L per 10 kg)")
    # Table 2 'Peripheral volume of distribution, V P,M (L)' = 241 (169-322)

    lq_desacetylrbn <- log(44.4)
    label("Des-rifabutin intercompartmental clearance Q_M (L/h per 10 kg)")
    # Table 2 'Intercompartmental clearance, Q M (L/h)' = 44.4 (31.5-63.3)

    # =========================================================================
    # Covariate effects (Table 2 'Covariates' block). All are encoded as
    # fractional changes, so the multiplier is (1 + effect) for the binary
    # covariates. The +158 percent bioavailability effect reproduces the
    # paper's own wording exactly: 1 + 1.58 = 2.58-fold.
    # =========================================================================
    e_conmed_lpv_fdepot <- 1.58
    label("Fractional change in rifabutin bioavailability on LPV/r co-treatment (unitless; +158 percent = 2.58-fold)")
    # Table 2 'LPV/r effect on bioavailability (%)' = +158 (+105 to +239);
    # Results 'Covariate model' states the 2.58-fold (1.93-3.46) increase.

    e_conmed_lpv_cl <- fixed(-1)
    label("Fractional change in the inhibitable CYP3A4 rifabutin clearance on LPV/r co-treatment (unitless; -100 percent = pathway fully inhibited)")
    # Table 2 'LPV/r effect on rifabutin clearance (inhibitable CYP3A4
    # pathway) (%)' = -100 fixed. Results 'Structural model' describes this
    # arm as "considered as fully inhibited by LPV/r".

    e_conmed_lpv_cl_desacetylrbn <- -0.766
    label("Fractional change in des-rifabutin clearance on LPV/r co-treatment (unitless; -76.6 percent)")
    # Table 2 'LPV/r effect on des-rifabutin clearance (%)' = -76.6
    # (-79.8 to -73.1); Results 'Covariate model' quotes -76.6 (74.4-78.3).

    e_waz_fdepot <- -0.26
    label("Fractional change in rifabutin bioavailability per weight-for-age z-score unit below -3 (unitless; -26.0 percent per unit)")
    # Table 2 'ZWFA effect (each point below -3) on bioavailability (%)' =
    # -26.0 (-34.1 to -18.1), footnote d "ZWFA effect on bioavailability per
    # unit decrease in children with a ZWFA of <-3"; Results 'Covariate
    # model'.
    #
    # FUNCTIONAL FORM IS AN EXTRACTION-SIDE READING. The paper reports the
    # coefficient and the -3 breakpoint but never prints the equation, and the
    # supplement (initialization text + Figures S1-S9) does not contain it
    # either. Encoded here as the COMPOUNDING (power) form
    #   fdepot *= (1 + e_waz_fdepot)^max(0, -3 - WAZ)
    # i.e. 0.74 per z-score unit below -3, because (a) that is what "a
    # decrease of 26.0% in bioavailability for each unit decrease in ZWFA"
    # says when read literally, and (b) it keeps bioavailability strictly
    # positive, whereas the linear alternative
    #   fdepot *= (1 + e_waz_fdepot * max(0, -3 - WAZ))
    # crosses zero at WAZ = -6.85 and goes negative below it -- reachable in
    # a severely malnourished African paediatric population and in the
    # paper's own 22,500-child simulation. The two forms agree at the
    # breakpoint and diverge with depth: at the cohort minimum WAZ = -5.15
    # they give 0.516 (power) versus 0.441 (linear), a 17 percent relative
    # difference. See the vignette 'Assumptions and deviations' section.

    e_age_le3_ka <- -0.723
    label("Fractional change in rifabutin absorption rate constant for children aged 3 years or younger (unitless; -72.3 percent)")
    # Table 2 'Age effect (<=3 years old) on Ka (%)' = -72.3 (-83.3 to -57.1);
    # Results 'Covariate model' quotes a 72.3 percent (48.5-85.7) slower rate
    # of absorption in the cohort aged 3 years or younger.

    # =========================================================================
    # Between-subject variability (BSV). Table 2 footnote b: BSV and BOV "were
    # obtained using sqrt(e^(OM^2) - 1) and reported as approximate % CV", so
    # the internal log-scale variance is omega^2 = log(1 + CV^2).
    #
    # Variability was deliberately sparse ("variability parameters were
    # included parsimoniously, where most relevant and essential"). Results
    # 'Structural model' is explicit that there is ONE COMMON BSV random
    # effect shared by the two rifabutin clearance pathways and a SEPARATE
    # BSV on des-rifabutin clearance -- so `etalcl` below is added to BOTH
    # `cl` and `cl_form_desacetylrbn` in model(). No BSV was estimated on any
    # volume or intercompartmental clearance; the blank Table 2 variability
    # cells on those rows are genuinely blank (confirmed against the
    # publisher PDF layout, where the merged BSV cell spans only the header
    # row) and the single reported shrinkage value per block is consistent
    # with one eta per block.
    # =========================================================================
    etalcl ~ 0.041956
    # Table 2 rifabutin block 'BSV: 20.7 (15.6-26.0)' %CV, shrinkage 14
    # percent; omega^2 = log(1 + 0.207^2) = 0.041956. Shared by the
    # inhibitable CYP3A4 clearance and the AADAC conversion clearance.

    etalcl_desacetylrbn ~ 0.090630
    # Table 2 des-rifabutin block 'BSV: 30.8 (23.2-39.2)' %CV, shrinkage 8
    # percent; omega^2 = log(1 + 0.308^2) = 0.090630.

    # =========================================================================
    # Between-occasion variability (BOV) on the three absorption parameters,
    # encoded as inter-occasion variability over four sampling occasions. One
    # magnitude per parameter is reported and shared across occasions (NONMEM
    # $OMEGA BLOCK(1) SAME); nlmixr2 has no SAME shortcut, so occasion 1
    # carries the estimated variance and occasions 2-4 fix it to that value,
    # per the registered idiom.
    # =========================================================================
    etaiov_fdepot_1 ~ 0.416709
    # Table 2 'Bioavailability, F' -> 'BOV: 71.9 (59.5-85.1)' %CV, shrinkage
    # 3 percent; omega^2 = log(1 + 0.719^2) = 0.416709 (estimated, occasion 1)
    etaiov_fdepot_2 ~ fix(0.416709)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_fdepot_3 ~ fix(0.416709)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_fdepot_4 ~ fix(0.416709)  # SAME-equivalent: equal to the occasion-1 BOV variance

    etaiov_ka_1 ~ 0.872297
    # Table 2 'Absorption rate constant, Ka (1/h)' -> 'BOV: 118 (84-175)'
    # %CV, shrinkage 11 percent; omega^2 = log(1 + 1.18^2) = 0.872297
    etaiov_ka_2 ~ fix(0.872297)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_ka_3 ~ fix(0.872297)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_ka_4 ~ fix(0.872297)  # SAME-equivalent: equal to the occasion-1 BOV variance

    etaiov_tlag_1 ~ 0.444368
    # Table 2 'Absorption lag time, Lag (h)' -> 'BOV: 74.8 (44.4-105.0)'
    # %CV, shrinkage 42 percent; omega^2 = log(1 + 0.748^2) = 0.444368
    etaiov_tlag_2 ~ fix(0.444368)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_tlag_3 ~ fix(0.444368)  # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_tlag_4 ~ fix(0.444368)  # SAME-equivalent: equal to the occasion-1 BOV variance

    # =========================================================================
    # Residual unexplained variability: combined additive + proportional on
    # each of the two outputs (Methods 'Population pharmacokinetic analysis':
    # "Combined error model (i.e., additive and proportional error) was used
    # to describe the unexplained residual variability"). Additive terms are
    # in ug/L, matching the assay and the `units$concentration` above.
    #
    # NOT ENCODED: the paper also estimated a 28.2 percent (14.6-42.1)
    # correlation coefficient between the rifabutin and des-rifabutin residual
    # error terms via the NONMEM L2 method, because the two analytes were
    # measured from the same sample. nlmixr2 has no construct for correlated
    # residual error across two outputs. Omitting it leaves every typical
    # value and every marginal residual variance unchanged and affects only
    # the joint parent/metabolite residual draw within a sample. See the
    # vignette 'Assumptions and deviations' section.
    # =========================================================================
    propSd <- 0.188
    label("Rifabutin proportional residual error (fraction)")
    # Table 2 'Rifabutin: proportional error (%)' = 18.8 (16.3-21.7)

    addSd <- 10.1
    label("Rifabutin additive residual error (ug/L)")
    # Table 2 'Rifabutin: additive error (ug/L)' = 10.1 (5.24-15.9)

    propSd_desacetylrbn <- 0.108
    label("Des-rifabutin proportional residual error (fraction)")
    # Table 2 'Des-rifabutin: proportional error (%)' = 10.8 (8.41-13.6)

    addSd_desacetylrbn <- 11.6
    label("Des-rifabutin additive residual error (ug/L)")
    # Table 2 'Des-rifabutin: additive error (ug/L)' = 11.6 (9.75-14.1)
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Derived covariate terms.
    # -----------------------------------------------------------------------
    # Allometric factors, reference 10 kg (the cohort median weight), with
    # exponents fixed a priori at 0.75 for clearances and 1.0 for volumes.
    wt_cl <- (WT / 10)^0.75
    wt_v  <- (WT / 10)

    # Lopinavir/ritonavir co-treatment multipliers. CONMED_LPV is time-varying
    # in the source design (the two younger cohorts crossed over from TB-only
    # treatment to LPV/r co-treatment mid-study).
    lpv_fdepot_factor          <- 1 + e_conmed_lpv_fdepot * CONMED_LPV
    lpv_cl_factor              <- 1 + e_conmed_lpv_cl * CONMED_LPV
    lpv_cl_desacetylrbn_factor <- 1 + e_conmed_lpv_cl_desacetylrbn * CONMED_LPV

    # Weight-for-age z-score effect on bioavailability. Active only below the
    # WHO severe-underweight cut-off of -3; `waz_below3` is the number of
    # z-score units below the breakpoint and is 0 at or above it, so children
    # with WAZ >= -3 get a factor of exactly 1.
    waz_below3    <- max(0, -3 - WAZ)
    waz_fdepot_factor <- (1 + e_waz_fdepot)^waz_below3

    # Age dichotomized at the paper's 3-year cut-off for the absorption-rate
    # effect.
    age_le3      <- (AGE <= 3)
    age_ka_factor <- 1 + e_age_le3_ka * age_le3

    # Occasion indicators multiplexing the between-occasion-variability etas
    # on the three absorption parameters. Pass OCC = 1 for single-occasion
    # simulations.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 +
      oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 +
      oc3 * etaiov_tlag_3 + oc4 * etaiov_tlag_4

    # -----------------------------------------------------------------------
    # 2. Individual parameters. `etalcl` is deliberately shared between the
    #    two rifabutin clearance pathways -- that is the paper's "common
    #    between-subject variability random effect on the two clearance
    #    pathways of rifabutin".
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wt_cl * lpv_cl_factor
    cl_form_desacetylrbn <- exp(lcl_form_desacetylrbn + etalcl) * wt_cl
    vc <- exp(lvc) * wt_v
    vp <- exp(lvp) * wt_v
    q  <- exp(lq)  * wt_cl

    ka     <- exp(lka + iov_ka) * age_ka_factor
    tlag   <- exp(ltlag + iov_tlag)
    fdepot <- exp(lfdepot + iov_fdepot) * lpv_fdepot_factor * waz_fdepot_factor

    cl_desacetylrbn <- exp(lcl_desacetylrbn + etalcl_desacetylrbn) * wt_cl *
      lpv_cl_desacetylrbn_factor
    vc_desacetylrbn <- exp(lvc_desacetylrbn) * wt_v
    vp_desacetylrbn <- exp(lvp_desacetylrbn) * wt_v
    q_desacetylrbn  <- exp(lq_desacetylrbn)  * wt_cl

    # -----------------------------------------------------------------------
    # 3. Molar stoichiometry of the parent-to-metabolite conversion.
    #    Methods 'Population pharmacokinetic analysis': "Molar conversion of
    #    the rifabutin dose and the concentrations of rifabutin and
    #    des-rifabutin was utilized to adjust for the difference in molecular
    #    weight between rifabutin (847.02 g/mol) and des-rifabutin (805
    #    g/mol)", and Results 'Structural model': "100% of the rifabutin
    #    eliminated by the clearance conversion is transformed into
    #    des-rifabutin". The ODE states here hold MASS (mg), so the 1:1 MOLAR
    #    transfer becomes a mass transfer scaled by the molecular-weight
    #    ratio 805 / 847.02.
    # -----------------------------------------------------------------------
    mwRatioDesacetylrbn <- 805 / 847.02

    # -----------------------------------------------------------------------
    # 4. ODE system (Figure S4 model schematic). Rifabutin leaves `central`
    #    by two parallel arms; the AADAC arm is the sole source of
    #    des-rifabutin.
    # -----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot

    d/dt(central) <- ka * depot -
      (cl / vc) * central -
      (cl_form_desacetylrbn / vc) * central -
      (q / vc) * central +
      (q / vp) * peripheral1

    d/dt(peripheral1) <- (q / vc) * central -
      (q / vp) * peripheral1

    d/dt(central_desacetylrbn) <-
      (cl_form_desacetylrbn / vc) * central * mwRatioDesacetylrbn -
      (cl_desacetylrbn / vc_desacetylrbn) * central_desacetylrbn -
      (q_desacetylrbn / vc_desacetylrbn) * central_desacetylrbn +
      (q_desacetylrbn / vp_desacetylrbn) * peripheral1_desacetylrbn

    d/dt(peripheral1_desacetylrbn) <-
      (q_desacetylrbn / vc_desacetylrbn) * central_desacetylrbn -
      (q_desacetylrbn / vp_desacetylrbn) * peripheral1_desacetylrbn

    # -----------------------------------------------------------------------
    # 5. Bioavailability and absorption lag on the depot.
    # -----------------------------------------------------------------------
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # -----------------------------------------------------------------------
    # 6. Observations. States are in mg and volumes in L, so amount / volume
    #    is mg/L; multiplying by 1000 gives ug/L, the unit in which the assay
    #    calibration ranges and both additive residual-error terms are
    #    reported.
    # -----------------------------------------------------------------------
    Cc                <- (central / vc) * 1000
    Cc_desacetylrbn   <- (central_desacetylrbn / vc_desacetylrbn) * 1000

    Cc              ~ add(addSd) + prop(propSd)
    Cc_desacetylrbn ~ add(addSd_desacetylrbn) + prop(propSd_desacetylrbn)
  })
}
