deAlwis_1998_ondansetron <- function() {
  description <- "Two-compartment population PK model with zero-order intravenous-infusion input for ondansetron in pooled paediatric, young-adult, elderly, and aged subjects (de Alwis 1998). The paper uses an empirical additive linear-regression covariate model in the 1990s NONMEM tradition (Maitre 1991 three-step approach): clearance CL and inter-compartmental clearance CLd are sex-stratified with separate male and female intercepts and slopes; the central volume V1 has a body-weight slope only; the steady-state volume Vss has body-weight and age slopes; the peripheral volume Vp is derived as Vss - V1. Inter-individual variability is diagonal log-normal on CL, V1, Vss, and CLd. Proportional residual error is stratified across five paper-defined study sub-populations (young healthy volunteers 18-41 y, elderly healthy volunteers 61-75 y, aged healthy volunteers >= 75 y, paediatric cancer patients receiving chemotherapy, paediatric patients receiving general anaesthesia), switched at runtime via the canonical AGE / DIS_HEALTHY / DIS_CANCER_PED covariates."
  reference <- "de Alwis DP, Aarons L, Palmer JL. Population pharmacokinetics of ondansetron: a covariate analysis. Br J Clin Pharmacol. 1998. doi:10.1046/j.1365-2125.1998.00756.x"
  vignette <- "deAlwis_1998_ondansetron"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  paper_specific_etas <- c("etalcl", "etalvc", "etalvss", "etalq")
  paper_specific_residual_sds <- c(
    "propSd_young_vol", "propSd_elderly_vol", "propSd_aged_vol",
    "propSd_paed_chemo", "propSd_paed_anaes"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters every structural parameter (CL, V1, Vss, CLd) as a linear slope in the additive linear-regression covariate model (Table 3, Three Step Full data column). Training-cohort range 10.2-95.8 kg (Methods, Data paragraph).",
      source_name        = "wt"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters CL (sex-specific slope) and Vss (single slope) in the linear-regression covariate model. Training-cohort range 2-82 y (Methods, Data paragraph). Also drives the residual-error stratification thresholds (volunteer sub-population: young AGE < 45, elderly 45 <= AGE < 75, aged AGE >= 75; chosen to separate the paper's three healthy-volunteer age bands, study 1 / study 2).",
      source_name        = "age"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. de Alwis 1998 encodes sex as 'gender 1 for females and 0 for males' (Table 2 footnote # and Stepwise regression analysis paragraph), matching the canonical SEXF orientation directly. CL and CLd carry sex-specific intercepts AND sex-specific slopes (Table 3 'CL male' / 'CL female' / 'CLd male' / 'CLd female' rows). V1 and Vss are NOT sex-stratified in the final model. Training-set sex breakdown 31 females / 68 males (Methods, Data paragraph).",
      source_name        = "gender"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy volunteer (studies 1 and 2), 0 = paediatric patient (studies 3 and 4).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (paediatric patient cohort; studies 3 and 4)",
      notes              = "Time-fixed per subject. Used together with DIS_CANCER_PED and AGE to assign the per-subject proportional-residual-error magnitude across the five paper-defined sub-populations (Table 1 footnote: *1 young volunteers, *2 elderly volunteers, *3 aged volunteers, *4 paediatric patients on chemotherapy, *5 paediatrics on anaesthesia). DIS_HEALTHY = 1 routes the subject into one of the three volunteer strata (selected by AGE); DIS_HEALTHY = 0 routes the subject into one of the two paediatric strata (selected by DIS_CANCER_PED).",
      source_name        = "(derived from study identifier; studies 1-2 = healthy volunteer, studies 3-4 = paediatric patient)"
    ),
    DIS_CANCER_PED = list(
      description        = "Paediatric oncology cohort indicator: 1 = paediatric cancer patient receiving chemotherapy (study 3), 0 = paediatric patient receiving general anaesthesia (study 4) or any non-paediatric subject (studies 1 and 2).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (paediatric general-anaesthesia cohort or non-paediatric subject)",
      notes              = "Time-fixed per subject. Discriminates the paediatric chemotherapy sub-population (study 3, n=14, AGE 2-13 y, BSA-banded i.v. doses) from the paediatric general-anaesthesia sub-population (study 4, n=19, AGE 3-11 y, prophylactic 2-4 mg i.v. before anaesthetic induction). Only used in the residual-error switch within model(); does not enter the structural covariate model.",
      source_name        = "(derived from study identifier; study 3 = paediatric chemo, study 4 = paediatric anaesthesia)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_studies      = 4L,
    age_range      = "2-82 years",
    weight_range   = "10.2-95.8 kg",
    sex_female_pct = 31.3,
    disease_state  = "Pooled healthy-volunteer + paediatric-patient cohort drawn from four phase I/II ondansetron i.v. studies. Study 1 (n=32): adult male healthy volunteers, 16 young (18-41 y) and 16 elderly (65-75 y), single 8 mg i.v. 15-min infusion. Study 2 (n=34): age-and-gender stratified healthy volunteers, young (21-38 y, 6 M + 5 F), elderly (61-74 y, 6 M + 6 F), aged (75-82 y, 5 M + 6 F), single 0.15 mg/kg i.v. 15-min infusion. Study 3 (n=14): paediatric cancer patients (2-13 y, 3 F + 11 M) receiving cancer chemotherapy, single BSA-banded i.v. 15-min infusion (5 mg/m^2 if BSA <= 1.2 m^2; 8 mg fixed otherwise). Study 4 (n=19): paediatric patients (3-11 y, 11 F + 8 M) undergoing general anaesthesia, single 5-min i.v. infusion before anaesthetic induction (2 mg if 3-7 y; 4 mg if 8-11 y).",
    dose_range     = "Single i.v. 15-min infusion (5-min in study 4): 8 mg fixed (studies 1, 5, 6), 0.15 mg/kg (study 2), 5 mg/m^2 or 8 mg BSA-banded (study 3), 2 mg or 4 mg age-banded (study 4). Repeated 0.15 mg/kg every 4 h on a single day in study 7.",
    regions        = "United Kingdom (GlaxoWellcome Research and Development sponsored studies).",
    bioanalysis    = "Plasma ondansetron quantified by validated solid-phase extraction HPLC with UV detection at 305 nm; LLOQ 1 ng/mL; linear calibration 1-250 ng/mL; CV < 10%. Plasma samples for the different studies were analysed by different centres, which is one source of the residual-error heterogeneity across the five sub-populations (Methods, Data analysis paragraph).",
    test_set       = "An independent test set of 54 subjects from studies 5-7 (young healthy male volunteers 20-35 y, n=16; adult cancer patients on chemotherapy 39-78 y, 18 M + 2 F, n=20; paediatric cancer patients 4-18 y, n=18) was used for the standardised-prediction-error validation. Weight range 16.2-100.0 kg; gender 7 F + 47 M.",
    notes          = "Total 1506 plasma concentrations from 99 model-building subjects (studies 1-4) and 607 plasma concentrations from 54 test-set subjects (studies 5-7). de Alwis 1998 Methods, Data paragraph. Demographics available for modelling were limited to age, weight, and gender; no race, no laboratory chemistry, no concomitant-medication, and no organ-function covariates were tested. The aged subgroup (study 2, AGE 75-82 y) is recorded as a distinct stratum from the elderly subgroup (AGE 61-74 y) only in the residual-error structure; the structural covariate model treats AGE as continuous."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters: empirical additive linear-regression
    # covariate model (de Alwis 1998 Table 3, Three Step / Full data
    # column; the paper's reference fitting method on the full 99-subject
    # training cohort). Sex-stratified equations for CL and CLd;
    # non-stratified equations for V1 and Vss. All clearances in L/h;
    # all volumes in L. Vp (peripheral) is derived in model() as
    # max(Vss - V1, 0.01) so the difference stays non-negative under
    # extreme covariate combinations.
    #
    # Paper equations (Table 3, Three Step Full data):
    #   CL_male   = 5.72 + 0.36 * WT - 0.17 * AGE      (L/h)
    #   CL_female = 7.78 + 0.20 * WT - 0.04 * AGE      (L/h)
    #   V1        = 8.08 + 0.31 * WT                   (L)
    #   Vss       = 4.79 + 1.88 * WT + 0.133 * AGE     (L)
    #   CLd_male   = -22.3 + 3.89 * WT                 (L/h)
    #   CLd_female = -76.9 + 7.49 * WT                 (L/h)
    #
    # The regression coefficients below use a `th_<param>_<term>` naming
    # scheme (paper-specific linear-regression intercepts and slopes that
    # do not map onto the canonical lcl / lvc / lq / e_<cov>_<param>
    # patterns); see the Assumptions and deviations section of the
    # validation vignette.

    # CL (male equation): linear regression in WT and AGE.
    th_cl_int_m <-  5.72;  label("CL intercept, males (L/h)")            # Table 3, Three Step Full data, CL male row
    th_cl_wt_m  <-  0.36;  label("CL slope on WT, males (L/h/kg)")       # Table 3, Three Step Full data, CL male row
    th_cl_age_m <- -0.17;  label("CL slope on AGE, males (L/h/year)")    # Table 3, Three Step Full data, CL male row

    # CL (female equation): linear regression in WT and AGE.
    th_cl_int_f <-  7.78;  label("CL intercept, females (L/h)")          # Table 3, Three Step Full data, CL female row
    th_cl_wt_f  <-  0.20;  label("CL slope on WT, females (L/h/kg)")     # Table 3, Three Step Full data, CL female row
    th_cl_age_f <- -0.04;  label("CL slope on AGE, females (L/h/year)")  # Table 3, Three Step Full data, CL female row

    # V1 (central volume): linear regression in WT only (no sex stratification).
    th_v1_int   <-  8.08;  label("V1 intercept (L)")                     # Table 3, Three Step Full data, V1 row
    th_v1_wt    <-  0.31;  label("V1 slope on WT (L/kg)")                # Table 3, Three Step Full data, V1 row

    # Vss (steady-state volume): linear regression in WT and AGE (no sex stratification).
    th_vss_int  <-  4.79;  label("Vss intercept (L)")                    # Table 3, Three Step Full data, Vss row
    th_vss_wt   <-  1.88;  label("Vss slope on WT (L/kg)")               # Table 3, Three Step Full data, Vss row
    th_vss_age  <-  0.133; label("Vss slope on AGE (L/year)")            # Table 3, Three Step Full data, Vss row

    # CLd (inter-compartmental clearance, male equation): linear regression in WT only.
    th_cld_int_m <- -22.3; label("CLd intercept, males (L/h)")           # Table 3, Three Step Full data, CLd male row
    th_cld_wt_m  <-   3.89; label("CLd slope on WT, males (L/h/kg)")     # Table 3, Three Step Full data, CLd male row

    # CLd (female equation): linear regression in WT only.
    th_cld_int_f <- -76.9; label("CLd intercept, females (L/h)")         # Table 3, Three Step Full data, CLd female row
    th_cld_wt_f  <-   7.49; label("CLd slope on WT, females (L/h/kg)")   # Table 3, Three Step Full data, CLd female row

    # -----------------------------------------------------------------
    # Inter-individual variability: diagonal log-normal IIV on CL, V1,
    # Vss, and CLd. Paper reports CV(%) on the natural-parameter scale
    # (multiplicative error). omega^2 = log(CV^2 + 1) with CV expressed
    # as a fraction.
    #   CL:   CV 29.6%  -> log(1 + 0.296^2) = log(1.0876) = 0.0840
    #   V1:   CV 71.1%  -> log(1 + 0.711^2) = log(1.5055) = 0.4090
    #   Vss:  CV 36.7%  -> log(1 + 0.367^2) = log(1.1347) = 0.1264
    #   CLd:  CV 37.0%  -> log(1 + 0.370^2) = log(1.1369) = 0.1283
    # All four etas (etalcl, etalvc, etalvss, etalq) are declared in
    # `paper_specific_etas` because the typical-value parameters cl_typ,
    # v1_typ, vss_typ, cld_typ are computed inside model() as additive
    # linear regressions of WT / AGE / SEXF rather than as canonical
    # lcl / lvc / lq fixed-effect parameters in ini(). The etas are
    # still log-normal multiplicative on the resulting typical values.
    etalcl  ~ 0.0840   # Table 1, Final covariate model, All data, CL CV 29.6%
    etalvc  ~ 0.4090   # Table 1, Final covariate model, All data, V1 CV 71.1% (V1 maps onto canonical vc)
    etalvss ~ 0.1264   # Table 1, Final covariate model, All data, Vss CV 36.7% (paper-specific eta on Vss; see paper_specific_etas above)
    etalq   ~ 0.1283   # Table 1, Final covariate model, All data, CLd CV 37.0% (CLd maps onto canonical q)

    # -----------------------------------------------------------------
    # Proportional residual error, split across five paper-defined
    # study sub-populations (Table 1, Final covariate model, All data,
    # Intra-subject variance group rows *1..*5). Switched per subject
    # in model() using AGE thresholds + DIS_HEALTHY + DIS_CANCER_PED.
    # All five names are declared in `paper_specific_residual_sds`
    # because the canonical residual-error matcher recognises only
    # bare `propSd` / `addSd` / `expSd` or per-output `propSd_<obs>`
    # suffixes -- not per-cohort suffixes.
    propSd_young_vol   <- 0.125; label("Proportional residual SD, young healthy volunteers (fraction)")   # Table 1, Final covariate model, All data, group *1 (young volunteers) CV 12.5%
    propSd_elderly_vol <- 0.133; label("Proportional residual SD, elderly healthy volunteers (fraction)") # Table 1, Final covariate model, All data, group *2 (elderly volunteers) CV 13.3%
    propSd_aged_vol    <- 0.169; label("Proportional residual SD, aged healthy volunteers (fraction)")    # Table 1, Final covariate model, All data, group *3 (aged volunteers) CV 16.9%
    propSd_paed_chemo  <- 0.178; label("Proportional residual SD, paediatric cancer chemotherapy (fraction)") # Table 1, Final covariate model, All data, group *4 (paediatric patients on chemotherapy) CV 17.8%
    propSd_paed_anaes  <- 0.145; label("Proportional residual SD, paediatric general anaesthesia (fraction)") # Table 1, Final covariate model, All data, group *5 (paediatrics on anaesthesia) CV 14.5%
  })

  model({
    # -----------------------------------------------------------------
    # Typical-value structural parameters from the linear-regression
    # covariate model. Sex is dispatched via SEXF (1 = female, 0 = male)
    # consistent with the paper's coding ("gender 1 for females and 0
    # for males", Table 2 footnote #).

    # CL: separate male and female equations, then dispatched by SEXF.
    cl_typ_male   <- th_cl_int_m + th_cl_wt_m * WT + th_cl_age_m * AGE
    cl_typ_female <- th_cl_int_f + th_cl_wt_f * WT + th_cl_age_f * AGE
    cl_typ        <- (1 - SEXF) * cl_typ_male + SEXF * cl_typ_female

    # V1 / Vss: single equations (no sex stratification).
    v1_typ  <- th_v1_int  + th_v1_wt  * WT
    vss_typ <- th_vss_int + th_vss_wt * WT + th_vss_age * AGE

    # CLd: separate male and female equations, then dispatched by SEXF.
    cld_typ_male   <- th_cld_int_m + th_cld_wt_m * WT
    cld_typ_female <- th_cld_int_f + th_cld_wt_f * WT
    cld_typ        <- (1 - SEXF) * cld_typ_male + SEXF * cld_typ_female

    # Individual PK parameters: multiplicative log-normal IIV on each.
    cl  <- cl_typ  * exp(etalcl)
    vc  <- v1_typ  * exp(etalvc)
    vss <- vss_typ * exp(etalvss)
    q   <- cld_typ * exp(etalq)

    # Vp (peripheral volume) derived from Vss - V1 = Vss - Vc. Floor at
    # 0.01 L so the peripheral-compartment k12 / k21 micro-constants
    # remain finite under simulation draws where exp(etalvc) shifts Vc
    # above the simulated Vss. The training cohort never hits this
    # regime; the floor only matters for outlier IIV draws.
    vp  <- max(vss - vc, 0.01)

    # Micro-constants for the explicit 2-compartment ODE.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: zero-order i.v. infusion delivered directly into the
    # central compartment via the data event table (cmt = "central",
    # amt = total dose, dur = infusion duration). No depot is needed
    # because the paper uses i.v. dosing only.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration in ng/mL. central is in mg (cumulative integral of
    # an infusion that delivers a dose in mg); vc is in L; mg/L = ug/mL;
    # multiply by 1000 to express as ng/mL, matching the paper's HPLC
    # assay units (LLOQ 1 ng/mL, linear calibration 1-250 ng/mL).
    Cc <- central / vc * 1000

    # -----------------------------------------------------------------
    # Proportional residual error: switched between five sub-populations
    # using existing canonical covariates AGE / DIS_HEALTHY /
    # DIS_CANCER_PED. The thresholds 18 y, 45 y, and 75 y separate the
    # paper's age bands cleanly (study 1 young 18-41 y / study 2 young
    # 21-38 y < 45 y; study 1 elderly 65-75 y / study 2 elderly 61-74 y
    # in [45, 75); study 2 aged 75-82 y >= 75 y; both paediatric studies
    # < 18 y). At inference time, exactly one of the five indicator
    # products is 1 for any given subject and the corresponding
    # propSd_<group> magnitude is selected. Edge cases at the threshold
    # boundaries (a healthy volunteer at exactly AGE = 45 or AGE = 75)
    # default to the higher-age stratum via the strict-less-than form
    # used here -- the training cohort has no subjects at the boundaries.
    is_paed         <- (1 - DIS_HEALTHY)
    is_paed_chemo   <- is_paed * DIS_CANCER_PED
    is_paed_anaes   <- is_paed * (1 - DIS_CANCER_PED)
    is_young_vol    <- DIS_HEALTHY * (AGE < 45)
    is_elderly_vol  <- DIS_HEALTHY * (AGE >= 45) * (AGE < 75)
    is_aged_vol     <- DIS_HEALTHY * (AGE >= 75)

    propSd <- propSd_young_vol   * is_young_vol    +
              propSd_elderly_vol * is_elderly_vol  +
              propSd_aged_vol    * is_aged_vol     +
              propSd_paed_chemo  * is_paed_chemo   +
              propSd_paed_anaes  * is_paed_anaes

    Cc ~ prop(propSd)
  })
}
