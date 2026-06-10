PerezRuixo_2006_tipifarnib <- function() {
  description <- "Three-compartment population PK model for oral and IV tipifarnib in healthy subjects and adult cancer patients (Perez-Ruixo 2006). Sequential zero-order release into the depot (duration D1) followed by first-order absorption (Ka) into the central compartment, with absorption lag time, linear elimination, two peripheral compartments, and bioavailability fixed at 26.7 percent. Covariate effects retained in the final model are total bilirubin on CL (power exponent -0.103 centred at 9 umol/L) and body weight on V2 (linear scaling, exponent fixed at 1, centred at 70 kg); healthy-vs-cancer cohort multipliers apply to CL, V2, Q4, V4, and Ka; a solution-vs-solid formulation indicator scales D1, Ka, and tlag. The mixture-model lag-time subpopulation (71.7 percent subpop 1 vs 28.3 percent subpop 2) is collapsed to the typical subpop-1 lag time for library simulation use; correlated IIVs with paper-reported correlation 1 (Q3-V3, CL-Q4, CL-V4) are encoded as derived etas via the published variance-expansion factors."
  reference   <- "Perez-Ruixo JJ, Piotrovskij V, Zhang S, Hayes S, De Porre P, Zannikos P. Population pharmacokinetics of tipifarnib in healthy subjects and adult cancer patients. Br J Clin Pharmacol. 2006;62(1):81-96. doi:10.1111/j.1365-2125.2006.02615.x"
  vignette    <- "PerezRuixo_2006_tipifarnib"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on V2 (central volume of distribution) with reference 70 kg (population denominator used by Perez-Ruixo 2006 Table 4 footnote and Methods covariate equation). The power coefficient on body weight was estimated and not statistically different from 1; the paper therefore fixed the exponent at 1 (Perez-Ruixo 2006 Results, second paragraph after backward elimination, and Table 3 footnote 2). Combined-data-set median body weight 70.0 kg (range 34.0-145).",
      source_name        = "WGT"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on systemic clearance with reference 9 umol/L (combined-data-set median; Perez-Ruixo 2006 Table 4 footnote a). The covariate equation is CL = 21.9 * (TBILI / 9)^-0.103. A 2-fold increase in TBILI corresponds to a 6.9 percent decrease in CL (Discussion). Combined-data-set median 10.0 umol/L (range 2.0-116).",
      source_name        = "TBIL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult cancer patient; reference cohort is the pooled advanced-cancer cohort across phase 1, 2, and 3 studies)",
      notes              = "Multiplicative ratio-form effects on five structural parameters: CL multiplier 1.21, V2 multiplier 0.55, Q4 multiplier 8.83, V4 multiplier 2.66, and Ka multiplier 2.31; Q3, V3, D1, tlag, and F are equal between cancer patients and healthy subjects (Perez-Ruixo 2006 Table 4 'Ratio healthy subjects : cancer subjects' column). The reference complement is the cancer-patient cohort, so DIS_HEALTHY = 0 leaves all ratio multipliers inactive.",
      source_name        = "HEALTHY"
    ),
    FORM_SOLUTION = list(
      description        = "Oral solution formulation indicator (paper's 'liquid' formulation)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (solid oral formulation: capsule or tablet; the two solid forms had statistically indistinguishable absorption parameters in Perez-Ruixo 2006 and are pooled as the reference)",
      notes              = "Per-dose-record indicator. Multiplicative ratio-form effects relative to the solid reference: D1 multiplier 0.348 (faster zero-order release from the solution), Ka multiplier 2.07, and tlag multiplier 0.183 (Perez-Ruixo 2006 Table 4 footnotes i, j, k). Bioavailability (F) is identical between solid and solution formulations (Discussion, 'Tipifarnib oral bioavailability did not differ between formulations'). The IV-administration routes (1 h, 2 h, 24 h infusions) are encoded by directing the dose to the central compartment with the paper's reported infusion duration; FORM_SOLUTION applies only to oral dose records.",
      source_name        = "FORM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1083L,
    n_studies      = 15L,
    age_range      = "18-89 years (combined-data-set range)",
    age_median     = "60 years",
    weight_range   = "34.0-145 kg (combined-data-set range)",
    weight_median  = "70.0 kg",
    sex_female_pct = 44.3,
    race_ethnicity = c(Caucasian = 93.6, African_American = 1.8, Other = 4.5),
    disease_state  = "Pooled cohort: 1035 adult cancer patients with advanced solid tumors or haematological malignancies (advanced breast, small-cell lung, urothelial transitional-cell, superficial bladder, advanced colorectal, advanced pancreatic, and acute myeloid leukaemia) and 48 healthy subjects.",
    dose_range     = "Oral 25-1300 mg single dose or twice-daily; IV 1, 2, or 24 h infusions of 50-500 mg (combined data set: 7339 plasma concentrations).",
    regions        = "Multinational (15 pooled phase 1, 2, and 3 studies).",
    n_observations = 7339L,
    formulations   = "Three oral formulations (solution, capsule, tablet) plus 1-h, 2-h, and 24-h IV infusions. The solid forms (capsule and tablet) shared a common absorption profile in the final model; the solution differed in D1, Ka, and tlag.",
    index_test_split = "Index data set: 7 phase 1 studies, 166 subjects (154 cancer + 12 healthy), 3445 concentrations; test data set: 5 phase 2 + 2 phase 3 + 1 phase 1 study, 917 subjects (881 cancer + 36 healthy), 3894 concentrations. The final model was fit to the combined data set (Table 4).",
    notes          = "Demographic counts and ranges from Perez-Ruixo 2006 Table 2 (combined data set). Sex distribution: 603 male, 480 female. Per Table 2 the missing-covariate percentage for sex is 0; for the race composite the combined data set lists 1014 Caucasian, 20 African-American, and 49 Other (with 0 percent missing). The 7339 plasma concentrations include 3445 from the index development data set and 3894 from the test data set (final-model fit used the combined merged data set after dropping two outliers from each data set, total 7283 retained)."
  )

  ini({
    # Structural parameters at the reference subject:
    #   DIS_HEALTHY = 0 (cancer patient), FORM_SOLUTION = 0 (solid oral form),
    #   WT = 70 kg, TBILI = 9 umol/L. All apparent / true clearances in L/h,
    #   volumes in L, ka in 1/h, durations in h.
    lcl  <- log(21.9);  label("Systemic clearance (CL, L/h) for cancer patients at TBILI = 9 umol/L")         # Perez-Ruixo 2006 Table 4: CL cancer = 21.9 L/h (footnote a normalises to TBIL = 9 umol/L)
    lvc  <- log(54.9);  label("Central volume of distribution (V2, L per 70 kg) for cancer patients")          # Perez-Ruixo 2006 Table 4: V2 cancer = 54.9 L / 70 kg
    lq   <- log(4.11);  label("Inter-compartmental clearance Q3 to peripheral1 (L/h)")                         # Perez-Ruixo 2006 Table 4: Q3 = 4.11 L/h (cancer = healthy)
    lvp  <- log(92.4);  label("First peripheral volume of distribution (V3, L)")                                # Perez-Ruixo 2006 Table 4: V3 = 92.4 L (cancer = healthy)
    lq2  <- log(14.8);  label("Inter-compartmental clearance Q4 to peripheral2 (L/h) for cancer patients")     # Perez-Ruixo 2006 Table 4: Q4 cancer = 14.8 L/h
    lvp2 <- log(21.4);  label("Second peripheral volume of distribution (V4, L) for cancer patients")          # Perez-Ruixo 2006 Table 4: V4 cancer = 21.4 L
    ld1  <- log(1.20);  label("Duration of zero-order release into depot D1 (h) for the solid formulation")    # Perez-Ruixo 2006 Table 4: D1 solid = 1.20 h
    lka  <- log(0.71);  label("First-order absorption rate constant Ka (1/h) for cancer patients, solid form") # Perez-Ruixo 2006 Table 4: Ka cancer solid = 0.71 1/h
    ltlag <- log(0.11); label("Absorption lag time (h) for the solid formulation, subpopulation 1")            # Perez-Ruixo 2006 Table 4: tlag solid subpop 1 = 0.11 h (representing 71.7 percent of subjects)

    # Bioavailability anchor. Perez-Ruixo 2006 used a logit transformation
    # constrained to [0, 1]; the typical value back-transforms to F = 0.267.
    # Stored here as logit(F) so the IIV (paper-reported logit-domain SD
    # 0.74) applies on the same scale.
    logitfdepot <- logit(0.267); label("Logit of absolute bioavailability into depot (F = 0.267)")             # Perez-Ruixo 2006 Table 4: Fabs = 26.7 percent

    # Healthy-vs-cancer multipliers (Perez-Ruixo 2006 Table 4 'Ratio healthy
    # subjects : cancer subjects' column). The reference cohort is cancer
    # patients (DIS_HEALTHY = 0); healthy subjects (DIS_HEALTHY = 1) receive
    # the indicated multiplicative ratio on each structural parameter.
    e_healthy_cl  <- 1.21; label("Healthy-vs-cancer CL ratio (multiplicative; applied as ratio^DIS_HEALTHY)")  # Perez-Ruixo 2006 Table 4 ratio CL = 1.21
    e_healthy_vc  <- 0.55; label("Healthy-vs-cancer V2 ratio (multiplicative; applied as ratio^DIS_HEALTHY)")  # Perez-Ruixo 2006 Table 4 ratio V2 = 0.55
    e_healthy_q2  <- 8.83; label("Healthy-vs-cancer Q4 ratio (multiplicative; applied as ratio^DIS_HEALTHY)")  # Perez-Ruixo 2006 Table 4 ratio Q4 = 8.83
    e_healthy_vp2 <- 2.66; label("Healthy-vs-cancer V4 ratio (multiplicative; applied as ratio^DIS_HEALTHY)")  # Perez-Ruixo 2006 Table 4 ratio V4 = 2.66
    e_healthy_ka  <- 2.31; label("Healthy-vs-cancer Ka ratio (multiplicative; applied as ratio^DIS_HEALTHY)")  # Perez-Ruixo 2006 Table 4 ratio Ka = 2.31

    # Solution-vs-solid formulation multipliers (Perez-Ruixo 2006 Table 4
    # footnotes i, j, k). The reference (FORM_SOLUTION = 0) is the pooled
    # solid oral cohort (capsule or tablet; both had statistically
    # indistinguishable absorption profiles per Results).
    e_solution_d1   <- 0.348; label("Solution-vs-solid D1 ratio (multiplicative; applied as ratio^FORM_SOLUTION)") # Perez-Ruixo 2006 Table 4 footnote i
    e_solution_ka   <- 2.07;  label("Solution-vs-solid Ka ratio (multiplicative; applied as ratio^FORM_SOLUTION)") # Perez-Ruixo 2006 Table 4 footnote j
    e_solution_tlag <- 0.183; label("Solution-vs-solid tlag ratio (multiplicative; applied as ratio^FORM_SOLUTION)") # Perez-Ruixo 2006 Table 4 footnote k

    # Continuous-covariate exponents.
    allo_vc      <- fixed(1);     label("WT exponent on V2 (unitless; fixed at 1 per Table 3 footnote 2)")    # Perez-Ruixo 2006 Table 3 footnote 2: exponent not different from 1, set to 1
    e_tbili_cl   <- -0.103;       label("TBILI power exponent on CL (unitless; reference 9 umol/L)")          # Perez-Ruixo 2006 Table 4 footnote a: theta_TBIL = -0.103

    # Inter-individual variability (cancer-patient values; Perez-Ruixo 2006
    # Table 4 IIV column). Lognormal variances computed as
    # omega^2 = log(1 + CV^2):
    #   CL  24.9 percent -> log(1 + 0.249^2) = 0.06037
    #   V2  20.3 percent -> log(1 + 0.203^2) = 0.04031
    #   Q3  74.0 percent -> log(1 + 0.74^2)  = 0.43669
    #   D1  52.7 percent -> log(1 + 0.527^2) = 0.24527
    #   Ka  86.1 percent -> log(1 + 0.861^2) = 0.55478
    # The F (Fabs) IIV is reported as a logit-domain standard deviation of
    # 0.74 (Table 4 footnote g), so the stored variance is 0.74^2 = 0.5476
    # on the logit scale rather than on the log scale.
    etalcl         ~ 0.06037   # Perez-Ruixo 2006 Table 4 IIV CL 24.9 percent CV
    etalvc         ~ 0.04031   # Perez-Ruixo 2006 Table 4 IIV V2 20.3 percent CV
    etalq          ~ 0.43669   # Perez-Ruixo 2006 Table 4 IIV Q3 74.0 percent CV
    etald1         ~ 0.24527   # Perez-Ruixo 2006 Table 4 IIV D1 52.7 percent CV
    etalka         ~ 0.55478   # Perez-Ruixo 2006 Table 4 IIV Ka 86.1 percent CV
    etalogitfdepot ~ 0.54760   # Perez-Ruixo 2006 Table 4 IIV Fabs SD = 0.74 (logit domain; variance = 0.74^2)

    # Residual error. Perez-Ruixo 2006 reports a stratified residual model:
    # 24.5 percent CV for full PK profiles (Table 4 footnote h), with larger
    # magnitudes for isolated samples in phase 1 (43.8 percent CV) and
    # phase 2/3 (72.3 percent CV) studies. The library uses the full-profile
    # value as the representative analytical-method residual error; the
    # vignette Assumptions and deviations section documents the alternative
    # arms.
    propSd <- 0.245; label("Proportional residual error (fraction; full PK profile)")                          # Perez-Ruixo 2006 Table 4 footnote h: full PK profiles = 24.5 percent CV
  })

  model({
    # Reference covariate values for the typical-value structural parameters.
    ref_wt   <- 70    # Perez-Ruixo 2006 Table 4: V2 normalised to 70 kg
    ref_tbil <- 9     # Perez-Ruixo 2006 Table 4 footnote a: CL normalised to TBILI = 9 umol/L

    # Healthy-vs-cancer ratio-form effects on five structural parameters.
    # DIS_HEALTHY = 0 (cancer reference) leaves every ratio inactive; = 1
    # applies the per-parameter healthy:cancer multiplier.
    f_healthy_cl  <- e_healthy_cl  ^ DIS_HEALTHY
    f_healthy_vc  <- e_healthy_vc  ^ DIS_HEALTHY
    f_healthy_q2  <- e_healthy_q2  ^ DIS_HEALTHY
    f_healthy_vp2 <- e_healthy_vp2 ^ DIS_HEALTHY
    f_healthy_ka  <- e_healthy_ka  ^ DIS_HEALTHY

    # Solution-vs-solid formulation effects on absorption.
    f_solution_d1   <- e_solution_d1   ^ FORM_SOLUTION
    f_solution_ka   <- e_solution_ka   ^ FORM_SOLUTION
    f_solution_tlag <- e_solution_tlag ^ FORM_SOLUTION

    # Correlated-IIV expansion factors (Perez-Ruixo 2006 Table 4 footnotes
    # c, d, e: 'Correlation between IIV of <X> and <Y> set to 1. Expansion
    # factor of <Y> is <f>'). The expansion factor is reported on the
    # variance scale, so the perfectly correlated derived eta is the square
    # root of the factor times the primary eta.
    eta_lvp_corr  <- sqrt(1.21) * etalq    # V3 IIV perfectly correlated to Q3, var(V3) = 1.21 * var(Q3)
    eta_lq2_corr  <- sqrt(2.08) * etalcl   # Q4 IIV perfectly correlated to CL, var(Q4) = 2.08 * var(CL)
    eta_lvp2_corr <- sqrt(0.95) * etalcl   # V4 IIV perfectly correlated to CL, var(V4) = 0.95 * var(CL)

    # Individual structural PK parameters with the Perez-Ruixo 2006
    # covariate equations applied. Bilirubin enters as a centred power law
    # on CL; body weight enters as a linear (exponent 1) ratio on V2.
    cl  <- exp(lcl + etalcl) * (TBILI / ref_tbil) ^ e_tbili_cl * f_healthy_cl
    vc  <- exp(lvc + etalvc) * (WT / ref_wt) ^ allo_vc          * f_healthy_vc
    q   <- exp(lq  + etalq)
    vp  <- exp(lvp + eta_lvp_corr)
    q2  <- exp(lq2 + eta_lq2_corr)                              * f_healthy_q2
    vp2 <- exp(lvp2 + eta_lvp2_corr)                            * f_healthy_vp2
    d1  <- exp(ld1 + etald1)                                    * f_solution_d1
    ka  <- exp(lka + etalka)                                    * f_healthy_ka  * f_solution_ka
    tlag <- exp(ltlag)                                          * f_solution_tlag
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # Micro-constants for the explicit three-compartment ODE system.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ODE system: depot -> central -> peripheral1, peripheral2.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Sequential zero-order release into depot followed by first-order
    # absorption: the dose enters `depot` over a zero-order window of
    # duration D1, then is absorbed first-order at rate Ka into `central`.
    dur(depot)  <- d1
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Concentration: dose in mg, volumes in L -> central / vc in mg/L.
    # Multiply by 1000 to express Cc in ng/mL to match Perez-Ruixo 2006
    # Methods (LC-MS/MS LLOQ 0.5 ng/mL; HPLC-UV LLOQ 1.0 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
