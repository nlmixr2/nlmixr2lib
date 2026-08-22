Mo_2025_everolimus <- function() {
  description <- "Two-compartment population PK model for a newly developed dispersible everolimus tablet (SVG-101) in 26 healthy adult Korean males from a randomized, open-label, single-dose, 2x2 crossover bioequivalence study against the reference Afinitor 5 mg tablet (Mo 2025). Only the test (dispersible) arm was used for model development. Absorption is parallel dual-input: fraction Fr_zo = 0.62 of the dose enters the central compartment directly by zero-order input over duration Tk0 = 1.12 h, and the remaining (1 - Fr_zo) = 0.38 enters the intestinal (depot) compartment and is absorbed first-order (ka = 9.48 1/h) after a lag time Tlag = 0.63 h. IIV on all eight structural parameters (log-normal, except the bounded zero-order fraction, which is logit-normal); proportional residual error only. No covariates were retained: age, weight, height, and BMI were screened and none was significant."
  reference <- paste(
    "Mo KH, Chae D, Jung YS, Jin BH, Keum DH, Choi MK, Cha JS, Park MS,",
    "Kim CO (2025). Application of population pharmacokinetic modeling of",
    "SVG-101 to evaluate proper dose selection.",
    "Transl Clin Pharmacol 33(3):156-167. doi:10.12793/tcp.2025.33.e14.",
    sep = " "
  )
  vignette <- "Mo_2025_everolimus"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Everolimus is assayed in whole blood (the everolimus
  # therapeutic-range convention of 5-15 ng/mL that Mo 2025 targets in the
  # dose simulation is a whole-blood trough range).
  compartmentData <- list(
    # `depot` is the intestinal compartment G of Mo 2025 Fig. 1.
    depot       = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "everolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "everolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  # Mo 2025 Methods 'Population pharmacokinetic modeling': "we performed a
  # comprehensive covariate exploration, which included examining a variety
  # of continuous covariates, such as age, weight, height, and BMI", and
  # Results 'Model development': "No significant covariates were identified
  # in this analysis." These are therefore screened-but-not-retained
  # covariates: documented for provenance, absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at screening",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the continuous-covariate exploration but not retained in the final model (Mo 2025 Methods 'Population pharmacokinetic modeling'; Results 'Model development'). Cohort value 31.0 +/- 10.8 years (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained (Mo 2025 Results 'Model development'). Cohort value 67.3 +/- 8.7 kg (Table 1); inclusion required at least 55 kg. Mo 2025 Discussion notes that weight HAS been identified as a significant covariate in pediatric everolimus populations (their reference 19) and attributes its absence here to the narrow healthy-adult-male weight range."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Mo 2025 Results 'Model development'). Cohort value 171.8 +/- 6.6 cm (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Mo 2025 Results 'Model development'). Cohort value 22.7 +/- 2.3 kg/m^2 (Table 1); inclusion required 18.5-27.0 kg/m^2."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    age_range      = "19-55 years by protocol inclusion criterion; cohort 31.0 +/- 10.8 years (mean +/- SD, Mo 2025 Table 1)",
    age_median     = "31.0 years (mean)",
    weight_range   = "at least 55 kg by protocol inclusion criterion; cohort 67.3 +/- 8.7 kg (mean +/- SD, Mo 2025 Table 1)",
    weight_median  = "67.3 kg (mean)",
    height_range   = "171.8 +/- 6.6 cm (mean +/- SD, Mo 2025 Table 1); full range not tabulated",
    height_median  = "171.8 cm (mean)",
    bmi_range      = "18.5-27.0 kg/m^2 by protocol inclusion criterion; cohort 22.7 +/- 2.3 kg/m^2 (mean +/- SD, Mo 2025 Table 1)",
    sex_female_pct = 0,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adult Korean male volunteers; participants with clinically significant hepatic, renal, neurological, immunological, respiratory, psychiatric, or gastrointestinal disorders were excluded, as were those with a history of hypersensitivity to everolimus or other rapamycin-related compounds. Concomitant medication was prohibited unless deemed necessary by the principal investigator.",
    dose_range     = "Single oral 5 mg dose in each of two periods of a 2x2 crossover: everolimus 5 mg dispersible tablet (SVG-101, test) and everolimus 5 mg tablet (Afinitor, Novartis Korea, reference), separated by a washout of at least 10 days. Only the 5 mg test-formulation occasion contributed to the population PK model.",
    regions        = "Republic of Korea (Severance Hospital, Yonsei University Health System, Seoul)",
    notes          = "26 participants enrolled and randomized 1:1 to the R-T (n = 13) and T-R (n = 13) sequences; 842 concentration samples from the test formulation were used for the population PK analysis (all 26 subjects). One subject withdrew after completing the PK schedule through day 4 but their test-formulation data were retained for model development, so the NCA in Tables 2 and 3 is based on n = 25 while the model is based on n = 26. Sampling at pre-dose and 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 24, 32, 48, 72, 96, 120, and 144 h post-dose. Estimation in Monolix 2023R1 (SAEM), with below-the-limit-of-quantification records handled by left censoring (censoring/limit columns) rather than exclusion. Dose linearity over the studied range was assumed, so CL/F and Vd/F are treated as dose-independent."
  )

  ini({
    # ====================================================================
    # Final population PK parameter estimates -- Mo 2025 Table 4
    # 'Population pharmacokinetic parameter estimates of the model
    # (n = 26)'. Estimates are printed as value (standard error).
    #
    # Structural model (Mo 2025 Fig. 1 schematic + the three displayed
    # ODEs on p. 161): a 2-compartment disposition with linear
    # elimination from the central compartment, fed by two PARALLEL
    # absorption arms that split the administered dose:
    #   * a zero-order arm carrying fraction Fr_zo of the dose directly
    #     into the central compartment over duration Tk0, and
    #   * a first-order arm carrying the remaining (1 - Fr_zo) through
    #     the intestinal compartment G at rate ka, starting after a lag
    #     time Tlag.
    # The paper's displayed ODEs contain only the first-order arm
    # (dG/dt = -ka*G, with ka*G entering dAc/dt) because the zero-order
    # arm is an administration property in Monolix, not an ODE term --
    # it is encoded below via f(central) / dur(central).
    #
    # All parameters are apparent (CL/F, Vc/F, ...): the study is oral
    # only, so absolute bioavailability is not identifiable and is
    # absorbed into the volume and clearance terms. The paper prints
    # them without the /F suffix; the labels below restore it.
    # ====================================================================
    lka <- log(9.48)   ; label("First-order absorption rate constant ka (1/h)")                          # Mo 2025 Table 4: ka = 9.48 (SE 2.99), %RSE 31.6
    lcl <- log(7.11)   ; label("Apparent clearance CL/F (L/h)")                                          # Mo 2025 Table 4: CL = 7.11 (SE 0.40), %RSE 5.6
    lvc <- log(47.66)  ; label("Apparent central volume of distribution Vc/F (L)")                        # Mo 2025 Table 4: Vc = 47.66 (SE 2.61), %RSE 5.5
    lq  <- log(14.32)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                        # Mo 2025 Table 4: Q = 14.32 (SE 0.90), %RSE 6.3
    lvp <- log(205.94) ; label("Apparent peripheral volume of distribution Vp/F (L)")                     # Mo 2025 Table 4: Vp = 205.94 (SE 10.89), %RSE 5.3

    ltk0  <- log(1.12) ; label("Duration of the zero-order absorption arm Tk0 (h)")                       # Mo 2025 Table 4: Tk0 = 1.12 (SE 0.14), %RSE 12.9 -- Table 4 labels Tk0 a "rate constant", but the Fig. 1 caption defines it as "zero-order absorption duration (h)"; the duration reading is used (see the vignette Errata)
    ltlag <- log(0.63) ; label("Lag time of the first-order absorption arm Tlag (h)")                     # Mo 2025 Table 4: Tlag = 0.63 (SE 0.15), %RSE 24.2

    # Fraction of the dose taken up by the zero-order arm. Mo 2025 Table 4
    # describes F as the "Proportion of the dose absorbed via the zero-order
    # process" -- it is a bounded proportion, NOT a bioavailability, so it is
    # carried on the logit scale (canonical fr_zo, founding example
    # Horita_2018_rifampicin.R) so that simulated values respect 0 < F < 1.
    logitfr_zo <- qlogis(0.62) ; label("Logit of the fraction of dose absorbed via the zero-order arm (unitless)")  # Mo 2025 Table 4: F = 0.62 (SE 0.05), %RSE 8.0

    # ====================================================================
    # Inter-individual variability -- Mo 2025 Table 4 'Random effects'.
    # The table prints omega as the standard deviation of the random
    # effect on its own transformed scale, alongside a "C.V. (%)" column.
    # nlmixr2's ini() takes the VARIANCE, so each entry below is omega^2.
    #
    # The omega scale is proven row by row against the printed C.V. (%)
    # column via the log-normal identity CV = sqrt(exp(omega^2) - 1):
    #   ka   0.94 -> 119.15% vs 118.88 printed
    #   Tk0  0.48 ->  50.90% vs  50.82 printed
    #   Tlag 0.78 ->  91.51% vs  91.25 printed
    #   Vc   0.26 ->  26.45% vs  26.66 printed
    #   Vp   0.26 ->  26.45% vs  26.38 printed
    #   Q    0.31 ->  31.76% vs  31.30 printed
    #   CL   0.28 ->  28.56% vs  28.76 printed
    # Every row agrees to within rounding of the two-decimal omega.
    #
    # F is the sole exception and it falsifies the paper's blanket prose
    # claim ("We assumed a lognormal distribution for all the model
    # parameters"): omega_F = 0.76 under a log-normal would give
    # CV = 88.4%, but Table 4 prints 27.05%. A LOGIT-normal F with
    # omega = 0.76 about logit(0.62) reproduces CV = 26.8% (4e6-draw
    # Monte Carlo), and a probit-normal gives 39.9%. The printed C.V.
    # column therefore pins F to the logit-normal scale -- which is also
    # Monolix's default distribution for a bounded proportion. See the
    # vignette Errata.
    # ====================================================================
    etalka        ~ 0.8836  # Mo 2025 Table 4: omega_ka   = 0.94 (SE 0.32), %RSE 34.5, C.V. 118.88% -- 0.94^2 = 0.8836
    etalcl        ~ 0.0784  # Mo 2025 Table 4: omega_CL   = 0.28 (SE 0.04), %RSE 14.2, C.V.  28.76% -- 0.28^2 = 0.0784
    etalvc        ~ 0.0676  # Mo 2025 Table 4: omega_Vc   = 0.26 (SE 0.04), %RSE 15.4, C.V.  26.66% -- 0.26^2 = 0.0676
    etalq         ~ 0.0961  # Mo 2025 Table 4: omega_Q    = 0.31 (SE 0.05), %RSE 15.3, C.V.  31.30% -- 0.31^2 = 0.0961
    etalvp        ~ 0.0676  # Mo 2025 Table 4: omega_Vp   = 0.26 (SE 0.04), %RSE 15.0, C.V.  26.38% -- 0.26^2 = 0.0676
    etaltk0       ~ 0.2304  # Mo 2025 Table 4: omega_Tk0  = 0.48 (SE 0.10), %RSE 20.9, C.V.  50.82% -- 0.48^2 = 0.2304
    etaltlag      ~ 0.6084  # Mo 2025 Table 4: omega_Tlag = 0.78 (SE 0.15), %RSE 18.8, C.V.  91.25% -- 0.78^2 = 0.6084
    etalogitfr_zo ~ 0.5776  # Mo 2025 Table 4: omega_F    = 0.76 (SE 0.19), %RSE 24.5, C.V.  27.05% -- 0.76^2 = 0.5776 (logit scale; see the note above)

    # Residual error -- proportional only. Mo 2025 Table 4 'Error model
    # parameter' lists a single row, sigma_prop; no additive or combined
    # component is reported, consistent with the Methods statement that
    # proportional, additive, and combined error models were compared and
    # the best-fitting one retained.
    propSd <- 0.1 ; label("Proportional residual error on Cc (fraction)")  # Mo 2025 Table 4: sigma_prop = 0.1 (SE 0.0056), %RSE 5.3
  })

  model({
    # Individual PK parameters. Mo 2025 Methods 'Population pharmacokinetic
    # modeling' places a log-normal random effect on the structural
    # parameters, i.e. P_i = P_TV * exp(eta_i); F is instead logit-normal
    # (see the ini() note), i.e. F_i = expit(logit(F_TV) + eta_i).
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq  + etalq)
    vp  <- exp(lvp + etalvp)

    tk0   <- exp(ltk0  + etaltk0)
    tlag  <- exp(ltlag + etaltlag)
    fr_zo <- expit(logitfr_zo + etalogitfr_zo)

    # Micro-constants for the 2-compartment disposition, matching the
    # displayed ODEs of Mo 2025 p. 161 term for term.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Mo 2025 p. 161, the three displayed differential equations:
    #   dG/dt   = -ka * G
    #   dAc/dt  =  ka * G + (Q/Vp) * Ap - (Q/Vc) * Ac - (CL/Vc) * Ac
    #   dAp/dt  =  (Q/Vc) * Ac - (Q/Vp) * Ap
    # G (the intestinal compartment of Fig. 1) maps to `depot`, Ac to
    # `central`, and Ap to `peripheral1`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Parallel dual absorption (Mo 2025 Fig. 1: "Some of the drug enters
    # the central compartment directly through zero-order absorption,
    # while the rest enters the gut compartment and then transitions into
    # the central compartment through first-order absorption"). The USER
    # encodes each oral administration as TWO parallel dose records in the
    # event table at the same time -- one with cmt = "central" (the
    # zero-order arm) and one with cmt = "depot" (the first-order arm),
    # each carrying the full nominal dose amount. The f() multipliers
    # below then split that nominal amount into the Fr_zo / (1 - Fr_zo)
    # fractions, so the two records together deliver exactly one dose.
    #
    # dur(central) turns the central-compartment record into a zero-order
    # input of duration Tk0 starting at t = 0. alag(depot) delays the
    # first-order arm by Tlag: the abstract and Discussion both describe
    # the model as "zero-order kinetics, FOLLOWED BY first-order kinetics
    # WITH LAG TIME", which attaches the single reported Tlag to the
    # first-order arm, not the zero-order arm. The vignette simulates both
    # readings and shows only this one reproduces the observed median
    # Tmax of 1.50 h.
    f(central)   <- fr_zo
    dur(central) <- tk0
    f(depot)     <- 1 - fr_zo
    alag(depot)  <- tlag

    # Whole-blood everolimus concentration. Doses are in mg and Vc/F is in
    # L, so central/vc is in mg/L; the factor 1000 converts to the ng/mL
    # units in which Mo 2025 reports every concentration (Tables 2-3, the
    # 5-15 ng/mL target trough range, and Figs. 2-4).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
