Choi_2018_metformin <- function() {
  description <- "Two-compartment population PK model for oral metformin in 36 healthy adult Korean men from a phase I single-dose 2-way crossover bioequivalence study comparing a single-agent metformin tablet against a metformin-containing fixed-dose combination (FDC) tablet (Choi 2018). The absorption process is parallel mixed-input: fraction F1 of the dose is absorbed first-order from the depot compartment (rate Ka), and fraction (1-F1) is absorbed zero-order directly into the central compartment over duration D2 with lag time ALAG2. Formulation enters as a binary covariate (FORM_FDC) with multiplicative power-style effects on Ka (Ka_FDC = 0.83 * Ka_single-agent) and on relative bioavailability F (F_FDC = 0.94 * F_single-agent = 0.94). IIV on CL/F, Vc/F (correlated, rho 0.225), and Ka; proportional residual error only."
  reference <- paste(
    "Choi S, Jeon S, Han S (2018). Population pharmacokinetic analysis of",
    "metformin administered as fixed-dose combination in Korean healthy",
    "adults. Transl Clin Pharmacol 26(1):25-31.",
    "doi:10.12793/tcp.2018.26.1.25.",
    sep = " "
  )
  vignette <- "Choi_2018_metformin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FORM_FDC = list(
      description        = "Fixed-dose-combination tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-agent metformin tablet; the typical-value reference for the structural Ka and the implicit relative bioavailability F = 1.0 in Choi 2018 Table 3)",
      notes              = "Per-subject (per-occasion) binary covariate. 1 = subject received the metformin-containing fixed-dose-combination (FDC) tablet during this occasion; 0 = subject received the single-agent metformin tablet. The paper does not name the specific co-formulant drug; typical Korean metformin FDC products co-formulate metformin with sitagliptin, glimepiride, vildagliptin, or dapagliflozin. The 2-way crossover design assigns each of the 36 subjects to both formulation arms across two periods with a 1-week wash-out, so FORM_FDC is per-occasion (not per-subject). Formulation effects enter the model multiplicatively as power-style coefficients (Choi 2018 Methods Eq. 2: theta_test = theta_ref * X^formulation): ka = exp(lka + etalka) * (e_form_fdc_ka^FORM_FDC) shrinks Ka to 83.0% of its single-agent value when FORM_FDC = 1; f_rel = (e_form_fdc_f^FORM_FDC) shrinks the relative bioavailability to 94.0% of 1 when FORM_FDC = 1. Reference category orientation here (single-agent = 0 = ref) is the opposite of the Wilkins 2008 antitubercular precedent (where FDC = 1 = ref); the canonical column value semantics (1 = FDC, 0 = single-drug tablet) are preserved across both papers.",
      source_name        = "formulation"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "20.0-42.0 years (mean 23.9, SD 5.0; Choi 2018 Table 1 + Results 'Dataset' prose)",
    age_median     = "23.9 years (mean)",
    weight_range   = "70.9 kg mean +/- 7.9 kg SD; full range not tabulated (Choi 2018 Results 'Dataset' prose -- the 'cm' unit on weight printed in the paper text 'Mean height and weight were 176.0 +/- 3.5 cm and 70.9 +/- 7.9 cm' is a typo; weight is in kg as confirmed by the Table 1 column header 'Weight (kg)')",
    weight_median  = "70.9 kg (mean)",
    height_range   = "169.1-183.5 cm (mean 176.0, SD 3.5; Choi 2018 Table 1 + Results 'Dataset' prose)",
    height_median  = "176.0 cm (mean)",
    sex_female_pct = 0,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adult Korean male volunteers recruited into a single-dose 2-way crossover bioequivalence study; no T2DM",
    dose_range     = "Single oral metformin dose (specific mg amount not reported in the paper). Each subject received both a single-agent metformin tablet (reference) and a metformin-containing FDC tablet (test) across two periods of a 2-way crossover with a 1-week wash-out, administered with 150 mL water after 10 h of fasting.",
    regions        = "Republic of Korea (Seoul St. Mary's Hospital, Catholic University of Korea, Seoul)",
    notes          = "IRB approval KC14MDSF0913. Sampling at 0 (predose), 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, and 24 h post-dose. Plasma metformin assayed by LC-MS/MS. Choi 2018 Table 1 contains a row-label typo: the values listed against 'Weight (kg)' (176.0; 169.1-183.5) are the height values, and the values listed against 'Height (cm)' (23.9; 20.0-42.0) are the age values; the correct mean +/- SD pairs are reproduced from the prose narrative. Covariates of age, weight, height, serum creatinine, and creatinine clearance were screened in stepwise selection but none was retained at p < 0.05; only formulation was retained as a final-model covariate (Choi 2018 Results 'Covariate Analysis and Formulation Difference')."
  )

  ini({
    # ====================================================================
    # Final population PK parameter estimates -- Choi 2018 Table 3
    # 'Summary of final population PK parameter estimates' (page 28-29).
    # The structural model is a 2-compartment disposition with first-
    # order elimination and parallel mixed-input absorption: a first-
    # order arm via depot (rate Ka) carrying fraction F1 of the dose,
    # and a zero-order arm directly into central (duration D2, lag
    # ALAG2) carrying fraction (1-F1). Reference values are for the
    # single-agent formulation (FORM_FDC = 0).
    # ====================================================================
    lcl <- log(76.7)  ; label("Apparent oral clearance CL/F (L/h)")                                      # Choi 2018 Table 3: CL/F = 76.7 L/h, RSE 3.62%
    lvc <- log(180)   ; label("Apparent central volume Vc/F (L)")                                        # Choi 2018 Table 3: Vc/F = 180 L, RSE 10.2%
    lq  <- log(21.3)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                        # Choi 2018 Table 3: Q/F = 21.3 L/h, RSE 14.1%
    lvp <- log(109)   ; label("Apparent peripheral volume Vp/F (L)")                                     # Choi 2018 Table 3: Vp/F = 109 L, RSE 6.50%
    lka <- log(1.19)  ; label("First-order absorption rate constant Ka (1/h), single-agent reference")   # Choi 2018 Table 3: Ka1 = 1.19 1/h, RSE 15.3%

    ld2    <- log(4.49)  ; label("Duration of zero-order absorption D2 (h)")                              # Choi 2018 Table 3: D2 = 4.49 h, RSE 3.30%
    lalag2 <- log(0.250) ; label("Lag time of zero-order absorption ALAG2 (h)")                           # Choi 2018 Table 3: ALAG2 = 0.250 h, RSE 0.160%
    logitf1 <- qlogis(0.289); label("Logit of fraction F1 absorbed by the first-order arm (unitless)")    # Choi 2018 Table 3: F1 = 0.289, RSE 20.5% -- bounded to [0,1] via logit so the simulated fraction respects 0 <= F1 <= 1

    # Formulation covariate effects. Choi 2018 Methods Eq. 2 parameterizes
    # the covariate as theta_test = theta_ref * X^formulation, so the
    # multiplicative shift X is applied to FORM_FDC = 1 (FDC arm) and the
    # single-agent reference (FORM_FDC = 0) sees no perturbation
    # (X^0 = 1).
    e_form_fdc_ka <- 0.830 ; label("Multiplicative formulation effect on Ka: ka_FDC = ka_ref * e_form_fdc_ka")             # Choi 2018 Table 3: Influence of formulation on Ka = 0.830, RSE 13.1%
    e_form_fdc_f  <- 0.940 ; label("Multiplicative formulation effect on relative bioavailability F: F_FDC = e_form_fdc_f") # Choi 2018 Table 3: Influence of formulation on F = 0.940, RSE 2.02%

    # ====================================================================
    # IIV (Choi 2018 Table 3). CV% on the log-normal scale; the internal
    # variance is omega^2 = log(1 + CV^2). etalcl and etalvc are
    # correlated (rho = 0.225 in Table 3); the off-diagonal of the
    # 2x2 block is cov = rho * sqrt(omega^2_lcl * omega^2_lvc).
    #   omega^2_lcl = log(1 + 0.198^2) ~= 0.03845
    #   omega^2_lvc = log(1 + 0.328^2) ~= 0.10215
    #   cov_lcl_lvc = 0.225 * sqrt(0.03845 * 0.10215) ~= 0.01411
    #   omega^2_lka = log(1 + 0.636^2) ~= 0.33967
    # No IIV is reported on F1, D2, ALAG2, or the formulation
    # coefficients (Choi 2018 Table 3 'Interindividual variability'
    # block lists only CL/F, Vc/F, and Ka).
    # ====================================================================
    etalcl + etalvc ~ c(0.03845, 0.01411, 0.10215)  # Choi 2018 Table 3: CV%(CL/F) = 19.8, CV%(Vc/F) = 32.8, rho = 0.225
    etalka          ~ 0.33967                       # Choi 2018 Table 3: CV%(Ka) = 63.6

    # Residual error -- proportional only (Choi 2018 Table 3 'Residual
    # error' block lists only sigma_prop; no additive component is
    # reported).
    propSd <- 0.259 ; label("Proportional residual error on Cc (fraction)")                              # Choi 2018 Table 3: sigma_prop = 0.259, RSE 2.73%
  })

  model({
    # Individual PK parameters. Choi 2018 Methods Eq. 1 places the
    # interindividual random effects as additive on the log-scale typical
    # value: P_i = P_TV * exp(eta_i). Translated here as exp(lparam +
    # etalparam). Formulation enters Ka multiplicatively as a power-style
    # covariate per Choi 2018 Methods Eq. 2 (theta_test = theta_ref *
    # X^formulation); the single-agent reference (FORM_FDC = 0) recovers
    # the bare typical Ka via X^0 = 1, and the FDC arm (FORM_FDC = 1)
    # multiplies it by e_form_fdc_ka = 0.83.
    ka <- exp(lka + etalka) * (e_form_fdc_ka ^ FORM_FDC)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    d2    <- exp(ld2)
    alag2 <- exp(lalag2)
    f1    <- expit(logitf1)

    # Formulation-dependent relative bioavailability. Choi 2018 Table 3
    # 'Influence of formulation on bioavailability' is F = 0.940 (RSE
    # 2.02%) applied as theta_test = theta_ref * X^formulation, with the
    # single-agent F implicitly = 1. f_rel = X^FORM_FDC is therefore 1.0
    # at FORM_FDC = 0 and 0.94 at FORM_FDC = 1.
    f_rel <- e_form_fdc_f ^ FORM_FDC

    # Micro-constants for the 2-compartment disposition
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 2-compartment ODE system with parallel mixed absorption. The depot
    # compartment receives the first-order absorption fraction of the
    # dose; the central compartment receives the zero-order absorption
    # fraction directly (with a duration-D2 infusion that starts after
    # lag ALAG2). Each oral administration is encoded by the USER as
    # two parallel dose records in the event table: cmt = depot (first-
    # order arm) and cmt = central (zero-order arm). The f() multipliers
    # split the dose into the F1 / (1 - F1) fractions and scale both by
    # f_rel; dur(central) imposes the zero-order infusion duration and
    # lag(central) imposes the zero-order lag time.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)     <- f_rel * f1
    f(central)   <- f_rel * (1 - f1)
    dur(central) <- d2
    lag(central) <- alag2

    # Plasma metformin concentration. Dose in mg, Vc/F in L gives
    # Cc = central / vc in mg/L (equivalent to ug/mL), matching the
    # typical metformin assay calibration range (0.05-3 mg/L plasma
    # after 500-1000 mg oral dosing).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
