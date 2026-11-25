Fidler_2018_olopatadine <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) differential-odds model of the conjunctival allergen challenge (CAC) ocular itching score (0-4 in 0.5-unit steps) after a single bilateral drop of olopatadine 0.1%, 0.2% or 0.7% ophthalmic solution or its vehicle in adults with allergic conjunctivitis (Fidler 2018; pooled CAC trials C-10-126 and C-12-053). The dose strength (%) enters a virtual one-compartment K-PD state that is eliminated first-order (half-life 18.46 h) and drives an Emax drug effect on the cumulative logits of the itching score. A vehicle effect that is linear in time since instillation and scaled by a power of the eye's baseline itching score acts on post-treatment records. Both the vehicle and the drug effect are spread across the score thresholds by differential-odds factors (inverse-logit of a per-threshold parameter, multiplied cumulatively), so they shrink at the higher scores. The eye's baseline itching score also shifts every cumulative logit linearly, with an extra shift in the 0.7% arm, and one between-subject random effect on the logit scale is shared by all thresholds."
  reference <- paste(
    "Fidler ML, Ogundele A, Covert D, Sarangapani R.",
    "Projected 24-hour post-dose ocular itching scores post-treatment with",
    "olopatadine 0.7% versus 0.2%.",
    "J Pharmacokinet Pharmacodyn. 2018;45(4):593-605.",
    "doi:10.1007/s10928-018-9588-7"
  )
  vignette <- "Fidler_2018_olopatadine"
  units <- list(
    time = "h",
    concentration = "%",
    dosing = "% (olopatadine ophthalmic solution strength, w/v; one drop per eye)",
    outcome = "CAC ocular itching score (ordinal, 0-4 in 0.5-unit steps; 9 categories)"
  )

  dosing <- "depot_kpd"

  covariateData <- list(
    OCULAR_ITCH_BL = list(
      description = "Baseline CAC ocular itching score of the eye: the mean of all itching scores recorded for that eye at the screening allergen challenge(s), 3, 5 and 7 min post-challenge",
      units = "(score, 0-4)",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-eye, time-fixed. Used as a surrogate of the patient's allergen sensitivity (Fidler 2018 Methods 'Covariate analysis'). Centred at the median baseline 2.8 (Table 3, overall median; the paper prints the centring as 'Median Covariate' without a number). Enters as a power on the vehicle shift Vb, as a linear shift of every cumulative logit, and (in the 0.7% arm) as an additional linear shift. Analysis-set range 1.5-4.0 (Table 3).",
      source_name = "baseline itching"
    ),
    TRT_PHASE = list(
      description = "Post-treatment record indicator: 1 = the record follows a bilateral instillation of vehicle or olopatadine, 0 = screening (baseline) challenge",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (screening baseline record)",
      notes = "The paper's vehicle indicator Vi ('0 for baseline observations, and 1 with drug or vehicle treatment', Methods). Gates the vehicle effect and the 0.7%-specific baseline shift; the drug effect is gated by the K-PD amount itself, which is zero before the first olopatadine drop.",
      source_name = "Vi"
    ),
    DOSE_HIGH = list(
      description = "Olopatadine 0.7% arm indicator: 1 = the eye was treated with olopatadine 0.7%, 0 = vehicle, olopatadine 0.1% or 0.2%",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (vehicle, 0.1% or 0.2%)",
      notes = "The 0.7% strength is the highest dose in the pooled trials. It gates the additional baseline-itching shift of the score distribution retained in the final model (Table 4 'Olopatadine 0.7% baseline on distribution', Table 5 'psi Olopatadine 0.7 % shift on baseline, distribution'), applied on post-treatment records only.",
      source_name = "olopatadine 0.7% treatment"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on Emax, on the vehicle shift Vb and on the score distribution; every test raised the objective function and none was retained (Fidler 2018 Table 4).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Gender (female indicator I_F) was tested on Emax, Vb and the score distribution and not retained (Fidler 2018 Table 4).",
      source_name = "gender (I_F)"
    )
  )

  compartmentData <- list(
    depot_kpd = list(analyte = "olopatadine", units = "% (dose strength)", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 547L,
    n_studies = 2L,
    age_range = "18-77 years (mean 40.9, SD 13.0)",
    sex_female_pct = 60.2,
    disease_state = "Adults with a history of allergic conjunctivitis and a positive bilateral conjunctival allergen challenge response at screening (itching >= 2 and redness >= 2 in 2 of 3 vessel beds).",
    dose_range = "Single bilateral drop of olopatadine 0.1%, 0.2% or 0.7% or vehicle before each CAC; challenge 27 min (onset), 16 h or 24 h post-instillation.",
    regions = "United States (multicenter).",
    baseline_itching = "Mean 2.8 (SD 0.5), median 2.8 (range 1.5-4.0) (Table 3).",
    allergens = "Grass 42.0%, dust mites 17.4%, ragweed 16.5%, trees 12.8%, cat dander 8.4%, dog dander 1.5%, cockroach 1.3% (Table 3).",
    notes = "Pooled ITT data from trials C-10-126 (NCT01479374) and C-12-053 (NCT01743027): 547 patients with 10,759 itching observations (3269 baseline; vehicle 117 patients, 0.1% 99, 0.2% 167, 0.7% 164). Itching scored 0-4 in 0.5-unit steps by each eye at 3, 5 and 7 min after each challenge. NONMEM 7.3, Laplacian estimation."
  )

  ini({
    # Baseline cumulative log-odds (Table 5, column 'alpha_s'). The paper
    # models logit P(Y >= s), s = 0.5..4; the library convention is
    # logit P(grade <= k) with grade k = 2 * score, so logite0 = -alpha_0.5
    # and logitstep<k> = -alpha_s with k = 2 * s - 1.
    logite0 <- -8.79
    label("Cumulative logit of P(grade <= 0), i.e. no itching, at baseline for the typical eye (unitless)") # Table 5: alpha 0.5 = 8.79 (RSE 2.50%), sign-flipped
    logitstep1 <- 1.47
    label("Increment to the cumulative logit of P(grade <= 1), grade = 2 x itching score (unitless)") # Table 5: alpha 1 = -1.47 (RSE 4.22%), sign-flipped
    logitstep2 <- 1.03
    label("Increment to the cumulative logit of P(grade <= 2), grade = 2 x itching score (unitless)") # Table 5: alpha 1.5 = -1.03 (RSE 4.61%), sign-flipped
    logitstep3 <- 0.85
    label("Increment to the cumulative logit of P(grade <= 3), grade = 2 x itching score (unitless)") # Table 5: alpha 2 = -0.85 (RSE 12.34%), sign-flipped
    logitstep4 <- 1.88
    label("Increment to the cumulative logit of P(grade <= 4), grade = 2 x itching score (unitless)") # Table 5: alpha 2.5 = -1.88 (RSE 5.61%), sign-flipped
    logitstep5 <- 2.88
    label("Increment to the cumulative logit of P(grade <= 5), grade = 2 x itching score (unitless)") # Table 5: alpha 3 = -2.88 (RSE 4.60%), sign-flipped
    logitstep6 <- 3.02
    label("Increment to the cumulative logit of P(grade <= 6), grade = 2 x itching score (unitless)") # Table 5: alpha 3.5 = -3.02 (RSE 3.82%), sign-flipped
    logitstep7 <- 2.14
    label("Increment to the cumulative logit of P(grade <= 7), grade = 2 x itching score (unitless)") # Table 5: alpha 4 = -2.14 (RSE 7.33%), sign-flipped

    # Vehicle effect Ve = Vi * (Vb + Vm * t) on logit P(Y >= s) (Methods
    # equation for Ve(S); Table 5 'Vehicle parameters').
    shift_placebo <- -4.98
    label("Vehicle shift Vb of logit P(itching score >= s) at the instillation time, typical eye (unitless)") # Table 5: Vehicle linear shift (Vb) = -4.98 (RSE 5.65%)
    slope_placebo <- 0.06
    label("Vehicle time slope Vm of logit P(itching score >= s) (1/h)") # Table 5: Vehicle time shift (Vm) = 0.06 (RSE 14.07%)
    e_ocular_itch_bl_shift_placebo <- 2.16
    label("Power exponent of the baseline itching score (OCULAR_ITCH_BL / 2.8) on Vb (unitless)") # Table 5: psi baseline on Vb = 2.16 (RSE 5.04%)

    # Vehicle differential-odds parameters v_s (Table 5 column 'vs'),
    # logit-scale; the factor applied at score s is the product of
    # expit(v_m) for m = 1..s. Index k = 2 * s - 1.
    logitfdiff_placebo1 <- fixed(10)
    label("Logit of the vehicle differential-odds factor at grade 2 (unitless)") # Table 5: v 1 = 10.00 FIXED
    logitfdiff_placebo2 <- fixed(10)
    label("Logit of the vehicle differential-odds factor at grade 3 (unitless)") # Table 5: v 1.5 = 10.00 FIXED
    logitfdiff_placebo3 <- fixed(10)
    label("Logit of the vehicle differential-odds factor at grade 4 (unitless)") # Table 5: v 2 = 10.00 FIXED
    logitfdiff_placebo4 <- 2.38
    label("Logit of the vehicle differential-odds factor at grade 5 (unitless)") # Table 5: v 2.5 = 2.38 (RSE 11.27%)
    logitfdiff_placebo5 <- 1.07
    label("Logit of the vehicle differential-odds factor at grade 6 (unitless)") # Table 5: v 3 = 1.07 (RSE 11.90%)
    logitfdiff_placebo6 <- fixed(10)
    label("Logit of the vehicle differential-odds factor at grade 7 (unitless)") # Table 5: v 3.5 = 10.00 FIXED
    logitfdiff_placebo7 <- fixed(10)
    label("Logit of the vehicle differential-odds factor at grade 8 (unitless)") # Table 5: v 4 = 10.00 FIXED

    # Olopatadine K-PD and Emax (Table 5 'Olopatadine parameters').
    lkel <- log(log(2) / 18.46)
    label("First-order elimination rate of the K-PD virtual compartment, log(2) / t1/2 (1/h)") # Table 5: KPD drug half-life t1/2 = 18.46 h (RSE 98.37%)
    lec50 <- log(0.03)
    label("K-PD amount (in dose-strength units) giving half the maximum drug effect EC50 (%)") # Table 5: EC50 = 0.03 % (RSE 96.36%)
    lemax <- log(3.40)
    label("Maximum olopatadine reduction of logit P(itching score >= s) Emax (unitless)") # Table 5: Maximum effect, Emax = 3.40 (RSE 8.70%)

    # Drug differential-odds parameters delta_s (Table 5 column 'delta s').
    logitfdiff_drug1 <- fixed(10)
    label("Logit of the drug differential-odds factor at grade 2 (unitless)") # Table 5: delta 1 = 10.00 FIXED
    logitfdiff_drug2 <- fixed(10)
    label("Logit of the drug differential-odds factor at grade 3 (unitless)") # Table 5: delta 1.5 = 10.00 FIXED
    logitfdiff_drug3 <- 3.64
    label("Logit of the drug differential-odds factor at grade 4 (unitless)") # Table 5: delta 2 = 3.64 (RSE 48.17%)
    logitfdiff_drug4 <- 1.56
    label("Logit of the drug differential-odds factor at grade 5 (unitless)") # Table 5: delta 2.5 = 1.56 (RSE 17.60%)
    logitfdiff_drug5 <- fixed(10)
    label("Logit of the drug differential-odds factor at grade 6 (unitless)") # Table 5: delta 3 = 10.00 FIXED
    logitfdiff_drug6 <- fixed(-10)
    label("Logit of the drug differential-odds factor at grade 7 (unitless)") # Table 5: delta 3.5 = -10.00 FIXED
    logitfdiff_drug7 <- fixed(10)
    label("Logit of the drug differential-odds factor at grade 8 (unitless)") # Table 5: delta 4 = 10.00 FIXED

    # Covariate shifts of the score distribution (Table 5 'Distribution
    # parameters'; Methods linear-shift equation).
    e_ocular_itch_bl_ocular_itch <- 5.37
    label("Linear shift of logit P(itching score >= s) per unit of OCULAR_ITCH_BL above 2.8 (unitless)") # Table 5: psi baseline, distribution = 5.37 (RSE 3.96%)
    e_dose_high_ocular_itch_bl_ocular_itch <- -0.64
    label("Additional linear shift per unit of OCULAR_ITCH_BL above 2.8 on post-treatment records of the 0.7% arm (unitless)") # Table 5: psi Olopatadine 0.7 % shift on baseline, distribution = -0.64 (RSE 54.17%)

    # Between-subject random effect on the logit scale, shared by every
    # threshold (Table 5: 'Log odds SD (w)' = 1.38 (RSE 4.76%), variance = SD^2).
    etalogite0 ~ 1.9044
  })

  model({
    # Differential-odds factors: the effect at threshold k (score
    # s = (k + 1) / 2) is the product of expit(v_m) over the thresholds up
    # to k; the factor at score 0.5 is 1 (Methods equations De(S), Ve(S)).
    fp1 <- expit(logitfdiff_placebo1)
    fp2 <- fp1 * expit(logitfdiff_placebo2)
    fp3 <- fp2 * expit(logitfdiff_placebo3)
    fp4 <- fp3 * expit(logitfdiff_placebo4)
    fp5 <- fp4 * expit(logitfdiff_placebo5)
    fp6 <- fp5 * expit(logitfdiff_placebo6)
    fp7 <- fp6 * expit(logitfdiff_placebo7)
    fd1 <- expit(logitfdiff_drug1)
    fd2 <- fd1 * expit(logitfdiff_drug2)
    fd3 <- fd2 * expit(logitfdiff_drug3)
    fd4 <- fd3 * expit(logitfdiff_drug4)
    fd5 <- fd4 * expit(logitfdiff_drug5)
    fd6 <- fd5 * expit(logitfdiff_drug6)
    fd7 <- fd6 * expit(logitfdiff_drug7)

    # K-PD virtual compartment holding the instilled dose strength (%).
    kel <- exp(lkel)
    d/dt(depot_kpd) <- -kel * depot_kpd
    eff_drug <- exp(lemax) * depot_kpd / (exp(lec50) + depot_kpd)

    # Vehicle effect on post-treatment records; t is the time since the
    # last instillation (h). Vehicle-arm instillations are amt = 0 doses
    # into depot_kpd so that tad() is defined.
    tveh <- 0
    if (TRT_PHASE == 1) tveh <- tad(depot_kpd)
    shift_veh <- shift_placebo * (OCULAR_ITCH_BL / 2.8)^e_ocular_itch_bl_shift_placebo
    eff_veh <- TRT_PHASE * (shift_veh + slope_placebo * tveh)

    # Linear shift of the whole distribution by the eye's baseline itching,
    # steeper in the 0.7% arm after treatment.
    shift_bl <- (e_ocular_itch_bl_ocular_itch + e_dose_high_ocular_itch_bl_ocular_itch * DOSE_HIGH * TRT_PHASE) * (OCULAR_ITCH_BL - 2.8)

    # logit P(grade <= k) = -logit P(score >= (k + 1) / 2). The paper's
    # vehicle and baseline terms enter logit P(score >= s) with a plus sign,
    # so they are subtracted here; the drug effect reduces logit
    # P(score >= s), so it is added here.
    lp <- logite0 + etalogite0 - shift_bl
    pcum0 <- expit(lp - eff_veh + eff_drug)
    pcum1 <- expit(lp + logitstep1 - eff_veh * fp1 + eff_drug * fd1)
    pcum2 <- expit(lp + logitstep1 + logitstep2 - eff_veh * fp2 + eff_drug * fd2)
    pcum3 <- expit(lp + logitstep1 + logitstep2 + logitstep3 - eff_veh * fp3 + eff_drug * fd3)
    pcum4 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4 - eff_veh * fp4 + eff_drug * fd4)
    pcum5 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4 + logitstep5 - eff_veh * fp5 + eff_drug * fd5)
    pcum6 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4 + logitstep5 + logitstep6 - eff_veh * fp6 + eff_drug * fd6)
    pcum7 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4 + logitstep5 + logitstep6 + logitstep7 - eff_veh * fp7 + eff_drug * fd7)

    # Per-grade probabilities; grade k is itching score k / 2.
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- pcum2 - pcum1
    p3 <- pcum3 - pcum2
    p4 <- pcum4 - pcum3
    p5 <- pcum5 - pcum4
    p6 <- pcum6 - pcum5
    p7 <- pcum7 - pcum6
    # Expected itching score (sum over s of P(score >= s) times 0.5).
    itch_mean <- 0.5 * (8 - pcum0 - pcum1 - pcum2 - pcum3 - pcum4 - pcum5 - pcum6 - pcum7)

    # rxode2 ordinal endpoint: simulated category k is itching score (k - 1) / 2.
    ocular_itch ~ c(p0, p1, p2, p3, p4, p5, p6, p7)
  })
}
