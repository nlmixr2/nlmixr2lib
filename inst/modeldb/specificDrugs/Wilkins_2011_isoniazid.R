Wilkins_2011_isoniazid <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for oral",
    "isoniazid in South African pulmonary tuberculosis patients",
    "(Wilkins 2011; 235 patients, 2352 plasma concentrations).",
    "First-order absorption with an absorption lag time, first-order",
    "elimination, and allometric scaling on all clearance and volume",
    "terms (WT exponent 0.75 on CL and Q, exponent 1 on Vc and Vp,",
    "reference weight 70 kg). A two-class mixture model on apparent",
    "clearance characterises the bimodal isoniazid elimination",
    "phenotype that arises from N-acetyltransferase-2 (NAT2)",
    "polymorphism: typical CL/F is 21.6 L/h in fast eliminators",
    "(13.2 % of subjects) and 9.70 L/h in slow eliminators",
    "(86.8 %). Two covariate effects were retained: female sex",
    "reduces Vc/F by 10.3 % and HIV-positive comorbidity reduces",
    "CL/F by 17.4 %. Inter-individual variability is reported on",
    "CL/F, Vc/F, Q/F, relative bioavailability F, and lag time;",
    "inter-occasion variability on ka (90.1 %) and F (8.4 %) is not",
    "propagated -- see the validation vignette Assumptions and",
    "deviations section for the single-occasion approximation."
  )
  reference <- paste(
    "Wilkins JJ, Langdon G, McIlleron H, Pillai G, Smith PJ, Simonsson USH.",
    "Variability in the population pharmacokinetics of isoniazid in",
    "South African tuberculosis patients.",
    "Br J Clin Pharmacol. 2011;72(1):51-62.",
    "doi:10.1111/j.1365-2125.2011.03940.x."
  )
  vignette <- "Wilkins_2011_isoniazid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = paste(
        "Total body weight at study entry. Time-fixed per subject in the",
        "Wilkins 2011 cohort (a single baseline weight was recorded per",
        "patient and carried across the observation window)."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling uses a 70 kg reference weight per Anderson",
        "& Holford 2008 (Wilkins 2011 Methods, equations 1 and 2):",
        "CL/F and Q/F scale with (WT/70)^0.75; Vc/F and Vp/F scale",
        "with (WT/70)^1. The studied cohort weight range is 33.7-68.0",
        "kg (Table 1 'Combined' column)."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = paste(
        "1 = female, 0 = male. Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Linear multiplicative effect on apparent central volume of",
        "distribution: Vc/F = Vc_typ * (1 + e_sexf_vc * SEXF) with",
        "e_sexf_vc = -0.103 (Wilkins 2011 Table 2 row 'Linear effect",
        "of being female on Vc/F (theta_Vc,gender,F)' = -0.103, RSE",
        "25.7 %). Female patients have Vc/F 10.3 % lower than male",
        "patients of the same body weight. Female proportion in the",
        "studied cohort is 43.4 % (102 / 235; Table 1)."
      ),
      source_name        = "SEX"
    ),
    HIV_POS = list(
      description        = paste(
        "1 = HIV-1 antibody positive at study entry, 0 = HIV-negative.",
        "Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = paste(
        "Linear multiplicative effect on apparent oral clearance:",
        "CL/F = CL_typ * (1 + e_hiv_pos_cl * HIV_POS) with",
        "e_hiv_pos_cl = -0.174 (Wilkins 2011 Table 2 row 'Linear",
        "effect of positive HIV status on CL/F (theta_CL,HIV)' =",
        "-0.174, RSE 22.2 %). HIV-positive patients have CL/F 17.4 %",
        "lower than HIV-negative patients regardless of eliminator",
        "phenotype. HIV prevalence in the studied cohort is 15.2 %",
        "(35 / 230 tested; Table 1). Five patients declined HIV",
        "testing; the paper does not describe how their HIV_POS",
        "value was imputed for the NONMEM fit."
      ),
      source_name        = "HIV"
    ),
    MIX_FAST_ELIM = list(
      description        = paste(
        "Per-subject latent mixture-model class indicator for the",
        "isoniazid eliminator phenotype: 1 = subject classified to the",
        "fast eliminator subpopulation (typical CL/F = 21.6 L/h at 70",
        "kg); 0 = subject classified to the slow eliminator",
        "subpopulation (typical CL/F = 9.70 L/h at 70 kg). Population",
        "probability of MIX_FAST_ELIM = 1 is 0.132 (Wilkins 2011 Table",
        "2 row 'Proportion of fast eliminators in population (P_fast)'",
        "= 0.132, RSE 23.2 %)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (slow eliminator subpopulation, majority class at 86.8 %)",
      notes              = paste(
        "Not a measured clinical covariate -- the mixture assignment",
        "is the per-subject posterior latent-class index from the",
        "NONMEM mixture model on apparent clearance described in",
        "Wilkins 2011 Methods ('The known trimodality of isoniazid",
        "elimination was investigated through the use of a mixture",
        "model for apparent clearance ... could only be parameterized",
        "in terms of two subpopulations rather than three, owing to",
        "insufficient information to differentiate between estimates",
        "of apparent clearance in intermediate acetylators and in fast",
        "acetylators'). For typical-value simulation set",
        "MIX_FAST_ELIM = 0 (dominant slow-eliminator phenotype). For",
        "population simulation, draw MIX_FAST_ELIM ~ Bernoulli(0.132)",
        "per subject. The mixture is gated structurally in model()",
        "via the power-form covariate effect e_mix_fast_elim_cl on",
        "CL/F, so CL/F = exp(lcl + e_mix_fast_elim_cl * MIX_FAST_ELIM",
        "+ etalcl) with e_mix_fast_elim_cl = log(21.6/9.70) recovers",
        "9.70 L/h at MIX_FAST_ELIM = 0 and 21.6 L/h at MIX_FAST_ELIM",
        "= 1 (both at the 70 kg reference weight)."
      ),
      source_name        = "$MIX class assignment"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 235L,
    n_studies       = 2L,
    n_observations  = 2352L,
    age_range       = "20-60 years; median 36 (Table 1 'Combined' column)",
    weight_range    = "33.7-68.0 kg; median 48.0 (Table 1 'Combined' column)",
    sex_female_pct  = 43.4,
    race_ethnicity  = c(
      Coloured  = 81.7,
      Black     = 17.4,
      Caucasian = 0.9
    ),
    disease_state   = paste(
      "Hospitalised pulmonary tuberculosis patients enrolled at two",
      "South African treatment centres -- the DP Marais SANTA Centre",
      "(DPM) near Cape Town (n = 91) and Brewelskloof Hospital (BKH)",
      "in the Breede River Valley (n = 144). Patients were treated",
      "with isoniazid in combination with rifampicin and, as",
      "indicated, pyrazinamide, ethambutol and streptomycin per the",
      "WHO DOTS strategy and South African national guidelines at the",
      "time. Males and non-pregnant females over the age of 18 years.",
      "HIV prevalence in the combined cohort was 15.2 %."
    ),
    dose_range      = paste(
      "Oral isoniazid 100-450 mg daily. DPM: 100 mg (n = 1), 225 mg",
      "(n = 8), 240 mg (n = 29), 300 mg (n = 53), 400 mg (n = 2)",
      "given Monday-Friday for at least 2 weeks. BKH: 200 mg",
      "(n = 1), 300 mg (n = 142), 450 mg (n = 1) given 7 days per",
      "week at the end of the 2-month intensive phase. Per-dose range",
      "spans 3.98-9.59 mg/kg with median 5.88 mg/kg (Table 1)."
    ),
    regions         = "South Africa (Western Cape)",
    notes           = paste(
      "Isoniazid was always co-administered with rifampicin. The",
      "DPM cohort contributed 3 samples per profile twice weekly at",
      "random times 0-12 h post-dose; the BKH cohort contributed a",
      "single-day profile with samples pre-dose and at 0.5, 1.0,",
      "1.5, 2.0, 2.5, 3.0, 4.0, 6.0 and 8.0 h post-dose. Plasma",
      "isoniazid was measured by reversed-phase HPLC-UV (LOQ 0.2",
      "mg/L); 52 of 2352 observations (1.6 %) were below the LOQ",
      "and were retained at their reported values per Wilkins 2011",
      "Methods 'Pharmacokinetic data analysis'. Patients were",
      "fasted overnight (from 22:00 the evening before sampling)."
    )
  )

  ini({
    # ============================================================
    # Structural PK -- two-compartment, first-order absorption with
    # an absorption lag time, first-order elimination. Point
    # estimates from Wilkins 2011 Table 2 ('Population mean' column,
    # final model). All apparent clearance and volume terms are
    # typical values at the 70 kg reference weight; the typical CL/F
    # entered here is the slow-eliminator value, and the fast-vs-slow
    # ratio is encoded as the e_mix_fast_elim_cl power-form effect
    # below so that the MIX_FAST_ELIM = 0 reference matches Table 2.
    # ============================================================
    lcl   <- log(9.70);  label("Apparent oral clearance, slow eliminator at 70 kg (CL_slow/F, L/h)")  # Wilkins 2011 Table 2 row "Typical apparent clearance, slow eliminators (CL_slow/F)" = 9.70 L/h (RSE 3.05 %)
    lvc   <- log(57.7);  label("Apparent central volume of distribution at 70 kg, male (Vc/F, L)")    # Wilkins 2011 Table 2 row "Typical apparent central volume of distribution (Vc/F)" = 57.7 L (RSE 2.81 %)
    lvp   <- log(1730);  label("Apparent peripheral volume of distribution at 70 kg (Vp/F, L)")       # Wilkins 2011 Table 2 row "Typical apparent peripheral volume of distribution (Vp/F)" = 1730 L (RSE 14.7 %)
    lq    <- log(3.34);  label("Apparent inter-compartmental clearance at 70 kg (Q/F, L/h)")          # Wilkins 2011 Table 2 row "Typical apparent intercompartmental clearance (Q/F)" = 3.34 L/h (RSE 8.52 %)
    lka   <- log(1.85);  label("First-order absorption rate constant (ka, 1/h)")                      # Wilkins 2011 Table 2 row "Typical absorption rate constant (ka)" = 1.85 1/h (RSE 4.01 %)
    ltlag <- log(0.180); label("Absorption lag time (tlag, h)")                                       # Wilkins 2011 Table 2 row "Typical absorption lag time (tlag)" = 0.180 h (RSE 7.63 %)

    # Bioavailability anchor. Oral-only data: absolute F is not
    # identifiable; the typical value is fixed at 1.0 and IIV on F
    # is propagated through etalfdepot below.
    lfdepot <- fixed(log(1)); label("Relative bioavailability anchor (F, fixed at 1)")

    # ============================================================
    # Allometric scaling (fixed at the Anderson & Holford 2008
    # canonical values). Wilkins 2011 Methods: "A priori scaling of
    # clearance and volume parameters was tested ... using a
    # reference weight of 70 kg. Typical values of clearance terms
    # were scaled as in the example in Equation (1), and typical
    # volume of distribution terms were scaled as in the example in
    # Equation (2)." Confirmed in Table 2 footnote:
    #   CL/F = CL * (WT/70)^0.75
    #   Vc/F = 57.7 * (WT/70)
    #   Vp/F = 1730 * (WT/70)
    #   Q/F  = 3.34 * (WT/70)^0.75
    # Wilkins 2011 reports the exponents 0.75 and 1 without
    # uncertainty (canonical Anderson-Holford values).
    # ============================================================
    e_wt_cl_q  <- fixed(0.75); label("Shared body-weight allometric exponent on CL/F and Q/F (unitless)")  # Wilkins 2011 Methods equations 1-2 and Table 2 footnote (Anderson & Holford 2008 canonical)
    e_wt_vc_vp <- fixed(1.0);  label("Shared body-weight allometric exponent on Vc/F and Vp/F (unitless)")  # Wilkins 2011 Methods equations 1-2 and Table 2 footnote (Anderson & Holford 2008 canonical)

    # ============================================================
    # Covariate effects on CL/F. Wilkins 2011 Table 2 footnote:
    #   CL/F = CL_typ * (WT/70)^0.75 * (1 - 0.174 * HIV),
    # where CL_typ = 21.6 (fast) or 9.70 (slow) and HIV is the
    # 0/1 indicator. The mixture effect is encoded structurally via
    # the power form below so the reference category MIX_FAST_ELIM
    # = 0 matches the slow-eliminator typical value entered as lcl.
    # ============================================================
    e_mix_fast_elim_cl <- log(21.6 / 9.70); label("Power-form effect of MIX_FAST_ELIM on CL/F (log fast/slow ratio, unitless)")  # Wilkins 2011 Table 2 rows "Typical apparent clearance, fast eliminators" (21.6) and "Typical apparent clearance, slow eliminators" (9.70)
    e_hiv_pos_cl       <- -0.174;            label("Linear fractional effect of HIV_POS on CL/F (unitless)")                       # Wilkins 2011 Table 2 row "Linear effect of positive HIV status on CL/F (theta_CL,HIV)" = -0.174 (RSE 22.2 %)

    # ============================================================
    # Covariate effects on Vc/F. Wilkins 2011 Table 2 footnote:
    #   Vc/F = 57.7 * (WT/70) * (1 - 0.103 * SEX),
    # where SEX = 1 for females and 0 for males (canonical SEXF).
    # ============================================================
    e_sexf_vc <- -0.103; label("Linear fractional effect of SEXF on Vc/F (unitless)")  # Wilkins 2011 Table 2 row "Linear effect of being female on Vc/F (theta_Vc,gender,F)" = -0.103 (RSE 25.7 %)

    # ============================================================
    # Population probability of MIX_FAST_ELIM = 1 (the fast-
    # eliminator subpopulation fraction) is 0.132 (Wilkins 2011
    # Table 2 row "Proportion of fast eliminators in population
    # (P_fast)" = 0.132, RSE 23.2 %). It is not declared as an
    # ini() parameter because the mixture is gated through the
    # per-subject MIX_FAST_ELIM covariate column; downstream users
    # constructing a virtual cohort draw the indicator as
    # MIX_FAST_ELIM ~ Bernoulli(0.132). The mixture-class
    # documentation lives in covariateData[[MIX_FAST_ELIM]]$notes.
    # ============================================================

    # ============================================================
    # Inter-individual variability. Wilkins 2011 Table 2 reports CV%
    # values; convert to log-normal variance via omega^2 = log(1 +
    # CV^2). All IIV terms are diagonal (no correlation block
    # estimated in the final model).
    # ============================================================
    etalcl     ~ log(1 + 0.184^2)
    # Wilkins 2011 Table 2 row "Apparent clearance (omega_CL^2)" = 18.4 % CV; omega^2 = log(1 + 0.184^2) approx 0.0333

    etalvc     ~ log(1 + 0.165^2)
    # Wilkins 2011 Table 2 row "Apparent central volume of distribution (omega_Vc^2)" = 16.5 % CV; omega^2 = log(1 + 0.165^2) approx 0.0269

    etalq      ~ log(1 + 0.931^2)
    # Wilkins 2011 Table 2 row "Apparent intercompartmental clearance (omega_Q^2)" = 93.1 % CV; omega^2 = log(1 + 0.931^2) approx 0.6249

    etalfdepot ~ log(1 + 0.262^2)
    # Wilkins 2011 Table 2 row "Relative bioavailability (omega_F^2)" = 26.2 % CV; omega^2 = log(1 + 0.262^2) approx 0.0665

    etaltlag   ~ log(1 + 0.884^2)
    # Wilkins 2011 Table 2 row "Absorption lag time (omega_tlag^2)" = 88.4 % CV; omega^2 = log(1 + 0.884^2) approx 0.5775

    etalka     ~ log(1 + 0.901^2)
    # Wilkins 2011 Table 2 row "Absorption rate constant (kappa_ka^2)" = 90.1 % CV (interoccasion).
    # The source reports this variability as IOV; a single per-subject random effect on lka is the
    # single-occasion approximation propagated to this model (see the validation vignette
    # Assumptions and deviations section). omega^2 = log(1 + 0.901^2) approx 0.5941.

    # ============================================================
    # Residual unexplained variability. Wilkins 2011 Methods:
    # "Residual variability, arising from unspecified within-subject
    # variability, model misspecification and experimental error,
    # was normally distributed with mean 0 and variance sigma^2; it
    # was applied as an additive model on the log scale." This is
    # the canonical lognormal residual: log(Y_obs) = log(Y_pred) +
    # epsilon with epsilon ~ N(0, sigma^2). The reported value is
    # the standard deviation on the log scale (0.205, RSE 2.74 %)
    # for the DPM cohort; the same parameter is used for the
    # combined dataset in the published final model.
    # ============================================================
    expSd <- 0.205; label("Residual error standard deviation on the log scale (unitless)")  # Wilkins 2011 Table 2 row "Additive variability for DPM (sigma, log scale, mg/L)" = 0.205 (RSE 2.74 %)
  })

  model({
    # Mixture-gated typical clearance. The MIX_FAST_ELIM = 0 case
    # collapses to exp(lcl) = 9.70 L/h (slow eliminator typical
    # value); MIX_FAST_ELIM = 1 multiplies by exp(e_mix_fast_elim_cl)
    # = 21.6 / 9.70 = 2.227 to yield 21.6 L/h (fast eliminator
    # typical value). HIV co-infection applies the linear -17.4 %
    # multiplicative shift to both mixture classes per the Table 2
    # footnote.
    cl <- exp(lcl + e_mix_fast_elim_cl * MIX_FAST_ELIM + etalcl) *
          (WT / 70) ^ e_wt_cl_q *
          (1 + e_hiv_pos_cl * HIV_POS)

    # Central volume: linear -10.3 % shift in female subjects.
    vc <- exp(lvc + etalvc) *
          (WT / 70) ^ e_wt_vc_vp *
          (1 + e_sexf_vc * SEXF)

    # Peripheral volume: allometric only.
    vp <- exp(lvp) * (WT / 70) ^ e_wt_vc_vp

    # Inter-compartmental clearance: allometric only.
    q  <- exp(lq + etalq) * (WT / 70) ^ e_wt_cl_q

    # First-order absorption rate constant.
    ka <- exp(lka + etalka)

    # Absorption lag time.
    tlag <- exp(ltlag + etaltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system: oral depot + central + peripheral1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                  k12 * central - k21 * peripheral1
    alag(depot)       <- tlag

    # Bioavailability (anchor at exp(lfdepot) = 1; IIV via
    # etalfdepot carries the 26.2 % CV reported in Table 2).
    f(depot) <- exp(lfdepot + etalfdepot)

    # Plasma concentration of isoniazid (mg/L).
    Cc <- central / vc

    # Lognormal residual error -- "additive on the natural
    # logarithmic scale" per Wilkins 2011 Methods.
    Cc ~ lnorm(expSd)
  })
}
