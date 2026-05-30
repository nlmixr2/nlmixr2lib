Shoji_2017_fosdagrocorat_p1np <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) model for serum amino-terminal propeptide of type I collagen (P1NP) bone-formation biomarker following once-daily oral fosdagrocorat (PF-04171327, a dissociated agonist of the glucocorticoid receptor) or oral prednisone comparator in adults with rheumatoid arthritis on background methotrexate (Shoji 2017). A virtual K-PD depot for the drug (zero-order Input mg/week, first-order elimination KDE) feeds a sigmoid Emax inhibition of biomarker synthesis (Hill coefficient fixed to 1); the synthesis rate carries an empirical dose-and-time-dependent rebound multiplier and an additive linear placebo-period slope captures the methotrexate-only time trend."
  reference <- "Shoji S, Suzuki A, Conrado DJ, Peterson MC, Hey-Hadavi J, McCabe D, Rojo R, Tammara BK. Dissociated Agonist of Glucocorticoid Receptor or Prednisone for Active Rheumatoid Arthritis: Effects on P1NP and Osteocalcin Pharmacodynamics. CPT Pharmacometrics Syst Pharmacol. 2017;6(7):439-448. doi:10.1002/psp4.12201"
  vignette <- "Shoji_2017_fosdagrocorat_p1np"
  units <- list(time = "week", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject assigned once-daily oral dose (mg) of fosdagrocorat or prednisone driving the dose-and-time-dependent rebound term on biomarker synthesis. Set to 0 for placebo (methotrexate only).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject (time-fixed) assigned q.d. dose level in milligrams. Enters the model only via the empirical rebound multiplier on synthesis 1 + RBmax * DOSE * t / (T50 + t); the drug input to the K-PD depot is supplied separately via dosing events with rate 7 * DOSE mg/week (Shoji 2017 Methods, K-PD model). Cohort levels in the source trial: 0 (placebo), 1, 5, 10, 15 (fosdagrocorat q.d. mg) and 5, 10 (prednisone q.d. mg).",
      source_name        = "DOSE"
    ),
    DRUG_PRED = list(
      description        = "1 = subject is in the prednisone comparator arm; 0 = subject is in the fosdagrocorat arm or placebo. Per-subject (time-fixed) categorical indicator switching the drug-elimination rate KDE and the sigmoid-Emax inhibition parameters (Imax, EDK50) between the two drugs.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fosdagrocorat or placebo)",
      notes              = "Shoji 2017 Table 2 reports separate KDE, Imax, and EDK50 estimates for fosdagrocorat and prednisone; the rebound parameters (RBmax, T50) and the response-side parameters (Kd, BL, SLP) are shared between the two drugs. The model selects the active parameter set by adding a log-ratio offset (dlkel_pred / dledk50_pred) and a logit-difference (dlogitimax_pred) when DRUG_PRED = 1, recovering the published Prednisone typical values exactly. For placebo subjects (no dosing events into the K-PD depot) DRUG_PRED is informational only.",
      source_name        = "DRUG_PRED"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 321L,
    n_observations   = 4837L,
    n_studies        = 1L,
    age_range        = "18-80 years (median 56)",
    weight_range     = "36.6-144 kg (median 71.0)",
    bmi_range        = "15.9-51.0 kg/m^2 (median 27.1)",
    sex_female_pct   = 80,
    race_ethnicity   = c(White = 87.5, Black = 2.2, Asian = 7.5, Other = 2.8),
    disease_state    = "Adults with active rheumatoid arthritis on background methotrexate (clinicaltrials.gov NCT01393639, phase II randomized double-blind parallel-group trial).",
    dose_range       = "Fosdagrocorat 1, 5, 10, or 15 mg q.d. orally; prednisone 5 or 10 mg q.d. orally; placebo q.d. All for 8 weeks followed by a 4-week taper (weeks 9-10 every other day at reduced dose; weeks 11-12 every 3 days) per protocol.",
    regions          = "Global multi-regional Phase II trial.",
    bl_p1np_ng_per_mL = "Median 48.5 ng/mL [min 10.9, max 184]",
    bl_oc_ng_per_mL   = "Median 22.4 ng/mL [min 4.85, max 85.9]",
    notes            = "Baseline demographics from Shoji 2017 Table 1 (full intent-to-treat n = 323; n = 321 used in this K-PD analysis after exclusion of 2 patients with missing baseline biomarker concentrations). Race breakdown is the all-treatments column (white/black/Asian/other = 281/7/24/9). Background therapy: methotrexate in all arms."
  )

  ini({
    # ===================================================================
    # Drug PK (virtual K-PD depot) -- Shoji 2017 Table 2 (P1NP)
    # KDE differs by drug. Reparameterized as base (fosdagrocorat KDE) plus
    # a log-ratio offset for prednisone so that a single etalkel pairs with
    # lkel and reproduces the published values for both arms:
    #   typical KDE(fos)  = exp(lkel)             = 0.597 /week
    #   typical KDE(pred) = exp(lkel + dlkel_pred) = 0.535 /week
    # ===================================================================
    lkel       <- log(0.597);          label("Effect-compartment elimination rate KDE for fosdagrocorat (1/week)")        # Shoji 2017 Table 2 P1NP: KDE Fosdagrocorat = 0.597 /week (RSE 17.8%)
    dlkel_pred <- log(0.535 / 0.597);  label("Log-ratio offset for prednisone KDE: log(KDE_pred / KDE_fos) (unitless)")    # Shoji 2017 Table 2 P1NP: KDE Prednisone   = 0.535 /week (RSE 28.9%)

    # ===================================================================
    # Biomarker (response compartment) -- Shoji 2017 Table 2 (P1NP)
    # ===================================================================
    lkd  <- log(0.609); label("First-order degradation rate Kd of P1NP (1/week)")  # Shoji 2017 Table 2 P1NP: Kd = 0.609 /week (RSE 17.5%)
    lrbase  <- log(47.0);  label("Baseline P1NP concentration BL (ng/mL)")             # Shoji 2017 Table 2 P1NP: BL = 47.0 ng/mL (RSE 2.83%)

    # Sigmoid Emax inhibition (Hill coefficient fixed to 1 -- paper Discussion:
    # c = 0.920 led to estimation instability, c fixed = 1.5-3 underpredicted).
    # Imax bounded in (0, 1) -> encoded on the logit scale; reparameterized as
    # base (fosdagrocorat) plus a logit-difference offset for prednisone.
    logitimax       <- log(0.751 / (1 - 0.751));                                  label("Logit of maximum fractional inhibition Imax for fosdagrocorat (unitless)") # Shoji 2017 Table 2 P1NP: Imax Fosdagrocorat = 0.751 (RSE 4.85%)
    dlogitimax_pred <- log(0.754 / (1 - 0.754)) - log(0.751 / (1 - 0.751));        label("Logit-difference offset for prednisone Imax (unitless)")                    # Shoji 2017 Table 2 P1NP: Imax Prednisone   = 0.754 (RSE 8.89%)

    # EDK50 reparameterized as base (fosdagrocorat) plus log-ratio offset for
    # prednisone so etaledk50 pairs with ledk50.
    ledk50       <- log(40.1);         label("Effect-compartment infusion rate IR producing 50% of Imax, fosdagrocorat (mg/week)") # Shoji 2017 Table 2 P1NP: EDK50 Fosdagrocorat = 40.1 mg/week (RSE 17.8%)
    dledk50_pred <- log(45.9 / 40.1);  label("Log-ratio offset for prednisone EDK50: log(EDK50_pred / EDK50_fos) (unitless)")      # Shoji 2017 Table 2 P1NP: EDK50 Prednisone   = 45.9 mg/week (RSE 30.2%)

    hill <- fixed(1); label("Hill coefficient gamma for the sigmoid Emax inhibition (unitless)")  # Shoji 2017 Table 2 P1NP: c FIX = 1 (paper Discussion)

    # Empirical dose-and-time-dependent rebound multiplier on synthesis,
    # common for fosdagrocorat and prednisone at the same q.d. dose.
    lrbmax <- log(0.0479); label("Maximum rebound effect on synthesis rate (1/mg)")  # Shoji 2017 Table 2 P1NP: RBmax = 0.0479 /mg (RSE 17.6%)
    lt50   <- log(1.13);   label("Time to 50% of maximum rebound effect (weeks)")    # Shoji 2017 Table 2 P1NP: T50   = 1.13 weeks (RSE 42.6%)

    # Additive linear placebo-period slope (observation-level, additive to R(t))
    slp <- 0.162; label("Linear placebo-period slope for P1NP (ng/mL per week)")  # Shoji 2017 Table 2 P1NP: SLP = 0.162 ng/mL/week (RSE 67.3%)

    # ===================================================================
    # IIV -- Shoji 2017 Table 2 (P1NP)
    # Reported "%CV [g_X]" follows the footnote convention %CV = sqrt(omega^2) * 100
    # (the paper's short-hand: omega^2 = (CV/100)^2). Off-diagonals derived from
    # the reported correlations q[g_i, g_j] = cov / (sd_i * sd_j):
    #   omega^2[g_KDE]    = 0.95^2  = 0.9025  (95.0 %CV)
    #   omega^2[g_EDK50]  = 0.655^2 = 0.4290  (65.5 %CV)
    #   omega^2[g_BL]     = 0.466^2 = 0.2172  (46.6 %CV)
    #   cov(g_KDE,  g_EDK50) = -0.312 * sqrt(0.9025 * 0.4290) = -0.1941
    #   cov(g_BL,   g_KDE)   = -0.316 * sqrt(0.2172 * 0.9025) = -0.1401
    #   cov(g_BL,   g_EDK50) = -0.410 * sqrt(0.2172 * 0.4290) = -0.1251
    # The same etalkel is added to whichever drug's KDE is active (single eta
    # on the K-PD depot per subject); analogous for etaledk50.
    # ===================================================================
    etalkel + etaledk50 + etalrbase ~ c(
       0.9025,
      -0.1941,  0.4290,
      -0.1401, -0.1251,  0.2172
    )

    # Additive (linear-scale) eta on SLP -- Shoji 2017 Methods: "For SLP only,
    # the parameter for the i th individual is described with an additive error
    # model SLP_i = theta_SLP + eta_SLP, eta_SLP ~ N(0, omega^2_SLP).
    # Table 2 reports IIV SD = 0.928 ng/mL/week, so omega^2 = 0.928^2 = 0.8612.
    etaslp ~ 0.8612

    # ===================================================================
    # Residual variability -- proportional in linear space; Shoji 2017 Table 2
    # reports "Residual variability %CV [e] = 15.2", with the footnote
    # %CV = sqrt(sigma^2) * 100, so propSd = sigma = 0.152.
    # ===================================================================
    propSd <- 0.152; label("Proportional residual error on P1NP (fraction)")  # Shoji 2017 Table 2 P1NP: e = 15.2 %CV (RSE 0.962%)
  })

  model({
    # ---- Drug switching: select the active KDE / logit-Imax / EDK50 by arm. ----
    # DRUG_PRED = 0 -> fosdagrocorat; DRUG_PRED = 1 -> prednisone. Placebo
    # subjects have no dose entering the K-PD depot so the choice is
    # informational. Single etalkel / etaledk50 add to whichever active value.
    lkel_active      <- lkel      + dlkel_pred      * DRUG_PRED
    logitimax_active <- logitimax + dlogitimax_pred * DRUG_PRED
    ledk50_active    <- ledk50    + dledk50_pred    * DRUG_PRED

    # ---- Individual parameters ----
    kel   <- exp(lkel_active + etalkel)
    imax  <- 1 / (1 + exp(-logitimax_active))
    edk50 <- exp(ledk50_active + etaledk50)
    kd    <- exp(lkd)
    rbase    <- exp(lrbase + etalrbase)
    rbmax <- exp(lrbmax)
    t50   <- exp(lt50)
    slp_i <- slp + etaslp

    # ---- Steady-state synthesis at baseline (no drug): Ks = BL * Kd ----
    ks <- rbase * kd

    # ---- K-PD depot (drug). Drug input enters as zero-order infusion via
    # dosing events (rate = 7 * q.d. dose in mg/week for continuous q.d.
    # dosing per Shoji 2017 Methods); first-order elimination at rate KDE.
    d/dt(depot_kpd) <- -kel * depot_kpd
    ir <- kel * depot_kpd

    # ---- Empirical rebound multiplier on synthesis (Shoji 2017 Results,
    # rebound-effect equation). DOSE is the per-subject q.d. dose covariate
    # (mg); t is time after the first dose (rxode2 simulation time).
    rebound <- 1 + rbmax * DOSE * t / (t50 + t)

    # ---- Sigmoid Emax inhibition (Hill coefficient fixed to 1) ----
    inhibition <- 1 - imax * ir^hill / (edk50^hill + ir^hill)

    # ---- Effect compartment (P1NP biomarker); initial condition = BL ----
    d/dt(effect) <- ks * rebound * inhibition - kd * effect
    effect(0)    <- rbase

    # ---- Observation: P1NP = R(t) + SLP_i * t (placebo-period linear trend),
    # with proportional residual error in linear space (Shoji 2017 Methods:
    # Y_ij = F_ij with proportional intra-individual variability e_ij).
    P1NP <- effect + slp_i * t
    P1NP ~ prop(propSd)
  })
}
