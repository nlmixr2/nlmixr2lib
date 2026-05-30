Shoji_2017_fosdagrocorat_oc <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) model for serum osteocalcin (OC) bone-formation biomarker following once-daily oral fosdagrocorat (PF-04171327, a dissociated agonist of the glucocorticoid receptor) or oral prednisone comparator in adults with rheumatoid arthritis on background methotrexate (Shoji 2017). Sister model to Shoji_2017_fosdagrocorat_p1np: identical K-PD structure (virtual K-PD depot with zero-order Input mg/week and first-order KDE; sigmoid Emax inhibition of biomarker synthesis with Hill coefficient fixed to 1; empirical dose-and-time-dependent rebound multiplier; additive placebo-period slope). For the OC fit Shoji 2017 fixed KDE to the P1NP-derived estimates and fixed Imax to 1 for both drugs, and used independent (not block) IIV on KDE, EDK50, and BL."
  reference <- "Shoji S, Suzuki A, Conrado DJ, Peterson MC, Hey-Hadavi J, McCabe D, Rojo R, Tammara BK. Dissociated Agonist of Glucocorticoid Receptor or Prednisone for Active Rheumatoid Arthritis: Effects on P1NP and Osteocalcin Pharmacodynamics. CPT Pharmacometrics Syst Pharmacol. 2017;6(7):439-448. doi:10.1002/psp4.12201"
  vignette <- "Shoji_2017_fosdagrocorat_oc"
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
      description        = "1 = subject is in the prednisone comparator arm; 0 = subject is in the fosdagrocorat arm or placebo. Per-subject (time-fixed) categorical indicator switching the drug-elimination rate KDE and the sigmoid-Emax EDK50 between the two drugs.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fosdagrocorat or placebo)",
      notes              = "Shoji 2017 Table 2 reports separate KDE and EDK50 estimates for fosdagrocorat and prednisone; Imax is fixed to 1 for both drugs in the OC fit, and the rebound parameters (RBmax, T50) and response-side parameters (Kd, BL, SLP) are shared between the drugs. The model selects the active parameter set by adding a fixed log-ratio offset (dlkel_pred / dledk50_pred) when DRUG_PRED = 1.",
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
    notes            = "Baseline demographics from Shoji 2017 Table 1 (full intent-to-treat n = 323; n = 321 used in this K-PD analysis after exclusion of 2 patients with missing baseline biomarker concentrations). Race breakdown is the all-treatments column (white/black/Asian/other = 281/7/24/9). Background therapy: methotrexate in all arms. The total observation count (4837) is across both biomarkers; OC was modelled on the same set of subjects as P1NP using the developed P1NP structural model."
  )

  ini({
    # ===================================================================
    # Drug PK (virtual K-PD depot) -- Shoji 2017 Table 2 (OC).
    # KDE for both drugs were FIXED to the P1NP-fit estimates because the
    # OC-only fit was unstable (high RSE). Reparameterized as base
    # (fosdagrocorat) plus a log-ratio offset for prednisone so etalkel
    # pairs with lkel.
    # ===================================================================
    lkel       <- fixed(log(0.597));         label("Effect-compartment elimination rate KDE for fosdagrocorat (1/week)")        # Shoji 2017 Table 2 OC: KDE Fosdagrocorat = 0.597 /week FIX (from P1NP)
    dlkel_pred <- fixed(log(0.535 / 0.597)); label("Log-ratio offset for prednisone KDE: log(KDE_pred / KDE_fos) (unitless)")    # Shoji 2017 Table 2 OC: KDE Prednisone   = 0.535 /week FIX (from P1NP)

    # ===================================================================
    # Biomarker (response compartment) -- Shoji 2017 Table 2 (OC).
    # ===================================================================
    lkd <- log(0.939); label("First-order degradation rate Kd of osteocalcin (1/week)")  # Shoji 2017 Table 2 OC: Kd = 0.939 /week (RSE 24.6%)
    lrbase <- log(22.2);  label("Baseline osteocalcin concentration BL (ng/mL)")             # Shoji 2017 Table 2 OC: BL = 22.2 ng/mL (RSE 2.58%)

    # Sigmoid Emax inhibition -- Imax fixed to 1 for both drugs because the
    # estimate was close to 1 in unconstrained OC fits (paper Discussion).
    # A single shared scalar imax (no drug switching needed because both
    # drugs share the value 1) and a per-drug EDK50 with log-ratio offset.
    imax <- fixed(1); label("Maximum fractional inhibition Imax (unitless, fixed to 1 for both drugs)")  # Shoji 2017 Table 2 OC: Imax Fosdagrocorat = Imax Prednisone = 1 FIX

    ledk50       <- log(148);          label("Effect-compartment infusion rate IR producing 50% of Imax, fosdagrocorat (mg/week)") # Shoji 2017 Table 2 OC: EDK50 Fosdagrocorat = 148 mg/week (RSE 6.68%)
    dledk50_pred <- log(122 / 148);    label("Log-ratio offset for prednisone EDK50: log(EDK50_pred / EDK50_fos) (unitless)")      # Shoji 2017 Table 2 OC: EDK50 Prednisone   = 122 mg/week (RSE 11.3%)

    hill <- fixed(1); label("Hill coefficient gamma for the sigmoid Emax inhibition (unitless)")  # Shoji 2017 Table 2 OC: c FIX = 1 (paper Discussion; same as P1NP fit)

    # Empirical dose-and-time-dependent rebound multiplier on synthesis,
    # common for fosdagrocorat and prednisone at the same q.d. dose.
    lrbmax <- log(0.0276); label("Maximum rebound effect on synthesis rate (1/mg)")  # Shoji 2017 Table 2 OC: RBmax = 0.0276 /mg (RSE 21.5%)
    lt50   <- log(2.24);   label("Time to 50% of maximum rebound effect (weeks)")    # Shoji 2017 Table 2 OC: T50   = 2.24 weeks (RSE 54.8%)

    # Additive linear placebo-period slope (observation-level, additive to R(t))
    slp <- 0.0675; label("Linear placebo-period slope for osteocalcin (ng/mL per week)")  # Shoji 2017 Table 2 OC: SLP = 0.0675 ng/mL/week (RSE 53.6%)

    # ===================================================================
    # IIV -- Shoji 2017 Table 2 (OC). Independent (not block) IIVs were used
    # because the joint estimation was unstable (paper Methods). Reported
    # "%CV [g_X]" follows the same footnote convention as P1NP
    # (%CV = sqrt(omega^2) * 100, i.e. omega^2 = (CV/100)^2):
    #   omega^2[g_KDE]   = 1.23^2  = 1.5129  (123  %CV)
    #   omega^2[g_EDK50] = 0.255^2 = 0.0650  (25.5 %CV)
    #   omega^2[g_BL]    = 0.436^2 = 0.1901  (43.6 %CV)
    # The same etalkel is added to whichever drug's KDE is active.
    # ===================================================================
    etalkel   ~ 1.5129
    etaledk50 ~ 0.0650
    etalrbase    ~ 0.1901

    # Additive (linear-scale) eta on SLP -- Shoji 2017 Methods (same
    # convention as P1NP). Table 2 reports IIV SD = 0.338 ng/mL/week,
    # so omega^2 = 0.338^2 = 0.1142.
    etaslp ~ 0.1142

    # ===================================================================
    # Residual variability -- proportional in linear space; Shoji 2017 Table 2
    # OC reports "Residual variability %CV [e] = 14.1", with the same
    # footnote convention as P1NP, so propSd = sigma = 0.141.
    # ===================================================================
    propSd <- 0.141; label("Proportional residual error on osteocalcin (fraction)")  # Shoji 2017 Table 2 OC: e = 14.1 %CV (RSE 0.887%)
  })

  model({
    # ---- Drug switching: select the active KDE / EDK50 by arm. Imax is
    # shared (= 1 for both drugs in the OC fit). ----
    lkel_active   <- lkel   + dlkel_pred   * DRUG_PRED
    ledk50_active <- ledk50 + dledk50_pred * DRUG_PRED

    # ---- Individual parameters ----
    kel   <- exp(lkel_active + etalkel)
    edk50 <- exp(ledk50_active + etaledk50)
    kd    <- exp(lkd)
    rbase    <- exp(lrbase + etalrbase)
    rbmax <- exp(lrbmax)
    t50   <- exp(lt50)
    slp_i <- slp + etaslp

    # ---- Steady-state synthesis at baseline (no drug): Ks = BL * Kd ----
    ks <- rbase * kd

    # ---- K-PD depot (drug). Drug input enters as zero-order infusion via
    # dosing events (rate = 7 * q.d. dose in mg/week); first-order
    # elimination at rate KDE.
    d/dt(depot_kpd) <- -kel * depot_kpd
    ir <- kel * depot_kpd

    # ---- Empirical rebound multiplier on synthesis ----
    rebound <- 1 + rbmax * DOSE * t / (t50 + t)

    # ---- Sigmoid Emax inhibition (Hill coefficient fixed to 1) ----
    inhibition <- 1 - imax * ir^hill / (edk50^hill + ir^hill)

    # ---- Effect compartment (OC biomarker); initial condition = BL ----
    d/dt(effect) <- ks * rebound * inhibition - kd * effect
    effect(0)    <- rbase

    # ---- Observation: OC = R(t) + SLP_i * t with proportional residual error ----
    OC <- effect + slp_i * t
    OC ~ prop(propSd)
  })
}
