Morris_2011_telapristone <- function() {
  description <- paste(
    "Population PK model for telapristone (CDB-4124, a selective",
    "progesterone-receptor antagonist developed for endometriosis and",
    "uterine fibroids) and its active monodemethylated metabolite",
    "CDB-4453 (Morris 2011). Parent is a two-compartment model with",
    "first-order oral absorption (no lag); metabolite is a one-compartment",
    "model with apparent volume V3/F fixed to 1 L for identifiability",
    "(Fmet is not separately identifiable from V3, so the estimated",
    "fmetest is interpreted as the ratio Fmet / V3 in 1/L). A NONMEM",
    "$MIXTURE block splits parent CL/F into a high-CL fast-eliminator",
    "subpopulation (CL/F = 11.6 L/h, population fraction P = 0.251) and",
    "a low-CL slow-eliminator subpopulation (CL/F = 3.34 L/h, P = 0.749);",
    "the mechanism is hypothesized in the Discussion to be polymorphic",
    "CYP3A5 activity but not directly tested. The mixture assignment is",
    "supplied as the binary covariate MIX_FAST_ELIM (1 = fast eliminator,",
    "0 = slow eliminator) drawn per subject from a Bernoulli(0.251). The",
    "only retained clinical covariate is moderate renal impairment",
    "(RENALIMP_MOD), which produces a 74% proportional decrease in the",
    "telapristone absorption rate constant Ka relative to the",
    "healthy / mild-renal-impaired reference cohort."
  )
  reference <- paste(
    "Morris D, Podolski J, Kirsch A, Wiehle R, Fleckenstein L. (2011).",
    "Population Pharmacokinetics of Telapristone (CDB-4124) and its Active",
    "Monodemethylated Metabolite CDB-4453, with a Mixture Model for Total",
    "Clearance. The AAPS Journal 13(4):665-673.",
    "doi:10.1208/s12248-011-9304-7."
  )
  vignette <- "Morris_2011_telapristone"
  paper_specific_etas <- "etalfmetest"

  units <- list(
    time          = "hour",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list(
    MIX_FAST_ELIM = list(
      description        = paste(
        "Per-subject binary mixture-model class indicator for the",
        "telapristone CL/F mixture. 1 = subject classified to the high-CL",
        "fast-eliminator subpopulation (typical CL/F = 11.6 L/h, ~12 h",
        "elimination half-life); 0 = subject classified to the low-CL",
        "slow-eliminator subpopulation (typical CL/F = 3.34 L/h, ~35 h",
        "elimination half-life)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (low-CL slow-eliminator subpopulation; the majority class at 74.9 % of the source cohort)",
      notes              = paste(
        "Not a measured patient covariate. Per-subject latent-class",
        "assignment from the NONMEM $MIXTURE block (Morris 2011 Methods",
        "'Mixture Model' subsection and Equation for the individual",
        "probability IPk per Carlsson et al.). Population probability of",
        "MIX_FAST_ELIM = 1 is the estimated population mixture fraction",
        "Ppop_high = 0.251 (Morris 2011 Table II row 'Probability', RSE",
        "61.0 %; bootstrap 95 % CI 0.038-0.534). The Discussion attributes",
        "the bimodal CL distribution to polymorphic CYP3A5 (functional",
        "CYP3A5 present in 10-40 % of Caucasians, ~50 % of African",
        "Americans, ~33 % of Asians) but does not test the genotype",
        "directly; the mixture indicator is a latent-class label, not a",
        "CYP3A5 genotype. For typical-value simulation set MIX_FAST_ELIM",
        "= 1 (fast eliminator, lower steady-state exposure) or 0 (slow",
        "eliminator, higher exposure; dominates the cohort). For",
        "population simulation, draw MIX_FAST_ELIM ~ Bernoulli(0.251) per",
        "subject. The reference category (= 0 = slow) is chosen so the",
        "binary numerically matches the paper's mixture-indicator",
        "orientation."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment)"
    ),
    RENALIMP_MOD = list(
      description        = paste(
        "Per-subject binary indicator of moderate renal impairment.",
        "1 = Cockcroft-Gault creatinine clearance 30-50 mL/min (the",
        "ZP-006 moderate-renal-impairment stratum); 0 = healthy renal",
        "function or mild renal impairment (CrCl >= 50 mL/min, the",
        "ZP-006 healthy + mild-renal-impairment reference stratum, and",
        "the ZP-005 healthy + hepatic-impairment cohort)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy renal function or mild renal impairment)",
      notes              = paste(
        "Morris 2011 Methods 'Covariate Analysis' and Results: the only",
        "retained covariate from a stepwise forward addition / backward",
        "elimination screen of age, weight, height, BMI, AST, ALT, CrCl,",
        "renal impairment status, and hepatic impairment status. Encoded",
        "as a proportional fractional-change effect on Ka:",
        "Ka_i = Ka_typ * (1 + e_renalimp_mod_ka * RENALIMP_MOD), giving a",
        "74 % decrease in Ka in moderate-renal-impaired subjects. The",
        "Discussion notes that the effect is not directly attributed to",
        "decreased GFR (creatinine clearance itself was screened and did",
        "not meet the inclusion threshold as a continuous covariate) but",
        "to other pathophysiological states associated with chronic renal",
        "impairment, e.g. delayed gastric emptying. The source cohort had",
        "n = 6 moderate-renal-impaired subjects (ZP-006); the remaining",
        "n = 26 subjects share the reference category. Severe renal",
        "impairment / end-stage renal disease subjects were not enrolled."
      ),
      source_name        = "RENAL (Morris 2011 Methods covariate model)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 2L,
    age_range      = "35-66 years (median 52.5)",
    weight_range   = "52.3-102.3 kg (median 71.65)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Adult women across two phase I/II open-label single-dose",
      "telapristone-acetate PK studies. ZP-005 (n = 11; NCT00741273):",
      "n = 4 healthy volunteers + n = 7 moderate hepatic impairment",
      "(Child-Pugh class B) -- each subject received 25 mg followed by",
      "50 mg after 14-day washout. ZP-006 (n = 21; NCT00787618): n = 9",
      "healthy + n = 6 mild renal impairment (CrCl 50-80 mL/min) + n = 6",
      "moderate renal impairment (CrCl 30-50 mL/min) -- single 50 mg",
      "dose."
    ),
    dose_range     = paste(
      "25 mg or 50 mg telapristone acetate, single oral dose under",
      "fasted conditions. Dose-to-amount conversion uses telapristone",
      "molecular weight 505.6 g/mol per Morris 2011 Methods 'Data",
      "Handling' (which cites the AMA USAN telapristone-acetate",
      "datasheet); a 25 mg dose maps to 49,447 nmol and 50 mg to 98,894",
      "nmol of telapristone equivalents in the depot."
    ),
    regions        = "United States (clinical trial sites under 21 CFR Part 56 IRB approval)",
    notes          = paste(
      "Demographics from Morris 2011 Table I. Plasma sampling 0-48 h",
      "post-dose at 22 time points per subject. Of 1,805 concentration",
      "measurements across both studies, 5.0 % (n = 90) were below the",
      "quantification limit of 5 ng/mL and were excluded from the",
      "analysis. Race / ethnicity not reported in Table I; the Discussion",
      "cites the racial distribution of functional CYP3A5 in the",
      "background population but does not stratify the cohort by race."
    )
  )

  ini({
    # Structural PK parameters. Final estimates from Morris 2011 Table II
    # ("Pharmacokinetic Parameter Estimates, Standard Errors, and
    # Bootstrap Results for Final Model"). All `l<base>` values stored as
    # natural log of the linear-scale typical value.

    # Mixture-selected parent clearance. The NONMEM $MIXTURE block
    # estimated two typical CL/F values; both share the single IIV
    # `etalcl` reported in Table II row 'IIV-CL/F'. The subject-level
    # selection is driven by the covariate MIX_FAST_ELIM.
    lcl_pop1 <- log(11.6)
    label("Parent CL/F -- high-CL fast-eliminator subpopulation (L/h)")  # Table II 'CL/F (L/h) - high clearance' (Estimate 11.6, RSE 36.0 %)
    lcl_pop2 <- log(3.34)
    label("Parent CL/F -- low-CL slow-eliminator subpopulation (L/h)")   # Table II 'CL/F (L/h) - low clearance' (Estimate 3.34, RSE 19.1 %)

    # Parent distribution and absorption
    lvc <- log(37.4); label("Parent V2/F -- central volume of distribution (L)")            # Table II 'V2/F (L)' (Estimate 37.4, RSE 10.2 %)
    lvp <- log(120);  label("Parent V4/F -- peripheral volume of distribution (L)")         # Table II 'V4/F (L)' (Estimate 120, RSE 12.4 %)
    lq  <- log(21.9); label("Parent Q/F -- inter-compartmental clearance (L/h)")            # Table II 'Q/F (L/h)' (Estimate 21.9, RSE 14.4 %)
    lka <- log(1.26); label("Parent Ka -- first-order absorption rate constant (1/h)")      # Table II 'Ka (h^-1)' (Estimate 1.26, RSE 7.21 %)

    # Covariate effect: moderate renal impairment on Ka, additive on
    # the log-transformed THETA per Morris 2011 Equation
    # Ka = theta1 * (1 + theta2 * RENAL); theta2 = -0.744 corresponds to
    # a 74 % decrease in absorption rate in moderate-renal-impaired
    # subjects vs the healthy / mild-renal reference.
    e_renalimp_mod_ka <- -0.744
    label("Proportional Ka effect of moderate renal impairment (fractional change)")        # Table II 'Moderate renal impairment on Ka' (Estimate -0.744, RSE 11.6 %)

    # Metabolite (CDB-4453) parameters
    lvc_cdb4453 <- fixed(log(1))
    label("Metabolite V3/F -- apparent volume (L), FIXED for identifiability")              # Table II 'V3/F (L)' = 1 (fixed); Methods 'V3/F was fixed to 1 L'
    lcl_cdb4453 <- log(2.43)
    label("Metabolite CLM/F -- apparent clearance (L/h)")                                   # Table II 'CLM/F (L/h)' (Estimate 2.43, RSE 15.1 %)
    lfmetest    <- log(0.201)
    label("Fmetest = Fmet / V3_metab (1/L)")                                                # Table II 'Fmet est (L^-1)' (Estimate 0.201, RSE 16.8 %)

    # Inter-individual variability. Morris 2011 Methods 'Base
    # Pharmacokinetic Model' states "IIV ... modelled assuming a
    # log-normal distribution ... The magnitude of IIV was expressed as
    # coefficient of variation (%CV)." Table II IIV column lists each
    # entry as "variance (%CV)". The variances in Table II are the
    # omega^2 values reported by NONMEM, which equal log(1 + CV^2) when
    # the log-normal CV is expressed as a fraction. Quick check on Ka:
    # reported variance 0.280, reported CV 52.9 % -> log(1 + 0.529^2) =
    # 0.249 vs reported 0.280; the small discrepancy reflects that the
    # paper reports approximate CV%, so we use the variance column
    # directly as the omega^2 entered into nlmixr2.
    etalcl       ~ 0.200       # Table II 'IIV-CL/F' variance (CV 44.7 %); shared across both mixture subpopulations per Table II row
    etalvc       ~ 0.443       # Table II 'IIV-V2/F' variance (CV 66.6 %)
    etalvp       ~ 0.426       # Table II 'IIV-V4/F' variance (CV 65.3 %)
    etalq        ~ 0.457       # Table II 'IIV-Q/F'  variance (CV 67.6 %)
    etalka       ~ 0.280       # Table II 'IIV-Ka'   variance (CV 52.9 %)
    etalfmetest  ~ 0.196       # Table II 'IIV-Fmet est' variance (CV 44.3 %)

    # Residual variability. Morris 2011 Methods 'Base Pharmacokinetic
    # Model' states "Residual variability was modelled using a log error
    # model": log(Cij) = log(Cpred,ij) + epsilon. In the linear space
    # this is equivalent to a proportional residual error structure
    # (Cij ~ Cpred,ij * exp(epsilon)), so the NONMEM "additive on
    # log-scale" sigma maps onto nlmixr2's prop() residual SD.
    propSd          <- 0.378
    label("Telapristone proportional residual SD (fraction; log-error model)")              # Table II 'Telapristone (nmol/L)' RV (Estimate 0.378, RSE 15.8 %)
    propSd_cdb4453  <- 0.109
    label("CDB-4453 proportional residual SD (fraction; log-error model)")                  # Table II 'CDB-4453 (nmol/L)' RV (Estimate 0.109, RSE 18.6 %)
  })

  model({
    # -----------------------------------------------------------------
    # 1. Mixture-class indicators (per-subject latent class)
    # -----------------------------------------------------------------
    # mix_fast / mix_slow are mutually-exclusive 0/1 selectors derived
    # from the per-subject covariate MIX_FAST_ELIM (~ Bernoulli(0.251)
    # at simulation time).
    mix_fast <- MIX_FAST_ELIM
    mix_slow <- 1 - mix_fast

    # -----------------------------------------------------------------
    # 2. Individual PK parameters
    # -----------------------------------------------------------------
    # Mixture-selected typical clearance on the log scale; the single
    # shared etalcl applies to whichever subpopulation the subject is in.
    lcl_select <- mix_fast * lcl_pop1 + mix_slow * lcl_pop2
    cl <- exp(lcl_select + etalcl)

    # Renal-impairment effect on Ka: proportional fractional change
    # encoded as `Ka = Ka_typ * (1 + theta_renal * RENALIMP_MOD)` per
    # Morris 2011 covariate-analysis equation. The reference Ka in
    # healthy / mild-renal subjects (RENALIMP_MOD = 0) reduces to
    # exp(lka + etalka); moderate-renal subjects (RENALIMP_MOD = 1) carry
    # the additional (1 - 0.744) = 0.256 factor (a 74 % decrease).
    ka <- exp(lka + etalka) * (1 + e_renalimp_mod_ka * RENALIMP_MOD)

    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Metabolite parameters. V3/F = 1 L by fixed-anchor identifiability;
    # no IIV on the metabolite CL or V per Morris 2011 Table II.
    vc_cdb4453 <- exp(lvc_cdb4453)
    cl_cdb4453 <- exp(lcl_cdb4453)
    fmetest    <- exp(lfmetest + etalfmetest)

    # -----------------------------------------------------------------
    # 3. Micro-constants
    # -----------------------------------------------------------------
    kel  <- cl  / vc
    k12  <- q   / vc
    k21  <- q   / vp
    kelm <- cl_cdb4453 / vc_cdb4453

    # -----------------------------------------------------------------
    # 4. ODE system
    # -----------------------------------------------------------------
    # Parent telapristone: depot -> central <-> peripheral1, first-order
    # elimination from central. Metabolite CDB-4453: 1-compartment with
    # formation rate Fmet * CL_parent * Cc and elimination CLM * Cc_metab.
    # Since fmetest = Fmet / V3 and V3 is fixed to 1 L, the numeric
    # equivalence Fmet = fmetest * V3 = fmetest holds here.
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(central_cdb4453)  <-  fmetest * vc_cdb4453 * cl * (central / vc) - kelm * central_cdb4453

    # -----------------------------------------------------------------
    # 5. Observation variables (nmol/L)
    # -----------------------------------------------------------------
    Cc          <- central          / vc
    Cc_cdb4453  <- central_cdb4453  / vc_cdb4453

    # -----------------------------------------------------------------
    # 6. Residual error -- proportional in linear space (log-additive in NONMEM)
    # -----------------------------------------------------------------
    Cc         ~ prop(propSd)
    Cc_cdb4453 ~ prop(propSd_cdb4453)
  })
}
