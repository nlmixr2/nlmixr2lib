Pierre_2017_morphine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for IV morphine and its",
    "primary glucuronide metabolite morphine-3-glucuronide (M3G) in 14",
    "healthy adults and 7 patients with biopsy-confirmed nonalcoholic",
    "steatohepatitis (NASH) following a single 5 mg morphine sulfate IV",
    "infusion (Pierre 2017). Morphine is described by a three-compartment",
    "disposition (central + two peripherals) with parallel renal (CL_M_R)",
    "and non-renal (CL_M_NR) clearances; the entire non-renal clearance is",
    "assumed to lead to M3G formation via a single liver transit compartment",
    "with first-order rate constant k_trans. M3G is described by a",
    "one-compartment model with a single total clearance (CL_M3G).",
    "Cumulative urinary morphine and M3G amounts are tracked as",
    "elimination-amount compartments. Total body weight enters all CL/Q",
    "and V parameters a priori with fixed allometric exponents (0.75 and",
    "1, respectively) referenced to 70 kg. The NASH severity score (NASF;",
    "combined NAFLD activity score and fibrosis staging, 0-12) is the only",
    "additional covariate retained in the final model; it acts on M3G",
    "clearance through a linear effect on the natural logarithm of",
    "(NASF / 4) for NASF >= 4 and is identically zero for NASF < 4 so that",
    "healthy and benign-NAFLD subjects (NASF < 5) recover the typical CL_M3G."
  )
  reference <- paste(
    "Pierre V, Johnston CK, Ferslew BC, Brouwer KLR, Gonzalez D.",
    "Population Pharmacokinetics of Morphine in Patients With",
    "Nonalcoholic Steatohepatitis (NASH) and Healthy Adults.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(5):331-339.",
    "doi:10.1002/psp4.12185.",
    sep = " "
  )
  vignette <- "Pierre_2017_morphine"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters all morphine and M3G CL/Q and V",
        "parameters a priori via the allometric relationships CL_i =",
        "CL_70kg * (WT / 70)^0.75 and V_i = V_70kg * (WT / 70)^1.0 (paper",
        "Methods, 'Covariate analysis', Eqs. 1 and 2). Reference weight 70 kg.",
        "Allometric exponents are fixed at the canonical values 0.75 (CL/Q)",
        "and 1 (V) and are not estimated."
      ),
      source_name        = "weight (paper Methods 'Covariate analysis', Eqs. 1-2)"
    ),
    NASF = list(
      description        = paste(
        "NASH severity score (NAFLD activity score plus fibrosis",
        "staging; integer 0-12)."
      ),
      units              = "(count, 0-12)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. The NAFLD activity score (NAS) sums",
        "steatosis (0-3), hepatocyte ballooning (0-2), and lobular",
        "inflammation (0-3) for 0-8 points; the fibrosis staging score",
        "(0 = absent, 1 = perisinusoidal/pericellular, 2 = periportal,",
        "3 = bridging, 4 = cirrhosis) adds 0-4 points; total NASF is 0-12.",
        "Healthy subjects in Pierre 2017 were assigned NASF = 0; observed",
        "NASH NASF values were 4, 5, 7, and 8. Used with a linear effect on",
        "the natural logarithm of (NASF / 4) for NASF >= 4 (paper Methods",
        "'Covariate analysis', Eq. 5); for NASF < 4 the covariate",
        "contribution is set to zero and individual CL_M3G equals the",
        "typical value. Source: paper Methods 'Covariate analysis' and",
        "references 31, 32."
      ),
      source_name        = "NASF (paper Methods 'Covariate analysis')"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21,
    n_studies      = 1,
    age_range      = "20-63 years",
    age_median     = "45 years",
    weight_range   = "52-128 kg",
    weight_median  = "78 kg",
    sex_female_pct = 52.4,
    race_ethnicity = NULL,
    disease_state  = paste(
      "14 healthy adults and 7 adults with biopsy-confirmed",
      "nonalcoholic steatohepatitis (NASH). All NASH subjects had",
      "total body weight > 70 kg and 4/7 had BMI > 30 kg/m^2;",
      "NASH NASF severity scores were 4 (n=1), 5 (n=2), 7 (n=3),",
      "and 8 (n=1). Healthy subjects were assigned NASF = 0. None of",
      "the participants exhibited overt renal dysfunction (median",
      "creatinine clearance 118 mL/min in healthy, 141 mL/min in NASH)."
    ),
    dose_range     = paste(
      "Single 5 mg morphine sulfate IV infused over 5 min, administered",
      "2 hr after a standardized 23.9 g fat meal consumed over 30 min.",
      "Doses were converted to nanomoles of morphine free base for",
      "modeling. Simulation cohort uses 10 mg morphine sulfate",
      "(approximately 13,178 nmol free base) infused over 10 min every",
      "4 hr for 24 hr."
    ),
    regions        = "United States (University of North Carolina at Chapel Hill)",
    notes          = paste(
      "Demographics from Pierre 2017 Table 1. NONMEM 7.3 with ADVAN13",
      "(stiff ODE solver), iterative two-stage followed by Monte Carlo",
      "importance sampling with Laplacian g-eta interaction. Serum samples",
      "(315 total) at 15 timepoints between predose and 480 min; two urine",
      "pools per subject over 4-hr intervals during the 8-hr post-dose",
      "sampling window. M3 method used for BLQ handling (23% of morphine",
      "and <1% of M3G serum samples were BLQ). Lower limits of",
      "quantification: morphine serum 8.24 nM, M3G serum 5.42 nM,",
      "morphine urine 82.4 nM, M3G urine 542 nM. M6G was measured but not",
      "included in the model because 42% of M6G samples were BLQ."
    )
  )

  ini({
    # ============================================================
    # Allometric exponents on weight -- fixed a priori at canonical
    # 0.75 (CL/Q) and 1 (V) values per Methods 'Covariate analysis',
    # Eqs. 1-2. Applied to all morphine and M3G CL/Q and V parameters.
    # ============================================================
    e_wt_cl_q  <- fixed(0.75)
    label("Allometric exponent on all CL/Q parameters (unitless, fixed)")
    # Paper Methods 'Covariate analysis': power exponent 0.75 fixed for
    # all CL and Q parameters (CL_M_NR, CL_M_R, Q_P1, Q_P2, CL_M3G).
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on all V parameters (unitless, fixed)")
    # Paper Methods 'Covariate analysis': power exponent 1 fixed for all
    # V parameters (V_M, V_P1, V_P2, V_M3G).

    # ============================================================
    # Morphine structural disposition -- paper Table 2 'Final model'.
    # Reference values are for a 70 kg subject.
    # ============================================================
    lcl_nonren <- log(44.1)
    label("Morphine non-renal clearance at 70 kg (L/h)")
    # Table 2: CL_M,NR = 44.1 L/h (RSE 9%); leads exclusively to M3G
    # formation (paper Results: 'CL_M,NR was assumed to lead exclusively
    # to the formation of M3G').

    lcl_renal <- log(6.32)
    label("Morphine renal clearance at 70 kg (L/h)")
    # Table 2: CL_M,R = 6.32 L/h (RSE 11%).

    lvc <- log(9.41)
    label("Morphine central volume V_M at 70 kg (L)")
    # Table 2: V_M = 9.41 L (RSE 13%).

    lvp <- log(108)
    label("Morphine first peripheral volume V_P1 at 70 kg (L)")
    # Table 2: V_P1 = 108 L (RSE 37%).

    lq <- log(67.1)
    label("Morphine first intercompartmental clearance Q_P1 at 70 kg (L/h)")
    # Table 2: Q_P1 = 67.1 L/h (RSE 17%).

    lvp2 <- log(50.7)
    label("Morphine second peripheral volume V_P2 at 70 kg (L)")
    # Table 2: V_P2 = 50.7 L (RSE 21%).

    lq2 <- log(83.4)
    label("Morphine second intercompartmental clearance Q_P2 at 70 kg (L/h)")
    # Table 2: Q_P2 = 83.4 L/h (RSE 15%).

    lktrans <- log(14.4)
    label("Liver transit rate constant k_trans (1/h)")
    # Table 2: k_trans = 14.4 1/h (RSE 14%); rate constant from the
    # liver transit compartment into the M3G central compartment.

    # ============================================================
    # M3G structural disposition -- paper Table 2 'Final model'.
    # ============================================================
    lcl_m3g <- log(7.32)
    label("M3G total clearance at 70 kg, NASF < 4 (L/h)")
    # Table 2: CL_M3G = 7.32 L/h (RSE 10%); reported as 'Total clearance
    # of M3G'.

    lvc_m3g <- log(9.51)
    label("M3G central volume at 70 kg (L)")
    # Table 2: V_M3G = 9.51 L (RSE 15%).

    # ============================================================
    # NASF covariate effect on M3G clearance -- paper Eq. 5.
    # For NASF >= 4: CL_M3G_i = CL_M3G_pop * (1 + e_nasf_cl_m3g * ln(NASF / 4))
    # For NASF <  4: CL_M3G_i = CL_M3G_pop (no covariate contribution).
    # ============================================================
    e_nasf_cl_m3g <- -0.628
    label("NASF linear effect on log(NASF/4) for CL_M3G (unitless)")
    # Table 2: NASF on CL_M3G = -0.628 (RSE 26%); applied for NASF >= 4
    # per Eq. 5. Inverse association: higher NASF reduces CL_M3G.

    # ============================================================
    # IIV -- log-normal variance parameterized as omega^2 = log(1 + CV^2).
    # CV% values from paper Table 2 'IIV, %CV (RSE, %)' column. The
    # correlation between IIV on CL_M3G and IIV on V_M3G is 0.751.
    # No IIV was retained for CL_M_R, V_P1, Q_P1, V_P2, Q_P2, or k_trans
    # (Table 2 shows '-' in those rows).
    # ============================================================
    etalcl_nonren ~ log(1 + 0.316^2)
    # Table 2 IIV on CL_M_NR = 31.6% CV (RSE 80%); omega^2 = log(1 + 0.316^2).

    etalvc ~ log(1 + 0.423^2)
    # Table 2 IIV on V_M = 42.3% CV (RSE 55%); omega^2 = log(1 + 0.423^2).

    etalcl_m3g + etalvc_m3g ~ c(
      log(1 + 0.345^2),
      0.751 * sqrt(log(1 + 0.345^2) * log(1 + 0.566^2)),
      log(1 + 0.566^2)
    )
    # Table 2: IIV CL_M3G = 34.5% CV (RSE 38), IIV V_M3G = 56.6% CV
    # (RSE 36); correlation Corr(CL_M3G, V_M3G) = 0.751 (RSE 20).
    # Off-diagonal covariance = rho * sqrt(var_cl * var_vc).

    # ============================================================
    # Residual variability -- paper Table 2 'Residual variability (%)'.
    # The paper modeled log-transformed serum concentrations with an
    # exponential residual error model; the reported %CV values are
    # interpreted as proportional residual SDs in linear space
    # (propSd = value / 100). Urine residual errors (62.1% and 66.5%
    # for morphine and M3G respectively) are documented in the
    # vignette source-trace but not declared here because this model
    # file does not emit urine-concentration observations; urinary
    # recovery is exposed only as the cumulative-amount compartments
    # urine_morphine and urine_m3g.
    # ============================================================
    propSd <- 0.272
    label("Proportional residual SD for morphine in serum (fraction)")
    # Table 2: sigma^2_M = 27.2% (RSE 6) for proportional residual error
    # of morphine in serum.

    propSd_m3g <- 0.384
    label("Proportional residual SD for M3G in serum (fraction)")
    # Table 2: sigma^2_M3G = 38.4% (RSE 10) for proportional residual
    # error of M3G in serum.
  })

  model({
    # ------------------------------------------------------------
    # Reference covariate values (paper Methods 'Covariate analysis').
    # ------------------------------------------------------------
    wt_ref   <- 70   # reference body weight (kg) for allometric scaling
    nasf_ref <- 4    # NASF cutoff below which the covariate effect is zero

    # ------------------------------------------------------------
    # Allometric scaling factors -- applied multiplicatively to all
    # CL/Q and V parameters per paper Methods Eqs. 1-2.
    # ------------------------------------------------------------
    allo_cl_factor <- (WT / wt_ref)^e_wt_cl_q
    allo_v_factor  <- (WT / wt_ref)^e_wt_vc_vp

    # ------------------------------------------------------------
    # NASF effect on M3G clearance. Eq. 5 defines a linear effect on
    # ln(NASF / nasf_ref) only when NASF >= nasf_ref; for NASF < 4 the
    # contribution is identically zero (paper Methods 'Covariate
    # analysis'). The clamp `nasf_safe <- max(NASF, nasf_ref)` avoids
    # evaluating log() of values < 4 (including log(0) for healthy
    # subjects with NASF = 0) and reduces to ln(1) = 0 below the cutoff,
    # so the gating multiplier (NASF >= nasf_ref) is also zero in that
    # regime; the two factors are redundant but each makes the cutoff
    # behavior explicit.
    # ------------------------------------------------------------
    nasf_safe   <- ifelse(NASF >= nasf_ref, NASF, nasf_ref)
    nasf_factor <- 1 + e_nasf_cl_m3g * log(nasf_safe / nasf_ref) *
                       (NASF >= nasf_ref)

    # ------------------------------------------------------------
    # Individual morphine PK parameters. CL_M_NR carries IIV; CL_M_R
    # does not (Table 2). All scale by weight via allometric factors.
    # ------------------------------------------------------------
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * allo_cl_factor
    cl_renal  <- exp(lcl_renal)                  * allo_cl_factor
    vc  <- exp(lvc + etalvc) * allo_v_factor
    vp  <- exp(lvp)          * allo_v_factor
    q   <- exp(lq)           * allo_cl_factor
    vp2 <- exp(lvp2)         * allo_v_factor
    q2  <- exp(lq2)          * allo_cl_factor

    # Liver transit rate constant (no IIV; no allometric scaling -- it
    # is a first-order rate constant rather than a CL or V).
    ktrans <- exp(lktrans)

    # ------------------------------------------------------------
    # Individual M3G PK parameters. CL_M3G picks up the NASF covariate
    # factor in addition to allometric scaling. The covariance block
    # captures the high (rho = 0.751) IIV correlation between CL_M3G
    # and V_M3G.
    # ------------------------------------------------------------
    cl_m3g <- exp(lcl_m3g + etalcl_m3g) * allo_cl_factor * nasf_factor
    vc_m3g <- exp(lvc_m3g + etalvc_m3g) * allo_v_factor

    # ------------------------------------------------------------
    # Micro-constants. Working in nanomoles of morphine free base for
    # every compartment so that the mole-for-mole conversion of
    # morphine to M3G across the liver transit compartment is
    # mass-balanced without an explicit stoichiometric factor.
    # ------------------------------------------------------------
    kel_m_nonren <- cl_nonren / vc   # rate from morphine central -> liver transit
    kel_m_renal  <- cl_renal  / vc   # rate from morphine central -> urine_morphine
    k12          <- q  / vc          # rate from morphine central -> peripheral1
    k21          <- q  / vp          # rate from peripheral1 -> morphine central
    k13          <- q2 / vc          # rate from morphine central -> peripheral2
    k31          <- q2 / vp2         # rate from peripheral2 -> morphine central
    kel_m3g      <- cl_m3g / vc_m3g  # rate from M3G central -> urine_m3g

    # ------------------------------------------------------------
    # ODE system. Morphine has a three-compartment disposition with
    # parallel renal and non-renal clearances. The non-renal flux feeds
    # a single liver transit compartment (paper Results: 'A
    # three-compartment model with an additional transit compartment
    # for metabolite conversion'); transit then drains at rate k_trans
    # into the M3G central compartment as mole-for-mole M3G formation
    # (paper Results: 'the fraction metabolized to M3G was assumed to
    # be unity'). Urine compartments accumulate the cumulative renal
    # morphine and total M3G elimination respectively.
    # ------------------------------------------------------------
    d/dt(central) <-
      -kel_m_nonren * central -
       kel_m_renal  * central -
       k12 * central + k21 * peripheral1 -
       k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    d/dt(transit1) <- kel_m_nonren * central - ktrans * transit1

    d/dt(central_m3g) <- ktrans * transit1 - kel_m3g * central_m3g

    d/dt(urine_morphine) <- kel_m_renal * central
    d/dt(urine_m3g)      <- kel_m3g * central_m3g

    # ------------------------------------------------------------
    # Observations. Serum concentrations are reported in nmol/L
    # (equivalent to nM). Cumulative urine amounts (nmol) are
    # available as the urine_morphine and urine_m3g compartments
    # but are not declared as observation variables because urine
    # concentration depends on the per-interval urine volume, which
    # is not a model parameter.
    # ------------------------------------------------------------
    Cc     <- central     / vc
    Cc_m3g <- central_m3g / vc_m3g

    Cc     ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
  })
}
