Grzesk_2016_m3M3FBS <- function() {
  description <- paste(
    "Preclinical (rat tail artery, ex vivo).",
    "Sigmoidal Emax concentration-response (CRC) model of perfusion pressure",
    "in the isolated, perfused male Wistar rat tail artery, parameterised",
    "for four mechanistically distinct vasoactive agonists (phenylephrine,",
    "arg-vasopressin, mastoparan-7, Bay K8644) in the absence or presence",
    "of the phospholipase-C activator m-3M3FBS (1e-5 M/L pretreatment).",
    "AGONIST_CODE (1..4) selects the agonist; M3M3FBS_PRESENT (0/1) selects",
    "the control vs +m-3M3FBS (EC50, Emax) sub-pair. CONC_AGONIST_M is the",
    "applied agonist concentration (M/L). The model is static (no time",
    "evolution, no PK), encoding only the per-condition concentration-",
    "response relationship reported in Tables I and II.",
    "OUT-OF-SCOPE caveats (see vignette Errata): the source paper is an",
    "ex-vivo concentration-response study that does NOT report a Hill",
    "coefficient, IIV / between-subject variability in the nlmixr2 sense,",
    "or a residual-error structure. The Hill coefficient is fixed to 1",
    "as the structural minimum, and the residual SD is a placeholder",
    "(propSd fixed at 0.10) so the model is nlmixr2-fit-compatible."
  )

  reference <- paste(
    "Grzesk E, Szadujkis-Szadurska K, Wicinski M, Malinowski B, Sinjab TA,",
    "Tejza B, Pujanek M, Janiszewska E, Kopczynska A, Grzesk G.",
    "Effect of 2,4,6-trimethyl-N-[3-(trifluoromethyl)phenyl]benzenesulfonamide",
    "on calcium influx in three contraction models.",
    "Biomed Rep. 2016;4(1):117-121.",
    "doi:10.3892/br.2015.543."
  )

  vignette <- "Grzesk_2016_m3M3FBS"

  units <- list(
    time          = "(none -- static concentration-response with no time evolution)",
    dosing        = "M (mol/L) applied to the perfusate; not a PK dose",
    concentration = "mmHg (perfusion pressure observation)"
  )

  covariateData <- list(
    AGONIST_CODE = list(
      description = paste(
        "Integer 1..4 selecting which vasoactive agonist occupies the",
        "agonist slot of the sigmoidal Emax CRC at simulation time.",
        "Mapping: 1 = phenylephrine (PHE, alpha1-adrenergic);",
        "2 = arg-vasopressin (AVP, V1 vasopressin);",
        "3 = mastoparan-7 (heterotrimeric G-protein direct activator);",
        "4 = Bay K8644 (L-type voltage-gated calcium channel agonist)."
      ),
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Selector covariate; out-of-range values cause every dispatch",
        "indicator to evaluate to 0 and the sigmoidal Emax expression",
        "to evaluate to 0 (no signal) rather than raise an error.",
        "Registered as the canonical AGONIST_CODE in",
        "inst/references/covariate-columns.md."
      ),
      source_name        = "(Grzesk 2016 Table I row labels: PHE / AVP / mastoparan-7 / Bay K8644)"
    ),
    M3M3FBS_PRESENT = list(
      description = paste(
        "Binary indicator selecting the m-3M3FBS-shifted (EC50, Emax)",
        "sub-pair within each agonist. 1 = the artery was pretreated",
        "with the phospholipase-C activator m-3M3FBS at 1e-5 M/L",
        "prior to the agonist CRC titration; 0 = control arm without",
        "m-3M3FBS pretreatment."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (control artery without m-3M3FBS pretreatment)",
      notes              = paste(
        "Time-fixed per artery / CRC. Registered as the canonical",
        "M3M3FBS_PRESENT in inst/references/covariate-columns.md."
      ),
      source_name        = "(Grzesk 2016 Tables I and II: 'controls' vs '+m-3M3FBS' arms)"
    ),
    CONC_AGONIST_M = list(
      description = paste(
        "Applied concentration (M, mol/L) of the agonist selected by",
        "AGONIST_CODE. Polymorphic in agonist identity: a single column",
        "carries phenylephrine M on PHE records, arg-vasopressin M on",
        "AVP records, etc. CONC_AGONIST_M is the dose-response x-axis;",
        "the model has no PK and no time dynamics."
      ),
      units              = "M (mol/L)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the sigmoidal Emax expression",
        "effect = emax * CONC_AGONIST_M^hill / (ec50^hill + CONC_AGONIST_M^hill).",
        "CONC_AGONIST_M = 0 yields effect = 0 (no agonist contraction).",
        "Registered as the canonical CONC_AGONIST_M in",
        "inst/references/covariate-columns.md."
      ),
      source_name        = "(Grzesk 2016 Table I 'EC50 [M/l]' column reports agonist concentrations in mol/L)"
    )
  )

  population <- list(
    species         = "rat (male Wistar) -- ex vivo isolated, perfused tail artery",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    weight_range    = "250-350 g",
    disease_state   = paste(
      "Healthy male Wistar rats; tail artery dissected, cannulated, and",
      "perfused in an isolated-organ bath at 37 C with oxygenated Krebs",
      "solution at 1 mL/min constant flow. Anaesthesia: intraperitoneal",
      "urethane 120 mg/kg followed by cervical dislocation. The study",
      "protocol was approved by the Local Ethics Committee of the",
      "University of Science and Technology (Krakow, Poland)."
    ),
    dose_range      = paste(
      "Concentration-response curves for four agonists",
      "(phenylephrine, arg-vasopressin, mastoparan-7, Bay K8644),",
      "in the absence and presence of m-3M3FBS 1e-5 M/L pretreatment.",
      "Per-CRC sample sizes (n = number of independent CRCs per",
      "condition): PHE 30 / +m-3M3FBS 16; AVP 25 / +m-3M3FBS 16;",
      "mastoparan-7 16 / +m-3M3FBS 16; Bay K8644 16 / +m-3M3FBS 16",
      "(Grzesk 2016 Table I)."
    ),
    regions         = "Poland (Collegium Medicum, Nicolaus Copernicus University, Bydgoszcz)",
    notes           = paste(
      "Ex-vivo isolated-organ-bath preparation; n in the source is the",
      "number of independent CRCs per condition, not the number of",
      "rats (each artery typically yields one CRC). The CRCs were",
      "summarised by the classical pharmacometric van Rossum method",
      "(Grzesk 2016 Methods, 'Data analysis and statistical procedures';",
      "cited refs 21 and 22). The paper does NOT report a Hill",
      "coefficient, BSV / IIV, or residual-error structure -- statistical",
      "comparisons are by Shapiro-Wilk + ANOVA + Newman-Keuls.",
      "Operator-authorised CRC extraction (sidecar response 001 to task",
      "frompeople-731): the Hill coefficient is fabricated at 1 and the",
      "residual SD is a placeholder; see vignette Errata."
    )
  )

  ini({
    # ====================================================================
    # Hill coefficient. Grzesk 2016 does NOT report a Hill coefficient
    # ('classical pharmacometric van Rossum method' summarises CRCs by
    # EC50 / Emax / pD2 without estimating a slope parameter). Fabricated
    # at 1 (no cooperativity, standard simple Emax CRC) and fixed; see
    # vignette Errata for the rationale and a sensitivity discussion.
    # ====================================================================
    lhill <- fixed(log(1))
    label("Hill coefficient (fabricated at 1, fixed; not reported in source)")

    # ====================================================================
    # Per-agonist, per-arm maximal perfusion pressure Emax (mmHg).
    # Source: Grzesk 2016 Table II 'Extracellular calcium phase 2' column
    # (the full-titration maximum-response perfusion pressure with the
    # extracellular calcium pool replete). Each value is the mean of n
    # independent CRCs per condition (n in Table II).
    # ====================================================================
    lemax_phe <- log(94.2)
    label("Phenylephrine Emax, control (mmHg)")
    # Grzesk 2016 Table II row PHE, Phase 2 column: 94.2 +/- 7.9 mmHg (n=30).

    lemax_phe_act <- log(112.2)
    label("Phenylephrine Emax, +m-3M3FBS (mmHg)")
    # Grzesk 2016 Table II row PHE+m-3M3FBS, Phase 2 column: 112.2 +/- 7.1 mmHg (n=16).

    lemax_avp <- log(103.4)
    label("Arg-vasopressin Emax, control (mmHg)")
    # Grzesk 2016 Table II row AVP, Phase 2 column: 103.4 +/- 5.9 mmHg (n=32).

    lemax_avp_act <- log(118.2)
    label("Arg-vasopressin Emax, +m-3M3FBS (mmHg)")
    # Grzesk 2016 Table II row AVP+m-3M3FBS, Phase 2 column: 118.2 +/- 7.5 mmHg (n=16).

    lemax_mas <- log(27.2)
    label("Mastoparan-7 Emax, control (mmHg)")
    # Grzesk 2016 Table II row mastoparan-7, Phase 2 column: 27.2 +/- 5.7 mmHg (n=16).

    lemax_mas_act <- log(42.0)
    label("Mastoparan-7 Emax, +m-3M3FBS (mmHg)")
    # Grzesk 2016 Table II row mastoparan-7+m-3M3FBS, Phase 2 column: 42.0 +/- 6.0 mmHg (n=16).

    lemax_bay <- log(75.2)
    label("Bay K8644 Emax, control (mmHg)")
    # Grzesk 2016 Table II row Bay K8644, Phase 2 column: 75.2 +/- 6.2 mmHg (n=16).

    lemax_bay_act <- log(75.3)
    label("Bay K8644 Emax, +m-3M3FBS (mmHg)")
    # Grzesk 2016 Table II row Bay K8644+m-3M3FBS, Phase 2 column: 75.3 +/- 4.1 mmHg (n=16).

    # ====================================================================
    # Per-agonist, per-arm EC50 (M/L).
    # Source: Grzesk 2016 Table I 'EC50 [M/l]' column. Each value is the
    # mean of n independent CRCs per condition (n in Table I).
    # ====================================================================
    lec50_phe <- log(7.50e-8)
    label("Phenylephrine EC50, control (M/L)")
    # Grzesk 2016 Table I row PHE: 7.50 +/- 0.98 x 10^-8 M/L (n=30); pD2 = 7.12.

    lec50_phe_act <- log(6.45e-8)
    label("Phenylephrine EC50, +m-3M3FBS (M/L)")
    # Grzesk 2016 Table I row PHE+m-3M3FBS: 6.45 +/- 2.10 x 10^-8 M/L (n=16);
    # pD2 = 7.19; relative potency 1.163 vs control (P = 0.0182).

    lec50_avp <- log(1.84e-8)
    label("Arg-vasopressin EC50, control (M/L)")
    # Grzesk 2016 Table I row AVP: 1.84 +/- 0.62 x 10^-8 M/L (n=25); pD2 = 7.74.

    lec50_avp_act <- log(1.42e-8)
    label("Arg-vasopressin EC50, +m-3M3FBS (M/L)")
    # Grzesk 2016 Table I row AVP+m-3M3FBS: 1.42 +/- 0.45 x 10^-8 M/L (n=16);
    # pD2 = 7.85; relative potency 1.296 vs control (P = 0.0071).

    lec50_mas <- log(4.48e-8)
    label("Mastoparan-7 EC50, control (M/L)")
    # Grzesk 2016 Table I row mastoparan-7: 4.48 +/- 2.36 x 10^-8 M/L (n=16); pD2 = 7.34.

    lec50_mas_act <- log(2.55e-8)
    label("Mastoparan-7 EC50, +m-3M3FBS (M/L)")
    # Grzesk 2016 Table I row mastoparan-7+m-3M3FBS: 2.55 +/- 1.52 x 10^-8 M/L (n=16);
    # pD2 = 7.59; relative potency 1.757 vs control (P = 0.0112).

    lec50_bay <- log(1.96e-6)
    label("Bay K8644 EC50, control (M/L)")
    # Grzesk 2016 Table I row Bay K8644: 1.96 +/- 0.26 x 10^-6 M/L (n=16); pD2 = 5.71.

    lec50_bay_act <- log(2.05e-6)
    label("Bay K8644 EC50, +m-3M3FBS (M/L)")
    # Grzesk 2016 Table I row Bay K8644+m-3M3FBS: 2.05 +/- 0.22 x 10^-6 M/L (n=16);
    # pD2 = 5.69; relative potency 0.956 vs control (P = 0.1824, not significant).

    # ====================================================================
    # Residual error. Grzesk 2016 does NOT report a residual-error
    # structure -- the CRCs are summarised by aggregate EC50 / Emax / pD2
    # with between-CRC SDs of those summary statistics (not within-CRC
    # residual error). Placeholder proportional SD of 0.10 (10%) fixed,
    # so the model is nlmixr2-fit-compatible. See vignette Errata.
    # ====================================================================
    propSd <- fixed(0.10)
    label("Placeholder proportional residual SD on perfusion pressure (unreported in source)")
  })

  model({
    # ====================================================================
    # Back-transform from log scale.
    # ====================================================================
    hill <- exp(lhill)

    emax_phe_arm <- exp(lemax_phe)
    emax_phe_act_arm <- exp(lemax_phe_act)
    emax_avp_arm <- exp(lemax_avp)
    emax_avp_act_arm <- exp(lemax_avp_act)
    emax_mas_arm <- exp(lemax_mas)
    emax_mas_act_arm <- exp(lemax_mas_act)
    emax_bay_arm <- exp(lemax_bay)
    emax_bay_act_arm <- exp(lemax_bay_act)

    ec50_phe_arm <- exp(lec50_phe)
    ec50_phe_act_arm <- exp(lec50_phe_act)
    ec50_avp_arm <- exp(lec50_avp)
    ec50_avp_act_arm <- exp(lec50_avp_act)
    ec50_mas_arm <- exp(lec50_mas)
    ec50_mas_act_arm <- exp(lec50_mas_act)
    ec50_bay_arm <- exp(lec50_bay)
    ec50_bay_act_arm <- exp(lec50_bay_act)

    # ====================================================================
    # Within each agonist, select the (Emax, EC50) sub-pair via
    # M3M3FBS_PRESENT (0 = control sub-pair, 1 = +m-3M3FBS sub-pair).
    # ====================================================================
    emax_phe_sel <- emax_phe_arm * (1 - M3M3FBS_PRESENT) + emax_phe_act_arm * M3M3FBS_PRESENT
    emax_avp_sel <- emax_avp_arm * (1 - M3M3FBS_PRESENT) + emax_avp_act_arm * M3M3FBS_PRESENT
    emax_mas_sel <- emax_mas_arm * (1 - M3M3FBS_PRESENT) + emax_mas_act_arm * M3M3FBS_PRESENT
    emax_bay_sel <- emax_bay_arm * (1 - M3M3FBS_PRESENT) + emax_bay_act_arm * M3M3FBS_PRESENT

    ec50_phe_sel <- ec50_phe_arm * (1 - M3M3FBS_PRESENT) + ec50_phe_act_arm * M3M3FBS_PRESENT
    ec50_avp_sel <- ec50_avp_arm * (1 - M3M3FBS_PRESENT) + ec50_avp_act_arm * M3M3FBS_PRESENT
    ec50_mas_sel <- ec50_mas_arm * (1 - M3M3FBS_PRESENT) + ec50_mas_act_arm * M3M3FBS_PRESENT
    ec50_bay_sel <- ec50_bay_arm * (1 - M3M3FBS_PRESENT) + ec50_bay_act_arm * M3M3FBS_PRESENT

    # ====================================================================
    # Dispatch on AGONIST_CODE (1 = PHE, 2 = AVP, 3 = mastoparan-7, 4 = Bay K8644).
    # Out-of-range AGONIST_CODE leaves all four indicators at 0, which
    # zeroes the sigmoidal expression (effect = 0) rather than raising.
    # ====================================================================
    is_phe <- (AGONIST_CODE == 1)
    is_avp <- (AGONIST_CODE == 2)
    is_mas <- (AGONIST_CODE == 3)
    is_bay <- (AGONIST_CODE == 4)

    emax <- emax_phe_sel * is_phe + emax_avp_sel * is_avp + emax_mas_sel * is_mas + emax_bay_sel * is_bay
    ec50 <- ec50_phe_sel * is_phe + ec50_avp_sel * is_avp + ec50_mas_sel * is_mas + ec50_bay_sel * is_bay

    # ====================================================================
    # Sigmoidal Emax concentration-response (Grzesk 2016 Methods,
    # 'classical pharmacometric van Rossum method'). The paper does not
    # write out the parametric form; the canonical sigmoidal Emax is used
    # with hill fixed at 1 (see vignette Errata for the choice).
    # ====================================================================
    effect <- emax * CONC_AGONIST_M^hill / (ec50^hill + CONC_AGONIST_M^hill)

    # ====================================================================
    # Observation model. Placeholder proportional residual (paper does
    # not report a residual-error structure); 10% on effect.
    # ====================================================================
    effect ~ prop(propSd)
  })
}
