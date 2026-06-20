vonHentig_2009_saquinavir <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for oral",
    "ritonavir-boosted saquinavir (1000/100 mg BID) in 136 HIV-1-infected",
    "adults including 13 pregnant women. Apparent oral clearance CL/F is",
    "modulated by two retained covariates: a binary atazanavir-coadministration",
    "indicator (CONMED_ATAZANAVIR; 49 of 136 patients on ATV 300 mg QD) as a",
    "power-of-binary multiplier 0.703^CONMED_ATAZANAVIR (30% CL reduction",
    "when atazanavir is coadministered), and the per-subject ritonavir 12 h",
    "AUC (CONMED_RTV_AUC_12h, cohort median 6.70355 mg*h/L) as a normalised",
    "power form (CONMED_RTV_AUC_12h / 6.70355)^(-0.403). Saquinavir formulation",
    "(Invirase hard gel vs Fortovase soft gel) was tested and not retained.",
    "Inter-individual variability is estimated on CL/F (53.1% CV) and V/F",
    "(54.8% CV); IIV on ka was rejected during model building. Residual error",
    "was reported as an additive-error model but the additive SD value is not",
    "reported anywhere in the paper -- addSd is encoded as fixed(0) and the",
    "vignette Errata documents the omission (von Hentig & Loetsch 2009)."
  )
  reference <- paste(
    "von Hentig N, Loetsch J. Cytochrome P450 3A inhibition by atazanavir and",
    "ritonavir, but not demography or drug formulation, influences saquinavir",
    "population pharmacokinetics in human immunodeficiency virus type",
    "1-infected adults. Antimicrob Agents Chemother. 2009 Aug;53(8):3524-3527.",
    "doi:10.1128/AAC.00025-09."
  )
  vignette <- "vonHentig_2009_saquinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_ATAZANAVIR = list(
      description        = "Concomitant atazanavir coadministration indicator (1 = on ATV 300 mg QD, 0 = no ATV)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atazanavir; nucleosidic reverse-transcriptase inhibitors backbone)",
      notes              = paste(
        "Per-subject time-fixed in von Hentig 2009. 49 of 136 patients",
        "(36%) received atazanavir 300 mg once daily as part of a boosted",
        "double-protease-inhibitor regimen because reverse-transcriptase",
        "inhibitors had had toxic effects or were ineffective due to viral",
        "resistance; the remaining 87 patients received NRTIs. Enters",
        "saquinavir CL/F as a power-of-binary multiplier:",
        "cl = exp(lcl) * e_atazanavir_cl^CONMED_ATAZANAVIR * ...",
        "with e_atazanavir_cl = 0.703 (Table 2 final model) so the on-ATV",
        "arm has CL/F reduced to 70.3% of the no-ATV reference. ATV is a",
        "CYP3A4 inhibitor and the reduction is attributed to CYP3A",
        "competition with saquinavir metabolism."
      ),
      source_name        = "atazanavir"
    ),
    CONMED_RTV_AUC_12h = list(
      description        = "Per-subject ritonavir AUC over the 12 h dosing interval (BID ritonavir 100 mg)",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject ritonavir AUC over the 12 h dosing interval, computed",
        "by the source authors from the observed ritonavir concentration",
        "profile (predose and 1, 2, 4, 6, 9, 12 h) using the",
        "log-trapezoidal rule. Cohort median 6703.55 ng/ml*h =",
        "6.70355 mg*h/L (Table 1 footnote). Enters saquinavir CL/F via a",
        "centred power form: cl = ... * (CONMED_RTV_AUC_12h / 6.70355)^",
        "e_rtv_auc_12h_cl with e_rtv_auc_12h_cl = -0.403 (Table 2 final",
        "model). At the cohort median the power term evaluates to 1 and",
        "CL/F equals the typical 60.4 L/h (without atazanavir)."
      ),
      source_name        = "AUCritonavir"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Median 41.5 years for men (range 20-71), 32.5 years for women",
        "(range 19-64) per Methods. Screened on CL/F and V/F per Table 1;",
        "did not meet the dOFV thresholds (delta-2LL < 6.63 to enter, or",
        "< 10.83 to retain) and was not in the final model."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Median 74 kg for men (range 51-116), 68.5 kg for women",
        "(range 42-100) per Methods. Screened on CL/F and V/F per Table 1;",
        "did not meet retention thresholds and was not in the final model."
      )
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "32 of 136 patients (24%) female per Methods. Screened on CL/F and",
        "V/F per Table 1 as (theta_sex)^sex with sex = 0 male / 1 female;",
        "did not meet retention thresholds and was not in the final model."
      )
    ),
    PREG = list(
      description = "Pregnancy status indicator (1 = pregnant at PK sample)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "13 of 32 women pregnant (mean gestational age 32 weeks + 4 days,",
        "range 24 wk 3 d to 36 wk 5 d) per Methods. Screened on CL/F and",
        "V/F per Table 1; pregnancy met the entry threshold for V/F",
        "(delta-2LL -11.34) but on backward elimination its removal",
        "produced delta-2LL +5.2 which is below the retention threshold",
        "of +10.83, so the covariate was dropped (Discussion)."
      )
    ),
    SAQ_FORMULATION = list(
      description = "Saquinavir formulation (0 = Invirase hard gel, 1 = Fortovase soft gel)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "84 of 136 patients received Invirase (hard gel), 52 received",
        "Fortovase (soft gel). Screened on CL/F and V/F per Table 1; on",
        "backward elimination the V/F effect produced delta-2LL +4.72,",
        "below the +10.83 retention threshold, so formulation was dropped.",
        "Provides confidence that product switching does not substantially",
        "jeopardize saquinavir therapy (Discussion)."
      )
    ),
    CONMED_RTV_AUC_12h_AUCATAZ = list(
      description = "Per-subject atazanavir AUC over the 24 h dosing interval (when CONMED_ATAZANAVIR = 1)",
      units       = "ng/ml*h",
      type        = "continuous",
      notes       = paste(
        "Cohort median 24029.6 ng/ml*h per Table 1 footnote. Tested as a",
        "continuous centred power-form covariate on CL/F and V/F in",
        "addition to (or instead of) the binary CONMED_ATAZANAVIR. Did not",
        "meet the entry threshold (delta-2LL -3.4 on CL/F, -1.26 on V/F",
        "vs the +6.63 entry threshold) and is not in the final model. The",
        "binary CONMED_ATAZANAVIR indicator was retained instead.",
        "Recorded here for provenance only; not a column the simulation",
        "needs."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 136L,
    n_studies       = 1L,
    n_observations  = NA_integer_,
    age_range       = "19-71 years (median 41.5 men, 32.5 women)",
    weight_range    = "42-116 kg (median 74 men, 68.5 women)",
    sex_female_pct  = 23.5,
    disease_state   = "HIV-1 infection; Child-Pugh class B/C and CYP3A4-modulating non-antiretroviral cotherapies excluded; 13 of 32 women pregnant (mean GA 32 wk 4 d)",
    dose_range      = "Saquinavir 1000 mg BID (Invirase hard gel n=84 or Fortovase soft gel n=52) + ritonavir 100 mg BID; 49 of 136 patients additionally received atazanavir 300 mg QD; 87 received nucleosidic reverse-transcriptase inhibitors instead of atazanavir",
    regions         = "Germany (Goethe-University, Frankfurt am Main)",
    notes           = paste(
      "Single 12-h sparse PK profile per patient (predose and 1, 2, 4, 6,",
      "9, 12 h post-dose) taken between the 9th and 2303rd",
      "saquinavir-ritonavir dose (median 61st) at steady state, after an",
      "overnight fast and a standardized 595 kcal breakfast (21% fat)",
      "served after drug administration. Concentrations measured by",
      "HPLC-MS/MS (LLOQ 20 ng/ml, linear to 20000 ng/ml, CV <20% across",
      "the range). NONMEM VI level 2 FOCE with eta-eps interaction. A",
      "two-compartment structural model did not converge so the",
      "one-compartment form was retained (Results paragraph 4)."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- von Hentig 2009 Table 2 full
    # (final) model column.
    # ============================================================
    lcl <- log(60.4)
    label("Apparent oral clearance CL/F at no-atazanavir + median CONMED_RTV_AUC_12h (L/h)")  # Table 2 full: CL/F = 60.4 L/h (95% CI 52.7-69)
    lvc <- log(126)
    label("Apparent volume of distribution V/F (L)")                                          # Table 2 full: V/F = 126 L (95% CI 105-147)
    lka <- log(0.21)
    label("First-order absorption rate ka (1/h)")                                             # Table 2 full: ka = 0.21 /h (95% CI 0.19-0.23)

    # ============================================================
    # Covariate effects on saquinavir CL/F.
    # Final-model CL/F equation (Table 2 full ODE column):
    #   CL = exp(lcl) * theta1_ATV^atazanavir
    #            * (AUC_ritonavir / 6703.55 ng/ml*h)^theta2_RTV
    # theta1_ATV = 0.703  (a power-of-binary multiplier; CL is 70.3%
    # of the no-ATV typical value when atazanavir is coadministered).
    # theta2_RTV = -0.403 (a centred power-form exponent; CL falls as
    # ritonavir AUC rises, consistent with stronger CYP3A inhibition).
    # 6.70355 mg*h/L = 6703.55 ng/ml*h is the cohort median of the
    # individual 12 h ritonavir AUCs per Table 1 footnote / Table 2 row d.
    # ============================================================
    e_atazanavir_cl <- 0.703
    label("Power-of-binary multiplier of atazanavir coadministration on CL/F (unitless; CL = ... * e_atazanavir_cl^CONMED_ATAZANAVIR)")  # Table 2 full: theta1_atazanavir = 0.703 (95% CI 0.58-0.87)
    e_rtv_auc_12h_cl <- -0.403
    label("Power exponent of normalised ritonavir 12 h AUC on CL/F (unitless; centred at 6.70355 mg*h/L)")                               # Table 2 full: theta2_ritonavir = -0.403 (95% CI -0.59 to -0.23)

    # ============================================================
    # Inter-individual variability -- log-normal on CL/F and V/F
    # (NONMEM exponential-eta model). Paper reports IIV as %CV and
    # uses omega^2 = log(1 + CV^2) for the internal log-normal
    # variance. IIV on ka was rejected during model building
    # (delta-2LL only -0.61 vs the +6.63 entry threshold; Results
    # paragraph 4).
    # ============================================================
    etalcl ~ log(1 + 0.531^2)  # Table 2 full: IIV CL/F = 53.1% CV (95% CI 44.8-61.1); omega^2 = log(1 + 0.531^2) = 0.2485
    etalvc ~ log(1 + 0.548^2)  # Table 2 full: IIV V/F = 54.8% CV (95% CI 41.2-65.9); omega^2 = log(1 + 0.548^2) = 0.2627

    # ============================================================
    # Residual error -- additive model per Results paragraph 4
    # ("an additive-error model were found to provide best fits").
    # The additive SD VALUE is NOT reported anywhere in the paper
    # (not in Table 1, not in Table 2, not in prose). Per the
    # general rule for unreported RUV magnitudes the value is
    # encoded as fixed(0); the vignette Errata documents the gap.
    # Downstream simulation users who need a prediction-interval
    # band can override the value when calling rxSolve.
    # ============================================================
    addSd <- fixed(0)
    label("Additive residual SD (mg/L); not reported in paper -- fixed(0) placeholder, see vignette Errata")
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters.
    # ------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) *
      e_atazanavir_cl^CONMED_ATAZANAVIR *
      (CONMED_RTV_AUC_12h / 6.70355)^e_rtv_auc_12h_cl

    # Micro-constant.
    kel <- cl / vc

    # ------------------------------------------------------------
    # ODE system -- one compartment with first-order absorption
    # from a depot. Verbatim from Table 2 column 1 (full model):
    #   dA(0)/dt = F * dose - ka * A(0)
    #   dA(1)/dt = ka * A(0) - CL * A(1) / V
    # The bioavailability term F is the structural anchor F = 1 in
    # NONMEM and is left implicit in the depot ODE (CL/F and V/F
    # are the apparent parameters).
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------
    # Observation and additive residual error. addSd is fixed at 0
    # because the paper does not report the additive magnitude.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
