Zingmark_2003_clomethiazole <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for",
    "clomethiazole (Zingmark 2003) in 774 adult acute-stroke",
    "patients dosed with a three-phase IV infusion of",
    "clomethiazole edisilate over 24 h (6 mg/kg over 0.25 h",
    "then 31 mg/kg over 0.25-8 h then 31 mg/kg over 8-24 h,",
    "total 68 mg/kg edisilate). The structural model is",
    "parameterized in CL/V1/Q/V2 with body weight as a linear",
    "covariate on V1 and V2 and a piecewise-linear covariate",
    "on CL (linear up to WT50 = 100 kg, constant above) plus",
    "a multiplicative effect of concomitant liver-enzyme-",
    "inducing drugs (carbamazepine, phenytoin, rifampicin) on",
    "CL. IIV uses parameter-specific etas combined with a",
    "shared eta common to all four PK parameters (paper text:",
    "attributed to clomethiazole adsorption to the infusion",
    "tubing) -- the joint structure induces a single pairwise",
    "correlation among the structural parameters. The paper",
    "also reports a proportional-odds sedation-score PD model",
    "with a sensitive/non-sensitive mixture component; that",
    "PD layer is not encoded here -- it requires a NONMEM",
    "MIXNUM-style mixture construct that is not naturally",
    "expressed in nlmixr2 / rxode2 model files, and the NIH",
    "stroke-scale covariate is not yet in the canonical",
    "covariate register. See the validation vignette's",
    "Assumptions and deviations section."
  )
  reference <- paste(
    "Zingmark P-H, Ekblom M, Odergren T, Ashwood T, Lyden P,",
    "Karlsson MO, Jonsson EN. (2003). Population",
    "pharmacokinetics of clomethiazole and its effect on the",
    "natural course of sedation in acute stroke patients.",
    "British Journal of Clinical Pharmacology 56(2):173-183.",
    "doi:10.1046/j.0306-5251.2003.01850.x.",
    sep = " "
  )
  vignette <- "Zingmark_2003_clomethiazole"
  units    <- list(
    time          = "hour",
    dosing        = "mg (clomethiazole free base)",
    concentration = "umol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline weight. Acts linearly on V1 and V2",
        "(1.3%/kg and 1.4%/kg respectively, centred at 75 kg) and",
        "piecewise-linearly on CL (0.9%/kg in non-inducer patients,",
        "1.3%/kg in inducer patients, with a smoothed step at",
        "WT50 = 100 kg above which CL is held constant per Zingmark",
        "2003 Methods 'Concomitant medications' / Table 3). The",
        "piecewise relationship is encoded inside model() via a",
        "Hill-style smooth-step function F = WT^50 / (WT^50 +",
        "WT50^50) following the paper's text (Pharmacokinetics",
        "section, Eq. 13). Cohort: median 75 kg (range 31-157 kg)",
        "per Zingmark 2003 Table 2."
      ),
      source_name        = "WT"
    ),
    CYP3A4_IND = list(
      description        = paste(
        "Concomitant 'liver enzyme inducing' drug indicator",
        "(binary; 1 = at least one of carbamazepine, phenytoin,",
        "or rifampicin coadministered at study entry, 0 = none)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no enzyme inducer)",
      notes              = paste(
        "Zingmark 2003 defines the indicator as patients on any of",
        "carbamazepine, phenytoin, or rifampicin (Pharmacokinetics",
        "section, paragraph describing the piecewise WT-CL",
        "relationship). All three drugs are broad-spectrum CYP /",
        "UGT inducers; clomethiazole is metabolised primarily by",
        "CYP2E1 but the pooled-inducer effect captures induction",
        "across the multiple CYP isoforms that contribute to its",
        "clearance. The canonical CYP3A4_IND name is used here",
        "because (a) all three of the paper's listed inducers are",
        "CYP3A4 inducers, (b) the per-paper convention of pooling",
        "broad-spectrum inducers under a single binary indicator",
        "matches the canonical's documented usage pattern, and (c)",
        "the canonical's notes explicitly anticipate per-paper",
        "documentation of which inducer set is pooled. The effect",
        "form combines an additive +39.9% shift on the population-",
        "typical CL at WT = 75 kg AND a different (steeper) WT-on-",
        "CL slope (1.3%/kg vs 0.9%/kg in non-inducer patients);",
        "both effects are encoded inside model(). Prevalence in",
        "the cohort: 100/1546 patients (~6.5%) per Zingmark 2003",
        "Table 2 'Concomitant medication' row 'Liver enzyme",
        "inducers'."
      ),
      source_name        = "Liver enzyme inducers (Table 2)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 774L,
    n_studies       = 3L,
    age_range       = "19-90 years",
    age_median      = "74 years",
    weight_range    = "31-157 kg",
    weight_median   = "75 kg",
    sex_female_pct  = 50.8,
    race_ethnicity  = c(
      Caucasian = 83.6,
      Black     = 10.9,
      Oriental  =  2.6,
      Hispanic  =  2.2,
      Other     =  0.7
    ),
    disease_state   = paste(
      "Acute stroke within 12 h of onset; three phase III",
      "randomized double-blind placebo-controlled trials:",
      "CLASS-I (acute ischaemic stroke with limb weakness,",
      "higher cortical dysfunction, and visual-field",
      "disturbance), CLASS-H (intracerebral haemorrhage), and",
      "CLASS-T (ischaemic stroke treated with t-PA). NIH",
      "stroke-scale score median 16 (range 1-34)."
    ),
    dose_range      = paste(
      "Three-phase IV infusion of clomethiazole edisilate over",
      "24 h: 6 mg/kg over 0.25 h, then 31 mg/kg from 0.25 to 8",
      "h, then 31 mg/kg from 8 to 24 h (total 68 mg/kg",
      "edisilate equivalent to approximately 42.8 mg/kg",
      "clomethiazole free base using the 1:2 base-per-",
      "edisilate stoichiometry). Target steady-state",
      "concentration approximately 10 umol/L per Methods",
      "'Treatment'. The salt-to-base conversion must be",
      "applied at data-assembly time -- the model parameters",
      "describe free-base clomethiazole in plasma."
    ),
    regions         = "United States and Canada (166 centres)",
    co_medication   = paste(
      "Concomitant medications coded as present/absent at",
      "study entry per Zingmark 2003 Table 2 (liver enzyme",
      "inducers 100/1546; CYP2E1 substrates 590/1546;",
      "dicumarol group 212/1546; heparin group 671/1546;",
      "antiplatelet 825/1546; beta-blockers 608/1546;",
      "thiazides 184/1546; loop diuretics 437/1546; calcium",
      "antagonists 428/1546; ACE inhibitors 508/1546; cardiac",
      "glycosides 424/1546; etc.). Only the liver-enzyme-",
      "inducer pool reached the backward-elimination",
      "significance criterion in the PK covariate model."
    ),
    notes           = paste(
      "Of the 1546 patients enrolled across CLASS-I (n=1200",
      "target), CLASS-H (n=200), and CLASS-T (n=200), 774",
      "were randomised to active clomethiazole and",
      "contributed to the PK analysis (the remainder",
      "contributed placebo arms to the PD sedation analysis",
      "only). Total PK observations 2288, of which 111 were",
      "excluded (concentrations > 150 umol/L deemed",
      "unrealistic from phase I prior; missing data;",
      "sampling-vs-dosing inconsistencies). Three plasma",
      "samples per patient per protocol (at 15 min, between 1",
      "and 2 h, and at infusion end or 24 h) with revised",
      "sampling times in the late CLASS-I phase (between 2-10",
      "h, 10-23 h, and 26-36 h). NONMEM VIb FOCE-I",
      "estimation, with the final model re-fitted under",
      "NONMEM V for cross-validation. Bioanalysis: reversed-",
      "phase LC with UV detection at 255 nm at AstraZeneca",
      "R&D Sodertalje, LLOQ 0.2 umol/L, intra-assay",
      "precision 1.4-5.4%. Demographic counts (sex, race)",
      "are reported across all 1546 patients in Table 2;",
      "the 774-subject PK subset is assumed to follow the",
      "same demographic proportions (the source paper does",
      "not stratify Table 2 by treatment arm). The 50.8%",
      "female and race percentages above are derived from",
      "Table 2 by dividing the per-category counts by 1546",
      "(e.g. female 785/1546 = 50.8%, Caucasian 1293/1546 =",
      "83.6%)."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Structural PK parameters at the typical 75 kg patient not on
    # liver enzyme inducing drugs. Zingmark 2003 Table 3 'Final
    # covariate model' column, page 178.
    # ------------------------------------------------------------
    lcl <- log(52.7); label("Clearance (L/h)")                    # Table 3: CL = 52.7 L/h at WT=75 kg, no inducer (RSE 2.9%)
    lvc <- log(82.5); label("Central volume V1 (L)")              # Table 3: V1 = 82.5 L at WT=75 kg (RSE 5.8%)
    lq  <- log(167);  label("Inter-compartmental clearance (L/h)") # Table 3: Q = 167 L/h (RSE 9.2%)
    lvp <- log(335);  label("Peripheral volume V2 (L)")           # Table 3: V2 = 335 L at WT=75 kg (RSE 9.2%)

    # ------------------------------------------------------------
    # WT covariate effects (proportional change per kg, centred at
    # 75 kg). For CL the slope differs between non-inducer and
    # inducer patients (per the paper text describing the piecewise
    # model). The piecewise CL-WT relationship caps the linear
    # increase at WT50 = 100 kg via a Hill-style smooth-step
    # function (also estimated, value below).
    # ------------------------------------------------------------
    e_wt_cl_noind <- 0.009; label("WT effect on CL in non-inducer patients (fraction per kg)") # Table 3: WT on CL = 0.009 (95% CI 0.006-0.011)
    e_wt_cl_ind   <- 0.013; label("WT effect on CL in inducer patients (fraction per kg)")    # Table 3: WT on CL = 0.013 (95% CI 0.005-0.021)
    e_wt_vc       <- 0.013; label("WT effect on V1 (fraction per kg)")                        # Table 3: WT on V1 = 0.013 (95% CI 0.009-0.017)
    e_wt_vp       <- 0.014; label("WT effect on V2 (fraction per kg)")                        # Table 3: WT on V2 = 0.014 (95% CI 0.010-0.017)

    # ------------------------------------------------------------
    # WT50: the body weight above which the linear WT-CL increase
    # is held constant. Encoded as a positive parameter (NOT log-
    # transformed -- the source value is mid-range so the encoded
    # form matches the paper's reported point estimate directly).
    # ------------------------------------------------------------
    wt50 <- 100; label("Body weight cap for the piecewise WT-CL relationship (kg)") # Table 3: WT50 = 100 kg (95% CI 89.3-110.7)

    # ------------------------------------------------------------
    # Liver-enzyme-inducer effect on the population-typical CL at
    # WT = 75 kg (proportional shift, additive on the (1 + .) scale).
    # ------------------------------------------------------------
    e_cyp3a4_ind_cl <- 0.399; label("Multiplicative inducer effect on CL at WT=75 kg (fraction)") # Table 3: Inducer on CL = 0.399 (95% CI 0.076-0.722)

    # ------------------------------------------------------------
    # Interindividual variability. The Zingmark 2003 model uses
    # parameter-specific etas PLUS a single shared eta added to
    # ALL log-parameters (paper text: attributed to clomethiazole
    # adsorption to the infusion tubing, which would correlate
    # apparent CL/V1/Q/V2 systematically across subjects). That
    # additive structure is mathematically equivalent to a single
    # 4x4 block omega matrix whose marginal variance is the sum
    # var_indiv + var_shared and whose off-diagonal covariance
    # equals var_shared (so the pairwise correlation between any
    # two PK parameters is var_shared / sqrt(var_X * var_Y) and
    # equals 22% for the CL-V1 pair Zingmark 2003 reports).
    # The block form is preferred here because nlmixr2's mu-
    # referencing parser rejects (lcl + etalcl_indiv + eta_shared)
    # ("theta + 2 etas") on a single line. Diagonal entries:
    #   var_CL = log(1 + 0.43^2) = 0.169
    #   var_V1 = log(1 + 0.48^2) = 0.207
    #   var_Q  = log(1 + 0.36^2) = 0.122
    #   var_V2 = log(1 + 0.49^2) = 0.215
    # Off-diagonal (constant since the paper reports a single
    # 22% correlation): rho * sqrt(var_CL * var_V1) = 0.22 *
    # sqrt(0.169 * 0.207) = 0.0412.
    # ------------------------------------------------------------
    etalcl + etalvc + etalq + etalvp ~ c(
      0.169,
      0.0412, 0.207,
      0.0412, 0.0412, 0.122,
      0.0412, 0.0412, 0.0412, 0.215
    )                                                                                                                  # Table 3: marginal CVs 43%/48%/36%/49%, pairwise correlation 22% via shared eta

    # ------------------------------------------------------------
    # Residual error. Paper text: "An additive residual error
    # model on the log-transformed data was sufficient". In nlmixr2
    # linear space this is a proportional error model. Table 3:
    # residual CV = 44% (RSE 12%). Using CV approx propSd for the
    # small-CV approximation typically used in popPK summaries.
    # ------------------------------------------------------------
    propSd <- 0.44; label("Proportional residual error (fraction)") # Table 3: residual CV = 44% (RSE 12%) -- log-additive in paper, proportional in linear space
  })

  model({
    # ------------------------------------------------------------
    # 1. Piecewise WT-CL relationship. Encode the Hill-style smooth
    #    step function F = WT^50 / (WT^50 + WT50^50) per the paper
    #    text. F approx 0 below WT50 and approx 1 above WT50.
    #    The capped effective weight inside the WT-CL slope is then
    #    WT_capped = WT * (1 - F) + WT50 * F, so the WT-on-CL
    #    fractional change is fixed at WT_capped above WT50.
    # ------------------------------------------------------------
    f_wt     <- WT^50 / (WT^50 + wt50^50)
    wt_eff   <- WT * (1 - f_wt) + wt50 * f_wt

    # ------------------------------------------------------------
    # 2. WT-CL slope is itself piecewise in the inducer indicator:
    #    non-inducer patients have slope e_wt_cl_noind (0.9%/kg),
    #    inducer patients have slope e_wt_cl_ind (1.3%/kg).
    # ------------------------------------------------------------
    slope_wt_cl <- e_wt_cl_noind * (1 - CYP3A4_IND) + e_wt_cl_ind * CYP3A4_IND

    # ------------------------------------------------------------
    # 3. Individual PK parameters. Log-normal IIV via parameter-
    #    specific etas plus a single shared eta on all log
    #    parameters (the paper attributes the shared term to
    #    clomethiazole adsorption to the infusion tubing).
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (1 + e_cyp3a4_ind_cl * CYP3A4_IND) *
      (1 + slope_wt_cl * (wt_eff - 75))

    vc <- exp(lvc + etalvc) *
      (1 + e_wt_vc * (WT - 75))

    q  <- exp(lq  + etalq)

    vp <- exp(lvp + etalvp) *
      (1 + e_wt_vp * (WT - 75))

    # ------------------------------------------------------------
    # 4. Two-compartment IV PK. Clomethiazole is dosed by IV
    #    infusion (no first-order absorption), so the dose enters
    #    the central compartment directly. Infusion duration is
    #    encoded on the dose record at the event-table level.
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------
    # 5. Observation. central is the amount of clomethiazole free
    #    base (mg) in the central compartment; vc is L; central /
    #    vc is mg/L. Convert mg/L to umol/L by dividing by the
    #    molecular weight of the free base 161.66 g/mol (0.16166
    #    mg/umol). Clomethiazole MW source: clomethiazole free base
    #    C6H8ClNS = 161.65 g/mol (PubChem CID 10109; ChEBI:46955).
    # ------------------------------------------------------------
    Cc <- (central / vc) / 0.16166
    Cc ~ prop(propSd)
  })
}
