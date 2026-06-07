# Population PK model for iloperidone and its two major plasma metabolites
# P-88 (M1) and P-95 (M2) in Chinese schizophrenia patients, quantifying the
# influence of the CYP2D6*10 (rs1065852) polymorphism on the metabolite
# formation rate constants (Pei 2016 Acta Pharmacol Sin 37(11):1499-1508;
# doi:10.1038/aps.2016.96).

Pei_2016_iloperidone <- function() {
  description <- paste(
    "Population PK model for iloperidone and its two major plasma metabolites",
    "P-88 (M1, contributes to the therapeutic profile via D2 / 5-HT2A binding",
    "affinity comparable to the parent) and P-95 (M2, CYP2D6-mediated",
    "hydroxylation metabolite, pharmacologically less active) in 70 Chinese",
    "patients with schizophrenia or schizoaffective disorder receiving oral",
    "iloperidone 12-24 mg/day twice daily (Pei 2016).",
    "One-compartment first-order absorption (Ka FIXED at 2.26 1/h, estimated",
    "in a separate forward analysis of healthy-volunteer concentration-time",
    "data digitised from Pei 2016 ref [23] and fixed for the patient model to",
    "stabilize absorption identification under sparse sampling) and",
    "first-order parallel-pathway elimination of iloperidone via three rate",
    "constants: K20 (other elimination pathways), K23 (formation of M1), and",
    "K24 (formation of M2). Each metabolite then occupies its own",
    "one-compartment model with FIXED apparent volume V3 = V4 = 10 L (the",
    "fractions of iloperidone converted to each metabolite are not",
    "identifiable from the cohort because no co-administered tracer was",
    "available, so K23, K30 and K24, K40 are estimated against the FIXED",
    "apparent metabolite volume per Methods) and first-order elimination",
    "(K30 for M1, K40 for M2). Inter-occasion variability was retained on",
    "K20 in the published final model but is NOT carried as a separate eta",
    "in this nlmixr2lib extraction (see Errata in the validation vignette).",
    "Mass units (mg) rather than molar units are used per the source's",
    "convention because the molecular weights of iloperidone (427.3 g/mol),",
    "M1 (429.4 g/mol), and M2 (429.2 g/mol) are within 0.5%.",
    "CYP2D6*10 (rs1065852) polymorphism affects both metabolite formation",
    "rate constants: T/T homozygotes have K23 1.34-fold the C/C + C/T",
    "pooled reference; C/T heterozygotes and T/T homozygotes have K24",
    "reduced to 0.693 and 0.492 of the C/C wild-type reference",
    "respectively.",
    sep = " "
  )
  reference <- paste(
    "Pei Q, Huang L, Huang J, Gu J-K, Kuang Y, Zuo X-C, Ding J-J, Tan H-Y,",
    "Guo C-X, Liu S-K, Yang G-P (2016). Influences of CYP2D6*10 polymorphisms",
    "on the pharmacokinetics of iloperidone and its metabolites in Chinese",
    "patients with schizophrenia: a population pharmacokinetic analysis.",
    "Acta Pharmacol Sin 37(11):1499-1508. doi:10.1038/aps.2016.96.",
    sep = " "
  )
  vignette <- "Pei_2016_iloperidone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_etas <- c("etalkel_form_p88", "etalkel_form_p95")

  covariateData <- list(
    CYP2D6_STAR10_HET = list(
      description        = "CYP2D6*10 (rs1065852) heterozygote indicator (C/T genotype)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type at rs1065852; paired with CYP2D6_STAR10_HOM)",
      notes              = paste(
        "1 = subject is CYP2D6*10 C/T heterozygote at rs1065852, 0 otherwise.",
        "Paired with CYP2D6_STAR10_HOM to encode the three-level C/C / C/T / T/T",
        "genotype using two binary indicators on the CYP3A5_STAR1_HET /",
        "CYP3A5_STAR1_HOM pattern with C/C wild-type as the implicit reference",
        "(both indicators = 0). Used by the K24 (M2 formation rate) covariate",
        "in the final model: C/T subjects multiply the typical K24 by 0.693",
        "relative to the C/C reference (Pei 2016 Results 'Population",
        "pharmacokinetic (PPK) models': `(K24)_i = theta * 0.00649 * exp(eta_i)`",
        "with theta = 0.693 for C/T). K23 (M1 formation) does NOT depend on",
        "CYP2D6_STAR10_HET in the final model because the paper pooled C/C and",
        "C/T together as a single reference group for K23. Time-fixed per",
        "subject (germline genotype). Cohort distribution per Pei 2016 Table 1:",
        "C/C 11/70 (15.7%), C/T 42/70 (60.0%), T/T 17/70 (24.3%)."
      ),
      source_name        = "CYP2D6*10 genotype (rs1065852); the C/T = 2 stratum in the paper's GENOTYPE = 1/2/3 encoding maps to CYP2D6_STAR10_HET = 1"
    ),
    CYP2D6_STAR10_HOM = list(
      description        = "CYP2D6*10 (rs1065852) homozygous-mutant indicator (T/T genotype)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type at rs1065852; paired with CYP2D6_STAR10_HET)",
      notes              = paste(
        "1 = subject is CYP2D6*10 T/T homozygote at rs1065852, 0 otherwise.",
        "Paired with CYP2D6_STAR10_HET; both indicators = 0 indicates the C/C",
        "wild-type reference. Used by both metabolite-formation rate constants",
        "in the final model: K23 (M1 formation) is multiplied by 1.34 in T/T",
        "subjects relative to the C/C + C/T pooled reference (Pei 2016",
        "Results: `(K23)_i = theta * 0.00451 * exp(eta_i)` with theta = 1.34",
        "for T/T); K24 (M2 formation) is multiplied by 0.492 in T/T subjects",
        "relative to the C/C reference (`(K24)_i = theta * 0.00649 * exp(eta_i)`",
        "with theta = 0.492 for T/T). Time-fixed per subject (germline",
        "genotype). Cohort distribution per Pei 2016 Table 1: C/C 11/70",
        "(15.7%), C/T 42/70 (60.0%), T/T 17/70 (24.3%)."
      ),
      source_name        = "CYP2D6*10 genotype (rs1065852); the T/T = 3 stratum in the paper's GENOTYPE = 1/2/3 encoding maps to CYP2D6_STAR10_HOM = 1"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 70L,
    n_observations      = 804L,
    n_studies           = 1L,
    age_range           = "18-65 years (mean 34 +/- 12; Pei 2016 Table 1)",
    weight_range        = "mean 62.2 +/- 10.5 kg (Pei 2016 Table 1)",
    sex_female_pct      = 60.0,
    race_ethnicity      = c(Asian = 100),
    disease_state       = "Schizophrenia or schizoaffective disorder (DSM-IV); acute psychotic episode at enrolment.",
    dose_range          = "Oral iloperidone 12-24 mg/day twice daily (12, 16, 20, or 24 mg/day) after individual dose titration on days 5-14 and fixed-dose maintenance on days 15-28; titration regimen 2 mg/day on day 1, 4 mg/day on day 2, 8 mg/day on day 3, 12 mg/day on day 4.",
    regions             = "China (multi-center: 13 hospitals)",
    cyp2d6_distribution = "CYP2D6*10 (rs1065852) genotype: C/C 11/70 (15.7%), C/T 42/70 (60.0%), T/T 17/70 (24.3%) per Pei 2016 Table 1. Allele frequencies for the *10 variant in Chinese populations are 48-70% per Introduction.",
    notes               = paste(
      "Sparse sampling. Four plasma samples per patient: C15-0 (predose on",
      "day 15, the first day of fixed-dose maintenance), C28-0 (predose on",
      "day 28), C28-4 (4 h post-AM-dose on day 28), and C28-12 (12 h",
      "post-AM-dose on day 28). 804 total concentration measurements across",
      "the three analytes: 266 iloperidone, 268 M1, 270 M2. BQL fractions",
      "were 1.48% (4/270) for iloperidone, 0.74% (2/270) for M1, and 0%",
      "(0/270) for M2; BQL records were excluded from the analysis. Plasma",
      "concentrations of iloperidone, M1, and M2 were simultaneously",
      "measured by validated HPLC-MS (linear ranges 1-100 ng/mL for parent",
      "and 3-120 ng/mL for metabolites)."
    )
  )

  ini({
    # ====================================================================
    # Structural rate constants and volumes - Pei 2016 Table 3 (final model)
    # ====================================================================
    # Ka FIXED at 2.26 1/h per Methods. Pei 2016 estimated Ka in a
    # separate forward analysis of healthy-volunteer iloperidone
    # concentration-time data digitised from ref [23] (Getdata Graph
    # Digitizer + two-compartment fit, which produced Ka = 1.68 1/h) and
    # subsequently fixed Ka at 2.26 1/h for the patient model to stabilize
    # absorption identification under sparse sampling. The RSE column in
    # Table 3 is reported as 'FIX'.
    lka              <- fixed(log(2.26))
    label("Iloperidone absorption rate constant Ka (1/h); FIXED per Pei 2016 Methods and Table 3")  # Pei 2016 Table 3: Ka = 2.26 (FIX)

    # K20 = iloperidone elimination via 'other' (non-M1, non-M2) pathways.
    # Inter-occasion variability (27.3% CV) is reported in Table 3 in
    # addition to the BSV (19.5% CV) but is NOT carried as a separate
    # IOV eta in this extraction; see vignette Errata for the rationale.
    lkel             <- log(0.067)
    label("Iloperidone non-metabolic elimination rate constant K20 (1/h)")  # Pei 2016 Table 3: K20 = 0.067 (RSE 13.3%; bootstrap 95% CI 0.051-0.082)

    # K23 = iloperidone -> M1 (P-88) formation rate constant. In the
    # final model only the CYP2D6*10 T/T (HOM) stratum has a distinct
    # typical value; C/C and C/T are pooled as the reference. The eta
    # name `etalkel_form_p88` is declared in `paper_specific_etas`
    # because `lkel_form_p88` is a paper-mechanistic parent-to-metabolite
    # formation rate rather than a canonical `lkel`-suffixed bare PK
    # parameter (the canonical `lkel_p88` is reserved for the metabolite
    # elimination rate, K30; see below).
    lkel_form_p88    <- log(0.00451)
    label("Iloperidone -> M1 (P-88) formation rate constant K23 in C/C or C/T reference subjects (1/h)")  # Pei 2016 Table 3: K23 = 0.00451 (RSE 27.7%; bootstrap 95% CI 0.00276-0.00627)

    # K24 = iloperidone -> M2 (P-95) formation rate constant. Both the
    # CYP2D6*10 C/T (HET) and T/T (HOM) strata have distinct typical
    # values relative to the C/C reference.
    lkel_form_p95    <- log(0.00649)
    label("Iloperidone -> M2 (P-95) formation rate constant K24 in C/C reference subjects (1/h)")  # Pei 2016 Table 3: K24 = 0.00649 (RSE 31.1%; bootstrap 95% CI 0.00292-0.0113)

    # K30 = M1 elimination rate constant; K40 = M2 elimination rate
    # constant. No covariates retained on K30 or K40 in the final model.
    # No IIV reported on K30 or K40 in Table 3.
    lkel_p88         <- log(0.160)
    label("M1 (P-88) elimination rate constant K30 (1/h)")  # Pei 2016 Table 3: K30 = 0.160 (RSE 30.6%; bootstrap 95% CI 0.101-0.223)
    lkel_p95         <- log(0.106)
    label("M2 (P-95) elimination rate constant K40 (1/h)")  # Pei 2016 Table 3: K40 = 0.106 (RSE 29.6%; bootstrap 95% CI 0.057-0.173)

    # Apparent volume of distribution of iloperidone V2 (estimated).
    lvc              <- log(679)
    label("Iloperidone apparent volume of distribution V2 (L)")  # Pei 2016 Table 3: V2 = 679 L (RSE 13.4%; bootstrap 95% CI 495-979)

    # Apparent volumes of M1 and M2 distribution V3, V4 are FIXED at 10 L
    # per Methods because the fractions of iloperidone converted to each
    # metabolite are not identifiable from the cohort (no co-administered
    # tracer); K23, K30 and K24, K40 are estimated against this fixed
    # apparent volume. The RSE column in Table 3 is reported as 'FIX'.
    lvc_p88          <- fixed(log(10))
    label("M1 (P-88) apparent volume of distribution V3 (L); FIXED per Pei 2016 Methods")  # Pei 2016 Table 3: V3 = 10 (FIX)
    lvc_p95          <- fixed(log(10))
    label("M2 (P-95) apparent volume of distribution V4 (L); FIXED per Pei 2016 Methods")  # Pei 2016 Table 3: V4 = 10 (FIX)

    # ====================================================================
    # CYP2D6*10 covariate effects on the metabolite formation rate
    # constants. The paper expresses the effect as a multiplicative
    # theta on the typical rate, e.g., (K23)_i = theta_TT * 0.00451 *
    # exp(eta) with theta_TT = 1.34. Equivalent log-scale shift on the
    # typical-value parameter: log(theta) added to log(K23_typ) before
    # exponentiation. The covariate-effect names end in the metabolite
    # suffix `_p88` / `_p95` so they match the canonical
    # `e_<cov>_<param>_<metab>` three-token covariate-effect form.
    # ====================================================================
    e_cyp2d6_star10_hom_kel_form_p88 <- log(1.34)
    label("Log fold-change in K23 (iloperidone -> M1 formation) for CYP2D6*10 T/T (HOM) relative to C/C + C/T pooled reference (unitless)")  # Pei 2016 Table 3: CYP2D6*10 T/T on K23 = 1.34 (RSE 12.4%; bootstrap 95% CI 1.15-1.55)
    e_cyp2d6_star10_het_kel_form_p95 <- log(0.693)
    label("Log fold-change in K24 (iloperidone -> M2 formation) for CYP2D6*10 C/T (HET) relative to C/C wild-type reference (unitless)")  # Pei 2016 Table 3: CYP2D6*10 C/T on K24 = 0.693 (RSE 18.9%; bootstrap 95% CI 0.518-0.919)
    e_cyp2d6_star10_hom_kel_form_p95 <- log(0.492)
    label("Log fold-change in K24 (iloperidone -> M2 formation) for CYP2D6*10 T/T (HOM) relative to C/C wild-type reference (unitless)")  # Pei 2016 Table 3: CYP2D6*10 T/T on K24 = 0.492 (RSE 24.2%; bootstrap 95% CI 0.362-0.652)

    # ====================================================================
    # Inter-individual variability (BSV). Pei 2016 Table 3 reports CV%
    # values for the lognormal IIV on K20, V2, K23, and K24. Convert to
    # log-scale variance via omega^2 = log(1 + CV^2). No IIV reported on
    # Ka (fixed), K30, K40, V3, V4 (fixed). The paper also reports an
    # inter-occasion variability of 27.3% CV on K20 (Table 3); IOV is
    # NOT carried as a separate eta in this extraction -- a model that
    # needs occasion-specific simulation should add it via the rxode2
    # IOV mechanism. See vignette Errata.
    # ====================================================================
    etalkel          ~ log(1 + 0.195^2)  # Pei 2016 Table 3: BSV K20 = 19.5% CV (RSE 34%; bootstrap 95% CI 6.7-29.8%)
    etalvc           ~ log(1 + 0.358^2)  # Pei 2016 Table 3: BSV V2 = 35.8% CV (RSE 19%; bootstrap 95% CI 25.6-45.6%)
    etalkel_form_p88 ~ log(1 + 0.196^2)  # Pei 2016 Table 3: BSV K23 = 19.6% CV (RSE 26%; bootstrap 95% CI 8.5-28.1%)
    etalkel_form_p95 ~ log(1 + 0.409^2)  # Pei 2016 Table 3: BSV K24 = 40.9% CV (RSE 13%; bootstrap 95% CI 30.5-48.2%)

    # ====================================================================
    # Residual variability. Pei 2016 Table 3 reports proportional CV%
    # for each of the three analytes (no additive term reported).
    # nlmixr2's propSd carries the proportional residual SD on a
    # fractional scale (SD of multiplicative noise), so the CV% values
    # are entered as fractions directly.
    # ====================================================================
    propSd           <- 0.377
    label("Proportional residual error on iloperidone Cc (fraction)")  # Pei 2016 Table 3: CV1 = 37.7% (RSE 15%; bootstrap 95% CI 32.8-42.6%)
    propSd_p88       <- 0.236
    label("Proportional residual error on M1 (P-88) Cc_p88 (fraction)")  # Pei 2016 Table 3: CV2 = 23.6% (RSE 9%; bootstrap 95% CI 19.1-29.0%)
    propSd_p95       <- 0.221
    label("Proportional residual error on M2 (P-95) Cc_p95 (fraction)")  # Pei 2016 Table 3: CV3 = 22.1% (RSE 7%; bootstrap 95% CI 17.2-26.5%)
  })

  model({
    # ====================================================================
    # Individual rate constants and volumes.
    # ====================================================================
    ka              <- exp(lka)
    kel             <- exp(lkel + etalkel)
    # K23 (iloperidone -> M1 formation): typical-value times exp(log-shift
    # for CYP2D6*10 T/T) times exp(eta). C/C and C/T are pooled in the
    # reference (the HOM shift gates on `CYP2D6_STAR10_HOM` alone).
    kel_form_p88    <- exp(lkel_form_p88 + etalkel_form_p88 +
                           e_cyp2d6_star10_hom_kel_form_p88 * CYP2D6_STAR10_HOM)
    # K24 (iloperidone -> M2 formation): typical-value with separate
    # log-shifts for C/T and T/T relative to the C/C wild-type reference.
    kel_form_p95    <- exp(lkel_form_p95 + etalkel_form_p95 +
                           e_cyp2d6_star10_het_kel_form_p95 * CYP2D6_STAR10_HET +
                           e_cyp2d6_star10_hom_kel_form_p95 * CYP2D6_STAR10_HOM)
    kel_p88         <- exp(lkel_p88)
    kel_p95         <- exp(lkel_p95)
    vc              <- exp(lvc + etalvc)
    vc_p88          <- exp(lvc_p88)
    vc_p95          <- exp(lvc_p95)

    # ====================================================================
    # ODE system: depot -> central -> (M1 / M2 in parallel)
    # Iloperidone in `central` eliminates via three parallel first-order
    # pathways: K20 (`kel`, 'other'), K23 (`kel_form_p88`, to M1), and
    # K24 (`kel_form_p95`, to M2). Each metabolite then eliminates via
    # its own first-order rate constant (K30 = `kel_p88` for M1, K40 =
    # `kel_p95` for M2). Mass units (mg) carry through because the
    # molecular weights of iloperidone (427.3 g/mol), M1 (429.4), and M2
    # (429.2) are within 0.5% (Pei 2016 Methods 'Population
    # pharmacokinetic (PPK) modeling').
    # ====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (kel + kel_form_p88 + kel_form_p95) * central
    d/dt(central_p88) <-  kel_form_p88 * central - kel_p88 * central_p88
    d/dt(central_p95) <-  kel_form_p95 * central - kel_p95 * central_p95

    # Plasma concentrations in ng/mL: dose is administered in mg, volumes
    # are in L; central/vc has units mg/L = ug/mL; multiply by 1000 ->
    # ng/mL, matching Pei 2016 Table 1 / Figure 2 observed C/D ratios
    # (e.g. typical iloperidone C28-0/D ~ 0.7-1.0 ng/mL per mg).
    Cc      <- 1000 * central     / vc
    Cc_p88  <- 1000 * central_p88 / vc_p88
    Cc_p95  <- 1000 * central_p95 / vc_p95

    # Proportional residual error per output (Pei 2016 Table 3 - CV1, CV2, CV3).
    Cc      ~ prop(propSd)
    Cc_p88  ~ prop(propSd_p88)
    Cc_p95  ~ prop(propSd_p95)
  })
}
