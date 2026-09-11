Ikawa_2008_imipenem <- function() {
  description <- paste(
    "Three-compartment IV population PK model for imipenem in 10 Japanese",
    "adults with intraabdominal infections (Ikawa 2008). A single 500 mg",
    "infusion into a central compartment with first-order elimination and",
    "linear distribution to two peripheral compartments. No covariates were",
    "retained. Inter-individual variability is exponential on clearance and",
    "all three volumes; residual error is additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Ikawa 2008 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Ikawa K, Morikawa N, Ikeda K, Ohge H, Sueda T.",
    "Development of breakpoints of carbapenems for intraabdominal",
    "infections based on pharmacokinetics and pharmacodynamics in",
    "peritoneal fluid.",
    "J Infect Chemother. 2008;14(4):330-332. doi:10.1007/s10156-008-0624-1.",
    "Parameters transcribed from Zhang P, Zhao Y, Zhu J, Yang Y, Liang G,",
    "Wang X, Yu Z. Population pharmacokinetics of imipenem in different",
    "populations for individualized dosing: a systematic review.",
    "Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055",
    "(Tables 1-3, Supplementary Table S1).",
    sep = " "
  )
  vignette <- "Zhang_2025_imipenem_model_review"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE because the primary publication is
  # not on disk; the review does not describe the assayed matrix beyond
  # "blood samples" and HPLC-UV (Supplementary Table S1).
  compartmentData <- list(
    central     = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  # No covariates. Zhang 2025 Table 3 records "NR" (no record) in every
  # column for this study -- covariate analysis method, covariates screened,
  # covariates incorporated, and formulation are all unreported by the
  # review. This is an absence of reporting in the secondary source rather
  # than positive evidence that the primary screened nothing.
  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 10L,
    n_studies        = 1L,
    age_mean         = "43.7 +/- 14.9 years (mean +/- SD)",
    weight_mean      = "56.7 +/- 10.5 kg (mean +/- SD)",
    sex_female_pct   = NA_real_,
    race_ethnicity   = NULL,
    disease_state    = "Adults with intraabdominal infections",
    dose_range       = paste(
      "500 mg imipenem as a single intravenous dose (Zhang 2025",
      "Supplementary Table S1). The infusion duration is not reported by",
      "the review."
    ),
    regions          = "Japan",
    n_concentrations = NA_integer_,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 1). Sample size (number",
      "of concentration records) and sex split are recorded as 'NR' in the",
      "review. Plasma samples were drawn at 0.5, 1, 2, 3, 4, 5 and 6 h after",
      "administration and assayed by HPLC-UV (Zhang 2025 Supplementary Table",
      "S1). Fitted in NONMEM with Crystal Ball 2000 (Zhang 2025 Table 2).",
      "The review records no model-evaluation method for this study.",
      "UNRESOLVED: IS V3 PERITONEAL FLUID? Two lines of evidence point that",
      "way. Zhang 2025's Discussion says the three-compartment structure",
      "matters for 'studies investigating imipenem PK in nonplasma",
      "compartments, for example, peritoneal fluid (Ikawa et al., 2008)';",
      "and the primary's own title, recovered from the review's reference",
      "list, is 'Development of breakpoints of carbapenems for",
      "intraabdominal infections based on pharmacokinetics and",
      "pharmacodynamics IN PERITONEAL FLUID'. Against that, the only",
      "definition the review gives of the symbol is its Table 2 footnote,",
      "which reads 'V3, volume of peripheral compartment 2' -- a generic",
      "distribution compartment -- and the review reports a single additive",
      "residual error rather than the two a plasma-plus-peritoneal-fluid",
      "model would need. Both peripheral compartments are therefore encoded",
      "as ordinary distribution compartments, which is what the review",
      "actually states, and neither is identified with peritoneal fluid.",
      "If the primary turns out to fit a peritoneal-fluid compartment with",
      "its own observable, this model is structurally incomplete rather",
      "than merely unverified, and the third state would need the canonical",
      "role sorted out before re-extraction. This is the single strongest",
      "reason to obtain the primary for this particular study.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Ikawa 2008 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Ikawa et al. (2008) row =====
    # The review gives the typical values only; no relative standard errors,
    # confidence intervals or bootstrap results are reproduced.
    #
    # Compartment mapping. The Zhang 2025 Table 2 footnote defines V1 as the
    # 'volume of central compartment', V2 as 'volume of peripheral
    # compartment/peripheral compartment 1', V3 as 'volume of peripheral
    # compartment 2', Q2 as 'intercompartmental clearance (peripheral
    # compartment 1)' and Q3 as 'intercompartmental clearance (peripheral
    # compartment 2)'. So the paper's (V2, Q2) pair maps to the canonical
    # (vp, q) and the (V3, Q3) pair maps to (vp2, q2).
    lcl  <- log(9.42); label("Clearance (L/h)")                                 # Zhang 2025 Table 2 (Ikawa 2008): CL = 9.42 L/h
    lvc  <- log(4.66); label("Central volume of distribution V1 (L)")           # Zhang 2025 Table 2 (Ikawa 2008): V1 = 4.66 L
    lvp  <- log(5.08); label("First peripheral volume of distribution V2 (L)")  # Zhang 2025 Table 2 (Ikawa 2008): V2 = 5.08 L
    lvp2 <- log(4.57); label("Second peripheral volume of distribution V3 (L)") # Zhang 2025 Table 2 (Ikawa 2008): V3 = 4.57 L
    lq   <- log(5.74); label("Intercompartmental clearance to peripheral1 Q2 (L/h)") # Zhang 2025 Table 2 (Ikawa 2008): Q2 = 5.74 L/h
    lq2  <- log(31.6); label("Intercompartmental clearance to peripheral2 Q3 (L/h)") # Zhang 2025 Table 2 (Ikawa 2008): Q3 = 31.6 L/h

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 heads this column
    # 'Inter-individual variability (IIV)' and prints a percentage per
    # parameter, without stating whether the percentage is an apparent
    # coefficient of variation of the log-normal distribution or an omega
    # standard deviation expressed as a percent. The review reproduces each
    # primary's own reported quantity without harmonising conventions --
    # measured directly against the four constituent models already
    # extracted in this package from their primaries, the column turns out
    # to mix at least three conventions (see the vignette's 'Assumptions and
    # deviations' section for the full audit). Every model transcribed from
    # this review therefore adopts ONE uniform, documented reading: the
    # printed percentage is treated as an apparent CV of a log-normal
    # random effect, so omega^2 = log(1 + CV^2). This is exact for a source
    # reporting NONMEM-style CV%, and understates omega^2 modestly for a
    # source reporting an omega SD as a percent (at CV = 30% the two
    # readings differ by 4% in variance; the gap widens at large CV). It is
    # recorded here, in every sibling model, and in the vignette Errata so
    # that a later re-extraction from the primary can correct it.
    etalcl  ~ log(1 + 0.269^2)  # Zhang 2025 Table 2 (Ikawa 2008): IIV CL = 26.9%, read as an apparent CV
    etalvc  ~ log(1 + 0.475^2)  # Zhang 2025 Table 2 (Ikawa 2008): IIV V1 = 47.5%, read as an apparent CV
    etalvp  ~ log(1 + 0.371^2)  # Zhang 2025 Table 2 (Ikawa 2008): IIV V2 = 37.1%, read as an apparent CV
    etalvp2 ~ log(1 + 0.862^2)  # Zhang 2025 Table 2 (Ikawa 2008): IIV V3 = 86.2%, read as an apparent CV

    # ===== Residual error =====
    # Additive-only. The review prints the additive term in mg/L, i.e. on
    # the concentration scale, so it is taken as a standard deviation --
    # which is what nlmixr2's add() expects. No proportional term is
    # reported for this study.
    addSd <- 1.13; label("Additive residual error (mg/L)")  # Zhang 2025 Table 2 (Ikawa 2008): Additive = 1.13 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # No covariates (Zhang 2025 Table 3 records 'NR' throughout for this
    # study), so each parameter is its typical value times an exponential
    # random effect.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q   <- exp(lq)
    q2  <- exp(lq2)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column. No absorption compartment.
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ----- Output -----
    # Dose in mg, vc in L -> mg/L, the units the review reports the additive
    # residual error in.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
