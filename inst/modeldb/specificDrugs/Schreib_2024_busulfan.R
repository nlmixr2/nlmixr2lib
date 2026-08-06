Schreib_2024_busulfan <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model with intravenous",
    "infusion and a time-varying elimination rate constant for busulfan in",
    "124 pediatric patients undergoing hematopoietic stem cell",
    "transplantation (HSCT) at the University Children's Hospital Zurich,",
    "October 2010 to February 2020, receiving twice-daily (q12h) 2-hour to",
    "4-hour intravenous busulfan infusions as conditioning. Schreib 2024 fit",
    "ln(V) and ln(k) directly rather than CL and V, because 'building the",
    "model on V and CL would result in two parameters with similar",
    "covariates ... potentially masking the covariates for k'. Two further",
    "parameters describe an exponential decline of the elimination rate",
    "constant over the course of therapy: kel(t) = kel(0) * (1 +",
    "kel_exp_famp * (1 - exp(-kel_exp_kdes * t))), with kel_exp_famp = -0.167",
    "(a 16.7% average fall in k and CL at steady state) and a 13.4 h",
    "half-life for the change. The fractional amplitude carries its own",
    "inter-individual variability and is doubled to -0.312 (a 31% fall) in",
    "the hemophagocytic lymphohistiocytosis / X-linked lymphoproliferative",
    "disease group. Volume scales near-proportionally with calculated total",
    "body water; the elimination rate constant depends on total body water, a",
    "postmenstrual-age maturation function, serum albumin, an acute",
    "lymphoblastic leukemia indicator, and the infusion duration. Total body",
    "water and the maturation function are derived inside model() from",
    "weight, height, age, and sex, so no separate columns are required.",
    "IMPORTANT: simulate this model with rxSolve(..., useLinCmt = FALSE).",
    "It is a one-compartment linear-elimination model, so rxSolve()'s default",
    "ODE-to-linCmt() auto-conversion replaces the ODE with a closed-form",
    "solution that holds the elimination rate constant at its t = 0 value and",
    "therefore silently discards the time dependence that is the entire point",
    "of this paper. The error is large and one-sided: exposure per dosing",
    "interval is under-predicted by about 17% at steady state (45% in the",
    "HLH/XLP group). See the validation vignette for the demonstration."
  )
  reference <- paste(
    "Schreib KM, Braem DS, Zeilhofer UB, Mueller D, Guengoer T, Kraemer SD,",
    "Hauri-Hohl MM. Population Pharmacokinetic Modeling for Twice-Daily",
    "Intravenous Busulfan in a Large Cohort of Pediatric Patients Undergoing",
    "Hematopoietic Stem Cell Transplantation-A 10-Year Single-Center",
    "Experience. Pharmaceutics. 2024;16(1):13.",
    "doi:10.3390/pharmaceutics16010013.",
    "Author R code deposit (saemix objective function BuFunD.R and the",
    "time-dependence explainer kchange.R):",
    "https://gitlab.ethz.ch/skraemer/busulfan_2022.git"
  )
  vignette <- "Schreib_2024_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters only through calculated total body water (TBW, Equation 7 of",
        "Schreib 2024, after Wells 2005). Not a covariate in its own right in",
        "the final model: Table 2 step 2a shows ln(TBW) on ln(V) outperformed",
        "ln(W), ln(BSA), ln(H), and ln(FFM). Population median 17.2 kg",
        "(range 4.3-85.0 kg, Table 1)."
      ),
      source_name        = "W"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters only through calculated total body water (Equation 7).",
        "Schreib 2024 back-calculated height from the recorded body surface",
        "area and weight as H = BSA^2 * 3600 / W (the Mosteller relation,",
        "Equation 6, inverted); the deposited script BuSaemix.R computes",
        "Height that way. Median BSA 0.71 m2 (range 0.25-2.06, Table 1)."
      ),
      source_name        = "H"
    ),
    AGE = list(
      description        = "Postnatal age at transplantation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used twice inside model(): as a term of the total body water",
        "equation (Equation 7) and to build postmenstrual age for the",
        "maturation function Fmat (Equation 9), PMA[weeks] = 52 * AGE + 40.",
        "The 40-week offset is stated in the Methods; the 52-weeks-per-year",
        "conversion is the one used in the deposited BuSaemix.R",
        "(Fmat = 1/(1 + ((AGE + 40/52)/(46/52))^-2.3)). Median 4.3 years",
        "(range 0.2-27.0 years, Table 1); 18% of patients were under 1 year."
      ),
      source_name        = "age at HSCT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Enters only through the total body water equation (Equation 7),",
        "where the female term contributes -0.047 on the natural-log scale.",
        "The deposited BuSaemix.R codes Sex01 = 1 for female and 0 for male,",
        "matching the canonical SEXF orientation with no value transformation.",
        "Sex had no additional effect on busulfan PK beyond TBW (Section 3.2,",
        "Figures S17-S19). 35 of 124 patients (28%) were female (Table 1)."
      ),
      source_name        = "Sex01"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Covariate on ln(k) as (ln(ALB) - ln(30 g/L)); reference value 30 g/L",
        "(Table 3), a rounded value close to the population median. The paper",
        "reports albumin in g/L, matching the canonical unit, so no conversion",
        "is needed. Range of the covariate term in the population -0.15 to",
        "0.16 (Table 3)."
      ),
      source_name        = "Alb"
    ),
    TINF = list(
      description        = "Duration of the intravenous busulfan infusion",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Covariate on both ln(V) (+0.226) and ln(k) (-0.161), entering as the",
        "linear deviation (TINF - 3 h); reference 3 h (Table 3). The center",
        "infused over 4 h from October 2010 to September 2014 and over 3 h",
        "from October 2014; only these two values occur, so the covariate term",
        "is 0 or 1 in this cohort. The opposite",
        "signs mean CL = k * V is virtually independent of TINF; the authors",
        "hypothesise the 3 h infusion uncovers a distribution phase that a",
        "one-compartment model absorbs into apparent k and V (Section 3.2).",
        "TINF must be supplied as a data column in addition to setting the",
        "physical infusion duration of the dose record, because rxode2 takes",
        "the infusion duration from the event table rather than from a",
        "covariate."
      ),
      source_name        = "Tinf"
    ),
    DIS_ALL = list(
      description        = "Acute lymphoblastic leukemia indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other HSCT indication in the cohort)",
      notes              = paste(
        "Covariate on ln(k), theta_k5 = -0.210, i.e. k and CL are about 19%",
        "lower in patients transplanted for acute lymphoblastic leukemia",
        "(Table 3). 13 of 124 patients (10%, Table 1). The reference group is",
        "all other HSCT indications in the cohort, not a matched comparator."
      ),
      source_name        = "DiaGroup == 'ALL'"
    ),
    DIS_HLHXLP = list(
      description        = paste(
        "Hemophagocytic lymphohistiocytosis / X-linked lymphoproliferative",
        "disease indicator"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other HSCT indication in the cohort)",
      notes              = paste(
        "Covariate on the fractional amplitude of the change in the",
        "elimination rate constant (kel_exp_famp), theta_dk2 = -0.145, which",
        "takes the amplitude from -0.167 to -0.312, i.e. from a 17% to a 31%",
        "fall in k and CL over the course of therapy (Table 3, Figure 4B).",
        "It does not affect the elimination rate constant at t = 0. Schreib",
        "2024 pools HLH and XLP into a single 14-patient diagnosis group",
        "(11%, Table 1) and estimates one coefficient for it."
      ),
      source_name        = "DiaGroup == 'HLH/XLP'"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 124L,
    n_studies      = 1L,
    age_range      = "0.2-27.0 years",
    age_median     = "4.3 years",
    weight_range   = "4.3-85.0 kg",
    weight_median  = "17.2 kg",
    sex_female_pct = 28,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Pediatric and young-adult patients conditioned with intravenous",
      "busulfan before allogeneic hematopoietic stem cell transplantation.",
      "42 of 124 (34%) had malignant disease (ALL 13, AML 12, neuroblastoma",
      "6, other 11) and 82 (66%) non-malignant disease (chronic",
      "granulomatous disease 32, HLH or XLP 14, primary immunodeficiencies",
      "14, hemoglobinopathies 12, metabolic diseases 8, thrombocytopenia 2)."
    ),
    dose_range     = paste(
      "Twice-daily (q12h) intravenous busulfan infusions over 4 h (October",
      "2010 to September 2014) or 3 h (from October 2014), four to ten doses",
      "over two to five consecutive days, dosed by weight-based",
      "nomogram and then adjusted by therapeutic drug monitoring toward a",
      "cumulative AUC target. Conditioning regimens were busulfan combined",
      "with fludarabine (66%), fludarabine plus thiotepa (10%), clofarabine",
      "(8%), melphalan (6%), cyclophosphamide plus melphalan (5%), or other."
    ),
    regions        = "Switzerland (single center, University Children's Hospital Zurich).",
    notes          = paste(
      "Baseline demographics from Table 1 of Schreib 2024. Plasma busulfan",
      "was measured by LC-MS/MS after the morning infusion only, pre-dose and",
      "at 0, 30, 60, 120, 240, and 360 min after the end of that infusion.",
      "18% of patients were under 1 year of age, 77% between 1 and 18 years,",
      "and 5% over 18 years. Age is postnatal age at transplantation."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Volume of distribution. Schreib 2024 Table 3, gray equation box:
    #   ln(V) = theta_V1
    #         + theta_V2 * (ln(TBW) - ln(10 L))
    #         + theta_V3 * (Tinf - 3 h)
    # Reference values are TBW = 10 L and Tinf = 3 h.
    # ------------------------------------------------------------------
    lvc       <- 2.459;  label("log of central volume of distribution at the reference covariate values (V, L)")   # Table 3, theta_V1 = 2.459 (SE 0.0189); untransformed 11.70 L
    e_tbw_vc  <- 0.931;  label("Exponent of total body water on the central volume (unitless)")                    # Table 3, theta_V2 = 0.931 (SE 0.0233, p < 0.001)
    e_tinf_vc <- 0.226;  label("Effect of infusion duration above 3 h on log central volume (per h)")              # Table 3, theta_V3 = 0.226 (SE 0.0356, p < 0.001)

    # ------------------------------------------------------------------
    # Elimination rate constant at the start of therapy. Table 3, gray box:
    #   ln(k) = theta_k1
    #         + theta_k2 * (ln(TBW) - ln(10 L))
    #         + theta_k3 * ln(Fmat)
    #         + theta_k4 * (ln(Alb) - ln(30 g/L))
    #         + theta_k5 * ALL
    #         + theta_k6 * (Tinf - 3 h)
    #
    # NOTE ON A SOURCE DEFECT. The multiplicative rewrite printed on the same
    # line of Table 3 reads "k = e^theta_k1 * Fmat^theta_k2 *
    # (TBW/(10 L))^theta_k3 * ...", i.e. it transposes the theta_k2 and
    # theta_k3 exponents relative to the additive form directly above it.
    # The additive form is the correct one; the multiplicative rewrite is a
    # typesetting error. Three independent checks agree: (a) Table 3's own row
    # labels read "theta_k2, effect of dln(TBW, L)" with reference value 10 L
    # and "theta_k3, effect of ln(Fmat)" with reference value 0; (b) Table 2
    # assigns theta_k2 to model 15 ("ln(TBW) on ln(k)") and theta_k3 to model
    # 16** ("ln(Fmat) on ln(k)"); (c) the tabulated population range of the
    # theta_k3 term is -0.40 to 0, which is only possible for a non-negative
    # coefficient on ln(Fmat), since Fmat <= 1 makes ln(Fmat) <= 0.
    # ------------------------------------------------------------------
    lkel       <- -1.007; label("log of the elimination rate constant at the start of therapy at the reference covariate values (1/h)") # Table 3, theta_k1 = -1.007 (SE 0.0329); untransformed 0.365 1/h
    e_tbw_kel  <- -0.189; label("Exponent of total body water on the elimination rate constant (unitless)")        # Table 3, theta_k2 = -0.189 (SE 0.0411, p < 0.001)
    e_fmat_kel <-  0.697; label("Exponent of the age maturation function on the elimination rate constant (unitless)") # Table 3, theta_k3 = 0.697 (SE 0.2001, p < 0.001)
    e_alb_kel  <-  0.331; label("Exponent of serum albumin on the elimination rate constant (unitless)")           # Table 3, theta_k4 = 0.331 (SE 0.1073, p = 0.001)
    e_all_kel  <- -0.210; label("Effect of acute lymphoblastic leukemia on the log elimination rate constant (unitless)") # Table 3, theta_k5 = -0.210 (SE 0.0592, p < 0.001)
    e_tinf_kel <- -0.161; label("Effect of infusion duration above 3 h on the log elimination rate constant (per h)")     # Table 3, theta_k6 = -0.161 (SE 0.0420, p < 0.001)

    # ------------------------------------------------------------------
    # Time dependence of the elimination rate constant. Schreib 2024
    # Equation (2):
    #   k*(t) = k * (-dk * (exp(-kappa_k * t) - 1) + 1)
    #         = k * (1 +  dk * (1 - exp(-kappa_k * t)))
    # so k falls (dk < 0) from k at t = 0 to k * (1 + dk) as t -> infinity.
    # kel_exp_famp is the paper's dk, the dimensionless fractional amplitude
    # of the change; lkel_exp_kdes is the paper's ln(kappa_k).
    #
    # Naming ratified by the operator on 2026-08-05 (task oare_PMC11154452
    # sidecar request-001 / response-001, question q1 answer A): mirror the
    # registered cl_exp_ time-varying-clearance family onto kel, founding the
    # new `_famp` fractional-amplitude role token.
    # ------------------------------------------------------------------
    kel_exp_famp          <- -0.167; label("Fractional amplitude of the change in the elimination rate constant over therapy (unitless)") # Table 3, theta_dk1 = -0.167 (SE 0.0191)
    e_hlhxlp_kel_exp_famp <- -0.145; label("Additional fractional amplitude in the HLH/XLP diagnosis group (unitless)")                   # Table 3, theta_dk2 = -0.145 (SE 0.0540, p = 0.0035)
    lkel_exp_kdes         <- -2.965; label("log of the rate constant of the exponential change in the elimination rate constant (1/h)")   # Table 3, theta_kappa_k = -2.965 (SE 0.1094); untransformed 0.0516 1/h, half-life 13.4 h

    # ------------------------------------------------------------------
    # Inter-individual variability, reported as variances Omega^2 in Table 3.
    # Random effects were estimated on ln(V), ln(k), and dk only; random
    # effects on ln(kappa_k) were dropped because their shrinkage exceeded
    # 60%, indicating overparameterization (Section 3.2 and Table 2, model 2
    # vs model 3*). No covariances between random effects are reported, so
    # the matrix is diagonal. The random effect on kel_exp_famp is additive
    # on the natural scale, because dk is estimated untransformed and is
    # negative.
    # ------------------------------------------------------------------
    etalvc          ~ 0.029; label("IIV on the log central volume (variance)")                        # Table 3, Omega^2 of ln(V) = 0.029 (Omega 0.169, SE 0.0039, shrinkage 5.5%)
    etalkel         ~ 0.033; label("IIV on the log elimination rate constant (variance)")             # Table 3, Omega^2 of ln(k) = 0.033 (Omega 0.182, SE 0.0048, shrinkage 12.5%)
    etakel_exp_famp ~ 0.032; label("IIV on the fractional amplitude of the change in kel (variance)") # Table 3, Omega^2 of dk = 0.032 (Omega 0.180, SE 0.0047, shrinkage 12.5%)

    # ------------------------------------------------------------------
    # Residual error. Schreib 2024 used an additive (constant) error model
    # throughout; a proportional model failed to compute -2LL and gave an
    # asymmetric residual distribution (Sections 2.5 and 3.3).
    # ------------------------------------------------------------------
    addSd <- 0.076; label("Additive residual standard deviation (mg/L)")  # Table 3, residual error = 0.076 (SE 0.0008)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariates
    # ------------------------------------------------------------------
    # Total body water (L), Schreib 2024 Equation (7) after Wells 2005. The
    # female term is -0.047 and SEXF is 1 for female, matching the deposited
    # BuSaemix.R (Sex01 = 1 for female).
    tbw <- exp(-2.952 + 0.551 * log(WT) + 0.796 * log(HT) + 0.008 * AGE - 0.047 * SEXF)

    # Age maturation function for busulfan elimination, Equation (9) after
    # Hassine 2021. Postmenstrual age is postnatal age plus a fixed 40 weeks
    # (Section 2.7); the deposited BuSaemix.R uses 52 weeks per year, so
    # PMA[weeks] = 52 * AGE + 40. Fmat is 0.5 at 46 weeks PMA and approaches
    # 1 in older children.
    pma  <- 52 * AGE + 40
    fmat <- 1 / (1 + (pma / 46)^(-2.3))

    # ------------------------------------------------------------------
    # 2. Individual parameters
    # ------------------------------------------------------------------
    # Central volume. Reference total body water 10 L, reference infusion
    # duration 3 h (Table 3).
    vc <- exp(lvc + etalvc) * (tbw / 10)^e_tbw_vc * exp(e_tinf_vc * (TINF - 3))

    # Elimination rate constant at the start of therapy (t = 0). Reference
    # values: TBW 10 L, Fmat 1, albumin 30 g/L, no ALL, infusion duration 3 h.
    kel <- exp(lkel + etalkel) * (tbw / 10)^e_tbw_kel * fmat^e_fmat_kel *
      (ALB / 30)^e_alb_kel * exp(e_all_kel * DIS_ALL + e_tinf_kel * (TINF - 3))

    # Individual fractional amplitude of the change in the elimination rate
    # constant, dk in Table 3. Additive on the natural scale, additive
    # HLH/XLP effect, additive random effect.
    kel_exp_famp_i <- kel_exp_famp + etakel_exp_famp +
      (e_hlhxlp_kel_exp_famp * DIS_HLHXLP)

    # Rate constant of the exponential change in kel (1/h).
    kel_exp_kdes <- exp(lkel_exp_kdes)

    # ------------------------------------------------------------------
    # 3. Time-varying elimination rate constant, Equation (2)
    # ------------------------------------------------------------------
    # kel_exp_total is what the ODE consumes; kel above is its value at
    # t = 0 and kel * (1 + kel_exp_famp_i) is its asymptote.
    #
    # Schreib 2024 prints two different time-varying rate constants and only
    # one of them is the model. Equation (2) is the instantaneous rate
    # constant k*(t) used here. Equation (3) defines k'(t), the running
    # time-average of k*(t) obtained by integration; it exists only because
    # saemix needed a closed-form superposition likelihood with no ODE
    # solver, and it is materially different (k*/k = 0.881 vs k'/k = 0.929 at
    # 24 h). Figure 4B, which plots simulated CL(t) normalized to CL(t = 0),
    # settles the choice: its solid "mean" line reads about 0.87 at 24 h,
    # 0.84 at 48 h, and 0.837 at 60 h, matching Equation (2) (0.881, 0.847,
    # 0.841) and not Equation (3) (0.929, 0.895, 0.885); its dashed HLH/XLP
    # line reads about 0.775, 0.72, and 0.70, matching Equation (2) at an
    # amplitude of -0.312 (0.779, 0.714, 0.702). The deposited kchange.R
    # calls the stepwise integration of k*(t) "the reference".
    #
    # `time` is continuous here. The deposited fit function BuFunD.R
    # discretizes it to the dosing-interval midpoint
    # (TimeAv <- (ceiling(Time/12) - 0.5) * 12), which would make this curve
    # start at 0.956 rather than 1 at t = 0+; Figure 4B starts at exactly 1
    # and is smooth, so the published curve uses continuous time. Continuous
    # time is also regimen-independent, whereas a hard-coded 12 h interval
    # would silently corrupt any non-q12h simulation.
    kel_exp_total <- kel * (1 + kel_exp_famp_i * (1 - exp(-kel_exp_kdes * time)))

    # ------------------------------------------------------------------
    # 4. ODE system. One compartment with intravenous infusion; the infusion
    # duration comes from the event table (see the TINF covariate note).
    #
    # SOLVE THIS WITH rxSolve(..., useLinCmt = FALSE). This is a
    # one-compartment linear-elimination system, so rxSolve()'s default
    # useLinCmt = TRUE rewrites it as a closed-form linCmt() solution, which
    # can only carry a CONSTANT rate constant and silently substitutes
    # kel_exp_total's t = 0 value for all time. Verified against a hand
    # integration: with the default flag, the AUC of every dosing interval
    # collapses to exactly Dose / (kel(0) * vc) and the entire time
    # dependence disappears. No way of writing the ODE avoids the
    # conversion -- an explicit rate variable and a clearance-form
    # denominator were both tested and are converted identically -- so the
    # solve-time flag is the only mitigation.
    # ------------------------------------------------------------------
    d/dt(central) <- -kel_exp_total * central

    # ------------------------------------------------------------------
    # 5. Observation and error
    # ------------------------------------------------------------------
    # central is in mg and vc in L, so Cc is in mg/L, the units the paper
    # and the deposited scripts use for busulfan plasma concentrations.
    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
