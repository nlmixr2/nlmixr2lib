Sulaiman_2026_piperacillin_tazobactam <- function() {
  description <- paste(
    "Joint one-compartment population PK model for piperacillin and",
    "tazobactam in 45 critically ill adults with sepsis or septic shock",
    "in three Malaysian intensive care units (Sulaiman 2026).",
    "Each drug has its own first-order-elimination central compartment",
    "with clearance and volume of distribution estimated separately;",
    "Cockcroft-Gault creatinine clearance enters both clearances as a",
    "power term normalised to the population median of 71 mL/min.",
    "Between-subject AND between-occasion variability are carried on the",
    "clearance and the volume of both drugs, the two occasions being",
    "day 1 and day 3 of therapy. Residual error is combined additive plus",
    "proportional for each drug. Piperacillin uses the unsuffixed",
    "canonical compartment / parameter set; tazobactam carries the",
    "sibling-drug suffix _taz throughout."
  )
  reference <- paste(
    "Sulaiman H, Wolky SA, Rozali MA, Adiraju SKS, Hasan MS,",
    "Hernandez-Mitre MP, Liu X, Mat-Nor MB, Mazlan MZ, Salmuna ZN,",
    "Wallis SC, Xie J, Roberts JA, Abdul-Aziz MH. A multicentre",
    "evaluation of pharmacokinetic/pharmacodynamic target attainment of",
    "piperacillin and tazobactam and the association with clinical",
    "outcomes in critically ill patients with sepsis and septic shock.",
    "J Antimicrob Chemother. 2026. doi:10.1093/jac/dkag199.",
    sep = " "
  )
  vignette <- "Sulaiman_2026_piperacillin_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    central_taz = list(analyte = "tazobactam",   units = "mg", specimen = "plasma", verified = TRUE)
  )

  dosing <- c("central", "central_taz")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Enters the",
        "clearance of BOTH drugs as a power term normalised to the",
        "population median, exactly as printed in the two display",
        "equations of Sulaiman 2026 Results ('Pharmacokinetic model",
        "building'):",
        "CL_pip = CLpop_pip * (CLcr / 71)^0.87 and",
        "CL_taz = CLpop_taz * (CLcr / 71)^0.90, with CLcr in mL/min.",
        "Adding the term dropped the objective function value by 103.19",
        "points and cut between-subject variability on clearance from 70%",
        "to 24% for piperacillin and from 75% to 30% for tazobactam.",
        "CENTRING-CONSTANT CAVEAT: the Results text describes the",
        "normalisation as 'the population median value of 71 mL/min' and",
        "the printed equations divide by 71, but Table 1 reports the",
        "cohort median Cockcroft-Gault CLcr as 70 mL/min (mean 75, SD 40,",
        "range 30-161). The 71 printed in the equations is used here",
        "because it is the value actually implemented in the covariate",
        "model; the 1 mL/min discrepancy against Table 1 changes a typical",
        "clearance by only 1.2% and is recorded in the vignette Errata.",
        "Cockcroft-Gault was computed from the serum creatinine on ICU",
        "admission (median 81.0 umol/L, range 38.0-196.0). Patients with",
        "pre-existing renal impairment -- renal replacement therapy, CLcr",
        "below 30 mL/min, or plasma creatinine above 200 umol/L -- were",
        "excluded, so the model carries no information below 30 mL/min",
        "and the paper's own dosing simulations span 30-130 mL/min.",
        "Stored under canonical CRCL with raw mL/min, matching the",
        "AbdulAziz_2016_doripenem.R precedent from the same senior",
        "authors and the same Malaysian ICU network."
      ),
      source_name        = "CLcr"
    ),
    OCC = list(
      description        = "Sampling-occasion index for between-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Two occasions: OCC = 1 is day 1 of piperacillin/tazobactam",
        "therapy and OCC = 2 is day 3 (Sulaiman 2026 Methods, 'Drug dosing",
        "and sample collection'; up to seven blood samples per patient per",
        "occasion). Between-occasion variability was estimated on the",
        "clearance and the volume of distribution of BOTH drugs (Table 2),",
        "and the Discussion singles the term out as capturing 'the dynamic",
        "physiological changes during critical illness'. Decomposed inside",
        "model() into the mutually exclusive indicators oc1 / oc2 that",
        "select the per-occasion eta slots, because rxode2 parses but",
        "cannot simulate the `eta ~ var | occ` multi-level syntax; the",
        "second-occasion variances are `fixed()` to the first, which is the",
        "same shared-magnitude structure NONMEM writes as",
        "$OMEGA BLOCK(1) SAME (precedents: Chen_2023_nemonoxacin.R,",
        "Jonsson_2011_ethambutol.R).",
        "Of 45 patients, 26 were sampled on BOTH occasions and 19 on one",
        "occasion only (Results, 'Plasma concentrations'), 381 samples in",
        "total."
      ),
      source_name        = "Occasion"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded at baseline (mean 64 kg, SD 16; median 65, range 31-114;",
        "Table 1) and screened during covariate testing, but NOT retained:",
        "Sulaiman 2026 reports Cockcroft-Gault CLcr on clearance as the",
        "only covariate that improved the fit. The volumes of distribution",
        "therefore carry no allometric term, which the Discussion flags",
        "explicitly ('the V was comparable despite the lower average body",
        "weight of the present cohort'). Documented for provenance only;",
        "not referenced in model()."
      )
    ),
    ALB = list(
      description        = "Serum albumin on ICU admission",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Collected as part of the baseline biochemistry (Methods, 'Data",
        "collection'); mean 24 g/L (SD 5), median 24, range 15-36, with",
        "hypoalbuminaemia (< 25 g/L) in 24/45 (53%) patients (Table 1).",
        "Not retained in the final model despite being a recognised driver",
        "of beta-lactam unbound fraction and distribution volume in",
        "critical illness. Documented for provenance only; not referenced",
        "in model()."
      )
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mean 55 years (SD 21), median 61, range 18-88 (Table 1). Age",
        "separates the target-attainment groups strongly -- patients with",
        "therapeutic concentrations on occasion 1 had a median age of 66",
        "versus 28 in the sub-therapeutic group (Table 3, P = 0.004) --",
        "but that association runs through Cockcroft-Gault CLcr, of which",
        "age is an input, and age itself was not retained as a covariate",
        "in the final PK model. Documented for provenance only; not",
        "referenced in model()."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = "18-88 years (median 61, IQR 36-72; mean 55, SD 21)",
    age_median     = "61 years",
    weight_range   = "31-114 kg (median 65; mean 64, SD 16)",
    weight_median  = "65 kg",
    sex_female_pct = 44,
    race_ethnicity = "Malaysian; the source paper reports no race / ethnicity breakdown",
    disease_state  = paste(
      "Sepsis or septic shock by the Sepsis-3 criteria with confirmed or",
      "suspected bacterial infection, in adults admitted to an intensive",
      "care unit. Primary infection site was lung in 32 (71%), blood in 7",
      "(16%), intra-abdominal in 3 (7%), skin in 2 (4%) and urinary in 1",
      "(2%). Median APACHE II 14 (range 2-35) and median SOFA 6 (range",
      "1-16) on ICU admission; 33 (73%) were mechanically ventilated.",
      "A causative pathogen was identified in 29 (64%), most often",
      "Pseudomonas aeruginosa (n = 7) and Klebsiella pneumoniae (n = 7)",
      "(Table 1, Results)."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance mean 75 mL/min (SD 40),",
      "median 70, range 30-161; serum creatinine on ICU admission median",
      "81.0 umol/L, range 38.0-196.0 (Table 1). Patients receiving renal",
      "replacement therapy, or with CLcr below 30 mL/min or plasma",
      "creatinine above 200 umol/L, were EXCLUDED, so the cohort spans",
      "preserved to augmented renal function and the model carries no",
      "information about renal impairment or dialysis."
    ),
    dose_range     = paste(
      "Piperacillin/tazobactam at the discretion of the treating",
      "physician: 4.5 g every 6 h in 37 (82%), 4.5 g every 8 h in 7 (16%)",
      "and 2.25 g every 6 h in 1 (2%). Administration was by intermittent",
      "infusion in 14 (31%), extended infusion of 3-4 h in 13 (29%) and",
      "continuous infusion in 18 (40%), the latter preceded by a 4.5 g",
      "loading dose infused over 30-60 min. A 4.5 g dose is 4 g",
      "piperacillin plus 0.5 g tazobactam (the clinical 8:1 ratio), which",
      "is how the two central compartments must be dosed (Table 1)."
    ),
    regions        = paste(
      "Malaysia -- three intensive care units: Sultan Ahmad Shah Medical",
      "Centre, Hospital Universiti Sains Malaysia and University of Malaya",
      "Medical Centre"
    ),
    notes          = paste(
      "Prospective multicentre pharmacokinetic study conducted March 2017",
      "to March 2018 (NMRR-16237231675). 381 plasma samples were collected",
      "across two sampling occasions (day 1 and day 3 of therapy), 26",
      "patients contributing both occasions and 19 a single occasion, with",
      "up to seven samples per occasion. Total plasma piperacillin and",
      "tazobactam were quantified by validated UHPLC-MS/MS; observed",
      "piperacillin spanned 0.93-400.42 mg/L and tazobactam 0.63-49.44",
      "mg/L. Estimation used Monolix 2021R2 and the Monte Carlo dosing",
      "simulations Simulx 2024R1.",
      "The PK/PD analysis of the source paper -- which this model file does",
      "NOT encode, because the authors fitted the PK to TOTAL plasma",
      "concentrations -- converts to unbound concentrations by assuming 30%",
      "protein binding for both drugs (i.e. an unbound fraction of 0.70,",
      "Methods 'Dosing simulations and probability of target attainment')",
      "and scores 100% fT > 16 mg/L for piperacillin combined with 85%",
      "fT > 2 mg/L for tazobactam, subject to a piperacillin toxicity",
      "threshold of a total trough or steady-state concentration of",
      "160 mg/L. Applying that free fraction to Cc / Cc_taz reproduces the",
      "paper's target-attainment analysis; see the vignette."
    )
  )

  ini({
    # =====================================================================
    # Structural fixed effects -- Table 2, "Fixed effects" block. Monolix
    # estimated CL and V for each drug separately on the log scale
    # (Results, "Pharmacokinetic model building": "All parameters were log
    # transformed"), so the typical values below are wrapped in log().
    # The reference individual has a Cockcroft-Gault CLcr of 71 mL/min.
    # =====================================================================
    lcl <- log(8.78);  label("Piperacillin clearance at CLcr = 71 mL/min (L/h)")            # Table 2 CLpip = 8.78 L/h (%RSE 5.69; bootstrap median 8.73, 95% CI 7.78-9.92)
    lvc <- log(25.53); label("Piperacillin volume of distribution (L)")                     # Table 2 Vpip = 25.53 L (%RSE 9.52; bootstrap median 25.73, 95% CI 2.17-31.59)

    lcl_taz <- log(8.74);  label("Tazobactam clearance at CLcr = 71 mL/min (L/h)")          # Table 2 CLtaz = 8.74 L/h (%RSE 6.39; bootstrap median 8.7, 95% CI 7.65-9.87)
    lvc_taz <- log(27.15); label("Tazobactam volume of distribution (L)")                   # Table 2 Vtaz = 27.15 L (%RSE 9.77; bootstrap median 27.88, 95% CI 15.95-34.47)

    # =====================================================================
    # Cockcroft-Gault creatinine-clearance effect on clearance, as printed
    # in the two display equations of Results, "Pharmacokinetic model
    # building":
    #   For piperacillin, CL = CLpop x (CLcr / 71)^0.87
    #   For tazobactam,   CL = CLpop x (CLcr / 71)^0.90
    # The 71 mL/min normaliser is a published constant, not an estimated
    # parameter, so it is fixed() here.
    # =====================================================================
    e_crcl_cl     <- 0.87; label("Power exponent of (CRCL / 71) on piperacillin clearance (unitless)") # Table 2 "CLcr effect on CLpip" = 0.87 (%RSE 10.2; bootstrap median 0.89, 95% CI 0.74-1.09); exponent form from the printed equation
    e_crcl_cl_taz <- 0.90; label("Power exponent of (CRCL / 71) on tazobactam clearance (unitless)")   # Table 2 "CLcr effect on CLtaz" = 0.9 (%RSE 10.8; bootstrap median 0.96, 95% CI 0.77-1.18); exponent form from the printed equation
    crcl_ref_cl <- fixed(71); label("Reference Cockcroft-Gault creatinine clearance for the covariate normalisation (mL/min)") # Results: "normalized to the population median value of 71 mL/min"; the printed equations divide by 71 (Table 1 reports the cohort median as 70 -- see the CRCL covariateData note)

    # =====================================================================
    # Between-subject variability -- Table 2, "Between-subject variability"
    # block. Those rows are reported as CV PERCENTAGES, and the variances
    # below reproduce them exactly through the log-normal identity
    # CV% = sqrt(exp(omega^2) - 1) x 100:
    #   omega 0.24 -> 24.35%   omega 0.41 -> 42.78%
    #   omega 0.30 -> 30.69%   omega 0.41 -> 42.78%
    # All four published CV% values invert to omegas that are round to two
    # decimal places (max residual 0.005 percentage points), which is what
    # identifies the "%" column as a CV rather than a raw omega; the same
    # inversion holds for all four between-occasion rows below.
    # =====================================================================
    etalcl     ~ 0.0576; label("BSV on piperacillin clearance (log-scale variance)")            # Table 2 BSV CLpip = 24.35% (%RSE 36.3) -> omega = 0.24, variance 0.24^2 = 0.0576
    etalvc     ~ 0.1681; label("BSV on piperacillin volume of distribution (log-scale variance)") # Table 2 BSV Vpip = 42.78% (%RSE 20.0) -> omega = 0.41, variance 0.41^2 = 0.1681
    etalcl_taz ~ 0.0900; label("BSV on tazobactam clearance (log-scale variance)")              # Table 2 BSV CLtaz = 30.69% (%RSE 29.0) -> omega = 0.30, variance 0.30^2 = 0.09
    etalvc_taz ~ 0.1681; label("BSV on tazobactam volume of distribution (log-scale variance)") # Table 2 BSV Vtaz = 42.78% (%RSE 19.4) -> omega = 0.41, variance 0.41^2 = 0.1681

    # =====================================================================
    # Between-occasion variability -- Table 2, "Between-occasion
    # variability" block, on the clearance and the volume of BOTH drugs.
    # There are exactly two occasions (day 1 and day 3 of therapy), and a
    # single magnitude is reported per parameter, so occasion 2 is fixed()
    # to the occasion-1 variance: the shared-magnitude structure NONMEM
    # writes as $OMEGA BLOCK(1) SAME. Same CV% -> omega inversion as above.
    # =====================================================================
    etaiov_cl_1 ~ 0.0961;        label("BOV on piperacillin clearance, occasion 1 (log-scale variance)")   # Table 2 BOV CLpip = 31.76% (%RSE 19.3) -> omega = 0.31, variance 0.31^2 = 0.0961
    etaiov_cl_2 ~ fixed(0.0961); label("BOV on piperacillin clearance, occasion 2 (log-scale variance)")   # shared magnitude with occasion 1 (one BOV value reported per parameter)
    etaiov_vc_1 ~ 0.0729;        label("BOV on piperacillin volume, occasion 1 (log-scale variance)")      # Table 2 BOV Vpip = 27.50% (%RSE 20.6) -> omega = 0.27, variance 0.27^2 = 0.0729
    etaiov_vc_2 ~ fixed(0.0729); label("BOV on piperacillin volume, occasion 2 (log-scale variance)")      # shared magnitude with occasion 1

    etaiov_cl_taz_1 ~ 0.1024;        label("BOV on tazobactam clearance, occasion 1 (log-scale variance)") # Table 2 BOV CLtaz = 32.84% (%RSE 19.1) -> omega = 0.32, variance 0.32^2 = 0.1024
    etaiov_cl_taz_2 ~ fixed(0.1024); label("BOV on tazobactam clearance, occasion 2 (log-scale variance)") # shared magnitude with occasion 1
    etaiov_vc_taz_1 ~ 0.0576;        label("BOV on tazobactam volume, occasion 1 (log-scale variance)")    # Table 2 BOV Vtaz = 24.35% (%RSE 20.9) -> omega = 0.24, variance 0.24^2 = 0.0576
    etaiov_vc_taz_2 ~ fixed(0.0576); label("BOV on tazobactam volume, occasion 2 (log-scale variance)")    # shared magnitude with occasion 1

    # =====================================================================
    # Residual unexplained variability -- Table 2, "Random error" block.
    # "Residual unexplained variability was described by a combined error
    # model" (Results, "Pharmacokinetic model building"). The paper does
    # not state WHICH of Monolix's two combined forms was used --
    # combined1, SD = a + b*f, or combined2, SD = sqrt(a^2 + (b*f)^2) --
    # and Text S4, which might, is not among the supplementary items
    # released with the open-access PDF. The nlmixr2 default add() + prop()
    # (= combined2) is used and the ambiguity is recorded in the vignette
    # Errata.
    # =====================================================================
    addSd  <- 2.43; label("Piperacillin additive residual SD (mg/L)")     # Table 2 Additive pip = 2.43 mg/L (%RSE 26.6; bootstrap median 2.54, 95% CI 0.71-4.78)
    propSd <- 0.17; label("Piperacillin proportional residual SD (fraction)") # Table 2 Proportional pip = 17% (%RSE 7.72; bootstrap median 17%, 95% CI 12-22)

    addSd_taz  <- 0.40; label("Tazobactam additive residual SD (mg/L)")   # Table 2 Additive taz = 0.4 mg/L (%RSE 24.2; bootstrap median 0.43, 95% CI 0.2-0.8)
    propSd_taz <- 0.16; label("Tazobactam proportional residual SD (fraction)") # Table 2 Proportional taz = 16% (%RSE 8.38; bootstrap median 16%, 95% CI 1-20)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Occasion indicators for the between-occasion variability. Two
    #    occasions only: day 1 (OCC = 1) and day 3 (OCC = 2) of therapy.
    #    Records carrying any other OCC value zero out both indicators and
    #    so receive no BOV, which is the correct behaviour for a
    #    simulation grid that extends past the sampled occasions.
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_cl     <- oc1 * etaiov_cl_1     + oc2 * etaiov_cl_2
    iov_vc     <- oc1 * etaiov_vc_1     + oc2 * etaiov_vc_2
    iov_cl_taz <- oc1 * etaiov_cl_taz_1 + oc2 * etaiov_cl_taz_2
    iov_vc_taz <- oc1 * etaiov_vc_taz_1 + oc2 * etaiov_vc_taz_2

    # ------------------------------------------------------------------
    # 2. Individual PK parameters. Clearance carries the published power
    #    term on Cockcroft-Gault creatinine clearance; the volumes carry
    #    no covariate (Results: CLcr on CL was the only retained effect).
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / crcl_ref_cl)^e_crcl_cl
    vc <- exp(lvc + etalvc + iov_vc)

    cl_taz <- exp(lcl_taz + etalcl_taz + iov_cl_taz) * (CRCL / crcl_ref_cl)^e_crcl_cl_taz
    vc_taz <- exp(lvc_taz + etalvc_taz + iov_vc_taz)

    # ------------------------------------------------------------------
    # 3. Micro-constants.
    # ------------------------------------------------------------------
    kel     <- cl     / vc
    kel_taz <- cl_taz / vc_taz

    # ------------------------------------------------------------------
    # 4. One-compartment IV disposition with first-order elimination per
    #    drug (Results: "A one-compartment model with first-order
    #    elimination best described the pharmacokinetics for both drugs").
    #    The two drugs do not interconvert; each is dosed into its own
    #    central compartment in the clinical 8:1 ratio (4 g piperacillin
    #    with 0.5 g tazobactam per 4.5 g dose).
    # ------------------------------------------------------------------
    d/dt(central)     <- -kel     * central
    d/dt(central_taz) <- -kel_taz * central_taz

    # ------------------------------------------------------------------
    # 5. Observations. Dose in mg with volumes in L gives mg/L = ug/mL,
    #    the units in which the source paper reports every concentration.
    #    These are TOTAL plasma concentrations; the source paper's PK/PD
    #    target attainment multiplies them by an assumed unbound fraction
    #    of 0.70 (30% protein binding) before scoring fT > MIC.
    # ------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_taz <- central_taz / vc_taz

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_taz ~ add(addSd_taz) + prop(propSd_taz)
  })
}
