Vonk_2025_ivacaftor <- function() {
  description <- paste(
    "Joint parent + two-metabolite population pharmacokinetic model for",
    "oral ivacaftor (CFTR potentiator, given as the tezacaftor-ivacaftor",
    "combination Symkevi / Symdeko) and its main metabolites M1 and M6 in",
    "21 children with cystic fibrosis aged 6-17 years, from the real-world",
    "prospective SYM-CF study (Vonk 2025). Ivacaftor is described by a",
    "two-compartment model whose depot receives the oral dose as a",
    "zero-order input of duration D1 and then transfers to the central",
    "compartment at first-order rate KA. M1 and M6 are one-compartment",
    "analytes formed from the parent central compartment by first-order",
    "rate constants, with the fractions metabolised fixed at 22% and 43%",
    "respectively; no metabolite model existed in the literature, so their",
    "apparent central volumes were fixed at 0.1 times the apparent",
    "ivacaftor central volume. Apparent clearances and volumes are CL/F,",
    "Q/F and V/F for the parent and CL/(F*fm) and V/(F*fm) for the",
    "metabolites; metabolite concentrations were expressed as parent",
    "equivalents using the molecular weight. Body weight is the only",
    "covariate, applied as fixed allometric scaling (exponent 0.75 on CL",
    "and Q, 1 on Vc and Vp, reference 70 kg) to parent and metabolites",
    "alike. Separate proportional residual errors were estimated for the",
    "M1 and M6 observations drawn as venous plasma and as dried blood",
    "spots. Because the paediatric data were sparse, the model was fitted",
    "with the NONMEM PRIOR subroutine using adolescent/adult priors from",
    "the Symdeko registration document; CL and its IIV were estimated",
    "without a prior.")
  reference <- "Vonk SEM, Terheggen-Lagro SWJ, Haarman EG, Janssens HM, Maitland-van der Zee AH, Kemper EM, Mathot RAA. Real-world population pharmacokinetics of tezacaftor-ivacaftor in children with cystic fibrosis: The SYM-CF study. Br J Clin Pharmacol. 2025;91(10):2969-2978. doi:10.1002/bcp.70131"
  vignette <- "Vonk_2025_tezacaftor_ivacaftor"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Residual-SD names carry a two-token suffix (metabolite then sample
  # matrix) because the source estimates one proportional error per
  # metabolite per collection matrix; the canonical propSd_<output>
  # matcher recognises only the single-token metabolite form.
  paper_specific_residual_sds <- c(
    "propSd_m1_plasma", "propSd_m1_dbs",
    "propSd_m6_plasma", "propSd_m6_dbs"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only structural covariate retained in the final model. Applied as fixed allometric scaling to the parent and to both metabolites: CL and Q scale with (WT/70)^0.75 and Vc and Vp with (WT/70)^1 (Vonk 2025 Table 2 footnote a). Weight scaling was predefined rather than selected by the covariate search; the covariates actually screened on CL (age, adherence, CF mutation) showed no relationship (Vonk 2025 Section 3.2.2). Study weight range 23.6-69.8 kg, median 43.5 kg (Table 1).",
      source_name        = "BW"
    ),
    SAMPLE_CAPILLARY = list(
      description        = "Dried blood spot (capillary) sampling indicator for the observation record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (venous plasma sample)",
      notes              = "1 = the ivacaftor-M1 / ivacaftor-M6 concentration came from a finger-prick dried blood spot (converted to an estimated plasma concentration by the Passing-Bablok regression of Vonk's earlier method paper); 0 = the concentration was measured directly in venous plasma. Per-observation indicator; a subject contributes both kinds of record. Used only to switch the proportional residual-error magnitude for the two metabolite outputs (Vonk 2025 Section 2.3: 'For ivacaftor-M1 and M6 separate proportional error models for plasma and DBS samples were implemented'); the parent ivacaftor observation has a single residual error covering both matrices. Of the 97 study samples, 13 (13%) were plasma and 84 (87%) were DBS (Table 1). Set to 1 for the whole simulation when reproducing the paper's DBS-dominated real-world design.",
      source_name        = "DBS"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "ivacaftor",    units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ivacaftor",    units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ivacaftor",    units = "mg", specimen = "plasma", verified = TRUE),
    central_m1  = list(analyte = "ivacaftor-M1", units = "mg", specimen = "plasma", verified = TRUE),
    central_m6  = list(analyte = "ivacaftor-M6", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "6 to 17 years; median 12 years (Vonk 2025 Table 1)",
    age_median     = "12 years",
    weight_range   = "23.6 to 69.8 kg (Vonk 2025 Table 1)",
    weight_median  = "43.5 kg",
    height_range   = "122 to 191 cm; median 153 cm (Vonk 2025 Table 1)",
    sex_female_pct = 52,
    disease_state  = "Children with cystic fibrosis carrying at least one F508del mutation: 16 (76%) homozygous F508del, 5 (24%) heterozygous F508del (4 with A455E, 1 with 3849+10kbC>T). Exocrine pancreatic insufficiency in 20 (95%), distal intestinal obstruction syndrome in 4 (19%), no CF-related diabetes.",
    dose_range     = "Ivacaftor 75 mg twice daily (6-11 years, <30 kg; 150 mg/day) or 150 mg twice daily (6-11 years >=30 kg, and 12-17 years; 300 mg/day), per the Symkevi product information; given with the once-daily tezacaftor component.",
    regions        = "The Netherlands (three Dutch hospitals; enrolment May 2021 to August 2022)",
    dosing_groups  = "6-11 years <30 kg, n = 3 (14%); 6-11 years >=30 kg, n = 7 (33%); 12-17 years, n = 11 (52%) (Vonk 2025 Table 1)",
    n_observations = "97 PK samples in total across all five analytes (13 plasma, 84 dried blood spot); median 5 samples per patient (range 2-7). Three samples (3%) were excluded for incorrect DBS sampling or missing dosing information, and three ivacaftor-M6 samples (3%) were below the LLOQ and excluded.",
    notes          = "Prospective, observational, multicentre real-world PK study (SYM-CF), METC Amsterdam UMC ABR NL75811.018.21. Concentrations were quantified by LC-MS/MS (LLOQ 0.01 mg/L, ULOQ 10 mg/L). Dried blood spot concentrations were converted to estimated plasma concentrations by a Passing-Bablok regression. Fitted in NONMEM 7.5.1 with FOCE-I using the PRIOR subroutine; adolescent/adult prior information came from the Symdeko registration document, which contained no ivacaftor metabolite model. Patients had used tezacaftor-ivacaftor for at least two weeks before inclusion, so all data are at steady state."
  )

  ini({
    # ------------------------------------------------------------------
    # IVACAFTOR (parent) STRUCTURAL PARAMETERS
    # Values are the "Estimates Value (RSE)" column of Vonk 2025 Table 2,
    # normalised to a body weight of 70 kg. Prior types recorded in the
    # comments come from the same table (V = vague 10^5, M = moderately
    # informative RSE 30%, I = informative RSE 10%); the prior weight
    # governed estimation and is provenance, not an encodable quantity.
    # ------------------------------------------------------------------

    lcl <- log(15.9)
    label("Ivacaftor apparent clearance CL/F at 70 kg (L/h)")                     # Table 2 ivacaftor CL = 15.9 (RSE 9%); no prior

    lvc <- log(178)
    label("Ivacaftor apparent central volume Vc/F at 70 kg (L)")                  # Table 2 ivacaftor Vc = 178 (RSE 20%); vague prior

    lq <- log(13.2)
    label("Ivacaftor apparent inter-compartmental clearance Q/F at 70 kg (L/h)")  # Table 2 ivacaftor Q = 13.2 (RSE 25%); moderately informative prior

    lvp <- log(106)
    label("Ivacaftor apparent peripheral volume Vp/F at 70 kg (L)")               # Table 2 ivacaftor Vp = 106 (RSE 26%); moderately informative prior

    lka <- log(0.506)
    label("Ivacaftor first-order absorption rate constant Ka (1/h)")              # Table 2 ivacaftor Ka = 0.506 (RSE 10%); informative prior

    ld1 <- log(2.59)
    label("Ivacaftor zero-order absorption duration D1 into the depot (h)")       # Table 2 ivacaftor D1 = 2.59 (RSE 10%); informative prior

    # ------------------------------------------------------------------
    # IVACAFTOR-M1 AND IVACAFTOR-M6 (metabolites)
    # Apparent values, i.e. CL/(F*fm) and V/(F*fm), on the parent-
    # equivalent molar basis the paper used for the metabolite assays.
    # The metabolite volumes were not estimable; both were fixed at
    # 0.1 * Vc,ivacaftor, and the arithmetic is kept inline so the rule
    # and the parent volume it was applied to stay visible.
    # ------------------------------------------------------------------

    lcl_m1 <- log(2.10)
    label("Ivacaftor-M1 apparent clearance CL/(F*fm) at 70 kg (L/h)")             # Table 2 ivacaftor-M1 CL = 2.10 (RSE 10%); no prior

    lvc_m1 <- fixed(log(0.1 * 178))
    label("Ivacaftor-M1 apparent central volume Vc/(F*fm) at 70 kg (L)")          # Table 2 ivacaftor-M1 Vc = 0.1 * Viva; Section 3.2.2 "Values for Viva-M1 and Viva-M6 were fixed at 0.1*Viva"

    lcl_m6 <- log(12.2)
    label("Ivacaftor-M6 apparent clearance CL/(F*fm) at 70 kg (L/h)")             # Table 2 ivacaftor-M6 CL = 12.2 (RSE 16%); no prior

    lvc_m6 <- fixed(log(0.1 * 178))
    label("Ivacaftor-M6 apparent central volume Vc/(F*fm) at 70 kg (L)")          # Table 2 ivacaftor-M6 Vc = 0.1 * Viva; Section 3.2.2 "Values for Viva-M1 and Viva-M6 were fixed at 0.1*Viva"

    # ------------------------------------------------------------------
    # FRACTIONS METABOLISED
    # ------------------------------------------------------------------

    fm_m1 <- fixed(0.22)
    label("Fraction of ivacaftor clearance forming M1 (unitless)")                # Vonk 2025 Section 2.3: "For ivacaftor the fraction parent drug metabolized into the metabolites were fixed to 22 and 43% for fm,M1 and fm,M6"

    fm_m6 <- fixed(0.43)
    label("Fraction of ivacaftor clearance forming M6 (unitless)")                # Vonk 2025 Section 2.3: fm,M6 = 43%

    # ------------------------------------------------------------------
    # ALLOMETRIC SCALING (shared by parent and metabolites)
    # ------------------------------------------------------------------

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of (WT/70) on CL and Q (unitless)")                # Table 2 footnote a: CL = thetaCL * (weight/70)^0.75 and Q = thetaQ * (weight/70)^0.75

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of (WT/70) on Vc and Vp (unitless)")               # Table 2 footnote a: Vc/p = thetaV * (weight/70)^1

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # Table 2 reports IIV on CL as a percent CV for the exponential model
    # CL = thetaCL * exp(eta_CL) (Vonk 2025 Equation 3), so the variance
    # on the eta scale is omega^2 = log(CV^2 + 1). See the vignette
    # Assumptions and deviations section for the arithmetic that
    # discriminates this reading from omega = CV.
    # ------------------------------------------------------------------

    etalcl ~ log(1 + 0.40^2)                                                      # Table 2 ivacaftor IIV CL = 40 CV% (RSE 37%, shrinkage 5%)

    etalcl_m1 ~ log(1 + 0.44^2)                                                   # Table 2 ivacaftor-M1 IIV CL = 44 CV% (RSE 38%, shrinkage 3%)

    etalcl_m6 ~ log(1 + 0.76^2)                                                   # Table 2 ivacaftor-M6 IIV CL = 76 CV% (RSE 37%, shrinkage 3%)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # Vonk 2025 Equation 4: Y = IPRED * (1 + theta_prop) + theta_add.
    # Only the proportional term was retained. The parent has a single
    # error; each metabolite has one error per collection matrix.
    # ------------------------------------------------------------------

    propSd <- 0.34
    label("Ivacaftor proportional residual error (fraction)")                     # Table 2 ivacaftor prop. error = 0.34 (RSE 10%)

    propSd_m1_plasma <- 0.37
    label("Ivacaftor-M1 proportional residual error, venous plasma samples (fraction)")  # Table 2 ivacaftor-M1 prop. error plasma = 0.37 (RSE 24%)

    propSd_m1_dbs <- 0.36
    label("Ivacaftor-M1 proportional residual error, dried blood spot samples (fraction)")  # Table 2 ivacaftor-M1 prop. error DBS = 0.36 (RSE 10%)

    propSd_m6_plasma <- 0.98
    label("Ivacaftor-M6 proportional residual error, venous plasma samples (fraction)")  # Table 2 ivacaftor-M6 prop. error plasma = 0.98 (RSE 25%)

    propSd_m6_dbs <- 0.50
    label("Ivacaftor-M6 proportional residual error, dried blood spot samples (fraction)")  # Table 2 ivacaftor-M6 prop. error DBS = 0.50 (RSE 10%)
  })

  model({
    # 1. Allometric size terms, shared by the parent and both metabolites
    #    (Vonk 2025 Table 2 footnote a, Equations 1 and 2).
    ref_wt   <- 70  # kg
    allom_cl <- (WT / ref_wt)^e_wt_cl_q
    allom_v  <- (WT / ref_wt)^e_wt_vc_vp

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc)          * allom_v
    q  <- exp(lq)           * allom_cl
    vp <- exp(lvp)          * allom_v
    ka <- exp(lka)
    d1 <- exp(ld1)

    cl_m1 <- exp(lcl_m1 + etalcl_m1) * allom_cl
    vc_m1 <- exp(lvc_m1)             * allom_v
    cl_m6 <- exp(lcl_m6 + etalcl_m6) * allom_cl
    vc_m6 <- exp(lvc_m6)             * allom_v

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_m1 <- cl_m1 / vc_m1
    kel_m6 <- cl_m6 / vc_m6

    # 4. ODE system (Vonk 2025 Figure 2). The oral dose enters `depot` as
    #    a zero-order input over D1 and leaves it at first-order rate KA.
    #    Each metabolite is formed from the parent central compartment at
    #    the corresponding fraction of the apparent parent elimination
    #    flux and is eliminated by its own first-order clearance.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                           k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_m1)  <-  fm_m1 * kel * central - kel_m1 * central_m1
    d/dt(central_m6)  <-  fm_m6 * kel * central - kel_m6 * central_m6

    # 5. Absorption input function
    dur(depot) <- d1

    # 6. Observations. Volumes are in L and doses in mg, so all three
    #    concentrations come out in mg/L, the unit the paper reports.
    Cc    <- central    / vc
    Cc_m1 <- central_m1 / vc_m1
    Cc_m6 <- central_m6 / vc_m6

    # Per-record residual magnitude: the metabolite errors switch on the
    # collection matrix (Vonk 2025 Section 2.3 and Table 2).
    w_m1 <- propSd_m1_plasma * (1 - SAMPLE_CAPILLARY) +
            propSd_m1_dbs    * SAMPLE_CAPILLARY
    w_m6 <- propSd_m6_plasma * (1 - SAMPLE_CAPILLARY) +
            propSd_m6_dbs    * SAMPLE_CAPILLARY

    Cc    ~ prop(propSd)
    Cc_m1 ~ prop(w_m1)
    Cc_m6 ~ prop(w_m6)
  })
}
