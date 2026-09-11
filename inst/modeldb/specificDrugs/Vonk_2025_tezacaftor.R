Vonk_2025_tezacaftor <- function() {
  description <- paste(
    "Joint parent + metabolite population pharmacokinetic model for oral",
    "tezacaftor (CFTR corrector, given as the tezacaftor-ivacaftor",
    "combination Symkevi / Symdeko) and its main metabolite M1 in 21",
    "children with cystic fibrosis aged 6-17 years, from the real-world",
    "prospective SYM-CF study (Vonk 2025). Both analytes are described by",
    "two-compartment disposition; the parent depot receives the oral dose",
    "as a zero-order input of duration D1 and then transfers to the",
    "central compartment at first-order rate KA. All of the apparent",
    "tezacaftor clearance forms M1 (fm fixed at 1), so the M1 compartment",
    "is driven by the whole parent elimination flux. Apparent clearances",
    "and volumes are CL/F, Q/F and V/F for the parent and CL/(F*fm),",
    "Q/(F*fm) and V/(F*fm) for the metabolite; metabolite concentrations",
    "were expressed as parent equivalents using the molecular weight.",
    "Body weight is the only covariate, applied as fixed allometric",
    "scaling (exponent 0.75 on CL and Q, 1 on Vc and Vp, reference 70 kg)",
    "to both parent and metabolite. Because the paediatric data were",
    "sparse, the model was fitted with the NONMEM PRIOR subroutine using",
    "adolescent/adult priors from the Symdeko registration document; CL",
    "and its IIV were estimated without a prior.")
  reference <- "Vonk SEM, Terheggen-Lagro SWJ, Haarman EG, Janssens HM, Maitland-van der Zee AH, Kemper EM, Mathot RAA. Real-world population pharmacokinetics of tezacaftor-ivacaftor in children with cystic fibrosis: The SYM-CF study. Br J Clin Pharmacol. 2025;91(10):2969-2978. doi:10.1002/bcp.70131"
  vignette <- "Vonk_2025_tezacaftor_ivacaftor"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Applied as fixed allometric scaling to the parent and the metabolite alike: CL and Q scale with (WT/70)^0.75 and Vc and Vp with (WT/70)^1 (Vonk 2025 Table 2 footnote a). Weight scaling was predefined rather than selected by the covariate search; the covariates actually screened on CL (age, adherence, CF mutation) showed no relationship (Vonk 2025 Section 3.2.1). Study weight range 23.6-69.8 kg, median 43.5 kg (Table 1).",
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    depot           = list(analyte = "tezacaftor",    units = "mg", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "tezacaftor",    units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "tezacaftor",    units = "mg", specimen = "plasma", verified = TRUE),
    central_m1      = list(analyte = "tezacaftor-M1", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m1  = list(analyte = "tezacaftor-M1", units = "mg", specimen = "plasma", verified = TRUE)
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
    dose_range     = "Tezacaftor 50 mg once daily (6-11 years, <30 kg) or 100 mg once daily (6-11 years >=30 kg, and 12-17 years), per the Symkevi product information; given with the twice-daily ivacaftor component.",
    regions        = "The Netherlands (three Dutch hospitals; enrolment May 2021 to August 2022)",
    dosing_groups  = "6-11 years <30 kg, n = 3 (14%); 6-11 years >=30 kg, n = 7 (33%); 12-17 years, n = 11 (52%) (Vonk 2025 Table 1)",
    n_observations = "97 PK samples in total across all five analytes (13 plasma, 84 dried blood spot); median 5 samples per patient (range 2-7). Three samples (3%) were excluded for incorrect DBS sampling or missing dosing information.",
    notes          = "Prospective, observational, multicentre real-world PK study (SYM-CF), METC Amsterdam UMC ABR NL75811.018.21. Concentrations were quantified by LC-MS/MS (LLOQ 0.01 mg/L, ULOQ 10 mg/L; tezacaftor-M1 LLOQ 0.025 mg/L, ULOQ 12.5 mg/L). Dried blood spot concentrations were converted to estimated plasma concentrations by a Passing-Bablok regression. Fitted in NONMEM 7.5.1 with FOCE-I using the PRIOR subroutine; adolescent/adult prior information came from the Symdeko registration document. Patients had used tezacaftor-ivacaftor for at least two weeks before inclusion, so all data are at steady state."
  )

  ini({
    # ------------------------------------------------------------------
    # TEZACAFTOR (parent) STRUCTURAL PARAMETERS
    # Values are the "Estimates Value (RSE)" column of Vonk 2025 Table 2,
    # normalised to a body weight of 70 kg. Prior types recorded in the
    # comments come from the same table (V = vague 10^5, M = moderately
    # informative RSE 30%, I = informative RSE 10%); the prior weight
    # governed estimation and is provenance, not an encodable quantity.
    # ------------------------------------------------------------------

    lcl <- log(1.95)
    label("Tezacaftor apparent clearance CL/F at 70 kg (L/h)")                    # Table 2 tezacaftor CL = 1.95 (RSE 6%); no prior

    lvc <- log(38.4)
    label("Tezacaftor apparent central volume Vc/F at 70 kg (L)")                 # Table 2 tezacaftor Vc = 38.4 (RSE 7%); vague prior

    lq <- log(0.19)
    label("Tezacaftor apparent inter-compartmental clearance Q/F at 70 kg (L/h)") # Table 2 tezacaftor Q = 0.19 (RSE 29%); moderately informative prior

    lvp <- log(36.4)
    label("Tezacaftor apparent peripheral volume Vp/F at 70 kg (L)")              # Table 2 tezacaftor Vp = 36.4 (RSE 29%); moderately informative prior

    lka <- log(2.95)
    label("Tezacaftor first-order absorption rate constant Ka (1/h)")             # Table 2 tezacaftor Ka = 2.95 (RSE 10%); informative prior

    ld1 <- log(1.06)
    label("Tezacaftor zero-order absorption duration D1 into the depot (h)")      # Table 2 tezacaftor D1 = 1.06 (RSE 28%); moderately informative prior

    # ------------------------------------------------------------------
    # TEZACAFTOR-M1 (metabolite) STRUCTURAL PARAMETERS
    # Apparent values, i.e. CL/(F*fm) and V/(F*fm), on the parent-
    # equivalent molar basis the paper used for the metabolite assay.
    # ------------------------------------------------------------------

    lcl_m1 <- log(1.01)
    label("Tezacaftor-M1 apparent clearance CL/(F*fm) at 70 kg (L/h)")            # Table 2 tezacaftor-M1 CL = 1.01 (RSE 6%); no prior

    lvc_m1 <- log(4.86)
    label("Tezacaftor-M1 apparent central volume Vc/(F*fm) at 70 kg (L)")         # Table 2 tezacaftor-M1 Vc = 4.86 (RSE 28%); informative prior

    lq_m1 <- log(3.70)
    label("Tezacaftor-M1 apparent inter-compartmental clearance Q/(F*fm) at 70 kg (L/h)")  # Table 2 tezacaftor-M1 Q = 3.70 (RSE 19%); moderately informative prior

    lvp_m1 <- log(37.5)
    label("Tezacaftor-M1 apparent peripheral volume Vp/(F*fm) at 70 kg (L)")      # Table 2 tezacaftor-M1 Vp = 37.5 (RSE 10%); informative prior

    # ------------------------------------------------------------------
    # FRACTION METABOLISED
    # ------------------------------------------------------------------

    fm <- fixed(1)
    label("Fraction of tezacaftor clearance forming M1 (unitless)")               # Vonk 2025 Section 2.3: "For tezacaftor the fraction parent drug metabolized into the metabolite was fixed to 1 for fm,M1"

    # ------------------------------------------------------------------
    # ALLOMETRIC SCALING (shared by parent and metabolite)
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

    etalcl ~ log(1 + 0.26^2)                                                      # Table 2 tezacaftor IIV CL = 26 CV% (RSE 19%, shrinkage 6%)

    etalcl_m1 ~ log(1 + 0.24^2)                                                   # Table 2 tezacaftor-M1 IIV CL = 24 CV% (RSE 19%, shrinkage 5%)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # Vonk 2025 Equation 4: Y = IPRED * (1 + theta_prop) + theta_add.
    # Only the proportional term was retained for both analytes.
    # ------------------------------------------------------------------

    propSd <- 0.26
    label("Tezacaftor proportional residual error (fraction)")                    # Table 2 tezacaftor prop. error = 0.26 (RSE 9%)

    propSd_m1 <- 0.20
    label("Tezacaftor-M1 proportional residual error (fraction)")                 # Table 2 tezacaftor-M1 prop. error = 0.20 (RSE 9%)
  })

  model({
    # 1. Allometric size terms, shared by the parent and the metabolite
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
    q_m1  <- exp(lq_m1)              * allom_cl
    vp_m1 <- exp(lvp_m1)             * allom_v

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_m1 <- cl_m1 / vc_m1
    k12_m1 <- q_m1  / vc_m1
    k21_m1 <- q_m1  / vp_m1

    # 4. ODE system (Vonk 2025 Figure 1). The oral dose enters `depot` as
    #    a zero-order input over D1 and leaves it at first-order rate KA.
    #    The whole apparent parent elimination flux forms M1 because fm is
    #    fixed at 1; the fm term is carried explicitly so the structure
    #    reads the same as the ivacaftor sibling model, where fm < 1.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - kel * central -
                              k12 * central + k21 * peripheral1
    d/dt(peripheral1)    <-  k12 * central - k21 * peripheral1

    d/dt(central_m1)     <-  fm * kel * central - kel_m1 * central_m1 -
                              k12_m1 * central_m1 + k21_m1 * peripheral1_m1
    d/dt(peripheral1_m1) <-  k12_m1 * central_m1 - k21_m1 * peripheral1_m1

    # 5. Absorption input function
    dur(depot) <- d1

    # 6. Observations. Volumes are in L and doses in mg, so both
    #    concentrations come out in mg/L, the unit the paper reports.
    Cc    <- central    / vc
    Cc_m1 <- central_m1 / vc_m1

    Cc    ~ prop(propSd)
    Cc_m1 ~ prop(propSd_m1)
  })
}
