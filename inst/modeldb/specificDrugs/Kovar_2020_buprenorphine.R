Kovar_2020_buprenorphine <- function() {
  description <- "Three-compartment intravenous PK model for buprenorphine in adults (typical values fitted in NONMEM to digitised mean literature profiles), scaled to children and preterm neonates by fixed-exponent allometry on body weight (Tod et al. 2008 approach). This is the classical allometric comparator of the Kovar 2020 paper; the paper's PK-Sim whole-body PBPK model is not included."
  reference <- "Kovar L, Schrapel C, Selzer D, Kohl Y, Bals R, Schwab M, Lehr T. Physiologically-Based Pharmacokinetic (PBPK) Modeling of Buprenorphine in Adults, Children and Preterm Neonates. Pharmaceutics. 2020;12(6):578. doi:10.3390/pharmaceutics12060578"
  vignette <- "Kovar_2020_buprenorphine"
  units <- list(time = "min", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size scaling of every clearance and volume, reference weight 71 kg (adult internal dataset; Supplementary Materials Section 3).",
      source_name = "BW"
    )
  )

  compartmentData <- list(
    central = list(analyte = "buprenorphine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "buprenorphine", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "buprenorphine", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 21,
    n_studies = 4,
    age_range = "adults: mean 64.2 years (Bullingham 1982), 32-39 years (Huestis 2013), 27-40 years (Kuhlman 1996), not reported (Everhart 1999); scaled to children 4.6-7.5 years and preterm neonates 27-34 weeks postmenstrual age",
    weight_range = "adults: mean 66.4 kg (Bullingham 1982), 62.1-82.6 kg (Huestis 2013), 62.6-72.7 kg (Kuhlman 1996); reference weight 71 kg",
    sex_female_pct = NA_real_,
    disease_state = "Adult volunteers and patients receiving intravenous buprenorphine (internal training dataset of the PBPK analysis); paediatric predictions for postoperative children (Olkkola 1989) and ventilated preterm neonates (Barrett 1993).",
    dose_range = "0.3-16 mg intravenous (1-60 min infusions) in adults; 3 ug/kg in children; 3 ug/kg loading plus 0.72-2.16 ug/kg/h for 11-118 h in preterm neonates",
    notes = "The three-compartment parameters were estimated in NONMEM 7.4.3 on the internal (training) dataset of the PBPK analysis: digitised mean profiles from Bullingham 1982 (2) (0.3 mg), Everhart 1999 (1 mg), Huestis 2013 (1) and (5) (2 and 16 mg) and Kuhlman 1996 (1.2 mg); main-text Table 1. n_subjects = 5 + 6 + 5 + 5 unique participants (the two Huestis 2013 arms share the same 5 participants). Arterial and venous samples were pooled. No inter-individual or residual variability was reported for this fit."
  )

  ini({
    # Supplementary Materials Table S4, row 'Adults (internal dataset)'. The
    # paper reports clearances in mL/min; the values are divided by 1000 here to
    # give L/min, consistent with volumes in L and time in min.
    lcl <- log(982.0 / 1000); label("Elimination clearance, 71 kg adult (L/min)") # Table S4: CL 982.0 mL/min
    lvc <- log(29.6); label("Central volume of distribution, 71 kg adult (L)") # Table S4: Vc 29.6 L
    lq <- log(2980.0 / 1000); label("Intercompartmental clearance central-peripheral1, 71 kg adult (L/min)") # Table S4: Q2 2980.0 mL/min
    lvp <- log(105.0); label("Peripheral1 volume of distribution, 71 kg adult (L)") # Table S4: V2 105.0 L
    lq2 <- log(554.0 / 1000); label("Intercompartmental clearance central-peripheral2, 71 kg adult (L/min)") # Table S4: Q3 554.0 mL/min
    lvp2 <- log(676.0); label("Peripheral2 volume of distribution, 71 kg adult (L)") # Table S4: V3 676.0 L

    # Allometric exponents, fixed by construction (Supplementary Materials
    # Equations S9-S14). For preterm neonates the paper also scaled CL with the
    # age-dependent exponent 1.2 of Mahmood and Tegenge (Equation S15); set
    # e_wt_cl = 1.2 to reproduce that variant.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Eq S9
    e_wt_q_q2 <- fixed(0.75); label("Allometric exponent on Q and Q2 (unitless)") # Eq S10-S11
    e_wt_vc_vp_vp2 <- fixed(1); label("Allometric exponent on Vc, Vp and Vp2 (unitless)") # Eq S12-S14
  })

  model({
    wt_ratio <- WT / 71

    cl <- exp(lcl) * wt_ratio^e_wt_cl
    q <- exp(lq) * wt_ratio^e_wt_q_q2
    q2 <- exp(lq2) * wt_ratio^e_wt_q_q2
    vc <- exp(lvc) * wt_ratio^e_wt_vc_vp_vp2
    vp <- exp(lvp) * wt_ratio^e_wt_vc_vp_vp2
    vp2 <- exp(lvp2) * wt_ratio^e_wt_vc_vp_vp2

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ug / L = ng/mL
    Cc <- central / vc
  })
}
