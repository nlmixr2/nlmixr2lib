Danielak_2020_treosulfan_rat <- function() {
  description <- "Preclinical (rat). Joint parent-metabolite population PK model for treosulfan (TREO) and its active monoepoxide (S,S)-1,2-epoxybutane-3,4-diol-4-methanesulfonate (EBDM) in plasma and brain of Wistar rats after a single 500 mg/kg intraperitoneal dose (Danielak 2020). First-order absorption into a one-compartment TREO plasma model; irreversible first-order conversion of TREO to EBDM (fixed rate constants in plasma and in brain); one-compartment EBDM plasma model; bidirectional blood-brain barrier transport of TREO and EBDM parameterised by an influx clearance and the influx/efflux clearance ratio; a second (deep) brain compartment for TREO. All clearances and volumes are apparent (divided by F) and per kg body weight. Male sex lowers TREO plasma clearance by 14.6%. IIV on ka and TREO CL; proportional residual error fixed to 15% (plasma) and 20% (brain)."
  reference <- "Danielak D, Romanski M, Kasprzyk A, Tezyk A, Glowka F. Population pharmacokinetic approach for evaluation of treosulfan and its active monoepoxide disposition in plasma and brain on the basis of a rat model. Pharmacol Rep. 2020;72(5):1297-1309. doi:10.1007/s43440-020-00115-0. Correction: Pharmacol Rep. 2020;72:1443. doi:10.1007/s43440-020-00144-9 (units of the first-order rate constants corrected from L/h to 1/h; no values changed)"
  vignette <- "Danielak_2020_treosulfan_rat"
  units <- list(time = "h", dosing = "umol/kg", concentration = "umol/L")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "female (SEXF = 1)",
      notes = "The paper's covariate is a male indicator: CL1-MALE/F = CL1/F + CL1/F * COV_CL-MALE (Table 1 footnote a), with the typical CL1/F (0.419 L/h/kg) referring to females. Encoded as cl <- tvcl * (1 + e_sex_cl * (1 - SEXF)) so the published coefficient (-0.146) keeps its sign. Typical male CL1/F = 0.419 * (1 - 0.146) = 0.358 L/h/kg, matching the male median in Figure 6.",
      source_name = "MALE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "treosulfan", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "treosulfan", units = "umol/kg", specimen = "plasma", verified = TRUE),
    central_ebdm = list(
      analyte = "EBDM (treosulfan monoepoxide)",
      units = "umol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    brain_extravascular = list(analyte = "treosulfan", units = "umol/kg", specimen = "tissue", verified = TRUE),
    brain_deep = list(analyte = "treosulfan", units = "umol/kg", specimen = "tissue", verified = TRUE),
    brain_extravascular_ebdm = list(
      analyte = "EBDM (treosulfan monoepoxide)",
      units = "umol/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 96L,
    n_studies = 1L,
    age_range = "10 weeks",
    weight_mean = "306 +/- 25 g (males); 188 +/- 15 g (females)",
    sex_female_pct = 50,
    disease_state = "Healthy animals",
    dose_range = "Single 500 mg/kg (1797 umol/kg) treosulfan intraperitoneal injection",
    regions = "Poland (Poznan)",
    notes = "One-animal-per-sample (destructive) design: 48 male and 48 female 10-week-old Wistar rats; plasma and brain sampled predose and at 0.25, 0.5, 1, 2, 4, 6 and 24 h (Methods 'Animals', 'Sample collection', Fig. 2). One female with an outlying 0.25 h plasma TREO concentration was excluded and the 24 h time point (all below LLOQ) was dropped (Results, first paragraph). Unbound fraction of TREO and EBDM in plasma and brain homogenate >= 0.94, so measured concentrations were treated as unbound. Brain concentrations were converted from homogenate supernatant to tissue concentrations (1 uM supernatant = 6 umol/kg tissue; brain density 1.04 g/mL)."
  )

  ini({
    # Absorption (intraperitoneal, first-order) -- Table 1
    lka <- log(5.12); label("First-order absorption rate constant k12 (1/h)") # Table 1: k12 = 5.12 1/h (RSE 17%)

    # Treosulfan plasma disposition -- Table 1 (apparent, per kg)
    lcl <- log(0.419); label("Treosulfan plasma clearance CL1/F in females (L/h/kg)") # Table 1: CL1/F = 0.419 L/h/kg (RSE 5%)
    lvc <- fixed(log(1.03)); label("Treosulfan plasma volume V2/F (L/kg)") # Table 1: V2/F = 1.03 L/kg (fixed)

    # Conversion of treosulfan to EBDM
    lkmet <- fixed(log(0.451)); label("Treosulfan-to-EBDM conversion rate constant in plasma k23 (1/h)") # Table 1: k23 = 0.451 1/h (fixed); Methods assumption (4), in vitro value at pH 7.40
    lkmet_brain <- fixed(log(0.271)); label("Treosulfan-to-EBDM conversion rate constant in brain k45 (1/h)") # Table 1: k45 = 0.271 1/h (fixed); Methods assumption (5), log10(kf) = -7.479 + 0.960 * pH at pH 7.2

    # EBDM plasma disposition -- Table 1
    lcl_ebdm <- log(6.24); label("EBDM plasma clearance CL2/F (L/h/kg)") # Table 1: CL2/F = 6.24 L/h/kg (RSE 3%)
    lvc_ebdm <- fixed(log(0.914)); label("EBDM plasma volume V3/F (L/kg)") # Table 1: V3/F = 0.914 L/kg (fixed)

    # Treosulfan brain -- Table 1
    lclin <- log(0.0233); label("Treosulfan blood-brain barrier influx clearance Q1/F = CLin (L/h/kg)") # Table 1: Q1/F = 0.0233 L/h/kg (RSE 14%)
    lkp_brain <- log(0.120); label("Treosulfan BBB1 = CLin/CLout, influx/efflux clearance ratio (unitless)") # Table 1: BBB1 = 0.120 (RSE 5%)
    lv_brain_extravascular <- fixed(log(6.56e-3)); label("Treosulfan central brain volume V4/F (L/kg)") # Table 1: V4/F = 6.56 x 10^-3 L/kg (fixed to rat brain extravascular space)
    lq_brain_deep <- log(0.304); label("Treosulfan central-to-peripheral brain clearance Q2/F (L/h/kg)") # Table 1: Q2/F = 0.304 L/h/kg (RSE 19%)
    lv_brain_deep <- fixed(log(0.364)); label("Treosulfan peripheral brain volume V6/F (L/kg)") # Table 1: V6/F = 0.364 L/kg (fixed)

    # EBDM brain -- Table 1
    lclin_ebdm <- log(0.0122); label("EBDM blood-brain barrier influx clearance Q3/F = CLin (L/h/kg)") # Table 1: Q3/F = 0.0122 L/h/kg (RSE 10%)
    lkp_brain_ebdm <- log(0.317); label("EBDM BBB2 = CLin/CLout, influx/efflux clearance ratio (unitless)") # Table 1: BBB2 = 0.317 (RSE 5%)
    lv_brain_extravascular_ebdm <- fixed(log(6.56e-3)); label("EBDM brain volume V5/F (L/kg)") # Table 1: V5/F = 6.56 x 10^-3 L/kg (fixed)

    # Sex covariate on treosulfan plasma clearance
    e_sex_cl <- -0.146; label("Fractional change in CL1/F for males (unitless)") # Table 1: COV_CL-MALE = -0.146; footnote a: CL1-MALE/F = CL1/F + CL1/F * COV_CL-MALE

    # IIV -- Table 1 reports %CV = sqrt(exp(OMEGA) - 1) * 100 (footnote b);
    # OMEGA = log(1 + (CV/100)^2)
    etalka ~ 0.38941 # Table 1: omega k12 = 69% CV (RSE 17%); log(1 + 0.69^2) = 0.38941
    etalcl ~ 0.00995 # Table 1: omega CL = 10% CV (RSE 35%); log(1 + 0.10^2) = 0.00995

    # Residual error -- proportional, fixed to the assay error (Methods,
    # 'Population pharmacokinetic analysis'; Table 1)
    propSd <- fixed(0.15); label("Proportional residual SD, treosulfan plasma (fraction)") # Table 1: 0.15 (fixed) for TREO and EBDM in plasma
    propSd_ebdm <- fixed(0.15); label("Proportional residual SD, EBDM plasma (fraction)") # Table 1: 0.15 (fixed) for TREO and EBDM in plasma
    propSd_Cbrain <- fixed(0.20); label("Proportional residual SD, treosulfan brain (fraction)") # Table 1: 0.20 (fixed) for TREO and EBDM in brain tissue
    propSd_Cbrain_ebdm <- fixed(0.20); label("Proportional residual SD, EBDM brain (fraction)") # Table 1: 0.20 (fixed) for TREO and EBDM in brain tissue
  })

  model({
    # Individual parameters (Methods: theta_ij = theta_j * exp(eta_ij))
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (1 + e_sex_cl * (1 - SEXF))
    vc <- exp(lvc)
    kmet <- exp(lkmet)
    kmet_brain <- exp(lkmet_brain)
    cl_ebdm <- exp(lcl_ebdm)
    vc_ebdm <- exp(lvc_ebdm)
    clin <- exp(lclin)
    kp_brain <- exp(lkp_brain)
    v_brain_extravascular <- exp(lv_brain_extravascular)
    q_brain_deep <- exp(lq_brain_deep)
    v_brain_deep <- exp(lv_brain_deep)
    clin_ebdm <- exp(lclin_ebdm)
    kp_brain_ebdm <- exp(lkp_brain_ebdm)
    v_brain_extravascular_ebdm <- exp(lv_brain_extravascular_ebdm)

    # Blood-brain barrier efflux clearances (Methods assumption (7):
    # BBB = CLin / CLout)
    clef <- clin / kp_brain
    clef_ebdm <- clin_ebdm / kp_brain_ebdm

    # Concentrations (umol/L)
    Cc <- central / vc
    Cc_ebdm <- central_ebdm / vc_ebdm
    Cbrain <- brain_extravascular / v_brain_extravascular
    Cbrain_deep <- brain_deep / v_brain_deep
    Cbrain_ebdm <- brain_extravascular_ebdm / v_brain_extravascular_ebdm

    # ODEs (Figure 3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - kmet * central - clin * Cc + clef * Cbrain
    d/dt(central_ebdm) <- kmet * central - cl_ebdm * Cc_ebdm - clin_ebdm * Cc_ebdm + clef_ebdm * Cbrain_ebdm
    d/dt(brain_extravascular) <- clin * Cc - clef * Cbrain - q_brain_deep * (Cbrain - Cbrain_deep) - kmet_brain * brain_extravascular
    d/dt(brain_deep) <- q_brain_deep * (Cbrain - Cbrain_deep)
    d/dt(brain_extravascular_ebdm) <- kmet_brain * brain_extravascular + clin_ebdm * Cc_ebdm - clef_ebdm * Cbrain_ebdm

    Cc ~ prop(propSd)
    Cc_ebdm ~ prop(propSd_ebdm)
    Cbrain ~ prop(propSd_Cbrain)
    Cbrain_ebdm ~ prop(propSd_Cbrain_ebdm)
  })
}
