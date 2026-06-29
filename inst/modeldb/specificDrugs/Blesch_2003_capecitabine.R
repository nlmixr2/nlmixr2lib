Blesch_2003_capecitabine <- function() {
  description <- paste(
    "Population PK model 'CAP7440' for oral capecitabine, evaluated by",
    "simultaneously fitting plasma concentrations of three sequential",
    "metabolites: 5'-DFUR (5'-deoxy-5-fluorouridine), 5-FU (5-fluorouracil),",
    "and FBAL (alpha-fluoro-beta-alanine). Parent capecitabine and the",
    "first metabolite 5'-DFCR are not modelled as compartments; the dose",
    "enters the 5'-DFUR pool directly through a first-order absorption rate",
    "constant KA with absorption lag TLAG. Sequential first-order kinetics",
    "carry mass through 5'-DFUR -> 5-FU -> FBAL using apparent oral",
    "clearances CL/F and volumes V/F (CL1/V1 for 5'-DFUR, CL2/V2 for 5-FU,",
    "CL3/V3 for FBAL). 5-FU volume V2 is fixed at 17.8 L from literature",
    "(Heggie 1987) because it was not sensitive to the dataset. Three",
    "retained covariate effects in the final population PK model: alkaline",
    "phosphatase on 5-FU clearance (CL2), creatinine clearance on FBAL",
    "clearance (CL3) and volume (V3), and body surface area on FBAL volume",
    "(V3) -- all multiplicative power forms. Fit to 481 patients with",
    "advanced or metastatic colorectal cancer (Phase III studies) plus 24",
    "patients with extensive sampling from a Phase I bioequivalence study."
  )
  reference <- "Blesch KS, Gieschke R, Tsukamoto Y, Reigner BG, Burger HU, Steimer JL. Clinical pharmacokinetic/pharmacodynamic and physiologically based pharmacokinetic modeling in new drug development: the capecitabine experience. Invest New Drugs. 2003;21(2):195-223. doi:10.1023/A:1023525513696. PMID: 12889740."
  vignette  <- "Blesch_2003_capecitabine"
  units     <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    ALP = list(
      description        = "Baseline serum alkaline phosphatase activity, multiplicative power covariate on apparent 5-FU clearance CL2/F. The source paper does not state the reference value used for centering; this model file uses 100 U/L as a clinically reasonable median for an advanced-colorectal-cancer cohort, applied as (ALP / 100)^e_alp_cl_5fu. Document the chosen reference in vignette Errata; effect coefficient -0.169 is paper-derived (Table 1).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper reports the effect as multiplicative power: a doubling of baseline ALP (ALP x 2) is associated with a 11% decrease in CL2 and a 12% increase in 5-FU AUC (Table 2). The reference (centering) value is NOT stated in Blesch 2003; the underlying Phase III dataset's median ALP is presumably in references [10] and [11]. The 100 U/L default chosen here is a non-paper provenance value (see vignette Errata).",
      source_name        = "ALP"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance, multiplicative power covariate on apparent FBAL clearance CL3/F and apparent FBAL volume V3/F. The paper does not state the assay (Cockcroft-Gault is the era-typical default) or the reference value; this model file uses 80 mL/min/1.73 m^2 as a clinically reasonable median for the advanced-cancer cohort.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper uses 'CLCR' as the column name with no explicit unit or assay statement (Cockcroft-Gault, raw mL/min is the era-typical assumption). Effect coefficients +0.615 on CL3 and +0.394 on V3 are paper-derived (Table 1); the 80 mL/min/1.73 m^2 reference is a non-paper provenance value (see vignette Errata). A 50% reduction in CLCR (CLCR x 0.5) gives 35% decrease CL3, 24% decrease V3, 53% increase FBAL AUC, 41% increase FBAL Cmax (Table 2).",
      source_name        = "CLCR"
    ),
    BSA = list(
      description        = "Body surface area, multiplicative power covariate on apparent FBAL volume V3/F. The source paper does not state the BSA computation formula (DuBois or Mosteller per era convention) nor the reference value; this model file uses 1.73 m^2 as the standard adult anchor.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Effect coefficient +0.812 on V3 is paper-derived (Table 1); the 1.73 m^2 reference is a non-paper provenance value (see vignette Errata). BSA x 1.3 (a 30% increase) gives 24% increase V3, 19% decrease FBAL Cmax (Table 2). BSA underlies the dosing (1250 mg/m^2 BID in Phase III) but only enters the PK structure as a power covariate on V3.",
      source_name        = "BSA"
    )
  )

  population <- list(
    species            = "human",
    n_subjects         = 505L,
    n_subjects_phase1  = 24L,
    n_subjects_phase3  = 481L,
    n_studies          = 3L,
    age_range          = "Adults with advanced or metastatic solid tumours; specific age range not tabulated in the review.",
    disease_state      = "Advanced or metastatic colorectal cancer (Phase III studies, 481 patients) plus advanced cancer (Phase I bioequivalence study, 24 patients with extensive sampling).",
    dose_range         = "Phase III: oral capecitabine 1250 mg/m^2 twice daily in 3-week treatment cycles (2 weeks on, 1 week off). Sparsely-sampled PK on the first day of cycles 2 and 4. Phase I bioequivalence: full-profile densely sampled data carried forward to stabilise the structural model.",
    regions            = "Phase III in advanced / metastatic colorectal cancer (multi-national); Phase I bioequivalence per Reigner et al. 1998.",
    notes              = "Patient-level demographics (gender, race, body weight, BSA, CLCR, ALP, ALT, AST, total bilirubin, albumin, liver-metastasis status, Karnofsky performance status) were screened as candidate covariates; only ALP, CLCR, and BSA met the inclusion criteria (p < 0.05 univariate; >10% individuals outside 0.8-1.25 reference range; backwards-deletion p < 0.001). The dataset combined sparse Phase III samples with the extensively-sampled Phase I bioequivalence cohort to anchor the structural model. NONMEM objective function value: 2844.378. Source: Blesch et al. 2003 Table 1 (citing Phase II in metastatic breast cancer per Reigner 1999 [20] for the structural model and Phase III colorectal per Hoff 2001 [10] and Van Cutsem 2001 [11] for the Phase III data; numeric demographic distributions live in those source publications, not in the Blesch 2003 review itself)."
  )

  ini({
    # ----------------------------------------------------------------------
    # Absorption (Blesch 2003 Table 1)
    # ----------------------------------------------------------------------
    # KA: first-order absorption into the 5'-DFUR pool. The dose is
    # capecitabine but no capecitabine or 5'-DFCR compartment is carried;
    # KA represents the lumped capecitabine -> 5'-DFCR -> 5'-DFUR formation
    # rate. ISV 70% CV. The source also reports a 70% CV inter-occasion
    # variability on KA estimated against the Phase II breast-cancer cycle
    # data (Reigner 1999); IOV is NOT represented in this static model
    # file (see vignette Errata).
    lka <- log(1.09)
    label("Lumped absorption rate constant KA into 5'-DFUR (1/h)")  # Blesch 2003 Table 1: TV KA = 1.09 1/h (SE 0.166)
    # TLAG: absorption lag time. The published typical value is 5.52E-4 h
    # (essentially zero) with a 49,498% CV inter-subject variability that
    # is effectively unidentified -- see vignette Errata. The IOV on
    # TLAG (52,915% CV) is similarly extreme and is NOT represented.
    ltlag <- log(5.52e-4)
    label("Capecitabine absorption lag TLAG (h)")  # Blesch 2003 Table 1: TV TLAG = 5.52E-4 h (SE 2.02E-4)

    # ----------------------------------------------------------------------
    # 5'-DFUR (compartment central_dfur, parameters V1/CL1)
    # ----------------------------------------------------------------------
    # V1: apparent 5'-DFUR volume V1/F. ISV on V1 was fixed at 30% CV in
    # the source model (Table 1 footnote "ISV fixed").
    lvc_dfur <- log(90.6)
    label("Apparent 5'-DFUR volume V1/F (L)")  # Blesch 2003 Table 1: TV V1 = 90.6 L (SE 14.1)
    lcl_dfur <- log(75.8)
    label("Apparent 5'-DFUR clearance CL1/F (L/h)")  # Blesch 2003 Table 1: TV CL1 = 75.8 L/h (SE 1.8)

    # ----------------------------------------------------------------------
    # 5-FU (compartment central_5fu, parameters V2/CL2)
    # ----------------------------------------------------------------------
    # V2: fixed to 17.8 L because sensitivity analysis showed V2 was not a
    # sensitive parameter for the model (Blesch 2003 page describing model
    # development; literature value from Heggie 1987 [reference 2 in
    # source]). No ISV estimated for V2.
    lvc_5fu <- fixed(log(17.8))
    label("Apparent 5-FU volume V2/F (L) -- FIXED to Heggie 1987 literature value")  # Blesch 2003 Table 1: V2 fixed at 17.8 L; "TV fixed, no ISV"
    lcl_5fu <- log(1190)
    label("Apparent 5-FU clearance CL2/F at ALP = 100 U/L (L/h)")  # Blesch 2003 Table 1: TV CL2 = 1190 L/h (SE 39.3)

    # ----------------------------------------------------------------------
    # FBAL (compartment central_fbal, parameters V3/CL3)
    # ----------------------------------------------------------------------
    lvc_fbal <- log(73.6)
    label("Apparent FBAL volume V3/F at CRCL = 80 mL/min/1.73 m^2 and BSA = 1.73 m^2 (L)")  # Blesch 2003 Table 1: TV V3 = 73.6 L (SE 2.35)
    lcl_fbal <- log(27.5)
    label("Apparent FBAL clearance CL3/F at CRCL = 80 mL/min/1.73 m^2 (L/h)")  # Blesch 2003 Table 1: TV CL3 = 27.5 L/h (SE 0.864)

    # ----------------------------------------------------------------------
    # Covariate effects (Blesch 2003 Table 1 / Table 2; power form)
    # Reference (centering) values for ALP, CRCL, and BSA are NOT stated
    # in Blesch 2003; clinically reasonable medians for an advanced-cancer
    # cohort are used here as non-paper provenance (vignette Errata).
    # ----------------------------------------------------------------------
    e_alp_cl_5fu  <- -0.169
    label("Power exponent of (ALP / 100) on apparent 5-FU clearance CL2/F (unitless)")  # Blesch 2003 Table 1: APHCL2 = -0.169 (SE 0.0586)
    e_crcl_cl_fbal <- 0.615
    label("Power exponent of (CRCL / 80) on apparent FBAL clearance CL3/F (unitless)")  # Blesch 2003 Table 1: CLRCL3 = 0.615 (SE 0.0769)
    e_crcl_vc_fbal <- 0.394
    label("Power exponent of (CRCL / 80) on apparent FBAL volume V3/F (unitless)")  # Blesch 2003 Table 1: CLRV3 = 0.394 (SE 0.109)
    e_bsa_vc_fbal <- 0.812
    label("Power exponent of (BSA / 1.73) on apparent FBAL volume V3/F (unitless)")  # Blesch 2003 Table 1: BSAV3 = 0.812 (SE 0.189)

    # ----------------------------------------------------------------------
    # Inter-individual variability (ISV; %CV in source -> log-normal omega^2
    # via omega^2 = log(1 + CV^2)). Inter-occasion variability on KA and
    # TLAG estimated in the source paper is NOT represented; see vignette
    # Errata. The 70% CV ISV on KA is reported alongside a 70% CV IOV on
    # KA; the 49,498% CV ISV on TLAG is reported alongside a 52,915% CV
    # IOV on TLAG -- the TLAG IIVs are effectively unidentified parameters
    # (see vignette Errata for the simulation implications).
    # ----------------------------------------------------------------------
    etalka       ~ 0.39878  # Blesch 2003 Table 1: ISV KA 70% CV -> log(1 + 0.70^2)
    etaltlag     ~ 12.41    # Blesch 2003 Table 1: ISV TLAG 49,498% CV -> log(1 + 494.98^2); effectively unidentified
    etalvc_dfur  ~ fixed(0.08618)  # Blesch 2003 Table 1: ISV V1 30% CV (FIXED) -> log(1 + 0.30^2)
    etalcl_dfur  ~ 0.05598  # Blesch 2003 Table 1: ISV CL1 24% CV -> log(1 + 0.24^2)
    etalcl_5fu   ~ 0.10336  # Blesch 2003 Table 1: ISV CL2 33% CV -> log(1 + 0.33^2)
    etalvc_fbal  ~ 0.06537  # Blesch 2003 Table 1: ISV V3 26% CV -> log(1 + 0.26^2)
    etalcl_fbal  ~ 0.09751  # Blesch 2003 Table 1: ISV CL3 32% CV -> log(1 + 0.32^2)

    # ----------------------------------------------------------------------
    # Residual error. Source reports proportional (%CV) residual variances
    # separately for extensive-sampling and sparse-sampling data; the values
    # here use the EXTENSIVE-sampling estimates (the typical-richness case).
    # The 5'-DFUR / 5-FU residual correlation of 0.77 (SE 0.0435) is NOT
    # encoded -- nlmixr2's residual-error syntax does not directly support
    # cross-output residual covariances; see vignette Errata. The FBAL
    # residual variance was estimated via sigma_3^2 = (tau * sigma_1)^2 with
    # tau = 0.767 (SE 0.121) fixed-effects parameter (Blesch 2003 Table 1
    # footnote); the resulting 34% CV is reported directly in Table 1 and
    # encoded as-is here.
    # ----------------------------------------------------------------------
    propSd_dfur <- 0.39
    label("Proportional residual SD for 5'-DFUR (extensive sampling)")  # Blesch 2003 Table 1: Res. Error 5'-DFUR 39% extensive / 51% sparse
    propSd_5fu  <- 0.58
    label("Proportional residual SD for 5-FU (extensive sampling)")    # Blesch 2003 Table 1: Res. Error 5-FU 58% extensive / 76% sparse
    propSd_fbal <- 0.34
    label("Proportional residual SD for FBAL (extensive sampling)")    # Blesch 2003 Table 1: Res. Error FBAL 34% extensive / 44% sparse (derived from tau = 0.767 x sigma_DFUR)
  })

  model({
    # 1. Individual PK parameters. KA, TLAG, V1, V2, CL1, CL2, V3, CL3
    # are exponentiated log-parameters with their respective etas. Power
    # covariates on CL2 (ALP), CL3 (CRCL), V3 (CRCL and BSA).
    ka       <- exp(lka + etalka)
    tlag     <- exp(ltlag + etaltlag)
    vc_dfur  <- exp(lvc_dfur + etalvc_dfur)
    cl_dfur  <- exp(lcl_dfur + etalcl_dfur)
    vc_5fu   <- exp(lvc_5fu)
    cl_5fu   <- exp(lcl_5fu + etalcl_5fu) * (ALP / 100)^e_alp_cl_5fu
    vc_fbal  <- exp(lvc_fbal + etalvc_fbal) *
                  (CRCL / 80)^e_crcl_vc_fbal *
                  (BSA / 1.73)^e_bsa_vc_fbal
    cl_fbal  <- exp(lcl_fbal + etalcl_fbal) *
                  (CRCL / 80)^e_crcl_cl_fbal

    # 2. Sequential metabolite micro-constants (each metabolite is the
    # input flux of the next). Apparent-clearance convention absorbs both
    # the molar conversion fraction and the prior-step bioavailability.
    kel_dfur <- cl_dfur / vc_dfur
    kel_5fu  <- cl_5fu  / vc_5fu
    kel_fbal <- cl_fbal / vc_fbal

    # 3. ODE system (Blesch 2003 page describing the differential
    # equations and Figure 6 structural diagram). Mass enters the 5'-DFUR
    # compartment from depot via KA; 5'-DFUR is sequentially metabolised
    # to 5-FU, then 5-FU to FBAL, then FBAL is eliminated.
    d/dt(depot)        <- -ka * depot
    d/dt(central_dfur) <-  ka       * depot         - kel_dfur * central_dfur
    d/dt(central_5fu)  <-  kel_dfur * central_dfur  - kel_5fu  * central_5fu
    d/dt(central_fbal) <-  kel_5fu  * central_5fu   - kel_fbal * central_fbal

    # 4. Absorption lag time.
    lag(depot) <- tlag

    # 5. Plasma concentration outputs.
    Cc_dfur <- central_dfur / vc_dfur
    Cc_5fu  <- central_5fu  / vc_5fu
    Cc_fbal <- central_fbal / vc_fbal

    Cc_dfur ~ prop(propSd_dfur)
    Cc_5fu  ~ prop(propSd_5fu)
    Cc_fbal ~ prop(propSd_fbal)
  })
}
