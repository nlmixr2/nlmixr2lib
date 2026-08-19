Karakitsios_2025_bedaquiline_human_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2023R1).",
    "Human extrapolation of bedaquiline distribution into healthy/uninvolved and",
    "tuberculosis-infected lung tissue in patients with pulmonary tuberculosis. The published",
    "Svensson 2016 three-compartment plasma model with two-transit-compartment oral absorption",
    "drives a permeability-limited lung model split into pulmonary blood, extracellular water,",
    "intracellular water and lysosomes, with pH-dependent ion trapping. The infected-lung",
    "extension adds a cellular-lesion (macrophage) compartment and a caseous-granuloma arm in",
    "which an outer caseum of foamy macrophages is in instantaneous equilibrium with a",
    "six-compartment catenary chain representing 300-micrometre rings of the necrotic core, with",
    "anomalous (fractal) diffusion described by a time-dependent rate coefficient. All lung,",
    "lesion and caseum parameters were estimated in mice and carried over unchanged, with the",
    "lysosome:cytosol partition rescaled to the human intracellular pH of 6.69. The lung states",
    "hold CONCENTRATIONS, not amounts. No human lung data exist, so this layer is an unvalidated",
    "extrapolation.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Della Pasqua O, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of bedaquiline across species using",
    "physiologically-based pharmacokinetic modelling.",
    "Br J Clin Pharmacol. 2025;91(11):3167-3178.",
    "doi:10.1002/bcp.70163.",
    "The plasma layer is the bedaquiline (parent) component of",
    "Svensson E. M., Dosne A.-G., Karlsson M. O. (2016).",
    "Population Pharmacokinetics of Bedaquiline and Metabolite M2 in Patients With",
    "Drug-Resistant Tuberculosis. CPT Pharmacometrics Syst Pharmacol 5(12):682-691.",
    "doi:10.1002/psp4.12147; see modellib('Svensson_2016_bedaquiline').",
    sep = " "
  )
  vignette <- "Karakitsios_2025_bedaquiline_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # cas1-cas6 are the paper's spatial binning of the necrotic caseous core into
  # 300-micrometre rings measured from the outer caseum edge (Supporting
  # Information, "PBPK models for the infected lung"). There is no canonical
  # chain prefix for this role; declared here pending operator ratification of
  # a `caseum<n>` canonical.
  paper_specific_compartments <- c(
    "caseum1", "caseum2", "caseum3", "caseum4", "caseum5", "caseum6"
  )

  compartmentData <- list(
    depot       = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion      = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum1     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum2     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum3     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum4     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum5     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum6     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  # The four covariates belong to the inherited Svensson 2016 plasma layer.
  # Karakitsios 2025 simulated "a mean typical patient", i.e. WT = 70,
  # ALB = 4.04, AGE = 32, RACE_BLACK = 0, which makes every factor below unity.
  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying in the source model; supplied per observation)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power on the plasma-layer clearances and volumes around a 70 kg reference:",
        "`cl = cl_typ * (WT/70)^e_wt_cl` (estimated exponent 0.181) and",
        "`vc = vc_typ * (WT/70)^e_wt_vc` (exponent fixed to 1.0). Inherited unchanged from",
        "Svensson 2016 Table 3; see modellib('Svensson_2016_bedaquiline') for the full",
        "parent-plus-metabolite model and its time-varying weight trajectory. The LUNG layer",
        "does NOT scale with WT: the Table S3 human lung blood flow, sub-volumes and cell",
        "surface area are fixed reference-human quantities.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin concentration (time-varying in the source model)",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect around the typical population steady-state value Ass = 4.04 g/dL:",
        "`cl = cl_typ * (ALB/4.04)^e_alb_cl` with estimated exponent 1.64, and an additional",
        "fixed-coefficient unbound-fraction adjustment `(4.04/ALB)^e_alb_vc` applied to the",
        "apparent volumes in the observation equation only. Inherited unchanged from",
        "Svensson 2016 Table 3.",
        sep = " "
      ),
      source_name        = "ALB"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear deviation around the cohort median of 32 years applied to clearance:",
        "`cl = cl_typ * (1 + e_age_cl * (32 - AGE))` with estimated coefficient 0.00881 1/year.",
        "Inherited unchanged from Svensson 2016 Table 3.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    RACE_BLACK = list(
      description        = "Black-race indicator (1 = Black, 0 = non-Black)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; the reference category includes White, Asian, Hispanic and Other in the Svensson 2016 cohort).",
      notes              = paste(
        "Multiplicative effect on clearance: `cl = cl_typ * (1 + e_race_black_cl * RACE_BLACK)`",
        "with estimated coefficient 0.84, i.e. about 84% higher clearance in Black patients.",
        "Inherited unchanged from Svensson 2016 Table 3.",
        sep = " "
      ),
      source_name        = "RACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 335L,
    n_studies      = 2L,
    age_range      = "18-68 years (Svensson 2016 Table 1)",
    age_median     = "32 years (Svensson 2016 Table 1)",
    weight_range   = "30-113 kg (Svensson 2016 Table 1)",
    weight_median  = "55-57 kg (study-stage dependent; Svensson 2016 Table 1)",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Adults with pulmonary multidrug-resistant, pre-extensively or extensively drug-resistant",
      "tuberculosis. The lung layer represents TB-infected lung containing both cellular",
      "(non-necrotising) lesions and caseous necrotic granulomas.",
      sep = " "
    ),
    dose_range     = paste(
      "The approved regimen simulated in Methods 2.4: oral bedaquiline 400 mg once daily for a",
      "2-week loading phase, then 200 mg three times weekly through week 24.",
      sep = " "
    ),
    regions        = "Multicentre international (C208 and C209 phase II trials).",
    notes          = paste(
      "The plasma layer's population is that of Svensson 2016 (pooled C208 and C209). The lung",
      "layer has NO human data of any kind: Methods 2.1 states that human lung data were",
      "unavailable for validation, and the Discussion acknowledges that the lack of human lung",
      "data, particularly in TB-infected lung, prevents full model validation. Confidence rests",
      "entirely on the successful healthy-lung extrapolation to rats and dogs. Simulations",
      "assumed no pharmacokinetic drug interactions with companion drugs (Methods 2.4). Because",
      "of bedaquiline's very long half-life, steady state was not reached even after 14 weeks of",
      "treatment (Results 3.5).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer -- the bedaquiline (parent) component of the
    #    Svensson 2016 three-compartment model with a two-transit-compartment
    #    absorption chain, selected in the Supporting Information as the human
    #    plasma driver. Values inherited from Svensson 2016 Table 3 via
    #    modellib('Svensson_2016_bedaquiline'); the M2 metabolite is omitted
    #    because it does not enter the lung model. CL and V are apparent
    #    (F-relative) values, so F is fixed to 1.
    # =====================================================================
    lmat      <- log(0.6620 * 6)               ; label("Mean absorption time MAT (h, log scale)")             # inherited from Svensson 2016 Table 3: MAT = 0.66 as a fraction of 6 h
    logitfmat <- log(0.4664 / (1 - 0.4664))    ; label("Fraction of MAT in the transit delay FR (logit)")     # inherited from Svensson 2016 Table 3: FR = 0.47
    lcl       <- log(2.616)                    ; label("Apparent clearance CL/F (L/h, log scale)")            # inherited from Svensson 2016 Table 3: CL/F = 2.62 L/h
    lvc       <- log(198.34)                   ; label("Apparent central volume Vc/F (L, log scale)")         # inherited from Svensson 2016 Table 3: V/F = 198 L per 70 kg
    lq        <- log(3.658)                    ; label("Apparent inter-compartmental clearance Q1/F (L/h)")   # inherited from Svensson 2016 Table 3: Q1/F = 3.66 L/h
    lvp       <- log(8549.06)                  ; label("Apparent first peripheral volume Vp1/F (L)")          # inherited from Svensson 2016 Table 3: VP1/F = 8550 L per 70 kg
    lq2       <- log(7.335)                    ; label("Apparent inter-compartmental clearance Q2/F (L/h)")   # inherited from Svensson 2016 Table 3: Q2/F = 7.34 L/h
    lvp2      <- log(2690.91)                  ; label("Apparent second peripheral volume Vp2/F (L)")         # inherited from Svensson 2016 Table 3: VP2/F = 2690 L per 70 kg
    lfdepot   <- fixed(log(1))                 ; label("Bioavailability F (CL and V are apparent, F-relative values)")

    e_wt_cl         <- 0.1809      ; label("Allometric body-weight exponent on clearances (70 kg reference)")   # inherited from Svensson 2016 Table 3
    e_wt_vc         <- fixed(1.0)  ; label("Allometric body-weight exponent on volumes (70 kg reference)")      # inherited from Svensson 2016 Table 3, coefficient fixed to 1
    e_alb_vc        <- fixed(1.0)  ; label("Albumin unbound-fraction exponent on apparent volumes (Ass 4.04)")  # inherited from Svensson 2016 Table 3, fixed
    e_alb_cl        <- 1.640       ; label("Albumin power exponent on clearance (Ass = 4.04 g/dL)")             # inherited from Svensson 2016 Table 3
    e_race_black_cl <- 0.8387      ; label("Multiplicative Black-race effect on clearance")                     # inherited from Svensson 2016 Table 3
    e_age_cl        <- 0.008808    ; label("Linear age effect on clearance (1/year, around 32 years)")          # inherited from Svensson 2016 Table 3

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S2).
    # =====================================================================
    pka1 <- fixed(8.5) ; label("First (strong) basic pKa (unitless)")   # Table S2 pKa1 = 8.5 (base), Janssen R&D data
    pka2 <- fixed(1.1) ; label("Second (weak) basic pKa (unitless)")    # Table S2 pKa2 = 1.1 (base), ACD predicted

    # =====================================================================
    # 3. Species-related lung physiology for humans (Table S3).
    # =====================================================================
    qlung        <- fixed(356000) ; label("Lung blood flow rate Qlung (mL/h)")                   # Table S3 human blood flow rate = 356000 mL/h
    sa_cells     <- fixed(960000) ; label("Surface area of lung cells SAcells (cm^2)")           # Table S3 human surface area lung cells = 960000 cm^2 (Stone 1992)
    v_pb         <- fixed(101.2)  ; label("Volume of pulmonary blood VPB (mL)")                  # Table S3 human volume of pulmonary blood = 101.2 mL (Gill 2016; 80.7 kg human)
    v_ew         <- fixed(89.2)   ; label("Volume of lung extracellular water VEW (mL)")         # Table S3 human volume of EW = 89.2 mL
    v_lysosome   <- fixed(5.47)   ; label("Volume of lung lysosomes VLYS (mL)")                  # Table S3 human volume of lysosomes = 5.47 mL (1% of a 547 mL total lung volume)
    v_iw         <- fixed(351)    ; label("Volume of lung intracellular water VIW (mL)")         # Table S3 human volume of IW = 351 mL
    ph_plasma    <- fixed(7.4)    ; label("Plasma pH (unitless)")                                # Table S3 human pH_plasma = 7.4
    ph_ew        <- fixed(7.4)    ; label("Lung extracellular-water pH (unitless)")              # Table S3 human pH_EW = 7.4
    ph_iw        <- fixed(6.69)   ; label("Lung intracellular-water pH (unitless)")              # Table S3 human pH_IW = 6.69 (Simcyp 2023)
    ph_lysosome  <- fixed(5.0)    ; label("Lysosomal pH (unitless)")                             # Table S3 human pH_Lysosomes = 5.0; Methods 2.2 "constant pH of 5.0 across species"
    ph_lesion_iw <- fixed(5.84)   ; label("Cytosolic pH of macrophages and foamy macrophages")   # Methods 2.3: "A pH of 5.84 was used for foamy macrophages" (Kempker 2017); applied to both cell types per the shared-lesion-pH statement in Methods 2.3

    # =====================================================================
    # 4. Drug- AND species-related parameters (Table S4).
    # =====================================================================
    bp        <- fixed(1.0)    ; label("Blood-to-plasma concentration ratio B/P (unitless)")     # Table S4 human blood to plasma ratio = 1.0
    fu_plasma <- fixed(5.0e-4) ; label("Fraction unbound in plasma fupls (unitless)")            # Table S4 human fupls = 5.0E-4; the anchor value used to derive KaHSA in equations (S5)-(S6)
    fu_ew     <- fixed(8.2e-4) ; label("Fraction unbound in lung extracellular water fuEW")      # Table S4 human fuEW = 8.2E-4, from equation (S7)
    fu_iw     <- fixed(1.0)    ; label("Fraction unbound in lung intracellular water fuIW")      # Supporting Information: "a fixed value of 1 was used for fuIW"

    # =====================================================================
    # 5. Lung parameters carried over unchanged from the mouse fit
    #    (Tables S8, S11, S13) -- all FIXED here because none was estimated
    #    on human data. Equation (S9) rescales KpLYS:IW/fuIW to pH_IW = 6.69
    #    in model(), giving about 7.86e5 versus 3.48e6 in rodents.
    # =====================================================================
    lperm_cells      <- fixed(log(55.7))      ; label("Permeability of lung cells PermCells (cm/h, log scale)")   # Table S8 Permcells = 5.57E+1 cm/h, estimated in mice and held constant across species
    lfulys_fnilys    <- fixed(log(1.9411e-8)) ; label("Product fuLYS*fniLYS in lysosomes (unitless, log scale)")  # derived from Table S8 KpLYS:IW/fuIW = 3.48E+6 in mice via equations (S8)-(S9); "assumed to be equal in all" species
    lrate_macrophage <- fixed(log(143e3))     ; label("Macrophage permeation rate constant rate_mphi (1/h, log)") # Table S11 rate_mphi = 143E+3 1/h, fitted in BALB/c mice and kept constant across species
    lkp_macrophage   <- fixed(log(568e4))     ; label("Lysosome:cytosol partition in macrophages Kpmphi (log)")   # Table S11 Kpmphi/(fuIW*fniIW) = 568E+4, with fuIW = fniIW = 1
    lrate_foamy      <- fixed(log(5.78e3))    ; label("Foamy-macrophage permeation rate constant rate_fmphi (1/h, log)") # Table S13 rate_fmphi = 5.78E+3 1/h, fitted in Kramnik mice and kept constant across species
    lkp_foamy        <- fixed(log(1.44e6))    ; label("Lysosome:cytosol partition in foamy macrophages Kpfmphi (log)")   # Table S13 Kpfmphi/(fuIW*fniIW) = 1.44E+6, with fuIW = fniIW = 1
    lkp_caseum       <- fixed(log(0.71))      ; label("Caseum-1 : outer-caseum partition KpCaseum (unitless, log)")      # Table S13 KpCaseum = 0.71
    lkcas            <- fixed(log(0.25))      ; label("Anomalous-diffusion rate constant kcas (h^(h-1), log scale)")     # Table S13 kcas = 0.25
    h_caseum         <- fixed(0.75)           ; label("Heterogeneity exponent h of the fractal caseum diffusion")        # Table S13 h = 0.75; equation (S15) requires 0 <= h <= 1

    fu_iw_macrophage  <- fixed(1.0) ; label("Unbound fraction in macrophage cytosol fuIW(mphi)")     # Supporting Information: "a fixed value of 1 was used for fuIW(mphi), fniIW(mphi), fuIW(fmphi) and fniIW(fmphi)"
    fni_iw_macrophage <- fixed(1.0) ; label("Non-ionised fraction in macrophage cytosol fniIW(mphi)") # Supporting Information, same sentence

    # Equation (S15) tau: "fixed to a very small positive constant to avoid
    # infinities for t = 0"; the value is not reported. See vignette Errata.
    tau_caseum  <- fixed(1e-3)  ; label("Onset-time regulariser tau of equation (S15) (h)")
    flys_cell   <- fixed(0.078) ; label("Lysosomal contribution to cellular volume fLys(AM) (unitless)")  # Supporting Information: 7.8% in rat alveolar macrophages (Ufuk 2017), "assumed that the lysosomal contribution to the cellular volume in humans is the same as in rats"
    ratio_cas1  <- fixed(1.67)  ; label("Volume ratio outer caseum : caseum-1 (ratio1, unitless)")        # Supporting Information: "The same value of 1.67 was used for ratio1 in humans as well"
    ratio_inner <- fixed(1.0)   ; label("Volume ratio caseum(n-1) : caseum(n) for n = 2..6 (unitless)")   # Supporting Information: "a value equal to 1 was used for all other ratios, namely ratio2 to ratio6"

    # =====================================================================
    # 6. Random effects on the inherited plasma layer (Svensson 2016). The
    #    lung layer has no random effects: nothing in it was estimated on
    #    human data.
    # =====================================================================
    etalcl     ~ 0.1528  # inherited from Svensson 2016: BSV CL = 40.7 %CV (the CL/CLM2 correlation block is dropped with the M2 compartment)
    etalvc     ~ 0.1719  # inherited from Svensson 2016: BSV V = 43.3 %CV
    etalq      ~ 0.1812  # inherited from Svensson 2016: BSV Q1 = 44.5 %CV
    etalfdepot ~ 0.0803  # inherited from Svensson 2016: BSV F = 28.9 %CV
    etalmat    ~ 1.1620  # inherited from Svensson 2016: BOV MAT = 148 %CV, used here as the MAT random effect

    # =====================================================================
    # 7. Residual error. Only the plasma layer has one; the lung layer was
    #    never fitted to human data.
    # =====================================================================
    propSd <- sqrt(0.05182) ; label("Plasma proportional residual error (fraction)")  # inherited from Svensson 2016: proportional residual error BDQ = 23.1 %CV
  })

  model({
    # -------------------------------------------------------------------
    # 1. Covariate factors on the inherited plasma layer.
    # -------------------------------------------------------------------
    blackcl <- 1 + e_race_black_cl * RACE_BLACK
    agecl   <- 1 + e_age_cl * (32 - AGE)
    allcl   <- (WT / 70)^e_wt_cl
    allv    <- (WT / 70)^e_wt_vc
    albcl   <- (ALB / 4.04)^e_alb_cl
    albvfu  <- (4.04 / ALB)^e_alb_vc

    # -------------------------------------------------------------------
    # 2. Individual plasma parameters.
    # -------------------------------------------------------------------
    mat <- exp(lmat + etalmat)
    fr  <- 1 / (1 + exp(-logitfmat))
    mtt <- mat * fr
    ka  <- log(2) / (mat * (1 - fr) / 3.3)
    ktr <- 2 / mtt

    cl  <- exp(lcl + etalcl) * blackcl * agecl * albcl * allcl
    vc  <- exp(lvc + etalvc) * allv
    q   <- exp(lq  + etalq)  * allcl
    vp  <- exp(lvp)          * allv
    q2  <- exp(lq2)          * allcl
    vp2 <- exp(lvp2)         * allv

    # The albumin-driven free-fraction correction enters the OBSERVATION
    # volume only, matching the source model's $ERROR bookkeeping.
    vc_obs <- vc * albvfu

    perm_cells      <- exp(lperm_cells)
    rate_macrophage <- exp(lrate_macrophage)
    kp_macrophage   <- exp(lkp_macrophage)
    rate_foamy      <- exp(lrate_foamy)
    kp_foamy        <- exp(lkp_foamy)
    kp_caseum       <- exp(lkp_caseum)
    kcas            <- exp(lkcas)

    # -------------------------------------------------------------------
    # 3. Ionisation. Bedaquiline is a diprotic BASE (Table S2, Figure S8),
    #    so the four-term Henderson-Hasselbalch factor printed in the
    #    Supporting Information is applied in its product form
    #    (1+10^(pKa1-pH))*(1+10^(pKa2-pH)). See vignette Errata.
    # -------------------------------------------------------------------
    x_plasma <- (1 + 10^(pka1 - ph_plasma)) * (1 + 10^(pka2 - ph_plasma))
    x_ew     <- (1 + 10^(pka1 - ph_ew))     * (1 + 10^(pka2 - ph_ew))
    x_iw     <- (1 + 10^(pka1 - ph_iw))     * (1 + 10^(pka2 - ph_iw))

    fi_plasma <- 1 / x_plasma
    fi_ew     <- 1 / x_ew
    fi_iw     <- 1 / x_iw

    kp_lysosome <- fi_iw / exp(lfulys_fnilys)   # equation (S9)

    # -------------------------------------------------------------------
    # 4. Derived lung constants.
    # -------------------------------------------------------------------
    kewb     <- (fu_plasma / fu_ew) * (x_ew / x_plasma)
    ps       <- perm_cells * sa_cells
    v_ew_eff <- (v_ew + v_pb / kewb) / fu_ew
    v_iw_eff <- v_iw + kp_lysosome * v_lysosome
    f_iw     <- v_iw / (v_pb + v_ew + v_iw + v_lysosome)
    r_lys    <- flys_cell / f_iw

    cblood <- (central / vc_obs) * bp
    jperm  <- ps * (lung_ew * fi_ew - fu_iw * lung_iw * fi_iw)

    # -------------------------------------------------------------------
    # 5. Plasma ODEs -- two-transit-compartment absorption chain.
    # -------------------------------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ka * transit2
    d/dt(central)     <-  ka * transit2 - cl * central / vc -
                          q * central / vc + q * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    f(depot) <- exp(lfdepot + etalfdepot)

    # -------------------------------------------------------------------
    # 6. Healthy/uninvolved lung -- equations (S1) and (S4).
    # -------------------------------------------------------------------
    d/dt(lung_ew) <- (qlung * cblood - (qlung / fu_ew) * (lung_ew / kewb) - jperm) / v_ew_eff
    d/dt(lung_iw) <- jperm / v_iw_eff

    cu_ew_ni <- lung_ew * fi_ew

    # -------------------------------------------------------------------
    # 7. Infected lung -- cellular lesion, equation (S11).
    # -------------------------------------------------------------------
    denom_macrophage <- 1 / (r_lys + 1) + kp_macrophage * r_lys / (r_lys + 1)
    d/dt(lesion) <- rate_macrophage *
      (cu_ew_ni - fu_iw_macrophage * fni_iw_macrophage / denom_macrophage * lesion)

    # -------------------------------------------------------------------
    # 8. Infected lung -- caseous granuloma, equations (S12)-(S15).
    # -------------------------------------------------------------------
    kcaseum   <- kcas / (t + tau_caseum)^h_caseum
    x_foamy   <- (r_lys + 1) / (1 + kp_foamy * r_lys)
    vfac_cas1 <- 1 + ratio_cas1 / kp_caseum

    d/dt(caseum1) <- (rate_foamy * ratio_cas1 *
                        (cu_ew_ni - x_foamy * caseum1 / kp_caseum) -
                      kcaseum * caseum1 + kcaseum * caseum2 / ratio_inner) / vfac_cas1
    d/dt(caseum2) <- kcaseum * caseum1 * ratio_inner - 2 * kcaseum * caseum2 +
                     kcaseum * caseum3 / ratio_inner
    d/dt(caseum3) <- kcaseum * caseum2 * ratio_inner - 2 * kcaseum * caseum3 +
                     kcaseum * caseum4 / ratio_inner
    d/dt(caseum4) <- kcaseum * caseum3 * ratio_inner - 2 * kcaseum * caseum4 +
                     kcaseum * caseum5 / ratio_inner
    d/dt(caseum5) <- kcaseum * caseum4 * ratio_inner - 2 * kcaseum * caseum5 +
                     kcaseum * caseum6 / ratio_inner
    d/dt(caseum6) <- kcaseum * caseum5 * ratio_inner - kcaseum * caseum6

    # -------------------------------------------------------------------
    # 9. Observations.
    # -------------------------------------------------------------------
    Cc <- central / vc_obs

    c_lung_ew  <- lung_ew / fu_ew
    c_lung_pb  <- c_lung_ew / kewb
    c_lysosome <- kp_lysosome * lung_iw
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew +
               lung_iw * v_iw + c_lysosome * v_lysosome) /
              (v_pb + v_ew + v_iw + v_lysosome)

    c_lesion       <- lesion
    c_outer_caseum <- caseum1 / kp_caseum

    # Equations (S16)-(S19): total lysosomal concentration, then the unbound
    # cytosolic concentration, then the unbound lysosomal concentration.
    #
    # Tables S11 and S13 report the RATIO Kp/(fuIW*fniIW), estimated with
    # fuIW and fniIW fixed to 1, so `kp_macrophage` and `kp_foamy` above are
    # that ratio and not Kp itself. Equations (S16)-(S19) are statements about
    # the true Kp = CLYS/CIW, so the ratio is un-normalised here with the
    # prediction-time fniIW evaluated at the cellular-lesion pH of 5.84
    # (Methods 2.3). Skipping this step leaves the unbound lysosomal
    # concentrations about 460-fold low and contradicts the paper's own
    # conclusion that they exceed the MIC50; see vignette Errata. The ODE
    # states are unaffected -- equation (S11) is written and was fitted in
    # terms of the ratio, and is left exactly as published.
    x_lesion_iw <- (1 + 10^(pka1 - ph_lesion_iw)) * (1 + 10^(pka2 - ph_lesion_iw))
    x_lysosome  <- (1 + 10^(pka1 - ph_lysosome))  * (1 + 10^(pka2 - ph_lysosome))
    fni_lesion_iw <- 1 / x_lesion_iw
    fni_lysosome  <- 1 / x_lysosome

    kp_macrophage_true <- kp_macrophage * fu_iw_macrophage * fni_lesion_iw
    kp_foamy_true      <- kp_foamy      * fu_iw_macrophage * fni_lesion_iw

    c_lesion_lys  <- c_lesion * (r_lys + 1) / (1 / kp_macrophage_true + r_lys)   # (S17)
    cu_lesion_iw  <- fu_iw_macrophage * c_lesion_lys / kp_macrophage_true        # (S18)
    cu_lesion_lys <- cu_lesion_iw * (fni_lesion_iw / fni_lysosome)               # (S19)

    c_caseum_lys  <- c_outer_caseum * (r_lys + 1) / (1 / kp_foamy_true + r_lys)  # (S17)
    cu_caseum_iw  <- fu_iw_macrophage * c_caseum_lys / kp_foamy_true             # (S18)
    cu_caseum_lys <- cu_caseum_iw * (fni_lesion_iw / fni_lysosome)               # (S19)

    cu_plasma <- Cc * fu_plasma

    Cc ~ prop(propSd)
  })
}
