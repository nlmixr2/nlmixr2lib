Zheng_2020_CNTO5048_mouse_mpbpk <- function() {
  description <- "Preclinical (mouse, CD45RB-high T-cell transfer colitis). mPBPK. Minimal PBPK model with quasi-equilibrium TMDD in serum and in a two-segment colon interstitium for CNTO 5048, an anti-murine TNF surrogate mAb of golimumab, and its soluble TNF target in IBD and non-IBD SCID mice (Zheng 2020 mPBPK/TE model)"
  reference <- "Zheng S, Niu J, Geist B, Fink D, Xu Z, Zhou H, Wang W. A minimal physiologically based pharmacokinetic model to characterize colon TNF suppression and treatment effects of an anti-TNF monoclonal antibody in a mouse inflammatory bowel disease model. mAbs. 2020;12(1):1813962. doi:10.1080/19420862.2020.1813962"
  vignette <- "Zheng_2020_CNTO5048_mouse_mpbpk"
  units <- list(
    time = "h",
    dosing = "pmol (convert mg via dose_pmol = dose_mg * 1e9 / 150000; 150 kDa IgG)",
    concentration = "nM"
  )

  # The two colon interstitial-fluid segments and their total-TNF pools are
  # paper-anatomical states of this model (Zheng 2020 Figure 3(b)-(d)). They
  # are deliberately not the luminal GI-transit `colon` compartment.
  paper_specific_compartments <- c(
    "colon1",
    "colon2",
    "total_target_colon1",
    "total_target_colon2"
  )

  compartmentData <- list(
    depot = list(analyte = "CNTO 5048", units = "pmol", specimen = "administration site", verified = TRUE),
    plasma = list(analyte = "CNTO 5048 (total)", units = "pmol", specimen = "serum", verified = TRUE),
    tight = list(analyte = "CNTO 5048", units = "pmol", specimen = "tissue", verified = TRUE),
    leaky = list(analyte = "CNTO 5048", units = "pmol", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "CNTO 5048", units = "pmol", specimen = "lymph", verified = TRUE),
    colon1 = list(analyte = "CNTO 5048 (total)", units = "pmol", specimen = "tissue", verified = TRUE),
    colon2 = list(analyte = "CNTO 5048 (total)", units = "pmol", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "TNF (total: free + CNTO 5048-bound)",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    total_target_colon1 = list(
      analyte = "TNF (total: free + CNTO 5048-bound)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    total_target_colon2 = list(
      analyte = "TNF (total: free + CNTO 5048-bound)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    DIS_TCT_COLITIS = list(
      description = "T-cell-transfer colitis (mouse IBD model) indicator: 1 = IBD mouse, 0 = non-IBD control mouse",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-IBD SCID control mouse that did not receive the CD45RB-high T-cell transfer)",
      notes = paste(
        "Zheng 2020 Table 2 estimates Vmax, sigma_colon and CL_colon separately for non-IBD and",
        "IBD mice, and fixes V_colon and L_colon at separate physiological values for each (eq. 6",
        "and 7). The leaky-tissue ISF volume and lymph flow are reduced by the colon values, so",
        "they also switch with disease status. Soluble TNF was below the assay limit in every",
        "non-IBD mouse, and TMDD was evaluated only in IBD mice ('TMDD was only evaluated for the",
        "IBD animals', Discussion); the model therefore sets the serum and colon TNF baselines to",
        "zero when DIS_TCT_COLITIS = 0, which makes every binding term vanish."
      ),
      source_name = "IBD vs non-IBD group (Table 1)"
    )
  )

  population <- list(
    species = "mouse (female SCID, Fox Chase C.B-17)",
    n_subjects = 144L,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "about 20 g (mean body weight used by the authors for per-kg conversions; Discussion)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state = paste(
      "114 SCID mice given an intraperitoneal injection of CD45RB-high T cells from female Balb/C",
      "donors on study day 0 to induce colitis ('IBD mice', Groups 2-6) and 30 non-IBD SCID controls",
      "(Group 1) that did not receive the T-cell transfer."
    ),
    dose_range = paste(
      "First dose on study day 21 (model time 0). Group 1 (non-IBD) and Group 3: single 10 mg/kg IV;",
      "Group 4: 10 mg/kg IV + 9 x 0.3 mg/kg IP every 3 days (study days 24-48); Group 5: 1.4 mg/kg IV +",
      "9 x 1.4 mg/kg IP Q3D; Group 6: 0.3 mg/kg IV + 9 x 0.3 mg/kg IP Q3D; Group 2: isotype control",
      "CNTO 1322 10 mg/kg IV + 9 x 0.3 mg/kg IP Q3D (Table 1)."
    ),
    regions = "Bolder BioPATH, Boulder, CO, USA",
    notes = paste(
      "Sparse destructive sampling (serum by retro-orbital or cardiac puncture, colon homogenate at",
      "necropsy). Naive-pooled fit in Monolix 2019R1 with the omega matrix fixed to zero, so the model",
      "has no between-animal variability. Serum CNTO 5048 had a proportional residual error and the",
      "colon CNTO 5048, serum TNF and colon TNF outputs additive residual errors, but the paper does not",
      "report their magnitudes."
    )
  )

  ini({
    # ---- Physiological constants (Table 2, footnote a: from Cao 2013, Chen 2018, Shah 2012) ----
    lvc <- fixed(log(0.85)); label("Serum volume Vs (mL)") # Table 2: Vs = 0.85 mL (fixed)
    lvisf <- fixed(log(4.35)); label("Total interstitial fluid volume ISF (mL)") # Table 2: ISF = 4.35 mL (fixed)
    lvlymph <- fixed(log(1.60)); label("Lymph volume Vlymph (mL)") # Table 2: Vlymph = 1.60 mL (fixed)
    llymphflow <- fixed(log(0.12)); label("Total lymph flow rate L (mL/h)") # Table 2: L = 0.12 mL/hr (fixed)
    sigma_l <- fixed(0.20); label("Lymphatic capillary reflection coefficient sigma_L (unitless)") # Table 2: sigma_L = 0.20 (fixed)
    kp <- fixed(0.8); label("Available fraction of ISF for IgG distribution Kp (unitless)") # Methods Step I text: Kp = .8

    # ---- Serum mAb elimination and distribution (Table 2) ----
    lvmax_healthy <- log(3.39); label("Serum mAb elimination capacity Vmax, non-IBD mice (pmol/h)") # Table 2: Vmax,nonIBD = 3.39 pmol/hr (RSE 13.5%)
    lvmax_ibd <- log(4.53); label("Serum mAb elimination capacity Vmax, IBD mice (pmol/h)") # Table 2: Vmax,IBD = 4.53 pmol/hr (RSE 43.7%)
    lkm <- log(221); label("Apparent affinity for mAb elimination Km (nM)") # Table 2: Km = 221 nM (RSE 18.9%)
    sigma_tight <- 0.955; label("Vascular reflection coefficient for tight tissues (unitless)") # Table 2: sigma_tight = 0.955 (RSE 1.37%)
    sigma_leaky <- 0.201; label("Vascular reflection coefficient for leaky tissues (unitless)") # Table 2: sigma_leaky = 0.201 (RSE 21.5%)
    lka <- log(0.104); label("First-order intraperitoneal absorption rate constant ka (1/h)") # Table 2: ka = 0.104 hr-1 (RSE 71.8%)

    # ---- Colon distribution (Table 2; eq. 6-7 for the fixed physiological values) ----
    sigma_colon_healthy <- 0.978; label("Reflection coefficient for colon, non-IBD mice (unitless)") # Table 2: sigma_colon,nonIBD = 0.978 (RSE 1.57%)
    sigma_colon_ibd <- 0.391; label("Reflection coefficient for colon, IBD mice (unitless)") # Table 2: sigma_colon,IBD = 0.391 (RSE 16.9%)
    lcl_colon_healthy <- log(0.000106); label("Clearance of mAb from colon ISF, non-IBD mice (mL/h)") # Table 2: CL_colon,nonIBD = 0.000106 mL/hr (RSE 59.6%)
    lcl_colon_ibd <- log(0.00217); label("Clearance of mAb from colon ISF, IBD mice (mL/h)") # Table 2: CL_colon,IBD = 0.00217 mL/hr (RSE 20.1%)
    lvcolon_healthy <- fixed(log(0.026)); label("Colon ISF volume, non-IBD mice (mL)") # Table 2 / eq. 6: V_colon,nonIBD = 0.026 mL (fixed)
    lvcolon_ibd <- fixed(log(0.039)); label("Colon ISF volume, IBD mice (mL)") # Table 2 / eq. 6: V_colon,IBD = 0.039 mL (fixed)
    llymphflow_colon_healthy <- fixed(log(0.00137)); label("Colon lymph flow rate, non-IBD mice (mL/h)") # Table 2 / eq. 7: L_colon,nonIBD = 0.00137 mL/hr (fixed)
    llymphflow_colon_ibd <- fixed(log(0.00205)); label("Colon lymph flow rate, IBD mice (mL/h)") # Table 2 / eq. 7: L_colon,IBD = 0.00205 mL/hr (fixed)
    ratio_isf <- fixed(0.175); label("Colon ISF volume per total colon tissue weight (mL/g)") # Methods Step II: ISF-to-tissue ratio in mouse large intestine = 17.5% (Shah 2012)
    ratio_serum <- fixed(0.0159); label("Residual serum volume per total colon tissue weight (mL/g)") # Methods Step II: residual serum-to-tissue ratio = 1.59% (Shah 2012)

    # ---- Serum TNF target engagement (Table 2; eq. 16-21) ----
    lkdeg <- log(12.8); label("Degradation rate constant of free TNF in serum kdeg (1/h)") # Table 2: kdeg = 12.8 hr-1 (RSE 17.9%)
    lkint <- log(0.866); label("Elimination rate constant of CNTO 5048-TNF complex in serum kint (1/h)") # Table 2: kint = 0.866 hr-1 (RSE 0.554%)
    lr0 <- fixed(log(0.00222)); label("Baseline TNF concentration in serum of IBD mice R0 (nM)") # Table 2: R0 = 2.22 pM (fixed) = 0.00222 nM
    lkss <- log(1.88); label("Quasi-equilibrium binding constant of CNTO 5048 and TNF in serum Kss (nM)") # Table 2: Kss = 1.88 nM (RSE 34.1%)

    # ---- Colon TNF target engagement (Table 2; eq. 22-35) ----
    lkdeg_colon <- log(1.17); label("Degradation rate constant of free TNF in colon ISF (1/h)") # Table 2: kdeg,colon = 1.17 hr-1 (RSE 33.1%)
    lkint_colon <- log(0.0492); label("Elimination rate constant of CNTO 5048-TNF complex in colon ISF (1/h)") # Table 2: kint,colon = 0.0492 hr-1 (RSE 50.8%)
    lr0_colon <- fixed(log(0.420)); label("Baseline TNF concentration in colon ISF of IBD mice (nM)") # Table 2: R0,colon = 420 pM (fixed) = 0.420 nM
    lkss_colon <- log(1.32); label("Quasi-equilibrium binding constant of CNTO 5048 and TNF in colon ISF (nM)") # Table 2: Kss,colon = 1.32 nM (RSE 0.09%)

    # ---- Residual error (Methods, Model fitting): forms stated, magnitudes not reported ----
    propSd <- fixed(0); label("Proportional residual error, serum total CNTO 5048 (fraction; not reported)") # Methods: proportional error for serum CNTO 5048; magnitude not reported
    addSd_Ccolon <- fixed(0); label("Additive residual error, colon homogenate CNTO 5048 (nM; not reported)") # Methods: constant error for colon CNTO 5048; magnitude not reported
    addSd_freeTnf <- fixed(0); label("Additive residual error, serum free TNF (nM; not reported)") # Methods: constant error for serum TNF; magnitude not reported
    addSd_freeTnf_colon <- fixed(0); label("Additive residual error, colon homogenate free TNF (nM; not reported)") # Methods: constant error for colon TNF; magnitude not reported
  })

  model({
    # Disease-status switches (Table 2 non-IBD / IBD parameter pairs)
    vmax <- DIS_TCT_COLITIS * exp(lvmax_ibd) + (1 - DIS_TCT_COLITIS) * exp(lvmax_healthy)
    sigma_colon <- DIS_TCT_COLITIS * sigma_colon_ibd + (1 - DIS_TCT_COLITIS) * sigma_colon_healthy
    cl_colon <- DIS_TCT_COLITIS * exp(lcl_colon_ibd) + (1 - DIS_TCT_COLITIS) * exp(lcl_colon_healthy)
    vcolon <- DIS_TCT_COLITIS * exp(lvcolon_ibd) + (1 - DIS_TCT_COLITIS) * exp(lvcolon_healthy)
    lymphflow_colon <- DIS_TCT_COLITIS * exp(llymphflow_colon_ibd) + (1 - DIS_TCT_COLITIS) * exp(llymphflow_colon_healthy)

    vc <- exp(lvc)
    visf <- exp(lvisf)
    vlymph <- exp(lvlymph)
    lymphflow <- exp(llymphflow)
    km <- exp(lkm)
    ka <- exp(lka)

    # Tissue volumes and lymph flows (Methods Step I-II): the colon is carved
    # out of the leaky tissue, Vleaky = 0.35*ISF*Kp - Vcolon, Lleaky = 2/3*L - Lcolon.
    vtight <- 0.65 * visf * kp
    vleaky <- 0.35 * visf * kp - vcolon
    lymphflow_tight <- lymphflow / 3
    lymphflow_leaky <- 2 * lymphflow / 3 - lymphflow_colon
    vcolon_seg <- 0.5 * vcolon

    # TNF turnover (eq. 17 and 28); no measurable TNF in non-IBD mice
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    kss <- exp(lkss)
    r0 <- DIS_TCT_COLITIS * exp(lr0)
    ksyn <- kdeg * r0
    kdeg_colon <- exp(lkdeg_colon)
    kint_colon <- exp(lkint_colon)
    kss_colon <- exp(lkss_colon)
    r0_colon <- DIS_TCT_COLITIS * exp(lr0_colon)
    ksyn_colon <- kdeg_colon * r0_colon

    total_target(0) <- r0
    total_target_colon1(0) <- r0_colon
    total_target_colon2(0) <- r0_colon

    # Total drug concentrations (nM = pmol/mL)
    cs <- plasma / vc
    ctight <- tight / vtight
    cleaky <- leaky / vleaky
    clymph <- lymph / vlymph
    ccolon1 <- colon1 / vcolon_seg
    ccolon2 <- colon2 / vcolon_seg

    # Quasi-equilibrium free drug and complex in serum (eq. 18, 20, 21)
    ds <- cs - kss - total_target
    cfree_s <- 0.5 * (ds + sqrt(ds * ds + 4 * kss * cs))
    ar <- total_target * cfree_s / (cfree_s + kss)
    rfree_s <- total_target * kss / (cfree_s + kss)

    # Quasi-equilibrium free drug and complex in each colon segment (eq. 22-23, 30-31, 33-34)
    d1 <- ccolon1 - kss_colon - total_target_colon1
    cfree_colon1 <- 0.5 * (d1 + sqrt(d1 * d1 + 4 * kss_colon * ccolon1))
    ar_colon1 <- total_target_colon1 * cfree_colon1 / (cfree_colon1 + kss_colon)
    rfree_colon1 <- total_target_colon1 * kss_colon / (cfree_colon1 + kss_colon)
    d2 <- ccolon2 - kss_colon - total_target_colon2
    cfree_colon2 <- 0.5 * (d2 + sqrt(d2 * d2 + 4 * kss_colon * ccolon2))
    ar_colon2 <- total_target_colon2 * cfree_colon2 / (cfree_colon2 + kss_colon)
    rfree_colon2 <- total_target_colon2 * kss_colon / (cfree_colon2 + kss_colon)

    # ODEs in amount form (pmol): each printed concentration equation multiplied by its volume.
    # Eq. 15: intraperitoneal absorption site
    d/dt(depot) <- -ka * depot
    # Eq. 11: serum (Michaelis-Menten elimination of free drug with total drug in the denominator, as printed)
    d/dt(plasma) <- clymph * lymphflow -
      cfree_s * lymphflow_tight * (1 - sigma_tight) -
      cfree_s * lymphflow_leaky * (1 - sigma_leaky) -
      cfree_s * vmax / (km + cs) -
      kint * ar * vc
    # Eq. 12-13: tight and leaky tissue ISF
    d/dt(tight) <- cfree_s * lymphflow_tight * (1 - sigma_tight) - ctight * lymphflow_tight * (1 - sigma_l)
    d/dt(leaky) <- cfree_s * lymphflow_leaky * (1 - sigma_leaky) - cleaky * lymphflow_leaky * (1 - sigma_l)
    # Eq. 14: lymph
    d/dt(lymph) <- ctight * lymphflow_tight * (1 - sigma_l) +
      cleaky * lymphflow_leaky * (1 - sigma_l) -
      clymph * lymphflow +
      ka * depot
    # Eq. 24-25: two sequential colon ISF segments
    d/dt(colon1) <- cfree_s * lymphflow_colon * (1 - sigma_colon) -
      cfree_colon1 * lymphflow_colon * (1 - sigma_colon) -
      cfree_colon1 * cl_colon -
      kint_colon * ar_colon1 * vcolon_seg
    d/dt(colon2) <- cfree_colon1 * lymphflow_colon * (1 - sigma_colon) -
      cfree_colon2 * lymphflow_colon * (1 - sigma_l) -
      cfree_colon2 * cl_colon -
      kint_colon * ar_colon2 * vcolon_seg
    # Eq. 16, 26-27: total TNF (nM)
    d/dt(total_target) <- ksyn - kdeg * (total_target - ar) - kint * ar
    d/dt(total_target_colon1) <- ksyn_colon - kdeg_colon * (total_target_colon1 - ar_colon1) - kint_colon * ar_colon1
    d/dt(total_target_colon2) <- ksyn_colon - kdeg_colon * (total_target_colon2 - ar_colon2) - kint_colon * ar_colon2

    # Observations (nM). Colon values are per gram of wet colon tissue (density 1 g/mL),
    # i.e. homogenate results divided by the 0.2 g/mL homogenisation ratio (Methods).
    Cc <- cs
    Ccolon <- ratio_isf * (ccolon1 + ccolon2) / 2 + ratio_serum * cs # eq. 10 (total CNTO 5048)
    Cfree_colon <- ratio_isf * (cfree_colon1 + cfree_colon2) / 2 + ratio_serum * cfree_s # eq. 32
    freeTnf <- rfree_s # eq. 21
    freeTnf_colon <- ratio_isf * (rfree_colon1 + rfree_colon2) / 2 + ratio_serum * rfree_s # eq. 35

    Cc ~ prop(propSd)
    Ccolon ~ add(addSd_Ccolon)
    freeTnf ~ add(addSd_freeTnf)
    freeTnf_colon ~ add(addSd_freeTnf_colon)
  })
}
