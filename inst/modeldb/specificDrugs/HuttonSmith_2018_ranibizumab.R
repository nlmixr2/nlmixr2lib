HuttonSmith_2018_ranibizumab <- function() {
  description <- "QSP. Three-compartment (retina, vitreous, aqueous) ocular PK/PD model of intravitreal ranibizumab and free VEGF suppression in wet AMD, with bivalent VEGF-ranibizumab binding (V, R, VR, RVR species per compartment; 12 ODEs). Retina-vitreous transfer uses hydrodynamic-radius-dependent ILM and RPE permeabilities (Hutton-Smith 2017 rabbit-PK power laws); vitreous-aqueous elimination k_el is derived from the ranibizumab ocular half-life via the mass-balance of Section S2. Parameters K_D = 19000 pM and log-normal population distributions of ocular half-life and VEGF production rate were re-estimated on 31 wet-AMD patients from Saunders 2015."
  reference <- "Hutton-Smith LA, Gaffney EA, Byrne HM, Caruso A, Maini PK, Mazer NA. Theoretical Insights into the Retinal Dynamics of Vascular Endothelial Growth Factor in Patients Treated with Ranibizumab, Based on an Ocular Pharmacokinetic/Pharmacodynamic Model. Mol Pharm. 2018;15(7):2770-2784. doi:10.1021/acs.molpharmaceut.8b00280. PMID: 29734810."
  vignette <- "HuttonSmith_2018_ranibizumab"

  paper_specific_compartments <- c(
    "vegf_ret",  "ranib_ret",  "vr_ret",  "rvr_ret",
    "vegf_vit",  "ranib_vit",  "vr_vit",  "rvr_vit",
    "vegf_aq",   "ranib_aq",   "vr_aq",   "rvr_aq"
  )

  units <- list(time = "day", dosing = "mg", concentration = "pM (pmol/L)")

  # Issue #482: the ocular matrices are the reason this model is not
  # comparable with a serum-clearance model. Verified against Hutton-Smith
  # 2018 Figure 1 (three-compartment eye: retina, vitreous, aqueous humour).
  compartmentData <- list(
    vegf_ret  = list(analyte = "VEGF (free)",               units = "pmol", specimen = "retina",         verified = TRUE),
    ranib_ret = list(analyte = "ranibizumab (free)",        units = "pmol", specimen = "retina",         verified = TRUE),
    vr_ret    = list(analyte = "VEGF-ranibizumab complex",  units = "pmol", specimen = "retina",         verified = TRUE),
    rvr_ret   = list(analyte = "ranibizumab-VEGF-ranibizumab complex", units = "pmol", specimen = "retina", verified = TRUE),
    vegf_vit  = list(analyte = "VEGF (free)",               units = "pmol", specimen = "vitreous",       verified = TRUE),
    ranib_vit = list(analyte = "ranibizumab (free)",        units = "pmol", specimen = "vitreous",       verified = TRUE),
    vr_vit    = list(analyte = "VEGF-ranibizumab complex",  units = "pmol", specimen = "vitreous",       verified = TRUE),
    rvr_vit   = list(analyte = "ranibizumab-VEGF-ranibizumab complex", units = "pmol", specimen = "vitreous", verified = TRUE),
    vegf_aq   = list(analyte = "VEGF (free)",               units = "pmol", specimen = "aqueous humour", verified = TRUE),
    ranib_aq  = list(analyte = "ranibizumab (free)",        units = "pmol", specimen = "aqueous humour", verified = TRUE),
    vr_aq     = list(analyte = "VEGF-ranibizumab complex",  units = "pmol", specimen = "aqueous humour", verified = TRUE),
    rvr_aq    = list(analyte = "ranibizumab-VEGF-ranibizumab complex", units = "pmol", specimen = "aqueous humour", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human (with rabbit-derived retinal permeabilities carried over from Hutton-Smith 2017)",
    n_subjects     = 31L,
    n_studies      = 1L,
    age_range      = NULL,
    weight_range   = NULL,
    sex_female_pct = NA_real_,
    disease_state  = "Neovascular (wet) age-related macular degeneration (wet AMD).",
    dose_range     = "Single intravitreal injection of 0.5 mg ranibizumab (the sensitivity analyses in the paper additionally simulate 1, 2, and 4 mg per Figure 11).",
    regions        = NULL,
    notes          = "Data set from Saunders DJ et al., Br J Ophthalmol 2015;99(11):1554-1559 (31 wet-AMD patients with serial aqueous humor free VEGF samples after 0.5 mg IVT ranibizumab). Individual fits used MATLAB lsqnonlin on log-transformed aqueous VEGF; K_D was estimated jointly by finding the value that reproduced the reference ocular half-life of ranibizumab (7.9 +/- 1.74 days, compiled in Hutton-Smith 2016). Retinal permeabilities and hydrodynamic radii are inherited from Hutton-Smith 2017 (rabbit intravitreal PK of therapeutic antibodies)."
  )

  ini({
    # ------------------------------------------------------------------------
    # Population log-normal parameters (Hutton-Smith 2018 Results, Figures 5
    # and 6). Individual fits gave t_1/2^(r) and V_in distributed log-normally
    # across the 31-patient sample; the paper reports the log-normal mu / sigma
    # with 95% CIs. Log-transformed here so the typical value (median) equals
    # exp(mu) and the IIV variance equals sigma^2. The paper's reported means
    # 7.9 days and 17.5 fmol/day correspond to arithmetic means of the
    # log-normal (mean = exp(mu + sigma^2/2)), which the parameterisation below
    # reproduces on simulation.
    # ------------------------------------------------------------------------
    lthalfr <- 2.03; label("Log ocular half-life of ranibizumab (log days)")  # Hutton-Smith 2018 Figure 6: mu = 2.03 (95% CI 1.94, 2.12); corresponds to a median 7.6 days and an arithmetic mean 7.9 days
    lVin    <- 2.76; label("Log VEGF production rate (log fmol/day)")         # Hutton-Smith 2018 Figure 5: mu = 2.76 (95% CI 2.65, 2.93); corresponds to a median 15.8 fmol/day and an arithmetic mean 17.5 fmol/day

    # Fixed structural parameters. K_D was estimated jointly by matching the
    # simulated t_1/2^(r) distribution to a literature reference value; the
    # paper reports a single point value with no uncertainty statistic, so it
    # is encoded as fixed here. k_off is fixed at the in-vitro estimate from
    # Yang et al. 2014 (Hutton-Smith 2018 reference 14).
    lKD   <- fixed(log(19000));  label("Log VEGF-ranibizumab dissociation constant (log pM)")             # Hutton-Smith 2018 Results (VEGF Production Rate) and Table 2: K_D = 19,000 pM (estimated jointly to reproduce t_1/2^(r) = 7.9 +/- 1.74 days)
    lkoff <- fixed(log(0.864));  label("Log single-site VEGF-ranibizumab off rate (log 1/day)")           # Hutton-Smith 2018 Table 2: k_off = 0.864 /day (Yang 2014 in vitro estimate; ref 14)

    # ------------------------------------------------------------------------
    # Inter-individual variability. Log-normal sigmas from Figures 5 and 6 map
    # directly to nlmixr2 eta variances (variance = sigma^2). No off-diagonal
    # covariance was reported.
    # ------------------------------------------------------------------------
    etalthalfr ~ 0.0625  # (0.25)^2; Hutton-Smith 2018 Figure 6: sigma = 0.25 (95% CI 0.20, 0.33)
    etalVin    ~ 0.1444  # (0.38)^2; Hutton-Smith 2018 Figure 5: sigma = 0.38 (95% CI 0.30, 0.51)

    # ------------------------------------------------------------------------
    # Residual error on the observed endpoint (aqueous humor free VEGF, pM).
    # The paper fit log-transformed aqueous VEGF (Methods, Fitting Protocol:
    # "penalizing the relative mean square error between the logarithms of the
    # numerical solution and the data") and reports that the residuals are
    # normally distributed with mean 0 and SD approximately 0.26 (Results,
    # Fitting Results and Figure 3b/c). Encoded here as a log-normal residual
    # so the SD applies on the log scale.
    # ------------------------------------------------------------------------
    expSd_Cc_vegf_aq <- 0.26; label("Log-normal residual SD on aqueous free VEGF (log pM)")  # Hutton-Smith 2018 Results (Fitting Results) and Figure 3b/c
  })

  model({
    # ------------------------------------------------------------------------
    # Anatomic constants (Hutton-Smith 2018 Table 2). All time units in days,
    # concentrations in pM (pmol/L), volumes in mL, areas in cm^2. 1 pM x 1 mL
    # equals 1 fmol, which is the amount unit used in the ODEs below.
    # ------------------------------------------------------------------------
    vol_ret <- 0.22   # mL   # Hutton-Smith 2018 Table 2 footnote b: retinal volume derived from Missel 2012 anatomy
    vol_vit <- 4.5    # mL   # Hutton-Smith 2018 Table 2: vitreal volume (refs 16, 17)
    vol_aq  <- 0.16   # mL   # Hutton-Smith 2018 Table 2: aqueous volume (refs 18, 19)
    S_ret   <- 9.71   # cm2  # Hutton-Smith 2018 Table 2 footnote b: retinal surface area derived from Missel 2012 anatomy
    CL      <- 3.6    # mL/day  # Hutton-Smith 2018 Table 2: aqueous humor outflow (ref 15)

    # Molecular weight of ranibizumab (Hutton-Smith 2018 Methods, paragraph
    # after Table 1). Used to convert IVT dose (mg) to amount (fmol) via
    # bioavailability scaling below.
    mw_ranib <- 48350  # g/mol  # Hutton-Smith 2018 ref 11: 48.35 kDa

    # ------------------------------------------------------------------------
    # Hydrodynamic radii (nm) and retinal permeabilities (cm/day) for each
    # chemical species (Hutton-Smith 2018 Table 3). Permeabilities in Table 3
    # are given in cm/s (x 10^-7); converted to cm/day via x 86400. Values
    # match eqs 15-16 evaluated at each species' R_h.
    # ------------------------------------------------------------------------
    pILM_v <- 1.89e-7 * 86400   # 0.01633 cm/day  # Hutton-Smith 2018 Table 3: p_ILM(V), R_h(V) = 2.39 nm
    pILM_r <- 1.89e-7 * 86400   # 0.01633 cm/day  # Hutton-Smith 2018 Table 3: p_ILM(R), R_h(R) = 2.45 nm
    pILM_c <- 1.79e-7 * 86400   # 0.01547 cm/day  # Hutton-Smith 2018 Table 3: p_ILM(VR), R_h(VR) = 3.29 nm
    pILM_h <- 1.73e-7 * 86400   # 0.01495 cm/day  # Hutton-Smith 2018 Table 3: p_ILM(RVR), R_h(RVR) = 4.07 nm

    pRPE_v <- 2.66e-7 * 86400   # 0.02298 cm/day  # Hutton-Smith 2018 Table 3: p_RPE(V)
    pRPE_r <- 2.63e-7 * 86400   # 0.02272 cm/day  # Hutton-Smith 2018 Table 3: p_RPE(R)
    pRPE_c <- 2.28e-7 * 86400   # 0.01970 cm/day  # Hutton-Smith 2018 Table 3: p_RPE(VR)
    pRPE_h <- 2.06e-7 * 86400   # 0.01780 cm/day  # Hutton-Smith 2018 Table 3: p_RPE(RVR)

    Rh_v <- 2.39   # nm   # Hutton-Smith 2018 Table 3: hydrodynamic radius of VEGF
    Rh_r <- 2.45   # nm   # Hutton-Smith 2018 Table 3: hydrodynamic radius of ranibizumab
    Rh_c <- 3.29   # nm   # Hutton-Smith 2018 Table 3: hydrodynamic radius of VR complex
    Rh_h <- 4.07   # nm   # Hutton-Smith 2018 Table 3: hydrodynamic radius of RVR complex

    # ------------------------------------------------------------------------
    # Individual estimated parameters.
    # ------------------------------------------------------------------------
    thalf_r <- exp(lthalfr + etalthalfr)   # days
    Vin     <- exp(lVin    + etalVin)      # fmol/day
    koff    <- exp(lkoff)                  # 1/day
    KD      <- exp(lKD)                    # pM
    kon     <- koff / KD                   # 1/(pM*day)  (Hutton-Smith 2018 Methods: k_on = k_off / K_D)

    # ------------------------------------------------------------------------
    # Species-specific half-lives and decay rates.
    # Eq 17: t_1/2^(i) proportional to R_h^(i); anchored to t_1/2^(r).
    # Eq 18: lambda^(i) = log(2) / t_1/2^(i).
    # ------------------------------------------------------------------------
    thalf_v  <- thalf_r * Rh_v / Rh_r
    thalf_c  <- thalf_r * Rh_c / Rh_r
    thalf_h  <- thalf_r * Rh_h / Rh_r

    lambda_v <- log(2) / thalf_v
    lambda_r <- log(2) / thalf_r
    lambda_c <- log(2) / thalf_c
    lambda_h <- log(2) / thalf_h

    # ------------------------------------------------------------------------
    # Vitreous-to-aqueous elimination rate k_el^(i), derived from the free-
    # species mass balance (Hutton-Smith 2018 Section S2 of the Supporting
    # Information, referenced but not reproduced in the main text). At the
    # long-time eigenmode where every ocular compartment for species i decays
    # exponentially at rate lambda^(i), the retinal-to-vitreous concentration
    # ratio is alpha^(i) = p_ILM^(i) / [(p_ILM^(i) + p_RPE^(i)) - lambda^(i) *
    # vol_ret / S_ret], and the vitreous mass balance gives k_el^(i) =
    # lambda^(i) + (S_ret / vol_vit) * p_ILM^(i) * (alpha^(i) - 1). For
    # ranibizumab (R_h = 2.45 nm, t_1/2^(r) = 7.5 days) this reproduces
    # k_el^(r) ~ 0.073 /day, which yields the vitreous:aqueous:retina
    # concentration ratios (10:1:5) reported in Results (Fitting Results).
    # ------------------------------------------------------------------------
    alpha_v <- pILM_v / ((pILM_v + pRPE_v) - lambda_v * vol_ret / S_ret)
    alpha_r <- pILM_r / ((pILM_r + pRPE_r) - lambda_r * vol_ret / S_ret)
    alpha_c <- pILM_c / ((pILM_c + pRPE_c) - lambda_c * vol_ret / S_ret)
    alpha_h <- pILM_h / ((pILM_h + pRPE_h) - lambda_h * vol_ret / S_ret)

    kel_v <- lambda_v + (S_ret / vol_vit) * pILM_v * (alpha_v - 1)
    kel_r <- lambda_r + (S_ret / vol_vit) * pILM_r * (alpha_r - 1)
    kel_c <- lambda_c + (S_ret / vol_vit) * pILM_c * (alpha_c - 1)
    kel_h <- lambda_h + (S_ret / vol_vit) * pILM_h * (alpha_h - 1)

    # ------------------------------------------------------------------------
    # Steady-state VEGF distribution before the intravitreal dose (all R, VR,
    # RVR pools at zero). Closed-form solution of the retina/vitreous/aqueous
    # steady-state balance for a single species with zero-order production
    # V_in in the retina, ILM diffusion between retina and vitreous, RPE
    # elimination from the retina, k_el^(v) transport vitreous -> aqueous, and
    # C_L outflow from the aqueous.
    # ------------------------------------------------------------------------
    A_v <- S_ret * pILM_v
    B_v <- S_ret * pRPE_v
    vegf_vit_ss <- Vin * A_v / (A_v * B_v + (A_v + B_v) * kel_v * vol_vit)
    vegf_ret_ss <- (A_v * vegf_vit_ss + Vin) / (A_v + B_v)
    vegf_aq_ss  <- kel_v * vol_vit * vegf_vit_ss / CL

    # ------------------------------------------------------------------------
    # Concentrations (pM) from the state amounts (fmol) via compartment
    # volumes (mL). 1 fmol / 1 mL = 1 pmol/L = 1 pM.
    # ------------------------------------------------------------------------
    vegf_ret_c  <- vegf_ret  / vol_ret
    ranib_ret_c <- ranib_ret / vol_ret
    vr_ret_c    <- vr_ret    / vol_ret
    rvr_ret_c   <- rvr_ret   / vol_ret

    vegf_vit_c  <- vegf_vit  / vol_vit
    ranib_vit_c <- ranib_vit / vol_vit
    vr_vit_c    <- vr_vit    / vol_vit
    rvr_vit_c   <- rvr_vit   / vol_vit

    vegf_aq_c   <- vegf_aq   / vol_aq
    ranib_aq_c  <- ranib_aq  / vol_aq
    vr_aq_c     <- vr_aq     / vol_aq
    rvr_aq_c    <- rvr_aq    / vol_aq

    # ------------------------------------------------------------------------
    # Net forward binding-reaction rates per compartment (pM/day). Positive
    # values represent net progression toward the bound species.
    #   r1: V + R  <-> VR   (2*k_on forward due to two V binding sites)
    #   r2: R + VR <-> RVR  (2*k_off reverse due to two dissociation paths)
    # Hutton-Smith 2018 Methods (Model Description) eqs 1-2 and Figure 1.
    # ------------------------------------------------------------------------
    r1_ret <- 2 * kon * vegf_ret_c * ranib_ret_c - koff * vr_ret_c
    r2_ret <- kon * ranib_ret_c * vr_ret_c - 2 * koff * rvr_ret_c

    r1_vit <- 2 * kon * vegf_vit_c * ranib_vit_c - koff * vr_vit_c
    r2_vit <- kon * ranib_vit_c * vr_vit_c - 2 * koff * rvr_vit_c

    r1_aq  <- 2 * kon * vegf_aq_c  * ranib_aq_c  - koff * vr_aq_c
    r2_aq  <- kon * ranib_aq_c * vr_aq_c - 2 * koff * rvr_aq_c

    # ------------------------------------------------------------------------
    # ODE system (Hutton-Smith 2018 eqs 3-14). State amounts in fmol; LHS in
    # fmol/day. Each concentration-form paper equation is multiplied through
    # by the appropriate compartment volume to give the amount-form LHS.
    # ------------------------------------------------------------------------

    # Retina (eqs 3-6)
    d/dt(vegf_ret)  <- -vol_ret * r1_ret            - S_ret * (pILM_v + pRPE_v) * vegf_ret_c  + S_ret * pILM_v * vegf_vit_c  + Vin
    d/dt(ranib_ret) <- -vol_ret * (r1_ret + r2_ret) - S_ret * (pILM_r + pRPE_r) * ranib_ret_c + S_ret * pILM_r * ranib_vit_c
    d/dt(vr_ret)    <-  vol_ret * (r1_ret - r2_ret) - S_ret * (pILM_c + pRPE_c) * vr_ret_c    + S_ret * pILM_c * vr_vit_c
    d/dt(rvr_ret)   <-  vol_ret * r2_ret            - S_ret * (pILM_h + pRPE_h) * rvr_ret_c   + S_ret * pILM_h * rvr_vit_c

    # Vitreous (eqs 7-10)
    d/dt(vegf_vit)  <- -vol_vit * r1_vit            + S_ret * pILM_v * (vegf_ret_c  - vegf_vit_c)  - kel_v * vegf_vit
    d/dt(ranib_vit) <- -vol_vit * (r1_vit + r2_vit) + S_ret * pILM_r * (ranib_ret_c - ranib_vit_c) - kel_r * ranib_vit
    d/dt(vr_vit)    <-  vol_vit * (r1_vit - r2_vit) + S_ret * pILM_c * (vr_ret_c    - vr_vit_c)    - kel_c * vr_vit
    d/dt(rvr_vit)   <-  vol_vit * r2_vit            + S_ret * pILM_h * (rvr_ret_c   - rvr_vit_c)   - kel_h * rvr_vit

    # Aqueous (eqs 11-14)
    d/dt(vegf_aq)   <- -vol_aq * r1_aq            + kel_v * vegf_vit  - CL * vegf_aq_c
    d/dt(ranib_aq)  <- -vol_aq * (r1_aq + r2_aq) + kel_r * ranib_vit - CL * ranib_aq_c
    d/dt(vr_aq)     <-  vol_aq * (r1_aq - r2_aq) + kel_c * vr_vit    - CL * vr_aq_c
    d/dt(rvr_aq)    <-  vol_aq * r2_aq           + kel_h * rvr_vit   - CL * rvr_aq_c

    # ------------------------------------------------------------------------
    # Initial conditions. Pre-dose steady state for VEGF; ranibizumab and its
    # complexes are zero until the intravitreal dose is administered.
    # ------------------------------------------------------------------------
    vegf_ret(0) <- vegf_ret_ss * vol_ret
    vegf_vit(0) <- vegf_vit_ss * vol_vit
    vegf_aq(0)  <- vegf_aq_ss  * vol_aq

    # ------------------------------------------------------------------------
    # Bioavailability-style unit conversion for the intravitreal dose. The
    # event-table amt is expressed in mg (the clinical dose unit); the state
    # amount is fmol, so the delivered amount is amt (mg) * 1e-3 g/mg * 1e15
    # fmol/mol / mw_ranib (g/mol) = amt * 1e12 / 48350 fmol. A dose of amt =
    # 0.5 gives 1.034e7 fmol delivered, matching r_vit(0) = 2.30e6 pM at
    # vol_vit = 4.5 mL (Hutton-Smith 2018 Methods, paragraph after Table 1).
    # ------------------------------------------------------------------------
    f(ranib_vit) <- 1e12 / mw_ranib

    # ------------------------------------------------------------------------
    # Observation variables (concentrations, pM). Cc_vegf_aq is the fitted
    # endpoint (Saunders 2015 aqueous humor VEGF); the other outputs are for
    # simulation and were not fit to observations.
    # ------------------------------------------------------------------------
    Cc_vegf_aq   <- vegf_aq_c
    Cc_vegf_vit  <- vegf_vit_c
    Cc_vegf_ret  <- vegf_ret_c
    Cc_ranib_ret <- ranib_ret_c
    Cc_ranib_vit <- ranib_vit_c
    Cc_ranib_aq  <- ranib_aq_c

    # Residual error on the fitted endpoint.
    Cc_vegf_aq ~ lnorm(expSd_Cc_vegf_aq)
  })
}
