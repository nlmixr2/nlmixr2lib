Yuan_2019_concizumab <- function() {
  description <- "QSP. Systems PK/PD model for concizumab (humanized anti-TFPI IgG4) describing binding to both membrane-bound TFPI (mTFPI; non-linear clearance via receptor-mediated endocytosis) and soluble TFPI (sTFPI; linear clearance via FcRn-recycled pinocytosis) in a minimal physiologically-based PK framework with two nested endothelial endosome compartments. Parameter values for 70 kg adult humans (Yuan 2019 Tables 1-2); the paper also tabulates monkey and rabbit parameter sets."
  reference   <- "Yuan D, Rode F, Cao Y. A systems pharmacokinetic/pharmacodynamic model for concizumab to explore the potential of anti-TFPI recycling antibodies. Eur J Pharm Sci. 2019 Oct 1;138:105032. doi:10.1016/j.ejps.2019.105032. PMID 31374317. mTFPI baseline, kdegm = kint, and koff optimized in the reduced PK/PD model using human PK/PD data digitized from Chowdary 2015 (J Thromb Haemost 13:743-754). Linear-clearance endosome parameters (CLup, CLe, krec, k1on, k1off, FcRn_b) inherited from Yuan 2018 (J Pharmacokinet Pharmacodyn 45:851-864), calibrated using adalimumab."
  vignette    <- "Yuan_2019_concizumab"
  units       <- list(time = "day", dosing = "nmol", concentration = "nM")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24,
    n_studies      = 1,
    age_range      = "adult (per Chowdary 2015 source study)",
    weight_range   = "70 kg reference body weight (Yuan 2019 Table 2)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "Healthy adult males (Chowdary 2015 first-in-human single ascending dose study of concizumab)",
    dose_range     = "IV 0.25, 1, 3, 9 mg/kg single (used for parameter optimisation); SC 1, 3 mg/kg single (external validation cohort).",
    regions        = NA,
    notes          = paste0(
      "Model is a typical-value mechanistic (systems QSP) simulator; ",
      "no IIV or residual error is reported for the systems model itself. ",
      "Residual-error terms in Yuan 2019 Table 1 were ADAPT 5 maximum-likelihood ",
      "fit residuals from the reduced PK/PD model and apply to fitting only, not ",
      "to simulation. Parameters can be replaced with monkey or rabbit values ",
      "(Yuan 2019 Table 2) to simulate those species. ",
      "Concizumab clinical PK/PD source: Chowdary P et al. J Thromb Haemost. ",
      "2015 May;13(5):743-754 (PMID 25641556)."
    ),
    model_class    = "QSP / systems PK-PD (mPBPK base with two nested endothelial endosome compartments, explicit TMDD for mTFPI and sTFPI, and explicit FcRn recycling)",
    n_states       = 23,
    upstream_model = "Yuan D, Krzyzanski W, Cao Y. A general pharmacokinetic model for the disposition of monoclonal antibodies against soluble antigens. J Pharmacokinet Pharmacodyn. 2018 Dec;45(6):851-864 (linear-clearance endosome calibration with adalimumab)."
  )

  ini({
    # Physiological parameters (Yuan 2019 Table 2, human column)
    lvp    <- fixed(log(2.6)); label("Plasma volume (L)")                                                                       # Yuan 2019 Table 2 (Cao 2013)
    visf   <- fixed(15.6);    label("Total interstitial fluid volume (L)")                                                      # Yuan 2019 Table 2 (Cao 2013)
    vlymph <- fixed(5.2);     label("Lymph volume (L)")                                                                         # Yuan 2019 Table 2 (Cao 2013)
    ltot   <- fixed(2.904);   label("Total lymph flow (L/day)")                                                                 # Yuan 2019 Table 2 (Cao 2013)
    sigma_tight <- fixed(0.945);   label("Vascular reflection coefficient for tight tissues (unitless)")                             # Yuan 2019 Table 2 (Cao 2013)
    sigma_leaky <- fixed(0.697);   label("Vascular reflection coefficient for leaky tissues (unitless)")                             # Yuan 2019 Table 2 (Cao 2013)
    sigmal <- fixed(0.2);     label("Lymphatic vascular reflection coefficient (unitless)")                                     # Yuan 2019 Table 2 (Cao 2013)
    kp     <- fixed(0.4);     label("Fraction of interstitial space available for IgG4 antibody distribution (unitless)")       # Yuan 2019 Table 2 (Cao 2013); 0.4 for IgG4, 0.8 for IgG1

    # Target (mTFPI / sTFPI) parameters (Yuan 2019 Tables 1-2, human)
    kon_per_nm_day  <- fixed(387.072);  label("Antibody-TFPI association rate constant in plasma (1/(nM*day))")                 # Yuan 2019 Table 2: 4.48e6 1/(M*s) (Hilden 2012) * 86400 s/day / 1e9 nM/M = 387.072 1/(nM*day)
    koff_per_day    <- fixed(150.336);  label("Antibody-TFPI complex dissociation rate constant in plasma (1/day)")             # Yuan 2019 Table 1 (optimised in reduced model): 1.74e-3 1/s * 86400 s/day = 150.336 1/day; CV 56.12%
    keon_per_nm_day <- fixed(387.072);  label("Antibody-TFPI association rate constant in endosomes (1/(nM*day))")              # Yuan 2019 Section 2.4: assumed equal to plasma kon for concizumab (not engineered for pH-dependent release)
    keoff_per_day   <- fixed(150.336);  label("Antibody-TFPI complex dissociation rate constant in endosomes (1/day)")          # Yuan 2019 Section 2.4: assumed equal to plasma koff for concizumab (not engineered for pH-dependent release)
    stfpi_b         <- fixed(1.6);      label("Plasma soluble TFPI baseline concentration (nM)")                                # Yuan 2019 Table 2 (Hansen 2014)
    kdegs           <- fixed(1.18);     label("Total degradation rate constant of soluble TFPI (1/day)")                        # Yuan 2019 Table 2 (Farrokhi 2018, clinical pulse-chase)
    mtfpi_b         <- fixed(20.03);    label("Plasma membrane-bound TFPI baseline concentration (nM)")                         # Yuan 2019 Table 1 (optimised in reduced model); CV 13.45%
    kdegm           <- fixed(1.155);    label("Membrane-bound TFPI internalisation rate constant kdegm = kint (1/day)")         # Yuan 2019 Table 1 (optimised in reduced model with constraint kdegm = kint); CV 26.76%

    # Endosome / FcRn parameters (Yuan 2019 Table 2, calibrated in Yuan 2018 with adalimumab)
    clup             <- fixed(1.48);    label("Endothelial nonspecific pinocytosis clearance into endosome 2 (L/day)")          # Yuan 2019 Table 2 (Yuan 2018)
    cle              <- fixed(5.75);    label("Lysosomal catabolism clearance from endosomes (L/day)")                          # Yuan 2019 Table 2 (Yuan 2018)
    krec             <- fixed(124.7);   label("FcRn-bound antibody recycling rate from endosomes to plasma (1/day)")            # Yuan 2019 Table 2 (Hopkins and Trowbridge 1983)
    ve1              <- fixed(0.0035);  label("Endosome volume for mTFPI-mediated endocytosis (Ve1, L)")                        # Yuan 2019 Table 2 (Li 2014): 0.005% of 70 kg body weight = 0.0035 L
    ve2              <- fixed(0.0035);  label("Endosome volume for nonspecific pinocytosis (Ve2, L)")                           # Yuan 2019 Table 2 (Li 2014): 0.005% of 70 kg body weight = 0.0035 L
    fcrn_b           <- fixed(49800);   label("FcRn baseline concentration in endosomes (nM)")                                  # Yuan 2019 Table 2 (Shah and Betts 2012)
    k1on_per_nm_day  <- fixed(20.8224); label("Antibody-FcRn association rate constant in endosomes (1/(nM*day))")              # Yuan 2019 Table 2: 2.41e5 1/(M*s) (adalimumab-human FcRn, Suzuki 2010) * 86400 s/day / 1e9 nM/M = 20.8224 1/(nM*day)
    k1off_per_day    <- fixed(13996.8); label("Antibody-FcRn complex dissociation rate constant in endosomes (1/day)")          # Yuan 2019 Table 2: 0.162 1/s (adalimumab-human FcRn, Suzuki 2010) * 86400 s/day = 13996.8 1/day

    # SC absorption (Yuan 2019 Table 2 -- fixed to monkey values; human SC absorption not optimised)
    fsc <- fixed(0.93);       label("SC bioavailability of concizumab (unitless)")                                              # Yuan 2019 Table 2 (Agerso 2014; monkey value assumed for human)
    lka <- fixed(log(0.231)); label("SC first-order absorption rate constant (1/day)")                                          # Yuan 2019 Table 2 (Agerso 2014; monkey value assumed for human)
  })

  model({
    # Back-transform log-scale parameters
    vp <- exp(lvp)
    ka <- exp(lka)

    # Tissue partitioning constants (Cao 2013 mPBPK)
    v1 <- 0.65 * visf
    v2 <- 0.35 * visf
    l1 <- 0.33 * ltot
    l2 <- 0.67 * ltot

    # Available distribution volumes
    v1eff <- v1 * kp
    v2eff <- v2 * kp

    # Endogenous baseline AMOUNTS (Yuan 2019 Appendix II ICs)
    bl_stfpi_p    <- stfpi_b * vp
    bl_mtfpi_p    <- mtfpi_b * vp
    bl_mtfpi_e1   <- mtfpi_b * kdegm * vp * ve1 / cle
    bl_stfpi_e2   <- stfpi_b * clup * ve2 / cle
    bl_fcrn_e1    <- fcrn_b * ve1
    bl_fcrn_e2    <- fcrn_b * ve2

    # Plasma synthesis rate of sTFPI (nmol/day) and pinocytosis-independent
    # degradation rate constant (Appendix II Eq. (3))
    ksyns_amt    <- kdegs * stfpi_b * vp
    kdegs_prime  <- kdegs - clup / vp

    # Plasma synthesis rate of mTFPI (nmol/day) (Appendix II Eq. (5))
    ksynm_amt    <- kdegm * mtfpi_b * vp

    # Concentrations from amounts (nM)
    c_ap          <- a_p / vp
    c_astfpi_p    <- astfpi_p / vp
    c_amtfpi_p    <- amtfpi_p / vp
    c_stfpi_p    <- stfpi_p / vp
    c_mtfpi_p    <- mtfpi_p / vp
    c_a_e1        <- a_e1 / ve1
    c_mtfpi_e1   <- mtfpi_e1 / ve1
    c_amtfpi_e1  <- amtfpi_e1 / ve1
    c_fcrn_e1    <- fcrn_e1 / ve1
    c_fcrna_e1   <- fcrna_e1 / ve1
    c_a_e2        <- a_e2 / ve2
    c_stfpi_e2   <- stfpi_e2 / ve2
    c_astfpi_e2  <- astfpi_e2 / ve2
    c_fcrn_e2    <- fcrn_e2 / ve2
    c_fcrna_e2   <- fcrna_e2 / ve2
    c_fcrnastfpi_e2 <- fcrnastfpi_e2 / ve2
    c_a_t         <- a_t / v1eff
    c_astfpi_t   <- astfpi_t / v1eff
    c_a_lk        <- a_lk / v2eff
    c_astfpi_lk  <- astfpi_lk / v2eff
    c_a_lm        <- a_lm / vlymph
    c_astfpi_lm  <- astfpi_lm / vlymph

    # Plasma binding fluxes (nmol/day): rate * concentration * concentration * volume
    rsp_bind_plasma <- kon_per_nm_day * c_ap * c_stfpi_p * vp - koff_per_day * astfpi_p
    rmp_bind_plasma <- kon_per_nm_day * c_ap * c_mtfpi_p * vp - koff_per_day * amtfpi_p

    # Endosome 1 (mTFPI-mediated endocytosis pathway) fluxes
    re1_target_bind <- keon_per_nm_day * c_a_e1 * c_mtfpi_e1 * ve1 - keoff_per_day * amtfpi_e1
    re1_fcrn_bind   <- k1on_per_nm_day * c_a_e1 * c_fcrn_e1 * ve1 - k1off_per_day * fcrna_e1

    # Endosome 2 (nonspecific pinocytosis pathway) fluxes
    re2_target_bind        <- keon_per_nm_day * c_a_e2 * c_stfpi_e2 * ve2 - keoff_per_day * astfpi_e2
    re2_fcrn_free          <- k1on_per_nm_day * c_a_e2 * c_fcrn_e2 * ve2 - k1off_per_day * fcrna_e2
    re2_fcrn_astfpi        <- k1on_per_nm_day * c_astfpi_e2 * c_fcrn_e2 * ve2 - k1off_per_day * fcrnastfpi_e2
    re2_target_on_fcrnab   <- keon_per_nm_day * c_fcrna_e2 * c_stfpi_e2 * ve2 - keoff_per_day * fcrnastfpi_e2

    # Eq. (1) SC absorption depot
    d/dt(depot) <- -ka * depot

    # Eq. (2) Free antibody in plasma (a_p, amount in nmol)
    d/dt(a_p) <- ka * depot -
                 rsp_bind_plasma -
                 rmp_bind_plasma -
                 (1 - sigma_tight) * l1 * c_ap -
                 (1 - sigma_leaky) * l2 * c_ap +
                 ltot * c_a_lm -
                 clup * c_ap +
                 krec * fcrna_e1 +
                 krec * fcrna_e2

    # Eq. (3) Free soluble TFPI in plasma
    d/dt(stfpi_p) <- -rsp_bind_plasma +
                     ksyns_amt -
                     kdegs_prime * stfpi_p -
                     clup * c_stfpi_p

    # Eq. (4) Antibody-sTFPI complex in plasma
    d/dt(astfpi_p) <- rsp_bind_plasma -
                      (1 - sigma_tight) * l1 * c_astfpi_p -
                      (1 - sigma_leaky) * l2 * c_astfpi_p +
                      ltot * c_astfpi_lm -
                      clup * c_astfpi_p +
                      krec * fcrnastfpi_e2

    # Eq. (5) Free membrane-bound TFPI in plasma
    d/dt(mtfpi_p) <- -rmp_bind_plasma +
                     ksynm_amt -
                     kdegm * mtfpi_p

    # Eq. (6) Antibody-mTFPI complex in plasma
    d/dt(amtfpi_p) <- rmp_bind_plasma -
                      kdegm * amtfpi_p

    # Eq. (7) Free antibody in endosome 1
    d/dt(a_e1) <- -re1_target_bind -
                  re1_fcrn_bind -
                  cle * c_a_e1

    # Eq. (8) Free mTFPI in endosome 1
    d/dt(mtfpi_e1) <- -re1_target_bind +
                      kdegm * mtfpi_p -
                      cle * c_mtfpi_e1

    # Eq. (9) Antibody-mTFPI complex in endosome 1
    d/dt(amtfpi_e1) <- re1_target_bind +
                       kdegm * amtfpi_p -
                       cle * c_amtfpi_e1

    # Eq. (10) Free FcRn in endosome 1
    d/dt(fcrn_e1) <- -re1_fcrn_bind +
                     krec * fcrna_e1

    # Eq. (11) FcRn-bound antibody in endosome 1
    d/dt(fcrna_e1) <- re1_fcrn_bind -
                      krec * fcrna_e1

    # Eq. (12) Free antibody in endosome 2
    d/dt(a_e2) <- -re2_target_bind -
                  re2_fcrn_free +
                  clup * c_ap -
                  cle * c_a_e2

    # Eq. (13) Free sTFPI in endosome 2
    d/dt(stfpi_e2) <- -re2_target_bind -
                      re2_target_on_fcrnab +
                      clup * c_stfpi_p -
                      cle * c_stfpi_e2

    # Eq. (14) Antibody-sTFPI complex in endosome 2
    d/dt(astfpi_e2) <- re2_target_bind -
                       re2_fcrn_astfpi +
                       clup * c_astfpi_p -
                       cle * c_astfpi_e2

    # Eq. (15) Free FcRn in endosome 2
    d/dt(fcrn_e2) <- -re2_fcrn_free -
                     re2_fcrn_astfpi +
                     krec * (fcrna_e2 + fcrnastfpi_e2)

    # Eq. (16) FcRn-bound antibody in endosome 2
    d/dt(fcrna_e2) <- re2_fcrn_free -
                      re2_target_on_fcrnab -
                      krec * fcrna_e2

    # Eq. (17) FcRn-bound antibody-sTFPI complex in endosome 2
    d/dt(fcrnastfpi_e2) <- re2_fcrn_astfpi +
                           re2_target_on_fcrnab -
                           krec * fcrnastfpi_e2

    # Eq. (18) Free antibody in tight tissue
    d/dt(a_t) <- (1 - sigma_tight) * l1 * c_ap -
                 (1 - sigmal) * l1 * c_a_t

    # Eq. (19) Antibody-sTFPI complex in tight tissue
    d/dt(astfpi_t) <- (1 - sigma_tight) * l1 * c_astfpi_p -
                      (1 - sigmal) * l1 * c_astfpi_t

    # Eq. (20) Free antibody in leaky tissue
    d/dt(a_lk) <- (1 - sigma_leaky) * l2 * c_ap -
                  (1 - sigmal) * l2 * c_a_lk

    # Eq. (21) Antibody-sTFPI complex in leaky tissue
    d/dt(astfpi_lk) <- (1 - sigma_leaky) * l2 * c_astfpi_p -
                       (1 - sigmal) * l2 * c_astfpi_lk

    # Eq. (22) Free antibody in lymph
    d/dt(a_lm) <- (1 - sigmal) * l1 * c_a_t +
                  (1 - sigmal) * l2 * c_a_lk -
                  ltot * c_a_lm

    # Eq. (23) Antibody-sTFPI complex in lymph
    d/dt(astfpi_lm) <- (1 - sigmal) * l1 * c_astfpi_t +
                       (1 - sigmal) * l2 * c_astfpi_lk -
                       ltot * c_astfpi_lm

    # SC bioavailability (applied to dose amounts targeting depot)
    f(depot) <- fsc

    # Endogenous initial conditions (Yuan 2019 Appendix II)
    stfpi_p(0)  <- bl_stfpi_p
    mtfpi_p(0)  <- bl_mtfpi_p
    mtfpi_e1(0) <- bl_mtfpi_e1
    stfpi_e2(0) <- bl_stfpi_e2
    fcrn_e1(0)  <- bl_fcrn_e1
    fcrn_e2(0)  <- bl_fcrn_e2

    # Outputs (concentrations in nM)
    # Total antibody = free + sTFPI-complex + mTFPI-complex in plasma
    Cc      <- c_ap + c_astfpi_p + c_amtfpi_p
    # Free soluble TFPI in plasma (PD marker)
    sTFPI   <- c_stfpi_p
  })
}
