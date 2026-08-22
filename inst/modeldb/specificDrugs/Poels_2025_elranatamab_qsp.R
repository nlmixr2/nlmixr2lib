Poels_2025_elranatamab_qsp <- function() {
  description <- "QSP. Quantitative systems pharmacology model of the BCMA x CD3 bispecific antibody elranatamab in relapsed/refractory multiple myeloma, calibrated to MagnetisMM-1 (NCT03269136) and MagnetisMM-3 (NCT04649359). Three compartments (central, peripheral, bone marrow / site of action) carry mass-action BsAb binding to membrane BCMA, soluble BCMA (drug sink) and T-cell CD3, forming BsAb-CD3-BCMA trimers that drive Simeoni-type myeloma-cell killing with time-dependent resistance, IL-6-like cytokine release with dose-to-dose attenuation, cytokine- and trimer-dependent T-cell retention in marrow, and M-protein / free light chain paraprotein read-outs. Builds on the Betts 2019 CD3-bispecific framework; see modellib('Betts_2019_pf_06671008_qsp')."
  reference <- paste(
    "Poels KE, Elmeliegy M, Hibma J, Wang D, Musante CJ, Shtylla B.",
    "Leveraging quantitative systems pharmacology modeling for elranatamab",
    "regimen optimization in relapsed or refractory multiple myeloma.",
    "npj Syst Biol Appl. 2025;11:102. doi:10.1038/s41540-025-00585-z.",
    "Structural framework adapted from Betts et al. (2019) AAPS J 21:66",
    "(see modellib('Betts_2019_pf_06671008_qsp')) and the cytokine-release",
    "framework of Chen et al. (2019) Clin Transl Sci 12:600-608.",
    sep = " "
  )
  vignette <- "Poels_2025_elranatamab"

  # Time is hours throughout: every rate constant in Supplementary Table 3 is
  # reported per hour. `depot` holds a subcutaneous AMOUNT (pmol) because
  # Supplementary Eq 1 divides the absorption flux by Vc; every other drug
  # state holds a CONCENTRATION (pM), exactly as printed in the supplement.
  # Doses are therefore given in pmol to `depot` for SC administration, or as
  # the central-compartment concentration equivalent (dose_pmol / vc) to
  # `central` for IV administration (the same convention as
  # Betts_2019_pf_06671008_qsp). Elranatamab MW = 148.5 kDa
  # (ELREXFIO US prescribing information, Description section, on disk as
  # PMC12402305_aux_ELREXFIO_USPI_dailymed.txt), so 76 mg = 5.1178e5 pmol.
  units <- list(time = "h", dosing = "pmol", concentration = "pM")

  # States whose role has no canonical entry in inst/references/compartment-names.md.
  # Canonical states used unchanged: depot, central, peripheral1, tumor,
  # drug_cd3_tumor, trimer, cycling_cells, damaged_cells1-3. `drug_bcma_tumor`
  # is registered as a new canonical in this PR, following the per-antigen
  # extension path documented under `drug_pcad_tumor`.
  paper_specific_compartments <- c(
    "sbcma_central", "sbcma_tumor",
    "drug_sbcma_central", "drug_sbcma_tumor", "drug_cd3_central",
    "tcell_central", "tcell_tumor",
    "mprotein", "flc",
    "cytokine_central", "cytokine_tumor",
    "cytokine_transit1", "cytokine_transit2", "cytokine_transit3",
    "cytokine_transit4", "cytokine_transit5",
    "cytokine_auc"
  )

  compartmentData <- list(
    depot              = list(analyte = "elranatamab",                     units = "pmol", specimen = "administration site", verified = TRUE),
    central            = list(analyte = "elranatamab",                     units = "pM",   specimen = "serum",               verified = TRUE),
    peripheral1        = list(analyte = "elranatamab",                     units = "pM",   specimen = "serum",               verified = TRUE),
    sbcma_central      = list(analyte = "soluble BCMA",                    units = "pM",   specimen = "serum",               verified = TRUE),
    drug_sbcma_central = list(analyte = "elranatamab-soluble BCMA dimer",  units = "pM",   specimen = "serum",               verified = TRUE),
    drug_cd3_central   = list(analyte = "elranatamab-CD3 dimer",           units = "pM",   specimen = "serum",               verified = TRUE),
    tcell_central      = list(analyte = "CD3+ T cells",                    units = "cells/uL", specimen = "whole blood",           verified = TRUE),
    cytokine_central   = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "serum",             verified = TRUE),
    mprotein           = list(analyte = "monoclonal M-protein",            units = "g/L",  specimen = "serum",               verified = TRUE),
    flc                = list(analyte = "involved serum free light chain", units = "mg/L", specimen = "serum",               verified = TRUE),
    tumor              = list(analyte = "elranatamab",                     units = "pM",   specimen = "tumor",         verified = TRUE),
    sbcma_tumor        = list(analyte = "soluble BCMA",                    units = "pM",   specimen = "tumor",         verified = TRUE),
    drug_bcma_tumor    = list(analyte = "elranatamab-membrane BCMA dimer", units = "pM",   specimen = "tumor",         verified = TRUE),
    drug_sbcma_tumor   = list(analyte = "elranatamab-soluble BCMA dimer",  units = "pM",   specimen = "tumor",         verified = TRUE),
    drug_cd3_tumor     = list(analyte = "elranatamab-CD3 dimer",           units = "pM",   specimen = "tumor",         verified = TRUE),
    trimer             = list(analyte = "elranatamab-CD3-BCMA trimer",     units = "pM",   specimen = "tumor",         verified = TRUE),
    cycling_cells      = list(analyte = "multiple myeloma cells",          units = "cells/uL", specimen = "tumor",     verified = TRUE),
    damaged_cells1     = list(analyte = "multiple myeloma cells",          units = "cells/uL", specimen = "tumor",     verified = TRUE),
    damaged_cells2     = list(analyte = "multiple myeloma cells",          units = "cells/uL", specimen = "tumor",     verified = TRUE),
    damaged_cells3     = list(analyte = "multiple myeloma cells",          units = "cells/uL", specimen = "tumor",     verified = TRUE),
    tcell_tumor        = list(analyte = "CD3+ T cells",                    units = "cells/uL", specimen = "tumor",     verified = TRUE),
    cytokine_tumor     = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "tumor",       verified = TRUE),
    cytokine_transit1  = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "not applicable",    verified = TRUE),
    cytokine_transit2  = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "not applicable",    verified = TRUE),
    cytokine_transit3  = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "not applicable",    verified = TRUE),
    cytokine_transit4  = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "not applicable",    verified = TRUE),
    cytokine_transit5  = list(analyte = "pro-inflammatory cytokine (IL-6)", units = "pg/mL", specimen = "not applicable",    verified = TRUE),
    cytokine_auc       = list(analyte = "cumulative cytokine exposure",    units = "pg/mL*h", specimen = "not applicable",   verified = TRUE)
  )

  # The published model carries no covariate columns: patient-to-patient
  # heterogeneity is expressed entirely through sampled baseline states
  # (ini() parameters ending in `0`) and the nine perturbed structural
  # parameters of Supplementary Table 1. Baseline sBCMA -- the paper's
  # stratifying biomarker -- is a model STATE, not a data column, so it is
  # documented here rather than in covariateData.
  covariatesDataExcluded <- list(
    SBCMA = list(
      description = "Baseline soluble BCMA concentration used to stratify response (high >= 100 ng/mL vs low < 100 ng/mL).",
      units       = "ng/mL",
      type        = "continuous",
      notes       = paste(
        "Not a data covariate in this model: baseline sBCMA enters as the",
        "initial condition of the `sbcma_central` state via the",
        "`sbcma0_central` parameter (pM). Convert with the paper's own scale,",
        "185 pM per ng/mL, so the published 100 ng/mL cut-off is 18500 pM",
        "(Supplementary Table 2 range 185-153520 pM maps onto the 1-830 ng/mL",
        "plausible-patient range of Figure 2b)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = "120 virtual patients per virtual population; 10 virtual populations calibrated (1200 virtual patients total), drawn from >4000 plausible patients out of 10000 sampled parametrizations",
    n_studies      = 2,
    disease_state  = "Relapsed or refractory multiple myeloma (RRMM); triple-class refractory and BCMA-treatment naive",
    dose_range     = "Subcutaneous 80, 130, 215, 360, 600 and 1000 ug/kg QW (MagnetisMM-1 Part 1); 600 then 1000 ug/kg QW or Q2W (Part 1.1); 44 mg then 76 mg QW (Part 2A); 12/32/76 mg on C1D1/C1D4/C1D8 then 76 mg QW, stepping down to 76 mg Q2W after cycle 6 and 76 mg Q4W after cycle 12 (MagnetisMM-3 Cohort A). Simulated dose-escalation arms 16, 28, 44, 76 and 152 mg QW plus 76 mg Q2W.",
    regions        = "Multicentre international (MagnetisMM-1 C1071001; MagnetisMM-3 C1071003)",
    notes          = paste(
      "Calibration data are summarised in Supplementary Table 4: MagnetisMM-1",
      "(C1071001) Part 1 IV n = 23, Part 1 SC n = 30, Part 1 C-D n = 13,",
      "Part 1.1 QW n = 7, Part 1.1 Q2W n = 13, Part 2A n = 15; MagnetisMM-3",
      "(C1071003) Part A n = 120. The paper reports no age, weight, sex or",
      "race distribution for the calibration cohorts, so those fields are",
      "omitted rather than guessed. Virtual-patient baselines were sampled",
      "from distributions fitted to pooled MagnetisMM-1 and MagnetisMM-3",
      "Cohort A baseline data (Supplementary Table 2, Figure 2b); ~70% of the",
      "virtual population had baseline sBCMA < 100 ng/mL and ~30% >= 100 ng/mL."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Transform policy: parameters the paper reports as estimated or
    # calibrated are log-transformed so a refit stays positive; parameters
    # the paper fixed from nonclinical data, literature or an upstream
    # popPK are carried bare (or as fixed(log(.))) inside fixed().
    # ------------------------------------------------------------------

    # --- BsAb binding kinetics (Supplementary Table 3, Sect. 1.1) --------
    # Paper notation: kon1/koff1 = BCMA arm, kon2/koff2 = CD3 arm.
    kon_bcma  <- fixed(0.0054);  label("Association rate of elranatamab to BCMA (1/pM/h)")          # Suppl Table 3 Sect 1.1: kon1, fixed internal nonclinical data
    koff_bcma <- fixed(0.162);   label("Dissociation rate of elranatamab from BCMA (1/h)")          # Suppl Table 3 Sect 1.1: koff1, fixed internal nonclinical data
    kon_cd3   <- fixed(0.00198); label("Association rate of elranatamab to CD3 (1/pM/h)")           # Suppl Table 3 Sect 1.1: kon2, fixed internal nonclinical data
    koff_cd3  <- fixed(82.17);   label("Dissociation rate of elranatamab from CD3 (1/h)")           # Suppl Table 3 Sect 1.1: koff2, fixed internal nonclinical data

    # --- Soluble BCMA turnover and transport ----------------------------
    kdeg_sbcma      <- fixed(0.01925); label("Degradation rate of soluble BCMA (1/h)")              # Suppl Table 3 Sect 1.1: kdeg_sBCMA, fixed internal nonclinical data
    ktr_sbcma_bm_c  <- fixed(0.1);     label("Transport rate of soluble BCMA and its drug dimer, bone marrow to central (1/h)")  # Suppl Table 3 Sect 1.1: ktr_sBCMA_BM_c, assumed equal to ktr_c_BM
    lktr_ab_c_bm    <- log(0.1);       label("Elranatamab transport rate, central to bone marrow (1/h)")   # Suppl Table 3 Sect 1.1: ktr_c_BM, estimated with clinical data
    lktr_ab_bm_c    <- log(0.3);       label("Elranatamab transport rate, bone marrow to central (1/h)")   # Suppl Table 3 Sect 1.1: ktr_BM_c, estimated with clinical data

    # --- Elranatamab PK, carried unchanged from the preliminary popPK ----
    # Supplementary Table 3 reports CL in mL/h (13.9583) but Q in L/h
    # (0.008416); both are transcribed on the units printed for their own row.
    # CL = 13.9583 mL/h = 0.335 L/day cross-checks against the ELREXFIO USPI
    # value of 0.324 L/day.
    lka      <- fixed(log(6.167e-3));  label("Subcutaneous absorption rate constant (1/h)")         # Suppl Table 3 Sect 1.1: ka, preliminary popPK
    lfdepot  <- fixed(log(0.511));     label("Subcutaneous bioavailability (fraction)")             # Suppl Table 3 Sect 1.1: F, preliminary popPK
    lvc      <- fixed(log(4.25));      label("Central volume of distribution (L)")                  # Suppl Table 3 Sect 1.1: Vc [V1], preliminary popPK
    lvp      <- fixed(log(7.8));       label("Peripheral volume of distribution (L)")               # Suppl Table 3 Sect 1.1: Vp [V2], preliminary popPK
    lq       <- fixed(log(0.008416));  label("Inter-compartmental clearance (L/h)")                 # Suppl Table 3 Sect 1.1: Q, preliminary popPK
    lcl      <- fixed(log(0.0139583)); label("Elimination clearance (L/h; reported as 13.9583 mL/h)")  # Suppl Table 3 Sect 1.1: CL, preliminary popPK; kel = CL/Vc

    # --- T cells and cytokines in the central compartment ---------------
    lktc_c_bm <- log(0.2);        label("T-cell transport rate, central to bone marrow (1/h)")      # Suppl Table 3 Sect 1.2: k_TC_c_BM, estimated with clinical data and literature
    lktc_bm_c <- log(0.3893);     label("T-cell transport rate, bone marrow to central (1/h)")      # Suppl Table 3 Sect 1.2: k_TC_BM_c, estimated with clinical data and literature
    kel_tc    <- fixed(0.1046);   label("T-cell elimination rate in central (1/h)")                 # Suppl Table 3 Sect 1.2: kel_T, fixed from literature [1] (Betts 2019)
    lbeta_prod <- log(1);         label("Basal cytokine production rate (pg/mL/h)")                 # Suppl Table 3 Sect 1.2: beta_prod, estimated from literature [14]
    lkdeg_cyt  <- log(0.5);       label("Cytokine degradation rate (1/h)")                          # Suppl Table 3 Sect 1.2: kdeg_cyt, estimated with clinical data and literature
    lktr_cyt   <- log(0.625);     label("Cytokine transit rate, bone marrow chain to central (1/h)")  # Suppl Table 3 Sect 1.2: ktr_cyt, estimated with clinical data and literature

    # --- Paraprotein turnover in the central compartment ----------------
    kdeg_mprotein <- fixed(0.0021); label("Degradation rate of serum M-protein (1/h)")              # Suppl Table 3 Sect 1.3: kdeg_Mp, fixed from literature [12] (Mills 2017)
    kdeg_flc      <- fixed(0.1733); label("Degradation rate of serum free light chain (1/h)")       # Suppl Table 3 Sect 1.3: kdeg_FLC, fixed from literature [13] (Tosi 2013)

    # --- Bone marrow / tumour compartment --------------------------------
    vbm         <- fixed(1.75);       label("Volume of the bone marrow compartment (L)")            # Suppl Table 3 Sect 2.2: V_BM, fixed from literature [11,10]
    lkg0        <- log(3.32e-4);      label("Exponential-phase myeloma growth rate (1/h)")          # Suppl Table 3 Sect 2.2: kg0, calibrated with clinical data and literature
    lkg1        <- log(1500);         label("Linear-phase myeloma growth rate (cells/uL/h)")        # Suppl Table 3 Sect 2.2: kg1, calibrated with clinical data and literature
    psi_switch  <- fixed(20);         label("Exponential-to-linear growth switch exponent (unitless)")  # Suppl Table 3 Sect 2.2: psi, fixed from literature [9] (Simeoni 2004)
    lalpha_kill <- log(0.02);         label("Maximum trimer-driven myeloma kill rate (1/h)")        # Suppl Table 3 Sect 2.2: alpha_kill, calibrated with clinical data
    lkm_kill    <- log(4.25e-4);      label("Trimer:myeloma-cell ratio at half-maximum kill (pM per cells/uL)")  # Suppl Table 3 Sect 2.2: km, estimated with clinical data and literature
    ln_kill     <- log(0.8);          label("Hill coefficient of the trimer-driven kill term (unitless)")  # Suppl Table 3 Sect 2.2: n_kill, calibrated with clinical data
    ltau_mm     <- log(24);           label("Transit time through the myeloma-cell death chain (h)")  # Suppl Table 3 Sect 2.2: tau_M [tau TA], calibrated with clinical data and literature
    lmm_max     <- log(3e6);          label("Maximum myeloma-cell burden in bone marrow (cells/uL)")  # Suppl Table 3 Sect 2.2: MM_max, estimated with clinical data
    lalpha_resis <- log(0.5);         label("Maximum multiplicative increase in half-maximal kill from resistance (unitless)")  # Suppl Table 3 Sect 2.2: alpha_resis, calibrated with clinical data
    lbeta_resis  <- log(0.1);         label("Rate at which resistance reaches its maximum effect (1/h)")  # Suppl Table 3 Sect 2.2: beta_resis, calibrated with clinical data
    ltau_resis   <- log(600);         label("Time at which resistance begins (h)")                  # Suppl Table 3 Sect 2.2: tau_resis, estimated with clinical data
    lbcma_density <- log(12590);      label("Membrane BCMA receptor density (receptors/cell)")             # Suppl Table 3 Sect 2.2: BCMA density, calibrated with clinical data and literature [3]
    lcd3_density  <- log(6e4);        label("CD3 receptor density (receptors/cell)")                             # Suppl Table 3 Sect 2.3: CD3 density, calibrated with clinical data and literature [5]

    # --- Cytokine release in bone marrow (Chen 2019 framework, ref [2]) ---
    cr_imax  <- fixed(1);       label("Maximum inhibition of cytokine release by cumulative exposure (fraction)")  # Suppl Table 3 Sect 2.3: Imax, fixed from literature [2]
    cr_n_ih  <- fixed(2.5);     label("Hill coefficient for cytokine-release inhibition (unitless)")   # Suppl Table 3 Sect 2.3: n1C, fixed from literature [2]
    cr_ic50  <- fixed(10000);   label("Cumulative cytokine exposure giving 50% release inhibition (pg/mL*h)")  # Suppl Table 3 Sect 2.3: IC50, from literature [2]
    cr_emax  <- fixed(80590);   label("Maximum cytokine release rate (pg/mL/h)")                    # Suppl Table 3 Sect 2.3: Emax, from literature [2]
    cr_ec50  <- fixed(0.5);     label("Trimer concentration giving 50% of maximum cytokine release (pM)")  # Suppl Table 3 Sect 2.3: EC50, from literature [2]
    cr_n     <- fixed(1);       label("Hill coefficient for cytokine release (unitless)")           # Suppl Table 3 Sect 2.3: n2C, from literature [2]

    # --- Cytokine- and trimer-dependent T-cell retention in bone marrow ---
    # Supplementary Table 3 lists Imax_TC as "Estimated" with a value of
    # exactly 1, the upper bound of a maximum-inhibition fraction; it is
    # therefore encoded as fixed rather than log-transformed (see vignette
    # Assumptions and deviations).
    imax_tc          <- fixed(1);  label("Maximum inhibition of T-cell egress from bone marrow (fraction)")  # Suppl Table 3 Sect 2.3: Imax_TC, estimated (value at the bound)
    lksat_tc         <- log(150);  label("Cytokine concentration at 50% of maximum T-cell egress inhibition (pg/mL)")  # Suppl Table 3 Sect 2.3: K1, estimated
    lksat_trimer_tc  <- log(0.02); label("Trimer concentration at 50% of maximum T-cell egress inhibition (pM)")  # Suppl Table 3 Sect 2.3: K2, estimated
    ln_cytokine_tc   <- log(2);    label("Hill coefficient of the cytokine arm of T-cell egress inhibition (unitless)")  # Suppl Table 3 Sect 2.3: n_cyt, estimated
    ln_trimer_tc     <- log(1.5);  label("Hill coefficient of the trimer arm of T-cell egress inhibition (unitless)")  # Suppl Table 3 Sect 2.3: n_tri, estimated

    # --- Numerical guards ------------------------------------------------
    # Supplementary Sect 1.2.2 introduces epsilon as "a small value introduced
    # for numerical stability" in the trimer:myeloma-cell ratio and gives no
    # value. It is a solver guard, not a calibrated quantity; the default is
    # 12 orders of magnitude below the smallest baseline myeloma burden of
    # Supplementary Table 2. See vignette Assumptions and deviations.
    eps_kill <- fixed(1e-6); label("Numerical stabiliser in the trimer:myeloma-cell ratio (cells/uL)")  # Suppl Sect 1.2.2, value not reported

    # --- Baseline states (Supplementary Table 2 sampling ranges) ---------
    # The paper publishes sampling RANGES, not point estimates: every result
    # comes from a 120-patient virtual population with independently sampled
    # baselines, so the model has no published "typical" patient. Each default
    # below is the geometric midpoint of the published range (the distributions
    # in Figure 2b are strongly right-skewed, which rules out the arithmetic
    # midpoint), and every one is flagged in the vignette Assumptions and
    # deviations. Vary them to reproduce the paper's virtual populations.
    lsbcma0_central   <- log(5329);   label("Baseline free soluble BCMA in central (pM; published range 185-153520 = 1-830 ng/mL)")   # Suppl Table 2: sBCMA_c range; default = geometric midpoint, no point estimate published
    ltcell0_central   <- log(203);    label("Baseline CD3+ T cells in central (cells/uL; published range 16-2580)")                   # Suppl Table 2: Tc_c range; default = geometric midpoint, no point estimate published
    lmm0              <- log(2.25e5); label("Baseline myeloma-cell burden in bone marrow (cells/uL; published range 0-450000)")       # Suppl Table 2: MM range; default = arithmetic midpoint (range starts at 0, so no geometric midpoint exists)
    lmprotein0        <- log(18.7);   label("Baseline serum M-protein (g/L; published range 0-70)")                                   # Suppl Table 2: M_P range; default = geometric midpoint of 5-70 g/L, 5 g/L being the IMWG measurability floor the paper cites
    lflc0             <- log(400.9);  label("Baseline involved serum free light chain (mg/L; published range 0.92-174680)")           # Suppl Table 2: FLC range; default = geometric midpoint, no point estimate published
    cytokine0_central <- fixed(2.857); label("Baseline cytokine (IL-6) concentration in central (pg/mL)")                             # Suppl Table 2: C_c initial value, reported explicitly

    # --- Unit conversion --------------------------------------------------
    # Supplementary Eqs 8, 10, 20 and 21 convert (receptors/cell x cells/uL)
    # to pM as `L_to_uL / av_num * 1e12`, i.e. 1e6 uL/L / Avogadro x 1e12.
    sbcma_pm_per_ngml <- fixed(185); label("Soluble BCMA scale: pM per ng/mL")  # derived from the paper's own two scales: Suppl Table 2 range 185-153520 pM equals the 1-830 ng/mL plausible-patient range of Figure 2b

    # --- Residual error ---------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # the source is a simulation-only QSP model and reports no residual-error model
  })

  model({
    # =================================================================
    # 1. Unit conversions and derived constants
    # =================================================================
    l_to_ul <- 1e6                       # uL per L
    av_num  <- 6.02214076e23             # Avogadro constant, 1/mol (SI defined)
    # receptors/cell x cells/uL -> pM (Suppl Eqs 8, 20, 21)
    recpt_to_pm <- l_to_ul / av_num * 1e12

    # =================================================================
    # 2. Individual parameters
    # =================================================================
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    cl     <- exp(lcl)

    ktr_ab_c_bm <- exp(lktr_ab_c_bm)
    ktr_ab_bm_c <- exp(lktr_ab_bm_c)
    ktc_c_bm    <- exp(lktc_c_bm)
    ktc_bm_c    <- exp(lktc_bm_c)
    beta_prod   <- exp(lbeta_prod)
    kdeg_cyt    <- exp(lkdeg_cyt)
    ktr_cyt     <- exp(lktr_cyt)

    kg0         <- exp(lkg0)
    kg1         <- exp(lkg1)
    alpha_kill  <- exp(lalpha_kill)
    km_kill     <- exp(lkm_kill)
    n_kill      <- exp(ln_kill)
    tau_mm      <- exp(ltau_mm)
    mm_max      <- exp(lmm_max)
    alpha_resis <- exp(lalpha_resis)
    beta_resis  <- exp(lbeta_resis)
    tau_resis   <- exp(ltau_resis)
    bcma_density <- exp(lbcma_density)
    cd3_density  <- exp(lcd3_density)

    ksat_tc        <- exp(lksat_tc)
    ksat_trimer_tc <- exp(lksat_trimer_tc)
    n_cytokine_tc  <- exp(ln_cytokine_tc)
    n_trimer_tc    <- exp(ln_trimer_tc)

    sbcma0_central <- exp(lsbcma0_central)
    tcell0_central <- exp(ltcell0_central)
    mm0            <- exp(lmm0)
    mprotein0      <- exp(lmprotein0)
    flc0           <- exp(lflc0)

    # =================================================================
    # 3. Micro-constants. Suppl Eq 1 writes the peripheral return as
    #    k21 * (Vp/Vc) * [Ab_p] with k21 = Q/Vp, i.e. Q/Vc * [Ab_p];
    #    the standard concentration-based two-compartment form.
    # =================================================================
    kel <- cl / vc                       # Suppl Table 3: kel = CL/Vc
    k12 <- q / vc
    k21 <- q / vp

    # =================================================================
    # 4. Steady-state baselines in the bone marrow.
    #    Supplementary Sects 1.1.1 and 1.2.1 both assert that the system is
    #    at steady state at t = 0. Supplementary Table 2 lists the marrow
    #    T-cell and marrow sBCMA states as separately sampled, but their ODEs
    #    (Eqs 14 and 28) are only stationary at t = 0 when the marrow value
    #    equals the central value scaled by the transport ratio below. The
    #    steady-state relations are used so an untreated simulation holds its
    #    baseline; see vignette Assumptions and deviations.
    # =================================================================
    tcell0_tumor <- ktc_c_bm * vc / (ktc_bm_c * vbm) * tcell0_central          # from Suppl Eq 28 at t = 0
    sbcma0_tumor <- kdeg_sbcma * vc / (ktr_sbcma_bm_c * vbm) * sbcma0_central  # from Suppl Eq 14 at t = 0

    # =================================================================
    # 5. Receptor pools (Suppl Eqs 8-9, 20-23)
    # =================================================================
    mm_tot <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    cd3_tot_central  <- tcell_central * cd3_density * recpt_to_pm             # Suppl Eq 8
    cd3_free_central <- cd3_tot_central - drug_cd3_central                    # Suppl Eq 9

    cd3_tot_tumor  <- tcell_tumor * cd3_density * recpt_to_pm                 # Suppl Eq 20
    bcma_tot_tumor <- mm_tot * bcma_density * recpt_to_pm                     # Suppl Eq 21
    bcma_free_tumor <- bcma_tot_tumor - drug_bcma_tumor - trimer              # Suppl Eq 22
    cd3_free_tumor  <- cd3_tot_tumor - drug_cd3_tumor - trimer                # Suppl Eq 23

    # =================================================================
    # 6. Trimer-driven kill with time-dependent resistance
    #    (Suppl Sect 1.2.2). The printed resistance equation
    #    km_kill * (1 + alpha_resis * (1 - exp(-beta_resis * (t - tau_resis))))
    #    diverges to large NEGATIVE half-maximal values for t < tau_resis
    #    (exp(+0.1 * 600) ~ 1e26), which is plainly not the stated intent
    #    ("tau_resis is the time it takes for the resistance to occur").
    #    The elapsed time is therefore floored at zero so resistance is
    #    absent until t = tau_resis and then rises to alpha_resis. Flagged
    #    in the vignette Assumptions and deviations.
    # =================================================================
    resis_time  <- max(0, t - tau_resis)
    km_kill_t   <- km_kill * (1 + alpha_resis * (1 - exp(-beta_resis * resis_time)))

    trimer_pos   <- max(trimer, 0)
    trimer_ratio <- trimer_pos / (eps_kill + cycling_cells)
    kmkill_rate  <- alpha_kill * trimer_ratio^n_kill /
      (km_kill_t^n_kill + trimer_ratio^n_kill)

    # =================================================================
    # 7. Cytokine release and T-cell retention (Suppl Eqs 29, 38)
    # =================================================================
    cyt_bar     <- (cytokine_tumor + cytokine_transit1 + cytokine_transit2 +
                      cytokine_transit3 + cytokine_transit4 + cytokine_transit5) / 6
    cyt_bar_pos <- max(cyt_bar, 0)
    cauc_pos    <- max(cytokine_auc, 0)

    t_inh <- imax_tc *
      (cyt_bar_pos^n_cytokine_tc / (ksat_tc^n_cytokine_tc + cyt_bar_pos^n_cytokine_tc)) *
      (trimer_pos^n_trimer_tc / (ksat_trimer_tc^n_trimer_tc + trimer_pos^n_trimer_tc))   # Suppl Eq 29

    ih    <- cr_imax * cauc_pos^cr_n_ih / (cr_ic50^cr_n_ih + cauc_pos^cr_n_ih)           # Suppl Eq 38
    r_syn <- cr_emax * trimer_pos^cr_n / (cr_ec50^cr_n + trimer_pos^cr_n)                # Suppl Eq 38

    # =================================================================
    # 8. ODE system
    # =================================================================
    # --- Central compartment (Suppl Eqs 1-5) -------------------------
    d/dt(depot) <- -ka * depot                                                # amount, pmol

    d/dt(central) <- ka * depot / vc - k12 * central + k21 * vp / vc * peripheral1 -
      kon_bcma * central * sbcma_central + koff_bcma * drug_sbcma_central -
      kon_cd3 * central * cd3_free_central + koff_cd3 * drug_cd3_central -
      ktr_ab_c_bm * central + ktr_ab_bm_c * vbm / vc * tumor -
      kel * central                                                           # Suppl Eq 1

    d/dt(peripheral1) <- k12 * vc / vp * central - k21 * peripheral1          # implied by Suppl Eq 1 (not printed separately)

    d/dt(sbcma_central) <-
      kdeg_sbcma * (sbcma0_central / sbcma0_tumor * sbcma_tumor - sbcma_central) -
      kon_bcma * central * sbcma_central + koff_bcma * drug_sbcma_central     # Suppl Eqs 2-3

    d/dt(drug_sbcma_central) <-
      kon_bcma * central * sbcma_central - koff_bcma * drug_sbcma_central +
      ktr_sbcma_bm_c * vbm / vc * drug_sbcma_tumor - kel * drug_sbcma_central # Suppl Eq 4

    d/dt(drug_cd3_central) <-
      kon_cd3 * central * cd3_free_central - koff_cd3 * drug_cd3_central      # Suppl Eq 5

    d/dt(tcell_central) <-
      kel_tc * tcell0_central - (kel_tc + ktc_c_bm) * tcell_central +
      ktc_bm_c * tcell_tumor * vbm / vc * (1 - t_inh)                         # Suppl Eq 6

    d/dt(cytokine_central) <-
      beta_prod + ktr_cyt * cytokine_transit5 * vbm / vc -
      kdeg_cyt * cytokine_central                                             # Suppl Eq 7

    # --- Paraproteins (Suppl Eqs 11-12) ------------------------------
    d/dt(mprotein) <- kdeg_mprotein * (mprotein0 / mm0) * mm_tot - kdeg_mprotein * mprotein  # Suppl Eq 11
    d/dt(flc)      <- kdeg_flc * (flc0 / mm0) * mm_tot - kdeg_flc * flc                      # Suppl Eq 12

    # --- Bone marrow drug and sBCMA (Suppl Eqs 13-14) ----------------
    d/dt(tumor) <-
      koff_bcma * drug_sbcma_tumor - kon_bcma * tumor * sbcma_tumor +
      koff_cd3 * drug_cd3_tumor - kon_cd3 * tumor * cd3_free_tumor +
      koff_bcma * drug_bcma_tumor - kon_bcma * tumor * bcma_free_tumor -
      ktr_ab_bm_c * tumor + ktr_ab_c_bm * vc / vbm * central                  # Suppl Eq 13

    d/dt(sbcma_tumor) <-
      kdeg_sbcma * vc / vbm * (sbcma0_central / mm0) * mm_tot -
      ktr_sbcma_bm_c * sbcma_tumor -
      kon_bcma * tumor * sbcma_tumor + koff_bcma * drug_sbcma_tumor           # Suppl Eq 14

    # --- Dimers and trimer in bone marrow (Suppl Eqs 15-19) ----------
    d/dt(drug_bcma_tumor) <-
      kon_bcma * tumor * bcma_free_tumor - koff_bcma * drug_bcma_tumor -
      kon_cd3 * drug_bcma_tumor * cd3_free_tumor + koff_cd3 * trimer          # Suppl Eq 15

    # Suppl Eq 17 puts the factor VBM/Vc on the loss of the marrow drug-sBCMA
    # dimer as well as on its gain in central (Eq 4). The free-drug pair
    # (Eqs 1 and 13) carries the factor on the gain term only, so the printed
    # dimer pair is not mass-conserving. Transcribed as printed; flagged in
    # the vignette Assumptions and deviations.
    d/dt(drug_sbcma_tumor) <-
      kon_bcma * tumor * sbcma_tumor - koff_bcma * drug_sbcma_tumor -
      ktr_sbcma_bm_c * vbm / vc * drug_sbcma_tumor                            # Suppl Eqs 16-17

    d/dt(drug_cd3_tumor) <-
      kon_cd3 * tumor * cd3_free_tumor - koff_cd3 * drug_cd3_tumor -
      kon_bcma * drug_cd3_tumor * bcma_free_tumor + koff_bcma * trimer        # Suppl Eq 18

    d/dt(trimer) <-
      kon_bcma * drug_cd3_tumor * bcma_free_tumor - koff_bcma * trimer +
      kon_cd3 * drug_bcma_tumor * cd3_free_tumor - koff_cd3 * trimer          # Suppl Eq 19

    # --- Myeloma cells: Simeoni growth with transduction chain --------
    growth_num <- kg0 * (1 - mm_tot / mm_max)
    growth_den <- (1 + (kg0 / kg1 * mm_tot)^psi_switch)^(1 / psi_switch)
    d/dt(cycling_cells)  <- growth_num / growth_den * cycling_cells - kmkill_rate * cycling_cells  # Suppl Eq 24
    d/dt(damaged_cells1) <- kmkill_rate * cycling_cells - 4 / tau_mm * damaged_cells1              # Suppl Eq 25
    d/dt(damaged_cells2) <- 4 / tau_mm * (damaged_cells1 - damaged_cells2)                         # Suppl Eq 26
    d/dt(damaged_cells3) <- 4 / tau_mm * (damaged_cells2 - damaged_cells3)                         # Suppl Eq 27

    # --- T cells and cytokines in bone marrow (Suppl Eqs 28, 30-36) ---
    d/dt(tcell_tumor) <-
      ktc_c_bm * vc / vbm * tcell_central - ktc_bm_c * tcell_tumor * (1 - t_inh)  # Suppl Eq 28

    d/dt(cytokine_tumor)    <- r_syn * (1 - ih) - (ktr_cyt + kdeg_cyt) * cytokine_tumor              # Suppl Eq 30
    d/dt(cytokine_transit1) <- ktr_cyt * cytokine_tumor    - (ktr_cyt + kdeg_cyt) * cytokine_transit1  # Suppl Eq 31
    d/dt(cytokine_transit2) <- ktr_cyt * cytokine_transit1 - (ktr_cyt + kdeg_cyt) * cytokine_transit2  # Suppl Eq 32
    d/dt(cytokine_transit3) <- ktr_cyt * cytokine_transit2 - (ktr_cyt + kdeg_cyt) * cytokine_transit3  # Suppl Eq 33
    d/dt(cytokine_transit4) <- ktr_cyt * cytokine_transit3 - (ktr_cyt + kdeg_cyt) * cytokine_transit4  # Suppl Eq 34
    d/dt(cytokine_transit5) <- ktr_cyt * cytokine_transit4 - (ktr_cyt + kdeg_cyt) * cytokine_transit5  # Suppl Eq 35
    d/dt(cytokine_auc)      <- cytokine_tumor - beta_prod                                             # Suppl Eq 36

    # =================================================================
    # 9. Dose-level adjustments
    # =================================================================
    f(depot) <- fdepot

    # Supplementary Eq 37 rescales the cumulative cytokine exposure at the
    # start of each dosing interval:
    #
    #   Cauc_N(0) = Cauc_(N-1)(tau) / ( 5 * (1 - (N+1)^2 / (1.3^2 + (N+1)^2)) )
    #
    # with N the number of doses given. This is a DISCRETE state reset, which
    # an rxode2 ODE cannot express, so it is delivered from the event table as
    # an `evid = 6` (multiply) record on `cytokine_auc` at each drug-dose time,
    # carrying the reciprocal of that divisor as its `amt`.
    # nlmixr2lib::Poels_2025_elranatamab_qsp_events() builds a compliant event
    # table and is the supported way to simulate this model; the multiplier is
    # computed there rather than here so that it never depends on rxode2's
    # internal ordering of same-time dose records. Simulating WITHOUT those
    # records simply omits the dose-to-dose attenuation of cytokine release.

    # =================================================================
    # 10. Initial conditions (Suppl Table 2 and Suppl Eq 10)
    # =================================================================
    depot(0)              <- 0
    central(0)            <- 0
    peripheral1(0)        <- 0
    tumor(0)              <- 0
    sbcma_central(0)      <- sbcma0_central
    sbcma_tumor(0)        <- sbcma0_tumor
    drug_sbcma_central(0) <- 0
    drug_cd3_central(0)   <- 0
    drug_bcma_tumor(0)    <- 0
    drug_sbcma_tumor(0)   <- 0
    drug_cd3_tumor(0)     <- 0
    trimer(0)             <- 0
    tcell_central(0)      <- tcell0_central
    tcell_tumor(0)        <- tcell0_tumor
    cycling_cells(0)      <- mm0
    damaged_cells1(0)     <- 0
    damaged_cells2(0)     <- 0
    damaged_cells3(0)     <- 0
    mprotein(0)           <- mprotein0
    flc(0)                <- flc0
    cytokine_tumor(0)     <- 0
    cytokine_transit1(0)  <- 0
    cytokine_transit2(0)  <- 0
    cytokine_transit3(0)  <- 0
    cytokine_transit4(0)  <- 0
    cytokine_transit5(0)  <- 0
    cytokine_central(0)   <- cytokine0_central
    cytokine_auc(0)       <- 0

    # =================================================================
    # 11. Observations
    # =================================================================
    Cc          <- central                                                   # free elranatamab in serum (pM)
    CcTotal     <- central + drug_sbcma_central + drug_cd3_central           # total-analyte assay equivalent (pM)
    sbcmaFree   <- sbcma_central / sbcma_pm_per_ngml                         # free serum sBCMA (ng/mL)
    mProtein    <- mprotein                                                  # serum M-protein (g/L)
    flcSerum    <- flc                                                       # involved serum FLC (mg/L)
    cytokineSerum <- cytokine_central                                        # serum IL-6-like cytokine (pg/mL)
    tumorBurden <- mm_tot                                                    # total myeloma burden in marrow (cells/uL)
    # "Effective binding ratio" of the Results text: trimers per BCMA receptor.
    bindingRatio <- trimer / (bcma_tot_tumor + eps_kill)
    # R_trimer of Supplementary Eq 39: trimer per BCMA-bound drug complex.
    trimerFraction <- trimer / (trimer + drug_bcma_tumor + eps_kill)

    Cc ~ prop(propSd)
  })
}
