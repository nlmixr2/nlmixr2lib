Zheng_2014_azithromycin <- function() {
  description <- "Semi-mechanistic tissue distribution population PK model for oral azithromycin in healthy adults (Zheng 2014). Three-compartment plasma PK (depot with absorption lag time and first-order absorption, central, two peripheral compartments) with concentration-dependent fraction unbound in plasma (equation 1). Three tissue distribution compartments (muscle interstitial space fluid, subcutaneous adipose tissue interstitial space fluid, polymorphonuclear-leukocyte (PML) cytosol) each driven by free unbound (or, for PML cytosol, free unionized) plasma drug via first-order rate constants kin and kout, with tissue-specific distribution factors df_muscle, df_adipose, df_pmn that scale the steady-state tissue:plasma free-unbound ratio. Each tissue compartment also exchanges with a deep nonspecific phospholipid-binding compartment via shared kon and koff (Methods equations 1-13)."
  reference <- "Zheng S, Matzneller P, Zeitlinger M, Schmidt S. Development of a population pharmacokinetic model characterizing the tissue distribution of azithromycin in healthy subjects. Antimicrob Agents Chemother. 2014 Nov;58(11):6675-6684. doi:10.1128/AAC.02904-14"
  vignette <- "Zheng_2014_azithromycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Reported in Methods (77.68 +/- 8.56 kg) but no covariate effect was retained: 'Due to the limited number of subjects (n = 6) and the relatively homogeneous study population ... a rigorous covariate analysis was deemed meaningless for this study.'"
    ),
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Reported in Methods (29.0 +/- 9.63 years) but not screened or retained for the same reason as WT (homogeneous cohort)."
    ),
    HT = list(
      description = "Body height.",
      units       = "cm",
      type        = "continuous",
      notes       = "Reported in Methods (184.17 +/- 6.74 cm) but not screened or retained for the same reason as WT (homogeneous cohort)."
    ),
    BMI = list(
      description = "Body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Reported in Methods (22.83 +/- 1.39 kg/m^2) but not screened or retained for the same reason as WT (homogeneous cohort)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Fixed at 0 for this cohort: all 6 enrolled subjects were healthy adult males. No covariate effect was screened."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 6,
    n_studies       = 1,
    age_range       = "29.0 +/- 9.63 years (mean +/- SD)",
    weight_range    = "77.68 +/- 8.56 kg (mean +/- SD)",
    height_range    = "184.17 +/- 6.74 cm (mean +/- SD)",
    bmi_range       = "22.83 +/- 1.39 kg/m^2 (mean +/- SD)",
    sex_female_pct  = 0,
    disease_state   = "Healthy male volunteers.",
    dose_range      = "500 mg once daily oral for 3 days.",
    regions         = "Austria (Medical University of Vienna).",
    n_observations  = "Total plasma azithromycin sampled at 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 6, 8 h on days 1 and 3, and three timepoints each on days 5 and 10. White-blood-cell total concentrations sampled at 0, 2, 6, 10 h on days 1 and 3, and one timepoint each on days 5 and 10. Muscle ISF and subcutis ISF unbound concentrations sampled via microdialysis at prespecified timepoints on the same study days (Matzneller et al. 2013, ref. 4).",
    notes           = "All demographics from Zheng 2014 Methods 'Subjects and pharmacokinetic study' (which cites Matzneller et al. 2013 for the upstream clinical microdialysis dataset). Underlying clinical study: Matzneller P, Krasniqi S, Kinzig M, Sorgel F, Huttner S, Lackner E, Muller M, Zeitlinger M. Pharmacokinetics of polymorphonuclear leukocyte and plasma concentrations of azithromycin in healthy volunteers. Antimicrob Agents Chemother. 2013;57(4):1736-42 (ref. 4)."
  )

  paper_specific_compartments <- c(
    "muscle_deep",
    "adipose_deep",
    "pmn_deep"
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma 3-compartment PK parameters
    # ------------------------------------------------------------------
    # Tlag, ka, and CL/F were FIXED at the values shown in Table 1 ('fixed')
    # during fitting of the full tissue distribution model. Footnote c to
    # Table 1: 'The bootstrap results for the Tlag, ka, and CL/F fixed
    # parameter estimates were based on bootstrap (n=1,000) for unbound
    # plasma only.' i.e., the bootstrap re-estimated them in an unbound-
    # plasma-only model, but the tissue model held them fixed.
    ltlag <- fixed(log(1.45))
    label("Absorption lag time (h, FIXED)")                                  # Zheng 2014 Table 1: 1.45 h (fixed)

    lka <- fixed(log(0.88))
    label("First-order absorption rate constant (1/h, FIXED)")               # Zheng 2014 Table 1: 0.88 1/h (fixed)

    lcl <- fixed(log(258))
    label("Apparent oral clearance CL/F (L/h, FIXED)")                       # Zheng 2014 Table 1: 258 L/h (fixed)

    lvc <- log(160)
    label("Apparent central volume Vc/F (L)")                                # Zheng 2014 Table 1: 160 L

    lvp <- log(1190)
    label("Apparent fast peripheral volume Vp1/F (L)")                       # Zheng 2014 Table 1: 1190 L

    lq <- log(207)
    label("Apparent fast inter-compartmental clearance Qp1/F (L/h)")         # Zheng 2014 Table 1: 207 L/h

    lvp2 <- log(9721)
    label("Apparent slow peripheral volume Vp2/F (L)")                       # Zheng 2014 Table 1: 9721 L

    lq2 <- log(101)
    label("Apparent slow inter-compartmental clearance Qp2/F (L/h)")         # Zheng 2014 Table 1: 101 L/h

    # ------------------------------------------------------------------
    # Saturable plasma protein binding (Equation 1, derived from
    # Bouvier d'Yvoire, Dresco, Tulkens 2001 (ref. 20) by Zheng et al.
    # via digitisation in GraphPad Prism 5). All three constants are
    # FIXED structural inputs to the Zheng 2014 model.
    #   fu_p(Cp) = fu_base + fu_emax * Cp / (fu_ec50 + Cp)
    # Note: fu is an INCREASING function of total plasma concentration
    # because the high-affinity binding sites saturate at higher Cp.
    # ------------------------------------------------------------------
    fu_base <- fixed(0.4984)
    label("Baseline fraction unbound in plasma at low Cp (unitless, FIXED)")  # Zheng 2014 Equation 1

    fu_emax <- fixed(0.5339)
    label("Maximum incremental fraction unbound at high Cp (unitless, FIXED)") # Zheng 2014 Equation 1

    fu_ec50 <- fixed(230.9)
    label("Cp at which fu_emax is half-saturated (ng/mL, FIXED)")             # Zheng 2014 Equation 1

    # ------------------------------------------------------------------
    # Unionized fractions (Equation 6) at the cited compartment pH.
    # f_unionized = 1 / (1 + 10^(pKa1 - pH) + 10^(pKa1 - pH) + 10^(pKa2 - pH))
    # AZM is a diprotic base with pKa1 = 8.1 and pKa2 = 8.8.
    # ------------------------------------------------------------------
    funi_plasma <- fixed(0.0076)
    label("Fraction unionized in plasma at pH 7.4 (unitless, FIXED)")          # Zheng 2014 Results: f_unionized in plasma/ISF = 0.0076

    funi_pmn_cyto <- fixed(0.0012)
    label("Fraction unionized in PML cytosol at pH ~7 (unitless, FIXED)")      # Zheng 2014 Results: f_unionized in PML cytosol = 0.0012

    # ------------------------------------------------------------------
    # Tissue distribution rate constants (Equation 2). kin and kout are
    # SHARED across the three tissues; kon and koff (nonspecific deep
    # binding) are also shared across the three tissues per Methods:
    # 'with the same kon (0.56 h-1) and koff (0.05 h-1) values for the
    # different tissue sites'.
    # ------------------------------------------------------------------
    lkin <- log(0.16)
    label("Plasma-to-tissue first-order distribution rate constant (1/h)")     # Zheng 2014 Table 1: 0.16 1/h

    lkout <- log(0.15)
    label("Tissue-to-plasma first-order distribution rate constant (1/h)")     # Zheng 2014 Table 1: 0.15 1/h

    lkon <- log(0.56)
    label("Nonspecific binding on rate (1/h, shared across muscle/adipose/PML)") # Zheng 2014 Table 1: 0.56 1/h

    lkoff <- log(0.05)
    label("Nonspecific binding off rate (1/h, shared across muscle/adipose/PML)") # Zheng 2014 Table 1: 0.05 1/h

    # ------------------------------------------------------------------
    # Tissue distribution factors (Equations 4-5, 7). One per tissue.
    # df = ratio of free-unbound (free-unionized for PML cytosol)
    # concentration in tissue vs plasma at steady state, modulo the
    # kin/kout ratio.
    # ------------------------------------------------------------------
    ldf_muscle <- log(0.55)
    label("Tissue distribution factor for muscle ISF (unitless)")              # Zheng 2014 Table 1: DF_muscle = 0.55

    ldf_adipose <- log(0.25)
    label("Tissue distribution factor for subcutis (subcutaneous adipose) ISF (unitless)") # Zheng 2014 Table 1: DF_subcutis = 0.25

    ldf_pmn <- log(52)
    label("Tissue distribution factor for PML cytosol unionized (unitless)")   # Zheng 2014 Table 1: DF_PML(cytosol) = 52

    # ------------------------------------------------------------------
    # Interindividual variability (Zheng 2014 Table 1 'Interindividual
    # variability (%)'). Exponential IIV model (Equation 12). Stored as
    # variance on the log scale: omega^2 = log((CV/100)^2 + 1).
    # The very small values for etalkin and etaldf_pmn (0.22%) likely
    # reflect convergence to a boundary in the single-run fit; the
    # nonparametric bootstrap (medians 22.4% and 9.0% respectively;
    # Table 1) supports more meaningful IIV but the model file uses
    # the single-run point estimates that the paper reports.
    # ------------------------------------------------------------------
    etaltlag      ~ 0.030559   # 17.6%  CV (Tlag);              log(0.176^2 + 1) = 0.030559
    etalcl        ~ 0.082382   # 29.3%  CV (CL/F);              log(0.293^2 + 1) = 0.082382
    etalvc        ~ 1.343243   # 168.3% CV (V_central/F);       log(1.683^2 + 1) = 1.343243
    etalkin       ~ 4.84e-06   # 0.22%  CV (kin);               log(0.0022^2 + 1) = 4.84e-06 (boundary)
    etaldf_muscle ~ 0.069794   # 26.9%  CV (DF_muscle);         log(0.269^2 + 1) = 0.069794
    etaldf_adipose ~ 0.094540  # 31.5%  CV (DF_subcutis);       log(0.315^2 + 1) = 0.094540
    etaldf_pmn    ~ 4.84e-06   # 0.22%  CV (DF_PML(cytosol));   log(0.0022^2 + 1) = 4.84e-06 (boundary)

    # ------------------------------------------------------------------
    # Residual error (Zheng 2014 Table 1 'Residual variability').
    # Equation 13: y_obs = y_pred * (1 + eps_prop) + eps_add, i.e.,
    # combined additive + proportional. The paper reports per-output
    # additive in ng/mL (e.g., plasma 35.2 ng/mL); the additive values
    # below are CONVERTED to mg/L (divide by 1000) to match the
    # concentration units declared at the top of this file. The two
    # near-zero additive values (1e-6 ng/mL = 1e-9 mg/L) represent
    # essentially pure-proportional fits in the paper.
    # ------------------------------------------------------------------
    propSd          <- 0.14
    label("Proportional residual error on free plasma Cc (fraction)")            # Zheng 2014 Table 1: Plasma proportional = 0.14
    addSd           <- 0.0352
    label("Additive residual error on free plasma Cc (mg/L)")                    # Zheng 2014 Table 1: Plasma additive = 35.2 ng/mL = 0.0352 mg/L

    propSd_Cmuscle  <- 0.14
    label("Proportional residual error on muscle ISF (fraction)")                # Zheng 2014 Table 1: Muscle ISF proportional = 0.14
    addSd_Cmuscle   <- 0.00051
    label("Additive residual error on muscle ISF (mg/L)")                        # Zheng 2014 Table 1: Muscle ISF additive = 0.51 ng/mL = 0.00051 mg/L

    propSd_Cadipose <- 0.34
    label("Proportional residual error on subcutis ISF (fraction)")              # Zheng 2014 Table 1: Subcutis ISF proportional = 0.34
    addSd_Cadipose  <- 1e-9
    label("Additive residual error on subcutis ISF (mg/L)")                      # Zheng 2014 Table 1: Subcutis ISF additive = 1e-6 ng/mL = 1e-9 mg/L

    propSd_Cpmncyto <- 0.23
    label("Proportional residual error on PML cytosol unionized (fraction)")     # Zheng 2014 Table 1: PML cytosol proportional = 0.23
    addSd_Cpmncyto  <- 1e-9
    label("Additive residual error on PML cytosol unionized (mg/L)")             # Zheng 2014 Table 1: PML cytosol additive = 1e-6 ng/mL = 1e-9 mg/L
  })

  model({
    # ----------------------------------------------------------------
    # Individual structural parameters
    # ----------------------------------------------------------------
    tlag       <- exp(ltlag + etaltlag)
    ka         <- exp(lka)
    cl         <- exp(lcl + etalcl)
    vc         <- exp(lvc + etalvc)
    vp         <- exp(lvp)
    q          <- exp(lq)
    vp2        <- exp(lvp2)
    q2         <- exp(lq2)

    kin        <- exp(lkin + etalkin)
    kout       <- exp(lkout)
    kon        <- exp(lkon)
    koff       <- exp(lkoff)

    df_muscle  <- exp(ldf_muscle  + etaldf_muscle)
    df_adipose <- exp(ldf_adipose + etaldf_adipose)
    df_pmn     <- exp(ldf_pmn     + etaldf_pmn)

    # Apparent tissue volumes are absorbed into df via the steady-state
    # algebra of Equations 3-5 (rearranged for V_tissue):
    #   V_tissue = V_plasma * (kin / kout) / df_tissue
    # This lets the tissue ODE be written in amounts (Equation 2) with
    # kin * (unbound plasma amount) as the influx term, and yields the
    # observed concentration as A_tissue / V_tissue. For PML cytosol the
    # 'plasma amount' is the unionized amount, so DF approximates the
    # tissue:plasma UNIONIZED ratio (Equation 7).
    v_muscle  <- vc * (kin / kout) / df_muscle
    v_adipose <- vc * (kin / kout) / df_adipose
    v_pmn     <- vc * (kin / kout) / df_pmn

    # ----------------------------------------------------------------
    # Saturable plasma protein binding (Equation 1)
    # ----------------------------------------------------------------
    # The paper's Equation 1 reports Cp in ng/mL; the central
    # concentration central/vc here is in mg/L. Convert internally:
    #   Cp_ngml = 1000 * central / vc
    # so fu_ec50 = 230.9 (ng/mL) stays in the paper's native units.
    cp_ngml <- 1000 * central / vc
    fu_p    <- fu_base + fu_emax * cp_ngml / (fu_ec50 + cp_ngml)

    # ----------------------------------------------------------------
    # ODE system
    # Plasma 3-compartment with first-order oral absorption.
    # central holds drug AMOUNT (mg). Inter-compartmental flux uses
    # k_ij = Q_ij / V_origin micro-constants.
    # ----------------------------------------------------------------
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (kel + k12 + k13) * central +
                          k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Absorption lag time on the depot compartment.
    alag(depot) <- tlag

    # Tissue distribution (Equation 2) -- for each tissue:
    #   d/dt(A_tissue) = kin * A_plasma_inflow
    #                  - kout * A_tissue
    #                  + koff * A_tissue_deep
    #                  - kon  * A_tissue
    #   d/dt(A_tissue_deep) = kon * A_tissue - koff * A_tissue_deep
    # For muscle and adipose ISF, A_plasma_inflow = central * fu_p
    # (unbound amount in plasma). For PML cytosol, A_plasma_inflow =
    # central * fu_p * funi_plasma (unionized amount in plasma); the
    # output concentration is the unionized concentration in the PML
    # cytosol (matches Figure 2B 'PML (cytosolic unionized)').
    d/dt(muscle)       <-  kin * fu_p * central               - kout * muscle  + koff * muscle_deep  - kon * muscle
    d/dt(muscle_deep)  <-  kon * muscle  - koff * muscle_deep

    d/dt(adipose)      <-  kin * fu_p * central               - kout * adipose + koff * adipose_deep - kon * adipose
    d/dt(adipose_deep) <-  kon * adipose - koff * adipose_deep

    d/dt(pmn)          <-  kin * fu_p * funi_plasma * central - kout * pmn     + koff * pmn_deep     - kon * pmn
    d/dt(pmn_deep)     <-  kon * pmn     - koff * pmn_deep

    # ----------------------------------------------------------------
    # Observation variables
    # Cc is the FREE UNBOUND plasma concentration (Figure 2A; the
    # paper fits the unbound plasma derived as Cp_total * fu_p). For
    # users who want total plasma, divide Cc by fu_p (also exposed
    # below as Cplasma_total for convenience).
    # ----------------------------------------------------------------
    Cplasma_total <- central / vc
    Cc            <- Cplasma_total * fu_p
    Cmuscle       <- muscle  / v_muscle
    Cadipose      <- adipose / v_adipose
    Cpmncyto      <- pmn     / v_pmn

    # Derived total PML concentration (Equations 8-9) for users who want
    # to compare to measured total PML drug. Equation 9 with volume ratio
    # V_PML(lysosome) / V_PML(cytosol) = 5 / 95 (macrophages / fibroblasts):
    #   C_PML(cytosol,total) = (1/3) * (1 + 5/95) * C_PML,total ~= 0.351 * C_PML,total
    # Inverting and converting unionized cytosol (Cpmncyto) -> total cytosol -> total PML:
    #   C_PML,total = Cpmncyto / (funi_pmn_cyto * 0.351)
    Cpml_total    <- Cpmncyto / (funi_pmn_cyto * 0.351)

    Cc        ~ add(addSd)          + prop(propSd)
    Cmuscle   ~ add(addSd_Cmuscle)  + prop(propSd_Cmuscle)
    Cadipose  ~ add(addSd_Cadipose) + prop(propSd_Cadipose)
    Cpmncyto  ~ add(addSd_Cpmncyto) + prop(propSd_Cpmncyto)
  })
}
