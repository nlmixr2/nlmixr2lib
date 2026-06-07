Brain_2008_ifosfamide <- function() {
  description <- "Joint population PK / PD model for ifosfamide in adults with advanced solid tumours (Brain 2008, n=17, single-agent ifosfamide 9 g/m^2 per cycle by either 3 h x 3 daily or 72 h continuous infusion, n=1 randomised crossover, NONMEM VI FOCE INTERACTION). One-compartment ifosfamide PK with Kerbusch 2000-style autoinduction of clearance via a relative enzyme-pool state (drug inhibits enzyme degradation), three coupled apparent-volume metabolite states (4-hydroxy-ifosfamide, 3-dechloroethyl-ifosfamide, 2-dechloroethyl-ifosfamide), an indirect-response model for urinary beta-2-microglobulin (BMG, renal tubular toxicity) with linear stimulation of production by parent ifosfamide concentration, and a five-compartment Friberg-style myelosuppression chain for absolute neutrophil count (ANC) with linear inhibition of proliferation by parent ifosfamide concentration and (CIRC0 / circ)^gamma feedback. No covariates were retained in the final model (one outlier patient on carbamazepine was excluded prior to the final analysis)."
  reference <- paste(
    "Brain EGC, Rezai K, Lokiec F, Gutierrez M, Urien S (2008).",
    "Population pharmacokinetics and exploratory pharmacodynamics of",
    "ifosfamide according to continuous or short infusion schedules:",
    "an n = 1 randomized study.",
    "Br J Clin Pharmacol 65(4):607-610.",
    "doi:10.1111/j.1365-2125.2007.03095.x.",
    "Autoinduction structure carried from Kerbusch T et al. (2000)",
    "Br J Clin Pharmacol 49(6):555-561",
    "doi:10.1046/j.1365-2125.2000.00217.x;",
    "myelosuppression structure from Friberg LE et al. (2002)",
    "J Clin Oncol 20(24):4713-4721 doi:10.1200/JCO.2002.02.140.",
    sep = " "
  )
  vignette <- "Brain_2008_ifosfamide"
  paper_specific_compartments <- c("ohif", "decloro3", "decloro2", "bmg")
  paper_specific_etas <- c(
    "etalfovv_ohif", "etalkm_ohif",
    "etalfovv_decloro3", "etalfovv_decloro2",
    "etalslope_bmg", "etalrbase_bmg",
    "etalmtt_anc", "etalrbase_anc"
  )
  paper_specific_residual_sds <- c(
    "expSd_ohif", "expSd_decloro3", "expSd_decloro2",
    "propSd_bmg", "addSd_bmg",
    "propSd_anc", "addSd_anc"
  )
  units <- list(
    time          = "h",
    dosing        = "umol",
    concentration = "umol/L for Cc (ifosfamide), Cohif, Cdecloro3, Cdecloro2 (metabolites); mg/L for BMG; cells/mm^3 for ANC"
  )

  # No covariates were retained in the final Brain 2008 model. The one
  # carbamazepine-coadministration observation (eight-fold CLINIT increase, OFV
  # drop of 106 units) was excluded prior to the final fit and is documented in
  # vignette Assumptions and deviations rather than as a model covariate.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "35-68 years",
    height_range   = "1.39-1.89 m",
    sex_female_pct = 50.0,
    disease_state  = "Adults with advanced solid tumours receiving single-agent ifosfamide",
    dose_range     = "Ifosfamide 3 g/m^2 per day for 3 days (total 9 g/m^2 per cycle) administered as either three successive 3 h daily IV infusions or a single 72 h continuous IV infusion; total per-cycle dose 12.6-16.8 g. Each subject received both schedules in an n = 1 randomised crossover separated by 3 weeks.",
    regions        = "France (single centre, Rene Huguenin Cancer Centre, Saint-Cloud)",
    notes          = paste(
      "Original enrollment was 17 patients (9 male, 8 female); 12 had two",
      "complete pharmacokinetic evaluations. One patient receiving",
      "carbamazepine co-administration (eight-fold increase of CLINIT, OFV drop",
      "of 106 units) was excluded from the final analysis, leaving 16 patients",
      "for the reported parameter estimates. Pharmacokinetic dataset contained",
      "572, 513 and 572 concentrations for ifosfamide, 4-hydroxy-ifosfamide",
      "and the (combined) dechloroethylated metabolites, respectively."
    )
  )

  ini({
    # =========================================================================
    # Ifosfamide one-compartment PK with Kerbusch 2000-style autoinduction
    # (Brain 2008 page 608 column 2 "The auto-induction model described
    # ifosfamide pharmacokinetics, whatever the infusion duration. ... The
    # final parameter estimates were V 46 l (6%), CLINIT 3.44 l h-1 (4%),
    # EC50 22 mmol l-1 (4%) and MTT 62 h (6%)").
    # =========================================================================
    lvc      <- log(46);    label("Ifosfamide central volume of distribution V (L)")                                                # Brain 2008 Results, p.608
    lcl      <- log(3.44);  label("Initial (baseline, enzyme = 1) ifosfamide clearance CLINIT (L/h)")                              # Brain 2008 Results, p.608
    lec50    <- log(22);    label("Ifosfamide concentration at half-maximal autoinduction effect EC50 (umol/L)")                   # Brain 2008 Results, p.608
    lmtt_ind <- log(62);    label("Mean transit time of the autoinduction enzyme pool MTT_ind (h); ktr_ind = 2 / MTT_ind")          # Brain 2008 Results, p.608

    # =========================================================================
    # Metabolite PK parameters (Brain 2008 page 609 column 1 and Equation set
    # for metabolites on p.608: V_m * dC_m/dt = f_m * CL_IF * C_IF - CL_m * C_m,
    # only f_m / V_m and K_m = CL_m / V_m identifiable). Order in the source
    # results paragraph is 4-hydroxy-ifosfamide, 3-dechloroethyl-ifosfamide,
    # 2-dechloroethyl-ifosfamide.
    # =========================================================================
    lfovv_ohif     <- log(0.0019);  label("Apparent f_m / V_m for 4-hydroxy-ifosfamide (1/L)")                                      # Brain 2008 Results, p.609
    lkm_ohif       <- log(0.14);    label("Metabolite elimination rate K_m for 4-hydroxy-ifosfamide (1/h)")                         # Brain 2008 Results, p.609
    lfovv_decloro3 <- log(0.0063);  label("Apparent f_m / V_m for 3-dechloroethyl-ifosfamide (1/L)")                                # Brain 2008 Results, p.609
    lkm_decloro3   <- log(0.020);   label("Metabolite elimination rate K_m for 3-dechloroethyl-ifosfamide (1/h)")                   # Brain 2008 Results, p.609
    lfovv_decloro2 <- log(0.0043);  label("Apparent f_m / V_m for 2-dechloroethyl-ifosfamide (1/L)")                                # Brain 2008 Results, p.609
    lkm_decloro2   <- log(0.036);   label("Metabolite elimination rate K_m for 2-dechloroethyl-ifosfamide (1/h)")                   # Brain 2008 Results, p.609

    # =========================================================================
    # Urinary beta-2-microglobulin (BMG) indirect-response PD model
    # (Brain 2008 page 608 column 2 equation "dR/dt = K_TR (1 + CIF*SLOPEIF)
    # - K_TR * R" with baseline scaling, and page 609 column 2 parameter
    # estimates "SLOPEIF 0.32 l mg-1 (34%), MTT 243 h (27%) and baseline
    # 0.05 mg l-1 (28%)"). Drug effect uses ifosfamide concentration converted
    # from umol/L to mg/L via the ifosfamide molecular weight; see the model()
    # block. Brain's MTT convention is MTT = (n_state + 1) / K_TR, so for the
    # single-state BMG model ktr_bmg = 2 / MTT_bmg.
    # =========================================================================
    lslope_bmg <- log(0.32);  label("Linear drug-effect slope SLOPEIF on BMG production, applied to Cif (mg/L) (L/mg)")            # Brain 2008 Results, p.609
    lmtt_bmg   <- log(243);   label("BMG indirect-response mean transit time MTT_bmg (h); ktr_bmg = 2 / MTT_bmg")                   # Brain 2008 Results, p.609
    lrbase_bmg <- log(0.05);  label("Baseline urinary beta-2-microglobulin concentration (mg/L)")                                   # Brain 2008 Results, p.609

    # =========================================================================
    # Friberg-style absolute neutrophil count (ANC) myelosuppression model
    # (Brain 2008 page 608 column 2 / page 609 column 2). Five-compartment
    # chain: PROL (= precursor1) + 3 transit (precursor2..precursor4) + CIRC
    # (circ; observation). Linear inhibition of proliferation by parent
    # ifosfamide concentration. Brain's MTT convention is MTT = (n_state + 1) /
    # K_TR, so for the five-state ANC chain ktr_anc = 6 / MTT_anc; SLOPEIF
    # multiplies Cif in umol/L directly.
    # =========================================================================
    lslope_anc <- log(0.014);   label("Linear drug-effect slope SLOPEIF on ANC proliferation rate (L/umol)")                        # Brain 2008 Results, p.609
    lmtt_anc   <- log(150);     label("ANC chain mean transit time MTT_anc (h); ktr_anc = 6 / MTT_anc")                              # Brain 2008 Results, p.609
    lrbase_anc <- log(4490);    label("Baseline circulating ANC count CIRC0 (cells/mm^3)")                                          # Brain 2008 Results, p.609
    lgamma     <- log(0.16);    label("Friberg feedback exponent gamma on (CIRC0 / circ) (unitless)")                                # Brain 2008 Results, p.609

    # =========================================================================
    # IIVs. Brain 2008 page 608 column 2 "Variabilities were expressed as the
    # square root of the variances, w^2 or s^2." -> reported values are SDs on
    # the natural-log scale; variance for the nlmixr2 omega is (reported SD)^2.
    # The "(NN%)" tail in the paper is the RSE of the SD estimate, not the
    # variability magnitude itself. Parameters without a reported ISV are not
    # declared as etas (per paper: "When ISV was not given for a parameter, it
    # meant that it was not statistically significant and its deletion did
    # not alter the fit and OFV").
    # =========================================================================
    etalvc            ~ 0.14^2                                                # ISV V = 0.14 SD on log scale       (Brain 2008 p.608)
    etalcl            ~ 0.18^2                                                # ISV CLINIT = 0.18 SD on log scale  (Brain 2008 p.608)
    etalfovv_ohif     ~ 0.42^2                                                # ISV fm/Vm 4-OH = 0.42 SD           (Brain 2008 p.609)
    etalkm_ohif       ~ 0.30^2                                                # ISV Km 4-OH = 0.30 SD              (Brain 2008 p.609)
    etalfovv_decloro3 ~ 0.21^2                                                # ISV fm/Vm 3-DCE = 0.21 SD          (Brain 2008 p.609)
    etalfovv_decloro2 ~ 0.26^2                                                # ISV fm/Vm 2-DCE = 0.26 SD          (Brain 2008 p.609)
    etalslope_bmg     ~ 3.4^2                                                 # ISV SLOPEIF (BMG) = 3.4 SD         (Brain 2008 p.609); large IIV consistent with the noisy urinary biomarker
    etalrbase_bmg     ~ 1.1^2                                                 # ISV baseline (BMG) = 1.1 SD        (Brain 2008 p.609)
    etalmtt_anc       ~ 0.30^2                                                # ISV MTT (ANC) = 0.30 SD            (Brain 2008 p.609)
    etalrbase_anc     ~ 0.46^2                                                # ISV CIRC0 = 0.46 SD                (Brain 2008 p.609)

    # =========================================================================
    # Residual error. Brain 2008 page 608 column 2 "ISV and residual
    # variabilities were modelled as exponential errors" -> exponential
    # (log-normal) residual on parent and metabolite Cc. BMG and ANC
    # residuals were reported as proportional plus additive components and
    # are encoded as add + prop on the linear scale.
    # =========================================================================
    expSd          <- 0.22;   label("Log-normal residual SD on ifosfamide concentration Cc (log scale)")                            # Brain 2008 p.608
    expSd_ohif     <- 0.71;   label("Log-normal residual SD on 4-hydroxy-ifosfamide concentration (log scale)")                     # Brain 2008 p.609
    expSd_decloro3 <- 0.36;   label("Log-normal residual SD on 3-dechloroethyl-ifosfamide concentration (log scale)")               # Brain 2008 p.609
    expSd_decloro2 <- 0.36;   label("Log-normal residual SD on 2-dechloroethyl-ifosfamide concentration (log scale)")               # Brain 2008 p.609
    propSd_bmg     <- 0.73;   label("Proportional residual SD on BMG (fraction)")                                                   # Brain 2008 p.609
    addSd_bmg      <- 0.02;   label("Additive residual SD on BMG (mg/L)")                                                           # Brain 2008 p.609
    propSd_anc     <- 0.34;   label("Proportional residual SD on ANC (fraction)")                                                   # Brain 2008 p.609
    addSd_anc      <- 650;    label("Additive residual SD on ANC (cells/mm^3)")                                                     # Brain 2008 p.609
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual parameters.
    # -----------------------------------------------------------------------
    vc      <- exp(lvc      + etalvc)
    cl      <- exp(lcl      + etalcl)                # CLINIT (baseline ifosfamide clearance)
    ec50    <- exp(lec50)
    mtt_ind <- exp(lmtt_ind)
    ktr_ind <- 2 / mtt_ind                            # Brain MTT convention: MTT = (n + 1) / K_TR; n = 1 here

    fovv_ohif     <- exp(lfovv_ohif     + etalfovv_ohif)
    km_ohif       <- exp(lkm_ohif       + etalkm_ohif)
    fovv_decloro3 <- exp(lfovv_decloro3 + etalfovv_decloro3)
    km_decloro3   <- exp(lkm_decloro3)
    fovv_decloro2 <- exp(lfovv_decloro2 + etalfovv_decloro2)
    km_decloro2   <- exp(lkm_decloro2)

    slope_bmg <- exp(lslope_bmg + etalslope_bmg)
    mtt_bmg   <- exp(lmtt_bmg)
    rbase_bmg <- exp(lrbase_bmg + etalrbase_bmg)
    ktr_bmg   <- 2 / mtt_bmg                          # n = 1 single-state IDR; Brain convention

    slope_anc <- exp(lslope_anc)
    mtt_anc   <- exp(lmtt_anc + etalmtt_anc)
    circ0     <- exp(lrbase_anc + etalrbase_anc)
    gamma     <- exp(lgamma)
    ktr_anc   <- 6 / mtt_anc                          # n = 5 Friberg chain; Brain "MTT of the system is 6/K_TR"

    # -----------------------------------------------------------------------
    # 2. Ifosfamide PK with Kerbusch 2000 autoinduction.
    # The relative enzyme amount A2 (compartment `enzyme`) follows the
    # Kerbusch 2000 form dA2/dt = K_enz - K_enz * A2 * (1 - C_p / (C_p + IC50))
    # with A2(0) = 1, which is what Brain 2008 carries forward as their
    # auto-induction model (Brain p.608 equation set, Brain "indirect response
    # model previously described [4]" -> Kerbusch 2000 (BJCP) and the
    # Kerbusch 2001 (Drug Metab Dispos) replicate). Effective clearance is
    # CL_app = CLINIT * enzyme.
    # -----------------------------------------------------------------------
    Cc <- central / vc                                # ifosfamide concentration (umol/L)

    d/dt(enzyme)  <- ktr_ind - ktr_ind * enzyme * (1 - Cc / (Cc + ec50))
    cl_app        <- cl * enzyme                      # apparent (induced) ifosfamide clearance
    d/dt(central) <- -cl_app * Cc

    enzyme(0)  <- 1                                   # baseline relative enzyme amount

    # -----------------------------------------------------------------------
    # 3. Metabolite concentrations (each state holds the metabolite plasma
    # concentration in umol/L; see paper equation
    #   V_m * dC_m/dt = f_m * CL_IF * C_IF - CL_m * C_m
    # divided through by V_m yields
    #   dC_m/dt = (f_m / V_m) * CL_IF * C_IF - K_m * C_m
    # with only f_m / V_m and K_m = CL_m / V_m identifiable per Brain 2008).
    # The instantaneous ifosfamide clearance cl_app is used (the auto-induced
    # value), consistent with Brain's sequential POSTHOC-feeding approach
    # where Bayesian individual PK estimates drove the metabolite dataset.
    # -----------------------------------------------------------------------
    d/dt(ohif)     <- fovv_ohif     * cl_app * Cc - km_ohif     * ohif
    d/dt(decloro3) <- fovv_decloro3 * cl_app * Cc - km_decloro3 * decloro3
    d/dt(decloro2) <- fovv_decloro2 * cl_app * Cc - km_decloro2 * decloro2

    Cohif     <- ohif                                 # 4-hydroxy-ifosfamide (umol/L)
    Cdecloro3 <- decloro3                             # 3-dechloroethyl-ifosfamide (umol/L)
    Cdecloro2 <- decloro2                             # 2-dechloroethyl-ifosfamide (umol/L)

    # -----------------------------------------------------------------------
    # 4. BMG indirect-response PD. Linear stimulation of production by the
    # parent ifosfamide concentration converted to mg/L via the ifosfamide
    # molecular weight (Brain reports SLOPEIF in L/mg; the underlying PK is
    # carried in umol/L per the paper's molar-concentration convention).
    # -----------------------------------------------------------------------
    mw_if      <- 261.1                               # ifosfamide molecular weight (g/mol)
    cif_mg     <- Cc * mw_if / 1000                   # Cc (umol/L) -> mg/L for the BMG drug-effect term
    edrug_bmg  <- slope_bmg * cif_mg
    d/dt(bmg)  <- ktr_bmg * rbase_bmg * (1 + edrug_bmg) - ktr_bmg * bmg
    bmg(0)     <- rbase_bmg
    BMG        <- bmg                                 # urinary beta-2-microglobulin (mg/L)

    # -----------------------------------------------------------------------
    # 5. Friberg myelosuppression PD on ANC. Five-compartment chain:
    # precursor1 = PROL (proliferation), precursor2..precursor4 = R2/R3/R4
    # (transit), circ = R5 = CIRC (circulating). Linear inhibition of
    # proliferation by parent ifosfamide concentration (in umol/L; SLOPEIF
    # carries units L/umol) and (CIRC0 / circ)^gamma feedback.
    # -----------------------------------------------------------------------
    edrug_anc <- 1 - slope_anc * Cc
    feed      <- (circ0 / circ)^gamma

    d/dt(precursor1) <- ktr_anc * precursor1 * edrug_anc * feed - ktr_anc * precursor1
    d/dt(precursor2) <- ktr_anc * precursor1 - ktr_anc * precursor2
    d/dt(precursor3) <- ktr_anc * precursor2 - ktr_anc * precursor3
    d/dt(precursor4) <- ktr_anc * precursor3 - ktr_anc * precursor4
    d/dt(circ)       <- ktr_anc * precursor4 - ktr_anc * circ

    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0
    circ(0)       <- circ0

    ANC <- circ                                       # absolute neutrophil count (cells/mm^3)

    # -----------------------------------------------------------------------
    # 6. Observation models. Parent and metabolites: log-normal residual
    # ("exponential errors" per Brain 2008). BMG and ANC: combined
    # proportional + additive on the linear scale.
    # -----------------------------------------------------------------------
    Cc        ~ lnorm(expSd)
    Cohif     ~ lnorm(expSd_ohif)
    Cdecloro3 ~ lnorm(expSd_decloro3)
    Cdecloro2 ~ lnorm(expSd_decloro2)
    BMG       ~ add(addSd_bmg) + prop(propSd_bmg)
    ANC       ~ add(addSd_anc) + prop(propSd_anc)
  })
}
