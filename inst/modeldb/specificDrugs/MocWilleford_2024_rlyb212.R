MocWilleford_2024_rlyb212 <- function() {
  description <- "Target-mediated drug disposition (TMDD) model of RLYB212, a recombinant human anti-HPA-1a IgG1 monoclonal antibody, with simultaneous fit of RLYB212 pharmacokinetics and HPA-1a-positive platelet dynamics in HPA-1b/b healthy volunteers (Moc Willeford 2024). Two-compartment SC PK for RLYB212 (central Vc, peripheral Vp, ka, CL) shares its central volume with the HPA-1a-positive platelet distribution (target central), which itself has a peripheral compartment (Vp_target) and first-order natural degradation kdeg. RLYB212 and free receptor form a drug-receptor complex in the central compartment via second-order kon and first-order koff (both fixed at internal / literature values). A novel threshold-gated phagocytic elimination pathway removes coated platelets (both free receptor and complex) once receptor occupancy exceeds a fixed threshold thres = 10%. Parameters and structure per Moc Willeford 2024 Table 2 and the supplement's differential-equation appendix."
  reference <- "Moc Willeford C, Shetty K, Sheridan D, Engler F. Informing pregnancy dose via target-mediated drug disposition modeling and simulations for a recombinant human monoclonal antibody. CPT Pharmacometrics Syst Pharmacol. 2024;13(11):2002-2015. doi:10.1002/psp4.13250"
  vignette <- "MocWilleford_2024_rlyb212"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("target_peripheral")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot             = list(analyte = "RLYB212", units = "mg", specimen = "administration site", verified = FALSE),
    central           = list(analyte = "RLYB212, free receptor, drug-receptor complex", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1       = list(analyte = "RLYB212", units = "mg", specimen = "plasma", verified = FALSE),
    target            = list(analyte = "HPA-1a-positive platelets", units = "mg", specimen = "blood cell", verified = FALSE),
    target_peripheral = list(analyte = "HPA-1a-positive platelets", units = "mg", specimen = "plasma", verified = FALSE),
    complex           = list(analyte = "drug-receptor complex", units = "mg", specimen = "blood cell", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Simple allometric scaling on all drug clearances (CL, Q, QP) with exponent 0.75 and on all volumes (Vc, Vp, Vp_target) with exponent 1, referenced to a 70 kg healthy adult. Applied uniformly in Moc Willeford 2024 Methods (simulation-side scaling equation): (BW_pw/BW_hv)^exp. Time-varying via gestational-age-based scaling in the simulation vignette; the model file itself carries WT as a static per-subject continuous covariate.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_pk_subjects  = 21L,
    n_pd_subjects  = 11L,
    n_studies      = 2L,
    study_names    = c("IPA2001 (Cohorts 1 and 2)", "IPA2003"),
    sex_female_pct = 0,
    sex_note       = "Both IPA2001 and IPA2003 phase 1 study populations were male HPA-1b/b healthy volunteers (Moc Willeford 2024 Table S1).",
    disease_state  = "Healthy adult HPA-1b/b (HPA-1a-negative) volunteers.",
    dose_range     = "IPA2001 Cohort 1: 0.21 mg SC single dose. IPA2001 Cohort 2: 0.29 mg SC on Day 1, then 0.1 mg SC every 2 weeks x 10 weeks. IPA2003: 0.09 mg or 0.29 mg SC single dose on Day 1; on Day 8 an IV transfusion of 10 x 10^9 HPA-1a-positive platelets was given.",
    notes          = "The final PK analysis dataset included 295 measurable observations from 21 subjects across IPA2001 and IPA2003. The final PD analysis dataset included 117 measurable HPA-1a-positive platelet observations from 11 subjects in IPA2003 (Moc Willeford 2024 Tables S3-S4). The model uses molar transformations of both drug and platelet-derived signal: RLYB212 concentrations in the paper are reported in ng/mL; HPA-1a-positive platelet counts are converted to molar receptor equivalents assuming 40,000 HPA-1a receptors per transfused heterozygous HPA-1a/b platelet."
  )

  ini({
    # -----------------------------------------------------------------------------------------------
    # Structural PK / TMDD parameter estimates from Moc Willeford 2024 Table 2 (final TMDD model).
    # Rate constants in 1/h (paper reports them in 1/h). Volumes in L. Clearances in L/h.
    # Concentrations use nM (equivalent to nmol/L) for the receptor / complex sub-model so that
    # kon's paper-reported units 1/(nM*h) are directly consistent.
    # -----------------------------------------------------------------------------------------------
    lka       <- log(0.0058);   label("SC absorption rate constant Ka (1/h)")                                       # Moc Willeford 2024 Table 2: Ka = 0.0058 1/h, RSE 0.05%
    lcl       <- log(0.0119);   label("Apparent RLYB212 clearance CL/F for a 70 kg healthy adult (L/h)")            # Moc Willeford 2024 Table 2: CL/F = 0.0119 L/h, RSE 0.3%
    lvc       <- log(4.61);     label("Apparent RLYB212 central volume V1/F (shared with target central) (L)")      # Moc Willeford 2024 Table 2: V1/F = 4.61 L, RSE 0.4%
    lvp       <- log(5.61);     label("Apparent RLYB212 peripheral volume V2/F (L)")                                # Moc Willeford 2024 Table 2: V2/F = 5.61 L, RSE 0.06%
    lq        <- log(0.693);    label("Apparent RLYB212 inter-compartmental clearance Q/F (L/h)")                   # Moc Willeford 2024 Table 2: Q/F = 0.693 L/h, RSE 3.8%
    lkdeg     <- log(0.0123);   label("First-order HPA-1a-positive platelet natural degradation rate kdeg (1/h)")   # Moc Willeford 2024 Table 2: kdeg = 0.0123 1/h, RSE 0.4%
    lqp       <- log(2.45);     label("Target inter-compartmental clearance QP (L/h)")                              # Moc Willeford 2024 Table 2: QP = 2.45 L/h, RSE 1.7%
    lvp_target <- log(0.523);   label("Target peripheral volume V3 (L)")                                            # Moc Willeford 2024 Table 2: V3 = 0.523 L, RSE 2.6%
    lkphag    <- log(0.197);    label("Threshold-gated phagocytic elimination rate of coated platelets (1/h)")      # Moc Willeford 2024 Table 2: kphagocytosis = 0.197 1/h, RSE 2.3%

    # Fixed binding constants and threshold from Moc Willeford 2024 Table 2 (FIX flag).
    # kon internal-data value 1/(nM*h); koff bivalent-avidity minimal-dissociation value 1/h.
    lkon    <- fixed(log(1.43));     label("Association rate kon (1/(nM*h)); from internal data")             # Moc Willeford 2024 Table 2: kon = 1.43 1/(nM*h) FIXED (internal data on soluble monomeric antigen)
    lkoff   <- fixed(log(0.0001));   label("Dissociation rate koff (1/h); for bivalent minimal dissociation") # Moc Willeford 2024 Table 2: koff = 0.0001 1/h FIXED (avidity-driven bivalent binding value 1.0e-4 /h)
    thres   <- fixed(10);            label("Receptor occupancy threshold for phagocytic elimination (percent)")     # Moc Willeford 2024 Table 2: THRES = 10 percent FIXED (paper Results: estimation ranged 3-22 percent across runs so THRES was fixed at 10)

    # Allometric scaling anchors used in simulation (fixed exponents; not estimated in Moc Willeford 2024).
    e_wt_cl    <- fixed(0.75);  label("Allometric exponent on CL (unitless)")                                       # Moc Willeford 2024 Table 1 / Methods: exponent 0.75 applied to CL FIXED
    e_wt_q     <- fixed(0.75);  label("Allometric exponent on Q (unitless)")                                        # Moc Willeford 2024 Table 1 / Methods: exponent 0.75 applied to intercompartmental clearance Q FIXED
    e_wt_qp    <- fixed(0.75);  label("Allometric exponent on QP (unitless)")                                       # Moc Willeford 2024 Table 1 / Methods: exponent 0.75 applied to target intercompartmental clearance QP FIXED
    e_wt_vc    <- fixed(1);     label("Allometric exponent on V1 (unitless)")                                       # Moc Willeford 2024 Table 1 / Methods: exponent 1 applied to central volume V1 FIXED
    e_wt_vp    <- fixed(1);     label("Allometric exponent on V2 (unitless)")                                       # Moc Willeford 2024 Table 1 / Methods: exponent 1 applied to peripheral volume V2 FIXED
    e_wt_vp_target <- fixed(1); label("Allometric exponent on V3 (unitless)")                                       # Moc Willeford 2024 Table 1 / Methods: exponent 1 applied to target peripheral volume V3 FIXED

    # Log-normal IIV, diagonal per Moc Willeford 2024 Table 2 (between-subject variability estimates).
    # Reported "Estimate" values are the OMEGA variances on the log scale (standard NONMEM omega**2);
    # the paper's "CV (%)" column is the RSE of the OMEGA estimate rather than a CV derived from it
    # (kphagocytosis: omega**2 = 0.401 corresponds to CV = 100*sqrt(exp(0.401)-1) = 68.6%, not the
    # tabulated 41.9%; the 41.9% is the RSE of the OMEGA estimate). Shrinkage values (6.5-37.9%)
    # are recorded in the vignette. Between-subject variability was estimated only for Ka, V1, CL,
    # and kphagocytosis.
    etalka   ~ 0.634                                                                                                # Moc Willeford 2024 Table 2: OMEGA Ka = 0.634 (RSE 25.7%, shrinkage 6.5%)
    etalvc   ~ 0.137                                                                                                # Moc Willeford 2024 Table 2: OMEGA V1 = 0.137 (RSE 66.8%, shrinkage 16.6%)
    etalcl   ~ 0.0388                                                                                               # Moc Willeford 2024 Table 2: OMEGA CL = 0.0388 (RSE 53.6%, shrinkage 9.7%)
    etalkphag ~ 0.401                                                                                               # Moc Willeford 2024 Table 2: OMEGA kphagocytosis = 0.401 (RSE 41.9%, shrinkage 37.9%)

    # Independent proportional residual error for the two outputs (Moc Willeford 2024 Table 2).
    # The paper labels these as "Proportional residual error of drug" and "Proportional residual
    # error of platelet"; both are on the linear scale of the respective observation.
    propSd           <- 0.152;   label("Proportional residual error on RLYB212 concentration (fraction)")               # Moc Willeford 2024 Table 2: prop residual error of drug = 0.152, RSE 1.2%
    propSd_platelet  <- 0.148;   label("Proportional residual error on HPA-1a-positive platelet concentration (fraction)") # Moc Willeford 2024 Table 2: prop residual error of platelet = 0.148, RSE 0.8%
  })

  model({
    # -------------------------------------------------------------------------------------
    # Unit convention.
    # Drug amounts (depot, central, peripheral1) are in mg so that the event-table `amt`
    # column can be given directly in mg. Target-species amounts (target, target_peripheral,
    # complex) are in nmol. Concentrations for the binding kinetics are in nM = nmol/L,
    # matching the paper's fixed kon (1/(nM*h)) and koff (1/h). The molecular weight of
    # RLYB212 (a human IgG1) is taken as 150 kDa per the paper's ~10 ng/mL <=> ~0.067 nM
    # conversion; this is the mwdrug constant used to bridge drug mass (mg) and molar
    # concentration (nM) at the binding interface.
    # -------------------------------------------------------------------------------------
    mwdrug <- 0.150   # mg/nmol; RLYB212 IgG1 ~150 kDa (Moc Willeford 2024: 10 ng/mL <=> 0.067 nM implies MW ~= 150 kDa)

    # Individual PK / TMDD parameters (allometric scaling to a 70 kg reference; WT is per-subject).
    ka        <- exp(lka + etalka)
    cl        <- exp(lcl + etalcl)      * (WT / 70)^e_wt_cl
    vc        <- exp(lvc + etalvc)      * (WT / 70)^e_wt_vc
    vp        <- exp(lvp)               * (WT / 70)^e_wt_vp
    q         <- exp(lq)                * (WT / 70)^e_wt_q
    qp        <- exp(lqp)               * (WT / 70)^e_wt_qp
    vp_target <- exp(lvp_target)        * (WT / 70)^e_wt_vp_target
    kdeg      <- exp(lkdeg)
    kphag     <- exp(lkphag + etalkphag)
    kon       <- exp(lkon)
    koff      <- exp(lkoff)

    # Micro-rate constants for the drug's linear distribution / elimination.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    # Target compartment rate constants (target central shares vc with drug; peripheral is vp_target).
    ktc12 <- qp / vc
    ktc21 <- qp / vp_target

    # Receptor occupancy (percent of platelet HPA-1a sites bound by RLYB212 in the central compartment).
    # complex and target are amounts in nmol; the numerator/denominator ratio is unit-free.
    ro  <- 100 * complex / (target + complex + 1e-30)
    swi <- (ro > thres)

    # ------- Drug ODEs (mg / h) -------
    # Binding rate at the drug side (mg / h):
    #   rate_bind_drug = kon [1/(nM*h)] * (central/(vc*mwdrug)) [nM] * (target/vc) [nM] * vc [L] * mwdrug [mg/nmol]
    #                  = kon * central * target / vc.
    # Dissociation gives drug back at:
    #   rate_dissoc_drug = koff [1/h] * complex [nmol] * mwdrug [mg/nmol] = koff * complex * mwdrug.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 -
                          (kon * central * target / vc) + koff * complex * mwdrug
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------- Target (HPA-1a-positive free receptor) ODEs (nmol / h) -------
    # Binding rate at the target side (nmol / h):
    #   rate_bind_target = kon [1/(nM*h)] * (central/(vc*mwdrug)) [nM] * (target/vc) [nM] * vc [L]
    #                    = kon * central * target / (vc * mwdrug).
    d/dt(target)          <- -kdeg * target - kphag * target * swi
                             - ktc12 * target + ktc21 * target_peripheral
                             - (kon * central * target / (vc * mwdrug)) + koff * complex
    d/dt(target_peripheral) <-  ktc12 * target - ktc21 * target_peripheral

    # ------- Drug-receptor complex ODE (nmol / h) -------
    d/dt(complex) <-  (kon * central * target / (vc * mwdrug)) - koff * complex - kphag * complex * swi

    # -------------------------------------------------------------------------------------
    # Observation and error model.
    # Drug: Cc in ng/mL from central [mg] and vc [L]:  Cc = 1000 * central / vc (mg/L to ng/mL).
    # Platelet: platelet in nM as free-receptor concentration = target / vc (the supplement
    # labels C4 as "concentration of free receptor in plasma compartment"; the observed
    # platelet-derived signal after the paper transformation Pi = PS * Nt * 40000 / NA
    # is on the same molar scale). Named `platelet` (rather than `Cplt`) so the checker does
    # not treat the observation as a PK metabolite of the drug.
    # -------------------------------------------------------------------------------------
    Cc       <- 1000 * central / vc
    platelet <- target / vc

    Cc       ~ prop(propSd)
    platelet ~ prop(propSd_platelet)
  })
}
