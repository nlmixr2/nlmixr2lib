Groenendaal_2007_morphine_brain_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Non-linear blood-brain barrier (BBB)",
    "distribution model for morphine published by Groenendaal et al. (2007,",
    "Br J Pharmacol 151(4):701-712). A three-compartment blood disposition",
    "model (paper Table 2, NONMEM ADVAN11 TRANS4; linear body-weight effects",
    "on CL and on the first peripheral volume V2) drives a single brain",
    "extracellular-fluid (ECF) compartment sampled by intracerebral",
    "striatal microdialysis. Mass exchange across the BBB is split into",
    "three explicit terms (paper equations 9-11): bidirectional passive",
    "diffusion kdiff * (Cblood - Cecf), a saturable active influx",
    "N*max * Cblood / (C50 + Cblood) that is already near-saturated at the",
    "lowest studied dose, and an active P-glycoprotein-mediated efflux",
    "keff * Cecf. Co-infusion of the Pgp inhibitor GF120918 (elacridar)",
    "lowers keff from 0.0195 to 0.0113 1/min (a 42% reduction) and leaves",
    "kdiff, N*max and C50 unchanged; the blood disposition is unaffected by",
    "GF120918. Inter-animal variability is on CL, V2, kdiff and keff, with a",
    "kdiff-keff covariance block (paper Table 3). Because the paper reports",
    "the volume-aggregated brain parameters (kdiff = Qdiff/Vecf,",
    "N*max = Nmax/Vecf, keff = Qeff/Vecf) and never identifies Vecf itself,",
    "the brain_ecf state carries a concentration (ng/mL) rather than an",
    "amount; see the validation vignette for the dimensional argument.",
    sep = " "
  )
  reference <- paste(
    "Groenendaal D, Freijer J, de Mik D, Bouw MR, Danhof M, de Lange ECM.",
    "Population pharmacokinetic modelling of non-linear brain distribution",
    "of morphine: influence of active saturable influx and P-glycoprotein",
    "mediated efflux.",
    "Br J Pharmacol. 2007;151(4):701-712.",
    "doi:10.1038/sj.bjp.0707257.",
    sep = " "
  )
  vignette <- "Groenendaal_2007_morphine_brain_rat"

  # `lkdiff_bbb` extends the `k<process>_bbb` BBB rate-constant family founded
  # by Pardridge_2023_propranolol_pbpk.R / Pardridge_2023_imipramine_pbpk.R
  # (`lkinflux_bbb`, `lkefflux_bbb`); `lkefflux_bbb` itself is reused
  # unchanged. `lNstarMax` follows Geldof_2008_fluvoxamine_rat.R, which
  # encodes the same Leiden group's volume-aggregated N-star saturable
  # transport capacity under that name. Declared here so
  # checkModelConventions() records them as deliberate.
  paper_specific_etas <- c("etalkdiff_bbb", "etalkefflux_bbb")
  paper_specific_residual_sds <- c("propSd_Cbrain_ecf")

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "morphine", units = "ng", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "ng", specimen = "whole blood", verified = TRUE),
    peripheral2 = list(analyte = "morphine", units = "ng", specimen = "whole blood", verified = TRUE),
    # NOTE the units on brain_ecf: this state is a CONCENTRATION, not an
    # amount. Paper equations 10-11 divide the ECF mass balance through by
    # the (never-identified) brain ECF distribution volume Vecf, so every
    # published brain parameter is already volume-aggregated and the state
    # variable that they act on is Cecf in ng/mL. Encoding an amount here
    # would require inventing Vecf, which the paper does not report.
    brain_ecf   = list(analyte = "morphine", units = "ng/mL (concentration state; see note)", specimen = "brain ISF", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at the start of the experiment.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per animal. Enters CL and V2 (peripheral 1) through the",
        "paper's centred LINEAR covariate form (equation 7),",
        "P_i = theta1 * (1 + theta2 * (BW_i - median BW)), NOT a power /",
        "allometric form. The paper does not print the numeric value of",
        "`median BW`; 0.300 kg is used here (see the vignette Errata). Rats",
        "weighed 250-350 g, the per-group means in Table 1 span 0.260-0.306",
        "kg, and the N = 20 largest group has a mean of exactly 0.300 kg.",
        "Keep simulated weights inside 0.25-0.35 kg: the linear form drives",
        "V2 negative below BW = 0.182 kg and CL negative below 0.113 kg."
      ),
      source_name        = "BW"
    ),
    CONMED_ELACRIDAR = list(
      description        = "Indicator for co-infusion of the P-glycoprotein inhibitor GF120918 (elacridar): 1 = elacridar arm, 0 = vehicle arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (vehicle co-infusion, no elacridar)",
      notes              = paste(
        "Groenendaal 2007 regimen: a 1 min bolus infusion of 6 mg/kg",
        "GF120918 in dimethyl sulphoxide followed by a continuous infusion",
        "of 25 ng/min in DMSO / glucose / cyclodextrine 5/5/10% in saline,",
        "started at t = -120 min and maintained for the whole experiment",
        "(paper Experimental procedures). The measured mean steady-state",
        "GF120918 blood concentration was 214 +/- 16 ng/mL, described by the",
        "sponsor as sufficient to block Pgp in vivo. Treated as time-fixed",
        "per animal here because the inhibitor infusion begins 2 h before",
        "the morphine dose and runs past the last morphine observation. Only",
        "the 4 mg/kg dose level was studied with elacridar; the paper found",
        "NO effect of GF120918 on the blood disposition, so the covariate",
        "acts on keff alone."
      ),
      source_name        = "(paper treatment-arm indicator: 'GF120918' vs 'Vehicle' in Table 1; paper equation 8 defines the factor as 1 when GF120918 is co-infused and 0 when vehicle is co-infused)"
    )
  )

  population <- list(
    species        = "rat (male Wistar, Charles River, Maastricht, The Netherlands)",
    n_subjects     = 71L,
    n_studies      = 1L,
    age_range      = "adult (age not reported; animals housed at least 7 days after arrival, then 10 days of recovery after electrode / cannula implantation)",
    weight_range   = "250-350 g body weight; per-group means 0.260-0.306 kg (paper Table 1)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy rats prepared with four indwelling cannulas (right femoral",
      "artery for serial arterial blood sampling; both left jugular veins",
      "for morphine and midazolam; right femoral vein for GF120918 /",
      "vehicle and vecuronium bromide). The microdialysis animals also",
      "carried a CMA/12 guide cannula with a 4 mm probe in the right",
      "striatum (AP +0.5, L +2.7 mm from bregma, V -3.5 mm from skull); the",
      "EEG animals carried four cortical screw electrodes. All rats received",
      "a continuous midazolam infusion (5.5 mg/kg/h after a Wagner loading",
      "scheme; measured steady state 937 +/- 56 ng/mL) to suppress",
      "opioid-induced seizure activity. Animals in the 40 mg/kg arm",
      "developed severe respiratory depression and muscle rigidity and were",
      "artificially ventilated from about 5 min after the start of the",
      "morphine infusion, with intravenous vecuronium bromide (0.15 mg,",
      "then 0.10 mg as needed) for muscle relaxation."
    ),
    dose_range     = paste(
      "Single 10 min zero-order intravenous infusion of morphine",
      "hydrochloride in saline at 4, 10 or 40 mg/kg. The 4 mg/kg level was",
      "run both with vehicle and with the Pgp inhibitor GF120918",
      "(elacridar); the 10 and 40 mg/kg levels were vehicle only."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "Two model layers, fitted sequentially in NONMEM V level 1.1 with",
      "FOCE INTERACTION. (1) The blood three-compartment model (Table 2) was",
      "fitted to the arterial blood morphine concentrations of all 71 rats",
      "in Table 1 (15 samples per animal over 360 min post-dose), pooling",
      "the EEG-only, EEG/MD and MD groups across 4, 10 and 40 mg/kg and",
      "across vehicle and GF120918; co-infusion of GF120918 had no effect on",
      "the blood PK. (2) The brain ECF transport model (Table 3, ADVAN6)",
      "was fitted to the microdialysate profiles with the blood model as the",
      "input function. Only the blank-perfusate (0 ng/mL) microdialysis rats",
      "contributed brain ECF observations: the 50 and 500 ng/mL",
      "retrodialysis / DNNF animals yielded no usable post-dose ECF data and",
      "were omitted (paper Microdialysis probe recovery). From Table 1 that",
      "leaves 10 (4 mg/kg vehicle EEG/MD) + 3 (4 mg/kg vehicle MD) + 10",
      "(4 mg/kg GF120918 EEG/MD) + 9 (40 mg/kg vehicle EEG/MD) = 32 rats",
      "with 25-30 dialysate fractions each over 4-6 h. Microdialysate",
      "concentrations were converted to true ECF concentrations with a fixed",
      "group in-vivo recovery of 16.1% for the two 4 mg/kg arms and 20.3%",
      "for the 40 mg/kg arm (individual recovery values were not estimable",
      "from the design). Assay LLOQ was 25 ng/mL in blood (50 uL sample) and",
      "0.5 ng/mL in dialysate (40 uL sample). The dose-normalised brain ECF",
      "AUC(0-165 min) fell from 6810 +/- 1890 (4 mg/kg) and 8460 +/- 2790",
      "(4 mg/kg + GF120918) to 3990 +/- 2180 ng.h/mL (40 mg/kg), the",
      "non-linearity this model was built to explain. The companion",
      "publication (Groenendaal et al. 2007, Br J Pharmacol 151(4):713-720,",
      "doi:10.1038/sj.bjp.0707258) reports the EEG PK-PD biophase model from",
      "the same experiments and is a separate model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Blood disposition -- Groenendaal 2007 Table 2 (p707), three-
    # compartment model, NONMEM PREDPP ADVAN11 TRANS4. Volumes in mL,
    # clearances in mL/min, absolute (whole-rat) rather than per kg.
    #
    # Body weight enters CL and V2 through paper equation 7:
    #     P_i = theta1 * (1 + theta2 * (BW_i - median BW))
    # a CENTRED LINEAR (not power / allometric) form. `theta1` is the
    # table's "Intercept" row, `theta2` its "Slope factor" row (units
    # 1/kg). Table 2 transposes the two Slope factor confidence
    # intervals -- each point estimate falls outside its printed CI and
    # inside the other's, and estimate +/- 1.96 * CV% * estimate
    # reproduces them swapped (5.35 -> 2.71-7.99, printed against V2;
    # 8.50 -> 5.65-11.35, printed against Cl). Point estimates are
    # unaffected. See vignette Errata.
    # ------------------------------------------------------------------
    lcl     <- log(20.0);  label("Blood clearance CL at the reference body weight (mL/min)")                        # Table 2, Cl Intercept = 20.0 mL/min (CV 5.6%)
    e_wt_cl <- 5.35;       label("Linear body-weight slope factor on CL (1/kg), applied as CL * (1 + e_wt_cl * (WT - 0.300))") # Table 2, Cl Slope factor = 5.35 (CV 25.2%); paper equation 7
    lvc     <- log(68.1);  label("Central blood volume of distribution V1 (mL)")                                    # Table 2: V1 = 68.1 mL (CV 16.7%)
    lq      <- log(15.5);  label("Inter-compartmental clearance central <-> peripheral1, Q2 (mL/min)")              # Table 2: Q2 = 15.5 mL/min (CV 11.3%)
    lvp     <- log(739);   label("Peripheral 1 volume of distribution V2 at the reference body weight (mL)")        # Table 2, V2 Intercept = 739 mL (CV 7.6%)
    e_wt_vp <- 8.50;       label("Linear body-weight slope factor on V2 (1/kg), applied as V2 * (1 + e_wt_vp * (WT - 0.300))") # Table 2, V2 Slope factor = 8.50 (CV 17.1%); paper equation 7
    lq2     <- log(17.8);  label("Inter-compartmental clearance central <-> peripheral2, Q3 (mL/min)")              # Table 2: Q3 = 17.8 mL/min (CV 18.4%)
    lvp2    <- log(133);   label("Peripheral 2 volume of distribution V3 (mL)")                                     # Table 2: V3 = 133 mL (CV 15.9%)

    # ------------------------------------------------------------------
    # Brain ECF transport -- Groenendaal 2007 Table 3 (p708). The three
    # BBB terms of paper equation 9 are reported after the volume
    # aggregation of equation 11 (kdiff = Qdiff/Vecf, N*max = Nmax/Vecf,
    # keff = Qeff/Vecf), so they act directly on concentrations.
    #
    # Table 3 prints "N*max (ng min-1)", but equation 10 puts Nmax/Vecf
    # additively against kdiff * (Cblood - Cecf) inside dCecf/dt, so the
    # aggregated parameter must carry ng/mL/min. The abstract states the
    # units correctly as "(Nmax/Vecf) of 0.66 ng min-1 ml-1". The Table 3
    # header is the typo. See vignette Errata.
    # ------------------------------------------------------------------
    lkdiff_bbb   <- log(0.0014); label("Passive-diffusion rate constant across the BBB, kdiff (1/min)")            # Table 3: kdiff = 0.0014 /min (CV 12.6%)
    lkefflux_bbb <- log(0.0195); label("Active (Pgp-mediated) efflux rate constant out of brain ECF, keff (1/min), vehicle arm") # Table 3: keff, -GF120918 = 0.0195 /min (CV 12.2%)
    lNstarMax    <- log(0.658);  label("Volume-aggregated maximal active influx N*max = Nmax/Vecf (ng/mL/min)")    # Table 3: N*max = 0.658 (CV 26.1%); units per abstract, ng min-1 ml-1
    lc50         <- log(9.92);   label("Blood morphine concentration giving half-maximal active influx, C50 (ng/mL)") # Table 3: C50 = 9.92 ng/mL (CV 71.5%)

    # Paper equation 8 encodes the GF120918 effect as an arm switch,
    #     P_i = theta3 * (1 - GF120918_i) + theta4 * GF120918_i,
    # i.e. two independently estimated keff values (0.0195 vehicle,
    # 0.0113 elacridar). The canonical multiplicative-exponential form
    # below is numerically identical: CONMED_ELACRIDAR = 0 gives 0.0195
    # and CONMED_ELACRIDAR = 1 gives 0.0195 * (0.0113/0.0195) = 0.0113.
    # Same encoding as Xie_2000_m3g_rat.R's e_conmed_probenecid_cluin.
    e_conmed_elacridar_kefflux_bbb <- log(0.0113 / 0.0195); label("Exponential effect of CONMED_ELACRIDAR on keff (log(0.0113/0.0195) = -0.5455; the 42% reduction reported in the abstract)") # Table 3: keff, +GF120918 = 0.0113 /min (CV 25.4%) vs vehicle reference 0.0195 /min

    # ------------------------------------------------------------------
    # Inter-animal variability. Paper equations 3-4 give the log-normal
    # form P_i = P_typ * exp(eta_i), eta_i ~ N(0, omega^2), and the
    # "Interanimal variability" rows of Tables 2 and 3 are the omega^2
    # VARIANCES (the row labels are literally "o2 Cl", "o2 kdiff", ...).
    # Corroboration: every printed confidence interval in both tables is
    # reproduced by estimate +/- 1.96 * (CV%/100) * estimate on the
    # tabulated number itself, so the tabulated number is the estimated
    # variance and not its square root. IIV on all other parameters was
    # fixed to zero because it could not be estimated (paper Results).
    # ------------------------------------------------------------------
    etalcl ~ 0.129  # Table 2: omega^2 Cl = 0.129 (CV 17.2%, CI 0.085-0.173)
    etalvp ~ 0.099  # Table 2: omega^2 V2 = 0.099 (CV 24.7%, CI 0.051-0.147)

    # Table 3 reports a kdiff-keff covariance block. The third row's label
    # is mangled in the published PDF ("o2keffB o2keff"), but the Results
    # text states a covariance block was added because the post-hoc
    # individual kdiff and keff estimates correlated, and the value is a
    # covariance rather than a correlation: 0.059 +/- 1.96 * 0.850 * 0.059
    # = (-0.039, 0.158) reproduces the printed "0.039-0.158" exactly once
    # the dropped minus sign is restored, which is the same PDF symbol-font
    # defect that swallowed the minus in the kdiff (-0.025) and C50 (-3.98)
    # lower limits. The implied correlation is
    # 0.059 / sqrt(0.238 * 0.080) = 0.428, comfortably inside +/-1.
    etalkdiff_bbb + etalkefflux_bbb ~ c(0.238, 0.059, 0.080)  # Table 3: omega^2 kdiff = 0.238 (CV 56.3%), cov(kdiff, keff) = 0.059 (CV 85.0%), omega^2 keff = 0.080 (CV 36.8%)

    # ------------------------------------------------------------------
    # Residual error. Paper equations 5-6:
    #   Cobs_ij = Cpred_ij * (1 + eps_ij),  eps_ij ~ N(0, sigma^2)
    # so the "Proportional error" rows of Tables 2 and 3 are sigma^2 and
    # the proportional SD is its square root. The blood and brain ECF
    # fits are separate NONMEM runs with their own sigma^2. Same reading
    # as the sibling Geldof_2008_fluvoxamine_rat.R (shared author J
    # Freijer, same LAP&P/Leiden group and notation), where sigma^2 =
    # 0.042 maps to propSd = 0.2049.
    # ------------------------------------------------------------------
    propSd            <- 0.2720; label("Proportional residual SD on blood morphine concentrations Cc (fraction)")             # Table 2: sigma^2 = 0.074 (CV 10.2%); SD = sqrt(0.074) = 0.27203
    propSd_Cbrain_ecf <- 0.3066; label("Proportional residual SD on brain ECF morphine concentrations Cbrain_ecf (fraction)") # Table 3: sigma^2 = 0.094 (CV 21.4%); SD = sqrt(0.094) = 0.30659
  })

  model({
    # ---- Individual blood disposition parameters ----
    # Paper equation 7 is a centred LINEAR body-weight model, so the
    # covariate factor multiplies the exponentiated typical value rather
    # than entering the exponent. The reference weight 0.300 kg stands in
    # for the paper's unprinted `median BW` (vignette Errata).
    cl  <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 0.300))
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp + etalvp) * (1 + e_wt_vp * (WT - 0.300))
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- Individual BBB transport parameters ----
    kdiff_bbb   <- exp(lkdiff_bbb + etalkdiff_bbb)
    kefflux_bbb <- exp(lkefflux_bbb + e_conmed_elacridar_kefflux_bbb * CONMED_ELACRIDAR + etalkefflux_bbb)
    NstarMax    <- exp(lNstarMax)
    c50         <- exp(lc50)

    # ---- Blood concentration, the BBB driving force ----
    Cc <- central / vc

    # ---- Three-compartment blood disposition (paper Table 2) ----
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- Brain ECF (paper equations 9-11) ----
    # Equation 9 is the mass balance on the ECF amount Aecf; equation 10
    # divides through by Vecf and equation 11 names the aggregated
    # constants. The result acts on the ECF CONCENTRATION, which is why
    # brain_ecf is a concentration state here (see compartmentData).
    #   dCecf/dt = kdiff * (Cblood - Cecf)          passive diffusion
    #            + N*max * Cblood / (C50 + Cblood)  saturable active influx
    #            - keff * Cecf                      active Pgp efflux
    # Brain uptake is negligible relative to systemic disposition, so
    # there is no back-coupling term into d/dt(central); the paper
    # likewise uses the blood model purely as an input function.
    d/dt(brain_ecf) <- kdiff_bbb * (Cc - brain_ecf) +
                       NstarMax * Cc / (c50 + Cc) -
                       kefflux_bbb * brain_ecf

    # ---- Observations ----
    Cbrain_ecf <- brain_ecf  # unbound morphine in brain ECF, recovery-corrected microdialysate

    Cc         ~ prop(propSd)
    Cbrain_ecf ~ prop(propSd_Cbrain_ecf)
  })
}
