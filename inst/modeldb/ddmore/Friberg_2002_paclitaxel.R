Friberg_2002_paclitaxel <- function() {
  description <- "Semi-mechanistic Friberg-style myelosuppression PK/PD model for paclitaxel in adult cancer patients (Friberg 2002, leukocyte arm of DDMODEL00000186). Paclitaxel exposure is driven by per-subject empirical-Bayes PK estimates supplied as data columns (CL_INDIV, VC_INDIV, VP_INDIV) with intercompartmental clearance Q fixed at 204 L/h. Leukocyte response is described by a self-renewing proliferating pool plus three transit compartments and a circulating compartment, with a linear drug effect (1 - SLOPU * Cc) on proliferation and a feedback term (CIRC0 / circ)^GAMMA. Output is total circulating leukocytes in 10^9 cells/L."
  reference <- paste(
    "Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO. (2002).",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 20(24):4713-4721. doi:10.1200/JCO.2002.02.140 (PMID 12488418).",
    "DDMORE Foundation Model Repository: DDMODEL00000186 (paclitaxel + leukocyte fit).",
    sep = " "
  )
  vignette <- "Friberg_2002_paclitaxel"
  units <- list(time = "hour", dosing = "umol", concentration = "umol/L", leukocyte = "10^9/L")
  ddmore_id    <- "DDMODEL00000186"
  replicate_of <- NULL

  covariateData <- list(
    CL_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel clearance",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel clearance EBE supplied per-subject as a data column. Friberg 2002 fixed the paclitaxel popPK structure (3-h IV infusion 2-compartment model) and parameters from a previously-published popPK analysis and used POSTHOC EBEs as inputs to the myelosuppression model rather than re-estimating the PK and PD jointly. The DDMORE bundle's NM-TRAN data file ships these as column `CLI`. Reference values: median ~285 L/h, range ~160-540 L/h across the 46 virtual subjects in `Simulated_WBC_pacl_ddmore.csv`.",
      source_name        = "CLI"
    ),
    VC_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel central volume of distribution",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel central volume EBE supplied per-subject as a data column. Bundle column `V1I`. Reference values: median ~290 L in the bundle's simulated dataset.",
      source_name        = "V1I"
    ),
    VP_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel peripheral volume of distribution",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel peripheral volume EBE supplied per-subject as a data column. Bundle column `V2I`. Reference values: median ~995 L in the bundle's simulated dataset.",
      source_name        = "V2I"
    )
  )

  population <- list(
    n_subjects     = 46L,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adult cancer patients receiving paclitaxel chemotherapy. Friberg 2002 develops the myelosuppression model on six anticancer drugs (docetaxel, paclitaxel, etoposide, DMDC, CPT-11, vinflunine) for both neutrophils and leukocytes; the DDMORE bundle for DDMODEL00000186 implements only the paclitaxel + leukocyte fit on a subset of the paclitaxel cohort.",
    dose_range     = "Intravenous paclitaxel as a 3-hour infusion. The bundle's `Simulated_WBC_pacl_ddmore.csv` records doses ~190-530 umol per cycle (~160-450 mg of paclitaxel using MW 853.9 g/mol) with multiple cycles per subject over up to ~360 days follow-up.",
    regions        = NA_character_,
    notes          = "Population demographic detail (age, weight, sex, race) is not reproduced in the DDMORE bundle for DDMODEL00000186, and the original Friberg 2002 publication is not on disk in this worktree. n_subjects = 46 is the count of distinct IDs in the bundle's simulated dataset (`Simulated_WBC_pacl_ddmore.csv`); the bundle's NMTRAN re-fit on those data reports 45 individuals contributing to the eta shrinkage statistics, consistent with one subject having no usable observations after EVID filtering. The Friberg 2002 paper itself characterises the model on six drugs for both neutrophils and leukocytes; this DDMORE entry implements only the paclitaxel + leukocyte fit. See the validation vignette's Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # Structural PD parameters - Executable_Myelosuppression.{mdl,xml} STRUCTURAL block
    # mirrored as the rendered NMTRAN $THETA in Output_simulated_SEE_NONMEM.lst lines 99-104.
    # The bundle ships only an Output_simulated_*.lst (no Output_real_*.lst); however
    # re-fitting the .mod against the shipped simulated dataset reaches MINIMIZATION SUCCESSFUL
    # and recovers the same point values to three significant figures (Output_simulated_SEE_NONMEM.lst
    # lines 410-412 after MINIMIZATION SUCCESSFUL: TH1=7.21, TH2=123, TH3=0.238, TH4=28.7, TH5=0.286).
    # The values used here are therefore the publication-derived point estimates that the bundle
    # encodes, not refitted final estimates from a real-data run; see the vignette Errata.
    lcirc0 <- log(7.21);  label("Baseline circulating leukocyte count CIRC0 (10^9 cells/L)")             # POP_CIRC0 in .mdl PARAMETERS / STRUCTURAL
    lmtt   <- log(124);   label("Mean transit time MTT through the proliferation->circulation chain (h)") # POP_MTT
    lslopu <- log(28.9);  label("Linear drug-effect slope SLOPU on paclitaxel concentration (L/umol)")   # POP_SLOPU
    gamma  <- 0.239;      label("Feedback exponent gamma on (CIRC0 / circ) (unitless)")                  # POP_GAMMA - no IIV in source

    # Inter-individual variability - diagonal $OMEGA in the rendered NMTRAN
    # (Output_simulated_SEE_NONMEM.lst lines 106-109). Each PPV is reported as a variance on
    # the log-eta scale (`type is var` in .mdl VARIABILITY); no transformation needed.
    etalcirc0 ~ 0.107   # POP_OMEGA_CIRC0 in .mdl VARIABILITY (variance form)
    etalmtt   ~ 0.0296  # POP_OMEGA_MTT
    etalslopu ~ 0.176   # POP_OMEGA_SLOPU

    # Residual error - the .mdl uses proportionalError(proportional=PROP_ERROR, prediction=CIRC, eps=eps_ERROR)
    # with var(eps) FIXED at 1 (Output_simulated_SEE_NONMEM.lst $SIGMA 1.0 FIX line 112).
    # Rendered NMTRAN: W = PROP_ERROR * IPRED ; Y = IPRED + W * EPS(1) - i.e., proportional residual
    # error whose SD is reported as a fraction of CIRC.
    propSd <- 0.286;      label("Proportional residual error on circulating leukocyte count (fraction)")  # PROP_ERROR
  })

  model({
    # Individual myelosuppression parameters (PD layer - the only layer estimated in the bundle)
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt   + etalmtt)
    slopu <- exp(lslopu + etalslopu)

    # PK driver: paclitaxel pharmacokinetics is supplied via per-subject EBE PK estimates
    # (CL_INDIV, VC_INDIV, VP_INDIV) carried as data columns. Intercompartmental clearance Q is
    # hard-coded at 204 L/h per the .mod (`Q = 204` literal in $DES). This reproduces the DDMORE
    # encoding where the upstream paclitaxel popPK was held fixed and only the PD parameters were
    # estimated.
    q  <- 204
    Cc <- central / VC_INDIV  # paclitaxel central concentration (umol/L)

    # Drug effect on proliferation rate and feedback from circulating cells
    edrug <- 1 - slopu * Cc
    feed  <- (circ0 / circ)^gamma

    # Transit-rate constant; the chain has 3 transit compartments (NN = 3 in the .mod) plus a
    # single proliferation compartment, so KTR = (NN + 1) / MTT = 4 / MTT.
    ktr <- 4 / mtt

    # Compartment ordering matches the bundle's NMTRAN $MODEL block so that the bundle's
    # `Simulated_WBC_pacl_ddmore.csv` (CMT = 1 doses, CMT = 3 leukocyte observations) maps
    # directly onto rxode2's positional compartment numbering: central = 1, peripheral1 = 2,
    # circ = 3, precursor1 (PROL) = 4, precursor2..precursor4 (TRANSIT1..3) = 5..7.
    # Paclitaxel two-compartment IV PK driven by the per-subject EBE columns
    d/dt(central)     <- -q / VC_INDIV * central - CL_INDIV / VC_INDIV * central + q / VP_INDIV * peripheral1
    d/dt(peripheral1) <-  q / VC_INDIV * central - q / VP_INDIV * peripheral1
    # Friberg myelosuppression chain - circ before precursors so CMT = 3 lands on circ.
    d/dt(circ)        <-  ktr * precursor4 - ktr * circ
    d/dt(precursor1)  <-  ktr * precursor1 * edrug * feed - ktr * precursor1
    d/dt(precursor2)  <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3)  <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4)  <-  ktr * precursor3 - ktr * precursor4

    # Initial conditions: all five myelosuppression compartments start at the individual baseline.
    circ(0)       <- circ0
    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0

    # Observation: total circulating leukocytes (10^9 cells/L)
    WBC <- circ
    WBC ~ prop(propSd)
  })
}
