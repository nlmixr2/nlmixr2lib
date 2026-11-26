Radulescu_2019_lead_qsp <- function() {
  description <- paste(
    "QSP. Phenomenological systems model of lead (Pb) and calcium (Ca)",
    "competitive toxicokinetics, focused on neurotoxicity. Ten ODE states",
    "carry the lead and calcium content of five compartments: blood (S),",
    "bone (O), brain (B), kidney (K) and other soft tissue (T). Transfer",
    "between compartments is a saturable competitive flux: a shared maximum",
    "rate A, split between Pb and Ca by their fractional abundance, times a",
    "sigmoidal 'crowdedness' nonlinearity N(x)=1/(1+exp(-b(x-theta))) on the",
    "joint Pb+Ca pool. Three mechanisms add feedback: gastrointestinal",
    "absorption is gated down by blood calcium (1-N_{a,tau}(S_Ca)); bone",
    "resorption is scaled by z (an age / physiological-state marker); renal",
    "excretion is blocked by accumulated kidney lead (exp(-k*K_Pb)); and the",
    "blood-brain barrier gates the blood-to-brain flux through",
    "xi=1-B*exp(-s*B_Pb), where B is the BBB filtering efficiency and s its",
    "sensitivity to lead neurotoxicity. Deterministic (no IIV, no residual",
    "error); driven by constant lead and calcium intake rates p and c and by",
    "the initial compartment contents, which determine which steady state is",
    "reached inside the model's bistability window. All transport rates A are",
    "normalised to 1. Time in weeks; Pb/Ca content in mg/kg (in blood,",
    "1 mg/kg is approximately 1 mg/L, so blood lead in mg/L x 100 = ug/dl)."
  )
  reference <- paste(
    "Radulescu A, Lundgren S. A pharmacokinetic model of lead absorption and",
    "calcium competitive dynamics. Sci Rep. 2019;9(1):14225.",
    "doi:10.1038/s41598-019-50654-7.",
    sep = " "
  )
  vignette <- "Radulescu_2019_lead_qsp"

  paper_specific_compartments <- c(
    "blood_Pb",
    "blood_Ca",
    "bone_Pb",
    "bone_Ca",
    "brain_Pb",
    "brain_Ca",
    "kidney_Pb",
    "kidney_Ca",
    "soft_Pb",
    "soft_Ca"
  )

  units <- list(time = "week", dosing = "mg/kg", concentration = "mg/kg")

  # Each ODE state holds the lead (Pb) or calcium (Ca) content of one
  # compartment, in mg/kg (Methods, 'unit convention'). Bone, brain, kidney and
  # soft tissue are solid matrices mapped to 'tissue'; blood lead/calcium is the
  # clinically measured whole-blood content. verified = TRUE: the paper names
  # each compartment's Pb and Ca content explicitly (Fig. 2, Key Compartments).
  compartmentData <- list(
    blood_Pb = list(analyte = "Lead", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    blood_Ca = list(analyte = "Calcium", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    bone_Pb = list(analyte = "Lead", units = "mg/kg", specimen = "tissue", verified = TRUE),
    bone_Ca = list(analyte = "Calcium", units = "mg/kg", specimen = "tissue", verified = TRUE),
    brain_Pb = list(analyte = "Lead", units = "mg/kg", specimen = "tissue", verified = TRUE),
    brain_Ca = list(analyte = "Calcium", units = "mg/kg", specimen = "tissue", verified = TRUE),
    kidney_Pb = list(analyte = "Lead", units = "mg/kg", specimen = "tissue", verified = TRUE),
    kidney_Ca = list(analyte = "Calcium", units = "mg/kg", specimen = "tissue", verified = TRUE),
    soft_Pb = list(analyte = "Lead", units = "mg/kg", specimen = "tissue", verified = TRUE),
    soft_Ca = list(analyte = "Calcium", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "Environmental lead exposure; neurotoxicity focus",
    dose_range = "Lead intake p and calcium intake c varied (p ~ 0-2, c ~ 0-20 mg/kg/week, with c > p)",
    notes = paste(
      "Phenomenological / theoretical model, not fitted to an individual",
      "cohort. Parameters were loosely tuned so the state variables lie within",
      "empirically observed human ranges (Methods, Results); the authors cite",
      "single-subject and group lead-exposure datasets (Moore, Sherlock,",
      "Williams, Cools, O'Flaherty) for qualitative and quantitative",
      "comparison. Age / developmental stage enters through the bone-resorption",
      "marker z and the blood-brain-barrier efficiency B and sensitivity s:",
      "an immature (child) BBB is low B / high s, a mature (adult) BBB is high",
      "B / low s. Baseline scenario in this file is a healthy adult",
      "(B = 0.8, s = 1, z = 1, k = 0.05) at medium lead intake p = 0.5 and",
      "calcium intake c = 10."
    )
  )

  ini({
    # All values are the baseline scenario fixed throughout Radulescu 2019
    # (Methods and the Fig. 3-5 captions). Nothing is estimated: this is a
    # deterministic phenomenological model, so every parameter is fixed().

    # Intake rates (environmental source -> blood), mg/kg/week.
    pbIntake <- fixed(0.5)
    label("Lead intake p, medium exposure (mg/kg/week)") # Fig. 3-5 captions
    caIntake <- fixed(10)
    label("Calcium intake c (mg/kg/week)") # Fig. 3-5 captions

    # Maximum inter-compartmental transport rate A_XY, mg/kg/week. The paper
    # normalises every A_XY to 1 (Methods, 'unit convention').
    aMaxRate <- fixed(1)
    label("Maximum inter-compartmental transport rate A (mg/kg/week)") # Methods, normalised to 1

    # Transport 'crowdedness' sigmoidal N_{b,theta}(x) = expit(b*(x-theta)).
    bSens <- fixed(0.6)
    label("Transport sigmoidal sensitivity b") # Methods; Appendix A robustness
    thetaSig <- fixed(6)
    label("Transport sigmoidal centre theta (mg/kg)") # Methods; Appendix A robustness

    # Absorption calcium-gating sigmoidal, 1 - N_{a,tau}(S_Ca).
    aGate <- fixed(0.6)
    label("Absorption gating sensitivity a") # Absorption into blood; Appendix A
    tauGate <- fixed(4)
    label("Absorption gating centre tau on blood calcium (mg/kg)") # Absorption into blood; Appendix A

    # Bone resorption multiplier z (age / physiological-state marker).
    zResorp <- fixed(1)
    label("Bone resorption multiplier z, healthy adult") # Fig. 3-5 captions

    # Renal sensitivity to lead toxicity k in the excretion feedback
    # exp(-k*K_Pb).
    kRenalSens <- fixed(0.05)
    label("Kidney sensitivity to lead toxicity k (kg/mg)") # Fig. 3-5 captions

    # Blood-brain-barrier gate xi = 1 - B*exp(-s*B_Pb).
    bbbEff <- fixed(0.8)
    label("BBB filtering efficiency B (0 to 1)") # Fig. 3-5 captions
    sBbbSens <- fixed(1)
    label("BBB sensitivity to lead neurotoxicity s (kg/mg)") # Fig. 3-5 captions

    # Initial compartment content used for every state, mg/kg. In the
    # bistability window the chosen steady state depends on these initial
    # conditions; the baseline 'low start' scenario (Fig. 5 top) sets all
    # contents to 0.1. Override individual states (e.g. an elevated initial
    # blood lead) via rxSolve(inits=) to explore the other basin (Fig. 5
    # bottom raises blood_Pb to 4).
    initContent <- fixed(0.1)
    label("Low-start initial content, all compartments (mg/kg)") # Fig. 5 top
  })

  model({
    # Joint Pb+Ca pools per compartment (the sigmoidal 'crowdedness' argument).
    sSum <- blood_Pb + blood_Ca
    tSum <- soft_Pb + soft_Ca
    oSum <- bone_Pb + bone_Ca
    kSum <- kidney_Pb + kidney_Ca
    bSumm <- brain_Pb + brain_Ca
    pcSum <- pbIntake + caIntake

    # Sigmoidal crowdedness nonlinearity N_{b,theta}(x) = expit(b*(x-theta)),
    # evaluated on each joint pool.
    nS <- expit(bSens * (sSum - thetaSig))
    nT <- expit(bSens * (tSum - thetaSig))
    nO <- expit(bSens * (oSum - thetaSig))
    nK <- expit(bSens * (kSum - thetaSig))
    nB <- expit(bSens * (bSumm - thetaSig))
    nPC <- expit(bSens * (pcSum - thetaSig))

    # Absorption calcium gate 1 - N_{a,tau}(S_Ca): high blood calcium slows
    # lead + calcium absorption from the gut.
    absGate <- 1 - expit(aGate * (blood_Ca - tauGate))

    # Renal excretion feedback phi = exp(-k*K_Pb): accumulated kidney lead
    # blocks excretion of both elements.
    phiRenal <- exp(-kRenalSens * kidney_Pb)

    # Blood-brain-barrier gate xi = 1 - B*exp(-s*B_Pb): brain lead erodes the
    # barrier's ability to discriminate lead from calcium.
    xiBBB <- 1 - bbbEff * exp(-sBbbSens * brain_Pb)
    denomB <- xiBBB * blood_Pb + (1 - xiBBB) * blood_Ca

    # Environment -> blood absorption (gated by blood calcium).
    fES_p <- aMaxRate * pbIntake / pcSum * nPC * absGate
    fES_c <- aMaxRate * caIntake / pcSum * nPC * absGate

    # Blood <-> other soft tissue.
    fST_p <- aMaxRate * blood_Pb / sSum * nS
    fST_c <- aMaxRate * blood_Ca / sSum * nS
    fTS_p <- aMaxRate * soft_Pb / tSum * nT
    fTS_c <- aMaxRate * soft_Ca / tSum * nT

    # Blood -> bone (deposit) and bone -> blood (resorption, scaled by z).
    fSO_p <- aMaxRate * blood_Pb / sSum * nS
    fSO_c <- aMaxRate * blood_Ca / sSum * nS
    fOS_p <- aMaxRate * zResorp * bone_Pb / oSum * nO
    fOS_c <- aMaxRate * zResorp * bone_Ca / oSum * nO

    # Blood <-> kidney, plus renal excretion with the exp(-k*K_Pb) feedback.
    fSK_p <- aMaxRate * blood_Pb / sSum * nS
    fSK_c <- aMaxRate * blood_Ca / sSum * nS
    fKS_p <- aMaxRate * kidney_Pb / kSum * nK
    fKS_c <- aMaxRate * kidney_Ca / kSum * nK
    fKE_p <- aMaxRate * phiRenal * kidney_Pb / kSum * nK
    fKE_c <- aMaxRate * phiRenal * kidney_Ca / kSum * nK

    # Blood -> brain (gated by the BBB) and brain -> blood (ungated).
    fSB_p <- aMaxRate * xiBBB * blood_Pb / denomB * nS
    fSB_c <- aMaxRate * (1 - xiBBB) * blood_Ca / denomB * nS
    fBS_p <- aMaxRate * brain_Pb / bSumm * nB
    fBS_c <- aMaxRate * brain_Ca / bSumm * nB

    # Lead-Calcium model (Radulescu 2019, 'Lead-Calcium model' system).
    d/dt(blood_Pb) <- fES_p + fTS_p + fKS_p + fOS_p + fBS_p -
      (fST_p + fSK_p + fSO_p + fSB_p)
    d/dt(blood_Ca) <- fES_c + fTS_c + fKS_c + fOS_c + fBS_c -
      (fST_c + fSK_c + fSO_c + fSB_c)
    d/dt(soft_Pb) <- fST_p - fTS_p
    d/dt(soft_Ca) <- fST_c - fTS_c
    d/dt(bone_Pb) <- fSO_p - fOS_p
    d/dt(bone_Ca) <- fSO_c - fOS_c
    d/dt(kidney_Pb) <- fSK_p - fKS_p - fKE_p
    d/dt(kidney_Ca) <- fSK_c - fKS_c - fKE_c
    d/dt(brain_Pb) <- fSB_p - fBS_p
    d/dt(brain_Ca) <- fSB_c - fBS_c

    # Low-start initial conditions (Fig. 5 top); override via rxSolve(inits=).
    blood_Pb(0) <- initContent
    blood_Ca(0) <- initContent
    soft_Pb(0) <- initContent
    soft_Ca(0) <- initContent
    bone_Pb(0) <- initContent
    bone_Ca(0) <- initContent
    kidney_Pb(0) <- initContent
    kidney_Ca(0) <- initContent
    brain_Pb(0) <- initContent
    brain_Ca(0) <- initContent

    # Clinically measured readout: whole-blood lead content (mg/kg approx mg/L;
    # x 100 gives ug/dl). The brain lead/calcium states are the paper's primary
    # neurotoxicity endpoints and are returned directly by rxSolve.
    Cc <- blood_Pb
  })
}
