Stein_2018_romosozumab <- function() {
  description <- "Two-compartment QSS TMDD typical-value fit for romosozumab (anti-sclerostin mAb) used to illustrate the critical concentration (Ccrit) for nonlinear PK (Stein and Peletier 2018 Table 1)"
  reference <- "Stein AM, Peletier LA. Predicting the onset of nonlinear pharmacokinetics. CPT Pharmacometrics Syst Pharmacol. 2018 Oct;7(10):670-677. doi:10.1002/psp4.12316. Underlying romosozumab PK data from Padhi D, Jang G, Stouch B, Fang L, Posvar E. J Bone Miner Res. 2011;26(1):19-26."
  vignette <- "Stein_2018_mAb_nonlinear_PK"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = 1,
    disease_state = "Healthy adults (romosozumab single-dose, placebo-controlled phase 1 study; Padhi 2011).",
    dose_range    = "Stein and Peletier 2018 simulated a 10 mg/kg single IV bolus (467 nmol assuming 70 kg patient and 150 kDa antibody; Stein and Peletier 2018 page 672).",
    regions       = NA_character_,
    notes         = paste(
      "Typical-value fit reproduced from Stein and Peletier 2018 Table 1. Parameters were estimated by",
      "Stein and Peletier (using the two-compartment quasi-steady-state TMDD model in their Eq. 7) by",
      "fitting the romosozumab PK profiles from Padhi D et al. J Bone Miner Res 2011;26:19-26.",
      "Stein and Peletier 2018 note (page 672, Model fit) that the large values for ke(R) and ke(CR)",
      "for romosozumab reflect the practical unidentifiability of these parameters in this model.",
      "No IIV or residual error was reported in the article; the model is provided as a typical-value",
      "(deterministic) simulation tool for the C_crit illustration in Stein and Peletier 2018 Figure 1."
    )
  )

  ini({
    # Drug disposition (Stein and Peletier 2018 Table 1, Romosozumab column)
    lvc   <- fixed(log(2.4));  label("Central volume of distribution (Vc, L)")                  # Stein and Peletier 2018 Table 1: Vc = 2.4 L
    lvp   <- fixed(log(2.6));  label("Peripheral volume of distribution (Vp, L)")               # Stein and Peletier 2018 Table 1: Vp = 2.6 L
    lcl   <- fixed(log(0.25)); label("Linear (nonspecific) clearance (CL, L/day)")              # Stein and Peletier 2018 Table 1: CL = 0.25 L/d
    lq    <- fixed(log(0.54)); label("Intercompartmental clearance (Q, L/day)")                 # Stein and Peletier 2018 Table 1: Q = 0.54 L/d

    # Target turnover and binding (QSS approximation; paper Eq. 7)
    lksyn <- fixed(log(6.1));  label("Target synthesis rate (ksyn = Vmax/Vc, nM/day)")          # Stein and Peletier 2018 Table 1: ksyn = Vmax/Vc = 6.1 nM/d
    lKss  <- fixed(log(12));   label("QSS binding constant (Kss = KM, nM)")                     # Stein and Peletier 2018 Table 1: Kss = KM = 12 nM
    lkdeg <- fixed(log(860));  label("Free target elimination rate (ke(R), 1/day)")             # Stein and Peletier 2018 Table 1: ke(R) = 860 1/d (practically unidentifiable; Stein and Peletier 2018 page 672)
    lkint <- fixed(log(860));  label("Drug-target complex internalization rate (ke(CR), 1/day)") # Stein and Peletier 2018 Table 1: ke(CR) = 860 1/d (assumed equal to ke(R); Stein and Peletier 2018 page 672)
  })

  model({
    # Typical-value parameters (no IIV reported in the source)
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    cl   <- exp(lcl)
    q    <- exp(lq)
    ksyn <- exp(lksyn)
    kss  <- exp(lKss)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)

    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    t0   <- ksyn / kdeg

    total_target(0) <- t0

    # QSS approximation in the central compartment (paper Eq. 5).
    ctot    <- central / vc
    disc    <- ctot - total_target - kss
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # ODEs (paper Eq. 7; central and peripheral1 are nmol amounts; total_target is nM).
    d/dt(central)      <- -kel * cfree * vc - kint * complex * vc - k12 * cfree * vc + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * cfree * vc - k21 * peripheral1
    d/dt(total_target) <-  ksyn - kdeg * (total_target - complex) - kint * complex

    # Stein and Peletier 2018 Figure 1 reports total drug concentration on the y-axis.
    Cc <- ctot
  })
}
