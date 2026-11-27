Calvier_2018_iv1cmt_dispersion <- function() {
  description <- paste(
    "Methodology reference. Two-drug population PK design model from",
    "Calvier 2018, the optimal-design study of whether the two hepatic",
    "PBPK parameters whole-liver unbound intrinsic clearance (CLint,u,WL)",
    "and hepatic blood flow (Qh) can be estimated from clinical PK data in",
    "children. There is no real molecule and no real patients: drugs A and",
    "B are hypothetical substrates of one isoenzyme, each given as a",
    "constant-rate IV infusion to its own arm of 1-year-old children, and",
    "every value is an author-chosen design constant encoded with fixed().",
    "Each drug follows a one-compartment model whose plasma clearance is",
    "the hepatic dispersion model (dispersion number 0.17). Both drugs",
    "share one Qh and one CLint,u,WL with their between-subject",
    "variability; drug B's intrinsic clearance is drug A's times a fixed",
    "in-vitro ratio (CLintratio), which is what makes Qh and CLint,u,WL",
    "jointly identifiable. Typical values are the stochastic",
    "simulation-and-estimation pair of the paper (drug A extraction ratio",
    "0.94, drug B 0.20) at the Simcyp-derived 1-year-old Qh of 20.3 L/h;",
    "the omegas are the deposited PFIM scenario (0.1 on Qh and on",
    "CLint,u,WL, 0.25 on each volume).",
    sep = " "
  )
  reference <- paste(
    "Calvier EAM, Nguyen TT, Johnson TN, Rostami-Hodjegan A, Tibboel D,",
    "Krekels EHJ, Knibbe CAJ. Can Population Modelling Principles be Used",
    "to Identify Key PBPK Parameters for Paediatric Clearance Predictions?",
    "An Innovative Application of Optimal Design Theory. Pharm Res.",
    "2018;35(11):209. doi:10.1007/s11095-018-2487-1. PMCID PMC6156772.",
    "Structure and design constants follow the authors' deposited PFIM 4.0",
    "model file (Supplementary material 2) and input file (Supplementary",
    "material 3, scenario 'OMEGA QH = 0.1 & OMEGA CLINT = 0.1'), which run",
    "in minutes and are converted to hours here; the two hypothetical drugs",
    "are rows of Supplementary table 1 with extraction ratios 0.94 and 0.20.",
    sep = " "
  )
  vignette <- "Calvier_2018_iv1cmt_dispersion"
  dosing <- c("central", "central_drugb")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Both drugs are hypothetical, so the analyte names carry the paper's
  # A / B labels. PFIM simulates and observes total plasma concentrations
  # (CLp is "the total (bound and unbound) plasma clearance", equation 1).
  compartmentData <- list(
    central = list(
      analyte = "hypothetical drug A",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_drugb = list(
      analyte = "hypothetical drug B",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "None (methodology paper; hypothetical drugs and a virtual paediatric design population, no fitted estimates).",
    n_subjects = 90L,
    n_studies = 1L,
    age_range = "1 year (the design population; Qh taken from 500 male and 500 female Simcyp V15 1-year-olds)",
    disease_state = "N/A (optimal-design and stochastic simulation-estimation study; not a fit of any real molecule).",
    dose_range = paste(
      "Constant-rate IV infusion started at t = 0 with the rate set to reach",
      "a steady-state plasma concentration of 70 mg/L (rate = CLp * 70);",
      "drug A to one arm and drug B to the other."
    ),
    regions = "N/A",
    scope_note = paste(
      "Filed under inst/modeldb/pharmacokinetics/ (not specificDrugs/)",
      "because there is no drug: the paper is a methodology reference whose",
      "PK content is a pair of author-defined hypothetical drugs."
    ),
    notes = paste(
      "Design (Clinical Trial Design; Fig. 2): two arms of 45 patients each,",
      "seven samples per patient taken every half-life from the first to the",
      "seventh half-life after the start of the infusion. fu and B:P are 1",
      "for both drugs ('In order to reduce computational costs ... fu and",
      "B:P were taken to be 1'). Step 2 swept 99 hypothetical drugs (ER",
      "0.01-0.99, Supplementary table 1) over four omega scenarios",
      "(omega^2 CLint,u,WL 0.1 or 0.8; omega^2 Qh 0.1 or 0.6); the",
      "deposited PFIM input carries the 0.1 / 0.1 scenario encoded here.",
      "Step 3 (stochastic simulation and estimation, NONMEM 7.3) used one",
      "drug with ER 0.94 and one with ER 0.2 in 1-year-olds, the pair",
      "encoded here."
    )
  )

  ini({
    lqh <- fixed(log(20.30106)); label("Hepatic blood flow Qh (L/h)") # Supp. material 3 beta QH = 0.338351 L/min * 60 = 20.30 L/h; Table I 'Theta Qh 20.3 L/h'; Table III 1 year 20.3
    lclint <- fixed(log(80.31137)); label("Whole-liver unbound intrinsic clearance of drug A, CLint,u,WL (L/h)") # Supp. table 1 row ER 0.940000027: clint 1.338522914 L/min * 60; within Table I range 0.10-156 L/h
    clint_ratio <- fixed(0.05818729); label("CLintratio = CLint,u,WL of drug B / CLint,u,WL of drug A (unitless)") # eq 7; Supp. table 1 clint 0.077884914 (ER 0.199999899) / 1.338522914 (ER 0.940000027); within Table I range 0.00065-0.97
    dn <- fixed(0.17); label("Dispersion number DN of the hepatic dispersion model (unitless)") # Methods text after eq 6 'taken to be 0.17'; Table I Theta DN 0.17 Fixed; Supp. material 2 Dn <- 0.17
    fu <- fixed(1); label("Fraction unbound in plasma of drug A (unitless)") # Table I Theta fu 1 Fixed; Supp. material 2 FuA <- 1
    fu_drugb <- fixed(1); label("Fraction unbound in plasma of drug B (unitless)") # Table I Theta fu 1 Fixed; Supp. material 2 FuB <- 1
    bp <- fixed(1); label("Blood-to-plasma ratio B:P of drug A (unitless)") # Table I Theta B:P 1 Fixed; Supp. material 2 BPA <- 1
    bp_drugb <- fixed(1); label("Blood-to-plasma ratio B:P of drug B (unitless)") # Table I Theta B:P 1 Fixed; Supp. material 2 BPB <- 1
    lvc <- fixed(log(4)); label("Volume of distribution of drug A, VA (L)") # Table I Theta VA 4 L; Supp. material 3 beta VA = 4
    lvc_drugb <- fixed(log(4)); label("Volume of distribution of drug B, VB (L)") # Table I Theta VB 4 L; Supp. material 3 beta VB = 4

    etalqh ~ fixed(0.1) # Supp. material 3 omega diag QH = 0.1 (project 'OMEGA QH = 0.1 & OMEGA CLINT = 0.1'); Table I range 0.1-0.6
    etalclint ~ fixed(0.1) # Supp. material 3 omega diag CLINTA = 0.1; Table I range 0.1-0.8
    etalvc ~ fixed(0.25) # Table I omega^2 VA 0.25; Supp. material 3 omega diag VA = 0.25
    etalvc_drugb ~ fixed(0.25) # Table I omega^2 VB 0.25; Supp. material 3 omega diag VB = 0.25

    propSd <- fixed(0.1); label("Proportional residual error SD, drug A (fraction)") # Table I sigma A 0.1; Supp. material 3 sig.slopeA = 0.1, sig.interA = 0
    propSd_drugb <- fixed(0.1); label("Proportional residual error SD, drug B (fraction)") # Table I sigma B 0.1; Supp. material 3 sig.slopeB = 0.1, sig.interB = 0
  })

  model({
    # System-specific Qh and the drug-A intrinsic clearance carry the
    # between-subject variability; drug B's intrinsic clearance is the same
    # individual value times the fixed in-vitro ratio (equations 8 and 9).
    qh <- exp(lqh + etalqh)
    clint <- exp(lclint + etalclint)
    clint_drugb <- clint * clint_ratio
    vc <- exp(lvc + etalvc)
    vc_drugb <- exp(lvc_drugb + etalvc_drugb)

    # Dispersion model (equations 1-6): efficiency number rn, hepatic
    # availability fh, extraction ratio eh, total plasma clearance
    # CLp = Qh * ER * B:P. Written exactly as the deposited PFIM model.
    rn <- (fu / bp) * clint / qh
    a_disp <- sqrt(1 + 4 * rn * dn)
    fh <- 4 * a_disp / ((1 + a_disp)^2 * exp((a_disp - 1) / (2 * dn)) - (1 - a_disp)^2 * exp(-(a_disp + 1) / (2 * dn)))
    eh <- 1 - fh
    cl <- qh * eh * bp

    rn_drugb <- (fu_drugb / bp_drugb) * clint_drugb / qh
    a_disp_drugb <- sqrt(1 + 4 * rn_drugb * dn)
    fh_drugb <- 4 * a_disp_drugb / ((1 + a_disp_drugb)^2 * exp((a_disp_drugb - 1) / (2 * dn)) - (1 - a_disp_drugb)^2 * exp(-(a_disp_drugb + 1) / (2 * dn)))
    eh_drugb <- 1 - fh_drugb
    cl_drugb <- qh * eh_drugb * bp_drugb

    kel <- cl / vc
    kel_drugb <- cl_drugb / vc_drugb

    d/dt(central) <- -kel * central
    d/dt(central_drugb) <- -kel_drugb * central_drugb

    Cc <- central / vc
    Cc_drugb <- central_drugb / vc_drugb

    Cc ~ prop(propSd)
    Cc_drugb ~ prop(propSd_drugb)
  })
}
