vanOs_2024_temocillin_isolmic16 <- function() {
  description <- paste0(
    "In vitro (Escherichia coli ISOL MIC16, temocillin MIC 16 mg/L). ",
    "Semi-mechanistic pharmacokinetic/pharmacodynamic model of temocillin bacterial killing and ",
    "resistance amplification, fitted jointly to static time-kill and hollow-fibre infection model data ",
    "and coupled, for simulation, to the two-compartment population pharmacokinetic model of temocillin ",
    "in critically ill patients of Laterre et al. 2015. Three bacterial subpopulations grow logistically ",
    "towards a shared maximum density and are killed by the unbound plasma concentration: a susceptible ",
    "and a less susceptible subpopulation (sigmoidal Emax kill, together making up the count on drug-free ",
    "agar) and a resistant subpopulation growing on agar containing 32 mg/L temocillin (linear kill), ",
    "which does not contribute to the total count. Sibling models fitted to the other three strains of ",
    "the same paper: vanOs_2024_temocillin_atcc25922, vanOs_2024_temocillin_isolmic8, vanOs_2024_temocillin_isolmic4."
  )
  reference <- paste0(
    "van Os W, Nussbaumer-Proll A, Pham AD, Wijnant GJ, Ngougni Pokem P, Van Bambeke F, ",
    "van Hasselt JGC, Zeitlinger M. Pharmacokinetic/pharmacodynamic model-based optimization of ",
    "temocillin dosing strategies for the treatment of systemic infections. ",
    "J Antimicrob Chemother. 2024;79(10):2484-2492. doi:10.1093/jac/dkae243. ",
    "Structural equations: main-text Equations 1-3. ",
    "Pharmacodynamic parameter estimates: main-text Table 3, column ISOL MIC16. ",
    "Strain characteristics: main-text Table 1. ",
    "Population pharmacokinetic parameters: main-text Table 2, reproduced from ",
    "Laterre PF, Wittebole X, Van de Velde S, Muller AE, Mouton JW, Carryn S, Tulkens PM, Dugernier T. ",
    "Temocillin (6 g daily) in critically ill patients: continuous infusion versus three times daily ",
    "administration. J Antimicrob Chemother. 2015;70(3):891-898. doi:10.1093/jac/dku465. ",
    "Structure and coupling verified against the Supplementary data file NONMEM_code.txt."
  )
  vignette <- "vanOs_2024_temocillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The three bacterial subpopulation states are paper-mechanistic and follow the
  # registry-wide "^bact_" paper-specific compartment family.
  paper_specific_compartment_pattern <- "^bact_"

  compartmentData <- list(
    central = list(analyte = "temocillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "temocillin", units = "mg", specimen = "plasma", verified = TRUE),
    bact_susceptible = list(
      analyte = "Escherichia coli ISOL MIC16 (susceptible subpopulation)",
      units = "cfu/mL", specimen = "not applicable", verified = TRUE
    ),
    bact_less_susceptible = list(
      analyte = "Escherichia coli ISOL MIC16 (less susceptible subpopulation)",
      units = "cfu/mL", specimen = "not applicable", verified = TRUE
    ),
    bact_resistant = list(
      analyte = "Escherichia coli ISOL MIC16 (resistant subpopulation growing on agar with 32 mg/L temocillin)",
      units = "cfu/mL", specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Escherichia coli ISOL MIC16); pharmacokinetics from critically ill adults",
    n_subjects = 1L,
    n_studies = 1L,
    strain = paste0(
      "ISOL MIC16 (catheter urine; ESBL genes: blaCTX-M-15, blaOXA-1; temocillin MIC 16 mg/L, determined in triplicate ",
      "by broth microdilution in cation-adjusted Mueller-Hinton broth following CLSI guidelines) (Table 1)"
    ),
    disease_state = paste0(
      "Escherichia coli selected to span the less susceptible end of the wild-type MIC distribution for ",
      "the species (4-16 mg/L). The pharmacokinetic layer describes critically ill patients, from the ",
      "population pharmacokinetic study of Laterre et al. 2015 (n = 11)."
    ),
    model_system = paste0(
      "Static time-kill (STK) experiments in 5 mL cation-adjusted Mueller-Hinton broth at 37 degrees C, ",
      "in triplicate, inoculated at 1.5e6 cfu/mL, with temocillin at 0.125-8x MIC in 2-fold steps plus a ",
      "growth control, sampled over 24 h; and a hollow-fibre infection model (HFIM; FX paed dialysis ",
      "cartridge, 50 mL extracapillary space) in which the unbound plasma concentration-time profiles of ",
      "four intravenous regimens were replicated over 72 h. Total counts were read on drug-free Columbia ",
      "agar and the resistant subpopulation on Mueller-Hinton agar containing 32 mg/L temocillin. The ",
      "theoretical limit of detection was 50 cfu/mL; observations below it were handled with the M3 method."
    ),
    dose_range = paste0(
      "HFIM regimens replicated over 72 h: 2 g q12h and 2 g q8h intermittent infusion, and 4 g/day and ",
      "6 g/day continuous infusion each with a 2 g loading dose; all infusions over 30 min. Published ",
      "simulations spanned 2 g q12h up to 12 g/day continuous infusion."
    ),
    initial_inoculum = paste0(
      "The published simulations set the total bacterial population at time zero to 1e6 cfu/mL and scaled ",
      "each subpopulation by the same factor, preserving the fitted subpopulation ratios (Methods, PK/PD ",
      "simulations). That rescaling is reproduced inside model() from log10cfu0_total, so the Table 3 ",
      "estimates below remain the traceable fitted values."
    ),
    notes = paste0(
      "Estimated in NONMEM 7.4 with Laplacian estimation, PsN 5.3.0 and Pirana 21.11.1. The model was ",
      "built in three steps: STK total counts first, then STK and HFIM data jointly, then the resistant ",
      "subpopulation observed on agar containing 32 mg/L temocillin. The maximum bacterial density was ",
      "estimated separately for the two experimental systems (Bmax,STK 8.86 and Bmax,HFIM ",
      "9.94 log10 cfu/mL for this strain); log10bmax below carries the HFIM value because that ",
      "is the value the published simulations used. No inter-individual variability was estimated on the ",
      "pharmacodynamic parameters; the only random effects are the pharmacokinetic ones inherited from ",
      "Laterre et al. 2015. The published Monte Carlo simulations (n = 1000 patients per regimen) ",
      "propagated that pharmacokinetic inter-individual variability but not the residual variability of ",
      "either model."
    )
  )

  ini({
    # ---- Population pharmacokinetics, critically ill patients ---------------
    # Table 2, reproduced from Laterre et al. 2015. Not estimated in this paper,
    # hence fixed(); confirmed as hard-coded constants in NONMEM_code.txt $PK.
    lcl <- fixed(log(3.69)); label("Clearance (L/h)")                       # Table 2: CL 3.69 L/h
    lvc <- fixed(log(14.0)); label("Central volume of distribution (L)")    # Table 2: Vc 14.0 L
    lq  <- fixed(log(8.45)); label("Intercompartmental clearance (L/h)")    # Table 2: Q 8.45 L/h
    lvp <- fixed(log(21.7)); label("Peripheral volume of distribution (L)") # Table 2: Vp 21.7 L
    fu  <- fixed(0.41);      label("Unbound fraction of temocillin in plasma") # Table 2: fu 0.41, the mean unbound fraction observed in the PK study; footnote b gives 0.25 and 0.57 as mean -/+ 1 SD

    # ---- Bacterial growth ---------------------------------------------------
    kg_s   <- 1.39; label("Growth rate constant of the susceptible subpopulation (1/h)")      # Table 3 (ISOL MIC16): k_g,S 1.39, %RSE 7.0
    kg_ls  <- 0.795; label("Growth rate constant of the less susceptible subpopulation (1/h)") # Table 3 (ISOL MIC16): k_g,LS 0.795, %RSE 11.3
    kg_res <- 0.589; label("Growth rate constant of the resistant subpopulation (1/h)")       # Table 3 (ISOL MIC16): k_g,RES 0.589, %RSE 4.2

    log10bmax <- 9.94; label("Maximum bacterial density (log10 cfu/mL)") # Table 3 (ISOL MIC16): Bmax,HFIM 9.94, %RSE 1.0; the static time-kill value Bmax,STK was 8.86 (%RSE 0.9)

    # ---- Drug effect --------------------------------------------------------
    # Sigmoidal Emax on the susceptible and less susceptible subpopulations
    # (Equations 1-2); linear on the resistant subpopulation (Equation 3).
    emax_s  <- 2.38; label("Maximum drug effect rate constant, susceptible subpopulation (1/h)")             # Table 3 (ISOL MIC16): Emax,S 2.38, %RSE 4.3
    ec50_s  <- 8.31; label("Unbound concentration at half-maximal effect, susceptible subpopulation (mg/L)") # Table 3 (ISOL MIC16): EC50,S 8.31, %RSE 4.7
    hill_s  <- 2.45; label("Hill coefficient, susceptible subpopulation")                                 # Table 3 (ISOL MIC16): H,S 2.45, %RSE 9.0

    emax_ls <- 1.09; label("Maximum drug effect rate constant, less susceptible subpopulation (1/h)")             # Table 3 (ISOL MIC16): Emax,LS 1.09, %RSE 10.1
    ec50_ls <- 37.0; label("Unbound concentration at half-maximal effect, less susceptible subpopulation (mg/L)") # Table 3 (ISOL MIC16): EC50,LS 37.0, %RSE 43.5
    hill_ls <- fixed(1); label("Hill coefficient, less susceptible subpopulation")                                 # Table 3 (ISOL MIC16): H,LS 1, FIX

    klin_res <- 0.00940; label("Linear drug effect rate constant, resistant subpopulation (L/(mg*h))") # Table 3 (ISOL MIC16): k_lin,RES 0.00940, %RSE 8.8

    # ---- Initial bacterial population sizes ---------------------------------
    log10cfu0_s   <- 5.83; label("Fitted initial size of the susceptible subpopulation (log10 cfu/mL)")      # Table 3 (ISOL MIC16): cfu_t0,S 5.83, %RSE 1.3
    log10cfu0_ls  <- 0.837; label("Fitted initial size of the less susceptible subpopulation (log10 cfu/mL)") # Table 3 (ISOL MIC16): cfu_t0,LS 0.837, %RSE 31.2
    log10cfu0_res <- 0.528; label("Fitted initial size of the resistant subpopulation (log10 cfu/mL)")       # Table 3 (ISOL MIC16): cfu_t0,RES 0.528, %RSE 21.4

    # Methods, PK/PD simulations: "The initial size of the total bacterial
    # population was set to 10^6 cfu/mL. The initial size of each subpopulation
    # was scaled accordingly, based on the PK/PD model estimates."
    log10cfu0_total <- fixed(6); label("Total bacterial population at time zero (log10 cfu/mL)") # Methods, PK/PD simulations; NONMEM_code.txt THETA(1)

    # ---- Pharmacokinetic inter-individual variability -----------------------
    # Table 2 footnote a: %CV = sqrt(exp(omega^2) - 1) x 100, i.e.
    # omega^2 = log(CV^2 + 1). log(0.36^2 + 1) = 0.1218636 and
    # log(0.58^2 + 1) = 0.2899794, matching NONMEM_code.txt $OMEGA exactly.
    # Inherited from Laterre et al. rather than estimated here, hence fixed().
    etalcl ~ fixed(0.1218636) # Table 2: IIV on CL, 36 %CV
    etalvc ~ fixed(0.2899794) # Table 2: IIV on Vc, 58 %CV

    # ---- Residual unexplained variability -----------------------------------
    # Methods, PK/PD model: "RUV was described with an additive error model on
    # the log10 scale."
    addSd_CFUtotal <- 0.507; label("Additive residual SD on the log10 total bacterial count (log10 cfu/mL)")     # Table 3 (ISOL MIC16): RUV total 0.507, %RSE 4.6
    addSd_CFUres   <- 0.711; label("Additive residual SD on the log10 resistant bacterial count (log10 cfu/mL)") # Table 3 (ISOL MIC16): RUV RES 0.711, %RSE 6.8
  })

  model({
    # ---- 0. Guard against log10() of a non-positive solver excursion --------
    eps <- 1e-30

    # ---- 1. Individual pharmacokinetic parameters ---------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 2. Two-compartment plasma pharmacokinetics -------------------------
    d/dt(central)     <- -kel * central + k21 * peripheral1 - k12 * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Total and unbound plasma concentration (mg/L). The unbound concentration
    # drives all three drug effects (NONMEM_code.txt $DES: TEMO = A(1)/VC*FU).
    Cc <- central / vc
    Cu <- Cc * fu

    # ---- 3. Back-transformed bacterial parameters ---------------------------
    bmax <- 10^log10bmax

    # Rescale the fitted inocula so the total population starts at the
    # simulation inoculum while preserving the fitted subpopulation ratios.
    # For this strain the shift is
    # log10(1e6 / (10^5.83 + 10^0.837)) log10 cfu/mL.
    cfu0scale <- log10cfu0_total - log10(10^log10cfu0_s + 10^log10cfu0_ls)

    # ---- 4. Drug effect (Equations 1-3) -------------------------------------
    kdrug_s   <- emax_s  * Cu^hill_s  / (ec50_s^hill_s   + Cu^hill_s)
    kdrug_ls  <- emax_ls * Cu^hill_ls / (ec50_ls^hill_ls + Cu^hill_ls)
    kdrug_res <- klin_res * Cu

    # ---- 5. Bacterial subpopulations (Equations 1-3) ------------------------
    # Equation 3 caps the growth of the resistant subpopulation on the count of
    # the drug-free-agar population (S + LS), not on its own size, and the
    # resistant subpopulation does not contribute to that total count
    # (Methods, PK/PD model).
    cfu <- bact_susceptible + bact_less_susceptible
    cap <- 1 - cfu / bmax

    d/dt(bact_susceptible)      <- kg_s   * cap * bact_susceptible      - kdrug_s   * bact_susceptible
    d/dt(bact_less_susceptible) <- kg_ls  * cap * bact_less_susceptible - kdrug_ls  * bact_less_susceptible
    d/dt(bact_resistant)        <- kg_res * cap * bact_resistant        - kdrug_res * bact_resistant

    bact_susceptible(0)      <- 10^(log10cfu0_s   + cfu0scale)
    bact_less_susceptible(0) <- 10^(log10cfu0_ls  + cfu0scale)
    bact_resistant(0)        <- 10^(log10cfu0_res + cfu0scale)

    # ---- 6. Observations (log10 cfu/mL) -------------------------------------
    # NONMEM_code.txt $ERROR constrains the count read on temocillin-containing
    # agar to the count read on drug-free agar before it is reported.
    # max(., eps) floors the argument of log10() so that a small negative
    # solver excursion, which happens once a subpopulation is driven far below
    # one organism in the whole system, yields the floor value rather than NaN.
    resobs <- min(bact_resistant, cfu)

    CFUtotal <- log10(max(cfu, eps))
    CFUres   <- log10(max(resobs, eps))

    CFUtotal ~ add(addSd_CFUtotal)
    CFUres   ~ add(addSd_CFUres)
  })
}

