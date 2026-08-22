Abouelhassan_2024_sulbactam_human <- function() {
  description <- paste(
    "Three-compartment plasma-plus-epithelial-lining-fluid (ELF) population PK",
    "model for intravenous sulbactam in healthy adults, co-modelled from the",
    "mean plasma and ELF concentrations of 30 healthy volunteers given",
    "sulbactam 1 g IV q6h as a 3 h infusion (Abouelhassan 2024, using the",
    "Rodvold 2018 intrapulmonary data). A two-compartment plasma model was",
    "built first and selected over one compartment on AIC (44.1 versus 54.6),",
    "then the ELF concentrations were added as a third compartment with its own",
    "volume and its own asymmetric first-order exchange with plasma, which is",
    "what generates the sub-unity ELF penetration. The model was fitted by",
    "nonparametric adaptive grid (NPAG) with adaptive gamma in Pmetrics. Only",
    "mean concentrations were available, so no random effects were estimated",
    "and no covariates were assessed; for the 5000-subject Monte Carlo",
    "simulation underlying the paper's probability-of-target-attainment",
    "analysis the authors artificially inflated the dispersion of every",
    "parameter to 40% CV to mimic patient variability, and that assumption is",
    "carried here as fixed 40% CV lognormal between-subject variability on all",
    "seven parameters. Residual error was not reported and is fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Abouelhassan Y, Kuti JL, Nicolau DP, Abdelraouf K.",
    "Pharmacokinetic/pharmacodynamic analysis of sulbactam against",
    "Acinetobacter baumannii pneumonia: establishing in vivo efficacy targets",
    "in the epithelial lining fluid. JAC Antimicrob Resist. 2024;6(6):dlae203.",
    "doi:10.1093/jacamr/dlae203.",
    "Final parameter estimates from the Results section 'Population",
    "pharmacokinetics modelling in healthy volunteers'; Monte Carlo simulation",
    "settings from Methods 'Monte Carlo simulation and PTA estimation'.",
    "The plasma and ELF concentrations that were modelled are from",
    "Rodvold KA, Gotfried MH, Isaacs RD et al. Plasma and intrapulmonary",
    "concentrations of ETX2514 and sulbactam following intravenous",
    "administration of ETX2514SUL to healthy adult subjects. Antimicrob Agents",
    "Chemother. 2018;62(10):e01089-18. doi:10.1128/AAC.01089-18.",
    sep = " "
  )
  vignette <- "Abouelhassan_2024_sulbactam_elf_pneumonia"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(
      analyte = "sulbactam", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "sulbactam", units = "mg",
      specimen = "tissue", verified = TRUE
    ),
    elf = list(
      analyte = "sulbactam", units = "mg",
      specimen = "epithelial lining fluid", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30,
    n_studies      = 1,
    disease_state  = "healthy adult volunteers",
    dose_range     = "sulbactam 1 g IV q6h as a 3 h infusion",
    notes          = paste(
      "Methods, 'Population pharmacokinetics modelling in healthy volunteers':",
      "mean sulbactam plasma and ELF concentrations at 1, 2.5, 3.25, 4 and 6 h",
      "from 30 healthy volunteers (Rodvold 2018) were co-modelled. Because a",
      "single data set of mean concentrations from a homogeneous healthy",
      "volunteer population was used, no covariates were assessed and no",
      "individual-level variability was estimable. The Results section cites",
      "reference 23 for these data, which is a CLSI standards document; the",
      "Methods section cites reference 25 (Rodvold 2018), which is the actual",
      "source. See the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Results, 'Population pharmacokinetics modelling
    # in healthy volunteers': "clearance (CL), 13.69 L/h; volume of central
    # compartment (Vc), 2.49 L; intercompartmental transfer rate constants K12,
    # 10.65 h-1; K21, 21.80 h-1; K13, 7.47 h-1; K31, 1.42 h-1 and volume of the
    # ELF compartment (VELF), 2.44 L".
    #
    # WHICH PRINTED PAIR DRIVES THE ELF COMPARTMENT. The prose says the ELF
    # "was the third compartment", which reads as ELF <-> plasma being K13/K31.
    # That reading is arithmetically impossible and the K12/K21 pair is the ELF
    # pair. Three independent checks, all internal to the paper:
    #   (i) For a linear system the steady-state / AUC ratio is
    #       AUC(Celf)/AUC(Cc) = (k_in / k_out) * (Vc / VELF). Taking K13/K31 as
    #       the ELF pair gives (7.47 / 1.42) * (2.49 / 2.44) = 5.37, i.e. ELF
    #       concentrations 5.4-fold ABOVE plasma. The Discussion states the mean
    #       sulbactam ELF AUC to FREE plasma AUC is 0.81, and free can never
    #       exceed total, so the ELF-to-total-plasma ratio must be <= 0.81.
    #   (ii) Taking K12/K21 as the ELF pair gives
    #       (10.65 / 21.80) * (2.49 / 2.44) = 0.499, which with the reported
    #       0.81 ELF-to-free-plasma ratio implies an unbound plasma fraction of
    #       0.499 / 0.81 = 0.62 -- the accepted sulbactam value (~38% bound).
    #   (iii) Simulating 1 g q6h 0.5 h infusion with 40% CV on all parameters
    #       reproduces the Table 4 PTA column only under the K12/K21 reading
    #       (17% fT>MIC target: 100 / 100 / 100 / 99.7 / 98.1 / 93.1 / 79.4 /
    #       52.8 versus the published 100 / 100 / 99 / 98 / 96 / 92 / 81 / 62 at
    #       MIC 0.0625-8 mg/L). The K13/K31 reading returns 100 / 100 / 100 /
    #       100 / 100 / 100 / 100 / 100, and 83% versus a published 12% at
    #       MIC 16 for the 48% target.
    # Vc*(K12/K21 + K13/K31) is invariant to the swap, so steady-state volume
    # (16.8 L) does not discriminate; only the ELF ratio does. Each parameter
    # below therefore carries the paper's own printed label in its comment.
    lcl <- log(13.69)
    label("Plasma clearance CL (L/h)")
    # Results: CL = 13.69 L/h

    lvc <- log(2.49)
    label("Central (plasma) volume of distribution Vc (L)")
    # Results: Vc = 2.49 L

    lk12 <- log(7.47)
    label("Central-to-peripheral1 transfer rate constant, the paper's K13 (1/h)")
    # Results: K13 = 7.47 1/h -- the plasma peripheral pair (see the header note)

    lk21 <- log(1.42)
    label("Peripheral1-to-central transfer rate constant, the paper's K31 (1/h)")
    # Results: K31 = 1.42 1/h -- the plasma peripheral pair (see the header note)

    lk_central_elf <- log(10.65)
    label("Central-to-ELF transfer rate constant, the paper's K12 (1/h)")
    # Results: K12 = 10.65 1/h -- the ELF pair (see the header note)

    lk_elf_central <- log(21.80)
    label("ELF-to-central transfer rate constant, the paper's K21 (1/h)")
    # Results: K21 = 21.80 1/h -- the ELF pair (see the header note)

    lvelf <- log(2.44)
    label("Epithelial lining fluid compartment volume VELF (L)")
    # Results: VELF = 2.44 L

    # Between-subject variability. NOT estimated: the model was fitted to a
    # single set of MEAN plasma and ELF concentrations, so no random effects
    # were identifiable. Methods, 'Monte Carlo simulation and PTA estimation':
    # "The CV was artificially inflated to 40% for all parameters to be
    # consistent with CVs observed in patients", restated in the Discussion as
    # "the dispersion around the mean for all sulbactam parameter estimates was
    # inflated to 40%CV during simulation". Encoded as lognormal variances
    # fixed at log(1 + 0.40^2) = 0.148420, which reproduces a 40% CV exactly.
    # These are a simulation assumption, not an estimate -- see the vignette
    # Errata before re-fitting.
    etalcl            ~ fixed(0.148420)
    etalvc            ~ fixed(0.148420)
    etalk12           ~ fixed(0.148420)
    etalk21           ~ fixed(0.148420)
    etalk_central_elf ~ fixed(0.148420)
    etalk_elf_central ~ fixed(0.148420)
    etalvelf          ~ fixed(0.148420)

    # Residual error was not reported -- the model was fitted to mean
    # concentrations and the paper gives no error model. Fixed at zero per the
    # standing policy for unreported RUV; see the vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual error on plasma concentration (fraction; ZERO - not reported in source)")

    propSd_Celf <- fixed(0)
    label("Proportional residual error on ELF concentration (fraction; ZERO - not reported in source)")
  })

  model({
    # 1. Individual parameters
    cl            <- exp(lcl + etalcl)
    vc            <- exp(lvc + etalvc)
    k12           <- exp(lk12 + etalk12)
    k21           <- exp(lk21 + etalk21)
    k_central_elf <- exp(lk_central_elf + etalk_central_elf)
    k_elf_central <- exp(lk_elf_central + etalk_elf_central)
    velf          <- exp(lvelf + etalvelf)

    # 2. Micro-constants
    kel <- cl / vc

    # 3. ODE system. Dose is an intravenous infusion into `central`. The ELF
    #    exchange is asymmetric (k_central_elf < k_elf_central), which is what
    #    produces an ELF-to-plasma AUC ratio below one.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k_central_elf * central + k_elf_central * elf
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf)         <- k_central_elf * central - k_elf_central * elf

    # 4. Observations. ELF concentrations are unbound by construction, so
    #    %fT>MIC in ELF is computed directly from Celf without a protein
    #    binding correction.
    Cc   <- central / vc
    Celf <- elf / velf

    Cc   ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
