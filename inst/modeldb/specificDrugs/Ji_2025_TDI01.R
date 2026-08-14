Ji_2025_TDI01 <- function() {
  description <- paste0(
    "Population PK model for orally administered TDI01, a novel selective ROCK2 inhibitor developed ",
    "for acute lung injury / acute respiratory distress syndrome (ALI/ARDS), in healthy Chinese adults ",
    "(Ji 2025, final Pop-PK model). Disposition is one-compartment with first-order absorption plus a ",
    "gallbladder compartment that reproduces the double-peak hepato-enteral (enterohepatic) circulation ",
    "seen after single dosing: drug leaves the central compartment continuously into the gallbladder at ",
    "the fixed rate k2G, and the gallbladder discharges back into the absorption compartment at the fixed ",
    "rate kG1 only while the authors' FLAG time switch is open, which they set to the 8-10 h post-dose ",
    "window. All disposition parameters are apparent (per unit bioavailability, /F). No covariate model ",
    "was developed because the phase 1 population was demographically homogeneous. Note that the source ",
    "describes the final model as a one-compartment model 'with dosage effect', reporting that relative ",
    "bioavailability was higher at 200-400 mg than at 800-1,200 mg, but it tabulates no dose-effect ",
    "parameter and writes no dose term into its published equations; the packaged model therefore carries ",
    "F = 1 at every dose. See the vignette 'Assumptions and deviations' section."
  )
  reference <- paste(
    "Ji X, Cui Y.",
    "Model-informed drug development of novel ROCK2 inhibitor TDI01:",
    "population pharmacokinetic study and simulation.",
    "Front Pharmacol. 2025;16:1477607.",
    "doi:10.3389/fphar.2025.1477607.",
    sep = " "
  )
  vignette <- "Ji_2025_TDI01"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Dose enters `depot`. Declared explicitly because buildModelDb()'s dosing
  # heuristic only inspects whether "depot" / "central" are present.
  dosing <- "depot"

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "TDI01", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "TDI01", units = "mg", specimen = "plasma", verified = FALSE),
    gallbladder = list(analyte = "TDI01", units = "mg", specimen = "bile", verified = FALSE)
  )

  # No covariate model was developed. Discussion: "The covariate model was not
  # developed in this study for the following reasons: (1) The PK data was
  # obtained from healthy subjects, and their demographic and physiological
  # characteristics were similar, so it is difficult to select significant
  # covariates; (2) The purpose of this study was to evaluate the in vivo
  # exposure level of TDI01 ... through simulation." No covariate was screened
  # and rejected either, so there is nothing to record in
  # covariatesDataExcluded.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_observations = 776L,
    n_studies      = 1L,
    disease_state  = "Healthy volunteers.",
    dose_range     = paste(
      "Oral TDI01. Single-dose regimens: 400 mg, 800 mg, and 1,200 mg.",
      "Multiple-dose regimens: 200 mg QD for 7 days and 400 mg QD for 7 days."
    ),
    regions        = "China (Institute of Clinical Pharmacology, Peking University First Hospital, Beijing; sponsor Beijing Tide Pharmaceutical Co., Ltd.).",
    sampling_window = paste(
      "Single-dose groups: 0 (pre-dose), 1, 2, 3, 3.5, 4, 5, 6, 8, 10, 12, 24, 36, 48, 72 h post-dose.",
      "Multiple-dose groups: 0, 1, 2, 3, 3.5, 4, 5, 6, 8, 10, 12, 24, 48, 72, 96, 120, 144, 145, 146,",
      "147, 147.5, 148, 149, 150, 152, 154, 156, 168, 180, 192, 216 h."
    ),
    notes = paste(
      "The source reports NO baseline demographic table: age, body weight, sex, and body-mass-index",
      "distributions are not published, and the only characterisation given is that the 39 healthy",
      "subjects had a 'homogeneous demographic profile' (Discussion, Limitations). Those fields are",
      "therefore deliberately absent here rather than guessed. Data handling: 10 of 786 observations",
      "(1.3%) were below the limit of quantification; because that is under the paper's 10% threshold",
      "the BQL records were simply excluded, leaving 776 quantifiable concentrations. Observations with",
      "|CWRES| > 5 were treated as outliers and omitted during model building, then re-introduced for a",
      "sensitivity analysis. Ethics approval: People's Liberation Army General Hospital Medical Ethical",
      "Committee. Estimation was FOCEI in NONMEM 7.5, with diagnostics in R 3.5.3, Xpose and PsN 4.8.0."
    )
  )

  # Implementation notes (see the vignette 'Assumptions and deviations'
  # section for the full justification):
  # * Compartment indexing in the source equations. Ji 2025 Equations 3-5 use
  #   DADT(1) / DADT(2) / DADT(3) on the left-hand side but A / A(1) / A(2) on
  #   the right-hand side, and the accompanying legend states "A, A (1) and
  #   A (2) represent absorption, central and gallbladder compartments,
  #   respectively". Under that legend the three equations are internally
  #   consistent and map unambiguously onto depot / central / gallbladder:
  #     Eq 3  DADT(1) = FLAG*KG1*A(2) - Ka*A          -> d(depot)/dt
  #     Eq 4  DADT(2) = Ka*A - K2G*A(1) - K20*A(1)    -> d(central)/dt
  #     Eq 5  DADT(3) = K2G*A(1) - FLAG*KG1*A(2)      -> d(gallbladder)/dt
  #   Mass balance closes exactly (each gallbladder loss term is the matching
  #   depot gain term, and the central loss to bile is the gallbladder gain),
  #   so the DADT-vs-A index offset is a typographic artifact, not an
  #   ambiguity.
  # * Rate-constant parameterisation. The gallbladder is mechanism-defined (a
  #   delayed-release container), not a volume-defined peripheral compartment,
  #   so the canonical lq / lvp form does not apply. The canonical `kbm`
  #   (biliary excretion out of central) and `kehc` (gated enterohepatic
  #   release) names are used, matching Courlet_2023_cabamiquine whose source
  #   supplement uses the identical k2g / kg1 notation, and Jain_2011_sorafenib.
  # * FLAG time switch. The paper reports "the starting and ending time of
  #   hepato-enteral circulation was set as 8-10 h" and simulates 7-day QD and
  #   BID regimens in which the recirculation recurs. The switch is therefore
  #   encoded on time-after-dose, `tad(depot)`, so the 2 h release window
  #   reopens after every dose rather than firing once on absolute time. The
  #   `if (tdose > 0)` guard is required because tad() returns NA before the
  #   first dose record and an NA would propagate through the whole
  #   trajectory (see Courlet_2023_cabamiquine for the same guard). For the
  #   q12h regimen the window (8-10 h post-dose) closes before the next dose,
  #   so successive windows never overlap.
  # * IIV_k12. Table 1's third variance component is labelled "IIV_k12", a
  #   name that appears nowhere in Equations 3-5. It is the IIV on ka: k12 is
  #   the NONMEM micro-constant for transfer from compartment 1 (depot) to
  #   compartment 2 (central), i.e. ka; and ka is the only estimated
  #   structural parameter left once IIV_CL/F and IIV_V/F are accounted for
  #   (k2G and kG1 are both FIX and so cannot carry a variance).
  # * Unquantified dosage effect. The Results call the final model a
  #   "one-compartment model with dosage effect", and the Discussion reports
  #   that dose-normalised exposure implied higher bioavailability at
  #   200-400 mg than at 800-1,200 mg, "which need more clinical evidence to
  #   support this conclusion". No dose-effect parameter appears in Table 1
  #   and no dose term appears in Equations 3-5, so none is invented here;
  #   f(depot) is left at unity and every parameter stays apparent (/F). The
  #   vignette quantifies the consequence against the paper's own Table 3.
  ini({
    # Structural disposition parameters (Ji 2025 Table 1, "Final model"
    # column). All are apparent parameters (reported as CL/F, Vc/F, ka/F).
    lcl <- log(95);    label("Apparent clearance CL/F (L/h)")                                  # Table 1 CL/F = 95 L/h (RSE 10%; bootstrap median 95.15, 95% CI 82.05-109.52)
    lvc <- log(1400);  label("Apparent central volume of distribution Vc/F (L)")               # Table 1 Vc/F = 1,400 L (RSE 8%; bootstrap median 1,395.22, 95% CI 1,182.17-1,586.25)
    lka <- log(0.345); label("First-order absorption rate constant ka/F (1/h)")                # Table 1 ka/F = 0.345 1/h (RSE 10%; bootstrap median 0.345, 95% CI 0.264-0.424)

    # Hepato-enteral (enterohepatic) circulation rate constants. Both were
    # held FIX by the authors, so both carry fixed().
    lkbm  <- fixed(log(0.023)); label("Central -> gallbladder biliary transfer rate constant k2G/F (1/h)")            # Table 1 k2G/F = 0.023 (FIX)
    lkehc <- fixed(log(2));     label("Gallbladder -> depot release rate constant kG1/F (1/h), open only while FLAG = 1") # Table 1 kG1/F = 2 (FIX)

    # FLAG time switch bounds. Results: "the starting and ending time of
    # hepato-enteral circulation was set as 8-10 h". Set by the authors rather
    # than estimated (they carry no estimate, RSE, or bootstrap CI in Table 1),
    # so both are fixed().
    ltgb <- fixed(log(8)); label("Hepato-enteral circulation start time, FLAG onset (h post-dose)")   # Results: hepato-enteral circulation "set as 8-10 h" -> onset 8 h
    ldgb <- fixed(log(2)); label("Hepato-enteral circulation window duration, FLAG open span (h)")    # Results: hepato-enteral circulation "set as 8-10 h" -> duration 10 - 8 = 2 h

    # Inter-individual variability. Methods Equation 1 is an exponential IIV
    # model (Pi = theta * exp(eta_i)), so the Table 1 percentages are CVs of a
    # log-normal and convert to log-scale variance via omega^2 = log(1 + CV^2):
    #   CL/F  CV 46.8% -> log(1 + 0.468^2) = 0.198051
    #   Vc/F  CV 32.2% -> log(1 + 0.322^2) = 0.098654
    #   ka    CV 47.3% -> log(1 + 0.473^2) = 0.201903
    etalcl ~ 0.198051   # Table 1 IIV_CL/F = 46.8% CV (RSE 12%, shrinkage 1%;  bootstrap median 45.9%, 95% CI 0.354-0.558)
    etalvc ~ 0.098654   # Table 1 IIV_V/F  = 32.2% CV (RSE 18%, shrinkage 10%; bootstrap median 31.1%, 95% CI 0.199-0.411)
    etalka ~ 0.201903   # Table 1 IIV_k12  = 47.3% CV (RSE 16%, shrinkage 19%; bootstrap median 45.9%, 95% CI 0.209-0.660) -- k12 = depot-to-central rate = ka

    # Residual error. Methods Equation 2 is Yij = Cij * (1 + eps1_ij) + eps2_ij,
    # i.e. eps1 multiplies the prediction (proportional) and eps2 is additive.
    # The surrounding prose labels eps1 "additive" and eps2 "proportional",
    # which contradicts its own equation; the equation governs, and the Table 1
    # magnitudes confirm it (34.2% is a fraction, 1.77 is on the ng/mL scale).
    propSd <- 0.342; label("Proportional residual error (fraction)")  # Table 1 Prop.error = 34.2% (RSE 3%;  bootstrap median 34.3%, 95% CI 0.322-0.363)
    addSd  <- 1.77;  label("Additive residual error (ng/mL)")         # Table 1 Add.error  = 1.77   (RSE 47%; bootstrap median 1.795, 95% CI 1.362-2.271)
  })
  model({
    # Individual parameters (exponential IIV, Methods Equation 1).
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    kbm  <- exp(lkbm)
    kehc <- exp(lkehc)
    tgb  <- exp(ltgb)
    dgb  <- exp(ldgb)

    # K20 = CL/Vc, the elimination rate constant of the central compartment
    # (Ji 2025, line preceding Equation 3).
    kel <- cl / vc

    # FLAG: the hepato-enteral circulation time switch, open on the
    # [8 h, 10 h] window after each dose. tad() returns NA until the first
    # dose record, so the branch guard keeps that NA out of the ODEs (a
    # comparison against NaN is false in C, so the else branch is taken).
    tdose <- tad(depot)
    if (tdose > 0) {
      flag <- (tdose >= tgb) * (tdose <= tgb + dgb)
    } else {
      flag <- 0
    }

    # ODE system, Ji 2025 Equations 3-5 (see the compartment-indexing note
    # above for the mapping of the source's A / A(1) / A(2) onto these states).
    d/dt(depot)       <- -ka * depot + flag * kehc * gallbladder
    d/dt(central)     <-  ka * depot - kel * central - kbm * central
    d/dt(gallbladder) <-  kbm * central - flag * kehc * gallbladder

    # TDI01 plasma concentration. `central` is in mg (dose units) and vc is in
    # L, so central/vc is mg/L; 1 mg/L = 1,000 ng/mL, hence the factor 1000 to
    # reach the source's ng/mL bioanalytical scale.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
