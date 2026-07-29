Anderson_1998_paracetamol <- function() {
  description <- paste(
    "One-compartment oral PK model for paracetamol (acetaminophen) with an",
    "explicit cerebrospinal-fluid (CSF) equilibration compartment in nine",
    "ventilator-dependent children (5 months to 12 years) with indwelling",
    "ventricular drains for raised intracranial pressure (Anderson 1998",
    "NONMEM fit, Table 3). First-order absorption, single nasogastric dose",
    "of 40 mg/kg paracetamol elixir, plasma + CSF sampled hourly for 4 h",
    "and 2-hourly through 10 h. CSF concentration follows the plasma",
    "concentration with first-order equilibration rate keq = ln(2)/teq and",
    "steady-state ratio PC = Ccsf/Cc. Parameters are standardized to a",
    "70 kg adult using fixed allometric exponents (0.75 on CL, 1 on V,",
    "0.25 on the equilibration half-time teq; keq therefore scales with",
    "exponent -0.25). The published equation 2 for residual error",
    "var = SF^2 * (C^PWR + V) is unconventional and the NONMEM PWR and V",
    "terms are not reported; placeholder additive residual SDs are used so",
    "the model simulates plausibly (see vignette Errata)."
  )
  reference <- paste(
    "Anderson BJ, Holford NHG, Woollard GA, Chan PLS (1998).",
    "Paracetamol plasma and cerebrospinal fluid pharmacokinetics in children.",
    "British Journal of Clinical Pharmacology 46(3):237-243.",
    "doi:10.1046/j.1365-2125.1998.00780.x.",
    sep = " "
  )
  vignette <- "Anderson_1998_paracetamol"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate; reference weight 70 kg per",
        "Anderson 1998 Methods (Pharmacokinetic modelling, equation 'Pi =",
        "Pstd * (Wi/Wstd)^b'). Children studied weighed 8-50 kg (Table 1)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_range      = "5 months to 12 years",
    age_median     = "5 years (Table 1)",
    weight_range   = "8 to 50 kg",
    weight_median  = "20 kg (Table 1)",
    sex_female_pct = 100 * 3 / 9,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Ventilator-dependent paediatric intensive-care patients with",
      "external ventricular drains placed for raised intracranial",
      "pressure (>20 mmHg). Diagnoses: closed head injury (n=3),",
      "subdural / extradural / compressed-skull haematoma (n=3),",
      "brainstem encephalitis (n=1), Dandy-Walker cyst (n=1),",
      "posterior fossa tumour (n=1). Six of nine children had suffered",
      "trauma (Anderson 1998 Methods, Table 1)."
    ),
    dose_range     = paste(
      "Single nasogastric dose of paracetamol elixir 40 mg/kg",
      "(Wellcome New Zealand, 250 mg / 5 mL); approximate range",
      "320-2000 mg given the 8-50 kg weight range."
    ),
    sampling       = paste(
      "Arterial blood and CSF (from the external ventricular drain)",
      "sampled hourly for the first 4 h and 2-hourly through 10 h",
      "after the dose (Anderson 1998 Methods)."
    ),
    regions        = "New Zealand (Auckland Children's Hospital)",
    notes          = paste(
      "Patient 1 was studied on two separate occasions five days apart;",
      "the two occasions are reported as records 1a and 1b in Tables 2",
      "and 3 of Anderson 1998. NONMEM final parameter estimates use the",
      "FOCE method with subroutine ADVAN6. Children with CSF red blood",
      "cell counts above 4 x 10^10 / L (a 100-fold dilution of the",
      "normal RBC count) were excluded. CSF drainage rates were",
      "3-10 mL/h, which Anderson 1998 Discussion notes is minor relative",
      "to the ~10 L/h plasma clearance and so not expected to alter the",
      "CSF kinetics."
    )
  )

  ini({
    # Structural parameters: Anderson 1998 Table 3 'Geometric mean' row,
    # standardized to a 70 kg adult; allometric exponents below apply
    # individual-level scaling. The NONMEM analysis used FOCE with
    # subroutine ADVAN6 (Methods, Population parameter estimations (b)).

    lka  <- log(0.77)
    label("Absorption rate constant ka (1/h)")  # Table 3 geometric mean: ka = 0.77 1/h (CV 49%)
    lcl  <- log(10.2)
    label("Apparent oral clearance CL/F (L/h)")  # Table 3 geometric mean: CL = 10.2 L/h (CV 47%); CL/F per Methods
    lvc  <- log(67.1)
    label("Apparent central volume of distribution V/F (L)")  # Table 3 geometric mean: V = 67.1 L (CV 58%); V/F per Methods
    lkeq <- log(log(2) / 0.72)
    label("Plasma-to-CSF equilibration rate constant keq (1/h)")  # Table 3 geometric mean: teq = 0.72 h (CV 117%); keq = ln(2)/teq = 0.963 1/h
    lpc  <- log(1.18)
    label("CSF / plasma partition coefficient PC (unitless)")  # Table 3 geometric mean: PC = 1.18 (CV 8%)

    # Allometric exponents: fixed per Anderson 1998 Methods,
    # "The allometric exponential (b) was assumed to be 0.75 for
    # clearance, 1 for volume of distribution and 0.25 for equilibration
    # half-time (teq)." keq = ln(2)/teq scales with exponent -0.25.
    allo_cl  <- fixed(0.75)
    label("Allometric exponent on CL (unitless, fixed)")  # Methods: fixed at 0.75
    allo_vc  <- fixed(1)
    label("Allometric exponent on V (unitless, fixed)")  # Methods: fixed at 1
    allo_keq <- fixed(-0.25)
    label("Allometric exponent on keq (unitless, fixed; equals -0.25 because teq scales at 0.25)")  # Methods: teq fixed at 0.25 -> keq at -0.25

    # IIV: paper states "interindividual variability in model parameters
    # was modelled by an exponential variance model" (Methods, NONMEM
    # paragraph), i.e. log-normal random effects. Table 3 reports %CV per
    # parameter; the internal variance is omega^2 = log(1 + CV^2).
    # Anderson 1998 also notes "A parameter covariance matrix was
    # incorporated into the structural model"; the off-diagonal covariance
    # terms are not reported in the paper, so this model uses a
    # diagonal-only IIV (see vignette Errata).
    etalka  ~ 0.2153  # CV(ka)  = 49 %  -> omega^2 = log(1 + 0.49^2) = 0.2153 (Table 3)
    etalcl  ~ 0.1996  # CV(CL)  = 47 %  -> omega^2 = log(1 + 0.47^2) = 0.1996 (Table 3)
    etalvc  ~ 0.2900  # CV(V)   = 58 %  -> omega^2 = log(1 + 0.58^2) = 0.2900 (Table 3)
    etalkeq ~ 0.8623  # CV(teq) = 117 % -> omega^2 = log(1 + 1.17^2) = 0.8623 (Table 3); same %CV applies to keq = ln(2)/teq because the transform is multiplicative
    etalpc  ~ 0.006379 # CV(PC) = 8 %   -> omega^2 = log(1 + 0.08^2) = 0.006379 (Table 3)

    # Residual error. Anderson 1998 Methods (NONMEM paragraph) states:
    # "An additive term characterized the residual error." Table 3 reports
    # SF_C = 2.9 mmol/L (plasma) and SF_CCSF = 2.1 mmol/L (CSF), but the
    # paper's residual-error equation 2 is the MKMODEL form
    # var = SF^2 * (C^PWR + V), and PWR / V are NOT reported for the
    # NONMEM fit. The reported magnitudes (2.9 and 2.1 mmol/L) are
    # 20-50x larger than the observed paracetamol plasma concentrations
    # (~0.05-0.15 mmol/L per Figures 2-3) which is incompatible with
    # a standard NONMEM additive residual SD on the concentration scale.
    # Following the precedent of Park 2001 ketoprofen (which used a
    # placeholder when residual error was not interpretable from the
    # source), this model uses a small additive residual SD on each
    # output so that simulations produce plausible noise; the paper's
    # reported SF values are recorded in the comments and in the vignette
    # Errata for traceability. Units are mg/L (matching units$concentration).
    addSd      <- 1.5
    label("Additive residual error on plasma Cc (mg/L) - placeholder; see vignette Errata")   # paper reports SF_C = 2.9 mmol/L = 438 mg/L (Table 3, NONMEM); placeholder 1.5 mg/L (~ 10 percent of peak)
    addSd_Ccsf <- 1.5
    label("Additive residual error on CSF Ccsf (mg/L) - placeholder; see vignette Errata")    # paper reports SF_CCSF = 2.1 mmol/L = 317 mg/L (Table 3, NONMEM); placeholder 1.5 mg/L
  })

  model({
    # Individual PK parameters, allometrically scaled to each subject's
    # body weight against the 70 kg reference. Anderson 1998 Methods.
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^allo_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^allo_vc
    keq <- exp(lkeq + etalkeq) * (WT / 70)^allo_keq
    pc  <- exp(lpc  + etalpc)

    # Micro-constant
    kel <- cl / vc

    # ODE system (Anderson 1998 Methods, equation 1):
    #   dAgut / dt = -ka * Agut
    #   dC    / dt = (Agut * ka - C * CL) / V
    #   dCcsf / dt = ln(2)/teq * (C * PC - Ccsf)
    # Here `central` is the amount in the central compartment (mg) and
    # `brain_csf` is the CSF concentration (mg/L); the latter is an
    # effect-compartment-style state whose time derivative is in
    # concentration units.
    d/dt(depot)     <- -ka * depot
    d/dt(central)   <-  ka * depot - kel * central
    d/dt(brain_csf) <-  keq * (pc * (central / vc) - brain_csf)

    # Bioavailability: paracetamol was administered as nasogastric elixir
    # and the NONMEM analysis identifies CL/F and V/F (Methods), so F is
    # absorbed into the apparent CL and V parameters and is left at the
    # rxode2 default of 1 for the depot.

    Cc   <- central / vc
    Ccsf <- brain_csf
    Cc   ~ add(addSd)
    Ccsf ~ add(addSd_Ccsf)
  })
}
