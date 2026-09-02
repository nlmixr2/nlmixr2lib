Gao_2012_exenatide_rat <- function() {
  description <- paste(
    "Preclinical (rat). Target-mediated drug disposition (TMDD) model for",
    "exendin-4 (exenatide) in male Sprague-Dawley rats (Gao 2012). Free drug",
    "in the central compartment binds glucagon-like peptide 1 receptor",
    "(GLP-1R) with second-order rate constant kon to form a drug-receptor",
    "complex that either dissociates (koff) or is internalised and degraded",
    "(kint); free drug additionally distributes to a nonspecific tissue",
    "compartment (k12 / k21) and is eliminated linearly from plasma (kel,",
    "physiologically renal filtration: CLc = kel * Vc = 3.62 mL/min, close to",
    "the reported renal clearance of 3.44 mL/min). The total receptor pool",
    "Rtot is held constant with no synthesis or degradation. Receptor binding",
    "is the source of the observed dose-dependent (nonlinear) clearance.",
    "Subcutaneous input is a first-order depot; the absorption rate constant",
    "was estimated separately for each of the three subcutaneous dose groups",
    "and DECREASES with dose (0.00820, 0.00579 and 0.00273 1/min at 0.5, 5",
    "and 50 nmol) -- the value carried in ini() is the 0.5 nmol estimate and",
    "must be overridden to simulate the higher dose groups. Naive-pooled fit",
    "to mean profiles in ADAPT II; there is no between-subject variability",
    "and the residual-variance parameters were not reported."
  )
  reference <- paste(
    "Gao W, Jusko WJ (2012).",
    "Target-mediated Pharmacokinetic and Pharmacodynamic Model of Exendin-4",
    "in Rats, Monkeys, and Humans.",
    "Drug Metabolism and Disposition 40(5):990-997.",
    "doi:10.1124/dmd.111.042291.",
    "PMCID PMC3336795.",
    sep = " "
  )
  vignette <- "Gao_2012_exenatide"

  units <- list(
    time          = "min",
    dosing        = "nmol",
    concentration = "nmol/L (nM); the paper converted all doses to moles using a molecular weight of 4186.6 g/mol for exendin-4"
  )

  compartmentData <- list(
    depot       = list(analyte = "exenatide", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "exenatide", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "exenatide", units = "nmol", specimen = "tissue", verified = TRUE),
    complex     = list(analyte = "exenatide/GLP-1R complex", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species      = "rat (male Sprague-Dawley)",
    n_subjects   = 7L,
    n_studies    = 1L,
    weight_range = "350-370 g",
    disease_state = "Healthy (normoglycaemic) male Sprague-Dawley rats.",
    dose_range   = paste(
      "Single intravenous bolus and single subcutaneous bolus at 0.5, 5 and",
      "50 nmol, plus continuous intravenous infusion at 0.5, 5 and 50 nmol/h.",
      "n = 4-7 per route and dose group."
    ),
    notes = paste(
      "Gao 2012 Materials and Methods, 'In the rat PK study'. Concentration",
      "data were supplied by Amylin Pharmaceuticals and assayed by a two-site",
      "sandwich assay with a minimum detectable concentration of 15 pM.",
      "Model fitted to MEAN profiles by naive pooling in ADAPT II with the",
      "maximum-likelihood method, so n_subjects records the largest per-group",
      "n rather than a pooled population size."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural TMDD parameters -- Gao 2012 Table 2 (rat), CV% in
    # parentheses in the source. All rates are 1/min; the rat model was
    # run on the nanomolar concentration scale (kon in 1/(nM*min) and
    # Rtot in nmol/L), which is confirmed by koff/kon = 0.0153/0.0207 =
    # 0.74 nM, the KD quoted in Gao 2012 Results 'Rat PK'.
    # ------------------------------------------------------------------
    lkel  <- log(0.0839)  ; label("Linear elimination rate constant from the central compartment (kel, 1/min)")     # Table 2: kel = 0.0839 (CV 10%)
    lk12  <- log(0.0282)  ; label("Transfer rate constant central -> peripheral1 (kpt, 1/min)")                     # Table 2: kpt = 0.0282 (CV 15%)
    lk21  <- log(0.0213)  ; label("Transfer rate constant peripheral1 -> central (ktp, 1/min)")                     # Table 2: ktp = 0.0213 (CV 5%)
    lvc   <- log(0.0432)  ; label("Central volume of distribution (Vc, L)")                                         # Table 2: Vc = 43.2 mL = 0.0432 L (CV 12%); CLc = kel * Vc = 3.62 mL/min per Results 'Rat PK'
    lkon  <- log(0.0207)  ; label("Second-order association rate constant of exendin-4 with GLP-1R (kon, 1/(nM*min))") # Table 2: kon = 0.0207 (CV 42%)
    lkoff <- log(0.0153)  ; label("First-order dissociation rate constant of the drug-receptor complex (koff, 1/min)") # Table 2: koff = 0.0153 (CV 206%)
    lkint <- log(0.0966)  ; label("Internalisation / degradation rate constant of the drug-receptor complex (kint, 1/min)") # Table 2: kint = 0.0966 (CV 38%)
    lrtot <- log(5.21)    ; label("Total GLP-1R concentration, held constant (Rtot, nmol/L)")                       # Table 2: Rtot = 5.21 nmol/L (CV 5%)

    # ------------------------------------------------------------------
    # Subcutaneous absorption. Gao 2012 estimated a SEPARATE ka for each
    # subcutaneous dose group and reported that it decreases with dose
    # (Table 2 and Results 'Rat PK': "The absorption rate constant in rats
    # decreased with increasing dose"); Fig. 7A plots the same three
    # values against dose. There is no published functional form linking
    # ka to dose, so a single ka is carried here and the value for the
    # dose group being simulated must be substituted:
    #   0.5 nmol  ka1 = 0.00820 1/min  (CV  9%)   <- the value below
    #   5  nmol   ka2 = 0.00579 1/min  (CV 11%)
    #   50 nmol   ka3 = 0.00273 1/min  (CV 11%)
    # ------------------------------------------------------------------
    lka     <- log(0.00820)     ; label("First-order subcutaneous absorption rate constant at the 0.5 nmol dose (ka1, 1/min)") # Table 2: ka1 = 0.00820 (CV 9%)
    lfdepot <- fixed(log(1))    ; label("Absolute subcutaneous bioavailability (F, fraction)")                      # Table 2: F = 1 (fixed); Results 'Rat PK': "Bioavailability was estimated as close to 1 and then fixed as 1 in the final model"

    # ------------------------------------------------------------------
    # Residual error. Gao 2012 Methods states the variance model
    # Vi = (sigma1 + sigma2 * Y)^2 -- i.e. a combined additive plus
    # proportional standard deviation -- but does NOT report sigma1 or
    # sigma2 anywhere in the paper. Both terms are therefore encoded at
    # zero rather than invented. See vignette Assumptions and deviations.
    # ------------------------------------------------------------------
    addSd  <- fixed(0) ; label("Additive residual SD (nmol/L; ZERO - sigma1 not reported in the source)")     # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma1 not reported
    propSd <- fixed(0) ; label("Proportional residual SD (fraction; ZERO - sigma2 not reported in the source)") # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma2 not reported

    # No inter-individual variability: the analysis is a naive-pooled fit
    # to mean concentration-time profiles, so no OMEGA block exists.
  })

  model({
    # 1. Individual parameters (no IIV; naive-pooled mean-data fit).
    kel  <- exp(lkel)
    k12  <- exp(lk12)
    k21  <- exp(lk21)
    vc   <- exp(lvc)
    kon  <- exp(lkon)
    koff <- exp(lkoff)
    kint <- exp(lkint)
    rtot <- exp(lrtot)
    ka   <- exp(lka)

    # 2. Free (unoccupied) receptor concentration. Gao 2012 eq. 3 writes
    # the binding term as kon * (Rtot - RC) * C, i.e. the receptor pool is
    # constant and free receptor is the difference between the total pool
    # and the occupied fraction. RC is a CONCENTRATION in the source; the
    # `complex` state below is the corresponding AMOUNT, so it is divided
    # by vc here.
    rfree <- rtot - complex / vc

    # 3. ODE system -- Gao 2012 eqs. 1-4, rewritten from the published
    # concentration form (dC/dt) to the amount form used by rxode2 by
    # multiplying eq. 1 through by Vc. The subcutaneous input function of
    # eq. 4, input(t) = ka * F * Dose * exp(-ka * t) / Vc, is exactly a
    # first-order depot with bioavailability F, so it is implemented as a
    # depot compartment rather than as an explicit exponential term.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (kel + k12) * central + k21 * peripheral1 -
                          kon * rfree * central + koff * complex
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(complex)     <-  kon * rfree * central - (koff + kint) * complex

    # 4. Bioavailability on the subcutaneous depot (F in eq. 4).
    f(depot) <- exp(lfdepot)

    # 5. Observation: free exendin-4 concentration in plasma (nmol/L).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
