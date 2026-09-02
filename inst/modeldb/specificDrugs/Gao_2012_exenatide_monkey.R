Gao_2012_exenatide_monkey <- function() {
  description <- paste(
    "Preclinical (rhesus monkey). Target-mediated drug disposition (TMDD)",
    "model for exendin-4 (exenatide) in male rhesus monkeys (Gao 2012). Same",
    "structure as the companion rat and human models: free drug in the",
    "central compartment binds glucagon-like peptide 1 receptor (GLP-1R) with",
    "second-order rate constant kon to form a drug-receptor complex that",
    "either dissociates (koff) or is internalised and degraded (kint); free",
    "drug additionally distributes to a nonspecific tissue compartment",
    "(k12 / k21) and is eliminated linearly from plasma (kel). The resulting",
    "CLc = kel * Vc = 2.39 mL/min/kg falls inside the reported glomerular",
    "filtration rate of healthy monkeys (2.2-3.6 mL/min/kg). The central",
    "volume is reported per kilogram, so body weight (WT) is a required",
    "covariate. Unlike the rat and human fits, subcutaneous bioavailability",
    "was ESTIMATED (0.688) rather than fixed at 1, and the binding parameters",
    "kon, koff and kint are very poorly identified (CV 179-2990%) because the",
    "monkey concentrations were digitised from a single published figure;",
    "the implied KD of 0.12 pM is roughly 5000-fold lower than in rats and",
    "should not be treated as a reliable affinity estimate. The absorption",
    "rate constant was estimated separately per subcutaneous dose group and",
    "decreases with dose. Naive-pooled fit to mean profiles in ADAPT II; no",
    "between-subject variability and no reported residual-variance values."
  )
  reference <- paste(
    "Gao W, Jusko WJ (2012).",
    "Target-mediated Pharmacokinetic and Pharmacodynamic Model of Exendin-4",
    "in Rats, Monkeys, and Humans.",
    "Drug Metabolism and Disposition 40(5):990-997.",
    "doi:10.1124/dmd.111.042291.",
    "PMCID PMC3336795.",
    "Monkey concentration-time data were digitised by Gao and Jusko from",
    "Ai G, Chen Z, Shan C, Che J, Hou Y, Cheng Y (2008),",
    "Eur J Pharmacol 588(2-3):297-303.",
    sep = " "
  )
  vignette <- "Gao_2012_exenatide"

  units <- list(
    time          = "min",
    dosing        = "pmol",
    concentration = "pmol/L (pM); the paper converted all doses to moles using a molecular weight of 4186.6 g/mol for exendin-4"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight. Required because Gao 2012 Table 3 reports the monkey central volume of distribution per kilogram (Vc = 69.3 mL/kg); the model multiplies it by WT to obtain an absolute volume in litres.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 4.3 kg, the mean body weight of the three male rhesus monkeys in Ai 2008 as reported in Gao 2012 Materials and Methods (4.3 +/- 0.7 kg). Time-fixed. Volume scales linearly with WT (exponent 1), which is how the source reports it, not as an estimated allometric exponent.",
      source_name        = "b.wt."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "exenatide", units = "pmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "exenatide", units = "pmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "exenatide", units = "pmol", specimen = "tissue", verified = TRUE),
    complex     = list(analyte = "exenatide/GLP-1R complex", units = "pmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species       = "rhesus monkey (male)",
    n_subjects    = 3L,
    n_studies     = 1L,
    weight_range  = "4.3 +/- 0.7 kg (mean +/- SD)",
    sex_female_pct = 0,
    disease_state = "Healthy male rhesus monkeys.",
    dose_range    = "Single subcutaneous injection of 1, 3 or 10 ug/kg, or a single intravenous injection of 3 ug/kg. n = 3.",
    notes = paste(
      "Gao 2012 Materials and Methods, 'In the monkey PK study'. Serum",
      "exendin-4 was measured by radioimmunoassay with a linear range of",
      "25-2000 pg/mL and a limit of quantitation of 25 pg/mL. The mean",
      "concentrations analysed by Gao and Jusko were captured by computer",
      "digitisation from Ai 2008, so this is a re-analysis of published",
      "figure data rather than of the original observations. Model fitted to",
      "MEAN profiles by naive pooling in ADAPT II."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural TMDD parameters -- Gao 2012 Table 3, 'Monkeys' column
    # (CV% in parentheses in the source). All rates are 1/min. The monkey
    # model was run on the PICOMOLAR concentration scale (kon in
    # 1/(pM*min) and Rtot in pM), confirmed by koff/kon = 0.0326/0.272 =
    # 0.12 pM, the KD quoted in Gao 2012 Results 'Monkey PK'.
    # ------------------------------------------------------------------
    lkel  <- log(0.0346)   ; label("Linear elimination rate constant from the central compartment (kel, 1/min)")      # Table 3 Monkeys: kel = 0.0346 (CV 10%)
    lk12  <- log(0.0143)   ; label("Transfer rate constant central -> peripheral1 (kpt, 1/min)")                      # Table 3 Monkeys: kpt = 0.0143 (CV 17%)
    lk21  <- log(0.00593)  ; label("Transfer rate constant peripheral1 -> central (ktp, 1/min)")                      # Table 3 Monkeys: ktp = 0.00593 (CV 17%)
    lvc   <- log(0.0693)   ; label("Central volume of distribution per kilogram (Vc, L/kg)")                          # Table 3 Monkeys: Vc = 69.3 mL/kg = 0.0693 L/kg (CV 12%); CLc = kel * Vc = 2.39 mL/min/kg per Results 'Monkey PK'
    lkon  <- log(0.272)    ; label("Second-order association rate constant of exendin-4 with GLP-1R (kon, 1/(pM*min))") # Table 3 Monkeys: kon = 0.272 (CV 607%)
    lkoff <- log(0.0326)   ; label("First-order dissociation rate constant of the drug-receptor complex (koff, 1/min)") # Table 3 Monkeys: koff = 0.0326 (CV 2990%)
    lkint <- log(0.00211)  ; label("Internalisation / degradation rate constant of the drug-receptor complex (kint, 1/min)") # Table 3 Monkeys: kint = 0.00211 (CV 179%)
    lrtot <- log(60.6)     ; label("Total GLP-1R concentration, held constant (Rtot, pmol/L)")                        # Table 3 Monkeys: Rtot = 60.6 pM (CV 33%)

    # ------------------------------------------------------------------
    # Subcutaneous absorption. Gao 2012 Table 3 reports only the RANGE of
    # the per-dose-group estimates for monkeys, "ka = 0.0244-0.0142 1/min
    # (CV 10-14%)", with Results 'Monkey PK' stating that ka "decreased
    # with dose, as observed in rats". Figure 7A plots the individual
    # points against dose; reading them off that figure gives
    #   1  ug/kg   ka = 0.0244 1/min   (the upper end of the Table 3 range;
    #                                   the value below)
    #   3  ug/kg   ka = 0.0205 1/min   (DIGITISED from Fig. 7A, +/- ~5%;
    #                                   this value is not printed anywhere
    #                                   in the text or tables)
    #   10 ug/kg   ka = 0.0142 1/min   (the lower end of the Table 3 range)
    # Substitute the value for the dose group being simulated.
    # ------------------------------------------------------------------
    lka     <- log(0.0244) ; label("First-order subcutaneous absorption rate constant at the 1 ug/kg dose (ka, 1/min)") # Table 3 Monkeys: ka range 0.0244-0.0142 (CV 10-14%); upper end assigned to the lowest dose per Fig. 7A
    lfdepot <- log(0.688)  ; label("Absolute subcutaneous bioavailability (F, fraction)")                             # Table 3 Monkeys: F = 0.688 (CV 8%) -- ESTIMATED, unlike the rat and human fits where F was fixed at 1

    # ------------------------------------------------------------------
    # Residual error. Gao 2012 Methods states the variance model
    # Vi = (sigma1 + sigma2 * Y)^2 -- a combined additive plus
    # proportional standard deviation -- but neither sigma1 nor sigma2 is
    # reported. Both are encoded at zero rather than invented.
    # ------------------------------------------------------------------
    addSd  <- fixed(0) ; label("Additive residual SD (pmol/L; ZERO - sigma1 not reported in the source)")     # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma1 not reported
    propSd <- fixed(0) ; label("Proportional residual SD (fraction; ZERO - sigma2 not reported in the source)") # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma2 not reported

    # No inter-individual variability: naive-pooled fit to mean profiles.
  })

  model({
    # 1. Individual parameters (no IIV; naive-pooled mean-data fit).
    # Vc is reported per kilogram in Table 3, so it is scaled linearly by
    # body weight to give an absolute volume in litres. This is the
    # source's own normalisation, not an estimated allometric exponent.
    kel  <- exp(lkel)
    k12  <- exp(lk12)
    k21  <- exp(lk21)
    vc   <- exp(lvc) * WT
    kon  <- exp(lkon)
    koff <- exp(lkoff)
    kint <- exp(lkint)
    rtot <- exp(lrtot)
    ka   <- exp(lka)

    # 2. Free (unoccupied) receptor concentration -- Gao 2012 eq. 3 binding
    # term kon * (Rtot - RC) * C. RC is a CONCENTRATION in the source; the
    # `complex` state is the corresponding AMOUNT, hence the division by vc.
    rfree <- rtot - complex / vc

    # 3. ODE system -- Gao 2012 eqs. 1-4, rewritten from the published
    # concentration form to the amount form used by rxode2 by multiplying
    # eq. 1 through by Vc. The subcutaneous input of eq. 4,
    # input(t) = ka * F * Dose * exp(-ka * t) / Vc, is exactly a
    # first-order depot with bioavailability F.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (kel + k12) * central + k21 * peripheral1 -
                          kon * rfree * central + koff * complex
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(complex)     <-  kon * rfree * central - (koff + kint) * complex

    # 4. Bioavailability on the subcutaneous depot (F in eq. 4).
    f(depot) <- exp(lfdepot)

    # 5. Observation: free exendin-4 concentration in serum (pmol/L).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
