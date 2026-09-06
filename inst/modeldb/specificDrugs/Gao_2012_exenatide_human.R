Gao_2012_exenatide_human <- function() {
  description <- paste(
    "Target-mediated drug disposition (TMDD) model for exendin-4 (exenatide)",
    "in healthy human subjects (Gao 2012). Same structure as the companion",
    "rat and monkey models: free drug in the central compartment binds",
    "glucagon-like peptide 1 receptor (GLP-1R) with second-order rate",
    "constant kon to form a drug-receptor complex that either dissociates",
    "(koff) or is internalised and degraded (kint); free drug additionally",
    "distributes to a nonspecific tissue compartment (k12 / k21) and is",
    "eliminated linearly from plasma (kel). The resulting CLc = kel * Vc =",
    "1.48 mL/min/kg is almost identical to glomerular filtration rate",
    "(125 mL/min in a 70 kg adult). The central volume is reported per",
    "kilogram, so body weight (WT) is a required covariate. The human data",
    "cover only subcutaneous doses plus one continuous intravenous infusion",
    "over a narrow dose range, so the total receptor pool Rtot could not be",
    "estimated and was FIXED at 1240 pmol/L, and every parameter is imprecise",
    "(CV 168-587%); the binding parameters in particular should not be read",
    "as reliable affinity estimates. The absorption rate constant was",
    "estimated separately per subcutaneous dose group and decreases with",
    "dose. Naive-pooled fit to mean profiles in ADAPT II; there is no",
    "between-subject variability and the residual-variance parameters were",
    "not reported. The authors themselves note that a population modelling",
    "approach with more extensive data would improve these estimates."
  )
  reference <- paste(
    "Gao W, Jusko WJ (2012).",
    "Target-mediated Pharmacokinetic and Pharmacodynamic Model of Exendin-4",
    "in Rats, Monkeys, and Humans.",
    "Drug Metabolism and Disposition 40(5):990-997.",
    "doi:10.1124/dmd.111.042291.",
    "PMCID PMC3336795.",
    "Human concentration-time data are from Kolterman OG et al. (2005)",
    "Am J Health Syst Pharm 62(2):173-181 (studies A and B) and",
    "Degn KB et al. (2004) Diabetes 53(9):2397-2403 (study C).",
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
      description        = "Body weight. Required because Gao 2012 Table 3 reports the human central volume of distribution per kilogram (Vc = 111 mL/kg); the model multiplies it by WT to obtain an absolute volume in litres.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 88.5 kg, the mean body weight of study A in Gao 2012 Materials and Methods (88.5 +/- 9.4 kg; study B was 88.8 +/- 12.1 kg and study C reported body mass index 21-29 kg/m^2 rather than weight). Time-fixed. Volume scales linearly with WT (exponent 1), which is how the source reports it, not as an estimated allometric exponent.",
      source_name        = "b.wt."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "exenatide", units = "pmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "exenatide", units = "pmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "exenatide", units = "pmol", specimen = "tissue", verified = TRUE),
    complex     = list(analyte = "exenatide/GLP-1R complex", units = "pmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 27L,
    n_studies     = 3L,
    weight_range  = "88.5 +/- 9.4 kg (study A) and 88.8 +/- 12.1 kg (study B); study C reported body mass index 21-29 kg/m^2 and no body weight",
    disease_state = "Healthy subjects.",
    dose_range    = paste(
      "Study A: single subcutaneous doses of 0.1, 0.2, 0.3 and 0.4 ug/kg",
      "(n = 8). Study B: single subcutaneous doses of 0.02, 0.05 and",
      "0.1 ug/kg (n = 8). Study C: continuous intravenous infusion of",
      "0.066 pmol/kg/min (0.276 ng/kg/min) for 360 min (n = 11)."
    ),
    notes = paste(
      "Gao 2012 Materials and Methods, 'Data from three human studies were",
      "included in the current analysis'. Studies A and B are from Kolterman",
      "2005 and study C from Degn 2004. Plasma exendin-4 was measured by",
      "Amylin Pharmaceuticals using an immunoenzymetric assay. Model fitted",
      "to MEAN profiles by naive pooling in ADAPT II, so n_subjects is the",
      "sum of the three study sizes rather than a modelled population."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural TMDD parameters -- Gao 2012 Table 3, 'Humans' column
    # (CV% in parentheses in the source). All rates are 1/min. The human
    # model was run on the PICOMOLAR concentration scale (kon in
    # 1/(pM*min) and Rtot in pM), confirmed by koff/kon = 0.566/0.000411 =
    # 1377 pM = 1.38 nM, the KD quoted in Gao 2012 Discussion
    # 'Pharmacokinetics'.
    #
    # Note on rounding: Table 3 prints kel to two significant figures
    # (0.013), and 0.013 * 111 mL/kg = 1.44 mL/min/kg, whereas Results
    # 'Human PK' quotes CLc = 1.48 mL/min/kg (which implies an unrounded
    # kel of 0.01333). The printed table value is used here; the ~3%
    # difference is pure display rounding, not a discrepancy.
    # ------------------------------------------------------------------
    lkel  <- log(0.013)     ; label("Linear elimination rate constant from the central compartment (kel, 1/min)")     # Table 3 Humans: kel = 0.013 (CV 198%)
    lk12  <- log(0.0685)    ; label("Transfer rate constant central -> peripheral1 (kpt, 1/min)")                     # Table 3 Humans: kpt = 0.0685 (CV 262%)
    lk21  <- log(0.0846)    ; label("Transfer rate constant peripheral1 -> central (ktp, 1/min)")                     # Table 3 Humans: ktp = 0.0846 (CV 168%)
    lvc   <- log(0.111)     ; label("Central volume of distribution per kilogram (Vc, L/kg)")                         # Table 3 Humans: Vc = 111 mL/kg = 0.111 L/kg (CV 168%)
    lkon  <- log(0.000411)  ; label("Second-order association rate constant of exendin-4 with GLP-1R (kon, 1/(pM*min))") # Table 3 Humans: kon = 0.000411 (CV 351%)
    lkoff <- log(0.566)     ; label("First-order dissociation rate constant of the drug-receptor complex (koff, 1/min)") # Table 3 Humans: koff = 0.566 (CV 399%)
    lkint <- log(0.00342)   ; label("Internalisation / degradation rate constant of the drug-receptor complex (kint, 1/min)") # Table 3 Humans: kint = 0.00342 (CV 587%)
    lrtot <- fixed(log(1240)); label("Total GLP-1R concentration, held constant (Rtot, pmol/L)")                      # Table 3 Humans: Rtot = 1240 pM (FIXED); Results 'Human PK': Rtot "was unable to be estimated; thus, this parameter was fixed as a computer generalized value of 1.24 nM"

    # ------------------------------------------------------------------
    # Subcutaneous absorption. Gao 2012 Table 3 reports only the RANGE of
    # the per-dose-group estimates for humans, "ka = 0.00550-0.0148 1/min
    # (CV 14-17%)", with Results 'Human PK' stating that "the absorption
    # rate constants were lower at higher doses, as was also found in other
    # species". Figure 7A plots the individual points against dose. The
    # value below is the upper end of the printed range, which Fig. 7A
    # assigns to the two lowest subcutaneous doses (0.02 and 0.05 ug/kg);
    # the lower end (0.00550) corresponds to the two highest doses
    # (0.3 and 0.4 ug/kg). Substitute the value for the dose group being
    # simulated.
    # ------------------------------------------------------------------
    lka     <- log(0.0148)   ; label("First-order subcutaneous absorption rate constant at the lowest doses (ka, 1/min)") # Table 3 Humans: ka range 0.00550-0.0148 (CV 14-17%); upper end assigned to the lowest doses per Fig. 7A
    lfdepot <- fixed(log(1)) ; label("Absolute subcutaneous bioavailability (F, fraction)")                           # Table 3 Humans: F = 1 (fixed); Results 'Human PK': "Bioavailability was fixed at 1"

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

    # 5. Observation: free exendin-4 concentration in plasma (pmol/L).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
