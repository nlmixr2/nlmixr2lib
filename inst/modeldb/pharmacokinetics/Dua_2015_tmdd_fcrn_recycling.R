Dua_2015_tmdd_fcrn_recycling <- function() {
  description <- "QSP. Full (explicit-binding) target-mediated drug disposition (TMDD) model with FcRn-mediated endosomal recycling (after Xiao 2010, as implemented in Dua 2015 Supplementary Model 7). Free drug binds a turning-over target in the central compartment to form a complex; both free drug and complex distribute to a tissue compartment and are taken up into an endosomal space, where each binds a constant pool of FcRn. FcRn-bound drug and FcRn-bound complex are recycled back to the central compartment, while unbound drug and complex in the endosome are eliminated. Nine ODE states, all concentrations in nM; the dose is given as an initial central concentration. Defaults are the tutorial's illustrative parameter set for a generic monoclonal antibody, not fitted to any dataset, so no drug and no real population is attached. Sibling in this library: modellib('Dua_2015_tmdd_multitarget_qss')."
  reference <- paste(
    "Dua P, Hawkins E, van der Graaf PH. A Tutorial on Target-Mediated Drug",
    "Disposition (TMDD) Models. CPT Pharmacometrics Syst Pharmacol.",
    "2015;4(6):324-337. doi:10.1002/psp4.41.",
    "Model is described in the main text ('FcRn recycling'); the equations and",
    "parameter values are Supplementary Materials 'Berkeley Madonna TMDD model",
    "library', Model 7. Model structure after Xiao JJ, et al. Pharmacokinetics of",
    "anti-hepcidin monoclonal antibody Ab 12B9m and hepcidin in cynomolgus",
    "monkeys. AAPS J. 2010;12(4):646-657. doi:10.1208/s12248-010-9222-0.",
    sep = " "
  )
  vignette <- "Dua_2015_tmdd_tutorial"

  # The four endosomal states have no canonical name: the registered endosomal
  # families (eu_/eb_/fr_<organ>, and the pH-resolved endo*_<organ> family) are
  # whole-body PBPK sub-compartments that require an organ suffix, whereas this
  # model has a single lumped endosomal space shared by free drug and complex.
  paper_specific_compartments <- c("endosome", "endosome_fcrn", "complex_endosome", "complex_endosome_fcrn")

  # Time is days and every state is a concentration in nM. Supplement Model 7
  # carries no volume: transfers between compartments are symmetric first-order
  # rate constants (e.g. Lt' = kpt*Lc - ktp*Lt), so every state is referenced to
  # the central volume, and the dose is the initial central concentration
  # Lc0 = 200 nM. A dose of amt = 200 into `central` reproduces it.
  units <- list(time = "day", dosing = "nM", concentration = "nM")

  compartmentData <- list(
    central = list(analyte = "drug (free)", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "drug (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target = list(analyte = "target (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex = list(analyte = "drug-target complex", units = "nM", specimen = "plasma", verified = TRUE),
    complex_peripheral1 = list(analyte = "drug-target complex", units = "nM", specimen = "tissue", verified = TRUE),
    endosome = list(analyte = "drug (free, endosomal)", units = "nM", specimen = "not applicable", verified = TRUE),
    endosome_fcrn = list(
      analyte = "drug-FcRn complex (endosomal)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_endosome = list(
      analyte = "drug-target complex (endosomal, not FcRn-bound)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_endosome_fcrn = list(
      analyte = "drug-target complex bound to FcRn (endosomal)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  dosing <- "central"

  covariateData <- list()

  population <- list(
    species = "not applicable (theoretical model)",
    n_subjects = NA_integer_,
    disease_state = "none; a teaching implementation of a generic monoclonal antibody with target binding and FcRn recycling",
    dose_range = "single i.v. bolus giving an initial central concentration of 200 nM (Supplement Model 7: Lc0 = 200 nM)",
    notes = paste(
      "TUTORIAL model. Dua 2015 is a review/tutorial of TMDD model structures; Supplementary",
      "Model 7 is the authors' runnable Berkeley Madonna implementation of a full TMDD model",
      "extended with FcRn-mediated endosomal recycling, after Xiao 2010. The rate constants",
      "are an illustrative set (the target-binding constants are those of Supplement Model 1,",
      "the distribution constants those of Supplement Models 2-6) and are not fitted to any",
      "data. Every value is printed in the supplement and is encoded with fixed(). There is",
      "no inter-individual variability and no residual error in the source; propSd is fixed(0)."
    )
  )

  ini({
    # ---- Linear disposition (Supplement Model 7, 'Values of Rate Constants') --
    lkel <- fixed(log(0.024)); label("Elimination rate constant of free drug from central (1/day)") # Supplement Model 7: kelc = 0.024 1/day
    lkel_endo <- fixed(log(0.024)); label("Elimination rate constant of non-FcRn-bound drug and complex from the endosome (1/day)") # Supplement Model 7: kele = 0.024 1/day
    lk12 <- fixed(log(0.902)); label("Central-to-tissue transfer rate constant of drug and complex (1/day)") # Supplement Model 7: kpt = 0.902 1/day
    lk21 <- fixed(log(0.725)); label("Tissue-to-central transfer rate constant of drug and complex (1/day)") # Supplement Model 7: ktp = 0.725 1/day
    lkup <- fixed(log(0.36)); label("Central-to-endosome uptake rate constant of drug and complex (1/day)") # Supplement Model 7: kup = 0.36 1/day

    # ---- Target turnover and binding ------------------------------------------
    lrbase_target <- fixed(log(0.03)); label("Baseline free target concentration R0 (nM)") # Supplement Model 7: R0 = 0.03 nM
    lkdeg <- fixed(log(0.823)); label("First-order degradation rate constant of free target (1/day)") # Supplement Model 7: kout = 0.823 1/day
    lkint <- fixed(log(0.201)); label("Elimination rate constant of the drug-target complex from central (1/day)") # Supplement Model 7: kepc = 0.201 1/day
    lkon <- fixed(log(0.592)); label("Association rate constant, drug to target (1/(nM*day))") # Supplement Model 7: kon = 0.592 1/(nM day)
    lkoff <- fixed(log(0.9)); label("Dissociation rate constant of the drug-target complex (1/day)") # Supplement Model 7: koff = 0.9 1/day

    # ---- FcRn binding and recycling -------------------------------------------
    lrtot_fcrn <- fixed(log(0.05)); label("Endosomal FcRn concentration, held constant Fc0 (nM)") # Supplement Model 7: Fc0 = 0.05 nM
    lkon_fcrn <- fixed(log(4.83)); label("Association rate constant, drug or complex to FcRn (1/(nM*day))") # Supplement Model 7: konA = 4.83 1/(nM day)
    lkoff_fcrn <- fixed(log(2)); label("Dissociation rate constant of the FcRn complexes (1/day)") # Supplement Model 7: koffA = 2 1/day
    lkrec_fcrn <- fixed(log(0.4)); label("Recycling rate constant of FcRn complexes back to central (1/day)") # Supplement Model 7: kr = 0.4 1/day

    # The source is a deterministic teaching implementation with no residual
    # error; fixed(0) keeps the model simulable without inventing a value.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - none in the source)")
  })

  model({
    kel <- exp(lkel)
    kel_endo <- exp(lkel_endo)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kup <- exp(lkup)

    rbase <- exp(lrbase_target)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ksyn <- kdeg * rbase # Supplement Model 7: kin = kout * R0

    fcrn <- exp(lrtot_fcrn)
    kon_fcrn <- exp(lkon_fcrn)
    koff_fcrn <- exp(lkoff_fcrn)
    krec_fcrn <- exp(lkrec_fcrn)

    target(0) <- rbase

    # Supplement Model 7 'Ordinary Differential Equations', in the printed order:
    # Lc, Lt, R, Pc, Pt, Le, FcLe, Pe, FcPe.
    d/dt(central) <- -(kel + kup + k12) * central - kon * target * central + koff * complex +
      krec_fcrn * endosome_fcrn + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(target) <- ksyn - kdeg * target - kon * target * central + koff * complex
    d/dt(complex) <- kon * target * central - (koff + kint + kup + k12) * complex +
      k21 * complex_peripheral1 + krec_fcrn * complex_endosome_fcrn
    d/dt(complex_peripheral1) <- k12 * complex - k21 * complex_peripheral1
    d/dt(endosome) <- kup * central - kel_endo * endosome - kon_fcrn * fcrn * endosome +
      koff_fcrn * endosome_fcrn
    d/dt(endosome_fcrn) <- kon_fcrn * fcrn * endosome - (koff_fcrn + krec_fcrn) * endosome_fcrn
    # The supplement eliminates endosomal complex with kele (the printed Pe' line),
    # not with the separately listed kepe, which no equation uses.
    d/dt(complex_endosome) <- kup * complex - kel_endo * complex_endosome -
      kon_fcrn * fcrn * complex_endosome + koff_fcrn * complex_endosome_fcrn
    d/dt(complex_endosome_fcrn) <- kon_fcrn * fcrn * complex_endosome -
      (koff_fcrn + krec_fcrn) * complex_endosome_fcrn

    # Total drug in the central compartment (free + target-bound).
    ctot <- central + complex

    # Free drug concentration in the central compartment (supplement Lc).
    Cc <- central
    Cc ~ prop(propSd)
  })
}
