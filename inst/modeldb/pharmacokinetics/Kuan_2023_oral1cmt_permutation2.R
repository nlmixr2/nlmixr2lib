Kuan_2023_oral1cmt_permutation2 <- function() {
  description <- paste(
    "Methodology reference. Permutation 2 of the one-compartment,",
    "first-order-input / first-order-output model used by Kuan 2023 to",
    "demonstrate flip-flop pharmacokinetics as a LOCAL IDENTIFIABILITY",
    "problem. There is no drug, no patients and no fitted estimates: the",
    "authors chose F = 1, D = 1, ka = 0.1, k = 0.5 and V = 2 (Appendix S1",
    "Table S1) purely to show that two distinct parameter sets produce the",
    "SAME input-output profile, so every value is encoded with fixed().",
    "This permutation is the ka < k branch (absorption slower than",
    "elimination, i.e. absorption-limited elimination). Its companion",
    "Kuan_2023_oral1cmt_permutation1 carries the ka > k branch and predicts",
    "an identical concentration-time curve. Clearance is the flip-flop",
    "INVARIANT: CL = V * k = 1 in both permutations, which is why",
    "AUC = Dose / CL and noncompartmental analysis are unaffected by whether",
    "a system is in the flip or the flop state. Parameters are stored in the",
    "canonical (CL, V, ka) parameterization; the paper's elimination rate",
    "constant k is recovered as kel = cl / vc = 0.5.",
    sep = " "
  )
  reference <- paste(
    "Kuan IHS, Wright DFB, Duffull SB. The influence of flip-flop in",
    "population pharmacokinetic analyses.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(3):285-287.",
    "doi:10.1002/psp4.12909. PMID 36647235. PMCID PMC10014047.",
    "Parameter values transcribed from Supplement Appendix S1 Table S1",
    "(Permutation 2 column); the permutation algebra is main-text Table 1.",
    sep = " "
  )
  vignette <- "Kuan_2023_flipflop"

  # The source simulation is dimensionless: Table S1 gives D = 1 and V = 2
  # with no units, and t = 0:100 with no time unit. The generic placeholders
  # are therefore the honest encoding (see Beal_2001_iv1cmt_bql, which does
  # the same for the other methodology toy model in the library).
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "dose_unit/vol_unit")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE because the source CANNOT confirm it:
  # Appendix S1 is an abstract mathematical demonstration that names no
  # molecule and no biological matrix. The analyte is nominal.
  compartmentData <- list(
    depot = list(
      analyte = "hypothetical drug", units = "dose_unit",
      specimen = "administration site", verified = FALSE
    ),
    central = list(
      analyte = "hypothetical drug", units = "dose_unit",
      specimen = "plasma", verified = FALSE
    )
  )

  covariateData <- list()

  # Creatinine clearance was screened as a covariate in the paper's MOTIVATING
  # metformin analysis (main text; Appendix S3), not in this Appendix S1 toy
  # model, and the paper reports NO coefficient for it -- only the direction of
  # the objective-function change. Documented here so the provenance of the
  # covariate screen survives without creating an unused-covariate warning.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault estimate, raw mL/min and NOT BSA-normalized. Kuan 2023 Appendix S3 states that IDEAL body weight was used as the body-size metric in the Cockcroft-Gault equation.",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened on ka, on CL, and on both in the metformin analysis of the",
        "main text, under both the (CL, V, ka) and (k, V, ka)",
        "parameterizations, using a 3.84-unit chi-square drop in objective",
        "function value as the significance test (Appendix S3). The paper's",
        "finding is that WHICH parameter CLcr appears significant on depends",
        "on the rank order of the INITIAL estimates of ka and k, i.e. it is an",
        "artefact of local identifiability. No point estimate, standard error",
        "or confidence interval is reported for any covariate coefficient",
        "anywhere in the paper or its supplement, so the metformin model is",
        "not encodable; see the vignette Errata.",
        sep = " "
      ),
      source_name = "CLcr_CG"
    )
  )

  population <- list(
    species        = "None (methodology paper; a single deterministic simulated profile of a hypothetical drug, not a fit of any real molecule).",
    n_subjects     = 1L,
    disease_state  = "N/A (algebraic demonstration of local identifiability).",
    dose_range     = "Single unit dose, D = 1, with F = 1 (Appendix S1 Table S1).",
    regions        = "N/A",
    scope_note     = paste(
      "Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/",
      "because there is no drug, following the precedent set by",
      "Beal_2001_iv1cmt_bql and the Schoning_2026_oral1cmt_* family. The file",
      "stem uses the structural descriptor oral1cmt in the slot where a drug",
      "name would normally go, again as Beal_2001 does with iv1cmt.",
      "IMPORTANT SCOPE LIMIT: this file encodes ONLY the Appendix S1",
      "one-compartment permutation demonstration, which is the sole part of",
      "the paper that reports numeric parameter values. The paper's",
      "motivating metformin population analysis (55 participants, 426 plasma",
      "concentrations pooled from three published studies) is NOT encodable:",
      "no fixed-effect estimate, no between-subject variance and no residual",
      "error value is reported for it anywhere in the paper or its",
      "supplement. The two-compartment permutations of main-text Table 1 are",
      "likewise not encodable -- Table 1 and Appendix S2 give the",
      "permutation algebra symbolically, with no numeric values.",
      sep = " "
    ),
    notes          = paste(
      "Appendix S1: 'Simulations using all the possible permutations for a",
      "one-compartment model was performed to illustrate that the same",
      "input-output profile could be produced with different sets of",
      "parameter values.' Simulated in R 3.5.3 over t = 0:100. The paper",
      "reports no between-subject variability and no residual error for this",
      "simulation, so neither is encoded here; the demonstration is purely",
      "deterministic. For context only, the separate metformin analysis",
      "pooled 55 participants (10 female / 45 male, median age 61.0 years,",
      "median weight 84.6 kg, median CLcr 29.0 mL/min) contributing 426",
      "plasma concentrations from Kuan 2021, Dissanayake 2017 and",
      "Pentikainen 1979 (Appendix S3 Table S2); that cohort is NOT the",
      "population of this model.",
      sep = " "
    )
  )

  ini({
    # Every value is an author-chosen demonstration constant from Appendix S1
    # Table S1, so all are fixed() rather than estimated.
    lka     <- fixed(log(0.1));  label("Absorption rate constant (1/time_unit)")            # Appendix S1 Table S1, Permutation 2: ka = 0.1
    lvc     <- fixed(log(2));    label("Central volume of distribution (vol_unit)")         # Appendix S1 Table S1, Permutation 2: V = 2
    lcl     <- fixed(log(1));    label("Clearance (vol_unit/time_unit)")                    # Appendix S1 Table S1, Permutation 2: CL = V * k = 2 * 0.5 = 1 (main text: "CL remains unchanged in both permutations and is therefore invariant to flip-flop")
    lfdepot <- fixed(log(1));    label("Bioavailability (fraction)")                        # Appendix S1 Table S1, Permutation 2: F = 1
  })

  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    cl  <- exp(lcl)
    kel <- cl / vc  # recovers the paper's elimination rate constant k = 1 / 2 = 0.5 (Appendix S1 Table S1, Permutation 2)

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    Cc <- central / vc
  })
}
