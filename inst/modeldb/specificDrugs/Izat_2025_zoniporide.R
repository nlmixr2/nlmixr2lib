Izat_2025_zoniporide <- function() {
  description <- paste(
    "One-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "full-PBPK model for the sodium-hydrogen exchanger 1 inhibitor",
    "zoniporide in healthy adults (Izat 2025). The source model was built",
    "in the Simcyp simulator (version 21) and its whole-body mass-balance",
    "equations are not published, so the platform model itself cannot be",
    "encoded here. What IS fully reported is the zoniporide compound",
    "layer: Table 1 gives the steady-state volume of distribution",
    "recovered by the fitted global tissue-partition scalar, the observed",
    "total intravenous plasma clearance, and the observed renal clearance,",
    "so the elimination is carried as an explicit renal plus non-renal",
    "sum. Zoniporide is the one compound in the paper that is not a CYP3A4",
    "substrate: its clearance is essentially all aldehyde oxidase plus",
    "renal excretion, and its clearance exceeds hepatic blood flow, which",
    "is why the paper predicts an extrahepatic (renal) aldehyde oxidase",
    "contribution of up to 18 percent. No parameter is fitted here and",
    "none is imported from a Simcyp population file: every value is either",
    "a Table 1 entry or an arithmetic consequence of one at the assumed",
    "70 kg reference body weight. Because the source used a full-PBPK",
    "distribution model, only a single lumped steady-state volume is",
    "reported and this reduction is mono-exponential, which overpredicts",
    "the observed peak concentration by about 14 percent; see the vignette",
    "for the quantitative comparison. This is a typical-value simulation",
    "model: the source reports no inter-individual variance components and",
    "no residual-error model, so there are no etas and propSd is fixed at",
    "zero. The fraction metabolized by aldehyde oxidase that is the",
    "paper's main contribution is a property of the unpublished platform",
    "layer and is NOT reproducible from this model; see the vignette for",
    "the full list of deviations.",
    sep = " "
  )
  reference <- paste(
    "Izat N, Bolleddula J, Carione P, Huertas Valentin L, Jones RS,",
    "Kulkarni P, Moss D, Peterkin VC, Tian DD, Treyer A, Venkatakrishnan K,",
    "Zientek MA, Barber J, Houston JB, Galetin A, Scotcher D. (2025).",
    "Establishing a physiologically based pharmacokinetic framework for",
    "aldehyde oxidase and dual aldehyde oxidase-CYP substrates.",
    "CPT Pharmacometrics Syst Pharmacol 14(1):164-178.",
    "doi:10.1002/psp4.13255.",
    sep = " "
  )
  vignette <- "Izat_2025_aldehyde_oxidase_substrates"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what the single ODE state holds, in what amount units, in
  # what biological matrix. Verified against Izat 2025 Table 1 (zoniporide
  # column, "Distribution model: Full PBPK") and Appendix S1 Supplement 2,
  # which states that the global Kp scalar was fitted "to recover the
  # observed volume of distribution at steady state (Vss)" against observed
  # plasma concentrations.
  compartmentData <- list(
    central = list(
      analyte = "zoniporide", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  # Body weight is implicit in the volume input but is never printed.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Izat 2025 Table 1 expresses Vss in L/kg, so the",
        "Simcyp model does scale distribution volume with weight. This",
        "reduction fixes a 70 kg reference weight instead of carrying a",
        "weight term, because the corresponding weight scaling of",
        "clearance lives in the unpublished Simcyp Healthy Volunteers",
        "population file. Scaling volume with weight while holding",
        "clearance fixed would produce an internally inconsistent model,",
        "so neither is scaled."
      ),
      units = "kg",
      type  = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. No body weight",
        "is printed anywhere in the paper or Appendix S1 -- Table S3",
        "reports only age range and female fraction. The 70 kg value is",
        "the standing rounded-standard assumption; see the vignette",
        "Errata."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4L,
    n_studies      = 1L,
    age_range      = "18-55 years",
    weight_median  = "70 kg (assumed reference weight for the L/kg volume input; not reported)",
    sex_female_pct = 0,
    disease_state  = "Healthy adult volunteers.",
    dose_range     = "Single 80 mg intravenous infusion.",
    route          = "intravenous",
    studies        = paste(
      "Dalvie D, Zhang C, Chen W, Smolarek T, Obach RS, Loi CM.",
      "Cross-species comparison of the metabolism and excretion of",
      "zoniporide: contribution of aldehyde oxidase to interspecies",
      "differences. Drug Metab Dispos 2010;38(4):641-654, cited as",
      "reference 13 of Appendix S1 Table S3. Four healthy men aged 18-55",
      "years received a single 80 mg intravenous dose in a combined",
      "pharmacokinetic and mass-balance study; this single study supplies",
      "both the disposition parameters reproduced here and the observed",
      "fractions metabolized."
    ),
    notes          = paste(
      "This is a PBPK analysis rather than a population-PK fit, so there",
      "is no pooled analysis dataset and no estimated variance",
      "components. The percent-CV figures in Table 1 and Appendix S1",
      "Table S10 are the spread of a Simcyp virtual population driven by",
      "unpublished population files, not estimated omegas.",
      "Observed dispositions from the human mass-balance study (Appendix",
      "S1 Table S2) are fmAO 0.52-0.69, fmother 0.13-0.30, renal",
      "excretion 0.17 and biliary excretion 0.01; zoniporide is the one",
      "compound in the paper with no CYP3A4 contribution. Those fractions",
      "are recorded here for provenance and their renal component is",
      "encoded as lcl_renal, but the metabolic pathway split itself is",
      "not encoded because it does not affect the plasma",
      "concentration-time profile and is a property of the unpublished",
      "platform layer."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Izat 2025 Table 1
    # entries or arithmetic consequences of them. The one assumption is
    # the 70 kg reference body weight needed to turn the L/kg volume
    # input into litres; no body weight is printed anywhere in the paper
    # (see covariatesDataExcluded$WT and the vignette Errata).
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Izat 2025 Table 1 classifies
    # zoniporide under "Distribution model: Full PBPK" with a fitted
    # global Kp scalar of 0.64. In a full-PBPK layout the multi-tissue
    # structure IS the distribution model, and the only distribution
    # quantity Table 1 reports is the lumped steady-state volume that the
    # Kp scalar was tuned to recover; no individual tissue volumes,
    # blood flows or partition coefficients are printed, and Table 1
    # gives neither a Vsac nor a Q for this compound. The reduction is
    # therefore mono-exponential. This IS an approximation of a genuinely
    # multi-exponential profile, and it costs about 14 percent on the
    # peak concentration. That deviation is reported in the vignette
    # rather than tuned away.
    # ------------------------------------------------------------------

    # Central volume. Table 1 gives the steady-state volume of
    # distribution Vss = 1.70 L/kg. At the 70 kg reference weight:
    #   1.70 L/kg * 70 kg = 119.0 L
    # There is no Vsac to subtract for this compound.
    lvc <- fixed(log(119.0))
    label("Central compartment volume vc (L) at the 70 kg reference weight")
    # Derived from Izat 2025 Table 1 (Vss 1.70 L/kg, recovered by the
    # optimized global Kp scalar of 0.64).

    # Renal arm of clearance. Table 1: CLRenal = 16.2 L/h. This is the
    # one compound in the paper with a non-negligible renal clearance,
    # so the elimination is carried as the canonical two-component
    # lcl_renal + lcl_nonren sum rather than a single lcl.
    lcl_renal <- fixed(log(16.2))
    label("Renal plasma clearance CL_renal (L/h)")
    # Izat 2025 Table 1, zoniporide column, "CLRenal (L/h)": 16.2.
    # Cross-check: 16.2 / 95.35 = 0.170, which reproduces the observed
    # renal excretion fraction ferenal of 0.17 in Appendix S1 Table S2 --
    # confirming that the renal clearance is a component of the reported
    # total rather than an additional term.

    # Non-renal arm of clearance, by difference from the reported total
    # intravenous plasma clearance:
    #   95.35 L/h total - 16.2 L/h renal = 79.15 L/h
    # Table 1: CLiv = 95.35 L/h [5.7% CV] (i.v.).
    lcl_nonren <- fixed(log(79.15))
    label("Non-renal plasma clearance CL_nonren (L/h)")
    # Derived from Izat 2025 Table 1, zoniporide column: CLiv 95.35 L/h
    # minus CLRenal 16.2 L/h. Cross-checked against Appendix S1 Table S10,
    # which reports observed CLiv 95.4 (5.71) and observed AUCinf 839
    # ng.h/mL after 80 mg, consistent with 80 mg / 95.35 L/h = 839
    # ng.h/mL.

    # Izat 2025 is a PBPK simulation analysis, not a population-PK fit.
    # It reports no residual-error model and no inter-individual variance
    # components. Rather than invent a variance, the residual error is
    # fixed at zero, which makes this a deterministic typical-value
    # simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    #    Total plasma clearance is the sum of the two reported arms and
    #    reproduces the Table 1 total of 95.35 L/h by construction.
    # ------------------------------------------------------------------
    vc         <- exp(lvc)
    cl_renal   <- exp(lcl_renal)
    cl_nonren  <- exp(lcl_nonren)
    cl         <- cl_renal + cl_nonren

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system. A single compartment, as justified in the ini()
    #    header. Amounts are in mg. Intravenous administration only --
    #    zoniporide has no oral arm anywhere in the paper (Table 1 lists
    #    "n/a" for its absorption model, fa, ka and tlag).
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is
    #    in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    #    throughout Izat 2025 Appendix S1 Table S10.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
