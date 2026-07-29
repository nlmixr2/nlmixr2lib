Lamoth_2009_imipenem <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in adult febrile",
    "neutropenic patients with hematological malignancies (Lamoth 2009).",
    "Total clearance is the additive sum of a non-renal arm and a renal arm",
    "linear in Cockcroft-Gault GFR; the central volume of distribution",
    "scales linearly with total body weight referenced to 70 kg. A single",
    "log-normal inter-individual variability term is applied multiplicatively",
    "to the total clearance (TVCL = CL_nonren + CL_renal * GFR / 100), and",
    "residual error is proportional."
  )
  reference <- paste(
    "Lamoth F, Buclin T, Csajka C, Pascual A, Calandra T, Marchetti O.",
    "Reassessment of recommended imipenem doses in febrile neutropenic",
    "patients with hematological malignancies. Antimicrob Agents Chemother.",
    "2009;53(2):785-787. doi:10.1128/AAC.00891-08.",
    sep = " "
  )
  vignette <- "Lamoth_2009_imipenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in Lamoth 2009 (single covariate snapshot",
        "around the imipenem sampling occasion). Enters central volume of",
        "distribution via a linear scaling V_i = V_70kg * (WT / 70) with",
        "fixed exponent 1 (paper Table 2 reports 'V (liters per 70 kg BW)",
        "= 33.5'). Reference weight 70 kg.",
        "Population median 73 kg (range 41-135 kg; paper Table 1)."
      ),
      source_name        = "body weight (paper Table 1; covariate model in Results paragraph 3)"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'GFR' in Lamoth 2009; computed by the Cockcroft-Gault",
        "equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2; the",
        "paper Methods explicitly cites the Cockcroft-Gault formula derived",
        "from age, sex, total body weight, and serum creatinine). Stored",
        "under the canonical CRCL column per",
        "inst/references/covariate-columns.md (CRCL accepts raw",
        "Cockcroft-Gault mL/min when the source paper does not BSA-normalize,",
        "with the per-model description recording the assay form).",
        "Reference value 100 mL/min for the linear renal-CL slope:",
        "cl_renal = exp(lcl_renal) * (CRCL / 100). Population median",
        "105 mL/min (range 38-285 mL/min; paper Table 1)."
      ),
      source_name        = "GFR (paper Methods and Table 1; calculated by the Cockcroft-Gault formula)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 1L,
    age_range      = "17-78 years",
    age_median     = "58 years",
    weight_range   = "41-135 kg",
    weight_median  = "73 kg",
    sex_female_pct = 22.8,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Febrile neutropenic adults with hematological malignancies",
      "(64.9% acute myeloid leukemia, 5.3% acute lymphoblastic leukemia,",
      "7% multiple myeloma, 8.8% lymphoma, 14% other). 47.8% of",
      "chemotherapy courses were induction for acute leukemia, 27.5%",
      "consolidation, and 14.5% autologous stem cell transplantation."
    ),
    dose_range     = paste(
      "Recommended schedule 500 mg imipenem IV infused over 30 min every",
      "6 h (2 g/day total), adjusted to calculated GFR per local renal",
      "dose-adjustment guidance. Observed daily doses 0.75-4 g/day",
      "(median 2 g/day, paper Table 1)."
    ),
    regions        = "Single centre: Centre Hospitalier Universitaire Vaudois and University of Lausanne, Lausanne, Switzerland",
    renal_function = "Cockcroft-Gault GFR median 105 mL/min (range 38-285 mL/min; raw mL/min, not BSA-normalized; paper Table 1)",
    notes          = paste(
      "159 plasma imipenem concentrations (86 troughs + 73 peaks) drawn",
      "around a single dose at steady state (median 3 days after start of",
      "therapy or the last dosing change, range 1-9 days). Trough samples",
      "10 min before the dose; peak samples median 2 h (range 0.5-4 h)",
      "after the start of the 30-min infusion. Free-circulating imipenem",
      "measured by HPLC (analytical range 0.25-200 mg/L; intra-/interassay",
      "accuracy and precision < 5%). 69 chemotherapy courses across the 57",
      "patients (median 2 samples per patient, range 1-10). Baseline",
      "demographics and dosing summary from paper Table 1."
    )
  )

  ini({
    # Structural parameters -- typical values for a 70 kg patient at the
    # paper's reference Cockcroft-Gault GFR of 100 mL/min. All four point
    # estimates come from Lamoth 2009 Table 2 ("Final population
    # pharmacokinetic model for imipenem in febrile neutropenic patients").
    lcl_nonren <- log(10.7); label("Non-renal clearance (L/h)")                              # Lamoth 2009 Table 2: Nonrenal CL = 10.7 +/- 1.7 L/h
    lcl_renal  <- log(4.79); label("Renal CL at CRCL = 100 mL/min (L/h per 100 mL/min GFR)") # Lamoth 2009 Table 2: Renal CL = 4.79 +/- 1.4 L/h per 100 mL/min GFR
    lvc        <- log(33.5); label("Central volume of distribution at 70 kg (L)")            # Lamoth 2009 Table 2: V = 33.5 +/- 4.0 L per 70 kg body weight

    # Linear (allometric exponent 1) body-weight scaling on Vc. The paper
    # Results paragraph 3 states "The definition of V as a multiple of body
    # weight improved the model" and Table 2 reports V in "liters per 70 kg
    # BW" -- i.e. a linear scaling V_i = V_70kg * (WT / 70). Encoded as a
    # fixed exponent so the structural assumption is visible in metadata.
    e_wt_vc <- fixed(1); label("Linear WT exponent on Vc (fixed at 1)")                      # Lamoth 2009 Results paragraph 3 and Table 2 ("V (liters per 70 kg BW)")

    # Inter-individual variability. Lamoth 2009 reports a single IIV
    # magnitude ("17% +/- 6%") in Table 2 alongside the CL parameters; the
    # Results paragraph above Table 2 states the simpler 1-compartment
    # covariate-free model had "26% proportional interpatient variability"
    # on total CL, and that the GFR covariate explained roughly the
    # difference. The 17% is therefore the residual unexplained IIV on
    # total CL after the GFR covariate. Applied here as an exponential ETA
    # on the total CL = (CL_nonren + CL_renal * GFR / 100) -- the standard
    # NONMEM idiom for an additive split-CL covariate model with a single
    # OMEGA on CL. omega^2 = log(1 + CV^2) for log-normal eta.
    etalcl ~ log(1 + 0.17^2)                                                                 # Lamoth 2009 Table 2: BSV on CL = 17% +/- 6% (proportional IIV)

    # Residual error -- proportional model. Lamoth 2009 Table 2 reports
    # "Residual error (59% +/- 7%)" as a proportional CV.
    propSd <- 0.59; label("Proportional residual error (fraction)")                          # Lamoth 2009 Table 2: Residual error = 59% +/- 7% (proportional)
  })
  model({
    # Typical non-renal and renal clearance components. The renal arm is a
    # deterministic linear function of CRCL with reference 100 mL/min; the
    # non-renal arm is constant. The single etalcl term is then applied
    # multiplicatively to total CL = (CL_nonren + CL_renal * CRCL / 100).
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal) * (CRCL / 100)
    cl        <- (cl_nonren + cl_renal) * exp(etalcl)

    # Central volume of distribution: linear (exponent 1) WT scaling at the
    # 70 kg reference. No IIV reported in Table 2.
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment IV-infusion disposition. Dosing into central via the
    # event-table rate column (30-min infusions in Lamoth 2009).
    d/dt(central) <- -kel * central

    # Imipenem plasma concentration; dose in mg, vc in L -> mg/L (matches
    # the HPLC assay reporting units in the paper Methods).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
