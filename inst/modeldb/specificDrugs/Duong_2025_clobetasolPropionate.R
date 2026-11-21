Duong_2025_clobetasolPropionate <- function() {
  description <- paste(
    "One-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "(version 22) full-PBPK model for the topical corticosteroid clobetasol",
    "propionate in healthy adults (Duong 2025). The source model couples a",
    "full-PBPK distribution layer to Simcyp's multi-phase multi-layer",
    "MechDermA dermal absorption model; neither the whole-body mass-balance",
    "equations nor the skin-physiology database values behind MechDermA are",
    "published, so the platform model itself cannot be encoded here. What IS",
    "reported is the intravenous disposition layer: the in vivo systemic",
    "clearance input of 27 L/h, its printed renal component of 0.12 L/h and",
    "the final model's steady-state volume of distribution of 2.34 L/kg",
    "(Table 6). At the 70 kg reference weight that volume is the 164 L Study",
    "A value the Kp scalar was optimised to reproduce. Elimination is carried",
    "as an explicit renal plus non-renal sum. No parameter is fitted here and",
    "none is imported from a Simcyp population file. Because the source used",
    "a full-PBPK distribution model, only a lumped steady-state volume is",
    "reported and this reduction is mono-exponential: it reproduces the",
    "paper's own predicted AUC0-30h and AUC0-inf after a 2 mg intravenous",
    "dose (Table 7) to within 8 percent and its terminal half-life to within",
    "20 percent. Only the intravenous arm is reproduced: the topical cream",
    "and ointment arms, the atopic-dermatitis skin-barrier-impairment",
    "scenarios and the oral arm rest on platform layers the paper does not",
    "print; see the vignette for the quantitative reasons. This is a",
    "typical-value simulation model: the source reports no inter-individual",
    "variance components and no residual-error model, so there are no etas",
    "and propSd is fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Duong JK, van Dijkman S, Ong G, Marta A, Ceci A, Bonifazi E,",
    "Della Pasqua O. (2025). Characterizing local and systemic exposure to",
    "clobetasol propionate in healthy subjects and patients with atopic",
    "dermatitis. Br J Clin Pharmacol 91(11):3150-3166.",
    "doi:10.1002/bcp.70102.",
    sep = " "
  )
  vignette <- "Duong_2025_clobetasolPropionate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What the single ODE state holds, in what amount units, in what matrix.
  # Duong 2025 'Intravenous administration' (Methods): 'A full PBPK
  # distribution model was used ... the tissue-plasma partition coefficient
  # (Kp) was predicted using the Rodgers and Roland Method'; Table 6 prints
  # only the lumped Vss and the Kp scalar.
  compartmentData <- list(
    central = list(
      analyte = "clobetasol propionate",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # Body weight is implicit in the L/kg volume input but is never printed.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Duong 2025 Table 6 expresses Vss in L/kg, so the",
        "Simcyp model scales distribution volume with weight. This",
        "reduction fixes a 70 kg reference weight instead of carrying a",
        "weight term, because the printed clearance of 27 L/h is an",
        "absolute value whose weight scaling (if any) lives in the",
        "unpublished Simcyp Healthy Volunteer population file. Scaling",
        "volume with weight while holding clearance fixed would produce an",
        "internally inconsistent model, so neither is scaled."
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. No body weight is",
        "printed for Study A anywhere in the paper or its supplement. The",
        "70 kg value is the standing rounded standard and is corroborated",
        "by the paper itself: the Kp scalar was optimised to match the",
        "observed Study A Vss of 164 L, and 164 L / 2.34 L/kg = 70.1 kg.",
        "See the vignette Assumptions and deviations."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "Not reported for Study A (healthy adult males; Table 1)",
    weight_median = "70 kg (assumed reference weight for the L/kg volume input; not reported, but implied by 164 L / 2.34 L/kg)",
    sex_female_pct = 0,
    race_ethnicity = "Not reported",
    disease_state = "Healthy adult male volunteers (intravenous Study A).",
    dose_range = "Single 2 mg intravenous dose of clobetasol propionate (Study A, Table 1).",
    route = "intravenous",
    regions = "United Kingdom (Glaxo Group Research report GMH/83/006, Greenford, 1983).",
    studies = paste(
      "The intravenous layer encoded here comes from Study A (Table 1;",
      "reference 15, an unpublished 1983 Glaxo bioavailability report): six",
      "healthy males given 2 mg clobetasol propionate intravenously and",
      "orally, sampled to 30 h with urine collected to 30 h. The paper as a",
      "whole used 12 clinical studies: Study A, five topical cream or",
      "ointment studies in healthy males (Studies B-D and Harding 1985), a",
      "stratum-corneum tape-stripping study (Au 2010, n = 30) and three",
      "studies in atopic dermatitis (Hehir 1983, Sparidans 2010, van Velsen",
      "2012; Tables 1 and 2)."
    ),
    notes = paste(
      "This is a PBPK analysis rather than a population-PK fit, so there is",
      "no pooled analysis dataset and no estimated variance components; the",
      "5th-95th percentile bands in Figures 2-6 are the spread of a Simcyp",
      "virtual population (trial designs in Supplementary Tables S4-S6).",
      "Compound layer (Table 6): MW 467 g/mol, log P 3.5, neutral, blood:",
      "plasma ratio 1.2 and fu 0.07 (both QSAR predicted), Kp scalar 0.22,",
      "predicted Vss 2.34 L/kg, in vivo CL 27 L/h, renal CL 0.12 L/h. The",
      "renal arm is encoded as lcl_renal, the rest as lcl_nonren. Reported",
      "but NOT reproducible from this model: the MechDermA topical",
      "absorption of the 0.05% cream and ointment (Tables 3, 4 and 6:",
      "Pcell 1e-7 cm/h, KpSClip:v 1.77, DSC,lip 3.03e-4 cm2/h for cream and",
      "7.0e-3 cm2/h for ointment or under occlusion, Kpdis:cont 20, cream",
      "evaporation profile), the atopic-dermatitis skin-barrier-impairment",
      "model (Table 5: 80% / 26% reduction in skin thickness for lesional /",
      "non-lesional skin, with hydration, pH and corneocyte-size changes),",
      "the stratum-corneum concentrations, and the oral first-order",
      "absorption arm (Supplementary Table S1: fa 0.5, ka 0.86 1/h)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Duong 2025 entries or
    # arithmetic consequences of them. The one assumption is the 70 kg
    # reference body weight needed to turn the L/kg volume into litres;
    # no body weight is printed (see covariatesDataExcluded$WT).
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Methods ('Intravenous
    # administration') state that a full PBPK distribution model with
    # Rodgers and Rowland partition coefficients was used and the Kp
    # scalar was optimised to the Study A Vss. Table 6 prints the Kp
    # scalar (0.22) and the predicted Vss but no tissue Kp, organ volume
    # or blood flow, so the multi-tissue structure cannot be rebuilt; the
    # only lumped distribution quantity is the Vss. The reduction is
    # therefore mono-exponential. Its agreement with the paper's own
    # predicted Study A AUC and half-life (Table 7) is in the vignette and
    # nothing was tuned to reach it.
    # ------------------------------------------------------------------

    # Central volume. Table 6: 'Vss (L/kg)' = 2.34. At the 70 kg
    # reference weight: 2.34 L/kg * 70 kg = 163.8 L, which matches the
    # observed Study A Vss of 164 L that the Kp scalar was optimised to
    # (Methods, 'Intravenous administration').
    lvc <- fixed(log(163.8))
    label("Central compartment volume vc (L) at the 70 kg reference weight")

    # Renal arm of clearance. Table 6: 'CLR (L/h)' = 0.12 (Study A); also
    # the abstract and Table 7 observed CLR 0.12 L/h.
    lcl_renal <- fixed(log(0.12))
    label("Renal plasma clearance CL_renal (L/h)")

    # Non-renal arm of clearance, by difference from Table 6
    # 'CLIV (L/h)' = 27 ('In vivo Clearance' type; abstract: mean clearance
    # 27 L/h): 27 - 0.12 = 26.88 L/h.
    lcl_nonren <- fixed(log(26.88))
    label("Non-renal plasma clearance CL_nonren (L/h)")

    # Duong 2025 is a PBPK simulation analysis, not a population-PK fit. It
    # reports no residual-error model and no inter-individual variance
    # components. Rather than invent a variance, the residual error is
    # fixed at zero, which makes this a deterministic typical-value
    # simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    #    Total plasma clearance is the sum of the two arms and reproduces
    #    the 27 L/h of Table 6 by construction.
    # ------------------------------------------------------------------
    vc <- exp(lvc)
    cl_renal <- exp(lcl_renal)
    cl_nonren <- exp(lcl_nonren)
    cl <- cl_renal + cl_nonren

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system. A single compartment, as justified in the ini()
    #    header. Amounts are in mg. Intravenous administration only.
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is
    #    in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units of
    #    Duong 2025 Figure 2 and Table 7.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
