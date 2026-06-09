Gatti_1996_dapsone <- function() {
  description <- paste(
    "One-compartment population PK model with first-order oral absorption",
    "and first-order elimination for dapsone 100 mg twice weekly oral",
    "Pneumocystis carinii pneumonia prophylaxis in 53 HIV-infected adults",
    "(Gatti 1996). Apparent clearance CL/F and apparent central volume",
    "V/F are scaled multiplicatively by concomitant rifampin",
    "co-administration (shared 69.6% increase on both parameters,",
    "reflecting a first-pass / bioavailability effect). Apparent",
    "absorption rate constant Ka is scaled multiplicatively by total",
    "serum bilirubin (per-mg/dL fractional decrease). IIV on CL/F (35%",
    "CV) and Ka (85% CV); V/F inter-individual variability was found",
    "non-significant after covariate inclusion and dropped from the",
    "final model. Residual-error magnitudes were not reported in the",
    "publication; propSd and addSd are FIXED at 0 in this packaged",
    "model so users must supply their own residual error to run any",
    "stochastic VPC -- see the validation vignette's Errata section."
  )
  reference <- paste(
    "Gatti G, Merighi M, Hossein J, Travaini S, Casazza R, Karlsson M,",
    "Cruciani M, Bassetti D.",
    "Population pharmacokinetics of dapsone administered biweekly to",
    "human immunodeficiency virus-infected patients.",
    "Antimicrob Agents Chemother. 1996;40(12):2743-2748.",
    "doi:10.1128/aac.40.12.2743.",
    sep = " "
  )
  vignette <- "Gatti_1996_dapsone"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    CONMED_RIF = list(
      description        = paste(
        "Concomitant rifampin co-administration indicator (1 = on",
        "rifampin for at least 2 weeks at the time of blood sampling,",
        "0 = not on rifampin)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant rifampin)",
      notes              = paste(
        "Time-fixed per subject in the Gatti 1996 cohort (7 of 53",
        "patients on rifampin at study entry, having received rifampin",
        "for at least 2 weeks before blood sampling so the chronic",
        "induction effect on CYP3A4 / first-pass metabolism was at",
        "equilibrium; Materials and Methods paragraph 3). Multiplicative",
        "fractional effect shared between CL/F and V/F: a single theta",
        "(theta4 = 0.696 in Table 3) increases both CL/F and V/F by",
        "69.6% in rifampin co-administered subjects, reflecting a",
        "primarily first-pass / bioavailability effect. The paper",
        "rejected a non-shared parameterisation (dOFV 1.23, P > 0.05;",
        "Results paragraph 3)."
      ),
      source_name        = "R (paper Results paragraph 3 covariate equation)"
    ),
    TBILI = list(
      description        = paste(
        "Total serum bilirubin at the time of pharmacokinetic sampling",
        "(mg/dL). 7 of 53 patients had TBILI > 1.2 mg/dL (the paper's",
        "upper limit of the normal range); the cohort median was 0.7",
        "mg/dL."
      ),
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in the Gatti 1996 cohort. Range 0.3-8.0",
        "mg/dL, median 0.7 mg/dL (Table 1). Multiplicative fractional",
        "effect on Ka analogous to the paper's explicit rifampin",
        "formula CL/F = theta1 * (1 + theta4 * R): Ka = theta3 * (1 +",
        "theta5 * TBILI) with theta5 = -0.119 per mg/dL bilirubin. At",
        "TBILI = 0.7 mg/dL this gives Ka = 1.04 * (1 - 0.119 * 0.7) =",
        "0.953 1/h; the paper's simulation used Ka = 0.957 1/h",
        "(Discussion paragraph 7), with the 0.4% discrepancy attributed",
        "to rounding of theta3 from a precise estimate near 1.043 down",
        "to 1.04 in Table 3 display. The form was confirmed by operator",
        "sidecar response 2026-05-30 (request-001 q1=A)."
      ),
      source_name        = "Bilirubin (paper Materials and Methods Table 1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 53L,
    n_studies      = 1L,
    age_range      = "27-46 years",
    age_median     = "33 years",
    weight_range   = "40-83 kg",
    weight_median  = "62 kg",
    sex_female_pct = 9.4,
    disease_state  = paste(
      "HIV-infected adults receiving dapsone for Pneumocystis carinii",
      "pneumonia (PCP) prophylaxis. CD4 lymphocyte count: median 25",
      "cells/uL (range 0-389); 24 of 53 were p24 antigen negative. 48",
      "primary prophylaxis, 5 secondary. 33 slow acetylators, 20 fast",
      "acetylators (acetylation ratio cutoff 0.35). 40 with intravenous",
      "drug use as HIV-acquisition risk factor; 13 other risk factors.",
      "42 tobacco smokers. 17 patients concomitantly on zidovudine",
      "(AZT); 17 on didanosine (DDI); 7 on rifampin. ALT median 54",
      "IU/L (range 23-508); total bilirubin median 0.7 mg/dL (range",
      "0.3-8.0)."
    ),
    dose_range     = paste(
      "100 mg dapsone twice weekly orally (biweekly: alternating 72-h",
      "and 96-h intervals between doses) for at least 1 month before",
      "the PK study. Two formulations: center 1 used commercial",
      "Farmitalia-Carlo Erba 100 mg tablets (11 patients); center 2",
      "used hospital-pharmacy-compounded tablets from Roussel-supplied",
      "powder (42 patients). Study site (formulation) was tested as a",
      "covariate and found non-significant on all PK parameters."
    ),
    regions        = "Italy (Genoa, Verona) -- 2 medical centers",
    notes          = paste(
      "53 subjects, 218 plasma dapsone concentrations. Sampling target:",
      "random within each of the time intervals 0-3, 3-24, 24-48, 48-96",
      "h after dosing (Materials and Methods). Concentrations measured",
      "by HPLC with an LLOQ of 31.2 ng/mL (0.031 mg/L); inter- and",
      "intra-day RSD < 10% at 62.5, 250, 1000, 5000 ng/mL QC levels.",
      "NONMEM IV (double precision, level 2.0) with FO approximation;",
      "FOCE used to confirm final-model parameters (Materials and",
      "Methods paragraph 5). ADVAN2 TRANS2 parameterisation (Ka, V/F,",
      "CL/F). Demographics from Table 1; final-model parameters from",
      "Table 3. Co-medications and covariate tests in Table 2."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Gatti 1996 Table 3 ('Mean' row).
    # Reference subject: CONMED_RIF = 0 (no rifampin), TBILI = 0 mg/dL.
    # All three structural THETAs are reported as point estimates with
    # 95% CIs in Table 3; none are fixed.
    # ============================================================
    lcl <- log(1.83); label("Apparent oral clearance CL/F for no-rifampin reference (L/h)")
    # Table 3 theta1: CL/F = 1.83 L/h (95% CI 1.57, 2.09)
    lvc <- log(69.6); label("Apparent volume V/F for no-rifampin reference (L)")
    # Table 3 theta2: V/F = 69.6 L (95% CI 57.4, 81.8)
    lka <- log(1.04); label("Absorption rate Ka for TBILI = 0 mg/dL reference (1/h)")
    # Table 3 theta3: Ka = 1.04 1/h (95% CI 0.72, 1.36)

    # ============================================================
    # Covariate effects -- Gatti 1996 Table 3.
    # The rifampin form is written explicitly in the Results paragraph 3:
    #   CL/F = theta1 + theta1*theta4*R = theta1 * (1 + theta4 * R)
    #   V/F  = theta2 + theta2*theta4*R = theta2 * (1 + theta4 * R)
    # with the SAME theta4 enforced for CL/F and V/F (dOFV 1.23, P > 0.05).
    # The bilirubin form on Ka is NOT written as an equation in the paper;
    # operator sidecar response (2026-05-30 request-001 q1=A) selected the
    # multiplicative / fractional form by analogy to the rifampin equation:
    #   Ka = theta3 * (1 + theta5 * TBILI)
    # ============================================================
    e_rif_cl_vc <- 0.696
    label("Rifampin fractional effect shared between CL/F and V/F (unitless)")
    # Table 3 theta4 (= theta6 by likelihood-ratio shared-effect test):
    # 0.696 (95% CI 0.318, 1.074); 69.6% increase in CL/F and V/F
    e_tbili_ka <- -0.119
    label("Bilirubin fractional effect on Ka (per mg/dL TBILI)")
    # Table 3 theta5: -0.119 (95% CI -0.08, -0.158); ~12% decrease in Ka
    # per mg/dL TBILI

    # ============================================================
    # Inter-individual variability -- Gatti 1996 Table 3 ('Mean' CV row)
    # and Results paragraph 1 ('interpatient variability was modeled
    # with a constant coefficient of variation' -- standard log-normal /
    # exponential form, so omega^2 = log(1 + CV^2)).
    #
    # V/F IIV was significant in the basic model (19% CV) but became
    # non-significant after covariate inclusion and was dropped from
    # the final model (Results paragraph 4). No etalvc here.
    # ============================================================
    etalcl ~ log(1 + 0.35^2)
    # Table 3 CV(CL/F) = 35% (95% CI 20, 46); omega^2 = log(1 + 0.35^2)
    etalka ~ log(1 + 0.85^2)
    # Table 3 CV(Ka) = 85% (95% CI 45, 111); omega^2 = log(1 + 0.85^2)

    # ============================================================
    # Residual error -- form is reported (Results paragraph 1:
    # 'proportional-plus-constant-error model') but the numerical
    # magnitudes for sigma_prop and sigma_add are NOT reported anywhere
    # in the publication. Operator sidecar response (2026-05-30
    # request-001 q2=A modified) instructs that both be FIXED at 0 here
    # and the gap be documented prominently in the validation vignette
    # Errata section. Users running stochastic VPCs must supply their
    # own residual error magnitudes; the assay specs (LLOQ 0.031 mg/L,
    # HPLC RSD < 10%) provide a defensible lower bound.
    # ============================================================
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; FIXED at 0 -- not reported in paper)")
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; FIXED at 0 -- not reported in paper)")
  })

  model({
    # ============================================================
    # Individual PK parameters with multiplicative covariate effects.
    # The rifampin form matches the explicit equation in Gatti 1996
    # Results paragraph 3 (shared theta on CL/F and V/F); the
    # bilirubin form matches by analogy per operator sidecar (q1=A).
    # ============================================================
    cl <- exp(lcl + etalcl) * (1 + e_rif_cl_vc * CONMED_RIF)
    vc <- exp(lvc) * (1 + e_rif_cl_vc * CONMED_RIF)
    ka <- exp(lka + etalka) * (1 + e_tbili_ka * TBILI)

    kel <- cl / vc

    # ============================================================
    # One-compartment open model with first-order absorption from the
    # depot (oral) and first-order elimination from the central
    # compartment (Gatti 1996 Results paragraph 1; NONMEM ADVAN2
    # TRANS2 parameterisation). Bioavailability F is structurally
    # absorbed into the apparent CL/F and V/F.
    # ============================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration in mg/L: dose in mg, vc in L, central in mg.
    Cc <- central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
