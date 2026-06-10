McLachlan_1996_fluconazole <- function() {
  description <- "One-compartment population PK model for fluconazole in adults with HIV/AIDS, fit to 770 plasma concentrations from 113 male subjects pooled across an intensive-sampling sub-study (Study 1, n=13, 12-17 samples per dose) and a sparse routine-care sub-study (Study 2, n=100, single sample per subject). Oral capsules (Diflucan, 50-800 mg) and 50 mg per 15 min IV infusions are described by a single linear central compartment with first-order absorption from a depot and zero-order input during IV infusion. The final NONMEM clearance model is an additive intercept-plus-slopes regression on Cockcroft-Gault creatinine clearance (raw, not BSA-normalized) and absolute CD4+ T-lymphocyte count: CL (L/h) = 0.25 + 0.0057 * CLcr (mL/min) + 0.00068 * CD4 (cells/mm^3); volume of distribution, absorption rate, and bioavailability are not modulated by covariates."
  reference <- "McLachlan AJ, Tett SE. Pharmacokinetics of fluconazole in people with HIV infection: a population analysis. Br J Clin Pharmacol. 1996;41(4):291-298. doi:10.1046/j.1365-2125.1996.03085.x"
  vignette <- "McLachlan_1996_fluconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column 'CLcr' in the paper text. Estimated from serum creatinine, age and weight via Cockcroft-Gault (paper Methods 'Study population', citing Cockcroft & Gault 1976); raw mL/min, not BSA-normalized. Cohort mean 68 mL/min (range 36-138 mL/min) per the Methods text. Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source does not apply BSA normalisation; precedent in Delattre_2010_amikacin.R). Enters CL as a raw additive linear slope (no centring, no divisive normalization): CL = exp(lcl) + e_crcl_cl * CRCL + e_cd4_abs_cl * CD4_ABS, per the final NONMEM equation in the abstract and paper Discussion p.296.",
      source_name        = "CLcr"
    ),
    CD4_ABS = list(
      description        = "Absolute CD4+ T-lymphocyte count (baseline)",
      units              = "cells/mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column 'CD4 cell count' in the paper text. Cohort mean 69 cells/mm^3 per the Methods (severely immunosuppressed HIV/AIDS cohort); the categorical-split analysis (Table 3) recorded 97 of 109 covariate-evaluable subjects with CD4 < 200 cells/mm^3 and 12 of 109 with CD4 >= 200 cells/mm^3. The paper does not tabulate the absolute-count range. Used as a surrogate disease-severity marker in HIV/AIDS; the paper hypothesises that lower CD4 (more advanced immunosuppression) co-occurs with reduced renal tubular reabsorption of fluconazole and altered drug-metabolising enzyme activity (paper Discussion p.297). Enters CL as a raw additive linear slope (cells/mm^3 in absolute count, not z-score): CL = exp(lcl) + e_crcl_cl * CRCL + e_cd4_abs_cl * CD4_ABS.",
      source_name        = "CD4 cell count"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 113L,
    n_studies      = 2L,
    age_range      = "23-60 years",
    age_median     = "38 years (mean)",
    weight_range   = "42-88 kg",
    weight_median  = "63 kg (mean)",
    sex_female_pct = 0,
    race_ethnicity = "(not reported in the source; single-site cohort recruited at St Vincent's Hospital, Darlinghurst NSW, Australia)",
    disease_state  = "HIV infection / AIDS; majority had a prior or current AIDS-defining illness (21 Kaposi's sarcoma, 35 Pneumocystis carinii pneumonia, 25 CMV retinitis, 12 HSV infection per the Methods). Mean baseline CD4+ T-lymphocyte count 69 cells/mm^3 (severe immunosuppression). Cockcroft-Gault creatinine clearance mean 68 mL/min (range 36-138). None of the subjects were receiving other drugs known to affect fluconazole disposition. Most were on antiretroviral nucleosides (zidovudine, didanosine, zalcitabine) and co-trimoxazole; the paper Results confirm these did not affect CL or V.",
    dose_range     = "Study 1 (n=13, intensive): 50, 100, or 400 mg single doses, each given as both an oral capsule (Diflucan) and an IV infusion at 50 mg per 15 min; doses separated by at least 2 weeks, 12-17 samples per dose at nominal times 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 24, 32, 48, 72, 96, 120, 168 h. Study 2 (n=100, sparse): 50-800 mg oral fluconazole at routine clinical frequency (twice daily to once weekly; majority once daily), all subjects at steady state (>= 7 days on the current regimen), one sample per subject at a precisely recorded post-dose time (mean 20.9 h, range 0.42-51.6 h).",
    regions        = "Australia (St Vincent's Hospital, Darlinghurst NSW)",
    n_observations = 770L,
    notes          = "Baseline demographics and study design from the Methods sections 'Study population' and 'Dose administration and plasma sampling'. All 113 subjects were male; the source recruited only male subjects. The combined dataset (Studies 1 + 2) was used for the final popPK fit reported in Tables 2-4; covariate analysis was run on the n=109 subjects for whom covariate data were available."
  )

  ini({
    # Structural parameters and CL covariate model from the final NONMEM regression
    # (paper Abstract and Discussion p.296):
    #   CL (L/h) = 0.25 (33% CV) + 0.0057 (32% CV) * CLcr (mL/min)
    #                            + 0.00068 (10% CV) * CD4 (cells/mm^3)
    # Volume, absorption rate, and bioavailability are not covariate-modulated;
    # the paper p.295 confirms V, ka, and F were unchanged after covariate
    # adjustment on CL. V, ka, F, and the residual error are taken from Table 2
    # (combined-dataset NONMEM column, n=113, 770 observations).

    lcl <- log(0.25)
    label("CL covariate-model intercept (L/h)")
    # McLachlan 1996 Abstract and Discussion p.296: final NONMEM CL regression
    # intercept = 0.25 L/h. Wrapped in log() so the structural intercept is
    # positive-by-construction inside model() per the operator-resolved sidecar
    # (request-001, option A).

    e_crcl_cl <- 0.0057
    label("Additive CL slope on Cockcroft-Gault CRCL (L/h per mL/min)")
    # McLachlan 1996 Abstract and Discussion p.296: CLcr slope = 0.0057 L/h per
    # mL/min in the final NONMEM additive-linear CL covariate model.

    e_cd4_abs_cl <- 0.00068
    label("Additive CL slope on absolute CD4 count (L/h per cell/mm^3)")
    # McLachlan 1996 Abstract and Discussion p.296: CD4 slope = 0.00068 L/h per
    # cell/mm^3 in the final NONMEM additive-linear CL covariate model.

    lvc <- log(47.6)
    label("Central volume of distribution (L)")
    # McLachlan 1996 Table 2 NONMEM column: V = 47.6 L (combined dataset,
    # 113 subjects, 770 observations).

    lka <- log(5.02)
    label("First-order absorption rate constant (1/h)")
    # McLachlan 1996 Table 2 NONMEM column: ka = 5.02 1/h (absorption half-life
    # ln(2)/5.02 = 0.138 h ~= 8 min, matching the paper's 'absorption half-life
    # 8 to 14 min' statement).

    lfdepot <- log(0.99)
    label("Bioavailability of the oral depot dose (unitless)")
    # McLachlan 1996 Table 2 NONMEM column: F = 0.99 (extent of absorption
    # for the oral capsule; IV infusion doses are 100% bioavailable and bypass
    # the depot).

    # Inter-individual variability (Table 2 NONMEM column, %CV). Per the
    # operator-resolved sidecar (request-001, option A), encoded as log-normal
    # etas: omega^2 = log(CV^2 + 1). The paper's Equation 1 uses a
    # multiplicative-on-linear-scale form (p_j = p_pop * (1 + g_pj)); the
    # log-normal form is the rxode2 / nlmixr2 idiom for positive-constrained
    # structural parameters and is documented as a deviation in the vignette.

    etalcl ~ 0.15532
    # log(0.41^2 + 1); 41% CV from McLachlan 1996 Table 2 base (no-covariate)
    # NONMEM CL IIV. The paper does not report a separate residual CL IIV
    # after covariate adjustment, so the base-model value is carried per the
    # operator-resolved sidecar (conservative upper bound; documented in the
    # vignette Assumptions and deviations section).

    etalvc ~ 0.006380
    # log(0.08^2 + 1); 8% CV from McLachlan 1996 Table 2 NONMEM V IIV.

    etalka ~ 1.86783
    # log(2.34^2 + 1); 234% CV from McLachlan 1996 Table 2 NONMEM ka IIV. The
    # paper Discussion p.295 explicitly states 'NONMEM had difficulty
    # estimating the intersubject variability for the absorption rate constant
    # (ka), probably because of the very few fluconazole concentrations taken
    # prior to the peak concentration'; the P-PHARM ka IIV from the same Table
    # is 41% CV, which the paper presents as more credible. The NONMEM value
    # is retained here because Table 2 NONMEM is the primary analysis used for
    # the other parameters; the vignette Assumptions and deviations section
    # documents the discrepancy.

    etalfdepot ~ 0.003594
    # log(0.06^2 + 1); 6% CV from McLachlan 1996 Table 2 NONMEM F IIV.

    addSd <- 0.52
    label("Additive residual error (mg/L)")
    # McLachlan 1996 Table 2 NONMEM 'Residual error' row: 0.52 mg/L. The paper
    # Discussion p.293 states NONMEM used an additive residual error model
    # (Equation 2: C_ij(t) = f(p_j, t_ij) + e_ij), whereas P-PHARM used a
    # heteroscedastic proportional error model with magnitude 0.07 mg/L; the
    # additive NONMEM form is retained because Table 2 NONMEM is the primary
    # analysis used here.
  })

  model({
    # ----- Individual PK parameters -----
    # CL: additive intercept-plus-linear-slopes regression on CRCL and CD4_ABS,
    # wrapped in exp(etalcl) to give log-normal IIV on the resulting CL value
    # (operator-resolved sidecar option A; same idiom as Delattre_2010_amikacin
    # for an additive-linear renal-function term).
    cl <- (exp(lcl) + e_crcl_cl * CRCL + e_cd4_abs_cl * CD4_ABS) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)

    # First-order elimination micro-constant.
    kel <- cl / vc

    # ----- ODEs -----
    # Oral dose enters via depot (first-order absorption with rate ka and
    # bioavailability F = exp(lfdepot + etalfdepot)). IV infusion doses bypass
    # the depot and target central directly; rate / duration are supplied via
    # the event dataset's rate column.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability on the oral depot dose only.
    f(depot) <- exp(lfdepot + etalfdepot)

    # ----- Observation and residual error -----
    # Concentration in central compartment. Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
