Por_2021_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 23 US adult",
    "burn patients with (n = 12) and without (n = 11) continuous venovenous",
    "haemofiltration (Por 2021). Body clearance has two branches: patients",
    "not on CVVH carry power effects of Cockcroft-Gault creatinine",
    "clearance and body weight, while patients on CVVH carry a 10% lower",
    "clearance (categorical CVVH effect) with a body-weight effect only,",
    "plus their individual measured haemofilter clearance supplied as a",
    "covariate. Both volumes carry an inverse power effect of serum albumin",
    "and the central volume also scales with body weight. The creatinine",
    "clearance and weight exponents, Q and the peripheral volume were fixed",
    "from the literature. Inter-individual variability is exponential on",
    "clearance and the central volume, and residual error is proportional.",
    sep = " "
  )
  reference <- paste(
    "Por ED, Akers KS, Chung KK, Livezey JR, Selig DJ.",
    "Population pharmacokinetic modeling and simulations of imipenem in",
    "burn patients with and without continuous venovenous hemofiltration",
    "in the military health system.",
    "J Clin Pharmacol. 2021;61(9):1182-1194. doi:10.1002/jcph.1865.",
    "This model is also catalogued (as study 12) by Zhang P, Zhao Y, Zhu J,",
    "Yang Y, Liang G, Wang X, Yu Z. Population pharmacokinetics of imipenem",
    "in different populations for individualized dosing: a systematic",
    "review. Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055.",
    sep = " "
  )
  vignette <- "Por_2021_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 99.5 kg, the cohort median (range 57.9-150.8 kg; Por",
        "2021 Results, Patient Demographics). Enters the central volume",
        "with exponent 0.74, the no-CVVH clearance branch with exponent",
        "0.33 (both fixed from Bhagunde 2019, Table 2 footnote c), and the",
        "CVVH clearance branch with exponent 0.75 (fixed; Table 2).",
        "Mean weight 105.06 +/- 28.66 kg without CVVH and 89.6 +/- 22.38 kg",
        "with CVVH (Table 1)."
      ),
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "UNIT CONVERSION. Por 2021 models albumin in g/dL (Table 1 'Albumin",
        "(g/dL)'; reference 2.7 g/dL, the cohort median, range 1.5-3.5",
        "g/dL; Results, Patient Demographics). The canonical ALB column is",
        "SI g/L per inst/references/covariate-columns.md, so model()",
        "converts with alb_gdL <- ALB * 0.1 before forming the ratio.",
        "Enters Vc with exponent -1.17 and Vp with exponent -3.68 (Table",
        "2). The Vp exponent is steep: halving albumin multiplies Vp by",
        "2^3.68 = 12.8. The paper uses albumin as a surrogate for burn",
        "severity (Figure 1: albumin = 3.68 - 0.021 * TBSA%), and the",
        "observed albumin range was only 1.5-3.5 g/dL, so predictions for",
        "healthy albumin (4 g/dL, Table 3) are extrapolations."
      ),
      source_name = "ALBUM"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance by the Cockcroft-Gault equation (Por 2021",
        "Methods, Covariate Model), raw mL/min, not normalised to body",
        "surface area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 145.83 mL/min, the median in the no-CVVH patients only",
        "(range 88.08-253.95 mL/min; Results, Patient Demographics). Enters",
        "the NO-CVVH clearance branch only, as (CRCL/145.83)^0.46 with the",
        "exponent fixed from Bhagunde 2019 (Table 2 footnote c, Equation",
        "10). It is multiplied by (1 - RRT_CRRT_STATUS) in model(), so its",
        "value is ignored for a patient on CVVH -- but it must still be a",
        "finite positive number there, because a missing value would",
        "propagate through the product."
      ),
      source_name = "CRCL"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous venovenous haemofiltration during imipenem therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not receiving CVVH",
      notes = paste(
        "Selects between the two published clearance equations (Equations",
        "10 and 11). It also carries the categorical CVVH effect of Table",
        "2, -0.1 in the form tvCL * (1 + theta * CVVH) of Equation 7, which",
        "is why the CVVH-branch intercept printed in Equation 11 is 13.78 =",
        "15.31 * 0.9 L/h. The effect was not statistically significant",
        "(RSE 130.86%, 95% CI -0.36 to 0.16) and was retained on",
        "physiological grounds (Discussion; Table S1 run 18). CVVH is a",
        "continuous renal replacement modality, so the canonical",
        "RRT_CRRT_STATUS column applies. It does NOT gate QEFF (see there)."
      ),
      source_name = "CVVH"
    ),
    QEFF = list(
      description = paste(
        "Individual imipenem clearance by the CVVH haemofilter (L/h),",
        "CL_CVVH = Qf * Sc * CF (Por 2021 Equations 1-3): ultrafiltrate",
        "flow rate times the measured sieving coefficient",
        "Sc = C_filter / ((C_pre + C_post) / 2), times the prefilter",
        "replacement-fluid correction CF = Qb / (Qb + Qrep)."
      ),
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "A per-patient data column, not an estimated parameter: it is added",
        "to the eta-bearing body clearance after the random effect",
        "(Equations 5 and 11). Mean 1.56 +/- 0.7 L/h in the 12 CVVH",
        "patients, with sieving coefficient 0.67 +/- 0.33, correction factor",
        "0.86 +/- 0.04 and effluent flow 30.13 +/- 6.45 mL/kg/h (Table 1).",
        "The paper's simulations use 3 and 5 L/h (Figure 5) and 3.27 L/h",
        "(Figure 3, the Boucher 2016 cohort). NOT gated by RRT_CRRT_STATUS:",
        "Equation 5 adds CL_CVVH unconditionally, so QEFF must be 0 for a",
        "patient not on CVVH. Leaving it ungated also expresses the paper's",
        "Figure 5 scenario, in which a patient with preserved renal function",
        "is placed on CVVH: RRT_CRRT_STATUS = 0 (the CrCl-driven branch) with",
        "QEFF = 3 or 5 L/h reproduces that figure, whereas the CVVH branch,",
        "which has no CrCl term, cannot separate its NRF and ARC curves."
      ),
      source_name = "CLCVVH"
    )
  )

  # Screened in the internal covariate search and not retained (Por 2021
  # Methods, Covariate Model; Table S1 runs 5-16). Lean body weight on CL
  # was significant (P = .01) but weight was preferred for comparability
  # with the literature (Table S1 footnote 3).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained. Mean 51.09 +/- 19.03 years without CVVH, 55 +/- 19.99 with CVVH (Table 1)."
    ),
    LBW = list(
      description = "Lean body weight by the Janmahasatian formula",
      units = "kg",
      type = "continuous",
      notes = "Significant on CL alone (Table S1 run 13, P = .01) and screened on Vc (run 6); not retained, total weight preferred."
    ),
    URINE_VOL_24H = list(
      description = "Urine output",
      units = "mL/24h",
      type = "continuous",
      notes = "Screened on CL (Table S1 run 14), not retained. Mean 3144 mL without CVVH, 949 mL with CVVH (Table 1)."
    ),
    TBSA = list(
      description = "Total burned body surface area, with total second- and third-degree burn areas screened separately",
      units = "%",
      type = "continuous",
      notes = "Screened on Vc, Vp and CL (Table S1 runs 7, 8, 15, 16), not retained; albumin was preferred as the physiological surrogate (Figure 1 relates the two). Documentation only, never referenced in model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    age_mean = "51.09 +/- 19.03 years without CVVH; 55 +/- 19.99 years with CVVH (mean +/- SD)",
    weight_mean = "105.06 +/- 28.66 kg without CVVH; 89.6 +/- 22.38 kg with CVVH (mean +/- SD)",
    weight_range = "57.9-150.8 kg (median 99.5 kg)",
    sex_female_pct = 26.1,
    race_ethnicity = NULL,
    disease_state = paste(
      "Adult patients with severe burns (mean total burn surface area 40%",
      "without and 45% with CVVH) at the US Army Institute of Surgical",
      "Research Burn Center, 12 of them on continuous venovenous",
      "haemofiltration."
    ),
    dose_range = paste(
      "500 mg imipenem every 6 h infused over 30 min or 1 h in 21 patients;",
      "one patient 1000 mg every 6 h over 1 h and one 250 mg every 6 h over",
      "30 min. Sampled at steady state."
    ),
    regions = "United States of America",
    n_concentrations = 81L,
    notes = paste(
      "23 patients enrolled (11 without and 12 with CVVH; 1 woman and 10",
      "men without, 5 women and 7 men with); one had no post-dose sample",
      "and was excluded, leaving 22 in the analysis. 81 prefilter plasma concentrations: a",
      "pre-dose trough plus up to four samples 0.5-8 h post-dose, assayed",
      "by HPLC-UV (linear 0.5-25 ug/mL). Median albumin 2.7 g/dL (range",
      "1.5-3.5) and median CrCl 145.83 mL/min in the no-CVVH group (range",
      "88.08-253.95). Mean CVVH clearance 1.56 +/- 0.7 L/h (Table 1).",
      "Fitted by FOCEI in Pumas 1.0.5 (Table S1 run 21)."
    )
  )

  ini({
    # ===== Structural PK -- Por 2021 Table 2 (final model) and Equations
    # 10-13. Reference subject: 99.5 kg, albumin 2.7 g/dL (27 g/L),
    # Cockcroft-Gault CrCl 145.83 mL/min, not on CVVH. =====
    lcl <- log(15.31); label("Body clearance at the reference subject, not on CVVH (L/h)") # Table 2 'CL (L/h)' 15.31 (RSE 13%); Equation 10 intercept
    lvc <- log(32.67); label("Central volume of distribution at the reference subject (L)") # Table 2 'Vc (L)' 32.67 (RSE 27.18%); Equation 12 intercept
    lq <- fixed(log(11)); label("Intercompartmental clearance (L/h)") # Table 2 'Q (L/h)' 11 fixed (footnote a, literature)
    lvp <- fixed(log(41.23)); label("Peripheral volume of distribution at albumin 2.7 g/dL (L)") # Table 2 'Vp (L)' 41.23 fixed (footnote b, literature); Equation 13 intercept

    # ===== Covariate effects on CL -- Table 2 =====
    e_rrt_crrt_status_cl <- -0.1; label("Fractional change in body clearance on CVVH (unitless)") # Table 2 'CVVH (categorical)' -0.1 (RSE 130.86%); Equation 7 form; 15.31 * 0.9 = 13.78, the Equation 11 intercept
    e_crcl_cl <- fixed(0.46); label("Power exponent on CRCL/145.83 for CL, no-CVVH branch (unitless)") # Table 2 'CrCL (power)' 0.46 fixed (footnote c, Bhagunde 2019); Equation 10
    e_wt_cl <- fixed(0.33); label("Power exponent on WT/99.5 for CL, no-CVVH branch (unitless)") # Table 2 'Weight no CVVH (power)' 0.33 fixed (footnote c); Equation 10
    e_wt_cl_cvvh <- fixed(0.75); label("Power exponent on WT/99.5 for CL, CVVH branch (unitless)") # Table 2 'Weight CVVH (power)' 0.75 fixed; Equation 11

    # ===== Covariate effects on volumes -- Table 2 =====
    e_wt_vc <- fixed(0.74); label("Power exponent on WT/99.5 for Vc (unitless)") # Table 2 'Covariates on Vc, Weight (power)' 0.74 fixed (footnote c); Equation 12
    e_alb_vc <- -1.17; label("Power exponent on albumin/2.7 g/dL for Vc (unitless)") # Table 2 'Covariates on Vc, Albumin (power)' -1.17 (RSE 42.84%); Equation 12
    e_alb_vp <- -3.68; label("Power exponent on albumin/2.7 g/dL for Vp (unitless)") # Table 2 'Covariates on Vp, Albumin (power)' -3.68 (RSE 17%); Equation 13

    # ===== Inter-individual variability -- Table 2 prints the VARIANCES
    # omega^2 directly (eta ~ N(0, omega^2), Equation 4). No IIV on Q or
    # Vp. Pearson correlation between the two etas was 0.05, so no
    # covariance is estimated. =====
    etalcl ~ 0.093 # Table 2 'omega2 CL' 0.093 (RSE 28.18%); eta-shrinkage 8.6%
    etalvc ~ 0.13 # Table 2 'omega2 Vc' 0.13 (RSE 36.45%); eta-shrinkage 45.41%

    # ===== Residual error =====
    propSd <- 0.3; label("Proportional residual error (fraction)") # Table 2 'Proportional error' 0.3 (RSE 22%), read as an SD; see the vignette's Assumptions section
  })

  model({
    # The albumin effects were fitted in g/dL; the canonical ALB column is
    # SI g/L.
    alb_gdL <- ALB * 0.1

    # ----- Body clearance (Equations 10 and 11) -----
    # No CVVH: 15.31 * (CRCL/145.83)^0.46 * (WT/99.5)^0.33 * exp(eta)
    # CVVH:    15.31 * (1 - 0.1) * (WT/99.5)^0.75 * exp(eta) + CL_CVVH
    cl_body <- exp(lcl + etalcl) * (1 + e_rrt_crrt_status_cl * RRT_CRRT_STATUS) *
      ((1 - RRT_CRRT_STATUS) * (CRCL / 145.83)^e_crcl_cl * (WT / 99.5)^e_wt_cl +
        RRT_CRRT_STATUS * (WT / 99.5)^e_wt_cl_cvvh)
    # Haemofilter clearance is a per-patient data column (0 off CVVH),
    # added after the random effect (Equation 5).
    cl_crrt <- QEFF
    # 'cl' must be the TOTAL clearance: rxode2 can solve analytically from
    # a cl/vc pair.
    cl <- cl_body + cl_crrt

    vc <- exp(lvc + etalvc) * (WT / 99.5)^e_wt_vc * (alb_gdL / 2.7)^e_alb_vc
    vp <- exp(lvp) * (alb_gdL / 2.7)^e_alb_vp
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Imipenem-cilastatin IV infusion into central; infusion duration comes
    # from the event table.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
