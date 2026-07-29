Akbar_2025_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order elimination for intravenous voriconazole in adult and pediatric Pakistani cancer patients receiving therapeutic drug monitoring (Akbar 2025); creatinine clearance and primary cancer diagnosis are covariates on clearance"
  reference <- "Akbar Z, Usman M, Aamir M, Saleem Z, Khan MR, Alamri A, Alharbi MS, Osman GEM. Population pharmacokinetic analysis of intravenous voriconazole in cancer patients. PLoS ONE. 2025;20(3):e0318883. doi:10.1371/journal.pone.0318883"
  vignette <- "Akbar_2025_voriconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centered-linear effect on CL: CL_typical *= (1 + 0.0077 * (CRCL - 120)). Reference 120 mL/min (Akbar 2025 Eq 7). Akbar 2025 does not state which formula was used to derive CRCL or whether the value is BSA-normalized; the inclusion of pediatric subjects (down to 3 years / 5.3 kg) suggests the Schwartz formula was likely used for the pediatric stratum (which yields BSA-normalized mL/min/1.73 m^2) while Cockcroft-Gault was likely used for adults (raw mL/min). The slope 0.0077 was estimated against whatever scale was used in the source data; downstream simulation should supply CRCL on the same scale. The canonical CRCL register entry is documented as BSA-normalized but is reused here because the source paper simply names the column 'creatinine clearance' without further qualification.",
      source_name        = "CLCR"
    ),
    TUMTP_LYMPH = list(
      description        = "Lymphoma cancer-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-lymphoma; when paired with all other Akbar 2025 TUMTP_* indicators = 0, the implicit reference is leukemia)",
      notes              = "Additive-fractional +1.91% effect on CL (Eq 2 of Akbar 2025). Source category 'Lymphoma' (n = 14, 15.9% of cohort). The 95% CI on the coefficient spans zero (-0.382 to 0.777 per Table 2 bootstrap), so the lymphoma effect on CL is not statistically distinguishable from leukemia in this dataset.",
      source_name        = "DISEASE == 'Lymphoma'"
    ),
    TUMTP_SARC = list(
      description        = "Sarcoma cancer-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-sarcoma; when paired with all other Akbar 2025 TUMTP_* indicators = 0, the implicit reference is leukemia)",
      notes              = "Additive-fractional +18.5% effect on CL (Eq 3 of Akbar 2025). Source category 'Sarcoma' (n = 9, 10.2% of cohort). The 95% CI on the coefficient spans zero (-0.403 to 2.01 per Table 2 bootstrap row labelled CL-DISEASE 3, where Table 2 footnotes 'd' and 'e' are mislabelled relative to the equation text -- see vignette Errata).",
      source_name        = "DISEASE == 'Sarcoma'"
    ),
    TUMTP_BREAST = list(
      description        = "Breast-cancer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-breast-cancer; when paired with all other Akbar 2025 TUMTP_* indicators = 0, the implicit reference is leukemia)",
      notes              = "Additive-fractional -81.7% effect on CL (Eq 4 of Akbar 2025); breast cancer subjects had ~5x lower CL than leukemia at the same CRCL. Source category 'Breast cancer' (n = 8, 9.1% of cohort). Effect is statistically significant (95% CI -0.919 to -0.593 per Table 2 bootstrap row labelled CL-DISEASE 4) and is consistent with Fig 1 of Akbar 2025 which shows breast cancer subjects with the lowest median CL (~2 L/h vs the leukemia reference 6.17 L/h).",
      source_name        = "DISEASE == 'Breast cancer'"
    ),
    TUMTP_MYELO = list(
      description        = "Multiple myeloma cancer-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-myeloma; when paired with all other Akbar 2025 TUMTP_* indicators = 0, the implicit reference is leukemia)",
      notes              = "Additive-fractional -2.33% effect on CL (Eq 5 of Akbar 2025). Source category 'Myeloma' (n = 4, 4.5% of cohort). The 95% CI on the coefficient spans zero (-0.365 to 0.553 per Table 2 bootstrap), so the myeloma effect on CL is not statistically distinguishable from leukemia.",
      source_name        = "DISEASE == 'Myeloma'"
    ),
    TUMTP_GLIO = list(
      description        = "Glioma cancer-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-glioma; when paired with all other Akbar 2025 TUMTP_* indicators = 0, the implicit reference is leukemia)",
      notes              = "Additive-fractional +8.81% effect on CL (Eq 6 of Akbar 2025). Source category 'Glioma' (n = 3, 3.4% of cohort). The 95% CI on the coefficient spans zero (-0.233 to 0.830 per Table 2 bootstrap), so the glioma effect on CL is not statistically distinguishable from leukemia in this small cohort.",
      source_name        = "DISEASE == 'Glioma'"
    )
  )

  population <- list(
    n_subjects     = 88L,
    n_studies      = 1L,
    n_observations = 88L,
    age_range      = "3-90 years",
    age_median     = "35 years (38 years per Abstract; Results section reports median 38 with range 3-90; Table 1 reports median 35 with range 3-67 -- internal inconsistency in the source paper)",
    weight_range   = "5.3-95.6 kg",
    weight_median  = "30.85 kg",
    sex_female_pct = 46.6,
    age_group      = c(adult_pct = 53.4, pediatric_pct = 46.6),
    cancer_type    = c(Leukemia_pct = 56.8, Lymphoma_pct = 15.9, Sarcoma_pct = 10.2, BreastCancer_pct = 9.1, Myeloma_pct = 4.5, Glioma_pct = 3.4),
    fungal_infection = c(Suspected_pct = 54.5, Aspergillosis_pct = 35.2, FungalPneumonia_pct = 10.2),
    disease_state  = "Adult and pediatric cancer patients (mixed solid and hematologic malignancies) receiving intravenous voriconazole for systemic fungal infections; majority empirical therapy for suspected fungal infection in immunocompromised cancer patients including febrile neutropenia.",
    dose_range     = "Intravenous voriconazole loading 6 mg/kg q12h for 24 h, then maintenance 4 mg/kg q12h. Administered doses 56-400 mg per dose in the dataset. Trough samples drawn on day 5 before the morning dose per hospital protocol.",
    regions        = "Single center: Shaukat Khanum Memorial Cancer Hospital and Research Centre, Lahore, Pakistan.",
    notes          = "Retrospective therapeutic drug monitoring (TDM) data collected 1 January 2023 - 31 December 2023 from electronic medical records. 488 admitted cancer patients with systemic fungal infections were screened, 112 had voriconazole TDM data, and 88 with complete information were retained. Trough concentrations spanned 0.10-21.0 ug/mL. Bioanalytical method: homogeneous enzyme immunoassay on Siemens Atellica CH-930. Single observation per subject (sparse TDM design). Baseline demographics per Akbar 2025 Table 1; final-model parameter estimates per Akbar 2025 Table 2."
  )

  ini({
    # Structural parameters -- typical-value reference is a leukemia patient
    # at CRCL = 120 mL/min (Akbar 2025 Eq 1 + Eq 7 reference). Voriconazole
    # is administered intravenously, so there is no absorption stage.
    lcl <- log(6.17); label("Clearance for the leukemia reference at CRCL 120 mL/min (L/h)")     # Akbar 2025 Table 2 (final estimate) and Eq 1
    lvc <- log(55.9); label("Volume of distribution (L)")                                        # Akbar 2025 Table 2 (final estimate)

    # Covariate effects on CL.
    # CRCL: centered-linear, 0.0077 per (mL/min - 120).
    e_crcl_cl  <-  0.0077;  label("Linear coefficient for centered CRCL on CL (per mL/min)")     # Akbar 2025 Table 2 (CRCL-CL) and Eq 7
    # Cancer-type indicators: additive fractional change from leukemia.
    # Values come from the Eq 1-6 typed-out equations in Akbar 2025
    # Results section ("The Equations 1-6 describes the effect of type of
    # cancer on CL"). Note: Akbar 2025 Table 2 footnotes 'd' / 'e' / 'f'
    # mislabel the diseases for the rows CL-DISEASE 3 and CL-DISEASE 4
    # (Table footnotes call them "Breast cancer" and "sarcoma" but the
    # equations and Fig 1 box-plot show them in the opposite order).
    # The numerical values are unambiguous; the equation text is internally
    # consistent with itself and with Fig 1 (breast cancer has the lowest
    # observed CL in Fig 1, matching the -0.817 fractional change), so the
    # equation-text mapping is taken as authoritative here.
    e_lymph_cl <-  0.0191;  label("Additive-fractional effect of lymphoma vs leukemia on CL (unitless)")        # Akbar 2025 Eq 2
    e_sarc_cl  <-  0.185;   label("Additive-fractional effect of sarcoma vs leukemia on CL (unitless)")         # Akbar 2025 Eq 3
    e_bc_cl    <- -0.817;   label("Additive-fractional effect of breast cancer vs leukemia on CL (unitless)")   # Akbar 2025 Eq 4
    e_myelo_cl <- -0.0233;  label("Additive-fractional effect of myeloma vs leukemia on CL (unitless)")         # Akbar 2025 Eq 5
    e_glio_cl  <-  0.0881;  label("Additive-fractional effect of glioma vs leukemia on CL (unitless)")          # Akbar 2025 Eq 6

    # IIV. Akbar 2025 reports IIV on CL only (no IIV on Vd). NONMEM
    # exponential IIV is `CL = TVCL * exp(eta)` with var(eta) = omega^2;
    # the source paper reports IIV-CL = 83.72% as a CV%. Using the standard
    # NONMEM convention CV% ~= sqrt(omega^2) * 100, omega^2 = 0.8372^2 =
    # 0.7009. The exact log-normal variance log(1 + CV^2) = log(1 + 0.7009)
    # = 0.5311 is an alternative interpretation; the approximate form is
    # used here as it is the more common reporting convention in NONMEM
    # popPK papers and aligns with the bootstrap CI (48.41% to 104.86%
    # in Table 2) being symmetric on the CV%-as-omega scale.
    etalcl ~ 0.7009                                           # Akbar 2025 Table 2 (IIV-CL 83.72%); see comment above for omega^2 derivation

    # Residual error. Akbar 2025 reports a proportional error model with
    # final-model proportional error 0.44.
    propSd <- 0.44; label("Proportional residual error (fraction)")  # Akbar 2025 Table 2 (Proportional error 0.44)
  })

  model({
    # Cancer-type effect on CL is encoded as the sum of fractional effects
    # over a one-hot decomposition of the cancer-type categorical (Eqs 1-6
    # of Akbar 2025). Exactly one of the indicator columns is 1 for any
    # given subject (or all are 0 when the subject is in the leukemia
    # reference category).
    disease_factor <- 1 +
      e_lymph_cl * TUMTP_LYMPH +
      e_sarc_cl  * TUMTP_SARC +
      e_bc_cl    * TUMTP_BREAST +
      e_myelo_cl * TUMTP_MYELO +
      e_glio_cl  * TUMTP_GLIO

    # Centered-linear CRCL effect on CL (Akbar 2025 Eq 7).
    crcl_factor <- 1 + e_crcl_cl * (CRCL - 120)

    # Individual CL combines disease_factor, crcl_factor, and exponential IIV.
    cl <- exp(lcl + etalcl) * disease_factor * crcl_factor
    vc <- exp(lvc)

    # 1-compartment IV PK with first-order elimination -- linCmt() resolves
    # to the analytical solution given cl and vc. NONMEM ADVAN1 TRANS2.
    Cc <- linCmt()

    # Concentration units: dose mg / volume L = mg/L = ug/mL, matching the
    # source TDM trough-concentration units.
    Cc ~ prop(propSd)
  })
}
