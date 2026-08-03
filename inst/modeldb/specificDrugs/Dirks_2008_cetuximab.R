Dirks_2008_cetuximab <- function() {
  description <- "Two-compartment population PK model for intravenous cetuximab (anti-EGFR chimeric IgG1) in adults with recurrent and/or metastatic squamous cell carcinoma of the head and neck (SCCHN), with Michaelis-Menten (target-mediated) elimination from the central compartment. Ideal body weight and white blood cell count are linear-deviation covariates on Vmax; total body weight is a linear-deviation covariate on V1 (Dirks 2008 J Clin Pharmacol; Chapter 3 of the Dirks 2010 UTHSC PhD dissertation)."
  reference <- "Dirks NL, Nolting A, Kovar A, Meibohm B. Population pharmacokinetics of cetuximab in patients with squamous cell carcinoma of the head and neck. J Clin Pharmacol. 2008;48(3):267-278. doi:10.1177/0091270007313393. See also Chapter 3 (pp. 45-64) and Appendix A (NONMEM control stream, pp. 141-149) of Dirks NL, Population Pharmacokinetics of Therapeutic Monoclonal Antibodies: Examples and Estimation Method Performance Differences, PhD dissertation, University of Tennessee Health Science Center, May 2010, ETD paper 63, doi:10.21007/etd.cghs.2010.0072."
  vignette <- "Dirks_2008_cetuximab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "cetuximab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "cetuximab", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    IBW = list(
      description        = "Ideal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only. Enters as an additive linear-deviation term on Vmax: Vmax = TVVmax * (1 + e_ibw_vmax * (IBW - 64)) with reference 64 kg (Dirks 2008 Table 3-2 in thesis Chapter 3, p. 53). The paper does not identify which Devine-family IBW variant was used by the source dataset; if only WT and HT are available the user must compute IBW externally with a chosen formula and record that choice.",
      source_name        = "IBW"
    ),
    WBC = list(
      description        = "Total white blood cell count at baseline",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only. Enters as an additive linear-deviation term on Vmax: Vmax = TVVmax * (1 + e_wbc_vmax * (WBC - 6.8)) with reference 6.8 x 10^9/L (Dirks 2008 Table 3-2 in thesis Chapter 3, p. 53). Interpreted mechanistically as a surrogate for EGFR-bearing leukocyte burden (lymphocytes, monocytes, macrophages, neutrophils all express EGFR); a 10^9/L rise in WBC increases Vmax by 2.2 percent.",
      source_name        = "WBC"
    ),
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only. Enters as an additive linear-deviation term on V1: V1 = TVV1 * (1 + e_wt_v1 * (WT - 60)) with reference 60 kg (Dirks 2008 Table 3-2 in thesis Chapter 3, p. 53). WGT in the NONMEM control stream (Appendix A) maps to WT here.",
      source_name        = "WGT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 143L,
    n_studies      = 2L,
    age_range      = "23-77 years",
    age_median     = "56 years",
    weight_range   = "34-113 kg",
    weight_median  = "60 kg",
    sex_female_pct = 16.1,
    race_ethnicity = c(Caucasian = 91.6, Unknown = 8.4),
    disease_state  = "Recurrent and/or metastatic squamous cell carcinoma of the head and neck (SCCHN). 138 of 143 patients had EGFR-expressing tumours by immunohistochemistry; human anti-chimeric antibodies were not detected in either study.",
    dose_range     = "Cetuximab IV: 400 mg/m^2 loading over 2 hours followed by 250 mg/m^2 weekly maintenance over 1 hour (currently approved regimen); duration of therapy 1-54 weeks in the modelled cohort (median 6 weeks). Study A n=47; Study B n=96. Concomitant 5-fluorouracil in 47/143; concomitant platinum in 89/143.",
    regions        = "Two phase I/II studies of cetuximab in SCCHN (institutions not specified in the thesis).",
    notes          = "Median patient characteristics from Table 3-1 (p. 51). Reference covariates for the typical-value equations (Table 3-2, p. 53): IBW 64 kg, WBC 6.8 x 10^9/L, WT (WGT) 60 kg. 912 concentrations available for model building (530 from study A, 382 from study B); Study A patients contributed a median of 13 samples each and Study B patients a median of 4 samples each. IBW median in Table 3-1 was 64.2 kg; the model uses 64 kg as the centring value."
  )

  ini({
    # Structural typical-value PK parameters at the reference covariates
    # (Dirks 2008 Table 3-2, thesis Chapter 3, p. 53; reference IBW 64 kg,
    # WBC 6.8 x 10^9/L, WT 60 kg).
    lvmax <- log(4.38);  label("Michaelis-Menten Vmax at reference covariates (mg/h)") # Dirks 2008 Table 3-2: theta1 = 4.38 mg/hr, 90% CI 3.40-6.64
    lkm   <- log(74);    label("Michaelis-Menten Km (ug/mL)")                              # Dirks 2008 Table 3-2: theta4 = 74 ug/mL, 90% CI 38.2-163.3
    lvc   <- log(2.83);  label("Central volume of distribution V1 at reference WT (L)")   # Dirks 2008 Table 3-2: theta5 = 2.83 L, 90% CI 2.69-2.96
    lvp   <- log(2.43);  label("Peripheral volume of distribution V2 (L)")                 # Dirks 2008 Table 3-2: theta7 = 2.43 L, 90% CI 1.95-2.85
    lq    <- log(0.103); label("Intercompartmental clearance Q (L/h)")                  # Dirks 2008 Table 3-2: theta8 = 0.103 L/hr, 90% CI 0.062-0.191

    # Covariate effects. Additive linear-deviation form on continuous
    # covariates as parameterised in Dirks 2008 Table 3-2:
    #   Vmax = theta1 * (1 + theta2 * (IBW - 64) + theta3 * (WBC - 6.8))
    #   V1   = theta5 * (1 + theta6 * (WT  - 60))
    e_ibw_vmax <- 0.0108; label("Additive linear-deviation coefficient of (IBW - 64) on Vmax (per kg)")            # Dirks 2008 Table 3-2: theta2 = 0.0108, 90% CI 0.0077-0.0139
    e_wbc_vmax <- 0.0216; label("Additive linear-deviation coefficient of (WBC - 6.8) on Vmax (per 10^9 cells/L)") # Dirks 2008 Table 3-2: theta3 = 0.0216, 90% CI 0.0169-0.0296
    e_wt_v1    <- 0.0083; label("Additive linear-deviation coefficient of (WT - 60) on V1 (per kg)")               # Dirks 2008 Table 3-2: theta6 = 0.0083, 90% CI 0.0057-0.0115

    # Inter-individual variability (Dirks 2008 Table 3-2, diagonal Omega).
    # Paper reports BSV as %CV; convert to log-normal variance via
    # omega^2 = log(CV^2 + 1).
    # Vmax  15.4% -> log(0.154^2 + 1) = 0.02348
    # V1    18.6% -> log(0.186^2 + 1) = 0.03410
    # V2    56.4% -> log(0.564^2 + 1) = 0.27443
    # Q     97.2% -> log(0.972^2 + 1) = 0.63959
    # BSV on Km was intentionally not estimated in the final model
    # (Chapter 3 Structural Model, p. 50: "K_m parameter was treated as
    # a fixed value across all subjects, which resulted in an OFV increase
    # of only 0.3 points"). etalkm is held at fixed(0) to preserve the
    # published diagonal structure without a spurious Km random effect.
    etalvmax ~ 0.02348    # Dirks 2008 Table 3-2: omega(Vmax) 15.4% CV, 90% CI 12.0-19.1
    etalkm   ~ fixed(0)   # Dirks 2008 Table 3-2: Km BSV not estimated (see p. 50 Structural Model)
    etalvc   ~ 0.03410    # Dirks 2008 Table 3-2: omega(V1)   18.6% CV, 90% CI 12.5-22.2
    etalvp   ~ 0.27443    # Dirks 2008 Table 3-2: omega(V2)   56.4% CV, 90% CI 18.0-72.8
    etalq    ~ 0.63959    # Dirks 2008 Table 3-2: omega(Q)    97.2% CV, 90% CI 40.2-133

    # Residual error. Dirks 2008 fit a log-additive residual (equivalent
    # to a proportional error in linear space) separately for each study
    # (Study A: 14.6% CV; Study B: 21.2% CV; Table 3-2 and p. 61). Only
    # one residual is exported here because the library model has no
    # STDY covariate; the Study A value (14.6%) is used because study A
    # contributed the denser sampling (median 13 observations per patient
    # vs 4 per patient in Study B, so it is more informative for the
    # per-observation residual). Study B is documented in the vignette
    # source-trace and in Assumptions and deviations.
    propSd <- 0.146; label("Proportional residual error (fraction)") # Dirks 2008 Table 3-2: Study A 14.6% CV; Study B 21.2% CV not encoded
  })
  model({
    # Individual PK parameters with covariate effects (Dirks 2008 final
    # covariate model; Table 3-2 equations, thesis Chapter 3, p. 53):
    #   Vmax = 4.38  * (1 + 0.0108 * (IBW - 64) + 0.0216 * (WBC - 6.8))
    #   Km   = 74
    #   V1   = 2.83  * (1 + 0.0083 * (WT  - 60))
    #   V2   = 2.43
    #   Q    = 0.103
    vmax <- exp(lvmax + etalvmax) *
              (1 + e_ibw_vmax * (IBW - 64) + e_wbc_vmax * (WBC - 6.8))
    km   <- exp(lkm   + etalkm)
    vc   <- exp(lvc   + etalvc) *
              (1 + e_wt_v1 * (WT - 60))
    vp   <- exp(lvp   + etalvp)
    q    <- exp(lq    + etalq)

    # Two-compartment IV-input PK with Michaelis-Menten elimination from
    # the central compartment (Dirks 2008 Chapter 3 Appendix A NM-TRAN,
    # thesis p. 141-149; ADVAN6 TRANS1). A parallel first-order
    # elimination pathway was tested but produced no OFV improvement and
    # was dropped in favour of the more parsimonious purely nonlinear
    # elimination model (p. 50 Structural Model). Vmax is in mg/hour and
    # depends on Cc = central / vc (concentration in ug/mL = mg/L because
    # doses are in mg and volumes in L).
    Cc <- central / vc

    d/dt(central)     <- -vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ prop(propSd)
  })
}
