McLaughlin_2026_axicabtageneCiloleucel_brexucabtageneAutoleucel <- function() {
  description <- "Mechanism-based population PK model of anti-CD19 CAR-T-cell kinetics for axicabtagene ciloleucel and brexucabtagene autoleucel in relapsed/refractory B-cell lymphomas: central and peripheral CAR-T-cell compartments with tumour-burden-driven first-order expansion after a lag, and a latent tumour (SPD) compartment depleted by second-order CAR-T killing (Mc Laughlin 2026)"
  reference <- paste(
    "Mc Laughlin AM, Bergstrand M, Ruiz-Garcia A, Filosto S, Shen R.",
    "Model-Based Analysis with Mechanistic Insights to CAR-T-Cell Therapy",
    "Kinetics: Case Study with Axicabtagene Ciloleucel and Brexucabtagene",
    "Autoleucel. Clin Pharmacol Ther. 2026;120(4):1027-1040.",
    "doi:10.1002/cpt.70395.",
    sep = " "
  )
  vignette <- "McLaughlin_2026_axicabtageneCiloleucel_brexucabtageneAutoleucel"
  units <- list(
    time = "day",
    dosing = "10^6 CAR-T cells",
    concentration = "cells/uL"
  )

  compartmentData <- list(
    central = list(
      analyte = "CAR-T cells (axicabtagene ciloleucel or brexucabtagene autoleucel)",
      units = "10^6 cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "CAR-T cells (axicabtagene ciloleucel or brexucabtagene autoleucel)",
      units = "10^6 cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tumor = list(
      analyte = "latent tumour burden (sum of the product of perpendicular diameters, SPD)",
      units = "mm^2",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUMSZ = list(
      description = "Baseline tumour burden: sum of the product of perpendicular diameters (SPD) of target lesions",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Initial condition of the latent tumour compartment, tumor(0) = TUMSZ; the expansion driver is the ratio tumor / TUMSZ, so the predicted CAR-T kinetics do not depend on the TUMSZ value itself (the killing term is first order in tumour burden). Overall median 2,244 mm^2 (Table 1); Figure 4 uses 2,244 mm^2 for the reference patient.",
      source_name = "SPD (baseline SPD, SPDBL)"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Log-transformed exponential effect centred on the overall median 82.7 kg, i.e. kel * (WT / 82.7)^0.781 (Table 2 footnote a; Methods 'Covariate analysis').",
      source_name = "Baseline body weight"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on the maximum expansion rate centred on the overall median 40 g/L: exp(0.0163 * (ALB - 40)).",
      source_name = "Baseline albumin"
    ),
    CART_TEFF_PCT = list(
      description = "Effector T cells (CD45RA+ CCR7-) as a percentage of viable CD3+ cells in the infused CAR-T product",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on the maximum expansion rate centred on the overall median 24.7%: exp(-0.00619 * (CART_TEFF_PCT - 24.7)).",
      source_name = "Product TEFF (% of viable CD3+ cells)"
    ),
    PRIOR_BENDAMUSTINE = list(
      description = "Prior bendamustine therapy before CAR-T-cell infusion",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior bendamustine)",
      notes = "Fractional-change effect on the maximum expansion rate: (1 - 0.241 * PRIOR_BENDAMUSTINE). 22/406 axi-cel and 36/67 brexu-cel patients (Table 1).",
      source_name = "Prior bendamustine treatment"
    ),
    TUMTP_MCL = list(
      description = "Brexucabtagene autoleucel product type / mantle cell lymphoma disease type composite indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (axicabtagene ciloleucel / non-MCL lymphoma)",
      notes = "Every brexu-cel patient had MCL and no axi-cel patient did, so the paper estimated product type and disease type as one composite covariate. Fractional-change effect on the killing rate constant: (1 - 0.617 * TUMTP_MCL). Set to 1 to simulate brexu-cel in MCL.",
      source_name = "Brexu-cel product type/MCL disease type"
    ),
    STUDY_ZUMA7 = list(
      description = "ZUMA-7 study indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ZUMA-1, ZUMA-2 or ZUMA-5)",
      notes = "Fractional-change effect on the residual-error SD: expSd * (1 + 0.719 * STUDY_ZUMA7). Affects only residual noise, not the typical kinetics.",
      source_name = "Study ZUMA-7"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 473,
    n_studies = 4,
    n_observations = 2785,
    age_range = "21-80 years (median 60)",
    weight_range = "43.1-166 kg (median 82.7)",
    sex_female_pct = 35.5,
    race_ethnicity = c(White = 87.5, Black = 4.02, Asian = 3.17, NativeHawaiian = 0.423, Other = 4.86),
    disease_state = "Relapsed or refractory B-cell lymphoma: large B-cell lymphoma (DLBCL, PMBCL, TFL and other), follicular lymphoma and marginal zone lymphoma (axicabtagene ciloleucel, n = 406), and mantle cell lymphoma (brexucabtagene autoleucel, n = 67).",
    dose_range = "Single 30-min intravenous infusion at a target 2 x 10^6 CAR-T cells/kg after cyclophosphamide/fludarabine lymphodepletion; actual infused CAR+ cells 51.8-201 x 10^6 (median 165 x 10^6).",
    regions = "Multinational (ZUMA-1, ZUMA-2, ZUMA-5, ZUMA-7)",
    studies = "ZUMA-1 (n = 99), ZUMA-2 (n = 67), ZUMA-5 (n = 146), ZUMA-7 (n = 161)",
    notes = "Peripheral-blood CAR-T-cell concentrations from qPCR converted to cells/uL; 14.6% of samples reported as zero were handled with the M3 method at study-specific censoring thresholds (ZUMA-1 0.00744, ZUMA-2 0.00612, ZUMA-5 0.00273273, ZUMA-7 0.0001624 cells/uL). Baseline demographics from Table 1."
  )

  ini({
    lkexp_max <- log(2.70); label("Maximum first-order CAR-T-cell expansion rate constant rho (1/day)") # Table 2 'rho' = 2.70 day^-1 (RSE 4.74%)
    lkel <- log(0.120); label("CAR-T-cell elimination rate constant kel (1/day)") # Table 2 'k el' = 0.120 day^-1 (RSE 8.06%)
    lk12 <- log(0.0465); label("Central-to-peripheral CAR-T-cell distribution rate constant k12 (1/day)") # Table 2 'k 12' = 0.0465 day^-1 (RSE 15.1%)
    lk21 <- log(0.0108); label("Peripheral-to-central CAR-T-cell distribution rate constant k21 (1/day)") # Table 2 'k 21' = 0.0108 day^-1 (RSE 15.4%)
    lkkill <- log(1.33e-6); label("Second-order tumour killing rate constant of CAR-T cells kkill (1/(day*10^6 cells))") # Table 2 'k kill' = 1.33 x 10^-6 (day*10^6 cells)^-1 (RSE 16.9%)
    ltlag <- log(2.47); label("Lag time before CAR-T-cell expansion starts (day)") # Table 2 't lag before expansion' = 2.47 days (RSE 5.24%)
    lvc <- log(3.14e4); label("Scaling factor linking central CAR-T amount to blood concentration (L)") # Table 2 'Scaling factor' = 3.14 x 10^4 L (RSE 12.0%)

    e_wt_kel <- 0.781; label("Exponent of log(WT/82.7) on kel (unitless)") # Table 2 'Baseline body weight on k el' = 0.781 (RSE 27.6%; footnote a log-transformed)
    e_cart_teff_pct_kexp_max <- -0.00619; label("Exponential effect of product %TEFF on rho (per %)") # Table 2 'Product % TEFF on rho' = -0.00619 (RSE 24.2%)
    e_prior_bendamustine_kexp_max <- -0.241; label("Fractional change in rho with prior bendamustine (unitless)") # Table 2 'Prior bendamustine therapy on rho' = -0.241 (RSE 18.7%)
    e_alb_kexp_max <- 0.0163; label("Exponential effect of baseline albumin on rho (per g/L)") # Table 2 'Baseline albumin on rho' = 0.0163 (RSE 26.3%)
    e_tumtp_mcl_kkill <- -0.617; label("Fractional change in kkill for brexu-cel product type / MCL disease type (unitless)") # Table 2 'Brexu-cel product type/MCL disease type on k kill' = -0.617 (RSE 12.3%)
    e_study_zuma7_expSd <- 0.719; label("Fractional change in the residual SD for study ZUMA-7 (unitless)") # Table 2 'Study ZUMA-7 on RUV' = 0.719 (RSE 18.1%)

    # IIV reported as approximate CV = sqrt(omega^2) (PsN convention; the Table 2
    # footnote reports RSEs on the approximate SD scale), so omega^2 = CV^2.
    etalkexp_max ~ 0.069696 # Table 2 'IIV rho' (CV) = 0.264; 0.264^2
    etalkel ~ 0.471969 # Table 2 'IIV k el' (CV) = 0.687; 0.687^2
    etalk12 ~ 1.8225 # Table 2 'IIV k 12' (CV) = 1.35; 1.35^2
    etalk21 ~ 3.8809 # Table 2 'IIV k 21' (CV) = 1.97; 1.97^2
    etalkkill ~ 0.912025 # Table 2 'IIV k kill' (CV) = 0.955; 0.955^2
    etaexpSd ~ 0.297025 # Table 2 'IIV RUV' (CV) = 0.545; 0.545^2

    expSd <- 0.969; label("Additive residual SD on the log scale (log cells/uL)") # Table 2 'RUV' (CV) = 0.969 (RSE 4.94%)
  })
  model({
    # Covariates: continuous effects exponential and centred on the Table 1
    # overall medians (Methods 'Covariate analysis'); body weight entered on
    # the log scale; categorical effects are fractional changes.
    kexp_max <- exp(lkexp_max + etalkexp_max) *
      exp(e_cart_teff_pct_kexp_max * (CART_TEFF_PCT - 24.7)) *
      exp(e_alb_kexp_max * (ALB - 40)) *
      (1 + e_prior_bendamustine_kexp_max * PRIOR_BENDAMUSTINE)
    kel <- exp(lkel + etalkel) * exp(e_wt_kel * (log(WT) - log(82.7)))
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kkill <- exp(lkkill + etalkkill) * (1 + e_tumtp_mcl_kkill * TUMTP_MCL)
    tlag <- exp(ltlag)
    vc <- exp(lvc)

    # Expansion switches on once the time since infusion exceeds tlag
    # (Equation 1, tlag,IND); the dose is the single infusion at time 0.
    tlag_ind <- 0
    if (t > tlag) tlag_ind <- 1

    tumor(0) <- TUMSZ

    d/dt(central) <- kexp_max * central * tlag_ind * tumor / TUMSZ -
      k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(tumor) <- -kkill * central * tumor

    # 10^6 cells / L = 1 cell/uL
    Cc <- central / vc

    # Additive error on log-transformed data with an eta on the residual SD
    expSdi <- expSd * exp(etaexpSd) * (1 + e_study_zuma7_expSd * STUDY_ZUMA7)
    Cc ~ lnorm(expSdi)
  })
}
