Mo_2018_olaratumab <- function() {
  description <- "Two-compartment population PK model with linear clearance for olaratumab in patients with advanced or metastatic cancer (Mo 2018)"
  reference <- "Mo G, Baldwin JR, Luffer-Atlas D, et al. Population Pharmacokinetic Modeling of Olaratumab, an Anti-PDGFRα Human Monoclonal Antibody, in Patients with Advanced and/or Metastatic Cancer. Clin Pharmacokinet. 2018;57(3):355-365. doi:10.1007/s40262-017-0562-0"
  vignette <- "Mo_2018_olaratumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed); Mo 2018 states covariates were held at their values at the time of initial assessment. Power scaling on CL and V1 with reference 79.7 kg (population median WTE per Mo 2018 Table 2). Renamed from source column WTE to the canonical WT per covariate-columns.md.",
      source_name        = "WTE"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size (sum of diameters of target lesions per RECIST v1.1)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear deviation from population median 86.5 mm on CL (Mo 2018 Table 3 footnote a; Table 2 median). Patients without measurable tumor data (e.g., GBM substudy) are represented by the median value so the deviation term is zero. Renamed from source column TUMR to the canonical TUMSZ per covariate-columns.md.",
      source_name        = "TUMR"
    )
  )

  population <- list(
    n_subjects     = 171,
    n_studies      = 4,
    age_range      = "22-82 years",
    age_median     = "57 years",
    weight_range   = "37.3-151 kg",
    weight_median  = "79.7 kg",
    sex_female_pct = 51,
    race_ethnicity = c(
      Caucasian = 86,
      `African descent` = 8.8,
      `Asian (east/southeast)` = 1.2,
      `Asian (western)` = 0.6,
      `Native Hawaiian/Other Pacific Islander` = 0.6,
      Other = 2.9
    ),
    disease_state  = "Advanced or metastatic cancer across four phase II studies: soft tissue sarcoma (STS, n=95), nonsmall cell lung cancer (NSCLC, n=50), gastrointestinal stromal tumor (GIST, n=19), glioblastoma multiforme (GBM, n=7).",
    dose_range     = "15 mg/kg IV over 60 min on days 1 and 8 of a 21-day cycle (STS, NSCLC; monotherapy or combined with doxorubicin or paclitaxel/carboplatin) or 20 mg/kg IV over 60-90 min every 14 days (GBM, GIST; monotherapy).",
    regions        = "Global (four phase II studies; regions not tabulated by the publication).",
    tumor_size_range = "12-571 mm (median 86.5 mm, n=164 with RECIST data; no tumor size measured in the 7 GBM patients)",
    notes          = "Baseline demographics per Mo 2018 Table 2. 1501 serum concentration observations from 171 patients analyzed. Treatment combination: monotherapy 31%, doxorubicin 43%, paclitaxel/carboplatin 26%. Treatment-emergent ADA incidence 5% with no effect on olaratumab PK."
  )

  ini({
    # Structural parameters — typical values (final model, Mo 2018 Table 3).
    # Reference covariate values: WT 79.7 kg, TUMSZ 86.5 mm (population medians, Table 2).
    lcl <- log(0.0233); label("Clearance at reference covariates (CL, L/h)")              # Mo 2018 Table 3, Final PK model
    lvc <- log(4.16);   label("Central volume of distribution at reference (V1, L)")      # Mo 2018 Table 3, Final PK model
    lvp <- log(3.58);   label("Peripheral volume of distribution (V2, L)")                # Mo 2018 Table 3, Final PK model
    lq  <- log(0.0315); label("Intercompartmental clearance (Q, L/h)")                    # Mo 2018 Table 3, Final PK model

    # Covariate effects (Mo 2018 Table 3 footnotes a and b)
    e_wt_cl    <- 0.431;   label("Power exponent of WT on CL (unitless)")                 # Mo 2018 Table 3, WTE_CL
    e_wt_vc    <- 0.610;   label("Power exponent of WT on V1 (unitless)")                 # Mo 2018 Table 3, WTE_V1
    e_tumsz_cl <- 0.00158; label("Linear coefficient of TUMSZ deviation on CL (per mm)")  # Mo 2018 Table 3, TUMR_CL

    # IIV on log-scale: omega^2 = log(CV^2 + 1).
    # Mo 2018 reports "no significant correlation between the IPV of V1 and CL" -> diagonal omega.
    etalcl ~ 0.1052  # 33.3% CV on CL (Mo 2018 Table 3 Final PK model, omega^2 = log(0.333^2 + 1))
    etalvc ~ 0.0240  # 15.6% CV on V1 (Mo 2018 Table 3 Final PK model, omega^2 = log(0.156^2 + 1))

    # Combined additive + proportional residual error (Mo 2018 Table 3 Final PK model)
    addSd  <- 10.1;  label("Additive residual error (ug/mL)")         # Mo 2018 Table 3, Final PK model
    propSd <- 0.225; label("Proportional residual error (fraction)")  # Mo 2018 Table 3, Final PK model
  })
  model({
    # Individual PK parameters (Mo 2018 Table 3 footnotes):
    #   CL_ind = CL * (WT / median(WT))^e_wt_cl * (1 + e_tumsz_cl * (TUMSZ - median(TUMSZ)))
    #   V1_ind = V1 * (WT / median(WT))^e_wt_vc
    # Reference (median) values: WT = 79.7 kg, TUMSZ = 86.5 mm (Mo 2018 Table 2).
    cl <- exp(lcl + etalcl) * (WT / 79.7)^e_wt_cl * (1 + e_tumsz_cl * (TUMSZ - 86.5))
    vc <- exp(lvc + etalvc) * (WT / 79.7)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
