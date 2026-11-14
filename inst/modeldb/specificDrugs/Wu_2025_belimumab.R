Wu_2025_belimumab <- function() {
  description <- "Two-compartment population PK model for intravenous and subcutaneous belimumab in healthy adults and adult and paediatric patients with systemic lupus erythematosus from nine studies, including Mainland Chinese participants (Wu 2025). First-order SC absorption with a lag time and logit-scale bioavailability. Fat-free mass scales all disposition parameters; baseline albumin and IgG act on CL; a saturable age term acts on Vc; a Mainland China indicator raises CL by 28% and Vc by 20%. Combined proportional + additive residual error with between-subject variability on the residual magnitude."
  reference <- paste(
    "Wu J, Zhou X, Dimelow R, Marshall S. Inter-regional pharmacokinetics",
    "and exposure-response analyses of belimumab in patients with system",
    "lupus erythematosus. Br J Clin Pharmacol. 2025;91(2):374-385.",
    "doi:10.1111/bcp.16263.",
    "Final PK parameters and equations: Supporting Information Table S2.",
    sep = " "
  )
  vignette <- "Wu_2025_belimumab"
  paper_specific_etas <- c("etaRUV")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  compartmentData <- list(
    depot = list(analyte = "belimumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "belimumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "belimumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Baseline fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column FFMBL. Computed in the source from baseline body weight, BMI and sex with the Janmahasatian 2005 formula (Wu 2025 Methods, reference 16), for adults and children alike. Power effects on CL and Q (shared exponent 0.739) and on Vc and Vp (shared exponent 0.711), reference 40.486 kg (Table S2 equations; the median of the PK population, Table S1).",
      source_name = "FFMBL"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column ALBBL. Power effect on CL with reference 41 g/L (Table S2). The main-text Results print the reference as '41.0 g dL-1'; Table S1 reports albumin in g/L (median 41.0, range 18-55), so the unit is g/L.",
      source_name = "ALBBL"
    ),
    IGG = list(
      description = "Baseline serum immunoglobulin G",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column IgGBL. Power effect on CL with reference 14.7 g/L (Table S2).",
      source_name = "IgGBL"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Saturable effect on Vc, AGE / (AGE + 2.69). The factor is 0.69 at 6 years, 0.84 at 14 years and 0.94 at 40 years. The paediatric study enrolled children aged 6-18 years; no information below 6 years.",
      source_name = "Age"
    ),
    REGION_CHINA = list(
      description = "Mainland China participant indicator (1 = enrolled in Mainland China, 0 = elsewhere)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other regions, including Hong Kong, Taiwan, Japan and Korea)",
      notes = "Source column MLC. Multiplicative factors 1.28 on CL and 1.20 on Vc (Table S2). Participants from Hong Kong and Taiwan are 0. The paper attributes the effect most plausibly to an inter-laboratory bioanalytical difference (Mainland Chinese samples were assayed in a Shanghai laboratory) rather than to ethnicity (Discussion).",
      source_name = "MLC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2002L,
    n_studies = 9L,
    n_observations = 12061L,
    age_range = "6-77 years (median 36.0)",
    weight_range = "17.0-165.4 kg (median 63.5)",
    sex_female_pct = 89.8,
    race_ethnicity = "Region strata (Table S1): Caucasian 919, Mainland Chinese 98, other North East Asian (Hong Kong, Taiwan, Japan, Korea) 219, others 713, paediatric 53 (27 Caucasian, 8 Asian, 3 Black, 15 Native American/Hawaiian).",
    disease_state = "Systemic lupus erythematosus (1832, 91.5%) and healthy volunteers (170, 8.5%).",
    dose_range = "IV 1 or 10 mg/kg every 4 weeks (after loading doses on days 0, 14 and 28) and single IV doses of 200-240 mg; SC 200 mg weekly and single SC doses of 120 x 2, 200 or 240 mg.",
    regions = "Global (BEL110751, BEL110752, BEL112341, BEL114055), North East Asia (BEL113750), USA (BEL114448), Japan (BEL116119) and Mainland China (200909, 209629).",
    fat_free_mass_range = "12.6-77.4 kg (median 40.5)",
    albumin_range = "18-55 g/L (median 41.0)",
    igg_range = "4.1-61.9 g/L (median 14.7)",
    notes = "Pooled analysis of nine GSK studies (Table 1); demographics from Supporting Information Table S1. Concentrations below the 100 ng/mL LLOQ were excluded."
  )

  ini({
    # Supporting Information Table S2 (final population PK model, NONMEM
    # ADVAN4 TRANS4, FOCEI). Published in mL/day and mL; divided by 1000 here
    # to give L/day and L, so central / vc is in mg/L = ug/mL.
    lcl <- log(0.223); label("Clearance for the reference participant (L/day)") # Table S2 theta1 = 223 mL/day (RSE 1%)
    lvc <- log(2.660); label("Central volume of distribution for the reference participant (L)") # Table S2 theta2 = 2660 mL (RSE 1.2%)
    lq <- log(0.494); label("Intercompartmental clearance for the reference participant (L/day)") # Table S2 theta4 = 494 mL/day (RSE 2.7%)
    lvp <- log(2.360); label("Peripheral volume of distribution for the reference participant (L)") # Table S2 theta5 = 2360 mL (RSE 2.1%)
    lka <- log(0.226); label("First-order SC absorption rate constant (1/day)") # Table S2 theta6 = 0.226 /day (RSE 3.8%)
    ltlag <- log(0.170); label("SC absorption lag time (day)") # Table S2 theta7 ALAG = 0.170 day (RSE 2.9%)
    # Table S2 equation F = 1/(1 + exp(theta3)), i.e. F = expit(-theta3)
    logitfdepot <- 1.40; label("Logit of SC bioavailability (unitless logit)") # Table S2 theta3 = -1.40 (logit, sign reversed by F = 1/(1+exp(theta3))); F = 80.2%

    e_ffm_cl <- 0.739; label("Power exponent of fat-free mass on CL and Q (unitless)") # Table S2 theta10 = 0.739 (shared by CL and Q)
    e_ffm_vc <- 0.711; label("Power exponent of fat-free mass on Vc and Vp (unitless)") # Table S2 theta11 = 0.711 (shared by Vc and Vp)
    e_alb_cl <- -1.11; label("Power exponent of baseline albumin on CL (unitless)") # Table S2 theta12 = -1.11 (CL block)
    e_igg_cl <- 0.313; label("Power exponent of baseline IgG on CL (unitless)") # Table S2 theta13 = 0.313
    e_region_china_cl <- 1.28; label("Multiplicative factor on CL for Mainland China participants (unitless)") # Table S2 theta14 = 1.28
    e_region_china_vc <- 1.20; label("Multiplicative factor on Vc for Mainland China participants (unitless)") # Table S2 theta15 = 1.20
    age50_vc <- 2.69; label("Age at half-maximal Vc in the AGE/(AGE + age50_vc) term (years)") # Table S2 theta16 Age row: printed 0.311 is the SE; 95% CI 2.08-3.30 midpoint 2.69, RSE 0.311/2.69 = 11.6%

    # IIV (Table S2 omega^2; footnote e Cov(CL, Vc) = 0.0168). omega3 (F) and
    # omega4 (Q) are 0 in the source and are omitted.
    etalcl + etalvc ~ c(
      0.0834,
      0.0168, 0.0076
    )
    etalvp ~ 0.0917 # Table S2 omega5^2 (CV 30.3%)
    etalka ~ 0.218 # Table S2 omega6^2 (CV 46.7%)
    etaltlag ~ 0.0359 # Table S2 omega7^2 (CV 18.9%)
    etaRUV ~ 0.272 # Table S2 omega8^2 on the residual error magnitude (CV 52.2%)

    # Residual error (Uppsala approach: SDs estimated as thetas, SIGMA 1 FIX)
    propSd <- 0.205; label("Proportional residual SD (fraction)") # Table S2 theta8 epsilon1 = 20.5%
    addSd <- 0.605; label("Additive residual SD (ug/mL)") # Table S2 theta9 epsilon2 = 0.605 mg/L
  })
  model({
    # Covariate effects (Table S2 equations, reference participant FFM
    # 40.486 kg, albumin 41 g/L, IgG 14.7 g/L, not Mainland China)
    ffm_cl <- (FFM / 40.486)^e_ffm_cl
    ffm_vc <- (FFM / 40.486)^e_ffm_vc
    alb_cl <- (ALB / 41)^e_alb_cl
    igg_cl <- (IGG / 14.7)^e_igg_cl
    china_cl <- e_region_china_cl^REGION_CHINA
    china_vc <- e_region_china_vc^REGION_CHINA
    age_vc <- AGE / (AGE + age50_vc)

    cl <- exp(lcl + etalcl) * ffm_cl * alb_cl * igg_cl * china_cl
    vc <- exp(lvc + etalvc) * ffm_vc * china_vc * age_vc
    q <- exp(lq) * ffm_cl
    vp <- exp(lvp + etalvp) * ffm_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- expit(logitfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # dose mg, vc L -> mg/L = ug/mL
    Cc <- central / vc

    # Table S2: W = sqrt(eps1^2 * IPRED^2 + eps2^2) * exp(eta), the eta being
    # the omega8 residual variability
    ruv_scale <- exp(etaRUV)
    propSd_i <- propSd * ruv_scale
    addSd_i <- addSd * ruv_scale
    Cc ~ prop(propSd_i) + add(addSd_i)
  })
}
