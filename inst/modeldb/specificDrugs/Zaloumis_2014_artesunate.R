Zaloumis_2014_artesunate <- function() {
  description <- paste(
    "One-compartment population PK model of dihydroartemisinin (DHA), the active",
    "metabolite of intravenous artesunate (ARS), in 266 adults and children with",
    "severe or moderately severe Plasmodium falciparum malaria pooled from six",
    "studies in Africa and Asia (Zaloumis 2014, final model fitted to all patients",
    "including period two of the crossover trials). ARS is assumed to convert",
    "instantly and completely to DHA, so the ARS dose enters the DHA compartment as",
    "an IV bolus scaled by the DHA/ARS molecular-weight ratio 284/384. Clearance and",
    "volume scale allometrically with body weight (fixed exponents 0.75 and 1,",
    "reference 15 kg); clearance carries a fixed exponential age-maturation function",
    "(value at birth 1/3.1 of the mature value, rate 0.082 per month) and volume is",
    "14% lower in males. Residual error is additive on the log scale with a",
    "separate SD for each of the six studies."
  )
  reference <- paste(
    "Zaloumis SG, Tarning J, Krishna S, Price RN, White NJ, Davis TME,",
    "McCaw JM, Olliaro P, Maude RJ, Kremsner P, Dondorp A, Gomes M, Barnes K,",
    "Simpson JA. Population pharmacokinetics of intravenous artesunate: a pooled",
    "analysis of individual data from patients with severe malaria.",
    "CPT Pharmacometrics Syst Pharmacol. 2014;3(11):e145.",
    "doi:10.1038/psp.2014.43"
  )
  vignette <- "Zaloumis_2014_artesunate"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at admission",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL (exponent 0.75) and V (exponent 1), both fixed, centred on the pooled median 15 kg (Supplementary Information 4; NONMEM code WGTmed = 15). Time-fixed in the source data set.",
      source_name = "WGT"
    ),
    AGE = list(
      description = "Age at admission",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the fixed exponential maturation of CL. The NONMEM code converts age to months (AGEmth = AGE * 12) before applying the 0.082 per-month rate; model() does the same conversion.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "Source column SEX is coded 1 = male, 0 = female; SEXF = 1 - SEX. Male volume is V * (1 + e_male_vc), so the typical values exp(lcl) and exp(lvc) are for a female patient (Table 3 footnote a).",
      source_name = "SEX"
    ),
    STUDY_KRISHNA2001 = list(
      description = "1 = record from the Krishna et al. 2001 Ghanaian paediatric crossover trial (study 2, conducted 1996)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the study-specific residual SD only. Reference (all STUDY_* = 0) is Kremsner et al. 2012 (study 1). Source column STDY (code 4 in the NONMEM data set).",
      source_name = "STDY"
    ),
    STUDY_NEALON2002 = list(
      description = "1 = record from the Nealon et al. 2002 Gabonese paediatric crossover trial (study 3, conducted 2000-2001)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the study-specific residual SD only. Source column STDY (code 5 in the NONMEM data set).",
      source_name = "STDY"
    ),
    STUDY_MAUDE2009 = list(
      description = "1 = record from the Maude et al. 2009 Bangladeshi adult study (study 4, conducted 2006-2007)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the study-specific residual SD only. Source column STDY (code 6 in the NONMEM data set).",
      source_name = "STDY"
    ),
    STUDY_WHO_SIMPSON2006 = list(
      description = "1 = record from the WHO crossover trial in Thai adults with moderately severe malaria (study 5, conducted 1996; subset published in Simpson et al. 2006)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the study-specific residual SD only. Source column STDY (code 8 in the NONMEM data set).",
      source_name = "STDY"
    ),
    STUDY_DAVIS2001 = list(
      description = "1 = record from the Davis et al. 2001 Vietnamese adult study (study 6, conducted 1997)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the study-specific residual SD only. Source column STDY (code 9 in the NONMEM data set).",
      source_name = "STDY"
    )
  )

  compartmentData <- list(
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 266L,
    n_studies = 6L,
    n_observations = 1014L,
    age_range = "6 months to 65 years (children median 3.47 years, range 0.57-9.26; adults median 24 years, range 16-65)",
    weight_range = "6-75 kg (children median 13 kg, range 6-29; adults median 50 kg, range 33-75)",
    weight_typical = "15 kg (pooled median; allometric reference)",
    sex_female_pct = 46,
    race_ethnicity = "African children (Gabon, Malawi, Ghana) and South and South-East Asian adults (Bangladesh, Thailand, Vietnam)",
    disease_state = "Severe or moderately severe Plasmodium falciparum malaria",
    dose_range = "IV artesunate 2.4 mg/kg (most patients), 4 mg/kg (one Kremsner arm) or 1.2 mg/kg (Nealon period two); some crossover patients received IR or IM artesunate in the other period",
    regions = "Gabon, Malawi, Ghana, Bangladesh, Thailand, Vietnam",
    notes = "195 children and 71 adults (Table 2, subgroup 1; 52% of children and 58% of adults male, giving about 46% female overall). Samples collected within 12 h of the IV dose; 3.6% of DHA concentrations below the LLOQ were handled by the M3 method in the source analysis (Supplementary Table S4)."
  )

  ini({
    # Typical values are for a female patient weighing 15 kg and fully
    # matured (Table 3 footnote a).
    lcl <- log(13.0)
    label("DHA clearance, CL, for a 15-kg female (L/h)") # Table 3 'All severe malaria patients': population mean CL = 13.0 L/h (%RSE 4.42)
    lvc <- log(11.56)
    label("DHA volume of distribution, V, for a 15-kg female (L)") # Table 3 'All severe malaria patients': population mean V = 11.56 L (%RSE 5.05)

    # Fixed allometric exponents and maturation constants (Supplementary
    # NONMEM code: CL = (TVCL/3.1)*(WGT/15)**0.75*EXP(-0.082*AGEmth) +
    # TVCL*(WGT/15)**0.75*(1 - EXP(-0.082*AGEmth)); V = TVV*(WGT/15)).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL (unitless)") # Supplementary NONMEM code: (WGT/WGTmed)**0.75
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V (unitless)") # Supplementary NONMEM code: TVV*(WGT/WGTmed)
    cl_pna_famp <- fixed(1 / 3.1 - 1)
    label("Fractional difference of CL at age 0 from mature CL (unitless)") # Supplementary NONMEM code: TVCL/3.1 at age 0 -> 1/3.1 - 1 = -0.677
    lcl_pna_kdes <- fixed(log(0.082))
    label("Rate of CL maturation with age (1/month)") # Supplementary NONMEM code: EXP(-0.082*AGEmth), AGEmth = AGE*12

    e_male_vc <- -0.14
    label("Fractional change in V for males vs females (unitless)") # Table 3 'All severe malaria patients': beta V,MALE = -0.14 (%RSE 39.1)

    # Table 3 reports SDs and the correlation: var = SD^2, cov = cor*SD_CL*SD_V
    # = 0.61 * 0.59 * 0.44 = 0.158356.
    etalcl + etalvc ~ c(0.3481, 0.158356, 0.1936) # Table 3 'All severe malaria patients': SD eta_CL = 0.59, SD eta_V = 0.44, cor = 0.61

    # Residual SD on the natural-log scale, one per study (Table 3 footnote:
    # sigma_i for study i = 1 Kremsner, 2 Krishna, 3 Nealon, 4 Maude, 5 WHO,
    # 6 Davis). expSd is the reference study (Kremsner et al. 2012).
    expSd <- 0.34
    label("Log-scale residual SD, Kremsner et al. 2012 study (unitless)") # Table 3 'All severe malaria patients': sigma_1 (2007-2008) = 0.34
    expSdKrishna2001 <- 0.36
    label("Log-scale residual SD, Krishna et al. 2001 study (unitless)") # Table 3 'All severe malaria patients': sigma_2 (1996) = 0.36
    expSdNealon2002 <- 0.43
    label("Log-scale residual SD, Nealon et al. 2002 study (unitless)") # Table 3 'All severe malaria patients': sigma_3 (2000-2001) = 0.43
    expSdMaude2009 <- 0.36
    label("Log-scale residual SD, Maude et al. 2009 study (unitless)") # Table 3 'All severe malaria patients': sigma_4 (2006-2007) = 0.36
    expSdWhoSimpson2006 <- 0.30
    label("Log-scale residual SD, WHO study (unitless)") # Table 3 'All severe malaria patients': sigma_5 (1996) = 0.30
    expSdDavis2001 <- 0.27
    label("Log-scale residual SD, Davis et al. 2001 study (unitless)") # Table 3 'All severe malaria patients': sigma_6 (1997) = 0.27
  })

  model({
    # Age maturation of CL, age converted to months as in the NONMEM code.
    age_mo <- AGE * 12
    mat_cl <- 1 + cl_pna_famp * exp(-exp(lcl_pna_kdes) * age_mo)

    cl <- exp(lcl + etalcl) * (WT / 15)^e_wt_cl * mat_cl
    vc <- exp(lvc + etalvc) * (WT / 15)^e_wt_vc * (1 + e_male_vc * (1 - SEXF))

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose is artesunate (mg); complete instantaneous conversion to DHA is
    # applied as the molecular-weight ratio 284/384 used by the authors.
    f(central) <- 284 / 384

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc

    expSdStudy <- expSd *
      (1 -
        STUDY_KRISHNA2001 -
        STUDY_NEALON2002 -
        STUDY_MAUDE2009 -
        STUDY_WHO_SIMPSON2006 -
        STUDY_DAVIS2001) +
      expSdKrishna2001 * STUDY_KRISHNA2001 +
      expSdNealon2002 * STUDY_NEALON2002 +
      expSdMaude2009 * STUDY_MAUDE2009 +
      expSdWhoSimpson2006 * STUDY_WHO_SIMPSON2006 +
      expSdDavis2001 * STUDY_DAVIS2001
    Cc ~ lnorm(expSdStudy)
  })
}
