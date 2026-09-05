Feng_2023_cyclosporine <- function() {
  description <- "One-compartment intravenous population PK model for cyclosporine A in Chinese paediatric patients undergoing allogeneic haematopoietic stem cell transplantation (Feng 2023; final covariate model with body weight on CL and body weight plus haematocrit on V)"
  reference <- "Feng H, Wang X, Zheng W, Liu S, Jiang H, Lin Y, Qiu H, Chan TF, Huang M, Li Y, Mo X, Li J. Initial dosage optimisation of cyclosporine in Chinese paediatric patients undergoing allogeneic haematopoietic stem cell transplantation based on population pharmacokinetics: a retrospective study. BMJ Paediatr Open. 2023;7(1):e002003. doi:10.1136/bmjpo-2023-002003"
  vignette <- "Feng_2023_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    # Cyclosporine A was measured in whole blood by an enzyme-multiplied
    # immunoassay technique (EMIT) on a Viva-E analyser (Feng 2023 Methods,
    # "Participants and data collection"), so the modelled matrix is whole
    # blood rather than plasma.
    central = list(analyte = "cyclosporine", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling on both CL and V, normalised to the 16.5 kg cohort",
        "median (Feng 2023 Eq. 7 and Eq. 8; Table 1 reports body weight",
        "median 16.5 kg, range 8.0-64.5 kg). Baseline body weight; the paper",
        "does not describe a time-varying weight record. Body weight was the",
        "only covariate retained on CL. The Discussion states that age and",
        "body weight were collinear and that body weight alone was kept in",
        "the final model."
      ),
      source_name        = "BW"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling on V only, normalised to 28.8 % (Feng 2023 Eq. 8).",
        "The paper prints haematocrit on the percent scale throughout",
        "(Table 1 median 27.0 %, range 5.0-37.5 %; Figure 5 simulates 10 %,",
        "30 % and 50 %), so no rescaling from a volume fraction is needed",
        "for the canonical percent column. Note that the 28.8 % normalising",
        "constant in Eq. 8 is not the 27.0 % median printed in Table 1 for",
        "the full 251-patient cohort; Eq. 8 was fitted on the 176-patient",
        "training set, whose demographics are in the online supplemental",
        "table S3 (not on disk). The exponent is negative: lower haematocrit",
        "leaves less red-cell mass for this highly lipophilic,",
        "erythrocyte-partitioned drug to bind, increasing distribution into",
        "fat and hence the apparent volume (Feng 2023 Discussion)."
      ),
      source_name        = "HCT"
    )
  )

  # Covariates that Feng 2023 screened but did not retain in the final model.
  # The paper is explicit that "none of the genetic polymorphisms in ABCB1,
  # CYP3A4, CYP3A5, POR and NR1I3 were significant covariates in the PK of CsA
  # in our study" (Discussion) and that, apart from body weight and
  # haematocrit, "other covariates were not found to have statistical
  # significance on PK parameters" (Results). These entries carry that
  # negative finding forward; they are documentation only and are not
  # referenced in model().
  covariatesDataExcluded <- list(
    AGE   = list(description = "Age", units = "years", type = "continuous",
                 notes = "Table 1 median 6 years (range 1-17). Screened; dropped for collinearity with body weight (Discussion)."),
    SEXF  = list(description = "Female sex indicator", units = "(binary)", type = "binary",
                 reference_category = "male", notes = "Table 1: 89/251 female (35.5 %). Screened, not significant."),
    ALB   = list(description = "Serum albumin", units = "g/L", type = "continuous",
                 notes = "Table 1 median 36.9 g/L (28.1-44.1). Screened, not significant."),
    TBILI = list(description = "Total bilirubin", units = "umol/L", type = "continuous",
                 notes = "Table 1 median 16.0 umol/L (3.6-105.1). Screened, not significant."),
    DBIL  = list(description = "Direct bilirubin", units = "umol/L", type = "continuous",
                 notes = "Table 1 median 6.1 umol/L (0.9-82.7). Screened, not significant."),
    ALT   = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                 notes = "Table 1 median 51 U/L (3-1120). Screened, not significant."),
    AST   = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
                 notes = "Table 1 median 26 U/L (1-606). Screened, not significant."),
    ALP   = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                 notes = "Table 1 median 119 U/L (46-291). Screened, not significant."),
    CRP   = list(description = "Hypersensitive C-reactive protein", units = "mg/L", type = "continuous",
                 notes = "Table 1 median 7.00 mg/L (0.20-263.40). Screened, not significant."),
    HGB   = list(description = "Haemoglobin", units = "g/L", type = "continuous",
                 notes = "Table 1 median 94 g/L (60-124). Screened, not significant; collinear with the retained HCT."),
    RBC   = list(description = "Red blood cell count", units = "10^12/L", type = "continuous",
                 notes = "Table 1 median 3.32 (2.12-5.28). Screened, not significant; collinear with the retained HCT."),
    PLT   = list(description = "Platelet count", units = "10^9/L", type = "continuous",
                 notes = "Table 1 median 24 (2-180). Screened, not significant."),
    SNP_ABCB1_RS1045642  = list(description = "ABCB1 rs1045642 (3435C>T) genotype", units = "(genotype)", type = "categorical",
                                reference_category = "GG", notes = "Table 2 genotype frequencies AA 0.15 / GA 0.43 / GG 0.42. Screened, not significant."),
    SNP_ABCB1_RS1128503  = list(description = "ABCB1 rs1128503 (1236C>T) genotype", units = "(genotype)", type = "categorical",
                                reference_category = "GG", notes = "Table 2 genotype frequencies AA 0.34 / AG 0.49 / GG 0.17. Screened, not significant."),
    SNP_ABCB1_RS34800935 = list(description = "ABCB1 rs34800935 genotype", units = "(genotype)", type = "categorical",
                                reference_category = "CC", notes = "Table 2 genotype frequencies CC 0.15 / TC 0.52 / TT 0.33. Screened, not significant."),
    SNP_ABCB1_RS3842     = list(description = "ABCB1 rs3842 genotype", units = "(genotype)", type = "categorical",
                                reference_category = "CC", notes = "Table 2 genotype frequencies CC 0.09 / TC 0.44 / TT 0.47. Screened, not significant."),
    SNP_CYP3A4_RS2242480 = list(description = "CYP3A4*1G (rs2242480) genotype", units = "(genotype)", type = "categorical",
                                reference_category = "CC", notes = "Table 2 genotype frequencies CC 0.50 / CT 0.43 / TT 0.08. Screened, not significant."),
    SNP_CYP3A5_RS776746  = list(description = "CYP3A5*3 (rs776746) genotype", units = "(genotype)", type = "categorical",
                                reference_category = "CC", notes = "Table 2 genotype frequencies CC 0.52 / CT 0.38 / TT 0.10. Screened, not significant."),
    SNP_POR_RS17685      = list(description = "POR rs17685 genotype", units = "(genotype)", type = "categorical",
                                reference_category = "GG", notes = "Table 2 genotype frequencies AA 0.12 / GA 0.52 / GG 0.37. Screened, not significant."),
    SNP_NR1I3_RS2307424  = list(description = "NR1I3 (CAR) rs2307424 genotype", units = "(genotype)", type = "categorical",
                                reference_category = "GG", notes = "Table 2 genotype frequencies AA 0.27 / AG 0.47 / GG 0.26. Screened, not significant.")
  )

  population <- list(
    species        = "human",
    n_subjects     = 251L,
    n_studies      = 1L,
    age_range      = "1-17 years",
    age_median     = "6 years",
    weight_range   = "8.0-64.5 kg",
    weight_median  = "16.5 kg",
    sex_female_pct = 35.5,
    race_ethnicity = "Chinese (single-centre cohort, Guangzhou)",
    disease_state  = "Paediatric recipients of allogeneic haematopoietic stem cell transplantation receiving cyclosporine A for acute graft-versus-host disease prophylaxis; 183/251 (72.9 %) beta-thalassaemia, 68/251 (27.1 %) other indications",
    dose_range     = "Intravenous infusion, usually 3 mg/kg/day divided every 12 h starting on day 1 of transplantation; beta-thalassaemia patients commonly started 1.5 mg/kg/day every 12 h from 10 days before transplantation and were escalated to 3 mg/kg/day on day 1. Doses were then adjusted by therapeutic drug monitoring to a trough of 150-200 ng/mL.",
    regions        = "Guangzhou, Guangdong, China (Guangzhou Women and Children's Medical Center)",
    n_observations = 865L,
    notes          = paste(
      "Retrospective single-centre analysis, January 2016 to December 2020,",
      "registered as ChiCTR2000040561. 865 whole-blood cyclosporine trough",
      "concentrations (Table 1 median C0 117.2 ng/mL, range 46.3-445.4)",
      "measured by enzyme-multiplied immunoassay technique on a Viva-E",
      "analyser. The 251 patients were split 7:3 by transplant date into a",
      "176-patient model-building set and a 75-patient external-validation",
      "set; the parameters in this file are the final model fitted to the",
      "176-patient training set (Table 3). Estimation used first-order",
      "conditional estimation with extended least squares in Phoenix NLME",
      "v7.0 (Certara). Baseline demographics are Feng 2023 Table 1;",
      "genotype frequencies are Table 2."
    )
  )

  ini({
    # Structural parameters at the paper's reference covariates
    # (WT = 16.5 kg, HCT = 28.8 %). Feng 2023 Table 3, "Final model" column.
    #
    # The data are trough-only (865 C0 samples), so V is supported almost
    # entirely by the rate of accumulation across days rather than by any
    # within-interval concentration decline. The resulting terminal
    # half-life, ln(2) * 2033.53 / 14.47 = 97 h, is much longer than the
    # 6-8 h usually quoted for cyclosporine; it is reproduced here as
    # published. This is what makes haematocrit -- which enters V only --
    # move the day-7 trough at all in Feng 2023 Figure 5: at true steady
    # state a volume covariate cannot change a trough concentration.
    lcl <- log(14.47);   label("Typical clearance CL at reference body weight (L/h)")                 # Feng 2023 Table 3: CL = 14.47 L/h (RSE 4.61 %, bootstrap median 14.48, 95% CI 13.00-15.73)
    lvc <- log(2033.53); label("Typical volume of distribution V at reference body weight and haematocrit (L)")  # Feng 2023 Table 3: V = 2033.53 L (RSE 8.78 %, bootstrap median 2024.48, 95% CI 1657.30-2416.52)

    # Covariate effects (power form; Feng 2023 Eq. 7 and Eq. 8):
    #   CL (L/h) = 14.47   * (BW / 16.5)^0.99
    #   V  (L)   = 2033.53 * (BW / 16.5)^1.00 * (HCT / 28.8)^-0.39
    # All three exponents were estimated (each is reported with an RSE and a
    # bootstrap 95% CI in Table 3), so none is wrapped in fixed().
    e_wt_cl  <-  0.99; label("Power exponent of (WT / 16.5 kg) on CL (unitless)")    # Feng 2023 Table 3: theta_BW,CL = 0.99 (RSE 9.54 %, bootstrap median 1.00, 95% CI 0.82-1.21); Eq. 7
    e_wt_vc  <-  1.00; label("Power exponent of (WT / 16.5 kg) on V (unitless)")     # Feng 2023 Table 3: theta_BW,V = 1.00 (RSE 21.04 %, bootstrap median 0.99, 95% CI 0.61-1.43); Eq. 8
    e_hct_vc <- -0.39; label("Power exponent of (HCT / 28.8 %) on V (unitless)")     # Feng 2023 Table 3: theta_HCT,V = -0.39 (RSE -18.94 %, bootstrap median -0.39, 95% CI -0.62 to -0.27); Eq. 8

    # Inter-individual variability, exponential model (Feng 2023 Eq. 1:
    # P_i = P_pop * exp(eta_i), eta ~ N(0, omega^2)). Table 3 prints these
    # two rows as "omega^2 CL" and "omega^2 V", i.e. on the VARIANCE scale,
    # in contrast with the residual row printed as "sigma". The variances
    # below are therefore transcribed as printed:
    #   omega^2 CL = 0.10 -> SD 0.316, CV 32.4 %
    #   omega^2 V  = 0.22 -> SD 0.469, CV 49.4 %
    # No off-diagonal covariance is reported, so the two etas are
    # independent here.
    etalcl ~ 0.10  # Feng 2023 Table 3: omega^2 CL = 0.10 (RSE 16.59 %, bootstrap median 0.10)
    etalvc ~ 0.22  # Feng 2023 Table 3: omega^2 V  = 0.22 (RSE 23.52 %, bootstrap median 0.21)

    # Residual error, proportional (Feng 2023 Eq. 2: Cobs = C * (1 + eps),
    # eps ~ N(0, sigma^2)). Table 3 prints this row as "sigma proportional",
    # on the SD scale, so it is used directly as propSd.
    propSd <- 0.24; label("Proportional residual error (fraction)")  # Feng 2023 Table 3: sigma proportional = 0.24 (RSE 4.20 %, bootstrap median 0.24, 95% CI 0.22-0.26)
  })
  model({
    # 1. Individual parameters with covariate effects (Feng 2023 Eq. 7, Eq. 8)
    cl <- exp(lcl + etalcl) * (WT / 16.5)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 16.5)^e_wt_vc * (HCT / 28.8)^e_hct_vc

    # 2. Micro-constant
    kel <- cl / vc

    # 3. ODE system. Cyclosporine A was given only by intravenous infusion in
    # this cohort, so doses enter the central compartment directly and there
    # is no absorption parameter (Table 3 reports none).
    d/dt(central) <- -kel * central

    # 4. Observation and error. central is in mg and vc in L, so central / vc
    # is mg/L = ug/mL; the factor 1000 converts to the ng/mL (= ug/L) that
    # Feng 2023 reports throughout.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
