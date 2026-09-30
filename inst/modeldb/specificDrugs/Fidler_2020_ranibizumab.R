Fidler_2020_ranibizumab <- function() {
  description <- "One-compartment population PK model with first-order absorption for serum ranibizumab after bilateral intravitreal injection in preterm infants with retinopathy of prematurity (RAINBOW trial; Fidler et al. 2020, TVST). The vitreous acts as a first-order depot (Ka = ocular elimination rate, flip-flop kinetics) into a systemic central compartment. Typical values are expressed for a 70-kg adult and scaled to the infant by body weight (fixed allometric exponent 0.75 on CL/F, linear on V/F); CL/F also carries a fixed creatinine-clearance adjustment evaluated at the study-median infant CrCl (54.82 mL/min, modified Schwartz) relative to the adult median (65.22 mL/min) with the adult exponent 0.266. Correlated log-normal IIV on CL/F, V/F and Ka; log-normal residual error."
  reference <- "Fidler M, Fleck BW, Stahl A, Marlow N, Chastain JE, Li J, Lepore D, Reynolds JD, Chiang MF, Fielder AR; RAINBOW study group. Ranibizumab population pharmacokinetics and free VEGF pharmacodynamics in preterm infants with retinopathy of prematurity in the RAINBOW trial. Transl Vis Sci Technol. 2020;9(8):43. doi:10.1167/tvst.9.8.43"
  vignette <- "Fidler_2020_ranibizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    vitreous = list(analyte = "ranibizumab", units = "mg", specimen = "vitreous", verified = TRUE),
    central = list(analyte = "ranibizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at treatment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to a 70-kg adult reference: (WT/70)^0.75 on CL/F and (WT/70)^1 on V/F, both exponents fixed (Supplementary Appendix, 'Final model equations'). RAINBOW baseline weight median 1.7-1.8 kg, range 0.8-4.2 kg (Table 1).",
      source_name = "w"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 95L,
    n_studies = 1L,
    age_range = "Preterm infants; postnatal age at baseline mean 10.8-11.1 weeks (SD 3.9-4.5) by arm; gestational age at birth mean 25.8-26.5 weeks",
    weight_range = "0.8-4.2 kg at baseline (median 1.7-1.8 kg by arm)",
    sex_female_pct = 46.3,
    race_ethnicity = c(Caucasian = 59.1, Asian = 31.5, Black = 2.7, Other = 6.7),
    disease_state = "Retinopathy of prematurity (ROP) requiring treatment",
    dose_range = "Single bilateral intravitreal injection of ranibizumab 0.1 mg or 0.2 mg per eye (0.2 or 0.4 mg total per infant); a laser-therapy arm (n = 69) contributed free-VEGF data only.",
    regions = "26 countries; investigator sites (Supplementary Table S3) in the USA, Mexico, Japan, Taiwan, India, Malaysia, Saudi Arabia, Egypt, Turkey, Russia and Europe.",
    notes = "RAINBOW trial. PK sampling in odd-numbered ranibizumab-treated infants at day 1, day 15 (7-21) and day 29 (22-28); 115 (0.1 mg) + 209 (0.2 mg) serum concentrations, 10 excluded after laser rescue (Results). Demographics summarised over the 149 ranibizumab-treated infants (76 + 73) of Table 1; percentages above are pooled over those two arms. PK analysis: 95 infants with serum ranibizumab observations (45 in 0.1 mg, 50 in 0.2 mg; Table 1). BLQ samples (LLOQ 0.015 ng/mL) were treated as missing. Plasma free VEGF showed no relationship with ranibizumab exposure, so no PD model was developed."
  )

  ini({
    # Typical values are for a 70-kg adult (Final Model: Parameter Estimates,
    # 'the model estimated clearances and volumes for the 70-kg adult values
    # instead of the infant values'). The Supplementary Appendix quotes
    # CLa = 23.8 L/day and Va = 2.97 L -- those are the prior adult-model
    # starting values; Table 2 holds the re-estimated final values used here.
    lcl <- log(28.48); label("Apparent clearance CL/F for a 70-kg adult (L/day)") # Table 2, CL/F = 28.48 L/day (RSE 3.96%)
    lvc <- log(27.58); label("Apparent volume of distribution V/F for a 70-kg adult (L)") # Table 2, V/F = 27.58 L (RSE 3.35%)
    lka <- log(0.12); label("First-order absorption rate from the vitreous = ocular elimination rate Ka (1/day)") # Table 2, Ka = 0.12 1/day (RSE 2.56%)

    # Fixed covariate terms (Supplementary Appendix, 'Final model equations')
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Supplementary Appendix, CLi = CLa * (w/70)^(3/4) * ...
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Supplementary Appendix, Vi = Va * (w/70) * exp(eta)
    e_crcl_cl <- fixed(0.266); label("Power exponent of creatinine clearance on CL/F, carried from the adult model (unitless)") # Supplementary Appendix, 'thetaCRCL = 0.266 ... assumed to be the same as that of the adult population'

    # IIV: Table 2 reports BSV as log-normal %CV = 100*sqrt(exp(omega^2) - 1),
    # so omega^2 = log(1 + CV^2). Covariances: the final model (Supplementary
    # Table S1, step 5) carries Cov(V, CL) and Cov(Ka, V) but not Cov(Ka, CL).
    # Final-model correlations are not tabulated; the correlation estimated
    # when each covariance was added (Table S1) is used:
    #   rho(CL, V) = 0.491 (step 3); rho(V, Ka) = -0.55 (step 5 = final model).
    # cov(CL, V) = 0.491 * sqrt(0.11506 * 2.24522) = 0.24956
    # cov(V, Ka) = -0.55 * sqrt(2.24522 * 0.026159) = -0.13329
    etalcl + etalvc + etalka ~ c(
      0.11506,
      0.24956, 2.24522,
      0, -0.13329, 0.026159
    ) # Table 2 BSV CV CL/F 34.92%, V/F 290.56%, Ka 16.28%; Supplementary Table S1 correlations 0.491 and -0.55

    expSd <- 0.60; label("Log-normal residual SD (log ng/mL)") # Table 2, 'Lognormal SD' = 0.60 (RSE 7.09%)
  })

  model({
    # Creatinine-clearance adjustment. The paper applies the STUDY-MEDIAN
    # infant CrCl (54.82 mL/min, modified Schwartz) for every infant, not
    # individual values, because site-measured serum creatinine was
    # unreliable (Supplementary Appendix). The adult reference is the median
    # adult Cockcroft-Gault CrCl of 65.22 mL/min. The term is therefore a
    # constant multiplier on CL/F (0.9548), kept explicit for traceability.
    crcl_adj <- (54.82 / 65.22)^e_crcl_cl

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * crcl_adj
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    # vitreous: ranibizumab amount in the (pooled bilateral) vitreous depot;
    # the dose is the total amount injected into both eyes.
    d/dt(vitreous) <- -ka * vitreous
    d/dt(central) <- ka * vitreous - kel * central

    # central amount in mg, vc in L: mg/L * 1000 = ng/mL
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
