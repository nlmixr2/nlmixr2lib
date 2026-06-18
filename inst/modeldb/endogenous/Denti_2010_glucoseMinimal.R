Denti_2010_glucoseMinimal <- function() {
  description <- "Bergman glucose minimal model fitted to insulin-modified IVGTT data from 204 healthy adults aged 18-87 years (Mayo Clinic cohort, Basu 2003/2006); population NLME implementation with mean-centred covariate effects of age, sex, visceral abdominal fat (CT-derived), percent total body fat (DEXA-derived), basal glucose, and basal insulin on the four mechanistic parameters glucose effectiveness (SG), insulin sensitivity (SI), insulin-action rate constant (P2), and apparent glucose distribution volume per kg body mass (VOL). Plasma insulin is supplied as a time-varying error-free forcing function (regressor INS)."
  reference <- paste(
    "Denti P, Bertoldo A, Vicini P, Cobelli C. (2010).",
    "IVGTT glucose minimal model covariate selection by nonlinear mixed-effects approach.",
    "Am J Physiol Endocrinol Metab 298(5):E950-E960.",
    "doi:10.1152/ajpendo.00656.2009.",
    sep = " "
  )
  vignette <- "Denti_2010_glucoseMinimal"
  paper_specific_compartments <- c("insulin_action")

  units <- list(
    time = "min",
    dosing = "mg/kg (glucose)",
    concentration = "mg/dL"
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean-centred at 55.53 years (Denti 2010 Table 1 pooled-population mean). Retained on lvd, lsi, and lp2 in the final covariate model (Table 5).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Denti 2010 reported two typical values for VOL (1.70 dL/kg male reference, 1.65 dL/kg female). Encoded here as a male-reference intercept lvd plus a log-scale female adjustment e_sexf_vd = log(1.65 / 1.70) so VOL_female = exp(lvd + e_sexf_vd) = 1.65 dL/kg exactly. Equivalent to the paper's two-typical-values form.",
      source_name = "SEX"
    ),
    VISCERAL_ABDOMINAL_FAT = list(
      description = "Visceral abdominal fat area from a single CT slice at the L2 / L3 vertebral level",
      units = "cm^2 / CT slice",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean-centred at 141.8 cm^2 / CT slice (Denti 2010 Table 1 pooled-population mean). Retained on lsi only (Table 5). Three subjects had missing values in the source paper and were imputed to the pooled mean; the mean-centred form sets their covariate deviation to zero so no covariate adjustment is applied for those subjects.",
      source_name = "VAF"
    ),
    BODYFAT_PCT = list(
      description = "Percent total body fat (DEXA-derived)",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean-centred at 32.39 % (Denti 2010 Table 1 pooled-population mean). Retained on lvd only (Table 5). DEXA-derived in the founding paper (Basu 2003 protocol). Same three subjects with missing VISCERAL_ABDOMINAL_FAT also missing BODYFAT_PCT; imputed to the pooled mean per the paper.",
      source_name = "%TBF"
    ),
    FPG = list(
      description = "Basal (fasting, predose) plasma glucose concentration; supplies Gb to the structural model",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean-centred at 91.34 mg/dL (Denti 2010 Table 1 pooled-population mean). Two roles in this model: (a) covariate on lvd (Table 5 e_fpg_vd = -0.00312 per mg/dL), and (b) Gb in the structural ODE -- the steady-state production rate sg * FPG * vd holds plasma glucose at FPG before the IVGTT bolus and the bolus initial-condition G(0) = FPG + dose / vd is reproduced automatically once the dose is applied. Per-subject baseline value, time-fixed. Units mg/dL per the source paper Table 1 (canonical FPG entry permits per-model mg/dL via covariateData[[FPG]]$units).",
      source_name = "GBSL"
    ),
    INS_BL = list(
      description = "Basal (fasting, predose) plasma insulin concentration; supplies Ib to the structural model",
      units = "pmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean-centred at 26.98 pmol/L (Denti 2010 Table 1 pooled-population mean). NOTE on units: Denti 2010 Table 1 labels the IBSL column 'pmol/ml' but the reported mean 26.98 is biologically pmol/L -- 26.98 pmol/mL would equal 26,980 pmol/L, far above any plausible fasting human insulin value (typical fasting insulin in healthy adults is 30-100 pmol/L). The 'pmol/ml' label in Table 1 is therefore an apparent typo for pmol/L and is interpreted as pmol/L throughout this model; see the vignette Errata section. Two roles in this model: (a) covariate on lsi (Table 5 e_insbl_si = -0.0282 per pmol/L) and lp2 (Table 5 e_insbl_p2 = -0.0150 per pmol/L), and (b) Ib in the structural insulin-action ODE p2 * (si * (INS - INS_BL) - insulin_action). Per-subject baseline value, time-fixed.",
      source_name = "IBSL"
    ),
    INS = list(
      description = "Plasma insulin concentration time course (forcing function)",
      units = "pmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying regressor supplied at every observation / event row in the dataset; linearly interpolated between rows via `linear(INS)` declared in model(). Denti 2010 Methods p. E951: 'I(t) is assumed as a known, error-free input (forcing) function' -- the minimal model does not estimate plasma-insulin kinetics; users must supply the observed or simulated insulin profile from the IVGTT (typical full insulin-modified IVGTT samples insulin at 0, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 19, 22, 25, 30, 40, 50, 60, 90, 120, 180, 240 min).",
      source_name = "I"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height (screened, not retained in final model)",
      units = "cm",
      type = "continuous",
      notes = "Screened in Denti 2010 multivariate regression (Table 3 'best covariates' lists for log(SG) and log(SI)) but not retained in the final population model (Table 5). Recorded here as documentation of the full screening set. Source name BH."
    ),
    WT = list(
      description = "Body weight (screened, not retained)",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 3) but not retained. Source name BW. Denti 2010 Table 1 mean 77.94 kg, range 53-127 kg."
    ),
    BMI = list(
      description = "Body mass index (screened, not retained)",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Table 3) but not retained. Derived covariate from BH and BW; collinearity with BW / BSA prevented retention even where signal was present (Discussion p. E955)."
    ),
    BSA = list(
      description = "Body surface area (screened, not retained)",
      units = "m^2",
      type = "continuous",
      notes = "Screened (Table 3) but not retained. Derived covariate from BH and BW."
    ),
    LBM = list(
      description = "Lean body mass (screened, not retained)",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 3) but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 204L,
    n_studies = 1L,
    age_range = "18-87 years",
    age_median = "56 years (mean; median 65 years from Table 1)",
    weight_range = "53-127 kg",
    weight_median = "77.94 kg (mean; median 79 kg from Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Healthy adults (no diagnosis of glucose-metabolic disorders); BMI range 20-35 kg/m^2",
    dose_range = "Insulin-modified IVGTT: glucose bolus 0.30 g/kg IV at t = 0 followed by a 5-minute insulin infusion starting at t = 20 min (Steil 1993 / Basu 2003 protocol). Full sampling schedule, 21 samples over 240 min.",
    regions = "United States (Mayo Clinic, Rochester MN)",
    notes = "204 healthy subjects pooled from a previously published Mayo Clinic dataset (Basu et al. 2003 J Clin Endocrinol Metab 88:6068; Basu et al. 2006 Diabetes 55:2001). Body composition assayed by dual-energy X-ray absorptiometry (DEXA) and single-slice computed tomography at L2 / L3. Three subjects had missing VISCERAL_ABDOMINAL_FAT, total-abdominal-fat, and BODYFAT_PCT values; the source paper imputed these to the pooled-population means so the mean-centred covariate deviations were zero (Methods p. E951). Baseline demographics in Table 1, individual-parameter empirical distribution in Table 2, covariate-model parameter estimates in Table 5."
  )

  ini({
    # ----- Structural typical values (Denti 2010 Table 5 'Covariate Model' column) -----
    # The four Bergman minimal-model parameters are estimated on the log scale for
    # log-normal IIV (paper assumption, Methods p. E951). lvd uses male as the
    # reference category; the female typical value is recovered via e_sexf_vd.

    lvd <- log(1.70)         ; label("Apparent glucose volume per body mass, male reference (dL/kg)")    # Denti 2010 Table 5 (theta_VOL = 1.70 dL/kg for males; 1% RSE)
    lsg <- log(0.0192)       ; label("Glucose effectiveness reference (1/min)")                          # Denti 2010 Table 5 (theta_SG = 0.0192 1/min; 2% RSE)
    lsi <- log(5.83e-5)      ; label("Insulin sensitivity reference (L per min per pmol/L)")             # Denti 2010 Table 5 (theta_SI = 5.83e-5 L/(min*pmol); 5% RSE)
    lp2 <- log(0.0254)       ; label("Insulin-action rate-constant reference (1/min)")                   # Denti 2010 Table 5 (theta_P2 = 0.0254 1/min; 5% RSE)

    # ----- Covariate effects (Denti 2010 Table 5) -----
    # All effects are log-scale slopes applied to mean-centred covariate deviations:
    # log(parameter) shifts by coefficient * (covariate - pooled-population mean).
    # Pooled-population means come from Denti 2010 Table 1 (the paper used
    # sex-stratified means for VOL but did not publish them, so the pooled means
    # are used here for both sexes; documented in the vignette Errata).

    # VOL covariate effects: sex, age, percent total body fat, basal glucose.
    e_sexf_vd    <- log(1.65 / 1.70)  ; label("Effect of SEXF (female vs male) on lvd, log-scale")       # Denti 2010 Table 5 (theta_VOL_FEMALE = 1.65 dL/kg, theta_VOL_MALE = 1.70 dL/kg)
    e_age_vd     <- 0.00181           ; label("Effect of (AGE - 55.53 yr) on lvd")                       # Denti 2010 Table 5 (theta_VOL~AGE = 0.00181; 18% RSE)
    e_bodyfat_vd <- -0.0101           ; label("Effect of (BODYFAT_PCT - 32.39 %) on lvd")                # Denti 2010 Table 5 (theta_VOL~%TBF = -0.0101; 7% RSE)
    e_fpg_vd     <- -0.00312          ; label("Effect of (FPG - 91.34 mg/dL) on lvd")                    # Denti 2010 Table 5 (theta_VOL~GBSL = -0.00312; 23% RSE; sign confirmed via Discussion p. E956 and bootstrap CI -0.00484 to -0.00153 in Table 5 right-hand column)

    # SI covariate effects: age, visceral abdominal fat, basal insulin.
    e_age_si     <- -0.00810  ; label("Effect of (AGE - 55.53 yr) on lsi")                                                       # Denti 2010 Table 5 (theta_SI~AGE = -0.00810; 7% RSE)
    e_vaf_si     <- -0.00208  ; label("Effect of (VISCERAL_ABDOMINAL_FAT - 141.8 cm^2/slice) on lsi")                            # Denti 2010 Table 5 (theta_SI~VAF = -0.00208; 20% RSE)
    e_insbl_si   <- -0.0282   ; label("Effect of (INS_BL - 26.98 pmol/L) on lsi")                                                # Denti 2010 Table 5 (theta_SI~IBSL = -0.0282; 11% RSE)

    # P2 covariate effects: age, basal insulin.
    e_age_p2     <- -0.0110   ; label("Effect of (AGE - 55.53 yr) on lp2")                               # Denti 2010 Table 5 (theta_P2~AGE = -0.0110; 17% RSE)
    e_insbl_p2   <- -0.0150   ; label("Effect of (INS_BL - 26.98 pmol/L) on lp2")                        # Denti 2010 Table 5 (theta_P2~IBSL = -0.0150; 25% RSE)

    # ----- Inter-individual variability (Denti 2010 Table 5 covariate-model omegas) -----
    # Paper reports omega as approximate CV% on the original scale; the internal
    # variance scale is var = log(CV^2 + 1) for log-normal IIV.
    # omega_SG  = 21.0%  -> var_sg  = log(0.210^2 + 1) = 0.04313
    # omega_VOL = 10.4%  -> var_vd  = log(0.104^2 + 1) = 0.01075
    # omega_SI  = 47.5%  -> var_si  = log(0.475^2 + 1) = 0.20337
    # omega_P2  = 37.9%  -> var_p2  = log(0.379^2 + 1) = 0.13399
    # Two correlated blocks (Denti 2010 Methods p. E951):
    #   rho_SG_VOL = -0.779 -> cov = -0.779 * sqrt(0.04313 * 0.01075) = -0.01666
    #   rho_SI_P2  =  0.876 -> cov =  0.876 * sqrt(0.20337 * 0.13399) =  0.14458

    etalsg + etalvd ~ c(0.04313, -0.01666, 0.01075)   # Denti 2010 Table 5: omega_SG = 21.0% CV, omega_VOL = 10.4% CV, rho_SG_VOL = -0.779
    etalsi + etalp2 ~ c(0.20337,  0.14458, 0.13399)   # Denti 2010 Table 5: omega_SI = 47.5% CV, omega_P2 = 37.9% CV, rho_SI_P2 = 0.876

    # ----- Residual error (Denti 2010 Table 5 covariate model) -----
    # Combined proportional + additive form per Methods p. E952:
    # y_obs = g(p) * (1 + eps_prop) + eps_add (proportional and additive
    # residuals independent).
    propSd <- 0.0227   ; label("Proportional residual SD on plasma glucose (fraction)")   # Denti 2010 Table 5 (sigma_prop = 2.27%; 3% RSE)
    addSd  <- 4.28     ; label("Additive residual SD on plasma glucose (mg/dL)")          # Denti 2010 Table 5 (sigma_add = 4.28 mg/dL; 2% RSE)
  })

  model({
    # Plasma insulin time course is supplied as an error-free linear regressor.
    # Denti 2010 Methods p. E951: "I(t) is assumed as a known, error-free input
    # (forcing) function." The minimal model does not estimate insulin kinetics;
    # the user supplies the observed or simulated insulin profile at every row.
    linear(INS)

    # ----- Individual parameters -----
    # Log-normal IIV around the population typical values, with covariate effects
    # entering as additive shifts on the log scale. Pooled-population means for
    # mean-centring (Denti 2010 Table 1):
    #   AGE 55.53 yr; BODYFAT_PCT 32.39 %; FPG 91.34 mg/dL;
    #   VISCERAL_ABDOMINAL_FAT 141.8 cm^2 / CT-slice; INS_BL 26.98 pmol/L.

    vd <- exp(lvd + etalvd
              + e_sexf_vd    * SEXF
              + e_age_vd     * (AGE - 55.53)
              + e_bodyfat_vd * (BODYFAT_PCT - 32.39)
              + e_fpg_vd     * (FPG - 91.34))

    sg <- exp(lsg + etalsg)

    si <- exp(lsi + etalsi
              + e_age_si   * (AGE - 55.53)
              + e_vaf_si   * (VISCERAL_ABDOMINAL_FAT - 141.8)
              + e_insbl_si * (INS_BL - 26.98))

    p2 <- exp(lp2 + etalp2
              + e_age_p2   * (AGE - 55.53)
              + e_insbl_p2 * (INS_BL - 26.98))

    # ----- Bergman glucose minimal model -----
    # Denti 2010 Methods p. E951, equations for the glucose minimal model
    # (Bergman 1979/1981; Pacini 1986):
    #   dG/dt = -(SG + X) * G + SG * Gb            ; G(0) = Gb + D / VOL
    #   dX/dt = -P2 * X + P2 * SI * (I - Ib)        ; X(0) = 0
    #
    # State `central` carries glucose mass per body weight (mg / kg). Observation
    # Cc = central / vd is plasma glucose concentration in mg/dL. The basal-glucose
    # self-regulation arm sg * FPG * vd holds the system at Cc = FPG when
    # insulin_action = 0 (no insulin action above basal). The IVGTT bolus dose D
    # (mg / kg of glucose) is delivered by an evid=1, cmt=central row in the
    # dataset; central(0+) = FPG * vd + D reproduces the paper's
    # G(0) = Gb + D / VOL exactly.
    #
    # State `insulin_action` carries X(t) in 1/min (paper-mechanistic state;
    # declared via paper_specific_compartments above). It drives the dynamic
    # insulin-dependent glucose clearance (sg + insulin_action) * central.

    d/dt(central)        <- sg * FPG * vd - (sg + insulin_action) * central
    d/dt(insulin_action) <- -p2 * insulin_action + p2 * si * (INS - INS_BL)

    # Steady-state initial conditions before the IVGTT bolus.
    central(0)        <- FPG * vd
    insulin_action(0) <- 0

    # ----- Observation and residual error -----
    # Plasma glucose concentration in mg/dL with combined proportional + additive
    # residual error matching Denti 2010 Methods p. E952.
    Cc <- central / vd
    Cc ~ add(addSd) + prop(propSd)
  })
}
