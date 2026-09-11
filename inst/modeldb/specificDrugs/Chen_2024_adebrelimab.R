Chen_2024_adebrelimab <- function() {
  description <- "Two-compartment population PK model for adebrelimab (anti-PD-L1 IgG4) with empirical sigmoid time-varying clearance, fitted to 263 Chinese patients with extensive-stage small-cell lung cancer or advanced solid tumours from the phase I SHR-1316-I-101 and phase III CAPSTONE-1 (SHR-1316-III-301) studies (Chen 2024)"
  reference <- "Chen P, Zhang Y, Wang Y, Ma K, Shi W, Djebli N, Shen K. Population pharmacokinetics of adebrelimab - Support of alternative flat dose regimen in extensive-stage small-cell lung cancer. CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1238-1251. doi:10.1002/psp4.13155"
  vignette <- "Chen_2024_adebrelimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Adebrelimab was quantified by ELISA in human SERUM
  # (Chen 2024 Methods, "Pharmacokinetic sampling and bioanalytical methods":
  # "The lower limit of quantification (LLOQ) of the adebrelimab assay in
  # human serum was 0.4 ug/mL"), so the specimen is serum rather than plasma.
  compartmentData <- list(
    central     = list(analyte = "adebrelimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "adebrelimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (not time-varying) body weight. Power scaling (WT / 64)^exponent on CL (0.710), on Vc (0.577) and on Vp (1.80). Reference 64 kg is the population median body weight (Chen 2024 Table 1: 64.0 kg, range 38.1-97.0 kg) and is the value hard-coded in the Data S1 NONMEM control stream (CLBW / V1BW / V2BW blocks all divide by 64).",
      source_name        = "BW"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum albumin in g/L, matching the ALB canonical SI unit -- no unit conversion is required for this model. Power scaling (ALB / 41.4)^-0.861 on CL. Reference 41.4 g/L is the total-population median (Chen 2024 Table 1) and is hard-coded in the Data S1 control stream (CLALB = (ALB/41.4)**THETA(9)). Chen 2024 Figure 1 quotes the 5th and 95th percentiles as 33.06 and 48.29 g/L.",
      source_name        = "ALB"
    ),
    NEUT = list(
      description        = "Baseline absolute neutrophil count",
      units              = "cells/mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (NEUT / 4150)^0.159 on CL. Chen 2024 reports neutrophil count in 10^9/L (Table 1 median 4.2 x 10^9/L; Figure 1 5th and 95th percentiles 2.381 and 7.781 x 10^9/L) and the Data S1 control stream normalises by 4.15 x 10^9/L. This model file carries the covariate in the NEUT canonical unit of cells/mm^3, so the reference is written as 4150 cells/mm^3 (= 4.15 x 10^9/L, since 1 L = 10^6 mm^3). The effect enters only as the ratio (NEUT / reference), so the exponent is numerically identical under either unit provided the data column and the reference share a unit.",
      source_name        = "NEUT"
    ),
    TUM_SLD = list(
      description        = "Baseline sum of the longest diameters of all target lesions (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (TUM_SLD / 90)^0.103 on CL. Reference 90 is hard-coded in the Data S1 control stream (CLSLD = (SLD/90)**THETA(12)); Chen 2024 Table 1 does not tabulate SLD, so the control stream is the only source for the normaliser. The unit is mm, established from Chen 2024 Figure 1, which labels the SLD 5th and 95th percentiles as 28.03 mm and 176.75 mm -- a reference of 90 mm sits between them as the population median. The effect enters only as the ratio (TUM_SLD / reference).",
      source_name        = "SLD"
    ),
    ADA_POS = list(
      description        = "Treatment-emergent anti-drug antibody positivity, at the subject level",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Subject-level (not time-varying) indicator: 1 for the whole treatment period if the patient experienced treatment-induced or treatment-boosted ADA positivity at any point, 0 otherwise (Chen 2024 Table 2 footnote a). Applied as the linear form (1 + 0.185 * ADA_POS) per Chen 2024 Equation 4 and the Data S1 control stream (IF(ADA.EQ.0) CLADA=1; IF(ADA.EQ.1) CLADA=(1+THETA(8))), so ADA-positive patients have 18.5% higher CL. 27.0% of the analysis population were ADA-positive (Chen 2024 Table 1).",
      source_name        = "ADA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 263L,
    n_studies      = 2L,
    age_range      = "18-73 years",
    age_median     = "61 years",
    weight_range   = "38.1-97.0 kg",
    weight_median  = "64.0 kg",
    sex_female_pct = 25.1,
    race_ethnicity = c(Asian_Chinese = 100),
    disease_state  = "Extensive-stage small-cell lung cancer (222 patients from the phase III CAPSTONE-1 study, all receiving adebrelimab with carboplatin and etoposide) pooled with 41 patients with advanced solid tumours from the phase I dose-escalation/expansion study. 87.5% lung cancer; 100% metastatic; 96.2% clinical stage IV; 87.5% ECOG performance status 1.",
    dose_range     = "3 mg/kg q3w (n = 3), 10 mg/kg q2w (n = 12), 10 mg/kg q3w (n = 13), and 20 mg/kg q3w (n = 235) as intravenous infusions",
    regions        = "China",
    ada_positive_pct = 27.0,
    albumin_median   = "41.4 g/L (range 26.6-53.2)",
    neutrophil_median = "4.2 x 10^9/L (range 1.7-11.4)",
    hepatic_impairment_pct = c(none = 84.0, mild = 15.2, moderate = 0.8, severe = 0),
    renal_impairment_pct   = c(none = 51.0, mild = 39.5, moderate = 9.5, severe = 0),
    notes          = "Baseline demographics from Chen 2024 Table 1 (total population column, N = 263). Studies: NCT03474289 (SHR-1316-I-101, phase I, advanced tumours) and NCT03711305 (SHR-1316-III-301, CAPSTONE-1, phase III, ES-SCLC). Reference patient for the structural parameters is the covariate median: 64 kg, albumin 41.4 g/L, neutrophils 4.15 x 10^9/L, SLD 90 mm, ADA-negative."
  )

  ini({
    # Structural PK parameters for the reference (median-covariate,
    # ADA-negative) patient at time zero. Values are the final-model point
    # estimates from the Data S1 NONMEM control stream $THETA block, which
    # carries more significant figures than the rounded Chen 2024 Table 2
    # values quoted in each trailing comment. CL and Q are already reported
    # in L/day, so the day time unit needs no conversion.
    lcl <- log(0.237784); label("Baseline clearance at t = 0 for the reference patient (L/day)")      # Chen 2024 Data S1 $THETA(1); Table 2 final model CL = 0.238 L/day
    lvc <- log(3.22559);  label("Central volume of distribution for the reference patient (L)")       # Chen 2024 Data S1 $THETA(2); Table 2 final model V1 = 3.23 L
    lq  <- log(0.703417); label("Intercompartmental clearance (L/day)")                                # Chen 2024 Data S1 $THETA(3); Table 2 final model Q = 0.703 L/day
    lvp <- log(1.69298);  label("Peripheral volume of distribution for the reference patient (L)")    # Chen 2024 Data S1 $THETA(4); Table 2 final model V2 = 1.69 L

    # Empirical time-varying clearance (sigmoid maximal-change function of
    # time since the first dose; Chen 2024 Equation 1):
    #   CL(t) = TVCL * exp(Imax * t^HILL / (TC50^HILL + t^HILL))
    # Imax is the maximal change in log-CL and is NEGATIVE, i.e. clearance
    # falls with time on treatment. Unlike the log-normal |Imax| idiom used by
    # Kuchimanchi_2024_dostarlimab, Chen 2024 gives Imax a NORMAL distribution
    # and an ADDITIVE eta (Table 2 footnote b, "The additive model for
    # interindividual variability of Imax was used"; Data S1 IMAX =
    # TVIMAX + ETA(4)), so Imax is carried here on the natural scale and an
    # individual Imax may cross zero. That is the authors' structure and is
    # reproduced faithfully. The asymptotic CL reduction at t >> TC50 is
    # 1 - exp(-0.34925) = 29.5%, consistent with the reported drop in
    # geometric-mean CL from 0.25 L/day at baseline to 0.177 L/day at steady
    # state (ratio 0.708 = exp(-0.345)).
    cl_time_max   <- -0.34925;     label("Imax; maximal change in log-CL at t >> TC50 (unitless)")   # Chen 2024 Data S1 $THETA(5); Table 2 final model Imax = -0.349
    lcl_t50       <- log(74.4482); label("log TC50; time at half of the maximal change in CL (days)") # Chen 2024 Data S1 $THETA(6); Table 2 final model TC50 = 74.4 days
    lcl_time_hill <- log(1.96911); label("log HILL; sigmoid shape coefficient in time (unitless)")    # Chen 2024 Data S1 $THETA(7); Table 2 final model HILL = 1.97

    # Covariate effects. Continuous covariates use the power form
    # (Cov / Cov_median)^theta (Chen 2024 Equation 3) and ADA uses the linear
    # form (1 + theta) for the positive group (Chen 2024 Equation 4). The
    # median normalisers are hard-coded in the Data S1 control stream:
    # BW/64, ALB/41.4, NEUT/4.15 (10^9/L, i.e. 4150 cells/mm^3), SLD/90.
    # All five act multiplicatively on CL (Data S1: CLCOV =
    # CLADA*CLALB*CLBW*CLNEUT*CLSLD); body weight additionally acts on both
    # volumes with two separately estimated exponents.
    e_ada_pos_cl <- 0.185112;  label("ADA-positive fractional change on CL vs the ADA-negative reference (fraction)")  # Chen 2024 Data S1 $THETA(8); Table 2 final model ADA on CL = 0.185
    e_alb_cl     <- -0.861333; label("ALB power exponent on CL, ALB/41.4 scaling (unitless)")                  # Chen 2024 Data S1 $THETA(9); Table 2 final model Albumin on CL = -0.861
    e_wt_cl      <- 0.70974;   label("WT power exponent on CL, WT/64 scaling (unitless)")                      # Chen 2024 Data S1 $THETA(10); Table 2 final model Body weight on CL = 0.710
    e_neut_cl    <- 0.158825;  label("NEUT power exponent on CL, NEUT/4150 scaling (unitless)")                # Chen 2024 Data S1 $THETA(11); Table 2 final model Neutrophil count on CL = 0.159
    e_tum_sld_cl <- 0.102909;  label("TUM_SLD power exponent on CL, TUM_SLD/90 scaling (unitless)")            # Chen 2024 Data S1 $THETA(12); Table 2 final model SLD on CL = 0.103
    e_wt_vc      <- 0.576703;  label("WT power exponent on Vc, WT/64 scaling (unitless)")                      # Chen 2024 Data S1 $THETA(13); Table 2 final model Body weight on V1 = 0.577
    e_wt_vp      <- 1.79986;   label("WT power exponent on Vp, WT/64 scaling (unitless)")                      # Chen 2024 Data S1 $THETA(14); Table 2 final model Body weight on V2 = 1.8

    # IIV. Chen 2024 estimated a full 3x3 covariance block across CL, V1 and
    # V2 (Table 2 "Cov" rows; Data S1 $OMEGA BLOCK(3)) plus an independent
    # eta on Imax. The block is given in NONMEM lower-triangular row-major
    # order, which is the same order nlmixr2 expects:
    #   var(CL); cov(CL,V1) var(V1); cov(CL,V2) cov(V1,V2) var(V2).
    # Q carries no IIV in the final model (Data S1: Q = TVQ).
    etalcl + etalvc + etalvp ~ c(0.0415871,
                                 0.0104212, 0.0296031,
                                 0.0370675, 0.0308242, 0.128194)  # Chen 2024 Data S1 $OMEGA BLOCK(3); Table 2: omega^2_CL 0.0416, Cov(CL,V1) 0.0104, omega^2_V1 0.0296, Cov(CL,V2) 0.0371, Cov(V1,V2) 0.0308, omega^2_V2 0.128
    # Additive eta on the natural-scale Imax (see cl_time_max above).
    etacl_time_max ~ 0.0533643  # Chen 2024 Data S1 $OMEGA (second block); Table 2: omega^2_Imax = 0.0534

    # Residual error, combined proportional and additive (Data S1 $ERROR:
    # Y = F*(1+EPS(1)) + EPS(2)). Chen 2024 Table 2 and the $SIGMA block
    # report VARIANCES; nlmixr2 wants standard deviations, so each value is
    # the square root of the reported sigma^2.
    propSd <- 0.149617; label("Proportional residual error SD (fraction)")  # sqrt(0.0223852); Chen 2024 Data S1 $SIGMA(1), Table 2 sigma^2_Prop = 0.0224
    addSd  <- 5.71370;  label("Additive residual error SD (ug/mL)")         # sqrt(32.6463);  Chen 2024 Data S1 $SIGMA(2), Table 2 sigma^2_add = 32.6 (ug/mL)^2
  })

  model({
    # 1. Covariate effects. Continuous covariates enter as the power of a
    #    median-normalised ratio; ADA enters linearly on the positive group.
    ada_cl     <- 1 + e_ada_pos_cl * ADA_POS
    alb_cl     <- (ALB / 41.4)^e_alb_cl
    wt_cl      <- (WT / 64)^e_wt_cl
    neut_cl    <- (NEUT / 4150)^e_neut_cl
    tum_sld_cl <- (TUM_SLD / 90)^e_tum_sld_cl

    wt_vc <- (WT / 64)^e_wt_vc
    wt_vp <- (WT / 64)^e_wt_vp

    # 2. Time-varying clearance multiplier (Chen 2024 Equation 1). t is time
    #    since the first dose, in days. At t = 0 the multiplier is exactly 1,
    #    so lcl is the baseline clearance; it decays towards exp(Imax) as
    #    t grows past TC50. The eta on Imax is additive, matching the paper's
    #    normal distribution for that parameter.
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_t50        <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)
    td_cl <- exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # 3. Individual PK parameters.
    cl_base <- exp(lcl + etalcl) * ada_cl * alb_cl * wt_cl * neut_cl * tum_sld_cl
    cl      <- cl_base * td_cl
    vc      <- exp(lvc + etalvc) * wt_vc
    vp      <- exp(lvp + etalvp) * wt_vp
    q       <- exp(lq)

    # 4. Two-compartment micro-constants (NONMEM ADVAN3 TRANS4).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation. Dose in mg / volume in L = mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
