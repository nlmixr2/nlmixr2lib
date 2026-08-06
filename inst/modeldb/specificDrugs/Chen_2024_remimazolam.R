Chen_2024_remimazolam <- function() {
  description <- "Joint population PK/PD model for remimazolam, its inactive metabolite CNS 7054, and the bispectral index (BIS) in healthy Chinese adult volunteers (Chen 2024). Remimazolam is described by a three-compartment model with first-order elimination; the whole of parent clearance feeds a single transit compartment that delays the appearance of CNS 7054, which is itself described by a two-compartment model. Sedation is described by an effect compartment equilibrating with remimazolam plasma concentration and driving an inhibitory sigmoid Imax model on BIS. Body weight enters every clearance and volume term by allometric scaling with fixed exponents (0.75 and 1) and a 60 kg reference weight; no other covariate was retained on either the PK or the PD."
  reference <- "Chen Y, Gong C, Liu F, Jiao Z, Zheng X. Toward Model-Informed Precision Dosing for Remimazolam: A Population Pharmacokinetic-Pharmacodynamic Analysis. Pharmaceutics. 2024;16(9):1122. doi:10.3390/pharmaceutics16091122"
  vignette <- "Chen_2024_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Chen 2024 Figure 1 (schematic of the
  # combined PK/PD model), which labels V1 remimazolam central, V2 "peripheral
  # shallow", V3 "peripheral deep", the transit compartment for metabolite
  # formation, V5 CNS 7054 central, V6 CNS 7054 peripheral, and V4 the effect
  # compartment.
  compartmentData <- list(
    central              = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2          = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    transit1             = list(analyte = "remimazolam-derived CNS 7054 precursor", units = "mg", specimen = "administration site", verified = TRUE),
    central_cns7054      = list(analyte = "CNS 7054", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cns7054  = list(analyte = "CNS 7054", units = "mg", specimen = "plasma", verified = TRUE),
    effect               = list(analyte = "remimazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on every clearance and volume term of both the remimazolam and the CNS 7054 disposition model, with the exponent fixed at 0.75 for clearance terms (CL, Q2, Q3, CLm, Q4) and 1 for volume terms (V1, V2, V3, V5, V6). Reference weight is 60 kg, taken from the denominator printed in Chen 2024 Equation (5); the cohort median weight is 62.5 kg (Table 2). Chen 2024 applies allometry to 'the PK models' only, so the transit rate constant Ktr and the effect-compartment rate constant ke0 -- neither of which is a clearance or a volume -- are left unscaled.",
      source_name        = "WT"
    )
  )

  # Screened by stepwise forward inclusion / backward elimination on both the
  # PK and the PD model but NOT retained in either final model ("No covariates
  # were found to affect the PK of either remimazolam or CNS 7054", Chen 2024
  # Section 3.2.1; "No covariates were identified in the PD model", Section
  # 3.2.2). No point estimates are reported for these, so they are documented
  # here rather than in covariateData.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a median-normalised power function (Chen 2024 Equation 6) on both PK and PD parameters; not retained in either final model."
    ),
    HT = list(
      description = "Body height at baseline.",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a median-normalised power function (Chen 2024 Equation 6) on both PK and PD parameters; not retained in either final model."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model (Chen 2024 Equation 7, coded 0 = male / 1 = female, which matches the SEXF polarity) on both PK and PD parameters; not retained in either final model. The authors attribute the absence of covariate effects to the homogeneity of the healthy-volunteer cohort (Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55,
    n_studies      = 2,
    age_range      = "19-43 years",
    age_median     = "28 years",
    weight_range   = "52-75 kg",
    weight_median  = "62.5 kg",
    height_range   = "151-185 cm",
    bmi_range      = "19.3-24.0 kg/m2",
    sex_female_pct = 27.3,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult volunteers (no organ disease; Mallampati score 1-2; BMI 19-24 kg/m2; weight at least 50 kg).",
    dose_range     = "Single IV bolus of 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, or 0.4 mg/kg (n = 46); or an IV bolus of 0.2 mg/kg given over 1 min followed by a 1 mg/kg/h IV infusion for 2 h (n = 9).",
    regions        = "China (single centre; conducted at Peking University First Hospital)",
    n_observations = "1113 remimazolam plasma concentrations, 1206 CNS 7054 plasma concentrations, and 1026 BIS observations. Assay calibration range 2-2000 ng/mL for both analytes.",
    notes          = "Secondary analysis of a single-centre, placebo-controlled, randomised, dose-escalation clinical pharmacology study (ChiCTR1800015185 and ChiCTR1800015186). Demographics reproduced from Chen 2024 Table 2 (Total column, N = 55); 40 men and 15 women. Sedation depth was measured as the bispectral index (BIS) by continuous EEG (BIS EEG VISTA) for 60 min after the bolus dose and for 3 h in the bolus-plus-infusion arm."
  )

  ini({
    # -----------------------------------------------------------------
    # Remimazolam (parent) disposition -- typical values at 60 kg.
    # Chen 2024 Table 3, "Typical value" block, left-hand column.
    # -----------------------------------------------------------------
    lcl  <- log(1.21);   label("Remimazolam clearance, CLp (L/min)")                          # Chen 2024 Table 3: CL for remimazolam = 1.21 L/min (RSE 3%)
    lvc  <- log(16);     label("Remimazolam central volume, V1 (L)")                          # Chen 2024 Table 3: V1 for remimazolam = 16 L (RSE 8%)
    lvp  <- log(22.6);   label("Remimazolam shallow peripheral volume, V2 (L)")               # Chen 2024 Table 3: V2 for remimazolam = 22.6 L (RSE 3%)
    lq   <- log(2.61);   label("Remimazolam intercompartmental clearance to the shallow peripheral compartment, Q2 (L/min)")  # Chen 2024 Table 3: Q2 for remimazolam = 2.61 L/min (RSE 4%)
    lvp2 <- log(23.5);   label("Remimazolam deep peripheral volume, V3 (L)")                  # Chen 2024 Table 3: V3 for remimazolam = 23.5 L (RSE 6%)
    lq2  <- log(0.227);  label("Remimazolam intercompartmental clearance to the deep peripheral compartment, Q3 (L/min)")     # Chen 2024 Table 3: Q3 for remimazolam = 0.227 L/min (RSE 14%)

    # Transit compartment delaying the appearance of CNS 7054. The whole of
    # remimazolam clearance feeds the transit compartment (Chen 2024 Figure 1
    # routes CLp from V1 into "transit", and Section 2.2.1 states that all
    # remimazolam was assumed to be converted into CNS 7054), so no separate
    # fraction-metabolised parameter is estimated.
    lktr <- log(0.447);  label("Transit rate constant for CNS 7054 formation, Ktr (1/min)")   # Chen 2024 Table 3: Ktr for remimazolam = 0.447 1/min (RSE 8%)

    # -----------------------------------------------------------------
    # CNS 7054 (inactive metabolite) disposition -- typical values at 60 kg.
    # Chen 2024 Table 3, "Typical value" block, right-hand column.
    # -----------------------------------------------------------------
    lcl_cns7054 <- log(0.0637); label("CNS 7054 clearance, CLm (L/min)")                      # Chen 2024 Table 3: CL for CNS7054 = 0.0637 L/min (RSE 3%)
    lvc_cns7054 <- log(3.72);   label("CNS 7054 central volume, V5 (L)")                      # Chen 2024 Table 3: V5 for CNS7054 = 3.72 L (RSE 5%)
    lq_cns7054  <- log(0.166);  label("CNS 7054 intercompartmental clearance, Q4 (L/min)")    # Chen 2024 Table 3: Q4 for CNS7054 = 0.166 L/min (RSE 5%)
    lvp_cns7054 <- log(5.15);   label("CNS 7054 peripheral volume, V6 (L)")                   # Chen 2024 Table 3: V6 for CNS7054 = 5.15 L (RSE 4%)

    # -----------------------------------------------------------------
    # Allometric scaling on body weight, Chen 2024 Equation (5):
    #   P_i = P_hat * (WT_i / 60)^theta1
    # with theta1 "set to" 0.75 for CL and 1 for volume -- i.e. imposed
    # rather than estimated, hence fixed() and no RSE in Table 3.
    # -----------------------------------------------------------------
    e_wt_cl_q   <- fixed(0.75); label("Allometric exponent on body weight for all clearance terms (unitless)")  # Chen 2024 Equation (5): exponent set to 0.75 for CL
    e_wt_vc_vp  <- fixed(1);    label("Allometric exponent on body weight for all volume terms (unitless)")     # Chen 2024 Equation (5): exponent set to 1 for volume

    # -----------------------------------------------------------------
    # Pharmacodynamics -- inhibitory sigmoid Imax on BIS driven by an
    # effect-compartment concentration. Chen 2024 Table 4, "Typical value".
    # -----------------------------------------------------------------
    lrbase <- log(92.5); label("Baseline bispectral index before remimazolam, BIS_baseline (BIS units)")        # Chen 2024 Table 4: BIS_baseline = 92.5 (RSE 1%)
    limax  <- log(54.5); label("Maximum reduction in BIS attributable to remimazolam, Imax (BIS units)")        # Chen 2024 Table 4: Imax = 54.5 (RSE 6%)
    lec50  <- log(504);  label("Effect-compartment concentration producing half of Imax, IC50 (ng/mL)")         # Chen 2024 Table 4: IC50 = 504 ng/mL (RSE 9%)
    lke0   <- log(1.38); label("Effect-compartment equilibration rate constant, ke0 (1/min)")                   # Chen 2024 Table 4: ke0 = 1.38 1/min (RSE 26%)
    lhill  <- log(1.44); label("Hill coefficient of the sigmoid Imax model (unitless)")                         # Chen 2024 Table 4: Hill coefficient = 1.44 (RSE 10%)

    # -----------------------------------------------------------------
    # Inter-individual variability. Chen 2024 reports BSV as a percentage
    # for an exponential BSV model (Equation 1). The percentages are read
    # as coefficients of variation and converted to the log-scale variance
    # used by ini() via omega^2 = log(CV^2 + 1); see the vignette
    # "Assumptions and deviations" section for this reading.
    # -----------------------------------------------------------------
    etalcl  ~ 0.0392207   # Chen 2024 Table 3: IIV_CL for remimazolam = 20% (RSE 10%, shrinkage 2%); log(0.20^2 + 1)
    etalvc  ~ 0.2642855   # Chen 2024 Table 3: IIV_V1 for remimazolam = 55% (RSE 10%, shrinkage 3%); log(0.55^2 + 1)
    etalq   ~ 0.0573713   # Chen 2024 Table 3: IIV_Q2 for remimazolam = 24.3% (RSE 26%, shrinkage 30%); log(0.243^2 + 1)
    etalvp  ~ 0.0980709   # Chen 2024 Table 3: IIV_V2 for remimazolam = 32.1% (RSE 13%, shrinkage 9%); log(0.321^2 + 1)
    etalktr ~ 0.1491103   # Chen 2024 Table 3: IIV_Ktr for remimazolam = 40.1% (RSE 14%, shrinkage 17%); log(0.401^2 + 1)

    etalcl_cns7054 ~ 0.0476857  # Chen 2024 Table 3: IIV_CL for CNS7054 = 22.1% (RSE 10%, shrinkage 2%); log(0.221^2 + 1)
    etalvc_cns7054 ~ 0.0867289  # Chen 2024 Table 3: IIV_V5 for CNS7054 = 30.1% (RSE 12%, shrinkage 7%); log(0.301^2 + 1)
    etalvp_cns7054 ~ 0.0301654  # Chen 2024 Table 3: IIV_V6 for CNS7054 = 17.5% (RSE 15%, shrinkage 22%); log(0.175^2 + 1)

    etalimax ~ 0.0713751  # Chen 2024 Table 4 Cont.: IIV_Imax = 27.2% (RSE 16%, shrinkage 28%); log(0.272^2 + 1)
    etalhill ~ 0.2136135  # Chen 2024 Table 4 Cont.: IIV_HILL = 48.8% (RSE 13%, shrinkage 17%); log(0.488^2 + 1)

    # -----------------------------------------------------------------
    # Residual unexplained variability. Proportional for remimazolam and
    # for BIS; combined proportional-plus-additive for CNS 7054. The
    # reported residual correlation between the two analytes (Table 3,
    # "COR (remimazolam _CNS7054)" = 0.0066) cannot be expressed in
    # nlmixr2's per-endpoint error model and is documented in the vignette.
    # -----------------------------------------------------------------
    propSd          <- 0.232;  label("Proportional residual error on remimazolam (fraction)")   # Chen 2024 Table 3: prop RUV for remimazolam = 23.2% (RSE 1%, shrinkage 7%)
    propSd_cns7054  <- 0.064;  label("Proportional residual error on CNS 7054 (fraction)")      # Chen 2024 Table 3: prop RUV for CNS7054 = 6.4% (RSE 0%, shrinkage 14%)
    addSd_cns7054   <- 43.13;  label("Additive residual error on CNS 7054 (ng/mL)")             # Chen 2024 Table 3: add RUV for CNS 7054 = 43.13 ng/mL (RSE 8%, shrinkage 8%)
    propSd_BIS      <- 0.112;  label("Proportional residual error on BIS (fraction)")           # Chen 2024 Table 4 Cont.: prop RUV pd = 11.2% (RSE 2.4%, shrinkage 4%)
  })

  model({
    # Reference body weight -- the denominator printed in Chen 2024
    # Equation (5), not the cohort median (62.5 kg, Table 2).
    ref_wt <- 60

    # Amounts are carried in mg and volumes in L, so central/vc is mg/L;
    # Chen 2024 reports every concentration (IC50, the CNS 7054 additive
    # residual error, and the 2-2000 ng/mL assay calibration range) in
    # ng/mL, so concentrations are converted with 1 mg/L = 1000 ng/mL.
    mgL_to_ngmL <- 1000

    # Allometric multipliers, Chen 2024 Equation (5).
    wt_cl <- (WT / ref_wt)^e_wt_cl_q
    wt_v  <- (WT / ref_wt)^e_wt_vc_vp

    # Remimazolam (parent) individual parameters.
    cl  <- exp(lcl  + etalcl)  * wt_cl
    vc  <- exp(lvc  + etalvc)  * wt_v
    q   <- exp(lq   + etalq)   * wt_cl
    vp  <- exp(lvp  + etalvp)  * wt_v
    q2  <- exp(lq2)            * wt_cl
    vp2 <- exp(lvp2)           * wt_v
    ktr <- exp(lktr + etalktr)

    # CNS 7054 (metabolite) individual parameters.
    cl_cns7054 <- exp(lcl_cns7054 + etalcl_cns7054) * wt_cl
    vc_cns7054 <- exp(lvc_cns7054 + etalvc_cns7054) * wt_v
    q_cns7054  <- exp(lq_cns7054)                   * wt_cl
    vp_cns7054 <- exp(lvp_cns7054 + etalvp_cns7054) * wt_v

    # Pharmacodynamic individual parameters.
    ke0   <- exp(lke0)
    imax  <- exp(limax + etalimax)
    ec50  <- exp(lec50)
    hill  <- exp(lhill + etalhill)
    rbase <- exp(lrbase)

    # Micro-constants.
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    k13  <- q2 / vc
    k31  <- q2 / vp2
    kelm <- cl_cns7054 / vc_cns7054
    k56  <- q_cns7054  / vc_cns7054
    k65  <- q_cns7054  / vp_cns7054

    # Three-compartment remimazolam disposition with IV input into the
    # central compartment (Chen 2024 Figure 1).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Metabolite formation. All remimazolam clearance is assumed to be
    # conversion to CNS 7054 (Chen 2024 Section 2.2.1), delayed by one
    # transit compartment. Chen 2024 states no stoichiometric conversion
    # factor, so the amount transfer is 1:1 and the fitted CNS 7054
    # volumes and clearance are apparent values on a remimazolam-mass
    # basis; see the vignette "Assumptions and deviations".
    d/dt(transit1)            <-  kel * central - ktr * transit1
    d/dt(central_cns7054)     <-  ktr * transit1 - kelm * central_cns7054 -
                                  k56 * central_cns7054 + k65 * peripheral1_cns7054
    d/dt(peripheral1_cns7054) <-  k56 * central_cns7054 - k65 * peripheral1_cns7054

    # Plasma concentrations (ng/mL).
    Cc         <- central         / vc         * mgL_to_ngmL
    Cc_cns7054 <- central_cns7054 / vc_cns7054 * mgL_to_ngmL

    # Effect compartment. Chen 2024 Figure 1 carries k1e into and ke0 out
    # of the effect compartment and Table 4 reports only ke0, so the
    # standard massless effect compartment (k1e = ke0) applies. The state
    # is a concentration in ng/mL, not an amount.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Inhibitory sigmoid Imax model on BIS, Chen 2024 Equation (8). The
    # published denominator reads "IC50^Hill x CE^Hill"; the multiplication
    # is a typesetting error for the sum, because the printed form is
    # independent of CE and would contradict the paper's own definition of
    # IC50 as "the concentration at half-maximum effect". See the vignette
    # "Assumptions and deviations".
    BIS <- rbase - imax * effect^hill / (ec50^hill + effect^hill)

    Cc         ~ prop(propSd)
    Cc_cns7054 ~ prop(propSd_cns7054) + add(addSd_cns7054)
    BIS        ~ prop(propSd_BIS)
  })
}
