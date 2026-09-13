Assmus_2026_dndi6148 <- function() {
  description <- paste(
    "One-compartment population PK model with first-order oral absorption for the",
    "benzoxaborole antileishmanial / antichagasic DNDI-6148 in healthy adult men,",
    "fitted to the first-in-human single-ascending-dose study (10-380 mg, eight",
    "cohorts, 48 subjects). The drug's non-linear pharmacokinetics are carried by",
    "two exponential effects of the weight-normalised dose, centred on the study",
    "median of 1.7 mg/kg: one on relative oral bioavailability and one on apparent",
    "clearance. The two act in opposite directions on exposure, which is why Cmax",
    "rises less than dose-proportionally while AUCinf stays approximately",
    "dose-linear. Allometric body-weight scaling (0.75 on clearance, 1 on volume,",
    "both fixed) sits underneath. Inter-individual variability on apparent volume",
    "was estimated near zero and fixed to zero by the authors, so no eta is carried",
    "on it here. See modellib('Henninger_2026_dndi6148_mouse') and",
    "modellib('Henninger_2026_dndi6148_human') for the independent murine",
    "target-site PK/PD model of the same molecule and its allometric human",
    "projection; this model is the first fit to actual human data."
  )
  reference <- paste(
    "Assmus F, Adehin A, Hoglund RM, Mowbray CE, Gillon JY, Blesson S,",
    "Braillard S, Chatelain E, Scandale I, Tarning J. (2026). Population",
    "pharmacokinetics of DNDI-6148 in healthy adults. PLoS Negl Trop Dis",
    "20(4):e0014220. doi:10.1371/journal.pntd.0014220.",
    "Final NONMEM control stream: S1 Code of that paper."
  )
  vignette <- "Assmus_2026_dndi6148"

  # Dose records carry milligrams; central/vc is therefore mg/L == ug/mL, which
  # is the unit Table 3 uses for AUCinf (ug x h/mL). Table 3 reports Cmax in
  # ng/mL, i.e. 1000 x the model's Cc.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on apparent clearance and apparent volume, standardised",
        "to 70 kg with exponents fixed a priori at 0.75 and 1 (Methods,",
        "'(iii) Covariate model'; S1 Code $PK, '((WT/70)**0.75)' and",
        "'((WT/70)**1.00)'). Time-fixed: this is a single-dose study with baseline",
        "weight only. Median 72.3 kg, range 56.9-96.5 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    DOSE_DNDI6148_MGKG = list(
      description        = "Administered DNDI-6148 dose level, weight-normalised",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a CENTRED exponential term on both relative bioavailability and",
        "apparent clearance, centred on the study median dose of 1.7 mg/kg (Results,",
        "'Population pharmacokinetic model'; Table 2 footnote c; S1 Code $PK,",
        "'EXP(THETA(5)*(DOSE_P_KG - 1.7))' and 'EXP(THETA(6)*(DOSE_P_KG - 1.7))').",
        "Both effects are therefore exactly 1 at 1.7 mg/kg. Required as a data",
        "column because rxode2 model code cannot read the amt of the dose record it",
        "is scaling; the absolute dose in mg that amt carries is",
        "DOSE_DNDI6148_MGKG * WT. Time-fixed per subject in this single-dose study.",
        "Pooled median 1.70 mg/kg, range 0.12-5.56 mg/kg (Table 1)."
      ),
      source_name        = "DOSE_P_KG"
    )
  )

  # Covariates screened by the stepwise covariate model (Methods '(iii) Covariate
  # model') but not retained after backward elimination (Results, 'No other
  # covariates tested during model development were retained after backward
  # elimination'). Baseline values for all of these are tabulated in Table 1; no
  # point estimate is reported for any of them because none entered the final
  # model, so they are documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Median 35 years, range 18-50 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Liver-function marker screened in the SCM; not retained. Median 19 IU/L (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Liver-function marker screened in the SCM; not retained. Median 19 IU/L (Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Liver-function marker screened in the SCM; not retained. Median 66.5 IU/L (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Liver-function marker screened in the SCM; not retained. Median 10 umol/L (Table 1)."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Liver-function marker screened in the SCM; not retained. Median 15.5 IU/L (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Kidney-function marker screened in the SCM; not retained. Median 113 mL/min,",
        "range 81.5-161 (Table 1). Consistent with the Discussion, which reports",
        "< 0.2% of unchanged drug recovered in urine."
      )
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Median 43.0%, range 37.1-48.4 (Table 1)."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "DNDI-6148", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "DNDI-6148", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 48,
    n_studies      = 1,
    age_range      = "18-50 years",
    age_median     = "35 years",
    weight_range   = "56.9-96.5 kg",
    weight_median  = "72.3 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = "healthy volunteers",
    dose_range     = paste(
      "10-380 mg single oral dose (free acid equivalent), eight ascending cohorts of",
      "6 active subjects each: 10, 20, 40, 80, 160, 220, 300, 380 mg"
    ),
    regions        = "France (single centre, Gieres)",
    formulation    = paste(
      "DNDI-6148 arginine monohydrate powder for suspension, reconstituted in",
      "ORA-Sweet vehicle, administered under fasting conditions"
    ),
    notes          = paste(
      "Phase 1 first-in-human single ascending dose study (EudraCT 2018-004023-37;",
      "ISRCTN54981564), 2018-2022. 64 healthy White men enrolled (48 active,",
      "16 placebo); the PK analysis used the 48 active subjects contributing 684",
      "plasma samples, all post-dose samples above the 1 ng/mL LLOQ. Baseline",
      "demographics and laboratory characteristics are in Table 1."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Table 2 reports these for 'an adult weighting
    # 70 kg' and, through the centred dose terms below, at the median dose of
    # 1.7 mg/kg. The S1 Code $THETA block carries the identical numbers
    # (2.55, 69.9, 0.576, 1 FIX, -0.123, -0.15), so the stream reproduces the
    # FINAL estimates rather than initial values.
    # -----------------------------------------------------------------------
    lka <- log(0.576)
    label("Absorption rate constant (1/h)")  # Table 2 'Absorption rate constant, K A (h-1)' = 0.576 (7.7% RSE); S1 Code $THETA 3

    lcl <- log(2.55)
    label("Apparent clearance at 70 kg and 1.7 mg/kg (L/h)")  # Table 2 'Apparent clearance, CL/F (L/h)' = 2.55 (5.9% RSE); S1 Code $THETA 1

    lvc <- log(69.9)
    label("Apparent central volume of distribution at 70 kg (L)")  # Table 2 'Apparent volume of distribution, V/F (L)' = 69.9 (2.8% RSE); S1 Code $THETA 2

    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability at 1.7 mg/kg (fraction)")  # Table 2 'Relative oral bioavailability, F' = '1 fixed'; S1 Code $THETA 4 '(1) FIX'

    # -----------------------------------------------------------------------
    # Allometric exponents, fixed a priori (not estimated): Methods
    # '(iii) Covariate model' -- 'Body weight (standardized to 70 kg) was
    # implemented a priori as an allometric function on all clearance
    # (exponent 0.75) and volume (exponent 1) parameters.' Neither appears in
    # the S1 Code $THETA block; both are hardcoded in $PK.
    # -----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on apparent clearance (unitless)")  # Methods '(iii) Covariate model'; S1 Code $PK '((WT/70)**0.75)'

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on apparent volume (unitless)")  # Methods '(iii) Covariate model'; S1 Code $PK '((WT/70)**1.00)'

    # -----------------------------------------------------------------------
    # Dose effects. Exponential functions centred on the median dose of
    # 1.7 mg/kg were selected over power functions, which fitted better but
    # 'produced steep relationships at low doses' (Results). The coefficient
    # multiplies a dose difference in mg/kg, so its unit is kg/mg.
    # -----------------------------------------------------------------------
    e_dose_fdepot <- -0.123
    label("Exponential effect of weight-normalised dose on relative bioavailability (kg/mg)")  # Table 2 'theta Dose_F (dose effect on F)' = -0.123 (13.5% RSE); S1 Code $THETA 5

    e_dose_cl <- -0.150
    label("Exponential effect of weight-normalised dose on apparent clearance (kg/mg)")  # Table 2 'theta Dose_CL (dose effect on CL/F)' = -0.150 (12.8% RSE); S1 Code $THETA 6 '-0.15'

    # -----------------------------------------------------------------------
    # Inter-individual variability. S1 Code $OMEGA carries the variances
    # directly; Table 2 reports the same quantities as %CV via its footnote a,
    # CV% = 100 x sqrt(exp(omega^2) - 1), which reproduces each printed value:
    #   0.0905 -> 30.8% (Table 2 CL/F IIV 30.8)
    #   0.274  -> 56.1% (Table 2 K A  IIV 56.1)
    #   0.0284 -> 17.0% (Table 2 F    IIV 17.0)
    # The fourth omega, on apparent volume, is '0 FIX' in S1 Code and has no
    # IIV row in Table 2 -- Results: 'Inter-individual variability on apparent
    # volume of distribution was estimated to be close to zero and was
    # therefore fixed to zero in the final model.' It is OMITTED here rather
    # than written as `etalvc ~ fixed(0)`, because a zero-variance diagonal
    # makes OMEGA singular and breaks the Cholesky sampler rxSolve uses.
    # -----------------------------------------------------------------------
    etalcl     ~ 0.0905  # S1 Code '$OMEGA 0.0905; 1 IIV_CL/F'; Table 2 IIV 30.8% CV (12.8% RSE)
    etalka     ~ 0.274   # S1 Code '$OMEGA 0.274; 3 IIV_KA';    Table 2 IIV 56.1% CV (10.9% RSE)
    etalfdepot ~ 0.0284  # S1 Code '$OMEGA 0.0284; 4 IIV_F';    Table 2 IIV 17.0% CV (13.6% RSE)

    # -----------------------------------------------------------------------
    # Residual unexplained variability. Methods '(ii) Structural and stochastic
    # model': 'Residual unexplained variability was implemented as an additive
    # error on log-transformed concentrations, equivalent to an exponential
    # error on the arithmetic scale.' S1 Code $ERROR confirms it:
    # IPRED = LOG(CP); Y = IPRED + EPS(1). That maps to `~ lnorm(expSd)` with
    # expSd the log-scale SD, i.e. sqrt of the reported $SIGMA VARIANCE:
    # sqrt(0.0361) = 0.19 exactly.
    # -----------------------------------------------------------------------
    expSd <- 0.19
    label("Log-normal residual error, SD on the natural-log scale (unitless)")  # Table 2 'Variance of residual error, sigma' = 0.0361 (13.4% RSE); S1 Code '$SIGMA 0.0361'; SD = sqrt(0.0361)
  })

  model({
    # ---- 1. Derived covariate terms --------------------------------------
    # Table 2 footnote c:
    #   F_i    = F    * exp(theta_Dose_F  * (Dose_i - Dose_median))
    #   CL_i/F = CL/F * exp(theta_Dose_CL * (Dose_i - Dose_median))
    # with Dose_median = 1.7 mg/kg (Results, 'Exponential functions centered on
    # the median dose (Dose_median = 1.7 mg/kg)'). S1 Code $PK writes the same
    # two terms as COV and COV2.
    doseEffectFdepot <- exp(e_dose_fdepot * (DOSE_DNDI6148_MGKG - 1.7))
    doseEffectCl     <- exp(e_dose_cl     * (DOSE_DNDI6148_MGKG - 1.7))

    # Reference weight 70 kg: Methods '(iii) Covariate model'.
    wtRatio <- WT / 70

    # ---- 2. Individual parameters ----------------------------------------
    # S1 Code $PK:
    #   TVCL = THETA(1) * ((WT/70)**0.75) * COV2 ; CL = TVCL * EXP(ETA(1))
    #   TVV2 = THETA(2) * ((WT/70)**1.00)        ; V2 = TVV2 * EXP(ETA(2))
    #   TVKA = THETA(3)                          ; KA = TVKA * EXP(ETA(3))
    #   TVF1 = THETA(4) * COV                    ; F1 = TVF1 * EXP(ETA(4))
    # ETA(2) is dropped here because its omega is 0 FIX (see ini()).
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl) * wtRatio^e_wt_cl * doseEffectCl
    vc     <- exp(lvc)          * wtRatio^e_wt_vc
    fdepot <- exp(lfdepot + etalfdepot) * doseEffectFdepot

    # ---- 3. Micro-constants ----------------------------------------------
    # S1 Code $PK: K12 = KA; K20 = CL/V2; S2 = V2.
    kel <- cl / vc

    # ---- 4. ODE system ----------------------------------------------------
    # S1 Code: $SUBROUTINE ADVAN5 TRANS1 with COMP=(1) absorption and
    # COMP=(2) central, linked by K12 = KA and eliminated by K20.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 5. Bioavailability -----------------------------------------------
    # F1 in NONMEM scales the dose into the absorption compartment.
    f(depot) <- fdepot

    # ---- 6. Observation and error -----------------------------------------
    # S1 Code $ERROR: CP = A(2)/S2 with S2 = V2.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
