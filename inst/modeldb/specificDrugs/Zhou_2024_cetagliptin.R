Zhou_2024_cetagliptin <- function() {
  description <- "Two-compartment population PK model with first-order absorption and saturable Michaelis-Menten elimination for cetagliptin, coupled by a direct-effect sigmoid Emax model to plasma DPP-4 inhibition, in Chinese patients with type 2 diabetes mellitus. Total bilirubin is a power covariate on the peripheral volume of distribution."
  reference   <- "Zhou C, Zhou S, Wang J, Xie L, Lv Z, Zhao Y, Wang L, Luo H, Xie D, Shao F. Safety, tolerability, pharmacokinetics and pharmacokinetic-pharmacodynamic modeling of cetagliptin in patients with type 2 diabetes mellitus. Front Endocrinol (Lausanne). 2024;15:1359407. doi:10.3389/fendo.2024.1359407"
  vignette    <- "Zhou_2024_cetagliptin"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the Zhou 2024 final population PK model.",
        "Supplementary Table 1 reports the coefficient dV2dTBIL = 0.3723 (RSE 3.066%)",
        "under the Phoenix NLME naming convention dPARAMdCOVARIATE, and Discussion",
        "page 10 states 'only TBIL had a significant effect on V2 ... V2 increases with",
        "the increase of TBIL'. Phoenix NLME's default continuous-covariate",
        "transformation is median-centred log-ratio,",
        "V2 = tvV2 * exp(dV2dTBIL * log(TBIL / median(TBIL))), which is algebraically",
        "the power form V2 = tvV2 * (TBIL / ref)^0.3723 encoded here. The paper reports",
        "NEITHER the covariate equation NOR the population median TBIL (total bilirubin",
        "does not appear in the Table 1 demographics), so the centring value is an",
        "assumption: 10 umol/L, a rounded mid-normal adult total bilirubin, chosen so",
        "that the Supplementary Table 1 typical value tvV2 = 558.0 L is the typical",
        "peripheral volume of the study population. Results page 6 states that no",
        "subject had clinically significant abnormal liver function, so all observed",
        "TBILI values lay inside the normal 5-21 umol/L range. See vignette Errata.",
        "Chinese clinical-chemistry laboratories report total bilirubin in umol/L, which",
        "is also the canonical unit for this column."
      ),
      source_name        = "TBIL"
    )
  )

  # Covariates the paper screened but did not retain. Discussion page 10 lists
  # the full stepwise forward-inclusion / backward-elimination candidate set:
  # "we investigated the effects of gender, body weight, glutamic-pyruvic
  # transaminase, total bilirubin, triglyceride, low-density lipoprotein
  # cholesterol, glucose, urea, and creatinine on pharmacokinetic parameters,
  # which finally proved that only TBIL had a significant effect on V2". Only
  # TBIL survived, so the other eight are recorded here as documentation of the
  # covariate screen; no effect size is reported for any of them anywhere in the
  # paper or supplement, and none is referenced in model().
  screened_not_retained <- paste(
    "Screened in the Zhou 2024 stepwise covariate search (Methods section",
    "2.7.1; candidate set listed in Discussion page 10) but not retained in the",
    "final population PK model. The paper reports no point estimate, standard",
    "error or confidence interval for this effect, so nothing is encodable."
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female (screened, not retained)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(screened_not_retained, "Source term: 'gender'."),
      source_name        = "gender"
    ),
    WT = list(
      description        = "Body weight (screened, not retained)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(screened_not_retained, "Source term: 'body weight'."),
      source_name        = "body weight"
    ),
    ALT = list(
      description        = "Alanine aminotransferase (screened, not retained)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'glutamic-pyruvic transaminase', the older name for ALT."
      ),
      source_name        = "glutamic-pyruvic transaminase"
    ),
    TRIG = list(
      description        = "Triglycerides (screened, not retained)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'triglyceride'. Units assumed mmol/L, the SI convention",
        "used by Chinese clinical-chemistry laboratories; the paper does not",
        "state the unit because no effect was retained."
      ),
      source_name        = "triglyceride"
    ),
    LDLC = list(
      description        = "Low-density lipoprotein cholesterol (screened, not retained)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'low-density lipoprotein cholesterol'. Units assumed",
        "mmol/L per the SI convention; not stated in the paper."
      ),
      source_name        = "low-density lipoprotein cholesterol"
    ),
    FPG = list(
      description        = "Fasting plasma glucose at baseline (screened, not retained)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'glucose'. Recorded as FPG rather than the time-varying",
        "regressor canonical GLU because the study measured fasting plasma",
        "glucose as a baseline clinical-chemistry value (Table 1 reports",
        "baseline FPG in mmol/L), not as a within-subject glucose input."
      ),
      source_name        = "glucose"
    ),
    BUN = list(
      description        = "Blood urea (screened, not retained)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'urea'. Chinese clinical-chemistry laboratories report",
        "urea in mmol/L rather than as mg/dL of urea nitrogen; 1 mmol/L urea",
        "is about 2.80 mg/dL BUN."
      ),
      source_name        = "urea"
    ),
    CREAT = list(
      description        = "Serum creatinine (screened, not retained)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        screened_not_retained,
        "Source term: 'creatinine'. Recorded as CREAT (serum creatinine)",
        "rather than CRCL: the paper screened the measured serum concentration,",
        "not a derived creatinine clearance. Units assumed umol/L per the SI",
        "convention; not stated in the paper."
      ),
      source_name        = "creatinine"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "cetagliptin", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "cetagliptin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cetagliptin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = 560L,
    age_range      = "18-65 years (protocol inclusion criterion)",
    age_median     = "mean 47.80 (SD 4.32) years in the 50 mg arm and 45.20 (SD 9.83) years in the 100 mg arm",
    weight_range   = "male >= 50.0 kg, female >= 45.0 kg (protocol inclusion criterion)",
    weight_median  = "mean 69.55 (SD 7.63) kg in the 50 mg arm and 72.56 (SD 7.62) kg in the 100 mg arm",
    sex_female_pct = 20,
    race_ethnicity = "Chinese (single-centre study at the First Affiliated Hospital with Nanjing Medical University)",
    disease_state  = "Type 2 diabetes mellitus (WHO 1999 criteria), newly diagnosed and untreated or controlled by diet and exercise with no hypoglycaemic drugs in the preceding 12 weeks. Inclusion required 6.5% <= HbA1c < 9%, fasting blood glucose < 13.4 mmol/L, and BMI 19.00-30.00 kg/m2. Baseline HbA1c was 8.21 (SD 0.66)% in the 50 mg arm and 7.79 (SD 0.53)% in the 100 mg arm; baseline fasting plasma glucose was 7.96 (SD 1.29) and 6.87 (SD 1.40) mmol/L (Table 1).",
    dose_range     = "50 mg or 100 mg cetagliptin orally once daily under fasting conditions for 14 consecutive days",
    regions        = "China (Nanjing; single centre, CTR20190599)",
    hepatic_function = "No subject had clinically significant abnormal liver function (Results section 3.2).",
    notes          = paste(
      "Randomised, double-blind, placebo- and positive-controlled single- and",
      "multiple-dose study. 32 Chinese adults with T2DM were enrolled in two dose",
      "groups of 16, randomised 10:2:4 within each group to cetagliptin (50 or",
      "100 mg), matching placebo, or open-label sitagliptin 100 mg; 20 subjects",
      "therefore received cetagliptin. n_subjects is recorded as those 20",
      "cetagliptin recipients because only they contribute cetagliptin",
      "concentrations. NOTE: Methods section 2.7.1 states that the 560 cetagliptin",
      "concentrations used for the population PK analysis came from '32 patients",
      "with T2DM', which cannot be reconciled with the 20 cetagliptin recipients",
      "reported in Table 1 and Table 2; the discrepancy is documented in the",
      "vignette Errata. The population PK/PD model was fit to 554 concentrations",
      "after excluding 6 records with |CWRES| > 5 (Methods section 2.7.2).",
      "Modelling was performed in Phoenix NLME by a sequential two-step approach:",
      "the population PK model was fit first and its parameters were then fixed",
      "when fitting the PK/PD model (Methods section 2.7)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Population PK -- Supplementary Table 1 ("Parameter estimates and
    # bootstrap results of the final population pharmacokinetic model"),
    # block "Population pharmacokinetic parameters".
    #
    # The supplement reports amounts in ug and concentrations in ug/L
    # (= ng/mL). Every value below is therefore the verbatim
    # Supplementary Table 1 estimate; the mg -> ug dose conversion is
    # applied once in model() via f(depot).
    #
    # All parameters are apparent (oral, /F): the paper never estimated
    # bioavailability and administered cetagliptin only by mouth.
    # ------------------------------------------------------------------
    lka   <- log(0.06521); label("First-order absorption rate constant ka (1/h)")                       # Supplementary Table 1: tvKa 0.06521 1/h (RSE 1.623%)
    lvc   <- log(8.668);   label("Apparent central volume of distribution V/F (L)")                     # Supplementary Table 1: tvV 8.668 L (RSE 2.234%)
    lvp   <- log(558.0);   label("Apparent peripheral volume of distribution V2/F (L) at TBILI 10 umol/L")  # Supplementary Table 1: tvV2 558.0 L (RSE 2.023%)
    lq    <- log(9.671);   label("Apparent intercompartmental clearance Cl2/F (L/h)")                   # Supplementary Table 1: tvCl2 9.671 L/h (RSE 2.050%)
    lkm   <- log(171.5);   label("Michaelis-Menten constant Km (ng/mL)")                                # Supplementary Table 1: tvKm 171.5 ug/L = 171.5 ng/mL (RSE 6.694%)
    lvmax <- log(9373);    label("Apparent maximum elimination rate Vmax/F (ug/h)")                     # Supplementary Table 1: tvVmax 9373 ug/h (RSE 5.044%)

    # Covariate effect. Power form on the peripheral volume, centred at
    # the assumed reference TBILI of 10 umol/L (see covariateData$TBILI
    # notes and the vignette Errata -- the paper reports neither the
    # equation nor the population median). Estimated, not fixed: the
    # supplement gives RSE 3.066% and a bootstrap 95% CI of 0.1417-0.5453.
    e_tbili_vp <- 0.3723;  label("Power exponent of total bilirubin on peripheral volume V2 (unitless)")  # Supplementary Table 1: dV2dTBIL 0.3723 (RSE 3.066%)

    # ------------------------------------------------------------------
    # Population PK/PD -- Supplementary Table 1, block "Population
    # pharmacokinetic/pharmacodynamic parameters". Direct-effect sigmoid
    # Emax model on plasma DPP-4 inhibition (Methods section 2.7.2).
    # ------------------------------------------------------------------
    lemax <- log(91.78);   label("Maximum DPP-4 inhibition Emax (%)")                                   # Supplementary Table 1: tvEmax 91.78 % (RSE 0.6484%)
    lec50 <- log(5.120);   label("Cetagliptin concentration giving half-maximal DPP-4 inhibition EC50 (ng/mL)")  # Supplementary Table 1: tvEC50 5.120 ug/L (RSE 3.850%)
    lhill <- log(1.008);   label("Hill coefficient on the DPP-4 inhibition concentration-response curve (unitless)")  # Supplementary Table 1: tvGam 1.008 (RSE 4.296%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- Supplementary Table 1, block
    # "Inter-individual variability". Methods Equation 1 specifies the
    # exponential (log-normal) form Pij = Pj * exp(eta_ij), so the
    # tabulated omega^2 values are already log-scale variances and are
    # carried over unchanged. Shrinkage from the same table.
    # Cl2 carries no IIV in the published model.
    # ------------------------------------------------------------------
    etalvc   ~ 0.8112     # Supplementary Table 1: omega2 V 0.8112 (RSE 8.340%, shrinkage 18.24%)
    etalvmax ~ 0.0260     # Supplementary Table 1: omega2 Vmax 0.0260 (RSE 4.482%, shrinkage 8.642%)
    etalka   ~ 0.01462    # Supplementary Table 1: omega2 Ka 0.01462 (RSE 2.673%, shrinkage 11.31%)
    etalvp   ~ 0.02874    # Supplementary Table 1: omega2 V2 0.02874 (RSE 1.115%, shrinkage 12.68%)
    etalkm   ~ 0.01278    # Supplementary Table 1: omega2 Km 0.01278 (RSE 1.622%, shrinkage 43.52%)
    etalec50 ~ 0.01259    # Supplementary Table 1: omega2 EC50 0.01259 (RSE 57.07%, shrinkage 18.13%)
    etalemax ~ 7.262e-05  # Supplementary Table 1: omega2 Emax 7.262E-05 (RSE 58.50%, shrinkage 17.58%)

    # ------------------------------------------------------------------
    # Residual variability -- Supplementary Table 1, final two blocks.
    # ------------------------------------------------------------------
    # PK: "Multiplicative residual variability PK(sigma), Stdev0 0.2241".
    # Phoenix's Multiplicative observation model is CObs = C * (1 + CEps),
    # i.e. proportional error in nlmixr2's linear space, with SD 0.2241.
    propSd <- 0.2241; label("Proportional residual error on cetagliptin plasma concentration (fraction)")  # Supplementary Table 1: PK sigma Stdev0 0.2241 (RSE 5.429%)

    # PD: "MixRatio residual variability PD(sigma), Stdev0 11.62". Phoenix's
    # Mixed-Ratio observation model carries an additive standard deviation
    # (Stdev0) plus a separate MixRatio coefficient scaling a proportional
    # component. Only Stdev0 is reported; the MixRatio coefficient appears
    # nowhere in the paper or supplement, so only the additive component is
    # encoded here and no proportional term is invented. 11.62 is on the
    # DPP-4 inhibition percentage-point scale -- it cannot be a proportional
    # fraction, which would imply a 1162% CV. See vignette Errata.
    addSd_dpp4Inhibition <- 11.62; label("Additive residual error on DPP-4 inhibition (percentage points)")  # Supplementary Table 1: PD sigma Stdev0 11.62 (RSE 22.63%)
  })

  model({
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    km   <- exp(lkm + etalkm)
    vmax <- exp(lvmax + etalvmax)

    # Total bilirubin raises the peripheral volume (Discussion page 10:
    # "V2 increases with the increase of TBIL").
    vp <- exp(lvp + etalvp) * (TBILI / 10)^e_tbili_vp

    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    # States hold ug and volumes are in L, so state/volume is ug/L = ng/mL.
    Cc <- central / vc
    Cp <- peripheral1 / vp

    # Two-compartment disposition with first-order oral absorption and
    # saturable Michaelis-Menten elimination from the central compartment
    # (Results section 3.6.1 selects the two-compartment structural model;
    # Supplementary Table 1 parameterises elimination as Vmax/Km rather
    # than as a clearance).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - vmax * Cc / (km + Cc) - q * (Cc - Cp)
    d/dt(peripheral1) <-  q * (Cc - Cp)

    # Unit conversion only, NOT a bioavailability estimate: doses are
    # entered in mg while the Supplementary Table 1 parameterisation works
    # in ug (Vmax in ug/h, Km in ug/L). All parameters remain apparent.
    f(depot) <- 1000

    # Direct-effect sigmoid Emax model on plasma DPP-4 inhibition (%),
    # Methods section 2.7.2. The equation is typeset in the PDF as
    # Emax*C^g/(EC50 + C^g), which drops the exponent on EC50 and is
    # dimensionally inconsistent with the paper's own definition of EC50 as
    # "plasma concentration of cetagliptin that achieves 50% of the maximum
    # drug effect". The canonical Phoenix Sigmoid-Emax form named in
    # Results section 3.6.2 is used instead; with the estimated Hill
    # coefficient of 1.008 the two forms differ by under 0.1%.
    dpp4Inhibition <- emax * Cc^hill / (ec50^hill + Cc^hill)

    Cc ~ prop(propSd)
    dpp4Inhibition ~ add(addSd_dpp4Inhibition)
  })
}
