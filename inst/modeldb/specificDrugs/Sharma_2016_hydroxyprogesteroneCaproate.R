Sharma_2016_hydroxyprogesteroneCaproate <- function() {
  description <- "Population PK model for 17alpha-hydroxyprogesterone caproate (17-OHPC) in pregnant women with singleton gestation receiving 250 mg weekly IM injections for prevention of recurrent preterm birth (Sharma 2016). Structural model is a maternal central compartment with first-order IM absorption (ka fixed to 3 /day) and first-order elimination, linked by reversible first-order rate constants kMF / kFM to a fetal compartment whose amount is tracked as a dynamic state. Allometric (power) scaling on apparent CL/F (exponent 0.80) and apparent Vmaternal/F (exponent 0.84) around the cohort median weight of 68 kg. Inter-individual variability is encoded as a shared random effect on log(CL) with an estimated scale factor (Theta6 = 1.90) applied to the same eta on log(Vmaternal), preserving the > 0.9 ETA correlation reported in the source. Inter-occasion variability (PK1 = 20-24 weeks gestation, PK2 = 31-34 weeks gestation) on CL and Vmaternal is encoded as paired per-occasion etas keyed to OCC. Residual error is combined proportional + additive on linear-scale ng/mL concentrations."
  reference <- paste(
    "Sharma S, Caritis S, Hankins G, Miodovnik M, Hebert MF, Mattison D,",
    "Venkataramanan R.",
    "Population pharmacokinetics of 17alpha-hydroxyprogesterone caproate",
    "in singleton gestation.",
    "Br J Clin Pharmacol. 2016 Nov;82(5):1084-1093.",
    "doi:10.1111/bcp.12990.",
    "ClinicalTrials.gov NCT00409825.",
    sep = " "
  )
  vignette <- "Sharma_2016_hydroxyprogesteroneCaproate"
  units    <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("fetal")

  covariateData <- list(
    WT = list(
      description        = "Baseline maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; baseline value recorded before the first 17-OHPC dose per Sharma 2016 Methods 'Covariate data'. Power-function covariate effect normalised to the cohort median of 68 kg per Sharma 2016 Results 'Data' (median weight at baseline = 68 kg, range 47-141 kg). Both CL/F (exponent 0.80, Theta_CL-WT in Table 2) and Vmaternal/F (exponent 0.84, Theta_Vmaternal-WT in Table 2) carry a (WT/68)^exp scaling per equations (9) and (10).",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued PK-sampling-occasion indicator (1 = PK1 visit at 20-24 weeks gestation, 2 = PK2 visit at 31-34 weeks gestation)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Time-varying within subject; constant within a PK sampling window. Sharma 2016 design (Methods 'Pharmacokinetic sampling schedule and sample analysis') uses two 7-day intensive PK occasions per subject: PK1 at 20-24 weeks gestation and PK2 at 31-34 weeks gestation. IOV on CL (CV 26.5%) and Vmaternal (CV 31.6%) is reported per Table 2 ('IOV - CL (interoccasion)' and 'IOV - Vmaternal/F (interoccasion)'). Decomposed inside model() into binary indicators oc1 / oc2 multiplexing per-occasion etas (etaiov_cl_1 / etaiov_cl_2 on CL, etaiov_vc_1 / etaiov_vc_2 on Vmaternal). The occasion-2 etas are fix()'d to the occasion-1 variance to preserve the NONMEM $OMEGA BLOCK(1) SAME idiom (single IOV variance per parameter; nlmixr2 has no SAME shortcut).",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the forward-inclusion screen on CL (Delta-OFV = -11, P < 0.001) and on Vmaternal (Delta-OFV = -2, P > 0.2) but not retained in the final model after the body-weight effect was included (Sharma 2016 Table 1).",
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Maternal age at enrolment",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the forward-inclusion screen on CL and Vmaternal (both Delta-OFV ~ 0, P > 0.2); not retained in the final model (Sharma 2016 Table 1). Median age 27 years (range 19-42 years) per Results 'Data'.",
      source_name        = "AGE"
    ),
    GA = list(
      description        = "Gestational age at first 17-OHPC injection",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the forward-inclusion screen on CL and Vmaternal (both Delta-OFV ~ 0, P > 0.2); not retained in the final model (Sharma 2016 Table 1). Median GA at first injection 17 weeks (range 16-21 weeks) per Results 'Data'. Recorded once per subject at enrolment.",
      source_name        = "GA"
    ),
    RACE_BLACK = list(
      description        = "Race indicator -- Black / African-American (1 = yes, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Tested as one of four race-and-ethnicity categories (White / Black / Hispanic / Others) in the forward-inclusion screen on CL and Vmaternal (Delta-OFV ~ 0, P > 0.2); not retained in the final model despite a ~32% higher median CL/F in African-Americans (Sharma 2016 Discussion 'Effect of race and ethnicity'). Cohort composition per Results 'Data': White 51%, Black 19%, Hispanic 27%, Other 3% (n = 30, 11, 16, 2 respectively in Figure 4 caption).",
      source_name        = "RACE_BLACK"
    ),
    RACE_HISPANIC = list(
      description        = "Race / ethnicity indicator -- Hispanic (1 = yes, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Tested with the Black indicator as part of the four-category race-and-ethnicity covariate; not retained (Delta-OFV ~ 0, P > 0.2; Sharma 2016 Table 1). ~10% higher CL/F vs Caucasians per Discussion 'Effect of race and ethnicity'.",
      source_name        = "RACE_HISPANIC"
    ),
    PROGESTERONE = list(
      description        = "Baseline plasma progesterone concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a continuous covariate on CL (Delta-OFV ~ 0, P > 0.2); not retained (Sharma 2016 Table 1). Median baseline progesterone 52.4 ng/mL per Discussion 'Effect of hormones'. The weak inverse correlation seen in Figure 4 was not statistically significant.",
      source_name        = "PROGESTERONE"
    ),
    HYDROXYPROGESTERONE = list(
      description        = "Baseline plasma 17-hydroxyprogesterone (endogenous) concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a continuous covariate on CL (Delta-OFV ~ 0, P > 0.2); not retained (Sharma 2016 Table 1). Median baseline hydroxyprogesterone 2.7 ng/mL per Discussion 'Effect of hormones'. Endogenous 17-OHP is structurally similar to 17-OHPC and shares CYP3A metabolism, motivating its inclusion in the screen.",
      source_name        = "HYDROXYPROGESTERONE"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 59L,
    n_studies       = 1L,
    n_observations  = 1514L,
    age_range       = "19-42 years",
    age_median      = "27 years",
    weight_range    = "47-141 kg",
    weight_median   = "68 kg",
    sex_female_pct  = 100,
    race_ethnicity  = c(White = 51, Black = 19, Hispanic = 27, Other = 3),
    disease_state   = "Pregnant women with singleton gestation and a history of prior preterm birth, receiving prophylactic 17-OHPC for prevention of recurrent preterm birth",
    ga_range        = "16-21 weeks gestation at first injection; treatment continued to 35 weeks gestation or delivery",
    dose_range      = "250 mg 17-OHPC in 1 mL castor oil, intramuscular injection, once weekly",
    regions         = "United States (multi-centre clinical study, NICHD-funded Obstetrical-Fetal Pharmacology Research Network)",
    notes           = "Demographics from Sharma 2016 Results 'Data'. 61 women enrolled (Methods 'Patients and drug administration'); 2 excluded for missing dosing information; 59 contributed to the analysis. 1514 plasma 17-OHPC concentrations including 18 cord-blood concentrations contributed to the population PK fit. Sampling design: two 7-day intensive PK occasions (PK1 at 20-24 weeks, PK2 at 31-34 weeks) per subject, with a 28-day extended-sampling cohort of 18 women drawn on days 9, 11, 14, 17, 20, 24, 28 after the start of PK2; cord blood collected when possible at delivery from a maternal vein and umbilical cord (artery or vein). Bioanalytical assay: HPLC-MS/MS, range 1-200 ng/mL, LLOQ 1 ng/mL, inter-/intra-assay variability 7.9%/5.2% at 10 ng/mL. Estimation method: FOCEI in NONMEM VII (ICON Development Solutions); bootstrap n = 500 (490 converged) per Results 'Final model evaluation'. ClinicalTrials.gov NCT00409825."
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Sharma 2016 Table 2 ('Population
    # estimate (%RSE)' column). Apparent values (absorbing F) for a
    # reference subject at the cohort median weight of 68 kg per
    # equations (9) and (10) (Sharma 2016, page 1087).
    # ============================================================
    lka <- fixed(log(3.0))
    label("Absorption rate constant ka (1/day; FIXED)")
    # Sharma 2016 Table 2: ka = 3.0 (Fixed); Methods 'Base model' (page 1086):
    # 'Due to lack of adequate sampling time points during the absorption phase
    # (Day 0 and 1), complex absorption profiles could not be characterized
    # and a simpler first-order absorption model (ka fixed to 3) was chosen.
    # A value of ka = 3 was estimated based on the model but subsequently
    # fixed due to low precision of the estimate.'

    lcl <- log(1797)
    label("Apparent maternal clearance CL/F at WT = 68 kg (L/day)")
    # Sharma 2016 Table 2: CL/F = 1797 L/day (%RSE 4)

    lvc <- log(32610)
    label("Apparent maternal volume of distribution Vmaternal/F at WT = 68 kg (L)")
    # Sharma 2016 Table 2: Vmaternal/F = 32 610 L (%RSE 7)

    lkmf <- log(0.005)
    label("Maternal-to-fetal first-order transfer rate constant kMF (1/day)")
    # Sharma 2016 Table 2: kmaternal to fetal = 0.005 1/day (%RSE 39);
    # equations (6) and (7) (page 1087) define it as an amount-on-amount
    # rate constant operating on Amaternal.

    lkfm <- log(0.56)
    label("Fetal-to-maternal first-order transfer rate constant kFM (1/day)")
    # Sharma 2016 Table 2: kfetal to maternal = 0.56 1/day (%RSE 40);
    # equations (6) and (7) (page 1087) define it as an amount-on-amount
    # rate constant operating on Afetal.

    # ------------------------------------------------------------
    # Covariate-effect parameters: power exponents on baseline maternal
    # body weight, normalised to the cohort median 68 kg. Equations (9)
    # and (10) (Sharma 2016 page 1087):
    #   CL/F        = Theta1 * (WT / Median(WT))^Theta7 * exp(eta1)
    #   Vmaternal/F = Theta2 * (WT / Median(WT))^Theta8 * exp(Theta6 * eta1)
    # with Median(WT) = 68 kg per Results 'Data'.
    # ------------------------------------------------------------
    e_wt_cl <- 0.80
    label("Power exponent of (WT/68) on CL/F (unitless)")
    # Sharma 2016 Table 2: Theta_CL-WT = 0.80 (%RSE 21)

    e_wt_vc <- 0.84
    label("Power exponent of (WT/68) on Vmaternal/F (unitless)")
    # Sharma 2016 Table 2: Theta_Vmaternal-WT = 0.84 (%RSE 29)

    scale_etalvc <- 1.90
    label("Scale factor on the shared IIV random effect for Vmaternal (unitless; Theta6)")
    # Sharma 2016 Table 2: Theta6 (Shared ETA Scale Factor) = 1.90 (%RSE 24).
    # Methods 'Base model' (page 1086): 'Individual ETAs for CL/F and
    # Vmaternal/F were highly correlated (r > 0.9), so a shared ETA was
    # used with an estimated scale factor applied to volume.' Encoded
    # following the Leshinsky_2017_caspofungin_cat 'scale_etalvc' precedent.

    # ============================================================
    # Inter-individual variability -- Sharma 2016 Table 2. The shared-eta
    # encoding (Methods 'Base model', page 1086) means a single random
    # effect etalcl drives BOTH CL and Vmaternal: log(CL_i) = log(CL_pop)
    # + etalcl_i, log(Vmaternal_i) = log(Vmaternal_pop) + scale_etalvc *
    # etalcl_i. Table 2 reports both %IIV - CL and %IIV - Vmaternal/F as
    # 20.0% CV (the variance of the underlying shared eta); the effective
    # IIV on Vmaternal is amplified by scale_etalvc = 1.90 to give
    # ~38% CV (Discussion paragraph 4: 'The volume for the maternal
    # compartment was estimated to be 32 610 L (IIV ~38% and IOV ~32%).').
    # Internal variance: omega^2 = log(1 + CV^2) = log(1 + 0.20^2) = 0.0392.
    # ============================================================
    etalcl ~ log(1 + 0.20^2)
    # Sharma 2016 Table 2: %IIV - CL = 20.0% CV (%RSE 17); same shared eta
    # also drives Vmaternal via scale_etalvc.

    # ============================================================
    # Inter-occasion variability -- Sharma 2016 Table 2 (%IOV - CL = 26.5%
    # CV, %IOV - Vmaternal/F = 31.6% CV). Two PK occasions per subject
    # (PK1 = 20-24 weeks, PK2 = 31-34 weeks; see Methods 'Pharmacokinetic
    # sampling schedule and sample analysis'). Encoded as paired
    # per-occasion etas with the occasion-2 variance fix()'d equal to
    # occasion-1 -- the NONMEM $OMEGA BLOCK(1) SAME idiom (single IOV
    # variance per parameter). IOV is independent for CL and Vmaternal
    # (the shared-eta structure applies only to IIV); evidence: IOV CV
    # values differ (26.5% vs 31.6%) and are not related by the
    # scale_etalvc factor.
    # ============================================================
    etaiov_cl_1 ~ log(1 + 0.265^2)
    # Sharma 2016 Table 2: %IOV - CL = 26.5% CV (%RSE 15) -> omega^2 = log(1 + 0.265^2)
    etaiov_cl_2 ~ fixed(log(1 + 0.265^2))
    # Occasion-2 variance fix()'d equal to occasion-1 (NONMEM $OMEGA BLOCK(1) SAME).

    etaiov_vc_1 ~ log(1 + 0.316^2)
    # Sharma 2016 Table 2: %IOV - Vmaternal/F = 31.6% CV (%RSE 14) -> omega^2 = log(1 + 0.316^2)
    etaiov_vc_2 ~ fixed(log(1 + 0.316^2))
    # Occasion-2 variance fix()'d equal to occasion-1 (NONMEM $OMEGA BLOCK(1) SAME).

    # ============================================================
    # Residual error -- Sharma 2016 Table 2. Combined proportional +
    # additive error model per Methods 'Base model development' (page
    # 1086) and Results 'Base model' (page 1087: 'The model with combined
    # proportional and additive residual errors was superior to the
    # proportional error or additive error model').
    # ============================================================
    propSd <- 0.173
    label("Proportional residual error (fraction)")
    # Sharma 2016 Table 2: % Residual error - proportional = 17.3% (%RSE 5);
    # NONMEM reports the proportional component as %CV, which is numerically
    # equal to the proportional-residual SD.

    addSd <- 1.0
    label("Additive residual error SD (ng/mL)")
    # Sharma 2016 Table 2: Residual error - additive = 1.0 ng/mL (%RSE 23).
  })

  model({
    # ------------------------------------------------------------
    # Occasion-indicator decomposition for IOV. OCC takes integer values
    # 1 (PK1, 20-24 weeks gestation) or 2 (PK2, 31-34 weeks gestation);
    # see covariateData$OCC$notes. Following the deWit_2016_everolimus
    # and Jonsson_2011_ethambutol IOV idiom (single shared variance,
    # per-occasion eta multiplexed by binary indicators).
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # ------------------------------------------------------------
    # Individual PK parameters. Sharma 2016 equations (9) and (10):
    #   CL/F        = Theta1 * (WT / 68)^e_wt_cl * exp(etalcl + iov_cl)
    #   Vmaternal/F = Theta2 * (WT / 68)^e_wt_vc * exp(scale_etalvc * etalcl + iov_vc)
    # The shared etalcl drives both CL and Vmaternal; scale_etalvc = 1.90
    # amplifies the effective IIV on Vmaternal. IOV is independent per
    # parameter (different reported CV; see comment on etaiov_*).
    # Operand parens guard against the rxode2 mu-ref false-positive
    # triggered by `theta + theta` adjacency.
    # ------------------------------------------------------------
    ka  <- exp(lka)
    cl  <- exp((lcl) + (etalcl) + (iov_cl)) * (WT / 68)^e_wt_cl
    vc  <- exp((lvc) + (scale_etalvc * etalcl) + (iov_vc)) * (WT / 68)^e_wt_vc
    kmf <- exp(lkmf)
    kfm <- exp(lkfm)
    kel <- cl / vc

    # ------------------------------------------------------------
    # Structural ODE system -- Sharma 2016 equations (5), (6), (7), and
    # (8) (page 1087). All states carry amount (dose units; mg). The
    # 'central' compartment is the maternal central compartment; the
    # paper-specific 'fetal' compartment is the fetal / cord-blood
    # circulation. Maternal-to-fetal exchange is reversible via
    # first-order rate constants kMF (slow) and kFM (fast).
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central - kmf * central + kfm * fetal
    d/dt(fetal)   <-  kmf * central - kfm * fetal

    # ------------------------------------------------------------
    # Observation. Dose in mg; vc in L -> central / vc gives maternal
    # plasma concentration in mg/L; multiplying by 1000 converts to
    # ng/mL (the paper's bioanalytical scale; LLOQ 1 ng/mL,
    # Methods 'Pharmacokinetic sampling schedule and sample analysis').
    # Combined proportional + additive residual error is applied to
    # the maternal-compartment Cc only.
    #
    # NOTE on the fetal compartment observable: Sharma 2016 equations
    # (5)-(8) define kMF and kFM as amount-on-amount rate constants and
    # do NOT report a fetal apparent volume of distribution. The paper
    # does fit the model against cord-blood concentrations (Figure 3
    # right VPC) and the Discussion (page 1089) reports a 'cord:maternal
    # (model predicted) ~ 0.28' ratio, but the supplemental NONMEM
    # control stream (Table S1) which encodes the cord-blood
    # observation prediction is not on disk. The fetal-compartment
    # amount 'fetal' is exposed as a dynamic state for downstream
    # interpretation, but no Ccord observable / residual error is
    # declared here. See the vignette Assumptions and deviations
    # section.
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
