Dave_2019_venetoclax <- function() {
  description <- "Integrated population PK / semimechanistic PD model for venetoclax and circulating B-lymphocytes in healthy female subjects. PK: two-compartment model with first-order absorption, an absorption lag time under fed conditions only, and first-order elimination; relative bioavailability depends on food (fasting / low-fat reference / moderate-or-high-fat), azithromycin coadministration, Chinese ethnicity, rifampin coadministration and a fixed power function of dose, and rifampin also raises apparent clearance. PD: Friberg-type lymphocyte model (proliferating pool, three transit compartments, circulating B-lymphocytes) with a (baseline / circulating)^gamma feedback on proliferation and a linear venetoclax effect increasing the first-order loss of circulating B-lymphocytes. The individual observed pre-dose B-lymphocyte count is the PD baseline. Final IIV and residual-error magnitudes are not reported by the source, so every random effect and residual SD is fixed(0): the model reproduces typical-value predictions only."
  reference <- paste(
    "Dave N, Gopalakrishnan S, Mensing S, Salem AH.",
    "Model-Informed Dosing of Venetoclax in Healthy Subjects:",
    "An Exposure-Response Analysis.",
    "Clin Transl Sci. 2019;12(6):625-632.",
    "doi:10.1111/cts.12665.",
    "Parameter estimates from Table 2; model code from the Supporting",
    "Information 'PK Model code' (CTS-12-625-s003.pdf) and",
    "'PD Model code' (CTS-12-625-s002.pdf) NONMEM control streams.",
    "Lymphocyte model backbone: Friberg LE et al. J Clin Oncol",
    "2002;20(24):4713-4721, doi:10.1200/JCO.2002.02.140.",
    sep = " "
  )
  vignette <- "Dave_2019_venetoclax"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L (venetoclax Cc); circulating B-lymphocytes (circ) in cells/uL"
  )

  covariateData <- list(
    FED = list(
      description = "Fed-vs-fasted indicator for the dose record: 1 = dosed within 30 minutes after a low-, moderate- or high-fat breakfast, 0 = dosed after a 10-hour overnight fast.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 with FED_LOWFAT = 1 (low-fat meal) is the model reference for bioavailability",
      notes = "Source NONMEM column FOOD (0 = fasting, 1 = low fat, >= 2 = moderate fat / high fat / any meal). FED = 1 when FOOD >= 1. Controls (i) the absorption lag time, which is 0 when fasted and exp(ltlag) when fed ('IF(FOOD.GE.1) ALAG1 = THETA(8)'), and (ii) together with FED_LOWFAT, the food effect on relative bioavailability.",
      source_name = "FOOD"
    ),
    FED_LOWFAT = list(
      description = "Low-fat-meal indicator for the dose record: 1 = dosed after a low-fat breakfast, 0 = fasting or moderate/high-fat breakfast.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (low-fat meal; the bioavailability reference category of the source model)",
      notes = "Source NONMEM column FOOD: FED_LOWFAT = 1 when FOOD == 1. With FED this reproduces the three food strata of the source: fasting (FED = 0) -> F x 0.29; low fat (FED = 1, FED_LOWFAT = 1) -> F x 1; moderate/high fat (FED = 1, FED_LOWFAT = 0, source FOOD >= 2) -> F x 1.47. Meal definitions are per the pooled studies' protocols and not restated in Dave 2019.",
      source_name = "FOOD"
    ),
    CONMED_AZITHROMYCIN = list(
      description = "Concomitant azithromycin (P-gp inhibitor) coadministration indicator for the venetoclax dose record: 1 = coadministered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no azithromycin)",
      notes = "Source NONMEM column AZIFL. Study V (Agarwal 2018, Adv Ther 35:2015) DDI arm. Multiplies relative bioavailability by 0.65.",
      source_name = "AZIFL"
    ),
    RACE_CHINESE = list(
      description = "Chinese-ethnicity indicator: 1 = subject from the study in healthy Chinese subjects, 0 = other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Chinese; the rest of the pooled healthy female cohort)",
      notes = "The source codes this by study, 'IF(STDY.EQ.14752) TVFA = THETA(12) ; Chinese Ethnicity' (Study I, Cheung 2018, Clin Pharmacol Drug Dev 7:435; 12 subjects). Multiplies relative bioavailability by 1.53.",
      source_name = "STDY == 14752"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampin (rifampicin) coadministration indicator for the dose / observation record: 1 = rifampin on board, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampin)",
      notes = "The source codes this by study and time, 'IF(STDY.EQ.14497.AND.TIME.GE.8)' (Study X, Agarwal 2016, J Clin Pharmacol 56:1335: venetoclax alone on day 1, then with single-dose and multiple-dose rifampin from day 8). The source pools the single-dose (transporter-inhibition) and multiple-dose (CYP3A-induction) phases into ONE indicator with two effects: apparent clearance x 2.67 and relative bioavailability x 4.91. Time-varying: set to 1 on dose AND observation records while rifampin is on board so the clearance effect applies throughout.",
      source_name = "STDY == 14497 & TIME >= 8"
    ),
    DOSE = list(
      description = "Administered venetoclax dose level (mg) of the current dose record.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a) of the DOSE register entry: drives the fixed dose nonlinearity in relative bioavailability, F1 x (DOSE / 400)^-0.178 (Table 2 'Dose nonlinearity on F1' = -0.178, Fixed; control stream '(DOSE/400)**(THETA(14))'). Reference dose 400 mg. Pooled doses 10-400 mg.",
      source_name = "DOSE"
    ),
    BLBCELL = list(
      description = "Individual baseline circulating CD19+ B-lymphocyte count (cells/uL): the median of the subject's observed B-lymphocyte counts before venetoclax administration.",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Source NONMEM column MIBTCD19 ('Lymbase = MIBTCD19'; 'Baseline B-Lymphocytes (10^6 cells per L)', numerically equal to cells/uL). Methods: 'The median of observed B-lymphocyte counts prior to venetoclax administration was used as the baseline, with an associated interindividual variability.' Required input for the PD layer; it sets the initial condition of every PD state. The Discussion recommends enrolling only subjects with a baseline >= 150 cells/uL.",
      source_name = "MIBTCD19"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(
      analyte = "proliferating B-lymphocyte precursors (e.g. bone marrow)",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "maturing B-lymphocytes (transit stage 1)",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing B-lymphocytes (transit stage 2)",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing B-lymphocytes (transit stage 3)",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(
      analyte = "circulating CD19+ B-lymphocytes",
      units = "cells/uL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 203L,
    n_studies = 10L,
    age_range = "21-65 years (median 40)",
    weight_range = "not reported; body mass index 18.6-30.1 kg/m^2 (inclusion criterion)",
    sex_female_pct = 100,
    race_ethnicity = "12 of 203 subjects were Chinese (Study I); the remaining race distribution is not reported",
    disease_state = "healthy volunteers",
    dose_range = "venetoclax 10-400 mg single oral doses, given 1-4 times per subject separated by >= 7-day washouts, fasting or after a low-, moderate- or high-fat breakfast",
    regions = "not reported",
    co_medication = "ritonavir (Study III), digoxin (IV), azithromycin (V), rifampin (X) and warfarin (II) DDI studies were pooled",
    biomarkers = "CD19+ B-lymphocyte counts at baseline and after dosing in 7 of the 10 studies (Studies III, VIII and IX had no B-lymphocyte data)",
    notes = "Pooled from 10 healthy-volunteer clinical pharmacology studies (Table 1). Only female subjects were enrolled because 4-week dog toxicology suggested venetoclax may compromise male fertility. The PK model was fit to all 203 subjects; the PD model to the 7 studies with B-lymphocyte data, using individual post hoc PK parameters (the number of PD subjects is not reported). In the PD control stream the IGNORE list removes studies 14253, 14497 and 15101 and two IDs."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters. Dave 2019 Table 2 ('PK model'); control
    # stream CTS-12-625-s003.pdf ($SUBROUTINE ADVAN4 TRANS4). Time unit is
    # day (CL/F in L/day, KA in 1/day, ALAG in day); dose in mg so Cc is in
    # mg/L. The $THETA block of the supplement lists INITIAL values; every
    # value below is the FINAL estimate from Table 2.
    # ------------------------------------------------------------------
    lcl <- log(449); label("Apparent clearance CL/F (L/day)") # Table 2 'CL/F (L/day)' = 449 (95% CI 392-506)
    lvc <- log(99); label("Apparent central volume of distribution V2/F (L)") # Table 2 'V2 (L)' = 99 (95% CI 82-116)
    lq <- log(130); label("Apparent intercompartmental clearance Q/F (L/day)") # Table 2 'Q (L/day)' = 130 (95% CI 116-144)
    lvp <- log(147); label("Apparent peripheral volume of distribution V3/F (L)") # Table 2 'V3 (L)' = 147 (95% CI 135-159)
    lka <- log(3.83); label("First-order absorption rate constant KA (1/day)") # Table 2 'KA (1/day)' = 3.83 (95% CI 3.48-4.18)
    ltlag <- log(0.04); label("Absorption lag time under fed conditions ALAG (day)") # Table 2 'ALAG (day)' = 0.04 (95% CI 0.0399-0.0401); applied only when FED = 1

    # Rifampin effect on apparent clearance (multiplicative factor).
    e_conmed_rifampicin_cl <- 2.67; label("Multiplicative factor on CL/F with rifampin coadministration (unitless)") # Table 2 'Rifampin on CL' = 2.67 (95% CI 1.59-3.75); control stream CL = THETA(1)*THETA(15)

    # Relative bioavailability. Reference = low-fat meal, 400 mg, no
    # interacting co-medication, non-Chinese (F1 = 1, Table 2 'Low fat on
    # F1' = 1). The food factors and the co-medication / ethnicity factors
    # are separate multiplicative blocks (TVF1 * TVFA in the control stream).
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 of the reference condition (unitless)") # Table 2 'Low fat on F1' = 1; control stream 'TVF1=1' and 'TVFA=1'
    e_fasted_fdepot <- 0.29; label("Multiplicative factor on F under fasting vs low-fat meal (unitless)") # Table 2 'Fasting on F1' = 0.29 (95% CI 0.28-0.29)
    e_modhighfat_fdepot <- 1.47; label("Multiplicative factor on F with a moderate- or high-fat meal vs low-fat meal (unitless)") # Table 2 'Moderate/high-fat on F1' = 1.47 (95% CI 1.43-1.51)
    e_conmed_azithromycin_fdepot <- 0.65; label("Multiplicative factor on F with azithromycin coadministration (unitless)") # Table 2 'Azithromycin on F1' = 0.65 (95% CI 0.61-0.70)
    e_race_chinese_fdepot <- 1.53; label("Multiplicative factor on F for Chinese subjects (unitless)") # Table 2 'Chinese on F1' = 1.53 (95% CI 0.81-2.47)
    e_conmed_rifampicin_fdepot <- 4.91; label("Multiplicative factor on F with rifampin coadministration (unitless)") # Table 2 'Rifampin on F1' = 4.91 (95% CI 2.38-7.44)
    e_dose_fdepot <- fixed(-0.178); label("Power exponent of (DOSE / 400 mg) on F (unitless)") # Table 2 'Dose nonlinearity on F1' = -0.178, Fixed; control stream '-0.178 FIX'

    # ------------------------------------------------------------------
    # PD structural parameters. Dave 2019 Table 2 ('PD model'); control
    # stream CTS-12-625-s002.pdf ($SUBROUTINE ADVAN13). t1/2 and gamma were
    # estimated with informative NWPRI priors ($THETAP 38.5 and 0.1,
    # $THETAPV 3 and 0.00005); they are ESTIMATED, not fixed. KIN and the
    # drug-effect slope are estimated on the log scale (KIN = EXP(THETA(3)),
    # SLB = EXP(THETA(4))); Table 2 reports them back-transformed.
    # ------------------------------------------------------------------
    lthalf <- log(37.5); label("Half-life of circulating B-lymphocytes (day)") # Table 2 't 1/2 (days)' = 37.5 (95% CI 34.0-41.0); control stream DELTA = LOG(2)/THETA(1)
    lktr <- log(0.1); label("Proliferation and transit maturation rate constant kin = ktr (1/day)") # Table 2 'k in, proliferation (and transit maturation) rate (day-1)' = 0.1 (95% CI 0.02-0.61); control stream KTR = KIN
    lgamma <- log(0.1); label("Feedback exponent gamma on (baseline / circulating B-lymphocytes) (unitless)") # Table 2 'Feedback exponent' = 0.1 (95% CI 0.09-0.11)
    lslope <- log(20.9); label("Slope of the linear venetoclax effect on circulating B-lymphocyte loss (L/mg)") # Table 2 'Slope of drug effect' = 20.9 (95% CI 18.4-23.7); control stream EFFB = 1 + SLB*CP with CP in mg/L

    # ------------------------------------------------------------------
    # Inter-individual variability. The control streams declare exponential
    # IIV on CL/F and V2/F (a correlated $OMEGA BLOCK(2)), on F1, and on the
    # PD baseline; IIV on KA is '0 FIX', and the PD code's second ETA (on
    # both the half-life and KIN) is '0 FIX'. Dave 2019 reports NO final
    # omega estimates (Table 2 lists fixed effects only; the supplement
    # lists initial values only), so the declared etas are fixed(0) and the
    # zero-fixed ones are omitted. The CL-V2 covariance is likewise
    # unreported and not encoded. See the vignette Assumptions and
    # deviations.
    # ------------------------------------------------------------------
    etalcl ~ fixed(0) # control stream '$OMEGA BLOCK(2)' ETA(1) on CL; final value not reported
    etalvc ~ fixed(0) # control stream '$OMEGA BLOCK(2)' ETA(2) on V2; final value not reported
    etalfdepot ~ fixed(0) # control stream '$OMEGA' ETA(4) on F1; final value not reported
    etalcirc0 ~ fixed(0) # PD control stream '$OMEGA' ETA(1) on baseline B-lymphocytes; final value not reported

    # ------------------------------------------------------------------
    # Residual error. PK: combined proportional + additive on the linear
    # scale ('Y = F*(1+EPS(1))+EPS(2)'); PD: additive on log-transformed
    # B-lymphocyte counts ('IPRED = LOG(A(8)); Y = IPRED + EPS(1)'), i.e.
    # log-normal. Final sigma estimates are not reported (the $SIGMA block
    # holds initial values 0.2, 3E-7 and 0.02), so all are fixed(0).
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD for venetoclax Cc (fraction); magnitude not reported") # PK control stream 'Y = F*(1+EPS(1))+EPS(2)'; final value not reported
    addSd <- fixed(0); label("Additive residual SD for venetoclax Cc (mg/L); magnitude not reported") # PK control stream 'Y = F*(1+EPS(1))+EPS(2)'; final value not reported
    expSd_circ <- fixed(0); label("Log-scale additive residual SD for circulating B-lymphocytes (log units); magnitude not reported") # PD control stream 'Y = IPRED + EPS(1)' on log counts; final value not reported
  })

  model({
    # ---- PK individual parameters (control stream $PK) ----
    cl <- exp(lcl + etalcl) * e_conmed_rifampicin_cl^CONMED_RIFAMPICIN
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    # Food effects on KA are 'THETA(6) = 1 FIX' and 'THETA(7) = 1 FIX' and
    # IIV on KA is '0 FIX' in the control stream, so KA is a single value.
    ka <- exp(lka)
    # 'ALAG1 = 0; IF(FOOD.GE.1) ALAG1 = THETA(8)': lag under fed conditions only.
    tlag <- exp(ltlag) * FED

    # Food factor on F (TVF1): fasting 0.29, low fat 1 (reference),
    # moderate/high fat 1.47.
    f_food <- e_fasted_fdepot * (1 - FED) + FED * (FED_LOWFAT + (1 - FED_LOWFAT) * e_modhighfat_fdepot)
    # Co-medication / ethnicity factor on F (TVFA). The control stream
    # OVERWRITES rather than multiplies, in the order azithromycin ->
    # Chinese -> rifampin; the three groups come from different studies
    # and never overlap in the data, but the override order is kept.
    f_other <- 1
    if (CONMED_AZITHROMYCIN == 1) f_other <- e_conmed_azithromycin_fdepot
    if (RACE_CHINESE == 1) f_other <- e_race_chinese_fdepot
    if (CONMED_RIFAMPICIN == 1) f_other <- e_conmed_rifampicin_fdepot
    fdepot <- exp(lfdepot + etalfdepot) * f_food * f_other * (DOSE / 400)^e_dose_fdepot

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PD individual parameters (PD control stream $PK) ----
    kout <- log(2) / exp(lthalf) # DELTA, first-order loss of circulating cells
    ktr <- exp(lktr) # KIN = KTR
    gamma <- exp(lgamma)
    slope <- exp(lslope)
    # BASEB = EXP(LOG(Lymbase) + ETA(1)): the baseline is the subject's
    # observed pre-dose count, not an estimated theta, so etalcirc0 has no
    # lcirc0 partner (a deliberate, documented convention warning).
    circ0 <- exp(log(BLBCELL) + etalcirc0)

    # ---- ODEs ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    # Linear drug effect on the loss of circulating B-lymphocytes (EFFB).
    edrug <- 1 + slope * Cc

    # The 1E-6 inside the feedback term is kept from the control stream
    # ('(BASEB/A(8) + 1E-6)**GAM'); it guards against division by zero and
    # perturbs the baseline steady state by < 1e-7 relative.
    d/dt(prol) <- ktr * prol * (circ0 / circ + 1e-6)^gamma - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - kout * edrug * circ

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Baseline steady state ('A_0(4..7) = DELTA*BASEB/KTR; A_0(8) = BASEB').
    prol(0) <- kout * circ0 / ktr
    transit1(0) <- kout * circ0 / ktr
    transit2(0) <- kout * circ0 / ktr
    transit3(0) <- kout * circ0 / ktr
    circ(0) <- circ0

    Cc ~ add(addSd) + prop(propSd)
    circ ~ lnorm(expSd_circ)
  })
}
