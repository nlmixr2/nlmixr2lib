Koh_2025_aspirin <- function() {
  description <- "Joint population PK-PD model for enteric-coated aspirin in healthy Korean adults (Koh 2025; two phase I studies, 100 mg once daily for 5 days). Parallel dual absorption -- a zero-order arm with a lag time plus a first-order arm -- delivers acetylsalicylic acid (ASA) into a pre-systemic compartment that drains by two competing first-order routes: intact ASA into the one-compartment ASA disposition, and gut-hydrolysed salicylic acid (SA) straight into the two-compartment SA disposition. Systemic ASA is eliminated entirely by conversion to SA (kmet), so apparent ASA clearance is kmet * V3/F = 69.8 L/h. Body weight enters kmet and SA clearance as power functions normalised to the 68.35 kg cohort median, and the enteric-coated capsule absorbs about four times faster than the enteric-coated tablet (ka 0.22 vs 0.053 1/h). A turnover model with an Imax function on ASA concentration inhibits serum thromboxane B2 production."
  reference   <- "Koh JE, Khwarg J, Yu KS, Lee S, Jang IJ, Lee S. Population Pharmacokinetic and Pharmacodynamic Modeling of Enteric-Coated Aspirin Capsule and Tablet Formulations in Healthy Subjects. Drug Des Devel Ther. 2025;19:7853-7863. doi:10.2147/DDDT.S533428"
  vignette    <- "Koh_2025_aspirin"
  units       <- list(time = "h", dosing = "umol", concentration = "umol/L")
  # Both parallel absorption arms are dosed explicitly: the first-order arm
  # through `depot` (fraction fdepot) and the zero-order arm straight into
  # `presystemic` (fraction 1 - fdepot, over `dur` = d1, after `alag` = tlag).
  # Declared because the registry's auto-detection only recognises depot and
  # central and would otherwise record `central`, which is never dosed here.
  dosing      <- c("depot", "presystemic")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight. Enters kmet and cl_sa as power functions normalised to the 68.35 kg cohort median (Koh 2025 Table S1 'Total' column, median [min-max] 68.35 [55.6-89.5] kg); the median normalisation form is Koh 2025 Covariate Analysis equation. Enrolment restricted weight to 50.0-90.0 kg and BMI to 18.0-27.0 kg/m2, so the model is not informed outside that range.",
      source_name        = "weight"
    ),
    FORM_CAPSULE = list(
      description        = "Enteric-coated capsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = enteric-coated tablet (Aspirin Protect, Bayer, Germany)",
      notes              = "1 = enteric-coated capsule (Astrix, Boryungbio, Korea); 0 = enteric-coated tablet (Aspirin Protect, Bayer, Germany). Unlike the bioavailability-oriented uses of this canonical elsewhere in the register, here the formulation acts only on the first-order absorption rate constant: ka = 0.22 1/h for the capsule (estimated, RSE 21.8%) versus 0.053 1/h for the tablet (fixed during covariate analysis to stabilise the model, Koh 2025 Covariate Analysis). Relative bioavailability is not estimated -- every parameter is apparent (/F) and F is common to both arms. The tablet is the reference category because it carries the fixed value. Koh 2025 Results: formulation on ka was retained with a correlation test statistic of -3.13 (p < 0.005); Tk0, Lag0 and fr were tested for the formulation difference and rejected.",
      source_name        = "formulation"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Koh 2025 Covariate Analysis) but not retained in the final model. Table S1 total median [min-max] 29 [20-49] years."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained. Table S1 total median [min-max] 174.55 [157.7-185.5] cm."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m2",
      type        = "continuous",
      notes       = "Screened but not retained. Table S1 total median [min-max] 22.65 [19.5-26.9] kg/m2."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "acetylsalicylic acid", units = "umol",
      specimen = "administration site", verified = TRUE
    ),
    presystemic = list(
      analyte = "acetylsalicylic acid", units = "umol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "acetylsalicylic acid", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    central_sa = list(
      analyte = "salicylic acid", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_sa = list(
      analyte = "salicylic acid", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    txb2 = list(
      analyte = "thromboxane B2", units = "ug/L",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44,
    n_studies      = 2,
    age_range      = "20-49 years",
    age_median     = "29 years",
    weight_range   = "55.6-89.5 kg",
    weight_median  = "68.35 kg",
    sex_female_pct = 6.8,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy volunteers",
    dose_range     = "aspirin 100 mg orally once daily for 5 days (enteric-coated capsule or enteric-coated tablet)",
    regions        = "Republic of Korea (Seoul National University Hospital)",
    notes          = "Pooled period-1 data from two open-label, two-period, one-sequence crossover phase I aspirin-rabeprazole interaction studies (NCT05481307, NCT05699070). Only period 1 -- aspirin alone, without the rabeprazole perpetrator -- was used to build the model. 21 subjects received the enteric-coated capsule and 23 the enteric-coated tablet; 41 male / 3 female. 779 observations: 317 ASA plasma, 379 SA plasma, 83 TXB2 serum. PK sampled at 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 12 and 24 h on days 1 and 5; TXB2 only pre-dose on day 1 and 24 h after the day-5 dose. LLOQ 10.0 ng/mL (ASA) and 200 ng/mL (SA); BLQ handled by the M1 method. Estimated by SAEM in Monolix 2021R1. Demographics from Koh 2025 Table S1."
  )

  ini({
    # -----------------------------------------------------------------
    # Absorption. Koh 2025 Figure 1 / Appendix 1: the dose splits between
    # a zero-order arm (fraction fr, lag time lag0, duration Tk0) and a
    # first-order arm (fraction 1 - fr, rate constant ka). BOTH arms
    # deliver into the pre-systemic compartment.
    #
    # The paper reports fr = 0.69 as the ZERO-order absorbed fraction. The
    # library canonical `fdepot` is the fraction absorbed via the depot,
    # i.e. the FIRST-order arm, so the stored value is the complement
    # 1 - 0.69 = 0.31 and the zero-order arm takes 1 - fdepot. IIV on fr
    # is logit-normal (Koh 2025 Methods equation); because
    # logit(1 - p) = -logit(p), the variance is unchanged by the
    # complement and omega carries over untouched.
    logitfdepot <- log(0.31 / 0.69); label("Logit of the fraction of dose absorbed by the first-order arm (unitless)")  # Koh 2025 Table 1: fr = 0.69 (zero-order fraction), RSE 4.99% -> logit(1 - 0.69)

    lka_capsule <- log(0.22); label("First-order absorption rate constant, enteric-coated capsule (1/h)")               # Koh 2025 Table 1: ka capsule = 0.22 1/h, RSE 21.8%
    lka_tablet  <- fixed(log(0.053)); label("First-order absorption rate constant, enteric-coated tablet (1/h)")        # Koh 2025 Table 1: ka tablet = 0.053 1/h, no RSE (fixed during covariate analysis)
    ld1         <- log(1.58); label("Duration of the zero-order absorption arm Tk0 (h)")                                # Koh 2025 Table 1: Tk0 = 1.58 h, RSE 15.8%
    ltlag       <- log(2.81); label("Lag time before the zero-order absorption arm starts, Lag0 (h)")                   # Koh 2025 Table 1: Lag0 = 2.81 h, RSE 8.26%

    # -----------------------------------------------------------------
    # Pre-systemic compartment exits. Koh 2025 Appendix 1:
    #   ddt_PSYS = -k23 * PSYS - k24 * PSYS
    # k23 carries intact ASA onward to the ASA central compartment;
    # k24 is pre-systemic (gut) hydrolysis delivering SA directly into
    # the SA central compartment.
    lk_presystemic_central    <- log(2.32); label("Pre-systemic to ASA-central first-order rate constant k23 (1/h)")               # Koh 2025 Table 1: k23 = 2.32 1/h, RSE 4.11%
    lk_presystemic_central_sa <- fixed(log(0.57)); label("Pre-systemic to SA-central first-order rate constant k24 (1/h)")         # Koh 2025 Table 1: k24 = 0.57 1/h, no RSE (fixed from Dings & Lehr, Koh 2025 ref 14, to aid convergence)

    # -----------------------------------------------------------------
    # ASA disposition (one compartment). ASA is eliminated ONLY by
    # conversion to SA, so apparent ASA clearance is kmet * vc.
    lkmet     <- log(2.97); label("ASA-to-SA systemic metabolism rate constant k34 (1/h)")   # Koh 2025 Table 1: k34 = 2.97 1/h, RSE 11.7%
    lvc       <- log(23.51); label("Apparent ASA central volume V3/F (L)")                   # Koh 2025 Table 1: V3/F = 23.51 L, RSE 12.3%
    e_wt_kmet <- 1.31; label("Power exponent of WT/68.35 on kmet (unitless)")                # Koh 2025 Table 1: beta weight_k34 = 1.31, RSE 30.9%

    # -----------------------------------------------------------------
    # SA disposition (two compartments).
    lcl_sa     <- log(2.76); label("Apparent SA clearance CLm/F (L/h)")                      # Koh 2025 Table 1: CLm/F = 2.76 L/h, RSE 3.86%
    e_wt_cl_sa <- 1.42; label("Power exponent of WT/68.35 on cl_sa (unitless)")              # Koh 2025 Table 1: beta weight_CLm = 1.42, RSE 22.8%
    lvc_sa     <- log(7.5); label("Apparent SA central volume V4/F (L)")                     # Koh 2025 Table 1: V4/F = 7.5 L, RSE 2.6%
    lq_sa      <- fixed(log(0.08)); label("Apparent SA inter-compartmental clearance Q/F (L/h)")  # Koh 2025 Table 1: Q/F = 0.08 L/h, no RSE (fixed from Dings & Lehr, ref 14)
    lvp_sa     <- fixed(log(1.98)); label("Apparent SA peripheral volume V5/F (L)")               # Koh 2025 Table 1: V5/F = 1.98 L, no RSE (fixed from Dings & Lehr, ref 14)

    # -----------------------------------------------------------------
    # Thromboxane B2 turnover with Imax inhibition of production.
    # Koh 2025 Figure 1 / Appendix 1:
    #   dR/dt = kin * (1 - Imax * Cp^gam / (Cp^gam + IC50^gam)) - kout * R
    #   kin = R0 * kout ; R(0) = R0
    # Cp is the ASA (parent) plasma concentration.
    lrbase <- log(26.4); label("Baseline serum thromboxane B2 R0 (ug/L)")                    # Koh 2025 Table 1: R0 = 26.4 ug/L, RSE 7.67%
    lkout  <- log(0.023); label("First-order thromboxane B2 elimination rate kout (1/h)")    # Koh 2025 Table 1: kout = 0.023 1/h, RSE 5.51%
    limax  <- fixed(log(1)); label("Maximum fractional inhibition of TXB2 production Imax (unitless)")  # Koh 2025 Table 1: Imax = 1, no RSE (fixed)
    # NOTE ON UNITS: Koh 2025 Table 1 prints "IC50 (mol/L) 0.0036". That unit
    # label is wrong by a factor of 1e6 -- the correct unit is umol/L, i.e. the
    # "u" was dropped. The printed NUMERAL 0.0036 is preserved unchanged; only
    # the unit is corrected, and it is already in this model's umol/L base.
    # Read literally as mol/L the model predicts essentially ZERO TXB2
    # inhibition, contradicting the paper's abstract ("exceeding 80%") and its
    # own Figure 3.
    #
    # The concentration base is pinned independently of the PD model by Koh 2025
    # Figure 2, whose axes are in nmol/L: simulated in this model's umol/L base,
    # the median ASA and SA peaks land on the published Figure 2A/2B/2E/2F
    # curves in both height and time. With that base fixed, scoring the
    # candidate readings of the IC50 unit against the twelve published Figure 3D
    # percent-inhibition bars puts the minimum at exactly 0.0036 umol/L, and it
    # is the only reading that preserves the published spread between the
    # formulation and dose arms. Both checks are computed and asserted in the
    # vignette (Figure 2 section and Errata); see there for the numbers.
    lic50  <- fixed(log(0.0036)); label("ASA concentration giving half-maximal TXB2 inhibition IC50 (umol/L)")  # Koh 2025 Table 1: IC50 = 0.0036, no RSE (fixed from Kimura 2014, ref 19); printed unit "(mol/L)" corrected to umol/L -- see vignette Errata
    lhill  <- fixed(log(1)); label("Hill coefficient on the ASA-vs-TXB2 sigmoid (unitless)")                    # Koh 2025 Table 1: Gamma = 1, no RSE (fixed from Kimura 2014, ref 19)

    # -----------------------------------------------------------------
    # Inter-individual variability. Koh 2025 Table 1 footnote a: "Random
    # effect is presented as standard deviations", so each variance below
    # is the published omega squared. IIV is exponential (log-normal) on
    # every parameter except fr, which is logit-normal.
    etalogitfdepot ~ 0.0841   # Koh 2025 Table 1: omega fr    = 0.29 (logit scale), RSE 16.5% -> 0.29^2
    etalka         ~ 1.1881   # Koh 2025 Table 1: omega ka    = 1.09, RSE 14.7% -> 1.09^2
    etald1         ~ 0.8281   # Koh 2025 Table 1: omega Tk0   = 0.91, RSE 13.7% -> 0.91^2
    etaltlag       ~ 0.2601   # Koh 2025 Table 1: omega lag0  = 0.51, RSE 12.6% -> 0.51^2
    etalkmet       ~ 0.0729   # Koh 2025 Table 1: omega k34   = 0.27, RSE 15.0% -> 0.27^2
    etalcl_sa      ~ 0.0576   # Koh 2025 Table 1: omega CLm/F = 0.24, RSE 11.9% -> 0.24^2
    etalrbase      ~ 0.2401   # Koh 2025 Table 1: omega R0    = 0.49, RSE 12.1% -> 0.49^2

    # -----------------------------------------------------------------
    # Residual unexplained variability: proportional for both PK analytes,
    # additive for TXB2 (Koh 2025 Development of Structural Model).
    propSd     <- 0.41; label("ASA proportional residual SD (fraction)")   # Koh 2025 Table 1: proportional error (ASA) = 0.41, RSE 5.52%
    propSd_sa  <- 0.17; label("SA proportional residual SD (fraction)")    # Koh 2025 Table 1: proportional error (SA)  = 0.17, RSE 5.2%
    addSd_TXB2 <- 2.58; label("TXB2 additive residual SD (ug/L)")          # Koh 2025 Table 1: additive error (TXB2)    = 2.58, RSE 12.1%
  })

  model({
    # 1. Formulation effect on the first-order absorption rate constant.
    #    The tablet (FORM_CAPSULE = 0) is the reference; the capsule
    #    typical value is estimated separately rather than as a shift, so
    #    that each arm keeps its published fixed / estimated status.
    lka_form <- (1 - FORM_CAPSULE) * lka_tablet + FORM_CAPSULE * lka_capsule

    # 2. Individual parameters.
    ka     <- exp(lka_form + etalka)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    d1     <- exp(ld1 + etald1)
    tlag   <- exp(ltlag + etaltlag)

    k_presystemic_central    <- exp(lk_presystemic_central)
    k_presystemic_central_sa <- exp(lk_presystemic_central_sa)

    kmet  <- exp(lkmet + etalkmet) * (WT / 68.35)^e_wt_kmet
    vc    <- exp(lvc)

    cl_sa <- exp(lcl_sa + etalcl_sa) * (WT / 68.35)^e_wt_cl_sa
    vc_sa <- exp(lvc_sa)
    q_sa  <- exp(lq_sa)
    vp_sa <- exp(lvp_sa)

    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout)
    kin   <- rbase * kout
    imax  <- exp(limax)
    ic50  <- exp(lic50)
    hill  <- exp(lhill)

    # 3. Micro-constants for the SA sub-model (Koh 2025 Appendix 1:
    #    k40 = CLm/V4, k45 = Q/V4, k54 = Q/V5).
    k40 <- cl_sa / vc_sa
    k45 <- q_sa / vc_sa
    k54 <- q_sa / vp_sa

    # 4. ODE system, matching Koh 2025 Appendix 1 term for term.
    d/dt(depot)       <- -ka * depot
    d/dt(presystemic) <-  ka * depot -
      (k_presystemic_central + k_presystemic_central_sa) * presystemic
    d/dt(central)     <-  k_presystemic_central * presystemic - kmet * central
    d/dt(central_sa)  <-  k_presystemic_central_sa * presystemic + kmet * central -
      k40 * central_sa - k45 * central_sa + k54 * peripheral1_sa
    d/dt(peripheral1_sa) <- k45 * central_sa - k54 * peripheral1_sa

    # 5. Dual absorption. The first-order arm carries fraction fdepot
    #    through the depot; the zero-order arm delivers the complement
    #    directly into the pre-systemic compartment over d1 hours,
    #    starting after the tlag lag time.
    f(depot)         <- fdepot
    f(presystemic)   <- 1 - fdepot
    dur(presystemic) <- d1
    alag(presystemic) <- tlag

    # 6. Observations.
    Cc    <- central / vc
    Cc_sa <- central_sa / vc_sa

    inh <- imax * max(Cc, 0)^hill / (max(Cc, 0)^hill + ic50^hill)
    d/dt(txb2) <- kin * (1 - inh) - kout * txb2
    txb2(0)    <- rbase
    TXB2       <- txb2

    Cc    ~ prop(propSd)
    Cc_sa ~ prop(propSd_sa)
    TXB2  ~ add(addSd_TXB2)
  })
}
