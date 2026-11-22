Yu_2026_dexamethasone_horse_mpbpk <- function() {
  description <- paste(
    "Veterinary (horse). mPBPK-PD. Minimal physiologically based PK model of",
    "dexamethasone (DEX) in horses (blood, liver with saturable DEX binding,",
    "lumped remainder tissue, renal and hepatic clearance) with route- and",
    "prodrug-specific inputs (IV and IM DEX alcohol, DEX sodium phosphate and",
    "DEX isonicotinate; oral solid and liquid DEX absorbed into the liver;",
    "intra-articular DEX sodium phosphate into healthy or LPS-inflamed joints),",
    "linked to an indirect-response model of adrenal suppression (inhibition of",
    "a circadian cortisol production rate, plus a dose-triggered stress surge)",
    "and of DEX-stimulated plasma glucose production (Yu 2026). Fitted with",
    "ADAPT 5 to mean data pooled from 12 published studies; typical-value model",
    "with no between-subject or between-study variability. Model time is clock",
    "time in hours after midnight of the first study day."
  )
  reference <- paste(
    "Yu R, Toutain P-L, Ekstrand C, Jusko WJ. Meta-Analysis and",
    "Physiologically-Based Modeling of the Pharmacokinetics and",
    "Pharmacodynamics of Dexamethasone in Horses. Pharm Res. 2026.",
    "doi:10.1007/s11095-026-04120-5"
  )
  vignette <- "Yu_2026_dexamethasone_horse_mpbpk"
  dosing <- c(
    "a_blood",
    "depot_im",
    "blood_pho",
    "depot_im_pho",
    "blood_iso",
    "depot_im_iso",
    "depot_oral_solid",
    "depot_oral_liquid",
    "joint_pho",
    "joint_inflamed1",
    "joint_inflamed2",
    "joint_inflamed3"
  )
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  # Dosing compartments for the prodrug pools, the route-specific depots and
  # the joint spaces (Yu 2026 Figure 1, Eqs. 1-7), plus the drug-induced
  # deviation of plasma cortisol from its circadian baseline (Eq. 16 rewritten
  # in deviation form; see model()).
  paper_specific_compartments <- c(
    "blood_pho",
    "blood_iso",
    "depot_im_pho",
    "depot_im_iso",
    "depot_oral_solid",
    "depot_oral_liquid",
    "joint_pho",
    "joint",
    "joint_inflamed1",
    "joint_inflamed2",
    "joint_inflamed3",
    "cortisol_dev"
  )

  compartmentData <- list(
    blood_pho = list(
      analyte = "dexamethasone sodium phosphate (prodrug)",
      units = "ug",
      specimen = "whole blood",
      verified = TRUE
    ),
    blood_iso = list(
      analyte = "dexamethasone isonicotinate (prodrug)",
      units = "ug",
      specimen = "whole blood",
      verified = TRUE
    ),
    depot_im = list(analyte = "dexamethasone", units = "ug", specimen = "administration site", verified = TRUE),
    depot_im_pho = list(
      analyte = "dexamethasone sodium phosphate (prodrug)",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_im_iso = list(
      analyte = "dexamethasone isonicotinate (prodrug)",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_oral_solid = list(analyte = "dexamethasone", units = "ug", specimen = "administration site", verified = TRUE),
    depot_oral_liquid = list(
      analyte = "dexamethasone",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    joint_pho = list(
      analyte = "dexamethasone sodium phosphate (prodrug)",
      units = "ug",
      specimen = "synovial fluid",
      verified = TRUE
    ),
    joint = list(analyte = "dexamethasone", units = "ug", specimen = "synovial fluid", verified = TRUE),
    joint_inflamed1 = list(analyte = "dexamethasone", units = "ug", specimen = "synovial fluid", verified = TRUE),
    joint_inflamed2 = list(analyte = "dexamethasone", units = "ug", specimen = "synovial fluid", verified = TRUE),
    joint_inflamed3 = list(analyte = "dexamethasone", units = "ug", specimen = "synovial fluid", verified = TRUE),
    a_blood = list(analyte = "dexamethasone", units = "ug", specimen = "whole blood", verified = TRUE),
    a_liver = list(analyte = "dexamethasone", units = "ug", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "dexamethasone", units = "ug", specimen = "tissue", verified = TRUE),
    cortisol_dev = list(
      analyte = "cortisol (deviation from circadian baseline)",
      units = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    glucose = list(analyte = "glucose", units = "mmol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "horse",
    n_subjects = NA_integer_,
    n_studies = 12L,
    age_range = "2-22 years",
    weight_range = "385-700 kg",
    sex_female_pct = NA_real_,
    breeds = "Saddlebred, Standardbred, Thoroughbred, Quarter Horse, warmblood (Table I)",
    disease_state = paste(
      "Healthy horses; one study (Ekstrand 2019) in horses with",
      "lipopolysaccharide-induced radiocarpal synovitis receiving",
      "intra-articular doses"
    ),
    dose_range = paste(
      "DEX-equivalent 0.1 ug/kg IV bolus + infusion up to 0.056 mg/kg IV;",
      "0.024-0.05 mg/kg IM; 0.0225-0.05 mg/kg oral (single and 3-day",
      "daily); 0.05 mg/kg IA into healthy joints and 0.01-3 mg IA into",
      "inflamed joints; 0.05 mg/kg IM daily for 7 days (glucose only)"
    ),
    regions = "Europe, North America (published studies 1984-2020)",
    notes = paste(
      "Model-based meta-analysis of mean concentration-time data, either",
      "digitised from the published figures or the authors' original",
      "measurements (Table I). At least 6 horses per study except the LPS",
      "synovitis study (one horse per IA dose level). Females, stallions and",
      "geldings. Physiological volumes and flows are those of a typical",
      "550-kg horse (Table S1)."
    )
  )

  ini({
    # Structure: Yu 2026 Eqs. 1-17 and the ADAPT 5 code in Supplementary
    # Materials sections 1-3. Values: Table II (PK) and Table III (PD), except
    # where noted. Physiological constants: Table S1 and the ADAPT code (which
    # carries the unrounded values used in the fit).

    # ---- Physiology of a typical 550-kg horse (Table S1; ADAPT code) ----
    lvblood <- fixed(log(38.896)); label("Blood volume Vb (L)") # ADAPT code Vp = 38896 mL; Table S1 38.90 L
    lvliver <- fixed(log(7.15)); label("Liver volume Vl (L)") # ADAPT code V1 = 7150 mL; Table S1 7.15 L
    lvremainder <- fixed(log(503.954)); label("Remainder-tissue volume Vr (L)") # ADAPT code V2 = 503954 mL; Table S1 504 L
    lqliver <- fixed(log(651.55)); label("Total liver blood flow Ql (L/h)") # ADAPT code Q1 = 651550 mL/h; Table S1 652 L/h
    lqremainder <- fixed(log(1311.99)); label("Remainder-tissue blood flow Qr (L/h)") # ADAPT code Q2 = 1311990 mL/h; Table S1 1312 L/h
    bpr <- fixed(0.6945); label("Blood-to-plasma concentration ratio Rb (unitless)") # ADAPT code Rb = 0.6945; Methods text 0.69 (rat value)
    lqurine <- fixed(log(0.494)); label("Mean urine flow Vu (L/h)") # Methods Eq. 12: 494 mL/h

    # ---- Shared DEX disposition (Table II) ----
    lbmax_liver <- log(140.1); label("DEX binding capacity in liver Bmax,l (ng/mL)") # Table II: 140.1 (CV 12%)
    lkd_liver <- log(0.269); label("Equilibrium dissociation constant of DEX in liver Kd (ng/mL)") # Table II: 0.269 (CV 16%)
    lkp_remainder <- log(1.09); label("Remainder-tissue partition coefficient Kp,r (unitless)") # Table II: 1.09 (CV 5%)
    lclint <- log(473.6); label("Intrinsic hepatic clearance CLh,int per kg body weight (mL/h/kg)") # Table II: 473.6 (CV 5%)
    lclr <- log(5.69); label("Renal clearance CLr per kg body weight (mL/h/kg)") # Table II: 5.69 (CV 10%)

    # ---- Prodrug conversion in blood after IV (and IM) prodrug doses (Table II) ----
    lkf_pho <- log(16.67); label("Formation rate constant of DEX from DEX sodium phosphate in blood kf,VM_PHO (1/h)") # Table II: 16.67 (CV 34%)
    lkf_iso <- log(20.59); label("Formation rate constant of DEX from DEX isonicotinate in blood kf,VM_ISO (1/h)") # Table II: 20.59 (CV 27%)
    lfconv_pho <- fixed(log(1)); label("Prodrug conversion bioavailability of DEX sodium phosphate Fc_PHO (fraction)") # Table II: 100% (fixed)
    lfconv_iso <- log(0.695); label("Prodrug conversion bioavailability of DEX isonicotinate Fc_ISO (fraction)") # Table II: 69.5% (CV 8%)

    # ---- IM absorption (Table II) ----
    lka_im <- log(0.131); label("IM absorption rate constant of DEX alcohol ka,IM_DEX (1/h)") # Table II: 0.131 (CV 19%)
    lka_im_pho <- log(2.95); label("IM absorption rate constant of DEX sodium phosphate ka,IM_PHO (1/h)") # Table II: 2.95 (CV 56%)
    lka_im_iso <- log(0.0105); label("IM absorption rate constant of DEX isonicotinate ka,IM_ISO (1/h)") # Table II: 0.0105 (CV 5.5%)
    lfdepot_im <- fixed(log(1)); label("IM bioavailability of DEX alcohol (fraction)") # ADAPT PD code Fml = 1; Results: 'negligible absorption loss'

    # ---- Oral absorption into the liver (Table II) ----
    lka_oral_solid <- log(0.421); label("Oral absorption rate constant, powder or tablet ka,PO (1/h)") # Table II: 0.421 (CV 11%)
    lka_oral_liquid <- log(0.529); label("Oral absorption rate constant, suspension or solution ka,PO (1/h)") # Table II: 0.529 (CV 11%)
    lfgi_solid <- log(0.521); label("Gastrointestinal bioavailability, powder or tablet FGI (fraction)") # Table II: 52.1% (CV 8%)
    lfgi_liquid <- log(0.682); label("Gastrointestinal bioavailability, suspension or solution FGI (fraction)") # Table II: 68.2% (CV 9%)

    # ---- Intra-articular dosing (Table II; Results text for the inflamed joint) ----
    lkf_ia_pho <- fixed(log(30)); label("Formation rate constant of DEX from DEX sodium phosphate in the joint kf,IA_PHO (1/h)") # Table II: 30 (fixed)
    lka_ia <- log(2.61); label("Absorption rate constant of DEX from healthy joints ka,IA (1/h)") # Table II: 2.61 (CV 24%)
    lfia <- log(0.830); label("Bioavailability after IA injection into healthy joints FIA (fraction)") # Table II: 83.0% (CV 11%)
    lka_ia_inflamed <- log(5.13); label("Absorption rate constant of DEX from inflamed joints ka,IA (1/h)") # Table II: 5.13 (CV 7%)
    lfia_inflamed <- fixed(log(1)); label("Bioavailability after IA injection into inflamed joints FIA (fraction)") # Table II: 100% (fixed)
    ffrac_sf1 <- 0.9783; label("Fractional intercept f1 of the inflamed-joint synovial-fluid DEX disposition (fraction)") # Results (IA paragraph), Eq. 7
    ffrac_sf2 <- 0.02098; label("Fractional intercept f2 of the inflamed-joint synovial-fluid DEX disposition (fraction)") # Results (IA paragraph), Eq. 7
    ffrac_sf3 <- 0.0007512; label("Fractional intercept f3 of the inflamed-joint synovial-fluid DEX disposition (fraction)") # Results (IA paragraph), Eq. 7
    lalpha_sf <- log(6.354); label("First disposition slope alpha of inflamed-joint synovial-fluid DEX (1/h)") # Results (IA paragraph), Eq. 7
    lbeta_sf <- log(0.5977); label("Second disposition slope beta of inflamed-joint synovial-fluid DEX (1/h)") # Results (IA paragraph), Eq. 7
    lgamma_sf <- log(0.07056); label("Third disposition slope gamma of inflamed-joint synovial-fluid DEX (1/h)") # Results (IA paragraph), Eq. 7
    lvsf <- fixed(log(0.0169)); label("Synovial fluid volume of the radiocarpal joint plus injected volume VSF (L)") # Methods: 16.9 mL

    # ---- Cortisol: circadian baseline and adrenal suppression (Table III) ----
    lrbase_cortisol <- log(62.25); label("Mesor (mean baseline) of plasma cortisol Rm (ng/mL)") # Table III: 62.25 (CV 4.8%), Ekstrand 2019; used for the Fig. 8 prediction
    lamp_cortisol <- log(22.01); label("Amplitude of the plasma cortisol rhythm Ra (ng/mL)") # Table III: 22.01 (CV 14%), Ekstrand 2019; used for the Fig. 8 prediction
    tacro_cortisol <- 8.94; label("Acrophase of the cortisol rhythm as clock time tp (h after midnight)") # ADAPT PD code tp1 = 0.94 h after an 8 AM start; Table III 8:56 AM (CV 61%)
    lkout_cortisol <- log(0.405); label("First-order loss rate constant of cortisol kout_CTS (1/h)") # Table III: 0.405 (CV 4.9%)
    lic50_cortisol <- log(0.0383); label("Plasma DEX concentration for 50% inhibition of cortisol production IC50 (ng/mL)") # Table III: 0.0383 (CV 8.9%)
    imax_cortisol <- fixed(1); label("Maximal inhibition of cortisol production Imax (fraction)") # Table III: 1.0 (fixed)
    lrmax_stress <- log(50.53); label("Stress-induced maximal cortisol surge after dosing Rmax,stress (ng/mL)") # Results (Fig. 8 paragraph): mean 50.53 of the three study estimates in Table III
    ltau_stress <- fixed(log(0.34)); label("Duration of the stress-induced cortisol surge tau (h)") # Table III: 0.34 (fixed)

    # ---- Glucose: DEX-stimulated production (Table III) ----
    # The ADAPT glucose code labels P(2) 'EC50' and P(3) 'Emax' (Psym) but
    # uses Emax = P(2) and EC50 = P(3) in DIFFEQ, so the fitted Emax is the
    # value Table III prints as EC50 and vice versa. With this assignment the
    # model reproduces all four Figure 9 profiles; with the Table III labels
    # it overshoots every one of them (see the vignette).
    lrbase_glucose <- log(5.14); label("Baseline plasma glucose E0 (mmol/L)") # Results: maximal GLU 11.24 mM = E0 * (1 + 1.187) -> E0 = 5.14; Table III range 4.65-5.47
    lkout_glucose <- log(0.0485); label("First-order loss rate constant of glucose kout_GLU (1/h)") # Table III: 0.0485 (CV 17.8%)
    lemax_glucose <- log(0.799); label("Maximal stimulation of glucose production Emax (fraction)") # Table III value printed as EC50 0.799 (CV 21%); ADAPT DIFFEQ Emax = P(2)
    lec50_glucose <- log(1.187); label("Plasma DEX concentration for 50% maximal stimulation EC50 (ng/mL)") # Table III value printed as Emax 1.187 (CV 51%); ADAPT DIFFEQ EC50 = P(3)

    # ---- Residual error ----
    # ADAPT variance model V = (sigma1 + sigma2 * Y)^2 (Eq. 18) per output group;
    # the variance-parameter estimates are not reported.
    propSd <- fixed(0); label("Proportional residual error, plasma DEX (fraction; not reported)") # Eq. 18; value not reported
    addSd <- fixed(0); label("Additive residual error, plasma DEX (ng/mL; not reported)") # Eq. 18; value not reported
    propSd_Curine <- fixed(0); label("Proportional residual error, urine DEX (fraction; not reported)") # Eq. 18; value not reported
    addSd_Curine <- fixed(0); label("Additive residual error, urine DEX (ng/mL; not reported)") # Eq. 18; value not reported
    propSd_Csynovial <- fixed(0); label("Proportional residual error, synovial fluid DEX (fraction; not reported)") # Eq. 18; value not reported
    addSd_Csynovial <- fixed(0); label("Additive residual error, synovial fluid DEX (ng/mL; not reported)") # Eq. 18; value not reported
    propSd_cortisol <- fixed(0); label("Proportional residual error, plasma cortisol (fraction; not reported)") # Eq. 18; value not reported
    addSd_cortisol <- fixed(0); label("Additive residual error, plasma cortisol (ng/mL; not reported)") # Eq. 18; value not reported
    propSd_glucose <- fixed(0); label("Proportional residual error, plasma glucose (fraction; not reported)") # Eq. 18; value not reported
    addSd_glucose <- fixed(0); label("Additive residual error, plasma glucose (mmol/L; not reported)") # Eq. 18; value not reported
  })

  model({
    # Molecular weights used to convert prodrug mass to DEX mass (Methods)
    mw_dex <- 392.46
    mw_pho <- 516.41
    mw_iso <- 492.54
    # Typical body weight behind Table S1 and the per-kg clearances (kg)
    bw_typ <- 550

    vblood <- exp(lvblood)
    vliver <- exp(lvliver)
    vremainder <- exp(lvremainder)
    qliver <- exp(lqliver)
    qremainder <- exp(lqremainder)
    qco <- qliver + qremainder
    qurine <- exp(lqurine)

    bmax_liver <- exp(lbmax_liver)
    kd_liver <- exp(lkd_liver)
    kp_remainder <- exp(lkp_remainder)
    # mL/h/kg -> L/h for the 550-kg horse (ADAPT code CL10 = 260500 mL/h, CLr = 3131 mL/h)
    clint <- exp(lclint) * bw_typ / 1000
    clr <- exp(lclr) * bw_typ / 1000

    kf_pho <- exp(lkf_pho)
    kf_iso <- exp(lkf_iso)
    fconv_pho <- exp(lfconv_pho)
    fconv_iso <- exp(lfconv_iso)
    ka_im <- exp(lka_im)
    ka_im_pho <- exp(lka_im_pho)
    ka_im_iso <- exp(lka_im_iso)
    fdepot_im <- exp(lfdepot_im)
    ka_oral_solid <- exp(lka_oral_solid)
    ka_oral_liquid <- exp(lka_oral_liquid)
    fgi_solid <- exp(lfgi_solid)
    fgi_liquid <- exp(lfgi_liquid)
    kf_ia_pho <- exp(lkf_ia_pho)
    ka_ia <- exp(lka_ia)
    fia <- exp(lfia)
    ka_ia_inflamed <- exp(lka_ia_inflamed)
    fia_inflamed <- exp(lfia_inflamed)
    alpha_sf <- exp(lalpha_sf)
    beta_sf <- exp(lbeta_sf)
    gamma_sf <- exp(lgamma_sf)
    vsf <- exp(lvsf)

    # Plasma DEX concentration (ng/mL = ug/L). The ADAPT code integrates the
    # plasma concentration with the blood equation divided by Vb * Rb, so the
    # blood amount is Cp * Vb * Rb and every blood-flow and renal term acts on
    # Cp (DIFFEQ XP(2): '- Qco*X(2) - CLr*X(2))/Vp/Rb').
    Cc <- a_blood / (vblood * bpr)
    cliver <- a_liver / vliver
    cremainder <- a_remainder / vremainder

    # Unbound fraction of DEX in liver from saturable binding (Eq. 11; ADAPT
    # Kp1 = 1/fut), written in the cancellation-free form of each sign branch
    # of a = Cl - Kd - Bmax so the low-concentration limit Kd / (Kd + Bmax)
    # is exact.
    a_bind <- cliver - kd_liver - bmax_liver
    sq_bind <- sqrt(a_bind * a_bind + 4 * kd_liver * cliver)
    if (a_bind >= 0) {
      fut <- 0.5 * (a_bind + sq_bind) / cliver
    } else {
      fut <- 2 * kd_liver / (sq_bind - a_bind)
    }

    # Inputs to blood (Eq. 8) and liver (Eq. 9)
    input_iv_im <- kf_pho * fconv_pho * blood_pho * mw_dex / mw_pho +
      kf_iso * fconv_iso * blood_iso * mw_dex / mw_iso +
      ka_im * fdepot_im * depot_im
    input_ia <- fia * ka_ia * joint +
      fia_inflamed * ka_ia_inflamed * (joint_inflamed1 + joint_inflamed2 + joint_inflamed3)
    input_po <- fgi_solid * ka_oral_solid * depot_oral_solid +
      fgi_liquid * ka_oral_liquid * depot_oral_liquid

    # Prodrugs in blood after IV or IM dosing (Eqs. 1 and 3)
    d/dt(blood_pho) <- ka_im_pho * depot_im_pho - kf_pho * blood_pho
    d/dt(blood_iso) <- ka_im_iso * depot_im_iso - kf_iso * blood_iso
    # IM depots (Eq. 2; DEX alcohol from the ADAPT code, Dose 14)
    d/dt(depot_im) <- -ka_im * depot_im
    d/dt(depot_im_pho) <- -ka_im_pho * depot_im_pho
    d/dt(depot_im_iso) <- -ka_im_iso * depot_im_iso
    # Oral depots (Eq. 4)
    d/dt(depot_oral_solid) <- -ka_oral_solid * depot_oral_solid
    d/dt(depot_oral_liquid) <- -ka_oral_liquid * depot_oral_liquid
    # Healthy joint: prodrug conversion then DEX absorption (Eqs. 5-6)
    d/dt(joint_pho) <- -kf_ia_pho * joint_pho
    d/dt(joint) <- kf_ia_pho * joint_pho * mw_dex / mw_pho - ka_ia * joint
    # Inflamed joint: DEX amount in synovial fluid follows the tri-exponential
    # Eq. 7, carried as three parallel first-order pools whose bioavailability
    # fractions are the intercepts f1-f3 (prodrug dose -> DEX-equivalent dose).
    d/dt(joint_inflamed1) <- -alpha_sf * joint_inflamed1
    d/dt(joint_inflamed2) <- -beta_sf * joint_inflamed2
    d/dt(joint_inflamed3) <- -gamma_sf * joint_inflamed3
    f(joint_inflamed1) <- ffrac_sf1 * mw_dex / mw_pho
    f(joint_inflamed2) <- ffrac_sf2 * mw_dex / mw_pho
    f(joint_inflamed3) <- ffrac_sf3 * mw_dex / mw_pho

    # Blood, liver and remainder tissue (Eqs. 8-10 as coded in ADAPT)
    d/dt(a_blood) <- input_iv_im + input_ia +
      qliver * cliver * fut + qremainder * cremainder / kp_remainder -
      qco * Cc - clr * Cc
    d/dt(a_liver) <- input_po + qliver * (Cc - cliver * fut) - clint * cliver * fut
    d/dt(a_remainder) <- qremainder * (Cc - cremainder / kp_remainder)

    # Urine concentration at collection midpoints (Eq. 12)
    Curine <- clr * Cc / qurine
    # Synovial-fluid DEX concentration (ng/mL = ug/L)
    Csynovial <- (joint + joint_inflamed1 + joint_inflamed2 + joint_inflamed3) / vsf

    # ---- Cortisol (Eqs. 13-16) ----
    rbase_cortisol <- exp(lrbase_cortisol)
    amp_cortisol <- exp(lamp_cortisol)
    kout_cortisol <- exp(lkout_cortisol)
    ic50_cortisol <- exp(lic50_cortisol)
    rmax_stress <- exp(lrmax_stress)
    tau_stress <- exp(ltau_stress)
    omega_circ <- 2 * pi / 24
    # Circadian baseline (Eq. 13) and the production rate that sustains it (Eq. 15)
    cortisol_bl <- rbase_cortisol + amp_cortisol * cos(omega_circ * (t - tacro_cortisol))
    kin_cortisol <- kout_cortisol * rbase_cortisol +
      kout_cortisol * amp_cortisol * cos(omega_circ * (t - tacro_cortisol)) -
      omega_circ * amp_cortisol * sin(omega_circ * (t - tacro_cortisol))
    # Stress surge Rmax,stress / tau for tau hours after a dose (Eq. 16);
    # tad() is NA before the first dose, so the surge is off until then.
    tad_dose <- tad()
    stress <- 0
    if (tad_dose >= 0 && tad_dose < tau_stress) {
      stress <- rmax_stress / tau_stress
    }
    inh_cortisol <- imax_cortisol * Cc / (Cc + ic50_cortisol)
    # Eq. 16 with R = cortisol_bl + cortisol_dev: substituting Eq. 15 gives
    # d(cortisol_dev)/dt = stress - kin(t) * I(Cp) - kout * cortisol_dev, with
    # cortisol_dev(0) = 0 (R starts on its baseline, Eq. 16 initial condition).
    d/dt(cortisol_dev) <- stress - kin_cortisol * inh_cortisol - kout_cortisol * cortisol_dev
    cortisol <- cortisol_bl + cortisol_dev

    # ---- Glucose (Eq. 17) ----
    rbase_glucose <- exp(lrbase_glucose)
    kout_glucose <- exp(lkout_glucose)
    emax_glucose <- exp(lemax_glucose)
    ec50_glucose <- exp(lec50_glucose)
    kin_glucose <- rbase_glucose * kout_glucose
    glucose(0) <- rbase_glucose
    d/dt(glucose) <- kin_glucose * (1 + emax_glucose * Cc / (Cc + ec50_glucose)) - kout_glucose * glucose

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Curine ~ add(addSd_Curine) + prop(propSd_Curine) + combined1()
    Csynovial ~ add(addSd_Csynovial) + prop(propSd_Csynovial) + combined1()
    cortisol ~ add(addSd_cortisol) + prop(propSd_cortisol) + combined1()
    glucose ~ add(addSd_glucose) + prop(propSd_glucose) + combined1()
  })
}
