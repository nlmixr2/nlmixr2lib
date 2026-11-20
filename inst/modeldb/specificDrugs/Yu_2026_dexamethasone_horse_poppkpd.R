Yu_2026_dexamethasone_horse_poppkpd <- function() {
  description <- paste(
    "Veterinary (horse). mPBPK-PD. Population minimal physiologically based",
    "PK model of dexamethasone (DEX) in horses (blood, liver with saturable",
    "DEX binding, lumped remainder tissue, renal and hepatic clearance, and a",
    "well-stirred bladder for urine concentrations) with route- and",
    "prodrug-specific inputs (IV DEX alcohol, IV DEX isonicotinate and DEX",
    "sodium phosphate, IM DEX alcohol and DEX isonicotinate, intra-articular",
    "DEX sodium phosphate), linked to an indirect-response model of adrenal",
    "suppression (Hill-type inhibition of a circadian cortisol production",
    "rate plus a post-dose stress surge). Prior-informed FOCEI (NONMEM",
    "FOCEI-Priors) fit to individual data from 32 horses in five studies",
    "(Yu 2026). Physiological volumes and flows scale linearly with body",
    "weight; hepatic and renal clearance scale with (WT/500)^0.75, and",
    "hepatic clearance is 297 L/h in mares versus 242.8 L/h in geldings at",
    "500 kg. Model time is",
    "clock time in hours after midnight of the first study day."
  )
  reference <- paste(
    "Yu R, Toutain P-L, Ekstrand C, Jusko WJ. Prior-informed population",
    "pharmacokinetic-pharmacodynamic modeling of dexamethasone in horses.",
    "J Pharmacokinet Pharmacodyn. 2026. doi:10.1007/s10928-026-10062-7"
  )
  vignette <- "Yu_2026_dexamethasone_horse_poppkpd"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  # Prodrug pools, route-specific depots, the joint space and the bladder
  # (Yu 2026 Figure 1, Eqs. 1-5 and 10), plus the drug-induced deviation of
  # plasma cortisol from its circadian baseline (Eq. 14 rewritten in
  # deviation form; see model()).
  paper_specific_compartments <- c(
    "blood_pho",
    "blood_iso",
    "depot_im_iso",
    "joint_pho",
    "joint",
    "bladder",
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
    depot_im_iso = list(
      analyte = "dexamethasone isonicotinate (prodrug)",
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
    a_blood = list(analyte = "dexamethasone", units = "ug", specimen = "whole blood", verified = TRUE),
    a_liver = list(analyte = "dexamethasone", units = "ug", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "dexamethasone", units = "ug", specimen = "tissue", verified = TRUE),
    bladder = list(analyte = "dexamethasone", units = "ug", specimen = "urine", verified = TRUE),
    cortisol_dev = list(
      analyte = "cortisol (deviation from circadian baseline)",
      units = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every physiological volume and flow linearly (NONMEM $PK, per-kg",
        "constants), and hepatic and renal clearance as (WT/500)^0.75. Observed",
        "range 430-700 kg (Table 1)."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female (mare), 0 = gelding",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the category-specific hepatic clearance: 297 L/h for mares and",
        "242.8 L/h for geldings at 500 kg (Table 2). Source column SEX with",
        "0 = female and 2 = gelding: SEXF = as.integer(SEX == 0). Intact males",
        "(stallions) were not in the data; the gelding value applies to SEXF = 0."
      ),
      source_name = "SEX"
    ),
    STUDY_TOUTAIN1984 = list(
      description = paste(
        "1 = horse from Study 1 (Toutain 1984; IV DEX and IV DEX isonicotinate",
        "0.05 mg/kg, conducted 1981), 0 = any other study"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Switches on the stress-induced cortisol surge Rmax,stress / tau for",
        "tau = 0.34 h after each dose; the PD control stream applies the surge",
        "only when STUDY_ID = 1. Source column STUDY_ID:",
        "STUDY_TOUTAIN1984 = as.integer(STUDY_ID == 1)."
      ),
      source_name = "STUDY_ID"
    ),
    STUDY_EKSTRAND2015 = list(
      description = paste(
        "1 = horse from Study 2 (Ekstrand 2015; IM DEX isonicotinate",
        "0.0237 mg/kg, conducted 2010), 0 = any other study"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the study-specific cortisol mesor Rm1 = 93 ng/mL instead of",
        "Rm2 = 47 ng/mL (Table 3). Source column STUDY_ID:",
        "STUDY_EKSTRAND2015 = as.integer(STUDY_ID == 2)."
      ),
      source_name = "STUDY_ID"
    )
  )

  population <- list(
    species = "horse",
    n_subjects = 32L,
    n_studies = 5L,
    age_range = "3-20 years",
    weight_range = "430-700 kg",
    sex_female_pct = 31.25,
    breeds = "Saddlebred, Standardbred (Table 1)",
    disease_state = "Healthy adult horses (mares and geldings)",
    dose_range = paste(
      "DEX-equivalent IV bolus and/or 3-h infusion 0.17 ug/kg to 0.05 mg/kg",
      "(DEX alcohol, DEX isonicotinate, DEX sodium phosphate); IM DEX",
      "isonicotinate 0.0237 mg/kg; IA DEX sodium phosphate 0.01-3 mg"
    ),
    regions = "France, Sweden (studies conducted 1981-2013)",
    notes = paste(
      "Individual data from five studies by one research network (Table 1;",
      "refs 8-12): 623 plasma DEX, 95 urine DEX and 1063 plasma cortisol",
      "concentrations. 10 mares and 22 geldings. Sequential PK then PD",
      "estimation in NONMEM 7.5 with informative priors from the authors'",
      "naive-pooled meta-analysis (Yu 2026, Pharm Res)."
    )
  )

  ini({
    # Structure: Yu 2026 Eqs. 1-14 and the NONMEM control streams in the
    # Supplementary Materials (DEXPKC3 for PK, DEXPDC8 for PD). Values:
    # Table 2 (PK) and Table 3 (PD), FOCEI-Priors column. Physiological
    # constants: the per-kg values in the NONMEM $PK block (Table S1 lists
    # the same values for a 550-kg horse).

    # ---- Physiology per kg body weight (NONMEM $PK) ----
    lvblood <- fixed(log(0.07072)); label("Blood volume per kg body weight Vb/BW (L/kg)") # $PK Vp = 70.72*BW mL; Table S1 38.90 L at 550 kg
    lvliver <- fixed(log(0.013045)); label("Liver volume per kg body weight Vl/BW (L/kg)") # $PK V1 = 13.045*BW mL; Table S1 7.15 L at 550 kg
    lvremainder <- fixed(log(0.91628)); label("Remainder-tissue volume per kg body weight Vr/BW (L/kg)") # $PK V2 = 916.28*BW mL; Table S1 504 L at 550 kg
    lqco <- fixed(log(3.570)); label("Cardiac output per kg body weight Qco/BW (L/h/kg)") # $PK Qco = 3570*BW mL/h; Table S1 1964 L/h at 550 kg
    lqliver <- fixed(log(1.18463)); label("Liver blood flow per kg body weight Ql/BW (L/h/kg)") # $PK Q1 = 1184.63*BW mL/h; Table S1 652 L/h at 550 kg
    lqremainder <- fixed(log(2.39544)); label("Remainder-tissue blood flow per kg body weight Qr/BW (L/h/kg)") # $PK Q2 = 2395.44*BW mL/h (Table S1 1312 L/h at 550 kg is 2385.45 mL/h/kg)
    lqurine <- fixed(log(0.00124)); label("Urine flow per kg body weight CLu/BW (L/h/kg)") # $PK CLu = 1.24*BW mL/h; Table S1 0.682 L/h at 550 kg
    lvbladder <- fixed(log(4)); label("Constant effective bladder volume Vbc (L)") # $PK Vbc = 4000 mL; Table S1 4.0 L

    # ---- DEX disposition (Table 2) ----
    lclint <- log(242.8); label("Hepatic (intrinsic) clearance in a 500-kg gelding CLh_gelding (L/h)") # Table 2: 242.8 L/h (RSE 0.2%)
    lclint_female <- log(297); label("Hepatic (intrinsic) clearance in a 500-kg mare CLh_female (L/h)") # Table 2: 297 L/h (RSE 0.2%)
    lclr <- log(3.430); label("Renal clearance in a 500-kg horse CLr (L/h)") # Table 2: 3.430 L/h (RSE 0.7%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on hepatic and renal clearance (unitless)") # Table 2: 0.75 FIX
    lbmax_liver <- log(121); label("DEX binding capacity in liver Bmax_l (ng/mL)") # Table 2: 121 ng/mL (RSE 8.8%)
    lkd_liver <- log(0.387); label("Equilibrium dissociation constant of DEX in liver Kd_l (ng/mL)") # Table 2: 0.387 ng/mL (RSE 8.3%)
    lkp_remainder <- log(1.04); label("Remainder-tissue partition coefficient Kp_r (unitless)") # Table 2: 1.04 (RSE 1.9%)
    bpr <- 0.852; label("Blood-to-plasma concentration ratio Rb (unitless)") # Table 2: 0.852 (RSE 1.3%)

    # ---- Prodrug conversion in blood (Table 2) ----
    lkf_iso <- log(22.8); label("Formation rate constant of DEX from DEX isonicotinate in blood kf,b_ISO (1/h)") # Table 2: 22.8 (RSE 16.8%)
    lfconv_iso <- log(0.812); label("Conversion bioavailability of DEX isonicotinate in blood Ff,b_ISO (fraction)") # Table 2: 81.2% (RSE 10.3%)
    lkf_pho <- fixed(log(16.7)); label("Formation rate constant of DEX from DEX sodium phosphate in blood kf,b_PHO (1/h)") # Table 2: 16.7 FIX
    lfconv_pho <- fixed(log(1)); label("Conversion bioavailability of DEX sodium phosphate in blood Ff,b_PHO (fraction)") # Table 2: 100% FIX

    # ---- IM absorption (Table 2) ----
    lka_im_iso <- log(0.0115); label("IM absorption rate constant of DEX isonicotinate ka,m_ISO (1/h)") # Table 2: 0.0115 (RSE 2.2%)
    lka_im <- log(0.131); label("IM absorption rate constant of DEX alcohol ka,m_DEX (1/h)") # Table 2: 0.131 (RSE 19%)

    # ---- Intra-articular dosing (Table 2) ----
    lkf_ia_pho <- fixed(log(30)); label("Formation rate constant of DEX from DEX sodium phosphate in the joint kf,j_PHO (1/h)") # Table 2: 30 FIX
    lka_ia <- log(4.94); label("Absorption rate constant of DEX from the joint ka,j_DEX (1/h)") # Table 2: 4.94 (RSE 6.9%)
    lfia <- fixed(log(1)); label("Bioavailability of DEX from the joint Fj_DEX (fraction)") # Table 2: 100% FIX

    # ---- Cortisol: circadian baseline and adrenal suppression (Table 3) ----
    lrbase_cortisol <- log(47); label("Cortisol mesor Rm2, studies other than Ekstrand 2015 (ng/mL)") # Table 3: Rm2 = 47 ng/mL (RSE 1.7%)
    lrbase_cortisol_ekstrand2015 <- log(93); label("Cortisol mesor Rm1, Ekstrand 2015 IM study (ng/mL)") # Table 3: Rm1 = 93 ng/mL (RSE 0.5%)
    lamp_cortisol <- log(17.6); label("Amplitude of the plasma cortisol rhythm Ra (ng/mL)") # Table 3: 17.6 ng/mL (RSE 4.6%)
    tacro_cortisol <- 8.992; label("Acrophase of the cortisol rhythm as clock time Tp (h after midnight)") # Table 3: Tp = 0.992 h after 8 AM (RSE 13.8%); PD $DES peak at tp + 8
    lkout_cortisol <- log(0.246); label("First-order loss rate constant of cortisol kout (1/h)") # Table 3: 0.246 (RSE 3.3%)
    lic50_cortisol <- log(0.0412); label("Plasma DEX concentration for 50% inhibition of cortisol production IC50 (ng/mL)") # Table 3: 0.0412 ng/mL (RSE 4.7%)
    imax_cortisol <- fixed(1); label("Maximal inhibition of cortisol production Imax (fraction)") # Table 3: 1 FIX
    hill_cortisol <- 1.19; label("Hill coefficient of cortisol inhibition gamma (unitless)") # Table 3 and Table S3: 1.19 (RSE 5%)
    lrmax_stress <- log(24.7); label("Stress-induced maximal cortisol surge after dosing Rmax,stress (ng/mL)") # Table 3: 24.7 ng/mL (RSE 5.7%)
    ltau_stress <- fixed(log(0.34)); label("Duration of the stress-induced cortisol surge tau (h)") # Table 3: 0.34 FIX

    # ---- Between-subject variability ----
    # The PK control stream applies ETA(1) to both CLh and Kd
    # ('Kd = TVKd* EXP(ETA(1))'); Table 2 reports it as the BSV in CLh.
    etalclint ~ 0.094 # Table 2: omega2 CLh = 0.094 (RSE 15.3%, shrinkage 13%)
    etalrbase_cortisol ~ 0.12 # Table 3: omega2 Rm = 0.12 (RSE 14.4%, shrinkage 5.1%)

    # ---- Residual error (tables report sqrt(SIGMA)) ----
    propSd <- 0.45; label("Proportional residual error, plasma DEX (fraction)") # Table 2: 0.45 (Results: 45%)
    propSd_Curine <- 5.44; label("Proportional residual error, urine DEX (fraction)") # Table 2: 5.44 (Results: 544%)
    addSd_cortisol <- 20; label("Additive residual error, plasma cortisol (ng/mL)") # Table 3: 20 ng/mL
  })

  model({
    # Molecular weights used to convert prodrug mass to DEX mass ($PK MWd, MWo)
    mw_dex <- 392.46
    mw_pho <- 516.41
    mw_iso <- 492.54

    # Physiology scaled linearly with body weight ($PK)
    vblood <- exp(lvblood) * WT
    vliver <- exp(lvliver) * WT
    vremainder <- exp(lvremainder) * WT
    qco <- exp(lqco) * WT
    qliver <- exp(lqliver) * WT
    qremainder <- exp(lqremainder) * WT
    qurine <- exp(lqurine) * WT
    vbladder <- exp(lvbladder)

    # Sex-specific hepatic clearance with power-law body weight (Methods
    # equation; $PK IF (SEX.EQ.0) ... CLh = TVCLhF * (BW/500)**0.75 * EXP(ETA(1)))
    clint <- exp((1 - SEXF) * lclint + SEXF * lclint_female + etalclint) * (WT / 500)^e_wt_cl
    clr <- exp(lclr) * (WT / 500)^e_wt_cl
    bmax_liver <- exp(lbmax_liver)
    kd_liver <- exp(lkd_liver + etalclint)
    kp_remainder <- exp(lkp_remainder)

    kf_iso <- exp(lkf_iso)
    fconv_iso <- exp(lfconv_iso)
    kf_pho <- exp(lkf_pho)
    fconv_pho <- exp(lfconv_pho)
    ka_im_iso <- exp(lka_im_iso)
    ka_im <- exp(lka_im)
    kf_ia_pho <- exp(lkf_ia_pho)
    ka_ia <- exp(lka_ia)
    fia <- exp(lfia)

    # Plasma DEX concentration (ng/mL = ug/L). The NONMEM code integrates the
    # plasma concentration with the blood equation divided by Vp * Rb, so the
    # blood amount is Cp * Vb * Rb and every blood-flow and renal term acts on
    # Cp ('- Qco * A(3) - CLr * A(3)) / Vp / Rb').
    Cc <- a_blood / (vblood * bpr)
    cliver <- a_liver / vliver
    cremainder <- a_remainder / vremainder

    # Unbound fraction of DEX in liver from saturable binding (Eq. 9; $DES
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

    # Inputs to blood (Eq. 6 input terms)
    input_iv_im <- kf_iso * fconv_iso * blood_iso * mw_dex / mw_iso +
      kf_pho * fconv_pho * blood_pho * mw_dex / mw_pho +
      ka_im * depot_im
    input_ia <- fia * ka_ia * joint

    # Prodrugs in blood after IV or IM dosing (Eqs. 1 and 3)
    d/dt(blood_iso) <- ka_im_iso * depot_im_iso - kf_iso * blood_iso
    d/dt(blood_pho) <- -kf_pho * blood_pho
    # IM depots (Eq. 2); DEX alcohol is absorbed straight into blood ($PK DOSE 4)
    d/dt(depot_im_iso) <- -ka_im_iso * depot_im_iso
    d/dt(depot_im) <- -ka_im * depot_im
    # Joint: prodrug conversion then DEX absorption (Eqs. 4-5)
    d/dt(joint_pho) <- -kf_ia_pho * joint_pho
    d/dt(joint) <- kf_ia_pho * joint_pho * mw_dex / mw_pho - ka_ia * joint

    # Blood, liver and remainder tissue (Eqs. 6-8 as coded in $DES)
    d/dt(a_blood) <- input_iv_im + input_ia +
      qliver * cliver * fut + qremainder * cremainder / kp_remainder -
      qco * Cc - clr * Cc
    d/dt(a_liver) <- qliver * (Cc - cliver * fut) - clint * cliver * fut
    d/dt(a_remainder) <- qremainder * (Cc - cremainder / kp_remainder)

    # Bladder as a well-stirred compartment of constant volume (Eq. 10)
    d/dt(bladder) <- clr * Cc - qurine * bladder / vbladder
    Curine <- bladder / vbladder

    # ---- Cortisol (Eqs. 11-14) ----
    rbase_cortisol <- exp((1 - STUDY_EKSTRAND2015) * lrbase_cortisol +
      STUDY_EKSTRAND2015 * lrbase_cortisol_ekstrand2015 + etalrbase_cortisol)
    amp_cortisol <- exp(lamp_cortisol)
    kout_cortisol <- exp(lkout_cortisol)
    ic50_cortisol <- exp(lic50_cortisol)
    rmax_stress <- exp(lrmax_stress)
    tau_stress <- exp(ltau_stress)
    omega_circ <- 2 * pi / 24
    # Circadian baseline (Eq. 11) and the production rate that sustains it (Eq. 13)
    cortisol_bl <- rbase_cortisol + amp_cortisol * cos(omega_circ * (t - tacro_cortisol))
    kin_cortisol <- kout_cortisol * rbase_cortisol +
      kout_cortisol * amp_cortisol * cos(omega_circ * (t - tacro_cortisol)) -
      omega_circ * amp_cortisol * sin(omega_circ * (t - tacro_cortisol))
    # Stress surge Rmax,stress / tau for tau hours after a dose, in the
    # Toutain 1984 study only ($DES IF (STUDY_ID .EQ. 1 .AND. TAD .LT. Tau));
    # tad() is NA before the first dose, so the surge is off until then.
    tad_dose <- tad()
    stress <- 0
    if (STUDY_TOUTAIN1984 == 1 && tad_dose >= 0 && tad_dose < tau_stress) {
      stress <- rmax_stress / tau_stress
    }
    cc_pos <- max(Cc, 0)
    inh_cortisol <- imax_cortisol * cc_pos^hill_cortisol /
      (cc_pos^hill_cortisol + ic50_cortisol^hill_cortisol)
    # Eq. 14 with R = cortisol_bl + cortisol_dev: substituting Eq. 13 gives
    # d(cortisol_dev)/dt = stress - kin(t) * I(Cp) - kout * cortisol_dev, with
    # cortisol_dev(0) = 0 (R starts on its baseline, $PK A_0(7)).
    d/dt(cortisol_dev) <- stress - kin_cortisol * inh_cortisol - kout_cortisol * cortisol_dev
    cortisol <- cortisol_bl + cortisol_dev

    Cc ~ prop(propSd)
    Curine ~ prop(propSd_Curine)
    cortisol ~ add(addSd_cortisol)
  })
}
