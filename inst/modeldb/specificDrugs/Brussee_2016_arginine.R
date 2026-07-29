Brussee_2016_arginine <- function() {
  description <- paste(
    "Two-compartment population PKPD model for intravenous L-arginine adjunctive",
    "therapy in 73 adults with moderately severe falciparum malaria. Exogenous",
    "L-arginine PK is two-compartment IV with allometric scaling on CL and V1 and",
    "a multiplicative Papuan-ethnicity effect on CL. Endogenous L-arginine",
    "concentration follows a second-order polynomial recovery function indexed",
    "from approximately two days before presentation (start of symptoms);",
    "lognormal between-subject variability multiplies the typical polynomial.",
    "Pharmacodynamic output 1 (exhaled NO, ppb) is linear in the exogenous",
    "arginine concentration with a per-subject baseline. Pharmacodynamic output 2",
    "(reactive-hyperemia peripheral-arterial-tonometry index, RH-PAT) is linear in",
    "the predicted NO minus its baseline, i.e., linear in the exogenous arginine",
    "concentration."
  )
  reference <- paste(
    "Brussee JM, Yeo TW, Lampah DA, Anstey NM, Duffull SB (2016).",
    "Pharmacokinetic-Pharmacodynamic Model for the Effect of L-Arginine on",
    "Endothelial Function in Patients with Moderately Severe Falciparum Malaria.",
    "Antimicrob Agents Chemother 60(1):198-205. doi:10.1128/AAC.01479-15.",
    "PK structure adapted from Yeo TW, Lampah DA, Gitawati R, Tjitra E, Kenangalem E,",
    "Price RN, Anstey NM, Duffull SB (2008). Pharmacokinetics of L-arginine in",
    "adults with moderately severe malaria. Antimicrob Agents Chemother 52(12):4381-4387.",
    "doi:10.1128/AAC.00421-08.",
    sep = " "
  )
  vignette <- "Brussee_2016_arginine"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL (estimated exponent e_wt_cl = 2.47) and V1",
        "(estimated exponent e_wt_vc = 0.757), normalised to the reference",
        "weight of 60 kg. Reference weight is the standardising centring value",
        "used in the paper's equations for CL_i and V1_i (Brussee 2016 Table 2",
        "footnotes a and b, in the Yeo 2008 lineage). Median observed WT in the",
        "73-subject cohort was 57 kg (range 42-77)."
      ),
      source_name        = "WT"
    ),
    RACE_PAPUAN = list(
      description        = "Papuan / indigenous Melanesian heritage indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Papuan; in the Brussee 2016 cohort, mainland-Indonesian transmigrants)",
      notes              = paste(
        "Multiplicative effect on CL: CL_i = CL * (WT/60)^e_wt_cl * f^RACE_PAPUAN",
        "with f = 1.9 (Brussee 2016 Table 2). Source paper's binary covariate is",
        "named 'Papuan' (1 = indigenous Melanesian Papuan, 0 = non-Papuan).",
        "Cohort was 86% Papuan (63 of 73 subjects); the non-Papuan reference",
        "group is predominantly mainland-Indonesian transmigrants enrolled at",
        "Mitra Masyarakat Hospital in Timika, Papua, Indonesia. PK structural",
        "model was developed in Yeo 2008 (52(12):4381) with the same",
        "Papuan-vs-non-Papuan stratification."
      ),
      source_name        = "Papuan"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 73L,
    n_studies      = 2L,
    age_range      = "18-56 years",
    age_median     = "27 years",
    weight_range   = "42-77 kg",
    weight_median  = "57 kg",
    sex_female_pct = 34,
    race_ethnicity = c(Papuan = 86, NonPapuan = 14),
    disease_state  = paste(
      "Moderately severe Plasmodium falciparum malaria (defined by parasitaemia",
      "and the absence of severe-malaria criteria). All subjects enrolled at",
      "Mitra Masyarakat Hospital in Timika, Papua, Indonesia. 30 subjects",
      "received intravenous L-arginine as adjunctive therapy in three dose",
      "groups (3, 6, or 12 g over 0.5 h); 43 subjects received saline placebo",
      "in an observational cohort. Both arms received standard antimalarial",
      "treatment (artesunate-based)."
    ),
    dose_range     = paste(
      "3, 6, or 12 g L-arginine as a 0.5-hour intravenous infusion (n = 10 per",
      "dose group); saline placebo (n = 43) for the observational arm. Single",
      "dose at study start."
    ),
    regions        = "Indonesia (Papua province, Timika)",
    notes          = paste(
      "Pooled cohorts from two prior studies: a natural-history observational",
      "study of 48 subjects (Yeo 2008 ref 16; 43 retained after exclusion for",
      "missing PD data) and a non-randomised intervention trial of 30 subjects",
      "(Yeo 2008 PK paper ref 17, AAC 52(12):4381). Patient characteristics in",
      "Brussee 2016 Table 1; cohort demographics are not significantly different",
      "between the L-arginine and saline arms. Sex distribution from Table 1",
      "(48 of 73 = 66% male, therefore 34% female)."
    )
  )

  ini({
    # =====================================================================
    # PK MODEL (exogenous L-arginine)
    # Two-compartment IV model. Reference values are Brussee 2016 Table 2
    # point estimates. Allometric exponents on CL (e_wt_cl, paper's k1) and V1
    # (e_wt_vc, paper's k2) were estimated, not fixed (no FIX flag in Table 2).
    # Papuan-vs-non-Papuan multiplicative factor f on CL was estimated.
    # =====================================================================
    lcl     <- log(30.3)     ; label("L-arginine clearance (L/h)")                                        # Brussee 2016 Table 2 (CL = 30.3, RSE 15.4%)
    lvc     <- log(26.6)     ; label("L-arginine central volume of distribution (L)")                     # Brussee 2016 Table 2 (V1 = 26.6, RSE 12.0%)
    lvp     <- log(80.1)     ; label("L-arginine peripheral volume of distribution (L)")                  # Brussee 2016 Table 2 (V2 = 80.1, RSE 27.5%)
    lq      <- log(22)       ; label("L-arginine inter-compartmental clearance (L/h)")                    # Brussee 2016 Table 2 (Q = 22, RSE 18.6%)

    # Allometric WT scaling (estimated, not fixed at 0.75 / 1)
    e_wt_cl <- 2.47          ; label("Allometric exponent of body weight on CL (unitless)")               # Brussee 2016 Table 2 (k1 = 2.47, RSE 21.9%)
    e_wt_vc <- 0.757         ; label("Allometric exponent of body weight on V1 (unitless)")               # Brussee 2016 Table 2 (k2 = 0.757, RSE 70.4%)

    # Papuan ethnicity effect on CL: multiplicative factor applied as
    # e_papuan_cl^RACE_PAPUAN (same form as Goel 2016 sonidegib, Wang 2017
    # benralizumab, Lu 2019 polatuzumab Asian effects). Paper's covariate is
    # named 'Papuan' (1 = Papuan, 0 = non-Papuan).
    e_papuan_cl <- 1.9       ; label("Multiplicative Papuan-vs-non-Papuan CL ratio (unitless)")           # Brussee 2016 Table 2 (f = 1.9, RSE 16.9%)

    # =====================================================================
    # ENDOGENOUS L-ARGININE BASELINE (paper Eq. 1)
    # Second-order polynomial in time indexed to approximately 2 days prior to
    # presentation (paper Methods + Eq. 1 narrative). t = 0 in the polynomial
    # is the start of symptoms; in the study clock that is approximately
    # t_study = -48 h, so during the 48 h observation window the polynomial is
    # evaluated at t_paper = study_time + 48 h. Per-subject IIV multiplies the
    # entire typical polynomial via exp(eta), so a single eta (etalrbase_arg)
    # scales the whole baseline trajectory.
    # =====================================================================
    lrbase_arg <- log(6.07)  ; label("Typical baseline arginine at symptom onset (mg/L)")                 # Brussee 2016 Table 2 (Arg_0 = 6.07 mg/L, RSE 5.7%)
    theta_t1   <- 0.0365     ; label("First-order recovery coefficient of arginine baseline (mg/L/h)")    # Brussee 2016 Table 2 (theta_t1 = 0.0365, RSE 50.4%)
    theta_t2   <- -0.0000348 ; label("Second-order recovery coefficient of arginine baseline (mg/L/h^2)") # Brussee 2016 Table 2 (theta_t2 = -3.48e-5, RSE 53.4%)

    # =====================================================================
    # PD MODEL: NITRIC OXIDE (paper Eq. 2)
    # Linear model: NO_pred = BL_NO + SL_NO * (E[Arg] - BL_Arg).
    # E[Arg] - BL_Arg = central/V1, i.e., only the exogenous arginine
    # contribution drives NO. The endogenous polynomial cancels because both
    # E[Arg] and BL_Arg evolve with it identically.
    # =====================================================================
    lrbase_no <- log(18.2)   ; label("Typical exhaled NO baseline (ppb)")                                 # Brussee 2016 Table 3 (BL_NO = 18.2 ppb, RSE 7.5%)
    lsl_no    <- log(0.0243) ; label("Log-slope of NO on exogenous arginine concentration (ppb per mg/L)")# Brussee 2016 Table 3 (SL_NO = 0.0243, RSE 56.8%)

    # =====================================================================
    # PD MODEL: ENDOTHELIAL FUNCTION / RH-PAT INDEX (paper Eq. 3)
    # Linear model: EF_pred = BL_EF + SL_EF * (NO_pred - BL_NO), where
    # NO_pred is evaluated at epsilon = 0 (typical NO trajectory).
    # Substituting gives EF_pred = BL_EF + SL_EF * SL_NO * central/V1.
    # =====================================================================
    lrbase_ef <- log(1.86)   ; label("Typical RH-PAT (endothelial-function) baseline (unitless)")         # Brussee 2016 Table 3 (BL_EF = 1.86, RSE 1.9%)
    lsl_ef    <- log(0.145)  ; label("Log-slope of RH-PAT on (NO - NO_baseline) (per ppb)")               # Brussee 2016 Table 3 (SL_EF = 0.145, RSE 73.1%)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY (paper Tables 2 and 3 omega %CV)
    # Reported as % coefficient of variation; converted to log-scale variance
    # via omega^2 = log(1 + CV^2).
    # =====================================================================
    etalcl        ~ 0.10878   # Brussee 2016 Table 2 (omega_CL = 33.9% CV; omega^2 = log(1 + 0.339^2) = 0.10878)
    etalvc        ~ 0.03732   # Brussee 2016 Table 2 (omega_V1 = 19.5% CV; omega^2 = log(1 + 0.195^2) = 0.03732)
    etalrbase_arg ~ 0.18067   # Brussee 2016 Table 2 (omega_Arg_0 = 44.5% CV; omega^2 = log(1 + 0.445^2) = 0.18067)
    etalrbase_no  ~ 0.29957   # Brussee 2016 Table 3 (omega_BL_NO = 59.1% CV; omega^2 = log(1 + 0.591^2) = 0.29957)
    etalrbase_ef  ~ 0.01550   # Brussee 2016 Table 3 (omega_BL_EF = 12.5% CV; omega^2 = log(1 + 0.125^2) = 0.01550)

    # =====================================================================
    # RESIDUAL ERROR (Tables 2 and 3)
    # Paper reports a combined additive + proportional model with sigma_add
    # fixed at approximately zero ("~= 0 FIX") for all three outputs and the
    # proportional term as % CV. Encoded as proportional-only since the
    # additive term carries no estimated information in the source.
    # =====================================================================
    propSd    <- 0.286       ; label("Arginine proportional residual error (fraction)")                   # Brussee 2016 Table 2 (sigma_prop = 28.6% CV, RSE 20.5%); sigma_add ~= 0 FIX -> dropped
    propSd_NO <- 0.355       ; label("NO proportional residual error (fraction)")                         # Brussee 2016 Table 3 (sigma_prop = 35.5% CV, RSE 22.7%); sigma_add ~= 0 FIX -> dropped
    propSd_EF <- 0.213       ; label("RH-PAT proportional residual error (fraction)")                     # Brussee 2016 Table 3 (sigma_prop = 21.3% CV, RSE 10.0%); sigma_add ~= 0 FIX -> dropped
  })
  model({
    # ===== Individual PK parameters =====
    # WT scaling normalised to the paper's 60 kg reference. Papuan effect is a
    # multiplicative factor f = 1.9, applied to CL via the same f^RACE_PAPUAN
    # form used in Goel 2016 sonidegib (e_race_japanese_cl^RACE_JAPANESE).
    cl <- exp(lcl + etalcl) * (WT / 60)^e_wt_cl * e_papuan_cl^RACE_PAPUAN
    vc <- exp(lvc + etalvc) * (WT / 60)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ===== Endogenous baseline polynomial =====
    # t = 0 in the paper polynomial is the start of symptoms; the study clock
    # starts approximately 48 h later (paper Eq. 1 narrative). At study time
    # 0 the polynomial is evaluated at t_paper = 48 h; over the 48 h
    # observation window t_paper goes from 48 to ~96 h.
    t_sym      <- time + 48
    arg_0      <- exp(lrbase_arg)
    bl_arg_typ <- arg_0 + theta_t1 * t_sym + theta_t2 * t_sym * t_sym
    bl_arg     <- bl_arg_typ * exp(etalrbase_arg)

    # ===== PD baselines (per-subject) =====
    bl_no <- exp(lrbase_no + etalrbase_no)
    bl_ef <- exp(lrbase_ef + etalrbase_ef)
    sl_no <- exp(lsl_no)
    sl_ef <- exp(lsl_ef)

    # ===== ODE system =====
    # central holds exogenous arginine mass (mg); peripheral1 likewise. IV
    # doses (3, 6, or 12 g over 0.5 h) target central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ===== Observations =====
    # Cc: total measured arginine = exogenous PK contribution (central/V1) +
    # endogenous polynomial baseline. This is what the HPLC assay measures.
    Cc <- central / vc + bl_arg

    # NO: linear in exogenous arginine (BL_Arg cancels in NO_pred - BL_NO).
    NO <- bl_no + sl_no * (central / vc)

    # EF (RH-PAT index): linear in NO_pred - BL_NO_i = sl_no * central/vc.
    EF <- bl_ef + sl_ef * sl_no * (central / vc)

    Cc ~ prop(propSd)
    NO ~ prop(propSd_NO)
    EF ~ prop(propSd_EF)
  })
}
