Yamada_2025_zolbetuximab <- function() {
  description <- "Two-compartment population PK model of zolbetuximab (anti-CLDN18.2 IgG1 mAb) with zero-order IV input and time-dependent clearance in patients with locally advanced unresectable or metastatic gastric/gastroesophageal junction (G/GEJ) adenocarcinoma (Yamada 2025)"
  reference <- "Yamada A, Takeuchi M, Komatsu K, Bonate PL, Poondru S, Yang J. Population PK and Exposure-Response Analyses of Zolbetuximab in Patients With Locally Advanced Unresectable or Metastatic G/GEJ Adenocarcinoma. Clinical and Translational Science. 2025;18(7):e70280. doi:10.1111/cts.70280"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CLss, CLT, Q (exponent 1.06) and V1, V2 (exponent 0.968); reference 1.70 m^2 per the Yamada 2025 Figure 1 reference population. BSA computation method (DuBois / Mosteller / Haycock) is not specified in the paper.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CLss (exponent -0.535) and Kdecay (exponent 1.48); reference 39.1 g/L per the Yamada 2025 Figure 1 reference population.",
      source_name        = "ALB"
    ),
    GAST = list(
      description        = "Prior gastrectomy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior gastrectomy)",
      notes              = "Time-fixed per subject. Fractional-change (dummy-variable) effects on CLss (-0.182), CLT (-0.495), and V1 (+0.103) in the final time-dependent-clearance model (Yamada 2025 Table 1).",
      source_name        = "GAST"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source paper uses the column label SEX (1 = female per the 'SEX on CLss (if female)' row of Table 1); renamed to the canonical SEXF per covariate-columns.md. Fractional-change effects on CLss (-0.195) and V1 (-0.108).",
      source_name        = "SEX"
    ),
    HGB = list(
      description        = "Baseline hemoglobin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on V1 (exponent -0.374); reference 118 g/L per the Yamada 2025 Figure 1 reference population.",
      source_name        = "HGB"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on V1 (exponent 0.0347); reference 0.38 mg/dL per the Yamada 2025 Figure 1 reference population.",
      source_name        = "TBILI"
    ),
    COMB_EOX = list(
      description        = "Concomitant EOX (epirubicin + oxaliplatin + capecitabine) chemotherapy backbone indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-EOX backbone; includes mFOLFOX6, CAPOX, and single-agent settings)",
      notes              = "Source paper uses 'COMB' as the categorical column with EOX as the non-reference level; renamed to the canonical COMB_EOX per covariate-columns.md to preserve the semantic meaning of the 1-level. Fractional-change effect on V1 (+0.466) in the time-dependent-clearance model (Yamada 2025 Table 1).",
      source_name        = "COMB"
    )
  )

  population <- list(
    n_subjects         = 714L,
    n_studies          = 8L,
    phase_mix          = "3 phase 1, 3 phase 2 (MONO, FAST, ILUSTRO), and 2 phase 3 studies (SPOTLIGHT, GLOW)",
    age_range          = "Adult (>= 18 years); full range not reported in the main paper (see Table S2).",
    weight_range       = "Not reported in the main paper (see Table S2); reference BSA 1.70 m^2 approximates a 70 kg adult with average height.",
    sex_female_pct     = NA_real_,
    race_ethnicity     = "Multi-regional enrollment: White, Asian (including Chinese, Japanese, Korean subgroups), and others. Race/ethnicity had no clinically important effect on zolbetuximab PK in the final model.",
    disease_state      = "Locally advanced unresectable or metastatic gastric (n = 540) or gastroesophageal junction (n = 174) adenocarcinoma. Prior gastrectomy status included as a covariate.",
    dose_range         = "33-1000 mg/m^2 IV infusion. Clinical regimen: 800 mg/m^2 loading dose followed by 600 mg/m^2 every 3 weeks (phase 3), also evaluated as 800/400 mg/m^2 every 2 weeks.",
    regions            = "Global / multi-regional.",
    n_observations     = 5066L,
    reference_subject  = "Male, BSA 1.70 m^2, ALB 39.1 g/L, HGB 118 g/L, TBILI 0.38 mg/dL, no prior gastrectomy, non-EOX chemotherapy backbone (per Yamada 2025 Figure 1 caption).",
    notes              = "Immunogenicity was not retained as a covariate because of low ADA incidence; race and mild/moderate renal impairment and mild hepatic impairment also showed no clinically important effect. Severe renal impairment and moderate/severe hepatic impairment were under-represented and not evaluated. The final model selected in the paper (Table 1, footnote a) is the time-dependent-clearance (TDC) variant."
  )

  ini({
    # Structural parameters (for the reference subject). Paper reports clearances
    # in L/h; we use time in days here to align with Kdecay (1/day), so L/h is
    # multiplied by 24 to give L/day. Values are from Yamada 2025 Table 1
    # (final TDC model; footnote a identifies it as the final model).
    lclss   <- log(0.0117 * 24); label("Steady-state clearance component (CLss, L/day)")    # Yamada 2025 Table 1 (0.0117 L/h * 24)
    lclt    <- log(0.0159 * 24); label("Time-decaying clearance component (CLT, L/day)")    # Yamada 2025 Table 1 (0.0159 L/h * 24)
    lkdecay <- log(0.0209);      label("First-order decay rate of CLT (Kdecay, 1/day)")     # Yamada 2025 Table 1
    lvc     <- log(3.04);        label("Central volume of distribution (V1, L)")            # Yamada 2025 Table 1
    lq      <- log(0.0235 * 24); label("Intercompartmental clearance (Q, L/day)")           # Yamada 2025 Table 1 (0.0235 L/h * 24)
    lvp     <- log(2.49);        label("Peripheral volume of distribution (V2, L)")         # Yamada 2025 Table 1

    # Covariate-effect parameters (Yamada 2025 Table 1, TDC model / final model).
    # Continuous covariates enter as power models, normalized to the Figure 1
    # reference values (BSA 1.70 m^2, ALB 39.1 g/L, HGB 118 g/L, TBILI 0.38 mg/dL).
    # Categorical covariates enter via the NONMEM dummy-variable form
    # (param = typical * (1 + theta * I)).
    e_bsa_cl     <-  1.06;   label("Power exponent of BSA on CLss, CLT, Q (unitless)")               # Yamada 2025 Table 1
    e_bsa_v      <-  0.968;  label("Power exponent of BSA on V1, V2 (unitless)")                     # Yamada 2025 Table 1
    e_alb_clss   <- -0.535;  label("Power exponent of albumin on CLss (unitless)")                   # Yamada 2025 Table 1
    e_alb_kdecay <-  1.48;   label("Power exponent of albumin on Kdecay (unitless)")                 # Yamada 2025 Table 1
    e_hgb_vc     <- -0.374;  label("Power exponent of hemoglobin on V1 (unitless)")                  # Yamada 2025 Table 1
    e_tbili_vc   <-  0.0347; label("Power exponent of total bilirubin on V1 (unitless)")             # Yamada 2025 Table 1
    e_gast_clss  <- -0.182;  label("Fractional change in CLss for prior gastrectomy (unitless)")      # Yamada 2025 Table 1
    e_gast_clt   <- -0.495;  label("Fractional change in CLT for prior gastrectomy (unitless)")       # Yamada 2025 Table 1
    e_gast_vc    <-  0.103;  label("Fractional change in V1 for prior gastrectomy (unitless)")        # Yamada 2025 Table 1
    e_sex_clss   <- -0.195;  label("Fractional change in CLss for females (unitless)")                # Yamada 2025 Table 1
    e_sex_vc     <- -0.108;  label("Fractional change in V1 for females (unitless)")                  # Yamada 2025 Table 1
    e_eox_vc     <-  0.466;  label("Fractional change in V1 for EOX chemotherapy backbone (unitless)")# Yamada 2025 Table 1

    # Inter-individual variability. The paper reports %CV on log-normal
    # parameters; the stored variance follows omega^2 = log(CV^2 + 1).
    etalclss   ~ 0.0669  # 26.3% CV; Yamada 2025 Table 1
    etalclt    ~ 0.4569  # 76.1% CV; Yamada 2025 Table 1
    etalkdecay ~ 0.4685  # 77.3% CV; Yamada 2025 Table 1
    etalvc     ~ 0.0396  # 20.1% CV; Yamada 2025 Table 1
    etalq      ~ 0.3424  # 63.9% CV; Yamada 2025 Table 1
    etalvp     ~ 0.0724  # 27.4% CV; Yamada 2025 Table 1

    # Residual error. Yamada 2025 Table 1 reports the NONMEM $SIGMA estimate
    # 0.169 for the proportional error and 4.03 ug/mL for the additive error
    # (final TDC model). Following the Thakre 2022 convention in this repo,
    # the proportional value is interpreted as the NONMEM variance so
    # propSd = sqrt(variance); the additive value is interpreted as an SD in
    # ug/mL (consistent with the column-header unit).
    propSd <- sqrt(0.169); label("Proportional residual error (SD, fraction)")   # Yamada 2025 Table 1
    addSd  <- 4.03;        label("Additive residual error (ug/mL)")              # Yamada 2025 Table 1
  })
  model({
    # Individual PK parameters. Reference subject (Yamada 2025 Figure 1):
    # BSA = 1.70 m^2, ALB = 39.1 g/L, HGB = 118 g/L, TBILI = 0.38 mg/dL, male,
    # no prior gastrectomy, non-EOX chemotherapy backbone.
    clss <- exp(lclss + etalclss) *
      (BSA / 1.70)^e_bsa_cl *
      (ALB / 39.1)^e_alb_clss *
      (1 + e_gast_clss * GAST) *
      (1 + e_sex_clss * SEXF)

    clt <- exp(lclt + etalclt) *
      (BSA / 1.70)^e_bsa_cl *
      (1 + e_gast_clt * GAST)

    kdecay <- exp(lkdecay + etalkdecay) *
      (ALB / 39.1)^e_alb_kdecay

    vc <- exp(lvc + etalvc) *
      (BSA / 1.70)^e_bsa_v *
      (HGB / 118)^e_hgb_vc *
      (TBILI / 0.38)^e_tbili_vc *
      (1 + e_gast_vc * GAST) *
      (1 + e_sex_vc * SEXF) *
      (1 + e_eox_vc * COMB_EOX)

    q <- exp(lq + etalq) *
      (BSA / 1.70)^e_bsa_cl

    vp <- exp(lvp + etalvp) *
      (BSA / 1.70)^e_bsa_v

    # Time-dependent total clearance (Yamada 2025 Equation 1):
    #   CL(t) = CLss + CLT * exp(-Kdecay * t)
    # `time` is the internal integration time in days, which corresponds to
    # time from the first dose for the event datasets this model expects.
    cl <- clss + clt * exp(-kdecay * time)

    # Two-compartment model with zero-order IV input (infusion rate supplied
    # via the `rate` column on dose events).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
