Wu_2012_modafinil <- function() {
  description <- "Joint parent + metabolite population pharmacokinetic model for oral modafinil and its principal carboxylic-acid metabolite modafinil acid (2-[(diphenylmethyl)sulfonyl]acetic acid) in 49 healthy volunteers from five major ethnic groups of China (Han, Mongolian, Korean, Uygur, Hui) under a single 200 mg oral dose (Wu 2012). Four-compartment NONMEM ADVAN6 / L2 structure: GI depot with first-order absorption (ka), two-compartment modafinil disposition (apparent CL/F, Vc/F, Q/F, Vp/F), and a one-compartment modafinil acid disposition (apparent CL3/F1F2, V3/F1F2). All modafinil elimination is treated as forming modafinil acid at the apparent-parameter level because F2 (the absolute modafinil-to-acid metabolic-conversion fraction) is not identifiable from oral plasma data alone; F2 is absorbed into the apparent acid parameters. Sex acts on CL/F, Q/F, and Vp/F of modafinil; ethnicity acts on Vc/F of modafinil (Korean and Hui share a single composite multiplier; Mongolian and Uygur each have their own) and on CL3/F1F2 of modafinil acid (Han and Mongolian share the reference; Korean has its own multiplier; Uygur and Hui share a composite multiplier)."
  reference <- paste(
    "Wu KH, Guo T, Deng CH, Guan Z, Li L, Zhou TY, Lu W.",
    "Population pharmacokinetics of modafinil acid and estimation of",
    "the metabolic conversion of modafinil into modafinil acid in 5",
    "major ethnic groups of China. Acta Pharmacol Sin.",
    "2012 Nov;33(11):1401-1408. doi:10.1038/aps.2012.124.",
    sep = " "
  )
  vignette <- "Wu_2012_modafinil"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Acts multiplicatively on CL/F, Q/F, and Vp/F of modafinil. Source paper reports sex-specific typical values in Wu 2012 Table 2 ('Male' / 'Female' rows for CL1/F1, CL2/F1, V2/F1); the canonical SEXF encoding is sex = female = 1, so the log-ratio effect e_sexf_<param> equals log(female-typical / male-typical).",
      source_name        = "Sex (male/female)"
    ),
    RACE_CN_MONGOLIAN = list(
      description        = "PRC-internal Mongolian ethnic-minority indicator, 1 = Mongolian (one of the 56 officially recognised PRC ethnic groups), 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Han is the typical-value reference for V1; Han + Mongolian are the joint typical-value reference for CL3, encoded by setting both RACE_CN_MONGOLIAN = 0 and the four other indicators = 0 = Han)",
      notes              = "Wu 2012 Table 1 'Ethnicity (male/female) [classification]' column codes Mongolian as 2. Acts multiplicatively on Vc/F of modafinil (factor 1.65 vs Han reference; Wu 2012 Table 2 'theta_COV-V1: Mongolian = 1.65'). Does NOT affect CL3/F1F2 of modafinil acid (Mongolian shares the Han reference for the acid clearance per Wu 2012 Table 2 footnote 'theta_COV-CL3 is 1 for Han and Mongolian groups').",
      source_name        = "Ethnicity (integer 2 = Mongolian)"
    ),
    RACE_CN_KOREAN = list(
      description        = "PRC-internal Korean ethnic-minority indicator, 1 = Chinese-Korean (one of the 56 officially recognised PRC ethnic groups; primarily resident in northeast China), 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Han is the typical-value reference for V1; Han + Mongolian for CL3)",
      notes              = "Wu 2012 Table 1 'Ethnicity (male/female) [classification]' column codes Korean as 3. Korean shares the Vc/F composite multiplier (0.86 vs Han) with Hui per Wu 2012 Table 2 ('theta_COV-V1: Korean or Hui = 0.86'); the model() block forms the composite via OR-logic (RACE_CN_KOREAN + RACE_CN_HUI). Korean has its own CL3/F1F2 multiplier (1.25 vs Han / Mongolian) per Wu 2012 Table 2 ('theta_COV-CL3: Korean = 1.25').",
      source_name        = "Ethnicity (integer 3 = Korean)"
    ),
    RACE_CN_UYGUR = list(
      description        = "PRC-internal Uygur (Uyghur) ethnic-minority indicator, 1 = Uygur (a Turkic Muslim ethnic group primarily resident in the Xinjiang Uygur Autonomous Region; one of the 56 officially recognised PRC ethnic groups), 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Han is the typical-value reference for V1; Han + Mongolian for CL3)",
      notes              = "Wu 2012 Table 1 'Ethnicity (male/female) [classification]' column codes Uygur as 4. Uygur has its own Vc/F multiplier (1.33 vs Han) per Wu 2012 Table 2 ('theta_COV-V1: Uygur = 1.33'). Uygur shares the CL3/F1F2 composite multiplier (1.15 vs Han / Mongolian) with Hui per Wu 2012 Table 2 ('theta_COV-CL3: Uygur or Hui = 1.15'); the model() block forms the composite via OR-logic (RACE_CN_UYGUR + RACE_CN_HUI).",
      source_name        = "Ethnicity (integer 4 = Uygur)"
    ),
    RACE_CN_HUI = list(
      description        = "PRC-internal Hui ethnic-minority indicator, 1 = Hui (a Sinophone Muslim ethnic group widely distributed across northwestern China with the largest concentration in the Ningxia Hui Autonomous Region; one of the 56 officially recognised PRC ethnic groups), 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Han is the typical-value reference for V1; Han + Mongolian for CL3)",
      notes              = "Wu 2012 Table 1 'Ethnicity (male/female) [classification]' column codes Hui as 5. Hui shares the Vc/F composite multiplier (0.86 vs Han) with Korean per Wu 2012 Table 2 ('theta_COV-V1: Korean or Hui = 0.86') and the CL3/F1F2 composite multiplier (1.15 vs Han / Mongolian) with Uygur per Wu 2012 Table 2 ('theta_COV-CL3: Uygur or Hui = 1.15'); both composites are formed via OR-logic in model().",
      source_name        = "Ethnicity (integer 5 = Hui)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49L,
    n_studies      = 1L,
    age_range      = "18-26 years (mean 22.4, SD 1.7)",
    age_median     = "approximately 22 years (mean)",
    weight_range   = "44-88 kg (mean 60.9, SD 10.6)",
    weight_median  = "60.9 kg (mean)",
    height_range   = "150-184 cm (mean 166.6, SD 8.3)",
    bmi_range      = "16.9-31.6 kg/m^2 (mean 21.9, SD 3.1)",
    sex_female_pct = 51.0,
    race_ethnicity = c(Han = 20.4, Mongolian = 20.4, Korean = 18.4, Uygur = 20.4, Hui = 20.4),
    disease_state  = "Healthy young volunteers; hepatic function and routine blood / biochemical parameters within normal range; no medications for at least 72 h before dosing; alcohol and smoking forbidden for at least 72 h before dosing and during sampling; females studied during luteal phase of the menstrual cycle.",
    dose_range     = "Single 200 mg oral modafinil dose (two 100 mg tablets, Jiangzhong Pharmaceutical Co Ltd, China) with 200 mL of water, after at least 8 h fasting.",
    regions        = "China (single-centre study, Shenyang Northern Hospital).",
    notes          = "Ethnicity counts (Table 1): Han 10 (5 male / 5 female), Mongolian 10 (5/5), Korean 9 (4/5), Uygur 10 (5/5), Hui 10 (5/5); minority families were single-ethnicity for three generations. Demographics from Wu 2012 Table 1; PK parameter estimates from Wu 2012 Table 2 'Final model'. The 49 subjects contributed 637 plasma concentration observations for each analyte. PopPK fitted with NONMEM VII (ICON) FOCE-I with interaction, ADVAN6 differential-equation system, L2 option for paired modafinil and acid observations."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification of each item):
  # * Apparent-parameter framework. Wu 2012 estimates CL/F, Vc/F, Q/F, Vp/F
  #   for modafinil and CL3/(F1*F2), V3/(F1*F2) for modafinil acid because
  #   bioavailability F1 (modafinil) and metabolic-conversion fraction F2
  #   (modafinil -> modafinil acid) are not separately identifiable from
  #   oral plasma data without a parenteral or urinary reference arm. The
  #   simulator runs in apparent units, treating all of modafinil's
  #   clearance flux (CL/Vc * central) as forming modafinil acid; F2 is
  #   absorbed into the apparent V3 / CL3, so the simulated acid
  #   concentration reproduces the observed plasma modafinil acid time
  #   course in micromolar units.
  # * Sex effect. Wu 2012 Table 2 reports separate male / female typical
  #   values for CL1/F1, CL2/F1, V2/F1 of modafinil ('Final model'
  #   columns). The model encodes the male values as the typical lcl,
  #   lq, lvp and applies a log-ratio shift e_sexf_<param> = log(female /
  #   male) on SEXF = 1 to recover the female-specific value.
  # * Ethnicity composite groupings. Wu 2012 estimates separate
  #   coefficients per ethnic-minority group on V1 and on CL3 of acid.
  #   Wu 2012 Table 2 ('Final model' column) shows the V1 grouping
  #   {Han, Uygur, Mongolian, (Korean or Hui)} and the CL3 grouping
  #   {(Han or Mongolian), Korean, (Uygur or Hui)}. The Korean+Hui V1
  #   composite and the Uygur+Hui CL3 composite are formed inside the
  #   model() block via OR-logic across the subject-level binary
  #   indicators; the data preserves the five-way ethnicity dimension so
  #   either composite could be reconstructed at simulation time.
  # * Residual error. Wu 2012 Table 2 reports a combined proportional +
  #   additive error model for modafinil (epsilon1 proportional = 13.4 %,
  #   epsilon2 additive = 0.001 umol/L) and a proportional-only error
  #   model for modafinil acid (epsilon3 proportional = 8.94 %).
  # * Concentration units. Wu 2012 Table 2 reports the additive residual
  #   error magnitude as 0.001 umol/L, so plasma concentrations were fit
  #   in micromolar units. The dose record is in mg of modafinil
  #   (modafinil molecular weight 273.35 g/mol); the model converts at
  #   the observation step via 1000 / 273.35 to produce micromolar plasma
  #   concentrations for both the parent and the acid, the latter
  #   transferred mole-for-mole through the apparent metabolite
  #   formation flux.
  ini({
    # ----- Structural PK -- modafinil parent (Wu 2012 Table 2 'Final model') -----
    # Typical values reported as the male / Han reference (see Implementation notes).
    lka  <- log(0.755); label("First-order absorption rate constant ka (1/h)")                                    # Wu 2012 Table 2 'Final model' ka = 0.755 (RSE 11.1 %)
    lcl  <- log(3.51);  label("Apparent modafinil clearance CL1/F1, male typical (L/h)")                          # Wu 2012 Table 2 'Final model' CL1/F1 (Male) = 3.51 (RSE 7.18 %)
    lvc  <- log(7.74);  label("Apparent modafinil central volume V1/F1, Han typical (L)")                         # Wu 2012 Table 2 'Final model' V1/F1 theta_2 = 7.74 (RSE 5.18 %); Han is the typical-value reference (theta_COV-V1 = 1 for Han)
    lq   <- log(7.02);  label("Apparent modafinil inter-compartmental clearance CL2/F1, male typical (L/h)")      # Wu 2012 Table 2 'Final model' CL2/F1 (Male) = 7.02 (RSE 9.82 %)
    lvp  <- log(35.0);  label("Apparent modafinil peripheral volume V2/F1, male typical (L)")                     # Wu 2012 Table 2 'Final model' V2/F1 (Male) = 35.0 (RSE 8.34 %)

    # ----- Structural PK -- modafinil acid metabolite (Wu 2012 Table 2 'Final model') -----
    # Typical values reported as the Han / Mongolian reference; F1*F2 absorbed into apparent CL3, V3.
    lcl_mfa <- log(4.94); label("Apparent modafinil acid clearance CL3/(F1*F2), Han / Mongolian typical (L/h)")    # Wu 2012 Table 2 'Final model' CL3/(F1*F2) theta_5 = 4.94 (RSE 7.65 %); theta_COV-CL3 = 1 for Han or Mongolian
    lvc_mfa <- log(2.73); label("Apparent modafinil acid central volume V3/(F1*F2) (L)")                           # Wu 2012 Table 2 'Final model' V3/(F1*F2) = 2.73 (RSE 10.5 %); no ethnicity covariate

    # ----- Sex covariate effects (log of female / male ratio, applied on SEXF = 1) -----
    # Females exhibit lower apparent CL/F, Q/F, Vp/F of modafinil per Wu 2012 Table 2 'Final model'.
    e_sexf_cl <- -0.111380; label("Sex (female-vs-male) log-ratio effect on apparent modafinil clearance CL1/F1 (unitless)") # Wu 2012 Table 2 'Final model' CL1/F1 (Female) = 3.14, (Male) = 3.51; log(3.14/3.51) = -0.111380
    e_sexf_q  <- -0.371726; label("Sex (female-vs-male) log-ratio effect on apparent modafinil inter-compartmental clearance CL2/F1 (unitless)") # Wu 2012 Table 2 'Final model' CL2/F1 (Female) = 4.84, (Male) = 7.02; log(4.84/7.02) = -0.371726
    e_sexf_vp <- -0.539189; label("Sex (female-vs-male) log-ratio effect on apparent modafinil peripheral volume V2/F1 (unitless)") # Wu 2012 Table 2 'Final model' V2/F1 (Female) = 20.41, (Male) = 35.0; log(20.41/35.0) = -0.539189

    # ----- Ethnicity covariate effects on apparent modafinil central volume Vc/F (log of multiplier vs Han) -----
    e_race_cn_uygur_vc     <- 0.285179; label("Ethnicity (Uygur-vs-Han) log-multiplier on apparent modafinil central volume V1/F1 (unitless)")          # Wu 2012 Table 2 'Final model' theta_COV-V1: Uygur = 1.33; log(1.33) = 0.285179
    e_race_cn_mongolian_vc <- 0.500775; label("Ethnicity (Mongolian-vs-Han) log-multiplier on apparent modafinil central volume V1/F1 (unitless)")      # Wu 2012 Table 2 'Final model' theta_COV-V1: Mongolian = 1.65; log(1.65) = 0.500775
    e_race_cn_korhui_vc    <- -0.150823; label("Ethnicity (Korean or Hui composite, vs Han) log-multiplier on apparent modafinil central volume V1/F1 (unitless)") # Wu 2012 Table 2 'Final model' theta_COV-V1: Korean or Hui = 0.86; log(0.86) = -0.150823

    # ----- Ethnicity covariate effects on apparent modafinil acid clearance CL3/F1F2 (log of multiplier vs Han / Mongolian) -----
    e_race_cn_korean_cl_mfa <- 0.223144; label("Ethnicity (Korean-vs-Han / Mongolian) log-multiplier on apparent modafinil acid clearance CL3/(F1*F2) (unitless)")              # Wu 2012 Table 2 'Final model' theta_COV-CL3: Korean = 1.25; log(1.25) = 0.223144
    e_race_cn_uyghui_cl_mfa <- 0.139762; label("Ethnicity (Uygur or Hui composite, vs Han / Mongolian) log-multiplier on apparent modafinil acid clearance CL3/(F1*F2) (unitless)") # Wu 2012 Table 2 'Final model' theta_COV-CL3: Uygur or Hui = 1.15; log(1.15) = 0.139762

    # ----- Inter-individual variability (log-normal; omega^2 = log(CV^2 + 1)) -----
    # CV% values from Wu 2012 Table 2 'Final model -- Inter-individual variability [shrinkage]'.
    etalcl     ~ 0.052424   # CL1/F1 23.2 % CV; log(1 + 0.232^2) = 0.052424
    etalvc     ~ 0.233618   # V1/F1 51.3 % CV; log(1 + 0.513^2) = 0.233618
    etalq      ~ 0.047683   # CL2/F1 22.1 % CV; log(1 + 0.221^2) = 0.047683
    etalvp     ~ 0.034733   # V2/F1 18.8 % CV; log(1 + 0.188^2) = 0.034733
    etalcl_mfa ~ 0.033653   # CL3/(F1*F2) 18.5 % CV; log(1 + 0.185^2) = 0.033653
    etalvc_mfa ~ 0.477274   # V3/(F1*F2) 78.2 % CV; log(1 + 0.782^2) = 0.477274
    # Note: Wu 2012 reports no IIV on ka (the Table 2 IIV block lists only CL1, V1, CL2, V2, CL3, V3).

    # ----- Residual error (Wu 2012 Table 2 'Final model -- Intra-individual variability') -----
    propSd     <- 0.134;  label("Proportional residual error on modafinil plasma concentration (fraction)")               # Wu 2012 Table 2 epsilon1 (Modafinil proportional) = 13.4 % [shrinkage 9.14 %]
    addSd      <- 0.001;  label("Additive residual error on modafinil plasma concentration (umol/L)")                     # Wu 2012 Table 2 epsilon2 (Modafinil additive) = 0.001 umol/L [shrinkage 9.14 %]
    propSd_mfa <- 0.0894; label("Proportional residual error on modafinil acid plasma concentration (fraction)")          # Wu 2012 Table 2 epsilon3 (Modafinil acid proportional) = 8.94 % [shrinkage 3.05 %]
  })
  model({
    # Modafinil molecular weight (g/mol) for the mg dose -> umol/L observation
    # conversion; modafinil acid is formed mole-for-mole from modafinil so the
    # acid central compartment accumulates in 'mg of modafinil equivalent' and
    # the same conversion factor applies at the acid observation step.
    mw_modafinil <- 273.35

    # ----- Ethnicity composite multipliers on V1 and CL3 -----
    # Each subject has at most one of the four binary indicators set;
    # the Korean+Hui sum and Uygur+Hui sum each take a value of 0 or 1.
    fethn_vc <- exp(
        e_race_cn_uygur_vc     * RACE_CN_UYGUR
      + e_race_cn_mongolian_vc * RACE_CN_MONGOLIAN
      + e_race_cn_korhui_vc    * (RACE_CN_KOREAN + RACE_CN_HUI)
    )
    fethn_cl_mfa <- exp(
        e_race_cn_korean_cl_mfa * RACE_CN_KOREAN
      + e_race_cn_uyghui_cl_mfa * (RACE_CN_UYGUR + RACE_CN_HUI)
    )

    # ----- Individual PK parameters -----
    ka      <- exp(lka)
    cl      <- exp(lcl     + etalcl     + e_sexf_cl * SEXF)
    vc      <- exp(lvc     + etalvc) * fethn_vc
    q       <- exp(lq      + etalq      + e_sexf_q  * SEXF)
    vp      <- exp(lvp     + etalvp     + e_sexf_vp * SEXF)
    cl_mfa  <- exp(lcl_mfa + etalcl_mfa) * fethn_cl_mfa
    vc_mfa  <- exp(lvc_mfa + etalvc_mfa)

    # ----- Micro rate constants -----
    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_mfa <- cl_mfa / vc_mfa

    # ----- ODE system (Wu 2012 Figure 1 4-compartment GI / 2-cmt modafinil / 1-cmt acid structure) -----
    # All modafinil elimination flux is treated as forming modafinil acid at the
    # apparent-parameter level; F2 is absorbed into apparent V3 / CL3 so the
    # simulator reproduces the observed acid time course (see Implementation notes).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1
    d/dt(central_mfa) <-  kel * central - kel_mfa * central_mfa

    # ----- Observations (umol/L; central in mg of modafinil, vc in L) -----
    # Cc:     central [mg of modafinil] / vc [L] * 1000 / MW_modafinil = umol of modafinil / L.
    # Cc_mfa: central_mfa [mg of modafinil equivalent, 1:1 molar with modafinil acid] / vc_mfa [L] * 1000 / MW_modafinil
    #         = umol of modafinil acid / L.
    Cc     <- central     / vc     * 1000 / mw_modafinil
    Cc_mfa <- central_mfa / vc_mfa * 1000 / mw_modafinil

    Cc     ~ add(addSd) + prop(propSd)
    Cc_mfa ~ prop(propSd_mfa)
  })
}
