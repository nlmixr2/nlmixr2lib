Zhang_2025_dupilumab_fev1 <- function() {
  description <- "Semi-mechanistic population PK/PD model for the effect of dupilumab on pre-bronchodilator FEV1 in adult and adolescent patients with uncontrolled moderate-to-severe asthma (Zhang 2025), combining the Zhang 2021 two-compartment asthma popPK layer (first-order SC absorption, parallel linear plus Michaelis-Menten elimination) with a direct-response Emax drug effect, an empirical exponential-onset placebo effect, and an additive baseline FEV1."
  reference <- paste(
    "Zhang L, Davis JD, Kanamaluru V, Xu C.",
    "Semi-mechanistic population pharmacokinetic/pharmacodynamic (PK/PD) modeling of dupilumab",
    "on pre-bronchodilator forced expiratory volume in 1 second (FEV1) in uncontrolled",
    "moderate-to-severe asthma.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1370-1380. doi:10.1002/psp4.70057.",
    "PK layer fixed from Zhang L, Gao Y, Li M, et al.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(8):941-952; doi:10.1002/psp4.12667;",
    "see modellib('Zhang_2021_dupilumab').",
    sep = " "
  )
  vignette <- "Zhang_2025_dupilumab_fev1"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dupilumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "dupilumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the model twice, at two different reference weights, because the PK and PD layers",
        "were fitted against different datasets. PK layer (fixed from Zhang 2021): power effects on",
        "Ke, V2 and Vmax normalised to 78 kg. PD layer (Zhang 2025 Equation 13): power effect",
        "(WT / 77)^0.259 on baseline FEV1, 77 kg being the median of the Zhang 2025 pooled Phase 2b +",
        "Phase 3 PK/PD dataset (Zhang 2025 Table 2). The NONMEM source column for the PD-layer effect",
        "is BLWT (baseline weight); the PK-layer weight is time-varying."
      ),
      source_name        = "BLWT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect (AGE / 50)^-0.423 on baseline FEV1 (Zhang 2025 Equation 13); 50 years is the",
        "median of the pooled dataset (Zhang 2025 Table 2). Age had no significant effect on the",
        "dupilumab treatment effect Emax. Source column AGEY."
      ),
      source_name        = "AGEY"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Selects the sex-specific typical baseline FEV1: 1.93 L in men, 1.54 L in women",
        "(Zhang 2025 Table 3 / Equation 13). The NONMEM source column SEX is 1 for female, and the",
        "control stream derives FEMALE = 1 when SEX == 1, so the source coding already matches the",
        "canonical SEXF orientation (no value inversion needed)."
      ),
      source_name        = "SEX"
    ),
    NEXAC12M = list(
      description        = "Number of severe asthma exacerbations in the 12 months before study entry",
      units              = "count",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Power effect (NEXAC12M / 1.0)^-0.0411 on baseline FEV1 (Zhang 2025 Equation 13 and Table 3);",
        "the normalising value 1.0 is the population median (Zhang 2025 Table 2; range 1-50). Used as a",
        "continuous count inside a power function, so it cannot be decomposed into band indicators.",
        "Zhang 2025 Table 4 shows the effect is not clinically meaningful (-7% at the 95th percentile of",
        "6 exacerbations). Source column PREEXAC."
      ),
      source_name        = "PREEXAC"
    ),
    FENO = list(
      description        = "Baseline fractional exhaled nitric oxide concentration",
      units              = "ppb",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Type-2 inflammation biomarker. Power effect (FENO / 25)^0.682 on the dupilumab maximum",
        "treatment effect Emax (Zhang 2025 Equation 12 and Table 3); 25 ppb is the median of the pooled",
        "dataset (Zhang 2025 Table 2; range 3-387 ppb). Baseline (pre-first-dose) value, time-fixed per",
        "subject. Source column BFENO."
      ),
      source_name        = "BFENO"
    ),
    EOS = list(
      description        = "Baseline blood eosinophil count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Type-2 inflammation biomarker. Power effect (EOS / 260)^0.334 on the dupilumab maximum",
        "treatment effect Emax. Zhang 2025 Equation 12 writes the same ratio as (EOS / 0.26) with EOS in",
        "10^9/L; the supplementary NONMEM control stream uses (BEOS / 260), i.e. the dataset column is in",
        "cells/uL. 260 cells/uL == 0.26 x 10^9/L, so the two forms are numerically identical; the",
        "cells/uL form is used here to match the register canonical. Baseline value, time-fixed per",
        "subject. Source column BEOS."
      ),
      source_name        = "BEOS"
    ),
    ALB = list(
      description        = "Serum albumin (baseline)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "PK-layer covariate only, inherited fixed from Zhang 2021: power effect (ALB / 44)^-0.484 on V2.",
        "Zhang 2025 uses 44 g/L as the reference-patient albumin in its covariate simulations",
        "(Zhang 2025 Section 2.4)."
      ),
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Creatinine clearance normalized to body surface area (Cockcroft-Gault, BSA-normalized as 1.73 * CrCl / BSA)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "PK-layer covariate only, inherited fixed from Zhang 2021: power effect (CRCL / 111)^0.217 on Ke.",
        "Zhang 2025 uses 111 mL/min/1.73 m^2 as the reference-patient value in its covariate simulations",
        "(Zhang 2025 Section 2.4). Source column CLCRN."
      ),
      source_name        = "CLCRN"
    ),
    ADA_POS = list(
      description        = "Stationary anti-drug antibody (ADA) positivity indicator (positive at any time on study)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative; typical patient)",
      notes              = paste(
        "PK-layer covariate only, inherited fixed from Zhang 2021: Ke * (1 + 0.191 * ADA_POS).",
        "Zhang 2025 tested stationary ADA status as a covariate on the PD parameters and found no",
        "significant effect on Emax (Zhang 2025 Results / Figure S2); the reference patient in the",
        "Zhang 2025 covariate simulations is ADA-negative."
      ),
      source_name        = "ADA2"
    )
  )

  # Covariates that Zhang 2025 screened on the PD parameters but did not retain in
  # the final model. Documented for provenance only; they are deliberately absent
  # from model().
  covariatesDataExcluded <- list(
    RACE = list(
      description = "Race group",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    REGION = list(
      description = "Geographic enrollment region",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    SMOKE = list(
      description = "Smoking history indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    FEV1 = list(
      description = "Baseline pre-bronchodilator FEV1 supplied as a covariate column",
      units       = "L",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate on Emax; not significant. Baseline FEV1 is a MODEL PARAMETER here",
        "(lrbase_male / lrbase_female with IIV etalrbase), not a covariate input."
      )
    ),
    IGE = list(
      description = "Baseline serum total immunoglobulin E",
      units       = "IU/mL",
      type        = "continuous",
      notes       = "Type-2 biomarker screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    TARC = list(
      description = "Baseline thymus and activation-regulated chemokine",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "Type-2 biomarker screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    POSTN = list(
      description = "Baseline serum periostin",
      units       = "ng/mL",
      type        = "continuous",
      notes       = "Type-2 biomarker screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    SCORE_ACQ5 = list(
      description = "Baseline 5-item asthma control questionnaire score",
      units       = "score",
      type        = "continuous",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    ICSD = list(
      description = "Background inhaled-corticosteroid dose level at randomization",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    AGEOS = list(
      description = "Age at onset of asthma",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    ),
    ATOP = list(
      description = "Atopic medical condition indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on Emax; not significant (Zhang 2025 Results, Figure S2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2654L,
    n_studies      = 2L,
    age_range      = "12-87 years",
    age_median     = "50 years",
    weight_range   = "30-227 kg",
    weight_median  = "77.0 kg",
    sex_female_pct = 62.9,
    race_ethnicity = "Race and geographic region were screened as covariates and were not significant; the race breakdown is not tabulated in the main text.",
    disease_state  = "Uncontrolled, persistent moderate-to-severe asthma treated with a medium-to-high dose of inhaled corticosteroid plus up to two long-acting beta2-agonists; dupilumab or placebo given as add-on maintenance therapy.",
    dose_range     = "Placebo (n = 794), or dupilumab 200 mg SC (400 mg loading dose) or 300 mg SC (600 mg loading dose) every 2 or every 4 weeks (n = 1860).",
    regions        = "Multi-regional Phase 2b (NCT01854047 / DRI12544, 24 weeks) and Phase 3 (NCT02414854 / EFC13579, 52 weeks) placebo-controlled pivotal studies.",
    notes          = paste(
      "Baseline demographics from Zhang 2025 Table 2 (pooled N = 2654; Phase 2b N = 761, Phase 3",
      "N = 1893). Adolescents (12 to <18 years) N = 107 (4.0%), all from the Phase 3 study. Median",
      "(range) baseline blood eosinophil count 0.26 (0-8.75) x 10^9/L = 260 (0-8750) cells/uL; median",
      "(range) FeNO 25 (3-387) ppb; median (range) number of exacerbations in the past year 1 (1-50).",
      "Pre-bronchodilator FEV1 was measured by spirometry per ATS/ERS standards in the morning after",
      "withholding short-acting bronchodilators for >= 6 h and ICS/LABA for 12 h. Serum dupilumab was",
      "measured by validated ELISA with LLOQ 0.078 mg/L. The individual PK parameters that drive the",
      "Emax term were post-hoc estimates from the Zhang 2021 asthma popPK model, read into the PK/PD",
      "dataset as the columns IKEL / IVC / IKCP / IKPC / IVMAX / IKM / IKA / IF1 (Zhang 2025",
      "Supporting Information, $INPUT and $PK blocks)."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # PK layer. Zhang 2025 fitted the PD parameters conditional on per-subject
    # post-hoc PK parameters carried in from the Zhang 2021 asthma popPK model
    # (Zhang 2025 Methods 2.3 and Supporting Information $PK: KEL = IKEL,
    # VC = IVC, ... read from the dataset). Every PK value below is therefore
    # FIXED from that upstream publication rather than estimated here; the
    # numbers are the Zhang 2021 Table 3 final-model estimates, reproduced from
    # modellib("Zhang_2021_dupilumab"). Reference covariate values: 78 kg,
    # 44 g/L albumin, 111 mL/min/1.73 m^2 creatinine clearance, ADA-negative.
    # ------------------------------------------------------------------------
    lka     <- fixed(log(0.263));  label("Absorption rate Ka (1/day)")                                 # Zhang 2021 Table 3, Ka row (upstream, ref [33] of Zhang 2025)
    lkel    <- fixed(log(0.0418)); label("Linear elimination rate Ke (1/day)")                         # Zhang 2021 Table 3, Ke row
    lvc     <- fixed(log(2.76));   label("Central compartment volume V2 (L)")                          # Zhang 2021 Table 3, V2 row
    lk12    <- fixed(log(0.0952)); label("Central-to-peripheral rate K23 (1/day)")                     # Zhang 2021 Table 3, K23 row
    lk21    <- fixed(log(0.163));  label("Peripheral-to-central rate K32 (1/day)")                     # Zhang 2021 Table 3, K32 row
    lvmax   <- fixed(log(1.39));   label("Maximum target-mediated elimination rate Vmax (mg/L/day)")   # Zhang 2021 Table 3, Vmax row
    lkm     <- fixed(log(2.08));   label("Michaelis-Menten constant Km (mg/L)")                        # Zhang 2021 Table 3, Km row
    lfdepot <- fixed(log(0.609));  label("Subcutaneous bioavailability Fsc (fraction)")                # Zhang 2021 Table 3, Fsc row

    e_wt_kel   <- fixed(0.222);  label("Power exponent of WT/78 on Ke (unitless)")                     # Zhang 2021 Table 3, weight effect on Ke
    e_wt_vc    <- fixed(0.667);  label("Power exponent of WT/78 on V2 (unitless)")                     # Zhang 2021 Table 3, weight effect on V2
    e_wt_vmax  <- fixed(0.224);  label("Power exponent of WT/78 on Vmax (unitless)")                   # Zhang 2021 Table 3, weight effect on Vmax
    e_alb_vc   <- fixed(-0.484); label("Power exponent of ALB/44 on V2 (unitless)")                    # Zhang 2021 Table 3, albumin effect on V2
    e_crcl_kel <- fixed(0.217);  label("Power exponent of CRCL/111 on Ke (unitless)")                  # Zhang 2021 Table 3, CLCRN effect on Ke
    e_ada_kel  <- fixed(0.191);  label("Proportional multiplicative effect of ADA-positive on Ke (unitless)")  # Zhang 2021 Table 3, ADA effect on Ke; coded as Ke * (1 + 0.191 * ADA_POS)

    etalkel    ~ fixed(0.0385)    # Zhang 2021 Table 3: Ke IIV omega^2 (CV ~19.6%)
    etalvc     ~ fixed(0.00834)   # Zhang 2021 Table 3: V2 IIV omega^2 (CV ~9.13%)
    etalvmax   ~ fixed(0.0589)    # Zhang 2021 Table 3: Vmax IIV omega^2 (CV ~24.3%)
    etalka     ~ fixed(0.243)     # Zhang 2021 Table 3: Ka IIV omega^2 (CV ~49.2%)
    etalfdepot ~ fixed(0.132)     # Zhang 2021 Table 3: Fsc IIV omega^2 (CV ~36.3%)

    # ------------------------------------------------------------------------
    # PD layer: pre-bronchodilator FEV1 (L). All values are Zhang 2025 Table 3
    # final-model estimates. Structure follows Figure 1 and Equations 1, 2, 12
    # and 13, confirmed line-by-line against the NONMEM control stream in the
    # Supporting Information ("Pop PK/PD model for dupilumab treatment arm").
    # ------------------------------------------------------------------------

    # Baseline FEV1 is reported separately for each sex, so both strata carry an
    # explicit suffix (parameter-names.md, "Stratum-suffixed parameters").
    lrbase_male   <- log(1.93);  label("Typical baseline pre-bronchodilator FEV1 in men (L)")    # Zhang 2025 Table 3, "Typical value of base for men (theta1, L)" = 1.93
    lrbase_female <- log(1.54);  label("Typical baseline pre-bronchodilator FEV1 in women (L)")  # Zhang 2025 Table 3, "Typical value of base for women (L)" = 1.54

    # Placebo effect. Pmax carries a NORMAL (additive) inter-individual random
    # effect in the source control stream -- PMAX = TVPMAX + ETA(2) -- so it is
    # kept on the natural scale rather than log-transformed.
    plbmax <- 0.172;         label("Maximum placebo effect on FEV1 (L)")                         # Zhang 2025 Table 3, "Typical value of Pmax (theta2, L)" = 0.172
    lkplb  <- log(0.0322);   label("First-order rate constant of placebo onset (1/day)")         # Zhang 2025 Table 3, "Typical value of Kplb (theta3, 1/day)" = 0.0322

    # Dupilumab direct-response Emax effect.
    lemax <- log(0.104);  label("Maximum dupilumab treatment effect on FEV1 (L)")                # Zhang 2025 Table 3, "Typical value of Emax (theta5, L)" = 0.104
    lec50 <- log(0.713);  label("Serum dupilumab concentration giving half-maximal effect EC50 (mg/L)")  # Zhang 2025 Table 3, "Typical value of EC50 (theta4, mg/L)" = 0.713

    # Covariate effects on baseline FEV1 (Zhang 2025 Equation 13).
    e_age_rbase      <- -0.423;   label("Power exponent of AGE/50 on baseline FEV1 (unitless)")                       # Zhang 2025 Table 3, "Power coefficient of AGEY on Base"
    e_wt_rbase       <-  0.259;   label("Power exponent of WT/77 on baseline FEV1 (unitless)")                        # Zhang 2025 Table 3, "Power coefficient of WT on Base"
    e_nexac12m_rbase <- -0.0411;  label("Power exponent of NEXAC12M/1 on baseline FEV1 (unitless)")                   # Zhang 2025 Table 3, "Power coefficient of PREEXAC on Base" = -0.0411 (Equation 13 prints the rounded -0.041)

    # Covariate effects on the dupilumab Emax (Zhang 2025 Equation 12).
    e_feno_emax <- 0.682;  label("Power exponent of FENO/25 on Emax (unitless)")                 # Zhang 2025 Table 3, "Power coefficient of FeNO on Emax"
    e_eos_emax  <- 0.334;  label("Power exponent of EOS/260 on Emax (unitless)")                 # Zhang 2025 Table 3, "Power coefficient of EOS on Emax"

    # Inter-individual variability. Zhang 2025 Table 3 reports the variance in the
    # "Estimate" column with the authors' small-variance approximation
    # CV(%) ~ 100 * sqrt(variance) in parentheses; each of the five rows
    # round-trips exactly under that convention (e.g. sqrt(0.0588) = 0.242 -> 24.2%),
    # which confirms the tabulated numbers are variances, not SDs.
    etalrbase ~ 0.0588   # Zhang 2025 Table 3: baseline FEV1 IIV omega^2 (reported 24.2%, shrinkage 5.48%); lognormal, BASE = TVBASE * EXP(ETA(1))
    etaplbmax ~ 0.0924   # Zhang 2025 Table 3: Pmax IIV omega^2 in L^2 (reported 30.4%, shrinkage 8.49%); ADDITIVE eta, PMAX = TVPMAX + ETA(2)
    etalkplb  ~ 1.64     # Zhang 2025 Table 3: Kplb IIV omega^2 (reported 128%, shrinkage 48.7%); lognormal
    etalec50  ~ 3.72     # Zhang 2025 Table 3: EC50 IIV omega^2 (reported 193%, shrinkage 80.3%); lognormal
    etalemax  ~ 0.710    # Zhang 2025 Table 3: Emax IIV omega^2 (reported 84.3%, shrinkage 52.2%); lognormal

    # Residual error on FEV1: combined proportional + additive
    # (Y = IPRED * (1 + ERR(1)) + ERR(2) for CMT 4, Zhang 2025 Supporting
    # Information $ERROR block). Table 3 reports the variances.
    propSd_FEV1 <- 0.0647;  label("Proportional residual error on FEV1 (fraction)")              # Zhang 2025 Table 3: proportional sigma^2 = 0.00419, SD = sqrt(0.00419) = 0.0647 (the tabulated 6.47% is the same number as a percent)
    addSd_FEV1  <- 0.120;   label("Additive residual error on FEV1 (L)")                         # Zhang 2025 Table 3: additive sigma^2 = 0.0144 L^2, SD = sqrt(0.0144) = 0.12 L (the tabulated "(0.12)" is that SD)

    # Residual error on serum dupilumab concentration. Zhang 2025 held both PK
    # sigmas fixed at the Zhang 2021 values while fitting the PD layer
    # ($SIGMA 0.04 FIXED ; PROP FOR PK  /  $SIGMA 2.98 FIXED ; ADD FOR PK).
    propSd <- fixed(0.200);  label("Proportional residual error on serum dupilumab (fraction)")  # Zhang 2025 Supporting Information $SIGMA: PK proportional variance 0.04 FIXED, SD = sqrt(0.04) = 0.200
    addSd  <- fixed(1.73);   label("Additive residual error on serum dupilumab (mg/L)")          # Zhang 2025 Supporting Information $SIGMA: PK additive variance 2.98 FIXED, SD = sqrt(2.98) = 1.73
  })
  model({
    # --- PK layer (Zhang 2021 covariate model, reference patient: 78 kg,
    # 44 g/L albumin, 111 mL/min/1.73 m^2 creatinine clearance, ADA-negative) ---
    kel <- exp(lkel + etalkel) *
           (1 + e_ada_kel * ADA_POS) *
           (WT / 78)^e_wt_kel *
           (CRCL / 111)^e_crcl_kel
    vc  <- exp(lvc + etalvc) *
           (WT / 78)^e_wt_vc *
           (ALB / 44)^e_alb_vc
    vmax <- exp(lvmax + etalvmax) *
            (WT / 78)^e_wt_vmax
    km <- exp(lkm)
    ka <- exp(lka + etalka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot <- exp(lfdepot + etalfdepot)

    # Serum dupilumab concentration Cp in Zhang 2025 Equation 1 / Figure 1.
    Cc <- central / vc

    # Two-compartment PK with first-order SC absorption and parallel linear plus
    # Michaelis-Menten elimination from the central compartment. The source
    # $DES writes the saturable term as A(2) * VMAX / (KM + A(2)/VC), i.e.
    # central * vmax / (km + Cc) with Vmax in mg/L/day.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          kel * central -
                          vmax * central / (km + Cc) -
                          k12 * central +
                          k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # --- PD layer: pre-bronchodilator FEV1 (L) ---
    # Baseline (Zhang 2025 Equation 13). The sex term selects between the two
    # published typical values; the three power terms are normalised to the
    # pooled-dataset medians (50 years, 77 kg, 1 prior exacerbation).
    rbasesex <- exp(lrbase_male) * (1 - SEXF) + exp(lrbase_female) * SEXF
    rbase <- rbasesex *
             (AGE / 50)^e_age_rbase *
             (WT / 77)^e_wt_rbase *
             (NEXAC12M / 1.0)^e_nexac12m_rbase *
             exp(etalrbase)

    # Placebo effect (Zhang 2025 Equation 2). Pmax uses an additive eta.
    plbmaxi <- plbmax + etaplbmax
    kplb <- exp(lkplb + etalkplb)
    plbeff <- plbmaxi * (1 - exp(-kplb * t))

    # Dupilumab direct-response Emax effect (Zhang 2025 Equations 1 and 12).
    # With no dupilumab dosing Cc is 0 and drugeff is 0, which reproduces the
    # source control stream's IF (NDOSE .EQ. 0) EDRUG = 0 branch, so the same
    # model serves the placebo and the dupilumab arms.
    emax <- exp(lemax + etalemax) *
            (FENO / 25)^e_feno_emax *
            (EOS / 260)^e_eos_emax
    ec50 <- exp(lec50 + etalec50)
    drugeff <- emax * Cc / (ec50 + Cc)

    FEV1 <- rbase + plbeff + drugeff

    Cc ~ prop(propSd) + add(addSd)
    FEV1 ~ prop(propSd_FEV1) + add(addSd_FEV1)
  })
}
