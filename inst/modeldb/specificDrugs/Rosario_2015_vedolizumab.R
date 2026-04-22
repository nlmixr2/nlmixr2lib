Rosario_2015_vedolizumab <- function() {
  description <- "Two-compartment population PK model for vedolizumab (humanised anti-alpha4-beta7 integrin IgG1 monoclonal antibody) with parallel linear and Michaelis-Menten elimination in adults with moderately-to-severely active ulcerative colitis or Crohn's disease and healthy volunteers (Rosario 2015)."
  reference   <- "Rosario M, Dirks NL, Gastonguay MR, Fasanmade AA, Wyant T, Parikh A, Sandborn WJ, Feagan BG, Reinisch W, Fox I. Population pharmacokinetics-pharmacodynamics of vedolizumab in patients with ulcerative colitis and Crohn's disease. Aliment Pharmacol Ther. 2015;42(2):188-202. doi:10.1111/apt.13243 (PMID 25996351). A corrigendum (doi:10.1111/apt.15571; PMC6885991) corrects a unit typo in the text (ng/mL -> ug/mL) and does not change any parameter value."
  vignette <- "Rosario_2015_vedolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on CLL (exponent 0.362, estimated) and Vc (exponent 0.467, estimated); fixed allometric exponents on Vp (1), Vmax (0.75), and Q (0.75). Reference 70 kg per Rosario 2015 Table 2 footnote (reference patient).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CLL: (ALB / 4)^(-1.18). Reference 4 g/dL is the reference-patient albumin (Rosario 2015 Table 2 footnote and Figure 5 caption). US-convention g/dL matches the source paper.",
      source_name        = "ALB"
    ),
    CALPRO = list(
      description        = "Baseline faecal calprotectin",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CLL: (CALPRO / 700)^0.0310. Reference 700 mg/kg per Rosario 2015 Table 2 footnote (reference patient).",
      source_name        = "CALPRO"
    ),
    CDAI = list(
      description        = "Crohn's Disease Activity Index (CD patients only)",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CLL for CD patients only: (CDAI / 300)^(-0.0515 * IBD_CD). Reference 300 per Rosario 2015 Table 2 footnote (reference CD patient). Gated by IBD_CD so the effect is identically 1 for UC patients; for CD patients supply the observed CDAI.",
      source_name        = "CDAI"
    ),
    PMAYO = list(
      description        = "Partial Mayo score (UC patients only)",
      units              = "(score 0-9)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CLL for UC patients only: (PMAYO / 6)^(0.0408 * (1 - IBD_CD)). Reference 6 per Rosario 2015 Table 2 footnote (reference UC patient). Gated by (1 - IBD_CD) so the effect is identically 1 for CD patients.",
      source_name        = "PMAYO"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CLL: (AGE / 40)^(-0.0346). Reference 40 years per Rosario 2015 Table 2 footnote (reference patient).",
      source_name        = "AGE"
    ),
    IBD_CD = list(
      description        = "IBD diagnosis indicator (Crohn's disease vs ulcerative colitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ulcerative colitis)",
      notes              = "Switches the typical CLL between UC (0.159 L/day) and CD (0.155 L/day), gates the disease-activity covariates (PMAYO applies only when IBD_CD = 0; CDAI applies only when IBD_CD = 1), and drives a small multiplicative effect on Vc (1.01^IBD_CD). The reference patient for Vc is UC (IBD_CD = 0) per Rosario 2015 Table 2 footnote.",
      source_name        = "DX"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF-alpha antagonist therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (TNF-naive)",
      notes              = "Multiplicative effect on CLL of the form 1.04^PRIOR_TNF (Rosario 2015 Table S4; null effect = 1).",
      source_name        = "PRIOR_TNF"
    ),
    ADA_POS = list(
      description        = "Anti-drug-antibody positivity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Multiplicative effect on CLL of the form 1.12^ADA_POS (Rosario 2015 Table S4; null effect = 1). Rosario 2015 did not find a statistically significant ADA-titre effect and used the binary positivity indicator in the final model.",
      source_name        = "ADA"
    ),
    CONMED_AZA = list(
      description        = "Concomitant azathioprine use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azathioprine)",
      notes              = "Multiplicative effect on CLL of the form 0.998^CONMED_AZA (Rosario 2015 Table S4; null effect = 1).",
      source_name        = "AZA"
    ),
    CONMED_MP = list(
      description        = "Concomitant 6-mercaptopurine use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant 6-MP)",
      notes              = "Multiplicative effect on CLL of the form 1.04^CONMED_MP (Rosario 2015 Table S4; null effect = 1).",
      source_name        = "MP"
    ),
    CONMED_MTX = list(
      description        = "Concomitant methotrexate use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant methotrexate)",
      notes              = "Multiplicative effect on CLL of the form 0.983^CONMED_MTX (Rosario 2015 Table S4; null effect = 1).",
      source_name        = "MTX"
    ),
    CONMED_AMINO = list(
      description        = "Concomitant aminosalicylate therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant aminosalicylate)",
      notes              = "Multiplicative effect on CLL of the form 1.02^CONMED_AMINO (Rosario 2015 Table S4; null effect = 1). Covers 5-ASA and the broader aminosalicylate class per the paper's 'AMINO' label.",
      source_name        = "AMINO"
    )
  )

  population <- list(
    n_subjects     = 2700L,
    n_studies      = 6L,
    age_range      = "18-78 years",
    age_median     = "39 years (IBD patients)",
    weight_range   = "29.4-156 kg",
    weight_median  = "70 kg (reference for allometric and power-form scaling)",
    sex_female_pct = 48,
    race_ethnicity = c(White = 87, Black = 2, Asian = 6, Other = 5),
    disease_state  = "Moderately-to-severely active ulcerative colitis or Crohn's disease (IBD cohort); healthy adult volunteers in the phase 1 study contributed to the linear-CL characterisation.",
    dose_range     = "Single IV doses 0.2-10 mg/kg; multiple IV doses 2.0/6.0/10.0 mg/kg or fixed 300 mg every 4 or 8 weeks (induction wks 0,2; maintenance through week 52 in phase 3).",
    regions        = "North America, Western/Northern Europe, Central Europe, Eastern Europe, Asia/Australia/Africa (phase 3 trials); United States (phase 1); Canada and Russia (phase 2).",
    ada_positive_pct = 4,
    prior_tnf_pct    = 50,
    cd_pct           = 51,
    uc_pct           = 43,
    healthy_pct      = 6,
    notes          = "Pooled phase 1/2/3 vedolizumab dataset: C13009 (phase 1, healthy volunteers), C13002 (phase 2 UC), GEMINI 1 (phase 3 UC), GEMINI 2 (phase 3 CD), GEMINI 3 (phase 3 CD), C13004 (phase 2 CD). Demographics from Rosario 2015 Table 1 (pooled population used for population PK-PD analyses). Approximate n_subjects 2700 is the total PK-evaluable pool reported in Table S1 of Appendix S1."
  )

  ini({
    # Structural parameters -- Rosario 2015 Table 2 (final-model Bayesian medians
    # from 4 MCMC chains). Reference patient: 70 kg, 40 y, albumin 4 g/dL,
    # faecal calprotectin 700 mg/kg, CDAI 300 (CD), partial Mayo 6 (UC),
    # UC diagnosis for Vc, no concomitant therapy, ADA-negative, TNF-naive.
    lcl    <- log(0.159);  label("Linear clearance CLL for the reference UC patient (L/day)")  # Table 2: CLL UC = 0.159 L/day
    lcl_cd <- log(0.155);  label("Linear clearance CLL for the reference CD patient (L/day)")  # Table 2: CLL CD = 0.155 L/day
    lvc    <- log(3.19);   label("Central volume of distribution Vc for the reference UC patient (L)")  # Table 2: Vc = 3.19 L
    lvp    <- log(1.65);   label("Peripheral volume of distribution Vp for 70 kg (L)")          # Table 2: Vp = 1.65 L
    lq     <- log(0.12);   label("Intercompartmental clearance Q for 70 kg (L/day)")            # Table 2: Q = 0.12 L/day
    lvmax  <- log(0.265);  label("Maximum elimination rate of the Michaelis-Menten pathway for 70 kg (Vmax, mg/day)")  # Table 2: Vmax = 0.265 mg/day
    lkm    <- log(0.964);  label("Michaelis-Menten constant (Km, ug/mL)")                       # Table 2: Km = 0.964 ug/mL

    # Continuous covariate effects (power-form; Table S4 "NULL effect = 0" entries).
    e_wt_cl     <-  0.362;  label("Weight exponent on CLL (unitless; reference 70 kg)")          # Table S4: weight on CLL = 0.362
    e_alb_cl    <- -1.18;   label("Albumin exponent on CLL (unitless; reference 4 g/dL)")        # Table S4: albumin on CLL = -1.18
    e_calpro_cl <-  0.0310; label("Faecal calprotectin exponent on CLL (unitless; reference 700 mg/kg)")  # Table S4: calprotectin on CLL = 0.0310
    e_cdai_cl   <- -0.0515; label("CDAI exponent on CLL, CD patients only (unitless; reference 300)")      # Table S4: CDAI on CLL = -0.0515
    e_pmayo_cl  <-  0.0408; label("Partial Mayo exponent on CLL, UC patients only (unitless; reference 6)")  # Table S4: pMayo on CLL = 0.0408
    e_age_cl    <- -0.0346; label("Age exponent on CLL (unitless; reference 40 years)")          # Table S4: age on CLL = -0.0346
    e_wt_vc     <-  0.467;  label("Weight exponent on Vc (unitless; reference 70 kg)")           # Table S4: weight on Vc = 0.467

    # Fixed allometric exponents (Table S4: "1 Fixed", "0.75 Fixed").
    allo_wt_vp   <- fixed(1);    label("Allometric exponent of WT on Vp (fixed)")                # Table S4: weight on Vp = 1 Fixed
    allo_wt_vmax <- fixed(0.75); label("Allometric exponent of WT on Vmax (fixed)")              # Table S4: weight on Vmax = 0.75 Fixed
    allo_wt_q    <- fixed(0.75); label("Allometric exponent of WT on Q (fixed)")                 # Table S4: weight on Q = 0.75 Fixed

    # Categorical covariate multipliers (Table S4 "NULL effect = 1").
    e_priortnf_cl   <- 1.04;  label("Prior TNF-alpha antagonist multiplier on CLL (power form: CLL * 1.04^PRIOR_TNF)")   # Table S4
    e_ada_cl        <- 1.12;  label("ADA-positive multiplier on CLL (power form: CLL * 1.12^ADA_POS)")                    # Table S4
    e_conmed_aza_cl <- 0.998; label("Concomitant azathioprine multiplier on CLL (power form: CLL * 0.998^CONMED_AZA)")    # Table S4
    e_conmed_mp_cl  <- 1.04;  label("Concomitant 6-MP multiplier on CLL (power form: CLL * 1.04^CONMED_MP)")              # Table S4
    e_conmed_mtx_cl <- 0.983; label("Concomitant methotrexate multiplier on CLL (power form: CLL * 0.983^CONMED_MTX)")    # Table S4
    e_conmed_amino_cl <- 1.02; label("Concomitant aminosalicylate multiplier on CLL (power form: CLL * 1.02^CONMED_AMINO)")  # Table S4
    e_ibd_cd_vc     <- 1.01;  label("Crohn's-vs-UC multiplier on Vc (power form: Vc * 1.01^IBD_CD)")                      # Table S4

    # Interindividual variability -- Rosario 2015 Table S2 and Table S3.
    # Table S2 reports %CV = 100 * omega (lognormal SD); therefore:
    #   omega2_CLL  = 0.346^2 = 0.119716  (35% CV per Table 2 / Table S2)
    #   omega2_Vc   = 0.191^2 = 0.036481  (19% CV per Table 2 / Table S2)
    #   omega2_Vmax = 1.05^2  = 1.1025    (105% CV per Table S2)
    # Correlations from Table S2: corr(CLL, Vc) = +0.566; corr(CLL, Vmax) = -0.192;
    # corr(Vc, Vmax) = -0.267. Block covariances:
    #   cov(CLL, Vc)   = 0.566 * 0.346 * 0.191 =  0.037402
    #   cov(CLL, Vmax) = -0.192 * 0.346 * 1.05 = -0.069754
    #   cov(Vc, Vmax)  = -0.267 * 0.191 * 1.05 = -0.053541
    # IIV on Vp, Q, and Km was fixed to 0 per Table S2.
    etalcl + etalvc + etalvmax ~ c(0.119716,
                                   0.037402,  0.036481,
                                  -0.069754, -0.053541, 1.1025)  # Rosario 2015 Tables S2/S3

    # Residual error -- Rosario 2015 Table 2: sigma^2_prop = 0.0554 (%CV = 23.5).
    propSd <- sqrt(0.0554); label("Proportional residual error on vedolizumab concentration (fraction)")  # Table 2
  })

  model({
    # Two typical values for linear clearance, switched by the IBD diagnosis indicator.
    cl_typ <- exp(lcl) * (1 - IBD_CD) + exp(lcl_cd) * IBD_CD

    # Individual PK parameters. Reference 70 kg for allometric / power scaling on WT;
    # reference 4 g/dL for ALB; reference 700 mg/kg for CALPRO; reference 300 for CDAI
    # (CD only); reference 6 for PMAYO (UC only); reference 40 y for AGE.
    cl <- cl_typ * exp(etalcl) *
      (WT  / 70  )^e_wt_cl *
      (ALB / 4   )^e_alb_cl *
      (CALPRO / 700)^e_calpro_cl *
      (CDAI   / 300)^(e_cdai_cl  *      IBD_CD ) *
      (PMAYO  / 6  )^(e_pmayo_cl * (1 - IBD_CD)) *
      (AGE / 40  )^e_age_cl *
      e_priortnf_cl^PRIOR_TNF *
      e_ada_cl^ADA_POS *
      e_conmed_aza_cl^CONMED_AZA *
      e_conmed_mp_cl^CONMED_MP *
      e_conmed_mtx_cl^CONMED_MTX *
      e_conmed_amino_cl^CONMED_AMINO

    vc   <- exp(lvc   + etalvc)   * (WT / 70)^e_wt_vc * e_ibd_cd_vc^IBD_CD
    vp   <- exp(lvp)              * (WT / 70)^allo_wt_vp
    q    <- exp(lq)               * (WT / 70)^allo_wt_q
    vmax <- exp(lvmax + etalvmax) * (WT / 70)^allo_wt_vmax
    km   <- exp(lkm)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment model, IV infusion into the central compartment, parallel
    # linear + Michaelis-Menten elimination (Rosario 2015 Figure 2).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 -
                         vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
