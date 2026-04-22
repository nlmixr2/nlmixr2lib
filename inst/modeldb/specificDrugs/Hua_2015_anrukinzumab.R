Hua_2015_anrukinzumab <- function() {
  description <- "Two-compartment population PK model for anrukinzumab (anti-IL-13 IgG1 monoclonal antibody) with first-order SC absorption and linear elimination, pooling healthy volunteers, mild-to-moderate asthma, moderate-to-severe asthma, and ulcerative colitis patients (Hua 2015)"
  reference <- "Hua F, Ribbing J, Reinisch W, Cataldi F, Martin S. A pharmacokinetic comparison of anrukinzumab, an anti-IL-13 monoclonal antibody, among healthy volunteers, asthma and ulcerative colitis patients. Br J Clin Pharmacol. 2015;80(1):101-109. doi:10.1111/bcp.12589"
  vignette <- "Hua_2015_anrukinzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling (WT/75)^0.75 on CL (exponent fixed per Hua 2015 base-model decision) and (WT/75)^0.688 on both Vc and Vp (estimated shared exponent). Reference 75 kg per Table 3.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL: (ALB/4.3)^theta_albumin with theta_albumin = -1.07 (Hua 2015 Table 3). Reference 4.3 g/dL per Table 3. Albumin units inferred as g/dL from Figure 1 covariate-plot axis range (3.0-5.0) which is consistent with US-convention reporting.",
      source_name        = "ALB"
    ),
    DIS_UC = list(
      description        = "Ulcerative colitis patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-UC: healthy volunteers or asthma patients in the Hua 2015 pooled cohort)",
      notes              = "Multiplicative fractional effect on CL: CL_UC = CL_nonUC * (1 + theta_UC * DIS_UC), theta_UC = 0.728 (Hua 2015 Table 3). UC patients in Hua 2015 received IV doses only (study 5), so F is not defined separately for UC.",
      source_name        = "UC"
    ),
    DIS_SASTHMA = list(
      description        = "Moderate-to-severe asthma patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer, mild-to-moderate asthma, or UC)",
      notes              = "Multiplicative fractional effect on SC bioavailability: F = F_pop * (1 + theta_sA * DIS_SASTHMA), theta_sA = -0.309 (Hua 2015 Table 3). Moderate-to-severe asthma defined per Hua 2015 study 4 inclusion: FEV1 55-80% and ACQ-5 >= 2. Mild-to-moderate asthma (study 1) used FEV1 > 70% and ACQ-5 <= 1. The reduced F in moderate-to-severe asthma is discussed as potentially reflecting sparse-sampling / SC-dosing-accuracy effects in the phase-2 multicenter protocol rather than a disease-biology effect.",
      source_name        = "sAsthma"
    )
  )

  population <- list(
    n_subjects     = 255L,
    n_studies      = 5L,
    age_median     = "37 years (mean 38, SD 13)",
    weight_median  = "81.3 kg (mean 82.6, SD 18.7)",
    sex_female_pct = 35,
    race_ethnicity = c(White = 73, Black = 12, Asian = 12, Other = 3),
    disease_state  = "Healthy volunteers (17%), mild-to-moderate asthma (20%), moderate-to-severe asthma (38%), ulcerative colitis (25%).",
    dose_range     = "0.2-600 mg total dose across 5 studies; SC (72% of subjects) and IV (28%) routes.",
    regions        = "Multinational; study 2 enrolled Japanese and non-Asian healthy volunteers separately.",
    ada_status     = "No anti-anrukinzumab ADA reported in any of the five studies.",
    notes          = "Baseline demographics pooled from Hua 2015 Table 2 (n = 255 across all 5 studies). Study breakdown from Table 1: (1) single-dose escalation in mild-to-moderate asthma, n=37, SC/IV 0.3-4 mg/kg; (2) single-dose SC escalation in healthy Japanese and non-Asian volunteers, n=44, 0.3-4 mg/kg; (3) allergen-challenge study in mild asthma, n=14, 2 mg/kg SC x 2 doses; (4) moderate-to-severe asthma efficacy, n=97, SC 0.2-2 mg/kg or 200 mg fixed-dose Q2-4W; (5) UC biomarker study, n=63, IV 200-600 mg Q2-4W x 5 doses. Reference covariate values: WT = 75 kg, ALB = 4.3 g/dL, non-UC, non-moderate-to-severe-asthma."
  )

  ini({
    # Structural PK parameters -- reference values for a 75-kg non-UC subject
    # with baseline albumin 4.3 g/dL (Hua 2015 Table 3, Final model column).
    lka     <- log(0.0119); label("First-order SC absorption rate constant (Ka, 1/h)")                         # Hua 2015 Table 3: Ka,pop = 0.0119 /h
    lcl     <- log(0.00732); label("Systemic clearance for 75-kg non-UC subject at ALB 4.3 g/dL (CL, L/h)")    # Hua 2015 Table 3: CL,pop = 0.00732 L/h
    lvc     <- log(3.81);   label("Central volume of distribution for 75-kg subject (Vc, L)")                  # Hua 2015 Table 3: Vc,pop = 3.81 L
    lvp     <- log(2.17);   label("Peripheral volume of distribution for 75-kg subject (Vp, L)")               # Hua 2015 Table 3: Vp,pop = 2.17 L
    lq      <- log(0.0224); label("Inter-compartmental clearance (Q, L/h)")                                    # Hua 2015 Table 3: Q,pop = 0.0224 L/h
    lfdepot <- log(0.973);  label("SC bioavailability in non-moderate-to-severe-asthma subjects (F, log-scale)")  # Hua 2015 Table 3: F,pop = 0.973

    # Allometric exponents on body weight (reference 75 kg; Hua 2015 Table 3).
    allo_cl <- fixed(0.75); label("Allometric exponent of WT on CL (fixed at 0.75)")                           # Hua 2015 Table 3 / Methods: fixed to 0.75
    allo_v  <- 0.688;       label("Allometric exponent of WT on Vc and Vp (estimated, shared)")                # Hua 2015 Table 3: theta_WT (Vc,Vp) = 0.688

    # Covariate effects on CL and F (Hua 2015 Table 3 final model).
    e_alb_cl <- -1.07;  label("Power exponent of baseline albumin on CL ((ALB/4.3)^e_alb_cl)")                 # Hua 2015 Table 3: theta_albumin = -1.07
    e_uc_cl  <-  0.728; label("Fractional increase in CL for UC patients vs. non-UC (unitless)")               # Hua 2015 Table 3: f0_UC = 0.728
    e_sa_f   <- -0.309; label("Fractional change in SC bioavailability for moderate-to-severe asthma (unitless)")  # Hua 2015 Table 3: theta_sA = -0.309

    # IIV (log-normal on the exponential-IIV estimation scale).
    # omega^2 = log(1 + CV^2) from reported CV% (Hua 2015 Table 3).
    # Hua 2015: "Vc and Vp were assigned the same eta values and assumed full correlation",
    # so a single eta (etalvc) drives both Vc and Vp in model() below. The block pairs
    # it with etalcl with correlation 0.727.
    #   var_lcl = log(1 + 0.316^2) = 0.0952
    #   var_lvc = log(1 + 0.265^2) = 0.0679
    #   cov     = 0.727 * sqrt(0.0952 * 0.0679) = 0.0584
    #   var_lka = log(1 + 0.540^2) = 0.2560
    etalcl + etalvc ~ c(0.0952,
                        0.0584, 0.0679)   # Hua 2015 Table 3: CV_CL=31.6%, CV_Vc,Vp=26.5%, corr=0.727
    etalka ~ 0.2560                       # Hua 2015 Table 3: CV_Ka=54.0%

    # Residual error -- additive on log-transformed concentrations per Hua 2015
    # Methods, which corresponds to proportional error in nlmixr2 linear-space syntax.
    propSd <- 0.235; label("Proportional residual error (fraction)")                                           # Hua 2015 Table 3: residual variability = 23.5%
  })
  model({
    # Covariate model (Hua 2015 Table 3 Final model equations, reproduced in full):
    #   CL = CL_pop * (WT/75)^0.75 * (ALB/4.3)^theta_albumin * (1 + theta_UC * DIS_UC)
    #   Vc = Vc_pop * (WT/75)^theta_WT_V
    #   Vp = Vp_pop * (WT/75)^theta_WT_V          # same eta as Vc (full correlation)
    #   Q  = Q_pop
    #   Ka = Ka_pop
    #   F  = F_pop * (1 + theta_sA * DIS_SASTHMA) # applied at depot (SC); UC was IV-only
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT / 75)^allo_cl *
      (ALB / 4.3)^e_alb_cl *
      (1 + e_uc_cl * DIS_UC)
    vc <- exp(lvc + etalvc) * (WT / 75)^allo_v
    vp <- exp(lvp + etalvc) * (WT / 75)^allo_v
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot) * (1 + e_sa_f * DIS_SASTHMA)

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
