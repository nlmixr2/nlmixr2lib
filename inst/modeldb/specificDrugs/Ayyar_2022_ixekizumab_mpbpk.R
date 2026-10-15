Ayyar_2022_ixekizumab_mpbpk <- function() {
  description <- "mPBPK. Minimal physiologically based PK model of ixekizumab (anti-IL-17A monoclonal antibody) in adults with moderate-to-severe plaque psoriasis, with quasi-equilibrium TMDD binding of IL-17A in serum and in skin interstitial fluid, plus the TE-based MBMA linking the 12-week average free skin IL-17A (% baseline) to placebo-adjusted PASI75 / PASI90 response (Ayyar 2022). Naive-pooled fit to published mean data; no between-subject variability."
  reference <- "Ayyar VS, Lee JB, Wang W, Pryor M, Zhuang Y, Wilde T, Vermeulen A. Minimal Physiologically-Based Pharmacokinetic (mPBPK) Metamodeling of Target Engagement in Skin Informs Anti-IL17A Drug Development in Psoriasis. Front Pharmacol. 2022;13:862291. doi:10.3389/fphar.2022.862291"
  vignette <- "Ayyar_2022_il17a_target_engagement"
  units <- list(
    time = "day",
    dosing = "pmol (convert mg via dose_pmol = dose_mg * 1e9 / 150000; nominal 150 kDa IgG -- the paper does not state the molecular weight it used)",
    concentration = "pM"
  )

  # The skin-ISF total-target pool is a paper-anatomical state of this model
  # (Supplementary Eq. 8). The 12-week average free skin IL-17A that drives
  # the TE-based MBMA (Results 'Target Engagement Model-Based Meta-Analysis')
  # is read from a cumulative accumulator of free skin IL-17A (% baseline).
  paper_specific_compartments <- c("total_target_skin", "auc_free_target_skin")

  compartmentData <- list(
    depot = list(analyte = "ixekizumab", units = "pmol", specimen = "administration site", verified = TRUE),
    plasma = list(
      analyte = "ixekizumab (total: free + IL-17A-bound)",
      units = "pmol",
      specimen = "serum",
      verified = TRUE
    ),
    is_muscle = list(analyte = "ixekizumab", units = "pmol", specimen = "tissue", verified = TRUE),
    leaky = list(analyte = "ixekizumab", units = "pmol", specimen = "tissue", verified = TRUE),
    is_skin = list(
      analyte = "ixekizumab (total: free + IL-17A-bound)",
      units = "pmol",
      specimen = "tissue",
      verified = TRUE
    ),
    lymph = list(analyte = "ixekizumab", units = "pmol", specimen = "lymph", verified = TRUE),
    total_target = list(
      analyte = "IL-17A (total: free + ixekizumab-bound)",
      units = "pM",
      specimen = "serum",
      verified = TRUE
    ),
    total_target_skin = list(
      analyte = "IL-17A (total: free + ixekizumab-bound)",
      units = "pM",
      specimen = "tissue",
      verified = TRUE
    ),
    auc_free_target_skin = list(
      analyte = "free IL-17A in skin ISF, cumulative",
      units = "%*day",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    age_range = "adults (not reported at the summary level used)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "Moderate-to-severe plaque psoriasis",
    dose_range = "15 mg IV and 5-150 mg SC at weeks 0, 2, 4 (Phase 1); 10-150 mg SC at weeks 0, 2, 4, 8 (Phase 2); 160 mg SC then 80 mg SC q2w or q4w (Phase 3) (Table 1, Figure 4)",
    regions = "international (Phase 1-3 trials)",
    notes = paste(
      "Model fitted to group-mean serum ixekizumab concentrations digitised from the",
      "Phase 1 dose-escalation study (FDA 2016, N = 46) and the Phase 3 UNCOVER-1/2/3",
      "studies (Papp 2018) (Figure 4); the Phase 2 dose-ranging study (Leonardi 2012,",
      "N = 141) contributed clinical response to the MBMA. The IL-17A system parameters,",
      "the serum complex elimination rate and the skin reflection coefficient were taken",
      "from the secukinumab fit (Table 2 footnote b), which used serum total IL-17A and",
      "skin drug data that are not available for ixekizumab. Naive-pooled fit in NONMEM 7.4",
      "(FOCE-I) with the omega matrix fixed to zero, so the model has no between-subject",
      "variability."
    )
  )

  ini({
    # Structure: Supplementary Material DataSheet1 Eq. 1-12 and the NONMEM control
    # stream in DataSheet2 ($DES block). Values: Table 2 'Ixekizumab' column (the system parameters are shared
    # with secukinumab);
    # where the control stream carries an unrounded fixed value it is used and
    # the Table 2 rounding is noted.

    # ---- Drug-specific PK (Table 2) ----
    lcl <- log(0.379); label("Linear clearance of free ixekizumab from serum CL (L/day)") # Table 2: CL = 0.379 L/day (RSE 8.0%)
    lfdepot <- fixed(log(0.81)); label("SC bioavailability F (fraction)") # Table 2: F = 0.81 (fixed; FDA 2016 / Jackson 2022 popPK)
    lka <- fixed(log(0.24)); label("First-order SC absorption rate constant ka (1/day)") # Table 2: ka = 0.24 1/day (fixed; FDA 2016 / Jackson 2022 popPK)

    # ---- Physiological flows and volumes (Table 2; Cao 2013, Shah and Betts 2012) ----
    llymphflow <- fixed(log(2.9)); label("Total lymph flow rate L (L/day)") # Table 2: L, total = 2.9 L/day (fixed)
    llymphflow_skin <- fixed(log(0.247)); label("Skin lymph flow rate Ls (L/day)") # Table 2: Ls, skin = 0.247 L/day (fixed)
    llymphflow_muscle <- fixed(log(0.71)); label("Muscle lymph flow rate L1 (L/day)") # Table 2: L1, muscle = 0.71 L/day (fixed)
    llymphflow_leaky <- fixed(log(1.943)); label("Leaky-tissue lymph flow rate L2 (L/day)") # Table 2: L2, leaky tissue = 1.943 L/day (fixed)
    lvc <- fixed(log(2.6)); label("Plasma volume Vp (L)") # Table 2: Vp, plasma = 2.6 L (fixed)
    lvleaky <- fixed(log(4.368)); label("Leaky-tissue ISF volume V2 (L)") # DataSheet2 THETA(9) V2 = 4.368 FIX; Table 2 prints 4.37 (rounded)
    lvskin <- fixed(log(1.81)); label("Skin ISF volume Vs (L)") # Table 2: Vs, skin = 1.81 L (fixed)
    lvmuscle <- fixed(log(6.3)); label("Muscle ISF volume V1 (L)") # Table 2: V1, muscle = 6.3 L (fixed)
    lvlymph <- fixed(log(2.6)); label("Lymph volume VL (L)") # Table 2: VL, lymph = 2.6 L (fixed)

    # ---- Vascular and lymphatic reflection coefficients (Table 2) ----
    sigma_skin <- fixed(0.63); label("Vascular reflection coefficient for skin sigma_s (unitless)") # Table 2: sigma_s, skin = 0.63, footnote b (assumed the same as for secukinumab)
    sigma_muscle <- fixed(0.95); label("Vascular reflection coefficient for muscle sigma_1 (unitless)") # Table 2: sigma_1, muscle = 0.95 (fixed)
    sigma_leaky <- 0.524; label("Vascular reflection coefficient for leaky tissues sigma_2 (unitless)") # Table 2: sigma_2, leaky = 0.524 (RSE 17.4%)
    sigma_l <- fixed(0.2); label("Lymphatic capillary reflection coefficient sigma_L (unitless)") # Table 2: sigma_L, lymph = 0.2 (fixed)

    # ---- IL-17A turnover and binding (Table 2) ----
    lkdeg <- fixed(log(45.5)); label("Elimination rate constant of free IL-17A in serum kdeg,p (1/day)") # Table 2: kdeg IL-17A Plasma = 45.5 1/day (fixed; Zheng 2020b)
    lkdeg_skin <- fixed(log(2.44)); label("Elimination rate constant of free IL-17A in skin kdeg,sk (1/day)") # Table 2: kdeg IL-17A Skin = 2.44 1/day (= ksyn / baseline skin)
    lkint <- fixed(log(1.24)); label("Elimination rate constant of the ixekizumab-IL-17A complex in serum kint,p (1/day)") # Table 2: kint complex Plasma = 1.24, footnote b (assumed the same as for secukinumab)
    lkint_skin <- fixed(log(0.34)); label("Elimination rate constant of the ixekizumab-IL-17A complex in skin kint,sk (1/day)") # Table 2: kint complex Skin = 0.34 1/day (= 2.5 x Ls/Vs)
    lr0 <- fixed(log(0.015)); label("Baseline IL-17A concentration in serum (pM)") # Table 2: baseline IL-17A Plasma = 0.015 pM (fixed; Dragatin 2016)
    lr0_skin <- fixed(log(0.28)); label("Baseline IL-17A concentration in skin ISF (pM)") # Table 2: baseline IL-17A Skin = 0.28 pM (fixed; Dragatin 2016)
    lkd <- fixed(log(1.8)); label("Equilibrium dissociation constant of ixekizumab for IL-17A KD (pM)") # Table 2: KD = 1.8 pM (fixed; Adams 2020)

    # ---- TE-based MBMA (Methods Eq. 1; Results 'Target Engagement Model-Based Meta-Analysis') ----
    # x = 12-week average predicted free IL-17A in skin (% baseline). One curve was
    # fitted jointly to secukinumab and ixekizumab arms, so these values are shared
    # with Ayyar_2022_secukinumab_mpbpk. The paper prints no values: each set was
    # digitised by the maintainers from the Figure 6 trend line (Figure 6 PDF image,
    # 500+ points on the plotted range 0.1-23% baseline) and Eq. 1 refitted by nls;
    # the refit reproduces the digitised line with residual SD 0.14 (PASI75) and
    # 0.18 (PASI90) percentage points.
    e0_pasi75 <- 78.11; label("TE-MBMA PASI75 response as free skin IL-17A tends to 0 E0 (% placebo-adjusted)") # digitised Figure 6 left panel, Eq. 1 refit
    emax_pasi75 <- -26.94; label("TE-MBMA PASI75 asymptote at high free skin IL-17A Emax (% placebo-adjusted)") # digitised Figure 6 left panel, Eq. 1 refit
    lec50_pasi75 <- log(8.851); label("TE-MBMA PASI75 free skin IL-17A at half-maximal change E50 (% baseline)") # digitised Figure 6 left panel, Eq. 1 refit
    hill_pasi75 <- 1.094; label("TE-MBMA PASI75 Hill coefficient (unitless)") # digitised Figure 6 left panel, Eq. 1 refit
    e0_pasi90 <- 66.46; label("TE-MBMA PASI90 response as free skin IL-17A tends to 0 E0 (% placebo-adjusted)") # digitised Figure 6 right panel, Eq. 1 refit
    emax_pasi90 <- -6.839; label("TE-MBMA PASI90 asymptote at high free skin IL-17A Emax (% placebo-adjusted)") # digitised Figure 6 right panel, Eq. 1 refit
    lec50_pasi90 <- log(3.062); label("TE-MBMA PASI90 free skin IL-17A at half-maximal change E50 (% baseline)") # digitised Figure 6 right panel, Eq. 1 refit
    hill_pasi90 <- 1.262; label("TE-MBMA PASI90 Hill coefficient (unitless)") # digitised Figure 6 right panel, Eq. 1 refit

    # ---- Residual error (Methods 'Data Analysis and Software') ----
    # Only serum ixekizumab concentrations were available (Figure 4); the fit used
    # log-transformed concentrations ('proportional error model', as in the
    # secukinumab control stream DataSheet2). The magnitude is not reported.
    expSd <- fixed(0); label("Log-scale residual error, serum ixekizumab (not reported)") # Methods: proportional error model; final value not reported
  })

  model({
    cl <- exp(lcl)
    ka <- exp(lka)
    lymphflow <- exp(llymphflow)
    lymphflow_skin <- exp(llymphflow_skin)
    lymphflow_muscle <- exp(llymphflow_muscle)
    lymphflow_leaky <- exp(llymphflow_leaky)
    vc <- exp(lvc)
    vleaky <- exp(lvleaky)
    vskin <- exp(lvskin)
    vmuscle <- exp(lvmuscle)
    vlymph <- exp(lvlymph)

    kdeg <- exp(lkdeg)
    kdeg_skin <- exp(lkdeg_skin)
    kint <- exp(lkint)
    kint_skin <- exp(lkint_skin)
    r0 <- exp(lr0)
    r0_skin <- exp(lr0_skin)
    kd <- exp(lkd)
    # DataSheet2 $PK: ksynSE = BSse * kdegSE, ksynSK = BSsk * kdegSK (both 0.683 pM/day, Table 2)
    ksyn <- kdeg * r0
    ksyn_skin <- kdeg_skin * r0_skin

    total_target(0) <- r0
    total_target_skin(0) <- r0_skin

    # Total drug concentrations (pM = pmol/L)
    cp <- plasma / vc
    cmuscle <- is_muscle / vmuscle
    cleaky <- leaky / vleaky
    cskin <- is_skin / vskin
    clymph <- lymph / vlymph

    # Quasi-equilibrium free drug (Supplementary Eq. 9-10). Written in the
    # algebraically identical form that avoids cancellation in each sign branch
    # of a = C - KD - Rtot (the printed 0.5 * (a + sqrt(a^2 + 4 KD C)) loses all
    # precision when a is large and negative, i.e. at low drug concentration).
    a_p <- cp - kd - total_target
    sq_p <- sqrt(a_p * a_p + 4 * kd * cp)
    if (a_p >= 0) {
      cfree <- 0.5 * (a_p + sq_p)
    } else {
      cfree <- 2 * kd * cp / (sq_p - a_p)
    }
    a_sk <- cskin - kd - total_target_skin
    sq_sk <- sqrt(a_sk * a_sk + 4 * kd * cskin)
    if (a_sk >= 0) {
      cfree_skin <- 0.5 * (a_sk + sq_sk)
    } else {
      cfree_skin <- 2 * kd * cskin / (sq_sk - a_sk)
    }

    # Drug-IL-17A complex and free IL-17A (Supplementary Eq. 11-12)
    ar <- total_target * cfree / (kd + cfree)
    ar_skin <- total_target_skin * cfree_skin / (kd + cfree_skin)
    free_target <- total_target - ar
    free_target_skin <- total_target_skin - ar_skin

    # ODEs in amount form (pmol): DataSheet2 $DES DADT(1)-DADT(6), each multiplied by its volume
    d/dt(depot) <- -ka * depot
    d/dt(plasma) <- ka * depot + lymphflow * clymph -
      cfree * lymphflow_muscle * (1 - sigma_muscle) -
      cfree * lymphflow_leaky * (1 - sigma_leaky) -
      cfree * lymphflow_skin * (1 - sigma_skin) -
      cl * cfree -
      kint * ar * vc
    d/dt(is_muscle) <- cfree * lymphflow_muscle * (1 - sigma_muscle) - cmuscle * lymphflow_muscle * (1 - sigma_l)
    d/dt(leaky) <- cfree * lymphflow_leaky * (1 - sigma_leaky) - cleaky * lymphflow_leaky * (1 - sigma_l)
    d/dt(is_skin) <- cfree * lymphflow_skin * (1 - sigma_skin) -
      cfree_skin * lymphflow_skin * (1 - sigma_l) -
      kint_skin * ar_skin * vskin
    d/dt(lymph) <- cfree_skin * lymphflow_skin * (1 - sigma_l) +
      cmuscle * lymphflow_muscle * (1 - sigma_l) +
      cleaky * lymphflow_leaky * (1 - sigma_l) -
      lymphflow * clymph
    # Total IL-17A (pM): DataSheet2 DADT(7)-DADT(8), Supplementary Eq. 7-8
    d/dt(total_target) <- ksyn - kdeg * (total_target - ar) - kint * ar
    d/dt(total_target_skin) <- ksyn_skin - kdeg_skin * (total_target_skin - ar_skin) - kint_skin * ar_skin
    # Cumulative free skin IL-17A (% baseline x day) for the 12-week TE average
    d/dt(auc_free_target_skin) <- 100 * free_target_skin / r0_skin

    f(depot) <- exp(lfdepot)

    # Average free skin IL-17A over [0, t] (% baseline). The paper's TE metric is
    # this value at t = 84 days (12 weeks after the first dose); TE = 100 - avg.
    if (t > 0) {
      free_target_skin_avg_pct <- auc_free_target_skin / t
    } else {
      free_target_skin_avg_pct <- 100
    }
    te_skin_avg_pct <- 100 - free_target_skin_avg_pct

    # TE-based MBMA (Methods Eq. 1): placebo-adjusted responder fraction, defined
    # for the week-12 read-out (t = 84 days) over the plotted range 0.1-23% baseline
    xte <- free_target_skin_avg_pct
    ec50_pasi75 <- exp(lec50_pasi75)
    ec50_pasi90 <- exp(lec50_pasi90)
    prob_pasi75_pbo_adj <- (e0_pasi75 + xte^hill_pasi75 * (emax_pasi75 - e0_pasi75) / (xte^hill_pasi75 + ec50_pasi75^hill_pasi75)) / 100
    prob_pasi90_pbo_adj <- (e0_pasi90 + xte^hill_pasi90 * (emax_pasi90 - e0_pasi90) / (xte^hill_pasi90 + ec50_pasi90^hill_pasi90)) / 100

    # Observations (pM). Only serum total ixekizumab was fitted (Figure 4); the
    # skin and IL-17A outputs are model predictions (Figure 5).
    Cc <- cp
    Cis_skin <- cskin
    Ctotal_target <- total_target
    Ctotal_target_skin <- total_target_skin
    Cfree_target <- free_target
    Cfree_target_skin <- free_target_skin

    Cc ~ lnorm(expSd)
  })
}
