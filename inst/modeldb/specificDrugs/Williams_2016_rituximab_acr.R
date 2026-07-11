Williams_2016_rituximab_acr <- function() {
  description <- paste(
    "Longitudinal ordered-categorical PD model of American College of",
    "Rheumatology (ACR) response (ACR20 / ACR50 / ACR70) over time for",
    "rituximab (both the proposed biosimilar PF-05280586 and the two rituximab",
    "reference products sourced from the EU and US) in adults with rheumatoid",
    "arthritis on background methotrexate with prior inadequate response to",
    "one or more TNF-antagonist therapies (Williams 2016). The model is a",
    "cumulative-probit latent-variable formulation with three thresholds",
    "(ACR20, ACR50, ACR70 as an ordered categorical score in {0, 1, 2, 3})",
    "and an exponential-onset time course for the population-mean latent",
    "variable. Unlike the companion DAS28cfb model, no rituximab exposure",
    "effect is included (Williams 2016 Results: 'rituximab exposure effects",
    "could not be supported in the ACR model'). The covariate model includes",
    "seven baseline disease-activity components (TJ68, SJ66, log CRP+1, PGA_PT,",
    "PhGA, PAIN, HAQ-DI) and two treatment-arm indicators (PF-05280586 vs",
    "rituximab-EU, rituximab-US vs rituximab-EU) additive on the log scale on",
    "both the maximum-effect parameter PMAX and the onset half-life parameter,",
    "with rituximab-EU as the reference arm. Outputs three continuous",
    "probabilities pACR20, pACR50, pACR70 = probability of achieving at least",
    "ACR20/50/70 response at the observation time; a single subject-level",
    "additive probit-scale eta captures between-subject variability.",
    sep = " "
  )

  reference <- paste(
    "Williams JH, Hutmacher MM, Zierhut ML, Becker JC, Gumbiner B,",
    "Spencer-Green G, Melia LA, Liao KH, Suster M, Yin D, Li R, Meng X.",
    "Comparative assessment of clinical response in patients with rheumatoid",
    "arthritis between PF-05280586, a proposed rituximab biosimilar, and",
    "rituximab. Br J Clin Pharmacol. 2016;82(6):1568-1579.",
    "doi:10.1111/bcp.13094. PMID: 27530379.",
    "Companion DAS28cfb model from the same paper:",
    "modellib('Williams_2016_rituximab_das28cfb').",
    "Latent-variable ordered-categorical framework adapted from",
    "Hutmacher MM, Krishnaswami S, Kowalski KG.",
    "Exposure-response modeling using latent variables for the efficacy of",
    "a JAK3 inhibitor administered to rheumatoid arthritis patients.",
    "J Pharmacokinet Pharmacodyn 2008;35:139-157.",
    "doi:10.1007/s10928-007-9080-2.",
    sep = " "
  )

  vignette <- "Williams_2016_rituximab_biosimilar"

  paper_specific_etas <- c("etaACR")
  paper_specific_residual_sds <- c(
    "addSd_pacr20", "addSd_pacr50", "addSd_pacr70"
  )

  units <- list(
    time          = "week",
    dosing        = "(none; PD-only ordered-categorical responder-rate model with no PK compartment or exposure driver)",
    concentration = "(observation is per-time probability of achieving ACR20/ACR50/ACR70 on the 0-1 scale; latent variable is on the probit scale)"
  )

  covariateData <- list(
    TEND_68JOINT = list(
      description        = "Baseline tender joint count on the 68-joint (extended) scale (integer 0-68).",
      units              = "count (0-68)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the covariate model as an additive",
        "log-scale effect (TEND_68JOINT - 24) on PMAX and on the onset",
        "half-life (Williams 2016 Supplemental Methods; reference value 24 =",
        "paper's declared median of the ACR dataset). Distinct from",
        "TEND_28JOINT which uses the 28-joint DAS28 subscale.",
        sep = " "
      ),
      source_name        = "TJ68"
    ),
    SWOL_66JOINT = list(
      description        = "Baseline swollen joint count on the 66-joint (extended) scale (integer 0-66).",
      units              = "count (0-66)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the covariate model as an additive",
        "log-scale effect (SWOL_66JOINT - 16) on PMAX and onset half-life",
        "(Williams 2016 Supplemental Methods; reference value 16 = paper's",
        "declared median of the ACR dataset). Distinct from SWOL_28JOINT",
        "which uses the 28-joint DAS28 subscale.",
        sep = " "
      ),
      source_name        = "SJ66"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein (BCRP).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline value). Enters the covariate model",
        "on the log(CRP + 1) scale to accommodate the highly skewed",
        "distribution and BCRP = 0 observations (Williams 2016 Supplemental",
        "Methods): lBCRP_i = log(CRP + 1). The reference value on that",
        "log-transformed scale is 2.2 log-units. Applied as (log(CRP + 1)",
        "- 2.2) additively on PMAX and onset half-life.",
        sep = " "
      ),
      source_name        = "BCRP"
    ),
    PGA_PT = list(
      description        = "Baseline patient's global assessment of arthritis (100-mm visual analogue scale).",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Additive log-scale effect (PGA_PT - 70)",
        "on PMAX and onset half-life; reference 70 mm is the paper's",
        "declared median of the ACR dataset. Distinct from the physician's",
        "global assessment (BLPHYVAS) and from patient-reported pain (PAIN).",
        sep = " "
      ),
      source_name        = "PGA"
    ),
    BLPHYVAS = list(
      description        = "Baseline physician's global assessment of disease activity (100-mm visual analogue scale).",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Additive log-scale effect (BLPHYVAS - 68) on",
        "PMAX and onset half-life (Williams 2016 Supplemental Methods;",
        "reference 68 = paper's declared median of the ACR dataset).",
        "Distinct from the patient's own global assessment (PGA_PT).",
        sep = " "
      ),
      source_name        = "PhGA"
    ),
    PAIN = list(
      description        = "Baseline patient-reported global arthritis pain (100-mm visual analogue scale).",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Additive log-scale effect (PAIN - 70) on",
        "PMAX and onset half-life; reference 70 = paper's declared median.",
        "The PAIN VAS is the patient's own pain rating; distinct from PGA_PT",
        "(patient's global disease-activity assessment) and BLPHYVAS",
        "(physician's global assessment).",
        sep = " "
      ),
      source_name        = "PAIN"
    ),
    BLHAQ = list(
      description        = "Baseline Health Assessment Questionnaire Disability Index (HAQ-DI; 0-3 composite score).",
      units              = "unitless (0-3 composite)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Additive log-scale effect (BLHAQ - 1.75) on",
        "PMAX and onset half-life; reference 1.75 = paper's declared median",
        "of the ACR dataset.",
        sep = " "
      ),
      source_name        = "HAQ-DI"
    ),
    TRT = list(
      description        = "Treatment-arm integer indicator: 0 = rituximab-EU (reference), 1 = PF-05280586 (proposed biosimilar), 2 = rituximab-US.",
      units              = "(categorical / integer-coded)",
      type               = "categorical",
      reference_category = "0 (rituximab-EU)",
      notes              = paste(
        "Time-fixed per subject. Same integer coding as the companion",
        "DAS28cfb model file. Encodes the trial's three treatment arms.",
        "Two derived indicators are computed inside model() as",
        "TRT_PF = 1*(TRT == 1) and TRT_US = 1*(TRT == 2), each entering",
        "the covariate function additively on the log scale.",
        sep = " "
      ),
      source_name        = "TRT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 214L,
    n_studies        = 1L,
    n_observations   = 1402L,
    age_range        = "adults (>=18 years); mean (SD) 54.8 (11.7) years (PF-05280586), 55.7 (10.2) years (rituximab-EU), 53.8 (11.8) years (rituximab-US)",
    weight_range     = "mean (SD) 86.2 (22.0) kg (PF-05280586), 82.6 (19.8) kg (rituximab-EU), 80.4 (21.6) kg (rituximab-US)",
    sex_female_pct   = 77.6,
    disease_state    = "Active rheumatoid arthritis on background methotrexate with inadequate response to one or more TNF-antagonist therapies (see companion DAS28cfb model for full baseline demographics).",
    dose_range       = "1000 mg IV on days 1 and 15 (standard rituximab RA induction course). All subjects received 100 mg IV methylprednisolone premedication.",
    regions          = "Multi-regional biosimilar development trial (ClinicalTrials.gov NCT01526057).",
    baseline_disease = list(
      TEND_68JOINT_by_arm = "PF-05280586 22.7 (12.6); rituximab-EU 23.7 (13.2); rituximab-US 29.7 (15.0)",
      SWOL_66JOINT_by_arm = "PF-05280586 15.4 (8.8); rituximab-EU 17.9 (10.6); rituximab-US 18.9 (8.4)",
      CRP_by_arm_mg_L     = "PF-05280586 12.4 (14.9); rituximab-EU 14.7 (17.6); rituximab-US 18.2 (25.1)",
      PGA_by_arm          = "PF-05280586 67.4 (16.8); rituximab-EU 67.7 (20.9); rituximab-US 74.8 (16.0)",
      PhGA_by_arm         = "PF-05280586 64.6 (15.3); rituximab-EU 66.1 (15.5); rituximab-US 70.1 (15.6)",
      PAIN_by_arm         = "PF-05280586 65.6 (17.8); rituximab-EU 66.1 (21.0); rituximab-US 72.1 (18.5)",
      HAQ_DI_by_arm       = "PF-05280586 1.67 (0.56); rituximab-EU 1.61 (0.53); rituximab-US 1.74 (0.62)"
    ),
    notes            = paste(
      "Baseline demographics from Williams 2016 Table 1. The ACR responder",
      "dataset included 1402 observations",
      "(Williams 2016 Results 'Population PK/PD models'). Composite",
      "responder-rate outcomes (ACR20, ACR50, ACR70) were simultaneously",
      "modeled via the paper's cumulative-probit latent-variable formulation",
      "(Williams 2016 Supp Methods 'ACR Responder Model').",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------------
    # Cumulative-probit thresholds (Williams 2016 Supp Table S1, With Baseline
    # Covariates column). The threshold for achieving at least ACR20 is
    # probitACR20 = theta_1; the additional gap to ACR50 is exp(lgap_acr50);
    # the additional gap to ACR70 is exp(lgap_acr70). Sum of exp(gap) values
    # gives the total offset from ACR20 to ACR70 on the probit scale.
    # Transformed estimates given in comments for reader convenience.
    # --------------------------------------------------------------------------
    probitACR20 <- -1.91;    label("Probit-scale intercept for the ACR20 (score >= 1) threshold")                # Williams 2016 Supp Table S1: theta_1 estimate -1.91 (SE 0.529)
    lgap_acr50  <-  0.306;   label("Log gap between the ACR20 and ACR50 thresholds on the probit scale (exp(0.306) = 1.36)") # Williams 2016 Supp Table S1: theta_2 estimate 0.306 (SE 0.0567)
    lgap_acr70  <- -0.196;   label("Log gap between the ACR50 and ACR70 thresholds on the probit scale (exp(-0.196) = 0.822)") # Williams 2016 Supp Table S1: theta_3 estimate -0.196 (SE 0.0822)

    # --------------------------------------------------------------------------
    # Latent-variable time-course parameters (Williams 2016 Supp Methods
    # ACR Responder Model; Supp Table S1 With Baseline Covariates column).
    # Both are on log-scale by design so that covariate effects and etas
    # enter additively.
    # --------------------------------------------------------------------------
    lpmax_acr     <- -2.73;   label("Typical maximum latent-variable effect for the ACR responder-rate exponential-onset time course (probit units)")  # Williams 2016 Supp Table S1: PMAX estimate -2.73 (SE 0.506)
    lthalfrec_acr <-  0.958;  label("Log typical onset half-life of the ACR latent-variable time course, giving kp = log(2)/exp(lthalfrec_acr) (weeks)") # Williams 2016 Supp Table S1: Onset (kp) estimate 0.958 (SE 0.442)

    # --------------------------------------------------------------------------
    # Covariate effects on PMAX (Williams 2016 Supp Table S1; additive on
    # log-probit scale of the maximum latent-variable effect).
    # --------------------------------------------------------------------------
    e_tj68_pmax_acr    <-  0.0422;   label("Log-scale additive effect of TEND_68JOINT - 24 on PMAX (per count)")  # Williams 2016 Supp Table S1: TJ68 on PMAX 0.0422 (SE 0.0122)
    e_sj66_pmax_acr    <- -0.0401;   label("Log-scale additive effect of SWOL_66JOINT - 16 on PMAX (per count)")   # Williams 2016 Supp Table S1: SJ66 on PMAX -0.0401 (SE 0.0190)
    e_lbcrp_pmax_acr   <-  0.0641;   label("Log-scale additive effect of (log(CRP+1) - 2.2) on PMAX (per log unit)") # Williams 2016 Supp Table S1: lBCRP on PMAX 0.0641 (SE 0.144)
    e_pga_pmax_acr     <- -0.00833;  label("Log-scale additive effect of PGA_PT - 70 on PMAX (per mm VAS)")        # Williams 2016 Supp Table S1: PGA on PMAX -0.00833 (SE 0.00925)
    e_phga_pmax_acr    <-  0.00434;  label("Log-scale additive effect of BLPHYVAS - 68 on PMAX (per mm VAS)")      # Williams 2016 Supp Table S1: PhGA on PMAX 0.00434 (SE 0.00746)
    e_pain_pmax_acr    <- -0.0158;   label("Log-scale additive effect of PAIN - 70 on PMAX (per mm VAS)")           # Williams 2016 Supp Table S1: PAIN on PMAX -0.0158 (SE 0.00983)
    e_haqdi_pmax_acr   <-  0.0418;   label("Log-scale additive effect of BLHAQ - 1.75 on PMAX (per HAQ-DI unit)")   # Williams 2016 Supp Table S1: HAQ-DI on PMAX 0.0418 (SE 0.222)
    e_trt_pf_pmax_acr  <-  0.244;    label("Log-scale additive effect on PMAX for TRT = 1 (PF-05280586)")          # Williams 2016 Supp Table S1: PF-05280586 on PMAX 0.244 (SE 0.284)
    e_trt_us_pmax_acr  <- -0.293;    label("Log-scale additive effect on PMAX for TRT = 2 (rituximab-US)")         # Williams 2016 Supp Table S1: Ritux-US on PMAX -0.293 (SE 0.274)

    # --------------------------------------------------------------------------
    # Covariate effects on the onset half-life (Williams 2016 Supp Table S1;
    # positive value means longer half-life = slower onset).
    # --------------------------------------------------------------------------
    e_tj68_thalfrec_acr    <-  0.0252;   label("Log-scale additive effect of TEND_68JOINT - 24 on lthalfrec_acr")   # Williams 2016 Supp Table S1: TJ68 on kp 0.0252 (SE 0.0158)
    e_sj66_thalfrec_acr    <-  0.00292;  label("Log-scale additive effect of SWOL_66JOINT - 16 on lthalfrec_acr")    # Williams 2016 Supp Table S1: SJ66 on kp 0.00292 (SE 0.0176)
    e_lbcrp_thalfrec_acr   <-  0.0676;   label("Log-scale additive effect of (log(CRP+1) - 2.2) on lthalfrec_acr")   # Williams 2016 Supp Table S1: lBCRP on kp 0.0676 (SE 0.131)
    e_pga_thalfrec_acr     <-  0.000573; label("Log-scale additive effect of PGA_PT - 70 on lthalfrec_acr")          # Williams 2016 Supp Table S1: PGA on kp 0.000573 (SE 0.0153)
    e_phga_thalfrec_acr    <-  0.00801;  label("Log-scale additive effect of BLPHYVAS - 68 on lthalfrec_acr")        # Williams 2016 Supp Table S1: PhGA on kp 0.00801 (SE 0.00664)
    e_pain_thalfrec_acr    <-  0.00745;  label("Log-scale additive effect of PAIN - 70 on lthalfrec_acr")             # Williams 2016 Supp Table S1: PAIN on kp 0.00745 (SE 0.0147)
    e_haqdi_thalfrec_acr   <- -0.413;    label("Log-scale additive effect of BLHAQ - 1.75 on lthalfrec_acr")          # Williams 2016 Supp Table S1: HAQ-DI on kp -0.413 (SE 0.255)
    e_trt_pf_thalfrec_acr  <- -0.368;    label("Log-scale additive effect on lthalfrec_acr for TRT = 1 (PF-05280586)") # Williams 2016 Supp Table S1: PF-05280586 on kp -0.368 (SE 0.361)
    e_trt_us_thalfrec_acr  <- -0.137;    label("Log-scale additive effect on lthalfrec_acr for TRT = 2 (rituximab-US)") # Williams 2016 Supp Table S1: Ritux-US on kp -0.137 (SE 0.346)

    # --------------------------------------------------------------------------
    # Subject-level additive random effect on the probit scale (Williams 2016
    # Supp Table S1: SD(eta1) log-scale estimate 0.192 -> transformed 1.21
    # probit units, so variance = 1.21^2 = 1.4641). Paper-specific eta (not
    # tied to a log-transformed structural parameter).
    # --------------------------------------------------------------------------
    etaACR ~ 1.4641

    # --------------------------------------------------------------------------
    # Placeholder additive residual errors on the three cumulative probability
    # outputs. The source paper's likelihood is Bernoulli (categorical) with
    # Laplace estimation for the ACR responder categorical response; nlmixr2's
    # residual-error grammar does not currently expose a categorical /
    # Bernoulli residual for a probability output declared in ini() / model().
    # These small additive residual SDs are structural placeholders for
    # simulation-only use; downstream fitting workflows should replace them
    # with the appropriate categorical likelihood.
    # --------------------------------------------------------------------------
    addSd_pacr20 <- 0.01;  label("Placeholder additive residual SD on pACR20 (probability units; see model comment)") # placeholder for Bernoulli likelihood
    addSd_pacr50 <- 0.01;  label("Placeholder additive residual SD on pACR50 (probability units; see model comment)") # placeholder for Bernoulli likelihood
    addSd_pacr70 <- 0.01;  label("Placeholder additive residual SD on pACR70 (probability units; see model comment)") # placeholder for Bernoulli likelihood
  })

  model({
    # Treatment-arm indicators (per Williams 2016 Supp Methods 'ACR Responder
    # Model'; g_i^(A) covariate function includes I(TRT_i = PF-05280586) and
    # I(TRT_i = Ritux-US) with rituximab-EU as the reference arm).
    TRT_PF <- (TRT == 1)
    TRT_US <- (TRT == 2)

    # Log-transformed baseline CRP (Williams 2016 Supp Methods:
    # lBCRP = log(BCRP + 1); reference median 2.2 log-units).
    lBCRP <- log(CRP + 1)

    # Individual log-scale structural parameters with covariate effects
    # centred on the paper's declared reference values (TJ68 = 24, SJ66 = 16,
    # lBCRP = 2.2, PGA = 70, PhGA = 68, PAIN = 70, HAQ-DI = 1.75,
    # TRT = 0 = rituximab-EU).
    pmax_i <- lpmax_acr +
              e_tj68_pmax_acr    * (TEND_68JOINT - 24) +
              e_sj66_pmax_acr    * (SWOL_66JOINT - 16) +
              e_lbcrp_pmax_acr   * (lBCRP - 2.2) +
              e_pga_pmax_acr     * (PGA_PT - 70) +
              e_phga_pmax_acr    * (BLPHYVAS - 68) +
              e_pain_pmax_acr    * (PAIN - 70) +
              e_haqdi_pmax_acr   * (BLHAQ - 1.75) +
              e_trt_pf_pmax_acr  * TRT_PF +
              e_trt_us_pmax_acr  * TRT_US

    lthalfrec_i <- lthalfrec_acr +
                    e_tj68_thalfrec_acr    * (TEND_68JOINT - 24) +
                    e_sj66_thalfrec_acr    * (SWOL_66JOINT - 16) +
                    e_lbcrp_thalfrec_acr   * (lBCRP - 2.2) +
                    e_pga_thalfrec_acr     * (PGA_PT - 70) +
                    e_phga_thalfrec_acr    * (BLPHYVAS - 68) +
                    e_pain_thalfrec_acr    * (PAIN - 70) +
                    e_haqdi_thalfrec_acr   * (BLHAQ - 1.75) +
                    e_trt_pf_thalfrec_acr  * TRT_PF +
                    e_trt_us_thalfrec_acr  * TRT_US
    kp_i <- log(2) / exp(lthalfrec_i)

    # Model equations (Williams 2016 Supp Methods 'ACR Responder Model'):
    #   fnon-C(t) = (theta_PMAX + g_i^(PMAX)(theta)) * (1 - exp(-kp_i * t))
    #   Phi^-1(P[score >= m | eta]) = theta_1 - exp(theta_2) * I_ACR50
    #                                 - exp(theta_3) * I_ACR70 + fnon-C(t) + eta_i
    # Cumulative thresholds on the probit scale (score >= 1 for ACR20, etc.):
    #   P[score >= 1] = P(achieve at least ACR20) uses only theta_1
    #   P[score >= 2] = P(achieve at least ACR50) subtracts exp(theta_2)
    #   P[score >= 3] = P(achieve ACR70)          subtracts exp(theta_2) + exp(theta_3)
    fnonC <- pmax_i * (1 - exp(-kp_i * time))
    lv    <- fnonC + etaACR

    pACR20 <- pnorm(probitACR20 + lv)
    pACR50 <- pnorm(probitACR20 - exp(lgap_acr50) + lv)
    pACR70 <- pnorm(probitACR20 - exp(lgap_acr50) - exp(lgap_acr70) + lv)

    pACR20 ~ add(addSd_pacr20)
    pACR50 ~ add(addSd_pacr50)
    pACR70 ~ add(addSd_pacr70)
  })
}
