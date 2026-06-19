Li_2018_PF04236921 <- function() {
  description <- "Integrated population PK and indirect-response PK/PD model for the anti-interleukin-6 monoclonal antibody PF-04236921 in healthy volunteers and adults with rheumatoid arthritis, Crohn's disease, or systemic lupus erythematosus (Li 2018). Two-compartment IV/SC PK with first-order absorption and linear elimination from the central compartment; disease-stratified linear clearance and PD parameters; PF-04236921 inhibits the zero-order CRP synthesis rate of an indirect-response model."
  reference <- "Li C, Shoji S, Beebe J. Pharmacokinetics and C-reactive protein modelling of anti-interleukin-6 antibody (PF-04236921) in healthy volunteers and patients with autoimmune disease. Br J Clin Pharmacol. 2018;84(9):2059-2074. doi:10.1111/bcp.13641"
  vignette <- "Li_2018_PF04236921"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL", CRP = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, Vc, Q and Vp with reference 72 kg per the Li 2018 final-model equations (paper Eq. for CL/Vc/Q/Vp on p2064). The reference 72 kg is the pooled-cohort median across the five studies (Table 2, p2065 median values 79/82/61/70/73 kg for B0151001/B0151004/B0151002/B0151003/B0151006). Treated as baseline weight; the source paper does not describe time-varying weight handling.",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on linear CL with form theta_SEX_on_CL^SEXF (per Li 2018 paper equation on p2064 with SEX=1 male, SEX=2 female; the canonical SEXF coding 0=male/1=female is identical to (SEX_paper - 1)). The estimate 0.862 (Table 3, p2068) means females have 14 percent lower CL than males, consistent with the Discussion (p2070).",
      source_name        = "SEX (1 = male, 2 = female; canonical SEXF = SEX - 1)"
    ),
    ALB = list(
      description        = "Serum albumin concentration (baseline)",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 4.0 g/dL per the Li 2018 final-model equation on p2064 (alongside the BWT/ALB/CRCL/CRP/SEX covariate block). Also a power effect on baseline CRP and on the logit-scale Imax parameter with the same reference 4.0 g/dL (paper Eq. for BLCRP and Imax_prime on p2066). The reference 4.0 g/dL is the pooled-cohort median (Table 2 medians 4.5/4.5/4.2/4.0/3.9 g/dL for HV_IV/HV_SC/RA/CD/SLE).",
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw Cockcroft-Gault, not BSA-normalized)",
      units              = "mL/min (raw, not BSA-normalized)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 113 mL/min per the Li 2018 final-model equation on p2064. The paper Table 2 lists CLcr in 'ml min-1' (without the /1.73 m^2 BSA normalisation) so this model carries the source unit raw mL/min rather than the canonical mL/min/1.73 m^2 BSA-normalised form (see covariate-columns.md notes on per-model unit documentation). The reference 113 mL/min is approximately the pooled-cohort median (Table 2 medians 116/123/93.3/110/119 mL/min).",
      source_name        = "CLCR"
    ),
    CRP = list(
      description        = "C-reactive protein concentration (baseline)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 7.6 mg/L per the Li 2018 final-model equation on p2064. The assay is the high-sensitivity CRP assay used for the PD endpoint (Methods, p2062; LLOQ 0.01-0.02 mg/dL = 0.1-0.2 mg/L). The reference 7.6 mg/L is approximately the pooled-cohort median across the five studies (Table 2 medians 0.9/0.8/7.4/19.4/2.9 mg/L for HV_IV/HV_SC/RA/CD/SLE). CRP also enters the PD model as the state variable but as a covariate here it is the baseline pre-dose value.",
      source_name        = "CRP"
    ),
    DIS_RA = list(
      description        = "Rheumatoid arthritis cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HV / CD / SLE)",
      notes              = "One of three orthogonal binary indicators (DIS_RA, DIS_CD, DIS_SLE) decomposing the four-level Li 2018 disease-cohort categorical. Healthy volunteers (B0151001 IV + B0151004 SC) are the reference (all three indicators = 0). Disease-stratified parameters: typical CL (Table 3A theta_CL,HV/RA/CD/SLE = 0.00546/0.00588/0.00946/0.00643 L/h), typical baseline CRP, typical IC50, typical Imax (Table 3B). Coded here as log-shift covariate effects on the HV reference value.",
      source_name        = "POPULATION (paper categorical: HV / RA / CD / SLE; encoded as 0 for non-RA)"
    ),
    DIS_CD = list(
      description        = "Crohn's disease cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HV / RA / SLE)",
      notes              = "Paired with DIS_RA and DIS_SLE; HV reference when all three indicators = 0. CD-cohort typical CL is approximately 60 percent higher than the other cohorts (Discussion, p2070; Table 3A theta_CL,CD = 0.00946 L/h vs typical 0.0058 L/h elsewhere). Distinct from IBD_CD (which is CD vs UC in a pooled IBD cohort) -- here the complement is the heterogeneous HV/RA/SLE pool, not UC.",
      source_name        = "POPULATION (paper categorical: HV / RA / CD / SLE; encoded as 0 for non-CD)"
    ),
    DIS_SLE = list(
      description        = "Systemic lupus erythematosus cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HV / RA / CD)",
      notes              = "Paired with DIS_RA and DIS_CD; HV reference when all three indicators = 0. SLE-cohort drives a non-unit Hill coefficient gamma (1.55 vs 1 in the HV/RA/CD reference; Table 3B theta_gamma_SLE = 1.55) as well as separate typical Imax / IC50 / baseline CRP values.",
      source_name        = "POPULATION (paper categorical: HV / RA / CD / SLE; encoded as 0 for non-SLE)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 392L,
    n_studies      = 5L,
    age_range      = "18-73 years",
    age_median     = "33 years (HV IV) / 43 (HV SC) / 55 (RA) / 37 (CD) / 39 (SLE)",
    weight_range   = "30-157 kg",
    weight_median  = "72 kg (pooled across cohorts; per-cohort medians 79/82/61/70/73 kg)",
    sex_female_pct = 47.5,
    race_ethnicity = c(White = 76.5, Black = 12.8, Asian = 7.7, Other = 3.1),
    disease_state  = "Pooled: healthy volunteers (n = 36 IV in B0151001 + 10 SC in B0151004), rheumatoid arthritis on background methotrexate (n = 31 in B0151002), moderate-to-severe Crohn's disease refractory to anti-TNF (n = 178 in B0151003), active generalised systemic lupus erythematosus (n = 138 in B0151006).",
    dose_range     = "PF-04236921 7-700 mg IV single dose (B0151001 HV); 1, 10, 30, 100, 250 mg IV Q4W x 3 doses (B0151002 RA); 200 mg SC single dose (B0151004 HV); 10, 50, 200 mg SC at Day 1 and Week 4 (B0151003 CD); 10, 50, 200 mg SC at Day 1, Week 8 and Week 16 (B0151006 SLE).",
    regions        = "Multi-regional phase 1 / phase 2 studies (NCT00838565, NCT01287897, NCT01166555, NCT01405196).",
    notes          = paste(
      "Baseline demographics from Li 2018 Table 2 (p2065). PK assay LLOQ 100 ng/mL ELISA; CRP assay LLOQ 0.01-0.02 mg/dL (high-sensitivity particle-enhanced immunoturbidimetry).",
      "8.5 percent of post-dose PK samples were below quantitation and excluded; 9.4 percent of CRP samples were below the LLOQ and imputed as LLOQ/2.",
      "One RA patient with CWRES > 6 was excluded from PK refinement.",
      "Race (W/B/A/O) per Table 2: 422/50/30/25 across the five studies pooled (denominators include placebo arms).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural PK parameters - Li 2018 Table 3A final-model estimates
    # (p2068). Reference covariate values for the typical subject:
    #   BWT 72 kg, ALB 4.0 g/dL, CRCL 113 mL/min, CRP 7.6 mg/L, male (SEXF 0),
    #   disease cohort HV (DIS_RA = DIS_CD = DIS_SLE = 0).
    # Concentrations in ng/mL; doses in mg; central/peripheral state amounts
    # in mg; ka in 1/h; CL/Q in L/h; Vc/Vp in L.
    # =====================================================================
    lcl <- log(0.00546); label("Linear clearance CL (L/h) in the healthy-volunteer reference cohort")  # Li 2018 Table 3A theta_CL,HV
    lvc <- log(3.03);    label("Central volume of distribution Vc (L)")                                  # Li 2018 Table 3A theta_Vc
    lq  <- log(0.0245);  label("Intercompartmental clearance Q (L/h)")                                   # Li 2018 Table 3A theta_Q
    lvp <- log(3.58);    label("Peripheral volume of distribution Vp (L)")                               # Li 2018 Table 3A theta_Vp
    lka <- log(0.00607); label("First-order SC absorption rate ka (1/h)")                                # Li 2018 Table 3A theta_ka

    # SC bioavailability fixed to 1 (Results "Final PK model" p2063, p2071).
    lfdepot <- fixed(log(1)); label("SC absolute bioavailability F (fraction)")                          # Li 2018 Results p2063 "Bioavailability was estimated to be 100 percent and was fixed to 100 percent"

    # ---------------------------------------------------------------------
    # Disease-cohort covariate effects on CL (log-shift relative to HV).
    # Li 2018 Table 3A theta_CL by cohort: HV 0.00546, RA 0.00588, CD 0.00946,
    # SLE 0.00643 (all L/h). Encoded as log(theta_CL_cohort / theta_CL_HV) so
    # cl = exp(lcl + e_ra_cl*DIS_RA + e_cd_cl*DIS_CD + e_sle_cl*DIS_SLE + ...)
    # reproduces the per-cohort typical value.
    # ---------------------------------------------------------------------
    e_ra_cl  <- 0.07414; label("Log-shift on CL for RA cohort vs HV reference (unitless)")               # log(0.00588 / 0.00546); Li 2018 Table 3A
    e_cd_cl  <- 0.54988; label("Log-shift on CL for CD cohort vs HV reference (unitless)")               # log(0.00946 / 0.00546); Li 2018 Table 3A
    e_sle_cl <- 0.16323; label("Log-shift on CL for SLE cohort vs HV reference (unitless)")              # log(0.00643 / 0.00546); Li 2018 Table 3A

    # ---------------------------------------------------------------------
    # Continuous covariate effects on PK parameters (Li 2018 Table 3A).
    # All applied as power-of-ratio: (cov/ref)^exponent. References:
    # BWT 72 kg, ALB 4.0 g/dL, CRCL 113 mL/min, CRP 7.6 mg/L.
    # ---------------------------------------------------------------------
    e_wt_cl   <- 0.298;  label("Power exponent of WT/72 on linear CL (unitless)")                        # Li 2018 Table 3A theta_BWT_on_CL
    e_alb_cl  <- -1.53;  label("Power exponent of ALB/4.0 on linear CL (unitless)")                      # Li 2018 Table 3A theta_ALB_on_CL
    e_crcl_cl <- 0.287;  label("Power exponent of CRCL/113 on linear CL (unitless)")                     # Li 2018 Table 3A theta_CLcr_on_CL
    e_crp_cl  <- 0.0470; label("Power exponent of CRP/7.6 on linear CL (unitless)")                      # Li 2018 Table 3A theta_CRP_on_CL
    e_sex_cl  <- 0.862;  label("Multiplier on CL when SEXF = 1 (female-vs-male ratio, unitless)")        # Li 2018 Table 3A theta_SEX_on_CL; cl *= e_sex_cl^SEXF

    e_wt_vc <- 0.550; label("Power exponent of WT/72 on Vc (unitless)")                                  # Li 2018 Table 3A theta_BWT_on_Vc
    e_wt_q  <- 0.845; label("Power exponent of WT/72 on Q (unitless)")                                   # Li 2018 Table 3A theta_BWT_on_Q
    e_wt_vp <- 0.840; label("Power exponent of WT/72 on Vp (unitless)")                                  # Li 2018 Table 3A theta_BWT_on_Vp

    # Q has no independent eta; per Li 2018 Table 3A "theta_IIV_on_Q" and the
    # paper Eq. for Q on p2064, eta_Q,i = q_vc_scale * eta_Vc,i. q_vc_scale
    # is therefore a fixed-effect scalar that maps the Vc IIV onto Q.
    q_vc_scale <- 2.73; label("Scaling factor mapping eta_Vc onto eta_Q (unitless, paper Eq.)")          # Li 2018 Table 3A theta_IIV_on_Q

    # ---------------------------------------------------------------------
    # Inter-individual variability for PK (Li 2018 Table 3A CV% reported with
    # the explicit convention "%CV = sqrt(omega^2) x 100" per Table 3
    # footnote b on p2069). The underlying omega^2 that NONMEM estimated is
    # therefore (CV/100)^2 (not log(1 + CV^2)):
    #   CL  CV 36.5%  -> omega^2 = 0.365^2 = 0.13323
    #   Vc  CV 21.5%  -> omega^2 = 0.215^2 = 0.04623
    #   Vp  CV 38.6%  -> omega^2 = 0.386^2 = 0.14900
    #   ka  CV 21.2%  -> omega^2 = 0.212^2 = 0.04494
    # Off-diagonals from the reported correlations (Table 3A rho rows):
    #   rho(CL, Vc) = 0.710 -> cov = 0.710 * sqrt(0.13323 * 0.04623) = 0.05573
    #   rho(CL, Vp) = 0.606 -> cov = 0.606 * sqrt(0.13323 * 0.14900) = 0.08541
    #   rho(Vc, Vp) = 0.907 -> cov = 0.907 * sqrt(0.04623 * 0.14900) = 0.07531
    # ---------------------------------------------------------------------
    etalcl + etalvc + etalvp ~ c(0.13323,
                                  0.05573, 0.04623,
                                  0.08541, 0.07531, 0.14900)                                              # Li 2018 Table 3A IIV block
    etalka ~ 0.04494                                                                                      # Li 2018 Table 3A CV(eta_ka)

    # ---------------------------------------------------------------------
    # PK residual error. Li 2018 used a log-additive residual with study-
    # stratified variance (HV.IV 10.7 percent CV, HV.SC 28.3 percent CV,
    # pooled RA/CD/SLE 17.9 percent CV; Table 3A) and an additional IIV
    # (CV(eta_error) = 40.1 percent) multiplying the per-subject residual
    # SD. This model file collapses the three study-specific SDs to the
    # pooled patient-cohort value 17.9 percent CV and does not encode the
    # eta_error IIV on the residual SD; see vignette "Assumptions and
    # deviations" for the fidelity trade-off discussion.
    # ---------------------------------------------------------------------
    expSd <- 0.179; label("Log-normal residual SD on PF-04236921 concentration (fraction)")              # Li 2018 Table 3A pooled patient CV (RA/CD/SLE)

    # =====================================================================
    # PK/PD parameters - Li 2018 Table 3B final-model estimates (p2068).
    # Indirect-response model: PF-04236921 inhibits the zero-order CRP
    # synthesis rate.
    #   d/dt(effect) = Kin * (1 - Imax * Cc^gamma / (IC50^gamma + Cc^gamma))
    #                - Kout * effect
    #   effect(0)    = BLCRP = Kin / Kout
    # =====================================================================
    # Baseline CRP - typical for the HV reference cohort at the reference
    # ALB = 4.0 g/dL. Disease cohorts get log-shifts; ALB enters as a power
    # effect (Li 2018 paper Eq. for BLCRP on p2066).
    lbase <- log(1.16); label("Typical baseline CRP in HV reference cohort (mg/L)")                      # Li 2018 Table 3B theta_BLCRP,HV

    e_ra_base  <- 1.90651; label("Log-shift on baseline CRP for RA cohort vs HV reference (unitless)")   # log(7.81 / 1.16); Li 2018 Table 3B
    e_cd_base  <- 2.59226; label("Log-shift on baseline CRP for CD cohort vs HV reference (unitless)")   # log(15.5 / 1.16); Li 2018 Table 3B
    e_sle_base <- 0.73553; label("Log-shift on baseline CRP for SLE cohort vs HV reference (unitless)")  # log(2.42 / 1.16); Li 2018 Table 3B

    e_alb_base <- -3.69;   label("Power exponent of ALB/4.0 on baseline CRP (unitless)")                 # Li 2018 Table 3B theta_ALB_on_BLCRP

    # IC50 - typical for the HV reference cohort; no continuous covariates.
    # Disease cohorts get log-shifts (Li 2018 Table 3B).
    lic50 <- log(391); label("Typical IC50 in HV reference cohort (ng/mL)")                              # Li 2018 Table 3B IC50,HV = 391 ng/mL

    e_ra_ic50  <- -2.23494; label("Log-shift on IC50 for RA cohort vs HV reference (unitless)")          # log(41.8 / 391); Li 2018 Table 3B
    e_cd_ic50  <- -0.70110; label("Log-shift on IC50 for CD cohort vs HV reference (unitless)")          # log(194 / 391); Li 2018 Table 3B
    e_sle_ic50 <- -0.00770; label("Log-shift on IC50 for SLE cohort vs HV reference (unitless)")         # log(388 / 391); Li 2018 Table 3B

    # Imax - typical logit-scale value for HV reference at ALB = 4.0 g/dL.
    # Imax (fraction) = exp(theta_prime_Imax) / (1 + exp(theta_prime_Imax)).
    # Disease cohorts get log-shifts on the logit. Per Li 2018 paper Eq.
    # (p2066) ALB enters multiplicatively on the logit: theta_prime_Imax_i =
    # theta_Imax_cohort * (ALB_i / 4.0)^e_alb_logitimax. We encode the
    # multiplier exactly as published. logit(0.881) = log(0.881/0.119) =
    # 1.9974.
    logitimax <- log(0.881 / (1 - 0.881)); label("Logit-scale Imax in HV reference cohort at ALB 4.0 g/dL (unitless)")  # Li 2018 Table 3B Imax,HV = 88.1 percent

    # Log-ratio shifts on the LOGIT scale: e_X_logitimax =
    # log(logit(Imax_X) / logit(Imax_HV)), so the disease-cohort logit equals
    # logitimax * exp(e_X_logitimax) -- consistent with the paper's
    # multiplicative form theta_prime_Imax = theta_Imax * (ALB/4)^k.
    # logit values (computed at high precision):
    #   HV  log(0.881/0.119) = 2.00197
    #   RA  log(0.970/0.030) = 3.47610
    #   CD  log(0.989/0.011) = 4.49872
    #   SLE log(0.932/0.068) = 2.61785
    e_ra_logitimax  <- 0.55196;  label("Log-shift on logit-Imax for RA cohort vs HV (unitless)")          # log(3.47610 / 2.00197); Li 2018 Table 3B
    e_cd_logitimax  <- 0.80963;  label("Log-shift on logit-Imax for CD cohort vs HV (unitless)")          # log(4.49872 / 2.00197); Li 2018 Table 3B
    e_sle_logitimax <- 0.26824;  label("Log-shift on logit-Imax for SLE cohort vs HV (unitless)")         # log(2.61785 / 2.00197); Li 2018 Table 3B

    e_alb_logitimax <- -1.53; label("Power exponent of ALB/4.0 on logit-Imax (paper multiplicative form, unitless)")  # Li 2018 Table 3B theta_ALB_on_Imax

    # Hill coefficient. Fixed at 1 for HV/RA/CD per paper (Table 3B
    # "theta_gamma,HV,RA,CD = 1 fix"); estimated as 1.55 for SLE only.
    lgamma <- fixed(log(1));    label("Hill coefficient gamma in HV/RA/CD reference (fixed at 1)")        # Li 2018 Table 3B theta_gamma,HV,RA,CD = 1 fix
    e_sle_gamma <- log(1.55);   label("Log-shift on Hill gamma for SLE cohort vs reference (unitless)")   # Li 2018 Table 3B theta_gamma,SLE = 1.55

    # First-order CRP elimination rate (single value across all cohorts; the
    # paper found per-disease Kout estimation did not significantly improve
    # the model fit, dOFV = -3.142 with df = 3).
    lkout <- log(0.0238); label("First-order CRP elimination rate Kout (1/h)")                            # Li 2018 Table 3B theta_Kout

    # ---------------------------------------------------------------------
    # PD inter-individual variability (Li 2018 Table 3B; same convention as
    # the PK IIV: "%CV = sqrt(omega^2) x 100" per Table 3 footnote d on
    # p2069 so omega^2 = (CV/100)^2):
    #   BLCRP  CV 95.9% -> omega^2 = 0.959^2 = 0.91969
    #   IC50   CV 120%  -> omega^2 = 1.200^2 = 1.44000
    #   Imax   CV 111%  -> omega^2 = 1.110^2 = 1.23210 (on the logit scale,
    #                       since eta_Imax enters logit-Imax additively)
    # Block (BLCRP, Imax) correlation rho = 0.852 (Table 3B):
    #   cov(eta_base, eta_logitimax) = 0.852 * sqrt(0.91969 * 1.23210)
    #                                = 0.852 * 1.06504 = 0.90742
    # IC50 has independent eta.
    # ---------------------------------------------------------------------
    etalbase + etalogitimax ~ c(0.91969,
                                 0.90742, 1.23210)                                                        # Li 2018 Table 3B IIV block (BLCRP, Imax)
    etalic50 ~ 1.44000                                                                                    # Li 2018 Table 3B CV(eta_IC50)

    # PD residual: Li 2018 reports CV(eps) = 45.3 percent for the PD residual
    # plus CV(eta_error_PD) = 53.4 percent IIV on the residual SD. Same
    # simplification as the PK residual: encode a single log-normal SD
    # without the eta-on-residual layer.
    expSd_CRP_pred <- 0.453; label("Log-normal residual SD on CRP concentration (fraction)")                  # Li 2018 Table 3B residual CV
  })

  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # =====================================================================
    # Disease-cohort log-shift adjustments. Healthy-volunteer cohort is the
    # reference (all three DIS_* indicators = 0).
    # =====================================================================
    lcl_cohort   <- lcl + e_ra_cl  * DIS_RA + e_cd_cl  * DIS_CD + e_sle_cl  * DIS_SLE
    lbase_cohort <- lbase + e_ra_base * DIS_RA + e_cd_base * DIS_CD + e_sle_base * DIS_SLE
    lic50_cohort <- lic50 + e_ra_ic50 * DIS_RA + e_cd_ic50 * DIS_CD + e_sle_ic50 * DIS_SLE
    lgamma_cohort <- lgamma + e_sle_gamma * DIS_SLE

    # Logit-Imax cohort-typical value built multiplicatively per Li 2018
    # paper Eq. (p2066): theta_prime_Imax = theta_Imax * (alb_gdL/4)^theta_ALB.
    # Disease shifts enter as exponentials so logitimax_cohort scales the
    # HV reference logit by the cohort-specific ratio.
    logitimax_cohort <- logitimax *
      exp(e_ra_logitimax  * DIS_RA +
          e_cd_logitimax  * DIS_CD +
          e_sle_logitimax * DIS_SLE)

    # =====================================================================
    # Individual PK parameters (Li 2018 paper Eq. for CL/Vc/Q/Vp on p2064).
    # Reference covariate values for the typical subject already documented
    # in the ini() block comments.
    # =====================================================================
    cl <- exp(lcl_cohort + etalcl) *
          (WT / 72)^e_wt_cl *
          (alb_gdL / 4.0)^e_alb_cl *
          (CRCL / 113)^e_crcl_cl *
          (CRP / 7.6)^e_crp_cl *
          (e_sex_cl)^SEXF
    vc <- exp(lvc + etalvc) * (WT / 72)^e_wt_vc
    q  <- exp(lq + q_vc_scale * etalvc) * (WT / 72)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 72)^e_wt_vp
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # =====================================================================
    # 2-compartment PK with first-order SC absorption + linear elimination
    # from central. IV doses can be administered directly to the central
    # compartment via the cmt column on the event record; SC doses go to the
    # depot compartment.
    # =====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Cc in ng/mL: doses are in mg, central state amount is in mg, vc in L,
    # so central/vc is in mg/L. 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # =====================================================================
    # PD: indirect response on CRP. PF-04236921 inhibits the synthesis rate.
    #   d/dt(effect) = Kin * (1 - u(t)) - Kout * effect
    #   u(t)         = Imax * Cc^gamma / (IC50^gamma + Cc^gamma)
    #   effect(0)    = BLCRP, Kin = Kout * BLCRP (steady-state pre-dose)
    # The state `effect` carries CRP in mg/L (paper Eq. on p2062-2066).
    # =====================================================================
    base <- exp(lbase_cohort + etalbase) * (alb_gdL / 4.0)^e_alb_base
    logitimax_prime <- logitimax_cohort * (alb_gdL / 4.0)^e_alb_logitimax
    imax <- exp(logitimax_prime + etalogitimax) /
            (1 + exp(logitimax_prime + etalogitimax))
    ic50 <- exp(lic50_cohort + etalic50)
    gamma <- exp(lgamma_cohort)
    kout <- exp(lkout)
    kin <- kout * base

    inh <- imax * Cc^gamma / (ic50^gamma + Cc^gamma)

    d/dt(effect) <- kin * (1 - inh) - kout * effect
    effect(0)    <- base

    CRP_pred <- effect

    Cc       ~ lnorm(expSd)
    CRP_pred ~ lnorm(expSd_CRP_pred)
  })
}
