Frey_2013_tocilizumab <- function() {
  description <- "Indirect-response PK/PD model of tocilizumab on the 28-joint Disease Activity Score (DAS28) in adults with rheumatoid arthritis (Levi/Grange/Frey 2013, OPTION + TOWARD phase III pool, n = 1703 patients with 12,618 DAS28 observations). Tocilizumab inhibits the DAS28 production rate kin via a sigmoid emax function whose driving concentration is the sum of circulating tocilizumab and a constant DMARD background term expressed in tocilizumab concentration units. The PK driver is the two-compartment, parallel linear + Michaelis-Menten model of Frey 2010 (PMID 20097931), reused unchanged for the exposure-response analysis."
  reference <- "Levi M, Grange S, Frey N. Exposure-response relationship of tocilizumab, an anti-IL-6 receptor monoclonal antibody, in a large population of patients with rheumatoid arthritis. J Clin Pharmacol. 2013;53(2):151-159. doi:10.1177/0091270012437585. PMID 23436260. PK backbone from Frey N, Grange S, Woodworth T. Population pharmacokinetic analysis of tocilizumab in patients with rheumatoid arthritis. J Clin Pharmacol. 2010;50(7):754-766. doi:10.1177/0091270009350623."
  vignette <- "Frey_2013_tocilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", response = "DAS28 (unitless 0-10 score)")

  covariateData <- list(
    IL6 = list(
      description        = "Baseline serum interleukin-6 (IL-6) concentration",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Frey 2013 enters IL-6 on the natural-log scale through the dimensionless ratio (log(IL-6 * 1000) / 9.9), where 9.9 = log(20000) corresponds to a reference IL-6 of ~20 pg/mL (the OPTION/TOWARD median of 19.9-22 pg/mL per Supplementary Table S1). The same log-IL-6 ratio enters three different parameters in the final model (Table 2): ec50 with exponent -4.4, BASE with exponent +0.13, and the DMARD background-effect parameter with exponent -6.4. The canonical column carries the raw IL-6 in pg/mL; the log transform is applied inside model().",
      source_name        = "IL-6"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Frey 2013 Table 2 reports the sex effect on emax with female as the reference (Emax_female = emax * 1.0; Emax_male = emax * 1.1, +11% in males). The canonical SEXF column is 1 = female, 0 = male, so the model applies the equation emax = Emax_typ * (1 + 0.11 * (1 - SEXF)) which preserves the female-as-reference NONMEM parameterization. The +11% male offset is below the known DAS28 measurement error (0.6 units) and is not clinically significant per the paper's Discussion.",
      source_name        = "SEX"
    ),
    RACE_ASIAN_AMIND_OTH = list(
      description        = "Composite race indicator: 1 = Asian, American Indian / Alaska Native, or Other; 0 = White or Black",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Black)",
      notes              = "Frey 2013 pools the smaller-N race groups (Asian, American Indian/Alaska Native, Other) into a single Asian-and-others composite and uses White + Black as the reference. The composite covers ~21-27% of the OPTION/TOWARD pool per Supplementary Table S1. Multiplicative effect on kout: kout = Kout_typ * (1 + (-0.25) * RACE_ASIAN_AMIND_OTH), i.e., kout is 25% lower in the Asian/AmInd/Other composite than in the White+Black reference. Paper-defined composite grouping; see the canonical RACE_ASIAN_AMIND_OTH register entry for the rationale.",
      source_name        = "RACE"
    ),
    BLHAQ = list(
      description        = "Baseline Health Assessment Questionnaire Disability Index (HAQ-DI; 0-3 score)",
      units              = "unitless (0-3 composite)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on the indirect-response BASE parameter: BASE = 6.8 * (BLHAQ/1.6)^0.043 (Frey 2013 Table 2; reference 1.6 is approximately the OPTION/TOWARD pooled median per Supplementary Table S1). Frey 2013 reports a paper-side HAQ floor of 0.010 in Table 2's covariate-range column to keep the power form well-defined when HAQ = 0; the model file applies the same floor inside model().",
      source_name        = "HAQ"
    ),
    PAIN = list(
      description        = "Baseline patient-reported global pain on a 100-mm visual analogue scale",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on BASE: BASE = 6.8 * (PAIN/60)^0.062 (Frey 2013 Table 2; reference 60 is the OPTION/TOWARD pooled median per Supplementary Table S1). Distinct from BLPHYVAS (the *physician*'s global VAS). Frey 2013 reports a paper-side PAIN floor of 0.010 to keep the power form well-defined when PAIN = 0; the model file applies the same floor inside model().",
      source_name        = "PAIN"
    ),
    BLPHYVAS = list(
      description        = "Baseline Physician's Global Assessment of Disease Activity (100-mm visual analogue scale)",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on BASE: BASE = 6.8 * (BLPHYVAS/65)^0.13 (Frey 2013 Table 2; reference 65 is the OPTION/TOWARD pooled median per Supplementary Table S1).",
      source_name        = "VASP"
    )
  )

  population <- list(
    n_subjects     = 1703L,
    n_observations = 12618L,
    n_studies      = 2L,
    age_range      = "18-89 years; OPTION median 52, TOWARD median 54",
    weight_range   = "36-148 kg; OPTION median 68, TOWARD median 70.9",
    sex_female_pct = 82,
    race_ethnicity = c(White = 73.5, "American Indian or Alaska Native" = 9.0, Asian = 10.4, Other = 3.9, Black = 3.6),
    disease_state  = "Moderate-to-severe rheumatoid arthritis (adults). OPTION: methotrexate-inadequate responders. TOWARD: traditional-DMARD-inadequate responders.",
    dose_range     = "Tocilizumab 4 or 8 mg/kg by 1-hour IV infusion every 4 weeks for up to 24 weeks, plus a placebo-on-DMARD-background arm. OPTION patients additionally received methotrexate; TOWARD patients additionally received methotrexate or other traditional DMARDs.",
    regions        = "International multi-regional (OPTION + TOWARD phase III studies pooled).",
    baseline_biomarkers = list(
      IL6_median_pg_mL_OPTION = 22,
      IL6_median_pg_mL_TOWARD = 19.9,
      sIL6R_median_ng_mL_OPTION = 36,
      sIL6R_median_ng_mL_TOWARD = 43.1,
      BLHAQ_median_OPTION  = 1.6,
      BLHAQ_median_TOWARD  = 1.5,
      PAIN_median_OPTION   = 62,
      PAIN_median_TOWARD   = 60,
      BLPHYVAS_median_OPTION = 65,
      BLPHYVAS_median_TOWARD = 65,
      RheumatoidFactor_median_U_mL_OPTION = 79,
      RheumatoidFactor_median_U_mL_TOWARD = 112
    ),
    notes          = "Baseline demographics from Frey 2013 Supplementary Table S1 (OPTION + TOWARD; OPTION N = 572 in the demographics table per the male+female totals of 105 + 467, TOWARD N = 1131 per 200 + 931 -- the paper's PKPD analysis subset n = 1703 reflects the union after restricting to subjects with at least one DAS28 observation and available individual PK estimates). Approximately 80% female and predominantly White; the smaller race groups (Asian, American Indian/Alaska Native, Other) are pooled by the paper into the Asian-and-others composite that drives the RACE_ASIAN_AMIND_OTH covariate. The PD analysis used DAS28-ESR (not DAS28-CRP); the residual error of 0.68 DAS28 units is consistent with the published DAS28 measurement error of 0.60."
  )

  ini({
    # ------------------------------------------------------------------
    # PK backbone (Frey 2010 Table II final-model estimates) - reused
    # unchanged in Frey 2013 (the 2013 paper does not refit PK; it uses
    # Frey 2010's individual empirical Bayes PK to drive the PD model).
    # The structural PK constants below correspond to Frey 2010's typical
    # subject (BSA = 1.8 m^2, HDL-C = 54 mg/dL, log(RF) = 4.7,
    # total protein = 74 g/L, albumin = 38 g/L, creatinine clearance =
    # 106 mL/min, non-smoker, male). Frey 2010 covariate effects on PK
    # (BSA / SEX / HDL-C / log-RF on CL; total protein / albumin on Vc;
    # albumin / creatinine clearance / smoking on Vm) are intentionally
    # omitted in this exposure-response file because they are not part
    # of the Frey 2013 PD model; see the vignette's Assumptions and
    # deviations section, and the standalone Frey_2010_tocilizumab.R
    # model for the full PK structure.
    # ------------------------------------------------------------------
    lcl <- log(0.3); label("Linear clearance CL (L/day) -- Frey 2010 typical PK reference")             # Frey 2010 Table II, CL
    lvc <- log(3.5); label("Central volume of distribution V1 (L) -- Frey 2010 typical PK reference")   # Frey 2010 Table II, V1
    lq  <- log(0.2); label("Inter-compartmental clearance Q (L/day) -- Frey 2010 typical PK reference") # Frey 2010 Table II, Q
    lvp <- log(2.9); label("Peripheral volume of distribution V2 (L) -- Frey 2010 typical PK reference") # Frey 2010 Table II, V2
    lvmax <- log(7.5); label("Maximum Michaelis-Menten elimination rate Vmax (mg/day) -- Frey 2010 typical PK reference") # Frey 2010 Table II, VM
    lkm <- log(2.7); label("Michaelis-Menten constant Km (ug/mL) -- Frey 2010 typical PK reference")    # Frey 2010 Table II, KM

    # ------------------------------------------------------------------
    # DAS28 indirect-response PD model -- Frey 2013 Table 1 typical values.
    # Reference covariate values for the typical subject are the OPTION /
    # TOWARD pooled medians: IL6 = 20 pg/mL (i.e., log(IL6 * 1000) = 9.9),
    # BLHAQ = 1.6, PAIN = 60, BLPHYVAS = 65, SEXF = 1 (female), and
    # RACE_ASIAN_AMIND_OTH = 0 (White or Black).
    # ------------------------------------------------------------------
    lec50  <- log(3.7);   label("Tocilizumab concentration at 50% of emax (ug/mL)")             # Frey 2013 Table 1, ec50
    lemax  <- log(0.73);  label("Maximum tocilizumab effect on DAS28 production rate kin (fraction)") # Frey 2013 Table 1, emax
    lkout  <- log(0.038); label("First-order DAS28 'loss' rate kout (1/day)")                   # Frey 2013 Table 1, kout
    lhill <- log(0.64);  label("Sigmoidicity (Hill) coefficient (unitless)")                   # Frey 2013 Table 1, GAMMA
    lrbase  <- log(6.8);   label("Typical baseline DAS28 score (unitless 0-10)")                 # Frey 2013 Table 1, Baseline DAS28
    lDMARD <- log(0.30);  label("DMARD background effect, in tocilizumab concentration units (ug/mL)") # Frey 2013 Table 1, DMARD effect

    # ------------------------------------------------------------------
    # Covariate effects on PD parameters (Frey 2013 Table 2 formulas).
    # ------------------------------------------------------------------
    e_lil6_ec50    <- -4.4;   label("Power exponent of (log(IL6 * 1000)/9.9) on ec50 (unitless)")           # Frey 2013 Table 2 ec50 row
    e_sexm_emax    <-  0.11;  label("Fractional increase in emax for males (unitless)")                    # Frey 2013 Table 2 SEX row (+11% male)
    e_race_amind_oth_kout    <- -0.25;  label("Fractional change in kout for RACE_ASIAN_AMIND_OTH = 1 (unitless)")         # Frey 2013 Table 2 RACE row (-25%)
    e_blhaq_base   <-  0.043; label("Power exponent of (BLHAQ/1.6) on BASE (unitless)")                    # Frey 2013 Table 2 HAQ row
    e_lil6_base    <-  0.13;  label("Power exponent of (log(IL6 * 1000)/9.9) on BASE (unitless)")          # Frey 2013 Table 2 log-IL-6 on BASE row
    e_pain_base    <-  0.062; label("Power exponent of (PAIN/60) on BASE (unitless)")                      # Frey 2013 Table 2 PAIN row
    e_blphyvas_base <- 0.13;  label("Power exponent of (BLPHYVAS/65) on BASE (unitless)")                  # Frey 2013 Table 2 VASP row
    e_lil6_dmard   <- -6.4;   label("Power exponent of (log(IL6 * 1000)/9.9) on DMARD effect (unitless)")  # Frey 2013 Table 2 log-IL-6 on DMARD row

    # ------------------------------------------------------------------
    # IIV - PD parameters (Frey 2013 Table 1 CV%).
    # Convert each CV to NONMEM-style log-normal variance:
    #   omega^2 = log(CV^2 + 1)
    #     ec50   CV 170%  -> omega^2 = log(1.70^2 + 1) = log(3.89)   = 1.358
    #     emax   CV  11%  -> omega^2 = log(0.11^2 + 1) = log(1.0121) = 0.01203
    #     kout   CV  60%  -> omega^2 = log(0.60^2 + 1) = log(1.36)   = 0.3075
    #     BASE   CV   9.4% -> omega^2 = log(0.094^2 + 1) = log(1.00884) = 0.008805
    #     DMARD  CV 193%  -> omega^2 = log(1.93^2 + 1) = log(4.7249) = 1.553
    # Correlation ec50-emax = 0.44 (Table 1, "Correlation ec50-emax" row)
    #   cov(ec50, emax) = 0.44 * sqrt(1.358) * sqrt(0.01203)
    #                   = 0.44 * 1.166 * 0.1097 = 0.05626
    # ------------------------------------------------------------------
    etalec50 + etalemax ~ c(1.358,
                            0.05626, 0.01203)   # Frey 2013 Table 1: ec50 IIV 170% CV, emax IIV 11% CV, ec50-emax correlation 0.44
    etalkout  ~ 0.3075     # Frey 2013 Table 1: kout IIV 60% CV
    etalrbase  ~ 0.008805   # Frey 2013 Table 1: Baseline IIV 9.4% CV
    etalDMARD ~ 1.553      # Frey 2013 Table 1: DMARD effect IIV 193% CV

    # ------------------------------------------------------------------
    # Residual error -- Frey 2013 Table 1 ("Error model" section).
    # Additive residual error on the DAS28 score (which is unitless,
    # 0-10).  Standard deviation 0.68 DAS28 units.
    # ------------------------------------------------------------------
    addSd <- 0.68; label("Additive residual error on DAS28 (DAS28 units)") # Frey 2013 Table 1: additive residual 0.68
  })
  model({
    # ------------------------------------------------------------------
    # 1. Tocilizumab PK (Frey 2010 backbone, typical values; no PK IIV
    #    or PK covariates here -- this file is the exposure-response
    #    overlay).
    # ------------------------------------------------------------------
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # ------------------------------------------------------------------
    # 2. Individual DAS28 PD parameters (Frey 2013 Table 2 covariate
    #    formulas).  Floor PAIN and BLHAQ at 0.010 to keep the power
    #    form well-defined when the covariate is zero (as documented in
    #    Frey 2013 Table 2's covariate-range column).
    # ------------------------------------------------------------------
    lil6_ratio <- log(IL6 * 1000) / 9.9         # Frey 2013 Table 2: log-IL-6 power-term basis; reference 20 pg/mL gives 1.0
    haq_floored  <- max(0.010, BLHAQ) / 1.6     # Frey 2013 Table 2: HAQ floored at 0.010 to keep the power form well-defined when HAQ = 0
    pain_floored <- max(0.010, PAIN ) / 60      # Frey 2013 Table 2: PAIN floored at 0.010 to keep the power form well-defined when PAIN = 0

    ec50  <- exp(lec50  + etalec50 ) * lil6_ratio^e_lil6_ec50
    emax  <- exp(lemax  + etalemax ) * (1 + e_sexm_emax * (1 - SEXF))
    kout  <- exp(lkout  + etalkout ) * (1 + e_race_amind_oth_kout * RACE_ASIAN_AMIND_OTH)
    hill <- exp(lhill)
    rbase  <- exp(lrbase  + etalrbase ) *
             haq_floored^e_blhaq_base *
             lil6_ratio^e_lil6_base *
             pain_floored^e_pain_base *
             (BLPHYVAS / 65)^e_blphyvas_base
    DMARD <- exp(lDMARD + etalDMARD) * lil6_ratio^e_lil6_dmard
    kin   <- kout * rbase

    # ------------------------------------------------------------------
    # 3. Two-compartment tocilizumab PK with parallel linear and
    #    Michaelis-Menten clearance from the central compartment
    #    (Frey 2010 final structural model).  Concentration is in
    #    ug/mL; central is in mg (1 mg/L = 1 ug/mL).
    # ------------------------------------------------------------------
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    # ------------------------------------------------------------------
    # 4. Indirect-response DAS28 with sigmoid-emax inhibition of kin.
    #    Effective drug concentration includes the DMARD background
    #    expressed in tocilizumab concentration units.
    # ------------------------------------------------------------------
    CeffP <- Cc + DMARD
    Eff   <- emax * CeffP^hill / (ec50^hill + CeffP^hill)

    das28(0)    <- rbase
    d/dt(das28) <- kin * (1 - Eff) - kout * das28

    # ------------------------------------------------------------------
    # 5. Observation and error model.  Tocilizumab Cc is solved as a
    #    typical-value driver only (no PK residual error here -- see
    #    Frey_2010_tocilizumab.R for the PK residual-error model).
    # ------------------------------------------------------------------
    das28 ~ add(addSd)
  })
}
