Storset_2014_tacrolimus <- function() {
  description <- "Theory-based two-compartment population pharmacokinetic model for oral tacrolimus in adult kidney-transplant recipients (Storset 2014): plasma-based disposition with first-order absorption and a lag time, allometric scaling on fat-free mass, CYP3A5-expresser effects on plasma clearance and oral bioavailability, a sigmoid-Emax prednisolone-driven reduction in bioavailability, a first-day-post-transplant bioavailability spike with subject-level random effect, and a saturable haematocrit-dependent red-blood-cell-binding equation that maps plasma concentration to whole-blood concentration."
  reference <- "Storset E, Holford N, Hennig S, Bergmann TK, Bergan S, Bremer S, Asberg A, Midtvedt K, Staatz CE. Improved prediction of tacrolimus concentrations early after kidney transplantation using theory-based pharmacokinetic modelling. Br J Clin Pharmacol. 2014;78(3):509-523. doi:10.1111/bcp.12361"
  vignette <- "Storset_2014_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    FFM = list(
      description        = "Predicted fat-free mass (Janmahasatian 2005 formula from total body weight, height, and sex)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on plasma CL/F, V1/F, Q/F, V2/F with reference 60 kg (Storset 2014 Methods Equation 2; Table 2 footnote: parameters reported at FFM 60 kg). Theory-based exponents fixed at 3/4 on clearances (CL/F, Q/F) and 1 on volumes (V1/F, V2/F). Cohort median 59 kg, range 35-80 kg (Storset 2014 Table 1 model-development column).",
      source_name        = "FFM"
    ),
    HCT = list(
      description        = "Haematocrit -- packed red-blood-cell volume fraction",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (changes substantially in the first weeks post-transplant). Used inside the saturable RBC-binding equation Cwb = Cp * (1 + Bmax * (HCT/100) / (Cp + KD)) (Storset 2014 Methods Equation 3) to map plasma to whole-blood concentration. Bmax = 418 ug/L erythrocytes and KD = 3.8 ug/L plasma are fixed at literature values (Storset 2014 Methods, citing reference 35; Jusko 1995 Clin Pharmacol Ther). Cohort medians 33% Brisbane / 36% Oslo, range 24-45% (Storset 2014 Table 1).",
      source_name        = "Haematocrit"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 functional-expression indicator (rs776746 / CYP3A5*3 polymorphism)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 non-expresser)",
      notes              = "1 = expresser (homozygous *1/*1 wild-type or heterozygous *1/*3 -- pooled because Storset 2014 had only n = 3 *1/*1 subjects); 0 = non-expresser (homozygous *3/*3). Time-fixed per subject. Effects: multiplicative on plasma CL/F (factor 1.30 for expressers; Storset 2014 Table 2) and on oral bioavailability F (factor 0.82 for expressers). Cohort distribution 3 / 53 / 205 (*1/*1, *1/*3, *3/*3 = 22.6% expressers; Storset 2014 Table 1 Hardy-Weinberg-equilibrium).",
      source_name        = "CYP3A5 expresser"
    ),
    PRED_DOSE = list(
      description        = "Concomitant oral prednisolone daily dose",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying with the post-transplant conmed_steroid taper (Storset 2014 Brisbane: 0.3 mg/kg ideal body weight per day initial, max 30 mg/day; Oslo: 20 mg/day initial, up to 80 mg/day in selected patients). Effect on F is a sigmoid-Emax fractional reduction (1 - Pred_max * PRED_DOSE / (Pred_50 + PRED_DOSE)) with Pred_max = 0.67 (95% CI 41%, 89%) and Pred_50 = 35 mg/day (95% CI 7, 50) -- Storset 2014 Table 2 and Methods Equations 4 + 6 with Hill = 1. Cohort medians 6.0 mg/day Brisbane / 0 mg/day Oslo (Storset 2014 Table 1).",
      source_name        = "Prednisolone dose"
    ),
    POSTTX_DAY1 = list(
      description        = "First-24-hours-post-transplant indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any observation outside the first 24 hours post-transplant)",
      notes              = "Time-varying: 1 only on the first day post-transplant (typically the day of the transplant surgery / first oral tacrolimus dose), 0 thereafter. Multiplicative ~2.68-fold increase in oral bioavailability on day 1 (Storset 2014 Table 2; 95% CI 2.28, 3.09; OFV decrease 209), with subject-level IIV of 57% CV on the day-1 multiplier. Storset 2014 Discussion attributes the effect to candidate mechanisms (methylprednisolone-bolus inhibition of intestinal CYP3A / P-glycoprotein, surgery-related inflammation, anaesthesia / opioid effects on gut motility, reduced food intake) without a single mechanistic explanation.",
      source_name        = "first day post-transplant"
    )
  )

  population <- list(
    n_subjects     = 242L,
    n_studies      = 2L,
    age_range      = "23-71 years",
    age_median     = "48 years",
    weight_range   = "51-121 kg total body weight; 35-80 kg fat-free mass",
    weight_median  = "80 kg total body weight; 59 kg fat-free mass",
    sex_female_pct = 31.8,
    race_ethnicity = "Not reported (Brisbane and Oslo cohorts; predominantly Caucasian inferred from CYP3A5 *3/*3 frequency 84.7%)",
    disease_state  = "Adult kidney-transplant recipients receiving oral tacrolimus, mycophenolate mofetil, and tapering prednisolone immunosuppression after high-dose IV methylprednisolone + basiliximab induction",
    dose_range     = "Tacrolimus oral starting dose 0.075 mg/kg twice daily (Brisbane) or 0.04 mg/kg twice daily (Oslo); subsequently adjusted to whole-blood trough targets 7-8 ug/L (Brisbane) or 3-7 ug/L (Oslo) over the first 3 months post-transplant; reported total daily doses median 8.4 mg (range 4.2-35.9 mg/day across both centres, Storset 2014 Table 1)",
    regions        = "Australia (Princess Alexandra Hospital, Brisbane: 173 subjects, 1546 concentrations) and Norway (Oslo University Hospital Rikshospitalet: 69 subjects, 1554 concentrations); all model-development data",
    cyp3a5_distribution = "*1/*1 = 1.2% (3/241), *1/*3 = 22.0% (53/241), *3/*3 = 84.7% (205/241); 22.6% expressers overall (Storset 2014 Table 1; one *2 subject excluded from genotype tally)",
    haematocrit_distribution = "Brisbane median 33% (range 25-43%); Oslo median 36% (range 24-45%) -- time-varying within subject across the first weeks post-transplant",
    sampling_design = "3100 whole-blood tacrolimus concentrations (median 11 per subject, range 4-52); 791 (26%) full PK profiles (>8 samples per occasion), 1277 (41%) limited profiles (4-5 samples per occasion), 1032 (33%) routine trough concentrations; sampling time median 20 days post-transplant, range 5-2591 days",
    notes          = "Pooled analysis of two cohorts previously analysed independently (Bergmann 2014 Brisbane / Asberg-Storset 2014 Oslo). External evaluation cohort: 72 additional Oslo subjects with 837 trough samples in the first 3 weeks post-transplant. Brisbane assays: LC-MS/MS (LLOQ 0.5 ug/L). Oslo assays: chemiluminescent microparticle immunoassay (80%) / LC-MS/MS (11%) / microparticle enzyme immunoassay (9%); immunoassay values converted to LC-MS/MS equivalents via C_LC-MS/MS = 0.80 * C_immunoassay + 0.19 (Storset 2014 Equation 1). Baseline demographics from Storset 2014 Table 1."
  )

  ini({
    # ----- Structural plasma-based PK parameters -----
    # Disposition is parameterized in plasma terms (Storset 2014 Methods: 'whole
    # blood concentrations were measured the model was parameterized in terms
    # of plasma concentration based disposition parameters'). Reference values
    # correspond to the typical subject FFM = 60 kg, CYP3A5 non-expresser.
    # The 'original model was as follows' block in the Storset 2014 Table 2
    # footnote gives the plasma-based equations:
    #   CLp/F = 811 * (FFM/60)^(3/4) * 1.30^CYP3A5_expresser   L/h
    #   V1p/F = 6290 * (FFM/60)                                L
    #   Qp/F  = 1200 * (FFM/60)^(3/4)                          L/h
    #   V2p/F = 32100 * (FFM/60)                               L
    lka  <- log(1.01);   label("Absorption rate constant (Ka, 1/h)")                               # Storset 2014 Table 2 final theory-based model (RSE 9%)
    ltlag <- log(0.41);  label("Absorption lag time (Tlag, h)")                                    # Storset 2014 Table 2 final theory-based model (RSE 8%)
    lcl  <- log(811);    label("Apparent plasma clearance at FFM 60 kg, CYP3A5 non-expresser (CLp/F, L/h)")    # Storset 2014 Table 2 footnote (original model equation); RSE 10% on the haematocrit-standardized estimate
    lvc  <- log(6290);   label("Apparent plasma central volume at FFM 60 kg (V1p/F, L)")           # Storset 2014 Table 2 footnote (original model equation); RSE 11%
    lq   <- log(1200);   label("Apparent plasma intercompartmental clearance at FFM 60 kg (Qp/F, L/h)") # Storset 2014 Table 2 footnote (original model equation); RSE 13%
    lvp  <- log(32100);  label("Apparent plasma peripheral volume at FFM 60 kg (V2p/F, L)")        # Storset 2014 Table 2 footnote (original model equation); RSE 16%

    # ----- Allometric exponents on FFM (theory-based, fixed) -----
    # Storset 2014 Methods Equation 2: 'P = Pstd * (SIZE/SIZE_std)^b' with b
    # 'fixed to theory-based values of 3/4 for clearances and 1 for volumes
    # of distribution'.
    e_ffm_cl <- fixed(0.75); label("Allometric exponent on plasma CL/F and Q/F (unitless, theory-based)") # Storset 2014 Methods Equation 2; theory-based 3/4 fixed
    e_ffm_vc <- fixed(1.00); label("Allometric exponent on plasma V1/F and V2/F (unitless, theory-based)") # Storset 2014 Methods Equation 2; theory-based 1 fixed

    # ----- Covariate effects on plasma CL -----
    # Storset 2014 Table 2 footnote: 'CLp/F = 811 * (FFM/60)^(3/4) * 1.30 (if
    # CYP3A5 expresser)'. Implemented as a multiplicative power-form factor
    # cl *= 1.30^CYP3A5_EXPR, so a non-expresser (CYP3A5_EXPR = 0) gets factor 1
    # and an expresser (CYP3A5_EXPR = 1) gets factor 1.30.
    e_cyp3a5_exp_cl <- log(1.30); label("Log-multiplicative effect of CYP3A5 expresser on plasma CL/F (unitless)")    # Storset 2014 Table 2 (CYP3A5 expresser CL factor 1.30; 95% CI 1.13, 1.46)

    # ----- Covariate effects on bioavailability (F) -----
    # Storset 2014 Table 2 footnote: 'F = 1 * [1 - (0.67 * Predisone)/(35 +
    # Predisone)] * 0.82 (if CYP3A5 expresser) * 2.68 (if first day
    # post-transplant)'.
    pred_max <- 0.67; label("Maximum fractional reduction in F due to prednisolone induction (unitless)")      # Storset 2014 Table 2 Predmax = -67%; 95% CI -41%, -89% (RSE 19%)
    pred_50  <- 35;   label("Prednisolone daily dose causing half-maximum reduction in F (Pred50, mg/day)")    # Storset 2014 Table 2 Pred50 = 35 mg; 95% CI 7, 50 (RSE 40%)
    e_cyp3a5_exp_fdepot <- log(0.82); label("Log-multiplicative effect of CYP3A5 expresser on F (unitless)")    # Storset 2014 Table 2 (CYP3A5 expresser F factor 0.82; 95% CI 0.71, 0.98)
    lfday1 <- log(2.68); label("Log-multiplicative effect of first day post-transplant on F (unitless)")        # Storset 2014 Table 2 (Fday1 factor 2.68; 95% CI 2.28, 3.09; RSE 8%)

    # ----- IIV (Storset 2014 Table 2 final theory-based model) -----
    # Log-normal random effects with omega^2 = log(1 + CV^2). Reported BSV CVs:
    # CLp/F = 40%, V1p/F = 54%, Qp/F = 63%, Fday1 (the day-1 multiplier) = 57%.
    # Correlations from Table 2: corr(CL, V1) = 0.43, corr(CL, Q) = 0.62. The
    # V1-Q correlation is not reported and is treated as 0 (partial omega
    # block: V1-Q covariance fixed at 0).
    #   omega^2: CL = log(1+0.40^2) = 0.14842; V1 = log(1+0.54^2) = 0.25612;
    #             Q = log(1+0.63^2) = 0.33429; Fday1 = log(1+0.57^2) = 0.28150.
    #   covariances:
    #     cov(CL, V1) = 0.43 * sqrt(0.14842 * 0.25612) = 0.08385
    #     cov(CL, Q)  = 0.62 * sqrt(0.14842 * 0.33429) = 0.13809
    #     cov(V1, Q)  = 0 (not reported; treated as 0)
    etalcl + etalvc + etalq ~ c(0.14842,
                                0.08385, 0.25612,
                                0.13809, 0,       0.33429)
    # IIV on the day-1 bioavailability multiplier; subject's eta only
    # contributes when POSTTX_DAY1 = 1 (model() block applies it as
    # exp((lfday1 + etalfday1) * POSTTX_DAY1)).
    etalfday1 ~ 0.28150  # 57% CV (Storset 2014 Table 2 BSV on the day-1 F multiplier)

    # ----- Residual error -----
    # Proportional on whole-blood Cc (Storset 2014 Table 2 final-model
    # residual error 14.9%, RSE 4%). The additive component was approached
    # zero in the full model and was removed.
    propSd <- 0.149; label("Proportional residual error on whole-blood Cc (fraction)")    # Storset 2014 Table 2 final theory-based model (14.9%; RSE 4%)
  })

  model({
    # ----- Individual PK parameters -----
    # Plasma-based clearance and volumes scaled by FFM (theory-based
    # allometry, fixed exponents) and CYP3A5 expresser status. Reference
    # values: typical subject FFM = 60 kg, CYP3A5 non-expresser
    # (CYP3A5_EXPR = 0).
    ka   <- exp(lka)
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl) * (FFM / 60)^e_ffm_cl * exp(e_cyp3a5_exp_cl * CYP3A5_EXPR)
    vc   <- exp(lvc + etalvc) * (FFM / 60)^e_ffm_vc
    q    <- exp(lq  + etalq)  * (FFM / 60)^e_ffm_cl
    vp   <- exp(lvp)          * (FFM / 60)^e_ffm_vc

    # Oral bioavailability F = baseline * F_pred * F_cyp * F_day1.
    # F_pred:  sigmoid-Emax fractional reduction by prednisolone daily dose.
    # F_cyp:   CYP3A5-expresser factor (0.82 for expressers, 1 otherwise).
    # F_day1:  ~2.68-fold spike on day 1 post-transplant; subject-level eta
    #          applied multiplicatively on the log scale (only kicks in when
    #          POSTTX_DAY1 = 1; on later days the term is 1 and the eta has
    #          no effect on F).
    fpred  <- 1 - pred_max * PRED_DOSE / (pred_50 + PRED_DOSE)
    fcyp   <- exp(e_cyp3a5_exp_fdepot * CYP3A5_EXPR)
    fday1  <- exp((lfday1 + etalfday1) * POSTTX_DAY1)
    fdepot <- fpred * fcyp * fday1

    # ----- Micro-rate constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (two-compartment with first-order absorption + lag) -----
    # central holds the apparent plasma amount (in plasma-equivalent units);
    # Cp = central / vc is the plasma concentration that drives the saturable
    # RBC-binding equation below.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Bioavailability and lag time on the depot compartment.
    f(depot)   <- fdepot
    lag(depot) <- tlag

    # ----- Plasma concentration -----
    # The ODE state `central` carries the apparent plasma amount in dose
    # units (mg). Plasma concentration in mg/L = central / vc; multiply by
    # 1000 to convert to ug/L (= ng/mL), the unit Storset 2014 reports
    # tacrolimus concentrations in throughout (Table 1, Figures 1-4) and
    # the unit in which the saturable RBC-binding constants Bmax and KD
    # are given.
    Cp <- 1000 * central / vc

    # ----- Whole-blood concentration via saturable RBC binding -----
    # Storset 2014 Methods Equation 3: Cwb = Cp * (1 + Bmax * fHCT / (Cp + KD))
    # with Bmax = 418 ug/L erythrocytes and KD = 3.8 ug/L plasma fixed at
    # literature values (Storset 2014 Methods, citing reference 35; Jusko
    # 1995). HCT is supplied as a percentage and converted to fraction here.
    # Cc is the observation (whole-blood tacrolimus concentration, ug/L).
    Bmax_rbc <- 418    # ug/L erythrocytes (Storset 2014 ref [35])
    KD_rbc   <- 3.8    # ug/L plasma (Storset 2014 ref [35])
    Cc <- Cp * (1 + Bmax_rbc * (HCT / 100) / (Cp + KD_rbc))
    Cc ~ prop(propSd)
  })
}
