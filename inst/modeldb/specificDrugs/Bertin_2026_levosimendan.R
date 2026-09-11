Bertin_2026_levosimendan <- function() {
  description <- "Joint parent-metabolite population PK model for intravenous levosimendan and its metabolites OR-1855 (inactive) and OR-1896 (active, long-lasting) in critically ill adults, neonates and infants supported by veno-arterial extracorporeal membrane oxygenation (VA-ECMO). Levosimendan disposition is two-compartment with first-order elimination; a single transit compartment carries the delayed formation of OR-1855, with the rate constant for entering the transit compartment set equal to the rate constant for leaving it. OR-1855 is either eliminated or acetylated to OR-1896, which is in turn eliminated or deacetylated back to OR-1855. Both metabolites are assumed to distribute into the levosimendan central volume (V3 = V4 = V1), which is what makes the metabolite rate constants identifiable. Allometric body weight on levosimendan clearance (exponent fixed at 0.75) and on the central volume (exponent estimated at 0.574), both referenced at 70 kg, plus body weight on the OR-1855 elimination rate constant (exponent -0.611) and a childhood (age 1 year or younger) effect that makes OR-1896 formation 3.7-fold slower in neonates and infants than in adults. The three metabolite rate constants that the 72-hour sampling window could not inform were fixed from the literature. All amounts are molar (umol) because the source analysis converted doses and concentrations to molar units so the three analytes' differing molecular weights would not distort the parent-metabolite mass transfers."
  reference <- paste(
    "Bertin S, Guidi M, Haefliger D, Thoueille P, Bardinet C, Decosterd LA,",
    "Perez MH, Giraud R, Assouline B, Schneider A, Buclin T, Livio F.",
    "Population pharmacokinetics of levosimendan and its metabolites OR-1855",
    "and OR-1896 in critically ill adults, neonates and infants on",
    "veno-arterial ECMO.",
    "Clin Pharmacokinet. 2026;65:XX. doi:10.1007/s40262-025-01591-4.",
    "Parameter values are the full-precision final estimates taken from the",
    "NONMEM control stream reproduced in Electronic Supplementary Material",
    "'Supplementary 10: NONMEM code for final model', cross-checked against",
    "Table 2 of the main paper (which rounds several of them).",
    sep = " "
  )
  vignette <- "Bertin_2026_levosimendan"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All amounts are umol: the source analysis states
  # "Levosimendan drug doses and parent-metabolite observed concentrations
  # were converted to umol/h and nmol/mL, respectively, to avoid
  # discrepancies related to the drugs' molecular weights" (Sect. 2.3.1).
  # umol/L is numerically identical to the source's nmol/mL.
  # verified = TRUE: analyte and specimen were read off the $MODEL block of
  # the source control stream (COMP = (CENTRAL), (PERIPH), (TRANSIT),
  # (METABO1), (METABO2)) together with Fig. 2's structural schematic.
  compartmentData <- list(
    central         = list(analyte = "levosimendan", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "levosimendan", units = "umol", specimen = "plasma", verified = TRUE),
    transit1        = list(analyte = "levosimendan", units = "umol", specimen = "not applicable", verified = TRUE),
    central_or1855  = list(analyte = "OR-1855", units = "umol", specimen = "plasma", verified = TRUE),
    central_or1896  = list(analyte = "OR-1896", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight is 70 kg, written explicitly into the final-model equations of Sect. 3.2 (CLi = CL * (BWi/70)^0.75, V1i = V1 * (BWi/70)^thetaBW, keM1i = keM1 * (BWi/70)^thetaBW) and into the control stream $PK block ('MWT = 70 ; MWT = MEAN BODY WEIGHT (KG)', 'RWT = WT/MWT'). It is NOT the cohort median: the cohort spans 2.7-5.8 kg (neonates/infants) and 52-125 kg (adults). Three separate exponents apply -- 0.75 fixed on CL, 0.574 estimated on V1, and -0.611 estimated on the OR-1855 elimination rate constant. Recorded as a baseline weight; the source collected weight once, on the first day of sampling (Sect. 2.1). Source column WT.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model ONLY through the binary childhood indicator that the control stream derives from it: 'IF(AGE.LE.1) Q1=1' / 'IF(AGE.GT.1) Q1=0', i.e. age 1 year or younger marks the neonate/infant group whose OR-1896 formation rate constant carries the e_child_kmet_or1896 effect. The main paper describes this covariate only as 'childhood' or 'neonates/infants vs adults' and never states the cut-off; the 1-year threshold is recoverable solely from the control stream. Age is otherwise not a covariate in the final model -- a linear age effect on V2 and on ktransit was significant in the univariate screen (ESM Supplementary 2, Step 1) but was not retained. The study's own groups are consistent with the cut-off: adults were 18-75 years and the neonates/infants 13-164 days (Table 1). Source column AGE."
    )
  )

  # Screened in the ESM covariate analysis but NOT retained in the final
  # model, so they are documented rather than declared as live covariates.
  # ESM Supplementary 2 Step 1 found each of these statistically significant
  # in the univariate screen; most were dropped for imprecision (RSE > 50%),
  # and CRRT was dropped in Step 2 despite a large dOFV because its effect
  # "could not be estimated precisely enough because of the small number of
  # patients involved" (Sect. 4).
  covariatesDataExcluded <- list(
    CRRT = list(
      description = "Continuous renal replacement therapy (continuous veno-venous haemodialysis)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on keM1 (dOFV -14.71, p < 0.001) and keM2 (dOFV -4.24, p < 0.05) in ESM Supplementary 2 Step 1, and retained into Step 2 where it gave dOFV -14.047 on keM1 but an effect estimate of 7.42 whose RSE could not be calculated. Not in the final model. Four of 21 patients were on CRRT and showed low metabolite concentrations; the paper warns clinically that CRRT may reduce OR-1855 and OR-1896 exposure but declines to quantify it (Sect. 4)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a power function on keM1 (dOFV -9.40) and keM2 (dOFV -10.28) in ESM Supplementary 2 Step 1; both effects had RSE > 50% and neither was retained. Body weight carried the size effect in the final model."
    ),
    SEXF = list(
      description = "Female sex",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as a categorical effect on keM2 (dOFV -9.31, p < 0.01) in ESM Supplementary 2 Step 1 with RSE 91%; not retained. The paper lists sex among the covariates whose high RSEs 'precluded any definitive conclusions' (Sect. 5)."
    ),
    ECMO_FLOW = list(
      description = "Veno-arterial ECMO circuit blood flow rate",
      units       = "L/min",
      type        = "continuous",
      notes       = "Screened as a linear effect on keM2 (dOFV -9.96, p < 0.01, RSE 60%) in ESM Supplementary 2 Step 1; not retained. Median 2.7 L/min (1.1-4.6) in adults and 0.5 L/min (0.2-1.1) in neonates/infants (Table 1). Named in Sect. 5 as one of the covariates the cohort size could not resolve."
    ),
    N_COMED = list(
      description = "Number of concomitant medications",
      units       = "(count)",
      type        = "continuous",
      notes       = "Screened as a linear effect on kM2 (dOFV -19.37, p < 0.0001, RSE 98%) and on kM2-M1 (dOFV -4.31, p < 0.05, RSE 63%) in ESM Supplementary 2 Step 1; neither was retained, both being far too imprecise. This was the single largest univariate dOFV in the whole screen, which the paper does not comment on."
    ),
    GFR = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min/1.73m2",
      type        = "continuous",
      notes       = "Tested (Sect. 2.3.2 lists it among the covariates screened, restricted to patients not on CRRT) but not significant and absent from the ESM Supplementary 2 Step 1 table. Sect. 4 explains why: renal impairment was only moderate in this cohort (adult median 78 mL/min/1.73m2), 'preventing GFR to be retained as a significant covariate on metabolite elimination', even though severe renal failure is known to raise OR-1855 and OR-1896 exposure."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested (Sect. 2.3.2) but not significant and absent from the ESM Supplementary 2 Step 1 table. Most patients were hypoalbuminaemic (median 28 g/L in adults, 30 g/L in neonates/infants; Table 1), which the paper notes would raise the free fraction of this highly protein-bound drug without changing total-concentration PK (Sect. 4)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested (Sect. 2.3.2) but not significant and absent from the ESM Supplementary 2 Step 1 table. Median 17 umol/L (8-116) in adults and 13 umol/L (6-38) in neonates/infants (Table 1)."
    ),
    ABX = list(
      description = "Concomitant antibiotic administration",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Tested (Sect. 2.3.2) but not significant: 'All the other covariates, including the use of antibiotics, were not associated with levosimendan PK' (ESM Supplementary 2). Notable because OR-1855 is formed by gut microbiota, so the paper hypothesises in Sect. 4 that antibiotic-driven dysbiosis may depress OR-1855 synthesis even though this cohort (67% of adults, 100% of neonates/infants on antibiotics) could not demonstrate it."
    ),
    T_ECMO = list(
      description = "Time since ECMO initiation",
      units       = "h",
      type        = "continuous",
      notes       = "Tested (Sect. 2.3.2) but not significant and absent from the ESM Supplementary 2 Step 1 table. Relevant to the saturable-sequestration hypothesis the paper raises in Sect. 4."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21,
    n_studies      = 1,
    n_samples      = 155,
    age_range      = "adults 18-75 years; neonates/infants 13-164 days",
    age_median     = "adults 62 years; neonates/infants 24 days",
    weight_range   = "adults 52-125 kg; neonates/infants 2.7-5.8 kg",
    weight_median  = "adults 79 kg; neonates/infants 3.4 kg",
    sex_female_pct = 33.3,
    disease_state  = "Critically ill intensive-care patients supported by veno-arterial extracorporeal membrane oxygenation (VA-ECMO), all 21 of them, for cardiac arrest (10), cardiogenic shock (7) or failure to wean from cardiopulmonary bypass after cardiac surgery (4). Four patients were on continuous veno-venous haemodialysis. Renal function was only moderately impaired (adult median estimated GFR 78 mL/min/1.73m2 by CKD-EPI; neonate/infant median 24 mL/min/1.73m2 by the Schwartz formula or 24-hour urine collection). Most were hypoalbuminaemic (median 28 and 30 g/L). ICU mortality was 33% in adults and 83% in neonates/infants.",
    dose_range     = "Adults: levosimendan started at 0.05 ug/kg/min for 1-4 h then increased to a maintenance rate of 0.1 (n = 9), 0.15 (n = 3) or 0.2 (n = 3) ug/kg/min, infused for about 24 h (23.5-29 h) in 14 of 15 and 6 h in one. Neonates/infants: a continuous 0.1 ug/kg/min infusion for 48 h. Two paediatric patients received two separate infusions 6 and 7 days apart.",
    regions        = "Switzerland -- Lausanne University Hospital (CHUV) adult and paediatric intensive care units and Geneva University Hospitals (HUG) intensive care, bicentric prospective observational study approved December 2022 (project-ID 2022-01262).",
    notes          = "Baseline demographics from Table 1. Sampling: adults at 1, 2, 4, 24, 25, 26, 28 and 48 h after the start of the infusion; paediatric patients at 1, 2, 4, 24, 48, 49, 52 and 72 h, a maximum of eight samples per patient and occasion (median 8, range 4-16). Assay LLOQ was 0.1 ng/mL for all three analytes by UHPLC-MS/MS. Levosimendan, OR-1855 and OR-1896 were below that limit in 5 (3%), 54 (35%) and 68 (44%) of samples respectively, handled during estimation by the M3 likelihood method. OR-1855 concentrations from one neonate were excluded from the analysis because residual drug from an earlier infusion made them decline throughout the sampling period. Because sampling stopped 24 h after the end of the infusion, metabolite elimination was never observed, which is why keM1, keM2 and kM2-M1 are fixed from the literature rather than estimated."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Values are the FINAL estimates read from the
    # $THETA block of the control stream in ESM Supplementary 10, which
    # carries more significant figures than Table 2 of the main paper.
    # Where the two disagree it is only by rounding, EXCEPT for ktransit:
    # Table 2 prints '0.01' where the control stream has 0.0127, a 27%
    # difference that matters because ktransit sets the whole timescale of
    # metabolite appearance. The control stream value is used.
    #
    # Clearances and volumes refer to a 70-kg patient.
    # ------------------------------------------------------------------
    lcl <- log(13.9);   label("Levosimendan elimination clearance CL at WT = 70 kg (L/h)")                    # $THETA 1 '13.9 ; CL'; Table 2 'CL (L/h) 14 (21%)'
    lvc <- log(15.9);   label("Levosimendan central volume of distribution V1 at WT = 70 kg (L)")             # $THETA 3 '15.9 ; V1'; Table 2 'V 1 (L) 16 (26%)'
    lq  <- log(0.501);  label("Levosimendan intercompartmental clearance Q (L/h)")                            # $THETA 5 '0.501 ;Q'; Table 2 'Q (L/h) 0.50 (36%)'
    lvp <- log(5.75);   label("Levosimendan peripheral volume of distribution V2 (L)")                        # $THETA 6 '5.75 ; V2'; Table 2 'V 2 (L) 5.8 (44%)'

    # The transit compartment that produces the delayed appearance of
    # OR-1855. One transit compartment beat direct conversion by
    # dAIC = +96, two transit compartments by +157, and a Savic / Stirling
    # approximation by +103 (Sect. 3.2). The SAME rate constant governs
    # entry into and exit from the compartment -- the control stream sets
    # 'K34 = K13' -- so this is one parameter used twice, not two.
    lktr <- log(0.0127); label("Rate constant for levosimendan entering and leaving the metabolite-formation transit compartment, ktransit (1/h)")  # $THETA 7 '0.0127 ; K13' with $PK 'K34 = K13'; Table 2 rounds to 'k transit (h -1 ) 0.01 (21%)'

    # ------------------------------------------------------------------
    # Metabolite rate constants. The study sampled only to 24 h after the
    # end of the infusion, so it never observed metabolite elimination;
    # Sect. 2.3.1 fixes keM1 and keM2 at 0.01 1/h, "corresponding to their
    # half-life of about 70 h", and fixes the OR-1896 to OR-1855 back
    # conversion at 0.012 1/h. Wrapped in fixed() because the control
    # stream marks all three 'FIX'.
    #
    # ESM Supplementary 1 derives the back-conversion constant from
    # Puttonen et al. 2007, in which 45% of an administered OR-1896 dose
    # was recovered unchanged in urine, so the fraction back-converted was
    # taken as FM2-M1 = 0.55 and kM2-M1 = keM2 * FM2-M1 / (1 - keM2) =
    # 0.01 * 0.55 / (1 - 0.55) = 0.012 1/h. Reproduced here exactly as the
    # ESM prints it, including its own inconsistency: the denominator is
    # written symbolically as (1 - keM2) but evaluated as (1 - 0.55), i.e.
    # (1 - FM2-M1). The arithmetic the ESM actually performs, 0.01 * 0.55 /
    # 0.45 = 0.0122, is what rounds to the 0.012 that the $THETA block
    # fixes, so the fixed value is unambiguous regardless of which
    # denominator was intended. Note that Table 2's "Final
    # model estimate" column prints this as '0.01 FIX' while its own
    # bootstrap column, the main text and the control stream all say
    # 0.012; 0.012 is used.
    # ------------------------------------------------------------------
    lkel_or1855  <- fixed(log(0.01));  label("OR-1855 elimination rate constant keM1 at WT = 70 kg (1/h)")                             # $THETA 8 '0.01 FIX ; K40'; Table 2 'k eM1 (h -1 ) 0.01 FIX'
    lkmet_or1896 <- log(0.0722);       label("Rate constant for acetylation of OR-1855 to OR-1896 in adults, kM2 (1/h)")               # $THETA 9 '0.0722 ; K45 adult'; Table 2 'k M2 (h -1 ) 0.07 (44%)'
    lkel_or1896  <- fixed(log(0.01));  label("OR-1896 elimination rate constant keM2 (1/h)")                                           # $THETA 10 '0.01 FIX ; K50'; Table 2 'k eM2 (h -1 ) 0.01 FIX'
    lkicv_or1855 <- fixed(log(0.012)); label("Rate constant for back conversion of OR-1896 to OR-1855, kM2-M1 (1/h)")                  # $THETA 11 '0.012 FIX ; K54'; ESM Supplementary 1 derivation; Sect. 2.3.1

    # ------------------------------------------------------------------
    # Covariate effects. All three body-weight terms are power functions
    # on WT/70; the childhood term is the categorical form
    # Param * (1 + COV * theta) given in Sect. 2.3.1.
    # ------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on levosimendan clearance (unitless)")                # $THETA 2 '0.75 FIX ; WT on CL'; Table 2 'theta BW 0.75 FIX'. Sect. 3.2: fixed at the allometric-theory value because the estimate's CI contained it
    e_wt_vc <- 0.574;       label("Exponent of body weight on the levosimendan central volume (unitless)")                  # $THETA 4 '0.574 ; WT on V1'; Table 2 'theta BW 0.57 (32%)'
    e_wt_kel_or1855 <- -0.611; label("Exponent of body weight on the OR-1855 elimination rate constant (unitless)")         # $THETA 13 '-0.611 ; WT on K40'; Table 2 'theta BW - 0.61 (22%)'
    e_child_kmet_or1896 <- -0.732; label("Fractional change in the OR-1855 to OR-1896 acetylation rate constant for patients aged 1 year or younger (unitless)")  # $THETA 12 '-0.732 ; K45 children'; Table 2 'theta neonates/infants - 0.73 (30%)'

    # ------------------------------------------------------------------
    # Between-subject variability. The $OMEGA block holds VARIANCES
    # (NONMEM definition), so they are used here as written. Cross-check
    # against Table 2, whose footnote b gives the conversion
    # CV% = sqrt(exp(omega^2) - 1); every one of the five reproduces the
    # published CV to two significant figures.
    #
    # The control stream carries nine etas, four of which are '0 FIX' --
    # on Q, keM1, keM2 and kM2-M1. Sect. 3.2 explains that adding BSV to
    # any of those four gave runs that did not converge. They are omitted
    # here rather than written as ~ fixed(0), because a zero on the
    # variance diagonal makes OMEGA singular and breaks the Cholesky
    # sampler that rxSolve uses to draw the cohort.
    # ------------------------------------------------------------------
    etalcl ~ 0.0998            # $OMEGA 1 'ETA CL';  sqrt(exp(0.0998)-1) = 32.4% CV = Table 2 'omega CL 32 (31%)'
    etalvc ~ 0.235             # $OMEGA 2 'ETA V1';  sqrt(exp(0.235)-1)  = 51.5% CV = Table 2 'omega V1 52 (48%)'
    etalvp ~ 0.684             # $OMEGA 4 'ETA V2';  sqrt(exp(0.684)-1)  = 99.1% CV = Table 2 'omega V2 99 (46%)'
    etalktr ~ 0.128            # $OMEGA 5 'ETA K13'; sqrt(exp(0.128)-1)  = 37.0% CV = Table 2 'omega ktransit 37 (71%)'
    etalkmet_or1896 ~ 0.599    # $OMEGA 7 'ETA K45'; sqrt(exp(0.599)-1)  = 90.6% CV = Table 2 'omega kM2 91 (88%)'

    # ------------------------------------------------------------------
    # Residual error. The source $ERROR block forms IPRED on the natural
    # log scale and adds the error there -- 'IPRED1 = LOG(A(1)/S1)'
    # followed by 'Y1 = IPRED1+ERR(1)' -- which is exactly nlmixr2's
    # ~ lnorm(expSd), the same translation used by the existing
    # Wattanakul_2024_primaquine.R extraction of an identically shaped
    # control stream. Sect. 2.3.1 describes it as "an additive error model
    # in the log-transformed domain for each compound, corresponding to a
    # proportional error in the raw concentration scale"; lnorm() is the
    # exact form and prop() the small-sigma approximation to it.
    #
    # $SIGMA holds VARIANCES, so each SD is written as sqrt() of the
    # published variance to keep the source number visible. Table 2's
    # "CV%" row for each analyte is sqrt(sigma^2) to two figures, which
    # confirms the variance reading.
    # ------------------------------------------------------------------
    expSd        <- sqrt(0.0952); label("Residual SD on the natural-log scale for levosimendan (log units)")  # $SIGMA 1 '0.0952 ; Prop Err Levosimendan'; sqrt = 0.309 = Table 2 'sigma levosimendan 31 (16%)'
    expSd_or1855 <- sqrt(0.135);  label("Residual SD on the natural-log scale for OR-1855 (log units)")       # $SIGMA 2 '0.135 ; Prop Err OR-1855'; sqrt = 0.367 = Table 2 'sigma M1 37 (32%)'
    expSd_or1896 <- sqrt(0.0916); label("Residual SD on the natural-log scale for OR-1896 (log units)")       # $SIGMA 3 '0.0916 ; Prop Err OR-1896'; sqrt = 0.303 = Table 2 'sigma M2 30 (42%)'
  })

  model({
    # ---- 1. Derived covariate terms -----------------------------------
    # The control stream derives the childhood indicator from age:
    #   IF(AGE.LE.1) Q1=1
    #   IF(AGE.GT.1) Q1=0
    # so the group boundary is one year, which only the control stream
    # states. Every neonate/infant in the study was 13-164 days old and
    # every adult 18 years or older, so no subject sat near the cut-off.
    child <- (AGE <= 1)

    # ---- 2. Individual parameters -------------------------------------
    # $PK: TVCL = THETA(1)*RWT**THETA(2), CL = TVCL*EXP(ETA(1)), with
    # RWT = WT/70. Sect. 3.2 equations CLi, V1i, keM1i, kM2i.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)                                   # ETA(3) on Q is '0 FIX'
    vp <- exp(lvp + etalvp)
    ktr <- exp(lktr + etalktr)

    kel_or1855  <- exp(lkel_or1855) * (WT / 70)^e_wt_kel_or1855        # ETA(6) on K40 is '0 FIX'
    kmet_or1896 <- exp(lkmet_or1896 + etalkmet_or1896) * (1 + child * e_child_kmet_or1896)
    kel_or1896  <- exp(lkel_or1896)                                    # ETA(8) on K50 is '0 FIX'
    kicv_or1855 <- exp(lkicv_or1855)                                   # ETA(9) on K54 is '0 FIX'

    # ---- 3. Micro-constants -------------------------------------------
    # $PK: K10 = CL/V1-K13, K12 = Q/V1, K21 = Q/V2.
    #
    # The subtraction is load-bearing and easy to get backwards. The flux
    # into the transit compartment is carved OUT of the total clearance
    # rather than added on top of it, so levosimendan still leaves the
    # central compartment at exactly CL/V1 in total, of which the ktr
    # share goes on to become OR-1855 and the remainder is the 95% of the
    # dose that Sect. 1 says is lost to the glutathione pathway. Encoding
    # it additively instead would inflate total elimination by ktr and
    # shift the published distribution half-life from 0.76 h to 0.75 h;
    # the vignette checks both half-lives against Sect. 3.2.
    kel <- cl / vc - ktr
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ------------------------------------------------
    # Transcribed from $DES. Metabolite concentrations are formed on the
    # levosimendan central volume because the control stream sets
    # V4 = V1 and V5 = V1 (S4 = S1, S5 = S1): Sect. 2.3.1 explains that
    # the metabolites' own volumes are not identifiable from
    # parent-metabolite data, so they were assumed equal to V1.
    d/dt(central)        <- -k12 * central + k21 * peripheral1 - ktr * central - kel * central
    d/dt(peripheral1)    <-  k12 * central - k21 * peripheral1
    d/dt(transit1)       <-  ktr * central - ktr * transit1
    d/dt(central_or1855) <-  ktr * transit1 - kmet_or1896 * central_or1855 -
                             kel_or1855 * central_or1855 + kicv_or1855 * central_or1896
    d/dt(central_or1896) <-  kmet_or1896 * central_or1855 - kel_or1896 * central_or1896 -
                             kicv_or1855 * central_or1896

    # ---- 5. Observations and error ------------------------------------
    Cc        <- central / vc
    Cc_or1855 <- central_or1855 / vc
    Cc_or1896 <- central_or1896 / vc

    Cc        ~ lnorm(expSd)
    Cc_or1855 ~ lnorm(expSd_or1855)
    Cc_or1896 ~ lnorm(expSd_or1896)
  })
}
