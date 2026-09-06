Pan_2026_adalimumab <- function() {
  description <- paste(
    "Population PK-PD model for subcutaneous adalimumab in adults with",
    "moderate-to-severe plaque psoriasis, developed on a UK real-world",
    "cohort (BSTOP / PSORT-D) to evaluate a proactive therapeutic drug",
    "monitoring strategy.",
    "PK: one-compartment disposition with first-order absorption and",
    "elimination, in APPARENT terms (F was fixed to 1 in the control",
    "stream, so CL and V are CL/F and V/F); ka 0.268 /day, CL/F 0.386",
    "L/day at the reference covariate vector, V/F fixed at 10.8 L from",
    "Ternant 2015. Clearance carries five covariates: allometric weight",
    "(exponent 0.75 fixed), a power effect of anti-drug-antibody level",
    "((ADA/76.03)^0.368), a power effect of waist circumference",
    "((WAIST/101)^0.888), and proportional increases for female sex",
    "(+21.6%) and hypertension (+17.7%); volume carries allometric weight",
    "with the exponent fixed at 1.",
    "PD: an indirect-response turnover model of the Psoriasis Area and",
    "Severity Index (PASI) in which adalimumab inhibits lesion formation",
    "through an Imax relationship, kin = rbase * kout with Imax fixed at 1",
    "and IC50 0.95 ug/mL; PASI starts at its baseline 14.3 and is lost at",
    "kout 0.04 /day, a lesion turnover half-life of 17.3 days.",
    "The PK layer was fit to all 543 patients and the PD layer to the 367",
    "with baseline PASI >= 10; the two are published as one coupled PK-PD",
    "system (Figure 1 and the Data S1 Code S2 $DES), which is how they are",
    "encoded here."
  )
  reference <- paste(
    "Pan S, Tsakok T, Wei R, Dand N, Loeff FC, Bloem K, de Vries A,",
    "Baudry D, Duckworth M, Pushpa-Rajah A, Russell A, Alsharqi A,",
    "Becher G, Murphy R, Wahie S, Wright A, Griffiths CEM, Reynolds NJ,",
    "Barker J, Warren RB, Burden AD, Rispens T, Mahil SK, Standing JF,",
    "Smith CH; BADBIR Study Group, BSTOP Study Group, PSORT Consortium.",
    "Evaluation of a Therapeutic Drug Monitoring Strategy for Adalimumab",
    "in Psoriasis: A Prospective Pharmacokinetic-Pharmacodynamic Study.",
    "Clin Transl Sci. 2026;19(4):e70563. doi:10.1111/cts.70563.",
    "PK parameter values from Table 2; PD parameter values from Table 3;",
    "structural equations and the covariate model from the Data S1",
    "supplement (Code S1 for PK, Code S2 for PK-PD).",
    sep = " "
  )
  vignette <- "Pan_2026_adalimumab"

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on both CL/F and V/F, normalised to 70 kg, with",
        "both exponents fixed a priori (0.75 on CL, 1 on V) rather than",
        "estimated (Pan 2026 Methods 2.2; Table 2 footnote b). Allowing the",
        "CL exponent to be estimated returned 0.67 (95% CI 0.35-1.00), which",
        "the authors judged a negligible change and so retained the fixed",
        "0.75 (Results 3.2.1). Note that 70 kg is the allometric reference,",
        "NOT the cohort median, which was 88.7 kg (Table 1); a median-weight",
        "patient therefore has CL/F about 19% above the reference value.",
        "2.1% missing in the source cohort."
      ),
      source_name        = "WEIGHT"
    ),
    SEXF = list(
      description        = "1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "The source dataset codes the opposite orientation: its GENDER column",
        "is 1 = male / 0 = female, and the Code S1 $PK block makes GENDER == 1",
        "the reference ('Most common') with the coefficient applied to",
        "GENDER == 0. The cohort was 63.2% male (Table 1), which confirms",
        "that GENDER == 1 is indeed male and therefore the most common level.",
        "Encoded here as SEXF = 1 - GENDER so the canonical orientation is",
        "preserved; the coefficient is unchanged and is applied to SEXF == 1,",
        "matching Table 2's own label 'coeff female on CL'. Females had the",
        "HIGHER clearance in this psoriasis cohort (+21.6%), which is the",
        "opposite direction to the rheumatoid-arthritis study of Ternant 2015",
        "and is flagged as such in Pan 2026 Discussion 4.2."
      ),
      source_name        = "GENDER"
    ),
    WAIST = list(
      description        = "Waist circumference measured at baseline.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL/F normalised to 101 cm: (WAIST/101)^0.888",
        "(Code S1 CLWAIST block; Table 2 'coeff waist on CL' = 0.888).",
        "101 cm is the cohort median (Table 1, range 46-165 cm), so this is a",
        "median-centred normalisation. The strongest single covariate on",
        "clearance after body weight. 8.1% of the source cohort had no waist",
        "measurement; the control stream codes those as -99 and sets the",
        "whole CLWAIST term to 1, which is exactly equivalent to supplying",
        "WAIST = 101 (the reference). Supply 101 for missing values rather",
        "than a sentinel -- this model has no -99 branch, and a literal -99",
        "would raise a negative number to a fractional power."
      ),
      source_name        = "WAIST"
    ),
    CONC_ADA_AUML = list(
      description        = paste(
        "Anti-drug-antibody level against adalimumab, measured by",
        "antigen-binding radioimmunoassay and reported on the assay's",
        "arbitrary-unit concentration scale."
      ),
      units              = "AU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL/F normalised to 76.03 AU/mL:",
        "(CONC_ADA_AUML/76.03)^0.368 (Code S1 CLARIA block; Table 2 'coeff",
        "ADA on CL' = 0.368). This is a quantitative concentration, NOT a",
        "positivity flag and NOT a dilution titer: the CLARIA block carries",
        "no IF branch at all -- unlike the CLWAIST block in the same stream,",
        "which does have a missing-value branch -- so every subject must",
        "carry a strictly positive value, and a 0 would make the power form",
        "undefined. Positivity, where Pan 2026 reports it, is a derived",
        "dichotomy at the assay cutoff of 30 AU/mL (Methods 2.1.3); do not",
        "dichotomise this column onto ADA_POS, which would discard an",
        "estimated parameter. The normalising constant 76.03 is not printed",
        "anywhere in the article or supplement; it is transcribed verbatim",
        "from the control stream (see the vignette Errata -- it is presumed",
        "to be the cohort median, matching the WAIST/101 median-centring",
        "convention used in the same $PK block, but this is not stated).",
        "Assay units are specific to this radioimmunoassay and are NOT",
        "convertible to ng/mL without the assay's own standard."
      ),
      source_name        = "ARIA"
    ),
    DIS_HYPERT = list(
      description        = "1 = hypertension recorded in the medical history, 0 = not recorded.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hypertension)",
      notes              = paste(
        "Proportional increase in CL/F of 17.7% (Code S1 CLHYPE block, which",
        "makes HYPE == 0 the reference 'Most common'; Table 2 'coeff",
        "hypertension on CL' = 0.177). 23.2% of the cohort was hypertensive",
        "(Table 1, 'HT'). Pan 2026 Discussion 4.2 reads this as likely",
        "reflecting underlying obesity rather than a hypertension mechanism",
        "per se, noting that it survived backward elimination alongside both",
        "body weight and waist circumference."
      ),
      source_name        = "HYPE"
    )
  )

  # Covariates that Pan 2026 screened in the stepwise covariate model but did
  # not retain in the final model. Listed for provenance only; none is
  # referenced in model(). The screened set is given in Methods 2.2.1 and the
  # $INPUT line of Code S1 (most are marked =DROP in the final PK run).
  #
  # Only the screened covariates that map onto an EXISTING canonical column
  # are listed here. The remainder of the screened set has no canonical entry,
  # and minting register names for covariates that no model actually uses
  # would pollute the register, so they are recorded in prose instead:
  # alcohol consumption (81.2% yes, 2.5% missing; Code S1 ALCO=DROP -- the
  # registered ALCOHOL_ABUSE is chronic abuse, a different concept),
  # dyslipidaemia (8.6%, DYSL=DROP), asthma (10.6%, ASTH=DROP -- the
  # registered DIS_SASTHMA is specifically moderate-to-severe asthma),
  # major depressive disorder (20.3%, DEPR=DROP), liver disease (8.7%,
  # LIVE=DROP), inflammatory/psoriatic arthritis (23.6%, 2.3% missing,
  # ANY_IA=DROP), psoriasis disease duration (median 21 years, range 2-63,
  # 6.1% missing), palmoplantar involvement (17.2%, 4.5% missing,
  # PALMS=DROP), and ethnicity (89.5% white, ETHNICITY=DROP).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline", units = "year", type = "continuous",
      notes = "Screened (Methods 2.2.1 demographics); not retained. Median 44.3 years, range 17.4-80.4 (Table 1). Marked AGE=DROP in the Code S1 $INPUT."
    ),
    BMI = list(
      description = "Body mass index at baseline", units = "kg/m^2", type = "continuous",
      notes = "Screened (Methods 2.2.1 anthropometrics); not retained. Median 29.7, range 16.6-63, 13.4% missing (Table 1). Marked BMI=DROP in the Code S1 $INPUT. Correlated with WAIST, but the stepwise covariate model retained waist circumference on CL and rejected BMI -- which is why this model carries WAIST and not BMI."
    ),
    HT = list(
      description = "Body height at baseline", units = "cm", type = "continuous",
      notes = paste(
        "Read into the Code S1 $INPUT as HEIGHT but not used in any $PK",
        "expression; not retained. NOTE THE TOKEN COLLISION: this is the",
        "canonical column for body HEIGHT, whereas Pan 2026's Table 1 uses",
        "'HT' as its abbreviation for HYPERTENSION. The hypertension",
        "covariate, which IS retained, is DIS_HYPERT above; the source",
        "dataset calls it HYPE, not HT."
      )
    ),
    CREAT = list(
      description = "Serum creatinine at baseline", units = "umol/L", type = "continuous",
      notes = "Screened (Methods 2.2.1); not retained. Median 76 umol/L, range 42-149, 2.7% missing (Table 1). Marked CREATININE=DROP in the Code S1 $INPUT. The source quantity is serum creatinine, not a computed clearance, hence CREAT rather than CRCL."
    ),
    SMOKE = list(
      description = "Current smoker", units = "(binary)", type = "binary",
      notes = "Screened (Methods 2.2.1 lifestyle); not retained. 56.3% smokers, 2.5% missing (Table 1)."
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus comorbidity", units = "(binary)", type = "binary",
      notes = "Screened (Methods 2.2.1 comorbidities); not retained. 17% of cohort (Table 1, 'DM'). Marked DIAB=DROP in the Code S1 $INPUT."
    ),
    PRIOR_BIO = list(
      description = "Prior biologic exposure before adalimumab initiation", units = "(binary)", type = "binary",
      notes = "Screened; not retained. Table 1 reports the complement, 69.2% biologic-naive, so PRIOR_BIO = 1 for the remaining 30.8%. Marked BIO_NAIVE=DROP in the Code S1 $INPUT."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "adalimumab",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "adalimumab",
      units    = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    pasi = list(
      analyte  = "none",
      units    = "PASI units (0-72 clinical score)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 543L,
    n_studies      = 2L,
    age_range      = "17.4-80.4 years",
    age_median     = "44.3 years",
    weight_range   = "42.6-170 kg",
    weight_median  = "88.7 kg",
    sex_female_pct = 36.8,
    race_ethnicity = c(White = 89.5, Other = 10.5),
    disease_state  = "moderate-to-severe chronic plaque psoriasis",
    dose_range     = paste(
      "Adalimumab 40 mg subcutaneously. Where administration dates were",
      "missing the standard label regimen was assumed: 80 mg loading dose,",
      "then 40 mg every 2 weeks starting 1 week after loading, with full",
      "adherence. No patient changed dose or dosing interval during the",
      "first treatment year."
    ),
    regions        = "United Kingdom (60 participating centers)",
    notes          = paste(
      "Baseline demographics from Pan 2026 Table 1 (n = 544 met the",
      "inclusion criteria; 543 provided 946 PK samples and 539 provided",
      "1700 PASI measurements within the first treatment year).",
      "Real-world observational data pooled from two prospective UK studies",
      "under the PSORT consortium: BSTOP (Biomarkers of Systemic Treatment",
      "Outcomes in Psoriasis, 79 centers) and PSORT-Discovery (11 BSTOP",
      "centers), with eligible adults drawn from the BADBIR registry.",
      "Samples were collected June 2009 to December 2016 during routine",
      "care WITHOUT regard to time since dose, so the dataset contains both",
      "trough and non-trough samples; 56 of 946 samples (5.9%) were below",
      "the 0.01 ug/mL assay limit and were replaced by half that limit.",
      "The PK model used all 543 patients (Table 2). The PD model was",
      "restricted to the 367 patients with baseline PASI >= 10, the",
      "threshold for biologic eligibility (Table 3); PD parameters refit to",
      "all 539 PASI-evaluable patients are given in Table S1 and are",
      "reproduced in the vignette as a sensitivity comparison.",
      "Missing covariates were imputed with medians or modal categories",
      "(Discussion 4.5). NONMEM 7.6, FOCE with INTERACTION; diagnostics in",
      "R 4.4.1."
    )
  )

  ini({
    # ================================================================
    # Pharmacokinetics. Final estimates from Pan 2026 Table 2 (n = 543).
    #
    # NOTE ON PROVENANCE: the Data S1 Code S1 control stream is printed
    # with its INITIAL estimates, which differ from the published final
    # estimates (e.g. $THETA(7) CLWAIST1 = 0.771436 initial vs 0.888
    # final in Table 2). Structure and the covariate functional forms
    # are therefore taken from Code S1, but every VALUE below comes
    # from Table 2.
    #
    # F1 was fixed to 1 in Code S1, so clearance and volume are the
    # APPARENT quantities CL/F and V/F. Adalimumab's absolute
    # subcutaneous bioavailability of about 64% (Introduction) is NOT
    # separated out by this model and must not be applied on top.
    # ================================================================
    lka <- log(0.268)
    label("First-order absorption rate constant (1/day)")
    # Table 2: ka = 0.268 /day (RSE 11.7%). Code S1 KA = THETA(1).

    lcl <- log(0.386)
    label("Apparent clearance CL/F at the reference covariate vector (L/day)")
    # Table 2: CL/F = 0.386 L/day (RSE 3.5%). Code S1 TVCL = THETA(2).
    # Reference vector is WT 70 kg, WAIST 101 cm, ADA 76.03 AU/mL,
    # male, no hypertension. A sensitivity analysis restricted to the
    # 349 samples with complete dosing records gave 0.358 L/day
    # (95% CI 0.319-0.397), overlapping this estimate (Results 3.2.1).

    lvc <- fixed(log(10.8))
    label("Apparent central volume V/F at 70 kg (L)")
    # Table 2: V/F = 10.8 L [fix], footnote a "Fixed on the value from
    # Ternant et al. [30]". Code S1 $THETA "10.8 FIX". Fixed a priori
    # because most samples were at steady state, which limits
    # estimation in the presence of flip-flop kinetics (Discussion 4.2).

    # ---- Covariate effects on CL/F and V/F (Code S1 $PK) ----
    # TVCL = THETA(2) * (WEIGHT/70)^0.75
    #        * (ARIA/76.03)^THETA(4)
    #        * (1 + THETA(5))  if female
    #        * (1 + THETA(6))  if hypertensive
    #        * (WAIST/101)^THETA(7)
    # TVV  = THETA(3) * (WEIGHT/70)^1

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")
    # Table 2 "coeff weight on CL" = 0.75 [fix], footnote b "Fixed a
    # priori using allometric scaling [25]". Hardcoded as 0.75 in the
    # Code S1 TVCL line rather than carried as a THETA.

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V/F (unitless)")
    # Table 2 "coeff weight on V" = 1 [fix], footnote b. Hardcoded as
    # the exponent 1 in the Code S1 TVV line.

    e_conc_ada_cl <- 0.368
    label("Power exponent of ADA level on CL/F ((CONC_ADA_AUML/76.03)^e_conc_ada_cl)")
    # Table 2 "coeff ADA on CL" = 0.368 (RSE 6.6%). Code S1 THETA(4),
    # CLARIA = (ARIA/76.03)**THETA(4). Table S3 confirms 0.368 for the
    # final model, against 0.366 and 0.382 for the two two-compartment
    # sensitivity fits.

    e_sexf_cl <- 0.216
    label("Proportional increase in CL/F for female sex (fraction)")
    # Table 2 "coeff female on CL" = 0.216 (RSE 25.1%). Code S1
    # THETA(5), CLGENDER = (1 + THETA(5)) when GENDER == 0 (female);
    # GENDER == 1 (male) is the reference. See covariateData[[SEXF]] --
    # this model uses SEXF = 1 - GENDER, so the coefficient multiplies
    # SEXF directly with no sign change.

    e_waist_cl <- 0.888
    label("Power exponent of waist circumference on CL/F ((WAIST/101)^e_waist_cl)")
    # Table 2 "coeff waist on CL" = 0.888 (RSE 18.2%). Code S1 THETA(7),
    # CLWAIST = (WAIST/101)**THETA(7).

    e_dis_hypert_cl <- 0.177
    label("Proportional increase in CL/F for hypertension (fraction)")
    # Table 2 "coeff hypertension on CL" = 0.177 (RSE 35.1%). Code S1
    # THETA(6), CLHYPE = (1 + THETA(6)) when HYPE == 1.

    # ---- PK between-subject variability ----
    # Table 2 reports BSV as a percentage; Code S1 $OMEGA carries
    # variances on the log scale (exponential ETAs, CL = TVCL*EXP(ETA1)).
    # Converted with the log-normal identity omega^2 = log(CV^2 + 1).
    # See the vignette Errata for why that identity, rather than
    # omega = CV, is the reading taken here.
    etalcl ~ 0.102774
    # Table 2 "BSV on CL (%)" = 32.9 (RSE 17.9%); log(0.329^2 + 1).

    etalvc ~ 0.463623
    # Table 2 "BSV on V (%)" = 76.8 (RSE 16.5%); log(0.768^2 + 1).

    # ---- PK residual error ----
    # Code S1 $ERROR: Y = IPRED*(1 + EPS(1)) + EPS(2), i.e. combined
    # proportional and additive on the linear scale.
    propSd <- 0.195
    label("Proportional residual error for serum adalimumab (fraction)")
    # Table 2 "Proportional error (%)" = 19.5 (RSE 25.9%).

    addSd <- 1.78
    label("Additive residual error for serum adalimumab (ug/mL)")
    # Table 2 "Additive error (SD)" = 1.78 (RSE 8.8%), reported directly
    # as an SD, so it is used as-is rather than square-rooted.

    # ================================================================
    # Pharmacodynamics: indirect-response PASI turnover with Imax
    # inhibition of lesion formation. Final estimates from Pan 2026
    # Table 3 (n = 367, the baseline PASI >= 10 subset). Structure from
    # Data S1 Code S2 $PK / $DES.
    # ================================================================
    lrbase <- log(14.3)
    label("Baseline PASI score (PASI units)")
    # Table 3 "Baseline PASI" = 14.3 (RSE 4.6%). Code S2 BSL = THETA(1),
    # used both as the PASI initial condition A_0(3) and to pin
    # KIN = BSL*KOUT so the drug-free system sits exactly at baseline.

    lkout <- log(0.04)
    label("First-order rate constant of psoriatic lesion loss (1/day)")
    # Table 3 "kout (/day)" = 0.04 (RSE 7.7%). Code S2 TVKOUT = THETA(2).
    # Consistency check: log(2)/0.04 = 17.33 days, matching the lesion
    # turnover half-life of "17.3 days" quoted in Results 3.2.2.

    limax <- fixed(log(1))
    label("Maximum fractional inhibition of lesion formation (unitless)")
    # Table 3 "Emax" = 1 [fix]. Code S2 $THETA "1 FIX". With Imax fixed
    # at 1 the drug can suppress lesion formation completely at
    # saturating concentrations.

    lic50 <- log(0.95)
    label("Serum concentration at half-maximal inhibition (ug/mL)")
    # Table 3 "EC50 (ug/mL)" = 0.95 (RSE 13.7%). Code S2 TVIC50 =
    # THETA(4) and Table S4 both name this IC50; the quantity is the
    # half-maximal INHIBITORY concentration in the Imax expression
    # DG = EMAX*CONC/(IC50+CONC), hence lic50 here. Table 3's "EC50"
    # and Table S4's "IC50" are the same parameter -- both are glossed
    # "concentration at 50% of maximum inhibition on TNF-alpha".

    # ---- PD between-subject variability ----
    # Table 3 percentages, converted with omega^2 = log(CV^2 + 1) as
    # for the PK etas. All three PD ETAs are exponential in Code S2.
    etalrbase ~ 0.138889
    # Table 3 "BSV on baseline (%)" = 38.6 (RSE 16.0%); log(0.386^2 + 1).

    etalkout ~ 0.858473
    # Table 3 "BSV on kout (%)" = 116.6 (RSE 19.7%); log(1.166^2 + 1).

    etalic50 ~ 0.664151
    # Table 3 "BSV on EC50 (%)" = 97.1 (RSE 22.0%); log(0.971^2 + 1).
    # Pan 2026 found no covariate explaining the large kout or IC50
    # variability, and a two-subgroup mixture on IC50 did not improve
    # the fit (Discussion 4.3).

    # ---- PD residual error ----
    # Code S2 $ERROR: Y = IPRED + EPS(1), additive only, on the PASI
    # score itself (IPRED = A(3), the PASI compartment).
    addSd_pasi <- 3.2
    label("Additive residual error for PASI score (PASI units)")
    # Table 3 "Additive error (SD)" = 3.2 (RSE 7.8%), reported as an SD.
  })

  model({
    # ------------------------------------------------------------
    # 1. Individual PK parameters (Code S1 $PK).
    #    The five covariate terms multiply the typical clearance in
    #    the order Code S1 assembles them:
    #      CLCOV = CLARIA * CLGENDER * CLHYPE * CLWAIST
    #      TVCL  = THETA(2) * (WEIGHT/70)^0.75 ; TVCL = CLCOV * TVCL
    #    Missing WAIST is handled by supplying the reference 101 cm,
    #    which reproduces the stream's CLWAIST = 1 branch exactly.
    # ------------------------------------------------------------
    ka <- exp(lka)

    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (CONC_ADA_AUML / 76.03)^e_conc_ada_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_dis_hypert_cl * DIS_HYPERT) *
      (WAIST / 101)^e_waist_cl

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # 2. Micro-constant. Code S2 KE = CL/V.
    kel <- cl / vc

    # ------------------------------------------------------------
    # 3. Individual PD parameters (Code S2 $PK). kin is not estimated:
    #    it is pinned to KIN = BSL*KOUT so that, with no drug on
    #    board, d/dt(pasi) = rbase*kout - kout*pasi is exactly zero at
    #    pasi = rbase and the untreated PASI trajectory is flat. The
    #    paper has no placebo arm, so there is no placebo-decay term
    #    (Discussion 4.5: "The lack of a placebo arm precluded
    #    assessment of the untreated disease trajectory").
    # ------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout + etalkout)
    imax  <- exp(limax)
    ic50  <- exp(lic50 + etalic50)
    kin   <- rbase * kout

    # ------------------------------------------------------------
    # 4. ODE system. Code S2 $MODEL declares COMP=(DOSE), COMP=(CENTRAL)
    #    and COMP=(PASI) in that order, and $DES gives:
    #      DADT(1) = -KA*A(1)
    #      DADT(2) =  KA*A(1) - KE*A(2)
    #      CONC    =  A(2)/V
    #      DG      =  EMAX*CONC/(IC50+CONC)
    #      DADT(3) =  KIN*(1-DG) - KOUT*A(3)
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    # Fractional inhibition of lesion formation by adalimumab.
    dg <- imax * Cc / (ic50 + Cc)

    d/dt(pasi) <- kin * (1 - dg) - kout * pasi

    # Code S2 A_0(3) = BSL: the PASI state starts at the individual
    # baseline rather than at zero.
    pasi(0) <- rbase

    # ------------------------------------------------------------
    # 5. Observations. No f(depot): Code S1 and Code S2 both set
    #    F1 = 1, so bioavailability is folded into CL/F and V/F.
    # ------------------------------------------------------------
    Cc   ~ add(addSd) + prop(propSd)
    pasi ~ add(addSd_pasi)
  })
}
