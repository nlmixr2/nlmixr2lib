Nava_2018_busulfan <- function() {
  description <- paste(
    "One-compartment IV population PK model with first-order elimination for",
    "intravenous busulfan in children and adolescents (0.1-20 years)",
    "undergoing conditioning before haematopoietic stem cell transplantation.",
    "Theory-based allometric scaling of actual body weight (exponent 0.75 on",
    "CL, 1 on V) referenced to 70 kg; a sigmoid Emax maturation function of",
    "postmenstrual age on clearance carried over unchanged from McCune 2014",
    "(TM50 = 46 weeks, Hill = 2.3); a power effect of postmenstrual age on",
    "volume; and a three-level GSTA1 diplotype effect on clearance (rapid",
    "metabolizers 7% faster, poor metabolizers 12% slower than normal",
    "metabolizers). Between-subject variability on CL and V, between-occasion",
    "variability on CL and V across the four days of busulfan dosing, and a",
    "proportional residual error. This is the first paediatric busulfan popPK",
    "model to incorporate a pharmacogenetic covariate (Nava 2018)."
  )
  reference <- paste(
    "Nava T, Kassir N, Rezgui MA, Uppugunduri CRS, Huezo-Diaz Curtis P,",
    "Duval M, Theoret Y, Daudt LE, Litalien C, Ansari M, Krajinovic M,",
    "Bittencourt H. Incorporation of GSTA1 genetic variations into a",
    "population pharmacokinetic model for IV busulfan in paediatric",
    "hematopoietic stem cell transplantation. Br J Clin Pharmacol.",
    "2018;84(7):1494-1504. doi:10.1111/bcp.13566.",
    "The maturation function is Equation D of the paper's Supplemental",
    "Material Table S2, which reproduces McCune JS, Bemer MJ, Barrett JS,",
    "Scott Baker K, Gamis AS, Holford NH. Clin Cancer Res 2014;20(3):754-763.",
    "doi:10.1158/1078-0432.CCR-13-1960."
  )
  vignette <- "Nava_2018_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the paper: "Plasma concentrations of Bu
  # were determined using a high-performance liquid chromatographic (HPLC)
  # assay with ultraviolet detection" (Methods, "PK analysis and genotyping"),
  # and busulfan is given as a constant-rate IV infusion straight into the
  # systemic circulation (Methods, "Population pharmacokinetic analysis").
  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Theory-based ('theoretical') allometric scaling with reference weight",
        "70 kg and exponents fixed at 0.75 on CL and 1 on V (Nava 2018 Results",
        "'Base model' and the displayed clearance equation on p. 1497).",
        "Treated as baseline over the 4-day busulfan course. The paper tested",
        "body surface area and a normal-fat-mass size descriptor (free fat mass",
        "plus 51% of fat mass for CL, 20% for V, per McCune 2014) against actual",
        "body weight and reports 'no advantage in comparison to ABW'",
        "(Discussion), so actual body weight is the descriptor in the final",
        "model. Nava 2018 Table 1 reports no weight summary; only weight",
        "adequacy (12.5% overweight, 6.3% obese among children older than 2)."
      ),
      source_name        = "ABW"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried in WEEKS, not the register-default months, because the source",
        "equations are written in weeks. Supplemental Material Table S2",
        "equation E defines PMA[weeks] = 52 * AGE[years] + 40; gestational age",
        "is set to 40 weeks for every subject because birth gestational age was",
        "not recorded (Methods 'Anthropometric data descriptors' and Discussion:",
        "'The PMA is calculated from the first day of the last menstruation",
        "period of the mother prior to the birth, which is set to 40 weeks due",
        "to missing information'). PAGE enters the model twice: (1) the sigmoid",
        "Emax maturation function on CL, Fmat = 1 / (1 + (PMA/46)^-2.3)",
        "(Supplemental Table S2 equation D; TM50 = 46 weeks and Hill = 2.3 both",
        "carried unchanged from McCune 2014); and (2) a power effect on V with",
        "exponent -0.06 (Table 2, 'PMA on V'). Cohort range 45.2-1080 weeks",
        "(age 0.1-20 years by equation E)."
      ),
      source_name        = "PMA"
    ),
    GSTA1_RM = list(
      description        = "GSTA1 rapid-metabolizer diplotype group (group G1) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer group G2, or poor metabolizer group G3)",
      notes              = paste(
        "1 = subject's GSTA1 diplotype places them in group G1 (rapid",
        "metabolizers, higher predicted promoter activity). Paired with",
        "GSTA1_PM so that both indicators = 0 selects the G2 normal-metabolizer",
        "reference group. Diplotypes were built from four promoter loci -- -69",
        "(rs3957357), -513 (rs11964968), -631 (rs4715333) and -1142",
        "(rs58912740) -- and grouped by predicted gene expression (Methods 'PK",
        "analysis and genotyping'; Supplemental Material Table S1 and Figures",
        "S1-S2). 16 of 112 subjects (14.3%) were G1 (Nava 2018 Table 1)."
      ),
      source_name        = "GSTA1 group G1"
    ),
    GSTA1_PM = list(
      description        = "GSTA1 poor-metabolizer diplotype group (group G3) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer group G2, or rapid metabolizer group G1)",
      notes              = paste(
        "1 = subject's GSTA1 diplotype places them in group G3 (poor",
        "metabolizers, lower predicted promoter activity). Paired with",
        "GSTA1_RM so that both indicators = 0 selects the G2 normal-metabolizer",
        "reference group. 16 of 112 subjects (14.3%) were G3 (Nava 2018",
        "Table 1)."
      ),
      source_name        = "GSTA1 group G3"
    ),
    OCC = list(
      description        = "Busulfan dosing-day occasion index",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer 1-4 identifying which of the four consecutive busulfan dosing",
        "days an observation belongs to. Decomposed inside model() into binary",
        "indicators multiplying per-occasion eta slots. The paper reports",
        "between-occasion variability (its 'BOV') on both CL and V (Table 2)",
        "but does not state the occasion count; four is used because busulfan",
        "was 'administered intravenously (IV) for four consecutive days'",
        "(Methods 'Treatment regimen') and the once-daily cohort contributed a",
        "mean of 3.7 PK profiles per conditioning regimen (Results), i.e. up to",
        "one profile per dosing day. One shared BOV magnitude is used across",
        "all four occasions, so occasions 2-4 are fixed to the occasion-1",
        "variance."
      ),
      source_name        = "occasion"
    )
  )

  # Screened during covariate analysis but not retained in the final model.
  # Documented here so the paper's covariate screen is preserved without
  # carrying "declared but not referenced" convention warnings.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested as a potential covariate (Methods 'Covariate analysis and",
        "sources of variability') and not retained by the stepwise forward",
        "(p = 0.05) / backward (p = 0.01) procedure. 59 of 112 subjects",
        "(52.7%) were female (Table 1). Sex is also an input to the free-fat-",
        "mass equation (Supplemental Table S2 equation A) of the normal-fat-",
        "mass size descriptor, which was itself rejected in favour of actual",
        "body weight."
      )
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Post-natal age was tested directly as a covariate and lost to",
        "postmenstrual age: the maturation factor Fmat 'was a more significant",
        "covariate on clearance than post-natal age' (Discussion). Age still",
        "reaches the model, but only through PAGE (= 52 * AGE + 40 weeks)."
      )
    ),
    DIS_MALIGNANT = list(
      description = "Malignant baseline disease indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Baseline disease (malignant vs non-malignant) was tested and not",
        "retained (Methods 'Covariate analysis and sources of variability').",
        "74 of 112 subjects (66.1%) had a malignant diagnosis (Table 1). Name",
        "is descriptive only -- no canonical register entry was created,",
        "because the covariate never enters model()."
      )
    ),
    REGIMEN_BUCY = list(
      description = "BuCy conditioning-regimen indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Conditioning regimen (BuCy vs others) was tested and not retained",
        "(Methods 'Covariate analysis and sources of variability'). 91 of 112",
        "subjects (81.3%) received BuCy (Table 1). Name is descriptive only --",
        "no canonical register entry was created, because the covariate never",
        "enters model()."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112L,
    n_studies      = 1L,
    age_range      = "0.1-20 years",
    age_median     = "5.4 years",
    sex_female_pct = 52.7,
    race_ethnicity = c(Caucasian = 83.9, African = 10.7, Other = 5.4),
    disease_state  = paste(
      "Children and adolescents receiving intravenous busulfan as part of the",
      "conditioning regimen before autologous or allogeneic haematopoietic",
      "stem cell transplantation. 66.1% malignant (AML 30.4%, MDS 17.0%,",
      "neuroblastoma 10.7%, ALL 7.1%, MPS 0.9%) and 33.9% non-malignant",
      "(haemoglobinopathy 14.3%, immunodeficiency 8.9%, metabolic disease",
      "6.3%, haemophagocytic syndrome 4.5%). Conditioning regimens: BuCy",
      "81.3%, BuMel 12.5%, BuCyVP16 7.1%, BuCyMel 0.9%, BuMelAraC 0.9%."
    ),
    dose_range     = paste(
      "Intravenous busulfan for four consecutive days on one of two schedules.",
      "Bu6 (April 2002 - April 2012): 2-h infusion every 6 h, initial dose",
      "16 mg/m^2 in infants <= 3 months, 0.8 mg/kg in children > 3 months and",
      "< 1 year or >= 4 years, and 1.0 mg/kg in children >= 1 and < 4 years.",
      "Bu24 (from May 2012): 3-h infusion every 24 h at four times the",
      "corresponding Bu6 dose. PK-guided dose adjustment was allowed from the",
      "2nd (Bu24) or 5th (Bu6) dose, targeting a cumulative 4-day AUC of",
      "18 000 uM.min (i.e. AUC24h 3600-6000 uM.min or AUC6h 900-1500 uM.min)."
    ),
    regions        = "Canada (single centre: CHU Sainte-Justine, Montreal, Quebec)",
    notes          = paste(
      "Retrospective chart review of paediatric transplants performed between",
      "April 2002 and August 2016; part of NCT01257854. 199 PK profiles and",
      "1735 plasma busulfan concentrations from 112 subjects; 114 profiles",
      "(57%) came from 31 once-daily conditioning regimens and the remainder",
      "from 81 four-times-daily regimens, where 99% of profiles were",
      "first-dose. Sampling: pre-dose and 120, 135, 150, 180, 240, 300, 360 min",
      "after infusion start for Bu6; 180, 195, 240, 300, 360, 480 min for Bu24",
      "(median 9 samples per Bu6 patient, range 7-18; median 28 per Bu24",
      "patient, range 14-28). Assay: HPLC with ultraviolet detection.",
      "Estimation: Phoenix NLME 6.4 (Certara), FOCE with interaction.",
      "GSTA1 diplotype groups: G1 14.3%, G2 71.4%, G3 14.3% (Table 1).",
      "Regimens containing another chemotherapeutic drug before or during the",
      "busulfan infusion days were excluded so that fludarabine-associated",
      "clearance changes could not bias the fit (Discussion)."
    )
  )

  ini({
    # ---- Structural parameters, reported for a 70 kg patient ---------------
    # Table 2 labels these "CL 70 kg" and "V 70 kg", and Results states the
    # final-model estimates are given "for a patient with a median body weight
    # of 70 kg". Because the maturation factor Fmat -> 1 at adult PMA and the
    # GSTA1 factor is 1 in the reference G2 group, 13.70 L/h is the clearance
    # of a fully mature 70 kg G2 subject; the printed clearance equation on
    # p. 1497 makes this explicit.
    lcl <- log(13.70); label("Typical busulfan clearance for a mature 70 kg normal-metabolizer (L/h)")  # Nava 2018 Table 2 final model (RSE 2.43%; bootstrap median 13.69, 95% CI 13.12-14.16); reproduced in the p. 1497 clearance equation
    lvc <- log(49.57); label("Typical busulfan volume of distribution at 70 kg and the reference PMA (L)")  # Nava 2018 Table 2 final model (RSE 1.15%; bootstrap median 49.57)

    # ---- Allometric exponents ---------------------------------------------
    # "Theoretical allometric scaling of the actual body weight (ABW) [was]
    # included in the base model" (Results, "Base model"). Theory-based
    # allometry fixes 0.75 on clearances and 1 on volumes; the 0.75 on CL is
    # printed directly in the p. 1497 clearance equation, CL = 13.7 *
    # (ABW/70)^0.75 * Fmat * F_GSTA1, and both are printed in Supplemental
    # Table S2 equation C, whose row label reads "Fsize, where pwr = 1 for V
    # and 0.75 for CL" over Fsize = (NFM/70)^pwr. (Nava substitutes ABW for
    # McCune's NFM -- the normal-fat-mass descriptor was tested and rejected,
    # Discussion -- but keeps the same theory-based exponents.) Neither
    # exponent is given an RSE in Table 2 (they do not appear there at all),
    # which is the paper's way of saying they were not estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")  # Nava 2018 p. 1497 clearance equation; Supplemental Table S2 eq. C row label ("0.75 for CL")
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on V (unitless)")   # Nava 2018 Supplemental Table S2 eq. C row label ("pwr = 1 for V"); Results 'Base model' ("theoretical allometric scaling")

    # ---- Maturation of busulfan metabolism on postmenstrual age -----------
    # Supplemental Material Table S2 ("Equations used in McCune's models"),
    # whose row label prints both constants and whose equation D prints the
    # form, verbatim:
    #     "Fmat, where TM50=46; Hill's coefficient=2.3"
    #     Fmat = 1 / (1 + (PMA / TM50)^-Hill)                             (D)
    # Both constants are carried unchanged from McCune 2014 -- the paper
    # neither re-estimates them nor reports an RSE for them, so both are
    # fixed(). The main text independently states the 46-week value:
    # "an Emax sigmoid curve to explain the Bu maturation, which would reach
    # 50% of the adult metabolic rate at 46 weeks of PMA" (Discussion).
    tm50_mat <- fixed(46);  label("Postmenstrual age at 50% of mature busulfan metabolism (weeks)")  # Nava 2018 Supplemental Table S2 row label "TM50=46" + eq. D; restated in Discussion p. 1499
    hill_mat <- fixed(2.3); label("Hill coefficient of the busulfan-metabolism maturation function (unitless)")  # Nava 2018 Supplemental Table S2 row label "Hill's coefficient=2.3" + eq. D

    # ---- Postmenstrual-age effect on volume of distribution ---------------
    # Table 2 reports a single "PMA on V" coefficient of -0.06 (RSE 23.40%;
    # bootstrap median -0.06, 95% CI -0.18 to 0.09) but never prints the
    # equation it sits in. It is a POWER effect -- see the vignette's
    # "Assumptions and deviations" section for the full argument; in brief, an
    # exponential or linear form on PMA in weeks is arithmetically impossible
    # over the cohort's 45-1080 week range (exp(-0.06 * 1080) underflows; a
    # linear 1 - 0.06 * PMA goes negative above PMA = 17 weeks, below every
    # subject's PMA), whereas a dimensionless power exponent of -0.06 gives a
    # 13% higher V per kg in a neonate than at the cohort's central PMA, which
    # is the expected direction and magnitude for total body water.
    e_page_vc <- -0.06; label("Power exponent of postmenstrual age on V (unitless)")  # Nava 2018 Table 2 final model 'PMA on V' (RSE 23.40%)

    # Reference (centring) postmenstrual age for the power effect above. NOT
    # PRINTED IN THE PAPER; derived here from the paper's own numbers. Table 1
    # gives a cohort median age of 5.4 years, and Supplemental Table S2
    # equation E gives PMA[weeks] = 52 * AGE[years] + 40, so the cohort median
    # PMA is 52 * 5.4 + 40 = 320.8 weeks. That the centring is at the cohort's
    # central PMA (rather than at an adult PMA) is forced by Table 2 itself:
    # adding the PMA term moved V70 only from 50.05 (base model) to 49.57
    # (final model), a 1.0% shift. A centring at an adult PMA of ~2000 weeks
    # would have required V70 to fall by (320.8/2000)^0.06 = 11% to keep the
    # same cohort predictions, and it did not.
    pma_ref <- fixed(320.8); label("Reference postmenstrual age for the effect on V (weeks)")  # derived from Nava 2018 Table 1 median age 5.4 y via Supplemental Table S2 eq. E; centring value not printed in the source

    # ---- GSTA1 diplotype-group effect on clearance ------------------------
    # Table 2 ("GSTA1-group on CL (Reference G2)") and the sentence beneath
    # the p. 1497 clearance equation, "where F_GSTA1 = 1.07 for G1, 1 for G2
    # and 0.88 for G3 patients". These are the multiplicative factors
    # themselves, not log- or fractional-shift coefficients. G2 (normal
    # metabolizers) is the reference and carries a factor of exactly 1.
    e_gsta1_rm_cl <- 1.07; label("Multiplicative CL factor for GSTA1 rapid metabolizers (group G1) vs normal (G2)")  # Nava 2018 Table 2 final model G1 (RSE 4.40%; bootstrap median 1.07, 95% CI 1.05-1.21) and p. 1497
    e_gsta1_pm_cl <- 0.88; label("Multiplicative CL factor for GSTA1 poor metabolizers (group G3) vs normal (G2)")   # Nava 2018 Table 2 final model G3 (RSE 0.40%; bootstrap median 0.89, 95% CI 0.82-0.96) and p. 1497

    # ---- Between-subject variability --------------------------------------
    # Table 2 reports BSV as percentages under an exponential random-effect
    # model with log-normally distributed individual parameters (Methods,
    # "Population pharmacokinetic analysis"), i.e. as %CV. Converted to the
    # internal log-scale variance with omega^2 = log(CV^2 + 1):
    #   BSV CL 13.30% -> log(0.1330^2 + 1) = 0.0175561
    #   BSV V   7.00% -> log(0.0700^2 + 1) = 0.0048879
    # Both final-model values are independently confirmed in the Discussion:
    # "as much as 13% and 7% of the variability between subjects (BSV) on CL
    # and V, respectively, remained unexplained", and "the incorporation of
    # the GSTA1 diplotypes seemed to reduce the overall BSV of CL by 27% (from
    # 18.3% to 13.3%)". The Table 2 bootstrap column reports 18.80% for BSV CL,
    # which matches the BASE model (18.30%) rather than the final one; the
    # prose settles it in favour of 13.30%.
    etalcl ~ 0.0175561  # Nava 2018 Table 2 final model BSV CL 13.30% (RSE 8.80%, shrinkage 6.00%); confirmed in Discussion p. 1500
    etalvc ~ 0.0048879  # Nava 2018 Table 2 final model BSV V 7.00% (RSE 19.20%, shrinkage 24.60%); confirmed in Discussion p. 1500

    # ---- Between-occasion variability on the four busulfan dosing days ----
    # Table 2 "Between-occasion variability (BOV)": BOV CL 7.30%, BOV V 9.60%.
    # Same %CV -> log-scale variance conversion as the BSV terms:
    #   BOV CL 7.30% -> log(0.0730^2 + 1) = 0.0053246
    #   BOV V  9.60% -> log(0.0960^2 + 1) = 0.0091775
    # One shared magnitude per parameter across occasions, so occasions 2-4
    # are fixed to the occasion-1 variance (the nlmixr2 spelling of NONMEM's
    # $OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.0053246          # Nava 2018 Table 2 final model BOV CL 7.30% (RSE 8.10%, shrinkage 46.60%)
    etaiov_cl_2 ~ fixed(0.0053246)   # same magnitude, occasion 2
    etaiov_cl_3 ~ fixed(0.0053246)   # same magnitude, occasion 3
    etaiov_cl_4 ~ fixed(0.0053246)   # same magnitude, occasion 4
    etaiov_vc_1 ~ 0.0091775          # Nava 2018 Table 2 final model BOV V 9.60% (RSE 16.40%, shrinkage 9.30%)
    etaiov_vc_2 ~ fixed(0.0091775)   # same magnitude, occasion 2
    etaiov_vc_3 ~ fixed(0.0091775)   # same magnitude, occasion 3
    etaiov_vc_4 ~ fixed(0.0091775)   # same magnitude, occasion 4

    # ---- Residual error ----------------------------------------------------
    # Table 2 "Random residual variability -- sigma [Proportional error (%)]":
    # 7.40% in BOTH the base and the final model column (RSE 7.90%). The
    # bootstrap column of the same row reads 27.20% (95% CI 25.00-29.00), which
    # cannot be reconciled with the rest of the table -- see the vignette's
    # "Assumptions and deviations" section. The 7.40% value is used: it is what
    # the final-model Estimate column reports, and it is the only value
    # consistent with the reported shrinkage. With 9-28 samples per profile the
    # standard error of an individual CL is roughly sigma / sqrt(n); at
    # sigma = 7.4% that is ~2.5%, giving eta-shrinkage of about
    # 0.025^2 / (0.025^2 + 0.133^2) = 3%, close to the 6.0% Table 2 reports,
    # whereas sigma = 27.2% would imply ~31% shrinkage -- five times what was
    # observed.
    propSd <- 0.074; label("Proportional residual error (fraction)")  # Nava 2018 Table 2 final model (7.40%, RSE 7.90%)
  })

  model({
    # ---- 1. Derived covariate terms ---------------------------------------
    # Sigmoid Emax maturation of busulfan metabolism on postmenstrual age.
    # Supplemental Table S2 equation D: Fmat = 1 / (1 + (PMA / TM50)^-Hill).
    # PAGE is carried in weeks for this model (see covariateData$PAGE).
    fmat <- 1 / (1 + (PAGE / tm50_mat)^(-hill_mat))

    # GSTA1 diplotype-group factor on clearance. The paper states F_GSTA1 =
    # 1.07 for G1, 1 for G2 and 0.88 for G3, so the reference G2 group (both
    # indicators zero) contributes exactly 1.
    fgsta1 <- 1 +
      (e_gsta1_rm_cl - 1) * GSTA1_RM +
      (e_gsta1_pm_cl - 1) * GSTA1_PM

    # Power effect of postmenstrual age on volume of distribution, centred at
    # the cohort's median PMA. See the ini() comment on pma_ref for why the
    # power form and this centring are the only reading consistent with the
    # paper's own numbers.
    fpage_vc <- (PAGE / pma_ref)^e_page_vc

    # Between-occasion variability across the four busulfan dosing days.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 + oc4 * etaiov_vc_4

    # ---- 2. Individual parameters -----------------------------------------
    # Nava 2018 p. 1497: CL = 13.7 * (ABW/70)^0.75 * Fmat * F_GSTA1.
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * fmat * fgsta1
    vc <- exp(lvc + etalvc + iov_vc) * (WT / 70)^e_wt_vc * fpage_vc

    # ---- 3. Micro-constants -----------------------------------------------
    kel <- cl / vc

    # ---- 4. ODE system ----------------------------------------------------
    # One compartment with first-order elimination (Results, "Base model":
    # "Bu concentration-time data was best described by a one-compartment
    # model with a first-order elimination"). Busulfan is infused directly
    # into the systemic circulation, so there is no depot and no
    # bioavailability term.
    d/dt(central) <- -kel * central

    # ---- 5. Observation and error -----------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
