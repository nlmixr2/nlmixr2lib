Marianski_2026_ceftriaxone <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "intravenous ceftriaxone in critically ill children with multiple organ",
    "dysfunction syndrome (Marianski 2026; 44 patients aged 2 months to 17",
    "years, up to 15 whole-blood microsamples each over 3 days). Clearance is",
    "0.91 L/h and central volume 8.22 L at the 27.4 kg / 120 mL/min/1.73 m^2",
    "reference, with body weight entering clearance and central volume",
    "through fixed allometric exponents (0.75 on CL, 1 on V1) and estimated",
    "glomerular filtration rate entering clearance as a through-origin linear",
    "ratio (eGFR / 120). Intercompartmental clearance and peripheral volume",
    "carry no covariates. Interindividual variability is diagonal on all four",
    "disposition parameters with a proportional residual error. Concentrations",
    "are WHOLE BLOOD, not plasma. Conference abstract: the complete parameter",
    "set is published, but only as a raster table image.",
    sep = " "
  )
  reference <- paste(
    "Marianski S, Amajor V, Shiau J, Bwint A, Sharova A, Hall M, Rhodes NJ,",
    "Downes KJ, Scheetz MH; PALISI study investigators (2026). P-1254.",
    "Evaluation of Ceftriaxone Population Pharmacokinetics (PK) and",
    "Pharmacodynamics (PD) in Critically Ill Pediatric Population. Open Forum",
    "Infectious Diseases 13(Suppl 1):S818-S819, abstract citation ID",
    "ofaf695.1445. doi:10.1093/ofid/ofaf695.1445. PMC12792806. IDWeek 2025",
    "poster abstract, Session 148 (PK/PD Studies), 21 October 2025. All final",
    "estimates and the covariate parameterisation come from Table 1, which is",
    "published as an EMBEDDED RASTER IMAGE that no text extraction recovers;",
    "the values were read from the decoded image (page 1, first embedded",
    "975x1010 JPEG). Table 1's 'Model parameterized as:' block supplies the",
    "four covariate equations verbatim.",
    sep = " "
  )
  vignette <- "Marianski_2026_ceftriaxone"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: ceftriaxone was quantified in WHOLE BLOOD collected by
  # volumetric absorptive microsampling (VAMS, 20 uL/sample) using a validated
  # LC-MS/MS assay (Marianski 2026 Methods). This model was fitted to those
  # whole-blood concentrations, so Cc is a whole-blood concentration and is NOT
  # interchangeable with the plasma concentrations that most ceftriaxone popPK
  # models report. See population$protein_binding and the vignette Assumptions
  # and deviations before comparing this model against a plasma-based model or
  # against a plasma-referenced susceptibility breakpoint.
  compartmentData <- list(
    central     = list(analyte = "ceftriaxone", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "ceftriaxone", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance and central volume as the allometric ratio",
        "(WT / 27.4) raised to a FIXED exponent: 0.75 on CL and 1 on V1.",
        "Table 1's 'Model parameterized as:' block prints both exponents",
        "explicitly ('CL_i = CL_pop * (WT/27.4)^0.75 * (GFR/120)' and",
        "'V1_i = V1_pop * (WT/27.4)^1'), and Methods states the scaling was",
        "imposed rather than fitted ('Fixed allometric scaling of clearance",
        "(CL) and central volume of distribution (V1) defined the base",
        "models'). The exponents carry no standard error, %RSE or confidence",
        "interval anywhere in Table 1 -- the table's Fixed Effects block has",
        "rows only for V1, Q, V2 and CL -- which is the second, independent",
        "signal that they were fixed; they are therefore wrapped in fixed() in",
        "ini(). NOTE that weight deliberately does NOT enter Q or V2: Table 1",
        "prints 'Q_i = Q_pop' and 'V2_i = V2_pop' with no covariate term at",
        "all, which is unusual for an allometric model and is reproduced here",
        "faithfully. The abstract does NOT report the cohort weight",
        "distribution, so it is not possible to confirm that 27.4 kg is the",
        "cohort median; it is simply the normalisation constant the authors",
        "printed. See the vignette Assumptions and deviations.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalized, by the CKiD U25 equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BSA-normalized paediatric eGFR in mL/min/1.73 m^2. Marianski 2026",
        "Methods names the estimating equation: 'Chronic Kidney Disease in",
        "Children Under 25 eGFR or U25', i.e. the CKiD U25 serum-creatinine",
        "equation, which is age- and sex-specific and reports in",
        "mL/min/1.73 m^2. The units are confirmed by Results, which gives the",
        "cohort range as '25 to 225 mL/min/1.73 m 2 (median 100)'. Enters",
        "clearance as the through-origin LINEAR ratio (CRCL / 120) -- Table 1",
        "prints '(GFR/120)' with NO exponent, in contrast to the two weight",
        "terms in the same block which both print their exponent explicitly,",
        "so the renal effect is strictly proportional (exponent 1) and no",
        "exponent parameter is introduced here. Body size is carried",
        "separately by the allometric WT term and, because eGFR is already",
        "BSA-normalized, the two terms are not double-counting size.",
        "IMPORTANT: the 120 reference is NOT the cohort median. Results gives",
        "the median eGFR as 100 mL/min/1.73 m^2, so 120 is an external",
        "normal-renal-function anchor and the typical patient in this cohort",
        "has CL about 17% BELOW the printed 0.91 L/h. Must be strictly",
        "positive. Because the term is linear through the origin, a patient",
        "with eGFR near zero is assigned essentially zero clearance, which",
        "extrapolates well outside the fitted 25-225 range; patients on",
        "extracorporeal support were excluded at enrollment, so the model",
        "carries no information about renal replacement therapy or ECMO.",
        sep = " "
      ),
      source_name        = "GFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    age_range      = "2 months to 17 years (Results). No median or mean age is reported.",
    weight_range   = "NOT REPORTED. The abstract gives no cohort weight distribution; the only weight in the paper is the 27.4 kg allometric normalisation constant in Table 1, which the authors do not identify as a median, mean or standard.",
    disease_state  = "Critically ill children under 18 years with multiple organ dysfunction syndrome (MODS), defined as two or more organ failures, prescribed ceftriaxone as standard of care in the pediatric intensive care unit. Patients receiving extracorporeal support (e.g. ECMO, renal replacement therapy) were EXCLUDED at enrollment, so the model carries no information about those settings.",
    renal_function = "Estimated GFR by the CKiD U25 equation ranged from 25 to 225 mL/min/1.73 m^2 with a median of 100 (Results). The range spans moderate renal impairment through frank hyperfiltration / augmented renal clearance, and the Conclusion identifies residual variability after accounting for renal function as the study's central finding.",
    dose_range     = "NOT REPORTED. Methods states only that patients were 'prescribed CRO' as standard of care and that first-24-hour exposures were computed 'from exact dosing and covariate histories'; no dose, frequency or infusion duration is given anywhere in the abstract.",
    regions        = "United States (multi-center; the PALISI network, with author affiliations at Children's Hospital of Philadelphia, Nationwide Children's Hospital and Midwestern University)",
    protein_binding = "Not fitted. This model outputs WHOLE-BLOOD ceftriaxone as Cc. The fT>MIC analysis converted those concentrations to free concentrations using a literature-based free fraction of 10% (i.e. 90% protein binding), taken from the literature and not estimated here (Methods; Figure 2 caption, '90% protein binding was assumed for all subjects'). Multiply Cc by 0.10 to obtain the free concentration the paper's fT>MIC targets are evaluated against. NOTE the approximation the authors made: a protein-binding fraction is a PLASMA property, and applying it directly to a whole-blood concentration assumes whole-blood and plasma ceftriaxone are equivalent. The same group published a whole-blood-to-plasma translation for VAMS antibiotic assays (Ther Drug Monit 2026, PMC13366314), but this abstract does not cite it and does not state that any translation was applied before modelling.",
    sampling       = "Up to 15 PK samples per patient collected over 3 days by volumetric absorptive microsampling (VAMS), 20 uL per sample, quantified in whole blood by a validated LC-MS/MS assay. The total number of concentrations is not reported; Figure 1 shows roughly 250-300 points. Concentrations in Figure 1 span about 0 to 400 mg/L.",
    notes          = "CONFERENCE ABSTRACT (IDWeek 2025 poster P-1254), not a peer-reviewed full paper. It is nonetheless completely encodable: Table 1 publishes all four structural estimates, all four interindividual variances (as BOTH log-scale SD and %CV, which pins the variance convention with no ambiguity), the residual-error coefficient, and the four covariate equations in full. Fitted in Monolix 2024R1 by SAEM ('Stoch. Approx.' heads the uncertainty columns of Table 1). One- and two-compartment models were tested and the two-compartment model retained; covariates were selected on the corrected Bayesian Information Criterion, reduction in between-subject variability, and physiological relevance. As of the extraction date no peer-reviewed full publication of this ceftriaxone model exists: a EuropePMC search on the funding grant (R01HD103755) returns this abstract plus the cefepime arm of the same PALISI/VAMS study (Antimicrob Agents Chemother 2026, PMC13321836), a different drug. Treat the covariate model as provisional in the way any conference abstract's is -- the cefepime companion of this same study retained a different renal covariate form. No supplementary material exists (EuropePMC hasSuppl = N) and none is cited."
  )

  ini({
    # ---------------------------------------------------------------------
    # All final estimates are from Marianski 2026 Table 1, published as an
    # embedded raster image on page 1 of the PDF (the first of three JPEGs,
    # 975x1010 px). Its caption reads: 'The population PK parameter fixed
    # effects estimates and random effects (i.e. between subject variability)
    # for central volume of distribution(V1), clearance (CL), beta
    # coefficients on weight adjusted (individual weight/27.4kg) V1 and CL,
    # peripheral volume of distribution (V2), and intercompartmental
    # bi-directional flow (Q).'
    #
    # Table 1 has four columns: Parameter | Value | S.E. | R.S.E.(%), the last
    # two headed 'Stoch. Approx.' (the SAEM stochastic approximation). Only
    # the point estimate is carried here.
    #
    # PANEL ATTRIBUTION. The PDF page holds two unrelated abstracts in two
    # columns, so panel ownership was proved rather than assumed, on three
    # independent matches: (1) Table 1's random-effects %CV of 43.11 and 45.45
    # reproduce the Results sentence 'Between subject variation (CV%) for V1
    # and CL was 43.1% and 45.5%, respectively (Table 1)' exactly; (2)
    # Figure 2's caption names 'CRO MIC' and the 90% protein binding that this
    # abstract's Methods introduce; (3) Table 1's covariate block uses the
    # 27.4 kg weight normalisation named in its own caption. The neighbouring
    # abstract's text (hyperfiltration, BMI, sex and 'nature of the
    # hematologic disease' not correlated) belongs to a DIFFERENT study and
    # none of it is attributed to this model.
    # ---------------------------------------------------------------------

    lcl <- log(0.91)
    label("Clearance at WT = 27.4 kg and CRCL = 120 mL/min/1.73 m^2 (L/h)")
    # Table 1 Fixed Effects, row 'CL (L/hr)': 0.91, S.E. 0.083, R.S.E. 9.18%.

    lvc <- log(8.22)
    label("Central volume of distribution at WT = 27.4 kg (L)")
    # Table 1 Fixed Effects, row 'V1 (L)': 8.22, S.E. 0.74, R.S.E. 9.01%. This
    # is the paper's V1.

    lq <- log(0.75)
    label("Intercompartmental clearance (L/h)")
    # Table 1 Fixed Effects, row 'Q (L/hr)': 0.75, S.E. 0.2, R.S.E. 26.3%.
    # Carries NO covariate: Table 1 prints 'Q_i = Q_pop'.

    lvp <- log(13.74)
    label("Peripheral volume of distribution (L)")
    # Table 1 Fixed Effects, row 'V2 (L)': 13.74, S.E. 4.56, R.S.E. 33.2%.
    # This is the paper's V2, and it is the least precisely estimated of the
    # four structural parameters. Carries NO covariate: Table 1 prints
    # 'V2_i = V2_pop'. Note that the peripheral volume exceeds the central
    # volume by about 1.7-fold.

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT / 27.4) for CL (unitless)")
    # Table 1 'Model parameterized as:' block prints
    # 'CL_i = CL_pop * (WT/27.4)^0.75 * (GFR/120)'. Fixed, not estimated:
    # Methods says 'Fixed allometric scaling of clearance (CL) and central
    # volume of distribution (V1) defined the base models', and Table 1's
    # Fixed Effects block has no exponent row and reports no uncertainty for
    # either exponent.

    e_wt_vc <- fixed(1)
    label("Allometric exponent on (WT / 27.4) for V1 (unitless)")
    # Table 1 'Model parameterized as:' block prints
    # 'V1_i = V1_pop * (WT/27.4)^1'. Fixed for the same two reasons as
    # e_wt_cl. The exponent of exactly 1 is the conventional allometric value
    # for a volume.

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 1's Random Effects block reports each
    # random effect TWICE -- once as 'SD' and once as 'C.V.(%)' -- and that
    # redundancy pins the variance convention with no ambiguity, which is
    # unusual and worth stating explicitly.
    #
    # The printed pairs are (SD, CV%): V1 (0.41, 43.11), Q (0.93, 117.84),
    # V2 (1.29, 205.6), CL (0.43, 45.45). The SD is the log-scale omega, and
    # the CV% is the log-normal CV derived from it as
    # sqrt(exp(omega^2) - 1) * 100. Back-transforming the CV% column via
    # omega^2 = log(CV^2 + 1) gives omega = 0.4129 / 0.9331 / 1.2860 / 0.4333,
    # each of which rounds to the printed 2-dp SD exactly, and re-deriving the
    # CV from those reproduces 43.11 / 117.84 / 205.60 / 45.45 to every digit
    # printed.
    #
    # The competing reading -- that 'C.V.(%)' is omega * 100, as some tools
    # print -- is FALSIFIED by the same two columns: it would require
    # omega = 0.4311 for V1 against a printed SD of 0.41, and omega = 2.056
    # for V2 against a printed SD of 1.29.
    #
    # The variances below are back-transformed from the CV% column rather than
    # squared from the SD column, because the CV% is printed to four
    # significant figures where the SD is printed to two decimals. The R.S.E.
    # column independently supports the higher-precision values: for CL,
    # S.E. 0.062 / omega 0.4333 = 14.31% against the printed 14.3%, whereas
    # 0.062 / 0.43 = 14.42% would print as 14.4%.
    #
    # The matrix is diagonal -- Table 1 reports no correlations or
    # off-diagonal terms. Every one of the four disposition parameters carries
    # a random effect, including Q and V2, whose 118% and 206% CVs are
    # enormous; both are reported with acceptable precision (R.S.E. 24.3% and
    # 19.3%) and are carried as published.
    # ---------------------------------------------------------------------

    etalvc ~ 0.170457
    # omega^2 for V1. Table 1 Random Effects: SD 0.41, 'C.V.(%)' 43.11,
    # S.E. 0.064, R.S.E. 15.6%.

    etalq ~ 0.870719
    # omega^2 for Q. Table 1 Random Effects: SD 0.93, 'C.V.(%)' 117.84,
    # S.E. 0.23, R.S.E. 24.3%. A CV above 100% on an intercompartmental
    # clearance means the log-normal distribution has a long right tail; the
    # vignette checks that the sampled cohort stays numerically well behaved.

    etalvp ~ 1.653864
    # omega^2 for V2. Table 1 Random Effects: SD 1.29, 'C.V.(%)' 205.6,
    # S.E. 0.25, R.S.E. 19.3%. This is the largest random effect in the model
    # by a wide margin and it sits on the least precisely estimated structural
    # parameter; see the vignette Assumptions and deviations.

    etalcl ~ 0.187782
    # omega^2 for CL. Table 1 Random Effects: SD 0.43, 'C.V.(%)' 45.45,
    # S.E. 0.062, R.S.E. 14.3%. This is the 45.5% figure quoted in Results.

    propSd <- 0.2
    label("Proportional residual error (fraction)")
    # Table 1 'Error Model Parameters', row 'b': 0.2, S.E. 0.0091,
    # R.S.E. 4.54%. In Monolix (used here, 2024R1) 'b' is the PROPORTIONAL
    # coefficient of the error model; 'a' is the additive one. Table 1 reports
    # a 'b' row and NO 'a' row, which identifies a pure proportional error
    # model, y = f * (1 + b * e). No additive term is included here.
  })

  model({
    # Covariate reference values, both printed in the Table 1 'Model
    # parameterized as:' block. 27.4 kg is the weight normalisation named in
    # the Table 1 caption; 120 mL/min/1.73 m^2 is an external normal-renal-
    # function anchor and is NOT the cohort median of 100 (see
    # covariateData$CRCL notes).
    wt_ref   <- 27.4 # kg
    crcl_ref <- 120  # mL/min/1.73 m^2

    # Allometric size factors. Named wt_<param> rather than allo_<param>:
    # allo_cl / allo_vc are RETIRED parameter names that an enumerating test
    # in test-checkModelConventions.R forbids anywhere in inst/modeldb.
    wt_cl <- (WT / wt_ref)^e_wt_cl
    wt_vc <- (WT / wt_ref)^e_wt_vc

    # CL_i = CL_pop * (WT/27.4)^0.75 * (GFR/120) * exp(eta_CL_i)
    # The renal term is a through-origin linear ratio: Table 1 prints
    # '(GFR/120)' with no exponent.
    cl <- exp(lcl + etalcl) * wt_cl * (CRCL / crcl_ref)

    # V1_i = V1_pop * (WT/27.4)^1 * exp(eta_V1_i)
    vc <- exp(lvc + etalvc) * wt_vc

    # Q_i = Q_pop * exp(eta_Q_i)      -- no covariate, per Table 1
    q <- exp(lq + etalq)

    # V2_i = V2_pop * exp(eta_V2_i)   -- no covariate, per Table 1
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # WHOLE-BLOOD ceftriaxone, matching the assayed quantity (VAMS whole blood
    # by LC-MS/MS). The free fraction of 0.10 used for the paper's fT>MIC
    # targets is a literature constant applied outside the PK model and is
    # recorded in population$protein_binding rather than in ini().
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
