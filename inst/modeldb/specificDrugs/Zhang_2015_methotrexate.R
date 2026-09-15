Zhang_2015_methotrexate <- function() {
  description <- "Two-compartment population PK model for intravenous high-dose methotrexate (8-12 g/m^2 infused over 4-6 h) in Chinese osteosarcoma patients at a single Beijing institution; central clearance carries a linear effect of the number of prior methotrexate chemotherapy cycles and of pre-dose creatinine clearance, and intercompartmental clearance and peripheral volume each carry a linear effect of body surface area (Zhang 2015)"
  reference <- paste(
    "Zhang W, Zhang Q, Tian X, Zhao H, Lu W, Zhen J, Niu X. (2015).",
    "Population Pharmacokinetics of High-dose Methotrexate After Intravenous",
    "Administration in Chinese Osteosarcoma Patients from a Single Institution.",
    "Chin Med J (Engl) 128(1):111-118.",
    "doi:10.4103/0366-6999.147829.",
    sep = " "
  )
  vignette <- "Zhang_2015_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Plasma methotrexate was quantified by fluorescence
  # polarization immunoassay (TDX, Abbott) with a quantification limit of
  # 0.01 umol/L, and every concentration in Zhang 2015 is expressed in umol/L
  # (Methods, "Methotrexate assay" and "Blood collection"). The amount unit
  # that pairs with the published V1 (L) and CL1 (L/h) is therefore umol.
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CYCLE = list(
      description        = "Index of the current high-dose methotrexate chemotherapy course (1 = first course)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Source column MTXNUM, glossed by Zhang 2015 as 'the number of methotrexate chemotherapy cycles before",
        "MTX infusion' (Abstract) and 'the time of chemotherapy using MTX before chemotherapy' (Results).",
        "Table 1 gives median 2, range 1-12, n = 270. The minimum of 1 -- not 0 -- in a dataset that must contain",
        "first courses identifies the column as a 1-based counter of the CURRENT course rather than a count of",
        "strictly prior courses, which is exactly the canonical CYCLE semantics.",
        "This is the ONE covariate in the final model that is NOT mean-centered: it enters as the bare count",
        "(1 - 0.0183 * CYCLE), so the tabulated typical CL1 of 6.20 L/h is the value extrapolated to CYCLE = 0,",
        "which never occurs in the data. At the median CYCLE of 2 the multiplier is 0.963.",
        "Zhang 2015 also screened NUM, the number of chemotherapy cycles of ANY drug (Table 1, median 8,",
        "range 1-28). NUM and MTXNUM were strongly correlated (Figure 2) and the paper retained MTXNUM only,",
        "explicitly to avoid an over-complicated model; see covariatesDataExcluded.",
        "Mechanistic reading offered by the paper (Discussion): repeated courses of high-dose methotrexate",
        "progressively impair renal function, so clearance falls as the course count rises.",
        sep = " "
      ),
      source_name        = "MTXNUM"
    ),
    CRCL = list(
      description        = "Pre-dose creatinine clearance (raw, NOT BSA-normalized), on the source's own numeric scale",
      units              = "mL/s",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CrCl1, 'creatinine clearance before administration'. Table 1 gives median 1.88,",
        "minimum 0.94, maximum 4.64, n = 245, and the final-model equation centers on 1.89.",
        "UNITS HAZARD -- READ BEFORE SUPPLYING THIS COLUMN. Zhang 2015 Table 1 labels this row 'ml/min', which",
        "cannot be right: a creatinine clearance of 1.88 mL/min is anuric, whereas the cohort is a hyperhydrated,",
        "renally intact group of adolescents and young adults. The tabulated magnitude is instead consistent with",
        "the SI unit mL/s (1.88 mL/s = 113 mL/min; the 0.94-4.64 range = 56-278 mL/min), which is a routine",
        "reporting convention in Chinese clinical laboratories. The same table independently mislabels serum",
        "creatinine as 'mg/dl' with median 49 and range 19-81 -- values that are only sane as umol/L -- so this",
        "table's unit column is demonstrably unreliable and the SI reading is corroborated twice over.",
        "A weight-normalized reading (mL/min/kg; 1.88 x 58 kg = 109 mL/min) reproduces essentially the same",
        "absolute clearance for the median patient and cannot be excluded from the paper alone; mL/s is recorded",
        "here as the better-supported of the two because the companion creatinine row is likewise SI.",
        "WHAT MATTERS FOR USING THIS MODEL: the coefficient 0.0416 is calibrated per unit ON THE SOURCE'S OWN",
        "NUMERIC SCALE and is centered at 1.89 on that same scale. Drive this model with creatinine clearance",
        "expressed on that scale (median 1.88, range 0.94-4.64) -- i.e. divide a conventional mL/min value by 60 --",
        "whatever the correct unit label turns out to be. Supplying a raw mL/min value would multiply the",
        "covariate deviation by about 60 and is the single largest misuse risk this model carries.",
        "Renal function was very likely MEASURED from a 24 h urine collection rather than estimated: the protocol",
        "recorded 24 h urine volume and Table 1 carries a separate urine-volume row. A Cockcroft-Gault estimate",
        "for the median patient (age 17, 58 kg, creatinine 49 umol/L) would give about 178 mL/min, roughly 1.6x",
        "the tabulated median, which is the expected direction of Cockcroft-Gault bias when creatinine production",
        "is low.",
        "SIGN: entered here as (1 + 0.0416 * (CRCL - 1.89)), i.e. clearance RISES with renal function. See the",
        "extended note in ini() -- the sign printed in the paper's own final-model equation is the opposite, and",
        "is contradicted by the paper's Methods template, its Discussion, and the renal elimination of",
        "methotrexate.",
        sep = " "
      ),
      source_name        = "CrCl1"
    ),
    BSA = list(
      description        = "Body surface area at the time of the course",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BODYAREA. Table 1 gives median 1.63 m^2, range 0.62-2.21, n = 270; the Subjects section",
        "reports 1.63 +/- 0.27 m^2. The final-model equations center on 1.62 m^2, consistent with the Methods",
        "covariate template, which centers each covariate on its MEAN rather than its median.",
        "The BSA computation formula (Du Bois / Mosteller / Haycock) is not stated by Zhang 2015; record it as",
        "unspecified. The cohort spans children through adults (weight 20-97 kg), which is why the BSA range",
        "reaches down to 0.62 m^2.",
        "BSA carries effects on BOTH the intercompartmental clearance CL2 and the peripheral volume V2, and on",
        "neither central parameter. It is also the dosing covariate: high-dose methotrexate is prescribed in",
        "g/m^2, so BSA sets the delivered dose as well as two disposition parameters.",
        "Zhang 2015 notes that weight, BMI and BSA were strongly mutually correlated (Figure 2) and that BSA was",
        "retained in preference to the other two; see covariatesDataExcluded.",
        "SIGN: entered here as (1 + 0.880 * (BSA - 1.62)) on CL2 and (1 + 0.874 * (BSA - 1.62)) on V2, i.e. both",
        "RISE with body size. See the extended note in ini() -- the sign printed in the paper's own final-model",
        "equations is the opposite of this and of the paper's own Discussion.",
        sep = " "
      ),
      source_name        = "BODYAREA"
    )
  )

  # Covariates that Zhang 2015 tabulated and screened but did NOT retain in the
  # final model. Documented here so the provenance of the paper's covariate
  # screen survives without triggering "declared but not referenced" warnings.
  # Table 1 additionally lists urinary pH, body temperature, urine volume, the
  # volume of infusion the day before dosing, serum sodium and serum chloride,
  # none of which have canonical register entries; those are recorded in
  # population$notes instead.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Table 1 median 17 years, range 6-49, n = 274; Subjects section 17.00 +/- 7.06 years. Screened, not retained.",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Table 1 median 58 kg, range 20-97, n = 268. Strongly correlated with BMI and BSA (Figure 2); BSA was retained in preference. Screened, not retained.",
      source_name = "Weight"
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Table 1 median 19.95 kg/m^2, range 11.32-40.17, n = 267. Strongly correlated with weight and BSA (Figure 2); BSA was retained in preference. Screened, not retained.",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "194 of 274 courses were in males and 80 in females, i.e. 29% female (Subjects section). Screened as 'Gender' (Table 1, Figure 2), not retained.",
      source_name = "Gender"
    ),
    CREAT = list(
      description = "Pre-dose serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Source column Cr1. Table 1 median 49, range 19-81, n = 250, labelled 'mg/dl' -- a label that cannot be",
        "correct at that magnitude and is only sane as umol/L (49 umol/L = 0.55 mg/dL). This mislabelling is the",
        "corroborating evidence that the companion CrCl1 row is likewise in SI units; see covariateData$CRCL.",
        "Screened, not retained -- the derived creatinine CLEARANCE was retained instead.",
        sep = " "
      ),
      source_name = "Cr1"
    ),
    ALB = list(
      description = "Pre-dose serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Table 1 median 41.2 g/L, range 24.4-51.7, n = 250. Screened, not retained.",
      source_name = "Albumin1"
    ),
    ALP = list(
      description = "Pre-dose alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Table 1 median 94 U/L, range 43-479, n = 238. Often raised in osteosarcoma. Screened, not retained.",
      source_name = "AKP1"
    ),
    HCT = list(
      description = "Pre-dose hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Table 1 median 36.5%, range 20.4-48.9, n = 262. Screened, not retained.",
      source_name = "HCT1"
    ),
    RBC = list(
      description = "Pre-dose erythrocyte count",
      units       = "10^12/L",
      type        = "continuous",
      notes       = "Table 1 median 4.06, range 2.28-5.70, n = 264, printed as '1 x 109/L' but at that magnitude the values are 10^12/L (a conventional red-cell count). Screened, not retained.",
      source_name = "RBC1"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 148L,
    n_studies         = 1L,
    n_administrations = 274L,
    age_range         = "6-49 years (Table 1); 17.00 +/- 7.06 years (Subjects section)",
    age_median        = "17 years",
    weight_range      = "20-97 kg (Table 1); 58.00 +/- 18.28 kg (Subjects section)",
    weight_median     = "58 kg",
    bsa_median        = "1.63 m^2 (range 0.62-2.21); 1.63 +/- 0.27 m^2 (Subjects section); BSA formula unspecified",
    height_mean       = "166.00 +/- 12.44 cm",
    sex_female_pct    = 29,
    race_ethnicity    = "Chinese; race was not evaluated as a covariate",
    disease_state     = "Osteosarcoma receiving neoadjuvant or adjuvant high-dose methotrexate with leucovorin rescue",
    dose_range        = "8-12 g/m^2 methotrexate in 500 mL 5% glucose infused intravenously over 4-6 h in darkness (Methods, 'High-dose methotrexate administration'); the Abstract quotes the regimen as 10 g/m^2",
    renal_function    = "Pre-dose creatinine clearance median 1.88 on the source's own numeric scale (range 0.94-4.64; see covariateData$CRCL for the units hazard), pre-dose serum creatinine median 49 umol/L (range 19-81). No renal-impairment stratification is reported and no exclusion for renal dysfunction is stated.",
    co_medication     = "Protocolized hyperhydration and urine alkalinization before, during and after the infusion (5% glucose, 5% glucose-saline, potassium chloride and 5% sodium bicarbonate), plus 2 mg vincristine and 5 mg tropisetron. Oral sodium bicarbonate 1.0 g three times daily and allopurinol 200 mg three times daily throughout. Leucovorin rescue 12 mg every 6 h beginning 6-8 h after the end of the infusion and continued until plasma methotrexate fell below 0.05 umol/L, escalated per protocol if elimination was delayed.",
    regions           = "Single center: Beijing Jishuitan Hospital, Beijing, China",
    notes             = paste(
      "148 patients contributed 274 high-dose methotrexate courses between August 2009 and August 2010.",
      "IMPORTANT STRUCTURAL CAVEAT: repeat courses in the same patient were treated as INDEPENDENT individuals",
      "during model building. Zhang 2015 says so explicitly (Discussion): 'There were some patients who received",
      "HD-MTX multiple times, and we treated them as totally separate individuals in the modeling process to",
      "gather more data, but that also meant that we ignored the internal correlation within patients.'",
      "Consequently the reported inter-individual variances are BETWEEN-COURSE variances that confound true",
      "inter-individual variability with between-occasion variability, and the n of the random-effects model is",
      "274 courses rather than 148 subjects. The paper's own abstract frames the analysis as exploring",
      "'between-occasion variability', but no separate IOV level was estimated.",
      "Plasma methotrexate was sampled at 0, 6, 12, 24, 48 and 72 h after the start of the infusion, with",
      "additional 24-hourly samples until the concentration fell below 0.05 umol/L, and assayed by fluorescence",
      "polarization immunoassay (TDX, Abbott) with a 0.01 umol/L quantification limit, 90-110% recovery and",
      "interday precision below 10%.",
      "Covariate screening used a sample-size-independent effect-size statistic with a 0.5 cut-off to drop",
      "laboratory indices that changed materially over the course, then stepwise forward inclusion at dOFV > 3.84",
      "and backward elimination at dOFV > 10.83. Table 1 additionally tabulates urinary pH before dosing",
      "(median 8, range 5-9), body temperature (median 36.3 C, range 35.3-37.9), urine volume the day before",
      "dosing (median 1820 mL, range 420-5670), infusion volume the day before dosing (median 1520 mL, range",
      "710-4240), serum sodium (median 142.5 mmol/L) and serum chloride (median 103 mmol/L); none were retained.",
      "The paper attributes the null result for urinary pH to poor measurement (patients self-read pH paper to",
      "one digit) and the null result for hydration to the uniformly high infusion volumes used at this center.",
      "NONMEM Version V Level 1.1, first-order conditional estimation. The objective function fell from 487.856",
      "in the base model to 373.294 in the final model (dOFV = 114.562). Internal validation was a 1000-replicate",
      "bootstrap. The paper reports its own goodness-of-fit diagnostics candidly: the model fits poorly at C0,",
      "and 'a portion of the WRESs was beyond the -4-4 range, with a few >20' even in the final model.",
      "No external validation was performed.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters: Zhang 2015 Table 2, "Final model standard value"
    # column. Two-compartment intravenous disposition, parameterized by the
    # paper as CL1 / V1 / CL2 / V2.
    #
    # CL2 is the INTERCOMPARTMENTAL clearance, not a second elimination
    # pathway. Table 2 and the Results text call it "clearance of the
    # peripheral compartment", which is ambiguous, but the Discussion resolves
    # it: "intercompartmental clearance (0.0172 L/h) and peripheral
    # distribution volume (0.515 L) values were in the same range as those
    # previously published by Aquerreta et al. (CLD1 0.053 L/h, V2 1.82 L)" --
    # CLD1 being an intercompartmental clearance. It is therefore mapped to
    # the canonical lq.
    # --------------------------------------------------------------------
    lcl <- log(6.20)   ; label("Central clearance, CL1 (L/h)")                    # Zhang 2015 Table 2 final model: CL1 6.20 L/h, RSE 4.87% (base model 5.81, RSE 2.07%)
    lvc <- log(19.6)   ; label("Central volume of distribution, V1 (L)")          # Zhang 2015 Table 2 final model: V1 19.6 L, RSE 4.39% (base model 19.2, RSE 2.34%)
    lq  <- log(0.0172) ; label("Intercompartmental clearance, CL2 (L/h)")         # Zhang 2015 Table 2 final model: CL2 0.0172 L/h, RSE 14.9% (base model 0.0154, RSE 5.61%)
    lvp <- log(0.515)  ; label("Peripheral volume of distribution, V2 (L)")       # Zhang 2015 Table 2 final model: V2 0.515 L, RSE 9.92% (base model 0.471, RSE 4.88%)

    # --------------------------------------------------------------------
    # COVARIATE-EFFECT SIGNS -- the single most important note in this file.
    #
    # Zhang 2015 prints its final model twice (Abstract, and Results "Final
    # regression model equation and parameter values"). Both printings put a
    # MINUS in all four covariate brackets:
    #
    #   CL1_i = CL1_tv * [1 - t_MTXNUM * MTXNUM]
    #                  * [1 - t_CrCl   * (CrCl1 - 1.89)]     * exp(eta_CL1)
    #   CL2_i = CL2_tv * [1 - t_BSA    * (BODYAREA - 1.62)]  * exp(eta_CL2)
    #   V2_i  = V2_tv  * [1 - t_BSA    * (BODYAREA - 1.62)]  * exp(eta_V2)
    #
    # Three of those four minus signs are wrong. This model uses a PLUS for
    # the CrCl term and for both BSA terms, and keeps the MINUS for MTXNUM.
    # The evidence, in descending order of force:
    #
    # 1. The paper's OWN generic covariate template (Methods, "Fixed effect
    #    model") is   P_TVij = P_TVj * [1 + theta_jk * (COVR_ik - COVR_k)]
    #    -- a PLUS, with the covariate centered on its mean. pdftotext renders
    #    that plus and the adjacent (COVR_ik - COVR_k) minus as two DIFFERENT
    #    glyphs in the same line, so the minus signs in the final-model
    #    equations are genuine characters, not a decoding artifact, and the
    #    template genuinely disagrees with them.
    #    This is the decisive structural tell: the three disputed terms are
    #    written in exactly the template's mean-centered form, (CrCl1 - 1.89)
    #    and (BODYAREA - 1.62), so they are instances of the template and
    #    inherit its plus. The ONE term that departs from the template -- the
    #    uncentered bare MTXNUM -- is the one whose minus is corroborated.
    #
    # 2. The Discussion states all four directions in words, and matches the
    #    mixed reading exactly: "The clearance rate of MTX decreased with
    #    increased times of MTX chemotherapy or a decreased creatinine
    #    clearance rate, while body surface area had a positive correlation
    #    with the peripheral clearance rate and the apparent volume of
    #    distribution of the peripheral compartment."
    #
    # 3. Physiology. Methotrexate is predominantly excreted unchanged by the
    #    kidney, so CL1 must RISE with creatinine clearance; the printed minus
    #    makes it fall. And under the printed minus the peripheral volume
    #    SHRINKS with body size -- V2 would be 0.97 L for a 0.62 m^2 child
    #    against 0.25 L for a 2.21 m^2 adult, a four-fold inversion of the
    #    usual size-volume relationship.
    #
    # The most economical explanation is that a single minus-bracket template
    # was copied across all four terms when the equations were typeset.
    # Applied under the standing "covariate equation absurd as printed ->
    # sensible centered interpretation" policy; flagged in the vignette's
    # Errata as the first thing a reviewer should check against the source.
    #
    # The centering constants 1.89 and 1.62 are taken verbatim from the
    # printed equations and are deliberately NOT the Table 1 medians of 1.88
    # and 1.63 -- the Methods template centers on the MEAN, and the equation
    # constant is the value the tabulated typical CL1, CL2 and V2 are
    # conditioned on.
    # --------------------------------------------------------------------
    e_cycle_cl <- 0.0183 ; label("Linear coefficient of methotrexate course number on CL1 (per course)")   # Zhang 2015 Table 2 final model: theta CL1-MTXNUM 0.0183, RSE 35.6%; enters as (1 - e_cycle_cl * CYCLE), UNCENTERED
    e_crcl_cl  <- 0.0416 ; label("Linear coefficient of pre-dose creatinine clearance on CL1 (per mL/s)")  # Zhang 2015 Table 2 final model: theta CL1-CrCl 0.0416, RSE 32.2%; enters as (1 + e_crcl_cl * (CRCL - 1.89)); sign flipped from the printed equation, see the block comment above
    e_bsa_q    <- 0.880  ; label("Linear coefficient of body surface area on CL2 (per m^2)")               # Zhang 2015 Table 2 final model: theta CL2-BODYAREA 0.880, RSE 28.3%; enters as (1 + e_bsa_q * (BSA - 1.62)); sign flipped from the printed equation, see the block comment above
    e_bsa_vp   <- 0.874  ; label("Linear coefficient of body surface area on V2 (per m^2)")                # Zhang 2015 Table 2 final model: theta V2-BODYAREA 0.874, RSE 21.2%; enters as (1 + e_bsa_vp * (BSA - 1.62)); sign flipped from the printed equation, see the block comment above

    # --------------------------------------------------------------------
    # Inter-individual variability. Zhang 2015 Table 2 reports an
    # "Inter-individual RSD %" column, and the Statistical model section
    # defines the random effect as P_ij = P_TVj * exp(eta_ij) with eta
    # normally distributed with variance omega^2. The RSD percentages are
    # therefore read as 100 * omega, i.e. the SD of eta on the log scale, and
    # the variances below are their squares:
    #   CL1  8.48% -> 0.0848^2 = 0.00719
    #   CL2 50.9%  -> 0.509^2  = 0.259
    #   V2  39.1%  -> 0.391^2  = 0.153
    # For an exponential random effect this reading and the exact lognormal
    # one, omega^2 = log(1 + CV^2), coincide to three digits at 8.48% and
    # differ by about 11% of the variance at 50.9% (0.230 vs 0.259); the
    # direct reading is used because the paper names the quantity an RSD of
    # eta rather than a CV of the parameter.
    #
    # NOTE ON WHAT THESE VARIANCES MEAN: because repeat courses in the same
    # patient were modelled as separate individuals (see population$notes),
    # these are BETWEEN-COURSE variances and confound inter-individual with
    # between-occasion variability.
    #
    # V1 carries NO random effect here. The printed equation shows
    # V1_i = V1_tv * exp(eta_V1), but Table 2 leaves the inter-individual RSD
    # cell blank for V1 in BOTH the base and the final model, and the Abstract
    # describes the V1 variability as "extremely small". No value is reported
    # anywhere in the paper, so the eta is omitted rather than invented.
    # --------------------------------------------------------------------
    etalcl ~ 0.00719  # Zhang 2015 Table 2 final model, inter-individual RSD 8.48% on CL1; variance = 0.0848^2 (base model 8.93%)
    etalq  ~ 0.259    # Zhang 2015 Table 2 final model, inter-individual RSD 50.9% on CL2; variance = 0.509^2 (base model 55.0%)
    etalvp ~ 0.153    # Zhang 2015 Table 2 final model, inter-individual RSD 39.1% on V2; variance = 0.391^2 (base model 47.3%)

    # --------------------------------------------------------------------
    # Residual unexplained variability. Zhang 2015 specifies the FORM
    # exactly (Methods, "Residual random effect model"):
    #   C_obs = C_pred * (1 + eps1) + eps2
    # i.e. a combined proportional-plus-additive error, with eps1 the
    # proportional and eps2 the additive component. The paper reports the
    # form but never reports the MAGNITUDE of either component: Table 2
    # carries no sigma rows, the text quotes no sigma values, and there is no
    # supplement. Both are therefore encoded at fixed(0) rather than invented,
    # following the standing "unreported RUV -> fixed(0) plus erratum" policy
    # and the precedent in Setiawan_2023_sulbactam.R.
    #
    # Practical consequence: simulations from this model carry structural and
    # inter-individual variability but NO residual error, so simulated
    # concentrations are individual predictions. This is flagged in the
    # vignette Errata.
    # --------------------------------------------------------------------
    propSd <- fixed(0) ; label("Proportional residual error (fraction; 0 -- form given but magnitude not reported in the source)")  # Zhang 2015 Methods 'Residual random effect model' gives C_obs = C_pred * (1 + eps1) + eps2 but no sigma values anywhere
    addSd  <- fixed(0) ; label("Additive residual error (umol/L; 0 -- form given but magnitude not reported in the source)")        # Zhang 2015 Methods 'Residual random effect model' gives C_obs = C_pred * (1 + eps1) + eps2 but no sigma values anywhere
  })

  model({
    # Individual parameters. Zhang 2015 final-model equations, with the three
    # sign corrections documented at length in ini(). CYCLE enters UNCENTERED
    # (bare count); CRCL and BSA enter mean-centered on 1.89 and 1.62.
    #
    # CRCL must be supplied on the source's own numeric scale (median 1.88) --
    # see covariateData$CRCL. Supplying conventional mL/min here would inflate
    # the covariate deviation roughly 60-fold.
    cl <- exp(lcl + etalcl) * (1 - e_cycle_cl * CYCLE) * (1 + e_crcl_cl * (CRCL - 1.89))
    vc <- exp(lvc)
    q  <- exp(lq  + etalq)  * (1 + e_bsa_q  * (BSA - 1.62))
    vp <- exp(lvp + etalvp) * (1 + e_bsa_vp * (BSA - 1.62))

    # Two-compartment intravenous disposition with first-order elimination
    # from the central compartment. The dose is delivered into `central` over
    # the 4-6 h infusion via the event table's rate / dur column; there is no
    # absorption compartment.
    d/dt(central)     <- q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <- q / vc * central     - q / vp * peripheral1

    # Plasma methotrexate concentration. The state holds umol and vc is in L,
    # so Cc is in umol/L, matching every concentration reported in Zhang 2015.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
