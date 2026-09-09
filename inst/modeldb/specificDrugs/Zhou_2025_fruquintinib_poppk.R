Zhou_2025_fruquintinib_poppk <- function() {
  description <- paste(
    "Joint parent plus metabolite population PK model for fruquintinib, a",
    "selective oral inhibitor of VEGFR-1, -2 and -3, and its major",
    "circulating metabolite M11, pooled from 557 subjects across five phase",
    "I/Ib studies and the global phase III FRESCO-2 trial in previously",
    "treated metastatic colorectal cancer. Fruquintinib disposition is",
    "one-compartment with first-order absorption, an absorption lag time and",
    "linear elimination; M11 is one-compartment with linear elimination, its",
    "formation flux being the fraction of the fruquintinib elimination flux",
    "fixed at 7.25% from a human mass-balance study. Apparent clearance and",
    "apparent volume of distribution of both analytes carry estimated (not",
    "fixed) body-weight power exponents centered at the 73 kg cohort median;",
    "the absorption rate constant is 60.7% lower with concurrent proton-pump",
    "inhibitor use, and fruquintinib apparent volume is 9.08% lower in",
    "healthy volunteers than in patients with cancer. Between-subject",
    "variability is carried on both clearances, both volumes and the",
    "absorption rate constant, with a full correlation block over the four",
    "disposition parameters and an uncorrelated absorption-rate eta.",
    sep = " "
  )
  reference <- paste(
    "Zhou X., Yang X., Grinshpun B., Taylor A., Strong L., Dasari A.,",
    "Wang-Gillam A., Li J., Xu R.-H., Gupta N., Chien C. (2025).",
    "Population pharmacokinetics of fruquintinib, a selective oral inhibitor",
    "of VEGFR-1, -2, and -3, in patients with refractory metastatic",
    "colorectal cancer.",
    "The Journal of Clinical Pharmacology 65(7):873-884.",
    "doi:10.1002/jcph.70001.",
    sep = " "
  )
  vignette <- "Zhou_2025_fruquintinib_poppk"

  # Zhou 2025 reports plasma concentrations in ng/mL: the bioanalytical
  # dynamic range is "1-750 ng/mL" (Methods, Analytical Methods), the additive
  # residual errors in Table 2 are 4.50 and 0.551 ng/mL, and the Discussion
  # quotes a mean steady-state trough of 228 ng/mL against a 176 ng/mL VEGFR-2
  # target. Clearances are L/h and volumes are L, so an amount in mg over a
  # volume in L gives mg/L; the model therefore carries an explicit factor of
  # 1000 on both observation lines (1 mg/L = 1000 ng/mL).
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every disposition parameter in the paper is APPARENT --
  # the paper writes CL/F and V/F for fruquintinib and CLM/f and VM/f for M11
  # -- so the state amounts are apparent amounts. The M11 state is additionally
  # in fruquintinib-mass-equivalents: the model multiplies the parent
  # elimination flux (mg of fruquintinib per hour) by the unitless fraction
  # metabolized 0.0725, so any fruquintinib-to-M11 molecular-weight ratio is
  # absorbed into the apparent VM/f, exactly as in the source.
  compartmentData <- list(
    depot = list(
      analyte = "fruquintinib", units = "mg (apparent, i.e. amount/F)",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "fruquintinib", units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma", verified = TRUE
    ),
    central_m11 = list(
      analyte = "M11", units = "mg fruquintinib-equivalents (apparent, i.e. amount/f)",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, the only continuous covariate retained in the final model; it acts on all four disposition parameters (fruquintinib CL/F and V/F, M11 CLM/f and VM/f)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered as a power function centered at 73 kg on each of the four",
        "disposition parameters (Zhou 2025 Results, 'Final Model' equation",
        "block: 'CL/F (L/h) = 0.808 x (BW/73)^0.43', 'V/F (L) = 50.4 x",
        "(BW/73)^0.924', 'CLM/f (L/h) = 0.161 x (BW/73)^1.06', 'VM/f (L) =",
        "11.7 x (BW/73)^1.28'). Methods, 'Covariate Selection': 'Continuous",
        "covariates were incorporated into the population model using a",
        "scaled structure based on the median value of the covariate in the",
        "population'. The centering value 73 kg is the overall cohort median",
        "body weight in Table 1 (73.0 kg, range 36.0-158). All four exponents",
        "are ESTIMATED, not fixed at allometric 0.75 / 1: the Discussion",
        "states 'The final model estimated the body weight-based allometry",
        "exponents for the clearances and volumes of distribution' and that",
        "both a fixed-allometry and a matched-allometry alternative performed",
        "worse. Table 2 gives %RSE 12.1, 3.9, 9.0 and 11.0 respectively, none",
        "of which is compatible with a fixed value. The centering and the",
        "exponents are jointly confirmed by the paper's own 70 kg typical",
        "values (Results, 'Final Model'): 0.808*(70/73)^0.430 = 0.7936 vs the",
        "published 0.794; 50.4*(70/73)^0.924 = 48.48 vs 48.5;",
        "0.161*(70/73)^1.06 = 0.1540 vs 0.154; 11.7*(70/73)^1.28 = 11.09 vs",
        "11.1. Baseline weight; the paper gives no indication that weight was",
        "carried as time-varying."
      ),
      source_name        = "BW"
    ),
    CONMED_PPI = list(
      description        = "Concurrent proton-pump inhibitor use, a covariate on the fruquintinib absorption rate constant",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concurrent proton-pump inhibitor)",
      notes              = paste(
        "1 = fruquintinib coadministered with a proton-pump inhibitor, 0 =",
        "not. Applied as a fractional change on the typical Ka (Zhou 2025",
        "Results, 'Final Model': 'Ka (1/h) = 2.52 (x [1 - 0.607] if",
        "coadministered with PPI)'), which is the categorical form described",
        "in Methods, 'Covariate Selection': 'Categorical covariates were",
        "incorporated into the population model using a proportional",
        "structure with either the most common level of the covariate being",
        "the reference'. Table 2 reports the coefficient as 'Ka ~ PPI",
        "fractional change = -0.607' (%RSE 11.0). The arithmetic is confirmed",
        "in the Discussion: 'when fruquintinib was given with a PPI, there was",
        "an estimated fractional change of -0.607, resulting in a Ka of 0.99",
        "per h' -- 2.52 * (1 - 0.607) = 0.990. Fruquintinib is a BCS class 2",
        "weak base with pH-dependent solubility (Introduction), which is the",
        "mechanistic rationale. The effect is on absorption RATE only, not",
        "extent: Figure 3 shows PPI coadministration produces a negligible",
        "change in steady-state exposure. Rate of PPI use in the cohort is",
        "not tabulated in Table 1; concurrent acid-reducing agents were",
        "tested 'to provide supportive data to complement the findings from",
        "dedicated drug-drug interaction studies' (Methods)."
      ),
      source_name        = "PPI"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer indicator, the paper's 'health status' covariate on fruquintinib apparent volume of distribution",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with cancer; the pooled non-healthy cohort of 515 subjects with colorectal or other solid tumors)",
      notes              = paste(
        "1 = healthy volunteer, 0 = patient with cancer. Applied as a",
        "fractional change on the typical V/F (Zhou 2025 Results, 'Final",
        "Model': 'V/F (L) = 50.4 x (BW/73)^0.924 (x [1 - 0.0908] if",
        "healthy)'); Table 2 reports 'V/F ~ Healthy subject fractional change",
        "= -0.0908' (%RSE 26.7, 95% CI -0.138 to -0.0434). Restated in the",
        "Results text as '9.08% lower V/F in healthy subjects compared with",
        "patients with cancer'. Table 1 gives 42 healthy subjects (7.5%,",
        "the two phase I studies NCT04645940 and NCT04557397) against 515",
        "patients (92.5%), so the patient group is both the reference and the",
        "most common level, consistent with Methods, 'Covariate Selection'.",
        "The paper's own conclusion is that the effect is not clinically",
        "meaningful (Figure 3: <= 2% change in fruquintinib CmaxSS and CminSS",
        "and no change in AUCSS relative to the reference)."
      ),
      source_name        = "HS"
    )
  )

  # Covariates that Zhou 2025 screened but did NOT retain in the final model.
  # These are documentation only: they are not referenced in model() and
  # checkModelConventions() does not require them to be. Recording them
  # preserves the provenance of the paper's covariate screen -- in particular
  # the fact that the paper explicitly evaluated and rejected renal function,
  # hepatic function, sex, age, race, ethnicity, country and ECOG status, which
  # is the paper's central clinical claim (no dose adjustment required).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Table 1 median 61.0 years, range 18.0-82.0. Screened; not retained.",
        "Conclusion: 'age (18.0-82.0 years) ... had no clinically meaningful",
        "impact on fruquintinib PK' and Figure 4a plots individual AUCSS",
        "against age with no trend."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "244 female / 313 male (43.8% female, Table 1). Sex was NOT retained",
        "in the final model. Two distinct traces exist and both were",
        "excluded. (1) An intermediate model from the stepwise covariate",
        "search carried 'sex on fruquintinib bioavailability' (Results,",
        "'Covariate Search'), but that term does not appear in the final",
        "model equations or in Table 2 -- it was dropped in the 'further",
        "model refinement' step that produced the final model. (2) A separate",
        "SENSITIVITY analysis estimated that 'Female subjects were estimated",
        "to have 17% lower values of both CL/F and CLM/f than male subjects,",
        "which is not considered to be clinically meaningful' (Discussion).",
        "The Discussion explains the exclusion mechanistically: 'The current",
        "model did not identify an effect of sex on fruquintinib CL/F due to",
        "the selection of body weight over sex, with which it was highly",
        "correlated.' Neither the bioavailability coefficient nor the",
        "sensitivity-analysis coefficients are printed, so no value could be",
        "encoded even if one wanted the sensitivity variant."
      )
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Retained on M11 apparent clearance in the INTERMEDIATE model coming",
        "out of the stepwise covariate search ('albumin on CLM/f', Results,",
        "'Covariate Search') but dropped in the subsequent 'further model",
        "refinement' step: no albumin term appears in the final model",
        "equations or in Table 2. No coefficient or centering value is",
        "printed anywhere in the paper, and albumin is not in Table 1, so the",
        "intermediate model cannot be reconstructed."
      )
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by the Cockcroft-Gault equation, the paper's renal-function covariate",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Table 1 median 97.9 mL/min, range 32.6-293.0. Tested because it was",
        "'deemed to be of clinical interest' (Methods) rather than on",
        "univariate significance; not retained. Figure 4f and the Results",
        "text report a 'relatively flat relationship between fruquintinib and",
        "M11 AUCSS and CrCl'."
      )
    ),
    RENALIMP_MILD = list(
      description = "Mild renal impairment indicator (CrCl category)",
      units       = "(binary)",
      type        = "binary",
      notes       = "177 subjects (31.8%, Table 1). Screened; not retained."
    ),
    RENALIMP_MOD = list(
      description = "Moderate renal impairment indicator (CrCl category)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "42 subjects (7.5%, Table 1). Screened; not retained. Conclusion:",
        "'mild to moderate renal impairment ... had no clinically meaningful",
        "impact on fruquintinib PK'."
      )
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator by NCI Organ Dysfunction Working Group criteria",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "133 subjects (23.9%, Table 1). Screened as clinically-of-interest;",
        "not retained. Discussion: 'mild hepatic impairment (total bilirubin",
        "<= upper limit of normal [ULN] and aspartate aminotransferase [AST]",
        "> ULN, or total bilirubin > 1.0 to 1.5 x ULN, with any AST) had no",
        "clinically meaningful impact on fruquintinib exposure'."
      )
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator by NCI Organ Dysfunction Working Group criteria",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Only 2 subjects (0.4%, Table 1). Results: 'There were only two",
        "subjects with moderate hepatic impairment, preventing an adequate",
        "assessment of this category.'"
      )
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Table 1 median 9.30 umol/L, range 2.57-38.0. Part of the hepatic",
        "function screen; not retained."
      )
    ),
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Table 1 median 25.2 kg/m^2, range 16.0-56.7. Part of the body-size",
        "screen; body weight was the retained body-size covariate."
      )
    ),
    ECOG_GE1 = list(
      description = "Eastern Cooperative Oncology Group performance-status >= 1 indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "310 subjects with ECOG 1 vs 247 with ECOG 0 (Table 1); no subject",
        "had ECOG >= 2, so the >= 1 indicator fully encodes the paper's",
        "'disease severity' covariate. Screened; not retained."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "140 subjects (25.1%, Table 1). Screened; not retained (Figure 4c)."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "29 subjects (5.2%, Table 1). Screened; not retained (Figure 4c)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "359 subjects (64.5%, Table 1). Screened; not retained (Figure 4c)."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "28 subjects (5.0%, Table 1). Screened; not retained (Figure 4d)."
    ),
    REGION_CHINA = list(
      description = "China enrollment-country indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "80 subjects (14.4%, Table 1; the two China-only studies",
        "NCT01645215 and NCT01975077). Screened; not retained. Discussion:",
        "'race (Asian, Black, and White), ethnicity (Hispanic and",
        "non-Hispanic) and country were not identified as significant",
        "covariates on fruquintinib CL/F.'"
      )
    ),
    REGION_JAPAN = list(
      description = "Japan enrollment-country indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "45 subjects (8.1%, Table 1; the Japanese arm of FRESCO-2). Screened; not retained."
    ),
    TUMTP_CRC = list(
      description = "Colorectal cancer tumor-type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "470 subjects with colorectal cancer, 45 with another tumor type and",
        "42 with no tumor (Table 1). Screened; not retained. Abstract: 'tumor",
        "type ... had no clinically meaningful impact on fruquintinib or M11",
        "PK'."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 557L,
    n_studies      = 6L,
    n_observations = "6668 post-treatment fruquintinib concentrations from 557 subjects; 4136 post-treatment M11 concentrations from 460 subjects (Results, 'Summary of Analysis Dataset')",
    age_range      = "18.0-82.0 years",
    age_median     = "61.0 years",
    weight_range   = "36.0-158 kg",
    weight_median  = "73.0 kg",
    bmi_median     = "25.2 kg/m^2 (range 16.0-56.7)",
    sex_female_pct = 43.8,
    race_ethnicity = c(
      White = 64.5, Black = 5.2, Asian = 25.1,
      HawaiianPacificIslander = 0.2, Multiple = 0.5, Other = 1.1, Missing = 3.4,
      Hispanic = 5.0
    ),
    disease_state  = "previously treated metastatic colorectal cancer and other advanced solid tumors (515 patients, 92.5%) plus healthy volunteers (42 subjects, 7.5%)",
    renal_function = "creatinine clearance 32.6-293.0 mL/min (median 97.9); 337 normal (60.5%), 177 mild (31.8%), 42 moderate (7.5%) impairment",
    hepatic_function = "NCI Organ Dysfunction Working Group category: 421 normal (75.6%), 133 mild (23.9%), 2 moderate (0.4%)",
    performance_status = "ECOG 0 in 247 subjects (44.3%), ECOG 1 in 310 (55.7%); no subject had ECOG >= 2",
    dose_range     = "oral fruquintinib 1-6 mg once daily across the phase I/Ib dose-ranging studies; 5 mg once daily for 21 days of each 28-day cycle in FRESCO-2",
    regions        = "China (14.4%), Japan (8.1%), rest of world (77.6%)",
    notes          = paste(
      "Baseline demographics are Table 1 of Zhou 2025, reported overall and",
      "by study. The six pooled studies are NCT01645215 (n = 40, China),",
      "NCT01975077 (n = 40, China), US1/NCT03251378 (n = 101),",
      "FRESCO-2/NCT04322539 (n = 334, global phase III), NCT04645940 (n = 14,",
      "healthy volunteers) and NCT04557397 (n = 28, healthy volunteers);",
      "Table S1 (not required for the model) summarizes them. Subjects in the",
      "FRESCO-2 placebo arm were excluded from the analysis dataset, as were",
      "169 (2.5%) fruquintinib and 1063 (20.4%) M11 concentrations below the",
      "assay lower limit of quantification -- the M3 method for M11 BLQ data",
      "was tested but failed to achieve a stable minimization and covariance",
      "step, so M11 BLQ observations were treated as missing."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # All values are the FINAL population estimates of Zhou 2025 Table 2
    # ("Parameter Estimates for the Final Population PK Model"), cross-checked
    # against the equation block printed in Results, "Final Model". Table 2
    # footnote: "Typical values correspond to a patient with cancer weighing
    # 73 kg and not taking any PPIs" -- so every typical value below is the
    # 73 kg, no-PPI, cancer-patient reference.
    #
    # The reference weight and all four body-weight exponents are confirmed
    # arithmetically by the paper's OWN derived 70 kg typical values quoted in
    # the same Results paragraph (see the WT covariateData note); every one of
    # the four reproduces to three significant figures.
    # -----------------------------------------------------------------------

    # -- Fruquintinib absorption. First order with a lag time. Ka is the
    #    no-PPI typical value; the lag time carries no covariate and no BSV
    #    ("BSV was applied to all structural model parameters except Tlag",
    #    Results, "Final Model").
    lka   <- log(2.52)  ; label("Fruquintinib first-order absorption rate constant without concurrent PPI (1/h, log scale)") # Zhou 2025 Table 2 'Ka (1/h) = 2.52' (%RSE 7.6, 95% CI 2.15-2.90); absorption half-life log(2)/2.52 = 0.275 h, matching the published 0.275 h
    ltlag <- log(0.463) ; label("Fruquintinib absorption lag time (h, log scale)")                                            # Zhou 2025 Table 2 'Lag time (h) = 0.463' (%RSE 1.1, 95% CI 0.453-0.473); restated in Results 'Final Model' as 'Tlag (h) = 0.463'

    # -- Fruquintinib disposition, one compartment. Both are APPARENT values
    #    (divided by the oral bioavailability F, which the paper never
    #    separates out) at the 73 kg reference weight.
    lcl <- log(0.808) ; label("Fruquintinib apparent clearance CL/F at 73 kg in patients with cancer (L/h, log scale)")               # Zhou 2025 Table 2 'CL/F (L/h) = 0.808' (%RSE 1.2, 95% CI 0.789-0.827); elimination half-life log(2)*50.4/0.808 = 43.2 h, matching the published 43.2 h
    lvc <- log(50.4)  ; label("Fruquintinib apparent volume of distribution V/F at 73 kg in patients with cancer (L, log scale)")     # Zhou 2025 Table 2 'V/F (L) = 50.4' (%RSE 0.9, 95% CI 49.5-51.3)

    # -- M11 disposition, one compartment. Both are apparent values divided by
    #    the paper's lowercase "f"; because the formation flux below is
    #    explicitly scaled by the fixed fraction metabolized 0.0725, the
    #    residual apparent scaling absorbed into these two parameters is the
    #    parent bioavailability together with any fruquintinib-to-M11
    #    molecular-weight ratio.
    lcl_m11 <- log(0.161) ; label("M11 apparent clearance CLM/f at 73 kg (L/h, log scale)")               # Zhou 2025 Table 2 'CLM/f (L/h) = 0.161' (%RSE 2.2, 95% CI 0.154-0.169)
    lvc_m11 <- log(11.7)  ; label("M11 apparent volume of distribution VM/f at 73 kg (L, log scale)")     # Zhou 2025 Table 2 'VM/f (L) = 11.7' (%RSE 3.3, 95% CI 10.9-12.4)

    # -- Fraction of the fruquintinib elimination flux forming M11. FIXED, not
    #    estimated: it is an external constant carried in from a dedicated
    #    human mass-balance study and therefore appears in neither Table 2 nor
    #    the final-model equation block.
    fm <- fixed(0.0725) ; label("Fraction of fruquintinib metabolized to M11, from the NCT02689752 mass-balance study (unitless)") # Zhou 2025 Results, 'Base Model Development': "Based on results from a phase I human mass balance study (2015-013-00CH2; NCT02689752), it was estimated that approximately 7.25% of the fruquintinib dose was metabolized to M11, and this value was used in the model as the fraction of fruquintinib metabolized to M11"; restated in the Discussion

    # -- Body-weight power exponents, all ESTIMATED (Discussion: "The final
    #    model estimated the body weight-based allometry exponents for the
    #    clearances and volumes of distribution"; fixed-allometry and
    #    matched-allometry alternatives were tested and performed worse).
    e_wt_cl     <- 0.430 ; label("Power exponent on (WT/73) for fruquintinib CL/F (unitless)")   # Zhou 2025 Table 2 'CL/F ~ Body weight exponent = 0.430' (%RSE 12.1, 95% CI 0.328-0.532)
    e_wt_vc     <- 0.924 ; label("Power exponent on (WT/73) for fruquintinib V/F (unitless)")    # Zhou 2025 Table 2 'V/F ~ Body weight exponent = 0.924' (%RSE 3.9, 95% CI 0.853-0.995)
    e_wt_cl_m11 <- 1.06  ; label("Power exponent on (WT/73) for M11 CLM/f (unitless)")           # Zhou 2025 Table 2 'CLM/f ~ Body weight exponent = 1.06' (%RSE 9.0, 95% CI 0.873-1.25)
    e_wt_vc_m11 <- 1.28  ; label("Power exponent on (WT/73) for M11 VM/f (unitless)")            # Zhou 2025 Table 2 'VM/f ~ Body weight exponent = 1.28' (%RSE 11.0, 95% CI 1.00-1.56)

    # -- Categorical covariate effects, both entered as FRACTIONAL CHANGES on
    #    the reference typical value, i.e. the parameter is multiplied by
    #    (1 + coefficient * indicator). This is the "proportional structure"
    #    of Methods, "Covariate Selection", and is exactly how the two effects
    #    are printed in the Results equation block.
    e_conmed_ppi_ka  <- -0.607  ; label("Fractional change in fruquintinib Ka with concurrent PPI use (unitless)")        # Zhou 2025 Table 2 'Ka ~ PPI fractional change = -0.607' (%RSE 11.0); 2.52 * (1 - 0.607) = 0.990 1/h, matching the Discussion's "resulting in a Ka of 0.99 per h"
    e_dis_healthy_vc <- -0.0908 ; label("Fractional change in fruquintinib V/F in healthy volunteers vs patients with cancer (unitless)") # Zhou 2025 Table 2 'V/F ~ Healthy subject fractional change = -0.0908' (%RSE 26.7, 95% CI -0.138 to -0.0434); Results restates it as "9.08% lower V/F in healthy subjects"

    # -----------------------------------------------------------------------
    # Between-subject variability. Methods, "Structural and Statistical Model
    # Development": IIV "was modeled assuming a log-normal distribution",
    # theta_ki = theta_k * exp(eta_ki) with eta ~ N(0, omega_k^2), which is
    # the exponential form used in model() below.
    #
    # SCALE. Table 2's BSV rows are headed "%CV" and the table's own footnote
    # gives the formula: "BSV %CV is calculated as sqrt(exp(Omega -1)), where
    # Omega is the BSV variance". The stray parenthesis is a typesetting slip
    # for the standard log-normal identity CV = sqrt(exp(omega^2) - 1); read
    # literally as written, sqrt(exp(Omega) - 1) with Omega already the
    # variance, which is the same thing. Either way the printed footnote
    # settles the scale outright -- these are NOT plain 100*omega values --
    # so every variance below is
    #
    #     omega^2 = log(1 + (%CV / 100)^2)
    #
    # and NOT (%CV/100)^2. The distinction is immaterial for V/F (0.02622 vs
    # 0.02657) but decisive for Ka, where 247% gives omega^2 = 1.960 rather
    # than 6.101 -- a factor of 3.1 in the variance.
    #
    # STRUCTURE. Results, "Final Model": "BSV was applied to all structural
    # model parameters except Tlag, with correlations between all BSV
    # components except for BSV on Ka." So the four disposition etas form one
    # full 4x4 block and the absorption eta is diagonal. The block below is
    # in the order CL/F, V/F, CLM/f, VM/f and is written lower-triangle
    # row-major; each off-diagonal is corr * omega_i * omega_j using the
    # correlation coefficients printed in Table 2. The resulting correlation
    # matrix is positive definite (smallest eigenvalue 0.177).
    #
    # Shrinkage, reported in Table 2 for every BSV term, is 8.3% (CL/F),
    # 19.4% (V/F), 10.6% (CLM/f), 12.5% (VM/f) and 21.0% (Ka), all below the
    # 30% the paper set as its acceptance criterion.
    # SOURCE TRACE FOR THE BLOCK BELOW. Every number is stated here rather
    # than as a trailing comment inside the `c(...)`: rxode2 rewrites trailing
    # comments in ini() into label() calls, and doing that inside a
    # multi-line correlated-eta block produces unparseable code.
    #
    # Diagonal (variance = log(1 + (%CV/100)^2)), Zhou 2025 Table 2:
    #   var CL/F  = log(1 + 0.262^2) = 0.06639056  'BSV CL/F %CV = 26.2'
    #                                              (%RSE 3.5, 95% CI 24.3-28.0)
    #   var V/F   = log(1 + 0.163^2) = 0.02622217  'BSV V/F %CV = 16.3'
    #                                              (%RSE 4.9, 95% CI 14.6-17.8)
    #   var CLM/f = log(1 + 0.493^2) = 0.21756720  'BSV CLM/f %CV = 49.3'
    #                                              (%RSE 3.5, 95% CI 45.5-53.0)
    #   var VM/f  = log(1 + 0.753^2) = 0.44916870  'BSV VM/f %CV = 75.3'
    #                                              (%RSE 3.7, 95% CI 68.4-82.0)
    #
    # Off-diagonals = correlation * omega_i * omega_j, with the correlations
    # printed in Table 2 and omega_i = sqrt(variance) above:
    #   cov(CL/F, V/F)    = 0.295 * 0.257664 * 0.161933 = 0.01230862  (%RSE 10.9)
    #   cov(CL/F, CLM/f)  = 0.563 * 0.257664 * 0.466441 = 0.06766410  (%RSE  4.9)
    #   cov(V/F,  CLM/f)  = 0.167 * 0.161933 * 0.466441 = 0.01261385  (%RSE 19.7)
    #   cov(CL/F, VM/f)   = 0.413 * 0.257664 * 0.670200 = 0.07131944  (%RSE  6.8)
    #   cov(V/F,  VM/f)   = 0.495 * 0.161933 * 0.670200 = 0.05372103  (%RSE  6.4)
    #   cov(CLM/f, VM/f)  = 0.708 * 0.466441 * 0.670200 = 0.22132720  (%RSE  4.4)
    #
    # Written lower-triangle row-major in the order CL/F, V/F, CLM/f, VM/f.
    # -----------------------------------------------------------------------
    etalcl + etalvc + etalcl_m11 + etalvc_m11 ~ c(
      0.06639056,
      0.01230862, 0.02622217,
      0.06766410, 0.01261385, 0.21756720,
      0.07131944, 0.05372103, 0.22132720, 0.44916870
    )
    # var Ka = log(1 + 2.47^2) ; Zhou 2025 Table 2 'BSV Ka %CV = 247'
    # (%RSE 4.4, 95% CI 201-299; shrinkage 21.0%). Uncorrelated with the four
    # disposition etas per Results, 'Final Model'.
    etalka ~ 1.96022200

    # -----------------------------------------------------------------------
    # Residual unexplained variability. Methods: a combined proportional plus
    # additive model, Y_ij = C_ij * (1 + eps_1ij) + eps_2ij. Results, "Final
    # Model": "A proportional and additive error model was used for both
    # fruquintinib and M11", with "separate components for each analyte".
    # The proportional terms are printed as %CV and the additive terms as an
    # SD in ng/mL, so propSd is the percentage over 100 and addSd is taken
    # directly. Fruquintinib residual shrinkage is 8.5% (Table 2).
    # -----------------------------------------------------------------------
    propSd     <- 0.199 ; label("Fruquintinib proportional residual error (fraction)")   # Zhou 2025 Table 2 'Fruquintinib proportional error %CV = 19.9' (%RSE 1.2, 95% CI 19.4-20.4; shrinkage 8.5%)
    addSd      <- 4.50  ; label("Fruquintinib additive residual error (ng/mL)")          # Zhou 2025 Table 2 'Fruquintinib additive error SD (ng/mL) = 4.50' (%RSE 4.1, 95% CI 4.14-4.86)
    propSd_m11 <- 0.195 ; label("M11 proportional residual error (fraction)")            # Zhou 2025 Table 2 'M11 proportional error %CV = 19.5' (%RSE 1.4, 95% CI 18.9-20.0)
    addSd_m11  <- 0.551 ; label("M11 additive residual error (ng/mL)")                   # Zhou 2025 Table 2 'M11 additive error SD (ng/mL) = 0.551' (%RSE 4.0, 95% CI 0.509-0.594)
  })

  model({
    # 1. Unit bridge. Doses are in mg and volumes in L, so an amount over a
    #    volume is mg/L; the paper's concentrations and additive residual
    #    errors are in ng/mL, and 1 mg/L = 1000 ng/mL.
    ugPerMg <- 1000

    # 2. Individual fruquintinib absorption parameters. The PPI effect is a
    #    fractional change, so CONMED_PPI = 0 recovers the typical 2.52 1/h
    #    exactly and CONMED_PPI = 1 gives 0.990 1/h. The lag time has neither
    #    a covariate nor an eta.
    ka   <- exp(lka + etalka) * (1 + e_conmed_ppi_ka * CONMED_PPI)
    tlag <- exp(ltlag)

    # 3. Individual fruquintinib disposition parameters. Body weight enters as
    #    a power function centered at the 73 kg cohort median; health status
    #    enters V/F as a fractional change, so DIS_HEALTHY = 0 (a patient with
    #    cancer, the reference) recovers the typical value exactly.
    cl <- exp(lcl + etalcl) * (WT / 73)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 73)^e_wt_vc *
      (1 + e_dis_healthy_vc * DIS_HEALTHY)

    # 4. Individual M11 disposition parameters. Body weight only.
    cl_m11 <- exp(lcl_m11 + etalcl_m11) * (WT / 73)^e_wt_cl_m11
    vc_m11 <- exp(lvc_m11 + etalvc_m11) * (WT / 73)^e_wt_vc_m11

    # 5. ODE system. Fruquintinib is one-compartment with first-order
    #    absorption and a lag; M11 is one-compartment, formed as the fraction
    #    `fm` of the fruquintinib elimination flux and cleared linearly. Note
    #    that `cl` is the WHOLE fruquintinib elimination -- the M11 route is a
    #    sub-fraction of it, not a parallel arm added on top -- so `cl / vc`
    #    is genuinely the parent elimination rate constant and the parent
    #    profile is unaffected by `fm`.
    #
    #    The steady-state consequence of writing the formation flux with an
    #    explicit `fm` is a metabolite-to-parent AUC ratio of
    #    fm * (CL/F) / (CLM/f) = 0.0725 * 0.808 / 0.161 = 0.364, against the
    #    "mean metabolite-to-parent area under the plasma concentration-time
    #    curve (AUC) ratio ... of 0.3 at steady state" reported in the
    #    Introduction from the separate mass-balance study. Dropping the `fm`
    #    factor would give 5.02 instead, so the placement is falsifiable and
    #    confirmed. See the vignette's steady-state mass-balance check.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * central / vc
    d/dt(central_m11) <-  fm * cl * central / vc - cl_m11 * central_m11 / vc_m11

    # 6. Absorption lag on the depot.
    alag(depot) <- tlag

    # 7. Observations, both in ng/mL, each with its own combined residual
    #    error.
    Cc     <- ugPerMg * central     / vc
    Cc_m11 <- ugPerMg * central_m11 / vc_m11

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_m11 ~ add(addSd_m11) + prop(propSd_m11)
  })
}
