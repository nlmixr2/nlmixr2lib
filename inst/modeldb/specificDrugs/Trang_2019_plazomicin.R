Trang_2019_plazomicin <- function() {
  description <- "Three-compartment population PK model for intravenous plazomicin in healthy adults and adults with complicated urinary tract infection, acute pyelonephritis, bloodstream infection or hospital-/ventilator-acquired bacterial pneumonia, with a sigmoidal Hill relationship between renal clearance and creatinine clearance and an additive continuous-renal-replacement-therapy clearance arm"
  reference <- paste(
    "Trang M, Seroogy JD, Van Wart SA, Bhavnani SM, Kim A, Gibbons JA,",
    "Ambrose PG, Rubino CM. 2019. Population pharmacokinetic analyses for",
    "plazomicin using pooled data from phase 1, 2, and 3 clinical studies.",
    "Antimicrob Agents Chemother 63:e02329-18. doi:10.1128/AAC.02329-18.",
    sep = " "
  )
  vignette <- "Trang_2019_plazomicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters CL as (WT / 75)^0.529 and Vp2 as (WT / 75)^1.62. The paper ",
        "does not print the normalising weight; 75 kg is used because the ",
        "Results sentence that reports the typical CL of 4.57 liters/h ",
        "describes 'a typical cUTI or HABP/VABP patient with a BW of 75 kg, ",
        "a BSA of 1.73 m2, and a CLCR of 90 ml/min', and that CL is ",
        "reproduced exactly by the Hill term alone (0.491 + 4.80 * ",
        "90^2.49 / (45.3^2.49 + 90^2.49) = 4.56 liters/h), i.e. with the ",
        "weight factor equal to 1. 75 kg is also the Table 1 pooled median ",
        "body weight. Range 40.5-165 kg (Table 1)."
      ),
      source_name = "BW"
    ),
    BSA = list(
      description = "Body surface area by the Du Bois and Du Bois equation",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters Vc as (BSA / 1.73)^1.23 and Vp1 as (BSA / 1.73)^1.17. The ",
        "1.73 m2 reference is the value named in the Results typical-patient ",
        "sentence. Materials and Methods: BSA = BW(kg)^0.425 * ",
        "height(cm)^0.725 * 0.007184. Pooled median 1.86 m2, range ",
        "1.29-2.58 m2 (Table 1)."
      ),
      source_name = "BSA"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters CLd2 (the second distributional clearance, q2 here) as ",
        "(HT / 170)^3.38. The paper does not print the normalising height; ",
        "170 cm is the Table 1 pooled median. Range 142-194 cm."
      ),
      source_name = "height"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters Vp1 as the linear proportional term (1 + 0.00954 * ",
        "(AGE - 39)). Table 2 labels the coefficient a 'Vp1-age slope', ",
        "which is why a linear rather than a power form is used; the ",
        "centring age is not printed and 39 years is the Table 1 pooled ",
        "median. Range 18-90 years. The term stays positive across the ",
        "whole observed age range (0.80 at 18 years, 1.49 at 90 years)."
      ),
      source_name = "age"
    ),
    CRCL = list(
      description = paste0(
        "Cockcroft-Gault creatinine clearance normalised to a body surface ",
        "area of 1.73 m^2; time-varying"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Materials and Methods 'Subject characteristics': CLCR (mL/min) for ",
        "males = (140 - age) * BW / (72 * SCr), for females = male CLCR * ",
        "0.85, with serum creatinine floored at 0.50 mg/dL, then normalised ",
        "to 1.73 m2 of Du Bois BSA. Recalculated and carried forward as a ",
        "time-varying covariate on every day a central-laboratory serum ",
        "creatinine was measured, with linear interpolation between measured ",
        "values. Drives the renal-clearance arm through the sigmoidal Hill ",
        "term. Pooled median 90.2, range 7.37-226 mL/min/1.73 m2 (Table 1). ",
        "The Discussion cautions that data were sparse above 150 ",
        "mL/min/1.73 m2."
      ),
      source_name = "CLCR"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection infection-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, when every DIS_* indicator is 0)",
      notes = paste0(
        "One of five mutually exclusive infection-type strata (healthy, ",
        "cUTI, AP, BSI, HABP/VABP; Table 1). cUTI and AP share a single ",
        "estimated proportional shift on Vc, CLd1, Vp1 and CLd2 (Results ",
        "'Final covariate model refinement': shifts with similar estimates ",
        "were combined and tested as one group). The CL:cUTI shift was ",
        "dropped during refinement, so cUTI patients take the reference CL. ",
        "213/564 subjects (37.8%)."
      ),
      source_name = "Infection type = cUTI"
    ),
    DIS_AP = list(
      description = "Acute pyelonephritis infection-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, when every DIS_* indicator is 0)",
      notes = paste0(
        "Grouped with DIS_CUTI for the Vc, CLd1, Vp1 and CLd2 shifts, but ",
        "carries its own +13.0% shift on CL (Table 2, 'Proportional increase ",
        "for AP patients'). 164/564 subjects (29.1%)."
      ),
      source_name = "Infection type = AP"
    ),
    DIS_BACTEREMIA = list(
      description = "Bloodstream-infection (bacteremia) infection-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, when every DIS_* indicator is 0)",
      notes = paste0(
        "Grouped with DIS_HABP and DIS_VABP for the Vc and CLd2 shifts, but ",
        "carries its own -18.9% shift on CL (Table 2, 'Proportional increase ",
        "for BSI patients'). The CLd1:BSI and Vp1:BSI relationships were ",
        "dropped after the bootstrap for poor precision, so BSI patients ",
        "take the reference CLd1 and the cUTI/AP-free Vp1. 29/564 subjects ",
        "(5.14%), all from the phase 3 CRE study 007."
      ),
      source_name = "Infection type = BSI"
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia infection-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, when every DIS_* indicator is 0)",
      notes = paste0(
        "Grouped with DIS_VABP and DIS_BACTEREMIA for the Vc and CLd2 ",
        "shifts. The CL:HABP/VABP shift was dropped during refinement, so ",
        "these patients take the reference CL, and CLd1:HABP/VABP and ",
        "Vp1:HABP/VABP were dropped after the bootstrap. The source reports ",
        "HABP and VABP as a single pooled stratum (15/564 subjects, 2.66%) ",
        "and estimates no separate HABP and VABP coefficients; the two ",
        "canonical columns are therefore always used together and always ",
        "carry the same coefficient."
      ),
      source_name = "Infection type = HABP/VABP"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia infection-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, when every DIS_* indicator is 0)",
      notes = paste0(
        "Pooled with DIS_HABP by the source into one HABP/VABP stratum; see ",
        "the DIS_HABP notes. Set both to 0 and DIS_HABP to 1, or both to 1, ",
        "to select the stratum -- the model uses their sum together with ",
        "DIS_BACTEREMIA, so exactly one of the five infection-type ",
        "indicators should be non-zero for a given subject."
      ),
      source_name = "Infection type = HABP/VABP"
    ),
    CONMED_INOTROPE = list(
      description = "Vasopressor / inotrope use at any time during the study",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no vasopressor)",
      notes = paste0(
        "Table 1 footnote d: administration of adrenaline, dobutamine, ",
        "etilefrine, isoproterenol, noradrenaline or norepinephrine at any ",
        "time during the study, handled as a time-fixed yes/no covariate. ",
        "Retained only on Vp2, where it raises the second peripheral volume ",
        "by 390% (Table 2); the Vp1:vasopressor relationship was dropped ",
        "after the bootstrap for poor precision. 24/564 subjects (4.26%), ",
        "all in the phase 3 CRE study 007."
      ),
      source_name = "Vasopressor use"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = paste0(
        "Total continuous-renal-replacement-therapy effluent flow rate, ",
        "i.e. the sum of the ultrafiltrate and dialysate flow rates"
      ),
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Materials and Methods: 'Clearance due to CRRT was set to the sum of ",
        "the actual patient-specific DFR and UFR, and multiplied by an ",
        "estimate of the sieving coefficient'. DFR and UFR were fixed from ",
        "the source data and only the sieving coefficient was estimated, so ",
        "unlike the Zurawska 2026 / Butragueno-Laiseca 2022 precedent the ",
        "sieving coefficient is NOT absorbed into an estimated clearance ",
        "parameter here -- it is the separate model parameter sc_crrt. ",
        "Table 2 gives the observed sum of UFR and DFR as 1.14-1.8 ",
        "liters/h, i.e. 1140-1800 mL/h; the model divides this column by ",
        "1000 to reach the liters/h in which clearance is expressed. Only 9 ",
        "of 564 subjects (study 007) received CRRT."
      ),
      source_name = "UFR + DFR"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Continuous-renal-replacement-therapy active indicator (time-varying gate)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CRRT not operative)",
      notes = paste0(
        "Results: CRRT clearance 'was estimated only during those periods ",
        "when CRRT was operative', and Materials and Methods assign the ",
        "timing of CRRT from the source data, so this is a within-subject ",
        "time-varying gate rather than a subject-level flag. Set to 1 only ",
        "while the circuit is running; the CRRT clearance arm is then added ",
        "to the residual clearance."
      ),
      source_name = "CRRT timing"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2 study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pooled phase 1 studies when STUDY_PHASE3 is also 0)",
      notes = paste0(
        "Selects the residual-error magnitude only; it touches no structural ",
        "or covariate parameter. Table 2 tabulates one constant-coefficient- ",
        "of-variation component per study phase. Set to 1 for the phase 2 ",
        "cUTI/AP study 002 (92 subjects)."
      ),
      source_name = "Study phase = 2"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pooled phase 1 studies when STUDY_PHASE2 is also 0)",
      notes = paste0(
        "Selects the residual-error magnitude only. Set to 1 for the phase 3 ",
        "studies 007 and 009 (329 subjects), which used sparse sampling; the ",
        "phase 1 reference stratum carries the rich single- and multiple- ",
        "dose profiles and the smallest residual error."
      ),
      source_name = "Study phase = 3"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste0(
        "Results '(ii) Covariate analysis': 'the occasions were categorized ",
        "as days 1 to 2, 3 to 6, 6 to 9, and >9'. Encode OCC = 1 for ",
        "0 <= t < 48 h, 2 for 48 <= t < 144 h, 3 for 144 <= t < 216 h and 4 ",
        "for t >= 216 h. The paper's day-6 boundary is printed as belonging ",
        "to both the second and the third occasion; the upper-exclusive ",
        "reading used here is the only one that makes the four occasions ",
        "mutually exclusive. IOV was retained only on CL and is small ",
        "(3.59% CV)."
      ),
      source_name = "occasion"
    )
  )

  compartmentData <- list(
    central = list(analyte = "plazomicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "plazomicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "plazomicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 564,
    n_studies = 7,
    age_range = "18-90 years",
    age_median = "39 years",
    weight_range = "40.5-165 kg",
    weight_median = "75.0 kg",
    height_range = "142-194 cm",
    bsa_range = "1.29-2.58 m^2",
    sex_female_pct = 52.8,
    race_ethnicity = c(
      White = 78.5,
      Black = 9.22,
      Asian = 4.26,
      `American Indian/Alaskan Native` = 7.09,
      Other = 0.887
    ),
    disease_state = paste0(
      "143 healthy adults (25.4%) plus 421 adults with complicated urinary ",
      "tract infection (213, 37.8%), acute pyelonephritis (164, 29.1%), ",
      "bloodstream infection (29, 5.14%) or hospital-acquired / ",
      "ventilator-associated bacterial pneumonia (15, 2.66%)"
    ),
    renal_function = paste0(
      "Creatinine clearance 7.37-226 mL/min/1.73 m^2 (median 90.2). Normal ",
      "renal function 237 (42.0%), mild impairment 186 (33.0%), moderate ",
      "128 (22.7%), severe 13 (2.30%). Nine subjects in study 007 received ",
      "continuous renal replacement therapy during treatment."
    ),
    co_medication = "Vasopressors in 24 subjects (4.26%); positive-pressure ventilation in 24 (4.26%); diabetes in 63 (11.2%)",
    dose_range = paste0(
      "1-20 mg/kg as 30-min intravenous infusions; single doses and q12h, ",
      "q24h or q48h multiple doses, with phase 3 doses adjusted for ",
      "baseline creatinine clearance and, when body weight was >=125% of ",
      "ideal body weight, based on adjusted body weight"
    ),
    n_observations = 4990,
    notes = paste0(
      "Baseline characteristics are Table 1 of Trang 2019, pooled across ",
      "four phase 1 studies (001, 003, 004, 006), one phase 2 study (002) ",
      "and two phase 3 studies (007, 009). Per-study retention is Table S1 ",
      "of the supplement: 573 subjects and 5,142 samples at start, 564 and ",
      "4,990 retained after excluding 106 outliers and 46 ",
      "below-limit-of-quantification samples."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Clearance. Total CL is the sum of a constant non-renal arm and a
    # saturable renal arm driven by creatinine clearance:
    #
    #   CL = (CL_nonrenal + CL_Rmax * CRCL^hill / (CRCL50^hill + CRCL^hill))
    #        * (WT / 75)^e_wt_cl * infection-type shift
    #
    # Results '(iv) Final population PK model': 'The relationship between
    # renal clearance (CLR) and CLCR was described using a sigmoidal
    # Hill-type function, and the relationship between total CL of
    # plazomicin and CLCR included an intercept to represent nonrenal
    # clearance. The other parameters ... were CLR maximum, CLCR50, and a
    # Hill coefficient.'
    #
    # This reproduces the paper's own typical values: at CRCL = 90 the
    # Hill fraction is 0.8468, so CL_R = 4.80 * 0.8468 = 4.06 liters/h
    # (paper: 4.08) and total CL = 0.491 + 4.06 = 4.56 liters/h (paper:
    # 4.57), with the weight factor equal to 1 at the 75 kg reference.
    # -----------------------------------------------------------------
    lcl_nonren <- log(0.491)
    label("Non-renal clearance (L/h)") # Table 2, 'Nonrenal CL (liters/h)' = 0.491 (%SEE 24.9)
    lcl_renal_max <- log(4.80)
    label("Maximum (asymptotic) renal clearance (L/h)") # Table 2, 'CLR maximum (liters/h)' = 4.80 (%SEE 7.39)
    lcrcl50 <- log(45.3)
    label("Creatinine clearance giving half-maximal renal clearance (mL/min/1.73 m^2)") # Table 2, 'Baseline CLCR50 (ml/min/1.73 m2)' = 45.3 (%SEE 5.27)
    lhill <- log(2.49)
    label("Hill coefficient of the renal-clearance / creatinine-clearance relationship (unitless)") # Table 2, 'Hill coefficient' = 2.49 (%SEE 13.9)
    e_wt_cl <- 0.529
    label("Body-weight power on clearance (unitless)") # Table 2, 'CL-weight power' = 0.529 (%SEE 14.0)
    e_ap_cl <- 0.130
    label("Proportional change in clearance for acute-pyelonephritis patients (fraction)") # Table 2, 'Proportional increase for AP patients' = 0.130 (%SEE 22.9)
    e_bsi_cl <- -0.189
    label("Proportional change in clearance for bloodstream-infection patients (fraction)") # Table 2, 'Proportional increase for BSI patients' = -0.189 (%SEE 41.0)

    # Central volume
    lvc <- log(9.10)
    label("Central volume at BSA = 1.73 m^2 in a healthy subject (L)") # Table 2, Vc 'Coefficient' = 9.10 (%SEE 4.07)
    e_bsa_vc <- 1.23
    label("Body-surface-area power on the central volume (unitless)") # Table 2, 'Vc-BSA power' = 1.23 (%SEE 17.5)
    e_cutiap_vc <- 1.05
    label("Proportional change in the central volume for cUTI and AP patients (fraction)") # Table 2, Vc 'Proportional increase for cUTI and AP patients' = 1.05 (%SEE 10.6)
    e_bsihabp_vc <- 1.55
    label("Proportional change in the central volume for BSI and HABP/VABP patients (fraction)") # Table 2, Vc 'Proportional increase for BSI and HABP/VABP patients' = 1.55 (%SEE 17.0)

    # First distributional clearance (the paper's CLd1)
    lq <- log(8.05)
    label("Distributional clearance to peripheral1 in a healthy subject (L/h)") # Table 2, CLd1 'Coefficient' = 8.05 (%SEE 7.97)
    e_cutiap_q <- -0.831
    label("Proportional change in the first distributional clearance for cUTI and AP patients (fraction)") # Table 2, CLd1 'Proportional increase for cUTI and AP patients' = -0.831 (%SEE 4.85)

    # First peripheral volume (the paper's Vp1)
    lvp <- log(8.71)
    label("First peripheral volume at BSA = 1.73 m^2 and age 39 y in a healthy subject (L)") # Table 2, Vp1 'Coefficient' = 8.71 (%SEE 3.97)
    e_bsa_vp <- 1.17
    label("Body-surface-area power on the first peripheral volume (unitless)") # Table 2, 'Vp1-BSA power' = 1.17 (%SEE 22.2)
    e_age_vp <- 0.00954
    label("Linear age slope on the first peripheral volume (fraction per year)") # Table 2, 'Vp1-age slope' = 0.00954 (%SEE 11.1)
    e_cutiap_vp <- -0.437
    label("Proportional change in the first peripheral volume for cUTI and AP patients (fraction)") # Table 2, Vp1 'Proportional increase for cUTI and AP patients' = -0.437 (%SEE 14.6)

    # Second distributional clearance (the paper's CLd2)
    lq2 <- log(0.199)
    label("Distributional clearance to peripheral2 at height 170 cm in a healthy subject (L/h)") # Table 2, CLd2 'Coefficient' = 0.199 (%SEE 3.64)
    e_ht_q2 <- 3.38
    label("Height power on the second distributional clearance (unitless)") # Table 2, 'CLd2-height power' = 3.38 (%SEE 17.5)
    e_cutiap_q2 <- -0.299
    label("Proportional change in the second distributional clearance for cUTI and AP patients (fraction)") # Table 2, CLd2 'Proportional increase for cUTI and AP patients' = -0.299 (%SEE 46.0)
    e_bsihabp_q2 <- 2.86
    label("Proportional change in the second distributional clearance for BSI and HABP/VABP patients (fraction)") # Table 2, CLd2 'Proportional increase for BSI and HABP/VABP patients' = 2.86 (%SEE 31.7)

    # Second peripheral volume (the paper's Vp2)
    lvp2 <- log(6.98)
    label("Second peripheral volume at WT = 75 kg without vasopressors (L)") # Table 2, Vp2 'Coefficient' = 6.98 (%SEE 9.21)
    e_wt_vp2 <- 1.62
    label("Body-weight power on the second peripheral volume (unitless)") # Table 2, 'Vp2-weight power' = 1.62 (%SEE 20.8)
    e_inotrope_vp2 <- 3.90
    label("Proportional change in the second peripheral volume for vasopressor use (fraction)") # Table 2, Vp2 'Proportional increase for vasopressor use' = 3.90 (%SEE 36.0)

    # Continuous renal replacement therapy
    sc_crrt <- 0.734
    label("Plazomicin sieving coefficient across the CRRT membrane (unitless)") # Table 2, CLCRRT 'Sieving coefficient' = 0.734 (%SEE 94.7); re-estimated during final refinement from the base-model value of 0.926

    # -----------------------------------------------------------------
    # Inter-individual variability. Table 2 reports omega-squared values
    # directly, and Materials and Methods states that IIV used 'an
    # exponential error model that assumed that these parameters are
    # log-normally distributed', so the tabulated values are variances on
    # the log scale and are used verbatim.
    #
    # Table 2 estimates exactly three covariances: CL-Vc, CL-Vp1 and
    # Vc-Vp1. The Results prose instead names the pairs eta-CL/eta-Vc,
    # eta-CL/eta-CLd1 and eta-Vc/eta-CLd1; the table is the final
    # parameter register and is followed here. Both readings give a
    # positive-definite block (the table's implied correlations are 0.632,
    # 0.878 and 0.543).
    # -----------------------------------------------------------------
    etalcl + etalvc + etalvp ~ c(
      0.103,
      0.0931, 0.211,
      0.0734, 0.0649, 0.0678
    ) # Table 2 omega-squared for CL / Vc / Vp1 and the three estimated covariances
    etalq ~ 0.0661 # Table 2, 'omega-squared for CLd1' = 0.0661 (%SEE 47.8)
    etalq2 ~ 0.0350 # Table 2, 'omega-squared for CLd2' = 0.0350 (%SEE 24.2)
    etalvp2 ~ 0.170 # Table 2, 'omega-squared for Vp2' = 0.170 (%SEE 34.9)

    # Inter-occasion variability on clearance. One variance is reported
    # and it is shared by all four occasions, so occasions 2-4 repeat the
    # occasion-1 value as fixed() -- the same encoding used for a NONMEM
    # $OMEGA BLOCK(1) SAME.
    etaiov_cl_1 ~ 0.00129 # Table 2, 'IOV on CL' = 0.00129 (%SEE 61.1); 3.59% CV
    etaiov_cl_2 ~ fixed(0.00129) # Table 2, 'IOV on CL'; shared across occasions
    etaiov_cl_3 ~ fixed(0.00129) # Table 2, 'IOV on CL'; shared across occasions
    etaiov_cl_4 ~ fixed(0.00129) # Table 2, 'IOV on CL'; shared across occasions

    # -----------------------------------------------------------------
    # Residual variability: an additive plus constant-coefficient-of-
    # variation (proportional) model with a separate proportional
    # component per study phase. Table 2 reports sigma-squared, so each
    # standard deviation is the square root of the tabulated variance.
    # -----------------------------------------------------------------
    addSd <- sqrt(0.0000414)
    label("Additive residual error (mg/L)") # Table 2, sigma-squared 'Additive component' = 0.0000414 (%SEE 51.5); sqrt = 0.006434 mg/L
    propSdPhase1 <- sqrt(0.0297)
    label("Proportional residual error, phase 1 studies (fraction)") # Table 2, 'CCV component for phase 1 studies' = 0.0297 (%SEE 8.99); sqrt = 0.1723
    propSdPhase2 <- sqrt(0.168)
    label("Proportional residual error, phase 2 study (fraction)") # Table 2, 'CCV component for phase 2 studies' = 0.168 (%SEE 14.9); sqrt = 0.4099
    propSdPhase3 <- sqrt(0.0846)
    label("Proportional residual error, phase 3 studies (fraction)") # Table 2, 'CCV component for phase 3 studies' = 0.0846 (%SEE 8.78); sqrt = 0.2909
  })

  model({
    # 1. Derived covariate terms.
    #
    # Infection-type groupings. Results '(iii) Final covariate model
    # refinement': proportional shifts with similar estimates were
    # combined, pairing cUTI with AP and BSI with HABP/VABP. Healthy
    # subjects are the reference stratum for Vc, CLd1, Vp1 and CLd2; for
    # CL the reference additionally absorbs cUTI and HABP/VABP, whose
    # shifts were dropped as negligible.
    cutiap <- DIS_CUTI + DIS_AP
    bsihabp <- DIS_BACTEREMIA + DIS_HABP + DIS_VABP

    # Occasion indicators for the inter-occasion variability on CL
    # (days 1-2, 3-6, 6-9, >9).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 +
      oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Sigmoidal Hill relationship between the renal clearance arm and
    # creatinine clearance.
    hill <- exp(lhill)
    crcl50 <- exp(lcrcl50)
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal_max) * CRCL^hill / (crcl50^hill + CRCL^hill)

    # Clearance contributed by the extracorporeal circuit while CRRT is
    # running. Materials and Methods: 'The total CL for patients while on
    # CRRT was equal to the sum of residual CL and the CRRT CL', so this
    # arm is additive on top of the individual clearance and carries no
    # IIV of its own. The effluent-flow column is in mL/h and clearance
    # in L/h, hence the division by 1000.
    cl_crrt <- sc_crrt * (RRT_CRRT_EFFLUENT_FLOW / 1000) * RRT_CRRT_ACTIVE

    # 2. Individual PK parameters.
    cl <- (cl_nonren + cl_renal) *
      (WT / 75)^e_wt_cl *
      (1 + e_ap_cl * DIS_AP + e_bsi_cl * DIS_BACTEREMIA) *
      exp(etalcl + iov_cl) +
      cl_crrt
    vc <- exp(lvc + etalvc) *
      (BSA / 1.73)^e_bsa_vc *
      (1 + e_cutiap_vc * cutiap + e_bsihabp_vc * bsihabp)
    q <- exp(lq + etalq) *
      (1 + e_cutiap_q * cutiap)
    vp <- exp(lvp + etalvp) *
      (BSA / 1.73)^e_bsa_vp *
      (1 + e_age_vp * (AGE - 39)) *
      (1 + e_cutiap_vp * cutiap)
    q2 <- exp(lq2 + etalq2) *
      (HT / 170)^e_ht_q2 *
      (1 + e_cutiap_q2 * cutiap + e_bsihabp_q2 * bsihabp)
    vp2 <- exp(lvp2 + etalvp2) *
      (WT / 75)^e_wt_vp2 *
      (1 + e_inotrope_vp2 * CONMED_INOTROPE)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. ODE system. Dosing is a zero-order 30-min intravenous infusion
    # into the central compartment, supplied through the event table's
    # infusion rate; the model itself has no absorption compartment.
    d/dt(central) <- -(kel + k12 + k13) * central +
      k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 5. Observation and residual error. The proportional component is
    # selected by study phase; phase 1 is the reference stratum.
    Cc <- central / vc
    propSdSel <- propSdPhase1 * (1 - STUDY_PHASE2 - STUDY_PHASE3) +
      propSdPhase2 * STUDY_PHASE2 +
      propSdPhase3 * STUDY_PHASE3
    Cc ~ add(addSd) + prop(propSdSel)
  })
}
