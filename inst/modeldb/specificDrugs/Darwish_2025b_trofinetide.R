Darwish_2025b_trofinetide <- function() {
  description <- "Updated population PK model for oral trofinetide in Rett syndrome (Darwish 2025b, DAFFODIL): two-compartment with first-order absorption and linear elimination, re-estimated after adding pediatric data from girls aged 2-4 years to the 13-study pool behind Darwish_2025a_trofinetide."
  reference <- "Darwish M, Passarell J, Maxwell K, Bradley H, Bishop KM, Youakim JM. Population Pharmacokinetics of Trofinetide in a Pediatric Population Aged 2 to 4 Years with Rett Syndrome. Advances in Therapy. 2025;42(2):1009-1025. doi:10.1007/s12325-024-03058-7"
  vignette <- "Darwish_2025b_trofinetide"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Trofinetide concentrations were quantified in
  # lithium-heparinized WHOLE BLOOD (Methods, "Study Design and Populations"),
  # not plasma.
  compartmentData <- list(
    depot       = list(analyte = "trofinetide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "trofinetide", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "trofinetide", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance; reference 58 kg, the median body weight of the analysis population (Darwish 2025b Table 2 footnote, 'WTKG/58'). Analysis-population range 9.8-140 kg (mean 56.5). The exponent rose from 0.443 in Darwish_2025a_trofinetide to 0.486 here, the largest covariate change apart from the two Rett syndrome shifts; Darwish 2025b Discussion attributes this to the addition of the 2-4 year old DAFFODIL cohort (mean weight 13.4 kg), which extends the weight range downward by 3.2 kg.",
      source_name        = "WTKG"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume; reference 22.4 years, the median age of the analysis population (Darwish 2025b Table 2 footnote, 'AGE/22.4'). Analysis-population range 2-64 years (mean 21.8). Note that the reference age is unchanged from Darwish_2025a_trofinetide even though the age range now extends down to 2 years.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Glomerular filtration rate, BSA-normalized (creatinine-based estimate)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance; reference 124 mL/min/1.73 m^2, the median GFR of the analysis population (Darwish 2025b Table 2 footnote, 'GFR/124'). The source column is named GFR and maps to the canonical general-scope CRCL covariate, which covers BSA-normalized renal function from either a creatinine-based estimate or a tracer-measured GFR. Darwish 2025b Methods states that GFR was estimated with the Schwartz equation for subjects under 18 years of age but does not name the equation used for adults. Creatinine clearance was screened as a separate covariate in the earlier model and was not retained; GFR is the renal descriptor carried into the final model.",
      source_name        = "GFR"
    ),
    DIS_RETT = list(
      description        = "Rett syndrome disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-Rett cohort)",
      notes              = "1 = patient with Rett syndrome, 0 = all other subjects (Darwish 2025b Table 2 footnote, definition of Rett_i). Decreases CL by 14.6% and increases Vp by 80.5%. Both shifts moved materially when the DAFFODIL data were added: Darwish 2025b Results reports them as the only covariate effects changing by more than 10% relative to Darwish_2025a_trofinetide, by +13.6% (CL, -0.169 to -0.146) and +30.7% (Vp, 0.616 to 0.805). Also selects the disease-cohort residual-error magnitude, which Darwish 2025b pooled across Rett syndrome, fragile X syndrome, and traumatic brain injury. All Rett syndrome participants were female.",
      source_name        = "Rett"
    ),
    DIS_TBI = list(
      description        = "Traumatic brain injury disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-TBI cohort)",
      notes              = "1 = patient with traumatic brain injury, 0 = all other subjects (Darwish 2025b Table 2 footnote, definition of TBI_i). Increases CL by 22.9% and decreases Vp by 75.3%. NOTE: distinct from the register's DIS_BURN_RECENT canonical, whose documented source alias 'TBI' denotes recent burn injury rather than traumatic brain injury. Also selects the disease-cohort residual-error magnitude.",
      source_name        = "TBI"
    ),
    DIS_FXS = list(
      description        = "Fragile X syndrome disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-FXS cohort)",
      notes              = "1 = patient with fragile X syndrome, 0 = all other subjects (Darwish 2025b Table 2 footnote, definition of FXS_i). Increases Vc by 116%, the largest single covariate effect in the model. Also selects the disease-cohort residual-error magnitude. All fragile X syndrome participants were male.",
      source_name        = "FXS"
    ),
    FED = list(
      description        = "Fed-vs-fasted state at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "1 = dose administered in the fed state, 0 = fasted (Darwish 2025b Table 2 footnote, definition of Fed_i). Decreases ka by 9.69% and F1 by 13.3%. Estimated from the dedicated phase 1 food-effect study ACP-2566-006, in which the fed arm followed a high-fat meal (Supplementary Material Table S2); Darwish 2025b does not report the meal composition in the main text, and the effect is pooled with the fed dosing in the phase 3 LAVENDER study, so the general FED canonical applies rather than FED_HIGHFAT.",
      source_name        = "FED"
    ),
    DOSE_18G = list(
      description        = "18 g trofinetide dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any dose level other than 18 g, including the 2-12 g therapeutic weight-banded doses)",
      notes              = "1 = subject received the 18 g supratherapeutic dose, 0 = all other subjects (Darwish 2025b Table 2 footnote, definition of DoseGrp1_i). Decreases F1 by 13.2%. Arises from the thorough-QTc study ACP-2566-008, which gave single 12, 18, or 24 g oral doses (Supplementary Material Table S1); the effect captures the less-than-proportional rise in exposure above the therapeutic range. Mutually exclusive with DOSE_24G. Unchanged from Darwish_2025a_trofinetide.",
      source_name        = "DoseGrp1"
    ),
    DOSE_24G = list(
      description        = "24 g trofinetide dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any dose level other than 24 g, including the 2-12 g therapeutic weight-banded doses)",
      notes              = "1 = subject received the 24 g supratherapeutic dose, 0 = all other subjects (Darwish 2025b Table 2 footnote, definition of DoseGrp2_i). Decreases F1 by 28.4%, roughly twice the 18 g reduction. Arises from the thorough-QTc study ACP-2566-008. Mutually exclusive with DOSE_18G. Unchanged from Darwish_2025a_trofinetide.",
      source_name        = "DoseGrp2"
    ),
    AE_DIARRHEA = list(
      description        = "Concurrent treatment-emergent diarrhea indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diarrhea at the time of the dose record)",
      notes              = "1 = subject experiencing diarrhea, 0 = otherwise. The Darwish 2025b Table 2 footnote explicitly describes Diar_i as a TIME-VARYING indicator, so the value is carried per dose record rather than per subject. Decreases F1 by 15.7%, up from 14.8% in Darwish_2025a_trofinetide. Diarrhea is the most common trofinetide adverse event, which is why the sponsor carried it structurally rather than screening it out.",
      source_name        = "Diar"
    )
  )

  # Covariates that the Darwish trofinetide analyses screened in the stepwise
  # covariate search (Darwish 2025b Methods, "Population Pharmacokinetic Model
  # Development") but did NOT retain in the final model. Documented here so the
  # provenance of the paper's screen survives without declaring covariates that
  # `model()` never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the composite 'sex/disease state' covariate, which Darwish 2025b Methods enumerates as six categories: male healthy volunteers, female healthy volunteers, patients with Rett syndrome (female only), patients with fragile X syndrome (male only), male patients with TBI, and female patients with TBI. Only the disease-state contrasts survived backward elimination, so sex itself is not in the final model; the retained members of that composite are DIS_RETT, DIS_TBI, and DIS_FXS. The analysis population was 56.9% female."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Listed among the covariates re-evaluated in Darwish 2025b Methods but not retained; body weight was the size descriptor carried into the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025b Methods) but not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025b Methods) but not retained in the final model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025b Methods) but not retained in the final model. Darwish 2025b does not report the units in which it was screened."
    ),
    ROUTE_NGT = list(
      description = "Gastric-tube-vs-oral enteral administration route indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Darwish 2025b Methods lists 'route of enteral administration (oral or gastric tube)' among the covariates screened; it was not retained, so oral and gastric-tube doses share the same absorption and bioavailability parameters here. Gastric-tube dosing was common in the Rett syndrome cohorts -- 31-41% of participants in the three Neu-2566-Rett studies and 4 of 13 (30.8%) in DAFFODIL (Darwish 2025b Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 455,
    n_observations = 5709,
    n_studies      = 14,
    age_range      = "2-64 years (mean 21.8)",
    weight_range   = "9.8-140 kg (mean 56.5)",
    sex_female_pct = 56.9,
    gfr_reference  = "124 mL/min/1.73 m^2 (analysis-population median)",
    disease_state  = "Pooled analysis of 156 healthy volunteers, 198 patients with Rett syndrome (female only), 57 patients with traumatic brain injury (male only), and 44 patients with fragile X syndrome (male only).",
    dose_range     = "Oral, gastric-tube, and intravenous (bolus and infusion) trofinetide; oral doses spanned the 2-12 g therapeutic weight-banded range plus supratherapeutic single 18 g and 24 g dose levels.",
    regions        = "Not reported",
    notes          = "Darwish 2025b Results, 'Population Pharmacokinetic Model'. This is the 14-study update of the 13-study model packaged as Darwish_2025a_trofinetide: it adds 114 concentrations from 13 girls aged 2-4 years with Rett syndrome enrolled in the phase 2/3 DAFFODIL study (ACP-2566-009), whose mean age was 3 years (range 2-4) and mean baseline body weight 13.4 kg (range 9.8-18.1). The DAFFODIL subgroup was 92.3% white and 7.7% Asian; race is not reported for the pooled population. Assay: LC-MS/MS in lithium-heparinized whole blood, LLOQ 0.100 ug/mL, ULOQ 100 ug/mL. The target steady-state exposure range used to confirm the weight-banded regimens was AUC0-12 of 800-1200 ug*h/mL. DAFFODIL dosing escalated from 2 g BID to 4 g BID at week 2, then to weight-banded 5 g BID (>= 9 to < 12 kg) or 6 g BID (>= 12 to < 20 kg) at week 4 (Supplementary Material Table S2)."
  )

  ini({
    # Structural parameters and covariate effects -- Darwish 2025b Table 2,
    # "Population mean estimate" column, which holds the UPDATED 14-study
    # estimates. Reference values are the analysis-population medians: WT 58 kg,
    # AGE 22.4 years, GFR 124 mL/min/1.73 m^2 (Table 2 footnote).
    #
    # PROVENANCE WARNING: the "Model equations:" block printed in the Table 2
    # footnote of this paper is a verbatim copy of the earlier 13-study paper's
    # equations and carries that model's coefficients (F1 0.828, ka 0.391,
    # CL 11.8, weight exponent 0.443, Vc 24.9, Vp 35.3, Rett shifts -0.169 and
    # 0.616, diarrhea -0.148). The footnote therefore contradicts the table it
    # annotates. The table is authoritative: Results states the primary
    # parameters differ by < 1.5% between the two models and that the Rett
    # shifts on CL and Vp increased by 13.6% and 30.7%, and only the table
    # column reproduces those figures exactly ((0.169 - 0.146)/0.169 = 13.6%;
    # (0.805 - 0.616)/0.616 = 30.7%). Only the FUNCTIONAL FORM is taken from
    # the footnote equations; every coefficient below comes from the table.
    lcl          <- log(11.7);   label("Clearance at the reference covariate values (L/h)")                          # Table 2: CL = 11.7 L/h (RSE 1.99%)
    e_wt_cl      <- 0.486;       label("Power exponent on (WT/58) for clearance (unitless)")                         # Table 2: Covariate exponent of weight on CL = 0.486 (RSE 7.40%)
    e_crcl_cl    <- 0.272;       label("Power exponent on (CRCL/124) for clearance (unitless)")                      # Table 2: Covariate exponent of GFR on CL = 0.272 (RSE 19.9%)
    e_rett_cl    <- -0.146;      label("Proportional shift in clearance for Rett syndrome (fraction)")               # Table 2: Shift in CL for RTT = 1 is -0.146 (RSE 29.4%)
    e_tbi_cl     <- 0.229;       label("Proportional shift in clearance for traumatic brain injury (fraction)")      # Table 2: Shift in CL for TBI = 1 is 0.229 (RSE 17.8%)

    lvc          <- log(25.0);   label("Central volume of distribution at the reference age (L)")                    # Table 2: Vc = 25.0 L (RSE 3.78%)
    e_age_vc     <- 0.556;       label("Power exponent on (AGE/22.4) for central volume (unitless)")                 # Table 2: Covariate exponent of age on Vc = 0.556 (RSE 7.98%)
    e_fxs_vc     <- 1.16;        label("Proportional shift in central volume for fragile X syndrome (fraction)")     # Table 2: Shift in Vc for FXS = 1 is 1.16 (RSE 20.0%)

    lq           <- log(1.42);   label("Intercompartmental clearance (L/h)")                                         # Table 2: Q = 1.42 L/h (RSE 5.97%)

    lvp          <- log(35.4);   label("Peripheral volume of distribution (L)")                                      # Table 2: Vp = 35.4 L (RSE 5.69%)
    e_rett_vp    <- 0.805;       label("Proportional shift in peripheral volume for Rett syndrome (fraction)")       # Table 2: Shift in Vp for RTT = 1 is 0.805 (RSE 23.4%)
    e_tbi_vp     <- -0.753;      label("Proportional shift in peripheral volume for traumatic brain injury (fraction)") # Table 2: Shift in Vp for TBI = 1 is -0.753 (RSE 4.43%)

    lka          <- log(0.394);  label("First-order absorption rate constant in the fasted state (1/h)")             # Table 2: ka = 0.394 1/h (RSE 4.00%)
    e_fed_ka     <- -0.0969;     label("Proportional shift in ka for the fed state (fraction)")                      # Table 2: Shift in ka for FED = -0.0969 (RSE 23.6%)

    lfdepot      <- log(0.832);  label("Oral bioavailability in the fasted, therapeutic-dose, diarrhea-free reference state (fraction)") # Table 2: F1 = 0.832 (RSE 3.59%)
    e_fed_f      <- -0.133;      label("Proportional shift in bioavailability for the fed state (fraction)")         # Table 2: Shift in F1 for FED = -0.133 (RSE 7.28%)
    e_dose18g_f  <- -0.132;      label("Proportional shift in bioavailability for the 18 g dose level (fraction)")   # Table 2: Shift in F1 for 18 g dose = -0.132 (RSE 20.0%)
    e_dose24g_f  <- -0.284;      label("Proportional shift in bioavailability for the 24 g dose level (fraction)")   # Table 2: Shift in F1 for 24 g dose = -0.284 (RSE 5.65%)
    e_diarrhea_f <- -0.157;      label("Proportional shift in bioavailability during diarrhea (fraction)")           # Table 2: Shift in F1 for diarrhea = -0.157 (RSE 19.5%)

    # IIV -- Darwish 2025b Table 2 reports interindividual variability twice:
    # the "Interindividual variability" block near the bottom holds the NONMEM
    # $OMEGA VARIANCES, and the "Variability / Estimate" column beside each
    # structural parameter holds the corresponding %CV. Methods Eq. (1) defines
    # %CV = sqrt(exp(omega^2) - 1) * 100, and inverting each variance
    # reproduces the printed %CV exactly (sqrt(exp(0.0182) - 1) = 13.6%;
    # 0.0906 -> 30.8%; 0.355 -> 65.3%; 0.0874 -> 30.2%; 0.04 -> 20.2%), which
    # confirms the block is variances rather than SDs. The variances are used
    # directly here because they are the primary reported quantity.
    # No IIV was estimated on ka ("NE" in Table 2), so ka carries no eta.
    etalcl     ~ 0.0182   # Table 2 IIV block: omega^2 for CL = 0.0182, printed as 13.6 %CV (RSE 17.4%; eta-shrinkage 37.5%)
    etalvc     ~ 0.0906   # Table 2 IIV block: omega^2 for Vc = 0.0906, printed as 30.8 %CV (RSE 12.5%; eta-shrinkage 31.9%)
    etalq      ~ 0.355    # Table 2 IIV block: omega^2 for Q = 0.355, printed as 65.3 %CV (RSE 17.7%; eta-shrinkage 9.8%)
    etalvp     ~ 0.0874   # Table 2 IIV block: omega^2 for Vp = 0.0874, printed as 30.2 %CV (RSE 14.3%; eta-shrinkage 41.0%)
    etalfdepot ~ 0.04     # Table 2 IIV block: omega^2 for F1 = 0.04, printed as 20.2 %CV (RSE 16.0%; eta-shrinkage 54.8%)

    # Residual error -- Darwish 2025b Methods Eq. (2) specifies a
    # log/exponential residual error model, equivalent to nlmixr2's `lnorm()`.
    # Two separate magnitudes were estimated, one for healthy subjects and one
    # for subjects with Rett syndrome, FXS, or TBI. The "Population mean
    # estimate" column holds the NONMEM $SIGMA VARIANCE and the adjacent
    # column holds the log-scale SD as a %CV: sqrt(0.0788) = 0.281 -> 28.1 %CV
    # and sqrt(0.140) = 0.374 -> 37.4 %CV, both exact. `lnorm()` takes an SD,
    # so each is entered as sqrt(sigma^2).
    expSdHealthy <- sqrt(0.0788); label("Log-scale residual SD, healthy subjects (log units)")                        # Table 2: Residual variability, healthy subjects, sigma^2 = 0.0788 (RSE 1.53%), printed as 28.1 %CV
    expSdDisease <- sqrt(0.140);  label("Log-scale residual SD, subjects with Rett syndrome, TBI, or FXS (log units)") # Table 2: Residual variability, patients with RTT, TBI or FXS, sigma^2 = 0.140 (RSE 3.19%), printed as 37.4 %CV
  })

  model({
    # Individual PK parameters. The functional form of every covariate term is
    # transcribed from the Darwish 2025b Table 2 footnote equations; each
    # categorical effect enters as a proportional shift of the form
    # (1 + theta * indicator), with theta carrying the sign printed in the
    # Table 2 "Population mean estimate" column.
    cl <- exp(lcl + etalcl) * (WT / 58)^e_wt_cl * (CRCL / 124)^e_crcl_cl *
      (1 + e_tbi_cl * DIS_TBI) * (1 + e_rett_cl * DIS_RETT)
    vc <- exp(lvc + etalvc) * (AGE / 22.4)^e_age_vc * (1 + e_fxs_vc * DIS_FXS)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (1 + e_rett_vp * DIS_RETT) * (1 + e_tbi_vp * DIS_TBI)
    ka <- exp(lka) * (1 + e_fed_ka * FED)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # Bioavailability carries the fed-state, supratherapeutic-dose, and
    # diarrhea effects. Applied to the depot so intravenous dosing (directly
    # into `central`) is unaffected, matching the source dataset in which F1
    # was identifiable only because IV data were pooled with oral data.
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_fed_f * FED) *
      (1 + e_dose18g_f * DOSE_18G) * (1 + e_dose24g_f * DOSE_24G) *
      (1 + e_diarrhea_f * AE_DIARRHEA)
    f(depot) <- fdepot

    # `central` is in mg and `vc` in L, so central/vc is mg/L = ug/mL, the
    # units in which Darwish 2025b reports trofinetide whole-blood
    # concentrations.
    Cc <- central / vc

    # A single residual-error model switched by the disease-state indicators.
    # The three indicators are mutually exclusive, so their sum is 1 for any
    # patient cohort and 0 for healthy volunteers.
    expSd <- expSdDisease * (DIS_RETT + DIS_TBI + DIS_FXS) +
      expSdHealthy * (1 - DIS_RETT - DIS_TBI - DIS_FXS)
    Cc ~ lnorm(expSd)
  })
}
