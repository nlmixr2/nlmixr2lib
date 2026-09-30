Chan_2021_pregabalin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "an absorption lag time for oral pregabalin in pooled pediatric",
    "(3 months to 16 years) and adult data (Chan 2021): healthy adults,",
    "adults with various degrees of renal function, and adult and",
    "pediatric patients with focal onset seizures (10 studies). CL/F is",
    "proportional to body-surface-area-normalised creatinine clearance",
    "(CRCL, mL/min/1.73 m^2) up to an estimated breakpoint of 96.4 and",
    "constant above it, with estimated allometric weight exponents on",
    "CL/F (0.52) and V/F (0.70) and female-sex multipliers on both. ka",
    "is estimated as a multiple of the individual elimination rate",
    "constant CL/V (to avoid flip-flop), with fed and unknown-food-status",
    "effects on ka and a fed effect on the lag time. Residual error is",
    "combined proportional + additive with separate magnitudes for phase",
    "I adult, phase III adult, phase I pediatric (A0081074) and phase III",
    "pediatric (A0081041, PERIWINKLE) studies. Individual predicted",
    "average steady-state concentrations from this model drive the",
    "exposure-response model Chan_2021_pregabalin_lsr28."
  )
  reference <- paste(
    "Chan PLS, Marshall SF, McFadyen L, Liu J.",
    "Pregabalin Population Pharmacokinetic and Exposure-Response Analyses",
    "for Focal Onset Seizures in Children (4-16 years) and Adults, to",
    "Support Dose Recommendations in Children.",
    "Clin Pharmacol Ther. 2021;110(1):132-140. doi:10.1002/cpt.2132.",
    "Covariate equations in Supplementary Information Appendix Equations I;",
    "final NONMEM control stream (run8.mod) in the second supplementary",
    "file."
  )
  vignette <- "Chan_2021_pregabalin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function normalised to 70 kg with estimated exponents on",
        "CL/F (0.52) and V/F (0.70) (Table 2 footnote d; control stream",
        "NORMBWT = 70)."
      ),
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Male is the reference (Table 2 footnote e). The control stream",
        "defines its own local SEXF with the OPPOSITE coding (SEXF = 1 for",
        "males, SEX == 2 -> 0 for females) and applies THETA(13) / THETA(7)",
        "when that local SEXF is 0; the canonical SEXF here is 1 for",
        "females, so the multipliers enter as e_sexf_cl^SEXF and",
        "e_sexf_vc^SEXF."
      ),
      source_name = "SEX"
    ),
    CRCL = list(
      description = paste(
        "Body-surface-area-normalised creatinine clearance (BSA-NCLcr)"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value. Cockcroft-Gault for subjects >= 13 years and the",
        "modified Schwartz equation (K = 0.55; K = 0.45 below 1 year) for",
        "younger children, normalised to 1.73 m^2 of body surface area",
        "(Methods, Population PK model; Table 1 footnote). CL/F is",
        "proportional to CRCL up to the breakpoint crcl_hinge = 96.4",
        "mL/min/1.73 m^2 and constant above it."
      ),
      source_name = "NBCCL"
    ),
    FED = list(
      description = "Fed state at dosing (1 = fed, 0 = fasted or unknown)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Control stream FOOD == 1. Reduces ka (all studies); in phase I",
        "studies (STUDY_PHASE3 = 0) it also lengthens the lag time and",
        "switches the ka random effect to etalka_fed. Known fed data came",
        "from one phase I adult study (1008-003; Table 2 footnote h and",
        "Discussion). Set FED = 0 and FED_MISSING = 0 for the fasted",
        "reference used in the paper's simulations."
      ),
      source_name = "FOOD"
    ),
    FED_MISSING = list(
      description = "Food status at dosing not recorded (1 = unknown)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Control stream FOOD == 2. All phase III adult studies were",
        "collected under unknown food status (Table 2 footnote h). Must be",
        "0 whenever FED = 1."
      ),
      source_name = "FOOD"
    ),
    STUDY_PHASE3 = list(
      description = "Phase III study stratum (1 = phase III, 0 = phase I)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Complement of the control-stream PTS flag (PTS = 1 when",
        "FLAGPHASE == 1). Selects the phase I vs phase III residual-error",
        "magnitudes for adult studies and gates the fed effect on the lag",
        "time and the fed ka random effect (both applied only when FED = 1",
        "and PTS = 1). Phase I = four healthy-adult studies, the renal",
        "impairment study and the pediatric PK study A0081074; phase III =",
        "the three adult focal-onset-seizure studies and PERIWINKLE",
        "(A0081041)."
      ),
      source_name = "FLAGPHASE"
    ),
    STUDY_A0081074 = list(
      description = "Phase I pediatric PK study A0081074 (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Control stream PROT == 1074; overrides the phase residual error",
        "with its own proportional SD (29.8%) and an additive SD fixed to",
        "0 ($SIGMA 8 0 FIX). Pediatric patients aged 1 month to 16 years",
        "with focal onset seizures (reference 10 of the paper)."
      ),
      source_name = "PROT"
    ),
    STUDY_A0081041 = list(
      description = "Phase III pediatric study A0081041 (PERIWINKLE) (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Control stream PROT == 1041; overrides the phase residual error",
        "with its own proportional (35.0%) and additive (0.68 ug/mL) SDs.",
        "Pediatric patients aged 4-16 years with focal onset seizures",
        "(NCT01389596)."
      ),
      source_name = "PROT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened on CL/F and V/F as a power function centred at 32 years",
        "(Appendix Equations I; control stream MEDAGE = 32) but not",
        "retained (Table S1); THETA(11) and THETA(12) are 0 FIX in the",
        "final control stream."
      )
    ),
    RACE = list(
      description = "Race (White reference; Black, Asian, Other)",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Screened on CL/F and V/F as fractional changes (3 degrees of",
        "freedom) but not retained (Table S1); THETA(14) to THETA(19) are",
        "0 FIX in the final control stream."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "pregabalin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "pregabalin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 979L,
    n_studies = 10L,
    age_range = "3 months to 75 years",
    age_median = "10 years (children); 38 years (adults)",
    weight_range = "6.6-180 kg",
    weight_median = "32.9 kg (children); 75.5 kg (adults)",
    sex_female_pct = 49.8,
    race_ethnicity = c(White = 80.2, Black = 5.2, Asian = 6.9, Other = 7.7),
    disease_state = paste(
      "Healthy adults, adults with various degrees of renal function, and",
      "adult and pediatric patients with focal onset seizures"
    ),
    dose_range = paste(
      "Oral pregabalin; adults 150-600 mg/day b.i.d. or t.i.d.; children",
      "2.5-10 mg/kg/day (>= 30 kg) or 3.5-14 mg/kg/day (< 30 kg)"
    ),
    renal_function = paste(
      "CLcr 42.2-261 mL/min in adults (excluding the renal study) and",
      "15.5-293 mL/min in children; NCLcr median 149 mL/min/1.73 m^2 in",
      "children and 101 in adults"
    ),
    regions = "Multinational",
    notes = paste(
      "724 adults and 255 pediatric patients (162 aged 3 months to < 12",
      "years, 93 aged 12-16 years); 5,258 PK samples (Chan 2021 Table 1).",
      "Two pediatric studies, three adult phase III studies, four phase I",
      "healthy-adult studies and one renal-impairment study."
    )
  )

  # Stratum-specific residual SDs, combined into one proportional and one
  # additive symbol inside model() from the study indicators (the
  # Rich_2026_momelotinib.R / Willmann_2019_moxifloxacin.R pattern).
  paper_specific_residual_sds <- c(
    "propSdPh1Adult",
    "addSdPh1Adult",
    "propSdPh3Adult",
    "addSdPh3Adult",
    "propSdPh1Ped",
    "propSdPh3Ped",
    "addSdPh3Ped"
  )
  # Fed-state ka random effect that REPLACES etalka (it is not added to it)
  # for fed phase I records, per the control stream
  # 'IF (FOOD.EQ.1.AND.PTS.EQ.1) KA = TVKA*EXP(ETA(4))'.
  paper_specific_etas <- c("etalka_fed")

  ini({
    # Chan 2021 Table 2 reports every structural row as the typical value
    # for the reference subject (70 kg male, CRCL at or above the
    # breakpoint, fasted), in natural units. The underlying NONMEM THETAs
    # are proportionality factors (Table 2 footnotes c and f), which the
    # model() block reconstructs.
    lcl <- log(4.96)
    label("CL/F for a 70 kg male with CRCL at or above the breakpoint (L/h)") # Table 2 'CL/F' 4.96 [1.78] L/hr (95% CI 4.71-5.20); plateau = THETA(1) * breakpoint
    lcrcl_hinge <- log(96.4)
    label("CRCL breakpoint above which CL/F is constant (mL/min/1.73 m^2)") # Table 2 'CLcr breakpoint' 96.4 [1.91] (95% CI 90.7-104)
    lvc <- log(39.8)
    label("V/F for a 70 kg male (L)") # Table 2 'V/F' 39.8 [1.62] L (95% CI 38.6-41.0)
    lka <- log(10.0)
    label("Fasted ka for the reference subject (1/h)") # Table 2 'ka fasted' 10.0 [16.2]/hr (95% CI 7.44-15.1); ka = THETA(3) * CL/V (footnote f)
    ltlag <- log(0.32)
    label("Absorption lag time, fasted or unknown food status (h)") # Table 2 'Tlag' 0.32 [1.52] hr (95% CI 0.31-0.32)

    e_wt_cl <- 0.52
    label("Allometric exponent of body weight on CL/F (unitless)") # Table 2 'Body weight on CL/F' 0.52 [4.72] (95% CI 0.48-0.57)
    e_wt_vc <- 0.70
    label("Allometric exponent of body weight on V/F (unitless)") # Table 2 'Body weight on V/F' 0.70 [4.59] (95% CI 0.64-0.76)
    e_sexf_cl <- 0.92
    label("CL/F multiplier for females vs males (unitless)") # Table 2 'Sex on CL/F' 0.92 [2.00] (95% CI 0.88-0.95); Results '8% lower CL/F in female patients'
    e_sexf_vc <- 0.83
    label("V/F multiplier for females vs males (unitless)") # Table 2 'Sex on V/F' 0.83 [2.48] (95% CI 0.79-0.87)
    # Table 2 prints the fed and unknown-food ka VALUES (0.71 and 1.22 per
    # hour for the reference subject, the same units as the fasted 10.0
    # row); the control stream estimates them as fractional changes on the
    # fasted ka, FKA*(1 + THETA(6)) and FKA*(1 + THETA(9)). The fractional
    # changes are therefore the printed ratios minus one.
    e_fed_ka <- 0.71 / 10.0 - 1
    label("Fractional change in ka when fed (unitless)") # Table 2 'Food: fed' (ka) 0.71 [2.39] (95% CI 0.33-5.14), as ka_fed / ka_fasted - 1 = -0.929
    e_fed_missing_ka <- 1.22 / 10.0 - 1
    label("Fractional change in ka when food status is unknown (unitless)") # Table 2 'Food: unknown' 1.22 [3.26] (95% CI 0.49-3.57), as ka_unknown / ka_fasted - 1 = -0.878
    e_fed_tlag <- 0.43
    label("Fractional change in Tlag when fed, phase I studies (unitless)") # Table 2 'Food: fed' (Tlag) 0.43 [10.5] (95% CI 0.34-0.83); footnote h

    # IIV: Table 2 percentages are sqrt(omega^2) x 100 (the same convention
    # as the residual rows, whose printed SDs are sqrt($SIGMA) of the
    # control stream), so omega^2 = (CV/100)^2. IIV on Tlag is 0 FIX
    # (Results; $OMEGA 5 0 FIX) and is omitted.
    etalcl ~ 0.040804 # Table 2 IIV 'CL/F' 20.2% [18.7], 0.202^2; shrinkage 27.4%
    etalvc ~ 0.016384 # Table 2 IIV 'V/F' 12.8% [21.7], 0.128^2; shrinkage 60.7%
    etalka ~ 1.3689 # Table 2 IIV 'ka' 117% [13.9], 1.17^2; shrinkage 46.2%
    etalka_fed ~ 0.335241 # Table 2 IIV 'ka: fed' 57.9% [74.4], 0.579^2; shrinkage 89.6%; replaces etalka for fed phase I records

    # Residual error: Y = F*exp(eps_prop) + eps_add ($ERROR), i.e. combined
    # proportional + additive. SDs as printed in Table 2.
    propSdPh1Adult <- 0.166
    label("Proportional residual SD, phase I adult studies (fraction)") # Table 2 'Proportional error, Phase I adult' 16.6% [10.1]
    addSdPh1Adult <- 0.021
    label("Additive residual SD, phase I adult studies (ug/mL)") # Table 2 'Additive error, Phase I adult' 0.021 [66.8] ug/mL
    propSdPh3Adult <- 0.289
    label("Proportional residual SD, phase III adult studies (fraction)") # Table 2 'Proportional error, Phase III adult' 28.9% [7.49]
    addSdPh3Adult <- 0.047
    label("Additive residual SD, phase III adult studies (ug/mL)") # Table 2 'Additive error, Phase III adult' 0.047 [77.3] ug/mL
    propSdPh1Ped <- 0.298
    label("Proportional residual SD, phase I pediatric study A0081074 (fraction)") # Table 2 'Proportional error, Phase I pediatric' 29.8% [22.7]; its additive SD is 0 FIX ($SIGMA 8)
    propSdPh3Ped <- 0.350
    label("Proportional residual SD, phase III pediatric study A0081041 (fraction)") # Table 2 'Proportional error, Phase III pediatric' 35.0% [21.0]
    addSdPh3Ped <- 0.68
    label("Additive residual SD, phase III pediatric study A0081041 (ug/mL)") # Table 2 'Additive error, Phase III pediatric' 0.68 [67.6] ug/mL
  })

  model({
    # ----- Covariate terms (Appendix Equations I; control stream $PK) -----
    crcl_hinge <- exp(lcrcl_hinge)
    # CL/F proportional to CRCL up to the breakpoint, constant above it
    # (TVCL = THETA(1)*CLCR if CLCR <= THETA(5), else THETA(1)*THETA(5)).
    # exp(lcl) is the plateau THETA(1)*THETA(5), so the ratio below is
    # the same piecewise-linear relationship.
    renal_cl <- min(CRCL, crcl_hinge) / crcl_hinge
    # Fed effects on Tlag and the ka random effect apply only to fed
    # records in phase I studies (FOOD == 1 and PTS == 1).
    fed_ph1 <- FED * (1 - STUDY_PHASE3)

    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl) * renal_cl * (WT / 70)^e_wt_cl * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * e_sexf_vc^SEXF
    kel <- cl / vc

    # ka is a multiple of the INDIVIDUAL elimination rate constant
    # (KA = CL/V * THETA(3) * food factor * exp(eta)). Table 2 prints ka
    # for the reference subject, so the multiple is exp(lka) / kel_ref.
    kel_ref <- exp(lcl) / exp(lvc)
    ka <- exp(lka + etalka * (1 - fed_ph1) + etalka_fed * fed_ph1) *
      (1 + e_fed_ka * FED) *
      (1 + e_fed_missing_ka * FED_MISSING) *
      kel / kel_ref
    tlag <- exp(ltlag) * (1 + e_fed_tlag * fed_ph1)

    # ----- ODEs -----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # ----- Observation and residual error -----
    # Pediatric study overrides (PROT == 1074 / 1041) take precedence over
    # the phase I / phase III adult strata, as in the $ERROR block.
    ped <- STUDY_A0081074 + STUDY_A0081041
    propSdCc <- (1 - ped) *
      (STUDY_PHASE3 * propSdPh3Adult + (1 - STUDY_PHASE3) * propSdPh1Adult) +
      STUDY_A0081074 * propSdPh1Ped +
      STUDY_A0081041 * propSdPh3Ped
    addSdCc <- (1 - ped) *
      (STUDY_PHASE3 * addSdPh3Adult + (1 - STUDY_PHASE3) * addSdPh1Adult) +
      STUDY_A0081041 * addSdPh3Ped
    Cc <- central / vc
    Cc ~ add(addSdCc) + prop(propSdCc)
  })
}
