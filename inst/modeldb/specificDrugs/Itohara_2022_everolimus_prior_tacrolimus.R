Itohara_2022_everolimus_prior_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral everolimus in Japanese adult liver transplant recipients who were on tacrolimus before starting everolimus (Itohara 2022, model with the prior tacrolimus dose-per-trough ratio). Identical in form to Itohara_2022_everolimus (allometric body weight, total daily everolimus dose, eGFR, fluconazole and female sex on apparent clearance; body weight on the apparent volume) with an added power effect of the pre-everolimus tacrolimus dose-per-trough (D/C) ratio on apparent clearance, all parameters re-estimated. Everolimus is assayed in whole blood and concentrations are in ng/mL."
  reference <- paste(
    "Itohara K, Yano I, Nakagawa S, Sugimoto M, Hirai M, Yonezawa A, Imai S,",
    "Nakagawa T, Hira D, Ito T, Hata K, Hatano E, Terada T, Matsubara K.",
    "Population pharmacokinetics of everolimus in adult liver transplant",
    "patients: Comparison to tacrolimus disposition and extrapolation to",
    "pediatrics. Clin Transl Sci. 2022;15(11):2652-2662.",
    "doi:10.1111/cts.13389.",
    sep = " "
  )
  vignette <- "Itohara_2022_everolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "everolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL/F (fixed exponent 0.75) and Vd/F (fixed",
        "exponent 1), normalised to 57 kg, the adult cohort median (Table 1)."
      ),
      source_name = "BW"
    ),
    DOSE_EVEROLIMUS_MGD = list(
      description = "Current total daily everolimus dose.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-record covariate (TDM effect), as in Itohara_2022_everolimus.",
        "Enters CL/F as (TDD / 1.0)^0.338 (Results, equation with the",
        "tacrolimus D/C ratio). For a b.i.d. regimen amt = DOSE_EVEROLIMUS_MGD / 2."
      ),
      source_name = "TDD"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole use (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluconazole)",
      notes = "CL/F multiplied by 0.803^FLCZ (Results, equation with the tacrolimus D/C ratio).",
      source_name = "FLCZ"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (Japanese equation).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Japanese Society of Nephrology creatinine eGFR (Methods). Enters",
        "CL/F as (eGFR / 55.2)^0.293 (Results, equation with the tacrolimus",
        "D/C ratio)."
      ),
      source_name = "eGFR"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "CL/F multiplied by 1.23^Female; paper polarity matches SEXF.",
      source_name = "Female"
    ),
    TAC_DOSE_TROUGH_RATIO_BL = list(
      description = "Tacrolimus dose-per-trough-concentration (D/C) ratio just before everolimus was started.",
      units = "mg/(ng/mL)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tacrolimus dose divided by the tacrolimus whole-blood trough",
        "concentration observed just before the first everolimus dose",
        "(Methods 'Correlation between estimated individual CL/F ...'), in",
        "mg/(ng/mL) as printed in the Discussion. The paper does not say",
        "whether the dose is the total daily dose or a single administration;",
        "a mean of 0.49 is consistent with the total daily dose.",
        "Baseline value, constant within a patient. Enters CL/F as",
        "(D/C_Tac / 0.44)^0.230 (Results); 0.44 is the normalising value",
        "printed in the equation. The 80 patients analysed had mean 0.49",
        "(SD 0.30) mg/(ng/mL) before everolimus (Discussion)."
      ),
      source_name = "D/C Tac"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "20-77 years (whole adult cohort)",
    weight_range = "34.5-91.0 kg (whole adult cohort)",
    race_ethnicity = "Japanese (single-centre cohort, Kyoto University Hospital).",
    disease_state = paste(
      "Adult liver transplant recipients who were receiving tacrolimus at",
      "the start of everolimus (the 80-patient subset of the 87-patient adult",
      "cohort used for Figure 4)."
    ),
    dose_range = "Oral everolimus, TDM-guided to a 3-8 ng/mL trough.",
    regions = "Japan.",
    notes = paste(
      "Subset demographics are not tabulated separately; Table 1 gives the",
      "whole 87-patient adult cohort. Tacrolimus D/C ratio mean 0.49 (SD 0.30)",
      "mg/(ng/mL) before and 0.51 (SD 0.32) after everolimus start",
      "(Discussion)."
    )
  )

  ini({
    # Itohara 2022 Results 'Relationship between estimated individual CL/F of
    # everolimus and D/C ratio of tacrolimus':
    #   CL/F (L/h) = 7.04 * (BW/57)^0.75 * (TDD/1.0)^0.338 * 0.803^FLCZ *
    #                (eGFR/55.2)^0.293 * (D/C_Tac/0.44)^0.230 * 1.23^Female
    #   Vd/F (L)   = 323 * (BW/57)
    # No standard errors are reported for this model.
    lka <- fixed(log(6.07)); label("Absorption rate constant ka (1/h)") # Methods: ka fixed at 6.07 /h from literature ref. 8 (as in the Table 2 model)
    lcl <- log(7.04); label("Apparent oral clearance CL/F for a 57-kg male, TDD 1 mg/day, eGFR 55.2, tacrolimus D/C 0.44, no fluconazole (L/h)") # Results equation: 7.04
    lvc <- log(323); label("Apparent volume of distribution Vd/F for a 57-kg patient (L)") # Results equation: 323

    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/57) on CL/F (unitless)") # Results equation: (BW/57)^0.75; Methods: fixed
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/57) on Vd/F (unitless)") # Results equation: (BW/57); Methods: fixed

    e_dose_everolimus_cl <- 0.338; label("Power exponent of (DOSE_EVEROLIMUS_MGD / 1.0 mg/day) on CL/F (unitless)") # Results equation: (TDD/1.0)^0.338
    e_conmed_fluconazole_cl <- 0.803; label("Multiplicative factor on CL/F with concomitant fluconazole (unitless)") # Results equation: 0.803^FLCZ
    e_crcl_cl <- 0.293; label("Power exponent of (CRCL / 55.2) on CL/F (unitless)") # Results equation: (eGFR/55.2)^0.293
    e_tac_dose_trough_ratio_bl_cl <- 0.230; label("Power exponent of (TAC_DOSE_TROUGH_RATIO_BL / 0.44) on CL/F (unitless)") # Results equation: (D/C Tac/0.44)^0.230
    e_sexf_cl <- 1.23; label("Multiplicative factor on CL/F for female sex (unitless)") # Results equation: 1.23^Female

    etalcl ~ 0.05325 # Results: IIV for CL/F decreased to 23.4%; omega^2 = log(0.234^2 + 1) = 0.05325

    propSd <- 0.182; label("Proportional residual error (fraction)") # Results: proportional error 18.2%
    addSd <- 0.510; label("Additive residual error (ng/mL)") # Results: additive error 0.510 ng/mL
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT / 57)^e_wt_cl *
      (DOSE_EVEROLIMUS_MGD / 1.0)^e_dose_everolimus_cl *
      e_conmed_fluconazole_cl^CONMED_FLUCONAZOLE *
      (CRCL / 55.2)^e_crcl_cl *
      (TAC_DOSE_TROUGH_RATIO_BL / 0.44)^e_tac_dose_trough_ratio_bl_cl *
      e_sexf_cl^SEXF
    vc <- exp(lvc) * (WT / 57)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
