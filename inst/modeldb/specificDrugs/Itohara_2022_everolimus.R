Itohara_2022_everolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral everolimus in Japanese adult liver transplant recipients (Itohara 2022 final model). Apparent oral clearance scales allometrically with body weight (fixed exponent 0.75) and carries power effects of the total daily everolimus dose (a therapeutic-drug-monitoring artefact, not nonlinearity) and of eGFR, plus multiplicative effects of concomitant fluconazole and female sex; the apparent volume scales linearly with body weight. The absorption rate constant is fixed from the literature. Everolimus is assayed in whole blood and concentrations are in ng/mL."
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

  # Everolimus was measured in whole blood by electrochemiluminescence
  # immunoassay (Methods 'Patients and data collection'), so `central` is a
  # whole-blood compartment.
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
        "exponent 1), normalised to 57 kg, the adult cohort median (Table 1:",
        "57.0 kg, range 34.5-91.0). Methods 'PPK modeling' fixes both",
        "exponents from the literature."
      ),
      source_name = "BW"
    ),
    DOSE_EVEROLIMUS_MGD = list(
      description = "Current total daily everolimus dose.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-record covariate: the total daily everolimus dose in force,",
        "updated whenever the dose is changed (0.5 mg b.i.d. -> 1.0). Enters",
        "CL/F as the power term (TDD / 1.0)^0.310 (Table 2). The authors",
        "call this a 'TDM effect': doses were titrated to a 3-8 ng/mL trough,",
        "so fast clearers were prescribed more drug, and the association runs",
        "from clearance to dose -- it does not imply nonlinear PK (Discussion).",
        "Must be kept consistent with the dosing records when simulating:",
        "for a b.i.d. regimen amt = DOSE_EVEROLIMUS_MGD / 2. Adult median",
        "1.0 mg/day, range 0.25-4.0 (Table 1)."
      ),
      source_name = "TDD"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole use (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluconazole)",
      notes = paste(
        "CL/F multiplied by 0.823^FLCZ, i.e. 17.7% lower with fluconazole",
        "(Table 2 theta3). 21 of 87 adults (24.1%) received fluconazole",
        "during the observation period (Table 1)."
      ),
      source_name = "FLCZ"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (Japanese equation).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Japanese Society of Nephrology creatinine eGFR,",
        "eGFR = 194 * Scr^-1.094 * Age^-0.287 (* 0.739 if female) with Scr in",
        "mg/dL (Methods 'Patients and data collection'). Enters CL/F as",
        "(eGFR / 55.2)^0.258 (Table 2), 55.2 mL/min/1.73 m^2 being the adult",
        "cohort median (Table 1, range 11.4-223.7)."
      ),
      source_name = "eGFR"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "CL/F multiplied by 1.23^Female (Table 2 theta5). The paper codes",
        "Female = 1 for women, which matches the canonical SEXF polarity, so",
        "no transformation is needed. 46 of 87 adults were female (Table 1)."
      ),
      source_name = "Female"
    )
  )

  # Screened in the stepwise covariate search (Methods 'PPK modeling'; Table S1)
  # but not retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Adult median 60 years, range 20-77 (Table 1).",
      source_name = "Age"
    ),
    HCT = list(
      description = "Hematocrit.",
      units = "%",
      type = "continuous",
      notes = "Screened on CL/F and Vd/F; not retained. Adult median 33.0%, range 18.8-51.4 (Table 1).",
      source_name = "Hematocrit"
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Adult median 26 U/L, range 10-177 (Table 1).",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Adult median 23 U/L, range 4-380 (Table 1).",
      source_name = "ALT"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL/F and Vd/F; not retained. Adult median 3.8 g/dL, range 1.6-5.1 (Table 1).",
      source_name = "Albumin"
    ),
    CONMED_CSA = list(
      description = "Concomitant cyclosporine use.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL/F; not retained. Only 5 of 87 adults received",
        "cyclosporine, which the authors give as the reason the known",
        "cyclosporine interaction could not be confirmed (Discussion)."
      ),
      source_name = "Cyclosporine"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 1L,
    n_observations = "937 everolimus whole-blood concentrations (Results 'Patient characteristics').",
    age_range = "20-77 years",
    age_median = "60 years",
    weight_range = "34.5-91.0 kg",
    weight_median = "57.0 kg",
    sex_female_pct = 52.9,
    race_ethnicity = "Japanese (single-centre cohort, Kyoto University Hospital).",
    disease_state = paste(
      "Adult (> 16 years) liver transplant recipients started on everolimus,",
      "median 1865 post-operative days (range 29-8092) at the first sample.",
      "Most also received a calcineurin inhibitor (tacrolimus or",
      "cyclosporine) and mycophenolate mofetil."
    ),
    renal_function = "eGFR median 55.2 mL/min/1.73 m^2, range 11.4-223.7 (Table 1).",
    co_medication = "Fluconazole 21 (24.1%), prednisolone 40 (46.0%), nifedipine 6 (6.9%), cyclosporine 5 (5.75%) (Table 1).",
    dose_range = paste(
      "Oral everolimus, usually 0.5 or 1.0 mg twice daily at start, then",
      "TDM-guided to a 3-8 ng/mL trough; total daily dose median 1.0 mg/day,",
      "range 0.25-4.0 (Table 1)."
    ),
    regions = "Japan.",
    notes = paste(
      "Retrospective routine-TDM data, March 2018 to December 2020 (Table 1",
      "gives baseline demographics). Almost all samples were troughs: no",
      "sample was taken within 10 h of a dose (Results). Dosing was assumed",
      "at 9:00/21:00. A further 13 children (373 concentrations) were used",
      "only for external evaluation of the model."
    )
  )

  ini({
    # Final model, Itohara 2022 Table 2:
    #   CL/F = theta1 * (BW/57)^0.75 * (TDD/1.0)^theta2 * theta3^FLCZ *
    #          (eGFR/55.2)^theta4 * theta5^Female
    #   Vd/F = theta6 * (BW/57)
    lka <- fixed(log(6.07)); label("Absorption rate constant ka (1/h)") # Table 2 ka = 6.07 (Fixed); Methods: fixed from literature ref. 8
    lcl <- log(6.89); label("Apparent oral clearance CL/F for a 57-kg male, TDD 1 mg/day, eGFR 55.2, no fluconazole (L/h)") # Table 2 theta1 = 6.89 L/h (RSE 4.63%)
    lvc <- log(312); label("Apparent volume of distribution Vd/F for a 57-kg patient (L)") # Table 2 theta6 = 312 L (RSE 9.65%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/57) on CL/F (unitless)") # Table 2 equation; Methods: fixed at 0.75 from literature ref. 14
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/57) on Vd/F (unitless)") # Table 2 equation; Methods: fixed at 1 from literature ref. 14

    e_dose_everolimus_cl <- 0.310; label("Power exponent of (DOSE_EVEROLIMUS_MGD / 1.0 mg/day) on CL/F (unitless)") # Table 2 theta2 (TDM effect) = 0.310 (RSE 15.6%)
    e_conmed_fluconazole_cl <- 0.823; label("Multiplicative factor on CL/F with concomitant fluconazole (unitless)") # Table 2 theta3 = 0.823 (RSE 5.13%)
    e_crcl_cl <- 0.258; label("Power exponent of (CRCL / 55.2) on CL/F (unitless)") # Table 2 theta4 = 0.258 (RSE 24.7%)
    e_sexf_cl <- 1.23; label("Multiplicative factor on CL/F for female sex (unitless)") # Table 2 theta5 = 1.23 (RSE 6.24%)

    # IIV on CL/F only; IIV on Vd/F was dropped for 49% shrinkage (Results).
    etalcl ~ 0.07085 # Table 2 'IIV for CL' 27.1 CV%; omega^2 = log(0.271^2 + 1) = 0.07085

    propSd <- 0.189; label("Proportional residual error (fraction)") # Table 2 proportional error 18.9 CV%
    addSd <- 0.495; label("Additive residual error (ng/mL)") # Table 2 additive error 0.495 ng/mL
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT / 57)^e_wt_cl *
      (DOSE_EVEROLIMUS_MGD / 1.0)^e_dose_everolimus_cl *
      e_conmed_fluconazole_cl^CONMED_FLUCONAZOLE *
      (CRCL / 55.2)^e_crcl_cl *
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
