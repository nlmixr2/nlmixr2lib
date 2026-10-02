Sakamoto_2021_fluconazole <- function() {
  description <- paste(
    "One-compartment population PK model with first-order oral absorption and first-order",
    "elimination for prophylactic oral fluconazole in Japanese adults with hematological",
    "malignancy receiving chemotherapy or hematopoietic stem cell transplantation. Apparent",
    "clearance scales as a power of Cockcroft-Gault creatinine clearance (normalized to 5.2 L/h)",
    "and apparent volume as a power of body weight (normalized to 57.6 kg); log-normal",
    "between-subject variability on clearance only.",
    sep = " "
  )
  reference <- paste(
    "Sakamoto Y, Isono H, Enoki Y, Taguchi K, Miyazaki T, Kunimoto H, Koike H, Hagihara M,",
    "Matsumoto K, Nakajima H, Sahashi Y, Matsumoto K. Population Pharmacokinetic Analysis and",
    "Dosing Optimization of Prophylactic Fluconazole in Japanese Patients with Hematological",
    "Malignancy. J Fungi (Basel). 2021;7(11):975. doi:10.3390/jof7110975.",
    sep = " "
  )
  vignette <- "Sakamoto_2021_fluconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods 2.3: CLcr 'calculated using ... the Cockcroft-Gault formula'; Table 1 row",
        "'CG-CLcr (mL/min)' median 87.1 (range 31.1-193.6). No BSA normalization is applied (the",
        "paper tabulates the BSA-normalized eGFRcre separately). Stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md, which accepts raw Cockcroft-Gault mL/min",
        "when the source does not BSA-normalize. The paper enters it into CL/F in L/h, normalized",
        "to 5.2 L/h (Results 3.3: 'normalized to the population median of 5.2 L/h' = 87.1 mL/min",
        "x 0.06); the model converts mL/min to L/h with the factor 0.06. Collected on the first day",
        "of fluconazole administration (Methods 2.3). The Japanese creatinine eGFR (eGFRcre,",
        "Matsuo equation; median 69.2 mL/min/1.73 m^2, range 31.6-157.5) was screened as a",
        "competing renal-function covariate and not retained; it is not a model input.",
        sep = " "
      ),
      source_name = "CLcr"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales V/F as (WT/57.6)^e_wt_vc; 57.6 kg is the cohort median (Table 1; Results 3.3",
        "'body weight (normalized to 57.6 kg)'). Range 39.8-99.1 kg. Recorded on the first day",
        "of fluconazole administration.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Screened in the stepwise forward-inclusion / backward-elimination search (Methods 2.5)
  # and not retained in the final model (Results 3.3 keeps only CLcr on CL/F and body weight
  # on V/F). No point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained. Median 53 years, range 20-77 (Table 1)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained. Median 162.6 cm, range 143.6-187.7 (Table 1)."
    ),
    BSA = list(
      description = "Body surface area (Du Bois)",
      units = "m^2",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained. Median 1.62 m^2, range 1.31-2.23 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained. Median 21.7 kg/m^2, range 15.8-29.1 (Table 1)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 1L,
    n_observations = 119L,
    age_range = "20-77 years",
    age_median = "53 years",
    weight_range = "39.8-99.1 kg",
    weight_median = "57.6 kg",
    sex_female_pct = 100 * 23 / 54,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Hematological malignancy receiving chemotherapy (38) or hematopoietic stem cell",
      "transplantation (autologous PBSCT 9, BMT 2, cord blood 2, allogeneic PBSCT 1, other 2).",
      "Diagnoses: non-Hodgkin lymphoma 28, AML/MDS 10, ALL 6, multiple myeloma 3, Hodgkin",
      "lymphoma 2, other 5 (Table 1).",
      sep = " "
    ),
    dose_range = "200 mg oral fluconazole once daily (prophylaxis); all patients",
    regions = "Japan (Yokohama City University Hospital, single centre)",
    renal_function = paste(
      "Cockcroft-Gault CLcr median 87.1 mL/min (range 31.1-193.6); eGFRcre median 69.2",
      "mL/min/1.73 m^2 (range 31.6-157.5) (Table 1).",
      sep = " "
    ),
    notes = paste(
      "Demographics from Table 1 (31 male / 23 female). Enrolment November 2018 to March 2020;",
      "age >= 16 years; critically ill patients excluded (Methods 2.2). 1-4 samples per patient",
      "at trough (1 h pre-dose), 2, 4 or 12 h post-dose; 125 samples collected, 119 analysed",
      "after excluding 6 where renal function changed substantially (Results 3.3). HPLC-UV",
      "assay, calibration 1.56-50 ug/mL. Estimation in Phoenix NLME (FOCE-ELS).",
      sep = " "
    )
  )

  ini({
    # Final-model estimates, Table 2 'Final Model / Estimate' column. Bootstrap medians
    # (n = 1000) in the same table agree to within 20% and are not used.
    lka <- log(0.34); label("Absorption rate constant ka (1/h)") # Table 2, ka = 0.34 /h (SE 0.12, CV 35.1%)
    lcl <- log(1.03); label("Apparent clearance CL/F coefficient theta1 at CLcr = 5.2 L/h, before the printed exp(0.16) factor (L/h)") # Table 2, theta1 = 1.03 L/h (SE 0.08, CV 7.38%)
    lvc <- log(62.3); label("Apparent volume V/F at 57.6 kg (L)") # Table 2, theta3 = 62.3 L (SE 8.87, CV 14.2%)

    # Table 2 prints the unit '(L/h)' against theta2 and '(L)' against theta4, but both are
    # power exponents (CL/F = theta1 x (CLcr/5.2)^theta2; Vd/F = theta3 x (BW/57.6)^theta4) and
    # are therefore unitless.
    e_crcl_cl <- 1.05; label("Power exponent of (CLcr [L/h] / 5.2) on CL/F (unitless)") # Table 2, theta2 = 1.05 (SE 0.14, CV 13.7%)
    e_wt_vc <- 1.06; label("Power exponent of (WT / 57.6) on V/F (unitless)") # Table 2, theta4 = 1.06 (SE 0.36, CV 34.4%)

    # IIV on CL/F only; Table 2 carries no variability term on V/F or ka.
    # Not printed as a variance: back-solved by the maintainers from the paper's own Monte
    # Carlo PTA curves (Figure 4, steady state, BW 60 kg, six CLcr panels; fitted omega^2 =
    # 0.147, fitted log median shift 0.170) and cross-checked on Figure 3 (Day 1). The value
    # 0.16 is the number printed in the CL/F equation of Table 2 ('x e^0.16'); see the
    # vignette 'Assumptions and deviations' for the reading.
    etalcl ~ 0.16

    # Residual error: Methods 2.5 tested additive, multiplicative and combined models, but the
    # final model and its magnitude are not reported (Table 2). Both components are encoded
    # and fixed to zero so that simulations are IPRED-only, like the paper's PTA analysis.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
    addSd <- fixed(0); label("Additive residual SD (mg/L; not reported in the source)")
  })

  model({
    # Table 2: CL/F (L/h) = theta1 x (CLcr/5.2)^theta2 x e^0.16, with CLcr in L/h (footnote
    # 'CLcr, estimated CLcr using the Cockcroft-Gault equation (L/h)'); CRCL is supplied in
    # mL/min, so CRCL * 0.06 converts it to L/h. The printed exp(0.16) factor is kept
    # verbatim: it reproduces the Results 3.3 'median ... CL/F = 1.2 L/h' (1.03 x exp(0.16) =
    # 1.21) and the median clearance implied by the paper's own PTA simulations (Figures 3-4,
    # Table 3).
    cl <- exp(lcl + etalcl) * exp(0.16) * (CRCL * 0.06 / 5.2)^e_crcl_cl
    # Table 2: Vd/F (L) = theta3 x (BW/57.6)^theta4
    vc <- exp(lvc) * (WT / 57.6)^e_wt_vc
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Plasma fluconazole in mg/L (= ug/mL, the paper's unit)
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
