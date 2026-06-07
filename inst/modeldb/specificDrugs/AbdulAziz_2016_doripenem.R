AbdulAziz_2016_doripenem <- function() {
  description <- "Two-compartment IV population PK model for doripenem in 12 Malaysian critically ill adults with sepsis receiving 500 mg as a 1-hour infusion every 8 hours (Abdul-Aziz 2016). Reported on free (unbound) doripenem; observed total concentrations were corrected by multiplying by 0.90 to account for ~10% protein binding. Body-weight allometric scaling is fixed (0.75 on CL/Q, 1 on V1/V2, reference 70 kg); Cockcroft-Gault creatinine clearance has an exponential effect on CL centred at the cohort mean 82.5 mL/min."
  reference <- paste(
    "Abdul-Aziz MH, Abd Rahman AN, Mat-Nor M-B, Sulaiman H,",
    "Wallis SC, Lipman J, Roberts JA, Staatz CE.",
    "Population pharmacokinetics of doripenem in critically ill patients",
    "with sepsis in a Malaysian intensive care unit.",
    "Antimicrob Agents Chemother. 2016;60(1):206-214.",
    "doi:10.1128/AAC.01543-15.",
    sep = " "
  )
  vignette <- "AbdulAziz_2016_doripenem"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling fixed a priori at 0.75 on clearances (CL, Q)",
        "and 1 on volumes (V1, V2), all standardised to a body weight of",
        "70 kg (Abdul-Aziz 2016 Methods, citing Anderson and Holford 2008).",
        "Source column not explicitly named in the paper; total body",
        "weight is referred to as WT throughout. BMI range 16.7-29.7 kg/m^2",
        "(Table 1) implies approximately 45-90 kg for the recruited cohort."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Estimated on each sampling occasion using the Cockcroft-Gault",
        "formula (Abdul-Aziz 2016 Methods). Cohort range 30-161 mL/min,",
        "median 85 mL/min, IQR 49-118 mL/min (Table 1). Exponential",
        "effect on CL centred at the population mean 82.5 mL/min:",
        "CL = theta_CLpop * exp(theta_CLCR * (CLCR - 82.5)) with",
        "theta_CLCR = 0.014 (Equation 1; Table 2 final model). A 30 mL/min",
        "increase in CLCR therefore increases CL by exp(0.014*30)-1 = 52%.",
        "Stored under canonical CRCL with raw mL/min (the source-aliases",
        "entry in inst/references/covariate-columns.md permits the raw",
        "Cockcroft-Gault form when the source paper does not apply",
        "BSA-normalisation)."
      ),
      source_name        = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 12L,
    n_studies        = 1L,
    age_range        = "19-74 years",
    age_median       = "52 years (IQR 33-60)",
    weight_range     = "Not directly reported; BMI range 16.7-29.7 kg/m^2 (Table 1)",
    weight_median    = "Not directly reported; BMI median 22.9 kg/m^2 (Table 1)",
    sex_female_pct   = 17,
    race_ethnicity   = "Malaysian (single-centre ICU in Pahang, Malaysia)",
    disease_state    = paste(
      "Critically ill adult ICU patients (>=18 years) with sepsis,",
      "two thirds with ventilator-associated pneumonia and the remainder",
      "with intra-abdominal sepsis; all required invasive mechanical",
      "ventilation and vasopressor support. Patients on extracorporeal",
      "renal support were excluded."
    ),
    dose_range       = "500 mg doripenem IV infusion over 1 hour, every 8 hours",
    regions          = "Malaysia (single tertiary-hospital adult ICU in Pahang)",
    n_concentrations = 140L,
    apache_ii_median = "19 (IQR 15-22)",
    sofa_median      = "6 (IQR 5-7)",
    crcl_range       = "30-161 mL/min (Cockcroft-Gault; median 85, IQR 49-118; cohort mean 82.5 mL/min used as the model reference)",
    notes            = paste(
      "Prospective open-label PK study (Abdul-Aziz 2016 Table 1), November",
      "2012 to October 2013. 140 plasma doripenem concentrations across",
      "two sampling occasions (day 1 and day 3 of therapy); samples drawn",
      "before dose then at 0.5, 0.75, 1, 1.5, 2, 4 and 8 h after the start",
      "of the infusion. Doripenem assay: validated HPLC-UV (LLOQ 0.2 mg/L,",
      "linear 0.2-100 mg/L). Observed total doripenem concentrations were",
      "corrected for ~10% protein binding by multiplying by 0.90 prior to",
      "modelling, so the reported parameters describe free doripenem.",
      "Causative pathogens: Acinetobacter baumannii (4), Pseudomonas",
      "aeruginosa (1), nil identified (7). Three patients died during the",
      "ICU stay (deaths attributed to preexisting comorbidities). Model",
      "fitted with NONMEM 7.3 using FOCE-I; allometric scaling fixed a",
      "priori (Methods)."
    )
  )

  ini({
    # ===== Structural PK (Abdul-Aziz 2016 Table 2 'Final model' column) =====
    # Reference subject: WT = 70 kg, CLCR = 82.5 mL/min.
    lcl <- log(10.1); label("Typical CL at WT=70 kg, CLCR=82.5 mL/min (L/h)")  # Abdul-Aziz 2016 Table 2 final model: CL = 10.1 L/h/70 kg (bootstrap median 9.9, 95% CI 8.9-10.9)
    lvc <- log(15.5); label("Typical central volume V1 at WT=70 kg (L)")        # Abdul-Aziz 2016 Table 2 final model: V1 = 15.5 L/70 kg (bootstrap median 15.8, 95% CI 10.5-22.0)
    lvp <- log(17.7); label("Typical peripheral volume V2 at WT=70 kg (L)")     # Abdul-Aziz 2016 Table 2 final model: V2 = 17.7 L/70 kg (bootstrap median 18.1, 95% CI 12.2-26.6)
    lq  <- log(36.3); label("Typical intercompartmental clearance Q at WT=70 kg (L/h)")  # Abdul-Aziz 2016 Table 2 final model: Q = 36.3 L/h/70 kg (bootstrap median 36.8, 95% CI 28.4-41.5)

    # ===== Allometric / covariate fixed coefficients =====
    # Allometric exponents fixed a priori to canonical values (Methods,
    # citing Anderson and Holford 2008): 0.75 on clearances, 1 on volumes,
    # standardised to 70 kg.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless, fixed)")  # Abdul-Aziz 2016 Methods, allometric scaling
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on (WT/70) for V1 and V2 (unitless, fixed)") # Abdul-Aziz 2016 Methods, allometric scaling

    # Exponential effect of Cockcroft-Gault CLCR on CL, centred at the
    # population mean 82.5 mL/min (Abdul-Aziz 2016 Equation 1, Table 2):
    #   CL = theta_CLpop * (WT/70)^0.75 * exp(theta_CLCR * (CLCR - 82.5))
    e_crcl_cl   <- 0.014; label("Exponential coefficient on (CRCL - 82.5) for CL (per mL/min)")  # Abdul-Aziz 2016 Table 2 final model: theta_CLCR = 0.014 (bootstrap median 0.014, 95% CI 0.012-0.017)
    crcl_ref_cl <- 82.5;  label("Reference Cockcroft-Gault CLCR for the exponential CL covariate (mL/min, population mean)")  # Abdul-Aziz 2016 Results, "normalized to the population mean value of 82.5 ml/min"

    # ===== Inter-individual variability (Abdul-Aziz 2016 Table 2 final model) =====
    # Exponential (log-normal) IIV reported as %CV; omega^2 = log(1 + CV^2).
    # The paper additionally reports between-occasion variability (BOV) on
    # CL of 22.2%. nlmixr2lib has no idiomatic encoding for BOV separate
    # from BSV; per the convention used in Bienczak_2016_nevirapine.R and
    # Svensson_2018_bedaquiline.R, BOV is dropped when a BSV term is
    # reported on the same parameter. The dropped BOV is documented in the
    # vignette's Assumptions and deviations section.
    etalcl ~ 0.01076  # Abdul-Aziz 2016 Table 2 final model: BSV CL = 10.4%; log(1 + 0.104^2) = 0.01076
    etalvc ~ 0.30491  # Abdul-Aziz 2016 Table 2 final model: BSV V1 = 59.7%; log(1 + 0.597^2) = 0.30491
    etalvp ~ 0.43503  # Abdul-Aziz 2016 Table 2 final model: BSV V2 = 73.8%; log(1 + 0.738^2) = 0.43503

    # ===== Residual error (Abdul-Aziz 2016 Table 2 final model) =====
    # Combined additive + proportional residual error model.
    propSd <- 0.090; label("Proportional residual error (fraction)")  # Abdul-Aziz 2016 Table 2 final model: 9.0% CV proportional
    addSd  <- 1.25;  label("Additive residual error (mg/L)")           # Abdul-Aziz 2016 Table 2 final model: additive SD = 1.25 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # CL: allometric on WT plus exponential effect of Cockcroft-Gault CLCR
    # centred at the cohort mean 82.5 mL/min (Equation 1).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * exp(e_crcl_cl * (CRCL - crcl_ref_cl))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # IV doripenem infusion (administered as a 1-hour infusion in the
    # source study); dose is delivered directly into the central
    # compartment via NONMEM ADVAN3 in the source.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    # Plasma free doripenem concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
