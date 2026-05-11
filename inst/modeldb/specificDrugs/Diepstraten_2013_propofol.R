Diepstraten_2013_propofol <- function() {
  description <- "Three-compartment intravenous population PK model for propofol in morbidly obese and nonobese adults, adolescents, and children (Diepstraten 2013 meta-analysis of five previously published studies; N = 94 patients, TBW 37-184 kg, age 9-79 years). Final model E in Table 3: total body weight scales clearance allometrically with an estimated exponent and scales the slow inter-compartmental clearance Q3 linearly; age modifies clearance via a bilinear function centered at 41 years with separate slopes below and above the breakpoint. Inter-individual variability on CL, V1, V3, and Q3 (log-normal) and proportional intra-individual error on log-transformed concentrations."
  reference <- paste(
    "Diepstraten J, Chidambaran V, Sadhasivam S, Blusse van Oud-Alblas HJ,",
    "Inge T, van Ramshorst B, van Dongen EPA, Vinks AA, Knibbe CAJ. (2013).",
    "An integrated population pharmacokinetic meta-analysis of propofol in",
    "morbidly obese and nonobese adults, adolescents, and children.",
    "CPT: Pharmacometrics & Systems Pharmacology 2:e73.",
    "doi:10.1038/psp.2013.47.",
    sep = " "
  )
  vignette <- "Diepstraten_2013_propofol"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as both an allometric scaler on CL with estimated exponent z = 0.77 (reference 70 kg) and a linear scaler on Q3 with exponent fixed at 1 (reference 70 kg) per Diepstraten 2013 Eq. 1 and Table 2 footnote c. Source range 37-184 kg pooled across the five studies (Table 1).",
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Subject age in years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters clearance through a bilinear factor centered at the cohort-median age of 41 years (Diepstraten 2013 Eq. 1 / Eq. 5). Slope below 41 y is b = 0.0103 per year (positive); slope above 41 y is c = -0.00539 per year (negative). At age = 41 y the factor equals 1. Source range 9-79 years.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 94L,
    n_studies      = 5L,
    age_range      = "9-79 years (pooled; cohort median 41 years)",
    weight_range   = "37-184 kg (total body weight; pooled mean 94 kg, SD 35)",
    sex_female_pct = 68,
    disease_state  = "Morbidly obese and nonobese adults, adolescents, and children. Adult cohorts: 20 morbidly obese patients scheduled for bariatric surgery and 40 nonobese patients (24 elective general-surgery patients receiving an induction bolus and 20 ICU patients receiving 2-5 days of continuous sedation). Pediatric / adolescent cohorts: 20 morbidly obese adolescents scheduled for bariatric surgery and 14 nonobese adolescents undergoing scoliosis surgery.",
    dose_range     = "Adults: induction bolus 200 or 350 mg or 2.5 mg/kg, followed by maintenance infusion (initially 10 mg/kg/h titrated to Bispectral Index 40-60) or 2-5 days of ICU sedation by Ramsay scale. Adolescents/children: induction bolus (4 mg/kg in nonobese, dosing-weight calculation per Servin in obese) followed by maintenance infusion 2-10 mg/kg/h.",
    regions        = "Pooled meta-analysis spanning Netherlands and USA (St. Antonius Hospital Nieuwegein; Erasmus Medical Centre Rotterdam; Cincinnati Children's Hospital Medical Center; Leiden / Amsterdam Center for Drug Research).",
    notes          = "1,652 propofol whole-blood concentration measurements pooled across the five source studies (references 6, 7, 30, 31, 32 in Diepstraten 2013). Demographics in Table 1; gender 30 M / 64 F across all cohorts. Concentrations were log-transformed prior to NONMEM VI ADVAN11 / TRANS4 fitting (Methods, Pharmacokinetic model)."
  )

  ini({
    # Final model E (Diepstraten 2013 Table 3 'Final model' column). The
    # paper reports CL and inter-compartmental clearances in L/min and the
    # volumes in L; time units are minutes and predicted concentration is
    # mg/L (= ug/mL). Population mean CL is for the reference subject of
    # 70 kg total body weight and 41 years of age (Eq. 1).
    lcl  <- log(2.34)  ; label("Clearance for a 70 kg, 41 y reference subject (L/min)")   # Diepstraten 2013 Table 3 final model: CL70 kg, 41 y = 2.34 (CV% 4.3)
    lvc  <- log(3.17)  ; label("Central volume of distribution V1 (L)")                   # Diepstraten 2013 Table 3 final model: V1 = 3.17 (CV% 11.3)
    lq   <- log(1.60)  ; label("Inter-compartmental clearance Q2 (L/min)")                # Diepstraten 2013 Table 3 final model: Q2 = 1.60 (CV% 11.7)
    lvp  <- log(5.89)  ; label("Peripheral volume of distribution V2 (L)")                # Diepstraten 2013 Table 3 final model: V2 = 5.89 (CV% 15.0)
    lq2  <- log(1.50)  ; label("Inter-compartmental clearance Q3 at 70 kg (L/min)")       # Diepstraten 2013 Table 3 final model: Q3,70 kg = 1.50 (CV% 6.2)
    lvp2 <- log(116)   ; label("Peripheral volume of distribution V3 (L)")                # Diepstraten 2013 Table 3 final model: V3 = 116 (CV% 7.5)

    # Allometric exponent of TBW on CL (estimated, not fixed at 0.75;
    # Diepstraten 2013 Discussion compares this 0.77 against the
    # theoretical allometric 0.75 and the previously published 0.72 / 0.80).
    e_wt_cl  <- 0.77      ; label("Allometric exponent of TBW/70 on CL (unitless)")       # Diepstraten 2013 Table 3 final model: z = 0.77 (CV% 6.9)

    # Bilinear age slopes on CL with breakpoint at the cohort median age of
    # 41 years (Diepstraten 2013 Eq. 1 and Methods Eq. 5).
    e_age_le41 <-  0.0103   ; label("Slope of (AGE - 41) on age factor for AGE <= 41 y (per year)")  # Diepstraten 2013 Table 3 final model: b = 0.0103 (CV% 13.5)
    e_age_gt41 <- -0.00539  ; label("Slope of (AGE - 41) on age factor for AGE >  41 y (per year)")  # Diepstraten 2013 Table 3 final model: c = -0.00539 (CV% -33.8)

    # Inter-individual variability (log-normal eta on log-scale parameters).
    # Diepstraten 2013 Table 3 final model reports IIV as CV%; the internal
    # variance is omega^2 = log(CV^2 + 1). IIV was retained on CL, V1, V3
    # and Q3 only; no IIV on V2 or Q2 (Table 3 entries blank for those).
    etalcl  ~ 0.030171  # Diepstraten 2013 Table 3 final model: IIV CL 17.5% (CV% 13.9), variance = log(0.175^2 + 1)
    etalvc  ~ 0.228065  # Diepstraten 2013 Table 3 final model: IIV V1 50.6% (CV% 41.3), variance = log(0.506^2 + 1)
    etalvp2 ~ 0.123131  # Diepstraten 2013 Table 3 final model: IIV V3 36.2% (CV% 34.9), variance = log(0.362^2 + 1)
    etalq2  ~ 0.151250  # Diepstraten 2013 Table 3 final model: IIV Q3 40.4% (CV% 37.5), variance = log(0.404^2 + 1)

    # Residual error. Diepstraten 2013 fit on log-transformed concentrations:
    # Yij = log(c_pred,ij) + eps_ij with eps ~ N(0, sigma^2) (Methods Eq. 3).
    # For small errors this is equivalent to a proportional error on the
    # linear scale; carried over here as propSd.
    propSd  <- 0.243   ; label("Proportional residual error (fraction)")                  # Diepstraten 2013 Table 3 final model: proportional intra-individual error 24.3% (CV% 10.3)
  })

  model({
    # Bilinear age factor on CL centered at 41 y (cohort median; Eq. 1).
    #   Fage = 1 + b * (AGE - 41)   for AGE <= 41
    #   Fage = 1 + c * (AGE - 41)   for AGE >  41
    # Implemented with min(0, AGE-41) / max(0, AGE-41) so the two halves
    # combine into a single smooth expression that returns 1 at AGE = 41.
    fage <- 1 + e_age_le41 * min(0, AGE - 41) + e_age_gt41 * max(0, AGE - 41)

    # Individual PK parameters. Allometric WT/70 on CL with estimated
    # exponent and bilinear age on CL; linear WT/70 on Q3 with exponent
    # fixed at 1 (Diepstraten 2013 Table 2 / Eq. 1). No covariates on V1,
    # V2, V3 or Q2 in the final model.
    cl  <- exp(lcl  + etalcl)  * (WT / 70) ^ e_wt_cl * fage
    vc  <- exp(lvc  + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2  + etalq2)  * (WT / 70)
    vp2 <- exp(lvp2 + etalvp2)

    # Three-compartment IV PK with first-order elimination from the central
    # compartment (NONMEM ADVAN11 / TRANS4 mass-balance ODEs; Methods
    # 'Pharmacokinetic model'). Dose lands in `central` via the cmt column
    # of the user data set.
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # Propofol whole-blood concentration in the central compartment.
    # Dose units mg, Vc units L -> Cc units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
