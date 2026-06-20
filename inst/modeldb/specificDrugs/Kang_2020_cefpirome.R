Kang_2020_cefpirome <- function() {
  description <- "Two-compartment IV-bolus population PK model for cefpirome in critically ill adults on venoarterial extracorporeal membrane oxygenation (VA-ECMO) (Kang 2020). Final-model covariates: power-form serum creatinine on CL (reference 1.6 mg/dL), and a binary ECMO-active treatment-status indicator on CL (1.41-fold higher when ECMO-ON) and V1 (4.22-fold higher when ECMO-ON)."
  reference <- paste(
    "Kang S, Jang JY, Hahn J, Kim D, Lee JY, Min KL, Yang S, Wi J, Chang MJ.",
    "Dose optimization of cefpirome based on population pharmacokinetics and",
    "target attainment during extracorporeal membrane oxygenation.",
    "Antimicrob Agents Chemother. 2020;64(5):e00249-20.",
    "doi:10.1128/AAC.00249-20.",
    sep = " "
  )
  vignette <- "Kang_2020_cefpirome"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine concentration (time-varying during sampling)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column SCr. Reference value 1.6 mg/dL (Kang 2020 Results, final-model equation 'CL = 5.71 * 0.487^(SCr/1.6) * 1.41^ECMO'; the 1.6 mg/dL centering point appears in the model expression itself, not in Table 1 demographics which report a cohort median of 1.58 mg/dL during ECMO and 1.83 mg/dL after weaning). Treated as time-varying -- Kang 2020 Methods (Covariate model development) states 'All data were recorded during sampling and tested as time-varying covariates'. Effect on CL is a power form: CL = 5.71 * 0.487^(SCr/1.6); SCr appears in the exponent. Cohort range during sampling 0.40-3.41 mg/dL (Table 1).",
      source_name        = "SCr"
    ),
    ECMO_STATUS = list(
      description        = "Binary indicator: 1 = patient currently on VA-ECMO support, 0 = patient off VA-ECMO (post-weaning) (time-varying within subject)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column ECMO. ECMO modality: venoarterial (VA-ECMO) only; the Kang 2020 cohort used a Capiox SP-101 centrifugal pump with a Capiox EBS X-coated conduit (Terumo) and Sechrist air-oxygen mixer (Methods, Extracorporeal membrane oxygenation system). Time-varying within subject: PK samples were drawn during ECMO (ECMO-ON, 14/15 patients) and after successful weaning (ECMO-OFF, 8/15 patients) -- the per-record indicator switches when the patient is decannulated. Circuit changes within an ECMO run were not modelled. Effect on CL is a power-form multiplier: CL = 5.71 * 0.487^(CREAT/1.6) * 1.41^ECMO_STATUS (i.e., 1.41-fold higher CL when ECMO is on). Effect on V1 is V1 = 2.74 * 4.22^ECMO_STATUS (i.e., 4.22-fold higher central volume when ECMO is on). Mechanistic discussion (Kang 2020 Discussion): the V1 expansion reflects the extra circulating volume contributed by the ECMO circuit; the CL increase reflects hyperdynamic vasodilatory state and altered drug-metabolizing-pathway activity in critical illness.",
      source_name        = "ECMO"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 15L,
    n_studies       = 1L,
    age_range       = "27-82 years",
    age_median      = "63 years (IQR 51.5-70.5)",
    weight_range    = "49.2-98.3 kg on ECMO; 48.4-84.9 kg off ECMO",
    weight_median   = "70.3 kg on ECMO (IQR 60.8-76.0); 62.25 kg off ECMO (IQR 57.3-70.2)",
    sex_female_pct  = 20,
    race_ethnicity  = "Not reported (single-centre South Korean cohort at Severance Hospital, Seoul; presumed predominantly Korean)",
    disease_state   = "Critically ill adults receiving venoarterial ECMO (VA-ECMO) for profound cardiogenic shock; indications were acute myocardial infarction (n = 11), arrhythmogenic right ventricular dysplasia (n = 1), myocarditis (n = 1), pulmonary thromboembolism (n = 1), and one additional AMI variant. Cefpirome administered as per hospital infection-prophylaxis protocol.",
    dose_range      = "Cefpirome 2 g IV bolus q12h per hospital protocol; 50% dose reduction (1 g IV bolus q12h) for patients with eGFR < 50 mL/min/1.73 m^2 (MDRD); occasional dose adjustments not specified per patient.",
    regions         = "South Korea (single-centre prospective cohort at Severance Hospital, Yonsei University College of Medicine, Seoul; cardiac intensive care unit; January 2018 - January 2019; IRB 4-2014-0919; ClinicalTrials.gov NCT02581280).",
    renal_function  = "Serum creatinine median 1.58 mg/dL during ECMO (IQR 0.99-2.36) and 1.84 mg/dL after weaning (IQR 1.48-2.27). 5/15 patients received continuous renal replacement therapy (CRRT) during ECMO; CRRT was screened as a covariate but dropped during stepwise selection (Kang 2020 Discussion).",
    apache2_score   = "Median 32 (IQR 29-36)",
    ecmo_duration   = "Median 166.1 h (IQR 124.2-254.0; range 34.6-720.2)",
    notes           = "Prospective cohort study (STROBE-conformant). 152 plasma samples (no BLQ) drawn at predose and 0.5-1, 2-3, 4-6, 8-10, and 12 h after cefpirome administration during ECMO (ECMO-ON, 14 patients sampled) and after ECMO weaning (ECMO-OFF, 8 patients sampled). LC-MS/MS quantification (validated 1.0-64.0 mg/L; LLOQ 1.0 mg/L). NONMEM 7.4.1 (FOCE-I) with Pirana 2.9.7 and Xpose4 v4.6.1 in R 3.5.3. Model validation by sampling-importance-resampling (SIR; 5000/1000, 5 iterations) and prediction-corrected VPC (n = 5000)."
  )

  ini({
    # Structural parameters. The 'typical-subject' anchor in Kang 2020 is the
    # ECMO-OFF arm with SCr = 0; physiological reference is the ECMO-OFF arm
    # with SCr = 1.6 mg/dL (the in-equation centering point that gives CL = 2.78
    # L/h and V1 = 2.74 L). Kang 2020 Table 2 final-model column.
    lcl <- log(5.71); label("Clearance at SCr = 0, ECMO-OFF (L/h)")                 # Kang 2020 Table 2: CL = 5.71 L/h (final-model fixed effect theta_CL)
    lvc <- log(2.74); label("Central volume of distribution V1 at ECMO-OFF (L)")     # Kang 2020 Table 2: V1 = 2.74 L (final-model fixed effect theta_V1)
    lvp <- log(16.7); label("Peripheral volume of distribution V2 (L)")              # Kang 2020 Table 2: V2 = 16.7 L
    lq  <- log(9.43); label("Intercompartmental clearance Q (L/h)")                  # Kang 2020 Table 2: Q  = 9.43 L/h

    # Covariate effects.
    #   CL  = exp(lcl + etalcl) * theta_screatcl^(CREAT / 1.6)
    #                            * theta_ecmocl^ECMO_STATUS
    #   V1  = exp(lvc + etalvc) * theta_ecmovc^ECMO_STATUS
    # The exponents are themselves estimated parameters (not fixed), so they are
    # encoded as bare `ini()` entries and read from Table 2.
    theta_screatcl <- 0.487; label("Multiplicative factor for (CREAT / 1.6 mg/dL) exponent on CL")     # Kang 2020 Table 2: theta_SCr/1.6 on CL = 0.487
    theta_ecmocl   <- 1.41;  label("Multiplicative factor for ECMO_STATUS exponent on CL")             # Kang 2020 Table 2: theta_ECMO on CL  = 1.41
    theta_ecmovc   <- 4.22;  label("Multiplicative factor for ECMO_STATUS exponent on V1")             # Kang 2020 Table 2: theta_ECMO on V1  = 4.22

    # Inter-individual variability (Kang 2020 Table 2 final-model %CV);
    # omega^2 = log(CV^2 + 1) for log-normal etas (exponential variance model
    # per Methods, Base model development).
    etalcl ~ 0.11618 # log(0.351^2 + 1); 35.1% CV on CL (Kang 2020 Table 2)
    etalvc ~ 0.13092 # log(0.374^2 + 1); 37.4% CV on V1 (Kang 2020 Table 2)
    etalvp ~ 0.20345 # log(0.475^2 + 1); 47.5% CV on V2 (Kang 2020 Table 2)

    # Proportional residual variability (Kang 2020 Table 2 final-model row
    # 'Proportional residual variability (sigma)'); the paper reports residual
    # variability on a %CV scale, and a proportional-only error model is
    # equivalent in linear DV space to propSd = sigma_prop_fraction.
    propSd <- 0.217; label("Proportional residual error (fraction)")                 # Kang 2020 Table 2: sigma = 21.7% CV (proportional-only)
  })
  model({
    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * theta_screatcl^(CREAT / 1.6) * theta_ecmocl^ECMO_STATUS
    vc <- exp(lvc + etalvc) * theta_ecmovc^ECMO_STATUS
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
