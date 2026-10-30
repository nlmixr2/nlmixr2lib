Chamzas_2024_tacrolimus <- function() {
  description <- "Two-compartment population PK model with Erlang-distributed transit absorption (3 transit compartments sharing one rate constant ktr) for twice-daily oral tacrolimus in self-identified Hispanic adult renal transplant recipients during the first 30 days post-transplant (Chamzas 2024). The structural parameters (ktr, Vc/F, Vp/F, Q/F) and their between-subject variabilities were estimated from digitised rich profiles of a separate bioequivalence trial and fixed; only apparent clearance CL/F and its variability were estimated from the sparse real-world trough data. CL/F, Q/F scale allometrically with body weight (exponent 0.75) and Vc/F, Vp/F linearly (exponent 1), all referenced to 70 kg; CL/F additionally decreases with time-varying hematocrit as a power function referenced to 30 %. Additive residual error on whole-blood concentrations."
  reference <- paste(
    "Chamzas A, Tellez E, SyBing A, Gobburu JVS, Gopalakrishnan M.",
    "Optimizing tacrolimus dosing in Hispanic renal transplant patients:",
    "insights from real-world data.",
    "Front Pharmacol. 2024;15:1443988.",
    "doi:10.3389/fphar.2024.1443988.",
    sep = " "
  )
  vignette <- "Chamzas_2024_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg (Equations 5-8): exponent 0.75 on CL/F and Q/F, exponent 1 on Vc/F and Vp/F. Cohort mean (SD) 77.3 (20.1) kg (Table 1).",
      source_name = "Wt"
    ),
    HCT = list(
      description = "Haematocrit, packed red blood cell volume fraction times 100",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying over the first 30 days post-transplant (Abstract: 'time-varying hematocrit on apparent clearance'). Enters CL/F as (HCT/30)^-0.83 (Equation 5 and Table 2 footnote: CL = tvcl * (HCT/HCT_median)^-HCTcov * ...; HCTcov = 0.83). The reference 30 % is the typical patient named in the Abstract ('a typical patient weighing 70 kg and at 30% hematocrit') and in Supplementary Table S2, and equals the cohort mean (SD) 30.0 (4.9) % in Table 1. Higher hematocrit lowers whole-blood apparent clearance because tacrolimus partitions extensively into erythrocytes.",
      source_name = "HT / HCT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 121L,
    n_studies = 1L,
    n_observations = 2215L,
    age_mean_sd = "45.7 (12.9) years",
    weight_mean_sd = "77.3 (20.1) kg",
    bmi_mean_sd = "27.6 (5.2) kg/m^2",
    sex_female_pct = 39,
    race_ethnicity = "100% self-identified Hispanic.",
    disease_state = "Adult renal (kidney) allograft recipients during the first 30 days of tacrolimus therapy after transplantation; 87% first transplant, 71% deceased donor. Baseline serum creatinine mean 4.7 mg/dL, creatinine clearance mean 29.9 mL/min, albumin 3.6 g/dL, haematocrit 30.0% (Table 1).",
    dose_range = "Oral tacrolimus twice daily, initial dose 0.02-0.14 mg/kg/day (median 0.06 mg/kg/day) titrated to a morning trough target of 10-12 ng/mL; maximum single dose 13 mg.",
    co_medication = "Thymoglobulin induction, methylprednisolone tapered to prednisone, mycophenolate mofetil 1 g twice daily (Methods 2.2).",
    regions = "United States (St. Joseph Hospital, Orange, California), transplants 2008-2020.",
    notes = "Retrospective electronic-medical-record data, sparse morning troughs only (mean 18 per patient), whole blood by LC-MS/MS; 1.7% of samples below the LLOQ were set to missing. Structural parameters were taken from a fit to 35 stable renal transplant recipients in the Alloway 2017 bioequivalence trial (15 samples over a 12-h steady-state interval) and fixed. Fitted in Pumas with FOCE."
  )

  ini({
    # Structural parameters fixed from the rich-data (Alloway 2017
    # bioequivalence) fit; values from Table 2 and Results 3.1.
    lktr <- fixed(log(4.27)); label("Transit rate constant ktr, shared by all transit steps (1/h)") # Table 2 Ktr = 4.27 (fixed); Results 3.1 '4.27 h-1'
    lcl <- log(13.65); label("Apparent clearance CL/F at 70 kg and 30% hematocrit (L/h)") # Table 2 CL/F = 13.65 L/hr (RSE 4.32%)
    lvc <- fixed(log(171)); label("Apparent central volume Vc/F at 70 kg (L)") # Table 2 Vc/F = 171 (fixed)
    lvp <- fixed(log(325.6)); label("Apparent peripheral volume Vp/F at 70 kg (L)") # Results 3.1 '325.6 L'; Table 2 prints 325 (fixed)
    lq <- fixed(log(62.46)); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 2 Q/F = 62.46 (fixed); Results 3.1 '62.5 L/h'

    # Covariate effects (Equations 5-8)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F and Q/F (unitless)") # Equations 5 and 8: (Wt/70)^0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT/70 on Vc/F and Vp/F (unitless)") # Equations 6 and 7: (Wt/70)
    e_hct_cl <- -0.83; label("Power exponent of HCT/30 on CL/F (unitless)") # Table 2 HCTcov = 0.83 (RSE 8.46%), entering Equation 5 with a negative sign

    # Between-subject variability. Table 2 reports %BSV; the variances below
    # are log(1 + (%BSV/100)^2), which gives round values (0.47, 0.88, ...)
    # consistent with %BSV being the exponential-scale CV of round omega^2.
    etalcl ~ 0.4700 # Table 2 CL/F BSV 77.46%: log(1 + 0.7746^2) = 0.4700
    etalvc ~ fixed(0.8800) # Table 2 Vc/F BSV 118.78%, from the rich-data fit: log(1 + 1.1878^2) = 0.8800
    etalvp ~ fixed(0.7140) # Table 2 Vp/F BSV 102.09%, from the rich-data fit: log(1 + 1.0209^2) = 0.7140
    etalq ~ fixed(0.7550) # Table 2 Q/F BSV 106.19%, from the rich-data fit: log(1 + 1.0619^2) = 0.7550
    etalktr ~ fixed(0.5580) # Table 2 Ktr BSV 86.44%, from the rich-data fit: log(1 + 0.8644^2) = 0.5580

    addSd <- 2.62; label("Additive residual error (ng/mL)") # Table 2 RUV additive error = 2.62 ng/mL (RSE 12.85%)
  })

  model({
    ktr <- exp(lktr + etalktr)
    cl <- exp(lcl + etalcl) * (HCT / 30)^e_hct_cl * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Three transit compartments (depot, transit1, transit2) connected by
    # the single estimated ktr; the last transit feeds central at ktr.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose mg, volume L -> mg/L; x1000 for whole-blood ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
