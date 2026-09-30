Li_2021_tedizolid <- function() {
  description <- paste(
    "Two-compartment population PK model for tedizolid (the active moiety of",
    "the prodrug tedizolid phosphate) in adults, adolescents and children",
    "pooled from 16 trials (Li 2021), with linear elimination and sigmoidal",
    "oral absorption: a zero-order release of the oral dose into a depot",
    "followed by first-order absorption into the central compartment, with an",
    "absorption lag time and absolute bioavailability. Intravenous doses are",
    "infused into the central compartment over a modelled, fixed duration.",
    "Body weight (power model, reference 77.3 kg) acts on CL, Q, Vc and Vp;",
    "acute bacterial skin and skin-structure infection (ABSSSI) raises CL and",
    "Vc relative to healthy volunteers; diabetes lowers Vc. Residual error is",
    "exponential (log scale), with fold-increases in the SD for the patient",
    "trials (phase 2 study 104 and the phase 3 trials) and for records after",
    "oral dosing. Doses are in mg of tedizolid free-base equivalent",
    "(200 mg tedizolid phosphate = 164.5 mg tedizolid).",
    sep = " "
  )
  reference <- paste(
    "Li D, Sabato PE, Guiastrennec B, Ouerdani A, Feng H-P, Duval V,",
    "De Anda CS, Sears PS, Chou MZ, Hardalo C, Broyde N, Rizk ML (2021).",
    "Population pharmacokinetics, exposure-response, and probability of",
    "target attainment analyses for tedizolid in adolescent patients with",
    "acute bacterial skin and skin structure infections.",
    "Antimicrobial Agents and Chemotherapy 65(12):e00895-21.",
    "doi:10.1128/AAC.00895-21.",
    "Parameter estimates are from Table 1 (final model, SAEM); the covariate",
    "equations are from the Table 1 footnote b. The absorption structure is",
    "described in Methods ('PopPK model development') and in the backbone",
    "model it updates, Flanagan S et al. (2014) Antimicrob Agents Chemother",
    "58:6462-6470, doi:10.1128/AAC.03423-14.",
    sep = " "
  )
  vignette <- "Li_2021_tedizolid"
  units <- list(
    time = "h",
    dosing = "mg (tedizolid free-base equivalent; 200 mg tedizolid phosphate = 164.5 mg tedizolid)",
    concentration = "ug/mL (mg/L)"
  )

  compartmentData <- list(
    depot = list(analyte = "tedizolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tedizolid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tedizolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model on CL, Q (sharing the CL exponent), Vc and Vp, normalized",
        "to 77.3 kg (Li 2021 Table 1 footnote b equations and footnote g:",
        "'Volumes and CL are reported for a typical individual of 77.3 kg').",
        "Analysis range 12.6-226 kg, median 76.0 kg (supplement Table A2).",
        sep = " "
      ),
      source_name = "WT"
    ),
    DIS_CSSSI = list(
      description = "Acute bacterial skin and skin-structure infection (ABSSSI) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no ABSSSI (healthy volunteer)",
      notes = paste(
        "Source flag INFEC, 'taking the value of 1 in case of infection ... and",
        "0 if otherwise' (Li 2021 Table 1 footnote b). Linear fractional effect",
        "on CL (+22.0%) and Vc (+9.87%). The reference population for the",
        "disease effect was changed from ABSSSI patients to healthy volunteers",
        "in this update (Methods, 'PopPK model development'), and Results",
        "('PopPK analysis') names the covariate as ABSSSI. How the 41",
        "hospitalized children and adolescents with suspected (not confirmed)",
        "Gram-positive infection in trials PN013 and PN026 were coded is not",
        "reported.",
        sep = " "
      ),
      source_name = "INFEC"
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = nondiabetic",
      notes = paste(
        "Source flag DIAB (Li 2021 Table 1 footnote b). Linear fractional",
        "effect on Vc only (-14.3%). 101 of 1,312 participants (7.7%) were",
        "diabetic; diabetic status was missing (coded -99) for 232 participants",
        "(supplement Table A1); how the missing value entered the covariate",
        "equation is not reported. No adolescent had diabetes.",
        sep = " "
      ),
      source_name = "DIAB"
    ),
    STUDY_PHASE3 = list(
      description = "Patient-trial residual-error stratum indicator (phase 2 study 104 and the phase 3 trials)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-phase 3 trials (phase 1 studies in healthy volunteers and hospitalized children / adolescents)",
      notes = paste(
        "1 = record from 'study 104 and phase 3 trials', 0 = record from the",
        "'non-phase 3 trials' (Li 2021 Table 1 residual-variability rows).",
        "Study 104 is the phase 2 ABSSSI trial, which the source pools into the",
        "phase 3 stratum, so set STUDY_PHASE3 = 1 for it as well. Selects the",
        "residual-error magnitude only (4.92-fold SD); the typical-value",
        "prediction is unaffected.",
        sep = " "
      ),
      source_name = "not published"
    ),
    ROUTE_ORAL = list(
      description = "Oral-administration record indicator (reference = intravenous infusion)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = intravenous infusion",
      notes = paste(
        "1 = observation record following oral tedizolid phosphate ('RV for",
        "oral data', Li 2021 Table 1), 0 = intravenous. Selects the",
        "residual-error magnitude only (2.01-fold SD); the structural route",
        "difference is carried by the dose record's target compartment",
        "(depot for oral, central for the intravenous infusion).",
        sep = " "
      ),
      source_name = "not published"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "No additional covariate screening was conducted in this update; 'No clear trends or specific impact of age or creatinine clearance on exposure were predicted' (Li 2021 Results, 'PopPK analysis'; Discussion).",
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "not reported",
      type = "continuous",
      notes = "Not retained; 'No clear trends or specific impact of age or creatinine clearance on exposure were predicted' (Li 2021 Results, 'PopPK analysis').",
      source_name = "not published"
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "A covariate on CL and Vc in the earlier Flanagan 2014 model; removed during the interim model updates and replaced by total body weight (Li 2021 Methods, 'Data selection' and 'PopPK model development').",
      source_name = "IBW"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "A covariate on CL in the earlier Flanagan 2014 model; removed during the interim model updates (Li 2021 Methods, 'PopPK model development').",
      source_name = "BILI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1312,
    n_studies = 16,
    n_observations = 9756,
    age_range = "3-94 years (median 40.0; 132 participants < 18 years)",
    weight_range = "12.6-226 kg (median 76.0, mean 77.6)",
    sex_female_pct = 32.6,
    race_ethnicity = c(White = 70.9, Asian = 15.2, Black = 12.1, Other = 1.8),
    disease_state = paste(
      "945 adults with ABSSSI, 223 healthy participants, 41 hospitalized",
      "children and adolescents with suspected Gram-positive infection, and",
      "103 adolescents with ABSSSI (91 from the phase 3 trial PN012)."
    ),
    dose_range = "200 mg tedizolid phosphate once daily, oral or intravenous, in the ABSSSI trials; the pooled phase 1 studies are not itemized by dose",
    regions = "Not tabulated; five trials enrolled only Asian participants (PN005, PN006, BAY-16101, BAY-16102, BAY-16411; supplement Table A1), and the interim updates added Japanese and Chinese patients (Methods)",
    notes = paste(
      "Update of the Flanagan 2014 tedizolid model, adding the phase 3",
      "adolescent ABSSSI trial PN012 (12 to < 18 years) and the ongoing phase 1",
      "trial PN013 (2 to < 12 years). 5,146 oral and 4,647 intravenous PK",
      "samples; lower limit of quantification 5 ng/mL; 37 outlier",
      "observations (CWRES > 6) excluded. Diabetes present in 7.7%",
      "(Li 2021 Results 'Participants'; supplement Tables A1-A3)."
    )
  )

  ini({
    # Structural parameters -- Li 2021 Table 1, final model (SAEM), typical
    # individual of 77.3 kg (footnote g).
    ld2 <- fixed(log(0.810))
    label("Log of the modelled intravenous infusion duration into central (h)")
    # Table 1 'Infusion time, fixed (h)' = 0.810 FIXED

    lfdepot <- log(0.857)
    label("Log of the oral bioavailability (fraction)")
    # Table 1 'F1' = 0.857 (RSE 0.965%)

    ld1 <- log(0.175)
    label("Log of the zero-order release duration of the oral dose into depot (h)")
    # Table 1 'Zero-order duration (h)' = 0.175 (RSE 28.3%)

    lka <- log(1.47)
    label("Log of the first-order absorption rate constant from depot (1/h)")
    # Table 1 'Ka (h-1)' = 1.47 (RSE 9.25%)

    ltlag <- log(0.226)
    label("Log of the oral absorption lag time (h)")
    # Table 1 'Lag time (h)' = 0.226 (RSE 0.0376%)

    lcl <- log(5.39)
    label("Log of clearance for a 77.3 kg healthy participant (L/h)")
    # Table 1 'CL (liter/h)' = 5.39 (RSE 6.98%)

    e_wt_cl <- 0.408
    label("Power exponent of body weight on CL and Q (unitless)")
    # Table 1 CL 'wt (power model)' = 0.408 (RSE 10.2%); Q 'wt (power model)' = 'Same as for CL'

    e_csssi_cl <- 0.220
    label("Fractional change in CL with ABSSSI (unitless)")
    # Table 1 CL 'Infection (%) (linear model)' = 22.0 (RSE 39.3%)

    lvc <- log(58.5)
    label("Log of central volume for a 77.3 kg healthy nondiabetic participant (L)")
    # Table 1 'Vc (liter)' = 58.5 (RSE 3.36%)

    e_wt_vc <- 0.903
    label("Power exponent of body weight on Vc (unitless)")
    # Table 1 Vc 'wt (power model)' = 0.903 (RSE 3.53%)

    e_csssi_vc <- 0.0987
    label("Fractional change in Vc with ABSSSI (unitless)")
    # Table 1 Vc 'Infection (%) (linear model)' = 9.87 (RSE 37.4%)

    e_dis_diab_vc <- -0.143
    label("Fractional change in Vc with diabetes (unitless)")
    # Table 1 Vc 'Diabetes (%) (linear model)' = -14.3 (RSE 22.2%); the PDF minus sign renders as '2' (bootstrap -14.0, 95% CI -18.6 to -9.49)

    lq <- log(1.43)
    label("Log of intercompartmental clearance for a 77.3 kg participant (L/h)")
    # Table 1 'Q(L/h)' = 1.43 (RSE 4.09%)

    lvp <- log(15.6)
    label("Log of peripheral volume for a 77.3 kg participant (L)")
    # Table 1 'Vp (L)' = 15.6 (RSE 2.41%)

    e_wt_vp <- 0.678
    label("Power exponent of body weight on Vp (unitless)")
    # Table 1 Vp 'wt (power model)' = 0.678 (RSE 6.87%)

    # Interindividual variability -- Table 1 'IIV, %CV', converted to log-scale
    # variances as omega^2 = log(CV^2 + 1).
    etald2 ~ fixed(0.00680) # Table 1 'Infusion time' IIV 8.26% CV, not estimated; log(0.0826^2 + 1) = 0.00680
    etald1 ~ 2.036 # Table 1 'Zero-order duration' IIV 258% CV (RSE 7.39%); log(2.58^2 + 1) = 2.036
    etalka ~ 0.4656 # Table 1 'Ka' IIV 77.0% CV (RSE 8.36%); log(0.770^2 + 1) = 0.4656
    etaltlag ~ 0.6931 # Table 1 'Lag time' IIV 100% CV (RSE 5.47%); log(1.00^2 + 1) = 0.6931
    etalcl + etalvc ~ c(0.09865, 0.04840, 0.06157) # Table 1 CL IIV 32.2% CV -> 0.09865; Vc IIV 25.2% CV -> 0.06157; 'Correlation CL-Vc' 62.1% -> cov = 0.621 * sqrt(0.09865 * 0.06157) = 0.04840
    etalvp ~ 0.02466 # Table 1 'Vp' IIV 15.8% CV (RSE 8.54%); log(0.158^2 + 1) = 0.02466

    # Residual error -- exponential (log-scale additive) with fold-increases
    # of the SD ('Fold increase on the square root scale', footnote h).
    expSd <- 0.123
    label("Log-scale residual SD, non-phase 3 trials (fraction)")
    # Table 1 'RV for non-phase 3 trials (%)' = 12.3 (RSE 1.36%)

    e_study_phase3_expsd <- 4.92
    label("Fold increase of the residual SD in study 104 and the phase 3 trials (unitless)")
    # Table 1 'RV for study 104 and phase 3 trials (fold)' = 4.92 (RSE 4.63%)

    e_route_oral_expsd <- 2.01
    label("Fold increase of the residual SD for records after oral dosing (unitless)")
    # Table 1 'RV for oral data (fold)' = 2.01 (RSE 5.22%)
  })

  model({
    # Covariate equations -- Li 2021 Table 1 footnote b:
    #   CLi = TVCL * (WT/77.3)^thetaWT-CL * (1 + thetaInfec-CL * INFEC) * exp(etaCLi)
    #   Vci = TVVc * (WT/77.3)^thetaWT-Vc * (1 + thetaInfec-Vc * INFEC)
    #         * (1 + thetaDiab-Vc * DIAB) * exp(etaVci)
    #   Vpi = TVVp * (WT/77.3)^thetaWT-Vp * exp(etaVpi)
    # Q carries the CL weight exponent ('Same as for CL') and no IIV.
    wt_ratio <- WT / 77.3
    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl * (1 + e_csssi_cl * DIS_CSSSI)
    vc <- exp(lvc + etalvc) *
      wt_ratio^e_wt_vc *
      (1 + e_csssi_vc * DIS_CSSSI) *
      (1 + e_dis_diab_vc * DIS_DIAB)
    q <- exp(lq) * wt_ratio^e_wt_cl
    vp <- exp(lvp + etalvp) * wt_ratio^e_wt_vp

    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)
    d2 <- exp(ld2 + etald2)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Sigmoidal absorption: the oral dose is released into depot as a
    # zero-order input of duration d1 (after the lag time) and absorbed into
    # central by a first-order process (Li 2021 Methods; Flanagan 2014).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Oral dose records need rate = -2 so the modelled d1 is used; intravenous
    # dose records into central need rate = -2 so the modelled infusion time
    # d2 is used.
    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1
    dur(central) <- d2

    # Dose in mg and volume in L give mg/L = ug/mL.
    Cc <- central / vc

    expSdCc <- expSd *
      e_study_phase3_expsd^STUDY_PHASE3 *
      e_route_oral_expsd^ROUTE_ORAL
    Cc ~ lnorm(expSdCc)
  })
}
