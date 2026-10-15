Wong_2025_efanesoctocog <- function() {
  description <- "One-compartment population PK model for factor VIII (FVIII) activity after intravenous efanesoctocog alfa (BIVV001, a high-sustained FVIII replacement fused to Fc, the D'D3 domain of von Willebrand factor and two XTEN polypeptides) in adults, adolescents and children with severe hemophilia A (Wong 2025). Linear elimination; body weight enters CL and V as power functions centred at 78.3 kg, and Asian race lowers CL by 10.4%. Correlated log-normal IIV on CL and V; combined additive and proportional residual error on FVIII activity measured by the one-stage clotting assay. Sister bleed repeated-time-to-event models from the same paper: modellib('Wong_2025_efanesoctocog_bleed_adult') and modellib('Wong_2025_efanesoctocog_bleed_pediatric')."
  reference <- paste(
    "Wong N, Bhagunde P, Nyberg J, Katragadda S, Demissie M, Willemze A,",
    "Benson C, Macha S. Efanesoctocog Alfa Population Pharmacokinetics and",
    "Repeated Time-To-Event Analysis of Bleeds in Adults, Adolescents, and",
    "Children with Severe Hemophilia A. J Clin Pharmacol. 2025;65(7):860-872.",
    "doi:10.1002/jcph.70008. Open Access under CC BY-NC.",
    "The CL and V equations are in the main text (Results, PopPK Model",
    "Development); the final parameter estimates are in Table S4 of the",
    "supplement (Wiley file JCPH-65-860-s001.docx).",
    sep = " "
  )
  vignette <- "Wong_2025_efanesoctocog"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING in the source analysis ('Time-varying body weight effect was included to allow consideration of body weight changes as children age', Methods). Power effects on CL (exponent 0.677) and V (exponent 0.791), both centred at 78.3 kg, the median baseline body weight of XTEND-1 and the Phase 1/2a studies (Table 1 footnote b). Age was not tested once weight was in the model because of the weight-age correlation.",
      source_name = "Body weight"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator: 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Time-fixed per subject. Multiplicative factor 0.896 on CL, i.e. CL 10.4% lower in Asian patients (Results; Table S4 'Proportional coefficient of Asian race on CL'). 40 of 260 patients (15.4%) were Asian (Table S3).",
      source_name = "Asian"
    )
  )

  covariatesDataExcluded <- list(
    VWF = list(
      description = "Time-varying von Willebrand factor antigen.",
      units = "IU/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained (Results, Model Development). Consistent with the design of efanesoctocog alfa, whose D'D3 domain decouples it from endogenous VWF."
    ),
    HCT = list(
      description = "Time-varying hematocrit.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained (Results, Model Development)."
    ),
    BLOOD_TYPE = list(
      description = "ABO blood type.",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = "Screened but not retained (Results, Model Development)."
    ),
    HCV_POS = list(
      description = "Hepatitis C infection status.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (negative)",
      notes = "Baseline status screened but not retained (Results, Model Development)."
    ),
    HIV_POS = list(
      description = "HIV infection status.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (negative)",
      notes = "Baseline status screened but not retained (Results, Model Development). 24 of 260 patients (9.2%) were HIV-positive (Table S3)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "efanesoctocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 260L,
    n_studies = 5L,
    age_range = "1.4-72 years (median 30)",
    weight_range = "5th-95th percentile 17.5-113 kg (median 72.0)",
    sex_female_pct = 0.4,
    race_ethnicity = c(White = 68.1, Asian = 15.4, Black = 1.9, Other = 3.1, `Not reported` = 11.5),
    disease_state = "Previously treated patients with severe hemophilia A (<1 IU/dL endogenous FVIII activity).",
    dose_range = "Intravenous efanesoctocog alfa: single doses of 25 or 65 IU/kg (EXTEN-A); four once-weekly doses of 50 or 65 IU/kg (Phase 1 repeat-dose study); 50 IU/kg once weekly for 52 weeks (XTEND-1, XTEND-Kids, XTEND-ed); XTEND-1 Arm B received 50 IU/kg on demand for 26 weeks before switching to once-weekly prophylaxis.",
    regions = "multinational",
    notes = "Model-development dataset: 3054 post-dose FVIII activity records from 260 patients (199 adults/adolescents and 61 children <12 years; 259 male, 1 female), about 3% below the 1 IU/dL LLOQ, pooled from EXTEN-A (NCT03205163), the Phase 1 repeat-dose study (EudraCT 2018-001535-51), XTEND-1 (NCT04161495), XTEND-Kids (NCT04759131) and one patient from XTEND-ed (NCT04644575) (Results; Tables S1-S3). A later maximum a posteriori assessment on 4405 records from 277 patients, including de novo Chinese participants from XTEND-ed, found the model adequate without re-estimation. FVIII activity was measured with an aPTT-based one-stage clotting assay (Actin FSL)."
  )

  ini({
    lcl <- log(0.433); label("Clearance at 78.3 kg in a non-Asian patient (dL/h)") # Table S4 'Typical value of CL' = 0.433 dL/h (RSE 1.43%)
    lvc <- log(30.2); label("Central volume of distribution at 78.3 kg (dL)") # Table S4 'Typical value of V' = 30.2 dL (RSE 1.04%)

    e_wt_cl <- 0.677; label("Power exponent of body weight on CL, centred at 78.3 kg (unitless)") # Table S4 'Power of coefficient of body weight on CL' = 0.677 (RSE 3.91%)
    e_wt_vc <- 0.791; label("Power exponent of body weight on V, centred at 78.3 kg (unitless)") # Table S4 'Power of coefficient of body weight on V' = 0.791 (RSE 3.15%)
    e_race_asian_cl <- log(0.896); label("Log fold change in CL for Asian vs non-Asian patients (unitless)") # Table S4 'Proportional coefficient of Asian race on CL' = 0.896 (RSE 2.56%); main-text equation (0.896)^Asian

    # Table S4 IIV block. 0.0354 and 0.0209 are log-scale VARIANCES:
    # sqrt(exp(0.0354) - 1) = 19.0% and sqrt(exp(0.0209) - 1) = 14.5%
    # reproduce the printed 19.0 / 14.6 %CV. The row labelled 'Corr-CL-V'
    # = 0.0148 is the COVARIANCE (OMEGA(2,1)): it implies a correlation of
    # 0.0148 / sqrt(0.0354 * 0.0209) = 0.544, whereas read as a correlation
    # coefficient it would be ~0.015, which is incompatible with its 21.6%
    # RSE and with the base model being built 'with correlation between CL
    # and V'.
    etalcl + etalvc ~ c(
      0.0354,
      0.0148, 0.0209
    ) # Table S4: IIV CL 0.0354 (19.0 %CV); Corr-CL-V 0.0148; IIV V 0.0209 (14.6 %CV)

    propSd <- 0.182; label("Proportional residual error (fraction)") # Table S4 'Proportional term' = 0.182 (RSE 5.70%)
    addSd <- 0.837; label("Additive residual error (IU/dL)") # Table S4 'Additive term (IU/dL)' = 0.837 (RSE 39.4%)
  })
  model({
    # Main-text equations (Results, PopPK Model Development):
    #   CL = 0.433 * (WT/78.3)^0.677 * 0.896^Asian * exp(eta1)
    #   V  = 30.2  * (WT/78.3)^0.791 * exp(eta2)
    cl <- exp(lcl + etalcl + e_race_asian_cl * RACE_ASIAN) * (WT / 78.3)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 78.3)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # FVIII activity; dose in IU over V in dL gives IU/dL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
