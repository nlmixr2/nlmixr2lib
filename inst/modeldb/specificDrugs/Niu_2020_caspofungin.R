Niu_2020_caspofungin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination and",
    "proportional residual error for once-daily 1-hour intravenous",
    "caspofungin infusions (70 mg/m^2 loading dose on day 1, then 50 mg/m^2",
    "daily) in children aged 0.6-14 years undergoing allogeneic",
    "hematopoietic stem cell transplantation (Niu 2020). Clearance scales",
    "with body surface area (power 0.89) and with the natural log of",
    "aspartate aminotransferase (power -0.23 on ln(AST), centred on",
    "ln(AST) = 3.38); volume of distribution is proportional to body",
    "surface area. Both are centred on BSA = 0.79 m^2. Exponential",
    "inter-individual variability on CL and Vd."
  )
  reference <- paste(
    "Niu C-H, Xu H, Gao L-L, Nie Y-M, Xing L-P, Yu L-P, Wu S-L, Wang Y.",
    "Population Pharmacokinetics of Caspofungin and Dosing Optimization in",
    "Children With Allogeneic Hematopoietic Stem Cell Transplantation.",
    "Front Pharmacol. 2020;11:184. doi:10.3389/fphar.2020.00184."
  )
  vignette <- "Niu_2020_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area (Mosteller formula)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "BSA = sqrt(height [cm] * weight [kg] / 3600) (Mosteller; Niu 2020",
        "Methods, 'Population Pharmacokinetic Modeling of Caspofungin').",
        "Power covariate on CL (exponent 0.89) and proportional (exponent 1)",
        "on Vd, both centred on 0.79 m^2 (Results, final-model equations).",
        "Cohort mean 0.80 m^2 (SD 0.27, range 0.38-1.50; Table 1 and",
        "Results 'Study Population'). Doses are BSA-based (mg/m^2), so the",
        "same BSA also sets the administered mg amount."
      ),
      source_name = "BSA"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Entered as the natural log: CL is multiplied by",
        "(ln(AST) / 3.38)^-0.23 (Results, final-model equation), i.e. the",
        "reference is ln(AST) = 3.38, AST = exp(3.38) ~= 29.4 U/L. The",
        "form requires AST > 1 U/L so that ln(AST) is positive. Cohort",
        "median 26 U/L (Table 1); 37 children had AST <= 40 U/L, 9 had",
        "40-120 U/L and 2 had > 120 U/L (Results 'Study Population')."
      ),
      source_name = "AST (entered as lnAST)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL in univariate screening (Table 2) but not retained alongside the collinear BSA (Figure 2 correlation 0.97).",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL in univariate screening (Table 2, 'CR') but not retained in the final model.",
      source_name = "CR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 1L,
    n_observations = 139L,
    age_range = "0.61-14 years (mean 6.58, median 6.2)",
    age_median = "6.2 years",
    weight_range = "7.5-54 kg (mean 21.7, median 20)",
    weight_median = "20 kg",
    bsa_range = "0.38-1.50 m^2 (mean 0.80, median 0.8)",
    sex_female_pct = 35.4,
    disease_state = paste(
      "Children undergoing allogeneic hematopoietic stem cell",
      "transplantation receiving caspofungin as antifungal prophylaxis.",
      "37 normal hepatic function (AST <= 40 U/L), 9 mild hepatic",
      "dysfunction (40 < AST <= 120 U/L), 2 severe (AST > 120 U/L)."
    ),
    dose_range = paste(
      "Caspofungin 70 mg/m^2 loading dose on day 1 followed by",
      "50 mg/m^2 once daily, each as a 1-hour intravenous infusion."
    ),
    regions = "China (Wuhan Children's Hospital; July 2017 - March 2019)",
    notes = paste(
      "Prospective open-label study with opportunistic sampling at",
      "steady state (2.9 samples per patient; mean 17.5 h after the last",
      "dose; 104 samples between peak and trough and 35 troughs).",
      "Concentrations 4.5-17.4 mg/L by HPLC-UV (LLOQ 0.6 mg/L). Phoenix",
      "NLME 8.1; forward selection at P < 0.05, backward elimination at",
      "P < 0.01. A two-compartment model was tested but not retained."
    )
  )

  ini({
    # Structural typical values -- Niu 2020 Table 3 'Final model' column and
    # the final-model equations in Results ('Population Pharmacokinetic
    # Model Building'). Reference subject: BSA = 0.79 m^2, ln(AST) = 3.38.
    lcl <- log(0.14); label("Clearance at BSA 0.79 m^2 and ln(AST) 3.38 (L/h)")                  # Table 3 'theta CL' = 0.14 L/h (SE 8.50%); equation CL = 0.14 x ...
    lvc <- log(1.36); label("Volume of distribution at BSA 0.79 m^2 (L)")                        # Table 3 'theta Vd' = 1.36 L (SE 15.58%); equation Vd = 1.36 x (BSA/0.79)

    e_bsa_cl <- 0.89; label("Power exponent of BSA on CL (unitless)")                            # Table 3 'theta 1' = 0.89 (SE 11.36%); equation (BSA/0.79)^0.89 on CL
    e_ast_cl <- -0.23; label("Power exponent of ln(AST) on CL (unitless)")                       # Table 3 'theta 2' = -0.23 (SE 39.13%); equation (lnAST/3.38)^-0.23
    e_bsa_vc <- fixed(1); label("Power exponent of BSA on Vd (unitless)")                        # Results equation Vd = 1.36 x (BSA/0.79): exponent 1, not estimated (absent from Table 3)

    # Inter-individual variability. Table 3 footnote defines omega as the
    # 'square root of inter-individual variance', so variance = omega^2:
    #   Vd: 0.329^2 = 0.1082; CL: 0.333^2 = 0.1109. No covariance reported.
    etalcl ~ 0.1109                                                                              # Table 3 'omega CL (%)' = 33.3
    etalvc ~ 0.1082                                                                              # Table 3 'omega Vd (%)' = 32.9

    propSd <- 0.266; label("Proportional residual error (fraction)")                             # Table 3 'Residual variability sigma (%)' = 26.6
  })
  model({
    cl <- exp(lcl + etalcl) * (BSA / 0.79)^e_bsa_cl * (log(AST) / 3.38)^e_ast_cl
    vc <- exp(lvc + etalvc) * (BSA / 0.79)^e_bsa_vc

    kel <- cl / vc

    d / dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
