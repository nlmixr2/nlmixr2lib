Stockmann_2015_vancomycin <- function() {
  description <- paste(
    "One-compartment IV population PK model for vancomycin in neonates, as implemented and externally",
    "validated by Stockmann 2015. Clearance is driven by three covariates: allometric body weight",
    "(exponent 0.75, reference 2.9 kg), a sigmoidal postmenstrual-age maturation function",
    "(TM50 = 34.8 weeks, Hill = 4.53), and serum creatinine entering as (1/CREAT)^0.267 with CREAT in",
    "mg/dL on the Jaffe scale. Central volume scales linearly with weight (1.75 L at 2.9 kg). The",
    "structural and variance parameters were NOT estimated in Stockmann 2015; they are fixed priors",
    "carried verbatim from the originating model of Frymoyer 2014, which Stockmann 2015 re-implemented",
    "in NONMEM 7.2 and evaluated against an independent cohort of 243 neonates.",
    sep = " "
  )
  reference <- paste(
    "Stockmann C, Hersh AL, Roberts JK, Bhongsatiern J, Korgenski EK, Spigarelli MG, Sherwin CMT,",
    "Frymoyer A. Predictive Performance of a Vancomycin Population Pharmacokinetic Model in Neonates.",
    "Infect Dis Ther. 2015;4(2):187-198. doi:10.1007/s40121-015-0067-9",
    "(Methods, Model Evaluation: Equations 1 and 2 and the accompanying variance text).",
    "The structural model and every parameter estimate originate from Frymoyer A, Hersh AL,",
    "El-Komy MH, et al. Association Between Vancomycin Trough Concentration and AUC in Neonates.",
    "Antimicrob Agents Chemother. 2014 (cited as reference 7 of Stockmann 2015, whose reference list",
    "gives no volume, pages or DOI). Stockmann 2015 reprints the complete model, so no value here is",
    "taken from a source other than Stockmann 2015 itself.",
    sep = " "
  )
  vignette <- "Stockmann_2015_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Enters twice: allometrically on clearance as (WT/2.9)^0.75 (Equation 1) and",
        "linearly on central volume as (WT/2.9) (Equation 2). The 2.9 kg reference is the MEDIAN",
        "weight of the model development cohort in Table 1 ('Weight, kg', model development cohort",
        "column, median 2.9, range 0.5-6.3), which corroborates the transcription of both equations.",
        "The external validation cohort of Stockmann 2015 was lighter: median 1.6 kg (range 0.4-6.8).",
        sep = " "
      ),
      source_name        = "Weight"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WEEKS, not the register-default months. Equation 1 writes PMA explicitly as 'PMA_weeks' and",
        "its maturation constant TM50 = 34.8 is only meaningful on a week scale (see the PAGE entry",
        "in inst/references/covariate-columns.md, which provides for this). Time-varying. Drives the",
        "sigmoidal maturation fraction 1 / (1 + (PMA/34.8)^-4.53) on clearance, described in the text",
        "as 'an indicator of maturation'. Table 1: model development cohort median 39 weeks (range",
        "24-53); external validation cohort median 33 weeks (range 23-54). Study inclusion required",
        "< 54 weeks postmenstrual age.",
        sep = " "
      ),
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine, Jaffe-method scale",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Equation 1 writes the term as (1 / Cr_mg/dL)^0.267, i.e. the covariate is",
        "NOT normalised to a reference creatinine; the factor equals 1 at CREAT = 1 mg/dL, so the",
        "typical value 0.345 L/h is the clearance of a fully mature 2.9 kg neonate with a serum",
        "creatinine of 1 mg/dL. The unit mg/dL is pinned by the equation's own subscript.",
        "ASSAY SCALE IS LOAD-BEARING: Table 1 footnote c states the model derivation cohort's",
        "creatinine was measured by the JAFFE method. Stockmann 2015 measured its validation cohort",
        "by the enzymatic method and converted to a Jaffe-standardized equivalent with the Srivastava",
        "linear equation printed in Methods, Validation Cohort: enzymatic = 1.050 * Jaffe - 0.122,",
        "i.e. Jaffe = (enzymatic + 0.122) / 1.050. Supply Jaffe-scale values, or convert first.",
        "Table 1: model development cohort median 0.4 mg/dL (range 0.1-2.7); external validation",
        "cohort median 0.6 mg/dL (range 0.3-1.5), Jaffe-standardized.",
        sep = " "
      ),
      source_name        = "Cr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 249L,
    n_studies      = 1L,
    age_range      = "24-53 weeks postmenstrual age",
    age_median     = "39 weeks postmenstrual age",
    weight_range   = "0.5-6.3 kg",
    weight_median  = "2.9 kg",
    sex_female_pct = 49.0,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Neonates receiving intravenous vancomycin with therapeutic drug monitoring performed.",
      "Neonates with congenital kidney disease, major congenital heart disease (other than",
      "ventricular septal defect, atrial septal defect or patent ductus arteriosus), or",
      "extracorporeal membrane oxygenation during the vancomycin course were excluded from the",
      "Stockmann 2015 validation cohort.",
      sep = " "
    ),
    dose_range     = paste(
      "Intravenous vancomycin as 1-h infusions per routine clinical practice. In the Stockmann 2015",
      "validation cohort the median dose was 15.5 mg/kg (IQR 13.9-19.3) and the median dosing",
      "interval 11.5 h (IQR 8.0-12.5).",
      sep = " "
    ),
    regions        = "United States",
    ga_range       = "23-42 weeks gestational age (median 34)",
    renal_function = "Serum creatinine median 0.4 mg/dL (range 0.1-2.7), Jaffe method",
    notes          = paste(
      "DEVELOPMENT population, from the 'Model development cohort (n = 249)' column of Stockmann 2015",
      "Table 1, whose footnote a identifies it as the cohort used to develop the model in Frymoyer",
      "2014 (reference 7). This is the cohort the parameters were estimated from, and the 2.9 kg",
      "reference weight in Equations 1 and 2 is its median weight. Stockmann 2015 estimated nothing;",
      "it set every structural and variance parameter equal to the previously reported estimates.",
      sep = " "
    ),
    validation_cohort = paste(
      "Stockmann 2015's own contribution: 243 neonates with 734 vancomycin concentrations (mean 3.0",
      "+/- 1.8 per neonate) from five Intermountain Healthcare neonatal intensive care units,",
      "2006-2013. Table 1, external validation cohort column: 42% female; gestational age median 30",
      "weeks (range 22-41); birth weight median 1.3 kg (range 0.5-5.1); weight median 1.6 kg (range",
      "0.4-6.8); postnatal age median 12 days (range 0-196); postmenstrual age median 33 weeks (range",
      "23-54); serum creatinine median 0.6 mg/dL (range 0.3-1.5). Lower weight and age and higher",
      "serum creatinine than the development cohort. Predictive performance (Table 3, all",
      "concentrations): median prediction error -0.8 mg/L (95% CI -1.4 to -0.4), median absolute",
      "prediction error 3.0 mg/L (95% CI 2.7-3.5), median absolute percent prediction error 15.2%;",
      "NPDE mean 0.05, variance 0.96.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Every parameter is FIXED. Stockmann 2015 Methods, Model Evaluation: "The
    # published population pharmacokinetic model was implemented in ... NONMEM
    # 7.2 ... as previously described [7]" and the Abstract: "with the
    # structural and variance parameter values set equal to the estimates
    # reported previously". No parameter was estimated in this paper, so none
    # carries a standard error, RSE or confidence interval anywhere in it.
    #
    # Typical subject implied by Equation 1: a 2.9 kg neonate of fully mature
    # postmenstrual age with a Jaffe-method serum creatinine of 1 mg/dL.
    # ------------------------------------------------------------------------
    lcl <- fixed(log(0.345)); label("Clearance at WT = 2.9 kg, mature PMA and CREAT = 1 mg/dL (L/h)")
    # Equation 1, leading coefficient: CL (L/h) = 0.345 * (Weight/2.9 kg)^0.75 * ...
    lvc <- fixed(log(1.75));  label("Central volume of distribution at WT = 2.9 kg (L)")
    # Equation 2, leading coefficient: V (L) = 1.75 * (Weight/2.9 kg)

    # Covariate effects, all read directly off Equations 1 and 2.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/2.9) for CL (unitless)")
    # Equation 1, superscript on the (Weight/2.9 kg) factor
    # NOTE: Equation 2 carries NO exponent on its (Weight/2.9 kg) factor, i.e.
    # central volume is scaled linearly with weight. That unity exponent is
    # written literally as (WT / 2.9) in model() rather than as a fixed(1)
    # parameter, so the model body matches the printed equation term for term.

    pma_tm50 <- fixed(34.8); label("Postmenstrual age at 50% of mature clearance (weeks)")
    # Equation 1, denominator of the (PMA_weeks / 34.8) ratio inside the maturation term
    pma_hill <- fixed(4.53); label("Hill coefficient of the postmenstrual-age maturation function (unitless)")
    # Equation 1, magnitude of the NEGATIVE exponent -4.53 on (PMA_weeks / 34.8)

    e_creat_cl <- fixed(0.267); label("Power exponent on (1 mg/dL / CREAT) for CL (unitless)")
    # Equation 1, superscript on the (1 / Cr_mg/dL) factor

    # Inter-individual variability. Methods: "the remaining variation between
    # neonates was described by an exponential error model for both CL (%
    # coefficient of variation [% CV] 21.6%) and V (% CV 10.9%)". Exponential
    # IIV, so the reported %CV is the omega SD on the log scale under the usual
    # NONMEM reporting convention omega^2 = CV^2; the variances below are
    # 0.216^2 and 0.109^2. The alternative log-normal conversion
    # omega^2 = log(CV^2 + 1) would give 0.04560 and 0.01184, a 1% and 0.3%
    # difference on the SD scale -- immaterial here and discussed in the
    # vignette Errata. No correlation between the two etas is reported.
    # Source traces are on their own lines: rxode2 rewrites a trailing comment
    # on an `eta ~ ...` line into a label() call.
    # Methods, Model Evaluation: CL % CV 21.6% -> 0.216^2 = 0.046656
    etalcl ~ fixed(0.046656)
    # Methods, Model Evaluation: V % CV 10.9% -> 0.109^2 = 0.011881
    etalvc ~ fixed(0.011881)

    # Combined proportional + additive residual error. Methods: "Residual
    # variability ... was captured using a combined proportional (% CV 20.5%)
    # and additive error model (standard deviation [SD] +/- 1.3 mg/L)".
    propSd <- fixed(0.205); label("Proportional residual error (fraction)")
    # Methods, Model Evaluation: proportional % CV 20.5%
    addSd  <- fixed(1.3);   label("Additive residual error (mg/L)")
    # Methods, Model Evaluation: additive standard deviation 1.3 mg/L
  })
  model({
    # 1. Derived covariate terms.
    #
    #    Maturation, written exactly as Equation 1 prints it:
    #
    #                             1
    #      ----------------------------------------------
    #      1 + (PMA_weeks / 34.8) ^ (-4.53)
    #
    #    This is the standard sigmoidal Hill maturation function; it is
    #    algebraically identical to PMA^4.53 / (34.8^4.53 + PMA^4.53), equals
    #    0.5 at PMA = 34.8 weeks and approaches 1 at large PMA. The printed
    #    form is kept so the body can be checked against the source term by
    #    term. The exponent's minus sign was confirmed against the typeset
    #    equation in the PDF, not only the extracted text.
    maturation_cl <- 1 / (1 + (PAGE / pma_tm50)^(-pma_hill))

    #    Renal function, from the (1 / Cr_mg/dL)^0.267 factor of Equation 1.
    #    Not normalised to a reference creatinine in the source; the factor is
    #    1 at CREAT = 1 mg/dL.
    creat_cl <- (1 / CREAT)^e_creat_cl

    # 2. Individual PK parameters.
    cl <- exp(lcl + etalcl) * (WT / 2.9)^e_wt_cl * maturation_cl * creat_cl
    vc <- exp(lvc + etalvc) * (WT / 2.9)

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. ODE system. Methods: "a one compartment model with first-order
    #    elimination was used to describe vancomycin pharmacokinetics".
    #    Vancomycin is given intravenously, so there is no absorption
    #    compartment; dose enters `central` directly. Table 2's heading
    #    establishes that doses were administered as 1-h infusions, which is a
    #    property of the event table rather than of the model.
    d/dt(central) <- -kel * central

    # 5. Observation and error. Dose in mg over volume in L gives mg/L, which
    #    matches the mg/L unit of the additive residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
