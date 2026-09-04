Leeds_2013_tecovirimat_human <- function() {
  description <- "Two-compartment population PK model for oral tecovirimat (ST-246) in healthy adult volunteers (Leeds 2013), with first-order absorption preceded by a lag time, allometric body-weight scaling (0.75 on CL/F and Q/F, 1 on Vc/F and Vp/F) about a 78.4 kg mean weight, and a sex effect on Vc/F (29.5% larger apparent central volume in women). Companion cynomolgus-monkey model used for the animal-rule dose bridge: modellib('Leeds_2013_tecovirimat_cyno')."
  reference <- paste(
    "Leeds JM, Fenneteau F, Gosselin NH, Mouksassi MS, Kassir N, Marier JF,",
    "Chen Y, Grosenbach D, Frimm AE, Honeychurch KM, Chinsangaram J,",
    "Tyavanagimatt SR, Hruby DE, Jordan R. (2013). Pharmacokinetic and",
    "pharmacodynamic modeling to determine the dose of ST-246 to protect",
    "against smallpox in humans. Antimicrob Agents Chemother 57(3):1136-1143.",
    "doi:10.1128/AAC.00959-12.",
    sep = " "
  )
  vignette <- "Leeds_2013_tecovirimat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "tecovirimat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tecovirimat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tecovirimat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with theory-based fixed exponents about 78.4 kg, the normalisation constant printed inside every human covariate formula of Leeds 2013 Table 1. Results 'Human POP PK model development' names it explicitly: 'For a healthy human with a body weight equivalent to the mean weight from the clinical study used in development of the POP PK model, 78.4 kg, the typical values of CL/F would be 41.15 liters/h ... and Vc/F in females and males would be 281.51 liters and 217.44 liters, respectively'. The exponents are the same theory-based 0.75 / 1 pair used in the companion NHP model.",
      source_name        = "wt"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Leeds 2013 Table 1 prints the human Vc/F row as two parallel typical values, '281.51 * (wt/78.4)^1 (female); 217.44 * (wt/78.4)^1 (male)', so the effect is recovered here as the log ratio of the two, with male taken as the reference. Table 1 footnote a: 'in humans, gender was a covariate for Vc/F, but not in NHPs'. The clinical study enrolled approximately equal numbers of male and female subjects. The paper judges the effect clinically unimportant because it barely moves the terminal half-life -- Results 'Human POP PK model development': 'The terminal elimination half-lives in male and female subjects were similar (16.7 h in males versus 17.4 h in females). Consequently, no clinically significant sex effect was predicted for ST-246 exposure or the terminal elimination half-life.' Sex was not retained in the companion NHP model.",
      source_name        = "gender"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 91L,
    n_studies      = 1L,
    age_range      = "19-74 years",
    weight_mean    = "78.4 kg",
    sex_female_pct = 50,
    disease_state  = "Healthy volunteers.",
    dose_range     = "Oral ST-246 400 mg (45 subjects) or 600 mg (46 subjects) once daily for 14 days, administered in the nonfasted state.",
    notes          = "Double-blind, randomised, placebo-controlled, multicenter phase 1 safety / tolerability / PK trial; 107 volunteers were randomised in three cohorts, of whom 91 received ST-246 and 16 received placebo, so only the 91 active-treated subjects contribute concentrations. Subjects with fewer than two concentration-time data were excluded (3 subjects in total), as were outlier concentrations showing high troughs or double peaks at 24 h postdose; 1,395 plasma concentrations of ST-246 entered the human population PK analysis, none with |CWRES| > 4. Sex distribution is given as 'approximately equal numbers of male and female subjects'; the exact percentage is not tabulated. Below-quantitation-limit values were handled by Beal's method M6. Model fit in Phoenix NLME 6.2.0.416 with FOCE-ELS. Parameters from Leeds 2013 Table 1, human columns; population description from Materials and Methods 'In vivo study summaries (ii) Human clinical PK study'."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural PK parameters -- Leeds 2013 Table 1, 'Human PK parameter /
    # Population estimate' column. All human formulae in Table 1 are printed
    # normalised to wt = 78.4 kg, so the typical values below are the values
    # at WT = 78.4 kg. Vc/F is printed as separate female and male values;
    # male is taken as the reference here and the female shift is a
    # covariate-effect parameter below.
    # ---------------------------------------------------------------------
    lka   <- log(1.06)   ; label("Absorption rate constant Ka (1/h)")                                       # Table 1 row 'Ka (h-1)' Human Population estimate = 1.06
    ltlag <- log(1.46)   ; label("Absorption lag time Tlag (h)")                                            # Table 1 row 'Tlag (h)' Human Population estimate = 1.46
    lcl   <- log(41.15)  ; label("Apparent oral clearance CL/F at WT = 78.4 kg (L/h)")                      # Table 1 row 'CL/F (liters/h)' Human: 41.15 * (wt/78.4)^0.75
    lvc   <- log(217.44) ; label("Apparent central volume of distribution Vc/F in males at WT = 78.4 kg (L)") # Table 1 row 'Vc/F (liters)' Human: 217.44 * (wt/78.4)^1 (male)
    lq    <- log(36.79)  ; label("Apparent inter-compartmental clearance Q/F at WT = 78.4 kg (L/h)")        # Table 1 row 'Q/F (liters/h)' Human: 36.79 * (wt/78.4)^0.75
    lvp   <- log(413.53) ; label("Apparent peripheral volume of distribution Vp/F at WT = 78.4 kg (L)")     # Table 1 row 'Vp/F (liters)' Human: 413.53 * (wt/78.4)^1

    # ---------------------------------------------------------------------
    # Allometric exponents, shared across the clearance terms and across the
    # volume terms. Theory-based, printed as literal exponents inside the
    # Table 1 formulae and without uncertainty, and described in Results
    # 'NHP POP PK model development' as 'theta_eff equal to 0.75 for
    # clearance-related parameters and theta_eff equal to 1 for
    # volume-related parameters (24)'; Results 'Human POP PK model
    # development' states the human model used 'allometric factors on PK
    # parameters using wt as a covariate (Table 1), similar to the NHP
    # model'. Encoded as structural rather than estimated.
    # ---------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F with body weight (unitless)")  # Table 1 human CL/F and Q/F formulae
    e_wt_vc_vp <- fixed(1.0)  ; label("Allometric exponent on Vc/F and Vp/F with body weight (unitless)") # Table 1 human Vc/F and Vp/F formulae

    # ---------------------------------------------------------------------
    # Sex effect on the apparent central volume, recovered as the log ratio
    # of the two typical values Table 1 prints for the same row. The weight
    # term is identical in both, so the ratio is exactly the sex effect.
    # ---------------------------------------------------------------------
    e_sexf_vc <- log(281.51 / 217.44) ; label("Log-additive effect of female sex on Vc/F (unitless)") # Table 1 row 'Vc/F (liters)' Human: 281.51 (female) vs 217.44 (male), i.e. 29.5% larger in women

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 1 column 'IIV (%)', human side.
    # Read as the standard deviation of the log-scale random effect
    # expressed as a percentage, i.e. omega^2 = (IIV/100)^2 -- the same
    # reading used for the companion NHP model file and for the '% CV' IIV
    # column of Jonsson_2011_ethambutol.R. See vignette Errata.
    # ---------------------------------------------------------------------
    etalka   ~ 0.1681 # Table 1 'Ka (h-1)' Human IIV = 41%; omega^2 = 0.41^2
    etaltlag ~ 0.0289 # Table 1 'Tlag (h)' Human IIV = 17%; omega^2 = 0.17^2
    etalcl   ~ 0.0961 # Table 1 'CL/F' Human IIV      = 31%; omega^2 = 0.31^2
    etalvc   ~ 0.0784 # Table 1 'Vc/F' Human IIV      = 28%; omega^2 = 0.28^2
    etalq    ~ 0.2916 # Table 1 'Q/F' Human IIV       = 54%; omega^2 = 0.54^2
    etalvp   ~ 0.2916 # Table 1 'Vp/F' Human IIV      = 54%; omega^2 = 0.54^2

    # ---------------------------------------------------------------------
    # Combined residual error. Table 1 'Error model' block, human column.
    # The additive term is printed in ug/liter and is converted to the
    # model's mg/L concentration scale by dividing by 1000.
    # ---------------------------------------------------------------------
    propSd <- 0.27          ; label("Proportional residual error (fraction)") # Table 1 'Error model / Proportional (%)' Human = 27
    addSd  <- 10.92 / 1000  ; label("Additive residual error (mg/L)")         # Table 1 'Error model / Additive (ug/liter)' Human = 10.92 ug/L
  })

  model({
    # Individual PK parameters, with allometric weight terms normalised to
    # the 78.4 kg cohort mean exactly as Table 1 prints them.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 78.4)^e_wt_cl_q
    vc <- exp(lvc + e_sexf_vc * SEXF + etalvc) * (WT / 78.4)^e_wt_vc_vp
    q  <- exp(lq + etalq) * (WT / 78.4)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 78.4)^e_wt_vc_vp

    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order oral absorption from a depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag applied to the dosing compartment.
    alag(depot) <- tlag

    # Dose units mg, Vc/F units L, so Cc units mg/L (= ug/mL).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
