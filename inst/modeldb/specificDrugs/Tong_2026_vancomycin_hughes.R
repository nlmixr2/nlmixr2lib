Tong_2026_vancomycin_hughes <- function() {
  description <- paste(
    "Hughes two-compartment IV population PK model for vancomycin in adults with obesity as implemented",
    "for model-informed precision dosing (MIPD) in the InsightRX Nova clinical decision support",
    "software, and used as the POST-intervention default model for patients with BMI >= 40 kg/m2 in",
    "Tong 2026. Fat-free mass is computed inside the model by the Janmahasatian equations from body",
    "weight, height and sex; creatinine clearance is then computed by a Cockcroft-Gault form that",
    "substitutes fat-free mass for total body weight. Clearance scales as (CRCL/100)^0.851 and both",
    "volumes scale allometrically on (FFM/70) with a fixed linear exponent. Correlated",
    "inter-individual variability on CL, Vc and Vp. All population parameters are FIXED priors for MAP",
    "Bayesian estimation; none were estimated in Tong 2026.",
    sep = " "
  )
  reference <- paste(
    "Tong DMH, Brooks JT, Keizer RJ, Hughes JH. Vancomycin target attainment improved following",
    "population pharmacokinetic model switch: a large-scale quasi-experimental study of precision",
    "dosing. JAC Antimicrob Resist. 2026. doi:10.1093/jacamr/dlag016 (Supplementary data, Code",
    "section, \"Hughes model\" NONMEM control stream; Table S2).",
    "Structural model and parameter estimates originate from Hughes MSA, Hughes JH, Endicott J et al.",
    "Developing parametric and nonparametric models for model-informed precision dosing: a quality",
    "improvement effort in vancomycin for patients with obesity. Ther Drug Monit 2024;46:575-583.",
    "doi:10.1097/FTD.0000000000001214.",
    "The fat-free-mass equations are Janmahasatian S, Duffull SB, Ash S et al.",
    "Clin Pharmacokinet 2005;44:1051-1065.",
    sep = " "
  )
  vignette <- "Tong_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI >= 40 kg/m2 cohort: median 133.0 kg (range 43.8-318.0). Does NOT enter",
        "the disposition parameters directly. It is an input to the two internally derived quantities",
        "that do: body mass index (BMI = WT / (HT/100)^2) and, through that, Janmahasatian fat-free",
        "mass. Using FFM rather than total weight is the point of this model in an obese cohort.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI >= 40 kg/m2 cohort: median 167.6 cm (range 93.3-200.7). The control",
        "stream divides by 100 to obtain metres before squaring (BMI = WT / ((HT/100)**2)), which",
        "fixes the unit as cm. Enters only through BMI and hence fat-free mass.",
        sep = " "
      ),
      source_name        = "HT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI >= 40 kg/m2 cohort: median 60.0 years (range 18.2 to 90+; ages above",
        "90 are aggregated for de-identification). Used only inside the internal fat-free-mass-based",
        "Cockcroft-Gault creatinine-clearance calculation.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tong 2026 Table 1, BMI >= 40 kg/m2 cohort: 1272 male / 1692 female treatment courses (57.1%",
        "female). The supplement control stream uses the OPPOSITE polarity, SEX with 1 = male: the",
        "default fat-free-mass line uses the male Janmahasatian coefficients (6680, 216) and",
        "IF(SEX.EQ.0) overrides them with the female coefficients (8780, 244), and the",
        "Cockcroft-Gault term is written 0.85**(1-SEX) so that females receive the 0.85 factor.",
        "Converted to the canonical SEXF (1 = female) via SEXF = 1 - SEX, so the sex-specific",
        "fat-free-mass coefficients are selected by SEXF and the Cockcroft-Gault factor becomes",
        "0.85^SEXF.",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI >= 40 kg/m2 cohort: median 0.95 mg/dL (range 0.10-10.9). Used as the",
        "denominator of the Cockcroft-Gault equation with the 72 constant, which fixes the unit as",
        "mg/dL. NOTE: unlike the sibling modified-Goti stream, this stream applies NO cap to the",
        "resulting creatinine clearance -- the (CRCL/100)^0.8509 power term is applied to the raw",
        "computed value.",
        sep = " "
      ),
      source_name        = "CR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2709L,
    n_studies      = 19L,
    age_range      = "18.2 to 90+ years",
    age_median     = "60.0 years",
    weight_range   = "43.8-318.0 kg",
    weight_median  = "133.0 kg",
    sex_female_pct = 57.1,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Hospitalized adults (>= 18 years) with BMI >= 40 kg/m2 receiving intravenous vancomycin under",
      "routine model-informed precision dosing; at least two doses and at least one measured",
      "concentration required. Patients undergoing haemodialysis at any point during treatment were",
      "excluded, as were patients dosed with a model other than their site's default.",
      sep = " "
    ),
    dose_range     = "Intravenous vancomycin per routine clinical practice; initial doses selected a priori, subsequent doses adapted by MAP Bayesian posterior estimates",
    regions        = "United States (19 hospital systems, patients beginning treatment August 2022 to December 2024)",
    renal_function = "Serum creatinine median 0.95 mg/dL (range 0.10-10.9); haemodialysis patients excluded",
    n_concentrations = 6572L,
    notes          = paste(
      "APPLICATION population from Tong 2026 Table 1 (BMI >= 40 kg/m2 cohort: 2709 patients, 2964",
      "treatment courses, 6572 samples), i.e. the cohort this model was USED to dose as the",
      "post-intervention default -- NOT the cohort it was estimated from. The DEVELOPMENT population",
      "is Hughes 2024: 100 patients with obesity and 276 concentrations (Tong 2026 Table S2 rows",
      "'# patients in development population' = 100 and '# of samples' = 276). Tong 2026 Discussion",
      "contrasts this model favourably with the Carreno model it replaced -- 'built as a parametric",
      "model on a larger data set of 100 patients, with a simpler structure, fewer variability terms,",
      "and allometric scaling by fat-free mass' -- while noting that the switch did not significantly",
      "change any AUC metric. Every population parameter is a FIXED prior inherited from Hughes 2024;",
      "Tong 2026 estimated none of them.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # ALL population parameters are FIXED priors for MAP Bayesian estimation.
    # Source: Tong 2026 Supplementary data, Code section, "Hughes model" control
    # stream. NOTE: only THETA(7) carries an explicit FIX flag in that listing;
    # the remaining $THETA records and the $OMEGA BLOCK(3) do not, unlike the
    # sibling modified-Goti and modified-Thomson streams. They are nevertheless
    # fixed priors, because Tong 2026 estimated no population parameters
    # (Methods, Pharmacokinetic analysis: only individual parameters were fitted,
    # by MAP Bayesian estimation) -- the un-flagged records are the estimation
    # control stream of the upstream Hughes 2024 fit, carried over verbatim.
    # Encoded with fixed() accordingly -- see vignette Errata.
    # ------------------------------------------------------------------------
    lcl <- fixed(log(5.0627)); label("Clearance at CRCL = 100 mL/min (L/h)")
    # $THETA(1) 5.0627 ; CL. Table S2 rounds it: "CL = 5.06 * (CRCL/100)^0.8509"
    lvc <- fixed(log(64.339)); label("Central volume at FFM = 70 kg (L)")
    # $THETA(2) 64.339 ; V1. Table S2 "V = 64.339 * FFM/70"
    lq  <- fixed(log(6.2402)); label("Intercompartmental clearance Q (L/h)")
    # $THETA(3) 6.2402 ; Q. Table S2 "Q = 6.2402"
    lvp <- fixed(log(59.019)); label("Peripheral volume at FFM = 70 kg (L)")
    # $THETA(4) 59.019 ; V2. Table S2 "V2 = 59.019 * FFM/70"

    # Covariate effects:
    #   TVCL  = THETA(1) * (CRCL/100)**THETA(8)
    #   TVV1  = THETA(2) * (FFM/70)**THETA(7)
    #   TVV2  = THETA(4) * (FFM/70)**THETA(7)
    e_crcl_cl    <- fixed(0.8509); label("Power exponent on (CRCL/100) for CL (unitless)")
    # $THETA(8) 0.8509 ; CRCL ~ CL. Table S2 exponent in "(CRCL/100)^0.8509"
    e_ffm_vc_vp  <- fixed(1);      label("Allometric exponent on (FFM/70) for Vc and Vp (unitless)")
    # $THETA(7) 1 FIX ; "allo V1/V2" -- the one theta explicitly flagged FIX in the listing.
    # A linear exponent, consistent with Table S2 writing the volumes as plain
    # proportionalities "V = 64.339 * FFM/70" and "V2 = 59.019 * FFM/70".

    # Inter-individual variability, $OMEGA BLOCK(3), lower triangle by row:
    #   0.0473154
    #   0.00271032  0.0322236
    #   0.0132618   0.138593   0.634656
    # mapping to ETA(1) on CL, ETA(2) on V1 and ETA(3) on V2. The stream stores
    # variances on the omega^2 = CV^2 convention (NOT log(CV^2 + 1)): each
    # sqrt(diagonal) reproduces the Table S2 %CV row exactly. The Vc-Vp
    # correlation implied by the off-diagonals is 0.969, which is high but leaves
    # the matrix positive-definite (smallest eigenvalue 0.00187).
    etalcl + etalvc + etalvp ~ fixed(c(
      0.0473154,
      0.00271032, 0.0322236,
      0.0132618,  0.138593,  0.634656
    ))
    # sqrt(0.0473154) = 0.2175 -> Table S2 "IIV on CL (%CV) 21.8";
    # sqrt(0.0322236) = 0.1795 -> Table S2 "IIV on V (%CV) 18.0";
    # sqrt(0.634656)  = 0.7967 -> Table S2 "IIV on V2 (%CV) 79.7".
    # No eta on Q: the stream sets Q = TVQ with no EXP(ETA), and Table S2 reports
    # "IIV on Q (%CV)" as "-" for this model.

    # Combined additive + proportional residual error, from the $ERROR block
    # W = SQRT(IPRED**2 * PROP**2 + ADD**2), with PROP and ADD carried as thetas
    # rather than hardcoded (unlike the three sibling streams).
    addSd  <- fixed(0.001);  label("Additive residual error (mg/L)")
    # $THETA(6) 0.001 ; add error. Effectively nil, consistent with Table S2
    # reporting "Additive error (mg/L)" as "-" for this model.
    propSd <- fixed(0.1372); label("Proportional residual error (fraction)")
    # $THETA(5) 0.1372 ; prop error. Table S2 reports "Proportional error (%) 13.3";
    # the control-stream value is used as authoritative -- see vignette Errata.
  })
  model({
    # 1. Derived covariate terms, transcribed from $PK. Body mass index, then
    #    Janmahasatian fat-free mass with sex-specific coefficients, then a
    #    Cockcroft-Gault creatinine clearance that substitutes FFM for total body
    #    weight:
    #      BMI  = WT / ((HT/100)**2)
    #      FFM  = 9270 * WT / (6680 + 216*BMI)            [male,   SEX = 1]
    #      FFM  = 9270 * WT / (8780 + 244*BMI)            [female, SEX = 0]
    #      CRCL = (140-AGE) * FFM * 0.85**(1-SEX) / (72*CR)
    #    with SEXF = 1 - SEX, so SEXF selects the female branch and the
    #    Cockcroft-Gault factor becomes 0.85^SEXF.
    bmi_i    <- WT / (HT / 100)^2
    ffm_male <- 9270 * WT / (6680 + 216 * bmi_i)
    ffm_fem  <- 9270 * WT / (8780 + 244 * bmi_i)
    ffm_i    <- ffm_male + SEXF * (ffm_fem - ffm_male)
    crcl_i   <- (140 - AGE) * ffm_i * 0.85^SEXF / (72 * CREAT)

    # 2. Individual PK parameters. Note that CL scales on renal function only and
    #    the volumes on fat-free mass only; Q carries neither a covariate nor an
    #    eta.
    cl <- exp(lcl + etalcl) * (crcl_i / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (ffm_i / 70)^e_ffm_vc_vp
    vp <- exp(lvp + etalvp) * (ffm_i / 70)^e_ffm_vc_vp
    q  <- exp(lq)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. The source stream uses the closed-form $SUBROUTINE ADVAN3
    #    TRANS4 (two-compartment, IV) rather than an explicit $DES; the ODEs below
    #    are the standard ADVAN3 system and are written out for house
    #    consistency with the three sibling models. Unlike those siblings this
    #    stream has no CL_HEMO term and, being ADVAN3, no AUC accumulator state
    #    (its $ERROR block nevertheless assigns AUC = A(3), a leftover reference
    #    to a state ADVAN3 does not define -- see vignette Errata).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation and error. S1 = V1 in the stream, so dose in mg over volume
    #    in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
