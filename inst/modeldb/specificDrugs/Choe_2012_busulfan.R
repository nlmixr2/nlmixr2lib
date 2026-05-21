Choe_2012_busulfan <- function() {
  description <- "One-compartment IV PK model for intravenous busulfan in adult Korean hematopoietic stem cell transplant recipients, with allometric scaling on actual body weight (fixed exponent 0.5) on CL and Vd and a sex effect on Vd (Choe 2012)."
  reference <- "Choe S, Kim G, Lim H-S, et al. A simple dosing scheme for intravenous busulfan based on retrospective population pharmacokinetic analysis in Korean patients. Korean J Physiol Pharmacol. 2012;16(4):273-280. doi:10.4196/kjpp.2012.16.4.273"
  vignette <- "Choe_2012_busulfan"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in the source analysis. Allometric exponent on both CL and Vd was fixed at 0.5; reference weight 65 kg corresponds to the paper's typical-patient typical values CL 7.6 L/h and Vd 32.2 L (male) / 29.1 L (female) from Table 2 (Choe 2012 Results, p. 276).",
      source_name        = "ABW"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Canonical SEXF encoding (1 = female). The source paper used SEX with 1 = male, 0 = female (Choe 2012 Results, p. 276); the relationship is SEX_paper = 1 - SEXF. The published Vd structural value 3.610 * ABW^0.5 corresponds to the female reference, with male Vd a +10.5% deviation. The model preserves the female-reference structural value by applying the sex effect via (1 + e_sexf_vc * (1 - SEXF)).",
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60,
    n_studies      = 1,
    age_range      = "16-58 years",
    age_mean       = "36.5 years (SD 10.9)",
    weight_range   = "52.5-116 kg",
    weight_mean    = "66.5 kg (SD 11.3)",
    sex_female_pct = 38.3,
    race_ethnicity = c(Asian_Korean = 100),
    disease_state  = "Adult patients with hematologic malignancies (AML/acute mixed leukemia, ALL, CML, MDS, other) receiving intravenous busulfan as conditioning therapy prior to hematopoietic stem cell transplantation (BuCy, BuFluATG, or Bu-only regimens). Patients had adequate cardiac, hepatic and renal function and Karnofsky performance scores >=70.",
    dose_range     = "Either 0.8 mg/kg every 6 h over a 2 h IV infusion x 4 days (BU4 arm) or 3.2 mg/kg every 24 h over a 3 h IV infusion x 4 days (BU1 arm). Doses calculated using selected body weight (SBW: actual body weight if <= ideal body weight (IBW), IBW if ABW within 120% of IBW, adjusted IBW otherwise).",
    regions        = "Korea (single-center: Asan Medical Center, Seoul)",
    notes          = "Source: Choe 2012 Table 1 and Methods. 60 Korean adults (37 men, 23 women) enrolled 1:1 across BU4 and BU1 arms. 295 plasma busulfan concentrations measured by validated LC-MS/MS (LOQ 30 ng/mL); samples drawn at 2.5, 3, 4, 5, 6 h post-infusion start (BU4) or 3.5, 5, 6, 7, 22 h post-infusion start (BU1). Median AUC0-inf with the once-daily 3.2 mg/kg BU1 regimen was 6,378 umol/L*min = 26.18 mg/L*h. All subjects received concomitant phenytoin (15 mg/kg loading + maintenance) for seizure prophylaxis. NONMEM VI ADVAN1 TRANS2 with FOCEi."
  )

  ini({
    # Structural parameters at reference weight 65 kg. The paper's Table 3 reports
    # CL = theta1 * ABW^0.5 and Vd_female = theta2 * ABW^0.5, with theta1 = 0.947
    # and theta2 = 3.610. At ABW = 65 kg this gives CL = 7.6 L/h and Vd_female =
    # 29.1 L; the male Vd is +10.5% (= 32.2 L) via the sex effect below. The
    # reference-weight form (exp(lcl) at WT = 65) lets the size term collapse to 1
    # for a 65 kg subject so the typical values from Choe 2012 Table 2 read off
    # directly.
    lcl <- log(7.6)  ; label("Typical clearance at 65 kg (L/h)")              # Choe 2012 Table 2 + Table 3 (theta1 = 0.947; 0.947 * sqrt(65) = 7.6)
    lvc <- log(29.1) ; label("Typical Vd at 65 kg, female reference (L)")     # Choe 2012 Table 2 + Table 3 (theta2 = 3.610; 3.610 * sqrt(65) = 29.1)

    # Allometric exponents on actual body weight, both fixed at 0.5 per the
    # paper: "Simplified power model with the power terms fixed at 0.5 was
    # fitted to the data and estimation of power terms as unknown parameters
    # was found not to improve the model compared to fixing the power term."
    e_wt_cl <- fixed(0.5) ; label("Allometric exponent of WT on CL (unitless, fixed)")   # Choe 2012 Results, p. 276
    e_wt_vc <- fixed(0.5) ; label("Allometric exponent of WT on Vd (unitless, fixed)")   # Choe 2012 Results, p. 276

    # Sex effect on Vd. The paper parameterised Vd as theta2 * ABW^0.5 *
    # (1 + SEX * theta3) with SEX = 1 for male, 0 for female and theta3 =
    # 0.105 (estimated). The female reference value is therefore the
    # unmodified theta2 * ABW^0.5, and the male Vd is 10.5% larger. With the
    # canonical SEXF (1 = female, 0 = male), the same effect is applied as
    # (1 + e_sexf_vc * (1 - SEXF)) so the structural value preserves the
    # published female reference. The coefficient itself stays at 0.105.
    e_sexf_vc <- 0.105 ; label("Sex effect on Vd; applied on (1 - SEXF) to preserve female reference")  # Choe 2012 Table 3 (theta3 = 0.105)

    # Inter-individual variability. Reported as CV% in Table 2: 16% on CL and
    # 9% on Vd. Converted to the internal log-normal variance scale via
    # omega^2 = log(CV^2 + 1): log(0.16^2 + 1) = 0.025277, log(0.09^2 + 1) =
    # 0.008068. Shrinkage was 1% on etalcl and 9% on etalvc (Choe 2012
    # Results, p. 277). IIV reported as diagonal (no correlation given).
    etalcl ~ 0.025277  # Choe 2012 Table 2 (IIV CL 16%)
    etalvc ~ 0.008068  # Choe 2012 Table 2 (IIV Vd 9%)

    # Proportional residual error. Table 2 reports residual variability 6.3%
    # for the proportional error model Y = F * (1 + epsilon); propSd is the
    # standard deviation of epsilon (= 0.063).
    propSd <- 0.063 ; label("Proportional residual error (fraction)")  # Choe 2012 Table 2
  })

  model({
    # Size scaling on actual body weight, reference 65 kg.
    size_cl <- (WT / 65)^e_wt_cl
    size_vc <- (WT / 65)^e_wt_vc

    # Sex effect on Vd: female reference (SEXF = 1) -> multiplier = 1;
    # male (SEXF = 0) -> multiplier = 1 + 0.105 = 1.105.
    sex_vc <- 1 + e_sexf_vc * (1 - SEXF)

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * size_cl
    vc <- exp(lvc + etalvc) * size_vc * sex_vc

    # One-compartment IV (NONMEM ADVAN1 TRANS2 equivalent). Busulfan is
    # administered as an IV infusion directly into the central compartment;
    # there is no depot. Dose events should target `central` with rate set in
    # the event table (2 h infusion in the BU4 arm, 3 h infusion in the BU1
    # arm).
    kel <- cl / vc
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
