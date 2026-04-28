Lin_2024_casirivimab <- function() {
  description <- "Two-compartment population PK model for casirivimab in pediatric and adult subjects (non-infected, ambulatory or hospitalized SARS-CoV-2-infected, or household contacts) following IV or SC administration (Lin 2024, casirivimab arm of the joint casirivimab + imdevimab popPK model)"
  reference <- "Lin K-J, Turner MA, Pasoll D, et al. Population Pharmacokinetics of Casirivimab and Imdevimab in Pediatric and Adult Non-Infected Individuals, Pediatric and Adult Ambulatory or Hospitalized Patients or Household Contacts of Patients Infected with SARS-COV-2. Pharmaceutical Research. 2024;41(10):1933-1949. doi:10.1007/s11095-024-03764-5"
  vignette <- "Lin_2024_casirivimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline; allometric scaling on CL and Vc with reference weight 81.6 kg (population median, Lin 2024 Table 1). Pediatrics < 6 years use fixed allometric exponents 0.75 (CL) and 1.0 (Vc); subjects >= 6 years use the estimated exponents 0.7959 (CL) and 0.5392 (Vc).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline; power scaling on CL with reference age 45 years. Also used to derive the pediatric indicator (AGE < 6) that switches the allometric exponents on CL/Vc and the SC bioavailability term.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on CL and Vc relative to the male reference.",
      source_name        = "SEXF"
    ),
    RACE_WHITE = list(
      description        = "White race indicator (1 = White, 0 = non-White)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-White, pooling Black/African American, Asian, American Indian/Alaska Native, Native Hawaiian/Pacific Islander, Other, Not reported, Unknown per Table 1 of Lin 2024)",
      notes              = "Multiplicative effect on CL relative to the non-White reference. Renamed from source column RACE to canonical RACE_WHITE per inst/references/covariate-columns.md.",
      source_name        = "RACE"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; power scaling on CL and Vc with reference 43 g/L (Lin 2024 Table 1 population median).",
      source_name        = "ALB"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment per NCI ODWG criteria (1 = mild, 0 = others)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function or pooled moderate/severe; reference complement per Lin 2024)",
      notes              = "Multiplicative effect on CL. Renamed from source column HEPIMP to canonical HEPIMP_MILD per inst/references/covariate-columns.md.",
      source_name        = "HEPIMP"
    ),
    SARS_VLOAD = list(
      description        = "SARS-CoV-2 baseline viral load (RT-qPCR, nasopharyngeal swab)",
      units              = "log10 copies/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 6.4 log10 copies/mL (Lin 2024 Table 1 population median in COVID-positive subjects). Non-infected subjects encoded as 0 in the source dataset.",
      source_name        = "VIRAL"
    ),
    SARS_SEROPOS = list(
      description        = "SARS-CoV-2 baseline antibody serostatus positive (1 = positive, 0 = negative or other/unknown)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (seronegative; 'Other' / unknown pooled into reference per Lin 2024 analysis plan)",
      notes              = "Multiplicative effect on CL. Renamed from source column SERPOS to canonical SARS_SEROPOS per inst/references/covariate-columns.md.",
      source_name        = "SERPOS"
    ),
    CRP = list(
      description        = "C-reactive protein (standard assay)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; power scaling on CL with reference 5.48 mg/L (Lin 2024 Table 1 pooled median across studies that collected CRP).",
      source_name        = "CRP"
    ),
    NLR = list(
      description        = "Neutrophil-to-lymphocyte ratio (CBC differential-derived)",
      units              = "ratio",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; power scaling on CL with reference 2.11 (Lin 2024 Table 1 pooled median).",
      source_name        = "NLR"
    ),
    OXYSUP_LOW = list(
      description        = "Low-flow supplemental oxygen at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no supplemental oxygen at baseline)",
      notes              = "Multiplicative effect on CL. Renamed from source column OXYSTAT1 to canonical OXYSUP_LOW per inst/references/covariate-columns.md.",
      source_name        = "OXYSTAT1"
    ),
    OXYSUP_HIGH = list(
      description        = "High-flow supplemental oxygen or mechanical ventilation at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no supplemental oxygen at baseline)",
      notes              = "Multiplicative effect on CL. In Lin 2024 the small mechanical-ventilation subset was pooled into the high-flow indicator. Renamed from source column OXYSTAT2 to canonical OXYSUP_HIGH per inst/references/covariate-columns.md.",
      source_name        = "OXYSTAT2"
    )
  )

  population <- list(
    n_subjects     = 7598,
    n_studies      = 7,
    age_range      = "0-98 years (median 45)",
    age_median     = 45,
    weight_range   = "8.6-235 kg (median 81.6)",
    weight_median  = 81.6,
    sex_female_pct = 50.1,
    race_ethnicity = "White 81.8%, Black/African American 7.3%, Asian 3.3%, American Indian/Alaska Native 1.2%, Native Hawaiian/Pacific Islander 0.2%, Other/Not reported 3.5%, Unknown 2.7%, Missing < 0.1% (Lin 2024 Table 1)",
    disease_state  = "Pediatric and adult non-infected individuals, ambulatory or hospitalized SARS-CoV-2-infected patients, and household contacts of infected patients (pooled).",
    dose_range     = "300-8000 mg IV single dose; 600-1200 mg SC single dose; or 1200 mg SC every 4 weeks (Lin 2024 Methods).",
    regions        = "Multinational (seven Phase 1/2/3 trials: NCT04426695, NCT04425629, NCT04452318, NCT04519437, NCT04666441, NCT05092581, NCT04992273)",
    serostatus     = "Positive 31.0%, Negative 62.8%, Other 6.1%, Missing 0.1%",
    oxygen_status  = "No supplemental oxygen 64.8%, Low-flow 23.8%, High-flow 1.4%, Mechanical ventilation 0.3%, Missing 9.7%",
    notes          = "Joint popPK model fits casirivimab + imdevimab simultaneously; this nlmixr2lib entry implements the casirivimab arm only (independent ODE chain in the source model). Imdevimab can be added in a parallel file. Pediatric coverage extends to infants (one subject < 1 year), but >= 98% of pediatrics are 2+ years old."
  )

  ini({
    # Structural parameters for casirivimab. Reference subject for the typical
    # values (Lin 2024 page 1939, paragraph after Eq. 4): 45-year-old non-White
    # male, 81.6 kg, ALB = 43 g/L, baseline viral load = 6.4 log10 copies/mL,
    # CRP = 5.48 mg/L, NLR = 2.11, seronegative, no supplemental oxygen.
    lcl     <- log(0.1926); label("Casirivimab clearance for the typical reference subject (CL, L/day)")  # Lin 2024 Table 2: theta1 = 0.1926
    lvc     <- log(3.917);  label("Casirivimab central volume of distribution at 81.6 kg, ALB 43 g/L, male (Vc, L)")  # Lin 2024 Table 2: theta2 = 3.917
    lq      <- log(0.4131); label("Casirivimab intercompartmental clearance (Q, L/day)")  # Lin 2024 Table 2
    lvp     <- log(3.065);  label("Casirivimab peripheral volume of distribution (Vp, L)")  # Lin 2024 Table 2
    lka     <- log(0.2183); label("Casirivimab first-order absorption rate constant (Ka, 1/day)")  # Lin 2024 Table 2
    lfdepot     <- log(0.7200);  label("Casirivimab subcutaneous bioavailability for adults / pediatrics >= 6 years (F1, fraction)")  # Lin 2024 Table 2
    lfdepot_ped <- log(0.8788);  label("Casirivimab subcutaneous bioavailability for pediatrics < 6 years (F1_ped, fraction)")  # Lin 2024 Table 2

    # Allometric exponents on body weight (reference 81.6 kg). For pediatrics
    # < 6 years the exponents are FIXED to classical values (0.75 on CL, 1.0
    # on Vc) per Lin 2024 page 1939 paragraph "The final PK model... had a
    # separate bioavailability term... and fixed exponents... to classical
    # allometric exponents (0.75 for CL and 1 for Vc) for children < 6 years
    # of age." For adults / pediatrics >= 6 years the exponents are estimated.
    e_wt_cl <- 0.7959;  label("Body-weight exponent on CL for adults and pediatrics >= 6 years (unitless)")  # Lin 2024 Table 2: theta13 (Weight on CL)
    e_wt_vc <- 0.5392;  label("Body-weight exponent on Vc for adults and pediatrics >= 6 years (unitless)")  # Lin 2024 Table 2: theta14 (Weight on Vc)

    # Other covariate effects on CL (reference values: AGE 45, ALB 43 g/L,
    # SARS_VLOAD 6.4 log10 copies/mL, CRP 5.48 mg/L, NLR 2.11; SEXF=0,
    # RACE_WHITE=0, HEPIMP_MILD=0, SARS_SEROPOS=0, OXYSUP_LOW=OXYSUP_HIGH=0).
    e_age_cl       <-  0.07037;  label("Power exponent: age effect on CL, (AGE/45)^e_age_cl (unitless)")  # Lin 2024 Table 2: theta15 (Age on CL)
    e_sexf_cl      <- -0.08051;  label("Casirivimab fractional change in CL for female vs male reference (unitless)")  # Lin 2024 Table 2: theta16 (Sex on CL, casirivimab)
    e_white_cl     <- -0.09478;  label("Fractional change in CL for White vs non-White reference (unitless)")  # Lin 2024 Table 2: theta17 (Race on CL, shared)
    e_alb_cl       <- -1.078;    label("Casirivimab power exponent: albumin on CL, (ALB/43)^e_alb_cl (unitless)")  # Lin 2024 Table 2: theta18 (Albumin on CL, casirivimab)
    e_hepimp_cl    <-  0.06602;  label("Fractional change in CL for mild hepatic impairment vs reference (unitless)")  # Lin 2024 Table 2: theta19 (Hepatic impairment on CL, shared)
    e_vload_cl     <- -0.00754;  label("Power exponent: SARS-CoV-2 viral load on CL, (SARS_VLOAD/6.4)^e_vload_cl (unitless)")  # Lin 2024 Table 2: theta20 (Viral load on CL, shared)
    e_seropos_cl   <-  0.07315;  label("Fractional change in CL for SARS-CoV-2 seropositive vs seronegative (unitless)")  # Lin 2024 Table 2: theta21 (Serostatus on CL, shared)
    e_crp_cl       <-  0.02252;  label("Power exponent: CRP on CL, (CRP/5.48)^e_crp_cl (unitless)")  # Lin 2024 Table 2: theta22 (CRP on CL, shared)
    e_nlr_cl       <-  0.02883;  label("Power exponent: NLR on CL, (NLR/2.11)^e_nlr_cl (unitless)")  # Lin 2024 Table 2: theta24 (NLR on CL, shared)
    e_oxylow_cl    <-  0.1064;   label("Fractional change in CL for low-flow supplemental oxygen vs none (unitless)")  # Lin 2024 Table 2: theta25 (Low oxygen supply on CL, shared)
    e_oxyhigh_cl   <-  0.3802;   label("Fractional change in CL for high-flow supplemental oxygen / ventilation vs none (unitless)")  # Lin 2024 Table 2: theta26 (High oxygen supply on CL, shared)

    # Covariate effects on Vc.
    e_sexf_vc      <- -0.1092;   label("Casirivimab fractional change in Vc for female vs male reference (unitless)")  # Lin 2024 Table 2: theta28 (Sex on Vc, casirivimab)
    e_alb_vc       <- -0.4167;   label("Power exponent: albumin on Vc, (ALB/43)^e_alb_vc (unitless)")  # Lin 2024 Table 2: theta30 (Albumin on Vc, shared)

    # Inter-individual variability (omega^2 = log(CV^2 + 1) for log-normal
    # random effects). The source paper (Lin 2024) shares the same eta values
    # for casirivimab and imdevimab on CL, Vc, and KA; this casirivimab-only
    # extraction carries the same variance estimates.
    etalcl ~ 0.08642  # 30.04 percent CV, Lin 2024 Table 2 IIV in CL; omega^2 = log(0.3004^2 + 1)
    etalvc ~ 0.11297  # 34.58 percent CV, Lin 2024 Table 2 IIV in Vc; omega^2 = log(0.3458^2 + 1)
    etalka ~ 0.47791  # 78.28 percent CV, Lin 2024 Table 2 IIV in KA narrative plus 95 percent CI midpoint 78.05 to 78.51; the printed table point estimate 72.28 lies outside its own CI and is treated as a typesetting typo (omega^2 = log(0.7828^2 + 1))

    # Residual variability. NONMEM "additive on log-transformed data" with
    # estimate 23.52 corresponds to a proportional error in linear space
    # with sigma ~= 0.2352 (Lin 2024 page 1939 narrative + Methods Eq. for
    # ln(Y_ij) = ln(C_ij) + e_ij).
    propSd <- 0.2352; label("Casirivimab proportional residual error (fraction)")  # Lin 2024 Table 2: residual variability 23.52
  })
  model({
    # Pediatric switch: subjects < 6 years use the fixed allometric exponents
    # 0.75 (CL) / 1.0 (Vc) and the higher SC bioavailability term, per Lin
    # 2024 page 1939 final-model description.
    PED <- (AGE < 6)
    allo_cl <- e_wt_cl * (1 - PED) + 0.75 * PED
    allo_vc <- e_wt_vc * (1 - PED) + 1.0  * PED
    fdepot  <- exp(lfdepot) * (1 - PED) + exp(lfdepot_ped) * PED

    # Multiplicative covariate effects on CL (Lin 2024 Eq. 1 for casirivimab).
    cl_cov <- (AGE / 45)^e_age_cl *
              (1 + e_sexf_cl    * SEXF) *
              (1 + e_white_cl   * RACE_WHITE) *
              (ALB / 43)^e_alb_cl *
              (1 + e_hepimp_cl  * HEPIMP_MILD) *
              (SARS_VLOAD / 6.4)^e_vload_cl *
              (1 + e_seropos_cl * SARS_SEROPOS) *
              (CRP / 5.48)^e_crp_cl *
              (NLR / 2.11)^e_nlr_cl *
              (1 + e_oxylow_cl  * OXYSUP_LOW) *
              (1 + e_oxyhigh_cl * OXYSUP_HIGH)

    # Multiplicative covariate effects on Vc (Lin 2024 Eq. 3 for casirivimab).
    vc_cov <- (1 + e_sexf_vc * SEXF) * (ALB / 43)^e_alb_vc

    cl <- exp(lcl + etalcl) * (WT / 81.6)^allo_cl * cl_cov
    vc <- exp(lvc + etalvc) * (WT / 81.6)^allo_vc * vc_cov
    q  <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Concentration: dose in mg, volume in L -> mg/L (the bioanalytical assay
    # in Lin 2024 reports mg/L; LLOQ 0.156 mg/L for both mAbs).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
