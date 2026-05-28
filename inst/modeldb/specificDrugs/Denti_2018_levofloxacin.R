Denti_2018_levofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral levofloxacin in South",
    "African children with multidrug-resistant tuberculosis (MDR-TB)",
    "disease or exposure (Denti 2018; n = 109; median age 2.1 yr;",
    "median weight 12.4 kg). First-order absorption with an absorption",
    "lag time, allometric scaling fixed to 0.75 on CL / Q and 1 on Vc /",
    "Vp with the population-median 12 kg as the reference weight, and a",
    "Hill-type maturation function on CL driven by postmenstrual age",
    "(PMAGE_50 = 10.6 mo, gamma = 3.39; PMAGE = postnatal age + 9 mo",
    "assuming term gestation). Covariate effects: HIV-positive children",
    "have 15.9% lower CL; nasogastric-tube (NGT) administration shortens",
    "the absorption lag time by 85.6% relative to the oral reference.",
    "F is fixed at 1; the additive residual error is fixed at 20% of",
    "the LLOQ (0.0160 mg/L)."
  )
  reference <- paste(
    "Denti P, Garcia-Prats AJ, Draper HR, Wiesner L, Winckler J, Thee S,",
    "Dooley KE, Savic RM, McIlleron HM, Schaaf HS, Hesseling AC (2018).",
    "Levofloxacin population pharmacokinetics in South African children",
    "treated for multidrug-resistant tuberculosis.",
    "Antimicrobial Agents and Chemotherapy 62(2):e01521-17.",
    "doi:10.1128/AAC.01521-17.",
    sep = " "
  )
  vignette <- "Denti_2018_levofloxacin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the allometric scaling of CL and Q",
        "(exponent 0.75 fixed) and Vc and Vp (exponent 1 fixed) with",
        "12 kg as the reference weight, matching the typical-value",
        "convention used in Denti 2018 Table 2 footnote 'a' ('The",
        "typical values reported here refer to a 12-kg child aged 2",
        "years.')."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the Hill-type maturation function on",
        "CL with PMAGE_50 = 10.6 months and shape parameter gamma =",
        "3.39 (Denti 2018 Methods 'Population pharmacokinetic model",
        "development' final paragraph; Table 2; Figure 1). PMAGE is",
        "computed as postnatal age (months) + 9 months assumed",
        "gestation, since gestational age at birth was not recorded",
        "(Methods 'Population pharmacokinetic model development'",
        "final paragraph). The reference postmenstrual age for the",
        "Table 2 typical values is 33 months (2 yr postnatal + 9 mo",
        "assumed term gestation; maturation 97.9% complete per",
        "Table 2 footnote a)."
      ),
      source_name        = "PMAGE"
    ),
    HIV_POS = list(
      description        = "HIV-1 infection status indicator (1 = HIV-positive on ART, 0 = HIV-negative)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative; structural reference for the Table 2 typical CL of 4.70 L/h)",
      notes              = paste(
        "Time-fixed per subject. Multiplicative effect on CL:",
        "(1 + e_hiv_pos_cl * HIV_POS) with e_hiv_pos_cl = -0.159 per",
        "Denti 2018 Table 2 'HIV+ on CL (%)'. All 16 HIV-positive",
        "children in the cohort were on antiretroviral therapy",
        "(13 lopinavir-ritonavir, 3 efavirenz) per Results 'Study",
        "population and pharmacokinetic samples'; the paper did not",
        "have power to ascribe the effect to a particular ART",
        "regimen."
      ),
      source_name        = "HIV"
    ),
    ROUTE_NGT = list(
      description        = "Nasogastric-tube administration indicator (1 = NGT, 0 = oral)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral administration; structural reference for the Table 2 T_lag of 0.242 h)",
      notes              = paste(
        "Per-dose-record indicator. Multiplicative effect on T_lag:",
        "(1 + e_route_ngt_tlag * ROUTE_NGT) with e_route_ngt_tlag =",
        "-0.856 per Denti 2018 Table 2 'NGT on T_lag (%)'. NGT",
        "delivery shortens the absorption lag time to roughly 14.4%",
        "of its oral value. In the cohort 90 / 109 (82.6%) of",
        "children were dosed by crushed tablet via NGT, 12 / 109",
        "(11.0%) by crushed tablet orally, and 7 / 109 (6.4%) by",
        "whole tablet orally (Table 1)."
      ),
      source_name        = "NGT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 109L,
    n_studies      = 1L,
    age_range      = "0.32 to 8.65 years",
    age_median     = "2.1 years",
    weight_range   = "5.88 to 21.8 kg",
    weight_median  = "12.4 kg",
    sex_female_pct = 100 * (1 - 56 / 109),
    race_ethnicity = c(Black = 63.3, MixedRace = 36.7),
    disease_state  = paste(
      "Children routinely treated with levofloxacin for confirmed,",
      "probable, or possible multidrug-resistant tuberculosis (MDR-TB)",
      "disease (71 / 109; 65.1%) or preventive therapy following",
      "exposure to an infectious MDR-TB source case (38 / 109; 34.9%)",
      "in Cape Town, South Africa. 16 / 109 (14.7%) were HIV-1 infected",
      "and on antiretroviral therapy (lopinavir-ritonavir n = 13 or",
      "efavirenz n = 3)."
    ),
    dose_range     = paste(
      "Once-daily oral levofloxacin 10 to 15 mg/kg (2012 to 2013) or",
      "15 to 20 mg/kg (2013 to 2017) using the 250 mg adult tablet",
      "(Austel, South Africa). The dose was administered as a whole",
      "tablet swallowed (7 / 109), a crushed tablet swallowed orally",
      "(12 / 109), or a crushed tablet dispersed in water and",
      "administered via a nasogastric tube (90 / 109). Median",
      "individual dose 212 mg (range 88.5 to 435), median 15 mg/kg",
      "(range 10 to 21.4)."
    ),
    sampling       = paste(
      "Pre-dose plus 1, 2, 4, 6, and 8 hours post-dose. 662",
      "quantifiable levofloxacin concentrations were analysed (3",
      "participants contributed data from more than one sampling",
      "occasion). 36 / 662 (5.4%) below-LLOQ samples were all in the",
      "pre-dose record and were imputed at LLOQ / 2 per Beal's M6",
      "method."
    ),
    regions        = "South Africa (Cape Town).",
    notes          = paste(
      "Baseline demographics from Denti 2018 Table 1 and Results 'Study",
      "population and pharmacokinetic samples'. Gestational age at",
      "birth was not recorded; postmenstrual age was assumed to be",
      "postnatal age plus 9 months for the maturation function. The",
      "additive residual error was fixed at 20% of the LLOQ (0.0160",
      "mg/L) because the estimate hit the stipulated lower boundary",
      "(Table 2 footnote c). The 4.48-fold scaling of the F",
      "between-occasion variability for unobserved (predose) doses",
      "(Table 2) is a fitting nuisance specific to the predose records",
      "in this cohort and is not encoded for forward simulation."
    )
  )

  ini({
    # Structural parameters. Reference weight 12 kg and reference postmenstrual age
    # 33 months (= 2 yr postnatal + 9 mo assumed term gestation): Denti 2018 Table 2
    # footnote 'a' typical values refer to a 12-kg child aged 2 years.

    lcl     <- log(4.70);  label("Clearance (CL, L/h) at 12 kg / 33 mo PMAGE / HIV-")  # Table 2 row "CL"
    lvc     <- log(19.2);  label("Central volume of distribution (Vc, L) at 12 kg")    # Table 2 row "Vc"
    lq      <- log(0.796); label("Intercompartmental clearance (Q, L/h) at 12 kg")     # Table 2 row "Q"
    lvp     <- log(3.40);  label("Peripheral volume of distribution (Vp, L) at 12 kg") # Table 2 row "Vp"
    lka     <- log(1.61);  label("Absorption rate constant (ka, 1/h)")                 # Table 2 row "ka"
    ltlag   <- log(0.242); label("Absorption lag time (T_lag, h) for oral dosing")     # Table 2 row "T_lag" (for oral)
    lfdepot <- fixed(log(1)); label("Bioavailability anchor (F, fraction)")            # Table 2 row "F" = 1 (fixed)

    # Allometric exponents (fixed per Denti 2018 Methods: 0.75 on CL parameters and
    # 1 on volume parameters, citing Anderson and Holford 2008).
    allo_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Methods, p.7 col.2 paragraph 2
    allo_v  <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)") # Methods, p.7 col.2 paragraph 2

    # Maturation parameters (Hill function on postmenstrual age, applied to CL).
    pmage50   <- 10.6; label("PMAGE_50: postmenstrual age (months) at 50% maturation") # Table 2 row "PMAGE_50"
    gamma_mat <- 3.39; label("Maturation shape parameter (unitless)")                  # Table 2 row "gamma"

    # Covariate effects.
    e_hiv_pos_cl     <- -0.159; label("HIV+ multiplicative effect on CL (fraction)")        # Table 2 row "HIV+ on CL (%)"
    e_route_ngt_tlag <- -0.856; label("NGT multiplicative effect on T_lag (fraction)")      # Table 2 row "NGT on T_lag (%)"

    # Inter-individual variability. Denti 2018 Table 2 reports BSV on CL only and
    # BOV on T_lag, ka, and F. nlmixr2lib has no idiomatic encoding for BOV
    # separate from BSV (see precedent in Bienczak_2016_nevirapine.R citing
    # Svensson_2018_bedaquiline.R / Svensson_2016_rifampicin.R); the BOV terms are
    # folded in here as BSV-equivalent etas. Internal variance scale is
    # omega^2 = log(1 + CV^2).
    etalcl     ~ 0.02284  # BSV CL  = 15.2% CV  -> log(1 + 0.152^2) (Table 2; eta shrinkage 24%)
    etaltlag   ~ 0.98954  # BOV Tlag = 130% CV  -> log(1 + 1.30^2)  (Table 2; eta shrinkage 76%, folded as BSV-equivalent)
    etalka     ~ 0.35068  # BOV ka   = 64.8% CV -> log(1 + 0.648^2) (Table 2; eta shrinkage 43%, folded as BSV-equivalent)
    etalfdepot ~ 0.04641  # BOV F    = 21.8% CV -> log(1 + 0.218^2) (Table 2; eta shrinkage 10%, folded as BSV-equivalent)

    # Residual error: combined proportional + additive (additive fixed at 20% of LLOQ).
    propSd <- 0.116;           label("Proportional residual error (fraction)")  # Table 2 row "Proportional error (%)"
    addSd  <- fixed(0.0160);   label("Additive residual error (mg/L; fixed at 20% of LLOQ)") # Table 2 row "Additive error" (fixed)
  })

  model({
    # Maturation function on CL: Hill-type, normalized to the Table 2 reference
    # postmenstrual age (33 months) so that lcl directly equals the displayed
    # typical CL of 4.70 L/h at 12 kg, 2 yr, HIV-.
    pmage_ref      <- 33
    maturation_pma <- PAGE^gamma_mat        / (pmage50^gamma_mat + PAGE^gamma_mat)
    maturation_ref <- pmage_ref^gamma_mat   / (pmage50^gamma_mat + pmage_ref^gamma_mat)
    mat_cl         <- maturation_pma / maturation_ref

    # HIV multiplicative effect on CL.
    hiv_cl <- 1 + e_hiv_pos_cl * HIV_POS

    # NGT multiplicative effect on absorption lag time.
    ngt_tlag <- 1 + e_route_ngt_tlag * ROUTE_NGT

    # Individual PK parameters with allometric weight scaling (12 kg reference).
    cl <- exp(lcl + etalcl)         * (WT / 12)^allo_cl * mat_cl * hiv_cl
    vc <- exp(lvc)                  * (WT / 12)^allo_v
    q  <- exp(lq)                   * (WT / 12)^allo_cl
    vp <- exp(lvp)                  * (WT / 12)^allo_v
    ka <- exp(lka + etalka)
    tlag_central <- exp(ltlag + etaltlag) * ngt_tlag

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)   <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag_central

    # Concentration in mg/L (= ug/mL). Combined proportional + additive residual.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
