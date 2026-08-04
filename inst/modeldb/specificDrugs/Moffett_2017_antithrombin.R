Moffett_2017_antithrombin <- function() {
  description <- paste(
    "One-compartment population PK model with exponential residual error",
    "for intravenous human (plasma-derived) antithrombin (AT) concentrate",
    "in paediatric patients (Moffett 2017). Structural clearance uses an",
    "allometric weight arm with additive concurrent-UFH-dose potentiation;",
    "volume of distribution uses allometric weight scaled by a power of the",
    "per-subject baseline AT activity level. The observed AT activity level",
    "is the sum of dose-derived AT concentration plus the endogenous",
    "baseline (Dansirikul B3 handling of endogenous baseline)."
  )
  reference <- paste(
    "Moffett BS, Diaz R, Galati M, Mahoney D, Teruya J, Yee DL.",
    "Population pharmacokinetics of human antithrombin concentrate in",
    "paediatric patients.",
    "Br J Clin Pharmacol. 2017 Nov;83(11):2450-2456.",
    "doi:10.1111/bcp.13359.",
    sep = " "
  )
  vignette <- "Moffett_2017_antithrombin"
  units <- list(
    time          = "h",
    dosing        = "units",
    concentration = "units/dL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "antithrombin", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying in paediatric infants over the analysis window.",
        "Used for theoretical-allometric scaling with reference 70 kg and",
        "fixed exponents 0.75 on CL and 1.00 on VD (Moffett 2017 Table 3).",
        "Cohort median 5.2 kg (IQR 3.4-14.8 kg)."
      ),
      source_name        = "WT"
    ),
    DOSE_UFH_UH = list(
      description        = "Concurrent continuous-infusion unfractionated heparin dose rate",
      units              = "units/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Absolute infusion rate (NOT weight-normalized). Time-varying per",
        "observation as the UFH infusion is titrated. Cohort median",
        "absolute rate 173 units/h serves as the reference in the",
        "(DOSE_UFH_UH / 173) linear-additive covariate ratio on CL.",
        "Set to 0 units/h for observations when the patient is not",
        "receiving concurrent UFH; the additive covariate term then",
        "collapses to zero and CL reduces to the baseline weight-scaled",
        "arm. The cohort mean per-kg UFH dose was 34.1 +/- 22.7 units/kg/h",
        "(Moffett 2017 Results); the source column value is the absolute",
        "rate (= per-kg rate * body weight)."
      ),
      source_name        = "UFH"
    ),
    AT_BL_UDL = list(
      description        = "Per-subject baseline (predose) plasma antithrombin activity",
      units              = "units/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject: the single predose measurement before the",
        "first AT-concentrate dose. Cohort mean 59 +/- 17 units/dL",
        "(Moffett 2017 Table 2); reference 60 in the (AT_BL_UDL / 60)",
        "power-form covariate ratio on VD. Enters the model in two places:",
        "(1) as a power-form covariate on VD with exponent -0.389 per",
        "Moffett 2017 Table 3, and (2) as an additive endogenous-baseline",
        "term in the observation equation per the Dansirikul et al. 2008",
        "B3 method for endogenous baseline handling.",
        "The Stachrom chromogenic AT-activity assay reports the same",
        "numerical value in both 'units/dL' and '% of normal'; both are",
        "used interchangeably by the source paper. Normal healthy paediatric",
        "range 85-130 units/dL; the cohort baseline is depressed",
        "relative to normal because the enrolled patients had consumptive",
        "AT depletion secondary to critical illness / concurrent UFH."
      ),
      source_name        = "BASE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 184L,
    n_studies      = 1L,
    age_range      = "birth to <19 years",
    age_median     = "0.35 years (IQR 0.07-3.9)",
    weight_range   = "not reported by the paper (IQR 3.4-14.8 kg)",
    weight_median  = "5.2 kg (IQR 3.4-14.8)",
    sex_female_pct = 53.3,
    race_ethnicity = "not reported by the paper",
    disease_state  = paste(
      "Paediatric inpatients receiving human (plasma-derived) antithrombin",
      "concentrate (Thrombate, Grifols) for consumptive antithrombin",
      "depletion secondary to critical illness or concurrent",
      "unfractionated heparin therapy. 87.5% of the cohort received",
      "concurrent continuous-infusion UFH (mean 34.1 +/- 22.7 units/kg/h);",
      "8.2% for coagulopathy and 4.3% during enoxaparin therapy.",
      "Patients on mechanical circulatory support during AT administration",
      "or monitoring were excluded."
    ),
    dose_range     = paste(
      "Median 2 doses per patient (IQR 1-4); mean 46.3 +/- 13.6 units/kg",
      "per dose infused over 15 minutes. Median 3 postdose AT activity",
      "levels per patient (IQR 2-7) at a median 20.8 h (IQR 12.4-32.1)",
      "after a dose; mean postdose AT activity 85.9 +/- 20.7 units/dL."
    ),
    regions        = "single-centre, Texas Children's Hospital, Houston, TX, USA",
    age_groups     = paste(
      "neonates (<=30 d, n=50), infants (31 d - 2 y, n=85),",
      "children (3-12 y, n=30), adolescents (13-18 y, n=19)."
    ),
    validation_cohort = paste(
      "Prospective validation cohort of 30 patients meeting the same",
      "inclusion criteria during 1 January 2016 - 30 June 2016;",
      "median 1.75% relative prediction error (IQR -11.75 to 6.5),",
      "median 10.1% absolute prediction error (IQR 4.8-19.6) for the",
      "first postdose AT activity level."
    ),
    notes          = paste(
      "Baseline demographics from Moffett 2017 Table 1; baseline",
      "laboratory values from Table 2. NONMEM v.7.2 (Icon) with",
      "PDx-Pop 5.1 using FOCE-I; bootstrap 1000 replications (Table 4).",
      "AT activity measured on the STA-R analyser with the STA-Stachrom",
      "AT kit; the assay coefficient of variation was 9.35%. The paper",
      "identifies serum albumin, age, PCA, PLT, HGB, HCT, and SCR as",
      "covariates tested but not retained after backward elimination",
      "(Appendix); only UFH on CL and BASE on VD survived the p < 0.001",
      "retention threshold."
    )
  )

  ini({
    # Structural intercepts (reference: WT = 70 kg, DOSE_UFH_UH = 0, AT_BL_UDL = 60 units/dL).
    # From Moffett 2017 Table 3 and Table 4 (final-model point estimates).
    lcl <- log(0.917); label("Reference clearance intercept (dL/h) for WT=70 kg, UFH=0")   # Moffett 2017 Table 3 & Table 4: CL_ref = 0.917 dL/h (RSE 16.8%)
    lvc <- log(67.9);  label("Reference volume of distribution (dL) for WT=70 kg, BASE=60") # Moffett 2017 Table 3 & Table 4: VD_ref = 67.9 dL (RSE 4.9%)

    # Allometric exponents were fixed at theoretical values in the source paper
    # (Methods "Model building": exponents of 0.75 on CL and 1 on VD used for
    # the primary model to which covariates were then applied). Table 3 reports
    # no uncertainty for these exponents, consistent with them being fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless; theoretical)")     # Moffett 2017 Methods "Model building" & Table 3 (fixed at theoretical value)
    e_wt_vc <- fixed(1.00); label("Allometric exponent on VD (unitless; theoretical)")     # Moffett 2017 Methods "Model building" & Table 3 (fixed at theoretical value)

    # Covariate effects: UFH on CL is linear-additive (dL/h per unit of UFH/173 ratio);
    # baseline AT on VD is a power form with negative exponent.
    e_ufh_cl   <-  0.129;  label("Linear-additive UFH effect on CL (dL/h per (UFH/173) unit)") # Moffett 2017 Table 3 & Table 4: +0.129 dL/h per (UFH/173) (RSE 15.5%)
    e_at_bl_vc <- -0.389;  label("Power exponent of baseline AT activity on VD (unitless)")     # Moffett 2017 Table 3 & Table 4: -0.389 (RSE 37.5%)

    # Inter-individual variability (exponential IIV, "modelled exponentially" per
    # Methods "Model building"). Table 4 reports CV% for each eta; the underlying
    # variance is omega^2 = log(CV^2 + 1).
    #   omega^2(CL) = log(0.475^2 + 1) = log(1.2256) = 0.20353
    #   omega^2(VD) = log(0.282^2 + 1) = log(1.07952) = 0.07649
    # Off-diagonal correlation not reported; treat as independent block.
    etalcl ~ 0.20353    # Moffett 2017 Table 4: omega1 = 47.5% CV (RSE 26.5%)  ->  variance = log(0.475^2 + 1)
    etalvc ~ 0.07649    # Moffett 2017 Table 4: omega2 = 28.2% CV (RSE 31.2%)  ->  variance = log(0.282^2 + 1)

    # Residual error: "exponential error model best fit the data" (Methods
    # "Pharmacokinetic modelling"). The paper's residual variability is
    # reported as CV% (11.1%); for the log-normal / exponential residual
    # `Cc ~ lnorm(expSd)`, the SD on the log scale is
    # expSd = sqrt(log(CV^2 + 1)) = sqrt(log(1.01232)) = 0.11056.
    expSd <- 0.11056; label("Exponential (log-normal) residual SD (unitless log-scale)")   # Moffett 2017 Table 4: exponential error 11.1% CV (RSE 26.9%)  ->  sqrt(log(0.111^2 + 1))
  })
  model({
    # Individual PK parameters. CL uses the paper's additive-covariate structure
    # (Moffett 2017 Table 3): TVCL = CL_ref * (WT/70)^0.75 + slope_UFH * (UFH/173),
    # then CL_i = TVCL_i * exp(etalcl_i). VD is the paper's multiplicative-covariate
    # structure: TVVD = VD_ref * (WT/70)^1 * (AT_BL_UDL/60)^-0.389, VD_i = TVVD * exp(etalvc_i).
    tvcl <- exp(lcl) * (WT / 70)^e_wt_cl + e_ufh_cl * (DOSE_UFH_UH / 173)
    tvvc <- exp(lvc) * (WT / 70)^e_wt_vc * (AT_BL_UDL / 60)^e_at_bl_vc

    cl <- tvcl * exp(etalcl)
    vc <- tvvc * exp(etalvc)

    # First-order elimination one-compartment model. The paper uses a one-
    # compartment structural model (Results: "A one-compartment exponential
    # error model best fit the data").
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Observation. Following the Dansirikul et al. 2008 B3 method (Moffett 2017
    # Methods "Data handling": "Baseline activity levels ... were included in
    # the error model to account for the endogenous production of AT"), the
    # measured plasma AT activity is the sum of dose-derived antithrombin
    # (central / vc) and the endogenous baseline (AT_BL_UDL). The residual
    # error is log-normal (exponential-error in NONMEM's parlance).
    #
    # Dimensional check: `central` is in units (mass), `vc` in dL, so
    # `central / vc` is in units/dL, matching AT_BL_UDL's units.
    Cc <- central / vc + AT_BL_UDL
    Cc ~ lnorm(expSd)
  })
}
