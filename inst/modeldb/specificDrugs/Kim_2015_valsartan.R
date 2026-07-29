Kim_2015_valsartan <- function() {
  description <- "Two-compartment population PK model for valsartan with zero-order absorption in healthy adult Korean male volunteers (Kim 2015)"
  reference <- "Kim Y, Son H, Son M, Lee D, Heo YA, Park K. Assessment of statistical power for covariate effects in data from phase I clinical trials. Transl Clin Pharmacol. 2015;23(1):31-34. doi:10.12793/tcp.2015.23.1.31"
  vignette <- "Kim_2015_valsartan"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL, V1, Q, V2 with reference weight 70 kg. Source column WT in kg (Kim 2015 Eq. 1).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR in raw mL/min (NOT BSA-normalized). Stored under canonical CRCL per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 125.5 mL/min (Kim 2015 cohort mean; Table 1). Power-law effect on CL: (CRCL / 125.5)^e_crcl_cl.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48,
    n_studies      = 1,
    age_range      = "20-45 years",
    age_median     = "26.8 years (mean, SD 6)",
    weight_range   = "53.5-85 kg",
    weight_median  = "68.7 kg (mean, SD 7.7)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy volunteers",
    dose_range     = "single oral dose (FDC tablet of amlodipine/valsartan; specific valsartan strength reported in upstream Kim 2013 source data, Clin Ther 35:934-940)",
    regions        = "South Korea",
    notes          = "Underlying concentration data taken from Kim et al. 2013 (Clin Ther 35:934-940, doi:10.1016/j.clinthera.2013.05.021), a fixed-dose-combination bioequivalence study of amlodipine and valsartan in healthy male Korean volunteers. Demographic distributions reported in Kim 2015 Table 1. The 2-compartment popPK model with zero-order absorption was newly developed for Kim 2015 from those data."
  )

  ini({
    # Structural parameters -- reference values for a 70 kg adult with CRCL = 125.5 mL/min
    # Kim 2015 Eq. 1 (page 32): CL/F, V1/F, Q/F, V2/F are apparent oral parameters;
    # bioavailability F is folded into all CL/V terms (not estimated separately).
    lcl  <- log(6.18);  label("Apparent oral clearance for a 70 kg adult with CRCL = 125.5 mL/min (CL/F, L/h)")   # Kim 2015 Methods, Eq. 1 + parameter list page 32 (THETA(1))
    lvc  <- log(25.9);  label("Apparent central volume of distribution for a 70 kg adult (V1/F, L)")              # Kim 2015 Methods, Eq. 1 + parameter list page 32 (THETA(2))
    ld1  <- log(4.39);  label("Zero-order absorption duration into central compartment (D1, h)")                  # Kim 2015 Methods, Eq. 1 + parameter list page 32 (THETA(3))
    lq   <- log(2.01);  label("Apparent inter-compartmental clearance for a 70 kg adult (Q/F, L/h)")              # Kim 2015 Methods, Eq. 1 + parameter list page 32 (THETA(4))
    lvp  <- log(17.4);  label("Apparent peripheral volume of distribution for a 70 kg adult (V2/F, L)")           # Kim 2015 Methods, Eq. 1 + parameter list page 32 (THETA(5))

    # Allometric exponents -- fixed at canonical values in Kim 2015 Eq. 1 (no uncertainty reported)
    e_wt_cl_q  <- fixed(0.75);  label("Allometric exponent on CL and Q (unitless)")  # Kim 2015 Eq. 1, page 32 (literal **0.75)
    e_wt_vc_vp <- fixed(1);     label("Allometric exponent on V1 and V2 (unitless)") # Kim 2015 Eq. 1, page 32 (linear WT/70)

    # Covariate effect: CRCL power-law on CL
    e_crcl_cl <- 0.793; label("CRCL power-law exponent on CL (unitless)")  # Kim 2015 Methods page 32 (THETA(6) estimated value, P = 0.0048)

    # IIV on CL and V1 with off-diagonal covariance (Kim 2015 page 32:
    # omega^2_11 = 0.139 (CV 37.3%), omega^2_22 = 0.269 (CV 51.9%),
    # omega^2_12 = 0.094 (rho_12 = 0.49)).
    etalcl + etalvc ~ c(0.139,
                        0.094, 0.269)

    # Residual error -- NOT REPORTED in Kim 2015. Sigma was used in the
    # paper's SSE simulation but the numeric value is omitted from the
    # publication. Approximated at 20% proportional (sidecar option A,
    # 1501 task report): plausible default for a healthy-volunteer
    # Phase I small-molecule oral popPK; see vignette Assumptions and
    # deviations. Not from the paper.
    propSd <- 0.2; label("Proportional residual error (fraction) -- approximated; not reported in Kim 2015")
  })
  model({
    # Individual PK parameters with allometric weight scaling (reference 70 kg)
    # and CRCL power-law on CL (reference 125.5 mL/min).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * (CRCL / 125.5)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    d1 <- exp(ld1)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE with zero-order input directly into central
    # (NONMEM ADVAN3-style: dose targets the central compartment with
    # a zero-order infusion rate of amt/d1 over duration d1).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Zero-order absorption: duration d1 applies to doses with cmt = central.
    dur(central) <- d1

    # Concentration: dose in mg, volume in L -> mg/L; multiply by 1000 to report ng/mL
    # to match Kim 2015 Figure 1 axis units.
    Cc <- (central / vc) * 1000
    Cc ~ prop(propSd)
  })
}
