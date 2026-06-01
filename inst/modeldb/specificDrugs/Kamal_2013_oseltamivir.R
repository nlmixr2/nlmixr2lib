Kamal_2013_oseltamivir <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral oseltamivir",
    "(prodrug, OP) and its active metabolite oseltamivir carboxylate",
    "(OC) in 390 subjects aged 1 to 78 years pooled from 13 clinical",
    "trials (healthy adults, influenza-inoculated and naturally infected",
    "adults, healthy geriatric subjects, renally impaired adults, and",
    "healthy and infected pediatric subjects 1 to 18 years). Oseltamivir",
    "is described by a two-compartment model with first-order absorption",
    "and first-order conversion to OC (CLp/F treated as the OP-to-OC",
    "conversion clearance under the assumption of complete metabolism;",
    "<5% of prodrug is excreted unchanged renally). OC is described by",
    "a one-compartment model with first-order elimination. All clearance",
    "and volume terms are apparent (conditioned on oral bioavailability",
    "F; OC terms additionally on the fraction metabolized fm, assumed 1).",
    "Covariates: body weight as a power function on OP CLp/F, OC CLm/F,",
    "and OC Vcm/F (allometric-style exponents estimated, not fixed);",
    "creatinine clearance (BSA-normalized to 1.73 m^2) as a power function",
    "on OC CLm/F; and age as a linear (additive) term on OC Vcm/F.",
    "Inter-individual variability is exponential on all seven structural",
    "parameters, with two off-diagonal covariances (CLp/F with CLm/F,",
    "and Vp/F with Vcm/F). Residual error is proportional only for",
    "oseltamivir (40.5% CV reduced CCV model) and combined additive",
    "plus proportional for OC (14.0% CV proportional + 17.9 ng/mL",
    "additive SD).")
  reference <- "Kamal MA, Van Wart SA, Rayner CR, Subramoney V, Reynolds DK, Bulik CC, Smith PF, Bhavnani SM, Ambrose PG, Forrest A. Population pharmacokinetics of oseltamivir: pediatrics through geriatrics. Antimicrob Agents Chemother. 2013;57(8):3470-3477. doi:10.1128/AAC.02438-12"
  vignette <- "Kamal_2013_oseltamivir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on OP CLp/F (exponent 0.838), OC CLm/F (exponent 0.560), and OC Vcm/F (exponent 0.830), all centered at the 70 kg adult reference per Kamal 2013 Methods page 3471 ('Covariate evaluation'). The 70 kg reference was chosen by the authors for convenient comparison to a typical adult rather than the cohort median; cohort median was 64.5 kg and range 8-115 kg (Kamal 2013 Results page 3471). Weight is also discussed under the allometric approach in Kamal 2013 Discussion page 3473 with the estimated exponents close to but distinct from the canonical 0.75 / 1 allometric values, which the authors note as a feature of their empirical (data-driven) covariate search rather than a fixed-exponent approach.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (baseline).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (additive) effect on OC Vcm/F: Vcm/F = 238*(WT/70)^0.830 - 2.25*(AGE - 21). The age effect is subtracted after the WT power term, NOT multiplied as a percentage (Kamal 2013 Table 2 final-estimate row 'Vc m /F (liters) = 238(WT/70) 0.830 - 2.25(age - 21)'). The 21-year reference is the cohort median age; the cohort age range was 1 to 78 years (Kamal 2013 Results page 3471). The age slope was a relatively minor covariate (18-unit MVOF decrease) and was the smallest of the four retained covariate effects (Kamal 2013 Results page 3473). Outside the fitted age x weight range the additive form can yield negative Vcm/F; users simulating very low-weight elderly subjects should sanity-check the typical-value prediction.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Body-surface-area-normalized creatinine clearance.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on OC CLm/F (exponent 0.487) centered at the 95 mL/min/1.73 m^2 cohort median per Kamal 2013 Methods page 3471 and Table 2 final-estimate row 'CL m /F (liters/h) = 20.7(WT/70) 0.560 (CL CR /95) 0.487'. In adults >=18 years the source paper computed CrCl via Cockcroft-Gault using ideal body weight when actual weight exceeded IBW, normalized to a BSA of 1.73 m^2; SCr was capped at a lower bound of 0.7 mg/dL. In pediatric subjects 1-17 years the revised Schwartz equation was used (CRCL = 0.413*HTCM/SCr) with SCr capped at a lower bound of 0.2 mg/dL. Both forms yield mL/min/1.73 m^2 in this column. Note: the source column name CLCR is the same canonical concept as Delattre 2010 amikacin's raw-mL/min CLCR alias, but Kamal 2013 normalizes to BSA -- hence mapped to canonical CRCL (BSA-normalized) per inst/references/covariate-columns.md rather than the raw CLCR alias.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 390L,
    n_studies      = 13L,
    age_range      = "1-78 years",
    age_median     = "21 years",
    weight_range   = "8-115 kg",
    weight_median  = "64.5 kg",
    sex_female_pct = 38.2,
    race_ethnicity = "Not reported in the model parameter table; race was screened as a covariate but not retained in the final model.",
    disease_state  = "Pooled cohort of healthy and naturally / experimentally influenza-infected subjects (no PK difference between infected and non-infected was detected); also includes adults with mild / moderate / severe renal impairment (CRCL down to ~14 mL/min/1.73 m^2; only 1 subject had severe impairment, CRCL <30).",
    dose_range     = "Oral oseltamivir 20-1,000 mg as single dose or repeated-dose; standard influenza treatment regimen 75 mg twice daily (q12h) for 5 days; pediatric weight-based dosing 2 mg/kg BID for ages 1-12 yr and weight-band fixed doses (30 mg for 1-2 yr, 45 mg for 3-5 yr).",
    regions        = "Not specified in the source paper.",
    crcl_range     = "13.9-178 mL/min/1.73 m^2, median 95.1",
    crcl_categories = c(
      "normal (>=80)"         = 297L,
      "mild (50-80)"          = 73L,
      "moderate (30-49)"      = 19L,
      "severe (<30)"          = 1L
    ),
    n_observations = "3,881 oseltamivir concentrations + 4,402 OC concentrations; 10 of 13 studies used intensive sampling (>=5 samples/subject), 3 studies used sparse sampling (2-3 samples/subject).",
    notes          = "Cohort statistics from Kamal 2013 Results page 3471 (241 male / 149 female = 38.2% female). 13 studies tabulated in Kamal 2013 Table 1: WP15517, WP15525, PV15616, NP15717, WV15670, WV15671, WV15730, WP15647, WP15648, WV15758, NP15826, JV16284, PP16351. NONMEM 7 level 1.2 (ADVAN13) with FOCE-with-interaction estimation was used (Kamal 2013 Methods 'Structural model development', page 3471). Excluded from the source dataset: drug-drug interaction studies, end-stage renal disease patients on intermittent hemodialysis, and neonates / infants <1 year (a later cohort not yet available at the time of the original analysis). LLOQ values were 1 ng/mL (OP) and 8.8 ng/mL (OC); only 3 OP records below LLOQ were excluded, so no Beal M3 likelihood-based handling was needed."
  )

  ini({
    # Structural fixed effects (parent oseltamivir, OP) -- Kamal 2013 Table 2
    # final-estimate column. All clearance / volume terms are apparent
    # (conditioned on the fraction of OP absorbed F).

    lka  <- log(0.775); label("OP first-order absorption rate constant (ka, 1/h)")               # Kamal 2013 Table 2: ka = 0.775 1/h, %SEM 3.74
    lcl  <- log(519);   label("OP apparent clearance CLp/F at 70 kg adult (CLp/F, L/h)")          # Kamal 2013 Table 2: CLp/F coefficient = 519 L/h, %SEM 3.99
    lvc  <- log(421);   label("OP apparent central volume Vcp/F (Vcp/F, L)")                       # Kamal 2013 Table 2: Vcp/F = 421 L, %SEM 6.05
    lq   <- log(120);   label("OP apparent inter-compartmental clearance CLd/F (CLd/F, L/h)")     # Kamal 2013 Table 2: CLd/F = 120 L/h, %SEM 4.95
    lvp  <- log(2800);  label("OP apparent peripheral volume Vp/F (Vp/F, L)")                      # Kamal 2013 Table 2: Vp/F = 2,800 L, %SEM 7.46

    # Structural fixed effects (metabolite oseltamivir carboxylate, OC) -- Kamal 2013 Table 2.
    # Apparent clearance / volume conditioned on F (parent) and fm (assumed 1).

    lcl_oc <- log(20.7);  label("OC apparent clearance CLm/F at 70 kg adult, CRCL 95 (CLm/F, L/h)")    # Kamal 2013 Table 2: CLm/F coefficient = 20.7 L/h, %SEM 3.36
    lvc_oc <- log(238);   label("OC apparent central volume Vcm/F at 70 kg adult, age 21 (Vcm/F, L)")  # Kamal 2013 Table 2: Vcm/F coefficient = 238 L, %SEM 5.16

    # Covariate effects on apparent OP CLp/F.

    e_wt_cl <- 0.838; label("Power exponent of (WT/70) on OP CLp/F (unitless)")  # Kamal 2013 Table 2: power of WT on CLp/F = 0.838, %SEM 4.67

    # Covariate effects on apparent OC CLm/F.

    e_wt_cl_oc   <- 0.560; label("Power exponent of (WT/70) on OC CLm/F (unitless)")     # Kamal 2013 Table 2: power of WT on CLm/F = 0.560, %SEM 6.63
    e_crcl_cl_oc <- 0.487; label("Power exponent of (CRCL/95) on OC CLm/F (unitless)")    # Kamal 2013 Table 2: power of CLCR on CLm/F = 0.487, %SEM 9.01

    # Covariate effects on apparent OC Vcm/F.
    # The Vcm/F equation is ADDITIVE in age (paper text and Table 2):
    #   Vcm/F = 238*(WT/70)^0.830 - 2.25*(AGE - 21)
    # The age slope is in liters per year and is applied AFTER the WT
    # power scaling, not multiplied as a percentage. See vignette
    # Assumptions / deviations for caution about negative Vcm/F at
    # extreme low-weight elderly covariate combinations.

    e_wt_vc_oc  <- 0.830; label("Power exponent of (WT/70) on OC Vcm/F (unitless)")              # Kamal 2013 Table 2: power of WT on Vcm/F = 0.830, %SEM 18.7
    e_age_vc_oc <- -2.25; label("Linear additive slope of (AGE - 21) on OC Vcm/F (L/year)")       # Kamal 2013 Table 2: slope of age on Vcm/F = -2.25 L/year, %SEM 31.2

    # Inter-individual variability. Paper reports diagonal omega^2 as
    # lognormal % CV; the canonical conversion omega^2 = log(CV^2 + 1)
    # is used here. Two off-diagonal covariances are reported on the
    # eta (log) scale: omega(etalcl, etalcl_oc) = 0.0987 and
    # omega(etalvp, etalvc_oc) = 0.218 (Kamal 2013 Table 2 'Covariance'
    # rows). The reported r^2 values (0.372 and 0.274) are computed in
    # the paper using the small-omega approximation r ~= omega_xy /
    # (CV_x * CV_y) and therefore differ slightly from the implied
    # log-scale correlation; see the vignette Errata.

    etalcl + etalcl_oc ~ c(0.16317,
                           0.09870, 0.13688)                                            # Kamal 2013 Table 2: omega^2(CLp/F) = log(0.421^2+1) = 0.16317; omega(CLp,CLm) = 0.0987; omega^2(CLm/F) = log(0.383^2+1) = 0.13688
    etalvp + etalvc_oc ~ c(0.34064,
                           0.21800, 0.35515)                                            # Kamal 2013 Table 2: omega^2(Vp/F) = log(0.637^2+1) = 0.34064; omega(Vp,Vcm) = 0.218; omega^2(Vcm/F) = log(0.653^2+1) = 0.35515
    etalka ~ 0.08999                                                                     # Kamal 2013 Table 2: omega^2(ka) = log(0.307^2+1) = 0.08999 (30.7% CV)
    etalvc ~ 0.39231                                                                     # Kamal 2013 Table 2: omega^2(Vcp/F) = log(0.693^2+1) = 0.39231 (69.3% CV)
    etalq  ~ 0.32607                                                                     # Kamal 2013 Table 2: omega^2(CLd/F) = log(0.621^2+1) = 0.32607 (62.1% CV)

    # Residual error.
    # OP: reduced CCV (proportional-only), sigma_CCV = 40.5% CV.
    # OC: additive plus CCV (proportional + additive), sigma_CCV = 14.0%
    # CV proportional and sigma_ADD = 17.9 ng/mL = 0.0179 mg/L additive
    # SD on the linear concentration scale. The paper labels the table
    # rows 'sigma^2 CCV' and 'sigma^2 ADD' but the magnitudes are
    # reported as SDs (CV% for CCV; ng/mL for ADD); see the vignette
    # Errata for the small notation slip.

    propSd    <- 0.405;  label("OP proportional residual SD on linear concentration (fraction)")  # Kamal 2013 Table 2: sigma_CCV for oseltamivir = 40.5% CV
    propSd_oc <- 0.140;  label("OC proportional residual SD on linear concentration (fraction)")  # Kamal 2013 Table 2: sigma_CCV for OC = 14.0% CV
    addSd_oc  <- 0.0179; label("OC additive residual SD on linear concentration (mg/L)")          # Kamal 2013 Table 2: sigma_ADD for OC = 17.9 ng/mL = 0.0179 mg/L
  })

  model({
    # Reference covariate values for the WT power, CRCL power, and AGE
    # additive equations (Kamal 2013 Methods page 3471, Table 2 footer).
    ref_wt   <- 70
    ref_crcl <- 95
    ref_age  <- 21

    # Individual parameters (OP, parent).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Individual parameters (OC, metabolite). Vcm/F is additive in age
    # AFTER the WT power scaling: Vcm/F = 238*(WT/70)^0.830 -
    # 2.25*(AGE - 21). IIV is applied via exp(etalvc_oc) on the entire
    # covariate-adjusted typical value.
    cl_oc <- exp(lcl_oc + etalcl_oc) *
             (WT / ref_wt)^e_wt_cl_oc *
             (CRCL / ref_crcl)^e_crcl_cl_oc
    vc_oc <- (exp(lvc_oc) * (WT / ref_wt)^e_wt_vc_oc +
              e_age_vc_oc * (AGE - ref_age)) * exp(etalvc_oc)

    # Micro-constants. kel (= cl/vc) is the OP-to-OC first-order
    # conversion rate constant; OP elimination as unchanged prodrug
    # (< 5% of dose) is folded into kel under the paper's complete-
    # metabolism assumption (Kamal 2013 Discussion page 3473).
    kel    <- cl    / vc
    k12    <- q     / vc
    k21    <- q     / vp
    kel_oc <- cl_oc / vc_oc

    # ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_oc)  <-  kel * central - kel_oc * central_oc

    # Observations. Dose in mg / volume in L => concentration in mg/L.
    Cc    <- central    / vc
    Cc_oc <- central_oc / vc_oc

    Cc    ~ prop(propSd)
    Cc_oc ~ add(addSd_oc) + prop(propSd_oc)
  })
}
