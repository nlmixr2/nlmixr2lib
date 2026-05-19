Naik_2016_vortioxetine <- function() {
  description <- "Two-compartment population PK model for vortioxetine in adult patients with major depressive disorder or generalized anxiety disorder, with first-order oral absorption, region-specific oral clearance, and linear creatinine-clearance and height effects on CL/F (Naik 2016)"
  reference <- "Naik H, Chan S, Vakilynejad M, Chen G, Loft H, Mahableshwarkar AR, Areberg J. A Population Pharmacokinetic-Pharmacodynamic Meta-Analysis of Vortioxetine in Patients with Major Depressive Disorder. Basic Clin Pharmacol Toxicol. 2016;118(5):344-355. doi:10.1111/bcpt.12513"
  vignette <- "Naik_2016_vortioxetine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear additive effect: CL/F increases by 0.18 L/hr per (CRCL - 106) mL/min. Naik 2016 used a raw mL/min value (not BSA-normalized); the paper does not state the formula explicitly. Reference 106 mL/min is the population median (Table 2).",
      source_name        = "CrCL"
    ),
    HT = list(
      description        = "Body height at baseline",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear additive effect: CL/F increases by 0.40 L/hr per (HT - 167) cm. Reference 167 cm is the population median (Table 2). Height was retained over weight and BMI in stepwise selection because it produced the larger reduction in CL IIV.",
      source_name        = "HT"
    ),
    REGION_EUROPE = list(
      description        = "Indicator for study site in the European Union",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = USA (USA is the typical-value reference; combine with REGION_ROW = 0 for USA)",
      notes              = "Naik 2016 estimated three region-specific typical CL/F values: USA = 51 L/hr (reference), EU = 39 L/hr, RoW = 38 L/hr. EU CL/F was about 23 percent lower than USA. Encoded here as a log-multiplicative shift e_region_europe_cl = log(39/51).",
      source_name        = "REGION_EU"
    ),
    REGION_ROW = list(
      description        = "Indicator for study site in the Rest of World (sites in Canada, Australia, and Asia for this paper)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = USA (USA is the typical-value reference; combine with REGION_EUROPE = 0 for USA)",
      notes              = "Naik 2016 RoW comprised study sites in Canada, Australia, and Asia (Table 1). RoW CL/F was about 25 percent lower than USA. Encoded here as a log-multiplicative shift e_region_row_cl = log(38/51).",
      source_name        = "REGION_RoW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3160,
    n_studies      = 12,
    age_range      = "18-88 years",
    age_median     = "46 years",
    weight_range   = "39-173 kg",
    weight_median  = "74 kg",
    height_range   = "137-203 cm",
    height_median  = "167 cm",
    sex_female_pct = 67,
    crcl_range     = "26-322 mL/min",
    crcl_median    = "106 mL/min",
    disease_state  = "Major depressive disorder (10 studies) or generalized anxiety disorder (2 studies)",
    dose_range     = "1-20 mg once-daily oral",
    regions        = "USA, European Union, and Rest of World (Canada, Australia, Asia)",
    notes          = "Pooled Phase II/III patient PK data from 12 trials (Table 1 of Naik 2016). A total of 10498 plasma vortioxetine concentrations (excluding 8 percent LLOQ samples) from 3160 patients were used for population PK. Baseline demographics in Table 2; PK/efficacy analysis subset comprised 2537 patients with end-of-treatment MADRS data (sex 1709F : 828M, USA:non-USA = 765:1772)."
  )

  ini({
    # Structural PK parameters - Table 3 of Naik 2016
    # Reference subject for covariates: USA, height = 167 cm, CRCL = 106 mL/min
    lcl       <- log(51);   label("Typical CL/F in USA reference, log-scale (L/hr)")            # Table 3 ("CL/F for US" = 51 L/hr)
    lvc       <- log(2900); label("Central volume of distribution V2/F, log-scale (L)")          # Table 3 (V2/F = 2.9 x 10^3 L)

    # Parameters fixed in Naik 2016 from the upstream Phase I popPK (Areberg et al. 2014).
    # The Phase I model is not yet packaged in nlmixr2lib; the fixed values are listed
    # explicitly here so this file is self-contained.
    lq        <- fixed(log(23));   label("Intercompartmental clearance Q/F, log-scale (L/hr)")    # Table 3 (Q/F = 23 L/hr, fixed)
    lvp       <- fixed(log(670));  label("Peripheral volume of distribution V3/F, log-scale (L)") # Table 3 (V3/F = 6.7 x 10^2 L, fixed)
    lka       <- fixed(log(0.14)); label("First-order absorption rate constant, log-scale (1/hr)")# Table 3 (ka = 0.14 /hr, fixed; Table 3 unit "L/hr" is a typo)

    # Region effects on CL/F (log-multiplicative form recovering the paper's
    # additive intercepts TVCL_USA = 51, TVCL_EU = 39, TVCL_RoW = 38 L/hr).
    e_region_europe_cl <- log(39 / 51); label("Log multiplicative effect of EU region on CL/F (TVCL_EU / TVCL_USA = 39 / 51)")  # Table 3
    e_region_row_cl    <- log(38 / 51); label("Log multiplicative effect of RoW region on CL/F (TVCL_RoW / TVCL_USA = 38 / 51)") # Table 3

    # Covariate effects on CL/F (linear additive on L/hr, per equation 12 of Naik 2016)
    e_crcl_cl <- 0.18; label("Linear effect of creatinine clearance on CL/F (L/hr per (CRCL - 106) mL/min)") # Table 3 (CrCL on CL/F)
    e_ht_cl   <- 0.40; label("Linear effect of height on CL/F (L/hr per (HT - 167) cm)")                    # Table 3 (Height on CL/F)

    # IIV: Naik 2016 estimated separate variances per region (EU:USA:RoW = 0.38:0.90:0.62).
    # This file uses a single etalcl with the RoW value (0.62) as a pragmatic single-value
    # compromise; the original region-specific values are documented in the vignette
    # under "Assumptions and deviations".
    etalcl    ~ 0.62  # Table 3 (omega^2 for CL/F, RoW value used as single representative)
    etalvc    ~ 0.82  # Table 3 (omega^2 for V2/F)

    # Residual error: Naik 2016 used an additive error model for log-transformed
    # concentrations (NONMEM LTBS). expSd with ~ lnorm() is the matching nlmixr2
    # parameterization.
    expSd     <- 0.26; label("Residual standard deviation on log-transformed concentrations")  # Table 4 (Residual error)
  })

  model({
    # Region-specific typical CL/F via log-multiplicative shifts. The exp() form recovers
    # the paper's additive intercepts: USA = 51, EU = 39, RoW = 38 L/hr.
    tvcl_region <- exp(lcl + e_region_europe_cl * REGION_EUROPE + e_region_row_cl * REGION_ROW)

    # Equation 12: CL/F = TVCL_region + 0.18 * (CRCL - 106) + 0.40 * (HT - 167)
    # IIV applied multiplicatively: CL_i = CL_typical * exp(etalcl)
    cl <- (tvcl_region + e_crcl_cl * (CRCL - 106) + e_ht_cl * (HT - 167)) * exp(etalcl)

    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Convert mg / L to ng/mL (factor 1000) for paper-native units
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
