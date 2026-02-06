Clegg_2024_nirsevimab <- function() {
  description <- "Population pharmacokinetics of nirsevimab in Preterm and Term Infants"
  reference <- "Clegg L, Freshwater E, Leach A, Villafana T, Wählby Hamrén U. Population Pharmacokinetics of Nirsevimab in Preterm and Term Infants. J Clin Pharmacol. 2024;64(5):555-567. doi:10.1002/jcph.2401"
  notes <-
    c(
      "depot is intramuscular",
      "Race covariates can be dervied with Clegg_2024_nirsevimab_derive_race_indicators()"
    )
  units <-
    list(
      dosing = "mg",
      concentration = "µg/mL",
      time = "day"
    )
  covariateData <-
    list(
      WT = list(
        description = "Body weight",
        units       = "kg",
        type        = "continuous",
        notes       = "Can be time-varying. Used for allometric scaling of CL, Q, V2, V3."
      ),
      PAGE = list(
        description = "Postmenstrual age at time of dosing",
        units       = "months",
        type        = "continuous",
        derivation  = "gestational_age_weeks / 4.35 + postnatal_age_months",
        notes       = "Drives the CL maturation function."
      ),
      RACE_BLK_OTH = list(
        description = "Race indicator: Black or African American, or Other",
        type        = "binary",
        coding      = list(`0` = "Reference (White, Native Hawaiian/Pacific Islander, Asian, Am.Ind./Alaskan Native, Multiple)",
                           `1` = "Black or African American, or Other"),
        affects     = "CL"
      ),
      RACE_ASN_AMI_MUL = list(
        description = "Race indicator: Asian, American Indian or Alaskan Native, or Multiple races",
        type        = "binary",
        coding      = list(`0` = "Reference (White, Native Hawaiian/Pacific Islander, Black/Afr.Am., Other)",
                           `1` = "Asian, American Indian or Alaskan Native, or Multiple races"),
        affects     = c("CL", "V2")
      ),
      SEASON2 = list(
        description = "RSV season indicator",
        type        = "binary",
        coding      = list(`0` = "First RSV season",
                            `1` = "Second RSV season (200 mg dose)"),
        affects     = "CL"
      ),
      ADACAT = list(
        description = "Antidrug antibody status (categorical)",
        type        = "binary",
        coding      = list(`0` = "ADA negative",
                            `1` = "ADA positive"),
        affects     = "CL"
      )
    )
  ini({

    # ── Structural parameters (70 kg adult reference) ──────────────────────
    lcl  <- log(38.8)   ; label("Clearance (mL/day)")
    lv2  <- log(1980)   ; label("Central volume of distribution (mL)")
    lq   <- log(709)    ; label("Intercompartmental clearance (mL/day)")
    lv3  <- log(2400)   ; label("Peripheral volume of distribution (mL)")
    lka  <- log(0.401)  ; label("Absorption rate constant (1/day)")
    tf1  <- 0.839       ; label("Bioavailability (unitless)")

    # ── Covariate parameters ──────────────────────────────────────────────
    # Maturation (PAGE effect on CL)
    beta_cl <- 0.364    ; label("Fractional CL at term birth vs maturation (unitless)")
    t50_cl  <- 14.8     ; label("Maturation half-life (months)")

    # Allometric scaling exponents
    allo_cl <- 0.589    ; label("Allometric exponent on CL and Q (unitless)")
    allo_v  <- 0.84     ; label("Allometric exponent on V2 and V3 (unitless)")

    # Race effects – one-hot encoded binary indicators
    #   Reference = White / Native Hawaiian or Pacific Islander
    race_cl_blk_oth     <-  0.132  ; label("CL race effect: Black/Afr.Am. or Other (unitless)")
    race_cl_asn_ami_mul <- -0.0894 ; label("CL race effect: Asian/Am.Ind./Multiple (unitless)")
    race_v2_asn_ami_mul <- -0.226  ; label("V2 race effect: Asian/Am.Ind./Multiple (unitless)")

    # Season 2 effect on CL
    season2_cl <- -0.122 ; label("Season 2 effect on CL (unitless)")

    # ADA categorical effect on CL (positive vs negative)
    ada_cl <- 0.124      ; label("ADA-positive effect on CL (unitless)")

    # ── Between-subject variability ───────────────────────────────────────
    # CL–V2 BLOCK(2) with r = 0.785 (Table 2)
    # OMEGA values from NONMEM final estimates:
    bsv_cl + bsv_v2 ~ c(
      0.0709,
      0.0789, 0.1625
    ) ; label("IIV on CL", "IIV on V2")
    bsv_ka ~ 0.2138      ; label("IIV on Ka")

    # ── Residual unexplained variability ──────────────────────────────────
    propSd <- 0.21     ; label("Proportional residual error (unitless)")
  })

  model({
    # ── Maturation function on CL ─────────────────────────────────────────
    clmatf <- 1 - (1 - beta_cl) * exp(-(PAGE - (40 / 4.35)) * log(2) / t50_cl)

    # ── Allometric scaling ────────────────────────────────────────────────
    alcl <- (WT / 70)^allo_cl
    alv  <- (WT / 70)^allo_v

    # ── Race covariate factors (one-hot encoded) ──────────────────────────
    cl_race <- 1 + race_cl_blk_oth * RACE_BLK_OTH +
                   race_cl_asn_ami_mul * RACE_ASN_AMI_MUL

    v2_race <- 1 + race_v2_asn_ami_mul * RACE_ASN_AMI_MUL

    # ── Season 2 effect ───────────────────────────────────────────────────
    cl_season <- 1 + season2_cl * SEASON2

    # ── ADA effect ────────────────────────────────────────────────────────
    cl_ada <- 1 + ada_cl * ADACAT

    # ── Individual PK parameters ──────────────────────────────────────────
    cl <- exp(lcl + bsv_cl) * alcl * clmatf * cl_race * cl_season * cl_ada
    v2 <- exp(lv2 + bsv_v2) * alv  * v2_race
    q  <- exp(lq)           * alcl
    v3 <- exp(lv3)          * alv
    ka <- exp(lka + bsv_ka)

    # ── Micro-rate constants (2-cpt with first-order absorption) ──────────
    k   <- cl / v2
    k23 <- q  / v2
    k32 <- q  / v3

    # ── Differential equations ────────────────────────────────────────────
    d/dt(depot)  = -ka * depot
    d/dt(central)  =  ka * depot - k * central - k23 * central + k32 * periph
    d/dt(periph) =  k23 * central - k32 * periph

    # ── Bioavailability for depot compartment ─────────────────────────────
    f(depot) <- tf1

    # ── Concentration and error model ─────────────────────────────────────
    # Dose is in mg, V in mL → central/v2 is mg/mL.
    # Multiply by 1000 to convert to µg/mL (matching DV units).
    Cc <- 1000 * central / v2
    Cc ~ prop(propSd)
  })
}
