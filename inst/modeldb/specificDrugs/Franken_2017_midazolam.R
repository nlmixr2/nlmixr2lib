Franken_2017_midazolam <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for midazolam, its primary",
    "active metabolite 1-OH-midazolam (1-OH-M), and the secondary metabolite",
    "1-OH-midazolam-glucuronide (1-OH-MG) in 45 terminally ill adult",
    "palliative-care patients (Franken 2017). Midazolam: one-compartment",
    "disposition with two parallel first-order absorption routes (oral and",
    "subcutaneous bolus) using route-specific absorption rate constants",
    "fixed from literature (Ka oral = 5.5 1/h, Ka SC = 10 1/h); oral",
    "bioavailability F is estimated and SC F assumed = 1. 1-OH-M: one",
    "compartment, central volume fixed equal to midazolam V, clearance",
    "estimated. 1-OH-MG: one compartment, clearance and volume estimated.",
    "All inter-compartment fluxes carry the parent / metabolite signal in",
    "midazolam-equivalent mass units (concentrations were adjusted to",
    "midazolam equivalents via molecular weight per Methods). Midazolam",
    "clearance depends on serum albumin (power form, reference 25 g/L) and",
    "1-OH-MG clearance depends on eGFR (standard four-variable MDRD, power",
    "form, reference 104 mL/min/1.73 m^2). IIV on midazolam CL was",
    "correlated with oral F (rho fixed to unity per Results); other IIVs",
    "are independent. Residual variability is additive on log-transformed",
    "concentrations (LTBS) for all three analytes; a cross-output residual",
    "correlation noted in Methods is not encoded in this nlmixr2 port",
    "(see vignette Assumptions and deviations)."
  )
  reference <- paste(
    "Franken LG, Masman AD, de Winter BCM, Baar FPM, Tibboel D,",
    "van Gelder T, Koch BCP, Mathot RAA.",
    "Hypoalbuminaemia and decreased midazolam clearance in terminally",
    "ill adult patients, an inflammatory effect?",
    "Br J Clin Pharmacol. 2017;83(8):1701-1712.",
    "doi:10.1111/bcp.13259.",
    sep = " "
  )
  vignette <- "Franken_2017_midazolam"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying serum albumin. Used with power scaling",
        "(ALB / 25)^1.08 on midazolam clearance, normalised to the",
        "population median of 25 g/L (Table 1; range 13-39). Same",
        "covariate would not improve volume of distribution (delta-OFV",
        "= 0.014 per Discussion); only clearance was retained. The",
        "paper's Discussion attributes the albumin-clearance correlation",
        "to an underlying inflammatory / catabolic state rather than to",
        "protein binding (CRP showed a comparable effect; albumin won",
        "selection because it is more clinically stable). Units: g/L",
        "(multiply g/dL values by 10)."
      ),
      source_name        = "albumin (g/L)"
    ),
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate computed from the",
        "four-variable Modification of Diet in Renal Disease (MDRD)",
        "formula (age, sex, race, serum creatinine; the standard",
        "abbreviated MDRD equation)."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying renal-function covariate. Used with power scaling",
        "(CRCL / 104)^0.53 on 1-OH-MG clearance, normalised to the",
        "population median of 104 mL/min/1.73 m^2 (Table 1; standard",
        "four-variable MDRD; range 6-328). The paper also tested the",
        "original six-variable MDRD equation but retained the simpler",
        "four-variable formula in the final model because the difference",
        "was minimal. eGFR did not affect midazolam or 1-OH-M clearance;",
        "only 1-OH-MG (the renally cleared end metabolite) retained the",
        "covariate. Source column 'eGFR' in the paper text and table;",
        "canonical register name CRCL covers MDRD-derived eGFR."
      ),
      source_name        = "eGFR (standard four-variable MDRD)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45,
    n_studies      = 1,
    age_range      = "43-93 years",
    age_median     = "71 years",
    weight_range   = "not reported (weight available for ~53% of subjects, baseline-only)",
    weight_median  = "not reported",
    sex_female_pct = 51.1,
    race_ethnicity = c(Caucasian = 91.1, AfroCaribbean = 6.7, Unknown = 2.2),
    disease_state  = paste(
      "Terminally ill adult palliative-care patients admitted to a",
      "Dutch hospice with prognosis of 2 days to 3 months and",
      "advanced malignancy as primary diagnosis in 97.8% (one",
      "respiratory-disease patient). Sparse venous samples (median 2",
      "per subject, range 1-10) drawn in both pre-terminal and",
      "terminal phases (terminal = bed-bound, semi-comatose, unable to",
      "take oral medication). Median duration of admittance 29 days",
      "(range 7-457)."
    ),
    dose_range     = paste(
      "Midazolam 7.5 mg oral up to four times daily, or subcutaneous",
      "bolus / infusion 2.5-180 mg/day. Administered per Dutch",
      "national palliative guidelines for insomnia or palliative",
      "sedation. Intravenous route not used (typical of palliative",
      "care)."
    ),
    regions        = "The Netherlands (Laurens Cadenza palliative care centre, Rotterdam)",
    notes          = paste(
      "Demographics from Franken 2017 Table 1. NONMEM 7.2 + PsN 4.4.8",
      "with FOCE-I; ADVAN7 subroutine; data were log-transformed",
      "before fitting. 192 blood samples collected for parent and",
      "metabolites (139 quantifiable; BQL fractions: 14% midazolam,",
      "16% 1-OH-M, 10% 1-OH-MG). 4-OH-M (75% BQL) was not modelled.",
      "M1 (discard) BQL handling. Baseline blood chemistry medians:",
      "albumin 25 g/L (range 13-39), urea 7.6 mmol/L, bilirubin 9",
      "umol/L, GGT 62 U/L, ALP 118 U/L, ALT 14 U/L, AST 30 U/L, CRP",
      "92 mg/L, creatinine 67 umol/L, eGFR (standard four-variable",
      "MDRD) 104 mL/min/1.73 m^2 (range 6-328). 37.8% used",
      "dexamethasone; 2.2% used phenytoin (concomitant CYP inducers).",
      "Final model evaluated by 200-run bootstrap and normalised",
      "prediction distribution errors (NPDE; global adjusted p-values",
      "0.75, 0.20, 0.41 for midazolam, 1-OH-M, 1-OH-MG)."
    )
  )

  ini({
    # ================================================================
    # Midazolam absorption (route-specific Ka, both fixed from
    # literature per Methods 'Structural model': 'the absorption
    # constants (Ka) could not be estimated. They were therefore
    # derived from literature (5.5 h-1 for oral administration, 10
    # h-1 for subcutaneous injection)').
    # ================================================================
    lka_oral <- fixed(log(5.5))
    label("Ka oral route (1/h, FIXED literature)")              # Methods / Table 2 footnote: Ka oral = 5.5 1/h fixed (literature)
    lka_sc   <- fixed(log(10))
    label("Ka subcutaneous route (1/h, FIXED literature)")      # Methods / Table 2 footnote: Ka SC = 10 1/h fixed (literature)

    # ================================================================
    # Oral bioavailability - estimated. SC F is structurally fixed
    # at 1 per Methods 'Base model development' ('Bioavailability of
    # subcutaneous midazolam was assumed to be 100%').
    # ================================================================
    lfdepot <- log(0.279)
    label("Oral bioavailability F (log scale)")                 # Table 2 Final model: F = 0.279 (RSE 12.6%)

    # ================================================================
    # Midazolam structural disposition - Franken 2017 Table 2 Final
    # model. Reference values are at the population median ALB =
    # 25 g/L (so typical midazolam CL = 8.42 at the reference subject).
    # ================================================================
    lcl <- log(8.42)
    label("Midazolam clearance at population median albumin (L/h)")  # Table 2 Final: CL midazolam = 8.42 L/h (RSE 9.0%)
    lvc <- log(113)
    label("Midazolam central volume V (L)")                          # Table 2 Final: V midazolam = 113 L (RSE 13.1%); also used for 1-OH-M (V_1OH-M assumed equal to V midazolam per Methods)

    # ================================================================
    # 1-OH-M metabolite disposition - one-compartment, V assumed
    # equal to midazolam V (Methods 'Base model development': 'The
    # volume of distribution (V) of 1-OH-M was assumed to be equal
    # to the volume of distribution of midazolam.'). Only CL is
    # estimated for 1-OH-M.
    # ================================================================
    lcl_1ohm <- log(45.4)
    label("1-OH-midazolam clearance (L/h)")                          # Table 2 Final: CL 1-OH-M = 45.4 L/h (RSE 11.5%)

    # ================================================================
    # 1-OH-MG metabolite disposition - one-compartment, both CL and
    # V estimated. CL reference is at population median eGFR =
    # 104 mL/min/1.73 m^2.
    # ================================================================
    lcl_1ohmg <- log(5.10)
    label("1-OH-midazolam-glucuronide clearance at population median eGFR (L/h)")  # Table 2 Final: CL 1-OH-MG = 5.10 L/h (RSE 11.0%)
    lvc_1ohmg <- log(2.98)
    label("1-OH-midazolam-glucuronide central volume V (L)")         # Table 2 Final: V 1-OH-MG = 2.98 L (RSE 71.5%; high RSE noted in Discussion - patients without renal insufficiency eliminate 1-OH-MG more rapidly so V is hard to estimate)

    # ================================================================
    # Covariate effects (power form per Methods Eq. 1).
    #   midazolam CL = CL_pop * (ALB / 25)^e_alb_cl
    #   1-OH-MG CL  = CL_pop * (CRCL / 104)^e_crcl_cl_1ohmg
    # Verified against Discussion 'Simulations' section:
    #   ALB=35 -> CL = 8.42 * (35/25)^1.08 = 12.10 L/h (paper: "12.1")
    #   ALB=15 -> CL = 8.42 * (15/25)^1.08 =  4.87 L/h (paper: "4.8")
    #   eGFR=90  -> 1OH-MG CL = 5.10 * (90/104)^0.53 = 4.73 (paper: "4.7")
    #   eGFR=50  -> 1OH-MG CL = 5.10 * (50/104)^0.53 = 3.50 (paper: "3.5")
    # ================================================================
    e_alb_cl <- 1.08
    label("Albumin power exponent on midazolam CL (unitless)")       # Table 2 Final: Albumin covariate effect on midazolam clearance = 1.08 (RSE 21.2%)
    e_crcl_cl_1ohmg <- 0.53
    label("eGFR power exponent on 1-OH-MG CL (unitless)")            # Table 2 Final: eGFR covariate effect on 1-OH-MG clearance = 0.53 (RSE 20.7%)

    # ================================================================
    # Inter-individual variability.
    # IIV reported as CV% in Table 2 Final model column; converted to
    # the internal log-scale variance via omega^2 = log(1 + CV^2).
    #
    # 'Base model development': 'The correlation between IIV of
    # midazolam clearance and F of oral midazolam was high (0.93) and
    # therefore fixed to unity.' Encoded as a 2x2 BLOCK with the
    # off-diagonal equal to sqrt(var_F * var_CL), which forces rho = 1
    # (same mathematically as the published shared-eta-with-scaling
    # NONMEM construct).
    # ================================================================
    etalfdepot + etalcl ~ c(
      log(1 + 0.506^2),
      sqrt(log(1 + 0.506^2) * log(1 + 0.490^2)),
      log(1 + 0.490^2)
    )
    # Table 2 Final IIV: F = 50.6% CV (RSE 17.4%, shrinkage 12.8%),
    # midazolam CL = 49.0% CV (RSE 14.0%, shrinkage 12.8%); rho fixed to 1.

    etalvc ~ log(1 + 0.709^2)
    # Table 2 Final IIV: midazolam V = 70.9% CV (RSE 15.1%, shrinkage 16.6%)

    etalcl_1ohm ~ log(1 + 0.605^2)
    # Table 2 Final IIV: 1-OH-M CL = 60.5% CV (RSE 18.0%, shrinkage 12.2%)

    etalcl_1ohmg ~ log(1 + 0.499^2)
    # Table 2 Final IIV: 1-OH-MG CL = 49.9% CV (RSE 23.1%, shrinkage 23.0%)

    # ================================================================
    # Residual variability.
    # Methods 'Base model development': 'additive residual error on
    # logarithmic transformed concentrations' (LTBS pattern). Table 2
    # reports the residual variability as CV% for each analyte; the
    # log-scale SD supplied to lnorm() is sqrt(log(1 + CV^2)).
    #
    # NOTE: Methods also state 'a correlation between the residual
    # errors was incorporated in the model' (across midazolam / 1-OH-M
    # / 1-OH-MG because they share the same assay sample). nlmixr2 /
    # rxode2 do not support cross-output residual correlations
    # directly, so the residuals here are treated as independent.
    # See vignette Assumptions and deviations for the cited deviation.
    # ================================================================
    expSd       <- sqrt(log(1 + 0.268^2))
    label("Midazolam log-normal residual SD on log scale")           # Table 2 Final: midazolam residual = 26.8% CV (RSE 13.3%, shrinkage 20.4%)
    expSd_1ohm  <- sqrt(log(1 + 0.423^2))
    label("1-OH-M log-normal residual SD on log scale")              # Table 2 Final: 1-OH-M residual = 42.3% CV (RSE 21.6%, shrinkage 18.5%)
    expSd_1ohmg <- sqrt(log(1 + 0.464^2))
    label("1-OH-MG log-normal residual SD on log scale")             # Table 2 Final: 1-OH-MG residual = 46.4% CV (RSE 13.1%, shrinkage 18.6%)
  })

  model({
    # ------------------------------------------------------------
    # Reference covariate values (Table 1 population medians).
    # ------------------------------------------------------------
    alb_ref  <- 25    # Table 1: median albumin = 25 g/L
    crcl_ref <- 104   # Table 1: median standard MDRD eGFR = 104 mL/min/1.73 m^2

    # ------------------------------------------------------------
    # Route-specific absorption rate constants (both fixed).
    # ------------------------------------------------------------
    ka_oral <- exp(lka_oral)
    ka_sc   <- exp(lka_sc)

    # ------------------------------------------------------------
    # Oral bioavailability with IIV (correlated with midazolam CL
    # at rho = 1 via the OMEGA block above). SC bioavailability is
    # rxode2 default 1 (no f(depot) assignment for the SC depot).
    # ------------------------------------------------------------
    f_oral <- exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------
    # Individual structural parameters with covariate effects.
    # Midazolam CL: ALB power-form on the typical CL.
    # ------------------------------------------------------------
    alb_factor <- (ALB / alb_ref)^e_alb_cl
    cl <- exp(lcl + etalcl) * alb_factor
    vc <- exp(lvc + etalvc)

    # 1-OH-M: CL estimated; V structurally equal to midazolam V
    # (Methods 'Base model development'). The same IIV draw for
    # midazolam V is implicitly applied to 1-OH-M V via the shared
    # vc symbol.
    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm)
    vc_1ohm <- vc

    # 1-OH-MG: CL with eGFR power-form covariate; V independent.
    crcl_factor_1ohmg <- (CRCL / crcl_ref)^e_crcl_cl_1ohmg
    cl_1ohmg <- exp(lcl_1ohmg + etalcl_1ohmg) * crcl_factor_1ohmg
    vc_1ohmg <- exp(lvc_1ohmg)

    # ------------------------------------------------------------
    # ODE system. Two parallel depots feed the midazolam central
    # compartment (oral via depot, SC via depot2). The full
    # midazolam elimination flux (cl * Cc) becomes the 1-OH-M
    # formation flux (apparent values: per the paper, the true
    # fraction metabolised to 1-OH-M is unknown but treated as 1
    # given the data; concentrations were pre-adjusted to midazolam
    # equivalents using molecular weight, so 1 mg of midazolam
    # cleared corresponds to 1 mg-equivalent of 1-OH-M). The full
    # 1-OH-M elimination flux similarly becomes the 1-OH-MG
    # formation flux; 1-OH-MG is fully cleared (renally per
    # Discussion).
    # ------------------------------------------------------------
    d/dt(depot)         <- -ka_oral * depot
    d/dt(depot2)        <- -ka_sc   * depot2
    d/dt(central)       <-  ka_oral * depot +
                            ka_sc   * depot2 -
                            cl * central / vc
    d/dt(central_1ohm)  <-  cl * central / vc -
                            cl_1ohm * central_1ohm / vc_1ohm
    d/dt(central_1ohmg) <-  cl_1ohm * central_1ohm / vc_1ohm -
                            cl_1ohmg * central_1ohmg / vc_1ohmg

    # Oral bioavailability on the oral depot only. SC depot stays
    # at rxode2 default F = 1 (per paper assumption).
    f(depot) <- f_oral

    # ------------------------------------------------------------
    # Observations. Internal states are mg of midazolam-equivalents
    # (concentrations were adjusted to midazolam equivalents using
    # molecular weight). Plasma concentrations are reported in ug/L
    # (paper LLOQs: 4 ug/L midazolam, 2 ug/L 1-OH-M, 8 ug/L
    # 1-OH-MG). 1 mg/L = 1000 ug/L.
    # ------------------------------------------------------------
    Cc        <- central         / vc        * 1000
    Cc_1ohm   <- central_1ohm    / vc_1ohm   * 1000
    Cc_1ohmg  <- central_1ohmg   / vc_1ohmg  * 1000

    Cc        ~ lnorm(expSd)
    Cc_1ohm   ~ lnorm(expSd_1ohm)
    Cc_1ohmg  ~ lnorm(expSd_1ohmg)
  })
}
