Soto_2014_ampicillin_sulbactam <- function() {
  description <- paste(
    "Joint two-compartment population PK model for ampicillin and sulbactam",
    "in 47 Japanese adults with moderate or severe community-acquired",
    "pneumonia receiving 30-minute IV infusions of 3 g ampicillin/sulbactam",
    "(2:1) every 6 hours (Soto 2014). Both drugs are fitted simultaneously",
    "via the NONMEM L2 data item; a single common Cockcroft-Gault CLcr",
    "power effect (0.701) is applied to CL of both drugs, and CL random",
    "effects are correlated across drugs (rho = 0.858). Body weight is a",
    "fixed linear allometric scalar on peripheral volume V2 for both drugs.",
    "Ampicillin uses the unsuffixed canonical compartment / parameter set;",
    "sulbactam carries the sibling-drug suffix _sbt throughout."
  )
  reference <- paste(
    "Soto E, Shoji S, Muto C, Tomono Y, Marshall S.",
    "Population pharmacokinetics of ampicillin and sulbactam in patients",
    "with community-acquired pneumonia: evaluation of the impact of renal",
    "impairment. Br J Clin Pharmacol. 2014;77(3):509-521.",
    "doi:10.1111/bcp.12232.",
    sep = " "
  )
  vignette <- "Soto_2014_ampicillin_sulbactam"
  units    <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as an allometric scalar on peripheral volume V2 for both",
        "ampicillin and sulbactam (Soto 2014 final-model equation): V2_i(k) =",
        "theta_V2(k) * (WT / 51.2)^theta_BWT * exp(eta_V2_i(k)) with theta_BWT",
        "FIXED at 1.0 (Table 2 combined final-model column, '- Fix' entry).",
        "Reference 51.2 kg is the population median body weight (Table 1).",
        "The reduced-final stepwise covariate analysis identified WT on V2 as",
        "significant for ampicillin (forwards inclusion) but not in the",
        "backwards exclusion step for sulbactam; both drugs were nevertheless",
        "given a physiological linear allometric V2 scaling (exponent = 1.0",
        "fixed) for clinical plausibility (Soto 2014 Results, paragraph 3 of",
        "the Pharmacokinetic analyses subsection).",
        "Time-fixed per subject (Soto 2014 Table 1 reports baseline weight",
        "only)."
      ),
      source_name        = "BWT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source paper uses raw (NOT BSA-normalized) creatinine clearance",
        "estimated by the Cockcroft and Gault equation (Soto 2014 Methods,",
        "Population pharmacokinetic analysis subsection and the cited",
        "reference [15]). Drives a common power effect on CL for both",
        "ampicillin and sulbactam: CL_i(k) = theta_CL(k) * (CRCL / 71)^",
        "theta_CLcr * exp(eta_CL_i(k)) with the SAME theta_CLcr (= 0.701)",
        "shared across drugs (Soto 2014 Table 2 combined final-model column;",
        "shared because integrating into a single parameter increased OFV",
        "by only 0.054 points). Reference 71 mL/min is the population median",
        "CLcr (Table 1). The cohort spanned 34.6-176 mL/min; patients with",
        "CLcr < 30 mL/min were excluded by protocol. Time-fixed per subject",
        "in the analysis."
      ),
      source_name        = "CLcr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47,
    n_studies      = 1,
    age_range      = "28-85 years (median 67)",
    age_median     = "67 years",
    weight_range   = "31.3-78.7 kg (median 51.2)",
    weight_median  = "51.2 kg",
    sex_female_pct = 45,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Moderate or severe community-acquired pneumonia (CAP) requiring",
      "in-hospital antimicrobial treatment. Renal function ranged from",
      "normal (CLcr >= 90 mL/min, 36% of cohort) to moderate impairment",
      "(30 mL/min <= CLcr < 60, 43%); patients with severe renal",
      "impairment (CLcr < 30 mL/min) were excluded by protocol",
      "(Soto 2014 Methods and Patient characteristics)."
    ),
    dose_range     = paste(
      "Ampicillin/sulbactam 3 g (2:1; ampicillin 2 g + sulbactam 1 g)",
      "administered as a 30-minute intravenous infusion every 6 hours",
      "for 3 to 14 days depending on clinical condition (Soto 2014",
      "Methods, Clinical studies and assay methods)."
    ),
    regions        = "Japan (multicentre)",
    crcl_range     = "34.6-176 mL/min (median 71; cohort excludes CLcr < 30)",
    bmi_range      = "13.7-29.0 kg/m^2 (median 20.4)",
    notes          = paste(
      "47 Japanese inpatients with CAP from ClinicalTrials.gov NCT01189487",
      "(Soto 2014 Methods). 444 plasma samples (222 each for ampicillin",
      "and sulbactam) were available for the joint population PK fit.",
      "Elderly patients (>= 65 years) accounted for 57% of the cohort;",
      "14 patients (30%) had baseline weight <= 45 kg.",
      "Demographics are from Soto 2014 Table 1."
    )
  )

  ini({
    # =====================================================================
    # Ampicillin structural parameters (Soto 2014 Table 2, combined
    # final-model column). Two-compartment IV PK with no absorption phase;
    # parameters refer to the typical subject at CRCL = 71 mL/min and
    # WT = 51.2 kg.
    # =====================================================================
    lcl <- log(10.7); label("Ampicillin CL at CRCL = 71 mL/min (L/h)")        # Table 2 ampicillin combined final-model column: CL(1) = 10.7 L/h (RSE 3.39%)
    lvc <- log(9.97); label("Ampicillin central volume V1 (L)")               # Table 2 ampicillin combined final-model column: V1(1) = 9.97 L (RSE 6.07%)
    lq  <- log(4.14); label("Ampicillin inter-compartmental clearance Q (L/h)")   # Table 2 ampicillin combined final-model column: Q(1) = 4.14 L/h (RSE 21.8%)
    lvp <- log(4.48); label("Ampicillin peripheral volume V2 at WT = 51.2 kg (L)") # Table 2 ampicillin combined final-model column: V2(1) = 4.48 L (RSE 9.91%)

    # =====================================================================
    # Sulbactam structural parameters (Soto 2014 Table 2, combined
    # final-model column).
    # =====================================================================
    lcl_sbt <- log(10.4); label("Sulbactam CL at CRCL = 71 mL/min (L/h)")        # Table 2 sulbactam combined final-model column: CL(2) = 10.4 L/h (RSE 3.40%)
    lvc_sbt <- log(10.2); label("Sulbactam central volume V1 (L)")               # Table 2 sulbactam combined final-model column: V1(2) = 10.2 L (RSE 7.04%)
    lq_sbt  <- log(4.58); label("Sulbactam inter-compartmental clearance Q (L/h)")   # Table 2 sulbactam combined final-model column: Q(2) = 4.58 L/h (RSE 28.2%)
    lvp_sbt <- log(4.04); label("Sulbactam peripheral volume V2 at WT = 51.2 kg (L)") # Table 2 sulbactam combined final-model column: V2(2) = 4.04 L (RSE 12.1%)

    # =====================================================================
    # Shared covariate effects (Soto 2014 Table 2 combined final-model
    # column). theta_CLcr (= 0.701) is a SINGLE parameter applied to CL of
    # both ampicillin and sulbactam (combining the two drug-specific CLcr
    # exponents into one common value increased the OFV by only 0.054
    # points; Results, Pharmacokinetic analyses subsection paragraph 4).
    # theta_BWT on V2 is FIXED at 1.0 for both drugs as a physiological
    # linear allometric scaling on peripheral volume.
    # =====================================================================
    e_crcl_cl <- 0.701;     label("Shared CRCL power exponent on CL for ampicillin and sulbactam (unitless)") # Table 2: theta CLcr on CL = 0.701 (RSE 7.65%), common to both drugs in the combined final model
    e_wt_vp   <- fixed(1.0); label("Body-weight allometric exponent on peripheral V2 for both drugs (FIXED at 1.0, unitless)") # Table 2: theta BWT on V2 = 1.00 Fix (combined final model, shared across drugs)

    # =====================================================================
    # Inter-individual variability (Soto 2014 Table 2, combined final-model
    # column). The paper reports IIV in CV%; log-normal variance is
    # omega^2 = log(1 + CV^2).
    #
    # Cross-drug correlation between eta_CL (ampicillin) and eta_CL
    # (sulbactam) is rho = 0.858 (Table 2; Results, paragraph 5). This
    # appears as the cov off-diagonal of the 2x2 block on lcl + lcl_sbt:
    #   cov(eta_CL_amp, eta_CL_sbt) = rho * sqrt(var_amp * var_sbt).
    # The correlation between eta_V2 of ampicillin and sulbactam was NOT
    # statistically significant (delta OFV = -0.673; Results, paragraph 4)
    # so the two V2 etas are declared independently below.
    # =====================================================================
    etalcl + etalcl_sbt ~ c(
      log(1 + 0.148^2),
      0.858 * sqrt(log(1 + 0.148^2) * log(1 + 0.152^2)),
      log(1 + 0.152^2)
    )
    # Table 2 ampicillin: CV[eta_CL,i(1)] = 14.8% (RSE 15.5%) -> omega^2 = log(1 + 0.148^2);
    # Table 2 sulbactam: CV[eta_CL,i(2)] = 15.2% (RSE 15.2%) -> omega^2 = log(1 + 0.152^2);
    # Table 2: rho[eta_CL,i(1), eta_CL,i(2)] = 0.858 (RSE 34.8%)

    etalvp     ~ log(1 + 0.152^2)
    # Table 2 ampicillin: CV[eta_V2,i(1)] = 15.2% (RSE 36.2%) -> omega^2 = log(1 + 0.152^2)
    etalvp_sbt ~ log(1 + 0.148^2)
    # Table 2 sulbactam: CV[eta_V2,i(2)] = 14.8% (RSE 28.3%) -> omega^2 = log(1 + 0.148^2)

    # =====================================================================
    # Residual error (Soto 2014 Table 2, combined final-model column).
    # The source paper used a log-transform-both-sides additive error
    # model, which is equivalent to a proportional error model in linear
    # space; CV[eps] is the proportional SD (no separate additive term).
    # The cross-drug residual correlation rho[eps_ij(1), eps_ij(2)] = 0.946
    # (Table 2) is documented in the vignette's Assumptions and deviations
    # section -- rxode2 / nlmixr2 typical forward-simulation does not
    # encode cross-output residual correlation for proportional errors,
    # so the residuals here are declared independently per drug.
    # =====================================================================
    propSd     <- 0.242; label("Ampicillin proportional residual SD (fraction)") # Table 2 ampicillin: CV[eps_ij(1)] = 24.2% (RSE 13.5%) under LTBS additive == linear-space proportional
    propSd_sbt <- 0.233; label("Sulbactam proportional residual SD (fraction)")  # Table 2 sulbactam: CV[eps_ij(2)] = 23.3% (RSE 14.4%) under LTBS additive == linear-space proportional
  })

  model({
    # ------------------------------------------------------------
    # Ampicillin individual PK parameters.
    # CL_i = theta_CL * (CRCL / 71)^theta_CLcr * exp(eta_CL)
    # V2_i = theta_V2 * (WT   / 51.2)^theta_BWT * exp(eta_V2)
    # V1 and Q have no covariates.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 71.0)^e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 51.2)^e_wt_vp

    # ------------------------------------------------------------
    # Sulbactam individual PK parameters.
    # ------------------------------------------------------------
    cl_sbt <- exp(lcl_sbt + etalcl_sbt) * (CRCL / 71.0)^e_crcl_cl
    vc_sbt <- exp(lvc_sbt)
    q_sbt  <- exp(lq_sbt)
    vp_sbt <- exp(lvp_sbt + etalvp_sbt) * (WT / 51.2)^e_wt_vp

    # ------------------------------------------------------------
    # Micro-constants.
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_sbt <- cl_sbt / vc_sbt
    k12_sbt <- q_sbt  / vc_sbt
    k21_sbt <- q_sbt  / vp_sbt

    # ------------------------------------------------------------
    # Ampicillin two-compartment IV disposition.
    # Dose goes into central (no absorption phase).
    # ------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central          - k21 * peripheral1

    # ------------------------------------------------------------
    # Sulbactam two-compartment IV disposition (independent
    # compartments; the two drugs do not interconvert).
    # ------------------------------------------------------------
    d/dt(central_sbt)     <- -(kel_sbt + k12_sbt) * central_sbt + k21_sbt * peripheral1_sbt
    d/dt(peripheral1_sbt) <-  k12_sbt * central_sbt              - k21_sbt * peripheral1_sbt

    # ------------------------------------------------------------
    # Observations. dose in mg, vc in L -> Cc in mg/L = ug/mL.
    # ------------------------------------------------------------
    Cc     <- central     / vc
    Cc_sbt <- central_sbt / vc_sbt

    Cc     ~ prop(propSd)
    Cc_sbt ~ prop(propSd_sbt)
  })
}
