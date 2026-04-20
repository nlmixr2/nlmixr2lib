Chua_2025_mirikizumab <- function() {
  description <- "Two-compartment population PK model for mirikizumab (anti-IL-23p19 IgG4 mAb) in patients with moderately-to-severely active Crohn's disease (Chua 2025 VIVID-1 phase 3)"
  reference <- "Chua L, Otani Y, Lin Z, Friedrich S, Durand F, Zhang XC. Mirikizumab Pharmacokinetics and Exposure-Response in Patients With Moderately-To-Severely Active Crohn's Disease: Results From Two Randomized Studies. Clin Transl Sci. 2025;18(8):e70320. doi:10.1111/cts.70320"
  vignette <- "Chua_2025_mirikizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Chua 2025 VIVID-1; power scaling on CL, Q, Vc, Vp with reference weight 65 kg per Table 2 footnote c/d.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Time-varying serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effect on CL per Table 2 footnote c: multiplier = 1 + (ALB - 44.57) * (-0.02). Reference 44.57 g/L (VIVID-1 geometric mean).",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL with reference 7.41 mg/L per Table 2 footnote c. Standard CRP assay (not high-sensitivity); see covariate-columns.md for CRP vs hsCRP distinction.",
      source_name        = "CRP"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on logit of bioavailability: logit(F) = logit(F_pop) + (-0.0354) * (BMI - 24.75). The paper reports the slope (-0.0354, Table 2) and the population F estimate (38.8%) but does not explicitly state the BMI centering value for the VIVID-1 logit parameterization; the 24.75 kg/m^2 centering is adopted from the SERENITY analysis in the same paper (Table 2 footnote e) as a reasonable default.",
      source_name        = "BMI"
    )
  )

  population <- list(
    n_subjects     = 711L,
    n_studies      = 1L,
    age_range      = "18-74 years (phase 3 inclusion: 18-80 years)",
    age_mean       = "36.3 years",
    weight_range   = "21.1-154.0 kg",
    weight_mean    = "68.1 kg",
    sex_female_pct = 43.5,
    race_ethnicity = c(White = 70.7, Asian = 25.0, Other = 4.2),
    disease_state  = "Moderately-to-severely active Crohn's disease (SES-CD >= 7, or >= 4 with isolated ileal disease; inadequate response, loss of response, or intolerance to conventional or biologic therapies).",
    dose_range     = "Induction: 900 mg IV Q4W x 3 doses. Maintenance: 300 mg SC Q4W through Week 52. Also placebo and ustekinumab active-control arms (not used for PK modeling).",
    regions        = "Global, multi-regional",
    notes          = "VIVID-1 (NCT03926130) phase 3 treat-through study, July 2019 - October 2023. PK dataset: 5318 observations from 711 patients (628 mirikizumab-treated during induction + 83 placebo non-responders who crossed over after Week 12). Baseline mean SES-CD 12.9 (3.3-35.5); baseline mean CDAI 319 (92.4-726.3). 62.4% had prior biologic therapy, 50.5% had failed prior biologic therapy (sic from Table 1; see source). Demographics from Chua 2025 Table 1."
  )

  ini({
    # Structural parameters from Chua 2025 Table 2 (VIVID-1 column, final
    # popPK). Typical values are for the reference subject: 65 kg body
    # weight, 44.57 g/L serum albumin, 7.41 mg/L CRP, 24.75 kg/m^2 BMI.
    # Paper reports rates in h^-1; converted to day^-1 for nlmixr2lib
    # convention (x 24).
    lka     <- log(0.00693 * 24); label("First-order SC absorption rate (ka, 1/day)")     # Chua 2025 Table 2 VIVID-1 (0.00693 /h)
    lcl     <- log(0.0197  * 24); label("Clearance at reference covariates (CL, L/day)")  # Chua 2025 Table 2 VIVID-1 (0.0197 L/h)
    lvc     <- log(2.78);         label("Central volume of distribution (V2 / Vc, L)")    # Chua 2025 Table 2 VIVID-1
    lvp     <- log(1.61);         label("Peripheral volume of distribution (V3 / Vp, L)") # Chua 2025 Table 2 VIVID-1
    lq      <- log(0.00921 * 24); label("Intercompartmental clearance (Q, L/day)")        # Chua 2025 Table 2 VIVID-1 (0.00921 L/h)
    logitfdepot <- log(0.388 / (1 - 0.388)); label("Logit of SC bioavailability at reference BMI (unitless)")  # Chua 2025 Table 2 VIVID-1 (F_pop = 38.8%)

    # Allometric-style body-weight exponents (Chua 2025 VIVID-1 Table 2
    # footnote c/d). Shared across CL and Q for clearance and across Vc
    # and Vp for volume.
    e_wt_cl <- 0.268; label("Power exponent of body weight on CL and Q (unitless)")  # Chua 2025 Table 2 VIVID-1
    e_wt_v  <- 0.445; label("Power exponent of body weight on Vc and Vp (unitless)") # Chua 2025 Table 2 VIVID-1

    # Covariate effects on CL (Table 2 footnote c, VIVID-1):
    #   CL_i = CL * (WT/65)^0.268 * (1 + (ALB-44.57) * -0.02) * (CRP/7.41)^0.0853
    e_alb_cl <- -0.0200;  label("Linear-deviation coefficient of albumin on CL (per g/L, relative)") # Chua 2025 Table 2 VIVID-1
    e_crp_cl <-  0.0853;  label("Power exponent of baseline CRP on CL (unitless)")                  # Chua 2025 Table 2 VIVID-1

    # Covariate effect on logit(F) (Table 2 footnote e, VIVID-1):
    #   logit(F_i) = logit(F_pop) + -0.0354 * (BMI - 24.75)
    e_bmi_fdepot <- -0.0354; label("Linear coefficient of BMI on logit of SC bioavailability (per kg/m^2)") # Chua 2025 Table 2 VIVID-1

    # IIV. Table 2 reports IIV as %CV on the linear scale for log-normal
    # parameters: omega^2 = log(CV^2 + 1).
    # CL: 32.3% -> omega^2 = log(0.323^2 + 1) = 0.09927
    # Vc: 17.1% -> omega^2 = log(0.171^2 + 1) = 0.02882
    # Vp: 23.6% -> omega^2 = log(0.236^2 + 1) = 0.05421
    # F:  additive in logit domain; per Table 2 footnote b,
    #   %IIV_F = 100 * sqrt(OMEGA_F), so OMEGA_F = (0.517)^2 = 0.26729.
    # Ka:  no IIV reported in VIVID-1 (dash in Table 2).
    etalcl         ~ 0.09927   # Chua 2025 Table 2 VIVID-1 (32.3% CV)
    etalvc         ~ 0.02882   # Chua 2025 Table 2 VIVID-1 (17.1% CV)
    etalvp         ~ 0.05421   # Chua 2025 Table 2 VIVID-1 (23.6% CV)
    etalogitfdepot ~ 0.26729   # Chua 2025 Table 2 VIVID-1 (OMEGA_F on logit scale; 51.7% IIV)

    # Residual error. Chua 2025 does not tabulate residual-error estimates
    # (they are shown only in the supplementary model code, Figure S13,
    # which is not available as structured data). A 20% proportional
    # residual error with a small additive floor is adopted as a
    # reasonable default for an IgG mAb popPK model; see the vignette
    # Assumptions section for detail.
    propSd <- 0.20;  label("Proportional residual error (SD, fraction) -- not published; assumed for simulation")
    addSd  <- 0.10;  label("Additive residual error (SD, ug/mL) -- not published; assumed for simulation")
  })
  model({
    # Individual PK parameters. Reference covariates per Table 2
    # footnotes c (CL, Q), d (Vc, Vp), and e (F): WT_ref = 65 kg,
    # ALB_ref = 44.57 g/L, CRP_ref = 7.41 mg/L, BMI_ref = 24.75 kg/m^2.

    # Covariate multipliers
    wt_cl <- (WT / 65)^e_wt_cl
    wt_v  <- (WT / 65)^e_wt_v
    alb_cl <- 1 + (ALB - 44.57) * e_alb_cl
    crp_cl <- (CRP / 7.41)^e_crp_cl

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt_cl * alb_cl * crp_cl
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp + etalvp) * wt_v
    q  <- exp(lq)           * wt_cl

    # Bioavailability on logit scale with linear BMI effect.
    logit_f <- logitfdepot + etalogitfdepot + e_bmi_fdepot * (BMI - 24.75)
    fdepot  <- 1 / (1 + exp(-logit_f))

    # Two-compartment model with first-order SC absorption.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    # SC bioavailability applied to the depot. IV doses bypass depot and
    # go directly to central via cmt=central in the event table.
    f(depot) <- fdepot

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
