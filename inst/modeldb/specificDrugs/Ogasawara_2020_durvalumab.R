Ogasawara_2020_durvalumab <- function() {
  description <- "Two compartment PK model of durvalumab (anti-PD-L1) in patients with hematologic malignancies (Ogasawara 2020)"
  reference <- "Ogasawara K, Newhall K, Maxwell SE, et al. Population Pharmacokinetics of an Anti-PD-L1 Antibody, Durvalumab in Patients with Hematologic Malignancies. Clin Pharmacokinet. 2020;59(2):217-227. doi:10.1007/s40262-019-00804-x"
  vignette <- "Ogasawara_2020_durvalumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power allometric effects on CL (exponent 0.581) and Vc (exponent 0.451), normalized to study-population median 74.7 kg.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on CL (exponent -1.58) and Vc (exponent -0.566), normalized to study-population median 40 g/L. Note: g/L not g/dL; 40 g/L = 4.0 g/dL.",
      source_name        = "ALB"
    ),
    IGG = list(
      description        = "Baseline serum immunoglobulin G",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.258), normalized to study-population median 7.6 g/L. Time-varying in studies MEDI4736-MM-002 and -MM-005; baseline only in NHL-001 and MDS-001.",
      source_name        = "IGG"
    ),
    SPDL1 = list(
      description        = "Baseline soluble PD-L1",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.0617), normalized to study-population median 173.8 pg/mL. Time-varying; values below LOD (67.1 pg/mL) imputed as LOD/2 = 33.55 pg/mL.",
      source_name        = "SPDL1"
    ),
    LDH = list(
      description        = "Baseline lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.115), normalized to study-population median 216 U/L. Time-varying.",
      source_name        = "LDH"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects on CL (factor 0.791) and Vc (factor 0.790) for females. Source column 'SEX' with 0 = male, 1 = female encoding matches the canonical SEXF orientation directly.",
      source_name        = "SEX"
    ),
    MDSAML = list(
      description        = "MDS/AML disease-type indicator (1 = myelodysplastic syndrome or acute myeloid leukemia)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-MDS/AML; includes multiple myeloma and non-Hodgkin lymphoma)",
      notes              = "Multiplicative effect on CL (factor 1.26). Captures the combined MDS and AML subpopulation from Ogasawara 2020 Table 1 (studies MEDI4736-MM-002 and -MDS-001).",
      source_name        = "MDSAML"
    ),
    MM = list(
      description        = "Multiple myeloma disease indicator (1 = active multiple myeloma)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-MM; includes MDS, AML, and non-Hodgkin lymphoma)",
      notes              = "Multiplicative effect on Vc (factor 0.820). Active (non-smoldering) multiple myeloma from studies MEDI4736-MM-002 and -MM-005 in Ogasawara 2020.",
      source_name        = "MM"
    )
  )
  # Note: ADA (anti-drug antibodies) were NOT examined as a covariate in this analysis
  dosing <- "central"
  ini({
    # Structural parameters (time in hours)
    # From Table 3 of the paper
    lcl <- log(0.0107)   ; label("Clearance (CL, L/h)")
    lvc <- log(4.63)     ; label("Central volume of distribution (Vc, L)")
    lq  <- log(0.0376)   ; label("Intercompartmental clearance (Q, L/h)")
    lvp <- log(2.68)     ; label("Peripheral volume of distribution (Vp, L)")

    # Covariate effects on CL (continuous: power model; categorical: multiplicative)
    # From Table 3 footnote b
    e_alb_cl   <- -1.58  ; label("Power exponent of albumin on CL (unitless)")
    e_igg_cl   <- 0.258  ; label("Power exponent of IgG on CL (unitless)")
    e_spdl1_cl <- 0.0617 ; label("Power exponent of soluble PD-L1 on CL (unitless)")
    e_ldh_cl   <- 0.115  ; label("Power exponent of LDH on CL (unitless)")
    e_wt_cl    <- 0.581  ; label("Allometric exponent of body weight on CL (unitless)")
    e_sex_cl   <- 0.791  ; label("Multiplicative factor on CL for females (unitless)")
    e_mdsaml_cl <- 1.26  ; label("Multiplicative factor on CL for MDS/AML (unitless)")

    # Covariate effects on Vc
    # From Table 3 footnote c
    e_alb_vc   <- -0.566 ; label("Power exponent of albumin on Vc (unitless)")
    e_wt_vc    <- 0.451  ; label("Allometric exponent of body weight on Vc (unitless)")
    e_sex_vc   <- 0.790  ; label("Multiplicative factor on Vc for females (unitless)")
    e_mm_vc    <- 0.820  ; label("Multiplicative factor on Vc for multiple myeloma (unitless)")

    # Inter-individual variability (diagonal).  nlmixr2 `etaX ~ value` stores the
    # variance (omega^2).  For log-normal parameters the paper's %CV back-converts
    # via omega^2 = log(1 + CV^2); sqrt(exp(omega^2) - 1) is the %CV actually
    # implied by the stored variance.  sqrt(exp(0.0654) - 1) = 25.6%, not 25.8%.
    etalcl ~ 0.0654   # implies 25.6% CV; paper Table 3 reports 25.8% CV
    etalvc ~ 0.0599   # implies 24.5% CV; paper Table 3 reports 24.7% CV

    # Residual error: log-additive (concentrations were log-transformed;
    # additive error on log scale = proportional-like on original scale)
    addSd <- 0.198     ; label("Log-additive residual error (SD on log scale)")
  })
  model({
    # Covariate-adjusted PK parameters
    # Reference values are study population medians from Table 2:
    #   ALB = 40 g/L, IgG = 7.6 g/L, sPD-L1 = 173.8 pg/mL,
    #   LDH = 216 U/L, WT = 74.7 kg
    # Equations from Table 3 footnotes b and c

    cl <- exp(lcl + etalcl) *
      (ALB / 40)^e_alb_cl *
      (IGG / 7.6)^e_igg_cl *
      (SPDL1 / 173.8)^e_spdl1_cl *
      (LDH / 216)^e_ldh_cl *
      (WT / 74.7)^e_wt_cl *
      e_sex_cl^SEXF *
      e_mdsaml_cl^MDSAML

    vc <- exp(lvc + etalvc) *
      (ALB / 40)^e_alb_vc *
      (WT / 74.7)^e_wt_vc *
      e_sex_vc^SEXF *
      e_mm_vc^MM

    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment model (IV administration)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc   # mg/L

    Cc ~ lnorm(addSd)
  })
}
