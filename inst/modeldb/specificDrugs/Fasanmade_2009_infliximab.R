Fasanmade_2009_infliximab <- function() {
  description <- "Two compartment PK model of infliximab (anti-TNF-alpha) in patients with ulcerative colitis (Fasanmade 2009)"
  reference <- "Fasanmade AA, Adedokun OJ, Blank M, Zhou H, Davis HM. Population pharmacokinetic analysis of infliximab in patients with ulcerative colitis. Eur J Clin Pharmacol. 2009;65(12):1211-1228. doi:10.1007/s00228-009-0718-4"
  units <- list(time = "day", dosing = "mg")
  covariateData <- list(
    WT = "Body weight in kg",
    ALB = "Serum albumin in g/dL",
    ADA = "Anti-drug antibody status (0 = negative, 1 = positive; called ATI [antibodies to infliximab] in the original publication; time-invariant: positive if detected at any visit through week 42/54)",
    SEX = "Sex (0 = male, 1 = female)"
  )
  dosing <- "central"
  ini({
    # Structural parameters (typical values for 77-kg male,
    # albumin 4.1 g/dL, ADA-negative)
    lcl <- log(0.407)    ; label("Clearance (CL, L/day)")
    lv1 <- log(3.29)     ; label("Central volume of distribution (V1, L)")
    lv2 <- log(4.13)     ; label("Peripheral volume of distribution (V2, L)")
    lq  <- log(7.14)     ; label("Intercompartmental clearance (Q, L/day)")

    # Covariate effects
    e_alb_cl  <- -1.54   ; label("Power exponent of albumin on CL (unitless)")
    e_ada_cl  <- 0.471   ; label("Fractional change in CL for ADA-positive (unitless)")
    e_sex_cl  <- -0.236  ; label("Fractional change in CL for females (unitless)")
    e_wt_v1   <- 0.538   ; label("Allometric exponent of body weight on V1 (unitless)")
    e_sex_v1  <- -0.137  ; label("Fractional change in V1 for females (unitless)")

    # Inter-individual variability (diagonal)
    etacl ~ 0.131   # 37.68% CV
    etav1 ~ 0.048   # 22.11% CV

    # Residual error (combined)
    prop.err <- 0.403    ; label("Proportional residual error (fraction)")
    add.err  <- 0.0413   ; label("Additive residual error (ug/mL)")
  })
  model({
    # Covariate-adjusted PK parameters
    # Reference: 77-kg male, albumin 4.1 g/dL, ADA-negative
    cl <- exp(lcl + etacl) *
      (ALB / 4.1)^e_alb_cl *
      (1 + e_ada_cl * ADA) *
      (1 + e_sex_cl * SEX)

    v1 <- exp(lv1 + etav1) *
      (WT / 77)^e_wt_v1 *
      (1 + e_sex_v1 * SEX)

    v2 <- exp(lv2)
    q  <- exp(lq)

    # Two-compartment model (IV administration)
    kel  <- cl / v1
    k12  <- q / v1
    k21  <- q / v2

    d/dt(central)    <- -kel * central - k12 * central + k21 * peripheral
    d/dt(peripheral) <-                  k12 * central - k21 * peripheral

    Cc <- central / v1

    Cc ~ add(add.err) + prop(prop.err)
  })
}
