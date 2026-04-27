Wang_2017_benralizumab <- function() {
  description <- "Two compartment PK model of benralizumab (anti-IL-5Ralpha) in healthy volunteers and patients with asthma (Wang 2017)"
  reference <- "Wang B, Lau YY, Liang M, et al. Population Pharmacokinetics and Pharmacodynamics of Benralizumab in Healthy Volunteers and Patients With Asthma. CPT Pharmacometrics Syst Pharmacol. 2017;6(4):249-257. doi:10.1002/psp4.12160"
  vignette <- "Wang_2017_benralizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL (exponent fixed at 0.75), Vc (estimated exponent allovc), and Vp (estimated exponent allovp), normalized to study-population mean 77 kg (Table 2; median not published).",
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody status, high-titer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (negative or low-titer, i.e., titer < 400)",
      notes              = "1 = high-titer ADA with titer >= 400; time-varying (assessed at each predesignated sampling visit). Called 'ADAs' in the original publication. exp(1.52) = ~4.6-fold increase in CL. Source column name 'ADA'.",
      source_name        = "ADA"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese heritage indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Multiplicative factor on Vc (e_jpn_vc = 1.34). Wang 2017 enrolled a dedicated Japanese healthy-volunteer cohort alongside non-Japanese subjects; the indicator captures both the racial and cohort effect. Source column name 'JAPANESE_HV'.",
      source_name        = "JAPANESE_HV"
    )
  )
  dosing <- c("depot", "central")
  ini({
    # Structural parameters
    lcl <- log(0.323)    ; label("Clearance (CL, L/day)")
    lvc <- log(3.16)     ; label("Central volume of distribution (Vc, L)")
    lq  <- log(0.939)    ; label("Intercompartmental clearance (Q, L/day)")
    lvp <- log(2.83)     ; label("Peripheral volume of distribution (Vp, L)")
    lka <- log(0.252)    ; label("Absorption rate constant (ka, 1/day)")
    lfdepot <- log(0.526) ; label("Subcutaneous bioavailability (F, fraction)")

    # Allometric scaling: CL fixed at 0.75, volumes estimated
    # CL exponent = 0.75 (fixed, standard allometric)
    allovc  <- 0.651     ; label("Allometric exponent of body weight on Vc (unitless)")
    allovp  <- 0.576     ; label("Allometric exponent of body weight on Vp (unitless)")

    # Covariate effects
    # ADA effect uses exp() form per Eq. 6 in paper: CL = theta1 * exp(I_ADA * theta2)
    # Table 3 footnote c: "Natural exponent" — so 1.52 is the exponent, not a multiplier
    e_ada_cl   <- 1.52   ; label("Natural exponent of ADA on CL: exp(1.52) = 4.57-fold increase (unitless)")
    # Race effect is a multiplicative fraction per Table 3 (no footnote c)
    e_jpn_vc   <- 1.34   ; label("Multiplicative factor on Vc for Japanese healthy volunteers (unitless)")

    # Inter-individual variability (diagonal)
    etalcl ~ 0.0608   # 24.9% CV
    etalvc ~ 0.0524   # 23.1% CV
    etalq  ~ 0.0791   # 28.5% CV
    etalvp ~ 0.1029   # 32.7% CV
    etalka ~ 0.2444   # 52.7% CV (omega^2 = log(1 + 0.527^2))
    etalfdepot ~ 0.1013   # 32.4% CV

    # Residual error (combined, separate for SC and IV routes)
    # Using SC route residual error as the default
    propSd <- 0.142    ; label("Proportional residual error, SC (fraction)")
    addSd  <- 0.0356   ; label("Additive residual error, SC (mg/L)")
  })
  model({
    # Covariate-adjusted PK parameters
    # Allometric scaling: paper states "normalized to the population median"
    # (Methods, PK model description). The median body weight is not published;
    # Table 2 reports only mean 77.0 ± 19.0 kg. Using 77 kg as a best-available
    # proxy for the unpublished median. If the median is later recovered from
    # a supplement or regulatory review, update this reference value.
    cl <- exp(lcl + etalcl) * (WT / 77)^0.75 * exp(e_ada_cl * ADA_POS)
    vc <- exp(lvc + etalvc) * (WT / 77)^allovc * e_jpn_vc^RACE_JAPANESE
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 77)^allovp
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    # Two-compartment model with first-order absorption (SC) or IV
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc    # mg/L (= ug/mL)

    Cc ~ add(addSd) + prop(propSd)
  })
}
