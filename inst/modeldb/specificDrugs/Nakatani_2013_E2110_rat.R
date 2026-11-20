Nakatani_2013_E2110_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment oral PK with first-order absorption",
    "linked through an effect compartment to central serotonin 1A (5-HT1A)",
    "receptor occupancy (RO) measured by [11C]WAY-100635 PET in the medial",
    "prefrontal cortex (MPFC) and the dorsal raphe nucleus (DRN) for E2110, an",
    "investigational 5-HT1A receptor antagonist for overactive bladder, in",
    "male and female Sprague-Dawley rats (Nakatani 2013, time-course analysis).",
    "The PK model was fitted separately to the mean plasma profiles of each",
    "sex after 1 mg/kg orally; the PD model (hyperbolic Emax with Emax = 100",
    "percent, region-specific EC50 on the plasma-equivalent effect-site",
    "concentration, one equilibration rate constant ke0) was fitted in male",
    "rats and applied to female PK to simulate female DRN occupancy. The paper",
    "prints the PD parameters (Table 1 EC50, ke0 as an equilibration half-life)",
    "but not the PK parameters; those were read back from the authors' drawn",
    "fitted curves in Figure 6 (RMSE 0.0016 log10 units for males, 0.0015 for",
    "females) and are apparent oral values (bioavailability unknown). Doses",
    "are in mg/kg and volumes and clearances per kg body weight.",
    "Typical-value model: mean-data fits, so no interindividual variability",
    "and no residual error were reported (residual terms are set to zero).",
    sep = " "
  )
  reference <- paste(
    "Nakatani Y, Suzuki M, Tokunaga M, Maeda J, Sakai M, Ishihara H,",
    "Yoshinaga T, Takenaka O, Zhang MR, Suhara T, Higuchi M. A small-animal",
    "pharmacokinetic/pharmacodynamic PET study of central serotonin 1A",
    "receptor occupancy by a potential therapeutic agent for overactive",
    "bladder. PLoS One. 2013;8(9):e75040. doi:10.1371/journal.pone.0075040.",
    "PMCID PMC3781034, PMID 24086433.",
    sep = " "
  )
  vignette <- "Nakatani_2013_E2110_rat"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex, female indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "1 = female, 0 = male. Selects the female two-compartment PK fit of",
        "Figure 6 (open circles, dotted line) instead of the male fit (filled",
        "circles, solid line). The authors fitted each sex separately; the",
        "female set is encoded as log-scale offsets from the male set so that",
        "SEXF = 1 reproduces the female fit exactly. The PD parameters were",
        "estimated in males and applied unchanged to females (Figure 8).",
        sep = " "
      ),
      source_name = "male / female SD rats"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "E2110", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "E2110", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "E2110", units = "mg/kg", specimen = "tissue", verified = TRUE),
    # Plasma-equivalent effect-site concentration (Equation 5), not an amount.
    effect = list(analyte = "E2110", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    disease_state = paste(
      "Healthy rats (PK and PET studies). Isoflurane-anaesthetised during the",
      "PET scans.",
      sep = " "
    ),
    dose_range = paste(
      "E2110 orally: 1 mg/kg single dose (PK study, both sexes; PET time-course",
      "study, males); 0.3, 1, 3 and 10 mg/kg (PET dose-occupancy study, males);",
      "0.1 mg/kg (pilot PET study, females).",
      sep = " "
    ),
    regions = "Japan (National Institute of Radiological Sciences and Eisai Co., Ltd.)",
    notes = paste(
      "PK study: 4 male and 4 female rats, 1 mg/kg orally, plasma sampled at",
      "0.25, 0.5, 1, 2, 4, 6 and 8 h; the mean profile of each sex was fitted.",
      "PET time-course study: 4 male rats scanned 2, 4, 6 or 8 h after",
      "1 mg/kg (n = 4 per time point), RO in MPFC and DRN relative to each",
      "rat's own baseline scan. PET dose-occupancy study: 4 male rats per dose",
      "scanned 4 h post-dose (fitted by the separate direct Emax model,",
      "Nakatani_2013_E2110_direct_rat). Female pilot PET: 3 rats, 0.1 mg/kg,",
      "4 h. Plasma unbound fraction 6 percent, concentration-independent.",
      "PK and PD parameters were estimated with SCIENTIST 2.01.",
      sep = " "
    )
  )

  ini({
    # PK: male two-compartment fit (Figure 6 solid line). The paper prints no
    # PK parameter values; these were read back by digitising the drawn line
    # from the publisher's full-resolution figure and least-squares fitting
    # the same model to it (RMSE 0.0016 log10 units over 1437 points). The
    # line joins model predictions at 0.1, 0.25, 0.4, 0.5, 0.75 and 1-10 h.
    # Terminal half-life of the read-back curve 1.51 h (Table 2 male t1/2
    # 1.5 h).
    lka <- log(2.352); label("Absorption rate constant, males (1/h)") # Figure 6 male line, read back; not printed
    lvc <- log(9.628); label("Apparent central volume V1/F per kg, males (L/kg)") # Figure 6 male line, read back; not printed
    lcl <- log(15.16); label("Apparent clearance CL/F per kg, males (L/h/kg)") # Figure 6 male line, read back; not printed
    lq <- log(5.305); label("Apparent intercompartmental clearance Q/F per kg, males (L/h/kg)") # Figure 6 male line, read back; not printed
    lvp <- log(7.742); label("Apparent peripheral volume V2/F per kg, males (L/kg)") # Figure 6 male line, read back; not printed

    # PK: female two-compartment fit (Figure 6 dotted line), read back the
    # same way (RMSE 0.0015 log10 units over 756 points; terminal half-life
    # 2.65 h, Table 2 female t1/2 2.6 h). CL/F is well determined; ka, V1/F,
    # Q/F and V2/F trade off along a ridge (ka 1.2-1.8 1/h fits the line
    # equally well), so the least-squares optimum is used.
    e_sexf_ka <- log(1.344 / 2.352); label("Log-scale effect of female sex on ka (unitless)") # Figure 6 female line, read back: ka 1.344 1/h
    e_sexf_vc <- log(3.339 / 9.628); label("Log-scale effect of female sex on V1/F (unitless)") # Figure 6 female line, read back: V1/F 3.339 L/kg
    e_sexf_cl <- log(1.830 / 15.16); label("Log-scale effect of female sex on CL/F (unitless)") # Figure 6 female line, read back: CL/F 1.830 L/h/kg
    e_sexf_q <- log(2.161 / 5.305); label("Log-scale effect of female sex on Q/F (unitless)") # Figure 6 female line, read back: Q/F 2.161 L/h/kg
    e_sexf_vp <- log(2.537 / 7.742); label("Log-scale effect of female sex on V2/F (unitless)") # Figure 6 female line, read back: V2/F 2.537 L/kg

    # PD: effect-compartment model (Equation 5) with hyperbolic Emax
    # (Equation 4). Discussion: CNS equilibration half-life of 1.4 h.
    lke0 <- log(log(2) / 1.4); label("Effect-compartment equilibration rate constant ke0 (1/h)") # Discussion: equilibration t1/2 1.4 h, ke0 = ln2/1.4
    lemax <- fixed(log(100)); label("Maximum receptor occupancy Emax (percent)") # Results: ROmax held at 100 because 3 and 10 mg/kg gave full occupancy
    lec50_mpfc <- log(2.87); label("Plasma-equivalent concentration for half-maximal RO in MPFC (ng/mL)") # Table 1 time-course study, MPFC: 5.85 nM; Results text 2.87 ng/mL
    lec50_drn <- log(3.6); label("Plasma-equivalent concentration for half-maximal RO in DRN (ng/mL)") # Table 1 time-course study, DRN: 7.38 nM; Results text 3.6 ng/mL

    # Residual error: the paper fitted mean data and reports no error model.
    propSd <- fixed(0); label("Proportional residual SD on plasma concentration (fraction)") # not reported
    addSd_RO_MPFC <- fixed(0); label("Additive residual SD on MPFC receptor occupancy (percent)") # not reported
    addSd_RO_DRN <- fixed(0); label("Additive residual SD on DRN receptor occupancy (percent)") # not reported
  })

  model({
    # Individual parameters (sex-specific PK)
    ka <- exp(lka + e_sexf_ka * SEXF)
    vc <- exp(lvc + e_sexf_vc * SEXF)
    cl <- exp(lcl + e_sexf_cl * SEXF)
    q <- exp(lq + e_sexf_q * SEXF)
    vp <- exp(lvp + e_sexf_vp * SEXF)

    ke0 <- exp(lke0)
    emax <- exp(lemax)
    ec50_mpfc <- exp(lec50_mpfc)
    ec50_drn <- exp(lec50_drn)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg/kg and volume in L/kg give mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc

    # Effect compartment (Equation 5): dCe/dt = ke0 * (Cp - Ce)
    d/dt(effect) <- ke0 * (Cc - effect)

    # Receptor occupancy (percent), Equation 4 on the effect-site concentration
    RO_MPFC <- emax * effect / (ec50_mpfc + effect)
    RO_DRN <- emax * effect / (ec50_drn + effect)

    Cc ~ prop(propSd)
    RO_MPFC ~ add(addSd_RO_MPFC)
    RO_DRN ~ add(addSd_RO_DRN)
  })
}
