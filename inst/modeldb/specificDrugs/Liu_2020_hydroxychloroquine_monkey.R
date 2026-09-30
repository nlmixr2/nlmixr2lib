Liu_2020_hydroxychloroquine_monkey <- function() {
  description <- "Preclinical (cynomolgus macaque). Five-compartment population PK model of oral (intragastric) hydroxychloroquine fitted jointly to plasma, whole-blood and lung-tissue concentrations in 17 male cynomolgus macaques. First-order absorption from a depot into a central (plasma) compartment with first-order elimination (Ke) and a linear exchange with an unobserved peripheral compartment (Kcp / Kpc). Transfer from the central compartment into a red-blood-cell compartment and into a lung-tissue compartment is saturable in the central AMOUNT (Kmax * A1 / (A50 + A1), Acb50 = 7.05 mg and Acl50 = 0.498 mg), with first-order return to the central compartment (Kbc, Klc). The whole-blood concentration is the sum of the red-cell and plasma concentrations (A3 / Vb + A1 / Vc). Fitted by the MLEM algorithm in ADAPT 5 with log-normal between-animal variability on all 13 structural parameters and a proportional-plus-additive residual error per matrix. The printed parameters reproduce the observed single-dose (3 mg/kg) whole-blood and plasma NCA and the paper's own population predictions for lung (Figure 5B), but not its simulated Figure 9 / Table 5 profiles, whose regimen is unstated (see the vignette)."
  reference <- "Liu Q, Bi G, Chen G, Guo X, Tu S, Tong X, Xu M, Liu M, Wang B, Jiang H, Wang J, Li H, Wang K, Liu D, Song C. Time-Dependent Distribution of Hydroxychloroquine in Cynomolgus Macaques Using Population Pharmacokinetic Modeling Method. Front Pharmacol. 2021;11:602880. doi:10.3389/fphar.2020.602880. PMCID: PMC7841297. Structural model: Equations 1-8 and Figure 2 (Methods, 'Structural PK Model'). Parameter estimates and between-animal variability: Table 4. Error model and estimation method: Methods, 'Structural PK Model' and 'Software and Platform Used'."
  vignette <- "Liu_2020_hydroxychloroquine_monkey"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL (plasma and whole blood); ng/g (lung tissue)")

  compartmentData <- list(
    depot = list(analyte = "hydroxychloroquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxychloroquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hydroxychloroquine", units = "mg", specimen = "not applicable", verified = TRUE),
    # Carried as an AMOUNT (mg) with its own apparent volume Vb, exactly as the
    # paper's A3 -- not in the concentration units the rbc_<analyte> family
    # usually carries. Recorded here so the unit deviation is machine-readable.
    rbc_hcq = list(analyte = "hydroxychloroquine", units = "mg", specimen = "blood cell", verified = TRUE),
    lung = list(analyte = "hydroxychloroquine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus macaque (Macaca fascicularis)",
    n_subjects = 17L,
    n_studies = 1L,
    sex_female_pct = 0,
    weight_range = "4.13 +/- 0.43 kg (mean +/- SD; Methods, 'Experimental Animals')",
    disease_state = "Healthy animals (biodistribution study motivated by COVID-19 repurposing of hydroxychloroquine)",
    dose_range = "Hydroxychloroquine 1-21 mg/kg intragastrically, in six regimens: single 3 mg/kg; 3 mg/kg twice 4 h apart on days 1-2; 2 mg/kg twice on day 1 then 1 mg/kg twice daily on days 2-3; 6 mg/kg twice on day 1 then 2 mg/kg twice daily on days 2-3; 21 mg/kg twice on day 1; 21 mg/kg twice on day 1 then 7 mg/kg twice daily on days 2-5. Doses within a day were given 4 h apart.",
    regions = "China (Pharmaron Beijing Co., Ltd.; Peking University Third Hospital)",
    notes = paste(
      "Six groups (n = 3 each, except group F n = 2), all male. 141 plasma and 149 whole-blood samples (intensive sampling to 72 h in the low-dose groups, to 264 h in group E and to 120 h in group F), and one terminal lung sample in 14 of the 17 animals at 120-504 h (Methods, 'Subjects and Study Design').",
      "The group letters are inconsistent between the Figure 1 caption and the Methods text; the dose levels and the six regimen shapes are the same in both.",
      "LLOQ 2.00 ng/mL in plasma, 5.00 ng/mL in whole blood and 0.40 ng/mL in tissue homogenate (Results, 'Validation of HPLC-MS/MS').",
      "The paper does not state whether the mg doses are expressed as hydroxychloroquine sulfate (the administered material) or as base."
    )
  )

  ini({
    # All values: Table 4 ('Mean'), the typical values of the final model
    # (the abstract calls Ke = 0.236 1/h the value 'for a typical monkey').
    # Table 4 footnote: 'Due to limited sample size, the RSD% cannot be
    # estimated', so no standard errors are available.
    lka <- log(0.592); label("Absorption rate constant Ka (1/h)")                                         # Table 4: Ka = 0.592 1/h
    lkel <- log(0.236); label("Elimination rate constant from the central (plasma) compartment Ke (1/h)") # Table 4: Ke = 0.236 1/h
    lvc <- log(114); label("Apparent central (plasma) volume Vc/F (L)")                                   # Table 4: Vc/F = 114 L
    lk12 <- log(0.600); label("Central-to-peripheral rate constant Kcp (1/h)")                            # Table 4: Kcp = 0.600 1/h
    lk21 <- log(0.514); label("Peripheral-to-central rate constant Kpc (1/h)")                            # Table 4: Kpc = 0.514 1/h

    # Red-blood-cell compartment (paper A3, volume Vb). Kcbmax is printed in
    # Table 4 with no unit and as 'h-1' in the Results text; Equation 1 and the
    # Figure 2 arrow 'Kcbmax/(Acb50 + A1)' make the flux Kcbmax * A1 /
    # (Acb50 + A1) with A1 in mg, so Kcbmax is carried as an amount rate (mg/h).
    lvmax_rbc <- log(2.48); label("Maximum central-to-red-cell transfer rate Kcbmax (mg/h)")                    # Table 4: Kcbmax = 2.48
    lkm_rbc <- log(7.05); label("Central amount giving half-maximal central-to-red-cell transfer Acb50 (mg)")   # Table 4: Acb50 = 7.05 mg
    lkeff_rbc <- log(0.718); label("Red-cell-to-central rate constant Kbc (1/h)")                               # Table 4: Kbc = 0.718 1/h
    lv_rbc <- log(2.68); label("Apparent red-blood-cell compartment volume Vb/F (L)")                           # Table 4: Vb/F = 2.68 L

    # Lung-tissue compartment (paper A4, volume VL). Same unit reasoning for
    # Kclmax as for Kcbmax above.
    lvmax_lung <- log(1.92); label("Maximum central-to-lung transfer rate Kclmax (mg/h)")                       # Table 4: Kclmax = 1.92
    lkm_lung <- log(0.498); label("Central amount giving half-maximal central-to-lung transfer Acl50 (mg)")     # Table 4: Acl50 = 0.498 mg
    lk_lung_central <- log(0.159); label("Lung-to-central rate constant Klc (1/h)")                             # Table 4: Klc = 0.159 1/h
    lv_lung <- log(5.55); label("Apparent lung compartment volume VL/F (L)")                                    # Table 4: VL/F = 5.55 L

    # Between-animal variability: Table 4 'IIV CV%' of a log-normal parameter,
    # converted as omega^2 = log(1 + CV^2). ADAPT's MLEM fits a multivariate
    # log-normal distribution, but only the CVs are printed, so the matrix is
    # encoded diagonal.
    etalka ~ 0.08075             # Table 4: Ka IIV CV 29%
    etalkel ~ 0.25839            # Table 4: Ke IIV CV 54.3%
    etalvc ~ 0.27533             # Table 4: Vc/F IIV CV 56.3%
    etalk12 ~ 0.91164            # Table 4: Kcp IIV CV 122%
    etalk21 ~ 0.58041            # Table 4: Kpc IIV CV 88.7%
    etalvmax_rbc ~ 0.99920       # Table 4: Kcbmax IIV CV 131%
    etalkm_rbc ~ 0.89200         # Table 4: Acb50 IIV CV 120%
    etalkeff_rbc ~ 0.24426       # Table 4: Kbc IIV CV 52.6%
    etalv_rbc ~ 0.94098          # Table 4: Vb/F IIV CV 125%
    etalvmax_lung ~ 0.48398      # Table 4: Kclmax IIV CV 78.9%
    etalkm_lung ~ 1.77975        # Table 4: Acl50 IIV CV 222%
    etalk_lung_central ~ 0.63818 # Table 4: Klc IIV CV 94.5%
    etalv_lung ~ 0.55366         # Table 4: VL/F IIV CV 86%

    # Residual error: ADAPT proportional-plus-additive model, SD = SDinter +
    # SDslope * Y (combined1). The additive SDs were fixed at 1e-5.
    propSd <- 0.286; label("Proportional residual error, plasma (fraction)")                    # Table 4: PD pl = 0.286
    addSd <- fixed(0.00001); label("Additive residual error SD, plasma (ng/mL)")                # Table 4: SD pl = 0.100E-04 ng/mL, Fixed
    propSd_Cblood <- 0.118; label("Proportional residual error, whole blood (fraction)")        # Table 4: PD bl = 0.118
    addSd_Cblood <- fixed(0.00001); label("Additive residual error SD, whole blood (ng/mL)")    # Table 4: SD bl = 0.100E-04 ng/mL, Fixed
    propSd_Clung <- 0.0514; label("Proportional residual error, lung tissue (fraction)")        # Table 4: PD lu = 0.514E-01
    addSd_Clung <- fixed(0.00001); label("Additive residual error SD, lung tissue (ng/g)")      # Table 4: SD lu = 0.100E-04 ng/g, Fixed
  })

  model({
    ka <- exp(lka + etalka)
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    vmax_rbc <- exp(lvmax_rbc + etalvmax_rbc)
    km_rbc <- exp(lkm_rbc + etalkm_rbc)
    keff_rbc <- exp(lkeff_rbc + etalkeff_rbc)
    v_rbc <- exp(lv_rbc + etalv_rbc)
    vmax_lung <- exp(lvmax_lung + etalvmax_lung)
    km_lung <- exp(lkm_lung + etalkm_lung)
    k_lung_central <- exp(lk_lung_central + etalk_lung_central)
    v_lung <- exp(lv_lung + etalv_lung)

    # Saturable transfer out of the central compartment, driven by the
    # central AMOUNT (Equations 1, 3 and 4; Figure 2).
    flux_rbc <- vmax_rbc * central / (km_rbc + central)
    flux_lung <- vmax_lung * central / (km_lung + central)

    d/dt(depot) <- -ka * depot                                                              # Equation 2
    d/dt(central) <- ka * depot + keff_rbc * rbc_hcq - flux_rbc - kel * central -
      flux_lung + k_lung_central * lung - k12 * central + k21 * peripheral1                 # Equation 1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1                                  # Equation 5
    d/dt(rbc_hcq) <- flux_rbc - keff_rbc * rbc_hcq                                          # Equation 3
    d/dt(lung) <- flux_lung - k_lung_central * lung                                         # Equation 4

    # Amounts in mg and volumes in L give mg/L; x 1000 gives ng/mL (ng/g for
    # lung, taking tissue density as 1 g/mL as the paper's units imply).
    Cc <- 1000 * central / vc                          # Equation 6
    Cblood <- 1000 * (rbc_hcq / v_rbc + central / vc)  # Equation 7
    Clung <- 1000 * lung / v_lung                      # Equation 8

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Cblood ~ add(addSd_Cblood) + prop(propSd_Cblood) + combined1()
    Clung ~ add(addSd_Clung) + prop(propSd_Clung) + combined1()
  })
}
