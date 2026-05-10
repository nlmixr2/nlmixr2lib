Tikiso_2021_abacavir <- function() {
  description <- "Two-compartment population PK model for oral abacavir in HIV-infected African children (Tikiso 2021), with a Savic 2007-style analytical transit-compartment chain feeding a first-order absorption depot, allometric body-weight scaling on disposition (0.75 on CL/Q, 1 on Vc/Vp at 70 kg), sigmoidal Hill-type maturation of CL on postmenstrual age, and multiplicative covariate effects of efavirenz co-medication on CL, rifampicin + super-boosted lopinavir/ritonavir co-medication on F, fixed-dose-combination tablet formulation on MTT, and a time-decaying malnutrition effect on F and CL."
  reference <- "Tikiso T, McIlleron H, Burger D, Gibb D, Rabie H, Lee J, Lallemant M, Cotton MF, Archary M, Hennig S, Denti P. Abacavir pharmacokinetics in African children living with HIV: A pooled analysis describing the effects of age, malnutrition and common concomitant medications. Br J Clin Pharmacol. 2021;1-13. doi:10.1111/bcp.14984"
  vignette <- "Tikiso_2021_abacavir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used for allometric scaling with reference weight 70 kg (exponent 0.75 on CL and Q; exponent 1 on Vc and Vp); see Tikiso 2021 Section 2.4 and Table 3.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (postnatal age + assumed gestational age of 9 months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the sigmoidal Hill maturation function on CL (Tikiso 2021 Eq. 1). Tikiso 2021 lacked subject-level gestational ages and assumed 9 months for every subject; for new simulations, supply PAGE = postnatal_age_months + 9 unless a real gestational age is available.",
      source_name        = "PMAGE"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz indicator (1 = on EFV-based ART, 0 = on standard LPV/r 4:1 reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (standard LPV/r 4:1)",
      notes              = "Multiplicative effect on CL relative to the LPV/r 4:1 reference; +12% in Tikiso 2021 Table 4 (efavirenz is a known UGT inducer).",
      source_name        = "EFV"
    ),
    CONMED_RIF_LPVR4 = list(
      description        = "Concomitant rifampicin-based antitubercular treatment with super-boosted lopinavir/ritonavir 4:4 indicator (1 = RIF + super-boosted LPV/r 4:4, 0 = standard LPV/r 4:1 or EFV reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (standard LPV/r 4:1 or EFV)",
      notes              = "Multiplicative effect on F relative to the LPV/r 4:1 reference; -29.4% in Tikiso 2021 Table 4. The reduction is attributed to PXR-mediated UGT induction by rifampicin (and possibly extra ritonavir).",
      source_name        = "RIF"
    ),
    FORM_TABLET = list(
      description        = "Abacavir tablet formulation indicator (1 = abacavir + lamivudine fixed-dose-combination tablet, 0 = abacavir liquid solution)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (liquid solution)",
      notes              = "Multiplicative effect on the absorption mean transit time MTT relative to the liquid reference; +24.9% (slower absorption) in Tikiso 2021 Table 4. Distinct from the FDC canonical (Wilkins 2008 antitubercular FDC vs SDC); see covariate-columns.md.",
      source_name        = "FORM_TABLET"
    ),
    MAL_NOURISH = list(
      description        = "Malnutrition indicator at study entry (1 = malnourished per WHO definition, height-for-age and weight-for-age Z-scores both < -2.0; 0 = not malnourished)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not malnourished)",
      notes              = "Subject-level time-fixed covariate flagging malnutrition status at the start of nutritional supplementation. Combined with T_NUT_SUPP via mal_decay = MAL_NOURISH * exp(-T_NUT_SUPP * log(2) / e_mal_thalf) to drive the time-decaying malnutrition effect on F (+115% at t=0) and on CL (-64% at t=0); see Tikiso 2021 Eq. 2 and Table 4.",
      source_name        = "MAL"
    ),
    T_NUT_SUPP = list(
      description        = "Time on nutritional supplementation in days (0 at start of supplementation; increases with time)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the exponential decay of the malnutrition effect with half-life 12.2 days (Tikiso 2021 Eq. 2 and Table 4). For non-malnourished subjects (MAL_NOURISH = 0) the value is irrelevant because the malnutrition decay is gated by MAL_NOURISH; supply 0 by default. For fully-recovered malnourished subjects, supply a large value (e.g., 100) so the decay reaches near zero.",
      source_name        = "TNUTRI"
    )
  )

  population <- list(
    n_subjects     = 230L,
    n_studies      = 4L,
    age_range      = "0.1-12.8 years (postnatal)",
    age_median     = "2.1 years",
    weight_range   = "2.5-30.0 kg",
    weight_median  = "9.8 kg",
    sex_female_pct = 52.6,
    disease_state  = "HIV-infected African children on abacavir-containing combination antiretroviral therapy; 115 of 230 (50.0%) malnourished by WHO height-for-age and weight-for-age Z-score criteria, 104 of 230 (45.2%) co-infected with tuberculosis.",
    dose_range     = "Oral abacavir 8 mg/kg twice daily or 16 mg/kg once daily per WHO weight-band guidelines; total daily dose 16 mg/kg.",
    regions        = "South Africa, Uganda, Zambia, Zimbabwe",
    co_medication  = "154 of 230 on lopinavir/ritonavir 4:1 (standard LPV/r); 76 on efavirenz; 104 on rifampicin-based antitubercular treatment of whom 101 received super-boosted LPV/r 4:4 and 3 received EFV.",
    studies        = "Pooled analysis of ARROW (Uganda/Zimbabwe; n=41 children on abacavir + EFV), CHAPAS-3 (Uganda/Zambia; n=27 children on abacavir + EFV), DNDi (South Africa; n=87 children on abacavir + LPV/r 4:1 plus rifampicin-based TB treatment), and MATCH (South Africa; n=75 severely malnourished children on abacavir + LPV/r 4:1).",
    notes          = "Baseline demographics from Tikiso 2021 Tables 1 and 2. n_subjects = 230 children with 2760 plasma concentrations of which 285 (10.3%) were below the lower limit of quantification (mostly pre-dose). Lower limits of quantification were 0.0243 ug/mL (DNDi) and 0.0238 ug/mL (CHAPAS, MATCH); the ARROW LLOQ was reported separately as 0.0243 ug/mL (LC-MS/MS at the University of Cape Town; HPLC at GlaxoSmithKline for ARROW). Sex breakdown: 109 male and 121 female (52.6% female)."
  )

  ini({
    # Structural PK parameters. Reference values quoted at 70 kg adult equivalent
    # so allometric scaling reads naturally; the paper itself reports values for
    # a 9.8 kg child on standard LPV/r 4:1 at steady state (Tikiso 2021 Table 4).
    # Conversion: CL_70 = CL_9.8 * (70/9.8)^0.75 and Vc_70 = Vc_9.8 * (70/9.8).
    # Allometric reference 70 kg matches the paper's Table 3 "Allometry CL/F
    # (L/h/70 kg) = 46.8" reported alongside the original 9.8 kg value of 10.7 L/h.
    lcl  <- log(46.77);  label("Apparent oral clearance CL/F at 70 kg, fully-mature CL, LPV/r 4:1 reference, well-nourished, liquid formulation, post-induction (L/h)")  # Tikiso 2021 Table 3 allometry column = 46.8 L/h/70 kg (= Table 4 10.7 L/h at 9.8 kg with exponent 0.75)
    lvc  <- log(78.57);  label("Apparent central volume of distribution Vc/F at 70 kg (L)")                                                                                  # Tikiso 2021 Table 4 Vc = 11.0 L at 9.8 kg, exponent 1; (70/9.8) = 7.143; 11.0 * 7.143 = 78.57 L
    lvp  <- log(23.79);  label("Apparent peripheral volume of distribution Vp/F at 70 kg (L)")                                                                               # Tikiso 2021 Table 4 Vp = 3.33 L at 9.8 kg, exponent 1; 3.33 * 7.143 = 23.79 L
    lq   <- log(4.808);  label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")                                                                                  # Tikiso 2021 Table 4 Q = 1.10 L/h at 9.8 kg, exponent 0.75; 1.10 * 4.371 = 4.808 L/h
    lka  <- log(2.29);   label("First-order absorption rate constant from depot to central, ka (1/h)")                                                                       # Tikiso 2021 Table 4 ka = 2.29 1/h (not weight-scaled)
    lmtt <- log(0.104);  label("Mean absorption transit-chain time MTT for the liquid formulation reference (h)")                                                            # Tikiso 2021 Table 4 MTT = 6.24 min = 0.104 h (liquid reference)
    lnn  <- log(11.9);   label("Number of absorption transit compartments NN (continuous, dimensionless)")                                                                   # Tikiso 2021 Table 4 NN = 11.9

    # Allometric exponents (paper-fixed: 0.75 on CL/Q, 1 on Vc/Vp)
    e_wt_cl_q  <- 0.75; label("Allometric exponent on CL and Q (unitless)")                  # Tikiso 2021 Section 2.4 (fixed at 0.75 on disposition clearances)
    e_wt_vc_vp <- 1;    label("Allometric exponent on Vc and Vp (unitless)")                 # Tikiso 2021 Section 2.4 (fixed at 1 on volumes of distribution)

    # Maturation function on CL (sigmoidal Hill on PMAGE; Tikiso 2021 Eq. 1)
    pmage50_mat   <- 8.10; label("Postmenstrual age at 50% mature CL (months)")              # Tikiso 2021 Table 4 PMAGE_50 = 8.10 months
    e_page_cl_hill <- 2.57; label("Hill coefficient (gamma) for the CL maturation function") # Tikiso 2021 Table 4 y_maturation = 2.57

    # Co-medication and formulation effects (multiplicative)
    e_efv_cl              <-  0.120; label("Efavirenz co-medication effect on CL relative to LPV/r 4:1 (fraction)")                          # Tikiso 2021 Table 4 "Change in CL when on EFV (%)" = +12.0%
    e_rif_lpvr4_fdepot    <- -0.294; label("Rifampicin + super-boosted LPV/r 4:4 effect on F relative to LPV/r 4:1 (fraction)")              # Tikiso 2021 Table 4 "Change in F on rifampicin + super-boosted lopinavir (%)" = -29.4%
    e_tablet_mtt          <-  0.249; label("FDC tablet formulation effect on MTT relative to liquid (fraction longer)")                      # Tikiso 2021 Table 4 "Change in speed of absorption for fixed-dose combination tablets (%)" = -24.9% (24.9% slower => 24.9% longer MTT)

    # Time-decaying malnutrition effects (Tikiso 2021 Eq. 2; gated by MAL_NOURISH)
    e_mal_fdepot <-  1.15; label("Malnutrition effect on F at start of nutritional supplementation (fraction)")  # Tikiso 2021 Table 4 "Change in F of malnourished children at start of supplementation (%)" = +115%
    e_mal_cl     <- -0.640; label("Malnutrition effect on CL at start of nutritional supplementation (fraction)") # Tikiso 2021 Table 4 "Change in CL of malnourished children at start of supplementation (%)" = -64.0%
    e_mal_thalf  <- 12.2;  label("Half-life of the malnutrition-effect decay (days)")                              # Tikiso 2021 Table 4 "Malnutrition effect half-life" = 12.2 days

    # Inter-individual variability (between-subject; omega^2 = log(1 + CV^2) for log-normal)
    etalcl ~ 0.02081  # Tikiso 2021 Table 4 BSV CL = 14.5% CV; omega^2 = log(1 + 0.145^2) = 0.02081
    etalvp ~ 0.18137  # Tikiso 2021 Table 4 BSV Vp = 44.6% CV; omega^2 = log(1 + 0.446^2) = 0.18137

    # Residual error (combined additive + proportional on linear ug/mL scale).
    # Tikiso 2021 dataset concentrations are in ug/mL (LLOQ 0.0238-0.0243 ug/mL);
    # the additive error of 2.01 ug/L = 0.00201 ug/mL is converted accordingly.
    propSd <- 0.238;    label("Proportional residual error (fraction)")  # Tikiso 2021 Table 4 "Proportional error" = 23.8%
    addSd  <- 0.00201;  label("Additive residual error (ug/mL)")          # Tikiso 2021 Table 4 "Additive error" = 2.01 ug/L = 0.00201 ug/mL
  })

  model({
    # 1. Maturation function for CL, sigmoidal Hill on postmenstrual age (Tikiso 2021 Eq. 1).
    #    Both PAGE and pmage50_mat are in months; the function returns 0 at PMAGE = 0,
    #    0.5 at PMAGE = pmage50_mat (8.10 months), and asymptotes to 1 at large PMAGE.
    #    For a typical term newborn (PAGE = 9 months), maturation is ~0.57; at PAGE = 19
    #    months (10 months postnatal age) it reaches ~0.90 (paper Section 3.2).
    maturation_cl <- (PAGE^e_page_cl_hill) / (pmage50_mat^e_page_cl_hill + PAGE^e_page_cl_hill)

    # 2. Time-decaying malnutrition effect (Tikiso 2021 Eq. 2). Gated by MAL_NOURISH so
    #    well-nourished subjects (MAL_NOURISH = 0) carry no effect. T_NUT_SUPP is in days.
    mal_decay <- MAL_NOURISH * exp(-T_NUT_SUPP * log(2) / e_mal_thalf)

    # 3. Multiplicative covariate factors on disposition and absorption.
    cl_efv  <- 1 + e_efv_cl              * CONMED_EFV
    cl_mal  <- 1 + e_mal_cl              * mal_decay
    f_rif   <- 1 + e_rif_lpvr4_fdepot    * CONMED_RIF_LPVR4
    f_mal   <- 1 + e_mal_fdepot          * mal_decay
    mtt_tab <- 1 + e_tablet_mtt          * FORM_TABLET

    # 4. Individual PK parameters with allometric scaling on disposition.
    #    Reference weight 70 kg; ka, MTT, NN, and e_*_cl factors are not weight-scaled.
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * maturation_cl * cl_efv * cl_mal
    vc  <- exp(lvc)          * (WT / 70)^e_wt_vc_vp
    vp  <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q   <- exp(lq)           * (WT / 70)^e_wt_cl_q
    ka  <- exp(lka)
    mtt <- exp(lmtt) * mtt_tab
    nn  <- exp(lnn)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. ODE system. Two-compartment disposition with a Savic 2007 analytical
    #    transit-compartment chain (rxode2 transit() built-in) feeding the depot,
    #    which then absorbs into central at rate ka. The transit() function reads
    #    the dose amount from podo(depot); f(depot) <- 0 suppresses the bolus
    #    contribution so the chain delivers the full dose. The bioavailability
    #    modifiers (RIF + super-boosted LPV/r and malnutrition on F) enter via
    #    the bio argument of transit() because nominal F is fixed at 1.
    f_total <- f_rif * f_mal
    d/dt(depot)        <- transit(nn, mtt, f_total) - ka * depot
    d/dt(central)      <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <- k12 * central - k21 * peripheral1

    f(depot) <- 0  # suppress bolus into depot; transit() drives the input rate

    # 6. Plasma abacavir concentration (dose mg / volume L = mg/L = ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
