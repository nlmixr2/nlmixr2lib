deKock_2017_sulfadoxinePyrimethamine <- function() {
  description <- paste(
    "Joint popPK model for the antimalarial fixed-dose combination of",
    "sulfadoxine (1500 mg) and pyrimethamine (75 mg) as intermittent",
    "preventive treatment during pregnancy (IPTp) and after delivery in",
    "98 women from Mali, Mozambique, Sudan, and Zambia (de Kock 2017).",
    "Sulfadoxine has 2-compartment disposition with first-order",
    "absorption; pyrimethamine has 3-compartment disposition with",
    "first-order absorption. Apparent volumes and flow rates are",
    "allometrically scaled with total body weight (exponents 1 and",
    "0.75 respectively, reference WT = 60 kg). Whole-blood predictions",
    "are derived from plasma predictions using hematocrit and an",
    "estimated RBC-to-plasma partition ratio per drug. Pregnancy",
    "effects on apparent CL differ by drug: sulfadoxine uses a",
    "sigmoidal time-after-delivery effect (asymptotic -75.7%, T50 =",
    "6.35 weeks, gamma = 4.90), while pyrimethamine uses a step",
    "contrast (+21.2% postpartum). Pyrimethamine apparent CL is",
    "additionally -20.2% in the Mozambique site. Residual",
    "country-specific scaling on the observed whole-blood",
    "concentrations is fitted with Mali as the reference."
  )
  reference <- paste(
    "de Kock M, Tarning J, Workman L, Nyunt MM, Adam I, Barnes KI,",
    "Denti P.",
    "Pharmacokinetics of Sulfadoxine and Pyrimethamine for Intermittent",
    "Preventive Treatment of Malaria During Pregnancy and After Delivery.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(7):430-438.",
    "doi:10.1002/psp4.12181.",
    sep = " "
  )
  vignette <- "deKock_2017_sulfadoxinePyrimethamine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L for sulfadoxine (whole blood); ng/mL for pyrimethamine (whole blood)"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling at the population median 60 kg applied to all",
        "apparent volumes (Vc, Vp, Vp2) with exponent 1.0 and to all",
        "apparent flows (CL, Q, Q2) with exponent 0.75, for both",
        "sulfadoxine and pyrimethamine. Height was not collected so",
        "fat-free mass was not testable."
      ),
      source_name        = "WT"
    ),
    HCT = list(
      description        = "Hematocrit (volume fraction)",
      units              = "fraction (0-1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-visit hematocrit expressed as a fraction of 1. The source",
        "data column was baseline hemoglobin (HGB in g/dL); the paper",
        "converts to hematocrit using the Lee 2008 formula for",
        "Plasmodium-falciparum-infected subjects (HCT = 0.0375 * HGB +",
        "0.0079 + ... , Eq. 1 of the paper). The conversion is treated",
        "as a data-assembly step; the model expects HCT directly. Used",
        "to scale model-predicted plasma concentrations to observed",
        "whole-blood concentrations via",
        "CWB = CPL * (HCT * hRBC_PL + (1 - HCT)) where hRBC_PL is the",
        "drug-specific RBC-to-plasma partition ratio estimated in",
        "ini()."
      ),
      source_name        = "HCT"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnancy visit (second or third trimester); 0 = postpartum",
        "visit. Encoded per occasion. Pyrimethamine apparent CL is",
        "modelled as a step contrast on PREG so that postpartum CL is",
        "21.2% higher than the pregnancy CL (Table 2 row \"Change in CL",
        "when non-pregnant\"). For sulfadoxine the pregnancy effect is",
        "carried instead by the continuous TPP covariate so that the",
        "model captures the gradual postpartum recovery; PREG itself",
        "does not enter the sulfadoxine CL equation."
      ),
      source_name        = "PREG"
    ),
    TPP = list(
      description        = "Time after delivery",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Weeks since delivery for postpartum visits. Set to 0 for",
        "pregnancy visits so the sigmoidal postpartum-recovery effect",
        "on sulfadoxine CL evaluates to 0 during pregnancy. The",
        "sigmoid (TPP^gamma) / (TPP^gamma + T50^gamma) with T50 = 6.35",
        "weeks and gamma = 4.90 reaches its asymptote near 13 weeks",
        "postpartum (Figure 3 of the paper)."
      ),
      source_name        = "TPP"
    ),
    REGION_MOZAMBIQUE = list(
      description        = "Mozambique study-site indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Mozambique; Mali is the reference site)",
      notes              = paste(
        "1 = Mozambique site (31 pregnant / 22 postpartum subjects);",
        "0 otherwise. Multiplicatively reduces pyrimethamine apparent",
        "CL by 20.2% (Table 2 \"Difference in clearance in",
        "Mozambique\"). Also scales observed sulfadoxine concentrations",
        "by +21.2% and observed pyrimethamine concentrations by +57.6%",
        "(Table 2 \"Site effect (scaling on observations)\")."
      ),
      source_name        = "SITE = 'Mozambique'"
    ),
    REGION_SUDAN = list(
      description        = "Sudan study-site indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Sudan; Mali is the reference site)",
      notes              = paste(
        "1 = Sudan site (24 pregnant / 9 postpartum subjects); 0",
        "otherwise. Scales observed sulfadoxine concentrations by",
        "+15.5% and observed pyrimethamine concentrations by +33.2%",
        "(Table 2 \"Site effect (scaling on observations)\")."
      ),
      source_name        = "SITE = 'Sudan'"
    ),
    REGION_ZAMBIA = list(
      description        = "Zambia study-site indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Zambia; Mali is the reference site)",
      notes              = paste(
        "1 = Zambia site (25 pregnant / 18 postpartum subjects); 0",
        "otherwise. Scales observed sulfadoxine concentrations by",
        "-24.8% and observed pyrimethamine concentrations by -5.40%",
        "(Table 2 \"Site effect (scaling on observations)\")."
      ),
      source_name        = "SITE = 'Zambia'"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 98,
    n_subjects_postpartum = 77,
    n_studies      = 1,
    n_sites        = 4,
    age_range      = "Adult women of reproductive age (per-site median 24-31 years; Table 1)",
    weight_range   = "Adult women (population median 60 kg, the allometric reference; per-site median 60-66 kg, Table 1)",
    sex_female_pct = 100,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Pregnant women in the second or third trimester (median",
      "gestational age 27-28 weeks across sites) and the same women",
      "after delivery (postpartum sampling 6.4 to 46 weeks after",
      "delivery, depending on site); no concurrent Plasmodium",
      "falciparum infection required for enrolment."
    ),
    dose_range     = paste(
      "Single oral fixed-dose combination tablet containing 1,500 mg",
      "sulfadoxine and 75 mg pyrimethamine, administered as part of",
      "intermittent preventive treatment of malaria in pregnancy",
      "(IPTp). One pregnancy dose plus one matched ad-hoc postpartum",
      "dose per subject."
    ),
    regions        = "Sub-Saharan Africa: Mali, Mozambique, Sudan, Zambia",
    notes          = paste(
      "Demographics from de Kock 2017 Table 1 (per-site medians of",
      "age, weight, hemoglobin, and gestational age / time after",
      "delivery). Whole-blood samples collected from a single dried",
      "blood spot per timepoint; the entire blood spot was cut out",
      "and analysed to avoid hematocrit-dependent punch effects. Drug",
      "concentrations were measured in capillary whole blood; the",
      "model is parameterised in plasma and predicted whole-blood",
      "concentrations are derived via the hematocrit and the",
      "estimated RBC-to-plasma partition ratio."
    )
  )

  ini({
    # ============================================================
    # Sulfadoxine structural parameters
    # ----- de Kock 2017 Table 2 "Sulfadoxine" column -----
    # All apparent volumes and flow rates are reported referring to
    # plasma and after allometric scaling to the median body weight
    # 60 kg (Table 2 footnote b).
    # ============================================================
    lka     <- log(0.531)
    label("Sulfadoxine first-order absorption rate, ka (1/h)")              # Table 2 sulfa: ka = 0.531 /h
    lcl     <- log(0.0303)
    label("Sulfadoxine apparent CL during pregnancy at WT = 60 kg (L/h)")   # Table 2 sulfa: CL/F (during pregnancy) = 0.0303 L/h
    lvc     <- log(14.1)
    label("Sulfadoxine apparent central volume of distribution at WT = 60 kg (L)")  # Table 2 sulfa: Vc/F = 14.1 L
    lq      <- log(0.0252)
    label("Sulfadoxine apparent inter-compartmental CL at WT = 60 kg (L/h)") # Table 2 sulfa: Qp1/F = 0.0252 L/h
    lvp     <- log(179)
    label("Sulfadoxine apparent peripheral volume of distribution at WT = 60 kg (L)") # Table 2 sulfa: Vp1/F = 179 L
    lfdepot <- fixed(log(1))
    label("Sulfadoxine relative bioavailability F (unitless, FIXED at 1)")  # Table 2 sulfa: F = 1 FIXED

    # ============================================================
    # Pyrimethamine structural parameters
    # ----- de Kock 2017 Table 2 "Pyrimethamine" column -----
    # ============================================================
    lka_pyra     <- log(1.31)
    label("Pyrimethamine first-order absorption rate, ka (1/h)")                 # Table 2 pyra: ka = 1.31 /h
    lcl_pyra     <- log(1.35)
    label("Pyrimethamine apparent CL during pregnancy at WT = 60 kg (L/h)")      # Table 2 pyra: CL/F (during pregnancy) = 1.35 L/h
    lvc_pyra     <- log(163)
    label("Pyrimethamine apparent central volume of distribution at WT = 60 kg (L)") # Table 2 pyra: Vc/F = 163 L
    lq_pyra      <- log(1.45)
    label("Pyrimethamine apparent shallow Q at WT = 60 kg (L/h)")                # Table 2 pyra: Qp1/F = 1.45 L/h
    lvp_pyra     <- log(29.8)
    label("Pyrimethamine apparent shallow peripheral volume at WT = 60 kg (L)")  # Table 2 pyra: Vp1/F = 29.8 L
    lq2_pyra     <- log(0.122)
    label("Pyrimethamine apparent deep Q at WT = 60 kg (L/h)")                   # Table 2 pyra: Qp2/F = 0.122 L/h
    lvp2_pyra    <- log(251)
    label("Pyrimethamine apparent deep peripheral volume at WT = 60 kg (L)")     # Table 2 pyra: Vp2/F = 251 L
    lfdepot_pyra <- fixed(log(1))
    label("Pyrimethamine relative bioavailability F (unitless, FIXED at 1)")     # Table 2 pyra: F = 1 FIXED

    # ============================================================
    # Whole-blood / plasma partition (per drug)
    # Estimated under a log-normal Bayesian prior (typical values
    # 0.16 sulfa / 0.42 pyra with 30% uncertainty, Methods).
    # ============================================================
    hrbcpl      <- 0.155
    label("Sulfadoxine RBC-to-plasma ratio (fraction)")                     # Table 2 sulfa: hRBC/PL = 0.155
    hrbcpl_pyra <- 0.324
    label("Pyrimethamine RBC-to-plasma ratio (fraction)")                   # Table 2 pyra: hRBC/PL = 0.324

    # ============================================================
    # Sulfadoxine postpartum-recovery sigmoid on CL
    # CL_sulfa(TPP) = exp(lcl) * (1 + e_tpp_cl * sigmoid(TPP))
    # where sigmoid(t) = t^gamma / (t^gamma + t50^gamma).
    # At TPP = 0 (pregnancy) the sigmoid is 0; far postpartum the
    # sigmoid approaches 1 and CL is reduced by 75.7%.
    # ============================================================
    e_tpp_cl  <- -0.757
    label("Sulfadoxine asymptotic postpartum fractional change in CL (unitless)") # Table 2 sulfa: Change in CL when non-pregnant = -75.7%
    t50_tpp   <- 6.35
    label("Sulfadoxine postpartum recovery T50 (weeks)")                          # Table 2 sulfa: T50 = 6.35 weeks
    gamma_tpp <- 4.90
    label("Sulfadoxine postpartum recovery sigmoid shape parameter (unitless)")   # Table 2 sulfa: gamma = 4.90

    # ============================================================
    # Pyrimethamine pregnancy step effect on CL
    # CL_pyra = exp(lcl_pyra) * (1 + e_preg_cl_pyra * (1 - PREG))
    # so pregnancy (PREG = 1) reproduces the structural CL and
    # postpartum (PREG = 0) raises CL by 21.2%.
    # ============================================================
    e_preg_cl_pyra <- 0.212
    label("Pyrimethamine postpartum fractional change in CL applied via (1 - PREG)") # Table 2 pyra: Change in CL when non-pregnant = +21.2%

    # Pyrimethamine Mozambique site effect on CL.
    e_region_mozambique_cl_pyra <- -0.202
    label("Pyrimethamine multiplicative CL effect for Mozambique site (unitless)")   # Table 2 pyra: Difference in clearance in Mozambique = -20.2%

    # ============================================================
    # Site-specific multiplicative scaling on observed whole-blood
    # concentrations (Table 2 "Site effect (scaling on
    # observations)"). Mali is the implicit reference (all three
    # indicators = 0).
    # ============================================================
    e_region_mozambique_cc <- 0.212
    label("Sulfadoxine observation-scaling effect, Mozambique (fraction)")     # Table 2 sulfa: site effect on observations, Mozambique = +21.2%
    e_region_sudan_cc      <- 0.155
    label("Sulfadoxine observation-scaling effect, Sudan (fraction)")          # Table 2 sulfa: site effect on observations, Sudan = +15.5%
    e_region_zambia_cc     <- -0.248
    label("Sulfadoxine observation-scaling effect, Zambia (fraction)")         # Table 2 sulfa: site effect on observations, Zambia = -24.8%

    e_region_mozambique_cc_pyra <- 0.576
    label("Pyrimethamine observation-scaling effect, Mozambique (fraction)")   # Table 2 pyra: site effect on observations, Mozambique = +57.6%
    e_region_sudan_cc_pyra      <- 0.332
    label("Pyrimethamine observation-scaling effect, Sudan (fraction)")        # Table 2 pyra: site effect on observations, Sudan = +33.2%
    e_region_zambia_cc_pyra     <- -0.054
    label("Pyrimethamine observation-scaling effect, Zambia (fraction)")       # Table 2 pyra: site effect on observations, Zambia = -5.40%

    # ============================================================
    # IIV - log-normal; omega^2 = log(1 + CV^2). Reported CV% from
    # Table 2 Final column. The paper distinguishes between-subject
    # variability (BSV) and between-occasion variability (BOV); for
    # nlmixr2lib forward-simulation use each is carried as a
    # subject-level eta and the BOV provenance is documented inline.
    # ============================================================
    etalcl       ~ log(1 + 0.313^2)
    # Table 2 sulfa: BSV in CL = 31.3% CV  ->  omega^2 = log(1 + 0.313^2)
    etalka       ~ log(1 + 0.564^2)
    # Table 2 sulfa: BOV in ka = 56.4% CV  ->  omega^2 = log(1 + 0.564^2); the paper's BOV is carried here as a subject-level eta
    etalcl_pyra  ~ log(1 + 0.123^2)
    # Table 2 pyra: BSV in CL = 12.3% CV  ->  omega^2 = log(1 + 0.123^2)

    # Correlated BOV on F across the two drugs.
    # Table 2 reports BOV(F_sulfa) = 20.7% CV, BOV(F_pyra) = 17.6%
    # CV, and corr(F_sulfa, F_pyra) = 0.677. Encoded as a 2x2 block
    # using the log-normal variances and the corresponding
    # covariance: cov = rho * sqrt(var_sulfa * var_pyra).
    etalfdepot + etalfdepot_pyra ~ c(
      log(1 + 0.207^2),
      0.677 * sqrt(log(1 + 0.207^2) * log(1 + 0.176^2)),
      log(1 + 0.176^2)
    )
    # Table 2: BOV in F sulfa = 20.7% CV, BOV in F pyra = 17.6% CV; correlation in bioavailability across the two drugs = 67.9% (Table 2)

    # ============================================================
    # Residual error - combined additive + proportional per drug,
    # on whole-blood concentrations.
    # ============================================================
    propSd      <- 0.170
    label("Sulfadoxine proportional residual SD (fraction)")                  # Table 2 sulfa: proportional error = 17.0%
    addSd       <- 2.13
    label("Sulfadoxine additive residual SD (mg/L = ug/mL whole blood)")      # Table 2 sulfa: additive error = 2.13 ug/mL
    propSd_pyra <- 0.180
    label("Pyrimethamine proportional residual SD (fraction)")                # Table 2 pyra: proportional error = 18.0%
    addSd_pyra  <- 2.45
    label("Pyrimethamine additive residual SD (ng/mL whole blood)")           # Table 2 pyra: additive error = 2.45 ng/mL
  })

  model({
    # ------------------------------------------------------------
    # Allometric ratio (reference WT = 60 kg, the population
    # median; Table 2 footnote b). Volumes scale linearly with WT;
    # clearance and inter-compartmental flow scale with WT^0.75.
    # ------------------------------------------------------------
    wt_ratio <- WT / 60

    # ------------------------------------------------------------
    # Sigmoidal postpartum-recovery factor for sulfadoxine CL.
    # During pregnancy TPP = 0 so the sigmoid evaluates to 0 and
    # the structural pregnancy-state CL is recovered.
    # ------------------------------------------------------------
    tpp_pow   <- TPP^gamma_tpp
    sigmoid_tpp <- tpp_pow / (tpp_pow + t50_tpp^gamma_tpp)

    # ------------------------------------------------------------
    # Sulfadoxine individual PK parameters.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wt_ratio^0.75 *
          (1 + e_tpp_cl * sigmoid_tpp)
    vc <- exp(lvc)                  * wt_ratio
    q  <- exp(lq)                   * wt_ratio^0.75
    vp <- exp(lvp)                  * wt_ratio
    ka <- exp(lka + etalka)

    # ------------------------------------------------------------
    # Pyrimethamine individual PK parameters.
    # ------------------------------------------------------------
    cl_pyra  <- exp(lcl_pyra + etalcl_pyra) * wt_ratio^0.75 *
                (1 + e_preg_cl_pyra * (1 - PREG)) *
                (1 + e_region_mozambique_cl_pyra * REGION_MOZAMBIQUE)
    vc_pyra  <- exp(lvc_pyra)               * wt_ratio
    q_pyra   <- exp(lq_pyra)                * wt_ratio^0.75
    vp_pyra  <- exp(lvp_pyra)               * wt_ratio
    q2_pyra  <- exp(lq2_pyra)               * wt_ratio^0.75
    vp2_pyra <- exp(lvp2_pyra)              * wt_ratio
    ka_pyra  <- exp(lka_pyra)

    # ------------------------------------------------------------
    # Sulfadoxine 2-compartment disposition with first-order
    # absorption from a `depot` compartment.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl + q) * central / vc +
                          q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp

    # Bioavailability on the sulfadoxine depot. F = 1 with
    # log-normal BOV in F (correlated across drugs via the
    # etalfdepot + etalfdepot_pyra block).
    f(depot) <- exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------
    # Pyrimethamine 3-compartment disposition with first-order
    # absorption from a separate `depot_pyra` compartment.
    # ------------------------------------------------------------
    d/dt(depot_pyra)        <- -ka_pyra * depot_pyra
    d/dt(central_pyra)      <-  ka_pyra * depot_pyra -
                                (cl_pyra + q_pyra + q2_pyra) * central_pyra / vc_pyra +
                                q_pyra  * peripheral1_pyra / vp_pyra +
                                q2_pyra * peripheral2_pyra / vp2_pyra
    d/dt(peripheral1_pyra)  <- q_pyra  * central_pyra / vc_pyra -
                               q_pyra  * peripheral1_pyra / vp_pyra
    d/dt(peripheral2_pyra)  <- q2_pyra * central_pyra / vc_pyra -
                               q2_pyra * peripheral2_pyra / vp2_pyra

    f(depot_pyra) <- exp(lfdepot_pyra + etalfdepot_pyra)

    # ------------------------------------------------------------
    # Whole-blood predictions. The model is parameterised on
    # plasma; capillary whole-blood concentrations are obtained
    # via the standard partition formula
    #   CWB = CPL * (HCT * hRBC_PL + (1 - HCT))
    # using the drug-specific RBC-to-plasma ratio. Site-specific
    # multiplicative scaling on the observed concentration (Mali =
    # reference) captures residual between-site differences in
    # apparent bioavailability and dried-blood-spot sample
    # handling (Discussion).
    # ------------------------------------------------------------
    wb_factor      <- HCT * hrbcpl      + (1 - HCT)
    wb_factor_pyra <- HCT * hrbcpl_pyra + (1 - HCT)

    site_factor <- 1 +
      e_region_mozambique_cc * REGION_MOZAMBIQUE +
      e_region_sudan_cc      * REGION_SUDAN +
      e_region_zambia_cc     * REGION_ZAMBIA
    site_factor_pyra <- 1 +
      e_region_mozambique_cc_pyra * REGION_MOZAMBIQUE +
      e_region_sudan_cc_pyra      * REGION_SUDAN +
      e_region_zambia_cc_pyra     * REGION_ZAMBIA

    # Sulfadoxine whole-blood concentration in mg/L (equivalently
    # ug/mL); pyrimethamine whole-blood concentration in ng/mL via
    # the 1000-fold molar-mass-independent unit conversion from
    # mg/L of central / vc_pyra.
    Cc      <- (central      / vc)      * wb_factor      * site_factor
    Cc_pyra <- (central_pyra / vc_pyra) * 1000 * wb_factor_pyra * site_factor_pyra

    Cc      ~ add(addSd)      + prop(propSd)
    Cc_pyra ~ add(addSd_pyra) + prop(propSd_pyra)
  })
}
