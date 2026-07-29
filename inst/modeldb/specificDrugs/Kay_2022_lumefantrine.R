# Population pharmacokinetic model for oral lumefantrine in Ugandan
# children with and without HIV receiving artemether-lumefantrine
# (Coartem Dispersible) with concomitant antiretroviral therapy (Kay
# 2022, Clin Pharmacol Ther 113(3):660-669; doi:10.1002/cpt.2768).

Kay_2022_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine in 277 Ugandan children",
    "(186 HIV-uninfected, 178 HIV-infected on efavirenz-, nevirapine-,",
    "or lopinavir/ritonavir-based antiretroviral therapy plus daily",
    "trimethoprim-sulfamethoxazole prophylaxis) ages ~2 months to 8.6",
    "years treated with six-dose weight-based Coartem Dispersible (20 mg",
    "artemether + 120 mg lumefantrine per tablet) for uncomplicated",
    "Plasmodium falciparum malaria (Kay 2022). Two-compartment",
    "disposition with first-order absorption. Body-weight fixed effects",
    "scale all clearance and volume terms with a reference weight of",
    "15 kg; volumes use an allometric exponent of 1, clearances use an",
    "age-dependent exponent (1.2 for age <= 3 months, 1.0 for >3 to 24",
    "months, 0.9 for >24 to 60 months, 0.75 for >60 months) from",
    "Anderson & Holford 2009 (paper ref 34). A power-form age effect on",
    "relative bioavailability captures reduced lumefantrine",
    "bioavailability in young children (F = (age_months / 50)^0.204).",
    "Concomitant antiretroviral therapy enters as mutually-exclusive",
    "linear-deviation effects on apparent oral CL/F and ka: efavirenz",
    "increases CL/F by 98.2% and ka by 48.4%, lopinavir/ritonavir",
    "decreases CL/F by 51.4% and ka by 21.2%, nevirapine has no",
    "statistically significant effect (both CIs cross zero). Q/F IIV",
    "is fixed at 15.9% CV; the remaining four structural parameters",
    "carry estimated log-normal IIV (104%, 112%, 127%, and 16.9% CV for",
    "CL/F, V2/F, V3/F, and ka respectively). NONMEM proportional",
    "residual error on linear concentration (sigma^2 = 0.200, 44.7%",
    "CV).",
    sep = " "
  )
  reference <- paste(
    "Kay K, Goodwin J, Ehrlich H, Ou J, Freeman T, Wang K, Li F, Wade M,",
    "French J, Huang L, Aweeka F, Mwebaza N, Kajubi R, Riggs M,",
    "Ruiz-Garcia A, Parikh S (2023).",
    "Impact of Drug Exposure on Resistance Selection Following",
    "Artemether-Lumefantrine Treatment for Malaria in Children With and",
    "Without HIV in Uganda.",
    "Clinical Pharmacology & Therapeutics 113(3):660-669.",
    "doi:10.1002/cpt.2768.",
    sep = " "
  )
  vignette <- "Kay_2022_lumefantrine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at episode enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with reference weight 15 kg (approximate",
        "cohort median per Table 2 footnote: '50-month-old",
        "HIV-uninfected child weighing 15 kg'). Volumes use exponent 1;",
        "clearances use the age-dependent exponent (Anderson & Holford",
        "2009 / paper ref 34): 1.2 for age <= 3 months, 1.0 for >3 to",
        "24 months, 0.9 for >24 to 60 months, 0.75 for >60 months.",
        "Cohort weight range 7.65-30.0 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at episode enrolment",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives (a) the piecewise allometric exponent on apparent CL/F",
        "and Q/F (see WT covariate notes) and (b) the power-form age",
        "effect on relative bioavailability F: F = (age_months / 50)^",
        "0.204 with reference 50 months (= 4.17 years; Table 2 footnote",
        "'reference subject ... 50-month-old HIV-uninfected child').",
        "Age in years is converted inside model() to months via",
        "agemo <- AGE * 12. Cohort age range ~0.16-8.58 years (Table",
        "1; HIV-uninfected median 3.58 years, HIV-infected medians",
        "4.5-6.0 years across EFV / LPV/r / NVP arms)."
      ),
      source_name        = "AGE"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz-based antiretroviral therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child on an efavirenz-based ART regimen",
        "(daily oral EFV-containing combination ART, plus daily",
        "trimethoprim-sulfamethoxazole prophylaxis); 0 = HIV-uninfected",
        "child (the model reference) or HIV-infected child on a non-EFV",
        "ART regimen. Time-fixed per episode in the Kay 2022 cohort.",
        "Efavirenz is a CYP3A4 inducer; the indicator enters as a",
        "linear-deviation multiplier on apparent oral lumefantrine CL/F",
        "and on the first-order absorption rate constant ka:",
        "CL/F = TVCL * (1 + e_efv_cl * CONMED_EFV) and",
        "ka  = TVKA * (1 + e_efv_ka * CONMED_EFV) with e_efv_cl =",
        "+0.982 and e_efv_ka = +0.484 (Table 2 theta_7 EFV_CL/F and",
        "theta_10 EFV_KA; Discussion: 'Patients receiving EFV had CL/F",
        "and KA estimates ~198.2% and 148% those estimated in",
        "HIV-uninfected patients, respectively'). Mutually exclusive",
        "with CONMED_LPV and CONMED_NVP within the source cohort."
      ),
      source_name        = "EFV"
    ),
    CONMED_LPV = list(
      description        = "Concomitant lopinavir/ritonavir-based antiretroviral therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child on a lopinavir/ritonavir-based ART",
        "regimen (daily oral LPV/r-containing combination ART, plus",
        "daily trimethoprim-sulfamethoxazole prophylaxis); 0 =",
        "HIV-uninfected child or HIV-infected child on a non-LPV/r ART",
        "regimen. Time-fixed per episode. Lopinavir/ritonavir is a",
        "potent CYP3A4 inhibitor (the ritonavir booster being the",
        "principal perpetrator); the indicator enters as a",
        "linear-deviation multiplier on apparent oral lumefantrine CL/F",
        "and on the first-order absorption rate constant ka:",
        "CL/F = TVCL * (1 + e_lpv_cl * CONMED_LPV) and",
        "ka  = TVKA * (1 + e_lpv_ka * CONMED_LPV) with e_lpv_cl =",
        "-0.514 and e_lpv_ka = -0.212 (Table 2 theta_8 LPV/r_CL/F and",
        "theta_11 LPV/r_KA; Discussion: 'Patients receiving LPV/r had",
        "CL/F and KA estimates ~48.6% and 78.8% that of HIV-uninfected",
        "patients, respectively'). Mutually exclusive with CONMED_EFV",
        "and CONMED_NVP within the source cohort."
      ),
      source_name        = "LPV/r"
    ),
    CONMED_NVP = list(
      description        = "Concomitant nevirapine-based antiretroviral therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child on a nevirapine-based ART regimen",
        "(daily oral NVP-containing combination ART, plus daily",
        "trimethoprim-sulfamethoxazole prophylaxis); 0 = HIV-uninfected",
        "child or HIV-infected child on a non-NVP ART regimen.",
        "Time-fixed per episode. Nevirapine is a milder CYP3A4",
        "inducer than efavirenz; the Kay 2022 model retains an effect",
        "on apparent oral lumefantrine CL/F and ka but neither effect",
        "is statistically significant (both 95% CIs cross zero):",
        "CL/F = TVCL * (1 + e_nvp_cl * CONMED_NVP) and",
        "ka  = TVKA * (1 + e_nvp_ka * CONMED_NVP) with e_nvp_cl =",
        "+0.0191 (95% CI -0.324, 0.362) and e_nvp_ka = -0.0589 (95% CI",
        "-0.207, 0.0891) (Table 2 theta_9 NVP_CL/F and theta_12",
        "NVP_KA). The text reports the practical interpretation 'no",
        "statistically significant effect of NVP treatment on",
        "lumefantrine CL/F and KA ... equivalent with nevirapine and",
        "HIV-uninfected children' (Results). Mutually exclusive with",
        "CONMED_EFV and CONMED_LPV within the source cohort."
      ),
      source_name        = "NVP"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 277L,
    n_episodes           = 364L,
    n_uninfected_episodes = 186L,
    n_efv_episodes       = 48L,
    n_lpv_episodes       = 68L,
    n_nvp_episodes       = 62L,
    n_studies            = 1L,
    age_range            = "0.16-8.58 years (HIV-uninfected median 3.58, EFV median 6.00, LPV/r median 4.50, NVP median 4.50; Table 1)",
    weight_range         = "7.65-30.0 kg (HIV-uninfected median 14.1, EFV median 18.0, LPV/r median 15.4, NVP median 16.0; Table 1)",
    sex_female_pct       = NA_real_,
    disease_state        = paste(
      "Uncomplicated Plasmodium falciparum malaria in high-transmission",
      "Tororo, Uganda. HIV-uninfected children and HIV-infected children",
      "on continuous efavirenz-, nevirapine-, or lopinavir/ritonavir-",
      "based combination antiretroviral therapy with daily",
      "trimethoprim-sulfamethoxazole prophylaxis (96% of HIV-infected",
      "children). Children allowed to re-enrol for repeat clinical",
      "episodes within and beyond the 42-day follow-up window."
    ),
    dose_range           = paste(
      "Coartem Dispersible (Novartis Pharma AG, Basel): 20 mg artemether",
      "+ 120 mg lumefantrine per tablet, weight-based dosing per WHO",
      "guidelines (1 tablet per dose for 5-14 kg, 2 tablets for 15-24",
      "kg, 3 tablets for 25-34 kg). Six oral doses (oral) at 0, 8, 24,",
      "36, 48, and 60 hours, taken with milk or during breastfeeding to",
      "promote lumefantrine absorption."
    ),
    regions              = "Uganda (Tororo, high malaria transmission setting, 2011-2014)",
    notes                = paste(
      "Per-episode sex breakdown from Table 1: % male episodes were",
      "53.2% (HIV-uninfected), 33.3% (EFV), 35.3% (LPV/r), and 53.2%",
      "(NVP); aggregate sex_female_pct is therefore approximately",
      "100 - weighted-male = ~54% but is not reported as a single",
      "value and is left NA. Parent NCA paper (used for the previous",
      "noncompartmental analysis): Parikh et al. 2016 (paper ref 29).",
      "Intensive PK sampling cohort: pre-first-dose (day 0) and",
      "pre/post-sixth dose (7 venous samples on day 3 at 0, 0.5, 1, 2,",
      "3, 4, 8 h post-last dose); capillary sampling on days 4, 7, 14,",
      "21. Sparse PK cohort: capillary sampling on days 7, 14, 21",
      "only. Lumefantrine quantified by LC-MS/MS with LLOQ 50 ng/mL.",
      "Companion downstream PK/PD time-to-event recurrence models",
      "(Table 3) and per-arm post-treatment chemoprophylaxis durations",
      "(Supplementary Table S4) are described in the vignette but not",
      "encoded here -- see vignette Assumptions and deviations."
    )
  )

  ini({
    # Structural population mean parameters from Kay 2022 Table 2,
    # "Structural model parameters" rows. Reference subject: 50-month-old
    # (4.17-year-old) HIV-uninfected child weighing 15 kg (Table 2 footnote).
    # Parameters are apparent (relative to F = 1) and reported on the
    # linear scale; log() is applied here for the nlmixr2 internal log scale.

    lcl <- log(1.20)
    label("Apparent lumefantrine clearance CL/F (L/h)")
    # Kay 2022 Table 2 theta_1: CL/F = 1.20 L/h (95% CI 0.952, 1.45)

    lvc <- log(24.1)
    label("Apparent lumefantrine central volume of distribution V2/F (L)")
    # Kay 2022 Table 2 theta_2: V2/F = 24.1 L (95% CI 19.0, 29.2)

    lq <- log(0.380)
    label("Apparent lumefantrine intercompartmental clearance Q/F (L/h)")
    # Kay 2022 Table 2 theta_3: Q/F = 0.380 L/h (95% CI 0.258, 0.501)

    lvp <- log(767)
    label("Apparent lumefantrine peripheral volume of distribution V3/F (L)")
    # Kay 2022 Table 2 theta_4: V3/F = 767 L (95% CI 174, 1.36e+03)

    lka <- log(0.0215)
    label("First-order absorption rate constant ka (1/h)")
    # Kay 2022 Table 2 theta_5: KA = 0.0215 1/h (95% CI 0.0197, 0.0234)

    # Bioavailability anchor: fixed at 1 because the model only identifies
    # apparent (relative-to-F) quantities. The age effect on F is the
    # estimated theta_6 below; without F = 1 fixed, the model would be
    # non-identifiable.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless; structural anchor, fixed at 1)")
    # Kay 2022 Methods: structural anchor for apparent-PK parameterisation.

    # Power-form age effect on relative bioavailability:
    # F = (AGE_months / 50_months)^e_age_f with reference age 50 months
    # (Table 2 footnote: 'reference subject ... 50-month-old HIV-uninfected
    # child'). Discussion ('estimated exponent of the age effect on
    # bioavailability was 0.204') makes the power-form parameterisation
    # explicit. The Discussion's verbal anchor 'compared with 5-year-old
    # children' uses a different reference (60 months) than the model's
    # 50-month reference; the model's parameter applies relative to 50 months.
    e_age_f <- 0.204
    label("Power-form age exponent on relative bioavailability F (unitless)")
    # Kay 2022 Table 2 theta_6: AGE_F = 0.204 (95% CI -0.0586, 0.467)

    # ART effects on apparent oral lumefantrine clearance. The Kay 2022
    # parameterisation is the linear-deviation form CL/F = TVCL * (1 + theta
    # * CONMED), which exactly reproduces the Discussion's reported
    # multipliers: EFV ~198.2% (= 1 + 0.982), LPV/r ~48.6% (= 1 - 0.514),
    # NVP ~102% (= 1 + 0.019; statistically not significant).
    e_efv_cl <- 0.982
    label("Linear-deviation effect of efavirenz on lumefantrine CL/F (fractional)")
    # Kay 2022 Table 2 theta_7: EFV_CL/F = 0.982 (95% CI 0.163, 1.80)

    e_lpv_cl <- -0.514
    label("Linear-deviation effect of lopinavir/ritonavir on lumefantrine CL/F (fractional)")
    # Kay 2022 Table 2 theta_8: LPV/r_CL/F = -0.514 (95% CI -0.696, -0.332)

    e_nvp_cl <- 0.0191
    label("Linear-deviation effect of nevirapine on lumefantrine CL/F (fractional; ns)")
    # Kay 2022 Table 2 theta_9: NVP_CL/F = 0.0191 (95% CI -0.324, 0.362; not significant)

    # ART effects on first-order absorption rate constant ka. Same
    # linear-deviation form. Discussion: EFV ~148% (= 1 + 0.484), LPV/r
    # ~78.8% (= 1 - 0.212), NVP ~94% (= 1 - 0.059; not significant).
    e_efv_ka <- 0.484
    label("Linear-deviation effect of efavirenz on absorption rate ka (fractional)")
    # Kay 2022 Table 2 theta_10: EFV_KA = 0.484 (95% CI 0.282, 0.685)

    e_lpv_ka <- -0.212
    label("Linear-deviation effect of lopinavir/ritonavir on absorption rate ka (fractional)")
    # Kay 2022 Table 2 theta_11: LPV/r_KA = -0.212 (95% CI -0.305, -0.120)

    e_nvp_ka <- -0.0589
    label("Linear-deviation effect of nevirapine on absorption rate ka (fractional; ns)")
    # Kay 2022 Table 2 theta_12: NVP_KA = -0.0589 (95% CI -0.207, 0.0891; not significant)

    # IIV. Kay 2022 Table 2 footnote: "CV% of omegas = sqrt(exp(estimate)
    # - 1) * 100" -- i.e., the tabled Estimate is the log-normal variance
    # (omega^2) on the internal log scale. The %CV values printed in
    # square brackets are derived from this variance and are not
    # independently estimated, so the internal variance encoded here
    # equals the tabled "Estimate" column directly.
    etalcl ~ 0.735
    # Kay 2022 Table 2 Omega(1,1): IIV CL/F = 0.735 (95% CI 0.526, 0.945; CV% 104, shrinkage 5.86)

    etalvc ~ 0.813
    # Kay 2022 Table 2 Omega(2,2): IIV V2/F = 0.813 (95% CI 0.504, 1.12; CV% 112, shrinkage 32.6)

    # IIV-Q/F was FIXED during estimation per Table 2 95%-CI column
    # ("FIXED" for Omega(3,3)).
    etalq ~ fixed(0.0250)
    # Kay 2022 Table 2 Omega(3,3): IIV Q/F = 0.0250 FIXED (CV% 15.9, shrinkage 67.5)

    etalvp ~ 0.956
    # Kay 2022 Table 2 Omega(4,4): IIV V3/F = 0.956 (95% CI 0.590, 1.32; CV% 127, shrinkage 36.4)

    etalka ~ 0.0280
    # Kay 2022 Table 2 Omega(5,5): IIV ka = 0.0280 (95% CI 0.00705, 0.0490; CV% 16.9, shrinkage 58.2)

    # Residual error. Kay 2022 Table 2 Sigma(1,1) row "Proportional ...
    # Variance" = 0.200 (95% CI 0.173, 0.226). Footnote: "CV% of sigma =
    # sqrt(estimate) * 100" -- so the tabled value is the variance of a
    # proportional residual on the linear-concentration scale, and the
    # corresponding fractional SD is sqrt(0.200) = 0.4472, matching the
    # tabled CV% = 44.7.
    propSd <- sqrt(0.200)
    label("Proportional residual SD for lumefantrine plasma concentration (fraction)")
    # Kay 2022 Table 2 Sigma(1,1): proportional variance = 0.200 (95% CI 0.173, 0.226; CV% 44.7)
  })

  model({
    # Reference values for the Kay 2022 parameterisation (Table 2 footnote):
    # 50-month-old, 15 kg HIV-uninfected child.
    ref_wt   <- 15        # kg
    ref_agemo <- 50       # months

    # Age in months drives both the piecewise allometric exponent on
    # apparent oral clearance terms and the power-form age effect on
    # relative bioavailability. The canonical covariate column AGE is in
    # years, so convert to months here.
    agemo <- AGE * 12

    # Piecewise age-dependent allometric exponent on CL/F and Q/F
    # (Anderson & Holford 2009; paper ref 34). Volumes use exponent 1.
    # Results, Population pharmacokinetics: "Fixed effects on volume used
    # an allometric exponent of 1, whereas the fixed effects on clearance
    # used an exponent of 0.75, 0.9, 1.0, or 1.2 for children age > 60
    # months, > 24 to 60 months, > 3 to 24 months, and <= 3 months,
    # respectively."
    allo_cl <- 1.2  * (agemo <= 3)                +
               1.0  * (agemo >  3 & agemo <= 24)  +
               0.9  * (agemo > 24 & agemo <= 60)  +
               0.75 * (agemo > 60)

    # Individual PK parameters with allometric scaling and ART
    # linear-deviation effects. CL/F, V2/F, Q/F, V3/F all scale with WT
    # relative to the 15-kg reference; CL/F and Q/F use the age-dependent
    # exponent, V2/F and V3/F use exponent 1. CL/F and ka additionally
    # carry the three mutually-exclusive ART indicator effects. IIV is
    # log-normal multiplicative.
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^allo_cl *
            (1 + e_efv_cl * CONMED_EFV) *
            (1 + e_lpv_cl * CONMED_LPV) *
            (1 + e_nvp_cl * CONMED_NVP)
    vc <- exp(lvc + etalvc) * (WT / ref_wt)
    q  <- exp(lq  + etalq ) * (WT / ref_wt)^allo_cl
    vp <- exp(lvp + etalvp) * (WT / ref_wt)
    ka <- exp(lka + etalka) *
            (1 + e_efv_ka * CONMED_EFV) *
            (1 + e_lpv_ka * CONMED_LPV) *
            (1 + e_nvp_ka * CONMED_NVP)

    # Power-form age effect on relative bioavailability.
    # F = exp(lfdepot) * (agemo / ref_agemo)^e_age_f with lfdepot = log(1)
    # the structural anchor.
    fdepot <- exp(lfdepot) * (agemo / ref_agemo)^e_age_f

    # Disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order absorption from depot into a two-compartment
    # disposition.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Plasma concentration in ng/mL. Dose in mg, Vc in L gives
    # central / vc in mg/L = ug/mL; multiplying by 1000 yields ng/mL,
    # the unit Kay 2022 reports throughout (Figures, LLOQ 50 ng/mL).
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
