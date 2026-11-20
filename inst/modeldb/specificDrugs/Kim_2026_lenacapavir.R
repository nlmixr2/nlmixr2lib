Kim_2026_lenacapavir <- function() {
  description <- paste(
    "QSP (semi-mechanistic PK-PD with HIV-1 viral dynamics). Long-acting",
    "lenacapavir (LEN) for HIV-1 treatment and pre-exposure prophylaxis",
    "after oral tablets, subcutaneous (SC) or intramuscular (IM) injection,",
    "from Kim et al. 2026. Plasma LEN is a one-compartment model with",
    "first-order elimination fed by first-order oral absorption (oral",
    "apparent volume a power function of the oral dose) or, for an",
    "injection, by two parallel release processes: a direct first-order",
    "release of a fraction of the dose and an indirect release of the rest",
    "through an Erlang-distributed chain of one to three transit",
    "compartments. Release parameters are specific to each injectable",
    "formulation (SC aqueous suspension, SC PEG/water solution, IM",
    "ethanol/water solution with 5% or 10% w/w ethanol). LEN blocks the",
    "production of infectious virions (sigmoid Emax, Hill 2.1) inside an",
    "eight-state HIV-1 viral dynamics model (T cells and macrophages,",
    "uninfected / early infected / productively infected, plus infectious",
    "and non-infectious virus) that starts at its untreated steady state;",
    "outputs are plasma viral load and its log10 change from baseline.",
    "Also reports, for wild-type virus and seven capsid-resistant mutants,",
    "the reproduction number and the prophylactic efficacy that the current",
    "plasma LEN concentration would give if held constant. Fitted to",
    "published mean concentration-time profiles by least squares:",
    "typical-value model with no IIV and no residual error."
  )
  reference <- paste(
    "Kim HY, Liebenberg A, Zhang L, von Kleist M. Pharmacokinetics,",
    "pharmacodynamics, efficacy and drug resistance selection of injectable",
    "long-acting lenacapavir pre-exposure prophylaxis (PrEP) against HIV.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70254.",
    "doi:10.1002/psp4.70254. PMCID PMC13093070.",
    "PK structure from Methods Equations (1)-(3) and Figure 1A; PK",
    "parameters from Table 1 and its caption; Emax model from Equation (5);",
    "viral dynamics from Supplementary Equations (S1)-(S8) and Table S3;",
    "reproduction number and infection probability from Equations (6)-(10)",
    "and Supplementary Equations (S13)-(S18); mutant phenotypes from Table",
    "S4. The LEN effect on virion production in the full viral dynamics",
    "model and the plasma and interstitial volumes that convert total-body",
    "virus to plasma viral load are from the authors' code deposit",
    "https://github.com/KleistLab/LenPrEP (GPL-3.0, commit 3933f74)."
  )
  vignette <- "Kim_2026_lenacapavir"

  # Doses and the five dosing-site states are in mg. central is an amount in
  # mg on the parenteral volume scale (see the f(depot_oral) note in model());
  # Cc = central / vc is converted from mg/L to ng/mL. The viral-dynamics
  # states are cell and virion counts in the whole body.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Oral tablets dose depot_oral. Every injection is TWO dose records of the
  # full injected amount, one to depot (direct release) and one to transit1
  # (indirect release); f() splits the amount between them.
  dosing <- c("depot_oral", "depot", "transit1")

  paper_specific_compartments <- c(
    "tcell_uninfected",
    "tcell_early",
    "tcell_infected",
    "macro_uninfected",
    "macro_early",
    "macro_infected",
    "virus_noninfectious"
  )

  covariateData <- list(
    DOSE_LENACAPAVIR_PO_MG = list(
      description = paste(
        "Nominal lenacapavir oral tablet dose per administration, in mg."
      ),
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the dose-dependent apparent volume of oral LEN, Vd(dose) =",
        "508.16 * dose^0.61 L (Table 1 caption and Methods 'Parameter",
        "Estimation'), which the model applies as the relative input",
        "f(depot_oral) = Vd / Vd(dose). Only oral dose records use it; a",
        "regimen without oral doses may carry 0 (then f(depot_oral) = 1 and",
        "has no effect). Constant within a subject."
      ),
      source_name = "dose (oral, mg)"
    ),
    FORM_LENACAPAVIR_SCAQ = list(
      description = paste(
        "1 = injections are the subcutaneous aqueous suspension (Studies 1",
        "and 3), 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SC PEG/water solution when all formulation indicators are 0)",
      notes = paste(
        "Selects the Table 1 release parameters. The aqueous suspension was",
        "fitted twice, to Study 1 (healthy adults) and Study 3 (adults with",
        "HIV-1); STUDY_GSUS2004072 chooses between the two sets. At most one",
        "of FORM_LENACAPAVIR_SCAQ, FORM_LENACAPAVIR_IM_ETOH5 and",
        "FORM_LENACAPAVIR_IM_ETOH10 may be 1."
      ),
      source_name = "SC aqueous suspension (Table S1)"
    ),
    STUDY_GSUS2004072 = list(
      description = paste(
        "1 = the Study 3 (GS-US-200-4072, adults with untreated HIV-1)",
        "parameter set of the SC aqueous suspension, 0 = the Study 1",
        "(GS-US-200-4070, healthy adults) set."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Study 1 set)",
      notes = paste(
        "Only read when FORM_LENACAPAVIR_SCAQ = 1. Study 3 is the",
        "monotherapy study whose viral-load decline fixed the IC50, so the",
        "Figure 3 replication uses this set."
      ),
      source_name = "Study 3 (Table 1, Table S1)"
    ),
    FORM_LENACAPAVIR_IM_ETOH5 = list(
      description = "1 = injections are the IM ethanol/water solution with 5% w/w ethanol (Study 5, F1), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the Study 5 release rates and the F1 direct fraction (12.8%).",
      source_name = "Study 5 F1 (Table 1, Table S1)"
    ),
    FORM_LENACAPAVIR_IM_ETOH10 = list(
      description = "1 = injections are the IM ethanol/water solution with 10% w/w ethanol (Study 5, F2), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Selects the Study 5 release rates and the F2 direct fraction (43.6%).",
      source_name = "Study 5 F2 (Table 1, Table S1)"
    )
  )

  compartmentData <- list(
    depot_oral = list(analyte = "lenacapavir", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "lenacapavir", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "lenacapavir", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "lenacapavir", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "lenacapavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lenacapavir", units = "mg", specimen = "plasma", verified = TRUE),
    tcell_uninfected = list(
      analyte = "uninfected CD4+ T cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early = list(
      analyte = "early infected CD4+ T cells (pre-integration)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected = list(
      analyte = "productively infected CD4+ T cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_uninfected = list(
      analyte = "uninfected macrophages",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_early = list(
      analyte = "early infected macrophages (pre-integration)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_infected = list(
      analyte = "productively infected macrophages",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(analyte = "infectious HIV-1 virions", units = "virions", specimen = "plasma", verified = TRUE),
    virus_noninfectious = list(
      analyte = "non-infectious HIV-1 virions",
      units = "virions",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 136,
    n_studies = 5,
    disease_state = paste(
      "Healthy adults (Studies 1, 2, 4 and 5) and adults with untreated",
      "HIV-1 infection (Study 3, also the viral-load data)"
    ),
    dose_range = paste(
      "Single doses: oral tablets 300, 900 and 1800 mg; SC aqueous",
      "suspension 20-450 mg; SC PEG/water solution 309 and 927 mg; IM",
      "ethanol/water solution 5000 mg"
    ),
    regions = "Not reported (Gilead phase 1 and 1b studies)",
    notes = paste(
      "Model fitted to digitised published MEAN or MEDIAN profiles, not",
      "individual data (Supplementary Table S1): 15 single-dose datasets",
      "(3 oral, 10 SC, 2 IM) from GS-US-200-4070 (Study 1, n = 8 per",
      "cohort), GS-US-200-4071 (Study 2, n = 8 per cohort),",
      "GS-US-200-4072 (Study 3, n = 6 per cohort), GS-US-200-4358 (Study",
      "4, n = 8 per cohort) and the once-yearly IM study of Jogiraju 2025",
      "(Study 5, n = 20 per cohort). n_subjects sums the cohort sizes;",
      "demographics were not reported. IC50 was fitted to mean viral-load",
      "declines after single SC doses of 20, 50, 150 and 450 mg in Study 3.",
      "Inter-individual variability is represented only by fixed",
      "0.15x / 4.0x multiples of the median concentration taken from the",
      "PURPOSE 2 trial (Methods 'PK Variability')."
    )
  )

  ini({
    # ---- Plasma LEN PK (Table 1; rate constants in 1/h) -------------------
    lkel <- fixed(log(0.0026)); label("Elimination rate constant ke (1/h)") # Table 1, ke = 0.0026 1/h; graphically determined from oral data and fixed for the parenteral fits (Methods 'LEN Parenteral PK')
    lvc <- log(1618); label("Volume of distribution for parenteral dosing Vd (L)") # Table 1, Vd = 1618 L (one Vd estimated across all parenteral formulations)
    lka_oral <- fixed(log(1.08)); label("Oral first-order absorption rate constant ka (1/h)") # Table 1, ka = 1.08 1/h (Study 2), computed from Tmax with Equation (4), not optimised
    lvc_oral <- log(508.16); label("Oral apparent volume at a 1 mg oral dose, power-law coefficient (L)") # Table 1 caption, f(dose) = 508.16 * dose^0.61
    e_dose_vc_oral <- 0.61; label("Power-law exponent of the oral apparent volume on the oral dose (unitless)") # Table 1 caption, f(dose) = 508.16 * dose^0.61

    # ---- Injection release, by formulation (Table 1) ---------------------
    # SC PEG/water solution, Study 4 (reference formulation)
    lka_scpeg <- log(0.039); label("Direct-release rate constant kdirect, SC PEG/water (1/h)") # Table 1, kdirect Study 4 = 0.039
    lktr_scpeg <- log(6.3e-4); label("Indirect-release transit rate constant kindirect, SC PEG/water (1/h)") # Table 1, kindirect Study 4 = 6.3e-4
    logitfrel_scpeg <- log(0.015 / (1 - 0.015)); label("Direct-release fraction Frac, SC PEG/water (logit of fraction)") # Table 1, Frac Study 4 = 1.5 %
    ntr_scpeg <- fixed(2); label("Number of indirect-release transit compartments n, SC PEG/water (count)") # Table 1, n Study 4 = 2
    # SC aqueous suspension, Study 1 (healthy adults)
    lka_scaq1 <- log(0.1); label("Direct-release rate constant kdirect, SC aqueous Study 1 (1/h)") # Table 1, kdirect Study 1 = 0.1
    lktr_scaq1 <- log(6.5e-4); label("Indirect-release transit rate constant kindirect, SC aqueous Study 1 (1/h)") # Table 1, kindirect Study 1 = 6.5e-4
    logitfrel_scaq1 <- log(0.12 / (1 - 0.12)); label("Direct-release fraction Frac, SC aqueous Study 1 (logit of fraction)") # Table 1, Frac Study 1 = 12 %
    ntr_scaq1 <- fixed(1); label("Number of indirect-release transit compartments n, SC aqueous Study 1 (count)") # Table 1, n Study 1 = 1
    # SC aqueous suspension, Study 3 (adults with HIV-1)
    lka_scaq3 <- log(0.063); label("Direct-release rate constant kdirect, SC aqueous Study 3 (1/h)") # Table 1, kdirect Study 3 = 0.063
    lktr_scaq3 <- log(0.0052); label("Indirect-release transit rate constant kindirect, SC aqueous Study 3 (1/h)") # Table 1, kindirect Study 3 = 0.0052
    logitfrel_scaq3 <- log(0.16 / (1 - 0.16)); label("Direct-release fraction Frac, SC aqueous Study 3 (logit of fraction)") # Table 1, Frac Study 3 = 16 %
    ntr_scaq3 <- fixed(3); label("Number of indirect-release transit compartments n, SC aqueous Study 3 (count)") # Table 1, n Study 3 = 3
    # IM ethanol/water solution, Study 5 (rates and n shared by F1 and F2)
    lka_im <- log(2.6e-4); label("Direct-release rate constant kdirect, IM ethanol/water (1/h)") # Table 1, kdirect Study 5 = 2.6e-4
    lktr_im <- log(4.1e-4); label("Indirect-release transit rate constant kindirect, IM ethanol/water (1/h)") # Table 1, kindirect Study 5 = 4.1e-4
    logitfrel_im_etoh5 <- log(0.128 / (1 - 0.128)); label("Direct-release fraction Frac, IM 5% w/w ethanol F1 (logit of fraction)") # Table 1, Frac Study 5 = 12.8 (F1)
    logitfrel_im_etoh10 <- log(0.436 / (1 - 0.436)); label("Direct-release fraction Frac, IM 10% w/w ethanol F2 (logit of fraction)") # Table 1, Frac Study 5 = 43.6 (F2)
    ntr_im <- fixed(2); label("Number of indirect-release transit compartments n, IM ethanol/water (count)") # Table 1, n Study 5 = 2

    # ---- PK variability envelope (Methods 'PK Variability') ---------------
    pkvar_lower <- fixed(0.15); label("Lower PK-variability limit as a multiple of the median concentration (unitless)") # Methods 'PK Variability', 0.15 times the median
    pkvar_upper <- fixed(4.0); label("Upper PK-variability limit as a multiple of the median concentration (unitless)") # Methods 'PK Variability', 4.0 times the median

    # ---- Pharmacodynamics (Equation 5) -------------------------------------
    lic50 <- log(1.6); label("Plasma LEN concentration for 50% inhibition of infectious-virion production IC50 (ng/mL)") # Results, IC50 = 1.6 ng/mL; Table S3 caption
    hill <- fixed(2.1); label("Hill coefficient m of the Emax model (unitless)") # Methods 'Estimation of Antiviral Potency', m = 2.1 (set from Smith 2024)

    # ---- HIV-1 viral dynamics (Table S3; rates in 1/day) ------------------
    lambda_t <- fixed(2e9); label("Production of uninfected T cells lambdaT (cells/day)") # Table S3, lambdaT = 2 x 10^9
    delta_t <- fixed(0.02); label("Death rate of uninfected T cells deltaT (1/day)") # Table S3, deltaT = 0.02
    delta_t1 <- fixed(0.7); label("Death rate of early infected T cells deltaT1 (1/day)") # Table S3, deltaT1 = 0.7
    delta_t2 <- fixed(1); label("Death rate of productively infected T cells deltaT2 (1/day)") # Table S3, deltaT2 = 1
    delta_pic_t <- fixed(0.35); label("Intracellular destruction of the pre-integration complex in T cells deltaPIC,T (1/day)") # Table S3, deltaPICT = 0.35
    k_t <- fixed(0.35); label("Integration (early to productive) rate in T cells kT (1/day)") # Table S3, kT = 0.35
    beta_t <- fixed(1.8e-12); label("Infection rate of T cells betaT (1/(virion day))") # Table S3, betaT = 1.8 x 10^-12 (adjusted to a 4.5 log10 baseline viral load)
    n_t_total <- fixed(1000); label("Total virions released per productively infected T cell NhatT (virions/(cell day))") # Table S3, NhatT = 1000
    n_t <- fixed(670); label("Infectious virions released per productively infected T cell NT (virions/(cell day))") # Table S3, NT = 670
    lambda_m <- fixed(6.9e7); label("Production of uninfected macrophages lambdaM (cells/day)") # Table S3, lambdaM = 6.9 x 10^7
    delta_m <- fixed(0.0069); label("Death rate of uninfected macrophages deltaM (1/day)") # Table S3, deltaM = 0.0069
    delta_m1 <- fixed(0.0069); label("Death rate of early infected macrophages deltaM1 (1/day)") # Table S3, deltaM1 = 0.0069
    delta_m2 <- fixed(0.09); label("Death rate of productively infected macrophages deltaM2 (1/day)") # Table S3, deltaM2 = 0.09
    delta_pic_m <- fixed(0.0035); label("Intracellular destruction of the pre-integration complex in macrophages deltaPIC,M (1/day)") # Table S3, deltaPICM = 0.0035
    k_m <- fixed(0.07); label("Integration (early to productive) rate in macrophages kM (1/day)") # Table S3, kM = 0.07
    beta_m <- fixed(2e-14); label("Infection rate of macrophages betaM (1/(virion day))") # Table S3, betaM = 2 x 10^-14 (adjusted to a 4.5 log10 baseline viral load)
    n_m_total <- fixed(100); label("Total virions released per productively infected macrophage NhatM (virions/(cell day))") # Table S3, NhatM = 100
    n_m <- fixed(67); label("Infectious virions released per productively infected macrophage NM (virions/(cell day))") # Table S3, NM = 67
    cl_virus <- fixed(23); label("Clearance rate of free virus in infected individuals CL (1/day)") # Table S3 caption, CL = 23 (simulation of treatment effects)
    cl_virus_naive <- fixed(2.3); label("Clearance rate of free virus in virus-naive individuals, prophylaxis setting (1/day)") # Table S3 caption, CL = 2.3 (PrEP simulations)
    p_infect_success <- fixed(0.5); label("Probability that a virion-target-cell encounter infects; sets CLT = (1/p - 1) * beta (unitless)") # Table S3 caption, CLT/M = (1/0.5 - 1) * betaT/M

    # ---- Plasma viral load conversion (Figure 3 caption) --------------------
    k_int_plas <- fixed(50); label("Interstitial-to-plasma virus partition coefficient Kint:plas (unitless)") # Figure 3 caption, Kint:plas approx 50
    v_int <- fixed(9.6); label("Interstitial fluid volume Vint (L)") # not stated in the paper; authors' code deposit utils.py solve_pd_ode() uses 9600 mL (von Kleist 2010)
    v_plas <- fixed(3.1); label("Plasma volume Vplas (L)") # not stated in the paper; authors' code deposit utils.py solve_pd_ode() uses 3100 mL (von Kleist 2010)
    rna_per_virion <- fixed(2); label("HIV-1 RNA copies per virion (copies)") # Figure 3 caption, plasma VL contains on average two viral RNAs per virion

    # ---- Capsid-resistance mutant phenotypes (Table S4) ---------------------
    # FC = fold change in EC50 vs wild type (median); fit = replication
    # capacity relative to wild type (average), used as s(i) in Eqs. (7) and (S18).
    fc_q67h <- fixed(4.7); label("EC50 fold change, Q67H (unitless)") # Table S4, Q67H FC median 4.7
    fit_q67h <- fixed(0.825); label("Relative replication capacity s, Q67H (fraction)") # Table S4, Q67H replication capacity 82.5 %
    fc_n74d <- fixed(18.1); label("EC50 fold change, N74D (unitless)") # Table S4, N74D FC median 18.1
    fit_n74d <- fixed(0.48); label("Relative replication capacity s, N74D (fraction)") # Table S4, N74D replication capacity 48 %
    fc_q67h_n74d <- fixed(467); label("EC50 fold change, Q67H+N74D (unitless)") # Table S4, Q67H+N74D FC median 467
    fit_q67h_n74d <- fixed(0.29); label("Relative replication capacity s, Q67H+N74D (fraction)") # Table S4, Q67H+N74D replication capacity 29 %
    fc_q67h_t107n <- fixed(47.8); label("EC50 fold change, Q67H+T107N (unitless)") # Table S4, Q67H+T107N FC median 47.8
    fit_q67h_t107n <- fixed(0.41); label("Relative replication capacity s, Q67H+T107N (fraction)") # Table S4, Q67H+T107N replication capacity 41 %
    fc_m66i_t107a <- fixed(386.3); label("EC50 fold change, M66I+T107A (unitless)") # Table S4, M66I+T107A FC median 386.3
    fit_m66i_t107a <- fixed(0.154); label("Relative replication capacity s, M66I+T107A (fraction)") # Table S4, M66I+T107A replication capacity 15.4 %
    fc_q67h_k70r <- fixed(28.6); label("EC50 fold change, Q67H+K70R (unitless)") # Table S4, Q67H+K70R FC median 28.6
    fit_q67h_k70r <- fixed(0.6115); label("Relative replication capacity s, Q67H+K70R (fraction)") # Table S4, Q67H+K70R replication capacity 61.15 %
    fc_q67h_k70r_t107n <- fixed(259); label("EC50 fold change, Q67H+K70R+T107N (unitless)") # Table S4, Q67H+K70R+T107N FC median 259
    fit_q67h_k70r_t107n <- fixed(0.8305); label("Relative replication capacity s, Q67H+K70R+T107N (fraction)") # Table S4, Q67H+K70R+T107N replication capacity 83.05 %
  })

  model({
    # ---- 1. PK parameters (no IIV) -------------------------------------------
    kel <- exp(lkel)
    vc <- exp(lvc)
    ka_oral <- exp(lka_oral)
    vc_oral <- vc
    if (DOSE_LENACAPAVIR_PO_MG > 0) {
      vc_oral <- exp(lvc_oral) * DOSE_LENACAPAVIR_PO_MG^e_dose_vc_oral
    }

    # ---- 2. Formulation-specific release parameters (Table 1) -------------
    is_im <- FORM_LENACAPAVIR_IM_ETOH5 + FORM_LENACAPAVIR_IM_ETOH10
    is_scaq1 <- FORM_LENACAPAVIR_SCAQ * (1 - STUDY_GSUS2004072)
    is_scaq3 <- FORM_LENACAPAVIR_SCAQ * STUDY_GSUS2004072
    is_scpeg <- 1 - FORM_LENACAPAVIR_SCAQ - is_im

    ka_inj <- is_scpeg * exp(lka_scpeg) + is_scaq1 * exp(lka_scaq1) +
      is_scaq3 * exp(lka_scaq3) + is_im * exp(lka_im)
    ktr <- is_scpeg * exp(lktr_scpeg) + is_scaq1 * exp(lktr_scaq1) +
      is_scaq3 * exp(lktr_scaq3) + is_im * exp(lktr_im)
    ntr <- is_scpeg * ntr_scpeg + is_scaq1 * ntr_scaq1 +
      is_scaq3 * ntr_scaq3 + is_im * ntr_im
    logitfrel_ind <- is_scpeg * logitfrel_scpeg + is_scaq1 * logitfrel_scaq1 +
      is_scaq3 * logitfrel_scaq3 +
      FORM_LENACAPAVIR_IM_ETOH5 * logitfrel_im_etoh5 +
      FORM_LENACAPAVIR_IM_ETOH10 * logitfrel_im_etoh10
    frel <- expit(logitfrel_ind)

    # Equation (3) feeds central with the Erlang(n, kindirect) density times
    # the indirectly released amount, which is exactly the output of a chain
    # of n transit compartments of rate kindirect dosed in its first member.
    # The chain is three long and n selects which member empties into
    # central, so one dose compartment (transit1) serves every n.
    to2 <- 0
    to3 <- 0
    if (ntr > 1.5) to2 <- 1
    if (ntr > 2.5) to3 <- 1
    out1 <- 1 - to2
    out2 <- to2 - to3
    out3 <- to3

    # ---- 3. PK ODEs (Equations 1-3) -----------------------------------------
    d/dt(depot_oral) <- -ka_oral * depot_oral
    d/dt(depot) <- -ka_inj * depot
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- to2 * ktr * transit1 - ktr * transit2
    d/dt(transit3) <- to3 * ktr * transit2 - ktr * transit3
    d/dt(central) <- ka_oral * depot_oral + ka_inj * depot +
      ktr * (out1 * transit1 + out2 * transit2 + out3 * transit3) -
      kel * central

    # Equation (3) is written for the concentration with a route-specific
    # volume: oral input divided by Vd(dose), injected input by Vd. Keeping
    # central on the Vd scale, the oral input enters as Vd / Vd(dose).
    f(depot_oral) <- vc / vc_oral
    # PIdir(0) = Frac * dose and PIind(0) = (1 - Frac) * dose.
    f(depot) <- frel
    f(transit1) <- 1 - frel

    # mg/L -> ng/mL
    Cc <- central / vc * 1000
    Cc_lower <- pkvar_lower * Cc
    Cc_upper <- pkvar_upper * Cc

    # ---- 4. Drug effect (Equation 5) ------------------------------------------
    ic50 <- exp(lic50)
    eta_wt <- Cc^hill / (Cc^hill + ic50^hill)

    # ---- 5. HIV-1 viral dynamics (Equations S1-S8; rates in 1/day, / 24) ----
    # LEN converts infectious into non-infectious virion production: NT and
    # NM are scaled by (1 - eta) and the total output NhatT / NhatM is
    # unchanged (deposit utils.py ode_complex()).
    cl_t <- (1 / p_infect_success - 1) * beta_t
    cl_m <- (1 / p_infect_success - 1) * beta_m
    n_t_drug <- (1 - eta_wt) * n_t
    n_m_drug <- (1 - eta_wt) * n_m

    # Untreated steady state (the baseline the authors simulate to before
    # dosing), in closed form: with no drug, dT1 = dT2 = dM1 = dM2 = dV = 0
    # reduce to a quadratic in V.
    a_t <- delta_t1 + k_t + delta_pic_t
    a_m <- delta_m1 + k_m + delta_pic_m
    g_t <- beta_t * (delta_t1 + k_t) / a_t
    g_m <- beta_m * (delta_m1 + k_m) / a_m
    r_t <- n_t * k_t * beta_t / (a_t * delta_t2) - (beta_t + cl_t)
    r_m <- n_m * k_m * beta_m / (a_m * delta_m2) - (beta_m + cl_m)
    qa <- cl_virus * g_t * g_m
    qb <- cl_virus * (delta_t * g_m + delta_m * g_t) - r_t * lambda_t * g_m -
      r_m * lambda_m * g_t
    qc <- cl_virus * delta_t * delta_m - r_t * lambda_t * delta_m -
      r_m * lambda_m * delta_t
    v_ss <- (-qb + sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)
    tu_ss <- lambda_t / (delta_t + g_t * v_ss)
    mu_ss <- lambda_m / (delta_m + g_m * v_ss)
    t1_ss <- beta_t * v_ss * tu_ss / a_t
    m1_ss <- beta_m * v_ss * mu_ss / a_m
    t2_ss <- k_t * t1_ss / delta_t2
    m2_ss <- k_m * m1_ss / delta_m2
    vn_ss <- ((n_t_total - n_t) * t2_ss + (n_m_total - n_m) * m2_ss) / cl_virus

    tcell_uninfected(0) <- tu_ss
    tcell_early(0) <- t1_ss
    tcell_infected(0) <- t2_ss
    macro_uninfected(0) <- mu_ss
    macro_early(0) <- m1_ss
    macro_infected(0) <- m2_ss
    virus(0) <- v_ss
    virus_noninfectious(0) <- vn_ss

    d/dt(tcell_uninfected) <- (lambda_t - delta_t * tcell_uninfected -
      beta_t * virus * tcell_uninfected + delta_pic_t * tcell_early) / 24
    d/dt(macro_uninfected) <- (lambda_m - delta_m * macro_uninfected -
      beta_m * virus * macro_uninfected + delta_pic_m * macro_early) / 24
    d/dt(tcell_early) <- (beta_t * virus * tcell_uninfected -
      a_t * tcell_early) / 24
    d/dt(macro_early) <- (beta_m * virus * macro_uninfected -
      a_m * macro_early) / 24
    d/dt(tcell_infected) <- (k_t * tcell_early - delta_t2 * tcell_infected) / 24
    d/dt(macro_infected) <- (k_m * macro_early - delta_m2 * macro_infected) / 24
    d/dt(virus) <- (n_t_drug * tcell_infected + n_m_drug * macro_infected -
      virus * (cl_virus + (beta_t + cl_t) * tcell_uninfected +
        (beta_m + cl_m) * macro_uninfected)) / 24
    d/dt(virus_noninfectious) <- ((n_t_total - n_t_drug) * tcell_infected +
      (n_m_total - n_m_drug) * macro_infected -
      cl_virus * virus_noninfectious) / 24

    # Plasma HIV-1 RNA (copies/mL): total-body virus distributed over plasma
    # and interstitial space (Figure 3 caption), volumes in mL.
    v_dist_ml <- (k_int_plas * v_int + v_plas) * 1000
    viral_load <- rna_per_virion * (virus + virus_noninfectious) / v_dist_ml
    viral_load0 <- rna_per_virion * (v_ss + vn_ss) / v_dist_ml
    d_viral_load <- log10(viral_load / viral_load0)

    # ---- 6. Reproduction number and prophylactic efficacy (Eqs. 6-10) -----
    # Virus-naive setting: uninfected T cells at lambdaT / deltaT, free-virus
    # clearance 2.3 /day. Per-entity propensities of Eqs. (S13)-(S18).
    tu_naive <- lambda_t / delta_t
    pa1 <- cl_virus_naive + cl_t * tu_naive
    pa4 <- beta_t * tu_naive
    pa2 <- delta_pic_t + delta_t1
    pa5 <- k_t
    # Equation (6)
    lambda_inf <- pa4 / (pa1 + pa4) * pa5 / (pa2 + pa5)
    # Equation (7) for the drug-free wild type
    r0_wt <- lambda_inf * n_t / delta_t2
    pinf_none <- lambda_inf * (1 - 1 / r0_wt)

    # Equation (8): R_t(C, i) = R0 * s(i) * (1 - eta(C, FC(i)))
    rt_wt <- r0_wt * (1 - eta_wt)
    rt_q67h <- r0_wt * fit_q67h *
      (1 - Cc^hill / (Cc^hill + (fc_q67h * ic50)^hill))
    rt_n74d <- r0_wt * fit_n74d *
      (1 - Cc^hill / (Cc^hill + (fc_n74d * ic50)^hill))
    rt_q67h_n74d <- r0_wt * fit_q67h_n74d *
      (1 - Cc^hill / (Cc^hill + (fc_q67h_n74d * ic50)^hill))
    rt_q67h_t107n <- r0_wt * fit_q67h_t107n *
      (1 - Cc^hill / (Cc^hill + (fc_q67h_t107n * ic50)^hill))
    rt_m66i_t107a <- r0_wt * fit_m66i_t107a *
      (1 - Cc^hill / (Cc^hill + (fc_m66i_t107a * ic50)^hill))
    rt_q67h_k70r <- r0_wt * fit_q67h_k70r *
      (1 - Cc^hill / (Cc^hill + (fc_q67h_k70r * ic50)^hill))
    rt_q67h_k70r_t107n <- r0_wt * fit_q67h_k70r_t107n *
      (1 - Cc^hill / (Cc^hill + (fc_q67h_k70r_t107n * ic50)^hill))

    # Equation (9), P_inf = max(0, Lambda * (1 - 1 / R_t)), per transmitted
    # virion; Equation (10), phi = 1 - P_inf / P_inf(wild type, no drug).
    pinf_wt <- 0
    if (rt_wt > 1) pinf_wt <- lambda_inf * (1 - 1 / rt_wt)
    pinf_q67h <- 0
    if (rt_q67h > 1) pinf_q67h <- lambda_inf * (1 - 1 / rt_q67h)
    pinf_n74d <- 0
    if (rt_n74d > 1) pinf_n74d <- lambda_inf * (1 - 1 / rt_n74d)
    pinf_q67h_n74d <- 0
    if (rt_q67h_n74d > 1) pinf_q67h_n74d <- lambda_inf * (1 - 1 / rt_q67h_n74d)
    pinf_q67h_t107n <- 0
    if (rt_q67h_t107n > 1) pinf_q67h_t107n <- lambda_inf * (1 - 1 / rt_q67h_t107n)
    pinf_m66i_t107a <- 0
    if (rt_m66i_t107a > 1) pinf_m66i_t107a <- lambda_inf * (1 - 1 / rt_m66i_t107a)
    pinf_q67h_k70r <- 0
    if (rt_q67h_k70r > 1) pinf_q67h_k70r <- lambda_inf * (1 - 1 / rt_q67h_k70r)
    pinf_q67h_k70r_t107n <- 0
    if (rt_q67h_k70r_t107n > 1) {
      pinf_q67h_k70r_t107n <- lambda_inf * (1 - 1 / rt_q67h_k70r_t107n)
    }

    phi_wt <- 1 - pinf_wt / pinf_none
    phi_q67h <- 1 - pinf_q67h / pinf_none
    phi_n74d <- 1 - pinf_n74d / pinf_none
    phi_q67h_n74d <- 1 - pinf_q67h_n74d / pinf_none
    phi_q67h_t107n <- 1 - pinf_q67h_t107n / pinf_none
    phi_m66i_t107a <- 1 - pinf_m66i_t107a / pinf_none
    phi_q67h_k70r <- 1 - pinf_q67h_k70r / pinf_none
    phi_q67h_k70r_t107n <- 1 - pinf_q67h_k70r_t107n / pinf_none
  })
}
