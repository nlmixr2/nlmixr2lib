Ding_2024_ivermectin <- function() {
  description <- paste(
    "In-host dengue virus and non-structural protein 1 (NS1) kinetic model",
    "with ivermectin population PK-PD in Thai adults and children with acute",
    "dengue fever (Ding 2024). Ivermectin PK: two-compartment disposition",
    "with a dose compartment plus six transit absorption compartments",
    "(ktr = 7 / MTT), relative bioavailability fixed to 1 with IIV, and",
    "allometric body-weight scaling (estimated exponents, centred on 70 kg).",
    "Dengue viral load: target-cell-limited model (target cells, an",
    "infected-cell eclipse chain of four compartments, virus-producing",
    "cells, free virus) with innate clearance c and an antibody-driven",
    "neutralisation term kA * A(t) * V, where A(t) is the mean of the",
    "normalised sigmoidal-Emax time courses of dengue-specific IgM and IgG.",
    "NS1: secreted from producing cells through nine transit compartments",
    "with a time-dependent sigmoidal-Emax secretion rate and a first-order",
    "elimination rate that rises linearly with time when NS1-specific IgG",
    "exceeds a cut-off. Ivermectin inhibits NS1 secretion through an Emax",
    "function of its plasma concentration (Emax fixed to 0.99); it has no",
    "effect on viral load. Model time is days since infection (fever onset",
    "assumed on day 6).",
    sep = " "
  )
  reference <- paste(
    "Ding J, Mairiang D, Prayongkul D, Puttikhunt C, Noisakran S, Kaewjiw N,",
    "Songjaeng A, Prommool T, Tangthawornchaikul N, Angkasekwinai N,",
    "Suputtamongkol Y, Lapphra K, Chokephaibulkit K, White NJ, Avirutnan P,",
    "Tarning J. In-host modeling of dengue virus and non-structural protein 1",
    "and the effects of ivermectin in patients with acute dengue fever.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2196-2208.",
    "doi:10.1002/psp4.13233",
    sep = " "
  )
  vignette <- "Ding_2024_ivermectin"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Viral-dynamics states with no register entry: the four-compartment
  # infected-cell eclipse chain (non-producing infected cells), the nine-
  # compartment NS1 secretion delay chain and the plasma NS1 antigen pool.
  paper_specific_compartments <- c(
    "eclipse1",
    "eclipse2",
    "eclipse3",
    "eclipse4",
    "ns1_transit1",
    "ns1_transit2",
    "ns1_transit3",
    "ns1_transit4",
    "ns1_transit5",
    "ns1_transit6",
    "ns1_transit7",
    "ns1_transit8",
    "ns1_transit9",
    "ns1"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling of CL/F and Q/F (exponent e_wt_cl) and of Vc/F and Vp/F (exponent e_wt_vc), centred on 70 kg (Supplementary Methods 'Ivermectin PK modeling').",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear spline with a cut-off at 20 years on the death rate of virus-producing cells and on the baseline NS1 elimination rate; no effect above 20 years. The NONMEM $INPUT comment labels AGE as months, but the 20-year cut-off in the text and Tables 1-2 is in years.",
      source_name = "AGE"
    ),
    DENV_SEROTYPE1 = list(
      description = "Infecting dengue virus serotype 1 indicator (1 = DENV-1, 0 = other serotype)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (serotype 2 or 4; DENV_SEROTYPE1 = DENV_SEROTYPE3 = 0)",
      notes = "Derived from the source SERO column: DENV_SEROTYPE1 = (SERO == 1). Multiplies the log10 exponent of beta and p.",
      source_name = "SERO"
    ),
    DENV_SEROTYPE3 = list(
      description = "Infecting dengue virus serotype 3 indicator (1 = DENV-3, 0 = other serotype)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (serotype 2 or 4; DENV_SEROTYPE1 = DENV_SEROTYPE3 = 0)",
      notes = "Derived from the source SERO column: DENV_SEROTYPE3 = (SERO == 3). Multiplies the log10 exponent of beta and p.",
      source_name = "SERO"
    ),
    DENV_PRIMARY = list(
      description = "Primary dengue infection indicator (1 = primary, 0 = secondary infection)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (secondary infection)",
      notes = "Derived from the source TYPE column (1 = primary, 2 = secondary): DENV_PRIMARY = (TYPE == 1). Proportional effect on ET50 of NS1 secretion.",
      source_name = "TYPE"
    ),
    ABNS1_IGG = list(
      description = "NS1-specific IgG level measured on day 5 after fever onset",
      units = "relative units (normalised to pooled convalescent sera)",
      type = "continuous",
      reference_category = NULL,
      notes = "One measurement per patient (day 5 after fever onset), treated as time-fixed. The source dataset column IGGN carries the natural log; supply the linear level here (ABNS1_IGG = exp(IGGN)). Enters the NS1 elimination rate only above the cut-off ln(ABNS1_IGG) = 4, i.e. ABNS1_IGG = exp(4) = 54.6 relative units.",
      source_name = "IGGN"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ivermectin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ivermectin", units = "mg", specimen = "plasma", verified = TRUE),
    target = list(
      analyte = "uninfected target cells (monocytes)",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    eclipse1 = list(
      analyte = "infected, not yet virus-producing cells",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    eclipse2 = list(
      analyte = "infected, not yet virus-producing cells",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    eclipse3 = list(
      analyte = "infected, not yet virus-producing cells",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    eclipse4 = list(
      analyte = "infected, not yet virus-producing cells",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    infected = list(
      analyte = "virus-producing infected cells",
      units = "cells/mL",
      specimen = "blood cell",
      verified = TRUE
    ),
    virus = list(analyte = "dengue virus RNA", units = "copies/mL", specimen = "plasma", verified = TRUE),
    ns1_transit1 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit2 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit3 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit4 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit5 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit6 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit7 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit8 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1_transit9 = list(
      analyte = "NS1 secretion delay state",
      units = "cells*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ns1 = list(analyte = "dengue non-structural protein 1 (NS1)", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 248L,
    n_studies = 2L,
    age_range = "3.6-72.6 years (adults 15.1-72.6, children 3.6-15.1)",
    age_median = "23.2 years (adult ivermectin arm), 24.5 years (adult placebo arm), 10.1 years (children)",
    weight_range = "17.2-135 kg",
    weight_median = "60.0 kg (adults), 29.9 kg (children)",
    sex_female_pct = 42.0,
    race_ethnicity = "Thai",
    disease_state = "Acute dengue fever (laboratory-confirmed dengue virus infection), 87% secondary infections; serotypes 1-4",
    dose_range = "Oral ivermectin 400 ug/kg once daily for 2 or 3 days (adults) or 400 or 600 ug/kg once daily for 3 days (children); placebo arm in adults",
    regions = "Thailand (Siriraj Hospital, Bangkok; Loei Hospital)",
    notes = paste(
      "Pooled ESIDEN phase II/III placebo-controlled adult trial (NCT02045069;",
      "250 enrolled: 145 ivermectin, 105 placebo) and the PKIDEN paediatric",
      "PK-PD trial (NCT03432442; 24 children). Supplementary Table S2. Analysis",
      "sets: ivermectin PK 169 patients (145 adults, 24 children; 1,369",
      "concentrations); viral load 224 adults plus 24 children; dengue IgM/IgG",
      "248 patients; NS1 238 patients. Female percentage pooled across Table S2",
      "(115 of 274 enrolled)."
    )
  )

  ini({
    # Ivermectin PK (Supplementary Table S3; reference weight 70 kg)
    lcl <- log(14.0); label("Apparent clearance CL/F at 70 kg (L/h)") # Table S3 'CL/F (L/h)' = 14.0
    lvc <- log(356); label("Apparent central volume Vc/F at 70 kg (L)") # Table S3 'Vc/F (L)' = 356
    lq <- log(24.9); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table S3 'Q/F (L/h)' = 24.9
    lvp <- log(1140); label("Apparent peripheral volume Vp/F at 70 kg (L)") # Table S3 'Vp/F (L)' = 1,140
    lmtt <- log(2.87); label("Mean absorption transit time MTT (h)") # Table S3 'MTT (h)' = 2.87
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (fraction)") # Table S3 'F (%)' = 100 Fix
    e_wt_cl <- 0.947; label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Table S3 'Bodyweight on CL/F' = 0.947; NONMEM code applies WTCL to CL and Q
    e_wt_vc <- 1.33; label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Table S3 'Bodyweight on Vc/F' = 1.33; NONMEM code applies WTV to V2 and V3

    # Dengue viral kinetics (Table 1; NONMEM viral-kinetic code)
    lbeta <- log(10^-6.47); label("Infection rate constant beta (mL/copy/day)") # Table 1 'beta (mL/day/copy)' = 10^-6.47
    ldelta <- log(2.76); label("Death rate constant of virus-producing cells delta (1/day)") # Table 1 'delta (1/day)' = 2.76
    lp <- log(10^3.3); label("Virus production rate constant p (copies/cell/day)") # Table 1 'p (copy/mL/day)' = 10^3.3
    lcvir <- log(3.04); label("Innate-immunity virus elimination rate constant c (1/day)") # Table 1 'c (1/day)' = 3.04
    lkab <- log(4.58); label("Antibody-driven virus neutralisation rate constant kH (1/day)") # Table 1 'k H (1/day)' = 4.58
    lktr_ecl <- log(6.73); label("Transit rate constant from infected to virus-producing cells KTR (1/day)") # Table 1 'KTR (1/day)' = 6.73
    lt0 <- fixed(log(3.5e5)); label("Initial target-cell concentration T0 (cells/mL)") # Methods 'fixed at a value of 3.5 x 10^5 cells/mL'; NONMEM T0 = 3.5*(10**5) (Table 1 prints 3.5 x 10^6)
    lv0 <- fixed(log(1)); label("Initial viral load V0 (copies/mL)") # Table 1 'V0' = 1 (fixed); Methods '1 dengue RNA copy number/mL'
    e_sero1_beta <- 0.152; label("Serotype 1 effect on the log10 exponent of beta (fraction)") # Table 1 'Serotype 1 on beta (%)' = 15.2
    e_sero3_beta <- -0.128; label("Serotype 3 effect on the log10 exponent of beta (fraction)") # Table 1 'Serotype 3 on beta (%)' = -12.8
    e_sero1_p <- 0.356; label("Serotype 1 effect on the log10 exponent of p (fraction)") # Table 1 'Serotype 1 on p (%)' = 35.6
    e_sero3_p <- -0.290; label("Serotype 3 effect on the log10 exponent of p (fraction)") # Table 1 'Serotype 3 on p (%)' = -29.0
    e_age_delta <- 0.0252; label("Fractional change in delta per year of age below 20 years (1/year)") # Table 1 'Age on delta (%)' = 2.5; NONMEM THETA(13) = 0.0252; Supplementary Results '2.52% increase ... per year'

    # Dengue-specific IgM and IgG sigmoidal Emax time courses (Supplementary Table S4)
    lemax_igm <- log(99.5); label("Maximum dengue-specific IgM level (assay units)") # Table S4 IgM 'Emax' = 99.5
    let50_igm <- log(10.8); label("Time to half-maximal IgM ET50 (day)") # Table S4 IgM 'ET 50 (day)' = 10.8
    lhill_igm <- log(10.9); label("Hill coefficient of the IgM time course (unitless)") # Table S4 IgM 'Hill coefficient' = 10.9
    lemax_igg <- log(233); label("Maximum dengue-specific IgG level (assay units)") # Table S4 IgG 'Emax' = 233
    let50_igg <- log(11.5); label("Time to half-maximal IgG ET50 (day)") # Table S4 IgG 'ET 50 (day)' = 11.5
    lhill_igg <- log(8.8); label("Hill coefficient of the IgG time course (unitless)") # Table S4 IgG 'Hill coefficient' = 8.8

    # Dengue NS1 kinetics and ivermectin effect (Table 2)
    lktr_ns1 <- log(6.76); label("Transit rate constant of NS1 secretion KTR (1/day)") # Table 2 'KTR (1/day)' = 6.76
    lemax_ns1 <- log(1.04); label("Maximum NS1 secretion rate constant Emax,NS1 (1/day)") # Table 2 'E max,NS1' = 1.04
    let50_ns1 <- log(8.3); label("Time to half-maximal NS1 secretion ET50,NS1 (day)") # Table 2 'ET50,NS1 (days)' = 8.3
    lgam_ns1 <- log(18.0); label("Shape factor of the NS1 secretion time course (unitless)") # Table 2 'gamma NS1' = 18.0
    lkbase_ns1 <- log(0.479); label("Baseline NS1 elimination rate constant kbase (1/day)") # Table 2 'k base' = 0.479
    e_ns1igg_kb <- 0.112; label("Slope of NS1 elimination on time x (ln NS1-specific IgG - 4) (1/day^2)") # Table 2 'NS1-specific IgG on kB' = 0.112
    e_primary_et50 <- 0.231; label("Primary-infection effect on ET50,NS1 (fraction)") # Table 2 'Primary infection on ET50,NS1 (%)' = 23.1
    e_age_kbase <- -0.301; label("Fractional change in kbase per year of age below 20 years (1/year)") # Table 2 'Age on k base (%)' = -30.1
    emax_ivm <- fixed(0.99); label("Maximum inhibition of NS1 secretion by ivermectin (fraction)") # Table 2 'Maximum effect of ivermectin on E max,NS1' = 0.99 (fixed)
    lec50_ivm <- log(67.5); label("Ivermectin concentration giving half-maximal NS1 inhibition EC50 (ng/mL)") # Table 2 'EC50 of ivermectin concentration (ng/mL)' = 67.5

    # IIV: variances from the reported CV% as omega^2 = log(1 + CV^2)
    # (Table S3 / S4 notes: CV = 100 x sqrt(exp(variance) - 1))
    etalcl ~ log(1 + 0.246^2) # Table S3 CL/F CV 24.6%
    etalq ~ log(1 + 0.457^2) # Table S3 Q/F CV 45.7%
    etalmtt ~ log(1 + 0.316^2) # Table S3 MTT CV 31.6%
    etalfdepot ~ log(1 + 0.488^2) # Table S3 F CV 48.8%
    etalbeta ~ log(1 + 0.505^2) # Table 1 beta CV 50.5% (NONMEM $OMEGA 0.227)
    etaldelta ~ log(1 + 0.349^2) # Table 1 delta CV 34.9% (NONMEM $OMEGA 0.115)
    etalemax_igm ~ log(1 + 0.622^2) # Table S4 IgM Emax CV 62.2%
    etalet50_igm ~ log(1 + 0.089^2) # Table S4 IgM ET50 CV 8.9%
    etalhill_igm ~ log(1 + 0.209^2) # Table S4 IgM Hill CV 20.9%
    etalet50_igg ~ log(1 + 0.242^2) # Table S4 IgG ET50 CV 24.2%
    etalhill_igg ~ log(1 + 0.510^2) # Table S4 IgG Hill CV 51.0%
    etalktr_ns1 ~ log(1 + 0.379^2) # Table 2 KTR CV 37.9%
    etalet50_ns1 ~ log(1 + 0.134^2) # Table 2 ET50,NS1 CV 13.4%
    etalkbase_ns1 ~ log(1 + 0.903^2) # Table 2 k base CV 90.3%

    # Residual error: additive on log-transformed observations
    expSd <- 0.251; label("Residual error on log ivermectin concentration (SD, log scale)") # Table S3 'Residual error (%)' = 25.1
    expSd_viral_load <- sqrt(1.31); label("Residual error on log viral load (SD, log scale)") # Table 1 'Residual error' = 1.31 = NONMEM $SIGMA (variance); SD = sqrt(1.31)
    expSd_ns1 <- sqrt(1.15); label("Residual error on log NS1 (SD, log scale)") # Table 2 'Residual error' = 1.15 = NONMEM $SIGMA (variance); SD = sqrt(1.15)
    expSd_igm_denv <- 0.426; label("Residual error on log dengue-specific IgM (SD, log scale)") # Table S4 IgM 'RUV (%CV)' = 42.6
    expSd_igg_denv <- 0.273; label("Residual error on log dengue-specific IgG (SD, log scale)") # Table S4 IgG 'RUV (%CV)' = 27.3
  })

  model({
    # Ivermectin PK; parameters are reported per hour and converted to per day
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * 24
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl * 24
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    mtt <- exp(lmtt + etalmtt) / 24
    # Dose compartment plus six transit compartments: ktr = (NN + 1) / MTT with NN = 6
    ktr <- 7 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Dengue-specific IgM / IgG; t is time since infection (day)
    emax_igm <- exp(lemax_igm + etalemax_igm)
    et50_igm <- exp(let50_igm + etalet50_igm)
    hill_igm <- exp(lhill_igm + etalhill_igm)
    emax_igg <- exp(lemax_igg)
    et50_igg <- exp(let50_igg + etalet50_igg)
    hill_igg <- exp(lhill_igg + etalhill_igg)
    frac_igm <- t^hill_igm / (et50_igm^hill_igm + t^hill_igm)
    frac_igg <- t^hill_igg / (et50_igg^hill_igg + t^hill_igg)
    igm_denv <- emax_igm * frac_igm
    igg_denv <- emax_igg * frac_igg
    # Normalised antibody response A(t) = (IgM / Emax_IgM + IgG / Emax_IgG) / 2
    ab_norm <- (frac_igm + frac_igg) / 2

    # Viral kinetics; serotype effects multiply the log10 exponent of beta and
    # p, as in the NONMEM code (TVBETA = THETA(1) * CSERO1; BETA = 10**TVBETA)
    beta <- exp(lbeta * (1 + e_sero1_beta * DENV_SEROTYPE1 + e_sero3_beta * DENV_SEROTYPE3) + etalbeta)
    p <- exp(lp * (1 + e_sero1_p * DENV_SEROTYPE1 + e_sero3_p * DENV_SEROTYPE3))
    delta <- exp(ldelta + etaldelta) * (1 + e_age_delta * (min(AGE, 20) - 20))
    cvir <- exp(lcvir)
    kab <- exp(lkab)
    ktr_ecl <- exp(lktr_ecl)
    t0 <- exp(lt0)
    v0 <- exp(lv0)

    # NS1 kinetics
    ktr_ns1 <- exp(lktr_ns1 + etalktr_ns1)
    emax_ns1 <- exp(lemax_ns1)
    et50_ns1 <- exp(let50_ns1 + etalet50_ns1) * (1 + e_primary_et50 * DENV_PRIMARY)
    gam_ns1 <- exp(lgam_ns1)
    kp_ns1 <- emax_ns1 * t^gam_ns1 / (t^gam_ns1 + et50_ns1^gam_ns1)
    # Elimination rises linearly with time once ln(NS1-specific IgG) exceeds 4
    kb_ns1 <- (exp(lkbase_ns1) * (1 + e_age_kbase * (min(AGE, 20) - 20)) +
      e_ns1igg_kb * t * max(log(ABNS1_IGG) - 4, 0)) * exp(etalkbase_ns1)
    ec50_ivm <- exp(lec50_ivm)

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(central) <- ktr * transit6 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L * 1000 = ng/mL
    Cc <- central / vc * 1000
    inh_ivm <- 1 - emax_ivm * Cc / (Cc + ec50_ivm)

    d/dt(target) <- -beta * target * virus
    d/dt(eclipse1) <- beta * target * virus - ktr_ecl * eclipse1
    d/dt(eclipse2) <- ktr_ecl * eclipse1 - ktr_ecl * eclipse2
    d/dt(eclipse3) <- ktr_ecl * eclipse2 - ktr_ecl * eclipse3
    d/dt(eclipse4) <- ktr_ecl * eclipse3 - ktr_ecl * eclipse4
    d/dt(infected) <- ktr_ecl * eclipse4 - delta * infected
    d/dt(virus) <- p * infected - cvir * virus - kab * ab_norm * virus

    # NS1 delay chain is fed by (not drained from) the producing cells
    d/dt(ns1_transit1) <- infected - ktr_ns1 * ns1_transit1
    d/dt(ns1_transit2) <- ktr_ns1 * ns1_transit1 - ktr_ns1 * ns1_transit2
    d/dt(ns1_transit3) <- ktr_ns1 * ns1_transit2 - ktr_ns1 * ns1_transit3
    d/dt(ns1_transit4) <- ktr_ns1 * ns1_transit3 - ktr_ns1 * ns1_transit4
    d/dt(ns1_transit5) <- ktr_ns1 * ns1_transit4 - ktr_ns1 * ns1_transit5
    d/dt(ns1_transit6) <- ktr_ns1 * ns1_transit5 - ktr_ns1 * ns1_transit6
    d/dt(ns1_transit7) <- ktr_ns1 * ns1_transit6 - ktr_ns1 * ns1_transit7
    d/dt(ns1_transit8) <- ktr_ns1 * ns1_transit7 - ktr_ns1 * ns1_transit8
    d/dt(ns1_transit9) <- ktr_ns1 * ns1_transit8 - ktr_ns1 * ns1_transit9
    d/dt(ns1) <- kp_ns1 * inh_ivm * ns1_transit9 - kb_ns1 * ns1

    f(depot) <- exp(lfdepot + etalfdepot)
    target(0) <- t0
    virus(0) <- v0

    viral_load <- virus

    Cc ~ lnorm(expSd)
    viral_load ~ lnorm(expSd_viral_load)
    ns1 ~ lnorm(expSd_ns1)
    igm_denv ~ lnorm(expSd_igm_denv)
    igg_denv ~ lnorm(expSd_igg_denv)
  })
}
