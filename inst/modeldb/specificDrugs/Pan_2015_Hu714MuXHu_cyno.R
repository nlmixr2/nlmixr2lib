Pan_2015_Hu714MuXHu_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Population PK/PD model for Hu714MuXHu, a chimeric murine-human anti-IL-15 IgG1 monoclonal antibody (cynomolgus surrogate of AMG 714), fitted simultaneously to serum concentrations and peripheral NK cell counts (Pan 2015): two-compartment PK with linear elimination, first-order SC absorption and estimated SC bioavailability, linked to an indirect-response model in which serum Hu714MuXHu inhibits NK cell production (Emax). Parameters are body-weight normalised (per kg) and doses are in mg/kg."
  reference <- "Pan WJ, Li H, Xiao JJ, Horner MJ, Lebrec HN, Butz EA, Kaliyaperumal A, Cheah TC, Ortiz RC, Prokop SP, Buntich SA, Boren BM, Wolford ST, Tsuji WH, Wienkers LC, Kock K. Modeling the pharmacokinetic-pharmacodynamic relationship of the monoclonal anti-macaque-IL-15 antibody Hu714MuXHu in cynomolgus monkeys. Pharmacol Res Perspect. 2015;3(6):e00199. doi:10.1002/prp2.199. PMCID: PMC4777250."
  vignette <- "Pan_2015_Hu714MuXHu_cyno"
  dosing <- c("depot", "central")
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Peripheral-blood NK cell count is the PD turnover state; it has no
  # canonical register entry (Roepcke_2018_tak_079.R also whitelists it).
  paper_specific_compartments <- c("nkcell")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "Hu714MuXHu", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "Hu714MuXHu", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "Hu714MuXHu", units = "mg/kg", specimen = "tissue", verified = TRUE),
    # The NK state holds a COUNT per volume of blood (Methods: 'R is the number
    # of NK cells x 10^3/uL in blood'), with kin = R0 * kout and R(0) = R0.
    nkcell = list(analyte = "NK cells (CD3- CD16+)", units = "10^3 cells/uL", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis, Chinese origin)",
    n_subjects = 104L,
    n_studies = 3L,
    age_range = "3-6 years",
    weight_range = "2.0-4.6 kg",
    sex_female_pct = 47.1,
    disease_state = "Healthy, treatment-naive monkeys (single-dose PK/PD study and two GLP toxicology studies).",
    dose_range = "Single dose 0.1 or 1 mg/kg IV; 0 (vehicle), 30, 60 or 150 mg/kg IV or 150 mg/kg SC once weekly for 5 doses; 0 (vehicle), 5, 30 or 150 mg/kg SC once weekly for 14 doses.",
    regions = NA_character_,
    notes = "Methods (Animal husbandry): 55 male and 49 female monkeys. Table 1: single-dose study 6 males (3 per group); 1-month GLP study 25 males + 25 females (5 per sex per group, 2 per sex per group kept for recovery to Day 337); 3-month GLP study 24 males + 24 females (6 per sex per group, 2 per sex per group kept for recovery to Day 205). 1408 serum concentrations and 982 NK cell counts were fitted simultaneously in NONMEM 7 (ADVAN13, FOCEI then IMP). No anti-drug antibodies were detected."
  )

  ini({
    # Table 3 reports the PK parameters per kg of body weight in mL; they are
    # converted to L/kg here so that mg/kg doses give Cc in mg/L = ug/mL.
    lka <- log(0.522); label("First-order SC absorption rate constant, Ka (1/day)") # Table 3: Ka 0.522 /day (RSE 7.11%)
    lcl <- log(3.66 / 1000); label("Linear clearance, CL (L/day/kg)") # Table 3: CL 3.66 mL/day/kg (RSE 3.22%)
    lvc <- log(44.8 / 1000); label("Central volume, Vc (L/kg)") # Table 3: Vc 44.8 mL/kg (RSE 4.51%)
    lq <- log(5.07 / 1000); label("Intercompartmental clearance, Q (L/day/kg)") # Table 3: Q 5.07 mL/day/kg (RSE 17.9%)
    lvp <- log(24.3 / 1000); label("Peripheral volume, Vp (L/kg)") # Table 3: Vp 24.3 mL/kg (RSE 8.25%)
    lfdepot <- log(0.824); label("SC bioavailability, F1 (fraction)") # Table 3: F1 0.824 (RSE 3.65%), no BSV

    lrbase <- log(0.584); label("Baseline NK cell count, R0 (10^3 cells/uL)") # Table 3: R0 0.584 x 10^3 cells/uL (RSE 4.96%)
    lkout <- log(0.305); label("First-order NK cell loss rate constant, Kout (1/day)") # Table 3: Kout 0.305 /day (RSE 14.9%)
    lemax <- log(0.941); label("Maximal fractional inhibition of NK cell production, Emax (fraction)") # Table 3: Emax 0.941 (RSE 0.386%); unit column misprinted as 10^3 cells/uL
    lec50 <- log(0.0978); label("Hu714MuXHu concentration at half-maximal inhibition, EC50 (ug/mL)") # Table 3: EC50 0.0978 ug/mL (RSE 28.0%)

    # Table 3 reports BSV as a percentage without stating the formula;
    # converted with the log-normal CV relation omega^2 = log(1 + CV^2).
    # A CL-Vc correlation was 'considered' (Methods) but none is reported.
    etalka ~ log(1 + 0.237^2) # Table 3: Ka BSV 23.7%
    etalcl ~ log(1 + 0.171^2) # Table 3: CL BSV 17.1%
    etalvc ~ log(1 + 0.152^2) # Table 3: Vc BSV 15.2%
    etalq ~ log(1 + 0.721^2) # Table 3: Q BSV 72.1%
    etalvp ~ log(1 + 0.126^2) # Table 3: Vp BSV 12.6%
    etalrbase ~ log(1 + 0.404^2) # Table 3: R0 BSV 40.4%
    etalkout ~ log(1 + 0.632^2) # Table 3: Kout BSV 63.2%
    etalemax ~ log(1 + 0.0245^2) # Table 3: Emax BSV 2.45%
    etalec50 ~ log(1 + 0.907^2) # Table 3: EC50 BSV 90.7%

    # Residual error Y = F * (1 + eps1) + eps2; Table 3 prints the variances,
    # and both additive variances are fixed to 0 (Methods), so only the
    # proportional components are kept.
    propSd <- sqrt(0.0349); label("Proportional residual error for Hu714MuXHu (fraction)") # Table 3: sigma1^2 (prop. PK) 0.0349
    propSd_nkcell <- sqrt(0.217); label("Proportional residual error for NK cell count (fraction)") # Table 3: sigma3^2 (prop. PD) 0.217
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    fdepot <- exp(lfdepot)

    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    kin <- rbase * kout # steady state before the first dose

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # Figure 1: SC doses enter depot, IV doses enter central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Methods equation: dR/dt = Kin * (1 - Emax * C / (EC50 + C)) - Kout * R
    d/dt(nkcell) <- kin * (1 - emax * Cc / (ec50 + Cc)) - kout * nkcell
    nkcell(0) <- rbase
    f(depot) <- fdepot

    Cc ~ prop(propSd)
    nkcell ~ prop(propSd_nkcell)
  })
}
