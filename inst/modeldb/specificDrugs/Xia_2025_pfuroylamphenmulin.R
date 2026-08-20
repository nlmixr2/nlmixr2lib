Xia_2025_pfuroylamphenmulin <- function() {
  description <- "Preclinical (chicken, specific-pathogen-free, Mycoplasma gallisepticum-infected). Two-compartment first-order-absorption plasma PK model for the novel pleuromutilin derivative p-furoylamphenmulin after intramuscular injection in M. gallisepticum-infected chicks, coupled to the paper's inhibitory sigmoid Emax PK/PD-index model for the reduction in lung mycoplasma load. PK from Xia 2025 Table 1 (mean of the 5, 40 and 80 mg/kg dose-group WinNonlin 5.2.1 fits): T1/2ka 0.0961 h, T1/2kel 4.2281 h, V1/F 6.3427 L/kg, V2/F 2.2405 L/kg, Cl/F 3.8692 L/h/kg. Xia 2025 reports no intercompartmental clearance, but the five reported quantities determine the two-compartment system exactly, so lq is DERIVED algebraically from them (see the ini() comment and the vignette). PD from Xia 2025 Table 2, AUC24h/MIC column: the 72 h reduction in lung load is E = E(0) + (Emax - E(0)) * Ce^N / (EC50^N + Ce^N) with E(0) = 1.07 log10 CFU/mL at zero exposure, a maximum of 4.29 log10 CFU/mL, EC50 = 7526.81 h and Hill N = 4.22. NOTE the paper's Table 2 row labels are inverted relative to the usual convention: its 'Emax' row (1.07) is the zero-exposure effect and its 'E0' row (4.29) is the maximum effect, as Xia 2025 defines them in Materials and methods and confirms in Results. Substituting Table 2 into the equation at the paper's own 3 log10 CFU/mL target of AUC24h/MIC = 8288.29 h returns E = 3.003, confirming the role assignment. The PK/PD index is formed from the closed form AUC24h = dose / (Cl/F), which is how Xia 2025 itself defines Cl/F (Dose/AUC24h reproduces every Table 1 AUC to nine decimal places) and how its dose-calculation formula (2) inverts it to the recommended 62.64 mg/kg. The lung mycoplasma density bact (linear CFU/mL) is integrated as d/dt(bact) = -ln(10) * (E / 72) * bact so that log10(bact) falls by exactly E across the paper's 72 h treatment course, reproducing the endpoint model at the single time Xia 2025 actually counted bacteria. Xia 2025 reports neither between-subject variability nor a residual error magnitude for either endpoint, so no eta parameters are present and both residual SDs are FIXED at 0 for typical-value simulation. The companion Cmax/MIC parameterisation of Table 2 is not packaged as a separate file because a running Cmax is not an ODE quantity; it is reproduced in full in the vignette."
  reference <- paste(
    "Xia X, Zhao H, Li Y, Long X, Liu X, Bai M, Tang Y, Shen X, Ding H. (2025).",
    "Pharmacokinetic/pharmacodynamic relationship of a novel pleuromutilin",
    "derivative p-furoylamphenmulin against Mycoplasma gallisepticum in vivo",
    "in chickens.",
    "Poultry Science 104(8):105249.",
    "doi:10.1016/j.psj.2025.105249.",
    sep = " "
  )
  vignette <- "Xia_2025_pfuroylamphenmulin"
  units <- list(
    time = "h",
    dosing = "mg/kg (p-furoylamphenmulin, intramuscular)",
    concentration = "ug/mL (p-furoylamphenmulin, Cc); log10 CFU/mL (M. gallisepticum lung load, log10cfu)"
  )

  # auc_0_24 -- p-furoylamphenmulin plasma AUC accumulated over the FIRST 24 h
  #   dosing interval only; carried so the vignette can check the closed-form
  #   index used below against the integrated exposure. Windowed-accumulator
  #   idiom follows Beredaki_2023_micafungin_clsi / Kuchimanchi_2018_evolocumab_ldlc.
  # bact -- linear M. gallisepticum density in lung homogenate (CFU/mL).
  paper_specific_compartments <- c("auc_0_24", "bact")

  covariateData <- list()

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "p-furoylamphenmulin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "p-furoylamphenmulin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "p-furoylamphenmulin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    auc_0_24    = list(analyte = "p-furoylamphenmulin", units = "ug*h/mL", specimen = "plasma", verified = TRUE),
    bact        = list(analyte = "Mycoplasma gallisepticum strain S6 (ATCC 15302)", units = "CFU/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "chicken (specific-pathogen-free, 1 day old at purchase, 35-45 g)",
    n_subjects     = 368L,
    n_studies      = 1L,
    age_range      = "6 days old at the start of the pharmacokinetic study (3 days acclimatisation, then 3 daily intratracheal inoculations)",
    weight_range   = "35-45 g at one day of age",
    organism       = "Mycoplasma gallisepticum standard strain S6 (ATCC 15302; China Institute of Veterinary Drug Control). p-furoylamphenmulin MIC = 0.001953125 ug/mL by broth microdilution at a 10^6 CFU/mL inoculum; tiamulin MIC against the same strain = 0.0078125 ug/mL",
    system         = "In vivo M. gallisepticum respiratory-infection model: 0.2 mL of a 10^9 CFU/mL suspension inoculated intratracheally once daily for three consecutive days through a sterile gavage needle. Lung colonisation at 0 h was 5.88 +/- 0.03 log10 CFU/mL in the infection-model validation cohort and 6.13 log10 CFU/mL in the efficacy cohort",
    disease_state  = "Experimental chronic respiratory disease (pneumonia and air sacculitis) from intratracheal M. gallisepticum challenge; air-sac lesion score 2.25 +/- 0.25 in the untreated infected control",
    dose_range     = "PK: single intramuscular injections of 5, 40 and 80 mg/kg (288 infected chickens, n = 8 per time point per dose). PD: intramuscular injections of 0, 5, 10, 20, 30, 40, 50, 60, 70 or 80 mg/kg once daily for 3 consecutive days (80 infected chickens, 8 per group), plus tiamulin fumarate 40 mg/kg comparator arms given orally and intramuscularly",
    design         = "Vehicle 10% DMSO / 10% Tween 80 / 80% saline. Blood was drawn once per chicken only (0.5 mL by cardiac puncture), so the plasma PK is a destructive-sampling naive-pooled design with 8 chickens contributing at each nominal time",
    sampling       = "PK: plasma at 0.083, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 10, 12 and 24 h after the single dose (n = 8 per time point). PD: lungs collected, homogenised and plated 24 h after the last of three daily administrations, i.e. 72 h after the first dose",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Animal Ethics Committee of South China Agricultural University approval 2024c019. Xia 2025 fitted each dose group separately in WinNonlin 5.2.1 and reported a Mean +/- SEM column; the packaged values are that mean column, which is also the parameter set the paper's own dose calculation uses (Cl/F = 3.87 L/h/kg -> 62.64 mg/kg). The per-dose-group parameter sets, at the full precision of supplementary file mmc5.xlsx, are (5 / 40 / 80 mg/kg): T1/2kel 4.205003 / 3.607271 / 4.871964 h; T1/2ka 0.08902 / 0.124374 / 0.074923 h; V1/F 7.392317 / 4.855187 / 6.780696 L/kg; V2/F 1.488610 / 2.113411 / 3.119444 L/kg; Cl/F 3.354635 / 3.941475 / 4.311523 L/h/kg; giving derived Q/F of 0.27723 / 0.46947 / 0.51166 L/h/kg (mean 0.41945, against 0.42210 derived from the mean column). Cl/F rises and AUC24h therefore falls short of dose proportionality across this range, so a mean-parameter model over-predicts AUC24h at 80 mg/kg by about 11% and under-predicts it at 5 mg/kg by about 13%; see the vignette. Xia 2025 Table 2 also reports the same exposure-response fitted against Cmax/MIC (zero-exposure effect 1.09, maximum 4.24 log10 CFU/mL, EC50 3945.62, Hill N 5.09, R^2 0.84, 3 log10 CFU/mL target 4299.30); AUC24h/MIC (R^2 0.83) is the index the paper designates as most important and the one packaged here. Xia 2025 deliberately did not use %T>MIC because it was 89.53% at 5 mg/kg and 100% at every higher dose and so could not discriminate the observed effects."
  )

  ini({
    # =================================================================
    # Plasma pharmacokinetics -- Xia 2025 Table 1, "Mean +/- SEM" column
    # =================================================================
    # Xia 2025 Materials and methods, "Pharmacokinetic data analysis":
    # "the first-order absorption two-compartment model was used to
    # calculate the main PK parameters" in WinNonlin 5.2.1. The Table 1
    # mean column is reported to two decimal places; the values below use
    # the full precision of supplementary file mmc5.xlsx, whose summary
    # block is the source of the printed Table 1.
    lka <- log(log(2) / 0.0961056666666667)
    label("Log first-order absorption rate constant ka (1/h)")  # Xia 2025 Table 1, T1/2ka = 0.10 +/- 0.01 h (mmc5.xlsx mean 0.0961056666666667 h); ka = ln(2)/T1/2ka = 7.2123 1/h
    lcl <- log(3.86921107)
    label("Log apparent clearance Cl/F (L/h/kg)")  # Xia 2025 Table 1, Cl/F = 3.87 +/- 0.28 L/h/kg (mmc5.xlsx mean 3.86921107)
    lvc <- log(6.3427335)
    label("Log apparent central volume V1/F (L/kg)")  # Xia 2025 Table 1, V1/F = 6.34 +/- 0.76 L/kg (mmc5.xlsx mean 6.3427335)
    lvp <- log(2.240488147)
    label("Log apparent peripheral volume V2/F (L/kg)")  # Xia 2025 Table 1, V2/F = 2.24 +/- 0.48 L/kg (mmc5.xlsx mean 2.240488147)

    # DERIVED, not printed. Xia 2025 reports no intercompartmental
    # clearance (nor k12/k21), but the four quantities above plus the
    # reported terminal half-life determine the two-compartment system
    # exactly, so Q/F follows algebraically with no free choice:
    #
    #   k10 = Cl/V1 = 0.6100226 1/h
    #   beta = ln(2) / T1/2kel = ln(2) / 4.22807933 = 0.16393908 1/h
    #   beta is the smaller root of
    #     beta^2 - (k10 + k12 + k21) beta + k10 * k21 = 0
    #   with k12 = Q/V1 and k21 = Q/V2, hence
    #     Q = beta (k10 - beta) / (k10/V2 - beta (1/V1 + 1/V2))
    #       = 0.42209916 L/h/kg
    #
    # Three independent checks that T1/2kel is the TERMINAL (beta)
    # half-life, which is what makes the derivation well posed:
    #   (i) k10 = 0.6100 1/h exceeds ln(2)/4.2281 = 0.1639 1/h, and in a
    #       two-compartment model alpha > k10 > beta, so 4.23 h can be
    #       neither the alpha half-life nor ln(2)/k10;
    #  (ii) ln(2)/(Cl/V1) is 1.53 / 0.86 / 1.09 h for the three dose
    #       groups, nothing like the reported 4.21 / 3.61 / 4.87 h;
    # (iii) the log-linear slope of the observed mean concentrations
    #       between 12 and 24 h (mmc5.xlsx) gives 4.13 / 3.42 / 4.27 h,
    #       matching the reported T1/2kel per dose group.
    # Substituting Cl/F, V1/F, V2/F and this Q/F back into the
    # two-compartment eigenvalue returns a terminal half-life of
    # 4.228079 h, i.e. the reported value exactly. See the vignette
    # "Source trace" and "Assumptions and deviations".
    lq <- log(0.42209916)
    label("Log apparent intercompartmental clearance Q/F (L/h/kg)")  # DERIVED from Xia 2025 Table 1 Cl/F, V1/F, V2/F and T1/2kel = 4.2281 h; not printed in the paper

    # =================================================================
    # PK/PD index integration window
    # =================================================================
    # FIXED experimental design input. Xia 2025 indexes the effect on
    # AUC24h/MIC and dosed once every 24 h for three days.
    tau <- fixed(24)
    label("Dosing interval and PK/PD index integration window (h;, design)")  # Xia 2025 Abstract and Dose calculation: "intramuscular injection once every 24 h for 3 consecutive days"; the index is AUC24h/MIC

    # =================================================================
    # Inhibitory sigmoid Emax exposure-response
    # Xia 2025 Table 2, AUC24h/MIC column
    # =================================================================
    # Xia 2025 Materials and methods, "PK/PD fitting analysis" fits the
    # Inhibitory Sigmoid Emax model of WinNonlin to the relationship
    # between the PK/PD index Ce and the mycoplasma reduction E, where
    # "E = antimicrobial effect, the change in mycoplasma load in the
    # lungs during 72 h of treatment in each group" (a positive
    # REDUCTION magnitude, in log10 CFU/mL).
    #
    # Xia 2025 uses INVERTED row labels relative to the usual convention,
    # and says so explicitly: "Emax = mycoplasma load change in the
    # positive control group during 72 h; E0 = maximum antimicrobial
    # effect reached after 72 h of treatment". Results repeats it:
    # "When the dose was 0 mg/kg, the load of mycoplasma in chicken lungs
    # naturally decreased by 1.07-1.09 Log10CFU/mL; E0 was 4.29 and 4.24
    # Log10CFU/mL". So the Table 2 'Emax' row is the ZERO-EXPOSURE effect
    # and the Table 2 'E0' row is the MAXIMUM effect. The parameters
    # below are named for their ROLE (e0 = zero-exposure, emax = maximum),
    # matching Chen_2023_tilmicosin and Beredaki_2023_micafungin_clsi.
    #
    # Numerical confirmation of the role assignment and the equation
    # form: substituting the four values below at the paper's own
    # 3 log10 CFU/mL target of AUC24h/MIC = 8288.29 h returns E = 3.003,
    # and the Cmax/MIC parameterisation at its target of 4299.30 returns
    # E = 3.004. Neither holds under any other assignment of the rows.
    e0 <- 1.07
    label("Reduction in lung mycoplasma load over 72 h at zero exposure (log10 CFU/mL)")  # Xia 2025 Table 2, AUC24h/MIC column, row labelled "Emax" = 1.07; Results: "When the dose was 0 mg/kg, the load of mycoplasma in chicken lungs naturally decreased by 1.07-1.09 Log10CFU/mL"
    lemax <- log(4.29)
    label("Log maximum reduction in lung mycoplasma load over 72 h Emax (log10 CFU/mL)")  # Xia 2025 Table 2, AUC24h/MIC column, row labelled "E0" = 4.29; Materials and methods: "E0 = maximum antimicrobial effect reached after 72 h of treatment"
    lec50 <- log(7526.81)
    label("Log AUC24h/MIC producing 50% of the maximum effect EC50 (h)")  # Xia 2025 Table 2, AUC24h/MIC column, EC50 = 7526.81
    lhill <- log(4.22)
    label("Log Hill coefficient N defining the slope of the fitted curve (unitless)")  # Xia 2025 Table 2, AUC24h/MIC column, Hill's slope = 4.22

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # The PK/PD index driving the sigmoid is AUC24h/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of
    # the challenge strain, not an estimated parameter. Change it to
    # apply the model to an isolate with a different susceptibility;
    # isolates recovered from the 60-80 mg/kg groups had MICs of
    # 0.00390625 ug/mL, one two-fold dilution higher.
    mic <- fixed(0.001953125)
    label("p-furoylamphenmulin MIC against M. gallisepticum S6 (ug/mL;, measured)")  # Xia 2025 Results, "Susceptibility changes": "For the original strain S6, the MIC of p-furoylamphenmulin was 0.001953125 ug/mL"; supplementary file mmc3.xlsx

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # FIXED experimental design input, not an estimated parameter. This
    # is the 0 h lung load of the EFFICACY cohort, which is the baseline
    # every reported reduction E is measured against (supplementary file
    # mmc4.xlsx computes each chicken's reduction as 6.13035274875 minus
    # its own 72 h count). The 5.88 +/- 0.03 log10 CFU/mL quoted in
    # Results is the separate infection-model validation cohort
    # (supplementary file mmc2.xlsx).
    log10_cfu0 <- fixed(6.13035274875)
    label("Log10 starting M. gallisepticum density in lung homogenate (log10 CFU/mL)")  # Xia 2025 supplementary file mmc4.xlsx, mean 0 h count of the efficacy cohort = 6.13035274875 log10 CFU/mL (n = 8)

    # =================================================================
    # Residual error
    # =================================================================
    # Xia 2025 reports no residual standard deviation for either
    # endpoint. For the plasma concentrations it reports only the
    # bioanalytical performance (recoveries 94.47-110.17%, intra-assay
    # precision 3.56-13.82%, inter-assay precision 2.10-6.30%), which
    # characterises the assay rather than the fit; for the sigmoid it
    # reports only R^2 = 0.83. Both residual SDs are therefore held at
    # zero for deterministic typical-value simulation. See the vignette
    # "Assumptions and deviations".
    propSd <- fixed(0)
    label("Proportional residual error on p-furoylamphenmulin concentration (0; not reported in Xia 2025)")  # Xia 2025 reports bioanalytical recovery and precision only, no residual error for the PK fit
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Xia 2025)")  # Xia 2025 Table 2 reports R^2 = 0.83 only, no residual SD
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---------------------------------------------------------------
    # Plasma pharmacokinetics. Xia 2025 Materials and methods,
    # "Pharmacokinetic data analysis": first-order absorption
    # two-compartment model, intramuscular injection into the thigh.
    # Bioavailability is not identifiable from an extravascular-only
    # design, so every disposition parameter carries the paper's /F and
    # F is implicitly 1; dose the `depot` compartment in mg/kg.
    # ---------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc

    # p-furoylamphenmulin AUC accumulated over the FIRST dosing interval
    # only -- the integrated counterpart of the closed-form index below.
    # The t < tau gate freezes the state at AUC0-24 so the vignette can
    # compare the two.
    d/dt(auc_0_24) <- Cc * (t < tau)

    # ---------------------------------------------------------------
    # PK/PD index. Xia 2025 Table 1 defines Cl/F such that
    # Dose / (Cl/F) reproduces every reported AUC24h to nine decimal
    # places (5/3.354635278 = 1.490474995 = the reported AUC24h, and
    # likewise at 40 and 80 mg/kg), and its dose-calculation formula (2),
    # Dose = (AUC/MIC)breakpoint x MIC90 x Cl_per_hour / (fu x F),
    # inverts exactly this relation to reach the recommended
    # 62.64 mg/kg. The closed form is used rather than the integrated
    # state because Xia 2025 fitted a 72 h ENDPOINT against the completed
    # AUC24h, so the index must be available from t = 0 for the
    # trajectory below to land on the fitted endpoint. `podo(depot)` is
    # the amount of the most recent depot dose in mg/kg, so the index
    # follows the actual dosing record and cannot drift out of step with
    # it. NOTE: simulate the untreated control arm with an explicit
    # `amt = 0` depot dose record rather than with no dose record at
    # all -- `podo()` is NA until a dose has been seen, whereas an
    # amt = 0 record makes it 0 and the trajectory reduces to the
    # zero-exposure effect e0 exactly.
    # ---------------------------------------------------------------
    auc24 <- podo(depot) / cl
    aucmic <- auc24 / mic

    # ---------------------------------------------------------------
    # Inhibitory sigmoid Emax exposure-response (Xia 2025 Materials and
    # methods, "PK/PD fitting analysis"; parameters from Table 2).
    # kill72 is the paper's E: the REDUCTION in lung mycoplasma load
    # (log10 CFU/mL) accrued over the 72 h treatment course, rising from
    # e0 at zero exposure towards emax at saturating exposure.
    # ---------------------------------------------------------------
    kill72 <- e0 + (emax - e0) * aucmic^hill / (ec50^hill + aucmic^hill)

    # Xia 2025 counted lung mycoplasma once, 24 h after the last of three
    # daily doses, i.e. at 72 h. Spreading the fitted 72 h reduction
    # uniformly across that window makes log10(bact) fall by exactly
    # kill72 between t = 0 and t = 72 h, so the trajectory matches the
    # fitted endpoint model at the time the paper fitted it:
    #   d(log10 N)/dt = -kill72 / 72
    #   => d(N)/dt     = -ln(10) * (kill72 / 72) * N
    d/dt(bact) <- -log(10) * (kill72 / 72) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1-CFU/mL floor so the log stays finite if bact
    # is driven below 1 CFU/mL (matches the antibacterial-PD convention
    # used by Chen_2023_tilmicosin and Beredaki_2023_micafungin_clsi).
    log10cfu <- log10(bact + 1)

    Cc ~ prop(propSd)
    log10cfu ~ add(addSd_log10cfu)
  })
}
