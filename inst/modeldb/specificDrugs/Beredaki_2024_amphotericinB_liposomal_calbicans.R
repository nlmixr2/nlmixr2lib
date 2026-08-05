Beredaki_2024_amphotericinB_liposomal_calbicans <- function() {
  description <- "In vitro (two Candida albicans isolates; one-compartment dilution PK/PD model). Liposomal amphotericin B (L-AMB) in vitro PK/PD model for the Candida albicans arm that Beredaki 2024 used to VALIDATE the in vitro system against a published neutropenic-mouse model of disseminated candidiasis. Beredaki 2024 simulated MOUSE L-AMB pharmacokinetics in the internal compartment of a dilution model, re-spiking to a constant target peak once daily with an average half-life of 8 h (range 7-11 h; target 11 h). The 48 h change in log10 CFU/mL from the starting inoculum was then related to the PK/PD index Cmax/MIC with the sigmoidal variable-slope Emax model E = Emax * EI^n / (EI50^n + EI^n) (Methods, 'PK/PD analysis'), where E is the REDUCTION in log10 CFU/mL relative to the drug-free control and the MIC is the CLSI M27 amphotericin B deoxycholate MIC. Beredaki 2024 printed the equation and its goodness of fit (R^2 = 0.91) but no coefficients, so e0, lemax, lec50 and lhill were recovered by digitising the fitted curve in Figure 2; the recovered coefficients place the stasis exposure at 2.06 Cmax/MIC against the 2.1 (0.5-3.9) printed in Figure 2, which in turn brackets the 1.6-3.8 Cmax/MIC required for stasis in mouse kidneys. The fungal density bact (linear CFU/mL) is integrated as d/dt(bact) = ln(10) * ((e0 - kill48) / 48) * bact so that log10(bact) changes by exactly (e0 - kill48) across the paper's 48 h observation window, reproducing the endpoint model exactly at the time the paper actually fitted it. The Candida auris models this arm validates are Beredaki_2024_amphotericinB_liposomal_cauris48h and Beredaki_2024_amphotericinB_liposomal_cauris24h."
  reference <- paste(
    "Beredaki MI, Sanidopoulos I, Pournaras S, Meletiadis J. (2024).",
    "Defining optimal doses of liposomal amphotericin B against Candida auris:",
    "data from an in vitro pharmacokinetic/pharmacodynamic model.",
    "The Journal of Infectious Diseases 229(2):599-607.",
    "doi:10.1093/infdis/jiad583.",
    sep = " "
  )
  vignette <- "Beredaki_2024_amphotericinB_liposomal"
  units <- list(
    time = "h",
    dosing = "mg (L-AMB added to the internal compartment of the dilution model)",
    concentration = "mg/L (L-AMB, Cc); log10 CFU/mL (fungal density, log10cfu)"
  )

  # bact -- linear fungal density (CFU/mL) in the internal compartment.
  # Follows the in-vitro PD accumulator precedent in Chen_2023_tilmicosin.
  paper_specific_compartments <- c("bact")

  # specimen is "not applicable" for both states: the matrix is the RPMI-1640
  # growth medium in the internal compartment of the in vitro dilution model,
  # which has no entry in the biological-specimen vocabulary. The matrix is
  # recorded in the analyte strings and in population$system.
  compartmentData <- list(
    central = list(analyte = "liposomal amphotericin B in the internal compartment (RPMI-1640 growth medium)", units = "mg", specimen = "not applicable", verified = TRUE),
    bact    = list(analyte = "Candida albicans viable count in the internal compartment (RPMI-1640 growth medium)", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "in vitro (two Candida albicans isolates; one-compartment dilution PK/PD model)",
    n_subjects     = 2L,
    n_studies      = 1L,
    organism       = paste(
      "Two Candida albicans isolates (Beredaki 2024 Table 1, median MIC in",
      "mg/L, CLSI M27-A3, amphotericin B deoxycholate / liposomal",
      "amphotericin B): C. albicans K1, wild type, 0.25 / 0.25, kindly",
      "provided by David Andes (University of Wisconsin) and previously",
      "tested in a neutropenic mouse model of disseminated candidiasis;",
      "C. albicans SSI-2699, fks1 (S649P) and Erg2 (F105fs) amphotericin",
      "B-resistant, >16 / 16, kindly provided by Maiken C. Arendrup",
      "(University of Copenhagen). The exposure index is formed on the CLSI",
      "amphotericin B deoxycholate MIC (Figure 2 x axis",
      "'Cmax L-AMB / CLSI MIC'; the Figure 2 legend quotes MIC 0.25 mg/L for",
      "K1 and >16 mg/L for SSI-2699, i.e. the deoxycholate values of",
      "Table 1)."
    ),
    system         = paste(
      "Previously validated one-compartment in vitro PK/PD dilution model:",
      "a 250 mL conical glass culture vessel (internal compartment) holding an",
      "initial 5 mL of fresh RPMI-1640, connected to a peristaltic pump",
      "(Minipuls Evolution, Gilson) that adds fresh medium at a rate matching",
      "the clearance of L-AMB in human plasma. The vessel is covered with",
      "aluminium foil to minimise light exposure and held at 37 degrees C on a",
      "magnetic stirrer; its volume rises to approximately 100 mL by 48 h.",
      "Growth medium was RPMI 1640 with L-glutamine, without bicarbonate,",
      "buffered to pH 7.0 with 0.165 M MOPS and supplemented with 100 mg/L",
      "chloramphenicol. L-AMB was AmBisome (Gilead Sciences) reconstituted to",
      "4 mg/mL."
    ),
    disease_state  = "In vitro infection of the internal compartment at a starting inoculum of 10^4 CFU/mL",
    dose_range     = paste(
      "L-AMB target peak concentrations Cmax 0.125-128 mg/L, added once",
      "daily. Target average half-life 11 h, chosen to reproduce the L-AMB",
      "exposures of the neutropenic mouse model of disseminated candidiasis",
      "(intraperitoneal L-AMB 0.312-80 mg/kg once daily for 72 h) that this",
      "arm validates against; achieved in vitro half-life 8 h (range 7-11 h)",
      "with peak concentrations within 20% of target."
    ),
    sampling       = paste(
      "PK: samples from the internal compartment assayed by microbiological",
      "diffusion against a Paecilomyces variotii strain; lowest limit of",
      "detection 0.03 mg/L for the 80% partial growth inhibition zone, linear",
      "over 0.03-1 mg/L (Figure 3A, r^2 = 0.9373), so samples expected above",
      "1 mg/L were pre-diluted into the linear range.",
      "PD: 500 uL sampled at regular intervals up to 48 h, 10-fold serially",
      "diluted in normal saline, 20 uL subcultured on Sabouraud dextrose agar,",
      "incubated at 30 degrees C for 24 h; dilutions yielding 10-50 colonies",
      "were counted (Figure 1 shows 0, 3, 6, 24, 27 and 48 h)."
    ),
    regions        = "Greece (Clinical Microbiology Laboratory, Attikon University Hospital, Athens)",
    notes          = paste(
      "Two independent experiments were conducted. This arm exists to",
      "validate the in vitro dilution model: the C. albicans K1 isolate had",
      "previously been used in a published neutropenic mouse model of",
      "disseminated candidiasis, and the in vitro fungistatic target",
      "recovered here (2.1 Cmax/MIC) is compared against the 1.6-3.8",
      "Cmax/MIC required for stasis in mouse kidneys once the L-AMB",
      "concentration at the site of infection is accounted for (kidney AUC",
      "for stasis 10 mg*h/L at MIC 0.25, i.e. 40 AUC/MIC, divided by the",
      "10.55-24.37 kidney AUC/Cmax ratio). A reduction of more than 1 log10",
      "CFU/mL was observed at 48 h for Cmax at or above 8 mg/L, while at the",
      "lowest doses regrowth appeared as early as 6 h. No killing was seen",
      "for the resistant isolate SSI-2699 at any exposure except Cmax",
      "128 mg/L. Data were analysed in GraphPad Prism 5.0. Note that the",
      "y axis of Figure 1 is labelled 'Change in log10 CFU/mL' but plots",
      "absolute log10 CFU/mL (curves start near the 10^4 CFU/mL inoculum);",
      "Figure 2 plots the true change from the initial inoculum."
    )
  )

  ini({
    # ================================================================
    # In vitro pharmacokinetics
    # Beredaki 2024 Methods, "Validation of the in vitro PK/PD model":
    # "L-AMB exposures in mice with Cmax 0.125-128 mg/L and average
    # half-life of 11 hours were simulated in the in vitro model", added
    # once daily. Figure 3B (the C. auris panel of the same apparatus)
    # shows the target profile returning to the SAME peak at t = 24 h,
    # i.e. the compartment is re-spiked to a constant peak rather than
    # dosed additively.
    # ================================================================
    # Beredaki 2024 Results, "Validation of the in vitro PK/PD model with
    # C. albicans pharmacokinetics": "The calculated Cmax in the IC was
    # within 20% of the target Cmax of 0.125-128 mg/L with mean half-life
    # (t1/2) of 8 (range, 7-11) hours". The ACHIEVED half-life is used,
    # matching the system that generated the pharmacodynamic data; the
    # 11 h design target is reproduced in the vignette by overriding lkel.
    lkel <- log(log(2) / 8)
    label("Log first-order rate of L-AMB removal from the internal compartment (1/h)")  # Beredaki 2024 Results: achieved in vitro t1/2 = 8 h (range 7-11 h); Methods design target 11 h (mouse L-AMB)

    # Volume of the internal compartment at the time the peak occurs.
    # FIXED: a measured dimension of the apparatus, not an estimated
    # parameter. It converts the L-AMB amount added to the culture vessel
    # into a concentration. Beredaki 2024 reports the compartment volume
    # rising to about 100 mL by 48 h, and abstracts the resulting dilution
    # into the reported mono-exponential half-life above, so a constant
    # volume anchored at the initial 5 mL reproduces the paper's own
    # concentration-time description exactly. See the vignette
    # "Assumptions and deviations".
    lvc <- fixed(log(0.005))
    label("Log volume of the internal compartment at dosing (L; apparatus dimension)")  # Beredaki 2024 Methods, "In vitro PK/PD model": culture vessel "containing fresh RPMI-1640 medium to an initial volume of 5 mL" = 0.005 L

    # Target peak L-AMB concentration reached at the start of every 24 h
    # dosing interval, and the numerator of the paper's PK/PD index.
    # FIXED experimental design input; change it to simulate another arm.
    # The default is the lowest arm at which Beredaki 2024 observed a
    # reduction of more than 1 log10 CFU/mL for C. albicans K1.
    cmax <- fixed(8)
    label("Target peak L-AMB concentration per 24 h interval (mg/L; experimental design input)")  # Beredaki 2024 Results, "Pharmacodynamics": "A reduction (>1 log10 CFU/mL) ... was observed at 48 hours for the highest L-AMB doses with Cmax >= 8 mg/L"; the simulated arms spanned Cmax 0.125-128 mg/L

    # ================================================================
    # Sigmoidal variable-slope Emax exposure-response, 48 h endpoint
    # Beredaki 2024 Methods, "PK/PD analysis":
    #   E = Emax * EI^n / (EI50^n + EI^n)
    # where EI is the exposure index Cmax/MIC, Emax is the maximum growth
    # rate, EI50 is the index giving 50% of Emax and n is the Hill
    # coefficient. Read with E as the REDUCTION relative to the drug-free
    # control, the quantity plotted in Figure 2 is (e0 - E) -- the printed
    # form alone gives E = 0 at zero exposure whereas Figure 2 starts near
    # +2.3 log10 CFU/mL.
    #
    # PROVENANCE: Beredaki 2024 printed this equation and the fitted
    # curve (Figure 2, R^2 = 0.91) but no coefficients. The four values
    # below were recovered by digitising the Figure 2 curve from the
    # 600 dpi page render (1459 traced points spanning Cmax/MIC
    # 0.0011-66) and least-squares fitting the four-parameter
    # log-logistic form; residual RMSE was 0.010 log10 CFU/mL, i.e. the
    # recovered coefficients reproduce the drawn curve to within its line
    # width. See the vignette "Source trace" and "Assumptions and
    # deviations".
    # ================================================================
    # Top asymptote: the 48 h change in log10 CFU/mL at zero exposure.
    # Kept on the natural scale because it is a signed net change. This
    # is the fitted asymptote; the drug-free control of Figure 1 rises
    # from about 4.0 to about 6.0 log10 CFU/mL over 48 h, so the fitted
    # +2.28 is consistent with the plotted control growth.
    e0 <- 2.283
    label("48 h change in log10 CFU/mL at zero exposure (fitted top asymptote)")  # digitised from Beredaki 2024 Figure 2
    lemax <- log(5.215)
    label("Log maximum reduction Emax relative to the drug-free control (log10 CFU/mL)")  # digitised from Beredaki 2024 Figure 2; implies a bottom asymptote of 2.283 - 5.215 = -2.93 log10 CFU/mL
    lec50 <- log(2.695)
    label("Log Cmax/MIC producing 50% of Emax (EI50, unitless)")  # digitised from Beredaki 2024 Figure 2
    lhill <- log(0.9306)
    label("Log Hill coefficient n of the exposure-response relationship (unitless)")  # digitised from Beredaki 2024 Figure 2

    # ================================================================
    # Minimum inhibitory concentration
    # ================================================================
    # The PK/PD index driving the sigmoid is Cmax/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of the
    # isolate, not an estimated parameter. Beredaki 2024 indexes on the
    # CLSI M27 amphotericin B DEOXYCHOLATE MIC, not the L-AMB MIC.
    # Default is the wild-type reference isolate C. albicans K1; change it
    # to apply the model to another isolate or to an arbitrary MIC from
    # the CLSI distribution. The resistant isolate SSI-2699 is reported
    # only as ">16" mg/L, so it has no usable point value.
    mic <- fixed(0.25)
    label("Amphotericin B CLSI M27 MIC of the simulated isolate (mg/L; measured property of the isolate)")  # Beredaki 2024 Table 1: C. albicans K1 (wild type) CLSI amphotericin B deoxycholate median MIC 0.25 mg/L (range 0.25)

    # ================================================================
    # Starting fungal density
    # ================================================================
    # FIXED experimental design input. Beredaki 2024 quotes a measured
    # mean t = 0 density only for the C. auris arm, so the C. albicans
    # model starts from the nominal inoculum, which the paper states was
    # confirmed by quantitative culture. Figure 1 shows the C. albicans
    # curves starting at about 4.0-4.2 log10 CFU/mL, consistent with it.
    log10_cfu0 <- fixed(4)
    label("Log10 starting fungal density in the internal compartment (log10 CFU/mL)")  # Beredaki 2024 Methods, "Isolates": inoculum suspension "adjusted to a final inoculum of 10^4 colony-forming units (CFU)/mL. The CFU number was confirmed by quantitative cultures on SDA plates" = 4 log10 CFU/mL

    # ================================================================
    # Residual error
    # ================================================================
    # L-AMB concentration: Beredaki 2024 reports only that the achieved
    # peaks were "within 20% of the target Cmax" for this arm (Results,
    # "Validation of the in vitro PK/PD model with C. albicans
    # pharmacokinetics") and gives no assay CV point estimate, so the
    # proportional residual SD is FIXED at that reported bound.
    propSd <- fixed(0.20)
    label("Proportional residual error on L-AMB concentration (fraction; the reported bound on peak accuracy)")  # Beredaki 2024 Results: "The calculated Cmax in the IC was within 20% of the target Cmax of 0.125-128 mg/L"

    # Fungal density: Beredaki 2024 reported only the coefficient of
    # determination of the Emax fit (R^2 = 0.91) and no residual standard
    # deviation, so the log10-density residual SD is held at zero for
    # deterministic typical-value simulation. See the vignette
    # "Assumptions and deviations".
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Beredaki 2024)")  # Beredaki 2024 Figure 2 reports R^2 = 0.91 only, no residual SD
  })

  model({
    kel <- exp(lkel)
    vc <- exp(lvc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---------------------------------------------------------------
    # In vitro pharmacokinetics (Beredaki 2024 Methods, "Pharmacokinetic
    # analysis"). L-AMB is added to the internal compartment once daily
    # and decays mono-exponentially. Dose the `central` compartment in
    # mg; the vignette computes the amount that re-spikes the compartment
    # to `cmax` at the start of every interval, matching the constant
    # target peaks the apparatus produced (Figure 3B).
    # ---------------------------------------------------------------
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # ---------------------------------------------------------------
    # PK/PD index and exposure-response (Beredaki 2024 Methods, "PK/PD
    # analysis"). Cmax is a design input rather than a state because the
    # paper fitted a 48 h ENDPOINT against the per-interval target peak,
    # so the index must be available from t = 0 for the trajectory below
    # to land on the fitted endpoint. kill48 is the paper's E: the
    # reduction in log10 CFU/mL over the 48 h observation window relative
    # to the drug-free control, rising from 0 at zero exposure towards
    # emax at saturating exposure.
    # ---------------------------------------------------------------
    ei <- cmax / mic
    kill48 <- emax * ei^hill / (ei^hill + ec50^hill)

    # Beredaki 2024 fitted only the 48 h change in log10 CFU/mL.
    # Spreading that change uniformly over the window makes log10(bact)
    # move by exactly (e0 - kill48) between t = 0 and t = 48 h, so the
    # trajectory matches the fitted endpoint model at the time the paper
    # fitted it:
    #   d(log10 N)/dt = (e0 - kill48) / 48
    #   => d(N)/dt    = ln(10) * ((e0 - kill48) / 48) * N
    d/dt(bact) <- log(10) * ((e0 - kill48) / 48) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1-CFU/mL floor so the log stays finite if bact
    # is driven below 1 CFU/mL (matches the in-vitro PD convention used
    # by Chen_2023_tilmicosin).
    log10cfu <- log10(bact + 1)

    Cc ~ prop(propSd)
    log10cfu ~ add(addSd_log10cfu)
  })
}
