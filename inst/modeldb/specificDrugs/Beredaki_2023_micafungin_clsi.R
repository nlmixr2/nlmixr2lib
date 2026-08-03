Beredaki_2023_micafungin_clsi <- function() {
  description <- "In vitro (Candida albicans clinical isolates in RPMI-1640 with 10% pooled human serum; two-compartment closed dialysis/diffusion PK/PD model). Micafungin in vitro PK/PD model with the exposure-response relationship indexed on CLSI M27 MICs. Beredaki 2023 simulated q24h micafungin exposures in the internal compartment of a dialysis/diffusion model and described the concentration decay with the one-compartment mono-exponential Ct = C0 * exp(-k * t) (Methods, 'In vitro pharmacokinetics'), with an average half-life of 14 h (range 9-15 h) in the presence of 10% pooled human serum. The 72 h change in log10 CFU/mL relative to the starting inoculum was then related to the PK/PD index fAUC0-24/MIC with the sigmoidal variable-slope Emax model E = Emax * EI^n / (EI^n + EI50^n) (Methods, 'PK/PD analysis'), where E is the REDUCTION in log10 CFU/mL relative to the drug-free control. Beredaki 2023 printed the equation but not its coefficients, so e0, lemax, lec50 and lhill were recovered by digitising the fitted curve in Figure 5(a) (CLSI panel, R^2 = 0.92); the recovered coefficients reproduce the paper's own Table 2 CLSI stasis target (2.8 fAUC0-24/MIC) to +3% and its 1-log-kill target (9.2 fAUC0-24/MIC) to -13%, equivalent to 0.10 log10 CFU/mL on the effect axis. The fungal density bact (linear CFU/mL) is integrated as d/dt(bact) = ln(10) * ((e0 - kill72) / 72) * bact so that log10(bact) changes by exactly (e0 - kill72) across the paper's 72 h observation window, reproducing the endpoint model exactly at the time the paper actually fitted it. The companion EUCAST parameterisation is Beredaki_2023_micafungin_eucast. See the vignette for the reported no-serum targets, which are NOT packaged because Beredaki 2023 published no exposure-response curve for that condition."
  reference <- paste(
    "Beredaki MI, Arendrup MC, Andes D, Meletiadis J. (2023).",
    "Development of an in vitro pharmacokinetic/pharmacodynamic model in the",
    "presence of serum for studying micafungin activity against Candida albicans:",
    "a need for revision of CLSI susceptibility breakpoints.",
    "Journal of Antimicrobial Chemotherapy 78(6):1386-1394.",
    "doi:10.1093/jac/dkad096.",
    sep = " "
  )
  vignette <- "Beredaki_2023_micafungin"
  units <- list(
    time = "h",
    dosing = "mg (micafungin added to the 10 mL internal compartment)",
    concentration = "mg/L (micafungin, Cc); log10 CFU/mL (fungal density, log10cfu)"
  )

  # bact  -- linear fungal density (CFU/mL) in the internal compartment.
  # fauc_0_24 -- unbound micafungin AUC accumulated over the FIRST 24 h
  #   dosing interval only; the numerator of the paper's PK/PD index. Follows
  #   the windowed-AUC accumulator precedent in Kuchimanchi_2018_evolocumab_ldlc.
  paper_specific_compartments <- c("fauc_0_24", "bact")

  population <- list(
    species        = "in vitro (four Candida albicans clinical isolates; two-compartment closed dialysis/diffusion PK/PD model)",
    n_subjects     = 4L,
    n_studies      = 1L,
    organism       = paste(
      "Four Candida albicans isolates (Beredaki 2023 Table 1, median MIC in mg/L,",
      "EUCAST E.Def 7.3 / CLSI M27): CA 580 fks1 wild-type 0.016/0.008;",
      "CA 9817 fks1 wild-type 0.03/0.03; CA SSI-5318 fks1 F641L (weak resistance)",
      "0.03/0.06; CA SSI-6683 fks1 R647G (strong resistance) 0.5/0.5.",
      "CA 580 and CA 9817 had previously been used in a neutropenic murine",
      "candidiasis model. Quality control: Candida parapsilosis ATCC 22019 and",
      "Candida krusei ATCC 6258."
    ),
    system         = paste(
      "Previously described two-compartment closed diffusion/dialysis in vitro",
      "PK/PD model: an external compartment (conical flask on a peristaltic pump)",
      "connected to an internal compartment consisting of a 10 mL semipermeable",
      "cellulose dialysis tube (Spectra/Por Float-A-Lyzer G2, 300 kDa molecular",
      "weight cut-off) inoculated with the yeast suspension. Growth medium was",
      "RPMI 1640 with L-glutamine, without bicarbonate, buffered to pH 7.0 with",
      "0.165 M MOPS and supplemented with 100 mg/L chloramphenicol, containing",
      "10% heat-inactivated pooled human serum (56 degrees C for 45 min)."
    ),
    disease_state  = "In vitro infection of the internal compartment at a starting inoculum of 10^4 CFU/mL",
    dose_range     = paste(
      "Micafungin target peak concentrations Cmax 0.004-32 mg/L simulated q24h",
      "for 72 h; the serum arms shown in Figure 4 used total Cmax of 0.25, 1, 4,",
      "8, 16 and 32 mg/L plus a drug-free control. Target half-life 15 h",
      "(human micafungin); achieved in vitro half-life 14 h (range 9-15 h) with",
      "10% serum and 9 h (range 8-10 h) without serum."
    ),
    sampling       = paste(
      "PK: 100 uL from the internal compartment at repeated intervals",
      "(Figure 2 shows 0, 4, 8, 12 h within each 24 h interval) assayed by",
      "microbiological diffusion against Aspergillus fumigatus AZN8196;",
      "limit of detection 0.125 mg/L, interexperimental CV < 3%.",
      "PD: 200 uL for quantitative culture at regular intervals up to 72 h",
      "(Figure 4 shows 0, 8, 24, 48 and 72 h), plated on SGC2 and counted",
      "after 24 h at 30 degrees C."
    ),
    regions        = "Greece (Attikon University Hospital, Athens); isolates from Greece and Statens Serum Institut, Denmark",
    notes          = paste(
      "Two independent experiments were conducted; Beredaki 2023 Table 2 reports",
      "the mean (95% CI) PK/PD targets across them while Figure 5 shows a single",
      "pooled Emax fit, which is why the packaged coefficients reproduce Table 2",
      "to within roughly 0.1 log10 CFU/mL rather than exactly. A preceding static",
      "time-kill experiment found no pharmacodynamic difference between 10%, 50%",
      "and 100% pooled human serum (static effect at 120, 141 and 167 total",
      "Cmax/MIC respectively), so 10% serum was used throughout the dynamic model."
    )
  )

  ini({
    # ================================================================
    # In vitro pharmacokinetics
    # Beredaki 2023 Methods, "In vitro pharmacokinetics": the
    # time-concentration profiles were fitted by non-linear regression to
    # the one-compartment model Ct = C0 * exp(-k * t) (printed in the PDF
    # as "Ct = Coe-k/t", a typesetting error for exp(-k * t): the paper
    # immediately defines t1/2 = 0.693/k, which is only consistent with a
    # mono-exponential in -k * t).
    # ================================================================
    # Beredaki 2023 Results, "In vitro dynamic PK/PD model": "an average
    # t1/2 of 9 h (8-10 h) and 14 h (9-15 h) in the absence and presence
    # of 10% pooled human serum, respectively". This model is the 10%
    # serum condition, so k = 0.693/14 using the paper's own constant.
    lkel <- log(0.693 / 14)
    label("Log first-order rate of micafungin removal from the internal compartment (1/h)")  # Beredaki 2023 Results: in vitro t1/2 = 14 h (9-15 h) with 10% pooled human serum; Methods: t1/2 = 0.693/k

    # Volume of the internal compartment. FIXED: a measured dimension of
    # the apparatus, not an estimated parameter. Converts the micafungin
    # amount added to the dialysis tube into a concentration.
    lvc <- fixed(log(0.010))
    label("Log volume of the internal compartment (L;, apparatus dimension)")  # Beredaki 2023 Methods, "In vitro PK/PD model": "an internal compartment (IC) of a 10 mL-volume semipermeable cellulose dialysis tube" = 0.010 L

    # Unbound fraction of micafungin. FIXED: a measured physicochemical
    # property. Beredaki 2023 applied it to convert the total in vitro
    # exposure into the fAUC0-24 that indexes the exposure-response.
    fu <- fixed(0.0025)
    label("Unbound fraction of micafungin in serum (measured)")  # Beredaki 2023 Methods, "PK/PD analysis": "For calculation of fAUC/MIC in serum, a protein binding of 99.75% was taken into account"; Discussion: "99.75% for micafungin"

    # Dosing interval. FIXED experimental design input. Also defines the
    # integration window of the paper's PK/PD index (fAUC0-24).
    tau <- fixed(24)
    label("Dosing interval and PK/PD index integration window (h;, design)")  # Beredaki 2023 Methods: "Drug concentrations were added at the corresponding Cmax values in the in vitro model once daily"; the index is fAUC0-24/MIC

    # Target peak total micafungin concentration reached at the start of
    # every dosing interval. FIXED experimental design input; change it to
    # simulate a different arm. Default 8 mg/L is one of the six serum-arm
    # levels. Beredaki 2023 re-spiked the internal compartment to this
    # peak once daily, so the peak is constant across intervals
    # (Figure 2a: target 20 mg/L at 0, 24 and 48 h).
    tcmax <- fixed(8)
    label("Target peak total micafungin concentration per 24 h interval (mg/L;, design)")  # Beredaki 2023 Figure 4 legend: total Cmax arms of 0.25, 1, 4, 8, 16 and 32 mg/L with 10% serum; Methods: Cmax range 0.004-32 mg/L

    # ================================================================
    # Sigmoidal variable-slope Emax exposure-response
    # Beredaki 2023 Methods, "PK/PD analysis":
    #   E = Emax * EI^n / (EI^n + EI50^n)
    # where EI is the PK/PD index fAUC0-24/MIC, Emax "is the maximum
    # reduction in log10 cfu/mL", EI50 is the index giving 50% of Emax and
    # n is the Hill coefficient. E is therefore the REDUCTION relative to
    # the drug-free control, and the quantity plotted in Figure 5 is the
    # change from the starting inoculum, (e0 - E).
    #
    # PROVENANCE: Beredaki 2023 printed this equation and the fitted
    # curve (Figure 5a, R^2 = 0.92) but no coefficients. The four values
    # below were recovered by digitising the Figure 5(a) curve from the
    # publisher figure at 600 dpi (839 points, x = 0.108-65.4) and
    # least-squares fitting the four-parameter log-logistic form; the
    # residual RMSE was 0.005 log10 CFU/mL, i.e. the recovered
    # coefficients reproduce the drawn curve to within line width. See
    # the vignette "Source trace" and "Assumptions and deviations".
    # ================================================================
    # Top asymptote of the fitted curve: the 72 h change in log10 CFU/mL
    # at zero exposure. Kept on the natural scale because it is a signed
    # net change. Note this is the fitted asymptote, which lies to the
    # left of the observed drug-free controls and so exceeds the mean
    # measured control growth of 3.39 log10 CFU/mL (4.37 -> 7.76).
    e0 <- 4.10
    label("72 h change in log10 CFU/mL at zero exposure (fitted top asymptote)")  # digitised from Beredaki 2023 Figure 5(a); cross-checks against the reported drug-free control growth 4.37 -> 7.76 log10 CFU/mL at 72 h
    lemax <- log(6.61)
    label("Log maximum reduction Emax relative to the drug-free control (log10 CFU/mL)")  # digitised from Beredaki 2023 Figure 5(a); implies a bottom asymptote of 4.10 - 6.61 = -2.51 log10 CFU/mL
    lec50 <- log(1.45)
    label("Log fAUC0-24/MIC producing 50% of Emax (EI50, unitless)")  # digitised from Beredaki 2023 Figure 5(a)
    lhill <- log(0.710)
    label("Log Hill coefficient n of the exposure-response relationship (unitless)")  # digitised from Beredaki 2023 Figure 5(a)

    # ================================================================
    # Minimum inhibitory concentration
    # ================================================================
    # The PK/PD index driving the sigmoid is fAUC0-24/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of the
    # isolate, not an estimated parameter. Default is the CLSI M27 median
    # MIC of the fks1 wild-type reference isolate CA 580. Change it to
    # apply the model to another isolate; this file is indexed on CLSI
    # MICs, so use Beredaki_2023_micafungin_eucast for EUCAST MICs.
    mic <- fixed(0.008)
    label("Micafungin CLSI M27 MIC of the simulated isolate (mg/L;, measured)")  # Beredaki 2023 Table 1: C. albicans CA 580 (fks1 wild-type) CLSI median MIC 0.008 mg/L (range 0.008-0.015)

    # ================================================================
    # Starting fungal density
    # ================================================================
    # FIXED experimental design input. The nominal inoculum was
    # 10^4 CFU/mL; the measured starting density in the presence of serum
    # is used so the simulated trajectory starts where Figure 4 does.
    log10_cfu0 <- fixed(4.37)
    label("Log10 starting fungal density in the internal compartment (log10 CFU/mL)")  # Beredaki 2023 Results: "In presence of serum, C. albicans grew from a mean +/- SD of 4.37 +/- 0.24 log10CFU/mL at t = 0 h"; nominal inoculum 10^4 CFU/mL

    # ================================================================
    # Residual error
    # ================================================================
    # Micafungin concentration: Beredaki 2023 Methods, "In vitro
    # pharmacokinetics" reports "interexperimental CV <3% for all drug
    # concentrations tested" for the microbiological diffusion assay.
    # FIXED at that reported upper bound; the paper gives no point
    # estimate of an assay CV.
    propSd <- fixed(0.03)
    label("Proportional residual error on micafungin concentration (fraction; the reported assay CV upper bound)")  # Beredaki 2023 Methods: LOD 0.125 mg/L with interexperimental CV < 3%

    # Fungal density: Beredaki 2023 reported only the coefficient of
    # determination of the Emax fit (R^2 = 0.92, Figure 5a) and no
    # residual standard deviation, so the log10-density residual SD is
    # held at zero for deterministic typical-value simulation. See the
    # vignette "Assumptions and deviations".
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Beredaki 2023)")  # Beredaki 2023 reported R^2 = 0.92 only, no residual SD
  })

  model({
    kel <- exp(lkel)
    vc <- exp(lvc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---------------------------------------------------------------
    # In vitro pharmacokinetics (Beredaki 2023 Methods, "In vitro
    # pharmacokinetics"). Micafungin is added to the internal compartment
    # once daily and decays mono-exponentially. Dose the `central`
    # compartment in mg; the vignette computes the amount that re-spikes
    # the compartment to `tcmax` at the start of every interval, matching
    # the constant peaks in Figure 2.
    # ---------------------------------------------------------------
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Unbound micafungin AUC accumulated over the FIRST dosing interval
    # only -- the numerator of the paper's PK/PD index fAUC0-24/MIC. The
    # t < tau gate freezes the state at fAUC0-24 so it can be compared
    # against the closed form used below (see the vignette's
    # self-consistency check). Windowed-accumulator idiom follows
    # Kuchimanchi_2018_evolocumab_ldlc.
    d/dt(fauc_0_24) <- fu * Cc * (t < tau)

    # ---------------------------------------------------------------
    # PK/PD index. Closed form of the same integral for a mono-exponential
    # interval starting at tcmax:
    #   fAUC0-24 = fu * integral(0..tau) tcmax * exp(-kel * t) dt
    #            = fu * tcmax * (1 - exp(-kel * tau)) / kel
    # Used rather than the integrated state because Beredaki 2023 fitted a
    # 72 h ENDPOINT against the completed fAUC0-24, so the index must be
    # available from t = 0 for the trajectory below to land on the fitted
    # endpoint. It equals fauc_0_24 at t = tau by construction.
    # ---------------------------------------------------------------
    fauc24 <- fu * tcmax * (1 - exp(-kel * tau)) / kel
    ei <- fauc24 / mic

    # ---------------------------------------------------------------
    # Exposure-response (Beredaki 2023 Methods, "PK/PD analysis").
    # kill72 is the paper's E: the reduction in log10 CFU/mL over the 72 h
    # observation window relative to the drug-free control, rising from 0
    # at zero exposure towards emax at saturating exposure.
    # ---------------------------------------------------------------
    kill72 <- emax * ei^hill / (ei^hill + ec50^hill)

    # Beredaki 2023 fitted only the 72 h change in log10 CFU/mL. Spreading
    # that change uniformly over the window makes log10(bact) move by
    # exactly (e0 - kill72) between t = 0 and t = 72 h, so the trajectory
    # matches the fitted endpoint model at the time the paper fitted it:
    #   d(log10 N)/dt = (e0 - kill72) / 72
    #   => d(N)/dt     = ln(10) * ((e0 - kill72) / 72) * N
    d/dt(bact) <- log(10) * ((e0 - kill72) / 72) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1-CFU/mL floor so the log stays finite if bact
    # is driven below 1 CFU/mL (matches the in-vitro PD convention used by
    # Chen_2023_tilmicosin).
    log10cfu <- log10(bact + 1)

    Cc ~ prop(propSd)
    log10cfu ~ add(addSd_log10cfu)
  })
}
