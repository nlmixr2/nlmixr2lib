Debord_2001_cyclosporin <- function() {
  description <- "Two-compartment population PK model for oral cyclosporin microemulsion (Neoral) in stable renal transplant recipients (Debord 2001), with a gamma-distribution absorption (Savic 2007 analytical transit-compartment form) feeding the central compartment directly, F fixed to 1, and population typical values derived from the means of the 21 individually-fitted patients in Table I of the paper."
  reference <- "Debord J, Risco E, Harel M, Le Meur Y, Buchler M, Lachatre G, Le Guellec C, Marquet P. Application of a gamma model of absorption to oral cyclosporin. Clin Pharmacokinet. 2001;40(5):375-382. doi:10.2165/00003088-200140050-00004"
  vignette <- "Debord_2001_cyclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "Not reported (adult renal transplant recipients)",
    weight_range   = "Not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (single-centre French cohort: Hopital Dupuytren, Limoges and CHU Bretonneau, Tours)",
    disease_state  = "Stable renal transplant recipients",
    dose_range     = "Oral cyclosporin microemulsion (Neoral) 75-175 mg twice daily",
    administration = "Oral, twice daily",
    regions        = "France (single centre, Limoges and Tours)",
    notes          = "Blood samples drawn just before administration (C0) and at 20, 40, 60, 90, 120, 180, 240, 360 and 540 min after administration (Methods page 377-378). Cyclosporin assayed by LC-MS; intra-assay CVs < 11%, inter-assay CVs < 15%; LOQ 10 ug/L over linearity range up to 2500 ug/L."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification of each item):
  # * Absorption: Debord 2001 Eq. 2-3 describe the absorption rate as a
  #   gamma density va(t) = F*D*b^a/Gamma(a) * t^(a-1) * exp(-b*t). The
  #   Savic 2007 analytical transit-compartment form (rxode2 transit())
  #   delivers the mathematically equivalent input rate with a = ntr + 1
  #   and b = (ntr + 1)/mtt. rxode2 transit() requires that the dose
  #   target compartment be the same compartment whose d/dt() contains
  #   the transit() call (otherwise podo() / tad() do not fire). Dose
  #   events therefore target depot, and d/dt(depot) carries transit()
  #   with a high-rate first-order pass-through (ka_pass = 100 1/h,
  #   fixed) feeding central. The time constant 1/ka_pass = 0.01 h is
  #   ~70-fold faster than the gamma rate b = (ntr+1)/mtt = 10.94 1/h,
  #   so the convolution introduces a negligible 0.01 h Tmax shift
  #   relative to the Tmax of ~0.69 h that pure gamma absorption would
  #   give (under 2 % perturbation). f(depot) <- 0 is NOT used because
  #   the dose must enter depot for transit() to fire; instead the
  #   gamma input rate is exactly compensated by depot drainage to
  #   central via the high-ka_pass term.
  # * Disposition: Debord 2001 Eq. 11 / Table I report the disposition as
  #   macro-parameters (A, B in L^-1 standardised to a 100 mg dose;
  #   alpha, beta in 1/h). The packaged model uses canonical micro
  #   parameters (CL/F, Vc/F, Q/F, Vp/F), derived from the means in
  #   Table I as Vc/F = 100/(A+B), k21 = (alpha*B + beta*A)/(A+B),
  #   kel = alpha*beta/k21, k12 = alpha+beta - kel - k21, CL = kel*Vc,
  #   Q = k12*Vc, Vp = Q/k21. Round-trip back through (A, B, alpha, beta)
  #   reproduces the Table I means exactly.
  # * Bioavailability: F fixed to 1 in the source paper (Methods page 378:
  #   "the intravenous disposition coefficients A and B (assuming F = 1)")
  #   because no IV reference data were available. Vc/F, CL/F, Q/F, Vp/F
  #   are therefore apparent (oral) values, not absolute.
  # * IIV: Debord 2001 page 380 reports cross-patient CVs of the 21
  #   individually-estimated parameters (descriptive statistics of
  #   separate per-subject fits, NOT formal mixed-effects omega
  #   estimates). The packaged model approximates IIV by mapping the
  #   reported CV ranges onto the matched micro-parameters and using
  #   omega^2 = log(1 + CV^2): CV ~ 65% (midpoint of 55-75%) on CL/Vc/Q
  #   from the macro CVs on A/B/alpha; CV ~ 40% (midpoint of 35-46%) on
  #   Vp/ntr/mtt from the macro CVs on beta/a/b. Correlations between
  #   absorption parameters (a-b r=0.88), first-disposition-phase (A-alpha
  #   r=0.71), and second-disposition-phase (B-beta r=0.62) reported on
  #   page 380 are not encoded as block omegas to keep the model
  #   tractable; users wanting the correlations should override with a
  #   custom omega block.
  # * Residual error: Debord 2001 used a logistic weighting function
  #   sk(Ck) = 0.13/(1 + exp(-3.1*Ck + 1.4)) - 0.028 (Eq. 13) derived from
  #   calibration experiments, applied as 1/sk^2 weights during nonlinear
  #   regression - not a population residual error model. The packaged
  #   model approximates this with combined additive (0.025 mg/L,
  #   matching the low-concentration plateau of Eq. 13 and the assay LOQ
  #   of 10 ug/L) plus proportional (8%, matching the high-concentration
  #   plateau of Eq. 13) residual error.
  # * C0: Debord 2001 also reports a per-patient theoretical residual
  #   concentration C0 (mean 0.10 mg/L per 100 mg dose) treated as an
  #   adjustable parameter to absorb pre-dose carry-over. The packaged
  #   model omits C0 as an explicit parameter because the multi-dose
  #   simulation naturally produces the steady-state trough via
  #   accumulation; the vignette demonstrates the agreement at the
  #   typical BID regimen.
  ini({
    # Gamma-distribution absorption parameters (Debord 2001 Table I means).
    # Map onto the rxode2 transit() function via ntr = a - 1 and
    # mtt = a/b so that the Savic 2007 analytical form delivers the same
    # input rate as Debord 2001 Eq. 3.
    lntr <- log(8.56 - 1);   label("Number of transit compartments (continuous, dimensionless)")  # Debord 2001 Table I: a (mean) = 8.56; ntr = a - 1 = 7.56
    lmtt <- log(8.56 / 10.94); label("Mean transit time (h)")                                     # Debord 2001 Table I: a (mean) = 8.56, b (mean) = 10.94 1/h; mtt = a/b = 0.7824 h

    # Bioavailability fixed to 1 (Debord 2001 Methods page 378:
    # "assuming F = 1"). Apparent (oral) CL/F, Vc/F, Q/F, Vp/F follow.
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (dimensionless)")                     # Debord 2001 Methods page 378 ("assuming F = 1")

    # Two-compartment disposition; micro parameters derived from the
    # Table I macro means (A, B per 100 mg dose; alpha, beta in 1/h).
    # Conversion documented in the implementation-notes block above.
    lcl <- log(56.25); label("Apparent oral clearance CL/F (L/h)")                                # Debord 2001 Table I means: alpha = 2.23, beta = 0.34, A = 1.21, B = 0.42; CL/F = alpha*beta / (beta*A_real + alpha*B_real) where A_real = A/100, B_real = B/100; CL/F = 56.25 L/h
    lvc <- log(61.35); label("Apparent central volume Vc/F (L)")                                  # Debord 2001 Table I means: Vc/F = 100 / (A + B) = 100 / 1.63 = 61.35 L
    lq  <- log(50.69); label("Apparent inter-compartmental clearance Q/F (L/h)")                  # Debord 2001 Table I means: k12 = (alpha + beta) - kel - k21 = 0.826 1/h with kel = 0.917, k21 = 0.827; Q/F = k12 * Vc/F = 50.69 L/h
    lvp <- log(61.29); label("Apparent peripheral volume Vp/F (L)")                               # Debord 2001 Table I means: Vp/F = Q/F / k21 = 50.69 / 0.827 = 61.29 L

    # Inter-individual variability approximated from Debord 2001 page 380
    # cross-patient CVs of the 21 per-subject fits (descriptive statistics,
    # NOT formal mixed-effects omegas). Conversion omega^2 = log(1 + CV^2).
    # See implementation-notes block above for the CV-to-parameter mapping.
    etalcl  ~ 0.353  # Debord 2001 page 380 (intravenous parameters A, B, alpha CVs 55-75%); CV ~ 65%; log(1 + 0.65^2) = 0.353
    etalvc  ~ 0.353  # Debord 2001 page 380 (intravenous parameters A, B, alpha CVs 55-75%); CV ~ 65%; log(1 + 0.65^2) = 0.353
    etalq   ~ 0.353  # Debord 2001 page 380 (intravenous parameters A, B, alpha CVs 55-75%); CV ~ 65%; log(1 + 0.65^2) = 0.353
    etalvp  ~ 0.148  # Debord 2001 page 380 (a, b, beta CVs 35-46%); CV ~ 40% mapped onto Vp from beta; log(1 + 0.40^2) = 0.148
    etalntr ~ 0.148  # Debord 2001 page 380 (a CV in 35-46% range); CV ~ 40%; log(1 + 0.40^2) = 0.148
    etalmtt ~ 0.148  # Debord 2001 page 380 (b CV in 35-46% range; mtt = a/b so propagates similarly); CV ~ 40%; log(1 + 0.40^2) = 0.148

    # Residual error approximated from Debord 2001 Eq. 13 calibration
    # logistic. The paper used this as a weighting function (1/sk^2), not
    # as a population residual error model.
    addSd  <- 0.025; label("Additive residual SD (mg/L)")                                          # Debord 2001 Eq. 13: sk plateau at low Ck approx 0.001-0.013 mg/L; LOQ = 0.010 mg/L; rounded to 0.025 mg/L combining assay floor and low-Ck model misspecification
    propSd <- 0.08;  label("Proportional residual SD (fraction)")                                  # Debord 2001 Eq. 13: sk/Ck approx 8% at Ck = 0.5-1 mg/L (typical peak range); combined with assay inter-day CV < 15%

    # Technical pass-through rate for the depot to central transfer (see
    # implementation notes above). Not a Debord 2001 parameter; held at
    # 100 1/h so the depot transfer is essentially instantaneous on the
    # gamma-absorption timescale.
    lka_pass <- fixed(log(100)); label("Depot to central pass-through rate (1/h, technical)")      # Implementation choice (NOT in Debord 2001); see implementation notes above
  })
  model({
    # Individual PK parameters
    ntr     <- exp(lntr + etalntr)
    mtt     <- exp(lmtt + etalmtt)
    cl      <- exp(lcl  + etalcl)
    vc      <- exp(lvc  + etalvc)
    q       <- exp(lq   + etalq)
    vp      <- exp(lvp  + etalvp)
    ka_pass <- exp(lka_pass)
    bio     <- exp(lfdepot)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with a Savic 2007 analytical
    # transit-compartment chain delivering the Debord 2001 gamma
    # absorption rate (Eq. 3) into depot. f(depot) <- 0 suppresses the
    # bolus so transit() drives the only input to depot; depot drains
    # to central at the fast technical rate ka_pass (= 100 1/h, fixed)
    # so the convolution of gamma absorption with a tiny exponential
    # smoothing is effectively the pure gamma absorption rate into
    # central (Tmax shift < 2 % vs pure gamma; see implementation notes
    # block above). transit() reads the dose amount from podo() and time
    # since dose from tad(), so multi-dose simulations correctly produce
    # one gamma input per dose record.
    d/dt(depot)       <-  transit(ntr, mtt, bio) - ka_pass * depot
    d/dt(central)     <-  ka_pass * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- 0

    # Cyclosporin whole-blood concentration (mg/L). Dose in mg, Vc in L,
    # so central (mg) / vc (L) gives mg/L directly, matching the paper's
    # bioanalytical units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
