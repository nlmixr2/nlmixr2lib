Shao_2023_bedaquiline <- function() {
  description <- paste0(
    "Three-compartment population PK model for bedaquiline in adult Chinese ",
    "patients with multidrug-resistant tuberculosis (MDR-TB), with the dual ",
    "zero-order absorption input of McLeay 2014 that reproduces the double ",
    "absorption peak: a fixed fraction FR1 = 58.5% of the dose enters a depot ",
    "compartment as a zero-order infusion of duration DUR1 = 2.22 h beginning ",
    "after a lag of 1.00 h, and is then handed to the central compartment by a ",
    "near-instantaneous first-order step (Ka fixed at 1,000 1/h), while the ",
    "remaining 41.5% enters the central compartment directly as a zero-order ",
    "infusion of duration DUR2 = 1.48 h beginning after an additional lag of ",
    "2.76 h (total 3.76 h). Apparent clearance CL/F is scaled by body weight ",
    "with the allometric exponent fixed at 0.75 and by a steep power function ",
    "of serum albumin with an estimated exponent of 3.76; apparent central ",
    "volume Vc/F is scaled by body weight with the exponent fixed at 1. Both ",
    "covariates are normalised to the development-cohort medians of 67 kg and ",
    "37.1 g/L. Between-subject variability is carried on CL/F and Vc/F and the ",
    "residual error is proportional."
  )
  reference <- paste(
    "Shao G, Bao Z, Davies Forsman L, Paues J, Werngren J, Niward K,",
    "Schon T, Bruchfeld J, Alffenaar J-W, Hu Y (2023).",
    "Population pharmacokinetics and model-based dosing evaluation of",
    "bedaquiline in multidrug-resistant tuberculosis patients.",
    "Frontiers in Pharmacology 14:1022090. doi:10.3389/fphar.2023.1022090.",
    "The absorption sub-model constants (Ka fixed at 1,000 1/h, DUR1 = 2.22 h,",
    "DUR2 = 1.48 h and FR1 = 58.5%) were fixed by Shao 2023 from",
    "McLeay SC, Vis P, van Heeswijk RPG, Green B (2014),",
    "Population pharmacokinetics of bedaquiline (TMC207), a novel",
    "antituberculosis drug, Antimicrobial Agents and Chemotherapy",
    "58(9):5315-5324, doi:10.1128/AAC.01418-13; every one of those values is",
    "reprinted in Shao 2023 Methods section 2.2 and Table 2, so no value in",
    "this file is taken from the McLeay paper itself.",
    sep = " "
  )
  vignette <- "Shao_2023_bedaquiline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Baseline body weight, normalised to the development-cohort median of ",
        "67 kg (Shao 2023 Table 1, development cohort: median 67 kg, IQR ",
        "59-74 kg). Enters CL/F with the allometric exponent fixed at 0.75 and ",
        "Vc/F with the exponent fixed at 1 (Shao 2023 Results 3.2 and Table 2, ",
        "'Weight effect on CL' / 'Weight effect on Vc'). Shao 2023 lists as a ",
        "limitation that weight was carried as a time-fixed baseline value ",
        "even though it changes over MDR-TB treatment; the packaged model ",
        "therefore treats WT as time-fixed, but a user with longitudinal ",
        "weights may supply it as a time-varying column."
      ),
      source_name        = "Weight"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Baseline serum albumin in SI units (g/L), the units Shao 2023 reports ",
        "and calibrated the coefficient in, so no unit conversion is applied ",
        "inside model(). Normalised to the development-cohort median of ",
        "37.1 g/L (Shao 2023 Table 1, development cohort: median 37.1 g/L, IQR ",
        "32.0-40.6 g/L) and entered on CL/F as a power function with an ",
        "estimated exponent of 3.76 (Shao 2023 Eq. 5 and Table 2, 'Albumin ",
        "effect on CL', RSE 14.2%, 95% CI 2.71-4.80). The exponent is steep: ",
        "across the development cohort's albumin IQR it spans a 2.4-fold range ",
        "in CL/F. Shao 2023 attributes the albumin effect to the high plasma ",
        "protein binding of bedaquiline (Discussion) and lists the time-fixed ",
        "treatment of albumin as a model limitation."
      ),
      source_name        = "Albumin"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 1L,
    n_observations = paste0(
      "1,205 quantifiable plasma bedaquiline concentrations (0.04-5.96 mg/L); ",
      "5 further samples below the 10 ng/mL LLOQ were discarded ",
      "(Shao 2023 Results 3.2)."
    ),
    age_range      = "18-70 years by protocol; observed median 44 years (IQR 34-53) (Shao 2023 Table 1)",
    age_median     = "44 years",
    weight_range   = "IQR 59-74 kg (Shao 2023 Table 1); no full range reported",
    weight_median  = "67 kg",
    height_median  = "165 cm (IQR 158-173) (Shao 2023 Table 1)",
    albumin_median = "37.1 g/L (IQR 32.0-40.6) (Shao 2023 Table 1)",
    sex_female_pct = 29.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Adults with multidrug-resistant tuberculosis (MDR-TB) confirmed by ",
      "phenotypic drug-susceptibility testing; 16.4% pre-XDR, 16.4% with a ",
      "pulmonary cavity, 16.4% with type 2 diabetes mellitus. Patients with ",
      "abnormal liver or kidney function, pregnancy, or HIV / hepatitis B / ",
      "hepatitis C infection were excluded (Shao 2023 Methods 2.1, Table 1)."
    ),
    dose_range     = paste0(
      "Oral bedaquiline 400 mg once daily for the first 2 weeks followed by ",
      "200 mg three times weekly for the following 22 weeks (the WHO ",
      "recommended regimen), within a standardized bedaquiline-containing ",
      "MDR-TB regimen."
    ),
    regions        = "China (Guizhou, Henan, Jiangsu and Sichuan Provinces), June 2016 to June 2019",
    co_medication  = paste0(
      "Development cohort: moxifloxacin or levofloxacin plus linezolid plus a ",
      "background regimen; validation cohort: moxifloxacin, linezolid, ",
      "clofazimine and cycloserine. Shao 2023 lists the unmodelled ",
      "drug-drug interaction potential of the companion drugs (notably ",
      "clofazimine) as a limitation."
    ),
    notes          = paste0(
      "The 55 subjects of the richly-sampled development cohort are the ",
      "subjects this PK model was fit to; plasma was drawn predose and at 1, ",
      "2, 3, 4, 5, 6, 8, 12, 18 and 24 h after dosing on treatment weeks 2 and ",
      "4. A separate prospective validation cohort of 159 MDR-TB patients ",
      "(sparse sampling: predose and 2, 4 and 6 h on week 4) was used only for ",
      "Bayesian post-hoc exposure estimation and CART threshold derivation, ",
      "not for estimating the parameters packaged here (Shao 2023 Methods 2.1, ",
      "2.3). Model built in Phoenix NLME 8.0."
    )
  )

  ini({
    # ==================================================================
    # Structural disposition -- Shao 2023 Table 2, 'Typical value
    # parameter of population'. All values are apparent (per bioavail-
    # ability F); the paper reports CL/F, Vc/F, Vp1/F, CLp1/F, Vp2/F and
    # CLp2/F. The canonical mapping is Vp1 -> vp, CLp1 -> q,
    # Vp2 -> vp2, CLp2 -> q2.
    # ==================================================================
    lcl  <- log(3.57);     label("Apparent clearance CL/F at reference covariates (L/h)")                       # Shao 2023 Table 2: CL/F = 3.57 L/h (RSE 11.9%, 95% CI 2.73-4.40)
    lvc  <- log(336.97);   label("Apparent central volume of distribution Vc/F at reference weight (L)")        # Shao 2023 Table 2: Vc/F = 336.97 L (RSE 12.6%, 95% CI 253.86-420.08)
    lvp  <- log(2839.13);  label("Apparent first peripheral volume of distribution Vp1/F (L)")                  # Shao 2023 Table 2: Vp1/F = 2839.13 L (RSE 46.3%, 95% CI 258.19-5420.06)
    lq   <- log(2.97);     label("Apparent inter-compartmental clearance CLp1/F to peripheral1 (L/h)")          # Shao 2023 Table 2: CLp1/F = 2.97 L/h (RSE 41.40%, 95% CI 0.56-5.39)
    lvp2 <- log(1391.89);  label("Apparent second peripheral volume of distribution Vp2/F (L)")                 # Shao 2023 Table 2: Vp2/F = 1391.89 L (RSE 37.41%, 95% CI 370.16-2413.63)
    lq2  <- log(9.81);     label("Apparent inter-compartmental clearance CLp2/F to peripheral2 (L/h)")          # Shao 2023 Table 2: CLp2/F = 9.81 L/h (RSE 18.86%, 95% CI 6.18-13.44)

    # ==================================================================
    # Dual zero-order absorption input. Shao 2023 Methods 2.2 states that
    # the absorption rate constant was FIXED at 1,000 1/h and that DUR1,
    # DUR2 and FR1 were also fixed from McLeay 2014 "for a better fit";
    # Table 2 accordingly reports FR1, DUR1 and DUR2 with no RSE and no
    # confidence interval. Ka = 1,000 1/h makes the depot half-life
    # 0.0007 h, so the zero-order fill of the depot is passed on to the
    # central compartment essentially without delay -- which is exactly
    # the paper's stated intent ("resulting in the dose via the input
    # compartment describing an initial zero-order (rather than
    # first-order) input").
    # ==================================================================
    lka     <- fixed(log(1000));   label("First-order transfer rate constant out of the depot Ka (1/h)")            # Shao 2023 Methods 2.2: "we fixed the absorption rate as 1,000 h-1 based on a published population PK model (McLeay et al., 2014)"
    lfdepot <- fixed(log(0.585));  label("Fraction of the dose entering the depot compartment FR1 (fraction)")      # Shao 2023 Table 2: FR1 = 58.5% (no RSE / CI reported); Methods 2.2 "the fraction of dose into the depot compartment (58.5%) were also fixed from this reference"
    ld1     <- fixed(log(2.22));   label("Zero-order duration of the infusion into the depot DUR1 (h)")             # Shao 2023 Table 2: DUR1 = 2.22 h (no RSE / CI reported); Methods 2.2 "The duration of infusion into the depot compartment (2.22 h) ... were also fixed"
    ld2     <- fixed(log(1.48));   label("Zero-order duration of the infusion into the central compartment DUR2 (h)")  # Shao 2023 Table 2: DUR2 = 1.48 h (no RSE / CI reported); Methods 2.2 "the duration of infusion into Vc/F (1.48 h) ... were also fixed"

    # Tlag is the lag before absorption starts (it delays the depot arm);
    # Tlag,add is, per the Table 2 footnote, the ADDITIONAL lag "prior to
    # administration of the remaining dose into Vc/F", so the central arm
    # starts at Tlag + Tlag,add = 1.00 + 2.76 = 3.76 h. Because the depot
    # arm finishes at 1.00 + 2.22 = 3.22 h, the additive reading leaves a
    # 0.54 h gap with no input at all -- which is what produces the
    # double absorption peak the dual-input structure was introduced to
    # capture (Shao 2023 Results 3.2). Both lag times are estimated.
    ltlag     <- log(1.00);  label("Lag time before absorption starts, applied to the depot arm Tlag (h)")                 # Shao 2023 Table 2: Tlag = 1.00 h (RSE 9.1%, 95% CI 0.82-1.18)
    ltlag_add <- log(2.76);  label("Additional lag time before the direct zero-order input into Vc/F Tlag,add (h)")        # Shao 2023 Table 2: Tlag,add = 2.76 h (RSE 21.0%, 95% CI 1.62-3.90)

    # ==================================================================
    # Covariate effects -- Shao 2023 Eq. 4, Eq. 5 and Table 2
    # 'Covariable effect'. The two weight exponents were fixed to the
    # allometric values; the albumin exponent was estimated.
    # ==================================================================
    e_wt_cl  <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")     # Shao 2023 Table 2: 'Weight effect on CL' = 0.75 (no RSE / CI); Results 3.2 "the covariate effects of body weight on the central volume of distribution and clearance were fixed to 1 and 0.75, respectively, according to the principles of allometric scaling"
    e_wt_vc  <- fixed(1);    label("Allometric exponent of body weight on Vc/F (unitless)")     # Shao 2023 Table 2: 'Weight effect on Vc' = 1 (no RSE / CI); Results 3.2, same sentence as above
    e_alb_cl <- 3.76;        label("Power exponent of serum albumin on CL/F (unitless)")               # Shao 2023 Table 2: 'Albumin effect on CL' = 3.76 (RSE 14.2%, 95% CI 2.71-4.80); Eq. 5

    # ==================================================================
    # Between-subject variability -- Shao 2023 Table 2, 'BSV'. The paper
    # states an exponential IIV model, theta_i = theta_TV * exp(eta_i)
    # with eta_i ~ N(0, omega^2) (Eq. 1), and tabulates omega^2 directly,
    # so the values below are used as reported with no CV% conversion.
    # omega^2_CL = 1.33 is very large (167% CV) and is itself imprecisely
    # estimated (RSE 167%); it is used as reported.
    # ==================================================================
    etalvc ~ 0.38   # Shao 2023 Table 2: omega^2 V = 0.38 (RSE 67.99%) -> 65% CV on Vc/F
    etalcl ~ 1.33   # Shao 2023 Table 2: omega^2 CL = 1.33 (RSE 166.76%) -> 167% CV on CL/F

    # ==================================================================
    # Residual variability. Shao 2023 Results 3.2 states plainly that "A
    # proportional error model was used to evaluate the residual
    # variability", so the single tabulated residual variance of 0.23 is
    # encoded as a proportional error with SD sqrt(0.23) = 0.4796. The
    # Table 2 row label reads 'sigma^2 add'; see the vignette Errata for
    # why the prose is followed instead (an additive SD of 0.48 mg/L is
    # not compatible with the 0.04 mg/L lower end of the observed
    # concentration range).
    # ==================================================================
    propSd <- sqrt(0.23); label("Proportional residual error on plasma bedaquiline Cc (fraction)")  # Shao 2023 Table 2: residual variance = 0.23 (RSE 2.75%, 95% CI 0.22-0.24); Results 3.2 "A proportional error model was used"
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. Shao 2023 Eq. 4 and Eq. 5:
    #      Vc/F = 336.97 * (Weight / 67)^1    * exp(eta_V)
    #      CL/F = 3.57   * (Weight / 67)^0.75 * (Albumin / 37.1)^3.76
    #                                         * exp(eta_CL)
    #    The centring constants are the development-cohort medians named
    #    in the sentence immediately after Eq. 5 ("The median weight was
    #    67 kg. The median albumin was 37.1 g/L") and confirmed by
    #    Table 1. The Table 2 footnote instead says "The reference weight
    #    was 59 kg"; 59 kg is the median weight of the McLeay 2014
    #    cohort that Shao 2023's own Discussion contrasts against
    #    ("the higher weight (median: 67 vs. 59 kg) of our subjects"), so
    #    the footnote is treated as a carry-over error and the printed
    #    equations are followed. See vignette Errata.
    # ------------------------------------------------------------------
    cl  <- exp(lcl + etalcl) * (WT / 67)^e_wt_cl * (ALB / 37.1)^e_alb_cl
    vc  <- exp(lvc + etalvc) * (WT / 67)^e_wt_vc
    vp  <- exp(lvp)
    q   <- exp(lq)
    vp2 <- exp(lvp2)
    q2  <- exp(lq2)

    ka       <- exp(lka)
    fdepot   <- exp(lfdepot)
    d1       <- exp(ld1)
    d2       <- exp(ld2)
    tlag     <- exp(ltlag)
    tlag_add <- exp(ltlag_add)

    # 2. Micro-constants
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ------------------------------------------------------------------
    # 3. Three-compartment disposition with a depot feeding the central
    #    compartment. Compartment roles follow Shao 2023 Table 2:
    #    peripheral1 carries Vp1/F and exchanges via CLp1/F, peripheral2
    #    carries Vp2/F and exchanges via CLp2/F.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ------------------------------------------------------------------
    # 4. Dual zero-order input. Each oral administration is encoded by
    #    the USER as TWO dose records at the same nominal time, both
    #    carrying the full dose amount and both with rate = -2 so that
    #    rxode2 uses the modelled durations below:
    #      (a) cmt = "depot"   -> fraction FR1  over DUR1, lagged by Tlag
    #      (b) cmt = "central" -> fraction 1-FR1 over DUR2, lagged by
    #                             Tlag + Tlag,add
    #    The f() multipliers do the FR1 split, so the two records must
    #    each carry the whole dose, not half of it.
    # ------------------------------------------------------------------
    f(depot)      <- fdepot
    dur(depot)    <- d1
    alag(depot)   <- tlag

    f(central)    <- 1 - fdepot
    dur(central)  <- d2
    alag(central) <- tlag + tlag_add

    # 5. Observation. Dose in mg and vc in L gives mg/L directly, which
    #    is the unit Shao 2023 reports concentrations in.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
