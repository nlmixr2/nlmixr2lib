Courlet_2023_cabamiquine_pk <- function() {
  description <- paste(
    "Population PK model for the Plasmodium elongation factor 2 (eEF2)",
    "inhibitor cabamiquine (formerly DDD107498 / M5717) in healthy adult",
    "men, pooled from a first-in-human single-ascending-dose / induced",
    "blood stage malaria (IBSM) trial and a sporozoite-challenge (SpzCh)",
    "chemoprophylaxis trial (Courlet 2023 Supplementary Material 1).",
    "Three-compartment apparent disposition with linear elimination; a",
    "Savic transit-absorption chain feeds a depot that is absorbed into",
    "central at rate ka; and an enterohepatic recycling loop moves drug",
    "from central into a recycling (gallbladder / bile) compartment at",
    "rate k2g and releases it back into the absorption depot at rate kg1",
    "once the depot emptying time MTIME is reached, reproducing the",
    "secondary concentration peak observed 24-30 h post-dose. The observed",
    "lack of dose proportionality is captured empirically by a power effect",
    "of the administered dose on the apparent central volume V2/F",
    "(exponent -0.50); allometric body-weight scaling is fixed at 0.75 on",
    "the apparent clearances and 1 on the apparent volumes. This is the",
    "parsimonious PK model that was carried forward, with individual",
    "parameters fixed, into the parasitemia PK/PD model; see",
    "modellib('Courlet_2023_cabamiquine').",
    sep = " "
  )
  reference <- paste(
    "Courlet P, Wilkins JJ, Oeuvray C, Gao W, Khandelwal A.",
    "Semi-mechanistic population pharmacokinetic/pharmacodynamic modeling",
    "of a Plasmodium elongation factor 2 inhibitor cabamiquine for",
    "prevention and cure of malaria.",
    "Antimicrob Agents Chemother. 2023;67(12):e00891-23.",
    "doi:10.1128/aac.00891-23.",
    sep = " "
  )
  vignette <- "Courlet_2023_cabamiquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Courlet 2023 Fig. 2 (left panel,
  # "Population pharmacokinetic model"), which draws Depot, Central (V2,
  # Ccab), Peripheral 1 (V3), Peripheral 2 (V4), and Recycling.
  compartmentData <- list(
    depot       = list(analyte = "cabamiquine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral1 = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral2 = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    gallbladder = list(analyte = "cabamiquine", units = "mg", specimen = "bile",                verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with exponents fixed at 0.75 on the apparent",
        "clearances (CL/F, Q2/F, Q3/F) and 1 on the apparent volumes",
        "(V2/F, V3/F, V4/F) per Courlet 2023 Supplementary Material 1 rows",
        "'Weight on CL/F ... 0.75 Fixed' through 'Weight on Q3/F ... 0.75",
        "Fixed' and the Results statement 'Allometric scaling was applied",
        "to clearance and volume parameters (19)'. The paper does not state",
        "the reference (centering) weight; the packaged model uses the",
        "conventional 70 kg. See vignette Assumptions and deviations."
      ),
      source_name        = "WT"
    ),
    DOSE_CABAMIQUINE_MG = list(
      description        = "Administered cabamiquine single oral dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject single oral dose. Enters the model as an empirical",
        "power effect on the apparent central volume,",
        "V2/F = 2363 * (WT/70)^1 * DOSE_CABAMIQUINE_MG^e_dose_vc, which the",
        "authors introduced to describe the observed greater-than-dose-",
        "proportional exposure (Results 'Population PK modeling': 'it was",
        "ultimately necessary to include an empirical covariate effect of",
        "dose on V2/F in the model'). The paper never reports the reference",
        "(centering) dose of the power term; the packaged model uses the",
        "uncentred dose in mg (i.e. reference dose 1 mg), which is the only",
        "reading compatible with the reported 146-193 h terminal half-life.",
        "Doses are FREE BASE milligrams, and so is the `amt` on the dose",
        "record. Study 1 (IBSM) administered cabamiquine as the succinate",
        "salt and study 2 (SpzCh) as free base, so salt doses must be",
        "converted before use: 1 mg succinate salt = 0.797 mg free base",
        "(Methods 'Data set'), i.e. the 150 / 400 / 800 mg salt IBSM cohorts",
        "are 119.6 / 318.8 / 637.6 mg free base. The paper gives that",
        "conversion factor precisely so the two studies can be pooled on a",
        "common basis; the vignette shows that the free-base reading, not",
        "the salt reading, reproduces the published recrudescence rates. See",
        "vignette Assumptions and deviations."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 101L,
    n_studies      = 2L,
    age_range      = "18-55 years",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Healthy malaria-naive adult volunteers. Study 1 part 2 subjects were",
      "inoculated with Plasmodium falciparum-infected erythrocytes (strain",
      "3D7) 8 days before dosing (induced blood stage malaria, IBSM);",
      "study 2 subjects were inoculated with approximately 3,200",
      "sporozoites (strain NF54) 2 h or 96 h before dosing (sporozoite",
      "challenge, SpzCh)."
    ),
    dose_range     = paste(
      "Single oral doses across 14 dose levels. Study 1 (NCT03261401)",
      "administered cabamiquine succinate salt (1 mg salt = 0.797 mg free",
      "base), including the 150, 400, and 800 mg IBSM cohorts; study 2",
      "administered free base at 30, 60, 80, 100, and 200 mg."
    ),
    regions        = "Australia (QIMR Berghofer) and the Netherlands",
    trial_registration = "ClinicalTrials.gov NCT03261401 (study 1)",
    notes          = paste(
      "The PK data set comprised 1,823 evaluable cabamiquine blood",
      "concentrations from 101 healthy subjects across 14 doses (Results",
      "'PK/PD data'). Parasitemia data used for the PK/PD layer came from",
      "22 IBSM participants and 39 SpzCh participants; the Discussion notes",
      "the pooled parasitemia population was 61 healthy, malaria-naive",
      "adult Caucasian men. The paper does not tabulate baseline weight or",
      "age distributions, so weight_range / age_median are omitted rather",
      "than guessed."
    )
  )

  ini({
    # =========================================================================
    # Structural PK. All values are Courlet 2023 Supplementary Material 1
    # ("PK model parameter estimates"), which is the parameter table of
    # record for the population PK model -- the main-text tables carry only
    # the PK/PD parameters. Estimates are apparent (F-relative) values.
    # =========================================================================
    lcl    <- log(17.8);   label("Apparent cabamiquine clearance CL/F for a 70 kg subject (L/h)")                 # Suppl. Material 1 row 'Clearance (CL/F, L/h) = 17.8 (RSE 5.64%)'
    lvc    <- log(2363);   label("Apparent central volume V2/F at a 1 mg dose for a 70 kg subject (L)")           # Suppl. Material 1 row 'Central volume of distribution (V2/F, L) = 2363 (RSE 3.23%)'
    lvp    <- log(2051);   label("Apparent first peripheral volume V3/F for a 70 kg subject (L)")                 # Suppl. Material 1 row 'Peripheral volume of distribution 1 (V3/F, L) = 2051 (RSE 7.76%)'
    lq     <- log(6.37);   label("Apparent first intercompartmental clearance Q2/F for a 70 kg subject (L/h)")    # Suppl. Material 1 row 'Intercompartmental clearance 1 (Q2/F, L/h) = 6.37 (RSE 0.0589%)'
    lvp2   <- log(2548);   label("Apparent second peripheral volume V4/F for a 70 kg subject (L)")                # Suppl. Material 1 row 'Peripheral volume of distribution 2 (V4/F, L) = 2548 (RSE 3.23%)'
    lq2    <- log(58.5);   label("Apparent second intercompartmental clearance Q3/F for a 70 kg subject (L/h)")   # Suppl. Material 1 row 'Intercompartmental clearance 2 (Q3/F, L/h) = 58.5 (RSE 3.08%)'
    lka    <- log(8.27);   label("First-order absorption rate constant ka, depot to central (1/h)")               # Suppl. Material 1 row 'Absorption rate constant (ka /h) = 8.27 (RSE 21.7%)'

    # Savic transit-absorption chain. Both MTT and ktr are estimated, so the
    # number of transit compartments is derived as nn = mtt * ktr - 1 =
    # 0.21 * 13.7 - 1 = 1.877 (non-integer, as the Savic 2007 analytical
    # solution allows).
    lktr   <- log(13.7);   label("Transit rate between absorption compartments ktr (1/h)")                        # Suppl. Material 1 row 'Transit rate between compartments (ktr /h) = 13.7 (RSE 5.66%)'
    lmtt   <- log(0.21);   label("Mean transit time MTT through the transit chain (h)")                           # Suppl. Material 1 row 'Mean transit time (MTT, h) = 0.21 (RSE 11.6%)'

    # Enterohepatic recycling loop (Courlet 2023 Fig. 2 left panel:
    # "Recycling compartment empties into depot at time MTIME"). k2g moves
    # drug from central into the recycling (bile) compartment continuously;
    # kg1 releases it back into the absorption depot, gated on the depot
    # emptying time MTIME. Both rate constants were fixed during estimation.
    lkbm   <- fixed(log(0.0039)); label("Central to recycling (biliary) transfer rate constant k2g (1/h)")        # Suppl. Material 1 row 'Central to depot transfer rate constant (k2g /h) = 0.0039 Fixed'
    lkehc  <- fixed(log(15));     label("Recycling to depot release rate constant kg1 (1/h), gated at MTIME")     # Suppl. Material 1 row 'Depot to absorption transfer rate constant (kg1 /h) = 15 Fixed*' (footnote: 'Fixed to a high number; transfer is near-instant.')
    lmtime <- log(25.05);         label("Depot emptying time MTIME, onset of the recycling release (h)")          # Suppl. Material 1 row 'Depot emptying time (h) = 25.05 (RSE 13.0%)'

    # Relative oral bioavailability is anchored at F = 1 because CL and the
    # volumes are reported as apparent (CL/F, V/F) values; only its IIV was
    # estimated (Results: 'IIV was estimated on apparent clearance (CL/F),
    # apparent central volume of distribution (V2/F), bioavailability (F)
    # and absorption rate (ka), transit rate between compartments (ktr),
    # mean transit time (MTT), and central to depot transfer rate constant
    # (k2g)').
    lfdepot <- fixed(log(1)); label("Reference relative oral bioavailability F (fraction)")                       # Apparent-parameter anchor; Suppl. Material 1 reports IIV on F but no typical value

    # Fixed allometric exponents, applied to the apparent clearances and the
    # apparent volumes. The paper does not state the reference weight; the
    # conventional 70 kg is used (see vignette Assumptions and deviations).
    e_wt_cl   <- fixed(0.75); label("Allometric exponent on CL/F, Q2/F and Q3/F with body weight (unitless)")     # Suppl. Material 1 rows 'Weight on CL/F = 0.75 Fixed', 'Weight on Q2/F = 0.75 Fixed', 'Weight on Q3/F = 0.75 Fixed'
    e_wt_vc   <- fixed(1);    label("Allometric exponent on V2/F, V3/F and V4/F with body weight (unitless)")     # Suppl. Material 1 rows 'Weight on V2/F = 1 Fixed', 'Weight on V3/F = 1 Fixed', 'Weight on V4/F = 1 Fixed'

    # Empirical dose effect on the apparent central volume. The main text
    # quotes -0.530 for this coefficient; Supplementary Material 1 -- the
    # parameter table of record, and the artefact the main text points the
    # reader to -- reports -0.50 with RSE 5.18%. The supplement value is
    # used; the conflict is recorded in the vignette Errata.
    e_dose_vc <- -0.50;       label("Power exponent of administered dose on V2/F (unitless)")                     # Suppl. Material 1 row 'Dose on V2/F = -0.50 (RSE 5.18%)'

    # =========================================================================
    # Inter-individual variability. Supplementary Material 1 reports IIV as
    # a coefficient of variation in percent (table footnote: 'IIV:
    # inter-individual variability, reported as coefficient of variation
    # (%)'), converted to log-normal variance via omega^2 = log(1 + CV^2).
    # =========================================================================
    etalcl     ~ 0.07059  # Suppl. Material 1 'IIV on CL/F = 27%'  -> omega^2 = log(1 + 0.27^2) = 0.07059
    etalvc     ~ 0.06062  # Suppl. Material 1 'IIV on V2/F = 25%'  -> omega^2 = log(1 + 0.25^2) = 0.06062
    etalka     ~ 2.51231  # Suppl. Material 1 'IIV on ka = 340%'   -> omega^2 = log(1 + 3.40^2) = 2.51231
    etalkbm    ~ 0.50974  # Suppl. Material 1 'IIV on k2g = 82%'   -> omega^2 = log(1 + 0.82^2) = 0.50974
    etalktr    ~ 0.06062  # Suppl. Material 1 'IIV on ktr = 25%'   -> omega^2 = log(1 + 0.25^2) = 0.06062
    etalmtt    ~ 0.79470  # Suppl. Material 1 'IIV on MTT = 110%'  -> omega^2 = log(1 + 1.10^2) = 0.79470
    etalfdepot ~ 0.05161  # Suppl. Material 1 'IIV on F = 23%'     -> omega^2 = log(1 + 0.23^2) = 0.05161

    # =========================================================================
    # Residual error. Supplementary Material 1 reports a single 'Residual
    # error (%) = 17', i.e. a proportional residual on the linear blood
    # concentration scale.
    # =========================================================================
    propSd <- 0.17; label("Proportional residual error on blood cabamiquine concentration (fraction)")            # Suppl. Material 1 row 'Residual error (%) = 17 (RSE 1.88%)'
  })

  model({
    # =========================================================================
    # 1. Individual parameters. Allometric weight scaling applies to the
    #    apparent clearances (exponent 0.75) and the apparent volumes
    #    (exponent 1); the empirical dose power term applies to V2/F only.
    # =========================================================================
    allcl <- (WT / 70)^e_wt_cl
    allv  <- (WT / 70)^e_wt_vc

    cl    <- exp(lcl + etalcl) * allcl
    q     <- exp(lq)           * allcl
    q2    <- exp(lq2)          * allcl
    vc    <- exp(lvc + etalvc) * allv * DOSE_CABAMIQUINE_MG^e_dose_vc
    vp    <- exp(lvp)          * allv
    vp2   <- exp(lvp2)         * allv

    ka     <- exp(lka  + etalka)
    ktr    <- exp(lktr + etalktr)
    mtt    <- exp(lmtt + etalmtt)
    kbm    <- exp(lkbm + etalkbm)
    kehc   <- exp(lkehc)
    # Bare name is tmtime, not mtime: `mtime` is a reserved rxode2 keyword
    # (the modeled-event-time function) and cannot be used as a variable.
    tmtime <- exp(lmtime)
    fdepot <- exp(lfdepot + etalfdepot)

    # Number of transit compartments implied by the Savic parameterisation,
    # nn = mtt * ktr - 1. The typical value is 0.21 * 13.7 - 1 = 1.877.
    #
    # Both MTT and ktr carry IIV (110% and 25% CV), so nn is a random
    # quantity, and roughly a third of random draws put it below 1. The
    # Savic input is proportional to (ktr * t)^nn, whose derivative at the
    # dose time behaves like t^(nn - 1): it is unbounded for nn < 1,
    # discontinuous at nn = 0, and negative nn makes the gamma density
    # undefined altogether. Those draws are both unphysical (a transit
    # chain shorter than one compartment) and unsolvable -- they make the
    # ODE solver fail outright for the affected subjects.
    #
    # nn is therefore floored at 1, the smallest value for which the input
    # is continuously differentiable at the dose time. The typical-value
    # profile (nn = 1.877) is unaffected. See vignette Assumptions and
    # deviations, which quantifies how many draws hit the floor.
    nn <- mtt * ktr - 1
    nn <- max(nn, 1)

    # =========================================================================
    # 2. Absorption input. The Savic 2007 analytical transit-chain density is
    #    written out explicitly rather than via the rxode2 transit() macro:
    #    transit() rescales its internal dose lookup by bioavailability, so
    #    the transit() + f(depot) <- 0 idiom silently delivers zero dose in
    #    nlmixr2 UI-form models. podo(depot) is not F-adjusted and returns
    #    the full dose amount, so the explicit form below is used and the
    #    delivered mass is asserted in the vignette.
    # =========================================================================
    #    tad() and podo() both return NA until the first dose record, and that
    #    NA would propagate into every state, so the input and the recycling
    #    gate are evaluated only once a dose has been given. The `else` branch
    #    also covers a placebo arm, which has no dose record at all.
    #
    #    Recycling release gate: the recycling compartment fills continuously
    #    from central and empties into the absorption depot once the depot
    #    emptying time MTIME has elapsed since the dose (Fig. 2 left panel).
    #    This reproduces the secondary peak reported between 24 and 30 h.
    tdose <- tad(depot)
    if (tdose > 0) {
      transit_in <- exp(log(podo(depot)) + log(fdepot) + log(ktr) +
                          nn * log(ktr * tdose) - ktr * tdose - lgamma(nn + 1))
      gate <- (tdose >= tmtime)
    } else {
      transit_in <- 0
      gate <- 0
    }

    # =========================================================================
    # 3. ODE system (Courlet 2023 Fig. 2, left panel).
    # =========================================================================
    d/dt(depot)       <-  transit_in + gate * kehc * gallbladder - ka * depot
    d/dt(central)     <-  ka * depot -
                            cl * central / vc -
                            q  * central / vc + q  * peripheral1 / vp -
                            q2 * central / vc + q2 * peripheral2 / vp2 -
                            kbm * central
    d/dt(peripheral1) <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2
    d/dt(gallbladder) <-  kbm * central - gate * kehc * gallbladder

    # 4. The bolus into depot is suppressed; the transit chain above delivers
    #    the whole dose (scaled by the individual bioavailability fdepot).
    f(depot) <- 0

    # =========================================================================
    # 5. Observation. Amount (mg) / volume (L) = mg/L; multiply by 1000 to
    #    report ng/mL, the unit used for cabamiquine blood concentrations and
    #    for the EC50 / MIC values in the companion PK/PD model.
    # =========================================================================
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
