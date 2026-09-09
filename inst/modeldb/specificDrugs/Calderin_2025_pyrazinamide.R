Calderin_2025_pyrazinamide <- function() {
  description <- paste(
    "One-compartment population PK model for oral pyrazinamide in plasma and",
    "lumbar cerebrospinal fluid (CSF) in South African adults with",
    "HIV-associated tuberculous meningitis (LASER-TBM PK substudy, Calderin",
    "2025). Absorption is a Savic transit chain (mean transit time 0.291 h,",
    "estimated chain length 4.25) feeding a first-order depot (ka 2.5 1/h);",
    "elimination is first order from the central compartment. Clearance",
    "(4.19 L/h) and central volume (45.0 L) are allometrically scaled on",
    "fat-free mass with a 45 kg reference and fixed 0.75 / 1 exponents (FFM",
    "beat total body weight by dOFV -15.5 vs -4.5). Clearance is 30.2%",
    "higher at the day-28 PK visit than at the day-3 visit. CSF is a",
    "Sheiner-style effect compartment holding a concentration, equilibrating",
    "with plasma at 1.05 1/h (equilibration half-life 0.66 h) toward a",
    "CSF-to-plasma pseudo-partition coefficient of 1.05, i.e. CSF exposure",
    "matches plasma. Random effects are between-subject variability on",
    "clearance (18.5%) and five-occasion between-occasion variability on",
    "bioavailability (15.8%), absorption rate constant (87.3%) and mean",
    "transit time (102%); all reported percentages are the omega standard",
    "deviation on the log scale. Residual error is combined proportional",
    "plus additive, separately for plasma (8.33%, 0.04 mg/L) and CSF",
    "(11.4%, 0.0468 mg/L). High-dose rifampicin (35 mg/kg) did not affect",
    "pyrazinamide plasma PK or CSF penetration."
  )
  reference <- paste(
    "Calderin JM, Wasserman S, Resendiz-Galvan JE, Abdelgawad N, Davis A,",
    "Stek C, Wiesner L, Wilkinson RJ, Denti P (2025).",
    "Population pharmacokinetics of pyrazinamide and isoniazid in plasma and",
    "cerebrospinal fluid from South African adults with tuberculous",
    "meningitis. Antimicrob Agents Chemother 69(8):e00099-25.",
    "doi:10.1128/aac.00099-25",
    sep = " "
  )
  vignette <- "Calderin_2025_pyrazinamide_isoniazid_tbm"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state is an exception to the usual "states hold an amount"
  # rule. Supplementary S2 writes the effect compartment directly in
  # concentration units, dC_CSF/dt = k_Plasma-CSF * (PPC * C_Plasma - C_CSF),
  # and the S9 control stream integrates that same equation as
  # DADT(3) = KE0*(PPC*C2 - A(3)) with the CSF observation read as
  # CE = A(3) (not A(3)/V). So `csf` carries mg/L, not mg. Same convention as
  # the sibling LASER-TBM model Abdelgawad_2024_linezolid.R.
  compartmentData <- list(
    depot = list(
      analyte = "pyrazinamide", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "pyrazinamide", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    csf = list(
      analyte = "pyrazinamide", units = "mg/L",
      specimen = "CSF", verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = paste(
        "Fat-free mass, computed with the Janmahasatian et al. formula from",
        "total body weight, height and sex. Allometric scaling on FFM was",
        "retained over total body weight (dOFV -15.53 vs -4.47)."
      ),
      units = "kg",
      type = "continuous",
      source_name = "FFM",
      reference_value = 45,
      notes = paste(
        "Cohort median FFM 45 kg (range 30-59), Table 1. The S9 control",
        "stream hardcodes TVFFM = 45 and forms ALLMCL_FFM = (FFM/45)^0.75",
        "and ALLMV_FFM = (FFM/45). Missing heights needed for the",
        "Janmahasatian formula were imputed with the Johansson and Karlsson",
        "regression given in supplementary S3."
      )
    ),
    DAY28 = list(
      description = paste(
        "Day-28 PK visit indicator: 1 = the record belongs to the second",
        "(day 28) PK sampling visit, 0 = the first (day 3) visit."
      ),
      units = "(binary)",
      type = "binary",
      source_name = "PK_VISIT",
      reference_category = "0 (day-3 visit)",
      notes = paste(
        "The S9 control stream carries PK_VISIT with the literal values 3",
        "and 28 and selects PK_VISIT_CL = 1 or 1 + THETA(7) from it, so the",
        "column is a two-level visit indicator rather than a continuous",
        "treatment-duration clock. Sampling was on day 3 (+/- 2) and day 28",
        "(+/- 2) after study enrolment; TB treatment had started 1-3 days",
        "before enrolment, so the visits sit at a median of 4 and 30 days on",
        "rifampicin respectively (Table 1)."
      )
    ),
    OCC = list(
      description = paste(
        "Occasion index for between-occasion variability. An occasion is a",
        "dosing event and its subsequent observations."
      ),
      units = "(count)",
      type = "categorical",
      source_name = "OCC",
      notes = paste(
        "The S9 control stream multiplexes the BOV etas with IF (OCC==k)",
        "blocks over five occasions (ETA 13-17 bioavailability, ETA 18-22",
        "ka, ETA 23-27 MTT), each declared $OMEGA BLOCK(1) SAME after the",
        "first. The sibling LASER-TBM linezolid control stream maps",
        "occasions 1, 2 and 5 to the day-3 visit and occasions 3 and 4 to",
        "the day-28 visit; the same trial and sampling schedule apply here."
      )
    )
  )

  population <- list(
    n_subjects = 49,
    n_studies = 1,
    species = "human",
    age_range = "39 years (range 25-78)",
    weight_range = "60.0 kg (range 30.0-107)",
    sex_female_pct = 45,
    race_ethnicity = "South African; not further reported",
    disease_state = paste(
      "HIV-associated tuberculous meningitis; all participants received",
      "adjunctive dexamethasone"
    ),
    dose_range = paste(
      "Pyrazinamide 25 mg/kg once daily by WHO weight band, given as an oral",
      "fixed-dose combination"
    ),
    regions = "South Africa (Cape Town and Gqeberha)",
    notes = paste(
      "PK substudy nested in the open-label randomised phase 2A LASER-TBM",
      "trial. 414 plasma and 44 CSF pyrazinamide concentrations from 49",
      "participants at the day-3 visit and 34 at the day-28 visit. Arms:",
      "standard of care (rifampicin 10 mg/kg) versus high-dose rifampicin",
      "(35 mg/kg) plus linezolid, with or without aspirin; neither the",
      "rifampicin dose nor aspirin affected pyrazinamide PK. Intensive",
      "sampling on day 3 (pre-dose, 0.5, 1, 2, 3, 6, 8-10, 24 h) and sparse",
      "sampling on day 28 (pre-dose, 2, 4 h), with one lumbar CSF sample per",
      "visit randomised to a 1-3, 3-6, 6-10 or 24 h window. Plasma LLOQ",
      "0.200 mg/L (0.97% BLQ), CSF LLOQ 0.234 mg/L (2.3% BLQ). Median",
      "fat-free mass 45 kg (range 30-59). The measured unbound plasma",
      "fraction of pyrazinamide was 93.3%; the model is written for total",
      "concentrations."
    )
  )

  ini({
    # --- Disposition. Table 2 typical values are quoted for the typical
    # individual in the cohort, i.e. at fat-free mass 45 kg (footnote b).
    lcl <- log(4.19)
    label("Log of clearance at the day-3 visit and FFM 45 kg (L/h)")                              # Table 2 Clearance 4.19 (3.86-4.45)
    lvc <- log(45.0)
    label("Log of central volume of distribution at FFM 45 kg (L)")                               # Table 2 Central volume of distribution 45.0 (43.4-46.6)

    # Allometric exponents fixed a priori, not estimated: "the exponents for
    # clearance and volume were fixed to 0.75 and 1, respectively".
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on fat-free mass for clearance")                             # Methods, Pharmacokinetic modeling; S9 ALLMCL_FFM = (FFM/45)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent on fat-free mass for central volume")                        # Methods, Pharmacokinetic modeling; S9 ALLMV_FFM = (FFM/45)

    # --- Absorption. Savic transit chain feeding a first-order depot.
    lka <- log(2.5)
    label("Log of the first-order absorption rate constant (1/h)")                                # Table 2 First-order absorption rate constant 2.5 (2.29-2.68)
    lmtt <- log(0.291)
    label("Log of the mean absorption transit time (h)")                                          # Table 2 Mean absorption transit time 0.291 (0.184-0.380)
    lntr <- log(4.25)
    label("Log of the estimated number of absorption transit compartments (unitless)")                       # Table 2 Number of absorption transit compartments 4.25 (3.76-4.88)
    lfdepot <- fixed(log(1))
    label("Log of oral bioavailability (unitless, 1; all other parameters are relative to it)")    # Table 2 Bioavailability 1 Fixed; S9 $THETA 4 (1) FIX

    # --- Covariate effect. Clearance is 30.2% higher at the day-28 visit.
    e_day28_cl <- 0.302
    label("Fractional increase in clearance at the day-28 PK visit")                              # Table 2 Change in clearance on day 28 +30.2% (+23.8 to +37.3)

    # --- CSF effect compartment. Table 2 reports the equilibration
    # half-life, so ke0 is recovered as log(2) / HL rather than read off the
    # control stream, whose $THETA block holds initial estimates only.
    lke0 <- log(log(2) / 0.66)
    label("Log of the plasma-to-CSF equilibration rate constant (1/h)")                           # Table 2 Plasma-to-CSF equilibrium half-life 0.66 h (0.43-0.90) -> ke0 = log(2)/0.66 = 1.05 1/h
    lppc <- log(1.05)
    label("Log of the CSF-to-plasma pseudo-partition coefficient (unitless)")                     # Table 2 CSF-to-plasma pseudo-partition coefficient 1.05 (0.99-1.09)

    # --- Random effects. Table 2 reports every variance component as a
    # percentage that is the omega standard deviation on the log scale, not
    # a lognormal CV: the S9 $OMEGA initials reproduce the tabulated
    # percentages as sqrt(omega^2) exactly (0.0249 -> 15.8%, 0.763 -> 87.3%,
    # 1.04 -> 102%), which fixes the convention for this table.
    etalcl ~ 0.0342
    label("Between-subject variability in clearance (log-scale variance)")                        # Table 2 BSV in clearance 18.5% (15.0-23.0); 0.185^2 = 0.0342; S9 $OMEGA 1 0.0341
    # The control stream fixes BSV on V, ka, bioavailability, MTT, ke0 and
    # PPC to zero ($OMEGA 2-7 all "0 FIX"), so those etas are absent here.

    etaiov_fdepot_1 ~ 0.0250
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")     # Table 2 BOV in bioavailability 15.8% (11.4-19.2); 0.158^2 = 0.0250; S9 $OMEGA 13 0.0249
    etaiov_fdepot_2 ~ fixed(0.0250)                                                               # OMEGA 14 equal to OMEGA 13 per S9 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.0250)                                                               # OMEGA 15 equal to OMEGA 13 per S9 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.0250)                                                               # OMEGA 16 equal to OMEGA 13 per S9 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fixed(0.0250)                                                               # OMEGA 17 equal to OMEGA 13 per S9 $OMEGA BLOCK(1) SAME

    etaiov_ka_1 ~ 0.762
    label("Between-occasion variability in absorption rate constant, occasion 1 (log-scale variance)")  # Table 2 BOV in absorption rate constant 87.3% (70.6-103); 0.873^2 = 0.762; S9 $OMEGA 18 0.763
    etaiov_ka_2 ~ fixed(0.762)                                                                    # OMEGA 19 equal to OMEGA 18 per S9 $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.762)                                                                    # OMEGA 20 equal to OMEGA 18 per S9 $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fixed(0.762)                                                                    # OMEGA 21 equal to OMEGA 18 per S9 $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fixed(0.762)                                                                    # OMEGA 22 equal to OMEGA 18 per S9 $OMEGA BLOCK(1) SAME

    etaiov_mtt_1 ~ 1.04
    label("Between-occasion variability in mean absorption transit time, occasion 1 (log-scale variance)")  # Table 2 BOV in mean absorption transit time 102% (73.5-131); 1.02^2 = 1.04; S9 $OMEGA 23 1.04
    etaiov_mtt_2 ~ fixed(1.04)                                                                    # OMEGA 24 equal to OMEGA 23 per S9 $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(1.04)                                                                    # OMEGA 25 equal to OMEGA 23 per S9 $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fixed(1.04)                                                                    # OMEGA 26 equal to OMEGA 23 per S9 $OMEGA BLOCK(1) SAME
    etaiov_mtt_5 ~ fixed(1.04)                                                                    # OMEGA 27 equal to OMEGA 23 per S9 $OMEGA BLOCK(1) SAME

    # --- Residual error, one combined proportional-plus-additive model per
    # matrix. The S9 $ERROR builds each additive term as THETA + 0.2 * LLOQ
    # with both THETAs FIX 0, so the additive standard deviations are
    # exactly 20% of the matrix-specific LLOQ, which is what Table 2
    # footnote c describes. Plasma: 0.2 * 0.200 = 0.04 mg/L, matching the
    # tabulated value. CSF: 0.2 * 0.234 = 0.0468 mg/L; Table 2 prints 0.04
    # for the CSF row too, which cannot be 20% of the 0.234 mg/L CSF LLOQ
    # that the same table's footnote invokes. The rule and the control
    # stream agree, so 0.0468 is used. See the vignette Errata.
    propSd <- 0.0833
    label("Proportional residual error for plasma pyrazinamide (fraction)")                       # Table 2 Proportional error for plasma 8.33% (7.25-9.08)
    addSd <- fixed(0.04)
    label("Additive residual error for plasma pyrazinamide (mg/L)")                               # Table 2 Additive error for plasma 0.04 Fixed = 0.2 * 0.200 mg/L LLOQ (footnote c; S9 LLOQ_P = 0.2)
    propSd_Ccsf <- 0.114
    label("Proportional residual error for CSF pyrazinamide (fraction)")                          # Table 2 Proportional error for CSF 11.4% (7.54-16.3)
    addSd_Ccsf <- fixed(0.0468)
    label("Additive residual error for CSF pyrazinamide (mg/L)")                                  # Footnote c rule = 0.2 * 0.234 mg/L CSF LLOQ (S9 LLOQ_E = 0.234); Table 2 prints 0.04
  })

  model({
    # 1. Occasion indicators. The S9 control stream multiplexes the
    # between-occasion etas with IF (OCC==k) blocks over five occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5

    # 2. Individual parameters. Clearance and volume are allometrically
    # scaled on fat-free mass against the 45 kg cohort median, and clearance
    # additionally takes the day-28 visit step
    # (S9: TVCL = THETA(1)*ALLMCL_FFM*PK_VISIT_CL with PK_VISIT_CL = 1 at
    # day 3 and 1 + THETA(7) at day 28).
    cl <- exp(lcl + etalcl) * (FFM / 45)^e_ffm_cl * (1 + e_day28_cl * DAY28)
    vc <- exp(lvc) * (FFM / 45)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    ke0 <- exp(lke0)
    ppc <- exp(lppc)

    kel <- cl / vc

    # 3. Savic transit-compartment absorption, written out in closed form.
    # The S9 $DES computes the gamma-density input as
    #   KTR     = (NN + 1) / MTT
    #   PIZZA   = LOG(BIO*PD*KTR) - GAMLN(NN+1)
    #   TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO)
    #   DADT(1) = TRANSIT - KA*A(1)
    # with PD the most recent dose amount and TEMPO the time after that
    # dose, guarded by IF (PD > 0 AND TEMPO > 0). podo(depot) and tad(depot)
    # supply PD and TEMPO. The guard is reproduced because at TEMPO = 0 the
    # NN*LOG(0) term is -Inf.
    #
    # This is written out rather than delegated to rxode2's transit()
    # built-in: transit() combined with the f(depot) <- 0 that the control
    # stream's F1 = 0 requires evaluates to an identically zero input rate,
    # so the model would simulate flat zero concentrations. The closed form
    # below is unaffected because podo() and tad() are both still live under
    # f(depot) <- 0.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) +
                    ntr * log(ktt) - ktt)
    }

    # 4. ODE system.
    Cc <- central / vc

    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot - kel * central
    # The CSF state holds a concentration, not an amount (supplementary S2).
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    # The dose amount is delivered entirely through the transit density, so
    # the ordinary bolus into the depot is suppressed (S9 $PK sets F1 = 0).
    f(depot) <- 0

    # 5. Observations and residual error.
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
