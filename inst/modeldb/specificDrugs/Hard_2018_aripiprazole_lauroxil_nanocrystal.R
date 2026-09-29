Hard_2018_aripiprazole_lauroxil_nanocrystal <- function() {
  description <- "Two-compartment population PK model for aripiprazole with three parallel inputs -- a double-Weibull dissolution input from the intramuscular nanocrystal dispersion of aripiprazole lauroxil (ALNCD), a lagged zero-order intramuscular input from standard aripiprazole lauroxil, and first-order oral absorption -- in adults with schizophrenia or schizoaffective disorder"
  reference <- paste(
    "Hard ML, Wehr AY, Sadler BM, Mills RJ, von Moltke L (2018).",
    "Population Pharmacokinetic Analysis and Model-Based Simulations of Aripiprazole",
    "for a 1-Day Initiation Regimen for the Long-Acting Antipsychotic Aripiprazole Lauroxil.",
    "Eur J Drug Metab Pharmacokinet 43(4):461-469. doi:10.1007/s13318-018-0488-4.",
    "All parameter estimates are from Supplementary Table S3 (Online Resource 1);",
    "the model structure is from Supplementary Fig. S1 (Online Resource 3) and the",
    "Supplementary Text (Online Resource 2).",
    "This model supersedes the earlier aripiprazole lauroxil PopPK model of Hard ML,",
    "Mills RJ, Sadler BM, Wehr AY, Weiden PJ, von Moltke L (2017) CNS Drugs 31(7):617-624,",
    "doi:10.1007/s40263-017-0447-7; see modellib('Hard_2017_aripiprazole_lauroxil').",
    sep = " "
  )
  vignette <- "Hard_2018_aripiprazole_lauroxil_nanocrystal"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # buildModelDb()'s dosing heuristic recognises only states literally named
  # 'depot' and 'central', so without this field the registry would record
  # 'depot,central' and silently omit the two nanocrystal-dispersion Weibull
  # depots, which are genuine dose targets.
  dosing <- c("depot", "depot2", "depot3", "central")

  # Supplementary Table S3 flags with '*' the parameters NONMEM estimated on the
  # log scale; those are written here as log()-wrapped thetas carrying an
  # exponential eta. GAM1, GAM2 and FRAC carry no '*' and are therefore written
  # on their natural scale, with the exponential inter-individual variability
  # applied multiplicatively -- Table S3 footnote 'e' confirms those rows are
  # log-normal, because its CV column is reproduced by sqrt(exp(omega^2) - 1)
  # (e.g. GAM1 sqrt(exp(0.461) - 1) = 0.765 -> the tabulated 76.5%).
  paper_specific_residual_sds <- c("propSdStdyA105", "propSdStdyAlncd")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Subject-level. The single covariate effect surviving backward elimination:",
        "an allometric power effect on the central apparent volume Vc/F with the exponent",
        "held at 1.0 and the weight centred at 70 kg (Supplementary Table S3 row",
        "'Weight ON VC/F' = 1, footnote 'a' 'Fixed at 1.00' and footnote 'f'",
        "'Power effect = VC/F*(weight/70)1.0'; the 70 kg centring is stated in Sect. 2.1.2",
        "of the paper, 'fixed allometric exponents of 0.75 and 1, respectively, and scaled",
        "to 70 kg'). Weight on CL/F was carried in the full model at the allometric 0.75",
        "(Supplementary Table S2) but was removed in backward elimination, so clearance in",
        "this final model does NOT scale with weight (Online Resource 2, 'Covariate model",
        "development'). Estimating the Vc/F exponent rather than fixing it returned",
        "1.15 (95% CI 0.940, 1.36), whose interval contains 1.0, so the fixed value was kept.",
        "Cohort mean weight was 89.1 kg (SD 17.9) across the four studies (Table 1)."
      ),
      source_name = "WT"
    ),
    STUDY_A105 = list(
      description = "Phase I study ALK9072-A105 membership indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the three nanocrystal-dispersion studies ALK9072-B101, B102 and B103)",
      notes = paste(
        "Subject-level. Selects the proportional residual-error magnitude. A105 (Study 4 of",
        "this paper) is the earlier phase I study in which aripiprazole lauroxil was given",
        "alone; adding a separate proportional error term for it dropped the objective",
        "function by 411 points and is reported as a lower residual variability of 14.4%",
        "versus 18.9% for the three nanocrystal-dispersion studies (Online Resource 2,",
        "'Final PopPK model update'). Supplementary Table S3 'Residual variability' block:",
        "sigma^2 prop Studies 1, 2, and 3 = 0.0359; sigma^2 prop Study 4 = 0.0207.",
        "Same study, and the same covariate name, as in",
        "modellib('Hard_2017_aripiprazole_lauroxil')."
      ),
      source_name = "Study"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "aripiprazole", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "aripiprazole lauroxil", units = "mg", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "aripiprazole lauroxil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 343,
    n_studies = 4,
    age_range = "adults",
    age_median = "mean 45.2 years (SD 10.8)",
    weight_range = "not reported",
    weight_median = "mean 89.1 kg (SD 17.9)",
    sex_female_pct = 27,
    race_ethnicity = c(Black = 78, White = 21, Other = 1),
    disease_state = "schizophrenia or schizoaffective disorder, stable on a first-line antipsychotic other than aripiprazole",
    dose_range = paste(
      "aripiprazole lauroxil 441-1064 mg IM and nanocrystal dispersion 110-662 mg IM",
      "(modelled as aripiprazole equivalents of 75, 150, 300, 450, 600 and 724 mg for",
      "110, 221, 441, 662, 882 and 1064 mg respectively), plus oral aripiprazole 15 mg",
      "once daily or a single 30 mg dose"
    ),
    regions = "USA",
    cyp2d6 = "extensive and intermediate metabolizers plus inconclusive phenotypes; poor metabolizers were absent (excluded from studies 1-3), so the model carries no CYP2D6 term and does not apply to poor metabolizers",
    notes = paste(
      "12,768 plasma aripiprazole concentrations (351 [3%] below the lower limit of",
      "quantitation of 1 ng/mL, handled by the M3 method) from 343 patients in four phase I",
      "studies: ALK9072-B101 (Study 1, n = 41, single-dose ascending ALNCD, gluteal),",
      "ALK9072-B102 (Study 2, n = 161, ALNCD + 30 mg oral aripiprazole + AL 441 or 882 mg,",
      "or the 21-day oral initiation regimen), ALK9072-B103 (Study 3, n = 47, ALNCD deltoid",
      "versus gluteal) and ALK9072-A105 (Study 4, n = 94, AL alone q4wk/q6wk/q8wk).",
      "The data set held 2536 dosing records (1742 oral aripiprazole, 626 AL, 168 ALNCD).",
      "Demographics are Table 1 of the paper; ethnicity was 5% Hispanic or Latino.",
      "Only 3 patients (1%) were of a race other than Black or African American or White,",
      "and race and sex were not carried into the covariate analysis."
    )
  )

  ini({
    # ---- Disposition (Supplementary Table S3) --------------------------------
    lcl <- log(1.98)
    label("Apparent clearance of aripiprazole (CL/F, L/h)")
    # Supplementary Table S3 row 'CL/F (L/h)' = 1.98 (%RSE 2.54; 95% CI 1.88, 2.07)

    lvc <- log(327)
    label("Apparent central volume of distribution at 70 kg (VC/F, L)")
    # Supplementary Table S3 row 'VC/F (L)' = 327 (%RSE 4.51; 95% CI 298, 356)

    lvp <- log(1720)
    label("Apparent peripheral volume of distribution (VP/F, L)")
    # Supplementary Table S3 row 'VP/F (L)' = 1720 (%RSE 13.9; 95% CI 1251, 2188)

    lq <- log(0.102)
    label("Apparent inter-compartmental clearance (Q/F, L/h)")
    # Supplementary Table S3 row 'Q/F (L/h)' = 0.102 (%RSE 10.9; 95% CI 0.080, 0.124)

    # ---- Oral aripiprazole input ---------------------------------------------
    lka <- log(0.47)
    label("First-order absorption rate constant for oral aripiprazole (1/h)")
    # Supplementary Table S3 row 'Ka PO ARI (h-1)' = 0.47 (%RSE 14.2; 95% CI 0.339, 0.601)

    lfdepot <- fixed(log(1))
    label("Bioavailability of oral aripiprazole, the reference route (FPO, unitless)")
    # Supplementary Table S3 row 'FPO ARI' = 1, footnote 'a' 'Fixed at 1.00'. Online
    # Resource 2 adds that its inter-individual variability was also fixed to zero,
    # so no eta is carried on this parameter.

    # ---- Aripiprazole lauroxil (AL) intramuscular input -----------------------
    ld1 <- log(934)
    label("Duration of the zero-order aripiprazole input after an intramuscular aripiprazole lauroxil injection (D AL, h)")
    # Supplementary Table S3 row 'D AL (h)' = 934 (%RSE 3.86; 95% CI 864, 1005);
    # 934 h = 38.9 days.

    ltlag <- log(106)
    label("Lag time before aripiprazole appears in the central compartment after an intramuscular aripiprazole lauroxil injection (ALAG AL, h)")
    # Supplementary Table S3 row 'ALAG AL (h)' = 106 (%RSE 7.85; 95% CI 89.4, 122);
    # 106 h = 4.4 days.

    lfdepot_im <- fixed(log(0.571))
    label("Bioavailability of intramuscular aripiprazole lauroxil relative to oral aripiprazole (FIM AL, unitless)")
    # Supplementary Table S3 row 'FIM AL' = 0.571, footnote 'b' 'Fixed at 57.1% from
    # previous analysis'; that analysis is Hard 2017 CNS Drugs, whose Supplemental
    # Table 7 estimated FIM = 0.571 (95% CI 0.542, 0.599).

    # ---- ALNCD intramuscular input: double Weibull dissolution ----------------
    lfdepot_ncd <- log(1.12)
    label("Bioavailability of the intramuscular nanocrystal dispersion relative to intramuscular aripiprazole lauroxil (FIM ALNCD, unitless)")
    # Supplementary Table S3 row 'FIM ALNCD' = 1.12 (%RSE 3.56; 95% CI 1.04, 1.20),
    # footnote 'c' 'ALNCD F estimated relative to AL'. Multiplying by FIM AL gives the
    # bioavailability relative to oral that the table reports on the next row as 0.638
    # (footnote 'd'); 1.12 * 0.571 = 0.6395, and the table notes 'more decimal places
    # used in calculation than presented'.

    lwa1 <- log(596)
    label("Scale of the slow-dissolving Weibull for the nanocrystal dispersion (MDT1, h)")
    # Supplementary Table S3 row 'MDT1 (h)' = 596 (%RSE 3.29; 95% CI 557, 634);
    # 596 h = 24.8 days.

    wb1 <- 2.2
    label("Shape of the slow-dissolving Weibull for the nanocrystal dispersion (GAM1, unitless)")
    # Supplementary Table S3 row 'GAM1' = 2.2 (%RSE 2.85; 95% CI 2.08, 2.32)

    lwa2 <- log(76.7)
    label("Scale of the fast-dissolving Weibull for the nanocrystal dispersion (MDT2, h)")
    # Supplementary Table S3 row 'MDT2 (h)' = 76.7 (%RSE 4.47; 95% CI 70.0, 83.4);
    # 76.7 h = 3.2 days.

    wb2 <- 2.09
    label("Shape of the fast-dissolving Weibull for the nanocrystal dispersion (GAM2, unitless)")
    # Supplementary Table S3 row 'GAM2' = 2.09 (%RSE 2.36; 95% CI 1.99, 2.19)

    logitfrac <- 2.02
    label("Logit of the nanocrystal-dispersion dose fraction routed through the slow-dissolving Weibull (FRAC, unitless)")
    # Supplementary Table S3 row 'FRAC' = 2.02 (%RSE 3.88; 95% CI 1.87, 2.17).
    # The paper never writes the double-Weibull equation, and 2.02 cannot be a bare
    # fraction, so the transform is inferred: expit(2.02) = 0.883 of the dose through
    # the MDT1 Weibull. See the vignette 'Assumptions and deviations' section -- this
    # is the one structural inference in the file and it is verified against the
    # independently reported ALNCD-alone pharmacokinetics of study ALK9072-B103.

    # ---- Baseline -------------------------------------------------------------
    lc0 <- log(0.378)
    label("Pre-first-dose aripiprazole plasma concentration (ARI(0), ng/mL)")
    # Supplementary Table S3 row 'ARI(0) (ng/mL)' = 0.378 (%RSE 14.8; 95% CI 0.269,
    # 0.488). Well below the 1 ng/mL lower limit of quantitation, reflecting the few
    # patients with quantifiable pre-dose aripiprazole (Online Resource 2).

    # ---- Covariate effect ------------------------------------------------------
    e_wt_vc <- fixed(1)
    label("Allometric exponent on body weight for VC/F, centred at 70 kg (unitless)")
    # Supplementary Table S3 row 'Weight ON VC/F' = 1, footnote 'a' 'Fixed at 1.00'
    # and footnote 'f' 'Power effect = VC/F*(weight/70)1.0'.

    # ---- Inter-individual variability -----------------------------------------
    # Supplementary Table S3's 'Interindividual variability' Point Estimate column
    # holds the omega^2 VARIANCE. Confirmed against the table's own CV% column and
    # its footnote 'e': the flagged rows satisfy CV = sqrt(exp(omega^2) - 1) * 100
    # -- CL/F sqrt(exp(0.539) - 1) = 0.845 -> 84.5%, Q/F sqrt(exp(2.6) - 1) = 3.53
    # -> 353%. The two unflagged rows (GAM2, D AL) are the small-variance
    # approximation sqrt(omega^2): sqrt(0.129) = 0.359 -> 35.9%.
    # The paper states these 13 terms sat in a full Omega block but publishes only
    # the diagonal ('Only diagonal elements of the full Omega block are presented'),
    # so they are carried here as independent etas; see the vignette Errata.
    etalcl ~ 0.539 # Supplementary Table S3 CL/F IIV variance 0.539 (95% CI 0.425, 0.653)
    etalvc ~ 0.239 # Supplementary Table S3 VC/F IIV variance 0.239 (95% CI 0.160, 0.318)
    etalwa1 ~ 0.197 # Supplementary Table S3 MDT1 IIV variance 0.197 (95% CI 0.147, 0.247)
    etawb1 ~ 0.461 # Supplementary Table S3 GAM1 IIV variance 0.461 (95% CI 0.273, 0.649)
    etalwa2 ~ 0.274 # Supplementary Table S3 MDT2 IIV variance 0.274 (95% CI 0.193, 0.355)
    etawb2 ~ 0.129 # Supplementary Table S3 GAM2 IIV variance 0.129 (95% CI 0.0684, 0.190)
    etalogitfrac ~ 0.978 # Supplementary Table S3 FRAC IIV variance 0.978 (95% CI 0.635, 1.32)
    etalvp ~ 0.686 # Supplementary Table S3 VP/F IIV variance 0.686 (95% CI 0.435, 0.937)
    etalq ~ 2.6 # Supplementary Table S3 Q/F IIV variance 2.6 (95% CI 2.05, 3.15)
    etalka ~ 2.45 # Supplementary Table S3 Ka PO ARI IIV variance 2.45 (95% CI 1.67, 3.23)
    etald1 ~ 0.147 # Supplementary Table S3 D AL IIV variance 0.147 (95% CI 0.0845, 0.210)
    etaltlag ~ 0.895 # Supplementary Table S3 ALAG AL IIV variance 0.895 (95% CI 0.630, 1.16)
    etalfdepot_im ~ 0.505 # Supplementary Table S3 FIM AL IIV variance 0.505 (95% CI 0.399, 0.611); estimated even though FIM AL itself is fixed
    etalc0 ~ 4.34 # Supplementary Table S3 Ari(0) IIV variance 4.34 (95% CI 3.44, 5.24); estimated outside the Omega block

    # ---- Residual error --------------------------------------------------------
    propSdStdyAlncd <- 0.18947
    label("Proportional residual SD for the three nanocrystal-dispersion studies (fraction)")
    # Supplementary Table S3 'sigma2 prop Studies 1, 2, and 3' = 0.0359 (%RSE 3.93;
    # 95% CI 0.0331, 0.0387); sqrt(0.0359) = 0.18947, matching its 18.9% CV column.

    propSdStdyA105 <- 0.14387
    label("Proportional residual SD for phase I study ALK9072-A105 (fraction)")
    # Supplementary Table S3 'sigma2 prop Study 4' = 0.0207 (%RSE 2.70; 95% CI 0.0196,
    # 0.0218); sqrt(0.0207) = 0.14387, matching its 14.4% CV column.
  })

  model({
    # ---- 1. Individual parameters ---------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)
    c0 <- exp(lc0 + etalc0)
    fdepot <- exp(lfdepot)
    fdepot_im <- exp(lfdepot_im + etalfdepot_im)
    fdepot_ncd <- exp(lfdepot_ncd)
    wa1 <- exp(lwa1 + etalwa1)
    wa2 <- exp(lwa2 + etalwa2)
    # GAM1, GAM2 and FRAC are tabulated on their natural scale (no '*' in
    # Supplementary Table S3) but carry log-normal inter-individual variability,
    # so the eta is applied multiplicatively rather than inside an exp(l... + eta).
    wb1i <- wb1 * exp(etawb1)
    wb2i <- wb2 * exp(etawb2)
    fracw <- expit(logitfrac * exp(etalogitfrac))

    # ---- 2. Micro-constants -----------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. Weibull dissolution hazards -----------------------------------------
    # The nanocrystal dispersion releases aripiprazole by a double Weibull
    # (Supplementary Fig. S1, 'Described by double Weibull function'). The fraction
    # released by time t after the injection is
    #     FR(t) = fracw * (1 - exp(-(t/wa1)^wb1)) + (1 - fracw) * (1 - exp(-(t/wa2)^wb2))
    # which is encoded as two depots emptying at the corresponding Weibull hazards,
    # the same parameterisation nlmixr2lib's own addWeibullAbs() uses.
    #
    # The hazard of a Weibull with shape > 1 grows without bound in time, so long
    # after the depot is numerically empty lsoda still sees an arbitrarily fast
    # rate acting on an amount of order 1e-300 and reports 'h too small for
    # machine precision'. The hazard is therefore capped at 50/wa. That bound
    # cannot change any result: the Weibull survivor exp(-(t/wa)^wb) is already
    # below 1e-12 by the time the uncapped hazard reaches ~13/wa for these shape
    # values, so the cap only ever engages after the depot has released more than
    # 99.9999999999% of its dose.
    h1 <- min((wb1i / wa1) * (tad0(depot2) / wa1)^(wb1i - 1), 50 / wa1)
    h2 <- min((wb2i / wa2) * (tad0(depot3) / wa2)^(wb2i - 1), 50 / wa2)

    # ---- 4. ODE system ----------------------------------------------------------
    # Supplementary Fig. S1: three inputs feed one aripiprazole central compartment
    # that exchanges with one peripheral compartment and is cleared by CL/F.
    #   depot       - oral aripiprazole, first-order at ka
    #   depot2      - the slow-dissolving Weibull share of a nanocrystal injection
    #   depot3      - the fast-dissolving Weibull share of a nanocrystal injection
    #   central     - aripiprazole lauroxil enters here directly as a lagged
    #                 zero-order input of duration D AL. The paper describes seven
    #                 IM depots ('the model was expanded by adding an additional 6 IM
    #                 dosing depot for AL'), but they are NONMEM book-keeping for
    #                 overlapping injections and share one D AL, one ALAG and one
    #                 FIM AL, so they superpose exactly onto a single lagged
    #                 modelled-duration input. AL dose records therefore use
    #                 cmt = "central" with rate = -2, as in
    #                 modellib('Hard_2017_aripiprazole_lauroxil').
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -h1 * depot2
    d/dt(depot3) <- -h2 * depot3
    d/dt(central) <- ka * depot + h1 * depot2 + h2 * depot3 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Pre-dose aripiprazole carried in as an initial condition. c0 is ng/mL and vc
    # is L, so ng/mL * L = ug and the /1000 converts to the mg the states hold.
    central(0) <- c0 * vc / 1000

    # ---- 5. Route-specific input -------------------------------------------------
    # A nanocrystal-dispersion injection is entered as two dose records carrying the
    # same amt, one into depot2 and one into depot3; the f() factors apply the
    # bioavailability and split the dose between the two Weibull pathways.
    f(depot) <- fdepot
    f(depot2) <- fdepot_im * fdepot_ncd * fracw
    f(depot3) <- fdepot_im * fdepot_ncd * (1 - fracw)
    f(central) <- fdepot_im
    dur(central) <- d1
    alag(central) <- tlag

    # ---- 6. Observation and error -------------------------------------------------
    # central is in mg and vc in L, so central/vc is ug/mL; *1000 gives ng/mL.
    Cc <- 1000 * central / vc
    propSdStudy <- propSdStdyA105 * STUDY_A105 + propSdStdyAlncd * (1 - STUDY_A105)
    Cc ~ prop(propSdStudy)
  })
}
