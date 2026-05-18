Jin_2014_simvastatin <- function() {
  description <- "Joint two-compartment population PK model for orally administered simvastatin (lactone parent) and its active metabolite simvastatin acid (open beta-hydroxyacid), describing atypical multiple-peak absorption via three parallel mixed zero-and-first-order absorption processes and non-equilibrium reversible interconversion between the two species (Jin 2014). The simvastatin lactone is delivered into three depot compartments with fractional bioavailabilities F1, F2, F3 (sum = 1) parameterised through two relative-bioavailability constants BA1 and BA2, each depot has its own first-order absorption rate constant Ka1, Ka2, Ka3, zero-order infusion duration D1, D2, D3, and absorption lag-times ALAG1 = 0, ALAG2, ALAG3, the lactone disposes via a 2-compartment system with apparent clearance CL and inter-compartmental clearance Q, the fraction FM of total parent CL is converted to simvastatin acid (V_acid central fixed at 1 L for identifiability), and a reverse clearance Q64 returns acid to the parent central compartment. Age, body weight, and height were tested as covariates and not retained in the final model. The source publication analysed data in molar units; this packaged model preserves that choice -- doses are expressed in nmol and concentrations in nmol/L. The validation vignette demonstrates the standard milligram-to-nanomole conversion using the simvastatin lactone molecular weight."
  reference <- paste(
    "Jin SJ, Bae KS, Cho SH, Jung JA, Kim U, Choe S, Ghim JL,",
    "Noh YH, Park HJ, Kim Hs, Lim HS.",
    "Population Pharmacokinetic Analysis of Simvastatin and its Active",
    "Metabolite with the Characterization of Atypical Complex Absorption Kinetics.",
    "Pharm Res. 2014;31(7):1801-1812.",
    "doi:10.1007/s11095-013-1284-0.",
    sep = " "
  )
  vignette <- "Jin_2014_simvastatin"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    # Jin 2014 tested AGE, body weight (WT), and height as candidate
    # covariates on both the parent (simvastatin lactone) and the
    # metabolite (simvastatin acid) parameters; none reached the
    # significance criterion and the final model carries no covariate
    # effects (Results, parent and metabolite sections). The covariate
    # list is therefore empty.
  )

  population <- list(
    species        = "human",
    n_subjects     = 133L,
    n_studies      = 4L,
    n_observations = "2,182 simvastatin and 2,130 simvastatin acid plasma concentrations (Table I).",
    age_range      = "approx 22-30 years (mean +/- SD = 25.0 +/- 2.6)",
    weight_range   = "approx 49-91 kg (mean +/- SD = 69.8 +/- 6.7)",
    height_range   = "approx 160-185 cm (mean +/- SD = 174.8 +/- 5.4)",
    bmi_range      = "approx 19-27 kg/m^2 (mean +/- SD = 22.8 +/- 1.7)",
    sex_female_pct = 0,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adult Korean male volunteers across one PK crossover (Study I, 20 subjects, randomised 2-of-5 dose combinations from 10, 20, 30, 40, 80 mg simvastatin separated by 3-day washout) and three single-dose 60 mg bioequivalence studies (Studies II, III, IV, 30+30+53 = 113 subjects, 2x2 crossover with reference and test formulations, 7-day washout, reference Zocor only used for this analysis).",
    dose_range     = "Single oral doses of 10, 20, 30, 40, 60, or 80 mg simvastatin (Zocor, Merck, NJ, USA) after overnight fast with 240 mL water; standard meals 4 and 9 h post-dose. The crossover study (Study I) used dense PK sampling at pre-dose and 10, 20, 30, 40, 50 min and 1, 1.25, 1.5, 1.75, 2, 3, 4, 6, 8, 10, 12, 24 h post-dose; the bioequivalence studies (II, III, IV) sampled at pre-dose and 15, 30, 45 min and 1, 1.5, 2, 3, 4, 6, 8, 10, 12 h post-dose.",
    regions        = "Republic of Korea (Asan Medical Center Clinical Trial Center, Seoul).",
    bloq_handling  = "Concentrations below the lower limit of quantification (LLOQ = 0.2 ng/mL simvastatin, 0.1 ng/mL simvastatin acid) were imputed at half the LLOQ for the first sub-LLOQ value in the disposition phase; subsequent sub-LLOQ values were removed.",
    notes          = "Baseline demographics reproduced from Jin 2014 Table I (subject characteristics by study). No female subjects; all subjects of Korean ancestry. The final model carries no significant demographic covariates (Results, parent and metabolite sections)."
  )

  ini({
    # Notes on parameter values:
    # * Point estimates and CIs come from Jin 2014 Table III ('Final
    #   pharmacokinetic model' columns) -- the canonical final-model
    #   estimates. The Results narrative paragraph reports BA1 = 0.16
    #   and BA2 = 0.70 (giving F1, F2, F3 = 53.8 %, 8.6 %, 37.6 %),
    #   but those values are inconsistent with the Table III CI of
    #   0.02 to 1.25 reported for BA1 (which requires BA1 ~ 0.636
    #   with the listed RSE of 49.4 %). Table III is treated as
    #   authoritative; the narrative percentages are noted as a
    #   transcription error in the vignette Errata section.
    # * Concentrations and amounts are in molar units (nmol/L, nmol)
    #   per the source ('the data were analyzed after changing the
    #   units of amount and concentration into molar units'). Oral
    #   dose in mg is converted to nmol inside model() via the
    #   simvastatin (lactone) molecular weight.
    # * %CV inter-individual variability values from Table III are
    #   translated to log-scale variance via omega^2 = log(1 + CV^2);
    #   each line below carries the source CV in a comment.
    # * Inter-occasion variability (IOV) was estimated on BA1, BA2,
    #   D1, D2, D3 (Table III) but is not encoded as a separate
    #   variance component in this packaged model; see the vignette
    #   Assumptions and deviations section for the rationale.

    # Parent disposition (simvastatin lactone) -- Table III
    lcl  <- log(571);    label("Apparent simvastatin clearance CL/F (L/h)")                   # Jin 2014 Table III: CL = 571 L/h (RSE 17.1 %)
    lvc  <- log(199);    label("Apparent simvastatin central volume V4/F (L)")                # Jin 2014 Table III: Vcentral = 199 L (RSE 19.1 %)
    lvp  <- log(2710);   label("Apparent simvastatin peripheral volume V5/F (L)")             # Jin 2014 Table III: Vperipheral = 2710 L (RSE 17.3 %)
    lq   <- log(199);    label("Apparent simvastatin inter-compartmental clearance Q (L/h)")  # Jin 2014 Table III: Q = 199 L/h (RSE 23.9 %)

    # Parent absorption -- Table III
    lka1   <- log(0.00126); label("Absorption rate constant from depot 1 Ka1 (1/h)")          # Jin 2014 Table III: Ka1 = 0.00126 1/h (RSE 54.6 %)
    lka2   <- log(0.964);   label("Absorption rate constant from depot 2 Ka2 (1/h)")          # Jin 2014 Table III: Ka2 = 0.964 1/h (RSE 0.3 %)
    lka3   <- log(0.179);   label("Absorption rate constant from depot 3 Ka3 (1/h)")          # Jin 2014 Table III: Ka3 = 0.179 1/h (RSE 8.2 %)
    lba1   <- log(0.636);   label("Relative bioavailability parameter BA1 (unitless)")        # Jin 2014 Table III: BA1 = 0.636 (RSE 49.4 %); F1 = 1/(1+BA1+BA2)
    lba2   <- log(0.662);   label("Relative bioavailability parameter BA2 (unitless)")        # Jin 2014 Table III: BA2 = 0.662 (RSE 13.2 %); F2 = BA1/(1+BA1+BA2), F3 = BA2/(1+BA1+BA2)
    llag2  <- log(0.142);   label("Absorption lag time of depot 2 ALAG2 (h)")                 # Jin 2014 Table III: ALAG2 = 0.142 h (RSE 6.9 %)
    llag3d <- log(0.787);   label("Difference between ALAG3 and ALAG2 (h); ALAG3 = ALAG2 + this")  # Jin 2014 Table III: ALAG3 - ALAG2 = 0.787 h (RSE 8.4 %)
    ldur1  <- log(0.102);   label("Zero-order absorption duration of depot 1 D1 (h)")         # Jin 2014 Table III: D1 = 0.102 h (RSE 0.06 %)
    ldur2  <- log(0.502);   label("Zero-order absorption duration of depot 2 D2 (h)")         # Jin 2014 Table III: D2 = 0.502 h (RSE 9.1 %)
    ldur3  <- log(1.38);    label("Zero-order absorption duration of depot 3 D3 (h)")         # Jin 2014 Table III: D3 = 1.38 h (RSE 5.54 %)

    # Simvastatin acid metabolite -- Table III
    logitfm  <- qlogis(0.236); label("Logit of fraction of simvastatin CL forming simvastatin acid FM (unitless)")  # Jin 2014 Table III: FM = 0.236 (RSE 0.4 %)
    lq64     <- log(112);      label("Reverse-conversion clearance from simvastatin acid central to simvastatin central Q64 (L/h)")  # Jin 2014 Table III: Q64 = 112 L/h (RSE 0.2 %)
    lk67     <- log(252);      label("Rate constant from simvastatin acid central to peripheral K67 (1/h)")          # Jin 2014 Table III: K67 = 252 1/h (RSE 6.0 %)
    lk76     <- log(2.30);     label("Rate constant from simvastatin acid peripheral to central K76 (1/h)")          # Jin 2014 Table III: K76 = 2.30 1/h (RSE 9.0 %)
    lcl_acid <- log(0.035);    label("Apparent simvastatin acid clearance CLm/F (L/h); reported with V6_acid fixed at 1 L")  # Jin 2014 Table III: CLm = 0.035 L/h (RSE 41.4 %)

    # Inter-individual variability (lognormal). omega^2 = log(1 + CV^2)
    # from each Table III percent-CV entry. The parent V_central, CL,
    # and Q form a 3x3 correlated block (Table III footnote a: CORR
    # Vcentral-CL = 0.62, Vcentral-Q = -0.17, CL-Q = -0.38). The
    # metabolite K67 and K76 form a 2x2 correlated block (CORR
    # K67-K76 = -0.90).
    etalvc + etalcl + etalq ~ c(0.4849,
                                0.2613,  0.3664,
                               -0.0939, -0.1825, 0.6293)                            # Jin 2014 Table III: IIV %CV V_central 79.0 / CL 66.5 / Q 93.6; CORR V_central-CL 0.62, V_central-Q -0.17, CL-Q -0.38
    etalka1   ~ 0.2003                                                              # Jin 2014 Table III: IIV Ka1 = 47.1 % CV -> log(1 + 0.471^2)
    etalka2   ~ 0.1763                                                              # Jin 2014 Table III: IIV Ka2 = 43.9 % CV -> log(1 + 0.439^2)
    etalka3   ~ 0.2231                                                              # Jin 2014 Table III: IIV Ka3 = 50.0 % CV -> log(1 + 0.500^2)
    etalba1   ~ 0.2043                                                              # Jin 2014 Table III: IIV BA1 = 47.6 % CV -> log(1 + 0.476^2)
    etalba2   ~ 0.1874                                                              # Jin 2014 Table III: IIV BA2 = 45.4 % CV -> log(1 + 0.454^2)
    etallag2  ~ 0.1949                                                              # Jin 2014 Table III: IIV ALAG2 = 46.4 % CV -> log(1 + 0.464^2)
    etallag3d ~ 0.1858                                                              # Jin 2014 Table III: IIV ALAG3-ALAG2 = 45.2 % CV -> log(1 + 0.452^2)
    etaldur1  ~ 0.2003                                                              # Jin 2014 Table III: IIV D1 = 47.1 % CV -> log(1 + 0.471^2)
    etaldur2  ~ 0.1942                                                              # Jin 2014 Table III: IIV D2 = 46.3 % CV -> log(1 + 0.463^2)
    etaldur3  ~ 0.1927                                                              # Jin 2014 Table III: IIV D3 = 46.1 % CV -> log(1 + 0.461^2)
    etalq64   ~ 0.6582                                                              # Jin 2014 Table III: IIV Q64 = 96.5 % CV -> log(1 + 0.965^2)
    etalk67 + etalk76 ~ c(0.5795,
                         -0.6822, 0.9914)                                           # Jin 2014 Table III: IIV K67 = 88.6 % CV, K76 = 130.2 % CV, CORR K67-K76 = -0.90
    etalcl_acid ~ 0.1554                                                            # Jin 2014 Table III: IIV CLm = 41.0 % CV -> log(1 + 0.410^2)

    # Residual error -- Jin 2014 Table III ('Combined additive and
    # proportional model' per Methods). The additive component for both
    # parent and metabolite is reported on the molar concentration
    # scale (nmole/L); the proportional component is a unitless CV.
    addSd       <- 0.118; label("Additive residual error for simvastatin (nmol/L)")              # Jin 2014 Table III: epsilon1 (additive) = 0.118 nmole/L (RSE 6.0 %)
    propSd      <- 0.236; label("Proportional residual error for simvastatin (fraction)")        # Jin 2014 Table III: epsilon2 (proportional) = 0.236 (RSE 2.8 %)
    addSd_acid  <- 0.27;  label("Additive residual error for simvastatin acid (nmol/L)")         # Jin 2014 Table III: epsilon3 (additive) = 0.27 nmole/L (RSE 0.7 %)
    propSd_acid <- 0.43;  label("Proportional residual error for simvastatin acid (fraction)")   # Jin 2014 Table III: epsilon4 (proportional) = 0.43 (RSE 1.0 %)
  })

  model({
    # Volume of the simvastatin-acid central compartment is fixed at 1 L
    # per Jin 2014 Methods ('the central volume of distribution Vd of
    # simvastatin acid was fixed at 1 L, because of the identifiability
    # problem'). The model is faithful to that convention; FM and CLm
    # are interpretable only at V6_acid = 1 L.
    vc_acid <- 1

    # Individual parameters (lognormal IIV; correlated block for the
    # parent V_central + CL + Q is set in ini()).
    ka1  <- exp(lka1  + etalka1)
    ka2  <- exp(lka2  + etalka2)
    ka3  <- exp(lka3  + etalka3)
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq    + etalq)
    vp   <- exp(lvp)
    ba1  <- exp(lba1  + etalba1)
    ba2  <- exp(lba2  + etalba2)
    alag2  <- exp(llag2  + etallag2)
    alag3d <- exp(llag3d + etallag3d)
    alag3  <- alag2 + alag3d
    d1   <- exp(ldur1 + etaldur1)
    d2   <- exp(ldur2 + etaldur2)
    d3   <- exp(ldur3 + etaldur3)

    fm   <- expit(logitfm)
    q64  <- exp(lq64 + etalq64)
    k67  <- exp(lk67 + etalk67)
    k76  <- exp(lk76 + etalk76)
    cl_acid <- exp(lcl_acid + etalcl_acid)

    # Bioavailability fractions for the three parallel parent depots.
    # F1 = 1 / (1 + BA1 + BA2), F2 = BA1 / (1 + BA1 + BA2),
    # F3 = BA2 / (1 + BA1 + BA2). F1 + F2 + F3 = 1 by construction.
    bsum <- 1 + ba1 + ba2
    f1   <- 1   / bsum
    f2   <- ba1 / bsum
    f3   <- ba2 / bsum

    # Inter-compartmental flow rates (parent disposition).
    kel       <- cl / vc
    k_c_to_p1 <- q  / vc
    k_p1_to_c <- q  / vp

    # Conversion rate from acid central back to parent central. Paper
    # parameterises this leg as a clearance Q64 (L/h) acting on the
    # acid central concentration A6/V6_acid (with V6_acid = 1 L); the
    # equivalent rate constant on A6 is Q64 / V6_acid = q64.
    k_acid_to_parent <- q64 / vc_acid

    # ODE system. Each depot is dosed with a mixed zero-and-first-order
    # absorption: dur(depot_k) sets the zero-order infusion duration and
    # the first-order leg follows with rate constant ka_k. ALAG1 = 0
    # (depot 1 has no lag); ALAG2 and ALAG3 = ALAG2 + ALAG3-ALAG2 are
    # imposed via lag(). Amounts are tracked in nmol throughout; the
    # vignette converts oral milligram doses to nmol via the
    # simvastatin lactone molecular weight (418.57 g/mol) before
    # entering the event table.
    d/dt(depot1)            <- -ka1 * depot1
    d/dt(depot2)            <- -ka2 * depot2
    d/dt(depot3)            <- -ka3 * depot3
    d/dt(central)           <- ka1 * depot1 + ka2 * depot2 + ka3 * depot3 -
                               kel * central -
                               k_c_to_p1 * central + k_p1_to_c * peripheral1 +
                               k_acid_to_parent * central_acid
    d/dt(peripheral1)       <- k_c_to_p1 * central - k_p1_to_c * peripheral1
    d/dt(central_acid)      <- fm * kel * central -
                               k_acid_to_parent * central_acid -
                               (cl_acid / vc_acid) * central_acid -
                               k67 * central_acid + k76 * peripheral1_acid
    d/dt(peripheral1_acid)  <- k67 * central_acid - k76 * peripheral1_acid

    # Per-depot bioavailability fractions, zero-order infusion durations,
    # and absorption lag times. A single oral dose event in the data
    # records targets each depot in turn (cmt = depot1, depot2, depot3)
    # with the same amt; the f() multipliers split the dose into
    # F1, F2, F3 (sum = 1).
    f(depot1) <- f1
    f(depot2) <- f2
    f(depot3) <- f3
    dur(depot1) <- d1
    dur(depot2) <- d2
    dur(depot3) <- d3
    lag(depot2) <- alag2
    lag(depot3) <- alag3

    # Observation variables (nmol/L).
    Cc       <- central      / vc
    Cc_acid  <- central_acid / vc_acid

    Cc      ~ add(addSd)      + prop(propSd)
    Cc_acid ~ add(addSd_acid) + prop(propSd_acid)
  })
}
