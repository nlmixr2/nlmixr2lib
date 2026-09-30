Kurosawa_2021_bedaquiline <- function() {
  description <- paste0(
    "Four-compartment population PK model for oral bedaquiline with the dual ",
    "zero-order absorption input of McLeay 2014, updated by Kurosawa 2021 with ",
    "the effect of steady-state clarithromycin (500 mg every 12 hours) on ",
    "bedaquiline apparent clearance, estimated from a phase 1 crossover study in ",
    "16 healthy adults (NCT03800550). A fraction FR1 of each dose enters a depot ",
    "as a zero-order input over DUR1 after a formulation-dependent lag Tlag and ",
    "passes to the central compartment at a fixed rate of 1,000 1/h; the ",
    "remaining 1 - FR1 enters the central compartment directly as a zero-order ",
    "input over DUR2 after a lag of Tlag + Tlag,add. Every structural, ",
    "covariate, variability and residual parameter is fixed to the McLeay 2014 ",
    "final estimates (study on relative bioavailability, Black race and ",
    "healthy-volunteer / drug-sensitive-TB status on CL/F, female sex on Vc/F); ",
    "only the clarithromycin effect was estimated, CL/F x (1 - 0.37) with ",
    "clarithromycin co-administration."
  )
  reference <- paste(
    "Kurosawa K, Rossenu S, Biewenga J, Ouwerkerk-Mahadevan S, Willems W,",
    "Ernault E, Kambili C (2021). Population Pharmacokinetic Analysis of",
    "Bedaquiline-Clarithromycin for Dose Selection Against Pulmonary",
    "Nontuberculous Mycobacteria Based on a Phase 1, Randomized,",
    "Pharmacokinetic Study. Journal of Clinical Pharmacology 61(10):1344-1355.",
    "doi:10.1002/jcph.1887. All parameters other than the clarithromycin",
    "effect are fixed to the previously developed model reprinted in",
    "Kurosawa 2021 Supplemental Table S1: McLeay SC, Vis P, van Heeswijk RPG,",
    "Green B (2014). Population pharmacokinetics of bedaquiline (TMC207), a",
    "novel antituberculosis drug. Antimicrobial Agents and Chemotherapy",
    "58(9):5315-5324. doi:10.1128/AAC.01418-13 (structure from McLeay 2014",
    "Fig. 1 and Results; values cross-checked against McLeay 2014 Table 3).",
    sep = " "
  )
  vignette <- "Kurosawa_2021_bedaquiline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  paper_specific_residual_sds <- c("expSd", "expSdC208C209")

  compartmentData <- list(
    depot = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "bedaquiline",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "bedaquiline",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral3 = list(
      analyte = "bedaquiline",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONMED_CLARITHROMYCIN = list(
      description = "Clarithromycin co-administration at steady state (1 = on clarithromycin 500 mg every 12 hours, 0 = bedaquiline alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (bedaquiline monotherapy)",
      notes = paste0(
        "Kurosawa 2021 Methods equation CL/F = CLpop * (1 + theta)^CLRi with ",
        "CLRi = 0 for no coadministration and 1 for coadministration; theta = ",
        "-0.37 (RSE 11%, Table S1; Results 95% CI -45% to -29%). In the study, ",
        "clarithromycin was given from day 1 to day 14 of treatment B and ",
        "bedaquiline on day 5, so steady-state inhibition was reached before ",
        "the bedaquiline dose. May be supplied time-varying."
      ),
      source_name = "CLR"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = any other race)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "McLeay 2014 Table 3 / Kurosawa 2021 Table S1: CL/F 52.0% higher in Black subjects, CL/F * (1 + 0.520)^RACE_BLACK (McLeay 2014 equation 5). Kurosawa 2021 simulated non-Black subjects only.",
      source_name = "black race"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "McLeay 2014 Table 3 / Kurosawa 2021 Table S1: Vc/F 15.7% lower in females, Vc/F * (1 - 0.157)^SEXF (McLeay 2014 equation 5).",
      source_name = "sex"
    ),
    DIS_TB_MDR = list(
      description = "Multidrug-resistant tuberculosis patient indicator (1 = MDR-TB patient, 0 = healthy volunteer or drug-sensitive TB patient)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 in the source parameterisation (MDR-TB patients of phase IIb studies C208 / C209 are the reference group); 0 applies the +37.5% CL/F increase",
      notes = paste0(
        "McLeay 2014 Results: 'healthy volunteers and patients with DS-TB were ",
        "added as a covariate to CL/F, with the effect described as an ",
        "increase in CL/F compared to the phase IIb studies in patients with ",
        "MDR-TB'; CL/F * (1 + 0.375)^(1 - DIS_TB_MDR). In the McLeay 2014 ",
        "dataset every MDR-TB patient came from C208 / C209. The Kurosawa 2021 ",
        "healthy volunteers carry DIS_TB_MDR = 0. Kurosawa 2021 simulated ",
        "pulmonary nontuberculous-mycobacteria patients assuming a disease ",
        "status similar to MDR-TB, i.e. without the 37.5% increase (Table S1 ",
        "footnote a lists the covariates used in those simulations and does ",
        "not include this one); reproduce that by setting DIS_TB_MDR = 1."
      ),
      source_name = "subject status (healthy volunteer / DS-TB vs MDR-TB)"
    ),
    STUDY_BDQ_C208_C209 = list(
      description = "Bedaquiline phase IIb MDR-TB studies TMC207-TiDP13-C208 / C209 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = paste0(
        "Relative bioavailability reference group (F = 1) and the higher ",
        "residual error (27.7% vs 20.6%) of McLeay 2014 Table 3. Mutually ",
        "exclusive with STUDY_BDQ_CDE102_C104; a subject with both indicators ",
        "0 is in the 'other studies' group (F = 2.03), which is where the ",
        "Kurosawa 2021 healthy volunteers sit (Table S1 footnote b). The ",
        "Kurosawa 2021 Table 5 regimen simulations are reproduced only with ",
        "F = 1, i.e. STUDY_BDQ_C208_C209 = 1 (with F = 2.03 every simulated ",
        "exposure is about 2.03-fold too high), although Table S1 footnote a ",
        "lists 'Other studies on F' among the simulation covariates; see the ",
        "vignette."
      ),
      source_name = "study"
    ),
    STUDY_BDQ_CDE102_C104 = list(
      description = "Bedaquiline phase 1 oral-solution studies R207910-CDE102 / TMC207-TiDP13-C104 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = "Relative bioavailability 1.51 versus the C208 / C209 reference (McLeay 2014 Table 3; Kurosawa 2021 Table S1). Mutually exclusive with STUDY_BDQ_C208_C209.",
      source_name = "study"
    ),
    FORM_TABLET = list(
      description = "Bedaquiline tablet formulation indicator (1 = tablet, 0 = oral solution)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral solution)",
      notes = "Selects the absorption lag time: Tlag = 0.541 h for the oral solution and 0.917 h for the tablet (McLeay 2014 Table 3; Kurosawa 2021 Table S1 'ALAG1 solution' / 'ALAG1 tablet'). The Kurosawa 2021 study and simulations used the 100 mg commercial tablet.",
      source_name = "formulation"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "24-55 years (median 43)",
    weight_range = "56.8-100.0 kg (median 65.65)",
    bmi_range = "18.9-29.5 kg/m^2 (median 22.44)",
    sex_female_pct = 56.3,
    race_ethnicity = c(White = 100),
    disease_state = "Healthy adults",
    dose_range = paste0(
      "Single oral 100 mg bedaquiline tablet with a standardized breakfast, ",
      "alone (treatment A) or on day 5 of clarithromycin 500 mg every 12 hours ",
      "for 14 days (treatment B); 2-sequence crossover with a washout of at ",
      "least 28 days."
    ),
    regions = "Belgium (single centre), March-June 2019",
    notes = paste0(
      "Kurosawa 2021 Table 1. Only the clarithromycin effect was estimated on ",
      "these 16 subjects; every other parameter is fixed to McLeay 2014, which ",
      "was fit to 5,222 bedaquiline concentrations from 480 subjects (111 ",
      "healthy volunteers, 44 drug-sensitive TB and 325 MDR-TB patients) in 9 ",
      "phase I/II studies, weight 30-113 kg, age 18-68 years (McLeay 2014 ",
      "Table 1 and Results)."
    )
  )

  ini({
    # ==================================================================
    # Disposition. Kurosawa 2021 Table S1 'Previously developed model'
    # column (= McLeay 2014 Table 3). The updated model re-estimated none
    # of these (Table S1 'Updated model' column shows '-'), so all are
    # fixed. Canonical mapping: CLp1-3/F -> q, q2, q3; Vp1-3/F -> vp,
    # vp2, vp3. Reference individual: male, non-Black, MDR-TB patient of
    # study C208 / C209 (McLeay 2014 Discussion).
    # ==================================================================
    lcl  <- fixed(log(2.78)); label("Apparent clearance CL/F, reference individual (L/h)")        # Table S1: CL/F = 2.78 L/h (McLeay 2014 Table 3: 2.78, RSE 5.1%)
    lvc  <- fixed(log(164));  label("Apparent central volume Vc/F, reference individual (L)")     # Table S1: Vc/F = 164 L (McLeay 2014 Table 3: 164, RSE 5.0%)
    lq   <- fixed(log(11.8)); label("Apparent intercompartmental clearance CLp1/F (L/h)")         # Table S1: CLp1/F = 11.8 L/h
    lvp  <- fixed(log(178));  label("Apparent first peripheral volume Vp1/F (L)")                 # Table S1: Vp1/F = 178 L
    lq2  <- fixed(log(8.03)); label("Apparent intercompartmental clearance CLp2/F (L/h)")         # Table S1: CLp2/F = 8.03 L/h
    lvp2 <- fixed(log(3010)); label("Apparent second peripheral volume Vp2/F (L)")                # Table S1: Vp2/F = 3010 L
    lq3  <- fixed(log(3.58)); label("Apparent intercompartmental clearance CLp3/F (L/h)")         # Table S1: CLp3/F = 3.58 L/h
    lvp3 <- fixed(log(7350)); label("Apparent third peripheral volume Vp3/F (L)")                 # Table S1: Vp3/F = 7350 L

    # ==================================================================
    # Dual zero-order input (McLeay 2014 Fig. 1 and Results). The
    # transfer rate out of the input compartment was fixed to 1,000 1/h
    # so that the depot arm is effectively zero-order. Tlag delays the
    # depot arm; the central arm starts at Tlag + Tlag,add (McLeay 2014
    # Fig. 1 labels its lag 'Tlag + Tlag,add').
    # ==================================================================
    lka       <- fixed(log(1000));         label("Transfer rate constant from the input compartment to central (1/h)")      # McLeay 2014 Results / Fig. 1: 'the rate parameter ... was fixed to a high rate (1,000)'
    logitfrel <- fixed(qlogis(0.585));     label("Logit of the fraction FR1 of the bioavailable dose entering the depot arm (unitless)")  # Table S1: FR1 = 58.5%; McLeay 2014 Results: 'constrained to a value between 0 and 1 using a logit function'
    ld1       <- fixed(log(2.22));         label("Duration of the zero-order input into the depot DUR1 (h)")                  # Table S1: D1 = 2.22 h
    ld2       <- fixed(log(1.48));         label("Duration of the zero-order input into central DUR2 (h)")                    # Table S1: D2 = 1.48 h
    ltlag     <- fixed(log(0.541));        label("Absorption lag time Tlag, oral solution (h)")                                # Table S1: 'ALAG1 solution' = 0.541 h
    ltlag_tablet <- fixed(log(0.917));     label("Absorption lag time Tlag, tablet (h)")                                       # Table S1: 'ALAG1 tablet' = 0.917 h
    ltlag_add <- fixed(log(1.48));         label("Additional lag time before the direct input into central Tlag,add (h)")      # Table S1: TLAG = 1.48 h ('additional lag time for absorption for the second pathway')

    # ==================================================================
    # Relative bioavailability by study (McLeay 2014: C208 / C209 is the
    # reference with F = 1).
    # ==================================================================
    lfdepot                 <- fixed(log(1)); label("Relative bioavailability, reference studies C208 / C209 (fraction)")   # McLeay 2014 Results: phase IIb MDR-TB studies 'fixed as the reference group for both F and CL/F'
    e_study_bdq_cde102_c104_f <- fixed(1.51); label("Relative bioavailability in studies CDE102 / C104 (ratio to reference)") # Table S1: 'Study R207910-CDE102 or TiDP13-C104 on F' = 1.51
    e_study_other_f         <- fixed(2.03);   label("Relative bioavailability in all other studies (ratio to reference)")     # Table S1: 'Other studies on F' = 2.03

    # ==================================================================
    # Covariate effects on disposition, (1 + theta)^cov form (McLeay 2014
    # equation 5; Kurosawa 2021 Methods equation).
    # ==================================================================
    e_race_black_cl <- fixed(0.520);  label("Fractional increase in CL/F for Black race (fraction)")                                  # Table S1: 'Increase in CL with Black race (%)' = 52.0
    e_nonmdr_cl     <- fixed(0.375);  label("Fractional increase in CL/F for healthy volunteers or DS-TB patients (fraction)")        # Table S1: 'Increase in CL for healthy volunteers or C202 (%)' = 37.5
    e_sexf_vc       <- fixed(-0.157); label("Fractional change in Vc/F for female sex (fraction)")                                    # Table S1: 'Decrease in Vc with female sex (%)' = -15.7
    e_conmed_clarithromycin_cl <- -0.37; label("Fractional change in CL/F with clarithromycin co-administration (fraction)")          # Table S1: 'Effect CLR on CL/F' = -0.37 (RSE 11%); Results: -37% (95% CI -45% to -29%)

    # ==================================================================
    # Between-subject variability (Table S1 'BSV' column). The printed
    # value is read as the SD of eta (omega = BSV / 100), not as a
    # back-transformed CV: only that reading reproduces the SD / mean
    # ratios of the Kurosawa 2021 Table 5 simulated exposures (see the
    # vignette). The CL/Vc correlation is 0.407; FR1 variability is on
    # the logit scale.
    # ==================================================================
    etalcl + etalvc ~ fixed(c(
      0.254016,
      0.080205, 0.152881
    ))                                   # Table S1: BSV CL/F 50.4, Vc/F 39.1 -> 0.504^2, 0.391^2; 'Correlation CL/Vc' 0.407 -> cov 0.407 * 0.504 * 0.391
    etalogitfrel ~ fixed(1.2769)         # Table S1: BSV FR1 = 113 -> 1.13^2 (logit scale)
    etalfdepot   ~ fixed(0.156816)       # Table S1: 'Between-subject variability on F' = 39.6 -> 0.396^2

    # ==================================================================
    # Residual error: log-transform-both-sides additive error (McLeay 2014
    # equation 2) = lognormal error on the linear scale.
    # ==================================================================
    expSd         <- fixed(0.206); label("Log-scale residual SD, all studies except C208 / C209 (log units)")  # Table S1: RUV = 20.6 CV%
    expSdC208C209 <- fixed(0.277); label("Log-scale residual SD, studies C208 / C209 (log units)")            # Table S1: 'RUV on TiDP13-C208 or TiDP13-C209' = 27.7 CV%
  })
  model({
    # 1. Individual parameters (McLeay 2014 equations 1 and 5; Kurosawa
    #    2021 Methods: CL/F = CLpop * (1 + theta)^CLRi).
    cl <- exp(lcl + etalcl) *
      (1 + e_race_black_cl)^RACE_BLACK *
      (1 + e_nonmdr_cl)^(1 - DIS_TB_MDR) *
      (1 + e_conmed_clarithromycin_cl)^CONMED_CLARITHROMYCIN
    vc <- exp(lvc + etalvc) * (1 + e_sexf_vc)^SEXF
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    q3 <- exp(lq3)
    vp3 <- exp(lvp3)

    ka <- exp(lka)
    d1 <- exp(ld1)
    d2 <- exp(ld2)
    tlag <- exp(ltlag * (1 - FORM_TABLET) + ltlag_tablet * FORM_TABLET)
    tlag_add <- exp(ltlag_add)

    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    # Study group on relative F: C208 / C209 = 1 (reference), CDE102 /
    # C104 = 1.51, all other studies = 2.03.
    study_other <- 1 - STUDY_BDQ_C208_C209 - STUDY_BDQ_CDE102_C104
    fstudy <- e_study_bdq_cde102_c104_f^STUDY_BDQ_CDE102_C104 *
      e_study_other_f^study_other
    fbio <- exp(lfdepot + etalfdepot) * fstudy

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    # 3. Four-compartment disposition fed by the input compartment
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2 -
      k14 * central + k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3

    # 4. Dual zero-order input. Each oral administration is TWO dose
    #    records at the same time, each carrying the WHOLE dose amount and
    #    a zero-order duration (rate = -2, or an explicit dur of DUR1 /
    #    DUR2): one to "depot", one to "central". The f() terms split it.
    f(depot) <- fbio * frel
    dur(depot) <- d1
    alag(depot) <- tlag

    f(central) <- fbio * (1 - frel)
    dur(central) <- d2
    alag(central) <- tlag + tlag_add

    # 5. Observation: mg / L. Log-transform-both-sides residual error,
    #    larger in the phase IIb studies C208 / C209.
    Cc <- central / vc
    w_rv <- expSd * (1 - STUDY_BDQ_C208_C209) + expSdC208C209 * STUDY_BDQ_C208_C209
    Cc ~ lnorm(w_rv)
  })
}
