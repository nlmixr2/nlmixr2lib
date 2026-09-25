Ezuruike_2018_ethinylestradiol <- function() {
  description <- paste(
    "Two-compartment oral population pharmacokinetic reduction of the",
    "Simcyp minimal-PBPK-with-single-adjusting-compartment (SAC) model",
    "for ethinylestradiol (EE), the estrogen component of combined oral",
    "contraceptives, in healthy adult women (Ezuruike 2018). The source",
    "model was built in the Simcyp Population-based Simulator V17r1;",
    "the whole-body mass-balance equations and the virtual-population",
    "database are not published, so the platform model itself cannot be",
    "encoded. What IS fully published is the EE compound layer (Table 4)",
    "together with the model's own predicted clearance split (Figure 3b),",
    "and those are sufficient to rebuild the disposition as an ordinary",
    "compartmental model with NO fitted parameters: first-order",
    "absorption into a depot, exchange between a systemic compartment and",
    "the SAC (the paper's Kin / Kout, encoded as the canonical k12 / k21),",
    "linear elimination split into renal and hepatic arms, and gut-wall",
    "plus hepatic first-pass extraction on the dose.",
    "The reduction reproduces every one of the fifteen population-mean",
    "steady-state AUC cells of the paper's Table 4 [Table 3 in print]",
    "(three doses x five CYP-modulation scenarios) to within 2.6%,",
    "and the predicted steady-state Cmax after 35 ug once daily to 3.3%.",
    "Coadministered CYP modulators act through a single relative CYP3A4",
    "activity term that scales gut-wall intrinsic clearance, hepatic",
    "intrinsic clearance and systemic clearance together. Each",
    "modulator's activity was back-solved from its published AUC ratio",
    "alone; the published Cmax ratios were held out and are reproduced",
    "within 5.6% for all five arms that report one. The hypothetical",
    "pan-CYP inhibitor arm is a pure out-of-sample prediction (it reuses",
    "ketoconazole's activity, the paper having given both perpetrators",
    "the same Ki) and lands within 1.9% of the published AUC.",
    "This is a typical-value simulation model: the source reports no",
    "inter-individual variance components and no residual-error model,",
    "so there are no etas and propSd is fixed at zero. The paper's",
    "counts of virtual subjects crossing the 1000 and 1675 pg/mL.h",
    "breakthrough-bleeding and cardiovascular-risk thresholds therefore",
    "cannot be reproduced; the population means behind them can.",
    "See the vignette for the full list of deviations.",
    sep = " "
  )
  reference <- paste(
    "Ezuruike U, Humphries H, Dickins M, Neuhoff S, Gardner I,",
    "Rowland Yeo K. (2018). Risk-Benefit Assessment of Ethinylestradiol",
    "Using a Physiologically Based Pharmacokinetic Modeling Approach.",
    "Clin Pharmacol Ther 104(6):1229-1239. doi:10.1002/cpt.1085.",
    "Includes the publisher's correction of 11 May 2018 to Table 2.",
    sep = " "
  )
  vignette <- "Ezuruike_2018_ethinylestradiol"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Ezuruike 2018 Methods "PBPK model
  # development", which describes a minimal PBPK model with an
  # additional single adjusting compartment fitted to intravenous data.
  compartmentData <- list(
    depot = list(
      analyte = "ethinylestradiol",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ethinylestradiol",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ethinylestradiol",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # The coadministered CYP modulators the paper simulates. Every flag is
  # a mutually-exclusive binary indicator; all zero = EE alone. Each
  # coefficient in ini() is the relative CYP3A4 activity back-solved from
  # that arm's single published AUC ratio.
  covariateData <- list(
    CONMED_KETOCONAZOLE = list(
      description = "Concomitant ketoconazole (strong CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 4 [Table 3 in print], column 'CYP3A4",
        "inhibition': multiple doses of ketoconazole 200 mg twice daily,",
        "CYP3A4 Ki 0.015 umol/L. This is also the arm the authors used to",
        "refine fm CYP3A4 against the observed clinical interaction."
      ),
      source_name = "ketoconazole"
    ),
    CONMED_PANCYP_INH = list(
      description = paste(
        "Concomitant hypothetical pan-CYP inhibitor with potent",
        "inhibition of CYP1A2, CYP2C8, CYP2C9 and CYP3A4"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 4 [Table 3 in print] footnote b: a",
        "hypothetical compound with Ki 0.015 umol/L against all four",
        "CYPs, i.e. the same Ki the same table gives ketoconazole against",
        "CYP3A4. This model therefore gives it ketoconazole's relative",
        "CYP3A4 activity and applies that activity to CYP1A2, CYP2C8 and",
        "CYP2C9 as well, which makes this arm an out-of-sample",
        "prediction rather than a calibrated one."
      ),
      source_name = "hypothetical potent CYP1A2/2C8/2C9/3A4 inhibitor"
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole (strong CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 1, row 1: voriconazole 400 mg twice daily on",
        "day 18 then 200 mg twice daily on days 19-21, against EE 35 ug",
        "once daily for 21 days. The voriconazole compound model the",
        "authors built for this arm is tabulated in supplementary",
        "Table S2 but is not encoded here; see the vignette."
      ),
      source_name = "voriconazole"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole (moderate CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = "Ezuruike 2018 Table 1, row 2: fluconazole 300 mg single dose on day 7.",
      source_name = "fluconazole"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine (moderate CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 1, rows 3-4. The paper simulates two",
        "carbamazepine regimens (300 mg twice daily for 21 days with EE",
        "20 ug; and a titration to 300 mg twice daily with EE 35 ug) and",
        "predicts the SAME Cmax and AUC ratios for both, so a single",
        "indicator covers them. Only CYP3A4 induction was considered for",
        "carbamazepine (Ezuruike 2018 Discussion)."
      ),
      source_name = "carbamazepine (CBZ)"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz (moderate CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 4 [Table 3 in print] footnote d: efavirenz",
        "600 mg once daily, Indmax 9.9 and IndC50 3.8 umol/L. Simulated",
        "only; no clinical EE-efavirenz arm is reported."
      ),
      source_name = "efavirenz"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin (strong CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP modulator coadministered",
      notes = paste(
        "Ezuruike 2018 Table 1 rows 5-6 (300 mg and 600 mg once daily)",
        "and Table 4 [Table 3 in print] footnote c (600 mg once daily,",
        "Indmax 16, IndC50 0.32 umol/L). Pair with DOSE_RIFAMPICIN_MG,",
        "which selects between the two simulated dose levels. The source",
        "model also induced CYP2C9 in the rifampicin arms (Ezuruike 2018",
        "Discussion); here that is absorbed into the single back-solved",
        "CYP3A4 activity."
      ),
      source_name = "rifampicin"
    ),
    DOSE_RIFAMPICIN_MG = list(
      description = "Daily rifampicin dose driving the magnitude of CYP3A4 induction",
      units = "mg",
      type = "continuous",
      notes = paste(
        "Ezuruike 2018 Table 1: the only two simulated levels are 300 mg",
        "once daily (predicted EE AUC ratio 0.40) and 600 mg once daily",
        "(predicted EE AUC ratio 0.36). The model carries a linear",
        "deviation from a 600 mg reference, so those two levels are",
        "reproduced exactly and intermediate values interpolate. The",
        "column is only read when CONMED_RIFAMPICIN is 1; set it to 0",
        "otherwise."
      ),
      source_name = "Rifampicin 300 mg QD / Rifampicin 600 mg QD"
    )
  )

  # Reported by the source but deliberately not carried as covariates.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Ezuruike 2018 Table 4 expresses Vss and Vsac in",
        "L/kg, so the Simcyp model does scale distribution volume with",
        "weight, and every simulation matched its clinical study's",
        "demographics. This reduction fixes a 70 kg reference weight",
        "instead, because the weight distribution of the Simcyp healthy-",
        "volunteer population file is not published."
      ),
      units = "kg",
      type = "continuous",
      notes = "Simcyp healthy-volunteer population file; distribution not published."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 9L,
    age_range = "18-50 years (trial-matched; the reference control simulation used 20-50 years)",
    sex_female_pct = 100,
    disease_state = "Healthy female volunteers; no patient data were used.",
    dose_range = paste(
      "Ethinylestradiol 20, 35 and 50 ug once daily by mouth (and a",
      "single 50 ug intravenous dose for the distribution fit)."
    ),
    regions = "Simcyp healthy-volunteer virtual population, trial-matched demographics.",
    studies = paste(
      "The reference simulation is 10 trials x 10 healthy female",
      "volunteers aged 20-50 years (Ezuruike 2018 Figure 1a and Table 4",
      "[Table 3 in print]). The single adjusting compartment was fitted",
      "to an intravenous clinical study; fm CYP3A4 was refined against a",
      "ketoconazole interaction study; and six further drug-interaction",
      "trials (voriconazole n = 16; fluconazole n = 21; carbamazepine",
      "n = 10 and n = 10; rifampicin n = 22 and n = 12) were simulated",
      "with matched subject numbers and age ranges (supplementary",
      "Table S1)."
    ),
    notes = paste(
      "n_subjects records the 100 virtual subjects of the reference",
      "control simulation; n_studies counts the intravenous study, the",
      "ketoconazole refinement study, the six drug-interaction",
      "verification trials and the multiple-dose 35 ug control study.",
      "This is a PBPK analysis, not a population-PK fit: there is no",
      "pooled analysis dataset and no estimated variance components.",
      "The paper's own predicted mean split of systemic elimination",
      "(Figure 3b) is CYP3A4 22.19%, CYP2C9 10.32%, CYP1A2 7.22%,",
      "CYP2C8 1.06%, UGT1A1 5.24%, additional HLM 36.85% and renal",
      "17.26%; those seven numbers pin the clearance layer of this",
      "reduction and are carried as fm_* parameters below.",
      "Exposure thresholds used by the paper's risk-benefit assessment:",
      "a steady-state AUC below 1000 pg/mL.h was taken as the",
      "breakthrough-bleeding / contraceptive-failure threshold, and the",
      "population mean of 1675 pg/mL.h at 50 ug as the upper",
      "cardiovascular-risk threshold."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed: this is a typical-value simulation model.
    # Values are Ezuruike 2018 Table 4 inputs, arithmetic consequences of
    # them and of the Figure 3b clearance split, or (for the seven
    # perpetrator terms) relative CYP3A4 activities back-solved from the
    # paper's own published AUC ratios. Each case is labelled below.
    #
    # Shared quantities, all Ezuruike 2018 Table 4 unless noted:
    #   MW     = 296.4 g/mol
    #   B:P    = 1 (assumed by the authors), so blood == plasma throughout
    #   fu     = 0.015
    #   fa     = 0.948          fraction absorbed
    #   ka     = 1.103 1/h
    #   Qgut   = 11.74 L/h,  fu,gut = 1
    #   Vss    = 4.06 L/kg
    #   Kin    = 0.287,  Kout = 0.096   SAC transfer, fitted to i.v. data
    #   Vsac   = 2 L/kg
    #   CLR    = 2.079 L/h
    #   EG     = 0.44, EH = 0.25        Results, from Back et al. 1982
    #
    # Un-printed constant used in the derivations: a 70 kg reference body
    # weight, needed only to turn the L/kg volume input into litres.
    # ------------------------------------------------------------------

    # -- Absorption ----------------------------------------------------
    lka <- fixed(log(1.103))
    label("First-order absorption rate constant ka (1/h)")
    # Ezuruike 2018 Table 4, 'ka (1/h)' = 1.103, predicted by Simcyp from
    # polar surface area and hydrogen-bond donor count.

    fa <- fixed(0.948)
    label("Fraction of the oral dose absorbed from the gut lumen")
    # Ezuruike 2018 Table 4, 'fa' = 0.948, predicted from the same
    # physicochemical data. Figure 3a rounds it to 'Fa = 1'.

    # -- Distribution --------------------------------------------------
    # Ezuruike 2018 Table 4 gives a whole-body Vss of 4.06 L/kg taken
    # from the clinical literature, and SAC transfer constants Kin and
    # Kout that act on drug MASS, so they map exactly onto the canonical
    # k12 / k21. At steady state the SAC therefore holds Kin/Kout times
    # the systemic amount, and the systemic volume is
    #   vc = Vss * BW / (1 + Kin/Kout)
    #      = 4.06 * 70 / (1 + 0.287/0.096) = 71.2355 L
    # The implied peripheral volume, 212.96 L = 3.04 L/kg, is NOT the
    # Table 4 'Vsac' of 2 L/kg; under a mass-based Kin/Kout
    # parameterisation Vsac never enters the plasma prediction and is not
    # carried here. See the vignette Errata.
    lvc <- fixed(log(71.235509))
    label("Central (systemic) volume of distribution Vc (L)")

    lk12 <- fixed(log(0.287))
    label("Central-to-SAC transfer rate constant k12 (1/h)")
    # Ezuruike 2018 Table 4, 'K in' = 0.287. The table's unit column
    # reads 'L/h'; the quantity is a first-order rate constant in 1/h
    # (Simcyp's SAC inputs are kin and kout in 1/h), which is what makes
    # the printed Vss, Kin and Kout mutually consistent. See the
    # vignette Errata.

    lk21 <- fixed(log(0.096))
    label("SAC-to-central transfer rate constant k21 (1/h)")
    # Ezuruike 2018 Table 4, 'K out' = 0.096; same unit note as k12.

    # -- Systemic clearance --------------------------------------------
    lcl_renal <- fixed(log(2.079))
    label("Renal clearance of unchanged ethinylestradiol CLR (L/h)")
    # Ezuruike 2018 Table 4, 'CL R (L/h)' = 2.079. The Results derive it
    # as 6% of a mean oral clearance of 34.6 L/h.

    lcl_nonren <- fixed(log(9.983054))
    label("Hepatic metabolic clearance at baseline CLH (L/h)")
    # NOT printed directly; an arithmetic consequence of two printed
    # numbers. Figure 3b gives the optimized model's own predicted mean
    # split of systemic elimination, in which renal clearance is 17.26%
    # of the total. With CLR = 2.079 L/h that fixes the total systemic
    # clearance at 2.079 / 0.172359 = 12.0621 L/h (the Figure 3b slices
    # are renormalised to sum to 1; as printed they sum to 100.14%), and
    # the hepatic arm at 12.0621 - 2.079 = 9.9831 L/h.
    # NOTE this is NOT the 16.47 L/h intravenous clearance quoted in the
    # Results: that literature value was the TARGET of the retrograde
    # calculation that produced the initial hepatic CLint, before the
    # ketoconazole refinement. The value carried here is what the final
    # optimized model predicts, and it is the value that reproduces the
    # paper's own predicted exposures (16.47 L/h would under-predict
    # every Table 4 [Table 3 in print] cell by about 26%).

    eh <- fixed(0.25)
    label("Baseline hepatic extraction ratio")
    # Ezuruike 2018 Results, 'First-pass metabolism of EE in the gut and
    # liver': a gut extraction ratio (EG) of 0.44 and a liver extraction
    # ratio (EH) of 0.25, measured by simultaneous portal-vein and
    # systemic sampling (Back et al. 1982, the paper's reference 15).
    # Baseline hepatic availability is therefore 1 - 0.25 = 0.75.

    # -- Gut-wall first pass -------------------------------------------
    # Ezuruike 2018 Table 4 prints Qgut = 11.74 L/h and fu,gut = 1, and
    # the Results fix the gut availability FG at 0.56 (EG = 0.44). The
    # well-stirred 'Qgut' relation FG = Qgut / (Qgut + fu,gut * CLint,G)
    # therefore fixes the total unbound gut intrinsic clearance at
    #   CLint,G = 11.74 * 0.44 / 0.56 = 9.2243 L/h.
    # The Results apportion that total: about 70% to sulfation, about 20%
    # to hydroxylation and the rest to glucuronidation. The ketoconazole
    # refinement then doubled every CYP CLint input (fm CYP rose from
    # about 0.2 to about 0.4 while total clearance was held), so the
    # hydroxylation share of the final gut model is about 40%; the
    # intestinal CYPs are CYP3A4 and CYP2C9 only, split in the Figure 3b
    # ratio 22.19 : 10.32. That gives the three arms below, which sum to
    # 9.2243 L/h and so return EG = 0.44 exactly.
    # This 40% is the only quantity in the model derived by an argument
    # rather than transcribed, and it is corroborated twice: the
    # published voriconazole AUC ratio of 1.40 is unattainable unless it
    # exceeds about 0.22, and the held-out Cmax ratios are reproduced
    # best over 0.35-0.40. See the vignette sensitivity analysis.
    lcl_int_g_cyp3a4 <- fixed(log(2.518448))
    label("Gut-wall unbound intrinsic clearance via CYP3A4 (L/h)")

    lcl_int_g_cyp2c9 <- fixed(log(1.171266))
    label("Gut-wall unbound intrinsic clearance via CYP2C9 (L/h)")

    lcl_int_g_other <- fixed(log(5.534571))
    label("Gut-wall unbound intrinsic clearance via sulfation and glucuronidation (L/h)")

    qgut <- fixed(11.74)
    label("Qgut, the composite enterocytic blood flow / permeability term (L/h)")
    # Ezuruike 2018 Table 4, 'Q gut (L/h)' = 11.74.

    # -- Split of hepatic metabolic clearance by pathway ---------------
    # Ezuruike 2018 Figure 3b, the optimized model's own predicted mean
    # contributions to systemic elimination: CYP3A4 22.19%, CYP2C9
    # 10.32%, CYP1A2 7.22%, CYP2C8 1.06%, UGT1A1 5.24%, additional HLM
    # 36.85% and renal 17.26%. Renal is carried separately as lcl_renal,
    # so the six values below are each slice divided by the sum of the
    # six non-renal slices (82.88%), and they sum to exactly 1. The
    # CYP total of 0.4079 of systemic clearance and the CYP3A4 share of
    # 0.2219 are the 'increased fmCYP of about 0.4' and the 'increase in
    # fmCYP3A4 from 0.11 to 0.22' quoted in the Results.
    fm_cyp3a4 <- fixed(0.267736)
    label("Share of hepatic metabolic clearance mediated by CYP3A4")
    fm_cyp2c9 <- fixed(0.124517)
    label("Share of hepatic metabolic clearance mediated by CYP2C9")
    fm_cyp1a2 <- fixed(0.087114)
    label("Share of hepatic metabolic clearance mediated by CYP1A2")
    fm_cyp2c8 <- fixed(0.012790)
    label("Share of hepatic metabolic clearance mediated by CYP2C8")
    fm_ugt1a1 <- fixed(0.063224)
    label("Share of hepatic metabolic clearance mediated by UGT1A1")
    fm_other <- fixed(0.444619)
    label("Share of hepatic metabolic clearance via unassigned pathways (additional HLM)")

    # -- CYP modulator effects -----------------------------------------
    # Each coefficient is log(relative CYP3A4 activity) for that arm,
    # back-solved from the arm's published AUC ratio ALONE. The published
    # Cmax ratios were held out; they are reproduced as follows
    # (vignette Table 4):
    #   ketoconazole  AUCR 1.347 (Table 4 [Table 3 in print] columns)
    #   voriconazole  AUCR 1.40  (Table 1)   Cmax 1.28 -> 1.27  (-0.8%)
    #   fluconazole   AUCR 1.13  (Table 1)   Cmax 1.10 -> 1.09  (-0.8%)
    #   carbamazepine AUCR 0.61  (Table 1)   Cmax 0.66 -> 0.70  (+5.6%)
    #   rifampicin300 AUCR 0.40  (Table 1)   Cmax 0.51 -> 0.51  (-0.8%)
    #   rifampicin600 AUCR 0.36  (Table 1)   Cmax 0.49 -> 0.47  (-4.9%)
    #   efavirenz     AUCR 0.614 (Table 4 [Table 3 in print] columns)
    # The back-solved activities ladder monotonically with each
    # perpetrator's regulatory potency class -- voriconazole 0.131 <
    # ketoconazole 0.222 < fluconazole 0.666 < 1 < efavirenz 2.58 ~
    # carbamazepine 2.61 < rifampicin 300 mg 4.35 < rifampicin 600 mg
    # 4.85 -- which is a free consistency check on the whole set.
    e_conmed_ketoconazole_cyp3a4 <- fixed(log(0.2223))
    label("log relative CYP3A4 activity with ketoconazole 200 mg twice daily")

    e_conmed_voriconazole_cyp3a4 <- fixed(log(0.1314))
    label("log relative CYP3A4 activity with voriconazole")

    e_conmed_fluconazole_cyp3a4 <- fixed(log(0.6655))
    label("log relative CYP3A4 activity with fluconazole 300 mg single dose")

    e_conmed_cbz_cyp3a4 <- fixed(log(2.6055))
    label("log relative CYP3A4 activity with carbamazepine")

    e_conmed_efv_cyp3a4 <- fixed(log(2.5796))
    label("log relative CYP3A4 activity with efavirenz 600 mg once daily")

    e_conmed_rifampicin_cyp3a4 <- fixed(log(4.8546))
    label("log relative CYP3A4 activity with rifampicin 600 mg once daily")

    e_dose_rifampicin_cyp3a4 <- fixed(log(4.8546 / 4.3542))
    label("Change in log relative CYP3A4 activity per 300 mg of daily rifampicin")
    # The two simulated rifampicin levels give back-solved activities of
    # 4.3542 (300 mg once daily, Table 1 AUC ratio 0.40) and 4.8546
    # (600 mg once daily, Table 1 AUC ratio 0.36). Encoded as a linear
    # deviation from a 600 mg reference so both levels are exact.

    propSd <- fixed(0)
    label("Proportional residual error (none reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Typical-value parameters. No random effects.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    cl_renal <- exp(lcl_renal)
    cl_nonren <- exp(lcl_nonren)

    # ------------------------------------------------------------------
    # 2. Relative CYP activity contributed by a coadministered modulator.
    # All indicators zero (EE alone) gives a3a4 = 1 and aoth = 1. The
    # terms are additive on the log scale, so the flags are
    # multiplicative on activity.
    # The hypothetical pan-CYP inhibitor is given ketoconazole's CYP3A4
    # activity, because the source assigns both the same Ki of
    # 0.015 umol/L; it is the only arm that also moves the non-CYP3A4
    # CYPs, and the only arm whose AUC is not used to calibrate anything.
    # ------------------------------------------------------------------
    lact_rif <- e_conmed_rifampicin_cyp3a4 +
      e_dose_rifampicin_cyp3a4 * (DOSE_RIFAMPICIN_MG - 600) / 300

    a3a4 <- exp(
      e_conmed_ketoconazole_cyp3a4 * CONMED_KETOCONAZOLE +
        e_conmed_ketoconazole_cyp3a4 * CONMED_PANCYP_INH +
        e_conmed_voriconazole_cyp3a4 * CONMED_VORICONAZOLE +
        e_conmed_fluconazole_cyp3a4 * CONMED_FLUCONAZOLE +
        e_conmed_cbz_cyp3a4 * CONMED_CBZ +
        e_conmed_efv_cyp3a4 * CONMED_EFV +
        lact_rif * CONMED_RIFAMPICIN
    )
    aoth <- exp(e_conmed_ketoconazole_cyp3a4 * CONMED_PANCYP_INH)

    # ------------------------------------------------------------------
    # 3. Hepatic arm. Multiplying an intrinsic clearance by a factor A
    # takes an extraction ratio E to E*A / (1 - E + E*A), so hepatic
    # availability 1 - E is divided by hfac and the hepatic clearance is
    # multiplied by m / hfac. The UGT1A1 and unassigned pathways are not
    # modulated by any of the perpetrators the paper simulates.
    # ------------------------------------------------------------------
    fmcyp_oth <- fm_cyp1a2 + fm_cyp2c8 + fm_cyp2c9
    m <- fm_ugt1a1 + fm_other + fm_cyp3a4 * a3a4 + fmcyp_oth * aoth
    hfac <- 1 - eh + eh * m
    fhep <- (1 - eh) / hfac
    cl <- cl_renal + cl_nonren * m / hfac
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. Gut-wall arm, re-derived through the same well-stirred relation
    # so that the extraction ratio can never leave [0, 1]. fu,gut = 1
    # (Table 4), so the unbound and total gut intrinsic clearances are
    # the same number.
    # ------------------------------------------------------------------
    cl_int_g <- exp(lcl_int_g_cyp3a4) * a3a4 +
      exp(lcl_int_g_cyp2c9) * aoth +
      exp(lcl_int_g_other)
    egut <- cl_int_g / (qgut + cl_int_g)

    # ------------------------------------------------------------------
    # 5. ODE system. `peripheral1` is the paper's single adjusting
    # compartment; k12 and k21 act on masses (ug), matching the paper's
    # Kin / Kout definitions.
    # ------------------------------------------------------------------
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 6. Oral bioavailability, F = fa * (1 - EG) * (1 - EH). At baseline
    # this is 0.948 * 0.56 * 0.75 = 0.3982, consistent with the
    # 'F ~ 50%' of Figure 3a and with oral bioavailability being less
    # than 50% (Results).
    # ------------------------------------------------------------------
    fdepot <- fa * (1 - egut) * fhep
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 7. Observation. Doses are in ug and vc is in L, so central / vc is
    # in ug/L = ng/mL; multiply by 1000 to report pg/mL, the unit the
    # source uses in Tables 2 and 4 [Tables 2 and 3 in print].
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
