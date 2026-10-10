Liu_2022_nifedipine <- function() {
  description <- paste(
    "Direct (ordinary, hyperbolic) Emax pharmacodynamic model of the change in",
    "systolic blood pressure (SBP) produced by the nifedipine plasma",
    "concentration in adults, the PD layer of the Liu 2022 nifedipine",
    "PBPK/PD model used to evaluate the CYP3A4-mediated drug-drug",
    "interaction with apatinib (Liu 2022 equation 3). The concentration-effect",
    "relationship is direct, with no effect compartment and no baseline term:",
    "the output dsbp is the drug-induced change in SBP (mmHg), negative for a",
    "reduction. Only this PD layer is reproduced here: the PK layer is a",
    "Simcyp minimal-PBPK model whose clearance is entered as recombinant",
    "CYP3A4 / CYP3A5 Vmax and Km per pmol of enzyme and whose controlled-",
    "release absorption uses an ADAM model with an unprinted Weibull",
    "dissolution profile, so it cannot be reproduced from the publication",
    "(see the vignette 'Assumptions and deviations'). The nifedipine plasma",
    "concentration is therefore supplied by the user as the canonical PD-driver",
    "covariate CEFFECT (ng/mL). Emax and EC50 were set by the authors (EC50",
    "from the literature, Emax chosen as the best of a range of literature",
    "values) rather than estimated, and no inter-individual variability or",
    "residual-error model is reported, so the model is typical-value only.",
    sep = " "
  )
  reference <- paste(
    "Liu H, Yu Y, Liu L, Wang C, Guo N, Wang X, Xiang X, Han B. (2022).",
    "Application of physiologically-based pharmacokinetic/pharmacodynamic",
    "models to evaluate the interaction between nifedipine and apatinib.",
    "Frontiers in Pharmacology 13:970539. doi:10.3389/fphar.2022.970539.",
    sep = " "
  )
  vignette <- "Liu_2022_nifedipine"
  units <- list(
    time = "h",
    dosing = "(not applicable; nifedipine plasma concentration is supplied as the CEFFECT covariate, not as a dose record)",
    concentration = "ng/mL (nifedipine plasma concentration via CEFFECT; the PD output dsbp is in mmHg)"
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Nifedipine plasma concentration (ng/mL), supplied as the direct driver of the Emax reduction in systolic blood pressure",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Member of the canonical PD-driver family CEFFECT; here the biophase is",
        "plasma itself, because Liu 2022 relates the change in SBP directly to",
        "the nifedipine concentration through an ordinary Emax model (equation",
        "3) with no effect compartment, citing 'a close and direct relationship",
        "between circulating drug concentrations and antihypertensive effect'",
        "(Discussion). EC50 is given as 12.12 ng/mL, so CEFFECT must be in",
        "ng/mL; the paper does not say whether this is total or unbound",
        "concentration, and every nifedipine concentration it reports",
        "(Tables 3 and 4, Figures 1 and 2) is total plasma, so total plasma is",
        "assumed. In the source study CEFFECT is the plasma profile predicted",
        "by the authors' Simcyp v16 minimal-PBPK model (Table 1: fu 0.039,",
        "B/P 0.685, Vss 0.57 L/kg, IR first-order ka 4.6 1/h with fa 1, CR",
        "tablet via ADAM, rCYP3A4 Vmax 22 pmol/min/pmol and Km 10.95 uM,",
        "rCYP3A5 Vmax 3.5 pmol/min/pmol and Km 31.9 uM). Its predicted mean",
        "Cmax was 26.51 ng/mL after 60 mg CR and 144.59 ng/mL after 20 mg IR",
        "in healthy volunteers (Table 3), and 20.97 ng/mL (alone) versus",
        "29.51 ng/mL (with apatinib 750 mg once daily, competitive CYP3A4",
        "inhibition, Ki 0.12 uM, fu,mic 0.65) after 30 mg CR (Table 4). That",
        "model is not reproduced in this file and cannot be reproduced from",
        "the publication. Users supply CEFFECT from their own nifedipine PK",
        "model or from observed plasma concentrations; set it to 0 for",
        "drug-free records."
      ),
      source_name = "C (Liu 2022 equation 3)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "26-65 years (Simcyp virtual populations)",
    weight_range = "not reported",
    sex_female_pct = 50,
    disease_state = paste(
      "PD layer checked against published mean SBP changes in hypertensive",
      "patients (Figure 3; Toal 2012 for 60 mg CR, single 10 mg IR dose) and",
      "then applied to virtual cancer patients and patients with Child-Pugh",
      "A / B / C hepatic impairment receiving apatinib"
    ),
    dose_range = paste(
      "Nifedipine single oral doses: 60 mg and 30 mg controlled-release (CR)",
      "tablet, 10 mg, 20 mg and 30 mg immediate-release (IR) tablet; with or",
      "without apatinib 750 mg orally once daily for 8 days (nifedipine on",
      "day 6)."
    ),
    regions = "China (authors' institution; Sim-Chinese population for the DDI verification)",
    notes = paste(
      "No subject-level data were fitted. The PD parameters are literature",
      "values: EC50 12.12 ng/mL 'reported ... for nifedipine at the regular",
      "doses' (Hirasawa 1985; Levine 2003; Meredith and Elliott 2004;",
      "Niu 2021), and Emax -30 mmHg selected as the best-fitting value from",
      "the literature range against observed mean SBP changes after 60 mg CR",
      "and 10 mg IR nifedipine in hypertensive patients (Results",
      "'Verification of pharmacodynamic model for nifedipine'; Table 5).",
      "All simulations in the paper used Simcyp virtual populations of",
      "100 subjects aged 26-65 years, 50% female (24 subjects for the",
      "Sim-Chinese DDI verification), so n_subjects and n_studies are not",
      "applicable."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Liu 2022 equation 3 (Methods, 'Development and verification of
    # pharmacodynamic model for nifedipine'):
    #
    #          Emax * C
    #   E = -----------
    #        EC50 + C
    #
    # 'An ordinary Emax model (Eq. 3) was used to describe the relationship
    # between nifedipine concentration and the change of SBP ... The Emax
    # and EC50 were obtained from literature, which were set at -30 mmHg
    # and 12.12 ng/ml, respectively.'
    #
    # The paper's Emax is signed (-30 mmHg, a reduction). Here the magnitude
    # 30 mmHg is carried on the log scale and the sign is applied in
    # model(), so dsbp = -emax * C / (ec50 + C) reproduces equation 3.
    # ------------------------------------------------------------------

    lemax <- fixed(log(30))
    label("Log of the Emax magnitude, maximum nifedipine-induced reduction in systolic blood pressure (mmHg)")
    # Liu 2022 Methods PD paragraph: Emax 'set at -30 mmHg'. Results,
    # 'Verification of pharmacodynamic model for nifedipine': 'When the Emax
    # was set to -30 mmHg, the PD model fitted best.' Discussion: chosen from
    # the 'range of Emax values reported in the literature' by comparing
    # predicted and observed SBP changes. A selected value with no standard
    # error, so it is fixed.

    lec50 <- fixed(log(12.12))
    label("Log of EC50, nifedipine plasma concentration giving half the maximum SBP reduction (ng/mL)")
    # Liu 2022 Methods PD paragraph: EC50 'set at ... 12.12 ng/ml'.
    # Discussion: 'The EC50 was reported to be 12.12 ng/ml for nifedipine at
    # the regular doses ... So, the EC50 value was set to 12.12 ng/ml in this
    # study.' A literature value, not estimated, so it is fixed.

    addSd <- fixed(0)
    label("Additive residual SD on the change in systolic blood pressure (mmHg; zero, no residual-error model is reported)")
  })

  model({
    # 1. Nifedipine plasma concentration for this record, ng/mL, supplied
    #    through the canonical CEFFECT covariate column. The link is direct
    #    (no effect compartment), as in Liu 2022 equation 3.
    conc <- CEFFECT

    # 2. Typical-value PD parameters; no IIV is reported.
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # 3. Liu 2022 equation 3 with the signed Emax of -30 mmHg: the
    #    drug-induced change in SBP (mmHg), negative for a reduction.
    dsbp <- -emax * conc / (ec50 + conc)

    # 4. Observation. addSd is fixed at zero because the source reports no
    #    residual-error model; free it to refit to individual SBP data.
    dsbp ~ add(addSd)
  })
}
