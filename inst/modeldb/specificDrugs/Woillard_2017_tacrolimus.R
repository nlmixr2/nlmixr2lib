Woillard_2017_tacrolimus <- function() {
  description <- paste(
    "One-compartment population PK model for immediate-release oral tacrolimus",
    "in the early period after cadaveric kidney transplantation (Woillard 2017),",
    "with two parallel gamma-distributed absorption routes (double-gamma",
    "absorption), first-order elimination and a dose-proportional steady-state",
    "trough offset C0. A three-level CYP3A metabolizer cluster (poor /",
    "intermediate / extensive, reconstructed inside model() from the recipient",
    "CYP3A5 expresser status and the CYP3A4*22 rs35599367 carrier indicator)",
    "acts as an ordinal power-model multiplier on the whole predicted whole-blood",
    "concentration. Fitted non-parametrically in Pmetrics and evaluated here in",
    "the closed form published by the authors, so the model describes ONE",
    "steady-state dosing interval.",
    sep = " "
  )
  reference <- paste(
    "Woillard JB, Mourad M, Neely M, Capron A, van Schaik RH, van Gelder T,",
    "Lloberas N, Hesselink DA, Marquet P, Haufroid V, Elens L.",
    "Tacrolimus Updated Guidelines through popPK Modeling: How to Benefit More",
    "from CYP3A Pre-emptive Genotyping Prior to Kidney Transplantation.",
    "Front Pharmacol. 2017;8:358. doi:10.3389/fphar.2017.00358",
    sep = " "
  )
  vignette <- "Woillard_2017_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte = "tacrolimus",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = paste(
        "Dose-registration state only, and it is deliberately always empty.",
        "Woillard 2017 Supplemental Data 1 publishes C(t) in closed form (the",
        "analytic convolution of the double-gamma absorption rate with the",
        "one-compartment disposition function) and states that this expression",
        "is what was entered as the output equation of the Pmetrics model file,",
        "so this file evaluates the same closed form algebraically rather than",
        "integrating an equivalent ODE system. rxode2 nevertheless requires a",
        "dosed compartment with a defined d/dt() for podo() and tad() to",
        "resolve, which is what this state provides. f(depot) <- 0 keeps it at",
        "zero, so depot is NOT a meaningful amount; the whole-blood",
        "concentration is Cc."
      )
    )
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (rs776746, CYP3A5*3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = paste(
        "1 = at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3);",
        "0 = CYP3A5*3/*3. TaqMan allelic discrimination (Woillard 2017 Methods,",
        "'Genotyping Analysis'). Not used on its own in the final model: it is",
        "one of the two inputs from which the three-level CYP3A metabolizer",
        "cluster (PM / IM / EM) is reconstructed inside model(), following the",
        "two-binary-input convention documented in the SNP_CYP3A4_RS35599367",
        "and CYP3A5_EXPR_DONOR entries of",
        "inst/references/covariate-columns.md. Cohort distribution (Table 1):",
        "*1/*1 4 (6.8%), *1/*3 14 (23.7%), *3/*3 41 (69.5%), so",
        "CYP3A5_EXPR = 1 in 18 of 59 subjects -- exactly the 18 patients the",
        "paper classifies as extensive metabolizers.",
        sep = " "
      ),
      source_name = "CYP3A5*3"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) reduced-function allele carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A4*1/*1 wild-type homozygote)",
      notes = paste(
        "1 = carries at least one T (*22) allele (*1/*22 or *22/*22); 0 =",
        "*1/*1. TaqMan allelic discrimination (Woillard 2017 Methods,",
        "'Genotyping Analysis'). Second of the two inputs to the CYP3A cluster",
        "reconstruction in model(). Cohort distribution (Table 1): *1/*22",
        "5 (8.5%), *1/*1 54 (91.5%); no *22/*22 homozygote was observed, and",
        "all 5 carriers were also CYP3A5*3/*3, so the 5 carriers are exactly",
        "the 5 patients the paper classifies as poor metabolizers.",
        sep = " "
      ),
      source_name = "CYP3A4*22"
    )
  )

  # Covariates that Woillard 2017 screened but did NOT retain in the final
  # model. Documented for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened in the univariate analysis and also tested as an allometric",
        "scaling term; 'Allometric scaling of age and bodyweight did not",
        "significantly decrease -2LL' (Results, 'Development of the Structural",
        "Model'). Cohort mean 70.4 +/- 13.9 kg (Table 1).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened univariately and as an allometric term; not significant. The",
        "Discussion attributes the null result to the cohort's narrow, elderly",
        "age distribution (mean 51.9 +/- 13.4 years, Table 1).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened univariately; not associated with any Bayesian posterior PK",
        "parameter. Cohort 38 female (64.4%) / 21 male (35.6%) (Table 1).",
        sep = " "
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened univariately; not associated with any Bayesian posterior PK",
        "parameter. Cohort mean 60.1 +/- 20.0 mL/min at the PK course",
        "(Table 1).",
        sep = " "
      )
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = paste(
        "Significant in the univariate screen (p = 0.0011, Figure 1C) but the",
        "forward inclusion step INCREASED -2LL by 31, and after CYP3A cluster",
        "inclusion by a further 133, so it was not retained. The Discussion",
        "attributes this to the unusually narrow haematocrit variability in a",
        "closely monitored inpatient cohort (CV 15% overall; 13.7% in PM and",
        "13.6% in IM). Cohort mean 31.9 +/- 5.0% (Table 1). No point estimate",
        "for a haematocrit effect is published, so none can be encoded.",
        sep = " "
      )
    ),
    SNP_PPARA_RS4253728 = list(
      description = "PPARA rs4253728 G>A variant, recessive coding (A/A vs G/A + G/G)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significant in the univariate screen (p = 0.007, Figure 1B) and",
        "improved the structural model by -2LL = 4 on disjoint forward",
        "inclusion, but after CYP3A cluster inclusion it increased -2LL by 184",
        "and was dropped. Only 4 of 59 patients were A/A homozygotes, 2 of whom",
        "were also CYP3A PM (Discussion). No point estimate is published.",
        "Column name is the shape-validated SNP_<GENE>_RS<rsid> form; it is NOT",
        "registered in inst/references/covariate-columns.md because the model",
        "does not use it.",
        sep = " "
      )
    ),
    SNP_POR_RS1057868 = list(
      description = "POR*28 (rs1057868 C>T) variant",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened univariately; not associated with any Bayesian posterior PK",
        "parameter. Cohort *1/*1 36 (61.0%), *1/*28 20 (33.9%), *28/*28",
        "3 (5.1%) (Table 1). The Discussion attributes the failure to replicate",
        "published POR*28 effects to insufficient statistical power.",
        sep = " "
      )
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 3435C>T variant",
      units = "(binary)",
      type = "binary",
      notes = "Screened univariately (Covariate Selection); not significant."
    ),
    SNP_ABCB1_RS9282564 = list(
      description = "ABCB1 1199G>A variant",
      units = "(binary)",
      type = "binary",
      notes = "Screened univariately (Covariate Selection); not significant."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 59L,
    n_studies = 1L,
    age_range = "Mean 51.9 +/- 13.4 years (Table 1)",
    weight_range = "Mean 70.4 +/- 13.9 kg (Table 1)",
    sex_female_pct = 64.4,
    race_ethnicity = "Not reported (single-centre Belgian cohort)",
    disease_state = "De novo cadaveric renal transplant recipients, early post-transplant hospitalization period",
    renal_function = "Creatinine clearance (Cockcroft-Gault) 60.1 +/- 20.0 mL/min at the PK course (Table 1)",
    dose_range = paste(
      "Immediate-release oral tacrolimus, initial dose 0.10 mg/kg body weight",
      "twice daily, then adjusted to trough concentration; dose before the PK",
      "course 5.5 +/- 2.7 mg (Table 1)",
      sep = " "
    ),
    administration = "Oral, twice daily",
    co_medication = paste(
      "Mycophenolate mofetil (81%) or mycophenolate sodium (19%) plus steroids",
      "on a standard tapering schedule. 21% received a P-glycoprotein inhibitor",
      "(atorvastatin, proton-pump inhibitors) at reduced dosage; no CYP3A",
      "inducer or inhibitor was documented.",
      sep = " "
    ),
    genotype = paste(
      "CYP3A metabolizer clusters (Elens 2011 clustering, Results): 5 poor",
      "metabolizers (CYP3A5*3/*3 carrying CYP3A4*22), 36 intermediate",
      "metabolizers (CYP3A5*3/*3 not carrying CYP3A4*22) and 18 extensive",
      "metabolizers (CYP3A5 expressers not carrying CYP3A4*22).",
      sep = " "
    ),
    regions = "Belgium (single centre: Cliniques Universitaires Saint-Luc, Brussels)",
    notes = paste(
      "Prospectively recruited July 2007 - January 2009; the same cohort as",
      "Elens 2013 (Ther Drug Monit 35:608-616). Full 12-h PK profile in every",
      "patient before hospital discharge: pre-dose and 30 min, 1 h 30 min, 3, 4,",
      "8 and 12 h after the morning dose, plus daily troughs. Tacrolimus",
      "measured in whole blood by chemiluminescent microparticle immunoassay",
      "(Abbott Architect). Target trough 10-20 ng/mL in week 1, 10-15 ng/mL",
      "thereafter. Estimation was NON-PARAMETRIC (Pmetrics NPAG), so the",
      "published typical values are means of the non-parametric marginal",
      "support-point distributions, not the parameters of a 'typical subject'.",
      sep = " "
    )
  )

  # Implementation notes (the vignette 'Assumptions and deviations' section
  # carries the full justification for each item):
  #
  # * STRUCTURE. Woillard 2017 Supplemental Data 1 gives the absorption rate as
  #   a two-component gamma mixture,
  #     v_abs(t) = F*D*sum_i r_i*f_i(t),
  #     f_i(t)   = b_i^a_i / Gamma(a_i) * t^(a_i-1) * exp(-b_i*t),
  #   the disposition after a unit IV bolus as I(t) = A_IV*exp(-k*t), and the
  #   analytic convolution of the two as
  #     C(t) = C0 + F*D*A_IV*exp(-k*t) *
  #            sum_i r_i * (b_i/(b_i-k))^a_i * P(a_i, (b_i-k)*t),
  #   where P is the regularised lower incomplete gamma function. The
  #   supplement states this closed form is "the output (C(t)) equation in our
  #   model file", so it is what this file evaluates, using rxode2's
  #   gammap(a, z) for P(a, z). Deriving the same expression from
  #   d/dt(central) = v_abs(t) - k*central with central(0) = 0 is an exact
  #   convolution identity, so the ODE and the closed form are the same model.
  #   The vignette gates the rxode2 solve against an independent R pgamma()
  #   evaluation AND against a numerical convolution of the published
  #   absorption rate with the published disposition function.
  #
  # * PARAMETERISATION. The gamma components are stored in the canonical
  #   transit parameterisation (ntr, mtt) rather than the paper's (a, b),
  #   following Debord_2001_cyclosporin.R and Fromage_2025_mycophenolic_acid.R
  #   (the same Limoges gamma-absorption lineage): a_i = ntr_i + 1 and
  #   b_i = a_i/mtt_i, a bijection, so mtt_i = a_i/b_i is exactly the mean
  #   absorption time of component i (0.606 h and 3.008 h here).
  #
  # * F*A_IV -> vc. Table 2 reports the product F*A_IV = 24.52 (the
  #   bioavailability coefficient times the concentration reached after a bolus
  #   of the model's unit dose); no independent estimate of F was possible
  #   because no IV data were available (Supplemental Data 1). With
  #   concentrations in ng/mL (= ug/L) and the unit dose equal to 1 mg
  #   (= 1000 ug; see the C0 note below), the amplitude F*A_IV*D equals
  #   1000*D/vc with the apparent volume vc = 1000/24.52 = 40.78 L, which is a
  #   physiologically sensible whole-blood V/F for tacrolimus. vc is therefore
  #   an apparent (oral) volume, V/F, not an absolute one.
  #
  # * C0 -> lrbase, AND ITS DOSE SCALING. The Table 2 footnote defines C0 as
  #   "the model estimated Tac trough level for a theoretical dose of 1000 mg
  #   (the real trough level can be calculated by dividing this value by 1000
  #   and multiplying by the patient dose)", i.e. C0 is DOSE-PROPORTIONAL, and
  #   Supplemental Data 1 adds that "as the patients were already at
  #   pharmacokinetic steady state, C0 corresponds to the trough concentration
  #   before the input dose". The stated unit "1000 mg" is a units slip for
  #   1000 ug = 1 mg: read literally it puts the typical trough at
  #   2.94/1000 * 5.5 = 0.016 ng/mL, against observed troughs of
  #   11.3 +/- 4.2 ng/mL (Table 1) and against the paper's own Table 3, where a
  #   7.5 mg dose in a poor metabolizer attains a trough >= 10 ng/mL 84.9% of
  #   the time. Reading the reference dose as 1000 ug gives a typical trough of
  #   C0 * D[mg] = 2.94 * 5.5 = 16.2 ng/mL at the cohort mean dose and a
  #   typical-to-median ratio across every cell of Table 3 of ~1.25-1.4, which
  #   is exactly what a right-skewed non-parametric distribution with the
  #   reported 40-80% CVs gives. It is also the only reading consistent with
  #   F*A_IV: the same unit dose has to scale both terms, and F*A_IV per ug
  #   would put V/F at 41 mL. lrbase is therefore the trough concentration PER
  #   mg of dose (ng/mL per mg) and enters as rbase*podo().
  #
  # * CYP3A CLUSTER. Results ('Covariate Analysis'): "The theta_CYP3A parameter
  #   was ascribed to the final output of the model (i.e., the Tac blood
  #   concentrations) in the form C(t) = C(t)_TPV * (theta_CYP3A)^CYP3A", with
  #   the cluster "encoded as a dummy variable". The dummy is the ordinal
  #   PM = 0 / IM = 1 / EM = 2, which the Abstract's own effect sizes confirm
  #   twice over: 1 - 0.77 = 23% lower in IM with CI 1 - [0.80, 0.74] =
  #   [20%, 26%] (printed "IC95%[20-26%]"), and 1 - 0.77^2 = 40.7% lower in EM
  #   with CI 1 - [0.80, 0.74]^2 = [36%, 45%] (printed "IC95%[36-45%]"). Both
  #   intervals reproduce exactly. The Abstract's "33%" for IM is a typo for
  #   23%: its own quoted interval [20-26%] does not contain 33%.
  #   The cluster is reconstructed from the two canonical genotype columns
  #   rather than registered as a collapsed three-level column, per the
  #   explicit instruction in the SNP_CYP3A4_RS35599367 register entry and
  #   following MohammedAli_2023_tacrolimus.R, which reconstructs the same
  #   Elens 2011 cluster from the same two inputs.
  #
  # * NO IIV IS ENCODED. Results reports only that "Inter-patient variability
  #   in PK parameters was represented by coefficients of variation ranging
  #   from 40 to 80% whereas the correlation between parameters fluctuated from
  #   r = -0.497 to 0.410" -- a range across an unnamed set of parameters, with
  #   no per-parameter variance and no covariance matrix. NPAG estimates a
  #   discrete joint support-point distribution, which is not published either.
  #   Inventing per-parameter omegas from the 40-80% range is not permitted, so
  #   this file carries typical values only. The etas are OMITTED rather than
  #   written as `~ fixed(0)` because a zero-variance diagonal makes OMEGA
  #   singular and breaks the Cholesky sampler used by rxSolve (the same
  #   reasoning as Thoueille_2026_salmeterol.R).
  #
  # * RESIDUAL ERROR. Methods gives the Pmetrics assay error polynomial
  #   SD = 0.0001 + 0.0762*C(t) - 0.1433*C(t)^2 with a fitted multiplier
  #   gamma = 0.43 for the final model, total noise = gamma*SD. The quadratic
  #   coefficient cannot be used as printed: -0.1433*C^2 drives SD negative
  #   above C = 0.53 ng/mL, i.e. over essentially the whole observed range
  #   (troughs alone are ~11 ng/mL), so the printed coefficient has lost a
  #   power of ten that the paper does not supply. The encoded error therefore
  #   keeps the two terms that are usable -- proportional 0.43*0.0762 = 3.28%
  #   and additive 0.43*0.0001 = 4.3e-5 ng/mL -- and drops the quadratic. The
  #   result is corroborated by the paper's own fit diagnostics: mean bias
  #   -0.11 +/- 3.7% and RMSE 4.5% (Results), both consistent with a ~3.3%
  #   proportional residual and inconsistent with any materially larger one.
  #
  # * SINGLE-INTERVAL SCOPE. tad() and podo() refer to the most recent dose, so
  #   each dose restarts the absorption input and there is no superposition.
  #   This matches the authors' framing (C0 IS the steady-state trough carried
  #   in from previous doses), but it means the model describes ONE
  #   steady-state dosing interval and must not be used to build up
  #   accumulation across doses.
  ini({
    # Double-gamma absorption, component 1 (the fast route). Stored as
    # (ntr, mtt); a1 = ntr1 + 1 and b1 = a1/mtt1 recover Table 2 exactly.
    lntr1 <- log(12.33 - 1)
    label("Gamma component 1 shape minus 1 (dimensionless)") # Woillard 2017 Table 2 final model: a1 = 12.33 [CI95% 6.25-18.41]; ntr1 = a1 - 1 = 11.33
    lmtt1 <- log(12.33 / 20.36)
    label("Gamma component 1 mean absorption time MAT1 (h)") # Woillard 2017 Table 2 final model: a1 = 12.33, b1 = 20.36 [CI95% 7.37-33.35] 1/h; MAT1 = a1/b1 = 0.6056 h

    # Double-gamma absorption, component 2 (the slow route).
    lntr2 <- log(15.19 - 1)
    label("Gamma component 2 shape minus 1 (dimensionless)") # Woillard 2017 Table 2 final model: a2 = 15.19 [CI95% 9.46-20.91]; ntr2 = a2 - 1 = 14.19
    lmtt2 <- log(15.19 / 5.05)
    label("Gamma component 2 mean absorption time MAT2 (h)") # Woillard 2017 Table 2 final model: a2 = 15.19, b2 = 5.05 [CI95% 1.02-9.08] 1/h; MAT2 = a2/b2 = 3.008 h

    # Fraction of the dose absorbed through the fast (component 1) route.
    lfdepot <- log(0.46)
    label("Fraction of dose absorbed via the fast gamma route (dimensionless)") # Woillard 2017 Table 2 final model: r = 0.46 [CI95% 0.40-0.51]

    # Disposition.
    lvc <- log(1000 / 24.52)
    label("Apparent central volume V/F (L)") # Woillard 2017 Table 2 final model: F*AIV = 24.52 [CI95% 20.61-28.43] ng/mL per mg; V/F = 1000/24.52 = 40.78 L (see the 'F*A_IV -> vc' implementation note)
    lkel <- log(1.52)
    label("First-order elimination rate constant (1/h)") # Woillard 2017 Table 2 final model: alpha = 1.52 [CI95% 1.19-1.85] 1/h (table footnote: 'alpha = elimination parameter')

    # Dose-proportional steady-state trough offset (the paper's C0).
    lrbase <- log(2.94)
    label("Steady-state pre-dose trough concentration C0 per mg of dose (ng/mL per mg)") # Woillard 2017 Table 2 final model: C0 = 2.94 [CI95% 2.42-3.47] (see the 'C0 -> lrbase' implementation note for the reference-dose unit)

    # CYP3A metabolizer cluster effect on the model output, applied as a power
    # of the ordinal cluster level (PM = 0, IM = 1, EM = 2).
    e_cyp3a_cluster_cc <- 0.77
    label("CYP3A metabolizer cluster multiplier per cluster level on predicted concentration") # Woillard 2017 Table 2 final model: theta_CYP3A = 0.77 [CI95% 0.74-0.80]

    # Residual error: Pmetrics assay polynomial scaled by the fitted gamma
    # multiplier; see the 'RESIDUAL ERROR' implementation note.
    propSd <- 0.0328
    label("Proportional residual SD (fraction)") # Woillard 2017 Methods 'Pharmacokinetic Population Modeling' (assay SD polynomial linear term 0.0762) x Results gamma = 0.43; 0.43 * 0.0762 = 0.03277
    addSd <- 4.3e-5
    label("Additive residual SD (ng/mL)") # Woillard 2017 Methods 'Pharmacokinetic Population Modeling' (assay SD polynomial constant term 0.0001) x Results gamma = 0.43; 0.43 * 0.0001 = 4.3e-5
  })
  model({
    # Individual parameters (typical values only; see the 'NO IIV IS ENCODED'
    # implementation note).
    ntr1 <- exp(lntr1)
    mtt1 <- exp(lmtt1)
    ntr2 <- exp(lntr2)
    mtt2 <- exp(lmtt2)
    fdepot <- exp(lfdepot)
    vc <- exp(lvc)
    kel <- exp(lkel)
    rbase <- exp(lrbase)

    # Three-level CYP3A metabolizer cluster (Elens 2011 clustering, as used in
    # Woillard 2017 Results), reconstructed from the two canonical genotype
    # columns exactly as MohammedAli_2023_tacrolimus.R does:
    #   EM = CYP3A5 expresser not carrying CYP3A4*22
    #   PM = CYP3A5*3/*3 carrying CYP3A4*22
    #   IM = everything else (CYP3A5*3/*3 not carrying CYP3A4*22)
    # The paper's dummy variable is the ordinal level PM = 0, IM = 1, EM = 2;
    # written as a product of complements so it stays exact 0/1 arithmetic.
    isEm <- (1 - SNP_CYP3A4_RS35599367) * CYP3A5_EXPR
    isPm <- SNP_CYP3A4_RS35599367 * (1 - CYP3A5_EXPR)
    isIm <- 1 - isEm - isPm
    cyp3aCluster <- 2 * isEm + isIm

    # Gamma shape / rate parameters of the two absorption routes
    # (Woillard 2017 Supplemental Data 1). a_i = ntr_i + 1, b_i = a_i/mtt_i.
    ga1 <- ntr1 + 1
    gb1 <- ga1 / mtt1
    ga2 <- ntr2 + 1
    gb2 <- ga2 / mtt2

    # Dose-registration state; always empty (see compartmentData notes).
    d/dt(depot) <- -kel * depot
    f(depot) <- 0

    # Published closed form (Woillard 2017 Supplemental Data 1). gammap(a, z)
    # is the regularised lower incomplete gamma function P(a, z). ugPerMg
    # converts the mg dose to ug so that dose/vc lands in ug/L == ng/mL, the
    # units the paper reports concentrations in.
    ugPerMg <- 1000
    tad1 <- tad()
    dose1 <- podo()
    arm1 <- fdepot * (gb1 / (gb1 - kel))^ga1 * gammap(ga1, (gb1 - kel) * tad1)
    arm2 <- (1 - fdepot) * (gb2 / (gb2 - kel))^ga2 * gammap(ga2, (gb2 - kel) * tad1)

    # The CYP3A effect multiplies the whole output, C(t) = C(t)_TPV *
    # theta_CYP3A^CYP3A, not any individual PK parameter (Results, 'Covariate
    # Analysis').
    Cc <- (rbase * dose1 +
      ugPerMg * dose1 / vc * (arm1 + arm2) * exp(-kel * tad1)) *
      e_cyp3a_cluster_cc^cyp3aCluster
    Cc ~ add(addSd) + prop(propSd)
  })
}
