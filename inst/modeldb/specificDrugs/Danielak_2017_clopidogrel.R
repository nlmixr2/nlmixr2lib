Danielak_2017_clopidogrel <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral clopidogrel and",
    "its active thiol H4 metabolite (the antiplatelet-active diastereomer)",
    "in adult Caucasian patients undergoing elective coronarography or",
    "percutaneous coronary intervention on chronic clopidogrel 75 mg/day",
    "(Danielak 2017). Clopidogrel is described by a one-compartment model",
    "with first-order absorption (rate constant ka = source k12) and",
    "first-order elimination (CL/F = source CL/F, V/F = source V2/F).",
    "The H4 metabolite is described by a one-compartment model with",
    "irreversible first-order formation from clopidogrel central at the",
    "rate FM * CL/F * (clopidogrel central / Vc) and first-order",
    "elimination (CL_h4/F = source Q2/F, V_h4/F = source V3/F). FM was",
    "constrained to <= 20% in the source fit because clopidogrel undergoes",
    "extensive first-pass metabolism to the inactive carboxylic acid (the",
    "competing CES1 pathway accounts for ~85% of the absorbed dose); the",
    "final estimate is FM = 4.5%. H4 plasma concentrations were assayed",
    "after bromo-3'-methoxyacetophenone derivatisation of the labile thiol",
    "and were adjusted to the mass equivalent of clopidogrel, so the",
    "parent <-> H4 flux carries 1:1 molar / mass-equivalent stoichiometry.",
    "Inter-individual variability is reported on ka, V/F, CL/F, and FM",
    "with a covariance between ka and V/F. The only retained covariate is",
    "CYP2C19*2 carriage on FM (linear-deviation effect, e_cyp2c19_s2_fm =",
    "-0.45); carriers convert 45% less of the absorbed dose to the active",
    "H4 metabolite. Bioavailability F was assumed to be unity (typical",
    "value 1, not estimated because no IV clopidogrel data exist).",
    "Residual error is proportional on the linear-concentration scale for",
    "both observed analytes; M3-method handling was used for samples below",
    "the quantitation limit (0.25 ng/mL for both clopidogrel and H4).")
  reference <- "Danielak D, Karazniewicz-Lada M, Komosa A, Burchardt P, Lesiak M, Kruszyna L, Graczyk-Szuster A, Glowka F. Influence of genetic co-factors on the population pharmacokinetic model for clopidogrel and its active thiol metabolite. Eur J Clin Pharmacol. 2017;73(12):1623-1632. doi:10.1007/s00228-017-2334-z"
  vignette <- "Danielak_2017_clopidogrel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP2C19_S2_CARRIER = list(
      description        = "CYP2C19*2 loss-of-function allele carrier indicator: 1 = subject carries at least one *2 allele (heterozygous *1/*2 or homozygous *2/*2); 0 = no *2 allele.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP2C19*2 allele -- *1/*1, *1/*17, or *17/*17).",
      notes              = "Time-fixed germline genotype determined by PCR-RFLP for rs4244285 (Danielak 2017 Methods 'Determination of genetic polymorphisms'). Pooled het + hom because the Danielak 2017 cohort had no *2/*2 homozygous-poor-metabolizers (n = 0 PM in the n = 63 cohort; Table 1 phenotype distribution lists only UM, EM, and IM). Linear-deviation effect on FM via the source equation FM = TVFM * (1 + e_cyp2c19_s2_fm * CYP2C19_S2_CARRIER) -- carriers convert 45% less of the absorbed clopidogrel to the active H4 thiol metabolite, giving a 36.7% lower predicted AUC of H4 vs non-carriers (Danielak 2017 Results page 1628 and Figure 4b).",
      source_name        = "CYP2C19*2"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 63L,
    n_studies            = 1L,
    n_observations       = "155 clopidogrel and 158 H4 plasma concentrations; 17 BQL samples (10 clopidogrel, 7 H4) handled with the NONMEM M3 method (Beal 2001).",
    age_range            = "mean 65.4 +/- 10.5 years (Table 1).",
    weight_range         = "mean 79.1 +/- 14.0 kg (n = 55 with weight reported; missing weights for 1 female and 7 male subjects imputed with sex-specific means -- female mean 73.5 kg, male mean 82.3 kg per Danielak 2017 Results 'Patients' characteristics').",
    sex_female_pct       = 33.3,
    race_ethnicity       = c(White = 100),
    disease_state        = "Stable adult patients undergoing elective coronarography or percutaneous coronary intervention (PCI) on chronic oral clopidogrel 75 mg once daily for at least 7 days prior to the procedure. Exclusion criteria: acute myocardial infarction, treatment with glycoprotein IIb/IIIa antagonists / coumarin derivatives / antiplatelet drugs other than aspirin, platelet count < 100,000/uL, ongoing malignancy, liver dysfunction, impaired renal function (serum creatinine > 2 mg/dL). 95% on concomitant statin; 73% on concomitant PPI (mostly pantoprazole).",
    dose_range           = "Oral clopidogrel 75 mg once daily for at least 7 days before sampling. Plasma sampling on the day of the procedure: full PK profile (0.5, 1, 2, 3, 4 h post-dose) in n = 17; sparse sampling (0.5 + 2 h or 1 + 3 h post-dose) in n = 46. The short sampling window (<= 4 h) was chosen because both clopidogrel and H4 are rapidly eliminated and concentrations are mostly BQL beyond 4 h post-dose.",
    regions              = "Single-centre, Poznan, Poland.",
    genotype_distribution = "CYP2C19*2 carriers n = 20 (31.7%); CYP2C19*17 carriers n = 34 (54.0%); CYP3A4*1G carriers n = 10 (15.8%); ABCB1 3435TT homozygotes n = 22 (34.9%). Phenotype: UM (*1/*17 or *17/*17) n = 25 (39.7%), EM (*1/*1) n = 18 (28.6%), IM (*1/*2) n = 20 (31.7%); no PM (*2/*2) observed.",
    bmi_range            = "mean 28.40 +/- 4.89 kg/m^2; 37/55 (67.3%) classified as obese (BMI >= 30) -- Table 1.",
    co_medication        = "Statins 60/63 (95.2%); PPI 46/63 (73.0%); pantoprazole 43/46 (93.5% of PPI users), omeprazole 2/46, esomeprazole 1/46. Diabetes mellitus 22/63 (34.9%). PPI subtype not separately modelled because pantoprazole dominated the cohort.",
    notes                = "Tested but not retained covariates: age, weight, BMI, obesity (BMI >= 30), sex, diabetes mellitus, PPI coadministration, statin coadministration, CYP2C19*17 allele, CYP3A4*1G allele, ABCB1 3435TT genotype. In the initial step-wise forward-selection step, sex on ka (k12 in the source) and V/F, ABCB1 3435TT on CL/F, and CYP2C19*17 on FM and V/F were significant individually; all were dropped after CYP2C19*2 was included on FM in the backward elimination (Danielak 2017 Results 'PK model' page 1626). Clopidogrel and H4 concentrations were quantified by validated HPLC-MS/MS within 0.25-5 ng/mL (clopidogrel) and 0.25-50 ng/mL (derivatised H4); within- and between-day precision below 19.9% and assay relative error below 16%. Plasma stabilised with 25 uL of 500 mM bromo-3'-methoxyacetophenone per 5 mL whole blood immediately on collection (Takahashi et al. derivatisation protocol)."
  )

  ini({
    # ------------------------------------------------------------------
    # CLOPIDOGREL (PARENT) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Danielak 2017 Table 2 'Final model' column. The source NONMEM
    # ADVAN5 parameterisation uses rate constants and apparent volumes
    # (k12 = absorption rate, V2/F = clopidogrel apparent central
    # volume, CL/F = clopidogrel apparent total clearance, FM =
    # fraction metabolised to H4 (irreversible), V3/F = H4 apparent
    # central volume, Q2/F = H4 apparent clearance; "Q2/F" is the
    # source's NONMEM-block label for the H4 elimination rate
    # constant * V3/F, i.e. it is functionally an apparent clearance
    # for the metabolite -- the paper's AUC equation AUC_H4 =
    # (F * D * FM) / CL_H4 confirms this). Bioavailability F was
    # fixed at the typical value of unity (Danielak 2017 page 1626:
    # "The relative bioavailability (F) of clopidogrel ... was
    # assumed a typical value of unity.").

    lka <- log(0.592)
    label("Clopidogrel first-order absorption rate constant (1/h)")             # Danielak 2017 Table 2 final-model k12 = 0.592 1/h (the source labels this 'k12' = depot -> central clopidogrel rate constant under NONMEM ADVAN5)

    lvc <- log(7660)
    label("Clopidogrel apparent central volume V/F (L)")                        # Danielak 2017 Table 2 final-model V2/F = 7660 L

    lcl <- log(14500)
    label("Clopidogrel apparent total clearance CL/F (L/h)")                    # Danielak 2017 Table 2 final-model CL/F = 14,500 L/h (total apparent clearance; sum of conversion-to-H4 plus competing carboxylic-acid esterase pathway)

    lfm <- log(0.045)
    label("Typical fraction of clopidogrel converted to the active H4 metabolite (TVFM, unitless)")  # Danielak 2017 Table 2 final-model FM = 0.045 (constrained <= 0.20 during estimation because ~85% of the absorbed dose is inactivated by CES1 to the carboxylic acid; Discussion page 1628)

    lfdepot <- fixed(log(1))
    label("Clopidogrel bioavailability into depot (F1, fixed at 1 -- absolute F not identifiable without IV clopidogrel data)")  # Danielak 2017 Methods page 1626: "The relative bioavailability (F) of clopidogrel ... was assumed a typical value of unity."

    # ------------------------------------------------------------------
    # H4 METABOLITE STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # H4 concentrations were derivatised with bromo-3'-
    # methoxyacetophenone and adjusted to the mass equivalent of the
    # parent compound, so the metabolite "central_h4" amount is
    # carried in clopidogrel mass-equivalent units and the
    # FM * clopidogrel-elimination-rate flux is 1:1 in mass-equivalent
    # space.

    lvc_h4 <- log(4.89)
    label("H4 metabolite apparent central volume V3/F (L)")                     # Danielak 2017 Table 2 final-model V3/F = 4.89 L

    lcl_h4 <- log(252)
    label("H4 metabolite apparent total clearance CL_h4/F (L/h)")               # Danielak 2017 Table 2 final-model Q2/F = 252 L/h (source labels this Q2/F under the NONMEM ADVAN5 micro-constant convention, but the AUC_H4 = F*D*FM/CL_h4 relationship in Methods page 1626 confirms it is functionally the H4 apparent clearance)

    # ------------------------------------------------------------------
    # COVARIATE EFFECTS
    # ------------------------------------------------------------------
    # The only retained covariate is CYP2C19*2 carriage on FM, applied
    # as a linear-deviation multiplier: FM_i = TVFM * (1 + COV *
    # CYP2C19*2_i). Danielak 2017 Table 2 final-model footnote b:
    # "FM = TVFM * (1 + CYP2C19 * COV)" with COV = -0.45.

    e_cyp2c19_s2_fm <- -0.45
    label("Linear-deviation coefficient for CYP2C19*2 carriage on FM (unitless)")  # Danielak 2017 Table 2 final-model 'Effect of CYP2C19*2 on FM (COV)' = -0.45 (carriers convert 45% less of the absorbed clopidogrel to active H4)

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Danielak 2017 Table 2 'Final model' column reports IIV as
    # %CV via footnote a: %CV = (SQRT(EXP(OMEGA(N)) - 1)) * 100%.
    # Inverting this gives omega^2 = log(1 + (CV/100)^2):
    #   ka  (k12) CV 25.3% -> log(1 + 0.253^2) = 0.0620
    #   V/F (V2/F) CV 63.7% -> log(1 + 0.637^2) = 0.3407
    #   CL/F      CV 49.8% -> log(1 + 0.498^2) = 0.2215
    #   FM        CV 71.4% -> log(1 + 0.714^2) = 0.4121
    # The off-diagonal covariance between k12 and V2/F is reported
    # directly as the OMEGA matrix element (not as a correlation):
    #   cov(eta_k12, eta_V2/F) = -0.087.
    # All other off-diagonals are zero (the source only retained
    # the ka-V/F covariance during covariate selection).

    etalka + etalvc ~ c(0.0620,
                        -0.087, 0.3407)                                          # Danielak 2017 Table 2 final-model: IIV(k12) = 25.3% CV -> omega^2 = 0.0620; cov(k12, V2/F) = -0.087; IIV(V2/F) = 63.7% CV -> omega^2 = 0.3407
    etalcl ~ 0.2215                                                              # Danielak 2017 Table 2 final-model: IIV(CL/F) = 49.8% CV -> omega^2 = 0.2215
    etalfm ~ 0.4121                                                              # Danielak 2017 Table 2 final-model: IIV(FM) = 71.4% CV -> omega^2 = 0.4121

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Danielak 2017 reports proportional residual error for both
    # clopidogrel and H4 (Methods page 1626: "Additive, proportional,
    # and combined error models ... were tested" and Results page 1626:
    # "a proportional error model, with separate terms for both
    # entities, was applied"). Table 2 final-model column lists
    # 'Proportional residual error for CLP = -0.45' and
    # 'Proportional residual error for H4 = -0.66'. The negative
    # signs are a NONMEM parameterisation artifact -- the source
    # ADVAN5 control stream codes residual error as
    # Y = F * (1 + THETA(N) * EPS(1)) with EPS(1) ~ N(0, 1) and
    # THETA(N) unconstrained, so the optimizer converges to either
    # sign with identical likelihood (|THETA(N)| * EPS = -|THETA(N)|
    # * EPS in distribution when EPS is symmetric). The magnitude
    # |THETA(N)| is the proportional residual SD on the linear-
    # concentration scale: 45% CV for clopidogrel and 66% CV for H4
    # (see vignette Errata for the full sign-artifact rationale and
    # an alternative log-transformed-SD interpretation that gives a
    # similar magnitude). The reported %RSE on each parameter (4.5%
    # and 5.4%) is small, supporting the magnitudes as well-
    # identified despite the sign ambiguity.

    propSd <- 0.45
    label("Clopidogrel proportional residual SD on linear concentration (fraction)")     # Danielak 2017 Table 2 final-model: |Proportional residual error for CLP| = 0.45 (45% CV; reported as -0.45 in the source -- NONMEM THETA-on-EPS sign artifact, see vignette Errata)

    propSd_h4 <- 0.66
    label("H4 metabolite proportional residual SD on linear concentration (fraction)")   # Danielak 2017 Table 2 final-model: |Proportional residual error for H4| = 0.66 (66% CV; reported as -0.66 in the source -- NONMEM THETA-on-EPS sign artifact, see vignette Errata)
  })

  model({
    # ------------------------------------------------------------------
    # Individual PK parameters
    # ------------------------------------------------------------------
    # FM uses the linear-deviation covariate form from Danielak 2017
    # Table 2 footnote b: FM_i = TVFM_i * (1 + COV * CYP2C19_S2_CARRIER_i)
    # where TVFM_i = exp(lfm + etalfm) and COV = e_cyp2c19_s2_fm.

    ka       <- exp(lka + etalka)
    vc       <- exp(lvc + etalvc)
    cl       <- exp(lcl + etalcl)
    fm       <- exp(lfm + etalfm) *
                (1 + e_cyp2c19_s2_fm * CYP2C19_S2_CARRIER)
    vc_h4    <- exp(lvc_h4)
    cl_h4    <- exp(lcl_h4)

    # Micro-constants. Total clopidogrel elimination rate kel = CL/Vc
    # splits into (1) the conversion-to-H4 path at rate fm * kel and
    # (2) the competing carboxylic-acid path at rate (1 - fm) * kel.
    # Only the conversion-to-H4 flux enters central_h4; the loss to
    # the carboxylic acid is implicit in the parent's kel.
    kel      <- cl    / vc
    kel_h4   <- cl_h4 / vc_h4

    # ODE system. Clopidogrel is dosed into the `depot` compartment
    # (first-order absorption at rate ka). The H4 metabolite forms
    # irreversibly from central clopidogrel at rate fm * kel * central
    # and is eliminated linearly. H4 plasma concentrations were
    # measured in clopidogrel mass-equivalent units in the source
    # paper, so the 1:1 stoichiometric flux is exact in those units.
    d/dt(depot)      <- -ka     * depot
    d/dt(central)    <-  ka     * depot   - kel    * central
    d/dt(central_h4) <-  fm     * kel     * central - kel_h4 * central_h4

    # Bioavailability into the depot. F1 fixed at 1 because the
    # absolute parent F is non-identifiable without IV clopidogrel.
    f(depot) <- exp(lfdepot)

    # Observations. Dose units are mg and volumes are L, so
    # central / vc returns concentrations in mg/L; multiply by 1000
    # to express in the source paper's units of ng/mL (Danielak 2017
    # Table 2 row units 'V2/F [L]', 'CL/F [L/h]', and AUC units of
    # ng h/mL).
    Cc    <- central    / vc    * 1000
    Cc_h4 <- central_h4 / vc_h4 * 1000

    Cc    ~ prop(propSd)
    Cc_h4 ~ prop(propSd_h4)
  })
}
