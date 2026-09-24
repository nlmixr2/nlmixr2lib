Lee_2019_tramadol <- function() {
  description <- "Joint parent-and-metabolite population PK model for sustained-release oral tramadol and its active CYP2D6-derived metabolite O-desmethyltramadol (M1) in healthy Korean male volunteers (Lee 2019). One-compartment disposition for tramadol with two PARALLEL absorption inputs that together deliver one dose: a fraction Fr enters a depot and is absorbed first-order at ka, and the complementary fraction 1 - Fr enters the central compartment as a lagged zero-order input of duration D2, which is how the authors captured the bimodal absorption phase of the extended-release formulation. Tramadol leaves the central compartment by two parallel first-order routes -- a non-M1 elimination clearance CL/F and a formation clearance CLPM/F into a one-compartment M1 pool that is itself eliminated by CLM/F. The CYP2D6*10/*10 genotype lowers both the tramadol elimination clearance (by 35.1 percent) and the M1 formation clearance (by 52.8 percent) relative to the wild-type reference group."
  reference <- "Lee HM, Kim SH, Kwon KH, Kim SJ, Cho CK, Bae JW, Jang CG, Lee SY. Population pharmacokinetic analysis of tramadol and O-desmethyltramadol with genetic polymorphism of CYP2D6. Drug Des Devel Ther. 2019;13:1751-1761. doi:10.2147/DDDT.S199574"
  vignette <- "Lee_2019_tramadol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = paste(
        "Compartment 1 of the Figure 2 scheme ('Depot'). Receives the fraction",
        "Fr = 0.88 of each oral dose (the f(depot) <- ffo line) and empties into",
        "central at the first-order rate ka. This is the SLOW arm: ka = 0.095",
        "1/h corresponds to an absorption half-life of 7.3 h, consistent with the",
        "Tridol SR sustained-release formulation the subjects received.",
        sep = " "
      )
    ),
    central = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "Compartment 2 of the Figure 2 scheme ('Tramadol'). Receives the",
        "complementary fraction 1 - Fr = 0.12 of the dose DIRECTLY as a zero-order",
        "input (NONMEM D2 = 1.93 h) delayed by NONMEM ALAG2 = 1.63 h, plus the",
        "first-order inflow ka * depot. Loses drug by two parallel first-order",
        "routes: k20 = CL/F / Vp/F (non-M1 elimination) and k23 = CLPM/F / Vp/F",
        "(formation of M1). Cc = 1000 * central / vc is the plasma tramadol",
        "concentration in ng/mL.",
        sep = " "
      )
    ),
    central_m1 = list(
      analyte = "O-desmethyltramadol (M1)",
      units = "mg (tramadol-equivalent; see notes)",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "Compartment 3 of the Figure 2 scheme ('O-desmethyl tramadol'). Fed",
        "mole-for-mole -- more precisely, unit-for-unit -- by k23 * central and",
        "eliminated by k30 = CLM/F / Vm/F. The source NONMEM ADVAN6 model applies",
        "NO molecular-weight correction on the parent-to-metabolite transfer",
        "(tramadol 263.4 g/mol versus M1 249.4 g/mol), so the amount in this state",
        "is carried in tramadol mass units and the 14 g/mol discrepancy, together",
        "with the unknown bioavailability F and the unknown fraction of tramadol",
        "clearance that is CYP2D6-mediated, is absorbed into the APPARENT",
        "metabolite volume Vm/F and clearance CLM/F. Those two parameters are",
        "therefore not interpretable as physiological M1 quantities on their own;",
        "only their ratio (the M1 elimination rate constant) and the resulting",
        "concentration Cc_m1 = 1000 * central_m1 / vc_m1 are. This is a faithful",
        "reproduction of the published model, not a simplification -- see the",
        "vignette Assumptions and deviations.",
        sep = " "
      )
    )
  )

  covariateData <- list(
    CYP2D6_STAR10_HOM = list(
      description = "CYP2D6*10 (rs1065852) homozygous-mutant indicator; 1 = CYP2D6*10/*10, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2D6*wt/*wt, the wild-type group; n = 14 of the 22 modelled subjects)",
      notes = paste(
        "Time-fixed germline genotype, assayed by pyrosequencing of the CYP2D6*2 and",
        "CYP2D6*10 alleles plus a long-PCR assay for CYP2D6*5 and gene duplication",
        "(Lee 2019 Methods 'Genotype analysis'). The cohort held only two genotype",
        "strata after exclusions -- *wt/*wt (n = 14) and *10/*10 (n = 8) -- so a",
        "single binary indicator is sufficient and the paired canonical",
        "CYP2D6_STAR10_HET indicator is not used; no heterozygous *10 carrier was",
        "enrolled. The one CYP2D6*5/*5 subject was excluded from the modelling",
        "dataset, and subjects carrying *5 or a duplicated CYP2D6 gene were",
        "excluded from the study at screening.",
        "",
        "The indicator enters TWO parameters as a fractional reduction, exactly as",
        "printed in Lee 2019 Results 'Covariate analysis':",
        "CL/F = 16.9 * (1 - 0.351 * G) and CLPM/F = 4.11 * (1 - 0.528 * G).",
        "The Discussion confirms both arithmetically: CL/F 16.9 L/h (wild type)",
        "versus 11.0 L/h (*10/*10), and CLPM/F 4.11 L/h versus 1.94 L/h.",
        "Adding the two effects lowered the inter-individual variability of CL/F",
        "and CLPM/F by 11.5 percent and 27.9 percent relative to the base model",
        "(Lee 2019 Results 'Covariate analysis'; Table 3 model 4 is the final",
        "model). No other covariate -- age, body weight, height or BMI -- was",
        "retained.",
        sep = " "
      ),
      source_name = "G_CYP2D6*10/*10 (Lee 2019 Results 'Covariate analysis' equations); CYP2D6 genotype in Table 1"
    )
  )

  # Screened but not retained in the final model, so no coefficient exists to
  # encode. Lee 2019 Results 'Covariate analysis': "Other covariates including
  # age, body weight, height and BMI did not appear to affect any of the PK
  # parameters in the final model."
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by scatter plot and by generalized additive modelling in Xpose 4.3.2; not retained in the final model."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened by scatter plot and by generalized additive modelling in Xpose 4.3.2; not retained in the final model. The cohort weight range was narrow (57-90 kg)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened by scatter plot and by generalized additive modelling in Xpose 4.3.2; not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened by scatter plot and by generalized additive modelling in Xpose 4.3.2; not retained in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    age_range = "20-40 years (mean 24.8, SD 4.8; median 24)",
    age_median = "24 years",
    weight_range = "57-90 kg (mean 71.6, SD 8.9; median 72)",
    weight_median = "72 kg",
    sex_female_pct = 0,
    race_ethnicity = "Korean; all volunteers",
    disease_state = "Healthy male volunteers (open-label single-centre PK study; KNUH IRB 2016-08-005)",
    dose_range = "Tramadol hydrochloride sustained-release tablets (Tridol SR, Yuhan), 100 mg orally every 12 h for a total of 5 doses, taken with 150 mL of water",
    administration = "Oral, sustained-release tablet",
    regions = "Korea (single centre: Kyungpook National University Hospital Clinical Trial Centre, Daegu)",
    genotypes = "CYP2D6*wt/*wt 14 (60.9%), CYP2D6*10/*10 8 (34.8%), CYP2D6*5/*5 1 (4.3%) among the 23 enrolled; the single *5/*5 subject was excluded from the modelling dataset, leaving 22 subjects",
    notes = paste(
      "Height 163-187 cm (mean 177, SD 5.6); BMI 18.1-26.9 kg/m2 (mean 23.0,",
      "SD 2.5). Table 1 describes all 23 enrolled subjects; Table 2 reports the",
      "final model fitted to the 22 subjects that remain after excluding the one",
      "CYP2D6*5/*5 volunteer (Lee 2019 Methods 'Population PK model",
      "development'). Subjects carrying a CYP2D6*5 allele or a duplicated CYP2D6",
      "gene were excluded at screening.",
      "",
      "SAMPLING. Blood was drawn at 0 (pre-dose), 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8,",
      "10, 12, 24, 48 and 72 h after the LAST (fifth) administration, giving 328",
      "tramadol and 323 M1 plasma concentrations. Assay: LC-MS/MS with calibration",
      "ranges of 1-1000 ng/mL for tramadol and 1-500 ng/mL for M1.",
      "",
      "ESTIMATION. NONMEM 7.3 with the PREDPP ADVAN6 general non-linear",
      "subroutine, first-order conditional estimation with eta-epsilon",
      "interaction. Evaluated by goodness-of-fit plots (Figures 3 and 4), a",
      "1000-replicate non-parametric bootstrap and a 1000-replicate visual",
      "predictive check (Figure 5).",
      "",
      "APPARENT PARAMETERS. Only oral data were collected, so bioavailability is",
      "not identifiable and every clearance and volume in this model is an",
      "apparent (/F) quantity. The metabolite parameters are additionally",
      "conditioned on the unknown fraction of tramadol clearance that proceeds",
      "through CYP2D6 -- see the central_m1 compartmentData notes.",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Tramadol disposition. Lee 2019 Table 2 ('Population pharmacokinetic
    # parameters for tramadol and O-desmethyltramadol in 22 subjects including
    # the results of bootstrap validation (final model)'). All values are
    # apparent (/F).
    #
    # The Figure 2 scheme puts CL/F (k20) and CLPM/F (k23) as two PARALLEL
    # first-order routes out of the tramadol central compartment, so the TOTAL
    # apparent tramadol clearance is cl + cl_form_m1 = 16.9 + 4.11 = 21.0 L/h
    # in a wild-type subject. `lcl` is therefore the non-M1 elimination route,
    # not the total.
    # -----------------------------------------------------------------------
    lcl <- log(16.9); label("Apparent non-M1 elimination clearance of tramadol CL/F (L/h), CYP2D6*wt/*wt") # Table 2 row 'CL/F (L/hr)' = 16.9 (RSE 20.7%; bootstrap median 16.50, 95% CI 15.4-18.6)
    lvc <- log(59.9); label("Apparent central volume of distribution of tramadol Vp/F (L)") # Table 2 row 'Vp/F (L)' = 59.9 (RSE 42.4%; bootstrap median 58.50, 95% CI 46.9-72.9)

    # -----------------------------------------------------------------------
    # Absorption. Lee 2019 Results 'Pharmacokinetic analysis': "We applied
    # combined first- and zero-order absorption to catch the bimodal absorption
    # phase of the extended release formulation of tramadol." Figure 2 shows
    # the dose splitting at the input: F1 into the depot (first-order, ka) and
    # 1 - F1 into the tramadol compartment (zero-order over D2, delayed by
    # ALAG2). NONMEM's D2 and ALAG2 both name compartment 2, which the Figure 2
    # scheme identifies as the tramadol central compartment.
    # -----------------------------------------------------------------------
    lka <- log(0.095); label("First-order absorption rate constant, depot to central (1/h)") # Table 2 row 'ka (hr-1)' = 0.095 (RSE 25.6%; bootstrap median 0.09, 95% CI 0.083-0.109)
    ld1 <- log(1.93); label("Duration of the zero-order input into the central compartment, D2 (h)") # Table 2 row 'D2 (hr)' = 1.93 (RSE 48.8%; bootstrap median 1.94, 95% CI 1.51-2.41)
    ltlag <- log(1.63); label("Lag time before the zero-order input into the central compartment starts, ALAG2 (h)") # Table 2 row 'ALAG2 (hr)' = 1.63 (RSE 47.1%; bootstrap median 1.65, 95% CI 1.24-2.02)
    logitffo <- qlogis(0.88); label("Logit of the fraction of the dose absorbed by the first-order route, Fr (unitless)") # Table 2 row 'Fr' = 0.88 (RSE 7.2%; bootstrap median 0.83, 95% CI 0.84-0.94)

    # -----------------------------------------------------------------------
    # O-desmethyltramadol (M1). Lee 2019 Results: "The metabolite, M1, was well
    # described by the one-compartment model as an extension of the parent drug
    # (tramadol) model, with first-order elimination."
    # -----------------------------------------------------------------------
    lcl_form_m1 <- log(4.11); label("Apparent formation clearance of tramadol to M1 CLPM/F (L/h), CYP2D6*wt/*wt") # Table 2 row 'CLPM/F (L/hr)' = 4.11 (RSE 34%; bootstrap median 4.01, 95% CI 3.42-4.81)
    lcl_m1 <- log(15.8); label("Apparent elimination clearance of O-desmethyltramadol CLM/F (L/h)") # Table 2 row 'CLM/F (L/hr)' = 15.8 (RSE 66.4%; bootstrap median 15.30, 95% CI 9.4-21.2)
    lvc_m1 <- log(8.63); label("Apparent volume of distribution of O-desmethyltramadol Vm/F (L)") # Table 2 row 'Vm/F (L)' = 8.63 (RSE 67.4%; bootstrap median 8.61, 95% CI 6.12-11.53)

    # -----------------------------------------------------------------------
    # CYP2D6*10/*10 covariate effects. Both are FRACTIONAL REDUCTIONS applied
    # multiplicatively as (1 - theta * G), exactly as printed in Lee 2019
    # Results 'Covariate analysis'.
    # -----------------------------------------------------------------------
    e_cyp2d6_star10_hom_cl <- 0.351; label("Fractional reduction in CL/F for CYP2D6*10/*10 subjects (unitless)") # Table 2 row 'CL/F, CYP2D6*10/*10' = 0.351 (RSE 81.4%; bootstrap median 0.34, 95% CI 0.191-0.498); Results equation 'CL/F = 16.9 (1 - 0.351 G)'
    e_cyp2d6_star10_hom_cl_form_m1 <- 0.528; label("Fractional reduction in CLPM/F for CYP2D6*10/*10 subjects (unitless)") # Table 2 row 'CLPM/F, CYP2D6*10/*10' = 0.528 (RSE 34.1%; bootstrap median 0.52, 95% CI 0.439-0.610); Results equation 'CLPM/F = 4.11 (1 - 0.528 G)'

    # -----------------------------------------------------------------------
    # Inter-individual variability. Lee 2019 Methods states the exponential
    # error model P_i = P_TV * exp(eta_i), and Table 2 tabulates the OMEGA
    # VARIANCES directly ('Variance of ...'), so the values below are used
    # as-is with no CV-to-variance conversion. Table 2 has a single unsubscripted
    # 'Variance of V/F' row; the parent volume is the one the Methods section
    # names V/F (the metabolite volume is always written Vm/F), so the eta is
    # carried on the parent central volume. Reported eta shrinkages: CL/F 1%,
    # V/F 14%, CLPM/F 26%, D2 16%.
    # -----------------------------------------------------------------------
    etalcl ~ 0.059 # Table 2 row 'omega^2 CL/F' variance = 0.059 (RSE 18.2%, shrinkage 1%; bootstrap median 0.058, 95% CI 0.053-0.065)
    etalvc ~ 0.023 # Table 2 row 'omega^2 V/F' variance = 0.023 (RSE 84.3%, shrinkage 14%; bootstrap median 0.024, 95% CI 0.018-0.028)
    etalcl_form_m1 ~ 0.017 # Table 2 row 'omega^2 CLPM/F' variance = 0.017 (RSE 56.2%, shrinkage 26%; bootstrap median 0.019, 95% CI 0.014-0.020)
    etald1 ~ 0.161 # Table 2 row 'omega^2 D2' variance = 0.161 (RSE 66.7%, shrinkage 16%; bootstrap median 0.165, 95% CI 0.135-0.192)

    # -----------------------------------------------------------------------
    # Residual unexplained variability. Lee 2019 Methods gives the combined
    # form C_ij = Cpred_ij * (1 + eps_pro,ij) + eps_add,ij, which is nlmixr2's
    # prop() + add(). Table 2 tabulates the SDs directly ('SD of proportional
    # error', 'SD of additive error'), and Results confirms them in words:
    # "The residual error of parent drug was 13% and 2.99 ng/mL ... the
    # residual error for the metabolite was 10.9% and 1.06 ng/mL".
    # -----------------------------------------------------------------------
    propSd <- 0.13; label("Proportional residual error for tramadol (fraction)") # Table 2 row 'sigma pro,tra' = 0.13 (RSE 52.1%, shrinkage 6.7%; bootstrap median 0.14, 95% CI 0.09-0.157)
    addSd <- 2.99; label("Additive residual error for tramadol (ng/mL)") # Table 2 row 'sigma add,tra' = 2.99 (RSE 51.3%, shrinkage 6.7%; bootstrap median 3.01, 95% CI 2.18-3.63)
    propSd_m1 <- 0.109; label("Proportional residual error for O-desmethyltramadol (fraction)") # Table 2 row 'sigma pro,ODMT' = 0.109 (RSE 33.8%, shrinkage 6.0%; bootstrap median 0.107, 95% CI 0.09-0.119)
    addSd_m1 <- 1.06; label("Additive residual error for O-desmethyltramadol (ng/mL)") # Table 2 row 'sigma add,ODMT' = 1.06 (RSE 67.3%, shrinkage 6.0%; bootstrap median 1.05, 95% CI 0.686-1.425)
  })

  model({
    # --- 1. Absorption -----------------------------------------------------
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag)
    # Fr is bounded in (0, 1) and is reported on the natural scale; it is
    # carried on the logit scale here so that a re-fit cannot leave the
    # feasible region. No IIV was reported on Fr.
    ffo <- expit(logitffo)

    # --- 2. Disposition ----------------------------------------------------
    # Lee 2019 Results 'Covariate analysis':
    #   CL/F   = 16.9 * (1 - 0.351 * G_CYP2D6*10/*10)
    #   CLPM/F = 4.11 * (1 - 0.528 * G_CYP2D6*10/*10)
    cl <- exp(lcl + etalcl) * (1 - e_cyp2d6_star10_hom_cl * CYP2D6_STAR10_HOM)
    vc <- exp(lvc + etalvc)
    cl_form_m1 <- exp(lcl_form_m1 + etalcl_form_m1) *
      (1 - e_cyp2d6_star10_hom_cl_form_m1 * CYP2D6_STAR10_HOM)
    cl_m1 <- exp(lcl_m1)
    vc_m1 <- exp(lvc_m1)

    # --- 3. Micro-constants and ODEs ---------------------------------------
    # Figure 2: k20 = CL/F / Vp/F, k23 = CLPM/F / Vp/F, k30 = CLM/F / Vm/F.
    kel <- cl / vc
    k23 <- cl_form_m1 / vc
    kel_m1 <- cl_m1 / vc_m1

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k23 * central
    d/dt(central_m1) <- k23 * central - kel_m1 * central_m1

    # --- 4. Parallel input routing -----------------------------------------
    # One oral dose is recorded twice: the fraction `ffo` = Fr enters `depot`
    # as an un-lagged first-order input, and the complementary fraction
    # 1 - Fr enters `central` as a zero-order input of duration d1 starting
    # after tlag. The central dose record MUST carry rate = -2 for rxode2 to
    # honour the modelled duration; without it the record becomes a bolus and
    # dur(central) is silently ignored.
    f(depot) <- ffo
    f(central) <- 1 - ffo
    dur(central) <- d1
    alag(central) <- tlag

    # --- 5. Observation ----------------------------------------------------
    # Dose is in mg and the volumes are in L, so amount/volume is mg/L; x 1000
    # gives ng/mL (= ug/L), the unit Lee 2019 reports concentrations in.
    Cc <- 1000 * central / vc
    Cc_m1 <- 1000 * central_m1 / vc_m1
    Cc ~ prop(propSd) + add(addSd)
    Cc_m1 ~ prop(propSd_m1) + add(addSd_m1)
  })
}
