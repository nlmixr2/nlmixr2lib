# Joint one-compartment population PK model for oral aripiprazole and its
# active metabolite dehydroaripiprazole in 80 Korean psychiatric patients,
# quantifying the four-level CYP2D6 stratification (Group I/II/III/IV) on the
# parent apparent oral clearance (Kim 2008 Br J Clin Pharmacol 66(6):802-810;
# doi:10.1111/j.1365-2125.2008.03223.x).

Kim_2008_aripiprazole <- function() {
  description <- paste(
    "Joint one-compartment population PK model for oral aripiprazole and its",
    "active metabolite dehydroaripiprazole in 80 Korean psychiatric patients",
    "(Kim 2008). First-order absorption (Ka FIXED at 1.06 1/h per a prior",
    "popPK analysis; the sparse-sampling design could not identify Ka) into",
    "a single aripiprazole central compartment with first-order elimination,",
    "and a metabolite central compartment that receives the entire parent",
    "elimination flux (fm = 1 assumed for identifiability; metabolite CL and",
    "V are apparent values scaled by the unknown fm and reported as CL(m)/fm",
    "and V(m)/fm) with first-order elimination of dehydroaripiprazole.",
    "Covariate analysis retained CYP2D6 genetic polymorphisms as the only",
    "significant covariate on parent CL/F, with four genotype strata fit as",
    "independent typical-value clearances (Group I 3.15 L/h, Group II 2.66,",
    "Group III 2.27, Group IV 1.83); the paper rejected pooling Groups I+II",
    "+III into a single CYP2D6 extensive-metabolizer stratum (uniting them",
    "increased OFV by 15.8 points). Age, body weight, gender, and CYP3A5",
    "genetic polymorphisms were screened and not retained. Inter-individual",
    "variability is fit on CL/F (shared across strata), V/F, and CL(m)/fm",
    "with an estimated covariance between CL/F and CL(m)/fm (value not",
    "reported by the paper; see vignette Assumptions and deviations). A",
    "proportional residual-error model is used separately for aripiprazole",
    "and dehydroaripiprazole.",
    sep = " "
  )
  reference <- paste(
    "Kim J-R, Seo H-B, Cho J-Y, Kang D-H, Kim Y K, Bahk W-M, Yu K-S,",
    "Shin S-G, Kwon J S, Jang I-J (2008).",
    "Population pharmacokinetic modelling of aripiprazole and its active",
    "metabolite, dehydroaripiprazole, in psychiatric patients.",
    "Br J Clin Pharmacol 66(6):802-810.",
    "doi:10.1111/j.1365-2125.2008.03223.x.",
    sep = " "
  )
  vignette <- "Kim_2008_aripiprazole"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or poor metabolizer; both CYP2D6_PM and CYP2D6_EM = 0 indicates IM)",
      notes              = paste(
        "1 = subject is a CYP2D6 extensive metabolizer carrying at least one",
        "functional CYP2D6 allele (*1 or *2) without a multiplicated allele;",
        "0 otherwise. In Kim 2008 the EM stratum is further subdivided by the",
        "second allele's functional status into Groups I/II/III via the paired",
        "binary canonicals CYP2D6_EM_1FUNC_1PD and CYP2D6_EM_1FUNC_1NULL (both",
        "indicators 0 within EM identifies Group I, the EM-2-functional",
        "reference subset). The Kim 2008 cohort had no PMs and no",
        "ultra-rapid metabolizers, so CYP2D6_PM is implicitly 0 throughout.",
        "Cohort distribution per Kim 2008 Results 'Study subjects and",
        "samples': Group I 15/80 (18.8%), Group II 26/80 (32.5%),",
        "Group III 12/80 (15.0%), Group IV (IM) 27/80 (33.8%); the EM",
        "indicator is 1 for Groups I/II/III combined (53/80, 66.3%)."
      ),
      source_name        = "Group I, II, or III in the paper's CYP2D6 stratification (Kim 2008 Table 2)"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, intermediate, or ultrarapid metabolizer)",
      notes              = paste(
        "1 = subject is a CYP2D6 poor metabolizer; 0 otherwise. Carried",
        "alongside CYP2D6_EM so the paired-indicator three-level encoding",
        "(EM = 1; PM = 1; both = 0 -> IM) generalises to cohorts that",
        "contain PMs. The Kim 2008 cohort had no PMs (Results 'Study",
        "subjects and samples': 'There were no patients who carried a",
        "multiplicated CYP2D6 allele or only null alleles.'), so CYP2D6_PM",
        "is 0 for every subject and the typical CL/F for the (hypothetical)",
        "PM stratum cannot be estimated from this paper. In a simulation",
        "or re-fit cohort that includes PMs, the typical CL/F for the PM",
        "stratum must be supplied externally."
      ),
      source_name        = "Not used in the Kim 2008 cohort (no PM subjects); included for paired-indicator generality"
    ),
    CYP2D6_EM_1FUNC_1PD = list(
      description        = "CYP2D6 extensive-metabolizer with 1 functional + 1 partially-deficient allele indicator (Kim 2008 Group II)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Group-II stratum). When paired with CYP2D6_EM = 1, CYP2D6_PM = 0, CYP2D6_EM_1FUNC_1NULL = 0, identifies Group II.",
      notes              = paste(
        "1 = subject is a CYP2D6 extensive metabolizer carrying one fully",
        "functional allele (*1, *2) paired with one partially-deficient",
        "allele (*10, *41); 0 otherwise. The paired *10 (100 C>T, P34S) and",
        "*41 (-1584 C>G) alleles encode partially functional CYP2D6 in the",
        "Kim 2008 classification scheme. Cohort distribution per Kim 2008",
        "Table 2: Group II 26/80 subjects (32.5%) comprising *1/*10 (22),",
        "*1/*41 (1), and *2/*10 (3). Time-fixed germline genotype. In the",
        "final model this indicator multiplicatively selects the lcl_em_1f1pd",
        "typical-value clearance (2.66 L/h) instead of the Group I reference",
        "lcl_em_2func (3.15 L/h)."
      ),
      source_name        = "Group II (Kim 2008 Table 2)"
    ),
    CYP2D6_EM_1FUNC_1NULL = list(
      description        = "CYP2D6 extensive-metabolizer with 1 functional + 1 null allele indicator (Kim 2008 Group III)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Group-III stratum). When paired with CYP2D6_EM = 1, CYP2D6_PM = 0, CYP2D6_EM_1FUNC_1PD = 0, identifies Group III.",
      notes              = paste(
        "1 = subject is a CYP2D6 extensive metabolizer carrying one fully",
        "functional allele (*1, *2) paired with one null allele (*4, *5, *14,",
        "*36); 0 otherwise. The *4 (1846 G>A), *5 (gene deletion), *14",
        "(1758 G>A), and *36 (CYP2D6 -> CYP2D7P conversion in exon 9)",
        "alleles encode no CYP2D6 activity. Cohort distribution per Kim 2008",
        "Table 2: Group III 12/80 subjects (15.0%) comprising *1/*5 (6),",
        "*1/*14 (2), *1/*36 (1), *2/*5 (2), and *2/*36 (1). Time-fixed",
        "germline genotype. In the final model this indicator selects the",
        "lcl_em_1f1null typical-value clearance (2.27 L/h) instead of the",
        "Group I reference lcl_em_2func (3.15 L/h)."
      ),
      source_name        = "Group III (Kim 2008 Table 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened by stepwise forward addition (Kim 2008 Table 4) and not",
        "retained: dOFV = -0.23 on CL/F (p = 0.632) and -0.17 on CL(m)/fm",
        "(p = 0.680). Cohort range 17-61 years (median 34, mean 35.2 +/- 9.8;",
        "Kim 2008 Table 1, excluding one subject of unknown age)."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened by stepwise forward addition (Kim 2008 Table 4) on CL/F",
        "(dOFV = -1.75, p = 0.186), CL(m)/fm (dOFV = -2.84, p = 0.092), and",
        "V/F (dOFV = -1.33, p = 0.249) and not retained. Cohort range",
        "38.0-95.0 kg (median 61.9, mean 60.8 +/- 10.9; Kim 2008 Table 1).",
        "The paper notes (Discussion) that the relatively narrow weight",
        "range (CV 18%) limits the power to detect a weight effect; prior",
        "popPK reports for aripiprazole DID identify a body-weight effect."
      )
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened by stepwise forward addition (Kim 2008 Table 4) on CL/F",
        "(dOFV = -0.07, p = 0.791) and CL(m)/fm (dOFV = -0.29, p = 0.590)",
        "and not retained. Cohort distribution: 34 male / 46 female",
        "(42.5% male; Kim 2008 Table 1)."
      )
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (1 = carries at least one functional CYP3A5*1 allele)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened by stepwise forward addition (Kim 2008 Table 4) on CL/F",
        "(dOFV = -1.69, p = 0.194) and CL(m)/fm (dOFV = -0.55, p = 0.458)",
        "and not retained. The non-significant CYP3A5 finding is",
        "interpreted (Discussion) as evidence that aripiprazole is",
        "metabolized selectively by CYP3A4 rather than CYP3A5. Cohort",
        "distribution per Kim 2008 Table 2: 3 *1/*1 + 33 *1/*3 = 36",
        "expressers (45.0%) vs 44 *3/*3 non-expressers (55.0%)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 1L,
    n_sites        = 15L,
    n_observations = 141L,
    age_range      = "17-61 years (median 34, mean 35.2 +/- 9.8; Kim 2008 Table 1, excluding one subject of unknown age)",
    age_median     = "34 years",
    weight_range   = "38.0-95.0 kg (median 61.9, mean 60.8 +/- 10.9; Kim 2008 Table 1)",
    weight_median  = "61.9 kg",
    sex_female_pct = 57.5,
    race_ethnicity = c(Korean = 100),
    disease_state  = paste(
      "Psychiatric patients with schizophrenia, schizophreniform disorder,",
      "or schizoaffective disorder enrolled in a multicentre Phase IV",
      "clinical trial in Korea (Kim 2008 Methods 'Study population')."
    ),
    dose_range     = paste(
      "Oral aripiprazole 10-30 mg once daily at steady state (mean dose",
      "24.3 +/- 6.7 mg; Kim 2008 Table 1). Initial titration from 15 mg",
      "with up-titration over the first several days, then maintained at",
      "steady state."
    ),
    regions        = "Republic of Korea (30 sites enrolled, 15 sites contributed sampled patients).",
    cyp2d6_distribution = paste(
      "Group I 15/80 (18.8%, 2 functional alleles), Group II 26/80 (32.5%,",
      "1 functional + 1 partially deficient), Group III 12/80 (15.0%, 1",
      "functional + 1 null), Group IV 27/80 (33.8%, IM = 2 partially",
      "deficient or 1 partially deficient + 1 null). No CYP2D6 PMs (no",
      "subjects with 2 null alleles) and no ultra-rapid metabolizers (no",
      "subjects with a multiplicated CYP2D6 allele) were observed; Kim 2008",
      "Results 'Study subjects and samples' and Table 2."
    ),
    cyp3a5_distribution = paste(
      "CYP3A5 expressers (Group A): 36/80 (45.0%, comprising 3 *1/*1 and",
      "33 *1/*3). CYP3A5 non-expressers (Group B): 44/80 (55.0%, all",
      "*3/*3). Kim 2008 Table 2."
    ),
    notes          = paste(
      "Sparse-sampling pop PK design: 141 plasma samples from 80 patients",
      "(1-2 samples per patient) drawn at steady state, mostly trough and",
      "approximately 6 h post-dose (5 trough-only patients, 14",
      "post-dose-only patients, and the remaining patients contributing",
      "both). Mean observed aripiprazole concentration 445.6 +/- 220.9",
      "ng/mL; mean observed dehydroaripiprazole concentration 129.3 +/-",
      "60.4 ng/mL (Kim 2008 Table 1). Aripiprazole oral bioavailability",
      "F = 0.87 per the prescribing information (Kim 2008 reference [2]);",
      "the structural CL/F and V/F estimates are apparent values that",
      "embed this F."
    )
  )

  ini({
    # ========================================================================
    # Absorption: Ka FIXED at 1.06 1/h per Kim 2008 Pharmacokinetics section
    # ("Since the available data contained little information on the oral
    # absorption of aripiprazole, ka could not be estimated properly and was
    # fixed to 1.06 h-1 in accordance with previous population analysis [20]").
    # ========================================================================
    lka <- fixed(log(1.06))
    label("Absorption rate constant Ka (1/h, FIXED)")  # Kim 2008 Table 3 row 'k_a (h-1) = 1.06 (FIXED)'

    # ========================================================================
    # Apparent clearance of aripiprazole, CYP2D6 stratum-specific typical
    # values (Kim 2008 Table 3 final model). Inter-individual variability is
    # shared across strata via the single etalcl random effect. Group I (EM
    # with 2 functional alleles) is the reference; indicator-gated linear
    # combination in model() selects the active typical value per subject
    # using the paired binary canonicals CYP2D6_EM, CYP2D6_PM,
    # CYP2D6_EM_1FUNC_1PD, and CYP2D6_EM_1FUNC_1NULL.
    # ========================================================================
    lcl_em_2func   <- log(3.15)
    label("Apparent clearance CL/F of aripiprazole in CYP2D6 EMs with 2 functional alleles, Group I (L/h)")  # Kim 2008 Table 3 row 'CL/F I = 3.15 (SE 0.223)'
    lcl_em_1f1pd   <- log(2.66)
    label("Apparent clearance CL/F of aripiprazole in CYP2D6 EMs with 1 functional + 1 partially-deficient allele, Group II (L/h)")  # Kim 2008 Table 3 row 'CL/F II = 2.66 (SE 0.119)'
    lcl_em_1f1null <- log(2.27)
    label("Apparent clearance CL/F of aripiprazole in CYP2D6 EMs with 1 functional + 1 null allele, Group III (L/h)")  # Kim 2008 Table 3 row 'CL/F III = 2.27 (SE 0.123)'
    lcl_im         <- log(1.83)
    label("Apparent clearance CL/F of aripiprazole in CYP2D6 IMs, Group IV (L/h)")  # Kim 2008 Table 3 row 'CL/F IV = 1.83 (SE 0.0715)'

    # ========================================================================
    # Apparent volume of distribution of aripiprazole V/F. Single typical
    # value across all CYP2D6 strata; final model retained no covariate on
    # V/F (Kim 2008 Table 3 row 'V/F').
    # ========================================================================
    lvc <- log(193)
    label("Apparent volume of distribution V/F of aripiprazole (L)")  # Kim 2008 Table 3 row 'V/F (l) = 193 (SE 11.7)'

    # ========================================================================
    # Metabolite dehydroaripiprazole. fm (fraction of absorbed aripiprazole
    # converted to dehydroaripiprazole) is set to 1 for identifiability, so
    # CL(m) and V(m) are reported as the apparent ratios CL(m)/fm and
    # V(m)/fm (Kim 2008 Methods 'Pharmacokinetic model and data analysis').
    # ========================================================================
    lcl_dehyari <- log(8.02)
    label("Apparent metabolic clearance CL(m)/fm of dehydroaripiprazole (L/h)")  # Kim 2008 Table 3 row 'CL(m)/f_m (l h-1) = 8.02 (SE 0.317)'
    lvc_dehyari <- log(587)
    label("Apparent volume of distribution V(m)/fm of dehydroaripiprazole (L)")  # Kim 2008 Table 3 row 'V(m)/f_m (l) = 587 (SE 217)'; the paper notes (Discussion) that V(m)/fm in the final model is much smaller and less precise than the base-model estimate of 950 L and that a bootstrap would help refine it; see vignette Assumptions and deviations

    # ========================================================================
    # Inter-individual variability (NONMEM OMEGA, log-normal variance scale).
    # Kim 2008 Table 3 reports w^2 CL/F = 0.0930, w^2 V/F = 0.100, and
    # w^2 CL(m)/fm = 0.113 in the final model. The paper notes (Results
    # 'Pharmacokinetics') that the covariance between CL/F and CL(m)/fm was
    # included in the final model because the full var-cov matrix did not
    # converge, but the numeric covariance value is not reported. The
    # covariance is encoded here as 0 (independent etalcl and etalcl_dehyari)
    # with a deviation noted in the vignette Assumptions and deviations.
    # IIV on V(m)/fm and on Ka was not retained ("Sparse sampling design led
    # to the incorporation of interindividual variability (IIV) only for
    # CL/F, V/F, and CL(m)/fm").
    # ========================================================================
    etalcl         ~ 0.0930   # Kim 2008 Table 3 row 'w^2 CL/F = 0.0930 (SE 0.0174)'; CV% sqrt scale 30.5%
    etalvc         ~ 0.100    # Kim 2008 Table 3 row 'w^2 V/F = 0.100 (SE 0.0475)'; CV% sqrt scale 31.6%
    etalcl_dehyari ~ 0.113    # Kim 2008 Table 3 row 'w^2 CL(m)/f_m = 0.113 (SE 0.0220)'; CV% sqrt scale 33.6%

    # ========================================================================
    # Residual variability. Kim 2008 Methods 'Pharmacokinetic model and data
    # analysis' specifies a proportional error model on linear-scale plasma
    # concentrations: C_ij = C_hat_ij * (1 + eps_ij) with eps ~ N(0, sigma^2).
    # Table 3 reports sigma^2 (variance scale) for parent (s^2 = 0.00264) and
    # metabolite (s^2_m = 0.00507) in the final model. nlmixr2 propSd is on
    # the SD scale, so propSd = sqrt(sigma^2).
    # ========================================================================
    propSd         <- sqrt(0.00264)
    label("Proportional residual error for aripiprazole (fraction; SD = sqrt(NONMEM sigma^2))")  # Kim 2008 Table 3 row 's^2 = 0.00264 (SE 0.00139)'; CV% sqrt scale 5.1%; propSd = sqrt(0.00264) = 0.0514
    propSd_dehyari <- sqrt(0.00507)
    label("Proportional residual error for dehydroaripiprazole (fraction; SD = sqrt(NONMEM sigma^2))")  # Kim 2008 Table 3 row 's^2_m = 0.00507 (SE 0.000715)'; CV% sqrt scale 7.1%; propSd_dehyari = sqrt(0.00507) = 0.0712
  })

  model({
    # -----------------------------------------------------------------------
    # 1. CYP2D6 stratum-specific indicator gating. Group I (EM with 2
    # functional alleles) is the implicit reference inside the EM stratum;
    # Group IV (IM) is selected when both CYP2D6_EM and CYP2D6_PM are 0.
    # The paired binary canonicals CYP2D6_EM, CYP2D6_PM, CYP2D6_EM_1FUNC_1PD,
    # and CYP2D6_EM_1FUNC_1NULL encode all four Kim 2008 groups without loss:
    #   Group I:   EM=1, PM=0, EM_1FUNC_1PD=0, EM_1FUNC_1NULL=0
    #   Group II:  EM=1, PM=0, EM_1FUNC_1PD=1, EM_1FUNC_1NULL=0
    #   Group III: EM=1, PM=0, EM_1FUNC_1PD=0, EM_1FUNC_1NULL=1
    #   Group IV:  EM=0, PM=0  (IM, the implicit non-PM non-EM stratum)
    # The Kim 2008 cohort contained no PMs, so the PM branch is never
    # exercised by paper data; it is carried in the parameterisation for
    # completeness but is mapped to the Group IV (IM) typical CL/F so a
    # simulated PM subject does not divide by zero.
    # -----------------------------------------------------------------------
    ind_g1 <- CYP2D6_EM * (1 - CYP2D6_EM_1FUNC_1PD - CYP2D6_EM_1FUNC_1NULL)
    ind_g2 <- CYP2D6_EM * CYP2D6_EM_1FUNC_1PD
    ind_g3 <- CYP2D6_EM * CYP2D6_EM_1FUNC_1NULL
    ind_g4 <- 1 - CYP2D6_EM - CYP2D6_PM

    # -----------------------------------------------------------------------
    # 2. Individual PK parameters. The four CYP2D6 stratum-specific typical
    # CL/F values share a single inter-individual variability term etalcl;
    # the indicator-gated linear combination selects which typical-value
    # serves as the multiplicative anchor for the per-subject log-normal
    # eta. Only one indicator is 1 per subject, so the active typical CL/F
    # equals the corresponding stratum's exp(lcl_*).
    # -----------------------------------------------------------------------
    cl_typ <- exp(lcl_em_2func)   * ind_g1 +
              exp(lcl_em_1f1pd)   * ind_g2 +
              exp(lcl_em_1f1null) * ind_g3 +
              exp(lcl_im)         * ind_g4
    cl <- cl_typ * exp(etalcl)

    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    cl_dehyari <- exp(lcl_dehyari + etalcl_dehyari)
    vc_dehyari <- exp(lvc_dehyari)

    # -----------------------------------------------------------------------
    # 3. ODE system. dA(dose)/dt = -ka*A(dose); dA(p)/dt = ka*A(dose) -
    # (cl/vc)*A(p); dA(m)/dt = fm*(cl/vc)*A(p) - (cl_dehyari/vc_dehyari)*A(m)
    # with fm = 1 fixed (Kim 2008 Methods 'Pharmacokinetic model and data
    # analysis'). Mass-fraction conversion of aripiprazole to dehydroari-
    # piprazole is treated as 1:1 because the molecular weights of the parent
    # (448) and the dehydro metabolite (446) differ by only 0.4%.
    # -----------------------------------------------------------------------
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (cl / vc) * central
    d/dt(central_dehyari) <- (cl / vc) * central - (cl_dehyari / vc_dehyari) * central_dehyari

    # -----------------------------------------------------------------------
    # 4. Observations. Dose in mg, volumes in L: central/vc has units mg/L =
    # ug/mL; multiply by 1000 to obtain ng/mL to match Kim 2008 Table 1
    # reported observed concentrations (aripiprazole 55.5-1442 ng/mL,
    # dehydroaripiprazole 31.7-428 ng/mL).
    # -----------------------------------------------------------------------
    Cc         <- 1000 * central         / vc
    Cc_dehyari <- 1000 * central_dehyari / vc_dehyari

    Cc         ~ prop(propSd)
    Cc_dehyari ~ prop(propSd_dehyari)
  })
}
