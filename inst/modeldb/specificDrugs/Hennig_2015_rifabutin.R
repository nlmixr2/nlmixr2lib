Hennig_2015_rifabutin <- function() {
  description <- "Two-compartment population pharmacokinetic model for rifabutin with simultaneous two-compartment metabolite (25-O-desacetyl rifabutin) modelling in 44 African HIV-infected adults with pulmonary tuberculosis on 300 mg daily oral rifabutin (Hennig 2015). Body weight allometrically scaled (a priori; CL exponent 0.75, V exponent 1) on all rifabutin apparent clearances and apparent volumes; sex effect on rifabutin V/F (males 1.84-fold higher than females); SLCO1B1 rs11045819 heterozygous-AC genotype increases rifabutin bioavailability F by 30.4 percent relative to homozygous-CC reference. Des-rifabutin parameters are apparent (with respect to rifabutin F and metabolite-formation fraction) and were estimated without allometric scaling, with metabolite Q and peripheral V fixed."
  reference <- "Hennig S, Naiker S, Reddy T, Egan D, Kellerman T, Wiesner L, Owen A, McIlleron H, Pym A. The effect of SLCO1B1 polymorphisms on the pharmacokinetics of rifabutin in African HIV-infected patients with tuberculosis. Antimicrob Agents Chemother. 2016 Jan;60(1):617-20. doi:10.1128/AAC.01195-15"
  vignette <- "Hennig_2015_rifabutin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of rifabutin apparent CL/F, Q/F, Cle/F, V/F, and Vp/F at reference 70 kg. Hennig 2015 Methods paragraph 3 (lines 100-102): allometric scaling fixed a priori per Anderson and Holford (Annu Rev Pharmacol Toxicol 2008;48:303-32; reference 20). Clearance exponent 0.75; volume exponent 1.0. Des-rifabutin parameters are NOT weight-scaled (Methods restricts scaling to rifabutin parameters; Table 2 metabolite-row units omit /70 kg).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Hennig 2015 Discussion (lines 130-131; line 175-179) and Table 2 covariate-effects block: males have a 1.84-fold higher rifabutin central V/F than females. Cohort: 27 male / 17 female (61 percent male) per Results paragraph 1 (lines 114-115). Effect coded as `(1 + e_sex_vc * (1 - SEXF))` so the female-baseline V/F equals exp(lvc) and males get the 1.84x factor.",
      source_name        = "SEX (paper text uses a male-deviation factor; the source NONMEM column name is in the unrecovered AAC supplement)"
    ),
    SNP_SLCO1B1_RS11045819 = list(
      description        = "SLCO1B1 rs11045819 binary genotype indicator (1 = at least one mutant A allele; 0 = homozygous CC reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type CC)",
      notes              = "Hennig 2015 Table 1: 5 of 35 successfully genotyped patients are AC heterozygotes; 30 of 35 are CC homozygotes; no AA homozygotes were observed. Effect on rifabutin bioavailability F: AC carriers have 30.4 percent higher F than CC reference (Hennig 2015 Table 2; Results paragraph 4 lines 137-138; dOFV = -6.5).",
      source_name        = "SLCO1B1 rs11045819 genotype (paper text)"
    )
  )

  population <- list(
    n_subjects     = 44,
    n_studies      = 1,
    age_range      = "32.7 (5.9) years (mean (SD))",
    age_median     = NULL,
    weight_range   = "60.7 (8.7) kg (mean (SD))",
    weight_median  = NULL,
    sex_female_pct = 39,
    race_ethnicity = c(Black_African = 100),
    disease_state  = "HIV-infected adults with microbiologically confirmed pulmonary tuberculosis (CD4 lymphocyte count 50-200 cells/mm^3; Karnofsky score 100; no grade 3-4 clinical or laboratory findings).",
    dose_range     = "Rifabutin 300 mg orally once daily, given for the last 2 weeks of the intensive phase of antituberculosis treatment and the first 2 weeks of the continuation phase, with concomitant standard doses of isoniazid, pyrazinamide, and ethambutol (no antiretroviral therapy at the time of PK sampling).",
    regions        = "South Africa (Durban; ANRS 12150a trial, ClinicalTrials.gov NCT00640887)",
    sampling_design = "After 4 weeks of daily rifabutin without ART, blood samples drawn following an overnight fast at predose (24 h after the previous dose) and at 2, 3, 4, 5, 6, 8, 12, and 24 h after the dose of interest. Standard hospital breakfast served > 2 h post dose. Total of 780 PK observations across 44 patients.",
    assay          = "LC-MS/MS for both rifabutin and 25-desacetyl rifabutin (calibration ranges 3.91-1000 ng/mL for rifabutin and 0.780-200 ng/mL for des-rifabutin). Inter-batch accuracy 99.1-109.0 percent; precision CV < 9.2 percent.",
    height_mean    = "159.6 (7.7) cm",
    bmi_mean       = "22.8 (3.3) kg/m^2",
    cd4_mean       = "126.1 (44.0) cells/mm^3",
    notes          = "All patients were of Black African ethnicity. Genetic samples were unavailable for 7 of 44 patients; rs4149032 genotyping was unsuccessful in 2 further patients. The full SLCO1B1 panel (rs4149032, rs2306283, rs4149056, rs11045819) was tested for covariate effects; only rs11045819 entered the final model (rs4149056 was excluded a priori because only 1 patient was a carrier). The AAC paper supplement (referenced as 'Supplementary material' for Fig. S1, S2, S3 and Table S1) was not available on disk during extraction; structural-model and parameter values were taken from the main paper Methods and Table 2."
  )

  ini({
    # ============================================================
    # Rifabutin (parent) - 2-compartment with first-order absorption
    # and absorption lag time. Body-weight allometric scaling
    # (a priori, exponent 0.75 on clearances and 1.0 on volumes,
    # reference 70 kg) is applied to CL/F, Q/F, Cle/F, V/F, and
    # Vp/F. Apparent clearance Cl/F here represents the
    # non-metabolic-formation arm of rifabutin elimination from the
    # central compartment; the formation clearance to the
    # 25-desacetyl rifabutin metabolite (lcl_form_desrbn) is a
    # separate parallel rifabutin-elimination arm so the total
    # rifabutin elimination clearance is (cl + cl_form_desrbn).
    # This matches the JAC 2016 Hennig pooled-rifabutin DDI control
    # stream parameterization (DADT(2) = ... - K20*A2 - K24*A2 with
    # K20 = CL/V and K24 = CLe/V).
    # ============================================================
    lka  <- log(0.24);   label("Rifabutin first-order absorption rate constant (1/h)")              # Hennig 2015 Table 2 (Final model row 'Absorption rate constant ka (1/h)')
    lcl  <- log(116.5);  label("Rifabutin apparent non-metabolic clearance Cl/F (L/h per 70 kg)")    # Hennig 2015 Table 2 (Final model row 'Clearance Cl/F (L/h/70 kg)')
    lvc  <- log(117.8);  label("Rifabutin apparent central volume V/F (L per 70 kg)")                # Hennig 2015 Table 2 (Final model row 'Central volume of distribution V/F (L/70kg)')
    lq   <- log(123.8);  label("Rifabutin apparent inter-compartmental clearance Q/F (L/h per 70 kg)") # Hennig 2015 Table 2 (Final model row 'Q/F (L/h/70 kg)')
    lvp  <- log(4897.8); label("Rifabutin apparent peripheral volume Vp/F (L per 70 kg)")             # Hennig 2015 Table 2 (Final model row 'Vpe/F (L/70kg)')
    lcl_form_desrbn <- log(21.2); label("Rifabutin-to-des-rifabutin apparent formation clearance Cle/F (L/h per 70 kg)") # Hennig 2015 Table 2 (Final model row 'Cle/F (metabolism of RBN to des-RBN)') ; allometrically weight-scaled per Hennig 2015 Methods/Results lines 100-102 and 127-128.
    ltlag <- log(1.6);   label("Rifabutin absorption lag time (h)")                                   # Hennig 2015 Table 2 (Final model row 'Lag time (h)')
    lfdepot <- fixed(log(1)); label("Rifabutin bioavailability F (typical value FIXED to 1)")         # Hennig 2015 Table 2 (Final model row 'Bioavailability F (Fixed)') ; fixed at population value 1 because absolute F is not identifiable from oral-only data.

    # Sex covariate effect on rifabutin V/F (males 1.84x higher than
    # females). Coded as 1 + e_sex_vc * (1 - SEXF) so the female
    # baseline equals exp(lvc) and males get the 1.84x factor.
    e_sex_vc <- 0.84; label("Fractional increase in rifabutin V/F for males relative to females (unitless; males +84 percent gives 1.84x factor)") # Hennig 2015 Discussion (lines 130-131); Table 2 'Increase of V/F for males (factor) = 1.8' (paper text reports 1.84 with an additional decimal place).

    # SLCO1B1 rs11045819 covariate effect on rifabutin F.
    # AC heterozygotes get a 1.304x factor; CC homozygotes (the
    # reference) get 1.0.
    e_snp_slco1b1_rs11045819_fdepot <- 0.304; label("Fractional increase in rifabutin bioavailability F for SLCO1B1 rs11045819 AC carriers relative to CC reference (unitless)") # Hennig 2015 Table 2 'Increase in bioavailability F (percent) for rs11045819 genotype = 30.4'.

    # ============================================================
    # 25-O-desacetyl rifabutin (metabolite, suffix '_desrbn')
    # 2-compartment with first-order elimination from the metabolite
    # central compartment. Metabolite parameters are apparent (i.e.,
    # with respect to the unidentifiable parent-F and
    # metabolite-formation-fraction product) and were NOT weight-
    # scaled in the source model. Qm/F and Vmper/F were FIXED in
    # Hennig 2015 (Table 2 marks both as 'Fixed').
    # ============================================================
    lcl_desrbn <- log(196.7); label("Des-rifabutin apparent elimination clearance Clm/F (L/h)")        # Hennig 2015 Table 2 (Final model row 'Clm/F (L/h)' under 'des-Rifabutin parameters')
    lvc_desrbn <- log(3.9);   label("Des-rifabutin apparent central volume Vm/F (L)")                  # Hennig 2015 Table 2 (Final model row 'Vm/F (L)' under 'des-Rifabutin parameters')
    lq_desrbn  <- fixed(log(0.15));   label("Des-rifabutin apparent inter-compartmental clearance Qm/F (L/h; FIXED)")  # Hennig 2015 Table 2 'Qm/F (L/h) (Fixed)'
    lvp_desrbn <- fixed(log(536.8));  label("Des-rifabutin apparent peripheral volume Vm-per/F (L; FIXED)")           # Hennig 2015 Table 2 'Vm-per/F (L) (Fixed)'

    # ============================================================
    # Inter-individual variability (between-subject variability,
    # BSV) - all reported as CV percent in Hennig 2015 Table 2; the
    # internal-scale variance is omega^2 = log(1 + CV^2).
    # ============================================================
    etalcl  ~ 0.01430                                                          # BSV 12.0 percent on Cl/F (Hennig 2015 Table 2 Final-model BSV column); omega^2 = log(1 + 0.12^2) = 0.01430
    etalvc  ~ 0.21528                                                          # BSV 49.0 percent on V/F (Hennig 2015 Table 2); omega^2 = log(1 + 0.49^2) = 0.21528
    etalka  ~ 0.05557                                                          # BSV 23.9 percent on ka (Hennig 2015 Table 2); omega^2 = log(1 + 0.239^2) = 0.05557
    etaltlag ~ 0.05923                                                         # BSV 24.7 percent on Lag time (Hennig 2015 Table 2); omega^2 = log(1 + 0.247^2) = 0.05923
    etalfdepot ~ 0.10336                                                       # BSV 33.0 percent on bioavailability F (Hennig 2015 Table 2 'F (Fixed)' BSV column); omega^2 = log(1 + 0.33^2) = 0.10336
    etalcl_desrbn ~ 0.08618                                                    # BSV 30.0 percent on des-rifabutin Clm/F (Hennig 2015 Table 2 des-Rifabutin BSV); omega^2 = log(1 + 0.30^2) = 0.08618

    # ============================================================
    # Residual error - combined additive + proportional on each
    # observation (Hennig 2015 Table 2 'Residual error' block).
    # Proportional component reported as CV percent; additive
    # component reported in ng/mL (the assay units).
    # ============================================================
    propSd <- 0.346; label("Rifabutin proportional residual SD (fraction)")              # Hennig 2015 Table 2 'Proportional error rifabutin (percent) = 34.6'
    addSd  <- 14.0;  label("Rifabutin additive residual SD (ng/mL)")                     # Hennig 2015 Table 2 'Additive error rifabutin (ng/mL) = 14.0'
    propSd_desrbn <- 0.346; label("Des-rifabutin proportional residual SD (fraction)")   # Hennig 2015 Table 2 'Proportional error des-rifabutin (percent) = 34.6'
    addSd_desrbn  <- 1.2;   label("Des-rifabutin additive residual SD (ng/mL)")          # Hennig 2015 Table 2 'Additive error des-rifabutin (ng/mL) = 1.2'
  })

  model({
    # ----------------------------------------------------------------
    # Body-weight allometric factors (a priori per Anderson & Holford;
    # exponents are fixed constants, not estimated). Reference 70 kg.
    # ----------------------------------------------------------------
    bw_cl <- (WT / 70) ^ 0.75
    bw_v  <- (WT / 70)

    # ----------------------------------------------------------------
    # Sex factor on rifabutin V/F (males 1.84x; females baseline).
    # SEXF = 1 (female) -> factor = 1; SEXF = 0 (male) -> factor = 1.84.
    # ----------------------------------------------------------------
    sex_vc_factor <- 1 + e_sex_vc * (1 - SEXF)

    # ----------------------------------------------------------------
    # SLCO1B1 rs11045819 factor on rifabutin F (AC carriers 1.304x).
    # ----------------------------------------------------------------
    snp_fdepot_factor <- 1 + e_snp_slco1b1_rs11045819_fdepot * SNP_SLCO1B1_RS11045819

    # ----------------------------------------------------------------
    # Individual rifabutin parameters (allometric weight-scaled).
    # ----------------------------------------------------------------
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)            * bw_cl
    vc  <- exp(lvc  + etalvc)            * bw_v * sex_vc_factor
    q   <- exp(lq)                        * bw_cl
    vp  <- exp(lvp)                       * bw_v
    cl_form_desrbn <- exp(lcl_form_desrbn) * bw_cl
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot + etalfdepot) * snp_fdepot_factor

    # ----------------------------------------------------------------
    # Individual des-rifabutin parameters (NOT allometric weight-scaled
    # per Hennig 2015 Methods/Results; Qm/F and Vmper/F FIXED).
    # ----------------------------------------------------------------
    cl_desrbn <- exp(lcl_desrbn + etalcl_desrbn)
    vc_desrbn <- exp(lvc_desrbn)
    q_desrbn  <- exp(lq_desrbn)
    vp_desrbn <- exp(lvp_desrbn)

    # ----------------------------------------------------------------
    # ODE system. Compartments: depot (oral dosing), central +
    # peripheral1 (rifabutin), central_desrbn + peripheral1_desrbn
    # (25-desacetyl rifabutin metabolite). The metabolite formation
    # term cl_form_desrbn / vc * central is a parallel elimination
    # arm on rifabutin central and a generation source on
    # central_desrbn (mass moves from parent into metabolite at rate
    # CLe * Cp_parent; the apparent metabolite parameters absorb the
    # unidentifiable molecular-weight ratio and metabolite-F factors).
    # ----------------------------------------------------------------
    d/dt(depot)              <- -ka * depot
    d/dt(central)            <-  ka * depot -
                                 (cl / vc) * central -
                                 (cl_form_desrbn / vc) * central -
                                 (q  / vc) * central +
                                 (q  / vp) * peripheral1
    d/dt(peripheral1)        <-  (q  / vc) * central -
                                 (q  / vp) * peripheral1

    d/dt(central_desrbn)     <-  (cl_form_desrbn / vc) * central -
                                 (cl_desrbn / vc_desrbn) * central_desrbn -
                                 (q_desrbn  / vc_desrbn) * central_desrbn +
                                 (q_desrbn  / vp_desrbn) * peripheral1_desrbn
    d/dt(peripheral1_desrbn) <-  (q_desrbn  / vc_desrbn) * central_desrbn -
                                 (q_desrbn  / vp_desrbn) * peripheral1_desrbn

    # ----------------------------------------------------------------
    # Bioavailability and absorption lag on the depot compartment.
    # ----------------------------------------------------------------
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ----------------------------------------------------------------
    # Observations. Internally central / vc has units mg/L = ug/mL
    # when dose is in mg and vc is in L; multiply by 1000 to convert
    # to ng/mL so that the proportional and additive residual-error
    # parameters (reported in ng/mL by Hennig 2015 Table 2) match the
    # output Cc and Cc_desrbn directly.
    # ----------------------------------------------------------------
    Cc        <- (central        / vc)        * 1000
    Cc_desrbn <- (central_desrbn / vc_desrbn) * 1000

    Cc        ~ add(addSd)        + prop(propSd)
    Cc_desrbn ~ add(addSd_desrbn) + prop(propSd_desrbn)
  })
}
