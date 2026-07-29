Horita_2018_isoniazid <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption and linear elimination for oral isoniazid in Ghanaian children with active tuberculosis (Horita 2018); NAT2 slow-vs-nonslow acetylator phenotype on apparent oral clearance with separate typical-value clearances and separate IIV omegas; allometric weight scaling on CL/F and Q/F (fixed 0.75) and V1/F and V2/F (fixed 1.0) normalised to the cohort median 14.3 kg."
  reference <- "Horita Y, Alsultan A, Kwara A, Antwi S, Enimil A, Ortsin A, Dompreh A, Yang H, Wiesner L, Peloquin CA. Evaluation of the Adequacy of WHO Revised Dosages of the First-Line Antituberculosis Drugs in Children with Tuberculosis Using Population Pharmacokinetic Modeling and Simulations. Antimicrob Agents Chemother. 2018;62(9):e00008-18. doi:10.1128/AAC.00008-18"
  vignette <- "Horita_2018_isoniazid"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F and Q/F with fixed exponent 0.75 and on V1/F and V2/F with fixed exponent 1.0 (Horita 2018 Results 'INH' paragraph 1: 'The fixed exponents were 0.75 for CL/F and clearance between compartments (Q/F) and 1.0 for V1/F and V2/F'). Reference weight is the cohort median 14.3 kg; see the vignette Errata for the reference-weight derivation.",
      source_name        = "WT"
    ),
    NAT2_SLOW = list(
      description        = "NAT2 slow-acetylator phenotype indicator (1 = slow, 0 = intermediate or rapid (pooled as 'nonslow'))",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nonslow: intermediate or rapid acetylator pooled)",
      notes              = "Horita 2018 Results 'INH' paragraph 1: 'NAT2 genotype was detected as a significant covariate for clearance. In light of the findings obtained from NCA that there were no significant differences in t1/2, CL/F, and AUC0-8 between NAT2 fast and intermediate genotypes, those were combined and named the nonslow group.' Cohort distribution (Table 1): 51 slow, 50 intermediate, 12 rapid -- so 51/113 = 45.1% slow. NAT2 genotyping via TaqMan PCR for SNPs rs1801279, rs1801280, rs1799930, rs1799931 (Methods 'Arylamine N-acetyltransferase 2 genotyping'). The slow phenotype reduces typical-value CL/F (4.44 L/h slow vs 8.08 L/h nonslow); see model() for the selection mechanism.",
      source_name        = "NAT2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 113L,
    n_studies      = 1L,
    age_range      = "3 months to 14 years (median 5.00 years, IQR 2.17 to 8.25)",
    age_median     = "5.00 years",
    weight_range   = "5-30 kg (median 14.3, IQR 9.70 to 20.1)",
    weight_median  = "14.3 kg",
    sex_female_pct = 44.2,
    hiv_positive_pct = 52.2,
    nat2_distribution = c(slow = 51L, intermediate = 50L, rapid = 12L),
    nat2_slow_pct  = 45.1,
    disease_state  = "Ghanaian children with active tuberculosis (HIV-positive and HIV-negative). 21.2% under 2 years of age.",
    dose_range     = "Isoniazid 7-15 mg/kg orally daily (median 11.0 mg/kg, IQR 9.06-12.8). Administered as part of standard four-drug anti-TB regimen.",
    regions        = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes          = "Patients enrolled October 2012-August 2015. PK sampling after at least 4 weeks of anti-TB treatment (steady state). Blood samples at 0, 1, 2, 4, 8 h postdose. INH concentrations 0.0977-26 ug/mL by LC-MS/MS. ClinicalTrials.gov NCT01687504. Demographics from Horita 2018 Table 1; structural model and parameters from Table 2."
  )

  ini({
    # Structural PK parameters -- Horita 2018 Table 2 final population pharmacokinetic
    # model for isoniazid. Typical values are at the cohort median weight 14.3 kg.
    lka         <- log(4.23);  label("First-order absorption rate constant ka (1/h)")                       # Table 2: ka = 4.23 1/h (RSE 6%) -- 'capped with 6.0 when selecting the base model'
    lcl_slow    <- log(4.44);  label("Apparent oral clearance CL/F for NAT2 slow acetylators at WT = 14.3 kg (L/h)")     # Table 2: CL/F slow    = 4.44 L/h (RSE 5%)
    lcl_nonslow <- log(8.08);  label("Apparent oral clearance CL/F for NAT2 nonslow acetylators at WT = 14.3 kg (L/h)")  # Table 2: CL/F nonslow = 8.08 L/h (RSE 6%)
    lvc         <- log(16.6);  label("Apparent central volume V1/F at WT = 14.3 kg (L)")                    # Table 2: V1/F = 16.6 L (RSE 4%)
    lq          <- log(8.46);  label("Apparent inter-compartmental clearance Q/F at WT = 14.3 kg (L/h)")    # Table 2: Q/F  = 8.46 L/h (RSE 21%)
    lvp         <- log(1.07);  label("Apparent peripheral volume V2/F at WT = 14.3 kg (L)")                 # Table 2: V2/F = 1.07 L (RSE 44%)

    # Allometric exponents on body weight -- fixed at canonical theoretical values
    # for CL/F and Q/F (0.75) and V1/F and V2/F (1.0). Horita 2018 Results 'INH'
    # paragraph 1.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (fixed, unitless)")   # Horita 2018 Results 'INH' paragraph 1: fixed exponent
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on V1/F and V2/F (fixed, unitless)")  # Horita 2018 Results 'INH' paragraph 1: fixed exponent

    # Inter-individual variability. Table 2 IIV column reports 'omega (CV%)' on the
    # log scale (variance = omega^2). The source uses SEPARATE omegas for the two
    # NAT2 groups -- Monolix's 'different parameters per category' option -- so the
    # CL/F IIV is encoded as two independent etas selected via NAT2_SLOW in model().
    etalka         ~ 0.321   # Table 2: 0.567 (61.6% CV)  -- 0.567^2 = 0.321 (ka)
    etalcl_slow    ~ 0.105   # Table 2: 0.324 (33.3% CV)  -- 0.324^2 = 0.105 (CL/F slow acetylators)
    etalcl_nonslow ~ 0.230   # Table 2: 0.480 (50.9% CV)  -- 0.480^2 = 0.230 (CL/F nonslow acetylators)
    etalvc         ~ 0.0581  # Table 2: 0.241 (24.5% CV)  -- 0.241^2 = 0.0581 (V1/F)
    etalq          ~ 0.406   # Table 2: 0.637 (70.7% CV)  -- 0.637^2 = 0.406 (Q/F)
    etalvp         ~ 3.61    # Table 2: 1.90 (599.7% CV)  -- 1.90^2  = 3.61   (V2/F; very large IIV reflects the small typical V2/F estimate)

    # Combined residual error. Table 2: 'Constant a' = 0.0393 (RSE 19%), 'Slope b'
    # = 0.193 (RSE 6%). See the vignette Errata for the Monolix combined-1 vs
    # nlmixr2 combined-2 distinction.
    addSd  <- 0.0393; label("Additive residual SD (ug/mL)")                  # Table 2: constant a = 0.0393
    propSd <- 0.193;  label("Proportional residual SD (fraction)")           # Table 2: slope b    = 0.193
  })

  model({
    # NAT2-genotype-dependent typical-value clearance and IIV selection. The
    # NAT2_SLOW indicator picks the slow-acetylator typical-value lcl_slow with its
    # corresponding etalcl_slow, OR the nonslow lcl_nonslow with etalcl_nonslow.
    # Each subject only "uses" one of the two etas; the other is a draw that does
    # not enter the prediction.
    lcl_i  <- lcl_slow * NAT2_SLOW + lcl_nonslow * (1 - NAT2_SLOW)
    eta_cl <- etalcl_slow * NAT2_SLOW + etalcl_nonslow * (1 - NAT2_SLOW)

    # Individual PK parameters with allometric weight scaling (reference 14.3 kg).
    ka <- exp(lka + etalka)
    cl <- exp(lcl_i + eta_cl) * (WT / 14.3)^e_wt_cl
    vc <- exp(lvc + etalvc)   * (WT / 14.3)^e_wt_vc
    q  <- exp(lq  + etalq)    * (WT / 14.3)^e_wt_cl
    vp <- exp(lvp + etalvp)   * (WT / 14.3)^e_wt_vc

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1

    # Concentration: dose mg / V1 L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
