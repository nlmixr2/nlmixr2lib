Vinnard_2017_isoniazid <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for oral isoniazid in HIV/tuberculosis co-infected adults in Botswana (Vinnard 2017). NAT2 acetylator genotype (slow acetylator as the reference, with proportional shifts for the intermediate and rapid phenotypes) and systemic immune activation (percent CD38 and HLA-DR co-expression on CD8+ T cells, entered as a median-normalised power term) act on apparent oral clearance; between-subject variability on CL/F, V/F and the absorption lag time, and inter-occasion variability on CL/F across the pre-ART and post-ART pharmacokinetic visits."
  reference <- paste(
    "Vinnard C, Ravimohan S, Tamuhla N, Ivaturi V, Pasipanodya J,",
    "Srivastava S, Modongo C, Zetola NM, Weissman D, Gumbo T, Bisson GP.",
    "(2017). Isoniazid clearance is impaired among human immunodeficiency",
    "virus/tuberculosis patients with high levels of immune activation.",
    "Br J Clin Pharmacol 83(4):801-811. doi:10.1111/bcp.13172.",
    sep = " "
  )
  vignette <- "Vinnard_2017_isoniazid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Isoniazid was quantified in SERUM by stable-isotope
  # dilution LC-ESI-MS/MS (Methods 'Serum isoniazid concentrations'), so the
  # sampled matrix is serum rather than plasma.
  compartmentData <- list(
    depot       = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "isoniazid", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    NAT2_SLOW = list(
      description        = "NAT2 slow-acetylator phenotype indicator (1 = slow acetylator, 0 = intermediate or rapid). Paired with NAT2_RAPID to encode the three-level NAT2 phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (slow acetylator) is the REFERENCE level of the source paper's CL/F covariate model; the typical value lcl is the slow-acetylator clearance.",
      notes              = "Vinnard 2017 Methods 'NAT-2 genotype': whole-exome sequencing; individuals carrying any combination of NAT2 *4, *11, *12, *13 were classified rapid acetylators, those carrying two of *5, *6, *7, *14 slow acetylators, and those carrying one allele from each group intermediate acetylators. Cohort distribution at visit 1 (Table 1, journal page 805): 13 rapid (32.5%), 18 intermediate (45%), 7 slow (17.5%), 2 ambiguous (5%); 38 of 40 participants had a non-ambiguous genotype. Unlike Horita 2018 (which pooled intermediate with rapid), this paper estimated a separate proportional CL/F shift for each of the intermediate and rapid phenotypes against the slow reference, so NAT2_SLOW and NAT2_RAPID are BOTH required and the joint state (0, 0) denotes an intermediate acetylator. The two ambiguous-genotype participants were not assigned a phenotype by the paper; see population$notes.",
      source_name        = "NAT2"
    ),
    NAT2_RAPID = list(
      description        = "NAT2 rapid (fast) acetylator phenotype indicator (1 = rapid acetylator, 0 = intermediate or slow). Paired with NAT2_SLOW to encode the three-level NAT2 phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0; the joint state NAT2_SLOW = 0 and NAT2_RAPID = 0 denotes an intermediate acetylator, and NAT2_SLOW = 1 is the reference level of the source paper's CL/F covariate model.",
      notes              = "Same phenotyping source as NAT2_SLOW (Vinnard 2017 Methods 'NAT-2 genotype'). Carries the source paper's Theta(Rapid NAT2) = 1.65 proportional increase in CL/F relative to the slow-acetylator reference (Table 2, journal page 808).",
      source_name        = "NAT2"
    ),
    CD8_CD38DR_PCT = list(
      description        = "Percentage of circulating CD8+ T cells co-expressing CD38 and HLA-DR (%CD38+DR+CD8+), the systemic immune-activation marker retained on apparent oral clearance in the final model.",
      units              = "% of CD8+ T cells",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Vinnard 2017 Methods 'Immune activation markers': cryopreserved PBMCs stained for CD3, CD8, CD4, CD38 and HLA-DR and acquired on an LSR II flow cytometer. Enters CL/F as a power term normalised to 36.9%, the visit-1 cohort median (Table 1, journal page 805; IQR 27.7-45.7, n = 38); the visit-2 median was 24.8% (IQR 21.9-35.9, n = 23) and the individual values plotted in Figure 2B (journal page 808) span roughly 11% to 64%. The paper prints no covariate equation, so the median-normalised power form was recovered from the published anchors; see the vignette source-trace and Assumptions sections for the four independent checks that identify it.",
      source_name        = "%CD38+DR+CD8+"
    ),
    OCC = list(
      description        = "Integer-valued pharmacokinetic study-visit occasion used to multiplex the inter-occasion-variability etas on apparent oral clearance.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Two occasions in the Vinnard 2017 cohort (Methods 'Data collection'): OCC = 1 is the first PK visit, 5 to 28 days after starting anti-TB therapy and before ART initiation; OCC = 2 is the second PK visit, a median of 33 days after ART initiation. Decomposed inside model() into the binary indicators oc1 and oc2 that select the two IOV etas on log CL/F. Records with OCC outside {1, 2} receive zero IOV contribution.",
      source_name        = "OCC"
    )
  )

  # Screened by the paper but NOT retained in the final model, so they are
  # documented here rather than in covariateData (they are never referenced in
  # model()). Results paragraph 3: 'Covariate relationships that were examined,
  # but not retained in the model, included creatinine clearance, sex, and
  # weight as covariates on V, CL, the volume of the peripheral compartment,
  # and intercompartmental clearance.' Results paragraph 4: 'No additional
  # covariate relationships with other immunologic parameters were identified
  # in the nonlinear mixed effects model (including CD4+ T cell count, IL-6,
  # neopterin, and CRP)'. No point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL/F, V/F, V2/F and Q/F and not retained. Median 102.1 mL/min (IQR 92.5-114.1) at visit 1 and 102.8 (97.8-123.0) at visit 2 (Table 1, journal page 805). Enrolment excluded creatinine clearance below 50 mL/min. The Discussion notes that changes in creatinine clearance between visits were not measured, which the authors list as a limitation."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F, V/F, V2/F and Q/F and not retained. 18 of 40 (45%) women at visit 1 and 12 of 24 (50%) at visit 2 (Table 1, journal page 805)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL/F, V/F, V2/F and Q/F and not retained; no allometric scaling appears in the final model. Median 55.0 kg (IQR 49.3-59.3) at visit 1 and 56.6 kg (52.5-61.8) at visit 2 (Table 1, journal page 805). Weight nonetheless set each participant's WHO weight-band isoniazid dose."
    ),
    CD4_ABS = list(
      description = "Absolute CD4+ T-lymphocyte count",
      units       = "cells/uL",
      type        = "continuous",
      notes       = "Screened on CL/F as an immunologic covariate and not retained. Median 238 cells/uL (IQR 105-339, n = 40) at visit 1, rising to 308 (212-400, n = 22) at visit 2 (Table 1, journal page 805). The Discussion draws the contrast explicitly: 'CD4+ T cell counts did not exert a significant covariate effect on isoniazid clearance, suggesting that immune activation rather than immune suppression drives the variability in NAT-2 activity.'"
    ),
    IL6 = list(
      description = "Plasma interleukin-6 concentration",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "Screened on CL/F as a secondary systemic-inflammation marker and not retained. Median 14.3 pg/mL (IQR 6.6-26.7, n = 39) at visit 1 and 6.7 (3.1-14.2, n = 24) at visit 2 (Table 1, journal page 805)."
    ),
    CRP = list(
      description = "Plasma C-reactive protein concentration",
      units       = "ug/mL",
      type        = "continuous",
      notes       = "Screened on CL/F as a secondary systemic-inflammation marker and not retained. Median 10.0 ug/mL (IQR 4.0-19.2, n = 36) at visit 1 and 5.6 (2.1-16.5, n = 23) at visit 2 (Table 1, journal page 805)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "21 years and older (enrolment criterion)",
    age_median     = "32 years (IQR 27-43) at visit 1; 32 years (IQR 28-43) at visit 2",
    weight_range   = "Not tabulated as a range; median 55.0 kg (IQR 49.3-59.3) at visit 1 and 56.6 kg (IQR 52.5-61.8) at visit 2",
    weight_median  = "55.0 kg",
    sex_female_pct = 45,
    race_ethnicity = "Citizens of Botswana (sub-Saharan African); detailed ancestry not reported.",
    disease_state  = "ART-naive HIV-infected adults newly diagnosed with pulmonary TB and established on a standard WHO first-line antitubercular regimen under directly observed therapy.",
    dose_range     = "Oral isoniazid once daily as part of a first-line fixed-dose combination, dosed by WHO weight band (Methods 'Study population'). The paper does not tabulate the administered isoniazid doses; the slow-acetylator individual CL/F and AUC0-inf pairs in Figure 3 (journal page 809) imply 300 mg for that subgroup (see the vignette Assumptions section).",
    regions        = "Botswana (Gaborone; 22 public clinics and Princess Marina Hospital).",
    renal_function = "Creatinine clearance below 50 mL/min was an exclusion criterion; median CrCl 102.1 mL/min (IQR 92.5-114.1) at visit 1.",
    hepatic_function = "Alanine or aspartate transaminase above 3 times the upper limit of normal was an exclusion criterion.",
    co_medication  = "First-line antitubercular fixed-dose combination therapy at both visits; tenofovir/emtricitabine/efavirenz ART in all participants by visit 2.",
    notes          = "Prospective observational two-visit design. 61 patients were screened and 40 enrolled and sampled at visit 1 (median 20 days, range 7-65, after starting anti-TB therapy, ART-naive); 24 returned for visit 2 (median 74 days, range 33-118, of anti-TB therapy and a median 33 days, range 5-44, of ART). Serum sampled at 0, 0.3, 0.9, 2.2, 4.5 and 8 h post-dose by optimal-design sampling; LLOQ 0.16 mg/L. Baseline demographics are Table 1 (journal page 805) and the final parameter estimates Table 2 (journal page 808). NAT2 genotype was non-ambiguous in 38 of 40 participants (7 slow, 18 intermediate, 13 rapid); the paper does not state how the 2 ambiguous-genotype participants were handled in the covariate model. Eta shrinkage in the final model was 10.0% for CL/F, 12.8% for V/F and 43.0% for the absorption lag time (Results paragraph 4). Estimation used FOCE in Phoenix NLME 1.3."
  )

  ini({
    # Structural PK parameters. All values are the point estimates in the
    # 'Typical value (% RSE)' column of Vinnard 2017 Table 2 (journal page
    # 808), 'Final model parameter estimates for adult population
    # pharmacokinetics of isoniazid in human immunodeficiency
    # virus/tuberculosis patients'. Dosing was oral throughout and no
    # bioavailability term was estimated, so every clearance and volume is an
    # APPARENT (per-F) quantity.
    lka   <- log(0.88);  label("First-order absorption rate constant ka (1/h)")                           # Table 2: Ka = 0.88 1/h (RSE 12.45%)
    lcl   <- log(10.99); label("Apparent oral clearance CL/F, NAT2 slow acetylator at 36.9% CD38+DR+CD8+ (L/h)") # Table 2: CL = 10.99 L/h (RSE 9.86%); the reference level of the CL/F covariate model
    lvc   <- log(5.55);  label("Apparent central volume of distribution V/F (L)")                         # Table 2: V = 5.55 L (RSE 48.11%)
    lvp   <- log(17.54); label("Apparent peripheral volume of distribution V2/F (L)")                     # Table 2: V2 = 17.54 L (RSE 16.29%)
    lq    <- log(11.61); label("Apparent intercompartmental clearance Q/F (L/h)")                         # Table 2: Q = 11.61 L/h (RSE 34.06%)
    ltlag <- log(0.25);  label("Absorption lag time Tlag (h)")                                            # Table 2: Tlag = 0.25 h (RSE 11.67%)

    # Covariate effects on apparent oral clearance. The paper prints no
    # covariate equation; the two functional forms below were identified from
    # the published anchors and are quantitatively verified against four
    # independent checks in the vignette source-trace section. In brief:
    # Figure 3A (journal page 809) plots individual CL/F by NAT2 genotype with
    # per-genotype pre-ART medians near 10, 20 and 27 L/h, which the
    # proportional (1 + theta) form reproduces as 10.99, 17.9 and 29.1 L/h and
    # which an exponential exp(theta) form cannot (it would put the rapid
    # typical value at 57 L/h, off the panel's axis).
    e_nat2_int_cl       <-  0.63; label("Proportional increase in CL/F, NAT2 intermediate vs slow acetylator (unitless)")  # Table 2: Theta(Intermediate NAT2) = 0.63 (RSE 30.62%); 95% CI 0.30 to 1.06
    e_nat2_rapid_cl     <-  1.65; label("Proportional increase in CL/F, NAT2 rapid vs slow acetylator (unitless)")         # Table 2: Theta(Rapid NAT2) = 1.65 (RSE 18.22%); 95% CI 1.14 to 2.25
    e_cd8_cd38dr_pct_cl <- -0.31; label("Power exponent on median-normalised %CD38+DR+CD8+ for CL/F (unitless)")           # Table 2: Theta(%CD38 + HLA-DR + |CD8) = -0.31 (RSE 40.06%); 95% CI -0.57 to -0.07

    # Between-subject variability. Table 2 reports BSV on the percent scale;
    # converted to log-normal variances with the exact identity
    # omega^2 = log(1 + CV^2). Table 2: BSV-CL 29.3% CV (RSE 2.4%).
    etalcl   ~ 0.08237
    # Table 2: BSV-V 139.3% CV (RSE 86.1%). This is the paper's estimate as
    # printed; V/F itself carries an RSE of 48.11% and a bootstrap 95% CI of
    # 2.01 to 12.31 L, so the central volume is only weakly identified.
    etalvc   ~ 1.07858
    # Table 2: BSV-Tlag 9.7% CV (RSE 1.7%). Eta shrinkage on the lag time was
    # 43.0% (Results paragraph 4).
    etaltlag ~ 0.00937

    # Inter-occasion variability on apparent oral clearance across the two PK
    # visits (OCC 1 = pre-ART, OCC 2 = post-ART). Table 2: BSV-IOV 4.5% CV
    # (RSE 0.1%), a single variance shared by both occasions. nlmixr2 has no
    # NONMEM-style 'SAME' shortcut, so the second occasion's variance is held
    # equal to the first (the pattern used by Wilkins_2008_rifampicin.R and by
    # the companion Vinnard_2017_rifampicin.R).
    etaiov_cl_1 ~ 0.00202
    etaiov_cl_2 ~ fix(0.00202)

    # Combined additive plus proportional residual error (Table 2).
    # Concentrations are in mg/L (equivalently ug/mL).
    addSd  <- 0.06; label("Additive residual error (mg/L)")             # Table 2: additive error SD = 0.06 mg/L (RSE 16.68%); the assay LLOQ was 0.16 mg/L
    propSd <- 0.25; label("Proportional residual error (fraction, CV)") # Table 2: proportional error = 0.25 (RSE 12.27%)
  })

  model({
    # 1. Derived covariate terms.

    #    NAT2 phenotype. The source paper's reference level is the SLOW
    #    acetylator, so lcl is the slow-acetylator clearance and both printed
    #    thetas raise it. The canonical register supplies binary indicators for
    #    the slow and rapid phenotypes; the intermediate phenotype is the joint
    #    (0, 0) state and is derived here rather than carried as a third
    #    column. Proportional form: Figure 3A's per-genotype CL/F medians.
    nat2_int <- (1 - NAT2_SLOW) * (1 - NAT2_RAPID)
    nat2_cl  <- 1 + e_nat2_int_cl * nat2_int + e_nat2_rapid_cl * NAT2_RAPID

    #    Systemic immune activation. Power term normalised to 36.9%, the
    #    visit-1 cohort median of %CD38+DR+CD8+ (Table 1, journal page 805).
    #    The normalisation makes lcl the clearance of a slow acetylator at the
    #    cohort-median activation level, and it is what makes the printed
    #    Theta = -0.31 reproduce the Figure 2B regression slope of about
    #    -0.009 per percentage point (theta / 36.9 = -0.0084).
    immune_cl <- (CD8_CD38DR_PCT / 36.9)^e_cd8_cd38dr_pct_cl

    #    Occasion decomposition for the IOV etas on log CL/F.
    oc1    <- (OCC == 1)
    oc2    <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # 2. Individual parameters. BSV was retained on CL/F, V/F and Tlag only
    #    (Results paragraph 3); ka, V2/F and Q/F carry no IIV in the published
    #    model.
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl + iov_cl) * nat2_cl * immune_cl
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    tlag <- exp(ltlag + etaltlag)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment system with first-order absorption from the depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Absorption lag. Results paragraph 3: 'Fit of the concentrations during
    #    the absorptive phase was improved with the use of a time-lag
    #    absorption model'.
    alag(depot) <- tlag

    # 6. Serum isoniazid concentration. Dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
