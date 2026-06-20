Chigutsa_2011_rifampicin <- function() {
  description <- "Population pharmacokinetic model for oral rifampicin in adults with sputum-positive pulmonary tuberculosis in South Africa (Cape Town). One-compartment disposition with a fixed-length Erlang transit-absorption chain (NN = 19 fixed) feeding the central compartment via first-order ka. Allometric scaling of CL/F and V/F to a 70 kg reference body weight with canonical Anderson and Holford (2008) exponents (0.75 on CL, 1.0 on V; cited as Chigutsa 2011 Methods reference 3 for the allometric model). Covariate effects: female sex on V/F (-30%) and on the mean transit time MTT (+30% per Results body text page 4124 -- women have a 30% LONGER absorption delay than men; Table 2 Final-model row prints -30% with a CI bit-identical to the V/F row immediately above, which is the canonical signature of a typesetting row-duplication error; per the operator sidecar request-001 directive the body text +30% is the source of truth); high-dose-band effect on MTT (-27% for daily doses >= 600 mg vs the 450 mg reference); SLCO1B1 rs4149032 genotype-dependent oral bioavailability F (heterozygous carriers -18%; homozygous variant carriers -28%; relative to the homozygous-common-allele wild-type reference). Between-subject variability (BSV) is carried on F, CL, and MTT with the CL-MTT correlation block 0.86 from Table 2; within-subject (WSV / IOV) variability reported in Table 2 is NOT carried (forward-simulation users do not need the second-occasion IOV layer; see vignette Errata). Combined additive + proportional residual error."
  reference <- "Chigutsa E, Visser ME, Swart EC, Denti P, Pushpakom S, Egan D, Holford NHG, Smith PJ, Maartens G, Owen A, McIlleron H. (2011). The SLCO1B1 rs4149032 polymorphism is highly prevalent in South Africans and is associated with reduced rifampin concentrations: dosing implications. Antimicrob Agents Chemother 55(9):4122-4127. doi:10.1128/AAC.01833-10"
  vignette <- "Chigutsa_2011_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the Chigutsa 2011 cohort (weight-banded dosing only; no within-treatment body weight remodelling). Used for allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1.0) with reference weight 70 kg per Anderson and Holford (Annu Rev Pharmacol Toxicol 2008;48:303-32; cited as Chigutsa 2011 Methods reference 3).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Effects on apparent central volume V/F (-30%; Chigutsa 2011 Table 2 row 'Effect of female sex on V/F (%) -30 (-21, -35)') and on mean transit time MTT (+30%, encoded per Results body text page 4124 stating 'women had a 30% LONGER mean transit time, showing that women have a longer absorption delay than men.' Table 2 row 'Effect of female sex on MTT (%)' prints -30 (-21, -35), but the CI bounds are bit-identical to the V/F row immediately above which is the canonical signature of a typesetting row-duplication error; per sidecar request-001 (2026-06-17) the body text +30% is the source of truth. See vignette Errata for the full discussion).",
      source_name        = "SEXF"
    ),
    DOSE = list(
      description        = "Per-record administered rifampicin daily dose (mg).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to derive a binary high-dose indicator (DOSE >= 600 mg) for the dose effect on mean transit time MTT (-27% relative to the 450 mg reference). The Chigutsa 2011 cohort received weight-banded doses: 38-54 kg -> 450 mg (3 Rifafour tablets), 55-70 kg -> 600 mg (4 tablets), >70 kg -> 750 mg (5 tablets). Per-record because rxode2 routes the simulated dose amount via the standard amt column.",
      source_name        = "DOSE"
    ),
    SNP_SLCO1B1_RS4149032_COUNT = list(
      description        = "Count of SLCO1B1 rs4149032 variant alleles per subject (0 = homozygous common allele wild-type, 1 = heterozygous one variant allele, 2 = homozygous two variant alleles). rs4149032 is an intronic SLCO1B1 SNP; Chigutsa 2011 Discussion notes it is in linkage disequilibrium with the functional SLCO1B1 variant C463A (rs11045819) in Caucasians but NOT in Africans, where the haplotype patterns differ due to greater genetic diversity.",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (germline genotype). Distribution in the Chigutsa 2011 cohort (Results paragraph 'The SLCO1B1 rs4149032 polymorphism...'): 31 of 60 homozygous variant (52%), 22 heterozygous (37%), 7 homozygous common (12%); overall variant allele frequency 0.70 (0.93 in the black African subcohort, 0.59 in the mixed-race subcohort). Used in model() as a decomposition: het_4149032 <- (SNP_SLCO1B1_RS4149032_COUNT == 1) and hom_4149032 <- (SNP_SLCO1B1_RS4149032_COUNT == 2). Independently estimated non-additive multiplicative bioavailability effects: -18% in HET, -28% in HOM variant; relative to WT homozygous reference. Follows the Schipani 2011 nevirapine count-decomposition precedent for non-additive het / hom shifts.",
      source_name        = "rs4149032 genotype"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 1L,
    age_range      = "18 - 55 years (median 29)",
    age_median     = "29 years",
    weight_range   = "41 - 69 kg (median 52)",
    weight_median  = "52 kg",
    sex_female_pct = 40,
    race_ethnicity = "40% black African (mainly Xhosa); 60% mixed ancestry (mainly Caucasian and African origins)",
    disease_state  = "Adults with sputum smear-positive pulmonary tuberculosis (Cape Town, South Africa). 16% co-infected with HIV.",
    dose_range     = "Oral rifampicin 450 mg (38-54 kg), 600 mg (55-70 kg), or 750 mg (>70 kg) once daily, administered as fixed-dose-combination Rifafour tablets (each containing 150 mg rifampicin, 75 mg isoniazid, 400 mg pyrazinamide, 275 mg ethambutol) under directly observed administration after at least 1 month of treatment.",
    regions        = "South Africa (Cape Town; Delft Community Health Centre)",
    sampling_design = "Sparse steady-state sampling: 4-8 plasma samples per patient collected randomly over a 7-h window after the morning dose; 24 of 57 patients had a second occasion sampled approximately 1 month after the first. Total 437 plasma concentrations. BLQ data (12% of observations, all in the absorption phase < 2 h post-dose) handled by Beal's M3 method in the original NONMEM fit.",
    assay          = "Plasma rifampicin by HPLC-MS/MS; LLOQ 0.08 mg/L.",
    notes          = "Cohort participants were enrolled in a parallel randomised micronutrient (vitamin A and zinc vs placebo) substudy; the micronutrient arm assignment had no significant effect on rifampicin PK and is not encoded as a covariate. The full Chigutsa 2011 Table 2 reports BSV on F, CL, and MTT (CL-MTT block correlation 0.86) and WSV on F, CL, V, and MTT (V-MTT block correlation -0.40); this model file carries the BSV layer but omits the WSV / inter-occasion-variability layer (forward simulation users do not need IOV for a single-occasion simulation; see vignette Errata)."
  )

  ini({
    # ============================================================
    # Structural PK at 70 kg reference body weight (Chigutsa 2011
    # Table 2 'Final model' block).
    # ============================================================
    lcl <- log(11)
    label("Apparent oral clearance CL/F at 70 kg (L/h)")
    # Chigutsa 2011 Table 2 row 'CL/F (liters/h/70 kg) = 11 (10, 13)'.

    lvc <- log(50)
    label("Apparent central volume of distribution V/F at 70 kg (L)")
    # Chigutsa 2011 Table 2 row 'V/F (liters/70 kg) = 50 (41, 53)'.

    lka <- log(1.1)
    label("First-order absorption rate constant ka (1/h)")
    # Chigutsa 2011 Table 2 row 'ka (1/h) = 1.1 (0.9, 1.4)'.

    lmtt <- log(1.6)
    label("Mean transit time MTT (h)")
    # Chigutsa 2011 Table 2 row 'MTT (h) = 1.6 (1.3, 1.8)'.

    lnn <- fixed(log(19))
    label("Erlang transit-chain shape NN (unitless; fixed)")
    # Chigutsa 2011 Table 2 row 'NN = 19 (fixed)'. The Savic transit
    # chain shape parameter is fixed at 19 transit compartments;
    # rxode2's transit(n, mtt, bio) closed-form gamma-density input
    # accepts non-integer n but Chigutsa 2011 reports an integer 19.

    lfdepot <- fixed(log(1))
    label("Bioavailability F (typical-value anchor; fixed)")
    # Chigutsa 2011 (Table 2 'BSV of F' and Results body text) treats F
    # as the bioavailability anchor relative to which the SLCO1B1
    # rs4149032 covariate effect is reported. Population value fixed
    # at 1 because absolute F is not identifiable from oral-only data;
    # the CL/F and V/F estimates are reported as the apparent
    # F-relative values per the Table 2 row labels.

    # ============================================================
    # Allometric exponents fixed at canonical Anderson and Holford
    # (2008) values per Chigutsa 2011 Methods reference 3.
    # ============================================================
    allo_cl <- fixed(0.75)
    label("Allometric exponent on CL (unitless; fixed)")
    # Chigutsa 2011 Methods reference 3 (Anderson and Holford 2008);
    # standard allometric scaling for clearance.

    allo_v <- fixed(1.0)
    label("Allometric exponent on V (unitless; fixed)")
    # Chigutsa 2011 Methods reference 3 (Anderson and Holford 2008);
    # standard allometric scaling for volume.

    # ============================================================
    # Covariate effects (fractional, additive in the 1 + e_*
    # multiplicative form). All values from Chigutsa 2011 Table 2.
    # ============================================================
    e_sexf_vc <- -0.30
    label("Fractional effect of female sex on apparent central volume V/F")
    # Chigutsa 2011 Table 2 'Effect of female sex on V/F (%) = -30
    # (-21, -35)' and Results page 4124 body text: 'Women had a V/F
    # that was 30% lower than that of men.' Sign and magnitude
    # internally consistent between table and text.

    e_sexf_mtt <- 0.30
    label("Fractional effect of female sex on mean transit time MTT")
    # Chigutsa 2011 Results page 4124 body text (source of truth per
    # operator sidecar request-001 directive 2026-06-17): 'women had
    # a 30% LONGER mean transit time, showing that women have a longer
    # absorption delay than men.' Table 2 'Effect of female sex on
    # MTT (%) = -30 (-21, -35)' is bit-identical to the V/F row
    # immediately above (typesetting row-duplication artifact); the
    # body-text sign +30% is encoded. See vignette Errata.

    e_dose_high_mtt <- -0.27
    label("Fractional effect of daily dose >= 600 mg on mean transit time MTT")
    # Chigutsa 2011 Table 2 'Effect of dose on MTT (%) = -27 (22, 36)'
    # and Results page 4124 body text: 'both men and women given
    # higher doses (600 or 750 mg of rifampin daily) were found to
    # have a 27% shorter absorption delay than those given 450 mg
    # daily.' Encoded with a binary indicator (DOSE >= 600). Table 2
    # CI bounds (22, 36) are missing minus signs in the PDF render
    # (rendering artifact -- not a sign error); the body-text 'shorter'
    # confirms negative direction.

    e_snp_slco1b1_rs4149032_het_fdepot <- -0.18
    label("Fractional effect of SLCO1B1 rs4149032 heterozygous (1 variant allele) on bioavailability F")
    # Chigutsa 2011 Table 2 'Effect of SLCO1B1 rs41490932 on F in
    # heterozygotes (%) = -18 (-15, -25)'. (Note: the paper writes
    # 'rs41490932' for the polymorphism in Table 2 and elsewhere in
    # Results; an extra '9' typo for rs4149032 -- the canonical dbSNP
    # rsid in the Methods text and Title. Body text page 4124:
    # 'Heterozygotes had 18% lower bioavailability ... than
    # homozygotes for the SLCO1B1 rs4149032 allele' [text inverts the
    # comparison wording but the Table 2 reference is unambiguously
    # the homozygous-common-allele wild-type subjects].)

    e_snp_slco1b1_rs4149032_hom_fdepot <- -0.28
    label("Fractional effect of SLCO1B1 rs4149032 homozygous variant (2 variant alleles) on bioavailability F")
    # Chigutsa 2011 Table 2 'Effect of SLCO1B1 rs41490932 on F in
    # variant homozygotes (%) = -28 (-19, -34)'. Body text page 4124:
    # 'homozygotes had 28% lower bioavailability than homozygotes for
    # the SLCO1B1 rs4149032 allele' [reference is unambiguously the
    # homozygous-common-allele wild-type subjects per Table 2 column
    # header 'variant homozygotes'].

    # ============================================================
    # Between-subject variability (BSV; NONMEM omega^2 scale).
    # Chigutsa 2011 Table 2 reports the BSV omega^2 directly (not as
    # CV%): BSV of F = 0.15; BSV of CL = 0.20; BSV of MTT = 0.52;
    # Correlation BSV(CL)-BSV(MTT) = 0.86.
    # Block off-diagonal: cov(CL, MTT) = corr * sqrt(var_cl * var_mtt)
    # = 0.86 * sqrt(0.20 * 0.52) = 0.27734.
    # WSV (within-subject / inter-occasion) rows from Table 2 are NOT
    # carried; see vignette Errata for the rationale.
    # ============================================================
    etalcl + etalmtt ~ c(0.20, 0.27734, 0.52)
    # Chigutsa 2011 Table 2 'BSV of CL = 0.20 (0.17, 0.32)',
    # 'BSV of MTT = 0.52 (0.45, 0.84)',
    # 'Correlation between BSV of CL and MTT = 0.86 (0.81, 0.96)';
    # off-diagonal = 0.86 * sqrt(0.20 * 0.52) = 0.27734.

    etalfdepot ~ 0.15
    # Chigutsa 2011 Table 2 'BSV of F = 0.15 (0.10, 0.16)'.

    # ============================================================
    # Residual error (Chigutsa 2011 Table 2 'Additive error' and
    # 'Proportional error' rows). Combined additive + proportional.
    # ============================================================
    addSd <- 0.03
    label("Additive residual error SD on rifampicin Cc (mg/L)")
    # Chigutsa 2011 Table 2 'Additive error (mg/liter) = 0.03
    # (0.02, 0.04)'.

    propSd <- 0.30
    label("Proportional residual error SD on rifampicin Cc (fraction)")
    # Chigutsa 2011 Table 2 'Proportional error = 0.30 (0.24, 0.33)'.
  })

  model({
    # ------------------------------------------------------------
    # 1. Body-weight allometric factors (Anderson and Holford 2008,
    #    Chigutsa 2011 Methods reference 3); reference 70 kg.
    # ------------------------------------------------------------
    bw_cl <- (WT / 70) ^ allo_cl
    bw_v  <- (WT / 70) ^ allo_v

    # ------------------------------------------------------------
    # 2. SLCO1B1 rs4149032 genotype indicators decomposed from the
    #    variant-allele count (non-additive HET / HOM effects on F).
    # ------------------------------------------------------------
    het_4149032 <- (SNP_SLCO1B1_RS4149032_COUNT == 1)
    hom_4149032 <- (SNP_SLCO1B1_RS4149032_COUNT == 2)
    snp_fdepot_factor <- 1 +
                         e_snp_slco1b1_rs4149032_het_fdepot * het_4149032 +
                         e_snp_slco1b1_rs4149032_hom_fdepot * hom_4149032

    # ------------------------------------------------------------
    # 3. High-dose indicator (>= 600 mg/day) for the MTT covariate.
    # ------------------------------------------------------------
    dose_high <- (DOSE >= 600)

    # ------------------------------------------------------------
    # 4. Individual PK parameters with allometric weight scaling,
    #    female-sex effects on V/F and MTT, and high-dose effect on
    #    MTT. SNP effect on F is multiplicative on the fixed F = 1
    #    anchor; BSV on F enters in log space via etalfdepot.
    # ------------------------------------------------------------
    cl  <- exp(lcl  + etalcl)  * bw_cl
    vc  <- exp(lvc)             * bw_v * (1 + e_sexf_vc * SEXF)
    ka  <- exp(lka)
    mtt <- exp(lmtt + etalmtt) *
           (1 + e_sexf_mtt * SEXF) *
           (1 + e_dose_high_mtt * dose_high)
    nn  <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot) * snp_fdepot_factor

    # ------------------------------------------------------------
    # 5. Micro-constant.
    # ------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------
    # 6. ODE system. Erlang transit-chain absorption (closed form)
    #    feeds depot via the gamma-density input from the most recent
    #    dose; ka drives depot -> central; first-order elimination
    #    from central.
    # ------------------------------------------------------------
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # transit() is the only depot input; disable the normal F-based
    # dose accumulation so the gamma-density flux is the unique
    # source term.
    f(depot) <- 0

    # ------------------------------------------------------------
    # 7. Observation: plasma rifampicin concentration (mg/L) with
    #    combined additive + proportional residual error.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
