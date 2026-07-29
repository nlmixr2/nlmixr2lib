deWit_2016_everolimus <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption",
    "for everolimus 10 mg once-daily in 40 adult patients with advanced",
    "thyroid carcinoma (de Wit 2016). Bioavailability F is structurally",
    "fixed at 1 (absolute F unknown), so reported CL, V1, Q, and V2 are",
    "apparent (oral / F). Allometric scaling on apparent clearance",
    "(exponent 0.75) and apparent central volume (exponent 1.0) using a",
    "70 kg reference weight per the Anderson and Holford theory cited by",
    "the paper. Apparent peripheral volume V2/F was held fixed at 400 L",
    "in the final model. Bioavailability is multiplied by 0.792 in",
    "subjects who carry at least one ABCB1 TTT haplotype (CYP3A / P-gp",
    "efflux marker). Inter-occasion variability on F captures the",
    "day-1-vs-day-15 sampling occasion contrast (CV 19.2%)."
  )
  reference <- paste(
    "de Wit D, Schneider TC, Moes DJAR, Roozen CFM, den Hartigh J,",
    "Gelderblom H, Guchelaar HJ, van der Hoeven JJ, Links TP,",
    "Kapiteijn E, van Erp NP. Everolimus pharmacokinetics and its",
    "exposure-toxicity relationship in patients with thyroid cancer.",
    "Cancer Chemother Pharmacol. 2016;78(1):63-71.",
    "doi:10.1007/s00280-016-3050-6.",
    "ClinicalTrials.gov NCT01118065.",
    sep = " "
  )
  vignette <- "deWit_2016_everolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed body weight used for allometric scaling of apparent",
        "clearance (exponent 0.75) and apparent central volume",
        "(exponent 1.0) at a 70 kg reference per Anderson and Holford",
        "(2008), Annu Rev Pharmacol Toxicol 48:303-332 -- the reference",
        "cited as [12] in the source paper for the allometric form. The",
        "paper does not state the reference weight explicitly; 70 kg is",
        "the canonical Anderson and Holford adult reference. Median",
        "weight in the cohort was 75 kg (range 45-105 kg) per Table 1."
      ),
      source_name        = "Weight (kg)"
    ),
    ABCB1_HAP_TTT = list(
      description        = paste(
        "ABCB1 TTT haplotype carrier indicator (1 = subject carries at",
        "least one TTT haplotype across the rs1128503 / rs2032582 /",
        "rs1045642 ABCB1 SNP block, 0 = no TTT haplotype). Time-fixed",
        "per subject (germline haplotype)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no TTT haplotype)",
      notes              = paste(
        "Carriers (1+ TTT) and non-carriers were the only contrast",
        "retained after backward elimination (dOFV = 9.6, P < 0.01,",
        "Results paragraph 2). Heterozygote and homozygote TTT carriers",
        "were pooled because the homozygote frequency was < 0.1 (Methods",
        "section 'Covariate analysis'; haplotype certainty > 0.97 in",
        "gPLINK, Methods section 'Pharmacogenetic analysis'). The TTT",
        "haplotype is associated with enhanced P-gp efflux activity",
        "(Discussion paragraph 5), producing a 21% reduction in apparent",
        "bioavailability via the multiplier theta_TTT = 0.792 on F",
        "(Table 2 final model)."
      ),
      source_name        = "ABCB1 TTT haplotype carrier (Methods 'Pharmacogenetic analysis', Table 2 final model)"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued occasion indicator (1 = day 1 sampling, 2 = day",
        "15 sampling). Multi-day steady-state contrast."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Day-1 and day-15 are the only two PK sampling occasions in the",
        "study (Methods, 'Pharmacokinetic sample collection and",
        "analysis'). Decomposed inside model() into binary indicators",
        "oc1 and oc2 multiplexing per-occasion etas on the log-",
        "bioavailability lfdepot. The IOV variance is single-valued",
        "(CV 19.2% in Table 2 final-model 'Inter-occasion variability'",
        "row) so the occasion-2 eta variance is fix()'d to match",
        "occasion-1 (the NONMEM $OMEGA BLOCK(1) SAME pattern; nlmixr2",
        "has no SAME shortcut)."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "40-80 years",
    age_median     = "63 years",
    weight_range   = "45-105 kg",
    weight_median  = "75 kg",
    sex_female_pct = 47.5,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Adults with advanced thyroid carcinoma enrolled in a phase II",
      "trial (NCT01118065): 26/40 (65%) differentiated, 7/40 (17.5%)",
      "undifferentiated (anaplastic), 7/40 (17.5%) medullary. Two of the",
      "originally enrolled 42 patients were excluded from the PK",
      "analysis (no samples collected in one; no measurable everolimus",
      "levels in the other)."
    ),
    dose_range     = paste(
      "Everolimus 10 mg orally once daily continuous dosing until tumor",
      "progression / unacceptable toxicity / death / withdrawal. First",
      "dose reduction (when needed) was to 5 mg once daily; second",
      "reduction to 5 mg every other day."
    ),
    regions        = "The Netherlands (Leiden University Medical Center + University Medical Center Groningen)",
    notes          = paste(
      "Demographics from de Wit 2016 Table 1 (median, range). PK",
      "sampling occurred on study days 1 and 15: predose and 1, 2, 3 h",
      "post-dose for everyone (sparse schedule), with an optional",
      "extension to 4, 5, 6, 7, 8 h (extensive schedule). 30 of 40",
      "patients completed the extensive schedule; 669 whole-blood",
      "everolimus concentrations contributed to the population model.",
      "NONMEM 7.2 with FOCE-INTER; bootstrap n=1000 (Table 2 'Bootstrap",
      "runs' column). Other genotypes / haplotypes tested but not",
      "retained (Methods 'Pharmacogenetic analysis', Results paragraph",
      "2): CYP2C8 haplotype, single SNPs in CYP3A4, CYP3A5, ABCB1,",
      "NR1I2 (PXR)."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- de Wit 2016 Table 2 ('Final model')
    # column. Apparent values reported with F structurally fixed at 1
    # (Methods 'Base model', paragraph 3). Allometric scaling on CL/F
    # and V1/F at a 70 kg reference per Anderson and Holford [ref 12].
    # ============================================================
    lcl <- log(17.4)
    label("Apparent oral clearance CL/F at WT = 70 kg (L/h)")            # Table 2 Final: Cl/F = 17.4 L/h (RSE 8.4%)
    lvc <- log(25.2)
    label("Apparent central volume V1/F at WT = 70 kg (L)")              # Table 2 Final: V1/F = 25.2 L (RSE 17.8%)
    lka <- log(0.647)
    label("First-order oral absorption rate constant ka (1/h)")          # Table 2 Final: ka = 0.647 1/h (RSE 6.2%)
    lq  <- log(51.1)
    label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table 2 Final: Q = 51.1 L/h (RSE 7.3%)
    lvp <- fixed(log(400))
    label("Apparent peripheral volume V2/F (L; FIXED)")                  # Table 2 Final: V2 = 400 L (no RSE in 'Final model' column, no bootstrap CI in 'Bootstrap runs' column -- held fixed in the final model)

    # Bioavailability anchor (F = 1 structurally) and ABCB1 TTT
    # multiplier (Methods 'Base model' paragraph 3 and Table 2 Final
    # 'theta TTT on F').
    lfdepot <- fixed(log(1))
    label("Reference bioavailability F (unitless, FIXED at 1)")          # Methods 'Base model' paragraph 3: F was fixed at 1 because absolute bioavailability is unknown
    e_abcb1_hap_ttt_fdepot <- 0.792
    label("Multiplicative effect of ABCB1 TTT carrier status on F (unitless)") # Table 2 Final: theta TTT on F = 0.792 (RSE 6.5%) -- 20.8% lower F in TTT carriers

    # Allometric exponents fixed to the canonical Anderson and Holford
    # values (paper cites [12] for the allometric form but does not
    # restate the exponents); reference weight 70 kg.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F (unitless, FIXED)")               # Methods 'Base model': allometric scaling per Anderson and Holford [12]; canonical 0.75 on CL
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on V1/F (unitless, FIXED)")               # Methods 'Base model': allometric scaling per Anderson and Holford [12]; canonical 1.0 on V

    # ============================================================
    # Inter-individual variability -- de Wit 2016 Table 2 ('Final
    # model') 'Inter-individual variability' rows. CV% reported in the
    # paper is converted to log-normal omega^2 = log(1 + CV^2).
    # ============================================================
    etalcl ~ log(1 + 0.351^2)
    # Table 2 Final IIV: CL/F = 35.1% CV (RSE 30.5%) -> omega^2 = log(1 + 0.351^2)
    etalvc ~ log(1 + 0.864^2)
    # Table 2 Final IIV: V1/F = 86.4% CV (RSE 35.3%) -> omega^2 = log(1 + 0.864^2)

    # ============================================================
    # Inter-occasion variability on F -- de Wit 2016 Table 2 ('Final
    # model') 'Inter-occasion variability F' row, CV 19.2%. Two
    # occasions (day 1, day 15) share a single estimated variance; the
    # second occasion is fix()'d to the first (NONMEM $OMEGA BLOCK(1)
    # SAME pattern, matching the Jonsson_2011_ethambutol.R and
    # Xie_2019_agomelatine.R idiom for IOV in nlmixr2).
    # ============================================================
    etaiov_fdepot_1 ~ log(1 + 0.192^2)
    # Table 2 Final IOV: F = 19.2% CV (RSE 38.1%) -> omega^2 = log(1 + 0.192^2)
    etaiov_fdepot_2 ~ fixed(log(1 + 0.192^2))
    # Occasion-2 variance fixed equal to occasion-1 (single 'Inter-occasion variability F' row in Table 2; NONMEM $OMEGA BLOCK(1) SAME)

    # ============================================================
    # Residual error -- de Wit 2016 Table 2 Final 'Residual variability
    # sigma (proportional error)' row, CV 27.3%.
    # ============================================================
    propSd <- 0.273
    label("Proportional residual error (fraction)")                      # Table 2 Final: sigma_prop = 27.3% (RSE 20.8%)
  })

  model({
    # ------------------------------------------------------------
    # Occasion indicators for IOV on bioavailability. OCC is the
    # integer-valued canonical (oc1 / oc2 binary decomposition; see
    # inst/references/covariate-columns.md OCC entry).
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # ------------------------------------------------------------
    # Individual PK parameters with allometric scaling on a 70 kg
    # reference (Methods 'Base model', paragraph 3; Anderson and
    # Holford theory-based exponents 0.75 / 1.0).
    # ------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------
    # Bioavailability. F = F_typical * theta_TTT^ABCB1_HAP_TTT *
    # exp(iov_fdepot), with theta_TTT = 0.792 (Table 2 Final 'theta
    # TTT on F') applied multiplicatively to carriers (ABCB1_HAP_TTT =
    # 1) and unchanged for non-carriers (ABCB1_HAP_TTT = 0). The IOV
    # term swaps in the occasion-1 vs occasion-2 eta via oc1/oc2.
    # ------------------------------------------------------------
    fdepot <- exp(lfdepot) *
              e_abcb1_hap_ttt_fdepot^ABCB1_HAP_TTT *
              exp(iov_fdepot)

    # ------------------------------------------------------------
    # Two-compartment oral PK (depot -> central <-> peripheral1) with
    # first-order absorption and first-order elimination (Results
    # paragraph 1: 'two-compartmental model with first-order
    # absorption and first-order elimination from the central
    # compartment').
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ------------------------------------------------------------
    # Observation. Dose units mg; vc units L -> compartment amount in
    # mg, central / vc in mg/L. Convert mg/L to ug/L (paper's unit) by
    # multiplying by 1000.
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
