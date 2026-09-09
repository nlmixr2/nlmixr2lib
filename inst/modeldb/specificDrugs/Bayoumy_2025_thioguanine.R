Bayoumy_2025_thioguanine <- function() {
  description <- "One-compartment population PK model with first-order absorption for erythrocyte 6-thioguanine nucleotide (6-TGN) concentrations after oral thioguanine in adults with inflammatory bowel disease, with fixed allometric body-weight scaling and a concomitant-aminosalicylate (5-ASA) effect on apparent clearance (Bayoumy 2025)"
  reference <- "Bayoumy AB, de Boer NKH, Keizer RJ, Derijks LJJ. Population pharmacokinetics model of thioguanine in patients with inflammatory bowel disease. Clin Pharmacokinet. 2025;64:1255-1262. doi:10.1007/s40262-025-01532-1"
  vignette <- "Bayoumy_2025_thioguanine"
  # Time base is DAYS: the final NONMEM control stream (ESM 1) estimates
  # Ka = 8 /day and CL/F = 110.405 L/day, which the paper's Table 2 reports
  # per hour (0.33 /h and 4.6 L/h). Concentrations are the historical 6-TGN
  # assay unit, picomole per 8e8 red blood cells, which the authors treat as
  # picomole per 200 microlitre of blood (Methods 2.3).
  units <- list(time = "day", dosing = "mg", concentration = "pmol/8e8 RBC")

  compartmentData <- list(
    # The depot holds the administered prodrug thioguanine; the central
    # compartment holds its active nucleotide metabolites (6-TGN), which are
    # what the HPLC-UV assay measures inside erythrocytes (Methods 2.1). The
    # apparent one-compartment disposition therefore lumps prodrug
    # absorption, intracellular anabolism to 6-TGN, and 6-TGN elimination.
    # Amounts are nanomoles: f(depot) converts the mg dose to nmol using the
    # thioguanine molecular weight (ESM 1 $PK, F1 = 1E6/167.19).
    depot   = list(analyte = "thioguanine", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "6-thioguanine nucleotides (6-TGN)", units = "nmol", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed allometric scaling on both disposition parameters, reference 70 kg: exponent 0.75 on CL/F and 1.0 on V/F (Bayoumy 2025 Results 3.3 and the final model equations; ESM 1 $THETA 5 and 6, both FIX). Median 60 kg (IQR 56-74.25) in the analysis population (Table 1). Height was not recorded, so fat-free mass and other body-size descriptors could not be evaluated.",
      source_name        = "WT"
    ),
    CONMED_AMINO = list(
      description        = "Concomitant aminosalicylate (5-ASA) use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant aminosalicylate)",
      notes              = "Power-form multiplier on CL/F: CL/F * 1.58165^CONMED_AMINO, i.e. 58% higher apparent clearance on concomitant 5-ASA (Bayoumy 2025 Results 3.3; ESM 1 $THETA 7 and $PK MU_3 = ... + LOG(THETA(7))*USEASA). Source column name USEASA. Binary only: the paper could not estimate a relationship between 5-ASA dose level and CL/F given the available data. 23/28 patients (82%) were on concomitant 5-ASA (Table 1).",
      source_name        = "USEASA"
    )
  )

  covariatesDataExcluded <- list(
    SNP_TPMT = list(
      description = "Thiopurine S-methyltransferase (TPMT) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Genotyped in 17/28 patients (61%), all of whom were TPMT *1/*1 wild type, so a TPMT effect could not be estimated (Bayoumy 2025 Table 1 and Limitations 4.2). The authors name TPMT genotype as the leading candidate covariate for a future refinement of this model."
    ),
    CONMED_AMINO_DOSE = list(
      description = "Daily dose of concomitant aminosalicylate (5-ASA)",
      units       = "mg",
      type        = "continuous",
      notes       = "Recorded (median 3000 mg, range 1500-3200 mg; Bayoumy 2025 Table 1) but not modelled: 'It was not possible to estimate a relationship (linear or otherwise) between 5-ASA dose level and CL due to the limited amount of data' (Results 3.3). The retained covariate is the binary CONMED_AMINO indicator."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 1L,
    n_observations = 131L,
    age_range      = "19-70 years (mean 38)",
    weight_median  = "60 kg (IQR 56-74.25)",
    sex_female_pct = 78.6,
    disease_state  = "Inflammatory bowel disease: Crohn's disease (n = 16) or ulcerative colitis (n = 12). Thioguanine was started for azathioprine intolerance (n = 9), mercaptopurine intolerance (n = 8), both (n = 8), disease activity (n = 1) or other reasons (n = 2).",
    dose_range     = "Oral thioguanine once daily; mean 0.32 mg/kg/day (95% CI 0.29-0.34), i.e. roughly 20 mg/day at the 60 kg population median weight.",
    regions        = "The Netherlands",
    co_medication  = "Concomitant aminosalicylate (5-ASA) in 23/28 patients (82%), median 5-ASA dose 3000 mg/day (range 1500-3200).",
    notes          = "Re-analysis of the cohort of Derijks et al. (Aliment Pharmacol Ther 2004), reference [9] of the source paper. Two data sets were pooled: (i) 28 patients sampled for trough 6-TGN at days 7, 14, 21, 28 and 56 of daily dosing, and (ii) three patients with a densely sampled 24-h steady-state profile (duplicate samples at 0, 2, 4, 8, 12 and 24 h). Because a joint fit destabilised the model and between-occasion variability could not be estimated from three subjects, the 24-h profiles were entered as separate subjects (Limitations 4.2). 131 measurements were available and 7 (5.0%) were removed as outliers -- 3 for |CWRES| > 3 and 4 for suspected non-adherence. 6-TGN was measured in erythrocytes by reversed-phase HPLC with UV detection (run-to-run CV 6.6%, LLOQ 30 pmol per 8e8 RBC). Estimation used SAEM followed by importance sampling in NONMEM 7.4. This is the first published population PK model of thioguanine in IBD; it was developed to support model-informed precision dosing."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters, in the control stream's DAY time base
    # (ESM 1 "Supplement 1. Final model code"). Reference covariates:
    # WT = 70 kg, CONMED_AMINO = 0.
    #
    # Table 2 of the paper reports the same estimates per hour, and its Unit
    # column is transposed between the two rows: it prints "CL/F ... L" and
    # "V/F ... L/h". The control stream is unambiguous -- CL/F is a clearance
    # (L/day) and V/F is a volume (L) -- and is used here.
    # ---------------------------------------------------------------------
    lka <- fixed(log(8)); label("Absorption rate constant Ka (1/day)")                 # ESM 1 $THETA 1: (0.1, 8) FIX "; 1 Ka /day"; 8/day = 0.33/h, Table 2 "0.33 (f)"
    lcl <- log(110.405);  label("Apparent oral 6-TGN clearance CL/F (L/day)")                 # ESM 1 $THETA 3: (0, 110.405) "; 3 CL/F L/day"; 110.405/24 = 4.60 L/h, Table 2 "4.6"
    # V/F carries the authors' x5 unit-conversion factor. ESM 1 $PK sets
    #   V_FACT = 5
    #   MU_2   = LOG(THETA(2) * V_FACT) + LOG(WT/70)*THETA(5)
    #   V      = EXP(MU_2 + ETA(2)); S2 = V
    # so the volume that actually scales the central amount -- and that sets
    # the elimination rate kel = CL/V -- is 270.46 * 5 = 1352.3 L, not the
    # 270.5 L printed in Table 2. The factor 5 converts nmol/L to the assay's
    # pmol per 200 microlitre (1 nmol/L = 0.2 pmol/200 uL); see Methods 2.3.
    lvc <- log(270.46 * 5); label("Apparent 6-TGN volume of distribution V/F (L), including the x5 assay unit-conversion factor")  # ESM 1 $THETA 2: (0, 270.46) x V_FACT 5 ($PK); Table 2 prints the unscaled 270.5

    # Allometric exponents, both fixed at the canonical values (Results 3.3:
    # "Estimating both exponents was not better than fixing to 0.75 and 1.0,
    # respectively. Therefore, the exponents were fixed in the final model.")
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless)")  # ESM 1 $THETA 6: (0, 0.75) FIX "; 6 WT on CL"
    e_wt_vc <- fixed(1.0);  label("Allometric exponent of WT on V/F (unitless)")   # ESM 1 $THETA 5: (0, 1.0) FIX "; 5 WT on V"

    # Concomitant aminosalicylate effect on CL/F, power form on the binary
    # indicator: CL/F * 1.58165^CONMED_AMINO (58% higher CL/F on 5-ASA).
    e_conmed_amino_cl <- 1.58165; label("Concomitant aminosalicylate multiplier on CL/F (power form: CL/F * 1.58165^CONMED_AMINO)")  # ESM 1 $THETA 7: (0, 1.58165) "; 7 USE of 5-ASA on CL"; Results 3.3 "58% higher CL", RSE 15%

    # ---------------------------------------------------------------------
    # Between-subject variability (omega^2 on the log scale).
    #
    # ESM 1 assigns ETA(1) to Ka, ETA(2) to V and ETA(3) to CL, and declares
    #   $OMEGA  1  FIX      ; BSV Ka
    #   $OMEGA  BLOCK(2)
    #    0.271878           ; BSV V
    #    0.119927 0.150194  ; BSV CL
    # so omega_V  = sqrt(0.271878) = 0.521 (52% CV) and
    #    omega_CL = sqrt(0.150194) = 0.388 (39% CV),
    # with correlation 0.119927 / sqrt(0.271878 * 0.150194) = 0.59.
    #
    # NOTE -- the paper's Table 2 BSV column assigns 52% to CL/F and 39% to
    # V/F, i.e. the OPPOSITE of the control stream, and Sections 3.4 and 4
    # repeat that assignment. The control stream is followed here because it
    # is the final model code, its ETA ordering and its own inline comments
    # agree with each other, and the adjacent Unit column of the same two
    # Table 2 rows is demonstrably transposed as well (see lcl / lvc above).
    # The variance magnitudes themselves are not in dispute -- only which
    # parameter each belongs to. See the vignette's Errata section.
    # ---------------------------------------------------------------------
    etalvc + etalcl ~ c(0.271878,
                        0.119927, 0.150194)  # ESM 1 $OMEGA BLOCK(2): omega^2 V = 0.271878, cov = 0.119927, omega^2 CL = 0.150194
    etalka ~ fixed(1)                        # ESM 1 $OMEGA 1 FIX "; BSV Ka"; Results 3.2: BSV on Ka was unidentifiable (>250%) and held at 100%
    etapropSd ~ 0.244129                     # ESM 1 $OMEGA 0.244129 "; eps" -> ETA(4), the BSV on the residual-error magnitude; sqrt = 0.494 (49% CV), Table 2 "BSV 49%"

    # ---------------------------------------------------------------------
    # Residual error. ESM 1 $ERROR uses log-transform-both-sides with
    #   PROP = THETA(4)*EXP(ETA(4));  W = PROP;  Y = LOG(F) + W*EPS(1)
    # and $SIGMA 1 FIX, i.e. a per-subject proportional error whose magnitude
    # is 0.145424 * exp(etapropSd). Results 3.2 states this log-transformed
    # error model is "equivalent to a proportional error model".
    # ---------------------------------------------------------------------
    propSd <- 0.145424; label("Proportional residual error scale (fraction; typical 14.5% before the BSV on residual-error magnitude)")  # ESM 1 $THETA 4: (0, 0.145424) "; 4 prop error"; Table 2 "Residual error magnitude 14.5%"
  })

  model({
    # Individual parameters. Reference covariates WT = 70 kg,
    # CONMED_AMINO = 0. Fixed allometric exponents applied inline; the
    # 5-ASA effect is a power of the binary indicator, exactly as
    # ESM 1 $PK writes LOG(THETA(7))*USEASA inside the log-scale MU.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * e_conmed_amino_cl^CONMED_AMINO
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Unit conversion carried on the depot, replicating ESM 1 $PK:
    #   F1 = 1E6/167.19   ; convert mg dose to nanomole, Mw 6-TG = 167.19 g/mol
    # Doses are therefore given in mg and amounts are carried in nmol.
    f(depot) <- 1e6 / 167.19

    # With amounts in nmol and vc in L (already including the x5 factor),
    # central / vc is directly in the assay's pmol per 8e8 RBC
    # (= pmol per 200 uL of blood; Methods 2.3).
    Cc <- central / vc

    # Per-subject proportional residual error: the paper estimates BSV on the
    # residual-error magnitude itself (ETA(4) in ESM 1 $ERROR).
    propSd_i <- propSd * exp(etapropSd)
    Cc ~ prop(propSd_i)
  })
}
