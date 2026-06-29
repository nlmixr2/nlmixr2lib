Rietveld_2025_docetaxel <- function() {
  description <- paste(
    "Integrated plasma + tumour population PK model for docetaxel (DTX)",
    "delivered as CPC634, a core-crosslinked polymeric micelle that",
    "covalently entraps DTX via a pH-responsive sulfone-ester linker",
    "(Rietveld 2025). Released / conventional DTX is described by a",
    "canonical three-compartment IV plasma model (central, peripheral1,",
    "peripheral2; CL, Vc, Q1, Vp1, Q2, Vp2). Unreleased DTX (still bound",
    "to CPC634) is described by a two-compartment plasma model on the",
    "paper-specific compartments entrapped + peripheral_entrapped with",
    "linear elimination CLcpc and intercompartmental clearance Qcpc.",
    "Release of DTX from CPC634 in plasma is time-dependent: six",
    "first-order release rates K122 / K123 / K124 / K125 / K126 / K12",
    "active in the post-dose time windows 0-0.5 / 0.5-1 / 1-2 / 2-6 /",
    "6-168 / 168+ hours (paper Methods 2.6 + Figure 1A). The two",
    "tumour-tissue states tumor_entrapped (unreleased DTX in tumour)",
    "and tumor_released (released DTX in tumour) are connected to plasma",
    "by an influx parameter Kbtn (unreleased DTX, one-way plasma ->",
    "tumour) and an in/efflux balance parameter KbtDTX (released DTX),",
    "and to each other by the tumour-local release rate KrelT; both",
    "share the same tumour distribution volume VcT (paper Methods 2.7).",
    "Inter-individual variability on CL (released-DTX clearance), Vcpc",
    "(CPC634 central volume), and K122 (first-window release rate).",
    "Additive-on-log residual error (equivalent to proportional in",
    "linear space) estimated separately for the four therapeutically",
    "relevant observation streams (released-DTX plasma, unreleased-DTX",
    "plasma, released-DTX tumour, total-DTX tumour). Pragmatic",
    "deviations from the published model: (1) the $MIXTURE on Qcpc",
    "(subpopulation 1, P = 0.69, Qcpc1 = 0.00122 L/h; subpopulation 2,",
    "Qcpc2 = 0.00769 L/h) is encoded as a typical value at the dominant",
    "subpopulation 1 (Qcpc = Qcpc1); the minority Qcpc2 = 0.00769 L/h",
    "is documented in the vignette Assumptions and deviations section",
    "and can be plugged in by overriding lqcpc <- log(0.00769) at",
    "simulation time. (2) The 89Zr-CPC634 radiotracer arm of the",
    "PICCOLO PET imaging study (compartments 6 and 7 in the supplement;",
    "fast-loss rate K002 = 0.336 1/h active in the first 2 h after the",
    "89Zr dose) is omitted -- the radiotracer arm shared the unreleased",
    "DTX disposition parameters (Vcpc / VPcpc / Qcpc / CLcpc) so it",
    "added no structural information beyond what entrapped already",
    "carries. See the vignette Assumptions and deviations section."
  )
  reference <- paste(
    "Rietveld PCS, Koolen SLW, Zeiser S, Rijcken CJF, van Noort M,",
    "van Eerden RAG, Atrafi F, Miedema IHC, Menke-van der Houven van",
    "Oordt CW, Koch BCP, Mathijssen RHJ, Snelder N, Sassen SDT.",
    "Drug release from docetaxel-entrapped core-crosslinked polymeric",
    "micelles: A population pharmacokinetic modelling approach based on",
    "clinical data. Biomed Pharmacother. 2025;185:118028.",
    "doi:10.1016/j.biopha.2025.118028."
  )
  vignette <- "Rietveld_2025_docetaxel"

  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c(
    "entrapped", "peripheral_entrapped",
    "tumor_entrapped", "tumor_released"
  )

  paper_specific_etas <- c("etalk_release_1")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 52L,
    n_studies       = 3L,
    n_observations  = 1811L,
    n_cycles        = "72 cycles CPC634 + 24 cycles conventional DTX + 15 cycles 89Zr-CPC634",
    age_range       = NA_character_,
    weight_range    = NA_character_,
    sex_female_pct  = NA_real_,
    disease_state   = paste(
      "Adults with advanced solid tumours pooled across three phase I /",
      "phase II / imaging trials of CPC634 (core-crosslinked polymeric",
      "micelle entrapping docetaxel via a pH-responsive sulfone-ester",
      "linker): NAPOLY phase I dose-escalation (n = 23; NCT02442531,",
      "15-100 mg/m^2 CPC634 IV every 2-3 weeks), CRITAX cross-over",
      "(n = 24; NL6299, 75 mg/m^2 CPC634 vs conventional DTX), and",
      "PICCOLO PET-imaging (n = 5; NCT03712423, 0.1-2 mg 89Zr-CPC634",
      "diagnostic dose + 60 mg/m^2 on-treatment CPC634). Tumour biopsies",
      "from CRITAX (n = 24 patients, two biopsies each at 24, 48, 72,",
      "96, 168, or 336 h after dose)."
    ),
    dose_range      = paste(
      "CPC634 15-100 mg/m^2 IV infusion (NAPOLY); 75 mg/m^2 CPC634 and",
      "75 mg/m^2 conventional DTX cross-over (CRITAX); 0.1-2 mg",
      "89Zr-CPC634 diagnostic dose followed two weeks later by 60 mg/m^2",
      "CPC634 + 89Zr-CPC634 on-treatment (PICCOLO). All routes IV",
      "infusion."
    ),
    regions         = "Netherlands (Erasmus MC, Amsterdam UMC VUmc).",
    notes           = paste(
      "Baseline demographics (age, weight, BSA, sex) are not tabulated",
      "in the source paper or supplement; the authors note in the",
      "Discussion that an exhaustive covariate screen (dose, body",
      "surface area, age) found no covariate explaining the Qcpc",
      "subpopulation split, but they do not publish the cohort-level",
      "demographics. NONMEM 7.5 with FOCE+I, Pirana 2.9.9, R 4.3.2;",
      "NONMEM2R 0.2.5 and Parid 1.0 for postprocessing. Nonparametric",
      "bootstrap n = 1000 for parameter precision (Methods 2.9). Six",
      "first-order release-rate constants and a $MIXTURE block on Qcpc",
      "drive the unreleased-DTX disposition; see paper Methods 2.6 +",
      "Table 2 and the supplement's $PK + $DES + $THETA blocks."
    )
  )

  ini({
    # ================================================================
    # Released / conventional DTX -- three-compartment IV plasma model
    # (paper Table 2 'Released/Conventional DTX plasma model' rows).
    # All values are typical final estimates.
    # ================================================================
    lcl  <- log(26.9)  ; label("Released DTX clearance, CL (L/h)")                                    # Table 2 row CL = 26.9, RSE 9%
    lvc  <- log(7.18)  ; label("Released DTX central volume, Vc (L)")                                 # Table 2 row Vc = 7.18, RSE 13%
    lq   <- log(16.5)  ; label("Released DTX intercompartmental clearance to peripheral1, Q1 (L/h)")  # Table 2 row Q1 = 16.5, RSE 11%
    lvp  <- log(1350)  ; label("Released DTX peripheral1 volume, Vp1 (L)")                            # Table 2 row Vp1 = 1350, RSE 12%
    lq2  <- log(9.34)  ; label("Released DTX intercompartmental clearance to peripheral2, Q2 (L/h)")  # Table 2 row Q2 = 9.34, RSE 11%
    lvp2 <- log(17.0)  ; label("Released DTX peripheral2 volume, Vp2 (L)")                            # Table 2 row Vp2 = 17, RSE 13%

    # ================================================================
    # CPC634 nanoparticle (unreleased DTX) -- two-compartment IV plasma
    # model. Paper Table 2 'Unreleased DTX plasma model' rows.
    # Qcpc value carried is the dominant (subpopulation 1) mixture
    # estimate (P = 0.69; Qcpc1 = 0.00122 L/h); the minority Qcpc2 is
    # preserved as a fixed anchor lqcpc_alt for documentation /
    # alternative-scenario simulation. See description and vignette
    # Assumptions and deviations.
    # ================================================================
    lvcpc  <- log(3.56)    ; label("CPC634 central volume, Vcpc (L)")                                 # Table 2 row Vcpc = 3.56, RSE 4%
    lvpcpc <- log(0.408)   ; label("CPC634 peripheral volume, VPcpc (L)")                             # Table 2 row VPcpc = 0.408, RSE 17%
    # Qcpc carried at the dominant (subpop 1) mixture estimate. The minority subpop 2 value
    # Qcpc2 = 0.00769 L/h (Table 2, P = 0.31, RSE 51%) is documented in description / vignette;
    # users who want subpop-2 behaviour can override lqcpc with log(0.00769) at simulation time.
    lqcpc  <- log(0.00122) ; label("CPC634 intercompartmental clearance Qcpc -- subpop 1 (L/h)")      # Table 2 row Qcpc1 = 0.00122, RSE 38% (P = 0.69)
    lclcpc <- log(0.0156)  ; label("CPC634 clearance, CLcpc (L/h)")                                   # Table 2 row CLcpc = 0.0156, RSE 19%

    # ================================================================
    # Time-dependent first-order release of DTX from CPC634 in plasma.
    # Six release-rate constants active in successive post-dose windows
    # (paper Methods 2.6 + supplement $DES; supplement encodes Krel via
    # IF(T.LT.0.5) ... IF(T.GE.168) ... cascade on time after dose).
    # IIV is estimated only on the first-window rate K122.
    # ================================================================
    lk_release_1 <- log(0.162)
    label("DTX release rate K122 (0-0.5 h post dose, 1/h)")                                            # Table 2 row K122 = 0.162, RSE 12%
    lk_release_2 <- log(0.0928)
    label("DTX release rate K123 (0.5-1 h post dose, 1/h)")                                            # Table 2 row K123 = 0.0928, RSE 7%
    lk_release_3 <- log(0.0699)
    label("DTX release rate K124 (1-2 h post dose, 1/h)")                                              # Table 2 row K124 = 0.0699, RSE 5%
    lk_release_4 <- log(0.0448)
    label("DTX release rate K125 (2-6 h post dose, 1/h)")                                              # Table 2 row K125 = 0.0448, RSE 5%
    lk_release_5 <- log(0.0199)
    label("DTX release rate K126 (6-168 h post dose, 1/h)")                                            # Table 2 row K126 = 0.0199, RSE 5%
    lk_release_terminal <- log(0.0157)
    label("DTX release rate K12 (>168 h post dose, 1/h)")                                              # Table 2 row K12 = 0.0157, RSE 18%

    # ================================================================
    # Tumour distribution and release. Paper Table 2 'Tumour model'.
    # All four free parameters are estimated; the joint plasma + tumour
    # fit uses the plasma parameters from above as informative priors
    # (supplement $PRIOR NWPRI). VcT shared between unreleased and
    # released DTX tumour states.
    # ================================================================
    lkbtn       <- log(0.00017) ; label("Plasma -> tumour influx for unreleased DTX, Kbtn (1/h; one-way)") # Table 2 row Kbtn = 0.00017, RSE 26%
    lkbt_dtx    <- log(0.0112)  ; label("Tumour <-> plasma exchange for released DTX, KbtDTX (1/h, balance)") # Table 2 row KbtDTX = 0.0112, RSE 18%
    lkrel_tumor <- log(0.00123) ; label("Tumour-local DTX release rate, KrelT (1/h)")                  # Table 2 row KrelT = 0.00123, RSE 41%
    lvct        <- log(220)     ; label("Tumour distribution volume for unreleased + released DTX, VcT (L)") # Table 2 row VcT = 220, RSE 14%

    # ================================================================
    # Inter-individual variability. NONMEM $OMEGA values from the
    # supplementary control stream's $OMEGA BLOCK(3) diagonal -- the
    # block is diagonal here (off-diagonals = 0), so the three IIVs are
    # independent. Table 2 reports the same parameters as approximate
    # CV%; the discrepancy with sqrt($OMEGA) is the standard
    # lognormal-CV-vs-log-scale approximation. omega^2 values are used
    # directly. Reported CVs: CL = 30.6%, Vcpc = 26.9%, K122 = 63.6%.
    # ================================================================
    etalcl          ~ 0.0966   # supplement $OMEGA BLOCK(3) row 1: IIV-CL = 0.0966
    etalvcpc        ~ 0.0729   # supplement $OMEGA BLOCK(3) row 2: IIV-Vcpc = 0.0729
    etalk_release_1 ~ 0.394    # supplement $OMEGA BLOCK(3) row 3: IIV-K122 = 0.394

    # ================================================================
    # Residual error. Supplement $ERROR carries IPRED = log(C) and
    # Y = IPRED + EPS, so the additive-on-log-scale error is
    # proportional in nlmixr2's linear space. SD = sqrt($SIGMA value).
    # Table 2 reports the $SIGMA variance entries as 'Add ERR CMT n';
    # propSd_<output> = sqrt(variance).
    # ================================================================
    propSd                <- sqrt(0.251)
    label("Released DTX plasma proportional residual SD (fraction)")                                   # supplement $SIGMA row 1 (CMT 2) = 0.251, RSE 5%; Table 2 Add ERR CMT 2
    propSd_Cc_entrapped   <- sqrt(0.247)
    label("Unreleased DTX (CPC634) plasma proportional residual SD (fraction)")                        # supplement $SIGMA row 2 (CMT 1) = 0.247, RSE 6%; Table 2 Add ERR CMT 1
    propSd_Cc_tumor_total <- sqrt(0.54)
    label("Total DTX (released + unreleased) tumour proportional residual SD (fraction)")              # supplement $SIGMA row 4 (CMT 7) = 0.54, RSE 35%; Table 2 Add ERR CMT 7
    propSd_Cc_tumor       <- sqrt(0.432)
    label("Released / conventional DTX tumour proportional residual SD (fraction)")                    # supplement $SIGMA row 5 (CMT 5) = 0.432, RSE 31%; Table 2 Add ERR CMT 5
  })

  model({
    # ------------------------------------------------------------
    # Individual structural parameters.
    # ------------------------------------------------------------
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    q2     <- exp(lq2)
    vp2    <- exp(lvp2)

    vcpc   <- exp(lvcpc + etalvcpc)
    vpcpc  <- exp(lvpcpc)
    qcpc   <- exp(lqcpc)
    clcpc  <- exp(lclcpc)

    kbtn       <- exp(lkbtn)
    kbt_dtx    <- exp(lkbt_dtx)
    krel_tumor <- exp(lkrel_tumor)
    vct        <- exp(lvct)

    # ------------------------------------------------------------
    # Micro-constants.
    # ------------------------------------------------------------
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    k13  <- q2 / vc
    k31  <- q2 / vp2

    kel_cpc <- clcpc / vcpc
    k_cpc_p <- qcpc / vcpc
    k_p_cpc <- qcpc / vpcpc

    # ------------------------------------------------------------
    # Time-dependent first-order release of DTX from CPC634 in plasma.
    # tadose = time elapsed since the most recent dose into the
    # entrapped (CPC634) compartment. The if/else cascade reproduces
    # the supplement $DES block.
    # ------------------------------------------------------------
    tadose <- t - tlast(entrapped)
    krel <- exp(lk_release_1 + etalk_release_1)
    if (tadose >= 0.5) {
      krel <- exp(lk_release_2)
    }
    if (tadose >= 1.0) {
      krel <- exp(lk_release_3)
    }
    if (tadose >= 2.0) {
      krel <- exp(lk_release_4)
    }
    if (tadose >= 6.0) {
      krel <- exp(lk_release_5)
    }
    if (tadose >= 168.0) {
      krel <- exp(lk_release_terminal)
    }

    # ------------------------------------------------------------
    # ODE system. Paper Figure 1 + supplement $DES.
    #   entrapped <-> peripheral_entrapped: CPC634 nanoparticle
    #     two-compartment distribution with linear elimination (CLcpc)
    #     and time-dependent first-order release of DTX into the
    #     released-DTX central compartment.
    #   central / peripheral1 / peripheral2: released-DTX three-cmt
    #     disposition (canonical 3-compartment).
    #   tumor_entrapped: unreleased DTX in tumour. Receives one-way
    #     plasma -> tumour influx (Kbtn) from the CPC634 plasma central
    #     compartment and feeds the tumour-released state via the
    #     intra-tumour release rate KrelT.
    #   tumor_released: released DTX in tumour. In/out balance with
    #     the released-DTX plasma central via KbtDTX (the source paper
    #     uses a single balance parameter for both directions); fed by
    #     the intra-tumour release rate KrelT from tumor_entrapped.
    # ------------------------------------------------------------
    d/dt(entrapped) <-
      -krel * entrapped -
      k_cpc_p * entrapped + k_p_cpc * peripheral_entrapped -
      kel_cpc * entrapped

    d/dt(peripheral_entrapped) <-
      k_cpc_p * entrapped - k_p_cpc * peripheral_entrapped

    d/dt(central) <-
      krel * entrapped -
      kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2

    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    d/dt(tumor_entrapped) <-
      kbtn * entrapped - krel_tumor * tumor_entrapped

    d/dt(tumor_released) <-
      kbt_dtx * central - kbt_dtx * tumor_released +
      krel_tumor * tumor_entrapped

    # ------------------------------------------------------------
    # Observation variables. Doses in mg, volumes in L give central /
    # vc in mg/L which equals ug/mL. Multiplying by 1000 yields ng/mL.
    # Tumour concentrations are reported in ng/mg by the source paper
    # (assuming soft-tissue density approx 1 g/mL: ng/mg = ng/uL =
    # ng/mL); we report tumour Cc in the same ng-equivalent scale
    # using vct in L, so the multiplier is also 1000.
    # ------------------------------------------------------------
    Cc              <- 1000 * central / vc                          # released / conventional DTX plasma (ng/mL)
    Cc_entrapped    <- 1000 * entrapped / vcpc                      # unreleased DTX (CPC634) plasma (ng/mL DTX-equivalent)
    Cc_tumor        <- 1000 * tumor_released / vct                  # released / conventional DTX tumour (ng/mg DTX-equivalent)
    Cc_tumor_total  <- 1000 * (tumor_entrapped + tumor_released) / vct # total DTX tumour (ng/mg DTX-equivalent)

    Cc              ~ prop(propSd)
    Cc_entrapped    ~ prop(propSd_Cc_entrapped)
    Cc_tumor        ~ prop(propSd_Cc_tumor)
    Cc_tumor_total  ~ prop(propSd_Cc_tumor_total)
  })
}
