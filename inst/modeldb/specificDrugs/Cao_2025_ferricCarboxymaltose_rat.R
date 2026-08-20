Cao_2025_ferricCarboxymaltose_rat <- function() {
  description <- paste(
    "QSP. Preclinical (rat, Sprague-Dawley).",
    "Mechanism-based PK/PD model of intravenous ferric carboxymaltose in",
    "iron-deficiency-anemia rats. A two-compartment serum-iron PK model with",
    "endogenous zero-order input drives a maturation-structured hematopoietic",
    "model in which the change in serum iron from baseline shifts commitment",
    "of hematopoietic stem and progenitor cells between the erythroid lineage",
    "(BFU-E -> CFU-E -> normoblast -> reticulocyte -> RBC) and the",
    "megakaryocytic lineage (10 MK-precursor aging compartments -> 10 platelet",
    "aging compartments), with a separate Emax stimulation of hemoglobin",
    "production. Three simultaneous outputs: RBC count, hemoglobin and",
    "platelet count.",
    sep = " "
  )
  reference <- paste(
    "Cao K, Fan X, Wong RSM, Yan X. Mechanism-Based Pharmacokinetic/",
    "Pharmacodynamic Modeling for Iron-Regulated Hematopoietic Stem and",
    "Progenitor Cells' Commitment toward Erythroid and Megakaryocytic",
    "Lineages. ACS Pharmacol Transl Sci. 2025;8(6):1711-1725.",
    "doi:10.1021/acsptsci.5c00097",
    sep = ""
  )
  vignette <- "Cao_2025_ferricCarboxymaltose"
  units <- list(time = "h", dosing = "mg", concentration = "ug/dL")

  # The PK model is parameterised per kilogram of body weight (Table S2 reports
  # V1/V2 in L/kg, Q and CL in L/h/kg, KIN_Iron in mg/h/kg), so amounts in
  # `central` / `peripheral1` are mg/kg and doses are entered as mg/kg.

  # Compartment names not in inst/references/compartment-names.md. `prol`
  # (HSPC pool), `hb`, `RBC` and `PLT` are canonical; the erythroid-lineage
  # states and the two aging chains are declared here pending operator
  # ratification as a canonical hematopoiesis family.
  paper_specific_compartments <- c("bfue", "cfue", "nor", "ret", "rbc")
  paper_specific_compartment_pattern <- "^(mk|plt)[0-9]+$"

  compartmentData <- list(
    central     = list(analyte = "iron", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "iron", units = "mg/kg", specimen = "tissue", verified = TRUE),
    prol        = list(analyte = "hematopoietic stem and progenitor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    bfue        = list(analyte = "burst-forming unit-erythroid cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    cfue        = list(analyte = "colony-forming unit-erythroid cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    nor         = list(analyte = "normoblasts", units = "cells/L", specimen = "tissue", verified = TRUE),
    ret         = list(analyte = "reticulocytes", units = "cells/L", specimen = "blood cell", verified = TRUE),
    rbc         = list(analyte = "mature red blood cells", units = "cells/L", specimen = "blood cell", verified = TRUE),
    hb          = list(analyte = "hemoglobin", units = "g/dL", specimen = "whole blood", verified = TRUE),
    mk1         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk2         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk3         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk4         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk5         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk6         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk7         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk8         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk9         = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    mk10        = list(analyte = "megakaryocytic precursor cells", units = "cells/L", specimen = "tissue", verified = TRUE),
    plt1        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt2        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt3        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt4        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt5        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt6        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt7        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt8        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt9        = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE),
    plt10       = list(analyte = "platelets", units = "cells/L", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species      = "rat (Sprague-Dawley)",
    n_subjects   = 24,
    n_studies    = 1,
    weight_range = "160-200 g at acquisition (Supporting Information, Animals)",
    sex_female_pct = 0,
    disease_state = paste(
      "Iron deficiency anemia induced in male Sprague-Dawley rats by a low-iron",
      "diet (<10 mg Fe/kg) for the whole experiment plus bi-weekly phlebotomy",
      "(1 mL per bleed) for the first three weeks followed by a two-week",
      "stabilization period; healthy controls received a 200 mg Fe/kg diet",
      "without phlebotomy",
      sep = " "
    ),
    dose_range = paste(
      "Ferric carboxymaltose 3, 15 or 90 mg/kg IV once weekly for 2 weeks,",
      "alone or combined with rHuEPO 450 IU/kg IV three times weekly;",
      "n = 3 per group",
      sep = " "
    ),
    regions = "Hong Kong SAR, China (Chinese University of Hong Kong)",
    notes = paste(
      "Naive pooled-data analysis in NONMEM 7.5.0: all individual data were",
      "treated as originating from a single individual, so the model carries",
      "NO between-subject variability. PK was estimated first (SAEM followed by",
      "importance-sampling EM) from sparse sampling in this study combined with",
      "intensive sampling from Funk 2022 (Eur J Pharm Biopharm 174:56-76) in an",
      "anemic rat model, then fixed while the PD model was fitted sequentially.",
      "Hematology was sampled on days 0, 2 and 4 of each week for 7 weeks after",
      "the first dose. The rHuEPO arms were used to demonstrate the",
      "iron-EPO interaction experimentally but no EPO term appears in the",
      "published model equations, so this model describes the iron effect only.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Serum-iron PK -- Cao 2025 Supporting Information Table S2
    # ---------------------------------------------------------------------
    lvc <- log(2.52)         ; label("Central volume of distribution of iron, V1 (L/kg)")            # Table S2 (RSE 9.5%)
    lvp <- log(0.687)        ; label("Peripheral volume of distribution of iron, V2 (L/kg)")         # Table S2 (RSE 31.4%)
    lq_cp <- log(0.264)      ; label("Central-to-peripheral distribution clearance, QCP (L/h/kg)")   # Table S2 (RSE 9.4%)
    lq_pc <- log(0.00132)    ; label("Peripheral-to-central distribution clearance, QPC (L/h/kg)")   # Table S2 (RSE 22.4%)
    lkin_iron <- log(0.000197); label("Endogenous zero-order iron input rate, KIN_Iron (mg/h/kg)")   # Table S2 (RSE 6.7%)
    # There is no regulated mechanism for iron excretion, so CL was fixed to a
    # very small number rather than estimated (Table S2 footnote a; Methods).
    lcl <- fixed(log(0.0001)); label("Clearance of iron, CL (L/h/kg)")                               # Table S2, fixed

    # ---------------------------------------------------------------------
    # Erythroid lineage -- Cao 2025 Table 1
    # ---------------------------------------------------------------------
    # Table 1 prints the reticulocyte row label as "T ERT"; the definition
    # column ("mean residence time for RETs") and every equation identify it
    # as T_RET.
    lt_ret <- log(39.92)     ; label("Mean residence time for reticulocytes, T_RET (h)")             # Table 1 (RSE 16.27%)
    lt_rbc <- log(161.8)     ; label("Mean residence time for red blood cells, T_RBC (h)")           # Table 1 (RSE 14.06%)
    lrbase_RBC <- log(7.141) ; label("Baseline RBC count, RBC0 (x10^12 cells/L)")                    # Table 1 (RSE 5.25%)
    lkdiff_ery <- log(166.5e-4); label("First-order rate constant for HSPC differentiation into BFU-E, KE (1/h)") # Table 1, 166.5 x 10^-4 /h (RSE 16.22%)

    # Disease-progression factor: a 0.87 (87%) reduction of the HSPC-to-BFU-E
    # differentiation rate relative to healthy conditions. Bounded on (0, 1),
    # so it is carried on the linear scale rather than log-transformed.
    dfprog <- 0.87               ; label("Disease progression factor on HSPC differentiation into BFU-E (fraction)") # Table 1 (RSE 4.94%)
    lsc50_df <- log(2062)    ; label("Change in serum iron giving half-maximal correction of the disease factor, SC50_DF (ug/dL)") # Table 1 (RSE 7%)
    # SmaxDF is the maximal correction of DF and was fixed to 1 (Methods, PD
    # Model), i.e. iron can fully abolish the disease effect but not overshoot.
    smax_df <- fixed(1)      ; label("Maximal correction of the disease progression factor, Smax_DF (fraction)")  # Methods, PD Model; fixed

    # Iron stimulation of HSPC commitment toward BFU-E, gated by an on/off
    # cutoff: only iron above CutoffIron promotes the erythroid lineage.
    lsmax_iron <- log(12.63) ; label("Maximal stimulation of HSPC differentiation into BFU-E, Smax_Iron (fraction)") # Table 1 (RSE 32.24%)
    lsc50_iron <- log(3.14)  ; label("Change in serum iron giving half-maximal stimulation of HSPC differentiation into BFU-E, SC50_Iron (ug/dL)") # Table 1 (RSE 38.93%)
    lcutoff_iron <- log(4.28); label("Change in serum iron above which HSPC differentiation into BFU-E is stimulated, Cutoff_Iron (ug/dL)") # Table 1 (RSE 16.72%)

    # ---------------------------------------------------------------------
    # Hemoglobin -- Cao 2025 Table 1
    # ---------------------------------------------------------------------
    lrbase_hb <- log(8.99)   ; label("Baseline hemoglobin concentration, HGB0 (g/dL)")               # Table 1 (RSE 4.94%)
    lsmax_hgb <- log(0.72)   ; label("Maximal stimulation of hemoglobin production, Smax_HGB (fraction)") # Table 1 (RSE 25.3%)
    lsc50_hgb <- log(16.24)  ; label("Change in serum iron giving half-maximal stimulation of hemoglobin production, SC50_HGB (ug/dL)") # Table 1 (RSE 37.73%)

    # ---------------------------------------------------------------------
    # Megakaryocytic lineage -- Cao 2025 Table 1
    # ---------------------------------------------------------------------
    lt_mp <- log(3.49)       ; label("Mean lifespan of megakaryocytic precursor cells, T_MP (h)")    # Table 1 (RSE 38.5%)
    lt_plt <- log(85.7)      ; label("Mean lifespan of platelets, T_PLT (h)")                        # Table 1 (RSE 18.79%)
    lrbase_PLT <- log(2.12)  ; label("Baseline platelet count, PLT0 (x10^12 cells/L)")               # Table 1 (RSE 18.67%)
    # Average number of platelets shed by a single MK cell; fixed at 4000 from
    # Krzyzanski 2013 (Pharm Res 30:655-669), reference 18.
    lcf <- fixed(log(4000))  ; label("Platelets produced per megakaryocyte, CF (platelets/cell)")    # Methods, PD Model; fixed from ref 18

    # ---------------------------------------------------------------------
    # Erythroid amplification exponents
    # ---------------------------------------------------------------------
    # 2^MCFU CFU-E are generated by one BFU-E and 2^MNOR normoblasts by one
    # CFU-E (Methods, PD Model; reference 19 Perez-Ruixo 2009). The paper does
    # not print MCFU or MNOR, but both are DETERMINED exactly by the
    # secondary-parameter table: Table S3 gives RET0 = 1.41, CFUE0 = 4.41e-2
    # and BFUE0 = 1.38e-3 (all x10^12 cells/L), and equations 24-25 make
    # RET0/CFUE0 = 2^MNOR and CFUE0/BFUE0 = 2^MCFU. Both ratios evaluate to
    # 32.0, so MCFU = MNOR = 5. Two independent ratios agree, so these are
    # determined rather than assumed.
    mcfu <- fixed(5)         ; label("log2 of the number of CFU-E generated per BFU-E, MCFU (unitless)") # Determined from Table S3 (see comment); fixed
    mnor <- fixed(5)         ; label("log2 of the number of normoblasts generated per CFU-E, MNOR (unitless)") # Determined from Table S3 (see comment); fixed

    # ---------------------------------------------------------------------
    # Residual error -- Cao 2025 Table 1
    # ---------------------------------------------------------------------
    # Table 1 labels these "sigma", not "sigma^2". Read as standard
    # deviations / CVs on the face of the table; the goodness-of-fit scatter in
    # Figure S9 supports that reading (see the vignette Errata section for the
    # magnitude check that rejects the variance interpretation).
    propSd_RBC <- 0.29       ; label("Proportional residual error for RBC (fraction)")               # Table 1, sigma_prop-RBC (RSE 7.13%)
    addSd_HGB <- 0.79         ; label("Additive residual error for hemoglobin (g/dL)")                # Table 1, sigma_add-HGB (RSE 5.11%)
    propSd_PLT <- 0.58       ; label("Proportional residual error for platelets (fraction)")         # Table 1, sigma_prop-PLT (RSE 6.80%)
  })

  model({
    # -------------------------------------------------------------------
    # Unit and structural constants
    # -------------------------------------------------------------------
    # Amounts are mg/kg and volumes L/kg, so amount/volume is mg/L; serum iron
    # is reported in ug/dL and 1 mg/L = 100 ug/dL.
    mgL_to_ugdL <- 100
    # Both aging chains use n = 10 compartments (Methods, PD Model: "MK_n,
    # where n = 10" and "PLT_n (where n = 10)"). The compartment count is
    # structural and therefore hard-coded rather than parameterised.
    nchain <- 10

    # -------------------------------------------------------------------
    # PK parameters and micro-constants (Cao 2025 equations 3-5)
    # -------------------------------------------------------------------
    vc <- exp(lvc)
    vp <- exp(lvp)
    q_cp <- exp(lq_cp)
    q_pc <- exp(lq_pc)
    cl <- exp(lcl)
    kin_iron <- exp(lkin_iron)

    kel <- cl / vc
    kcp <- q_cp / vc
    kpc <- q_pc / vp

    # -------------------------------------------------------------------
    # PD parameters
    # -------------------------------------------------------------------
    t_ret <- exp(lt_ret)
    t_rbc <- exp(lt_rbc)
    rbase_RBC <- exp(lrbase_RBC)
    kdiff_ery <- exp(lkdiff_ery)
    sc50_df <- exp(lsc50_df)
    smax_iron <- exp(lsmax_iron)
    sc50_iron <- exp(lsc50_iron)
    cutoff_iron <- exp(lcutoff_iron)
    rbase_hb <- exp(lrbase_hb)
    smax_hgb <- exp(lsmax_hgb)
    sc50_hgb <- exp(lsc50_hgb)
    t_mp <- exp(lt_mp)
    t_plt <- exp(lt_plt)
    rbase_PLT <- exp(lrbase_PLT)
    cf <- exp(lcf)

    # T_EP is the average time for a precursor to transition to the next cell
    # population. Equations 10-13 carry T_EP1, T_EP2 and T_EP3 separately, but
    # the Methods state "To simplify the model parameters, T_EP was assumed to
    # be equal to T_RET", so all three collapse to T_RET.
    t_ep <- t_ret

    # -------------------------------------------------------------------
    # Baseline (secondary) parameters -- Cao 2025 equations 22-29,
    # tabulated in Supporting Information Table S3
    # -------------------------------------------------------------------
    # The measured RBC count is the sum of the reticulocyte and mature-RBC
    # pools, which is what makes equation 22 divide by (T_RBC + T_RET) rather
    # than by T_RBC alone. Splitting RBC0 this way puts the erythroid chain in
    # exact flux balance at time zero: ret0 / T_RET == rbc0 / T_RBC.
    ret0 <- rbase_RBC * t_ret / (t_rbc + t_ret)
    rbc0 <- rbase_RBC * t_rbc / (t_rbc + t_ret)
    nor0 <- ret0 * t_ep / t_ret
    cfue0 <- ret0 * t_ep / (t_ep * 2^mnor)
    bfue0 <- ret0 * t_ep / (t_ret * 2^mnor * 2^mcfu)
    hspc0 <- bfue0 / (t_ep * kdiff_ery)

    # First-order rate constant for HSPC differentiation into the MK lineage,
    # back-calculated from the platelet baseline (equation 29). Table S3 lists
    # KM = 28.7e-4 /h; recomputing it here from the published equation keeps
    # the platelet chain in exact steady state at time zero (see vignette
    # Errata for the ~4% arithmetic discrepancy in the printed table).
    kdiff_mk <- rbase_PLT / (cf * t_plt * hspc0)
    # Zero-order HSPC production rate (equation 28), evaluated at the healthy
    # steady state, i.e. without the disease factor.
    kin_hspc <- hspc0 * (kdiff_mk + kdiff_ery)
    # Per-compartment MK precursor baseline (equation 27).
    mk0 <- t_mp * rbase_PLT / (cf * t_plt * nchain)

    # Hemoglobin turnover (equations 15-16). The elimination rate of HGB was
    # assumed to match that of RBC, so KOUT_HGB = 1 / T_RBC and
    # KIN_HGB = HGB0 / T_RBC.
    kin_hb <- rbase_hb / t_rbc
    kout_hb <- 1 / t_rbc

    # -------------------------------------------------------------------
    # Iron driver
    # -------------------------------------------------------------------
    # C_Iron is the CHANGE in serum iron concentration from its baseline, and
    # the baseline is the PK steady state A1(0)/V1 = KIN_Iron / CL
    # = 0.000197 / 0.0001 = 1.97 mg/L = 197 ug/dL. This is confirmed by the
    # Results text: an absolute serum iron of 201.284 ug/dL stimulates HSPCs,
    # and 201.284 - 197.0 = 4.28 ug/dL = Cutoff_Iron in Table 1.
    c_iron_base <- mgL_to_ugdL * kin_iron / cl
    c_iron <- mgL_to_ugdL * central / vc - c_iron_base

    # On/off gate (equation 9): only iron above the cutoff drives HSPCs toward
    # the erythroid lineage.
    fs <- (c_iron > cutoff_iron)

    # Disease progression, partially corrected by supplemental iron, and the
    # gated Emax stimulation of erythroid commitment (equations 8 and 10).
    dfprog_eff <- 1 - dfprog * (1 - smax_df * c_iron / (sc50_df + c_iron))
    stim_ery <- 1 + smax_iron * c_iron / (sc50_iron + c_iron) * fs
    flux_ery <- kdiff_ery * dfprog_eff * prol * stim_ery

    # -------------------------------------------------------------------
    # Serum-iron PK (equations 1-2)
    # -------------------------------------------------------------------
    d/dt(central) <- kin_iron - kcp * central + kpc * peripheral1 - kel * central
    d/dt(peripheral1) <- kcp * central - kpc * peripheral1

    # -------------------------------------------------------------------
    # Erythroid lineage (equations 8, 10-14)
    # -------------------------------------------------------------------
    d/dt(prol) <- kin_hspc - flux_ery - kdiff_mk * prol
    d/dt(bfue) <- flux_ery - bfue / t_ep
    d/dt(cfue) <- 2^mcfu * bfue / t_ep - cfue / t_ep
    d/dt(nor) <- 2^mnor * cfue / t_ep - nor / t_ep
    d/dt(ret) <- nor / t_ep - ret / t_ret
    d/dt(rbc) <- ret / t_ret - rbc / t_rbc

    # -------------------------------------------------------------------
    # Hemoglobin (equation 15)
    # -------------------------------------------------------------------
    d/dt(hb) <- kin_hb * (1 + smax_hgb * c_iron / (sc50_hgb + c_iron)) - kout_hb * hb

    # -------------------------------------------------------------------
    # Megakaryocytic precursor aging chain, n = 10 (equations 17-18)
    # -------------------------------------------------------------------
    d/dt(mk1) <- kdiff_mk * prol - (nchain / t_mp) * mk1
    d/dt(mk2) <- (nchain / t_mp) * (mk1 - mk2)
    d/dt(mk3) <- (nchain / t_mp) * (mk2 - mk3)
    d/dt(mk4) <- (nchain / t_mp) * (mk3 - mk4)
    d/dt(mk5) <- (nchain / t_mp) * (mk4 - mk5)
    d/dt(mk6) <- (nchain / t_mp) * (mk5 - mk6)
    d/dt(mk7) <- (nchain / t_mp) * (mk6 - mk7)
    d/dt(mk8) <- (nchain / t_mp) * (mk7 - mk8)
    d/dt(mk9) <- (nchain / t_mp) * (mk8 - mk9)
    d/dt(mk10) <- (nchain / t_mp) * (mk9 - mk10)

    # -------------------------------------------------------------------
    # Platelet aging chain, n = 10 (equations 19-20)
    # -------------------------------------------------------------------
    d/dt(plt1) <- cf * (nchain / t_mp) * mk10 - (nchain / t_plt) * plt1
    d/dt(plt2) <- (nchain / t_plt) * (plt1 - plt2)
    d/dt(plt3) <- (nchain / t_plt) * (plt2 - plt3)
    d/dt(plt4) <- (nchain / t_plt) * (plt3 - plt4)
    d/dt(plt5) <- (nchain / t_plt) * (plt4 - plt5)
    d/dt(plt6) <- (nchain / t_plt) * (plt5 - plt6)
    d/dt(plt7) <- (nchain / t_plt) * (plt6 - plt7)
    d/dt(plt8) <- (nchain / t_plt) * (plt7 - plt8)
    d/dt(plt9) <- (nchain / t_plt) * (plt8 - plt9)
    d/dt(plt10) <- (nchain / t_plt) * (plt9 - plt10)

    # -------------------------------------------------------------------
    # Initial conditions -- physiological steady state at time zero
    # (equations 6-7 for the PK; equations 22-29 for the PD)
    # -------------------------------------------------------------------
    central(0) <- kin_iron / kel
    peripheral1(0) <- (kin_iron / kel) * (kcp / kpc)
    prol(0) <- hspc0
    bfue(0) <- bfue0
    cfue(0) <- cfue0
    nor(0) <- nor0
    ret(0) <- ret0
    rbc(0) <- rbc0
    hb(0) <- rbase_hb
    mk1(0) <- mk0
    mk2(0) <- mk0
    mk3(0) <- mk0
    mk4(0) <- mk0
    mk5(0) <- mk0
    mk6(0) <- mk0
    mk7(0) <- mk0
    mk8(0) <- mk0
    mk9(0) <- mk0
    mk10(0) <- mk0
    plt1(0) <- rbase_PLT / nchain
    plt2(0) <- rbase_PLT / nchain
    plt3(0) <- rbase_PLT / nchain
    plt4(0) <- rbase_PLT / nchain
    plt5(0) <- rbase_PLT / nchain
    plt6(0) <- rbase_PLT / nchain
    plt7(0) <- rbase_PLT / nchain
    plt8(0) <- rbase_PLT / nchain
    plt9(0) <- rbase_PLT / nchain
    plt10(0) <- rbase_PLT / nchain

    # -------------------------------------------------------------------
    # Observations
    # -------------------------------------------------------------------
    # Serum iron concentration, reported in ug/dL.
    Cc <- mgL_to_ugdL * central / vc
    # The measured red-cell count includes reticulocytes (see the ret0 / rbc0
    # split above).
    RBC <- ret + rbc
    # Platelets are the total across all 10 aging compartments (equation 21).
    PLT <- plt1 + plt2 + plt3 + plt4 + plt5 + plt6 + plt7 + plt8 + plt9 + plt10
    HGB <- hb

    RBC ~ prop(propSd_RBC)
    HGB ~ add(addSd_HGB)
    PLT ~ prop(propSd_PLT)
  })
}
