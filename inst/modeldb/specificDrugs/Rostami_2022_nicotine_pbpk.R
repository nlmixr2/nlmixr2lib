Rostami_2022_nicotine_pbpk <- function() {
  description <- "PBPK (whole-body, MCSim/deSolve). Nicotine, cotinine and nicotine glucuronide disposition in humans, the parent nicotine PBPK model of Rostami 2022 in its originally published parameterization (Robinson / Teeguarden lineage). Carries eleven flow-limited nicotine tissues plus a nicotine-driven heart-rate feedback on cardiac output, a five-tissue cotinine sub-model and a one-compartment nicotine-glucuronide sub-model. Two administration routes are reproducible from the published sources and are implemented: an intravenous infusion directly into arterial blood (the paper's own primary validation, Figure 11) and a first-order oral / gastrointestinal input with bioavailability. Deterministic typical-value model: the sources report no between-subject variability and no residual-error model. The four product-specific routes of the paper (conventional cigarette, ENDS, vapor inhaler and smokeless tobacco) are NOT implemented because their respiratory-tract deposition fractions and buccal / airway permeation coefficients are unreported zero placeholders in the published MCSim listing and are tabulated in no available source; see the vignette Errata. This is the un-adjusted parent of the Salehi 2025 nicotine-pouch model, which halved the hepatic and urinary clearances and rebuilt the buccal front-end."
  reference <- paste(
    "Rostami AA, Campbell JL, Pithawalla YB, Pourhashem H, Muhammad-Kah RS,",
    "Sarkar MA, Liu J, McKinney WJ, Gentry R, Gogova M. (2022).",
    "A comprehensive physiologically based pharmacokinetic (PBPK) model for nicotine in humans",
    "from using nicotine-containing products with different routes of exposure.",
    "Sci Rep 12:1091. doi:10.1038/s41598-022-05108-y. PMCID PMC8776883.",
    "Author Correction Sci Rep 12:2436, doi:10.1038/s41598-022-06693-8 (bibliography fix only, no parameter impact).",
    "Whole-body disposition structure and chemical parameters trace to",
    "Robinson DE, Balter NJ, Schwartz SL (1992) J Pharmacokinet Biopharm 20:591-609 and",
    "Teeguarden JG et al. (2013) Regul Toxicol Pharmacol 65:12-28.",
    sep = " "
  )
  vignette <- "Rostami_2022_nicotine_pbpk"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    amount = "mg",
    weight = "kg"
  )

  # The intravenous validation arm doses arterial blood directly; the oral arm
  # doses the gastrointestinal compartment. Both sit outside the depot/central
  # default, so the dosing targets are declared explicitly.
  dosing <- c("a_gut", "a_arterial")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every tissue volume, blood flow and clearance.",
        "Rostami 2022 Table 1 uses a 73.0 kg reference adult;",
        "the intravenous validation study (Gourlay and Benowitz, cited by Rostami as reference 26)",
        "dosed 2 ug/kg/min, so the delivered milligram amount is body-weight dependent."
      ),
      source_name = "BW"
    )
  )

  # Every ODE state. `analyte` / `specimen` checked against the Rostami 2022
  # supplementary MCSim listing and Rostami 2022 Tables 1-2.
  compartmentData <- list(
    a_gut = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    a_arterial = list(analyte = "nicotine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_venous = list(analyte = "nicotine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_pulmonary = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_buccal = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_conducting_airway = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_transitional_airway = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_heart = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_skin = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_hepatic = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine = list(analyte = "nicotine", units = "mg", specimen = "urine", verified = TRUE),
    effect = list(analyte = "nicotine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    a_venous_cot = list(analyte = "cotinine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial_cot = list(analyte = "cotinine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_liver_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_hepatic_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine_cot = list(analyte = "cotinine", units = "mg", specimen = "urine", verified = TRUE),
    central_gluc = list(analyte = "nicotine glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    a_urine_gluc = list(analyte = "nicotine glucuronide", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    age_range = "adults",
    weight_median = "73.0 kg (Rostami 2022 Table 1 reference adult)",
    disease_state = "healthy adults; validation cohorts are habitual tobacco users",
    dose_range = paste(
      "intravenous nicotine 2 ug/kg/min over 30 min (Figure 11 validation arm);",
      "oral nicotine (first-order gastrointestinal uptake, FA = 0.67, KA = 1.34 /h)"
    ),
    regions = "United States",
    notes = paste(
      "The PBPK model is not fitted to individual-level data; it is a mechanistic",
      "structural model with parameters adopted from Robinson 1992, Teeguarden 2013,",
      "ICRP and Corley et al. The intravenous arm (Figure 11) is the paper's own",
      "self-contained validation and needs no dosimetry or permeation front-end.",
      "The four product-specific inhalation / oral-tobacco routes require CFD-derived",
      "deposition fractions and permeation coefficients that are unreported in every",
      "available source and are therefore not implemented; see the vignette Errata."
    )
  )

  ini({
    # =====================================================================
    # GASTROINTESTINAL / ORAL UPTAKE -- Rostami 2022 Table 2
    # A pseudo-physiological gut compartment with first-order uptake of the
    # bioavailable fraction directly into the liver.
    # =====================================================================
    KA1 <- fixed(1.34); label("Oral / gastrointestinal absorption rate constant (1/h)") # Rostami 2022 Table 2 KA = 1.34 (Teeguarden et al.)
    FA <- fixed(0.67); label("Oral bioavailability of swallowed nicotine (unitless)") # Rostami 2022 Table 2 FA = 0.67 (Teeguarden et al.)

    # =====================================================================
    # METABOLISM AND EXCRETION -- Rostami 2022 Table 2 (original,
    # un-adjusted values; Salehi 2025 later halved the four clearances).
    # =====================================================================
    CLMC <- fixed(2.70); label("Hepatic nicotine metabolic clearance constant (L/h/kg)") # Rostami 2022 Table 2 CLMC = 2.70 L/h/kg^0.75. NOTE: applied linearly in body weight, not BW^0.75 -- see model() and vignette Errata
    CLKC <- fixed(0.42); label("Urinary nicotine clearance constant (L/h/kg^0.75)") # Rostami 2022 Table 2 CLKC = 0.42
    CLLMC <- fixed(0.14); label("Hepatic cotinine metabolic clearance constant (L/h/kg^0.75)") # Rostami 2022 Table 2 CLLMC = 0.14
    CLKMC <- fixed(0.025); label("Urinary cotinine clearance constant (L/h/kg^0.75)") # Rostami 2022 Table 2 CLKMC = 0.025
    FNC <- fixed(0.80); label("Fraction of metabolized nicotine converted to cotinine (unitless)") # Rostami 2022 Table 2 FNC = 0.80 (the supplementary MCSim listing instead carries FNC = 0.72; the published table value is adopted -- see vignette Errata)
    FNG <- fixed(0.10); label("Fraction of metabolized nicotine converted to nicotine glucuronide (unitless)") # Rostami 2022 supplementary MCSim listing: FNG = 0.10
    VGDC <- fixed(1.0); label("Nicotine glucuronide distribution volume constant (L/kg)") # Rostami 2022 supplementary MCSim listing: VGDC = 1.0
    GCLC <- fixed(0.1); label("Urinary nicotine glucuronide clearance constant (L/h/kg^0.75)") # Rostami 2022 supplementary MCSim listing: GCLC = 0.1

    # =====================================================================
    # PHYSIOLOGY -- Rostami 2022 Table 1
    # =====================================================================
    QCC <- fixed(16); label("Cardiac output constant (L/h/kg^0.75)") # Rostami 2022 Table 1 QCC = 16 (ICRP)
    HRO <- fixed(61.1); label("Basal heart rate (beats/min)") # Rostami 2022 Table 1 HRO = 61.1 (Teeguarden et al.)

    VHC <- fixed(0.0044); label("Heart volume (fraction of body weight)") # Rostami 2022 Table 1 VHC
    VBC <- fixed(0.02); label("Brain volume (fraction of body weight)") # Rostami 2022 Table 1 VBC
    VFC <- fixed(0.258); label("Fat volume (fraction of body weight)") # Rostami 2022 Table 1 VFC
    VLC <- fixed(0.024); label("Liver volume (fraction of body weight)") # Rostami 2022 Table 1 VLC
    VSKC <- fixed(0.042); label("Skin volume (fraction of body weight)") # Rostami 2022 Table 1 VSKC
    VMC <- fixed(0.34); label("Muscle volume (fraction of body weight)") # Rostami 2022 Table 1 VMC
    VABC <- fixed(0.02); label("Arterial blood volume (fraction of body weight)") # Rostami 2022 Table 1 VABC
    VVBC <- fixed(0.05); label("Venous blood volume (fraction of body weight)") # Rostami 2022 Table 1 VVBC
    VRC <- fixed(0.03); label("Rapidly perfused volume (fraction of body weight)") # Rostami 2022 Table 1 VRC
    VSLOWC <- fixed(0.08); label("Slowly perfused volume (fraction of body weight)") # Rostami 2022 Table 1 VSLOWC

    QFC <- fixed(0.068); label("Fat blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QFC
    QBC <- fixed(0.12); label("Brain blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QBC
    QHC <- fixed(0.04); label("Heart blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QHC
    QSKC <- fixed(0.05); label("Skin blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QSKC
    QMC <- fixed(0.14); label("Muscle blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QMC
    QLC <- fixed(0.26); label("Liver blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QLC
    QRC <- fixed(0.19); label("Rapidly perfused blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QRC
    QSC <- fixed(0.08); label("Slowly perfused blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QSC
    QBUC <- fixed(0.0215); label("Buccal cavity blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QBUC (Corley et al.)
    QCAC <- fixed(0.025); label("Conducting airway blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QCAC (Corley et al.)
    QTAC <- fixed(0.007); label("Transitional airway blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QTAC (Corley et al.)

    SABU <- fixed(103.10); label("Buccal cavity surface area (cm^2)") # Rostami 2022 Table 1 SABU
    SACA <- fixed(199.50); label("Conducting airway surface area (cm^2)") # Rostami 2022 Table 1 SACA
    SATA <- fixed(163.60); label("Transitional airway surface area (cm^2)") # Rostami 2022 Table 1 SATA
    SAPUL <- fixed(540000); label("Pulmonary surface area (cm^2)") # Rostami 2022 Table 1 SAPUL
    WTBU <- fixed(0.0065); label("Buccal cavity epithelium width (cm)") # Rostami 2022 Table 1 WTBU
    WTCA <- fixed(0.0065); label("Conducting airway epithelium width (cm)") # Rostami 2022 Table 1 WTCA
    WTTA <- fixed(0.0065); label("Transitional airway epithelium width (cm)") # Rostami 2022 Table 1 WTTA
    WTPUL <- fixed(0.000036); label("Pulmonary epithelium width (cm)") # Rostami 2022 Table 1 WTPUL

    # =====================================================================
    # PARTITION COEFFICIENTS -- Rostami 2022 Table 2
    # =====================================================================
    PLU <- fixed(0.90); label("Nicotine lung:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PLU
    PF <- fixed(0.80); label("Nicotine fat:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PF
    PBR <- fixed(3.00); label("Nicotine brain:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PBR
    PL <- fixed(7.50); label("Nicotine liver:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PL
    PH <- fixed(1.60); label("Nicotine heart:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PH
    PSK <- fixed(1.50); label("Nicotine skin:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PSK
    PM <- fixed(1.50); label("Nicotine muscle:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PM
    PR <- fixed(7.50); label("Nicotine rapidly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PR (set to liver)
    PS <- fixed(1.50); label("Nicotine slowly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PS (set to muscle)
    PML <- fixed(2.00); label("Cotinine liver:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PML
    PMM <- fixed(1.50); label("Cotinine muscle:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PMM
    PMR <- fixed(1.50); label("Cotinine rapidly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PMR
    PMS <- fixed(1.00); label("Cotinine slowly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PMS
    PMF <- fixed(0.50); label("Cotinine fat:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PMF

    # =====================================================================
    # SATURABLE TISSUE BINDING -- Rostami 2022 supplementary MCSim listing.
    # Both affinities are zero in the published code, which degenerates the
    # Langmuir isotherm to a hard-threshold sink of capacity BM*V. The
    # capacities are ~1.6e-3 mg (heart) and ~5e-5 mg (lung), i.e. well under
    # 0.1% of a typical nicotine dose. Implemented as coded.
    # =====================================================================
    BMLURC <- fixed(0.00235); label("Lung nicotine binding capacity per unit tissue volume (mg/L)") # Rostami 2022 supplementary MCSim listing: BMLURC = 0.00235
    BMHRC <- fixed(0.00427); label("Heart nicotine binding capacity per unit tissue volume (mg/L)") # Rostami 2022 supplementary MCSim listing: BMHRC = 0.00427
    KBLU <- fixed(0); label("Lung nicotine binding affinity (mg/L)") # Rostami 2022 supplementary MCSim listing: KBLU = 0
    KBH <- fixed(0); label("Heart nicotine binding affinity (mg/L)") # Rostami 2022 supplementary MCSim listing: KBH = 0

    # =====================================================================
    # HEART-RATE / CARDIAC-OUTPUT FEEDBACK -- Rostami 2022 Table 2
    # =====================================================================
    SPD <- fixed(933.66); label("Nicotine concentration-heart rate slope (beats/min per mg/L)") # Rostami 2022 Table 2 S = 933.66 (Teeguarden et al.)
    KANT <- fixed(1.6617); label("First-order rate of loss of nicotine tolerance (1/h)") # Rostami 2022 Table 2 KANT = 1.6617 (Teeguarden et al.)
    CANT50 <- fixed(0.0152); label("Nicotine tolerance half-maximal concentration (mg/L)") # Rostami 2022 Table 2 CANT50 = 0.0152 (Teeguarden et al.)

    # =====================================================================
    # MOLECULAR WEIGHTS (physical constants; declared as 0 placeholders in
    # the published MCSim listing, whose .in input files were not published)
    # =====================================================================
    MW <- fixed(162.23); label("Nicotine molecular weight (g/mol)") # physical constant; PubChem CID 89594. Declared MW = 0 in the published MCSim listing
    MWCOT <- fixed(178); label("Cotinine molecular weight (g/mol)") # Rostami 2022 supplementary MCSim listing: hard-coded literal 178 in the cotinine formation term
    MWGLU <- fixed(338); label("Nicotine glucuronide molecular weight (g/mol)") # Rostami 2022 supplementary MCSim listing: hard-coded literal 338 in the glucuronide formation term (the MWG = 338.4 declaration is not the value used)
  })

  model({
    # =====================================================================
    # 1. TISSUE VOLUMES (L)
    #    Epithelial volumes are surface area x width, converted cm^3 -> L.
    #    Rostami 2022 renormalizes every body-weight fraction by VTOTALC so
    #    the explicit epithelial volumes do not double-count body mass.
    #    The nasal and parallel perfused-lung compartments of the parent
    #    model are omitted: their blood flows QNC and QLNGC are absent from
    #    Rostami 2022 Table 1 (whose flow fractions already sum to ~1.0) and
    #    are 0 in the published code, so both compartments are inert.
    # =====================================================================
    vca <- SACA * WTCA / 1000
    vta <- SATA * WTTA / 1000
    vpul <- SAPUL * WTPUL / 1000
    vbu <- SABU * WTBU / 1000
    vlu <- vca + vta + vpul
    vtotalc <- (vlu + vbu) / WT + VHC + VBC + VFC + VSKC + VMC + VABC + VVBC + VSLOWC + VRC + VLC
    vf <- VFC * WT / vtotalc
    vb <- VBC * WT / vtotalc
    vh <- VHC * WT / vtotalc
    vr <- VRC * WT / vtotalc
    vl <- VLC * WT / vtotalc
    vsk <- VSKC * WT / vtotalc
    vm <- VMC * WT / vtotalc
    vs <- VSLOWC * WT / vtotalc
    vab <- VABC * WT / vtotalc
    vvb <- VVBC * WT / vtotalc
    vgd <- VGDC * WT
    bmh <- BMHRC * vh
    bmlu <- BMLURC * vlu

    # =====================================================================
    # 2. SCALED CLEARANCES (L/h)
    #    CLM is assigned twice in the published MCSim listing -- first as
    #    CLMC*BW^0.75, then four lines later as CLMC*BW. MCSim's generated
    #    C keeps the LAST assignment, so hepatic metabolic clearance is
    #    linear in body weight while every other clearance is allometric.
    #    Reproduced as coded; see vignette Errata.
    # =====================================================================
    clk <- CLKC * WT^0.75
    clm <- CLMC * WT
    clkm <- CLKMC * WT^0.75
    cllm <- CLLMC * WT^0.75
    gcl <- GCLC * WT^0.75

    # =====================================================================
    # 3. CARDIAC OUTPUT WITH NICOTINE HEART-RATE FEEDBACK
    #    Rostami 2022: heart rate rises with venous nicotine and is damped
    #    by an acquired-tolerance state; cardiac output is heart rate x
    #    stroke volume, and every tissue flow is a fixed fraction of it.
    # =====================================================================
    qci <- QCC * WT^0.75
    sv <- qci / (HRO * 60)
    cvvb <- a_venous / vvb
    eout <- HRO + SPD * cvvb / (1 + effect / CANT50)
    qc <- eout * 60 * sv
    qtotalc <- QFC + QBC + QRC + QLC + QHC + QSKC + QMC + QSC + QBUC + QTAC + QCAC
    qf <- QFC * qc / qtotalc
    qb <- QBC * qc / qtotalc
    qr <- QRC * qc / qtotalc
    ql <- QLC * qc / qtotalc
    qh <- QHC * qc / qtotalc
    qsk <- QSKC * qc / qtotalc
    qm <- QMC * qc / qtotalc
    qs <- QSC * qc / qtotalc
    qbu <- QBUC * qc / qtotalc
    qca <- QCAC * qc / qtotalc
    qta <- QTAC * qc / qtotalc

    # Tolerance state: rises toward the current venous concentration and is
    # rectified so tolerance is never lost faster than it is gained
    # (published code: DCANT = (KANT*(CVVB-CANT) < 0) ? 0 : KANT*(CVVB-CANT)).
    d/dt(effect) <- max(KANT * (cvvb - effect), 0)

    # =====================================================================
    # 4. NICOTINE TISSUE CONCENTRATIONS (mg/L)
    #    The lung and heart venous-outflow concentrations invert a Langmuir
    #    binding isotherm; the remaining tissues are flow-limited. Fat uses
    #    a MULTIPLICATION by its partition coefficient where every other
    #    tissue divides -- reproduced as coded; see vignette Errata.
    #    The buccal / conducting-airway / transitional-airway compartments
    #    are perfused tissues carried for flow-balance fidelity; with no
    #    product dosimetry input they simply equilibrate with arterial blood.
    # =====================================================================
    ca <- a_arterial / vab
    ch <- a_heart / vh
    cbr <- a_brain / vb
    cr <- a_rapidly_perfused / vr
    cliv <- a_liver / vl
    csk <- a_skin / vsk
    cmu <- a_muscle / vm
    cslow <- a_slowly_perfused / vs
    cfa <- a_fat / vf

    cvh <- 0.5 * (sqrt(((bmh + vh * PH * KBH - a_heart) / (vh * PH))^2 +
      4 * a_heart * KBH / (vh * PH)) -
      (bmh + vh * PH * KBH - a_heart) / (vh * PH))
    cvpul <- 0.5 * (sqrt(((bmlu + vpul * PLU * KBLU - a_pulmonary) / (vpul * PLU))^2 +
      4 * a_pulmonary * KBLU / (vpul * PLU)) -
      (bmlu + vpul * PLU * KBLU - a_pulmonary) / (vpul * PLU))
    cvbr <- cbr / PBR
    cvr <- cr / PR
    cvl <- cliv / PL
    cvsk <- csk / PSK
    cvm <- cmu / PM
    cvfa <- cfa * PF
    cvs <- cslow / PS
    cvbus <- (a_buccal / vbu) / PLU
    cvcas <- (a_conducting_airway / vca) / PLU
    cvtas <- (a_transitional_airway / vta) / PLU

    # =====================================================================
    # 5. ORAL / GASTROINTESTINAL INPUT
    #    First-order uptake of the bioavailable fraction into the liver.
    #    Oral doses target `a_gut`; bioavailability FA is applied at the
    #    dose so only the absorbed fraction enters, matching the parent
    #    model's ODOSE handling (FA*ODOSE into the gut, all absorbed).
    # =====================================================================
    goral <- KA1 * a_gut
    f(a_gut) <- FA
    d/dt(a_gut) <- -goral

    # =====================================================================
    # 6. NICOTINE MASS BALANCE
    #    Intravenous doses target `a_arterial` directly (Figure 11).
    # =====================================================================
    rametl <- clm * cvl

    d/dt(a_venous) <- qm * cvm + qh * cvh + qf * cvfa + qb * cvbr + qs * cvs +
      qsk * cvsk + qr * cvr + ql * cvl + qta * cvtas + qbu * cvbus +
      qca * cvcas - qc * cvvb
    d/dt(a_arterial) <- qc * (cvpul - ca) - clk * ca
    d/dt(a_urine) <- clk * ca
    d/dt(a_pulmonary) <- qc * (cvvb - cvpul)
    d/dt(a_buccal) <- qbu * (ca - cvbus)
    d/dt(a_conducting_airway) <- qca * (ca - cvcas)
    d/dt(a_transitional_airway) <- qta * (ca - cvtas)
    d/dt(a_liver) <- ql * (ca - cvl) - rametl + goral
    d/dt(a_hepatic) <- rametl
    d/dt(a_heart) <- qh * (ca - cvh)
    d/dt(a_brain) <- qb * (ca - cvbr)
    d/dt(a_rapidly_perfused) <- qr * (ca - cvr)
    d/dt(a_skin) <- qsk * (ca - cvsk)
    d/dt(a_muscle) <- qm * (ca - cvm)
    d/dt(a_fat) <- qf * (ca - cvfa)
    d/dt(a_slowly_perfused) <- qs * (ca - cvs)

    # =====================================================================
    # 7. COTININE SUB-MODEL
    #    Rostami 2022 lumps cotinine into five tissue groups: rapidly
    #    perfused + brain, liver, muscle + heart, fat, and slowly perfused
    #    + skin + buccal + airways. Cotinine is formed in the liver from
    #    the metabolized-nicotine flux, mass-corrected by the cotinine /
    #    nicotine molecular-weight ratio.
    # =====================================================================
    cmr <- a_rapidly_perfused_cot / (vr + vb)
    cml <- a_liver_cot / vl
    cmm <- a_muscle_cot / (vm + vh)
    cmf <- a_fat_cot / vf
    cms <- a_slowly_perfused_cot / (vs + vsk + vbu + vca + vta)
    cvmr <- cmr / PMR
    cvml <- cml / PML
    cvmm <- cmm / PMM
    cvmf <- cmf / PMF
    cvms <- cms / PMS
    cvbm <- a_venous_cot / vvb
    cam <- a_arterial_cot / vab

    d/dt(a_venous_cot) <- (qr + qb) * cvmr + ql * cvml + (qm + qh) * cvmm +
      qf * cvmf + (qs + qsk + qbu + qca + qta) * cvms - qc * cvbm
    d/dt(a_arterial_cot) <- qc * (cvbm - cam) - clkm * cam
    d/dt(a_urine_cot) <- clkm * cam
    d/dt(a_rapidly_perfused_cot) <- (qr + qb) * (cam - cvmr)
    d/dt(a_liver_cot) <- ql * (cam - cvml) + FNC * rametl * (MWCOT / MW) - cllm * cvml
    d/dt(a_hepatic_cot) <- cllm * cvml
    d/dt(a_muscle_cot) <- (qm + qh) * (cam - cvmm)
    d/dt(a_fat_cot) <- qf * (cam - cvmf)
    d/dt(a_slowly_perfused_cot) <- (qs + qsk + qbu + qca + qta) * (cam - cvms)

    # =====================================================================
    # 8. NICOTINE GLUCURONIDE SUB-MODEL (one compartment, renal clearance)
    # =====================================================================
    cgb <- central_gluc / vgd
    d/dt(central_gluc) <- FNG * rametl * (MWGLU / MW) - gcl * cgb
    d/dt(a_urine_gluc) <- gcl * cgb

    # =====================================================================
    # 9. OUTPUTS -- mg/L converted to ng/mL. Rostami 2022 reports venous
    #    plasma nicotine; the arterial output reproduces the arterial /
    #    venous separation of Figure 11. No residual-error model is
    #    reported in any available source, so none is declared.
    # =====================================================================
    Cc <- cvvb * 1000
    Cart <- ca * 1000
    Cc_cot <- cvbm * 1000
    Cc_gluc <- cgb * 1000
    HR <- eout
  })
}
