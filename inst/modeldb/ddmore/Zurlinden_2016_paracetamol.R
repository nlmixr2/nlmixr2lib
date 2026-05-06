Zurlinden_2016_paracetamol <- function() {
  description <- "Whole-body physiologically-based pharmacokinetic (PBPK) model for paracetamol (acetaminophen, APAP) and its conjugated metabolites APAP-glucuronide (AG) and APAP-sulfate (AS) in healthy adults (Zurlinden & Reisfeld 2016, DDMODEL00000237 Scenario 4 = 1000 mg single oral dose). Each chemical is distributed across nine flow-limited tissue compartments (fat, kidney, muscle, rapidly perfused, slowly perfused, liver, arterial blood, venous blood) with a separate hepatic sub-compartment for the conjugates. Liver metabolism uses Michaelis-Menten kinetics with partial substrate inhibition for CYP-mediated NAPQI formation, sulfation by SULT (cofactor PAPS), and glucuronidation by UGT (cofactor UDP-glucuronic acid, GA); both cofactors are tracked as relative-amount states with zeroth-order resynthesis. Renal elimination is linear, scaled by body weight. Oral absorption is encoded as a bi-exponential gastric-emptying rate function (Tg, Tp) added directly to the liver compartment, with dose-dependent bioavailability fa = 0.0005*Dose_mg + 0.37 (Dose < 1000 mg) or 0.88 (>= 1000 mg). The model is deterministic typical-value (no IIV, no residual error) — the DDMORE bundle exposes only the Bayesian posterior-mean parameters from Forward_APAP1.in, not the per-individual variability or measurement-error distributions reported in the publication."
  reference <- paste(
    "Zurlinden TJ, Reisfeld B. (2016).",
    "Physiologically based modeling of the pharmacokinetics of acetaminophen and its major metabolites in humans using a Bayesian population approach.",
    "Eur J Drug Metab Pharmacokinet 41(3):267-280.",
    "doi:10.1007/s13318-015-0253-x. PMID 25636597.",
    "DDMORE Foundation Model Repository: DDMODEL00000237 (Scenario 4 = 1000 mg single oral dose).",
    sep = " "
  )
  vignette <- "Zurlinden_2016_paracetamol"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mcg/L",
    amount        = "mcmol",
    weight        = "kg"
  )
  ddmore_id    <- "DDMODEL00000237"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to scale all tissue volumes (linear in WT) and the cardiac output / metabolic Vmax (allometric exponent 0.75) and the renal clearance (linear in WT).  The DDMORE bundle's Forward_APAP1.in fixes BW = 70 kg; users can override via rxSolve(params = c(WT = ...)) or supply WT as a per-subject column in the event table.  The DDMORE bundle does NOT report covariates beyond body weight; the Bayesian fit was on a small adult cohort and the model is intended for typical-value adult simulation.",
      source_name        = "BW"
    )
  )
  # NB: `OralDose_APAP_mg` is declared in ini() (default = 1000 mg, the bundle
  # Scenario 4 dose). Users override at simulation time via
  # `rxSolve(params = c(OralDose_APAP_mg = ...))`. It is a typical-value scalar
  # rather than a per-subject data-column covariate, so it is intentionally NOT
  # listed in covariateData.

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "Adult (specific age distribution not in bundle)",
    weight_range   = "70 kg reference (Forward_APAP1.in)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy adult, no co-medication, single 1000 mg oral paracetamol dose. The Bayesian population fit in Zurlinden & Reisfeld (2016) used pooled human PK data from multiple published studies (per Methods of the publication, which is not on disk in this worktree); the DDMORE bundle does not redistribute the underlying dataset. The bundle's Real_APAP_data.csv is a digitisation of plasma concentrations from Jansen et al. (2004) J Pharm Biomed Anal 34:585-593 — a single-study reference dataset re-used by the authors as one validation source, not the full Bayesian inference dataset.",
    dose_range     = "1000 mg single oral dose (Scenario 4)",
    regions        = NA_character_,
    notes          = "Population demographic detail (n, age range, sex distribution, race) is NOT exposed by the DDMORE bundle (Forward_APAP1.in fixes BW=70 kg and reports only the 21 Bayesian posterior-mean parameter values). The original Zurlinden & Reisfeld 2016 publication is not on disk in this worktree, so no cross-check against the publication's Methods is possible. The model is intended for typical-value adult simulation under the Scenario-4 dosing regimen; downstream users wishing to characterise variability must consult the publication directly to obtain the Bayesian posterior distributions for each parameter (these are summarised in the paper but are not exposed by the dpastoor scrape of DDMODEL00000237)."
  )

  ini({
    # ---------------------------------------------------------------------
    # 21 Bayesian posterior-mean parameter estimates from Forward_APAP1.in
    # `# Mean parameters` block (lines 5-25). The bundle reports each on
    # the natural-log estimation scale; the back-transformed linear value
    # is stored in ini() via `exp(<ln value>)` so that model() can use
    # each parameter directly without additional combination expressions
    # (an extra back-transform pass inside model() trips rxode2's
    # mu-reference parser when many population parameters appear in a
    # single derived expression).  The trailing comment carries the bundle
    # line number AND the original ln value for source-trace audits.
    # Forward_APAP1.in does NOT expose IIV or residual-error distributions
    # — see the vignette Errata.
    # ---------------------------------------------------------------------

    # CYP-mediated metabolism (APAP -> NAPQI)
    CYP_VmaxC      <- exp(  0.0918)  ; label("CYP allometric Vmax constant for liver oxidation (mcmol/hr/kg^0.75)")           # Forward_APAP1.in line 14: lnCYP_VmaxC =  0.0918
    CYP_Km         <- exp(  4.8122)  ; label("CYP Michaelis constant (mcmol/L)")                                              # line 6:   lnCYP_Km       =  4.8122

    # SULT-mediated sulfation (APAP -> AS) with PAPS cofactor
    SULT_VmaxC     <- exp(  5.9494)  ; label("SULT allometric Vmax constant for liver sulfation (mcmol/hr/kg^0.75)")          # line 11:  lnSULT_VmaxC   =  5.9494
    SULT_Km_apap   <- exp(  6.7507)  ; label("SULT Km on APAP substrate (mcmol/L)")                                           # line 23:  lnSULT_Km_apap =  6.7507
    SULT_Km_paps   <- exp( -1.1493)  ; label("SULT Km on PAPS cofactor (relative amount, unitless)")                          # line 15:  lnSULT_Km_paps = -1.1493
    SULT_Ki        <- exp(  5.9984)  ; label("SULT partial substrate-inhibition constant (mcmol/L)")                          # line 12:  lnSULT_Ki      =  5.9984

    # UGT-mediated glucuronidation (APAP -> AG) with GA (UDP-glucuronic acid) cofactor
    UGT_VmaxC      <- exp(  8.2277)  ; label("UGT allometric Vmax constant for liver glucuronidation (mcmol/hr/kg^0.75)")     # line 5:   lnUGT_VmaxC    =  8.2277
    UGT_Km         <- exp(  8.2749)  ; label("UGT Km on APAP substrate (mcmol/L)")                                            # line 9:   lnUGT_Km       =  8.2749
    UGT_Km_GA      <- exp( -1.4898)  ; label("UGT Km on GA cofactor (relative amount, unitless)")                             # line 16:  lnUGT_Km_GA    = -1.4898
    UGT_Ki         <- exp( 10.7505)  ; label("UGT partial substrate-inhibition constant (mcmol/L)")                           # line 22:  lnUGT_Ki       = 10.7505

    # Hepatic-to-systemic transport of conjugates (Vmax, Km on amount in hepatocyte sub-compartment)
    Vmax_AG        <- exp( 10.9996)  ; label("Vmax for AG transport from hepatocyte to liver-blood (mcmol/hr)")               # line 18:  lnVmax_AG      = 10.9996
    Km_AG          <- exp(  9.6067)  ; label("Km for AG hepatocyte-to-liver-blood transport (mcmol)")                         # line 24:  lnKm_AG        =  9.6067
    Vmax_AS        <- exp( 13.6788)  ; label("Vmax for AS transport from hepatocyte to liver-blood (mcmol/hr)")               # line 21:  lnVmax_AS      = 13.6788
    Km_AS          <- exp(  9.7220)  ; label("Km for AS hepatocyte-to-liver-blood transport (mcmol)")                         # line 13:  lnKm_AS        =  9.7220

    # Cofactor zeroth-order resynthesis rates (relative-amount/hr)
    kPAPS_syn      <- exp(  7.9251)  ; label("PAPS zeroth-order resynthesis rate (1/hr)")                                     # line 7:   lnkPAPS_syn    =  7.9251
    kGA_syn        <- exp(  9.0430)  ; label("GA zeroth-order resynthesis rate (1/hr)")                                       # line 25:  lnkGA_syn      =  9.0430

    # Renal blood clearance constants (L/hr/kg) for APAP and conjugates
    CLC_APAP       <- exp( -4.6564)  ; label("APAP renal blood clearance per body weight (L/hr/kg)")                          # line 17:  lnCLC_APAP     = -4.6564
    CLC_AG         <- exp( -1.9876)  ; label("AG renal blood clearance per body weight (L/hr/kg)")                            # line 19:  lnCLC_AG       = -1.9876
    CLC_AS         <- exp( -2.0404)  ; label("AS renal blood clearance per body weight (L/hr/kg)")                            # line 8:   lnCLC_AS       = -2.0404

    # Bi-exponential gastric-emptying time constants (hours)
    Tg             <- exp( -1.1567)  ; label("Gastric-emptying time constant Tg (hr)")                                        # line 20:  lnTg           = -1.1567
    Tp             <- exp( -2.0559)  ; label("Gastric-emptying time constant Tp (hr)")                                        # line 10:  lnTp           = -2.0559

    # ---------------------------------------------------------------------
    # Physiological constants from Executable_APAP.model (NOT estimated;
    # carried verbatim from the literature defaults the bundle ships).
    # Line numbers refer to Executable_APAP.model.
    # ---------------------------------------------------------------------

    # Cardiac-output allometric constant
    QCC            <- 16.2    ; label("Cardiac-output allometric constant (L/hr/kg^0.75)")  # Executable_APAP.model line 24

    # Tissue volume fractions (fraction of body weight, dimensionless)
    VFC            <- 0.214   ; label("Fat tissue volume fraction (kg/kg)")              # line 27
    VKC            <- 0.0044  ; label("Kidney tissue volume fraction (kg/kg)")           # line 28
    VLC            <- 0.0257  ; label("Liver tissue volume fraction (kg/kg)")            # line 30
    VMC            <- 0.4     ; label("Muscle tissue volume fraction (kg/kg)")           # line 31
    VBLAC          <- 0.0243  ; label("Arterial-blood volume fraction (kg/kg)")          # line 32
    VBLVC          <- 0.0557  ; label("Venous-blood volume fraction (kg/kg)")            # line 33
    VSC            <- 0.185   ; label("Slowly-perfused tissue volume fraction (kg/kg)")  # line 34
    VRC            <- 0.0765  ; label("Rapidly-perfused tissue volume fraction (kg/kg)") # line 35
    # NB: bundle line 29 declares VGC = 0.0144 (gut) — used only inside the VTC
    # normalising sum which collapses to 1.0 exactly, so the gut volume fraction
    # is unused in this nlmixr2 translation and is not declared in this ini()
    # block to avoid a `parameter declared but not used` build-time error.

    # Tissue blood-flow fractions (fraction of cardiac output)
    QFC            <- 0.052   ; label("Fat blood-flow fraction (unitless)")              # line 39
    QKC            <- 0.175   ; label("Kidney blood-flow fraction (unitless)")           # line 40
    QGC            <- 0.181   ; label("Gut blood-flow fraction (unitless; gut joins the liver via QL = QG + QLB)")  # line 41
    QLBC           <- 0.046   ; label("Hepatic-artery blood-flow fraction (unitless)")   # line 42
    QMC            <- 0.191   ; label("Muscle blood-flow fraction (unitless)")           # line 43
    QSC            <- 0.14    ; label("Slowly-perfused blood-flow fraction (unitless)")  # line 44
    QRC            <- 0.215   ; label("Rapidly-perfused blood-flow fraction (unitless)") # line 45

    # Renal-clearance scalars (alpha) per chemical (unitless)
    alpha_APAP     <- 1.0     ; label("APAP renal-clearance scalar (unitless)")  # line 59
    alpha_AS       <- 1.0     ; label("AS renal-clearance scalar (unitless)")    # line 60
    alpha_AG       <- 1.0     ; label("AG renal-clearance scalar (unitless)")    # line 61

    # Blood:plasma ratio (also reused for AG and AS in the bundle's CalcOutputs)
    BP_APAP        <- 0.9     ; label("Blood:plasma ratio for APAP (reused for AG and AS by the bundle)")  # line 63

    # Tissue:blood partition coefficients for APAP (Rodgers/Lukacova predictions).
    # The gut compartment is not a state in this model (gut joins the liver via QL),
    # so PG_APAP from the bundle line 68 is omitted from this ini() block.
    PF_APAP        <- 0.447   ; label("APAP fat:blood partition coefficient (unitless)")              # line 67
    PK_APAP        <- 0.711   ; label("APAP kidney:blood partition coefficient (unitless)")           # line 69
    PL_APAP        <- 0.687   ; label("APAP liver:blood partition coefficient (unitless)")            # line 70
    PM_APAP        <- 0.687   ; label("APAP muscle:blood partition coefficient (unitless)")           # line 71
    PR_APAP        <- 0.676   ; label("APAP rapidly-perfused tissue:blood partition coefficient")     # line 73
    PS_APAP        <- 0.606   ; label("APAP slowly-perfused tissue:blood partition coefficient")      # line 74

    # Tissue:blood partition coefficients for AS (sulfate).  PG_AS line 78 is omitted (gut not a state).
    PF_AS          <- 0.088   ; label("AS fat:blood partition coefficient (unitless)")                # line 77
    PK_AS          <- 0.261   ; label("AS kidney:blood partition coefficient (unitless)")             # line 79
    PL_AS          <- 0.203   ; label("AS liver:blood partition coefficient (unitless)")              # line 80
    PM_AS          <- 0.199   ; label("AS muscle:blood partition coefficient (unitless)")             # line 81
    PR_AS          <- 0.207   ; label("AS rapidly-perfused tissue:blood partition coefficient")       # line 83
    PS_AS          <- 0.254   ; label("AS slowly-perfused tissue:blood partition coefficient")        # line 84

    # Tissue:blood partition coefficients for AG (glucuronide).  PG_AG line 88 omitted.
    # NOTE: bundle line 91 has the linear value `PF_AG = 0.336` paired with a
    # log-line `lnPM_AG = log(0.366)` — the values disagree.  This file uses
    # the linear value 0.336 (which is what the bundle's Initialize{} block
    # picks up via `PM_AG = exp(lnPM_AG)` ... unless the lnPM_AG override on
    # the same line is loaded; the order-of-execution is bundle-dependent).
    # Flagged in the vignette Errata.
    PF_AG          <- 0.128   ; label("AG fat:blood partition coefficient (unitless)")                # line 87
    PK_AG          <- 0.392   ; label("AG kidney:blood partition coefficient (unitless)")             # line 89
    PL_AG          <- 0.321   ; label("AG liver:blood partition coefficient (unitless)")              # line 90
    PM_AG          <- 0.336   ; label("AG muscle:blood partition coefficient (unitless; bundle .model line 91 is internally inconsistent — see Errata)")  # line 91
    PR_AG          <- 0.364   ; label("AG rapidly-perfused tissue:blood partition coefficient")       # line 93
    PS_AG          <- 0.351   ; label("AG slowly-perfused tissue:blood partition coefficient")        # line 94

    # Molecular weights (g/mol) — unit-conversion factors mcmol -> mcg
    MW_APAP        <- 151.17  ; label("APAP molecular weight (g/mol)")                                # line 15
    MW_AG          <- 327.28  ; label("APAP-glucuronide molecular weight (g/mol)")                    # line 16
    MW_AS          <- 231.22  ; label("APAP-sulfate molecular weight (g/mol)")                        # line 17

    # Default scenario input (Scenario 4 = 1000 mg single oral dose)
    OralDose_APAP_mg <- 1000  ; label("Default oral paracetamol dose (mg); user-overridable via rxSolve(params=...)")  # Forward_APAP1.in line 34

    # ---------------------------------------------------------------------
    # NO IIV (eta) and NO residual-error blocks. The DDMORE bundle's
    # Forward_APAP1.in exposes only the population-mean parameter values;
    # the publication's Bayesian posterior distributions are not exposed by
    # the bundle. This file is therefore a deterministic typical-value
    # mechanistic PBPK simulator. See the vignette Errata.
    # ---------------------------------------------------------------------
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Dose-dependent oral bioavailability fa(Dose).
    #    Bundle Initialize{} block, Executable_APAP.model lines 597-598:
    #      fa = (ODose_APAP_mg<1000 ? 0.0005*ODose_APAP_mg + 0.37 : 0.88)
    #    For Scenario 4 (1000 mg): condition is FALSE so fa = 0.88.
    # ---------------------------------------------------------------------
    fa <- ifelse(OralDose_APAP_mg < 1000,
                 0.0005 * OralDose_APAP_mg + 0.37,
                 0.88)

    # true_dose = fa * Dose_mg * 1000 / MW_APAP  -> total absorbed dose in mcmol
    # Bundle line 599: true_dose = (fa*ODose_APAP_mg)*(1000./MW_APAP)
    true_dose <- fa * OralDose_APAP_mg * 1000 / MW_APAP

    # ---------------------------------------------------------------------
    # 3. Compartment volumes (L). Bundle lines 608-619.
    #    The bundle divides each volume by VTC = sum of tissue fractions.
    #    Numerically, VTC = VFC + VKC + VGC + VLC + VMC + VBLAC + VBLVC +
    #    VRC + VSC = 0.214 + 0.0044 + 0.0144 + 0.0257 + 0.4 + 0.0243 +
    #    0.0557 + 0.0765 + 0.185 = 1.0 exactly, so the division is a
    #    mathematical no-op and is omitted here.  This avoids declaring a
    #    nine-term additive expression of bare population parameters,
    #    which trips rxode2's mu-reference parser
    #    (`.muRefHandlePlus` rejects 2+ bare population parameters in a
    #    single `+` expression as a missed-eta hint).  Verbatim numerical
    #    behaviour is preserved.
    # ---------------------------------------------------------------------
    VF   <- VFC   * WT
    VK   <- VKC   * WT
    VL   <- VLC   * WT
    VM   <- VMC   * WT
    VBLA <- VBLAC * WT
    VBLV <- VBLVC * WT
    VR   <- VRC   * WT
    VS   <- VSC   * WT

    # ---------------------------------------------------------------------
    # 4. Compartment blood flows (L/hr). Bundle lines 622-634.
    #    QC = QCC * BW^0.75 (allometric cardiac output).
    #    The bundle divides each tissue flow by QTC = sum of flow
    #    fractions = QFC + QKC + QGC + QLBC + QMC + QRC + QSC =
    #    0.052 + 0.175 + 0.181 + 0.046 + 0.191 + 0.14 + 0.215 = 1.0
    #    exactly, so the division is a mathematical no-op and is omitted
    #    here (same mu-ref-parser reason as for VTC above).
    #    QL = gut + hepatic-artery (the gut tissue is collapsed into the
    #    liver in the bundle's APAP model, which is why VG / QG / QGC
    #    are declared but the gut compartment is not integrated
    #    separately).
    # ---------------------------------------------------------------------
    QC  <- QCC * WT^0.75
    QF  <- QFC  * QC
    QK  <- QKC  * QC
    QG  <- QGC  * QC
    QM  <- QMC  * QC
    QLB <- QLBC * QC
    QR  <- QRC  * QC
    QS  <- QSC  * QC
    QL  <- QG + QLB

    # ---------------------------------------------------------------------
    # 5. Renal clearances (L/hr) and allometric metabolic Vmax (mcmol/hr).
    #    Bundle lines 638-649.
    # ---------------------------------------------------------------------
    CLR_APAP  <- alpha_APAP * CLC_APAP * WT
    CLR_AS    <- alpha_AS   * CLC_AS   * WT
    CLR_AG    <- alpha_AG   * CLC_AG   * WT
    CYP_Vmax  <- CYP_VmaxC  * WT^0.75
    UGT_Vmax  <- UGT_VmaxC  * WT^0.75
    SULT_Vmax <- SULT_VmaxC * WT^0.75

    # ---------------------------------------------------------------------
    # 6. Tissue concentrations (mcmol/L). Bundle Dynamics{} lines 670-707.
    #    Each tissue's blood-side concentration is the partitioned form
    #    Cv_<tissue> = C_<tissue> / P_<tissue>.
    # ---------------------------------------------------------------------
    # APAP tissues
    CF_APAP  <- a_fat_apap / VF
    CVF_APAP <- CF_APAP / PF_APAP
    CK_APAP  <- a_kid_apap / VK
    CVK_APAP <- CK_APAP / PK_APAP
    CM_APAP  <- a_mus_apap / VM
    CVM_APAP <- CM_APAP / PM_APAP
    CL_APAP  <- a_liv_apap / VL
    CVL_APAP <- CL_APAP / PL_APAP
    CR_APAP  <- a_rap_apap / VR
    CVR_APAP <- CR_APAP / PR_APAP
    CS_APAP  <- a_slo_apap / VS
    CVS_APAP <- CS_APAP / PS_APAP

    # AS tissues
    CF_AS    <- a_fat_as / VF
    CVF_AS   <- CF_AS / PF_AS
    CK_AS    <- a_kid_as / VK
    CVK_AS   <- CK_AS / PK_AS
    CM_AS    <- a_mus_as / VM
    CVM_AS   <- CM_AS / PM_AS
    CL_AS    <- a_liv_as / VL
    CVL_AS   <- CL_AS / PL_AS
    CR_AS    <- a_rap_as / VR
    CVR_AS   <- CR_AS / PR_AS
    CS_AS    <- a_slo_as / VS
    CVS_AS   <- CS_AS / PS_AS

    # AG tissues
    CF_AG    <- a_fat_ag / VF
    CVF_AG   <- CF_AG / PF_AG
    CK_AG    <- a_kid_ag / VK
    CVK_AG   <- CK_AG / PK_AG
    CM_AG    <- a_mus_ag / VM
    CVM_AG   <- CM_AG / PM_AG
    CL_AG    <- a_liv_ag / VL
    CVL_AG   <- CL_AG / PL_AG
    CR_AG    <- a_rap_ag / VR
    CVR_AG   <- CR_AG / PR_AG
    CS_AG    <- a_slo_ag / VS
    CVS_AG   <- CS_AG / PS_AG

    # ---------------------------------------------------------------------
    # 7. Blood and plasma concentrations (mcmol/L). Bundle lines 712-722.
    #    Note: the bundle's APAP arterial / venous formulas use the
    #    OPPOSITE-side blood volume (CA_APAP = ABLA / VBLV; CV_APAP = ABLV
    #    / VBLA); the AS and AG formulas use the same-side volumes
    #    correctly (CA_AS = ABLA_AS / VBLA, CV_AS = ABLV_AS / VBLV). This
    #    is preserved verbatim per the operator's `extract_verbatim`
    #    decision (see vignette Errata "Bundle quirks"); the Bayesian
    #    posterior means in Forward_APAP1.in were obtained against this
    #    same .model so changing the volumes here would invalidate the
    #    parameter values.
    # ---------------------------------------------------------------------
    CA_APAP    <- a_art_apap / VBLV   # bundle line 712 (APAP only — see Errata)
    CV_APAP    <- a_ven_apap / VBLA   # bundle line 713 (APAP only — see Errata)
    Cplasma_apap <- CV_APAP / BP_APAP
    CA_AS      <- a_art_as / VBLA     # bundle line 716
    CV_AS      <- a_ven_as / VBLV     # bundle line 717
    Cplasma_as <- CV_AS / BP_APAP     # bundle line 718 (BP_APAP reused for AS)
    CA_AG      <- a_art_ag / VBLA     # bundle line 720
    CV_AG      <- a_ven_ag / VBLV     # bundle line 721
    Cplasma_ag <- CV_AG / BP_APAP     # bundle line 722 (BP_APAP reused for AG)

    # ---------------------------------------------------------------------
    # 8. Liver metabolic rates (mcmol/hr). Bundle lines 732-734.
    #    Each rate has Michaelis-Menten kinetics with partial substrate
    #    inhibition (Ki) for SULT/UGT and cofactor dependence (PAPS, GA).
    # ---------------------------------------------------------------------
    r_napqi <- CYP_Vmax * CL_APAP / (CYP_Km + CL_APAP)
    r_sult  <- SULT_Vmax * CL_APAP * a_paps /
                ((SULT_Km_apap + CL_APAP + CL_APAP^2 / SULT_Ki) *
                 (SULT_Km_paps + a_paps))
    r_ugt   <- UGT_Vmax  * CL_APAP * a_ga /
                ((UGT_Km   + CL_APAP + CL_APAP^2 / UGT_Ki) *
                 (UGT_Km_GA + a_ga))

    # Hepatocyte-to-liver-blood transport rates for the conjugates
    # (mcmol/hr). Bundle lines 786, 819.
    r_hep_as <- Vmax_AS * a_hep_as / (Km_AS + a_hep_as)
    r_hep_ag <- Vmax_AG * a_hep_ag / (Km_AG + a_hep_ag)

    # ---------------------------------------------------------------------
    # 9. Renal elimination rates (mcmol/hr) and gastric-emptying input
    #    rate (mcmol/hr). Bundle lines 749, 770, 809, 842.
    #    Gastric emptying is encoded as a closed-form function of the
    #    simulation time t (NOT of stomach state), so this model is
    #    parameterised for a single oral dose at t = 0 only — see vignette
    #    Errata.
    # ---------------------------------------------------------------------
    r_renal_apap <- CLR_APAP * CA_APAP
    r_renal_as   <- CLR_AS   * CA_AS
    r_renal_ag   <- CLR_AG   * CA_AG
    gastric_in   <- true_dose * (exp(-t / Tg) - exp(-t / Tp)) / (Tg - Tp)

    # ---------------------------------------------------------------------
    # 10. APAP ODE system (9 states). Bundle lines 753-780.
    #     Liver: receives arterial inflow from QL, the gastric-emptying
    #     rate (collapsed gut/portal/hepatic-artery), and loses mass to
    #     three metabolic pathways and venous outflow.
    #     Other tissues: flow-limited (Q*(CA - Cv_tissue)).
    #     Kidney: flow-limited PLUS a renal-elimination sink (CLR * CA).
    #     Venous blood: sums all tissue venous outflows minus arterial
    #     output to lungs.
    #     Arterial blood: receives venous output, supplies tissues.
    #     Urine: cumulative mass excreted via the kidney.
    # ---------------------------------------------------------------------
    d/dt(a_liv_apap) <- QL * CA_APAP + gastric_in - QL * CVL_APAP - r_napqi - r_sult - r_ugt
    d/dt(a_fat_apap) <- QF * (CA_APAP - CVF_APAP)
    d/dt(a_mus_apap) <- QM * (CA_APAP - CVM_APAP)
    d/dt(a_kid_apap) <- QK * (CA_APAP - CVK_APAP) - r_renal_apap
    d/dt(a_rap_apap) <- QR * (CA_APAP - CVR_APAP)
    d/dt(a_slo_apap) <- QS * (CA_APAP - CVS_APAP)
    d/dt(a_ven_apap) <- QF * CVF_APAP + QM * CVM_APAP + QK * CVK_APAP +
                        QL * CVL_APAP + QR * CVR_APAP + QS * CVS_APAP -
                        QC * CV_APAP
    d/dt(a_art_apap) <- QC * (CV_APAP - CA_APAP)
    d/dt(a_uri_apap) <- r_renal_apap

    # ---------------------------------------------------------------------
    # 11. AS (sulfate) ODE system (10 states). Bundle lines 786-813.
    #     Hepatocyte sub-compartment a_hep_as accumulates AS from
    #     sulfation (r_sult) and empties to liver-blood via Vmax/Km
    #     (r_hep_as). Tissue + blood + urine equations mirror APAP.
    # ---------------------------------------------------------------------
    d/dt(a_hep_as) <- r_sult - r_hep_as
    d/dt(a_liv_as) <- QL * (CA_AS - CVL_AS) + r_hep_as
    d/dt(a_fat_as) <- QF * (CA_AS - CVF_AS)
    d/dt(a_mus_as) <- QM * (CA_AS - CVM_AS)
    d/dt(a_kid_as) <- QK * (CA_AS - CVK_AS) - r_renal_as
    d/dt(a_rap_as) <- QR * (CA_AS - CVR_AS)
    d/dt(a_slo_as) <- QS * (CA_AS - CVS_AS)
    d/dt(a_ven_as) <- QF * CVF_AS + QM * CVM_AS + QK * CVK_AS +
                      QL * CVL_AS + QR * CVR_AS + QS * CVS_AS -
                      QC * CV_AS
    d/dt(a_art_as) <- QC * (CV_AS - CA_AS)
    d/dt(a_uri_as) <- r_renal_as

    # ---------------------------------------------------------------------
    # 12. AG (glucuronide) ODE system (10 states). Bundle lines 819-846.
    # ---------------------------------------------------------------------
    d/dt(a_hep_ag) <- r_ugt - r_hep_ag
    d/dt(a_liv_ag) <- QL * (CA_AG - CVL_AG) + r_hep_ag
    d/dt(a_fat_ag) <- QF * (CA_AG - CVF_AG)
    d/dt(a_mus_ag) <- QM * (CA_AG - CVM_AG)
    d/dt(a_kid_ag) <- QK * (CA_AG - CVK_AG) - r_renal_ag
    d/dt(a_rap_ag) <- QR * (CA_AG - CVR_AG)
    d/dt(a_slo_ag) <- QS * (CA_AG - CVS_AG)
    d/dt(a_ven_ag) <- QF * CVF_AG + QM * CVM_AG + QK * CVK_AG +
                      QL * CVL_AG + QR * CVR_AG + QS * CVS_AG -
                      QC * CV_AG
    d/dt(a_art_ag) <- QC * (CV_AG - CA_AG)
    d/dt(a_uri_ag) <- r_renal_ag

    # ---------------------------------------------------------------------
    # 13. Cofactor dynamics (relative-amount states, 1.0 at baseline).
    #     Bundle lines 737-738. Each cofactor is depleted by its
    #     conjugation reaction and regenerated by zeroth-order synthesis
    #     proportional to the deficit (1 - a_cofactor).
    # ---------------------------------------------------------------------
    d/dt(a_paps) <- -r_sult + kPAPS_syn * (1 - a_paps)
    d/dt(a_ga)   <- -r_ugt  + kGA_syn   * (1 - a_ga)

    # Initial conditions: cofactors at baseline 1.0 (bundle Initialize{}
    # lines 654-655). All tissue and blood states default to 0.
    a_paps(0) <- 1
    a_ga(0)   <- 1

    # ---------------------------------------------------------------------
    # 14. Plasma-concentration outputs in mcg/L (= ng/mL). The bundle's
    #     CalcOutputs{} block emits CPL_<x>_mcgL = CPL_<x> (mcmol/L) *
    #     MW_<x> (g/mol). Bundle lines 858-861.
    # ---------------------------------------------------------------------
    Cplasma_apap_mcgL <- Cplasma_apap * MW_APAP
    Cplasma_ag_mcgL   <- Cplasma_ag   * MW_AG
    Cplasma_as_mcgL   <- Cplasma_as   * MW_AS
  })
}
