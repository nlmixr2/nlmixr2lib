Itohara_2019_tacrolimus_pbpk <- function() {
  description <- paste(
    "PBPK (minimal, Simcyp v17 parameterisation re-encoded in NONMEM 7.2).",
    "Tacrolimus whole-blood pharmacokinetics in adult Japanese living-donor",
    "liver transplant (LDLT) recipients during the first postoperative month.",
    "One-compartment disposition (Vd = 1.29 L/kg) with first-order oral",
    "absorption; clearance is a well-stirred hepatic extraction of an",
    "enzyme-abundance intrinsic clearance (CYP3A4 + CYP3A5, 13-O-desmethylation",
    "and 12-hydroxylation) scaled by a regenerating graft liver",
    "(60% of standard size at surgery, half-recovery in one postoperative day).",
    "Oral bioavailability is Fa x Fg x Fh with a Qgut gut-wall model whose",
    "CYP3A Michaelis-Menten extraction saturates on the portal-vein",
    "concentration. Hepatic CYP3A5 follows the DONOR (graft) genotype and",
    "intestinal CYP3A5 follows the RECIPIENT genotype. Deterministic given",
    "covariates (the source control stream fixes both etas to zero)."
  )
  reference <- paste(
    "Itohara K, Yano I, Tsuzuki T, Uesugi M, Nakagawa S, Yonezawa A,",
    "Okajima H, Kaido T, Uemoto S, Matsubara K. A Minimal Physiologically-Based",
    "Pharmacokinetic Model for Tacrolimus in Living-Donor Liver Transplantation:",
    "Perspectives Related to Liver Regeneration and the Cytochrome P450 3A5",
    "(CYP3A5) Genotype. CPT Pharmacometrics Syst Pharmacol. 2019;8(8):587-595.",
    "doi:10.1002/psp4.12420. Model equations from main-text Eq. 1 and",
    "Supplementary Information Eq. 1-7; drug and system parameters from",
    "Table S2; final NONMEM control stream from the Supplementary Information.",
    "The absorption rate constant (not reported) is taken from",
    "Ji E, Kim MG, Oh JM. Ther Clin Risk Manag. 2018;14:2119-2126,",
    "doi:10.2147/TCRM.S184376 (Ji_2018_tacrolimus.R)."
  )
  vignette <- "Itohara_2019_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales only the volume of distribution (Vd = 1.29 L/kg x WT,",
        "SI Eq. 6). Clearance is NOT weight-scaled in this encoding: the",
        "Simcyp standard liver volume and microsomal protein content enter",
        "only as their back-solved typical product (mppgl_vliver)."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the fraction of cardiac output that perfuses the liver",
        "(0.255 male / 0.280 female) and the portal vein (0.19 / 0.215).",
        "The control stream codes SEX = 1 for male and 0 for female, so",
        "SEXF = 1 - SEX. The effect on clearance is below 1% because",
        "hepatic extraction of tacrolimus is low."
      ),
      source_name = "SEX"
    ),
    POD = list(
      description = "Postoperative day after living-donor liver transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying, fractional days from surgery (POD = 0 at",
        "transplantation). Drives the liver-regeneration model of main-text",
        "Eq. 1, LW = (GR + (1 - GR) * POD / (POD + 1)) * LWstd, written here",
        "in the division-safe form of the paper's (1 - GR) / (1 + 1/POD)."
      ),
      source_name = "POD"
    ),
    CYP3A5_EXPR_DONOR = list(
      description = "Graft (donor) CYP3A5 expresser status, 1 = *1 carrier, 0 = *3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (donor CYP3A5*3/*3, 'L*3')",
      notes = paste(
        "Sets the hepatic CYP3A5 abundance: 17 pmol/mg microsomal protein",
        "for a *1-carrying graft ('L*1'), 0 for a *3/*3 graft ('L*3')."
      ),
      source_name = "L*1 / L*3 (L3A5 data item)"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (intestine), 1 = *1 carrier, 0 = *3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (recipient CYP3A5*3/*3, 'I*3')",
      notes = paste(
        "Sets the whole-gut CYP3A5 abundance (G3A5 data item): the",
        "back-solved expresser value for a *1-carrying recipient ('I*1'),",
        "0 for a *3/*3 recipient ('I*3')."
      ),
      source_name = "I*1 / I*3 (G3A5 data item)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 298L,
    n_studies = 2L,
    age_range = "adults (20 years and older)",
    age_mean = "53.5 years (model building), 52.2 years (verification)",
    weight_mean = "66.6 kg (model building), 60.1 kg (verification)",
    disease_state = "Adult living-donor liver transplant recipients during the first postoperative month",
    dose_range = paste(
      "Intravenous infusion 1.25 ug/h/kg for 12 h then 0.83 ug/h/kg,",
      "switched on POD 3 to oral 0.04 mg/kg/day (Kyoto University Hospital",
      "2013 protocol, used for both normalisation and simulation)"
    ),
    regions = "Japan (Kyoto University Hospital)",
    race_ethnicity = c(Asian = 100),
    cyp3a5_genotypes = paste(
      "Model building (n = 26): donor *3/*3 / *1/*3 / *1/*1 = 14 / 11 / 1,",
      "recipient 15 / 9 / 2. Verification (n = 272): donor 166 / 92 / 14,",
      "recipient 171 / 88 / 13 (Table S1)."
    ),
    notes = paste(
      "Model building: 26 LDLT patients (Dec 2013 - Jan 2016) with routine",
      "TDM concentrations during intravenous infusion and oral troughs.",
      "Verification: 272 LDLT patients (Jul 2004 - Jun 2011), trough",
      "concentrations over the first 28 postoperative days (Table S1).",
      "The model itself is a bottom-up Simcyp parameterisation adapted to",
      "these data (liver-regeneration form and Fa selected by fit); the",
      "paper's simulations draw 100 Simcyp Japanese virtual adults per",
      "genotype combination. Mean hematocrit 25.5%."
    )
  )

  ini({
    # ---- Absorption --------------------------------------------------------
    # Not reported by Itohara 2019 (the SKA data item is a Simcyp output).
    # Borrowed per maintainer ruling 2026-09-26 from Ji_2018_tacrolimus.R, a
    # literature-fixed ka for adult LDLT recipients; the paper's own data are
    # troughs and do not identify ka.
    lka <- fixed(log(4.48)); label("First-order absorption rate constant (1/h)") # NOT IN PAPER: Ji 2018 Table 2 ka fixed 4.48 1/h (maintainers' ruling 2026-09-26)
    fa <- fixed(0.6); label("Fraction of the oral dose absorbed (unitless)") # Results 'PBPK modeling' Fa = 0.6 selected (smallest RMSE 5.69); Results 'Simulation' Fa designated 0.6

    # ---- Distribution ------------------------------------------------------
    lvc <- fixed(log(1.29)); label("Blood volume of distribution per kg body weight (L/kg)") # Table S2 'Vb (L/kg) 1.29'; Methods 'set as 1.29 L/kg'; SI Eq. 6

    # ---- Blood binding -----------------------------------------------------
    fu <- fixed(0.013); label("Fraction unbound in plasma (unitless)") # Table S2 'fu 0.013'
    bpr <- fixed(35); label("Blood-to-plasma concentration ratio (unitless)") # Table S2 'BP ratio 35'
    fu_gut <- fixed(1); label("Fraction unbound in the gut wall (unitless)") # Table S2 'fugut 1'; Methods and Figure S1 sensitivity analysis

    # ---- Recombinant-CYP kinetics (Table S2, Dai 2006) ----------------------
    isef <- fixed(0.24); label("Inter-system extrapolation factor rhCYP to HLM (unitless)") # Table S2 'ISEF 0.24' (all four pathways)
    vmax_cyp3a4_13odm <- fixed(8); label("CYP3A4 Vmax, 13-O-desmethylation (pmol/min/pmol CYP)") # Table S2; control stream CLIL '8/0.21'
    km_cyp3a4_13odm <- fixed(0.21); label("CYP3A4 Km, 13-O-desmethylation (umol/L)") # Table S2; control stream CLIL '8/0.21'
    vmax_cyp3a4_12oh <- fixed(0.6); label("CYP3A4 Vmax, 12-hydroxylation (pmol/min/pmol CYP)") # Table S2; control stream CLIL '0.6/0.29'
    km_cyp3a4_12oh <- fixed(0.29); label("CYP3A4 Km, 12-hydroxylation (umol/L)") # Table S2; control stream CLIL '0.6/0.29'
    vmax_cyp3a5_13odm <- fixed(17); label("CYP3A5 Vmax, 13-O-desmethylation (pmol/min/pmol CYP)") # Table S2; control stream CLIL '17/0.21'
    km_cyp3a5_13odm <- fixed(0.21); label("CYP3A5 Km, 13-O-desmethylation (umol/L)") # Table S2; control stream CLIL '17/0.21'
    vmax_cyp3a5_12oh <- fixed(1.4); label("CYP3A5 Vmax, 12-hydroxylation (pmol/min/pmol CYP)") # Table S2; control stream CLIL '1.4/0.35'
    km_cyp3a5_12oh <- fixed(0.35); label("CYP3A5 Km, 12-hydroxylation (umol/L)") # Table S2; control stream CLIL '1.4/0.35'

    # ---- Liver -------------------------------------------------------------
    abund_cyp3a4_liver <- fixed(93); label("Hepatic CYP3A4 abundance (pmol/mg microsomal protein)") # Table S2 'Abundance of CYP3A4 in the liver 93'; Methods 'PBPK model development'
    abund_cyp3a5_liver <- fixed(17); label("Hepatic CYP3A5 abundance in a *1-carrying graft (pmol/mg microsomal protein)") # Table S2 'Abundance of CYP3A5 in the liver 17'; nonexpressor graft = 0 (back-solve, see vignette)
    mppgl_vliver <- fixed(56.04); label("Microsomal protein per g liver x standard liver volume (mg/g x L)") # NOT PRINTED: back-solved from Results/Discussion CL 1.59/1.19 (POD 30), 0.975/0.729 (POD 0) and Fh 0.979/0.985
    gr <- fixed(0.6); label("Graft-to-standard liver size ratio at transplantation (unitless)") # Methods 'GR was assumed to be 60%'; control stream LVD = LVR*0.6

    # ---- Blood flows -------------------------------------------------------
    qco <- fixed(251); label("Cardiac output (L/h)") # NOT PRINTED: back-solved from the same CL and Fh targets (hepatic blood flow 64.0 L/h = 0.255 x qco)
    fq_liver_male <- fixed(0.255); label("Hepatic blood flow as a fraction of cardiac output, male (unitless)") # control stream 'B=CO*0.255'
    fq_liver_female <- fixed(0.280); label("Hepatic blood flow as a fraction of cardiac output, female (unitless)") # control stream 'B=CO*0.280'
    fq_pv_male <- fixed(0.19); label("Portal-vein blood flow as a fraction of cardiac output, male (unitless)") # control stream 'C=CO*0.19'
    fq_pv_female <- fixed(0.215); label("Portal-vein blood flow as a fraction of cardiac output, female (unitless)") # control stream 'C=CO*0.215'

    # ---- Gut wall ----------------------------------------------------------
    qgut <- fixed(11.4); label("Qgut, flow for drug delivery to the gut wall (L/h)") # Table S2 'Qgut (L/h) 11.4'
    # The two gut abundances are Simcyp data items the paper does not print.
    # They are back-solved so that this model's saturable Fg at the last
    # standard-regimen oral dose (POD 27.5, 0.02 mg/kg q12h, 60 kg male,
    # averaged over the two graft genotypes) equals the printed Fg; the
    # control stream tabulates the saturable FGA, not the linear FG.
    abund_cyp3a4_gut <- fixed(76585); label("Whole-gut CYP3A4 abundance (pmol)") # NOT PRINTED: back-solved from Results 'Fg ... I*3 0.227' (vignette Errata)
    abund_cyp3a5_gut <- fixed(20283); label("Whole-gut CYP3A5 abundance in a *1-carrying recipient (pmol)") # NOT PRINTED: back-solved from Results 'Fg ... I*1 0.153' (vignette Errata)

    # ---- Variability -------------------------------------------------------
    # The control stream's '$OMEGA 0 FIXED 0 FIXED' sets both etas (on CL and
    # V) to zero, so they are omitted: the model is deterministic given the
    # covariates. The paper's between-subject spread comes from Simcyp
    # virtual-population draws of the physiological inputs, which it does not
    # report beyond the hepatic CYP3A4 / CYP3A5 abundance CVs (81% / 185%).
    expSd <- fixed(0.3162); label("Log-normal residual error (SD on log scale)") # control stream '$SIGMA 0.1' with Y=F*EXP(EPS(1)); sqrt(0.1)
  })

  model({
    # Tacrolimus molecular weight (g/mol), Table S2; converts ng/mL to umol/L.
    mw <- 804.02
    # Liver density (g/L) converting liver volume to weight, and the unit
    # factor 1e-6 * 60 (pmol/min per umol/L -> L/h); control stream CLIL.
    rho_liver <- 1080
    ufac <- 0.00006

    # ---- Individual parameters ---------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc) * WT

    # Hepatic CYP3A5 is donor-derived (graft), gut CYP3A5 recipient-derived.
    l3a5 <- abund_cyp3a5_liver * CYP3A5_EXPR_DONOR
    g3a5 <- abund_cyp3a5_gut * CYP3A5_EXPR

    # Liver regeneration, main-text Eq. 1 (selected model).
    f_regen <- gr + (1 - gr) * POD / (POD + 1)

    # Hepatic intrinsic clearance, SI Eq. 1 (L/h).
    cint_3a4 <- vmax_cyp3a4_13odm / km_cyp3a4_13odm + vmax_cyp3a4_12oh / km_cyp3a4_12oh
    cint_3a5 <- vmax_cyp3a5_13odm / km_cyp3a5_13odm + vmax_cyp3a5_12oh / km_cyp3a5_12oh
    clint_h <- rho_liver * isef * (cint_3a4 * abund_cyp3a4_liver + cint_3a5 * l3a5) *
      mppgl_vliver * f_regen * ufac
    clint_hb <- fu / bpr * clint_h

    qh <- qco * (fq_liver_male + SEXF * (fq_liver_female - fq_liver_male))
    qpv <- qco * (fq_pv_male + SEXF * (fq_pv_female - fq_pv_male))

    # Well-stirred liver, SI Eq. 4 and 7.
    cl <- qh * clint_hb / (qh + clint_hb)
    fh <- qh / (qh + clint_hb)

    # Gut wall, SI Eq. 2 and 5, with the control stream's two-step
    # fixed-point refinement of the saturable extraction: each step replaces
    # Vmax/Km by Vmax/(Km + C) at the portal-vein concentration C (umol/L),
    # which is the absorbed-flux term at the dose record plus the systemic
    # blood concentration. podo(depot) is the oral dose amount (mg). The
    # control stream's exp(-KA*TAD) factor on the absorbed-flux term is 1 at
    # the dose record (TAD = 0), where f(depot) is evaluated, so it is omitted.
    dose_po <- podo(depot)
    if (is.na(dose_po) || dose_po <= 0) {
      dose_po <- 0
    }
    c_sys <- central / vc * 1000 / mw
    clint_g <- isef * (cint_3a4 * abund_cyp3a4_gut + cint_3a5 * g3a5) * ufac
    fg_lin <- qgut / (qgut + fu_gut * clint_g)
    c_pv1 <- fa * fg_lin * dose_po * ka / qpv * 1000 / mw + c_sys
    clint_g1 <- isef * ((vmax_cyp3a4_13odm / (km_cyp3a4_13odm + c_pv1) +
      vmax_cyp3a4_12oh / (km_cyp3a4_12oh + c_pv1)) * abund_cyp3a4_gut +
      (vmax_cyp3a5_13odm / (km_cyp3a5_13odm + c_pv1) +
        vmax_cyp3a5_12oh / (km_cyp3a5_12oh + c_pv1)) * g3a5) * ufac
    fg1 <- qgut / (qgut + fu_gut * clint_g1)
    c_pv2 <- fa * fg1 * dose_po * ka / qpv * 1000 / mw + c_sys
    clint_g2 <- isef * ((vmax_cyp3a4_13odm / (km_cyp3a4_13odm + c_pv2) +
      vmax_cyp3a4_12oh / (km_cyp3a4_12oh + c_pv2)) * abund_cyp3a4_gut +
      (vmax_cyp3a5_13odm / (km_cyp3a5_13odm + c_pv2) +
        vmax_cyp3a5_12oh / (km_cyp3a5_12oh + c_pv2)) * g3a5) * ufac
    fg <- qgut / (qgut + fu_gut * clint_g2)

    # Oral bioavailability, F = Fa x Fg x Fh.
    f_oral <- fa * fg * fh

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- f_oral

    # Whole-blood concentration (ng/mL): dose mg, volume L.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
