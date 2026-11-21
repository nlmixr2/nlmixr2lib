Luu_2001_d4_human <- function() {
  description <- paste(
    "PBPK (whole-body, human). Octamethylcyclotetrasiloxane (D4), a",
    "low-molecular-weight cyclic siloxane used in personal-care products",
    "and silicone breast implants (Luu and Hutter 2001, Environ Health",
    "Perspect). The rat PBPK (Luu_2001_d4_rat) scaled to human physiology",
    "by power law (Y = a*M^0.75). Flow-limited lung, richly-perfused",
    "tissue, a two-subcompartment fat depot, liver with first-order",
    "metabolism, GI tract and kidney, plus blood D4 dissolved in plasma",
    "and weakly protein-bound. Unlike the rat, the human model has no",
    "strongly-bound plasma pool (Table 3 Vstrong = 0). Measured plasma",
    "radioactivity is total D4 (dissolved plus bound). The inhalation",
    "dose is the net absorbed amount delivered to the lung over the",
    "exposure (11.1 mg over 1 h in the validation scenario), with 0.8%",
    "preloaded to the weakly-bound blood pool. Validated against human",
    "inhalation plasma data (Utell 1998). Deterministic, no random",
    "effects."
  )
  reference <- paste(
    "Luu HMD, Hutter JC. Bioavailability of octamethylcyclotetrasiloxane",
    "(D4) after exposure to silicones by inhalation and implantation.",
    "Environ Health Perspect. 2001;109(11):1095-1101.",
    "doi:10.1289/ehp.109-1240468. Calibrated parameters from Table 3",
    "(human column); physiologic volumes and flows from Table 2 (human);",
    "material-balance ODEs from the Appendix."
  )
  vignette <- "Luu_2001_d4_siloxane"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  paper_specific_compartments <- c(
    "Clung",
    "Cr",
    "Cf1",
    "Cf2",
    "Cl",
    "Cgi",
    "Ck",
    "Cai",
    "Cwk",
    "Mmet"
  )

  compartmentData <- list(
    Clung = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cr = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cf1 = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cf2 = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cl = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cgi = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Ck = list(analyte = "D4", units = "ug/mL", specimen = "tissue", verified = TRUE),
    Cai = list(analyte = "D4", units = "ug/mL", specimen = "plasma", verified = TRUE),
    Cwk = list(analyte = "D4", units = "ug/mL", specimen = "plasma", verified = TRUE),
    Mmet = list(analyte = "D4 metabolites", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  # Body weight is not a model covariate: the tissue volumes and blood
  # flows are the fixed physiologic values of Table 2 (human), and the
  # dose is supplied directly as DOSE (ug). WT is documented here only to
  # record the 70 kg adult the Table 3 scale-up factors reference.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Reference 70 kg adult for the Table 3 scale-up factors. Not used",
        "in model(): physiology is the fixed Table 2 (human) values and",
        "the dose enters as DOSE (ug, net absorbed)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "adult",
    weight_range = "70 kg (Table 3 scale-up reference)",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Healthy adult human volunteers. Validation inhalation data are",
      "from 12 volunteers exposed to 10 ppm D4 vapour for 1 h (Utell",
      "1998). The model also simulates repeated inhalation and silicone",
      "breast-implant exposure, motivated by D4 accumulation in women's",
      "fatty tissues."
    ),
    dose_range = paste(
      "Inhalation: 11.1 mg net absorbed (single) to 77.7 mg (7-day",
      "repeated); breast implant: 10.4 mg leached."
    ),
    regions = NA_character_,
    notes = paste(
      "Toxicokinetic PBPK for an environmental/industrial siloxane, not",
      "a therapeutic drug. Human physiology scaled from the rat",
      "calibration by power law (Table 2 human column; metabolism and",
      "binding parameters from Table 3 human column). All parameters are",
      "FIXED; the model is deterministic."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Human D4 parameters. Source: Luu and Hutter 2001 Table 3 (human
    # column). The human model has NO strongly-bound plasma pool
    # (Vstrong = 0), so the strong-binding state is dropped relative to
    # the rat. kel scales up by the power-law factor (70/0.20)^0.75. All
    # FIXED; deterministic model.
    # -------------------------------------------------------------------
    lDr <- fixed(log(10)); label("Richly-perfused tissue:blood distribution coefficient (unitless)") # Table 3 'Dr' = 10
    lDf1 <- fixed(log(700)); label("Weakly-bound fat:blood distribution coefficient (unitless)") # Table 3 'Df1' = 700
    lDf2 <- fixed(log(2000)); label("Strongly-bound fat:blood distribution coefficient (unitless)") # Table 3 'Df2' = 2,000
    lDl <- fixed(log(7.4)); label("Liver:blood distribution coefficient (unitless)") # Table 3 'Dl' human = 7.4
    lDgi <- fixed(log(0.41)); label("GI tract:blood distribution coefficient (unitless)") # Table 3 'Dgi' = 0.41
    lDk <- fixed(log(40)); label("Kidney:blood distribution coefficient (unitless)") # Table 3 'Dk' = 40
    lHair <- fixed(log(20)); label("Blood:air partition constant (unitless)") # Table 3 'Hair' = 20.0

    lkel <- fixed(log(728)); label("Hepatic metabolism rate constant (1/h)") # Table 3 'kel' human = 728 /h (= 9 x (70/0.20)^0.75)

    lkwi <- fixed(log(0.2)); label("Forward rate constant, weak plasma protein binding (1/h)") # Table 3 'kwi' = 0.2 /h
    lkwo <- fixed(log(1.0)); label("Reverse rate constant, weak plasma protein binding (1/h)") # Table 3 'kwo' = 1.0 /h

    lVweak <- fixed(log(10.4)); label("Volume of weakly-bound plasma lipoprotein pool (mL)") # Table 3 'Vweak' human = 10.4 mL

    lQf2Qf1 <- fixed(log(0.059)); label("Deep-fat exchange flow as a fraction of fat blood flow (Qf2/Qf1)") # Table 3 'Qf2/Qf1' = 0.059
    lVf2Vf1 <- fixed(log(0.73)); label("Deep-fat volume as a fraction of fat volume (Vf2/Vf1)") # Table 3 'Vf2/Vf1' = 0.73

    fwb0 <- fixed(0.008); label("Fraction of absorbed dose preloaded to weakly-bound blood at t=0 (unitless)") # Appendix: human weak-bound D4 initialised to 0.8% of delivered dose

    # Dosing inputs; overridden per scenario in the vignette via rxSolve.
    DOSE <- fixed(11100); label("Net absorbed D4 dose delivered to the lung (ug)") # Appendix: 11.1 mg delivered to the body over the 1-h inhalation
    tinj <- fixed(1.0); label("Duration of the inhalation input (h)") # Utell 1998 1-h exposure

    # No residual-error model is reported for this deterministic PBPK.
    # FIXED placeholder, carries no information from the paper; see Errata.
    propSd <- fixed(0.1); label("Proportional residual error placeholder, blood (fraction)") # not reported in Luu 2001; placeholder only
  })

  model({
    # ===================================================================
    # 0. Back-transforms (Table 3 human parameters)
    # ===================================================================
    Dr <- exp(lDr)
    Df1 <- exp(lDf1)
    Df2 <- exp(lDf2)
    Dl <- exp(lDl)
    Dgi <- exp(lDgi)
    Dk <- exp(lDk)
    Hair <- exp(lHair)
    kel <- exp(lkel)
    kwi <- exp(lkwi)
    kwo <- exp(lkwo)
    Vweak <- exp(lVweak)
    Qf2Qf1 <- exp(lQf2Qf1)
    Vf2Vf1 <- exp(lVf2Vf1)

    # ===================================================================
    # 1. Fixed human physiology (Table 2). Flows mL/min -> mL/h (x 60).
    # ===================================================================
    Qt <- 5200 * 60 # Table 2 human total cardiac output 5,200 mL/min
    Qr <- 2054 * 60 # Table 2 human richly-perfused flow 2,054 mL/min
    Qf1 <- 442 * 60 # Table 2 human fat flow 442 mL/min
    Ql <- 1196 * 60 # Table 2 human liver flow 1,196 mL/min
    Qgi <- 728 * 60 # Table 2 human GI flow 728 mL/min
    Qk <- 780 * 60 # Table 2 human kidney flow 780 mL/min
    Qf2 <- Qf2Qf1 * Qf1 # deep-fat exchange flow (Table 3 ratio)

    Vlung <- 3800 # Table 2 human lung volume 3,800 mL
    Vr <- 40672 # Table 2 human richly-perfused volume 40,672 mL
    Vf1 <- 14562 # Table 2 human fat volume 14,562 mL
    Vf2 <- Vf2Vf1 * Vf1 # deep-fat volume (Table 3 ratio)
    Vl <- 1419 # Table 2 human liver volume 1,419 mL
    Vgi <- 948 # Table 2 human GI volume 948 mL
    Vk <- 243 # Table 2 human kidney volume 243 mL
    Vblood <- 2800 # Table 2 human plasma volume 2,800 mL

    # ===================================================================
    # 2. Arterial concentration and the inhalation input to the lung.
    #    The net absorbed dose DOSE is delivered to the lung over tinj;
    #    tissues receive D4 only via arterial perfusion (no IV forcing).
    # ===================================================================
    Cart <- Hair * Clung
    input <- (t < tinj) * (DOSE / tinj)

    # ===================================================================
    # 3. Material-balance ODEs (Appendix), human (no strong binding).
    #    The same three printed-equation typos corrected as in the rat.
    # ===================================================================
    d / dt(Clung) <- (Qt * Cai - Qt * Cart + input) / Vlung
    d / dt(Cr) <- (Qr * Cart - Qr * Cr / Dr) / Vr
    d / dt(Cf2) <- (Qf2 * Cf1 / Df1 - Qf2 * Cf2 / Df2) / Vf2
    d / dt(Cf1) <-
      (Qf1 * Cart + Qf2 * Cf2 / Df2 -
        Qf1 * Cf1 / Df1 - Qf2 * Cf1 / Df1) / Vf1
    d / dt(Cl) <- (Ql * Cart - Ql * Cl / Dl) / Vl - kel * Cl
    d / dt(Cgi) <- (Qgi * Cart - Qgi * Cgi / Dgi) / Vgi
    d / dt(Ck) <- (Qk * Cart - Qk * Ck / Dk) / Vk

    d / dt(Cai) <-
      (Qr * Cr / Dr + Qf1 * Cf1 / Df1 + Ql * Cl / Dl + Qgi * Cgi / Dgi +
        Qk * Ck / Dk + kwo * Cwk - kwi * Cai - Qt * Cai) / Vblood
    d / dt(Cwk) <- (kwi * Cai - kwo * Cwk) / Vweak

    d / dt(Mmet) <- kel * Cl * Vl

    # Weakly-bound blood is preloaded with 0.8% of the absorbed dose.
    Cwk(0) <- fwb0 * DOSE / Vweak

    # ===================================================================
    # 4. Observations. Measured plasma radioactivity is total D4
    #    (dissolved + weakly protein-bound), reported in ng/mL.
    # ===================================================================
    Cc <- 1000 * (Cai * Vblood + Cwk * Vweak) / Vblood
    Cfat <- 1000 * Cf1
    Cliver <- 1000 * Cl

    Cc ~ prop(propSd)
  })
}
