Luu_2001_d4_rat <- function() {
  description <- paste(
    "PBPK (whole-body, female Sprague-Dawley rat). Octamethyl-",
    "cyclotetrasiloxane (D4), a low-molecular-weight cyclic siloxane",
    "used in personal-care products and as a silicone-polymer precursor",
    "(Luu and Hutter 2001, Environ Health Perspect). Flow-limited",
    "lung, richly-perfused tissue, a two-subcompartment fat depot",
    "(weakly- and strongly-bound), liver with first-order metabolism,",
    "GI tract and kidney, plus a three-state blood model (D4 dissolved",
    "in plasma, weakly protein-bound, strongly protein-bound). The",
    "measured blood/tissue radioactivity is total D4 (dissolved plus",
    "bound). The intravenous dose enters as a flow-partitioned forcing",
    "function distributed to the tissues, with 5% of the dose preloaded",
    "to the weakly-bound blood pool. Calibrated to 14C-D4 IV disposition",
    "in female SD rats (Kirkpatrick 1995). Deterministic, no random",
    "effects. See the Luu_2001_d4_human companion for the human",
    "inhalation/implant scale-up."
  )
  reference <- paste(
    "Luu HMD, Hutter JC. Bioavailability of octamethylcyclotetrasiloxane",
    "(D4) after exposure to silicones by inhalation and implantation.",
    "Environ Health Perspect. 2001;109(11):1095-1101.",
    "doi:10.1289/ehp.109-1240468. Calibrated parameters from Table 3",
    "(rat / Sprague-Dawley column); physiologic volumes and flows from",
    "Table 2 (rat); material-balance ODEs from the Appendix."
  )
  vignette <- "Luu_2001_d4_siloxane"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  # Material-balance states are D4 concentrations (ug/mL); the metabolite
  # pool Mmet is an amount (ug). None map to a canonical compartment role.
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
    "Cstr",
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
    Cstr = list(analyte = "D4", units = "ug/mL", specimen = "plasma", verified = TRUE),
    Mmet = list(analyte = "D4 metabolites", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  # Body weight is not a model covariate: the tissue volumes and blood
  # flows are the fixed physiologic values of Table 2, and the dose is
  # supplied directly as DOSE (ug). WT is documented here only to record
  # the 0.20 kg rat the Table 3 scale-up factors reference.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Reference 0.20 kg female SD rat for the Table 3 scale-up factors.",
        "Not used in model(): physiology is the fixed Table 2 values and",
        "the dose enters as DOSE (ug). The 7 mg/kg scenario is DOSE = 1400",
        "ug (7 mg/kg x 0.20 kg)."
      )
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley, female)",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "adult",
    weight_range = "~0.20 kg (Table 3 scale-up reference)",
    sex_female_pct = 100,
    disease_state = paste(
      "Healthy adult female SD rats. Only the female data were used",
      "(females retain more D4 than males); the calibration study dosed",
      "14C-D4 intravenously at 7 and 70 mg/kg single and 7 mg/kg daily",
      "for 14 days (Kirkpatrick 1995)."
    ),
    dose_range = "7-70 mg/kg 14C-D4 intravenous",
    regions = NA_character_,
    notes = paste(
      "Toxicokinetic PBPK for an environmental/industrial siloxane, not",
      "a therapeutic drug. Physiologic parameters from Table 2 (rat);",
      "calibrated distribution coefficients, protein-binding and",
      "metabolism parameters from Table 3 (SD-rat column). All",
      "parameters are FIXED; the model is deterministic."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Calibrated D4 parameters, female SD rat. Source: Luu and Hutter
    # 2001 Table 3 (rat column, Sprague-Dawley values where the table
    # distinguishes SD from F344). All FIXED; deterministic model. The
    # 'd' / 'e' footnotes in Table 3 mark Sprague-Dawley vs F344; the IV
    # calibration is Sprague-Dawley, so the SD values are used.
    # -------------------------------------------------------------------
    lDr <- fixed(log(10)); label("Richly-perfused tissue:blood distribution coefficient (unitless)") # Table 3 'Dr' = 10
    lDf1 <- fixed(log(700)); label("Weakly-bound fat:blood distribution coefficient (unitless)") # Table 3 'Df1' = 700
    lDf2 <- fixed(log(2000)); label("Strongly-bound fat:blood distribution coefficient (unitless)") # Table 3 'Df2' = 2,000 (estimated from Log P = 5.1)
    lDl <- fixed(log(7.4)); label("Liver:blood distribution coefficient (unitless)") # Table 3 'Dl' SD-rat = 7.4
    lDgi <- fixed(log(0.41)); label("GI tract:blood distribution coefficient (unitless)") # Table 3 'Dgi' = 0.41
    lDk <- fixed(log(40)); label("Kidney:blood distribution coefficient (unitless)") # Table 3 'Dk' = 40
    lHair <- fixed(log(20)); label("Blood:air partition constant (unitless)") # Table 3 'Hair' = 20.0

    lkel <- fixed(log(9)); label("Hepatic metabolism rate constant (1/h)") # Table 3 'kel' SD-rat = 9 /h

    lksi <- fixed(log(0.4)); label("Forward rate constant, strong plasma protein binding (1/h)") # Table 3 'ksi' = 0.4 /h
    lkso <- fixed(log(0.07)); label("Reverse rate constant, strong plasma protein binding (1/h)") # Table 3 'kso' = 0.07 /h
    lkwi <- fixed(log(0.2)); label("Forward rate constant, weak plasma protein binding (1/h)") # Table 3 'kwi' = 0.2 /h
    lkwo <- fixed(log(1.0)); label("Reverse rate constant, weak plasma protein binding (1/h)") # Table 3 'kwo' = 1.0 /h

    lVweak <- fixed(log(1.0)); label("Volume of weakly-bound plasma lipoprotein pool (mL)") # Table 3 'Vweak' rat = 1.0 mL
    lVstrong <- fixed(log(1.0)); label("Volume of strongly-bound plasma lipoprotein pool (mL)") # Table 3 'Vstrong' rat = 1.0 mL

    lQf2Qf1 <- fixed(log(0.059)); label("Deep-fat exchange flow as a fraction of fat blood flow (Qf2/Qf1)") # Table 3 'Qf2/Qf1' = 0.059
    lVf2Vf1 <- fixed(log(0.73)); label("Deep-fat volume as a fraction of fat volume (Vf2/Vf1)") # Table 3 'Vf2/Vf1' = 0.73

    fwb0 <- fixed(0.05); label("Fraction of IV dose preloaded to weakly-bound blood at t=0 (unitless)") # Appendix: weak-bound D4 initialised to 5% of IV dose

    # Dosing inputs; overridden per scenario in the vignette via rxSolve.
    DOSE <- fixed(1400); label("Intravenous D4 dose (ug)") # 7 mg/kg x 0.20 kg = 1.4 mg = 1400 ug (default single-dose scenario)
    tinj <- fixed(0.1); label("Duration of the IV forcing function (h)") # near-bolus; paper reports tmax = 0 for the single IV dose

    # No residual-error model is reported for this deterministic PBPK.
    # FIXED placeholder, carries no information from the paper; see Errata.
    propSd <- fixed(0.1); label("Proportional residual error placeholder, blood (fraction)") # not reported in Luu 2001; placeholder only
  })

  model({
    # ===================================================================
    # 0. Back-transforms (calibrated Table 3 parameters)
    # ===================================================================
    Dr <- exp(lDr)
    Df1 <- exp(lDf1)
    Df2 <- exp(lDf2)
    Dl <- exp(lDl)
    Dgi <- exp(lDgi)
    Dk <- exp(lDk)
    Hair <- exp(lHair)
    kel <- exp(lkel)
    ksi <- exp(lksi)
    kso <- exp(lkso)
    kwi <- exp(lkwi)
    kwo <- exp(lkwo)
    Vweak <- exp(lVweak)
    Vstrong <- exp(lVstrong)
    Qf2Qf1 <- exp(lQf2Qf1)
    Vf2Vf1 <- exp(lVf2Vf1)

    # ===================================================================
    # 1. Fixed physiology, female SD rat (Table 2). Table 2 lists blood
    #    flows in mL/min; the Appendix rate constants are per hour, so
    #    the flows are converted to mL/h (x 60). Qr+Qf+Ql+Qgi+Qk = Qt
    #    exactly, which fixes Qf1 == tabulated fat flow Qf and makes the
    #    flow-partitioned IV forcing sum to the whole dose.
    # ===================================================================
    Qt <- 116.7 * 60 # Table 2 total cardiac output 116.7 mL/min
    Qr <- 44 * 60 # Table 2 richly-perfused flow 44 mL/min
    Qf1 <- 34.2 * 60 # Table 2 fat flow 34.2 mL/min (= weakly-bound-fat flow)
    Ql <- 15.2 * 60 # Table 2 liver flow 15.2 mL/min
    Qgi <- 11.6 * 60 # Table 2 GI flow 11.6 mL/min
    Qk <- 11.7 * 60 # Table 2 kidney flow 11.7 mL/min
    Qf2 <- Qf2Qf1 * Qf1 # deep-fat exchange flow (Table 3 ratio)

    Vlung <- 14 # Table 2 lung volume 14 mL
    Vr <- 146 # Table 2 richly-perfused volume 146 mL
    Vf1 <- 14.2 # Table 2 fat volume 14.2 mL (= weakly-bound-fat volume)
    Vf2 <- Vf2Vf1 * Vf1 # deep-fat volume (Table 3 ratio)
    Vl <- 8.0 # Table 2 liver volume 8.0 mL
    Vgi <- 6.0 # Table 2 GI volume 6.0 mL
    Vk <- 1.9 # Table 2 kidney volume 1.9 mL
    Vblood <- 8.5 # Table 2 plasma volume 8.5 mL

    # ===================================================================
    # 2. Arterial concentration and the flow-partitioned IV forcing
    #    function xi(t) = Qinj*Cinj = DOSE/tinj during the injection.
    # ===================================================================
    Cart <- Hair * Clung
    xi <- (t < tinj) * (DOSE / tinj)

    # ===================================================================
    # 3. Material-balance ODEs (Appendix). Each V*dC/dt is divided by the
    #    compartment volume. Three printed-equation typos are corrected
    #    to satisfy mass balance (noted inline).
    # ===================================================================
    d / dt(Clung) <- (Qt * Cai - Qt * Cart) / Vlung
    d / dt(Cr) <- (Qr * Cart + Qr * xi / Qt - Qr * Cr / Dr) / Vr
    # Deep-fat inflow: Appendix prints Qf2*Cf2/Df2 - Qf2*Cf2/Df2 (self-
    # cancelling); mass balance with the Cf1 equation requires the inflow
    # to be Qf2*Cf1/Df1 (exchange from the weakly-bound fat).
    d / dt(Cf2) <- (Qf2 * Cf1 / Df1 - Qf2 * Cf2 / Df2) / Vf2
    d / dt(Cf1) <-
      (Qf1 * Cart + Qf1 * xi / Qt + Qf2 * Cf2 / Df2 -
        Qf1 * Cf1 / Df1 - Qf2 * Cf1 / Df1) / Vf1
    # Liver venous return: Appendix prints Ql*Cl/Dr; the liver:blood
    # coefficient is Dl, so Dr is a typo corrected to Dl here.
    d / dt(Cl) <- (Ql * Cart + Ql * xi / Qt - Ql * Cl / Dl) / Vl - kel * Cl
    # GI arterial inflow: Appendix prints Ql*Clung*Hair; the GI flow is
    # Qgi, so Ql is a typo corrected to Qgi here.
    d / dt(Cgi) <- (Qgi * Cart + Qgi * xi / Qt - Qgi * Cgi / Dgi) / Vgi
    d / dt(Ck) <- (Qk * Cart + Qk * xi / Qt - Qk * Ck / Dk) / Vk

    d / dt(Cai) <-
      (Qr * Cr / Dr + Qf1 * Cf1 / Df1 + Ql * Cl / Dl + Qgi * Cgi / Dgi +
        Qk * Ck / Dk + kwo * Cwk - kwi * Cai - Qt * Cai) / Vblood
    d / dt(Cwk) <- (kwi * Cai + kso * Cstr - kwo * Cwk - ksi * Cwk) / Vweak
    d / dt(Cstr) <- (ksi * Cwk - kso * Cstr) / Vstrong

    d / dt(Mmet) <- kel * Cl * Vl

    # Weakly-bound blood is preloaded with 5% of the IV dose at t = 0.
    Cwk(0) <- fwb0 * DOSE / Vweak

    # ===================================================================
    # 4. Observations. Measured blood/tissue radioactivity is TOTAL D4
    #    (dissolved + weakly + strongly protein-bound).
    # ===================================================================
    Cc <- (Cai * Vblood + Cwk * Vweak + Cstr * Vstrong) / Vblood
    Cfat <- Cf1
    Cliver <- Cl

    Cc ~ prop(propSd)
  })
}
