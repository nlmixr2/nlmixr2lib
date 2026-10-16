Zhang_2022_dexketoprofen_pbpk <- function() {
  description <- paste(
    "PBPK (compartmental absorption and transit (CAT) model, Berkeley",
    "Madonna 8.3.18). Oral dexketoprofen (as dexketoprofen trometamol",
    "immediate-release tablet) in fasted healthy adults. The gastrointestinal",
    "tract is one stomach segment and seven equal small-intestinal transit",
    "segments, each split into undissolved and dissolved drug (Eqs 1-4,",
    "Figure 2). Undissolved drug dissolves first-order with a pH-specific",
    "rate constant (pH 1.2 in the stomach, pH 4.5 in segment 1, pH 6.8 in",
    "segments 2-7), empties from the stomach at the solid (undissolved) or",
    "liquid (dissolved) gastric emptying rate, and transits the intestine at",
    "a common rate; only dissolved drug in the intestinal segments is",
    "absorbed, at Ka = 2 * Peff / R (Eq 7). Absorbed drug enters a",
    "two-compartment disposition model fitted to intravenous data.",
    "Deterministic: the paper's between-subject variability is a Monte Carlo",
    "draw of the liquid gastric emptying rate from a migrating myoelectric",
    "complex (MMC) model, implemented in the validation vignette as a",
    "per-subject value of lkst_liquid."
  )
  reference <- paste(
    "Zhang X, Ye X, Hu K, Li W, Li W, Xiao Q, Chen L, Yang J (2022).",
    "A Physiologically Based Pharmacokinetic Model for Studying the",
    "Biowaiver Risk of Biopharmaceutics Classification System Class I Drugs",
    "With Rapid Elimination: Dexketoprofen Trometamol Case Study.",
    "Front Pharmacol 13:808456. doi:10.3389/fphar.2022.808456.",
    "MMC gastric-emptying inputs from Oberle RL et al. (1990) The influence",
    "of the interdigestive migrating myoelectric complex on the gastric",
    "emptying of liquids. Gastroenterology 99(5):1275-1282.",
    "doi:10.1016/0016-5085(90)91150-5.",
    sep = " "
  )
  vignette <- "Zhang_2022_dexketoprofen_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # CAT luminal states. Each segment carries an undissolved and a dissolved
  # pool, which no canonical GI compartment resolves; the seven intestinal
  # segments are equal-transit CAT compartments (Yu 1996), not anatomical
  # duodenum / jejunum / ileum.
  paper_specific_compartments <- c(
    "stomach_undissolved",
    "stomach_dissolved",
    paste0("intestine", 1:7, "_undissolved"),
    paste0("intestine", 1:7, "_dissolved")
  )

  compartmentData <- c(
    list(
      stomach_undissolved = list(
        analyte = "dexketoprofen",
        units = "mg",
        specimen = "administration site",
        verified = TRUE
      ),
      stomach_dissolved = list(
        analyte = "dexketoprofen",
        units = "mg",
        specimen = "administration site",
        verified = TRUE
      )
    ),
    stats::setNames(
      rep(list(list(analyte = "dexketoprofen", units = "mg", specimen = "administration site", verified = TRUE)), 14),
      c(paste0("intestine", 1:7, "_undissolved"), paste0("intestine", 1:7, "_dissolved"))
    ),
    list(
      central = list(analyte = "dexketoprofen", units = "mg", specimen = "plasma", verified = TRUE),
      peripheral1 = list(analyte = "dexketoprofen", units = "mg", specimen = "plasma", verified = TRUE)
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_studies = 3L,
    age_range = "adults",
    disease_state = "Healthy volunteers, fasted.",
    dose_range = paste(
      "Single oral 37 mg dexketoprofen trometamol (25 mg dexketoprofen) and",
      "18.5 mg dexketoprofen trometamol (12.5 mg dexketoprofen) tablets;",
      "single intravenous bolus 37 mg dexketoprofen trometamol for the",
      "disposition fit."
    ),
    notes = paste(
      "Not a population fit. Disposition (Vc, K10, K12, K21) was fitted in",
      "Phoenix WinNonlin 7.0 to the mean intravenous profile of Valles 2006",
      "(n = 12, Supplementary Figure S4); oral verification used the mean",
      "profiles of Barbanoj 1998 (25 mg, n = 18, Figure 3) and Mauleon 1996",
      "(12.5 mg). Dissolution rate constants were curve-fitted to the",
      "Garcia-Arieta 2015 USP apparatus 2 (75 rpm, 900 mL) profiles of the",
      "reference tablet (Supplementary Figure S1). Dexketoprofen is a BCS",
      "class I drug; plasma protein binding exceeds 99%."
    )
  )

  ini({
    # Dissolution (Table 1 'Kdis', curve fit to Supplementary Figure S1)
    lkdiss_stomach <- log(8.16); label("Dissolution rate constant in the stomach, pH 1.2 (1/h)") # Table 1 Kdis 8.16 (pH 1.2); Methods: pH 1.2 profile applied to the stomach
    lkdiss_duodenum <- log(6.59); label("Dissolution rate constant in intestinal segment 1, pH 4.5 (1/h)") # Table 1 Kdis 6.59 (pH 4.5); Methods: pH 4.5 profile applied to the first intestinal segment (duodenum)
    lkdiss_intestine <- log(8.29); label("Dissolution rate constant in intestinal segments 2-7, pH 6.8 (1/h)") # Table 1 Kdis 8.29 (pH 6.8); Methods: pH 6.8 profile applied to the other intestinal segments

    # Gastrointestinal transit
    lkst_solid <- fixed(log(1.65)); label("Solid (undissolved drug) gastric emptying rate constant Ktgs (1/h)") # Table 1 Ktgs 1.65 (Wilding 2003)
    lkst_liquid <- fixed(log(5.67)); label("Liquid (dissolved drug) gastric emptying rate constant Ktgl (1/h)") # Table 1 Ktgl 5.67, footnote a: mean of 10,000 runs of the gastric emptying variation model (Eqs 5-6)
    lktr <- fixed(log(2.10)); label("Intestinal transit rate constant between the 7 CAT segments Kt (1/h)") # Table 1 Kt 2.10 (Yu 1996; 7 segments over a 3.3 h small-intestinal transit time)

    # Absorption, Ka = 2 * Peff / R (Eq 7)
    lpeff <- fixed(log(3.37)); label("Human effective jejunal permeability Peff (cm/h)") # Table 1 Peff,man 3.37 (Simcyp 18 prediction)
    lr_intestine <- fixed(log(1.14)); label("Intestinal radius R (cm)") # Eq 7 R, not printed: back-solved by the maintainers from Supplementary Table S3 Cmax 2.90 mg/L (25 mg), cross-checked on Figure 7B segment-1 Fa ~70%

    # Two-compartment disposition (Table 1, fitted in Phoenix WinNonlin 7.0 to i.v. data)
    lvc <- log(3.55); label("Central volume of distribution Vc (L)") # Table 1 Vc 3.55
    lkel <- log(1.68); label("Elimination rate constant from central K10 (1/h)") # Table 1 K10 1.68
    lk12 <- log(1.14); label("Central -> peripheral rate constant K12 (1/h)") # Table 1 K12 1.14
    lk21 <- log(1.24); label("Peripheral -> central rate constant K21 (1/h)") # Table 1 K21 1.24

    # The paper reports no residual-error model; the PBPK model is a
    # forward predictor assessed by %PE of Cmax and AUCtlast (Eq 8).
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported by the paper")
  })

  model({
    kdiss_stomach <- exp(lkdiss_stomach)
    kdiss_duodenum <- exp(lkdiss_duodenum)
    kdiss_intestine <- exp(lkdiss_intestine)
    kst_solid <- exp(lkst_solid)
    kst_liquid <- exp(lkst_liquid)
    ktr <- exp(lktr)
    peff <- exp(lpeff)
    r_intestine <- exp(lr_intestine)
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Eq 7; Peff in cm/h and R in cm give Ka in 1/h
    ka <- 2 * peff / r_intestine

    # Eqs 1-2: stomach; dosed drug enters undissolved and is not absorbed here
    d/dt(stomach_undissolved) <- -kst_solid * stomach_undissolved - kdiss_stomach * stomach_undissolved
    d/dt(stomach_dissolved) <- -kst_liquid * stomach_dissolved + kdiss_stomach * stomach_undissolved

    # Eqs 3-4: intestinal segments. Segment 1 receives gastric emptying and
    # dissolves at pH 4.5; segments 2-7 dissolve at pH 6.8. Drug leaving
    # segment 7 passes to the colon and is not absorbed.
    d/dt(intestine1_undissolved) <- kst_solid * stomach_undissolved - ktr * intestine1_undissolved - kdiss_duodenum * intestine1_undissolved
    d/dt(intestine1_dissolved) <- kst_liquid * stomach_dissolved - ktr * intestine1_dissolved + kdiss_duodenum * intestine1_undissolved - ka * intestine1_dissolved
    d/dt(intestine2_undissolved) <- ktr * intestine1_undissolved - ktr * intestine2_undissolved - kdiss_intestine * intestine2_undissolved
    d/dt(intestine2_dissolved) <- ktr * intestine1_dissolved - ktr * intestine2_dissolved + kdiss_intestine * intestine2_undissolved - ka * intestine2_dissolved
    d/dt(intestine3_undissolved) <- ktr * intestine2_undissolved - ktr * intestine3_undissolved - kdiss_intestine * intestine3_undissolved
    d/dt(intestine3_dissolved) <- ktr * intestine2_dissolved - ktr * intestine3_dissolved + kdiss_intestine * intestine3_undissolved - ka * intestine3_dissolved
    d/dt(intestine4_undissolved) <- ktr * intestine3_undissolved - ktr * intestine4_undissolved - kdiss_intestine * intestine4_undissolved
    d/dt(intestine4_dissolved) <- ktr * intestine3_dissolved - ktr * intestine4_dissolved + kdiss_intestine * intestine4_undissolved - ka * intestine4_dissolved
    d/dt(intestine5_undissolved) <- ktr * intestine4_undissolved - ktr * intestine5_undissolved - kdiss_intestine * intestine5_undissolved
    d/dt(intestine5_dissolved) <- ktr * intestine4_dissolved - ktr * intestine5_dissolved + kdiss_intestine * intestine5_undissolved - ka * intestine5_dissolved
    d/dt(intestine6_undissolved) <- ktr * intestine5_undissolved - ktr * intestine6_undissolved - kdiss_intestine * intestine6_undissolved
    d/dt(intestine6_dissolved) <- ktr * intestine5_dissolved - ktr * intestine6_dissolved + kdiss_intestine * intestine6_undissolved - ka * intestine6_dissolved
    d/dt(intestine7_undissolved) <- ktr * intestine6_undissolved - ktr * intestine7_undissolved - kdiss_intestine * intestine7_undissolved
    d/dt(intestine7_dissolved) <- ktr * intestine6_dissolved - ktr * intestine7_dissolved + kdiss_intestine * intestine7_undissolved - ka * intestine7_dissolved

    absorption <- ka * (intestine1_dissolved + intestine2_dissolved + intestine3_dissolved +
      intestine4_dissolved + intestine5_dissolved + intestine6_dissolved + intestine7_dissolved)

    # Two-compartment disposition (Figure 2)
    d/dt(central) <- absorption - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
