Gasthuys_2023_paracetamol_human <- function() {
  description <- "One-compartment population PK model for oral paracetamol (acetaminophen) in healthy adult men studied fasted, after a reference meal and after infant formula (Gasthuys 2023). Absorption uses a dual input function: a fraction Bio of the dose is released by zero-order input of duration dT1 into a first dosing compartment and the remaining 1 - Bio by zero-order input of duration dT2 into a second dosing compartment that additionally carries a lag time T2; both compartments then deliver drug to the central compartment by a common first-order rate Ka. The second, lagged input reproduces the double absorption peak or shoulder seen in several individual profiles. CL and Vd are apparent (CL/F, Vd/F) because no intravenous arm was run, and both carry body weight as a power function with exponents fixed at 0.75 and 1. Fed state was screened as a covariate but no food effect was retained. NOTE: the published CL/F and Vd/F are internally consistent with each other (they imply a 2.2 h half-life) but are jointly about 2.6-fold too low to reproduce the paper's own non-compartmental AUC for the stated 1000 mg dose; see the validation vignette."
  reference <- "Gasthuys E, Sandra L, Statelova M, Vertzoni M, Vermeulen A. The Use of Population Pharmacokinetics to Extrapolate Food Effects from Human Adults and Beagle Dogs to the Pediatric Population Illustrated with Paracetamol as a Test Case. Pharmaceuticals. 2024;17(1):53. doi:10.3390/ph17010053 (published online 2023-12-28). Human adult study design and data originally reported in Statelova et al. 2019 (cited as ref [16])."
  vignette <- "Gasthuys_2023_paracetamol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Subject-level, time-fixed. Applied as a power function normalised to a 70 kg",
        "reference with exponents fixed at 0.75 for CL/F and 1 for Vd/F",
        "(Gasthuys 2023 Eq. 1 and Table 3, 'Covariate model' rows). Cohort median 81.5 kg,",
        "range 60-104 kg (Gasthuys 2023 Table 4)."
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "paracetamol", units = "mg", specimen = "administration site", verified = TRUE),
    depot2  = list(analyte = "paracetamol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 8L,
    n_studies      = 1L,
    n_observations = "360 plasma paracetamol concentrations from eight adults across three occasions (Gasthuys 2023 Results 2.1). Ten volunteers were enrolled and eight completed (Gasthuys 2023 Table 4 footnote). Concentrations below the 7.5 ng/mL limit of quantification (0.83% of human records) were excluded before fitting.",
    age_range      = "21-48 years (median 25)",
    age_median     = "25 years",
    weight_range   = "60-104 kg",
    weight_median  = "81.5 kg",
    height_range   = "1.69-1.92 m (median 1.85)",
    bmi_range      = "20.3-27.7 kg/m2 (median 23.8)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = "healthy adult male volunteers",
    dose_range     = "single oral 1000 mg paracetamol (42 mL of Panadol suspension at 24 mg/mL) per occasion",
    regions        = "Greece (Red Cross Hospital of Athens; protocol Nr. 4145/14-02-18)",
    notes          = paste(
      "Single-dose, open-label, randomised, crossover, three-period comparative oral",
      "bioavailability study (Gasthuys 2023 Figure 6): fasted (84 mL water plus the",
      "paracetamol suspension given twice over 1 min); reference meal (990 kcal: 240 mL whole",
      "cow's milk, two eggs, two strips of bacon, two slices of toast and 56 g of French",
      "fries) given 30 min before dosing; and infant formula (Noulac, 800 mL total, 520 kcal)",
      "given as 400 mL over 8 min before and 400 mL over 8 min after dosing. Fifteen PK",
      "samples over 0-10 h. A double absorption peak or shoulder was seen in four of eight",
      "subjects on the reference-meal occasion, two of eight fasted and one of eight on infant",
      "formula, which is what the dual input function was added to describe.",
      "The Figure 6 caption states '21 mL paracetamol (Panadol, 168 mg paracetamol)' per",
      "administration, which is inconsistent with Table 4's 42 mL / 24 mg/mL / 1000 mg dose",
      "(21 mL at 24 mg/mL is 504 mg, not 168 mg); 168 mg is the dose used in the companion",
      "beagle-dog study. See the vignette 'Assumptions and deviations' section."
    )
  )

  ini({
    # ======================================================================
    # Structural model -- Gasthuys 2023 Table 3, 'Estimate Human Adults'
    # column. All structural values are estimated point estimates reported
    # with an RSE%; none carries a FIX flag. CL and Vd are apparent
    # (CL/F, Vd/F) per the Table 3 footnote '*', because the human study
    # ran no intravenous arm.
    #
    # Disposition parameters are typical values at the 70 kg reference of
    # Gasthuys 2023 Eq. 1.
    # ======================================================================
    logitfrel <- qlogis(0.45); label("Fraction of the dose entering the first dosing compartment, Bio (unitless)")  # Gasthuys 2023 Table 3: Bio = 0.45, RSE 27.5% -- held on the logit scale so 1 - Bio, the fraction routed to the second dosing compartment, cannot go negative
    lka  <- log(1.79);  label("First-order absorption rate constant Ka from both dosing compartments to central (1/h)")  # Gasthuys 2023 Table 3: ka = 1.79 1/h, RSE 38.1%
    lvc  <- log(27.6);  label("Apparent central volume of distribution Vd/F (L)")                                    # Gasthuys 2023 Table 3: Vd = 27.6 L, RSE 2.42%, footnote * = Vd/F
    lcl  <- log(8.79);  label("Apparent total body clearance CL/F (L/h)")                                            # Gasthuys 2023 Table 3: CL = 8.79 L/h, RSE 1.99%, footnote * = CL/F
    ld1  <- log(0.22);  label("Duration of the zero-order input into the first dosing compartment dT1 (h)")          # Gasthuys 2023 Table 3: dT1 = 0.22 h, RSE 39.2%
    ld2  <- log(2.73);  label("Duration of the zero-order input into the second dosing compartment dT2 (h)")         # Gasthuys 2023 Table 3: dT2 = 2.73 h, RSE 34.2%
    ltlag <- log(0.97); label("Lag time before the second dosing compartment starts releasing, T2 (h)")              # Gasthuys 2023 Table 3: T2 = 0.97 h, RSE 31.8%; Results 2.1 "For the second dosing compartment (1-Bio), a lag time was implemented"

    # ======================================================================
    # Covariate model -- Gasthuys 2023 Table 3, 'Covariate model' rows.
    # Both exponents are reported as [FIX], and Methods 4.2 states the
    # exponent was fixed to 0.75 for clearance and 1 for volumes.
    # ======================================================================
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")  # Gasthuys 2023 Table 3: beta_WTonCL = 0.75 [FIX]
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vd/F (unitless)")  # Gasthuys 2023 Table 3: beta_WTonVd = 1 [FIX]

    # ======================================================================
    # Random effects -- Gasthuys 2023 Table 3.
    #
    # The human column reports between-subject variability on Vd and CL,
    # and inter-occasion variability on the five absorption parameters
    # (Bio, ka, dT1, dT2, T2). The two sets are disjoint. nlmixr2lib has no
    # idiomatic encoding for inter-occasion variability separate from
    # between-subject variability, so per the convention already used in
    # Bienczak_2016_efavirenz.R, Bienczak_2016_nevirapine.R and
    # Svensson_2018_bedaquiline.R, BSV is kept where reported and IOV on a
    # parameter carrying no BSV term is folded in as a BSV-equivalent.
    #
    # Monolix (version 2019R2, per Gasthuys 2023 Methods 4.5) reports omega
    # as the standard deviation of the random effect, so the tabulated
    # values are squared to give the log-scale (or logit-scale) variances
    # nlmixr2 expects.
    # ======================================================================
    etalvc ~ 0.098^2  # Gasthuys 2023 Table 3: BSV Vd = 0.098 (RSE 39.9%); variance = 0.098^2 = 0.009604
    etalcl ~ 0.13^2   # Gasthuys 2023 Table 3: BSV CL = 0.13 (RSE 26.0%); variance = 0.13^2 = 0.0169

    etalogitfrel ~ 0.69^2  # Gasthuys 2023 Table 3: IOV Bio = 0.69 (RSE 20.7%), folded as a BSV-equivalent; variance = 0.69^2 = 0.4761 on the logit scale
    etalka       ~ 0.53^2  # Gasthuys 2023 Table 3: IOV ka  = 0.53 (RSE 56.4%), folded as a BSV-equivalent; variance = 0.53^2 = 0.2809
    etald1       ~ 1.58^2  # Gasthuys 2023 Table 3: IOV dT1 = 1.58 (RSE 25.9%), folded as a BSV-equivalent; variance = 1.58^2 = 2.4964
    etald2       ~ 0.87^2  # Gasthuys 2023 Table 3: IOV dT2 = 0.87 (RSE 22.0%), folded as a BSV-equivalent; variance = 0.87^2 = 0.7569
    etaltlag     ~ 0.59^2  # Gasthuys 2023 Table 3: IOV T2  = 0.59 (RSE 38.3%), folded as a BSV-equivalent; variance = 0.59^2 = 0.3481

    # ======================================================================
    # Residual error -- Gasthuys 2023 Table 3, 'Error model parameters'
    # rows, footnote *1 "combined proportional and additive error model"
    # with a = additive term and b = proportional term. Both happen to be
    # 0.094 in the human fit.
    # ======================================================================
    addSd  <- 0.094; label("Additive residual error on plasma paracetamol (ug/mL)")  # Gasthuys 2023 Table 3: a = 0.094, RSE 23.9%
    propSd <- 0.094; label("Proportional residual error on plasma paracetamol (fraction)")  # Gasthuys 2023 Table 3: b = 0.094, RSE 11.2%
  })

  model({
    # Allometric reference body weight, from Gasthuys 2023 Eq. 1
    # ("theta_70pop: typical population value for the PK parameter in a
    # 70 kg individual").
    wtref <- 70  # kg

    # Dual input function. Bio is carried on the logit scale; keep the
    # fixed effect and eta on their own line so the term stays
    # mu-referenced.
    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    ka   <- exp(lka + etalka)
    d1   <- exp(ld1 + etald1)
    d2   <- exp(ld2 + etald2)
    tlag <- exp(ltlag + etaltlag)

    cl <- exp(lcl + etalcl) * (WT / wtref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wtref)^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition fed by two dosing compartments that
    # share a single first-order absorption rate constant (Gasthuys 2023
    # Results 2.1: "Following the zero-order release into the respective
    # dosing compartments, the first-order absorption (ka) into the
    # central compartment happened").
    d/dt(depot)   <- -ka * depot
    d/dt(depot2)  <- -ka * depot2
    d/dt(central) <-  ka * depot + ka * depot2 - kel * central

    # Each oral administration is encoded by the USER as two parallel dose
    # records in the event table, one targeting cmt = "depot" and one
    # targeting cmt = "depot2", each carrying the full dose amount. The
    # f() multipliers then split it into Bio and 1 - Bio.
    f(depot)     <- frel
    dur(depot)   <- d1
    f(depot2)    <- 1 - frel
    dur(depot2)  <- d2
    alag(depot2) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
