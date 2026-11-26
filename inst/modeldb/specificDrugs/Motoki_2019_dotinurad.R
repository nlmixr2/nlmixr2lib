Motoki_2019_dotinurad <- function() {
  description <- paste(
    "PK/PD model of the URAT1-selective uricosuric agent dotinurad and serum",
    "urate (SUA) in healthy Japanese men (Motoki 2019). Dotinurad PK: one",
    "compartment with first-order absorption and first-order elimination",
    "parameterised as V/F and Ke. Urate PD: a single-pool urate mass balance",
    "with a constant production rate (Synthesis = 3/2 x the mean pre-dose",
    "24-h urinary urate excretion), a constant non-renal (gut) clearance equal",
    "to half the mean pre-dose renal clearance, and renal clearance written as",
    "filtration minus fractional tubular reabsorption, CLRUA = GFR - freabs *",
    "GFR * (1 - Emax * Cp / (EC50 + Cp)), so dotinurad inhibits reabsorption",
    "with a simple Emax function of plasma concentration. GFR is creatinine",
    "clearance normalised to 1.73 m^2; it is the only subject-level input",
    "besides the pre-dose SUA, which sets the initial urate pool.",
    "Typical-value model: the paper reports no between-subject or residual",
    "variability."
  )
  reference <- paste(
    "Motoki K, Igarashi T, Omura K, Nakatani H, Iwanaga T, Tamai I, Ohashi T.",
    "Pharmacokinetic/pharmacodynamic modeling and simulation of dotinurad, a",
    "novel uricosuric agent, in healthy volunteers.",
    "Pharmacol Res Perspect. 2019;7(6):e00533.",
    "doi:10.1002/prp2.533"
  )
  vignette <- "Motoki_2019_dotinurad"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "dotinurad", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dotinurad", units = "mg", specimen = "plasma", verified = TRUE),
    urate = list(analyte = "urate", units = "mg", specimen = "serum", verified = TRUE),
    urate_urine = list(analyte = "urate", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance normalised to a body surface area of 1.73 m^2, used as the glomerular filtration rate (GFR) of urate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Motoki 2019 Methods (PK/PD modeling): 'creatinine clearance (Ccr), adjusted to body surface area of 1.73 m2, was used as GFR'. Enters Eq. 8 linearly as the filtration term of renal urate clearance, CLRUA = CRCL - freabs * CRCL * (1 - Emax * Cp / (EC50 + Cp)), with no reference value. Because urate production and non-renal clearance are constants, a lower CRCL lowers renal urate clearance and raises the urate steady state, which is the use for renal impairment proposed in the Motoki 2019 Discussion. The cohort mean is not printed; the pre-dose identity CLRUA = CRCL * (1 - freabs) (Eq. 9) with the printed mean baseline CLRUA 7.10 mL/min and freabs 0.943 implies 7.10 / 0.057 = 124.6 mL/min/1.73 m^2, which the vignette uses for the typical subject.",
      source_name = "Ccr (GFR)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 72L,
    n_studies = 2L,
    age_range = "20-35 years",
    weight_range = "Body mass index 18.5 to < 25.0 kg/m^2; mean body weight about 60 kg in the single-dose study (Motoki 2019 Discussion).",
    sex_female_pct = 0,
    race_ethnicity = "Japanese",
    disease_state = "Healthy men (normouricaemic). Mean (SD) pre-dose SUA 5.71 (0.75) mg/dL, 24-h urinary urate excretion 575 (100) mg/day and renal urate clearance 7.10 (1.51) mL/min across the single-dose study.",
    dose_range = "Single oral doses of 0.5, 1, 2, 5, 10 and 20 mg fasted (plus 5 mg fed); 2 or 5 mg once daily for 7 days fed.",
    regions = "Japan",
    notes = "Model building used the single ascending dose study (54 subjects, 6 active + 3 placebo per dose level; NCT02348307). The multiple ascending dose study (18 subjects, 6 active + 3 placebo per dose; NCT02348333) was held out for external validation of the SUA prediction (Motoki 2019 Figure 4). PK parameters were estimated with Phoenix WinNonlin; VdUA was the mean of individual calculations in 35 dotinurad-treated subjects (Eq. 4), freabs the mean in 18 placebo subjects (Eq. 9), and Emax / EC50 a WinNonlin simple Emax fit of interval-wise effect on CLRUA versus mean Cp."
  )

  ini({
    # Dotinurad PK (Motoki 2019 Table 3; fasted single-dose data)
    lka <- log(0.770); label("Dotinurad absorption rate constant Ka (1/h)") # Table 3: Ka 0.770 /h (14.2% CV)
    lvc <- log(10.585); label("Dotinurad apparent volume of distribution V/F (L)") # Table 3: V/F 10 585 mL (7.2% CV)
    lkel <- log(0.0795); label("Dotinurad elimination rate constant Ke (1/h)") # Table 3: Ke 0.0795 /h (9.3% CV)

    # Urate system (Motoki 2019 Table 3)
    lvd_urate <- log(16.0); label("Distribution volume of urate VdUA (L)") # Table 3: VdUA 16.0 L (36.6% CV), Eq. 4
    f_reab_urate <- 0.943; label("Fractional tubular reabsorption of filtered urate freabs (fraction)") # Table 3: freabs 0.943 (0.6% CV), Eq. 9

    # Urate system constants derived from the single-dose study's mean
    # pre-dose values (Results: SUA 5.71 mg/dL, XuUA 575 mg/day, CLRUA 7.10 mL/min)
    lksyn <- fixed(log(35.9375)); label("Urate production rate Synthesis (mg/h)") # Eq. 7: Synthesis = 3/2 x XuUA = 1.5 x 575 mg/day = 862.5 mg/day = 35.9375 mg/h
    lcl_nonren_urate <- fixed(log(0.213)); label("Non-renal (gut) clearance of urate CLGutUA (L/h)") # Eq. 6: CLGutUA = 1/2 x CLRUA = 0.5 x 7.10 mL/min = 3.55 mL/min = 0.213 L/h
    lrbase <- fixed(log(5.71)); label("Pre-dose serum urate, initial value of the urate pool (mg/dL)") # Results: mean baseline SUA 5.71 (0.75) mg/dL

    # Dotinurad effect on renal urate reabsorption (Motoki 2019 Table 3, Eq. 8)
    lemax <- log(0.51); label("Maximum fractional inhibition of urate reabsorption Emax (fraction)") # Table 3: Emax 0.51 (4.6% CV)
    lec50 <- log(196); label("Dotinurad plasma concentration at half-maximal effect EC50 (ng/mL)") # Table 3: EC50 196 ng/mL (11.9% CV)
  })

  model({
    # Dotinurad PK
    ka <- exp(lka)
    vc <- exp(lvc)
    kel <- exp(lkel)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # central (mg) / vc (L) = mg/L; x 1000 gives ng/mL, the unit of EC50
    Cc <- 1000 * central / vc

    # Urate system
    vd_urate <- exp(lvd_urate)
    ksyn <- exp(lksyn)
    cl_nonren_urate <- exp(lcl_nonren_urate)
    rbase <- exp(lrbase)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    # CRCL (mL/min/1.73 m^2) to L/h
    gfr <- CRCL * 0.06

    # Eq. 8: renal urate clearance = filtration - inhibitable reabsorption.
    # Dotinurad affects only renal clearance (Methods assumption 3).
    effect <- emax * Cc / (ec50 + Cc)
    cl_renal_urate <- gfr - f_reab_urate * gfr * (1 - effect)

    # Eq. 1 in differential form; urate in mg, vd_urate in L, so sua_mgl is mg/L
    sua_mgl <- urate / vd_urate
    d/dt(urate) <- ksyn - (cl_renal_urate + cl_nonren_urate) * sua_mgl
    d/dt(urate_urine) <- cl_renal_urate * sua_mgl
    # Initial pool from the observed pre-dose SUA (Model simulation section);
    # rbase is mg/dL, x 10 gives mg/L
    urate(0) <- rbase * 10 * vd_urate

    # Serum urate in mg/dL
    sUA <- sua_mgl / 10
  })
}
