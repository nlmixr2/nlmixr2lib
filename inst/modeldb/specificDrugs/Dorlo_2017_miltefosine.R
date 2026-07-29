Dorlo_2017_miltefosine <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption",
    "for miltefosine in 95 Eastern African adults and children (>=7 years)",
    "with visceral leishmaniasis (Dorlo 2017), enrolled across three",
    "treatment centres in Kenya and Sudan and randomised to either a",
    "28-day monotherapy regimen of oral miltefosine 2.5 mg/kg/day or a",
    "10-day oral miltefosine 2.5 mg/kg/day arm combined with a single",
    "10 mg/kg liposomal amphotericin B IV dose on day 1. CL/F, Q/F, Vc/F,",
    "and Vp/F are allometrically scaled on fat-free mass (exponents 0.75",
    "and 1.0; reference FFM 53 kg). Relative bioavailability is",
    "structurally fixed at 100% from the end of the initial reduced",
    "absorption window onwards, and reduced by a typical 74.3% during",
    "the window itself (0 < t <= 7 days for monotherapy, 0 < t <= 1 day",
    "for the combination arm); the duration is regimen-dependent via the",
    "MIL_REGIMEN indicator. The combined-error residual model is",
    "proportional (31.0%) with additive component fixed at 0.001 ug/mL."
  )
  reference <- paste(
    "Dorlo TPC, Kip AE, Younis BM, Ellis SJ, Alves F, Beijnen JH,",
    "Njenga S, Kirigi G, Hailu A, Olobo J, Musa AM, Balasegaram M,",
    "Wasunna M, Karlsson MO, Khalil EAG. Visceral leishmaniasis relapse",
    "hazard is linked to reduced miltefosine exposure in patients from",
    "Eastern Africa: a population pharmacokinetic/pharmacodynamic study.",
    "J Antimicrob Chemother. 2017;72(11):3131-3140.",
    "doi:10.1093/jac/dkx283. ClinicalTrials.gov NCT01067443.",
    sep = " "
  )
  vignette <- "Dorlo_2017_miltefosine"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    FFM = list(
      description        = "Baseline fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling on CL, Vc, Q, Vp with reference 53 kg",
        "(Dorlo 2017 Table 2 footnote b: 'Values reported are normalized",
        "to a standard fat-free mass of 53 kg'). Cohort median FFM was",
        "29.9 kg (range 15.2-55.9, Table 1). Computed per subject from",
        "body weight, height, and sex via the Janmahasatian 2005",
        "semi-mechanistic FFM formula (Dorlo 2017 cites Anderson and",
        "Holford [ref 19] for the allometric form)."
      ),
      source_name        = "Fat-free mass (kg)"
    ),
    MIL_REGIMEN = list(
      description        = paste(
        "Miltefosine dosing-regimen indicator. 1 = monotherapy arm (28",
        "days oral miltefosine 2.5 mg/kg/day, max 150 mg/day); 0 =",
        "combination arm (single IV liposomal amphotericin B 10 mg/kg",
        "on day 1 plus 10 days oral miltefosine 2.5 mg/kg/day, max 150",
        "mg/day). Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (combination arm: single LAmB IV + 10-day miltefosine)",
      notes              = paste(
        "Selects the duration of the initial reduced-bioavailability",
        "window inside model() via tlowf = 7 * MIL_REGIMEN + 1 * (1 -",
        "MIL_REGIMEN): 7 days for monotherapy, 1 day for the combination",
        "arm. Source: Dorlo 2017 Results 'Pharmacokinetics of",
        "miltefosine' paragraph 2 and Table 2 footnote d ('Duration of",
        "the reduction in F was empirically determined to be from",
        ">0 to <=7 days for the monotherapy regimen and from >0 to <=1",
        "day for the combination therapy regimen'). The faster recovery",
        "in the combination arm is hypothesised to reflect a quicker",
        "improvement in patients' physiological state after the LAmB",
        "infusion and possibly a transient miltefosine-LAmB membrane",
        "interaction (Dorlo 2017 Discussion)."
      ),
      source_name        = "Treatment arm (monotherapy MIL vs combination LAmB + MIL)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 95L,
    n_studies      = 1L,
    age_range      = "7-41 years",
    age_median     = "13 years",
    weight_range   = "15-65 kg",
    weight_median  = "32 kg",
    ffm_range      = "15.2-55.9 kg",
    ffm_median     = "29.9 kg",
    height_range   = "107-185 cm",
    sex_female_pct = 13.7,
    race_ethnicity = "Eastern African (Kenyan and Sudanese cohorts; no further breakdown reported)",
    disease_state  = paste(
      "Eastern African adults and children with parasitologically",
      "confirmed visceral leishmaniasis (kala-azar) starting first-line",
      "treatment; 41 of 95 (43.2%) paediatric (>=7 and <12 years), 49",
      "of 95 (51%) malnourished by BMI / BMI-for-age z-score criteria,",
      "15 of 95 (15.7%) relapsed during the 210-day follow-up."
    ),
    dose_range     = paste(
      "Oral miltefosine (Impavido) 2.5 mg/kg/day, maximum 150 mg/day,",
      "administered as 50 mg gelatin capsules with food and directly",
      "observed: 28 days in the monotherapy arm (median total 72.6",
      "mg/kg) and 10 days in the combination arm (median total 20.8",
      "mg/kg). The combination arm also received a single IV dose of",
      "liposomal amphotericin B (AmBisome) 10 mg/kg on day 1 (not",
      "included in this miltefosine PK model)."
    ),
    regions        = "Eastern Africa (Kenya: Kimalel; Sudan: Dooka, Kassab)",
    co_medication  = paste(
      "Combination arm received single IV liposomal amphotericin B 10",
      "mg/kg on day 1. The PK model describes miltefosine concentrations",
      "only; LAmB is captured in the published time-to-event PD model",
      "as a proportional reduction in baseline relapse hazard (not",
      "encoded here)."
    ),
    samples        = paste(
      "608 post-dose miltefosine plasma concentrations (excluding",
      "pre-dose samples) from 95 patients, of which 28 (4.6%) were",
      "below the LLOQ of 4.00 ng/mL and handled via the M3 likelihood",
      "method. Sampling schedule (combination arm, >=12 years): 4 h, 8",
      "h, day 2, day 5, day 9 post first dose plus days 60 and 210;",
      "monotherapy arm (>=12 years): 4 h, 8 h, days 2, 6, 13, 20, 27",
      "plus days 60 and 210. Children <12 years had a sparser schedule."
    ),
    notes          = paste(
      "Demographics from Dorlo 2017 Table 1 (median, range) pooled",
      "across the monotherapy (n=48) and combination (n=47) arms.",
      "Eastern African patients had lower miltefosine exposure (~30%",
      "lower AUC0-28d) than historical Indian and European cohorts at",
      "the same 2.5 mg/kg/day regimen, driven by the transient",
      "first-week bioavailability reduction modelled here. Children",
      "<12 years had ~30-40% lower AUC0-end and Time>EC90 than adults",
      "within each arm, reinforcing the case for paediatric dose",
      "adaptation (Dorlo 2017 Table 3, Discussion). NONMEM 7.3,",
      "FOCE-I, M3 method for BLQ; bootstrap n=1000."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Dorlo 2017 Table 2 'Pharmacokinetics'
    # block. Apparent values reported with F structurally fixed at 1
    # (absolute miltefosine bioavailability unknown, Results paragraph
    # 'Pharmacokinetics of miltefosine'). CL/F, Q/F, Vc/F, and Vp/F are
    # normalized to a 53 kg reference fat-free mass (Table 2 footnote
    # b); allometric exponents fixed at 0.75 on clearances and 1.0 on
    # volumes per Anderson and Holford theory (Table 2 footnote b and
    # Methods 'Population pharmacokinetic analysis' citing [19]).
    # The Discussion paragraph mentions a typical absorption rate of
    # 1.25 day^-1; the authoritative final-model estimate is Table 2's
    # 1.49 day^-1 and is used here (see Assumptions and deviations).
    # ============================================================
    lka <- log(1.49)
    label("First-order oral absorption rate constant ka (1/day)")          # Table 2 Pharmacokinetics row 'Absorption rate, k_a (/day)': 1.49 (RSE 17.7%)
    lcl <- log(4.29)
    label("Apparent oral clearance CL/F at FFM = 53 kg (L/day)")            # Table 2 Pharmacokinetics row 'Clearance, CL/F (L/day)': 4.29 (RSE 3.22%)
    lvc <- log(51.7)
    label("Apparent central volume Vc/F at FFM = 53 kg (L)")                # Table 2 Pharmacokinetics row 'Central volume of distribution, V_c/F (L)': 51.7 (RSE 4.33%)
    lq  <- log(0.0266)
    label("Apparent inter-compartmental clearance Q/F at FFM = 53 kg (L/day)") # Table 2 Pharmacokinetics row 'Intercompartmental clearance, Q/F (L/day)': 0.0266 (RSE 40.7%)
    lvp <- log(2.25)
    label("Apparent peripheral volume Vp/F at FFM = 53 kg (L)")              # Table 2 Pharmacokinetics row 'Peripheral volume of distribution, V_p/F (L)': 2.25 (RSE 14.1%)

    # Bioavailability anchor (F = 1 structurally at the end of treatment)
    # and the transient reduction during the initial absorption window.
    # Re-parameterised on the log scale as a multiplier <1 so the
    # individual-level bioavailability remains positive across the
    # log-normal BSV draw (see model() for the predicate that gates
    # this multiplier to t > 0 and t <= tlowf).
    lfdepot <- fixed(log(1))
    label("Reference bioavailability F at end of treatment (unitless, FIXED at 1)") # Table 2 'F (%, at end of treatment): 100 fixed'; Methods 'Population pharmacokinetic analysis'
    lfred_mult <- log(1 - 0.743)
    label("Log of the F multiplier during the reduced-bioavailability window (unitless)") # Table 2 'Reduction in F at baseline (% change from end of treatment): -74.3' -> multiplier 1 - 0.743 = 0.257 of normal F

    # Allometric exponents fixed at the canonical Anderson and Holford
    # values; reference FFM 53 kg matches Table 2 footnote b.
    e_ffm_cl_q  <- fixed(0.75)
    label("Allometric exponent shared by CL/F and Q/F (unitless, FIXED)")  # Table 2 footnote b: 'allometrically scaled (power exponents of 0.75 and 1, respectively) based on fat-free mass'
    e_ffm_vc_vp <- fixed(1.0)
    label("Allometric exponent shared by Vc/F and Vp/F (unitless, FIXED)") # Table 2 footnote b: see above; volumes get exponent 1

    # ============================================================
    # Inter-individual variability -- Dorlo 2017 Table 2 'Between-
    # subject variability (% RSE)' column. CV% reported in the paper
    # is converted to log-normal omega^2 = log(1 + CV^2). BSV could
    # only be identified for ka, CL/F, and the temporary reduction in
    # F (Results paragraph 'Pharmacokinetics of miltefosine'); the
    # 'NE' entries for Vc/F, Q/F, Vp/F, and F are not encoded as
    # estimated etas.
    # ============================================================
    etalka       ~ log(1 + 0.680^2)
    # Table 2 BSV row 'Absorption rate, k_a (/day)': 68.0% CV (RSE 24.5%) -> omega^2 = log(1 + 0.680^2)
    etalcl       ~ log(1 + 0.170^2)
    # Table 2 BSV row 'Clearance, CL/F (L/day)': 17.0% CV (RSE 23.3%)  -> omega^2 = log(1 + 0.170^2)
    etalfred_mult ~ log(1 + 0.964^2)
    # Table 2 BSV row 'Reduction in F at baseline': 96.4% CV (RSE 19.5%) -> omega^2 = log(1 + 0.964^2). The variability is applied multiplicatively on the log-scale F multiplier so F stays positive on every draw.

    # ============================================================
    # Residual error -- Dorlo 2017 Table 2 'Pharmacokinetics' block:
    # combined proportional (31.0% CV) + additive (0.001 ug/mL, FIXED).
    # ============================================================
    propSd <- 0.310
    label("Proportional residual error on Cc (fraction)")                  # Table 2 'Residual proportional error (%)': 31.0 (RSE 5.74%)
    addSd  <- fixed(0.001)
    label("Additive residual error on Cc (ug/mL, FIXED)")                  # Table 2 'Residual additive error (ug/mL)': 0.001 fixed
  })

  model({
    # ------------------------------------------------------------
    # Regimen-dependent duration of the initial reduced-F window.
    # MIL_REGIMEN = 1 -> monotherapy 28-day arm: tlowf = 7 days
    # MIL_REGIMEN = 0 -> combination 10-day + LAmB arm: tlowf = 1 day
    # (Dorlo 2017 Table 2 footnote d, Results 'Pharmacokinetics of
    # miltefosine' paragraph 2).
    # ------------------------------------------------------------
    tlowf <- 7 * MIL_REGIMEN + 1 * (1 - MIL_REGIMEN)

    # ------------------------------------------------------------
    # Individual PK parameters with allometric scaling on fat-free
    # mass at the 53 kg reference (Dorlo 2017 Table 2 footnote b).
    # Fixed allometric exponents (0.75 on CL/Q, 1 on Vc/Vp) per
    # Anderson and Holford theory.
    # ------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (FFM / 53)^e_ffm_cl_q
    vc <- exp(lvc)          * (FFM / 53)^e_ffm_vc_vp
    q  <- exp(lq)           * (FFM / 53)^e_ffm_cl_q
    vp <- exp(lvp)          * (FFM / 53)^e_ffm_vc_vp

    # ------------------------------------------------------------
    # Bioavailability. F = F_typical at the end of treatment
    # (lfdepot, FIXED at log(1)); during the initial window
    # 0 < t <= tlowf, F is reduced by an individual multiplier
    # exp(lfred_mult + etalfred_mult) with typical value
    # exp(log(0.257)) = 0.257 (Table 2 'Reduction in F at baseline').
    # The (t > 0) * (t <= tlowf) product implements the AND predicate
    # so the reduced-F factor is applied only inside the window.
    # ------------------------------------------------------------
    fred_active <- (t > 0) * (t <= tlowf)
    fred_mult_i <- exp(lfred_mult + etalfred_mult)
    fdepot      <- exp(lfdepot) * (fred_mult_i * fred_active + 1 * (1 - fred_active))

    # ------------------------------------------------------------
    # Two-compartment oral PK (depot -> central <-> peripheral1)
    # with first-order absorption and first-order elimination from
    # the central compartment (Results paragraph 'Pharmacokinetics
    # of miltefosine'; Dorlo 2008 [ref 11] base structural model).
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ------------------------------------------------------------
    # Observation. Dose units mg; Vc units L -> central / vc in
    # mg/L = ug/mL (the paper's reporting unit). Combined error.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
