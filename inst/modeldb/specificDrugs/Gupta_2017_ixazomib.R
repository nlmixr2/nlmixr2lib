Gupta_2017_ixazomib <- function() {
  description <- "Three-compartment population pharmacokinetic model for the oral proteasome inhibitor ixazomib (Ninlaro) in 755 adult patients with multiple myeloma, lymphoma, solid tumours, or light-chain amyloidosis pooled across ten phase I, I/II, and III trials including TOURMALINE-MM1 (Gupta 2017). First-order linear absorption with a 13 min lag time describes oral dosing; intravenous and oral data share the same disposition kinetics. Inter-individual variability is estimated on clearance, bioavailability F, and the second peripheral volume V4, with a strong (82%) correlation between log CL and log F. Body surface area on V4 (reference 1.87 m^2, exponent 2.06) is the only retained covariate; sex, age, race, mild/moderate renal impairment, mild hepatic impairment, smoking status, and CYP-modulatory concomitant medications had no clinically relevant effect on systemic exposure. Residual error is additive on log-transformed concentration with a time-after-dose-varying standard deviation declining exponentially from SD1 = 1.90 to SD0 = 0.46 with rate KSD = 0.84/h (Karlsson 1995 model 3)."
  reference <- paste(
    "Gupta N, Diderichsen PM, Hanley MJ, Berg D, van de Velde H, Harvey RD,",
    "Venkatakrishnan K. (2017). Population pharmacokinetic analysis of",
    "ixazomib, an oral proteasome inhibitor, including data from the phase III",
    "TOURMALINE-MM1 study to inform labelling.",
    "Clin Pharmacokinet 56(11):1355-1368.",
    "doi:10.1007/s40262-017-0526-4.",
    sep = " "
  )
  vignette <- "Gupta_2017_ixazomib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power covariate on the second peripheral volume V4 (canonical name `lvp2`) with reference BSA = 1.87 m^2 and exponent 2.06 (Gupta 2017 Table 3). The paper does not state the BSA formula (DuBois / Mosteller / Haycock); recorded as unspecified. The covariate-effect coefficient is encoded as `e_bsa_vp2` consistent with the e_<cov>_<param> convention. BSA was preferred over body weight in model selection because it gave a marginally larger drop in OFV and the V4 inter-individual variance (Gupta 2017 Results section 3.2).",
      source_name        = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 755L,
    n_studies      = 10L,
    age_range      = "23-91 years (median 65)",
    age_median     = "65 years",
    weight_range   = "36.7-151 kg (median 75.5)",
    weight_median  = "75.5 kg",
    sex_female_pct = 42.4,
    race_ethnicity = c(White = 79.9, Black = 5.56, Asian = 11.7, Other = 2.91),
    disease_state  = "Adult patients with advanced hematologic or non-hematologic malignancies: multiple myeloma (n = 632; relapsed/refractory or newly diagnosed), advanced solid tumours (n = 80), lymphoma (n = 28), and relapsed/refractory light-chain amyloidosis (n = 15).",
    dose_range     = "Oral capsule (n = 647) at 0.2-10.6 mg per dose or intravenous (n = 108) at 0.2-8.9 mg per dose. Once-weekly dosing (28-day cycle; days 1, 8, 15; n = 560 patients) or twice-weekly dosing (21-day cycle; days 1, 4, 8, 11; n = 195 patients). Approved label dose is 4 mg orally on days 1, 8, 15 of a 28-day cycle in combination with lenalidomide and dexamethasone.",
    regions        = "Global (TOURMALINE-MM1 phase III; plus seven phase I and two phase I/II studies including East-Asian and Japanese cohorts).",
    studies        = "C16001, C16002, C16003, C16004, C16005, C16007, C16008, C16010/TOURMALINE-MM1, C16013, TB-MC010034 (Gupta 2017 Table 1).",
    baseline_labs  = "Albumin 39 g/L (12-55); AST 22 U/L (4-127); total bilirubin 7 uM (1.71-39.3) = 0.4 mg/dL (0.1-2.3); creatinine clearance 86.8 mL/min (25.8-297, Cockcroft-Gault); hematocrit 0.35 (0.15-0.54); hemoglobin 11.6 g/dL (4.6-16.8). Median values with ranges per Gupta 2017 Table 2.",
    co_medication  = "29.9% single agent; 70.1% ixazomib + lenalidomide + dexamethasone. CYP1A2-modulatory and CYP3A4-modulatory concomitant medications were screened as time-dependent covariates and found not statistically significant.",
    smoking        = "Never 33.1%; Current 4.37%; Former 18.5%; Unknown 44.1% (Gupta 2017 Table 2).",
    notes          = "Demographics from Gupta 2017 Table 2. Patients receiving fixed-dose oral ixazomib (the approved 4 mg once-weekly regimen) made up the majority of the pooled dataset; BSA-based dosing was used in the early phase I studies before the switch to fixed dosing supported by Gupta 2015 (Br J Clin Pharmacol 79:789-800). Severe renal impairment and moderate/severe hepatic impairment were studied separately (label dose 3 mg) and are not represented in this analysis dataset."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters - three-compartment model with first-order
    # oral absorption and an absorption lag time. NONMEM theta names in
    # the paper (Ka, CL, V2, F, Q3, V3, Q4, V4, TLAG) map to nlmixr2
    # canonical names as: V2 -> central volume `lvc`; V3 -> first
    # peripheral volume `lvp`; V4 -> second peripheral volume `lvp2`;
    # Q3 -> inter-compartmental clearance to peripheral1 `lq`;
    # Q4 -> inter-compartmental clearance to peripheral2 `lq2`.
    # Estimates are the log-scale point values from Gupta 2017 Table 3
    # (NONMEM 7.2, FOCEI, the model was developed on log-transformed
    # concentrations). The untransformed values in Table 3 are the
    # exponentiated reference values quoted in the paper text.
    # ------------------------------------------------------------------
    lka  <- log(0.34)  ; label("First-order absorption rate constant (1/h)")           # Gupta 2017 Table 3: log Ka = -1.09 (RSE 8%); Ka = 0.34/h; t1/2 = 124 min
    lcl  <- log(1.86)  ; label("Systemic clearance (L/h)")                              # Gupta 2017 Table 3: log CL = 0.62 (RSE 7%); CL = 1.86 L/h
    lvc  <- log(13.7)  ; label("Central volume of distribution V2 (L)")                 # Gupta 2017 Table 3: log V2 = 2.62 (RSE 4%); V2 = 13.7 L
    lfdepot <- log(0.58) ; label("Absolute oral bioavailability F (fraction)")          # Gupta 2017 Table 3: log F = -0.55 (RSE 9%); F = 58%
    lq   <- log(5.18)  ; label("Inter-compartmental clearance Q3 (L/h)")                # Gupta 2017 Table 3: log Q3 = 1.65 (RSE 7%); Q3 = 5.18 L/h
    lvp  <- log(309)   ; label("Volume of first peripheral compartment V3 (L)")         # Gupta 2017 Table 3: log V3 = 5.73 (RSE 1%); V3 = 309 L
    lq2  <- log(26.1)  ; label("Inter-compartmental clearance Q4 (L/h)")                # Gupta 2017 Table 3: log Q4 = 3.26 (RSE 2%); Q4 = 26.1 L/h
    lvp2 <- log(205)   ; label("Volume of second peripheral compartment V4 (L)")        # Gupta 2017 Table 3: log V4 = 5.32 (RSE 1%); V4 = 205 L
    ltlag <- log(13 / 60) ; label("Absorption lag time (h)")                            # Gupta 2017 Table 3: log TLAG = -1.52 (RSE 0%); TLAG = 13 min = 0.2167 h

    # ------------------------------------------------------------------
    # Covariate effect: power-form BSA effect on V4 with reference BSA
    # = 1.87 m^2 (Gupta 2017 Table 3 V4[BSA] row). At the 5th and 95th
    # percentile of BSA (1.5 and 2.25 m^2) the implied factor on V4 is
    # -37% and +46% respectively, both stated in the same Table 3 row.
    # ------------------------------------------------------------------
    e_bsa_vp2 <- 2.06  ; label("Power exponent for BSA effect on V4 (unitless)")        # Gupta 2017 Table 3 V4[BSA] = 2.06 (RSE 18%)

    # ------------------------------------------------------------------
    # Inter-individual variability. The paper reports IIV as % CV (Table
    # 3 random-effect block); the internal log-scale variance is
    # omega^2 = log(1 + CV^2). IIV on log CL and log F is correlated
    # (rho = 0.82) and entered as a 2x2 block; IIV on log V4 is on its
    # own. No IIV on log Ka, log V2, log Q3, log V3, log Q4, or log
    # TLAG (Gupta 2017 Table 3 lists no IIV row for these structural
    # parameters).
    # Block lower-triangle order: var_lcl, cov(lcl, lfdepot), var_lfdepot.
    # cov = rho * sqrt(var_lcl * var_lfdepot).
    # ------------------------------------------------------------------
    etalcl + etalfdepot ~ c(0.17697, 0.22550, 0.42726)                                  # Gupta 2017 Table 3: IIV CL 44% -> var = log(1 + 0.44^2) = 0.17697; IIV F 73% -> var = log(1 + 0.73^2) = 0.42726; rho(CL,F) = 0.82 -> cov = 0.82 * sqrt(0.17697 * 0.42726) = 0.22550
    etalvp2 ~ 0.48492                                                                   # Gupta 2017 Table 3: IIV V4 79% -> var = log(1 + 0.79^2) = 0.48492

    # ------------------------------------------------------------------
    # Residual error: additive on log-transformed concentration
    # (Eq. 2 of Gupta 2017), which maps to a proportional error in
    # linear concentration space for nlmixr2. The residual standard
    # deviation is itself a function of time after dose (Eq. 3,
    # Karlsson 1995 model 3) declining exponentially from SD1 at
    # t = 0 toward the plateau SD0:
    #   sd(eps) = SD0 + (SD1 - SD0) * exp(-KSD * tad)
    # The decay half-life log(2)/KSD = log(2)/0.84 = 0.825 h (~50 min)
    # is quoted in Table 3 as "t1/2 = 50 min". The derived residual
    # SD (`errSd`) is computed in `model()` and passed to `prop()`.
    # All three components are estimated.
    # ------------------------------------------------------------------
    sd0  <- 0.46 ; label("Steady-state residual SD on log-concentration (Karlsson 1995 SD0)")  # Gupta 2017 Table 3 SD0 = 0.46 (RSE 3%)
    sd1  <- 1.90 ; label("Initial residual SD on log-concentration at t=0 post-dose (Karlsson 1995 SD1)") # Gupta 2017 Table 3 SD1 = 1.90 (RSE 11%)
    ksd  <- 0.84 ; label("Rate constant of residual-SD decline post-dose (1/h)")                # Gupta 2017 Table 3 KSD = 0.84 (RSE 22%); residual SD t1/2 = log(2)/0.84 = 50 min
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters. BSA scales V4 by (BSA / 1.87)^2.06
    # per Gupta 2017 Table 3. No covariate effects on any other
    # parameter were retained in the final model (Gupta 2017 Results
    # section 3.2 and Discussion). IIV is applied as a multiplicative
    # log-normal random effect on CL, F (depot bioavailability), and
    # V4 (`lvp2`).
    # ------------------------------------------------------------------
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    q2   <- exp(lq2)
    vp2  <- exp(lvp2 + etalvp2) * (BSA / 1.87) ^ e_bsa_vp2
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)

    # ------------------------------------------------------------------
    # 2. ODE system. depot -> central first-order absorption (rate ka);
    # central <-> peripheral1 (Q3) and central <-> peripheral2 (Q4)
    # distribution; first-order linear elimination (CL) from central.
    # IV dose lands directly in `central` via cmt column; oral dose
    # lands in `depot` and transits via ka after the lag time.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          (q  / vc) * central + (q  / vp ) * peripheral1 -
                          (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <-  (q  / vc) * central - (q  / vp ) * peripheral1
    d/dt(peripheral2) <-  (q2 / vc) * central - (q2 / vp2) * peripheral2

    # ------------------------------------------------------------------
    # 3. Bioavailability and absorption lag on the oral depot
    # compartment (Gupta 2017 Fig. 2 and Methods section 2.2).
    # ------------------------------------------------------------------
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 4. Concentration in ng/mL. Dose units are mg, vc is in L, so
    # central/vc has units mg/L = ug/mL. Multiply by 1000 to express
    # the observation in ng/mL, matching the assay reporting units
    # (Gupta 2017 Methods section 2.1 LC-MS/MS assay).
    # ------------------------------------------------------------------
    Cc <- (central / vc) * 1000

    # ------------------------------------------------------------------
    # 5. Time-varying residual SD on log-concentration. `tad()` returns
    # the time since the most recent dose; the source paper says the
    # residual SD "decreases over the first 2-4 h post-dose towards a
    # plateau" so the time-after-dose interpretation is consistent
    # with the paper's text. NONMEM additive-on-log-scale residual
    # error maps to nlmixr2 `prop()` in linear concentration space
    # (per references/naming-conventions.md NONMEM -> nlmixr2 syntax
    # table).
    # ------------------------------------------------------------------
    errSd <- sd0 + (sd1 - sd0) * exp(-ksd * tad())
    Cc ~ prop(errSd)
  })
}
