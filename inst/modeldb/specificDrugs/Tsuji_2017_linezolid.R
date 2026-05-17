Tsuji_2017_linezolid <- function() {
  description <- paste(
    "Population PK/PD model for linezolid in hospitalized adult and",
    "pediatric patients with MRSA or gram-positive cocci infections",
    "(Tsuji 2017). PK is a two-compartment model with first-order",
    "oral absorption and an additive renal-plus-non-renal clearance",
    "structure (CL = CL_nonren + CL_renal * RF, where RF = CrCl /",
    "100 mL/min/70 kg standardized to 70 kg by allometry); plasma",
    "total and unbound concentrations are modelled simultaneously",
    "with an estimated fraction-unbound (FU = 0.823) linking the two.",
    "PD is a Friberg-style semi-mechanistic platelet turnover model",
    "(one proliferating compartment, three transit compartments, one",
    "circulating compartment) with an empirical (PLTZERO/PLT)^gamma",
    "feedback term and a published mixture model of two thrombocytopenia",
    "mechanisms: linear inhibition of platelet synthesis (PDI, 97% of",
    "patients, SLOPE on RFORM) and saturable stimulation of platelet",
    "elimination (PDS, 3% of patients, Emax on Kcirc), selected per",
    "subject by the binary covariate MIX_PDI."
  )
  reference <- paste(
    "Tsuji Y, Holford NHG, Kasai H, Ogami C, Heo Y-A, Higashi Y,",
    "Mizoguchi A, To H, Yamamoto Y. (2017).",
    "Population pharmacokinetics and pharmacodynamics of",
    "linezolid-induced thrombocytopenia in hospitalized patients.",
    "Br J Clin Pharmacol 83(8):1758-1772.",
    "doi:10.1111/bcp.13262"
  )
  vignette <- "Tsuji_2017_linezolid"
  units <- list(
    time = "hour", dosing = "mg",
    concentration = "mg/L", platelet = "cells/uL"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL/Q (exponent 0.75) and VC/VP (exponent 1)",
        "with reference 70 kg per Tsuji 2017 Methods Equation 3. Also used",
        "to standardise the Cockcroft-Gault CrCl to 70 kg before the renal",
        "function ratio RF = (CrCl x (70/WT)^0.75) / 100."
      ),
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear fractional-change effect on non-renal CL centred at the",
        "cohort median 69 years per Tsuji 2017 Methods Equation 4:",
        "FAGE_CL = 1 + KAGECL * (AGE - 69), with KAGECL = -0.021/year",
        "(roughly -2% CL per year above 69). FAGE on Q, VC and VP are",
        "held at 1."
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance by the Cockcroft-Gault formula (raw mL/min,",
        "NOT BSA-normalized)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried as the raw Cockcroft-Gault value in mL/min, not the canonical",
        "BSA-normalized mL/min/1.73 m^2 form (the same deviation documented",
        "in Jonckheere_2019_cefepime.R and Delattre_2010_amikacin.R). Inside",
        "the model the value is first standardised to a 70 kg body weight via",
        "Holford-style allometry CrCl_70 = CrCl * (70/WT)^0.75 and then",
        "expressed as a ratio RF = CrCl_70 / 100 mL/min/70 kg (Tsuji 2017",
        "Methods Equation 5). Pediatric patients aged 1, 5, 8 and 13 years",
        "were assigned RF = 0.5 by the authors (Methods, paragraph after",
        "Equation 5). To reproduce that handling at simulation time, supply",
        "a CRCL value that yields RF_paper = 0.5 (i.e., CRCL = 50 * (WT/70)^0.75",
        "for those four pediatric IDs)."
      ),
      source_name        = "CLcr"
    ),
    MIX_PDI = list(
      description        = paste(
        "Latent mixture-model indicator: 1 = subject classified to the",
        "PDI (inhibition of platelet synthesis) sub-population, 0 = subject",
        "classified to the PDS (stimulation of platelet elimination)",
        "sub-population."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "PDS = 0",
      notes              = paste(
        "Not a measured patient covariate; this is the per-subject mixture",
        "assignment of the published model. Population probability of MIX_PDI = 1",
        "is the estimated mixture fraction FPOP_inhibit = 0.969 (Tsuji 2017",
        "Table 2; 95% CI 0.867-1.00). For typical-value simulation, set",
        "MIX_PDI = 1 to reproduce the dominant PDI mechanism (78/80 patients",
        "in the source dataset; Figure 5 left panel) or MIX_PDI = 0 to",
        "reproduce the rarer PDS mechanism (2/80 patients; Figure 5 right",
        "panel). For population simulation, draw MIX_PDI ~ Bernoulli(0.969)",
        "per subject. The reference category (MIX_PDI = 0 = PDS) is chosen so",
        "the binary numerically maps onto the paper's mixture indicator;",
        "the dominant class is the non-reference category by intent."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 81L,
    n_studies      = 2L,
    age_range      = "1-85 years (2.5-97.5% interval; median 69)",
    weight_range   = "21.0-99.5 kg (2.5-97.5% interval; median 53.2)",
    sex_female_pct = 37.0,
    race_ethnicity = "Japanese (single-country cohort)",
    disease_state  = paste(
      "Hospitalized adult and pediatric patients (n = 81) with gram-positive",
      "cocci (GPC) or methicillin-resistant Staphylococcus aureus (MRSA)",
      "infections, including sepsis (n = 26), wound / skin / soft-tissue",
      "infection (n = 25), pneumonia (n = 14), abscess (n = 8),",
      "osteomyelitis (n = 6) and undetermined (n = 2)."
    ),
    dose_range     = paste(
      "Linezolid 10 mg/kg three times daily (pediatric) or 300 mg once",
      "daily to 600 mg twice daily (adult), administered orally as",
      "Zyvox film-coated tablets and/or by 1-2 h intravenous infusion;",
      "of 81 patients 54 received only IV, 13 only PO, and 14 both."
    ),
    regions        = paste(
      "Two centres in Japan: Sasebo Chuo Hospital (Nagasaki) and",
      "Toyama University Hospital (Toyama). November 2008 - August 2015."
    ),
    notes          = paste(
      "Demographics from Tsuji 2017 Table 1. Renal function spans CrCl 5.6-188.4",
      "mL/min (median 59.6), with the four youngest pediatric patients (1, 5,",
      "8 and 13 years) assigned RF = 0.5 by the authors. Concentrations were",
      "493 total linezolid + 380 unbound linezolid samples; the platelet PD",
      "fit used 575 platelet counts from 80 patients (one PDS-classified",
      "patient was excluded a posteriori because they already had",
      "thrombocytopenia before starting linezolid)."
    )
  )

  ini({
    # PK structural parameters - reference 70 kg, age 69 years (cohort median),
    # CrCl 100 mL/min/70 kg. All values from Tsuji 2017 Table 2 (final model column).
    lcl_nonren <- log(1.86);  label("Non-renal CL at 70 kg, age 69 (L/h)")          # Tsuji 2017 Table 2 (CL_nonrenal)
    lcl_renal  <- log(1.44);  label("Renal CL at 70 kg, age 69, RF=1 (L/h)")        # Tsuji 2017 Table 2 (CL_renal)
    lvc        <- log(22.9);  label("Central volume at 70 kg (L)")                  # Tsuji 2017 Table 2 (VC)
    lvp        <- log(24.7);  label("Peripheral volume at 70 kg (L)")               # Tsuji 2017 Table 2 (VP)
    lq         <- log(10.9);  label("Inter-compartmental CL at 70 kg (L/h)")        # Tsuji 2017 Table 2 (Q)
    ltabs      <- log(3.61);  label("Absorption half-life (h); ka = ln(2)/Tabs")    # Tsuji 2017 Table 2 (Tabs)
    logitfdepot <- log(0.922 / (1 - 0.922));
    label("Logit-transformed oral bioavailability (F = 0.922)")                     # Tsuji 2017 Table 2 (F)
    logitfu     <- log(0.823 / (1 - 0.823));
    label("Logit-transformed fraction unbound (FU = 0.823; 18% protein binding)")   # Tsuji 2017 Table 2 (FU)

    # Covariate effects on PK
    e_age_cl_nonren <- -0.021;
    label("KAGECL: linear slope of non-renal CL on (AGE - 69), per year")           # Tsuji 2017 Table 2 (KAGECL)

    # Fixed allometric exponents per Tsuji 2017 Methods Equation 3 ("PWR fixed
    # to 0.75 for CL and Q, and 1 for VC and VP"); also fixed for the
    # 70 kg standardisation of CrCl inside RF.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Tsuji 2017 Methods Equation 3
    e_wt_vc_vp <- fixed(1.00); label("Allometric exponent on VC and VP (unitless)") # Tsuji 2017 Methods Equation 3

    # PD structural parameters (Friberg-style turnover; Tsuji 2017 Table 2).
    lmtt     <- log(113);     label("Mean transit time MTT (h); Ntr = 3, Ktr = 4/MTT") # Tsuji 2017 Table 2 (MTT)
    gamma_pd <- -0.187;       label("Empirical feedback exponent gamma (unitless)")  # Tsuji 2017 Table 2 (gamma)
    lpltzero <- log(206000);  label("Baseline circulating platelet count (cells/uL)") # Tsuji 2017 Table 2 (PLTZERO)
    lslope   <- log(0.00566); label("PDI linear slope on total Cc (1/(mg/L))")       # Tsuji 2017 Table 2 (SLOPE)
    lsmax    <- log(2.55);    label("PDS maximum stimulation Emax (unitless)")      # Tsuji 2017 Table 2 (SMAX)
    lsc50    <- log(0.00364); label("PDS SC50 on total Cc (mg/L)")                  # Tsuji 2017 Table 2 (SC50)

    # Inter-individual variability. Tsuji 2017 Table 2 footnote a defines
    # BSV = sqrt(NONMEM OMEGA), so the values reported in the table are
    # omega (the SD on the log-eta scale); omega^2 (the variance entered
    # into nlmixr2 ini()) is the squared value. Tabs, F, SMAX and SC50 are
    # reported as "0 FIXED" in the table and carry no eta here. The IIV on
    # CL is one eta on the composite total CL = (CL_nonren + CL_renal * RF) *
    # FAGE_CL * FSIZE_CL applied as an outer multiplicative factor; the eta
    # name pairs with lcl_nonren for convention-check compliance but does
    # not represent an IIV specific to the non-renal arm. See vignette
    # Assumptions and deviations.
    etalcl_nonren ~ 0.1362  # Tsuji 2017 Table 2 (BSV CL  = 0.369; omega^2 = 0.369^2)
    etalvc     ~ 2.0192  # Tsuji 2017 Table 2 (BSV VC  = 1.421; omega^2 = 1.421^2; high IIV reflects critically ill cohort)
    etalvp     ~ 0.0025  # Tsuji 2017 Table 2 (BSV VP  = 0.050; omega^2 = 0.050^2)
    etalq      ~ 3.3197  # Tsuji 2017 Table 2 (BSV Q   = 1.822; omega^2 = 1.822^2)
    etalmtt    ~ 0.0571  # Tsuji 2017 Table 2 (BSV MTT = 0.239; omega^2 = 0.239^2)
    etagamma_pd ~ 0.0942 # Tsuji 2017 Table 2 (BSV gamma = 0.307; omega^2 = 0.307^2)
    etalpltzero ~ 0.3249 # Tsuji 2017 Table 2 (BSV PLTZERO = 0.570; omega^2 = 0.570^2)
    etalslope  ~ 0.2237  # Tsuji 2017 Table 2 (BSV SLOPE = 0.473; omega^2 = 0.473^2)

    # Residual unidentified variability. Tsuji 2017 Table 2 footnote b notes
    # the RUVs were estimated as THETA-style proportional + additive SDs;
    # propSd / addSd map directly. Output suffixes per nlmixr2lib multi-output
    # convention: Cc (total) bare suffix-free, Cu (unbound) with _Cu suffix,
    # plt (circulating platelets) with _plt suffix.
    propSd     <- 0.318;  label("Proportional residual SD for total Cc (fraction)")       # Tsuji 2017 Table 2 (RUV PROP_TOTAL)
    addSd      <- 0.251;  label("Additive residual SD for total Cc (mg/L)")               # Tsuji 2017 Table 2 (RUV ADD_TOTAL)
    propSd_Cu  <- 0.319;  label("Proportional residual SD for unbound Cu (fraction)")     # Tsuji 2017 Table 2 (RUV PROP_UNBOUND)
    addSd_Cu   <- 0.034;  label("Additive residual SD for unbound Cu (mg/L)")             # Tsuji 2017 Table 2 (RUV ADD_UNBOUND)
    propSd_circ <- 0.234;  label("Proportional residual SD for platelet count (fraction)") # Tsuji 2017 Table 2 (RUV PROP_PLT)
  })

  model({
    # ---- Size factor (allometric on weight, reference 70 kg) ----
    fsize_cl <- (WT / 70)^e_wt_cl_q
    fsize_vc <- (WT / 70)^e_wt_vc_vp

    # ---- Age effect on non-renal CL (linear fractional change about 69 y) ----
    fage_cl <- 1 + e_age_cl_nonren * (AGE - 69)

    # ---- Renal function ratio (size-standardised CrCl normalised to 100 mL/min/70 kg) ----
    crcl_70 <- CRCL * (70 / WT)^e_wt_cl_q
    rf      <- crcl_70 / 100

    # ---- Individual PK parameters ----
    cl_nonren_i <- exp(lcl_nonren) * fage_cl
    cl_renal_i  <- exp(lcl_renal)  * rf
    cl <- (cl_nonren_i + cl_renal_i) * fsize_cl * exp(etalcl_nonren)

    vc <- exp(lvc + etalvc) * fsize_vc
    vp <- exp(lvp + etalvp) * fsize_vc
    q  <- exp(lq  + etalq)  * fsize_cl

    tabs <- exp(ltabs)
    ka   <- log(2) / tabs                        # ka = ln(2) / Tabs per Tsuji 2017 Methods
    fdepot <- 1 / (1 + exp(-logitfdepot))        # back-transform bioavailability
    fu     <- 1 / (1 + exp(-logitfu))            # back-transform fraction unbound

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- PD individual parameters ----
    mtt     <- exp(lmtt + etalmtt)
    ktr     <- 4 / mtt                           # 3 transit compartments => (Ntr + 1)/MTT = 4/MTT
    kcirc   <- ktr                               # paper assumes Kcirc = Ktr
    # gamma is reported on the natural (signed) scale (-0.187). Sign is preserved
    # across individuals by using a multiplicative log-normal eta on |gamma|;
    # additive iiv would let some subjects flip the feedback direction, which is
    # not physiologically meaningful. See vignette Errata.
    gamma   <- gamma_pd * exp(etagamma_pd)
    pltzero <- exp(lpltzero + etalpltzero)
    slope   <- exp(lslope   + etalslope)         # PDI linear slope on total Cc
    smax    <- exp(lsmax)                        # PDS Emax (no IIV per paper)
    sc50    <- exp(lsc50)                        # PDS SC50 (no IIV per paper)

    # ---- PK ODE system: oral depot + central + one peripheral compartment ----
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # ---- PK observations: total Cc and unbound Cu ----
    Cc <- central / vc                           # total linezolid plasma concentration (mg/L)
    Cu <- fu * Cc                                # unbound linezolid plasma concentration (mg/L)

    # ---- Drug-effect terms driving the PD chain ----
    # PDI: linear inhibition of platelet synthesis (drug acts on RFORM).
    # PDS: Emax stimulation of platelet elimination (drug acts on Kcirc).
    # MIX_PDI is a per-subject binary indicator (1 = PDI mechanism, 0 = PDS).
    edrug_syn  <- slope * Cc                                       # PDI linear effect
    edrug_elim <- smax  * Cc / (sc50 + Cc)                         # PDS Emax effect
    rform_factor <- 1 - MIX_PDI       * edrug_syn                  # PDI inhibits synthesis
    elim_factor  <- 1 + (1 - MIX_PDI) * edrug_elim                 # PDS stimulates elimination

    # ---- Feedback and platelet turnover ----
    # FBACK = (circ / PLTZERO)^gamma with gamma negative (-0.187) per Tsuji 2017
    # Table 2 and Discussion ("feedback parameter (gamma) with an absolute
    # value ... were 113 h, 0.187 and 206000 ul-1"). Same sign convention as
    # Sasaki (-0.203) and Boak (-1.02) in Table 4; gives FBACK > 1 (increased
    # formation) when circ < PLTZERO, i.e., compensatory feedback to platelet
    # depletion. Note that Friberg 2002 uses the inverse ratio with positive
    # gamma; the two parametrizations are mathematically equivalent.
    fback <- (circ / pltzero)^gamma

    # Friberg-style turnover chain (compartment names match Friberg_2002_paclitaxel):
    # precursor1 = proliferating-platelet (PLTFORM) compartment; precursor2..4 =
    # three maturation transit compartments; circ = circulating platelet count.
    d/dt(precursor1) <- ktr * precursor1 * rform_factor * fback - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr * precursor4 - kcirc * circ * elim_factor

    precursor1(0) <- pltzero
    precursor2(0) <- pltzero
    precursor3(0) <- pltzero
    precursor4(0) <- pltzero
    circ(0)       <- pltzero

    # ---- Residual error ----
    Cc   ~ add(addSd)    + prop(propSd)
    Cu   ~ add(addSd_Cu) + prop(propSd_Cu)
    circ ~ prop(propSd_circ)
  })
}
