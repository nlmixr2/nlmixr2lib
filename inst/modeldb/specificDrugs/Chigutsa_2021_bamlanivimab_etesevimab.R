Chigutsa_2021_bamlanivimab_etesevimab <- function() {
  description <- paste(
    "Sequential population PK/PD model of the anti-SARS-CoV-2 spike-protein",
    "neutralizing monoclonal antibodies bamlanivimab and etesevimab and of the",
    "SARS-CoV-2 viral-load time course in outpatients with mild to moderate",
    "COVID-19 (BLAZE-1 and BLAZE-4 phase II trials, Chigutsa 2021). Each",
    "antibody has its own two-compartment PK model with linear elimination",
    "after a single intravenous infusion, with body weight entering all",
    "clearances (fixed exponent 0.81) and volumes (fixed exponent 1)",
    "allometrically at a 70 kg reference; the two PK models were fitted",
    "separately and share no parameters. The viral load follows a",
    "target-cell-limited model (Baccam 2006) with uninfected target cells,",
    "productively infected cells and free virus, the model clock running from",
    "the onset of symptoms. Both antibodies increase the virus elimination",
    "rate through additive Emax terms driven by their serum concentrations,",
    "with a common Emax and an etesevimab EC50 fixed at three times the",
    "bamlanivimab EC50 (relative in vitro potency). No covariates other than",
    "body weight on PK were retained in either the PK or the viral-dynamic",
    "analysis.",
    sep = " "
  )
  reference <- paste(
    "Chigutsa E, O'Brien L, Ferguson-Sells L, Long A, Chien J.",
    "Population Pharmacokinetics and Pharmacodynamics of the Neutralizing",
    "Antibodies Bamlanivimab and Etesevimab in Patients With Mild to",
    "Moderate COVID-19 Infection.",
    "Clin Pharmacol Ther. 2021;110(5):1302-1310. doi:10.1002/cpt.2420.",
    "PK parameters are Supplementary Table S2, viral-dynamic parameters are",
    "Table 1, and the viral-dynamic NONMEM control stream is Supplementary",
    "Material S1.",
    sep = " "
  )
  vignette <- "Chigutsa_2021_bamlanivimab_etesevimab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = paste(
      "ug/mL (= mg/L) for the serum antibody concentrations Cc",
      "(bamlanivimab) and Cc_ete (etesevimab); log10 viral load for",
      "log10_viral_load"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "bamlanivimab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "bamlanivimab",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_ete = list(
      analyte = "etesevimab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1_ete = list(
      analyte = "etesevimab",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "free SARS-CoV-2 virus (measured by RT-PCR on nasopharyngeal swabs)",
      units = "1e6 copies/mL (log10(virus * 1e6) is the modelled log10 viral load)",
      specimen = "not applicable",
      verified = TRUE
    ),
    target = list(
      analyte = "uninfected target cells (ACE2-expressing type II pneumocytes)",
      units = "1e6 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    infected = list(
      analyte = "productively infected target cells",
      units = "1e6 cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power on the PK parameters of both antibodies with a",
        "70 kg reference (Supplementary Table S2 footnotes a and b):",
        "CL and Q scale as (WT/70)^0.81 (exponent fixed from Betts 2018),",
        "V1 and V2 as (WT/70)^1. Body weight did not affect the viral",
        "dynamics once its effect on PK was accounted for. The NONMEM",
        "dataset column is WTE (weight at entry)."
      ),
      source_name = "WTE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Tested on PK (stepwise covariate modelling) and on the viral dynamics; not retained."
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Tested on PK; not retained."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested on the viral-dynamic parameters; not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2970L,
    n_studies = 4L,
    age_range = "12-94 years (median 48 years in the bamlanivimab PK, 50 years in the etesevimab PK and viral-load populations)",
    weight_range = "41.4-220 kg (median 85.7 kg bamlanivimab PK, 87.8 kg etesevimab PK and viral-load populations)",
    sex_female_pct = 52.4,
    race_ethnicity = c(White = 85.9, Black = 7.2, Asian = 3.9, Other = 3.0),
    disease_state = paste(
      "Outpatients with mild to moderate COVID-19 (BLAZE-1, BLAZE-4) plus,",
      "for PK only, 18 hospitalized patients with severe COVID-19 (J2W-MC-PYAA)",
      "and 20 healthy participants (J2Z-MC-PGAA)."
    ),
    dose_range = paste(
      "Single IV dose: bamlanivimab 700, 2800 or 7000 mg alone, or",
      "bamlanivimab/etesevimab 175/350, 350/700, 700/1400 or 2800/2800 mg."
    ),
    regions = "Mainly United States (BLAZE-1, BLAZE-4).",
    n_observations = paste(
      "5915 bamlanivimab concentrations from 1899 participants;",
      "4961 etesevimab concentrations from 1498 participants;",
      "17805 viral-load measurements from 2970 participants (BLAZE-1 N = 2303,",
      "BLAZE-4 N = 667)."
    ),
    hepatic_function = "Normal 78.7%, mild impairment 20.2%, moderate 0.1% (viral-load population)",
    notes = paste(
      "Demographics are Supplementary Table S1; the viral-load population is",
      "summarised above (the PK populations are listed in that table",
      "separately). Median baseline viral load 5.30 log10 [(40 - Ct)/log2(10)]."
    )
  )

  ini({
    # ---- Bamlanivimab PK (Supplementary Table S2, bamlanivimab columns) ----
    lcl <- log(0.231); label("Bamlanivimab clearance CL at 70 kg (L/day)")                 # Table S2: CL = 0.231 L/d (%SEE 0.732)
    lvc <- log(2.68);  label("Bamlanivimab central volume V1 at 70 kg (L)")                # Table S2: V1 = 2.68 L (%SEE 0.888)
    lq  <- log(0.281); label("Bamlanivimab intercompartmental clearance Q at 70 kg (L/day)") # Table S2: Q = 0.281 L/d (%SEE 2.85)
    lvp <- log(2.68);  label("Bamlanivimab peripheral volume V2 at 70 kg (L)")             # Table S2: V2 = 2.68 L (%SEE 1.22)
    e_wt_cl_q  <- fixed(0.81); label("Allometric exponent of body weight on bamlanivimab CL and Q (unitless)")   # Table S2: 0.81 (Fixed), footnote a (Betts 2018)
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent of body weight on bamlanivimab V1 and V2 (unitless)") # Table S2: 1.00 (Fixed), footnote b

    # ---- Etesevimab PK (Supplementary Table S2, etesevimab columns) ----
    lcl_ete <- log(0.111); label("Etesevimab clearance CL at 70 kg (L/day)")                 # Table S2: CL = 0.111 L/d (%SEE 0.928)
    lvc_ete <- log(2.45);  label("Etesevimab central volume V1 at 70 kg (L)")                # Table S2: V1 = 2.45 L (%SEE 1.00)
    lq_ete  <- log(0.308); label("Etesevimab intercompartmental clearance Q at 70 kg (L/day)") # Table S2: Q = 0.308 L/d (%SEE 3.96)
    lvp_ete <- log(2.18);  label("Etesevimab peripheral volume V2 at 70 kg (L)")             # Table S2: V2 = 2.18 L (%SEE 1.40)
    e_wt_cl_q_ete  <- fixed(0.81); label("Allometric exponent of body weight on etesevimab CL and Q (unitless)")   # Table S2: 0.81 (Fixed), footnote a
    e_wt_vc_vp_ete <- fixed(1);    label("Allometric exponent of body weight on etesevimab V1 and V2 (unitless)") # Table S2: 1.00 (Fixed), footnote b

    # ---- Viral dynamics (Table 1; NONMEM stream in Supplementary Material S1) ----
    # The stream carries the virus and the cell pools in units of 1e6
    # (IPRED = LOG10(A(3) * 1000000); THETA(1) = 400 FIX for the 4 x 10^8
    # target-cell pool). The Table 1 rate constants are the stream's THETAs
    # and are used on that scale.
    ltarget0 <- fixed(log(400)); label("Uninfected target-cell pool at symptom onset (1e6 cells)")                  # Table 1: 4 x 10^8 (Fixed, Baccam 2006); stream THETA(1) = 400 FIX
    lvirus0  <- log(10^(7.60 - 6)); label("Viral load at symptom onset (1e6 copies/mL; log10 copies/mL = 7.60)")    # Table 1: log10 viral load at onset of symptoms = 7.60 (%RSE 9.05)
    lbeta    <- log(5.23e-7);    label("Infection rate constant beta ((1e6 copies/mL)^-1 day^-1, stream scale)")    # Table 1: beta = 5.23 x 10^-7 (%RSE 31.2)
    lp       <- log(0.0844);     label("Virus production rate PV (day^-1)")                                       # Table 1: PV = 0.0844 day^-1 (%RSE 26.1)
    lc       <- log(1.42);       label("Virus elimination rate CV (day^-1)")                                      # Table 1: CV = 1.42 day^-1 (%RSE 2.30)
    ldelta   <- log(0.290);      label("Death rate of infected cells DI (day^-1)")                                # Table 1: DI = 0.290 day^-1 (%RSE 3.04)
    lemax    <- log(0.462);      label("Maximum antibody-driven increase in the virus elimination rate (day^-1)") # Table 1: Emax = 0.462 (%RSE 3.90), footnote c
    lec50    <- log(0.467);      label("Bamlanivimab EC50 on virus elimination (ug/mL)")                          # Table 1: EC50 = 0.467 ug/mL (%RSE 16.3)
    ec50_ratio_ete <- fixed(3);  label("Etesevimab EC50 as a multiple of the bamlanivimab EC50 (unitless)")       # Results 'Viral dynamic modeling' and Table 1 footnote c: EC50 x 3; stream '(EC50*3)'

    # ---- Interindividual variability ----
    # Table 1 and Table S2 report IIV as %CV; Table 1's 2526% for beta is only
    # consistent with CV = sqrt(exp(omega^2) - 1), and its '15% Fixed' rows are
    # OMEGA 0.0225 FIX in the stream, so omega^2 = log(1 + CV^2) throughout.
    etalcl ~ 0.0611     # Table S2: bamlanivimab CL IIV 25.1% -> log(1 + 0.251^2)
    etalvc ~ 0.0802     # Table S2: bamlanivimab V1 IIV 28.9% -> log(1 + 0.289^2)
    etalcl_ete ~ 0.0776 # Table S2: etesevimab CL IIV 28.4% -> log(1 + 0.284^2)
    etalvc_ete ~ 0.0901 # Table S2: etesevimab V1 IIV 30.7% -> log(1 + 0.307^2)
    etalvirus0 ~ fixed(0.0225) # Table 1: 15% CV, set for SAEM efficiency (footnote b); stream OMEGA BLOCK(1) 0.0225
    # Table 1: IIV beta 2526%, PV 47.2%, CV 85.4% -> variances 6.460, 0.2011,
    # 0.5477; correlations beta-PV 0.939, beta-CV 0.344, PV-CV 0.359.
    etalbeta + etalp + etalc ~ c(
      6.460,
      1.0703, 0.2011,
      0.6471, 0.11916, 0.5477
    )
    etaldelta ~ fixed(0.0225) # Table 1: 15% CV, set for SAEM efficiency (footnote b); stream OMEGA 0.0225
    etalemax  ~ fixed(0.0225) # Table 1: 15% CV, set for SAEM efficiency (footnote b); stream OMEGA 0.0225
    etalec50  ~ fixed(0.0225) # Table 1: 15% CV, set for SAEM efficiency (footnote b); stream OMEGA 0.0225

    # ---- Residual error ----
    propSd     <- 0.190; label("Proportional residual error, bamlanivimab (fraction)") # Table S2: 19.0% (%SEE 2.86)
    propSd_ete <- 0.183; label("Proportional residual error, etesevimab (fraction)")   # Table S2: 18.3% (%SEE 2.42)
    addSd_log10_viral_load <- 0.939; label("Additive residual error on log10 viral load") # Table 1: additive error 0.939 (%RSE 0.656); stream Y = IPRED + THETA(7) * ERR(1), SIGMA 1 FIX
  })

  model({
    # ---- Bamlanivimab PK (Table S2 footnotes a and b) ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Etesevimab PK ----
    cl_ete <- exp(lcl_ete + etalcl_ete) * (WT / 70)^e_wt_cl_q_ete
    vc_ete <- exp(lvc_ete + etalvc_ete) * (WT / 70)^e_wt_vc_vp_ete
    q_ete  <- exp(lq_ete) * (WT / 70)^e_wt_cl_q_ete
    vp_ete <- exp(lvp_ete) * (WT / 70)^e_wt_vc_vp_ete

    kel_ete <- cl_ete / vc_ete
    k12_ete <- q_ete / vc_ete
    k21_ete <- q_ete / vp_ete

    # ---- Viral-dynamic parameters (stream $PK, mu-referenced) ----
    target0 <- exp(ltarget0)
    virus0  <- exp(lvirus0 + etalvirus0)
    beta    <- exp(lbeta + etalbeta)
    p       <- exp(lp + etalp)
    c       <- exp(lc + etalc)
    delta   <- exp(ldelta + etaldelta)
    emax    <- exp(lemax + etalemax)
    ec50    <- exp(lec50 + etalec50)
    ec50_ete <- ec50 * ec50_ratio_ete

    # Serum concentrations (mg / L = ug/mL)
    Cc     <- central / vc
    Cc_ete <- central_ete / vc_ete

    # Table 1 footnote c: CV_t = CV + Emax*conc1/(EC50 + conc1)
    #                            + Emax*conc2/(3*EC50 + conc2), Hill = 1
    eff     <- emax * Cc / (ec50 + Cc)
    eff_ete <- emax * Cc_ete / (ec50_ete + Cc_ete)

    d/dt(central)         <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1
    d/dt(central_ete)     <- -kel_ete * central_ete - k12_ete * central_ete + k21_ete * peripheral1_ete
    d/dt(peripheral1_ete) <-  k12_ete * central_ete - k21_ete * peripheral1_ete

    # Target-cell-limited viral dynamics (Methods equations; stream $DES)
    d/dt(virus)    <- p * infected - c * virus - eff * virus - eff_ete * virus
    d/dt(target)   <- -beta * target * virus
    d/dt(infected) <- beta * target * virus - delta * infected

    # Time 0 is the onset of symptoms; infected cells start at zero.
    virus(0)    <- virus0
    target(0)   <- target0
    infected(0) <- 0

    # Stream $ERROR: IPRED = LOG10(A(3) * 1e6) when A(3) * 1e6 > 1, else 1e-5
    vl_linear <- virus * 1e6
    log10_viral_load <- 0.00001
    if (vl_linear > 1) log10_viral_load <- log10(vl_linear)

    Cc ~ prop(propSd)
    Cc_ete ~ prop(propSd_ete)
    log10_viral_load ~ add(addSd_log10_viral_load)
  })
}
