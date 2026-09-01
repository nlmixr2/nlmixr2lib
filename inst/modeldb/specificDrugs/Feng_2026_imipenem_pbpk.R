Feng_2026_imipenem_pbpk <- function() {
  description <- paste(
    "PBPK-derived reduced one-compartment intravenous model for imipenem across",
    "healthy adults, adults with renal impairment (RI), and children aged 3-18",
    "years with normal renal function or mild / moderate / severe RI. The source",
    "paper built a whole-body PBPK model in GastroPlus 9.9, using the software's",
    "age-related physiology module to generate organ weights, blood flows, cardiac",
    "output and in vivo renal clearance; that whole-body structure is a GastroPlus",
    "platform model whose ODEs, organ volumes, blood flows and tissue-to-plasma",
    "partition coefficients are NOT written out in the publication and are",
    "therefore not reproduced here. What IS fully reported, and is what this file",
    "encodes, is the reduced disposition model the authors carried out of the PBPK",
    "and into their dose-adjustment and Monte Carlo target-attainment analysis:",
    "Supplementary Table S3 tabulates total clearance, steady-state volume of",
    "distribution, half-life and unbound fraction for every simulated population,",
    "and those triples satisfy t(1/2) = ln2 * Vss / CL to within 0.5% on every",
    "row, i.e. the paper's own summary of its PBPK is mono-exponential. Volume is",
    "weight-normalized and shared by every stratum (Vss / body weight is",
    "0.155-0.173 L/kg across all Table S3 rows with a published weight, spanning",
    "58.75-86.6 kg adults and, independently, an 8-year-old child). Clearance is",
    "stratum-specific: adult clearance is absolute and body-weight-independent",
    "(10.7-11.3 L/h over a 58.75-86.6 kg range, because it is driven by glomerular",
    "filtration), while pediatric clearance is weight-normalized because pediatric",
    "dosing is per kg. Intravenous infusion dosing; linear elimination; no",
    "absorption model. The PK/PD index is the percentage of the dosing interval",
    "during which the free concentration exceeds the MIC (%fT>MIC) with 70%fT>MIC",
    "as the target and a probability of target attainment of at least 80% as the",
    "acceptance criterion; %fT>MIC is computed from the simulated profile using the",
    "unbound fraction carried in the model rather than integrated as a model state,",
    "following the Yang_2024_meropenem_pbpk precedent. The paper ran a 200-subject",
    "population simulator and reports 90% confidence intervals graphically, but",
    "reports no interindividual variance magnitudes and no residual error model, so",
    "all etas are omitted and both residual error terms are encoded as fixed(0).",
    sep = " "
  )
  reference <- paste(
    "Feng C, Xiao P, Qu Y, Fan K, Wang Y, Zhang X, Wang X, Pan J, Deng Y, Yu Y.",
    "Physiologically based pharmacokinetic modeling and dose adjustment of",
    "imipenem in pediatric patients with renal impairment.",
    "Front Cell Infect Microbiol. 2026;16:1798911.",
    "doi:10.3389/fcimb.2026.1798911. PMCID PMC13194417.",
    "Clearance, steady-state volume of distribution, half-life, GFR and unbound",
    "fraction for every simulated population are Supplementary Table S3; the",
    "clinical studies that supplied the observed data and the per-study",
    "demographics are Supplementary Tables S1 (adults) and S2 (children).",
    "Predicted and observed AUC0-t and Cmax are Table 2 (adults) and Table 3",
    "(children). Pediatric AUC0-24h by WHO age band and renal-function stratum,",
    "the geometric mean ratios and the proposed doses are Table 4. The WHO age",
    "bands (preschool 3-6, school-age 6-12, adolescent 12-18 years) and the FDA",
    "renal-function classification (GFR >= 90, 60-89, 30-59, 15-29",
    "mL/min/1.73 m^2) are Methods. The PK/PD target 70%fT>MIC, the >= 80% PTA",
    "acceptance criterion and the infusion-duration comparison are Methods",
    "'Pharmacodynamic evaluation of imipenem'; the resulting PTA curves are",
    "Figure 6.",
    sep = " "
  )
  vignette <- "Feng_2026_imipenem_pbpk"
  # The paper tabulates concentrations as ug/mL and AUC as ug.h/mL; mg/L is the
  # numerically identical form and matches the sibling Yang_2024_meropenem_pbpk.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales the volume of distribution in every stratum and, in the pediatric",
        "strata only, the clearance. Adult clearance is deliberately NOT scaled by",
        "body weight: Supplementary Table S3 reports 10.726-11.262 L/h across adult",
        "cohorts weighing 58.75-86.6 kg, i.e. essentially flat, because the",
        "GastroPlus renal-clearance term is driven by glomerular filtration rate",
        "rather than by size. Adult weights are Supplementary Table S1 (58.75-86.6",
        "kg means). Supplementary Table S2 reports a mean weight only for the",
        "Bradley 2023 pediatric cohort (15.4 kg); see the vignette Errata for why",
        "that value is inconsistent with the subject the paper actually simulated.",
        sep = " "
      ),
      source_name        = "Weight, kg"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Selects the population stratum via the WHO age bands the paper adopts in",
        "Methods 'Prediction of imipenem exposures in pediatric patients with RI':",
        "preschool 3-6 years, school-age 6-12 years, adolescent 12-18 years, adult",
        "18 years and over. The representative simulated children were 3, 8 and 16",
        "years old. The paper states in Limitations that the model applies only to",
        "children over 3 years old with normal growth and development, and not to",
        "neonates, infants under 3 years, or children who are overweight or",
        "underweight; an AGE below 3 falls into no band and yields a clearance of",
        "zero rather than an error, so callers must screen for it.",
        sep = " "
      ),
      source_name        = "Age, years"
    ),
    RENALIMP_MILD = list(
      description        = "Mild renal impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = paste(
        "1 = mild renal impairment, 0 = normal renal function or any non-mild",
        "category. The paper classifies renal function by GFR per FDA guidance",
        "(Methods 'Model for adult patients with RI'): normal >= 90, mild 60-89,",
        "moderate 30-59 and severe 15-29 mL/min/1.73 m^2. Mutually exclusive with",
        "RENALIMP_MOD and RENALIMP_SEV; all three are 0 for normal renal function.",
        "The adult effect is transcribed from the Gibson 1985 renal-insufficiency",
        "cohort simulated in Supplementary Table S3; the pediatric effects are the",
        "software's built-in renal-impairment physiology and are, as the paper",
        "states in Limitations, an unvalidated extrapolation because no pediatric",
        "renal-impairment PK data exist.",
        sep = " "
      ),
      source_name        = "Mild RI"
    ),
    RENALIMP_MOD = list(
      description        = "Moderate renal impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = paste(
        "1 = moderate renal impairment (GFR 30-59 mL/min/1.73 m^2 per the FDA",
        "classification the paper adopts in Methods 'Model for adult patients with",
        "RI'), 0 otherwise. Mutually exclusive with RENALIMP_MILD and RENALIMP_SEV.",
        sep = " "
      ),
      source_name        = "Moderate RI"
    ),
    RENALIMP_SEV = list(
      description        = "Severe renal impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = paste(
        "1 = severe renal impairment (GFR 15-29 mL/min/1.73 m^2 per the FDA",
        "classification the paper adopts in Methods 'Model for adult patients with",
        "RI'), 0 otherwise. Mutually exclusive with RENALIMP_MILD and RENALIMP_MOD.",
        "Note that the Gibson 1985 severe-RI cohort the adult effect is transcribed",
        "from was simulated at GFR 0.1433 mL/s = 8.6 mL/min, which sits below the",
        "15-29 mL/min band the paper's own classification assigns to severe RI.",
        sep = " "
      ),
      source_name        = "Severe RI"
    )
  )

  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 76,
    n_studies      = 8,
    age_range      = paste(
      "Adults 18-68 years (study means 25-60.8 years, Supplementary Table S1);",
      "children 2-12 years (Supplementary Table S2). Simulated pediatric",
      "representatives were 3, 8 and 16 years old.",
      sep = " "
    ),
    weight_range   = paste(
      "Adult study means 58.75-86.6 kg (individual range 51-116 kg,",
      "Supplementary Table S1); Bradley 2023 pediatric cohort mean 15.4 kg",
      "(range 12-19), Claesson 1992 weights not reported (Supplementary Table S2).",
      sep = " "
    ),
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Healthy adult volunteers; adults with mild, moderate or severe renal",
      "insufficiency (Gibson 1985); children with peritonitis (Claesson 1992) and",
      "children with confirmed or suspected gram-negative bacterial infections",
      "(Bradley 2023). The pediatric renal-impairment populations are simulated",
      "only -- no pediatric renal-impairment PK data exist and the paper says so",
      "in Limitations.",
      sep = " "
    ),
    dose_range     = paste(
      "Adults 0.25-1.0 g as an intravenous bolus or a 0.15, 0.5 or 2.0 h infusion,",
      "single dose and q6h; children 15 and 25 mg/kg as a 0.5 or 1.0 h infusion.",
      "Proposed pediatric regimens are 15 mg/kg q6h (normal renal function and mild",
      "RI), 12 mg/kg q6h (moderate RI) and 7 mg/kg q6h (severe RI), each as a 3 h",
      "infusion (Table 4 and Results 'Pharmacodynamic evaluation of imipenem').",
      sep = " "
    ),
    regions        = "China, Thailand, Sweden, United States",
    renal_function = paste(
      "Normal (GFR >= 90), mild (60-89), moderate (30-59) and severe (15-29)",
      "renal impairment per the FDA classification, in mL/min/1.73 m^2",
      "(Methods 'Model for adult patients with RI').",
      sep = " "
    ),
    co_medication  = paste(
      "All included clinical studies administered imipenem with cilastatin or with",
      "another enzyme inhibitor, and all reported imipenem concentrations measured",
      "during that co-administration (Methods 'Model development and evaluation').",
      "The renal clearance the model carries is therefore the cilastatin-inhibited",
      "value, as Table 1's footnote states explicitly.",
      sep = " "
    ),
    notes          = paste(
      "Eight clinical studies pooled from the literature and digitized with GetData",
      "Graph Digitizer 2.26: six adult (Wang 2021, Jaruratanasirikul 2005, Nilsson",
      "1991, Norrby 1983, Drusano 1984, Gibson 1985) in Supplementary Table S1 and",
      "two pediatric (Claesson 1992, Bradley 2023) in Supplementary Table S2. The",
      "n_subjects total of 76 sums the per-study N of those tables. Population",
      "simulations used 200 virtual subjects per trial at a 1:1 male-to-female",
      "ratio (Methods 'Population simulation').",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Volume of distribution, weight-normalized and shared by every stratum.
    #
    # Supplementary Table S3 reports Vss under a "Vss(L/kg)" header, and the
    # Results text repeats it as "about 3.259 L/kg in children, compared with
    # 11.663 L/kg in adults". That unit is wrong: 11.663 L/kg would be 863 L
    # in a 74 kg adult and an 11 h half-life for a drug the same table gives a
    # 0.723 h half-life. The values are absolute litres, which is proved three
    # ways -- t(1/2) = ln2 * Vss / CL holds to within 0.5% on every row only
    # with Vss in L and CL in L/h; dose / Vss reproduces the paper's own
    # predicted bolus Cmax scale; and Vss / body weight is then a physiological
    # 0.155-0.173 L/kg. See the vignette Errata.
    #
    # 0.1636 is the mean of Vss / body weight over the seven Table S3 rows that
    # have a published weight in Table S1: 0.1730 (Wang 2021, 60.4 kg), 0.1555
    # (Jaruratanasirikul 2005, 58.75), 0.1576 (Nilsson 1991, 74), 0.1649
    # (Norrby 1983, 75), 0.1626 (Gibson mild, 86.6), 0.1605 (Gibson moderate,
    # 76.5), 0.1710 (Gibson severe, 71.3). Independently, solving Table 3's
    # four school-age (AUC0-t, Cmax) pairs for a 0.5 h infusion of 15 mg/kg
    # gives 0.1629 / 0.1635 / 0.1659 / 0.1652 L/kg for normal / mild / moderate
    # / severe RI -- an 8-year-old child and an 87 kg adult on the same value.
    # Renal impairment therefore carries no effect on volume; the 11.663 to
    # 14.081 L spread across the adult renal strata is body weight, not RI.
    # ---------------------------------------------------------------------
    lvcPerKg <- fixed(log(0.1636)) ; label("Volume of distribution per kg body weight (L/kg)")  # Supplementary Table S3 Vss divided by the Supplementary Table S1 weights; see comment above

    # ---------------------------------------------------------------------
    # Adult clearance, absolute and body-weight-independent.
    #
    # Supplementary Table S3 gives 10.726-11.262 L/h across five healthy adult
    # cohorts spanning 58.75-86.6 kg, so clearance does not track size; the
    # Discussion quotes "approximately 11.126 L/h" as the theoretical healthy
    # adult value. 11.176 L/h is the Table S3 row that pairs with the 11.663 L
    # volume (Wang 2021 q6h and Drusano 1984 q6h) and is used as the reference.
    # ---------------------------------------------------------------------
    lcl_adult <- fixed(log(11.176)) ; label("Clearance, adult with normal renal function (L/h)")  # Supplementary Table S3, Wang 2021 q6h / Drusano 1984 rows

    # Adult renal-impairment multipliers on clearance. Numerators are the
    # Gibson 1985 rows of Supplementary Table S3 (10.682 / 6.376 / 3.486 L/h);
    # the denominator is the 11.176 L/h reference above. Gibson's own
    # normal-renal control, implied by Table 2's predicted AUC0-t of 22.33
    # mg.h/mL after 0.25 g, is 11.196 L/h -- within 0.2% of the reference, so
    # the ratios are not sensitive to which healthy row is used.
    e_renalimp_mild_cl_adult <- fixed(0.95580) ; label("Clearance multiplier for mild renal impairment, adults (unitless)")      # Supplementary Table S3 Gibson mild CL 10.682 / 11.176
    e_renalimp_mod_cl_adult  <- fixed(0.57051) ; label("Clearance multiplier for moderate renal impairment, adults (unitless)")  # Supplementary Table S3 Gibson moderate CL 6.376 / 11.176
    e_renalimp_sev_cl_adult  <- fixed(0.31192) ; label("Clearance multiplier for severe renal impairment, adults (unitless)")    # Supplementary Table S3 Gibson severe CL 3.486 / 11.176

    # ---------------------------------------------------------------------
    # Pediatric clearance, weight-normalized, one value per WHO age band.
    #
    # Table 4 reports AUC0-24h after a single 15 mg/kg dose for each band and
    # renal stratum. The paper's half-lives are 0.29-0.70 h, so AUC0-24h is
    # AUC0-inf for a single dose and CL/WT = 15 / AUC0-24h exactly. The
    # school-age value below reproduces Table 3's school-age row (predicted
    # AUC0-t 42.82 mg.h/mL for the 8-year-old representative) identically, which
    # confirms Table 4 is single-dose rather than q6h steady state.
    # ---------------------------------------------------------------------
    lclPerKg_preschool  <- fixed(log(0.47835)) ; label("Clearance per kg, preschool children 3-6 years, normal renal function (L/h/kg)")   # Table 4 preschool normal, 15 / 31.358
    lclPerKg_schoolage  <- fixed(log(0.35032)) ; label("Clearance per kg, school-age children 6-12 years, normal renal function (L/h/kg)") # Table 4 school-age normal, 15 / 42.818
    lclPerKg_adolescent <- fixed(log(0.17280)) ; label("Clearance per kg, adolescents 12-18 years, normal renal function (L/h/kg)")        # Table 4 adolescent normal, 15 / 86.804

    # Pediatric renal-impairment multipliers on clearance, one set per band.
    # Each is the reciprocal of that band's AUC0-24h ratio in Table 4. The
    # ratios are computed from the tabulated AUC values rather than from the
    # tabulated geometric mean ratios because one printed GMR disagrees with
    # its own AUC pair: Table 4 gives 1.19 for preschool moderate RI, but
    # 38.928 / 31.358 = 1.2414 (4.3% apart). The other eight agree to within
    # 0.9%. See the vignette Errata.
    e_renalimp_mild_cl_preschool  <- fixed(0.96246) ; label("Clearance multiplier for mild renal impairment, preschool children (unitless)")      # Table 4, 31.358 / 32.581
    e_renalimp_mod_cl_preschool   <- fixed(0.80554) ; label("Clearance multiplier for moderate renal impairment, preschool children (unitless)")  # Table 4, 31.358 / 38.928
    e_renalimp_sev_cl_preschool   <- fixed(0.49650) ; label("Clearance multiplier for severe renal impairment, preschool children (unitless)")    # Table 4, 31.358 / 63.158
    e_renalimp_mild_cl_schoolage  <- fixed(0.95723) ; label("Clearance multiplier for mild renal impairment, school-age children (unitless)")     # Table 4, 42.818 / 44.731
    e_renalimp_mod_cl_schoolage   <- fixed(0.79639) ; label("Clearance multiplier for moderate renal impairment, school-age children (unitless)") # Table 4, 42.818 / 53.765
    e_renalimp_sev_cl_schoolage   <- fixed(0.46826) ; label("Clearance multiplier for severe renal impairment, school-age children (unitless)")   # Table 4, 42.818 / 91.441
    e_renalimp_mild_cl_adolescent <- fixed(0.94831) ; label("Clearance multiplier for mild renal impairment, adolescents (unitless)")             # Table 4, 86.804 / 91.535
    e_renalimp_mod_cl_adolescent  <- fixed(0.74580) ; label("Clearance multiplier for moderate renal impairment, adolescents (unitless)")         # Table 4, 86.804 / 116.39
    e_renalimp_sev_cl_adolescent  <- fixed(0.42434) ; label("Clearance multiplier for severe renal impairment, adolescents (unitless)")           # Table 4, 86.804 / 204.56

    # ---------------------------------------------------------------------
    # Unbound fraction in plasma, used to form the free concentration that
    # drives the paper's %fT>MIC endpoint. Supplementary Table S3 reports one
    # value per simulated population. The healthy value of 0.8 is the Table 1
    # literature input; the renal-impairment values are GastroPlus defaults
    # that rise with impairment, which the Discussion attributes to the fall in
    # serum albumin in renal disease. Supplementary Table S3 carries NO
    # pediatric renal-impairment row, so fup_ped is applied to children of
    # every renal stratum; see the vignette Errata.
    # ---------------------------------------------------------------------
    fup_adult      <- fixed(0.80000) ; label("Unbound fraction in plasma, adult with normal renal function (fraction)")  # Table 1 Fup 80%, and Supplementary Table S3 healthy adult rows
    fup_adult_mild <- fixed(0.80588) ; label("Unbound fraction in plasma, adult with mild renal impairment (fraction)")      # Supplementary Table S3 Gibson mild
    fup_adult_mod  <- fixed(0.82723) ; label("Unbound fraction in plasma, adult with moderate renal impairment (fraction)")  # Supplementary Table S3 Gibson moderate
    fup_adult_sev  <- fixed(0.83352) ; label("Unbound fraction in plasma, adult with severe renal impairment (fraction)")    # Supplementary Table S3 Gibson severe
    fup_ped        <- fixed(0.81217) ; label("Unbound fraction in plasma, children (fraction)")                              # Supplementary Table S3 Claesson 1992 row; the Bradley 2023 row gives 0.81182

    # ---------------------------------------------------------------------
    # Residual error. The paper ran a 200-subject population simulator and
    # displays 90% confidence intervals in Figures 2-4, but reports no
    # interindividual variance magnitudes anywhere in the paper or its
    # supplement, and describes no residual error model. Both terms are
    # therefore fixed at zero rather than invented, and no etas are declared;
    # see the vignette Errata.
    # ---------------------------------------------------------------------
    addSd  <- fixed(0) ; label("Additive residual error (mg/L; not published by Feng 2026)")
    propSd <- fixed(0) ; label("Proportional residual error (fraction; not published by Feng 2026)")
  })

  model({
    # -------------------------------------------------------------------
    # 1. Population strata, from the WHO age bands the paper adopts in
    #    Methods "Prediction of imipenem exposures in pediatric patients
    #    with RI" (preschool 3-6, school-age 6-12, adolescent 12-18 years)
    #    and the adult cohorts of Supplementary Table S1.
    #
    #    The bands are exclusive and exhaustive for AGE >= 3. An AGE below 3
    #    falls into no band, giving a clearance of zero and a profile that
    #    never declines, rather than an error -- the paper states in
    #    Limitations that it does not apply to neonates or children under 3
    #    years, so callers must screen for that themselves.
    # -------------------------------------------------------------------
    isPreschool  <- (AGE >= 3) * (AGE < 6)
    isSchoolage  <- (AGE >= 6) * (AGE < 12)
    isAdolescent <- (AGE >= 12) * (AGE < 18)
    isAdult      <- (AGE >= 18)

    # -------------------------------------------------------------------
    # 2. Renal-function stratum. The three indicators are mutually
    #    exclusive; normal renal function is all three set to 0. A data set
    #    that sets two of them to 1 makes isNormalRenal negative and
    #    produces nonsense rather than erroring, so assemblers must preserve
    #    exclusivity.
    # -------------------------------------------------------------------
    isNormalRenal <- 1 - RENALIMP_MILD - RENALIMP_MOD - RENALIMP_SEV

    riAdult <- isNormalRenal +
      RENALIMP_MILD * e_renalimp_mild_cl_adult +
      RENALIMP_MOD  * e_renalimp_mod_cl_adult +
      RENALIMP_SEV  * e_renalimp_sev_cl_adult
    riPreschool <- isNormalRenal +
      RENALIMP_MILD * e_renalimp_mild_cl_preschool +
      RENALIMP_MOD  * e_renalimp_mod_cl_preschool +
      RENALIMP_SEV  * e_renalimp_sev_cl_preschool
    riSchoolage <- isNormalRenal +
      RENALIMP_MILD * e_renalimp_mild_cl_schoolage +
      RENALIMP_MOD  * e_renalimp_mod_cl_schoolage +
      RENALIMP_SEV  * e_renalimp_sev_cl_schoolage
    riAdolescent <- isNormalRenal +
      RENALIMP_MILD * e_renalimp_mild_cl_adolescent +
      RENALIMP_MOD  * e_renalimp_mod_cl_adolescent +
      RENALIMP_SEV  * e_renalimp_sev_cl_adolescent

    # -------------------------------------------------------------------
    # 3. Clearance. Adults take an absolute, body-weight-independent value;
    #    children take a weight-normalized one, matching how each is
    #    reported (Supplementary Table S3 in L/h, Table 4 as AUC after a
    #    mg/kg dose).
    # -------------------------------------------------------------------
    clPerKgPed <-
      isPreschool  * exp(lclPerKg_preschool)  * riPreschool +
      isSchoolage  * exp(lclPerKg_schoolage)  * riSchoolage +
      isAdolescent * exp(lclPerKg_adolescent) * riAdolescent

    cl <- isAdult * exp(lcl_adult) * riAdult + clPerKgPed * WT

    # -------------------------------------------------------------------
    # 4. Volume. Shared by every stratum; renal impairment carries no
    #    effect on volume (see the ini() comment).
    # -------------------------------------------------------------------
    vc <- exp(lvcPerKg) * WT

    # -------------------------------------------------------------------
    # 5. Unbound fraction, for the free concentration that drives %fT>MIC.
    # -------------------------------------------------------------------
    fup <- isAdult * (
      isNormalRenal * fup_adult +
        RENALIMP_MILD * fup_adult_mild +
        RENALIMP_MOD  * fup_adult_mod +
        RENALIMP_SEV  * fup_adult_sev
    ) + (1 - isAdult) * fup_ped

    # -------------------------------------------------------------------
    # 6. One-compartment intravenous disposition. Imipenem is given only by
    #    intravenous bolus or infusion in every study the paper used, so
    #    there is no absorption model and no bioavailability term.
    # -------------------------------------------------------------------
    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc    <- central / vc
    Cfree <- Cc * fup

    Cc ~ add(addSd) + prop(propSd)
  })
}
