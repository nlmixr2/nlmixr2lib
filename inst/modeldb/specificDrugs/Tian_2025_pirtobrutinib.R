Tian_2025_pirtobrutinib <- function() {
  description <- paste(
    "Two-compartment oral pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK-with-single-adjusting-compartment (SAC) model for the",
    "reversible BTK inhibitor pirtobrutinib in healthy adults",
    "(Tian 2025). The source model was built in the Simcyp Simulator",
    "version 19 and its whole-body mass-balance equations are not",
    "published, so the platform model itself cannot be encoded here.",
    "What IS fully reported is the pirtobrutinib compound layer, and it",
    "is sufficient to rebuild the disposition as an ordinary",
    "compartmental model: first-order absorption into a depot with a lag",
    "time, distribution between a systemic compartment and the SAC (the",
    "Simcyp kin / kout, encoded as the canonical k12 / k21), and",
    "first-order elimination from the systemic compartment. Systemic",
    "clearance is carried as the three components the paper's disposition",
    "diagram prints separately, so a CYP3A4 perturbation can be applied",
    "to the CYP3A4-mediated arm alone: cl_3a4 (40% of systemic",
    "clearance), cl_other (non-CYP3A4 hepatic metabolism plus biliary",
    "excretion) and cl_renal. Hepatic availability is not a fitted",
    "constant but is derived inside the model from the well-stirred",
    "relationship the paper states as its Equation 2, so that oral",
    "bioavailability responds correctly when the CYP3A4 arm is changed.",
    "No parameter is fitted here: every value is either a Table 1 /",
    "Figure 2 / Data S3 input or an arithmetic consequence of one. The",
    "reduction reproduces the paper's own predicted median tmax exactly",
    "after both a single 200 mg dose (2.94 h) and 200 mg once daily",
    "(2.72 h), and its predicted Cmax and AUC to within 8.3%; the",
    "residual offset is the difference between a typical-value",
    "individual and the geometric mean of the paper's 100-subject",
    "virtual population (see the validation vignette).",
    "This is a typical-value simulation model: the source reports no",
    "estimated inter-individual variance components and no residual-error",
    "model, so there are no etas and propSd is fixed at zero.",
    "The drug-drug-interaction predictions that are the paper's main",
    "contribution depend on proprietary Simcyp perpetrator compound",
    "files and on a CYP3A4 degradation rate constant that is never",
    "printed, so they are NOT reproducible mechanistically from this",
    "model; the vignette instead shows the closed-form CYP3A4-fraction",
    "envelope this reduction does support and compares it against the",
    "paper's Figure 5.",
    sep = " "
  )
  reference <- paste(
    "Tian DD, Hall SD, Chapman SC, Posada MM. (2025).",
    "Application of physiologically-based pharmacokinetic modeling to",
    "support drug labeling: prediction of CYP3A4-mediated",
    "pirtobrutinib-drug interactions.",
    "CPT Pharmacometrics Syst Pharmacol 14:2221-2231.",
    "doi:10.1002/psp4.70134.",
    sep = " "
  )
  vignette <- "Tian_2025_pirtobrutinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Tian 2025 Figure 2 (the proposed
  # disposition diagram) and the Data S3 Simcyp input sheet rows
  # "Distribution Model | Minimal PBPK Model", "SAC kin (1/h)",
  # "SAC kout (1/h)" and "SAC Q (L/h)".
  compartmentData <- list(
    depot = list(
      analyte = "pirtobrutinib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "pirtobrutinib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "pirtobrutinib", units = "mg",
      specimen = "tissue", verified = TRUE
    )
  )

  # No covariates are carried. The Simcyp Sim-Healthy Volunteers
  # population behind Tian 2025 varies body weight, age and sex across
  # the 100 virtual subjects, but those act through Simcyp population
  # files that are not reported, so no covariate relationship in the
  # published paper is reproducible here. Recorded as
  # screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. The source expresses Vss (0.54 L/kg) and Vsac",
        "(0.17 L/kg) per kilogram in Table 1, so the Simcyp model does",
        "scale distribution volume with weight. This reduction does not",
        "need a weight term at all, because Data S3 reports the SAC",
        "rate constants kin and kout together with the",
        "inter-compartmental clearance Q, and vc = Q / kin is an",
        "absolute volume in litres. Carrying a weight exponent on vc",
        "would therefore be a guess about the platform rather than a",
        "reported relationship."
      ),
      units = "kg",
      type  = "continuous",
      notes = "Implicit in the L/kg volume inputs; not needed by this reduction."
    ),
    AGE = list(
      description = paste(
        "Age. Data S3 sets the Sim-Healthy Volunteers age range to",
        "20-55 years for the multiple-dose simulation, and Table S4",
        "records 19-55 years for the single-dose simulation, matched to",
        "the contributing clinical studies. Age drives liver weight and",
        "CYP3A4 abundance inside the Simcyp population files, but no",
        "age relationship is printed, so none is carried."
      ),
      units = "years",
      type  = "continuous",
      notes = "Sets the virtual-population sampling range only; no printed relationship."
    ),
    SEXF = list(
      description = paste(
        "Female sex indicator (1 = female). Data S3 sets the proportion",
        "of females to 0.25 for the multiple-dose simulation and Table",
        "S4 records 0.18 for the single-dose simulation. Sex drives",
        "organ weights and blood flows inside the Simcyp population",
        "files; no sex effect on any pirtobrutinib parameter is printed."
      ),
      units = "unitless",
      type  = "categorical",
      notes = "Sets the virtual-population sex split only; no printed relationship."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 9L,
    age_range      = "20-55 years for the 200 mg once-daily simulation (Data S3); 19-55 years for the single-dose simulation (Table S4)",
    weight_median  = "not reported; the Simcyp Sim-Healthy Volunteers default distribution was used",
    sex_female_pct = 25,
    disease_state  = "Healthy adult volunteers, fasted, CYP3A5 genotype unrestricted for the pirtobrutinib-alone simulations.",
    dose_range     = "Single oral dose of 200 mg, and 200 mg once daily for 13 doses (Data S3).",
    regions        = "Simcyp Sim-Healthy Volunteers virtual population (Data S3).",
    studies        = paste(
      "Nine clinical pharmacology studies supplied the observations and",
      "the compound-layer inputs (Supplemental Material 1).",
      "The absolute bioavailability study (n = 5 healthy males, age",
      "19-32, NCT06180954) is the study that anchors this reduction:",
      "it supplies the intravenous systemic clearance of 1.69 L/h, the",
      "renal clearance of 0.12 L/h, the steady-state volume of 0.54",
      "L/kg and the absolute oral bioavailability of 0.86.",
      "The mass balance study (n = 4 healthy males, age 35-43,",
      "NCT06180954) supplies the excretion split that fixes f_e,urine",
      "= 0.071 and f_e,bile = 0.14 and the metabolite recoveries that",
      "corroborate f_m,other = 0.39.",
      "The pirtobrutinib-itraconazole and pirtobrutinib-rifampin",
      "interaction studies (NCT05134337) are the studies ka, t_lag and",
      "f_m,CYP3A4 were fitted against.",
      "The pirtobrutinib-midazolam interaction study (n = 15, age",
      "22-55, NCT06180967) supplied the multiple-dose profiles the",
      "time-dependent-inhibition k_inact was fitted against.",
      "The pilot and formal food effect studies (NCT05134350,",
      "NCT06180980), the CYP cocktail study (NCT06215430) and the",
      "repaglinide study (NCT06165146) are model verification studies."
    ),
    notes          = paste(
      "n_subjects records the 100 virtual subjects (10 trials of 10)",
      "Data S3 simulated for pirtobrutinib alone, because this is a PBPK",
      "analysis rather than a population-PK fit: there is no pooled",
      "analysis dataset and no estimated variance components.",
      "n_studies counts the nine clinical studies listed above."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Tian 2025 Table 1 /
    # Figure 2 / Data S3 inputs, or arithmetic consequences of them.
    #
    # The three shared physiological constants, all main text:
    #   B:P = 0.79     blood-to-plasma ratio (Table 1, measured)
    #   Q_H = 90 L/h   hepatic blood flow (main text after Equation 2)
    #   Q   = 3.56 L/h SAC inter-compartmental clearance (Table 1)
    #
    # Figure 2 (the proposed disposition diagram) prints the whole
    # clearance decomposition, and it is internally exact:
    #   CL                        = 1.69 L/h
    #   CL_R                      = 0.12 L/h,  f_e,urine = 0.071
    #   CL_H = CL - CL_R          = 1.57 L/h            (Equation 1)
    #   CL_bile                   = 0.24 L/h,  f_e,bile = 0.14
    #   CL_hepatic metabolism     = 1.33 L/h
    #   f_m,CYP3A4 = 0.40, f_m,other = 0.39
    # so 0.40 * 1.69 = 0.676 L/h is CYP3A4-mediated, the non-CYP3A4
    # hepatic metabolism is 1.33 - 0.676 = 0.654 L/h, and the three
    # components below sum back to 1.69 L/h exactly.
    # ------------------------------------------------------------------

    lka <- fixed(log(0.82))
    label("First-order absorption rate constant ka (1/h)")
    # Tian 2025 Table 1: ka = 0.82 1/h, "manually fitted to match
    # observed tmax and Cmax after a single dose of pirtobrutinib alone
    # in the pirtobrutinib-itraconazole and pirtobrutinib-rifampin
    # interaction studies". Data S3 row "ka (1/h)" confirms 0.82.

    ltlag <- fixed(log(0.25))
    label("Absorption lag time t_lag (h)")
    # Tian 2025 Table 1: t_lag = 0.25 h, fitted alongside ka in the same
    # sensitivity analysis. Data S3 row "lag time (h)" confirms 0.25.

    lcl_3a4 <- fixed(log(0.676))
    label("CYP3A4-mediated systemic plasma clearance cl_3a4 (L/h)")
    # Derived from Tian 2025 Figure 2: f_m,CYP3A4 = 0.40 times the
    # systemic clearance CL = 1.69 L/h gives 0.676 L/h. f_m,CYP3A4 was
    # itself estimated by sensitivity analysis against the itraconazole
    # and rifampin interaction studies (Methods section 2.4; Figure S2
    # shows it is insensitive to the CL_R and CL_bile assumptions).
    # This is the arm a CYP3A4 inhibitor or inducer acts on.

    lcl_other <- fixed(log(0.894))
    label("Non-CYP3A4, non-renal systemic plasma clearance cl_other (L/h)")
    # Derived from Tian 2025 Figure 2 as the sum of the two remaining
    # hepatic routes: non-CYP3A4 hepatic metabolism, which is
    # CL_hepatic metabolism - cl_3a4 = 1.33 - 0.676 = 0.654 L/h
    # (UGT1A8 / UGT1A9 glucuronidation and amide hydrolysis, f_m,other
    # = 0.39), plus biliary excretion of unchanged drug,
    # CL_bile = 0.24 L/h (f_e,bile = 0.14). 0.654 + 0.240 = 0.894 L/h.
    # Kept as one parameter because neither route is a CYP3A4
    # drug-interaction target: the paper shows in the Discussion that
    # even complete UGT inhibition would raise exposure less than
    # 1.5-fold and that no clinical UGT1A8 / UGT1A9 inhibitor exists.

    lcl_renal <- fixed(log(0.12))
    label("Renal clearance cl_renal (L/h)")
    # Tian 2025 Table 1 and Figure 2: CL_R = 0.12 L/h, "measured from
    # intravenous administration of pirtobrutinib". Section 3.1 notes
    # that an oral solution gave 0.29 L/h instead, and that the choice
    # does not change f_m,CYP3A4 (Figure S2); the intravenous value is
    # the one the model used.

    # Systemic-compartment volume. Data S3 reports the two SAC first-order
    # rate constants together with the inter-compartmental clearance:
    #   SAC kin  = 0.12402043422142926 1/h
    #   SAC kout = 0.2594761264470896  1/h
    #   SAC Q    = 3.56 L/h
    # In the Simcyp minimal-PBPK layout kin = Q / V_systemic and
    # kout = Q / V_SAC, so
    #   vc      = 3.56 / 0.1240204 = 28.7049 L
    #   V_SAC   = 3.56 / 0.2594761 = 13.7200 L
    # and 13.7200 L / (0.17 L/kg) = 80.7 kg recovers the body weight of
    # the Simcyp population representative, which confirms the reading.
    # No reference body weight has to be assumed: vc is absolute.
    lvc <- fixed(log(28.7049))
    label("Systemic compartment volume vc (L)")
    # Derived from Data S3 rows "SAC Q (L/h)" = 3.56 and "SAC kin (1/h)"
    # = 0.12402043422142926 as Q / kin.

    # The Simcyp kin and kout are first-order rate constants acting on
    # the MASSES of drug in the systemic compartment and in the single
    # adjusting compartment respectively. That is exactly the definition
    # of the canonical central-to-peripheral1 and peripheral1-to-central
    # micro-constants, so they are encoded as k12 / k21. The SAC volume
    # of 0.17 L/kg does not enter the plasma prediction under a
    # mass-based parameterisation and is therefore not carried; the
    # peripheral volume implied by vc * k12 / k21 does equal it here
    # (28.7049 * 0.1240204 / 0.2594761 = 13.72 L), because both rate
    # constants share the same Q.
    lk12 <- fixed(log(0.1240204))
    label("Systemic-to-SAC transfer rate constant k12 (Simcyp kin, 1/h)")
    # Data S3 row "SAC kin (1/h)" = 0.12402043422142926. Table 1 records
    # that Q and V_SAC were "estimated from intravenous geometric mean
    # concentration-time profile using the Parameter Estimation module
    # within Simcyp".

    lk21 <- fixed(log(0.2594761))
    label("SAC-to-systemic transfer rate constant k21 (Simcyp kout, 1/h)")
    # Data S3 row "SAC kout (1/h)" = 0.2594761264470896, equivalently
    # Q / V_SAC = 3.56 / (0.17 L/kg * 80.7 kg).

    # Tian 2025 is a PBPK simulation analysis, not a population-PK fit.
    # It reports no residual-error model and no estimated
    # inter-individual variance components; the percent-CV figures in
    # Table 2 are the spread of a Simcyp virtual population driven by
    # unpublished population files, not estimated omegas. Data S3 does
    # print 30% input CVs on ka, lag time, Q_gut, the additional HLM
    # CLint and the biliary CLint, but total-clearance variability in
    # Simcyp is dominated by liver weight and CYP3A4 abundance, which
    # are population-database outputs. Rather than invent a variance,
    # the residual error is fixed at zero, which makes this a
    # deterministic typical-value simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    # ------------------------------------------------------------------
    ka        <- exp(lka)
    tlag      <- exp(ltlag)
    cl_3a4    <- exp(lcl_3a4)
    cl_other  <- exp(lcl_other)
    cl_renal  <- exp(lcl_renal)
    vc        <- exp(lvc)
    k12       <- exp(lk12)
    k21       <- exp(lk21)

    # ------------------------------------------------------------------
    # 2. Clearance assembly, Tian 2025 Equation 1 and Figure 2.
    # Hepatic clearance is everything except the renal route, and total
    # systemic clearance adds the renal route back:
    #   cl_h = 0.676 + 0.894 = 1.570 L/h
    #   cl   = 1.570 + 0.120 = 1.690 L/h
    # Scaling lcl_3a4 up or down therefore moves both the elimination
    # rate and, through fh below, the oral bioavailability, which is how
    # a CYP3A4 inhibitor or inducer acts in this model.
    # ------------------------------------------------------------------
    cl_h <- cl_3a4 + cl_other
    cl   <- cl_h + cl_renal
    kel  <- cl / vc

    # ------------------------------------------------------------------
    # 3. Oral bioavailability, Tian 2025 Equations 2 and 6.
    #
    # Equation 2, the well-stirred hepatic availability, with the
    # measured blood-to-plasma ratio of 0.79 and hepatic blood flow of
    # 90 L/h, so (B:P) * Q_H = 71.1 L/h:
    #   fh = 1 - cl_h / 71.1 = 1 - 1.570 / 71.1 = 0.977915
    # Figure 2 prints Fh = 0.98, which matches.
    #
    # Equation 6 defines the fraction absorbed from the measured
    # absolute oral bioavailability F = 0.86 (Figure 2, absolute
    # bioavailability study) and the intestinal availability Fg = 0.96
    # (Equation 5, from the 4% pirtobrutinib Cmax rise under
    # itraconazole and the remaining intestinal CYP3A4 activity
    # X = 0.009 of Table S3):
    #   fa = F / (Fg * fh) = 0.86 / (0.96 * 0.977915) = 0.91606
    # Table 1 prints this rounded to 0.91. The product fa * Fg is
    # therefore 0.91606 * 0.96 = 0.879419, and it is carried as that
    # single constant because neither factor changes when hepatic CYP3A4
    # activity changes. Multiplying by the derived fh recovers
    # F = 0.879419 * 0.977915 = 0.86000 at baseline.
    # ------------------------------------------------------------------
    fh     <- 1 - cl_h / 71.1
    fdepot <- 0.879419 * fh

    # ------------------------------------------------------------------
    # 4. ODE system, corresponding to Tian 2025 Figure 2 once the liver
    # and portal-vein compartments of the Simcyp minimal-PBPK layout are
    # lumped into the systemic compartment. `peripheral1` is the single
    # adjusting compartment: k12 acts on the mass in the systemic
    # compartment and moves drug into the SAC, k21 acts on the mass in
    # the SAC and returns it. Amounts are in mg.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Bioavailability and lag time on the oral depot.
    # ------------------------------------------------------------------
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    # throughout Tian 2025 Table 2 and Tables S10-S12.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
