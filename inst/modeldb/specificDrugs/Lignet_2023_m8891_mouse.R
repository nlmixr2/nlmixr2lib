Lignet_2023_m8891_mouse <- function() {
  description <- "Preclinical (mouse, Caki-1 renal-carcinoma xenograft). One-compartment oral PK of M8891, a selective and reversible methionine aminopeptidase 2 (MetAP2) inhibitor, with first-order absorption and elimination, linked through an effect compartment to a turnover model for the tumour target-engagement biomarker Met-EF1a (uncleaved methionine-elongation-factor-1-alpha). M8891 inhibits the first-order degradation of Met-EF1a, so the biomarker accumulates above its kin/kout baseline. Naive-pooled fit (Phoenix WinNonlin 6.4) to single-dose and 4-day repeated-dose PK/PD data at 10, 25, and 100 mg/kg p.o.; the authors report no inter-individual variability because a nonlinear mixed-effects fit did not converge on this dataset."
  reference <- paste(
    "Lignet F, Friese-Hamim M, Jaehrling F, El Bawab S, Rohdich F.",
    "Preclinical Pharmacokinetics and Translational",
    "Pharmacokinetic/Pharmacodynamic Modeling of M8891, a Potent and",
    "Reversible Inhibitor of Methionine Aminopeptidase 2.",
    "Pharm Res. 2023;40(12):3011-3023.",
    "doi:10.1007/s11095-023-03611-z.",
    sep = " "
  )
  vignette <- "Lignet_2023_m8891"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # Met-EF1a is the paper's target-engagement biomarker (uncleaved
  # methionine-elongation-factor-1-alpha in tumour tissue, ug per mg total
  # protein). It is a genuinely paper-mechanistic state with no canonical
  # analogue in inst/references/compartment-names.md, so it is declared here
  # rather than forced onto a canonical compartment name. The canonical
  # `effect` compartment is already used for the paper's effect compartment
  # Ce (Lignet 2023 Eq. 5 / Fig. 1).
  paper_specific_compartments <- c("metef1a")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "M8891", units = NA_character_, specimen = "administration site", verified = FALSE),
    central = list(analyte = "M8891", units = NA_character_, specimen = "plasma", verified = FALSE),
    effect  = list(analyte = "Met-EF1a", units = NA_character_, specimen = "not applicable", verified = FALSE),
    metef1a = list(analyte = "Met-EF1a", units = NA_character_, specimen = "tumor", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (female CD1 nu/nu, subcutaneous Caki-1 human renal-carcinoma xenograft)",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "5-6 weeks old at tumour-cell inoculation",
    weight_range   = NA_character_,
    sex_female_pct = 100,
    disease_state  = "Subcutaneous Caki-1 (ATCC HTB-46) human renal-cell-carcinoma xenograft in the right flank; treatment started when tumours reached 300-500 mm3.",
    dose_range     = "M8891 10, 25, or 100 mg/kg p.o. in 0.25% Methocel: a single administration (study PK/PD-01) or four daily administrations (study PK/PD-02). Plasma and tumour tissue were collected 1, 7, 24, 48, 72, and 96 h after the (last) dose.",
    regions        = NA_character_,
    notes          = "Lignet 2023 Methods, 'PK/PD Studies' and 'PK/PD Modeling of Caki-1 Xenograft Data'. The paper states 'mice were randomized into treatment groups (n = 5)' for each of the three dose levels in each of the two PK/PD studies, but does not tabulate a total animal count; because animals were euthanised at each sampling time (1, 7, 24, 48, 72, 96 h) the design is destructive/serial-sacrifice rather than serial-sampling, and it is not stated whether n = 5 is per group or per group-and-timepoint. n_subjects is therefore left NA rather than guessed. Modelling proceeded in three steps: (1) pool PK/PD-01 and PK/PD-02 plasma data to fit the PK model; (2) pool the Met-EF1a data and fit the PD model driven by PK-model-simulated plasma concentrations; (3) validate by simulating Met-EF1a at the doses of efficacy study EFF-01 (Fig. 3e). No time-dependent PK was observed between day 1 and day 4. Parameter estimates and their estimation CV% are Lignet 2023 Table IV."
  )

  ini({
    # --- PK, Lignet 2023 Table IV (mouse Caki-1 xenograft PK/PD model). ---
    # CL/F and V/F are apparent oral parameters normalised to body weight
    # (L/h/kg and L/kg); doses are given in mg/kg, so Cc is in mg/L
    # (= ug/mL = 1000 ng/mL, the unit used in Lignet 2023 Fig. 3a-b).
    lka <- fixed(log(2.7)); label("Absorption rate constant ka (1/h)")                            # Lignet 2023 Table IV: ka = 2.7 1/h, reported as "Fixed" in the "Estimate CV (%)" column
    lcl <- log(0.415); label("Apparent oral clearance CL/F (L/h/kg)")                             # Lignet 2023 Table IV: CL/F = 0.415 L/h/kg (estimation CV 13.8%)
    lvc <- log(1.034); label("Apparent oral volume of distribution V/F (L/kg)")                   # Lignet 2023 Table IV: V/F = 1.034 L/kg (estimation CV 15.7%)

    # --- PD, Lignet 2023 Table IV and Eqs. (5)-(6). ---
    lke0 <- log(0.0566); label("Effect-compartment equilibration rate constant ke0 (1/h)")        # Lignet 2023 Table IV: ke0 = 0.0566 1/h (estimation CV 7.3%)
    lkin <- log(29.1); label("Zero-order Met-EF1a synthesis rate kin (ug/mg protein/h)")          # Lignet 2023 Table IV: Kin = 29.1 ug/mg/h (estimation CV 22%)
    lkout <- log(1.45); label("First-order Met-EF1a degradation rate constant kout (1/h)")        # Lignet 2023 Table IV: Kout = 1.45 1/h (estimation CV 25.7%)
    limax <- log(0.91); label("Maximum fractional inhibition of Met-EF1a degradation Imax (unitless)") # Lignet 2023 Table IV: Imax = 0.91 (estimation CV 1.3%)
    lic50 <- log(0.340); label("Effect-site concentration producing half-maximal inhibition IC50 (mg/L)") # Lignet 2023 Table IV: IC50 = 340 ng/mL = 0.340 mg/L (estimation CV 17%); Results also quote 0.88 uM total / 28 nM free

    # --- Residual error. ---
    # Lignet 2023 Results ("PK/PD Modeling") states that the PK data were
    # fitted "with a multiplicative error model", i.e. proportional error in
    # nlmixr2's linear space, but the paper does not report its magnitude
    # anywhere (Table IV lists only the estimation CV% of the structural
    # parameters). No error model is described at all for the Met-EF1a data.
    # Both are therefore encoded with the correct FORM and a magnitude fixed
    # at zero, so the model reproduces the published typical-value curves
    # exactly; see the vignette's "Assumptions and deviations" section.
    propSd <- fixed(0); label("Proportional residual error on plasma concentration (fraction)")   # Lignet 2023 Results: "with a multiplicative error model"; magnitude not reported
    propSd_MetEF1a <- fixed(0); label("Proportional residual error on Met-EF1a (fraction)")       # not reported in Lignet 2023; form assumed to match the PK error model
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ke0 <- exp(lke0)
    kin <- exp(lkin)
    kout <- exp(lkout)
    imax <- exp(limax)
    ic50 <- exp(lic50)

    kel <- cl / vc

    # One-compartment oral PK with first-order absorption and elimination.
    # Lignet 2023 Eq. (3) gives the closed form
    #   C(t) = ka * Dose / ((ka - ke) * V/F) * (exp(-ke*t) - exp(-ka*t))
    # and Eq. (4) gives ke = CL/V. F is folded into the apparent CL/F and
    # V/F, so no separate f(depot) term is applied (Eq. (3) does not
    # multiply Dose by F).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Effect compartment, Lignet 2023 Eq. (5): dCe/dt = ke0 * (C - Ce).
    # The paper's Fig. 1 shows the same rate constant ke0 governing transfer
    # into and out of the effect compartment, so Ce equilibrates towards the
    # plasma concentration with no gain factor.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Met-EF1a turnover with inhibition of degradation, Lignet 2023 Eq. (6):
    #   dE/dt = kin - kout * E * (1 - Imax * Ce / (Ce + IC50))
    # The biomarker starts from its undrugged steady state kin/kout
    # (29.1 / 1.45 = 20.1 ug/mg protein, matching the pre-dose and 96-h
    # levels in Lignet 2023 Fig. 3c-d).
    d/dt(metef1a) <- kin - kout * metef1a * (1 - imax * effect / (effect + ic50))
    metef1a(0) <- kin / kout

    MetEF1a <- metef1a

    Cc ~ prop(propSd)
    MetEF1a ~ prop(propSd_MetEF1a)
  })
}
