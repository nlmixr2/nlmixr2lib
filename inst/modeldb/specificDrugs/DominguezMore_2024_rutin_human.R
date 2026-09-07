DominguezMore_2024_rutin_human <- function() {
  description <- paste(
    "One-compartment intravenous human PK model for rutin",
    "(quercetin-3-O-rutinoside) obtained by simple interspecies allometry from",
    "the rat and rabbit population PK fits of Dominguez More 2024 (Table 5).",
    "It is a forward projection, not a fit to human data: V and Cl are the",
    "only parameters the allometry supplies, so the two-compartment",
    "distribution of the animal models cannot be carried over, and no",
    "between-subject variability or residual error is reported. The",
    "P. peruviana calyx extract matrix is retained as a categorical covariate,",
    "the paper's central finding being that the increased volume of",
    "distribution and clearance seen with the extract in animals is preserved",
    "in the human predictions."
  )
  reference <- paste(
    "Dominguez More GP, Rey DP, Valderrama IH, Ospina LF, Aragon DM.",
    "Rutin and Physalis peruviana extract: population pharmacokinetics in",
    "New Zealand rabbits. Pharmaceutics. 2024;16(10):1241.",
    "doi:10.3390/pharmaceutics16101241.",
    "The rat parameters entering the allometry come from Dominguez-More GP,",
    "Sepulveda PM, Echeverry SM, Oliveira-Simoes CM, Aragon DM. Matrix effects",
    "of the hydroethanolic extract of calyces of Physalis peruviana L. on",
    "rutin pharmacokinetics in Wistar rats using population modeling.",
    "Pharmaceutics. 2021;13(4):535. doi:10.3390/pharmaceutics13040535",
    sep = " "
  )
  vignette <- "DominguezMore_2024_rutin_physalis"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "rutin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg -- the human weight Dominguez More 2024",
        "Section 2.3.3 states the allometry was evaluated at ('W is the human",
        "weight, typically set at 70 kg'). Both V and Cl scale as a power of",
        "body weight with an arm-specific exponent read from Table 5, so a",
        "prediction at a weight other than 70 kg diverges between the pure",
        "rutin and extract arms by more than the fixed 70 kg ratio."
      ),
      source_name        = "W"
    ),
    FORM_RUTIN_EXTRACT = list(
      description        = "Source of the administered rutin: within the Physalis peruviana calyx extract matrix versus the isolated pure compound",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pure rutin)",
      notes              = paste(
        "Selects the RUT (pure rutin) or EXT (rutin within the extract) row",
        "of Dominguez More 2024 Table 5. Because the allometric coefficient a",
        "AND the exponent b both differ between the arms, the covariate acts",
        "on the typical value and on the weight exponent, not on the typical",
        "value alone."
      ),
      source_name        = "source of rutin"
    )
  )

  population <- list(
    species        = "human (projected by allometry from rat and rabbit)",
    n_subjects     = 0L,
    n_studies      = 0L,
    weight_range   = "70 kg reference (Section 2.3.3)",
    disease_state  = "not applicable -- forward projection, no human subjects were studied",
    dose_range     = "not applicable -- no human dose was administered or proposed",
    regions        = "not applicable",
    notes          = paste(
      "Simple allometry, Y = a * W^b (Eq. 6), fitted across two species: the",
      "Wistar rat individual parameters of Dominguez-More 2021 and the New",
      "Zealand White rabbit individual parameters of the present paper. Table",
      "5 gives, per arm, the animal parameters, the coefficient a, the",
      "exponent b, and the 70 kg human prediction using both the fitted",
      "('experimental') b and the theoretical b (1.0 for V, 0.75 for Cl):",
      "",
      "  V,  pure rutin: rat 0.024 L,    rabbit 0.096 L,    a 0.057, b 0.8;",
      "                  human 1.410 L (experimental b), 3.988 L (theoretical b)",
      "  V,  extract:    rat 0.035 L,    rabbit 0.190 L,    a 0.102, b 0.9;",
      "                  human 4.592 L (experimental b), 7.128 L (theoretical b)",
      "  Cl, pure rutin: rat 0.031 L/h,  rabbit 0.188 L/h,  a 0.096, b 0.9;",
      "                  human 5.389 L/h (experimental b), 2.323 L/h (theoretical b)",
      "  Cl, extract:    rat 0.072 L/h,  rabbit 0.690 L/h,  a 0.296, b 1.2;",
      "                  human 48.938 L/h (experimental b), 7.160 L/h (theoretical b)",
      "",
      "This model file carries the EXPERIMENTAL-b column, which is the",
      "paper's own data-derived allometric relationship. The theoretical-b",
      "column is reproduced exactly by a * 70^b with b = 1.0 for V and 0.75",
      "for Cl (0.057 * 70 = 3.99; 0.102 * 70 = 7.14; 0.096 * 70^0.75 = 2.32;",
      "0.296 * 70^0.75 = 7.16), so a user who prefers the theoretical",
      "exponents can substitute those four numbers directly. Note that the",
      "printed b values are rounded to one decimal and do not exactly",
      "regenerate the experimental-b column from a: 0.057 * 70^0.8 = 1.71",
      "against the tabulated 1.410 for V of pure rutin (the unrounded",
      "exponent is 0.755); the other three rows agree within about 5%. The",
      "typical values below are therefore taken from the tabulated human",
      "column rather than recomputed, so the model returns Table 5 exactly at",
      "70 kg.",
      "",
      "Only two species entered the regression, so each exponent is",
      "determined by a single pair of points and carries no uncertainty",
      "estimate; the paper itself flags the extract clearance exponent of 1.2",
      "as positive allometry. The absence of a third species, of human data,",
      "and of any variability term means this model is a scaling prior and",
      "not a population model."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Human typical values at the 70 kg reference weight -- Table 5, "Human
    # Parameter (Experimental b)" column. Every value is fixed: none was
    # estimated from human data, they are arithmetic consequences of the
    # animal fits and Eq. 6.
    #
    # The two treatment arms carry SEPARATE allometric coefficients AND
    # exponents in Table 5, so both the typical value and the weight exponent
    # take a stratum suffix naming the arm; neither arm keeps the bare
    # canonical (see inst/references/parameter-names.md, stratum-suffixed
    # parameters).
    # ----------------------------------------------------------------------
    lvc_purerutin <- fixed(log(1.410))  ; label("Human volume of distribution at 70 kg, pure rutin (L)")        # Table 5: V RUT human, experimental b
    lvc_extract   <- fixed(log(4.592))  ; label("Human volume of distribution at 70 kg, rutin in extract (L)")  # Table 5: V EXT human, experimental b
    lcl_purerutin <- fixed(log(5.389))  ; label("Human clearance at 70 kg, pure rutin (L/h)")                   # Table 5: Cl RUT human, experimental b
    lcl_extract   <- fixed(log(48.938)) ; label("Human clearance at 70 kg, rutin in extract (L/h)")             # Table 5: Cl EXT human, experimental b

    # Allometric exponents b of Eq. 6, per arm and per parameter (Table 5,
    # "b Exponent" column). Fixed: read off a two-species regression, with no
    # reported uncertainty.
    e_wt_vc_purerutin <- fixed(0.8) ; label("Allometric weight exponent on V, pure rutin (unitless)")           # Table 5: b = 0.8 for V RUT
    e_wt_vc_extract   <- fixed(0.9) ; label("Allometric weight exponent on V, rutin in extract (unitless)")     # Table 5: b = 0.9 for V EXT
    e_wt_cl_purerutin <- fixed(0.9) ; label("Allometric weight exponent on Cl, pure rutin (unitless)")          # Table 5: b = 0.9 for Cl RUT
    e_wt_cl_extract   <- fixed(1.2) ; label("Allometric weight exponent on Cl, rutin in extract (unitless)")    # Table 5: b = 1.2 for Cl EXT

    # No between-subject variability and no residual error are reported for
    # the allometric projection -- it is a deterministic scaling of typical
    # values -- so the residual SD is fixed at zero rather than invented.
    propSd <- fixed(0) ; label("Proportional residual error (fraction; not reported)")                          # Dominguez More 2024 reports no residual-error model for the allometric projection
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Arm selection. FORM_RUTIN_EXTRACT picks the RUT or the EXT row of
    #    Table 5 for both the 70 kg typical value and the weight exponent.
    # ----------------------------------------------------------------------
    lvc70   <- lvc_purerutin     * (1 - FORM_RUTIN_EXTRACT) + lvc_extract     * FORM_RUTIN_EXTRACT
    lcl70   <- lcl_purerutin     * (1 - FORM_RUTIN_EXTRACT) + lcl_extract     * FORM_RUTIN_EXTRACT
    e_wt_vc <- e_wt_vc_purerutin * (1 - FORM_RUTIN_EXTRACT) + e_wt_vc_extract * FORM_RUTIN_EXTRACT
    e_wt_cl <- e_wt_cl_purerutin * (1 - FORM_RUTIN_EXTRACT) + e_wt_cl_extract * FORM_RUTIN_EXTRACT

    # ----------------------------------------------------------------------
    # 2. Allometry, Eq. 6 (Y = a * W^b), re-anchored at the 70 kg reference
    #    so that WT = 70 returns Table 5's human column exactly.
    # ----------------------------------------------------------------------
    vc <- exp(lvc70) * (WT / 70)^e_wt_vc
    cl <- exp(lcl70) * (WT / 70)^e_wt_cl

    kel <- cl / vc

    # ----------------------------------------------------------------------
    # 3. One-compartment intravenous disposition. The allometry supplies only
    #    V and Cl, so the peripheral compartment of the animal models is not
    #    carried forward.
    # ----------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ----------------------------------------------------------------------
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is in
    #    mg/L; the factor 1000 converts to ng/mL.
    # ----------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
