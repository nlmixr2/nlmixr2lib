An_2012_mitoxantrone_human_pbpk <- function() {
  description <- paste(
    "Human-scaled simulation. Semi-mechanistic whole-body PBPK model 3",
    "for mitoxantrone (Novantrone) in adult cancer patients after a",
    "single 12 mg/m^2 IV bolus, projected forward from the mouse fit in",
    "An_2012_mitoxantrone_mouse_pbpk (An and Morris 2012, AAPS J). Same",
    "topology as the mouse model: seven physiological tissue compartments",
    "(central plasma plus six perfusion-limited well-stirred organs - lung,",
    "heart, spleen, liver, kidney, brain) and a permeability-limited",
    "remainder compartment lumping muscle, fat, bone, and skin and",
    "resolving into an interstitial (is_remainder) and an intracellular",
    "(int_remainder) subspace coupled by a permeability-surface area",
    "product. Plasma unbound fraction fu = 0.2, DNA dissociation constant",
    "K_DNA = 0.0013 uM, protein-binding dissociation constant K_macro =",
    "1.44 uM, and per-organ T_macro are carried over from the mouse fit",
    "as cross-mammalian constants (paper Methods). Per-organ T_DNA values",
    "are replaced with the human DNA-content literature values reported",
    "in Table III: lung and spleen DNA use the literature rapidly-",
    "perfused-organ value (15 uM); brain T_DNA is the mouse-derived value",
    "(0.10 uM); the remainder T_DNA uses the literature slowly-perfused-",
    "organ value (4.5 uM). The remainder ISF/intracellular split is",
    "assumed to follow the mouse proportion 33/67 (Table I) applied to",
    "the human remainder volume V_other = 62 L. Hepatic and renal",
    "intrinsic clearances are derived from clinical CL_H = 19 L/h/m^2 and",
    "CL_R = 2.7 L/h/m^2 (Ehninger 1990 ref 6) via the well-stirred",
    "rearrangement Clint = Q_H * CL / (fu * Q_H - fu * CL) yielding",
    "Clint_H = 250 L/h and Clint_R = 27 L/h. PS_remainder is",
    "allometrically scaled from the mouse value PS = 1.44 mL/min via",
    "PS = A * M^0.75 (Kawai 1998 ref 28) giving PS = 31.1 L/h. The model",
    "is a typical-value forward simulation; there is no IIV and no",
    "residual error from the paper."
  )
  reference <- paste(
    "An G, Morris ME. A Physiologically Based Pharmacokinetic Model of",
    "Mitoxantrone in Mice and Scale-up to Humans: a Semi-Mechanistic",
    "Model Incorporating DNA and Protein Binding. AAPS J.",
    "2012;14(2):352-364. doi:10.1208/s12248-012-9344-7.",
    "PK structure and binding constants K_DNA / K_macro / T_macro",
    "inherited from the mouse fit; see modellib('An_2012_mitoxantrone_mouse_pbpk').",
    sep = " "
  )
  vignette <- "An_2012_mitoxantrone_pbpk"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # BSA is used externally as the IV-dose-mass multiplier (dose mg =
  # 12 mg/m^2 * BSA), but does NOT appear in model() because the An
  # and Morris 2012 human simulation fixes Q_i and V_i at the 70-kg
  # Table III reference values. Documented in covariatesDataExcluded
  # to preserve the dose-conversion provenance while passing the
  # checkModelConventions() unused-covariate check.
  covariatesDataExcluded <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "External dose-mass multiplier only (dose mg = 12 mg/m^2 *",
        "BSA). The An and Morris 2012 human PBPK simulation fixes",
        "Q_i and V_i at the 70-kg Davies and Morris 1993 (ref 25)",
        "reference physiology in Table III; BSA does not enter the",
        "ODE system. Users simulating across patients of different",
        "body sizes should adjust the physiological constants in",
        "their own simulation wrapper."
      ),
      source_name        = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "adult (Larson 1987 J Clin Oncol acute non-lymphocytic leukemia; Peng 1982 J Chromatogr breast cancer)",
    weight_range   = "Davies and Morris (1993) 70-kg reference subject for the physiology",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Adult cancer patients receiving 12 mg/m^2 IV mitoxantrone in two",
      "previously published clinical pharmacokinetic studies (Larson",
      "1987 J Clin Oncol 5:391-7 acute non-lymphocytic leukemia; Peng",
      "1982 J Chromatogr 233:235-47 advanced breast cancer). The An",
      "and Morris 2012 PBPK simulation used these two studies' plasma",
      "concentration vs. time observations as the human validation",
      "target (Fig 6)."
    ),
    dose_range     = "12 mg/m^2 single IV bolus mitoxantrone (typical clinical dose)",
    regions        = NA_character_,
    notes          = paste(
      "Pure forward simulation: no individual-level data were fit. The",
      "An and Morris 2012 paper reports comparisons of the model-",
      "predicted plasma concentration vs. mean digitised observations",
      "from Larson 1987 and Peng 1982 (Fig 6, gray and black symbols",
      "respectively). The model is intended for typical-value",
      "simulations of human plasma and tissue mitoxantrone disposition;",
      "tissue concentrations are not reported in the source clinical",
      "studies."
    )
  )

  ini({
    # All parameters in the human-scaled model are FIXED. Either inherited
    # from the An and Morris 2012 mouse fit (K_DNA, K_macro, per-organ
    # T_macro) per the paper's cross-mammalian assumption, taken from
    # Table III human-physiology references (per-organ T_DNA), or derived
    # by paper Eq 14 / Eq 13 from clinical CL_H, CL_R, and the mouse
    # PS (Clint_H, Clint_R, PS_remainder).

    # DNA binding capacities (uM) -- Table III values. Lung and spleen
    # use the human DNA value for rapidly perfused organs; remainder
    # uses the human DNA value for slowly perfused organs; brain uses
    # the mouse-fit value (no human brain DNA value available).
    lt_dna_lung   <- fixed(log(15));     label("Tissue DNA binding capacity, lung (uM)")            # An 2012 Table III (15, rapidly-perfused human DNA)
    lt_dna_heart  <- fixed(log(8.3));    label("Tissue DNA binding capacity, heart (uM)")           # An 2012 Table III (8.3, Gustafson 2002 ref 17)
    lt_dna_spleen <- fixed(log(15));     label("Tissue DNA binding capacity, spleen (uM)")          # An 2012 Table III (15, rapidly-perfused human DNA)
    lt_dna_liver  <- fixed(log(23.7));   label("Tissue DNA binding capacity, liver (uM)")           # An 2012 Table III (23.7, Gustafson 2002 ref 17)
    lt_dna_kidney <- fixed(log(16.2));   label("Tissue DNA binding capacity, kidney (uM)")          # An 2012 Table III (16.2, Gustafson 2002 ref 17)
    lt_dna_brain  <- fixed(log(0.1));    label("Tissue DNA binding capacity, brain (uM)")           # An 2012 Table III (0.1, mouse fit; no human brain value available)
    lt_dna_rem    <- fixed(log(4.5));    label("Tissue DNA binding capacity, remainder (uM)")       # An 2012 Table III (4.5, slowly-perfused human DNA)

    # Protein (macromolecule) binding capacities (uM) -- Table III; per
    # paper text, T_macro values are "assumed to be identical in
    # mammals" and use the An and Morris 2012 mouse Model 3 fit.
    lt_macro_lung   <- fixed(log(19));   label("Tissue protein binding capacity, lung (uM)")        # An 2012 Table III / mouse fit (19)
    lt_macro_heart  <- fixed(log(43.8)); label("Tissue protein binding capacity, heart (uM)")       # An 2012 Table III / mouse fit (43.8)
    lt_macro_spleen <- fixed(log(3.54)); label("Tissue protein binding capacity, spleen (uM)")      # An 2012 Table III / mouse fit (3.54)
    lt_macro_liver  <- fixed(log(514));  label("Tissue protein binding capacity, liver (uM)")       # An 2012 Table III / mouse fit (514)
    lt_macro_kidney <- fixed(log(52.3)); label("Tissue protein binding capacity, kidney (uM)")      # An 2012 Table III / mouse fit (52.3)
    lt_macro_rem    <- fixed(log(4.67)); label("Tissue protein binding capacity, remainder (uM)")   # An 2012 Table III / mouse fit (4.67)

    # Binding affinities -- inherited from the An and Morris 2012 mouse
    # Model 3 fit (paper: "assumed to be identical in mammals").
    lk_dna   <- fixed(log(0.0013));      label("DNA binding dissociation constant K_DNA (uM)")      # An 2012 Table II Model 3 mouse fit (0.0013)
    lk_macro <- fixed(log(1.44));        label("Protein binding dissociation constant K_macro (uM)") # An 2012 Table II Model 3 mouse fit (1.44)

    # Intrinsic clearances -- derived from clinical CL_H = 19 L/h/m^2 and
    # CL_R = 2.7 L/h/m^2 (Ehninger 1990 ref 6) via paper Eq 14:
    # Clint = Q_H * CL / (fu * Q_H - fu * CL); the paper Methods
    # section states the resulting values are Clint_H = 250 L/h and
    # Clint_R = 27 L/h for the human simulation.
    lcl_int_h <- fixed(log(250));        label("Hepatic intrinsic clearance Clint_H (L/h)")         # An 2012 Methods text ("Clint,H and Clint,R used in human simulation were 250 L/h and 27 L/h")
    lcl_int_r <- fixed(log(27));         label("Renal intrinsic clearance Clint_R (L/h)")           # An 2012 Methods text (same sentence)

    # Permeability-surface area product -- allometrically scaled from the
    # mouse PS = 1.44 mL/min via paper Eq 13 PS = A * M^B with B = 0.75
    # (fixed) and yielded PS_human = 31.1 L/h.
    lps_rem <- fixed(log(31.1));         label("Permeability-surface area product, remainder PS (L/h)") # An 2012 Methods text ("PSt used in the human simulation was 31.1 L/h")

    # Residual error -- the An and Morris 2012 human simulation is a
    # forward typical-value projection only; no residual error
    # estimate was fit to the human data. propSd is a fixed placeholder
    # for syntactic completeness; see vignette Assumptions and
    # deviations.
    propSd <- fixed(0.10);               label("Proportional residual error placeholder (fraction)") # not estimated for human simulation; placeholder
  })

  model({
    # ===== Drug physicochemistry =====
    mw_mitox <- 444.49   # mitoxantrone molecular weight (g/mol)

    # Plasma unbound fraction fixed to 0.2 (paper Methods: "fu was
    # fixed to 0.2" for both the mouse fit and the human simulation).
    fu <- 0.2

    # ===== Human physiology (Table III, 70-kg reference) =====
    # Organ volumes (L) directly from Table III (Davies and Morris 1993
    # 70-kg adult reference, ref 25).
    vp        <- 3.0       # plasma
    v_lung    <- 1.17
    v_heart   <- 0.31
    v_spleen  <- 0.192
    v_liver   <- 1.69
    v_kidney  <- 0.28
    v_brain   <- 1.45
    v_rem_tot <- 62.0      # remainder (muscle + fat + bone + skin lumped)

    # ISF / intracellular split for the remainder. Table III gives only
    # the total volume V_other = 62 L; the An and Morris 2012 paper
    # does not list a separate ISF / intracellular split for human.
    # Applied here the mouse Table I 33/67 split (V_isf / V_total =
    # 8.26 / 25.03 = 0.330, V_int / V_total = 16.77 / 25.03 = 0.670),
    # which matches standard human ECF / ICF physiology (ECF ~ 14 L,
    # ICF ~ 28 L per 70-kg subject; ratio ~ 1/2 -> 33% extracellular).
    v_isf_rem <- v_rem_tot * (8.26 / 25.03)    # = 20.46 L
    v_int_rem <- v_rem_tot * (16.77 / 25.03)   # = 41.54 L

    # Organ blood flow rates -- Table III lists Q in L/min; convert to
    # L/h by multiplying by 60.
    q_lung   <- 5.6   * 60   # = 336 L/h (cardiac output)
    q_heart  <- 0.24  * 60   # = 14.4 L/h
    q_spleen <- 0.077 * 60   # = 4.62 L/h
    q_liver  <- 1.45  * 60   # = 87 L/h
    q_kidney <- 1.24  * 60   # = 74.4 L/h
    q_brain  <- 0.7   * 60   # = 42 L/h
    q_rem    <- 1.893 * 60   # = 113.58 L/h

    # ===== Individual structural parameters (all fixed) =====
    t_dna_lung   <- exp(lt_dna_lung)
    t_dna_heart  <- exp(lt_dna_heart)
    t_dna_spleen <- exp(lt_dna_spleen)
    t_dna_liver  <- exp(lt_dna_liver)
    t_dna_kidney <- exp(lt_dna_kidney)
    t_dna_brain  <- exp(lt_dna_brain)
    t_dna_rem    <- exp(lt_dna_rem)

    t_macro_lung   <- exp(lt_macro_lung)
    t_macro_heart  <- exp(lt_macro_heart)
    t_macro_spleen <- exp(lt_macro_spleen)
    t_macro_liver  <- exp(lt_macro_liver)
    t_macro_kidney <- exp(lt_macro_kidney)
    t_macro_rem    <- exp(lt_macro_rem)

    # Brain T_macro is not reported in Table III ("Brain T_macro = -")
    # and was not used in the An and Morris 2012 human simulation; the
    # brain Kp_eff therefore reduces to the DNA-only form.
    t_macro_brain <- 0

    k_dna   <- exp(lk_dna)
    k_macro <- exp(lk_macro)

    cl_int_h <- exp(lcl_int_h)
    cl_int_r <- exp(lcl_int_r)
    ps_rem   <- exp(lps_rem)

    # ===== Concentrations (mg/L) and plasma uM conversion =====
    c_plasma <- central / vp
    cp_uM    <- c_plasma * 1000 / mw_mitox

    # ===== Cp-dependent effective tissue:plasma partition coefficients =====
    # An and Morris 2012 Eq 9 inline definition; see
    # An_2012_mitoxantrone_mouse_pbpk for the rationale. Same
    # functional form, just evaluated with the human T_DNA and T_macro
    # values above.
    denom_dna   <- k_dna   / fu + cp_uM
    denom_macro <- k_macro / fu + cp_uM

    kp_lung   <- t_dna_lung   / denom_dna + t_macro_lung   / denom_macro
    kp_heart  <- t_dna_heart  / denom_dna + t_macro_heart  / denom_macro
    kp_spleen <- t_dna_spleen / denom_dna + t_macro_spleen / denom_macro
    kp_liver  <- t_dna_liver  / denom_dna + t_macro_liver  / denom_macro
    kp_kidney <- t_dna_kidney / denom_dna + t_macro_kidney / denom_macro
    kp_brain  <- t_dna_brain  / denom_dna + t_macro_brain  / denom_macro
    kp_rem    <- t_dna_rem    / denom_dna + t_macro_rem    / denom_macro

    c_lung    <- lung          / v_lung
    c_heart   <- heart         / v_heart
    c_spleen  <- spleen        / v_spleen
    c_liver   <- liver         / v_liver
    c_kidney  <- kidney        / v_kidney
    c_brain   <- brain         / v_brain
    c_isf_rem <- is_remainder  / v_isf_rem
    c_int_rem <- int_remainder / v_int_rem

    # ===== ODEs (Eqs 9-11) =====
    d/dt(lung)   <- q_lung   * (c_plasma - c_lung   / kp_lung)
    d/dt(heart)  <- q_heart  * (c_plasma - c_heart  / kp_heart)
    d/dt(spleen) <- q_spleen * (c_plasma - c_spleen / kp_spleen)
    d/dt(brain)  <- q_brain  * (c_plasma - c_brain  / kp_brain)

    d/dt(liver)  <- q_liver  * (c_plasma - c_liver  / kp_liver) -
                    cl_int_h * fu * (c_liver  / kp_liver)
    d/dt(kidney) <- q_kidney * (c_plasma - c_kidney / kp_kidney) -
                    cl_int_r * fu * (c_kidney / kp_kidney)

    flux_perm <- ps_rem * fu * (c_isf_rem - c_int_rem / kp_rem)
    d/dt(is_remainder)  <- q_rem * (c_plasma - c_isf_rem) - flux_perm
    d/dt(int_remainder) <- flux_perm

    d/dt(central) <-
        q_lung   * (c_lung   / kp_lung   - c_plasma) +
        q_heart  * (c_heart  / kp_heart  - c_plasma) +
        q_spleen * (c_spleen / kp_spleen - c_plasma) +
        q_liver  * (c_liver  / kp_liver  - c_plasma) +
        q_kidney * (c_kidney / kp_kidney - c_plasma) +
        q_brain  * (c_brain  / kp_brain  - c_plasma) +
        q_rem    * (c_isf_rem            - c_plasma)

    Cc <- c_plasma
    Cc ~ prop(propSd)
  })
}
