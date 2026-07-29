An_2012_mitoxantrone_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). Semi-mechanistic whole-body PBPK model 3 for",
    "mitoxantrone (Novantrone) in male ND4 Swiss Webster mice (24-32 g)",
    "after a single 5 mg/kg intravenous bolus (penile vein) (An and",
    "Morris 2012, AAPS J). Seven physiological tissue compartments:",
    "central plasma plus six perfusion-limited well-stirred organs (lung,",
    "heart, spleen, liver, kidney, brain), with a permeability-limited",
    "remainder compartment that lumps muscle, fat, bone, intestine, and",
    "skin and resolves into an interstitial (is_remainder) and an",
    "intracellular (int_remainder) subspace coupled by a permeability-",
    "surface area product PS. Hepatic and renal elimination act on the",
    "unbound cellular concentration of liver and kidney via well-stirred",
    "intrinsic clearances Clint_H and Clint_R. Saturable tissue binding",
    "to DNA (capacity T_DNA, affinity K_DNA) and macromolecular protein",
    "(capacity T_macro, affinity K_macro) is encoded as a Cp-dependent",
    "effective tissue:plasma partition coefficient Kp_eff(Cp) that varies",
    "instantaneously with plasma concentration (Eqs. 9-11 of the paper).",
    "Plasma unbound fraction is fixed to 0.2; tissue binding affinities",
    "are shared across all organs. The model is intended for typical-",
    "value simulation of mouse plasma and tissue concentration-time",
    "profiles."
  )
  reference <- paste(
    "An G, Morris ME. A Physiologically Based Pharmacokinetic Model of",
    "Mitoxantrone in Mice and Scale-up to Humans: a Semi-Mechanistic",
    "Model Incorporating DNA and Protein Binding. AAPS J.",
    "2012;14(2):352-364. doi:10.1208/s12248-012-9344-7.",
    sep = " "
  )
  vignette <- "An_2012_mitoxantrone_pbpk"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mouse body weight. The An and Morris (2012) cohort weighed",
        "24-32 g with mean 27.4 g; physiological organ volumes and blood",
        "flows in Table I were scaled from the 25 g Davies and Morris",
        "(1993) reference mouse by the An and Morris BW1/BW2 ratio (Eq 1)",
        "to the actual mean 27.4 g. In model() the WT covariate scales",
        "every organ volume and blood flow linearly from the Table I",
        "values for a 27.4 g reference mouse (multiplier WT / 0.0274 kg),",
        "matching the An and Morris Eq 1 BW2/BW1 scaling. WT also enters",
        "as the IV bolus mass multiplier (dose mg = 5 mg/kg * WT)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "mouse (male ND4 Swiss Webster)",
    n_subjects     = 27L,
    n_studies      = 1L,
    age_range      = "adult",
    weight_range   = "24-32 g (mean 27.4 g)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy, non-tumour-bearing adult male ND4 Swiss Webster mice",
      "(Harlan Labs, Indianapolis, IN). Housed at 12 h light/dark with",
      "standard diet and water ad libitum. Destructive sampling at 8",
      "post-dose times (5, 15, 30 min and 1, 2, 4, 8, 48 h) with 3-4",
      "mice per time point. Animal study approved by the IACUC,",
      "University at Buffalo."
    ),
    dose_range     = paste(
      "5 mg/kg single intravenous bolus via penile vein. Mitoxantrone",
      "dissolved in saline at 1.25 mg/mL and delivered as 4 microL/g",
      "body weight."
    ),
    regions        = "University at Buffalo, New York, USA",
    notes          = paste(
      "Plasma and six tissue (lung, heart, spleen, liver, kidney, brain)",
      "concentrations were assayed by validated HPLC (LLOQ plasma",
      "5 ng/mL, tissues 12.5 ng/mL) and fit simultaneously by maximum",
      "likelihood in ADAPT 5 (Biomedical Simulations Resource) to a",
      "variance model Var(t) = (sigma1 + sigma2 * Y(t))^2 with sigma1,",
      "sigma2 not reported numerically; this model file therefore",
      "carries no IIV and uses a placeholder fixed proportional",
      "residual error for syntactic completeness only - see the",
      "vignette Assumptions and deviations section. The mouse fit",
      "anchored the K_DNA, K_macro and tissue T_macro estimates that",
      "are reused in An_2012_mitoxantrone_human_pbpk."
    )
  )

  ini({
    # Estimated parameters of An and Morris 2012 Model 3 (Table II).
    # The paper reports T_DNA and T_macro as concentrations in nmol/g
    # tissue, which equal micromolar (uM) under the standard
    # assumption that 1 g tissue ~ 1 mL. K_DNA and K_macro are
    # reported in uM. Reported CV% are in the trailing comment.

    # DNA binding capacities -- estimated where the experimental
    # mouse-tissue DNA content was not reported in the literature
    # (lung, spleen, brain) and for the lumped remainder.
    lt_dna_lung   <- log(16.2);   label("Tissue DNA binding capacity, lung (uM)")            # An 2012 Table II Model 3 (T_DNA lung = 16.2, RSE 13%)
    lt_dna_spleen <- log(18.4);   label("Tissue DNA binding capacity, spleen (uM)")          # An 2012 Table II Model 3 (T_DNA spleen = 18.4, RSE 18%)
    lt_dna_brain  <- log(0.10);   label("Tissue DNA binding capacity, brain (uM)")           # An 2012 Table II Model 3 (T_DNA brain = 0.10, RSE 11%)
    lt_dna_rem    <- log(4.82);   label("Tissue DNA binding capacity, remainder (uM)")       # An 2012 Table II Model 3 (T_DNA other = 4.82, RSE 43%)

    # Protein (macromolecule) binding capacities -- estimated where the
    # mouse-tissue cardiolipin content was not available (lung,
    # spleen, remainder) and in the liver (literature cardiolipin
    # value 44.6 uM did not reproduce the liver concentration profile;
    # see paper Results section).
    lt_macro_lung   <- log(19.0);  label("Tissue protein binding capacity, lung (uM)")         # An 2012 Table II Model 3 (T_macro lung = 19.0, RSE 87%)
    lt_macro_spleen <- log(3.54);  label("Tissue protein binding capacity, spleen (uM)")       # An 2012 Table II Model 3 (T_macro spleen = 3.54, RSE >300%)
    lt_macro_liver  <- log(514);   label("Tissue protein binding capacity, liver (uM)")        # An 2012 Table II Model 3 (T_macro liver = 514, RSE 45%)
    lt_macro_rem    <- log(4.67);  label("Tissue protein binding capacity, remainder (uM)")    # An 2012 Table II Model 3 (T_macro other = 4.67, RSE >300%)

    # Binding affinities -- one global value for DNA and one for the
    # macromolecular protein binding site, shared across all organs.
    lk_dna   <- log(0.0013); label("DNA binding dissociation constant K_DNA (uM)")           # An 2012 Table II Model 3 (K_DNA = 0.0013, RSE 29%)
    lk_macro <- log(1.44);   label("Protein binding dissociation constant K_macro (uM)")     # An 2012 Table II Model 3 (K_macro = 1.44, RSE 52%)

    # Intrinsic clearances (well-stirred liver / kidney). Reported in
    # mL/min in the paper; the model uses L/h, so 29.5 mL/min = 1.77 L/h
    # and 2.14 mL/min = 0.128 L/h.
    lcl_int_h <- log(1.77);  label("Hepatic intrinsic clearance Clint_H (L/h)")              # An 2012 Results section (Clint_H = 29.5 mL/min, RSE 22%; see vignette Errata for the 29.5 vs 25.2 mL/min Results-vs-Discussion discrepancy)
    lcl_int_r <- log(0.128); label("Renal intrinsic clearance Clint_R (L/h)")                # An 2012 Results section (Clint_R = 2.14 mL/min, RSE 54%; see vignette Errata for the 0.21 mL/min apparent-renal-clearance inconsistency)

    # Permeability-surface area product for the permeability-limited
    # remainder compartment. 1.44 mL/min = 0.0864 L/h.
    lps_rem   <- log(0.0864); label("Permeability-surface area product, remainder PS (L/h)") # An 2012 Results section (PS = 1.44 mL/min, RSE not reported)

    # Residual error -- An and Morris 2012 used an ADAPT 5 maximum-
    # likelihood variance model Var(t) = (sigma1 + sigma2 * Y(t))^2,
    # but the numeric sigma1, sigma2 values are not reported in the
    # paper. nlmixr2 model definitions require a residual-error term,
    # so propSd is set to a small fixed placeholder. The model is
    # intended for typical-value simulation (no IIV); this residual
    # value should not be used as an inferential estimate. See the
    # vignette Assumptions and deviations section.
    propSd <- fixed(0.10); label("Proportional residual error placeholder (fraction)")       # not reported in An 2012; placeholder for syntactic completeness only
  })

  model({
    # ===== Drug physicochemistry =====
    # Mitoxantrone molecular weight (g/mol). DrugBank DB01204 / PubChem
    # CID 4212 report 444.49 g/mol; the An and Morris 2012 paper uses
    # micromolar concentrations for K_DNA / K_macro / T_DNA / T_macro
    # consistent with this value.
    mw_mitox <- 444.49

    # Plasma unbound fraction. Fixed to 0.2 throughout the An and Morris
    # 2012 fit (paper Data Analysis section: "In the data fitting, fu
    # was fixed to 0.2").
    fu <- 0.2

    # ===== Mouse physiology (Table I, Model 3 columns) =====
    # An and Morris 2012 Table I reports organ volumes and blood flows
    # for the mean 27.4 g study mouse. Eq 1 of the paper scales these
    # values linearly with body weight (Q2 = BW2 * Q1 / BW1). The WT
    # covariate enters here as the linear scaler relative to the
    # 27.4 g reference, so a user simulating a 25 g or 32 g mouse gets
    # appropriately scaled Q_i and V_i.
    wt_ref <- 0.0274                # An and Morris 2012 mean study weight (kg)
    wt_scale <- WT / wt_ref         # Eq 1 linear scaling factor

    # Organ volumes (Table I Model 3 column, mL) scaled by WT, then
    # converted mL -> L for unit consistency with Clint, PS in L/h.
    vp        <- 1.2   / 1000 * wt_scale   # plasma
    v_lung    <- 0.136 / 1000 * wt_scale
    v_heart   <- 0.144 / 1000 * wt_scale
    v_spleen  <- 0.06  / 1000 * wt_scale
    v_liver   <- 1.18  / 1000 * wt_scale
    v_kidney  <- 0.489 / 1000 * wt_scale
    v_brain   <- 0.36  / 1000 * wt_scale
    v_isf_rem <- 8.26  / 1000 * wt_scale   # extracellular (ISF) of remainder
    v_int_rem <- 16.77 / 1000 * wt_scale   # intracellular of remainder

    # Organ blood flow rates (Q) -- Table I Model 3 column (mL/min)
    # scaled linearly by WT per Eq 1, converted mL/min -> L/h via
    # * 60 / 1000.
    q_lung   <- 9.00 * 60 / 1000 * wt_scale   # = 0.540 L/h at WT_ref
    q_heart  <- 0.50 * 60 / 1000 * wt_scale   # = 0.030 L/h
    q_spleen <- 0.05 * 60 / 1000 * wt_scale   # = 0.003 L/h
    q_liver  <- 1.21 * 60 / 1000 * wt_scale   # = 0.0726 L/h
    q_kidney <- 1.99 * 60 / 1000 * wt_scale   # = 0.1194 L/h
    q_brain  <- 1.00 * 60 / 1000 * wt_scale   # = 0.060 L/h
    q_rem    <- 4.32 * 60 / 1000 * wt_scale   # = 0.2592 L/h

    # ===== Fixed mouse tissue T_DNA and T_macro values (Table I) =====
    # Mouse tissue DNA content was measurable for liver, heart, kidney
    # only (Gustafson 2002 DMD ref 17); the other organs are estimated
    # in ini(). Cardiolipin content was used for heart and kidney
    # T_macro (Courtade 1967 ref 27).
    t_dna_heart    <- 20.6   # An 2012 Table I (T_DNA heart, mouse)
    t_dna_liver    <- 26.0   # An 2012 Table I (T_DNA liver, mouse)
    t_dna_kidney   <- 50.2   # An 2012 Table I (T_DNA kidney, mouse)
    t_macro_heart  <- 43.8   # An 2012 Table I (T_macro heart, mouse cardiolipin)
    t_macro_kidney <- 52.3   # An 2012 Table I (T_macro kidney, mouse cardiolipin)

    # The brain T_macro is not reported in Table I (Brain T_macro = "-")
    # and was not estimated in the final model 3; the paper implicitly
    # omits the protein-binding term for brain. Set to 0 so the brain
    # Kp_eff reduces to the DNA-only form.
    t_macro_brain <- 0

    # ===== Individual structural parameters =====
    # Estimated parameters back-transformed from the log scale.
    t_dna_lung   <- exp(lt_dna_lung)
    t_dna_spleen <- exp(lt_dna_spleen)
    t_dna_brain  <- exp(lt_dna_brain)
    t_dna_rem    <- exp(lt_dna_rem)

    t_macro_lung   <- exp(lt_macro_lung)
    t_macro_spleen <- exp(lt_macro_spleen)
    t_macro_liver  <- exp(lt_macro_liver)
    t_macro_rem    <- exp(lt_macro_rem)

    k_dna   <- exp(lk_dna)
    k_macro <- exp(lk_macro)

    cl_int_h <- exp(lcl_int_h)
    cl_int_r <- exp(lcl_int_r)
    ps_rem   <- exp(lps_rem)

    # ===== Concentrations (mg/L) and plasma uM conversion =====
    # State variables hold mass in mg (matching the units = "mg"
    # dosing convention). Volume in L gives concentrations in mg/L.
    c_plasma <- central / vp

    # Convert plasma mg/L to uM for the binding equation, which is
    # parameterised in uM throughout. 1 mg/L * 1000 ug/mg / MW g/mol =
    # umol/L = uM.
    cp_uM <- c_plasma * 1000 / mw_mitox

    # ===== Cp-dependent effective tissue:plasma partition coefficients =====
    # An and Morris 2012 Eq. 9 inline definition: at the prevailing
    # plasma concentration Cp, the effective Kp combines saturable DNA
    # and macromolecular protein binding,
    #   Kp_eff(Cp) = T_DNA  / (K_DNA  / fu + Cp)
    #              + T_macro/ (K_macro/ fu + Cp)
    # where Cp, T_DNA, T_macro, K_DNA, K_macro are all in uM. Kp_eff
    # is dimensionless. As Cp -> 0 the partition saturates at the
    # maximum tissue:plasma ratio; as Cp -> infinity Kp_eff drops as
    # 1/Cp, mimicking high-concentration saturation of the binding
    # sites.
    denom_dna   <- k_dna   / fu + cp_uM
    denom_macro <- k_macro / fu + cp_uM

    kp_lung   <- t_dna_lung   / denom_dna + t_macro_lung   / denom_macro
    kp_heart  <- t_dna_heart  / denom_dna + t_macro_heart  / denom_macro
    kp_spleen <- t_dna_spleen / denom_dna + t_macro_spleen / denom_macro
    kp_liver  <- t_dna_liver  / denom_dna + t_macro_liver  / denom_macro
    kp_kidney <- t_dna_kidney / denom_dna + t_macro_kidney / denom_macro
    kp_brain  <- t_dna_brain  / denom_dna + t_macro_brain  / denom_macro
    kp_rem    <- t_dna_rem    / denom_dna + t_macro_rem    / denom_macro

    # Tissue concentrations (mg/L) from amount states (mg) / volume (L).
    c_lung    <- lung          / v_lung
    c_heart   <- heart         / v_heart
    c_spleen  <- spleen        / v_spleen
    c_liver   <- liver         / v_liver
    c_kidney  <- kidney        / v_kidney
    c_brain   <- brain         / v_brain
    c_isf_rem <- is_remainder  / v_isf_rem
    c_int_rem <- int_remainder / v_int_rem

    # ===== ODEs (Eqs 9-11 of An and Morris 2012, with elimination) =====
    # Well-stirred non-eliminating organs (Eq 9):
    #   V_t * dC_t/dt = Q_t * (C_p - C_t / Kp_eff(C_p))
    # In amount form (state * = mass in mg):
    #   d/dt(organ) = Q_t * (C_p - C_t / Kp_eff)
    d/dt(lung)   <- q_lung   * (c_plasma - c_lung   / kp_lung)
    d/dt(heart)  <- q_heart  * (c_plasma - c_heart  / kp_heart)
    d/dt(spleen) <- q_spleen * (c_plasma - c_spleen / kp_spleen)
    d/dt(brain)  <- q_brain  * (c_plasma - c_brain  / kp_brain)

    # Well-stirred eliminating organs (liver, kidney) add the intrinsic
    # clearance term acting on the unbound concentration at the organ
    # outflow (Cp_equivalent = C_t / Kp_eff). For mg/h units the
    # elimination flux is Clint * fu * (C_t / Kp_eff).
    d/dt(liver)  <- q_liver  * (c_plasma - c_liver  / kp_liver) -
                    cl_int_h * fu * (c_liver  / kp_liver)
    d/dt(kidney) <- q_kidney * (c_plasma - c_kidney / kp_kidney) -
                    cl_int_r * fu * (c_kidney / kp_kidney)

    # Permeability-limited remainder (Eqs 10-11):
    #   V_ISF * dC_ISF/dt = Q_rem * (C_p - C_ISF)
    #                      - PS * fu * (C_ISF - C_int / Kp_eff(C_p))
    #   V_int * dC_int/dt = PS * fu * (C_ISF - C_int / Kp_eff(C_p))
    # In amount form:
    flux_perm <- ps_rem * fu * (c_isf_rem - c_int_rem / kp_rem)
    d/dt(is_remainder)  <- q_rem * (c_plasma - c_isf_rem) - flux_perm
    d/dt(int_remainder) <- flux_perm

    # Plasma mass balance: inflow = sum of venous returns; outflow =
    # total flow to all organs. The An and Morris 2012 schematic
    # (Fig 2) treats every tissue as exchanging in parallel with one
    # central plasma pool, so plasma sees both the lung flow (Q_lung =
    # 9 mL/min, equal to cardiac output) and the sum of all the other
    # parallel organ flows (= 9 mL/min). This is the paper's stated
    # topology (see vignette Assumptions and deviations).
    d/dt(central) <-
        q_lung   * (c_lung   / kp_lung   - c_plasma) +
        q_heart  * (c_heart  / kp_heart  - c_plasma) +
        q_spleen * (c_spleen / kp_spleen - c_plasma) +
        q_liver  * (c_liver  / kp_liver  - c_plasma) +
        q_kidney * (c_kidney / kp_kidney - c_plasma) +
        q_brain  * (c_brain  / kp_brain  - c_plasma) +
        q_rem    * (c_isf_rem            - c_plasma)

    # ===== Observation: plasma concentration (mg/L = ug/mL) =====
    Cc <- c_plasma
    Cc ~ prop(propSd)
  })
}
