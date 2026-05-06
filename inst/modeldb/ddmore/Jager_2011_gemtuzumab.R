Jager_2011_gemtuzumab <- function() {
  description <- "Mechanism-based PKPD model for the antibody-drug conjugate gemtuzumab ozogamicin (GO; Mylotarg) in patients with acute myeloid leukemia, as packaged in DDMORE Foundation Model Repository entry DDMODEL00000229. The model couples drug pharmacokinetics with explicit binding to the cell-surface antigen CD33: free drug in a central compartment binds free CD33 receptor to form a drug-receptor complex that is internalized; the toxic ozogamicin component then drives linear depletion of leukemic blast cells. The DDMORE entry extends the original Jager 2011 PK structure by adding a peripheral drug compartment and re-estimating all parameters simultaneously in Monolix, so it is not a literal reproduction of the published model — see the validation vignette for the comparison."
  reference <- paste(
    "Jager E, van der Velden VHJ, te Marvelde JG, Walter RB, Agur Z, Vainstein V (2011).",
    "Targeted drug delivery by gemtuzumab ozogamicin: mechanism-based mathematical model for treatment strategy improvement and therapy individualization.",
    "PLoS One 6(9):e24265.",
    "doi:10.1371/journal.pone.0024265.",
    "DDMORE Foundation Model Repository: DDMODEL00000229 (Monolix re-fit with an added peripheral compartment).",
    sep = " "
  )
  vignette <- "Jager_2011_gemtuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000229"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Patients with acute myeloid leukemia (AML) receiving the CD33-directed antibody-drug conjugate gemtuzumab ozogamicin (GO; Mylotarg).",
    dose_range     = "Intravenous gemtuzumab ozogamicin. Specific clinical doses, infusion durations, and study design are reported in the original Jager 2011 publication; the DDMORE bundle does not ship a dataset and the publication itself is not on disk in this worktree, so dose-level detail could not be extracted.",
    regions        = NA_character_,
    notes          = paste(
      "Demographic detail is not reproduced in the DDMORE bundle for DDMODEL00000229,",
      "and the Jager 2011 PLoS One paper itself is not on disk under the literature tree.",
      "The DDMORE bundle's `Model_Accommodations.txt` documents that the model parameters",
      "were re-evaluated in Monolix using the same data as the original publication,",
      "with the additional change that a peripheral PK compartment was introduced and all",
      "parameters were estimated simultaneously rather than step-wise as in Jager 2011.",
      "The original 2011 paper itself is a single-cohort population-PKPD analysis on AML",
      "patients; consult the publication for the complete population description."
    )
  )

  ini({
    # Parameter VALUES come from `GO_PK_model.mdl` parObj$STRUCTURAL and parObj$VARIABILITY.
    # The DDMORE bundle for DDMODEL00000229 ships only the MDL DSL source (`GO_PK_model.mdl`)
    # and its PharmML rendering (`GO_PK_model.xml`) — there is NO `Output_real_*.lst` listing.
    # `Model_Accommodations.txt` describes the parObj as the Monolix re-fit of Jager 2011 on
    # the same data as the publication; absent a `.lst` for cross-checking, the parObj values
    # are treated as the bundle's authoritative final estimates. See the vignette Errata.

    # Drug PK — central volume and rate-constant parameterization (the MDL keeps drug PK in
    # rate-constant form rather than CL/V form, with a peripheral compartment added on top of
    # the original 2011 structure).
    lvc    <- log(5.42);    label("Drug central volume of distribution V1 (L)")                                            # parObj POP_v1
    lk     <- log(0.0135);  label("Drug linear elimination rate k from central (1/h)")                                     # parObj POP_k
    lkm    <- log(0.00757); label("Drug central -> peripheral first-order rate km (1/h)")                                  # parObj POP_km
    lkn    <- log(0.0185);  label("Drug peripheral -> central first-order rate kn (1/h)")                                  # parObj POP_kn

    # CD33-directed binding kinetics — drug binds free receptor to form an internalized complex.
    lkb    <- log(9.24e5);  label("Drug-receptor association (binding) rate kb (MDL native rate units)")                   # parObj POP_kb
    lku    <- log(310);     label("Drug-receptor dissociation (unbinding) rate ku (1/h)")                                  # parObj POP_ku
    lki    <- log(0.624);   label("Drug-receptor complex internalization rate ki (1/h)")                                   # parObj POP_ki

    # Free CD33 receptor turnover (zero-order production rp, first-order elimination ke).
    lke    <- log(0.199);   label("Free receptor first-order elimination rate ke (1/h)")                                   # parObj POP_ke
    lrp    <- log(823);     label("Free receptor zero-order production rate rp (MDL native receptor units / h)")           # parObj POP_rp

    # Leukemic blast cell dynamics (linear depletion at typical-population rate alph).
    lalph  <- log(0.0755);  label("Leukemic blast first-order depletion rate alpha (1/h)")                                 # parObj POP_alph
    ln0    <- log(1.32e-5); label("Initial leukemic blast count n0 at t = 0 (MDL native cell-count units)")                # parObj POP_n0

    # Inter-individual variability — only ki, ke, and rp carry non-zero IIV in the parObj.
    # The MDL declares omega_X with `type is sd`, so the parObj value IS the standard deviation;
    # nlmixr2 expects a variance for the `~` form, hence the squared values below.
    # The MDL parObj also declares an OMEGA correlation block among ETA_k, ETA_alph, ETA_n0,
    # ETA_v1, ETA_kn, ETA_km, but every one of those etas has omega = 0 (sd) in the parObj
    # so the correlation block is mathematically vestigial and is dropped here. See the
    # vignette Errata.
    etalki ~ 0.258^2  # parObj omega_ki (sd 0.258 -> variance 0.0666)
    etalke ~ 0.281^2  # parObj omega_ke (sd 0.281 -> variance 0.0790)
    etalrp ~ 0.611^2  # parObj omega_rp (sd 0.611 -> variance 0.3733)

    # Residual error - additive on the predicted central concentration output1 = central / vc.
    # MDL: Y = additiveError(additive = a, eps = EPS_Y, prediction = output1) with EPS_Y ~ N(0, 1).
    addSd  <- 0.000583; label("Additive residual error on central concentration (mg/L per the dosing/Vc unit assumption)")  # parObj a
  })

  model({
    # Individual parameters - log-linear with random effects only where the parObj reports a
    # non-zero IIV sd (ki, ke, rp); typical-value-only otherwise.
    vc   <- exp(lvc)
    k    <- exp(lk)
    km   <- exp(lkm)
    kn   <- exp(lkn)
    kb   <- exp(lkb)
    ku   <- exp(lku)
    ki   <- exp(lki + etalki)
    ke   <- exp(lke + etalke)
    rp   <- exp(lrp + etalrp)
    alph <- exp(lalph)
    n0   <- exp(ln0)

    # MDL inline constants from the MODEL_PREDICTION block. nav is Avogadro's number divided
    # by 1e23 (so absolute amounts in `central` are scaled to a `molecules / 1e23` working
    # range), and vblood is the assumed blood volume used to convert receptor density to a
    # plasma-equivalent concentration in the drug-side equation for `central`.
    nav    <- 6.02214
    vblood <- 4.8

    # Initial conditions per the MDL DEQ block. The free receptor starts at the steady-state
    # baseline rp / ke; the drug-receptor complex starts empty; leukemic blasts start at n0;
    # both drug compartments start empty (dosing populates `central` via amt).
    target(0)      <- rp / ke
    complex(0)     <- 0
    cells(0)       <- n0
    central(0)     <- 0
    peripheral1(0) <- 0

    # ODE system - DDMODEL00000229 GO_PK_model.mdl `MODEL_PREDICTION$DEQ`. Compartment mapping:
    #   A1 -> central       drug in the central plasma compartment
    #   A2 -> target        free CD33 receptor (TMDD-canonical name)
    #   A3 -> complex       drug-receptor complex (TMDD-canonical name)
    #   A4 -> cells         leukemic blast cell population
    #   A5 -> peripheral1   drug peripheral compartment (Monolix re-fit addition)
    d/dt(central)     <- (-kb * central * target + vc * ku * complex) * cells / (vblood * nav) - k * central - km * central + kn * peripheral1
    d/dt(target)      <- -(kb / vc) * central * target + ku * complex + rp - ke * target
    d/dt(complex)     <-  (kb / vc) * central * target - ku * complex - ki * complex
    d/dt(cells)       <- -alph * cells
    d/dt(peripheral1) <-  km * central - kn * peripheral1

    # Observation - drug concentration in the central compartment (output1 = A1 / v1).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
