Hennig_2006_itraconazole <- function() {
  description <- "Population PK model for oral itraconazole and its active metabolite hydroxy-itraconazole in paediatric cystic-fibrosis and bone-marrow-transplant patients (Hennig 2006). One-compartment parent + one-compartment metabolite with first-order absorption, first-order metabolic conversion (fm fixed to 1), allometric weight scaling on parent CL/F (0.75) and Vd/F (1.0), and formulation-specific ka and relative bioavailability for capsule vs oral solution."
  reference   <- "Hennig S, Wainwright CE, Bell SC, Miller H, Friberg LE, Charles BG. Population pharmacokinetics of itraconazole and its active metabolite hydroxy-itraconazole in paediatric cystic fibrosis and bone marrow transplant patients. Clin Pharmacokinet. 2006;45(11):1099-1114. doi:10.2165/00003088-200645110-00004"
  vignette    <- "Hennig_2006_itraconazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; constant within an individual in the source dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed allometric scaling: (WT/70)^0.75 on CL/F of itraconazole and (WT/70)^1.0 on Vd/F of itraconazole. No weight scaling on metabolite parameters per the source NMTRAN control stream. Reference weight 70 kg.",
      source_name        = "WT"
    ),
    CAPSULE = list(
      description        = "Formulation indicator at the dose record: 1 = capsule (Sporanox), 0 = oral solution.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral solution; relative bioavailability fixed to 1).",
      notes              = "Per-record covariate (each subject may switch formulations during follow-up). Switches absorption rate between ka(capsule) = 0.09 h^-1 and ka(solution) = 0.96 h^-1, and gates the IIV-bearing relative-bioavailability term so etalfdepot only contributes when CAPSULE = 1.",
      source_name        = "PREP"
    )
  )

  population <- list(
    n_subjects     = 49,
    n_studies      = 2,
    age_range      = "0.4-30 years (median 8 years; 5 adult cystic-fibrosis patients aged 19-30 years pooled with 44 paediatric patients)",
    age_median     = "8 years",
    weight_range   = "6.8-83.5 kg",
    weight_median  = "29.3 kg",
    sex_female_pct = 38.8,
    race_ethnicity = c(White = 95.9, Other = 4.1),
    disease_state  = "Paediatric cystic fibrosis (n = 29; itraconazole prescribed for allergic bronchopulmonary aspergillosis) or paediatric bone-marrow transplant (n = 20; itraconazole prescribed for fungal-infection prophylaxis); 5 adult cystic-fibrosis subjects from The Prince Charles Hospital pooled with the paediatric cohort for comparison.",
    dose_range     = "Oral itraconazole 1.5-12.5 mg/kg/day (median 5.4 mg/kg/day) as Sporanox capsules or oral cyclodextrin solution; switching between formulations was allowed mid-study after at least three consecutive doses of the most recent formulation.",
    regions        = "Australia (The Royal Children's Hospital, Brisbane, QLD; The Prince Charles Hospital, Brisbane, QLD).",
    n_observations = "227 itraconazole observations (median 0.264 mg/L) and 192 hydroxy-itraconazole observations (median 0.531 mg/L); LLOQ 0.075 mg/L for both analytes (sub-LLOQ values imputed as half-LLOQ = 0.0375 mg/L per Hennig 2006 Drug Analysis section).",
    food_status    = "75% of capsule doses were taken with food, 25% on an empty stomach; 9% of oral-solution doses were taken with food, 91% on an empty stomach. Food and acidic-beverage effects were tested as covariates on Frel but did not reach significance.",
    notes          = "Demographic and dosing summary from Hennig 2006 Table I. Sampling was empirical (1-2 samples per subject per visit), with preferred windows of 0-6 h and 4 h before next dose; samples below LLOQ assigned to LLOQ/2 = 0.0375 mg/L."
  )

  ini({
    # Structural parameters - parent itraconazole, reference 70 kg.
    # Final estimates from Hennig 2006 Table II "Final model" column.
    lcl <- log(35.5);   label("Apparent clearance of itraconazole, CL/F, at 70 kg (L/h)")                       # Hennig 2006 Table II, final model: 35.5 L/h
    lvc <- log(672);    label("Apparent volume of distribution of itraconazole, Vd/F, at 70 kg (L)")            # Hennig 2006 Table II, final model: 672 L

    # Structural parameters - hydroxy-itraconazole metabolite.
    # Reported THETA values in Table II are CLm/(F*fm) and Vd(m)/(F*fm); the model
    # internally rescales by the metabolite/parent molecular-weight ratio
    # (722.64/705.64) inside model() per the supplied .ctl, which converts the
    # parent-mass-equivalent compartment amount into the measured metabolite
    # mass concentration.
    lcl_ohi <- log(10.6);  label("Apparent clearance of hydroxy-itraconazole, CLm/(F*fm) (L/h)")                 # Hennig 2006 Table II, final model: 10.6 L/h
    lvc_ohi <- log(5.29);  label("Apparent volume of distribution of hydroxy-itraconazole, Vd(m)/(F*fm) (L)")    # Hennig 2006 Table II, final model: 5.29 L

    # Absorption: oral-solution ka is the typical-value baseline, the capsule
    # arm shifts ka multiplicatively via a binary CAPSULE covariate effect.
    lka          <- log(0.96);                                                                                   # Hennig 2006 Table II, final model: ka(oral solution) = 0.96 h^-1
    label("First-order absorption rate constant for oral solution (1/h)")
    e_capsule_ka <- log(0.09 / 0.96);                                                                            # Hennig 2006 Table II, final model: ka(capsule) = 0.09 h^-1; log-shift = log(0.09/0.96) = -2.367
    label("Log-ratio shift on ka for capsule formulation (unitless)")

    # Lag time (same for capsule and oral solution per the .ctl).
    llag <- log(19.1 / 60);                                                                                      # Hennig 2006 Table II, final model: t_lag = 19.1 min = 0.3183 h
    label("Absorption lag time (h)")

    # Relative bioavailability: oral solution Frel fixed to 1 (reference);
    # capsule Frel = exp(lfdepot) = 0.55 with the etalfdepot IIV applied
    # only on the capsule arm via the (lfdepot + etalfdepot) * CAPSULE form
    # in model().
    lfdepot <- log(0.55);                                                                                         # Hennig 2006 Table II, final model: Frel(capsule:solution) = 0.55
    label("Log relative bioavailability of capsule vs oral solution (unitless)")

    # Allometric scaling exponents - fixed per the paper's Methods (Eq. 5 and
    # surrounding text) and the .ctl lines 20-21.
    allo_cl <- fixed(0.75);  label("Allometric exponent on CL/F of itraconazole (fixed)")                         # Hennig 2006 Methods (Eq. 5): "exponent x was fixed to 0.75 for clearance"
    allo_vc <- fixed(1.0);   label("Allometric exponent on Vd/F of itraconazole (fixed)")                         # Hennig 2006 Methods (Eq. 5): "and 1 for volume of distribution"

    # IIV. Hennig 2006 reports omega as %CV (Table II footer); convert to
    # log-normal variance via omega^2 = log(1 + CV^2). Block on CL_itra and
    # Vd_itra has correlation rho = 0.69 (Random Effects Model paragraph,
    # page 1105); cov = rho * sqrt(omega2_CL) * sqrt(omega2_VC).
    # CL_itra: log(1 + 0.688^2) = 0.38731
    # Vd_itra: log(1 + 0.758^2) = 0.45404
    # cov   : 0.69 * sqrt(0.38731) * sqrt(0.45404) = 0.28938
    etalcl + etalvc ~ c(0.38731, 0.28938, 0.45404)                                                                # Hennig 2006 Table II + Random Effects Model paragraph (rho = 0.69)
    # CLm   : log(1 + 0.734^2) = 0.43098
    etalcl_ohi ~ 0.43098                                                                                          # Hennig 2006 Table II, final model: CV 73.4% on CLm
    # Frel (capsule arm only; gated by CAPSULE inside model()):
    # log(1 + 0.611^2) = 0.31734
    etalfdepot ~ 0.31734                                                                                          # Hennig 2006 Table II, final model: CV 61.1% on Frel

    # Residual error. Hennig 2006 Eq. 3 specifies an exponential residual
    # model Y = Cpred * exp(eps), which the paper also calls "proportional"
    # (Random Effects Model paragraph). Translated to the proportional /
    # linear-space convention used by nlmixr2lib: propSd is the reported
    # sigma value (~CV in linear space for the magnitudes here). See
    # vignette Assumptions and deviations for the equivalence note.
    propSd     <- 0.499;  label("Itraconazole proportional residual error (fraction)")                          # Hennig 2006 Table II, final model: sigma_itraconazole = 49.9%
    propSd_ohi <- 0.470;  label("Hydroxy-itraconazole proportional residual error (fraction)")                  # Hennig 2006 Table II, final model: sigma_m = 47.0%
  })

  model({
    # Reference body weight for allometric scaling (BWT_std = 70 kg per
    # Hennig 2006 Methods Eq. 5 and .ctl lines 20-21).
    ref_wt <- 70

    # Parent/metabolite molecular-weight ratio (722.64 / 705.64 from .ctl
    # lines 32 and 35); applied to cl_ohi and vc_ohi so that the
    # observation A(central_ohi) / vc_ohi is in the measured-metabolite
    # mass concentration units (mg/L of hydroxy-itraconazole).
    mw_ratio <- 722.64 / 705.64

    # Individual parameters - parent itraconazole.
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^allo_vc

    # Individual parameters - hydroxy-itraconazole metabolite (no IIV on
    # vc_ohi per Table II; mw_ratio applied to both clearance and volume
    # to match the .ctl).
    cl_ohi <- exp(lcl_ohi + etalcl_ohi) * mw_ratio
    vc_ohi <- exp(lvc_ohi) * mw_ratio

    # Absorption: solution baseline + capsule log-shift.
    ka <- exp(lka + e_capsule_ka * CAPSULE)

    # Relative bioavailability on the depot:
    #   CAPSULE = 0 (solution): exp(0) = 1 (reference)
    #   CAPSULE = 1 (capsule):  exp(lfdepot + etalfdepot) = ~0.55 with IIV
    fdepot <- exp((lfdepot + etalfdepot) * CAPSULE)

    # Micro-constants. fm fixed to 1 per the paper (Pharmacokinetic
    # Analysis paragraph): all eliminated parent appears as metabolite.
    kpm <- cl / vc      # parent -> metabolite formation rate constant
    kel <- cl_ohi / vc_ohi  # metabolite elimination rate constant

    # ODE system (matches the .ctl ADVAN5 GUTT / PARENT / METABOLITE chain
    # with K20 = 0 and K32 = 0).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kpm * central
    d/dt(central_ohi) <-  kpm * central - kel * central_ohi

    # Bioavailability and lag on the depot (lag is the same for both
    # formulations per Table II).
    f(depot)   <- fdepot
    lag(depot) <- exp(llag)

    # Observations.
    Cc     <- central / vc
    Cc_ohi <- central_ohi / vc_ohi

    Cc     ~ prop(propSd)
    Cc_ohi ~ prop(propSd_ohi)
  })
}
