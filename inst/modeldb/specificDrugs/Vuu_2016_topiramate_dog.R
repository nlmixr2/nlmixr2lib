Vuu_2016_topiramate_dog <- function() {
  description <- paste(
    "Preclinical (dog). Population two-compartment intravenous PK model for",
    "topiramate (TPM) in dogs with naturally-occurring epilepsy (Vuu 2016).",
    "Stable-labelled TPM was given as a 5-min IV infusion at 10 mg/kg (n = 4)",
    "or 20 mg/kg (n = 3); pooled across the low- and high-dose data, a two-",
    "compartment model with first-order elimination from the central",
    "compartment described the disposition best. Concomitant phenobarbital",
    "(CONMED_PB) was identified as an enzyme-inducer covariate on systemic",
    "clearance via an exponential effect (Cl = tvCl * exp(dCl * CONMED_PB)),",
    "yielding a 5.64-fold higher CL in PB-coadministered dogs. Per-kg",
    "structural parameters (Vc, Vp, CL, Q) are scaled to absolute units by",
    "individual body weight (WT, kg) inside the model; the dose in the event",
    "table is therefore absolute mg (mg/kg dose times WT). IIV is exponential",
    "on Vc and CL; residual error is proportional (~15%).",
    sep = " "
  )
  reference   <- paste(
    "Vuu I, Coles LD, Maglalang P, Leppik IE, Worrell G, Crepeau D,",
    "Mishra U, Cloyd JC, Patterson EE.",
    "Intravenous Topiramate: Pharmacokinetics in Dogs with Naturally",
    "Occurring Epilepsy.",
    "Front Vet Sci. 2016 Dec 5;3:107.",
    "doi:10.3389/fvets.2016.00107.",
    sep = " "
  )
  vignette <- "Vuu_2016_topiramate_dog"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L (equivalent to ug/mL)",
    dosing_notes  = paste(
      "Doses in the source paper are mg/kg; the event-table AMT is the",
      "absolute mg dose (mg/kg dose times the dog's body weight). Per-kg",
      "structural PK parameters (Vc, Vp, CL, Q) are converted to absolute",
      "units inside model() by multiplying by the WT covariate (kg).",
      sep = " "
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Individual body weight at the time of TPM dosing.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Required to convert the per-kg structural PK parameters reported",
        "by Vuu 2016 (Vc, Vp in mL/kg; CL, Q in mL/(kg*min)) into absolute",
        "units (L, L/h). The cohort range is 15 - 35 kg (Table 1).",
        "Vuu 2016 did not estimate a separate body-weight covariate effect",
        "(the model is linearly per-kg by construction), so the scaling is",
        "an exact linear-in-WT mapping equivalent to allometric scaling",
        "with exponent 1.0 on Vc, Vp, CL, and Q.",
        sep = " "
      ),
      source_name = "Weight"
    ),
    CONMED_PB = list(
      description        = paste(
        "Concomitant phenobarbital (PB) coadministration indicator at the",
        "time of TPM dosing: 1 = dog on chronic PB maintenance therapy,",
        "0 = no concomitant PB.",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = paste(
        "Exponential effect on CL on the log scale:",
        "Cl_i = tvCl * exp(dCl * CONMED_PB) * exp(eta_Cl).",
        "With dCl = 1.73 (Table S1), CONMED_PB = 1 produces a",
        "exp(1.73) = 5.64-fold higher Cl than CONMED_PB = 0",
        "(consistent with the Results-text 5.6-fold figure).",
        "Vuu 2016 cohort: dogs 1, 2, and 5 were on PB (Table 1);",
        "dogs 3 and 4 were on no antiseizure maintenance therapy.",
        "Primidone is conventionally pooled with phenobarbital under",
        "CONMED_PB (no dog in this cohort was on primidone).",
        sep = " "
      ),
      source_name = "PB"
    )
  )

  population <- list(
    species        = "dog (mixed-breed and beagle, naturally occurring epilepsy)",
    n_subjects     = 5L,
    n_studies      = 1L,
    age_range      = "3 - 9 years",
    weight_range   = "15 - 35 kg",
    sex_female_pct = 20.0,
    disease_state  = paste(
      "Naturally-occurring canine epilepsy (focal with secondary",
      "generalisation, or in remission). Three dogs (1, 2, 5) on chronic",
      "antiseizure-drug maintenance including phenobarbital; two dogs",
      "(3, 4) on no concomitant antiseizure therapy.",
      sep = " "
    ),
    dose_range     = paste(
      "Low-dose IV TPM 10 mg/kg over a 5-min infusion (n = 4 dogs",
      "ID 1-4); high-dose IV TPM 20 mg/kg over a 5-min infusion",
      "(n = 3 dogs ID 3-5). One hour after the 10 mg/kg IV dose,",
      "dogs 1-4 also received a 5 mg/kg unlabelled oral tablet for",
      "the bioavailability arm (analysed by NCA only; not part of",
      "the popPK fit).",
      sep = " "
    ),
    regions        = "United States (University of Minnesota Veterinary College).",
    notes          = paste(
      "Stable-isotope-labelled TPM (six 13C) used for the IV formulation",
      "(10 mg/mL in 10% Captisol) to permit simultaneous IV / oral dosing",
      "while distinguishing the two formulations analytically. Pooled",
      "low-/high-dose IV data fit in Phoenix NLME 1.3 via first-order",
      "conditional estimation extended least squares; exponential IIV on",
      "Vc and CL, proportional residual error.",
      "iEEG and oral NCA arms are described in the paper but not part",
      "of the population PK structural model implemented here.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural PK parameters -- Table S1 of Vuu 2016 (Supplementary
    # Material, population two-compartment fit pooling the 10 and 20 mg/kg
    # IV doses). Per-kg values are converted to absolute units in model()
    # by multiplying by the WT covariate (kg). Volume mL/kg / 1000 = L/kg;
    # clearance mL/(kg*min) * 60 / 1000 = L/(kg*h).
    # =====================================================================
    lvc <- log(0.376)
    label("Per-kg central volume of distribution Vc (L/kg)")             # Table S1: tvV = 376 mL/kg (SE 72.4, CV% 19.2)
    lvp <- log(0.298)
    label("Per-kg peripheral volume of distribution Vp (L/kg)")          # Table S1: tvV2 = 298 mL/kg (SE 56.0, CV% 18.7)
    lcl <- log(1.84 * 60 / 1000)
    label("Per-kg systemic clearance CL (L/(kg*h))")                     # Table S1: tvCl = 1.84 mL/(kg*min) (SE 0.08, CV% 4.52); 1.84 mL/(kg*min) * 60 min/h / 1000 mL/L = 0.1104 L/(kg*h)
    lq  <- log(21.0 * 60 / 1000)
    label("Per-kg intercompartmental clearance Q (L/(kg*h))")            # Table S1: tvCl2 = 21.0 mL/(kg*min) (SE 9.37, CV% 44.7); 21.0 mL/(kg*min) * 60 min/h / 1000 mL/L = 1.26 L/(kg*h)

    # Exponential CYP-inducer effect on CL: cl = tvcl * exp(e_conmed_pb_cl
    # * CONMED_PB), with e_conmed_pb_cl on the log scale (paper symbol dCl).
    e_conmed_pb_cl <- 1.73
    label("Log-scale exponential effect of CONMED_PB on CL (unitless)")  # Table S1: dCl = 1.73 (SE 0.13, CV% 7.66); exp(1.73) = 5.64-fold higher Cl when on PB, matching Results-text "5.6-fold"

    # =====================================================================
    # IIV (exponential model; Phoenix NLME "BSV" estimates reported as
    # variances of the log-normal random effect, i.e. omega^2). Table S1
    # of Vuu 2016: BSV(V) = 0.08 (RSE 24.6%), BSV(Cl) = 0.02 (RSE 53.4%).
    # =====================================================================
    etalvc ~ 0.08  # Table S1: BSV(V) = 0.08, RSE 24.6%, shrinkage 9.3%; ~29% CV on Vc
    etalcl ~ 0.02  # Table S1: BSV(Cl) = 0.02, RSE 53.4%, shrinkage 9.18%; ~14% CV on Cl

    # =====================================================================
    # Residual error. Vuu 2016 used a "multiplicative" (proportional)
    # error model on plasma TPM. Table S1: CV% = 14.9 (SE 1.71, RSE 11.5%).
    # The main Results text quotes this as "~15%". Encoded as a fractional
    # proportional SD on Cc.
    # =====================================================================
    propSd <- 0.149
    label("Proportional residual error on plasma TPM (fraction)")        # Table S1: Residual error CV% = 14.9 (SE 1.71, RSE 11.5%)
  })

  model({
    # Individual structural parameters: per-kg values multiplied by the WT
    # covariate to get absolute units (L, L/h). The CONMED_PB exponential
    # effect on Cl is on the log scale (paper: cl = tvcl * exp(dCl * PB)).
    cl <- exp(lcl + etalcl + e_conmed_pb_cl * CONMED_PB) * WT
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq) * WT
    vp <- exp(lvp) * WT

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system, IV infusion direct to central (no depot).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma TPM concentration in mg/L (dose in mg, vc in L); mg/L is
    # numerically equivalent to ug/mL, matching the paper's reported units.
    Cc <- central / vc

    # Proportional residual error on plasma TPM.
    Cc ~ prop(propSd)
  })
}
