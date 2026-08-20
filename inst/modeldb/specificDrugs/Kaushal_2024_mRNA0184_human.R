Kaushal_2024_mRNA0184_human <- function() {
  description <- "Human-scaled allometric projection (NO human data were fitted). Translational semi-mechanistic PK/PD model for mRNA-0184, a lipid-nanoparticle-encapsulated mRNA encoding human relaxin-2 fused to a variable light chain kappa domain (Rel2-vlk), for heart failure with reduced ejection fraction. This is the cynomolgus-monkey model of Kaushal 2024 forward-projected to a 70 kg adult by allometric scaling (exponent 1 on volumes, 0.75 on Rel2-vlk mRNA clearances, 0.85 on Rel2-vlk protein clearances). PK: 3-compartment plasma1-tissue-plasma2 redistribution of Rel2-vlk mRNA (both plasma compartments share V1; elimination from the tissue compartment). PD: Rel2-vlk protein production linear in a hypothetical effect-compartment mRNA concentration, distributed in a 2-compartment system with a single symmetric intercompartmental clearance and first-order central elimination. Used to support the 0.025 mg/kg every-2-weeks first-in-human starting dose for NCT05659264."
  reference <- paste(
    "Kaushal N, Attarwala H, Iqbal MJ, Saini R, Van L, Liang M.",
    "Translational pharmacokinetic/pharmacodynamic model for mRNA-0184,",
    "an investigational therapeutic for the treatment of heart failure.",
    "Clin Transl Sci. 2024;17(8):e13894. doi:10.1111/cts.13894.",
    "Structural equations from Supporting Information Table S1;",
    "human-scaled parameter values from Table 2.",
    sep = " "
  )
  vignette <- "Kaushal_2024_mRNA0184"
  paper_specific_compartments <- c("rel2vlk", "rel2vlk_p")

  units <- list(
    time          = "h",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  # Table S1 names these A1, A2, A3, Ce, protein1 and protein2. Note that
  # `effect` carries a CONCENTRATION (the paper's Ce), not an amount --
  # Table S1 writes dCe/dt = Ke0 * (C(t) - Ce(t)) directly in
  # concentration units.
  compartmentData <- list(
    central     = list(analyte = "Rel2-vlk mRNA (A1, plasma-1)",                 units = "ng",    specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "Rel2-vlk mRNA (A2, tissue / target site)",     units = "ng",    specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "Rel2-vlk mRNA (A3, plasma-2)",                 units = "ng",    specimen = "plasma", verified = TRUE),
    effect      = list(analyte = "Rel2-vlk mRNA (Ce, hypothetical effect compartment)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    rel2vlk     = list(analyte = "Rel2-vlk protein (protein1, central)",         units = "ng",    specimen = "plasma", verified = TRUE),
    rel2vlk_p   = list(analyte = "Rel2-vlk protein (protein2, peripheral)",      units = "ng",    specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "human (allometric projection from cynomolgus monkey; no human PK or PD data were fitted)",
    n_subjects     = 0L,
    n_studies      = 0L,
    age_range      = "adults (target population of NCT05659264); not otherwise specified",
    weight_range   = "70 kg reference body weight (Table 2 column header)",
    sex_female_pct = NA_real_,
    disease_state  = "Intended population: adult patients with stable heart failure with reduced ejection fraction (first-in-human trial NCT05659264). The projection itself carries no human covariate or disease information.",
    dose_range     = "Simulated doses 0.01-1 mg/kg every 2 weeks; 0.025 mg/kg every 2 weeks was the selected first-in-human starting dose (Figure 5, Discussion).",
    regions        = "Not applicable (simulation only).",
    notes          = "Every parameter value is a deterministic allometric transform of the cynomolgus-monkey estimate in Table 1 using Equation 1, Y = a * (BW_human / BW_cyno)^beta with BW ratio 70 / 2.5 = 28: exponent 1 for all volumes and for Slope, 0.75 for the Rel2-vlk mRNA clearances (CL, CL2, CL3), and 0.85 for the Rel2-vlk protein clearances (Kprot, K50). The 0.85 protein-clearance exponent was chosen from a literature review of therapeutic proteins (paper references 19-21). Inter-individual variances, the effect-compartment rate constant Ke0, and both residual-error SDs are carried over unscaled (Table 2). No within-human body-weight scaling is included because the paper does not fit or report any; the values are for a 70 kg adult and doses were simulated in mg/kg. Because concentrations are in ng/mL and volumes in mL, the dose amount supplied to the model must be expressed in ng (e.g. 0.025 mg/kg x 70 kg = 1.75e6 ng). See the vignette Assumptions and deviations section for the clearance-vs-rate-constant reading of Kprot and K50, and for the small rounding differences between Table 2 as printed and an exact application of Equation 1."
  )

  ini({
    # =====================================================================
    # Rel2-vlk mRNA PK -- Table 2, "Scaled estimate (human 70kg)" column.
    # Structure identical to Kaushal_2024_mRNA0184_cyno (Table S1):
    #   central = plasma-1 (volume V1), peripheral1 = tissue (volume V2),
    #   peripheral2 = plasma-2 (ALSO volume V1, per Figure 1).
    # A single intercompartmental clearance CL3 serves both the
    # tissue -> plasma-2 and plasma-2 -> tissue directions.
    # Scaling checks against Equation 1 with BW ratio 70 / 2.5 = 28:
    #   V1  112  * 28^1    = 3136   -> Table 2 prints 3136 (exact)
    #   CL   258 * 28^0.75 = 3140   -> Table 2 prints 3144 (rounding)
    #   CL2 42.4 * 28^0.75 =  516.1 -> Table 2 prints  516 (exact)
    #   V2   160 * 28^1    = 4480   -> Table 2 prints 4490 (rounding)
    #   CL3 9.86 * 28^0.75 =  120.0 -> Table 2 prints  120 (exact)
    # The printed Table 2 values are used verbatim.
    # =====================================================================
    lvc   <- fixed(log(3136)); label("V1, Rel2-vlk mRNA volume of distribution in plasma, shared by both plasma compartments (mL)") # Table 2; cyno value was fixed from a prior estimation run
    lcl12 <- log(3144);        label("CL, plasma-1 -> tissue intercompartmental clearance of Rel2-vlk mRNA (mL/h)")                  # Table 2
    lcl20 <- log(516);         label("CL2, systemic (tissue) elimination clearance of Rel2-vlk mRNA (mL/h)")                         # Table 2
    lvp   <- log(4490);        label("V2, Rel2-vlk mRNA volume of distribution in tissue / target site (mL)")                        # Table 2
    lcl23 <- log(120);         label("CL3, tissue <-> plasma-2 intercompartmental clearance of Rel2-vlk mRNA, same value both directions (mL/h)") # Table 2

    # IIV carried over unscaled from the cynomolgus fit (Table 2, scaling
    # coefficient "-"). Phoenix NLME reports the Omega diagonal as
    # log-scale VARIANCES; see vignette Assumptions.
    etalcl12 ~ 0.266   # Table 2 etaCL
    etalvp   ~ 1.06    # Table 2 etaV2
    etalcl23 ~ 0.171   # Table 2 etaCL3

    propSd <- 0.557; label("Proportional residual SD on Rel2-vlk mRNA concentration (fraction)") # Table 2, carried over unscaled

    # =====================================================================
    # Rel2-vlk protein PK/PD -- Table 2, "Scaled estimate (human 70kg)".
    #   Ke0    0.193 unscaled                 -> Table 2 prints  0.193
    #   Vc     163  * 28^1    =  4564         -> Table 2 prints  4550 (rounding)
    #   Vp     364  * 28^1    = 10192         -> Table 2 prints 10182 (rounding)
    #   Kprot  25.5 * 28^0.85 =   433.3       -> Table 2 prints   433
    #   K50    4.58 * 28^0.85 =    77.8       -> Table 2 prints    77.7
    #   Slope  0.540 * 28^1   =    15.12      -> Table 2 prints    15.1
    #
    # As in the cynomolgus file, Kprot and K50 are labelled (h^-1) in
    # Table 1 but are implemented here as CLEARANCES (mL/h) acting on
    # concentrations. The 0.85 exponent Table 2 applies to them is the
    # paper's own allometric exponent for the CLEARANCE of therapeutic
    # proteins, and only the clearance reading reproduces the paper's
    # reported human Rel2-vlk protein half-life of "approximately 6-7
    # days" (Discussion): the clearance form gives 143 h = 5.96 days,
    # the literal rate-constant form gives 0.019 h. Full justification
    # is in the vignette Assumptions and deviations section.
    # =====================================================================
    lke0     <- log(0.193);        label("Ke0, equilibration rate constant of the hypothetical effect compartment with plasma Rel2-vlk mRNA (1/h)") # Table 2, carried over unscaled
    lvc_prot <- fixed(log(4550));  label("Vc, Rel2-vlk protein volume of distribution in plasma / central compartment (mL)")                        # Table 2; cyno value was fixed from a prior estimation run
    lvp_prot <- log(10182);        label("Vp, Rel2-vlk protein volume of distribution in tissue / peripheral compartment (mL)")                     # Table 2
    lkprot   <- log(433);          label("Kprot, symmetric central <-> peripheral intercompartmental clearance of Rel2-vlk protein (mL/h)")         # Table 2
    lk50     <- log(77.7);         label("K50, first-order elimination clearance of Rel2-vlk protein from the central compartment (mL/h)")          # Table 2
    lslope   <- log(15.1);         label("Slope, Rel2-vlk protein production per unit effect-compartment mRNA concentration (mL/h)")                # Table 2

    etalk50   ~ 0.073  # Table 2 etaK50
    etalslope ~ 0.417  # Table 2 etaSlope

    propSd_Rel2vlk <- 0.378; label("Proportional residual SD on Rel2-vlk protein concentration (fraction)") # Table 2, carried over unscaled
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters.
    # ------------------------------------------------------------------
    vc    <- exp(lvc)
    cl12  <- exp(lcl12 + etalcl12)
    cl20  <- exp(lcl20)
    vp    <- exp(lvp + etalvp)
    cl23  <- exp(lcl23 + etalcl23)

    ke0   <- exp(lke0)
    vcp   <- exp(lvc_prot)
    vpp   <- exp(lvp_prot)
    kprot <- exp(lkprot)
    k50   <- exp(lk50 + etalk50)
    slope <- exp(lslope + etalslope)

    # ------------------------------------------------------------------
    # Rel2-vlk mRNA: plasma-1 (central) -> tissue (peripheral1) <->
    # plasma-2 (peripheral2); elimination from the tissue compartment.
    # Both plasma compartments share vc (= V1). State expressions are
    # written inline inside d/dt() on purpose.
    # ------------------------------------------------------------------
    d/dt(central)     <- -cl12 * (central / vc)
    d/dt(peripheral1) <-  cl12 * (central / vc) -
      cl23 * (peripheral1 / vp) + cl23 * (peripheral2 / vc) -
      cl20 * (peripheral1 / vp)
    d/dt(peripheral2) <-  cl23 * (peripheral1 / vp) - cl23 * (peripheral2 / vc)

    # ------------------------------------------------------------------
    # Hypothetical effect compartment (carried as a CONCENTRATION, the
    # paper's Ce) driven toward total plasma mRNA C = (A1 + A3) / V1.
    # ------------------------------------------------------------------
    d/dt(effect) <- ke0 * ((central / vc + peripheral2 / vc) - effect)

    # ------------------------------------------------------------------
    # Rel2-vlk protein: production linear in the effect-compartment mRNA
    # concentration, symmetric central <-> peripheral distribution, and
    # first-order elimination from the central compartment.
    # ------------------------------------------------------------------
    d/dt(rel2vlk)   <- slope * effect -
      kprot * (rel2vlk / vcp - rel2vlk_p / vpp) - k50 * (rel2vlk / vcp)
    d/dt(rel2vlk_p) <- kprot * (rel2vlk / vcp - rel2vlk_p / vpp)

    # ------------------------------------------------------------------
    # Observations.
    # ------------------------------------------------------------------
    Cc      <- central / vc + peripheral2 / vc   # Table S1: C = (A1 + A3) / V1
    Rel2vlk <- rel2vlk / vcp                     # Table S1: Rel2-vlk = protein1 / Vc

    Cc      ~ prop(propSd)
    Rel2vlk ~ prop(propSd_Rel2vlk)
  })
}
