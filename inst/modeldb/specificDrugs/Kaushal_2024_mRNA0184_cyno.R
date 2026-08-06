Kaushal_2024_mRNA0184_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Translational semi-mechanistic PK/PD model for mRNA-0184, a lipid-nanoparticle-encapsulated mRNA encoding human relaxin-2 fused to a variable light chain kappa domain (Rel2-vlk), for heart failure. PK: 3-compartment plasma1-tissue-plasma2 redistribution of Rel2-vlk mRNA (both plasma compartments share V1; elimination is from the tissue compartment) reproducing the delayed second concentration peak characteristic of LNP modalities; observed mRNA concentration is the sum of the two plasma compartments. PD: Rel2-vlk protein production linear in the concentration of a hypothetical effect compartment equilibrating with plasma mRNA, distributed in a 2-compartment (central/peripheral) system with a single symmetric intercompartmental clearance and first-order elimination from the central compartment."
  reference <- paste(
    "Kaushal N, Attarwala H, Iqbal MJ, Saini R, Van L, Liang M.",
    "Translational pharmacokinetic/pharmacodynamic model for mRNA-0184,",
    "an investigational therapeutic for the treatment of heart failure.",
    "Clin Transl Sci. 2024;17(8):e13894. doi:10.1111/cts.13894.",
    "Structural equations from Supporting Information Table S1;",
    "parameter estimates from Table 1.",
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
    species       = "cynomolgus monkey (male, healthy)",
    n_subjects    = 12L,
    n_studies     = 1L,
    age_range     = "not reported",
    weight_range  = "2.5 kg reference body weight (Table 2 column header)",
    sex_female_pct = 0,
    disease_state = "Healthy (non-disease) non-human primates. The efficacious-exposure anchor used for human dose selection (AUEC over a 2-week interval = 486 ng/mL*h after weekly 0.15 mg/kg dosing) came from a SEPARATE study in aged, high-fat-diet-induced obese cynomolgus monkeys with naturally developed cardiovascular and metabolic disease; that study is reported as 'data not shown' and is NOT part of this model's fitted dataset.",
    dose_range    = "Single dose of mRNA-0184 at 0.15, 0.5, or 1 mg/kg administered as a 1-h intravenous infusion (N = 4 per dose level).",
    regions       = "Preclinical study sponsored by Moderna, Inc. (Cambridge, MA, USA).",
    notes         = "Rel2-vlk mRNA and Rel2-vlk protein plasma concentrations were measured at various timepoints up to 337 h post-dose. Assays: pre-qualified bDNA (mRNA) and ELISA (protein); LLOQ 0.050 ng/mL for Rel2-vlk mRNA and 20 pg/mL for Rel2-vlk protein (Methods, Data used for model development). Estimation was performed in Phoenix NLME 8.3.4.295 using first-order conditional estimation-extended least squares. V1 (mRNA plasma volume) and Vc (protein central volume) were fixed based on a prior estimation run (Table 1 footnote). Doses are stated in mg/kg; because concentrations are in ng/mL and volumes in mL, the dose amount supplied to the model must be expressed in ng (e.g. 0.15 mg/kg x 2.5 kg = 3.75e5 ng). See the vignette Assumptions and deviations section for the clearance-vs-rate-constant reading of Kprot and K50."
  )

  ini({
    # =====================================================================
    # Rel2-vlk mRNA PK (Table 1, "PK model parameter estimates"; structure
    # from Supporting Information Table S1).
    #
    # Table S1 writes the system in rate-constant form on amounts:
    #   dA1/dt = Input - K12 * A1
    #   dA2/dt = K12 * A1 - K23 * A2 + K32 * A3 - K20 * A2
    #   dA3/dt = K23 * A2 - K32 * A3
    #   C      = (A1 + A3) / V1
    # with K12 = CL/V1, K23 = CL3/V2, K32 = CL3/V1, K20 = CL2/V2.
    # Substituting those definitions gives the equivalent clearance x
    # concentration form implemented in model():
    #   A1 = plasma-1 (central, volume V1), A2 = tissue (peripheral1,
    #   volume V2), A3 = plasma-2 (peripheral2, ALSO volume V1 -- see
    #   Figure 1, where the third circle is labelled "Plasma V1").
    # A single intercompartmental clearance CL3 serves BOTH tissue ->
    # plasma-2 (K23) and plasma-2 -> tissue (K32) directions.
    # =====================================================================
    lvc   <- fixed(log(112)); label("V1, Rel2-vlk mRNA volume of distribution in plasma, shared by both plasma compartments (mL)") # Table 1, fixed from a prior estimation run
    lcl12 <- log(258);        label("CL, plasma-1 -> tissue intercompartmental clearance of Rel2-vlk mRNA (mL/h)")                  # Table 1 tvCL, SE 25.3%
    lcl20 <- log(42.4);       label("CL2, systemic (tissue) elimination clearance of Rel2-vlk mRNA (mL/h)")                         # Table 1 tvCL2, SE 5.85%
    lvp   <- log(160);        label("V2, Rel2-vlk mRNA volume of distribution in tissue / target site (mL)")                        # Table 1 tvV2, SE 34.1%
    lcl23 <- log(9.86);       label("CL3, tissue <-> plasma-2 intercompartmental clearance of Rel2-vlk mRNA, same value both directions (mL/h)") # Table 1 tvCL3, SE 1.22%

    # IIV on the mRNA PK parameters (Table 1 eta rows). Phoenix NLME
    # reports the Omega diagonal as VARIANCES on the log scale, which is
    # the scale nlmixr2 expects here; see vignette Assumptions.
    etalcl12 ~ 0.266   # Table 1 etaCL,  SE 0.106, shrinkage 36.8%
    etalvp   ~ 1.06    # Table 1 etaV2,  SE 0.253, shrinkage 23.9%
    etalcl23 ~ 0.171   # Table 1 etaCL3, SE 0.034, shrinkage 16.6%

    # Proportional residual SD on Rel2-vlk mRNA (Table 1, SE 0.080;
    # eps-shrinkage 12.3%).
    propSd <- 0.557; label("Proportional residual SD on Rel2-vlk mRNA concentration (fraction)") # Table 1

    # =====================================================================
    # Rel2-vlk protein PK/PD (Table 1, "mRNA-0184-Rel2-vlk PK/PD model
    # parameter estimates"; structure from Table S1).
    #
    # Table S1:
    #   dCe/dt        = Ke0 * (C(t) - Ce(t))
    #   dprotein1/dt  = slope * Ce - Kprot * [protein1 - protein2] - K50 * protein1
    #   dprotein2/dt  = Kprot * (protein1 - protein2)
    #   Rel2-vlk      = protein1 / Vc
    #
    # Kprot and K50 are labelled (h^-1) in Table 1 but are implemented
    # here as CLEARANCES (mL/h) acting on concentrations, i.e.
    #   Kprot * (protein1/Vc - protein2/Vp)   and   K50 * protein1/Vc,
    # exactly the substitution Table S1 itself performs for the mRNA
    # rate constants. Four independent checks select this reading; see
    # the vignette Assumptions and deviations section:
    #   (1) Vp is estimated (SE 42.7%) and allometrically scaled in
    #       Table 2, but never appears in the printed equations -- it is
    #       only identifiable in the clearance form.
    #   (2) Table 2 scales Kprot and K50 with the exponent 0.85, which
    #       the paper describes as the allometric exponent for
    #       CLEARANCE of therapeutic proteins; a rate constant would
    #       scale with 0.85 - 1 = -0.15.
    #   (3) The clearance reading gives a terminal protein half-life of
    #       143 h (5.96 days) for the human-scaled parameters, matching
    #       the paper's "approximately 6-7 days" (Discussion); the
    #       literal rate-constant reading gives 0.019 h.
    #   (4) The single symbol Kprot for both directions is only
    #       self-consistent as a symmetric intercompartmental clearance,
    #       mirroring CL3 on the mRNA side and the paper's assumption
    #       that "the distribution rate constant of Rel2-vlk protein is
    #       the same between the plasma and tissue compartments".
    # =====================================================================
    lke0     <- log(0.193);        label("Ke0, equilibration rate constant of the hypothetical effect compartment with plasma Rel2-vlk mRNA (1/h)") # Table 1 tvKe0, SE 0.050
    lvc_prot <- fixed(log(163));   label("Vc, Rel2-vlk protein volume of distribution in plasma / central compartment (mL)")                        # Table 1 tvVc, fixed from a prior estimation run
    lvp_prot <- log(364);          label("Vp, Rel2-vlk protein volume of distribution in tissue / peripheral compartment (mL)")                     # Table 1 tvVp, SE 42.7%
    lkprot   <- log(25.5);         label("Kprot, symmetric central <-> peripheral intercompartmental clearance of Rel2-vlk protein (mL/h)")         # Table 1 tvKprot, SE 5.21%
    lk50     <- log(4.58);         label("K50, first-order elimination clearance of Rel2-vlk protein from the central compartment (mL/h)")          # Table 1 tvK50, SE 0.491
    lslope   <- log(0.540);        label("Slope, Rel2-vlk protein production per unit effect-compartment mRNA concentration (mL/h)")                # Table 1 tvSlope, SE 0.106

    # IIV on the PD parameters (Table 1 eta rows, variances as above).
    etalk50   ~ 0.073  # Table 1 etaK50,   SE 0.065, shrinkage 14.5%
    etalslope ~ 0.417  # Table 1 etaSlope, SE 0.160, shrinkage 2.76%

    # Proportional residual SD on Rel2-vlk protein (Table 1, SE 0.027;
    # eps-shrinkage 9.02%).
    propSd_Rel2vlk <- 0.378; label("Proportional residual SD on Rel2-vlk protein concentration (fraction)") # Table 1
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
    # paper's Ce) driven toward the total plasma mRNA concentration
    # C = (A1 + A3) / V1.
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
