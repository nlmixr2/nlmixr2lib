Falkenhagen_2023_warfarin_qsp <- function() {
  description <- paste(
    "QSP-derived mechanism-based warfarin/INR pharmacodynamic model (Falkenhagen 2023).",
    "Obtained by systematically reducing the 62-ODE / 174-parameter Wajima 2009 blood-coagulation",
    "quantitative systems pharmacology model down to 6 ODEs and 11 structural parameters, while",
    "guaranteeing under 10% relative INR error for at least 95% of a diverse virtual population.",
    "One-compartment oral warfarin PK inhibits vitamin K hydroquinone (VKH2) synthesis through an",
    "Imax function; VKH2 in turn drives the synthesis of coagulation Factors II, VII, and X, each",
    "modelled as a turnover pool holding its own pre-stimulus steady state. The INR is recovered",
    "algebraically as a power law in the product of the three relative factor concentrations,",
    "INR = INR0 * (II/II0 * VII/VII0 * X/X0)^gamma with gamma = -0.1975. CYP2C9 *1/*2/*3 allele",
    "counts set warfarin clearance and VKORC1 -1639 G/A allele counts set IC50, both as per-allele",
    "sums. All parameter values are the Wajima 2009 reference parameterization carried through the",
    "reduction; no parameter was estimated from clinical data in this paper."
  )
  reference <- paste(
    "Falkenhagen U, Knoechel J, Kloft C, Huisinga W.",
    "Deriving mechanism-based pharmacodynamic models by reducing quantitative systems pharmacology",
    "models: An application to warfarin.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(4):432-443. doi:10.1002/psp4.12903.",
    "Reduced in vivo ODE system from Supplementary Material Section S5 (Eqs S27-S29);",
    "INR power law from Equation (13); genotype parameters from Supplementary Table S1;",
    "steady-state parameter interdependencies from Supplementary Material Section S2 (Eq S3).",
    "Numeric values of the reference parameterization are inherited from the underlying QSP model:",
    "Wajima T, Isbister GK, Duffull SB. A comprehensive model for the humoral coagulation network",
    "in humans. Clin Pharmacol Ther. 2009;86(3):290-298. doi:10.1038/clpt.2009.87,",
    "as implemented in the authors' own supplemental code archive doi:10.5281/zenodo.7417886."
  )
  vignette <- "Falkenhagen_2023_warfarin_qsp"
  paper_specific_compartments <- c("vkh2", "factor_ii", "factor_vii", "factor_x")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 (wild-type) alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "CYP2C9_S1_COUNT + CYP2C9_S2_COUNT + CYP2C9_S3_COUNT = 2. Each *1 allele contributes",
        "0.1000 L/h to warfarin CL (Falkenhagen 2023 Supplementary Table S1). The per-allele",
        "values are anchored so that the *1/*1 wild-type CL equals the Wajima 2009 reference",
        "CL of 0.2 L/h. Allele frequency in the reduction population: 0.815 (Table S1)."
      ),
      source_name        = "CYP2C9 genotype (allele pair a/b; CL^{a/b} = CL^a + CL^b, Eq S1)"
    ),
    CYP2C9_S2_COUNT = list(
      description        = "Count of CYP2C9*2 alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Each *2 allele contributes 0.0505 L/h to warfarin CL (Falkenhagen 2023 Supplementary",
        "Table S1). Allele frequency in the reduction population: 0.112 (Table S1)."
      ),
      source_name        = "CYP2C9 genotype (allele pair a/b; CL^{a/b} = CL^a + CL^b, Eq S1)"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Each *3 allele contributes 0.0243 L/h to warfarin CL (Falkenhagen 2023 Supplementary",
        "Table S1). Allele frequency in the reduction population: 0.073 (Table S1)."
      ),
      source_name        = "CYP2C9 genotype (allele pair a/b; CL^{a/b} = CL^a + CL^b, Eq S1)"
    ),
    VKORC1_1639G_COUNT = list(
      description        = paste(
        "Count of VKORC1 -1639G alleles per subject (0, 1, or 2).",
        "The complementary -1639A count is 2 - VKORC1_1639G_COUNT."
      ),
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Each -1639G allele contributes 0.1700 mg/L and each -1639A allele 0.0796 mg/L to the",
        "warfarin IC50 (Falkenhagen 2023 Supplementary Table S1, Eq S1). The per-allele values",
        "are anchored so that the GG genotype IC50 equals the Wajima 2009 reference of 0.34 mg/L.",
        "Allele frequencies in the reduction population: G 0.608, A 0.392 (Table S1), matching the",
        "Warfarin Genetics (WARG) study frequencies reported by Hamberg 2010."
      ),
      source_name        = "VKORC1 genotype (allele pair a/b; IC50^{ab} = IC50^a + IC50^b, Eq S1)"
    )
  )

  population <- list(
    species        = "human (in silico virtual population; no clinical data were fitted)",
    n_subjects     = 1000L,
    n_studies      = 0L,
    disease_state  = paste(
      "Virtual adults on chronic oral warfarin anticoagulation. The model was NOT fitted to",
      "patient data: it is a systematic mathematical reduction of the Wajima 2009 blood-coagulation",
      "QSP model, and the virtual population exists to define the parameter region over which the",
      "reduction is guaranteed accurate."
    ),
    dose_range     = paste(
      "4 mg orally once daily for 30 days; reduced to 1 mg once daily for virtual individuals whose",
      "steady-state INR under the 4 mg regimen would exceed 4, so as to keep INR in a clinically",
      "relevant range (Falkenhagen 2023, Workflow of reducing the scenarios considering IIV)."
    ),
    cyp2c9_freq    = "CYP2C9 allele frequencies *1 0.815, *2 0.112, *3 0.073 (Supplementary Table S1)",
    vkorc1_freq    = "VKORC1 -1639 allele frequencies G 0.608, A 0.392 (Supplementary Table S1)",
    iiv_assumption = paste(
      "Unexplained random IIV was assumed (not estimated) as an independent log-normal distribution",
      "on every parameter and initial value with a 40% coefficient of variation,",
      "q_j ~ logN(log(q_ref_j), 0.4^2) (Equation 3 / Supplementary Eq S2). Latin hypercube sampling",
      "was used to draw the 1000 virtual individuals. Parameter variability was assumed uncorrelated",
      "for lack of knowledge on the correlation structure (Discussion)."
    ),
    validation     = paste(
      "The reduced model reproduced the full QSP model's INR to within 10% relative error for more",
      "than 99% of the 1000 virtual individuals; the 5 excluded individuals had errors between 10%",
      "and 13% (Results / Discussion)."
    ),
    notes          = paste(
      "Genotypes were assigned deterministically so that allele frequencies matched those reported",
      "for the Warfarin Genetics (WARG) study in Hamberg 2010. The INR range of validity is below 4,",
      "and the model presumes the standard prothrombin-time test (Discussion)."
    )
  )

  ini({
    # =================================================================
    # Warfarin PK. Reference parameterization from the Wajima 2009 QSP
    # model, carried through the reduction unchanged. Falkenhagen 2023
    # Supplementary Eq S27 writes the PK on the concentration scale:
    #   dA_warf/dt = -ka_warf * A_warf
    #   dC_warf/dt = (ka_warf / Vd_warf) * A_warf - ke_warf * C_warf
    # model() below multiplies the second equation through by Vd_warf so
    # that `central` carries an AMOUNT (the nlmixr2 convention) and
    # Cc <- central / vc; the two forms are algebraically identical.
    # =================================================================
    lka <- fixed(log(1.0)); label("Warfarin absorption rate constant ka_warf (1/h)")  # Wajima 2009 reference value ka_Warf = 1.0 1/h (Zenodo code archive doi:10.5281/zenodo.7417886, Wajima2009BloodCoagulation_parameters.m); Supplementary S6 notes ka_warf is reasonable to fix
    lvc <- fixed(log(10));  label("Warfarin apparent volume of distribution Vd_warf (L)")  # Wajima 2009 reference value Vd_Warf = 10 L (Zenodo code archive, Wajima2009BloodCoagulation_parameters.m)

    # Supplementary Section S6 lists k_e,warf (not CL) as one of the 11
    # structural parameters, and the authors' population code applies the
    # CYP2C9 genotype as a multiplicative ratio on k_e,warf. The reference
    # rate below is the CYP2C9*1/*1 value, ke = CL/Vd = 0.2 / 10 = 0.02 1/h.
    lkel <- fixed(log(0.02)); label("Warfarin elimination rate constant ke_warf for the CYP2C9*1/*1 reference (1/h)")  # Wajima 2009 reference: Cl_Warf = 0.2 L/h and Vd_Warf = 10 L give ke_Warf = 0.02 1/h (Zenodo code archive, Wajima2009BloodCoagulation_parameters.m line par(I.ke_Warf) = par(I.Cl_Warf)/par(I.Vd_Warf))

    # Per-CYP2C9-allele CL contributions; subject CL is the sum over the
    # two alleles, CL^{a/b} = CL^a + CL^b (Supplementary Eq S1). The
    # *1/*1 sum, 2 * 0.1000 = 0.2 L/h, reproduces the Wajima 2009
    # reference clearance exactly, so model() applies the genotype as the
    # ratio cl_typ / cl_ref, which is 1 for the *1/*1 reference.
    lcl_cyp2c9_s1 <- fixed(log(0.1000)); label("Warfarin CL per CYP2C9*1 allele (L/h)")  # Falkenhagen 2023 Supplementary Table S1
    lcl_cyp2c9_s2 <- fixed(log(0.0505)); label("Warfarin CL per CYP2C9*2 allele (L/h)")  # Falkenhagen 2023 Supplementary Table S1
    lcl_cyp2c9_s3 <- fixed(log(0.0243)); label("Warfarin CL per CYP2C9*3 allele (L/h)")  # Falkenhagen 2023 Supplementary Table S1

    # =================================================================
    # Warfarin effect on the vitamin K cycle. Warfarin inhibits VKORC1
    # and thereby VKH2 synthesis through an Imax function of the
    # warfarin concentration (Supplementary Eq S27, inhibition term
    # 1 - Imax * C_warf / (IC50 + C_warf)).
    # Per-VKORC1-allele IC50 contributions sum over the two alleles,
    # IC50^{ab} = IC50^a + IC50^b (Supplementary Eq S1). The GG sum,
    # 2 * 0.1700 = 0.34 mg/L, reproduces the Wajima 2009 reference IC50.
    # =================================================================
    lic50_vkorc1_g <- fixed(log(0.1700)); label("Warfarin IC50 per VKORC1 -1639G allele (mg/L)")  # Falkenhagen 2023 Supplementary Table S1
    lic50_vkorc1_a <- fixed(log(0.0796)); label("Warfarin IC50 per VKORC1 -1639A allele (mg/L)")  # Falkenhagen 2023 Supplementary Table S1
    limax          <- fixed(log(1));      label("Maximum fractional inhibition of VKH2 synthesis Imax (unitless)")  # Wajima 2009 reference value lmax = 1 (Zenodo code archive); Supplementary S6 notes Imax is reasonable to fix, and the authors' population code excludes Imax from random IIV

    # =================================================================
    # Turnover rate constants (1/h). deg_VK2 is the vitamin-K to VKH2
    # conversion rate; deg_II / deg_VII / deg_X are the first-order
    # degradation rates of the three warfarin-dependent coagulation
    # factors retained by the reduction. The corresponding half-lives
    # log(2)/deg are 30.4 h (VK2), 69.3 h (II), 5.8 h (VII) and 38.5 h
    # (X); the main text quotes ~69 h, ~6 h and ~39 h for II, VII and X
    # respectively (Biomarker proposal section), independently
    # corroborating the three factor rate constants.
    # =================================================================
    lkdeg_vk2 <- fixed(log(0.0228)); label("Vitamin K to VKH2 conversion rate deg_VK2 (1/h)")  # Wajima 2009 reference value degVK2 = 0.0228 1/h (Zenodo code archive)
    lkdeg_ii  <- fixed(log(0.010));  label("Factor II degradation rate deg_II (1/h)")   # Wajima 2009 reference value degII = 0.010 1/h (Zenodo code archive); half-life 69.3 h matches the ~69 h quoted in Falkenhagen 2023 Results
    lkdeg_vii <- fixed(log(0.12));   label("Factor VII degradation rate deg_VII (1/h)")  # Wajima 2009 reference value degVII = 0.12 1/h (Zenodo code archive); half-life 5.8 h matches the ~6 h quoted in Falkenhagen 2023 Results
    lkdeg_x   <- fixed(log(0.018));  label("Factor X degradation rate deg_X (1/h)")     # Wajima 2009 reference value degX = 0.018 1/h (Zenodo code archive); half-life 38.5 h matches the ~39 h quoted in Falkenhagen 2023 Results

    # =================================================================
    # INR power law (Equation 13):
    #   INR(t)/INR(0) = (II(t)/II(0) * VII(t)/VII(0) * X(t)/X(0))^gamma
    # gamma was obtained as the slope of a linear fit in log-log space
    # over the clinically relevant INR range 2-3 (Figure 4). It is
    # negative, so it is NOT log-transformed.
    # =================================================================
    gamma <- fixed(-0.1975); label("Power-law exponent linking the relative coagulation-factor product to the INR (unitless)")  # Falkenhagen 2023 Equation (13) / Figure 4
    le0   <- fixed(log(1));  label("Pre-treatment (drug-free) baseline INR (unitless)")  # INR is defined as PT/PT_ref (Equation 1), so the untreated reference individual has INR = 1; Figure 6 INR-time profiles start at 1

    # =================================================================
    # Pre-stimulus (drug-free) steady-state concentrations. VK and VKH2
    # are in the QSP model's relative vitamin-K units; the coagulation
    # factors are in nmol/L. Supplementary Section S2 (Eq S3) requires
    # the synthesis rates to satisfy k_syn,j = k_deg,j * [x0]_j so that
    # the pre-stimulus steady state is sustained; model() therefore
    # DERIVES deg_VKH2 and the three factor synthesis terms from these
    # baselines rather than carrying them as free parameters.
    # =================================================================
    bl_vk   <- fixed(1.0);    label("Pre-stimulus vitamin K concentration VK(0) (relative units)")  # Wajima 2009 reference initial value X0(VK) = 1.0 (Zenodo code archive, Wajima2009BloodCoagulation_initialvalues.m)
    bl_vkh2 <- fixed(0.1);    label("Pre-stimulus vitamin K hydroquinone concentration VKH2(0) (relative units)")  # Wajima 2009 reference initial value X0(VKH2) = 0.1 (Zenodo code archive)
    bl_ii   <- fixed(1394.4); label("Pre-stimulus Factor II concentration II(0) (nmol/L)")   # Wajima 2009 reference initial value X0(II) = 1394.4 nmol/L (Zenodo code archive)
    bl_vii  <- fixed(10.0);   label("Pre-stimulus Factor VII concentration VII(0) (nmol/L)")  # Wajima 2009 reference initial value X0(VII) = 10.0 nmol/L (Zenodo code archive)
    bl_x    <- fixed(174.3);  label("Pre-stimulus Factor X concentration X(0) (nmol/L)")     # Wajima 2009 reference initial value X0(X) = 174.3 nmol/L (Zenodo code archive)

    # =================================================================
    # Assumed (not estimated) inter-individual variability. Equation (3)
    # / Supplementary Eq S2 place an independent log-normal distribution
    # with 40% CV on every parameter and initial value:
    #   q_j ~ logN(log(q_ref_j), 0.4^2)   =>  omega^2 = 0.4^2 = 0.16
    # All variances are fixed, not estimated -- no clinical data were
    # fitted. Two parameters are deliberately excluded from random IIV
    # in the authors' population code: Imax (held at 1) and the dose.
    # deg_VKH2 and the factor synthesis rates carry no independent eta
    # because they are derived from the sampled parameters and baselines
    # to preserve the pre-stimulus steady state (Supplementary S2).
    # =================================================================
    etalka       ~ fixed(0.16)
    etalvc       ~ fixed(0.16)
    etalkel      ~ fixed(0.16)
    etalic50     ~ fixed(0.16)
    etalkdeg_vk2 ~ fixed(0.16)
    etalkdeg_ii  ~ fixed(0.16)
    etalkdeg_vii ~ fixed(0.16)
    etalkdeg_x   ~ fixed(0.16)
    etabl_vk     ~ fixed(0.16)
    etabl_vkh2   ~ fixed(0.16)
    etabl_ii     ~ fixed(0.16)
    etabl_vii    ~ fixed(0.16)
    etabl_x      ~ fixed(0.16)

    # No residual error is reported: the model was never fitted to
    # observed INR data. Fixed at 0 per the standing convention for
    # unreported RUV (documented in the vignette Errata).
    addSd <- fixed(0); label("Additive residual error SD on INR (unitless)")
  })

  model({
    # -----------------------------------------------------------------
    # Individual warfarin PK parameters. The authors' population code
    # places independent log-normal IIV on ka_warf, Vd_warf and
    # ke_warf; the genotype acts as a multiplicative factor on ke_warf.
    # Reproduced here by forming the genotype-dependent typical
    # elimination rate from the TYPICAL volume (so the CYP2C9 ratio is
    # preserved exactly) and then applying the ke_warf eta.
    # -----------------------------------------------------------------
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)

    # Typical CL as the per-allele sum over the two CYP2C9 alleles
    # (Supplementary Eq S1 / Table S1), expressed as a ratio to the
    # CYP2C9*1/*1 reference so it scales the reference ke_warf.
    cl_typ <- CYP2C9_S1_COUNT * exp(lcl_cyp2c9_s1) +
      CYP2C9_S2_COUNT * exp(lcl_cyp2c9_s2) +
      CYP2C9_S3_COUNT * exp(lcl_cyp2c9_s3)
    cl_ref <- 2 * exp(lcl_cyp2c9_s1)
    kel <- exp(lkel + etalkel) * (cl_typ / cl_ref)

    # Typical IC50 as the per-allele sum over the two VKORC1 alleles
    # (Supplementary Eq S1 / Table S1).
    ic50 <- (VKORC1_1639G_COUNT * exp(lic50_vkorc1_g) +
      (2 - VKORC1_1639G_COUNT) * exp(lic50_vkorc1_a)) * exp(etalic50)

    imax <- exp(limax)

    # -----------------------------------------------------------------
    # Individual turnover rate constants and pre-stimulus baselines.
    # -----------------------------------------------------------------
    kdeg_vk2 <- exp(lkdeg_vk2 + etalkdeg_vk2)
    kdeg_ii  <- exp(lkdeg_ii + etalkdeg_ii)
    kdeg_vii <- exp(lkdeg_vii + etalkdeg_vii)
    kdeg_x   <- exp(lkdeg_x + etalkdeg_x)

    vk_base   <- bl_vk * exp(etabl_vk)
    vkh2_base <- bl_vkh2 * exp(etabl_vkh2)
    ii_base   <- bl_ii * exp(etabl_ii)
    vii_base  <- bl_vii * exp(etabl_vii)
    x_base    <- bl_x * exp(etabl_x)

    # deg_VKH2 is NOT a free parameter: Supplementary Section S2 fixes it
    # by the requirement that the vitamin K cycle be at steady state
    # before dosing, deg_VKH2 = deg_VK2 * VK(0) / VKH2(0). With the
    # reference values this gives 0.0228 * 1.0 / 0.1 = 0.228 1/h.
    kdeg_vkh2 <- kdeg_vk2 * vk_base / vkh2_base

    # -----------------------------------------------------------------
    # One-compartment oral warfarin PK (Supplementary Eq S27, rewritten
    # on the amount scale as noted in ini()).
    # -----------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # Fractional inhibition of VKH2 synthesis by warfarin. With Imax
    # fixed at 1 and no IIV on it, this is bounded in [0, 1).
    inh <- imax * Cc / (ic50 + Cc)

    # -----------------------------------------------------------------
    # Vitamin K hydroquinone. Vitamin K itself is an ENVIRONMENTAL state
    # in the reduced model: the reduction found that holding it constant
    # at its initial value VK(0) preserves the response, so it appears
    # here as the constant `vk_base` rather than as an ODE (Figure 5,
    # where vitamin K is marked "*"; Supplementary Eq S27).
    # -----------------------------------------------------------------
    d/dt(vkh2) <- kdeg_vk2 * vk_base * (1 - inh) - kdeg_vkh2 * vkh2

    # -----------------------------------------------------------------
    # Warfarin-dependent coagulation factors (Supplementary Eqs S27-S29).
    # Each factor is a turnover pool whose synthesis is proportional to
    # the relative VKH2 concentration, with the synthesis constant fixed
    # by the pre-stimulus steady state (Supplementary S2 Eq S3):
    #   d/dt(F) = deg_F * F(0) * VKH2/VKH2(0) - deg_F * F
    # so that d/dt(F) = 0 when VKH2 = VKH2(0) and F = F(0).
    # -----------------------------------------------------------------
    d/dt(factor_ii)  <- kdeg_ii * ii_base * (vkh2 / vkh2_base) - kdeg_ii * factor_ii
    d/dt(factor_vii) <- kdeg_vii * vii_base * (vkh2 / vkh2_base) - kdeg_vii * factor_vii
    d/dt(factor_x)   <- kdeg_x * x_base * (vkh2 / vkh2_base) - kdeg_x * factor_x

    vkh2(0)       <- vkh2_base
    factor_ii(0)  <- ii_base
    factor_vii(0) <- vii_base
    factor_x(0)   <- x_base

    # -----------------------------------------------------------------
    # Relative factor concentrations and the algebraic INR (Equation 13).
    # The reduction showed that the in vitro prothrombin-time scenario
    # collapses to a dependence on the PRODUCT of the three relative
    # factor concentrations only (Equation 12), which is what makes this
    # single power law sufficient.
    # -----------------------------------------------------------------
    rel_ii  <- factor_ii / ii_base
    rel_vii <- factor_vii / vii_base
    rel_x   <- factor_x / x_base

    INR <- exp(le0) * (rel_ii * rel_vii * rel_x)^gamma
    INR ~ add(addSd)
  })
}
