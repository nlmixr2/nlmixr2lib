Ayyar_2024_givosiran <- function() {
  description <- "Mechanistic translational PK model for the GalNAc-siRNA givosiran (Ayyar & Song 2024) parameterized for human (70 kg adult). 22-ODE system covering SC depot, central plasma (parent + AS(N-1)3' active metabolite), competitive ASGPR receptor binding (free target, parent-target complex, metabolite-target complex), receptor-mediated hepatocyte internalization, endolysosomal sequestration / degradation / endosomal escape, free cytoplasmic siRNA, RISC-loaded siRNA (combined parent + metabolite), kidney vascular and tissue distribution with a deep bound pool and GFR elimination - for parent and metabolite. Pharmacodynamic ALAS1 mRNA silencing (rat-only in the paper) is not included in the human parameterization."
  reference <- "Ayyar VS, Song D. Mechanistic Pharmacokinetics and Pharmacodynamics of GalNAc-siRNA: Translational Model Involving Competitive Receptor-Mediated Disposition and RISC-Dependent Gene Silencing Applied to Givosiran. J Pharm Sci. 2024;113(1):176-190. doi:10.1016/j.xphs.2023.10.026"
  vignette <- "Ayyar_2024_givosiran"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adult",
    weight_range   = "70 kg reference adult (model parameters scaled / calibrated to a 70 kg human)",
    sex_female_pct = NA_real_,
    species        = "Human (70 kg). The paper additionally fits a rat (0.25 kg) parameterization including the ALAS1 mRNA silencing PD layer and a cynomolgus monkey parameterization; both are summarized in the assumptions section of the validation vignette but are not packaged in this file.",
    disease_state  = "Adults with acute hepatic porphyria (AHP) - validation cohort is the published Phase 1 SC dose-ranging study (Ayyar 2024 reference 26).",
    dose_range     = "0.35, 1, 2.5, 5 mg/kg SC single dose (validation cohort, Fig 6). Model is structurally suitable for IV bolus / short infusion as well; SC bioavailability fixed F = 0.9.",
    regions        = "Multi-source - rat and monkey datasets digitized from published literature; human single-dose PK from a published Phase 1 study (ENVISION program for givosiran).",
    notes          = "Pooled time-course data fit with the Nelder-Mead simplex algorithm (without inter-individual variability) in Monolix 2020R1; final human predictions in Fig 6 use Monte Carlo simulations with assumed 20% CV on ka, Vc, and kint (Ayyar 2024 p.184-185, Methods). The packaged model retains those 20% CV IIV terms and otherwise operates as a typical-value mechanistic simulator.",
    scope_note     = "Mechanistic translational PK model: parameters are mechanistic constants. Human values are scaled allometrically from rat (or monkey) and calibrated within < ~20% to match observed Phase 1 plasma profiles. No residual error is fitted to the human cohort - the published prediction interval (Fig 6 shaded region) is generated from the IIV terms alone. PD (ALAS1 mRNA silencing in liver) was characterized in rat only and is not included in this human-parameter file."
  )

  ini({
    # =============================================================
    # ASGPR (asialoglycoprotein receptor) - Table 1, ASGPR section
    # Human values, scaled / calibrated as noted; see footnotes a, c.
    # =============================================================
    lrtot   <- log(340);     label("Baseline ASGPR density in liver (nM)")                     # Table 1, human: 340 (footnote a, see Ayyar 2021)
    lkdeg   <- log(1.9);     label("ASGPR turnover (degradation) rate constant (1/h)")          # Table 1, human (calibrated 1.9; initial scale-up 2.4); footnote b: kdeg = kint
    lkoff   <- log(1.32);    label("Tris-GalNAc-ASGPR dissociation rate constant (1/h)")        # Table 1, shared across species: 1.32 (Ayyar 2021)
    lkint   <- log(1.9);     label("Tris-GalNAc-ASGPR complex internalization rate (1/h)")      # Table 1, human (calibrated 1.9; initial scale-up 2.4)
    lkd     <- log(27.7);    label("Tris-GalNAc-ASGPR equilibrium dissociation constant (nM)")  # Table 1, shared rat/human: 27.7 (30% RSE in rat fit)

    # =============================================================
    # Central compartment - Table 1, Central compartment section
    # =============================================================
    lka     <- log(0.036);   label("SC absorption rate constant (1/h)")                         # Table 1, human (calibrated, footnote c; initial allometric 0.03 with babs=-0.27)
    lfdepot <- log(0.9);     label("SC bioavailability (fraction)")                             # Table 1, fixed in all species at 0.9
    lvc     <- log(2.8);     label("Central volume of distribution (L)")                        # Table 1, human (calibrated 2.8 ~ human plasma volume; initial scale-up 3.13)
    lclm    <- log(3.7);     label("Metabolite (AS(N-1)3') formation / elimination clearance (L/h)") # Table 1, human (footnote c, calibrated to human plasma PK)

    # =============================================================
    # Kidney - Table 1, Kidney section
    # All physiological values fixed; no estimation.
    # =============================================================
    lqkid   <- log(23.5);    label("Renal plasma flow Q_kid (L/h)")                             # Table 1, human (physiological)
    lvk     <- log(0.296);   label("Kidney interstitial and cellular volume V_k (L)")           # Table 1, human (physiological)
    lvkvas  <- log(0.0182);  label("Kidney plasma vascular volume V_k_vas (L)")                 # Table 1, human (physiological)
    lgfr    <- log(7.20);    label("Glomerular filtration rate GFR (L/h)")                      # Table 1, human (physiological)
    lfuk    <- log(0.04);    label("Kidney free fraction f_u_k (unitless)")                     # Table 1, shared (35% RSE rat fit)
    lfugfr  <- log(1);       label("Fraction renally filtered f_u_GFR (unitless)")              # Table 1, footnote g: free non-ASGPR-bound siRNA assumed subject to GFR
    lkassk  <- log(2.9);     label("Kidney binding association constant k_ass_k (1/h)")         # Table 1, shared rat/human: 2.9 (18% RSE)
    lkdisk  <- log(0.004);   label("Kidney binding dissociation constant k_dis_k (1/h)")        # Table 1, shared rat/human: 0.004 (8% RSE)

    # =============================================================
    # Liver (hepatocyte) - Table 1, Liver section
    # =============================================================
    lkassl   <- log(0.00035); label("Liver binding association constant k_ass_l (1/h)")          # Table 1, shared rat/human: 0.00035 (58% RSE)
    lkdisl   <- log(0.0036);  label("Liver binding dissociation constant k_dis_l (1/h)")         # Table 1, shared rat/human: 0.0036 (29% RSE)
    lvhep    <- log(0.69);    label("Liver hepatocyte cellular volume V_hep (L)")                # Table 1, human (allometric exponent b_vol=1.0; rat 0.00272 L)
    lfesc    <- log(0.01);    label("Escaped fraction of endosomal siRNA f_esc (unitless)")      # Table 1, fixed at 0.01 (1%); paper text p.180
    lkdege   <- log(0.007);   label("Endosomal degradation rate constant k_deg_e (1/h)")         # Table 1, human (allometric, b_deg=-0.24; rat 0.028, monkey 0.013)
    lkcle    <- log(1.32);    label("GalNAc cleavage rate from GalNAc-siRNA k_cle (1/h)")        # Table 1, shared across species: 1.32 (= k_off; Ayyar 2021)
    lkonapp  <- log(1.4e-5);  label("RISC association rate constant k_on_app (1/(nM h))")        # Table 1, shared rat/human: 1.4e-5 (36% RSE)
    lrisctot <- log(30);      label("Total liver RISC concentration RISC_tot (nM)")              # Table 1, fixed across species: 30 (paper Table 1)
    lkcomplx <- log(0.042);   label("RISC-loaded siRNA degradation rate k_complx (1/h)")         # Table 1, shared rat/human: 0.042 (31% RSE)
    lkdegc   <- log(0.01);    label("Cytosolic free-siRNA degradation rate k_deg_c (1/h)")       # Table 1, shared across species: 0.01

    # =============================================================
    # IIV - 20% CV assumed by paper for human Monte Carlo simulation
    # (Ayyar 2024 p.184: "Monte Carlo simulations, incorporating
    # variability and uncertainty in ka, Vc, and kint, were performed
    # to predict the plasma PK of givosiran and AS(N-1)3' givosiran in
    # 500 simulated subjects"; assumed 20% CV stated in same paragraph
    # and in the Discussion p.187).
    # omega^2 = log(CV^2 + 1) = log(1.04) = 0.0392.
    # =============================================================
    etalka   ~ 0.0392
    etalvc   ~ 0.0392
    etalkint ~ 0.0392
  })

  model({
    # 1. Individual parameters
    rtot    <- exp(lrtot)
    kdeg    <- exp(lkdeg)
    koff    <- exp(lkoff)
    kint    <- exp(lkint + etalkint)
    kd      <- exp(lkd)
    ka      <- exp(lka + etalka)
    vc      <- exp(lvc + etalvc)
    clm     <- exp(lclm)
    qkid    <- exp(lqkid)
    vk      <- exp(lvk)
    vkvas   <- exp(lvkvas)
    gfr     <- exp(lgfr)
    fuk     <- exp(lfuk)
    fugfr   <- exp(lfugfr)
    kassk   <- exp(lkassk)
    kdisk   <- exp(lkdisk)
    kassl   <- exp(lkassl)
    kdisl   <- exp(lkdisl)
    vhep    <- exp(lvhep)
    fesc    <- exp(lfesc)
    kdege   <- exp(lkdege)
    kcle    <- exp(lkcle)
    konapp  <- exp(lkonapp)
    risctot <- exp(lrisctot)
    kcomplx <- exp(lkcomplx)
    kdegc   <- exp(lkdegc)

    # 2. Secondary parameters (Table 1 footnotes)
    # k_syn = k_deg * R_tot (steady-state ASGPR turnover; nM/h)
    # k_on  = k_off / K_D   (1/(nM h))
    ksyn <- kdeg * rtot
    kon  <- koff / kd

    # 3. Concentrations derived from amount states
    # Plasma (parent + metabolite) live in central volume Vc:
    Cp     <- central / vc
    Cp_m   <- central_asn1 / vc
    # Free ASGPR and receptor complexes are also tracked in the central
    # volume Vc - the paper's simplifying assumption (text p.180:
    # "GalNAc-siRNA concentrations in plasma closely reflect those
    # available for binding ASGPR in liver interstitial fluid"):
    Rf     <- target / vc
    RC     <- complex / vc
    RCm    <- complex_asn1 / vc
    # Cytoplasmic and RISC concentrations live in V_hep:
    cytoC  <- cyto / vhep
    cytoCm <- cyto_asn1 / vhep
    RISC   <- risc / vhep
    # Kidney compartments:
    Ckvas   <- kid_vas / vkvas
    Ckvas_m <- kid_vas_asn1 / vkvas
    Ck      <- kid / vk
    Ck_m    <- kid_asn1 / vk

    # 4. Receptor-binding fluxes (nmol/h) - parent and metabolite
    flux_par_bind   <- kon * Cp   * Rf  * vc
    flux_par_unbind <- koff * complex
    flux_par_int    <- kint * complex
    flux_met_bind   <- kon * Cp_m * Rf  * vc
    flux_met_unbind <- koff * complex_asn1
    flux_met_int    <- kint * complex_asn1

    # RISC-loading fluxes (nmol/h)
    free_RISC <- risctot - RISC
    flux_par_risc <- konapp * cyto      * free_RISC
    flux_met_risc <- konapp * cyto_asn1 * free_RISC

    # 5. ODE system (paper Eqs 1-20 in amount form)
    # SC absorption depot
    d/dt(depot) <- -ka * depot

    # Eq 1: parent givosiran in central plasma (nmol)
    d/dt(central) <- ka * depot -
                     flux_par_bind + flux_par_unbind -
                     clm * Cp -
                     qkid * Cp + qkid * Ckvas

    # Eq 2: AS(N-1)3' metabolite in central plasma (nmol)
    d/dt(central_asn1) <- clm * Cp -
                          flux_met_bind + flux_met_unbind -
                          clm * Cp_m -
                          qkid * Cp_m + qkid * Ckvas_m

    # Eq 3: free ASGPR (target, nmol) - amount form of dRf/dt * Vc
    d/dt(target) <- ksyn * vc - kdeg * target -
                    flux_par_bind + flux_par_unbind -
                    flux_met_bind + flux_met_unbind

    # Eq 4: parent-ASGPR complex (nmol)
    d/dt(complex) <- flux_par_bind - flux_par_unbind - flux_par_int

    # Eq 5: metabolite-ASGPR complex (nmol)
    d/dt(complex_asn1) <- flux_met_bind - flux_met_unbind - flux_met_int

    # Eq 6: internalized parent-ASGPR complex (nmol)
    d/dt(liv) <- flux_par_int - kcle * liv

    # Eq 7: internalized metabolite-ASGPR complex (nmol)
    d/dt(liv_asn1) <- flux_met_int - kcle * liv_asn1

    # Eq 8: free endosomal parent siRNA (nmol)
    # Two kdeg_e terms: degradation and escape-to-cytosol; both retained
    # exactly as in paper Eq 8.
    d/dt(liv_endo) <- kcle * liv -
                      kdege * liv_endo - fesc * kdege * liv_endo -
                      kassl * liv_endo + kdisl * liv_deep

    # Eq 9: free endosomal metabolite siRNA (nmol)
    d/dt(liv_endo_asn1) <- kcle * liv_asn1 -
                           kdege * liv_endo_asn1 - fesc * kdege * liv_endo_asn1 -
                           kassl * liv_endo_asn1 + kdisl * liv_deep_asn1

    # Eq 10: liver endosomal parent bound pool (nmol)
    d/dt(liv_deep) <- kassl * liv_endo - kdisl * liv_deep

    # Eq 11: liver endosomal metabolite bound pool (nmol)
    d/dt(liv_deep_asn1) <- kassl * liv_endo_asn1 - kdisl * liv_deep_asn1

    # Eq 12: free cytoplasmic parent siRNA (nmol; paper concentration form * Vhep)
    d/dt(cyto) <- fesc * kdege * liv_endo -
                  kdegc * cyto -
                  flux_par_risc

    # Eq 13: free cytoplasmic metabolite siRNA (nmol)
    d/dt(cyto_asn1) <- fesc * kdege * liv_endo_asn1 -
                       kdegc * cyto_asn1 -
                       flux_met_risc

    # Eq 14: combined RISC-loaded siRNA (nmol)
    d/dt(risc) <- flux_par_risc + flux_met_risc - kcomplx * risc

    # Eq 15: parent in kidney vasculature (nmol)
    d/dt(kid_vas) <- qkid * Cp + qkid * fuk * Ck -
                     2 * qkid * Ckvas -
                     gfr * fugfr * Ckvas

    # Eq 16: metabolite in kidney vasculature (nmol)
    d/dt(kid_vas_asn1) <- qkid * Cp_m + qkid * fuk * Ck_m -
                          2 * qkid * Ckvas_m -
                          gfr * fugfr * Ckvas_m

    # Eq 17: parent in kidney tissue (nmol)
    d/dt(kid) <- qkid * Ckvas - qkid * fuk * Ck -
                 kassk * fuk * kid + kdisk * kid_deep

    # Eq 18: metabolite in kidney tissue (nmol)
    d/dt(kid_asn1) <- qkid * Ckvas_m - qkid * fuk * Ck_m -
                      kassk * fuk * kid_asn1 + kdisk * kid_deep_asn1

    # Eq 19: kidney bound pool parent (nmol)
    d/dt(kid_deep) <- kassk * fuk * kid - kdisk * kid_deep

    # Eq 20: kidney bound pool metabolite (nmol)
    d/dt(kid_deep_asn1) <- kassk * fuk * kid_asn1 - kdisk * kid_deep_asn1

    # 6. Bioavailability for SC depot (paper Eq 1: Input = ka * Dose * F)
    f(depot) <- exp(lfdepot)

    # 7. Initial conditions - free ASGPR at steady state.
    # State amount target(0) = R_tot * Vc (since Rf is treated as a
    # concentration in the central volume; see receptor section above).
    target(0) <- rtot * vc

    # 8. Observations - molar plasma concentrations (nmol/L = nM)
    # Convert to ng/mL externally with: ng_per_mL = nM * MW_giv / 1000
    # = nM * 16.3 (MW_givosiran = MW_AS(N-1)3' = 16,300 g/mol; paper p.179
    # Eq 23 / Eq 26).
    Cc      <- Cp
    Cc_asn1 <- Cp_m
  })
}
