Wang_2024_ionizableLipid_hela_qsp <- function() {
  description <- paste(
    "In vitro (HeLa cell line).",
    "QSP (cellular PBPK, MATLAB SimBiology 2022a) model of siRNA-lipid",
    "nanoparticle (LNP) trafficking inside HeLa cells, tracking the lipid and",
    "the siRNA cargo separately. Internalised LNP matures from the early",
    "endosome or macropinosome to the late endosome and then to the lysosome.",
    "From the late endosome a fraction of vesicles undergoes an endosomal",
    "escape (release) event at rate krel; of the siRNA in those vesicles a",
    "fraction frel reaches the cytoplasm and the remainder, together with all",
    "the lipid, is engulfed by an autophagosome that matures into an",
    "autolysosome. Late endosomes also egress material back out of the cell at",
    "rate keg. This is the 'reduced' structure of Figure 3, in which LNP",
    "disassembly is not simulated explicitly. Three ionizable lipids are",
    "covered: C12-200 (the reference) plus MC3 and L319, selected with",
    "FORM_LNP_MC3 / FORM_LNP_L319; they differ only in krel.",
    sep = " "
  )
  reference <- paste(
    "Wang W, Deng S, Lin J, Ouyang D (2024).",
    "Modeling on in vivo disposition and cellular transportation of RNA lipid",
    "nanoparticles via quantum mechanics/physiologically-based pharmacokinetic",
    "approaches. Acta Pharm Sin B 14(10):4591-4607.",
    "doi:10.1016/j.apsb.2024.06.011.",
    "Model structure from Figure 3 (right panel) and Eqs. (10)-(17); the full",
    "SimBiology ODE export is in the Supporting Information (mmc1.pdf) section",
    "4 ('The reduced cell model'); fitted parameters from Figures 8C and 9D.",
    "The underlying HeLa imaging data are from the three studies Wang 2024",
    "cites as references 8, 9 and 17.",
    sep = " "
  )
  vignette <- "Wang_2024_rnaLipidNanoparticle"
  units <- list(
    time          = "h",
    dosing        = "ng",
    concentration = "fraction of internalised siRNA (dimensionless)"
  )

  # Every state is a mass (ng) of lipid or siRNA in one intracellular vesicle
  # population, or in the culture medium. These are paper-mechanistic
  # endosomal-trafficking states with no analogue in the canonical
  # compartment register; the `lnp_` prefix marks material still associated
  # with an intact LNP and the bare prefix marks material that has been
  # released or egressed.
  paper_specific_compartments <- c(
    "lnp_lipid_ee", "lnp_rna_ee",
    "lnp_lipid_le", "lnp_rna_le",
    "lnp_lipid_ly", "lnp_rna_ly",
    "lipid_ap", "rna_ap",
    "lipid_al", "rna_al",
    "rna_cyto",
    "lipid_med", "rna_med"
  )

  compartmentData <- list(
    lnp_lipid_ee = list(analyte = "ionizable lipid (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lnp_rna_ee   = list(analyte = "siRNA (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lnp_lipid_le = list(analyte = "ionizable lipid (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lnp_rna_le   = list(analyte = "siRNA (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lnp_lipid_ly = list(analyte = "ionizable lipid (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lnp_rna_ly   = list(analyte = "siRNA (in LNP)", units = "ng", specimen = "endosome", verified = TRUE),
    lipid_ap     = list(analyte = "ionizable lipid (released vesicle)", units = "ng", specimen = "endosome", verified = TRUE),
    rna_ap       = list(analyte = "siRNA (released vesicle)", units = "ng", specimen = "endosome", verified = TRUE),
    lipid_al     = list(analyte = "ionizable lipid (released vesicle)", units = "ng", specimen = "endosome", verified = TRUE),
    rna_al       = list(analyte = "siRNA (released vesicle)", units = "ng", specimen = "endosome", verified = TRUE),
    rna_cyto     = list(analyte = "siRNA (free, cytoplasmic)", units = "ng", specimen = "not applicable", verified = TRUE),
    lipid_med    = list(analyte = "ionizable lipid (egressed)", units = "ng", specimen = "not applicable", verified = TRUE),
    rna_med      = list(analyte = "siRNA (egressed)", units = "ng", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    FORM_LNP_MC3 = list(
      description        = "1 = LNP formulated with the ionizable lipid DLin-MC3-DMA (MC3); 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C12-200, the formulation used to fit the shared trafficking rates)",
      notes              = paste(
        "Selects the MC3 column of Wang 2024 Figure 9D (krel = 0.0058 1/h).",
        "MC3's own frel was not measured; the paper assumes it equals the",
        "L319 value of 0.5 because the two lipids are structurally similar",
        "(Figure 9D footnote #, and section 3.2.2). Mutually exclusive with",
        "FORM_LNP_L319."
      ),
      source_name        = "MC3"
    ),
    FORM_LNP_L319 = list(
      description        = "1 = LNP formulated with the ionizable lipid L319; 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C12-200, the formulation used to fit the shared trafficking rates)",
      notes              = paste(
        "Selects the L319 column of Wang 2024 Figure 9D (krel = 0.0157 1/h).",
        "L319 is the only formulation for which frel was measured directly",
        "(0.5, from the fraction of galectin-8-positive vesicles that yielded",
        "a detectable cytoplasmic signal). Mutually exclusive with",
        "FORM_LNP_MC3."
      ),
      source_name        = "L319"
    )
  )

  population <- list(
    species       = "in vitro (HeLa cell line)",
    n_subjects    = NA_integer_,
    n_studies     = 3L,
    disease_state = "not applicable (cell culture)",
    dose_range    = paste(
      "LNP incubated with HeLa cells; 3 h exposure followed by a wash and 25 h",
      "of monitoring for C12-200, 6 h continuous exposure for MC3, and an",
      "uptake time course for L319"
    ),
    notes = paste(
      "Confocal and electron microscopy time courses of LNP uptake, vesicle",
      "co-localisation and siRNA release digitised from the three HeLa studies",
      "Wang 2024 cites as references 8, 9 and 17. The trafficking rates",
      "kEE_LE, kLE_LY (= kAP_AL) and keg were fitted simultaneously to the",
      "C12-200 and MC3 datasets and are treated as cell-physiological",
      "constants shared by all formulations; only krel (and frel) are",
      "formulation-specific. Fitted to pooled mean profiles, so the model",
      "carries no between-subject variability."
    )
  )

  ini({
    # ================================================================
    # Cell-physiological trafficking rates, shared across formulations
    # (Figure 8C, "Reduced model" column). kAP_AL is assumed equal to
    # kLE_LY (paper section 2.2.2: "In addition, kLE_LY is assumed to be
    # equal to kAP_AL"), so a single parameter carries both.
    # ================================================================
    lkee_le <- log(3.3820) ; label("Early-endosome to late-endosome transport rate (1/h)")        # Fig 8C, "kEE_LE", reduced model
    lkle_ly <- log(0.2094) ; label("Late-endosome to lysosome and autophagosome to autolysosome transport rate (1/h)") # Fig 8C, "kLE_LY and kAP_AL", reduced model
    lkeg    <- log(0.5222) ; label("Egress rate out of the cell from the late endosome (1/h)")    # Fig 8C, "keg", reduced model

    # ================================================================
    # Formulation-specific endosomal-escape (release) rate, Figure 9D.
    # ================================================================
    lkrel_c12200 <- log(0.0080) ; label("siRNA release rate, C12-200 (1/h)") # Fig 9D, "krel", C12-200 column (also Fig 8C)
    lkrel_mc3    <- log(0.0058) ; label("siRNA release rate, MC3 (1/h)")     # Fig 9D, "krel", MC3 column
    lkrel_l319   <- log(0.0157) ; label("siRNA release rate, L319 (1/h)")    # Fig 9D, "krel", L319 column

    # frel is a measured fraction, not an estimated parameter: section 3.2.2
    # reports that for L319 "only half of siRNA could eventually escape into
    # the cytoplasm (frel = 0.5) even if the release event happens to an
    # endosome, which was evidenced by the detected fluorescence signal", and
    # Figure 9D footnote # records that the same value was assumed for MC3.
    # It was not obtainable for C12-200 (Figure 9D "-, not applicable").
    frel <- fixed(0.5) ; label("Fraction of released siRNA reaching the cytoplasm, measured for L319 and assumed equal for MC3 (unitless)") # Fig 9D, "frel"
  })

  model({
    # ------------------------------------------------------------------
    # 1. Formulation selection. Both indicators zero selects C12-200.
    # ------------------------------------------------------------------
    krel <- exp(lkrel_c12200) * (1 - FORM_LNP_MC3) * (1 - FORM_LNP_L319) +
      exp(lkrel_mc3) * FORM_LNP_MC3 +
      exp(lkrel_l319) * FORM_LNP_L319

    kee_le <- exp(lkee_le)
    kle_ly <- exp(lkle_ly)
    kap_al <- kle_ly          # assumed equal, paper section 2.2.2
    keg    <- exp(lkeg)

    # ------------------------------------------------------------------
    # 2. ODE system, Supporting Information section 4.2 (reduced model),
    #    Eqs. (10)-(17). States are masses (ng).
    #
    #    LNP uptake into the early endosome is encoded as a dose rather than
    #    as an explicit input term. The authors modelled uptake as an
    #    empirical function of time fitted separately to each study's uptake
    #    curve and did not report its parameters; every quantity the paper
    #    reports from this model is a ratio of first-order fluxes downstream
    #    of the early endosome, so it is independent of the shape of the
    #    input (see the vignette, where the same asymptotes are reproduced
    #    from a bolus and from a slow infusion).
    #
    #    Metabolism of lipid and siRNA is deliberately absent: Table S1
    #    assumption 5 states that metabolism is not modelled in the cellular
    #    PBPK model, and the SimBiology export's kel term on the lysosomal
    #    siRNA pool is therefore zero.
    # ------------------------------------------------------------------
    # Early endosome / macropinosome
    d/dt(lnp_lipid_ee) <- -kee_le * lnp_lipid_ee
    d/dt(lnp_rna_ee)   <- -kee_le * lnp_rna_ee

    # Late endosome: matures to lysosome, egresses, or undergoes release
    d/dt(lnp_lipid_le) <- kee_le * lnp_lipid_ee -
      kle_ly * lnp_lipid_le - keg * lnp_lipid_le - krel * lnp_lipid_le
    # The two release branches for siRNA, krel * frel to the cytoplasm and
    # krel * (1 - frel) to the autophagosome, together remove krel from the
    # late-endosomal pool.
    d/dt(lnp_rna_le)   <- kee_le * lnp_rna_ee -
      kle_ly * lnp_rna_le - keg * lnp_rna_le - krel * lnp_rna_le

    # Lysosome (no metabolism, Table S1 assumption 5)
    d/dt(lnp_lipid_ly) <- kle_ly * lnp_lipid_le
    d/dt(lnp_rna_ly)   <- kle_ly * lnp_rna_le

    # Endosomal escape: siRNA to the cytoplasm, Eq. (14)
    d/dt(rna_cyto) <- frel * krel * lnp_rna_le

    # Damaged vesicles engulfed by autophagosomes, Eqs. (15)-(16), then
    # maturing to autolysosomes, Eq. (17). All the lipid of a vesicle that
    # released follows this route; only the (1 - frel) part of its siRNA does.
    d/dt(lipid_ap) <- krel * lnp_lipid_le - kap_al * lipid_ap
    d/dt(rna_ap)   <- (1 - frel) * krel * lnp_rna_le - kap_al * rna_ap
    d/dt(lipid_al) <- kap_al * lipid_ap
    d/dt(rna_al)   <- kap_al * rna_ap

    # Egress back into the culture medium, Eq. (13)
    d/dt(lipid_med) <- keg * lnp_lipid_le
    d/dt(rna_med)   <- keg * lnp_rna_le

    # ------------------------------------------------------------------
    # 3. Observations, from the SimBiology "repeated assignment" block
    #    (Supporting Information section 4.1). All are dimensionless
    #    fractions, so they do not depend on the dose amount.
    # ------------------------------------------------------------------
    # Total siRNA still inside the cell (`cellR` in the export).
    rnaCell <- lnp_rna_ee + lnp_rna_le + rna_ap + lnp_rna_ly + rna_al + rna_cyto

    # `occurfrac`: fraction of intracellular siRNA in vesicles that have
    # undergone a release event. Plotted for L319 in Figure 9B against the
    # observed 0.07 galectin-8 co-localisation fraction.
    fracReleaseEvent <- (rna_ap + rna_al + rna_cyto) / rnaCell

    # `cRfrac`: fraction of intracellular siRNA that actually reached the
    # cytoplasm. Plotted for MC3 in Figure 9C against the observed 0.0134.
    fracCytoplasm <- rna_cyto / rnaCell

    # `occurfracR`: probability that a molecule of internalised siRNA
    # undergoes a release event, with egressed siRNA kept in the denominator.
    # Plotted in Figures 9B and 9C against the observed 0.0228 for L319.
    probReleaseEvent <- (rna_cyto + rna_al) /
      (rna_cyto + rna_al + lnp_rna_ly + rna_med)

    # Vesicle-population fractions of siRNA, replicating Figure 8B.
    rnaEE <- lnp_rna_ee
    rnaLE <- lnp_rna_le + rna_ap
    rnaLY <- lnp_rna_ly + rna_al
    fracEE <- rnaEE / (rnaEE + rnaLE + rnaLY)
    fracLE <- rnaLE / (rnaEE + rnaLE + rnaLY)
    fracLY <- rnaLY / (rnaEE + rnaLE + rnaLY)

    # Lipid counterpart of rnaCell, for lipid mass-balance checks.
    lipidCell <- lnp_lipid_ee + lnp_lipid_le + lipid_ap + lnp_lipid_ly + lipid_al
  })
}
