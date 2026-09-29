Liu_2018_pneumocystis_mouse_qsp <- function() {
  description <- "QSP. Preclinical (mouse, C3H/HeN). Pharmacodynamic module of the Liu 2018 Pneumocystis quantitative systems pharmacology model: the untreated natural history of a Pneumocystis murina lung infection described as a two-stage life cycle of trophic forms and asci. Trophic forms proliferate first-order (KsTro), are lost by a second-order (logistic, crowding) death term (KdTro), and interconvert with asci through the forward transformation KTA and the reverse transformation KAT; asci do not proliferate and are lost by a first-order death term (KdAsci) that the paper sets to a near-zero value. The KsTro : KdTro ratio sets the carrying capacity, so the untreated organism burden plateaus at (KsTro - KTA + KAT*KTA/(KAT + KdAsci)) / KdTro = 1e7 organisms with the Table 4 basal values, reproducing the reported steady state of about 1e7 reached at roughly 35 days (Figure 3b, Figure 4c). Liu 2018 does not report the day-0 inoculum, so trophic0 / asci0 are exposed as parameters and defaulted to a single organism of each form, the value that reproduces the paper's own described time course (slow accumulation over the first two weeks, exponential growth from week three, plateau at the end of week five). This is the drug-free module that Figure 3 validates; the four integrated PK/PD models add a drug on top of exactly these equations -- see modellib('Liu_2018_anidulafungin_mouse_qsp'), modellib('Liu_2018_caspofungin_mouse_qsp'), modellib('Liu_2018_micafungin_mouse_qsp') and modellib('Liu_2018_sulfamethoxazole_mouse_qsp')."
  reference <- paste(
    "Liu GS, Ballweg R, Ashbaugh A, Zhang Y, Facciolo J, Cushion MT, Zhang T.",
    "(2018). A quantitative systems pharmacology (QSP) model for Pneumocystis",
    "treatment in mice. BMC Syst Biol 12(1):77.",
    "doi:10.1186/s12918-018-0603-9.",
    "Constraining organism-burden data (Liu 2018 reference 33) from",
    "Cushion MT, Linke MJ, Ashbaugh A, Sesterhenn T, Collins MS, Lynch K,",
    "Brubaker R, Walzer PD. (2010). Echinocandin treatment of pneumocystis",
    "pneumonia in rodent models depletes cysts leaving trophic burdens that",
    "cannot transmit the infection. PLoS One 5(1):e8524.",
    "doi:10.1371/journal.pone.0008524.",
    "The near-zero asci death rate is justified (Liu 2018 reference 44) by",
    "Icenhour CR, Kottom TJ, Limper AH. (2003). Evidence for a melanin cell",
    "wall component in Pneumocystis carinii. Infect Immun 71(9):5360-5363.",
    "doi:10.1128/IAI.71.9.5360-5363.2003.",
    "The two-state transformation structure follows the multistate",
    "tuberculosis model (Liu 2018 reference 32) of",
    "Clewe O, Aulin L, Hu Y, Coates AR, Simonsson US. (2016). A multistate",
    "tuberculosis pharmacometric model: a framework for studying",
    "anti-tubercular drug effects in vitro.",
    "J Antimicrob Chemother 71(4):964-974. doi:10.1093/jac/dkv416;",
    "see modellib('Clewe_2018_TB_MTP_GPDI_invitro')."
  )
  vignette <- "Liu_2018_pneumocystis_qsp"

  # The two Pneumocystis life-cycle stages are genuinely paper-mechanistic
  # states with no canonical analogue in inst/references/compartment-names.md:
  # `trophic` is the vegetative, actively proliferating trophic form and `asci`
  # the non-proliferating ascus (cyst) form that the echinocandins target. They
  # are NOT bacterial MTP subpopulations (`fbugs` / `sbugs` / `nbugs`), which
  # partition one organism by growth rate rather than by life-cycle stage.
  paper_specific_compartments <- c("trophic", "asci")

  units <- list(
    time = "day",
    dosing = "n/a (no drug is administered in the PD module)",
    concentration = "organisms per mouse lung (trophic, asci); log10 organisms for the log10_trophic / log10_asci observations"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE -- both entries were checked against
  # Liu 2018 Methods 'Construction of the PD module in mice' and Figure 2
  # (right panel).
  compartmentData <- list(
    trophic = list(
      analyte = "Pneumocystis murina trophic form",
      units = "organisms per mouse lung",
      specimen = "not applicable",
      verified = TRUE
    ),
    asci = list(
      analyte = "Pneumocystis murina asci (cyst form)",
      units = "organisms per mouse lung",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # No subject-level covariates. Liu 2018 induced variability by resampling the
  # five PD rate constants from a uniform distribution spanning 70-130% of the
  # Table 4 basal values, not by covariate effects; see the vignette for how to
  # reproduce that population.
  covariateData <- list()

  population <- list(
    species = "mouse (C3H/HeN)",
    n_subjects = "n = 2 or 3 immunosuppressed, Pneumocystis-infected mice per nuclei-count time point (Liu 2018 Figure 3b legend); the reported steady-state organism levels (log10 7.62 +/- 0.17 trophic forms, log10 7.79 +/- 0.13 asci) come from Cushion 2010 (Liu 2018 reference 33)",
    n_studies = 2L,
    age_range = "6 weeks old at the start of the study",
    sex_female_pct = 0,
    disease_state = "Pneumocystis murina lung infection in immunosuppressed mice (the model of Pneumocystis pneumonia, PCP)",
    dose_range = "not applicable (untreated natural history)",
    regions = "United States (University of Cincinnati; mice supplied by Charles River)",
    notes = "Liu 2018 Experimental methods 'Measuring pneumocystis numbers in mice': 6-week-old male C3H/HeN mice. Organism burden was quantified either by RT-qPCR of Pneumocystis mitochondrial large-subunit rRNA (total nuclei, cannot distinguish the two life-cycle stages) or by microscopic quantification with cresyl echt violet (asci) and a rapid Wright-Giemsa stain (all stages), which does distinguish them."
  )

  ini({
    # ==================================================================
    # PD module -- Liu 2018 Table 4 'Basal Parameter Values'. Every value
    # was chosen manually by the authors to recapture the observed data
    # ('all parameters sets and initial conditions were derived manually
    # using a trial and error method'), so none carries an uncertainty
    # estimate and all are encoded as fixed().
    # ==================================================================
    kstro <- fixed(1)
    label("Trophic form proliferation rate KsTro (1/day)")  # Table 4
    kdtro <- fixed(1e-7)
    label("Trophic form second-order death rate KdTro (1/day/organism)")  # Table 4
    kta <- fixed(0.1)
    label("Trophic form to asci transformation rate KTA (1/day)")  # Table 4
    kat <- fixed(0.1)
    label("Asci to trophic form transformation rate KAT (1/day)")  # Table 4
    kdasci <- fixed(2e-12)
    label("Asci death rate KdAsci (1/day)")  # Table 4

    # ==================================================================
    # Initial organism burden. NOT REPORTED BY LIU 2018 -- the paper states
    # only that 'the initial conditions for each drug were estimated from
    # the literature data when available' for the PK module and says
    # nothing about the day-0 Pneumocystis inoculum. These two values were
    # back-solved from the paper's own description of the untreated time
    # course in Results ('The initial growth of the organism was very slow
    # within the first two weeks, however, starting from the third week, an
    # exponential growth of Pneumocystis was observed which peaked at the
    # end of the fifth week') plus 'At about 35 days, the levels of both
    # trophic forms and asci reached a steady state of about 10^7'. A unit
    # inoculum of each form reproduces that timeline: the linearised net
    # growth rate of the drug-free system is 0.910/day, so a single
    # organism reaches 98% of the 1e7 plateau at day 35 and crosses half
    # of it at day 18. Override with rxSolve(inits = ...) for any other
    # infection scenario. See the vignette 'Assumptions and deviations'.
    # ==================================================================
    trophic0 <- fixed(1)
    label("Initial trophic form burden (organisms per lung)")  # back-solved; not reported
    asci0 <- fixed(1)
    label("Initial asci burden (organisms per lung)")  # back-solved; not reported
  })

  model({
    # ==================================================================
    # Two-stage Pneumocystis life cycle -- Liu 2018 Table 4 / Equations
    # (d) and (e).
    #
    #   dTro/dt  = KsTro*Tro - KdTro*Tro*Tro - KTA*Tro + KAT*Asci
    #   dAsci/dt = KTA*Tro - KAT*Asci - KdAsci*Asci
    #
    # The trophic loss term is deliberately SECOND order: 'Following
    # logistic growth models, the decay of the trophic form is a second
    # order reaction since the trophic forms actively proliferate and
    # compete for space and nutrients. On the contrary, the asci do not
    # actively proliferate but rather result from the transformation of
    # trophic forms. Hence, the decay of the asci is set to be a first
    # order reaction.' (Methods, 'Construction of the PD module in mice')
    # ==================================================================
    d/dt(trophic) <- kstro * trophic - kdtro * trophic * trophic -
      kta * trophic + kat * asci
    d/dt(asci) <- kta * trophic - kat * asci - kdasci * asci

    trophic(0) <- trophic0
    asci(0) <- asci0

    # Liu 2018 plots both stages on a log10 scale (Figures 3a, 4c, 4d) and
    # reports the observed burdens as log10 values (7.62 trophic forms,
    # 7.79 asci), so both log10 transforms are exposed as observations.
    log10_trophic <- log10(trophic)
    log10_asci <- log10(asci)
  })
}
