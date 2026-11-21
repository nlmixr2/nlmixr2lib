Tiley_2026_metformin_placentaperfusion <- function() {
  description <- "Ex vivo (human term placenta, dual-side recirculating single-cotyledon perfusion). Three-compartment model of metformin transfer across an isolated perfused cotyledon: a maternal reservoir, the syncytiotrophoblast (placental barrier) and a fetal reservoir, linked by bidirectional passive diffusion clearances at the maternal-placental (CLPDM) and placental-fetal (CLPDF) barriers. The system is closed (no metabolism or loss). Fitted in Phoenix WinNonlin by the source authors to the closed-system cotyledon data of Nanovskaya 2006 (Tiley 2026 Equations 2-5, 'Approach 3'); the fitted clearances are printed per mL of whole-placenta volume and are converted back to per-cotyledon clearances with the paper's own Equation 1. The paper's Simcyp pregnancy PBPK model, into which these clearances were imputed, is NOT reproduced here."
  reference <- paste(
    "Tiley JB, Hartauer ME, Whigham TL, De Sousa Mendes M, Brouwer KLR, Hebert MF.",
    "Comparison of Metformin PBPK Models Incorporating Placental Transfer to Predict",
    "Fetal and Maternal Exposure. CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70136.",
    "doi:10.1002/psp4.70136. PMCID: PMC12823317.",
    "Model structure and fixed volumes and dose: Methods, Equations 2-5;",
    "unit conversion: Equation 1; estimates: Methods (Approach 3) and Table 1;",
    "fitted data: Table S1 (digitised by the source authors from Figure 1 of",
    "Nanovskaya et al., Am J Obstet Gynecol 2006;195:1081-1085).",
    sep = " "
  )
  vignette <- "Tiley_2026_metformin_placentaperfusion"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  # Metformin is added to the maternal reservoir in the ex vivo experiment.
  dosing <- c("maternal_reservoir")

  # The three states are the perfusion-circuit reservoirs and the
  # syncytiotrophoblast of the perfused cotyledon (Equations 2-5); none is a
  # whole-body compartment. `placental_barrier` follows the Bukkems 2022
  # doravirine perfusion model for the same tissue layer.
  paper_specific_compartments <- c(
    "maternal_reservoir",
    "placental_barrier",
    "fetal_reservoir"
  )

  compartmentData <- list(
    maternal_reservoir = list(analyte = "metformin", units = "ug", specimen = "not applicable", verified = TRUE),
    placental_barrier = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    fetal_reservoir = list(analyte = "metformin", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "ex vivo (human term placenta, isolated single-cotyledon dual perfusion)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = paste(
      "Closed (recirculating) dual perfusion of one cotyledon; maternal and fetal",
      "reservoir volumes 250 and 150 mL (Methods, Approach 3, citing the Nanovskaya",
      "2002 buprenorphine perfusion set-up)."
    ),
    dose_range = "metformin 1250 ug into the maternal reservoir (final concentration 5 ug/mL in 250 mL)",
    sampling = "maternal artery and fetal vein concentrations from 0 to 240 min (Table S1)",
    disease_state = "not applicable (ex vivo placental tissue from term deliveries)",
    notes = paste(
      "The fitted data are the mean maternal and fetal curves of Nanovskaya 2006",
      "Figure 1, digitised by the source authors (Table S1); the number of perfused",
      "placentas is in that earlier report, not in this paper. The paper's in vivo",
      "data (16 women with gestational diabetes on metformin 500 mg every 12 h)",
      "evaluate the Simcyp PBPK model, not this ex vivo model, and are out of scope",
      "for this file."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Experimental constants (Methods, Approach 3). Not estimated.
    # ------------------------------------------------------------------
    v_mres <- fixed(250)
    label("Maternal reservoir volume, Vm (mL)") # Methods, Approach 3: 'The maternal and fetal compartment volumes were set as the reservoir volumes of 250 mL and 150 mL, respectively'
    v_fres <- fixed(150)
    label("Fetal reservoir volume, Vf (mL)") # Methods, Approach 3, same sentence
    v_pb <- fixed(3.395)
    label("Syncytiotrophoblast volume of one cotyledon, Vs (mL)") # Methods, Approach 3: 'The syncytiotrophoblast volume was calculated as 3.395 mL (67.9 mL placental volume/20 cotyledon per placenta)'

    # ------------------------------------------------------------------
    # Fitted passive diffusion clearances of one cotyledon (mL/min). The
    # paper prints them only after its Equation 1 conversion to L/h per mL of
    # placental volume (0.0023 and 0.039; Methods, Approach 3 and Table 1).
    # Inverting Equation 1 with its stated 20 cotyledons per placenta and
    # 659 mL placental volume: CL_cot (mL/min) = value x 659 / 20 x 1000 / 60.
    # ------------------------------------------------------------------
    lcl_pdm <- log(1.263083)
    label("Maternal-to-placental passive diffusion clearance of one cotyledon, CLPDM,cot (mL/min)") # Table 1 'CLPDM (L/h/mL placental volume)' Approach 3 = 0.0023; 0.0023 x 659 / 20 x 1000 / 60 = 1.263083 mL/min
    lcl_pdf <- log(21.4175)
    label("Placental-to-fetal passive diffusion clearance of one cotyledon, CLPDF,cot (mL/min)") # Table 1 'CLPDF (L/h/mL placental volume)' Approach 3 = 0.039; 0.039 x 659 / 20 x 1000 / 60 = 21.4175 mL/min

    # No between-experiment variability or residual error is reported for this
    # Phoenix fit to the mean digitised curves, so the file carries typical
    # values only.
  })

  model({
    cl_pdm <- exp(lcl_pdm)
    cl_pdf <- exp(lcl_pdf)

    # Equation 2: concentrations (ug/mL) from amounts (ug) and volumes (mL).
    c_mres <- maternal_reservoir / v_mres
    c_pb <- placental_barrier / v_pb
    c_fres <- fetal_reservoir / v_fres

    # Equations 3-5: closed system, bidirectional passive diffusion only.
    d/dt(maternal_reservoir) <- cl_pdm * (c_pb - c_mres)
    d/dt(placental_barrier) <- cl_pdm * (c_mres - c_pb) - cl_pdf * (c_pb - c_fres)
    d/dt(fetal_reservoir) <- -cl_pdf * (c_fres - c_pb)

    # Observed: metformin in the maternal and fetal reservoirs (Table S1).
    Cmaternal <- c_mres
    Cfetal <- c_fres
  })
}
