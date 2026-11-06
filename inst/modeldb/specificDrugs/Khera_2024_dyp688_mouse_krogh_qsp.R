Khera_2024_dyp688_mouse_krogh_qsp <- function() {
  description <- "QSP. Preclinical (mouse, 92.1 uveal melanoma xenograft). Mechanistic Krogh cylinder model of tumour delivery and intracellular target engagement of DYP688, a PMEL17-targeted antibody-drug conjugate carrying two FR900359 (Galpha-q/11 inhibitor) payloads that inactivate while still attached to the antibody (Khera 2024, Supplementary Data 'Krogh cylinder model equations' and Supplementary Table 1). The typical-value DYP688 mouse PK model (three ADC species by number of active payloads: dar2 / dar1 / dar0) supplies the plasma concentration of each species. Each species extravasates across the capillary wall (Robin boundary condition), diffuses radially through the tumour interstitium, binds competitively to cell-surface PMEL17, and is internalised; in the lysosome linker cleavage and antibody degradation release active FR900359 in proportion to each species' active-payload count. Released payload escapes to the cytoplasm, binds Galpha-q/11 reversibly (with Galpha-q/11 turnover), and can leave and re-enter cells and wash out through the capillary (plasma payload = 0). Inactive payload is inert and does not inactivate in the tumour. The authors solved the radial partial differential equations in MATLAB; here they are discretised by the method of lines into 10 equal-width annular shells from the capillary wall (8 um) to the Krogh radius (100 um), one ODE state per species per shell (states <species>_slab1 ... _slab10, shell 1 at the capillary). gq_free_frac is the radial average of free Galpha-q/11 as a fraction of baseline (Khera 2024 Figure 4B). No variability: the model is fully predictive and was not fitted to tumour data. Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr."
  reference <- paste(
    "Khera E, Dharmarajan L, Hainzl D, Engelhardt V, Vostiarova H, Davis J,",
    "Ebel N, Wuersch K, Romanet V, Sharaby S, Kearns JD. (2024). QSP modeling",
    "of a transiently inactivating antibody-drug conjugate highlights benefit",
    "of short antibody half life. J Pharmacokinet Pharmacodyn 52(1):7.",
    "doi:10.1007/s10928-024-09956-1."
  )
  vignette <- "Khera_2024_dyp688"
  units <- list(
    time = "h",
    dosing = "nmol (antibody; one record per ADC species, see description)",
    concentration = "nM for plasma ADC and for every tumour species (concentration per tumour volume as written in the Krogh equations)"
  )

  compartmentData <- list(
    dar2_central = list(
      analyte = "DYP688 with two active FR900359 payloads",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar1_central = list(
      analyte = "DYP688 with one active FR900359 payload",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar0_central = list(
      analyte = "DYP688 with no active FR900359 payload",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar2_peripheral1 = list(
      analyte = "DYP688 with two active FR900359 payloads",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar1_peripheral1 = list(
      analyte = "DYP688 with one active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar0_peripheral1 = list(
      analyte = "DYP688 with no active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    adc_free_twoactive_slab1 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab1 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab1 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab1 = list(
      analyte = "free cell-surface PMEL17, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab1 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab1 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab1 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab1 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab1 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab1 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab1 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab1 = list(
      analyte = "extracellular active FR900359, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab1 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab1 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab1 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 1",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab2 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab2 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab2 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab2 = list(
      analyte = "free cell-surface PMEL17, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab2 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab2 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab2 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab2 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab2 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab2 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab2 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab2 = list(
      analyte = "extracellular active FR900359, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab2 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab2 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab2 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 2",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab3 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab3 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab3 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab3 = list(
      analyte = "free cell-surface PMEL17, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab3 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab3 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab3 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab3 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab3 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab3 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab3 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab3 = list(
      analyte = "extracellular active FR900359, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab3 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab3 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab3 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 3",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab4 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab4 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab4 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab4 = list(
      analyte = "free cell-surface PMEL17, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab4 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab4 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab4 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab4 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab4 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab4 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab4 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab4 = list(
      analyte = "extracellular active FR900359, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab4 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab4 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab4 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 4",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab5 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab5 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab5 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab5 = list(
      analyte = "free cell-surface PMEL17, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab5 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab5 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab5 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab5 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab5 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab5 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab5 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab5 = list(
      analyte = "extracellular active FR900359, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab5 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab5 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab5 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 5",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab6 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab6 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab6 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab6 = list(
      analyte = "free cell-surface PMEL17, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab6 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab6 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab6 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab6 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab6 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab6 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab6 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab6 = list(
      analyte = "extracellular active FR900359, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab6 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab6 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab6 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 6",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab7 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab7 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab7 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab7 = list(
      analyte = "free cell-surface PMEL17, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab7 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab7 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab7 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab7 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab7 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab7 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab7 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab7 = list(
      analyte = "extracellular active FR900359, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab7 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab7 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab7 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 7",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab8 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab8 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab8 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab8 = list(
      analyte = "free cell-surface PMEL17, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab8 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab8 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab8 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab8 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab8 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab8 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab8 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab8 = list(
      analyte = "extracellular active FR900359, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab8 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab8 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab8 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 8",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab9 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab9 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab9 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab9 = list(
      analyte = "free cell-surface PMEL17, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab9 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab9 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab9 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab9 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab9 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab9 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab9 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab9 = list(
      analyte = "extracellular active FR900359, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab9 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab9 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab9 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 9",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_twoactive_slab10 = list(
      analyte = "free DYP688 with two active FR900359 payloads (tissue concentration), radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_oneactive_slab10 = list(
      analyte = "free DYP688 with one active FR900359 payload (tissue concentration), radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_free_noactive_slab10 = list(
      analyte = "free DYP688 with no active FR900359 payload (tissue concentration), radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_free_slab10 = list(
      analyte = "free cell-surface PMEL17, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_twoactive_slab10 = list(
      analyte = "PMEL17-bound DYP688 with two active payloads, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_oneactive_slab10 = list(
      analyte = "PMEL17-bound DYP688 with one active payload, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_noactive_slab10 = list(
      analyte = "PMEL17-bound DYP688 with no active payload, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_twoactive_slab10 = list(
      analyte = "internalised (lysosomal) DYP688 with two active payloads, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_oneactive_slab10 = list(
      analyte = "internalised (lysosomal) DYP688 with one active payload, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_internal_noactive_slab10 = list(
      analyte = "internalised (lysosomal) DYP688 with no active payload, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_lyso_slab10 = list(
      analyte = "released active FR900359 in the lysosome, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_ext_slab10 = list(
      analyte = "extracellular active FR900359, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_int_slab10 = list(
      analyte = "free intracellular (cytoplasmic) active FR900359, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    payload_bound_slab10 = list(
      analyte = "intracellular FR900359 bound to Galpha-q/11, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gq_free_slab10 = list(
      analyte = "free intracellular Galpha-q/11, radial shell 10",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female nude, 92.1 uveal melanoma xenograft); simulation model",
    n_subjects = NA_integer_,
    disease_state = "Subcutaneous 92.1 uveal melanoma cell-derived xenograft with moderate PMEL17 expression (approximately 100,000 receptors per cell)",
    dose_range = "Simulated at 1.5, 6 and 12.5 mg/kg single IV doses (Figure 4B) and the dose-fractionation regimens of the TGI study (Figure S9)",
    notes = "The Krogh cylinder model is not fitted to data (Khera 2024 Results: 'fully predictive and not fit to any data'). Plasma input is the DYP688 mouse PK model (Table 1). Predicted Galpha-q/11 inhibition was compared qualitatively with tumour pERK immunohistochemistry H-scores (Figure 4B)."
  )

  ini({
    # Plasma PK input: DYP688 mouse PK typical values (Table 1, full
    # precision from the Supplementary Material Monolix output).
    lkel <- log(0.0161)
    label("Antibody / ADC elimination rate constant kel (1/h)") # Supplement Monolix output kel_pop = 0.0161; Table 1 0.016
    lk12 <- log(0.0355)
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Supplement Monolix output k12_pop = 0.0355; Table 1 0.035
    lk21 <- log(0.0356)
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Supplement Monolix output k21_pop = 0.0356; Table 1 0.036
    lkdec <- log(0.0125)
    label("On-antibody inactivation rate constant of each active payload in plasma ki.a. (1/h)") # Supplement Monolix output kia1_pop = 0.0125; Table 1 0.013
    lvc <- log(1.53)
    label("Plasma volume of distribution Vc (mL)") # Supplement Monolix output Vc_pop = 1.53; Table 1 1.53
    fr <- 0.0274
    label("Fraction of inactive ADC in the dose fr (fraction)") # Supplement Monolix output fr_pop = 0.0274; Table 1 0.027
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # Supplement Monolix output xfr_pop = 0.99 (no SE)

    # Krogh cylinder geometry and antibody transport (Supplementary
    # Table 1, printed units; converted to um and h in model()).
    rkrogh <- 100
    label("Krogh cylinder radius, half the intercapillary distance R_Krogh (um)") # Supplementary Table 1, R_Krogh = 100 um
    rcap <- 8
    label("Capillary radius R_capillary (um)") # Supplementary Table 1, R_capillary = 8 um
    p_ab <- 3e-9
    label("Vascular permeability of antibody / ADC P_ADC (m/s)") # Supplementary Table 1, P_ADC = 3 x 10^-9 m/s
    d_ab <- 10
    label("Antibody / ADC effective diffusivity in tumour, same for all ADC species D_ADC (um^2/s)") # Supplementary Table 1, D_ADC = 10 um^2/s
    eps <- 0.24
    label("Tumour void fraction available to antibody epsilon (fraction)") # Supplementary Table 1, epsilon = 0.24

    # PMEL17 binding and trafficking (Supplementary Table 1).
    ag0 <- 83.3
    label("Initial cell-surface PMEL17 concentration [Ag]0, 100,000 receptors per cell (nM)") # Supplementary Table 1, [Ag]0 = 83.3 nM (estimated from IHC)
    kon_ab <- 4.07e5
    label("Antibody / ADC association rate constant k_on,Ab (1/(M*s))") # Supplementary Table 1, k_on,Ab = 4.07 x 10^5 M^-1 s^-1 (measured)
    koff_ab <- 9.7e-3
    label("Antibody / ADC dissociation rate constant k_off,Ab (1/s)") # Supplementary Table 1, k_off,Ab = 9.7 x 10^-3 s^-1 (k_on,Ab * K_D, K_D = 23.9 nM)
    ke_ag <- 7.7e-5
    label("Internalisation rate constant of free PMEL17 k_e (1/s)") # Supplementary Table 1, k_e = 7.7 x 10^-5 s^-1
    kint <- 7.7e-5
    label("Internalisation rate constant of PMEL17-bound ADC k_int (1/s)") # Supplementary Table 1, k_int = 7.7 x 10^-5 s^-1 (assumed = k_e)
    krel <- 1.925e-5
    label("Lysosomal valine-citrulline linker hydrolysis rate constant per payload k_rel (1/s)") # Supplementary Table 1, k_rel = 1.925 x 10^-5 s^-1
    kdeg_lyso <- 8.0e-6
    label("Lysosomal ADC (antibody) degradation rate constant k_deg (1/s)") # Supplementary Table 1, k_deg = 8.0 x 10^-6 s^-1

    # Active FR900359 payload transport and Galpha-q/11 binding
    # (Supplementary Table 1).
    kin_pl <- 1.44e-4
    label("Payload cell influx (and lysosomal escape) rate constant k_in,P (1/s)") # Supplementary Table 1, k_in,P = 1.44 x 10^-4 s^-1
    kout_pl <- 7.06e-5
    label("Payload cell efflux rate constant k_out,P (1/s)") # Supplementary Table 1, k_out,P = 7.06 x 10^-5 s^-1 (k_in,P / (1 + R))
    d_pl <- 7.63
    label("Payload effective diffusion coefficient in tumour D_P,eff (um^2/s)") # Supplementary Table 1, D_P,eff = 7.63 um^2/s (D_P,equil / (1 + R))
    p_pl <- 1e-6
    label("Vascular permeability of payload P_P (m/s)") # Supplementary Table 1, P_P = 1 x 10^-6 m/s (estimated)
    part_pl <- 1.03
    label("Payload partition coefficient R (unitless)") # Supplementary Table 1, R = 1.03 (calculated; printed unit 'm/s')
    kon_pl <- 7.25e4
    label("Payload-Galpha-q association rate constant k_on,P (1/(M*s))") # Supplementary Table 1, k_on,P = 7.25 x 10^4 M^-1 s^-1 (k_off,P / K_D,P, K_D,P = 0.59 nM)
    koff_pl <- 4.3e-5
    label("Payload-Galpha-q dissociation rate constant k_off,P (1/s)") # Supplementary Table 1, k_off,P = 4.3 x 10^-5 (printed unit 'M^-1 s^-1')
    eps_pl <- 0.44
    label("Cell void fraction epsilon_P (fraction)") # Supplementary Table 1, epsilon_P = 0.44 (estimated)
    gq0 <- 400
    label("Baseline intracellular Galpha-q/11 concentration Int_Gq (nM)") # Supplementary Table 1, Int_Gq = 400 nM (assumed)
    kdeg_gq <- 8.0e-6
    label("Galpha-q/11 turnover rate constant k_deg,Gq (1/s)") # Supplementary Table 1, k_deg,Gq = 8.0 x 10^-6 s^-1 (assumed)
  })

  model({
    # ---- Plasma PK (Supplement 'Mouse PK model structure file') ---------
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kdec <- exp(lkdec)
    vc <- exp(lvc)

    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - 2 * kdec * dar2_central - kel * dar2_central
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 + 2 * kdec * dar2_central - kdec * dar1_central - kel * dar1_central
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 + kdec * dar1_central - kel * dar0_central
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1 - 2 * kdec * dar2_peripheral1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1 + 2 * kdec * dar2_peripheral1 - kdec * dar1_peripheral1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1 + kdec * dar1_peripheral1

    f(dar2_central) <- 1 - fr
    f(dar1_central) <- (1 - xfr) * fr
    f(dar0_central) <- xfr * fr

    # Plasma concentration of each ADC species, nM (nmol states, vc in mL):
    # C_plasma,ADC in the capillary boundary condition of each species.
    cp_dar2 <- 1000 * dar2_central / vc
    cp_dar1 <- 1000 * dar1_central / vc
    cp_dar0 <- 1000 * dar0_central / vc
    Cc <- cp_dar2 + cp_dar1 + cp_dar0
    Cc_active_payload <- 2 * cp_dar2 + cp_dar1

    # ---- Unit conversions to um, nM and h -------------------------------
    p_ab_h <- p_ab * 1e6 * 3600
    d_ab_h <- d_ab * 3600
    p_pl_h <- p_pl * 1e6 * 3600
    d_pl_h <- d_pl * 3600
    kon_ab_h <- kon_ab * 1e-9 * 3600
    koff_ab_h <- koff_ab * 3600
    ke_ag_h <- ke_ag * 3600
    kint_h <- kint * 3600
    krel_h <- krel * 3600
    kdeg_lyso_h <- kdeg_lyso * 3600
    kin_pl_h <- kin_pl * 3600
    kout_pl_h <- kout_pl * 3600
    koff_pl_h <- koff_pl * 3600
    kdeg_gq_h <- kdeg_gq * 3600
    # Antigen recycling keeps surface PMEL17 at [Ag]0 (R_s = k_e * [Ag]0)
    # and Galpha-q/11 synthesis keeps it at Int_Gq (k_syn,Gq =
    # k_deg,Gq * Int_Gq) (Supplementary Table 1).
    rs <- ke_ag_h * ag0
    ksyn_gq <- kdeg_gq_h * gq0
    # Payload-Galpha-q binding enters as k_on,P / ((1 - epsilon_P) * R).
    kon_gq <- kon_pl * 1e-9 * 3600 / ((1 - eps_pl) * part_pl)

    # ---- Method-of-lines grid --------------------------------------------
    # Shell i spans r = rcap + (i - 1) * dr to rcap + i * dr; r_face<i> is
    # its outer face and area_slab<i> its area per radian, so the cylindrical
    # diffusion operator (1/r) d/dr (r dC/dr) becomes the net radial flux
    # through the two faces divided by the shell area. At the capillary face
    # the flux is the Robin condition P * (C_plasma - C / epsilon) (ADC) or
    # P_P * (C_plasma,P - C_ext,P / (epsilon * (1 + R))) with C_plasma,P = 0
    # (payload); at the Krogh radius the flux is zero (Neumann).
    #
    # Per shell, following the Supplementary Data equations: free ADC
    # species C_2P, C_1P, C_0P; free target T_free; bound B_2P, B_1P, B_0P;
    # internalised C_int,2P, C_int,1P, C_int,0P; lysosomal C_lyso,P;
    # extracellular C_ext,P; intracellular C_int,P; Galpha-q-bound
    # C_bound,P; free Int_Gq. Lysosomal payload is released at 2 k_rel and
    # k_rel (linker cleavage) and 2 k_deg and k_deg (antibody degradation)
    # from the two- and one-active-payload species and escapes to the
    # cytoplasm at k_in,P.
    dr <- (rkrogh - rcap) / 10
    r_face1 <- rcap + 1 * dr
    r_face2 <- rcap + 2 * dr
    r_face3 <- rcap + 3 * dr
    r_face4 <- rcap + 4 * dr
    r_face5 <- rcap + 5 * dr
    r_face6 <- rcap + 6 * dr
    r_face7 <- rcap + 7 * dr
    r_face8 <- rcap + 8 * dr
    r_face9 <- rcap + 9 * dr
    area_slab1 <- ((rcap + 1 * dr)^2 - (rcap + 0 * dr)^2) / 2
    area_slab2 <- ((rcap + 2 * dr)^2 - (rcap + 1 * dr)^2) / 2
    area_slab3 <- ((rcap + 3 * dr)^2 - (rcap + 2 * dr)^2) / 2
    area_slab4 <- ((rcap + 4 * dr)^2 - (rcap + 3 * dr)^2) / 2
    area_slab5 <- ((rcap + 5 * dr)^2 - (rcap + 4 * dr)^2) / 2
    area_slab6 <- ((rcap + 6 * dr)^2 - (rcap + 5 * dr)^2) / 2
    area_slab7 <- ((rcap + 7 * dr)^2 - (rcap + 6 * dr)^2) / 2
    area_slab8 <- ((rcap + 8 * dr)^2 - (rcap + 7 * dr)^2) / 2
    area_slab9 <- ((rcap + 9 * dr)^2 - (rcap + 8 * dr)^2) / 2
    area_slab10 <- ((rcap + 10 * dr)^2 - (rcap + 9 * dr)^2) / 2
    d/dt(adc_free_twoactive_slab1) <- (rcap * p_ab_h * (cp_dar2 - adc_free_twoactive_slab1 / eps) - (r_face1 * d_ab_h * (adc_free_twoactive_slab1 - adc_free_twoactive_slab2) / dr)) / area_slab1 - kon_ab_h * adc_free_twoactive_slab1 / eps * target_free_slab1 + koff_ab_h * adc_bound_twoactive_slab1
    d/dt(adc_free_oneactive_slab1) <- (rcap * p_ab_h * (cp_dar1 - adc_free_oneactive_slab1 / eps) - (r_face1 * d_ab_h * (adc_free_oneactive_slab1 - adc_free_oneactive_slab2) / dr)) / area_slab1 - kon_ab_h * adc_free_oneactive_slab1 / eps * target_free_slab1 + koff_ab_h * adc_bound_oneactive_slab1
    d/dt(adc_free_noactive_slab1) <- (rcap * p_ab_h * (cp_dar0 - adc_free_noactive_slab1 / eps) - (r_face1 * d_ab_h * (adc_free_noactive_slab1 - adc_free_noactive_slab2) / dr)) / area_slab1 - kon_ab_h * adc_free_noactive_slab1 / eps * target_free_slab1 + koff_ab_h * adc_bound_noactive_slab1
    d/dt(target_free_slab1) <- rs - kon_ab_h * (adc_free_twoactive_slab1 + adc_free_oneactive_slab1 + adc_free_noactive_slab1) / eps * target_free_slab1 + koff_ab_h * (adc_bound_twoactive_slab1 + adc_bound_oneactive_slab1 + adc_bound_noactive_slab1) - ke_ag_h * target_free_slab1
    d/dt(adc_bound_twoactive_slab1) <- kon_ab_h * adc_free_twoactive_slab1 / eps * target_free_slab1 - koff_ab_h * adc_bound_twoactive_slab1 - kint_h * adc_bound_twoactive_slab1
    d/dt(adc_bound_oneactive_slab1) <- kon_ab_h * adc_free_oneactive_slab1 / eps * target_free_slab1 - koff_ab_h * adc_bound_oneactive_slab1 - kint_h * adc_bound_oneactive_slab1
    d/dt(adc_bound_noactive_slab1) <- kon_ab_h * adc_free_noactive_slab1 / eps * target_free_slab1 - koff_ab_h * adc_bound_noactive_slab1 - kint_h * adc_bound_noactive_slab1
    d/dt(adc_internal_twoactive_slab1) <- kint_h * adc_bound_twoactive_slab1 - 2 * krel_h * adc_internal_twoactive_slab1 - kdeg_lyso_h * adc_internal_twoactive_slab1
    d/dt(adc_internal_oneactive_slab1) <- kint_h * adc_bound_oneactive_slab1 + 2 * krel_h * adc_internal_twoactive_slab1 - krel_h * adc_internal_oneactive_slab1 - kdeg_lyso_h * adc_internal_oneactive_slab1
    d/dt(adc_internal_noactive_slab1) <- kint_h * adc_bound_noactive_slab1 + krel_h * adc_internal_oneactive_slab1 - kdeg_lyso_h * adc_internal_noactive_slab1
    d/dt(payload_lyso_slab1) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab1 + kdeg_lyso_h * adc_internal_oneactive_slab1 + 2 * krel_h * adc_internal_twoactive_slab1 + krel_h * adc_internal_oneactive_slab1 - kin_pl_h * payload_lyso_slab1
    d/dt(payload_ext_slab1) <- (rcap * p_pl_h * (0 - payload_ext_slab1 / (eps * (1 + part_pl))) - (r_face1 * d_pl_h * (payload_ext_slab1 - payload_ext_slab2) / dr)) / area_slab1 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab1 + kout_pl_h * payload_int_slab1
    d/dt(payload_int_slab1) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab1 - kout_pl_h * payload_int_slab1 - kon_gq * gq_free_slab1 * payload_int_slab1 + koff_pl_h * payload_bound_slab1 + kin_pl_h * payload_lyso_slab1
    d/dt(payload_bound_slab1) <- kon_gq * gq_free_slab1 * payload_int_slab1 - koff_pl_h * payload_bound_slab1
    d/dt(gq_free_slab1) <- ksyn_gq - kdeg_gq_h * gq_free_slab1 - kon_gq * gq_free_slab1 * payload_int_slab1 + koff_pl_h * payload_bound_slab1
    target_free_slab1(0) <- ag0
    gq_free_slab1(0) <- gq0
    d/dt(adc_free_twoactive_slab2) <- (r_face1 * d_ab_h * (adc_free_twoactive_slab1 - adc_free_twoactive_slab2) / dr - (r_face2 * d_ab_h * (adc_free_twoactive_slab2 - adc_free_twoactive_slab3) / dr)) / area_slab2 - kon_ab_h * adc_free_twoactive_slab2 / eps * target_free_slab2 + koff_ab_h * adc_bound_twoactive_slab2
    d/dt(adc_free_oneactive_slab2) <- (r_face1 * d_ab_h * (adc_free_oneactive_slab1 - adc_free_oneactive_slab2) / dr - (r_face2 * d_ab_h * (adc_free_oneactive_slab2 - adc_free_oneactive_slab3) / dr)) / area_slab2 - kon_ab_h * adc_free_oneactive_slab2 / eps * target_free_slab2 + koff_ab_h * adc_bound_oneactive_slab2
    d/dt(adc_free_noactive_slab2) <- (r_face1 * d_ab_h * (adc_free_noactive_slab1 - adc_free_noactive_slab2) / dr - (r_face2 * d_ab_h * (adc_free_noactive_slab2 - adc_free_noactive_slab3) / dr)) / area_slab2 - kon_ab_h * adc_free_noactive_slab2 / eps * target_free_slab2 + koff_ab_h * adc_bound_noactive_slab2
    d/dt(target_free_slab2) <- rs - kon_ab_h * (adc_free_twoactive_slab2 + adc_free_oneactive_slab2 + adc_free_noactive_slab2) / eps * target_free_slab2 + koff_ab_h * (adc_bound_twoactive_slab2 + adc_bound_oneactive_slab2 + adc_bound_noactive_slab2) - ke_ag_h * target_free_slab2
    d/dt(adc_bound_twoactive_slab2) <- kon_ab_h * adc_free_twoactive_slab2 / eps * target_free_slab2 - koff_ab_h * adc_bound_twoactive_slab2 - kint_h * adc_bound_twoactive_slab2
    d/dt(adc_bound_oneactive_slab2) <- kon_ab_h * adc_free_oneactive_slab2 / eps * target_free_slab2 - koff_ab_h * adc_bound_oneactive_slab2 - kint_h * adc_bound_oneactive_slab2
    d/dt(adc_bound_noactive_slab2) <- kon_ab_h * adc_free_noactive_slab2 / eps * target_free_slab2 - koff_ab_h * adc_bound_noactive_slab2 - kint_h * adc_bound_noactive_slab2
    d/dt(adc_internal_twoactive_slab2) <- kint_h * adc_bound_twoactive_slab2 - 2 * krel_h * adc_internal_twoactive_slab2 - kdeg_lyso_h * adc_internal_twoactive_slab2
    d/dt(adc_internal_oneactive_slab2) <- kint_h * adc_bound_oneactive_slab2 + 2 * krel_h * adc_internal_twoactive_slab2 - krel_h * adc_internal_oneactive_slab2 - kdeg_lyso_h * adc_internal_oneactive_slab2
    d/dt(adc_internal_noactive_slab2) <- kint_h * adc_bound_noactive_slab2 + krel_h * adc_internal_oneactive_slab2 - kdeg_lyso_h * adc_internal_noactive_slab2
    d/dt(payload_lyso_slab2) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab2 + kdeg_lyso_h * adc_internal_oneactive_slab2 + 2 * krel_h * adc_internal_twoactive_slab2 + krel_h * adc_internal_oneactive_slab2 - kin_pl_h * payload_lyso_slab2
    d/dt(payload_ext_slab2) <- (r_face1 * d_pl_h * (payload_ext_slab1 - payload_ext_slab2) / dr - (r_face2 * d_pl_h * (payload_ext_slab2 - payload_ext_slab3) / dr)) / area_slab2 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab2 + kout_pl_h * payload_int_slab2
    d/dt(payload_int_slab2) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab2 - kout_pl_h * payload_int_slab2 - kon_gq * gq_free_slab2 * payload_int_slab2 + koff_pl_h * payload_bound_slab2 + kin_pl_h * payload_lyso_slab2
    d/dt(payload_bound_slab2) <- kon_gq * gq_free_slab2 * payload_int_slab2 - koff_pl_h * payload_bound_slab2
    d/dt(gq_free_slab2) <- ksyn_gq - kdeg_gq_h * gq_free_slab2 - kon_gq * gq_free_slab2 * payload_int_slab2 + koff_pl_h * payload_bound_slab2
    target_free_slab2(0) <- ag0
    gq_free_slab2(0) <- gq0
    d/dt(adc_free_twoactive_slab3) <- (r_face2 * d_ab_h * (adc_free_twoactive_slab2 - adc_free_twoactive_slab3) / dr - (r_face3 * d_ab_h * (adc_free_twoactive_slab3 - adc_free_twoactive_slab4) / dr)) / area_slab3 - kon_ab_h * adc_free_twoactive_slab3 / eps * target_free_slab3 + koff_ab_h * adc_bound_twoactive_slab3
    d/dt(adc_free_oneactive_slab3) <- (r_face2 * d_ab_h * (adc_free_oneactive_slab2 - adc_free_oneactive_slab3) / dr - (r_face3 * d_ab_h * (adc_free_oneactive_slab3 - adc_free_oneactive_slab4) / dr)) / area_slab3 - kon_ab_h * adc_free_oneactive_slab3 / eps * target_free_slab3 + koff_ab_h * adc_bound_oneactive_slab3
    d/dt(adc_free_noactive_slab3) <- (r_face2 * d_ab_h * (adc_free_noactive_slab2 - adc_free_noactive_slab3) / dr - (r_face3 * d_ab_h * (adc_free_noactive_slab3 - adc_free_noactive_slab4) / dr)) / area_slab3 - kon_ab_h * adc_free_noactive_slab3 / eps * target_free_slab3 + koff_ab_h * adc_bound_noactive_slab3
    d/dt(target_free_slab3) <- rs - kon_ab_h * (adc_free_twoactive_slab3 + adc_free_oneactive_slab3 + adc_free_noactive_slab3) / eps * target_free_slab3 + koff_ab_h * (adc_bound_twoactive_slab3 + adc_bound_oneactive_slab3 + adc_bound_noactive_slab3) - ke_ag_h * target_free_slab3
    d/dt(adc_bound_twoactive_slab3) <- kon_ab_h * adc_free_twoactive_slab3 / eps * target_free_slab3 - koff_ab_h * adc_bound_twoactive_slab3 - kint_h * adc_bound_twoactive_slab3
    d/dt(adc_bound_oneactive_slab3) <- kon_ab_h * adc_free_oneactive_slab3 / eps * target_free_slab3 - koff_ab_h * adc_bound_oneactive_slab3 - kint_h * adc_bound_oneactive_slab3
    d/dt(adc_bound_noactive_slab3) <- kon_ab_h * adc_free_noactive_slab3 / eps * target_free_slab3 - koff_ab_h * adc_bound_noactive_slab3 - kint_h * adc_bound_noactive_slab3
    d/dt(adc_internal_twoactive_slab3) <- kint_h * adc_bound_twoactive_slab3 - 2 * krel_h * adc_internal_twoactive_slab3 - kdeg_lyso_h * adc_internal_twoactive_slab3
    d/dt(adc_internal_oneactive_slab3) <- kint_h * adc_bound_oneactive_slab3 + 2 * krel_h * adc_internal_twoactive_slab3 - krel_h * adc_internal_oneactive_slab3 - kdeg_lyso_h * adc_internal_oneactive_slab3
    d/dt(adc_internal_noactive_slab3) <- kint_h * adc_bound_noactive_slab3 + krel_h * adc_internal_oneactive_slab3 - kdeg_lyso_h * adc_internal_noactive_slab3
    d/dt(payload_lyso_slab3) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab3 + kdeg_lyso_h * adc_internal_oneactive_slab3 + 2 * krel_h * adc_internal_twoactive_slab3 + krel_h * adc_internal_oneactive_slab3 - kin_pl_h * payload_lyso_slab3
    d/dt(payload_ext_slab3) <- (r_face2 * d_pl_h * (payload_ext_slab2 - payload_ext_slab3) / dr - (r_face3 * d_pl_h * (payload_ext_slab3 - payload_ext_slab4) / dr)) / area_slab3 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab3 + kout_pl_h * payload_int_slab3
    d/dt(payload_int_slab3) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab3 - kout_pl_h * payload_int_slab3 - kon_gq * gq_free_slab3 * payload_int_slab3 + koff_pl_h * payload_bound_slab3 + kin_pl_h * payload_lyso_slab3
    d/dt(payload_bound_slab3) <- kon_gq * gq_free_slab3 * payload_int_slab3 - koff_pl_h * payload_bound_slab3
    d/dt(gq_free_slab3) <- ksyn_gq - kdeg_gq_h * gq_free_slab3 - kon_gq * gq_free_slab3 * payload_int_slab3 + koff_pl_h * payload_bound_slab3
    target_free_slab3(0) <- ag0
    gq_free_slab3(0) <- gq0
    d/dt(adc_free_twoactive_slab4) <- (r_face3 * d_ab_h * (adc_free_twoactive_slab3 - adc_free_twoactive_slab4) / dr - (r_face4 * d_ab_h * (adc_free_twoactive_slab4 - adc_free_twoactive_slab5) / dr)) / area_slab4 - kon_ab_h * adc_free_twoactive_slab4 / eps * target_free_slab4 + koff_ab_h * adc_bound_twoactive_slab4
    d/dt(adc_free_oneactive_slab4) <- (r_face3 * d_ab_h * (adc_free_oneactive_slab3 - adc_free_oneactive_slab4) / dr - (r_face4 * d_ab_h * (adc_free_oneactive_slab4 - adc_free_oneactive_slab5) / dr)) / area_slab4 - kon_ab_h * adc_free_oneactive_slab4 / eps * target_free_slab4 + koff_ab_h * adc_bound_oneactive_slab4
    d/dt(adc_free_noactive_slab4) <- (r_face3 * d_ab_h * (adc_free_noactive_slab3 - adc_free_noactive_slab4) / dr - (r_face4 * d_ab_h * (adc_free_noactive_slab4 - adc_free_noactive_slab5) / dr)) / area_slab4 - kon_ab_h * adc_free_noactive_slab4 / eps * target_free_slab4 + koff_ab_h * adc_bound_noactive_slab4
    d/dt(target_free_slab4) <- rs - kon_ab_h * (adc_free_twoactive_slab4 + adc_free_oneactive_slab4 + adc_free_noactive_slab4) / eps * target_free_slab4 + koff_ab_h * (adc_bound_twoactive_slab4 + adc_bound_oneactive_slab4 + adc_bound_noactive_slab4) - ke_ag_h * target_free_slab4
    d/dt(adc_bound_twoactive_slab4) <- kon_ab_h * adc_free_twoactive_slab4 / eps * target_free_slab4 - koff_ab_h * adc_bound_twoactive_slab4 - kint_h * adc_bound_twoactive_slab4
    d/dt(adc_bound_oneactive_slab4) <- kon_ab_h * adc_free_oneactive_slab4 / eps * target_free_slab4 - koff_ab_h * adc_bound_oneactive_slab4 - kint_h * adc_bound_oneactive_slab4
    d/dt(adc_bound_noactive_slab4) <- kon_ab_h * adc_free_noactive_slab4 / eps * target_free_slab4 - koff_ab_h * adc_bound_noactive_slab4 - kint_h * adc_bound_noactive_slab4
    d/dt(adc_internal_twoactive_slab4) <- kint_h * adc_bound_twoactive_slab4 - 2 * krel_h * adc_internal_twoactive_slab4 - kdeg_lyso_h * adc_internal_twoactive_slab4
    d/dt(adc_internal_oneactive_slab4) <- kint_h * adc_bound_oneactive_slab4 + 2 * krel_h * adc_internal_twoactive_slab4 - krel_h * adc_internal_oneactive_slab4 - kdeg_lyso_h * adc_internal_oneactive_slab4
    d/dt(adc_internal_noactive_slab4) <- kint_h * adc_bound_noactive_slab4 + krel_h * adc_internal_oneactive_slab4 - kdeg_lyso_h * adc_internal_noactive_slab4
    d/dt(payload_lyso_slab4) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab4 + kdeg_lyso_h * adc_internal_oneactive_slab4 + 2 * krel_h * adc_internal_twoactive_slab4 + krel_h * adc_internal_oneactive_slab4 - kin_pl_h * payload_lyso_slab4
    d/dt(payload_ext_slab4) <- (r_face3 * d_pl_h * (payload_ext_slab3 - payload_ext_slab4) / dr - (r_face4 * d_pl_h * (payload_ext_slab4 - payload_ext_slab5) / dr)) / area_slab4 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab4 + kout_pl_h * payload_int_slab4
    d/dt(payload_int_slab4) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab4 - kout_pl_h * payload_int_slab4 - kon_gq * gq_free_slab4 * payload_int_slab4 + koff_pl_h * payload_bound_slab4 + kin_pl_h * payload_lyso_slab4
    d/dt(payload_bound_slab4) <- kon_gq * gq_free_slab4 * payload_int_slab4 - koff_pl_h * payload_bound_slab4
    d/dt(gq_free_slab4) <- ksyn_gq - kdeg_gq_h * gq_free_slab4 - kon_gq * gq_free_slab4 * payload_int_slab4 + koff_pl_h * payload_bound_slab4
    target_free_slab4(0) <- ag0
    gq_free_slab4(0) <- gq0
    d/dt(adc_free_twoactive_slab5) <- (r_face4 * d_ab_h * (adc_free_twoactive_slab4 - adc_free_twoactive_slab5) / dr - (r_face5 * d_ab_h * (adc_free_twoactive_slab5 - adc_free_twoactive_slab6) / dr)) / area_slab5 - kon_ab_h * adc_free_twoactive_slab5 / eps * target_free_slab5 + koff_ab_h * adc_bound_twoactive_slab5
    d/dt(adc_free_oneactive_slab5) <- (r_face4 * d_ab_h * (adc_free_oneactive_slab4 - adc_free_oneactive_slab5) / dr - (r_face5 * d_ab_h * (adc_free_oneactive_slab5 - adc_free_oneactive_slab6) / dr)) / area_slab5 - kon_ab_h * adc_free_oneactive_slab5 / eps * target_free_slab5 + koff_ab_h * adc_bound_oneactive_slab5
    d/dt(adc_free_noactive_slab5) <- (r_face4 * d_ab_h * (adc_free_noactive_slab4 - adc_free_noactive_slab5) / dr - (r_face5 * d_ab_h * (adc_free_noactive_slab5 - adc_free_noactive_slab6) / dr)) / area_slab5 - kon_ab_h * adc_free_noactive_slab5 / eps * target_free_slab5 + koff_ab_h * adc_bound_noactive_slab5
    d/dt(target_free_slab5) <- rs - kon_ab_h * (adc_free_twoactive_slab5 + adc_free_oneactive_slab5 + adc_free_noactive_slab5) / eps * target_free_slab5 + koff_ab_h * (adc_bound_twoactive_slab5 + adc_bound_oneactive_slab5 + adc_bound_noactive_slab5) - ke_ag_h * target_free_slab5
    d/dt(adc_bound_twoactive_slab5) <- kon_ab_h * adc_free_twoactive_slab5 / eps * target_free_slab5 - koff_ab_h * adc_bound_twoactive_slab5 - kint_h * adc_bound_twoactive_slab5
    d/dt(adc_bound_oneactive_slab5) <- kon_ab_h * adc_free_oneactive_slab5 / eps * target_free_slab5 - koff_ab_h * adc_bound_oneactive_slab5 - kint_h * adc_bound_oneactive_slab5
    d/dt(adc_bound_noactive_slab5) <- kon_ab_h * adc_free_noactive_slab5 / eps * target_free_slab5 - koff_ab_h * adc_bound_noactive_slab5 - kint_h * adc_bound_noactive_slab5
    d/dt(adc_internal_twoactive_slab5) <- kint_h * adc_bound_twoactive_slab5 - 2 * krel_h * adc_internal_twoactive_slab5 - kdeg_lyso_h * adc_internal_twoactive_slab5
    d/dt(adc_internal_oneactive_slab5) <- kint_h * adc_bound_oneactive_slab5 + 2 * krel_h * adc_internal_twoactive_slab5 - krel_h * adc_internal_oneactive_slab5 - kdeg_lyso_h * adc_internal_oneactive_slab5
    d/dt(adc_internal_noactive_slab5) <- kint_h * adc_bound_noactive_slab5 + krel_h * adc_internal_oneactive_slab5 - kdeg_lyso_h * adc_internal_noactive_slab5
    d/dt(payload_lyso_slab5) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab5 + kdeg_lyso_h * adc_internal_oneactive_slab5 + 2 * krel_h * adc_internal_twoactive_slab5 + krel_h * adc_internal_oneactive_slab5 - kin_pl_h * payload_lyso_slab5
    d/dt(payload_ext_slab5) <- (r_face4 * d_pl_h * (payload_ext_slab4 - payload_ext_slab5) / dr - (r_face5 * d_pl_h * (payload_ext_slab5 - payload_ext_slab6) / dr)) / area_slab5 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab5 + kout_pl_h * payload_int_slab5
    d/dt(payload_int_slab5) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab5 - kout_pl_h * payload_int_slab5 - kon_gq * gq_free_slab5 * payload_int_slab5 + koff_pl_h * payload_bound_slab5 + kin_pl_h * payload_lyso_slab5
    d/dt(payload_bound_slab5) <- kon_gq * gq_free_slab5 * payload_int_slab5 - koff_pl_h * payload_bound_slab5
    d/dt(gq_free_slab5) <- ksyn_gq - kdeg_gq_h * gq_free_slab5 - kon_gq * gq_free_slab5 * payload_int_slab5 + koff_pl_h * payload_bound_slab5
    target_free_slab5(0) <- ag0
    gq_free_slab5(0) <- gq0
    d/dt(adc_free_twoactive_slab6) <- (r_face5 * d_ab_h * (adc_free_twoactive_slab5 - adc_free_twoactive_slab6) / dr - (r_face6 * d_ab_h * (adc_free_twoactive_slab6 - adc_free_twoactive_slab7) / dr)) / area_slab6 - kon_ab_h * adc_free_twoactive_slab6 / eps * target_free_slab6 + koff_ab_h * adc_bound_twoactive_slab6
    d/dt(adc_free_oneactive_slab6) <- (r_face5 * d_ab_h * (adc_free_oneactive_slab5 - adc_free_oneactive_slab6) / dr - (r_face6 * d_ab_h * (adc_free_oneactive_slab6 - adc_free_oneactive_slab7) / dr)) / area_slab6 - kon_ab_h * adc_free_oneactive_slab6 / eps * target_free_slab6 + koff_ab_h * adc_bound_oneactive_slab6
    d/dt(adc_free_noactive_slab6) <- (r_face5 * d_ab_h * (adc_free_noactive_slab5 - adc_free_noactive_slab6) / dr - (r_face6 * d_ab_h * (adc_free_noactive_slab6 - adc_free_noactive_slab7) / dr)) / area_slab6 - kon_ab_h * adc_free_noactive_slab6 / eps * target_free_slab6 + koff_ab_h * adc_bound_noactive_slab6
    d/dt(target_free_slab6) <- rs - kon_ab_h * (adc_free_twoactive_slab6 + adc_free_oneactive_slab6 + adc_free_noactive_slab6) / eps * target_free_slab6 + koff_ab_h * (adc_bound_twoactive_slab6 + adc_bound_oneactive_slab6 + adc_bound_noactive_slab6) - ke_ag_h * target_free_slab6
    d/dt(adc_bound_twoactive_slab6) <- kon_ab_h * adc_free_twoactive_slab6 / eps * target_free_slab6 - koff_ab_h * adc_bound_twoactive_slab6 - kint_h * adc_bound_twoactive_slab6
    d/dt(adc_bound_oneactive_slab6) <- kon_ab_h * adc_free_oneactive_slab6 / eps * target_free_slab6 - koff_ab_h * adc_bound_oneactive_slab6 - kint_h * adc_bound_oneactive_slab6
    d/dt(adc_bound_noactive_slab6) <- kon_ab_h * adc_free_noactive_slab6 / eps * target_free_slab6 - koff_ab_h * adc_bound_noactive_slab6 - kint_h * adc_bound_noactive_slab6
    d/dt(adc_internal_twoactive_slab6) <- kint_h * adc_bound_twoactive_slab6 - 2 * krel_h * adc_internal_twoactive_slab6 - kdeg_lyso_h * adc_internal_twoactive_slab6
    d/dt(adc_internal_oneactive_slab6) <- kint_h * adc_bound_oneactive_slab6 + 2 * krel_h * adc_internal_twoactive_slab6 - krel_h * adc_internal_oneactive_slab6 - kdeg_lyso_h * adc_internal_oneactive_slab6
    d/dt(adc_internal_noactive_slab6) <- kint_h * adc_bound_noactive_slab6 + krel_h * adc_internal_oneactive_slab6 - kdeg_lyso_h * adc_internal_noactive_slab6
    d/dt(payload_lyso_slab6) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab6 + kdeg_lyso_h * adc_internal_oneactive_slab6 + 2 * krel_h * adc_internal_twoactive_slab6 + krel_h * adc_internal_oneactive_slab6 - kin_pl_h * payload_lyso_slab6
    d/dt(payload_ext_slab6) <- (r_face5 * d_pl_h * (payload_ext_slab5 - payload_ext_slab6) / dr - (r_face6 * d_pl_h * (payload_ext_slab6 - payload_ext_slab7) / dr)) / area_slab6 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab6 + kout_pl_h * payload_int_slab6
    d/dt(payload_int_slab6) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab6 - kout_pl_h * payload_int_slab6 - kon_gq * gq_free_slab6 * payload_int_slab6 + koff_pl_h * payload_bound_slab6 + kin_pl_h * payload_lyso_slab6
    d/dt(payload_bound_slab6) <- kon_gq * gq_free_slab6 * payload_int_slab6 - koff_pl_h * payload_bound_slab6
    d/dt(gq_free_slab6) <- ksyn_gq - kdeg_gq_h * gq_free_slab6 - kon_gq * gq_free_slab6 * payload_int_slab6 + koff_pl_h * payload_bound_slab6
    target_free_slab6(0) <- ag0
    gq_free_slab6(0) <- gq0
    d/dt(adc_free_twoactive_slab7) <- (r_face6 * d_ab_h * (adc_free_twoactive_slab6 - adc_free_twoactive_slab7) / dr - (r_face7 * d_ab_h * (adc_free_twoactive_slab7 - adc_free_twoactive_slab8) / dr)) / area_slab7 - kon_ab_h * adc_free_twoactive_slab7 / eps * target_free_slab7 + koff_ab_h * adc_bound_twoactive_slab7
    d/dt(adc_free_oneactive_slab7) <- (r_face6 * d_ab_h * (adc_free_oneactive_slab6 - adc_free_oneactive_slab7) / dr - (r_face7 * d_ab_h * (adc_free_oneactive_slab7 - adc_free_oneactive_slab8) / dr)) / area_slab7 - kon_ab_h * adc_free_oneactive_slab7 / eps * target_free_slab7 + koff_ab_h * adc_bound_oneactive_slab7
    d/dt(adc_free_noactive_slab7) <- (r_face6 * d_ab_h * (adc_free_noactive_slab6 - adc_free_noactive_slab7) / dr - (r_face7 * d_ab_h * (adc_free_noactive_slab7 - adc_free_noactive_slab8) / dr)) / area_slab7 - kon_ab_h * adc_free_noactive_slab7 / eps * target_free_slab7 + koff_ab_h * adc_bound_noactive_slab7
    d/dt(target_free_slab7) <- rs - kon_ab_h * (adc_free_twoactive_slab7 + adc_free_oneactive_slab7 + adc_free_noactive_slab7) / eps * target_free_slab7 + koff_ab_h * (adc_bound_twoactive_slab7 + adc_bound_oneactive_slab7 + adc_bound_noactive_slab7) - ke_ag_h * target_free_slab7
    d/dt(adc_bound_twoactive_slab7) <- kon_ab_h * adc_free_twoactive_slab7 / eps * target_free_slab7 - koff_ab_h * adc_bound_twoactive_slab7 - kint_h * adc_bound_twoactive_slab7
    d/dt(adc_bound_oneactive_slab7) <- kon_ab_h * adc_free_oneactive_slab7 / eps * target_free_slab7 - koff_ab_h * adc_bound_oneactive_slab7 - kint_h * adc_bound_oneactive_slab7
    d/dt(adc_bound_noactive_slab7) <- kon_ab_h * adc_free_noactive_slab7 / eps * target_free_slab7 - koff_ab_h * adc_bound_noactive_slab7 - kint_h * adc_bound_noactive_slab7
    d/dt(adc_internal_twoactive_slab7) <- kint_h * adc_bound_twoactive_slab7 - 2 * krel_h * adc_internal_twoactive_slab7 - kdeg_lyso_h * adc_internal_twoactive_slab7
    d/dt(adc_internal_oneactive_slab7) <- kint_h * adc_bound_oneactive_slab7 + 2 * krel_h * adc_internal_twoactive_slab7 - krel_h * adc_internal_oneactive_slab7 - kdeg_lyso_h * adc_internal_oneactive_slab7
    d/dt(adc_internal_noactive_slab7) <- kint_h * adc_bound_noactive_slab7 + krel_h * adc_internal_oneactive_slab7 - kdeg_lyso_h * adc_internal_noactive_slab7
    d/dt(payload_lyso_slab7) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab7 + kdeg_lyso_h * adc_internal_oneactive_slab7 + 2 * krel_h * adc_internal_twoactive_slab7 + krel_h * adc_internal_oneactive_slab7 - kin_pl_h * payload_lyso_slab7
    d/dt(payload_ext_slab7) <- (r_face6 * d_pl_h * (payload_ext_slab6 - payload_ext_slab7) / dr - (r_face7 * d_pl_h * (payload_ext_slab7 - payload_ext_slab8) / dr)) / area_slab7 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab7 + kout_pl_h * payload_int_slab7
    d/dt(payload_int_slab7) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab7 - kout_pl_h * payload_int_slab7 - kon_gq * gq_free_slab7 * payload_int_slab7 + koff_pl_h * payload_bound_slab7 + kin_pl_h * payload_lyso_slab7
    d/dt(payload_bound_slab7) <- kon_gq * gq_free_slab7 * payload_int_slab7 - koff_pl_h * payload_bound_slab7
    d/dt(gq_free_slab7) <- ksyn_gq - kdeg_gq_h * gq_free_slab7 - kon_gq * gq_free_slab7 * payload_int_slab7 + koff_pl_h * payload_bound_slab7
    target_free_slab7(0) <- ag0
    gq_free_slab7(0) <- gq0
    d/dt(adc_free_twoactive_slab8) <- (r_face7 * d_ab_h * (adc_free_twoactive_slab7 - adc_free_twoactive_slab8) / dr - (r_face8 * d_ab_h * (adc_free_twoactive_slab8 - adc_free_twoactive_slab9) / dr)) / area_slab8 - kon_ab_h * adc_free_twoactive_slab8 / eps * target_free_slab8 + koff_ab_h * adc_bound_twoactive_slab8
    d/dt(adc_free_oneactive_slab8) <- (r_face7 * d_ab_h * (adc_free_oneactive_slab7 - adc_free_oneactive_slab8) / dr - (r_face8 * d_ab_h * (adc_free_oneactive_slab8 - adc_free_oneactive_slab9) / dr)) / area_slab8 - kon_ab_h * adc_free_oneactive_slab8 / eps * target_free_slab8 + koff_ab_h * adc_bound_oneactive_slab8
    d/dt(adc_free_noactive_slab8) <- (r_face7 * d_ab_h * (adc_free_noactive_slab7 - adc_free_noactive_slab8) / dr - (r_face8 * d_ab_h * (adc_free_noactive_slab8 - adc_free_noactive_slab9) / dr)) / area_slab8 - kon_ab_h * adc_free_noactive_slab8 / eps * target_free_slab8 + koff_ab_h * adc_bound_noactive_slab8
    d/dt(target_free_slab8) <- rs - kon_ab_h * (adc_free_twoactive_slab8 + adc_free_oneactive_slab8 + adc_free_noactive_slab8) / eps * target_free_slab8 + koff_ab_h * (adc_bound_twoactive_slab8 + adc_bound_oneactive_slab8 + adc_bound_noactive_slab8) - ke_ag_h * target_free_slab8
    d/dt(adc_bound_twoactive_slab8) <- kon_ab_h * adc_free_twoactive_slab8 / eps * target_free_slab8 - koff_ab_h * adc_bound_twoactive_slab8 - kint_h * adc_bound_twoactive_slab8
    d/dt(adc_bound_oneactive_slab8) <- kon_ab_h * adc_free_oneactive_slab8 / eps * target_free_slab8 - koff_ab_h * adc_bound_oneactive_slab8 - kint_h * adc_bound_oneactive_slab8
    d/dt(adc_bound_noactive_slab8) <- kon_ab_h * adc_free_noactive_slab8 / eps * target_free_slab8 - koff_ab_h * adc_bound_noactive_slab8 - kint_h * adc_bound_noactive_slab8
    d/dt(adc_internal_twoactive_slab8) <- kint_h * adc_bound_twoactive_slab8 - 2 * krel_h * adc_internal_twoactive_slab8 - kdeg_lyso_h * adc_internal_twoactive_slab8
    d/dt(adc_internal_oneactive_slab8) <- kint_h * adc_bound_oneactive_slab8 + 2 * krel_h * adc_internal_twoactive_slab8 - krel_h * adc_internal_oneactive_slab8 - kdeg_lyso_h * adc_internal_oneactive_slab8
    d/dt(adc_internal_noactive_slab8) <- kint_h * adc_bound_noactive_slab8 + krel_h * adc_internal_oneactive_slab8 - kdeg_lyso_h * adc_internal_noactive_slab8
    d/dt(payload_lyso_slab8) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab8 + kdeg_lyso_h * adc_internal_oneactive_slab8 + 2 * krel_h * adc_internal_twoactive_slab8 + krel_h * adc_internal_oneactive_slab8 - kin_pl_h * payload_lyso_slab8
    d/dt(payload_ext_slab8) <- (r_face7 * d_pl_h * (payload_ext_slab7 - payload_ext_slab8) / dr - (r_face8 * d_pl_h * (payload_ext_slab8 - payload_ext_slab9) / dr)) / area_slab8 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab8 + kout_pl_h * payload_int_slab8
    d/dt(payload_int_slab8) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab8 - kout_pl_h * payload_int_slab8 - kon_gq * gq_free_slab8 * payload_int_slab8 + koff_pl_h * payload_bound_slab8 + kin_pl_h * payload_lyso_slab8
    d/dt(payload_bound_slab8) <- kon_gq * gq_free_slab8 * payload_int_slab8 - koff_pl_h * payload_bound_slab8
    d/dt(gq_free_slab8) <- ksyn_gq - kdeg_gq_h * gq_free_slab8 - kon_gq * gq_free_slab8 * payload_int_slab8 + koff_pl_h * payload_bound_slab8
    target_free_slab8(0) <- ag0
    gq_free_slab8(0) <- gq0
    d/dt(adc_free_twoactive_slab9) <- (r_face8 * d_ab_h * (adc_free_twoactive_slab8 - adc_free_twoactive_slab9) / dr - (r_face9 * d_ab_h * (adc_free_twoactive_slab9 - adc_free_twoactive_slab10) / dr)) / area_slab9 - kon_ab_h * adc_free_twoactive_slab9 / eps * target_free_slab9 + koff_ab_h * adc_bound_twoactive_slab9
    d/dt(adc_free_oneactive_slab9) <- (r_face8 * d_ab_h * (adc_free_oneactive_slab8 - adc_free_oneactive_slab9) / dr - (r_face9 * d_ab_h * (adc_free_oneactive_slab9 - adc_free_oneactive_slab10) / dr)) / area_slab9 - kon_ab_h * adc_free_oneactive_slab9 / eps * target_free_slab9 + koff_ab_h * adc_bound_oneactive_slab9
    d/dt(adc_free_noactive_slab9) <- (r_face8 * d_ab_h * (adc_free_noactive_slab8 - adc_free_noactive_slab9) / dr - (r_face9 * d_ab_h * (adc_free_noactive_slab9 - adc_free_noactive_slab10) / dr)) / area_slab9 - kon_ab_h * adc_free_noactive_slab9 / eps * target_free_slab9 + koff_ab_h * adc_bound_noactive_slab9
    d/dt(target_free_slab9) <- rs - kon_ab_h * (adc_free_twoactive_slab9 + adc_free_oneactive_slab9 + adc_free_noactive_slab9) / eps * target_free_slab9 + koff_ab_h * (adc_bound_twoactive_slab9 + adc_bound_oneactive_slab9 + adc_bound_noactive_slab9) - ke_ag_h * target_free_slab9
    d/dt(adc_bound_twoactive_slab9) <- kon_ab_h * adc_free_twoactive_slab9 / eps * target_free_slab9 - koff_ab_h * adc_bound_twoactive_slab9 - kint_h * adc_bound_twoactive_slab9
    d/dt(adc_bound_oneactive_slab9) <- kon_ab_h * adc_free_oneactive_slab9 / eps * target_free_slab9 - koff_ab_h * adc_bound_oneactive_slab9 - kint_h * adc_bound_oneactive_slab9
    d/dt(adc_bound_noactive_slab9) <- kon_ab_h * adc_free_noactive_slab9 / eps * target_free_slab9 - koff_ab_h * adc_bound_noactive_slab9 - kint_h * adc_bound_noactive_slab9
    d/dt(adc_internal_twoactive_slab9) <- kint_h * adc_bound_twoactive_slab9 - 2 * krel_h * adc_internal_twoactive_slab9 - kdeg_lyso_h * adc_internal_twoactive_slab9
    d/dt(adc_internal_oneactive_slab9) <- kint_h * adc_bound_oneactive_slab9 + 2 * krel_h * adc_internal_twoactive_slab9 - krel_h * adc_internal_oneactive_slab9 - kdeg_lyso_h * adc_internal_oneactive_slab9
    d/dt(adc_internal_noactive_slab9) <- kint_h * adc_bound_noactive_slab9 + krel_h * adc_internal_oneactive_slab9 - kdeg_lyso_h * adc_internal_noactive_slab9
    d/dt(payload_lyso_slab9) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab9 + kdeg_lyso_h * adc_internal_oneactive_slab9 + 2 * krel_h * adc_internal_twoactive_slab9 + krel_h * adc_internal_oneactive_slab9 - kin_pl_h * payload_lyso_slab9
    d/dt(payload_ext_slab9) <- (r_face8 * d_pl_h * (payload_ext_slab8 - payload_ext_slab9) / dr - (r_face9 * d_pl_h * (payload_ext_slab9 - payload_ext_slab10) / dr)) / area_slab9 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab9 + kout_pl_h * payload_int_slab9
    d/dt(payload_int_slab9) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab9 - kout_pl_h * payload_int_slab9 - kon_gq * gq_free_slab9 * payload_int_slab9 + koff_pl_h * payload_bound_slab9 + kin_pl_h * payload_lyso_slab9
    d/dt(payload_bound_slab9) <- kon_gq * gq_free_slab9 * payload_int_slab9 - koff_pl_h * payload_bound_slab9
    d/dt(gq_free_slab9) <- ksyn_gq - kdeg_gq_h * gq_free_slab9 - kon_gq * gq_free_slab9 * payload_int_slab9 + koff_pl_h * payload_bound_slab9
    target_free_slab9(0) <- ag0
    gq_free_slab9(0) <- gq0
    d/dt(adc_free_twoactive_slab10) <- (r_face9 * d_ab_h * (adc_free_twoactive_slab9 - adc_free_twoactive_slab10) / dr - (0)) / area_slab10 - kon_ab_h * adc_free_twoactive_slab10 / eps * target_free_slab10 + koff_ab_h * adc_bound_twoactive_slab10
    d/dt(adc_free_oneactive_slab10) <- (r_face9 * d_ab_h * (adc_free_oneactive_slab9 - adc_free_oneactive_slab10) / dr - (0)) / area_slab10 - kon_ab_h * adc_free_oneactive_slab10 / eps * target_free_slab10 + koff_ab_h * adc_bound_oneactive_slab10
    d/dt(adc_free_noactive_slab10) <- (r_face9 * d_ab_h * (adc_free_noactive_slab9 - adc_free_noactive_slab10) / dr - (0)) / area_slab10 - kon_ab_h * adc_free_noactive_slab10 / eps * target_free_slab10 + koff_ab_h * adc_bound_noactive_slab10
    d/dt(target_free_slab10) <- rs - kon_ab_h * (adc_free_twoactive_slab10 + adc_free_oneactive_slab10 + adc_free_noactive_slab10) / eps * target_free_slab10 + koff_ab_h * (adc_bound_twoactive_slab10 + adc_bound_oneactive_slab10 + adc_bound_noactive_slab10) - ke_ag_h * target_free_slab10
    d/dt(adc_bound_twoactive_slab10) <- kon_ab_h * adc_free_twoactive_slab10 / eps * target_free_slab10 - koff_ab_h * adc_bound_twoactive_slab10 - kint_h * adc_bound_twoactive_slab10
    d/dt(adc_bound_oneactive_slab10) <- kon_ab_h * adc_free_oneactive_slab10 / eps * target_free_slab10 - koff_ab_h * adc_bound_oneactive_slab10 - kint_h * adc_bound_oneactive_slab10
    d/dt(adc_bound_noactive_slab10) <- kon_ab_h * adc_free_noactive_slab10 / eps * target_free_slab10 - koff_ab_h * adc_bound_noactive_slab10 - kint_h * adc_bound_noactive_slab10
    d/dt(adc_internal_twoactive_slab10) <- kint_h * adc_bound_twoactive_slab10 - 2 * krel_h * adc_internal_twoactive_slab10 - kdeg_lyso_h * adc_internal_twoactive_slab10
    d/dt(adc_internal_oneactive_slab10) <- kint_h * adc_bound_oneactive_slab10 + 2 * krel_h * adc_internal_twoactive_slab10 - krel_h * adc_internal_oneactive_slab10 - kdeg_lyso_h * adc_internal_oneactive_slab10
    d/dt(adc_internal_noactive_slab10) <- kint_h * adc_bound_noactive_slab10 + krel_h * adc_internal_oneactive_slab10 - kdeg_lyso_h * adc_internal_noactive_slab10
    d/dt(payload_lyso_slab10) <- 2 * kdeg_lyso_h * adc_internal_twoactive_slab10 + kdeg_lyso_h * adc_internal_oneactive_slab10 + 2 * krel_h * adc_internal_twoactive_slab10 + krel_h * adc_internal_oneactive_slab10 - kin_pl_h * payload_lyso_slab10
    d/dt(payload_ext_slab10) <- (r_face9 * d_pl_h * (payload_ext_slab9 - payload_ext_slab10) / dr - (0)) / area_slab10 - kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab10 + kout_pl_h * payload_int_slab10
    d/dt(payload_int_slab10) <- kin_pl_h * (1 - eps_pl) / eps_pl * payload_ext_slab10 - kout_pl_h * payload_int_slab10 - kon_gq * gq_free_slab10 * payload_int_slab10 + koff_pl_h * payload_bound_slab10 + kin_pl_h * payload_lyso_slab10
    d/dt(payload_bound_slab10) <- kon_gq * gq_free_slab10 * payload_int_slab10 - koff_pl_h * payload_bound_slab10
    d/dt(gq_free_slab10) <- ksyn_gq - kdeg_gq_h * gq_free_slab10 - kon_gq * gq_free_slab10 * payload_int_slab10 + koff_pl_h * payload_bound_slab10
    target_free_slab10(0) <- ag0
    gq_free_slab10(0) <- gq0

    # ---- Outputs ------------------------------------------------------------
    # Radial averages over the equal-width shells: free Galpha-q/11 as a
    # fraction of baseline (Figure 4B, 'Predicted Galpha-q/11 free'),
    # PMEL17 occupancy by all ADC species as a fraction of [Ag]0, and total
    # intracellular active payload (free + Galpha-q-bound + lysosomal), nM.
    gq_free_frac <- (gq_free_slab1 + gq_free_slab2 + gq_free_slab3 + gq_free_slab4 + gq_free_slab5 + gq_free_slab6 + gq_free_slab7 + gq_free_slab8 + gq_free_slab9 + gq_free_slab10) / (10 * gq0)
    pmel17_occupancy <- (adc_bound_twoactive_slab1 + adc_bound_twoactive_slab2 + adc_bound_twoactive_slab3 + adc_bound_twoactive_slab4 + adc_bound_twoactive_slab5 + adc_bound_twoactive_slab6 + adc_bound_twoactive_slab7 + adc_bound_twoactive_slab8 + adc_bound_twoactive_slab9 + adc_bound_twoactive_slab10 + adc_bound_oneactive_slab1 + adc_bound_oneactive_slab2 + adc_bound_oneactive_slab3 + adc_bound_oneactive_slab4 + adc_bound_oneactive_slab5 + adc_bound_oneactive_slab6 + adc_bound_oneactive_slab7 + adc_bound_oneactive_slab8 + adc_bound_oneactive_slab9 + adc_bound_oneactive_slab10 + adc_bound_noactive_slab1 + adc_bound_noactive_slab2 + adc_bound_noactive_slab3 + adc_bound_noactive_slab4 + adc_bound_noactive_slab5 + adc_bound_noactive_slab6 + adc_bound_noactive_slab7 + adc_bound_noactive_slab8 + adc_bound_noactive_slab9 + adc_bound_noactive_slab10) / (10 * ag0)
    payload_cell <- (payload_int_slab1 + payload_int_slab2 + payload_int_slab3 + payload_int_slab4 + payload_int_slab5 + payload_int_slab6 + payload_int_slab7 + payload_int_slab8 + payload_int_slab9 + payload_int_slab10 + payload_bound_slab1 + payload_bound_slab2 + payload_bound_slab3 + payload_bound_slab4 + payload_bound_slab5 + payload_bound_slab6 + payload_bound_slab7 + payload_bound_slab8 + payload_bound_slab9 + payload_bound_slab10 + payload_lyso_slab1 + payload_lyso_slab2 + payload_lyso_slab3 + payload_lyso_slab4 + payload_lyso_slab5 + payload_lyso_slab6 + payload_lyso_slab7 + payload_lyso_slab8 + payload_lyso_slab9 + payload_lyso_slab10) / 10
  })
}
