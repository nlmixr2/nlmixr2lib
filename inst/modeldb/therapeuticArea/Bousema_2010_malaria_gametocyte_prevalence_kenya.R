Bousema_2010_malaria_gametocyte_prevalence_kenya <- function() {
  description <- paste(
    "Two-state gametocyte-prevalence model of Plasmodium falciparum",
    "gametocyte carriage after efficacious antimalarial treatment in Kenyan",
    "children with uncomplicated malaria (Bousema 2010, Kenya trial",
    "2003-2004, Pfs25 QT-NASBA data). A patient who is gametocyte-negative",
    "(E) becomes gametocyte-positive (I) at rate f = 1/(mean time to",
    "positivity) as sequestered gametocytes emerge, and a positive patient",
    "clears gametocytaemia at rate r = 1/(mean duration of gametocyte",
    "carriage). Time zero is the start of treatment (day 0); the chain",
    "starts from the patient's enrolment status GAMETOCYTE_POS_BASE. The",
    "duration of carriage differs between the non-ACT arm",
    "(sulfadoxine-pyrimethamine + amodiaquine; reference) and the ACT arm",
    "(sulfadoxine-pyrimethamine + artesunate, or artemether-lumefantrine),",
    "selected by TRT_ARTEMISININ. Observations are Bernoulli",
    "(gametocyte_status) with a Normal patient random effect on the logit",
    "of the probability of positivity; its variance is not reported and is",
    "fixed to 0, so simulations are typical-value. Sibling models: the",
    "Tanzania prevalence fit",
    "(Bousema_2010_malaria_gametocyte_prevalence_tanzania) and the density",
    "models (Bousema_2010_malaria_gametocyte_density_kenya,",
    "Bousema_2010_malaria_gametocyte_density_tanzania).",
    sep = " "
  )
  reference <- paste(
    "Bousema T, Okell L, Shekalaghe S, Griffin JT, Omar S, Sawa P,",
    "Sutherland C, Sauerwein R, Ghani AC, Drakeley C (2010). Revisiting the",
    "circulation time of Plasmodium falciparum gametocytes: molecular",
    "detection methods to estimate the duration of gametocyte carriage and",
    "the effect of gametocytocidal drugs. Malar J 9:136.",
    "doi:10.1186/1475-2875-9-136.",
    sep = " "
  )
  vignette <- "Bousema_2010_malaria_gametocyte"
  units <- list(
    time = "day (time zero = start of treatment)",
    dosing = "no dosing (treatment arm enters through TRT_ARTEMISININ)",
    concentration = "probability (unitless)"
  )

  # The patient random effect is additive on the logit of the probability of
  # gametocyte positivity; it has no typical-value theta of its own.
  paper_specific_etas <- c("etalogitgam")

  compartmentData <- list(
    s_gametocyteneg = list(
      analyte = "Probability that the patient is infected but gametocyte-negative by Pfs25 QT-NASBA (gametocytes not yet released into the circulation)",
      units = "(fraction)",
      specimen = "whole blood",
      verified = TRUE
    ),
    s_gametocytepos = list(
      analyte = "Probability that the patient is gametocyte-positive by Pfs25 QT-NASBA",
      units = "(fraction)",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    TRT_ARTEMISININ = list(
      description = "Artemisinin-based combination therapy (ACT) arm indicator: 1 = sulfadoxine-pyrimethamine + artesunate or artemether-lumefantrine, 0 = sulfadoxine-pyrimethamine + amodiaquine (non-ACT, reference)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-ACT (sulfadoxine-pyrimethamine + amodiaquine)",
      notes = paste(
        "Bousema 2010 Data analysis pools the SP+AS and AL arms as one ACT",
        "group because both had an indistinguishable effect on gametocyte",
        "carriage. Selects the arm-specific mean duration of carriage."
      ),
      source_name = "treatment group (non-ACT / ACT)"
    ),
    GAMETOCYTE_POS_BASE = list(
      description = "Gametocyte-positive by Pfs25 QT-NASBA at enrolment (day 0): 1 = positive, 0 = negative",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = gametocyte-negative at enrolment",
      notes = paste(
        "Sets the initial state of the two-state chain:",
        "s_gametocytepos(0) = GAMETOCYTE_POS_BASE (the paper's I0). Supplying",
        "a cohort proportion instead of a 0/1 value gives the paper's",
        "Equation for the expected prevalence I(t) of that cohort."
      ),
      source_name = "I0 (prevalence at day 0)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 126L,
    n_studies = 1L,
    age_range = "0.5-10 years (median 3, IQR 1-5; Table 1)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria, 500-100,000 asexual",
      "parasites/uL at enrolment, Hb >= 5 g/dL; successfully cured of",
      "asexual parasites, no re-infection, and Pfs25 QT-NASBA gametocyte",
      "positive at least once between day 0 and day 28."
    ),
    dose_range = paste(
      "Non-ACT: sulfadoxine-pyrimethamine (25 mg/kg S + 1.25 mg/kg P single",
      "dose) + amodiaquine 10 mg/kg once daily for 3 days. ACT:",
      "sulfadoxine-pyrimethamine + artesunate 4 mg/kg once daily for 3",
      "days, or artemether-lumefantrine (10 mg A + 60 mg L per 5 kg twice",
      "daily for 3 days)."
    ),
    regions = "Mbita, western Kenya (EIR ~70 infectious bites per person per year)",
    notes = paste(
      "36 non-ACT and 90 ACT children were included in the model fitting",
      "(Table 1); 4 and 12 of them were gametocyte-negative at enrolment",
      "(Table 2 footnote). Gametocyte status was measured by Pfs25 QT-NASBA",
      "on days 0, 3, 7, 14 and 28. Fitted with SAS PROC NLMIXED."
    )
  )

  ini({
    lmtpos <- log(5.60); label("Mean time until a gametocyte-negative patient becomes gametocyte-positive, 1/f (day)") # Table 2 '1/f' = 5.60 days (95% CI 1.06-29.4), estimated from the Kenya trial (footnote)

    lmdur_nonact <- log(55.6); label("Mean duration of gametocyte carriage, 1/r, non-ACT arm (day)") # Table 3 non-ACT 55.6 (28.7-107.7); Figure 4 caption
    lmdur_act <- log(13.4); label("Mean duration of gametocyte carriage, 1/r, ACT arm (day)") # Table 3 ACT 13.4 (10.2-17.5); Figure 4 caption

    etalogitgam ~ fixed(0) # Methods: Normal patient random effect on the logit of gametocyte positivity; variance not reported
  })

  model({
    f <- 1 / exp(lmtpos)
    r <- 1 / exp(lmdur_nonact * (1 - TRT_ARTEMISININ) + lmdur_act * TRT_ARTEMISININ)

    s_gametocyteneg(0) <- 1 - GAMETOCYTE_POS_BASE
    s_gametocytepos(0) <- GAMETOCYTE_POS_BASE

    d/dt(s_gametocyteneg) <- -f * s_gametocyteneg
    d/dt(s_gametocytepos) <- f * s_gametocyteneg - r * s_gametocytepos

    # expit(logit(I) + eta) written so that it stays finite at I = 0 or 1.
    p_gametocytepos <- s_gametocytepos * exp(etalogitgam) /
      (1 - s_gametocytepos + s_gametocytepos * exp(etalogitgam))
    p0_gametocyte <- 1 - p_gametocytepos
    gametocyte_status ~ c(p0_gametocyte)
  })
}
