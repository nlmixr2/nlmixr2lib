# nlmixr2lib

# development version

- Add Wang 2023 zoledronic acid ([doi:10.3389/fphar.2023.1089774](https://doi.org/10.3389/fphar.2023.1089774)) -- adults with primary osteoporosis, pooled from 10 published trials.

- Add Ousey 2026 plozasiran ([doi:10.1002/jcph.70190](https://doi.org/10.1002/jcph.70190)) -- adults with familial chylomicronemia syndrome.

- Add Ravva 2010 varenicline exposure-response models ([doi:10.1038/clpt.2009.282](https://doi.org/10.1038/clpt.2009.282)) -- adult cigarette smokers in smoking-cessation trials.

- Rename the sigmoidal-in-time clearance parameter family so the names state the
  structure rather than only the curve shape, across all 20 models that use it:
  `cl_hill_max` -> `cl_time_max`, `cl_hill_t50` -> `cl_t50`, `cl_hill_gamma` ->
  `cl_time_hill` (with matching `l`-prefixed, `eta`-prefixed and `_i` forms). The
  family now extends the registered `cl_time` stem, and `hill` survives only on the
  shape coefficient where it is meaningful. `checkModelConventions()` and its test
  were updated to recognise the new names. The trio had been in use for ~20 models
  with no register entry; all three are now documented in `parameter-names.md`.
  Unrelated maturation Hill coefficients (`e_page_cl_hill`, `e_pna_cl_hill`) are
  deliberately unchanged.

- Add Riggs 2012 CKD-mineral bone disorder multiscale QSP model
  ([doi:10.1177/0091270011412967](https://doi.org/10.1177/0091270011412967)) —
  a hypothetical adult over a 10-year course of chronic kidney disease.

- Add Riggs 2013 empagliflozin ([doi:10.1002/jcph.147](https://doi.org/10.1002/jcph.147)) - adults with type 2 diabetes mellitus.

- Add Salem 2013 drug-elimination-pathway ontogeny ([doi:10.1002/jcph.100](https://doi.org/10.1002/jcph.100)) -- healthy humans from birth to 20 years.

- Add Zhang 2016 burosumab ([doi:10.1002/jcph.611](https://doi.org/10.1002/jcph.611)) -- adults with X-linked hypophosphatemia.

- Add Kim 2026 zolpidem ([doi:10.1002/psp4.70208](https://doi.org/10.1002/psp4.70208)) -- healthy Korean volunteers; PK plus three sequential PK-PD models for DSST, choice reaction time and sedation VAS.

- Add Songvut 2026 andrographolide ([doi:10.3389/fphar.2026.1781740](https://doi.org/10.3389/fphar.2026.1781740)) -- adults with mild COVID-19.

- Add Dings 2026 cafedrine/theodrenaline and ephedrine ([doi:10.3390/pharmaceutics18030296](https://doi.org/10.3390/pharmaceutics18030296)) - parturients with spinal-anaesthesia-induced hypotension during caesarean section, plus a companion neonatal acid-base outcome model.

- Add Saporta 2026 meropenem ([doi:10.1128/aac.01788-25](https://doi.org/10.1128/aac.01788-25)) - neutropenic, intermediately suppressed and immunocompetent mice with Klebsiella pneumoniae lung infection.

- Add Bae 2026 MIT-001 ([doi:10.1002/jcph.70189](https://doi.org/10.1002/jcph.70189)) -- healthy Korean adults.

- Add Bamgboye 2026 intravenous topiramate ([doi:10.1002/jcph.70191](https://doi.org/10.1002/jcph.70191)) — adults with epilepsy or migraine.

- Add Bhagunde 2026 lecanemab plasma-biomarker PK/PD models ([doi:10.1002/trc2.70246](https://doi.org/10.1002/trc2.70246)) -- adults with early Alzheimer's disease.

- Add Cao 2026 lecanemab ([doi:10.1038/s41540-026-00677-4](https://doi.org/10.1038/s41540-026-00677-4)) — adults with early symptomatic Alzheimer's disease.

- Add Hu 2026 pemetrexed + osimertinib sequence-dependent synergy QSP model ([doi:10.3390/pharmaceutics18040408](https://doi.org/10.3390/pharmaceutics18040408)) -- HCC827 xenograft-bearing female BALB/c nude mice.

- Add Lim 2026 leuprolide acetate 3-month depot ([doi:10.1007/s40272-025-00733-2](https://doi.org/10.1007/s40272-025-00733-2)) -- children aged 1-10 years with central precocious puberty.

- Add Gafar 2026 rifampicin ([doi:10.1093/infdis/jiag052](https://doi.org/10.1093/infdis/jiag052)) -- adolescents and adults receiving tuberculosis preventive therapy in Canada, Indonesia and Vietnam.

- Add Wei 2026 enrofloxacin + cefquinome ([doi:10.1186/s12917-026-05301-5](https://doi.org/10.1186/s12917-026-05301-5)) -- chicks with experimentally induced *Klebsiella pneumoniae* pneumonia.

- Add Tang 2026 vixarelimab ([doi:10.1002/psp4.70230](https://doi.org/10.1002/psp4.70230)) -- healthy volunteers and patients with chronic pruritic conditions.

- Add Tseng 2026 piperacillin ([doi:10.2147/DDDT.S602835](https://doi.org/10.2147/DDDT.S602835)) -- hospitalized adults with low body weight (BMI at or below 18.5 kg/m^2).

- Add Ugolkov 2026 B cell immune response QSP model
  ([doi:10.3389/fimmu.2026.1745710](https://doi.org/10.3389/fimmu.2026.1745710))
  -- healthy mice, aggregated across 21 published studies.

- Add Ding 2026 vancomycin ([doi:10.1007/s11096-025-02077-w](https://doi.org/10.1007/s11096-025-02077-w)) -- adults undergoing cardiac surgery.

- Add Li 2026 alvespimycin (17-DMAG) NSGA-II Pareto-front models at NEP = 5, 9,
  11 and 16 ([doi:10.1007/s10928-026-10036-9](https://doi.org/10.1007/s10928-026-10036-9))
  -- adults with advanced solid tumors.

- Add Yang 2026 encorafenib ([doi:10.1007/s40262-025-01608-y](https://doi.org/10.1007/s40262-025-01608-y)) -- adults with BRAF V600-mutant solid tumors and healthy participants.

- Add Gaffney 2026 niraparib ([doi:10.1002/jcph.70210](https://doi.org/10.1002/jcph.70210)) -- adults with advanced solid tumours or ovarian cancer.

- Add Feng 2026 imipenem PBPK ([doi:10.3389/fcimb.2026.1798911](https://doi.org/10.3389/fcimb.2026.1798911)) -- healthy adults, adults with renal impairment, and children aged 3-18 years with normal or impaired renal function.

- Add Blackman 2026 high-dose methotrexate ([doi:10.1007/s00228-026-04080-0](https://doi.org/10.1007/s00228-026-04080-0)) -- adults with leukemia, lymphoma, or sarcoma.

- Add Yang 2026 lurasidone ([doi:10.3389/fphar.2026.1810528](https://doi.org/10.3389/fphar.2026.1810528)) -- Chinese psychiatric inpatients aged 13-70 years.

- Add Chung 2026 vancomycin ([doi:10.3390/children13050649](https://doi.org/10.3390/children13050649)) -- preterm and term neonates in the neonatal intensive care unit.

- Add Ren 2026 rivaroxaban PK/PD ([doi:10.1016/j.rpth.2026.106618](https://doi.org/10.1016/j.rpth.2026.106618)) -- Chinese adults treated for pulmonary embolism.

- Add Yin 2026 ceftazidime ([doi:10.1128/aac.01810-25](https://doi.org/10.1128/aac.01810-25)) — Chinese neonates sampled by quantitative dried blood spot.

- Add Sulaiman 2026 piperacillin/tazobactam ([doi:10.1093/jac/dkag199](https://doi.org/10.1093/jac/dkag199)) -- critically ill adults with sepsis or septic shock in Malaysian intensive care units.

- Add Ma 2026 colistin sulfate ([doi:10.2147/DDDT.S600942](https://doi.org/10.2147/DDDT.S600942)) -- critically ill adults with carbapenem-resistant infections.

- Add Simpson 2006 artesunate ([doi:10.1371/journal.pmed.0030444](https://doi.org/10.1371/journal.pmed.0030444)) -- adults and children with moderately severe falciparum malaria.

- Add van Steeg 2007 S(-)-atenolol and isoprenaline PK-PD models ([doi:10.1038/sj.bjp.0707234](https://doi.org/10.1038/sj.bjp.0707234)) -- conscious male Wistar-Kyoto rats, with and without isoprenaline-induced tachycardia.

- Add Groenendaal 2007 morphine ([doi:10.1038/sj.bjp.0707258](https://doi.org/10.1038/sj.bjp.0707258)) -- male Wistar rats, EEG effect with P-glycoprotein interaction.

- Add Groenendaal 2007 morphine brain distribution ([doi:10.1038/sj.bjp.0707257](https://doi.org/10.1038/sj.bjp.0707257)) -- male Wistar rats.

- Add Geldof 2008 fluvoxamine SERT occupancy PK-PD models ([doi:10.1038/bjp.2008.179](https://doi.org/10.1038/bjp.2008.179)) -- male Wistar rats, ex vivo 5-HT transporter occupancy in frontal cortex.

- Add Rohatagi 2008 apricoxib ([doi:10.1111/j.1365-2125.2008.03175.x](https://doi.org/10.1111/j.1365-2125.2008.03175.x)) -- adults with acute postoperative dental pain.

- Add Hyland 2008 maraviroc in vitro CYP3A4 N-dealkylation kinetics ([doi:10.1111/j.1365-2125.2008.03198.x](https://doi.org/10.1111/j.1365-2125.2008.03198.x)) -- pooled human liver microsomes (60 donors) and recombinant CYP3A4 Supersomes.

- Add Muller 2009 amoxicillin ([doi:10.1128/AAC.00119-08](https://doi.org/10.1128/AAC.00119-08)) -- women in labour, the venous umbilical cord and their neonates.

- Add Hyland 2009 midazolam N-glucuronidation ([doi:10.1111/j.1365-2125.2009.03386.x](https://doi.org/10.1111/j.1365-2125.2009.03386.x)) -- in vitro, pooled human liver microsomes and recombinant UGT1A4.

- Add Mandema 2005 gemcabene, statins and ezetimibe LDL-C dose-response MBMA ([doi:10.1208/aapsj070352](https://doi.org/10.1208/aapsj070352)) -- adults with hypercholesterolemia, plus healthy-volunteer, obese and low-HDL-C cohorts.

- Add Benson 2010 reboxetine hNET target-binding kinetics
  ([doi:10.1111/j.1476-5381.2010.00719.x](https://doi.org/10.1111/j.1476-5381.2010.00719.x))
  -- in vitro, HEK-293 membranes expressing human noradrenaline transporter.

- Add Feng 2011 escitalopram ([doi:10.2147/NDT.S15921](https://doi.org/10.2147/NDT.S15921)) -- adults with chronic psychiatric disorders in a MEMS-monitored depression cohort.

- Add Marier 2011 cenicriviroc (TBR-652) HIV-1 RNA and MCP-1 exposure-response models ([doi:10.1128/AAC.00713-10](https://doi.org/10.1128/AAC.00713-10)) -- HIV-1-infected, antiretroviral treatment-experienced, CCR5-antagonist-naive adults.

- Add Roiko 2012 gamma-hydroxybutyric acid (GHB) blood-brain-barrier uptake
  ([doi:10.1124/dmd.111.041749](https://doi.org/10.1124/dmd.111.041749)) --
  in vitro, rat (RBE4) and human (hCMEC/D3) brain capillary endothelial cells.

- Add Morse 2012 gamma-hydroxybutyrate red-cell transport
  ([doi:10.1124/dmd.111.041285](https://doi.org/10.1124/dmd.111.041285)) --
  in vitro rat Sprague-Dawley erythrocytes.

- Add Fang 2013 pramlintide ([doi:10.1208/s12248-012-9409-7](https://doi.org/10.1208/s12248-012-9409-7)) -- adults with type 1 diabetes mellitus.

- Add Leeds 2013 tecovirimat (ST-246) ([doi:10.1128/AAC.00959-12](https://doi.org/10.1128/AAC.00959-12)) -- two models: uninfected and monkeypox-virus-infected cynomolgus monkeys, and healthy adult volunteers.

- Add Bruno 2012 capecitabine + docetaxel tumor growth inhibition, overall
  survival and progression-free survival models
  ([doi:10.1038/psp.2012.20](https://doi.org/10.1038/psp.2012.20)) -- adults
  with second-line locally advanced or metastatic breast cancer.

- Add Reif 2013 ethinylestradiol and drospirenone ([doi:10.1136/jfprhc-2012-100397](https://doi.org/10.1136/jfprhc-2012-100397)) -- healthy young women using an extended-cycle combined oral contraceptive.

- Add Rayner 2013 oseltamivir exposure-response models ([doi:10.1128/AAC.02440-12](https://doi.org/10.1128/AAC.02440-12)) — healthy adults experimentally inoculated with influenza A or B in two phase 2 studies.

- Add Snelder 2013 drug-independent cardiovascular systems (MAP / CO / TPR) model ([doi:10.1111/bph.12190](https://doi.org/10.1111/bph.12190)) -- conscious male spontaneously hypertensive rats.

- Add Dodds 2013 psoriasis biologics dose-response meta-analysis ([doi:10.1038/psp.2013.32](https://doi.org/10.1038/psp.2013.32)) -- adults with moderate-to-severe plaque psoriasis.

- Add Jones 2013 generic perfusion-limited whole-body PBPK template
  ([doi:10.1038/psp.2013.41](https://doi.org/10.1038/psp.2013.41)) -- no drug;
  70 kg human reference physiology.

- Add Snelder 2014 extended drug-independent cardiovascular systems (MAP / CO / HR / SV / TPR) model ([doi:10.1111/bph.12824](https://doi.org/10.1111/bph.12824)) -- conscious spontaneously hypertensive and normotensive Wistar-Kyoto rats.

- Add Riggs 2014 empagliflozin exposure-response ([doi:10.1111/bcp.12453](https://doi.org/10.1111/bcp.12453)) -- adults with type 2 diabetes.

- Add Chetty 2014 efalizumab CD11a target engagement and PASI response
  ([doi:10.3389/fimmu.2014.00670](https://doi.org/10.3389/fimmu.2014.00670)) --
  adults with moderate-to-severe plaque psoriasis.

- Add Cortez 2015 long-acting injectable nevirapine ([doi:10.1128/AAC.03906-14](https://doi.org/10.1128/AAC.03906-14)) -- Sprague-Dawley rats for the two formulation fits, and breastfeeding infants for the HIV-prophylaxis projections.

- Add Chandorkar 2015 ceftolozane and tazobactam ([doi:10.1002/jcph.395](https://doi.org/10.1002/jcph.395)) -- healthy adults, adults with mild to severe renal impairment, and patients with complicated urinary tract or intra-abdominal infections.

- Add Almond 2016 CYP3A4 induction by rifampicin, carbamazepine, phenobarbital, phenytoin, efavirenz and nifedipine ([doi:10.1124/dmd.115.066845](https://doi.org/10.1124/dmd.115.066845)) — cryopreserved human hepatocytes from four donors, plus healthy white adults for the in-vivo rifampicin reference.

- Add Guo 2016 taurocholate ([doi:10.1124/jpet.116.231928](https://doi.org/10.1124/jpet.116.231928)) -- sandwich-cultured human hepatocytes from three donors.

- Add Edwards 2016 obeticholic acid ([doi:10.1111/cts.12421](https://doi.org/10.1111/cts.12421)) -- healthy adults and subjects with Child-Pugh A/B/C cirrhosis.

- Add Gibbs 2017 evolocumab ([doi:10.1002/jcph.840](https://doi.org/10.1002/jcph.840)) -- healthy subjects and statin-treated patients with hypercholesterolemia.

- Add Kamal 2017 oseltamivir influenza-pandemic SEIR model ([doi:10.1111/bcp.13229](https://doi.org/10.1111/bcp.13229)) -- a simulated population of 100 000 individuals susceptible to pandemic influenza.

- Add Schulthess 2017 frequency-domain response analysis PD archetypes
  ([doi:10.1002/psp4.12266](https://doi.org/10.1002/psp4.12266)) -- 15
  theoretical models (no data fitted, no subjects).

- Add Nava 2018 busulfan ([doi:10.1111/bcp.13566](https://doi.org/10.1111/bcp.13566)) -- children and adolescents undergoing haematopoietic stem cell transplantation.

- Add Vujkovic 2018 efavirenz ([doi:10.1038/s41397-018-0028-2](https://doi.org/10.1038/s41397-018-0028-2)) -- black African adults with HIV-1 initiating efavirenz-based therapy in Botswana.

- Add Bartko 2018 sutimlimab classical-complement-pathway inhibition model ([doi:10.1002/cpt.1111](https://doi.org/10.1002/cpt.1111)) -- healthy volunteers in a first-in-human single- and multiple-ascending-dose trial.

- Add Chotsiri 2019 piperaquine ([doi:10.1038/s41467-019-08297-9](https://doi.org/10.1038/s41467-019-08297-9)) — children aged 2-58 months receiving seasonal malaria chemoprevention in Burkina Faso.

- Add Troy 2020 arylsulfatase A / TAK-611 ([doi:10.1002/cpt.1752](https://doi.org/10.1002/cpt.1752)) -- children with metachromatic leukodystrophy receiving intrathecal enzyme replacement.

- Add Kang 2020 adalimumab ([doi:10.1111/bcp.14330](https://doi.org/10.1111/bcp.14330))
  as three models (`_phase1`, `_phase3_base`, `_phase3_extension`) -- healthy
  male volunteers and adults with active rheumatoid arthritis.

- Add Keunecke 2020 regorafenib ([doi:10.1111/bcp.14334](https://doi.org/10.1111/bcp.14334)) -- adults with advanced solid tumours, as a phase 1 and a phase 3 model.

- Add Downes 2022 tobramycin ([doi:10.1128/aac.02377-21](https://doi.org/10.1128/aac.02377-21)) -- hospitalized children under 5 years of age with cystic fibrosis.

- Add Pohl 2022 linzagolix population PK and PK-oestradiol models ([doi:10.1111/bcp.15171](https://doi.org/10.1111/bcp.15171)) — healthy women and women with endometriosis.

- Add Shchelokov 2023 nivolumab PD-1 receptor-occupancy QSP model ([doi:10.1080/19420862.2022.2156317](https://doi.org/10.1080/19420862.2022.2156317)) -- virtual patients with advanced solid tumors.

- Add Zhang 2023 aztreonam, amoxicillin and clavulanic acid ([doi:10.3390/pharmaceutics15010251](https://doi.org/10.3390/pharmaceutics15010251)) -- simulated adults with creatinine clearance 10-150 mL/min/1.73 m2.

- Add Li 2023 TQ-B3203 ([doi:10.3389/fphar.2023.1102244](https://doi.org/10.3389/fphar.2023.1102244)) -- Chinese adults with advanced solid tumors.

- Enrich `Lahu_2010_roflumilast` provenance from the FDA NDA 22-522 Clinical
  Pharmacology review, which reprints the sponsor's population PK study report
  114/2005 -- the same analysis Lahu published. All 40 shipped parent and
  N-oxide `ini()` values reproduce the review's "Final With Race" columns
  exactly; **no parameter changed**. The review supplies provenance the journal
  article omits: the phase I demographics table now fills
  `population$sex_female_pct`, `population$weight_range`,
  `population$age_range` and `population$race_ethnicity` (previously `NULL` or
  placeholder text); alcohol use and several other tested-and-dropped
  covariates are recorded in `covariateData` notes; and two vignette errata are
  downgraded from open questions to confirmed readings -- the female-reference
  intercept (the FDA reviewer states it three times, so the paper's prose
  rather than its equations is wrong) and the genuine absence of inter-
  individual variability on N-oxide `tlag`. The withheld tPDE4i and
  adverse-event logistic-regression layer is *not* recoverable from the review;
  that gap is now documented as closed to further search.

- Add van Schaick 2016 prucalopride ([doi:10.1002/prp2.236](https://doi.org/10.1002/prp2.236)) -- children aged 6 months to 18 years with functional constipation.

- Add Feng 2019 ipilimumab mixture tumor growth dynamics model ([doi:10.1002/psp4.12454](https://doi.org/10.1002/psp4.12454)) -- adults with advanced melanoma.

- Add Chatterjee 2016 pembrolizumab tumor-size exposure-response model ([doi:10.1093/annonc/mdw174](https://doi.org/10.1093/annonc/mdw174)) -- adults with advanced non-small-cell lung cancer (KEYNOTE-001).

- Add Chandralayam Ayyappa Menon 2026 ISB 2001 and teclistamab QSP models ([doi:10.1002/cpt.70319](https://doi.org/10.1002/cpt.70319)) -- adults with relapsed and/or refractory multiple myeloma.

- Add Chanu 2010 methoxy polyethylene glycol-epoetin beta (C.E.R.A.) ([doi:10.1177/0091270009343931](https://doi.org/10.1177/0091270009343931)) -- ESA-naive and ESA-treated adults with chronic kidney disease on dialysis.

- Add Claret 2014 motesanib tumor size and overall survival ([doi:10.1038/clpt.2014.11](https://doi.org/10.1038/clpt.2014.11)) -- adults with advanced nonsquamous non-small cell lung cancer in the MONET1 study.

- Add Claret 2014 week-8 disease control rate OS and PFS models ([doi:10.1002/jcph.191](https://doi.org/10.1002/jcph.191)) -- Western and Chinese adults with first-line non-small cell lung cancer.

- Add Codaccioni 2024 hepatic CYP3A4 ontogeny, Salem vs modified Upreti ([doi:10.1002/jcph.2452](https://doi.org/10.1002/jcph.2452)) -- paediatric system-parameter model, birth to adulthood.

- Add Darwish 2012 armodafinil and modafinil ([doi:10.1177/0091270011417825](https://doi.org/10.1177/0091270011417825)) -- patients with excessive sleepiness associated with shift work disorder.

- Add Gosselin 2015 motesanib and its active metabolite M4 ([doi:10.1002/cpdd.196](https://doi.org/10.1002/cpdd.196)) -- patients with advanced solid tumors.

- Add Hopkins 2024 nonracemic amisulpride (SEP-4199) ([doi:10.1002/cpt.3311](https://doi.org/10.1002/cpt.3311)) -- healthy adult volunteers.

- Add Idkaidek 2011 ibuprofen at normal gravity and under simulated microgravity ([doi:10.1177/0091270010388652](https://doi.org/10.1177/0091270010388652)) -- healthy adult men in an antiorthostatic bed-rest crossover.

- Add Kamal 2014 oseltamivir ([doi:10.1038/clpt.2014.120](https://doi.org/10.1038/clpt.2014.120)) -- infants aged 2 weeks through 11 months with influenza.

- Add Kastrissios 2006 apricoxib (CS-706) ([doi:10.1177/0091270006287122](https://doi.org/10.1177/0091270006287122)) -- healthy adult volunteers.

- Add Kastrissios 2012 managlinat dialanetil (CS-917) ([doi:10.1177/0091270010396373](https://doi.org/10.1177/0091270010396373)) -- adults with type 2 diabetes mellitus.

- Add Knebel 2008 epoetin delta ([doi:10.1177/0091270008318218](https://doi.org/10.1177/0091270008318218)) -- pediatric patients with chronic kidney disease.

- Add Knebel 2011 pantoprazole ([doi:10.1177/0091270010366146](https://doi.org/10.1177/0091270010366146)) -- pediatric patients from birth to 16 years with gastroesophageal reflux disease.

- Add Knebel 2012 istradefylline ([doi:10.1177/0091270011420566](https://doi.org/10.1177/0091270011420566)) -- adults with Parkinson disease and levodopa-related motor complications.

- Add Knebel 2013 atorvastatin and o-hydroxyatorvastatin ([doi:10.1002/jcph.66](https://doi.org/10.1002/jcph.66)) -- children and adolescents with heterozygous familial hypercholesterolemia.

- Add Machavaram 2013 interleukin-6 / CYP3A4 suppression ([doi:10.1038/clpt.2013.79](https://doi.org/10.1038/clpt.2013.79)) -- adults with rheumatoid arthritis, bone marrow transplant or surgical trauma.

- Add Mandema 2011 anticoagulants dose-response meta-analysis ([doi:10.1038/clpt.2011.232](https://doi.org/10.1038/clpt.2011.232)) -- adults undergoing total hip or knee replacement surgery.

- Add Mandema 2011 biologic DMARDs ([doi:10.1038/clpt.2011.256](https://doi.org/10.1038/clpt.2011.256)) -- adults with rheumatoid arthritis.

- Add Mao 2012 vernakalant ([doi:10.1177/0091270011408425](https://doi.org/10.1177/0091270011408425)) -- adults with atrial fibrillation or atrial flutter and healthy volunteers.

- Add Marier 2010 teduglutide ([doi:10.1177/0091270009342252](https://doi.org/10.1177/0091270009342252)) -- healthy adults and patients with short bowel syndrome, Crohn's disease, or moderate hepatic impairment.

- Add Marier 2014 dutogliptin ([doi:10.1002/cpdd.87](https://doi.org/10.1002/cpdd.87)) -- healthy subjects and patients with type 2 diabetes mellitus.

- Add Mascarenhas 2015 pentadecanoic and triheptadecanoic acid ([doi:10.1002/jcph.484](https://doi.org/10.1002/jcph.484)) -- children and adults with cystic fibrosis and pancreatic insufficiency, plus healthy comparison subjects.

- Add Mondick 2018 empagliflozin ([doi:10.1002/jcph.1051](https://doi.org/10.1002/jcph.1051)) -- adults with type 1 diabetes, plus the type 2 diabetes parameter set reported alongside it.

- Add Mouksassi 2015 thrombomodulin alfa ([doi:10.1002/cpdd.163](https://doi.org/10.1002/cpdd.163)) -- healthy adults and adults with sepsis and suspected disseminated intravascular coagulation, across normal to severe renal impairment.

- Add Nielsen 2015 vigabatrin ([doi:10.1002/jcph.378](https://doi.org/10.1002/jcph.378)) - adults and children with refractory complex partial seizures.

- Add Niu 2017 veliparib and its M8 metabolite ([doi:10.1002/jcph.892](https://doi.org/10.1002/jcph.892)) -- patients with BRCA 1/2-mutated cancer or PARP-sensitive tumor types.

- Add Panday 2025 SERT inhibitors tremor MBMA ([doi:10.1002/cpt.3696](https://doi.org/10.1002/cpt.3696)) -- 29,677 adults across 33 treatment arms and 20 serotonin reuptake transporter inhibitors.

- Add Reinecke 2018 levonorgestrel contraceptives, as eight models covering the intrauterine, oral and subdermal routes ([doi:10.1002/jcph.1288](https://doi.org/10.1002/jcph.1288)) -- healthy premenopausal women.

- Add Riccobene 2017 ceftaroline ([doi:10.1002/jcph.809](https://doi.org/10.1002/jcph.809)) -- children from birth to under 18 years pooled with adults.

- Add Riggs 2012 albinterferon alfa-2b ([doi:10.1177/0091270011399576](https://doi.org/10.1177/0091270011399576)) -- adults with chronic hepatitis C virus infection.

- Add Rini 2013 axitinib ([doi:10.1002/jcph.73](https://doi.org/10.1002/jcph.73)) -- pooled healthy volunteers and patients with metastatic renal cell carcinoma or other solid tumours.

- Add Siccardi 2012 efavirenz ([doi:10.1038/clpt.2012.61](https://doi.org/10.1038/clpt.2012.61)) -- European HIV-1-positive adults and healthy volunteers.

- Add Song 2013 olmesartan / amlodipine / hydrochlorothiazide (CS-8635) ([doi:10.1002/cpdd.17](https://doi.org/10.1002/cpdd.17)) -- adults with hypertension.

- Add Stringer 2013 sipoglitazar ([doi:10.1177/0091270012447121](https://doi.org/10.1177/0091270012447121)) -- healthy volunteers and adults with type 2 diabetes; the UGT2B15 genotype model and the parallel latent-subpopulation mixture model.

- Add Stringer 2014 sipoglitazar ([doi:10.1002/jcph.227](https://doi.org/10.1002/jcph.227)) -- drug-naive adults with type 2 diabetes mellitus.

- Add Zhu 2013 ganitumab ([doi:10.1002/cpdd.48](https://doi.org/10.1002/cpdd.48)) -- adults with metastatic pancreatic cancer or other advanced solid cancers.

- Add Butragueno-Laiseca 2022 piperacillin ([doi:10.1016/j.cmi.2022.03.031](https://doi.org/10.1016/j.cmi.2022.03.031)) -- critically ill children with and without continuous kidney replacement therapy.

- Add Butragueno-Laiseca 2024 meropenem ([doi:10.1128/aac.01729-23](https://doi.org/10.1128/aac.01729-23)) -- critically ill children in a paediatric ICU, with and without continuous kidney replacement therapy.

- Add Elmokadem 2023 mavoglurant whole-body PBPK ([doi:10.1002/psp4.12926](https://doi.org/10.1002/psp4.12926)) -- healthy adults receiving a single intravenous infusion.

- Add Yin 2023 soticlestat ([doi:10.1111/cts.13517](https://doi.org/10.1111/cts.13517)) -- healthy adults.

- Add Fukae 2024 valemetostat exposure-response ([doi:10.1002/psp4.13203](https://doi.org/10.1002/psp4.13203)) -- adults with relapsed/refractory adult T-cell leukemia/lymphoma and other non-Hodgkin lymphomas.

- Add Lin 2024 adverse-event grade proportional-odds Markov model ([doi:10.3389/fphar.2024.1487062](https://doi.org/10.3389/fphar.2024.1487062)) -- simulated subjects on a de-identified oral compound.

- Add Siemers 2025 sabirnetug ([doi:10.1016/j.tjpad.2024.100005](https://doi.org/10.1016/j.tjpad.2024.100005)) -- adults with mild cognitive impairment or mild dementia due to Alzheimer's disease.

- Add Cao 2025 busulfan ([doi:10.3389/fphar.2025.1632588](https://doi.org/10.3389/fphar.2025.1632588)) -- Chinese children undergoing allogeneic hematopoietic cell transplantation.

- Add Assmus 2025 benznidazole ([doi:10.1371/journal.pntd.0013522](https://doi.org/10.1371/journal.pntd.0013522)) -- adults with chronic indeterminate Chagas disease in the BENDITA trial; population PK plus a qPCR-positivity exposure-response model.

- Add Zhang 2025 anti-psoriatic "Drug A" ([doi:10.1002/psp4.70090](https://doi.org/10.1002/psp4.70090)) -- adults with plaque psoriasis in a simulated multi-regional phase 2 trial.

- Add Muhamad 2025 cholecalciferol PBPK ([doi:10.3390/nu17193028](https://doi.org/10.3390/nu17193028)) -- healthy Cape Town schoolchildren aged 6-11 years.

- Add Park 2025 efineptakin alfa (rhIL-7-hyFc) ([doi:10.2147/DDDT.S564085](https://doi.org/10.2147/DDDT.S564085)) -- adults with locally advanced or metastatic solid tumours.

- Add Zhou 2025 fruquintinib concentration-QTc ([doi:10.1002/jcph.70051](https://doi.org/10.1002/jcph.70051)) -- adults with previously treated metastatic colorectal cancer in the phase 3 FRESCO-2 ECG substudy.

- Add Yang 2025 matrine ([doi:10.3389/fvets.2025.1620161](https://doi.org/10.3389/fvets.2025.1620161)) -- pigs, intestinal-lumen exposure after oral dosing.

- Add Dohmann 2025 piperacillin ([doi:10.1002/bcp.70153](https://doi.org/10.1002/bcp.70153)) -- adults with renal impairment or on intermittent haemodialysis in a non-intensive-care ward setting.

- Add Karakitsios 2025 bedaquiline lung PBPK models ([doi:10.1002/bcp.70163](https://doi.org/10.1002/bcp.70163)) -- mouse, rat (Sprague-Dawley), beagle dog, and adults with pulmonary drug-resistant tuberculosis (ratifies the new numbered-chain compartment canonical `caseum<n>`, the concentric rings of the necrotic caseous core of a tuberculosis granuloma, numbered outward-to-inward from the outer caseum edge; the mouse and human models previously carried these states under the `paper_specific_compartments` escape hatch).

- Add Sun 2025 paclitaxel ([doi:10.1111/cts.70404](https://doi.org/10.1111/cts.70404)) -- adult women with breast cancer receiving weekly paclitaxel.

- Add Hu 2025 tepotinib ([doi:10.3389/fphar.2025.1685468](https://doi.org/10.3389/fphar.2025.1685468)) -- adults with MET-aberrant non-small cell lung cancer; the human PK submodule of a MET-pathway QSP model. Parameters are transcribed verbatim, but the source's own printed values do not reproduce its own published PK figure; read the vignette before using it to predict exposure.

- Add Perlstein 2026 extended-release injectable olanzapine TV-44749 and oral olanzapine ([doi:10.1002/jcph.70144](https://doi.org/10.1002/jcph.70144)) -- healthy adults and adults with schizophrenia or schizoaffective disorder.

- Add Echterhof 2026 Pseudomonas bacteriophage whole-body PBPK ([doi:10.1128/aac.01506-25](https://doi.org/10.1128/aac.01506-25)) -- three models: CD-1 mice (the fitted species), plus rat and human interspecies projections.

- Add Kim 2026 midazolam ([doi:10.1097/ALN.0000000000005811](https://doi.org/10.1097/ALN.0000000000005811)) -- adults on venoarterial ECMO and after decannulation.

- Add Yu 2026 tenofovir ([doi:10.1007/s40262-025-01589-y](https://doi.org/10.1007/s40262-025-01589-y)) -- non-pregnant, pregnant and postpartum women taking tenofovir disoproxil fumarate or tenofovir alafenamide.

- Add Stitt 2026 tranexamic acid ([doi:10.1111/trf.70047](https://doi.org/10.1111/trf.70047)) -- adults with severe traumatic injury, allometrically scaled for extrapolation to children with trauma-related bleeding.

- Add Tsuchitani 2026 apixaban ([doi:10.1002/psp4.70163](https://doi.org/10.1002/psp4.70163)) - healthy adults; whole-body PBPK with enterohepatic circulation.

- Add Thoueille 2026 inhaled salmeterol and alpha-hydroxysalmeterol ([doi:10.1002/psp4.70187](https://doi.org/10.1002/psp4.70187)) -- healthy adults, chronic asthmatics and athletes/endurance-trained individuals.

- Add Hanan 2026 pegylated-interferon-alfa HBsAg-loss model-based meta-analysis ([doi:10.1002/psp4.70164](https://doi.org/10.1002/psp4.70164)) -- adults with chronic hepatitis B virus infection.

- Add Storgaard 2026 delta-9-tetrahydrocannabinol ([doi:10.1002/bcp.70284](https://doi.org/10.1002/bcp.70284)) -- acutely hospitalized older medical patients with poor appetite.

- Add Ward 2026 magnesium bone implant PBPK ([doi:10.1021/acsomega.5c06910](https://doi.org/10.1021/acsomega.5c06910)) -- average healthy 70 kg adults.

- Add Mukker 2026 tuvusertib ([doi:10.1002/cpt.70029](https://doi.org/10.1002/cpt.70029)) -- adults with advanced solid tumors, plus a mouse ARID1A-mutant xenograft tumor-growth model.

- Add Zhang 2026 ribociclib intracranial CDK4/6 occupancy QSP ([doi:10.1186/s12885-026-15561-x](https://doi.org/10.1186/s12885-026-15561-x)) -- patients with breast cancer and brain metastasis.

- Add Mahadevan 2026 polymyxin B + meropenem + fosfomycin mechanism-based time-kill models ([doi:10.1128/aac.00782-25](https://doi.org/10.1128/aac.00782-25)) -- six carbapenem-resistant Klebsiella pneumoniae clinical isolates in vitro; one file per isolate.

- Add Zurawska 2026 piperacillin ([doi:10.1128/aac.01760-25](https://doi.org/10.1128/aac.01760-25)) -- critically ill adults with hospital-acquired pneumonia, including those on continuous renal replacement therapy.

- Add Mukker 2026 tuvusertib concentration-QTcF, concentration-heart-rate and in vitro hERG models ([doi:10.1111/cts.70496](https://doi.org/10.1111/cts.70496)) -- patients with advanced solid tumors.

- Add Klein 2026 inhaled Staccato alprazolam ([doi:10.1111/epi.18643](https://doi.org/10.1111/epi.18643)) -- adults and adolescents with epilepsy.

- Add Sun 2026 tilmicosin ([doi:10.1021/acs.jafc.5c11368](https://doi.org/10.1021/acs.jafc.5c11368)) -- healthy swine, with hollow-fiber PK/PD against *Pasteurella multocida*.

- Add Tong 2026 vancomycin ([doi:10.1093/jacamr/dlag016](https://doi.org/10.1093/jacamr/dlag016)) -- four MIPD models (modified Goti, capped Thomson, Carreno, Hughes) in hospitalized adults.

- Add Yang 2026 APTM, a novel pleuromutilin derivative ([doi:10.1016/j.psj.2026.106560](https://doi.org/10.1016/j.psj.2026.106560)) -- Mycoplasma gallisepticum in vitro and in specific-pathogen-free chickens.

- Add Bardhi 2026 ampicillin ([doi:10.1093/jvimsj/aalag021](https://doi.org/10.1093/jvimsj/aalag021)) -- hospitalized neonatal foals.

- Add Rolsma 2026 cefepime ([doi:10.1093/ofid/ofag069](https://doi.org/10.1093/ofid/ofag069)) -- critically ill children in intensive care.

- Add Hu 2026 utreloxastat ([doi:10.1002/jcph.70171](https://doi.org/10.1002/jcph.70171)) -- healthy adults in a first-in-human study.

- Add Subrtova 2026 daptomycin ([doi:10.1128/aac.01532-25](https://doi.org/10.1128/aac.01532-25)) -- adults with serious Gram-positive infections.

- Add Taylor 2026 methotrexate ([doi:10.1007/s40262-026-01618-4](https://doi.org/10.1007/s40262-026-01618-4)) -- adults with lymphoma receiving high-dose methotrexate.

- Add Lau 2026 paracetamol ([doi:10.1002/psp4.70168](https://doi.org/10.1002/psp4.70168)) - adults with and without obesity, oral and intravenous.

- Add Retout 2026 baloxavir ([doi:10.1002/cpt.70204](https://doi.org/10.1002/cpt.70204)) -- influenza patients aged 1 year and older, pooled pediatric, adolescent and adult.

- Add Zhao 2026 venetoclax ([doi:10.2147/DDDT.S583847](https://doi.org/10.2147/DDDT.S583847)) -- Chinese pediatric patients with hematological malignancy.

- Add Teixeira 2026 meloxicam ([doi:10.3390/ph19030412](https://doi.org/10.3390/ph19030412)) - female mixed-breed dogs undergoing ovariohysterectomy, free vs polymeric-nanocapsule oral formulations.

- Add Nemitz 2026 dapagliflozin whole-body PBPK/PD ([doi:10.3390/pharmaceutics18030287](https://doi.org/10.3390/pharmaceutics18030287)) -- healthy adults, type 1 and type 2 diabetes, and renal or hepatic impairment.

- Add Soria-Chacartegui 2026 tramadol ([doi:10.1007/s13318-026-00986-3](https://doi.org/10.1007/s13318-026-00986-3)) - European healthy volunteers genotyped for CYP2D6.

- Add Zhou 2026 tacrolimus ([doi:10.2147/DDDT.S576598](https://doi.org/10.2147/DDDT.S576598)) -- adult Chinese patients with nephrotic syndrome.

- Add Schoning 2026 artificial-patient popPK methodology models ([doi:10.1002/prp2.70241](https://doi.org/10.1002/prp2.70241)) - a hypothetical drug in 20 simulated adults, packaged as the author-defined ground truth plus the five WGAN-GP-augmented fitted arms.

- Add Yu 2026 lacosamide ([doi:10.1186/s40360-026-01114-2](https://doi.org/10.1186/s40360-026-01114-2)) -- adults with epilepsy.

- Add Henninger 2026 DNDI-6148 ([doi:10.1111/cts.70535](https://doi.org/10.1111/cts.70535)) -- two files: the murine target-site PK/PD fit in *Leishmania major*-infected BALB/c mice, and its allometrically scaled human translation for efficacious-dose prediction in cutaneous leishmaniasis.

- Add Lledo-Garcia 2022 rozanolixizumab ([doi:10.1002/psp4.12739](https://doi.org/10.1002/psp4.12739)) -- cynomolgus monkeys, and healthy adults in the first-in-human study.

- Rename `Hansson_2013_sunitinib_OS` to `Hansson_2013_sunitinib_svegfr3_os` and its
  vignette to match, resolving a case-only filename collision with the sibling
  `Hansson_2013_sunitinib_os` that broke `git clone` on case-insensitive
  filesystems (issue 492). The two are DIFFERENT models from companion papers
  published back-to-back (e84 doi:10.1038/psp.2013.61, sVEGFR-3-driven; e85
  doi:10.1038/psp.2013.62, adverse-effect-driven); both now cross-reference each
  other explicitly.
- Add Siebinga 2023 [68Ga]Ga- and [177Lu]Lu-HA-DOTATATE ([doi:10.1186/s40658-023-00565-4](https://doi.org/10.1186/s40658-023-00565-4)) - adults with neuroendocrine tumors.

- Add Gracia 2025 51Cr-EDTA and 99mTc-DTPA GFR tracers ([doi:10.1007/s00467-025-06828-9](https://doi.org/10.1007/s00467-025-06828-9)) -- oncopediatric children receiving cisplatin and/or ifosfamide.

- Add Koh 2025 enteric-coated aspirin ([doi:10.2147/DDDT.S533428](https://doi.org/10.2147/DDDT.S533428)) -- healthy Korean adults.

- Add Royston 2025 letermovir ([doi:10.1128/aac.00697-25](https://doi.org/10.1128/aac.00697-25)) -- adult allogeneic hematopoietic cell transplant recipients.

- Add Matcha 2025 amikacin ([doi:10.1038/s41390-025-04044-7](https://doi.org/10.1038/s41390-025-04044-7)) -- term neonates.

- Add Olsson Gisleskog 2025 ibrutinib ([doi:10.1002/psp4.70061](https://doi.org/10.1002/psp4.70061)) -- patients 65 years and older with previously untreated mantle cell lymphoma in the phase 3 SHINE study.

- Add Ujihira 2025 glycochenodeoxycholic acid 3-O-sulfate, an endogenous OATP1B3 / OAT3 biomarker ([doi:10.1002/cpt.70023](https://doi.org/10.1002/cpt.70023)) -- healthy adults.

- Add Chigutsa 2025 tirzepatide ([doi:10.1002/cpt.3750](https://doi.org/10.1002/cpt.3750)) -- adults with obesity or overweight without type 2 diabetes (SURMOUNT-1).

- Add Ibrahim 2025 ibrutinib ([doi:10.1002/psp4.70124](https://doi.org/10.1002/psp4.70124)) -- adults with treatment-naive or relapsed/refractory chronic lymphocytic leukemia.

- Add Boivin-Champeaux 2026 acute hepatitis B virus QSP model ([doi:10.1002/psp4.70172](https://doi.org/10.1002/psp4.70172)) -- adults with acute HBV infection pooled from eight published cohorts.

- Add Proctor 2026 durvalumab ([doi:10.1002/psp4.70185](https://doi.org/10.1002/psp4.70185)) -- adults with advanced cancer; a whole-body two-pore PBPK model of cachexia-driven time-dependent clearance plus the empirical popPK comparator.

- Add Almquist 2015 Mig1 nuclear-localization models ([doi:10.1371/journal.pone.0124050](https://doi.org/10.1371/journal.pone.0124050)) -- in vitro single yeast cells (Saccharomyces cerevisiae); four per-experiment fits.

- Add Yau 2023 diazepam and midazolam simplified PBPK models ([doi:10.1002/psp4.12915](https://doi.org/10.1002/psp4.12915)) -- rat, cynomolgus monkey and the human translation.

* Add Ketharanathan 2023 pentobarbital ([doi:10.1007/s40262-023-01249-z](https://doi.org/10.1007/s40262-023-01249-z)) -- critically ill children in the paediatric intensive care unit treated for refractory status epilepticus or severe traumatic brain injury.

- Add Xu 2023 busulfan ([doi:10.1002/psp4.13004](https://doi.org/10.1002/psp4.13004)) -- Chinese adults and children undergoing allogeneic hematopoietic stem cell transplantation.

- Add Bjornsson 2023 buprenorphine / CAM2038 depot ([doi:10.1007/s40262-023-01288-6](https://doi.org/10.1007/s40262-023-01288-6)) -- healthy adults and adults with opioid use disorder.

- Add Zhu 2024 pyrotinib ([doi:10.3389/fphar.2024.1432944](https://doi.org/10.3389/fphar.2024.1432944)) - Chinese patients with HER2-positive breast cancer.

- Add Kado 2023 benzathine benzylpenicillin G ([doi:10.1128/aac.00962-23](https://doi.org/10.1128/aac.00962-23)) - healthy adults receiving high-dose subcutaneous infusions.

- Add Roy 2023 erythropoiesis QSP model ([doi:10.3389/fphar.2023.1274490](https://doi.org/10.3389/fphar.2023.1274490)) - virtual patients with anemia due to chronic kidney disease.

- Add Yin 2024 soticlestat ([doi:10.1111/cts.13722](https://doi.org/10.1111/cts.13722)) -- healthy volunteers and patients with developmental and epileptic encephalopathies.

- Add Li 2024 ORIN1001 ([doi:10.3389/fphar.2024.1322557](https://doi.org/10.3389/fphar.2024.1322557)) - Chinese patients with advanced solid tumors.

- Add Winchell 2024 posaconazole ([doi:10.1128/aac.01197-23](https://doi.org/10.1128/aac.01197-23)) - pediatric patients aged 2 to 17 years with neutropenia.

- Register three new canonical parameter names in
  `inst/references/parameter-names.md`: `boxcox_<param>` for a Box-Cox shape
  parameter applied to an interindividual random effect (distinct from
  `lambda_<output>`, which transforms the *residual* error), `lkd_direct` /
  `lkd_delay` for linear drug-effect coefficients on a turnover loss rate
  (distinct from `kd`, a mechanistic dissociation rate), and `e_wt_k` for a
  single allometric exponent a paper applies to every rate constant at once.
  The `sf<n>` note that had pointed at the deprecated `allo_<param>` spelling
  now points at `e_wt_<param>`, which is what `checkModelConventions()`
  actually enforces.

- Add Verrest 2024 *Leishmania* blood parasite dynamics ([doi:10.1371/journal.pntd.0012078](https://doi.org/10.1371/journal.pntd.0012078)) - Eastern African visceral leishmaniasis patients treated with liposomal amphotericin B, sodium stibogluconate, miltefosine, or fexinidazole.

- Add Derippe 2024 BH3-mimetic apoptosis QSP, mouse venetoclax PK, and SU-DHL-4 xenograft growth models ([doi:10.1002/psp4.13158](https://doi.org/10.1002/psp4.13158)) -- SU-DHL-4 and KARPAS-422 lymphoma cell lines in vitro, with a mouse xenograft bridge.

- Add Kiriyama 2024 nifedipine and captopril ([doi:10.1002/prp2.1249](https://doi.org/10.1002/prp2.1249)) -- spontaneously hypertensive rats under urethane anaesthesia.

- Add Moein 2024 apitolisib ([doi:10.1007/s40268-024-00459-5](https://doi.org/10.1007/s40268-024-00459-5)) - 786-O xenograft mice and patients with advanced solid tumors.

- Add Alfosea-Cuadrado 2024 reserpine ([doi:10.3390/pharmaceutics16081101](https://doi.org/10.3390/pharmaceutics16081101)) - male Sprague-Dawley rats in the reserpine-induced myalgia model of fibromyalgia.

- Add Okada 2024 triazolam ([doi:10.1186/s40360-024-00777-z](https://doi.org/10.1186/s40360-024-00777-z)) - young and elderly adults.

- Add Henthorn 2024 delta-9-tetrahydrocannabinol ([doi:10.1097/FTD.0000000000001224](https://doi.org/10.1097/FTD.0000000000001224)) - adults inhaling commercial-market cannabis, occasional and daily users.

- Add Ohara 2014 S-warfarin PK/PD ([doi:10.1371/journal.pone.0105891](https://doi.org/10.1371/journal.pone.0105891)) -- Chinese adults in Taiwan starting warfarin induction therapy.

- Add van Os 2024 temocillin PK/PD ([doi:10.1093/jac/dkae243](https://doi.org/10.1093/jac/dkae243)) -- four in vitro *Escherichia coli* strains, with plasma pharmacokinetics from critically ill patients.

- Add Fan 2024 sirolimus ([doi:10.3389/fphar.2024.1457614](https://doi.org/10.3389/fphar.2024.1457614)) - Chinese children with vascular anomalies, aged 0.08-12 years.

- Add Khaowroongrueng 2024 sufentanil ([doi:10.1002/psp4.13205](https://doi.org/10.1002/psp4.13205)) - adult Korean patients undergoing cardiopulmonary bypass surgery.

- Add Yang 2024 naloxone auto-injector and its four mechanistic opioid-induced
  respiratory depression reversal models -- buprenorphine, morphine, fentanyl
  and carfentanil ([doi:10.1002/psp4.13215](https://doi.org/10.1002/psp4.13215))
  - healthy adults.

- Add Ivanova 2024 anti-alpha-synuclein immunotherapy ([doi:10.1002/psp4.13223](https://doi.org/10.1002/psp4.13223)) -- QSP model of alpha-synuclein pathology in Parkinson's-disease-like mice.

- Add Zhang 2024 F-53B ([doi:10.1021/acs.est.4c05405](https://doi.org/10.1021/acs.est.4c05405)) - gestational PBPK in pregnant mice, plus pre-pregnant and pregnant human extrapolations.

- Add Zhang 2024 tucatinib ([doi:10.1007/s40262-024-01412-0](https://doi.org/10.1007/s40262-024-01412-0)) - healthy participants and patients with HER2+ metastatic breast or colorectal cancer.

- Add Zuo 2024 apatinib ([doi:10.1186/s12885-024-13118-4](https://doi.org/10.1186/s12885-024-13118-4)) - Chinese adults with solid tumours.

- Add Duan 2024 linezolid ([doi:10.1128/aac.01148-24](https://doi.org/10.1128/aac.01148-24)) - Chinese premature neonates undergoing therapeutic drug monitoring.

- Add Saeheng 2024 Atractylodes lancea ([doi:10.1186/s12906-024-04618-8](https://doi.org/10.1186/s12906-024-04618-8)) - patients with advanced-stage intrahepatic cholangiocarcinoma.

- Add Bian 2024 lefamulin popPK plus epithelial lining fluid models ([doi:10.3389/fphar.2024.1456741](https://doi.org/10.3389/fphar.2024.1456741)) -- pooled Phase 1 healthy adults, Phase 2 skin-infection patients, and Phase 3 community-acquired bacterial pneumonia patients.

- Add Wang 2024 RNA lipid nanoparticles ([doi:10.1016/j.apsb.2024.06.011](https://doi.org/10.1016/j.apsb.2024.06.011)) - ionizable-lipid PBPK in Sprague-Dawley rats, CD-1 and C57BL/6 mice, and healthy volunteers given patisiran, plus a cellular trafficking model in HeLa cells.

- Add Abdullah-Koolmees 2024 voriconazole and posaconazole whole-body PBPK ([doi:10.1007/s13318-024-00916-1](https://doi.org/10.1007/s13318-024-00916-1)) -- a typical 73-kg adult, simulated as healthy, bacterially infected and ICU patient types with and without concomitant flucloxacillin.

- Add Iwama 2024 febuxostat ([doi:10.1002/prp2.70032](https://doi.org/10.1002/prp2.70032)) -- Japanese pediatric patients with hyperuricemia including gout and adults who were healthy or had renal dysfunction.

- Add Willemin 2024 interleukin-6 / CYP modulation after talquetamab ([doi:10.1007/s11523-024-01093-6](https://doi.org/10.1007/s11523-024-01093-6)) - adults with relapsed/refractory multiple myeloma and cytokine release syndrome.

- Add Chotsiri 2019 lumefantrine ([doi:10.1002/cpt.1531](https://doi.org/10.1002/cpt.1531)) -- children aged 6-59 months with uncomplicated malaria, with and without severe acute malnutrition.

- Add Xu 2024 linezolid ([doi:10.2147/DDDT.S474470](https://doi.org/10.2147/DDDT.S474470)) -- adults with hospital-acquired pneumonia and renal insufficiency.

- Add Li 2024 norepinephrine ([doi:10.1007/s40262-024-01430-y](https://doi.org/10.1007/s40262-024-01430-y)) -- healthy volunteers, awake and under propofol/remifentanil general anesthesia.

- Add Ashraf 2024 codeine ([doi:10.1007/s40262-024-01433-9](https://doi.org/10.1007/s40262-024-01433-9)) -- CYP2D6-genotyped adults undergoing ambulatory surgery.

- Add Choules 2024 enfortumab vedotin and brentuximab vedotin ([doi:10.1007/s10928-023-09877-5](https://doi.org/10.1007/s10928-023-09877-5)) -- adults with metastatic urothelial carcinoma or CD30-positive lymphoma; compartmental reductions of the published Simcyp ADC-module PBPK models, with the MMAE payload and its CYP3A4/P-gp drug interactions.

- Add Patidar 2024 platform monoclonal-antibody minimal-PBPK models ([doi:10.1007/s10928-023-09899-z](https://doi.org/10.1007/s10928-023-09899-z)) -- mouse (wild-type and FcRn knockout) and 70 kg adult human.

- Add Ojara 2024 lamivudine ([doi:10.1002/psp4.13274](https://doi.org/10.1002/psp4.13274)) -- breastfeeding Ugandan women living with HIV, with paired maternal plasma and breast-milk data.

- Add Sexton 2024 lanadelumab kallikrein-kinin system QSP model ([doi:10.1007/s10928-024-09919-6](https://doi.org/10.1007/s10928-024-09919-6)) -- adults with hereditary angioedema due to C1-inhibitor deficiency.

- Add Zhang 2024 valproic acid ([doi:10.3389/fphar.2024.1423411](https://doi.org/10.3389/fphar.2024.1423411)) -- Chinese children and adults with epilepsy or after neurosurgery.

- Add Tuey 2024 cholecalciferol ([doi:10.3390/ijms252212279](https://doi.org/10.3390/ijms252212279)) -- adults with chronic kidney disease and vitamin D insufficiency or deficiency.

- Add Vaddady 2024 quizartinib and its active metabolite AC886 ([doi:10.1111/cts.70074](https://doi.org/10.1111/cts.70074)) -- adults with newly diagnosed or relapsed/refractory FLT3-ITD-positive acute myeloid leukemia, healthy volunteers and subjects with hepatic impairment.

- Add Sadiq 2024 tozorakimab ([doi:10.1111/bcp.16195](https://doi.org/10.1111/bcp.16195)) -- healthy adults and patients with mild COPD.

- Add Zhang 2024 nedosiran ([doi:10.1111/bcp.16194](https://doi.org/10.1111/bcp.16194)) -- healthy adults and patients with primary hyperoxaluria type 1 or 2.

- Add Yang 2024 dabigatran ([doi:10.3389/fphar.2024.1454612](https://doi.org/10.3389/fphar.2024.1454612)) -- healthy Chinese adults.

- Add Na 2024 tovecimig / ABL001 / CTX-009 ([doi:10.1111/cas.16363](https://doi.org/10.1111/cas.16363)) -- adults with relapsed or refractory solid tumors.

- Add Mo 2024 shikimic acid ([doi:10.1021/acs.jafc.4c09250](https://doi.org/10.1021/acs.jafc.4c09250)) -- growing pigs (Landrace x Large White) dosed intravenously and intragastrically.

- Add Suzuki 2024 mycophenolic acid ([doi:10.1111/cts.70097](https://doi.org/10.1111/cts.70097)) -- adult kidney transplant recipients.

- Add Kurup 2024 DZIF-10c ([doi:10.1007/s10928-024-09947-2](https://doi.org/10.1007/s10928-024-09947-2)) -- pooled adults (SARS-CoV-2 infected and uninfected) and cynomolgus macaques.

- Add Khwarg 2024 proguanil ([doi:10.1111/cts.70103](https://doi.org/10.1111/cts.70103)) -- healthy Korean male adults genotyped for SLC22A1 (OCT1) 1022C>T.

- Add Yang 2024 zastaprazan ([doi:10.1002/psp4.13228](https://doi.org/10.1002/psp4.13228)) -- patients with erosive gastroesophageal reflux disease and healthy volunteers.

- Add Willmann 2024 elinzanetant ([doi:10.1002/psp4.13226](https://doi.org/10.1002/psp4.13226)) -- healthy volunteers and women with vasomotor symptoms associated with menopause.

- Add Crass 2024 pegcetacoplan population PK and hemoglobin / LDH PK/PD models ([doi:10.1007/s40268-024-00500-7](https://doi.org/10.1007/s40268-024-00500-7)) -- healthy adults and adults with paroxysmal nocturnal hemoglobinuria.

- Add Winter 2024 oxytetracycline ([doi:10.3389/fmicb.2024.1498219](https://doi.org/10.3389/fmicb.2024.1498219)) -- calves and adult cattle.

- Add Schouwenburg 2025 clavulanic acid ([doi:10.1002/cpt.3423](https://doi.org/10.1002/cpt.3423)) -- preterm and term neonates and infants up to 1 year of age.

- Add Melander 2025 cetirizine ([doi:10.1111/bcpt.14100](https://doi.org/10.1111/bcpt.14100)) -- breast milk of lactating women.

- Add Abouelhassan 2024 sulbactam ([doi:10.1093/jacamr/dlae203](https://doi.org/10.1093/jacamr/dlae203)) -- mouse (ICR/CD-1 neutropenic *Acinetobacter baumannii* pneumonia model) and healthy adults, both with an epithelial lining fluid compartment.

- Add Alqurain 2024 vancomycin ([doi:10.2147/DDDT.S496512](https://doi.org/10.2147/DDDT.S496512)) -- non-critical-care adults aged 40 years and older on medical wards in Saudi Arabia.

- Add Zhu 2024 SPT-07A (D-borneol) whole-body PBPK ([doi:10.3390/pharmaceutics16121596](https://doi.org/10.3390/pharmaceutics16121596)) -- three models: rats (Sprague-Dawley), beagle dogs, and healthy human adults.

- Add Panetta 2024 palbociclib ([doi:10.3390/pharmaceutics16121528](https://doi.org/10.3390/pharmaceutics16121528)) -- children and young adults with recurrent, progressive or refractory brain tumors.

- Add Liu 2024 Spatholobi Caulis constituents 3'-methoxydaidzein, 8-O-methylretusin, daidzin and isolariciresinol ([doi:10.3390/ph17121621](https://doi.org/10.3390/ph17121621)) -- rats given a single oral gavage of Spatholobi Caulis aqueous extract.

- Add Al-Zubaydi 2024 gabapentin ([doi:10.3390/pharmaceutics16121514](https://doi.org/10.3390/pharmaceutics16121514)) -- hospitalized adults with therapeutic drug monitoring concentrations.

- Add Lee 2024 L-serine (AST-001) PK-PD ([doi:10.3389/fphar.2024.1452526](https://doi.org/10.3389/fphar.2024.1452526)) -- children aged 2-11 years with autism spectrum disorder.

- Add Yan 2024 amisulpride ([doi:10.2147/DDDT.S469149](https://doi.org/10.2147/DDDT.S469149)) - Chinese adult inpatients with schizophrenia.

- Add Mao 2024 sirolimus ([doi:10.2147/DDDT.S503463](https://doi.org/10.2147/DDDT.S503463)) -- adult liver transplant recipients.

- Add Jiang 2025 filgrastim ([doi:10.1111/cts.70121](https://doi.org/10.1111/cts.70121)) -- healthy Korean adult men.

- Add Hernandez-Lozano 2025 apramycin translational PKPD models ([doi:10.1093/jac/dkae409](https://doi.org/10.1093/jac/dkae409)) -- in vitro (Escherichia coli EN591 and ATCC 700336), mouse (C3H/HeJ) complicated urinary tract infection, and a simulated human efficacy prediction.

- Add Comisar 2025 zavegepant ([doi:10.1002/psp4.13257](https://doi.org/10.1002/psp4.13257)) -- healthy adults and patients with migraine dosed intravenously, intranasally, or orally.

- Add Izat 2025 zaleplon ([doi:10.1002/psp4.13255](https://doi.org/10.1002/psp4.13255)) -- healthy adult women and men, intravenous.

- Add Izat 2025 ziprasidone ([doi:10.1002/psp4.13255](https://doi.org/10.1002/psp4.13255)) -- healthy adult men, intravenous.

- Add Izat 2025 zoniporide ([doi:10.1002/psp4.13255](https://doi.org/10.1002/psp4.13255)) -- healthy adult men, intravenous.

- Add Eaton 2025 dabigatran ([doi:10.1177/02676591231226291](https://doi.org/10.1177/02676591231226291)) -- anaesthetised sheep, with idarucizumab reversal.

- Add Cerqueira 2025 resveratrol ([doi:10.3390/nu17010181](https://doi.org/10.3390/nu17010181)) -- rat (Wistar, male); separate intravenous and oral models.

- Add Centanni 2025 sunitinib thrombocytopenia ([doi:10.1007/s40273-024-01438-z](https://doi.org/10.1007/s40273-024-01438-z)) -- adults with imatinib-resistant gastrointestinal stromal tumours (GIST).

- Add Zhang 2025 epoetin alfa erythroferrone PK/PD models ([doi:10.1021/acsptsci.4c00575](https://doi.org/10.1021/acsptsci.4c00575)) -- rat (Sprague-Dawley) models of adenine-induced CKD anemia and carboplatin-induced anemia.

- Add Marques 2025 oral salbutamol ([doi:10.3390/pharmaceutics17010039](https://doi.org/10.3390/pharmaceutics17010039)) -- virtual (PBPK-generated) patients aged 5-65 years.

- Add Gaspar 2025 fexofenadine ([doi:10.1007/s40262-024-01470-4](https://doi.org/10.1007/s40262-024-01470-4)) -- hospitalized older adult polymorbid patients pooled with healthy volunteers.

- Add Darwish 2025 trofinetide ([doi:10.1007/s12325-024-03056-9](https://doi.org/10.1007/s12325-024-03056-9)) -- healthy volunteers and patients with Rett syndrome, fragile X syndrome, or traumatic brain injury.

- Add Golzaryan 2025 [177Lu]Lu-PSMA I&T ([doi:10.1038/s41598-025-86159-9](https://doi.org/10.1038/s41598-025-86159-9)) - men with metastatic castration-resistant prostate cancer; a 21-compartment whole-body PBPK model with parallel labelled and unlabelled circulations.

- Add Lam 2025 ondansetron ([doi:10.1111/cts.70147](https://doi.org/10.1111/cts.70147)) -- neonates with neonatal opioid withdrawal syndrome.

- Add Kim 2025 evogliptin ([doi:10.1002/psp4.13263](https://doi.org/10.1002/psp4.13263)) -- adults spanning normal renal function through end-stage renal disease on hemodialysis.

- Add Braniff 2025 lorlatinib signaling / shell-and-core tumor QSP model ([doi:10.1002/psp4.13270](https://doi.org/10.1002/psp4.13270)) -- adults with ALK-positive advanced non-small cell lung cancer.

- Add Crass 2025 pegcetacoplan geographic-atrophy lesion-area models
  ([doi:10.1002/psp4.13264](https://doi.org/10.1002/psp4.13264))
  -- adults with geographic atrophy secondary to age-related macular degeneration.

- Add Balis 2025 busulfan ([doi:10.1007/s00280-025-04757-w](https://doi.org/10.1007/s00280-025-04757-w)) -- infants and children receiving busulfan-containing stem cell transplant conditioning.

- Add Park 2025 anti-CD19 CAR-T cell therapy in systemic lupus erythematosus ([doi:10.1111/cts.70146](https://doi.org/10.1111/cts.70146)) -- adults with severe refractory SLE.

- Add Brossard 2025 emapalumab ([doi:10.1111/cts.70163](https://doi.org/10.1111/cts.70163)) -- children and young adults with macrophage activation syndrome associated with Still's disease.

- Add Nakai 2025 tranexamic acid ([doi:10.1007/s00228-025-03802-0](https://doi.org/10.1007/s00228-025-03802-0)) -- adults undergoing cardiac surgery with cardiopulmonary bypass.

- Add Deferm 2025 postpartum maternal physiology and milk composition ([doi:10.3389/fphar.2025.1517069](https://doi.org/10.3389/fphar.2025.1517069)) -- healthy breastfeeding women from birth to 12 months postpartum.

- Add Nie 2025 epidural sufentanil ([doi:10.2147/DDDT.S500189](https://doi.org/10.2147/DDDT.S500189)) -- primiparous women receiving patient-controlled epidural labour analgesia.

- Add DeJongh 2025 AZD7648 and olaparib ([doi:10.1007/s10928-025-09962-x](https://doi.org/10.1007/s10928-025-09962-x)) -- mice (SCID and athymic nude, FaDu ATM-knockout xenograft).

- Add Rolsma 2025 cefepime, meropenem, piperacillin and tazobactam ([doi:10.1093/infdis/jiae451](https://doi.org/10.1093/infdis/jiae451)) -- children and adults with cystic fibrosis.

- Add Baiardi 2025 dalbavancin ([doi:10.3390/antibiotics14020190](https://doi.org/10.3390/antibiotics14020190)) -- adults receiving multidose long-term therapy for difficult-to-treat Gram-positive infections.

- Add Han 2025 midazolam, fentanyl, alfentanil and sufentanil age-dependent whole-body PBPK models ([doi:10.3390/pharmaceutics17020214](https://doi.org/10.3390/pharmaceutics17020214)) -- preterm neonates through adults to the oldest old.

- Add Olafuyi 2025 propylene glycol ([doi:10.1002/jcph.6150](https://doi.org/10.1002/jcph.6150)) - healthy adults and term neonates; compartmental reductions of the published Simcyp full-PBPK models, with saturable ADH-mediated clearance.

- Add Wu 2025 sitafloxacin ([doi:10.3389/fphar.2025.1476158](https://doi.org/10.3389/fphar.2025.1476158)) -- Japanese and Chinese healthy volunteers, subjects with renal impairment, and patients with respiratory-tract infection.

- Add Fan 2025 intravenous iron whole-body PBPK models ([doi:10.1007/s13346-024-01675-x](https://doi.org/10.1007/s13346-024-01675-x)) -- mice on iron-deficient, iron-adequate and iron-loaded diets, plus ferric carboxymaltose extrapolations to rat (iron-deficiency anaemia) and adults with iron-deficiency anaemia.

- Add Butragueno-Laiseca 2025 teicoplanin ([doi:10.1093/jac/dkaf012](https://doi.org/10.1093/jac/dkaf012)) -- critically ill children, including those on continuous kidney replacement therapy.

- Add Xie 2025 midazolam ([doi:10.2147/DDDT.S495647](https://doi.org/10.2147/DDDT.S495647)) -- mechanically ventilated Chinese adult ICU patients.

- Add Takada 2025 vancomycin ([doi:10.1186/s40780-025-00423-8](https://doi.org/10.1186/s40780-025-00423-8)) -- Japanese inpatients aged 75 years and older with a body mass index below 25 kg/m^2.

- Add Li 2025 doxycycline PK and doxycycline + florfenicol PK/PD models ([doi:10.1016/j.psj.2025.104922](https://doi.org/10.1016/j.psj.2025.104922)) -- Riemerella anatipestifer-infected ducks (Tadorna tadorna).

- Add Zhang 2025 abemaciclib CDK4/6 occupancy and pRB / TOPO-IIa biomarker QSP model ([doi:10.1021/acsomega.4c09472](https://doi.org/10.1021/acsomega.4c09472)) -- patients with metastatic breast cancer, including brain metastases.

- Add Winning 2025 certepetide ([doi:10.1002/cpdd.1502](https://doi.org/10.1002/cpdd.1502)) -- adults with metastatic pancreatic ductal adenocarcinoma.

- Add Ji 2025 TDI01 ([doi:10.3389/fphar.2025.1477607](https://doi.org/10.3389/fphar.2025.1477607)) -- healthy Chinese adults.

- Add Said 2025 imatinib ([doi:10.1002/psp4.13299](https://doi.org/10.1002/psp4.13299)) -- adults with COVID-19 ARDS pooled with CML/GIST oncology patients.

- Add Michelet 2025 BI 754111 anti-LAG-3 two-pore minimal-PBPK model ([doi:10.1002/psp4.13285](https://doi.org/10.1002/psp4.13285)) -- adults with anti-PD-1-refractory NSCLC or HNSCC in an 89Zr-immuno-PET biodistribution study.

- Add Boone 2025 vinyl chloride PBPK ([doi:10.3390/jox15020042](https://doi.org/10.3390/jox15020042)) -- 70 kg reference adult; community exposure reconstruction after the 2012 Paulsboro, New Jersey train derailment.

- Add Miano 2024 tacrolimus ([doi:10.1016/j.jhlto.2024.100134](https://doi.org/10.1016/j.jhlto.2024.100134)) -- adult lung transplant recipients during the first 14 postoperative days.

- Add Budiansah 2025 DOTATATE ([doi:10.1186/s40658-025-00726-7](https://doi.org/10.1186/s40658-025-00726-7)) -- adults with neuroendocrine tumours or meningiomas undergoing peptide receptor radionuclide therapy dosimetry.

- Add Warren 2025 apremilast and orismilast ([doi:10.1007/s13555-025-01371-9](https://doi.org/10.1007/s13555-025-01371-9)) -- adults with atopic dermatitis, plus an in vitro human whole-blood IL-13 inhibition layer for each drug.

- Add Zhang 2025 fluconazole ([doi:10.3389/fphar.2025.1564070](https://doi.org/10.3389/fphar.2025.1564070)) -- critically ill adults with acute renal failure receiving continuous renal replacement therapy.

- Add Visscher 2025 parathyroid hormone rhPTH(1-84) ([doi:10.1111/bcp.16342](https://doi.org/10.1111/bcp.16342)) -- a single adult with chronic postsurgical hypoparathyroidism.

- Add Wen 2025 salbutamol pulmonary PBPK ([doi:10.1002/psp4.13316](https://doi.org/10.1002/psp4.13316)) -- rat (male Wistar Han) after intratracheal instillation.

- Add Demin 2025 zanubrutinib, acalabrutinib and ibrutinib BTK-occupancy QSP model ([doi:10.1002/psp4.13307](https://doi.org/10.1002/psp4.13307)) -- virtual patients with B-cell malignancies.

- Add Zhang 2025 cefiderocol ([doi:10.1007/s12325-025-03147-1](https://doi.org/10.1007/s12325-025-03147-1)) - healthy Chinese adults.

- Add Adamiszak 2025 fluconazole ([doi:10.3390/pharmaceutics17040488](https://doi.org/10.3390/pharmaceutics17040488)) -- hemato-oncologic pediatric patients aged 7 months to 18 years.

- Add Na 2025 HOSU-53 (JBZ-001) PK/PD models ([doi:10.3390/pharmaceutics17040412](https://doi.org/10.3390/pharmaceutics17040412)) -- mice and beagle dogs.

- Add Wickramasinghe 2025 pamiparib ([doi:10.3390/pharmaceutics17040524](https://doi.org/10.3390/pharmaceutics17040524)) -- adults with newly diagnosed or recurrent glioblastoma.

- Add Gao 2025 cefquinome ([doi:10.3390/vetsci12040294](https://doi.org/10.3390/vetsci12040294)) -- Ili foals, with an ex vivo PK/PD-integration model against *Escherichia coli*.

- Add Fromage 2025 mycophenolic acid ([doi:10.1111/bcp.16374](https://doi.org/10.1111/bcp.16374)) -- solid organ transplant, haematopoietic cell transplant and autoimmune disease patients on enteric-coated mycophenolate sodium.

- Add Athanassa 2025 minocycline ([doi:10.1093/jac/dkaf090](https://doi.org/10.1093/jac/dkaf090)) -- critically ill adults with ventilator-associated pneumonia.

- Add Lacroix 2025 polymyxin B ([doi:10.1128/aac.01535-24](https://doi.org/10.1128/aac.01535-24)) -- in vitro time-kill against two multidrug-resistant Acinetobacter baumannii clinical isolates, with and without 1% mucin.

- Add Van Wart 2025 tobramycin ([doi:10.1128/aac.00908-24](https://doi.org/10.1128/aac.00908-24)) -- adults with pneumonia, with serum and epithelial lining fluid concentrations.

- Add Roberts 2025 remdesivir and GS-441524 ([doi:10.1007/s40262-025-01496-2](https://doi.org/10.1007/s40262-025-01496-2)) -- hospitalised adults with COVID-19.

- Add Mauro 2025 nilotinib ([doi:10.1007/s00280-025-04777-6](https://doi.org/10.1007/s00280-025-04777-6)) -- healthy adults receiving Danziten tablets or Tasigna capsules under four prandial conditions.

- Add Wang 2025 rivaroxaban ([doi:10.3389/fphar.2025.1574949](https://doi.org/10.3389/fphar.2025.1574949)) - Chinese adults with non-valvular atrial fibrillation.

- Add Song 2025 infliximab ([doi:10.5009/gnl240503](https://doi.org/10.5009/gnl240503)) -- Korean adults with inflammatory bowel disease on maintenance intravenous or subcutaneous therapy.

- Add Garcia 2025 garadacimab ([doi:10.1002/psp4.70009](https://doi.org/10.1002/psp4.70009)) -- healthy volunteers and adolescents and adults with hereditary angioedema.

- Add van der Gaag 2025 osimertinib target-site PBPK model ([doi:10.1002/psp4.70006](https://doi.org/10.1002/psp4.70006)) -- adults with advanced-stage EGFR-mutated non-small cell lung cancer.

- Add Wickramasinghe 2025 abemaciclib ([doi:10.1002/psp4.70026](https://doi.org/10.1002/psp4.70026)) -- 9-compartment spatial CNS and brain-tumor PBPK model in adults with glioblastoma.

- Add Patel 2025 eteplirsen ([doi:10.1002/psp4.70001](https://doi.org/10.1002/psp4.70001)) -- boys with Duchenne muscular dystrophy aged 6 months to 16 years.

- Add Assmus 2025 benznidazole ([doi:10.1371/journal.pntd.0012968](https://doi.org/10.1371/journal.pntd.0012968)) -- mouse (BALB/c), uninfected satellite PK cohort and chronic *Trypanosoma cruzi* infection efficacy cohort.

- Add McCann 2025 midazolam ([doi:10.1111/cts.70247](https://doi.org/10.1111/cts.70247)) -- children and young adults with and without obesity receiving standard-of-care IV midazolam.

- Add Kobuchi 2025 dapagliflozin ([doi:10.7150/ijms.111519](https://doi.org/10.7150/ijms.111519)) -- Japanese outpatients with type 2 diabetes mellitus treated for one year in routine practice.

- Add Zhao 2025 methotrexate ([doi:10.3389/fphar.2025.1548203](https://doi.org/10.3389/fphar.2025.1548203)) -- Chinese children and adults with intracranial germ cell tumors.

- Add Kir 2025 atenolol and metoprolol minimal-PBPK absorption models ([doi:10.1007/s13318-025-00943-6](https://doi.org/10.1007/s13318-025-00943-6)) -- non-malnourished and malnourished rats (Sprague-Dawley).

- Add Wolowich 2025 delta-9-tetrahydrocannabinol heart-rate PK/PD models ([doi:10.1007/s13318-025-00941-8](https://doi.org/10.1007/s13318-025-00941-8)) -- healthy volunteers given a single intravenous bolus.

- Add Sheiner 1979 d-tubocurarine ([doi:10.1002/cpt1979253358](https://doi.org/10.1002/cpt1979253358)) -- adults undergoing elective surgery with normal renal function or chronic end-stage renal failure.

- Add Perlstein 2025 TV-46000 long-acting subcutaneous risperidone ([doi:10.1007/s40120-025-00723-z](https://doi.org/10.1007/s40120-025-00723-z)) -- healthy volunteers and adults with schizophrenia or schizoaffective disorder.

- Add Lee 2025 levofloxacin ([doi:10.3390/ph18050621](https://doi.org/10.3390/ph18050621)) -- healthy Korean adults given a single 500 mg intravenous infusion.

- Add de Cacqueray 2022 cefepime ([doi:10.1016/j.cmi.2022.05.007](https://doi.org/10.1016/j.cmi.2022.05.007)) -- critically ill infants and children aged 1.1 months to 17.6 years.

- Add Zhao 2020 cefepime ([doi:10.3389/fphar.2020.00014](https://doi.org/10.3389/fphar.2020.00014)) -- Chinese neonates and young infants with postmenstrual age below 48 weeks.

- Add van den Maagdenberg 2025 A2AR / anti-PD-L1 tumour-microenvironment QSP model ([doi:10.1021/acs.jcim.5c00107](https://doi.org/10.1021/acs.jcim.5c00107)) -- mouse (MCA205 syngeneic tumour model).

- Add Onita 2025 sulbactam ([doi:10.1093/jpids/piaf043](https://doi.org/10.1093/jpids/piaf043)) -- pediatric patients from 4 weeks to 16 years, pooled from 23 published studies.

- Add van den Berg 2025 generic monoclonal antibody meta-analytic models ([doi:10.1080/19420862.2025.2512217](https://doi.org/10.1080/19420862.2025.2512217)) -- medians across 160 published population PK models of 69 marketed IgG mAbs in adults.

- Add Valadez 2025 cefepime ([doi:10.1128/aac.00102-25](https://doi.org/10.1128/aac.00102-25)) -- mechanically ventilated adults in the ICU with suspected hospital-acquired pneumonia, with and without ECMO.

- Add Okumura 2025 ergothioneine PBPK model ([doi:10.1002/fsn3.70382](https://doi.org/10.1002/fsn3.70382)) -- healthy Japanese adults taking daily oral ergothioneine supplements.

- Add Xia 2025 p-furoylamphenmulin ([doi:10.1016/j.psj.2025.105249](https://doi.org/10.1016/j.psj.2025.105249)) -- Mycoplasma gallisepticum-infected specific-pathogen-free chickens.

- Add Gu 2025 rivaroxaban ([doi:10.3389/fphar.2025.1562259](https://doi.org/10.3389/fphar.2025.1562259)) - Chinese healthy volunteers and patients after radiofrequency ablation for non-valvular atrial fibrillation.

- Add Hornik 2025 furosemide ([doi:10.1007/s40262-025-01515-2](https://doi.org/10.1007/s40262-025-01515-2)) -- adults with chronic heart failure and volume overload, allometrically scaled to adolescents.

- Add Shigetome 2025 paroxetine ([doi:10.1002/psp4.70032](https://doi.org/10.1002/psp4.70032)) -- Japanese adults with major depressive disorder; population PK plus a MADRS exposure-response model.

- Add Cao 2025 ferric carboxymaltose ([doi:10.1021/acsptsci.5c00097](https://doi.org/10.1021/acsptsci.5c00097)) -- rats (Sprague-Dawley) with iron deficiency anemia.

- Add Glatard 2025 octreotide ([doi:10.1007/s40262-025-01522-3](https://doi.org/10.1007/s40262-025-01522-3)) -- healthy adults and adults with acromegaly.

- Fix the dose basis in the Steichert 2025 enalapril + enalaprilat vignette
  ([doi:10.1007/s40262-025-01520-5](https://doi.org/10.1007/s40262-025-01520-5)) --
  doses are enalapril free base, so 0.25 mg enalapril maleate is 191.1 ug, not
  250 ug. The vignette now reproduces the published enalaprilat Cmax,1 instead of
  overshooting it by 31%.

- Add Shah 2025 clarithromycin ([doi:10.3390/antibiotics14060559](https://doi.org/10.3390/antibiotics14060559)) -- critically ill adults receiving intravenous clarithromycin.

- Add Liu 2025 avatrombopag ([doi:10.3390/ph18060903](https://doi.org/10.3390/ph18060903)) -- healthy Chinese adults dosed fasting or after a high-fat meal.

- Add Zhao 2025 paracetamol ([doi:10.1002/bcp.70028](https://doi.org/10.1002/bcp.70028)) - children and adults with spinal muscular atrophy and healthy controls.

- Add Van Wart 2025 telavancin ([doi:10.1128/aac.01382-24](https://doi.org/10.1128/aac.01382-24)) - healthy subjects and patients with complicated skin and skin-structure infection, hospital-acquired or ventilator-associated bacterial pneumonia, or uncomplicated bacteremia, spanning the full range of renal function including hemodialysis.

- Add Harada 2025 cyanide ([doi:10.1007/s11419-025-00713-8](https://doi.org/10.1007/s11419-025-00713-8)) -- adults in fire-related deaths examined at forensic autopsy.

- Add Dasti 2025 multiscale mRNA vaccine QSP models for BNT162b2 (general adult
  and over-60 populations), mRNA-1273, and the single-cell antigen-presentation
  molecular layer
  ([doi:10.1002/psp4.70041](https://doi.org/10.1002/psp4.70041))
  -- healthy adults receiving prophylactic COVID-19 mRNA vaccination; the early
  antigen-presenting-cell events are calibrated on rhesus macaque data.

- Add Toutain 2025 doxycycline ([doi:10.1111/jvp.13511](https://doi.org/10.1111/jvp.13511)) -- pigs of 8.5-100.6 kg dosed intravenously or orally in feed, drinking water or by stomach tube.

- Add Kong 2025 piperacillin/tazobactam ([doi:10.1007/s40262-025-01527-y](https://doi.org/10.1007/s40262-025-01527-y)) -- adults with end-stage kidney disease on intermittent haemodialysis.

- Add Fiandaca 2025 mRNA-encoded BiTE short-chain and long-chain multiscale PBPK models ([doi:10.1016/j.omtn.2025.102606](https://doi.org/10.1016/j.omtn.2025.102606)) -- mouse (immunodeficient, tumour-bearing).

- Add Chen 2025 iohexol + creatinine joint GFR / OCT2-MATE model ([doi:10.1002/cpt.3612](https://doi.org/10.1002/cpt.3612)) -- healthy adult volunteers.

- Add Li 2025 modakafusp alfa ([doi:10.1111/cts.70296](https://doi.org/10.1111/cts.70296)) - adults with relapsed or refractory multiple myeloma.

- Add Chen 2025 remimazolam ([doi:10.3389/fphar.2025.1526266](https://doi.org/10.3389/fphar.2025.1526266)) -- critically ill adults receiving continuous infusion for ICU sedation.

- Add Rognas 2025 bitopertin ([doi:10.1007/s10928-025-09990-7](https://doi.org/10.1007/s10928-025-09990-7)) -- healthy adults (ratifies seven new erythropoiesis compartment canonicals: `ret_imm_marrow`, `ret_mat_marrow`, `ret_imm_blood`, `ret_mat_blood`, and the numbered-chain prefixes `erythrocytes<n>`, `mch<n>` and `moderator<n>`).

- Add Ebihara 2025 tebipenem ([doi:10.3390/antibiotics14070648](https://doi.org/10.3390/antibiotics14070648)) -- Japanese adults, stratified by renal function.

- Add He 2025 lidocaine ([doi:10.2147/DDDT.S485389](https://doi.org/10.2147/DDDT.S485389)) -- Chinese adults undergoing partial hepatectomy, with the MEGX and GX active metabolites.

- Add Richardson 2025 automated model-search popPK models for ribociclib, camizestrant, osimertinib, olaparib and tezepelumab ([doi:10.1038/s43856-025-01054-8](https://doi.org/10.1038/s43856-025-01054-8)) -- one synthetic cohort and four pooled Phase 1 adult cohorts.

- Add Xie 2025 aztreonam-avibactam ([doi:10.1128/aac.01950-24](https://doi.org/10.1128/aac.01950-24)) -- healthy adults and hospitalized adults with complicated intra-abdominal infection, nosocomial pneumonia, complicated urinary tract infection, or bloodstream infection, spanning augmented renal clearance through end-stage renal disease.

- Add Kong 2025 sudapyridine ([doi:10.1016/j.ejps.2025.107160](https://doi.org/10.1016/j.ejps.2025.107160)) - Chinese healthy volunteers and drug-susceptible / multidrug-resistant tuberculosis patients.

- Add Comisar 2025 rimegepant ([doi:10.1002/psp4.70051](https://doi.org/10.1002/psp4.70051)) -- healthy adults, elderly adults with stable chronic illness, and adults with renal or hepatic impairment.

- Add Zhang 2025 dupilumab FEV1 PK/PD ([doi:10.1002/psp4.70057](https://doi.org/10.1002/psp4.70057)) -- adults and adolescents with uncontrolled moderate-to-severe asthma.

- Add Vucicevic 2025 efavirenz ([doi:10.7759/cureus.88533](https://doi.org/10.7759/cureus.88533)) -- Caucasian adults with HIV-1 in Serbia.

- Add Lee 2025 piperacillin and tazobactam ([doi:10.3390/ph18081124](https://doi.org/10.3390/ph18081124)) - healthy Korean adults.

- Add Rong 2019 and van Hest 2005 mycophenolic acid models, as re-implemented by Maizaud 2025 ([doi:10.1371/journal.pone.0330854](https://doi.org/10.1371/journal.pone.0330854)) - adult kidney transplant recipients on tacrolimus (Rong) or ciclosporin (van Hest).

- Add Poels 2025 elranatamab QSP model ([doi:10.1038/s41540-025-00585-z](https://doi.org/10.1038/s41540-025-00585-z)) -- adults with relapsed or refractory multiple myeloma.

- Add Randell 2024 metronidazole ([doi:10.1128/aac.01533-23](https://doi.org/10.1128/aac.01533-23)) -- critically ill preterm and term infants.

- Add Zavrelova 2025 linezolid ([doi:10.1111/cts.70346](https://doi.org/10.1111/cts.70346)) -- hematooncological adults with suspected or proven Gram-positive sepsis.

- Add Zhang 2025 nedosiran ([doi:10.1007/s40262-025-01540-1](https://doi.org/10.1007/s40262-025-01540-1)) -- healthy volunteers and children and adults with primary hyperoxaluria.

- Add Francke 2025 tacrolimus full and starting-dose models ([doi:10.1007/s40262-025-01533-0](https://doi.org/10.1007/s40262-025-01533-0)) - adult living and deceased donor kidney transplant recipients.

- Add Huang 2025 dexmedetomidine nasal spray ([doi:10.3389/fphar.2025.1662364](https://doi.org/10.3389/fphar.2025.1662364)) -- Chinese healthy volunteers and adults undergoing elective abdominal surgery.

- Add Yamada 2025 fluorouracil and oxaliplatin ([doi:10.1007/s00280-025-04808-2](https://doi.org/10.1007/s00280-025-04808-2)) -- adults with locally advanced unresectable or metastatic gastric or gastroesophageal junction adenocarcinoma receiving mFOLFOX6 with or without zolbetuximab.

- Add Xiang 2025 tacrolimus ([doi:10.2147/DDDT.S542786](https://doi.org/10.2147/DDDT.S542786)) -- adult renal transplant recipients.

- Add 14 published imatinib population PK models transcribed from the Yang 2025 external evaluation ([doi:10.1007/s11523-025-01172-2](https://doi.org/10.1007/s11523-025-01172-2)) -- Judson 2005, Schmidli 2005, Widmer 2006, Petain 2008, Demetri 2009, Menon-Andersen 2009, Yamakawa 2011, Eechoute 2012, Golabchifar 2014, Gotta 2014, Di Paolo 2014, Wang 2019, Shriyan 2022 and He 2023; adults with chronic myeloid leukemia or gastrointestinal stromal tumor, two of them also including children.

- Add Bouazza 2025 prednisolone ([doi:10.1002/bcp.70103](https://doi.org/10.1002/bcp.70103)) -- paediatric and adult patients with active systemic lupus erythematosus.

- Add Mohammed Ali 2025 tacrolimus ([doi:10.3390/pharmaceutics17091185](https://doi.org/10.3390/pharmaceutics17091185)) -- stable adult renal transplant recipients converted from IR-Tac to LCP-Tac.

- Add Shen 2025 voriconazole ([doi:10.3389/fphar.2025.1671652](https://doi.org/10.3389/fphar.2025.1671652)) -- immunocompromised children under 2 years of age.

- Add de Vries 2025 durvalumab ([doi:10.1007/s40262-025-01555-8](https://doi.org/10.1007/s40262-025-01555-8)) -- adults with non-small cell lung cancer.

- Add Boulanger 2025 trimethoprim + sulfadiazine / sulfadimethoxine / sulfamethoxazole ([doi:10.1080/01652176.2025.2565351](https://doi.org/10.1080/01652176.2025.2565351)) -- healthy growing pigs.

- Add Salehi 2025 nicotine ([doi:10.1002/jcph.70038](https://doi.org/10.1002/jcph.70038)) -- adults who use nicotine pouches or moist smokeless tobacco.

- Add Resendiz-Galvan 2025 cycloserine ([doi:10.1128/aac.00101-25](https://doi.org/10.1128/aac.00101-25)) - Indian adolescents and adults with multidrug-resistant tuberculosis.

- Add Cendros 2025 enflicoxib ([doi:10.3389/fvets.2025.1645857](https://doi.org/10.3389/fvets.2025.1645857)) -- dogs with naturally occurring osteoarthritis.

- Add Tao 2025 meropenem ([doi:10.3389/fphar.2025.1643553](https://doi.org/10.3389/fphar.2025.1643553)) -- critically ill adult and elderly surgical-ICU patients with Pseudomonas aeruginosa infections.

- Add Yao 2025 flurbiprofen enantiomers ([doi:10.2147/DDDT.S542722](https://doi.org/10.2147/DDDT.S542722)) -- Chinese adults with postoperative pain after joint replacement, with paired plasma and cerebrospinal-fluid sampling.

- Add Mo 2025 everolimus ([doi:10.12793/tcp.2025.33.e14](https://doi.org/10.12793/tcp.2025.33.e14)) -- healthy adult Korean males.

- Add Kim 2025 infliximab, a panel of eight externally validated population PK models ([doi:10.1002/psp4.70089](https://doi.org/10.1002/psp4.70089)) -- adults and children with inflammatory bowel disease.

- Add McBride 2025 recombinant ADAMTS13 ([doi:10.1002/psp4.70063](https://doi.org/10.1002/psp4.70063)) -- patients with congenital thrombotic thrombocytopenic purpura.

- Add Duffull 2025 monoclonal antibody TMDD models, full and target-saturated simplification ([doi:10.1002/psp4.70049](https://doi.org/10.1002/psp4.70049)) -- 80 simulated subjects, unnamed mAb.

- Add Comisar 2025 rimegepant ([doi:10.1111/cts.70360](https://doi.org/10.1111/cts.70360)) -- pooled adults and children 6 to <12 years with a history of migraine.

- Add Tsyplakova 2025 mycophenolate sodium and mycophenolate mofetil ([doi:10.3390/medsci13040235](https://doi.org/10.3390/medsci13040235)) -- adult renal transplant recipients.

- Add Karakitsios 2025 bedaquiline lung PBPK models ([doi:10.1002/bcp.70163](https://doi.org/10.1002/bcp.70163)) -- mouse, rat (Sprague-Dawley), beagle dog, and adults with pulmonary drug-resistant tuberculosis.

- Add Rosenborg 2025 inhaled fluticasone propionate and salmeterol ([doi:10.2147/DDDT.S480189](https://doi.org/10.2147/DDDT.S480189)) - healthy adults in three dry-powder-inhaler bioequivalence crossover studies.

- Add Jia 2025 TAK-071 ([doi:10.1002/cpdd.1579](https://doi.org/10.1002/cpdd.1579)) -- healthy adults and participants with Parkinson disease and cognitive impairment.

- Add Fan 2025 Nb457-NbHSA-Nb457 anti-CD4 trimeric nanobody and ibalizumab TMDD PK-PD models ([doi:10.1128/spectrum.00805-25](https://doi.org/10.1128/spectrum.00805-25)) -- HIV-1-infected humanized NDG-HuPBL mice, plus allometric human projections.

- Add Stoschus 2025 phenobarbital ([doi:10.1111/epi.18517](https://doi.org/10.1111/epi.18517)) -- critically ill adults with refractory and superrefractory status epilepticus.

- Add Parkinson 2025 balcinrenone ([doi:10.1007/s40262-025-01572-7](https://doi.org/10.1007/s40262-025-01572-7)) -- healthy participants, participants with renal impairment, and patients with heart failure and chronic kidney disease.

- Add Ozdin 2025 dexamethasone, erythrocyte-encapsulated ([doi:10.1002/psp4.70103](https://doi.org/10.1002/psp4.70103)) -- healthy adults and pediatric patients with ataxia telangiectasia.

- Add Nguyen 2025 valbenazine ([doi:10.1002/jcph.70092](https://doi.org/10.1002/jcph.70092)) -- healthy adults and patients with tardive dyskinesia or Huntington's disease chorea.

- Add Gu 2025 whole-body heart-failure PBPK models for digoxin, furosemide, bumetanide, torasemide, captopril, valsartan, felodipine and midazolam ([doi:10.3390/pharmaceutics17111394](https://doi.org/10.3390/pharmaceutics17111394)) -- healthy adults and chronic heart-failure patients in NYHA classes II-IV.

- Add Zhou 2025 tacrolimus ([doi:10.1007/s00228-025-03920-9](https://doi.org/10.1007/s00228-025-03920-9)) -- Chinese adult lung transplant recipients in the first three months post-transplantation.

- Add Jian 2025 peginterferon alfa-2b / Pegbing ([doi:10.1002/psp4.70104](https://doi.org/10.1002/psp4.70104)) -- healthy volunteers and chronic hepatitis B patients.

- Add Tian 2025 pirtobrutinib ([doi:10.1002/psp4.70134](https://doi.org/10.1002/psp4.70134)) - healthy adults; a compartmental reduction of the published Simcyp minimal-PBPK model.

- Add Xu 2025 aficamten ([doi:10.1002/psp4.70099](https://doi.org/10.1002/psp4.70099)) - healthy adults and adults with obstructive hypertrophic cardiomyopathy.

- Add Collins 2025 midazolam ([doi:10.1002/psp4.70116](https://doi.org/10.1002/psp4.70116)) -- healthy volunteers given oral midazolam as a CYP3A probe with single-dose and multiple-dose efavirenz.

- Add Varela-Gonzalez-Aller 2025 fludarabine ([doi:10.3390/pharmaceutics17121592](https://doi.org/10.3390/pharmaceutics17121592)) -- adults with relapsed/refractory large B-cell lymphoma receiving fludarabine-based lymphodepletion before CAR T-cell therapy.

- Add Trujillo 2025 proinsulin/glucose homeostasis QSP model ([doi:10.3390/pharmaceutics17121522](https://doi.org/10.3390/pharmaceutics17121522)) -- four virtual patients spanning healthy physiology and the type 2 diabetes spectrum.

- Add Zhang 2026 tebipenem ([doi:10.1111/cts.70453](https://doi.org/10.1111/cts.70453)) -- Japanese children aged 0.5-16 years, simulated for Bangladeshi children aged 24-59 months with shigellosis.

- Add Sinha 2026 oxcarbazepine ([doi:10.1007/s40262-025-01579-0](https://doi.org/10.1007/s40262-025-01579-0)) -- children and young adults aged 44 days to 20.9 years, 52% with obesity.

- Add Yu 2025 vancomycin ([doi:10.1007/s40268-025-00523-8](https://doi.org/10.1007/s40268-025-00523-8)) -- Chinese pediatric inpatients (birth to 15 years) receiving IV vancomycin.

- Add Cheng 2026 levamisole interspecies models ([doi:10.1007/s40005-025-00770-6](https://doi.org/10.1007/s40005-025-00770-6)) -- joint allometric two-compartment and minimal-PBPK fits across duck, rabbit, chicken, goat, dog, sheep, pig and human.

- Add Morath 2025 apixaban ([doi:10.1007/s40262-025-01534-z](https://doi.org/10.1007/s40262-025-01534-z)) -- adults with postoperative atrial fibrillation after cardiac surgery, with and without concomitant amiodarone.

- Add Decker 2026 baricitinib ([doi:10.1007/s40262-025-01563-8](https://doi.org/10.1007/s40262-025-01563-8)) -- pediatric patients aged 2 to <18 years with moderate-to-severe atopic dermatitis.

- Add Zakria 2026 voriconazole ([doi:10.1080/20523211.2025.2601420](https://doi.org/10.1080/20523211.2025.2601420)) -- adult Pakistani cancer patients, elderly (> 65 y) versus young.

- Add Li 2025 amikacin ([doi:10.1186/s12879-025-11747-z](https://doi.org/10.1186/s12879-025-11747-z)) -- Chinese premature infants receiving therapeutic drug monitoring.

- Add Zeng 2026 atezolizumab ([doi:10.1007/s00228-025-03974-9](https://doi.org/10.1007/s00228-025-03974-9)) -- adults with metastatic non-small cell lung cancer; all parameters fixed, transcribed from the FDA CDER review of BLA 761041Orig1s000.

- Add Kollo 2026 selumetinib ([doi:10.1002/psp4.70156](https://doi.org/10.1002/psp4.70156)) - children aged 5-16 years with inoperable neurofibromatosis type I or plexiform neurofibromas.

- Add Moein 2024 apitolisib translational PK/PD models ([doi:10.1007/s40268-024-00459-5](https://doi.org/10.1007/s40268-024-00459-5)) -- mouse (786-O renal cell adenocarcinoma xenograft) and adults with advanced solid tumors or non-Hodgkin's lymphoma.

- Add Decrane 2023 oxyfluorfen ([doi:10.1016/j.crtox.2023.100138](https://doi.org/10.1016/j.crtox.2023.100138)) -- rat (Sprague-Dawley) and extrapolated human.

- Add Courlet 2023 cabamiquine ([doi:10.1128/aac.00891-23](https://doi.org/10.1128/aac.00891-23)) - healthy adult men in induced blood stage and sporozoite malaria challenge studies.

- Add Hanley 2024 brigatinib ([doi:10.1002/psp4.13106](https://doi.org/10.1002/psp4.13106)) - healthy adults; a compartmental reduction of the published Simcyp minimal-PBPK model.

- Add Majid 2024 lenvatinib ([doi:10.1002/psp4.13130](https://doi.org/10.1002/psp4.13130)) - patients with radioiodine-refractory differentiated thyroid cancer.

- Add Aoki 2024 intra-target microdosing PBPK-PKRO model
  ([doi:10.3389/fphar.2024.1366160](https://doi.org/10.3389/fphar.2024.1366160))
  -- drug-agnostic simulated small-molecule compounds, no human subjects.

- Add Chen 2024 IL-6-mediated CYP3A suppression ([doi:10.1002/psp4.13073](https://doi.org/10.1002/psp4.13073)) -- adults with relapsed/refractory non-Hodgkin lymphoma receiving mosunetuzumab.

- Add Saleh 2023 LeiCNS-PK3.0 mouse CNS PBPK models for cyclophosphamide,
  quinidine, erlotinib, phenobarbital, colchicine, ribociclib, topotecan,
  cefadroxil, prexasertib and methotrexate
  ([doi:10.1007/s11095-023-03554-5](https://doi.org/10.1007/s11095-023-03554-5))
  -- laboratory mice (CD1 nude, NMRI, FVB, ICR, C57BL/6 Pept2+/+).

- Add Granda 2024 tenofovir, oseltamivir carboxylate and kynurenic acid ([doi:10.1111/cts.13678](https://doi.org/10.1111/cts.13678)) - adult outpatients spanning CKD stages 1-5.

- Add Wattanakul 2024 primaquine ([doi:10.1038/s41467-024-47908-y](https://doi.org/10.1038/s41467-024-47908-y)) -- lactating women with *Plasmodium vivax* infection and their breastfed infants.

- Add Pei 2023 tacrolimus PBPK and popPK models ([doi:10.3390/pharmaceutics15112580](https://doi.org/10.3390/pharmaceutics15112580)) -- adult heart transplant recipients.

- Add Wojciechowski 2023 ritlecitinib ([doi:10.1007/s40262-023-01318-3](https://doi.org/10.1007/s40262-023-01318-3)) - three model iterations covering healthy participants and patients with rheumatoid arthritis, ulcerative colitis, alopecia areata or vitiligo, plus moderate hepatic and severe renal impairment.

- Add Walsh 2024 buprenorphine ([doi:10.1038/s41386-023-01793-z](https://doi.org/10.1038/s41386-023-01793-z)) - non-treatment-seeking adults with moderate to severe opioid use disorder.

- Canonical unit spellings. The machine-readable `units` block wrote the same
  time unit three ways -- `"hour"` in 643 models, `"h"` in 208 and `"hr"` in 28
  -- so a consumer parsing `units$time` could not canonicalise it without
  carrying its own spelling table. Same for `"minute"`/`"min"` and
  `"microgram"`/`"ug"`.

  790 spellings in the `units` block and 204 unit hints in labels are
  normalised (`/hour` and `/hr` to `/h`, `pM.day` to `pM*day`, `mcmol` to
  `umol`). `conventions$timeUnitSpellings` and `$doseUnitSpellings` hold the
  map; `checkModelConventions()` errors on a non-canonical spelling and
  `buildModelDb()` aborts, so it cannot regrow. The extraction skill's template
  and checklist require it of new models.

  **Spelling is normalised; dimension is never converted.** `min` and `h` are
  both canonical and are never conflated -- rewriting one as the other would
  misstate every value. Generic dimensionless models (`PK_1cmt`, `PK_2cmt`, ...)
  keep their `"time_unit"` / `"dose_unit"` placeholders, and
  `Beal_2001_iv1cmt_bql` keeps time in half-lives; those are exempt by design.

  No model's numeric values changed.

  `kon` is deliberately **not** canonicalised: the prefix covers at least three
  different dimensionalities in this library (3D molar rates, QSP 2D on-rates
  carrying a length dimension, and mass-concentration forms), plus four
  parameters where `KON..` is the source paper's name for an EC50 or an Emax.
  The reasoning is recorded in `inst/references/parameter-names.md`.

- Add Foster 2023 enrofloxacin and ciprofloxacin ([doi:10.1111/jvim.16866](https://doi.org/10.1111/jvim.16866)) - client-owned cats with normal to severely reduced kidney function.

- Add Yang 2024 meropenem ([doi:10.1038/s41598-024-64223-0](https://doi.org/10.1038/s41598-024-64223-0)) - critically ill adult ICU patients with severe pneumonia.

- Add Schreib 2024 busulfan ([doi:10.3390/pharmaceutics16010013](https://doi.org/10.3390/pharmaceutics16010013)) - pediatric patients undergoing hematopoietic stem cell transplantation.

- Add Sharma 2023 nitrofurantoin whole-body PBPK ([doi:10.3390/pharmaceutics15092199](https://doi.org/10.3390/pharmaceutics15092199)) -- three models: rabbits, rats, and human adults.

- Canonical names for lactation transfer and tissue partitioning.
  `cmpr` / `lcmpr` is the estimated milk-to-plasma **concentration** ratio, for
  lactation popPK models with too few milk samples to support a milk compartment
  (`Cmilk <- cmpr * Cc`). Separately, the per-tissue partition-coefficient family
  `kp_<tissue>` / `lkp_<tissue>` -- de-facto used across six PBPK model files
  since before the register existed, but never registered -- is now documented,
  with `kp_milk` / `lkp_milk` added as its lactation member. The two are distinct
  concepts and coexist: `cmpr` is fitted to paired plasma and milk observations,
  while a `kp` is a predicted or literature-fixed partition constant. Both are in
  `inst/references/parameter-names.md`.

- Add Li 2023 ornidazole ([doi:10.3390/pharmaceutics15112524](https://doi.org/10.3390/pharmaceutics15112524)) - breastfeeding women after caesarean section, with colostrum concentrations.

- Add Adeojo 2024 levonorgestrel ([doi:10.3390/pharmaceutics16081050](https://doi.org/10.3390/pharmaceutics16081050)) -- adult women, pooled across four published trials with and without efavirenz.

- Add Gong 2023 pemigatinib ([doi:10.1002/psp4.13064](https://doi.org/10.1002/psp4.13064)) - healthy participants and patients with advanced solid tumors including cholangiocarcinoma.

- Add Qi 2024 vosoritide ([doi:10.1007/s40262-024-01371-6](https://doi.org/10.1007/s40262-024-01371-6)) - children with achondroplasia aged 0.95-15 years.

- Add Liang 2024 rituximab and anti-PLA2R titer models ([doi:10.3389/fphar.2024.1197651](https://doi.org/10.3389/fphar.2024.1197651)) - adults with primary membranous nephropathy.

- Add Shu 2024 posaconazole ([doi:10.1038/s41598-024-70955-w](https://doi.org/10.1038/s41598-024-70955-w)) — Chinese hematopoietic stem cell transplantation recipients on oral suspension.

- Add Mody 2023 doxorubicin + dexrazoxane ([doi:10.3389/fphar.2023.1239141](https://doi.org/10.3389/fphar.2023.1239141)) -- in vitro JIMT-1 and MDA-MB-468 human breast cancer cell lines, with clinical translation.

- Add Abdelgawad 2024 linezolid ([doi:10.1093/infdis/jiad413](https://doi.org/10.1093/infdis/jiad413)) -- adults with HIV-associated tuberculous meningitis, in plasma and cerebrospinal fluid.

- Add Jung 2024 clopidogrel ([doi:10.1002/psp4.13053](https://doi.org/10.1002/psp4.13053)) -- healthy Korean male adults stratified by CYP2C19 phenotype.

- Add Singu 2024 gentamicin ([doi:10.3390/children11080898](https://doi.org/10.3390/children11080898)) - Namibian neonates with suspected or confirmed sepsis.

- Add Goulooze 2022 finerenone UACR and eGFR dose-exposure-response models ([doi:10.1007/s40262-022-01124-3](https://doi.org/10.1007/s40262-022-01124-3)) -- adults with chronic kidney disease and type 2 diabetes (FIDELIO-DKD).

- Add Liang 2024 osimertinib ([doi:10.3389/fphar.2024.1363259](https://doi.org/10.3389/fphar.2024.1363259)) - simulated Caucasian, Japanese and Chinese NSCLC patients with EGFR T790M / L858R mutations.

- Add Dias 2024 quetiapine ([doi:10.1002/psp4.13107](https://doi.org/10.1002/psp4.13107)) - naive and schizophrenia phenotyped Wistar rats.

- Add Zazo 2024 ceftazidime-avibactam ([doi:10.3390/antibiotics13090861](https://doi.org/10.3390/antibiotics13090861)) — simulated non-ICU and ICU adults with renal impairment (two population-specific models).

- Canonical unit spellings. The machine-readable `units` block wrote the same
  time unit three ways -- `"hour"` in 643 models, `"h"` in 208 and `"hr"` in 28
  -- so a consumer parsing `units$time` could not canonicalise it without
  carrying its own spelling table. Same for `"minute"`/`"min"` and
  `"microgram"`/`"ug"`.

- Add Willemin 2024 interleukin-6 / CYP PBPK ([doi:10.1002/psp4.13144](https://doi.org/10.1002/psp4.13144)) - patients with relapsed/refractory multiple myeloma experiencing cytokine release syndrome after teclistamab.

- Add Shen 2024 vancomycin ([doi:10.1002/psp4.13151](https://doi.org/10.1002/psp4.13151)) - Southern Chinese children aged 1 month to 17 years on routine therapeutic drug monitoring.

- Add Marques 2024 salbutamol ([doi:10.3390/pharmaceutics16070881](https://doi.org/10.3390/pharmaceutics16070881)) - healthy adult volunteers given a single 600 ug dry-powder-inhaler dose.

- Add Wang 2024 saxagliptin and 5-hydroxy saxagliptin ([doi:10.1186/s40360-024-00757-3](https://doi.org/10.1186/s40360-024-00757-3)) - streptozotocin/high-fat-diet type 2 diabetic Sprague-Dawley rats.

- Add Lee 2024 eculizumab ([doi:10.1007/s00228-024-03703-8](https://doi.org/10.1007/s00228-024-03703-8)) - healthy adults and patients with paroxysmal nocturnal haemoglobinuria.

- Add Chen 2024 noscapine ([doi:10.1007/s40268-024-00466-6](https://doi.org/10.1007/s40268-024-00466-6)) - healthy adults genotyped for CYP2C9.

- Add Ait-Oudhia 2024 sotatercept ([doi:10.1002/psp4.13166](https://doi.org/10.1002/psp4.13166)) - healthy post-menopausal women and patients with pulmonary arterial hypertension.

- Add Khwarg 2024 donepezil ([doi:10.1007/s40120-024-00643-4](https://doi.org/10.1007/s40120-024-00643-4)) - healthy adult men given oral tablets or long-acting intramuscular GB-5001 injections.

- Add Kim 2024 meropenem ([doi:10.3390/antibiotics13090849](https://doi.org/10.3390/antibiotics13090849)) - healthy Korean adults with normal renal function.

- Add Kuroda 2023 cephalothin ([doi:10.1294/jes.34.111](https://doi.org/10.1294/jes.34.111)) - Thoroughbred horses given intramuscular and intravenous doses.

- Add Huppe 2023 fosfomycin ([doi:10.1038/s41598-023-45084-5](https://doi.org/10.1038/s41598-023-45084-5)) - critically ill adults with renal insufficiency during continuous venovenous hemodialysis.

- Add Zhang 2024 sertraline ([doi:10.1016/j.heliyon.2024.e25231](https://doi.org/10.1016/j.heliyon.2024.e25231)) - Chinese inpatients with psychiatric disorders, aged 11-79 years.

- Add Yates 2023 theoretical oncology dose-response models ([doi:10.1002/psp4.13020](https://doi.org/10.1002/psp4.13020)) - illustrative exponential, Mayneord and von Bertalanffy tumor-growth laws; no data fitted.

- Add Wu 2023 SPI-62 ([doi:10.1007/s40262-023-01278-8](https://doi.org/10.1007/s40262-023-01278-8)) - healthy adults.

- Add Wang 2024 amphenmulin ([doi:10.1128/spectrum.03675-23](https://doi.org/10.1128/spectrum.03675-23)) - broiler chickens and in-vitro *Mycoplasma gallisepticum*.

- Add Park 2023 mycophenolic acid ([doi:10.3390/pharmaceutics15122741](https://doi.org/10.3390/pharmaceutics15122741)) - paediatric haematopoietic stem cell transplant recipients.

- Disambiguated the overloaded `OC` name, which denoted five unrelated
  concepts. Osteocalcin `OC` -> `OSTCALC` (uppercase, matching the sibling
  biomarkers `P1NP` / `PSA` / `PLT` / `WBC` and the sister model
  `Shoji_2017_fosdagrocorat_p1np`); the oseltamivir-carboxylate metabolite
  suffix `_oc` -> `_oselcarb` (`central_oselcarb`, `lcl_oselcarb`,
  `Cc_oselcarb`, ...) across Chairat 2016, Kamal 2013 and Standing 2012;
  Hussein 1997's unregistered `OC` covariate column -> the existing
  canonical `CONMED_BIRTHCONTROL`. **Breaking for simulation code** that
  references the old names. Model ids and vignette filenames are unchanged.
- Merged the `CONMED_DIUR` covariate canonical into `CONMED_DIURETIC`. The two
  names denoted the same concept (concomitant diuretic use); the split was an
  artifact of independent extractions. `Wright_2016_allopurinol`,
  `Wright_2013_allopurinol` and `Stocker_2012_oxypurinol` now use
  `CONMED_DIURETIC` (and `e_conmed_diuretic_*` effect parameters). **Breaking for
  simulation code**: event tables / `keep=` vectors referencing `CONMED_DIUR`
  must be renamed. Per-paper diuretic class composition (which differs between
  these models) remains documented in each model's `covariateData` notes.
* Add Wada 2023 sparsentan ([doi:10.1002/psp4.12996](https://doi.org/10.1002/psp4.12996)) -- healthy volunteers, subjects with hepatic impairment, and patients with primary or genetic focal segmental glomerulosclerosis (ratifies new `CONMED_CYP3A4_INH_MOD` and `FORM_CRUSHED_TABLET` covariate canonicals).
* Add Beal 2001 one-compartment IV-bolus BQL methodology template ([doi:10.1023/a:1012299115260](https://doi.org/10.1023/a:1012299115260)) -- methodology reference (no drug, no patients); packages the SI1 generative model from Beal's M1-M7 below-quantification-limit paper as a teaching template with CL = 0.693 and Vd = 1 (time in half-lives).
* Add Luu 2017 nusinersen ([doi:10.1002/jcph.884](https://doi.org/10.1002/jcph.884)) -- pediatric patients with spinal muscular atrophy receiving intrathecal nusinersen.
* Add Gaohua 2012 pregnancy PBPK ([doi:10.1111/j.1365-2125.2012.04363.x](https://doi.org/10.1111/j.1365-2125.2012.04363.x)) -- healthy pregnant Caucasian women (14-compartment whole-body p-PBPK with GA-dependent maternal physiology, applied to caffeine [CYP1A2], metoprolol [CYP2D6], and midazolam [CYP3A4]); ratifies new canonical bare `skin` PBPK compartment.
* Add Morris 2011 telapristone ([doi:10.1208/s12248-011-9304-7](https://doi.org/10.1208/s12248-011-9304-7)) -- adult women with single-dose telapristone-acetate phase I/II PK studies; ratifies new `cdb4453` metabolite suffix, `MIX_FAST_ELIM`, and `RENALIMP_MOD` covariate canonicals.
* Add Padoin 1998 cephalexin ([doi:10.1128/aac.42.6.1463](https://doi.org/10.1128/aac.42.6.1463)) -- preclinical male Wistar rats with a cephalexin / quinapril oral coadministration DDI (ratifies new `CONMED_QPRL_ORAL` covariate canonical).
* Add van Rongen 2018 metformin ([doi:10.1007/s40272-018-0293-1](https://doi.org/10.1007/s40272-018-0293-1)) -- overweight and obese Caucasian adolescents.
* Add Tarning 2012 artemether and dihydroartemisinin ([doi:10.1186/1475-2875-11-293](https://doi.org/10.1186/1475-2875-11-293)) -- pregnant women with uncomplicated Plasmodium falciparum malaria in Uganda (joint parent + DHA popPK with zero-order dissolution and 6-compartment transit absorption).
* Add Lowe 2009 omalizumab ([doi:10.1111/j.1365-2125.2009.03401.x](https://doi.org/10.1111/j.1365-2125.2009.03401.x)) -- adults and adolescents (12-79 years) with severe persistent allergic asthma plus healthy atopic volunteers (mechanism-based binding popPK/PD for free omalizumab, free IgE, and the omalizumab-IgE complex; extends Hayashi 2007).
* Add Grover 2011 tacrolimus ([doi:10.1124/dmd.111.041350](https://doi.org/10.1124/dmd.111.041350)) -- adult Native American renal transplant recipients on stable twice-daily oral tacrolimus.
* Add Schoemaker 2017 brivaracetam ([doi:10.1007/s00228-017-2230-6](https://doi.org/10.1007/s00228-017-2230-6)) -- paediatric patients with epilepsy aged 1 month to 16 years; ratifies new `CONMED_PB`, `CONMED_CBZ`, and `CONMED_VPA` covariate canonicals.
* Add Nanga 2019 tacrolimus meta-model ([doi:10.1111/bcp.14110](https://doi.org/10.1111/bcp.14110)) -- pooled paediatric and adult solid-organ transplant recipients (n = 281; 201 liver and 80 kidney) on oral tacrolimus; ratifies new `TX_LIVER` and `FORM_SYRUP` covariate canonicals.
* Rewrite Fiedler-Kelly 2019 fremanezumab ([doi:10.1111/bcp.14096](https://doi.org/10.1111/bcp.14096)) to support both IV and SC routes -- the packaged model now carries the route-specific central volume (Vc,IV = 2.98 L FIXED vs Vc,SC = 1.88 L) and the route-specific residual-error structure (IV proportional-only vs SC additive+proportional) from Fiedler-Kelly 2019 Table 2, switched by the `ROUTE_IV` covariate column. Removes the previous spurious allometric weight effect on Vp (the source paper holds Vp FIXED with no weight effect).
* Add Yuan 2019 concizumab ([doi:10.1016/j.ejps.2019.105032](https://doi.org/10.1016/j.ejps.2019.105032)) -- healthy adult males.
* `_pkgdown.yml` is now refreshed in place using `# AUTOGEN:<name>:BEGIN`/`:END` signpost comments rather than a full `yaml::read_yaml` + `yaml::write_yaml` round-trip. The signpost-managed regions are the `specific_drugs` and `ddmore` navbar menus and the top-level `articles:` index; everything else in the file (template, navbar.structure, maintainer-added keys, comments, formatting) is preserved byte-for-byte. The modeldb gains two new cached columns -- `label` (human-readable navbar text like `Ustekinumab (Aguiar 2021)` or `DDMoRe: lidocaine`) and `category` (`specificDrugs` / `ddmore` / `other`) -- so the regeneration costs no extra work beyond the existing per-file md5 cache. Fixes the recurring pkgdown CI breakage where stale entries crept into the `articles:` index when the YAML round-trip captured a transient filesystem state.
* Add Quartino 2016 trastuzumab ([doi:10.1007/s00280-015-2922-5](https://doi.org/10.1007/s00280-015-2922-5)) -- women with HER2-positive early breast cancer (HannaH study; subcutaneous fixed 600 mg q3w and intravenous weight-based regimens).
* Add Diepstraten 2013 propofol ([doi:10.1038/psp.2013.47](https://doi.org/10.1038/psp.2013.47)) -- morbidly obese and nonobese adults, adolescents, and children (pooled meta-analysis across five previously published studies).
* List of models vignette is grouped into General / Specific Drugs / DDMoRe with capped-width tables and uses the navbar-style "Drug (Author Year)" name (`DDMoRe: <drug>` for the parameterless DDMoRe entries). Drug-specific articles in `vignettes/articles/` are now hidden from the auto-generated pkgdown Articles index and reachable only via the navbar dropdowns or direct URL. Validation vignettes that previously titled themselves with the machine basename now use the human form so the page title matches the link.
* Add Diep 2026 donidalorsen ([doi:10.1002/psp4.70206](https://doi.org/10.1002/psp4.70206)) -- pooled healthy volunteers and adult / adolescent patients with hereditary angioedema (ratifies new `DIS_HAE` covariate canonical).
* Add Attarwala 2023 mRNA-3927 ([doi:10.1089/nat.2022.0036](https://doi.org/10.1089/nat.2022.0036)) -- preclinical PCC-deficient mice, juvenile Sprague-Dawley rats, and cynomolgus monkeys (translational semi-mechanistic PK + PK/PD model for an LNP-encapsulated dual mRNA encoding PCCA / PCCB subunits of propionyl-CoA carboxylase; allometrically scalable to humans).
* Add Dunn 2025 tranexamic acid ([doi:10.1002/jcph.70031](https://doi.org/10.1002/jcph.70031)) -- pregnant individuals receiving IV, IM, or oral TXA for prevention or treatment of postpartum hemorrhage.
* Add Parhiz 2024 mRNA-LNP ([doi:10.1016/j.omtn.2024.102175](https://doi.org/10.1016/j.omtn.2024.102175)) -- preclinical mice (C57BL/6, 25 g); whole-body PBPK plus luciferase-expression model for systemically administered firefly-luciferase mRNA in lipid nanoparticles (bare-LNP homogenate-assay parameter set).
* Add Schoenmakers 2025 betamethasone ([doi:10.1002/bcp.70035](https://doi.org/10.1002/bcp.70035)) -- pregnant women with imminent preterm birth, including early-onset pre-eclampsia (ratifies new `DIS_EOPE` covariate canonical for early-onset pre-eclampsia indicator).
* Add Rambiritch 2016 glibenclamide ([doi:10.2147/CPAA.S102676](https://doi.org/10.2147/CPAA.S102676)) -- poorly controlled South African adults with type 2 diabetes mellitus.
* Add Lahu 2010 roflumilast ([doi:10.2165/11536600-000000000-00000](https://doi.org/10.2165/11536600-000000000-00000)) -- adult healthy volunteers and patients with moderate-to-severe COPD (joint parent + roflumilast N-oxide popPK; ratifies new `noxide` registered metabolite suffix and new `DIS_COPD` covariate canonical).
* Add Brill 2014 midazolam ([doi:10.1007/s40262-014-0166-x](https://doi.org/10.1007/s40262-014-0166-x)) -- 20 morbidly obese surgical patients and 12 non-obese healthy volunteers receiving semi-simultaneous oral and intravenous midazolam.
* Add Diao 2014 rFIXFc ([doi:10.1007/s40262-013-0129-7](https://doi.org/10.1007/s40262-013-0129-7)) -- adolescents and adults with severe to moderate haemophilia B (three-compartment IV popPK with estimated body-weight power exponents on CL and V1).
* Add Plock 2014 ferumoxytol ([doi:10.1007/s40262-014-0203-9](https://doi.org/10.1007/s40262-014-0203-9)) -- pooled adult healthy volunteers and CKD-stage-5D haemodialysis patients (two-compartment Michaelis-Menten IV PK with WT and sex covariates on V1; haemodialysis time-varying V1 mechanism documented in the vignette but not enabled in the model file).
* Add Qi 2014 sapropterin ([doi:10.1007/s40262-014-0196-4](https://doi.org/10.1007/s40262-014-0196-4)) -- pediatric and adult patients with phenylketonuria pooled from BioMarin studies PKU-015 (0-6 years) and PKU-004 (9-50 years).
* Add Retlich 2015 linagliptin ([doi:10.1007/s40262-014-0232-4](https://doi.org/10.1007/s40262-014-0232-4)) -- adults with type 2 diabetes mellitus (coupled two-compartment popPK with concentration-dependent DPP-4 binding in both compartments plus a sigmoid Emax PD model relating linagliptin to plasma DPP-4 activity; ratifies new `GGT`, `FPG`, `DPP4_BL_RFU`, `CONMED_METFORMIN`, and `FORM_LINAG_TAB1` covariate canonicals).
* Add Buil-Bruna 2015 lanreotide ([doi:10.1007/s40262-015-0329-4](https://doi.org/10.1007/s40262-015-0329-4)) -- adults with gastroenteropancreatic neuroendocrine tumors.
* Add van Rongen 2016 acetaminophen ([doi:10.1007/s40262-015-0357-0](https://doi.org/10.1007/s40262-015-0357-0)) -- morbidly obese and non-obese adults receiving 2 g IV acetaminophen (parent + glucuronide + sulphate + CYP2E1-oxidation metabolite popPK; ratifies new `cysmer` registered metabolite suffix for combined acetaminophen cysteine + mercapturate).
* Add Gupta 2017 ixazomib ([doi:10.1007/s40262-017-0526-4](https://doi.org/10.1007/s40262-017-0526-4)) -- adult patients with multiple myeloma, lymphoma, solid tumours, or light-chain amyloidosis (pooled phase I-III TOURMALINE-MM1 popPK; three-compartment oral + IV PK with BSA on V4 and time-after-dose-varying residual SD).
* Add Savic 2017 cladribine ([doi:10.1007/s40262-017-0516-6](https://doi.org/10.1007/s40262-017-0516-6)) -- adults with relapsing-remitting multiple sclerosis (three-compartment popPK with renal-clearance scaling on Cockcroft-Gault CRCL and a non-renal-clearance effect of concomitant subcutaneous IFN beta-1a; ratifies new `CONMED_IFNB1A` covariate canonical).
* Add Hendriksen 2013 artesunate ([doi:10.1038/clpt.2013.26](https://doi.org/10.1038/clpt.2013.26)) -- Tanzanian children (7 months to 11 years) with severe Plasmodium falciparum malaria receiving intramuscular artesunate.
* Add Ait-Oudhia 2012 canakinumab ([doi:10.1038/psp.2012.6](https://doi.org/10.1038/psp.2012.6)) -- adults with rheumatoid arthritis (integrated PK / quasi-equilibrium IL-1beta binding / CRP transduction / ACR latent-variable model with logit-transform mapping to ACR20/50/70 response probabilities).
* Add Wu 2012 sirolimus ([doi:10.1038/psp.2012.18](https://doi.org/10.1038/psp.2012.18)) -- adults with advanced solid tumors receiving 1-60 mg/week oral sirolimus (first popPK report describing the nonlinear absorption with a two-compartment + Michaelis-Menten saturable absorption structure; hematocrit power covariate on apparent oral clearance).
* Add de Hoogd 2017 morphine ([doi:10.1007/s40262-017-0544-2](https://doi.org/10.1007/s40262-017-0544-2)) -- 20 morbidly obese adults (post gastric-bypass) and 20 healthy adult volunteers (joint parent + M3G + M6G popPK).
* Add Doldan-Martelli 2013 EGF-IFN chimeric ligand ([doi:10.1038/psp.2013.2](https://doi.org/10.1038/psp.2013.2)) -- in-vitro Daudi / Daudi-EGFR cells (mechanistic two-receptor binding model for an EGF-IFNalpha-2a chimera with selectivity-by-affinity-mutant analysis).
* Add Tetschke 2018 erythropoiesis ([doi:10.3390/pr6090157](https://doi.org/10.3390/pr6090157)) -- healthy adult male volunteers undergoing a single 1-unit blood donation (three-compartment endogenous mixed-effects model for red blood cell regeneration after phlebotomy; ratifies new `THB_MASS` covariate canonical for total hemoglobin mass measured by CO-rebreathing).
* Add Ekhart 2008 carboplatin ([doi:10.1007/s00280-008-0856-x](https://doi.org/10.1007/s00280-008-0856-x)) -- adult cancer patients (underweight to obese) receiving combination chemotherapy.
* Add Aregbe 2012 alvespimycin / 17-DMAG ([doi:10.1007/s00280-012-1859-1](https://doi.org/10.1007/s00280-012-1859-1)) -- adult patients with advanced solid tumors.
* Add Ferron 2013 cabazitaxel ([doi:10.1007/s00280-012-2058-9](https://doi.org/10.1007/s00280-012-2058-9)) -- patients with advanced solid tumors (pooled five Phase I-III studies including TROPIC mCRPC).
* Add Panoilia 2015 bevacizumab ([doi:10.1007/s00280-015-2701-3](https://doi.org/10.1007/s00280-015-2701-3)) -- adults with stage IV metastatic colorectal cancer (TMDD QSS binding model for bevacizumab and free VEGF165 with fixed allometric body-weight scaling).
* Add Varatharajan 2016 daunorubicin ([doi:10.1007/s00280-016-3166-8](https://doi.org/10.1007/s00280-016-3166-8)) -- adult de novo acute myeloid leukaemia patients (joint daunorubicin + daunorubicinol popPK with two-compartment disposition for each component; ratifies new `dol` registered metabolite suffix).
* Add Ozawa 2007 docetaxel ([doi:10.1111/j.1349-7006.2007.00615.x](https://doi.org/10.1111/j.1349-7006.2007.00615.x)) -- Japanese adult cancer patients (62 subjects with breast, NSCLC, head-and-neck, or other solid tumours; coupled three-compartment IV PK and modified Friberg myelosuppression PD with a dexamethasone-induced ANC-bump compartment; AGP power-form effect on the drug-effect slope).
* Add Alsultan 2017 pyrazinamide ([doi:10.1128/AAC.02625-16](https://doi.org/10.1128/AAC.02625-16)) -- adults with drug-susceptible pulmonary tuberculosis (TBTC studies 27 and 28).
* Add Yukawa 1990 phenytoin ([doi:10.1248/cpb.38.1973](https://doi.org/10.1248/cpb.38.1973)) -- Japanese epileptic outpatients on chronic oral phenytoin maintenance therapy (ratifies `CONMED_AED`, `FORM_POWDER`, and `DOSE_PHT_MGKGD` covariate canonicals).
* Add Vezina 2010 valganciclovir ([doi:10.2147/CPAA.S8341](https://doi.org/10.2147/CPAA.S8341)) -- pediatric solid organ transplant recipients receiving valganciclovir prophylaxis for Epstein-Barr virus disease.
* Add Kunisawa 2014 olprinone ([doi:10.2147/CPAA.S50626](https://doi.org/10.2147/CPAA.S50626)) -- healthy adult Japanese male volunteers.
* Add Wang 2021 pertuzumab ([doi:10.1007/s00280-021-04296-0](https://doi.org/10.1007/s00280-021-04296-0)) -- women with HER2-positive early breast cancer (FeDeriCa study, NCT03493854; pooled IV and subcutaneous fixed-dose-combination pertuzumab+trastuzumab cohorts).
* Add Storset 2014 tacrolimus ([doi:10.1111/bcp.12361](https://doi.org/10.1111/bcp.12361)) -- adult kidney-transplant recipients (theory-based plasma-parameterised popPK with FFM allometry, CYP3A5-expresser effects, prednisolone-induced reduction in F, first-day-post-transplant F spike, and saturable haematocrit-dependent RBC binding linking plasma to whole-blood concentration).
* Add Overgaard 2019 semaglutide ([doi:10.1007/s13300-019-0581-y](https://doi.org/10.1007/s13300-019-0581-y)) -- pooled clinical pharmacology cohort of healthy volunteers and adults with type 2 diabetes.
* Add Rovei 1982 theophylline ([doi:10.1111/j.1365-2125.1982.tb02035.x](https://doi.org/10.1111/j.1365-2125.1982.tb02035.x)) -- healthy adult volunteers receiving single oral 125-500 mg theophylline tablets.
* Add Dunlap 2025 tacrolimus ([doi:10.1007/s40262-025-01529-w](https://doi.org/10.1007/s40262-025-01529-w)) -- adult allogeneic hematopoietic cell transplant recipients (CYP3A5 metabolizer phenotype and reduced-intensity conditioning effects on apparent clearance; ratifies new `HCT_COND_RIC` covariate canonical).
* Add Benkali 2010 tacrolimus (Clin Pharmacokinet 2010;49(10):683-92; DOI not on disk) -- stable adult renal transplant recipients switched to once-daily extended-release tacrolimus (Advagraf).
* Add Boer-Perez 2026 piperacillin ([doi:10.1128/aac.00998-25](https://doi.org/10.1128/aac.00998-25)) -- preterm and term neonates with severe infections.
* Add Jonckheere 2019 cefepime ([doi:10.1128/AAC.01552-19](https://doi.org/10.1128/AAC.01552-19)) -- critically ill ICU adults on continuous-infusion cefepime via target-controlled infusion.
* Add Wang 2015 rucaparib ([doi:10.1002/cpdd.176](https://doi.org/10.1002/cpdd.176)) -- adults with advanced solid tumors (Phase 1 first-in-patient study A4991002).
* Add Brown 2017 osimertinib ([doi:10.1111/bcp.13223](https://doi.org/10.1111/bcp.13223)) -- adults with advanced EGFR-mutation-positive non-small cell lung cancer pooled with healthy volunteers (joint osimertinib + AZ5104 metabolite popPK).
* Add Caldes 2009 ganciclovir ([doi:10.1128/AAC.00085-09](https://doi.org/10.1128/AAC.00085-09)) -- adult solid organ transplant recipients (kidney, liver, heart) with cytomegalovirus infection receiving IV ganciclovir followed by oral valganciclovir.
* Add Vezina 2014 valganciclovir ([doi:10.1111/bcp.12343](https://doi.org/10.1111/bcp.12343)) -- paediatric and adult solid organ transplant recipients receiving valganciclovir prophylaxis.
* Add Chen 2021 ganciclovir ([doi:10.1002/jcph.1735](https://doi.org/10.1002/jcph.1735)) -- adult Chinese renal allograft recipients receiving oral valganciclovir for CMV prophylaxis.
* Add Guo 2022 PF-06939999 ([doi:10.1002/psp4.12882](https://doi.org/10.1002/psp4.12882)) -- adults with advanced solid tumors (Phase I FIP study NCT03854227; coupled PK + plasma SDMA indirect-response + Friberg semi-mechanistic platelet model).
* Add Urien 2005 capecitabine ([doi:10.1007/s10928-005-0018-2](https://doi.org/10.1007/s10928-005-0018-2)) -- adults with metastatic cancer.
* Add Garg 2014 pertuzumab ([doi:10.1007/s00280-014-2560-3](https://doi.org/10.1007/s00280-014-2560-3)) -- patients with a variety of HER2-targeted solid tumors.
* Add Chien 2022 imatinib ([doi:10.1007/s00280-022-04454-y](https://doi.org/10.1007/s00280-022-04454-y)) -- healthy Caucasian volunteers receiving a single oral 400 mg dose.
* Add Akbar 2025 voriconazole ([doi:10.1371/journal.pone.0318883](https://doi.org/10.1371/journal.pone.0318883)) -- adult and pediatric Pakistani cancer patients receiving therapeutic drug monitoring of intravenous voriconazole.
* Add Ayyar 2024 givosiran ([doi:10.1016/j.xphs.2023.10.026](https://doi.org/10.1016/j.xphs.2023.10.026)) -- adults with acute hepatic porphyria (mechanistic translational PK with ASGPR-mediated hepatocyte uptake, endolysosomal kinetics, and cytoplasmic RISC loading; ratifies new `asn1` registered metabolite suffix for AS(N-1)3' GalNAc-siRNA).
* Add Salinger 2013 magnesium sulphate ([doi:10.1111/1471-0528.12222](https://doi.org/10.1111/1471-0528.12222)) -- pregnant women with pre-eclampsia receiving IV or IM MgSO4-7H2O for eclampsia prevention.
* Add Chakraborty 2012 canakinumab ([doi:10.2165/11599820-000000000-00000](https://doi.org/10.2165/11599820-000000000-00000)) -- adults and children with cryopyrin-associated periodic syndromes (CAPS) and four other inflammatory cohorts (population PK / IL-1b binding model with quasi-steady-state binding; ratifies new `il1b` registered metabolite suffix).
* Add Rodrigues 2017 oxcarbazepine ([doi:10.1111/bcp.13392](https://doi.org/10.1111/bcp.13392)) -- epileptic children aged 2-12 years (parent oxcarbazepine + active monohydroxy-derivative metabolite popPK with allometric weight scaling, MHD-to-OXC back-transformation, and concomitant enzyme-inducing antiepileptic drug effect; ratifies new `CONMED_EIAED` covariate canonical and `mhd` registered metabolite suffix).
* Add Li 2017 cediranib ([doi:10.1111/bcp.13266](https://doi.org/10.1111/bcp.13266)) -- adult cancer patients on oral cediranib (AZD2171).
* Add Tortorici 2017 alpha-1 proteinase inhibitor ([doi:10.1111/bcp.13358](https://doi.org/10.1111/bcp.13358)) -- adults with severe alpha-1 antitrypsin deficiency.
* Add Valitalo 2017 ketorolac ([doi:10.1111/bcp.13311](https://doi.org/10.1111/bcp.13311)) -- women at delivery, postpartum women, nonpregnant women, and men receiving single IV racemic ketorolac (joint R- / S-enantiomer popPK).
* Add Talke 2018 dexmedetomidine ([doi:10.1111/bcp.13571](https://doi.org/10.1111/bcp.13571)) -- healthy adult volunteers (3-cmt IV PK + effect-compartment sigmoid Emax PD for vasoconstriction).
* Add Huynh 2026 VRC07-523LS ([doi:10.1093/jac/dkaf449](https://doi.org/10.1093/jac/dkaf449)) -- HIV-exposed infants and healthy adults.
* Add Pu 2021 evinacumab ([doi:10.1002/psp4.12711](https://doi.org/10.1002/psp4.12711)) -- healthy volunteers and adult / pediatric patients with homozygous familial hypercholesterolemia.
* Add Zhang 2021 dupilumab ([doi:10.1002/psp4.12667](https://doi.org/10.1002/psp4.12667)) -- adult and adolescent patients with moderate-to-severe asthma (pooled with healthy adults).
* Add Hennig 2015 rifabutin ([doi:10.1128/AAC.01195-15](https://doi.org/10.1128/AAC.01195-15)) -- HIV-infected adults with pulmonary tuberculosis (joint parent + 25-O-desacetyl rifabutin metabolite popPK with allometric weight scaling, sex on V/F, and SLCO1B1 rs11045819 on bioavailability).
* Add Hennig 2013 tobramycin ([doi:10.1007/s40262-013-0036-y](https://doi.org/10.1007/s40262-013-0036-y)) -- adults and children with and without cystic fibrosis.
* Add Yassen 2025 asundexian ([doi:10.1002/psp4.70142](https://doi.org/10.1002/psp4.70142)) -- healthy volunteers and adults at risk for thromboembolic / cardiovascular events (PACIFIC-AF, PACIFIC-STROKE, PACIFIC-AMI Phase II + six Phase I studies; pooled n = 2914).
* Add Jelliffe 2014 digoxin ([doi:10.1097/FTD.0000000000000023](https://doi.org/10.1097/FTD.0000000000000023)) -- adults requiring digoxin therapy across normal-to-anephric renal function.
* Add Aksenov 2018 uricAcid ([doi:10.14814/phy2.13614](https://doi.org/10.14814/phy2.13614)) -- healthy adults and gout patients with hyperuricemia.
* Add Garmann 2017 BAY 81-8973 ([doi:10.1111/hae.13192](https://doi.org/10.1111/hae.13192)) -- patients aged 1-61 years with severe haemophilia A (LEOPOLD I/II/Kids).
* Add Shah 2012 mAb PBPK ([doi:10.1007/s10928-011-9232-2](https://doi.org/10.1007/s10928-011-9232-2)) -- platform PBPK model parameterized for the human (71 kg) reference subject; 15 anatomical tissues with FcRn-mediated recycling.
* Add Kunarajah 2017 doxorubicin ([doi:10.1007/s00280-017-3309-6](https://doi.org/10.1007/s00280-017-3309-6)) -- paediatric oncology patients (joint doxorubicin / doxorubicinol popPK with cardiac troponin I turnover sub-model; ratifies new `PRIOR_ANTHRACYCLINE_DOSE` covariate canonical and `doxol` registered metabolite suffix).
* Add Llanos-Paez 2017 gentamicin ([doi:10.1208/s12248-017-0173-6](https://doi.org/10.1208/s12248-017-0173-6)) -- pediatric oncology patients with febrile or fever-only neutropenia.
* Add Hennig 2015 phenytoin ([doi:10.1002/jcph.417](https://doi.org/10.1002/jcph.417)) -- critically ill children in a paediatric ICU (joint protein-unbound and protein-bound popPK with linear albumin partition).
* Add Bergmann 2014 tacrolimus ([doi:10.1097/FTD.0b013e31829f1ab8](https://doi.org/10.1097/FTD.0b013e31829f1ab8)) -- adult kidney transplant recipients.
* Add Harun 2019 cysticFibrosis ([doi:10.1136/thoraxjnl-2018-211550](https://doi.org/10.1136/thoraxjnl-2018-211550)) -- children with classic cystic fibrosis (sigmoid-Emax disease-progression model of FEV1% predicted from age 5 to 14 years).
* Add Hennig 2007 itraconazole ([doi:10.1111/j.1365-2125.2006.02778.x](https://doi.org/10.1111/j.1365-2125.2006.02778.x)) -- adult cystic fibrosis patients (D-optimal-designed cross-over single-dose 200 mg oral itraconazole with parallel parent / hydroxy-itraconazole metabolite popPK and per-formulation absorption / bioavailability).
* Add Hennig 2008 tobramycin ([doi:10.1111/j.1365-2125.2007.03045.x](https://doi.org/10.1111/j.1365-2125.2007.03045.x)) -- paediatric cystic fibrosis patients on once-daily intravenous tobramycin.
* Add Kirubakaran 2022 tacrolimus ([doi:10.1111/bcp.15566](https://doi.org/10.1111/bcp.15566)) -- adult heart transplant recipients with concomitant azole antifungal therapy.
* Add Hennig 2006 itraconazole ([doi:10.2165/00003088-200645110-00004](https://doi.org/10.2165/00003088-200645110-00004)) -- paediatric cystic-fibrosis and bone-marrow-transplant patients (parent + active metabolite popPK for oral itraconazole and hydroxy-itraconazole).
* Add Lawson 2022 busulfan ([doi:10.1002/psp4.12809](https://doi.org/10.1002/psp4.12809)) -- pediatric hematopoietic stem cell transplant recipients receiving once-daily IV busulfan.
* Add Archary 2019 abacavir ([doi:10.1111/bcp.13998](https://doi.org/10.1111/bcp.13998)) -- severely malnourished HIV-infected children.
* Add Bista 2015 fentanyl (manuscript, journal/DOI not on disk) -- adults with advanced cancer receiving Durogesic transdermal fentanyl matrix patches.
* Add Archary 2018 lopinavir ([doi:10.1097/INF.0000000000001867](https://doi.org/10.1097/INF.0000000000001867)) -- severely malnourished HIV-infected children (1-month to 12-year-olds) on twice-daily oral LPV/rtv with FFM allometric scaling and a total-cholesterol covariate effect on apparent clearance.
* Add Xu 2023 sabatolimab MBG453 ([doi:10.1002/psp4.12962](https://doi.org/10.1002/psp4.12962)) -- adults with advanced solid tumors or hematologic malignancies (AML, MDS, CMML).
* Add Goel 2016 sonidegib ([doi:10.1007/s00280-016-2982-1](https://doi.org/10.1007/s00280-016-2982-1)) -- healthy subjects and adults with advanced solid tumors or basal cell carcinoma.
* Add Tikiso 2021 abacavir ([doi:10.1111/bcp.14984](https://doi.org/10.1111/bcp.14984)) -- HIV-infected African children on abacavir-containing antiretroviral therapy.
* Add Stein 2019 tisagenlecleucel ([doi:10.1002/psp4.12388](https://doi.org/10.1002/psp4.12388)) -- pediatric and young adult patients with relapsed or refractory B-cell ALL (CAR-T cellular kinetic model: exponential expansion to Tmax followed by biexponential effector / memory-cell decline).
* Add van Wijk 2019 paracetamol ([doi:10.1038/s41598-019-38530-w](https://doi.org/10.1038/s41598-019-38530-w)) [DDMODEL00000294] -- preclinical zebrafish (Danio rerio) larvae 3-5 dpf under continuous 1 mM paracetamol bath exposure.
* Add Llanos-Paez 2017 gentamicin ([doi:10.1128/AAC.00205-17](https://doi.org/10.1128/AAC.00205-17)) -- pediatric oncology patients with febrile neutropenia.
* Add Llanos-Paez 2020 gentamicin ([doi:10.1128/AAC.01730-19](https://doi.org/10.1128/AAC.01730-19)) -- pediatric oncology and nononcology patients.
* Add Delattre 2010 amikacin ([doi:10.1097/FTD.0b013e3181f675c2](https://doi.org/10.1097/FTD.0b013e3181f675c2)) -- critically ill adults with severe sepsis or septic shock.
* Add Diep 2022 eplontersen ([doi:10.1111/bcp.15468](https://doi.org/10.1111/bcp.15468)) -- healthy volunteers in two phase 1 studies (two-compartment popPK with site-specific SC absorption and indirect-response PD on serum transthyretin).
* Add Laffont 2024 nalmefene ([doi:10.3389/fpsyt.2024.1399803](https://doi.org/10.3389/fpsyt.2024.1399803)) -- healthy adult volunteers receiving intranasal nalmefene HCl.
* Add Laffont 2024 naloxone ([doi:10.3389/fpsyt.2024.1399803](https://doi.org/10.3389/fpsyt.2024.1399803)) -- healthy adult volunteers receiving intranasal naloxone HCl.
* Add Cao 2013 mPBPK 12-mAb cohort ([doi:10.1007/s10928-013-9332-2](https://doi.org/10.1007/s10928-013-9332-2)) -- second-generation minimal PBPK structural model packaged as one entry per fit: 10 human-mAb entries (adecatumumab, mepolizumab, gevokizumab, GNbAC1, MEDI528, tefibazumab, PAmAb, PRO95780, siltuximab, visilizumab) under specificDrugs/, plus 2 preclinical mouse entries (mab7E3, mab8C2) under pharmacokinetics/.
* Add Sadouki 2025 meropenem / gentamicin / ciprofloxacin ([doi:10.1038/s41598-025-29354-y](https://doi.org/10.1038/s41598-025-29354-y)) -- in-vitro pharmacodynamic time-kill model for two- and three-way antibiotic combinations against Escherichia coli NCTC 12,241.
* Add Themans 2019 meropenem ([doi:10.1007/s40268-019-0268-x](https://doi.org/10.1007/s40268-019-0268-x)) [DDMODEL00000301] -- adults with severe pneumonia.
* Add Kovalenko 2016 dupilumab ([doi:10.1002/psp4.12136](https://doi.org/10.1002/psp4.12136)) [DDMODEL00000273] -- DDMORE-bundle replicate of the existing Kovalenko_2016_dupilumab specificDrugs entry; encodes the bundle's non-standard V2 weight covariate and `Output_simulated_*.lst` final estimates (no `Output_real_*.lst` shipped).
* Add Terranova 2018 paclitaxel ([doi:10.1016/j.jtbi.2018.04.012](https://doi.org/10.1016/j.jtbi.2018.04.012)) [DDMODEL00000274] -- preclinical xenograft mice (Dynamic Energy Budget tumor-growth-inhibition model coupling paclitaxel PK to tumor mass and host body weight / cachexia).
* Add Bajaj 2017 nivolumab ([doi:10.1002/psp4.12143](https://doi.org/10.1002/psp4.12143)) [DDMODEL00000284] -- adults with advanced solid tumors (DDMORE-source replicate of the paper-source Bajaj 2017 nivolumab popPK; time kept in hours to mirror the bundle).
* Add Birgersson 2019 artesunate ([doi:10.12688/wellcomeopenres.14849.2](https://doi.org/10.12688/wellcomeopenres.14849.2)) [DDMODEL00000297] -- pregnant and non-pregnant women with uncomplicated Plasmodium falciparum malaria in Burkina Faso (joint parent-metabolite popPK of artesunate and dihydroartemisinin with 3-compartment transit absorption and pregnancy / ALT / log-parasitaemia covariates).
* Add Dao 2020 sultiame ([doi:10.1002/prp2.558](https://doi.org/10.1002/prp2.558)) [DDMODEL00000298] -- healthy adult volunteers (4-compartment popPK with saturable plasma <-> erythrocyte binding plus cumulative urinary excretion).
* Add Conrado 2014 alzheimer ([doi:10.1007/s10928-014-9375-z](https://doi.org/10.1007/s10928-014-9375-z)) [DDMODEL00000290] -- adults with Alzheimer's disease (CAMD ADAS-Cog database; Richards three-parameter logistic disease-progression model with beta-regression residual).
* Add Schoemaker 2018 levetiracetam ([doi:10.1007/s40262-017-0597-2](https://doi.org/10.1007/s40262-017-0597-2)) [DDMODEL00000239] -- adults and children (4-16 years) with focal seizures (negative-binomial seizure-count PD model with mixture and Markovian dependence on previous-day count; LEV adult+pediatric fit used in the publication to scaffold a brivaracetam pediatric extrapolation).
* Add MPD6 Sutent sunitinib NSCLC PK/PD/tumor-growth model (no linked publication) [DDMODEL00000231] -- semi-mechanistic 15-state ODE model with parent + metabolite 2-compartment oral PK, four indirect-response PD biomarkers, sphere-volume tumor growth, and three resistance / memory chains; MDL-only deposit, validated by F.2 self-consistency only (no `.lst`, no companion paper, no simulated dataset).
* Add Leuppi-Taegtmeyer 2019 colistin ([doi:10.1128/AAC.01957-18](https://doi.org/10.1128/AAC.01957-18)) [DDMODEL00000295] -- critically ill adults receiving colistimethate sodium / colistin during continuous renal replacement therapy.
* Add Voller 2017 phenobarbital ([doi:10.1016/j.ejps.2017.05.026](https://doi.org/10.1016/j.ejps.2017.05.026)) [DDMODEL00000256] -- preterm and term newborns receiving a phenobarbital loading dose followed by oral maintenance.
* Add BAST 2017 PTTE four-event teaching library [DDMODEL00000243] -- 200 simulated patients, no linked publication; four parametric time-to-event hazard models packaged separately as `NA_NA_tte_gompertz` (Event 1, exponential / NEUT + AGE), `NA_NA_tte_gompertz_ev2` (Event 2, Gompertz / first-week AUC), `NA_NA_tte_lognormal` (Competing Event 1, log-normal / AGE), and `NA_NA_tte_loglogistic` (Competing Event 2, log-logistic / no covariate).
* Add Cook 2016 paracetamol ([doi:10.1007/s40262-016-0408-1](https://doi.org/10.1007/s40262-016-0408-1)) [DDMODEL00000271] -- term and preterm newborns receiving IV paracetamol (parent + glucuronide + sulphate plasma compartments with cumulative urinary excretion).
* Add Wang 2013 morphine ([doi:10.1007/s40261-013-0097-6](https://doi.org/10.1007/s40261-013-0097-6)) [DDMODEL00000269] -- neonates / infants / children / adolescents / adults across the entire paediatric age range (bodyweight-dependent allometric exponent on morphine clearance).
* Add Stevens 2012 remoxipride ([doi:10.1007/s10928-012-9262-4](https://doi.org/10.1007/s10928-012-9262-4)) [DDMODEL00000268] -- preclinical Wistar rats (mechanism-based PK/PD pool model for the prolactin response with brain-ECF-driven Emax stimulation and positive feedback on synthesis).
* Add Allegaert 2015 paracetamol ([doi:10.1186/s12871-015-0144-3](https://doi.org/10.1186/s12871-015-0144-3)) [DDMODEL00000267] -- young women across pregnancy / postpartum / non-pregnant-volunteer reproductive states with and without oral contraceptive use (8-compartment IV propacetamol PK with glucuronide and sulphate metabolites and cumulative-urine outputs).
* Update miridesap (CPHPC) [DDMODEL00000262] to the FINAL model from Sahota 2015 ([doi:10.1002/psp4.15](https://doi.org/10.1002/psp4.15)) -- healthy volunteers (CPH113776) and patients with systemic amyloidosis (CPH114527); two-compartment CPHPC PK plus two-compartment SAP turnover with CPHPC-SAP binding/internalization, now with the full final-model covariate panel (AMLIVER on Q4, AMLOAD on V4, SEXF on SAP baseline, CRCL on CL) and parameter values cross-checked against Sahota 2015 Table 2. Ratifies new `AMLOAD` and `AMLIVER` covariate canonicals and adds a `SEX 1=male/2=female` source-alias note to `SEXF`.
* Add Svensson 2018 rifampicin ([doi:10.1002/cpt.778](https://doi.org/10.1002/cpt.778)) [DDMODEL00000244] -- adult pulmonary tuberculosis patients on high-dose rifampicin (PanACEA HIGHRIF1 dose-escalation trial; auto-induction enzyme turnover with Michaelis-Menten clearance and dose-dependent bioavailability).
* Add Wilbaux 2015 prostate ([doi:10.1002/psp4.34](https://doi.org/10.1002/psp4.34)) [DDMODEL00000261] -- adults with metastatic castration-resistant prostate cancer on chemotherapy and/or hormonotherapy (joint K-PD model of CTC count and PSA kinetics with cell-lifespan delay and negative-binomial count likelihood).
* Add Clewe 2018 rifampicin ([doi:10.1093/jac/dkx380](https://doi.org/10.1093/jac/dkx380)) [DDMODEL00000259] -- in vitro M. tuberculosis B1585 time-kill under rifampicin + isoniazid + ethambutol triple-combination GPDI model (scenario 4).
* Add van Rongen 2018 midazolam ([doi:10.1007/s40262-017-0579-4](https://doi.org/10.1007/s40262-017-0579-4)) [DDMODEL00000250] -- 19 obese adolescents and 20 morbidly obese adults (CYP3A probe, oral + IV midazolam, sub-population-specific weight covariates on CL and Vp).
* Add Vet 2016 midazolam ([doi:10.1164/rccm.201510-2114OC](https://doi.org/10.1164/rccm.201510-2114OC)) [DDMODEL00000249] -- critically ill children in the paediatric ICU receiving continuous IV midazolam (per-stratum CL by number of failing organs, with CRP and body weight covariates).
* Add Knibbe 2009 morphine ([doi:10.2165/00003088-200948060-00003](https://doi.org/10.2165/00003088-200948060-00003)) [DDMODEL00000248] -- preterm neonates, term newborns, infants and toddlers <3 years (joint parent-metabolite popPK for morphine + M3G + M6G with body-weight allometric scaling and a postnatal-age-stratified glucuronidation step at PNA = 10 days).
* Add Zurlinden 2016 paracetamol ([doi:10.1007/s13318-015-0253-x](https://doi.org/10.1007/s13318-015-0253-x)) [DDMODEL00000237] -- healthy adults receiving a single 1000 mg oral dose (whole-body PBPK with Michaelis-Menten liver metabolism, cofactor depletion / resynthesis, and renal elimination of APAP and its conjugates AG and AS; first MCSim-source extraction in `inst/modeldb/ddmore/`).
* Add Zierhut 2008 osteoprotegerin ([doi:10.1007/s10928-008-9093-5](https://doi.org/10.1007/s10928-008-9093-5)) [DDMODEL00000233] -- healthy postmenopausal women receiving single IV or SC Fc-osteoprotegerin (two-peripheral-compartment PK with parallel linear and Michaelis-Menten elimination, logistic-style SC bioavailability F = FSC / (1 + FSC), and an indirect-response uNTX biomarker turnover model with route-conditional PK observation residual SD).
* Add Jonsson 2005 disufenton ([doi:10.2165/00003088-200544080-00007](https://doi.org/10.2165/00003088-200544080-00007)) [DDMODEL00000245] -- adults with acute ischaemic or haemorrhagic stroke receiving 72-h IV NXY-059 infusion across CRCL 20-143 mL/min.
* Add Wilkins 2008 rifampicin ([doi:10.1128/AAC.00461-07](https://doi.org/10.1128/AAC.00461-07)) [DDMODEL00000280] -- adult South African pulmonary tuberculosis patients.
* Add NA NA lidocaine ([DDMORE Foundation Model Repository: DDMODEL00000281](https://repository.ddmore.eu/model/DDMODEL00000281)) [DDMODEL00000281] -- population unspecified (no linked publication; lidocaine + MEGX + GX + 2,6-xylidide ADVAN5 parent-and-three-metabolites model).
* Add Laouenan 2015 ribavirin ([doi:10.1002/psp4.8](https://doi.org/10.1002/psp4.8)) [DDMODEL00000285] -- HCV genotype 1 cirrhotic patients on telaprevir- or boceprevir-based triple therapy (hemoglobin turnover PD model with upstream-PK ribavirin regressors).
* Add Clewe 2016 rifampicin ([doi:10.1093/jac/dkv416](https://doi.org/10.1093/jac/dkv416)) [DDMODEL00000240] -- in vitro Mycobacterium tuberculosis H37Rv multistate natural-growth scaffold (fast / slow / non-multiplying bacterial states; rifampicin exposure-response not encoded in the bundled scenario).
* Add Germovsek 2016 gentamicin ([doi:10.1128/AAC.00577-16](https://doi.org/10.1128/AAC.00577-16)) [DDMODEL00000238] -- neonates and infants receiving gentamicin.
* Add Novakovic 2017 cladribine ([doi:10.1208/s12248-016-9977-z](https://doi.org/10.1208/s12248-016-9977-z)) [DDMODEL00000223] -- adults with relapsing-remitting multiple sclerosis (8-item EDSS Item Response Theory disease-progression model with FREM covariates).
* Add Schindler 2016 sunitinib ([doi:10.1002/psp4.12057](https://doi.org/10.1002/psp4.12057)) [DDMODEL00000221] -- adults with imatinib-resistant or imatinib-intolerant advanced gastrointestinal stromal tumor.
* Add Jager 2011 gemtuzumab ([doi:10.1371/journal.pone.0024265](https://doi.org/10.1371/journal.pone.0024265)) [DDMODEL00000229] -- patients with acute myeloid leukemia (mechanism-based PKPD with explicit CD33 binding and leukemic-blast cell killing; DDMORE Monolix re-fit with added peripheral PK compartment).
* Add Khan 2015 ciprofloxacin ([doi:10.1093/jac/dkv233](https://doi.org/10.1093/jac/dkv233)) [DDMODEL00000225] -- in vitro time-kill experiments on E. coli K-12 wild-type and quinolone-resistant single-step mutants.
* Add Ter Heine 2014 tamoxifen ([doi:10.1111/bcp.12388](https://doi.org/10.1111/bcp.12388)) [DDMODEL00000212] -- adults with breast cancer on chronic 20 mg PO QD tamoxifen at steady state (joint parent-metabolite popPK with CYP2D6 and CYP3A4/5 covariates on endoxifen formation).
* Add Mohamed 2016 colistin + meropenem ([doi:10.1093/jac/dkv488](https://doi.org/10.1093/jac/dkv488)) [DDMODEL00000173] -- in vitro time-kill PK/PD against P. aeruginosa wild-type ATCC 27853 and meropenem-resistant ARU552.
* Add Bizzotto 2016 glucose ([doi:10.1152/ajpendo.00045.2016](https://doi.org/10.1152/ajpendo.00045.2016)) [DDMODEL00000227] -- adults across the glucose-tolerance spectrum (mechanistic glucose-tracer kinetics simulator).
* Add Netterberg 2017 docetaxel ([doi:10.1007/s00280-017-3366-x](https://doi.org/10.1007/s00280-017-3366-x)) [DDMODEL00000224] -- adults receiving docetaxel chemotherapy (Friberg-style myelosuppression PD model with Kloft 2006 parameter set).
* Add Zecchin 2016 survival ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000218] -- adults with advanced epithelial ovarian cancer.
* Add Hansson 2013 sunitinib fatigue Markov + proportional-odds model ([doi:10.1038/psp.2013.62](https://doi.org/10.1038/psp.2013.62)) [DDMODEL00000222] -- adults with imatinib-resistant gastrointestinal stromal tumors.
* Add Svensson 2016 bedaquiline ([doi:10.1002/psp4.12147](https://doi.org/10.1002/psp4.12147)) [DDMODEL00000219] -- adults with multidrug-resistant tuberculosis (parent + N-desmethyl M2 metabolite, time-varying weight and albumin).
* Add Jonsson 2011 ethambutol ([doi:10.1128/AAC.00274-11](https://doi.org/10.1128/AAC.00274-11)) [DDMODEL00000220] -- adult South African pulmonary tuberculosis patients.
* Add Zecchin 2016 tumorovarian ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000217] -- adults with advanced epithelial ovarian cancer receiving carboplatin monotherapy or carboplatin + gemcitabine combination chemotherapy.
* Add Girard 2012 pimasertib ([www.page-meeting.org/?abstract=2458](https://www.page-meeting.org/?abstract=2458)) [DDMODEL00000215] -- adults with advanced solid tumours and hematological malignancies (joint K-PD / cumulative-logit Markov ocular-AE-grade and Weibull-TTE dropout model).
* Add Lestini 2015 TGF-beta inhibitor ([doi:10.1007/s11095-015-1693-3](https://doi.org/10.1007/s11095-015-1693-3)) [DDMODEL00000192] -- simulated 50-subject oncology cohort (one-compartment first-order absorption PK + indirect-response biomarker turnover; first nlmixr2lib `inst/modeldb/ddmore/` entry).
* Add Li 2006 meropenem ([doi:10.1177/0091270006291035](https://doi.org/10.1177/0091270006291035)) [DDMODEL00000213] -- adult patients receiving meropenem.
* Add Friberg 2002 paclitaxel ([doi:10.1200/JCO.2002.02.140](https://doi.org/10.1200/JCO.2002.02.140)) [DDMODEL00000186] -- adult cancer patients receiving paclitaxel chemotherapy (semi-mechanistic myelosuppression PK/PD with leukocyte output).
* Add Hansson 2013 sunitinib ([doi:10.1038/psp.2013.61](https://doi.org/10.1038/psp.2013.61)) [DDMODEL00000197] -- adults with imatinib-resistant gastrointestinal stromal tumours (biomarker PD model for VEGF, sVEGFR-2, sVEGFR-3, sKIT).
* Add Henin 2009 capecitabine ([doi:10.1038/clpt.2008.220](https://doi.org/10.1038/clpt.2008.220)) [DDMODEL00000214] -- adults with metastatic colorectal or advanced/metastatic breast cancer (Markov-proportional-odds model for hand-and-foot syndrome grades 0-2 with K-PD capecitabine exposure).
* Add Plan 2012 pain ([doi:10.1038/clpt.2011.301](https://doi.org/10.1038/clpt.2011.301)) [DDMODEL00000194] -- adults in placebo arm of three Phase III neuropathic-pain trials (Markov Integer Model for daily 11-point Likert pain scores).
* Add Kloft 2004 sibrotuzumab ([doi:10.1023/B:DRUG.0000006173.72210.1c](https://doi.org/10.1023/B:DRUG.0000006173.72210.1c)) [DDMODEL00000195] -- adults with metastatic FAP-positive solid tumors.
* Change: Standardize Xie 2019 agomelatine parameter, compartment, IOV-eta, and residual-error naming. Compartments rename `DEPOT1` / `DEPOT2` / `LIVER` / `CENTPRNT` / `ALMTPERI` / `METB3OH` / `METB7DM` / `METB7DMPERI` to canonical `depot` / `depot2` / `liver` / `central` / `peripheral1` / `central_3oh` / `central_7dm` / `peripheral1_7dm`; the per-occasion IOV etas (formerly `e.IOV1`-`e.IOV5` and `eta17`-`eta35`) become descriptive `etaiov_<param>_<occ>` names paired with the structural parameter they apply to (`k13`, `alag2`, `k23`, `clint`, `fpop`); residual-error parameters `sdalmt` / `sd3oh` / `sd7dm` become canonical `addSd_lcalmt` / `addSd_lc3oh` / `addSd_lc7dm`; primary/secondary `1` suffixes are dropped where redundant (`fDepot1` -> `fDepot`, `alag1` -> `alag`, `ltvalag1` -> `ltvalag`, `etaltvalag1` -> `etaltvalag`, `F1` -> `fpop`, `etaF1` -> `etafpop`); per-metabolite molecular-weight ratios are renamed descriptively (`mpr1` / `mpr2` -> `mpr_3oh` / `mpr_7dm`). Convention infrastructure (`R/conventions.R`, `R/checkModelConventions.R`) extended to register the `liver` compartment, the `depot[0-9]+` numbered-depot pattern, the `3oh` / `7dm` metabolite suffixes, and the `etaiov_<param>_<occ>` IOV-eta naming pattern.
* Change: Standardize parameter, compartment, covariate-effect, and residual-error naming across the model library. Volumes use `vc` / `vp` / `vp2`; Michaelis-Menten Vmax uses `lvmax` / `vmax`; multi-component CL uses `lcl_ss` / `lcl_time`; shared exponents use `e_<cov>_<param1>_<param2>`; reversed-order covariate effects are flipped to `e_<cov>_<param>`. Parent-metabolite ADC models drop the parent `_adc` suffix and add canonical `<canonical>_<metab>` compartments (`central_mmae`, `central_dxd`, `central_sn38`, `central_tab`, `central_nab`, `central_dm4`, `central_medm4`); metabolite outputs become `Cc_<metab>`. Residual-error parameters now use the parameter-first form `propSd_<X>` / `addSd_<X>` (e.g. `propSd_dxd`, `addSd_tumorSize`) for every non-parent output, while the parent observation `Cc` keeps bare `propSd` / `addSd`. Convention infrastructure (`R/checkModelConventions.R`) extended to enforce the new canonical forms and flag deprecated names.
* Fix Robbie 2012 palivizumab population metadata after revisit against the full-text PDF: pediatric `n` corrected from 1,767 (derived) to 1,684 (paper PK dataset), `sex_female_pct` denominator from 1,660 to 1,684, and added a vignette note clarifying that the 2012 erratum's equation 10 correction (6,900 -> 16,900) is a NONMEM `$PRIOR` variance for adult Vp, not a final-model parameter. All `ini()` parameter values continue to match Table 2 exactly.
* Add de Vries Schultink 2020 zenocutuzumab (MCLA-128) ([doi:10.1007/s40262-020-00858-2](https://doi.org/10.1007/s40262-020-00858-2)) -- adults with advanced solid tumors (first bispecific antibody in nlmixr2lib).
* Add Almquist 2022 anifrolumab ([doi:10.1002/jcph.2055](https://doi.org/10.1002/jcph.2055)) -- adults with moderate-to-severe systemic lupus erythematosus and healthy volunteers (QSS-TMDD with dynamic IFNAR1 receptor pool and time-varying linear clearance).
* Fix Clegg 2024 nirsevimab vignette Figure 4: pass disjoint `id_offset` per cohort and carry `trial` via `rxSolve(keep = )` so the four panels no longer collapse into a single Frankenstein-subject simulation (predictions had been ~3x too high).
* Add Lu 2022 patritumab deruxtecan ([doi:10.1002/jcph.2137](https://doi.org/10.1002/jcph.2137)) -- adults with HER3-expressing advanced or metastatic solid tumors (NSCLC, breast cancer, colorectal cancer).
* Add Cheng 2026 immunoglobulin ([doi:10.1002/bcp.70420](https://doi.org/10.1002/bcp.70420)) -- children with primary immunodeficiency or secondary antibody deficiency receiving intravenous immunoglobulin replacement therapy.
* Add Frey 2010 tocilizumab ([doi:10.1177/0091270009350623](https://doi.org/10.1177/0091270009350623)) -- adults with moderate-to-severe rheumatoid arthritis.
* Add Hayashi 2007 omalizumab ([doi:10.1111/j.1365-2125.2006.02803.x](https://doi.org/10.1111/j.1365-2125.2006.02803.x)) -- Japanese adults with atopic asthma or seasonal allergic rhinitis.
* Add Hu 2014 bapineuzumab ([doi:10.1002/jcph.393](https://doi.org/10.1002/jcph.393)) -- adults with mild-to-moderate Alzheimer's disease.
* Add Huang 2017 VRC01 ([doi:10.1080/19420862.2017.1311435](https://doi.org/10.1080/19420862.2017.1311435)) -- HIV-uninfected healthy adults receiving IV or SC HIV-1 broadly neutralizing monoclonal antibody.
* Add Ide 2020 elotuzumab ([doi:10.1002/jcph.1698](https://doi.org/10.1002/jcph.1698)) -- Japanese and non-Japanese adults with multiple myeloma.
* Add Marquez-Megias 2023 adalimumab ([doi:10.3390/biomedicines11102822](https://doi.org/10.3390/biomedicines11102822)) -- adults with inflammatory bowel disease.
* Add Nestorov 2014 factor VIII Fc fusion protein ([doi:10.1002/cpdd.167](https://doi.org/10.1002/cpdd.167)) -- previously treated patients with severe hemophilia A.
* Add Perez-Ruixo 2025 posdinemab ([doi:10.1002/cpt.70173](https://doi.org/10.1002/cpt.70173)) -- healthy adults and adults with Alzheimer's disease.
* Add Petrov 2024 romiplostim ([doi:10.1002/cpdd.1494](https://doi.org/10.1002/cpdd.1494)) -- adults with chronic immune thrombocytopenia (ITP).
* Add Pouzin 2022 tusamitamab ravtansine ([doi:10.1007/s10928-021-09799-0](https://doi.org/10.1007/s10928-021-09799-0)) -- adults with advanced solid tumors expressing CEACAM5 (multi-analyte ADC + NAB + DM4 + MeDM4 model).
* Add Takeuchi 2023 ozoralizumab ([doi:10.1002/jcph.2380](https://doi.org/10.1002/jcph.2380)) -- Japanese adults with rheumatoid arthritis.
* Add Toukam 2025 BIIB107 ([doi:10.1002/jcph.70109](https://doi.org/10.1002/jcph.70109)) -- healthy adult volunteers (first-in-human study supporting MS dose optimization).
* Add Yao 2018 guselkumab ([doi:10.1002/jcph.1063](https://doi.org/10.1002/jcph.1063)) -- adults with moderate-to-severe plaque psoriasis.
* Add Frey 2013 tocilizumab DAS28 exposure-response ([doi:10.1177/0091270012437585](https://doi.org/10.1177/0091270012437585)) -- adults with moderate-to-severe rheumatoid arthritis (OPTION + TOWARD phase III pool).
* Add Hwang 2023 monalizumab ([doi:10.1002/jcph.2220](https://doi.org/10.1002/jcph.2220)) -- adults with advanced solid tumors or squamous cell carcinoma of the head and neck.
* Add Suri 2018 brentuximab vedotin ([doi:10.1002/cpt.1037](https://doi.org/10.1002/cpt.1037)) -- adults with CD30-positive Hodgkin lymphoma, systemic anaplastic large-cell lymphoma, mycosis fungoides, or primary cutaneous ALCL pooled across six clinical studies (including the phase III ALCANZA trial in CTCL).
* Add Zhong 2026 abatacept ([doi:10.1002/jcph.70156](https://doi.org/10.1002/jcph.70156)) -- pooled adults with rheumatoid arthritis, patients aged 2-17 years with polyarticular juvenile idiopathic arthritis, and patients aged 6+ years with hematologic malignancies receiving HLA-matched unrelated-donor HSCT.
* Add Mukai 2019 mogamulizumab ([doi:10.1002/jcph.1564](https://doi.org/10.1002/jcph.1564)) -- adults with cutaneous T-cell lymphoma or adult T-cell lymphoma.
* Add Lin 2024 casirivimab ([doi:10.1007/s11095-024-03764-5](https://doi.org/10.1007/s11095-024-03764-5)) -- pediatric and adult subjects, non-infected or with SARS-CoV-2 infection or household contacts.
* Add Faelens 2021 infliximab ([doi:10.3390/pharmaceutics13101623](https://doi.org/10.3390/pharmaceutics13101623)) -- adults with moderate-to-severe ulcerative colitis.
* Add Sathe 2024 sacituzumab govitecan ([doi:10.1007/s40262-024-01366-3](https://doi.org/10.1007/s40262-024-01366-3)) -- adults with metastatic triple-negative breast cancer and other solid tumors.
* Add Takahashi 2023 abatacept ([doi:10.1182/blood.2023020035](https://doi.org/10.1182/blood.2023020035)) -- pooled adult/pediatric rheumatoid-arthritis / pJIA patients and adult/pediatric hematopoietic-cell-transplant recipients in the ABA2 trial.
* Add Hong 2025 datopotamab deruxtecan ([doi:10.1002/psp4.70118](https://doi.org/10.1002/psp4.70118)) -- adults with advanced or metastatic NSCLC and breast cancer.
* Add Fau 2020 isatuximab ([doi:10.1002/psp4.12561](https://doi.org/10.1002/psp4.12561)) -- adults with relapsed/refractory multiple myeloma.
* Add Lu 2019 polatuzumab vedotin ([doi:10.1002/psp4.12482](https://doi.org/10.1002/psp4.12482)) -- adults with relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma.
* Add Okada 2025 rocatinlimab ([doi:10.1002/psp4.70069](https://doi.org/10.1002/psp4.70069)) -- adults with moderate-to-severe atopic dermatitis (pooled with plaque psoriasis, ulcerative colitis, and healthy volunteers).
* Add Zhou 2025 brentuximab vedotin ([doi:10.1002/cpt.3629](https://doi.org/10.1002/cpt.3629)) -- pediatric patients (5-18 years) with newly diagnosed Hodgkin lymphoma or relapsed/refractory systemic anaplastic large-cell lymphoma.
* Add Lin 2024 pozelimab ([doi:10.1007/s10928-024-09941-8](https://doi.org/10.1007/s10928-024-09941-8)) -- healthy adult volunteers, adults with paroxysmal nocturnal hemoglobinuria, and pediatric and adult patients with CHAPLE disease.
* Add Yin 2021 trastuzumab deruxtecan ([doi:10.1002/cpt.2096](https://doi.org/10.1002/cpt.2096)) -- adults with HER2-positive breast cancer or other HER2-expressing solid tumors.
* Add Hwang 2022 tremelimumab ([doi:10.1111/bcp.15622](https://doi.org/10.1111/bcp.15622)) -- adults with advanced solid tumours.
* Add Papathanasiou 2025 belantamab mafodotin ([doi:10.1007/s40262-025-01508-1](https://doi.org/10.1007/s40262-025-01508-1)) -- adults with relapsed/refractory multiple myeloma (ADC moiety only).
* Add Kuchimanchi 2024 dostarlimab ([doi:10.1111/bcp.16325](https://doi.org/10.1111/bcp.16325)) -- adults with primary advanced or recurrent endometrial cancer (RUBY Part 1) and advanced solid tumours (GARNET).
* Add Diao 2016 daclizumab CD25 occupancy ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Rewrite Diao 2016 daclizumab CD25 occupancy ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) as a kinetic-binding kon/koff ODE (replacing the previous sigmoidal-Emax-with-desaturation-only form that under-predicted Figure 1A saturation onset); kon and koff calibrated to reproduce both Figure 1A (~7 h saturation) and Figure 1B (~24-week baseline return).
* Add Diao 2016 daclizumab CD56 bright NK expansion ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Diao 2016 daclizumab Treg reduction ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Fiedler-Kelly 2020 fremanezumab exposure-response ([doi:10.1111/head.13855](https://doi.org/10.1111/head.13855)) -- adults with episodic migraine and adults with chronic migraine (two PD-only models: `FiedlerKelly_2020_fremanezumab_em` and `FiedlerKelly_2020_fremanezumab_cm`).
* Add Koopman 2023 factor IX-Fc ([doi:10.1111/bcp.15881](https://doi.org/10.1111/bcp.15881)) -- children, adolescents, and adults with haemophilia B (real-world data including patients aged < 12 years).
* Add Le Tilly 2021 trastuzumab ([doi:10.1002/cpt.2188](https://doi.org/10.1002/cpt.2188)) -- adults with HER2+ breast cancer leptomeningeal carcinomatosis receiving intrathecal and intravenous trastuzumab.
* Add Wang 2024 sugemalimab ([doi:10.1111/bcp.16276](https://doi.org/10.1111/bcp.16276)) -- adults with advanced solid tumours or lymphomas across nine Phase I-III sugemalimab trials.
* Add Yang 2024 axatilimab ([doi:10.1002/cpt.3503](https://doi.org/10.1002/cpt.3503)) -- pooled healthy adults, adults with advanced solid tumors, and adults / children with chronic graft-versus-host disease.
* Add Hood 2021 MEDI7836 ([doi:10.3390/pharmaceutics13050613](https://doi.org/10.3390/pharmaceutics13050613)) -- healthy adult males in a first-in-human single-ascending-dose trial.
* Amend Castro-Suarez 2020 nimotuzumab: V1 decreased 53% for 50 mg cohort per Figure 4 visual inspection (direction not stated in paper text; corresponding author contacted).
* Add Castro-Suarez 2020 nimotuzumab ([doi:10.3390/pharmaceutics12121147](https://doi.org/10.3390/pharmaceutics12121147)) -- adults with autosomal dominant polycystic kidney disease.
* Add Yang 2021 cemiplimab ([doi:10.1007/s10928-021-09739-y](https://doi.org/10.1007/s10928-021-09739-y)) -- adults with advanced solid tumors including cutaneous squamous cell carcinoma.
* Add Papachristos 2020 bevacizumab ([doi:10.3390/ijms21113753](https://doi.org/10.3390/ijms21113753)) -- adults with metastatic colorectal cancer (three co-equal final models: descriptive PK, binding QSS TMDD, and PK/PD Imax inhibition of free VEGF-A).
* Add Ngo 2020 HL2351 ([doi:10.1002/psp4.12552](https://doi.org/10.1002/psp4.12552)) -- healthy adult Korean men.
* Add Wojciechowski 2022 domagrozumab ([doi:10.1111/cts.13418](https://doi.org/10.1111/cts.13418)) -- healthy adult volunteers and pediatric patients with Duchenne muscular dystrophy.
* Add Yu 2022 ofatumumab ([doi:10.1007/s40263-021-00895-w](https://doi.org/10.1007/s40263-021-00895-w)) -- adults with relapsing multiple sclerosis.
* Add Melhem 2022 dostarlimab ([doi:10.1111/bcp.15339](https://doi.org/10.1111/bcp.15339)) -- adults with advanced solid tumours.
* Add Brillac 2025 isatuximab ([doi:10.1007/s00280-025-04832-2](https://doi.org/10.1007/s00280-025-04832-2)) -- pediatric and adult patients with relapsed/refractory acute leukemias.
* Add Wu 2024 inotuzumab ozogamicin ([doi:10.1007/s40262-024-01386-z](https://doi.org/10.1007/s40262-024-01386-z)) -- pediatric and adult patients with relapsed/refractory B-cell precursor acute lymphoblastic leukemia and adults with B-cell non-Hodgkin's lymphoma.
* Fix vignettes: derive concentration units in `labs()` from model `$units` metadata; replace inline trapezoidal NCA with PKNCA; add PKNCA sections to nalmefene and clesrovimab vignettes; add PKNCA treatment grouping to benralizumab and durvalumab vignettes.
* Add Chen 2022 guselkumab ([doi:10.1111/cts.13197](https://doi.org/10.1111/cts.13197)) -- adults with active psoriatic arthritis (DISCOVER-1 and DISCOVER-2 phase 3 trials).
* Add `concentration` field to `units` metadata for `Wang_2017_benralizumab` and `Ogasawara_2020_durvalumab` models.
* Add Chen 2020 luspatercept ([doi:10.1002/psp4.12515](https://doi.org/10.1002/psp4.12515)) -- adults with anemia due to lower-risk myelodysplastic syndromes.
* Add Martinez 2019 alirocumab ([doi:10.1007/s40262-018-0669-y](https://doi.org/10.1007/s40262-018-0669-y)) -- healthy volunteers and adults with hypercholesterolemia.
* Add Li 2017 brentuximab vedotin ([doi:10.1002/jcph.920](https://doi.org/10.1002/jcph.920)) -- adults with relapsed/refractory CD30-expressing hematologic malignancies (Hodgkin lymphoma and systemic anaplastic large cell lymphoma).
* Add Quartino 2019 trastuzumab ([doi:10.1007/s00280-018-3728-z](https://doi.org/10.1007/s00280-018-3728-z)) -- adults with metastatic breast cancer, early breast cancer, advanced gastric cancer, or other solid tumors.
* Add Suleiman 2019 risankizumab ([doi:10.1007/s40262-019-00759-z](https://doi.org/10.1007/s40262-019-00759-z)) -- healthy volunteers and adults with moderate-to-severe plaque psoriasis.
* Add Berends 2019 infliximab ([doi:10.1007/s10928-019-09652-5](https://doi.org/10.1007/s10928-019-09652-5)) -- adults with moderate-to-severe ulcerative colitis.
* Add Nikanjam 2019 siltuximab ([doi:10.1007/s00280-019-03939-7](https://doi.org/10.1007/s00280-019-03939-7)) -- adults pooled across healthy volunteers and oncology cohorts including Castleman's disease and smoldering multiple myeloma.
* Add Sanghavi 2020 ipilimumab ([doi:10.1002/psp4.12477](https://doi.org/10.1002/psp4.12477)) -- adults with advanced solid tumors receiving ipilimumab alone or in combination with nivolumab.
* Add Zhang 2019 nivolumab ([doi:10.1002/psp4.12476](https://doi.org/10.1002/psp4.12476)) -- adults with advanced solid tumors, alone or in combination with ipilimumab or chemotherapy.
* Add Wang 2020 ontamalimab ([doi:10.1002/jcph.1590](https://doi.org/10.1002/jcph.1590)) -- adults with moderate-to-severe ulcerative colitis or Crohn's disease.
* Add Hanzel 2021 infliximab CT-P13 ([doi:10.1111/apt.16609](https://doi.org/10.1111/apt.16609)) -- adults with Crohn's disease or ulcerative colitis.
* Add Zhou 2021 belimumab ([doi:10.1007/s40268-021-00363-2](https://doi.org/10.1007/s40268-021-00363-2)) -- adult and pediatric patients with systemic lupus erythematosus (Chinese and non-Chinese).
* Add Aguiar 2021 ustekinumab ([doi:10.3390/pharmaceutics13101587](https://doi.org/10.3390/pharmaceutics13101587)) -- adults with active Crohn's disease.
* Add Li 2019 abatacept ([doi:10.1002/jcph.1308](https://doi.org/10.1002/jcph.1308)) -- adults with rheumatoid arthritis.
* Add Gandhi 2021 abatacept ([doi:10.1002/jcph.1781](https://doi.org/10.1002/jcph.1781)) -- pooled adults with rheumatoid arthritis and patients aged 2-17 years with polyarticular juvenile idiopathic arthritis.
* Add Mulyukov 2018 ranibizumab ([doi:10.1002/psp4.12322](https://doi.org/10.1002/psp4.12322)) -- anti-VEGF-naive adults with neovascular age-related macular degeneration.
* Add Bajaj 2017 nivolumab ([doi:10.1002/psp4.12143](https://doi.org/10.1002/psp4.12143)) -- patients with advanced solid tumors (melanoma, NSCLC, RCC, other).
* Add Kielbasa 2020 galcanezumab ([doi:10.1002/jcph.1511](https://doi.org/10.1002/jcph.1511)) -- healthy adults and adults with episodic or chronic migraine.
* Add Othman 2014 daclizumab ([doi:10.1007/s40262-014-0159-9](https://doi.org/10.1007/s40262-014-0159-9)) -- healthy adult volunteers.
* Add Timmermann 2019 brodalumab ([doi:10.1111/bcpt.13202](https://doi.org/10.1111/bcpt.13202)) -- adults with moderate-to-severe plaque psoriasis.
* Add Kuchimanchi 2018 evolocumab ([doi:10.1007/s10928-018-9592-y](https://doi.org/10.1007/s10928-018-9592-y)) -- healthy adults and patients with hypercholesterolemia (including heterozygous familial hypercholesterolemia and statin-intolerant cohorts).
* Add Mo 2018 olaratumab ([doi:10.1007/s40262-017-0562-0](https://doi.org/10.1007/s40262-017-0562-0)) -- patients with advanced or metastatic cancer (soft tissue sarcoma, NSCLC, GIST, glioblastoma multiforme).
* Add Kovalenko 2016 dupilumab ([doi:10.1002/psp4.12136](https://doi.org/10.1002/psp4.12136)) -- healthy volunteers and adults with moderate-to-severe atopic dermatitis.
* Add Narwal 2013 sifalimumab ([doi:10.1007/s40262-013-0085-2](https://doi.org/10.1007/s40262-013-0085-2)) -- adults with moderately active systemic lupus erythematosus.
* Add Valenzuela 2025 nipocalimab ([doi:10.1002/psp4.70109](https://doi.org/10.1002/psp4.70109)) -- healthy adults and adults with generalized myasthenia gravis.
* Add Hua 2015 anrukinzumab ([doi:10.1111/bcp.12589](https://doi.org/10.1111/bcp.12589)) -- healthy volunteers, asthma, and ulcerative colitis patients.
* Add Robbie 2012 palivizumab ([doi:10.1128/AAC.06446-11](https://doi.org/10.1128/AAC.06446-11)) -- adults and pediatric patients at high risk of RSV lower respiratory tract disease.
* Add Bender 2014 trastuzumab emtansine (T-DM1, ADC) ([doi:10.1208/s12248-014-9618-3](https://doi.org/10.1208/s12248-014-9618-3)) -- preclinical: rats + cynomolgus monkeys; two model variants (reduced + mechanistic DAR chain).
* Add Rosario 2015 vedolizumab ([doi:10.1111/apt.13243](https://doi.org/10.1111/apt.13243)) -- adults with moderately-to-severely active ulcerative colitis or Crohn's disease.
* Add Long 2017 necitumumab ([doi:10.1007/s40262-016-0452-x](https://doi.org/10.1007/s40262-016-0452-x)) -- cancer patients with advanced solid tumors.
* Add Wade 2015 certolizumab ([doi:10.1002/jcph.491](https://doi.org/10.1002/jcph.491)) -- adults with Crohn's disease.
* Add Lon 2013 abatacept ([doi:10.1007/s10928-013-9341-1](https://doi.org/10.1007/s10928-013-9341-1)) -- male Lewis rats with collagen-induced arthritis (preclinical).
* Add Lu 2014 trastuzumab emtansine ([doi:10.1007/s00280-014-2500-2](https://doi.org/10.1007/s00280-014-2500-2)) -- adults with HER2-positive locally advanced or metastatic breast cancer.
* Add Gupta 2016 amatuximab ([doi:10.1007/s00280-016-2984-z](https://doi.org/10.1007/s00280-016-2984-z)) -- adults with advanced mesothelin-expressing cancers including malignant pleural mesothelioma.
* Add Zheng 2016 sifalimumab ([doi:10.1111/bcp.12864](https://doi.org/10.1111/bcp.12864)) -- adults with systemic lupus erythematosus.
* Add Mould 2007 alemtuzumab ([doi:10.1111/j.1365-2125.2007.02914.x](https://doi.org/10.1111/j.1365-2125.2007.02914.x)) -- adults with B-cell chronic lymphocytic leukaemia.
* Add Farrell 2012 farletuzumab ([doi:10.1007/s00280-012-1959-y](https://doi.org/10.1007/s00280-012-1959-y)) -- women with advanced epithelial ovarian cancer.
* Add Xu 2011 sirukumab ([doi:10.1111/j.1365-2125.2011.03964.x](https://doi.org/10.1111/j.1365-2125.2011.03964.x)) -- healthy adult volunteers in a first-in-human study.
* Add Yamada 2025 zolbetuximab ([doi:10.1111/cts.70280](https://doi.org/10.1111/cts.70280)) -- patients with locally advanced unresectable or metastatic gastric/gastroesophageal junction adenocarcinoma.
* Add Jackson 2022 ixekizumab ([doi:10.1111/bcp.15034](https://doi.org/10.1111/bcp.15034)) -- paediatric patients with moderate-to-severe plaque psoriasis.
* Add Moein 2022 etrolizumab ([doi:10.1002/psp4.12846](https://doi.org/10.1002/psp4.12846)) -- patients with moderately-to-severely active ulcerative colitis.
* Add Chua 2025 mirikizumab ([doi:10.1111/cts.70320](https://doi.org/10.1111/cts.70320)) -- patients with moderately-to-severely active Crohn's disease.
* Add Ma 2020 sarilumab ANC ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) -- adults with rheumatoid arthritis.
* Add Ma 2020 sarilumab DAS28-CRP ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) -- adults with rheumatoid arthritis.
* Add Tiraboschi 2025 amlitelimab ([doi:10.1002/psp4.70121](https://doi.org/10.1002/psp4.70121)) -- adults with moderate-to-severe atopic dermatitis.
* Add Budha 2023 tislelizumab ([doi:10.1002/psp4.12880](https://doi.org/10.1002/psp4.12880)) -- patients with advanced tumors.
* Add Masters 2022 avelumab ([doi:10.1002/psp4.12771](https://doi.org/10.1002/psp4.12771)) -- patients with advanced solid tumors.
* Add Xu 2019 sarilumab ([doi:10.1007/s40262-019-00765-1](https://doi.org/10.1007/s40262-019-00765-1)) -- adults with rheumatoid arthritis.
* `checkModelConventions()` -- new function that reports deviations from the package's parameter-naming, covariate, compartment, and metadata conventions for a single model or the entire `modeldb`. Called automatically during `buildModelDb()` so convention drift surfaces at package-build time (existing grandfathered deviations continue to build). Canonical standards (PK parameter prefixes, `eta`-prefixed IIV, `propSd`/`addSd` residual error, canonical covariate column register with aliases, canonical compartment vocabulary) are codified in an internal `.nlmixr2libConventions` list that mirrors the `extract-literature-model` skill references. Addresses issue #39.
* Added canonical TMDD archetype models under `inst/modeldb/pharmacokinetics/`: `PK_1cmt_tmdd_full` (Mager & Jusko 2001), `PK_1cmt_tmdd_qss`, `PK_1cmt_tmdd_mm`, `PK_2cmt_tmdd_qss`, and `PK_2cmt_tmdd_mm` (Gibiansky et al. 2008 QSS and MM approximations). Built to the `extract-literature-model` skill conventions with `reference` / `units` / `population` metadata, per-parameter source-trace comments, and CL/V parameterization with `kel` derived inside `model()` so later transforms can re-parameterize. Replaces the 38 draft models from PR #60 (#60).
* Registered `target`, `complex`, and `total_target` as canonical TMDD compartment names in the `extract-literature-model` skill naming conventions (per @iamstein's proposal on PR #60).
* Added `vignettes/tmdd_archetypes.Rmd` comparing the five TMDD archetypes with typical-value trajectories of drug, free target, and complex, and a regime-convergence check showing QSS/MM collapsing onto the full Mager & Jusko 2001 model in the fast-binding regime of Gibiansky 2008.
* Add Thakre 2022 risankizumab ([doi:10.1007/s40744-022-00495-0](https://doi.org/10.1007/s40744-022-00495-0)) -- patients with active psoriatic arthritis.
* `addEta()`, `addResErr()`, `addDepot()`, `removeDepot()`, `addTransit()`, and `removeTransit()` now accept `model` as a deprecated alias for `ui` (issue #84). Passing `model = ...` emits a deprecation warning; passing both `ui` and `model` is an error.
* `addDepot()` and `addTransit()` now work correctly when `d/dt(central)` or `d/dt(depot)` appears at the beginning or end of the model block, or when transit-compartment ODEs and residual-error (`~`) specs are interleaved with assignment lines. The newly introduced helper and ODE lines are inserted immediately adjacent to the modified ODE so that the relative order of every pre-existing model line is preserved (#77, #78).
* Markov modeling creation functions including `createMarkovModel()` were added
* Add Kotani 2022 astegolimab ([doi:10.1002/jcph.2021](https://doi.org/10.1002/jcph.2021)) -- adults with severe asthma.
* Add Fasanmade 2009 infliximab ([doi:10.1007/s00228-009-0718-4](https://doi.org/10.1007/s00228-009-0718-4)) -- adults with moderately-to-severely active ulcerative colitis.
* Add Fiedler-Kelly 2019 fremanezumab ([doi:10.1111/bcp.14096](https://doi.org/10.1111/bcp.14096)) -- adults with chronic or episodic migraine.
* Add Hu 2026 clesrovimab ([doi:10.1002/cpt.70199](https://doi.org/10.1002/cpt.70199)) -- preterm and full-term infants.
* Add Clegg 2024 nirsevimab ([doi:10.1002/jcph.2401](https://doi.org/10.1002/jcph.2401)) -- preterm and term infants.
* Add Chawla 2023 gefapixant ([doi:10.1002/psp4.12978](https://doi.org/10.1002/psp4.12978)) -- healthy adults and adults with refractory or unexplained chronic cough.
* Verified all published-literature specific-drug and mAb-consensus models against their source papers and fixed several parameter-encoding bugs that had been latent in the package since their original addition:
  - **CarlssonPetri 2021 liraglutide**: fixed categorical covariate encoding that was zeroing individual clearance for subjects not in the indexed group. `(1 - SEXF)^e_sex_cl` -> `e_sex_cl^(1 - SEXF)` (previously evaluated `0^1.12 = 0` for females); `CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl` -> `e_age_child_cl^CHILD * e_age_adolescent_cl^ADOLESCENT` (previously evaluated `0^1.11 * 0^1.06 = 0` for adults). IIV rewritten as `omega^2 = log(1 + CV^2)` per Table 3's explicit `%CV = sqrt(exp(omega^2) - 1) * 100` footnote.
  - **Zhu 2017 lebrikizumab**: fixed IIV variance-covariance block that was storing `sqrt(variance)` (SDs) instead of variances/covariances from Table 3.
  - **Soehoel 2022 tralokinumab**: fixed IIV block (SDs -> variance-covariance matrix from Table 2 footnote `IIV = sqrt(exp(omega^2) - 1)`, correlation 0.61 applied on variances not SDs); corrected body-weight exponent on V2/V3 from `0.793` to Table 2's `0.783`.
  - **Kovalenko 2020 dupilumab**: squared the five IIV values so they store variances. Paper Methods explicitly defines `omega` as "the standard deviation [SD] of between-subject variability", and nlmixr2 `etaX ~ value` stores the variance (omega^2), so SDs needed squaring.
  - **Davda 2014 mAb consensus (PK_2cmt_mAb_Davda_2014)**: all parameters verified; no changes required.
* Filled in previously-TODO `population` metadata blocks for the five verified models above with demographics from each paper's Table 1.
* Added validation vignettes for the five verified models above (`CarlssonPetri_2021_liraglutide.Rmd`, `Zhu_2017_lebrikizumab.Rmd`, `Soehoel_2022_tralokinumab.Rmd`, `Kovalenko_2020_dupilumab.Rmd`, `PK_2cmt_mAb_Davda_2014.Rmd`). Each vignette follows the `extract-literature-model` skill conventions with a population description, source trace, virtual cohort, simulation, figure replication, and PKNCA-based NCA validation.
* Retrofit Cirincione 2017 exenatide model to the `extract-literature-model` skill conventions and fix parameter encoding bugs: `ka_max` corrected from `0.0813` to paper value `12.8` /hr, `Km` rescaled to ng/mL so units are consistent with `Cc = central / vc`, and IIV variances rewritten as `log(1 + CV^2)` rather than the `log(1 + CV)` shortcut. Replaces the character-valued `DVID` covariate with `STUDY1` / `STUDY5` binary indicators and adds a companion validation vignette with PKNCA checks against the paper's Figure 5 typical values.
* Fix endogenous-model bugs surfaced during a parameter audit:
  - **Charbonneau 2021 phenylalanine**: fixed `vd_phe` typo on the `f_gut_plasma` line (the undefined symbol made `rxSolve()` fail) -- should read `vd`; removed a stray `vd *` factor from the `daily_phe_intake` augmentation output so the reported value is in mg/day rather than (L/kg)*mg/day. All Table 4 parameter values verified against the authors' Zenodo companion notebook.
  - **Kim 2006 IgG**: corrected `V1` units label from `(mg/kg)` to `(mL/kg)` (paper Table 1); closed a missing `)` in the `ljmax` label; removed a redundant `igg_0 <- css` line. Parameter values verified against Table 1.
* Extended the `extract-literature-model` skill to cover endogenous, mechanistic, and turnover models: added naming conventions for Vmax/Km/baseline/rate parameters, a new `references/endogenous-validation.md` reference covering steady-state / perturbation-recovery / mass-balance / dimensional-analysis patterns, and an explicit dimensional-analysis item in `references/verification-checklist.md` (driven by the Kim units-label and Charbonneau `daily_phe_intake` bugs, which were only catchable by dimensional analysis).
* Retrofit the remaining specific-drug models and the Davda 2014 mAb 2-compartment model to the `extract-literature-model` skill conventions (no parameter-value changes). Structured `covariateData` with `description` / `units` / `type` / `reference_category` / `notes` / `source_name`, canonical `units` list (time/dosing/concentration), and `population` metadata blocks added (with TODO placeholders where not yet sourced). IIV etas renamed to `eta` + transformed-parameter name (e.g., `etacl` -> `etalcl`, `iiv_lka` -> `etalka`, `bsv_fpla_*` -> `etalfpla_*`). Race columns renamed to the `RACE_` prefix (`BLACK` -> `RACE_BLACK`, `ASIAN` -> `RACE_ASIAN`, `MULTIRACIAL` -> `RACE_MULTI`, `BLACK_OTH` -> `RACE_BLACK_OTH`, `ASIAN_AMIND_MULTI` -> `RACE_ASIAN_AMIND_MULTI`); `ADA` -> `ADA_POS`; `SEXM` -> `SEXF` (value inverted to keep the original effect magnitude). Touched vignettes updated to use the canonical column names. Models touched: `CarlssonPetri_2021_liraglutide`, `Clegg_2024_nirsevimab`, `Grimm_2023_gantenerumab`, `Grimm_2023_trontinemab`, `Hu_2026_clesrovimab`, `Kovalenko_2020_dupilumab`, `Kyhl_2016_nalmefene`, `PK_2cmt_mAb_Davda_2014`, `Soehoel_2022_tralokinumab`, `Xie_2019_agomelatine`, `Zhu_2017_lebrikizumab`.

# Version 0.3.2

* Add Kim 2006 model for IgG metabolism
* Add Xie 2019 agomelatine PK model
* Drop `qs` since it will be archived
* Update tumor growth inhibition models
* `addResErr()` now works with multiple-endpoint models
* Additional testing

# Version 0.3.1

* Bug fix for replacement of multiplicative expressions in `nlmixr2lib`
* phenylalanine_charbonneau_2021 had its net protein breakdown parameter corrected
* Kyhl_2016_nalmefene model was added

# Version 0.3.0

* Added ability to choose style type when modifying models.  Currently
  supported styles are: "camel" for `variablesLikeThis`, "snake" for
  `variables_like_this`, "dot" for `variables.like.this` and "blank"
  for `variableslikethis`.  This can be selected with
  `setCombineType()`.

* With the new combination style, you can change how `eta` variables
  are constructed with the `option(nlmixr2lib.etaCombineType="camel")`
  or whatever you wish it to the variable style to be.

* Added new model building framework for building models

  - **PK model building functions**

     - `addTransit()`/`removeTransit()` which were present before, but now modified and
       made a bit more robust, more closely matching literature method
       of transit compartments.

     - `addDepot()`/`removeDepot()` which were present before, but
       modified to be a bit more robust.

     - `addWeibullAbs()` which adds a Weibull absorption to a PK model

     - `convertMM()` converts linear elimination to Michaelis-Menten elimination

     - `transPK()` converts the `cl` style parameter transformations
       to various other PK transformations like `k`, `aob`, `alpha`,
       `k12`

  - **PD model building functions**

   - `addIndirectLin()` -- this adds an indirect effect model to a PK
     model that has a concentration `Cc` in the model.  This purposely
     uses a simple linear effect of `Cc*Ek` or `Cc*Ik` so it will be
     easy to parse and turn into other functional forms (like `Emax`
     or `Hill`).  If the PK model is not present it will use `Cc` as a
     covariate in a purely PD models.

   - `addIndirect()` -- this builds on `addIndirectLin()` and adds
     `Emax` or `Hill` models to a PK model. You can also set `imax=1`
     or `emax=1` to drop these parameters from being estimated in the
     model.  Additionally `hill=TRUE` will add a Hill coefficient to
     the sigmoid model.

   - `addEffectCmtLin()` -- this adds an effect compartment based on
     the `Cc` in the model.  The linear effect can be modified into
     other function forms.

   - `addDirectLin()` -- this adds a direct effect model based on the
     `Cc` in the model.

   - **Changing functional forms of Effect models**

     - `convertEmax()` changes linear effect models to Emax models

     - `convertEmaxHill()` changes linear effect models to Hill models

     - `convertQuad()` changes linear effect models to quadratic models

     - `convertLogLin()` changes linear effect models to log-linear models

   - **Changing functional forms of Baselines in non-indirect response models**

     - `addBaselineConst()` changes the zero baseline to a estimated
       constant

     - `addBaselineLin()` changes the zero baseline to a estimated
       constant and a linear constant with respect to `time`.

     - `addBaselineExp()` changes the zero baseline to a exponential
       decay with respect to time

     - `addBaseline1exp()` -- the baseline effect is changed from zero
       to to an exponential approaching to a constant (with respect to
       time).

   - **Changing model properties** (all use `addCmtProp()`)

      - `addBioavailability()` adds bioavailability property to a
        compartment

      - `addRate()` adds a modeled rate to a compartment

      - `addDur()` adds modeled duration to a compartment

      - `addIni()` adds an initial value to a compartment

      - `addLag()` adds a lag time to the a compartment

* Add Carlsson Petri (2021) liraglutide PK model
* Add Cirincione (2017) exenatide immediate-release PK model
* Add a variety of indirect response models
* Add a variety of tumor growth inhibition models and move all oncology models
  into a new model database directory
* Add a variety of double-absorption PK models
* `cp` and related `cpddSd` and `cppropSd` were renamed to `Cc`, `CcAddSd` and
  `CcPropSd` (fix #70).
* Multiple-endpoint models will have the `DV` column in the modeldb separated by
  commas.

# Version 0.2.0

* Work with the new `rxode2` version 2.0.12 `model()` and `ini()` assignment
  methods.
* Therapeutic-area specific models have begun being added.
* Models can now give the user some additional information load via the
  `message` meta-data.
* Models can now be in different directories.  The change is for ease of
  maintaining the library, it is not a change that affects users.
* A regression where `addEta()` did not change the parameter, related to a
  change in `rxode2`, was fixed.
* `addEta()` detects where to add etas more robustly when covariates are on the
  parameter.

## Models added

* Add Davda (2014) mAb consensus model
* Add Liu (2017) time-dependent clearance model based on nivolumab
* Add Kovalenko (2020) dupilumab PK model
* Add Soehoel (2022) tralokinumab PK model
* Add Zhu (2017) lebrikizumab PK model

# Version 0.1.0

* Initial version
