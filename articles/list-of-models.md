# List of models

``` r
knitr::kable(modeldb[, c("name", "description")])
```

| name                                | description                                                                                                                                      |
|:------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------|
| PK_1cmt                             | One compartment PK model with linear clearance                                                                                                   |
| PK_1cmt_des                         | One compartment PK model with linear clearance using differential equations                                                                      |
| PK_2cmt                             | Two compartment PK model with linear clearance                                                                                                   |
| PK_2cmt_des                         | Two compartment PK model with linear clearance using differential equations                                                                      |
| PK_2cmt_no_depot                    | Two compartment PK model with linear clearance using differential equations                                                                      |
| PK_2cmt_tdcl_des                    | Two compartment PK model with time-dependent clearance using differential equations (structured like nivolumab PK model)                         |
| PK_3cmt                             | Three compartment PK model with linear clearance                                                                                                 |
| PK_3cmt_des                         | Three compartment PK model with linear clearance using differential equations                                                                    |
| igg_kim_2006                        | Immunoglobulin G (IgG) model for nonlinear metabolism in healthy subjects                                                                        |
| phenylalanine_charbonneau_2021      | Phenylalanine model for absorption and metabolism in healthy subjects and patients with PKU                                                      |
| indirect_0cpt_transitEx             | Two compartment PK model with Michealis-Menten clearance using differential equations                                                            |
| indirect_1cpt_inhi_kin              | One compartment indirect response model with inhibition of kin.                                                                                  |
| indirect_1cpt_inhi_kin_CLV          | One compartment indirect response model with inhibition of kin.                                                                                  |
| indirect_1cpt_inhi_kin_r0rmaxcrmax  | One compartment indirect response model with inhibition of kin.                                                                                  |
| indirect_1cpt_inhi_kout             | One compartment indirect response model with inhibition of kout.                                                                                 |
| indirect_1cpt_inhi_kout_CLV         | One compartment indirect response model with inhibition of kout.                                                                                 |
| indirect_1cpt_inhi_kout_r0rmaxcrmax | One compartment indirect response model with inhibition of kout.                                                                                 |
| indirect_1cpt_stim_kin              | One compartment indirect response model with stimulation of kin.Parameterized using rate cosntants                                               |
| indirect_1cpt_stim_kin_CLV          | One compartment indirect response model with stimulation of kin.                                                                                 |
| indirect_1cpt_stim_kin_r0rmaxcrmax  | One compartment indirect response model with stimulation of kin.                                                                                 |
| indirect_1cpt_stim_kout             | One compartment indirect response model with stimulation of kout.Parameterized using rate cosntants                                              |
| indirect_1cpt_stim_kout_CLV         | One compartment indirect response model with stimulation of kout.                                                                                |
| indirect_1cpt_stim_kout_r0rmaxcrmax | One compartment indirect response model with stimulation of kout.                                                                                |
| indirect_circ_1cpt_inhi_kin_kin_t   | One compartment indirect response model with inhibition of kin and circadian kin_t.                                                              |
| indirect_circ_1cpt_inhi_kin_kout_t  | One compartment indirect response model with inhibition of kin and circadian kin_t.                                                              |
| indirect_circ_1cpt_inhi_kout_kin_t  | One compartment indirect response model with inhibition of kout and circadian kin_t.                                                             |
| indirect_circ_1cpt_inhi_kout_kout_t | One compartment indirect response model with inhibition of kout and circadian kin_t.                                                             |
| indirect_circ_1cpt_stim_kin_kin_t   | One compartment indirect response model with stimulation of kin and circadian kin_t.Parameterized using rate cosntants                           |
| indirect_circ_1cpt_stim_kin_kout_t  | One compartment indirect response model with stimulation of kin and circadian kout_t.Parameterized using rate constants                          |
| indirect_circ_1cpt_stim_kout_kin_t  | One compartment indirect response model with stimulation of kout and circadian kin_t.Parameterized using rate cosntants                          |
| indirect_circ_1cpt_stim_kout_kout_t | One compartment indirect response model with stimulation of kout and circadian kout_t.Parameterized using rate cosntants                         |
| indirect_prec_1cpt_inhi_CLV         | One compartment precursor-dependent indirect response model with inhibition of drug response. Parameterized with clearance and volume. (effect). |
| indirect_prec_1cpt_inhi_r0rmaxcrmax | One compartment precursor-dependent indirect response model with inhibition of drug response (effect).                                           |
| indirect_prec_1cpt_stim_CLV         | One compartment precursor-dependent indirect response model with inhibition of drug response (effect). Parameterized with clearance and volume   |
| indirect_prec_1cpt_stim_r0rmaxcrmax | One compartment precursor-dependent indirect response model with inhibition of drug response (effect). Parameterized with clearance and volume   |
| PK_2cmt_mAb_Davda_2014              | Two compartment PK model with linear clearance for average monoclonal antibodies (Davda 2014)                                                    |
| PK_double_sim_01                    | PK double absorption model with simultaneous zero order and first order absorptions                                                              |
| PK_double_sim_10                    | PK double absorption model with simultaneous first order and zero order absorptions                                                              |
| PK_double_sim_11                    | PK double absorption model with simultaneous first order absorptions                                                                             |
| CarlssonPetri_2021_liraglutide      | Liraglutide PK model in adolescents (Carlsson Petri 2021)                                                                                        |
| Cirincione_2017_exenatide           | Exenatide immediate-release PK model (Cirincione 2017)                                                                                           |
| Grimm_2023_gantenerumab             | Gantenerumab PK model (Grimm 2017)                                                                                                               |
| Grimm_2023_trontinemab              | Trontinemab PK model (Grimm 2017)                                                                                                                |
| Kovalenko_2020_dupilumab            | Dupilumab PK model (Kovalenko 2020)                                                                                                              |
| Kyhl_2016_nalmefene                 | NA                                                                                                                                               |
| Soehoel_2022_tralokinumab           | Tralokinumab PK model (Soehoel 2022)                                                                                                             |
| Xie_2019_agomelatine                | A semiphysiological population pharmacokinetic model of agomelatine and its metabolites in Chinese healthy volunteers                            |
| Zhu_2017_lebrikizumab               | Lebrikizumab PK model (Zhu 2017)                                                                                                                 |
| oncology_sdm_lobo_2002              | Signal transduction model for delayed concentration effects on cancer cell growth                                                                |
| oncology_xenograft_simeoni_2004     | Oncology tumor growth model in xenograft models                                                                                                  |
| tgi_no_sat_Koch                     | One compartment TGI model with with exponential tumor growth, without saturation.                                                                |
| tgi_no_sat_expo                     | One compartment TGI model with with exponential tumor growth, without saturation.                                                                |
| tgi_no_sat_linear                   | One compartment TGI model with with linear tumor growth, without saturation.                                                                     |
| tgi_no_sat_powerLaw                 | One compartment TGI model with with exponential tumor growth, without saturation.                                                                |
| tgi_sat_VonBertalanffy              | One compartment TGI model where tumor growth is limited by a loss term, with saturation.                                                         |
| tgi_sat_genVonBertalanffy           | One compartment TGI model where tumor growth is limited by a loss term, with saturation.                                                         |
| tgi_sat_logistic                    | One compartment TGI model with with exponential tumor growth that decelerates linearly, with saturation.                                         |
