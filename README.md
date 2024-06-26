# High K diet
This code is for the K+ regulation model used in [Stadt & Layton, A modeling analysis of whole-body potassium regulation on a high potassium diet: Proximal tubule and tubuloglomerular feedback effects, AJP-Regulatory, Integrative and Comparative Physiology, 2024](https://doi.org/10.1152/ajpregu.00283.2023). 

Supplementary materials availabe in S1_text.pdf



## Figure replication
These are the files use for making the figures in the manuscript. The simulation data is stored in the directory "MultiDaySim/" and the sensitivity analysis data is stored in "Sensitivity/". 

**plot_urine_highK_GFRPT.m** Fig 3.1 and 3.3

**plot_highK_GFRPT.m** Fig 3.2

**plot_etaPTKreab_change.m** Fig 3.4

**plot_alphaTGF.m** Fig 3.5

**plot_alphaTGF_etaPT_combine.m** Fig 3.6

**postprocess_localsensitivity.m** Fig 3.7

**postprocess_MA.m** Fig 3.8 and Fig 3.9



## Key files
- The simulations for the data in "MultiDaySim" are done by using **driver_TGFeffect_multiday.m**
- Model code **kreg_eqns.m**
- Parameter values are set using **set_params.m**
- Files used to conduct the sensitivity are in "SA". Use **localsensitivity.r** for local sensitivity analysis and **MorrisAnalysis_r100.r** for the Morris analysis. Postprocess results with **postprocess_localsensitivity.r** and **postprocess_MorrisAnalysis.r** to convert to .csv for the figures.


## Notes
- The global sensitivity analysis takes several hours to run.

