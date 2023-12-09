# High K diet
This code is for the K+ regulation model under a high K+ diet. 

Please cite [Stadt & Layton, PLOS Comp Bio, 2022](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010607)

## Figure replication
These are the files use for making the figures in the manuscript. The simulation data is stored in the directory "MultiDaySim/" and the sensitivity analysis data is sotred in "Sensitivity/". This is used to make the plots.
**plot_urine_highK_GFRPT.m** Fig 3.1 and 3.3

**plot_highK_GFRPT.m** Fig 3.2 and A.1

**plot_etaPTKreab_change.m** Fig 3.4

**plot_alphaTGF.m** Fig 3.5

**postprocess_MA.m** Fig 3.6 and 3.7

**plot_MKX.m** Fig 3.8 and 3.9

**postprocess_localsensitivity.m** Fig A.2

## Key files
- The simulations for the data in "MultiDaySim" are done by using **driver_TGFeffect_multiday.m**
- Model code **kreg_eqns.m**
- Parameter values are set using **set_params.m**
- Files used to conduct the sensitivity are in "SA". Use **localsensitivity.r** for local sensitivity analysis and **MorrisAnalysis_r100.r** for the Morris Analysis. Postprocess results with **postprocess_localsensitivity** and **postprocess_MorrisAnalysis.r**

