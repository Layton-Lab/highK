# Set up parameter ranges for Morris Analysis

testpars <- c("kgut",
                "GFR_base",
                "eta_ptKreab_base",
                "eta_LoHKreab" ,
                "eta_ptKreab_highK",
                "dtKsec_eq",
                "A_dtKsec",
                "B_dtKsec",
                "alpha_TGF",
                "A_cdKreab",
                "ALD_eq",
                "m_K_ALDO",
                "FF",
                "A_insulin",
                "alpha_al",
                "beta_al")

parsinf <- c(0.5 * p$kgut,
            0.095, # GFR_base (95/1000)
            0.36, # eta_ptKreab_base
            0.20, # eta_LoHKreab
            0.36, # eta_ptKreab_highK
            0.75 * p$dtKsec_eq, 
            0.5 * p$A_dtKsec,
            0.5 * p$B_dtKsec,
            0.5 * p$alpha_TGF,
            0.6 * p$A_cdKreab,
            0.5 * p$ALD_eq,
            0.5 * p$m_K_ALDO,
            0.5 * p$FF,
            0.5 * p$A_insulin,
            0.5 * p$alpha_al,
            0.5 * p$beta_al
            )


parssup <- c(1.5 * p$kgut,
            0.130, # GFR_base (130/1000)
            0.67, # eta_ptKreab_base
            0.30, # eta_LoH-Kreab
            0.67, # eta_ptKreab_highK, Up to the largest eta_ptKreab_base
            1.75 * p$dtKsec_eq, # dt K sec should have a larger range
            1.5 * p$A_dtKsec,
            1.5 * p$B_dtKsec,
            1.5 * p$alpha_TGF,
            1.5 * p$A_cdKreab,
            1.5 * p$ALD_eq,
            1.5 * p$m_K_ALDO,
            1.5 * p$FF,
            1.5 * p$A_insulin,
            1.5 * p$alpha_al,
            1.5 * p$beta_al
            )