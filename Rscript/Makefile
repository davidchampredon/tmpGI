#####################################################
#####     
#####     GENERATION INTERVALS PROJECT
#####
#####     Create figures for manuscript
#####     Depends on C++ simulations run in
#####     ~/Dropbox/PROJECTS/Gillespie_SEmInR/
#####
#####     Created on 2015-07-10 by David Champredon
#####
#####################################################


EXPL_GI := explain_GIfwdbck_v6.R
SIM_VS_THEO := GI_Distribution_sim_vs_theory_4.R
MULTI_GI := GIfwdbck_multi_R0.R
FIT_GI = fit_model_bckGI_2.R

FCT_FILES := figures_ms.R GIfwdbck_multi_R0_FCT.R calc_theoretical_GI.R SEmInR_deterministic.R read_simul_GI_FCT.R fit_model_bckGI_FCT_2.R

DIR_CPP := ../Gillespie_SEmInR


###########################################################
### Re-create all figures (delete all existing!) 
###########################################################

allfig.out: $(EXPL_GI) $(SIM_VS_THEO) $(MULTI_GI) $(FCT_FILES) $(DIR_CPP)/OUT/*.out
	/bin/rm -f *.pdf
	Rscript $(EXPL_GI) &
	Rscript $(SIM_VS_THEO) &
	Rscript $(MULTI_GI) &
	Rscript $(FIT_GI) &
	cat $(DIR_CPP)/simul_param*.csv > $@


