
.PHONY: all

MDLS = sdd_exp_hyb sdd_exp_hyb_reg sdd_wei_hyb sdd_wei_hyb_reg
SPPS = BBE SSI ETB SND
    
all: 
	$(foreach x,$(MDLS),	\
		$(foreach y,$(SPPS),	\
			mkdir -p ./$(x)/$(y);	\
			cp -f ./$(x)/run.R ./$(x)/$(y)/run.R;	\
			cd ./$(x)/$(y);  Rscript run.R; cd ../..;	\
			rm ./$(x)/$(y)/run.R;	\
	))

