
.PHONY: all

MDLS = sdd_exp sdd_exp_reg sdd_wei sdd_wei_reg
SPPS = ETB BBE SND SSI
    
all: 
	$(foreach x,$(MDLS),	\
		$(foreach y,$(SPPS),	\
			mkdir -p ./$(x)/$(y);	\
			cp -f ./$(x)/run.R ./$(x)/$(y)/run.R;	\
			cd ./$(x)/$(y);  Rscript run.R; cd ../..;	\
			rm ./$(x)/$(y)/run.R;	\
	))

