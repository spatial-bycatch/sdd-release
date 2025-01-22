
.PHONY: all

MDLS = idd_exp idd_exp_reg idd_wei idd_wei_reg
SPPS = BBE SSI ETB SND

all:
	$(foreach x,$(MDLS),	\
		$(foreach y,$(SPPS),	\
			mkdir -p ./$(x)/$(y);	\
			cp -f ./$(x)/run.R ./$(x)/$(y)/run.R;	\
			cd ./$(x)/$(y);  Rscript run.R; cd ../..;	\
			rm ./$(x)/$(y)/run.R;	\
	))
	
test: 
	$(foreach x,$(MDLS),	\
		$(foreach y,$(SPPS),	\
			mkdir -p ./$(x)/$(y);	\
			cp -f ./$(x)/test.R ./$(x)/$(y)/test.R;	\
			cd ./$(x)/$(y);  Rscript test.R; cd ../..;	\
			rm ./$(x)/$(y)/test.R;	\
	))


