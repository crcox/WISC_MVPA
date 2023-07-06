MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
MATLABDIR ?= /localbuildsoftware/matlab/R2018b
MCC=$(MATLABDIR)/bin/mcc
MEX=$(MATLABDIR)/bin/mex
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
GLMNET=$(DEP)/glmnet
INCL= -I $(SRC) -I $(JSON) -I $(GLMNET)
.PHONEY: all clean-all clean-postbuild glmnet sdist

all: setup glmnet WISC_MVPA clean-postbuild

loadmodule:
	module load MATLAB/r2018b

setup: $(SRC) $(DEP)
$(SRC) $(DEP):
	tar xzf $(SRCTAR)

glmnet: $(GLMNET)/glmnetMex.mexa64
$(GLMNET)/glmnetMex.mexa64: $(GLMNET)/glmnetMex.F $(GLMNET)/GLMnet.f
	$(MEX) -fortran -outdir $(GLMNET) $^

WISC_MVPA: $(SRC)/WISC_MVPA.m
	$(MCC) -v $(MFLAGS) $(INCL) -o $@ $<

clean-postbuild:
	-rm *.dmr
	-rm mccExcludedFiles.log
	-rm readme.txt
	-rm run_WISC_MVPA.sh

sdist:
	tar czhf $(SRCTAR) src dependencies

clean-all:
	-rm WISC_MVPA
