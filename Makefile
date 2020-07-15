export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
export MCC="mcc"
export MEX="mex"
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
GLMNET=$(DEP)/glmnet
UTILS=$(DEP)/mri_coordinate_tools
SEARCHMIGHT=$(DEP)/Searchmight
LIBSVM=$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/libsvm
CVX=$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/cvx
INCL= -I $(SRC) -I $(JSON) -I $(UTILS)
#-I $(GLMNET) \
#      -I $(SEARCHMIGHT) \
#      -I $(SEARCHMIGHT)/CoreToolbox \
#      -I $(LIBSVM) \
#      -I $(CVX) \
#      -I $(CVX)/builtins \
#      -I $(CVX)/commands \
#      -I $(CVX)/functions \
#      -I $(CVX)/keywords \
#      -I $(CVX)/lib \
#      -I $(CVX)/matlab6 \
#      -I $(CVX)/sdpt3 \
#      -I $(CVX)/sdpt3/Examples \
#      -I $(CVX)/sdpt3/HSDSolver \
#      -I $(CVX)/sdpt3/Linsysolver \
#      -I $(CVX)/sdpt3/Linsysolver/spchol \
#      -I $(CVX)/sdpt3/Solver \
#      -I $(CVX)/sdpt3/Solver/Mexfun \
#      -I $(CVX)/sdpt3/dimacs \
#      -I $(CVX)/sdpt3/sdplib \
#      -I $(CVX)/sedumi \
#      -I $(CVX)/sedumi/conversion \
#      -I $(CVX)/sedumi/doc \
#      -I $(CVX)/sedumi/mexw32 \
#      -I $(CVX)/sets \
#      -I $(CVX)/structures
.PHONEY: all clean-all clean-libsvm clean-Searchmight clean-postbuild cvx glmnet libsvm Searchmight sdist

all: setup glmnet WISC_MVPA clean-postbuild

loadmodule:
	module load MATLAB/r2018b

setup: $(SRC) $(DEP) $(UTILS)
$(SRC) $(DEP):
	tar xzvf $(SRCTAR)
$(UTILS): $(DEP)
	git clone https://github.com/crcox/DefineCommonGrid.git $(DEP)/DefineCommonGrid

libsvm: $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/svmpredict.mexa64 $(LIBSVM)/svmtrain.mexa64
$(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/svmpredict.mexa64 $(LIBSVM)/svmtrain.mexa64:
	$(MAKE) -C $(LIBSVM)

cvx:
	$(MAKE) -C $(CVX)

Searchmight: $(SEARCHMIGHT)/searchmightGNB.mexa64
$(SEARCHMIGHT)/searchmightGNB.mexa64: $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/read_sparse.mexa64
	$(MAKE) -C $(SEARCHMIGHT)
	mv $(SEARCHMIGHT)/repmat.mex* $(SEARCHMIGHT)/private/

glmnet: $(GLMNET)/glmnetMex.mexa64
$(GLMNET)/glmnetMex.mexa64: $(GLMNET)/glmnetMex.F $(GLMNET)/GLMnet.f
	$(MEX) -fortran -outdir $(GLMNET) $^

# $(SEARCHMIGHT)/searchmightGNB.mexa64
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
	-rm $(SEARCHMIGHT)/searchmightGNB.mexa64
	-rm $(LIBSVM)/read_sparse.mexa64
	-rm $(LIBSVM)/svmpredict.mexa64
	-rm $(LIBSVM)/svmtrain.mexa64
