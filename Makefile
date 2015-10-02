export MATLABDIR = /usr/local/MATLAB/R2013b
export MCC=/usr/local/MATLAB/R2013b/bin/mcc
export MEX=/usr/local/MATLAB/R2013b/bin/mex
export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
UTILS=$(DEP)/DefineCommonGrid
SEARCHMIGHT=$(DEP)/Searchmight
LIBSVM=$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/libsvm
CVX=$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/cvx
INCL= -I $(SRC) -I $(JSON) -I $(UTILS) \
      -I $(SEARCHMIGHT) \
      -I $(SEARCHMIGHT)/CoreToolbox \
      -I $(LIBSVM) \
      -I $(CVX) \
      -I $(CVX)/builtins \
      -I $(CVX)/commands \
      -I $(CVX)/functions \
      -I $(CVX)/keywords  \
      -I $(CVX)/lib \
      -I $(CVX)/matlab6 \
      -I $(CVX)/sdpt3 \
      -I $(CVX)/sedumi \
      -I $(CVX)/sets \
      -I $(CVX)/structures
.PHONEY: all clean-all clean-libsvm clean-Searchmight clean-postbuild libsvm Searchmight sdist

all: setup libsvm Searchmight WholeBrain_MVPA clean-postbuild

setup: $(SRC) $(DEP) $(UTILS)
$(SRC) $(DEP):
	tar xzvf $(SRCTAR)
$(UTILS): $(DEP)
	git clone https://github.com/crcox/DefineCommonGrid.git $(DEP)/DefineCommonGrid

libsvm: $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/svmpredict.mexa64 $(LIBSVM)/svmtrain.mexa64
$(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/svmpredict.mexa64 $(LIBSVM)/svmtrain.mexa64:
	$(MAKE) -C $(LIBSVM)

Searchmight: $(SEARCHMIGHT)/searchmightGNB.mexa64
$(SEARCHMIGHT)/searchmightGNB.mexa64: $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/read_sparse.mexa64 $(LIBSVM)/read_sparse.mexa64
	$(MAKE) -C $(SEARCHMIGHT)
	mv repmat.mex* private/

WholeBrain_MVPA: $(SRC)/WholeBrain_MVPA.m $(SEARCHMIGHT)/searchmightGNB.mexa64
	$(MCC) $(MFLAGS) $(INCL) -o $@ $<

clean-postbuild:
	rm *.dmr
	rm mccExcludedFiles.log
	rm readme.txt
	rm run_WholeBrain_MVPA.sh

sdist:
	tar czhf $(SRCTAR) src dependencies

clean-all:
	-rm WholeBrain_MVPA
	-rm $(SEARCHMIGHT)/searchmightGNB.mexa64
	-rm $(LIBSVM)/read_sparse.mexa64
	-rm $(LIBSVM)/svmpredict.mexa64
	-rm $(LIBSVM)/svmtrain.mexa64
