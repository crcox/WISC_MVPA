MCC=/usr/local/MATLAB/R2013b/bin/mcc
MEX=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
SEARCHMIGHT=$(DEP)/Searchmight
INCL= -I $(SRC) -I $(JSON) -I $(SEARCHMIGHT) -I $(SEARCHMIGHT)/CoreToolbox \
      -I $(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/cvx \
      -I $(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/libsvm
.PHONEY: all clean-postbuild sdist

all: setup WholeBrain_MVPA clean-postbuild

setup:
	tar xzvf $(SRCTAR)

WholeBrain_MVPA: $(SRC)/WholeBrain_MVPA.m
	$(MCC) $(MFLAGS) $(INCL) -o $@ $^

clean-postbuild:
	rm *.dmr
	rm mccExcludedFiles.log
	rm readme.txt
	rm run_WholeBrain_MVPA.sh

sdist:
	tar czhf $(SRCTAR) src dependencies
