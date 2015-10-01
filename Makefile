MCC=/usr/local/MATLAB/R2013b/bin/mcc
MEX=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
SEARCHMIGHT=$(DEP)/Searchmight
INCL= -I$(SRC) -I$(JSON) -I$(SEARCHMIGHT) -I$(SEARCHMIGHT)/private \
      -I$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/cvx \
      -I$(SEARCHMIGHT)/CoreToolbox/ExternalPackages.Linux_x86_64/libsvm
.PHONEY: clean clean-all all source_code.tar.gz extract

all: WholeBrain_MVPA binaries.tar.gz

sdist:
	tar czhf source_code.tar.gz src dependencies

extract:
	-mkdir source_code/
	tar xzf source_code.tar.gz -C ./source_code/

WholeBrain_MVPA: $(SRCDIR)/WholeBrain_MVPA.m
	$(MCC) $(MFLAGS) $(IDIRS) -o $@ WholeBrain_MVPA.m

binaries.tar.gz: WholeBrain_MVPA run_WholeBrain_MVPA.sh
	mkdir bin/
	mv WholeBrain_MVPA WholeBrain_MVPA.sh bin/
	tar czvf $@ bin/

clean:
	-rm *.dmr
	-rm _condor_std???
	-rm readme.txt
	-rm mccExcludedFiles.log

clean-all: binaries.tar.gz
	rm -rf src/ bin/
