MCC=/usr/local/MATLAB/R2013b/bin/mcc
MEX=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRCDIR=src

.PHONEY: clean clean-all all

all: SOSLasso summarize_jobs binaries.tar.gz

SOSLasso: $(SRCDIR)/runSOSLasso.m
	$(MCC) $(MFLAGS) $(SRCDIR)/runSOSLasso.m -o $@

summarize_jobs: $(SRCDIR)/summarize_jobs.m
	$(MCC) $(MFLAGS) $(SRCDIR)/runSOSLasso.m -o $@

binaries.tar.gz: SOSLasso run_SOSLasso.sh summarize_jobs run_summarize_jobs.sh
	mkdir bin/
	mv SOSLasso run_SOSLasso.sh bin/
	mv summarize_jobs run_summarize_jobs.sh bin/
	tar czvf $@ bin/

clean:
	-rm *.dmr
	-rm _condor_std???
	-rm readme.txt
	-rm mccExcludedFiles.log

clean-all: binaries.tar.gz
	rm -rf src/ bin/
