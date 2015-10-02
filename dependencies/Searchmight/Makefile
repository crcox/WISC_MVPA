MATLABDIR ?= /usr/local/MATLAB/R2013b
INCL= -I$(MATLABDIR)/extern/include
all: searchmightGNB repmat

searchmightGNB:
	mex $(INCL) searchmightGNB.c

repmat: mexutil
	mex $(INCL) repmat.c mexutil.o

mexutil:
	mex $(INCL) -c mexutil.c

clean:
	rm searchmightGNB.mex* repmat.mex* mexutil.o
