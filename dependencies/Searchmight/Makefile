all: searchmightGNB repmat

searchmightGNB:
	mex searchmightGNB.c

repmat: mexutil
	mex repmat.c mexutil.o

mexutil:
	mex -c mexutil.c

clean:
	rm searchmightGNB.mex* repmat.mex* mexutil.o
