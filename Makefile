
SRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSEngine.c smacofSSMajorize.c smacofSSMonotone.c

%.o: %.c smacofSS.h
	clang -c $@

shlib: smacofSS.h $(SRC)
	R CMD SHLIB -o smacofSS.so $(SRC)

clean:
	rm -f *.o *.so

