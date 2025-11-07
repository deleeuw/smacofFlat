
SRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSEngine.c smacofSSMajorize.c smacofSSMonotone.c

%.o: %.c smacofSS.h
	clang -c $@

shlib: $(SRC)
	R CMD SHLIB -o smacofSS.so $(SRC)

walker: walker.c $(SRC)
	clang -o walker walker.c $(SRC)

runner: runner.c $(SRC)
	clang -o runner runner.c $(SRC)

clean:
	rm -f runner walker *.o *.so

