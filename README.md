# smacofFlat
The smacof at https://www.github.com/deleeuw/smacofFlat combines R and C, trying
to use the strengths of both languages. It also

- minimizes the use of matrices and uses a flat data structures with indices
- has the smacof upgrade formula in C
- uses Busing's monotone function in C for monotone regression
- implements the three tie approaches for monotone regression, also in C
- is a monolithic C program with R frontend for data entry
- smacofSS is for square symmetric MDS problems weighted/unweighted and ratio/ordinal
- it consists of four separate R programs for these four options (more will follow)
- the default initial is torgerson (if missing data, impute using average dissimilarity)
- the repository has the R and C code but not the compiled code
- I have been using it for ratio/ordinal versions of McGee Elastic and Sammon Mapping
- creates the MDS data structure (a matrix) from Delta or (Delta, W) and back
- only entries below the diagonal are used
- in the MDS data structure dissimilarities are already sorted
- smacofSS uses a single .C() call, so the C code does not use the R API
- the C code is supposed to be portable
- 5 to 30 times as fast as smacofSym from the smacof package
- R routines for data manipulation and plotting are included 
