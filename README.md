# pmnrg - Patrizia's minimalistic multivariate normal random generator

## Purpose and nature of 'pmnrg'

'pmnrg.f90' is a self-standing Fortran 2003 module for generating multivariate normal random deviates, given the mean vector and covariance matrix.

## Usage and functions synopsis

Access to 'pmnrg' is made through a use clause,

    use pmnrg
  
Once "used", the module gives access to a data type, 'MultiNormalGenerator', which can be instantiated in user code as in

    type(MultiNormlGenerator) :: tRand
    
Please remember to assign the type value after declaration, to ensure the type-bound constructor is invoked:

    tRand = MultiNormalGenerator()
    
Not doing so will cause the 'get' function to terminate in error.

The 'pmnrg' module contains an only public function, invoked as

        iRetCode = tRand % get(rvMean, rmCov, rmData)
        
where 'rvMean' is a length-N vector containing the means, 'rmCov' is the order-N covariance matrix, and 'rmData' a NxK matrix containing on successful exit the desired random deviates.

The value of N must be positive; the case N=1 is admitted, meaning univariate data series.

The 'get' function yields an integer return code, with the following values:

0.  Successful completion: multivariate normal deviates generated.
1.  Attempt made to use a non-initialized generator.
2.  N = size(rvMean) <= 0 (valid values should be positive).
3.  size(rmCov, dim=1) /= N (should be 'N').
4.  size(rmCov, dim=2) /= N (should be 'N').
5.  size(rmData, dim=1) /= N (should be 'N').
6.  size(rmData, dim=2) <= 0 (should be positive).
7.  Matrix 'rmCov' is not symmetric.
8.  Matrix 'rmCov' is not positive definite.

## Credits

As a module, "pmnrg" univariate normal generator is a translation of the ZIGNOR code by prof. J.A. Doornik, described in "An Improved Ziggurat Method to Generate Normal Random Samples", as I've found it in https://www.doornik.com/research/ziggurat.pdf.

Translation was performed quite "aggressively", employing a bit of Fortran 2003 and using some object-oriented programming in sake of ease of use and clarity. Nonetheless, the important variable names have been retained to permit readers to feel familiar with prof. Doornik's report.

## Rationale of having selected the Ziggurat method

Choice of Ziggurat method instead of, say, the selection procedure shown in D.E. Knuth, "The Art of Computer Programming", volume 2, "Seminumerical Algorithms", is due to Ziggurat's higher computational efficiency.

This "may" matter in the use case of interest to Patrizia, where possibly millions or billions deviates need to be generated: a tiny time saving may have then a deep impact on overall feasibility of the process at hand.

## Theory of operation

Function 'Ziggurat' is a univariate generator, and is not directly accessible to users within 'pmnrg' module. The really accessible generator is the multi-variate one, 'get', explained above.

The implementation is straightforward. Let's assume 'M' represents the means vector, and 'C' the covariance matrix ('M' and 'C' should be known before invoking the generator). Then the multivariate normal deviates are obtained as

    r = M + u*L
    
where 'u' is a vector of normally distributed values with mean 0 and standard deviation 1, and 'L' the lower matrix of C Cholesky decomposition (L Transpose(L) = C). Then it can be shown (prof. Doornik does) that 'r' is normally distributed with means 'M' and covariances 'C'.

The covariance matrix 'C' is positive definite and symmetric, so the Cholesky factorization 'L' is guaranteed to exist.

Cholesky decomposition has been coded using Algorithm 4.2.1 in G.H. Golub, C.F. van Loan, "Matrix Computations", 3rd edition, John Hopkins University Press, page 144.

In a sense, computing the Cholesky decomposition is the matrix equivalent of extracting a square root. The previous formula can then be seen as a matrix analogous of the univariate scaling

    r = mu + u*sigma
    
with 'sigma' the standard deviation, that is, the square root of variance.

## Advantages of Cholesky decomposition

The Cholesky decomposition is just one of the many decompositions available, the LU decomposition among these.

Before to begin coding, I've done some evaluations of competing methods, and Cholesky decomposition emerged as the clear winner - at least among the methods I've heard of when a college girl.

Indeed using Cholesky decomposition makes a lot of sense: as mentioned in the previous section the covariance matrix is positive definite and symmetric, so the Cholesky decomposition exists by necessity. Cholesky decomposition is also reputed for its numerical stability, and does not require pivoting steps.

The lack of pivoting steps has the interesting consequence that matrix rows are not rearranged as they would be using e.g. LU decomposition. I've considered the need to rearrange decomposition results before to rescale the deviates quite _obscuring_, a very good reason to avoid it.

Last, but not least, Cholesky decomposition is very efficient compared to other methods. This is not a critical point, the matrix decomposition being computed only once for each generation step. But it does not harm...

## Reference use cases

Module 'pmnrg' is purpose-built, to support the realization of Lagrangian particle dispersion models in two and three dimensions. That is to say, I have been quite brutal in my code writing, and paid no attention to large matrices. A larger scale implementation would require consideration of roundoff, maybe not severe because of the method stability, yet to consider: use of professionally-developed libraries like LAPACK would be advisable in case. And, I did not do...

That said, code is not formally restricted to 2x2 and 3x3 covariance matrices: it may be used in any case, at least in principle. I made no effort to test with cases larger than 3x3 however: use at your own risk!

## Why MIT open-source license?

A final word about open-source licensing. I decided to use the rather permissive MIT license in this case, because the code you are reading was written with me in the role of thesis supervisor. As such, it was intended for my cubs' and anyone else's use. I cannot of course guarantee my code is absolutely correct (I did my best to write and test it correctly), and assume anyone wishing to use it has enough intelligence and preparation to understand whether it fulfils their own necessities.

On the other side, I admit not being obsessed with open-source "purity" and similar things, looking quite ideological to my simple and down-to-earth mind. I'm sure I've missed something in the overall debate, and would be sincerely grateful to anyone taking their time to illuminate it. In the same time, I objectively have no need to be remembered as a developer-of-the-past when I'll be passed out. And, very strongly feel that writing code is (or _should_ be) an act of care, with no degree of contriction or any other form of violence built-in. Besides, there is more than enough of it in this rough time.

For these reasons I decided to _not_ use any variant of GPL license. That's positively fully intentional.

Thanks,
Patrizia

