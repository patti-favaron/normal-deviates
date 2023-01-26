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
        
where 'rvMean' is a length-N vector containing the means, 'rmCov' is the order-N covariance matrix, and 'rmData' a KxN matrix containing on successful exit the desired random deviates.

The value of N must be positive; the case N=1 is admitted, meaning univariate data series. K must be positive, with typical values 2 and 3 in the author's application.

The 'get' function yields an integer return code, with the following values:

0.  Successful completion: multivariate normal deviates generated.
1.  Attempt made to use a non-initialized generator.
2.  N = size(rvMean) <= 0 (valid values should be positive).
3.  size(rmCov, dim=1) /= N (should be 'N').
4.  size(rmCov, dim=2) /= N (should be 'N').
5.  size(rmData, dim=2) /= N (should be 'N').
6.  size(rmData, dim=1) <= 0 (should be positive).
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

## Why the MIT open-source license, and not something else?

A final word about open-source licensing. I decided to use the permissive MIT license in this case, because the code you are reading was written with me in the role of thesis co-supervisor. As such, this code was intended for my cubs' and anyone else's use.

Distributing it as open-source was then out of question.

But then, _which_ open-source license?

Possibilities abound, some very popular as the GPL.

I'm not a lawyer, and can't really appreciate their subtleties in full. Yet I made some research, but before I tried to clarify what is the final purpose of this work (the same considerations apply to many others I've done).

As far as I've understood, licenses like the GPL contain an obligation. Any user of GPL licensed code should release their work under the same GPL license. Commercial use is not allowed. I realize this may be a push to enlarge the open-source code base more and more.

But on the other side, I feel at unease forcing others to adopt any method (and a license is just a special case) I could like. I can guess this code is far from originality, but yet it's useful. And usefulness is its very reason to be: surely I will never be remembered as one of the "coders of the past" when I'll pass out, nor imagining to become feels relevant to my being.

Quite on the contrary, a code is an expression of a will in a specific moment (ultimately of love). Is, in itself, a form of care. Or, I wish it to be.

And so, building constrictions into it would be, well, _limiting_. One other crippling forms of little violence, one more attempt to control others' lifes and possibilities. This is definitely not what _I_ would wish to a son or a daughter, would I happen to have one. And sure, the World abounds in little and big violences: adding one more would have had been decidedly not an amelioration.

Reading the text of the GPL I got the impression it's imbued of an idea of open-source "purity", with quite ad ideological flavor to my simple and down-to-earth tastes. I'm sure I've missed something in the overall debate, and would be sincerely grateful to anyone taking their time to illuminate it. Yet I understand licensing this work under the GPL would not be in my cubs best interest.

For these reasons I decided to _not_ use any variant of GPL license, and opt for the MIT's.

That's positively fully intentional, and longly thought.

Thanks,
Patrizia

