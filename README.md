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

## Credits

As a module, "pmnrg" univariate normal generator is a Fortran translation of the ZIGNOR code by J.A. Doornik, in "An Improved Ziggurat Method to Generate Normal Random Samples", as I've found it in https://www.doornik.com/research/ziggurat.pdf.

Translation was performed quite "aggressively", employing a bit of Fortran 2003 and more, doing my best to ease the users' life up a little bit.

Choice of Ziggurat method instead of, say, the selection procedure shown in D.E. Knuth, "The Art of Computer Programming", volume 2, "Seminumerical Algorithms", is due to Ziggurat's higher computational efficiency. This "may" matter in the use case of interest to Patrizia, where possibly millions or billions deviates need to be generated: a tiny time saving may have then a deep impact on overall
feasibility of the process at hand.

What's really of use is not the unifariate generator, but the multi-variate. This has been coded in a straightforward way, using direct Cholesky decomposition, the latter calculated as shown in G.H. Golub, C.F. van Loan, "Matrix Computations", 3rd edition, John Hopkins University Press. Using Cholesky decomposition makes a lot of sense, the covariance matrix being positive definite, in which case the Cholesky decomposition is guaranteed to exist. Cholesky decomposition is also reputed for its numerical stability, and does not require pivoting steps.

Incidentally, "pmnrg" means "Patti's Minimalistic Normal Random Generator". It is purpose-built, supporting the realization of Lagrangian particle dispersion models in two and three dimensions. That is to say, I have been quite brutal in my code writing, and used no libraries - LAPACK being a first-idea-jumping-to-mind. No libraries means no dependencies, and easier to deploy code base...

I could have used LU decomposition instead, but in that case pivoting and permutators would have had made my code trickier, harder to understand, and (marginally) less computationally efficient. I'm decidedly _not_ a genius, and anyone like me will appreciate my quest for simplicity. Besides, using LU decomposition would had demanded a tricky nitty-gritty code I sincerely have no time to conceive and debug: use of library code would have then be absolutely advisable in that case as an error-preventing
approach.

That said, code is not formally restricted to 2x2 and 3x3 covariance matrices: it may be used in any case, at least in principle. I made no effort to test with cases larger than 3x3 however: use at your own risk!

## Why MIT open-source license?

A final word about open-source licensing. I decided to use the rather permissive MIT license in this case, because the code you are reading was written with me in the role of thesis supervisor. As such, it was intended for my cubs' and anyone else's use. I cannot of course guarantee my code is absolutely correct (I did my best to write and test it correctly), and assume anyone wishing to use it has enough intelligence and preparation to understand whether it fulfils their own necessities.

On the other side, I admit not being obsessed with open-source "purity" and similar things, looking quite ideological to my simple and down-to-earth mind. I'm sure I've missed something in the overall debate, and would be sincerely grateful to anyone taking their time to illuminate it. In the same time, I objectively have no need to be remembered as a developer-of-the-past when I'll be passed out. And, very strongly feel that writing code is (or _should_ be) an act of care, with no degree of contriction or any other form of violence built-in. Besides, there is more than enough of it in this rough time.

For these reasons I decided to _not_ use any variant of GPL license. That's positively fully intentional.

Thanks,
Patrizia

