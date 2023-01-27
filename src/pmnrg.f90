! pmnrg.f90 - Module providing simple 3D and 2D multivariate random generators.
!
! -------------------------------------------------------------------------------------------------
!
! Copyright 2023 Patrizia Favaron
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
! documentation files (the "Software"), to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
! and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions
! of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
! TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!
! -------------------------------------------------------------------------------------------------
!
module pmnrg
    
    implicit none
    
    private
    
    ! Public interface
    public  :: MultiNormalGenerator
    
    ! Internal constants
    integer, parameter  :: C = 128
    real(8), parameter  :: R = 3.442619855899d0
    real(8), parameter  :: V = 9.91256303526217d-3
    
    ! Generator data type
    type MultiNormalGenerator
        private
        logical                     :: lInitialized = .false.
        real(8), dimension(0:C)     :: s_adZigX
        real(8), dimension(0:C-1)   :: s_adZigR
    contains
        procedure, public   :: get
        procedure, private  :: Ziggurat
    end type MultiNormalGenerator
    
    interface MultiNormalGenerator
        procedure   :: InitZiggurat
    end interface MultiNormalGenerator
    
contains

    ! Automatic type constructor for 'MultiNormalGenerator' (notice the not-that-unusual
    ! syntax: it's a regular function returning a MultiNormalGenerator, aliased to MultiNormalGenerator
    ! in the module declarations part: this is Fortran 2003 to say "constructor"!).
    !
    ! The implementation, object stuff aside, follows almost identically professor Doornik's.
    !
    type(MultiNormalGenerator) function InitZiggurat()
        
        ! Routine arguments
        ! --none--
        
        ! Locals
        integer :: i
        real(8) :: f
        
        ! Initialize tables used when generating random data
        f = exp(-0.5d0 * R**2)
        InitZiggurat % s_adZigX(0)  = V / f
        InitZiggurat % s_adZigX(1)  = R
        InitZiggurat % s_adZigX(C) = 0.d0
        do i = 2, C-1
            InitZiggurat % s_adZigX(i) = sqrt(-2.d0 * log(V / InitZiggurat % s_adZigX(i - 1) + f))
            f = exp(-0.5d0 * InitZiggurat % s_adZigX(i)**2)
        end do
        do i = 0, C-1
            InitZiggurat % s_adZigR(i) = InitZiggurat % s_adZigX(i + 1) / InitZiggurat % s_adZigX(i)
        end do
        InitZiggurat % lInitialized = .true.

    end function InitZiggurat


    ! Public member, providing multivariate random numbers: here is where the real work is made.
    ! Use of this function MUST be preceded by an assigment having form
    !
    !   tRandom = MultiNormalGenerator()
    !
    ! where tRandom has been declared as type(MultiNormalGenerator). This is to force the
    ! type's constructor to start: would this step not be made, then routine 'get' would terminate
    ! with error code 1.
    !
    ! Return codes:
    !
    !   0   Successful completion: multivariate normal deviates generated
    !
    !   1   Attempt made to use a non-initialized generator (assign object the default
    !       MultiNormalGenerator() before)
    !
    !   2   n = size(rvMean) <= 0 (valid values should be positive)
    !
    !   3   size(rmCov, dim=1) /= n (should be 'n')
    !
    !   4   size(rmCov, dim=2) /= n (should be 'n')
    !
    !   5   size(rmData, dim=2) /= n (should be 'n')
    !
    !   6   size(rmData, dim=1) <= 0 (should be positive)
    !
    !   7   matrix 'rmCov' is not symmetric
    !
    !   8   matrix 'rmCov' is not positive definite
    !
    function get(this, rvMean, rmCov, rmData) result(iRetCode)
        
        ! Routine arguments
        class(MultiNormalGenerator), intent(in) :: this     ! This function's class
        real(8), dimension(:), intent(in)       :: rvMean   ! 'n'-dimensioned vector containing means
        real(8), dimension(:,:), intent(in)     :: rmCov    ! 'n' x 'n' matrix containing covariances
        real(8), dimension(:,:), intent(out)    :: rmData   ! 'm' x 'n' matrix containing random deviates
                                                            ! (row-wise, one data vector per column)
        integer                                 :: iRetCode ! Return code
        
        ! Locals
        integer :: i
        integer :: j
        integer :: n
        integer :: k
        integer :: iErrCode
        real(8), dimension(:,:), allocatable    :: rmL
        
        ! Steering constants
        real(8), parameter  :: TOL = 1.0d-6
        
        ! Assume success (will falsify on failure)
        iRetCode = 0
        
        ! Get and check dimensions
        if(.not. this % lInitialized) then
            iRetCode = 1
            return
        end if
        n = size(rvMean)
        if(n <= 0) then
            iRetCode = 2
            return
        end if
        if(n /= size(rmCov, dim=1)) then
            iRetCode = 3
            return
        end if
        if(n /= size(rmCov, dim=2)) then
            iRetCode = 4
            return
        end if
        if(size(rmData, dim=2) /= n) then
            iRetCode = 5
            return
        end if
        if(size(rmData, dim=1) <= 0) then
            iRetCode = 6
            return
        end if
        do i = 1, n-1
            do j = i+1, n
                if(abs(rmCov(i,j) - rmCov(j,i)) > TOL) then
                    iRetCode = 7
                    return
                end if
            end do
        end do
        
        ! Compute Cholesky decomposition of the covariance matrix
        allocate(rmL(n,n))
        rmL = rmCov
        iErrCode = Cholesky(rmL)
        if(iErrCode /= 0) then
            iRetCode = 8
            return
        end if
        
        ! Set all upper elements of L to 0 (they are not, because of the
        ! way "Cholesky" function is written)
        do i = 1, n-1
            rmL(i, (i+1):n) = 0.d0
        end do
        
        ! Extract random sample with 0 mean and unit standard deviation
        do i = 1, n
            do j = 1, k
                rmData(j,i) = this % Ziggurat()
            end do
        end do
        
        ! Apply transformation
        do j = 1,k
            rmData(j,:) = rvMean + matmul(rmL, rmData(j,:))
        end do
        
        ! Leave
        deallocate(rmL)
        
    end function get
    

    ! Single N(0,1) random generator
    function Ziggurat(this) result(rNormalDeviate)
        
        ! Routine arguments
        class(MultiNormalGenerator), intent(in) :: this
        real(8)                                 :: rNormalDeviate
        
        ! Locals
        integer     :: i
        real(8)     :: x
        real(8)     :: u
        real(8)     :: f0
        real(8)     :: f1
        real(8)     :: rRanU
        integer     :: iRanU
        
        ! Extract random normal deviate through Ziggurat method
        do
            
            call random_number(rRanU)
            u = 2.d0 * rRanU - 1.d0
            call random_number(rRanU)
            iRanU = nint(rRanU * 10000.d0)
            i = iand(iRanU, Z'7F')
            if(abs(u) < this % s_adZigR(i)) then
                rNormalDeviate = u * this % s_adZigX(i)
                return
            end if
            if(i == 0) then
                rNormalDeviate = DRanNormalTail(R, u < 0)
                return
            end if
            x = u * this % s_adZigX(i)
            f0 = exp(-0.5d0 * (this % s_adZigX(i)**2 - x**2))
            f1 = exp(-0.5d0 * (this % s_adZigX(i+1)**2 - x**2))
            call random_number(rRanU)
            if(f1 + rRanU * (f0 - f1) < 1.0d0) then
                rNormalDeviate = x
                return
            end if
            
        end do
        
    end function Ziggurat

    ! **********************
    ! * Internal functions *
    ! **********************
    
    ! In-place Cholesky decomposition according to Algorithm 4.2.1 of
    ! Golub - van Loan, ed.3, page 144
    function Cholesky(A) result(iRetCode)
        
        ! Routine arguments
        real(8), dimension(:,:), intent(inout)  :: A
        integer                                 :: iRetCode
        
        ! Locals
        integer     :: i, j, n
        real(8)     :: rDenom
        
        ! Assume success (will falsify on failure)
        iRetCode = 0
        
        ! Check dimensions: only non-null square matrices allowed
        n = size(A, dim=1)
        if(n <= 0) then
            iRetCode = 1
            return
        end if
        if(n /= size(A, dim=2)) then
            iRetCode = 2
            return
        end if
        
        ! Perform the Cholesky decomposition
        do j = 1, n
            if(j > 1) then
                do i = j, n
                    A(i,j) = A(i,j) - dot_product(A(i,1:(j-1)), A(j,1:(j-1)))
                end do
            end if
            rDenom = A(j,j)
            if(abs(rDenom) < 1.0d-12) then
                ! Matrix is not definite positive!
                iRetCode = 3
                return
            end if
            A(j:n,j) = A(j:n,j)/sqrt(rDenom)
        end do
            
    end function Cholesky
    
    
    ! Auxiliary function, used by "Ziggurat"
    function DRanNormalTail(rMin, lNegative) result(rSample)
        
        ! Routine arguments
        real(8), intent(in)     :: rMin
        logical, intent(in)     :: lNegative
        real(8)                 :: rSample
        
        ! Locals
        real(8)     :: rX
        real(8)     :: rY
        real(8)     :: rU
        real(8)     :: rV

        ! Compute the desired quantity
        do
            call random_number(rU)
            call random_number(rV)
            rX = log(rU) / rMin
            rY = log(rV)
            if(-2.d0*rY >= rX**2) exit
        end do
        
        ! Yield result
        if(lNegative) then
            rSample = rX - rMin
        else
            rSample = rMin - rX
        end if
        
    end function DRanNormalTail

end module pmnrg
