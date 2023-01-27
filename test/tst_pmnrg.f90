! test_pmnrg.f90 - Test driver for 'pmnrg' module
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
program t_pmnrg
    
    use pmnrg
    
    implicit none
    
    ! Local variables
    real(8), dimension(:), allocatable      :: rvMean
    real(8), dimension(:,:), allocatable    :: rmData
    real(8), dimension(:,:), allocatable    :: rmCov
    integer                                 :: iRetCode
    integer                                 :: i, j
    type(MultiNormalGenerator)              :: tRand
    
    ! Test 1: non-initialized value
    print *, 'Test 1 - Non-initialized generator'
    allocate(rvMean(3), rmCov(3,3), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 1)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Assign a sensible initial value to the generator (this
    ! actually causes its constructor to execute)
    tRand = MultiNormalGenerator()
    
    ! Test 2: null means vector
    print *, 'Test 2 - Null means vector'
    allocate(rvMean(0), rmCov(3,3), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 2)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 3: non-square covariance matrix
    print *, 'Test 3 - Non conformable covariances matrix, case 1'
    allocate(rvMean(3), rmCov(3,4), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 4)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 4: non-square covariance matrix
    print *, 'Test 4 - Non conformable covariances matrix, case 2'
    allocate(rvMean(3), rmCov(3,2), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 4)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 5: non-square covariance matrix
    print *, 'Test 5 - Non conformable covariances matrix, case 3'
    allocate(rvMean(3), rmCov(4,3), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 3)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 6: non-square covariance matrix
    print *, 'Test 6 - Non conformable covariances matrix, case 4'
    allocate(rvMean(3), rmCov(2,3), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 3)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 7: non-square covariance matrix
    print *, 'Test 7 - Non conformable covariances matrix, case 5'
    allocate(rvMean(3), rmCov(2,6), rmData(1024,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 3)'
    deallocate(rvMean, rmCov, rmData)
    
    ! Test 8: non-conformable data matrix
    print *, 'Test 8 - Non conformable data matrix, case 1'
    allocate(rvMean(3), rmCov(3,3), rmData(1024,2))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 5)'
    deallocate(rvMean, rmCov, rmData)

    ! Test 9: non-conformable data matrix
    print *, 'Test 9 - Non conformable data matrix, case 2'
    allocate(rvMean(3), rmCov(3,3), rmData(1024,4))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 5)'
    deallocate(rvMean, rmCov, rmData)

    ! Test 10: null data vectors
    print *, 'Test 10 - Non conformable data matrix, case 2'
    allocate(rvMean(3), rmCov(3,3), rmData(0,3))
    rvMean = 0.0d0
    rmCov  = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 6)'
    deallocate(rvMean, rmCov, rmData)

    ! Test 11: non-symmetric matrix
    print *, 'Test 11 - Nonsymmetric covariance matrix'
    allocate(rvMean(3), rmCov(3,3), rmData(1024,3))
    call random_number(rmData)
    do i = 1, 3
        rvMean(i) = sum(rmData(i,:)) / 1024.0d0
        do j = i, 3
            rmCov(i,j) = sum((rmData(:,i)-rvMean(i))*(rmData(:,j)-rvMean(j)))/1024.0d0
            if(i /= j) rmCov(j,i) = rmCov(i,j)
        end do
    end do
    rmCov(1,2) = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 7)'
    deallocate(rvMean, rmCov, rmData)

    ! Test 12: non positive-definite matrix
    print *, 'Test 12 - Non positive-definite covaiance matrix'
    allocate(rvMean(3), rmCov(3,3), rmData(1024,3))
    call random_number(rmData)
    do i = 1, 3
        rvMean(i) = sum(rmData(i,:)) / 1024.0d0
        do j = i, 3
            rmCov(i,j) = sum((rmData(:,i)-rvMean(i))*(rmData(:,j)-rvMean(j)))/1024.0d0
            if(i /= j) rmCov(j,i) = rmCov(i,j)
        end do
    end do
    rmCov(2,2) = 0.0d0
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 8)'
    deallocate(rvMean, rmCov, rmData)

    ! Test 13: non positive-definite matrix
    print *, 'Test 13 - Correct data'
    allocate(rvMean(3), rmCov(3,3), rmData(1024*1024,3))
    call random_number(rmData)
    do i = 1, 3
        rvMean(i) = sum(rmData(:,i)) / (1024.0d0 * 1024.0d0)
        do j = i, 3
            rmCov(i,j) = sum((rmData(:,i)-rvMean(i))*(rmData(:,j)-rvMean(j)))/(1024.0d0 * 1024.0d0)
            if(i /= j) rmCov(j,i) = rmCov(i,j)
        end do
    end do
    open(10,file="mean_cov.csv", status="unknown", action="write")
    write(10, "('M, C1, C2, C3')")
    do i = 1, 3
        write(10, "(f9.6, 3(',',f9.6))") rvMean(i), rmCov(i,1), rmCov(i,2), rmCov(i,3)
    end do
    close(10)
    iRetCode = tRand % get(rvMean, rmCov, rmData)
    print *, 'Return code = ', iRetCode, ' (Expected = 0)'
    open(10,file="data.csv", status="unknown", action="write")
    write(10, "('x, y, z')")
    do j = 1, 1024*1024
        write(10, "(f9.6, 2(',',f9.6))") rmData(j,1), rmData(j,2), rmData(j,3)
    end do
    close(10)
    deallocate(rvMean, rmCov, rmData)

end program t_pmnrg
