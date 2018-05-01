program mpi_diffusion
use cube_mem
implicit none
real(kind=8) :: h,hsqrd, lRoom, tStep, dcoef,tAccum,cMax,cMin,urms,tStart,tEnd,&
                timeSpent
integer :: rank, ierr, i, j, k, t0=0, t1=1,it, numThreads,tSec, tMin, tHr
integer, parameter:: N=20
logical :: partition
character :: selection
logical :: isPart = .false.


    allocate(cube(0:N+1,0:N+1,0:N+1), stat=ierr)
    allocate(mask(0:N+1,0:N+1,0:N+1), stat=ierr)
    allocate(cubeCopy(N,N,N), stat=ierr)
   
    write( *, "(A40)", advance = "no")"Would you like a partition? (y/n): "
    read(*,"(A)")selection
    if(selection .eq. 'y') then
        isPart = .true.
    endif
    print *, "How many threads? "
    read *,numThreads
    
    if(isPart) then 
        do i = 1, N
            do j = 1, N
                do k =1 ,N
                        mask(j,i,k) = 1
                enddo
            enddo
        enddo

            do j = N, (N/2), -1
                do k = 1, N
                    mask(j,(N/2),k) = 0 
                enddo
            enddo
     else
        do i = 1, N
            do j = 1, N
                do k =1 ,N
                        mask(j,i,k) = 1
                enddo
            enddo
        enddo
     endif

!        do i = 0, N+1
!            do j = 0, N+1
!                do k = 0 ,N+1
!                    write(*,'(i2)',advance="NO") mask(j,i,k)
!                    if(mod(k,N+2) .eq. N+1) then
!                        print *,""
!                    endif
!                enddo
!            enddo
!                print *,""
!        enddo

    call cpu_time(tStart)
    lRoom = 5.d0
    urms = 250.d0
    tStep = (lRoom/urms)/dble(N)
    h = lRoom/dble(N)
    hsqrd = h**2.d0

    cube(1,1,1) = 1e21
    dcoef = (0.175d0 * tstep) / hsqrd
    cMax = cube(1,1,1)
    cMin = 0
    tAccum = 0

    call omp_set_num_threads(numThreads)
    do while(cMin/cMax <= 0.99)
        tAccum = tAccum + tStep
        !$OMP PARALLEL DO private(i,j,k) 
        do i = 1, N
            do j = 1, N
                do k =1 ,N
                    cubeCopy(j,i,k) = cube(j,i,k) + (&
                    cube(j+1,i,k)*mask(j+1,i,k) + cube(j-1,i,k)*mask(j-1,i,k) + &
                    cube(j,i+1,k)*mask(j,i+1,k) + cube(j,i-1,k)*mask(j,i-1,k) + &
                    cube(j,i,k+1)*mask(j,i,k+1) + cube(j,i,k-1)*mask(j,i,k-1) - &
                    dble(mask(j+1,i,k)+mask(j-1,i,k)+mask(j,i-1,k)+mask(j,i+1,k) +&
                        mask(j,i,k-1)+mask(j,i,k+1))*(cube(j,i,k)))* 2 * dcoef 
                enddo
            enddo
        enddo
            cMax = maxval(cubeCopy)
            cMin = minval(cubeCopy)

        !$OMP PARALLEL DO private(i,j,k) 
        do i = 1, N
            do j = 1, N
                do k = 1 ,N
                    cube(j,i,k) = cubeCopy(j,i,k)
                enddo
            enddo
        enddo
    enddo
    call cpu_time(tEnd)
    timeSpent = tEnd-tStart
    tSec = mod(int(timeSpent), 60)
    tMin = timeSpent / 60
    tHr = timeSpent / 3600 

    print *, ""
    write( *, "(A17)", advance = "no")"Molecular Sum: "
    print "(es20.1)",sum(cube)
    print "(a4, a20, a20, a20, a20, a20, a20)", "Dim", "Max Init", "MinInit", "Max End","Min End", "T Step",&
          &"Sim Time"
    print "(i3.0, es20.1, es20.1, es20.1, es20.1, es20.1,f20.2,i5,a1,i2,a1,i2)", N, 1E21, 0.0, cMax,&
          & cMin, tStep, tAccum  

    deallocate(cube)
    deallocate(cubeCopy) 
    deallocate(mask) 
end program mpi_diffusion
