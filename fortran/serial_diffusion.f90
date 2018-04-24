program mpi_diffusion
use cube_mem
implicit none
real(kind=8) deltas, lRoom, tStep, dcoef
integer :: rank, size, ierr, i, j, k, t0=0, t1=1,it
integer, parameter :: N=10
logical :: partition
    
    allocate(cube(0:N+1,0:N+1,0:N+1), stat=ierr)
    allocate(cubeCopy(N,N,N), stat=ierr)
    lRoom = 5
    tStep = (lRoom/250.d0)/dble(N)
    deltas = lRoom/dble(N)
    deltas = deltas**2.d0

    cube(1,1,1) = 1e21
    dcoef = 0.175 * tstep / deltas

    do it=1, 1000
        do i = 1, N
            do j = 1, N
                do k =1 ,N
                    !Compare this to equation in book.  I also added
                    !the diffusion coefficient which you were missing.
                    cubeCopy(j,i,k) = cube(j,i,k) + (&
                    cube(j+1,i,k) + cube(j-1,i,k) + &
                    cube(j,i+1,k) + cube(j,i-1,k) + &
                    cube(j,i,k+1) + cube(j,i,k-1) - &
                    6.0D0*(cube(j,i,k))) * dcoef / deltas
                enddo
            enddo
        enddo
        do i = 1, N
            do j = 1, N
                do k = 1 ,N
                    cube(j,i,k) = cubeCopy(j,i,k)
                    print *, cube(j,i,k)
                enddo
            enddo
        enddo
        
    enddo
    print *, sum(cube) 
    deallocate(cube)
    deallocate(cubeCopy) 
end program mpi_diffusion
