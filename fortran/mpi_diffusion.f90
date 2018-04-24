program mpi_diffusion
USE MPI
implicit none
integer :: rank, size, ierr, N=0,i
logical :: partition

    call MPI_INIT(ierr)

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size, ierr)
    if(rank .eq. 0) then
        write(*,'(A)', ADVANCE = "NO") 'Please enter value number of dimensions: '
            read(*,*) N
        call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        print *, N
    else
        call MPI_RECV(N, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            print *, "worker " , N
    endif

    call MPI_FINALIZE(ierr)
    
end program mpi_diffusion
