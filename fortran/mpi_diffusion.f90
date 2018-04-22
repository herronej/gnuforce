program mpi_diffusion
USE MPI
implicit none
integer :: ierr

    call MPI_INIT(ierr)

    call MPI_FINALIZE(ierr)
    
end program mpi_diffusion
