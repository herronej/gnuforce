! CSC 330
! Assignment 2 - Diffusion
!
! Author: Devin McBryde

program diff

        !needed to suppress unwanted behavior in fortran
        implicit none


        !Determines the number of divisions used n each dimension of the room
        integer :: N
        
        !The declaration of all necessary variables needed for the program 
        real(kind=8) :: mTotal
        real(kind=8) :: tot
        real(kind=8) :: mSpeed
        real(kind=8) :: hval
        real(kind=8) :: D
        real(kind=8) :: conMax
        real(kind=8) :: conMin
        real(kind=8) :: tStep
        real(kind=8) :: time

        real         :: start, endTime

        !A 3 dimensional array that will operate as a rank 3 tensor used 
        !       to represent the room 

	real(kind=8), dimension(:,:,:), allocatable :: room 


        !This array will store the different values of 
        !       change in concentration between two cells
        !       The name means concentration difference
        real(kind=8), dimension(6) :: dCon

        !Every time we check to see the flux of gas between cells
        !       we would also need to multiply several values,
        !       slowing the speed of computation. By calculating the
        !       value once we only need to perform a single
        !       multiplication each time afterwards for each cell
        !       instead of several
        real         :: coefficient

        logical      :: partition

        !declaring the variables that will be used for looping through the room
        integer      :: i,j,k

        print*, "Enter Number of Divisions for Room "
        read*, N

        allocate( room(N,N,N))

        call cpu_time(start)

        !initalizing the time to zero
        time = 0

        !initalizing the variables that describe the state of the room and
        !       information about the molecules of gas
        mTotal = 1000000000000000000000.0
        mSpeed = 250.0
        hval   = 5.0/N
        D      = 0.175
        conMax = 0.0
        conMin = 0.0
        tStep  = hval/mSpeed

        !logical type variable that controls whether or not the program will 
        !       be run with a partition in the room
        partition = .false.

        coefficient = (tStep*D)/(hval*hval)

        !Following for loops will initialize the room tensor with 0 values 
        !       when partioning is turned off, otherwise locations that 
        !       represent the partion in the room will be initialized
        !       to the value -1
        do i = 1, N
                do j = 1, N
                        do k = 1, N
                                !negative values are used when partition is true and will
                                !       place them half way into the room (when j == (N/2)-1)
                                !       and half way up (when i >= (N/2)-1)
                                if (j==((N/2)) .and. i>=((N/2)) .and. partition) then
                                        room(i,j,k) = -1.0
                                else
                                        room(i,j,k) = 0.0
                                end if
                        end do
                end do
        end do

        !Initializing the dCon array
        do i = 1, 6
                dCon(i) = 0
        end do

        !This is where we put the gas in the room that will diffuse
        !       This can be thought of as the top corner of the room
        room(1,1,1) = mTotal

        !initializing these values so that the loop can start, some values would
        !       prevent this
        conMax = mTotal
        conMin = 1.0

!We want the simulation to stop when the room has become sufficiently
!       diffuse with the gas, thus we check if the ratio of lowest
!       concentration to highest is less than 0.99, and when it is 
!       higher we know the gas has diffused
do while ((conMin/conMax) .lt. 0.99)
        
        !every step of the program has the time tick by the tStep variable
        !       tStep is based on qualities of the diffuse material and the room
        time = time + tStep

        do i = 1, N
                do j = 1, N
                        do k = 1, N
                               
                                !calculate the difference in concentration from flux with each cube face
                                !       The 6 faces of the cube are represented with different address values
                                !       and an if is used to determine if it is safe to move molecules
                                !       if a value = N-1  or 0 then we have hit a face of the cube and do not calculate 

                                if (room(i,j,k) .ne. -1) then
                                       
                                        if (k==N .or. room(i,j,k+1) == -1) then
                                                dCon(1) = 0
                                        else
                                                dcon(1) = (room(i,j,k)-room(i,j,k+1)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(1)
                                                room(i,j,k+1) = room(i,j,k+1) + dCon(1)
                                        end if
        
                                        if (j==N .or. room(i,j+1,k) == -1) then
                                                dCon(2) = 0
                                        else
                                                dcon(2) = (room(i,j,k)-room(i,j+1,k)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(2)
                                                room(i,j+1,k) = room(i,j+1,k) + dCon(2)
                                        end if
        
                                        if (i==N .or. room(i+1,j,k) == -1) then
                                                dCon(3) = 0 
                                        else
                                                dcon(3) = (room(i,j,k)-room(i+1,j,k)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(3)
                                                room(i+1,j,k) = room(i+1,j,k) + dCon(3)
                                        end if
        
                                        if (k==1 .or. room(i,j,k-1) == -1) then
                                                dCon(4) = 0
                                        else
                                                dcon(4) = (room(i,j,k)-room(i,j,k-1)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(4)
                                                room(i,j,k-1) = room(i,j,k-1) + dCon(4)
                                        end if
        
                                        if (j==1 .or. room(i,j-1,k) == -1) then
                                                dCon(5) = 0
                                        else
                                                dcon(5) = (room(i,j,k)-room(i,j-1,k)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(5)
                                                room(i,j-1,k) = room(i,j-1,k) + dCon(5)
                                        end if
        
                                        if (i==1 .or. room(i-1,j,k) == -1) then
                                                dCon(6) = 0
                                        else
                                                dcon(6) = (room(i,j,k)-room(i-1,j,k)) * coefficient
                                                room(i,j,k) = room(i,j,k) - dCon(6)
                                                room(i-1,j,k) = room(i-1,j,k) + dCon(6)
                                        end if
                                end if
                        end do
                end do
        end do

        !Setting up these variables with some set value from the room
        !       to prime them for the max and min search
        conMin = room(1,1,1)
        conMax = room(1,1,1)

         do i = 1, N
                do j = 1, N
                        do k = 1, N

                        if (room(i,j,k) < conMin .and. room(i,j,k) .ne. -1) then
                                conMin = room(i,j,k)
                        end if

                        if (room(i,j,k) > conMax .and. room(i,j,k) .ne. -1) then
                                conMax = room(i,j,k)
                        end if

                        end do
                end do
        end do

        
end do

        !To check for matter consistency we check the total amount of molecules in the
        !       room after the smulation and output it with the starting value
        do i = 1, N
                do j = 1, N
                        do k = 1, N
                                tot = tot + room(i,j,k)
                        end do
                end do
        end do

        call cpu_time(endTime)

        !output of the simulation detailing 5 vaules
        !       How many molecules did we start with
        !       How many molecules did we end with
        !       The total amount of time it took for the room to become diffused
        !       The minimum concentration in the room
        !       the maximum concentration in the room
        print *, "Total molecules starting:", mTotal
        print *, "Total molecules left:", tot
        print *, "Time Simulated:", time
        print *, "min concentration:", conMin
        print *, "max concentration:", conMax
        print *, "Wall time:", endTime-start


end program diff
