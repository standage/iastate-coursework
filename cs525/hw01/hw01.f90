! Daniel S. Standage
! Iowa State University
! February 10, 2011
! HW 1, CS 525
use mpi
implicit none

! Variables and paraneters
integer, parameter :: dp=mpi_double_precision, comm = mpi_comm_world
integer, parameter :: n = 64
integer            :: i, ierror, p, rank, status(mpi_status_size)
double precision   :: A(n,n), B(n,n), C(n,n)
double precision   :: t1, t2, time, max_time, min_time, tmp
double precision, allocatable :: time_array(:)

! Initialize MPI
call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Allocate time_array and initialize A, B and C:
allocate(time_array(0:p-1))
time_array(0:p-1) = -1.d0
call random_number(A)
call random_number(B)
A = A + float(rank)
B = B - float (rank)
C = 0.d0

! Do the matrix multiplication
call mpi_barrier(comm, ierror)
t1 = mpi_wtime()
C = C + matmul(A,B)
t2 = mpi_wtime()
time = t2 - t1

! Find the maximum and minimum of all the times and put on process 0
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
call mpi_reduce(time, min_time, 1, dp, mpi_min, 0, comm, ierror)

! Send time to process 0 and receive into time_array
if(rank == 0) then
  time_array(0) = time
  do i = 1, p-1
    call mpi_recv( tmp, 1, dp, mpi_any_source, 0, comm, status, ierror )
    time_array( status(mpi_source) ) = tmp
  enddo
else
  call mpi_ssend( time, 1, dp, 0, 0, comm, ierror )
endif

! Print results
if(rank == 0) then
  do i = 0, p-1
    print*,'For rank = ',i,' time = ',time_array(i)
  enddo
  print*,' '
  print*,' maximum time = ', max_time,' seconds'
  print*,' minimum time = ', min_time,' seconds'
  print*,' '
  print*,'The value of n used was n = ', n
endif

call mpi_finalize(ierror)
end

!---------------------------------------------------------
! 4 processors
!---------------------------------------------------------
! For rank =            0  time =   2.064943313598633E-003
! For rank =            1  time =   2.057075500488281E-003
! For rank =            2  time =   2.061128616333008E-003
! For rank =            3  time =   2.047061920166016E-003
!  
!  maximum time =   2.064943313598633E-003  seconds
!  minimum time =   2.047061920166016E-003  seconds

!---------------------------------------------------------
! 8 processors
!---------------------------------------------------------
! The value of n used was n =           64
! For rank =            0  time =   2.053976058959961E-003
! For rank =            1  time =   2.053022384643555E-003
! For rank =            2  time =   2.046108245849609E-003
! For rank =            3  time =   2.052068710327148E-003
! For rank =            4  time =   2.058982849121094E-003
! For rank =            5  time =   2.063035964965820E-003
! For rank =            6  time =   2.072095870971680E-003
! For rank =            7  time =   2.053022384643555E-003
!  
!  maximum time =   2.072095870971680E-003  seconds
!  minimum time =   2.046108245849609E-003  seconds

!---------------------------------------------------------
! 16 processors
!---------------------------------------------------------
! The value of n used was n =           64
! For rank =            0  time =   2.054929733276367E-003
! For rank =            1  time =   2.058029174804688E-003
! For rank =            2  time =   2.053022384643555E-003
! For rank =            3  time =   2.063989639282227E-003
! For rank =            4  time =   2.087116241455078E-003
! For rank =            5  time =   2.057075500488281E-003
! For rank =            6  time =   2.068042755126953E-003
! For rank =            7  time =   2.055168151855469E-003
! For rank =            8  time =   2.052068710327148E-003
! For rank =            9  time =   2.058982849121094E-003
! For rank =           10  time =   2.053022384643555E-003
! For rank =           11  time =   2.048969268798828E-003
! For rank =           12  time =   2.074003219604492E-003
! For rank =           13  time =   2.054929733276367E-003
! For rank =           14  time =   2.054929733276367E-003
! For rank =           15  time =   2.061843872070312E-003
!  
!  maximum time =   2.087116241455078E-003  seconds
!  minimum time =   2.048969268798828E-003  seconds
!  
! The value of n used was n =           64
