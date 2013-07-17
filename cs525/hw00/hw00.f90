! Daniel S. Standage
! Iowa State University
! February 1, 2011
! HW 0, CS 525
use mpi
implicit none

! Parameters
integer, parameter :: comm = mpi_comm_world
integer, parameter :: dp   = mpi_double_precision
integer, parameter :: n    = 3

! Variables
double precision :: A(n,n)
integer          :: i,j
integer          :: ierror
integer          :: p, rank
integer          :: next, prev
integer          :: status(mpi_status_size)

! Initialize MPI
call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Initialize A
do i = 1,n
  do j = 1,n
    A(i,j) = i + j + rank
  enddo
enddo

! Bounds checking for next and previous processes
if( rank == 0 ) then
  prev = mpi_proc_null
else
  prev = rank - 1
endif
if( rank + 1 == p ) then
  next = mpi_proc_null
else
  next = rank + 1
endif

! Send and recieve data
if( next /= mpi_proc_null ) then
  call mpi_ssend( A(1,1), n*n, dp, next, rank, comm, ierror )
endif
if( prev /= mpi_proc_null ) then
  call mpi_recv( A(1,1), n*n, dp, prev, prev, comm, status, ierror )
endif

! Print first element of the new matrix
write(*, '(A, I0, A, F5.1)') "! A(1,1) on process ", rank, ":", A(1,1)

call mpi_finalize(ierror)
end

!==============================
! Program output
!==============================
!
! [Run test p=4] mpirun -np 4 hw00
! A(1,1) on process 3:  4.0
! A(1,1) on process 1:  2.0
! A(1,1) on process 0:  2.0
! A(1,1) on process 2:  3.0

! [Run test p=8] mpirun -np 8 hw00
! A(1,1) on process 7:  8.0
! A(1,1) on process 6:  7.0
! A(1,1) on process 5:  6.0
! A(1,1) on process 4:  5.0
! A(1,1) on process 3:  4.0
! A(1,1) on process 2:  3.0
! A(1,1) on process 1:  2.0
! A(1,1) on process 0:  2.0
