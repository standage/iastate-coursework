! Daniel S. Standage
! Iowa State University
! March 31, 2011
! HW 6.2, CS 525
use mpi
implicit none

! Parameters
integer, parameter :: comm = mpi_comm_world
integer, parameter :: dp   = mpi_double_precision
integer, parameter :: n    = 3500

! Variables
double precision :: A(n,n), temp(n,n)
double precision :: t1, t2, time, max_time
integer          :: i,j
integer          :: ierror
integer          :: p, rank
integer          :: right, left
integer          :: sendreq, recvreq
integer          :: status(mpi_status_size)

! Initialize MPI
call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Print dimensions of A
if( rank == 0 ) then
  write(*, '(A, I4, A, I2, A)', advance="no") "! [n=", n, ", p=", p, "] "
endif

! Initialize A
do j = 1,n
  do i = 1,n
    A(i,j) = i + j + rank
  enddo
enddo

! Bounds checking for next and previous processes
if( rank == 0 ) then
  left = mpi_proc_null
else
  left = rank - 1
endif
if( rank + 1 == p ) then
  right = mpi_proc_null
else
  right = rank + 1
endif

! Send and receive--version 1
t1 = mpi_wtime()
if( right /= mpi_proc_null ) then
  call mpi_ssend( A(1,1), n*n, dp, right, 0, comm, ierror )
endif
if( left /= mpi_proc_null ) then
  call mpi_recv( A(1,1), n*n, dp, left, 0, comm, mpi_status_ignore, ierror )
endif
t2 = mpi_wtime()
time = t2 - t1
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
if( rank == 0 ) then
  write(*, '(A, ES12.6)', advance="no") "   Test 1) ", max_time
endif

! Send and receive--version 2
t1 = mpi_wtime()
if( right /= mpi_proc_null ) then
  call mpi_send( A(1,1), n*n, dp, right, 0, comm, ierror )
endif
if( left /= mpi_proc_null ) then
  call mpi_recv( A(1,1), n*n, dp, left, 0, comm, mpi_status_ignore, ierror )
endif
t2 = mpi_wtime()
time = t2 - t1
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
if( rank == 0 ) then
  write(*, '(A, ES12.6)', advance="no") "    Test 2) ", max_time
endif

! Send and receive--version 3
t1 = mpi_wtime()
  call mpi_sendrecv_replace( A(1,1), n*n, dp, right, 0, left, 0, comm, mpi_status_ignore, ierror )
t2 = mpi_wtime()
time = t2 - t1
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
if( rank == 0 ) then
  write(*, '(A, ES12.6)', advance="no") "    Test 3) ", max_time
endif

! Send and receive--version 4
t1 = mpi_wtime()
if( left /= mpi_proc_null ) then
  call mpi_irecv( temp(1,1), n*n, dp, left, 0, comm, recvreq, ierror )
endif
if( right /= mpi_proc_null ) then
  call mpi_send( A(1,1), n*n, dp, right, 0, comm, ierror )
endif
if( left /= mpi_proc_null ) then
  call mpi_wait(recvreq, mpi_status_ignore, ierror)
  do i = 1,n
    do j = 1,n
      A(j,i) = temp(j,i)
    enddo
  enddo
endif
t2 = mpi_wtime()
time = t2 - t1
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
if( rank == 0 ) then
  write(*, '(A, ES12.6)', advance="no") "    Test 4) ", max_time
endif

if( rank == 0 ) then
  write(*, '(A)') ""
endif

call mpi_finalize(ierror)
end

! The new implementation (using mpi_irecv) performs similar to versions 2 and 3
! from the previous homework assignment, except for n=250, for which versions 2
! and 3 provide a consistent performance improvement (the new version is almost
! two times slower than version 2 for n=250).
!
!-------------------
! n=10
!-------------------
! [n=  10, p= 4]    Test 1) 8.538961E-03    Test 2) 3.190041E-04    Test 3) 3.149509E-04    Test 4) 2.160072E-04
! [n=  10, p= 8]    Test 1) 2.254009E-02    Test 2) 1.947880E-04    Test 3) 1.900196E-04    Test 4) 1.940727E-04
! [n=  10, p=16]    Test 1) 5.629206E-02    Test 2) 2.131462E-04    Test 3) 2.059937E-04    Test 4) 2.119541E-04
! [n=  10, p=32]    Test 1) 8.984995E-02    Test 2) 1.940727E-04    Test 3) 1.888275E-04    Test 4) 1.950264E-04
!-------------------
! n=250
!-------------------
! [n= 250, p= 4]    Test 1) 1.219106E-02    Test 2) 3.722906E-03    Test 3) 5.004883E-03    Test 4) 6.593227E-03
! [n= 250, p= 8]    Test 1) 3.136110E-02    Test 2) 3.692865E-03    Test 3) 5.844831E-03    Test 4) 6.947041E-03
! [n= 250, p=16]    Test 1) 7.064390E-02    Test 2) 3.722906E-03    Test 3) 5.760908E-03    Test 4) 7.803202E-03
! [n= 250, p=32]    Test 1) 1.583230E-01    Test 2) 3.757000E-03    Test 3) 5.759954E-03    Test 4) 7.580996E-03
!-------------------
! n=3500
!-------------------
! [n=3500, p= 4]    Test 1) 9.494319E-01    Test 2) 6.679900E-01    Test 3) 8.985791E-01    Test 4) 9.612288E-01
! [n=3500, p= 8]    Test 1) 2.390953E+00    Test 2) 9.196651E-01    Test 3) 1.059201E+00    Test 4) 1.074881E+00
! [n=3500, p=16]    Test 1) 5.017354E+00    Test 2) 1.288790E+00    Test 3) 1.059065E+00    Test 4) 1.070065E+00
! [n=3500, p=32]    Test 1) 1.021079E+01    Test 2) 1.807550E+00    Test 3) 1.062239E+00    Test 4) 1.073590E+00
