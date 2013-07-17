! Daniel S. Standage
! Iowa State University
! February 17, 2011
! HW 2, CS 525
use mpi
implicit none

! Parameters
integer, parameter :: comm = mpi_comm_world
integer, parameter :: dp   = mpi_double_precision
integer, parameter :: n    = 3500

! Variables
double precision :: A(n,n)
double precision :: t1, t2, time, max_time
integer          :: i,j
integer          :: ierror
integer          :: p, rank
integer          :: right, left
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

if( rank == 0 ) then
  write(*, '(A)') ""
endif
! Print first element of the new matrix
! write(*, '(A, I0, A, F5.1)') "! A(1,1) on process ", rank, ":", A(1,1)

call mpi_finalize(ierror)
end

! [[[ Output ]]]
!-------------------
! n=10
!-------------------
! [n=  10, p= 4]    Test 1) 1.036310E-02    Test 2) 3.130436E-04    Test 3) 3.080368E-04
! [n=  10, p= 4]    Test 1) 8.574963E-03    Test 2) 1.819134E-04    Test 3) 5.602837E-05
! [n=  10, p= 4]    Test 1) 8.589029E-03    Test 2) 1.780987E-04    Test 3) 5.793571E-05
! [n=  10, p= 8]    Test 1) 7.986307E-02    Test 2) 2.048016E-04    Test 3) 2.009869E-04
! [n=  10, p= 8]    Test 1) 3.316689E-02    Test 2) 2.613282E-02    Test 3) 1.423836E-03
! [n=  10, p= 8]    Test 1) 3.376007E-02    Test 2) 2.689886E-02    Test 3) 1.969337E-04
! [n=  10, p=16]    Test 1) 5.573297E-02    Test 2) 2.009869E-04    Test 3) 1.981258E-04
! [n=  10, p=16]    Test 1) 5.740499E-02    Test 2) 2.040863E-04    Test 3) 2.000332E-04
! [n=  10, p=16]    Test 1) 4.506993E-02    Test 2) 7.762194E-03    Test 3) 7.755041E-03
! [n=  10, p=32]    Test 1) 2.580500E-01    Test 2) 3.090143E-03    Test 3) 1.959801E-04
! [n=  10, p=32]    Test 1) 1.131599E-01    Test 2) 2.090931E-04    Test 3) 2.059937E-04
! [n=  10, p=32]    Test 1) 8.968115E-02    Test 2) 2.758503E-04    Test 3) 1.990795E-04
!-------------------
! n=250
!-------------------
! [n= 250, p= 4]    Test 1) 1.216006E-02    Test 2) 3.748178E-03    Test 3) 5.181789E-03
! [n= 250, p= 4]    Test 1) 1.224399E-02    Test 2) 3.676891E-03    Test 3) 5.119085E-03
! [n= 250, p= 4]    Test 1) 1.232386E-02    Test 2) 3.633976E-03    Test 3) 5.120039E-03
! [n= 250, p= 8]    Test 1) 2.426195E-02    Test 2) 3.745079E-03    Test 3) 5.535126E-03
! [n= 250, p= 8]    Test 1) 2.600789E-02    Test 2) 3.737926E-03    Test 3) 5.486012E-03
! [n= 250, p= 8]    Test 1) 2.587581E-02    Test 2) 3.741980E-03    Test 3) 5.547047E-03
! [n= 250, p=16]    Test 1) 6.209302E-02    Test 2) 3.748894E-03    Test 3) 5.475998E-03
! [n= 250, p=16]    Test 1) 6.172895E-02    Test 2) 3.770828E-03    Test 3) 5.526066E-03
! [n= 250, p=16]    Test 1) 6.205893E-02    Test 2) 3.782034E-03    Test 3) 5.621910E-03
! [n= 250, p=32]    Test 1) 1.458151E-01    Test 2) 3.751040E-03    Test 3) 5.538940E-03
! [n= 250, p=32]    Test 1) 1.425140E-01    Test 2) 3.773928E-03    Test 3) 5.445957E-03
! [n= 250, p=32]    Test 1) 3.750579E-01    Test 2) 2.330971E-01    Test 3) 2.345929E-01
!-------------------
! n=3500
!-------------------
! [n=3500, p= 4]    Test 1) 7.957239E-01    Test 2) 6.021209E-01    Test 3) 7.569580E-01
! [n=3500, p= 4]    Test 1) 7.952299E-01    Test 2) 5.985160E-01    Test 3) 7.606289E-01
! [n=3500, p= 4]    Test 1) 7.836478E-01    Test 2) 5.932260E-01    Test 3) 7.528610E-01
! [n=3500, p= 8]    Test 1) 2.396737E+00    Test 2) 9.270079E-01    Test 3) 1.058917E+00
! [n=3500, p= 8]    Test 1) 2.314695E+00    Test 2) 9.275730E-01    Test 3) 1.050240E+00
! [n=3500, p= 8]    Test 1) 2.299999E+00    Test 2) 9.280930E-01    Test 3) 1.056798E+00
! [n=3500, p=16]    Test 1) 5.087708E+00    Test 2) 1.213626E+00    Test 3) 1.049239E+00
! [n=3500, p=16]    Test 1) 5.079564E+00    Test 2) 1.209199E+00    Test 3) 1.048371E+00
! [n=3500, p=16]    Test 1) 5.103053E+00    Test 2) 1.205890E+00    Test 3) 1.044103E+00
! [n=3500, p=32]    Test 1) 1.046298E+01    Test 2) 1.781815E+00    Test 3) 1.046560E+00
! [n=3500, p=32]    Test 1) 1.042169E+01    Test 2) 1.781001E+00    Test 3) 1.048451E+00
! [n=3500, p=32]    Test 1) 1.045141E+01    Test 2) 1.766217E+00    Test 3) 1.049237E+00
!
! [[[ Questions ]]]
!
! 1. The right shift is executed serially with 'mpi_ssend' and for large messages with
!    'mpi_send'. The right shift is executed in parallel for small messages with 'mpi_send'
!    and for all messages with 'mpi_sendrecv_replace'.
! 2. Consequently, 'mpi_sendrecv_replace' gives the best overall performance.
