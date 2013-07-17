! Daniel S. Standage
! Iowa State University
! February 24, 2011
! HW 3, CS 525
use mpi
implicit none

integer, parameter :: dp=mpi_double_precision, comm = mpi_comm_world
integer :: i, j, ierror, p, rank, status(mpi_status_size)
integer, parameter :: m = 1 , n = 1 , nloops = 5000
!!integer, parameter :: m = 2500, n = 4000, nloops = 5
double precision :: A(m,n), B(m,n), t1, t2, time, bandwidth

call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

!--------------------------------------------------
! Part 1: ping-pong using 'mpi_ssend'
!--------------------------------------------------
if(rank==0) then
  write(*, '(A, I0, A, I0, A, I0)') '! m=',m,' n=',n,' nloops=', nloops
  write(*,'(A)') '!     Part 1 using "mpi_ssend"'
endif

! Initialize B and A
B(1:m,1:n) = 0.d0
if (rank==0) then
  call random_number(A)
end if

! Time send and receives, compute bandwidth
do i = 1, p-1
  call mpi_barrier(comm, ierror)
  ! Proc 0 sends and receives every time
  if( rank == 0 ) then
    t1 = mpi_wtime()
    do j = 1, nloops
      call mpi_ssend( A(1,1), m*n, dp, i, 0, comm, ierror )
      call mpi_recv( B(1,1), m*n, dp, i, 0, comm, mpi_status_ignore, ierror )
    enddo
    t2 = mpi_wtime()
    time = (t2 - t1)/nloops
    bandwidth = dble(8*m*n)/(time*(1024**2))
    write(*, '(A, I0)') '!         Processor i=', i
    write(*, '(A, ES12.6, A)') '!             time:      ', time, ' seconds'
    write(*, '(A, ES12.6, A)') '!             bandwidth: ', bandwidth, ' MBps'
  else
    ! Other procs only receive and send if it's their turn
    if( rank == i ) then
      do j = 1, nloops
        call mpi_recv( A(1,1), m*n, dp, 0, 0, comm, mpi_status_ignore, ierror )
        call mpi_ssend( B(1,1), m*n, dp, 0, 0, comm, ierror )
      enddo
    endif
  endif
enddo



!--------------------------------------------------
! Part 2: ping-pong using 'mpi_send'
!--------------------------------------------------
if(rank==0) then
  write(*,'(A)') '!     Part 2 using "mpi_send"'
endif

! Initialize B and A
B(1:m,1:n) = 0.d0
if (rank==0) then
  call random_number(A)
end if

! Time send and receives, compute bandwidth
do i = 1, p-1
  call mpi_barrier(comm, ierror)
  ! Proc 0 sends and receives every time
  if( rank == 0 ) then
    t1 = mpi_wtime()
    do j = 1, nloops
      call mpi_send( A(1,1), m*n, dp, i, 0, comm, ierror )
      call mpi_recv( B(1,1), m*n, dp, i, 0, comm, mpi_status_ignore, ierror )
    enddo
    t2 = mpi_wtime()
    time = (t2 - t1)/nloops
    bandwidth = dble(8*m*n)/(time*(1024**2))
    write(*, '(A, I0)') '!         Processor i=', i
    write(*, '(A, ES12.6, A)') '!             time:      ', time, ' seconds'
    write(*, '(A, ES12.6, A)') '!             bandwidth: ', bandwidth, ' MBps'
  else
    ! Other procs only receive and send if it's their turn
    if( rank == i ) then
      do j = 1, nloops
        call mpi_recv( A(1,1), m*n, dp, 0, 0, comm, mpi_status_ignore, ierror )
        call mpi_send( B(1,1), m*n, dp, 0, 0, comm, ierror )
      enddo
    endif
  endif
enddo



!--------------------------------------------------
! Part 3: ping-pong using 'mpi_sendrecv'
!--------------------------------------------------
if(rank==0) then
  write(*,'(A)') '!     Part 3 using "mpi_sendrecv"'
endif

! Initialize B and A
B(1:m,1:n) = 0.d0
if (rank==0) then
  call random_number(A)
end if

! Time send and receives, compute bandwidth
do i = 1, p-1
  call mpi_barrier(comm, ierror)
  ! Proc 0 sends and receives every time
  if( rank == 0 ) then
    t1 = mpi_wtime()
    do j = 1, nloops
      call mpi_sendrecv( A(1,1), m*n, dp, i, 0, B(1,1), m*n, dp, i, 0, comm, mpi_status_ignore, ierror )
    enddo
    t2 = mpi_wtime()
    time = (t2 - t1)/nloops
    bandwidth = dble(8*m*n)/(time*(1024**2))
    write(*, '(A, I0)') '!         Processor i=', i
    write(*, '(A, ES12.6, A)') '!             time:      ', time, ' seconds'
    write(*, '(A, ES12.6, A)') '!             bandwidth: ', bandwidth, ' MBps'
  else
    ! Other procs only receive and send if it's their turn
    if( rank == i ) then
      do j = 1, nloops
        call mpi_sendrecv( B(1,1), m*n, dp, 0, 0, A(1,1), m*n, dp, 0, 0, comm, mpi_status_ignore, ierror )
      enddo
    endif
  endif
enddo


call mpi_finalize(ierror)
end

! [[[ Output ]]]
!
! m=1 n=1 nloops=5000
!     Part 1 using "mpi_ssend"
!         Processor i=1
!             time:      5.817413E-06 seconds
!             bandwidth: 1.311475E+00 MBps
!         Processor i=2
!             time:      5.662503E-05 seconds
!             bandwidth: 1.347354E-01 MBps
!         Processor i=3
!             time:      5.661821E-05 seconds
!             bandwidth: 1.347516E-01 MBps
!         Processor i=4
!             time:      5.567183E-05 seconds
!             bandwidth: 1.370423E-01 MBps
!         Processor i=5
!             time:      5.567141E-05 seconds
!             bandwidth: 1.370433E-01 MBps
!         Processor i=6
!             time:      5.660219E-05 seconds
!             bandwidth: 1.347897E-01 MBps
!         Processor i=7
!             time:      5.672083E-05 seconds
!             bandwidth: 1.345078E-01 MBps
!     Part 2 using "mpi_send"
!         Processor i=1
!             time:      2.166224E-06 seconds
!             bandwidth: 3.521979E+00 MBps
!         Processor i=2
!             time:      1.561942E-05 seconds
!             bandwidth: 4.884557E-01 MBps
!         Processor i=3
!             time:      1.560221E-05 seconds
!             bandwidth: 4.889946E-01 MBps
!         Processor i=4
!             time:      1.526418E-05 seconds
!             bandwidth: 4.998235E-01 MBps
!         Processor i=5
!             time:      1.529679E-05 seconds
!             bandwidth: 4.987578E-01 MBps
!         Processor i=6
!             time:      1.560040E-05 seconds
!             bandwidth: 4.890514E-01 MBps
!         Processor i=7
!             time:      1.561298E-05 seconds
!             bandwidth: 4.886570E-01 MBps
!     Part 3 using "mpi_sendrecv"
!         Processor i=1
!             time:      2.114391E-06 seconds
!             bandwidth: 3.608317E+00 MBps
!         Processor i=2
!             time:      1.006818E-05 seconds
!             bandwidth: 7.577731E-01 MBps
!         Processor i=3
!             time:      1.013703E-05 seconds
!             bandwidth: 7.526260E-01 MBps
!         Processor i=4
!             time:      9.900999E-06 seconds
!             bandwidth: 7.705681E-01 MBps
!         Processor i=5
!             time:      9.926176E-06 seconds
!             bandwidth: 7.686137E-01 MBps
!         Processor i=6
!             time:      1.004581E-05 seconds
!             bandwidth: 7.594600E-01 MBps
!         Processor i=7
!             time:      1.007180E-05 seconds
!             bandwidth: 7.575004E-01 MBps
! m=2500 n=4000 nloops=5
!     Part 1 using "mpi_ssend"
!         Processor i=1
!             time:      4.187000E-01 seconds
!             bandwidth: 1.822163E+02 MBps
!         Processor i=2
!             time:      6.497050E-01 seconds
!             bandwidth: 1.174286E+02 MBps
!         Processor i=3
!             time:      6.485072E-01 seconds
!             bandwidth: 1.176455E+02 MBps
!         Processor i=4
!             time:      6.484068E-01 seconds
!             bandwidth: 1.176637E+02 MBps
!         Processor i=5
!             time:      6.482488E-01 seconds
!             bandwidth: 1.176924E+02 MBps
!         Processor i=6
!             time:      6.484690E-01 seconds
!             bandwidth: 1.176524E+02 MBps
!         Processor i=7
!             time:      6.488180E-01 seconds
!             bandwidth: 1.175891E+02 MBps
!     Part 2 using "mpi_send"
!         Processor i=1
!             time:      4.184922E-01 seconds
!             bandwidth: 1.823067E+02 MBps
!         Processor i=2
!             time:      6.480902E-01 seconds
!             bandwidth: 1.177212E+02 MBps
!         Processor i=3
!             time:      6.480648E-01 seconds
!             bandwidth: 1.177258E+02 MBps
!         Processor i=4
!             time:      6.478212E-01 seconds
!             bandwidth: 1.177701E+02 MBps
!         Processor i=5
!             time:      6.478466E-01 seconds
!             bandwidth: 1.177654E+02 MBps
!         Processor i=6
!             time:      6.480748E-01 seconds
!             bandwidth: 1.177240E+02 MBps
!         Processor i=7
!             time:      6.480920E-01 seconds
!             bandwidth: 1.177209E+02 MBps
!     Part 3 using "mpi_sendrecv"
!         Processor i=1
!             time:      4.162866E-01 seconds
!             bandwidth: 1.832726E+02 MBps
!         Processor i=2
!             time:      4.402070E-01 seconds
!             bandwidth: 1.733138E+02 MBps
!         Processor i=3
!             time:      4.420750E-01 seconds
!             bandwidth: 1.725814E+02 MBps
!         Processor i=4
!             time:      4.415510E-01 seconds
!             bandwidth: 1.727863E+02 MBps
!         Processor i=5
!             time:      4.422062E-01 seconds
!             bandwidth: 1.725303E+02 MBps
!         Processor i=6
!             time:      4.416720E-01 seconds
!             bandwidth: 1.727389E+02 MBps
!         Processor i=7
!             time:      4.420462E-01 seconds
!             bandwidth: 1.725927E+02 MBps
!
! [[[ Discussion ]]]
!
! I ran the ping-pong tests with both sets of parameters (m=1, n=1, nloops=5000
! versus m=2500, n=4000, nloops=5). Processor 1 had the best times in all
! cases, so I assume it was on the same node as processor 0. It was interesting
! to note that for processor 1, the send method used did not seem to have a
! significant effect on time or bandwidth.
!
! For all the other processors, each test offered an improvement in time and
! bandwidth. This improvement was pretty very pronounced in the latency test
! with small messages, whereas it was only a slight improvement when sending
! large messages. So while 'mpi_sendrecv' is clearly the best choice in terms of
! send time and bandwidth, the margin of improvement shrinks as the message size
! grows.
