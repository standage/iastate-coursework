! Daniel S. Standage
! Iowa State University
! March 10, 2011
! HW 5, CS 525
use mpi
implicit none

integer :: i,ierror
integer,parameter :: ntrial = 100
integer, parameter :: dp = mpi_double_precision, comm=mpi_comm_world
integer :: p,rank,status(mpi_status_size)
double precision :: t0, t1, t, max_time

call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Time mpi_barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call mpi_barrier(comm,ierror)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t *1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call mpi_barrier is         ', max_time,' microseconds'
endif

! 4. time the alltoall barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call barrier4(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call alltoall barrier is    ', max_time,' microseconds'
end if

! 5. time the isend/irecv barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call barrier5(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call isend/irecv barrier is ', max_time,' microseconds'
end if

! 6. time the central manager barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call barrier6(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call centman barrier is     ', max_time,' microseconds'
end if

call mpi_finalize(ierror)
end

subroutine barrier4(p, rank, comm) ! using mpi alltoall.
use mpi
implicit none
integer :: comm, p, rank, ierror
integer :: messages(p)

call mpi_alltoall( rank, 1, mpi_integer, messages, 1, mpi_integer, comm, ierror)

return
end subroutine barrier4

subroutine barrier5(p, rank, comm) ! using mpi isends and irecvs.
use mpi
implicit none
integer :: comm, p, rank, ierror, message, i, ind, sendreq, recvreq
integer :: array_of_requests(p-1), array_of_statuses(mpi_status_size, p-1)

message = 1
array_of_requests = 0


ind = 1
do i = 0, p-1
  if( i /= rank ) then
    call mpi_isend( message, 1, mpi_integer, i, 0, comm, sendreq, ierror )
    array_of_requests(ind) = sendreq
    ind = ind+1
  endif
enddo
do i = 1, p-1
  call mpi_irecv( message, 1, mpi_integer, mpi_any_source, 0, comm, recvreq, ierror )
enddo
call mpi_waitall( p-1, array_of_requests, array_of_statuses, ierror )

return
end subroutine barrier5

subroutine barrier6(p, rank, comm) ! using central manager.
use mpi
implicit none
integer :: comm, p, rank, ierror
integer :: i, message

if( rank == 0 ) then
  do i=1, p-1
    call mpi_recv( message, 1, mpi_integer, mpi_any_source, 0, comm, mpi_status_ignore, ierror )
  enddo
  do i=1, p-1
    call mpi_send( message, 1, mpi_integer, i, 0, comm, ierror )
  enddo
else
  call mpi_send( message, 1, mpi_integer, 0, 0, comm, ierror )
  call mpi_recv( message, 1, mpi_integer, 0, 0, comm, mpi_status_ignore, ierror )
endif

return
end subroutine barrier6

! ------ np=4 -----
! Time to call mpi_barrier is               24.1708755493164063 microseconds
! Time to call alltoall barrier is          39.8707389831542969 microseconds
! Time to call isend/irecv barrier is        6.4516067504882812 microseconds
! Time to call centman barrier is           26.7100334167480469 microseconds
! 
! ------ np=8 -----
! Time to call mpi_barrier is               42.4313545227050781 microseconds
! Time to call alltoall barrier is          61.3498687744140625 microseconds
! Time to call isend/irecv barrier is       13.5493278503417969 microseconds
! Time to call centman barrier is           66.5307044982910156 microseconds
! 
! ------ np=16 -----
! Time to call mpi_barrier is               63.9915466308593679 microseconds
! Time to call alltoall barrier is          87.8810882568359375 microseconds
! Time to call isend/irecv barrier is       28.5696983337402344 microseconds
! Time to call centman barrier is          145.6189155578613281 microseconds
! 
! ------ np=32 -----
! Time to call mpi_barrier is               84.9103927612304688 microseconds
! Time to call alltoall barrier is         123.3983039855957173 microseconds
! Time to call isend/irecv barrier is       60.5392456054687500 microseconds
! Time to call centman barrier is          304.8706054687500000 microseconds
