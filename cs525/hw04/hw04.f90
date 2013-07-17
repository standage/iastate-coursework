! Daniel S. Standage
! Iowa State University
! March 3, 2011
! HW 4, CS 525
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
  write(*, '(A, F, A)') 'Time to call mpi_barrier is       ', max_time,' microseconds'
endif

! 1. time the send - recv barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call barrier1(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call send/recv barrier is ', max_time,' microseconds'
end if

! 2. time the bcast barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime()
do i = 1, ntrial
  call barrier2(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call bcast barrier is     ', max_time,' microseconds'
end if

! 3. time the scatter barrier
call mpi_barrier(comm,ierror)
t0 = mpi_wtime() 
do i = 1, ntrial 
  call barrier3(p, rank, comm)
enddo
t1 = mpi_wtime()
t = (t1 - t0)/dble(ntrial)
t = t*1.d6 ! convert from seconds to microseconds
call mpi_reduce(t, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank.eq.0) then
  write(*, '(A, F, A)') 'Time to call scatter barrier is   ', max_time,' microseconds'
end if

call mpi_finalize(ierror)
end

subroutine barrier1(p, rank, comm) ! using mpi sends and recvs.
use mpi
implicit none
integer :: comm, p, rank, ierror, message, i

message = 0

do i = 0, p-1
  if( i /= rank ) then
    call mpi_send( message, 1, mpi_integer, i, 0, comm, ierror )
  endif
enddo
do i = 1, p-1
  call mpi_recv( message, 1, mpi_integer, mpi_any_source, 0, comm, mpi_status_ignore, ierror )
enddo

return
end subroutine barrier1

subroutine barrier2(p, rank, comm) ! using mpi broadcast.
use mpi
implicit none
integer :: comm, p, rank, ierror, message, i

message = 1
do i=0, p-1
  call mpi_bcast( message, 1, mpi_integer, i, comm, ierror )
enddo

return
end subroutine barrier2

subroutine barrier3(p, rank, comm) ! using mpi gather/scatter.
use mpi
implicit none
integer :: comm, p, rank, ierror, a(0:p-1), b

a(0:p-1) = 0
call mpi_gather(b, 1, mpi_integer, a(0), 1, mpi_integer, 0, comm, ierror)
call mpi_scatter(a(0), 1, mpi_integer, b, 1, mpi_integer, 0, comm, ierror)

return
end subroutine barrier3

! ------ np=4 -----
! Time to call mpi_barrier is             24.4021415710449219 microseconds
! Time to call send/recv barrier is       43.2395935058593750 microseconds
! Time to call bcast barrier is           48.7303733825683594 microseconds
! Time to call scatter barrier is         20.3108787536621094 microseconds
! 
! ------ np=8 -----
! Time to call mpi_barrier is             44.0406799316406250 microseconds
! Time to call send/recv barrier is      124.0301132202148437 microseconds
! Time to call bcast barrier is          152.0800590515136719 microseconds
! Time to call scatter barrier is         36.7188453674316406 microseconds
! 
! ------ np=16 -----
! Time to call mpi_barrier is             63.6005401611328125 microseconds
! Time to call send/recv barrier is      266.9501304626464844 microseconds
! Time to call bcast barrier is          396.5902328491210937 microseconds
! Time to call scatter barrier is         52.8001785278320313 microseconds
! 
! ------ np=32 -----
! Time to call mpi_barrier is             84.6600532531738281 microseconds
! Time to call send/recv barrier is      564.1698837280273437 microseconds
! Time to call bcast barrier is          986.8311882019041832 microseconds
! Time to call scatter barrier is         69.7302818298339844 microseconds
