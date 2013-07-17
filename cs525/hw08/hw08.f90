! Daniel S. Standage
! Iowa State University
! April 7, 2011
! HW 8, CS 525
use mpi
implicit none

integer, parameter :: n=4096, niter=50, tag=1
integer :: i, j, m, iter, ierror, left, right
double precision :: t, time, max_time
double precision,allocatable :: A(:,:), B(:,:)
integer, parameter :: dp=mpi_double_precision, comm = mpi_comm_world
integer :: p, rank, status(mpi_status_size)

call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Compute the size of the local blocks
m = n/p
if (rank < (n-p*m)) then
m = m + 1
endif ! if n = m*p + r, 0 < r < p, then m = n/p + 1 for
      ! rank = 0, 1, ..., r-1.
      
! Allocate the local arrays A and B
allocate (A(0:n+1,0:m+1), B(n,m))

! Initialize A
A = 0.d0
if (rank == 0) then
  A(1:n,0) = 1.d0
endif

! Use mpi_proc_null to determine the left and right neighbors when
! dest or source = mpi_proc_null, then the routine does nothing.
if (rank == 0) then
  left = mpi_proc_null
else
  left = rank - 1
endif
if (rank == p-1) then
  right = mpi_proc_null
else
  right = rank + 1
endif

! Parallel Jacobi Iteration
call mpi_barrier(comm, ierror)
t = mpi_wtime()
do iter = 1, niter
  do j = 1, m
    do i = 1, n
      B(i,j) = 0.25d0*(A(i-1,j)+A(i+1,j) + A(i,j-1) +A(i,j+1))
    enddo
  enddo
  do j = 1, m
    do i = 1, n
      A(i,j) = B(i,j)
    enddo
  enddo
  ! Perform the communication between processors
  call mpi_sendrecv( B(1,m), n, dp, right, 0, A(1,m+1), n, dp, right, 0, comm, mpi_status_ignore, ierror )
  call mpi_sendrecv( B(1,1), n, dp, left, 0, A(1,0), n, dp, left, 0, comm, mpi_status_ignore, ierror )
enddo
time = mpi_wtime() - t

call mpi_barrier(comm, ierror) ! not sure this is needed
call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
if (rank == 0) then
  !print *,'Jacobi Iteration using mpi_sendrecv'
  print *,'! n, niter, and p = ', n, niter, p
  print *,'! time for niter iterations = ',time,' seconds'
  print *,'! average time per iteration = ', max_time/dble(niter),' seconds'
  print *,'! prevent dead code elimination, A(2,2) = ', A(2,2)
  print *,'!'
endif

call mpi_finalize(ierror)
end

! [[[ Discussion ]]]
! The scalability of the iteration approaches the optimum: runtime cut in half
! when the number of processors doubles. This does not seem to change with the
! number of iterations: the same pattern (and almost the same exact runtimes per
! iteration) is seen whether niter is 50 or 300.
!
! [[[ Output ]]]
! n, niter, and p =         4096         300           4
! time for niter iterations =    85.8410379886627       seconds
! average time per iteration =   0.286136793295542       seconds
! prevent dead code elimination, A(2,2) =   0.491609681926248      
!
! n, niter, and p =         4096         300           8
! time for niter iterations =    43.1221950054169       seconds
! average time per iteration =   0.143740650018056       seconds
! prevent dead code elimination, A(2,2) =   0.491609681926248      
!
! n, niter, and p =         4096         300          16
! time for niter iterations =    21.7420301437378       seconds
! average time per iteration =   7.247343381245931E-002  seconds
! prevent dead code elimination, A(2,2) =   0.491609681926248      
!
! n, niter, and p =         4096         300          32
! time for niter iterations =    10.9497749805450       seconds
! average time per iteration =   3.649924993515015E-002  seconds
! prevent dead code elimination, A(2,2) =   0.491609681926248      
!
! n, niter, and p =         4096          50           4
! time for niter iterations =    14.3591160774231       seconds
! average time per iteration =   0.287182722091675       seconds
! prevent dead code elimination, A(2,2) =   0.452413635350912     
!
! n, niter, and p =         4096          50           8
! time for niter iterations =    7.18668198585510       seconds
! average time per iteration =   0.143733639717102       seconds
! prevent dead code elimination, A(2,2) =   0.452413635350912     
!
! n, niter, and p =         4096          50          16
! time for niter iterations =    3.63717508316040       seconds
! average time per iteration =   7.274350166320800E-002  seconds
! prevent dead code elimination, A(2,2) =   0.452413635350912     
!
! n, niter, and p =         4096          50          32
! time for niter iterations =    1.84469008445740       seconds
! average time per iteration =   3.689380168914795E-002  seconds
! prevent dead code elimination, A(2,2) =   0.452413635350912     
!
