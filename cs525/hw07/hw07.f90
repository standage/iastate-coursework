! Daniel S. Standage
! Iowa State University
! March 24, 2011
! HW 7, CS 525
use mpi
implicit none

! Parameters
integer, parameter :: n = 32
integer, parameter :: ltr = 528 ! # elements in lower tri portion = n(n+1)/2
!integer, parameter :: n = 1024
!integer, parameter :: ltr = 524800 ! # elements in lower tri portion = n(n+1)/2
integer, parameter :: ntrial = 3
integer, parameter :: dp = mpi_double_precision, comm=mpi_comm_world

! Data
double precision :: A(n,n)
double precision :: temp(ltr)
double precision :: t0, t1, t, maxtime
integer :: i,j,k,counter,ierror,p,rank,status(mpi_status_size)
integer :: blocklengths(n), displacements(n)
integer :: ltriangle

! Initialize MPI
call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

! Initialize data
A = 0
if(rank == 0) then
  call random_number(A)
  write(*, '(A, I0, A, I0)') "! p=", p, ', n=', n
endif


! Send data: user packing
do k = 1,ntrial
  call mpi_barrier(comm, ierror)
  t0 = mpi_wtime()
  if(rank == 0) then
    counter = 1;
    do i = 1,n
      do j = i,n
        temp(counter) = A(j, i)
        counter = counter + 1
      enddo
    enddo
  endif
  call mpi_bcast(temp, ltr, dp, 0, comm, ierror)
  counter = 1
  if(rank > 0) then
    do i = 1,n
      do j = i,n
        A(j, i) = temp(counter)
        counter = counter + 1
      enddo
    enddo
  endif
  t1 = mpi_wtime()
  t = t1 - t0
  call mpi_reduce(t, maxtime, 1, dp, mpi_max, 0, comm, ierror)
  if(rank == 0) then
    write(*, '(A, ES12.6, A)') "!    User packing: ", maxtime, " seconds"
  endif
enddo

! Send data: MPI derived types
do i = 0,n-1
  blocklengths(i+1) = n-i
  displacements(i+1) = i*(n+1)
enddo
call mpi_type_indexed(n, blocklengths, displacements, dp, ltriangle, ierror)
call mpi_type_commit(ltriangle, ierror)
do k = 1,ntrial
  call mpi_barrier(comm, ierror)
  t0 = mpi_wtime()
  call mpi_bcast(A(1,1), 1, ltriangle, 0, comm, ierror)
  t1 = mpi_wtime()
  t = t1 - t0
  call mpi_reduce(t, maxtime, 1, dp, mpi_max, 0, comm, ierror)
  if(rank == 0) then
    write(*, '(A, ES12.6, A)') "!   Derived types: ", maxtime, " seconds"
  endif
enddo

! Terminate
call mpi_type_free(ltriangle, ierror)
call mpi_finalize(ierror)
end

! p=8, n=32
!    User packing: 1.409054E-04 seconds
!    User packing: 1.020432E-04 seconds
!    User packing: 9.012222E-05 seconds
!   Derived types: 1.578331E-04 seconds
!   Derived types: 1.020432E-04 seconds
!   Derived types: 1.010895E-04 seconds
! p=16, n=32
!    User packing: 1.819134E-04 seconds
!    User packing: 1.468658E-04 seconds
!    User packing: 1.418591E-04 seconds
!   Derived types: 1.888275E-04 seconds
!   Derived types: 1.530647E-04 seconds
!   Derived types: 1.552105E-04 seconds
! p=32, n=32
!    User packing: 2.231598E-04 seconds
!    User packing: 2.119541E-04 seconds
!    User packing: 1.947880E-04 seconds
!   Derived types: 2.520084E-04 seconds
!   Derived types: 2.059937E-04 seconds
!   Derived types: 2.450943E-04 seconds
! p=8, n=1024
!    User packing: 6.397414E-02 seconds
!    User packing: 5.315089E-02 seconds
!    User packing: 5.315709E-02 seconds
!   Derived types: 8.330798E-02 seconds
!   Derived types: 8.160019E-02 seconds
!   Derived types: 8.077407E-02 seconds
! p=16, n=1024
!    User packing: 7.407188E-02 seconds
!    User packing: 5.927992E-02 seconds
!    User packing: 5.894709E-02 seconds
!   Derived types: 8.968592E-02 seconds
!   Derived types: 9.126091E-02 seconds
!   Derived types: 8.948207E-02 seconds
! p=32, n=1024
!    User packing: 8.381987E-02 seconds
!    User packing: 6.545711E-02 seconds
!    User packing: 6.472993E-02 seconds
!   Derived types: 1.014378E-01 seconds
!   Derived types: 1.011369E-01 seconds
!   Derived types: 1.008551E-01 seconds
