!########################################################################
!#   This is a parallel program for CPSC424/524 Assignment #3           #
!#   Task 2. Gravitational N-Body Simulation                            #
!#           with Collective Communication                              #
!#   Author: Xin Yan                                                    #
!#   Credit: Part of the program is adapted from a serial program       #
!#           by Andrew Sherman                                          #
!#   Date: 10/18/2014                                                   #
!########################################################################

!### MODULES ###
! A module that defines the type of each body
MODULE typedef
IMPLICIT NONE
  ! define a new type for each body
  TYPE :: body
    DOUBLE PRECISION :: p(3)
    DOUBLE PRECISION :: v(3)
    DOUBLE PRECISION :: m
  ENDTYPE
  ! define a new MPI_bd type for communication
CONTAINS
  SUBROUTINE mpistruct(MPI_bd)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
    INTEGER IERROR
    INTEGER :: MPI_bd ! new MPI type for the body
    INTEGER :: cnt(3) = (/3,3,1/) ! 
    INTEGER :: disp(3)
    INTEGER :: mpityp(3) = (/MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE/)
    INTEGER :: i, mpisz 
    disp(:) = 0
    DO i=1,2
      CALL MPI_TYPE_SIZE(mpityp(i),mpisz,IERROR)
      disp(i+1) = disp(i) + cnt(i)*mpisz 
    ENDDO
    CALL MPI_TYPE_STRUCT(3,cnt,disp,mpityp,MPI_bd,IERROR)
    CALL MPI_TYPE_COMMIT(MPI_bd,IERROR)
    RETURN
  END SUBROUTINE mpistruct 
END MODULE typedef

! A module that defines all the constants
MODULE const
IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: G = 1.0
  INTEGER, PARAMETER :: MXB = 100000 ! max # of pairwise interactions
  INTEGER, PARAMETER :: MXBCUT = 5000 ! max # of bodies that are within
                                      ! 5DU of another octant
  INTEGER, PARAMETER :: NOCT = 8
  DOUBLE PRECISION, PARAMETER :: CUT = 5.0 ! cutoff distance
END MODULE const

! A module that declares all vars for each worker
MODULE local
USE typedef
USE const
IMPLICIT NONE
  INTEGER :: MPI_bd
  INTEGER :: nts, nt
  DOUBLE PRECISION :: dt
  TYPE(body), ALLOCATABLE :: bdloc(:)
  TYPE(body), ALLOCATABLE :: sbuf(:) ! sendbuf for body update
  TYPE(body) :: avgloc ! store the local avg p, v and mtotloc
  TYPE(body) :: sbdcut(MXBCUT,0:NOCT-1) ! sendbuf for bodies within 5DU of another
  TYPE(body), ALLOCATABLE :: bdcut(:) ! buffer for bodies within 5DU of another 
  INTEGER :: nbloc, nbcut
  DOUBLE PRECISION :: dt2, r, r2
  DOUBLE PRECISION :: deltaf(3)
  DOUBLE PRECISION, ALLOCATABLE :: f(:,:)
  DOUBLE PRECISION :: ax, ay, az, vavgx, vavgy, vavgz, mx, my, mz
  DOUBLE PRECISION :: t0loc, t1loc, cptloc, wctloc ! local timing
  INTEGER :: i, j, thisbody, otherbody
  INTEGER :: IERROR
  INTEGER :: scnt(0:NOCT-1), rcnt(0:NOCT-1)
  INTEGER :: sdisp(0:NOCT-1), rdisp(0:NOCT-1)
END MODULE local

! A module that declares all vars for the master
MODULE mast
USE typedef
USE const
IMPLICIT NONE
  INTEGER :: nb
  TYPE(body), ALLOCATABLE :: bd(:) ! buffer for all bodies in the system
  DOUBLE PRECISION :: vavg(3), pcom(3) ! global average
  DOUBLE PRECISION :: time, time0, time1, cptime, wctime ! global timing
  DOUBLE PRECISION :: mtotal ! total mass of the system
  TYPE(body) :: avg(NOCT) ! recvbuf for local avg
  INTEGER :: INFILE = 5, OUTFILE = 6 ! in/out unit
  DOUBLE PRECISION :: wct(0:NOCT-1) !recvbuf for local timing
  INTEGER :: nbd(0:NOCT-1) ! recvbuf for local bdcnt
END MODULE mast

!### Main program ###
PROGRAM task2
IMPLICIT NONE 
INCLUDE 'mpif.h'
  INTEGER IERROR, iproc
  ! MPI Init and Branching
  CALL MPI_INIT(IERROR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, IERROR)
  IF (iproc .EQ. 0) THEN
    CALL master(iproc)
  ELSE 
    CALL worker(iproc)
  ENDIF
  ! MPI Finalize
  CALL MPI_FINALIZE(IERROR)
STOP
END PROGRAM task2

!### Subroutine for master ###
SUBROUTINE master(iproc)
USE typedef
USE const
USE local
USE mast
IMPLICIT NONE 
INCLUDE 'mpif.h'
  INTEGER iproc 
  INTEGER tmp
  ! Read data from stdin
  READ(INFILE,*) nb
  READ(INFILE,*) nts
  READ(INFILE,*) dt
  IF (nb .EQ. 0) THEN
     WRITE(OUTFILE,*) 'Nothing to do'
     STOP
  ENDIF
  WRITE(OUTFILE,*) 'N-Body Problem (Parallel Run) N = ', nb
  ALLOCATE(bd(nb))
  READ(INFILE,*) (bd(i)%m, i=1,nb)
  READ(INFILE,*) ((bd(i)%p(j), j=1,3),i=1,nb)
  READ(INFILE,*) ((bd(i)%v(j), j=1,3),i=1,nb)
  ! Calc mtotal
  mtotal = 0.
  DO i=1,nb
    mtotal = mtotal + bd(i)%m
  ENDDO
  ! Assign global data to workers
  CALL MPI_BCAST(nb,1,MPI_INT,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(nts,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  ! Sort init bd array according to oct index
  CALL semiqsort(bd,nb,1,nb,3) 
  CALL cntbd(nb,bd,scnt)
  CALL mpistruct(MPI_bd)
  CALL calcdisp(NOCT,scnt,sdisp)
  ! Reorder to make the octant index rotate counterclockwise
  tmp = scnt(2)
  scnt(2) = scnt(3)
  scnt(3) = tmp
  tmp = scnt(6)
  scnt(6) = scnt(7)
  scnt(7) = tmp
  tmp = sdisp(2)
  sdisp(2) = sdisp(3)
  sdisp(3) = tmp
  tmp = sdisp(6)
  sdisp(6) = sdisp(7)
  sdisp(7) = tmp
  ! Assign init data of the bodies to workers
  CALL MPI_SCATTER(scnt,1,MPI_INT,nbloc,1,MPI_INT,0,MPI_COMM_WORLD,IERROR) 
  ALLOCATE(bdloc(nbloc)) ! set up body recvbuf
  CALL MPI_SCATTERV(bd,scnt,sdisp,MPI_bd,bdloc,nbloc,MPI_bd,0, &
                MPI_COMM_WORLD,IERROR)
  ! Master frees the space occupied by bd
  DEALLOCATE(bd)
  ! Calculate total local mass
  avgloc%m = 0.
  DO i=1,nbloc
    avgloc%m = avgloc%m + bdloc(i)%m
  ENDDO
  dt2 = dt/2.
  ! Start global timing
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  CALL timing(time0,cptime)
  ! Start local timing
!  t0loc = time0
!  cptloc = cptime
  ! Timestep loop
  DO nt=0,nts
    ! Output if needed
    IF ((MOD(nt,128) .EQ. 0) .OR. (nt .EQ. nts)) THEN
      ! Calculate local average
      avgloc%p(:) = 0.
      avgloc%v(:) = 0.
      DO i=1,nbloc
        DO j=1,3
          avgloc%p(j) = avgloc%p(j) + bdloc(i)%m*bdloc(i)%p(j)
          avgloc%v(j) = avgloc%v(j) + bdloc(i)%v(j)
        ENDDO
      ENDDO
      ! Gather the local averages from all workers
      CALL MPI_GATHER(avgloc,1,MPI_bd,avg,1,MPI_bd,0,MPI_COMM_WORLD,IERROR)
      ! Calculate global average  
      pcom(:) = 0.
      vavg(:) = 0.
      DO i=1,NOCT
        DO j=1,3
          pcom(j) = pcom(j) + avg(i)%p(j)
          vavg(j) = vavg(j) + avg(i)%v(j)
        ENDDO
      ENDDO
      DO j=1,3
        pcom(j) = pcom(j)/mtotal
        vavg(j) = vavg(j)/nb
      ENDDO
      ! Output global average
      time = nt * dt
      IF(nt .EQ. 0) THEN
        WRITE(OUTFILE,990) time, pcom(1), pcom(2), pcom(3), &
                   vavg(1), vavg(2), vavg(3)
      ELSE 
        WRITE(OUTFILE,991) nt, time, pcom(1), pcom(2), pcom(3), &
                   vavg(1), vavg(2), vavg(3)
      ENDIF
 990  FORMAT(/' Initial Conditions (time = ', g16.8, '):'// &
            '     Center of Mass:   (', &
                     g16.8, ', ', g16.8, ', ', g16.8, ' )'/ &
            '     Average Velocity: (', &
                     g16.8, ', ', g16.8, ', ', g16.8, ' )' )
 991  FORMAT(/' Conditions after timestep ', i5, &
            ' (time = ', g16.8, '):'// &
            '     Center of Mass:   (',&
                     g16.8, ', ', g16.8, ', ', g16.8, ' )'/ &
            '     Average Velocity: (',&
                     g16.8, ', ', g16.8, ', ', g16.8, ' )' )
    ENDIF     
    ! Computation
    IF (nt .LT. nts) THEN
      ! Determines which bodies are within 5DU of other octant
      scnt(:) = MXBCUT
      CALL calcdisp(NOCT,scnt,sdisp)
      CALL chkcut(iproc,bdloc,nbloc,sbdcut,scnt)
      CALL MPI_ALLTOALL(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD, &
                      IERROR)
      nbcut = 0
      DO i=0,NOCT-1
        nbcut = nbcut + rcnt(i)
      ENDDO
      ALLOCATE(bdcut(nbcut))
      CALL calcdisp(NOCT,rcnt,rdisp)
      CALL MPI_ALLTOALLV(sbdcut,scnt,sdisp,MPI_bd,bdcut,rcnt,rdisp,MPI_bd, &
                        MPI_COMM_WORLD,IERROR)
    ! Initialize forces
      ALLOCATE(f(nbloc,3))
      DO thisbody=1,nbloc
        f(thisbody,:) = 0.
      ENDDO
      DO thisbody=1,nbloc
    ! Compute all pairwise interbody forces
        DO otherbody=thisbody+1,nbloc
          CALL force(nbloc,nbloc,thisbody,otherbody,deltaf,bdloc,bdloc)
          f(thisbody,1) = f(thisbody,1) + deltaf(1)
          f(thisbody,2) = f(thisbody,2) + deltaf(2)
          f(thisbody,3) = f(thisbody,3) + deltaf(3)
          f(otherbody,1) = f(otherbody,1) - deltaf(1)
          f(otherbody,2) = f(otherbody,2) - deltaf(2)
          f(otherbody,3) = f(otherbody,3) - deltaf(3)
        ENDDO
    ! Compute forces between bodies that belong to different octants
        DO otherbody=1,nbcut
          CALL force(nbloc,nbcut,thisbody,otherbody,deltaf,bdloc,bdcut)
          f(thisbody,1) = f(thisbody,1) + deltaf(1)
          f(thisbody,2) = f(thisbody,2) + deltaf(2)
          f(thisbody,3) = f(thisbody,3) + deltaf(3)
        ENDDO
      ENDDO
    ! Move the bodies
      DO thisbody=1,nbloc
        ax = f(thisbody,1)/bdloc(thisbody)%m
        ay = f(thisbody,2)/bdloc(thisbody)%m
        az = f(thisbody,3)/bdloc(thisbody)%m
        vavgx = bdloc(thisbody)%v(1) + dt2*ax
        vavgy = bdloc(thisbody)%v(2) + dt2*ay
        vavgz = bdloc(thisbody)%v(3) + dt2*az
        bdloc(thisbody)%p(1) = bdloc(thisbody)%p(1) + dt*vavgx
        bdloc(thisbody)%p(2) = bdloc(thisbody)%p(2) + dt*vavgy
        bdloc(thisbody)%p(3) = bdloc(thisbody)%p(3) + dt*vavgz
        bdloc(thisbody)%v(1) = bdloc(thisbody)%v(1) + dt*ax
        bdloc(thisbody)%v(2) = bdloc(thisbody)%v(2) + dt*ay
        bdloc(thisbody)%v(3) = bdloc(thisbody)%v(3) + dt*az
      ENDDO
      ! Update ownership list
      scnt(:) = 0
      CALL semiqsort(bdloc,nbloc,1,nbloc,3) 
      CALL cntbd(nbloc,bdloc,scnt)
      CALL calcdisp(NOCT,scnt,sdisp)
      ! Reorder to make the octant index rotate counterclockwise
      tmp = scnt(2)
      scnt(2) = scnt(3)
      scnt(3) = tmp
      tmp = scnt(6)
      scnt(6) = scnt(7)
      scnt(7) = tmp
      tmp = sdisp(2)
      sdisp(2) = sdisp(3)
      sdisp(3) = tmp
      tmp = sdisp(6)
      sdisp(6) = sdisp(7)
      sdisp(7) = tmp
      CALL MPI_ALLTOALL(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD, &
                  IERROR)
      CALL MPI_GATHER(nbloc,1,MPI_INT,nbd,1,MPI_INT,0,MPI_COMM_WORLD,IERROR)
      ALLOCATE(sbuf(nbloc))
      sbuf(1:nbloc) =  bdloc(1:nbloc)
      DEALLOCATE(bdloc)
      nbloc = 0
      DO i=0,NOCT-1
        nbloc = nbloc + rcnt(i)
      ENDDO
      ALLOCATE(bdloc(nbloc))
      CALL calcdisp(NOCT,rcnt,rdisp)
      CALL MPI_ALLTOALLV(sbuf,scnt,sdisp,MPI_bd,bdloc,rcnt,rdisp,MPI_bd, &
                        MPI_COMM_WORLD,IERROR)
      ! Deallocate bdcut
      DEALLOCATE(f)
      DEALLOCATE(bdcut)
      DEALLOCATE(sbuf)
    ENDIF
  ENDDO
  ! Stop timing
  ! calculate local elapsed time
!  CALL timing(t1loc,cptloc)
!  wctloc = t1loc - t0loc
!  CALL MPI_GATHER(wctloc,1,MPI_DOUBLE,wct,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  ! calculate global elapsed time
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  CALL timing(time1,cptime)
  wctime = time1 - time0
  ! Output timing result
  WRITE(OUTFILE,992) nts, nb, wctime
 992  FORMAT(/' Elapsed time for ', I5, ' timesteps with ', &
              I5, ' bodies: ', F9.4, ' seconds')
!  WRITE(OUTFILE,*) 'Summary for each processor:'
!  DO i=0,NOCT-1
!    WRITE(OUTFILE,993) i, nbd(i), wct(i)
!  ENDDO
! 993  FORMAT(/' Elapsed time for ', I5, ' processor with ', &
!              I5, ' bodies: ', F9.4, ' seconds')
  RETURN
END SUBROUTINE master

! ### Subroutine for workers ###
SUBROUTINE worker(iproc)
USE typedef
USE const
USE local
USE mast
IMPLICIT NONE
INCLUDE 'mpif.h'
  INTEGER iproc
  INTEGER tmp
  ! collect global data
  CALL MPI_BCAST(nb,1,MPI_INT,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(nts,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  ! collect init data of bodies from the master
  CALL MPI_SCATTER(scnt,1,MPI_INT,nbloc,1,MPI_INT,0,MPI_COMM_WORLD,IERROR) 
  ALLOCATE(bd(nb)) ! set up body sendbuf
  ALLOCATE(bdloc(nbloc)) ! set up body recvbuf
  CALL mpistruct(MPI_bd)
  CALL MPI_SCATTERV(bd,scnt,sdisp,MPI_bd,bdloc,nbloc,MPI_bd,0, &
                MPI_COMM_WORLD,IERROR)
  DEALLOCATE(bd)
  ! Calculate total local mass
  avgloc%m = 0.
  DO i=1,nbloc
    avgloc%m = avgloc%m + bdloc(i)%m
  ENDDO
  dt2 = dt/2.
  ! Start global timing
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  ! Start local timing 
!  CALL timing(t0loc,cptloc)
  ! Timestep loop
  DO nt=0,nts
    ! Output if needed
    IF ((MOD(nt,128) .EQ. 0) .OR. (nt .EQ. nts)) THEN
      ! Calculate local average
      avgloc%p(:) = 0.
      avgloc%v(:) = 0.
      DO i=1,nbloc
        DO j=1,3
          avgloc%p(j) = avgloc%p(j) + bdloc(i)%m*bdloc(i)%p(j)
          avgloc%v(j) = avgloc%v(j) + bdloc(i)%v(j)
        ENDDO
      ENDDO
      ! Gather the local averages from all workers
      CALL MPI_GATHER(avgloc,1,MPI_bd,avg,1,MPI_bd,0,MPI_COMM_WORLD,IERROR)
    ENDIF     
    ! Computation
    IF (nt .LT. nts) THEN
      ! Determines which bodies are within 5DU of other octant
      scnt(:) = MXBCUT
      CALL calcdisp(NOCT,scnt,sdisp)
      CALL chkcut(iproc,bdloc,nbloc,sbdcut,scnt)
      CALL MPI_ALLTOALL(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD, &
                      IERROR)
      nbcut = 0
      DO i=0,NOCT-1
        nbcut = nbcut + rcnt(i)
      ENDDO
      ALLOCATE(bdcut(nbcut))
      CALL calcdisp(NOCT,rcnt,rdisp)
      CALL MPI_ALLTOALLV(sbdcut,scnt,sdisp,MPI_bd,bdcut,rcnt,rdisp,MPI_bd, &
                        MPI_COMM_WORLD,IERROR)
    ! Initialize forces
      ALLOCATE(f(nbloc,3))
      DO thisbody=1,nbloc
        f(thisbody,:) = 0.
      ENDDO
      DO thisbody=1,nbloc
    ! Compute all pairwise interbody forces
        DO otherbody=thisbody+1,nbloc
          CALL force(nbloc,nbloc,thisbody,otherbody,deltaf,bdloc,bdloc)
          f(thisbody,1) = f(thisbody,1) + deltaf(1)
          f(thisbody,2) = f(thisbody,2) + deltaf(2)
          f(thisbody,3) = f(thisbody,3) + deltaf(3)
          f(otherbody,1) = f(otherbody,1) - deltaf(1)
          f(otherbody,2) = f(otherbody,2) - deltaf(2)
          f(otherbody,3) = f(otherbody,3) - deltaf(3)
        ENDDO
    ! Compute forces between bodies that belong to different octants
        DO otherbody=1,nbcut
          CALL force(nbloc,nbcut,thisbody,otherbody,deltaf,bdloc,bdcut)
          f(thisbody,1) = f(thisbody,1) + deltaf(1)
          f(thisbody,2) = f(thisbody,2) + deltaf(2)
          f(thisbody,3) = f(thisbody,3) + deltaf(3)
        ENDDO
      ENDDO
    ! Move the bodies
      DO thisbody=1,nbloc
        ax = f(thisbody,1)/bdloc(thisbody)%m
        ay = f(thisbody,2)/bdloc(thisbody)%m
        az = f(thisbody,3)/bdloc(thisbody)%m
        vavgx = bdloc(thisbody)%v(1) + dt2*ax
        vavgy = bdloc(thisbody)%v(2) + dt2*ay
        vavgz = bdloc(thisbody)%v(3) + dt2*az
        bdloc(thisbody)%p(1) = bdloc(thisbody)%p(1) + dt*vavgx
        bdloc(thisbody)%p(2) = bdloc(thisbody)%p(2) + dt*vavgy
        bdloc(thisbody)%p(3) = bdloc(thisbody)%p(3) + dt*vavgz
        bdloc(thisbody)%v(1) = bdloc(thisbody)%v(1) + dt*ax
        bdloc(thisbody)%v(2) = bdloc(thisbody)%v(2) + dt*ay
        bdloc(thisbody)%v(3) = bdloc(thisbody)%v(3) + dt*az
      ENDDO
      ! Update ownership list
      scnt(:) = 0
      CALL semiqsort(bdloc,nbloc,1,nbloc,3) 
      CALL cntbd(nbloc,bdloc,scnt)
      CALL calcdisp(NOCT,scnt,sdisp)
      ! Reorder to make the octant index rotate counterclockwise
      tmp = scnt(2)
      scnt(2) = scnt(3)
      scnt(3) = tmp
      tmp = scnt(6)
      scnt(6) = scnt(7)
      scnt(7) = tmp
      tmp = sdisp(2)
      sdisp(2) = sdisp(3)
      sdisp(3) = tmp
      tmp = sdisp(6)
      sdisp(6) = sdisp(7)
      sdisp(7) = tmp
      CALL MPI_ALLTOALL(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD, &
                  IERROR)
      CALL MPI_GATHER(nbloc,1,MPI_INT,nbd,1,MPI_INT,0,MPI_COMM_WORLD,IERROR)
      ALLOCATE(sbuf(nbloc))
      sbuf(1:nbloc) = bdloc(1:nbloc)
      DEALLOCATE(bdloc)
      nbloc = 0
      DO i=0,NOCT-1
        nbloc = nbloc + rcnt(i)
      ENDDO
      ALLOCATE(bdloc(nbloc))
      CALL calcdisp(NOCT,rcnt,rdisp)
      CALL MPI_ALLTOALLV(sbuf,scnt,sdisp,MPI_bd,bdloc,rcnt,rdisp,MPI_bd, &
                        MPI_COMM_WORLD,IERROR)
      ! Deallocate bdcut
      DEALLOCATE(f)
      DEALLOCATE(bdcut)
      DEALLOCATE(sbuf)
    ENDIF
  ENDDO
  ! Stop timing
  ! calculate local elapsed time
!  CALL timing(t1loc,cptloc)
!  wctloc = t1loc - t0loc
!  CALL MPI_GATHER(wctloc,1,MPI_DOUBLE,wct,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERROR)
  DEALLOCATE(bdloc)
  ! calculate global elapsed time
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  RETURN
END SUBROUTINE worker

! Subroutine for sorting the body array to different octants 
RECURSIVE SUBROUTINE semiqsort(A,N,S,E,id)
USE typedef
USE const
IMPLICIT NONE
  INTEGER :: N
  TYPE(body) :: A(N)
  INTEGER :: S, E, L, R
  INTEGER :: K
  TYPE(body) :: temp
  INTEGER :: id ! current sorting dimension (x=1,y=2,z=3)
  L = S - 1
  R = E + 1
  IF ((id .LE. 0).OR.(S.GE.E)) RETURN
  K = 0 !sorting by value of K
  DO WHILE(.TRUE.)
    DO WHILE(.TRUE.)
      L = L + 1
      IF ((A(L)%p(id) .GT. K) .OR. (L .GE. E)) EXIT 
    ENDDO
    DO WHILE(.TRUE.)
      R = R - 1
      IF ((A(R)%p(id) .LT. K) .OR. (R .LE. S)) EXIT
    ENDDO
    IF(R .LE. L) EXIT
    temp = A(L)
    A(L) = A(R)
    A(R) = temp
  ENDDO
  CALL semiqsort(A,N,S,R,id-1)
  CALL semiqsort(A,N,R+1,E,id-1)
  RETURN
END SUBROUTINE semiqsort

! Subroutine to count the number of bodies in each octant
SUBROUTINE cntbd(n,bdary,cnt)
USE typedef
USE const
IMPLICIT NONE
  INTEGER :: n, i
  TYPE(body) :: bdary(n)
  INTEGER :: cnt(0:NOCT-1)
  cnt(:) = 0
  DO i=1,n
    IF(bdary(i)%p(3) .le. 0) THEN
      IF(bdary(i)%p(2) .le. 0) THEN
        IF(bdary(i)%p(1) .le. 0) THEN
          cnt(0) = cnt(0) + 1
        ELSE
          cnt(1) = cnt(1) + 1
        ENDIF
      ELSE
        IF(bdary(i)%p(1) .le. 0) THEN
          cnt(2) = cnt(2) + 1
        ELSE
          cnt(3) = cnt(3) + 1
        ENDIF
      ENDIF
    ELSE
      IF(bdary(i)%p(2) .le. 0) THEN
        IF(bdary(i)%p(1) .le. 0) THEN
          cnt(4) = cnt(4) + 1
        ELSE
          cnt(5) = cnt(5) + 1
        ENDIF
      ELSE
        IF(bdary(i)%p(1) .le. 0) THEN
          cnt(6) = cnt(6) + 1
        ELSE
          cnt(7) = cnt(7) + 1
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  RETURN
END SUBROUTINE cntbd
  
! subroutine to calculate displacement
SUBROUTINE calcdisp(n,blen,disp)
IMPLICIT NONE
  INTEGER :: i, n
  INTEGER :: disp(n), blen(n)
  INTEGER IERROR
  disp(:) = 0
  DO i=1,n-1
    disp(i+1) = disp(i) + blen(i) 
  ENDDO
  RETURN
END SUBROUTINE calcdisp

! subroutine to return the body arrays that are within 5DU of octants
SUBROUTINE chkcut(rank,array,n,sbdcut,scnt)
USE typedef
USE const
IMPLICIT NONE
  INTEGER :: rank, ip, ibd, n, p(7)
  TYPE(body) :: sbdcut(MXBCUT,0:NOCT-1)
  INTEGER :: scnt(0:NOCT-1)
  DOUBLE PRECISION :: CUT2 = CUT*CUT 
  DOUBLE PRECISION :: x, y, z
  DOUBLE PRECISION :: dist(7)
  TYPE(body) :: array(n)
  scnt(:) = 0
  ! define pointers to other octants according to relative locations
  IF (MOD(rank,2) .EQ. 0) THEN
    p(1) = MOD((rank+1),4) + INT(rank/4)*4 !oct share y-z plane
    p(2) = MOD((rank+2),4) + INT(rank/4)*4 !oct share z axis
    p(3) = MOD((rank+3),4) + INT(rank/4)*4 !oct share x-z plane
    p(4) = MOD((rank+4),8) !oct share x-y plane
    p(5) = MOD((rank+5),8) !oct share y axis
    p(6) = 6 - rank !oct share origin
    p(7) = 7 - rank !oct share x axis
  ELSE 
    p(1) = MOD((rank+3),4) + INT(rank/4)*4
    p(2) = MOD((rank+2),4) + INT(rank/4)*4
    p(3) = MOD((rank+1),4) + INT(rank/4)*4
    p(4) = MOD((rank+4),8)
    p(5) = MOD((rank+3),8)
    p(6) = 8 - rank 
    p(7) = 7 - rank 
  ENDIF
!! Debug start
!  write(*,*) rank, 'p:', p(:)
!! debug end
  DO ibd=1,n
      x = array(ibd)%p(1)
      y = array(ibd)%p(2)
      z = array(ibd)%p(3)
      dist(1) = x*x !dist to y-z plane 
      dist(3) = y*y !dist to x-z plane
      dist(4) = z*z !dist to x-y plane
      dist(2) = dist(1) + dist(3) !dist to z axis
      dist(5) = dist(1) + dist(4) !dist to y axis
      dist(6) = dist(1) + dist(3) + dist(4) !dist to origin
      dist(7) = dist(3) + dist(4) !dist to x axis
    DO ip = 1,7
      IF(dist(ip) .LE. CUT2) THEN
        scnt(p(ip)) = scnt(p(ip)) + 1
        sbdcut(scnt(p(ip)),p(ip)) = array(ibd)
      ENDIF
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE chkcut

! Subroutine for force calculation    
SUBROUTINE force(N1, N2, body1, body2, deltaf, bdary1, bdary2)
USE typedef
USE const
IMPLICIT NONE
  INTEGER body1, body2
  INTEGER :: N1, N2
  DOUBLE PRECISION deltaf(3)
  TYPE(body) :: bdary1(N1), bdary2(N2)
  DOUBLE PRECISION :: gmmr3, r, r2, dx, dy, dz
  dx = bdary2(body2)%p(1) - bdary1(body1)%p(1)
  dy = bdary2(body2)%p(2) - bdary1(body1)%p(2)
  dz = bdary2(body2)%p(3) - bdary1(body1)%p(3)
  r2 = dx*dx + dy*dy + dz*dz
  r = SQRT(r2)
  IF (r .LE. CUT) THEN
     gmmr3 = G * bdary1(body1)%m * bdary2(body2)%m / (r2 * r)
     deltaf(1) = gmmr3 * dx
     deltaf(2) = gmmr3 * dy
     deltaf(3) = gmmr3 * dz
  ELSE
     deltaf(1) = 0.
     deltaf(2) = 0.
     deltaf(3) = 0.
  ENDIF
  RETURN
END SUBROUTINE force


