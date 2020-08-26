!  Copyright 1991-2020 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
!  Mumps Technologies, University of Bordeaux.
!
!  This is released under the CeCILL-C license:
!  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html
!
!  This file contains stub MPI library functions for
!  linking/running MUMPS on a platform where MPI is not installed.

module mpi

use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64
use, intrinsic :: iso_c_binding, only : c_double
implicit none
public

integer, parameter :: MPI_STATUS_SIZE=2
integer :: mpi_status_ignore(MPI_STATUS_SIZE), mpi_proc_null

INTEGER, parameter :: MPI_2DOUBLE_PRECISION=1
INTEGER, parameter ::  MPI_2INTEGER=2
INTEGER, parameter ::  MPI_2REAL=3
INTEGER, parameter ::  MPI_ANY_SOURCE=4
INTEGER, parameter ::  MPI_ANY_TAG=5
INTEGER, parameter ::  MPI_BYTE=6
INTEGER, parameter ::  MPI_CHARACTER=7
INTEGER, parameter ::  MPI_COMM_NULL=8
INTEGER, parameter ::  MPI_COMM_WORLD=9
INTEGER, parameter ::  MPI_COMPLEX=10
INTEGER, parameter ::  MPI_DOUBLE_COMPLEX=11
INTEGER, parameter ::  MPI_DOUBLE_PRECISION=12
INTEGER, parameter ::  MPI_INTEGER=13
INTEGER, parameter ::  MPI_LOGICAL=13
INTEGER, parameter ::  MPI_MAX=15
INTEGER, parameter ::  MPI_MAX_PROCESSOR_NAME=31
INTEGER, parameter ::  MPI_MAXLOC=16
INTEGER, parameter ::  MPI_MIN=17
INTEGER, parameter ::  MPI_MINLOC=18
INTEGER, parameter ::  MPI_PACKED=19
INTEGER, parameter ::  MPI_PROD=20
INTEGER, parameter ::  MPI_REAL=21
INTEGER, parameter ::  MPI_REPLACE=22
INTEGER, parameter ::  MPI_REQUEST_NULL=23
INTEGER, parameter ::  MPI_SOURCE=1
INTEGER, parameter ::  MPI_SUM=26
INTEGER, parameter ::  MPI_TAG=2
INTEGER, parameter ::  MPI_UNDEFINED=28
INTEGER, parameter ::  MPI_WTIME_IS_GLOBAL=30
INTEGER, parameter ::  MPI_LOR=31
INTEGER, parameter ::  MPI_LAND=32
INTEGER, parameter ::  MPI_INTEGER8=33
INTEGER, parameter ::  MPI_REAL8=34
INTEGER, parameter ::  MPI_BSEND_OVERHEAD=0

integer :: MPI_IN_PLACE

interface
subroutine mumps_elapse(val) bind(c)
import c_double
real(c_double), intent(inout) :: val
end subroutine mumps_elapse
end interface

contains

SUBROUTINE MPI_BSEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IERR )
INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
INTEGER BUF(*)
error stop 'MPI_BSEND should not be called.'
END SUBROUTINE MPI_BSEND

SUBROUTINE MPI_BUFFER_ATTACH(BUF, CNT,  IERR )
INTEGER CNT, IERR
INTEGER BUF(*)
IERR = 0
END SUBROUTINE MPI_BUFFER_ATTACH

SUBROUTINE MPI_BUFFER_DETACH(BUF, CNT,  IERR )
INTEGER CNT, IERR
INTEGER BUF(*)
IERR = 0
END SUBROUTINE MPI_BUFFER_DETACH

SUBROUTINE MPI_GATHER( SENDBUF, CNT, DATATYPE, RECVBUF, RECCNT, RECTYPE, ROOT, COMM, IERR )
INTEGER CNT, DATATYPE, RECCNT, RECTYPE, ROOT, COMM, IERR
class(*) :: SENDBUF(:), RECVBUF(:)

IF (RECCNT /= CNT) error stop 'ERROR in MPI_GATHER, RECCNT != CNT'

CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF)
IERR = 0
END SUBROUTINE MPI_GATHER

SUBROUTINE MPI_GATHERV( SENDBUF, CNT, DATATYPE, RECVBUF, RECCNT, DISPLS, RECTYPE, ROOT, COMM, IERR )
INTEGER CNT, DATATYPE, RECTYPE, ROOT, COMM, IERR
INTEGER RECCNT(1)
class(*) :: SENDBUF(:), RECVBUF(:)
INTEGER DISPLS(*)
!     Note that DISPLS is ignored in this version. One may
!     want to copy in reception buffer with a shift DISPLS(1).
!     This requires passing the offset DISPLS(1) to
!     "MUMPS_COPY_DATATYPE" routines.
IF (RECCNT(1) /= CNT) error stop 'ERROR in MPI_GATHERV, RECCNT(1) != CNT'
CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF )
IERR = 0
END SUBROUTINE MPI_GATHERV

SUBROUTINE MPI_ALLREDUCE( SENDBUF, RECVBUF, CNT, DATATYPE, OPERATION, COMM, IERR )
INTEGER CNT, DATATYPE, OPERATION, COMM, IERR
class(*) :: SENDBUF(:), RECVBUF(:)
IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, CNT))  CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF )
IERR = 0
END SUBROUTINE MPI_ALLREDUCE

SUBROUTINE MPI_REDUCE( SENDBUF, RECVBUF, CNT, DATATYPE, OP, ROOT, COMM, IERR )
INTEGER CNT, DATATYPE, OP, ROOT, COMM, IERR
class(*) :: SENDBUF(:), RECVBUF(:)
IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, CNT)) CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF )
IERR = 0
END SUBROUTINE MPI_REDUCE

SUBROUTINE MPI_REDUCE_SCATTER( SENDBUF, RECVBUF, RCVCNT,  DATATYPE, OP, COMM, IERR )
INTEGER RCVCNT, DATATYPE, OP, COMM, IERR
class(*) :: SENDBUF(:), RECVBUF(:)
IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, RCVCNT)) CALL MUMPS_COPY( RCVCNT, SENDBUF, RECVBUF )
IERR = 0
END SUBROUTINE MPI_REDUCE_SCATTER

SUBROUTINE MPI_ABORT( COMM, IERRCODE, IERR )
INTEGER COMM, IERRCODE, IERR
error stop "MPI_ABORT called"
END SUBROUTINE MPI_ABORT

SUBROUTINE MPI_ALLTOALL( SENDBUF, SENDCNT, SENDTYPE, RECVBUF, RECVCNT, RECVTYPE, COMM, IERR )
INTEGER SENDCNT, SENDTYPE, RECVCNT, RECVTYPE, COMM, IERR
class(*) :: SENDBUF(:), RECVBUF(:)
IF ( RECVCNT .NE. SENDCNT ) error stop 'ERROR in MPI_ALLTOALL, RECVCNT != SENDCNT'
if ( RECVTYPE .NE. SENDTYPE ) error stop 'ERROR in MPI_ALLTOALL, RECVTYPE != SENDTYPE'

CALL MUMPS_COPY( SENDCNT, SENDBUF, RECVBUF )
IERR = 0
END SUBROUTINE MPI_ALLTOALL

SUBROUTINE MPI_ATTR_PUT( COMM, KEY, VAL, IERR )
INTEGER COMM, KEY, VAL, IERR
END SUBROUTINE MPI_ATTR_PUT

SUBROUTINE MPI_BARRIER( COMM, IERR )
INTEGER COMM, IERR
IERR = 0
END SUBROUTINE MPI_BARRIER

SUBROUTINE MPI_GET_PROCESSOR_NAME( NAME, RESULTLEN, IERROR)
CHARACTER (LEN=*) NAME
INTEGER RESULTLEN,IERROR
RESULTLEN = 1
IERROR = 0
NAME = 'X'
END SUBROUTINE MPI_GET_PROCESSOR_NAME

SUBROUTINE MPI_BCAST( BUFFER, CNT, DATATYPE, ROOT, COMM, IERR )
INTEGER CNT, DATATYPE, ROOT, COMM, IERR
INTEGER BUFFER( * )
IERR = 0
END SUBROUTINE MPI_BCAST

SUBROUTINE MPI_CANCEL( IREQ, IERR )
INTEGER IREQ, IERR
IERR = 0
END SUBROUTINE MPI_CANCEL

SUBROUTINE MPI_COMM_CREATE( COMM, GROUP, COMM2, IERR )

INTEGER COMM, GROUP, COMM2, IERR
IERR = 0
END SUBROUTINE MPI_COMM_CREATE

SUBROUTINE MPI_COMM_DUP( COMM, COMM2, IERR )

INTEGER COMM, COMM2, IERR
IERR = 0
END SUBROUTINE MPI_COMM_DUP

SUBROUTINE MPI_COMM_FREE( COMM, IERR )

INTEGER COMM, IERR
IERR = 0
END SUBROUTINE MPI_COMM_FREE

SUBROUTINE MPI_COMM_GROUP( COMM, GROUP, IERR )

INTEGER COMM, GROUP, IERR
IERR = 0
END SUBROUTINE MPI_COMM_GROUP

SUBROUTINE MPI_COMM_RANK( COMM, RANK, IERR )

INTEGER COMM, RANK, IERR
RANK = 0
IERR = 0
END SUBROUTINE MPI_COMM_RANK

SUBROUTINE MPI_COMM_SIZE( COMM, SIZE, IERR )

INTEGER COMM, SIZE, IERR
SIZE = 1
IERR = 0
END SUBROUTINE MPI_COMM_SIZE

SUBROUTINE MPI_COMM_SPLIT( COMM, COLOR, KEY, COMM2, IERR )

INTEGER COMM, COLOR, KEY, COMM2, IERR
IERR = 0
END SUBROUTINE MPI_COMM_SPLIT

!     SUBROUTINE MPI_ERRHANDLER_SET( COMM, ERRHANDLER, IERR )
!     INTEGER COMM, ERRHANDLER, IERR
!     IERR = 0
!     END SUBROUTINE MPI_ERRHANDLER_SET

SUBROUTINE MPI_FINALIZE( IERR )

INTEGER IERR
IERR = 0
END SUBROUTINE MPI_FINALIZE

SUBROUTINE MPI_GET_COUNT( STATUS, DATATYPE, CNT, IERR )

INTEGER DATATYPE, CNT, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
error stop 'MPI_GET_CNT should not be called.'
END SUBROUTINE MPI_GET_COUNT

SUBROUTINE MPI_GROUP_FREE( GROUP, IERR )

INTEGER GROUP, IERR
IERR = 0
END SUBROUTINE MPI_GROUP_FREE

SUBROUTINE MPI_GROUP_RANGE_EXCL( GROUP, N, RANGES, GROUP2, IERR )

INTEGER GROUP, N, GROUP2, IERR
INTEGER RANGES(*)
IERR = 0
END SUBROUTINE MPI_GROUP_RANGE_EXCL

SUBROUTINE MPI_GROUP_SIZE( GROUP, SIZE, IERR )

INTEGER GROUP, SIZE, IERR
SIZE = 1 ! Or should it be zero ?
IERR = 0
END SUBROUTINE MPI_GROUP_SIZE

SUBROUTINE MPI_INIT(IERR)

INTEGER IERR
IERR = 0
END SUBROUTINE MPI_INIT

SUBROUTINE MPI_INITIALIZED( FLAG, IERR )

LOGICAL FLAG
INTEGER IERR
FLAG = .TRUE.
IERR = 0
END SUBROUTINE MPI_INITIALIZED

SUBROUTINE MPI_IPROBE( SOURCE, TAG, COMM, FLAG, STATUS, IERR )

INTEGER SOURCE, TAG, COMM, IERR
INTEGER STATUS(MPI_STATUS_SIZE)
LOGICAL FLAG
FLAG = .FALSE.
IERR = 0
END SUBROUTINE MPI_IPROBE

SUBROUTINE MPI_IRECV( BUF, CNT, DATATYPE, SOURCE, TAG, COMM, IREQ, IERR )
INTEGER CNT, DATATYPE, SOURCE, TAG, COMM, IREQ, IERR
class(*) :: BUF(..)
IERR = 0
END SUBROUTINE MPI_IRECV

SUBROUTINE MPI_ISEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IREQ, IERR )
INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR, IREQ
class(*) :: BUF(..)
error stop 'MPI_ISEND should not be called.'
END SUBROUTINE MPI_ISEND

SUBROUTINE MPI_TYPE_COMMIT( NEWTYP, IERR_MPI )
INTEGER NEWTYP, IERR_MPI
END SUBROUTINE MPI_TYPE_COMMIT

SUBROUTINE MPI_TYPE_FREE( NEWTYP, IERR_MPI )
INTEGER NEWTYP, IERR_MPI
END SUBROUTINE MPI_TYPE_FREE

SUBROUTINE MPI_TYPE_CONTIGUOUS( LENGTH, DATATYPE, NEWTYPE, IERR_MPI )
INTEGER LENGTH, DATATYPE, NEWTYPE, IERR_MPI
END SUBROUTINE MPI_TYPE_CONTIGUOUS

SUBROUTINE MPI_OP_CREATE( FUNC, COMMUTE, OP, IERR )
EXTERNAL FUNC
LOGICAL COMMUTE
INTEGER OP, IERR
OP = 0
END SUBROUTINE MPI_OP_CREATE

SUBROUTINE MPI_OP_FREE( OP, IERR )
INTEGER OP, IERR
END SUBROUTINE MPI_OP_FREE

SUBROUTINE MPI_PACK( INBUF, INCNT, DATATYPE, OUTBUF, OUTCNT, POSITION, COMM, IERR )
INTEGER INCNT, DATATYPE, OUTCNT, POSITION, COMM, IERR
class(*) :: INBUF(..), OUTBUF(..)
error stop 'MPI_PACKED should not be called.'
END SUBROUTINE MPI_PACK

SUBROUTINE MPI_PACK_SIZE( INCNT, DATATYPE, COMM, SIZE, IERR )
INTEGER INCNT, DATATYPE, COMM, SIZE, IERR
error stop 'MPI_PACK_SIZE should not be called.'
END SUBROUTINE MPI_PACK_SIZE

SUBROUTINE MPI_PROBE( SOURCE, TAG, COMM, STATUS, IERR )
INTEGER SOURCE, TAG, COMM, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
error stop 'MPI_PROBE should not be called.'
END SUBROUTINE MPI_PROBE

SUBROUTINE MPI_RECV( BUF, CNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERR )
INTEGER CNT, DATATYPE, SOURCE, TAG, COMM, IERR
class(*) :: BUF(..)
integer :: STATUS(MPI_STATUS_SIZE)
error stop 'MPI_RECV should not be called.'
END SUBROUTINE MPI_RECV

SUBROUTINE MPI_REQUEST_FREE( IREQ, IERR )
INTEGER IREQ, IERR
IERR = 0
END SUBROUTINE MPI_REQUEST_FREE

SUBROUTINE MPI_SEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IERR )
INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
class(*) :: BUF(..)
error stop 'MPI_SEND should not be called.'
END SUBROUTINE MPI_SEND

SUBROUTINE MPI_SSEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IERR)

INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
class(*) :: BUF(..)
error stop 'MPI_SSEND should not be called.'
END SUBROUTINE MPI_SSEND

SUBROUTINE MPI_TEST( IREQ, FLAG, STATUS, IERR )
INTEGER IREQ, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
LOGICAL FLAG
FLAG = .FALSE.
IERR = 0
END SUBROUTINE MPI_TEST

SUBROUTINE MPI_UNPACK( INBUF, INSIZE, POSITION, OUTBUF, OUTCNT,  DATATYPE, COMM, IERR )
INTEGER INSIZE, POSITION, OUTCNT, DATATYPE, COMM, IERR
class(*) :: INBUF(..), OUTBUF(..)
error stop 'MPI_UNPACK should not be called.'
END SUBROUTINE MPI_UNPACK

SUBROUTINE MPI_WAIT( IREQ, STATUS, IERR )
INTEGER IREQ, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
error stop 'MPI_WAIT should not be called.'
END SUBROUTINE MPI_WAIT

SUBROUTINE MPI_WAITALL( CNT, ARRAY_OF_REQUESTS, STATUS, IERR )
INTEGER CNT, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
INTEGER ARRAY_OF_REQUESTS( CNT )
error stop 'MPI_WAITALL should not be called.'
END SUBROUTINE MPI_WAITALL

SUBROUTINE MPI_WAITANY( CNT, ARRAY_OF_REQUESTS, INDEX, STATUS, IERR )
INTEGER CNT, INDEX, IERR
INTEGER STATUS( MPI_STATUS_SIZE )
INTEGER ARRAY_OF_REQUESTS( CNT )
error stop 'MPI_WAITANY should not be called.'
END SUBROUTINE MPI_WAITANY

real(real64) FUNCTION MPI_WTIME( )
!     elapsed time
real(real64) :: VAL
!     write(*,*) 'Entering MPI_WTIME'
CALL MUMPS_ELAPSE( VAL )
MPI_WTIME = VAL
!     write(*,*) 'Exiting MPI_WTIME'
END FUNCTION MPI_WTIME



!  Utilities to copy data

subroutine mumps_copy(N, S, R)
class(*), dimension(N), intent(in) :: R
class(*), dimension(N), intent(out) :: S
integer, intent(in) :: N

select type (S)
type is (real(real32))
  select type (R)
  type is (real(real32))
    R = S
  end select
type is (real(real64))
  select type (R)
  type is (real(real64))
    R = S
  end select
type is (integer(int32))
  select type (R)
  type is (integer(int32))
    R = S
  end select
end select
end subroutine mumps_copy

SUBROUTINE MUMPS_COPY_INTEGER( S, R, N )
INTEGER N
INTEGER S(N),R(N)
R = S
END SUBROUTINE MUMPS_COPY_INTEGER

SUBROUTINE MUMPS_COPY_INTEGER8( S, R, N )
INTEGER N
INTEGER(8) S(N),R(N)
R = S
END SUBROUTINE MUMPS_COPY_INTEGER8

SUBROUTINE MUMPS_COPY_LOGICAL( S, R, N )
INTEGER N
LOGICAL S(N),R(N)
R = S
END SUBROUTINE MUMPS_COPY_LOGICAL

SUBROUTINE MUMPS_COPY_2INTEGER( S, R, N )
INTEGER N
INTEGER S(N+N),R(N+N)
R = S
END SUBROUTINE MUMPS_COPY_2INTEGER

SUBROUTINE MUMPS_COPY_REAL( S, R, N )
INTEGER N
REAL S(N),R(N)
R = S
END

SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION( S, R, N )
INTEGER N
DOUBLE PRECISION S(N+N),R(N+N)
R = S
END SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION

SUBROUTINE MUMPS_COPY_DOUBLE_PRECISION( S, R, N )
INTEGER N
DOUBLE PRECISION S(N),R(N)
R = S
END

SUBROUTINE MUMPS_COPY_COMPLEX( S, R, N )
INTEGER N
COMPLEX S(N),R(N)
R = S
END SUBROUTINE MUMPS_COPY_COMPLEX

SUBROUTINE MUMPS_COPY_DOUBLE_COMPLEX( S, R, N )
INTEGER N
COMPLEX(kind=kind(0.0D0)) :: S(N),R(N)
R = S
END

LOGICAL FUNCTION MUMPS_IS_IN_PLACE( SENDBUF, CNT )
class(*) :: SENDBUF(:)
INTEGER, intent(in) :: CNT
INTEGER :: I

!! Check address using C code
! MUMPS_IS_IN_PLACE = .FALSE.
! IF ( CNT .GT. 0 ) THEN
!   CALL MUMPS_CHECKADDREQUAL(SENDBUF(1), MPI_IN_PLACE, I)
!   MUMPS_IS_IN_PLACE = I == 1
! ENDIF

! Begin old code which requires the MPI_IN_PLACE
! variable to have the F2003 attribute VOLATILE
IF ( CNT .GT. 0 ) THEN
  MPI_IN_PLACE = -1

  select type (sendbuf)
  type is (real(real32))
    if (SENDBUF(1) == MPI_IN_PLACE) then
      MPI_IN_PLACE = -9876543
      MUMPS_IS_IN_PLACE = SENDBUF(1) == MPI_IN_PLACE
    endif
  type is (real(real64))
    if (SENDBUF(1) == MPI_IN_PLACE) then
      MPI_IN_PLACE = -9876543
      MUMPS_IS_IN_PLACE = SENDBUF(1) == MPI_IN_PLACE
    endif
  type is (integer(int32))
    if (SENDBUF(1) == MPI_IN_PLACE) then
      MPI_IN_PLACE = -9876543
      MUMPS_IS_IN_PLACE = SENDBUF(1) == MPI_IN_PLACE
    endif
  class default
    error stop 'MUMPS_IS_IN_PLACE: unknown type'
  end select

endif

! End old code
END FUNCTION MUMPS_IS_IN_PLACE

! Begin old code
! LOGICAL FUNCTION MUMPS_CHECK_EQUAL(I,J)
! INTEGER :: I,J
! IF (I.EQ.J) THEN
!   MUMPS_CHECK_EQUAL = .TRUE.
! ELSE
!   MUMPS_CHECK_EQUAL = .FALSE.
! ENDIF
! END FUNCTION MUMPS_CHECK_EQUAL
! End old code


end module mpi
