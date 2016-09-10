********************************************************************************
** FICHE F.14.  ALGORITHM TO HANDLE INDICES IN CONSTANT MU VT MONTE CARLO     **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** INDEX-HANDLING IN GRAND CANONICAL MONTE CARLO SIMULATION.     **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
C    **    ADDS AN ATOM TO THE ARRAY LOCATE.                          **
C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
C    **    REMOVES AN ATOM FROM THE ARRAY LOCATE.                     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS BEFORE TRIAL      **
C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS           **
C    ** INTEGER IPULL               INDEX OF ATOM FOR REMOVAL         **
C    ** INTEGER NLOC                POSITION OF N IN LOCATE           **
C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING TRIAL      **
C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF AN ATOM  **
C    ** REAL    RX(NMAX), ETC.      POSITIONS OF ATOMS                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** ROUTINE ADD IS CALLED AFTER A SUCCESSFUL TRIAL ADDITION.      **
C    ** ROUTINE REMOVE IS CALLED AFTER A SUCCESSFUL TRIAL REMOVAL.    **
C    ** THE ARRAY LOCATE IS UPDATED IN EACH CASE.                     **
C    *******************************************************************



        SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** SUBROUTINE TO ADD AN ATOM TO THE ARRAY LOCATE.                **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE NEW ADDITION   **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )

        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        INTEGER     LOCATE(NMAX)

        REAL        RXNEW, RYNEW, RZNEW
        INTEGER     N, IPULL

        INTEGER     INEW, NTRIAL

C    *******************************************************************

        NTRIAL = N + 1
        INEW = LOCATE(NTRIAL)

        IF ( INEW .EQ. 0 ) THEN

C       ** ATOM REQUIRES A NEW NUMBER **

           LOCATE(NTRIAL) = NTRIAL
           INEW           = NTRIAL

        ENDIF

C    ** FIT NEW ATOM INTO THE ARRAY **

        RX(INEW) = RXNEW
        RY(INEW) = RYNEW
        RZ(INEW) = RZNEW

        RETURN
        END



        SUBROUTINE REMOVE ( NLOC, IPULL, N )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** SUBROUTINE TO REMOVE AN ATOM FROM THE ARRAY LOCATE.           **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE REMOVAL.       **
C    ** ELEMENT IPULL OF LOCATE IS TO BE DESTROYED.                   **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )

        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        INTEGER     LOCATE(NMAX)
        INTEGER     N, IPULL, NLOC

        INTEGER     K

C    *******************************************************************

        IF ( NLOC .LT. N ) THEN

C       ** CLOSE UP THE ARRAY LOCATE AFTER THE REMOVAL **

           DO 10 K = NLOC + 1, N

              LOCATE(K - 1) = LOCATE(K)

10         CONTINUE

C       ** PLACE THE GHOST ATOM IPULL JUST OUTSIDE THE ACTIVE **
C       ** RANGE OF THE ARRAY LOCATE FOR FUTURE USE           **

           LOCATE(N) = IPULL

        ENDIF

        RETURN
        END



