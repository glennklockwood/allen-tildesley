********************************************************************************
** FICHE F.26.  ROUTINES TO FOLD/UNFOLD TRAJECTORIES IN PERIODIC BOUNDARIES   **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        SUBROUTINE FOLD ( N, RX, RY, RZ )

C    *******************************************************************
C    ** SUBROUTINE TO FOLD TRAJECTORIES IN PERIODIC BOUNDARIES.       **
C    **                                                               **
C    ** THE FOLDING ROUTINE IS SIMPLY THE USUAL APPLICATION OF        **
C    ** BOUNDARY CONDITIONS.  WE TAKE THE UNIT CUBE AS AN EXAMPLE.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE IS CALLED FOR EVERY CONFIGURATION GENERATED IN A  **
C    ** SIMULATION.                                                   **
C    *******************************************************************

        INTEGER N
        REAL    RX(N), RY(N), RZ(N)

        INTEGER I

C    *******************************************************************

        DO 100 I = 1, N

           RX(I) = RX(I) - ANINT ( RX(I) )
           RY(I) = RY(I) - ANINT ( RY(I) )
           RZ(I) = RZ(I) - ANINT ( RZ(I) )

100     CONTINUE

        RETURN
        END



        SUBROUTINE UNFOLD ( N, RX, RY, RZ, RX0, RY0, RZ0 )

C    *******************************************************************
C    ** SUBROUTINE TO FOLD TRAJECTORIES IN PERIODIC BOUNDARIES.       **
C    **                                                               **
C    ** THE UNFOLDING ROUTINE UNDOES THE EFFECT OF FOLDING.           **
C    ** AGAIN WE TAKE THE UNIT CUBE AS AN EXAMPLE.                    **
C    ** THIS ROUTINE REQUIRES THAT COORDINATES FROM SUCCESSIVE STEPS  **
C    ** BE SUPPLIED AND ASSUMES THAT NATURAL MOVEMENT ACROSS HALF A   **
C    ** BOX LENGTH IN ONE STEP WILL NEVER OCCUR.                      **
C    ** THE RESULTING COORDINATES WOULD BE SUITABLE, FOR EXAMPLE, TO  **
C    ** CALCULATE THE DIFFUSION COEFFICIENT BY EINSTEIN'S RELATION.   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                    NUMBER OF MOLECULES              **
C    ** REAL    RX(N),RY(N),RZ(N)    MOLECULAR POSITIONS AT TIME T    **
C    ** REAL    RX0(N),RY0(N),RY0(N) MOLECULAR POSITIONS AT TIME T-DT **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** (1)  READ INITIAL COORDINATES INTO RX0,RY0,RZ0                **
C    ** (2)  WRITE RX0,RY0,RZ0 TO OUTPUT FILE (OR USE IMMEDIATELY)    **
C    ** (3)  READ NEXT STEP COORDINATES INTO RX,RY,RZ                 **
C    ** (4)  CALL UNFOLD ( N, RX, RY, RZ, RX0, RY0, RZ0 )             **
C    ** (5)  WRITE RX,RY,RZ TO OUTPUT FILE (OR USE IMMEDIATELY)       **
C    ** (6)  SET RX0(I)=RX(I), RY0(I)=RY(I), RZ0(I)=RZ(I), I = 1,N    **
C    ** (7)  UNLESS DATA EXHAUSTED, GO TO (3)                         **
C    *******************************************************************

        INTEGER N
        REAL    RX(N), RY(N), RZ(N), RX0(N), RY0(N), RZ0(N)

        INTEGER I
        REAL    DX, DY, DZ

C    *******************************************************************

        DO 100 I = 1, N

           DX = RX(I) - RX0(I)
           DY = RY(I) - RY0(I)
           DZ = RZ(I) - RZ0(I)
           DX = DX - ANINT ( DX )
           DY = DY - ANINT ( DY )
           DZ = DZ - ANINT ( DZ )
           RX(I) = RX0(I) + DX
           RY(I) = RY0(I) + DY
           RZ(I) = RZ0(I) + DZ

100     CONTINUE

        RETURN
        END



