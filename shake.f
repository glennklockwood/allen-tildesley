********************************************************************************
** FICHE F.8.  SHAKE ALGORITHM FOR CONSTRAINT DYNAMICS OF A CHAIN MOLECULE    **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        SUBROUTINE MOVE ( DT, TOL, MAXIT, NB, BOX, K, WC )

        COMMON / BLOCK1 / RX, RY, RZ, OX, OY, OZ, FX, FY, FZ
        COMMON / BLOCK2 / DSQ, M

C    *******************************************************************
C    ** CONSTRAINT DYNAMICS OF A CHAIN OF ATOMS USING SHAKE.          **
C    **                                                               **
C    ** THIS ROUTINE COMPUTES CONSTRAINT EFFECTS IN AN ITERATIVE WAY. **
C    ** WE APPLY BOND LENGTH CONSTRAINTS TO ADJACENT ATOMS ONLY IN A  **
C    ** CHAIN MOLECULE WHICH MAY BE CYCLIC.  THE GENERALIZATION TO    **
C    ** TO MORE COMPLICATED SYSTEMS IS STRAIGHTFORWARD.  THE          **
C    ** CONSTRAINT EQUATIONS ARE LINEARIZED, AND EACH CONSTRAINT IS   **
C    ** TREATED IN TURN, UNTIL BOND LENGTHS ARE SATISFIED TO WITHIN A **
C    ** SPECIFIED TOLERANCE.                                          **
C    ** IN THIS EXAMPLE WE TAKE A 6-ATOM CHAIN.                       **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** RYCKAERT ET AL., J. COMP. PHYS. 23, 327, 1977.                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                          NUMBER OF MOLECULES        **
C    ** INTEGER NA                         NUMBER OF ATOMS PER MOL.   **
C    ** REAL    DT                         TIMESTEP                   **
C    ** REAL    TOL                        BOND LENGTH TOLERANCE      **
C    ** INTEGER MAXIT                      MAXIMUM ALLOWED ITERATIONS **
C    ** INTEGER NB                         NUMBER OF BONDS            **
C    ** REAL    BOX                        BOX LENGTH                 **
C    ** REAL    K                          KINETIC ENERGY             **
C    ** REAL    WC                         CONSTRAINT VIRIAL          **
C    ** REAL    RX(N,NA),RY(N,NA),RZ(N,NA) ATOM POSITIONS AT TIME T   **
C    ** REAL    OX(N,NA),OY(N,NA),OZ(N,NA) OLD POSITIONS AT TIME T-DT **
C    ** REAL    FX(N,NA),FY(N,NA),FZ(N,NA) ATOM FORCES AT TIME T      **
C    ** REAL    DSQ(NA)                    SQUARED BOND LENGTHS       **
C    ** REAL    M(NA)                      ATOMIC MASSES              **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE IS CALLED AFTER COMPUTATION OF THE FORCES.        **
C    ** THE VERLET ALGORITHM IS USED TO ADVANCE THE POSITIONS FROM    **
C    ** RX,RY,RZ TO PX,PY,PZ, WITHOUT ANY CONSTRAINTS APPLIED.        **
C    ** THE ROUTINE THEN USES THE DESIRED SQUARED BOND LENGTHS IN DSQ **
C    ** TO CONSTRAIN THE NEW POSITIONS.                               **
C    ** DSQ(A) IS THE SQUARED BOND LENGTH BETWEEN ATOM A AND A+1.     **
C    ** NB IS THE NUMBER OF SUCH BONDS.  IF NB = NA - 1 THE MOLECULE  **
C    ** IS NON-CYCLIC WHILE IF NB = NA THE MOLECULE IS CYCLIC.        **
C    ** THE ROUTINE ALSO REQUIRES THE DESIRED TOLERANCE AND AN UPPER  **
C    ** LIMIT TO THE NUMBER OF ITERATIONS IN CASE OF NON-CONVERGENCE. **
C    ** THE ROUTINE USES TWO LOGICAL ARRAYS TO KEEP TRACK OF WHETHER  **
C    ** OR NOT WE HAVE MOVED (I.E. CORRECTED) THE ATOM POSITIONS:     **
C    ** MOVING(A) A=1,NA SAYS WHETHER WE ARE MOVING ATOM A THIS TIME  **
C    ** MOVED(A)  A=1,NA SAYS WHETHER WE MOVED ATOM A LAST TIME.      **
C    ** THIS IS SO THAT WE CAN STOP CORRECTING THE POSITIONS OF ATOMS **
C    ** WHENEVER POSSIBLE, SO AS TO CUT DOWN ON UNNECESSARY WORK.     **
C    ** THE ROUTINE RETURNS WITH NEW (T+DT) CONSTRAINED POSITIONS IN  **
C    ** RX,RY,RZ AND OLD (T) POSITIONS IN OX,OY,OZ.                   **
C    ** THE ROUTINE ALSO CALCULATES THE KINETIC ENERGY AT TIME T      **
C    ** AND THE CONSTRAINT CONTRIBUTION TO THE VIRIAL WC.             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     NA
        PARAMETER ( NA = 6 )

        INTEGER     MAXIT, NB
        REAL        RX(N,NA), RY(N,NA), RZ(N,NA)
        REAL        OX(N,NA), OY(N,NA), OZ(N,NA)
        REAL        FX(N,NA), FY(N,NA), FZ(N,NA)
        REAL        DSQ(NA)
        REAL        M(NA)
        REAL        TOL, DT, K, WC, BOX

        REAL        RXI(NA), RYI(NA), RZI(NA)
        REAL        PXI(NA), PYI(NA), PZI(NA)
        LOGICAL     MOVING(NA)
        LOGICAL     MOVED(NA)

        LOGICAL     DONE
        INTEGER     IT, A, B, I
        REAL        PXAB, PYAB, PZAB, PABSQ
        REAL        RXAB, RYAB, RZAB, RABSQ, DIFFSQ, RPAB
        REAL        GAB, DX, DY, DZ, TOL2, DTSQ, VXIA, VYIA, VZIA
        REAL        BOXINV, RPTOL, RMA, RMB
        PARAMETER ( RPTOL = 1.0E-6 )

C    *******************************************************************

        IF ( ( NB .NE. NA ) .AND. ( NB .NE. NA-1 ) ) STOP 'NB IN ERROR'

        BOXINV = 1.0 / BOX
        DTSQ   = DT ** 2
        TOL2   = 2.0 * TOL
        K      = 0.0
        WC     = 0.0

C    ** LOOP OVER ALL MOLECULES **

        DO 2000 I = 1, N

C       ** VERLET ALGORITHM **

           DO 100 A = 1, NA

              RXI(A) = RX(I,A)
              RYI(A) = RY(I,A)
              RZI(A) = RZ(I,A)
              PXI(A) = 2.0 * RX(I,A) - OX(I,A) + DTSQ * FX(I,A) / M(A)
              PYI(A) = 2.0 * RY(I,A) - OY(I,A) + DTSQ * FY(I,A) / M(A)
              PZI(A) = 2.0 * RZ(I,A) - OZ(I,A) + DTSQ * FZ(I,A) / M(A)

              MOVING(A) = .FALSE.
              MOVED(A)  = .TRUE.

100        CONTINUE

           IT = 0
           DONE = .FALSE.

C       ** BEGIN ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              DONE = .TRUE.

              DO 300 A = 1, NB

                 B = A + 1
                 IF ( B .GT. NA ) B = 1

                 IF ( MOVED(A) .OR. MOVED(B) ) THEN

                    PXAB = PXI(A) - PXI(B)
                    PXAB = PXAB - ANINT ( PXAB * BOXINV ) * BOX
                    PYAB = PYI(A) - PYI(B)
                    PYAB = PYAB - ANINT ( PYAB * BOXINV ) * BOX
                    PZAB = PZI(A) - PZI(B)
                    PZAB = PZAB - ANINT ( PZAB * BOXINV ) * BOX

                    PABSQ  = PXAB ** 2 + PYAB ** 2 + PZAB ** 2
                    RABSQ  = DSQ(A)
                    DIFFSQ = RABSQ - PABSQ

                    IF ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) THEN

                       RXAB = RXI(A) - RXI(B)
                       RXAB = RXAB - ANINT ( RXAB * BOXINV ) * BOX
                       RYAB = RYI(A) - RYI(B)
                       RYAB = RYAB - ANINT ( RYAB * BOXINV ) * BOX
                       RZAB = RZI(A) - RZI(B)
                       RZAB = RZAB - ANINT ( RZAB * BOXINV ) * BOX
                       RPAB = RXAB * PXAB + RYAB * PYAB + RZAB * PZAB

                       IF ( RPAB .LT. ( RABSQ * RPTOL ) ) THEN

                          STOP 'CONSTRAINT FAILURE '

                       ENDIF

                       RMA = 1.0 / M(A)
                       RMB = 1.0 / M(B)
                       GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB )
                       WC  = WC + GAB * RABSQ
                       DX  = RXAB * GAB
                       DY  = RYAB * GAB
                       DZ  = RZAB * GAB

                       PXI(A) = PXI(A) + RMA * DX
                       PYI(A) = PYI(A) + RMA * DY
                       PZI(A) = PZI(A) + RMA * DZ
                       PXI(B) = PXI(B) - RMB * DX
                       PYI(B) = PYI(B) - RMB * DY
                       PZI(B) = PZI(B) - RMB * DZ

                       MOVING(A) = .TRUE.
                       MOVING(B) = .TRUE.
                       DONE = .FALSE.

                    ENDIF

                 ENDIF

300           CONTINUE

              DO 400 A = 1, NA

                 MOVED(A) = MOVING(A)
                 MOVING(A) = .FALSE.

400           CONTINUE

              IT = IT + 1
              GOTO 1000

           ENDIF

C       ** END ITERATIVE LOOP **

           IF ( .NOT. DONE ) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
              WRITE(*,'('' MOLECULE '',I5)') I
              STOP

           ENDIF

           DO 500 A = 1, NA

              VXIA = 0.5 * ( PXI(A) - OX(I,A) ) / DT
              VYIA = 0.5 * ( PYI(A) - OY(I,A) ) / DT
              VZIA = 0.5 * ( PZI(A) - OZ(I,A) ) / DT
              K  = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(A)

              RX(I,A) = PXI(A)
              RY(I,A) = PYI(A)
              RZ(I,A) = PZI(A)
              OX(I,A) = RXI(A)
              OY(I,A) = RYI(A)
              OZ(I,A) = RZI(A)

500        CONTINUE

2000    CONTINUE

C    ** END LOOP OVER MOLECULES **

        WC = WC / DTSQ / 3.0
        K  = 0.5 * K

        RETURN
        END



