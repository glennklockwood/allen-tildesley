********************************************************************************
** FICHE F.9.  RATTLE ALGORITHM FOR CONSTRAINT DYNAMICS OF A CHAIN MOLECULE   **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** CONSTRAINT DYNAMICS OF A CHAIN OF ATOMS USING RATTLE.         **
C    **                                                               **
C    ** WE APPLY BOND LENGTH CONSTRAINTS TO ADJACENT ATOMS ONLY IN A  **
C    ** CHAIN MOLECULE WHICH MAY BE CYCLIC.  THE GENERALIZATION TO    **
C    ** MORE COMPLICATED SYSTEMS IS STRAIGHTFORWARD.  THE CONSTRAINT  **
C    ** EQUATIONS ARE LINEARIZED, AND EACH CONSTRAINT IS TREATED IN   **
C    ** TURN, UNTIL BOND LENGTHS ARE SATISFIED TO WITHIN A SPECIFIED  **
C    ** TOLERANCE.                                                    **
C    ** IN THIS EXAMPLE WE TAKE A 6-ATOM CHAIN.                       **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** HC ANDERSEN, J. COMPUT. PHYS. 52, 24, 1983.                   **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVEA ( DT, TOL, MAXIT, NB, BOX )                  **
C    **    ADVANCES POSITIONS AND HALF ADVANCES VELOCITIES WITH       **
C    **    APPLIED CONSTRAINTS                                        **
C    ** SUBROUTINE MOVEB ( DT, TOL, MAXIT, NB, BOX, K, WC )           **
C    **    COMPLETES VELOCITY MOVE AND CALCULATES NEW KINETIC ENERGY  **
C    **    AND CONSTRAINT CONTRIBUTION TO VIRIAL.                     **
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
C    ** REAL    VX(N,NA),VY(N,NA),VZ(N,NA) ATOM VELOCITIES            **
C    ** REAL    FX(N,NA),FY(N,NA),FZ(N,NA) ATOM FORCES                **
C    ** REAL    DSQ(NA)                    SQUARED BOND LENGTHS       **
C    ** REAL    M(NA)                      ATOMIC MASSES              **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THESE ROUTINES COMPUTE CONSTRAINT EFFECTS IN AN ITERATIVE WAY.**
C    ** POSITIONS, VELOCITIES, AND FORCES AT TIME T ARE SUPPLIED TO   **
C    ** THE FIRST ROUTINE MOVEA.                                      **
C    ** THE VELOCITY VERLET ALGORITHM IS USED TO ADVANCE THE          **
C    ** POSITIONS THROUGH A TIMESTEP T -> T+DT FROM RX,RY,RZ TO       **
C    ** PX,PY,PZ, AND THE VELOCITIES VX,VY,VZ THROUGH HALF A TIMESTEP **
C    ** T -> T+DT/2, WITHOUT ANY CONSTRAINTS APPLIED:                 **
C    ** PX(T+DT)   = RX(T) + VX(T)*DT + AX(T)*DT**2/2 ETC.            **
C    ** VX(T+DT/2) = VX(T) + AX(T)*DT/2 ETC.                          **
C    ** THE DESIRED SQUARED BOND LENGTHS AND ATOMIC MASSES ARE THEN   **
C    ** USED TO APPLY CONSTRAINTS TO POSITIONS AND HALF-STEP          **
C    ** VELOCITIES.                                                   **
C    ** DSQ(A) CONTAINS SQUARED BOND LENGTH BETWEEN ATOMS A AND A+1.  **
C    ** IF NB=NA THE MOLECULE IS CYCLIC, IF NB=NA-1 IT IS NOT.        **
C    ** THE ROUTINE ALSO REQUIRES THE DESIRED TOLERANCE AND AN UPPER  **
C    ** LIMIT TO THE NUMBER OF ITERATIONS IN CASE OF NON-CONVERGENCE. **
C    ** THE ROUTINE USES TWO LOGICAL ARRAYS TO KEEP TRACK OF WHETHER  **
C    ** OR NOT WE HAVE MOVED (I.E. CORRECTED) THE ATOM POSITIONS:     **
C    ** MOVING(A) A=1,NA SAYS WHETHER WE ARE MOVING ATOM A THIS TIME  **
C    ** MOVED(A)  A=1,NA SAYS WHETHER WE MOVED ATOM A LAST TIME.      **
C    ** THIS IS SO THAT WE CAN STOP CORRECTING THE POSITIONS OF ATOMS **
C    ** WHENEVER POSSIBLE, SO AS TO CUT DOWN ON UNNECESSARY WORK.     **
C    ** THE ROUTINE RETURNS FINAL VALUES IN RX,RY,RZ,VX,VY,VZ.        **
C    ** NEW FORCES ARE COMPUTED FROM THE POSITIONS IN A FORCE ROUTINE **
C    ** (NOT SUPPLIED HERE) AND THE SECOND ROUTINE MOVEB CALLED.      **
C    ** THIS ADVANCES THE VELOCITIES FROM T+DT/2 TO T+DT:             **
C    ** VX(T+DT) = VX(T+DT/2) + AX(T+DT)*DT/2 ETC.                    **
C    ** AND COMPLETES THE CONSTRAINT PROCEDURE ON VX,VY,VZ.           **
C    ** IT ALSO COMPUTES KINETIC ENERGY AND CONSTRAINT VIRIAL.        **
C    *******************************************************************

        SUBROUTINE MOVEA ( DT, TOL, MAXIT, NB, BOX )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / DSQ, M

C    *******************************************************************
C    ** FIRST PART OF VELOCITY VERLET ALGORITHM WITH CONSTRAINTS      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     NA
        PARAMETER ( NA = 6 )

        REAL        DT, TOL, BOX
        INTEGER     MAXIT, NB
        REAL        RX(N,NA), RY(N,NA), RZ(N,NA)
        REAL        VX(N,NA), VY(N,NA), VZ(N,NA)
        REAL        FX(N,NA), FY(N,NA), FZ(N,NA)
        REAL        DSQ(NA), M(NA)

        LOGICAL     DONE
        LOGICAL     MOVING(NA), MOVED(NA)
        REAL        RXI(NA), RYI(NA), RZI(NA)
        REAL        PXI(NA), PYI(NA), PZI(NA)
        REAL        VXI(NA), VYI(NA), VZI(NA)
        REAL        TOL2, PXAB, PYAB, PZAB, PABSQ, DT2, DTSQ2
        REAL        RABSQ, DIFFSQ, RXAB, RYAB, RZAB, RPAB, GAB
        REAL        DX, DY, DZ, RMA, RMB, BOXINV, RPTOL
        REAL        AXIA, AYIA, AZIA
        INTEGER     I, A, B, IT
        PARAMETER ( RPTOL = 1.0E-6 )

C    *******************************************************************

        IF ( ( NB .NE. NA ) .AND. ( NB .NE. NA-1 ) ) STOP 'NB IN ERROR'

        BOXINV = 1.0 / BOX
        TOL2   = 2.0 * TOL
        DT2    = DT / 2.0
        DTSQ2  = DT * DT2

C    ** LOOP OVER MOLECULES **

        DO 2000 I = 1, N

C       ** VELOCITY VERLET ALGORITHM PART A **

           DO 100 A = 1, NA

              AXIA = FX(I,A) / M(A)
              AYIA = FY(I,A) / M(A)
              AZIA = FZ(I,A) / M(A)

              RXI(A) = RX(I,A)
              RYI(A) = RY(I,A)
              RZI(A) = RZ(I,A)
              PXI(A) = RX(I,A) + DT * VX(I,A) + DTSQ2 * AXIA
              PYI(A) = RY(I,A) + DT * VY(I,A) + DTSQ2 * AYIA
              PZI(A) = RZ(I,A) + DT * VZ(I,A) + DTSQ2 * AZIA
              VXI(A) = VX(I,A) + DT2 * AXIA
              VYI(A) = VY(I,A) + DT2 * AYIA
              VZI(A) = VZ(I,A) + DT2 * AZIA

              MOVING(A) = .FALSE.
              MOVED(A)  = .TRUE.

100        CONTINUE

           IT = 0
           DONE = .FALSE.

C       ** START OF ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              DONE = .TRUE.

              DO 300 A = 1, NB

                 B = A + 1
                 IF ( B .GT. NA ) B = 1

                 IF ( MOVED(A) .OR. MOVED(B) ) THEN

                    PXAB = PXI(A) - PXI(B)
                    PYAB = PYI(A) - PYI(B)
                    PZAB = PZI(A) - PZI(B)
                    PXAB = PXAB - ANINT ( PXAB * BOXINV ) * BOX
                    PYAB = PYAB - ANINT ( PYAB * BOXINV ) * BOX
                    PZAB = PZAB - ANINT ( PZAB * BOXINV ) * BOX

                    PABSQ = PXAB ** 2 + PYAB ** 2 + PZAB ** 2
                    RABSQ = DSQ(A)
                    DIFFSQ = RABSQ - PABSQ

                    IF ( ABS ( DIFFSQ ) .GT. ( RABSQ * TOL2 ) ) THEN

                       RXAB = RXI(A) - RXI(B)
                       RYAB = RYI(A) - RYI(B)
                       RZAB = RZI(A) - RZI(B)
                       RXAB = RXAB - ANINT ( RXAB * BOXINV ) * BOX
                       RYAB = RYAB - ANINT ( RYAB * BOXINV ) * BOX
                       RZAB = RZAB - ANINT ( RZAB * BOXINV ) * BOX
                       RPAB = RXAB * PXAB + RYAB * PYAB + RZAB * PZAB

                       IF ( RPAB .LT. ( RABSQ * RPTOL ) ) THEN

                          WRITE(*,'('' CONSTRAINT FAILURE '')')
                          STOP

                       ENDIF

                       RMA = 1.0 / M(A)
                       RMB = 1.0 / M(B)
                       GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB )
                       DX  = RXAB * GAB
                       DY  = RYAB * GAB
                       DZ  = RZAB * GAB

                       PXI(A) = PXI(A) + RMA * DX
                       PYI(A) = PYI(A) + RMA * DY
                       PZI(A) = PZI(A) + RMA * DZ
                       PXI(B) = PXI(B) - RMB * DX
                       PYI(B) = PYI(B) - RMB * DY
                       PZI(B) = PZI(B) - RMB * DZ

                       DX = DX / DT
                       DY = DY / DT
                       DZ = DZ / DT

                       VXI(A) = VXI(A) + RMA * DX
                       VYI(A) = VYI(A) + RMA * DY
                       VZI(A) = VZI(A) + RMA * DZ
                       VXI(B) = VXI(B) - RMB * DX
                       VYI(B) = VYI(B) - RMB * DY
                       VZI(B) = VZI(B) - RMB * DZ

                       MOVING(A) = .TRUE.
                       MOVING(B) = .TRUE.
                       DONE = .FALSE.

                    ENDIF

                 ENDIF

300           CONTINUE

              DO 500 A = 1, NA

                 MOVED(A)  = MOVING(A)
                 MOVING(A) = .FALSE.

500           CONTINUE

              IT = IT + 1
              GOTO 1000

           ENDIF

C       ** END OF ITERATIVE LOOP **

           IF (.NOT. DONE) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS IN MOVEA '')')
              WRITE(*,'('' MOLECULE '',I5)') I
              STOP

           ENDIF

C       ** STORE AWAY NEW VALUES    **

           DO 600 A = 1, NA

              RX(I,A) = PXI(A)
              RY(I,A) = PYI(A)
              RZ(I,A) = PZI(A)
              VX(I,A) = VXI(A)
              VY(I,A) = VYI(A)
              VZ(I,A) = VZI(A)

600        CONTINUE

2000    CONTINUE

C    ** END OF LOOP OVER MOLECULES **

        RETURN
        END



        SUBROUTINE MOVEB ( DT, TOL, MAXIT, NB, BOX, K, WC )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / DSQ, M

C    *******************************************************************
C    ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     NA
        PARAMETER ( NA = 6 )

        REAL        DT, TOL, BOX, K, WC
        INTEGER     MAXIT, NB
        REAL        RX(N,NA), RY(N,NA), RZ(N,NA)
        REAL        VX(N,NA), VY(N,NA), VZ(N,NA)
        REAL        FX(N,NA), FY(N,NA), FZ(N,NA)
        REAL        DSQ(NA), M(NA)

        LOGICAL     DONE
        LOGICAL     MOVING(NA), MOVED(NA)
        REAL        RXI(NA), RYI(NA), RZI(NA)
        REAL        VXI(NA), VYI(NA), VZI(NA)
        REAL        RXAB, RYAB, RZAB, RVAB, GAB
        REAL        VXAB, VYAB, VZAB
        REAL        DX, DY, DZ, DT2, RMA, RMB, BOXINV
        INTEGER     I, A, B, IT

C    *******************************************************************

        BOXINV = 1.0 / BOX
        DT2    = DT / 2.0
        K      = 0.0
        WC     = 0.0

C    ** LOOP OVER ALL MOLECULES **

        DO 2000 I = 1, N

C       ** VELOCITY VERLET ALGORITHM PART B **

           DO 100 A = 1, NA

              RXI(A) = RX(I,A)
              RYI(A) = RY(I,A)
              RZI(A) = RZ(I,A)
              VXI(A) = VX(I,A) + DT2 * FX(I,A) / M(A)
              VYI(A) = VY(I,A) + DT2 * FY(I,A) / M(A)
              VZI(A) = VZ(I,A) + DT2 * FZ(I,A) / M(A)

              MOVING(A) = .FALSE.
              MOVED(A)  = .TRUE.

100        CONTINUE

C       ** START OF ITERATIVE LOOP **

           IT = 0
           DONE = .FALSE.

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              DONE = .TRUE.

              DO 300 A = 1, NB

                 B = A + 1
                 IF ( B .GT. NA ) B = 1

                 IF ( MOVED(A) .OR. MOVED(B) ) THEN

                    VXAB = VXI(A) - VXI(B)
                    VYAB = VYI(A) - VYI(B)
                    VZAB = VZI(A) - VZI(B)
                    RXAB = RXI(A) - RXI(B)
                    RYAB = RYI(A) - RYI(B)
                    RZAB = RZI(A) - RZI(B)
                    RXAB = RXAB - ANINT ( RXAB * BOXINV ) * BOX
                    RYAB = RYAB - ANINT ( RYAB * BOXINV ) * BOX
                    RZAB = RZAB - ANINT ( RZAB * BOXINV ) * BOX
                    RVAB = RXAB * VXAB + RYAB * VYAB + RZAB * VZAB
                    RMA  = 1.0 / M(A)
                    RMB  = 1.0 / M(B)
                    GAB  = -RVAB / ( ( RMA + RMB ) * DSQ(A) )

                    IF ( ABS ( GAB ) .GT. TOL ) THEN

                       WC = WC + GAB * DSQ(A)
                       DX = RXAB * GAB
                       DY = RYAB * GAB
                       DZ = RZAB * GAB

                       VXI(A) = VXI(A) + RMA * DX
                       VYI(A) = VYI(A) + RMA * DY
                       VZI(A) = VZI(A) + RMA * DZ
                       VXI(B) = VXI(B) - RMB * DX
                       VYI(B) = VYI(B) - RMB * DY
                       VZI(B) = VZI(B) - RMB * DZ

                       MOVING(A) = .TRUE.
                       MOVING(B) = .TRUE.
                       DONE = .FALSE.

                    ENDIF

                 ENDIF

300           CONTINUE

              DO 500 A = 1, NA

                 MOVED(A)  = MOVING(A)
                 MOVING(A) = .FALSE.

500           CONTINUE

              IT = IT + 1
              GOTO 1000

           ENDIF

C       ** END OF ITERATIVE LOOP **

           IF (.NOT. DONE) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS IN MOVEB '')')
              WRITE(*,'('' MOLECULE '',I5)') I
              STOP

           ENDIF

           DO 600 A = 1, NA

              VX(I,A) = VXI(A)
              VY(I,A) = VYI(A)
              VZ(I,A) = VZI(A)
              K = K + M(A) * ( VXI(A) ** 2 + VYI(A) ** 2 + VZI(A) ** 2 )

600        CONTINUE

2000    CONTINUE

C    ** END OF LOOP OVER MOLECULES **

        K  = K * 0.5
        WC = WC / DT2 / 3.0

        RETURN
        END



