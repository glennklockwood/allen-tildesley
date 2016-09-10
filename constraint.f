********************************************************************************
** FICHE F.7.  CONSTRAINT DYNAMICS FOR A NONLINEAR TRIATOMIC MOLECULE         **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO SEPARATE PARTS: NONRIGID AND RIGID TRIATOMIC MOLECULES.   **
C    *******************************************************************



C    *******************************************************************
C    ** FICHE F.7 - PART A                                            **
C    ** CONSTRAINT DYNAMICS FOR A NONRIGID TRIATOMIC MOLECULE.        **
C    *******************************************************************



        SUBROUTINE MOVE ( DT, TOL, MAXIT, BOX, K, WC )

        COMMON / BLOCK1 / RX, RY, RZ, OX, OY, OZ, FX, FY, FZ
        COMMON / BLOCK2 / DSQ, M

C    *******************************************************************
C    ** UPDATES ATOMIC POSITIONS WITH BOND CONSTRAINTS APPLIED.       **
C    **                                                               **
C    ** THIS ROUTINE USES THE VERLET ALGORITHM AND APPLIES BOND       **
C    ** LENGTH CONSTRAINTS TO TWO BOND LENGTHS ONLY, 1-2 AND 2-3.     **
C    ** WE SOLVE A SYSTEM OF QUADRATIC EQUATIONS ITERATIVELY, USING   **
C    ** MATRIX INVERSION TO SOLVE ASSOCIATED LINEAR EQUATIONS,        **
C    ** AND ITERATING TO CONVERGENCE.                                 **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** RYCKAERT ET AL., J. COMP. PHYS, 23, 327, 1977.                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                          NUMBER OF MOLECULES        **
C    ** INTEGER NA                         ATOMS PER MOL. (3 HERE)    **
C    ** REAL    DT                         TIMESTEP                   **
C    ** INTEGER MAXIT                      MAX NUMBER OF ITERATIONS   **
C    ** REAL    TOL                        PRESCRIBED TOLERANCE       **
C    ** REAL    BOX                        BOX LENGTH                 **
C    ** REAL    K                          KINETIC ENERGY             **
C    ** REAL    WC                         CONSTRAINT VIRIAL          **
C    ** REAL    RX(N,NA),RY(N,NA),RZ(N,NA) ATOMIC POSITIONS AT TIME T **
C    ** REAL    OX(N,NA),OY(N,NA),OZ(N,NA) OLD POSITIONS AT TIME T-DT **
C    ** REAL    FX(N,NA),FY(N,NA),FZ(N,NA) ATOMIC FORCES AT TIME T    **
C    ** REAL    M(NA)                      ATOMIC MASSES.             **
C    ** REAL    DSQ(NA)                    SQUARED BOND LENGTHS       **
C    ** DSQ(A) IS THE SQUARED BOND LENGTH BETWEEN ATOMS A AND A+1.    **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE SHOULD BE CALLED AFTER COMPUTING FORCES.          **
C    ** IT RETURNS THE NEW POSITIONS (TIME T+DT) IN RX,RY,RZ, AND THE **
C    ** CURRENT (TIME T) POSITIONS IN OX,OY,OZ.                       **
C    ** IT ALSO COMPUTES THE KINETIC ENERGY K AND THE CONSTRAINT      **
C    ** CONTRIBUTION TO THE ATOMIC VIRIAL WC.                         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     NA
        PARAMETER ( NA = 3 )

        INTEGER     MAXIT
        REAL        TOL, DT, K, WC, BOX
        REAL        RX(N,NA), RY(N,NA), RZ(N,NA)
        REAL        OX(N,NA), OY(N,NA), OZ(N,NA)
        REAL        FX(N,NA), FY(N,NA), FZ(N,NA)
        REAL        M(NA), DSQ(NA)

        INTEGER     I, IT
        REAL        RX1, RY1, RZ1, RX2, RY2, RZ2, RX3, RY3, RZ3
        REAL        PX1, PY1, PZ1, PX2, PY2, PZ2, PX3, PY3, PZ3
        REAL        RX12, RY12, RZ12, RX23, RY23, RZ23
        REAL        PX12, PY12, PZ12, PX23, PY23, PZ23
        REAL        VXIA, VYIA, VZIA, OM1, OM2, OM3, DTSQ
        REAL        R12SQ, R23SQ, R12R23
        REAL        P12SQ, P23SQ
        REAL        P12R12, P12R23, P23R12, P23R23
        REAL        L12, L23, L12NEW, L23NEW
        REAL        MAT11, MAT12, MAT21, MAT22
        REAL        INV11, INV12, INV21, INV22
        REAL        CONST1, CONST2, QUAD1, QUAD2, VEC1, VEC2
        REAL        Q111, Q122, Q112, Q211, Q222, Q212
        REAL        DETERM, LTOL, BOXINV
        LOGICAL     DONE

C    *******************************************************************

        BOXINV = 1.0 / BOX
        DTSQ   = DT ** 2
        LTOL   = TOL * DTSQ
        R12SQ  = DSQ(1)
        R23SQ  = DSQ(2)
        OM1    = 1.0 / M(1)
        OM2    = 1.0 / M(2)
        OM3    = 1.0 / M(3)

        K  = 0.0
        WC = 0.0

C    ** LOOP OVER MOLECULES **

        DO 2000 I = 1, N

C       ** VERLET ALGORITHM **

           RX1 = RX(I,1)
           RY1 = RY(I,1)
           RZ1 = RZ(I,1)
           PX1 = 2.0 * RX1 - OX(I,1) + DTSQ * FX(I,1) * OM1
           PY1 = 2.0 * RY1 - OY(I,1) + DTSQ * FY(I,1) * OM1
           PZ1 = 2.0 * RZ1 - OZ(I,1) + DTSQ * FZ(I,1) * OM1
           RX2 = RX(I,2)
           RY2 = RY(I,2)
           RZ2 = RZ(I,2)
           PX2 = 2.0 * RX2 - OX(I,2) + DTSQ * FX(I,2) * OM2
           PY2 = 2.0 * RY2 - OY(I,2) + DTSQ * FY(I,2) * OM2
           PZ2 = 2.0 * RZ2 - OZ(I,2) + DTSQ * FZ(I,2) * OM2
           RX3 = RX(I,3)
           RY3 = RY(I,3)
           RZ3 = RZ(I,3)
           PX3 = 2.0 * RX3 - OX(I,3) + DTSQ * FX(I,3) * OM3
           PY3 = 2.0 * RY3 - OY(I,3) + DTSQ * FY(I,3) * OM3
           PZ3 = 2.0 * RZ3 - OZ(I,3) + DTSQ * FZ(I,3) * OM3

C       ** CALCULATE RELATIVE VECTORS **

           RX12 = RX1 - RX2
           RX12 = RX12 - ANINT ( RX12 * BOXINV ) * BOX
           RY12 = RY1 - RY2
           RY12 = RY12 - ANINT ( RY12 * BOXINV ) * BOX
           RZ12 = RZ1 - RZ2
           RZ12 = RZ12 - ANINT ( RZ12 * BOXINV ) * BOX
           RX23 = RX2 - RX3
           RX23 = RX23 - ANINT ( RX23 * BOXINV ) * BOX
           RY23 = RY2 - RY3
           RY23 = RY23 - ANINT ( RY23 * BOXINV ) * BOX
           RZ23 = RZ2 - RZ3
           RZ23 = RZ23 - ANINT ( RZ23 * BOXINV ) * BOX
           PX12 = PX1 - PX2
           PX12 = PX12 - ANINT ( PX12 * BOXINV ) * BOX
           PY12 = PY1 - PY2
           PY12 = PY12 - ANINT ( PY12 * BOXINV ) * BOX
           PZ12 = PZ1 - PZ2
           PZ12 = PZ12 - ANINT ( PZ12 * BOXINV ) * BOX
           PX23 = PX2 - PX3
           PX23 = PX23 - ANINT ( PX23 * BOXINV ) * BOX
           PY23 = PY2 - PY3
           PY23 = PY23 - ANINT ( PY23 * BOXINV ) * BOX
           PZ23 = PZ2 - PZ3
           PZ23 = PZ23 - ANINT ( PZ23 * BOXINV ) * BOX

C       ** CALCULATE SCALAR PRODUCTS **

           R12R23 = RX12 * RX23 + RY12 * RY23 + RZ12 * RZ23
           P12SQ  = PX12 ** 2 + PY12 ** 2 + PZ12 ** 2
           P23SQ  = PX23 ** 2 + PY23 ** 2 + PZ23 ** 2
           P12R12 = PX12 * RX12 + PY12 * RY12 + PZ12 * RZ12
           P12R23 = PX12 * RX23 + PY12 * RY23 + PZ12 * RZ23
           P23R12 = PX23 * RX12 + PY23 * RY12 + PZ23 * RZ12
           P23R23 = PX23 * RX23 + PY23 * RY23 + PZ23 * RZ23
           CONST1 = R12SQ - P12SQ
           CONST2 = R23SQ - P23SQ

C       ** EVALUATE MATRIX ELEMENTS AND QUADRATIC COEFFICIENTS **

           MAT11 =   2.0 * P12R12 * ( OM1 + OM2 )
           MAT12 = - 2.0 * P12R23 *   OM2
           MAT21 = - 2.0 * P23R12 *   OM2
           MAT22 =   2.0 * P23R23 * ( OM2 + OM3 )

           Q111 = - R12SQ  * ( OM1 + OM2 ) ** 2
           Q122 = - R23SQ  *   OM2 ** 2
           Q112 = + 2.0 * R12R23 * ( OM1 + OM2 ) * OM2
           Q211 = - R12SQ  *   OM2 ** 2
           Q222 = - R23SQ  * ( OM2 + OM3 ) ** 2
           Q212 = + 2.0 * R12R23 * ( OM2 + OM3 ) * OM2

C       ** INVERT MATRIX **

           DETERM = 1.0 / ( MAT11 * MAT22 - MAT21 * MAT12 )
           INV11 =  MAT22 * DETERM
           INV12 = -MAT12 * DETERM
           INV21 = -MAT21 * DETERM
           INV22 =  MAT11 * DETERM

C       ** PREPARE FOR ITERATIVE LOOP **

           L12  = 0.0
           L23  = 0.0
           DONE = .FALSE.
           IT   = 0

C       ** BEGIN ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              QUAD1 = Q111 * L12 ** 2
     :              + Q122 * L23 ** 2
     :              + Q112 * L12 * L23
              QUAD2 = Q211 * L12 ** 2
     :              + Q222 * L23 ** 2
     :              + Q212 * L12 * L23

              VEC1 = CONST1 + QUAD1
              VEC2 = CONST2 + QUAD2

              L12NEW = INV11 * VEC1 + INV12 * VEC2
              L23NEW = INV21 * VEC1 + INV22 * VEC2

              DONE = ( ( ABS ( L12NEW - L12 ) .LE. LTOL ) .AND.
     :                 ( ABS ( L23NEW - L23 ) .LE. LTOL ) )

              L12 = L12NEW
              L23 = L23NEW
              IT  = IT + 1

              GOTO 1000

           ENDIF

C       ** END OF ITERATION **

           PX1 = PX1 + OM1 * ( L12 * RX12              )
           PY1 = PY1 + OM1 * ( L12 * RY12              )
           PZ1 = PZ1 + OM1 * ( L12 * RZ12              )
           PX2 = PX2 + OM2 * ( L23 * RX23 - L12 * RX12 )
           PY2 = PY2 + OM2 * ( L23 * RY23 - L12 * RY12 )
           PZ2 = PZ2 + OM2 * ( L23 * RZ23 - L12 * RZ12 )
           PX3 = PX3 + OM3 * (            - L23 * RX23 )
           PY3 = PY3 + OM3 * (            - L23 * RY23 )
           PZ3 = PZ3 + OM3 * (            - L23 * RZ23 )

           IF ( .NOT. DONE ) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
              WRITE(*,'('' MOLECULE '',I5)') I
              STOP

           ENDIF

C       ** CALCULATE KINETIC ENERGY CONTRIBUTION **

           VXIA = 0.5 * ( PX1 - OX(I,1) ) / DT
           VYIA = 0.5 * ( PY1 - OY(I,1) ) / DT
           VZIA = 0.5 * ( PZ1 - OZ(I,1) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(1)
           VXIA = 0.5 * ( PX2 - OX(I,2) ) / DT
           VYIA = 0.5 * ( PY2 - OY(I,2) ) / DT
           VZIA = 0.5 * ( PZ2 - OZ(I,2) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(2)
           VXIA = 0.5 * ( PX3 - OX(I,3) ) / DT
           VYIA = 0.5 * ( PY3 - OY(I,3) ) / DT
           VZIA = 0.5 * ( PZ3 - OZ(I,3) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(3)

C       ** CALCULATE CONSTRAINT VIRIAL CONTRIBUTION **

           WC = WC + L12 * R12SQ + L23 * R23SQ

C       ** STORE AWAY POSITIONS **

           OX(I,1) = RX1
           OY(I,1) = RY1
           OZ(I,1) = RZ1
           OX(I,2) = RX2
           OY(I,2) = RY2
           OZ(I,2) = RZ2
           OX(I,3) = RX3
           OY(I,3) = RY3
           OZ(I,3) = RZ3
           RX(I,1) = PX1
           RY(I,1) = PY1
           RZ(I,1) = PZ1
           RX(I,2) = PX2
           RY(I,2) = PY2
           RZ(I,2) = PZ2
           RX(I,3) = PX3
           RY(I,3) = PY3
           RZ(I,3) = PZ3

2000    CONTINUE

C    ** END OF LOOP OVER MOLECULES **

        K  = 0.5 * K
        WC = WC / DTSQ / 3.0

        RETURN
        END



C    *******************************************************************
C    ** FICHE F.7 - PART B                                            **
C    ** CONSTRAINT DYNAMICS FOR A RIGID NONLINEAR TRIATOMIC MOLECULE. **
C    *******************************************************************



        SUBROUTINE MOVE ( DT, TOL, MAXIT, BOX, K, WC )

        COMMON / BLOCK1 / RX, RY, RZ, OX, OY, OZ, FX, FY, FZ
        COMMON / BLOCK2 / DSQ, M

C    *******************************************************************
C    ** UPDATES ATOMIC POSITIONS WITH BOND CONSTRAINTS APPLIED.       **
C    **                                                               **
C    ** THIS ROUTINE USES THE VERLET ALGORITHM AND APPLIES BOND       **
C    ** LENGTH CONSTRAINTS TO ALL THREE BONDS, NAMELY 1-2, 2-3, 3-1.  **
C    ** WE SOLVE A SYSTEM OF QUADRATIC EQUATIONS ITERATIVELY, USING   **
C    ** MATRIX INVERSION TO SOLVE ASSOCIATED LINEAR EQUATIONS,        **
C    ** AND ITERATING TO CONVERGENCE.                                 **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** RYCKAERT ET AL., J. COMP. PHYS, 23, 327, 1977.                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                          NUMBER OF MOLECULES        **
C    ** INTEGER NA                         ATOMS PER MOL. (3 HERE)    **
C    ** REAL    DT                         TIMESTEP                   **
C    ** INTEGER MAXIT                      MAX NUMBER OF ITERATIONS   **
C    ** REAL    BOX                        BOX LENGTH                 **
C    ** REAL    TOL                        PRESCRIBED TOLERANCE       **
C    ** REAL    K                          KINETIC ENERGY             **
C    ** REAL    WC                         CONSTRAINT VIRIAL          **
C    ** REAL    RX(N,NA),RY(N,NA),RZ(N,NA) ATOMIC POSITIONS AT TIME T **
C    ** REAL    OX(N,NA),OY(N,NA),OZ(N,NA) OLD POSITIONS AT TIME T-DT **
C    ** REAL    FX(N,NA),FY(N,NA),FZ(N,NA) ATOMIC FORCES AT TIME T    **
C    ** REAL    M(NA)                      ATOMIC MASSES.             **
C    ** REAL    DSQ(NA)                    SQUARED BOND LENGTHS       **
C    ** DSQ(A) IS THE SQUARED BOND LENGTH BETWEEN ATOMS A AND A+1.    **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE SHOULD BE CALLED AFTER COMPUTING FORCES.          **
C    ** IT RETURNS THE NEW POSITIONS (TIME T+DT) IN RX,RY,RZ, AND THE **
C    ** CURRENT (TIME T) POSITIONS IN OX,OY,OZ.                       **
C    ** IT ALSO COMPUTES THE KINETIC ENERGY K AND THE CONSTRAINT      **
C    ** CONTRIBUTION TO THE VIRIAL WC.                                **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     NA
        PARAMETER ( NA = 3 )

        INTEGER     MAXIT
        REAL        TOL, DT, K, WC, BOX
        REAL        RX(N,NA), RY(N,NA), RZ(N,NA)
        REAL        OX(N,NA), OY(N,NA), OZ(N,NA)
        REAL        FX(N,NA), FY(N,NA), FZ(N,NA)
        REAL        M(NA), DSQ(NA)

        INTEGER     IT, I
        REAL        RX1, RY1, RZ1, RX2, RY2, RZ2, RX3, RY3, RZ3
        REAL        PX1, PY1, PZ1, PX2, PY2, PZ2, PX3, PY3, PZ3
        REAL        RX12, RY12, RZ12, RX23, RY23, RZ23, RX31, RY31, RZ31
        REAL        PX12, PY12, PZ12, PX23, PY23, PZ23, PX31, PY31, PZ31
        REAL        VXIA, VYIA, VZIA, OM1, OM2, OM3, DTSQ
        REAL        R12SQ, R23SQ, R31SQ, R12R23, R23R31, R31R12
        REAL        P12SQ, P23SQ, P31SQ
        REAL        P12R12, P12R23, P12R31
        REAL        P23R12, P23R23, P23R31
        REAL        P31R12, P31R23, P31R31
        REAL        L12, L23, L31, L12NEW, L23NEW, L31NEW, LTOL, BOXINV
        REAL        CONST1, CONST2, CONST3
        REAL        QUAD1, QUAD2, QUAD3, VEC1, VEC2, VEC3
        REAL        Q111, Q122, Q133, Q112, Q123, Q131
        REAL        Q211, Q222, Q233, Q212, Q223, Q231
        REAL        Q311, Q322, Q333, Q312, Q323, Q331
        REAL        MAT(3,3), INV(3,3)
        LOGICAL     DONE

C    *******************************************************************

        BOXINV = 1.0 / BOX
        DTSQ   = DT ** 2
        LTOL   = TOL * DTSQ
        R12SQ  = DSQ(1)
        R23SQ  = DSQ(2)
        R31SQ  = DSQ(3)
        OM1    = 1.0 / M(1)
        OM2    = 1.0 / M(2)
        OM3    = 1.0 / M(3)

        K  = 0.0
        WC = 0.0

C    ** LOOP OVER MOLECULES **

        DO 2000 I = 1, N

C       ** VERLET ALGORITHM **

           RX1 = RX(I,1)
           RY1 = RY(I,1)
           RZ1 = RZ(I,1)
           PX1 = 2.0 * RX1 - OX(I,1) + DTSQ * FX(I,1) * OM1
           PY1 = 2.0 * RY1 - OY(I,1) + DTSQ * FY(I,1) * OM1
           PZ1 = 2.0 * RZ1 - OZ(I,1) + DTSQ * FZ(I,1) * OM1
           RX2 = RX(I,2)
           RY2 = RY(I,2)
           RZ2 = RZ(I,2)
           PX2 = 2.0 * RX2 - OX(I,2) + DTSQ * FX(I,2) * OM2
           PY2 = 2.0 * RY2 - OY(I,2) + DTSQ * FY(I,2) * OM2
           PZ2 = 2.0 * RZ2 - OZ(I,2) + DTSQ * FZ(I,2) * OM2
           RX3 = RX(I,3)
           RY3 = RY(I,3)
           RZ3 = RZ(I,3)
           PX3 = 2.0 * RX3 - OX(I,3) + DTSQ * FX(I,3) * OM3
           PY3 = 2.0 * RY3 - OY(I,3) + DTSQ * FY(I,3) * OM3
           PZ3 = 2.0 * RZ3 - OZ(I,3) + DTSQ * FZ(I,3) * OM3

C       ** CALCULATE RELATIVE VECTORS **

           RX12 = RX1 - RX2
           RX12 = RX12 - ANINT ( RX12 * BOXINV ) * BOX
           RY12 = RY1 - RY2
           RY12 = RY12 - ANINT ( RY12 * BOXINV ) * BOX
           RZ12 = RZ1 - RZ2
           RZ12 = RZ12 - ANINT ( RZ12 * BOXINV ) * BOX
           RX23 = RX2 - RX3
           RX23 = RX23 - ANINT ( RX23 * BOXINV ) * BOX
           RY23 = RY2 - RY3
           RY23 = RY23 - ANINT ( RY23 * BOXINV ) * BOX
           RZ23 = RZ2 - RZ3
           RZ23 = RZ23 - ANINT ( RZ23 * BOXINV ) * BOX
           RX31 = RX3 - RX1
           RX31 = RX31 - ANINT ( RX31 * BOXINV ) * BOX
           RY31 = RY3 - RY1
           RY31 = RY31 - ANINT ( RY31 * BOXINV ) * BOX
           RZ31 = RZ3 - RZ1
           RZ31 = RZ31 - ANINT ( RZ31 * BOXINV ) * BOX

           PX12 = PX1 - PX2
           PX12 = PX12 - ANINT ( PX12 * BOXINV ) * BOX
           PY12 = PY1 - PY2
           PY12 = PY12 - ANINT ( PY12 * BOXINV ) * BOX
           PZ12 = PZ1 - PZ2
           PZ12 = PZ12 - ANINT ( PZ12 * BOXINV ) * BOX
           PX23 = PX2 - PX3
           PX23 = PX23 - ANINT ( PX23 * BOXINV ) * BOX
           PY23 = PY2 - PY3
           PY23 = PY23 - ANINT ( PY23 * BOXINV ) * BOX
           PZ23 = PZ2 - PZ3
           PZ23 = PZ23 - ANINT ( PZ23 * BOXINV ) * BOX
           PX31 = PX3 - PX1
           PX31 = PX31 - ANINT ( PX31 * BOXINV ) * BOX
           PY31 = PY3 - PY1
           PY31 = PY31 - ANINT ( PY31 * BOXINV ) * BOX
           PZ31 = PZ3 - PZ1
           PZ31 = PZ31 - ANINT ( PZ31 * BOXINV ) * BOX

C       ** CALCULATE SCALAR PRODUCTS **

           R12R23 = RX12 * RX23 + RY12 * RY23 + RZ12 * RZ23
           R23R31 = RX23 * RX31 + RY23 * RY31 + RZ23 * RZ31
           R31R12 = RX31 * RX12 + RY31 * RY12 + RZ31 * RZ12
           P12SQ  = PX12 ** 2 + PY12 ** 2 + PZ12 ** 2
           P23SQ  = PX23 ** 2 + PY23 ** 2 + PZ23 ** 2
           P31SQ  = PX31 ** 2 + PY31 ** 2 + PZ31 ** 2
           P12R12 = PX12 * RX12 + PY12 * RY12 + PZ12 * RZ12
           P12R23 = PX12 * RX23 + PY12 * RY23 + PZ12 * RZ23
           P12R31 = PX12 * RX31 + PY12 * RY31 + PZ12 * RZ31
           P23R12 = PX23 * RX12 + PY23 * RY12 + PZ23 * RZ12
           P23R23 = PX23 * RX23 + PY23 * RY23 + PZ23 * RZ23
           P23R31 = PX23 * RX31 + PY23 * RY31 + PZ23 * RZ31
           P31R12 = PX31 * RX12 + PY31 * RY12 + PZ31 * RZ12
           P31R23 = PX31 * RX23 + PY31 * RY23 + PZ31 * RZ23
           P31R31 = PX31 * RX31 + PY31 * RY31 + PZ31 * RZ31
           CONST1 = R12SQ - P12SQ
           CONST2 = R23SQ - P23SQ
           CONST3 = R31SQ - P31SQ

C       ** CALCULATE MATRIX AND QUADRATIC COEFFICIENTS **

           MAT(1,1) =  2.0 * ( OM1 + OM2 ) * P12R12
           MAT(1,2) = -2.0 *   OM2         * P12R23
           MAT(1,3) = -2.0 *   OM1         * P12R31
           MAT(2,1) = -2.0 *   OM2         * P23R12
           MAT(2,2) =  2.0 * ( OM2 + OM3 ) * P23R23
           MAT(2,3) = -2.0 *   OM3         * P23R31
           MAT(3,1) = -2.0 *   OM1         * P31R12
           MAT(3,2) = -2.0 *   OM3         * P31R23
           MAT(3,3) =  2.0 * ( OM1 + OM3 ) * P31R31

           Q111 = - R12SQ * ( OM1 + OM2 ) ** 2
           Q122 = - R23SQ * OM2 ** 2
           Q133 = - R31SQ * OM1 ** 2
           Q112 = + 2.0 * R12R23 * ( OM1 + OM2 ) * OM2
           Q123 = - 2.0 * R23R31 * OM1 * OM2
           Q131 = + 2.0 * R31R12 * ( OM1 + OM2 ) * OM1

           Q211 = - R12SQ * OM2 ** 2
           Q222 = - R23SQ * ( OM2 + OM3 ) ** 2
           Q233 = - R31SQ * OM3 ** 2
           Q212 = + 2.0 * R12R23 * ( OM2 + OM3 ) * OM2
           Q223 = + 2.0 * R23R31 * ( OM2 + OM3 ) * OM3
           Q231 = - 2.0 * R31R12 * OM2 * OM3

           Q311 = - R12SQ * OM1 ** 2
           Q322 = - R23SQ * OM3 ** 2
           Q333 = - R31SQ * ( OM1 + OM3 ) ** 2
           Q312 = - 2.0 * R12R23 * OM1 * OM3
           Q323 = + 2.0 * R23R31 * ( OM1 + OM3 ) * OM3
           Q331 = + 2.0 * R31R12 * ( OM1 + OM3 ) * OM1

C       ** NOW CALL ROUTINE TO INVERT THE MATRIX MAT **

           CALL MATINV ( MAT, INV )

C       ** PREPARE FOR ITERATIVE LOOP **

           DONE = .FALSE.
           IT   = 0
           L12  = 0.0
           L23  = 0.0
           L31  = 0.0

C       ** BEGIN ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              QUAD1 = Q111 * L12 ** 2 + Q112 * L12 * L23
     :              + Q122 * L23 ** 2 + Q123 * L23 * L31
     :              + Q133 * L31 ** 2 + Q131 * L31 * L12
              QUAD2 = Q211 * L12 ** 2 + Q212 * L12 * L23
     :              + Q222 * L23 ** 2 + Q223 * L23 * L31
     :              + Q233 * L31 ** 2 + Q231 * L31 * L12
              QUAD3 = Q311 * L12 ** 2 + Q312 * L12 * L23
     :              + Q322 * L23 ** 2 + Q323 * L23 * L31
     :              + Q333 * L31 ** 2 + Q331 * L31 * L12

              VEC1 = CONST1 + QUAD1
              VEC2 = CONST2 + QUAD2
              VEC3 = CONST3 + QUAD3

C          ** OBTAIN SOLUTIONS OF LINEARIZED EQUATION **

              L12NEW = INV(1,1) * VEC1 + INV(1,2) * VEC2
     :               + INV(1,3) * VEC3
              L23NEW = INV(2,1) * VEC1 + INV(2,2) * VEC2
     :               + INV(2,3) * VEC3
              L31NEW = INV(3,1) * VEC1 + INV(3,2) * VEC2
     :               + INV(3,3) * VEC3

              DONE = ( ( ABS ( L12NEW - L12 ) .LE. LTOL ) .AND.
     :                 ( ABS ( L23NEW - L23 ) .LE. LTOL ) .AND.
     :                 ( ABS ( L31NEW - L31 ) .LE. LTOL ) )

              L12 = L12NEW
              L23 = L23NEW
              L31 = L31NEW
              IT  = IT + 1

              GOTO 1000

           ENDIF

C       ** END OF ITERATION **

           PX1 = PX1 + OM1 * ( L12 * RX12 - L31 * RX31 )
           PY1 = PY1 + OM1 * ( L12 * RY12 - L31 * RY31 )
           PZ1 = PZ1 + OM1 * ( L12 * RZ12 - L31 * RZ31 )
           PX2 = PX2 + OM2 * ( L23 * RX23 - L12 * RX12 )
           PY2 = PY2 + OM2 * ( L23 * RY23 - L12 * RY12 )
           PZ2 = PZ2 + OM2 * ( L23 * RZ23 - L12 * RZ12 )
           PX3 = PX3 + OM3 * ( L31 * RX31 - L23 * RX23 )
           PY3 = PY3 + OM3 * ( L31 * RY31 - L23 * RY23 )
           PZ3 = PZ3 + OM3 * ( L31 * RZ31 - L23 * RZ23 )

           IF ( .NOT. DONE ) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
              WRITE(*,'('' MOLECULE '',I5)') I
              STOP

           ENDIF

C       ** CALCULATE KINETIC ENERGY CONTRIBUTION **

           VXIA = 0.5 * ( PX1 - OX(I,1) ) / DT
           VYIA = 0.5 * ( PY1 - OY(I,1) ) / DT
           VZIA = 0.5 * ( PZ1 - OZ(I,1) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(1)
           VXIA = 0.5 * ( PX2 - OX(I,2) ) / DT
           VYIA = 0.5 * ( PY2 - OY(I,2) ) / DT
           VZIA = 0.5 * ( PZ2 - OZ(I,2) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(2)
           VXIA = 0.5 * ( PX3 - OX(I,3) ) / DT
           VYIA = 0.5 * ( PY3 - OY(I,3) ) / DT
           VZIA = 0.5 * ( PZ3 - OZ(I,3) ) / DT
           K = K + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * M(3)

C       ** CALCULATE CONSTRAINT VIRIAL CONTRIBUTION **

           WC = WC + L12 * R12SQ + L23 * R23SQ + L31 * R31SQ

C       ** STORE RESULTS **

           OX(I,1) = RX1
           OY(I,1) = RY1
           OZ(I,1) = RZ1
           RX(I,1) = PX1
           RY(I,1) = PY1
           RZ(I,1) = PZ1
           OX(I,2) = RX2
           OY(I,2) = RY2
           OZ(I,2) = RZ2
           RX(I,2) = PX2
           RY(I,2) = PY2
           RZ(I,2) = PZ2
           OX(I,3) = RX3
           OY(I,3) = RY3
           OZ(I,3) = RZ3
           RX(I,3) = PX3
           RY(I,3) = PY3
           RZ(I,3) = PZ3

2000    CONTINUE

C    ** END OF LOOP OVER MOLECULES **

        K  = 0.5 * K
        WC = WC / DTSQ / 3.0

        RETURN
        END



        SUBROUTINE MATINV ( MAT, INV )

C    *******************************************************************
C    ** A SIMPLE ROUTINE TO INVERT THE 3X3 MATRIX MAT AND RETURN INV  **
C    *******************************************************************

        REAL        MAT(3,3), INV(3,3), DETERM
        REAL        TOL
        PARAMETER ( TOL = 1.E-7 )

C    *******************************************************************

C    ** GET ALL SIGNED COFACTORS, TRANSPOSE, AND PUT IN INV **

        INV(1,1) = MAT(2,2) * MAT(3,3) - MAT(2,3) * MAT(3,2)
        INV(2,1) = MAT(3,1) * MAT(2,3) - MAT(3,3) * MAT(2,1)
        INV(3,1) = MAT(2,1) * MAT(3,2) - MAT(2,2) * MAT(3,1)
        INV(1,2) = MAT(3,2) * MAT(1,3) - MAT(3,3) * MAT(1,2)
        INV(2,2) = MAT(1,1) * MAT(3,3) - MAT(1,3) * MAT(3,1)
        INV(3,2) = MAT(3,1) * MAT(1,2) - MAT(3,2) * MAT(1,1)
        INV(1,3) = MAT(1,2) * MAT(2,3) - MAT(1,3) * MAT(2,2)
        INV(2,3) = MAT(2,1) * MAT(1,3) - MAT(2,3) * MAT(1,1)
        INV(3,3) = MAT(1,1) * MAT(2,2) - MAT(1,2) * MAT(2,1)

C    ** GET DETERMINANT AND GUARD AGAINST ZERO **

        DETERM = MAT(1,1) * INV(1,1) + MAT(1,2) * INV(2,1)
     :         + MAT(1,3) * INV(3,1)

        IF ( ABS ( DETERM ) .LT. TOL ) STOP 'ZERO DETERM IN MATINV'

        INV(1,1) = INV(1,1) / DETERM
        INV(1,2) = INV(1,2) / DETERM
        INV(1,3) = INV(1,3) / DETERM
        INV(2,1) = INV(2,1) / DETERM
        INV(2,2) = INV(2,2) / DETERM
        INV(2,3) = INV(2,3) / DETERM
        INV(3,1) = INV(3,1) / DETERM
        INV(3,2) = INV(3,2) / DETERM
        INV(3,3) = INV(3,3) / DETERM

        RETURN
        END



