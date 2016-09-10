********************************************************************************
** FICHE F.5.  QUATERNION PARAMETER PREDICTOR-CORRECTOR ALGORITHM             **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** RIGID MOLECULE ROTATION USING QUATERNION PREDICTOR-CORRECTOR. **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** EVANS AND MURAD, MOLEC. PHYS. 34, 327, 1977.                  **
C    ** GEAR, NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY            **
C    ** DIFFERENTIAL EQUATIONS, PRENTICE-HALL, 1971.                  **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE PREDIC ( DT )                                      **
C    **    PREDICTS POSITIONS, VELOCITIES ETC. AT NEXT STEP.          **
C    ** SUBROUTINE MOLATM                                             **
C    **    CONVERTS MOLECULAR COORDINATES TO ATOM POSITIONS.          **
C    ** SUBROUTINE ATMMOL                                             **
C    **    CONVERTS ATOMIC FORCES TO MOLECULAR FORCES AND TORQUES.    **
C    ** SUBROUTINE CORREC ( DT, M, IXX, IYY, IZZ, K )                 **
C    **    CORRECTS POSITIONS, VELOCITIES ETC.                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** REAL    DT                            TIMESTEP                **
C    ** REAL    M                             MOLECULAR MASS          **
C    ** REAL    IXX,IYY,IZZ                   PRINCIPAL INERTIAS      **
C    ** REAL    K                             KINETIC ENERGY          **
C    ** REAL    RX (N),RY (N),RZ (N)          C-O-M POSITIONS         **
C    ** REAL    RX1(N),RY1(N),RZ1(N)          FIRST DERIVATIVES       **
C    ** REAL    RX2(N),RY2(N),RZ2(N)          SECOND DERIVATIVES      **
C    ** REAL    RX3(N),RY3(N),RZ3(N)          THIRD DERIVATIVES       **
C    ** REAL    FX (N),FY (N),FZ (N)          TOTAL FORCES            **
C    ** REAL    QW (N),QX (N),QY (N),QZ (N)   QUATERNION PARAMETERS   **
C    ** REAL    QW1(N),QX1(N),QY1(N),QZ1(N)   FIRST DERIVATIVES       **
C    ** REAL    QW2(N),QX2(N),QY2(N),QZ2(N)   SECOND DERIVATIVES      **
C    ** REAL    QW3(N),QX3(N),QY3(N),QZ3(N)   THIRD DERIVATIVES       **
C    ** REAL    QW4(N),QX4(N),QY4(N),QZ4(N)   FOURTH DERIVATIVES      **
C    ** REAL    OX (N),OY (N),OZ (N)          ANGULAR VELOCITIES      **
C    ** REAL    OX1(N),OY1(N),OZ1(N)          FIRST DERIVATIVES       **
C    ** REAL    OX2(N),OY2(N),OZ2(N)          SECOND DERIVATIVES      **
C    ** REAL    OX3(N),OY3(N),OZ3(N)          THIRD DERIVATIVES       **
C    ** REAL    OX4(N),OY4(N),OZ4(N)          FOURTH DERIVATIVES      **
C    ** REAL    TX (N),TY (N),TZ (N)          TOTAL TORQUES           **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE PREDICTOR ROUTINE IS CALLED TO ADVANCE THE POSITIONS,     **
C    ** QUATERNIONS, AND ANGULAR VELOCITIES GIVEN THE CURRENT VALUES  **
C    ** OF THESE QUANTITIES AND THEIR SUCCESSIVE TIME DERIVATIVES.    **
C    ** FOLLOWING THIS, SUBROUTINE MOLATM USES THE QUATERNIONS TO     **
C    ** OBTAIN THE ATOM OR SITE POSITIONS RSX,RSY,RSZ.  THESE ARE FED **
C    ** INTO A FORCE ROUTINE (NOT SUPPLIED HERE) WHICH GIVES THE      **
C    ** FORCE ACTING ON EACH SITE OR ATOM. IN TURN, THESE ARE         **
C    ** CONVERTED INTO THE TOTAL FORCE AND TORQUE ACTING ON EACH      **
C    ** MOLECULE, BY SUBROUTINE ATMMOL, AND THESE ARE USED IN THE     **
C    ** CORRECTOR STAGE.                                              **
C    *******************************************************************



        SUBROUTINE PREDIC ( DT )

        COMMON / BLOCK1 / RX , RY , RZ , RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX , FY , FZ
        COMMON / BLOCK2 / QW , QX , QY , QZ , QW1, QX1, QY1, QZ1,
     :                    QW2, QX2, QY2, QZ2, QW3, QX3, QY3, QZ3,
     :                    QW4, QX4, QY4, QZ4,
     :                    OX , OY , OZ , OX1, OY1, OZ1,
     :                    OX2, OY2, OZ2, OX3, OY3, OZ3,
     :                    OX4, OY4, OZ4, TX, TY, TZ

C    *******************************************************************
C    ** PREDICTOR ROUTINE                                             **
C    **                                                               **
C    ** WE ADOPT A 5-VALUE METHOD FOR REORIENTATIONAL VARIABLES       **
C    ** EMPLOYING BODY-FIXED ANGULAR VELOCITIES AND QUATERNIONS,      **
C    ** AND A 4-VALUE METHOD FOR CENTRE-OF-MASS (C-O-M) TRANSLATION.  **
C    ** THE PREDICTOR STAGE IS A SIMPLE TAYLOR SERIES.                **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT
        REAL        RX (N), RY (N), RZ (N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX (N), FY (N), FZ (N)
        REAL        QW (N), QX (N), QY (N), QZ (N)
        REAL        QW1(N), QX1(N), QY1(N), QZ1(N)
        REAL        QW2(N), QX2(N), QY2(N), QZ2(N)
        REAL        QW3(N), QX3(N), QY3(N), QZ3(N)
        REAL        QW4(N), QX4(N), QY4(N), QZ4(N)
        REAL        OX (N), OY (N), OZ (N)
        REAL        OX1(N), OY1(N), OZ1(N)
        REAL        OX2(N), OY2(N), OZ2(N)
        REAL        OX3(N), OY3(N), OZ3(N)
        REAL        OX4(N), OY4(N), OZ4(N)
        REAL        TX (N), TY (N), TZ (N)

        INTEGER     I
        REAL        C1, C2, C3, C4

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        DO 100 I = 1, N

           RX (I) = RX (I) + C1*RX1(I) + C2*RX2(I) + C3*RX3(I)
           RY (I) = RY (I) + C1*RY1(I) + C2*RY2(I) + C3*RY3(I)
           RZ (I) = RZ (I) + C1*RZ1(I) + C2*RZ2(I) + C3*RZ3(I)
           RX1(I) = RX1(I) + C1*RX2(I) + C2*RX3(I)
           RY1(I) = RY1(I) + C1*RY2(I) + C2*RY3(I)
           RZ1(I) = RZ1(I) + C1*RZ2(I) + C2*RZ3(I)
           RX2(I) = RX2(I) + C1*RX3(I)
           RY2(I) = RY2(I) + C1*RY3(I)
           RZ2(I) = RZ2(I) + C1*RZ3(I)

           QW(I) = QW(I) + C1*QW1(I) + C2*QW2(I) + C3*QW3(I) + C4*QW4(I)
           QX(I) = QX(I) + C1*QX1(I) + C2*QX2(I) + C3*QX3(I) + C4*QX4(I)
           QY(I) = QY(I) + C1*QY1(I) + C2*QY2(I) + C3*QY3(I) + C4*QY4(I)
           QZ(I) = QZ(I) + C1*QZ1(I) + C2*QZ2(I) + C3*QZ3(I) + C4*QZ4(I)
           QW1(I) = QW1(I) + C1*QW2(I) + C2*QW3(I) + C3*QW4(I)
           QX1(I) = QX1(I) + C1*QX2(I) + C2*QX3(I) + C3*QX4(I)
           QY1(I) = QY1(I) + C1*QY2(I) + C2*QY3(I) + C3*QY4(I)
           QZ1(I) = QZ1(I) + C1*QZ2(I) + C2*QZ3(I) + C3*QZ4(I)
           QW2(I) = QW2(I) + C1*QW3(I) + C2*QW4(I)
           QX2(I) = QX2(I) + C1*QX3(I) + C2*QX4(I)
           QY2(I) = QY2(I) + C1*QY3(I) + C2*QY4(I)
           QZ2(I) = QZ2(I) + C1*QZ3(I) + C2*QZ4(I)
           QW3(I) = QW3(I) + C1*QW4(I)
           QX3(I) = QX3(I) + C1*QX4(I)
           QY3(I) = QY3(I) + C1*QY4(I)
           QZ3(I) = QZ3(I) + C1*QZ4(I)

           OX(I) = OX(I) + C1*OX1(I) + C2*OX2(I) + C3*OX3(I) + C4*OX4(I)
           OY(I) = OY(I) + C1*OY1(I) + C2*OY2(I) + C3*OY3(I) + C4*OY4(I)
           OZ(I) = OZ(I) + C1*OZ1(I) + C2*OZ2(I) + C3*OZ3(I) + C4*OZ4(I)
           OX1(I) = OX1(I) + C1*OX2(I) + C2*OX3(I) + C3*OX4(I)
           OY1(I) = OY1(I) + C1*OY2(I) + C2*OY3(I) + C3*OY4(I)
           OZ1(I) = OZ1(I) + C1*OZ2(I) + C2*OZ3(I) + C3*OZ4(I)
           OX2(I) = OX2(I) + C1*OX3(I) + C2*OX4(I)
           OY2(I) = OY2(I) + C1*OY3(I) + C2*OY4(I)
           OZ2(I) = OZ2(I) + C1*OZ3(I) + C2*OZ4(I)
           OX3(I) = OX3(I) + C1*OX4(I)
           OY3(I) = OY3(I) + C1*OY4(I)
           OZ3(I) = OZ3(I) + C1*OZ4(I)

100     CONTINUE

        RETURN
        END



        SUBROUTINE MOLATM

        COMMON / BLOCK1 / RX , RY , RZ , RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX , FY , FZ
        COMMON / BLOCK2 / QW , QX , QY , QZ , QW1, QX1, QY1, QZ1,
     :                    QW2, QX2, QY2, QZ2, QW3, QX3, QY3, QZ3,
     :                    QW4, QX4, QY4, QZ4,
     :                    OX , OY , OZ , OX1, OY1, OZ1,
     :                    OX2, OY2, OZ2, OX3, OY3, OZ3,
     :                    OX4, OY4, OZ4, TX , TY , TZ
        COMMON / BLOCK3 / RSX, RSY, RSZ, FSX, FSY, FSZ
        COMMON / BLOCK4 / DX , DY , DZ

C    *******************************************************************
C    ** CONVERSION OF MOLECULAR COORDINATES TO ATOM OR SITE POSITIONS **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    ** REAL    DX(NA),DY(NA),DZ(NA)          ATOM POSITIONS IN MOLEC **
C    ** REAL    AXX,AXY,AXZ ETC.              ROTATION MATRIX         **
C    ** THE VARIABLES DX,DY,DZ ARE ACTUALLY THE POSITION VECTORS OF   **
C    ** EACH ATOM IN THE MOLECULE RELATIVE TO THE CENTRE OF MASS IN   **
C    ** THE UNROTATED, BODY-FIXED, AXIS SYSTEM.                       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE TRANSPOSE OF THE ROTATION MATRIX IS USED TO OBTAIN THE    **
C    ** POSITIONS OF EACH ATOM FROM THE CENTRE-OF-MASS POSITION AND   **
C    ** THE BODY-FIXED ATOM POSITION VECTORS (KNOWN FROM THE START).  **
C    ** THESE MAY THEN BE FED INTO THE FORCE ROUTINE.                 **
C    ** FOR THIS EXAMPLE WE TAKE (NONLINEAR) TRIATOMIC MOLECULES.     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX (N), RY (N), RZ (N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX (N), FY (N), FZ (N)
        REAL        QW (N), QX (N), QY (N), QZ (N)
        REAL        QW1(N), QX1(N), QY1(N), QZ1(N)
        REAL        QW2(N), QX2(N), QY2(N), QZ2(N)
        REAL        QW3(N), QX3(N), QY3(N), QZ3(N)
        REAL        QW4(N), QX4(N), QY4(N), QZ4(N)
        REAL        OX (N), OY (N), OZ (N)
        REAL        OX1(N), OY1(N), OZ1(N)
        REAL        OX2(N), OY2(N), OZ2(N)
        REAL        OX3(N), OY3(N), OZ3(N)
        REAL        OX4(N), OY4(N), OZ4(N)
        REAL        TX (N), TY (N), TZ (N)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)
        REAL        DX(NA), DY(NA), DZ(NA)

        INTEGER     I, A
        REAL        AXX, AXY, AXZ, AYX, AYY, AYZ, AZX, AZY, AZZ

C    *******************************************************************

C    ** LOOP OVER ALL MOLECULES **

        DO 200 I = 1, N

C       ** CALCULATE ROTATION MATRIX ELEMENTS **

           AXX = QW(I) ** 2 + QX(I) ** 2 - QY(I) ** 2 - QZ(I) ** 2
           AXY = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
           AXZ = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
           AYX = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
           AYY = QW(I) ** 2 - QX(I) ** 2 + QY(I) ** 2 - QZ(I) ** 2
           AYZ = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
           AZX = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
           AZY = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
           AZZ = QW(I) ** 2 - QX(I) ** 2 - QY(I) ** 2 + QZ(I) ** 2

C       ** LOOP OVER ALL SITES IN MOLECULE **

           DO 199 A = 1, NA

              RSX(I,A) = RX(I) + AXX * DX(A) + AYX * DY(A) + AZX * DZ(A)
              RSY(I,A) = RY(I) + AXY * DX(A) + AYY * DY(A) + AZY * DZ(A)
              RSZ(I,A) = RZ(I) + AXZ * DX(A) + AYZ * DY(A) + AZZ * DZ(A)

199        CONTINUE

200     CONTINUE

        RETURN
        END



        SUBROUTINE ATMMOL

        COMMON / BLOCK1 / RX , RY , RZ , RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX , FY , FZ
        COMMON / BLOCK2 / QW , QX , QY , QZ , QW1, QX1, QY1, QZ1,
     :                    QW2, QX2, QY2, QZ2, QW3, QX3, QY3, QZ3,
     :                    QW4, QX4, QY4, QZ4,
     :                    OX , OY , OZ , OX1, OY1, OZ1,
     :                    OX2, OY2, OZ2, OX3, OY3, OZ3,
     :                    OX4, OY4, OZ4, TX, TY, TZ
        COMMON / BLOCK3 / RSX, RSY, RSZ, FSX, FSY, FSZ

C    *******************************************************************
C    ** CONVERTS ATOMIC SITE FORCES TO MOLECULAR FORCE AND TORQUE.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    ** REAL    FSX(N,NA),FSY(N,NA),FSZ(N,NA) ATOM FORCES             **
C    ** REAL    AXX,AXY,AXZ ETC.              ROTATION MATRIX         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROTATION MATRIX IS USED TO CONVERT THE TORQUES FROM       **
C    ** SPACE-FIXED TO BODY-FIXED AXES PRIOR TO THE CORRECTOR STEP    **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX (N), RY (N), RZ (N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX (N), FY (N), FZ (N)
        REAL        QW (N), QX (N), QY (N), QZ (N)
        REAL        QW1(N), QX1(N), QY1(N), QZ1(N)
        REAL        QW2(N), QX2(N), QY2(N), QZ2(N)
        REAL        QW3(N), QX3(N), QY3(N), QZ3(N)
        REAL        QW4(N), QX4(N), QY4(N), QZ4(N)
        REAL        OX (N), OY (N), OZ (N)
        REAL        OX1(N), OY1(N), OZ1(N)
        REAL        OX2(N), OY2(N), OZ2(N)
        REAL        OX3(N), OY3(N), OZ3(N)
        REAL        OX4(N), OY4(N), OZ4(N)
        REAL        TX (N), TY (N), TZ (N)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)

        INTEGER     I, A
        REAL        AXX, AXY, AXZ, AYX, AYY, AYZ, AZX, AZY, AZZ
        REAL        FXI, FYI, FZI, TXI, TYI, TZI
        REAL        RXI, RYI, RZI, QWI, QXI, QYI, QZI
        REAL        RSXIA, RSYIA, RSZIA, FSXIA, FSYIA, FSZIA

C    *******************************************************************

C    ** LOOP OVER MOLECULES **

        DO 300 I = 1, N

           FXI = 0.0
           FYI = 0.0
           FZI = 0.0
           TXI = 0.0
           TYI = 0.0
           TZI = 0.0
           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           QWI = QW(I)
           QXI = QX(I)
           QYI = QY(I)
           QZI = QZ(I)

C       ** LOOP OVER SITES IN A MOLECULE **

           DO 299 A = 1, NA

              FSXIA = FSX(I,A)
              FSYIA = FSY(I,A)
              FSZIA = FSZ(I,A)
              RSXIA = RSX(I,A) - RXI
              RSYIA = RSY(I,A) - RYI
              RSZIA = RSZ(I,A) - RZI

C          ** TOTAL FORCE AND TORQUE CONTRIBUTIONS **

              FXI = FXI + FSXIA
              FYI = FYI + FSYIA
              FZI = FZI + FSZIA
              TXI = TXI + RSYIA * FSZIA - RSZIA * FSYIA
              TYI = TYI + RSZIA * FSXIA - RSXIA * FSZIA
              TZI = TZI + RSXIA * FSYIA - RSYIA * FSXIA

299        CONTINUE

C       ** STORE TOTAL FORCE **

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI

C       ** CALCULATE ROTATION MATRIX ELEMENTS **

           AXX = QWI ** 2 + QXI ** 2 - QYI ** 2 - QZI ** 2
           AXY = 2.0 * ( QXI * QYI + QWI * QZI )
           AXZ = 2.0 * ( QXI * QZI - QWI * QYI )
           AYX = 2.0 * ( QXI * QYI - QWI * QZI )
           AYY = QWI ** 2 - QXI ** 2 + QYI ** 2 - QZI ** 2
           AYZ = 2.0 * ( QYI * QZI + QWI * QXI )
           AZX = 2.0 * ( QXI * QZI + QWI * QYI )
           AZY = 2.0 * ( QYI * QZI - QWI * QXI )
           AZZ = QWI ** 2 - QXI ** 2 - QYI ** 2 + QZI ** 2

C       ** CONVERT TORQUE TO BODY-FIXED COORDINATES **

           TX(I) = AXX * TXI + AXY * TYI + AXZ * TZI
           TY(I) = AYX * TXI + AYY * TYI + AYZ * TZI
           TZ(I) = AZX * TXI + AZY * TYI + AZZ * TZI

300     CONTINUE

        RETURN
        END



        SUBROUTINE CORREC ( DT, M, IXX, IYY, IZZ, K )

        COMMON / BLOCK1 / RX , RY , RZ , RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX , FY , FZ
        COMMON / BLOCK2 / QW , QX , QY , QZ , QW1, QX1, QY1, QZ1,
     :                    QW2, QX2, QY2, QZ2, QW3, QX3, QY3, QZ3,
     :                    QW4, QX4, QY4, QZ4,
     :                    OX , OY , OZ , OX1, OY1, OZ1,
     :                    OX2, OY2, OZ2, OX3, OY3, OZ3,
     :                    OX4, OY4, OZ4, TX, TY, TZ

C    *******************************************************************
C    ** CORRECTS TRANSLATIONAL AND ROTATIONAL VARIABLES.              **
C    **                                                               **
C    ** THE CORRECTOR STAGE USES GEAR COEFFICIENTS (SEE REF ABOVE).   **
C    ** FOR TIMESTEP-SCALED VARIABLES THESE WOULD BE AS FOLLOWS.      **
C    ** FOR TRANSLATIONAL ALGORITHM, 4-VALUE METHOD, 2ND-ORDER D.E.   **
C    ** COEFFICIENTS ARE 1/6, 5/6, 1, 1/3                             **
C    ** FOR ROTATIONAL ALGORITHM, 5-VALUE METHOD, 1ST-ORDER D.E.      **
C    ** COEFFICIENTS ARE 251/720, 1, 11/12, 1/3, 1/24.                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL   GEART0, GEART1, GEART3         TRANSLATIONAL COEFFTS   **
C    ** REAL   GEARR0, GEARR2, GEARR3, GEARR4 ROTATIONAL COEFFTS      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS ROUTINE IS CALLED AFTER EVALUATION OF FORCES AND TORQUES **
C    ** AND THE CONVERSION OF TORQUES INTO BODY-FIXED AXES.           **
C    ** IT ALSO RETURNS THE KINETIC ENERGY.                           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT, M, IXX, IYY, IZZ, K
        REAL        RX (N), RY (N), RZ (N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX (N), FY (N), FZ (N)
        REAL        QW (N), QX (N), QY (N), QZ (N)
        REAL        QW1(N), QX1(N), QY1(N), QZ1(N)
        REAL        QW2(N), QX2(N), QY2(N), QZ2(N)
        REAL        QW3(N), QX3(N), QY3(N), QZ3(N)
        REAL        QW4(N), QX4(N), QY4(N), QZ4(N)
        REAL        OX (N), OY (N), OZ (N)
        REAL        OX1(N), OY1(N), OZ1(N)
        REAL        OX2(N), OY2(N), OZ2(N)
        REAL        OX3(N), OY3(N), OZ3(N)
        REAL        OX4(N), OY4(N), OZ4(N)
        REAL        TX (N), TY (N), TZ (N)

        INTEGER     I
        REAL        C1, C2, C3, C4
        REAL        CORRW, CORRX, CORRY, CORRZ
        REAL        RX2I, RY2I, RZ2I
        REAL        QW1I, QX1I, QY1I, QZ1I, OX1I, OY1I, OZ1I
        REAL        CTRAN0, CTRAN1, CTRAN3
        REAL        CROT0, CROT2, CROT3, CROT4

        REAL        GEART0, GEART1, GEART3

        PARAMETER ( GEART0 = 1.0 / 6.0,
     :              GEART1 = 5.0 / 6.0,
     :              GEART3 = 1.0 / 3.0 )

        REAL        GEARR0, GEARR2, GEARR3, GEARR4
        PARAMETER ( GEARR0 = 251.0 / 720.0,
     :              GEARR2 = 11.0  / 12.0,
     :              GEARR3 = 1.0   / 3.0,
     :              GEARR4 = 1.0   / 24.0 )

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        CTRAN0 = GEART0 * C2
        CTRAN1 = GEART1 * C2 / C1
        CTRAN3 = GEART3 * C2 / C3

        CROT0 = GEARR0 * C1
        CROT2 = GEARR2 * C1 / C2
        CROT3 = GEARR3 * C1 / C3
        CROT4 = GEARR4 * C1 / C4

        DO 400 I = 1, N

           RX2I = FX(I) / M
           RY2I = FY(I) / M
           RZ2I = FZ(I) / M
           CORRX = RX2I - RX2(I)
           CORRY = RY2I - RY2(I)
           CORRZ = RZ2I - RZ2(I)

           RX (I) = RX (I) + CTRAN0 * CORRX
           RY (I) = RY (I) + CTRAN0 * CORRY
           RZ (I) = RZ (I) + CTRAN0 * CORRZ
           RX1(I) = RX1(I) + CTRAN1 * CORRX
           RY1(I) = RY1(I) + CTRAN1 * CORRY
           RZ1(I) = RZ1(I) + CTRAN1 * CORRZ
           RX2(I) = RX2I
           RY2(I) = RY2I
           RZ2(I) = RZ2I
           RX3(I) = RX3(I) + CTRAN3 * CORRX
           RY3(I) = RY3(I) + CTRAN3 * CORRY
           RZ3(I) = RZ3(I) + CTRAN3 * CORRZ

           K = K + M * ( RX1(I) ** 2 + RY1(I) ** 2 + RZ1(I) ** 2 )

           QW1I = ( - QX(I)*OX(I) - QY(I)*OY(I) - QZ(I)*OZ(I) ) * 0.5
           QX1I = (   QW(I)*OX(I) - QZ(I)*OY(I) + QY(I)*OZ(I) ) * 0.5
           QY1I = (   QZ(I)*OX(I) + QW(I)*OY(I) - QX(I)*OZ(I) ) * 0.5
           QZ1I = ( - QY(I)*OX(I) + QX(I)*OY(I) + QW(I)*OZ(I) ) * 0.5

           CORRW = QW1I - QW1(I)
           CORRX = QX1I - QX1(I)
           CORRY = QY1I - QY1(I)
           CORRZ = QZ1I - QZ1(I)

           QW (I) = QW (I) + CROT0 * CORRW
           QX (I) = QX (I) + CROT0 * CORRX
           QY (I) = QY (I) + CROT0 * CORRY
           QZ (I) = QZ (I) + CROT0 * CORRZ
           QW1(I) = QW1I
           QX1(I) = QX1I
           QY1(I) = QY1I
           QZ1(I) = QZ1I
           QW2(I) = QW2(I) + CROT2 * CORRW
           QX2(I) = QX2(I) + CROT2 * CORRX
           QY2(I) = QY2(I) + CROT2 * CORRY
           QZ2(I) = QZ2(I) + CROT2 * CORRZ
           QW3(I) = QW3(I) + CROT3 * CORRW
           QX3(I) = QX3(I) + CROT3 * CORRX
           QY3(I) = QY3(I) + CROT3 * CORRY
           QZ3(I) = QZ3(I) + CROT3 * CORRZ
           QW4(I) = QW4(I) + CROT4 * CORRW
           QX4(I) = QX4(I) + CROT4 * CORRX
           QY4(I) = QY4(I) + CROT4 * CORRY
           QZ4(I) = QZ4(I) + CROT4 * CORRZ

           OX1I = ( TX(I) + OY(I) * OZ(I) * (IYY-IZZ) ) / IXX
           OY1I = ( TY(I) + OZ(I) * OX(I) * (IZZ-IXX) ) / IYY
           OZ1I = ( TZ(I) + OX(I) * OY(I) * (IXX-IYY) ) / IZZ

           CORRX = OX1I - OX1(I)
           CORRY = OY1I - OY1(I)
           CORRZ = OZ1I - OZ1(I)

           OX (I) = OX (I) + CROT0 * CORRX
           OY (I) = OY (I) + CROT0 * CORRY
           OZ (I) = OZ (I) + CROT0 * CORRZ
           OX1(I) = OX1I
           OY1(I) = OY1I
           OZ1(I) = OZ1I
           OX2(I) = OX2(I) + CROT2 * CORRX
           OY2(I) = OY2(I) + CROT2 * CORRY
           OZ2(I) = OZ2(I) + CROT2 * CORRZ
           OX3(I) = OX3(I) + CROT3 * CORRX
           OY3(I) = OY3(I) + CROT3 * CORRY
           OZ3(I) = OZ3(I) + CROT3 * CORRZ
           OX4(I) = OX4(I) + CROT4 * CORRX
           OY4(I) = OY4(I) + CROT4 * CORRY
           OZ4(I) = OZ4(I) + CROT4 * CORRZ

           K = K + IXX * OX(I) ** 2
     :           + IYY * OY(I) ** 2
     :           + IZZ * OZ(I) ** 2

400     CONTINUE

        K = 0.5 * K

        RETURN
        END



