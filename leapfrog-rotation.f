********************************************************************************
** FICHE F.6.  LEAPFROG ALGORITHMS FOR ROTATIONAL MOTION                      **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO SEPARATE PARTS: ROTATION OF LINEAR, NONLINEAR MOLECULES.  **
C    *******************************************************************



C    *******************************************************************
C    ** FICHE F.6 - PART A                                            **
C    ** LEAPFROG ALGORITHM FOR ROTATIONAL MOTION OF LINEAR MOLECULES. **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FINCHAM, CCP5 QUARTERLY 2, 6, 1981.                           **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVE ( DT, M, INERT, K )                           **
C    **    ADVANCES POSITIONS AND VELOCITIES                          **
C    ** SUBROUTINE MOLATM                                             **
C    **    CONVERTS MOLECULAR COORDINATES TO ATOMIC/SITE POSITIONS    **
C    ** SUBROUTINE ATMMOL                                             **
C    **    CONVERTS ATOMIC FORCES TO MOLECULAR FORCES AND "TORQUES"   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    DT                            TIMESTEP                **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    M                             MOLECULAR MASS          **
C    ** REAL    INERT                         MOMENT OF INERTIA       **
C    ** REAL    K                             KINETIC ENERGY          **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    VX(N),VY(N),VZ(N)             VELOCITIES AT TIME T    **
C    ** REAL    FX(N),FY(N),FZ(N)             C-O-M FORCES            **
C    ** REAL    EX(N),EY(N),EZ(N)             UNIT BOND VEC AT TIME T **
C    ** REAL    UX(N),UY(N),UZ(N)             TIME DERIV AT T-DT/2    **
C    ** REAL    GX(N),GY(N),GZ(N)             AUXILIARY TORQUE AT T   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SUBROUTINE MOLATM IS CALLED, TO OBTAIN ATOMIC SITE POSITIONS  **
C    ** WHICH ARE USED BY THE FORCE ROUTINE (NOT SUPPLIED HERE) TO    **
C    ** CALCULATE ATOMIC FORCES.  SUBROUTINE ATMMOL THEN CONVERTS     **
C    ** THESE INTO MOLECULAR FORCE AND MODIFIED TORQUE TERMS.         **
C    ** SUBROUTINE MOVE THEN ADVANCES THE POSITIONS ETC.              **
C    ** FOR THIS EXAMPLE WE TAKE A (LINEAR) TRIATOMIC MOLECULE.       **
C    *******************************************************************



        SUBROUTINE MOVE ( DT, M, INERT, K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / EX, EY, EZ, UX, UY, UZ, GX, GY, GZ

C    *******************************************************************
C    ** ADVANCES POSITIONS, BOND VECTORS, AND THEIR TIME DERIVATIVES. **
C    **                                                               **
C    ** THIS METHOD USES AN AUXILIARY VECTOR TO DESCRIBE THE TORQUE   **
C    ** AND THE BOND VECTOR DERIVATIVE INSTEAD OF ANGULAR VELOCITY.   **
C    ** EVERYTHING IS IN SPACE-FIXED AXES.                            **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        DT
        REAL        M
        REAL        INERT, K
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        EX(N), EY(N), EZ(N)
        REAL        UX(N), UY(N), UZ(N)
        REAL        GX(N), GY(N), GZ(N)

        INTEGER     I
        REAL        UXI, UYI, UZI, EXI, EYI, EZI, VXI, VYI, VZI, DOT

C    *******************************************************************

        K = 0.0

        DO 400 I = 1, N

C       ** MOVE BOND VECTOR DERIVATIVES         **
C       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY **

           UXI = UX(I)
           UYI = UY(I)
           UZI = UZ(I)
           EXI = EX(I)
           EYI = EY(I)
           EZI = EZ(I)
           DOT = 2.0 * ( UXI * EXI + UYI * EYI + UZI * EZI )
           UX(I) = UXI + DT * GX(I) / INERT - DOT * EXI
           UY(I) = UYI + DT * GY(I) / INERT - DOT * EYI
           UZ(I) = UZI + DT * GZ(I) / INERT - DOT * EZI
           UXI = 0.5 * ( UXI + UX(I) )
           UYI = 0.5 * ( UYI + UY(I) )
           UZI = 0.5 * ( UZI + UZ(I) )
           K = K + INERT * ( UXI ** 2 + UYI ** 2 + UZI ** 2 )

C       ** ADVANCE BOND VECTORS TO T+DT **

           EX(I) = EXI + DT * UX(I)
           EY(I) = EYI + DT * UY(I)
           EZ(I) = EZI + DT * UZ(I)

C       ** MOVE THE LINEAR VELOCITIES ALL THE WAY **
C       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY   **

           VXI = VX(I)
           VYI = VY(I)
           VZI = VZ(I)
           VX(I) = VXI + DT * FX(I) / M
           VY(I) = VYI + DT * FY(I) / M
           VZ(I) = VZI + DT * FZ(I) / M
           VXI = 0.5 * ( VXI + VX(I) )
           VYI = 0.5 * ( VYI + VY(I) )
           VZI = 0.5 * ( VZI + VZ(I) )
           K = K + M * ( VXI **2 + VYI ** 2 + VZI ** 2 )

C       ** ADVANCE POSITIONS TO T+DT **

           RX(I) = RX(I) + DT * VX(I)
           RY(I) = RY(I) + DT * VY(I)
           RZ(I) = RZ(I) + DT * VZ(I)

400     CONTINUE

        K = 0.5 * K

        RETURN
        END



        SUBROUTINE MOLATM

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / EX, EY, EZ, UX, UY, UZ, GX, GY, GZ
        COMMON / BLOCK3 / D, RSX, RSY, RSZ, FSX, FSY, FSZ

C    *******************************************************************
C    ** CONVERTS C-O-M COORDINATES AND BOND VECTOR TO SITE POSITIONS. **
C    **                                                               **
C    ** THE POSITION OF EACH ATOM IN THE MOLECULE IS DEFINED IN TERMS **
C    ** OF THE UNIT BOND VECTOR EX(I),EY(I),EZ(I) AND THE ATOM        **
C    ** POSITION VARIABLE D(A): RSX(I,A) = RX(I) + D(A)*EX(I) ETC.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    EX(N),EY(N),EZ(N)             UNIT BOND VEC AT TIME T **
C    ** REAL    D(NA)                         ATOM POSITIONS IN MOLEC **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        EX(N), EY(N), EZ(N)
        REAL        UX(N), UY(N), UZ(N)
        REAL        GX(N), GY(N), GZ(N)
        REAL        D(NA)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)

        INTEGER     I, A
        REAL        EXI, EYI, EZI

C    *******************************************************************

        DO 200 I = 1, N

           EXI = EX(I)
           EYI = EY(I)
           EZI = EZ(I)

           DO 199 A = 1, NA

              RSX(I,A) = RX(I) + D(A) * EXI
              RSY(I,A) = RY(I) + D(A) * EYI
              RSZ(I,A) = RZ(I) + D(A) * EZI

199        CONTINUE

200     CONTINUE

        RETURN
        END



        SUBROUTINE ATMMOL

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / EX, EY, EZ, UX, UY, UZ, GX, GY, GZ
        COMMON / BLOCK3 / D, RSX, RSY, RSZ, FSX, FSY, FSZ

C    *******************************************************************
C    ** CONVERT ATOM FORCES TO TOTAL FORCES AND AUXILIARY TORQUES.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    FX(N),FY(N),FZ(N)             C-O-M FORCES            **
C    ** REAL    EX(N),EY(N),EZ(N)             UNIT BOND VEC AT TIME T **
C    ** REAL    GX(N),GY(N),GZ(N)             AUXILIARY TORQUE AT T   **
C    ** REAL    D(NA)                         ATOM POSITIONS IN MOLEC **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    ** REAL    FSX(N,NA),FSY(N,NA),FSZ(N,NA) FORCES ON EACH ATOM     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        EX(N), EY(N), EZ(N)
        REAL        UX(N), UY(N), UZ(N)
        REAL        GX(N), GY(N), GZ(N)
        REAL        D(NA)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)

        INTEGER     I, A
        REAL        FXI, FYI, FZI, GXI, GYI, GZI
        REAL        RXI, RYI, RZI, EXI, EYI, EZI
        REAL        FSXIA, FSYIA, FSZIA, DOT

C    *******************************************************************

        DO 300 I = 1, N

           FXI = 0.0
           FYI = 0.0
           FZI = 0.0
           GXI = 0.0
           GYI = 0.0
           GZI = 0.0
           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           EXI = EX(I)
           EYI = EY(I)
           EZI = EZ(I)

           DO 299 A = 1, NA

              FSXIA = FSX(I,A)
              FSYIA = FSY(I,A)
              FSZIA = FSZ(I,A)
              FXI = FXI + FSXIA
              FYI = FYI + FSYIA
              FZI = FZI + FSZIA
              GXI = GXI + D(A) * FSXIA
              GYI = GYI + D(A) * FSYIA
              GZI = GZI + D(A) * FSZIA

299        CONTINUE

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI
           DOT   = GXI * EXI + GYI * EYI + GZI * EZI
           GX(I) = GXI - DOT * EXI
           GY(I) = GYI - DOT * EYI
           GZ(I) = GZI - DOT * EZI

300     CONTINUE

        RETURN
        END



C    *******************************************************************
C    ** FICHE F.6 - PART B                                            **
C    ** LEAPFROG ALGORITHM FOR ROTATIONAL MOTION, NONLINEAR MOLECULES.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FINCHAM, CCP5 QUARTERLY 12, 47, 1984.                         **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVE ( DT, M, IXX, IYY, IZZ, K )                   **
C    **    ADVANCES POSITIONS, ORIENTATIONS, AND TIME DERIVATIVES     **
C    ** SUBROUTINE MOLATM                                             **
C    **    CONVERTS MOLECULAR COORDINATES INTO ATOMIC SITE POSITIONS  **
C    ** SUBROUTINE ATMMOL                                             **
C    **    CONVERTS ATOMIC FORCES INTO MOLECULAR FORCES AND TORQUES   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    DT                            TIMESTEP                **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** REAL    M                             MOLECULAR MASS          **
C    ** REAL    IXX,IYY,IZZ                   PRINCIPAL INERTIAS      **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    VX(N),VY(N),VZ(N)             VELOCITIES AT TIME T    **
C    ** REAL    FX(N),FY(N),FZ(N)             C-O-M FORCES            **
C    ** REAL    QW(N),QX(N),QY(N),QZ(N)       QUATERNIONS AT TIME T   **
C    ** REAL    JX(N),JY(N),JZ(N)             ANGULAR MOM. AT T-DT/2  **
C    ** REAL    TX(N),TY(N),TZ(N)             TORQUE AT T             **
C    ** REAL    K                             KINETIC ENERGY          **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** WE USE QUATERNION PARAMETERS FOR THE ORIENTATION.             **
C    ** THIS METHOD USES AN AUXILIARY EQUATION TO OBTAIN ACCURATE     **
C    ** QUATERNIONS AND ROTATION MATRICES AT THE HALF-STEP TIME.      **
C    ** ANGULAR MOMENTUM AND TORQUE ARE IN SPACE-FIXED AXES.          **
C    ** WE ASSUME THAT WE ARE ALSO USING LEAPFROG FOR TRANSLATION     **
C    ** SUBROUTINE MOLATM IS CALLED, FOLLOWED BY THE FORCE ROUTINE    **
C    ** (NOT SUPPLIED HERE).  AFTER THIS, SUBROUTINE ATMMOL IS CALLED **
C    ** AND THEN SUBROUTINE MOVE ADVANCES THE CONFIGURATION.          **
C    ** FOR THIS EXAMPLE WE TAKE A (NONLINEAR) TRIATOMIC MOLECULE.    **
C    *******************************************************************



        SUBROUTINE MOVE ( DT, M, IXX, IYY, IZZ, K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / QW, QX, QY, QZ, JX, JY, JZ, TX, TY, TZ

C    *******************************************************************
C    ** ADVANCE THE CONFIGURATION AND CALCULATE KINETIC ENERGY        **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT
        REAL        M
        REAL        IXX, IYY, IZZ, K
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        QW(N), QX(N), QY(N), QZ(N)
        REAL        JX(N), JY(N), JZ(N)
        REAL        TX(N), TY(N), TZ(N)

        INTEGER     I
        REAL        DT2
        REAL        JXI, JYI, JZI, OXI, OYI, OZI, QWI, QXI, QYI, QZI
        REAL        QW1I, QX1I, QY1I, QZ1I, VXI, VYI, VZI
        REAL        AXX, AXY, AXZ, AYX, AYY, AYZ, AZX, AZY, AZZ

C    *******************************************************************

        K   = 0.0
        DT2 = DT / 2.0

        DO 400 I = 1, N

C       ** AUXILIARY EQUATION MOVES   **
C       ** ANGULAR MOMENTUM TO TIME T **

           JXI = JX(I) + DT2 * TX(I)
           JYI = JY(I) + DT2 * TY(I)
           JZI = JZ(I) + DT2 * TZ(I)

C       ** OBTAIN ROTATION MATRIX AT TIME T **

           AXX = QW(I) ** 2 + QX(I) ** 2 - QY(I) ** 2 - QZ(I) ** 2
           AXY = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
           AXZ = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
           AYX = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
           AYY = QW(I) ** 2 - QX(I) ** 2 + QY(I) ** 2 - QZ(I) ** 2
           AYZ = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
           AZX = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
           AZY = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
           AZZ = QW(I) ** 2 - QX(I) ** 2 - QY(I) ** 2 + QZ(I) ** 2

C       ** CONVERT ANGULAR MOMENTUM TO BODY-FIXED **
C       ** FORM AND HENCE TO ANGULAR VELOCITIES   **

           OXI = ( AXX * JXI + AXY * JYI + AXZ * JZI ) / IXX
           OYI = ( AYX * JXI + AYY * JYI + AYZ * JZI ) / IYY
           OZI = ( AZX * JXI + AZY * JYI + AZZ * JZI ) / IZZ

           K = K + IXX * OXI ** 2 + IYY * OYI ** 2 + IZZ * OZI ** 2

C       ** OBTAIN TIME-DERIVATIVES OF QUATERNIONS **
C       ** AND ADVANCE TO TIME T+DT/2             **

           QW1I = ( - QX(I) * OXI - QY(I) * OYI - QZ(I) * OZI ) * 0.5
           QX1I = (   QW(I) * OXI - QZ(I) * OYI + QY(I) * OZI ) * 0.5
           QY1I = (   QZ(I) * OXI + QW(I) * OYI - QX(I) * OZI ) * 0.5
           QZ1I = ( - QY(I) * OXI + QX(I) * OYI + QW(I) * OZI ) * 0.5
           QWI = QW(I) + DT2 * QW1I
           QXI = QX(I) + DT2 * QX1I
           QYI = QY(I) + DT2 * QY1I
           QZI = QZ(I) + DT2 * QZ1I

C       ** OBTAIN ROTATION MATRIX AT TIME T+DT/2 **

           AXX = QWI ** 2 + QXI ** 2 - QYI ** 2 - QZI ** 2
           AXY = 2.0 * ( QXI * QYI + QWI * QZI )
           AXZ = 2.0 * ( QXI * QZI - QWI * QYI )
           AYX = 2.0 * ( QXI * QYI - QWI * QZI )
           AYY = QWI ** 2 - QXI ** 2 + QYI ** 2 - QZI ** 2
           AYZ = 2.0 * ( QYI * QZI + QWI * QXI )
           AZX = 2.0 * ( QXI * QZI + QWI * QYI )
           AZY = 2.0 * ( QYI * QZI - QWI * QXI )
           AZZ = QWI ** 2 - QXI ** 2 - QYI ** 2 + QZI ** 2

C       ** MOVE THE ANGULAR MOMENTA ALL THE WAY     **
C       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY     **
C       ** CONVERT TO BODY-FIXED ANGULAR VELOCITIES **
C       ** AT TIME T+DT/2                           **

           JX(I) = JX(I) + DT * TX(I)
           JY(I) = JY(I) + DT * TY(I)
           JZ(I) = JZ(I) + DT * TZ(I)

           OXI = ( AXX * JX(I) + AXY * JY(I) + AXZ * JZ(I) ) / IXX
           OYI = ( AYX * JX(I) + AYY * JY(I) + AYZ * JZ(I) ) / IYY
           OZI = ( AZX * JX(I) + AZY * JY(I) + AZZ * JZ(I) ) / IZZ

C       ** OBTAIN TIME-DERIVATIVES OF QUATERNIONS **
C       ** AND ADVANCE TO T+DT                    **

           QW1I = ( - QXI * OXI - QYI * OYI - QZI * OZI ) * 0.5
           QX1I = (   QWI * OXI - QZI * OYI + QYI * OZI ) * 0.5
           QY1I = (   QZI * OXI + QWI * OYI - QXI * OZI ) * 0.5
           QZ1I = ( - QYI * OXI + QXI * OYI + QWI * OZI ) * 0.5
           QW(I) = QW(I) + DT * QW1I
           QX(I) = QX(I) + DT * QX1I
           QY(I) = QY(I) + DT * QY1I
           QZ(I) = QZ(I) + DT * QZ1I

C       ** MOVE THE LINEAR VELOCITIES ALL THE WAY   **
C       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY     **

           VXI = VX(I)
           VYI = VY(I)
           VZI = VZ(I)

           VX(I) = VXI + DT * FX(I) / M
           VY(I) = VYI + DT * FY(I) / M
           VZ(I) = VZI + DT * FZ(I) / M

           VXI = 0.5 * ( VXI + VX(I) )
           VYI = 0.5 * ( VYI + VY(I) )
           VZI = 0.5 * ( VZI + VZ(I) )

           K = K + M * ( VXI ** 2 + VYI ** 2 + VZI ** 2 )

C       ** ADVANCE POSITIONS TO T+DT **

           RX(I) = RX(I) + DT * VX(I)
           RY(I) = RY(I) + DT * VY(I)
           RZ(I) = RZ(I) + DT * VZ(I)

400     CONTINUE

        K = 0.5 * K

        RETURN
        END



        SUBROUTINE MOLATM

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / QW, QX, QY, QZ, JX, JY, JZ, TX, TY, TZ
        COMMON / BLOCK3 / DX, DY, DZ, RSX, RSY, RSZ, FSX, FSY, FSZ

C    *******************************************************************
C    ** COMPUTE ELEMENTS OF ROTATION MATRIX FOR EACH MOLECULE I.      **
C    **                                                               **
C    ** THE TRANSPOSE OF THE ROTATION MATRIX IS USED TO OBTAIN THE    **
C    ** POSITIONS OF EACH ATOM FROM THE CENTRE-OF-MASS POSITION AND   **
C    ** THE BODY-FIXED ATOM POSITION VECTORS (KNOWN FROM THE START).  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    QW(N),QX(N),QY(N),QZ(N)       QUATERNIONS AT TIME T   **
C    ** REAL    DX(NA),DY(NA),DZ(NA)          ATOM POSITIONS IN MOLEC **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        QW(N), QX(N), QY(N), QZ(N)
        REAL        JX(N), JY(N), JZ(N)
        REAL        TX(N), TY(N), TZ(N)
        REAL        DX(NA), DY(NA), DZ(NA)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)

        INTEGER     I, A
        REAL        AXX, AXY, AXZ, AYX, AYY, AYZ, AZX, AZY, AZZ

C    *******************************************************************

        DO 200 I = 1, N

           AXX = QW(I) ** 2 + QX(I) ** 2 - QY(I) ** 2 - QZ(I) ** 2
           AXY = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
           AXZ = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
           AYX = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
           AYY = QW(I) ** 2 - QX(I) ** 2 + QY(I) ** 2 - QZ(I) ** 2
           AYZ = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
           AZX = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
           AZY = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
           AZZ = QW(I) ** 2 - QX(I) ** 2 - QY(I) ** 2 + QZ(I) ** 2

           DO 199 A = 1, NA

              RSX(I,A) = RX(I) + AXX * DX(A) + AYX * DY(A) + AZX * DZ(A)
              RSY(I,A) = RY(I) + AXY * DX(A) + AYY * DY(A) + AZY * DZ(A)
              RSZ(I,A) = RZ(I) + AXZ * DX(A) + AYZ * DY(A) + AZZ * DZ(A)

199        CONTINUE

200     CONTINUE

        RETURN
        END



        SUBROUTINE ATMMOL

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / QW, QX, QY, QZ, JX, JY, JZ, TX, TY, TZ
        COMMON / BLOCK3 / DX, DY, DZ, RSX, RSY, RSZ, FSX, FSY, FSZ

C    *******************************************************************
C    ** CONVERT ATOM FORCES TO TOTAL FORCES AND TORQUES               **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C    ** REAL    FX(N),FY(N),FZ(N)             C-O-M FORCES            **
C    ** REAL    QW(N),QX(N),QY(N),QZ(N)       QUATERNIONS AT TIME T   **
C    ** REAL    TX(N),TY(N),TZ(N)             TORQUE AT T             **
C    ** REAL    DX(NA),DY(NA),DZ(NA)          ATOM POSITIONS IN MOLEC **
C    ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C    ** REAL    FSX(N,NA),FSY(N,NA),FSZ(N,NA) FORCES ON EACH ATOM     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        INTEGER     NA
        PARAMETER ( NA = 3 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        QW(N), QX(N), QY(N), QZ(N)
        REAL        JX(N), JY(N), JZ(N)
        REAL        TX(N), TY(N), TZ(N)
        REAL        DX(NA), DY(NA), DZ(NA)
        REAL        RSX(N,NA), RSY(N,NA), RSZ(N,NA)
        REAL        FSX(N,NA), FSY(N,NA), FSZ(N,NA)

        INTEGER     I, A
        REAL        RXI, RYI, RZI, FXI, FYI, FZI, TXI, TYI, TZI
        REAL        FSXIA, FSYIA, FSZIA, RSXIA, RSYIA, RSZIA

C    *******************************************************************

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

           DO 299 A = 1, NA

              FSXIA = FSX(I,A)
              FSYIA = FSY(I,A)
              FSZIA = FSZ(I,A)
              RSXIA = RSX(I,A) - RXI
              RSYIA = RSY(I,A) - RYI
              RSZIA = RSZ(I,A) - RZI
              FXI = FXI + FSXIA
              FYI = FYI + FSYIA
              FZI = FZI + FSZIA
              TXI = TXI + RSYIA * FSZIA - RSZIA * FSYIA
              TYI = TYI + RSZIA * FSXIA - RSXIA * FSZIA
              TZI = TZI + RSXIA * FSYIA - RSYIA * FSXIA

299        CONTINUE

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI
           TX(I) = TXI
           TY(I) = TYI
           TZ(I) = TZI

300     CONTINUE

        RETURN
        END



