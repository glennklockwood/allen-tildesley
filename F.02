********************************************************************************
** FICHE F.2.  GEAR 5-VALUE PREDICTOR-CORRECTOR ALGORITHM                     **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** GEAR 5-VALUE PREDICTOR-CORRECTOR ALGORITHM FOR TRANSLATION.   **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** GEAR, NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY            **
C    ** DIFFERENTIAL EQUATIONS (PRENTICE-HALL,1971).                  **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE PREDIC ( DT )                                      **
C    **    PREDICTS THE NEW POSITIONS, VELOCITIES, ETC.               **
C    ** SUBROUTINE CORREC ( DT, M, K )                                **
C    **    CORRECTS THE POSITIONS, VELOCITIES ETC. USING GEAR METHOD  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF MOLECULES                **
C    ** REAL    DT                 TIMESTEP                           **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    VX(N),VY(N),VZ(N)  VELOCITIES                         **
C    ** REAL    AX(N),AY(N),AZ(N)  ACCELERATIONS                      **
C    ** REAL    BX(N),BY(N),BZ(N)  THIRD DERIVATIVES                  **
C    ** REAL    CX(N),CY(N),CZ(N)  FOURTH DERIVATIVES                 **
C    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** AT EACH TIMESTEP, CALL PREDIC, FORCE, CORREC IN ORDER         **
C    ** FOLLOWED BY ACCUMULATION OF THERMODYNAMIC QUANTITIES.         **
C    ** THE FORCE ROUTINE (NOT SUPPLIED HERE: SEE F.17) CALCULATES    **
C    ** POTENTIAL ENERGY AND FORCES ON ALL ATOMS.                     **
C    *******************************************************************



        SUBROUTINE PREDIC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ

C    *******************************************************************
C    ** PREDICTOR ROUTINE                                             **
C    **                                                               **
C    ** IN TIMESTEP-SCALED VARIABLES THE PREDICTOR IS THE PASCAL      **
C    ** TRIANGLE MATRIX.  IN UNSCALED VARIABLES IT IS A TAYLOR SERIES **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** PREDIC IS CALLED TO ADVANCE THE COORDINATES, VELOCITIES ETC.  **
C    ** BY ONE TIMESTEP DT, PRIOR TO FORCE EVALUATION.                **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N)
        REAL        BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL        C1, C2, C3, C4

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        DO 100 I = 1, N

           RX(I) = RX(I) + C1*VX(I) + C2*AX(I) + C3*BX(I) + C4*CX(I)
           RY(I) = RY(I) + C1*VY(I) + C2*AY(I) + C3*BY(I) + C4*CY(I)
           RZ(I) = RZ(I) + C1*VZ(I) + C2*AZ(I) + C3*BZ(I) + C4*CZ(I)
           VX(I) = VX(I) + C1*AX(I) + C2*BX(I) + C3*CX(I)
           VY(I) = VY(I) + C1*AY(I) + C2*BY(I) + C3*CY(I)
           VZ(I) = VZ(I) + C1*AZ(I) + C2*BZ(I) + C3*CZ(I)
           AX(I) = AX(I) + C1*BX(I) + C2*CX(I)
           AY(I) = AY(I) + C1*BY(I) + C2*CY(I)
           AZ(I) = AZ(I) + C1*BZ(I) + C2*CZ(I)
           BX(I) = BX(I) + C1*CX(I)
           BY(I) = BY(I) + C1*CY(I)
           BZ(I) = BZ(I) + C1*CZ(I)

100     CONTINUE

        RETURN
        END



        SUBROUTINE CORREC ( DT, M, K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ

C    *******************************************************************
C    ** CORRECTOR ROUTINE                                             **
C    **                                                               **
C    ** CORRECTS POSITIONS, VELOCITIES ETC. AFTER FORCE EVALUATION.   **
C    ** IN TIMESTEP-SCALED VARIABLES THE NUMERICAL COEFFICIENTS ARE   **
C    ** GIVEN BY GEAR (REF ABOVE): 19/120, 3/4, 1, 1/2, 1/12.         **
C    ** IN UNSCALED FORM THESE MUST BE MULTIPLIED BY FACTORS          **
C    ** INVOLVING THE TIMESTEP AS SHOWN HERE.                         **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    M                        ATOMIC MASS                  **
C    ** REAL    K                        KINETIC ENERGY PER ATOM      **
C    ** REAL    GEAR0,GEAR1,GEAR2,GEAR3  GEAR COEFFICIENTS            **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** IT IS ASSUMED THAT INTERMOLECULAR FORCES HAVE BEEN CALCULATED **
C    ** AND STORED IN FX,FY,FZ. CORREC SIMPLY APPLIES THE CORRECTOR   **
C    ** EQUATIONS BASED ON THE DIFFERENCES BETWEEN PREDICTED AND      **
C    ** EVALUATED ACCELERATIONS.  IT ALSO CALCULATES KINETIC ENERGY.  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108)
        REAL        DT
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N)
        REAL        BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        M, K

        INTEGER     I
        REAL        AXI, AYI, AZI
        REAL        CORRX, CORRY, CORRZ
        REAL        C1, C2, C3, C4
        REAL        CR, CV, CB, CC
        REAL        GEAR0, GEAR1, GEAR3, GEAR4
        PARAMETER ( GEAR0 = 19.0 / 120.0, GEAR1 = 3.0 / 4.0,
     :              GEAR3 = 1.0 / 2.0,    GEAR4 = 1.0 / 12.0 )

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        CR = GEAR0 * C2
        CV = GEAR1 * C2 / C1
        CB = GEAR3 * C2 / C3
        CC = GEAR4 * C2 / C4

        K = 0.0

        DO 100 I = 1, N

           AXI = FX(I) / M
           AYI = FY(I) / M
           AZI = FZ(I) / M
           CORRX = AXI - AX(I)
           CORRY = AYI - AY(I)
           CORRZ = AZI - AZ(I)

           RX(I) = RX(I) + CR * CORRX
           RY(I) = RY(I) + CR * CORRY
           RZ(I) = RZ(I) + CR * CORRZ
           VX(I) = VX(I) + CV * CORRX
           VY(I) = VY(I) + CV * CORRY
           VZ(I) = VZ(I) + CV * CORRZ
           AX(I) = AXI
           AY(I) = AYI
           AZ(I) = AZI
           BX(I) = BX(I) + CB * CORRX
           BY(I) = BY(I) + CB * CORRY
           BZ(I) = BZ(I) + CB * CORRZ
           CX(I) = CX(I) + CC * CORRX
           CY(I) = CY(I) + CC * CORRY
           CZ(I) = CZ(I) + CC * CORRZ

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

100     CONTINUE

        K = 0.5 * M * K

        RETURN
        END



