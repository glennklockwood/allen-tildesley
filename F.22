********************************************************************************
** FICHE F.22.  ROUTINES TO PERFORM THE EWALD SUM                             **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** REAL-SPACE AND RECIPROCAL-SPACE PARTS OF EWALD SUM FOR IONS.  **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** WOODCOCK AND SINGER, TRANS. FARADAY SOC. 67, 12, 1971.        **
C    ** DE LEEUW ET AL., PROC. ROY. SOC. A 373, 27, 1980.             **
C    ** HEYES, J. CHEM. PHYS. 74, 1924, 1981.                         **
C    ** SEE ALSO FINCHAM, MDIONS, CCP5 PROGRAM LIBRARY.               **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE SETUP ( KAPPA )                                    **
C    **    SETS UP THE WAVEVECTORS FOR USE IN THE EWALD SUM           **
C    ** SUBROUTINE RWALD ( KAPPA, VR )                                **
C    **    CALCULATES THE R-SPACE PART OF THE SUM                     **
C    ** SUBROUTINE KWALD ( KAPPA, VK )                                **
C    **    CALCULATES THE K-SPACE PART OF THE SUM                     **
C    ** REAL FUNCTION ERFC ( X )                                      **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER  TOTK         THE TOTAL NUMBER OF K-VECTORS STORED    **
C    ** INTEGER  MAXK         MAXIMUM POSSIBLE NUMBER OF K-VECTORS    **
C    ** INTEGER  KMAX         MAX INTEGER COMPONENT OF THE K-VECTOR   **
C    ** INTEGER  KSQMAX       MAX SQUARE MOD OF THE K-VECTOR REQUIRED **
C    ** REAL     VR           ENERGY FROM R-SPACE SUM                 **
C    ** REAL     VK           ENERGY FROM K-SPACE SUM                 **
C    ** REAL     KVEC(MAXK)   ARRAY USED TO STORE K-VECTORS           **
C    ** REAL     KAPPA        WIDTH OF CANCELLING DISTRIBUTION        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SETUP IS CALLED ONCE AT THE BEGINNING OF THE SIMULATION       **
C    ** TO CALCULATE ALL THE K-VECTORS REQUIRED IN THE EWALD SUM.     **
C    ** THESE VECTORS ARE USED THROUGHOUT THE SIMULATION IN THE       **
C    ** SUBROUTINE KWALD TO CALCULATE THE K-SPACE CONTRIBUTION TO THE **
C    ** POTENTIAL ENERGY AT EACH CONFIGURATION. THE SELF TERM IS      **
C    ** SUBTRACTED FROM THE K-SPACE CONTRIBUTION IN KWALD.            **
C    ** THE SURFACE TERM FOR SIMULATIONS IN VACUUM IS NOT INCLUDED.   **
C    ** ROUTINE RWALD RETURNS THE R-SPACE CONTRIBUTION TO THE EWALD   **
C    ** SUM AND IS CALLED FOR EACH CONFIGURATION IN THE SIMULATION.   **
C    ** A CUBIC BOX AND UNIT BOX LENGTH ARE ASSUMED THROUGHOUT.       **
C    *******************************************************************



        SUBROUTINE SETUP ( KAPPA )

        COMMON / BLOCK2 / KVEC

C    *******************************************************************
C    ** ROUTINE TO SET UP THE WAVE-VECTORS FOR THE EWALD SUM.         **
C    **                                                               **
C    ** THE WAVEVECTORS MUST FIT INTO A BOX OF UNIT LENGTH.           **
C    ** IN THIS EXAMPLE WE ALLOW A MAXIMUM OF 1000 WAVEVECTORS.       **
C    *******************************************************************

        INTEGER     MAXK
        PARAMETER ( MAXK = 1000 )

        REAL        KVEC(MAXK), KAPPA

        INTEGER     KMAX, KSQMAX, KSQ, KX, KY, KZ, TOTK
        REAL        TWOPI, B, RKX, RKY, RKZ, RKSQ
        PARAMETER ( KMAX = 5, KSQMAX = 27 , TWOPI = 6.2831853 )

C    *******************************************************************

        B = 1.0 / 4.0 / KAPPA / KAPPA

C    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **

        TOTK = 0

        DO 100 KX = 0, KMAX

           RKX = TWOPI * REAL ( KX )

           DO 99 KY = -KMAX, KMAX

              RKY = TWOPI * REAL ( KY )

              DO 98 KZ = -KMAX, KMAX

                 RKZ = TWOPI * REAL ( KZ )

                 KSQ = KX * KX + KY * KY + KZ * KZ

                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN

                    TOTK = TOTK + 1

                    IF ( TOTK .GT. MAXK ) STOP 'KVEC IS TOO SMALL'

                    RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
                    KVEC(TOTK) = TWOPI * EXP ( -B * RKSQ ) / RKSQ

                 ENDIF

98            CONTINUE

99         CONTINUE

100     CONTINUE

        WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
        WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) TOTK

        RETURN
        END



        SUBROUTINE RWALD ( KAPPA, VR )

        COMMON / BLOCK1 / RX, RY, RZ, Z

C    *******************************************************************
C    ** CALCULATES R-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                     NUMBER OF IONS                  **
C    ** REAL    RX(N),RY(N),RZ(N)     POSITIONS OF IONS               **
C    ** REAL    Z(N)                  IONIC CHARGES                   **
C    ** REAL    VR                    R-SPACE POTENTIAL ENERGY        **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION ERFC ( X )                                      **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 216 )

        REAL        RX(N), RY(N), RZ(N), Z(N)
        REAL        KAPPA, VR

        REAL        RXI, RYI, RZI, ZI, RXIJ, RYIJ, RZIJ
        REAL        RIJSQ, RIJ, KRIJ, ERFC, VIJ

        INTEGER     I, J

C    *******************************************************************

        VR = 0.0

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           ZI  = Z(I)

           DO 99 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RZIJ = RZIJ - ANINT ( RZIJ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
              RIJ   = SQRT ( RIJSQ )
              KRIJ  = KAPPA * RIJ
              VIJ   = ZI * Z(J) * ERFC ( KRIJ ) / RIJ

              VR    = VR + VIJ

99         CONTINUE

100     CONTINUE

        RETURN
        END



        SUBROUTINE KWALD ( KAPPA, VK )

        COMMON / BLOCK1 / RX, RY, RZ, Z
        COMMON / BLOCK2 / KVEC

C    *******************************************************************
C    ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** THE SELF TERM IS SUBTRACTED.                                  **
C    ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
C    ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
C    ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
C    ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF IONS                    **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS OF IONS                 **
C    ** REAL    Z(N)                IONIC CHARGES                     **
C    ** REAL    VK                  K-SPACE POTENTIAL ENERGY          **
C    ** REAL    VKS                 SELF PART OF K-SPACE SUM          **
C    *******************************************************************

        INTEGER     MAXK, N
        PARAMETER ( MAXK = 1000, N = 216 )

        REAL        KVEC(MAXK), RX(N), RY(N), RZ(N), Z(N)
        REAL        KAPPA, VK
        INTEGER     TOTK

        INTEGER     KMAX, KX, KY, KZ, I, KSQMAX, KSQ
        REAL        TWOPI, FACTOR, VD, VS, RSQPI
        PARAMETER ( KMAX = 5, KSQMAX = 27 )
        PARAMETER ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX     EIKX(1:N, 0:KMAX)
        COMPLEX     EIKY(1:N, -KMAX:KMAX)
        COMPLEX     EIKZ(1:N, -KMAX:KMAX)
        COMPLEX     EIKR(N), SUM

C    *******************************************************************

C    ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **

C    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **

        DO 10 I = 1, N

           EIKX(I, 0) = (1.0, 0.0)
           EIKY(I, 0) = (1.0, 0.0)
           EIKZ(I, 0) = (1.0, 0.0)

           EIKX(I, 1) = CMPLX ( COS ( TWOPI * RX(I) ) ,
     :                          SIN ( TWOPI * RX(I) ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPI * RY(I) ) ,
     :                          SIN ( TWOPI * RY(I) ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPI * RZ(I) ) ,
     :                          SIN ( TWOPI * RZ(I) ) )

           EIKY(I, -1) = CONJG ( EIKY(I, 1) )
           EIKZ(I, -1) = CONJG ( EIKZ(I, 1) )

10      CONTINUE

C    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **

        DO 12 KX = 2, KMAX

           DO 11 I = 1, N

              EIKX(I, KX) = EIKX(I, KX-1) * EIKX(I, 1)

11         CONTINUE

12      CONTINUE

        DO 14 KY = 2, KMAX

           DO 13 I = 1, N

              EIKY(I,  KY) = EIKY(I, KY-1) * EIKY(I, 1)
              EIKY(I, -KY) = CONJG ( EIKY(I, KY) )

13         CONTINUE

14      CONTINUE

        DO 16 KZ = 2, KMAX

           DO 15 I = 1, N

              EIKZ(I,  KZ) = EIKZ(I, KZ-1) * EIKZ(I, 1)
              EIKZ(I, -KZ) = CONJG ( EIKZ(I, KZ) )

15         CONTINUE

16      CONTINUE

C    ** SUM OVER ALL VECTORS **

        VD   = 0.0
        TOTK = 0

        DO 24 KX = 0, KMAX

           IF ( KX .EQ. 0 ) THEN

              FACTOR = 1.0

           ELSE

              FACTOR = 2.0

           ENDIF

           DO 23 KY = -KMAX, KMAX

              DO 22 KZ = -KMAX, KMAX

                 KSQ = KX * KX + KY * KY + KZ * KZ

                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN

                    TOTK = TOTK + 1
                    SUM  = (0.0, 0.0)

                    DO 21 I = 1, N

                       EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ)
                       SUM     = SUM + Z(I) * EIKR(I)

21                  CONTINUE

                    VD = VD + FACTOR * KVEC(TOTK) * CONJG ( SUM ) * SUM

                 ENDIF

22            CONTINUE

23         CONTINUE

24      CONTINUE

C    ** CALCULATES SELF PART OF K-SPACE SUM **

        VS = 0.0

        DO 25 I = 1, N

           VS = VS + Z(I) * Z(I)

25      CONTINUE

        VS = RSQPI * KAPPA * VS

C    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **

        VK = VD - VS

        RETURN
        END



        REAL FUNCTION ERFC ( X )

C    *******************************************************************
C    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
C    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
C    *******************************************************************

        REAL        A1, A2, A3, A4, A5, P

        PARAMETER ( A1 = 0.254829592, A2 = -0.284496736 )
        PARAMETER ( A3 = 1.421413741, A4 = -1.453152027 )
        PARAMETER ( A5 = 1.061405429, P  =  0.3275911   )

        REAL        T, X, XSQ, TP

C    *******************************************************************

        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X

        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

        ERFC = TP * EXP ( -XSQ )

        RETURN
        END



