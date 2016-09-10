********************************************************************************
** FICHE F.37.  ROUTINES TO CALCULATE FOURIER TRANSFORMS.                     **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** THREE SEPARATE ROUTINES FOR DIFFERENT APPLICATIONS.           **
C    *******************************************************************



        SUBROUTINE FILONC ( DT, DOM, NMAX, C, CHAT )

C    *******************************************************************
C    ** CALCULATES THE FOURIER COSINE TRANSFORM BY FILON'S METHOD     **
C    **                                                               **
C    ** A CORRELATION FUNCTION, C(T), IN THE TIME DOMAIN, IS          **
C    ** TRANSFORMED TO A SPECTRUM CHAT(OMEGA) IN THE FREQUENCY DOMAIN.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FILON, PROC ROY SOC EDIN, 49 38, 1928.                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    C(NMAX)            THE CORRELATION FUNCTION.          **
C    ** REAL    CHAT(NMAX)         THE 1-D COSINE TRANSFORM.          **
C    ** REAL    DT                 TIME INTERVAL BETWEEN POINTS IN C. **
C    ** REAL    DOM                FREQUENCY INTERVAL FOR CHAT.       **
C    ** INTEGER NMAX               NO. OF INTERVALS ON THE TIME AXIS  **
C    ** REAL    OMEGA              THE FREQUENCY                      **
C    ** REAL    TMAX               MAXIMUM TIME IN CORRL. FUNCTION    **
C    ** REAL    ALPHA, BETA, GAMMA FILON PARAMETERS                   **
C    ** INTEGER TAU                TIME INDEX                         **
C    ** INTEGER NU                 FREQUENCY INDEX                    **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE REQUIRES THAT THE NUMBER OF INTERVALS, NMAX, IS   **
C    ** EVEN AND CHECKS FOR THIS CONDITION. THE FIRST VALUE OF C(T)   **
C    ** IS AT T=0. THE MAXIMUM TIME FOR THE CORRELATION FUNCTION IS   **
C    ** TMAX=DT*NMAX. FOR AN ACCURATE TRANSFORM C(TMAX)=0.            **
C    *******************************************************************

        INTEGER    NMAX
        REAL       DT, DOM, C(0:NMAX), CHAT(0:NMAX)

        REAL       TMAX, OMEGA, THETA, SINTH, COSTH, CE, CO
        REAL       SINSQ, COSSQ, THSQ, THCUB, ALPHA, BETA, GAMMA
        INTEGER    TAU, NU

C    *******************************************************************

C    ** CHECKS NMAX IS EVEN **

        IF ( MOD ( NMAX, 2 ) .NE. 0 ) THEN

           STOP ' NMAX SHOULD BE EVEN '

        ENDIF

        TMAX = REAL ( NMAX ) * DT

C    ** LOOP OVER OMEGA **

        DO 30 NU = 0, NMAX

           OMEGA = REAL ( NU ) * DOM
           THETA = OMEGA * DT

C       ** CALCULATE THE FILON PARAMETERS **

           SINTH = SIN ( THETA )
           COSTH = COS ( THETA )
           SINSQ = SINTH * SINTH
           COSSQ = COSTH * COSTH
           THSQ  = THETA * THETA
           THCUB = THSQ * THETA

           IF ( THETA. EQ. 0.0 ) THEN

              ALPHA = 0.0
              BETA  = 2.0 / 3.0
              GAMMA = 4.0 / 3.0

            ELSE

              ALPHA = ( 1.0 / THCUB )
     :                * ( THSQ + THETA * SINTH * COSTH - 2.0 * SINSQ )
              BETA  = ( 2.0 / THCUB )
     :                * ( THETA * ( 1.0 + COSSQ ) -2.0 * SINTH * COSTH )
              GAMMA = ( 4.0 / THCUB ) * ( SINTH - THETA * COSTH )

           ENDIF

C       ** DO THE SUM OVER THE EVEN ORDINATES **

           CE = 0.0

           DO 10 TAU = 0, NMAX, 2

              CE = CE + C(TAU) * COS ( THETA * REAL ( TAU ) )

10         CONTINUE

C       ** SUBTRACT HALF THE FIRST AND LAST TERMS **

           CE = CE - 0.5 * ( C(0) + C(NMAX) * COS ( OMEGA * TMAX ) )

C       ** DO THE SUM OVER THE ODD ORDINATES **

           CO = 0.0

           DO 20 TAU = 1, NMAX - 1, 2

              CO = CO + C(TAU) * COS ( THETA * REAL ( TAU ) )

20         CONTINUE

C       ** FACTOR OF TWO FOR THE REAL COSINE TRANSFORM **

           CHAT(NU) = 2.0 * ( ALPHA * C(NMAX) * SIN ( OMEGA * TMAX )
     :                       + BETA * CE + GAMMA * CO ) * DT

30      CONTINUE

        RETURN
        END



        SUBROUTINE LADO ( DT, NMAX, C, CHAT )

C    *******************************************************************
C    ** CALCULATES THE FOURIER COSINE TRANSFORM BY LADO'S METHOD      **
C    **                                                               **
C    ** A CORRELATION FUNCTION, C(T), IN THE TIME DOMAIN, IS          **
C    ** TRANSFORMED TO A SPECTRUM CHAT(OMEGA) IN THE FREQUENCY DOMAIN.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** LADO, J COMPUT PHYS, 8 417, 1971.                             **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    C(NMAX)     THE CORRELATION FUNCTION.                 **
C    ** REAL    CHAT(NMAX)  THE 1-D COSINE TRANSFORM.                 **
C    ** REAL    DT          TIME INTERVAL BETWEEN POINTS IN C.        **
C    ** REAL    DOM         FREQUENCY INTERVAL BETWEEN POINTS IN CHAT.**
C    ** INTEGER NMAX        NO. OF INTERVALS ON THE TIME AXIS         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE CORRELATION FUNCTION IS REQUIRED AT HALF INTEGER          **
C    ** INTERVALS, I.E. C(T), T=(TAU-0.5)*DT FOR TAU=1 .. NMAX.       **
C    ** THE COSINE TRANSFORM IS RETURNED AT HALF INTERVALS, I.E.      **
C    ** CHAT(OMEGA), OMEGA=(NU-0.5)*DOM FOR NU = 1 .. NMAX.           **
C    *******************************************************************

        INTEGER     NMAX
        REAL        DT, C(NMAX), CHAT(NMAX)

        INTEGER     TAU, NU
        REAL        TAUH, NUH, NMAXH, PI, SUM
        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

        NMAXH = REAL ( NMAX ) - 0.5

C    ** LOOP OVER OMEGA **

        DO 20 NU = 1, NMAX

           NUH = REAL ( NU ) - 0.5
           SUM = 0.0

C       ** LOOP OVER T **

           DO 10 TAU = 1, NMAX

              TAUH =  REAL ( TAU ) - 0.5
              SUM  = SUM + C(TAU) * COS ( TAUH * NUH * PI / NMAXH )

10         CONTINUE

C       ** FACTOR OF TWO FOR THE REAL COSINE TRANSFORM **

           CHAT(NU) = 2.0 * DT * SUM

20      CONTINUE

        RETURN
        END



        SUBROUTINE FILONS ( DR, DK, NMAX, H, HHAT )

C    *******************************************************************
C    ** FOURIER SINE TRANSFORM BY FILON'S METHOD                      **
C    **                                                               **
C    ** A SPATIAL CORRELATION FUNCTION, H(R), IS TRANSFORMED TO       **
C    ** HHAT(K) IN RECIPROCAL SPACE.                                  **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FILON, PROC ROY SOC EDIN, 49 38, 1928.                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    KVEC                THE WAVENUMBER                    **
C    ** REAL    RMAX                MAXIMUM DIST IN CORREL. FUNCTION  **
C    ** REAL    ALPHA, BETA, GAMMA  FILON PARAMETERS                  **
C    ** REAL    H(NMAX)             THE CORRELATION FUNCTION          **
C    ** REAL    HHAT(NMAX)          THE 3-D TRANSFORM                 **
C    ** REAL    DR                  INTERVAL BETWEEN POINTS IN H      **
C    ** REAL    DK                  INTERVAL BETWEEN POINTS IN HHAT   **
C    ** INTEGER NMAX                NO. OF INTERVALS                  **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE REQUIRES THAT THE NUMBER OF INTERVALS, NMAX, IS   **
C    ** EVEN AND CHECKS FOR THIS CONDITION. THE FIRST VALUE OF H(R)   **
C    ** IS AT R=0. THE MAXIMUM R FOR THE CORRELATION FUNCTION IS      **
C    ** RMAX=DR*NMAX. FOR AN ACCURATE TRANSFORM H(RMAX)=0.            **
C    *******************************************************************

        INTEGER     NMAX
        REAL        DR, DK, H(0:NMAX), HHAT(0:NMAX)

        REAL        RMAX, K, THETA, SINTH, COSTH
        REAL        SINSQ, COSSQ, THSQ, THCUB, ALPHA, BETA, GAMMA
        REAL        SE, SO, FOURPI, R
        INTEGER     IR, IK

C    *******************************************************************

C    ** CHECKS NMAX IS EVEN **

        IF ( MOD ( NMAX, 2 ) .NE. 0 ) THEN

           STOP ' NMAX SHOULD BE EVEN '

        ENDIF

        FOURPI = 16.0 * ATAN ( 1.0 )
        RMAX   = REAL ( NMAX ) * DR

C    ** LOOP OVER K **

        DO 30 IK = 0, NMAX

           K  = REAL ( IK ) * DK
           THETA = K * DR

C       ** CALCULATE THE FILON PARAMETERS **

           SINTH = SIN ( THETA )
           COSTH = COS ( THETA )
           SINSQ = SINTH * SINTH
           COSSQ = COSTH * COSTH
           THSQ  = THETA * THETA
           THCUB = THSQ * THETA

           IF ( THETA. EQ. 0.0 ) THEN

              ALPHA = 0.0
              BETA  = 2.0 / 3.0
              GAMMA = 4.0 / 3.0

            ELSE

              ALPHA = ( 1.0 / THCUB )
     :               * ( THSQ + THETA * SINTH * COSTH - 2.0 * SINSQ )
              BETA  = ( 2.0 / THCUB )
     :               * ( THETA * ( 1.0 + COSSQ ) -2.0 * SINTH * COSTH )
              GAMMA = ( 4.0 / THCUB ) * ( SINTH - THETA * COSTH )

           ENDIF

C       ** THE INTEGRAND IS H(R) * R FOR THE 3-D TRANSFORM **

C       ** DO THE SUM OVER THE EVEN ORDINATES **

           SE = 0.0

           DO 10 IR = 0, NMAX, 2

              R = REAL ( IR ) * DR
              SE = SE + H(IR) * R * SIN ( K * R )

10         CONTINUE

C       ** SUBTRACT HALF THE FIRST AND LAST TERMS **
C       ** HERE THE FIRST TERM IS ZERO            **

           SE = SE - 0.5 * ( H(NMAX) * RMAX * SIN ( K * RMAX ) )

C       ** DO THE SUM OVER THE ODD ORDINATES **

           SO = 0.0

           DO 20 IR = 1, NMAX - 1, 2

              R = REAL ( IR ) * DR
              SO = SO + H(IR) * R * SIN ( K * R )

20         CONTINUE

           HHAT(IK) = ( - ALPHA * H(NMAX) * RMAX * COS ( K * RMAX)
     :                 + BETA * SE + GAMMA * SO ) * DR

C       ** INCLUDE NORMALISING FACTOR **

           HHAT(IK) = FOURPI * HHAT(IK) / K

30      CONTINUE

        RETURN
        END

