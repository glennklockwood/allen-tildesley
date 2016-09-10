********************************************************************************
** FICHE F.24.  INITIAL VELOCITY DISTRIBUTION                                 **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** CENTRE OF MASS AND ANGULAR VELOCITIES FOR LINEAR MOLECULES    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   THE NUMBER OF MOLECULES           **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    EX(N),EY(N),EZ(N)   ORIENTATIONS                      **
C    ** REAL    OX(N),OY(N),OZ(N)   SPACE-FIXED ANGULAR VELOCITIES    **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    INERT               REDUCED MOMENT OF INERTIA         **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE COMVEL ( TEMP )                                    **
C    **    SETS THE CENTRE OF MASS VELOCITIES FOR A CONFIGURATION OF  **
C    **    LINEAR MOLECULES AT A GIVEN TEMPERATURE.                   **
C    ** SUBROUTINE ANGVEL ( TEMP, INERT )                             **
C    **    SETS THE ANGULAR VELOCITIES FOR A CONFIGURATION OF LINEAR  **
C    **    MOLECULES AT A GIVEN TEMPERATURE.                          **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** WE ASSUME UNIT MOLECULAR MASS AND EMPLOY LENNARD-JONES UNITS  **
C    **       PROPERTY                      UNITS                     **
C    **       RX, RY, RZ           (EPSILON/M)**(1.0/2.0)             **
C    **       OX, OY, OZ           (EPSILON/M*SIGMA**2)**(1.0/2.0)    **
C    **       INERT                 M*SIGMA**2                        **
C    *******************************************************************

        SUBROUTINE COMVEL ( TEMP )

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ,
     :                    VX, VY, VZ, OX, OY, OZ

C    *******************************************************************
C    ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
C    **                                                               **
C    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
C    ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
C    ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)
        REAL        VX(N), VY(N), VZ(N), OX(N), OY(N), OZ(N)
        REAL        TEMP

        REAL        RTEMP, SUMX, SUMY, SUMZ
        REAL        GAUSS, DUMMY
        INTEGER     I

C    *******************************************************************

        RTEMP = SQRT ( TEMP )

        DO 100 I = 1, N

           VX(I) = RTEMP * GAUSS ( DUMMY )
           VY(I) = RTEMP * GAUSS ( DUMMY )
           VZ(I) = RTEMP * GAUSS ( DUMMY )

100     CONTINUE

C    ** REMOVE NET MOMENTUM **

        SUMX = 0.0
        SUMY = 0.0
        SUMZ = 0.0

        DO 200 I = 1, N

           SUMX = SUMX + VX(I)
           SUMY = SUMY + VY(I)
           SUMZ = SUMZ + VZ(I)

200     CONTINUE

        SUMX = SUMX / REAL ( N )
        SUMY = SUMY / REAL ( N )
        SUMZ = SUMZ / REAL ( N )

        DO 300 I = 1, N

           VX(I) = VX(I) - SUMX
           VY(I) = VY(I) - SUMY
           VZ(I) = VZ(I) - SUMZ

300     CONTINUE

        RETURN
        END



        SUBROUTINE ANGVEL ( TEMP, INERT )

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ,
     :                    VX, VY, VZ, OX, OY, OZ

C    *******************************************************************
C    ** ANGULAR VELOCITIES FROM THE MAXWELL-BOLTZMANN DISTRIBUTION.   **
C    **                                                               **
C    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND INERTIA.    **
C    ** THIS ROUTINE IS SPECIFIC TO LINEAR MOLECULES.                 **
C    ** IT CHOOSES THE DIRECTION OF THE ANGULAR VELOCITY RANDOMLY BUT **
C    ** PERPENDICULAR TO THE MOLECULAR AXIS. THE SQUARE OF THE        **
C    ** MAGNITUDE OF THE ANGULAR VELOCITY IS CHOSEN FROM AN           **
C    ** EXPONENTIAL DISTRIBUTION. THERE IS NO ATTEMPT TO SET THE      **
C    ** TOTAL ANGULAR MOMENTUM TO ZERO.                               **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)
        REAL        VX(N), VY(N), VZ(N), OX(N), OY(N), OZ(N)
        REAL        TEMP, INERT

        REAL        NORM, DOT, OSQ, O, MEAN
        REAL        XISQ, XI1, XI2, XI
        REAL        RANF, DUMMY
        INTEGER     I

C       ****************************************************************

        MEAN = 2.0 * TEMP / INERT

C    ** SET DIRECTION OF THE ANGULAR VELOCITY **

        DO 100 I = 1, N

C       ** CHOOSE A RANDOM VECTOR IN SPACE **

           XISQ = 1.0

1000       IF ( XISQ .GE. 1.0 ) THEN

              XI1  = RANF ( DUMMY ) * 2.0 - 1.0
              XI2  = RANF ( DUMMY ) * 2.0 - 1.0
              XISQ = XI1 * XI1 + XI2 * XI2

              GO TO 1000

           ENDIF

           XI    = SQRT ( 1.0 - XISQ )
           OX(I) = 2.0 * XI1 * XI
           OY(I) = 2.0 * XI2 * XI
           OZ(I) = 1.0 - 2.0 * XISQ

C       ** CONSTRAIN THE VECTOR TO BE PERPENDICULAR TO THE MOLECULE **

           DOT   = OX(I) * EX(I) + OY(I) * EY(I) + OZ(I) * EZ(I)
           OX(I) = OX(I) - DOT * EX(I)
           OY(I) = OY(I) - DOT * EY(I)
           OZ(I) = OZ(I) - DOT * EZ(I)

C       ** RENORMALIZE **

           OSQ   = OX(I) * OX(I) + OY(I) * OY(I) + OZ(I) * OZ(I)
           NORM  = SQRT ( OSQ )
           OX(I) = OX(I) / NORM
           OY(I) = OY(I) / NORM
           OZ(I) = OZ(I) / NORM

C       ** CHOOSE THE MAGNITUDE OF THE ANGULAR VELOCITY **

           OSQ   = - MEAN * LOG ( RANF ( DUMMY ) )
           O     = SQRT ( OSQ )
           OX(I) = O * OX(I)
           OY(I) = O * OY(I)
           OZ(I) = O * OZ(I)

100     CONTINUE

        RETURN
        END



        REAL FUNCTION GAUSS ( DUMMY )

C    *******************************************************************
C    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
C    **                                                               **
C    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
C    **    ADDISON-WESLEY), 1978                                      **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    *******************************************************************

        REAL        A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138, A3 = 0.252408784 )
        PARAMETER ( A5 = 0.076542912, A7 = 0.008355968 )
        PARAMETER ( A9 = 0.029899776                   )

        REAL        SUM, R, R2
        REAL        RANF, DUMMY
        INTEGER     I

C    *******************************************************************

        SUM = 0.0

        DO 10 I = 1, 12

           SUM = SUM + RANF ( DUMMY )

10      CONTINUE

        R  = ( SUM - 6.0 ) / 4.0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 )
     :          * R

        RETURN
        END



        REAL FUNCTION RANF ( DUMMY )

C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    **                                                               **
C    **                 ***************                               **
C    **                 **  WARNING  **                               **
C    **                 ***************                               **
C    **                                                               **
C    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
C    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
C    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        REAL        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /

C    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = REAL ( SEED ) / M

        RETURN
        END



