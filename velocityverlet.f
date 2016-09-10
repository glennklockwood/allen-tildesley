********************************************************************************
** FICHE F.4.  VELOCITY VERSION OF VERLET ALGORITHM                           **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO ROUTINES THAT TOGETHER IMPLEMENT VELOCITY VERLET METHOD.  **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.                   **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVEA ( DT, M )                                    **
C    **    MOVES POSITIONS AND PARTIALLY UPDATES VELOCITIES.          **
C    ** SUBROUTINE MOVEB ( DT, M, K )                                 **
C    **    COMPLETES VELOCITY MOVE AND CALCULATES KINETIC ENERGY.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** REAL    DT                  TIMESTEP                          **
C    ** REAL    M                   ATOMIC MASS                       **
C    ** REAL    K                   KINETIC ENERGY                    **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    FX(N),FY(N),FZ(N)   FORCES                            **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE    **
C    ** POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE  **
C    ** ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH        **
C    ** COMPLETES THE ADVANCEMENT OF VELOCITIES.                      **
C    *******************************************************************



        SUBROUTINE MOVEA ( DT, M )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** FIRST PART OF VELOCITY VERLET ALGORITHM                       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH      **
C    ** ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM       **
C    ** T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT, M
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL        DT2, DTSQ2

C    *******************************************************************

        DT2   = DT / 2.0
        DTSQ2 = DT * DT2

        DO 100 I = 1, N

           RX(I) = RX(I) + DT * VX(I) + DTSQ2 * FX(I) / M
           RY(I) = RY(I) + DT * VY(I) + DTSQ2 * FY(I) / M
           RZ(I) = RZ(I) + DT * VZ(I) + DTSQ2 * FZ(I) / M
           VX(I) = VX(I) + DT2 * FX(I) / M
           VY(I) = VY(I) + DT2 * FY(I) / M
           VZ(I) = VZ(I) + DT2 * FZ(I) / M

100     CONTINUE

        RETURN
        END



        SUBROUTINE MOVEB ( DT, M, K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** SECOND PART OF VELOCITY VERLET ALGORITHM                      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM     **
C    ** T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN        **
C    ** COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.       **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT, M, K
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL        DT2

C    *******************************************************************

        DT2 = DT / 2.0

        K = 0.0

        DO 200 I = 1, N

           VX(I) = VX(I) + DT2 * FX(I) / M
           VY(I) = VY(I) + DT2 * FY(I) / M
           VZ(I) = VZ(I) + DT2 * FZ(I) / M

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

200     CONTINUE

        K = 0.5 * M * K

        RETURN
        END



