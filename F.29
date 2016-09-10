********************************************************************************
** FICHE F.29.  CONSTANT-NVT MOLECULAR DYNAMICS - CONSTRAINT METHOD           **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** CONSTANT-TEMPERATURE MOLECULAR DYNAMICS USING CONSTRAINT.     **
C    **                                                               **
C    ** THE METHOD EMPLOYED HERE IS A MODIFICATION OF THE LEAP-FROG   **
C    ** ALGORITHM (SEE ALSO FICHE F.3).                               **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** HOOVER, LADD, AND MORAN, PHYS REV LETT 48, 1818, 1982.        **
C    ** EVANS, J CHEM PHYS 78, 3297, 1983.                            **
C    ** BROWN AND CLARKE, MOL PHYS, 51, 1243, 1984.                   **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVEA ( DT, M, TEMPER )                            **
C    **    FIRST PART OF MOVE WITH VELOCITY CONSTRAINTS APPLIED.      **
C    ** SUBROUTINE MOVEB ( DT )                                       **
C    **    SECOND PART OF MOVE.                                       **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** REAL    DT                  TIMESTEP                          **
C    ** REAL    M                   ATOMIC MASS                       **
C    ** REAL    TEMPER              TEMPERATURE                       **
C    ** REAL    CHI                 SCALING PARAMETER                 **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    FX(N),FY(N),FZ(N)   FORCES                            **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE FORCE ROUTINE SHOULD BE CALLED FIRST.  THEN, WITHIN THE   **
C    ** MAIN LOOP, MOVEA ADVANCES VELOCITIES WITH THE CONSTRAINT OF   **
C    ** CONSTANT KINETIC ENERGY APPLIED.  AFTER THE ACCUMULATION      **
C    ** OF THERMODYNAMIC DATA, MOVEB MAKES THE POSITIONAL MOVE TO     **
C    ** TIME T+DT AND A NEW CALL TO THE FORCE ROUTINE MAY BE MADE.    **
C    ** THIS COMPLETES A STEP.                                        **
C    *******************************************************************



        SUBROUTINE MOVEA ( DT, M, TEMPER )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** FIRST PART OF THE CONSTANT TEMPERATURE ALGORITHM.             **
C    **                                                               **
C    ** THE FIRST PART OF THE ALGORITHM ADVANCES VELOCITIES FROM      **
C    ** T-DT/2 TO T, WITHOUT CONSTRAINT, AND THEN CALCULATES THE      **
C    ** SCALING FACTOR CHI.  THIS IS THEN USED TO PROPERLY ADVANCE    **
C    ** THE VELOCITIES FROM T-DT/2 TO T+DT/2 WITH CONSTRAINT APPLIED. **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT, M, TEMPER
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL        DT2, CHI, FREE, TEMP, K
        REAL        VXI, VYI, VZI

C    *******************************************************************

        FREE  = REAL ( ( N - 1 ) * 3 )
        DT2   = DT / 2.0
        K     = 0.0

C    ** CALCULATE THE UNCONSTRAINED VELOCITIES AT TIME T **

        DO 100 I = 1, N

           VXI = VX(I) + DT2 * FX(I) / M
           VYI = VY(I) + DT2 * FY(I) / M
           VZI = VZ(I) + DT2 * FZ(I) / M
           K   =  K + VXI * VXI + VYI * VYI + VZI * VZI

100     CONTINUE

C    ** CALCULATE THE SCALING FACTOR CHI **

        TEMP = M * K / FREE
        CHI = SQRT ( TEMPER / TEMP )

C    ** CALCULATE THE CONSTRAINED VELOCITIES AT TIME T+DT/2 **

        DO 200 I = 1, N

           VX(I) = VX(I) * ( 2.0 * CHI - 1.0 ) + CHI * DT * FX(I) / M
           VY(I) = VY(I) * ( 2.0 * CHI - 1.0 ) + CHI * DT * FY(I) / M
           VZ(I) = VZ(I) * ( 2.0 * CHI - 1.0 ) + CHI * DT * FZ(I) / M

200     CONTINUE

        RETURN
        END



        SUBROUTINE MOVEB ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** SECOND PART OF THE CONSTANT TEMPERATURE ALGORITHM             **
C    **                                                               **
C    ** THIS ADVANCES THE POSITIONS FROM T TO T + DT.                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        DT
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I

C    *******************************************************************

        DO 300 I = 1, N

           RX(I) = RX(I) + DT * VX(I)
           RY(I) = RY(I) + DT * VY(I)
           RZ(I) = RZ(I) + DT * VZ(I)

300     CONTINUE

        RETURN
        END



