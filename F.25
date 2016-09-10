********************************************************************************
** FICHE F.25.  ROUTINE TO CALCULATE TRANSLATIONAL ORDER PARAMETER            **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        SUBROUTINE ORDER ( KLATX, KLATY, KLATZ, PARAM )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** CALCULATION OF TRANSLATIONAL ORDER PARAMETER (MELTING FACTOR).**
C    **                                                               **
C    ** CLASSICALLY, THE ORDER PARAMETER IS A NORMALIZED SUM OF       **
C    ** COSINE TERMS WHICH SHOULD BE UNITY IN THE PERFECT LATTICE     **
C    ** AND FLUCTUATE AROUND ZERO FOR A DISORDERED SYSTEM.            **
C    ** HOWEVER, THIS IS NOT ORIGIN-INDEPENDENT: WITH AN UNSUITABLE   **
C    ** CHOICE OF ORIGIN IT COULD VANISH EVEN IN A PERFECT LATTICE.   **
C    ** ACCORDINGLY, WE CALCULATE HERE A QUANTITY THAT IS INDEPENDENT **
C    ** OF THE ORIGIN OF COORDINATES.                                 **
C    ** IT SHOULD BE UNITY IN A LATTICE FOR WHICH A RECIPROCAL VECTOR **
C    ** (KLATX,KLATY,KLATZ) IS SUPPLIED.                              **
C    ** IT SHOULD BE POSITIVE BUT SMALL, OF ORDER SQRT(N) IN A        **
C    ** DISORDERED SYSTEM.                                            **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR COORDINATES               **
C    ** REAL    VX(N),VY(N),VZ(N) MOLECULAR VELOCITIES (NOT USED)     **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES (NOT USED)         **
C    ** REAL    KLATX,KLATY,KLATZ RECIPROC. VECTOR OF INITIAL LATTICE **
C    ** REAL    PARAM             RESULT: ORDER PARAMETER             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        KLATX, KLATY, KLATZ, PARAM
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL        SINSUM, COSSUM

C    *******************************************************************

        SINSUM = 0.0
        COSSUM = 0.0

        DO 100 I = 1, N

           COSSUM = COSSUM + COS (  KLATX * RX(I)
     :                            + KLATY * RY(I)
     :                            + KLATZ * RZ(I) )
           SINSUM = SINSUM + SIN (  KLATX * RX(I)
     :                            + KLATY * RY(I)
     :                            + KLATZ * RZ(I) )

100     CONTINUE

        COSSUM = COSSUM / REAL ( N )
        SINSUM = SINSUM / REAL ( N )
        PARAM  = SQRT ( COSSUM ** 2 + SINSUM ** 2 )

        RETURN
        END



