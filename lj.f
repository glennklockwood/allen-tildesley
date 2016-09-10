********************************************************************************
** FICHE F.17.  A SIMPLE LENNARD-JONES FORCE ROUTINE                          **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        SUBROUTINE FORCE ( EPSLON, SIGMA, RCUT, BOX, V, VC, W )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** FORCE CALCULATION FOR LENNARD-JONES ATOMS.                    **
C    **                                                               **
C    ** IN THIS WE AIM TO SHOW HOW THE FORCES, POTENTIAL ENERGY AND   **
C    ** VIRIAL FUNCTION ARE CALCULATED IN A FAIRLY EFFICIENT WAY.     **
C    ** UNDOUBTEDLY FURTHER IMPROVEMENT WOULD BE POSSIBLE ON SPECIFIC **
C    ** MACHINES.                                                     **
C    ** THE POTENTIAL IS V(R) = 4*EPSLON*((SIGMA/R)**12-(SIGMA/R)**6) **
C    ** WE INCLUDE SPHERICAL CUTOFF AND MINIMUM IMAGING IN CUBIC BOX. **
C    ** THE BOX LENGTH IS BOX.  THE CUTOFF IS RCUT.                   **
C    ** THE ROUTINE ACTUALLY RETURNS TWO DIFFERENT POTENTIAL ENERGIES.**
C    ** V IS CALCULATED USING THE LENNARD-JONES POTENTIAL TO BE USED  **
C    ** FOR CALCULATING THE THERMODYNAMIC INTERNAL ENERGY.            **
C    ** LONG-RANGE CORRECTIONS SHOULD BE APPLIED TO THIS OUTSIDE THE  **
C    ** ROUTINE, IN THE FORM                                          **
C    **         SR3 = ( SIGMA / RCUT ) ** 3                           **
C    **         SR9 = SR3 ** 3                                        **
C    **         DENS = REAL(N) * ( SIGMA / BOX ) ** 3                 **
C    **         VLRC = ( 8.0 /9.0 ) * PI * EPSLON * DENS * REAL ( N ) **
C    **      :           * ( SR9 - 3.0 * SR3 )                        **
C    **         WLRC = ( 16.0 / 9.0 ) * PI * EPSLON * DENS * REAL( N )**
C    **      :           * ( 2.0 * SR9 - 3.0 * SR3 )                  **
C    **         V = V + VLRC                                          **
C    **         W = W + WLRC                                          **
C    ** VC IS CALCULATED USING THE SHIFTED LENNARD-JONES POTENTIAL,   **
C    ** WITH NO DISCONTINUITY AT THE CUTOFF, TO BE USED IN ASSESSING  **
C    ** ENERGY CONSERVATION.                                          **
C    ** NO REDUCED UNITS ARE USED: FOR THIS POTENTIAL WE COULD SET    **
C    ** EPSLON = 1 AND EITHER SIGMA = 1 OR BOX = 1 TO IMPROVE SPEED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    ** REAL    VX(N),VY(N),VZ(N) MOLECULAR VELOCITIES (NOT USED)     **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES                    **
C    ** REAL    SIGMA             PAIR POTENTIAL LENGTH PARAMETER     **
C    ** REAL    EPSLON            PAIR POTENTIAL ENERGY PARAMETER     **
C    ** REAL    RCUT              PAIR POTENTIAL CUTOFF               **
C    ** REAL    BOX               SIMULATION BOX LENGTH               **
C    ** REAL    V                 POTENTIAL ENERGY                    **
C    ** REAL    VC                SHIFTED POTENTIAL                   **
C    ** REAL    W                 VIRIAL FUNCTION                     **
C    ** REAL    VIJ               PAIR POTENTIAL BETWEEN I AND J      **
C    ** REAL    WIJ               NEGATIVE OF PAIR VIRIAL FUNCTION W  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        SIGMA, EPSLON, RCUT, BOX, V, VC, W
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)

        INTEGER     I, J, NCUT
        REAL        BOXINV, RCUTSQ, SIGSQ, EPS4, EPS24
        REAL        RXI, RYI, RZI, FXI, FYI, FZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL        SR2, SR6, SR12, VIJ, WIJ, FIJ

C    *******************************************************************

C    ** CALCULATE USEFUL QUANTITIES **

        BOXINV = 1.0 / BOX
        RCUTSQ = RCUT ** 2
        SIGSQ  = SIGMA ** 2
        EPS4   = EPSLON * 4.0
        EPS24  = EPSLON * 24.0

C    ** ZERO FORCES, POTENTIAL, VIRIAL **

        DO 100 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

100     CONTINUE

        NCUT = 0
        V    = 0.0
        W    = 0.0

C    ** OUTER LOOP BEGINS **

        DO 200 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           FXI = FX(I)
           FYI = FY(I)
           FZI = FZ(I)

C       ** INNER LOOP BEGINS **

           DO 199 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)
              RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
              RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
              RZIJ = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
              RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2

              IF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 SR12  = SR6 ** 2
                 VIJ   = SR12 - SR6
                 V     = V + VIJ
                 WIJ   = VIJ + SR12
                 W     = W + WIJ
                 FIJ   = WIJ / RIJSQ
                 FXIJ  = FIJ * RXIJ
                 FYIJ  = FIJ * RYIJ
                 FZIJ  = FIJ * RZIJ
                 FXI   = FXI + FXIJ
                 FYI   = FYI + FYIJ
                 FZI   = FZI + FZIJ
                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ
                 NCUT  = NCUT + 1

              ENDIF

199        CONTINUE

C       ** INNER LOOP ENDS **

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI

200     CONTINUE

C    ** OUTER LOOP ENDS **

C    ** CALCULATE SHIFTED POTENTIAL **

        SR2 = SIGSQ / RCUTSQ
        SR6 = SR2 * SR2 * SR2
        SR12 = SR6 * SR6
        VIJ = SR12 - SR6
        VC = V - REAL ( NCUT ) * VIJ

C    ** MULTIPLY RESULTS BY ENERGY FACTORS **

        DO 300 I = 1, N

           FX(I) = FX(I) * EPS24
           FY(I) = FY(I) * EPS24
           FZ(I) = FZ(I) * EPS24

300     CONTINUE

        V  = V  * EPS4
        VC = VC * EPS4
        W  = W  * EPS24 / 3.0

        RETURN
        END



