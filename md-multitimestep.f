********************************************************************************
** FICHE F.21.  MULTIPLE TIMESTEP MOLECULAR DYNAMICS                          **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

        SUBROUTINE FORMTS ( STEP, DT, RCUT, RPRIM, BOX, NTS, V, W )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / POINT, LIST
        COMMON / BLOCK3 / FXS, FYS, FZS, FX1, FY1, FZ1,
     :                    FX2, FY2, FZ2, FX3, FY3, FZ3,
     :                    VS, V1, V2, V3, WS, W1, W2, W3, START

C    *******************************************************************
C    ** CALCULATES THE FORCE ON AN ATOM USING THE MTS METHOD          **
C    **                                                               **
C    ** PRINCIPAL VARIABLES                                           **
C    **                                                               **
C    ** REAL     STEP                 NUMBER OF CURRENT TIME STEP     **
C    ** REAL     RCUT                 CUTOFF DISTANCE FOR THE FORCE   **
C    ** REAL     RPRIM                RADIUS OF THE PRIMARY LIST      **
C    ** REAL     BOX                  BOX LENGTH IN SIGMA             **
C    ** REAL     V                    POTENTIAL ENERGY                **
C    ** REAL     W                    VIRIAL                          **
C    ** REAL     RX(N),RY(N),RZ(N)    ATOM POSITIONS                  **
C    ** REAL     RX1(N),RY1(N),RZ1(N) VELOCITIES                      **
C    ** REAL     FX(N),FY(N),FZ(N)    FORCE ON AN ATOM                **
C    ** REAL     FXP(N),FYP(N),FZP(N) FORCE FROM PRIMARY ATOMS        **
C    ** REAL     FXS(N),FYS(N),FZS(N) FORCE FROM SECONDARY ATOMS      **
C    ** REAL     FX1(N),FY1(N),FZ1(N) 1ST DERIVATIVE SECONDARY FORCE  **
C    ** REAL     FX2(N),FY2(N),FZ2(N) 2ND DERIVATIVE SECONDARY FORCE  **
C    ** REAL     FY3(N),FZ3(N),FZ3(N) 3RD DERIVATIVE SECONDARY FORCE  **
C    ** REAL     VP, WP               PRIMARY ENERGY AND VIRIAL       **
C    ** REAL     VS, WS               SECONDARY ENERGY AND VIRIAL     **
C    ** REAL     V1, V2, V3           DERIVATIVES OF SECONDARY ENERGY **
C    ** REAL     W1, W2, W3           DERIVATIVES OF SECONDARY VIRIAL **
C    ** INTEGER  START                STEP ZERO FOR EXTRAPOLATION     **
C    ** INTEGER  POINT(N)             INDEX TO THE NEIGHBOUR LIST     **
C    ** INTEGER  LIST(MAXNAB)         LIST OF PRIMARY NEIGHBOURS      **
C    ** INTEGER  NTS                  STEPS BETWEEN RECALCULATION OF  **
C    **                               THE SECONDARY FORCE             **
C    ** LOGICAL  MTS                  TRUE FOR AN EXTRAPOLATED STEP   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FORMTS IS CALLED IN TWO MODES. IF MTS IS FALSE, THE PRIMARY   **
C    ** FORCES ARE CALCULATED EXPLICITLY. THE SECONDARY FORCES AND    **
C    ** THEIR THREE DERIVATIVES ARE CALCULATED AND STORED. A LIST OF  **
C    ** PRIMARY NEIGHBOURS IS COMPILED. WHILE MTS IS TRUE FOR THE     **
C    ** NEXT NTS-1 STEPS, THE PRIMARY FORCES ARE CALCULATED           **
C    ** EXPLICITLY AND THE SECONDARY FORCES ESTIMATED FROM A TAYLOR   **
C    ** SERIES APPROXIMATION. THE SECONDARY FORCES, THEIR DERIVATIVES **
C    ** AND AN EXTRAPOLATION STEP MARKER ARE STORED IN COMMON BLOCK3  **
C    ** FOR USE DURING THE EXTRAPOLATION STEPS.                       **
C    ** IN THIS EXAMPLE WE TAKE THE SHIFTED-FORCE LENNARD-JONES       **
C    ** PAIR POTENTIAL.  IN PRACTICE THE MTS METHOD IS MOST EFFECTIVE **
C    ** FOR MORE COMPLICATED MOLECULAR SYSTEMS.                       **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THIS ROUTINE USES LENNARD-JONES UNITS. THE BOX IS A CUBE      **
C    ** CENTRED AT THE ORIGIN.                                        **
C    *******************************************************************

        INTEGER     N, MAXNAB

        PARAMETER ( N = 108, MAXNAB = 25 * N )

        REAL        RX(N), RY(N), RZ(N), RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N), RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N), RCUT, BOX, V, W, RPRIM
        REAL        FXS(N), FYS(N), FZS(N), FX1(N), FY1(N), FZ1(N)
        REAL        FX2(N), FY2(N), FZ2(N), FX3(N), FY3(N), FZ3(N)
        REAL        VS, V1, V2, V3
        REAL        WS, W1, W2, W3
        INTEGER     POINT(N), LIST(MAXNAB), NTS, STEP, START

        REAL        DT
        REAL        RXI, RYI, RZI, RX1I, RY1I, RZ1I, RIJ
        REAL        RX2I, RY2I, RZ2I, RX3I, RY3I, RZ3I
        REAL        FXP(N), FYP(N), FZP(N), VP, WP
        REAL        FXPI, FYPI, FZPI, FXSI, FYSI, FZSI
        REAL        FX1I, FY1I, FZ1I, FX2I, FY2I, FZ2I
        REAL        FX3I, FY3I, FZ3I
        REAL        RIJSQ, WIJ, VIJ, FIJ
        REAL        SR, SR2, SR3, SR5, SR6, SR7, SR8, SR10, SR12, SR14
        REAL        BOXINV, RCUTSQ, RPRMSQ, SF1, SF2, SRCUT
        REAL        RXIJ, RYIJ, RZIJ, RX1IJ, RY1IJ, RZ1IJ
        REAL        RX2IJ, RY2IJ, RZ2IJ, RX3IJ, RY3IJ, RZ3IJ
        REAL        R0R1, R1R1, R0R2, R0R3, R1R2
        REAL        FXIJ, FYIJ, FZIJ, FX1IJ, FY1IJ, FZ1IJ
        REAL        FX2IJ, FY2IJ, FZ2IJ, FX3IJ, FY3IJ, FZ3IJ
        REAL        BB, BC, BD, X, Y, Z, A, B, C, D
        REAL        ESTEP, ET1, ET2, ET3
        INTEGER     I, J, NLIST, JBEG, JEND, JNAB
        LOGICAL     MTS

C    *******************************************************************

        RPRMSQ = RPRIM * RPRIM
        RCUTSQ = RCUT  * RCUT
        BOXINV = 1.0 / BOX

C    ** SHIFTED FORCE CONSTANTS **

        SRCUT = 1.0 / RCUT
        SF1   = 48.0 * SRCUT ** 13  - 24.0 * SRCUT **  7
        SF2   = 28.0 * SRCUT **  6  - 52.0 * SRCUT ** 12

C ** IDENTIFY FIRST STEP IN SEQUENCE **

        IF ( MOD ( ( STEP - 1 ), NTS ) .EQ. 0 ) THEN

           MTS   = .FALSE.
           START = STEP

        ELSE

           MTS   = .TRUE.

        ENDIF

        IF ( .NOT. MTS ) THEN

C       ** AN INITIAL STEP IN A GROUP OF NTS STEPS **
C       ** COMPLETE CALCULATION OF FORCES REQUIRED **

           NLIST = 0

           DO 100 I = 1, N

C          ** ZERO PRIMARY FORCES **

              FXP(I) = 0.0
              FYP(I) = 0.0
              FZP(I) = 0.0

C          ** ZERO SECONDARY FORCES **

              FXS(I) = 0.0
              FYS(I) = 0.0
              FZS(I) = 0.0

C          ** ZERO DERIVATIVES OF SECONDARY FORCES **

              FX1(I) = 0.0
              FY1(I) = 0.0
              FZ1(I) = 0.0
              FX2(I) = 0.0
              FY2(I) = 0.0
              FZ2(I) = 0.0
              FX3(I) = 0.0
              FY3(I) = 0.0
              FZ3(I) = 0.0

100        CONTINUE

C       ** ZERO PRIMARY AND SECONDARY ENERGY AND VIRIAL **

           VP = 0.0
           WP = 0.0
           VS = 0.0
           WS = 0.0

C      ** ZERO SECONDARY ENERGY AND VIRIAL DERIVATIVES **

           V1 = 0.0
           V2 = 0.0
           V3 = 0.0
           W1 = 0.0
           W2 = 0.0
           W3 = 0.0

C       ** BEGINNING OF DOUBLE LOOP OVER ATOMS **

           DO 102 I = 1, N - 1

C          ** PLACE ARRAY ELEMENTS IN LOCAL VARIABLES **

              RXI  = RX(I)
              RYI  = RY(I)
              RZI  = RZ(I)
              RX1I = RX1(I)
              RY1I = RY1(I)
              RZ1I = RZ1(I)
              RX2I = RX2(I)
              RY2I = RY2(I)
              RZ2I = RZ2(I)
              RX3I = RX3(I)
              RY3I = RY3(I)
              RZ3I = RZ3(I)

              FXPI = FXP(I)
              FYPI = FYP(I)
              FZPI = FZP(I)
              FXSI = FXS(I)
              FYSI = FYS(I)
              FZSI = FZS(I)

              FX1I = FX1(I)
              FY1I = FY1(I)
              FZ1I = FZ1(I)
              FX2I = FX2(I)
              FY2I = FY2(I)
              FZ2I = FZ2(I)
              FX3I = FX3(I)
              FY3I = FY3(I)
              FZ3I = FZ3(I)

              POINT (I) = NLIST + 1

              DO 101 J = I + 1, N

                 RXIJ = RXI - RX(J)
                 RYIJ = RYI - RY(J)
                 RZIJ = RZI - RZ(J)

                 RXIJ  = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
                 RYIJ  = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                 RZIJ  = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

C                ** PAIR INSIDE CUTOFF **

                    RIJ = SQRT ( RIJSQ )
                    SR  = 1.0 / RIJ
                    SR2 = SR * SR
                    SR6 = SR2 * SR2 * SR2
                    VIJ =  4.0 * SR6 * ( SR6 - 1.0 ) + SF1 * RIJ + SF2
                    WIJ = 48.0 * SR6 * ( SR6 - 0.5 ) - SF1 * RIJ
                    FIJ = SR2 * WIJ

                    FXIJ = RXIJ * FIJ
                    FYIJ = RYIJ * FIJ
                    FZIJ = RZIJ * FIJ

                    IF ( RIJSQ .LT. RPRMSQ ) THEN

C                   ** J IS A PRIMARY NEIGHBOUR OF I **

                       NLIST       = NLIST + 1
                       LIST(NLIST) = J

                       FXPI   = FXPI   + FXIJ
                       FYPI   = FYPI   + FYIJ
                       FZPI   = FZPI   + FZIJ
                       FXP(J) = FXP(J) - FXIJ
                       FYP(J) = FYP(J) - FYIJ
                       FZP(J) = FZP(J) - FZIJ
                       VP     = VP     + VIJ
                       WP     = WP     + WIJ

                    ELSE

C                   ** J IS A SECONDARY NEIGHBOUR OF I **

                       FXSI   = FXSI   + FXIJ
                       FYSI   = FYSI   + FYIJ
                       FZSI   = FZSI   + FZIJ
                       FXS(J) = FXS(J) - FXIJ
                       FYS(J) = FYS(J) - FYIJ
                       FZS(J) = FZS(J) - FZIJ
                       VS     = VS     + VIJ
                       WS     = WS     + WIJ

C                   ** CALCULATE THE FIRST THREE DERIVATIVES **

                       SR3   = SR2 * SR
                       SR5   = SR3 * SR2
                       SR7   = SR5 * SR2
                       SR8   = SR6 * SR2
                       SR10  = SR5 * SR5
                       SR12  = SR10 * SR2
                       SR14  = SR12 * SR2

                       A = FIJ
                       B =   192.0 * SR10 * ( 1.0 - 3.5 * SR6 )
     :                                              +        SR3 * SF1
                       C =  1920.0 * SR12 * ( 5.6 * SR6 - 1.0 )
     :                                              -  3.0 * SR5 * SF1
                       D = 23040.0 * SR14 * ( 1.0 - 8.4 * SR6 )
     :                                              + 15.0 * SR7 * SF1

C                   ** DERIVATIVES OF PAIR SEPARATIONS **

                       RX1IJ = RX1I - RX1(J)
                       RY1IJ = RY1I - RY1(J)
                       RZ1IJ = RZ1I - RZ1(J)
                       RX2IJ = RX2I - RX2(J)
                       RY2IJ = RY2I - RY2(J)
                       RZ2IJ = RZ2I - RZ2(J)
                       RX3IJ = RX3I - RX3(J)
                       RY3IJ = RY3I - RY3(J)
                       RZ3IJ = RZ3I - RZ3(J)

C                   ** APPROPRIATE DOT PRODUCTS **

                       R0R1 =  RXIJ*RX1IJ +  RYIJ*RY1IJ +  RZIJ*RZ1IJ
                       R0R2 =  RXIJ*RX2IJ +  RYIJ*RY2IJ +  RZIJ*RZ2IJ
                       R0R3 =  RXIJ*RX3IJ +  RYIJ*RY3IJ +  RZIJ*RZ3IJ
                       R1R1 = RX1IJ*RX1IJ + RY1IJ*RY1IJ + RZ1IJ*RZ1IJ
                       R1R2 = RX1IJ*RX2IJ + RY1IJ*RY2IJ + RZ1IJ*RZ2IJ

C                   ** FIRST DERIVATIVES OF SECONDARY FORCES **

                       BB    = B * R0R1
                       FX1IJ = A * RX1IJ + BB * RXIJ
                       FY1IJ = A * RY1IJ + BB * RYIJ
                       FZ1IJ = A * RZ1IJ + BB * RZIJ

                       FX1I   = FX1I   + FX1IJ
                       FY1I   = FY1I   + FY1IJ
                       FZ1I   = FZ1I   + FZ1IJ
                       FX1(J) = FX1(J) - FX1IJ
                       FY1(J) = FY1(J) - FY1IJ
                       FZ1(J) = FZ1(J) - FZ1IJ

C                   ** SECOND DERIVATIVES OF SECONDARY FORCES **

                       BC    = B * ( R0R2 + R1R1 ) + C * R0R1 * R0R1
                       FX2IJ = BC * RXIJ + 2.0 * BB * RX1IJ + A * RX2IJ
                       FY2IJ = BC * RYIJ + 2.0 * BB * RY1IJ + A * RY2IJ
                       FZ2IJ = BC * RZIJ + 2.0 * BB * RZ1IJ + A * RZ2IJ

                       FX2I   = FX2I   + FX2IJ
                       FY2I   = FY2I   + FY2IJ
                       FZ2I   = FZ2I   + FZ2IJ
                       FX2(J) = FX2(J) - FX2IJ
                       FY2(J) = FY2(J) - FY2IJ
                       FZ2(J) = FZ2(J) - FZ2IJ

C                   ** THIRD DERIVATIVES OF SECONDARY FORCES **

                       BD = B*(R0R3 + 3.0*R1R2) + 3.0*C*(R0R1*R0R2 +
     :                      R0R1*R1R1) + D*R0R1*R0R1*R0R1
                       FX3IJ =BD*RXIJ+3.0*BC*RX1IJ+3.0*BB*RX2IJ+A*RX3IJ
                       FY3IJ =BD*RYIJ+3.0*BC*RY1IJ+3.0*BB*RY2IJ+A*RY3IJ
                       FZ3IJ =BD*RZIJ+3.0*BC*RZ1IJ+3.0*BB*RZ2IJ+A*RZ3IJ

                       FX3I   = FX3I   + FX3IJ
                       FY3I   = FY3I   + FY3IJ
                       FZ3I   = FZ3I   + FZ3IJ
                       FX3(J) = FX3(J) - FX3IJ
                       FY3(J) = FY3(J) - FY3IJ
                       FZ3(J) = FZ3(J) - FZ3IJ

C                   ** ENERGY DERIVATIVES **

                       V1 = V1 - A * R0R1
                       V2 = V2 - B*R0R1*R0R1 - A*(R0R2+R1R1)
                       V3 = V3 - C*R0R1*R0R1*R0R1 - 3.0*B*R0R1*(R0R2+
     :                           R1R1) - A*(R0R3+3.0*R1R2)

C                   ** VIRIAL DERIVATIVES **

                       X =   144.0*SR8 *( 1.0 -  4.0*SR6) - SR *SF1
                       Y =  1152.0*SR10*( 7.0*SR6 -  1.0) + SR3*SF1
                       Z = 11520.0*SR12*( 1.0 - 11.2*SR6) - SR5*SF1*3.0

                       W1 = W1 + X * R0R1
                       W2 = W2 + Y*R0R1*R0R1 + X*(R0R2+R1R1)
                       W3 = W3 + Z*R0R1*R0R1*R0R1 + 3.0*Y*R0R1*(R0R2+
     :                           R1R1) + X*(R0R3+3.0*R1R2)

                    ENDIF

                 ENDIF

101           CONTINUE

C          ** COLLECT ACCUMULATORS FOR FORCES AND DERIVATIVES **

              FXP(I) = FXPI
              FYP(I) = FYPI
              FZP(I) = FZPI
              FXS(I) = FXSI
              FYS(I) = FYSI
              FZS(I) = FZSI
              FX1(I) = FX1I
              FY1(I) = FY1I
              FZ1(I) = FZ1I
              FX2(I) = FX2I
              FY2(I) = FY2I
              FZ2(I) = FZ2I
              FX3(I) = FX3I
              FY3(I) = FY3I
              FZ3(I) = FZ3I

102        CONTINUE

C       ** END OF DOUBLE LOOP OVER ATOMS **

           POINT(N) = NLIST + 1

C       ** ADD PRIMARY AND SECONDARY FORCES **

           DO 103 I = 1, N

              FX(I) = FXP(I) + FXS(I)
              FY(I) = FYP(I) + FYS(I)
              FZ(I) = FZP(I) + FZS(I)

103        CONTINUE

           V = VP + VS
           W = WP + WS

        ELSE

C       ** A MULTIPLE TIMESTEP EXTRAPOLATION STEP **
C       ** SECONDARY FORCES FROM TAYLOR SERIES    **

C       ** ZERO PRIMARY FORCES **

           DO 104 I = 1, N

              FXP(I) = 0.0
              FYP(I) = 0.0
              FZP(I) = 0.0

104        CONTINUE

           VP = 0.0
           WP = 0.0

C       ** BEGINNING OF DOUBLE LOOP OVER PRIMARY ATOMS USING LIST **

           DO 106 I = 1, N - 1

              JBEG = POINT(I)
              JEND = POINT(I+1) - 1

              IF ( JBEG .LE. JEND ) THEN

                 RXI  = RX(I)
                 RYI  = RY(I)
                 RZI  = RZ(I)
                 FXPI = FXP(I)
                 FYPI = FYP(I)
                 FZPI = FZP(I)

                 DO 105 JNAB = JBEG, JEND

                    J = LIST(JNAB)

                    RXIJ  = RXI - RX(J)
                    RYIJ  = RYI - RY(J)
                    RZIJ  = RZI - RZ(J)

                    RXIJ  = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
                    RYIJ  = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                    RZIJ  = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
                    RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                    IF ( RIJSQ .LT. RCUTSQ ) THEN

                       RIJ = SQRT ( RIJSQ )
                       SR  = 1.0 / RIJ
                       SR2 = SR * SR
                       SR6 = SR2 * SR2 * SR2
                       VIJ = 4.0 * SR6 * ( SR6 - 1.0 ) + SF1*RIJ + SF2
                       WIJ = 48.0 * SR6 * ( SR6 - 0.5 ) - SF1 * RIJ

                       FIJ  = WIJ * SR2
                       FXIJ = FIJ * RXIJ
                       FYIJ = FIJ * RYIJ
                       FZIJ = FIJ * RZIJ

                       FXPI   = FXPI   + FXIJ
                       FYPI   = FYPI   + FYIJ
                       FZPI   = FZPI   + FZIJ
                       FXP(J) = FXP(J) - FXIJ
                       FYP(J) = FYP(J) - FYIJ
                       FZP(J) = FZP(J) - FZIJ

                       VP     = VP + VIJ
                       WP     = WP + WIJ

                    ENDIF

105              CONTINUE

                 FXP(I) = FXPI
                 FYP(I) = FYPI
                 FZP(I) = FZPI

              ENDIF

106        CONTINUE

C       ** END OF DOUBLE LOOP OVER PRIMARY ATOMS USING LIST **

C       ** ESTIMATE SECONDARY FORCES **

           ESTEP = REAL ( STEP - START )
           ET1   = ESTEP * DT
           ET2   = 0.5 * ET1 * ET1
           ET3   = ( 1.0 / 3.0 ) * ET1 * ET2

C       ** TAYLOR SERIES EXPANSION **

           DO 107 I = 1, N

              FX(I) = FXP(I)+FXS(I)+ET1*FX1(I)+ET2*FX2(I)+ET3*FX3(I)
              FY(I) = FYP(I)+FYS(I)+ET1*FY1(I)+ET2*FY2(I)+ET3*FY3(I)
              FZ(I) = FZP(I)+FZS(I)+ET1*FZ1(I)+ET2*FZ2(I)+ET3*FZ3(I)

107        CONTINUE

           V = VP + VS + ET1 * V1 + ET2 * V2 + ET3 * V3
           W = WP + WS + ET1 * W1 + ET2 * W2 + ET3 * W3

        ENDIF

        W = W / 3.0

        RETURN
        END



