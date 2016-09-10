********************************************************************************
** FICHE F.30.  CONSTANT-NPH MOLECULAR DYNAMICS - EXTENDED SYSTEM METHOD      **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM ANDERS

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** CONSTANT-NPH MOLECULAR DYNAMICS BY ANDERSEN'S METHOD.         **
C    **                                                               **
C    ** THE MODIFIED EQUATIONS OF MOTION ARE AS FOLLOWS:              **
C    **       R2 = F/M + (R/3)[V2/V - (2/3)*(V1/V)**2]                **
C    **       V2 = ( PRES - PRESUR ) / MP                             **
C    ** WHERE R,R1,R2 ARE THE ATOM POSITIONS AND THEIR DERIVATIVES,   **
C    ** V,V1,V2 ARE THE VOLUME AND ITS DERIVATIVE, AND MP IS THE      **
C    ** MASS OF THE PISTON SURROUNDING THE BOX. PRES IS THE           **
C    ** CALCULATED PRESSURE, AND PRESUR THE REQUIRED PRESSURE.        **
C    ** WE SOLVE THESE EQUATIONS BY A GEAR 4-VALUE METHOD FOR SECOND  **
C    ** ORDER DIFFERENTIAL EQUATIONS. FOLLOWING BROWN AND CLARKE, WE  **
C    ** USE UNSCALED DISTANCE VARIABLES WHICH ARE REDUCED BY SIGMA.   **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** ANDERSEN, J. CHEM. PHYS. 72, 283, 1980.                       **
C    ** HAILE AND GRABEN, J. CHEM. PHYS., 73, 2412, 1980.             **
C    ** BROWN AND CLARKE, MOLEC. PHYS., 51, 1243, 1984.               **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** REAL    DT                            TIMESTEP                **
C    ** REAL    MP                            PISTON MASS             **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS               **
C    ** REAL    RX1(N),RY1(N),RZ1(N)          FIRST DERIVATIVES       **
C    ** REAL    RX2(N),RY2(N),RZ2(N)          SECOND DERIVATIVES      **
C    ** REAL    RX3(N),RY3(N),RZ3(N)          THIRD DERIVATIVES       **
C    ** REAL    VOL,VOL1,VOL2,VOL3            VOLUME AND DERIVATIVES  **
C    ** REAL    FX(N),FY(N),FZ(N)             TOTAL FORCES            **
C    **                                                               **
C    ** ROUTINES REFERENCED                                           **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN CONFIGURATION AND BOX VARIABLES                   **
C    ** SUBROUTINE FORCE ( RCUT, V, W )                               **
C    **    CALCULATES FORCES, POTENTIAL, AND VIRIAL                   **
C    ** SUBROUTINE KINET ( K )                                        **
C    **    CALCULATES KINETIC ENERGY                                  **
C    ** SUBROUTINE PREDIC ( DT )                                      **
C    **    PREDICTOR ROUTINE FOR CONFIGURATION AND BOX VARIABLES      **
C    ** SUBROUTINE CORREC ( DT )                                      **
C    **    CORRECTOR ROUTINE FOR CONFIGURATION AND BOX VARIABLES      **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT CONFIGURATION AND BOX VARIABLES                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES

        INTEGER     STEP, NSTEP, IPRINT
        REAL        ACV, ACK, ACE, ACP, ACT, ACH, ACD
        REAL        ACVSQ, ACKSQ, ACESQ, ACPSQ, ACTSQ, ACHSQ, ACDSQ
        REAL        AVV, AVK, AVE, AVP, AVT, AVH, AVD
        REAL        FLV, FLK, FLE, FLP, FLT, FLH, FLD
        REAL        DT, DENS, TEMP, RCUT, PRES, PRESUR, NORM
        REAL        K, V, W, E, H, HAM
        REAL        KN, VN, EN, HN, HAMN
        REAL        SR3, SR9, VLRC, WLRC, VLRC0, WLRC0, PI, MP
        CHARACTER   TITLE*80, CNFILE*30
        REAL        FREE

        PARAMETER ( FREE = 3.0 )
        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

        WRITE(*,'(1H1, '' **** PROGRAM ANDERS ****                 '')')
        WRITE(*,'(//1X,'' MOLECULAR DYNAMICS OF LENNARD-JONES ATOMS'')')
        WRITE(*,'(1X,  '' CONSTANT-NPH ALGORITHM OF ANDERSEN       '')')

C    ** BASIC SIMULATION PARAMETERS **

        WRITE(*,'('' ENTER RUN TITLE                               '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER NUMBER OF STEPS                         '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER INTERVAL BETWEEN PRINTS                 '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER THE FOLLOWING IN L-J REDUCED UNITS      '')')
        WRITE(*,'('' ENTER TIMESTEP                                '')')
        READ (*,*) DT
        WRITE(*,'('' ENTER POTENTIAL CUTOFF                        '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER DESIRED PRESSURE                        '')')
        READ (*,*) PRESUR
        WRITE(*,'('' ENTER PISTON MASS PARAMETER, M                '')')
        READ (*,*) MP

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE
        WRITE(*,'('' NUMBER OF STEPS  = '',I6   )') NSTEP
        WRITE(*,'('' PRINT INTERVAL   = '',I6   )') IPRINT
        WRITE(*,'('' TIMESTEP         = '',F10.5)') DT
        WRITE(*,'('' POTENTIAL CUTOFF = '',F10.5)') RCUT
        WRITE(*,'('' DESIRED PRES.    = '',F10.5)') PRESUR
        WRITE(*,'('' M PARAMETER      = '',F10.5)') MP

C    ** READCN MUST READ IN INITIAL CONFIGURATION    **
C    ** AND ASSIGN VALUES TO BOX AND ITS DERIVATIVES **

        CALL READCN ( CNFILE )
        DENS = REAL ( N ) / VOL

        WRITE(*,'('' INITIAL DENS.    = '',F10.5)') DENS

        IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1

C    ** PREPARE FACTORS FOR LONG-RANGE CORRECTIONS **
C    ** NB: SPECIFIC TO LENNARD-JONES POTENTIAL    **

        SR3 = ( 1.0 / RCUT ) ** 3
        SR9 = SR3 ** 3
        VLRC0 = 8.0 * PI * REAL ( N ) * SR9 / 9.0
     :        - 8.0 * PI * REAL ( N ) * SR3 / 3.0
        WLRC0 = 32.0 * PI * REAL ( N ) * SR9 / 9.0
     :         - 16.0 * PI * REAL ( N ) * SR3 / 3.0

C    ** ZERO ACCUMULATORS **

        ACV = 0.0
        ACK = 0.0
        ACE = 0.0
        ACP = 0.0
        ACT = 0.0
        ACH = 0.0
        ACD = 0.0

        ACVSQ = 0.0
        ACKSQ = 0.0
        ACESQ = 0.0
        ACPSQ = 0.0
        ACTSQ = 0.0
        ACHSQ = 0.0
        ACDSQ = 0.0

        FLV = 0.0
        FLK = 0.0
        FLE = 0.0
        FLP = 0.0
        FLT = 0.0
        FLH = 0.0
        FLD = 0.0

        WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
        WRITE(*,10001)

C    *******************************************************************
C    ** MAIN LOOP STARTS                                              **
C    *******************************************************************

        DO 1000 STEP = 1, NSTEP

C       ** IMPLEMENT ALGORITHM **

           CALL PREDIC ( DT )
           CALL FORCE ( RCUT, V, W )
           CALL KINET ( K )

C       ** INCLUDE LONG-RANGE CORRECTIONS IN ALGORITHM **

           KN    = K / REAL ( N )
           DENS  = REAL ( N ) / VOL
           WLRC  = WLRC0 * DENS
           TEMP  = 2.0 * KN / FREE
           PRES  = DENS * TEMP + ( W + WLRC ) / VOL
           DPRES = ( PRES - PRESUR ) / MP

           CALL CORREC ( DT )

C       ** CALCULATE INSTANTANEOUS VALUES   **
C       ** INCLUDING LONG-RANGE CORRECTIONS **

           DENS = REAL ( N ) / VOL
           VLRC = VLRC0 * DENS
           WLRC = WLRC0 * DENS
           V    = V + VLRC
           W    = W + WLRC
           E    = K + V
           VN   = V / REAL ( N )
           EN   = E / REAL ( N )
           TEMP = 2.0 * KN / FREE
           PRES = DENS * TEMP + W / VOL
           H    = E + PRES * VOL + 0.5 * MP * VOL1 ** 2
           HAM  = E + PRESUR * VOL + 0.5 * MP * VOL1 ** 2
           HN   = H / REAL ( N )
           HAMN = HAM / REAL ( N )

C       ** INCREMENT ACCUMULATORS **

           ACE = ACE + EN
           ACK = ACK + KN
           ACV = ACV + VN
           ACP = ACP + PRES
           ACH = ACH + HN
           ACD = ACD + DENS

           ACESQ = ACESQ + EN ** 2
           ACKSQ = ACKSQ + KN ** 2
           ACVSQ = ACVSQ + VN ** 2
           ACPSQ = ACPSQ + PRES ** 2
           ACHSQ = ACHSQ + HN ** 2
           ACDSQ = ACDSQ + DENS ** 2

C       ** OPTIONALLY PRINT INFORMATION **

           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN

              WRITE(*,'(1X,I8,9(2X,F10.5))')
     :                 STEP, EN, HN, KN, VN, PRES, TEMP, DENS, HAMN, VOL

           ENDIF

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

C    ** WRITE OUT FINAL CONFIGURATION **
C    ** INCLUDING VOL,VOL1,VOL2,VOL3  **

        CALL WRITCN ( CNFILE )

C    ** WRITE OUT FINAL AVERAGES **

        NORM = REAL ( NSTEP )
        AVE  = ACE / NORM
        AVK  = ACK / NORM
        AVV  = ACV / NORM
        AVP  = ACP / NORM
        AVH  = ACH / NORM
        AVD  = ACD / NORM

        ACESQ = ( ACESQ / NORM ) - AVE ** 2
        ACKSQ = ( ACKSQ / NORM ) - AVK ** 2
        ACVSQ = ( ACVSQ / NORM ) - AVV ** 2
        ACPSQ = ( ACPSQ / NORM ) - AVP ** 2
        ACHSQ = ( ACHSQ / NORM ) - AVH ** 2
        ACDSQ = ( ACDSQ / NORM ) - AVD ** 2

        IF ( ACESQ .GT. 0.0 ) FLE = SQRT ( ACESQ )
        IF ( ACKSQ .GT. 0.0 ) FLK = SQRT ( ACKSQ )
        IF ( ACVSQ .GT. 0.0 ) FLV = SQRT ( ACVSQ )
        IF ( ACPSQ .GT. 0.0 ) FLP = SQRT ( ACPSQ )
        IF ( ACHSQ .GT. 0.0 ) FLH = SQRT ( ACHSQ )
        IF ( ACDSQ .GT. 0.0 ) FLD = SQRT ( ACDSQ )

        AVT = AVK / 1.5
        FLT = FLK / 1.5

        WRITE(*,'('' AVERAGES'',7(2X,F10.5))')
     :             AVE, AVH, AVK, AVV, AVP, AVT, AVD
        WRITE(*,'('' FLUCTS  '',7(2X,F10.5))')
     :             FLE, FLH, FLK, FLV, FLP, FLT, FLD

        STOP

10001   FORMAT(//1X,'TIMESTEP  ..ENERGY..  .ENTHALPY.  ..KINETIC.',
     :          '  ..POTENT..  .PRESSURE.  ..TEMPER..  ..DENSITY.',
     :          '  ...HAMIL..  ..VOLUME..')
        END



        SUBROUTINE FORCE ( RCUT, V, W )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** LENNARD-JONES FORCE ROUTINE IN REDUCED UNITS                  **
C    **                                                               **
C    ** THE POTENTIAL IS V(R) = 4*((1/R)**12-(1/R)**6)                **
C    ** WE INCLUDE SPHERICAL CUTOFF AND MINIMUM IMAGING IN CUBIC BOX. **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES                    **
C    ** REAL    VOL               SIMULATION VOLUME                   **
C    ** REAL    BOX               SIMULATION BOX LENGTH               **
C    ** REAL    RCUT              PAIR POTENTIAL CUTOFF               **
C    ** REAL    V                 POTENTIAL ENERGY                    **
C    ** REAL    W                 VIRIAL FUNCTION                     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES

        INTEGER     I, J
        REAL        RCUT, V, W, BOX
        REAL        BOXINV, RCUTSQ
        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, RIJSQ
        REAL        FXI, FYI, FZI, FXIJ, FYIJ, FZIJ
        REAL        SR2, SR6, SR12, VIJ, WIJ, FIJ

C    *******************************************************************

C    ** USEFUL QUANTITIES **

        BOX    = VOL ** ( 1.0 / 3.0 )
        BOXINV = 1.0 / BOX
        RCUTSQ = RCUT ** 2

C    ** ZERO FORCES, POTENTIAL, VIRIAL **

        DO 100 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

100     CONTINUE

        V   = 0.0
        W   = 0.0

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
              RIJSQ  = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2

              IF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2   = 1.0 / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )
                 V     = V + VIJ
                 WIJ   = SR6 * ( SR6 - 0.5 )
                 W     = W + WIJ
                 FIJ   = WIJ * SR2
                 FXIJ  = FIJ * RXIJ
                 FYIJ  = FIJ * RYIJ
                 FZIJ  = FIJ * RZIJ
                 FXI   = FXI + FXIJ
                 FYI   = FYI + FYIJ
                 FZI   = FZI + FZIJ
                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ

              ENDIF

199        CONTINUE

C       ** INNER LOOP ENDS **

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI

200     CONTINUE

C    ** OUTER LOOP ENDS **

C    ** MULTIPLY RESULTS BY NUMERICAL FACTORS **

        DO 300 I = 1, N

           FX(I) = FX(I) * 48.0
           FY(I) = FY(I) * 48.0
           FZ(I) = FZ(I) * 48.0

300     CONTINUE

        V   = V * 4.0
        W   = W * 48.0 / 3.0

        RETURN
        END



        SUBROUTINE KINET ( K )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** ROUTINE TO COMPUTE KINETIC ENERGY.                            **
C    **                                                               **
C    ** MOMENTUM AND VELOCITY ARE RELATED THROUGH THE FOLLOWING       **
C    ** DIFFERENTIAL EQUATION: P = R1 - V1*R/3.0/V                    **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES

        REAL        K, V1V3, PX, PY, PZ
        INTEGER     I

C    *******************************************************************

        K    = 0.0
        V1V3 = VOL1 / VOL / 3.0

        DO 1000 I = 1, N

           PX = RX1(I) - RX(I) * V1V3
           PY = RY1(I) - RY(I) * V1V3
           PZ = RZ1(I) - RZ(I) * V1V3
           K = K + PX ** 2 + PY ** 2 + PZ ** 2

1000    CONTINUE

        K = 0.5 * K

        RETURN
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT, NN
        PARAMETER ( CNUNIT = 10 )

C     ******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) NN, VOL, VOL1, VOL2, VOL3
        IF ( NN .NE. N ) STOP 'INCORRECT VALUE OF N'
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) RX1, RY1, RZ1
        READ ( CNUNIT ) RX2, RY2, RZ2
        READ ( CNUNIT ) RX3, RY3, RZ3

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION TO UNIT 10           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N, VOL, VOL1, VOL2, VOL3
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) RX1, RY1, RZ1
        WRITE ( CNUNIT ) RX2, RY2, RZ2
        WRITE ( CNUNIT ) RX3, RY3, RZ3

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE PREDIC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** STANDARD TAYLOR SERIES PREDICTORS                             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES
        REAL        DT

        INTEGER     I
        REAL        C1, C2, C3

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0

        DO 100 I = 1, N

           RX(I)  = RX(I)  + C1*RX1(I) + C2*RX2(I) + C3*RX3(I)
           RY(I)  = RY(I)  + C1*RY1(I) + C2*RY2(I) + C3*RY3(I)
           RZ(I)  = RZ(I)  + C1*RZ1(I) + C2*RZ2(I) + C3*RZ3(I)
           RX1(I) = RX1(I) + C1*RX2(I) + C2*RX3(I)
           RY1(I) = RY1(I) + C1*RY2(I) + C2*RY3(I)
           RZ1(I) = RZ1(I) + C1*RZ2(I) + C2*RZ3(I)
           RX2(I) = RX2(I) + C1*RX3(I)
           RY2(I) = RY2(I) + C1*RY3(I)
           RZ2(I) = RZ2(I) + C1*RZ3(I)

100     CONTINUE

        VOL  = VOL  + C1*VOL1 + C2*VOL2 + C3*VOL3
        VOL1 = VOL1 + C1*VOL2 + C2*VOL3
        VOL2 = VOL2 + C1*VOL3

        RETURN
        END



        SUBROUTINE CORREC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / VOL, VOL1, VOL2, VOL3, DPRES

C    *******************************************************************
C    ** GEAR CORRECTOR ALGORITHM.                                     **
C    **                                                               **
C    ** FOR TIMESTEP-SCALED VARIABLES, GEAR COEFFICIENTS WOULD BE AS  **
C    ** FOLLOWS (4-VALUE METHOD, 2ND-ORDER D.E.): 1/6, 5/6, 1, 1/3.   **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        VOL, VOL1, VOL2, VOL3, DPRES
        REAL        DT

        INTEGER     I
        REAL        C1, C2, C3, COEFF0, COEFF1, COEFF3
        REAL        CORV, CORRX, CORRY, CORRZ, VFAC
        REAL        RX2I, RY2I, RZ2I

        REAL        GEAR0, GEAR1, GEAR3
        PARAMETER ( GEAR0 = 1.0 / 6.0,
     :              GEAR1 = 5.0 / 6.0,
     :              GEAR3 = 1.0 / 3.0 )

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0

        COEFF0 = GEAR0 * C2
        COEFF1 = GEAR1 * C2 / C1
        COEFF3 = GEAR3 * C2 / C3

        VFAC = ( ( VOL2 / VOL ) - 2.0 * ( VOL1 / VOL ) ** 2 / 3.0 )
     :         / 3.0

        DO 400 I = 1, N

           RX2I = FX(I) + VFAC * RX(I)
           RY2I = FY(I) + VFAC * RY(I)
           RZ2I = FZ(I) + VFAC * RZ(I)
           CORRX = RX2I - RX2(I)
           CORRY = RY2I - RY2(I)
           CORRZ = RZ2I - RZ2(I)

           RX(I)  = RX(I)  + COEFF0 * CORRX
           RY(I)  = RY(I)  + COEFF0 * CORRY
           RZ(I)  = RZ(I)  + COEFF0 * CORRZ
           RX1(I) = RX1(I) + COEFF1 * CORRX
           RY1(I) = RY1(I) + COEFF1 * CORRY
           RZ1(I) = RZ1(I) + COEFF1 * CORRZ
           RX2(I) = RX2I
           RY2(I) = RY2I
           RZ2(I) = RZ2I
           RX3(I) = RX3(I) + COEFF3 * CORRX
           RY3(I) = RY3(I) + COEFF3 * CORRY
           RZ3(I) = RZ3(I) + COEFF3 * CORRZ

400     CONTINUE

        CORV = DPRES - VOL2
        VOL  = VOL   + COEFF0 * CORV
        VOL1 = VOL1  + COEFF1 * CORV
        VOL2 = DPRES
        VOL3 = VOL3  + COEFF3 * CORV

        RETURN
        END



