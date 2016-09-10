********************************************************************************
** FICHE F.31.  CONSTANT-NPT MOLECULAR DYNAMICS - CONSTRAINT METHOD           **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM EVAMOR

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** CONSTANT-NPT MOLECULAR DYNAMICS USING CONSTRAINT ALGORITHM.   **
C    **                                                               **
C    ** THE MODIFIED EQUATIONS OF MOTION ARE AS FOLLOWS:              **
C    **       R1 = P/M + CHI*R                                        **
C    **       P1 = F - CHI*P - XI*P                                   **
C    **       V1 = 3*V*CHI                                            **
C    ** WHERE R1 IS THE TIME DERIVATIVE OF POSITION R                 **
C    ** P1 IS THE TIME DERIVATIVE OF MOMENTUM P                       **
C    ** V1 IS THE TIME DERIVATIVE OF VOLUME V                         **
C    ** AND CHI AND XI ARE LAGRANGE MULTIPLIERS                       **
C    **       CHI = ( - SUM( (R.P) X / (R.R) ) ) / ( SUM(X) + 9PV )   **
C    **       XI  = SUM(P.F) / SUM(P.P) - CHI                         **
C    ** HERE TERMS LIKE (P.F) ARE SCALAR PRODUCTS OVER PAIR TERMS AND **
C    ** PV STANDS FOR PRESSURE * VOLUME.                              **
C    ** WE SOLVE THESE EQUATIONS BY A GEAR 4-VALUE METHOD FOR FIRST   **
C    ** ORDER DIFFERENTIAL EQUATIONS                                  **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** EVANS AND MORRISS, CHEM PHYS 77, 63, 1983.                    **
C    ** EVANS AND MORRISS, COMPUT PHYS REP 1, 297, 1984.              **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** REAL    DT                            TIMESTEP                **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS               **
C    ** REAL    RX1(N),RY1(N),RZ1(N)          FIRST DERIVATIVES       **
C    ** REAL    RX2(N),RY2(N),RZ2(N)          SECOND DERIVATIVES      **
C    ** REAL    RX3(N),RY3(N),RZ3(N)          THIRD DERIVATIVES       **
C    ** REAL    PX(N),PY(N),PZ(N)             MOMENTA                 **
C    ** REAL    PX1(N),PY1(N),PZ1(N)          FIRST DERIVATIVES       **
C    ** REAL    PX2(N),PY2(N),PZ2(N)          SECOND DERIVATIVES      **
C    ** REAL    PX3(N),PY3(N),PZ3(N)          THIRD DERIVATIVES       **
C    ** REAL    BOX,BOX1,BOX2,BOX3            BOX LENGTH AND DERIVS   **
C    ** REAL    FX(N),FY(N),FZ(N)             TOTAL FORCES            **
C    **                                                               **
C    ** ROUTINES REFERENCED                                           **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN CONFIGURATION AND BOX VARIABLES                   **
C    ** SUBROUTINE FORCE ( RCUT, V, W, X, RPX, PF )                   **
C    **    CALCULATES FORCES, POTENTIAL, VIRIAL, HYPERVIRIAL ETC.     **
C    ** SUBROUTINE KINET ( K )                                        **
C    **    CALCULATES KINETIC ENERGY                                  **
C    ** SUBROUTINE PREDIC ( DT )                                      **
C    **    PREDICTOR ROUTINE FOR CONFIGURATION AND BOX VARIABLES      **
C    ** SUBROUTINE CORREC ( DT, CHI, CHIPXI )                         **
C    **    CORRECTOR ROUTINE FOR CONFIGURATION AND BOX VARIABLES      **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT CONFIGURATION AND BOX VARIABLES                 **
C    ** SUBROUTINE SCALE ( TEMPER, PRESUR, WLRC0, XLRC0, RCUT )       **
C    **    SCALES MOMENTA AND BOX SIZE TO GIVE DESIRED T AND P        **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        INTEGER     STEP, NSTEP, IPRINT, ISCALE
        REAL        ACV, ACK, ACE, ACP, ACT, ACH, ACD
        REAL        ACVSQ, ACKSQ, ACESQ, ACPSQ, ACTSQ, ACHSQ, ACDSQ
        REAL        AVV, AVK, AVE, AVP, AVT, AVH, AVD
        REAL        FLV, FLK, FLE, FLP, FLT, FLH, FLD
        REAL        DT, DENS, TEMP, RCUT, VOL, PRES, H, NORM
        REAL        TEMPER, PRESUR
        REAL        K, V, W, E
        REAL        KN, VN, EN, HN
        REAL        X, RPX, PF, PP, PV, CHI, CHIPXI
        REAL        SR3, SR9, VLRC, WLRC, XLRC, VLRC0, WLRC0, XLRC0, PI
        CHARACTER   TITLE*80, CNFILE*30
        REAL        FREE

        PARAMETER ( FREE = 3.0 )
        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM EVAMOR ****                  '')')
        WRITE(*,'(//1X,''MOLECULAR DYNAMICS OF LENNARD-JONES ATOMS '')')
        WRITE(*,'(1X,''CONSTANT-NPT ALGORITHM OF EVANS AND MORRISS '')')

C    ** BASIC SIMULATION PARAMETERS **

        WRITE(*,'('' ENTER RUN TITLE                               '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER NUMBER OF STEPS                         '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER INTERVAL BETWEEN PRINTS                 '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER INTERVAL BETWEEN BOX AND MOMENTUM SCALES'')')
        READ (*,*) ISCALE
        WRITE(*,'('' ENTER THE FOLLOWING IN L-J REDUCED UNITS      '')')
        WRITE(*,'('' ENTER TIMESTEP                                '')')
        READ (*,*) DT
        WRITE(*,'('' ENTER POTENTIAL CUTOFF                        '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER DESIRED TEMPERATURE                     '')')
        READ (*,*) TEMPER
        WRITE(*,'('' ENTER DESIRED PRESSURE                        '')')
        READ (*,*) PRESUR

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE
        WRITE(*,'('' NUMBER OF STEPS  = '',I6   )') NSTEP
        WRITE(*,'('' PRINT INTERVAL   = '',I6   )') IPRINT
        WRITE(*,'('' SCALE INTERVAL   = '',I6   )') ISCALE
        WRITE(*,'('' TIMESTEP         = '',F10.5)') DT
        WRITE(*,'('' POTENTIAL CUTOFF = '',F10.5)') RCUT
        WRITE(*,'('' DESIRED TEMP.    = '',F10.5)') TEMPER
        WRITE(*,'('' DESIRED PRES.    = '',F10.5)') PRESUR

C    ** READCN MUST READ IN INITIAL CONFIGURATION    **
C    ** AND ASSIGN VALUES TO BOX AND ITS DERIVATIVES **

        CALL READCN ( CNFILE )
        VOL = BOX ** 3
        DENS = REAL ( N ) / VOL

        WRITE(*,'('' INITIAL DENS.    = '',F10.5)') DENS

        IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1
        IF ( ISCALE .LE. 0 ) ISCALE = NSTEP + 1

C    ** PREPARE FACTORS FOR LONG-RANGE CORRECTIONS **
C    ** NB: SPECIFIC TO LENNARD-JONES POTENTIAL    **

        SR3 = ( 1.0 / RCUT ) ** 3
        SR9 = SR3 ** 3
        VLRC0 =   8.0 * PI * REAL ( N ) * SR9 / 9.0
     :          - 8.0 * PI * REAL ( N ) * SR3 / 3.0
        WLRC0 =  32.0 * PI * REAL ( N ) * SR9 / 9.0
     :         - 16.0 * PI * REAL ( N ) * SR3 / 3.0
        XLRC0 = 128.0 * PI * REAL ( N ) * SR9 / 9.0
     :         - 32.0 * PI * REAL ( N ) * SR3 / 3.0

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

           CALL FORCE ( RCUT, V, W, X, RPX, PF )
           CALL KINET ( K )

C       ** INCLUDE LONG-RANGE CORRECTIONS IN ALGORITHM **

           KN   = K / REAL ( N )
           VOL  = BOX ** 3
           DENS = REAL ( N ) / VOL
           WLRC = WLRC0 * DENS
           XLRC = XLRC0 * DENS
           TEMP = 2.0 * KN / FREE
           PRES = DENS * TEMP + ( W + WLRC ) / VOL
           PV   = PRES * VOL
           PP   = 2.0 * K

           CHI    = - RPX / 9.0 / ( X + XLRC + PV )
           CHIPXI = PF / PP

           CALL CORREC ( DT, CHI, CHIPXI )

C       ** CALCULATE INSTANTANEOUS VALUES   **
C       ** INCLUDING LONG-RANGE CORRECTIONS **

           VOL  = BOX ** 3
           DENS = REAL ( N ) / VOL
           VLRC = VLRC0 * DENS
           WLRC = WLRC0 * DENS
           XLRC = XLRC0 * DENS
           V    = V + VLRC
           W    = W + WLRC
           X    = X + XLRC
           E    = K + V
           VN   = V / REAL ( N )
           EN   = E / REAL ( N )
           TEMP = 2.0 * KN / FREE
           PRES = DENS * TEMP + W / VOL
           H    = E + PRES * VOL
           HN   = H / REAL ( N )

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

              WRITE(*,'(1X,I8,7(2X,F10.5))')
     :                 STEP, EN, HN, KN, VN, PRES, TEMP, DENS

           ENDIF

C       ** OPTIONALLY SCALE MOMENTA AND BOX SIZE **

           IF ( MOD ( STEP, ISCALE ) .EQ. 0 ) THEN

              CALL SCALE ( TEMPER, PRESUR, WLRC0, XLRC0, RCUT )

           ENDIF

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

C    ** WRITE OUT FINAL CONFIGURATION **
C    ** INCLUDING BOX,BOX1,BOX2,BOX3  **

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
     :          '  ..POTENT..  .PRESSURE.  ..TEMPER..  ..DENSITY.')
        END



        SUBROUTINE FORCE ( RCUT, V, W, X, RPX, PF )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

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
C    ** REAL    PX(N),PY(N),PZ(N) MOLECULAR MOMENTA                   **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES                    **
C    ** REAL    BOX               SIMULATION BOX LENGTH               **
C    ** REAL    RCUT              PAIR POTENTIAL CUTOFF               **
C    ** REAL    V                 POTENTIAL ENERGY                    **
C    ** REAL    W                 VIRIAL FUNCTION                     **
C    ** REAL    X                 HYPERVIRIAL FUNCTION                **
C    ** REAL    RPX               ADDITIONAL PAIR SUMMATION           **
C    ** REAL    PF                ANOTHER PAIR SUMMATION              **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        INTEGER     I, J
        REAL        RCUT, V, W, X, RPX, PF
        REAL        BOXINV, RCUTSQ
        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, RIJSQ
        REAL        PXI, PYI, PZI, PXIJ, PYIJ, PZIJ
        REAL        FXI, FYI, FZI, FXIJ, FYIJ, FZIJ
        REAL        SR2, SR6, SR12, VIJ, WIJ, FIJ, XIJ, RP

C    *******************************************************************

C    ** USEFUL QUANTITIES **

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
        X   = 0.0
        RPX = 0.0
        PF  = 0.0

C    ** OUTER LOOP BEGINS **

        DO 200 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           FXI = FX(I)
           FYI = FY(I)
           FZI = FZ(I)
           PXI = PX(I)
           PYI = PY(I)
           PZI = PZ(I)

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
                 SR12  = SR6 ** 2
                 VIJ   = SR12 - SR6
                 V     = V + VIJ
                 WIJ   = VIJ + SR12
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
                 XIJ   = WIJ + 2.0 * SR12
                 X     = X + XIJ

                 PXIJ  = PXI - PX(J)
                 PYIJ  = PYI - PY(J)
                 PZIJ  = PZI - PZ(J)
                 RP    = RXIJ * PXIJ + RYIJ * PYIJ + RZIJ * PZIJ
                 RPX   = RPX + RP * XIJ * SR2
                 PF    = PF + PXIJ * FXIJ + PYIJ * FYIJ + PZIJ * FZIJ

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

           FX(I) = FX(I) * 24.0
           FY(I) = FY(I) * 24.0
           FZ(I) = FZ(I) * 24.0

300     CONTINUE

        V   = V * 4.0
        W   = W * 24.0 / 3.0
        X   = X * 144.0 / 9.0
        RPX = RPX * 144.0
        PF  = PF * 24.0

        RETURN
        END



        SUBROUTINE KINET ( K )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** ROUTINE TO COMPUTE KINETIC ENERGY                             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        REAL        K
        INTEGER     I

C    *******************************************************************

        K = 0.0

        DO 1000 I = 1, N

           K = K + PX(I) ** 2 + PY(I) ** 2 + PZ(I) ** 2

1000    CONTINUE

        K = 0.5 * K

        RETURN
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT, NN
        PARAMETER ( CNUNIT = 10 )

C     ******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) NN, BOX, BOX1, BOX2, BOX3
        IF ( NN .NE. N ) STOP 'INCORRECT VALUE OF N'
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) RX1, RY1, RZ1
        READ ( CNUNIT ) RX2, RY2, RZ2
        READ ( CNUNIT ) RX3, RY3, RZ3
        READ ( CNUNIT ) PX, PY, PZ
        READ ( CNUNIT ) PX1, PY1, PZ1
        READ ( CNUNIT ) PX2, PY2, PZ2
        READ ( CNUNIT ) PX3, PY3, PZ3

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION TO UNIT 10           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N, BOX, BOX1, BOX2, BOX3
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) RX1, RY1, RZ1
        WRITE ( CNUNIT ) RX2, RY2, RZ2
        WRITE ( CNUNIT ) RX3, RY3, RZ3
        WRITE ( CNUNIT ) PX, PY, PZ
        WRITE ( CNUNIT ) PX1, PY1, PZ1
        WRITE ( CNUNIT ) PX2, PY2, PZ2
        WRITE ( CNUNIT ) PX3, PY3, PZ3

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE PREDIC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** STANDARD TAYLOR SERIES PREDICTORS                             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3
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

           PX(I)  = PX(I)  + C1*PX1(I) + C2*PX2(I) + C3*PX3(I)
           PY(I)  = PY(I)  + C1*PY1(I) + C2*PY2(I) + C3*PY3(I)
           PZ(I)  = PZ(I)  + C1*PZ1(I) + C2*PZ2(I) + C3*PZ3(I)
           PX1(I) = PX1(I) + C1*PX2(I) + C2*PX3(I)
           PY1(I) = PY1(I) + C1*PY2(I) + C2*PY3(I)
           PZ1(I) = PZ1(I) + C1*PZ2(I) + C2*PZ3(I)
           PX2(I) = PX2(I) + C1*PX3(I)
           PY2(I) = PY2(I) + C1*PY3(I)
           PZ2(I) = PZ2(I) + C1*PZ3(I)

100     CONTINUE

        BOX  = BOX  + C1*BOX1 + C2*BOX2 + C3*BOX3
        BOX1 = BOX1 + C1*BOX2 + C2*BOX3
        BOX2 = BOX2 + C1*BOX3

        RETURN
        END



        SUBROUTINE CORREC ( DT, CHI, CHIPXI )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** GEAR CORRECTOR ALGORITHM.                                     **
C    **                                                               **
C    ** FOR TIMESTEP-SCALED VARIABLES GEAR COEFFICIENTS WOULD BE AS   **
C    ** FOLLOWS (4-VALUE METHOD, 1ST-ORDER D.E.): 3/8, 1, 3/4, 1/6.   **
C    ** CHIPXI IS SHORT FOR ( CHI + XI ) (SEE CHAPTER 7).             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3
        REAL        DT, CHI, CHIPXI

        INTEGER     I
        REAL        C1, C2, C3, COEFF0, COEFF2, COEFF3
        REAL        CORR, CORRX, CORRY, CORRZ, CORPX, CORPY, CORPZ
        REAL        RX1I, RY1I, RZ1I, PX1I, PY1I, PZ1I

        REAL        GEAR0, GEAR2, GEAR3
        PARAMETER ( GEAR0 = 3.0 / 8.0,
     :              GEAR2 = 3.0 / 4.0,
     :              GEAR3 = 1.0 / 6.0 )

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0

        COEFF0 = GEAR0 * C1
        COEFF2 = GEAR2 * C1 / C2
        COEFF3 = GEAR3 * C1 / C3

        DO 400 I = 1, N

           RX1I = PX(I) + CHI * RX(I)
           RY1I = PY(I) + CHI * RY(I)
           RZ1I = PZ(I) + CHI * RZ(I)
           CORRX = RX1I - RX1(I)
           CORRY = RY1I - RY1(I)
           CORRZ = RZ1I - RZ1(I)

           RX(I)  = RX(I)  + COEFF0 * CORRX
           RY(I)  = RY(I)  + COEFF0 * CORRY
           RZ(I)  = RZ(I)  + COEFF0 * CORRZ
           RX1(I) = RX1I
           RY1(I) = RY1I
           RZ1(I) = RZ1I
           RX2(I) = RX2(I) + COEFF2 * CORRX
           RY2(I) = RY2(I) + COEFF2 * CORRY
           RZ2(I) = RZ2(I) + COEFF2 * CORRZ
           RX3(I) = RX3(I) + COEFF3 * CORRX
           RY3(I) = RY3(I) + COEFF3 * CORRY
           RZ3(I) = RZ3(I) + COEFF3 * CORRZ

           PX1I = FX(I) - CHIPXI * PX(I)
           PY1I = FY(I) - CHIPXI * PY(I)
           PZ1I = FZ(I) - CHIPXI * PZ(I)
           CORPX = PX1I - PX1(I)
           CORPY = PY1I - PY1(I)
           CORPZ = PZ1I - PZ1(I)

           PX(I)  = PX(I)  + COEFF0 * CORPX
           PY(I)  = PY(I)  + COEFF0 * CORPY
           PZ(I)  = PZ(I)  + COEFF0 * CORPZ
           PX1(I) = PX1I
           PY1(I) = PY1I
           PZ1(I) = PZ1I
           PX2(I) = PX2(I) + COEFF2 * CORPX
           PY2(I) = PY2(I) + COEFF2 * CORPY
           PZ2(I) = PZ2(I) + COEFF2 * CORPZ
           PX3(I) = PX3(I) + COEFF3 * CORPX
           PY3(I) = PY3(I) + COEFF3 * CORPY
           PZ3(I) = PZ3(I) + COEFF3 * CORPZ

400     CONTINUE

        CORR = CHI * BOX - BOX1
        BOX  = BOX  + COEFF0 * CORR
        BOX1 = CHI * BOX
        BOX2 = BOX2 + COEFF2 * CORR
        BOX3 = BOX3 + COEFF3 * CORR

        RETURN
        END



        SUBROUTINE SCALE ( TEMPER, PRESUR, WLRC0, XLRC0, RCUT )

        COMMON / BLOCK1 / RX, RY, RZ, RX1, RY1, RZ1,
     :                    RX2, RY2, RZ2, RX3, RY3, RZ3,
     :                    PX, PY, PZ, PX1, PY1, PZ1,
     :                    PX2, PY2, PZ2, PX3, PY3, PZ3,
     :                    FX, FY, FZ
        COMMON / BLOCK2 / BOX, BOX1, BOX2, BOX3

C    *******************************************************************
C    ** SCALES MOMENTA AND BOX SIZE TO GIVE DESIRED VALUES OF T AND P.**
C    **                                                               **
C    ** WE USE DIRECT MOMENTUM SCALING TO GIVE THE TEMPERATURE AND A  **
C    ** NEWTON-RAPHSON METHOD FOR THE BOX SIZE.                       **
C    ** WE USE THE PRESSURE DERIVATIVE FUNCTION                       **
C    ** BOX * DP/D(BOX) = 3 * VOL * DP/D(VOL) = 3 * ( - P - X / VOL ) **
C    ** WHERE P IS THE PRESSURE AND X THE HYPERVIRIAL FUNCTION.       **
C    ** LONG-RANGE CORRECTIONS ARE TAKEN INTO ACCOUNT.                **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RX1(N), RY1(N), RZ1(N)
        REAL        RX2(N), RY2(N), RZ2(N)
        REAL        RX3(N), RY3(N), RZ3(N)
        REAL        PX(N), PY(N), PZ(N)
        REAL        PX1(N), PY1(N), PZ1(N)
        REAL        PX2(N), PY2(N), PZ2(N)
        REAL        PX3(N), PY3(N), PZ3(N)
        REAL        FX(N), FY(N), FZ(N)
        REAL        BOX, BOX1, BOX2, BOX3

        INTEGER     I
        REAL        TEMPER, PRESUR, WLRC0, XLRC0, RCUT
        REAL        K, KN, V, W, X, RPX, PF, FACTOR
        REAL        VOL, DENS, WLRC, XLRC, TEMP, PRES
        REAL        FREE, TOL
        PARAMETER ( FREE = 3.0, TOL = 0.01 )

C    *******************************************************************

C    ** SCALE MOMENTA **

        CALL KINET ( K )
        KN = K / REAL ( N )
        TEMP = 2.0 * KN / FREE
        FACTOR = SQRT ( TEMPER / TEMP )

        DO 100 I = 1, N

           PX(I) = PX(I) * FACTOR
           PY(I) = PY(I) * FACTOR
           PZ(I) = PZ(I) * FACTOR

100     CONTINUE

        TEMP = TEMPER

C    ** SCALE PRESSURE USING NEWTON-RAPHSON PROCEDURE **

        CALL FORCE ( RCUT, V, W, X, RPX, PF )
        VOL  = BOX ** 3
        DENS = REAL ( N ) / VOL
        WLRC = WLRC0 * DENS
        XLRC = XLRC0 * DENS
        W    = W + WLRC
        X    = X + XLRC
        PRES = DENS * TEMP + W / VOL

1000    IF ( ABS ( PRES - PRESUR ) .GT. TOL ) THEN

           FACTOR = 1.0 + ( PRES - PRESUR ) / ( PRES + X / VOL ) / 3.0

           DO 200 I = 1, N

              RX(I) = RX(I) * FACTOR
              RY(I) = RY(I) * FACTOR
              RZ(I) = RZ(I) * FACTOR

200        CONTINUE

           BOX = BOX * FACTOR

           CALL FORCE ( RCUT, V, W, X, RPX, PF )
           VOL  = BOX ** 3
           DENS = REAL ( N ) / VOL
           WLRC = WLRC0 * DENS
           XLRC = XLRC0 * DENS
           W    = W + WLRC
           X    = X + XLRC
           PRES = DENS * TEMP + W / VOL

           GOTO 1000

        ENDIF

        RETURN
        END



