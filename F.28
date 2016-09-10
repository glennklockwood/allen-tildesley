********************************************************************************
** FICHE F.28.  CONSTANT-NVT MOLECULAR DYNAMICS - EXTENDED SYSTEM METHOD      **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM NOSE

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** CONSTANT-NVT MOLECULAR DYNAMICS USING AN EXTENDED SYSTEM.     **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** NOSE, MOLEC. PHYS. 52, 255, 1984.                             **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE FORCE ( BOX, RCUT, V, VC, W )                      **
C    **    CALCULATES ACCELERATIONS, POTENTIAL AND VIRIAL             **
C    ** SUBROUTINE KINET ( K )                                        **
C    **    CALCULATES KINETIC ENERGY                                  **
C    ** SUBROUTINE READCN ( CNFILE, BOX )                             **
C    **    READS CONFIGURATION                                        **
C    ** SUBROUTINE WRITCN ( CNFILE, BOX )                             **
C    **    WRITES CONFIGURATION                                       **
C    ** SUBROUTINE PREDIC ( DT )                                      **
C    **    PREDICTS POSITIONS AT TIME T + DT                          **
C    ** SUBROUTINE CORREC ( DT )                                      **
C    **    CORRECTS POSITIONS AT TIME T + DT                          **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF ATOMS         **
C    ** REAL    DT                            TIMESTEP                **
C    ** REAL    RX(N),RY(N),RZ(N)             ATOM POSITIONS          **
C    ** REAL    VX(N),VY(N),VZ(N)             FIRST DERIVATIVES       **
C    ** REAL    AX(N),AY(N),AZ(N)             SECOND DERIVATIVES      **
C    ** REAL    BX(N),BY(N),BZ(N)             THIRD DERIVATIVES       **
C    ** REAL    FX(N),FY(N),FZ(N)             TOTAL FORCES            **
C    ** REAL    S,SV,SA,SB,SC,SF              S AND DERIVATIVES       **
C    ** REAL    Q                             THERMAL INERTIA         **
C    ** REAL    ACV,ACK ETC.                  AV VALUE ACCUMULATORS   **
C    ** REAL    AVV,AVK ETC.                  AV VALUES               **
C    ** REAL    ACVSQ, ACKSQ ETC.             AV**2 VAL ACCUMULATORS  **
C    ** REAL    FLV, FLK ETC.                 FLUCT AVERAGES          **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE MODIFIED EQUATIONS OF MOTION ARE AS FOLLOWS:              **
C    ** A = F/(M*S**2) - 2*SV*V/S                                     **
C    ** SA = (SUM (M*S/Q)*V**2 ) - ((3N-3)+1)*KT/(Q*S)                **
C    ** WHERE R,V,A,.. ARE SUCCESSIVE DERIVATIVES OF POSITION         **
C    ** AND SV,SA,.... ARE SUCCESSIVE DERIVATIVES OF EXTRA VARIABLE S **
C    ** F REPRESENTS FORCES, M PARTICLE MASS, Q EXTRA VARIABLE MASS   **
C    ** KT IS TEMPERATURE*BOLTZMANN'S CONSTANT                        **
C    ** AND WE HAVE (3N-3) DEGREES OF FREEDOM IF MOMENTUM FIXED.      **
C    ** IN A SENSE, S IS A TIME-SCALING FACTOR.                       **
C    ** WE SOLVE THESE EQUATIONS BY A GEAR 5-VALUE METHOD FOR SECOND  **
C    ** ORDER DIFFERENTIAL EQUATIONS.                                 **
C    ** THIS PROGRAM USES UNIT 10 FOR CONFIGURATION INPUT AND OUTPUT  **
C    *******************************************************************

        INTEGER       N
        PARAMETER   ( N = 108 )

        REAL          RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL          AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL          CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL          S, SV, SA, SB, SC, SF

        INTEGER       STEP, NSTEP, IPRINT
        REAL          Q, QK, QV
        REAL          ACV, ACK, ACE, ACEC, ACP, ACT, ACH, ACS
        REAL          AVV, AVK, AVE, AVEC, AVP, AVT, AVH, AVS
        REAL          ACVSQ, ACKSQ, ACESQ, ACECSQ, ACPSQ, ACTSQ
        REAL          ACHSQ, ACSSQ
        REAL          FLV, FLK, FLE, FLEC, FLP, FLT, FLH, FLS
        REAL          DT, DENS, TEMPER, RCUT, BOX, PRES, NORM, TEMP
        REAL          K, V, VC, W, E, EC, FREE, FREEN, H, VOL
        REAL          KN, VN, EN, ECN, HN
        REAL          SR3, SR9, VLRC, WLRC, PI
        CHARACTER     TITLE*80, CNFILE*30

        PARAMETER   ( PI = 3.1415927 )
        PARAMETER   ( FREE = 3.0 )

C    *******************************************************************

C    ** READ IN INITIAL PARAMETERS **

        WRITE(*,'(1H1,'' **** PROGRAM NOSE ****                    '')')
        WRITE(*,'(//1X,''MOLECULAR DYNAMICS OF LENNARD-JONES ATOMS '')')
        WRITE(*,'('' CONSTANT-NVT DYNAMICS BY THE METHOD OF NOSE   '')')
        WRITE(*,'('' ENTER RUN TITLE                               '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF STEPS                         '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER INTERVAL BETWEEN PRINTS                 '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE FOLLOWING IN L-J REDUCED UNITS      '')')
        WRITE(*,'('' ENTER TIMESTEP                                '')')
        READ (*,*) DT
        WRITE(*,'('' ENTER TEMPERATURE                             '')')
        READ (*,*) TEMPER
        WRITE(*,'('' ENTER POTENTIAL CUTOFF                        '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER THERMAL INERTIA PARAMETER               '')')
        READ (*,*) Q

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' NUMBER OF ATOMS  = '',I6   )') N
        WRITE(*,'('' NUMBER OF STEPS  = '',I6   )') NSTEP
        WRITE(*,'('' PRINT INTERVAL   = '',I6   )') IPRINT
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE
        WRITE(*,'('' TIMESTEP         = '',F10.5)') DT
        WRITE(*,'('' TEMPERATURE      = '',F10.5)') TEMPER
        WRITE(*,'('' POTENTIAL CUTOFF = '',F10.5)') RCUT
        WRITE(*,'('' THERMAL INERTIA  = '',F10.2)') Q

        FREEN = FREE * REAL ( N )

C    ** READCN MUST READ IN INITIAL CONFIGURATION  **
C    ** AND ASSIGN VALUES TO S AND ITS DERIVATIVES **

        CALL READCN ( CNFILE, BOX )

        VOL  = BOX ** 3
        DENS = REAL ( N ) / VOL

        WRITE(*,'('' DENSITY          = '',F10.5)') DENS

C    ** CALCULATE LONG-RANGE CORRECTIONS **
C    ** NB: SPECIFIC TO LENNARD-JONES    **

        SR3  = ( 1.0 / RCUT ) ** 3
        SR9  = SR3 ** 3
        VLRC = ( 8.0 / 9.0  ) * PI * DENS * REAL ( N ) *
     :         ( SR9 - 3.0 * SR3 )
        WLRC = ( 16.0 / 9.0 ) * PI * DENS * REAL ( N ) *
     :         ( 2.0 * SR9 - 3.0 * SR3 )

C    ** ZERO ACCUMULATORS **

        ACV  = 0.0
        ACK  = 0.0
        ACE  = 0.0
        ACEC = 0.0
        ACP  = 0.0
        ACT  = 0.0
        ACH  = 0.0
        ACS  = 0.0

        ACVSQ  = 0.0
        ACKSQ  = 0.0
        ACESQ  = 0.0
        ACECSQ = 0.0
        ACPSQ  = 0.0
        ACTSQ  = 0.0
        ACHSQ  = 0.0
        ACSSQ  = 0.0

        FLV  = 0.0
        FLK  = 0.0
        FLE  = 0.0
        FLEC = 0.0
        FLP  = 0.0
        FLT  = 0.0
        FLH  = 0.0
        FLS  = 0.0

        IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1

        WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
        WRITE(*,10001)

C    *******************************************************************
C    ** MAIN LOOP BEGINS                                              **
C    *******************************************************************

        DO 1000 STEP = 1, NSTEP

C       ** IMPLEMENT ALGORITHM **

           CALL PREDIC ( DT )
           CALL FORCE  ( BOX, RCUT, V, VC, W )
           CALL KINET  ( K )

           SF = ( 2.0 * K - ( FREEN + 1.0 ) * TEMPER ) / ( S * Q )

           CALL CORREC ( DT )

           QK   = 0.5 * Q * ( SV ** 2 )
           QV   = ( FREEN + 1.0 ) * TEMPER * LOG ( S )
           V    = V + VLRC
           W    = W + WLRC
           E    = K + V
           EC   = K + VC
           H    = K + VC + QK + QV
           EN   = E / REAL ( N )
           VN   = V / REAL ( N )
           ECN  = EC / REAL ( N )
           HN   = H / REAL ( N )
           TEMP = KN * 2.0 / FREE
           PRES = DENS * TEMP + W / VOL

C       ** INCREMENT ACCUMULATORS **

           ACE  = ACE  + EN
           ACEC = ACEC + ECN
           ACK  = ACK  + KN
           ACV  = ACV  + VN
           ACP  = ACP  + PRES
           ACH  = ACH  + HN
           ACS  = ACS  + S

           ACESQ  = ACESQ  + EN ** 2
           ACECSQ = ACECSQ + ECN ** 2
           ACKSQ  = ACKSQ  + KN ** 2
           ACVSQ  = ACVSQ  + VN ** 2
           ACPSQ  = ACPSQ  + PRES ** 2
           ACHSQ  = ACHSQ  + HN ** 2
           ACSSQ  = ACSSQ  + S ** 2

C       ** OPTIONALLY PRINT INFORMATION **

           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN

              WRITE(*,'(1X,I8,8(2X,F10.4))')
     :           STEP, HN, EN, ECN, KN, VN, TEMP, PRES, S

           ENDIF

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

C    ** WRITE OUT FINAL AVERAGES **

        NORM = REAL ( NSTEP )

        AVE  = ACE  / NORM
        AVEC = ACEC / NORM
        AVK  = ACK  / NORM
        AVV  = ACV  / NORM
        AVP  = ACP  / NORM
        AVH  = ACH  / NORM
        AVS  = ACS  / NORM

        ACESQ  = ( ACESQ  / NORM ) - AVE  ** 2
        ACECSQ = ( ACECSQ / NORM ) - AVEC ** 2
        ACKSQ  = ( ACKSQ  / NORM ) - AVK  ** 2
        ACVSQ  = ( ACVSQ  / NORM ) - AVV  ** 2
        ACPSQ  = ( ACPSQ  / NORM ) - AVP  ** 2
        ACHSQ  = ( ACHSQ  / NORM ) - AVH  ** 2
        ACSSQ  = ( ACSSQ  / NORM ) - AVS  ** 2

        IF ( ACESQ  .GT. 0.0 ) FLE  = SQRT ( ACESQ  )
        IF ( ACECSQ .GT. 0.0 ) FLEC = SQRT ( ACECSQ )
        IF ( ACKSQ  .GT. 0.0 ) FLK  = SQRT ( ACKSQ  )
        IF ( ACVSQ  .GT. 0.0 ) FLV  = SQRT ( ACVSQ  )
        IF ( ACPSQ  .GT. 0.0 ) FLP  = SQRT ( ACPSQ  )
        IF ( ACHSQ  .GT. 0.0 ) FLH  = SQRT ( ACHSQ  )
        IF ( ACSSQ  .GT. 0.0 ) FLS  = SQRT ( ACSSQ  )

        AVT = AVK * 2.0 / FREE
        FLT = FLK * 2.0 / FREE

        WRITE(*,'('' AVERAGES'',8(2X,F10.5))')
     :     AVH, AVE, AVEC, AVK, AVV, AVT, AVP, AVS
        WRITE(*,'('' FLUCTS  '',8(2X,F10.5))')
     :     FLH, FLE, FLEC, FLK, FLV, FLT, FLP, FLS

C    ** WRITE OUT FINAL CONFIGURATION **

        CALL WRITCN ( CNFILE, BOX )

        STOP

10001   FORMAT(//1X,'TIMESTEP    ...HAMIL..  ..ENERGY..',
     :              '  CUTENERGY.  ..KINETIC.  ..POTENT..',
     :              '  ..TEMPER..  .PRESSURE.  ....S.....'/)
        END



        SUBROUTINE FORCE ( BOX, RCUT, V, VC, W )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ

C    *******************************************************************
C    ** LENNARD-JONES FORCE ROUTINE IN REDUCED UNITS                  **
C    **                                                               **
C    ** THE POTENTIAL IS V(R) = 4*((1/R)**12-(1/R)**6)                **
C    ** WE INCLUDE SPHERICAL CUTOFF AND MINIMUM IMAGING IN CUBIC BOX. **
C    ** TWO POTENTIAL ENERGIES ARE RETURNED.                          **
C    ** V IS CALCULATED USING THE UNSHIFTED POTENTIAL. WHEN LONG-     **
C    ** RANGE TAIL CORRECTIONS ARE ADDED, THIS MAY BE USED TO         **
C    ** CALCULATE THERMODYNAMIC INTERNAL ENERGY ETC.                  **
C    ** VC IS CALCULATED USING THE SHIFTED POTENTIAL WITH NO          **
C    ** DISCONTINUITY AT THE CUTOFF.  THIS MAY BE USED TO CHECK THE   **
C    ** CONSERVATION LAWS.                                            **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES                    **
C    ** REAL    BOX               SIMULATION BOX LENGTH               **
C    ** REAL    RCUT              PAIR POTENTIAL CUTOFF               **
C    ** REAL    V                 POTENTIAL ENERGY                    **
C    ** REAL    VC                SHIFTED POTENTIAL                   **
C    ** REAL    W                 VIRIAL FUNCTION                     **
C    ** REAL    VIJ               PAIR POTENTIAL BETWEEN I AND J      **
C    ** REAL    WIJ               NEGATIVE OF PAIR VIRIAL FUNCTION W  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RCUT, BOX, V, VC, W
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)

        INTEGER     I, J, NCUT
        REAL        BOXINV, RCUTSQ
        REAL        RXI, RYI, RZI, FXI, FYI, FZI
        REAL        RXIJ, RYIJ, RZIJ, FXIJ, FYIJ, FZIJ
        REAL        RIJSQ, SR2, SR6, SR12, VIJ, WIJ, FIJ

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

        SR2  = 1.0 / RCUTSQ
        SR6  = SR2 * SR2 * SR2
        SR12 = SR6 * SR6
        VIJ  = SR12 - SR6
        VC   = V - REAL ( NCUT ) * VIJ

C    ** MULTIPLY RESULTS BY NUMERICAL FACTORS **

        DO 300 I = 1, N

           FX(I) = FX(I) * 24.0
           FY(I) = FY(I) * 24.0
           FZ(I) = FZ(I) * 24.0

300     CONTINUE

        V  = V  * 4.0
        VC = VC * 4.0
        W  = W  * 24.0 / 3.0

        RETURN
        END



        SUBROUTINE KINET ( K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** ROUTINE TO COMPUTE KINETIC ENERGY IN NOSE METHOD              **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        K

        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL        S, SV, SA, SB, SC, SF

        INTEGER     I

C    *******************************************************************

        K = 0.0

        DO 1000 I = 1, N

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

1000    CONTINUE

        K = 0.5 * ( S ** 2 ) * K

        RETURN
        END




        SUBROUTINE READCN ( CNFILE, BOX )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        CHARACTER   CNFILE*(*)
        REAL        BOX
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL        S, SV, SA, SB, SC, SF

        INTEGER     CNUNIT, NN
        PARAMETER ( CNUNIT = 10 )

C     ******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS= 'OLD', FORM = 'UNFORMATTED' )

        READ  ( CNUNIT ) NN, BOX, S, SV, SA, SB, SC
        IF ( NN .NE. N ) STOP ' INCORRECT NUMBER OF ATOMS '
        READ  ( CNUNIT ) RX, RY, RZ
        READ  ( CNUNIT ) VX, VY, VZ
        READ  ( CNUNIT ) AX, AY, AZ
        READ  ( CNUNIT ) BX, BY, BZ
        READ  ( CNUNIT ) CX, CY, CZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE, BOX )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION                      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        CHARACTER   CNFILE*(*)
        REAL        BOX
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL        S, SV, SA, SB, SC, SF

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C    ****************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS='OLD', FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N, BOX, S, SV, SA, SB, SC
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) VX, VY, VZ
        WRITE ( CNUNIT ) AX, AY, AZ
        WRITE ( CNUNIT ) BX, BY, BZ
        WRITE ( CNUNIT ) CX, CY, CZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE PREDIC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** STANDARD TAYLOR SERIES PREDICTORS                             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL        S, SV, SA, SB, SC, SF
        REAL        DT

        INTEGER     I
        REAL        C1, C2, C3, C4

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        DO 100 I = 1, N

           RX(I) = RX(I) + C1*VX(I) + C2*AX(I) + C3*BX(I) + C4*CX(I)
           RY(I) = RY(I) + C1*VY(I) + C2*AY(I) + C3*BY(I) + C4*CY(I)
           RZ(I) = RZ(I) + C1*VZ(I) + C2*AZ(I) + C3*BZ(I) + C4*CZ(I)
           VX(I) = VX(I) + C1*AX(I) + C2*BX(I) + C3*CX(I)
           VY(I) = VY(I) + C1*AY(I) + C2*BY(I) + C3*CY(I)
           VZ(I) = VZ(I) + C1*AZ(I) + C2*BZ(I) + C3*CZ(I)
           AX(I) = AX(I) + C1*BX(I) + C2*CX(I)
           AY(I) = AY(I) + C1*BY(I) + C2*CY(I)
           AZ(I) = AZ(I) + C1*BZ(I) + C2*CZ(I)
           BX(I) = BX(I) + C1*CX(I)
           BY(I) = BY(I) + C1*CY(I)
           BZ(I) = BZ(I) + C1*CZ(I)

100     CONTINUE

        S  = S  + C1*SV + C2*SA + C3*SB + C4*SC
        SV = SV + C1*SA + C2*SB + C3*SC
        SA = SA + C1*SB + C2*SC
        SB = SB + C1*SC

        RETURN
        END



        SUBROUTINE CORREC ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, AX, AY, AZ,
     :                    BX, BY, BZ, CX, CY, CZ, FX, FY, FZ
        COMMON / BLOCK2 / S, SV, SA, SB, SC, SF

C    *******************************************************************
C    ** GEAR CORRECTOR ALGORITHM                                      **
C    **                                                               **
C    ** WE USE GEAR'S ALTERNATIVE COEFFICIENT GEAR0 WHICH APPLIES IN  **
C    ** THE CASE THAT FIRST TIME DERIVATIVES APPEAR ON THE RIGHT OF   **
C    ** THESE SECOND-ORDER DIFFERENTIAL EQUATIONS.                    **
C    ** FOR TIMESTEP-SCALED VARIABLES THE COEFFICIENTS WOULD BE       **
C    ** 19/90, 3/4, 1, 1/2, 1/12.                                     **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** GEAR, NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY            **
C    ** DIFFERENTIAL EQUATIONS (PRENTICE-HALL, 1971).                 **
C    ** GEAR, REPORT ANL 7126 (ARGONNE NATIONAL LAB., 1966).          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        AX(N), AY(N), AZ(N), BX(N), BY(N), BZ(N)
        REAL        CX(N), CY(N), CZ(N), FX(N), FY(N), FZ(N)
        REAL        S, SV, SA, SB, SC, SF
        REAL        DT

        INTEGER     I
        REAL        C1, C2, C3, C4
        REAL        CR, CV, CB, CC
        REAL        CORR, CORRX, CORRY, CORRZ
        REAL        AXI, AYI, AZI

        REAL        GEAR0, GEAR1, GEAR3, GEAR4
        PARAMETER ( GEAR0 = 19.0 / 90.0, GEAR1 = 3.0 /4.0,
     :              GEAR3 = 1.0 / 2.0,   GEAR4 = 1.0 /12.0 )

C    *******************************************************************

        C1 = DT
        C2 = C1 * DT / 2.0
        C3 = C2 * DT / 3.0
        C4 = C3 * DT / 4.0

        CR = GEAR0 * C2
        CV = GEAR1 * C2 / C1
        CB = GEAR3 * C2 / C3
        CC = GEAR4 * C2 / C4

        DO 400 I = 1, N

           AXI = FX(I) / ( S ** 2 ) - 2.0 * SV * VX(I) / S
           AYI = FY(I) / ( S ** 2 ) - 2.0 * SV * VY(I) / S
           AZI = FZ(I) / ( S ** 2 ) - 2.0 * SV * VZ(I) / S
           CORRX = AXI - AX(I)
           CORRY = AYI - AY(I)
           CORRZ = AZI - AZ(I)

           RX(I) = RX(I) + CR * CORRX
           RY(I) = RY(I) + CR * CORRY
           RZ(I) = RZ(I) + CR * CORRZ
           VX(I) = VX(I) + CV * CORRX
           VY(I) = VY(I) + CV * CORRY
           VZ(I) = VZ(I) + CV * CORRZ
           AX(I) = AXI
           AY(I) = AYI
           AZ(I) = AZI
           BX(I) = BX(I) + CB * CORRX
           BY(I) = BY(I) + CB * CORRY
           BZ(I) = BZ(I) + CB * CORRZ
           CX(I) = CX(I) + CC * CORRX
           CY(I) = CY(I) + CC * CORRY
           CZ(I) = CZ(I) + CC * CORRZ

400     CONTINUE

        CORR = SF - SA
        S  = S  + CR * CORR
        SV = SV + CV * CORR
        SA = SF
        SB = SB + CB * CORR
        SC = SC + CC * CORR

        RETURN
        END



