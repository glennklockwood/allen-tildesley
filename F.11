********************************************************************************
** FICHE F.11.  CONSTANT-NVT MONTE CARLO FOR LENNARD JONES ATOMS              **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM MCNVT

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** MONTE CARLO SIMULATION PROGRAM IN THE CONSTANT-NVT ENSEMBLE.  **
C    **                                                               **
C    ** THIS PROGRAM TAKES A CONFIGURATION OF LENNARD JONES ATOMS     **
C    ** AND PERFORMS A CONVENTIONAL NVT MC SIMULATION. THE BOX IS OF  **
C    ** UNIT LENGTH, -0.5 TO +0.5 AND THERE ARE NO LOOKUP TABLES.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** INTEGER NSTEP               MAXIMUM NUMBER OF CYCLES          **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    DENS                REDUCED DENSITY                   **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    SIGMA               REDUCED LJ DIAMETER               **
C    ** REAL    RMIN                MINIMUM REDUCED PAIR SEPARATION   **
C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
C    ** REAL    DRMAX               REDUCED MAXIMUM DISPLACEMENT      **
C    ** REAL    V                   THE POTENTIAL ENERGY              **
C    ** REAL    W                   THE VIRIAL                        **
C    ** REAL    PRES                THE PRESSURE                      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE PROGRAM TAKES IN A CONFIGURATION OF ATOMS                 **
C    ** AND RUNS A MONTE CARLO SIMULATION AT THE GIVEN TEMPERATURE    **
C    ** FOR THE SPECIFIED NUMBER OF CYCLES.                           **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM USES LENNARD-JONES UNITS FOR USER INPUT AND       **
C    ** OUTPUT BUT CONDUCTS THE SIMULATION IN A BOX OF UNIT LENGTH.   **
C    ** FOR EXAMPLE, FOR A BOXLENGTH L, AND LENNARD-JONES PARAMETERS  **
C    ** EPSILON AND SIGMA, THE UNITS ARE:                             **
C    **                                                               **
C    **     PROPERTY       LJ  UNITS            PROGRAM UNITS         **
C    **                                                               **
C    **     TEMP           EPSILON/K            EPSILON/K             **
C    **     PRES           EPSILON/SIGMA**3     EPSILON/L**3          **
C    **     V              EPSILON              EPSILON               **
C    **     DENS           1/SIGMA**3           1/L**3                **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, V, W )          **
C    **    CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION  **
C    ** SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, SIGMA, V, W )     **
C    **    CALCULATES THE POTENTIAL ENERGY OF ATOM I WITH ALL THE     **
C    **    OTHER ATOMS IN THE LIQUID                                  **
C    ** SUBROUTINE READCN (CNFILE )                                   **
C    **    READS IN A CONFIGURATION                                   **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT A CONFIGURATION                                 **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)

        REAL        DRMAX, DENS, TEMP, DENSLJ, SIGMA, RMIN, RCUT, BETA
        REAL        RANF, DUMMY, ACM, ACATMA, PI, RATIO, SR9, SR3
        REAL        V, VNEW, VOLD, VEND, VN, DELTV, DELTVB, VS
        REAL        W, WEND, WNEW, WOLD, PRES, DELTW, WS, PS
        REAL        VLRC, VLRC6, VLRC12, WLRC, WLRC6, WLRC12
        REAL        RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW
        REAL        AVV, AVP, AVW, ACV, ACP, ACVSQ, ACPSQ, FLV, FLP
        INTEGER     STEP, I, NSTEP, IPRINT, ISAVE, IRATIO
        LOGICAL     OVRLAP
        CHARACTER   TITLE*80, CNFILE*30

        PARAMETER ( PI = 3.1415927 )

C       ****************************************************************

C    ** READ INPUT DATA **

        WRITE(*,'(1H1,'' **** PROGRAM MCLJ ****                  ''/)')
        WRITE(*,'('' CONSTANT-NVT MONTE CARLO PROGRAM            '' )')
        WRITE(*,'('' FOR LENNARD JONES ATOMS                      '')')

        WRITE(*,'('' ENTER THE RUN TITLE                          '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF CYCLES                       '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT LINES   '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN DATA SAVES     '')')
        READ (*,*) ISAVE
        WRITE(*,'('' ENTER INTERVAL FOR UPDATE OF MAX. DISPL.     '')')
        READ (*,*) IRATIO
        WRITE(*,'('' ENTER THE CONFIGURATION FILE NAME            '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE FOLLOWING IN LENNARD-JONES UNITS '',/)')
        WRITE(*,'('' ENTER THE DENSITY                            '')')
        READ (*,*) DENS
        WRITE(*,'('' ENTER THE TEMPERATURE                        '')')
        READ (*,*) TEMP
        WRITE(*,'('' ENTER THE POTENTIAL CUTOFF DISTANCE          '')')
        READ (*,*) RCUT

C    ** WRITE INPUT DATA **

        WRITE(*,'(       //1X                    ,A     )') TITLE
        WRITE(*,'('' NUMBER OF ATOMS           '',I10   )') N
        WRITE(*,'('' NUMBER OF CYCLES          '',I10   )') NSTEP
        WRITE(*,'('' OUTPUT FREQUENCY          '',I10   )') IPRINT
        WRITE(*,'('' SAVE FREQUENCY            '',I10   )') ISAVE
        WRITE(*,'('' RATIO UPDATE FREQUENCY    '',I10   )') IRATIO
        WRITE(*,'('' CONFIGURATION FILE  NAME  '',A     )') CNFILE
        WRITE(*,'('' TEMPERATURE               '',F10.4 )') TEMP
        WRITE(*,'('' DENSITY                   '',F10.4 )') DENS
        WRITE(*,'('' POTENTIAL CUTOFF          '',F10.4 )') RCUT

C    ** READ INITIAL CONFIGURATION **

        CALL READCN ( CNFILE )

C    ** CONVERT INPUT DATA TO PROGRAM UNITS **

        BETA   = 1.0 / TEMP
        SIGMA  = ( DENS / REAL ( N ) ) ** ( 1.0 / 3.0 )
        RMIN   = 0.70 * SIGMA
        RCUT   = RCUT * SIGMA
        DRMAX  = 0.15 * SIGMA
        DENSLJ = DENS
        DENS   = DENS / ( SIGMA ** 3 )

        IF ( RCUT .GT. 0.5 ) STOP ' CUT-OFF TOO LARGE '

C    ** ZERO ACCUMULATORS **

        ACV    = 0.0
        ACVSQ  = 0.0
        ACP    = 0.0
        ACPSQ  = 0.0
        FLV    = 0.0
        FLP    = 0.0
        ACM    = 0.0
        ACATMA = 0.0

C    ** CALCULATE LONG RANGE CORRECTIONS    **
C    ** SPECIFIC TO THE LENNARD JONES FLUID **

        SR3 = ( SIGMA / RCUT ) ** 3
        SR9 = SR3 ** 3

        VLRC12 =   8.0 * PI * DENSLJ * REAL ( N ) * SR9 / 9.0
        VLRC6  = - 8.0 * PI * DENSLJ * REAL ( N ) * SR3 / 3.0
        VLRC   =   VLRC12 + VLRC6
        WLRC12 =   4.0  * VLRC12
        WLRC6  =   2.0  * VLRC6
        WLRC   =   WLRC12 + WLRC6

C    ** WRITE OUT SOME USEFUL INFORMATION **

        WRITE(*,'('' SIGMA/BOX              =  '',F10.4)')  SIGMA
        WRITE(*,'('' RMIN/BOX               =  '',F10.4)')  RMIN
        WRITE(*,'('' RCUT/BOX               =  '',F10.4)')  RCUT
        WRITE(*,'('' LRC FOR <V>            =  '',F10.4)')  VLRC
        WRITE(*,'('' LRC FOR <W>            =  '',F10.4)')  WLRC

C    ** CALCULATE INITIAL ENERGY AND CHECK FOR OVERLAPS **

        CALL SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, V, W )

        IF ( OVRLAP ) STOP ' OVERLAP IN INITIAL CONFIGURATION '

        VS = ( V + VLRC ) / REAL ( N )
        WS = ( W + WLRC ) / REAL ( N )
        PS = DENS * TEMP + W + WLRC
        PS = PS * SIGMA ** 3

        WRITE(*,'('' INITIAL V              =  '', F10.4 )' ) VS
        WRITE(*,'('' INITIAL W              =  '', F10.4 )' ) WS
        WRITE(*,'('' INITIAL P              =  '', F10.4 )' ) PS

        WRITE(*,'(//'' START OF MARKOV CHAIN               ''//)')
        WRITE(*,'(''  NMOVE     RATIO       V/N            P''/)')

C    *******************************************************************
C    ** LOOPS OVER ALL CYCLES AND ALL MOLECULES                       **
C    *******************************************************************

        DO 100 STEP = 1, NSTEP

           DO 99 I = 1, N

              RXIOLD = RX(I)
              RYIOLD = RY(I)
              RZIOLD = RZ(I)

C          ** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **

              CALL ENERGY ( RXIOLD, RYIOLD, RZIOLD, I, RCUT, SIGMA,
     :                      VOLD, WOLD                            )

C          ** MOVE I AND PICKUP THE CENTRAL IMAGE **

              RXINEW = RXIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RYINEW = RYIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RZINEW = RZIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX

              RXINEW = RXINEW - ANINT ( RXINEW )
              RYINEW = RYINEW - ANINT ( RYINEW )
              RZINEW = RZINEW - ANINT ( RZINEW )

C          ** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **

              CALL ENERGY ( RXINEW, RYINEW, RZINEW, I, RCUT, SIGMA,
     :                      VNEW, WNEW                            )

C          ** CHECK FOR ACCEPTANCE **

              DELTV  = VNEW - VOLD
              DELTW  = WNEW - WOLD
              DELTVB = BETA * DELTV

              IF ( DELTVB .LT. 75.0 ) THEN

                 IF ( DELTV .LE. 0.0 ) THEN

                    V      = V + DELTV
                    W      = W + DELTW
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ELSEIF ( EXP ( - DELTVB ) .GT. RANF ( DUMMY ) ) THEN

                    V      = V + DELTV
                    W      = W + DELTW
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ENDIF

              ENDIF

              ACM = ACM + 1.0

C          ** CALCULATE INSTANTANEOUS VALUES **

              VN     = ( V + VLRC ) / REAL ( N )
              PRES   = DENS * TEMP + W + WLRC

C          ** CONVERT PRESSURE TO LJ UNITS **

              PRES   = PRES * SIGMA ** 3

C          ** ACCUMULATE AVERAGES **

              ACV    = ACV   + VN
              ACP    = ACP   + PRES
              ACVSQ  = ACVSQ + VN * VN
              ACPSQ  = ACPSQ + PRES * PRES

C          *************************************************************
C          ** ENDS LOOP OVER ATOMS                                    **
C          *************************************************************

99         CONTINUE

C       ** PERFORM PERIODIC OPERATIONS  **

           IF ( MOD ( STEP, IRATIO ) .EQ. 0 ) THEN

C          ** ADJUST MAXIMUM DISPLACEMENT **

              RATIO = ACATMA / REAL ( N * IRATIO )

              IF ( RATIO .GT. 0.5 ) THEN

                 DRMAX  = DRMAX  * 1.05

              ELSE

                 DRMAX  = DRMAX  * 0.95

              ENDIF

              ACATMA = 0.0

           ENDIF

           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN

C          ** WRITE OUT RUNTIME INFORMATION **

              WRITE(*,'(I8,3F12.6)') INT(ACM), RATIO, VN, PRES

           ENDIF

           IF ( MOD ( STEP, ISAVE ) .EQ. 0 ) THEN

C          ** WRITE OUT THE CONFIGURATION AT INTERVALS **

              CALL WRITCN ( CNFILE )

           ENDIF

100     CONTINUE

C    *******************************************************************
C    ** ENDS THE LOOP OVER CYCLES                                     **
C    *******************************************************************

        WRITE(*,'(//'' END OF MARKOV CHAIN          ''//)')

C    ** CHECKS FINAL VALUE OF THE POTENTIAL ENERGY IS CONSISTENT **

        CALL SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, VEND, WEND )

        IF ( ABS ( VEND - V ) .GT. 1.0E-03 ) THEN

           WRITE(*,'('' PROBLEM WITH ENERGY,'')')
           WRITE(*,'('' VEND              = '', E20.6)') VEND
           WRITE(*,'('' V                 = '', E20.6)') V

        ENDIF

C    ** WRITE OUT THE FINAL CONFIGURATION FROM THE RUN **

        CALL WRITCN ( CNFILE )

C    ** CALCULATE AND WRITE OUT RUNNING AVERAGES **

        AVV   = ACV / ACM
        ACVSQ = ( ACVSQ / ACM ) - AVV ** 2
        AVP   = ACP / ACM
        ACPSQ = ( ACPSQ / ACM ) - AVP ** 2

C    ** CALCULATE FLUCTUATIONS **

        IF ( ACVSQ .GT. 0.0 ) FLV = SQRT ( ACVSQ )
        IF ( ACPSQ .GT. 0.0 ) FLP = SQRT ( ACPSQ )

        WRITE(*,'(/'' AVERAGES ''/ )')
        WRITE(*,'('' <V/N>   = '',F10.6)') AVV
        WRITE(*,'('' <P>     = '',F10.6)') AVP

        WRITE(*,'(/'' FLUCTUATIONS ''/)')

        WRITE(*,'('' FLUCTUATION IN <V/N> = '',F10.6)') FLV
        WRITE(*,'('' FLUCTUATION IN <P>   = '',F10.6)') FLP
        WRITE(*,'(/'' END OF SIMULATION '')')

        STOP
        END



        SUBROUTINE SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, V, W )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS                 **
C    ** REAL    RX(N(,RY(N),RZ(N) THE POSITIONS OF THE ATOMS          **
C    ** REAL    V                 THE POTENTIAL ENERGY                **
C    ** REAL    W                 THE VIRIAL                          **
C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SUBROUTINE RETURNS THE TOTAL POTENTIAL ENERGY AT THE      **
C    ** BEGINNING AND END OF THE RUN.                                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        SIGMA, RMIN, RCUT, V, W
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, SIGSQ, RXIJ, RYIJ, RZIJ
        REAL        RXI, RYI, RZI, VIJ, WIJ, SR2, SR6, RIJSQ
        INTEGER     I, J

C    *******************************************************************

        OVRLAP = .FALSE.
        RCUTSQ = RCUT * RCUT
        RMINSQ = RMIN * RMIN
        SIGSQ  = SIGMA * SIGMA

        V      = 0.0
        W      = 0.0

C    ** LOOP OVER ALL THE PAIRS IN THE LIQUID **

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)

           DO 99 J = I + 1, N

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)

C          ** MINIMUM IMAGE THE PAIR SEPARATIONS **

              RXIJ  = RXIJ - ANINT ( RXIJ )
              RYIJ  = RYIJ - ANINT ( RYIJ )
              RZIJ  = RZIJ - ANINT ( RZIJ )
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RMINSQ ) THEN

                 OVRLAP = .TRUE.
                 RETURN

              ELSEIF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2 = SIGSQ / RIJSQ
                 SR6 = SR2 * SR2 * SR2
                 VIJ = SR6 * ( SR6 - 1.0 )
                 WIJ = SR6 * ( SR6 - 0.5 )
                 V   = V + VIJ
                 W   = W + WIJ

              ENDIF

99         CONTINUE

100     CONTINUE

        V = 4.0 * V
        W = 48.0 * W / 3.0

        RETURN
        END



        SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, SIGMA, V, W )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY OF ATOM I WITH ALL OTHER ATOMS.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER I                 THE ATOM OF INTEREST                **
C    ** INTEGER N                 THE NUMBER OF ATOMS                 **
C    ** REAL    RX(N),RY(N),RZ(N) THE ATOM POSITIONS                  **
C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF ATOM I           **
C    ** REAL    V                 THE POTENTIAL ENERGY OF ATOM I      **
C    ** REAL    W                 THE VIRIAL OF ATOM I                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL MOVE OF ATOM I. IT IS CALLED BEFORE AND        **
C    ** AFTER THE RANDOM DISPLACEMENT OF I.                           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        RCUT, SIGMA, RXI, RYI, RZI, V, W
        INTEGER     I

        REAL        RCUTSQ, SIGSQ, SR2, SR6
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ
        INTEGER     J

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA

        V      = 0.0
        W      = 0.0

C    ** LOOP OVER ALL MOLECULES EXCEPT I  **

        DO 100 J = 1, N

           IF ( I .NE. J ) THEN

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)

              RXIJ  = RXIJ - ANINT ( RXIJ )
              RYIJ  = RYIJ - ANINT ( RYIJ )
              RZIJ  = RZIJ - ANINT ( RZIJ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2 = SIGSQ / RIJSQ
                 SR6 = SR2 * SR2 * SR2
                 VIJ = SR6 * ( SR6 - 1.0 )
                 WIJ = SR6 * ( SR6 - 0.5 )
                 V   = V + VIJ
                 W   = W + WIJ

              ENDIF

           ENDIF

100     CONTINUE

        V = 4.0 * V
        W = 48.0 * W / 3.0

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



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 10          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

        INTEGER     NN

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED'                        )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP 'N ERROR IN READCN'
        READ ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :         FORM = 'UNFORMATTED'                            )

        WRITE ( CNUNIT ) N
        WRITE ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



