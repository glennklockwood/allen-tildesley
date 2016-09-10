********************************************************************************
** FICHE F.33.  BROWNIAN DYNAMICS FOR  A LENNARD-JONES FLUID                  **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

        PROGRAM BROWND

        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK2 / D, XIC

C    *******************************************************************
C    ** BROWNIAN DYNAMICS WITH HYDRODYNAMIC INTERACTIONS              **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ERMAK AND MCCAMMON, J CHEM PHYS, 69, 1352, 1982.              **
C    **                                                               **
C    ** THIS PROGRAM TAKES A CONFIGURATION OF LENNARD JONES ATOMS     **
C    ** AND PERFORMS A BROWNIAN DYNAMICS SIMULATION.                  **
C    ** THE ALGORITHM, DUE TO ERMAK AND MCCAMMON INCLUDES THE         **
C    ** HYDRODYNAMIC INTERACTION THROUGH EITHER THE OSEEN OR THE      **
C    ** ROTNE-PRAGER TENSOR.  BEWARE! UNDER CERTAIN CONDITIONS THE    **
C    ** APPROXIMATE DIFFUSION TENSOR MAY NOT BE POSITIVE-DEFINITE     **
C    ** IN WHICH CASE THE PROGRAM WILL FAIL IN SUBROUTINE COVAR.      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS                   **
C    ** INTEGER N3                  NUMBER OF DEGREES OF FREEDOM      **
C    ** INTEGER MSTEP               MAXIMUM NUMBER OF STEPS           **
C    ** INTEGER ISAVE               STEPS BETWEEN DATA SAVE           **
C    ** INTEGER IPRINT              STEPS BETWEEN OUTPUT              **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    FX(N),FY(N),FZ(N)   FORCES                            **
C    ** REAL    XIC(I)              CORRELATED RANDOM NORMAL DEVIATES **
C    ** REAL    D(N3,N3)            THE DIFFUSION TENSOR              **
C    ** REAL    DENS                REDUCED DENSITY                   **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    DT                  REDUCED TIMESTEP                  **
C    ** REAL    SIGMA               REDUCED LJ DIAMETER               **
C    ** REAL    ETA                 REDUCED VISCOSITY                 **
C    ** REAL    CONSII              CONSTANT FOR THE DIFFUSION TENSOR **
C    ** REAL    CONSIJ              CONSTANT FOR THE DIFFUSION TENSOR **
C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
C    ** REAL    V                   THE CONFIGURATIONAL ENERGY        **
C    ** REAL    W                   THE VIRIAL                        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SIMULATION IS PERFORMED IN A BOX OF UNIT LENGTH CENTRED   **
C    ** AT THE ORIGIN.                                                **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM USES LENNARD-JONES UNITS FOR USER INPUT AND       **
C    ** OUTPUT BUT CONDUCTS THE SIMULATION IN A BOX OF UNIT LENGTH.   **
C    ** FOR EXAMPLE, FOR A BOXLENGTH L, THE UNITS ARE:                **
C    **                                                               **
C    **     PROPERTY       LJ  UNITS            PROGRAM UNITS         **
C    **                                                               **
C    **     TEMP           EPSILON/K            EPSILON/K             **
C    **     DENS           1/SIGMA**3           1/L**3                **
C    **     ETA            SQRT(M*EPSILON/      SQRT(M*EPSILON/       **
C    **                    SIGMA**4)            L**4)                 **
C    **     DT             SQRT(M*SIGMA**2/     SQRT(M*L**2/          **
C    **                    EPSILON)             EPSILON)              **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE FORCE ( SIGMA, RCUT, CONSII, CONSIJ, V, W )        **
C    **    CALCULATES THE DIFFUSION TENSOR AND THE SYSTEMATIC FORCE   **
C    **    ON EACH ATOM IN A PARTICULAR CONFIGURATION                 **
C    ** SUBROUTINE MOVE ( DT, TEMP )                                  **
C    **    MOVES THE ATOMS                                            **
C    ** SUBROUTINE READCN (CNFILE )                                   **
C    **    READS IN A CONFIGURATION                                   **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT A CONFIGURATION                                 **
C    ** SUBROUTINE COVAR ( DT )                                       **
C    **    CALCULATES 3N CORRELATED NORMAL RANDOM DEVIATES            **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    CALCULATES A NORMAL RANDOM VARIATE FROM A DISTRIBUTION     **
C    **    WITH ZERO MEAN AND UNIT VARIANCE                           **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
C    *******************************************************************

        INTEGER     N, N3
        PARAMETER ( N = 32, N3 = 3 * N )

        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        D(N3,N3), XIC(N3)

        REAL        DENS, TEMP, DENLJ, ETA, DT
        REAL        SIGMA, RCUT, CONSII, CONSIJ
        REAL        PI, ACV, ACP, ACVSQ, ACPSQ
        REAL        AVV, AVP, FLV, FLP
        REAL        VLRC, VLRCA, VLRCR, WLRC, WLRCA, WLRCR
        REAL        V, W, RADIUS, VN, PRES, RANF, GAUSS, DUMMY
        INTEGER     STEP, I, NSTEP, ISAVE, IPRINT
        CHARACTER   TITLE*80, CNFILE*30

        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

C    ** READ INPUT DATA **

        WRITE(*,'(1H1,'' **** PROGRAM BROWND ****                 '')')
        WRITE(*,'('' BROWNIAN DYNAMICS SIMULATION                 '')')
        WRITE(*,'('' WITH HYDRODYNAMIC INTERACTIONS               '')')
        WRITE(*,'('' ENTER THE RUN TITLE                          '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF STEPS                        '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN DATA SAVES     '')')
        READ (*,*) ISAVE
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT         '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER THE CONFIGURATION FILE NAME            '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'(/)')
        WRITE(*,'('' ENTER THE FOLLOWING IN LENNARD-JONES UNITS '',/)')
        WRITE(*,'('' ENTER THE DENSITY                            '')')
        READ (*,*) DENS
        WRITE(*,'('' ENTER THE TEMPERATURE                        '')')
        READ (*,*) TEMP
        WRITE(*,'('' ENTER THE VISCOSITY                          '')')
        READ (*,*) ETA
        WRITE(*,'('' ENTER THE POTENTIAL CUTOFF DISTANCE          '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER THE TIMESTEP                           '')')
        READ (*,*) DT

C    ** WRITE INPUT DATA **

        WRITE(*,'(       //1X                    ,A     )') TITLE
        WRITE(*,'('' NUMBER OF ATOMS           '',I10   )') N
        WRITE(*,'('' NUMBER OF STEPS           '',I10   )') NSTEP
        WRITE(*,'('' SAVE FREQUENCY            '',I10   )') ISAVE
        WRITE(*,'('' OUTPUT FREQUENCY          '',I10   )') IPRINT
        WRITE(*,'('' CONFIGURATION FILE  NAME  '',A     )') CNFILE
        WRITE(*,'('' TEMPERATURE               '',F10.4 )') TEMP
        WRITE(*,'('' DENSITY                   '',F10.4 )') DENS
        WRITE(*,'('' VISCOSITY                 '',F10.4 )') ETA
        WRITE(*,'('' POTENTIAL CUTOFF          '',F10.4 )') RCUT
        WRITE(*,'('' TIMESTEP                  '',F10.4 )') DT

C    ** READ IN INITIAL CONFIGURATION **

        CALL READCN ( CNFILE )

C    ** CONVERT INPUT DATA TO PROGRAM UNITS **

        SIGMA  = ( DENS / REAL(N) ) ** ( 1.0 / 3.0 )
        RCUT   = RCUT * SIGMA
        DENLJ  = DENS
        DENS   = DENS / ( SIGMA ** 3 )
        ETA    = ETA  / ( SIGMA ** 2 )
        RADIUS = SIGMA * 0.5
        CONSII = TEMP / 6.0 / PI / ETA / RADIUS
        CONSIJ = TEMP / 8.0 / PI / ETA
        DT     = DT * SIGMA

        IF ( RCUT .GT. 0.5 ) STOP 'CUT-OFF TOO LARGE'

C    ** LONG RANGE CORRECTIONS              **
C    ** SPECIFIC TO THE LENNARD JONES FLUID **

        VLRCR =  ( 8.0 * PI * DENLJ * ( SIGMA / RCUT ) ** 9 ) / 9.0
        VLRCA = -( 8.0 * PI * DENLJ * ( SIGMA / RCUT ) ** 3 ) / 3.0
        VLRC  = VLRCR + VLRCA
        WLRCR =  4.0 * REAL(N) * VLRCR
        WLRCA =  2.0 * REAL(N) * VLRCA
        WLRC  = WLRCR + WLRCA

C    ** ZERO ACCUMULATORS **

        ACV    = 0.0
        ACP    = 0.0
        ACVSQ  = 0.0
        ACPSQ  = 0.0
        FLV    = 0.0
        FLP    = 0.0

C    ** WRITE OUT SOME USEFUL INFORMATION **

        WRITE(*,'('' SIGMA/BOX              =  '',F10.4)')  SIGMA
        WRITE(*,'('' RCUT/BOX               =  '',F10.4)')  RCUT
        WRITE(*,'('' DT                     =  '',F10.4)')  DT

        WRITE(*,'(/'' ** BROWNIAN DYNAMICS BEGINS ** ''/// )')
        WRITE(*,'(''   STEP       V/N       P    ''/ )')

C    *******************************************************************
C    ** MAIN LOOP BEGINS                                              **
C    *******************************************************************

        DO 100 STEP = 1, NSTEP

C       ** CALCULATE THE DIFFUSION TENSOR AND SYSTEMATIC **
C       ** FORCES AT THE BEGINNING OF THE STEP           **

           CALL FORCE ( SIGMA, RCUT, CONSII, CONSIJ, V, W )

C       ** CALCULATE THE CORRELATED NORMAL VARIATES **

           CALL COVAR ( DT )

C       ** MOVE THE ATOMS **

           CALL MOVE ( DT, TEMP )

C       ** CALCULATE INSTANTANEOUS VALUES FOR PREVIOUS STEP **

           VN     = V / REAL ( N ) + VLRC
           PRES   = DENS * TEMP + W + WLRC

C       ** CONVERT PRESSURE TO LJ UNITS **

           PRES   = PRES * SIGMA ** 3

C       ** UPDATE ACCUMULATORS **

           ACV    = ACV   + VN
           ACP    = ACP   + PRES
           ACVSQ  = ACVSQ + VN * VN
           ACPSQ  = ACPSQ + PRES * PRES

C       ** WRITE OUT RUNTIME INFORMATION **

           IF( MOD( STEP, IPRINT ) .EQ. 0 ) THEN

              WRITE( *, '( I8, 3F12.6 )' ) STEP, VN, PRES

           ENDIF

C       ** WRITE OUT THE CONFIGURATION AT INTERVALS **

           IF ( MOD ( STEP, ISAVE ) .EQ. 0 ) THEN

              CALL WRITCN ( CNFILE )

           ENDIF

100     CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/'' ** BROWNIAN DYNAMICS ENDS  ** ''///)')

C    ** CALCULATE AND WRITE OUT RUNNING AVERAGES **

        AVV   = ACV / REAL ( NSTEP )
        AVP   = ACP / REAL ( NSTEP )
        ACVSQ = ( ACVSQ / REAL ( NSTEP ) )  - AVV ** 2
        ACPSQ = ( ACPSQ / REAL ( NSTEP ) )  - AVP ** 2

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

C    ** WRITE OUT THE FINAL CONFIGURATION FROM THE RUN **

        CALL WRITCN ( CNFILE )

        STOP
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 10          **
C    *******************************************************************

        INTEGER      N
        PARAMETER (  N = 32 )
        CHARACTER    CNFILE * ( * )
        REAL         RX(N), RY(N), RZ(N)

        INTEGER      CNUNIT, NN
        PARAMETER (  CNUNIT = 10 )

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :         FORM = 'UNFORMATTED'                            )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP 'PROBLEM WITH N IN READCN'
        READ ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
C    *******************************************************************

        INTEGER      N
        PARAMETER (  N = 32 )
        CHARACTER    CNFILE * ( * )
        REAL         RX(N), RY(N), RZ(N)

        INTEGER      CNUNIT
        PARAMETER (  CNUNIT = 10 )

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED'                        )

        WRITE ( CNUNIT ) N
        WRITE ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE FORCE ( SIGMA, RCUT, CONSII, CONSIJ, V, W )

        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK2 / D, XIC

C    *******************************************************************
C    ** ROUTINE TO COMPUTE SYSTEMATIC FORCES AND THE DIFFUSION TENSOR **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER N3                 NUMBER OF DEGREES OF FREEDOM       **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
C    ** REAL    D(N3,N3)           THE DIFFUSION TENSOR               **
C    ** REAL    XIC(N3)            CORRELATED RANDOM NORMAL DEVIATES  **
C    ** REAL    SIGMA              THE LJ LENGTH PARAMETER            **
C    ** REAL    RCUT               THE CUT-OFF DISTANCE               **
C    ** REAL    CONSII             CONSTANT IN THE DIFFUSION TENSOR   **
C    ** REAL    CONSIJ             CONSTANT IN THE DIFFUSION TENSOR   **
C    ** REAL    V                  THE POTENTIAL ENERGY               **
C    ** REAL    W                  THE VIRIAL                         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FORCE IS CALLED IN A BROWNIAN DYNAMICS PROGRAM TO CALCULATE   **
C    ** THE SYSTEMATIC FORCE ON EACH ATOM AND THE ELEMENTS OF THE     **
C    ** DIFFUSION TENSOR. A CUTOFF IS APPLIED TO THE SYSTEMATIC FORCE **
C    ** IT IS ASSUMED THAT THE LENNARD-JONES SIGMA IS ALSO THE ATOMIC **
C    ** DIAMETER.                                                     **
C    *******************************************************************

        INTEGER     N, N3
        PARAMETER ( N = 32, N3 = N * 3 )

        REAL        SIGMA, RCUT, CONSII, CONSIJ, V, W
        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        D(N3,N3), XIC(N3)

        INTEGER     IC, JC, I, J
        REAL        RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, FIJ, RCUTSQ
        REAL        SIGSQ, FXI, FYI, FZI, SR2, SR6, RIJ, RRIJSQ, SIGSQ6
        REAL        RIJSQ ,RXIJ, RYIJ, RZIJ, VIJ, WIJ, OIJ, RPIJ

C    *******************************************************************

        SIGSQ  = SIGMA ** 2
        RCUTSQ = RCUT ** 2
        SIGSQ6 = SIGSQ / 6.0

C    ** ZERO FORCES AND POTENTIAL **

        DO 10 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

10      CONTINUE

        V = 0.0
        W = 0.0

C       ** LOOP OVER ALL PAIRS OF ATOMS **

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           FXI = FX(I)
           FYI = FY(I)
           FZI = FZ(I)
           IC = 3 * ( I - 1) + 1

           DO 99 J = I + 1, N

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)
              RXIJ  = RXIJ - ANINT( RXIJ )
              RYIJ  = RYIJ - ANINT( RYIJ )
              RZIJ  = RZIJ - ANINT( RZIJ )
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

C          ** CALCULATE OFF-DIAGONAL BLOCKS OF DIFFUSION TENSOR **
C          ** HERE WE ASSUME THE ROTNE-PRAGER TENSOR FORM       **
C          ** TAKE RPIJ = 0 INSTEAD BELOW FOR OSEEN TENSOR      **

              JC     = ( J - 1 ) * 3 + 1
              RIJ    = SQRT ( RIJSQ )
              RRIJSQ = 1.0 / RIJSQ
              OIJ    = CONSIJ / RIJ
              RPIJ   = OIJ * SIGSQ6 * RRIJSQ

              D( IC  , JC   ) = OIJ + RPIJ
     :           + ( OIJ - 3.0 * RPIJ ) * RXIJ * RXIJ * RRIJSQ
              D( IC+1, JC+1 ) = OIJ + RPIJ
     :           + ( OIJ - 3.0 * RPIJ ) * RYIJ * RYIJ * RRIJSQ
              D( IC+2, JC+2 ) = OIJ + RPIJ
     :           + ( OIJ - 3.0 * RPIJ ) * RZIJ * RZIJ * RRIJSQ
              D( IC  , JC+1 ) =
     :             ( OIJ - 3.0 * RPIJ ) * RXIJ * RYIJ * RRIJSQ
              D( IC  , JC+2 ) =
     :             ( OIJ - 3.0 * RPIJ ) * RXIJ * RZIJ * RRIJSQ
              D( IC+1, JC+2 ) =
     :             ( OIJ - 3.0 * RPIJ ) * RYIJ * RZIJ * RRIJSQ
              D( IC+1, JC   ) = D( IC  , JC+1 )
              D( IC+2, JC   ) = D( IC  , JC+2 )
              D( IC+2, JC+1 ) = D( IC+1, JC+2 )

C          ** CALCULATE SYSTEMATIC FORCES **

              IF( RIJSQ .LT. RCUTSQ ) THEN

                 SR2   = SIGSQ * RRIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )
                 WIJ   = SR6 * ( SR6 - 0.5 )
                 FIJ   = WIJ * RRIJSQ
                 FXIJ  = FIJ * RXIJ
                 FYIJ  = FIJ * RYIJ
                 FZIJ  = FIJ * RZIJ
                 V     = V  + VIJ
                 W     = W  + WIJ
                 FXI   = FXI + FXIJ
                 FYI   = FYI + FYIJ
                 FZI   = FZI + FZIJ
                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ

              ENDIF

99         CONTINUE

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI

100     CONTINUE

C    ** INCORPORATE FACTORS **

        V = V * 4.0
        W = W * 48.0 / 3.0

        DO 50 I = 1, N

           FX(I) = FX(I) * 48.0
           FY(I) = FY(I) * 48.0
           FZ(I) = FZ(I) * 48.0

C       ** CALCULATE ON-DIAGONAL BLOCKS OF DIFFUSION TENSOR **

           IC = 3 * ( I - 1 ) + 1

           D( IC  , IC   ) = CONSII
           D( IC+1, IC+1 ) = CONSII
           D( IC+2, IC+2 ) = CONSII
           D( IC  , IC+1 ) = 0.0
           D( IC  , IC+2 ) = 0.0
           D( IC+1, IC+2 ) = 0.0

50      CONTINUE

C    ** FILL THE LOWER TRIANGLE OF THE DIFFUSION TENSOR **

        DO 70 IC = 1, N3 - 1

           DO 60 JC = IC + 1, N3

              D( JC, IC ) = D( IC, JC )

60         CONTINUE

70      CONTINUE

        RETURN
        END



        SUBROUTINE COVAR  ( DT )

        COMMON / BLOCK2 / D, XIC

C    *******************************************************************
C    ** ROUTINE TO COMPUTE 3N CORRELATED RANDOM NORMAL DEVIATES.      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                NUMBER OF ATOMS                      **
C    ** INTEGER N3               NUMBER OF DEGREES OF FREEDOM         **
C    ** REAL    D(N3,N3)         THE DIFFUSION TENSOR                 **
C    ** REAL    XIC(N3)          CORRELATED RANDOM NORMAL DEVIATES    **
C    ** REAL    XI(N3)           UNCORRELATED RANDOM NORMAL DEVIATES  **
C    ** REAL    L(N3,N3)         A LOWER TRIANGULAR MATRIX            **
C    ** REAL    DT               REDUCED TIMESTEP                     **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** COVAR IS CALLED IN A BROWNIAN DYNAMICS SIMULATION AFTER THE   **
C    ** THE DIFFUSION TENSOR HAS BEEN CONSTRUCTED IN FORCE. ON EXIT   **
C    ** THE ARRAY XIC CONTAINS THE CORRELATED GAUSSIAN DISPLACEMENTS. **
C    **                                                               **
C    **    *****************************************************      **
C    **    **               WARNING                           **      **
C    **    **                                                 **      **
C    **    ** THIS ROUTINE PERFORMS A STANDARD DECOMPOSITION  **      **
C    **    ** OF A POSITIVE DEFINITE MATRIX D INTO A PRODUCT  **      **
C    **    ** L * L(TRANSPOSE), WHERE L IS A LOWER TRIANGULAR **      **
C    **    ** MATRIX. THIS IS EXPENSIVE FOR A LARGE MATRIX    **      **
C    **    ** AND YOU MAY FIND A MORE EFFICIENT OR ACCURATE   **      **
C    **    ** MACHINE CODE ROUTINE IN THE COMMON SCIENTIFIC   **      **
C    **    ** LIBRARIES SUCH AS NAG OR IMSL.  IF THE MATRIX   **      **
C    **    ** IS NOT POSITIVE DEFINITE THE METHOD WILL FAIL.  **      **
C    **    *****************************************************      **
C    **                                                               **
C    *******************************************************************

        INTEGER     N, N3
        PARAMETER ( N = 32, N3 = N * 3 )

        REAL        D(N3,N3), XIC(N3)
        REAL        DT

        INTEGER     I, J, K, IC
        REAL        GAUSS, DUMMY, L(N3,N3), SUM, XI(N3)

C    *******************************************************************

C    ** CALCULATE THE LOWER TRIANGULAR MATRIX L **

        L(1, 1) = SQRT ( D(1, 1) )
        L(2, 1) = D(2, 1) / L(1, 1)
        L(2, 2) = SQRT ( D(2, 2) - L(2, 1) * L(2, 1) )

        DO 60 I = 3, N3

           L(I, 1) = D(I, 1) / L(1, 1)

           DO 40 J = 2, I - 1

              SUM = 0.0

              DO 30 K = 1, J - 1

                 SUM = SUM + L(I, K) * L(J, K)

30            CONTINUE

              L(I, J) = ( D(I, J) - SUM ) / L(J, J)

40         CONTINUE

           SUM = 0.0

           DO 50 K = 1, I - 1

              SUM = SUM + L(I, K) * L(I, K)

50         CONTINUE

           L(I, I) = SQRT ( D(I, I) - SUM )

60      CONTINUE

C    ** CALCULATE CORRELATED RANDOM DISPLACEMENTS **

        DO 80 I = 1, N3

C       ** CALCULATE UNCORRELATED RANDOM NORMAL DEVIATES **
C       ** WITH ZERO MEAN AND VARIANCE 2.0 * DT          **

           XI(I) = GAUSS ( DUMMY ) * SQRT ( 2.0 * DT )
           SUM = 0.0

           DO 70 J = 1, I

              SUM = SUM + L(I, J) * XI(J)

70         CONTINUE

           XIC(I) = SUM

80      CONTINUE

        RETURN
        END



        SUBROUTINE MOVE ( DT, TEMP )

        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK2 / D, XIC

C    *******************************************************************
C    ** ROUTINE TO MOVE THE ATOMS IN A BROWNIAN DYNAMICS SIMULATION   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER N3                 NUMBER OF DEGREES OF FREEDOM       **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
C    ** REAL    D(N3,N3)           THE DIFFUSION TENSOR               **
C    ** REAL    XIC(N3)            CORRELATED RANDOM NORMAL DEVIATES  **
C    ** REAL    DT                 REDUCED TIMESTEP                   **
C    ** REAL    TEMP               REDUCED TEMPERATURE                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** MOVE IS CALLED AFTER FORCE AND COVAR TO MOVE THE ATOMS.       **
C    *******************************************************************

        INTEGER     N, N3
        PARAMETER ( N = 32, N3 = N * 3 )

        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        D(N3,N3), XIC(N3)
        REAL        DT, TEMP

        REAL        F(N3), SUMX, SUMY, SUMZ
        INTEGER     I, J, IC, JC

C    *******************************************************************

C    ** PLACE FORCES IN A TEMPORARY ARRAY OF SIZE 3N **

        DO 10 I = 1, N

           IC      = ( I - 1 ) * 3 + 1
           F(IC)   = FX(I)
           F(IC+1) = FY(I)
           F(IC+2) = FZ(I)

10      CONTINUE

C    ** MOVE THE ATOMS **

        DO 30 I = 1, N

           IC   = ( I - 1 ) * 3  + 1
           SUMX = 0.0
           SUMY = 0.0
           SUMZ = 0.0

           DO 20 JC = 1, N3

              SUMX = SUMX + D( IC  , JC  ) * F(JC)
              SUMY = SUMY + D( IC+1, JC  ) * F(JC)
              SUMZ = SUMZ + D( IC+2, JC  ) * F(JC)

20         CONTINUE

           RX(I) = RX(I) + SUMX * DT / TEMP + XIC( IC )
           RY(I) = RY(I) + SUMY * DT / TEMP + XIC( IC + 1 )
           RZ(I) = RZ(I) + SUMZ * DT / TEMP + XIC( IC + 2 )

30      CONTINUE

        RETURN
        END



        REAL FUNCTION RANF ( DUMMY )

C    *******************************************************************
C    ** FUNCTION RANF RETURNS A UNIFORM RANDOM VARIATE BETWEEN 0 AND 1**
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

        END



        REAL FUNCTION GAUSS ( DUMMY )

C    *******************************************************************
C    ** FUNCTION GAUSS RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM   **
C    ** A DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.              **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
C    **          ADDISON-WESLEY), 1978.                               **
C    *******************************************************************

        REAL        A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138, A3 = 0.252408784 )
        PARAMETER ( A5 = 0.076542912, A7 = 0.008355968 )
        PARAMETER ( A9 = 0.029899776                   )

        REAL        SUM, R, R2
        INTEGER     I

C    *******************************************************************

        SUM = 0.0

        DO 10 I = 1, 12

           SUM = SUM + RANF ( DUMMY )

10      CONTINUE

        R  = ( SUM - 6.0 ) / 4.0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 + A1 )
     :          * R

        END



