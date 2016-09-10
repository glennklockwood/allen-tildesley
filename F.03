********************************************************************************
** FICHE F.3.  LOW-STORAGE MD PROGRAMS USING THE LEAPFROG VERLET ALGORITHM    **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO SEPARATE PARTS: FORTRAN AND BASIC VERSIONS.               **
C    *******************************************************************



C    *******************************************************************
C    ** FICHE F.3 - PART A                                            **
C    ** FORTRAN PROGRAM USING THE LEAPFROG ALGORITHM.                 **
C    *******************************************************************



        PROGRAM FROGGY

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** FORTRAN PROGRAM TO CONDUCT MOLECULAR DYNAMICS OF ATOMS.       **
C    **                                                               **
C    ** A SPECIAL LOW-STORAGE FORM OF THE LEAPFROG VERLET ALGORITHM   **
C    ** IS USED AND WE TAKE THE LENNARD-JONES POTENTIAL AS AN EXAMPLE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FINCHAM AND HEYES, CCP5 QUARTERLY, 6, 4, 1982.                **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN CONFIGURATION                                     **
C    ** SUBROUTINE FORCE ( DT, SIGMA, RCUT, NEWV, NEWVC, NEWW )       **
C    **    CALCULATES THE ACCELERATIONS AND ADVANCES THE VELOCITIES   **
C    **    FROM T - DT/2 TO T + DT/2. ALSO CALCULATES POTENTIAL       **
C    **    ENERGY AND VIRIAL AT TIMESTEP T.                           **
C    ** SUBROUTINE MOVE ( DT )                                        **
C    **    ADVANCES POSITIONS FROM T TO T + DT.                       **
C    ** SUBROUTINE KINET ( NEWK )                                     **
C    **    CALCULATES KINETIC ENERGY.                                 **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT CONFIGURATION                                   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF ATOMS                     **
C    ** REAL    DT                TIMESTEP                            **
C    ** REAL    RX(N),RY(N),RZ(N) ATOMIC POSITIONS                    **
C    ** REAL    VX(N),VY(N),VZ(N) ATOMIC VELOCITIES                   **
C    ** REAL    ACV,ACK ETC.      AVERAGE VALUE ACCUMULATORS          **
C    ** REAL    AVV,AVK ETC.      AVERAGE VALUES                      **
C    ** REAL    ACVSQ, ACKSQ ETC. AVERAGE SQUARED VALUE ACCUMULATORS  **
C    ** REAL    FLV,FLK ETC.      FLUCTUATION AVERAGES                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE LEAPFROG ALGORITHM IS OF THE FORM                         **
C    ** VX(T + DT/2) = VX(T - DT/2) + DT * AX(T)  (SIMILARLY FOR Y,Z) **
C    ** RX(T + DT)   = RX(T) + DT * VX(T + DT/2)  (SIMILARLY FOR Y,Z) **
C    ** TO SAVE STORAGE IN THIS PROGRAM WE ACCUMULATE VALUES AX(T)    **
C    ** DIRECTLY ONTO THE VELOCITIES VX(T - DT/2) IN THE FORCE LOOP.  **
C    ** THIS MEANS THAT AN APPROXIMATION MUST BE USED FOR THE KINETIC **
C    ** ENERGY AT EACH TIME STEP:                                     **
C    ** K = ( OLDK + NEWK ) / 2 + ( OLDV - 2 * V + NEWV ) / 8         **
C    ** WHERE K    = K( T        ), V    = V( T )                     **
C    **       OLDK = K( T - DT/2 ), OLDV = V( T - DT )                **
C    **       NEWK = K( T + DT/2 ), NEWV = V( T + DT )                **
C    ** AT THE START OF A STEP THE FOLLOWING VARIABLES ARE STORED:    **
C    ** R    : R(STEP)         V    : V(STEP+1/2)                     **
C    ** OLDV : V(STEP-2)       V    : V(STEP-1)        NEWV : V(STEP) **
C    ** OLDK : K(STEP-1/2)     NEWK : K(STEP+1/2)                     **
C    ** THIS PROGRAM USES UNIT 10 FOR CONFIGURATION INPUT AND OUTPUT  **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM USES LENNARD-JONES REDUCED UNITS FOR USER INPUT   **
C    ** AND OUTPUT BUT CONDUCTS SIMULATION IN A BOX OF UNIT LENGTH.   **
C    ** SUMMARY FOR BOX LENGTH L, ATOMIC MASS M, AND LENNARD-JONES    **
C    ** POTENTIAL PARAMETERS SIGMA AND EPSILON:                       **
C    **                OUR PROGRAM           LENNARD-JONES SYSTEM     **
C    ** LENGTH         L                     SIGMA                    **
C    ** MASS           M                     M                        **
C    ** ENERGY         EPSILON               EPSILON                  **
C    ** TIME           SQRT(M*L**2/EPSILON)  SQRT(M*SIGMA**2/EPSILON) **
C    ** VELOCITY       SQRT(EPSILON/M)       SQRT(EPSILON/M)          **
C    ** PRESSURE       EPSILON/L**3          EPSILON/SIGMA**3         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        INTEGER     STEP, NSTEP, IPRINT
        REAL        DENS, NORM, RCUT, DT, SIGMA
        REAL        V, K, E, W
        REAL        OLDK, NEWK, OLDV, NEWV, NEWW
        REAL        OLDVC, NEWVC, VC, KC, EC
        REAL        VN, KN, EN, ECN, PRES, TEMP
        REAL        ACV, ACK, ACE, ACEC, ACP, ACT
        REAL        AVV, AVK, AVE, AVEC, AVP, AVT
        REAL        ACVSQ, ACKSQ, ACESQ, ACECSQ, ACPSQ, ACTSQ
        REAL        FLV, FLK, FLE, FLEC, FLP, FLT
        REAL        SR3, SR9, VLRC, WLRC, PI, SIGCUB
        CHARACTER   CNFILE*30, TITLE*80
        REAL        FREE

        PARAMETER ( FREE = 3.0 )
        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

C    ** READ IN INITIAL PARAMETERS **

        WRITE(*,'(1H1,'' **** PROGRAM FROGGY ****                  '')')
        WRITE(*,'(//, '' MOLECULAR DYNAMICS OF LENNARD-JONES ATOMS '')')
        WRITE(*,'(    '' LEAPFROG ALGORITHM WITH MINIMAL STORAGE   '')')

        WRITE(*,'('' ENTER RUN TITLE                               '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF STEPS                         '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT LINES    '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE FOLLOWING IN LENNARD-JONES UNITS    '')')
        WRITE(*,'('' ENTER DENSITY                                 '')')
        READ (*,*) DENS
        WRITE(*,'('' ENTER POTENTIAL CUTOFF DISTANCE               '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER TIME STEP                               '')')
        READ (*,*) DT

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' NUMBER OF ATOMS  = '',I10  )') N
        WRITE(*,'('' NUMBER OF STEPS  = '',I10  )') NSTEP
        WRITE(*,'('' OUTPUT FREQUENCY = '',I10  )') IPRINT
        WRITE(*,'('' POTENTIAL CUTOFF = '',F10.4)') RCUT
        WRITE(*,'('' DENSITY          = '',F10.4)') DENS
        WRITE(*,'('' TIME STEP        = '',F10.6)') DT

C    ** READ CONFIGURATION INTO COMMON / BLOCK1 / VARIABLES **

        CALL READCN ( CNFILE )

C    ** CONVERT INPUT DATA TO PROGRAM UNITS **

        SIGMA = ( DENS / REAL ( N ) ) ** ( 1.0 / 3.0 )
        RCUT  = RCUT * SIGMA
        DT    = DT * SIGMA
        DENS  = DENS / ( SIGMA ** 3 )

C    ** CALCULATE LONG-RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        SIGCUB = SIGMA ** 3
        VLRC = ( 8.0 /9.0 ) * PI * DENS * SIGCUB * REAL ( N )
     :           * ( SR9 - 3.0 * SR3 )
        WLRC = ( 16.0 / 9.0 ) * PI * DENS * SIGCUB * REAL ( N )
     :           * ( 2.0 * SR9 - 3.0 * SR3 )

C    ** ZERO ACCUMULATORS **

        ACV  = 0.0
        ACK  = 0.0
        ACE  = 0.0
        ACEC = 0.0
        ACP  = 0.0
        ACT  = 0.0

        ACVSQ  = 0.0
        ACKSQ  = 0.0
        ACESQ  = 0.0
        ACECSQ = 0.0
        ACPSQ  = 0.0
        ACTSQ  = 0.0

        FLV  = 0.0
        FLK  = 0.0
        FLE  = 0.0
        FLEC = 0.0
        FLP  = 0.0
        FLT  = 0.0

C    ** CALCULATE INITIAL VALUES **

        CALL FORCE ( -DT, SIGMA, RCUT, NEWV, NEWVC, NEWW )
        CALL MOVE  ( -DT )
        CALL FORCE ( -DT, SIGMA, RCUT, V, VC, W )
        CALL FORCE (  DT, SIGMA, RCUT, V, VC, W )
        CALL KINET ( OLDK )
        CALL MOVE  (  DT )
        CALL FORCE (  DT, SIGMA, RCUT, NEWV, NEWVC, NEWW )
        CALL KINET ( NEWK )

C    ** INCLUDE LONG-RANGE CORRECTIONS **

        V    = V + VLRC
        W    = W + WLRC
        NEWV = NEWV + VLRC
        NEWW = NEWW + WLRC

        IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1

        WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
        WRITE(*,10001)

C    *******************************************************************
C    ** MAIN LOOP BEGINS                                              **
C    *******************************************************************

        DO 1000 STEP = 1, NSTEP

C       ** IMPLEMENT ALGORITHM **

           CALL MOVE ( DT )

           OLDV  = V
           V     = NEWV
           OLDVC = VC
           VC    = NEWVC
           W     = NEWW

           CALL FORCE ( DT, SIGMA, RCUT, NEWV, NEWVC, NEWW )

C       ** INCLUDE LONG-RANGE CORRECTIONS **

           NEWV = NEWV + VLRC
           NEWW = NEWW + WLRC

C       ** CALCULATE KINETIC ENERGY AT CURRENT STEP **

           K  = ( NEWK + OLDK ) / 2.0
     :        + ( NEWV  - 2.0 * V  + OLDV  ) / 8.0
           KC = ( NEWK + OLDK ) / 2.0
     :        + ( NEWVC - 2.0 * VC + OLDVC ) / 8.0
           OLDK = NEWK

           CALL KINET ( NEWK )

C       ** CALCULATE INSTANTANEOUS VALUES **

           E    = K + V
           EC   = KC + VC
           VN   = V  / REAL ( N )
           KN   = K  / REAL ( N )
           EN   = E  / REAL ( N )
           ECN  = EC / REAL ( N )
           TEMP = 2.0 * KN / FREE
           PRES = DENS * TEMP + W

C       ** CONVERT TO LENNARD-JONES UNITS **

           PRES = PRES * SIGMA ** 3

C       ** INCREMENT ACCUMULATORS **

           ACE  = ACE  + EN
           ACEC = ACEC + ECN
           ACK  = ACK  + KN
           ACV  = ACV  + VN
           ACP  = ACP  + PRES

           ACESQ  = ACESQ  + EN  ** 2
           ACECSQ = ACECSQ + ECN ** 2
           ACKSQ  = ACKSQ  + KN  ** 2
           ACVSQ  = ACVSQ  + VN  ** 2
           ACPSQ  = ACPSQ  + PRES ** 2

C       ** OPTIONALLY PRINT INFORMATION **

           IF ( MOD( STEP, IPRINT ) .EQ. 0 ) THEN

              WRITE(*,'(1X,I8,6(2X,F10.4))')
     :                 STEP, EN, ECN, KN, VN, PRES, TEMP
           ENDIF

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

C    ** OUTPUT RUN AVERAGES **

        NORM   = REAL ( NSTEP )

        AVE  = ACE  / NORM
        AVEC = ACEC / NORM
        AVK  = ACK  / NORM
        AVV  = ACV  / NORM
        AVP  = ACP  / NORM

        ACESQ  = ( ACESQ  / NORM ) - AVE  ** 2
        ACECSQ = ( ACECSQ / NORM ) - AVEC ** 2
        ACKSQ  = ( ACKSQ  / NORM ) - AVK  ** 2
        ACVSQ  = ( ACVSQ  / NORM ) - AVV  ** 2
        ACPSQ  = ( ACPSQ  / NORM ) - AVP  ** 2

        IF ( ACESQ  .GT. 0.0 ) FLE  = SQRT ( ACESQ  )
        IF ( ACECSQ .GT. 0.0 ) FLEC = SQRT ( ACECSQ )
        IF ( ACKSQ  .GT. 0.0 ) FLK  = SQRT ( ACKSQ  )
        IF ( ACVSQ  .GT. 0.0 ) FLV  = SQRT ( ACVSQ  )
        IF ( ACPSQ  .GT. 0.0 ) FLP  = SQRT ( ACPSQ  )

        AVT = AVK * 2.0 / FREE
        FLT = FLK * 2.0 / FREE

        WRITE(*,'('' AVERAGES'',6(2X,F10.5))')
     :            AVE, AVEC, AVK, AVV, AVP, AVT
        WRITE(*,'('' FLUCTS  '',6(2X,F10.5))')
     :            FLE, FLEC, FLK, FLV, FLP, FLT

C    ** WRITE OUT FINAL CONFIGURATION **

        CALL WRITCN ( CNFILE )

        STOP

10001   FORMAT(//1X,'TIMESTEP  ..ENERGY..  CUTENERGY.'
     :              '  ..KINETIC.  ..POTENT..',
     :              '  .PRESSURE.  ..TEMPER..'/)
        END



        SUBROUTINE FORCE ( DT, SIGMA, RCUT, V, VC, W )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** ROUTINE TO COMPUTE FORCES WITH CUTOFF & MINIMUM IMAGING.      **
C    **                                                               **
C    ** THE ROUTINE ACTUALLY RETURNS TWO DIFFERENT POTENTIAL ENERGIES.**
C    ** V IS CALCULATED USING THE UNSHIFTED LENNARD-JONES POTENTIAL.  **
C    ** WHEN LONG-RANGE TAIL CORRECTIONS ARE ADDED, THIS MAY BE USED  **
C    ** TO CALCULATE THERMODYNAMIC INTERNAL ENERGY ETC.               **
C    ** VC IS CALCULATED USING THE SHIFTED LENNARD-JONES POTENTIAL    **
C    ** WITH NO DISCONTINUITY AT THE CUTOFF.  THIS MAY BE USED TO     **
C    ** CHECK ENERGY CONSERVATION.                                    **
C    ** PROGRAM UNITS MAKE EPSILON = MASS = 1, BUT SIGMA IS NOT UNITY **
C    ** SINCE WE TAKE A BOX OF UNIT LENGTH.                           **
C    ** THE ROUTINE IS COMBINED WITH THE LEAPFROG ALGORITHM FOR       **
C    ** ADVANCING VELOCITIES, I.E. FORCES ARE ACCUMULATED DIRECTLY    **
C    ** ONTO VELOCITIES.                                              **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        V, RCUT, DT, SIGMA, W, VC

        INTEGER     I, J, NCUT
        REAL        RCUTSQ, VCUT, SIGSQ
        REAL        RXI, RYI, RZI, VXI, VYI, VZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ
        REAL        SR2, SR6, SR12
        REAL        VIJ, WIJ, VELIJ, DVX, DVY, DVZ

C    *******************************************************************

C    ** CALCULATE SQUARED DISTANCES FOR INNER LOOP **

        SIGSQ  = SIGMA ** 2
        RCUTSQ = RCUT ** 2

C    ** CONDUCT OUTER LOOP OVER ATOMS **

        V    = 0.0
        W    = 0.0
        NCUT = 0

        DO 1000 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           VXI = VX(I)
           VYI = VY(I)
           VZI = VZ(I)

C       ** CONDUCT INNER LOOP OVER ATOMS **

           DO 900 J = I + 1, N

              RXIJ = RXI - RX(J)
              RXIJ = RXIJ - ANINT ( RXIJ )

              IF ( ABS ( RXIJ ) .LT. RCUT ) THEN

                 RYIJ = RYI - RY(J)
                 RYIJ = RYIJ - ANINT ( RYIJ )
                 RIJSQ = RXIJ ** 2 + RYIJ ** 2

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    RZIJ = RZI - RZ(J)
                    RZIJ = RZIJ - ANINT ( RZIJ )
                    RIJSQ = RIJSQ + RZIJ ** 2

                    IF ( RIJSQ .LT. RCUTSQ ) THEN

                       SR2   = SIGSQ / RIJSQ
                       SR6   = SR2 ** 3
                       SR12  = SR6 ** 2
                       VIJ   = 4.0 * ( SR12 - SR6 )
                       WIJ   = 24.0 * ( 2.0 * SR12 - SR6 )
                       VELIJ = WIJ * DT / RIJSQ
                       DVX   = VELIJ * RXIJ
                       DVY   = VELIJ * RYIJ
                       DVZ   = VELIJ * RZIJ
                       VXI   = VXI + DVX
                       VYI   = VYI + DVY
                       VZI   = VZI + DVZ
                       VX(J) = VX(J) - DVX
                       VY(J) = VY(J) - DVY
                       VZ(J) = VZ(J) - DVZ
                       V     = V + VIJ
                       W     = W + WIJ
                       NCUT  = NCUT + 1

                    ENDIF

                 ENDIF

              ENDIF

900        CONTINUE

C       ** END OF INNER LOOP **

           VX(I) = VXI
           VY(I) = VYI
           VZ(I) = VZI

1000    CONTINUE

C    ** END OF OUTER LOOP **

C    ** CALCULATE POTENTIAL SHIFT **

        SR2  = SIGSQ / RCUTSQ
        SR6  = SR2 * SR2 * SR2
        SR12 = SR6 * SR6
        VIJ  = 4.0 * ( SR12 - SR6 )
        VC   = V - REAL ( NCUT ) * VIJ

        W    = W / 3.0

        RETURN
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)

        INTEGER     CNUNIT, NN
        PARAMETER ( CNUNIT = 10 )

C     ******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP ' INCORRECT NUMBER OF ATOMS '
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION TO UNIT 10           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C    ****************************************************************

        OPEN  ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :          FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE MOVE ( DT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** LEAPFROG VERLET ALGORITHM FOR ADVANCING POSITIONS             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        DT

        INTEGER     I

C    *******************************************************************

        DO 1000 I = 1, N

           RX(I) = RX(I) + VX(I) * DT
           RY(I) = RY(I) + VY(I) * DT
           RZ(I) = RZ(I) + VZ(I) * DT

1000    CONTINUE

        RETURN
        END



        SUBROUTINE KINET ( K )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** COMPUTES KINETIC ENERGY                                       **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        K

        INTEGER     I

C    *******************************************************************

        K  = 0.0

        DO 100 I = 1, N

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

100     CONTINUE

        K = K * 0.5

        RETURN
        END



C    *******************************************************************
C    ** FICHE F.3 - PART B                                            **
C    ** BASIC PROGRAM USING THE LEAPFROG VERLET ALGORITHM.            **
C    *******************************************************************

C    *******************************************************************
C    ** THE FOLLOWING PROGRAMS WERE WRITTEN ON AN ACORN MODEL B MICRO,**
C    ** USING BBC BASIC AND SOME MACHINE CODE FOR THE GRAPHICS.       **
C    ** CONSEQUENTLY MINOR ALTERATIONS MAY BE NEEDED FOR OTHER MICROS.**
C    ** THE FIRST PROGRAM RUNS MOLECULAR DYNAMICS FROM AN INITIAL     **
C    ** CONFIGURATION WHICH MAY BE GENERATED USING THE SECOND PROGRAM.**
C    ** DEFAULT INPUT PARAMETERS ARE PROVIDED IN THE MD PROGRAM:      **
C    ** THESE CAN BE ALTERED TO SUIT THE USER.                        **
C    *******************************************************************



       MOLECULAR DYNAMICS PROGRAM

  10 IF PAGE<>&1700 THEN PAGE=&1700:CHAIN"2.FROG1"
   20 MODE1
   30 CLS:PRINT''''TAB(12);"Program Froggy"
   40PRINT'TAB(10);"Molecular Dynamics"'
   50PRINTTAB(7);" of Lennard-Jones atoms"
   60 PRINT'TAB(1);"Leapfrog algorithm with minimal storage"
   70 PRINT'TAB(7);"Repulsive WCA potential"'
   80 PRINTTAB(9);"in two dimensions"
   90 PRINT'TAB(15);"Fiche 3b"
  100 INPUTTAB(10,31)"( press <RETURN> )"A$
  110 CLS
  120*KEY0"100"
  130*KEY1"2.DATA"
  140*KEY2"0.3"
  150*KEY3"0.02"
  160 N=10:FREE=2:x1=200:y1=200:x2=512:y2=512
  170 HIMEM=&2E00
  180DIM X(12),Y(12),RX(N),RY(N),VX(N),VY(N)
  190HIMEM=&2E00:GCOL 3,4:VDU 19,2,1,0,0,0
  200 *FX138,0,128
  210NSTEP%= FNgetnum( "Enter number of steps")
  220 *FX138,0,129
  230PRINT"Enter configuration filename ";:INPUTCNFILE$
  240PRINT'TAB(5);"IN Lennard-Jones units,"
  250 *FX138,0,130
  260DENS= FNgetnum("Enter density")
  270 *FX138,0,131
  280DT= FNgetnum("Enter time step")
  290PROCreadcn(CNFILE$):SIGMA=(DENS/N)^(1/2):RCUT=SIGMA*(2.0^(1.0/6.0))
  300DT=DT*SIGMA
  310CLS
  320R=(DENS*(x2-x1)*(y2-y1)/N)^0.5/2
  330PROCcircle(R)
  340PROCcode
  350PROCforce( -DT,SIGMA,RCUT )
  360PROCmove( -DT )
  370PROCforce( -DT,SIGMA,RCUT )
  380PROCforce(  DT,SIGMA,RCUT )
  390PROCmove( DT )
  400PROCforce(DT,SIGMA,RCUT )
  410GCOL 3,1
  420PROCbox(x1-2*R,y1,x1+x2-2*R,y1+y2)
  430GCOL 3,3
  440PROCgraphics(&2F00,FALSE)
  450PROCerase
  460FOR STP%=1 TO NSTEP%
  470PROCmove( DT )
  480PROCforce(DT,SIGMA,RCUT )
  490PROCgraphics(&2F04,TRUE)
  500NEXT STP%
  510END
  520DEF FNgetnum(S$)
  530PRINT S$;" ";:INPUT A
  540=A
  550DEF PROCreadcn( FILE$ )
  560LOCAL PORT
  570PORT=OPENUP( FILE$ )
  580IF PORT=0 THEN PRINT''CHR$(34);FILE$;CHR$(34);" - File not found"':END
  590FOR I%= 1 TO N:INPUT#PORT,RX(I%),RY(I%):NEXT I%
  600FOR I%= 1 TO N:INPUT#PORT,VX(I%),VY(I%):NEXT I%
  610CLOSE #PORT
  620DEF PROCmove( DT )
  630FOR I%= 1 TO N
  640RX(I%)= RX(I%)+VX(I%)*DT
  650IF RX(I%)>0.5 REPEAT:RX(I%)=RX(I%)-1:UNTIL RX(I%)<0.5
  660IF RX(I%)<-0.5 REPEAT:RX(I%)=RX(I%)+1:UNTIL RX(I%)>-0.5
  670RY(I%)= RY(I%)+VY(I%)*DT
  680IF RY(I%)>0.5 REPEAT:RY(I%)=RY(I%)-1:UNTIL RY(I%)<0.5
  690IF RY(I%)<-0.5 REPEAT:RY(I%)=RY(I%)+1:UNTIL RY(I%)>-0.5
  700NEXT I%
  710ENDPROC
  720DEF PROCforce( DT,SIGMA,RCUT )
  730A=SIGMA*SIGMA:B=RCUT*RCUT
  740FOR I%=1 TO N-1
  750C=RX(I%):D=RY(I%):H=VX(I%):I=VY(I%)
  760FOR J%= I%+1 TO N
  770E=C-RX(J%):E=E-INT(E+0.5)
  780IF ABS(E)>=RCUT THEN 830
  790F=D-RY(J%):F=F-INT(F+0.5):G=E*E+F*F
  800IF G>=B THEN 830
  810S=A/G:T=S*S*S:U=T*T:W=24*(U+U-S):J=W*DT/G:M=J*E:O=J*F:H=H+M:I=I+O
  820VX(J%)=VX(J%)-M:VY(J%)=VY(J%)-O
  830NEXT J%
  840VX(I%)=H:VY(I%)=I
  850NEXT I%
  860ENDPROC
  870DEF PROCcircle(R)
  880D=100:T=0:X(1)=R*SIN(T):Y(1)=SQR(R^2-X(1)^2)
  890FOR A=2 TO 12
  900X(A)=X(A-1)*COS(D)+Y(A-1)*SIN(D):Y(A)=Y(A-1)*COS(D)-X(A-1)*SIN(D)
  910NEXT
  920B%=&2E00
  930FOR A%=1 TO 12
  940?B%=25:B%?1=1:B%?2=X(A%) MOD 256:B%?4=Y(A%) MOD 256
  950IF X(A%)<0 B%?3=255-(X(A%) DIV 256) ELSE B%?3=X(A%) DIV 256
  960IF Y(A%)<0 B%?5=255-(Y(A%) DIV 256) ELSE B%?5=Y(A%) DIV 256
  970B%=B%+6
  980NEXT
  990?&2E17=0
 1000ENDPROC
 1010DEF PROCcode
 1020P%=&1100
 1030FOR I%=0 TO 2 STEP2
 1040[OPT I%
 1050 .circle LDY#0
 1060.loop LDA &2E00,Y
 1070JSR &FFEE
 1080INY
 1090CPY #72
 1100BNE loop
 1110RTS
 1120.code LDX #20
 1130 LDA #&00
 1140 STA &70
 1150 LDA #&2F
 1160 STA &71
 1170 .l2 JSR getnum
 1180 JSR pl2
 1190 JSR image
 1200 INC &70
 1210 INC &70
 1220 INC &70
 1230 INC &70
 1240 DEX
 1250 BNE l2
 1260 RTS
 1270.getnum LDY#3
 1280.getnum1 LDA(&70),Y
 1290 STA &72,Y
 1300 DEY
 1310 BPL getnum1
 1320 RTS
 1330.pl2 LDA #25
 1340 JSR &FFEE
 1350 LDA #4
 1360 JSR &FFEE
 1370 LDY #0
 1380.pl21 LDA &72,Y
 1390 JSR &FFEE
 1400 INY
 1410 CPY #4
 1420 BNE pl21
 1430 JSR circle
 1440 RTS
 1450
 1460.image
 1470 LDA &73
 1480 ADC #1
 1490 STA &73
 1500 CMP #3
 1510 BPL P%+5
 1520 JSR pl2
 1530 JSR getnum
 1540 .n7 LDA &75
 1550 ADC #1
 1560 STA &75
 1570 CMP #3
 1580 BPL P%+5
 1590 JSR pl2
 1600JSR getnum
 1610LDA &75
 1620SEC
 1630SBC #2
 1640STA &75
 1650CMP #0
 1660BNE P%+11
 1670 LDA &74
 1680 CMP #&80
 1690 BMI P%+5
 1700JSR pl2
 1710JSR getnum
 1720LDA &73
 1730SEC
 1740SBC #2
 1750STA &73
 1760CMP #0
 1770BNE P%+11
 1780  LDA &72
 1790 CMP #&80
 1800 BMI P%+5
 1810 JSR pl2
 1820 RTS
 1830.go JSR code
 1840 LDY #0
 1850 .l3 LDA &2F04,Y
 1860 STA &2F00,Y
 1870 INY
 1880 CPY #81
 1890 BNE l3
 1900 RTS
 1910]
 1920NEXT
 1930ENDPROC
 1940DEF PROCgraphics(B%,run)
 1950FOR A%=1 TO 10
 1960 D%=x1+(RX(A%)+0.5)*x2:C%=y1+(RY(A%)+0.5)*y2:?(B%+(A%-1)*8)=D% MOD 256
 1970?(B%+(A%-1)*8+1)=D% DIV 256:?(B%+(A%-1)*8+2)=C% MOD 256:
 1980?(B%+(A%-1)*8+3)=C% DIV 256
 1990 NEXT
 2000 IF run CALL go
 2010ENDPROC
 2020DEF PROCerase
 2030FOR A%=1 TO 10
 2040 X%=x1+(RX(A%)+0.5)*x2
 2050 Y%=y1+(RY(A%)+0.5)*y2
 2060MOVE x1+(RX(A%)+0.5)*x2,y1+(RY(A%)+0.5)*y2
 2070CALL circle
 2080 IF X%+512<256*3 MOVEX%+512,Y%:CALL circle
 2090 IF Y%+512<256*3 MOVEX%,Y%+512:CALL circle
 2100 IF (Y%-512<256) AND (Y%-512)>128 MOVEX%,Y%-512:CALL circle
 2110 IF X%-512<256 AND X%-512>128 MOVEX%-512,Y%:CALL circle
 2120NEXT
 2130ENDPROC
 2140DEF PROCbox(A%,B%,C%,D%)
 2150MOVE A%,B%:DRAW A%,D%:DRAW C%,D%:DRAW C%,B%:DRAW A%,B%
 2160ENDPROC



      CONFIGURATION GENERATOR

   10 MODE 1
   20REM PROGRAM SETUP
   30REM
   40N= 10: REM CONSTANT
   50DIM RX(N),RY(N)
   60DIM VX(N),VY(N)
   70 A1= 3.949846138: A3= 0.252408784: A5= 0.076542912
   80 A7= 0.008355968: A9= 0.029899776
   90PROCfcc
  100 *KEY0"1.75"
  110 *KEY1"2.DATA"
  120
  130 *FX138,0,128
  140PRINT'" Enter Temperature in LJ units ";:INPUT TEMP
  150PROCcomvel( TEMP )
  160
  170 *FX138,0,129
  180PRINT'"Enter Output Data File ";:INPUT CNFILE$
  190
  200PORT= OPENOUT( CNFILE$ )
  210IF PORT=0:PRINT''"Error in opening ";CHR$(34);CNFILE$;CHR$(34)'':END
  220
  230FOR I%=1 TO N:PRINT #PORT,RX(I%),RY(I%):NEXT I%
  240FOR I%=1 TO N:PRINT #PORT,VX(I%),VY(I%):NEXT I%
  250CLOSE #PORT
  260END
  270
  280 DEF PROCfcc
  290 LOCAL X,Y,M,NC,KL,KX
  300 NC=INT((N/4)^(1/2) +0.5)
  310 X=1/NC:Y=1
  320RX(1)=0:RY(1)=0
  330RX(2)=2*Y/6:RY(2)=0
  340RX(3)=4*Y/6:RY(3)=0
  350RX(4)=Y/6:RY(4)=Y/5
  360RX(5)=3*Y/6:RY(5)=Y/5
  370RX(6)=0:RY(6)=2*Y/5
  380RX(7)=2*Y/6:RY(7)=2*Y/5
  390RX(8)=4*Y/6:RY(8)=2*Y/5
  400RX(9)=Y/6:RY(9)=3*Y/5
  410RX(10)=3*Y/6:RY(10)=3*Y/5
  420 FOR A%=1 TO 10:RX(A%)=RX(A%)-0.5:RY(A%)=RY(A%)-0.5:NEXT
  430 ENDPROC
  440DEF PROCcomvel( TEMP )
  450
  460LOCAL RTEMP,SUMX,SUMY,GAUSS,DUMMY
  470
  480RTEMP= SQR( TEMP )
  490FOR I%= 1 TO N
  500VX(I%)= RTEMP*FNgauss
  510VY(I%)= RTEMP*FNgauss
  520NEXT I%
  530SUX= 0:SUMY= 0
  540FOR I%=1 TO N
  550SUMX= SUMX+VX(I%)
  560SUMY= SUMY+VY(I%)
  570NEXT I%
  580SUMX= SUMX/N: SUMY= SUMY/N
  590FOR I%=1 TO N
  600VX(I%)= VX(I%)-SUMX
  610VY(I%)= VY(I%)-SUMY
  620NEXT I%
  630ENDPROC
  640DEF FNgauss
  650LOCAL SUM,R,R2
  660FOR L%=1 TO 12
  670SUM= SUM+RND(1)
  680NEXT L%
  690R= (SUM-6)/4: R2= R*R
  700= ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R
  710END



