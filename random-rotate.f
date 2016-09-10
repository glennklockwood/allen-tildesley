********************************************************************************
** FICHE F.15.  ROUTINES TO RANDOMLY ROTATE MOLECULES.                        **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** THREE METHODS FOR UNIFORM RANDOM ROTATION OF LINEAR MOLECULE. **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE FLIP1 ( EXIOLD, EYIOLD, EZIOLD, DPHIMX, DCOSMX,    **
C    **                    EXINEW, EYINEW, EZINEW         )           **
C    ** SUBROUTINE FLIP2 ( EXIOLD, EYIOLD, EZIOLD, DGAMAX,            **
C    **                    EXINEW, EYINEW, EZINEW         )           **
C    ** SUBROUTINE FLIP3 ( EXIOLD, EYIOLD, EZIOLD, DOTMIN,            **
C    **                    EXINEW, EYINEW, EZINEW         )           **
C    **                                                               **
C    ** ROUTINE REQUIRED:                                             **
C    **                                                               **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    **    INCLUSIVE (SEE F.11).                                      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL      EXIOLD,EYIOLD,EZIOLD  OLD AXIAL VECTOR FOR MOL. I   **
C    ** REAL      EYINEW,EYINEW,EZINEW  TRIAL AXIAL VECTOR FOR MOL. I **
C    *******************************************************************



        SUBROUTINE FLIP1 ( EXIOLD, EYIOLD, EZIOLD, DPHIMX, DCOSMX,
     :                     EXINEW, EYINEW, EZINEW         )

C    *******************************************************************
C    ** MAKES A RANDOM CHANGE IN THE POLAR ANGLES PHI AND THETA.      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL      DCOSMX  MAXIMUM CHANGE IN COS(THETA)                **
C    ** REAL      DPHIMX  MAXIMUM CHANGE IN PHI                       **
C    ** REAL      PHIOLD  PHI IN THE OLD STATE                        **
C    ** REAL      PHINEW  PHI IN THE NEW TRIAL STATE                  **
C    ** REAL      COSOLD  COS(THETA) IN THE OLD STATE                 **
C    ** REAL      COSNEW  COS(THETA) IN THE NEW TRIAL STATE           **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FLIP1 MAKES A RANDOM CHANGE IN PHI AND COS(THETA).            **
C    ** THE MAXIMUM ALLOWED CHANGES IN THESE VARIABLES ARE CONTROLLED **
C    ** BY THE PARAMETERS DPHIMX AND DCOSMX RESPECTIVELY.             **
C    ** PHI AND THETA ARE THE EULER ANGLES DESCRIBING THE ORIENTATION **
C    ** OF THE AXIAL VECTOR. THIS METHOD CAN BE READILY EXTENDED TO   **
C    ** POLYATOMICS BY CHANGING THE THIRD EULER ANGLE , PSI, IN THE   **
C    ** RANGE ZERO TO TWOPI. IT WOULD BE FASTER TO PASS THE VARIABLES **
C    ** PHIOLD AND COSOLD DIRECTLY TO FLIP1 IF THEY ARE AVAILABLE IN  **
C    ** THE MAIN PROGRAM. SIMILARLY PHINEW AND COSNEW COULD BE        **
C    ** PASSED DIRECTLY BACK THROUGH THE SUBROUTINE HEADER            **
C    *******************************************************************

        REAL        EXIOLD, EYIOLD, EZIOLD, EXINEW, EYINEW, EZINEW
        REAL        DPHIMX, DCOSMX

        REAL        COSNEW, COSOLD, PHINEW, PHIOLD, SINNEW
        REAL        TWOPI, PI
        REAL        RANF, DUMMY

        PARAMETER ( TWOPI = 6.2831853 )

C    *******************************************************************

C    ** CONVERT THE AXIAL VECTOR TO THE EULER ANGLES **

        COSOLD = EZIOLD
        PHIOLD = ATAN2 ( EYIOLD, EXIOLD )

C    ** PERFORM THE DISPLACEMENTS **

        PHINEW = PHIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DPHIMX
        PHINEW = PHINEW - ANINT ( PHINEW / TWOPI ) * TWOPI
        COSNEW = COSOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DCOSMX
        COSNEW = COSNEW - ANINT ( COSNEW / 2.0 ) * 2.0
        SINNEW = SQRT ( 1.0 - COSNEW * COSNEW )

C    ** CONVERT THE EULER ANGLES TO AXIAL VECTORS **

        EXINEW = COS ( PHINEW ) * SINNEW
        EYINEW = SIN ( PHINEW ) * SINNEW
        EZINEW = COSNEW

        RETURN
        END



        SUBROUTINE FLIP2 ( EXIOLD, EYIOLD, EZIOLD, DGAMAX,
     :                     EXINEW, EYINEW, EZINEW         )

C    *******************************************************************
C    ** PERFORMS RANDOM ROTATION ABOUT SPACE-FIXED AXES.              **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** BARKER AND WATTS, CHEM PHYS LETTS 3, 144, 1969.               **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL      DGAMAX  MAXIMUM ANGULAR DISPLACEMENT IN RADIANS     **
C    ** INTEGER   IAXIS   SPACE FIXED AXIS FOR ROTATION               **
C    **                   (1 = X, 2 = Y, 3 = Z)                       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FLIP2 CHOOSES ONE OF THE THREE SPACE-FIXED AXES AT RANDOM     **
C    ** AND ROTATES THE MOLECULE AROUND THIS AXIS BY DGAMMA RADIANS.  **
C    ** THE MAXIMUM ANGULAR DISPLACEMENT IS DGAMAX. THIS METHOD CAN   **
C    ** READILY EXTENDED TO POLYATOMIC MOLECULES.                     **
C    *******************************************************************

        REAL    EXIOLD, EYIOLD, EZIOLD, EXINEW, EYINEW, EZINEW, DGAMAX

        REAL    COSDG, SINDG, DGAMMA
        REAL    RANF, DUMMY
        INTEGER IAXIS

C    *******************************************************************

C    ** CHOOSE A SPACE FIXED AXIS AT RANDOM **

        IAXIS = INT ( 3.0 * RANF ( DUMMY ) ) + 1

C    ** CHOOSE A RANDOM ROTATION **

        DGAMMA = ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DGAMAX

C    ** SET UP THE ROTATION MATRIX **

        COSDG = COS ( DGAMMA )
        SINDG = SIN ( DGAMMA )

C    ** PERFORM ROTATIONS **

        IF ( IAXIS .EQ. 1 ) THEN

           EXINEW = EXIOLD
           EYINEW = COSDG * EYIOLD + SINDG * EZIOLD
           EZINEW = COSDG * EZIOLD - SINDG * EYIOLD

        ELSE IF ( IAXIS .EQ. 2 ) THEN

           EXINEW = COSDG * EXIOLD - SINDG * EZIOLD
           EYINEW = EYIOLD
           EZINEW = COSDG * EZIOLD + SINDG * EXIOLD

        ELSE

           EXINEW = COSDG * EXIOLD + SINDG * EYIOLD
           EYINEW = COSDG * EYIOLD - SINDG * EXIOLD
           EZINEW = EZIOLD

        ENDIF

        RETURN
        END



        SUBROUTINE FLIP3 ( EXIOLD, EYIOLD, EZIOLD, DOTMIN,
     :                     EXINEW, EYINEW, EZINEW         )

C    *******************************************************************
C    ** CHOOSES A RANDOM DISPLACEMENT ON THE SURFACE OF A UNIT SPHERE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** MARSAGLIA, ANN MATHS STAT 43, 645, 1972.                      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL      DOT                   DOT PRODUCT OF OLD AND NEW    **
C    **                                 AXIAL VECTORS                 **
C    ** REAL      DOTMIN                PARAMETER TO ADJUST MAXIMUM   **
C    **                                 DISPLACEMENT. DOTMIN SHOULD   **
C    **                                 BE LESS THAN ONE              **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FLIP3 USES A REJECTION TECHNIQUE TO CREATE A TRIAL            **
C    ** ORIENTATION OF MOLECULE I SUBJECT TO THE CONSTRAINT THAT      **
C    ** THE COSINE OF THE ANGLE BETWEEN THE OLD AND NEW AXIAL         **
C    ** VECTORS IS GREATER THAN ( 1.0 - DOTMIN ).                     **
C    *******************************************************************

        REAL    EXIOLD, EYIOLD, EZIOLD, EXINEW, EYINEW, EZINEW, DOTMIN

        REAL    DOT, XI1, XI2, XI, XISQ
        REAL    RANF, DUMMY

C    *******************************************************************

C    ** INITIALISE DOT **

        DOT  = 0.0

C    ** ITERATIVE LOOP **

1000    IF ( ( 1.0 - DOT ) .GE. DOTMIN ) THEN

C       ** INITIALISE XISQ **

           XISQ = 1.0

C       ** INNER ITERATIVE LOOP **

2000       IF ( XISQ .GE. 1.0 ) THEN

              XI1  = RANF ( DUMMY ) * 2.0 - 1.0
              XI2  = RANF ( DUMMY ) * 2.0 - 1.0
              XISQ = XI1 * XI1 + XI2 * XI2

              GOTO 2000

           ENDIF

           XI = SQRT ( 1.0 - XISQ )
           EXINEW = 2.0 * XI1 * XI
           EYINEW = 2.0 * XI2 * XI
           EZINEW = 1.0 - 2.0 * XISQ
           DOT    = EXINEW * EXIOLD + EYINEW * EYIOLD + EZINEW * EZIOLD

           GOTO 1000

        ENDIF

        RETURN
        END



