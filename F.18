********************************************************************************
** FICHE F.18.  ALGORITHM FOR AVOIDING THE SQUARE ROOT OPERATION              **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        REAL FUNCTION SSQRT ( XSQ )

C    *******************************************************************
C    ** FUNCTION TO COMPUTE SQUARE ROOT FASTER THAN FORTRAN SQRT.     **
C    **                                                               **
C    ** THIS FUNCTION RETURNS THE SQUARE ROOT OF A REAL NUMBER XSQ    **
C    ** WITH 0.1 < XSQ < 1.0. THE METHOD USES A POLYNOMIAL            **
C    ** APPROXIMATION FOLLOWED BY NEWTON-RAPHSON ITERATION.           **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** SINGER, CCP5 QUARTERLY, 8, 47, 1983.                          **
C    ** POWLES, CCP5 QUARTERLY, 11, 39, 1984.                         **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    XSQ          THE NUMBER WHOSE SQUARE ROOT IS REQUIRED **
C    ** REAL    X            THE SQUARE ROOT                          **
C    ** REAL    C0,C1,C2,C3  COEFFICIENTS IN THE POLYNOMIAL APPROX.   **
C    *******************************************************************

        REAL        XSQ

        REAL        X
        REAL        C0, C1, C2, C3

        PARAMETER ( C0 =  0.188030699,  C1 = 1.48359853  )
        PARAMETER ( C2 = -1.0979059  ,  C3 = 0.430357353 )

C    *******************************************************************

C    ** POLYNOMIAL APPROXIMATION TO X **

        X = C0 + XSQ * ( C1 + XSQ * ( C2 + XSQ * C3 ) )

C    ** ITERATION OF APPROXIMATION **

        X =     X + 0.5 * ( XSQ / X - X )
        X =     X + 0.5 * ( XSQ / X - X )
        SSQRT = X + 0.5 * ( XSQ / X - X )

        RETURN
        END



