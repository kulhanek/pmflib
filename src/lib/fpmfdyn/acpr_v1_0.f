CACPRRANLUX.  RANLUX, A FORTRAN IMPLEMENTATION OF THE HIGH-QUALITY      ACPR0000
C1   PSEUDORANDOM NUMBER GENERATOR OF LUSCHER.  F. JAMES.               ACPR0000
CREF. IN COMP. PHYS. COMMUN. 79 (1994) 111                              ACPR0000
      SUBROUTINE RANLUX(RVEC,LENV)                                      ACPR0001
C         Subtract-and-borrow random number generator proposed by       ACPR0002
C         Marsaglia and Zaman, implemented by F. James with the name    ACPR0003
C         RCARRY in 1991, and later improved by Martin Luescher         ACPR0004
C         in 1993 to produce "Luxury Pseudorandom Numbers".             ACPR0005
C     Fortran 77 coded by F. James, 1993                                ACPR0006
C                                                                       ACPR0007
C   LUXURY LEVELS.                                                      ACPR0008
C   ------ ------      The available luxury levels are:                 ACPR0009
C                                                                       ACPR0010
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia      ACPR0011
C           and Zaman, very long period, but fails many tests.          ACPR0012
C  level 1  (p=48): considerable improvement in quality over level 0,   ACPR0013
C           now passes the gap test, but still fails spectral test.     ACPR0014
C  level 2  (p=97): passes all known tests, but theoretically still     ACPR0015
C           defective.                                                  ACPR0016
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible         ACPR0017
C           correlations have very small chance of being observed.      ACPR0018
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.      ACPR0019
C                                                                       ACPR0020
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ACPR0021
C!!!  Calling sequences for RANLUX:                                  ++ ACPR0022
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++ ACPR0023
C!!!                   32-bit random floating point numbers between  ++ ACPR0024
C!!!                   zero (not included) and one (also not incl.). ++ ACPR0025
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++ ACPR0026
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++ ACPR0027
C!!!               which is integer between zero and MAXLEV, or if   ++ ACPR0028
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++ ACPR0029
C!!!               should be set to zero unless restarting at a break++ ACPR0030
C!!!               point given by output of RLUXAT (see RLUXAT).     ++ ACPR0031
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++ ACPR0032
C!!!               which can be used to restart the RANLUX generator ++ ACPR0033
C!!!               at the current point by calling RLUXGO.  K1 and K2++ ACPR0034
C!!!               specify how many numbers were generated since the ++ ACPR0035
C!!!               initialization with LUX and INT.  The restarting  ++ ACPR0036
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++ ACPR0037
C!!!   A more efficient but less convenient way of restarting is by: ++ ACPR0038
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++ ACPR0039
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++ ACPR0040
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++ ACPR0041
C!!!                 32-bit integer seeds, to be used for restarting ++ ACPR0042
C!!!      ISVEC must be dimensioned 25 in the calling program        ++ ACPR0043
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ACPR0044
      DIMENSION RVEC(LENV)                                              ACPR0045
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)                       ACPR0046
      PARAMETER (MAXLEV=4, LXDFLT=3)                                    ACPR0047
      DIMENSION NDSKIP(0:MAXLEV)                                        ACPR0048
      DIMENSION NEXT(24)                                                ACPR0049
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)       ACPR0050
      PARAMETER (ITWO24=2**24, ICONS=2147483563)                        ACPR0051
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV       ACPR0052
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED             ACPR0053
      INTEGER LUXLEV                                                    ACPR0054
      LOGICAL NOTYET                                                    ACPR0055
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/  ACPR0056
      DATA I24,J24,CARRY/24,10,0./                                      ACPR0057
C                               default                                 ACPR0058
C  Luxury Level   0     1     2   *3*    4                              ACPR0059
      DATA NDSKIP/0,   24,   73,  199,  365 /                           ACPR0060
Corresponds to p=24    48    97   223   389                             ACPR0061
C     time factor 1     2     3     6    10   on slow workstation       ACPR0062
C                 1    1.5    2     3     5   on fast mainframe         ACPR0063
C                                                                       ACPR0064
C  NOTYET is .TRUE. if no initialization has been performed yet.        ACPR0065
C              Default Initialization by Multiplicative Congruential    ACPR0066
      IF (NOTYET) THEN                                                  ACPR0067
         NOTYET = .FALSE.                                               ACPR0068
         JSEED = JSDFLT                                                 ACPR0069
         INSEED = JSEED                                                 ACPR0070
C         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED    ACPR0071
         LUXLEV = LXDFLT                                                ACPR0072
         NSKIP = NDSKIP(LUXLEV)                                         ACPR0073
         LP = NSKIP + 24                                                ACPR0074
         IN24 = 0                                                       ACPR0075
         KOUNT = 0                                                      ACPR0076
         MKOUNT = 0                                                     ACPR0077
C         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',    ACPR0078
C     +        LUXLEV,'      p =',LP                                     ACPR0079
            TWOM24 = 1.                                                 ACPR0080
         DO 25 I= 1, 24                                                 ACPR0081
            TWOM24 = TWOM24 * 0.5                                       ACPR0082
         K = JSEED/53668                                                ACPR0083
         JSEED = 40014*(JSEED-K*53668) -K*12211                         ACPR0084
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                         ACPR0085
         ISEEDS(I) = MOD(JSEED,ITWO24)                                  ACPR0086
   25    CONTINUE                                                       ACPR0087
         TWOM12 = TWOM24 * 4096.                                        ACPR0088
         DO 50 I= 1,24                                                  ACPR0089
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24                              ACPR0090
         NEXT(I) = I-1                                                  ACPR0091
   50    CONTINUE                                                       ACPR0092
         NEXT(1) = 24                                                   ACPR0093
         I24 = 24                                                       ACPR0094
         J24 = 10                                                       ACPR0095
         CARRY = 0.                                                     ACPR0096
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24                          ACPR0097
      ENDIF                                                             ACPR0098
C                                                                       ACPR0099
C          The Generator proper: "Subtract-with-borrow",                ACPR0100
C          as proposed by Marsaglia and Zaman,                          ACPR0101
C          Florida State University, March, 1989                        ACPR0102
C                                                                       ACPR0103
      DO 100 IVEC= 1, LENV                                              ACPR0104
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY                             ACPR0105
      IF (UNI .LT. 0.)  THEN                                            ACPR0106
         UNI = UNI + 1.0                                                ACPR0107
         CARRY = TWOM24                                                 ACPR0108
      ELSE                                                              ACPR0109
         CARRY = 0.                                                     ACPR0110
      ENDIF                                                             ACPR0111
      SEEDS(I24) = UNI                                                  ACPR0112
      I24 = NEXT(I24)                                                   ACPR0113
      J24 = NEXT(J24)                                                   ACPR0114
      RVEC(IVEC) = UNI                                                  ACPR0115
C  small numbers (with less than 12 "significant" bits) are "padded".   ACPR0116
      IF (UNI .LT. TWOM12)  THEN                                        ACPR0117
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)                    ACPR0118
C        and zero is forbidden in case someone takes a logarithm        ACPR0119
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24            ACPR0120
      ENDIF                                                             ACPR0121
C        Skipping to luxury.  As proposed by Martin Luscher.            ACPR0122
      IN24 = IN24 + 1                                                   ACPR0123
      IF (IN24 .EQ. 24)  THEN                                           ACPR0124
         IN24 = 0                                                       ACPR0125
         KOUNT = KOUNT + NSKIP                                          ACPR0126
         DO 90 ISK= 1, NSKIP                                            ACPR0127
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY                          ACPR0128
         IF (UNI .LT. 0.)  THEN                                         ACPR0129
            UNI = UNI + 1.0                                             ACPR0130
            CARRY = TWOM24                                              ACPR0131
         ELSE                                                           ACPR0132
            CARRY = 0.                                                  ACPR0133
         ENDIF                                                          ACPR0134
         SEEDS(I24) = UNI                                               ACPR0135
         I24 = NEXT(I24)                                                ACPR0136
         J24 = NEXT(J24)                                                ACPR0137
   90    CONTINUE                                                       ACPR0138
      ENDIF                                                             ACPR0139
  100 CONTINUE                                                          ACPR0140
      KOUNT = KOUNT + LENV                                              ACPR0141
      IF (KOUNT .GE. IGIGA)  THEN                                       ACPR0142
         MKOUNT = MKOUNT + 1                                            ACPR0143
         KOUNT = KOUNT - IGIGA                                          ACPR0144
      ENDIF                                                             ACPR0145
      RETURN                                                            ACPR0146
C                                                                       ACPR0147
C           Entry to input and float integer seeds from previous run    ACPR0148
      ENTRY RLUXIN(ISDEXT)                                              ACPR0149
         TWOM24 = 1.                                                    ACPR0150
         DO 195 I= 1, 24                                                ACPR0151
         NEXT(I) = I-1                                                  ACPR0152
  195    TWOM24 = TWOM24 * 0.5                                          ACPR0153
         NEXT(1) = 24                                                   ACPR0154
         TWOM12 = TWOM24 * 4096.                                        ACPR0155
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:' ACPR0156
      WRITE(6,'(5X,5I12)') ISDEXT                                       ACPR0157
      DO 200 I= 1, 24                                                   ACPR0158
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24                                 ACPR0159
  200 CONTINUE                                                          ACPR0160
      CARRY = 0.                                                        ACPR0161
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24                            ACPR0162
      ISD = IABS(ISDEXT(25))                                            ACPR0163
      I24 = MOD(ISD,100)                                                ACPR0164
      ISD = ISD/100                                                     ACPR0165
      J24 = MOD(ISD,100)                                                ACPR0166
      ISD = ISD/100                                                     ACPR0167
      IN24 = MOD(ISD,100)                                               ACPR0168
      ISD = ISD/100                                                     ACPR0169
      LUXLEV = ISD                                                      ACPR0170
        IF (LUXLEV .LE. MAXLEV) THEN                                    ACPR0171
          NSKIP = NDSKIP(LUXLEV)                                        ACPR0172
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', ACPR0173
     +                         LUXLEV                                   ACPR0174
        ELSE  IF (LUXLEV .GE. 24) THEN                                  ACPR0175
          NSKIP = LUXLEV - 24                                           ACPR0176
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV ACPR0177
        ELSE                                                            ACPR0178
          NSKIP = NDSKIP(MAXLEV)                                        ACPR0179
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV   ACPR0180
          LUXLEV = MAXLEV                                               ACPR0181
        ENDIF                                                           ACPR0182
      INSEED = -1                                                       ACPR0183
      RETURN                                                            ACPR0184
C                                                                       ACPR0185
C                    Entry to ouput seeds as integers                   ACPR0186
      ENTRY RLUXUT(ISDEXT)                                              ACPR0187
      DO 300 I= 1, 24                                                   ACPR0188
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)                        ACPR0189
  300 CONTINUE                                                          ACPR0190
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV          ACPR0191
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)                      ACPR0192
      RETURN                                                            ACPR0193
C                                                                       ACPR0194
C                    Entry to output the "convenient" restart point     ACPR0195
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)                                    ACPR0196
      LOUT = LUXLEV                                                     ACPR0197
      INOUT = INSEED                                                    ACPR0198
      K1 = KOUNT                                                        ACPR0199
      K2 = MKOUNT                                                       ACPR0200
      RETURN                                                            ACPR0201
C                                                                       ACPR0202
C                    Entry to initialize from one or three integers     ACPR0203
      ENTRY RLUXGO(LUX,INS,K1,K2)                                       ACPR0204
         IF (LUX .LT. 0) THEN                                           ACPR0205
            LUXLEV = LXDFLT                                             ACPR0206
         ELSE IF (LUX .LE. MAXLEV) THEN                                 ACPR0207
            LUXLEV = LUX                                                ACPR0208
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN                  ACPR0209
            LUXLEV = MAXLEV                                             ACPR0210
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX    ACPR0211
         ELSE                                                           ACPR0212
            LUXLEV = LUX                                                ACPR0213
            DO 310 ILX= 0, MAXLEV                                       ACPR0214
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX                ACPR0215
  310       CONTINUE                                                    ACPR0216
         ENDIF                                                          ACPR0217
      IF (LUXLEV .LE. MAXLEV)  THEN                                     ACPR0218
         NSKIP = NDSKIP(LUXLEV)                                         ACPR0219
C         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :', ACPR0220
C     +        LUXLEV,'     P=', NSKIP+24                                ACPR0221
      ELSE                                                              ACPR0222
          NSKIP = LUXLEV - 24                                           ACPR0223
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV ACPR0224
      ENDIF                                                             ACPR0225
      IN24 = 0                                                          ACPR0226
      IF (INS .LT. 0)  WRITE (6,'(A)')                                  ACPR0227
     +   ' Illegal initialization by RLUXGO, negative input seed'       ACPR0228
      IF (INS .GT. 0)  THEN                                             ACPR0229
        JSEED = INS                                                     ACPR0230
C        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS', ACPR0231
C     +      JSEED, K1,K2                                                ACPR0232
      ELSE                                                              ACPR0233
        JSEED = JSDFLT                                                  ACPR0234
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED' ACPR0235
      ENDIF                                                             ACPR0236
      INSEED = JSEED                                                    ACPR0237
      NOTYET = .FALSE.                                                  ACPR0238
      TWOM24 = 1.                                                       ACPR0239
         DO 325 I= 1, 24                                                ACPR0240
           TWOM24 = TWOM24 * 0.5                                        ACPR0241
         K = JSEED/53668                                                ACPR0242
         JSEED = 40014*(JSEED-K*53668) -K*12211                         ACPR0243
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                         ACPR0244
         ISEEDS(I) = MOD(JSEED,ITWO24)                                  ACPR0245
  325    CONTINUE                                                       ACPR0246
      TWOM12 = TWOM24 * 4096.                                           ACPR0247
         DO 350 I= 1,24                                                 ACPR0248
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24                              ACPR0249
         NEXT(I) = I-1                                                  ACPR0250
  350    CONTINUE                                                       ACPR0251
      NEXT(1) = 24                                                      ACPR0252
      I24 = 24                                                          ACPR0253
      J24 = 10                                                          ACPR0254
      CARRY = 0.                                                        ACPR0255
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24                             ACPR0256
C        If restarting at a break point, skip K1 + IGIGA*K2             ACPR0257
C        Note that this is the number of numbers delivered to           ACPR0258
C        the user PLUS the number skipped (if luxury .GT. 0).           ACPR0259
      KOUNT = K1                                                        ACPR0260
      MKOUNT = K2                                                       ACPR0261
      IF (K1+K2 .NE. 0)  THEN                                           ACPR0262
        DO 500 IOUTER= 1, K2+1                                          ACPR0263
          INNER = IGIGA                                                 ACPR0264
          IF (IOUTER .EQ. K2+1)  INNER = K1                             ACPR0265
          DO 450 ISK= 1, INNER                                          ACPR0266
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY                       ACPR0267
            IF (UNI .LT. 0.)  THEN                                      ACPR0268
               UNI = UNI + 1.0                                          ACPR0269
               CARRY = TWOM24                                           ACPR0270
            ELSE                                                        ACPR0271
               CARRY = 0.                                               ACPR0272
            ENDIF                                                       ACPR0273
            SEEDS(I24) = UNI                                            ACPR0274
            I24 = NEXT(I24)                                             ACPR0275
            J24 = NEXT(J24)                                             ACPR0276
  450     CONTINUE                                                      ACPR0277
  500   CONTINUE                                                        ACPR0278
C         Get the right value of IN24 by direct calculation             ACPR0279
        IN24 = MOD(KOUNT, NSKIP+24)                                     ACPR0280
        IF (MKOUNT .GT. 0)  THEN                                        ACPR0281
           IZIP = MOD(IGIGA, NSKIP+24)                                  ACPR0282
           IZIP2 = MKOUNT*IZIP + IN24                                   ACPR0283
           IN24 = MOD(IZIP2, NSKIP+24)                                  ACPR0284
        ENDIF                                                           ACPR0285
C       Now IN24 had better be between zero and 23 inclusive            ACPR0286
        IF (IN24 .GT. 23) THEN                                          ACPR0287
           WRITE (6,'(A/A,3I11,A,I5)')                                  ACPR0288
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,     ACPR0289
     +     K1, K2, ' cannot occur at luxury level', LUXLEV              ACPR0290
           IN24 = 0                                                     ACPR0291
        ENDIF                                                           ACPR0292
      ENDIF                                                             ACPR0293
      RETURN                                                            ACPR0294
      END                                                               ACPR0295
