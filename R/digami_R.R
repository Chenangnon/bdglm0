#C     ALGORITHM AS 187  APPL. STATIST. (1982) VOL.31, NO.3
#C
#C     Computes derivatives of the incomplete gamma integral for positive
#C     parameters, X, P, using a series expansion if P > X or X <= 1, and
#C     a continued fraction expansion otherwise.
#C
#C     Calculation of D(4) in line 60 corrected 5 October 1993.
#C
#C     N.B. The user must input values of the incomplete gamma, digamma
#C          and trigamma functions.  These can be obtained using AS 239
#C          (or 32), AS 103 and AS 121 respectively.
#C
.digami_R <-
  function (X, P, GPLOG, GP1LOG, PSIP,
           PSIP1, PSIDP, PSIDP1,
           E = 1.e-6, OFLO = 1.E300,
           TMAX = 100, VSMALL = 1.E-30) {
    # C
    PN <- D <- DP <- DPP <- numeric(6)
    ZERO = 0
    ONE = 1
    TWO = 2

    IFAULT = 0
    #C
    #C     Derivatives with respect to X
    #C
    PM1 = P - ONE
    XLOG = logb(X)
    D[1] = exp(-GPLOG + PM1*XLOG - X)
    D[2] = D[1] * (PM1/X - ONE)
    D[5] = D[1] * (XLOG - PSIP)
    #C
    #C     Derivatives with respect to P
    #C
    #IF (X .GT. ONE .AND. X .GE. P) GO TO 30
    if (X <= ONE | X <= P) {
      #C
      #C     Series expansion
      #C
      F = exp(P*XLOG - GP1LOG - X)
      DFP = F * (XLOG - PSIP1)
      DFPP = DFP*DFP/F - F*PSIDP1
      #C
      TMAXP = TMAX + P
      C = S = ONE
      CP = CPP = DSP = DSPP = ZERO
      for (A in (P + 1):(TMAXP + 1)) {
        CPC = CP / C
        CP = CPC - ONE/A
        CPP = CPP/C - CPC*CPC + ONE/A^2
        C = C*X/A
        CP = CP*C
        CPP = CPP*C + CP*CP/C
        S = S + C
        DSP = DSP + CP
        DSPP = DSPP + CPP
        if (A > TMAXP) {
          IFAULT = 1
          break
        }
        if (C <= E*S)
          break
      }
      if (IFAULT)
        return(structure(D, fail = 1))
      D[6] = S*F
      D[3] = S*DFP + F*DSP
      D[4] = S*DFPP + TWO*DFP*DSP + F*DSPP
      return(structure(D, fail = 0))
    }
    #C
    #C     Continued fraction expansion
    #C
    F = exp(P * XLOG - GPLOG - X)
    DFP = F * (XLOG - PSIP)
    DFPP = DFP * DFP/F - F*PSIDP
    #C
    A = PM1
    B = X + ONE - A
    TERM = ZERO
    PN[1] = ONE
    PN[2] = X
    PN[3] = X + ONE
    PN[4] = X * B
    S0 = PN[3] / PN[4]
    # DO 31 I = 1, 4
    for (I in 1:4) { # REMOVE
      DP[I] = ZERO
      DPP[I] = ZERO
    }
    # 31 CONTINUE
    DP[4] = -X
    #C
    #32
    for (TERM in 1:(TMAX + 1)) {
      A = A - ONE
      B = B + TWO
      #TERM = TERM + ONE
      AN = A*TERM
      PN[5] = B*PN[3] + AN*PN[1]
      PN[6] = B*PN[4] + AN*PN[2]
      DP[5] = B*DP[3] - PN[3] + AN*DP[1] + PN[1] *TERM
      DP[6] = B*DP[4] - PN[4] + AN*DP[2] + PN[2] *TERM
      DPP[5] = B*DPP[3] + AN*DPP[1] + TWO*(TERM*DP[1] - DP[3])
      DPP[6] = B*DPP[4] + AN*DPP[2] + TWO*(TERM*DP[2] - DP[4])
      #C
      if (abs(PN[6]) >= VSMALL) {
        S = PN[5] / PN[6]
        C = abs(S - S0)
        if (C*P <= E & C < E*S)
          break
        #C
        #34
        S0 = S
      }
      # 35 DO 36 I = 1, 4
      for (I in 1:4) {
        I2 = I + 2
        DP[I] = DP[I2]
        DPP[I] = DPP[I2]
        PN[I] = PN[I2]
      }
      # 36 CONTINUE
      #C
      if (TERM > TMAX)
        return(structure(D, fail = 0))
      if (abs(PN[5]) >= OFLO) {
        for (I in 1:4) {
          DP[I] = DP[I] / OFLO
          DPP[I] = DPP[I] / OFLO
          PN[I] = PN[I] / OFLO
        }
      }
      #41 CONTINUE
      #GO TO 32
    }
    #C
    #42
    D[6] = ONE - F*S
    DSP = (DP[5] - S*DP[6]) / PN[6]
    DSPP = (DPP[5] - S*DPP[6] - TWO*DSP*DP[6]) / PN[6]
    D[3] = -F*DSP - S*DFP
    D[4] = -F*DSPP - TWO*DSP*DFP - S*DFPP
    return(structure(D, fail = 0))
    #C
    #C     Set fault indicator
    #C
    #1001 IFAULT = 1
    #RETURN
    #END
  }
