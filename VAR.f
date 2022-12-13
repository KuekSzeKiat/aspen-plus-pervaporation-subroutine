      INTEGER FUNCTION VAR(NVARS,ICOL,IROW)
      PARAMETER (IPARAM_NVARS=19)
      INTEGER IVARIABLES(8,IPARAM_NVARS), IVRSN
      DATA IVARIABLES/
     1  4HTHIC,  4HKNES,  4HS   ,  4H    ,  1,  -1,  0,  1,
     1  4HLENG,  4HTH  ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HWIDT,  4HH   ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHEIG,  4HHT_F,  4HEED ,  4H    ,  1,  -1,  0,  1,
     1  4HHEIG,  4HHT_P,  4HERME,  4HATE ,  1,  -1,  0,  1,
     1  4HREF_,  4HPERM,  4H_WAT,  4HER  ,  1,  -1,  0,  1,
     1  4HT_RE,  4HF_PE,  4HRM_W,  4HATER,  1,  -1,  0,  1,
     1  4HREF_,  4HPERM,  4H_ETH,  4H    ,  1,  -1,  0,  1,
     1  4HT_RE,  4HF_PE,  4HRM_E,  4HTH  ,  1,  -1,  0,  1,
     1  4HREF_,  4HPERM,  4H_ETL,  4H    ,  1,  -1,  0,  1,
     1  4HT_RE,  4HF_PE,  4HRM_E,  4HTL  ,  1,  -1,  0,  1,
     1  4HREF_,  4HPERM,  4H_LA ,  4H    ,  1,  -1,  0,  1,
     1  4HT_RE,  4HF_PE,  4HRM_L,  4HA   ,  1,  -1,  0,  1,
     1  4HE_AC,  4HT_W ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HE_AC,  4HT_ET,  4HH   ,  4H    ,  1,  -1,  0,  1,
     1  4HE_AC,  4HT_ET,  4HL   ,  4H    ,  1,  -1,  0,  1,
     1  4HE_AC,  4HT_LA,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HPRES,  4H_P  ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HNUM_,  4HSEGM,  4HENTS,  4H    ,  0,  -1,  0,  1/

      DATA IVRSN/1599185789/
      NVARS = IPARAM_NVARS
      IF(IROW .EQ. 0) THEN
          VAR=IVRSN
      ELSE IF(IROW .LE. NVARS) THEN
          VAR=IVARIABLES(ICOL,IROW)
      ENDIF
      RETURN
      END
