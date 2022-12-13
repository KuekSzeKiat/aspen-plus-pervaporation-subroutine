      SUBROUTINE PERVA (NMATI, SIN, NINFI, SINFI, NMATO,
     +                SOUT, NINFO, SINFO, IDSMI, IDSII,
     +                IDSMO, IDSIO, NTOT, NSUBS, IDXSUB,
     +                ITYPE, NINT, INT, NREAL, REAL,
     +                IDS, NPO, NBOPST, NIWORK, IWORK,
     +                NWORK, WORK, NSIZE, SIZE, INTSIZ,
     +                LD )
      IMPLICIT NONE
#include "ppexec_user.cmn"
#include "dms_plex.cmn"
      Real*8 B(1)
      Equivalence (B(1),IB(1))
#include "dms_ncomp.cmn"
C     Include files pass additional variables via COMMONs.
C     Ppexec_user.cmn passes USER_NHSTRY. (1) Dms_plex.cmn passes
C     arrays containing component data such as molecular weight. (2)
C      Dms_ncomp.cmn passes NCOMP_NCC. (3)
C ------------------- Declare arguments -------------------------------------
      INTEGER NMATI, NINFI, NMATO, NINFO, NTOT,
     +        NSUBS, NINT, NPO, NIWORK, NWORK,
     +        NSIZE, NREAL
      
      INTEGER IDSMI(2,NMATI), IDSII(2,NINFI),
     +        IDSMO(2,NMATO), IDSIO(2,NINFO),
     +        IDXSUB(NSUBS), ITYPE(NSUBS), INT(NINT),
     +        IDS(2,3), NBOPST(6,NPO),
     +        IWORK(NIWORK), INTSIZ(NSIZE), LD
      REAL*8  SIN(NTOT,NMATI), SINFI(NINFI),
     +        SOUT(NTOT,NMATO), SINFO(NINFO),
     +        WORK(NWORK), SIZE(NSIZE), REAL(NREAL)
      
C --------------- Declare Local Variables ---------------------------------
      
      INTEGER OFFSET, IERR, LDATA, KDIAG, IDX(10), NCP, I, J, K, INDEX,
     +        IPERM, IRET,IFAIL, KBASE, KH, num_segments
      REAL*8  THICK, LENG, WIDTH, PERM_REF_W, T_REF_W, 
     +        PERM_REF_ETH, T_REF_ETH, E_ACT_W, E_ACT_ETH, HEI_F, 
     +        HEI_P, X(10), Y(10), FLOW, MU_IN, DIFF_IN(10), MU_P, V(10)
     +        , SA(10), CC(10), Tile(10,10), ACT_COEF(10), RC(10), 
     +        PERM_REF(10), T_REF(10), PERM(10), VELO, D_h_f, RE_f,
     +        sh_a, sh_b, sh_c, sh_d, SC(10), SH(10), K_B(10), FLUX(10),
     +        HMX_F, DHMX_F, HMX_R, DHMX_R, HMX_P, DHMX_P, E_ACT(10),
     +        Tc_w, Tr_w, HVAP_W, Tc_eth, Tr_eth, HVAP_eth, HVAP_LA,
     +        HVAP_ETL, Q_F, Q_R, Q_P, Q_vap, T_R, D_h_p, Rho_P, u_p, 
     +        RE_P, P_p, X_IN(10), F_w, F_eth, F_etl, F_la, F_f, T_f, 
     +        P_f, Rho_f, MW, T_R_cal, Flux_w, Flux_eth, Flux_etl, 
     +        Flux_la, F_p, Psat_w, Psat_eth, Psat_etl, Psat_la, 
     +        Psat(10), PERM_REF_ETL, T_REF_ETL, 
     +        PERM_REF_LA, T_REF_LA, E_ACT_ETL, E_ACT_LA

C ----------------- Declare Functions --------------------------------------

      INTEGER USRUTL_GET_REAL_PARAM, ! These functions allow access to real
     +        USRUTL_GET_INT_PARAM ! and integer parameters using named references
      
      INTEGER DMS_IFCMNC !Determines offset to universal constant data. (5)
      REAL*8 DLOG, LOG, EXP !Standard Fortran function.
      
C ----------------- Begin Executable Code ----------------------------------
C ----------------- Get configured REAL variables from Aspen Plus. ---------
      IFAIL = 0
      INDEX = 0 !Used for passing a structure. (6)

C ----------------- Get Real Variables for mass balance --------------------
      
      IERR = USRUTL_GET_REAL_PARAM('THICKNESS', INDEX, THICK) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING MEMBRANE THICKNESS'
          IFAIL = 1
      END IF

      IERR = USRUTL_GET_REAL_PARAM('LENGTH', INDEX, LENG) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING MEMBRANE LENGTH'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('WIDTH', INDEX, WIDTH) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING MEMBRANE WIDTH'
          IFAIL = 1
      END IF
     
      IERR = USRUTL_GET_REAL_PARAM('HEIGHT_FEED', INDEX, HEI_F) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING MEMBRANE HEIGHT AT FEED'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('HEIGHT_PERMEATE', INDEX, HEI_P) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING MEMBRANE HEIGHT AT PERMEATE'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('REF_PERM_WATER', INDEX, 
     +        PERM_REF_W)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING PERMEABILITY REFERENCE 
     +                    OF WATER'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('T_REF_PERM_WATER', INDEX, 
     +        T_REF_W)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING TEMPERATURE REFERENCE 
     +                        OF WATER'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('REF_PERM_ETH', INDEX, 
     +        PERM_REF_ETH)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING PERMEABILITY REFERENCE 
     +                    OF ETHANOL'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('T_REF_PERM_ETH', INDEX, 
     +        T_REF_ETH)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING TEMPERATURE REFERENCE 
     +                        OF ETHANOL'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('REF_PERM_ETL', INDEX, 
     +        PERM_REF_ETL)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING PERMEABILITY REFERENCE 
     +                    OF ETL'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('T_REF_PERM_ETL', INDEX, 
     +        T_REF_ETL)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING TEMPERATURE REFERENCE 
     +                        OF ETL'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('REF_PERM_LA', INDEX, 
     +        PERM_REF_LA)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING PERMEABILITY REFERENCE 
     +                    OF LA'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('T_REF_PERM_LA', INDEX, 
     +        T_REF_LA)
      IF (IERR .NE. 0) THEN 
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING TEMPERATURE REFERENCE 
     +                        OF LA'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('E_ACT_W', INDEX, E_ACT_W)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING W ACTIVATION ENERGY'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('E_ACT_ETH', INDEX, E_ACT_ETH)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING ETH ACTIVATION ENERGY'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('E_ACT_ETL', INDEX, E_ACT_ETL)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING ETL ACTIVATION ENERGY'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('E_ACT_LA', INDEX, E_ACT_LA)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING LA ACTIVATION ENERGY'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_REAL_PARAM('PRES_P', INDEX, P_p)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING Permeate Pressure'
          IFAIL = 1
      END IF
      
      IERR = USRUTL_GET_INT_PARAM('NUM_SEGMENTS', INDEX, num_segments)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY,*) ' ERROR FETCHING number of iterations'
          IFAIL = 1
      END IF
      
C ---------------- Calculate viscosity -------------------------------------
      CALL SHS_CPACK(SIN(1,1), NCP, IDX, X, FLOW) !Pack stream data. (9)
      KDIAG = 4
      CALL PPMON_VISCL(SIN(NCOMP_NCC+2,1), SIN(NCOMP_NCC+3,1), X, NCP,
     +                IDX, NBOPST, KDIAG, MU_IN, IERR) !Calculate viscosity,
      
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY, *) ' ERROR EVALUATING VISCOSITY FOR FEED'
      IFAIL = 1
      END IF
      IF (IFAIL .EQ. 1) RETURN
      
C ---------------- Calculate diffusion coefficient -------------------
      CALL SHS_CPACK(SIN(1,1), NCP, IDX, X, FLOW) !Pack stream data. (9)
      KDIAG = 4
      CALL PPMON_DIFCOL(SIN(NCOMP_NCC+2,1), SIN(NCOMP_NCC+3,1), X, NCP,
     +                IDX, NBOPST, KDIAG, DIFF_IN, IERR) !Calculate diffusion coefficient,
      
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY, *) ' ERROR EVALUATING DIFFUSION COEFFICEINT 
     +                         FOR FEED'
      IFAIL = 1
      END IF
      IF (IFAIL .EQ. 1) RETURN
      
C ---------------- Calculate enthalpy at feed -------------------
      
      CALL SHS_CPACK(SIN(1,1), NCP, IDX, X, FLOW)
      KDIAG=4
      KBASE = 0
      KH = 3
     
      !DHMX = J/KMOL.K
      CALL PPMON_ENTHL (SIN(NCOMP_NCC+2,1), SIN(NCOMP_NCC+3,1), X, NCP, 
     +                 IDX, NBOPST, KDIAG, KBASE, KH,
     +                 HMX_F, DHMX_F, IERR)

C --------------- Assume PERMEATE stream is first, switch if not. ----------
      IPERM = 1
      IRET = 2
      IF (IDSMO(1,1) .EQ. 'RETE') THEN !IDSMO is an argument passed to the
          IPERM = 2 !subroutine. It contains the outlet
          IRET = 1 !stream ID’s. (12)
      END IF
      
C ---------------- Model Equations -----------------------------------------
      Flux_w = 0
      Flux_eth = 0
      Flux_la = 0
      Flux_etl = 0
      F_p = 0
      y(10) = 0
      
      PERM_REF(1) = PERM_REF_W
      T_REF(1) = T_REF_W
      E_ACT(1) = E_ACT_W
      
      PERM_REF(2) = PERM_REF_ETH     ! kmol/m.s
      T_REF(2) = T_REF_ETH
      E_ACT(2) = E_ACT_ETH
      
      PERM_REF(3) = PERM_REF_ETL
      T_REF(3) = T_REF_ETL
      E_ACT(3) = E_ACT_ETL
      
      PERM_REF(4) = PERM_REF_LA     ! kmol/m.s
      T_REF(4) = T_REF_LA
      E_ACT(4) = E_ACT_LA
      
      DO I=1, 4
          X_In(I) = SIN(I,1)/SIN(5,1)
      END DO
      
      F_w = sin(1,1)
      F_eth = sin(2,1)
      F_etl = sin(3,1)
      F_la = sin(4,1)
      F_f = SIN(NCOMP_NCC+1,1)
      T_f = SIN(NCOMP_NCC+2,1)
      P_f = SIN(NCOMP_NCC+3,1)
      Rho_f = SIN(NCOMP_NCC + 8,1)
      MW = SIN(NCOMP_NCC+9,1)
      
      DO K = 1, num_segments
      WRITE(USER_NHSTRY,*) 'Number of Segments = ', K
          
          
      V(1) = 0.92/(X_IN(2) * 2.1055 + X_IN(3) * 4.4555 + X_IN(4)*3.1648)
      SA(1) = 1.4 /(X_IN(2) * 1.9720+ X_IN(3) * 3.9280 + X_IN(4)*2.88)
      CC(1) = (1-V(1) + LOG(V(1))) - (10.0 /2) * 1.4 * (1 - V(1)/SA(1)
     +        + LOG(V(1)/SA(1)))
      ACT_COEF(1) = EXP(abs(CC(1)))
      
      V(2) = 2.1055/(X_IN(1) * 0.92+ X_IN(3) * 4.4555 + X_IN(4)*3.1648)
      SA(2) = 1.9720 /(X_IN(1) * 1.4+ X_IN(3) * 3.9280 + X_IN(4)*2.88)
      CC(2) = (1-V(2) + LOG(V(2))) - (10.0 /2) * 1.972 * (1 - V(2)/SA(2)
     +        + LOG(V(2)/SA(2))) 
      ACT_COEF(2) = EXP(abs(CC(2)))
      
      V(3) = 4.4555/(X_IN(2)*2.1055 + X_IN(1) * 0.92 + X_IN(4)*3.1648)
      SA(3) = 3.9280 /(X_IN(2) * 1.9720+ X_IN(1) * 1.4 + X_IN(4)*2.88)
      CC(3) = (1-V(3) + LOG(V(3))) - (10.0 /2) * 1.4 * (1 - V(3)/SA(3)
     +        + LOG(V(3)/SA(3)))
      ACT_COEF(3) = EXP(abs(CC(3)))
      
      V(4) = 3.1648/(X_IN(1) * 0.92+ X_IN(3) * 4.4555 + X_IN(2)*2.1055)
      SA(4) = 2.88 /(X_IN(1) * 1.4+ X_IN(3) * 3.9280 + X_IN(2)*1.972)
      CC(4) = (1-V(4) + LOG(V(4))) - (10.0 /2) * 1.972 * (1 - V(4)/SA(4)
     +        + LOG(V(4)/SA(4))) 
      ACT_COEF(4) = EXP(abs(CC(4)))
      
      DO I=1,2
          PERM(I) = PERM_REF(I) * EXP((E_ACT(I)/8314.32)*(
     +     (1.0 / T_REF(I)) - (1.0 / T_f)))
      END DO
      
      VELO = F_f * MW / 
     +       (Rho_f * HEI_F * WIDTH) ! m/s
      
      D_h_f = 4 * (HEI_F * WIDTH) / (2 * ( HEI_F + WIDTH))    !m
      
      RE_F = Rho_F * VELO * D_H_F / MU_IN      !unitless
      
      IF (RE_F < 2300) THEN
          sh_a = 1.615
          sh_b = 0.33
          sh_c = 0.33
          sh_d = 0.33 
      ELSE IF (RE_F >= 2300) THEN
          sh_a = 0.026
          sh_b = 0.80
          sh_c = 0.30
          sh_d = 0
      END IF
      
      Psat_w = (10**(7.0436 - (1636.91/ ((T_f-273.15)+224.92))))*1000
      Psat_eth = (10**(7.1688 - (1552.60/ ((T_f- 273.15)+222.42))))*1000
      Psat_etl = (10**(7.8269 - (2489.7/ ((T_f- 273.15)+273.15))))*1000
      Psat_la = (10**(7.2471 - (1968.21/ ((T_f- 273.15)+158.94))))*1000
      
      Psat(1) = Psat_w
      Psat(2) = Psat_eth
      Psat(3) = Psat_etl
      Psat(4) = Psat_la
      
      DO I=1, 4                                                      !1 IS WATER 2 IS OTHER
          SC(I) = MU_IN / (Rho_F * DIFF_IN(i))       !unitless
                                                                      !!! LEFT PERMEA AND ACTIVITY COEFF
          SH(i) =sh_a *(RE_F**sh_b) * (SC(i)**sh_c)*(D_H_F / 
     +             (LENG/num_segments))**sh_d
      
          K_B(i) = SH(i) * DIFF_IN(i) / D_H_F
      
          FLUX(i) = (1 ./ (( ACT_COEF(i) / (k_b(i) * 
     +     (Rho_f/ MW)))
     +     + (THICK / PERM(i)))) 
     +     * ((X(I)*ACT_COEF(I)) - y(I)*(P_p/Psat(I)))*
     +      (Width*(LENG/num_segments))
      END DO
      
      Write(USER_NHSTRY,*) ' flux =  ', flux(1), ',', flux(2),',', 
     +                        Flux(3),',', Flux(4)
      
      SOUT(1,1) = F_w - FLUX(1)
      SOUT(2,1) = F_eth - FLUX(2)
      SOUT(3,1) = F_etl - FLUX(3)
      SOUT(4,1) = F_la - FLUX(4)
      SOUT(5,1) = F_f - FLUX(1)-FLUX(2)- FLUX(3)-FLUX(4)
      
      SOUT(1,2) = Flux_w + FLUX(1)
      SOUT(2,2) = Flux_eth + FLUX(2)
      SOUT(3,2) = Flux_etl + FLUX(3)
      SOUT(4,2) = Flux_la + FLUX(4)
      SOUT(5,2) = F_p + FLUX(1) + FLUX(2)+ FLUX(3)+FLUX(4)
      
C
C      
C------------- ENERGY BALANCE----------------------------------------
      
C-------- HEAT OF VAPORIZATION (T in heat of vaporisation want loop to get T?)
      Tc_w= 647.13    !K      
      Tc_eth= 513.92    !K   
      
      T_r = T_f
      
      Write(USER_NHSTRY,*) ' t_r =  ', T_r
      Write(USER_NHSTRY,*) ' P_f =  ', P_f
      
      HVAP_LA = 6.91 * (10**7)    !J/kmol
      HVAP_ETL = 5.177 * (10**7)
      
      DO
          Tr_w = T_f / Tc_w
          HVAP_W = 5.2053 * (10**7) * (1-Tr_w)**(0.3199 - 0.212*Tr_w + 
     +         0.25795*Tr_w*Tr_w)
      
          Tr_eth = T_f / Tc_eth
          HVAP_ETH =5.69 * (10**7) * (1-Tr_eth)**(0.3359)
          
          CALL SHS_CPACK(SOUT(1,1), NCP, IDX, X, FLOW)
          KDIAG=4
          KBASE = 0
          KH = 3
          CALL PPMON_ENTHL(T_R, P_f, X, NCP, 
     +                 IDX, NBOPST, KDIAG, KBASE, KH,
     +                 HMX_R, DHMX_R, IERR)
          IF (IERR .NE. 0) THEN
          WRITE(USER_NHSTRY, *) ' ERROR EVALUATING enth FOR retentate'
          IFAIL = 1
          END IF
          IF (IFAIL .EQ. 1) RETURN
          
          CALL SHS_CPACK(SOUT(1,2), NCP, IDX, Y, FLOW)
          KDIAG=4
          KBASE = 0
          KH = 3
          CALL PPMON_ENTHV(T_r, P_p, Y, NCP, 
     +                 IDX, NBOPST, KDIAG, KBASE, KH,
     +                 HMX_P, DHMX_P, IERR)
          IF (IERR .NE. 0) THEN
          WRITE(USER_NHSTRY, *) ' ERROR EVALUATING enth FOR perm'
          IFAIL = 1
          END IF
          IF (IFAIL .EQ. 1) RETURN
          
          Q_F = F_f * DHMX_F * T_f  !J/s
          Q_R = SOUT(NCOMP_NCC+1,1) * DHMX_R !J/s.K
          Q_P = SOUT(NCOMP_NCC+1,2) * DHMX_P                      !J/s.K
          Q_vap = FLUX(1) * HVAP_W + FLUX(2) * HVAP_ETH + FLUX(3) 
     +             * HVAP_ETL + FLUX(4) * HVAP_LA      !J/s

          T_R_cal = (Q_F - Q_vap ) / (Q_R + Q_P)
      
          IF (abs((T_r-T_R_cal)/T_r_cal)<0.01) THEN
              SOUT(6,1) = T_R_cal
              SOUT(6,2) = T_R_cal
              Write(USER_NHSTRY,*) ' t_r1 =  ', T_r_cal
              EXIT
          Else
              T_r = T_r_cal
              Write(USER_NHSTRY,*) ' t_r2 =  ', T_r_cal
          END IF
          
      END DO
      
C
C ---------------- Calculate permeate viscosity ---------------------------
      CALL SHS_CPACK(SOUT(1,2), NCP, IDX, Y, FLOW) !Pack stream data. (9)
      KDIAG = 4
      CALL PPMON_VISCV(T_r, P_p, Y, NCP,
     +                IDX, NBOPST, KDIAG, MU_P, IERR)
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY, *) ' ERROR EVALUATING VISCOSITY FOR Permeate'
      IFAIL = 1
      END IF
      IF (IFAIL .EQ. 1) RETURN
      
C ---------------- Permeate Pressure ------------------------------     
      D_h_p = 4 * (HEI_P * WIDTH) / (2 * ( HEI_P + WIDTH))
      
      Rho_P = y(1) * 997 + y(2) * 789 + y(3) * 1030 + y(4) * 1206
      
      RE_P = Rho_P * u_p * D_H_P / MU_P
      
      IF (RE_P < 2320) THEN                                          
         P_p =  ((P_p)**2 - 19 * 8314.32 * 
     +          T_r * MU_P* SOUT(NCOMP_NCC+1,2) * 
     +        (LENG/num_segments) /(D_h_p**2 * WIDTH * HEI_P))**(0.5)
      ELSE IF (RE_P >= 2320) THEN
         P_p =  ((P_p)**2 - 0.444 * 8314.32 *
     +          T_r * MU_P * 
     +        (WIDTH*LENG)**0.748 * SOUT(NCOMP_NCC+1,2)**1.748  * 
     +        (LENG/num_segments) /
     +          (D_h_p**1.252 * WIDTH**1.748 * HEI_P**1.748 ))**(0.5) 
      END IF
      
      SOUT(7,1) = SIN(7,1)
      SOUT(7,2) = P_p
      
C-----------Data for iterations-----------------
      F_w = SOUT(1,1)
      F_eth = SOUT(2,1)
      F_etl = SOUT(3,1)
      F_la = SOUT(4,1)
      F_f = SOUT(5,1)
      
      Flux_w = SOUT(1,2)
      Flux_eth = SOUT(2,2)
      Flux_etl = SOUT(3,2)
      Flux_la = SOUT(4,2)
      F_p = SOUT(5,2)
      T_f = T_r_cal
      DO I=1, 4
          X_In(I) = SOUT(I,1)/SOUT(5,1)
      END DO
      
      CALL SHS_CPACK(SOUT(1,1), NCP, IDX, X, FLOW)
      KDIAG = 4
      CALL PPMON_VISCL(T_r, P_f, X, NCP,
     +                IDX, NBOPST, KDIAG, MU_IN, IERR)
      
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY, *) ' ERROR EVALUATING VISCOSITY FOR FEED'
      IFAIL = 1
      END IF
      IF (IFAIL .EQ. 1) RETURN
      
      CALL SHS_CPACK(SOUT(1,1), NCP, IDX, X, FLOW)
      KDIAG = 4
      CALL PPMON_DIFCOL(T_r, P_f, X, NCP,
     +                IDX, NBOPST, KDIAG, DIFF_IN, IERR) 
      IF (IERR .NE. 0) THEN
      WRITE(USER_NHSTRY, *) ' ERROR EVALUATING DIFFUSION COEFFICEINT 
     +                         FOR FEED'
      IFAIL = 1
      END IF
      IF (IFAIL .EQ. 1) RETURN
      
      CALL SHS_CPACK(SOUT(1,1), NCP, IDX, X, FLOW)
      KDIAG=4
      KBASE = 0
      KH = 3
     
      !DHMX = J/KMOL.K
      CALL PPMON_ENTHL (T_f,P_f, X, NCP, 
     +                 IDX, NBOPST, KDIAG, KBASE, KH,
     +                 HMX_F, DHMX_F, IERR)
      
      MW = X_in(1) * 18.01528 + X_in(2) * 46.07 +X_in(3) * 118.132 + 
     +    X_in(4) * 90.08
      Rho_f = X_in(1) * 997 + X_in(2) * 789 + x_in(3) * 1030 + 
     +        x_in(4) * 1206
      
      END DO
      
      RETURN
       END