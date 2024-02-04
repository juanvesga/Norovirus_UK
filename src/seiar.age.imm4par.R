## Definition of the time-step and output as "time"
dt <- user(1)
steps_per_week<- 7/dt
steps_per_year<- 365/dt
initial(time) <- 0
update(time) <- (step + 1) * dt


# index 1 is GI3 infection
# index 2 is other Gi infection
# index 3 is GII4 infection
# index 4 is other GII infection
# For states E, I and A the first number is the active infecyive strain
# the following numbers are the carrying immunity
# R compartments show carrying immunity in ascending order

## Equations for transitions between compartments by age group
###############################################################
# Not infected
update(M[]) <-  
  M[i] + 
  n_bM[i] - 
  n_ageoM[i] + 
  n_ageiM[i] - 
  n_MS[i] - 
  n_Mvacc[i] - 
  n_muM[i] 

update(G[]) <-   
  G[i] +
  n_Gvacc[i]+ 
  n_bG[i] - 
  n_ageoG[i] + 
  n_ageiG[i] - 
  n_muG[i] -
  n_Gvacc[i]

update(S[]) <-
  S[i] + 
  n_bS[i] - 
  n_ageoS[i] + 
  n_ageiS[i] + 
  n_MS[i] + 
  n_R1_S[i] + 
  n_R2_S[i] + 
  n_R3_S[i] +
  n_R4_S[i] +
  n_R13_wane[i]+
  n_R12_wane[i]+  
  n_R14_wane[i] +
  n_R23_wane[i]+
  n_R24_wane[i]+
  n_R34_wane[i]+
  n_R123_wane[i] +
  n_R134_wane[i]+
  n_R124_wane[i]+
  n_R234_wane[i]+
  n_R1234_wane[i] -
  n_S_E1[i] - 
  n_S_E2[i] -
  n_S_E3[i] - 
  n_S_E4[i] -
  n_Svacc[i] - 
  n_muS[i]


# infected(GI.3) (strain 1)
update(E1[]) <-
  E1[i]- 
  n_ageoE1[i] + 
  n_ageiE1[i] + 
  n_S_E1[i] - 
  n_E1_I1[i] - 
  n_muE1[i]

update(I1[]) <-   
  I1[i] - 
  n_ageoI1[i] + 
  n_ageiI1[i] + 
  n_E1_I1[i] -
  n_I1_A1[i] - 
  n_muI1[i] 

update(A1[]) <-  
  A1[i]- 
  n_ageoA1[i] + 
  n_ageiA1[i] + 
  n_I1_A1[i] + 
  n_R1_A1[i] - 
  n_A1_R1[i] - 
  n_muA1[i]

update(R1[]) <- 
  R1[i] - 
  n_ageoR1[i] + 
  n_ageiR1[i] + 
  n_A1_R1[i] - 
  n_R1_S[i] -
  n_R1_A1[i] -
  n_R1_E21[i] - 
  n_R1_E31[i] - 
  n_R1_E41[i] - 
  n_muR1[i] -
  n_R1vacc[i]


# infected(other GI) (strain 2)
update(E2[]) <-
  E2[i]- 
  n_ageoE2[i] + 
  n_ageiE2[i] + 
  n_S_E2[i] - 
  n_E2_I2[i] - 
  n_muE2[i]

update(I2[]) <-   
  I2[i] - 
  n_ageoI2[i] + 
  n_ageiI2[i] + 
  n_E2_I2[i] -
  n_I2_A2[i] - 
  n_muI2[i] 

update(A2[]) <-  
  A2[i]- 
  n_ageoA2[i] + 
  n_ageiA2[i] + 
  n_I2_A2[i] + 
  n_R2_A2[i] - 
  n_A2_R2[i] - 
  n_muA2[i]

update(R2[]) <- 
  R2[i] - 
  n_ageoR2[i] + 
  n_ageiR2[i] + 
  n_A2_R2[i]- 
  n_R2_S[i] -
  n_R2_A2[i] -
  n_R2_E12[i] - 
  n_R2_E32[i] - 
  n_R2_E42[i] - 
  n_muR2[i] -
  n_R2vacc[i]


# infected(GII.4) (strain 3)
update(E3[]) <-
  E3[i]- 
  n_ageoE3[i] + 
  n_ageiE3[i] + 
  n_S_E3[i] - 
  n_E3_I3[i] - 
  n_muE3[i]

update(I3[]) <-   
  I3[i] - 
  n_ageoI3[i] + 
  n_ageiI3[i] + 
  n_E3_I3[i] -
  n_I3_A3[i] - 
  n_muI3[i] 

update(A3[]) <-  
  A3[i]- 
  n_ageoA3[i] + 
  n_ageiA3[i] + 
  n_I3_A3[i] + 
  n_R3_A3[i] - 
  n_A3_R3[i] - 
  n_muA3[i]

update(R3[]) <- 
  R3[i] - 
  n_ageoR3[i] + 
  n_ageiR3[i] + 
  n_A3_R3[i] - 
  n_R3_S[i] -
  n_R3_A3[i] -
  n_R3_E13[i] - 
  n_R3_E23[i] - 
  n_R3_E43[i] - 
  n_muR3[i] -
  n_R3vacc[i]

# infected(other GII) (strain 4)
update(E4[]) <-
  E4[i]- 
  n_ageoE4[i] + 
  n_ageiE4[i] + 
  n_S_E4[i] - 
  n_E4_I4[i] - 
  n_muE4[i]

update(I4[]) <-   
  I4[i] - 
  n_ageoI4[i] + 
  n_ageiI4[i] + 
  n_E4_I4[i] -
  n_I4_A4[i] - 
  n_muI4[i] 

update(A4[]) <-  
  A4[i]- 
  n_ageoA4[i] + 
  n_ageiA4[i] + 
  n_I4_A4[i] + 
  n_R4_A4[i] - 
  n_A4_R4[i] - 
  n_muA4[i]

update(R4[]) <- 
  R4[i] - 
  n_ageoR4[i] + 
  n_ageiR4[i] + 
  n_A4_R4[i] - 
  n_R4_S[i] -
  n_R4_A4[i] -
  n_R4_E14[i] - 
  n_R4_E24[i] - 
  n_R4_E34[i] - 
  n_muR4[i] -
  n_R4vacc[i]



# infected(GII4) (strain 3) from s1
update(E31[]) <-
  E31[i]- 
  n_ageoE31[i] + 
  n_ageiE31[i] +
  n_R1_E31[i] - 
  n_E31_I31[i] - 
  n_muE31[i]

update(I31[]) <-   
  I31[i] - 
  n_ageoI31[i] + 
  n_ageiI31[i] + 
  n_E31_I31[i] -
  n_I31_A31[i] - 
  n_muI31[i] 

update(A31[]) <-  
  A31[i]- 
  n_ageoA31[i] + 
  n_ageiA31[i] + 
  n_I31_A31[i] + 
  n_R13_A31[i] - 
  n_A31_R13[i] - 
  n_muA31[i]

# infected(Other GI) (strain 2) from s1
update(E21[]) <-
  E21[i]- 
  n_ageoE21[i] + 
  n_ageiE21[i] +
  n_R1_E21[i]  - 
  n_E21_I21[i] - 
  n_muE21[i]

update(I21[]) <-   
  I21[i] - 
  n_ageoI21[i] + 
  n_ageiI21[i] + 
  n_E21_I21[i] -
  n_I21_A21[i] - 
  n_muI21[i] 

update(A21[]) <-  
  A21[i]- 
  n_ageoA21[i] + 
  n_ageiA21[i] + 
  n_I21_A21[i] + 
  n_R12_A21[i] - 
  n_A21_R12[i] - 
  n_muA21[i]

# infected(Other GII) (strain 4) from s1
update(E41[]) <-
  E41[i]- 
  n_ageoE41[i] + 
  n_ageiE41[i] +
  n_R1_E41[i]  - 
  n_E41_I41[i] - 
  n_muE41[i]

update(I41[]) <-   
  I41[i] - 
  n_ageoI41[i] + 
  n_ageiI41[i] + 
  n_E41_I41[i] -
  n_I41_A41[i] - 
  n_muI41[i] 

update(A41[]) <-  
  A41[i]- 
  n_ageoA41[i] + 
  n_ageiA41[i] + 
  n_I41_A41[i] + 
  n_R14_A41[i] - 
  n_A41_R14[i] - 
  n_muA41[i]


# infected(GI3) (strain 1) from s2
update(E12[]) <-
  E12[i]- 
  n_ageoE12[i] + 
  n_ageiE12[i] +
  n_R2_E12[i] - 
  n_E12_I12[i] - 
  n_muE12[i]

update(I12[]) <-   
  I12[i] - 
  n_ageoI12[i] + 
  n_ageiI12[i] + 
  n_E12_I12[i] -
  n_I12_A12[i] - 
  n_muI12[i] 

update(A12[]) <-  
  A12[i]- 
  n_ageoA12[i] + 
  n_ageiA12[i] + 
  n_I12_A12[i] + 
  n_R12_A12[i] - 
  n_A12_R12[i] - 
  n_muA12[i]

# infected(GII4) (strain 3) from s2
update(E32[]) <-
  E32[i]- 
  n_ageoE32[i] + 
  n_ageiE32[i] +
  n_R2_E32[i] - 
  n_E32_I32[i] - 
  n_muE32[i]

update(I32[]) <-   
  I32[i] - 
  n_ageoI32[i] + 
  n_ageiI32[i] + 
  n_E32_I32[i] -
  n_I32_A32[i] - 
  n_muI32[i] 

update(A32[]) <-  
  A32[i]- 
  n_ageoA32[i] + 
  n_ageiA32[i] + 
  n_I32_A32[i] + 
  n_R23_A32[i] - 
  n_A32_R23[i] - 
  n_muA32[i]

# infected(O GII) (strain 4) from s2
update(E42[]) <-
  E42[i]- 
  n_ageoE42[i] + 
  n_ageiE42[i] +
  n_R2_E42[i] - 
  n_E42_I42[i] - 
  n_muE42[i]

update(I42[]) <-   
  I42[i] - 
  n_ageoI42[i] + 
  n_ageiI42[i] + 
  n_E42_I42[i] -
  n_I42_A42[i] - 
  n_muI42[i] 

update(A42[]) <-  
  A42[i]- 
  n_ageoA42[i] + 
  n_ageiA42[i] + 
  n_I42_A42[i] + 
  n_R24_A42[i] - 
  n_A42_R24[i] - 
  n_muA42[i]


# infected(GI3) (strain 1) from s3
update(E13[]) <-
  E13[i]- 
  n_ageoE13[i] + 
  n_ageiE13[i] +
  n_R3_E13[i] - 
  n_E13_I13[i] - 
  n_muE13[i]

update(I13[]) <-   
  I13[i] - 
  n_ageoI13[i] + 
  n_ageiI13[i] + 
  n_E13_I13[i] -
  n_I13_A13[i] - 
  n_muI13[i] 

update(A13[]) <-  
  A13[i]- 
  n_ageoA13[i] + 
  n_ageiA13[i] + 
  n_I13_A13[i] + 
  n_R13_A13[i] - 
  n_A13_R13[i] - 
  n_muA13[i]

# infected(O-GI) (strain 2) from s3
update(E23[]) <-
  E23[i]- 
  n_ageoE23[i] + 
  n_ageiE23[i] +
  n_R3_E23[i] - 
  n_E23_I23[i] - 
  n_muE23[i]

update(I23[]) <-   
  I23[i] - 
  n_ageoI23[i] + 
  n_ageiI23[i] + 
  n_E23_I23[i] -
  n_I23_A23[i] - 
  n_muI23[i] 

update(A23[]) <-  
  A23[i]- 
  n_ageoA23[i] + 
  n_ageiA23[i] + 
  n_I23_A23[i] + 
  n_R23_A23[i] - 
  n_A23_R23[i] - 
  n_muA23[i]

# infected(O-GII) (strain 4) from s3
update(E43[]) <-
  E43[i]- 
  n_ageoE43[i] + 
  n_ageiE43[i] +
  n_R3_E43[i] - 
  n_E43_I43[i] - 
  n_muE43[i]

update(I43[]) <-   
  I43[i] - 
  n_ageoI43[i] + 
  n_ageiI43[i] + 
  n_E43_I43[i] -
  n_I43_A43[i] - 
  n_muI43[i] 

update(A43[]) <-  
  A43[i]- 
  n_ageoA43[i] + 
  n_ageiA43[i] + 
  n_I43_A43[i] + 
  n_R34_A43[i] - 
  n_A43_R34[i] - 
  n_muA43[i]


# infected(GI3) (strain 1) from s4
update(E14[]) <-
  E14[i]- 
  n_ageoE14[i] + 
  n_ageiE14[i] +
  n_R4_E14[i] - 
  n_E14_I14[i] - 
  n_muE14[i]

update(I14[]) <-   
  I14[i] - 
  n_ageoI14[i] + 
  n_ageiI14[i] + 
  n_E14_I14[i] -
  n_I14_A14[i] - 
  n_muI14[i] 

update(A14[]) <-  
  A14[i]- 
  n_ageoA14[i] + 
  n_ageiA14[i] + 
  n_I14_A14[i] + 
  n_R14_A14[i] - 
  n_A14_R14[i] - 
  n_muA14[i]

# infected(O-GI) (strain 2) from s4
update(E24[]) <-
  E24[i]- 
  n_ageoE24[i] + 
  n_ageiE24[i] +
  n_R4_E24[i] - 
  n_E24_I24[i] - 
  n_muE24[i]

update(I24[]) <-   
  I24[i] - 
  n_ageoI24[i] + 
  n_ageiI24[i] + 
  n_E24_I24[i] -
  n_I24_A24[i] - 
  n_muI24[i] 

update(A24[]) <-  
  A24[i]- 
  n_ageoA24[i] + 
  n_ageiA24[i] + 
  n_I24_A24[i] + 
  n_R24_A24[i] - 
  n_A24_R24[i] - 
  n_muA24[i]

# infected(GII4) (strain 3) from s4
update(E34[]) <-
  E34[i]- 
  n_ageoE34[i] + 
  n_ageiE34[i] +
  n_R4_E34[i] - 
  n_E34_I34[i] - 
  n_muE34[i]

update(I34[]) <-   
  I34[i] - 
  n_ageoI34[i] + 
  n_ageiI34[i] + 
  n_E34_I34[i] -
  n_I34_A34[i] - 
  n_muI34[i] 

update(A34[]) <-  
  A34[i]- 
  n_ageoA34[i] + 
  n_ageiA34[i] + 
  n_I34_A34[i] + 
  n_R34_A34[i] - 
  n_A34_R34[i] - 
  n_muA34[i]


###### 2 STRAIN IMMUNITY Rs 

update(R13[]) <- 
  R13[i] - 
  n_ageoR13[i] + 
  n_ageiR13[i] + 
  n_A13_R13[i] + 
  n_A31_R13[i] +
  n_toR13vacc[i] - 
  n_R13_wane[i] -
  n_R13_A13[i] -
  n_R13_A31[i] -
  n_R13_E213[i] - 
  n_R13_E413[i] - 
  n_muR13[i] -
  n_R13vacc[i]

update(R12[]) <- 
  R12[i] - 
  n_ageoR12[i] + 
  n_ageiR12[i] + 
  n_A12_R12[i] + 
  n_A21_R12[i] - 
  n_R12_wane[i] -
  n_R12_A12[i] -
  n_R12_A21[i] -
  n_R12_E312[i] - 
  n_R12_E412[i] - 
  n_muR12[i] -
  n_R12vacc[i]

update(R23[]) <- 
  R23[i] - 
  n_ageoR23[i] + 
  n_ageiR23[i] + 
  n_A23_R23[i] + 
  n_A32_R23[i]  - 
  n_R23_wane[i] -
  n_R23_A23[i] -
  n_R23_A32[i] -
  n_R23_E123[i] - 
  n_R23_E423[i] - 
  n_muR23[i] -
  n_R23vacc[i]

update(R14[]) <- 
  R14[i] - 
  n_ageoR14[i] + 
  n_ageiR14[i] + 
  n_A14_R14[i] + 
  n_A41_R14[i]  - 
  n_R14_wane[i] -
  n_R14_A14[i] -
  n_R14_A41[i] -
  n_R14_E214[i] - 
  n_R14_E314[i] - 
  n_muR14[i] -
  n_R14vacc[i]

update(R24[]) <- 
  R24[i] - 
  n_ageoR24[i] + 
  n_ageiR24[i] + 
  n_A24_R24[i] + 
  n_A42_R24[i]  - 
  n_R24_wane[i] -
  n_R24_A24[i] -
  n_R24_A42[i] -
  n_R24_E124[i] - 
  n_R24_E324[i] - 
  n_muR24[i]-
  n_R24vacc[i]

update(R34[]) <- 
  R34[i] - 
  n_ageoR34[i] + 
  n_ageiR34[i] + 
  n_A34_R34[i] + 
  n_A43_R34[i]  - 
  n_R34_wane[i] -
  n_R34_A34[i] -
  n_R34_A43[i] -
  n_R34_E134[i] - 
  n_R34_E234[i] - 
  n_muR34[i] -
  n_R34vacc[i]


# infected(O-GI) (strain 2) from r13

update(E213[]) <-
  E213[i]- 
  n_ageoE213[i] + 
  n_ageiE213[i] +
  n_R13_E213[i] - 
  n_E213_I213[i] - 
  n_muE213[i]

update(I213[]) <-   
  I213[i] - 
  n_ageoI213[i] + 
  n_ageiI213[i] + 
  n_E213_I213[i] -
  n_I213_A213[i] - 
  n_muI213[i] 

update(A213[]) <-  
  A213[i]- 
  n_ageoA213[i] + 
  n_ageiA213[i] + 
  n_I213_A213[i] + 
  n_R123_A213[i] - 
  n_A213_R123[i] - 
  n_muA213[i]


# infected(O-GII) (strain 4) from r13

update(E413[]) <-
  E413[i]- 
  n_ageoE413[i] + 
  n_ageiE413[i] +
  n_R13_E413[i] - 
  n_E413_I413[i] - 
  n_muE413[i]

update(I413[]) <-   
  I413[i] - 
  n_ageoI413[i] + 
  n_ageiI413[i] + 
  n_E413_I413[i] -
  n_I413_A413[i] - 
  n_muI413[i] 

update(A413[]) <-  
  A413[i]- 
  n_ageoA413[i] + 
  n_ageiA413[i] + 
  n_I413_A413[i] + 
  n_R134_A413[i] - 
  n_A413_R134[i] - 
  n_muA413[i]

# infected(GII4) (strain 3) from r12

update(E312[]) <-
  E312[i]- 
  n_ageoE312[i] + 
  n_ageiE312[i] +
  n_R12_E312[i] - 
  n_E312_I312[i] - 
  n_muE312[i]

update(I312[]) <-   
  I312[i] - 
  n_ageoI312[i] + 
  n_ageiI312[i] + 
  n_E312_I312[i] -
  n_I312_A312[i] - 
  n_muI312[i] 

update(A312[]) <-  
  A312[i]- 
  n_ageoA312[i] + 
  n_ageiA312[i] + 
  n_I312_A312[i] + 
  n_R123_A312[i] - 
  n_A312_R123[i] - 
  n_muA312[i]

# infected(O-GII) (strain 4) from r12

update(E412[]) <-
  E412[i]- 
  n_ageoE412[i] + 
  n_ageiE412[i] +
  n_R12_E412[i] - 
  n_E412_I412[i] - 
  n_muE412[i]

update(I412[]) <-   
  I412[i] - 
  n_ageoI412[i] + 
  n_ageiI412[i] + 
  n_E412_I412[i] -
  n_I412_A412[i] - 
  n_muI412[i] 

update(A412[]) <-  
  A412[i]- 
  n_ageoA412[i] + 
  n_ageiA412[i] + 
  n_I412_A412[i] + 
  n_R124_A412[i] - 
  n_A412_R124[i] - 
  n_muA412[i]

## infected(GI3) (strain 1) from r23

update(E123[]) <-
  E123[i]- 
  n_ageoE123[i] + 
  n_ageiE123[i] +
  n_R23_E123[i] - 
  n_E123_I123[i] - 
  n_muE123[i]

update(I123[]) <-   
  I123[i] - 
  n_ageoI123[i] + 
  n_ageiI123[i] + 
  n_E123_I123[i] -
  n_I123_A123[i] - 
  n_muI123[i] 

update(A123[]) <-  
  A123[i]- 
  n_ageoA123[i] + 
  n_ageiA123[i] + 
  n_I123_A123[i] + 
  n_R123_A123[i] - 
  n_A123_R123[i] - 
  n_muA123[i]


# infected(O-GII) (strain 4) from r23

update(E423[]) <-
  E423[i]- 
  n_ageoE423[i] + 
  n_ageiE423[i] +
  n_R23_E423[i] - 
  n_E423_I423[i] - 
  n_muE423[i]

update(I423[]) <-   
  I423[i] - 
  n_ageoI423[i] + 
  n_ageiI423[i] + 
  n_E423_I423[i] -
  n_I423_A423[i] - 
  n_muI423[i] 

update(A423[]) <-  
  A423[i]- 
  n_ageoA423[i] + 
  n_ageiA423[i] + 
  n_I423_A423[i] + 
  n_R234_A423[i] - 
  n_A423_R234[i] - 
  n_muA423[i]


# infected(GII4) (strain 3) from r14

update(E314[]) <-
  E314[i]- 
  n_ageoE314[i] + 
  n_ageiE314[i] +
  n_R14_E314[i] - 
  n_E314_I314[i] - 
  n_muE314[i]

update(I314[]) <-   
  I314[i] - 
  n_ageoI314[i] + 
  n_ageiI314[i] + 
  n_E314_I314[i] -
  n_I314_A314[i] - 
  n_muI314[i] 

update(A314[]) <-  
  A314[i]- 
  n_ageoA314[i] + 
  n_ageiA314[i] + 
  n_I314_A314[i] + 
  n_R134_A314[i] - 
  n_A314_R134[i] - 
  n_muA314[i]

# infected(O-GI) (strain 2) from r14

update(E214[]) <-
  E214[i]- 
  n_ageoE214[i] + 
  n_ageiE214[i] +
  n_R14_E214[i] - 
  n_E214_I214[i] - 
  n_muE214[i]

update(I214[]) <-   
  I214[i] - 
  n_ageoI214[i] + 
  n_ageiI214[i] + 
  n_E214_I214[i] -
  n_I214_A214[i] - 
  n_muI214[i] 

update(A214[]) <-  
  A214[i]- 
  n_ageoA214[i] + 
  n_ageiA214[i] + 
  n_I214_A214[i] + 
  n_R124_A214[i] - 
  n_A214_R124[i] - 
  n_muA214[i]


# infected(GI3) (strain 1) from r24

update(E124[]) <-
  E124[i]- 
  n_ageoE124[i] + 
  n_ageiE124[i] +
  n_R24_E124[i] - 
  n_E124_I124[i] - 
  n_muE124[i]

update(I124[]) <-   
  I124[i] - 
  n_ageoI124[i] + 
  n_ageiI124[i] + 
  n_E124_I124[i] -
  n_I124_A124[i] - 
  n_muI124[i] 

update(A124[]) <-  
  A124[i]- 
  n_ageoA124[i] + 
  n_ageiA124[i] + 
  n_I124_A124[i] + 
  n_R124_A124[i] - 
  n_A124_R124[i] - 
  n_muA124[i]

# infected(GII4) (strain 3) from r24

update(E324[]) <-
  E324[i]- 
  n_ageoE324[i] + 
  n_ageiE324[i] +
  n_R24_E324[i] - 
  n_E324_I324[i] - 
  n_muE324[i]

update(I324[]) <-   
  I324[i] - 
  n_ageoI324[i] + 
  n_ageiI324[i] + 
  n_E324_I324[i] -
  n_I324_A324[i] - 
  n_muI324[i] 

update(A324[]) <-  
  A324[i]- 
  n_ageoA324[i] + 
  n_ageiA324[i] + 
  n_I324_A324[i] + 
  n_R234_A324[i] - 
  n_A324_R234[i] - 
  n_muA324[i]

# infected(GI3) (strain 1) from r34

update(E134[]) <-
  E134[i]- 
  n_ageoE134[i] + 
  n_ageiE134[i] +
  n_R34_E134[i] - 
  n_E134_I134[i] - 
  n_muE134[i]

update(I134[]) <-   
  I134[i] - 
  n_ageoI134[i] + 
  n_ageiI134[i] + 
  n_E134_I134[i] -
  n_I134_A134[i] - 
  n_muI134[i] 

update(A134[]) <-  
  A134[i]- 
  n_ageoA134[i] + 
  n_ageiA134[i] + 
  n_I134_A134[i] + 
  n_R134_A134[i] - 
  n_A134_R134[i] - 
  n_muA134[i]

# infected(O-GI) (strain 2) from r34

update(E234[]) <-
  E234[i]- 
  n_ageoE234[i] + 
  n_ageiE234[i] +
  n_R34_E234[i] - 
  n_E234_I234[i] - 
  n_muE234[i]

update(I234[]) <-   
  I234[i] - 
  n_ageoI234[i] + 
  n_ageiI234[i] + 
  n_E234_I234[i] -
  n_I234_A234[i] - 
  n_muI234[i] 

update(A234[]) <-  
  A234[i]- 
  n_ageoA234[i] + 
  n_ageiA234[i] + 
  n_I234_A234[i] + 
  n_R234_A234[i] - 
  n_A234_R234[i] - 
  n_muA234[i]


###### 3 STRAIN IMMUNITY Rs 
###########################

update(R123[]) <- 
  R123[i] - 
  n_ageoR123[i] + 
  n_ageiR123[i] + 
  n_A123_R123[i] + 
  n_A213_R123[i] +
  n_A312_R123[i] +
  n_toR123vacc[i] -
  n_R123_wane[i] -
  n_R123_A123[i] -
  n_R123_A213[i] -
  n_R123_A312[i] -
  n_R123_E4123[i] - 
  n_muR123[i] -
  n_R123vacc[i]

update(R124[]) <- 
  R124[i] - 
  n_ageoR124[i] + 
  n_ageiR124[i] + 
  n_A124_R124[i] + 
  n_A214_R124[i] +
  n_A412_R124[i]  -
  n_R124_wane[i] -
  n_R124_A124[i] -
  n_R124_A214[i] -
  n_R124_A412[i] -
  n_R124_E3124[i] - 
  n_muR124[i] -
  n_R124vacc[i]

update(R134[]) <- 
  R134[i] - 
  n_ageoR134[i] + 
  n_ageiR134[i] + 
  n_A134_R134[i] + 
  n_A314_R134[i] +
  n_A413_R134[i] +
  n_toR134vacc[i] -
  n_R134_wane[i] -
  n_R134_A134[i] -
  n_R134_A314[i] -
  n_R134_A413[i] -
  n_R134_E2134[i] - 
  n_muR134[i] -
  n_R134vacc[i]


update(R234[]) <- 
  R234[i] - 
  n_ageoR234[i] + 
  n_ageiR234[i] + 
  n_A234_R234[i] + 
  n_A324_R234[i] +
  n_A423_R234[i]  -
  n_R234_wane[i] -
  n_R234_A234[i] -
  n_R234_A324[i] -
  n_R234_A423[i] -
  n_R234_E1234[i] - 
  n_muR234[i] -
  n_R234vacc[i]

# infected(GI3) (strain 1) from r234

update(E1234[]) <-
  E1234[i]- 
  n_ageoE1234[i] + 
  n_ageiE1234[i] +
  n_R234_E1234[i] - 
  n_E1234_I1234[i] - 
  n_muE1234[i]

update(I1234[]) <-   
  I1234[i] - 
  n_ageoI1234[i] + 
  n_ageiI1234[i] + 
  n_E1234_I1234[i] -
  n_I1234_A1234[i] - 
  n_muI1234[i] 

update(A1234[]) <-  
  A1234[i]- 
  n_ageoA1234[i] + 
  n_ageiA1234[i] + 
  n_I1234_A1234[i] + 
  n_R1234_A1234[i] - 
  n_A1234_R1234[i] - 
  n_muA1234[i]


# infected(O-GI) (strain 2) from r134

update(E2134[]) <-
  E2134[i]- 
  n_ageoE2134[i] + 
  n_ageiE2134[i] +
  n_R134_E2134[i] - 
  n_E2134_I2134[i] - 
  n_muE2134[i]

update(I2134[]) <-   
  I2134[i] - 
  n_ageoI2134[i] + 
  n_ageiI2134[i] + 
  n_E2134_I2134[i] -
  n_I2134_A2134[i] - 
  n_muI2134[i] 

update(A2134[]) <-  
  A2134[i]- 
  n_ageoA2134[i] + 
  n_ageiA2134[i] + 
  n_I2134_A2134[i] + 
  n_R1234_A2134[i] - 
  n_A2134_R1234[i] - 
  n_muA2134[i]

# infected(GII4) (strain 3) from r124

update(E3124[]) <-
  E3124[i]- 
  n_ageoE3124[i] + 
  n_ageiE3124[i] +
  n_R124_E3124[i] - 
  n_E3124_I3124[i] - 
  n_muE3124[i]

update(I3124[]) <-   
  I3124[i] - 
  n_ageoI3124[i] + 
  n_ageiI3124[i] + 
  n_E3124_I3124[i] -
  n_I3124_A3124[i] - 
  n_muI3124[i] 

update(A3124[]) <-  
  A3124[i]- 
  n_ageoA3124[i] + 
  n_ageiA3124[i] + 
  n_I3124_A3124[i] + 
  n_R1234_A3124[i] - 
  n_A3124_R1234[i] - 
  n_muA3124[i]

# infected(O-GII) (strain 4) from r123

update(E4123[]) <-
  E4123[i]- 
  n_ageoE4123[i] + 
  n_ageiE4123[i] +
  n_R123_E4123[i] - 
  n_E4123_I4123[i] - 
  n_muE4123[i]

update(I4123[]) <-   
  I4123[i] - 
  n_ageoI4123[i] + 
  n_ageiI4123[i] + 
  n_E4123_I4123[i] -
  n_I4123_A4123[i] - 
  n_muI4123[i] 

update(A4123[]) <-  
  A4123[i]- 
  n_ageoA4123[i] + 
  n_ageiA4123[i] + 
  n_I4123_A4123[i] + 
  n_R1234_A4123[i] - 
  n_A4123_R1234[i] - 
  n_muA4123[i]


### 4 STRAIN IMUNITY R 
#######################

update(R1234[]) <- 
  R1234[i] - 
  n_ageoR1234[i] + 
  n_ageiR1234[i] + 
  n_A1234_R1234[i] + 
  n_A2134_R1234[i] +
  n_A3124_R1234[i] +
  n_A4123_R1234[i] + 
  n_toR1234vacc[i] -
  n_R1234_wane[i] -
  n_R1234_A1234[i] -
  n_R1234_A2134[i] -
  n_R1234_A3124[i] -
  n_R1234_A4123[i] -
  n_muR1234[i] -
  n_R1234vacc[i]


####### Outputs
#########################
reported_fac[1]<-1/repfac_0
reported_fac[2]<-1/repfac_5
reported_fac[3]<-1/repfac_15
reported_fac[4]<-1/repfac_65p
dim(reported_fac)<-4

reported_fac_long[1:4]<-1/repfac_0
reported_fac_long[5:8]<-1/repfac_5
reported_fac_long[9:13]<-1/repfac_15
reported_fac_long[N_age]<-1/repfac_65p
dim(reported_fac_long)<-N_age


# Daily infections incidence
update(infections_day_gi3[]) <- 
  if (i==1)
    (infections_day_gi3[i]+
       sum(n_E1_I1[1:4]) * reported_fac[i] +
       sum(n_E12_I12[1:4]) * reported_fac[i] +
       sum(n_E13_I13[1:4]) * reported_fac[i] +
       sum(n_E14_I14[1:4]) * reported_fac[i] +
       sum(n_E123_I123[1:4]) * reported_fac[i] +
       sum(n_E124_I124[1:4]) * reported_fac[i] +
       sum(n_E134_I134[1:4]) * reported_fac[i] +
       sum(n_E1234_I1234[1:4]) * reported_fac[i] ) else
         if(i==2) 
           (infections_day_gi3[i]+
              sum(n_E1_I1[5:8]) * reported_fac[i] +
              sum(n_E12_I12[5:8]) * reported_fac[i] +
              sum(n_E13_I13[5:8]) * reported_fac[i] +
              sum(n_E14_I14[5:8]) * reported_fac[i] +
              sum(n_E123_I123[5:8]) * reported_fac[i] +
              sum(n_E124_I124[5:8]) * reported_fac[i] +
              sum(n_E134_I134[5:8]) * reported_fac[i] +
              sum(n_E1234_I1234[5:8]) * reported_fac[i] ) else
                if(i==3) 
                  (infections_day_gi3[i]+
                     sum(n_E1_I1[9:13]) * reported_fac[i] +
                     sum(n_E12_I12[9:13]) * reported_fac[i] +
                     sum(n_E13_I13[9:13]) * reported_fac[i] +
                     sum(n_E14_I14[9:13]) * reported_fac[i] +
                     sum(n_E123_I123[9:13]) * reported_fac[i] +
                     sum(n_E124_I124[9:13]) * reported_fac[i] +
                     sum(n_E134_I134[9:13]) * reported_fac[i] +
                     sum(n_E1234_I1234[9:13]) * reported_fac[i] ) else
                       (infections_day_gi3[i]+
                          (n_E1_I1[14]) * reported_fac[i] +
                          (n_E12_I12[14]) * reported_fac[i] +
                          (n_E13_I13[14]) * reported_fac[i] +
                          (n_E14_I14[14]) * reported_fac[i] +
                          (n_E123_I123[14]) * reported_fac[i] +
                          (n_E124_I124[14]) * reported_fac[i] +
                          (n_E134_I134[14]) * reported_fac[i] +
                          (n_E1234_I1234[14]) * reported_fac[i] )


update(infections_day_gi[]) <- 
  if (i==1)
    (infections_day_gi[i]+
       sum(n_E2_I2[1:4]) * reported_fac[i] +
       sum(n_E21_I21[1:4]) * reported_fac[i] +
       sum(n_E23_I23[1:4]) * reported_fac[i] +
       sum(n_E24_I24[1:4]) * reported_fac[i] +
       sum(n_E213_I213[1:4]) * reported_fac[i] +
       sum(n_E214_I214[1:4]) * reported_fac[i] +
       sum(n_E234_I234[1:4]) * reported_fac[i] +
       sum(n_E2134_I2134[1:4]) * reported_fac[i]) else
         if (i==2)
           (infections_day_gi[i]+
              sum(n_E2_I2[5:8]) * reported_fac[i] +
              sum(n_E21_I21[5:8]) * reported_fac[i] +
              sum(n_E23_I23[5:8]) * reported_fac[i] +
              sum(n_E24_I24[5:8]) * reported_fac[i] +
              sum(n_E213_I213[5:8]) * reported_fac[i] +
              sum(n_E214_I214[5:8]) * reported_fac[i] +
              sum(n_E234_I234[5:8]) * reported_fac[i] +
              sum(n_E2134_I2134[5:8]) * reported_fac[i]) else
                if (i==3)
                  (infections_day_gi[i]+
                     sum(n_E2_I2[9:13]) * reported_fac[i] +
                     sum(n_E21_I21[9:13]) * reported_fac[i] +
                     sum(n_E23_I23[9:13]) * reported_fac[i] +
                     sum(n_E24_I24[9:13]) * reported_fac[i] +
                     sum(n_E213_I213[9:13]) * reported_fac[i] +
                     sum(n_E214_I214[9:13]) * reported_fac[i] +
                     sum(n_E234_I234[9:13]) * reported_fac[i] +
                     sum(n_E2134_I2134[9:13]) * reported_fac[i]) else
                       (infections_day_gi[i]+
                          (n_E2_I2[14]) * reported_fac[i] +
                          (n_E21_I21[14]) * reported_fac[i] +
                          (n_E23_I23[14]) * reported_fac[i] +
                          (n_E24_I24[14]) * reported_fac[i] +
                          (n_E213_I213[14]) * reported_fac[i] +
                          (n_E214_I214[14]) * reported_fac[i] +
                          (n_E234_I234[14]) * reported_fac[i] +
                          (n_E2134_I2134[14]) * reported_fac[i])


update(infections_day_gii4[]) <- 
  if(i==1)
    (infections_day_gii4[i] +
       sum(n_E3_I3[1:4]) * reported_fac[i] +
       sum(n_E31_I31[1:4]) * reported_fac[i] +
       sum(n_E32_I32[1:4]) * reported_fac[i] +
       sum(n_E34_I34[1:4]) * reported_fac[i] +
       sum(n_E312_I312[1:4]) * reported_fac[i] +
       sum(n_E314_I314[1:4]) * reported_fac[i] +
       sum(n_E324_I324[1:4]) * reported_fac[i] +
       sum(n_E3124_I3124[1:4]) * reported_fac[i])else
         if(i==2)
           (infections_day_gii4[i] +
              sum(n_E3_I3[5:8]) * reported_fac[i] +
              sum(n_E31_I31[5:8]) * reported_fac[i] +
              sum(n_E32_I32[5:8]) * reported_fac[i] +
              sum(n_E34_I34[5:8]) * reported_fac[i] +
              sum(n_E312_I312[5:8]) * reported_fac[i] +
              sum(n_E314_I314[5:8]) * reported_fac[i] +
              sum(n_E324_I324[5:8]) * reported_fac[i] +
              sum(n_E3124_I3124[5:8]) * reported_fac[i])else
                if(i==3)
                  (infections_day_gii4[i] +
                     sum(n_E3_I3[9:13]) * reported_fac[i] +
                     sum(n_E31_I31[9:13]) * reported_fac[i] +
                     sum(n_E32_I32[9:13]) * reported_fac[i] +
                     sum(n_E34_I34[9:13]) * reported_fac[i] +
                     sum(n_E312_I312[9:13]) * reported_fac[i] +
                     sum(n_E314_I314[9:13]) * reported_fac[i] +
                     sum(n_E324_I324[9:13]) * reported_fac[i] +
                     sum(n_E3124_I3124[9:13]) * reported_fac[i])else
                       (infections_day_gii4[i] +
                          (n_E3_I3[14]) * reported_fac[i] +
                          (n_E31_I31[14]) * reported_fac[i] +
                          (n_E32_I32[14]) * reported_fac[i] +
                          (n_E34_I34[14]) * reported_fac[i] +
                          (n_E312_I312[14]) * reported_fac[i] +
                          (n_E314_I314[14]) * reported_fac[i] +
                          (n_E324_I324[14]) * reported_fac[i] +
                          (n_E3124_I3124[14]) * reported_fac[i])



update(infections_day_gii[]) <- 
  if(i==1)
    (infections_day_gii[i] +
       sum(n_E4_I4[1:4] ) * reported_fac[i] +
       sum(n_E41_I41[1:4] ) * reported_fac[i] +
       sum(n_E42_I42[1:4] ) * reported_fac[i] +
       sum(n_E43_I43[1:4] ) * reported_fac[i] +
       sum(n_E412_I412[1:4] ) * reported_fac[i] +
       sum(n_E413_I413[1:4] ) * reported_fac[i] +
       sum(n_E423_I423[1:4] ) * reported_fac[i] +
       sum(n_E4123_I4123[1:4] ) * reported_fac[i]) else
         if(i==2)
           (infections_day_gii[i] +
              sum(n_E4_I4[5:8] ) * reported_fac[i] +
              sum(n_E41_I41[5:8] ) * reported_fac[i] +
              sum(n_E42_I42[5:8] ) * reported_fac[i] +
              sum(n_E43_I43[5:8] ) * reported_fac[i] +
              sum(n_E412_I412[5:8] ) * reported_fac[i] +
              sum(n_E413_I413[5:8] ) * reported_fac[i] +
              sum(n_E423_I423[5:8] ) * reported_fac[i] +
              sum(n_E4123_I4123[5:8] ) * reported_fac[i]) else
                if(i==3)
                  (infections_day_gii[i] +
                     sum(n_E4_I4[9:13] ) * reported_fac[i] +
                     sum(n_E41_I41[9:13] ) * reported_fac[i] +
                     sum(n_E42_I42[9:13] ) * reported_fac[i] +
                     sum(n_E43_I43[9:13] ) * reported_fac[i] +
                     sum(n_E412_I412[9:13] ) * reported_fac[i] +
                     sum(n_E413_I413[9:13] ) * reported_fac[i] +
                     sum(n_E423_I423[9:13] ) * reported_fac[i] +
                     sum(n_E4123_I4123[9:13] ) * reported_fac[i]) else
                       (infections_day_gii[i] +
                          (n_E4_I4[14] ) * reported_fac[i] +
                          (n_E41_I41[14] ) * reported_fac[i] +
                          (n_E42_I42[14] ) * reported_fac[i] +
                          (n_E43_I43[14] ) * reported_fac[i] +
                          (n_E412_I412[14] ) * reported_fac[i] +
                          (n_E413_I413[14] ) * reported_fac[i] +
                          (n_E423_I423[14] ) * reported_fac[i] +
                          (n_E4123_I4123[14] ) * reported_fac[i])



# Number of vaccines per day

update(vaccines_perday) <- 
  sum(n_Mvacc) +
  sum(n_Gvacc) +
  sum(n_toR13vacc)+
  sum(n_toR123vacc)+
  sum(n_toR134vacc)+
  sum(n_toR1234vacc)


# Weekly reported cases (match sgss)

update(reported_wk[]) <- if (step %% steps_per_week == 0)
  (   n_E1_I1[i] * reported_fac_long[i] +
        n_E12_I12[i] * reported_fac_long[i] +
        n_E13_I13[i] * reported_fac_long[i] +
        n_E14_I14[i] * reported_fac_long[i] +
        n_E123_I123[i] * reported_fac_long[i] +
        n_E124_I124[i] * reported_fac_long[i] +
        n_E134_I134[i] * reported_fac_long[i] +
        n_E1234_I1234[i] * reported_fac_long[i] +
        n_E2_I2[i] * reported_fac_long[i] +
        n_E21_I21[i] * reported_fac_long[i] +
        n_E23_I23[i] * reported_fac_long[i] +
        n_E24_I24[i] * reported_fac_long[i] +
        n_E213_I213[i] * reported_fac_long[i] +
        n_E214_I214[i] * reported_fac_long[i] +
        n_E234_I234[i] * reported_fac_long[i] +
        n_E2134_I2134[i] * reported_fac_long[i] +
        n_E3_I3[i] * reported_fac_long[i] +
        n_E31_I31[i] * reported_fac_long[i] +
        n_E32_I32[i] * reported_fac_long[i] +
        n_E34_I34[i] * reported_fac_long[i] +
        n_E312_I312[i] * reported_fac_long[i] +
        n_E314_I314[i] * reported_fac_long[i] +
        n_E324_I324[i] * reported_fac_long[i] +
        n_E3124_I3124[i] * reported_fac_long[i] +
        n_E4_I4[i] * reported_fac_long[i] +
        n_E41_I41[i] * reported_fac_long[i] +
        n_E42_I42[i] * reported_fac_long[i] +
        n_E43_I43[i] * reported_fac_long[i] +
        n_E412_I412[i] * reported_fac_long[i] +
        n_E413_I413[i] * reported_fac_long[i] +
        n_E423_I423[i] * reported_fac_long[i] +
        n_E4123_I4123[i] * reported_fac_long[i]) else  
          ( reported_wk[i] + 
              ( n_E1_I1[i] * reported_fac_long[i] +
                  n_E12_I12[i] * reported_fac_long[i] +
                  n_E13_I13[i] * reported_fac_long[i] +
                  n_E14_I14[i] * reported_fac_long[i] +
                  n_E123_I123[i] * reported_fac_long[i] +
                  n_E124_I124[i] * reported_fac_long[i] +
                  n_E134_I134[i] * reported_fac_long[i] +
                  n_E1234_I1234[i]* reported_fac_long[i] +
                  n_E2_I2[i] * reported_fac_long[i] +
                  n_E21_I21[i] * reported_fac_long[i] +
                  n_E23_I23[i] * reported_fac_long[i] +
                  n_E24_I24[i] * reported_fac_long[i] +
                  n_E213_I213[i] * reported_fac_long[i] +
                  n_E214_I214[i] * reported_fac_long[i] +
                  n_E234_I234[i] * reported_fac_long[i] +
                  n_E2134_I2134[i]* reported_fac_long[i] +
                  n_E3_I3[i] * reported_fac_long[i] +
                  n_E31_I31[i] * reported_fac_long[i] +
                  n_E32_I32[i] * reported_fac_long[i] +
                  n_E34_I34[i] * reported_fac_long[i] +
                  n_E312_I312[i] * reported_fac_long[i] +
                  n_E314_I314[i] * reported_fac_long[i] +
                  n_E324_I324[i] * reported_fac_long[i] +
                  n_E3124_I3124[i]* reported_fac_long[i] +
                  n_E4_I4[i] * reported_fac_long[i] +
                  n_E41_I41[i] * reported_fac_long[i] +
                  n_E42_I42[i] * reported_fac_long[i] +
                  n_E43_I43[i] * reported_fac_long[i] +
                  n_E412_I412[i] * reported_fac_long[i] +
                  n_E413_I413[i] * reported_fac_long[i] +
                  n_E423_I423[i] * reported_fac_long[i] +
                  n_E4123_I4123[i]* reported_fac_long[i] ) ) 



update(inc_day_gi3[])<-
  (if(i==1) sum(n_E1_I1[1]) + sum(n_E12_I12[1]) +
     sum(n_E13_I13[1]) + sum(n_E14_I14[1]) + 
     sum(n_E123_I123[1]) + sum(n_E124_I124[1]) + 
     sum(n_E134_I134[1]) + sum(n_E1234_I1234[1])  else 
       if(i==2)  sum(n_E1_I1[2:4]) + sum(n_E12_I12[2:4]) +
     sum(n_E13_I13[2:4]) + sum(n_E14_I14[2:4]) + 
     sum(n_E123_I123[2:4]) + sum(n_E124_I124[2:4]) +
     sum(n_E134_I134[2:4]) + sum(n_E1234_I1234[2:4])  else 
       if(i==3)  
         sum(n_E1_I1[5:8]) + sum(n_E12_I12[5:8]) + sum(n_E13_I13[5:8]) +
     sum(n_E14_I14[5:8]) + sum(n_E123_I123[5:8]) + sum(n_E124_I124[5:8]) +
     sum(n_E134_I134[5:8]) + sum(n_E1234_I1234[5:8]) else 
       if (i==4) 
         sum(n_E1_I1[9:13]) + sum(n_E12_I12[9:13]) + sum(n_E13_I13[9:13]) +
     sum(n_E14_I14[9:13]) + sum(n_E123_I123[9:13]) + sum(n_E124_I124[9:13]) +
     sum(n_E134_I134[9:13]) + sum(n_E1234_I1234[9:13])  
   else
     sum(n_E1_I1[14]) + sum(n_E12_I12[14]) + sum(n_E13_I13[14]) +
     sum(n_E14_I14[14]) + sum(n_E123_I123[14]) + sum(n_E124_I124[14]) +
     sum(n_E134_I134[14]) + sum(n_E1234_I1234[14])     ) 


update(inc_day_gi[])<-
  (if(i==1) sum(n_E2_I2[1]) +
     sum(n_E21_I21[1]) +  sum(n_E23_I23[1]) +  sum(n_E24_I24[1]) +
     sum(n_E213_I213[1]) + sum(n_E214_I214[1]) + sum(n_E234_I234[1]) +
     sum(n_E2134_I2134[1])  else
       if(i==2)  sum(n_E2_I2[2:4]) +
     sum(n_E21_I21[2:4]) + sum(n_E23_I23[2:4]) +
     sum(n_E24_I24[2:4]) + sum(n_E213_I213[2:4]) +
     sum(n_E214_I214[2:4]) + sum(n_E234_I234[2:4]) +
     sum(n_E2134_I2134[2:4])  else 
       if(i==3)  sum(n_E2_I2[5:8]) +
     sum(n_E21_I21[5:8]) + sum(n_E23_I23[5:8]) +
     sum(n_E24_I24[5:8]) + sum(n_E213_I213[5:8]) +
     sum(n_E214_I214[5:8]) + sum(n_E234_I234[5:8]) +
     sum(n_E2134_I2134[5:8])  else 
       if (i==4) sum(n_E2_I2[9:13]) +
     sum(n_E21_I21[9:13]) + sum(n_E23_I23[9:13]) +
     sum(n_E24_I24[9:13]) + sum(n_E213_I213[9:13]) +
     sum(n_E214_I214[9:13]) + sum(n_E234_I234[9:13]) +
     sum(n_E2134_I2134[9:13])  else
       sum(n_E2_I2[14]) + sum(n_E21_I21[14]) +
     sum(n_E23_I23[14]) + sum(n_E24_I24[14]) +
     sum(n_E213_I213[14]) + sum(n_E214_I214[14]) +
     sum(n_E234_I234[14]) + sum(n_E2134_I2134[14]) ) 


update(inc_day_gii4[])<-
  (if(i==1)   sum(n_E3_I3[1]) +
     sum(n_E31_I31[1]) + sum(n_E32_I32[1]) +
     sum(n_E34_I34[1]) + sum(n_E312_I312[1]) +
     sum(n_E314_I314[1]) + sum(n_E324_I324[1]) +
     sum(n_E3124_I3124[1]) else
       if(i==2)    sum(n_E3_I3[2:4]) +
     sum(n_E31_I31[2:4]) + sum(n_E32_I32[2:4]) +
     sum(n_E34_I34[2:4]) + sum(n_E312_I312[2:4]) +
     sum(n_E314_I314[2:4]) + sum(n_E324_I324[2:4]) +
     sum(n_E3124_I3124[2:4])  else 
       if(i==3)    sum(n_E3_I3[5:8]) +
     sum(n_E31_I31[5:8]) + sum(n_E32_I32[5:8]) +
     sum(n_E34_I34[5:8]) + sum(n_E312_I312[5:8]) +
     sum(n_E314_I314[5:8]) + sum(n_E324_I324[5:8]) +
     sum(n_E3124_I3124[5:8])  else 
       if (i==4)   sum(n_E3_I3[9:13]) +
     sum(n_E31_I31[9:13]) + sum(n_E32_I32[9:13]) +
     sum(n_E34_I34[9:13]) + sum(n_E312_I312[9:13]) +
     sum(n_E314_I314[9:13]) + sum(n_E324_I324[9:13]) +
     sum(n_E3124_I3124[9:13]) else
       sum(n_E3_I3[14]) + sum(n_E31_I31[14]) +
     sum(n_E32_I32[14]) + sum(n_E34_I34[14]) +
     sum(n_E312_I312[14]) + sum(n_E314_I314[14]) +
     sum(n_E324_I324[14]) + sum(n_E3124_I3124[14]) ) 


update(inc_day_gii[])<-
  (if(i==1)    sum(n_E4_I4) +  sum(n_E41_I41[1]) +
     sum(n_E42_I42[1]) + sum(n_E43_I43[1]) +
     sum(n_E412_I412[1]) + sum(n_E413_I413[1]) +
     sum(n_E423_I423[1]) + sum(n_E4123_I4123[1])  else
       if(i==2)    
         sum(n_E4_I4[2:4]) + sum(n_E41_I41[2:4]) +
     sum(n_E42_I42[2:4]) + sum(n_E43_I43[2:4]) +
     sum(n_E412_I412[2:4]) + sum(n_E413_I413[2:4]) +
     sum(n_E423_I423[2:4]) + sum(n_E4123_I4123[2:4])   else 
       if(i==3)    
         sum(n_E4_I4[5:8]) + sum(n_E41_I41[5:8]) +
     sum(n_E42_I42[5:8]) + sum(n_E43_I43[5:8]) +
     sum(n_E412_I412[5:8]) + sum(n_E413_I413[5:8]) +
     sum(n_E423_I423[5:8]) + sum(n_E4123_I4123[5:8])  else 
       if (i==4)  
         sum(n_E4_I4[9:13]) + sum(n_E41_I41[9:13]) +
     sum(n_E42_I42[9:13]) + sum(n_E43_I43[9:13]) +
     sum(n_E412_I412[9:13]) + sum(n_E413_I413[9:13]) +
     sum(n_E423_I423[9:13]) + sum(n_E4123_I4123[9:13])   else
       sum(n_E4_I4[14]) +  sum(n_E41_I41[14]) +
     sum(n_E42_I42[14]) + sum(n_E43_I43[14]) +
     sum(n_E412_I412[14]) + sum(n_E413_I413[14]) +
     sum(n_E423_I423[14]) + sum(n_E4123_I4123[14]) ) 



dim(pop_forinc)<-5

pop_forinc[1]<-N_byage[1]
pop_forinc[2]<-sum(N_byage[2:4])
pop_forinc[3]<-sum(N_byage[5:8])
pop_forinc[4]<-sum(N_byage[9:13])
pop_forinc[5]<-N_byage[14]




update(inc_day_all[])<-1000*(inc_day_gi3[i]+inc_day_gi[i]+inc_day_gii4[i]+inc_day_gii[i])/pop_forinc[i]

update(new_cases)<- sum(inc_day_gi3)+sum(inc_day_gi)+sum(inc_day_gii4)+sum(inc_day_gii)

update(new_cases_week)<- if (step %% steps_per_week == 0)
  sum(inc_day_gi3)+sum(inc_day_gi)+sum(inc_day_gii4)+sum(inc_day_gii) else 
    new_cases_week +  sum(inc_day_gi3)+sum(inc_day_gi)+sum(inc_day_gii4)+sum(inc_day_gii)


update(new_cases_week_gi)<- if (step %% steps_per_week == 0)
  sum(inc_day_gi) else 
    new_cases_week_gi +  sum(inc_day_gi)

update(new_cases_week_gi3)<- if (step %% steps_per_week == 0)
  sum(inc_day_gi3) else 
    new_cases_week_gi3 +  sum(inc_day_gi3)

update(new_cases_week_gii)<- if (step %% steps_per_week == 0)
  sum(inc_day_gii) else 
    new_cases_week_gii +  sum(inc_day_gii)

update(new_cases_week_gii4)<- if (step %% steps_per_week == 0)
  sum(inc_day_gii4) else 
    new_cases_week_gii4 +  sum(inc_day_gii4)



update(seasonality) <- (1 + w1_1*cos((2*pi*time)/364 + (w2/12)*pi))


## Incidence by Year and strain
###############


update(inc_year_gi3[]) <- if (step %% steps_per_year==0)
  (if(i==1) 
    (n_E1_I1[1]) + 
     (n_E12_I12[1]) +
     (n_E13_I13[1]) + 
     (n_E14_I14[1]) + 
     (n_E123_I123[1]) + 
     (n_E124_I124[1]) + 
     (n_E134_I134[1]) + 
     (n_E1234_I1234[1])  else 
       if(i==2)  
         sum(n_E1_I1[2:4]) + 
     sum(n_E12_I12[2:4]) +
     sum(n_E13_I13[2:4]) + 
     sum(n_E14_I14[2:4]) + 
     sum(n_E123_I123[2:4]) + 
     sum(n_E124_I124[2:4]) +
     sum(n_E134_I134[2:4]) + 
     sum(n_E1234_I1234[2:4])  else 
       if(i==3)  
         sum(n_E1_I1[5:8]) + 
     sum(n_E12_I12[5:8]) + 
     sum(n_E13_I13[5:8]) +
     sum(n_E14_I14[5:8]) + 
     sum(n_E123_I123[5:8]) + 
     sum(n_E124_I124[5:8]) +
     sum(n_E134_I134[5:8]) + 
     sum(n_E1234_I1234[5:8]) else 
       if (i==4) 
         sum(n_E1_I1[9:13]) + 
     sum(n_E12_I12[9:13]) + 
     sum(n_E13_I13[9:13]) +
     sum(n_E14_I14[9:13]) + 
     sum(n_E123_I123[9:13]) + 
     sum(n_E124_I124[9:13]) +
     sum(n_E134_I134[9:13]) + 
     sum(n_E1234_I1234[9:13])  
   else
     (n_E1_I1[14]) + 
     (n_E12_I12[14]) + 
     (n_E13_I13[14]) +
     (n_E14_I14[14]) + 
     (n_E123_I123[14]) + 
     (n_E124_I124[14]) +
     (n_E134_I134[14]) + 
     (n_E1234_I1234[14])
  ) else ( 
    if(i==1) 
      inc_year_gi3[1] +
      sum(n_E1_I1[1]) + 
      sum(n_E12_I12[1]) +
      sum(n_E13_I13[1]) + 
      sum(n_E14_I14[1]) + 
      sum(n_E123_I123[1]) + 
      sum(n_E124_I124[1]) + 
      sum(n_E134_I134[1]) + 
      sum(n_E1234_I1234[1])  else 
        if(i==2) 
          inc_year_gi3[2] +  
      sum(n_E1_I1[2:4]) + 
      sum(n_E12_I12[2:4]) +
      sum(n_E13_I13[2:4]) + 
      sum(n_E14_I14[2:4]) + 
      sum(n_E123_I123[2:4]) +
      sum(n_E124_I124[2:4]) +
      sum(n_E134_I134[2:4]) + 
      sum(n_E1234_I1234[2:4])  else 
        if(i==3) 
          inc_year_gi3[3] +  
      sum(n_E1_I1[5:8]) + 
      sum(n_E12_I12[5:8]) + 
      sum(n_E13_I13[5:8]) +
      sum(n_E14_I14[5:8]) + 
      sum(n_E123_I123[5:8]) + 
      sum(n_E124_I124[5:8]) +
      sum(n_E134_I134[5:8]) + 
      sum(n_E1234_I1234[5:8]) else 
        if (i==4) 
          inc_year_gi3[4] + 
      sum(n_E1_I1[9:13]) + 
      sum(n_E12_I12[9:13]) + 
      sum(n_E13_I13[9:13]) +
      sum(n_E14_I14[9:13]) + 
      sum(n_E123_I123[9:13]) + 
      sum(n_E124_I124[9:13]) +
      sum(n_E134_I134[9:13]) + 
      sum(n_E1234_I1234[9:13])  
    else
      inc_year_gi3[5] +
      sum(n_E1_I1[14]) + 
      sum(n_E12_I12[14]) +
      sum(n_E13_I13[14]) +
      sum(n_E14_I14[14]) + 
      sum(n_E123_I123[14]) + 
      sum(n_E124_I124[14]) +
      sum(n_E134_I134[14]) + 
      sum(n_E1234_I1234[14])
  )




update(inc_year_gi[]) <-  if (step %% steps_per_year==0)
  ( if(i==1) 
    sum(n_E2_I2[1]) +
      sum(n_E21_I21[1]) +  
      sum(n_E23_I23[1]) +  
      sum(n_E24_I24[1]) +
      sum(n_E213_I213[1]) + 
      sum(n_E214_I214[1]) + 
      sum(n_E234_I234[1]) +
      sum(n_E2134_I2134[1])  else
        if(i==2)  
          sum(n_E2_I2[2:4]) +
      sum(n_E21_I21[2:4]) + 
      sum(n_E23_I23[2:4]) +
      sum(n_E24_I24[2:4]) + 
      sum(n_E213_I213[2:4]) +
      sum(n_E214_I214[2:4]) + 
      sum(n_E234_I234[2:4]) +
      sum(n_E2134_I2134[2:4])  else 
        if(i==3)  
          sum(n_E2_I2[5:8]) +
      sum(n_E21_I21[5:8]) + 
      sum(n_E23_I23[5:8]) +
      sum(n_E24_I24[5:8]) + 
      sum(n_E213_I213[5:8]) +
      sum(n_E214_I214[5:8]) + 
      sum(n_E234_I234[5:8]) +
      sum(n_E2134_I2134[5:8])  else 
        if (i==4) 
          sum(n_E2_I2[9:13]) +
      sum(n_E21_I21[9:13]) + 
      sum(n_E23_I23[9:13]) +
      sum(n_E24_I24[9:13]) + 
      sum(n_E213_I213[9:13]) +
      sum(n_E214_I214[9:13]) + 
      sum(n_E234_I234[9:13]) +
      sum(n_E2134_I2134[9:13])  else
        sum(n_E2_I2[14]) + 
      sum(n_E21_I21[14]) +
      sum(n_E23_I23[14]) + 
      sum(n_E24_I24[14]) +
      sum(n_E213_I213[14]) + 
      sum(n_E214_I214[14]) +
      sum(n_E234_I234[14]) + 
      sum(n_E2134_I2134[14])) else
        (if(i==1) 
          inc_year_gi[1] +
           sum(n_E2_I2[1]) +
           sum(n_E21_I21[1]) +  
           sum(n_E23_I23[1]) +  
           sum(n_E24_I24[1]) +
           sum(n_E213_I213[1]) + 
           sum(n_E214_I214[1]) + 
           sum(n_E234_I234[1]) +
           sum(n_E2134_I2134[1])  else
             if(i==2)
               inc_year_gi[2] +
           sum(n_E2_I2[2:4]) +
           sum(n_E21_I21[2:4]) + 
           sum(n_E23_I23[2:4]) +
           sum(n_E24_I24[2:4]) + 
           sum(n_E213_I213[2:4]) +
           sum(n_E214_I214[2:4]) + 
           sum(n_E234_I234[2:4]) +
           sum(n_E2134_I2134[2:4])  else 
             if(i==3)  
               inc_year_gi[3] +
           sum(n_E2_I2[5:8]) +
           sum(n_E21_I21[5:8]) +
           sum(n_E23_I23[5:8]) +
           sum(n_E24_I24[5:8]) + 
           sum(n_E213_I213[5:8]) +
           sum(n_E214_I214[5:8]) + 
           sum(n_E234_I234[5:8]) +
           sum(n_E2134_I2134[5:8])  else 
             if (i==4) 
               inc_year_gi[4] +
           sum(n_E2_I2[9:13]) +
           sum(n_E21_I21[9:13]) + 
           sum(n_E23_I23[9:13]) +
           sum(n_E24_I24[9:13]) + 
           sum(n_E213_I213[9:13]) +
           sum(n_E214_I214[9:13]) + 
           sum(n_E234_I234[9:13]) +
           sum(n_E2134_I2134[9:13])  
         else
           inc_year_gi[5] +
           sum(n_E2_I2[14]) + 
           sum(n_E21_I21[14]) +
           sum(n_E23_I23[14]) + 
           sum(n_E24_I24[14]) +
           sum(n_E213_I213[14]) + 
           sum(n_E214_I214[14]) +
           sum(n_E234_I234[14]) + 
           sum(n_E2134_I2134[14]) )


update(inc_year_gii4[]) <- if (step %% steps_per_year==0)
  (if(i==1)   
    sum(n_E3_I3[1]) +
     sum(n_E31_I31[1]) + 
     sum(n_E32_I32[1]) +
     sum(n_E34_I34[1]) + 
     sum(n_E312_I312[1]) +
     sum(n_E314_I314[1]) + 
     sum(n_E324_I324[1]) +
     sum(n_E3124_I3124[1]) else
       if(i==2)    
         sum(n_E3_I3[2:4]) +
     sum(n_E31_I31[2:4]) + 
     sum(n_E32_I32[2:4]) +
     sum(n_E34_I34[2:4]) + 
     sum(n_E312_I312[2:4]) +
     sum(n_E314_I314[2:4]) + 
     sum(n_E324_I324[2:4]) +
     sum(n_E3124_I3124[2:4])  else 
       if(i==3)    
         sum(n_E3_I3[5:8]) +
     sum(n_E31_I31[5:8]) + 
     sum(n_E32_I32[5:8]) +
     sum(n_E34_I34[5:8]) + 
     sum(n_E312_I312[5:8]) +
     sum(n_E314_I314[5:8]) + 
     sum(n_E324_I324[5:8]) +
     sum(n_E3124_I3124[5:8])  else 
       if (i==4)   
         sum(n_E3_I3[9:13]) +
     sum(n_E31_I31[9:13]) + 
     sum(n_E32_I32[9:13]) +
     sum(n_E34_I34[9:13]) + 
     sum(n_E312_I312[9:13]) +
     sum(n_E314_I314[9:13]) + 
     sum(n_E324_I324[9:13]) +
     sum(n_E3124_I3124[9:13]) else
       sum(n_E3_I3[14]) + 
     sum(n_E31_I31[14]) +
     sum(n_E32_I32[14]) + 
     sum(n_E34_I34[14]) +
     sum(n_E312_I312[14]) + 
     sum(n_E314_I314[14]) +
     sum(n_E324_I324[14]) + 
     sum(n_E3124_I3124[14]) ) else
       (if(i==1)   
         inc_year_gii4[1] +
          sum(n_E3_I3[1]) +
          sum(n_E31_I31[1]) + 
          sum(n_E32_I32[1]) +
          sum(n_E34_I34[1]) + 
          sum(n_E312_I312[1]) +
          sum(n_E314_I314[1]) + 
          sum(n_E324_I324[1]) +
          sum(n_E3124_I3124[1]) else
            if(i==2) 
              inc_year_gii4[2] +
          sum(n_E3_I3[2:4]) +
          sum(n_E31_I31[2:4]) + 
          sum(n_E32_I32[2:4]) +
          sum(n_E34_I34[2:4]) + 
          sum(n_E312_I312[2:4]) +
          sum(n_E314_I314[2:4]) + 
          sum(n_E324_I324[2:4]) +
          sum(n_E3124_I3124[2:4])  else 
            if(i==3)    
              inc_year_gii4[3] +
          sum(n_E3_I3[5:8]) +
          sum(n_E31_I31[5:8]) + 
          sum(n_E32_I32[5:8]) +
          sum(n_E34_I34[5:8]) + 
          sum(n_E312_I312[5:8]) +
          sum(n_E314_I314[5:8]) + 
          sum(n_E324_I324[5:8]) +
          sum(n_E3124_I3124[5:8])  else 
            if (i==4)   
              inc_year_gii4[4] +
          sum(n_E3_I3[9:13]) +
          sum(n_E31_I31[9:13]) + 
          sum(n_E32_I32[9:13]) +
          sum(n_E34_I34[9:13]) + 
          sum(n_E312_I312[9:13]) +
          sum(n_E314_I314[9:13]) + 
          sum(n_E324_I324[9:13]) +
          sum(n_E3124_I3124[9:13]) else
            inc_year_gii4[5] +
          sum(n_E3_I3[14]) + 
          sum(n_E31_I31[14]) +
          sum(n_E32_I32[14]) + 
          sum(n_E34_I34[14]) +
          sum(n_E312_I312[14]) + 
          sum(n_E314_I314[14]) +
          sum(n_E324_I324[14]) + 
          sum(n_E3124_I3124[14]) )

update(inc_year_gii[]) <-  if (step %% steps_per_year==0)
  (if(i==1)    
    sum(n_E4_I4[1]) +  
     sum(n_E41_I41[1]) +
     sum(n_E42_I42[1]) + 
     sum(n_E43_I43[1]) +
     sum(n_E412_I412[1]) + 
     sum(n_E413_I413[1]) +
     sum(n_E423_I423[1]) + 
     sum(n_E4123_I4123[1])  
   else
     if(i==2)    
       sum(n_E4_I4[2:4]) + 
     sum(n_E41_I41[2:4]) +
     sum(n_E42_I42[2:4]) + 
     sum(n_E43_I43[2:4]) +
     sum(n_E412_I412[2:4]) + 
     sum(n_E413_I413[2:4]) +
     sum(n_E423_I423[2:4]) + 
     sum(n_E4123_I4123[2:4])   
   else 
     if(i==3)    
       sum(n_E4_I4[5:8]) + 
     sum(n_E41_I41[5:8]) +
     sum(n_E42_I42[5:8]) + 
     sum(n_E43_I43[5:8]) +
     sum(n_E412_I412[5:8]) + 
     sum(n_E413_I413[5:8]) +
     sum(n_E423_I423[5:8]) + 
     sum(n_E4123_I4123[5:8])  
   else 
     if (i==4)  
       sum(n_E4_I4[9:13]) + 
     sum(n_E41_I41[9:13]) +
     sum(n_E42_I42[9:13]) + 
     sum(n_E43_I43[9:13]) +
     sum(n_E412_I412[9:13]) + 
     sum(n_E413_I413[9:13]) +
     sum(n_E423_I423[9:13]) + 
     sum(n_E4123_I4123[9:13])   
   else
     sum(n_E4_I4[14]) +  
     sum(n_E41_I41[14]) +
     sum(n_E42_I42[14]) + 
     sum(n_E43_I43[14]) +
     sum(n_E412_I412[14]) + 
     sum(n_E413_I413[14]) +
     sum(n_E423_I423[14]) + 
     sum(n_E4123_I4123[14]) ) else
       ( if(i==1)    
         inc_year_gii[1] +
           sum(n_E4_I4[1]) +  
           sum(n_E41_I41[1]) +
           sum(n_E42_I42[1]) + 
           sum(n_E43_I43[1]) +
           sum(n_E412_I412[1]) + 
           sum(n_E413_I413[1]) +
           sum(n_E423_I423[1]) + 
           sum(n_E4123_I4123[1])  
         else
           if(i==2)
             inc_year_gii[2] +
           sum(n_E4_I4[2:4]) + 
           sum(n_E41_I41[2:4]) +
           sum(n_E42_I42[2:4]) + 
           sum(n_E43_I43[2:4]) +
           sum(n_E412_I412[2:4]) + 
           sum(n_E413_I413[2:4]) +
           sum(n_E423_I423[2:4]) + 
           sum(n_E4123_I4123[2:4])  
         else 
           if(i==3)    
             inc_year_gii[3] +
           sum(n_E4_I4[5:8]) + 
           sum(n_E41_I41[5:8]) +
           sum(n_E42_I42[5:8]) + 
           sum(n_E43_I43[5:8]) +
           sum(n_E412_I412[5:8]) + 
           sum(n_E413_I413[5:8]) +
           sum(n_E423_I423[5:8]) + 
           sum(n_E4123_I4123[5:8])  
         else 
           if (i==4)  
             inc_year_gii[4] +
           sum(n_E4_I4[9:13]) + 
           sum(n_E41_I41[9:13]) +
           sum(n_E42_I42[9:13]) + 
           sum(n_E43_I43[9:13]) +
           sum(n_E412_I412[9:13]) + 
           sum(n_E413_I413[9:13]) +
           sum(n_E423_I423[9:13]) + 
           sum(n_E4123_I4123[9:13])   
         else
           inc_year_gii[5] +
           sum(n_E4_I4[14]) +  
           sum(n_E41_I41[14]) +
           sum(n_E42_I42[14]) + 
           sum(n_E43_I43[14]) +
           sum(n_E412_I412[14]) + 
           sum(n_E413_I413[14]) +
           sum(n_E423_I423[14]) + 
           sum(n_E4123_I4123[14]))


## compute at risk population ( match IDD2 data)
update(pop_by4age[]) <-
  (if(i==1)  sum(N_byage[1]) else
    if(i==2)  sum(N_byage[2:4])  else 
      if(i==3)  sum(N_byage[5:8])  else 
        if (i==4)  sum(N_byage[9:13])  else
          sum(N_byage[14]) ) 

## compute prevalence og GII4
update(seroprev_num[]) <- if (step %% steps_per_year==0)
  I3[i]+ 
  A3[i]+
  I13[i]+
  A13[i]+
  I23[i]+
  A23[i]+
  I31[i]+ 
  A31[i]+
  I34[i]+
  A34[i]+
  I43[i]+
  A43[i]+
  I123[i]+
  A123[i]+
  I134[i]+
  A134[i]+
  I213[i]+
  A213[i]+
  I234[i]+
  A234[i]+
  I312[i]+ 
  A312[i]+
  I314[i]+
  A314[i]+
  I324[i]+ 
  A324[i]+
  I413[i]+
  A413[i]+
  I423[i]+
  A423[i]+
  I1234[i]+
  A1234[i]+
  I2134[i]+ 
  A2134[i]+
  I3124[i]+
  A3124[i]+
  I4123[i]+
  A4123[i]+
  R123[i] +
  R1234[i] +
  R13[i] +
  R134[i] +
  R23[i] +
  R234[i] +
  R3[i] +
  R34[i]  else
    seroprev_num[i] + 
  I3[i]+ 
  A3[i]+
  I13[i]+
  A13[i]+
  I23[i]+
  A23[i]+
  I31[i]+ 
  A31[i]+
  I34[i]+
  A34[i]+
  I43[i]+
  A43[i]+
  I123[i]+
  A123[i]+
  I134[i]+
  A134[i]+
  I213[i]+
  A213[i]+
  I234[i]+
  A234[i]+
  I312[i]+ 
  A312[i]+
  I314[i]+
  A314[i]+
  I324[i]+ 
  A324[i]+
  I413[i]+
  A413[i]+
  I423[i]+
  A423[i]+
  I1234[i]+
  A1234[i]+
  I2134[i]+ 
  A2134[i]+
  I3124[i]+
  A3124[i]+
  I4123[i]+
  A4123[i]+
  R123[i] +
  R1234[i] +
  R13[i] +
  R134[i] +
  R23[i] +
  R234[i] +
  R3[i] +
  R34[i]


update(seroprev_den[]) <- if (step %% steps_per_year==0)
  N_byage[i] else
    seroprev_den[i] + N_byage[i] 

update(seroprev[]) <- prev_byage[i]/N_byage[i]


update(sus_gi3) <- 
  sum(S)+
  sum(R2)+
  sum(R3)+
  sum(R4)+
  sum(R23)+
  sum(R24)+
  sum(R34)+
  sum(R234)

update(sus_gi) <- 
  sum(S)+
  sum(R1)+
  sum(R3)+
  sum(R4)+
  sum(R13)+
  sum(R14)+
  sum(R34)+
  sum(R134)

update(sus_gii4) <- 
  sum(S)+
  sum(R1)+
  sum(R2)+
  sum(R4)+
  sum(R12)+
  sum(R14)+
  sum(R24)+
  sum(R124)

update(sus_gii) <- 
  sum(S)+
  sum(R1)+
  sum(R2)+
  sum(R3)+
  sum(R12)+
  sum(R13)+
  sum(R23)+
  sum(R123)


update(covid_timeline)<-comix_switcher


#### Individual probabilities of transition:
#######
p_mu[] <- 1 - exp(-mu[i] * dt) # mortality
p_aging[]<- 1 - exp(-aging_vec[i] * dt)
p_MS[] <- 1 - exp(-(1/maternalAB) * dt)  # M to S
p_SE1[] <-  1 - exp(-lambda_1[i] * dt) # S to E
p_SE2[] <- 1 - exp(-lambda_2[i] * dt) # S to E
p_SE3[] <- 1 - exp(-lambda_3[i] * dt) # S to E
p_SE4[] <- 1 - exp(-lambda_4[i] * dt) # S to E
p_R1A1[] <- 1 - exp(-lambda_1[i] * (1-alpha) * dt)
p_R2A2[] <- 1 - exp(-lambda_2[i] * (1-alpha) * dt)
p_R3A3[] <- 1 - exp(-lambda_3[i]  *(1-alpha) * dt)
p_R4A4[] <- 1 - exp(-lambda_4[i] * (1-alpha) * dt)
p_R1_to_E2_gi[] <- 1 - exp(-lambda_2[i] * (1-crossp_GI) * dt)
p_R2_to_E1_gi[] <- 1 - exp(-lambda_1[i] * (1-crossp_GI) * dt)
p_R3_to_E4_gii[] <- 1 - exp(-lambda_4[i] * (1-crossp_GII) * dt)
p_R4_to_E3_gii[] <- 1 - exp(-lambda_3[i] * (1-crossp_GII) * dt)
p_EI   <- 1 - exp(-(1/epsilon) * dt) # E to I
p_IA_5    <- 1 - exp(-(1/theta_5) * dt) # I to A
p_IA_5p   <- 1 - exp(-(1/theta_5p) * dt) # I to A
p_AR   <- 1 - exp(-(1/sigma) * dt) # A to R
p_RS   <- 1 - exp(- (1/(imm_yr*365)) * dt) # R to S
p_RS2   <- 1 - exp(- (1/(imm_yr*imm_fac1*365)) * dt) # R to S (1 infection)
p_RS3   <- 1 - exp(- (1/(imm_yr*imm_fac1*imm_fac2*365)) * dt) # R to S (2 infection)
p_RS4   <- 1 - exp(- (1/(imm_yr*imm_fac1*imm_fac2*imm_fac3*365)) * dt) # R to S (3 infection)
p_vacc[]<- 1 - exp(- (1/365) * vac_eff * (
  vac_camp_switch * vac_camp_cov[i]  + vac_sche_cov[i]) * dt) # R to S

vac_camp_switch <- if (as.integer(step) <= 365) 1L else 0L

N_byage[]<-( M[i] +  
               G[i] +
               S[i] + 
               E1[i]+
               I1[i]+
               A1[i]+
               R1[i]+
               E2[i]+
               I2[i]+ 
               A2[i]+
               R2[i]+
               E3[i]+
               I3[i]+ 
               A3[i]+
               R3[i]+
               E4[i]+ 
               I4[i]+ 
               A4[i]+
               R4[i]+
               E12[i]+
               I12[i]+ 
               A12[i]+ 
               E13[i]+
               I13[i]+
               A13[i]+
               E14[i]+
               I14[i]+
               A14[i]+
               E21[i]+ 
               I21[i]+ 
               A21[i]+
               E23[i]+
               I23[i]+
               A23[i]+
               E24[i]+
               I24[i]+
               A24[i]+
               E31[i]+
               I31[i]+ 
               A31[i]+
               E32[i]+ 
               I32[i]+
               A32[i]+
               E34[i]+
               I34[i]+
               A34[i]+
               E41[i]+ 
               I41[i]+
               A41[i]+
               E42[i]+ 
               I42[i]+
               A42[i]+
               E43[i]+
               I43[i]+
               A43[i]+
               R12[i]+
               R13[i]+
               R14[i]+
               R23[i]+
               R24[i]+
               R34[i]+
               E123[i]+
               I123[i]+
               A123[i]+
               E124[i]+
               I124[i]+
               A124[i]+
               E134[i]+
               I134[i]+
               A134[i]+
               E213[i]+
               I213[i]+
               A213[i]+
               E214[i]+
               I214[i]+
               A214[i]+
               E234[i]+ 
               I234[i]+
               A234[i]+
               E312[i]+ 
               I312[i]+ 
               A312[i]+
               E314[i]+ 
               I314[i]+
               A314[i]+
               E324[i]+
               I324[i]+ 
               A324[i]+
               E412[i]+ 
               I412[i]+
               A412[i]+
               E413[i]+
               I413[i]+
               A413[i]+
               E423[i]+
               I423[i]+
               A423[i]+
               R123[i]+
               R124[i]+
               R134[i]+
               R234[i]+
               E1234[i]+
               I1234[i]+
               A1234[i]+
               E2134[i]+
               I2134[i]+ 
               A2134[i]+
               E3124[i]+ 
               I3124[i]+
               A3124[i]+
               E4123[i]+
               I4123[i]+
               A4123[i]+
               R1234[i]) 


prev_byage[]<-(M[i] +
                 A3[i]+
                 A13[i]+
                 A23[i]+
                 A31[i]+
                 A34[i]+
                 A43[i]+
                 A123[i]+
                 A134[i]+
                 A213[i]+
                 A234[i]+
                 A312[i]+
                 A314[i]+
                 A324[i]+
                 A413[i]+
                 A423[i]+
                 A1234[i]+
                 A2134[i]+
                 A3124[i]+
                 A4123[i]+
                 R1[i] +
                 R12[i] +
                 R123[i] +
                 R1234[i] +
                 R124[i] +
                 R13[i] +
                 R134[i] +
                 R14[i] +
                 R2[i] +
                 R23[i] +
                 R234[i] +
                 R24[i] +
                 R3[i] +
                 R34[i] +
                 R4[i])

dim(N_byage)<-N_age
dim(prev_byage)<-N_age

N <- sum(N_byage)

prev   <- sum(prev_byage)/N  

p_birth[1]   <- (1-p_nonsecretor)*prev
p_birth[2]   <- p_nonsecretor
p_birth[3]   <- (1-p_nonsecretor)*(1-prev) 


## Force of infection
###########################
m[, ] <- user() # age-structured contact matrix
m_holi[,]<-user()
cmx_1[, ] <- user() 
cmx_2[, ] <- user() 
cmx_3[, ] <- user() 
cmx_4[, ] <- user() 
cmx_5[, ] <- user() 
cmx_6[, ] <- user() 
cmx_7[, ] <- user() 
cmx_8[, ] <- user() 
cmx_9[, ] <- user() 


school_switcher <- if (as.integer(step) >= n_school_steps)
  school_step[n_school_steps]*
  (1- (as.integer(step)-n_school_steps)) else 
    school_step[as.integer(step) + 1] 

comix_switcher <- if (as.integer(step) >= n_covid_steps)
  covid_step[n_covid_steps]*
  (1- (as.integer(step)-n_covid_steps)) else 
    covid_step[as.integer(step) + 1] 



contact_matrix[,]<- if (comix_switcher >0)( 
  if (comix_switcher==1)
    cmx_1[i,j] else 
      if (comix_switcher==2)
        cmx_2[i,j] else 
          if (comix_switcher==3)
            cmx_3[i,j] else 
              if (comix_switcher==4)
                cmx_4[i,j] else 
                  if (comix_switcher==5)
                    cmx_5[i,j] else 
                      if (comix_switcher==6)
                        cmx_6[i,j] else 
                          if (comix_switcher==7)
                            cmx_7[i,j] else 
                              if (comix_switcher==8)
                                cmx_8[i,j] else 
                                  cmx_9[i,j] 
) else (
  school_switcher*m[i, j]  + (1-school_switcher)*m_holi[i, j])







# infectives GI3
c1_ij[, ] <- if (j>4) (I1[j]+I12[j] + 
                         I13[j]+I14[j] +
                         I123[j]+I124[j] + 
                         I134[j]+I1234[j] + 
                         rr_inf_asymp * (A1[j]+A12[j] + 
                                           A13[j]+A14[j] +
                                           A123[j]+A124[j] + 
                                           A134[j]+A1234[j] +
                                           E1[j]+E12[j] + 
                                           E13[j]+E14[j] +
                                           E123[j]+E124[j] + 
                                           E134[j]+E1234[j]) ) * contact_matrix[i,j] * aduRR else
                                             (I1[j]+I12[j] + 
                                                I13[j]+I14[j] +
                                                I123[j]+I124[j] + 
                                                I134[j]+I1234[j] + 
                                                rr_inf_asymp * (A1[j]+A12[j] + 
                                                                  A13[j]+A14[j] +
                                                                  A123[j]+A124[j] + 
                                                                  A134[j]+A1234[j] +
                                                                  E1[j]+E12[j] + 
                                                                  E13[j]+E14[j] +
                                                                  E123[j]+E124[j] + 
                                                                  E134[j]+E1234[j]) ) * contact_matrix[i,j]



# infectives GI
c2_ij[, ] <- if (j>4) (I2[j]+I21[j]+
                         I23[j]+I24[j]+
                         I213[j]+I214[j]+
                         I234[j]+I2134[j]+
                         rr_inf_asymp * (A2[j]+A21[j]+
                                           A23[j]+A24[j]+
                                           A213[j]+A214[j]+
                                           A234[j]+A2134[j]+
                                           E2[j]+E21[j]+
                                           E23[j]+E24[j]+
                                           E213[j]+E214[j]+
                                           E234[j]+E2134[j])) * contact_matrix[i,j] * aduRR else
                                             (I2[j]+I21[j]+
                                                I23[j]+I24[j]+
                                                I213[j]+I214[j]+
                                                I234[j]+I2134[j]+
                                                rr_inf_asymp * (A2[j]+A21[j]+
                                                                  A23[j]+A24[j]+
                                                                  A213[j]+A214[j]+
                                                                  A234[j]+A2134[j]+
                                                                  E2[j]+E21[j]+
                                                                  E23[j]+E24[j]+
                                                                  E213[j]+E214[j]+
                                                                  E234[j]+E2134[j])) * contact_matrix[i,j]



# infectives GII4
c3_ij[, ] <- if (j>4 )(I3[j]+I31[j]+
                         I32[j]+I34[j]+
                         I312[j]+I314[j]+
                         I324[j]+I3124[j]+ 
                         rr_inf_asymp *(A3[j]+A31[j]+
                                          A32[j]+A34[j]+
                                          A312[j]+A314[j]+
                                          A324[j]+A3124[j]+
                                          E3[j]+E31[j]+
                                          E32[j]+E34[j]+
                                          E312[j]+E314[j]+
                                          E324[j]+E3124[j]) ) * contact_matrix[i,j] * aduRR else
                                            (I3[j]+I31[j]+
                                               I32[j]+I34[j]+
                                               I312[j]+I314[j]+
                                               I324[j]+I3124[j]+ 
                                               rr_inf_asymp *(A3[j]+A31[j]+
                                                                A32[j]+A34[j]+
                                                                A312[j]+A314[j]+
                                                                A324[j]+A3124[j]+
                                                                E3[j]+E31[j]+
                                                                E32[j]+E34[j]+
                                                                E312[j]+E314[j]+
                                                                E324[j]+E3124[j]) ) * contact_matrix[i,j]



# infectives GII
c4_ij[, ] <- if (j>4) (I4[j]+I41[j]+
                         I42[j]+I43[j]+
                         I412[j]+I413[j]+
                         I423[j]+I4123[j] +
                         rr_inf_asymp *(E4[j]+E41[j]+
                                          E42[j]+E43[j]+
                                          E412[j]+E413[j]+
                                          E423[j]+E4123[j]+
                                          A4[j]+A41[j]+
                                          A42[j]+A43[j]+
                                          A412[j]+A413[j]+
                                          A423[j]+A4123[j]) ) * contact_matrix[i,j] * aduRR else
                                            (I4[j]+I41[j]+
                                               I42[j]+I43[j]+
                                               I412[j]+I413[j]+
                                               I423[j]+I4123[j] +
                                               rr_inf_asymp *(E4[j]+E41[j]+
                                                                E42[j]+E43[j]+
                                                                E412[j]+E413[j]+
                                                                E423[j]+E4123[j]+
                                                                A4[j]+A41[j]+
                                                                A42[j]+A43[j]+
                                                                A412[j]+A413[j]+
                                                                A423[j]+A4123[j]) ) * contact_matrix[i,j]





# if (school[as.integer(step)]==1) m[i, j] else m_holi[i, j]

beta_1_t <- beta_1 *(1 + w1_1*cos((2*pi*time)/364 + (w2/12)*pi))
beta_2_t <- beta_2 *(1 + w1_1*cos((2*pi*time)/364 + (w2/12)*pi))
beta_3_t <- beta_3 *(1 + w1_1*cos((2*pi*time)/364 + (w2/12)*pi))
beta_4_t <- beta_4 *(1 + w1_1*cos((2*pi*time)/364 + (w2/12)*pi))


lambda_1[] <-  beta_1_t   * sum(c1_ij[i , ])
lambda_2[] <-  beta_2_t   * sum(c2_ij[i , ])
lambda_3[] <-  beta_3_t   * sum(c3_ij[i , ])
lambda_4[] <-  beta_4_t   * sum(c4_ij[i , ])



########### Aging numbers ::::::::::::::::::::::::
####################################################
#Age out
n_ageoM[] <- round(M[i] *  p_aging[i])
n_ageoG[] <- round(G[i] *  p_aging[i])
n_ageoS[] <- round(S[i] *  p_aging[i])
n_ageoE1[] <- round(E1[i] *  p_aging[i])
n_ageoI1[] <- round(I1[i] *  p_aging[i])
n_ageoA1[] <- round(A1[i] *  p_aging[i])
n_ageoR1[] <- round(R1[i] *  p_aging[i])
n_ageoE2[] <- round(E2[i] *  p_aging[i])
n_ageoI2[] <- round(I2[i] *  p_aging[i])
n_ageoA2[] <- round(A2[i] *  p_aging[i])
n_ageoR2[] <- round(R2[i] *  p_aging[i])
n_ageoE3[] <- round(E3[i] *  p_aging[i])
n_ageoI3[] <- round(I3[i] *  p_aging[i])
n_ageoA3[] <- round(A3[i] *  p_aging[i])
n_ageoR3[] <- round(R3[i] *  p_aging[i])
n_ageoE4[] <- round(E4[i] *  p_aging[i])
n_ageoI4[] <- round(I4[i] *  p_aging[i])
n_ageoA4[] <- round(A4[i] *  p_aging[i])
n_ageoR4[] <- round(R4[i] *  p_aging[i])
n_ageoE12[] <- round(E12[i] *  p_aging[i])
n_ageoI12[] <- round(I12[i] *  p_aging[i])
n_ageoA12[] <- round(A12[i] *  p_aging[i])
n_ageoE13[] <- round(E13[i] *  p_aging[i])
n_ageoI13[] <- round(I13[i] *  p_aging[i])
n_ageoA13[] <- round(A13[i] *  p_aging[i])
n_ageoE14[] <- round(E14[i] *  p_aging[i])
n_ageoI14[] <- round(I14[i] *  p_aging[i])
n_ageoA14[] <- round(A14[i] *  p_aging[i])
n_ageoE21[] <- round(E21[i] *  p_aging[i])
n_ageoI21[] <- round(I21[i] *  p_aging[i])
n_ageoA21[] <- round(A21[i] *  p_aging[i])
n_ageoE23[] <- round(E23[i] *  p_aging[i])
n_ageoI23[] <- round(I23[i] *  p_aging[i])
n_ageoA23[] <- round(A23[i] *  p_aging[i])
n_ageoE24[] <- round(E24[i] *  p_aging[i])
n_ageoI24[] <- round(I24[i] *  p_aging[i])
n_ageoA24[] <- round(A24[i] *  p_aging[i])
n_ageoE31[] <- round(E31[i] *  p_aging[i])
n_ageoI31[] <- round(I31[i] *  p_aging[i])
n_ageoA31[] <- round(A31[i] *  p_aging[i])
n_ageoE32[] <- round(E32[i] *  p_aging[i])
n_ageoI32[] <- round(I32[i] *  p_aging[i])
n_ageoA32[] <- round(A32[i] *  p_aging[i])
n_ageoE34[] <- round(E34[i] *  p_aging[i])
n_ageoI34[] <- round(I34[i] *  p_aging[i])
n_ageoA34[] <- round(A34[i] *  p_aging[i])
n_ageoE41[] <- round(E41[i] *  p_aging[i])
n_ageoI41[] <- round(I41[i] *  p_aging[i])
n_ageoA41[] <- round(A41[i] *  p_aging[i])
n_ageoE42[] <- round(E42[i] *  p_aging[i])
n_ageoI42[] <- round(I42[i] *  p_aging[i])
n_ageoA42[] <- round(A42[i] *  p_aging[i])
n_ageoE43[] <- round(E43[i] *  p_aging[i])
n_ageoI43[] <- round(I43[i] *  p_aging[i])
n_ageoA43[] <- round(A43[i] *  p_aging[i])
n_ageoR12[] <- round(R12[i] *  p_aging[i])
n_ageoR13[] <- round(R13[i] *  p_aging[i])
n_ageoR14[] <- round(R14[i] *  p_aging[i])
n_ageoR23[] <- round(R23[i] *  p_aging[i])
n_ageoR24[] <- round(R24[i] *  p_aging[i])
n_ageoR34[] <- round(R34[i] *  p_aging[i])
n_ageoE123[] <- round(E123[i] *  p_aging[i])
n_ageoI123[] <- round(I123[i] *  p_aging[i])
n_ageoA123[] <- round(A123[i] *  p_aging[i])
n_ageoE124[] <- round(E124[i] *  p_aging[i])
n_ageoI124[] <- round(I124[i] *  p_aging[i])
n_ageoA124[] <- round(A124[i] *  p_aging[i])
n_ageoE134[] <- round(E134[i] *  p_aging[i])
n_ageoI134[] <- round(I134[i] *  p_aging[i])
n_ageoA134[] <- round(A134[i] *  p_aging[i])
n_ageoE213[] <- round(E213[i] *  p_aging[i])
n_ageoI213[] <- round(I213[i] *  p_aging[i])
n_ageoA213[] <- round(A213[i] *  p_aging[i])
n_ageoE214[] <- round(E214[i] *  p_aging[i])
n_ageoI214[] <- round(I214[i] *  p_aging[i])
n_ageoA214[] <- round(A214[i] *  p_aging[i])
n_ageoE234[] <- round(E234[i] *  p_aging[i])
n_ageoI234[] <- round(I234[i] *  p_aging[i])
n_ageoA234[] <- round(A234[i] *  p_aging[i])
n_ageoE312[] <- round(E312[i] *  p_aging[i])
n_ageoI312[] <- round(I312[i] *  p_aging[i])
n_ageoA312[] <- round(A312[i] *  p_aging[i])
n_ageoE314[] <- round(E314[i] *  p_aging[i])
n_ageoI314[] <- round(I314[i] *  p_aging[i])
n_ageoA314[] <- round(A314[i] *  p_aging[i])
n_ageoE324[] <- round(E324[i] *  p_aging[i])
n_ageoI324[] <- round(I324[i] *  p_aging[i])
n_ageoA324[] <- round(A324[i] *  p_aging[i])
n_ageoE412[] <- round(E412[i] *  p_aging[i])
n_ageoI412[] <- round(I412[i] *  p_aging[i])
n_ageoA412[] <- round(A412[i] *  p_aging[i])
n_ageoE413[] <- round(E413[i] *  p_aging[i])
n_ageoI413[] <- round(I413[i] *  p_aging[i])
n_ageoA413[] <- round(A413[i] *  p_aging[i])
n_ageoE423[] <- round(E423[i] *  p_aging[i])
n_ageoI423[] <- round(I423[i] *  p_aging[i])
n_ageoA423[] <- round(A423[i] *  p_aging[i])
n_ageoR123[] <- round(R123[i] *  p_aging[i])
n_ageoR124[] <- round(R124[i] *  p_aging[i])
n_ageoR134[] <- round(R134[i] *  p_aging[i])
n_ageoR234[] <- round(R234[i] *  p_aging[i])
n_ageoE1234[] <- round(E1234[i] *  p_aging[i])
n_ageoI1234[] <- round(I1234[i] *  p_aging[i])
n_ageoA1234[] <- round(A1234[i] *  p_aging[i])
n_ageoE2134[] <- round(E2134[i] *  p_aging[i])
n_ageoI2134[] <- round(I2134[i] *  p_aging[i])
n_ageoA2134[] <- round(A2134[i] *  p_aging[i])
n_ageoE3124[] <- round(E3124[i] *  p_aging[i])
n_ageoI3124[] <- round(I3124[i] *  p_aging[i])
n_ageoA3124[] <- round(A3124[i] *  p_aging[i])
n_ageoE4123[] <- round(E4123[i] *  p_aging[i])
n_ageoI4123[] <- round(I4123[i] *  p_aging[i])
n_ageoA4123[] <- round(A4123[i] *  p_aging[i])
n_ageoR1234[] <- round(R1234[i] *  p_aging[i])


#Age in
n_ageiM[] <- if (i>1) round(M[i-1] * p_aging  [i-1]) else 0
n_ageiG[] <- if (i>1) round(G[i-1] * p_aging  [i-1]) else 0
n_ageiS[] <- if (i>1) round(S[i-1] * p_aging  [i-1]) else 0
n_ageiE1[] <- if (i>1) round(E1[i-1] * p_aging  [i-1]) else 0
n_ageiI1[] <- if (i>1) round(I1[i-1] * p_aging  [i-1]) else 0
n_ageiA1[] <- if (i>1) round(A1[i-1] * p_aging  [i-1]) else 0
n_ageiR1[] <- if (i>1) round(R1[i-1] * p_aging  [i-1]) else 0
n_ageiE2[] <- if (i>1) round(E2[i-1] * p_aging  [i-1]) else 0
n_ageiI2[] <- if (i>1) round(I2[i-1] * p_aging  [i-1]) else 0
n_ageiA2[] <- if (i>1) round(A2[i-1] * p_aging  [i-1]) else 0
n_ageiR2[] <- if (i>1) round(R2[i-1] * p_aging  [i-1]) else 0
n_ageiE3[] <- if (i>1) round(E3[i-1] * p_aging  [i-1]) else 0
n_ageiI3[] <- if (i>1) round(I3[i-1] * p_aging  [i-1]) else 0
n_ageiA3[] <- if (i>1) round(A3[i-1] * p_aging  [i-1]) else 0
n_ageiR3[] <- if (i>1) round(R3[i-1] * p_aging  [i-1]) else 0
n_ageiE4[] <- if (i>1) round(E4[i-1] * p_aging  [i-1]) else 0
n_ageiI4[] <- if (i>1) round(I4[i-1] * p_aging  [i-1]) else 0
n_ageiA4[] <- if (i>1) round(A4[i-1] * p_aging  [i-1]) else 0
n_ageiR4[] <- if (i>1) round(R4[i-1] * p_aging  [i-1]) else 0
n_ageiE12[] <- if (i>1) round(E12[i-1] * p_aging  [i-1]) else 0
n_ageiI12[] <- if (i>1) round(I12[i-1] * p_aging  [i-1]) else 0
n_ageiA12[] <- if (i>1) round(A12[i-1] * p_aging  [i-1]) else 0
n_ageiE13[] <- if (i>1) round(E13[i-1] * p_aging  [i-1]) else 0
n_ageiI13[] <- if (i>1) round(I13[i-1] * p_aging  [i-1]) else 0
n_ageiA13[] <- if (i>1) round(A13[i-1] * p_aging  [i-1]) else 0
n_ageiE14[] <- if (i>1) round(E14[i-1] * p_aging  [i-1]) else 0
n_ageiI14[] <- if (i>1) round(I14[i-1] * p_aging  [i-1]) else 0
n_ageiA14[] <- if (i>1) round(A14[i-1] * p_aging  [i-1]) else 0
n_ageiE21[] <- if (i>1) round(E21[i-1] * p_aging  [i-1]) else 0
n_ageiI21[] <- if (i>1) round(I21[i-1] * p_aging  [i-1]) else 0
n_ageiA21[] <- if (i>1) round(A21[i-1] * p_aging  [i-1]) else 0
n_ageiE23[] <- if (i>1) round(E23[i-1] * p_aging  [i-1]) else 0
n_ageiI23[] <- if (i>1) round(I23[i-1] * p_aging  [i-1]) else 0
n_ageiA23[] <- if (i>1) round(A23[i-1] * p_aging  [i-1]) else 0
n_ageiE24[] <- if (i>1) round(E24[i-1] * p_aging  [i-1]) else 0
n_ageiI24[] <- if (i>1) round(I24[i-1] * p_aging  [i-1]) else 0
n_ageiA24[] <- if (i>1) round(A24[i-1] * p_aging  [i-1]) else 0
n_ageiE31[] <- if (i>1) round(E31[i-1] * p_aging  [i-1]) else 0
n_ageiI31[] <- if (i>1) round(I31[i-1] * p_aging  [i-1]) else 0
n_ageiA31[] <- if (i>1) round(A31[i-1] * p_aging  [i-1]) else 0
n_ageiE32[] <- if (i>1) round(E32[i-1] * p_aging  [i-1]) else 0
n_ageiI32[] <- if (i>1) round(I32[i-1] * p_aging  [i-1]) else 0
n_ageiA32[] <- if (i>1) round(A32[i-1] * p_aging  [i-1]) else 0
n_ageiE34[] <- if (i>1) round(E34[i-1] * p_aging  [i-1]) else 0
n_ageiI34[] <- if (i>1) round(I34[i-1] * p_aging  [i-1]) else 0
n_ageiA34[] <- if (i>1) round(A34[i-1] * p_aging  [i-1]) else 0
n_ageiE41[] <- if (i>1) round(E41[i-1] * p_aging  [i-1]) else 0
n_ageiI41[] <- if (i>1) round(I41[i-1] * p_aging  [i-1]) else 0
n_ageiA41[] <- if (i>1) round(A41[i-1] * p_aging  [i-1]) else 0
n_ageiE42[] <- if (i>1) round(E42[i-1] * p_aging  [i-1]) else 0
n_ageiI42[] <- if (i>1) round(I42[i-1] * p_aging  [i-1]) else 0
n_ageiA42[] <- if (i>1) round(A42[i-1] * p_aging  [i-1]) else 0
n_ageiE43[] <- if (i>1) round(E43[i-1] * p_aging  [i-1]) else 0
n_ageiI43[] <- if (i>1) round(I43[i-1] * p_aging  [i-1]) else 0
n_ageiA43[] <- if (i>1) round(A43[i-1] * p_aging  [i-1]) else 0
n_ageiR12[] <- if (i>1) round(R12[i-1] * p_aging  [i-1]) else 0
n_ageiR13[] <- if (i>1) round(R13[i-1] * p_aging  [i-1]) else 0
n_ageiR14[] <- if (i>1) round(R14[i-1] * p_aging  [i-1]) else 0
n_ageiR23[] <- if (i>1) round(R23[i-1] * p_aging  [i-1]) else 0
n_ageiR24[] <- if (i>1) round(R24[i-1] * p_aging  [i-1]) else 0
n_ageiR34[] <- if (i>1) round(R34[i-1] * p_aging  [i-1]) else 0
n_ageiE123[] <- if (i>1) round(E123[i-1] * p_aging  [i-1]) else 0
n_ageiI123[] <- if (i>1) round(I123[i-1] * p_aging  [i-1]) else 0
n_ageiA123[] <- if (i>1) round(A123[i-1] * p_aging  [i-1]) else 0
n_ageiE124[] <- if (i>1) round(E124[i-1] * p_aging  [i-1]) else 0
n_ageiI124[] <- if (i>1) round(I124[i-1] * p_aging  [i-1]) else 0
n_ageiA124[] <- if (i>1) round(A124[i-1] * p_aging  [i-1]) else 0
n_ageiE134[] <- if (i>1) round(E134[i-1] * p_aging  [i-1]) else 0
n_ageiI134[] <- if (i>1) round(I134[i-1] * p_aging  [i-1]) else 0
n_ageiA134[] <- if (i>1) round(A134[i-1] * p_aging  [i-1]) else 0
n_ageiE213[] <- if (i>1) round(E213[i-1] * p_aging  [i-1]) else 0
n_ageiI213[] <- if (i>1) round(I213[i-1] * p_aging  [i-1]) else 0
n_ageiA213[] <- if (i>1) round(A213[i-1] * p_aging  [i-1]) else 0
n_ageiE214[] <- if (i>1) round(E214[i-1] * p_aging  [i-1]) else 0
n_ageiI214[] <- if (i>1) round(I214[i-1] * p_aging  [i-1]) else 0
n_ageiA214[] <- if (i>1) round(A214[i-1] * p_aging  [i-1]) else 0
n_ageiE234[] <- if (i>1) round(E234[i-1] * p_aging  [i-1]) else 0
n_ageiI234[] <- if (i>1) round(I234[i-1] * p_aging  [i-1]) else 0
n_ageiA234[] <- if (i>1) round(A234[i-1] * p_aging  [i-1]) else 0
n_ageiE312[] <- if (i>1) round(E312[i-1] * p_aging  [i-1]) else 0
n_ageiI312[] <- if (i>1) round(I312[i-1] * p_aging  [i-1]) else 0
n_ageiA312[] <- if (i>1) round(A312[i-1] * p_aging  [i-1]) else 0
n_ageiE314[] <- if (i>1) round(E314[i-1] * p_aging  [i-1]) else 0
n_ageiI314[] <- if (i>1) round(I314[i-1] * p_aging  [i-1]) else 0
n_ageiA314[] <- if (i>1) round(A314[i-1] * p_aging  [i-1]) else 0
n_ageiE324[] <- if (i>1) round(E324[i-1] * p_aging  [i-1]) else 0
n_ageiI324[] <- if (i>1) round(I324[i-1] * p_aging  [i-1]) else 0
n_ageiA324[] <- if (i>1) round(A324[i-1] * p_aging  [i-1]) else 0
n_ageiE412[] <- if (i>1) round(E412[i-1] * p_aging  [i-1]) else 0
n_ageiI412[] <- if (i>1) round(I412[i-1] * p_aging  [i-1]) else 0
n_ageiA412[] <- if (i>1) round(A412[i-1] * p_aging  [i-1]) else 0
n_ageiE413[] <- if (i>1) round(E413[i-1] * p_aging  [i-1]) else 0
n_ageiI413[] <- if (i>1) round(I413[i-1] * p_aging  [i-1]) else 0
n_ageiA413[] <- if (i>1) round(A413[i-1] * p_aging  [i-1]) else 0
n_ageiE423[] <- if (i>1) round(E423[i-1] * p_aging  [i-1]) else 0
n_ageiI423[] <- if (i>1) round(I423[i-1] * p_aging  [i-1]) else 0
n_ageiA423[] <- if (i>1) round(A423[i-1] * p_aging  [i-1]) else 0
n_ageiR123[] <- if (i>1) round(R123[i-1] * p_aging  [i-1]) else 0
n_ageiR124[] <- if (i>1) round(R124[i-1] * p_aging  [i-1]) else 0
n_ageiR134[] <- if (i>1) round(R134[i-1] * p_aging  [i-1]) else 0
n_ageiR234[] <- if (i>1) round(R234[i-1] * p_aging  [i-1]) else 0
n_ageiE1234[] <- if (i>1) round(E1234[i-1] * p_aging  [i-1]) else 0
n_ageiI1234[] <- if (i>1) round(I1234[i-1] * p_aging  [i-1]) else 0
n_ageiA1234[] <- if (i>1) round(A1234[i-1] * p_aging  [i-1]) else 0
n_ageiE2134[] <- if (i>1) round(E2134[i-1] * p_aging  [i-1]) else 0
n_ageiI2134[] <- if (i>1) round(I2134[i-1] * p_aging  [i-1]) else 0
n_ageiA2134[] <- if (i>1) round(A2134[i-1] * p_aging  [i-1]) else 0
n_ageiE3124[] <- if (i>1) round(E3124[i-1] * p_aging  [i-1]) else 0
n_ageiI3124[] <- if (i>1) round(I3124[i-1] * p_aging  [i-1]) else 0
n_ageiA3124[] <- if (i>1) round(A3124[i-1] * p_aging  [i-1]) else 0
n_ageiE4123[] <- if (i>1) round(E4123[i-1] * p_aging  [i-1]) else 0
n_ageiI4123[] <- if (i>1) round(I4123[i-1] * p_aging  [i-1]) else 0
n_ageiA4123[] <- if (i>1) round(A4123[i-1] * p_aging  [i-1]) else 0
n_ageiR1234[] <- if (i>1) round(R1234[i-1] * p_aging  [i-1]) else 0


########### Draws from binomial distributions for numbers changing between
## compartments:
###############################################

## Binomial draw for mortality
n_muM	[]<- rbinom(	M[i] -  n_ageoM[i]  ,p_mu[i])
n_muG	[]<- rbinom(G[i] -  n_ageoG[i]  ,p_mu[i])
n_muS	[]<- rbinom(S[i] -  n_ageoS[i]  ,p_mu[i])
n_muE1	[]<- rbinom(E1[i] -  n_ageoE1[i]  ,p_mu[i])
n_muI1	[]<- rbinom(I1[i] -  n_ageoI1[i]  ,p_mu[i])
n_muA1	[]<- rbinom(A1[i] -  n_ageoA1[i]  ,p_mu[i])
n_muR1	[]<- rbinom(R1[i] -  n_ageoR1[i]  ,p_mu[i])
n_muE2	[]<- rbinom(E2[i] -  n_ageoE2[i]  ,p_mu[i])
n_muI2	[]<- rbinom(I2[i] -  n_ageoI2[i]  ,p_mu[i])
n_muA2	[]<- rbinom(A2[i] -  n_ageoA2[i]  ,p_mu[i])
n_muR2	[]<- rbinom(R2[i] -  n_ageoR2[i]  ,p_mu[i])
n_muE3	[]<- rbinom(E3[i] -  n_ageoE3[i]  ,p_mu[i])
n_muI3	[]<- rbinom(I3[i] -  n_ageoI3[i]  ,p_mu[i])
n_muA3	[]<- rbinom(A3[i] -  n_ageoA3[i]  ,p_mu[i])
n_muR3	[]<- rbinom(R3[i] -  n_ageoR3[i]  ,p_mu[i])
n_muE4	[]<- rbinom(E4[i] -  n_ageoE4[i]  ,p_mu[i])
n_muI4	[]<- rbinom(I4[i] -  n_ageoI4[i]  ,p_mu[i])
n_muA4	[]<- rbinom(A4[i] -  n_ageoA4[i]  ,p_mu[i])
n_muR4	[]<- rbinom(R4[i] -  n_ageoR4[i]  ,p_mu[i])
n_muE12	[]<- rbinom(E12[i] -  n_ageoE12[i]  ,p_mu[i])
n_muI12	[]<- rbinom(I12[i] -  n_ageoI12[i]  ,p_mu[i])
n_muA12	[]<- rbinom(A12[i] -  n_ageoA12[i]  ,p_mu[i])
n_muE13	[]<- rbinom(E13[i] -  n_ageoE13[i]  ,p_mu[i])
n_muI13	[]<- rbinom(I13[i] -  n_ageoI13[i]  ,p_mu[i])
n_muA13	[]<- rbinom(A13[i] -  n_ageoA13[i]  ,p_mu[i])
n_muE14	[]<- rbinom(E14[i] -  n_ageoE14[i]  ,p_mu[i])
n_muI14	[]<- rbinom(I14[i] -  n_ageoI14[i]  ,p_mu[i])
n_muA14	[]<- rbinom(A14[i] -  n_ageoA14[i]  ,p_mu[i])
n_muE21	[]<- rbinom(E21[i] -  n_ageoE21[i]  ,p_mu[i])
n_muI21	[]<- rbinom(I21[i] -  n_ageoI21[i]  ,p_mu[i])
n_muA21	[]<- rbinom(A21[i] -  n_ageoA21[i]  ,p_mu[i])
n_muE23	[]<- rbinom(E23[i] -  n_ageoE23[i]  ,p_mu[i])
n_muI23	[]<- rbinom(I23[i] -  n_ageoI23[i]  ,p_mu[i])
n_muA23	[]<- rbinom(A23[i] -  n_ageoA23[i]  ,p_mu[i])
n_muE24	[]<- rbinom(E24[i] -  n_ageoE24[i]  ,p_mu[i])
n_muI24	[]<- rbinom(I24[i] -  n_ageoI24[i]  ,p_mu[i])
n_muA24	[]<- rbinom(A24[i] -  n_ageoA24[i]  ,p_mu[i])
n_muE31	[]<- rbinom(E31[i] -  n_ageoE31[i]  ,p_mu[i])
n_muI31	[]<- rbinom(I31[i] -  n_ageoI31[i]  ,p_mu[i])
n_muA31	[]<- rbinom(A31[i] -  n_ageoA31[i]  ,p_mu[i])
n_muE32	[]<- rbinom(E32[i] -  n_ageoE32[i]  ,p_mu[i])
n_muI32	[]<- rbinom(I32[i] -  n_ageoI32[i]  ,p_mu[i])
n_muA32	[]<- rbinom(A32[i] -  n_ageoA32[i]  ,p_mu[i])
n_muE34	[]<- rbinom(E34[i] -  n_ageoE34[i]  ,p_mu[i])
n_muI34	[]<- rbinom(I34[i] -  n_ageoI34[i]  ,p_mu[i])
n_muA34	[]<- rbinom(A34[i] -  n_ageoA34[i]  ,p_mu[i])
n_muE41	[]<- rbinom(E41[i] -  n_ageoE41[i]  ,p_mu[i])
n_muI41	[]<- rbinom(I41[i] -  n_ageoI41[i]  ,p_mu[i])
n_muA41	[]<- rbinom(A41[i] -  n_ageoA41[i]  ,p_mu[i])
n_muE42	[]<- rbinom(E42[i] -  n_ageoE42[i]  ,p_mu[i])
n_muI42	[]<- rbinom(I42[i] -  n_ageoI42[i]  ,p_mu[i])
n_muA42	[]<- rbinom(A42[i] -  n_ageoA42[i]  ,p_mu[i])
n_muE43	[]<- rbinom(E43[i] -  n_ageoE43[i]  ,p_mu[i])
n_muI43	[]<- rbinom(I43[i] -  n_ageoI43[i]  ,p_mu[i])
n_muA43	[]<- rbinom(A43[i] -  n_ageoA43[i]  ,p_mu[i])
n_muR12	[]<- rbinom(R12[i] -  n_ageoR12[i]  ,p_mu[i])
n_muR13	[]<- rbinom(R13[i] -  n_ageoR13[i]  ,p_mu[i])
n_muR14	[]<- rbinom(R14[i] -  n_ageoR14[i]  ,p_mu[i])
n_muR23	[]<- rbinom(R23[i] -  n_ageoR23[i]  ,p_mu[i])
n_muR24	[]<- rbinom(R24[i] -  n_ageoR24[i]  ,p_mu[i])
n_muR34	[]<- rbinom(R34[i] -  n_ageoR34[i]  ,p_mu[i])
n_muE123	[]<- rbinom(E123[i] -  n_ageoE123[i]  ,p_mu[i])
n_muI123	[]<- rbinom(I123[i] -  n_ageoI123[i]  ,p_mu[i])
n_muA123	[]<- rbinom(A123[i] -  n_ageoA123[i]  ,p_mu[i])
n_muE124	[]<- rbinom(E124[i] -  n_ageoE124[i]  ,p_mu[i])
n_muI124	[]<- rbinom(I124[i] -  n_ageoI124[i]  ,p_mu[i])
n_muA124	[]<- rbinom(A124[i] -  n_ageoA124[i]  ,p_mu[i])
n_muE134	[]<- rbinom(E134[i] -  n_ageoE134[i]  ,p_mu[i])
n_muI134	[]<- rbinom(I134[i] -  n_ageoI134[i]  ,p_mu[i])
n_muA134	[]<- rbinom(A134[i] -  n_ageoA134[i]  ,p_mu[i])
n_muE213	[]<- rbinom(E213[i] -  n_ageoE213[i]  ,p_mu[i])
n_muI213	[]<- rbinom(I213[i] -  n_ageoI213[i]  ,p_mu[i])
n_muA213	[]<- rbinom(A213[i] -  n_ageoA213[i]  ,p_mu[i])
n_muE214	[]<- rbinom(E214[i] -  n_ageoE214[i]  ,p_mu[i])
n_muI214	[]<- rbinom(I214[i] -  n_ageoI214[i]  ,p_mu[i])
n_muA214	[]<- rbinom(A214[i] -  n_ageoA214[i]  ,p_mu[i])
n_muE234	[]<- rbinom(E234[i] -  n_ageoE234[i]  ,p_mu[i])
n_muI234	[]<- rbinom(I234[i] -  n_ageoI234[i]  ,p_mu[i])
n_muA234	[]<- rbinom(A234[i] -  n_ageoA234[i]  ,p_mu[i])
n_muE312	[]<- rbinom(E312[i] -  n_ageoE312[i]  ,p_mu[i])
n_muI312	[]<- rbinom(I312[i] -  n_ageoI312[i]  ,p_mu[i])
n_muA312	[]<- rbinom(A312[i] -  n_ageoA312[i]  ,p_mu[i])
n_muE314	[]<- rbinom(E314[i] -  n_ageoE314[i]  ,p_mu[i])
n_muI314	[]<- rbinom(I314[i] -  n_ageoI314[i]  ,p_mu[i])
n_muA314	[]<- rbinom(A314[i] -  n_ageoA314[i]  ,p_mu[i])
n_muE324	[]<- rbinom(E324[i] -  n_ageoE324[i]  ,p_mu[i])
n_muI324	[]<- rbinom(I324[i] -  n_ageoI324[i]  ,p_mu[i])
n_muA324	[]<- rbinom(A324[i] -  n_ageoA324[i]  ,p_mu[i])
n_muE412	[]<- rbinom(E412[i] -  n_ageoE412[i]  ,p_mu[i])
n_muI412	[]<- rbinom(I412[i] -  n_ageoI412[i]  ,p_mu[i])
n_muA412	[]<- rbinom(A412[i] -  n_ageoA412[i]  ,p_mu[i])
n_muE413	[]<- rbinom(E413[i] -  n_ageoE413[i]  ,p_mu[i])
n_muI413	[]<- rbinom(I413[i] -  n_ageoI413[i]  ,p_mu[i])
n_muA413	[]<- rbinom(A413[i] -  n_ageoA413[i]  ,p_mu[i])
n_muE423	[]<- rbinom(E423[i] -  n_ageoE423[i]  ,p_mu[i])
n_muI423	[]<- rbinom(I423[i] -  n_ageoI423[i]  ,p_mu[i])
n_muA423	[]<- rbinom(A423[i] -  n_ageoA423[i]  ,p_mu[i])
n_muR123	[]<- rbinom(R123[i] -  n_ageoR123[i]  ,p_mu[i])
n_muR124	[]<- rbinom(R124[i] -  n_ageoR124[i]  ,p_mu[i])
n_muR134	[]<- rbinom(R134[i] -  n_ageoR134[i]  ,p_mu[i])
n_muR234	[]<- rbinom(R234[i] -  n_ageoR234[i]  ,p_mu[i])
n_muE1234	[]<- rbinom(E1234[i] -  n_ageoE1234[i]  ,p_mu[i])
n_muI1234	[]<- rbinom(I1234[i] -  n_ageoI1234[i]  ,p_mu[i])
n_muA1234	[]<- rbinom(A1234[i] -  n_ageoA1234[i]  ,p_mu[i])
n_muE2134	[]<- rbinom(E2134[i] -  n_ageoE2134[i]  ,p_mu[i])
n_muI2134	[]<- rbinom(I2134[i] -  n_ageoI2134[i]  ,p_mu[i])
n_muA2134	[]<- rbinom(A2134[i] -  n_ageoA2134[i]  ,p_mu[i])
n_muE3124	[]<- rbinom(E3124[i] -  n_ageoE3124[i]  ,p_mu[i])
n_muI3124	[]<- rbinom(I3124[i] -  n_ageoI3124[i]  ,p_mu[i])
n_muA3124	[]<- rbinom(A3124[i] -  n_ageoA3124[i]  ,p_mu[i])
n_muE4123	[]<- rbinom(E4123[i] -  n_ageoE4123[i]  ,p_mu[i])
n_muI4123	[]<- rbinom(I4123[i] -  n_ageoI4123[i]  ,p_mu[i])
n_muA4123	[]<- rbinom(A4123[i] -  n_ageoA4123[i]  ,p_mu[i])
n_muR1234	[]<- rbinom(R1234[i] -  n_ageoR1234[i]  ,p_mu[i])



## Initial infection and M to S and vaccinations

n_Gvacc	[]<- rbinom(G[i] -  n_ageoG[i] - n_muG[i]  , p_vacc[i])

n_MS[] <- rbinom(M[i]-n_ageoM[i]-n_muM[i], p_MS[i])
n_Mvacc[]<- rbinom(M[i]-n_ageoM[i]-n_muM[i]-n_MS[i], p_vacc[i])


n_S_E1[] <- rbinom(S[i]-
                     n_ageoS[i]-
                     n_muS[i], p_SE1[i])

n_S_E2[] <- rbinom(S[i]-
                     n_ageoS[i]-
                     n_muS[i]-
                     n_S_E1[i], p_SE2[i])

n_S_E3[] <- rbinom(S[i]-
                     n_ageoS[i]-
                     n_muS[i]-
                     n_S_E1[i]-
                     n_S_E2[i], p_SE3[i])

n_S_E4[] <- rbinom(S[i]-
                     n_ageoS[i]-
                     n_muS[i]-
                     n_S_E1[i]-
                     n_S_E2[i]-
                     n_S_E3[i], p_SE4[i])

n_Svacc[] <- rbinom(S[i]-
                      n_ageoS[i]-
                      n_muS[i]-
                      n_S_E1[i]-
                      n_S_E2[i]-
                      n_S_E3[i]-
                      n_S_E4[i], p_vacc[i])




## Transition E to I

n_E1_I1[]<- rbinom(E1[i]- n_ageoE1[i]- n_muE1[i], p_EI)
n_E12_I12[]<- rbinom(E12[i]- n_ageoE12[i]- n_muE12[i], p_EI)
n_E123_I123[]<- rbinom(E123[i]- n_ageoE123[i]- n_muE123[i], p_EI)
n_E1234_I1234[]<- rbinom(E1234[i]- n_ageoE1234[i]- n_muE1234[i], p_EI)
n_E124_I124[]<- rbinom(E124[i]- n_ageoE124[i]- n_muE124[i], p_EI)
n_E13_I13[]<- rbinom(E13[i]- n_ageoE13[i]- n_muE13[i], p_EI)
n_E134_I134[]<- rbinom(E134[i]- n_ageoE134[i]- n_muE134[i], p_EI)
n_E14_I14[]<- rbinom(E14[i]- n_ageoE14[i]- n_muE14[i], p_EI)
n_E2_I2[]<- rbinom(E2[i]- n_ageoE2[i]- n_muE2[i], p_EI)
n_E21_I21[]<- rbinom(E21[i]- n_ageoE21[i]- n_muE21[i], p_EI)
n_E213_I213[]<- rbinom(E213[i]- n_ageoE213[i]- n_muE213[i], p_EI)
n_E2134_I2134[]<- rbinom(E2134[i]- n_ageoE2134[i]- n_muE2134[i], p_EI)
n_E214_I214[]<- rbinom(E214[i]- n_ageoE214[i]- n_muE214[i], p_EI)
n_E23_I23[]<- rbinom(E23[i]- n_ageoE23[i]- n_muE23[i], p_EI)
n_E234_I234[]<- rbinom(E234[i]- n_ageoE234[i]- n_muE234[i], p_EI)
n_E24_I24[]<- rbinom(E24[i]- n_ageoE24[i]- n_muE24[i], p_EI)
n_E3_I3[]<- rbinom(E3[i]- n_ageoE3[i]- n_muE3[i], p_EI)
n_E31_I31[]<- rbinom(E31[i]- n_ageoE31[i]- n_muE31[i], p_EI)
n_E312_I312[]<- rbinom(E312[i]- n_ageoE312[i]- n_muE312[i], p_EI)
n_E3124_I3124[]<- rbinom(E3124[i]- n_ageoE3124[i]- n_muE3124[i], p_EI)
n_E314_I314[]<- rbinom(E314[i]- n_ageoE314[i]- n_muE314[i], p_EI)
n_E32_I32[]<- rbinom(E32[i]- n_ageoE32[i]- n_muE32[i], p_EI)
n_E324_I324[]<- rbinom(E324[i]- n_ageoE324[i]- n_muE324[i], p_EI)
n_E34_I34[]<- rbinom(E34[i]- n_ageoE34[i]- n_muE34[i], p_EI)
n_E4_I4[]<- rbinom(E4[i]- n_ageoE4[i]- n_muE4[i], p_EI)
n_E41_I41[]<- rbinom(E41[i]- n_ageoE41[i]- n_muE41[i], p_EI)
n_E412_I412[]<- rbinom(E412[i]- n_ageoE412[i]- n_muE412[i], p_EI)
n_E4123_I4123[]<- rbinom(E4123[i]- n_ageoE4123[i]- n_muE4123[i], p_EI)
n_E413_I413[]<- rbinom(E413[i]- n_ageoE413[i]- n_muE413[i], p_EI)
n_E42_I42[]<- rbinom(E42[i]- n_ageoE42[i]- n_muE42[i], p_EI)
n_E423_I423[]<- rbinom(E423[i]- n_ageoE423[i]- n_muE423[i], p_EI)
n_E43_I43[]<- rbinom(E43[i]- n_ageoE43[i]- n_muE43[i], p_EI)



# Transitions I to A

n_I1_A1[]<- rbinom(I1[i]- n_ageoI1[i]- n_muI1[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I12_A12[]<- rbinom(I12[i]- n_ageoI12[i]- n_muI12[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I123_A123[]<- rbinom(I123[i]- n_ageoI123[i]- n_muI123[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I1234_A1234[]<- rbinom(I1234[i]- n_ageoI1234[i]- n_muI1234[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I124_A124[]<- rbinom(I124[i]- n_ageoI124[i]- n_muI124[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I13_A13[]<- rbinom(I13[i]- n_ageoI13[i]- n_muI13[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I134_A134[]<- rbinom(I134[i]- n_ageoI134[i]- n_muI134[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I14_A14[]<- rbinom(I14[i]- n_ageoI14[i]- n_muI14[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I2_A2[]<- rbinom(I2[i]- n_ageoI2[i]- n_muI2[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I21_A21[]<- rbinom(I21[i]- n_ageoI21[i]- n_muI21[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I213_A213[]<- rbinom(I213[i]- n_ageoI213[i]- n_muI213[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I2134_A2134[]<- rbinom(I2134[i]- n_ageoI2134[i]- n_muI2134[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I214_A214[]<- rbinom(I214[i]- n_ageoI214[i]- n_muI214[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I23_A23[]<- rbinom(I23[i]- n_ageoI23[i]- n_muI23[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I234_A234[]<- rbinom(I234[i]- n_ageoI234[i]- n_muI234[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I24_A24[]<- rbinom(I24[i]- n_ageoI24[i]- n_muI24[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I3_A3[]<- rbinom(I3[i]- n_ageoI3[i]- n_muI3[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I31_A31[]<- rbinom(I31[i]- n_ageoI31[i]- n_muI31[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I312_A312[]<- rbinom(I312[i]- n_ageoI312[i]- n_muI312[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I3124_A3124[]<- rbinom(I3124[i]- n_ageoI3124[i]- n_muI3124[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I314_A314[]<- rbinom(I314[i]- n_ageoI314[i]- n_muI314[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I32_A32[]<- rbinom(I32[i]- n_ageoI32[i]- n_muI32[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I324_A324[]<- rbinom(I324[i]- n_ageoI324[i]- n_muI324[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I34_A34[]<- rbinom(I34[i]- n_ageoI34[i]- n_muI34[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I4_A4[]<- rbinom(I4[i]- n_ageoI4[i]- n_muI4[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I41_A41[]<- rbinom(I41[i]- n_ageoI41[i]- n_muI41[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I412_A412[]<- rbinom(I412[i]- n_ageoI412[i]- n_muI412[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I4123_A4123[]<- rbinom(I4123[i]- n_ageoI4123[i]- n_muI4123[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I413_A413[]<- rbinom(I413[i]- n_ageoI413[i]- n_muI413[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I42_A42[]<- rbinom(I42[i]- n_ageoI42[i]- n_muI42[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I423_A423[]<- rbinom(I423[i]- n_ageoI423[i]- n_muI423[i], if (i<=5)p_IA_5 else p_IA_5p )
n_I43_A43[]<- rbinom(I43[i]- n_ageoI43[i]- n_muI43[i], if (i<=5)p_IA_5 else p_IA_5p )



# A to R transitions

n_A1_R1[]<- rbinom(A1[i]- n_ageoA1[i]- n_muA1[i], p_AR)
n_A12_R12[]<- rbinom(A12[i]- n_ageoA12[i]- n_muA12[i], p_AR)
n_A123_R123[]<- rbinom(A123[i]- n_ageoA123[i]- n_muA123[i], p_AR)
n_A1234_R1234[]<- rbinom(A1234[i]- n_ageoA1234[i]- n_muA1234[i], p_AR)
n_A124_R124[]<- rbinom(A124[i]- n_ageoA124[i]- n_muA124[i], p_AR)
n_A13_R13[]<- rbinom(A13[i]- n_ageoA13[i]- n_muA13[i], p_AR)
n_A134_R134[]<- rbinom(A134[i]- n_ageoA134[i]- n_muA134[i], p_AR)
n_A14_R14[]<- rbinom(A14[i]- n_ageoA14[i]- n_muA14[i], p_AR)
n_A2_R2[]<- rbinom(A2[i]- n_ageoA2[i]- n_muA2[i], p_AR)
n_A21_R12[]<- rbinom(A21[i]- n_ageoA21[i]- n_muA21[i], p_AR)
n_A213_R123[]<- rbinom(A213[i]- n_ageoA213[i]- n_muA213[i], p_AR)
n_A2134_R1234[]<- rbinom(A2134[i]- n_ageoA2134[i]- n_muA2134[i], p_AR)
n_A214_R124[]<- rbinom(A214[i]- n_ageoA214[i]- n_muA214[i], p_AR)
n_A23_R23[]<- rbinom(A23[i]- n_ageoA23[i]- n_muA23[i], p_AR)
n_A234_R234[]<- rbinom(A234[i]- n_ageoA234[i]- n_muA234[i], p_AR)
n_A24_R24[]<- rbinom(A24[i]- n_ageoA24[i]- n_muA24[i], p_AR)
n_A3_R3[]<- rbinom(A3[i]- n_ageoA3[i]- n_muA3[i], p_AR)
n_A31_R13[]<- rbinom(A31[i]- n_ageoA31[i]- n_muA31[i], p_AR)
n_A312_R123[]<- rbinom(A312[i]- n_ageoA312[i]- n_muA312[i], p_AR)
n_A3124_R1234[]<- rbinom(A3124[i]- n_ageoA3124[i]- n_muA3124[i], p_AR)
n_A314_R134[]<- rbinom(A314[i]- n_ageoA314[i]- n_muA314[i], p_AR)
n_A32_R23[]<- rbinom(A32[i]- n_ageoA32[i]- n_muA32[i], p_AR)
n_A324_R234[]<- rbinom(A324[i]- n_ageoA324[i]- n_muA324[i], p_AR)
n_A34_R34[]<- rbinom(A34[i]- n_ageoA34[i]- n_muA34[i], p_AR)
n_A4_R4[]<- rbinom(A4[i]- n_ageoA4[i]- n_muA4[i], p_AR)
n_A41_R14[]<- rbinom(A41[i]- n_ageoA41[i]- n_muA41[i], p_AR)
n_A412_R124[]<- rbinom(A412[i]- n_ageoA412[i]- n_muA412[i], p_AR)
n_A4123_R1234[]<- rbinom(A4123[i]- n_ageoA4123[i]- n_muA4123[i], p_AR)
n_A413_R134[]<- rbinom(A413[i]- n_ageoA413[i]- n_muA413[i], p_AR)
n_A42_R24[]<- rbinom(A42[i]- n_ageoA42[i]- n_muA42[i], p_AR)
n_A423_R234[]<- rbinom(A423[i]- n_ageoA423[i]- n_muA423[i], p_AR)
n_A43_R34[]<- rbinom(A43[i]- n_ageoA43[i]- n_muA43[i], p_AR)



# Transitions R to A & E

n_R1_A1[]<- rbinom(	R1[i] - 
                      n_ageoR1[i]- 
                      n_muR1[i], p_R1A1[i])

n_R1_S[]<- rbinom(	R1[i] - 
                     n_ageoR1[i] - 
                     n_muR1[i]-
                     n_R1_A1[i], p_RS)

# Infection 
n_R1_E31[]<- rbinom(	R1[i] - 
                       n_ageoR1[i] - 
                       n_muR1[i]-
                       n_R1_A1[i]-
                       n_R1_S[i], p_SE3[i])

n_R1_E41[]<- rbinom(	R1[i] - 
                       n_ageoR1[i] - 
                       n_muR1[i]-
                       n_R1_A1[i]-
                       n_R1_S[i] -
                       n_R1_E31[i], p_SE4[i])

n_R1_E21[]<- rbinom(	R1[i] - 
                       n_ageoR1[i] - 
                       n_muR1[i]-
                       n_R1_A1[i]-
                       n_R1_E31[i] -
                       n_R1_E41[i] -
                       n_R1_S[i], p_R1_to_E2_gi[i])

n_R1vacc[]<- rbinom(	R1[i] - 
                       n_ageoR1[i] - 
                       n_muR1[i]-
                       n_R1_A1[i]-
                       n_R1_E31[i] -
                       n_R1_E41[i] -
                       n_R1_S[i] - n_R1_E21[i], p_vacc[i])


# n_R1_E31[]<-rbinom(	n_R1_E[i], p_SE3[i]/(p_R1_to_E2_gi[i]+p_SE3[i]+p_SE4[i]))
# n_R1_E41[]<-rbinom(	n_R1_E[i], p_SE4[i]/(p_R1_to_E2_gi[i]+p_SE3[i]+p_SE4[i]))
# n_R1_E21[]<-rbinom(	n_R1_E[i], p_R1_to_E2_gi[i]/(p_R1_to_E2_gi[i]+p_SE3[i]+p_SE4[i]))


n_R2_A2[]<- rbinom(	R2[i] - 
                      n_ageoR2[i]- 
                      n_muR2[i], p_R2A2[i])


n_R2_S[]<- rbinom(	R2[i] - 
                     n_ageoR2[i] - 
                     n_R2_A2[i] -
                     n_muR2[i], p_RS)

n_R2_E12[]<- rbinom(R2[i] - 
                      n_ageoR2[i] - 
                      n_R2_A2[i] -
                      n_muR2[i] -
                      n_R2_S[i], p_R2_to_E1_gi[i])

n_R2_E32[]<- rbinom(R2[i] - 
                      n_ageoR2[i] - 
                      n_R2_A2[i] -
                      n_muR2[i] -
                      n_R2_E12[i] -
                      n_R2_S[i], p_SE3[i])

n_R2_E42[]<- rbinom(R2[i] - 
                      n_ageoR2[i] - 
                      n_R2_A2[i] -
                      n_muR2[i] -
                      n_R2_E12[i] -
                      n_R2_E32[i] -
                      n_R2_S[i], p_SE4[i])

n_R2vacc[]<- rbinom(R2[i] - 
                      n_ageoR2[i] - 
                      n_R2_A2[i] -
                      n_muR2[i] -
                      n_R2_E12[i] -
                      n_R2_E32[i] -
                      n_R2_S[i] - n_R2_E42[i], p_vacc[i])




n_R3_A3[]<- rbinom(	R3[i] - 
                      n_ageoR3[i]- 
                      n_muR3[i], p_R3A3[i])

n_R3_S[]<- rbinom(	R3[i] - 
                     n_ageoR3[i] - 
                     n_R3_A3[i] -
                     n_muR3[i], p_RS)

n_R3_E13[]<- rbinom(R3[i] - 
                      n_ageoR3[i] - 
                      n_R3_A3[i] -
                      n_muR3[i] -
                      n_R3_S[i], p_SE1[i])

n_R3_E23[]<- rbinom(R3[i] - 
                      n_ageoR3[i] - 
                      n_R3_A3[i] -
                      n_muR3[i] -
                      n_R3_E13[i]-
                      n_R3_S[i], p_SE2[i])

n_R3_E43[]<- rbinom(R3[i] - 
                      n_ageoR3[i] - 
                      n_R3_A3[i] -
                      n_muR3[i] -
                      n_R3_E13[i]-
                      n_R3_E23[i] -
                      n_R3_S[i], p_R3_to_E4_gii[i])

n_R3vacc[]<- rbinom(R3[i] - 
                      n_ageoR3[i] - 
                      n_R3_A3[i] -
                      n_muR3[i] -
                      n_R3_E13[i]-
                      n_R3_E23[i] -
                      n_R3_S[i] - n_R3_E43[i], p_vacc[i])



n_R4_A4[]<- rbinom(	R4[i] - 
                      n_ageoR4[i]- 
                      n_muR4[i], p_R4A4[i])

n_R4_S[]<- rbinom(	R4[i] - 
                     n_ageoR4[i] -
                     n_R4_A4[i] -
                     n_muR4[i], p_RS)


n_R4_E14[]<- rbinom(R4[i] - 
                      n_ageoR4[i] -
                      n_R4_A4[i] -
                      n_muR4[i] -
                      n_R4_S[i], p_SE1[i])


n_R4_E24[]<- rbinom(R4[i] - 
                      n_ageoR4[i] -
                      n_R4_A4[i] -
                      n_muR4[i] -
                      n_R4_E14[i] -
                      n_R4_S[i], p_SE2[i])


n_R4_E34[]<- rbinom(R4[i] - 
                      n_ageoR4[i] -
                      n_R4_A4[i] -
                      n_muR4[i] -
                      n_R4_E14[i] -
                      n_R4_E24[i] -
                      n_R4_S[i], p_R4_to_E3_gii[i])

n_R4vacc[]<- rbinom(R4[i] - 
                      n_ageoR4[i] -
                      n_R4_A4[i] -
                      n_muR4[i] -
                      n_R4_E14[i] -
                      n_R4_E24[i] -
                      n_R4_S[i] - n_R4_E34[i], p_vacc[i])




n_R13_A13[]<- rbinom(	R13[i] - 
                        n_ageoR13[i]- 
                        n_muR13[i], p_R1A1[i])

n_R13_A31[]<- rbinom(	R13[i] - 
                        n_ageoR13[i]- 
                        n_R13_A13[i] -
                        n_muR13[i], p_R3A3[i])

n_R13_wane[]<- rbinom(	R13[i] - 
                         n_ageoR13[i] - 
                         n_R13_A13[i] -
                         n_R13_A31[i] -
                         n_muR13[i], p_RS2)

n_R13_E213[]<- rbinom(R13[i] - 
                        n_ageoR13[i] - 
                        n_R13_A13[i] -
                        n_R13_A31[i] -
                        n_R13_wane[i] -
                        n_muR13[i], p_R1_to_E2_gi[i])

n_R13_E413[]<- rbinom(R13[i] - 
                        n_ageoR13[i] - 
                        n_R13_A13[i] -
                        n_R13_A31[i] -
                        n_R13_wane[i] -
                        n_R13_E213[i]-
                        n_muR13[i], p_R3_to_E4_gii[i])

n_R13vacc[]<- rbinom(R13[i] - 
                       n_ageoR13[i] - 
                       n_R13_A13[i] -
                       n_R13_A31[i] -
                       n_R13_wane[i] -
                       n_R13_E213[i]-
                       n_muR13[i] - n_R13_E413[i], p_vacc[i])



n_R12_A12[]<- rbinom(R12[i] - 
                       n_ageoR12[i]- 
                       n_muR12[i], p_R1A1[i])

n_R12_A21[]<- rbinom(R12[i] - 
                       n_ageoR12[i]- 
                       n_muR12[i]-
                       n_R12_A12[i], p_R2A2[i])

n_R12_wane[]<- rbinom(R12[i] - 
                        n_ageoR12[i] -
                        n_muR12[i]-
                        n_R12_A12[i]-
                        n_R12_A21[i], p_RS2)


n_R12_E312[]<- rbinom(R12[i] - 
                        n_ageoR12[i] -
                        n_muR12[i]-
                        n_R12_A12[i]-
                        n_R12_wane[i] -
                        n_R12_A21[i], p_SE3[i])

n_R12_E412[]<- rbinom(R12[i] - 
                        n_ageoR12[i] -
                        n_muR12[i]-
                        n_R12_A12[i]-
                        n_R12_wane[i] -
                        n_R12_E312[i]-
                        n_R12_A21[i], p_SE4[i])

n_R12vacc[]<- rbinom(R12[i] - 
                       n_ageoR12[i] -
                       n_muR12[i]-
                       n_R12_A12[i]-
                       n_R12_wane[i] -
                       n_R12_E312[i]-
                       n_R12_A21[i] - n_R12_E412[i], p_vacc[i])



n_R23_A23[]<- rbinom(	R23[i] - 
                        n_ageoR23[i]- 
                        n_muR23[i], p_R2A2[i])

n_R23_A32[]<- rbinom(	R23[i] - 
                        n_ageoR23[i]- 
                        n_R23_A23[i]-
                        n_muR23[i], p_R3A3[i])

n_R23_wane[]<- rbinom(	R23[i] - 
                         n_ageoR23[i] - 
                         n_R23_A23[i] -
                         n_R23_A32[i] -
                         n_muR23[i], p_RS2)


n_R23_E123[]<- rbinom(R23[i] - 
                        n_ageoR23[i] - 
                        n_R23_A23[i] -
                        n_R23_A32[i] -
                        n_R23_wane[i] -
                        n_muR23[i], p_R2_to_E1_gi[i])

n_R23_E423[]<- rbinom(R23[i] - 
                        n_ageoR23[i] - 
                        n_R23_A23[i] -
                        n_R23_A32[i] -
                        n_R23_wane[i] -
                        n_R23_E123[i]-
                        n_muR23[i], p_R3_to_E4_gii[i])

n_R23vacc[]<- rbinom(R23[i] - 
                       n_ageoR23[i] - 
                       n_R23_A23[i] -
                       n_R23_A32[i] -
                       n_R23_wane[i] -
                       n_R23_E123[i]-
                       n_muR23[i] - n_R23_E423[i], p_vacc[i])



n_R14_A14[]<- rbinom(	R14[i] - 
                        n_ageoR14[i]- 
                        n_muR14[i], p_R1A1[i])

n_R14_A41[]<- rbinom(	R14[i] - 
                        n_ageoR14[i]- 
                        n_R14_A14[i]-
                        n_muR14[i], p_R4A4[i])

n_R14_wane[]<- rbinom(	R14[i] - 
                         n_ageoR14[i] - 
                         n_R14_A14[i] -
                         n_R14_A41[i] -
                         n_muR14[i], p_RS2)


n_R14_E214[]<- rbinom(R14[i] - 
                        n_ageoR14[i] - 
                        n_R14_A14[i] -
                        n_R14_A41[i] -
                        n_R14_wane[i]-
                        n_muR14[i], p_R1_to_E2_gi[i])

n_R14_E314[]<- rbinom(R14[i] - 
                        n_ageoR14[i] - 
                        n_R14_A14[i] -
                        n_R14_A41[i] -
                        n_R14_wane[i]-
                        n_R14_E214[i]-
                        n_muR14[i], p_R4_to_E3_gii[i])

n_R14vacc[]<- rbinom(R14[i] - 
                       n_ageoR14[i] - 
                       n_R14_A14[i] -
                       n_R14_A41[i] -
                       n_R14_wane[i]-
                       n_R14_E214[i]-
                       n_muR14[i] - n_R14_E314[i], p_vacc[i])



n_R24_A24[]<- rbinom(	R24[i] - 
                        n_ageoR24[i]- 
                        n_muR24[i], p_R2A2[i])

n_R24_A42[]<- rbinom(	R24[i] - 
                        n_ageoR24[i]- 
                        n_R24_A24[i]-
                        n_muR24[i], p_R4A4[i])

n_R24_wane[]<- rbinom(	R24[i] - 
                         n_ageoR24[i] - 
                         n_R24_A24[i] -
                         n_R24_A42[i] -
                         n_muR24[i], p_RS2)


n_R24_E124[]<- rbinom(R24[i] - 
                        n_ageoR24[i] - 
                        n_R24_A24[i] -
                        n_R24_A42[i] -
                        n_R24_wane[i]-
                        n_muR24[i], p_R2_to_E1_gi[i])

n_R24_E324[]<- rbinom(R24[i] - 
                        n_ageoR24[i] - 
                        n_R24_A24[i] -
                        n_R24_A42[i] -
                        n_R24_wane[i]-
                        n_R24_E124[i]-
                        n_muR24[i], p_R4_to_E3_gii[i])

n_R24vacc[]<- rbinom(R24[i] - 
                       n_ageoR24[i] - 
                       n_R24_A24[i] -
                       n_R24_A42[i] -
                       n_R24_wane[i]-
                       n_R24_E124[i]-
                       n_muR24[i] - n_R24_E324[i], p_vacc[i])


n_R34_A34[]<- rbinom(	R34[i] - 
                        n_ageoR34[i]- 
                        n_muR34[i], p_R1A1[i])

n_R34_A43[]<- rbinom(	R34[i] - 
                        n_ageoR34[i]- 
                        n_R34_A34[i]-
                        n_muR34[i], p_R1A1[i])

n_R34_wane[]<- rbinom(	R34[i] - 
                         n_ageoR34[i] - 
                         n_R34_A34[i] -
                         n_R34_A43[i] -
                         n_muR34[i], p_RS2)


n_R34_E134[]<- rbinom(R34[i] - 
                        n_ageoR34[i] - 
                        n_R34_A34[i] -
                        n_R34_A43[i] -
                        n_R34_wane[i] -
                        n_muR34[i], p_SE1[i])

n_R34_E234[]<- rbinom(R34[i] - 
                        n_ageoR34[i] - 
                        n_R34_A34[i] -
                        n_R34_A43[i] -
                        n_R34_wane[i] -
                        n_R34_E134[i]-
                        n_muR34[i], p_SE2[i])

n_R34vacc[]<- rbinom(R34[i] - 
                       n_ageoR34[i] - 
                       n_R34_A34[i] -
                       n_R34_A43[i] -
                       n_R34_wane[i] -
                       n_R34_E134[i]-
                       n_muR34[i] - n_R34_E234[i], p_vacc[i])




n_R123_A123[]<- rbinom(	R123[i] - 
                          n_ageoR123[i] - 
                          n_muR123[i], p_R1A1[i])

n_R123_A213[]<- rbinom(	R123[i] - 
                          n_ageoR123[i]- 
                          n_muR123[i]-
                          n_R123_A123[i], p_R2A2[i])

n_R123_A312[]<- rbinom(	R123[i] - 
                          n_ageoR123[i]- 
                          n_muR123[i]-
                          n_R123_A213[i]-
                          n_R123_A123[i],  p_R3A3[i])


n_R123_wane[]<- rbinom(	R123[i] - 
                          n_ageoR123[i] - 
                          n_muR123[i]-
                          n_R123_A123[i]-
                          n_R123_A213[i]-
                          n_R123_A312[i], p_RS3)


n_R123_E4123[]<- rbinom(R123[i] - 
                          n_ageoR123[i] - 
                          n_muR123[i]-
                          n_R123_A123[i]-
                          n_R123_A213[i]-
                          n_R123_A312[i]-
                          n_R123_wane[i] ,p_R3_to_E4_gii[i])

n_R123vacc[]<- rbinom(R123[i] - 
                        n_ageoR123[i] - 
                        n_muR123[i]-
                        n_R123_A123[i]-
                        n_R123_A213[i]-
                        n_R123_A312[i]-
                        n_R123_wane[i] - n_R123_E4123[i] ,p_vacc[i])




n_R124_A124[]<- rbinom(	R124[i] - 
                          n_ageoR124[i]- 
                          n_muR124[i], p_R1A1[i])

n_R124_A214[]<- rbinom(	R124[i] - 
                          n_ageoR124[i] - 
                          n_R124_A124[i] -
                          n_muR124[i], p_R2A2[i])

n_R124_A412[]<- rbinom(	R124[i] - 
                          n_ageoR124[i]- 
                          n_R124_A124[i] -
                          n_R124_A214[i] -
                          n_muR124[i], p_R4A4[i])

n_R124_wane[]<- rbinom(	R124[i] - 
                          n_ageoR124[i]- 
                          n_R124_A124[i] -
                          n_R124_A214[i] -
                          n_R124_A412[i] -
                          n_muR124[i], p_RS3)

n_R124_E3124[]<- rbinom(	R124[i] - 
                           n_ageoR124[i]- 
                           n_R124_A124[i] -
                           n_R124_A214[i] -
                           n_R124_A412[i] -
                           n_muR124[i]-
                           n_R124_wane[i],p_R4_to_E3_gii[i])

n_R124vacc[]<- rbinom(	R124[i] - 
                         n_ageoR124[i]- 
                         n_R124_A124[i] -
                         n_R124_A214[i] -
                         n_R124_A412[i] -
                         n_muR124[i]-
                         n_R124_wane[i]- n_R124_E3124[i],p_vacc[i])



n_R134_A134[]<- rbinom(	R134[i] - 
                          n_ageoR134[i]- 
                          n_muR134[i], p_R1A1[i])

n_R134_A413[]<- rbinom(	R134[i] - 
                          n_ageoR134[i]- 
                          n_R134_A134[i]-
                          n_muR134[i], p_R4A4[i])


n_R134_A314[]<- rbinom(	R134[i] - 
                          n_ageoR134[i]- 
                          n_R134_A134[i]-
                          n_R134_A413[i]-
                          n_muR134[i], p_R3A3[i])

n_R134_wane[]<- rbinom(	R134[i] - 
                          n_ageoR134[i] -
                          n_R134_A134[i] -
                          n_R134_A413[i] -
                          n_R134_A314[i] -
                          n_muR134[i], p_RS3)

n_R134_E2134[]<- rbinom(R134[i] - 
                          n_ageoR134[i] -
                          n_R134_A134[i] -
                          n_R134_A413[i] -
                          n_R134_A314[i] -
                          n_R134_wane[i]-
                          n_muR134[i], p_R1_to_E2_gi[i])

n_R134vacc[]<- rbinom(R134[i] - 
                        n_ageoR134[i] -
                        n_R134_A134[i] -
                        n_R134_A413[i] -
                        n_R134_A314[i] -
                        n_R134_wane[i] -
                        n_muR134[i] - n_R134_E2134[i], p_vacc[i])




n_R234_A234[]<- rbinom(	R234[i] - 
                          n_ageoR234[i]- 
                          n_muR234[i], p_R2A2[i])

n_R234_A423[]<- rbinom(	R234[i] - 
                          n_ageoR234[i]- 
                          n_R234_A234[i]-
                          n_muR234[i], p_R4A4[i])

n_R234_A324[]<- rbinom(	R234[i] - 
                          n_ageoR234[i]- 
                          n_R234_A234[i]-
                          n_R234_A423[i]-
                          n_muR234[i], p_R3A3[i])



n_R234_wane[]<- rbinom(	R234[i] - 
                          n_ageoR234[i]- 
                          n_R234_A234[i]-
                          n_R234_A423[i]-
                          n_R234_A324[i]-
                          n_muR234[i], p_RS3)


n_R234_E1234[]<- rbinom(R234[i] - 
                          n_ageoR234[i]- 
                          n_R234_A234[i]-
                          n_R234_A423[i]-
                          n_R234_A324[i]-
                          n_R234_wane[i] -
                          n_muR234[i], p_R2_to_E1_gi[i])

n_R234vacc[]<- rbinom(R234[i] - 
                        n_ageoR234[i]- 
                        n_R234_A234[i]-
                        n_R234_A423[i]-
                        n_R234_A324[i]-
                        n_R234_wane[i] -
                        n_muR234[i] - n_R234_E1234[i], p_vacc[i])




n_R1234_A1234[]<- rbinom(	R1234[i] - 
                            n_ageoR1234[i]- 
                            n_muR1234[i], p_R1A1[i])

n_R1234_A2134[]<- rbinom(	R1234[i] - 
                            n_ageoR1234[i]- 
                            n_R1234_A1234[i]-
                            n_muR1234[i], p_R2A2[i])

n_R1234_A3124[]<- rbinom(	R1234[i] - 
                            n_ageoR1234[i] - 
                            n_R1234_A1234[i]-
                            n_R1234_A2134[i]-
                            n_muR1234[i], p_R3A3[i])

n_R1234_A4123[]<- rbinom(	R1234[i] - 
                            n_ageoR1234[i] - 
                            n_R1234_A1234[i]-
                            n_R1234_A2134[i]-
                            n_R1234_A3124[i]-
                            n_muR1234[i], p_R4A4[i])

n_R1234_wane[]<- rbinom(	R1234[i] - 
                           n_ageoR1234[i] - 
                           n_R1234_A1234[i]-
                           n_R1234_A2134[i]-
                           n_R1234_A3124[i]-
                           n_R1234_A4123[i]-
                           n_muR1234[i], p_RS4)

n_R1234vacc[]<- rbinom(	R1234[i] - 
                          n_ageoR1234[i] - 
                          n_R1234_A1234[i]-
                          n_R1234_A2134[i]-
                          n_R1234_A3124[i]-
                          n_R1234_A4123[i]-
                          n_muR1234[i] - n_R1234_wane[i], p_vacc[i])


###############
######### Compund vaccine movements by immunity strain

n_toR13vacc[] <- n_Svacc[i] + n_R1vacc[i] + n_R3vacc[i] + n_R13vacc[i]

n_toR123vacc[] <- n_R2vacc[i] + n_R12vacc[i] + n_R23vacc[i] + n_R123vacc[i]

n_toR134vacc[] <- n_R4vacc[i] + n_R14vacc[i] + n_R34vacc[i] + n_R134vacc[i]  

n_toR1234vacc[] <- n_R24vacc[i] + n_R124vacc[i] + n_R234vacc[i] + n_R1234vacc[i]



##### Deaths 
###################################

n_allDeath<- 
  sum(	n_muM	)+
  sum(	n_muG	)+
  sum(	n_muS	)+
  sum(	n_muE1	)+
  sum(	n_muI1	)+
  sum(	n_muA1	)+
  sum(	n_muR1	)+
  sum(	n_muE2	)+
  sum(	n_muI2	)+
  sum(	n_muA2	)+
  sum(	n_muR2	)+
  sum(	n_muE3	)+
  sum(	n_muI3	)+
  sum(	n_muA3	)+
  sum(	n_muR3	)+
  sum(	n_muE4	)+
  sum(	n_muI4	)+
  sum(	n_muA4	)+
  sum(	n_muR4	)+
  sum(	n_muE12	)+
  sum(	n_muI12	)+
  sum(	n_muA12	)+
  sum(	n_muE13	)+
  sum(	n_muI13	)+
  sum(	n_muA13	)+
  sum(	n_muE14	)+
  sum(	n_muI14	)+
  sum(	n_muA14	)+
  sum(	n_muE21	)+
  sum(	n_muI21	)+
  sum(	n_muA21	)+
  sum(	n_muE23	)+
  sum(	n_muI23	)+
  sum(	n_muA23	)+
  sum(	n_muE24	)+
  sum(	n_muI24	)+
  sum(	n_muA24	)+
  sum(	n_muE31	)+
  sum(	n_muI31	)+
  sum(	n_muA31	)+
  sum(	n_muE32	)+
  sum(	n_muI32	)+
  sum(	n_muA32	)+
  sum(	n_muE34	)+
  sum(	n_muI34	)+
  sum(	n_muA34	)+
  sum(	n_muE41	)+
  sum(	n_muI41	)+
  sum(	n_muA41	)+
  sum(	n_muE42	)+
  sum(	n_muI42	)+
  sum(	n_muA42	)+
  sum(	n_muE43	)+
  sum(	n_muI43	)+
  sum(	n_muA43	)+
  sum(	n_muR12	)+
  sum(	n_muR13	)+
  sum(	n_muR14	)+
  sum(	n_muR23	)+
  sum(	n_muR24	)+
  sum(	n_muR34	)+
  sum(	n_muE123	)+
  sum(	n_muI123	)+
  sum(	n_muA123	)+
  sum(	n_muE124	)+
  sum(	n_muI124	)+
  sum(	n_muA124	)+
  sum(	n_muE134	)+
  sum(	n_muI134	)+
  sum(	n_muA134	)+
  sum(	n_muE213	)+
  sum(	n_muI213	)+
  sum(	n_muA213	)+
  sum(	n_muE214	)+
  sum(	n_muI214	)+
  sum(	n_muA214	)+
  sum(	n_muE234	)+
  sum(	n_muI234	)+
  sum(	n_muA234	)+
  sum(	n_muE312	)+
  sum(	n_muI312	)+
  sum(	n_muA312	)+
  sum(	n_muE314	)+
  sum(	n_muI314	)+
  sum(	n_muA314	)+
  sum(	n_muE324	)+
  sum(	n_muI324	)+
  sum(	n_muA324	)+
  sum(	n_muE412	)+
  sum(	n_muI412	)+
  sum(	n_muA412	)+
  sum(	n_muE413	)+
  sum(	n_muI413	)+
  sum(	n_muA413	)+
  sum(	n_muE423	)+
  sum(	n_muI423	)+
  sum(	n_muA423	)+
  sum(	n_muR123	)+
  sum(	n_muR124	)+
  sum(	n_muR134	)+
  sum(	n_muR234	)+
  sum(	n_muE1234	)+
  sum(	n_muI1234	)+
  sum(	n_muA1234	)+
  sum(	n_muE2134	)+
  sum(	n_muI2134	)+
  sum(	n_muA2134	)+
  sum(	n_muE3124	)+
  sum(	n_muI3124	)+
  sum(	n_muA3124	)+
  sum(	n_muE4123	)+
  sum(	n_muI4123	)+
  sum(	n_muA4123	)+
  sum(	n_muR1234	)





# Births to keep stable population equal to deaths
n_bM[] <- if (i==1) n_allDeath*p_birth[1] else 0 
n_bG[] <- if (i==1) n_allDeath*p_birth[2] else 0 
n_bS[] <- if (i==1) n_allDeath*p_birth[3] else 0 






## Initial states:

initial(M[])<- init[	1	,i]
initial(G[])<- init[	2	,i]
initial(S[])<- init[	3	,i]
initial(E1[])<- init[	4	,i]
initial(I1[])<- init[	5	,i]
initial(A1[])<- init[	6	,i]
initial(R1[])<- init[	7	,i]
initial(E2[])<- init[	8	,i]
initial(I2[])<- init[	9	,i]
initial(A2[])<- init[	10	,i]
initial(R2[])<- init[	11	,i]
initial(E3[])<- init[	12	,i]
initial(I3[])<- init[	13	,i]
initial(A3[])<- init[	14	,i]
initial(R3[])<- init[	15	,i]
initial(E4[])<- init[	16	,i]
initial(I4[])<- init[	17	,i]
initial(A4[])<- init[	18	,i]
initial(R4[])<- init[	19	,i]
initial(E12[])<- init[	20	,i]
initial(I12[])<- init[	21	,i]
initial(A12[])<- init[	22	,i]
initial(E13[])<- init[	23	,i]
initial(I13[])<- init[	24	,i]
initial(A13[])<- init[	25	,i]
initial(E14[])<- init[	26	,i]
initial(I14[])<- init[	27	,i]
initial(A14[])<- init[	28	,i]
initial(E21[])<- init[	29	,i]
initial(I21[])<- init[	30	,i]
initial(A21[])<- init[	31	,i]
initial(E23[])<- init[	32	,i]
initial(I23[])<- init[	33	,i]
initial(A23[])<- init[	34	,i]
initial(E24[])<- init[	35	,i]
initial(I24[])<- init[	36	,i]
initial(A24[])<- init[	37	,i]
initial(E31[])<- init[	38	,i]
initial(I31[])<- init[	39	,i]
initial(A31[])<- init[	40	,i]
initial(E32[])<- init[	41	,i]
initial(I32[])<- init[	42	,i]
initial(A32[])<- init[	43	,i]
initial(E34[])<- init[	44	,i]
initial(I34[])<- init[	45	,i]
initial(A34[])<- init[	46	,i]
initial(E41[])<- init[	47	,i]
initial(I41[])<- init[	48	,i]
initial(A41[])<- init[	49	,i]
initial(E42[])<- init[	50	,i]
initial(I42[])<- init[	51	,i]
initial(A42[])<- init[	52	,i]
initial(E43[])<- init[	53	,i]
initial(I43[])<- init[	54	,i]
initial(A43[])<- init[	55	,i]
initial(R12[])<- init[	56	,i]
initial(R13[])<- init[	57	,i]
initial(R14[])<- init[	58	,i]
initial(R23[])<- init[	59	,i]
initial(R24[])<- init[	60	,i]
initial(R34[])<- init[	61	,i]
initial(E123[])<- init[	62	,i]
initial(I123[])<- init[	63	,i]
initial(A123[])<- init[	64	,i]
initial(E124[])<- init[	65	,i]
initial(I124[])<- init[	66	,i]
initial(A124[])<- init[	67	,i]
initial(E134[])<- init[	68	,i]
initial(I134[])<- init[	69	,i]
initial(A134[])<- init[	70	,i]
initial(E213[])<- init[	71	,i]
initial(I213[])<- init[	72	,i]
initial(A213[])<- init[	73	,i]
initial(E214[])<- init[	74	,i]
initial(I214[])<- init[	75	,i]
initial(A214[])<- init[	76	,i]
initial(E234[])<- init[	77	,i]
initial(I234[])<- init[	78	,i]
initial(A234[])<- init[	79	,i]
initial(E312[])<- init[	80	,i]
initial(I312[])<- init[	81	,i]
initial(A312[])<- init[	82	,i]
initial(E314[])<- init[	83	,i]
initial(I314[])<- init[	84	,i]
initial(A314[])<- init[	85	,i]
initial(E324[])<- init[	86	,i]
initial(I324[])<- init[	87	,i]
initial(A324[])<- init[	88	,i]
initial(E412[])<- init[	89	,i]
initial(I412[])<- init[	90	,i]
initial(A412[])<- init[	91	,i]
initial(E413[])<- init[	92	,i]
initial(I413[])<- init[	93	,i]
initial(A413[])<- init[	94	,i]
initial(E423[])<- init[	95	,i]
initial(I423[])<- init[	96	,i]
initial(A423[])<- init[	97	,i]
initial(R123[])<- init[	98	,i]
initial(R124[])<- init[	99	,i]
initial(R134[])<- init[	100	,i]
initial(R234[])<- init[	101	,i]
initial(E1234[])<- init[	102	,i]
initial(I1234[])<- init[	103	,i]
initial(A1234[])<- init[	104	,i]
initial(E2134[])<- init[	105	,i]
initial(I2134[])<- init[	106	,i]
initial(A2134[])<- init[	107	,i]
initial(E3124[])<- init[	108	,i]
initial(I3124[])<- init[	109	,i]
initial(A3124[])<- init[	110	,i]
initial(E4123[])<- init[	111	,i]
initial(I4123[])<- init[	112	,i]
initial(A4123[])<- init[	113	,i]
initial(R1234[])<- init[	114	,i]



# Outputs
initial(inc_day_gii4[]) <- 0
initial(inc_day_gii[]) <- 0
initial(inc_day_gi3[]) <- 0
initial(inc_day_gi[]) <- 0
initial(inc_day_all[]) <- 0
initial(inc_year_gii4[]) <- 0
initial(inc_year_gii[]) <- 0
initial(inc_year_gi3[]) <- 0
initial(inc_year_gi[]) <- 0

initial(pop_by4age[]) <- 0
initial(seroprev[]) <-0
initial(seroprev_num[]) <-0
initial(seroprev_den[]) <-1
initial(infections_day_gi3[])<-0
initial(infections_day_gi[])<-0
initial(infections_day_gii4[])<-0
initial(infections_day_gii[])<-0
initial(reported_wk[])<-0
initial(vaccines_perday)<-0
initial(covid_timeline)<-0
initial(new_cases)<-0
initial(new_cases_week)<-0
initial(new_cases_week_gi)<-0
initial(new_cases_week_gii)<-0
initial(new_cases_week_gi3)<-0
initial(new_cases_week_gii4)<-0
initial(seasonality)<-0
initial(sus_gi)<-0
initial(sus_gii)<-0
initial(sus_gi3)<-0
initial(sus_gii4)<-0
#initial(infections_tot)<-0
# initial(reported_wk0_4)<-0
# initial(reported_wk5_65)<-0
# initial(reported_wk65_p)<-0

## User defined parameters - default in parentheses:
init[,] <- user()

beta_1 <- user(0.2)   # transm coefficient
beta_2 <- user(0.2)   # transm coefficient
beta_3 <- user(0.2)   # transm coefficient
beta_4 <- user(0.2)   # transm coefficient
aduRR<-user(0.1) # ratio adult to child infectiousness

repfac_0<- user(287)    # reported to community factor
repfac_5<- user(287)    # reported to community factor
repfac_15<- user(287)    # reported to community factor
repfac_65p<-user(287)
maternalAB <- user(180)  # maternal Ab decay (days)
epsilon <- user(1.0)   # incubation
theta_5 <- user(2.5)   # duration symptoms in under 5
theta_5p <- user(1.5)   # duration symptoms in over 5
sigma <- user(15) # duration asymp shedding
imm_yr   <- user(5.1)     # duration immunity
imm_fac1   <- user(2)     # duration immunity
imm_fac2   <- user(1)     # duration immunity
imm_fac3   <- user(1)     # duration immunity
rr_inf_asymp   <- user(0.05) # rel infect asymptomatic 
crossp_GI<-user(0.05) # cross_protection prob from j to k
crossp_GII<-user(0.05) # cross_protection prob from j to k
p_nonsecretor<-user(0.2) # Fraction immune genetically
w1_1 <- user(0.8) # sesonality
w2 <-  user(-2) # sesonality
pi <-user(3.141593)
alpha<-user(0)# rel susc in R 
mu[]  <- user()      # mortality rates 
aging_vec[]<-user() # aging transitions matrix
school_step[]<-user()
n_school_steps<-user()
covid_step[]<-user()
n_covid_steps<-user()

# Vaccinatons
vac_eff<- user(0.95)
vac_camp_cov[] <- user() 
vac_sche_cov[]<- user()



#
# dimensions of arrays
N_age <- user()
dim(school_step)<-user()
dim(covid_step)<-user()
dim(init) <-  c(114,N_age)
dim(aging_vec)<- N_age
dim(M	)<- N_age
dim(G	)<- N_age
dim(S	)<- N_age
dim(E1	)<- N_age
dim(I1	)<- N_age
dim(A1	)<- N_age
dim(R1	)<- N_age
dim(E2	)<- N_age
dim(I2	)<- N_age
dim(A2	)<- N_age
dim(R2	)<- N_age
dim(E3	)<- N_age
dim(I3	)<- N_age
dim(A3	)<- N_age
dim(R3	)<- N_age
dim(E4	)<- N_age
dim(I4	)<- N_age
dim(A4	)<- N_age
dim(R4	)<- N_age
dim(E12	)<- N_age
dim(I12	)<- N_age
dim(A12	)<- N_age
dim(E13	)<- N_age
dim(I13	)<- N_age
dim(A13	)<- N_age
dim(E14	)<- N_age
dim(I14	)<- N_age
dim(A14	)<- N_age
dim(E21	)<- N_age
dim(I21	)<- N_age
dim(A21	)<- N_age
dim(E23	)<- N_age
dim(I23	)<- N_age
dim(A23	)<- N_age
dim(E24	)<- N_age
dim(I24	)<- N_age
dim(A24	)<- N_age
dim(E31	)<- N_age
dim(I31	)<- N_age
dim(A31	)<- N_age
dim(E32	)<- N_age
dim(I32	)<- N_age
dim(A32	)<- N_age
dim(E34	)<- N_age
dim(I34	)<- N_age
dim(A34	)<- N_age
dim(E41	)<- N_age
dim(I41	)<- N_age
dim(A41	)<- N_age
dim(E42	)<- N_age
dim(I42	)<- N_age
dim(A42	)<- N_age
dim(E43	)<- N_age
dim(I43	)<- N_age
dim(A43	)<- N_age
dim(R12	)<- N_age
dim(R13	)<- N_age
dim(R14	)<- N_age
dim(R23	)<- N_age
dim(R24	)<- N_age
dim(R34	)<- N_age
dim(E123	)<- N_age
dim(I123	)<- N_age
dim(A123	)<- N_age
dim(E124	)<- N_age
dim(I124	)<- N_age
dim(A124	)<- N_age
dim(E134	)<- N_age
dim(I134	)<- N_age
dim(A134	)<- N_age
dim(E213	)<- N_age
dim(I213	)<- N_age
dim(A213	)<- N_age
dim(E214	)<- N_age
dim(I214	)<- N_age
dim(A214	)<- N_age
dim(E234	)<- N_age
dim(I234	)<- N_age
dim(A234	)<- N_age
dim(E312	)<- N_age
dim(I312	)<- N_age
dim(A312	)<- N_age
dim(E314	)<- N_age
dim(I314	)<- N_age
dim(A314	)<- N_age
dim(E324	)<- N_age
dim(I324	)<- N_age
dim(A324	)<- N_age
dim(E412	)<- N_age
dim(I412	)<- N_age
dim(A412	)<- N_age
dim(E413	)<- N_age
dim(I413	)<- N_age
dim(A413	)<- N_age
dim(E423	)<- N_age
dim(I423	)<- N_age
dim(A423	)<- N_age
dim(R123	)<- N_age
dim(R124	)<- N_age
dim(R134	)<- N_age
dim(R234	)<- N_age
dim(E1234	)<- N_age
dim(I1234	)<- N_age
dim(A1234	)<- N_age
dim(E2134	)<- N_age
dim(I2134	)<- N_age
dim(A2134	)<- N_age
dim(E3124	)<- N_age
dim(I3124	)<- N_age
dim(A3124	)<- N_age
dim(E4123	)<- N_age
dim(I4123	)<- N_age
dim(A4123	)<- N_age
dim(R1234	)<- N_age


dim(n_ageiM	)<- N_age
dim(n_ageiG	)<- N_age
dim(n_ageiS	)<- N_age
dim(n_ageiE1	)<- N_age
dim(n_ageiI1	)<- N_age
dim(n_ageiA1	)<- N_age
dim(n_ageiR1	)<- N_age
dim(n_ageiE2	)<- N_age
dim(n_ageiI2	)<- N_age
dim(n_ageiA2	)<- N_age
dim(n_ageiR2	)<- N_age
dim(n_ageiE3	)<- N_age
dim(n_ageiI3	)<- N_age
dim(n_ageiA3	)<- N_age
dim(n_ageiR3	)<- N_age
dim(n_ageiE4	)<- N_age
dim(n_ageiI4	)<- N_age
dim(n_ageiA4	)<- N_age
dim(n_ageiR4	)<- N_age
dim(n_ageiE12	)<- N_age
dim(n_ageiI12	)<- N_age
dim(n_ageiA12	)<- N_age
dim(n_ageiE13	)<- N_age
dim(n_ageiI13	)<- N_age
dim(n_ageiA13	)<- N_age
dim(n_ageiE14	)<- N_age
dim(n_ageiI14	)<- N_age
dim(n_ageiA14	)<- N_age
dim(n_ageiE21	)<- N_age
dim(n_ageiI21	)<- N_age
dim(n_ageiA21	)<- N_age
dim(n_ageiE23	)<- N_age
dim(n_ageiI23	)<- N_age
dim(n_ageiA23	)<- N_age
dim(n_ageiE24	)<- N_age
dim(n_ageiI24	)<- N_age
dim(n_ageiA24	)<- N_age
dim(n_ageiE31	)<- N_age
dim(n_ageiI31	)<- N_age
dim(n_ageiA31	)<- N_age
dim(n_ageiE32	)<- N_age
dim(n_ageiI32	)<- N_age
dim(n_ageiA32	)<- N_age
dim(n_ageiE34	)<- N_age
dim(n_ageiI34	)<- N_age
dim(n_ageiA34	)<- N_age
dim(n_ageiE41	)<- N_age
dim(n_ageiI41	)<- N_age
dim(n_ageiA41	)<- N_age
dim(n_ageiE42	)<- N_age
dim(n_ageiI42	)<- N_age
dim(n_ageiA42	)<- N_age
dim(n_ageiE43	)<- N_age
dim(n_ageiI43	)<- N_age
dim(n_ageiA43	)<- N_age
dim(n_ageiR12	)<- N_age
dim(n_ageiR13	)<- N_age
dim(n_ageiR14	)<- N_age
dim(n_ageiR23	)<- N_age
dim(n_ageiR24	)<- N_age
dim(n_ageiR34	)<- N_age
dim(n_ageiE123	)<- N_age
dim(n_ageiI123	)<- N_age
dim(n_ageiA123	)<- N_age
dim(n_ageiE124	)<- N_age
dim(n_ageiI124	)<- N_age
dim(n_ageiA124	)<- N_age
dim(n_ageiE134	)<- N_age
dim(n_ageiI134	)<- N_age
dim(n_ageiA134	)<- N_age
dim(n_ageiE213	)<- N_age
dim(n_ageiI213	)<- N_age
dim(n_ageiA213	)<- N_age
dim(n_ageiE214	)<- N_age
dim(n_ageiI214	)<- N_age
dim(n_ageiA214	)<- N_age
dim(n_ageiE234	)<- N_age
dim(n_ageiI234	)<- N_age
dim(n_ageiA234	)<- N_age
dim(n_ageiE312	)<- N_age
dim(n_ageiI312	)<- N_age
dim(n_ageiA312	)<- N_age
dim(n_ageiE314	)<- N_age
dim(n_ageiI314	)<- N_age
dim(n_ageiA314	)<- N_age
dim(n_ageiE324	)<- N_age
dim(n_ageiI324	)<- N_age
dim(n_ageiA324	)<- N_age
dim(n_ageiE412	)<- N_age
dim(n_ageiI412	)<- N_age
dim(n_ageiA412	)<- N_age
dim(n_ageiE413	)<- N_age
dim(n_ageiI413	)<- N_age
dim(n_ageiA413	)<- N_age
dim(n_ageiE423	)<- N_age
dim(n_ageiI423	)<- N_age
dim(n_ageiA423	)<- N_age
dim(n_ageiR123	)<- N_age
dim(n_ageiR124	)<- N_age
dim(n_ageiR134	)<- N_age
dim(n_ageiR234	)<- N_age
dim(n_ageiE1234	)<- N_age
dim(n_ageiI1234	)<- N_age
dim(n_ageiA1234	)<- N_age
dim(n_ageiE2134	)<- N_age
dim(n_ageiI2134	)<- N_age
dim(n_ageiA2134	)<- N_age
dim(n_ageiE3124	)<- N_age
dim(n_ageiI3124	)<- N_age
dim(n_ageiA3124	)<- N_age
dim(n_ageiE4123	)<- N_age
dim(n_ageiI4123	)<- N_age
dim(n_ageiA4123	)<- N_age
dim(n_ageiR1234	)<- N_age

dim(n_ageoM	)<- N_age
dim(n_ageoG	)<- N_age
dim(n_ageoS	)<- N_age
dim(n_ageoE1	)<- N_age
dim(n_ageoI1	)<- N_age
dim(n_ageoA1	)<- N_age
dim(n_ageoR1	)<- N_age
dim(n_ageoE2	)<- N_age
dim(n_ageoI2	)<- N_age
dim(n_ageoA2	)<- N_age
dim(n_ageoR2	)<- N_age
dim(n_ageoE3	)<- N_age
dim(n_ageoI3	)<- N_age
dim(n_ageoA3	)<- N_age
dim(n_ageoR3	)<- N_age
dim(n_ageoE4	)<- N_age
dim(n_ageoI4	)<- N_age
dim(n_ageoA4	)<- N_age
dim(n_ageoR4	)<- N_age
dim(n_ageoE12	)<- N_age
dim(n_ageoI12	)<- N_age
dim(n_ageoA12	)<- N_age
dim(n_ageoE13	)<- N_age
dim(n_ageoI13	)<- N_age
dim(n_ageoA13	)<- N_age
dim(n_ageoE14	)<- N_age
dim(n_ageoI14	)<- N_age
dim(n_ageoA14	)<- N_age
dim(n_ageoE21	)<- N_age
dim(n_ageoI21	)<- N_age
dim(n_ageoA21	)<- N_age
dim(n_ageoE23	)<- N_age
dim(n_ageoI23	)<- N_age
dim(n_ageoA23	)<- N_age
dim(n_ageoE24	)<- N_age
dim(n_ageoI24	)<- N_age
dim(n_ageoA24	)<- N_age
dim(n_ageoE31	)<- N_age
dim(n_ageoI31	)<- N_age
dim(n_ageoA31	)<- N_age
dim(n_ageoE32	)<- N_age
dim(n_ageoI32	)<- N_age
dim(n_ageoA32	)<- N_age
dim(n_ageoE34	)<- N_age
dim(n_ageoI34	)<- N_age
dim(n_ageoA34	)<- N_age
dim(n_ageoE41	)<- N_age
dim(n_ageoI41	)<- N_age
dim(n_ageoA41	)<- N_age
dim(n_ageoE42	)<- N_age
dim(n_ageoI42	)<- N_age
dim(n_ageoA42	)<- N_age
dim(n_ageoE43	)<- N_age
dim(n_ageoI43	)<- N_age
dim(n_ageoA43	)<- N_age
dim(n_ageoR12	)<- N_age
dim(n_ageoR13	)<- N_age
dim(n_ageoR14	)<- N_age
dim(n_ageoR23	)<- N_age
dim(n_ageoR24	)<- N_age
dim(n_ageoR34	)<- N_age
dim(n_ageoE123	)<- N_age
dim(n_ageoI123	)<- N_age
dim(n_ageoA123	)<- N_age
dim(n_ageoE124	)<- N_age
dim(n_ageoI124	)<- N_age
dim(n_ageoA124	)<- N_age
dim(n_ageoE134	)<- N_age
dim(n_ageoI134	)<- N_age
dim(n_ageoA134	)<- N_age
dim(n_ageoE213	)<- N_age
dim(n_ageoI213	)<- N_age
dim(n_ageoA213	)<- N_age
dim(n_ageoE214	)<- N_age
dim(n_ageoI214	)<- N_age
dim(n_ageoA214	)<- N_age
dim(n_ageoE234	)<- N_age
dim(n_ageoI234	)<- N_age
dim(n_ageoA234	)<- N_age
dim(n_ageoE312	)<- N_age
dim(n_ageoI312	)<- N_age
dim(n_ageoA312	)<- N_age
dim(n_ageoE314	)<- N_age
dim(n_ageoI314	)<- N_age
dim(n_ageoA314	)<- N_age
dim(n_ageoE324	)<- N_age
dim(n_ageoI324	)<- N_age
dim(n_ageoA324	)<- N_age
dim(n_ageoE412	)<- N_age
dim(n_ageoI412	)<- N_age
dim(n_ageoA412	)<- N_age
dim(n_ageoE413	)<- N_age
dim(n_ageoI413	)<- N_age
dim(n_ageoA413	)<- N_age
dim(n_ageoE423	)<- N_age
dim(n_ageoI423	)<- N_age
dim(n_ageoA423	)<- N_age
dim(n_ageoR123	)<- N_age
dim(n_ageoR124	)<- N_age
dim(n_ageoR134	)<- N_age
dim(n_ageoR234	)<- N_age
dim(n_ageoE1234	)<- N_age
dim(n_ageoI1234	)<- N_age
dim(n_ageoA1234	)<- N_age
dim(n_ageoE2134	)<- N_age
dim(n_ageoI2134	)<- N_age
dim(n_ageoA2134	)<- N_age
dim(n_ageoE3124	)<- N_age
dim(n_ageoI3124	)<- N_age
dim(n_ageoA3124	)<- N_age
dim(n_ageoE4123	)<- N_age
dim(n_ageoI4123	)<- N_age
dim(n_ageoA4123	)<- N_age
dim(n_ageoR1234	)<- N_age


#dim(cumu_inc) <- N_age
dim(seroprev) <- 7
dim(seroprev_num) <- 7
dim(seroprev_den) <- 7
dim(inc_day_gii4)<- 5
dim(inc_day_gii)<- 5
dim(inc_day_gi3)<- 5
dim(inc_day_gi)<- 5
dim(inc_day_all)<- 5
dim(inc_year_gii4)<- 5
dim(inc_year_gii)<- 5
dim(inc_year_gi3)<- 5
dim(inc_year_gi)<- 5
dim(pop_by4age)<- 5
#dim(n_risk) <- N_age
dim(infections_day_gii4)<- 4
dim(infections_day_gii)<- 4
dim(infections_day_gi3)<- 4
dim(infections_day_gi)<- 4
dim(reported_wk)<- N_age
dim(n_muM	)<- N_age
dim(n_muG	)<- N_age
dim(n_muS	)<- N_age
dim(n_muE1	)<- N_age
dim(n_muI1	)<- N_age
dim(n_muA1	)<- N_age
dim(n_muR1	)<- N_age
dim(n_muE2	)<- N_age
dim(n_muI2	)<- N_age
dim(n_muA2	)<- N_age
dim(n_muR2	)<- N_age
dim(n_muE3	)<- N_age
dim(n_muI3	)<- N_age
dim(n_muA3	)<- N_age
dim(n_muR3	)<- N_age
dim(n_muE4	)<- N_age
dim(n_muI4	)<- N_age
dim(n_muA4	)<- N_age
dim(n_muR4	)<- N_age
dim(n_muE12	)<- N_age
dim(n_muI12	)<- N_age
dim(n_muA12	)<- N_age
dim(n_muE13	)<- N_age
dim(n_muI13	)<- N_age
dim(n_muA13	)<- N_age
dim(n_muE14	)<- N_age
dim(n_muI14	)<- N_age
dim(n_muA14	)<- N_age
dim(n_muE21	)<- N_age
dim(n_muI21	)<- N_age
dim(n_muA21	)<- N_age
dim(n_muE23	)<- N_age
dim(n_muI23	)<- N_age
dim(n_muA23	)<- N_age
dim(n_muE24	)<- N_age
dim(n_muI24	)<- N_age
dim(n_muA24	)<- N_age
dim(n_muE31	)<- N_age
dim(n_muI31	)<- N_age
dim(n_muA31	)<- N_age
dim(n_muE32	)<- N_age
dim(n_muI32	)<- N_age
dim(n_muA32	)<- N_age
dim(n_muE34	)<- N_age
dim(n_muI34	)<- N_age
dim(n_muA34	)<- N_age
dim(n_muE41	)<- N_age
dim(n_muI41	)<- N_age
dim(n_muA41	)<- N_age
dim(n_muE42	)<- N_age
dim(n_muI42	)<- N_age
dim(n_muA42	)<- N_age
dim(n_muE43	)<- N_age
dim(n_muI43	)<- N_age
dim(n_muA43	)<- N_age
dim(n_muR12	)<- N_age
dim(n_muR13	)<- N_age
dim(n_muR14	)<- N_age
dim(n_muR23	)<- N_age
dim(n_muR24	)<- N_age
dim(n_muR34	)<- N_age
dim(n_muE123	)<- N_age
dim(n_muI123	)<- N_age
dim(n_muA123	)<- N_age
dim(n_muE124	)<- N_age
dim(n_muI124	)<- N_age
dim(n_muA124	)<- N_age
dim(n_muE134	)<- N_age
dim(n_muI134	)<- N_age
dim(n_muA134	)<- N_age
dim(n_muE213	)<- N_age
dim(n_muI213	)<- N_age
dim(n_muA213	)<- N_age
dim(n_muE214	)<- N_age
dim(n_muI214	)<- N_age
dim(n_muA214	)<- N_age
dim(n_muE234	)<- N_age
dim(n_muI234	)<- N_age
dim(n_muA234	)<- N_age
dim(n_muE312	)<- N_age
dim(n_muI312	)<- N_age
dim(n_muA312	)<- N_age
dim(n_muE314	)<- N_age
dim(n_muI314	)<- N_age
dim(n_muA314	)<- N_age
dim(n_muE324	)<- N_age
dim(n_muI324	)<- N_age
dim(n_muA324	)<- N_age
dim(n_muE412	)<- N_age
dim(n_muI412	)<- N_age
dim(n_muA412	)<- N_age
dim(n_muE413	)<- N_age
dim(n_muI413	)<- N_age
dim(n_muA413	)<- N_age
dim(n_muE423	)<- N_age
dim(n_muI423	)<- N_age
dim(n_muA423	)<- N_age
dim(n_muR123	)<- N_age
dim(n_muR124	)<- N_age
dim(n_muR134	)<- N_age
dim(n_muR234	)<- N_age
dim(n_muE1234	)<- N_age
dim(n_muI1234	)<- N_age
dim(n_muA1234	)<- N_age
dim(n_muE2134	)<- N_age
dim(n_muI2134	)<- N_age
dim(n_muA2134	)<- N_age
dim(n_muE3124	)<- N_age
dim(n_muI3124	)<- N_age
dim(n_muA3124	)<- N_age
dim(n_muE4123	)<- N_age
dim(n_muI4123	)<- N_age
dim(n_muA4123	)<- N_age
dim(n_muR1234	)<- N_age

dim(	n_MS	)<-  N_age
dim(	n_S_E1	)<-  N_age
dim(	n_S_E2	)<-  N_age
dim(	n_S_E3	)<-  N_age
dim(	n_S_E4	)<-  N_age


dim(	n_E1_I1	)<-N_age
dim(	n_E12_I12	)<-N_age
dim(	n_E123_I123	)<-N_age
dim(	n_E1234_I1234	)<-N_age
dim(	n_E124_I124	)<-N_age
dim(	n_E13_I13	)<-N_age
dim(	n_E134_I134	)<-N_age
dim(	n_E14_I14	)<-N_age
dim(	n_E2_I2	)<-N_age
dim(	n_E21_I21	)<-N_age
dim(	n_E213_I213	)<-N_age
dim(	n_E2134_I2134	)<-N_age
dim(	n_E214_I214	)<-N_age
dim(	n_E23_I23	)<-N_age
dim(	n_E234_I234	)<-N_age
dim(	n_E24_I24	)<-N_age
dim(	n_E3_I3	)<-N_age
dim(	n_E31_I31	)<-N_age
dim(	n_E312_I312	)<-N_age
dim(	n_E3124_I3124	)<-N_age
dim(	n_E314_I314	)<-N_age
dim(	n_E32_I32	)<-N_age
dim(	n_E324_I324	)<-N_age
dim(	n_E34_I34	)<-N_age
dim(	n_E4_I4	)<- N_age
dim(	n_E41_I41	)<-N_age
dim(	n_E412_I412	)<-N_age
dim(	n_E4123_I4123	)<-N_age
dim(	n_E413_I413	)<-N_age
dim(	n_E42_I42	)<-N_age
dim(	n_E423_I423	)<-N_age
dim(	n_E43_I43	)<-N_age
dim(	n_I1_A1	)<-N_age
dim(	n_I12_A12	)<-N_age
dim(	n_I123_A123	)<-N_age
dim(	n_I1234_A1234	)<-N_age
dim(	n_I124_A124	)<-N_age
dim(	n_I13_A13	)<-N_age
dim(	n_I134_A134	)<-N_age
dim(	n_I14_A14	)<-N_age
dim(	n_I2_A2	)<-N_age
dim(	n_I21_A21	)<-N_age
dim(	n_I213_A213	)<-N_age
dim(	n_I2134_A2134	)<-N_age
dim(	n_I214_A214	)<-N_age
dim(	n_I23_A23	)<-N_age
dim(	n_I234_A234	)<-N_age
dim(	n_I24_A24	)<-N_age
dim(	n_I3_A3	)<-N_age
dim(	n_I31_A31	)<-N_age
dim(	n_I312_A312	)<-N_age
dim(	n_I3124_A3124	)<-N_age
dim(	n_I314_A314	)<-N_age
dim(	n_I32_A32	)<-N_age
dim(	n_I324_A324	)<-N_age
dim(	n_I34_A34	)<-N_age
dim(	n_I4_A4	)<-N_age
dim(	n_I41_A41	)<-N_age
dim(	n_I412_A412	)<-N_age
dim(	n_I4123_A4123	)<-N_age
dim(	n_I413_A413	)<-N_age
dim(	n_I42_A42	)<-N_age
dim(	n_I423_A423	)<-N_age
dim(	n_I43_A43	)<-N_age
dim(	n_A1_R1	)<-N_age
dim(	n_A12_R12	)<-N_age
dim(	n_A123_R123	)<-N_age
dim(	n_A1234_R1234	)<-N_age
dim(	n_A124_R124	)<-N_age
dim(	n_A13_R13	)<-N_age
dim(	n_A134_R134	)<-N_age
dim(	n_A14_R14	)<-N_age
dim(	n_A2_R2	)<-N_age
dim(	n_A21_R12	)<-N_age
dim(	n_A213_R123	)<-N_age
dim(	n_A2134_R1234	)<-N_age
dim(	n_A214_R124	)<-N_age
dim(	n_A23_R23	)<-N_age
dim(	n_A234_R234	)<-N_age
dim(	n_A24_R24	)<-N_age
dim(	n_A3_R3	)<-N_age
dim(	n_A31_R13	)<-N_age
dim(	n_A312_R123	)<-N_age
dim(	n_A3124_R1234	)<-N_age
dim(	n_A314_R134	)<-N_age
dim(	n_A32_R23	)<-N_age
dim(	n_A324_R234	)<-N_age
dim(	n_A34_R34	)<-N_age
dim(	n_A4_R4	)<-N_age
dim(	n_A41_R14	)<-N_age
dim(	n_A412_R124	)<-N_age
dim(	n_A4123_R1234	)<-N_age
dim(	n_A413_R134	)<-N_age
dim(	n_A42_R24	)<-N_age
dim(	n_A423_R234	)<-N_age
dim(	n_A43_R34	)<-N_age



dim(	n_R1_A1	)<-  N_age
dim(	n_R1_S	)<-  N_age
dim(	n_R1_E31	)<-  N_age
dim(	n_R1_E41	)<-  N_age
dim(	n_R1_E21	)<-  N_age
dim(	n_R2_A2	)<-  N_age
dim(	n_R2_S	)<-  N_age
dim(	n_R2_E12	)<-  N_age
dim(	n_R2_E32	)<-  N_age
dim(	n_R2_E42	)<-  N_age
dim(	n_R3_A3	)<-  N_age
dim(	n_R3_S	)<-  N_age
dim(	n_R3_E13	)<-  N_age
dim(	n_R3_E23	)<-  N_age
dim(	n_R3_E43	)<-  N_age
dim(	n_R4_A4	)<-  N_age
dim(	n_R4_S	)<-  N_age
dim(	n_R4_E14	)<-  N_age
dim(	n_R4_E24	)<-  N_age
dim(	n_R4_E34	)<-  N_age
dim(	n_R13_A13	)<-  N_age
dim(	n_R13_A31	)<-  N_age
dim(	n_R13_wane	)<-  N_age
dim(	n_R13_E213	)<-  N_age
dim(	n_R13_E413	)<-  N_age
dim(	n_R12_A12	)<-  N_age
dim(	n_R12_A21	)<-  N_age
dim(	n_R12_wane	)<-  N_age
dim(	n_R12_E312	)<-  N_age
dim(	n_R12_E412	)<-  N_age
dim(	n_R23_A23	)<-  N_age
dim(	n_R23_A32	)<-  N_age
dim(	n_R23_wane	)<-  N_age
dim(	n_R23_E123	)<-  N_age
dim(	n_R23_E423	)<-  N_age
dim(	n_R14_A14	)<-  N_age
dim(	n_R14_A41	)<-  N_age
dim(	n_R14_wane	)<-  N_age
dim(	n_R14_E214	)<-  N_age
dim(	n_R14_E314	)<-  N_age
dim(	n_R24_A24	)<-  N_age
dim(	n_R24_A42	)<-  N_age
dim(	n_R24_wane	)<-  N_age
dim(	n_R24_E124	)<-  N_age
dim(	n_R24_E324	)<-  N_age
dim(	n_R34_A34	)<-  N_age
dim(	n_R34_A43	)<-  N_age
dim(	n_R34_wane	)<-  N_age
dim(	n_R34_E134	)<-  N_age
dim(	n_R34_E234	)<-  N_age
dim(	n_R123_A123	)<-  N_age
dim(	n_R123_A213	)<-  N_age
dim(	n_R123_A312	)<-  N_age
dim(	n_R123_wane	)<-  N_age
dim(	n_R123_E4123	)<-  N_age
dim(	n_R124_A124	)<-  N_age
dim(	n_R124_A214	)<-  N_age
dim(	n_R124_A412	)<-  N_age
dim(	n_R124_wane	)<-  N_age
dim(	n_R124_E3124	)<-  N_age
dim(	n_R134_A134	)<-  N_age
dim(	n_R134_A413	)<-  N_age
dim(	n_R134_A314	)<-  N_age
dim(	n_R134_wane	)<-  N_age
dim(	n_R134_E2134	)<-  N_age
dim(	n_R234_A234	)<-  N_age
dim(	n_R234_A423	)<-  N_age
dim(	n_R234_A324	)<-  N_age
dim(	n_R234_wane	)<-  N_age
dim(	n_R234_E1234	)<-  N_age
dim(	n_R1234_A1234	)<-  N_age
dim(	n_R1234_A2134	)<-  N_age
dim(	n_R1234_A3124	)<-  N_age
dim(	n_R1234_A4123	)<-  N_age
dim(	n_R1234_wane	)<-  N_age
dim(	n_R1vacc)<-  N_age
dim(	n_R2vacc)<-  N_age
dim(	n_R3vacc)<-  N_age
dim(	n_R4vacc)<-  N_age
dim(	n_R13vacc)<-  N_age
dim(	n_R12vacc)<-  N_age
dim(	n_R23vacc)<-  N_age
dim(	n_R24vacc)<-  N_age
dim(	n_R14vacc)<-  N_age
dim(	n_R34vacc)<-  N_age
dim(	n_R123vacc)<-  N_age
dim(	n_R124vacc)<-  N_age
dim(	n_R134vacc)<-  N_age
dim(	n_R234vacc)<-  N_age
dim(	n_R1234vacc)<-  N_age
dim(	n_toR13vacc)<-  N_age
dim(	n_toR123vacc)<-  N_age
dim(	n_toR134vacc)<-  N_age
dim(	n_toR1234vacc)<-  N_age

dim(n_Gvacc)<-N_age
dim(n_Mvacc)<-N_age
dim(n_Svacc)<-N_age

dim(n_bM) <- N_age
dim(n_bG) <- N_age
dim(n_bS) <- N_age
dim(p_birth)   <- 3
dim(p_mu) <- N_age
dim(p_aging) <- N_age
dim(p_MS) <- N_age
dim(p_SE1) <- N_age
dim(p_SE2) <- N_age
dim(p_R1A1) <- N_age
dim(p_R2A2) <- N_age
dim(p_R3A3) <- N_age
dim(p_R4A4) <- N_age
dim(p_SE3) <- N_age
dim(p_SE4) <- N_age
dim(p_R1_to_E2_gi) <- N_age
dim(p_R2_to_E1_gi) <- N_age
dim(p_R3_to_E4_gii) <- N_age
dim(p_R4_to_E3_gii) <- N_age
dim(p_vacc) <- N_age
dim(mu)<-N_age
dim(m) <- c(N_age, N_age)
dim(m_holi) <- c(N_age, N_age)
dim(cmx_1) <- c(N_age, N_age)
dim(cmx_2) <- c(N_age, N_age)
dim(cmx_3) <- c(N_age, N_age)
dim(cmx_4) <- c(N_age, N_age)
dim(cmx_5) <- c(N_age, N_age)
dim(cmx_6) <- c(N_age, N_age)
dim(cmx_7) <- c(N_age, N_age)
dim(cmx_8) <- c(N_age, N_age)
dim(cmx_9) <- c(N_age, N_age)
dim(c1_ij) <- c(N_age, N_age)
dim(c2_ij) <- c(N_age, N_age)
dim(c3_ij) <- c(N_age, N_age)
dim(c4_ij) <- c(N_age, N_age)
dim(contact_matrix) <- c(N_age, N_age)
dim(lambda_1) <- N_age
dim(lambda_2) <- N_age
dim(lambda_3) <- N_age
dim(lambda_4) <- N_age
dim(vac_camp_cov) <- N_age 
dim(vac_sche_cov) <- N_age