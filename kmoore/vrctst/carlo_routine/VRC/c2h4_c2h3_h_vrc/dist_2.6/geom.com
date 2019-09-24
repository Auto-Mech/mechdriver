 %mem=204MW                         
 %chk=tmp
 %NProcShared=          12
 # uwb97xd/6-311+g(d,p)  opt=internal                                   
 # int=ultrafine nosym guess(mix,always)                                

 geom            0

           0           1
  c1                                                         
  c2 c1 rcc1                                                 
  h3 c2 rcc2 c1 accc1                                        
  h4 c2 rcc3 c1 ahcc1 h3 b1                                  
  h5 c1 rch1 c2 ahcc2 h3 b2                                  
  h6 c1 rch2 c2 ahcc3 h5 b3                                  

 RCC1                             1.3212999999999999     
 RCC2                             1.0825000000000000     
 RCC3                             1.0825000000000000     
 RCH2                             1.0825000000000000     
 ACCC1                            121.62090000000001     
 AHCC1                            121.62090000000001     
 AHCC2                            121.62170000000000     
 AHCC3                            121.62170000000000     
 B1                               179.99390000000000     
 B2                               180.00129999999999     
 B3                               180.00319999999999     

 RCH1                             1.8000000000000000     

