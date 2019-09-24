 %mem=204MW                         
 %chk=tmp
 %NProcShared=          12
 # uwb97xd/6-311+g(d,p)  opt=internal                                   
 # int=ultrafine nosym guess(mix,always)                                

 geom            0

           0           1
  c1                                                         
  h2 c1 rch1                                                 
  h3 c1 rch2 h2 ahcc1                                        
  h4 c1 rch3 h2 ahcc2 h3 b1                                  
  h5 c1 rch4 h2 ahcc3 h3 b2                                  

 RCH1                             1.0870000000000000     
 RCH2                             1.0869000000000000     
 RCH3                             1.0870000000000000     
 AHCC1                            109.46850000000001     
 AHCC2                            109.45970000000000     
 AHCC3                            109.45970000000000     
 B1                               120.00170000000000     
 B2                              -120.00170000000000     

 RCH4                             1.8000000000000000     

