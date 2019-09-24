 %mem=204MW                         
 %chk=tmp
 %NProcShared=          12
 # um062x/6-311+g(d,p) opt=(noeig,intern,maxcyc=50)                     
 # int=ultrafine nosym guess=(mix,always)                               

 geom            6

           0           1
               c1                                            
               c2  c1 rcc1                                   
               h3  c2 rcc2 c1 accc1                          
               h4  c2 rcc3 c1 ahcc1 h3 b1                    
               h41 c2 rcc4 c1 ahcc5 h3 b5                    
               h5  c1 rch1 c2 ahcc2 h3 b2                    
               h6  c1 rch2 c2 ahcc3 h5 b3                    
               h7  c1 rch3 c2 ahcc4 h5 b4                    

 RCC2                             1.0821000000000001     
 RCC3                             1.0821000000000001     
 RCC4                             1.0822000000000001     
 RCH1                             1.0822000000000001     
 RCH2                             1.0821000000000001     
 RCH3                             1.0821000000000001     
 ACCC1                            101.67890000000000     
 AHCC1                            101.72190000000001     
 AHCC2                            101.64210000000000     
 AHCC3                            101.67570000000001     
 AHCC4                            101.66980000000000     
 AHCC5                            101.65670000000000     
 B1                               120.09970000000000     
 B5                              -119.90600000000001     
 B3                               120.03060000000001     
 B4                              -119.93989999999999     

 RCC1                             2.3999999999999999     
 B2                               20.069550725641612     

