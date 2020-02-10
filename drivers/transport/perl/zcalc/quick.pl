#!/usr/bin/perl

if($#ARGV < 3) {die("usage: quick.pl epCALC sigCALC epTAB sigTAB cm-1 and A")}

$ec=abs($ARGV[0]);  # calc
$sc=$ARGV[1];
$et=abs($ARGV[2]);  # tab
$st=$ARGV[3];

print " CALC: $ec cm-1 $sc A \n";
print " TAB : $et cm-1 $st A \n\n";

$s2=($sc/$st)**2;

#@tlist=(300,500,1000,1500,2000,2500,3000);
@tlist=(300,1000,2000);

print " T,K  T*(CALC)  T*(TAB) \n";
foreach $temp (@tlist) {
$tc=$temp*0.695016/$ec;
$tt=$temp*0.695016/$et;
print "$temp $tc $tt";
if ($tc < .1 || $tt < .1) {print " LOW ";}
if ($tc > 200 || $tt > 200) {print " HIGH ";}
print "\n";
}
print "\n";

print "CALC/TAB \n";
print "method    sig**2 omega(T1) z(T1) omega(T2) z(T2) etc... \n";

# 9-6
print " LJ9-6  ";
printf("% 10.4f",$s2);
foreach $temp (@tlist) {
$tc=$temp*0.695016/$ec;
$tt=$temp*0.695016/$et;
$ts=log($tc);
$tmp=-3.775786e-04*$ts**6+4.824077e-03*$ts**5-1.582658e-02*$ts**4-3.434983e-02*$ts**3+3.221786e-01*$ts**2-8.249583e-01*$ts+1.499615;
if ($ts > 5.2) { $tmp = -8.757856E-02*$ts+8.511083E-01}
$ts=log($tt);
$tmp2=-3.775786e-04*$ts**6+4.824077e-03*$ts**5-1.582658e-02*$ts**4-3.434983e-02*$ts**3+3.221786e-01*$ts**2-8.249583e-01*$ts+1.499615;
if ($ts > 5.2) { $tmp2 = -8.757856E-02*$ts+8.511083E-01}
$o=$tmp/$tmp2;
$z=$o*$s2;
printf("% 10.4f",$o);
printf("% 10.4f",$z);
}
print " \n";

# 12-6
print " LJ12-6 ";
printf("% 10.4f",$s2);
foreach $temp (@tlist) {
$tc=$temp*0.695016/$ec;
$tt=$temp*0.695016/$et;
$ts=log($tc);
$tmp=-2.375644e-04*$ts**6+3.285312e-03*$ts**5-1.230859e-02*$ts**4-1.979413e-02*$ts**3+2.444541e-01*$ts**2-6.922991e-01*$ts+1.450545;
if ($ts > 5.2) { $tmp = -8.085029E-02*$ts+8.899438E-01}
$ts=log($tt);
$tmp2=-2.375644e-04*$ts**6+3.285312e-03*$ts**5-1.230859e-02*$ts**4-1.979413e-02*$ts**3+2.444541e-01*$ts**2-6.922991e-01*$ts+1.450545;
if ($ts > 5.2) { $tmp2 = -8.085029E-02*$ts+8.899438E-01}
$o=$tmp/$tmp2;
$z=$o*$s2;
printf("% 10.4f",$o);
printf("% 10.4f",$z);
}
print " \n";

# Troe's 12-6 (2,2)* fit
print " Troe   ";
printf("% 10.4f",$s2);
foreach $temp (@tlist) {
$tc=$temp*0.695016/$ec;
$tt=$temp*0.695016/$et;
$ts=log($tc)/log(10.);
$tmp = 1./(0.7+0.52*$ts);
$ts=log($tt)/log(10.);
$tmp2 = 1./(0.7+0.52*$ts);
$o=$tmp/$tmp2;
$z=$o*$s2;
printf("% 10.4f",$o);
printf("% 10.4f",$z);
}
print " \n\n";

print "CALC 9-6/TAB 12-6 \n";
print "method    sig**2 omega(T1) z(T1) omega(T2) z(T2) etc... \n";
print " 9/12   ";
$s2=1.;
$s2=($sc/$st)**2;
printf("% 10.4f",$s2);
foreach $temp (@tlist) {
$tc=$temp*0.695016/$ec;
$ts=log($tc);
$tmp=-3.775786e-04*$ts**6+4.824077e-03*$ts**5-1.582658e-02*$ts**4-3.434983e-02*$ts**3+3.221786e-01*$ts**2-8.249583e-01*$ts+1.499615;
if ($ts > 5.2) { $tmp = -8.757856E-02*$ts+8.511083E-01}
$tt=$temp*0.695016/$et;
$ts=log($tt);
$tmp2=-2.375644e-04*$ts**6+3.285312e-03*$ts**5-1.230859e-02*$ts**4-1.979413e-02*$ts**3+2.444541e-01*$ts**2-6.922991e-01*$ts+1.450545;
if ($ts > 5.2) { $tmp2 = -8.085029E-02*$ts+8.899438E-01}
$o=$tmp/$tmp2;
$z=$o*$s2;
printf("% 10.4f",$o);
printf("% 10.4f",$z);
}
print " \n";



