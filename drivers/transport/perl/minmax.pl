#!/usr/bin/perl

open(IN,"<$ARGV[0]");
@dump = (<IN>);
close(IN);

$n = $dump[0];
chomp($n);

$a = 0;
$rrminmax = 1.e10;

while ($a*($n+2) < $#dump) {
$rrmax=0.;
for ($i=0;$i<$n;$i++) {
for ($j=$i+1;$j<$n;$j++) {
@x = split/ +/,$dump[$i+$a*($n+2)+2];
@y = split/ +/,$dump[$j+$a*($n+2)+2];
$rr = sqrt(($x[1]-$y[1])**2+($x[2]-$y[2])**2+($x[3]-$y[3])**2);
if ($rr > $rrmax) {$rrmax = $rr};
}
}
if ($rrmax < $rrminmax) {$rrminmax = $rrmax; $win = $a};
$a++;
}

open(OUT,">smallest.geo");
#print OUT "$n\n";
#print OUT "Conformer # $win out of $a \n";
for ($i=0;$i<$n;$i++) {
print OUT $dump[$i+$win*($n+2)+2];
}
