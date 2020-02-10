#!/usr/bin/perl

# usage: ./builddb.pl species.list
# The species list consists of SMILES strings obtained from RMG or elsewhere
# For each SMILES read from the species list, this script:
# 1) Uses openbabel to convert SMILES to canonical SMILES
# 2) If needed, creates directories of the type 
# 	./database/heavy-atom-list/heavy-atom-stoichiometry/stoichiometry/canonical-SMILES
# 3) Converts SMILES to a geometry guess using openbabel
# 4) [Missing: Deals with torsions/stereoisomers somehow.]
# 5) Executes "qcscript.pl QC-SCRIPT.dual" for the quantum chemistry runs specified in
#    QC-SCRIPT.qcscript. This script does the following:
# 	a) Looks for missing input/output QC files specified in QC-SCRIPT
# 	b) If output files are missing, QC jobs are sent to the queue
# 	c) If output files are available, they are read and summary files are written.

# allowed atom types & ordering
@atomtypes = ("C","O","N","Cl","H");

# flags
$canonicalsmiles = 1;     # convert to canonical (unique) smiles using obabel? 0 = no
$debug = 1;               # 0 = less output; 1 = more output

# summarize, could re-write to use these to control new jobs; right now these only control output
$summarize[$i++] = "1dmin/lj.dat";

# base filenames
$cwd = qx ! pwd !; 
chomp($cwd);
$cwdd = "$cwd/database";
#$cwdd = "$cwd/database2";
$cwdk = "$cwd/dictionary";
$cwdbin = "/lcrc/project/KTP/kmoore/auto1dmin/scripts";

# executables
$obabel = "$cwdbin/obabel";

open(IN,"<$ARGV[0]");  # Species list in RMG annotated format

foreach $line (<IN>) { ## MAIN LOOP ##
chomp($line);

if ($line =~ /END/) {last;}

$index++;

$line = " $line";
@dat = split/ +/,$line;
$name[$index]=$dat[1];
$smiles[$index]=$dat[3];
$multiplicity[$index]=$dat[4];

# SMILES
if ($canonicalsmiles) {
$tmp = qx ! $obabel -:"$smiles[$index]" -ocan 2> /dev/null !;  # convert to canonical SMILES
chomp ($tmp);
$smiles[$index] = $tmp;
}
$smiles[$index] =~ s/	//g;
$smiles[$index] =~ s/ //g;
print "Doing: $name[$index], canonical SMILES = $smiles[$index], multiplicity = $multiplicity[$index]\n";

# edit SMILES for directory names
$smilesdir[$index] = $smiles[$index];
$smilesdir[$index] =~ s/\(/p/g;
$smilesdir[$index] =~ s/\)/p/g;
$smilesdir[$index] =~ s/\//sf/g;
$smilesdir[$index] =~ s/\\/sb/g;
$smilesdir[$index] =~ s/#/h/g;
$smilesdir[$index] =~ s/@/a/g;
#$smilesdir[$index] =~ s/[/b/g;

# get stoichiometry using obabel's xyz
foreach $x (@atomtypes) { $atomcount{lc $x}=0 };
qx ! $obabel -:"$smiles[$index]" -O obabel-guess.xyz --gen3d 2> /dev/null !;
open(TMP,"<obabel-guess.xyz"); @dump = (<TMP>); close(TMP);
for ($i=2;$i<=$#dump;$i++) { 
	if ($dump[$i] =~ / *([a-zA-Z]+) /) { print($1); $atomcount{lc $1}++; } else { die("died at 1 in thermo.pl"); }
 print($dump[$i]);
}
foreach $x (@atomtypes) { 
  print($x);
  print("$atomcount{lc $x}\n");
	if ($atomcount{lc $x} > 0) { $stoichiometry[$index].=$x };
	if ($atomcount{lc $x} > 0) { $stoichiometry[$index].=$atomcount{lc $x} };
}
$natom[$index]=$#dump-1;

# FILE STRUCTURE
$level1 = $stoichiometry[$index];
$level1 =~ s/[0-9]+//g;
$level1 =~ s/H//;
if ($level1 eq "") {$level1="H"};

$level2 = $stoichiometry[$index];
$level2 =~ s/H[0-9]+//;
if ($level2 eq "") {$level2="H"};

$check = "$level1";
if (!(-e "$cwdd/$check")) {
if ($debug) {print "New heavy atom combination found. Adding $check.\n";}
mkdir "$cwdd/$check";
}

$check = "$level1/$level2";
if (!(-e "$cwdd/$check")) {
if ($debug) {print "New heavy atom stoichiometry found. Adding $check.\n";}
mkdir "$cwdd/$check";
}

$check = "$level1/$level2/$stoichiometry[$index]";
if (!(-e "$cwdd/$check")) {
if ($debug) {print "New heavy stoichiometry found. Adding $check.\n";}
mkdir "$cwdd/$check";
}

$check = "$level1/$level2/$stoichiometry[$index]/$smilesdir[$index]";
if (!(-e "$cwdd/$check")) {
if ($debug) {print "New SMILES found. Adding $check.\n";}
}
$path[$index] = $check;
if (!(-e "$cwdd/$path[$index]")) {mkdir "$cwdd/$path[$index]";}
if ($debug) {print "Database path: $path[$index]\n";}

# OBABEL GEOMETRY GUESSES
if (!(-e "$cwdd/xyz/$smilesdir[$index].xyz")) {
qx! cp -f "$cwdd/xyz/$smilesdir[$index].xyz" "./obabel-guess.xyz" !;
}
if (!(-e "$cwdd/$path[$index]/obabel-guess.xyz")) {
qx! mv -f "obabel-guess.xyz" "$cwdd/$path[$index]/obabel-guess.xyz" !;
} else {
qx! rm -f "obabel-guess.xyz" !;
}
if (!(-e "$cwdd/$path[$index]/obabel-guess.geo")) {
qx! tail -q -n +3 "$cwdd/$path[$index]/obabel-guess.xyz" > "$cwdd/$path[$index]/obabel-guess.geo" !;
}

# OBABEL CONFORMERS
if (!(-e "$cwdd/$path[$index]/obabel-confs.xyz")) {
qx! obabel $cwdd/$path[$index]/obabel-guess.xyz -O $cwdd/$path[$index]/obabel-confs.xyz --confab !;
}
open(TMP,"<$cwdd/$path[$index]/obabel-confs.xyz"); @dump = (<TMP>); close(TMP);
$nconf[$index] = ($#dump+1)/($natom[$index]+2);

if($nconf[$index] == 0) {
$nconf[$index]=1;
qx ! cp -f "$cwdd/$path[$index]/obabel-guess.xyz" "$cwdd/$path[$index]/obabel-confs.xyz" !;
open(TMP,"<$cwdd/$path[$index]/obabel-confs.xyz"); @dump = (<TMP>); close(TMP);

}

s/F/H/ for @dump;

print "Confomers: $nconf[$index]\n";

chdir "$cwdd/$path[$index]";
if (!(-e "1dmin")) { qx ! mkdir 1dmin ! };
chdir "$cwdd/$path[$index]/1dmin";
if (!(-e "lj.dat")) { 
  print "running $cwdbin/minmax.pl < ../obabel-confs.xyz \n";
  qx ! $cwdbin/minmax.pl ../obabel-confs.xyz !;
  $tmp = qx ! ls !;
  if ($tmp =~ /1DMIN/ ) {
    qx ! cat 1DMIN*/a*/output > output !;
    $tmp = qx ! $cwdk/lj.pl output !;
    open(OUT,">lj.dat");
    print OUT $tmp;
    close(OUT);
		} else {
			if (-e "q.x") {
				print " 1DMIN job in queue \n";
			} else {
				print " 1DMIN job lauched to queue \n";
				qx ! cp -R $cwdk/r0 .!;
				$spin=$multiplicity[$index]-1;  # MOLPRO SPIN = 2 * REAL SPIN
				open(INQCMOL,"<r0/qc.template");
				open(OUTQCMOL,">r0/qc.mol");
				foreach $tmp (<INQCMOL>) {
				$tmp =~ s/SPIN/$spin/;
				print OUTQCMOL $tmp;
				}
				close(INQCMOL);
				close(OUTQCMOL);
				qx ! cp $cwdk/batch-submit.sh .!;
				qx ! cp $cwdk/ingen.pl .!;
				qx ! cp smallest.geo r0/. !;
				qx ! ./batch-submit.sh !;
		}
	}
}
if (-e "lj.dat") { 
	open(LJIN,"<lj.dat");
	@ljdump = (<LJIN>);
	close(LJIN);
	@tmp = split/[ +|\n]+/,"@ljdump";
	$ep = $tmp[4];
	$si = $tmp[7];
	open(OUTLJ,">lj.table");
	print OUTLJ "$name[$index] $si $ep ! $smiles[$index] $multiplicity[$index] \n";
	close(OUTLJ);
}

print "\n";

} ## END MAIN LOOP ##

# OUTPUT
print "OUTPUT:\n";
for ($i=1;$i<=$index;$i++) {
print "Name: $name[$i]\nMultiplicity: $multiplicity[$i]\nSmiles: $smiles[$i]\nPath: $path[$i]\n";
foreach $k (@summarize) {
print "Output for $k\n";
$tmppath="$cwdd/$path[$i]/$k";
print " PATH: $tmppath \n";
if (-e "$tmppath") { $tmp = qx ! cat $tmppath !; print "$tmp"; } else { print " MISSING\n";}
}
print "\n";
}

# TABLE
for ($i=1;$i<=$index;$i++) {
$tmppath="$cwdd/$path[$i]/1dmin/lj.table";
if (-e $tmppath) { $tmp = qx ! cat $tmppath !; $savetable.=$tmp; };
}
print $savetable;

