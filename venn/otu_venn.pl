#!/usr/bin/perl
use strict;
use warnings;
use SVG;
use Getopt::Long;
use Math::Trig;

my %opts;
GetOptions (\%opts,"otu_tab=s","group=s","prefix=s","circle_r=i","ellipse_la=i","ellipse_ma=i","opacity=f","font_size=i");

my $usage = <<"USAGE";
	Program : $0
	Discription: Draw venn for a large number of samples
	Contact: 
	Usage:perl $0 [options]
		-otu_tab    <file>	OTU table file
		-group      <file>	sample group file
		-prefix     <str>	The prefix of output file, dafult otu_venn
		-circle_r   <num>	The radius of center circle, circle_r <50, dafult 20
		-ellipse_la <num>	The long axis of the ellipse, dafult 360
		-ellipse_ma <num>	The minor axis of the ellipse, dafult 80
		-opacity    <num>	The opacity of the ellipse's color,0=<opacity=<1,  dafult 0.5
		-font_size  <num>	The labels' font fize, dafult 21
	Example: perl $0 -otu_tab otu.xls -group group.txt -prefix venn -opacity 0.8
USAGE

die $usage if ( ! $opts{otu_tab} || ! $opts{group});
$opts{prefix} ||= "otu_venn";
$opts{circle_r} ||= 20;
$opts{ellipse_la} ||= 360;
$opts{ellipse_ma} ||= 80;
$opts{opacity} ||= 0.5;
$opts{font_size} ||= 20;

open G,"<$opts{group}" || die "Error: don't open group file: $opts{group}!\n";
my %group;my $sample_num=0; my %color;
while(<G>){
	chomp;
	next if(/^#/);
	my @l=split /\t/;
	$group{$l[0]}=$l[1];
	$color{$l[0]}=$l[2];
	$sample_num++;
}
close G;

open F,"<$opts{otu_tab}" || die "Error: don't open infile:$opts{otu_tab}!\n";
my %unique;my $common_num=0;
<F>;my $line2=<F>; chomp $line2; my @sample=split /\t/,$line2;
while(<F>){
	chomp;
	my @l=split /\t/;
	my $num;my %hash;
	for(my $i=1;$i<@l;$i++){
		if($l[$i] != 0){
			$num++;
			$hash{$sample[$i]}=$l[$i];
		}
	}
	$common_num++ if($num==$sample_num);
	if($num==1){
		(my $s)=keys %hash;
		$unique{$s}++;
	}
}
close F;

my $svg = SVG->new(width=>$opts{ellipse_la}*2+500, height=>$opts{ellipse_la}*2+500);
my $x0 = $opts{ellipse_la}+250; my $y0 = $opts{ellipse_la}+250;
my $rx = $opts{ellipse_la}/2; my $ry = $opts{ellipse_ma}/2;

my @labels=sort{$group{$a} cmp $group{$b}} keys %group;
my $n=$sample_num;
for(my $i=0;$i<$n;$i++){
	my ($cx,$cy,$angle,$text_x,$text_y);
	my $offset=35;my $sample_length = length($labels[$i]);
	if(360*$i/$n>=0 && 360*$i/$n<90){
		$cx = $x0+($rx-$opts{circle_r})*sin(6.28*$i/$n);
		$cy = $y0-($rx-$opts{circle_r})*cos(6.28*$i/$n);
		$angle = 360*$i/$n-90;
		$svg->ellipse(cx => $cx, cy => $cy, rx => $rx, ry => $ry, transform => "rotate($angle $cx $cy)", style=>{stroke=>'black','stroke-width',0,fill=>"$color{$labels[$i]}",'fill-opacity'=> "$opts{opacity}"});
		$text_x = $x0+($opts{ellipse_la}-3*$opts{circle_r})*sin(6.28*$i/$n);
		$text_y = $y0-($opts{ellipse_la}-3*$opts{circle_r})*cos(6.28*$i/$n);
		$offset=35*$sample_length/9 if((90-360*$i/$n)<=20);
		$svg->text(x => $text_x, y => $text_y, 'font-size'=>$opts{font_size}, 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', $unique{$labels[$i]});
		$svg->text(x => $x0+($opts{ellipse_la}+$offset)*sin(6.28*$i/$n), y => $y0-($opts{ellipse_la}+$offset)*cos(6.28*$i/$n), 'font-size'=>"$opts{font_size}", 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', "$labels[$i]");
	}elsif(360*$i/$n>=90 && 360*$i/$n<180){
		$cx = $x0+($rx-$opts{circle_r})*sin(3.14-6.28*$i/$n);
		$cy = $y0+($rx-$opts{circle_r})*cos(3.14-6.28*$i/$n);
		$angle = 90-(180-360*$i/$n);
		$svg->ellipse(cx => $cx, cy => $cy, rx => $rx, ry => $ry, transform => "rotate($angle $cx $cy)", style=>{stroke=>'black','stroke-width',0,fill=>"$color{$labels[$i]}",'fill-opacity'=> "$opts{opacity}"});
		$text_x = $x0+($opts{ellipse_la}-3*$opts{circle_r})*sin(3.14-6.28*$i/$n);
		$text_y = $y0+($opts{ellipse_la}-3*$opts{circle_r})*cos(3.14-6.28*$i/$n);
		$offset=35*$sample_length/9 if((360*$i/$n-90)<=20);
		$svg->text(x => $text_x, y => $text_y, 'font-size'=>$opts{font_size}, 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', $unique{$labels[$i]});
		$svg->text(x => $x0+($opts{ellipse_la}+$offset)*sin(3.14-6.28*$i/$n), y => $y0+($opts{ellipse_la}+$offset)*cos(3.14-6.28*$i/$n), 'font-size'=>"$opts{font_size}", 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', "$labels[$i]");
	}elsif(360*$i/$n>=180 && 360*$i/$n<270){
		$cx = $x0-($rx-$opts{circle_r})*sin(6.28*$i/$n-3.14);
		$cy = $y0+($rx-$opts{circle_r})*cos(6.28*$i/$n-3.14);
		$angle = 360*$i/$n-180-90;
		$svg->ellipse(cx => $cx, cy => $cy, rx => $rx, ry => $ry, transform => "rotate($angle $cx $cy)", style=>{stroke=>'black','stroke-width',0,fill=>"$color{$labels[$i]}",'fill-opacity'=> "$opts{opacity}"});
		$text_x = $x0-($opts{ellipse_la}-3*$opts{circle_r})*sin(6.28*$i/$n-3.14);
		$text_y = $y0+($opts{ellipse_la}-3*$opts{circle_r})*cos(6.28*$i/$n-3.14);
		$offset=35*$sample_length/9 if((270-360*$i/$n)<=20);
		$svg->text(x => $text_x, y => $text_y, 'font-size'=>$opts{font_size}, 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', $unique{$labels[$i]});
		$svg->text(x => $x0-($opts{ellipse_la}+$offset)*sin(6.28*$i/$n-3.14), y => $y0+($opts{ellipse_la}+$offset)*cos(6.28*$i/$n-3.14), 'font-size'=>"$opts{font_size}", 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', "$labels[$i]");
	}elsif(360*$i/$n>=270 && 360*$i/$n<360){
		$cx = $x0-($rx-$opts{circle_r})*sin(6.28-6.28*$i/$n);
		$cy = $y0-($rx-$opts{circle_r})*cos(6.28-6.28*$i/$n);
		$angle = 90-(360-360*$i/$n);
		$svg->ellipse(cx => $cx, cy => $cy, rx => $rx, ry => $ry, transform => "rotate($angle $cx $cy)", style=>{stroke=>'black','stroke-width',0,fill=>"$color{$labels[$i]}",'fill-opacity'=> "$opts{opacity}"});
		$text_x = $x0-($opts{ellipse_la}-3*$opts{circle_r})*sin(6.28-6.28*$i/$n);
		$text_y = $y0-($opts{ellipse_la}-3*$opts{circle_r})*cos(6.28-6.28*$i/$n);
		$offset=35*$sample_length/9 if((360*$i/$n-270)<=20);
		$svg->text(x => $text_x, y => $text_y, 'font-size'=>$opts{font_size}, 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', $unique{$labels[$i]});
		$svg->text(x => $x0-($opts{ellipse_la}+$offset)*sin(6.28-6.28*$i/$n), y => $y0-($opts{ellipse_la}+$offset)*cos(6.28-6.28*$i/$n), 'font-size'=>"$opts{font_size}", 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.5, '-cdata', "$labels[$i]");
	}
}
$svg->circle(cx => $x0, cy => $y0, r => $opts{circle_r}, style=>{stroke=>'black','stroke-width',0.5,fill=>"red",'fill-opacity'=> 0.92});
$svg->text(x => $x0, y => $y0+($opts{font_size}-3)*6/15, 'font-size'=>$opts{font_size}-3, 'text-anchor'=>'middle', 'stroke', 'black', 'stroke-width',0.3, '-cdata', "$common_num");

my $out = $svg->xmlify;
open SVGFILE,">$opts{prefix}_venn.svg" || die $!;
print SVGFILE $out;
close SVGFILE;
`/usr/bin/convert $opts{prefix}_venn.svg $opts{prefix}_venn.png`;
