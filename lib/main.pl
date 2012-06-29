#!/usr/bin/perl

use InSilicoSpectro::InSilico::MassCalculator;
InSilicoSpectro::InSilico::MassCalculator::init('insilicodef.xml');
use MSMSTheoSpectrum;

# get fragType hash
my %fragType = %{InSilicoSpectro::InSilico::MassCalculator::fragType};
my %series= %{InSilicoSpectro::InSilico::MassCalculator::series};
my $peptide = 'QCTIPADFK';
my @modif   = ('');
my $modif   = '';
print "$peptide mass is ", getPeptideMass( pept => $peptide, modif => \@modif ),
  "\n";
getFragmentMasses(
    pept      => $peptide,
    modif     => $modif,
    fragTypes => [ 'b', 'y', 'b++', 'y++', 'immo' ],
    spectrum  => \%spectrum
);
my $theoSpectrum = new MSMSTheoSpectrum(
    theoSpectrum => \%spectrum,
    massType     => getMassType()
);
my $len = $theoSpectrum->getPeptideLength();
print "Fragments of ", $theoSpectrum->getPeptide(), " (",
  modifToString( $theoSpectrum->getModif(), $len ), ", ",
  $theoSpectrum->getPeptideMass(), " Da):\n";

foreach ( $theoSpectrum->getTermIons(\%fragType,\%series) ) {
    print $theoSpectrum->toString($_)."\n";;
}
print "\n";
foreach ( $theoSpectrum->getInternIons() ) {
    print $theoSpectrum->toString($_)."\n";
}

################## output ##########################
#Fragments of QCTIPADFK (::::::::::, 1021.490255 Da):
#y: 1:147.11335  2:294.18176     3:409.2087      4:480.24581     5:577.29857     6:690.38263     7:791.43031     8:894.4395      9:      (z=1, mono)
#y++: 1:74.060315        2:147.59452     3:205.10799     4:240.626545    5:289.152925    6:345.694955    7:396.218795    8:447.72339     9:      (z=2, mono)
#b: 1:   2:232.075595    3:333.123275    4:446.207335    5:543.260095    6:614.297205    7:729.324145    8:876.392555    9:      (z=1, mono)
#b++: 1: 2:116.5414375   3:167.0652775   4:223.6073075   5:272.1336875   6:307.6522425   7:365.1657125   8:438.6999175   9:      (z=2, mono)
#immo: F:120.08131       I:86.09696      P:70.06566      (z=1, mono)
####################################################
