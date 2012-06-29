
=head1 NAME

InSilicoSpectro - Open source Perl library for proteomics

=head1 INSILICOSPECTRO PROJECT DESCRIPTION

This is the description of the entire InSilicoSpectro project; a description of the InSilicoSpectro.pm module
is provided hereunder.

InSilicoSpectro is a proteomics open-source project intended to cover common operations in
mass list file format conversions, protein sequence digestion, theoretical mass spectra computations, theoretical
and experimental MS data matching, text/graphic display, peptide retention time predictions, etc.

The problems of raw data processing, storage and database searching are not addressed by the
InSilicoSpectro project. InSilicoSpectro is released under the LGPL license and it is available from
a dedicated web site at http://insilicospectro.vital-it.ch.

The general design of the modules follows the object oriented programming (OOP) model and most of the
modules are class definitions actually. The module that implements most of the theoretical mass computation routines supports
a dual OOP and procedural programming model. InSilicoSpectro modules make use of some
Perl modules that are not part of the standard Perl distribution, such as Statistics:Regression, XML:Twig, GD, and
IA:NNFlex.

We have developed a simple and minimal hierarchy to represent protein sequences and peptides (as digestion product)
in a way that, on the one hand, fits the needs of the computations we perform and, on the other hand, stays relatively
neutral in its design. Thus it should be possible to combine the latter classes with existing
projects at users sites, e.g. via multiple inheritance, or to use them as the basis of more sophisticated objects.

InSilicoSpectro Perl code is documented mainly via pod and a wide collection of simple and focused examples. An
introductory explanation is provided here to guide new users and give them an understanding of the library that
should be sufficient such that pod and the examples are the only necessary documentation.

=head2 Installation

=head2 Library organization

InSilicoSpectro modules (lib/InSilicoSpectro) are organized according to their function. At the more general level there is
a module named InSilicoSpectro.pm (This one!!) that provides general functionalities for initializing all other modules.
More specialized modules are grouped in three folders:

=over 4

=item Spectra, for mass list-related;

=item InSilico, for computational modules;

=item Utils, for a few utility modules.

=back

In addition, illustrative examples can be found in three folders:

=over 4

=item scripts, which contains a set of tools implemented with InSilicoSpectro modules;

=item cgi, which contains scripts implementing a simple web-based set of tools;

=item t, which contains test programs that are examples as well.

=back

Now, by considering the main topics we cover in InSilicoSpectro one after another, we introduce the
main modules and examples the user should try and look at to gain autonomy with the whole library.

=head2 Mass list file format conversion

A general purpose conversion program, convertSpectra.pl in folder scripts, allows you to convert
one mass list format to another. A CGIzed version exists in the cgi folder: cgiConvertSpectra.pl.

convertSpectra.pl is a good starting point to see a high-level usage of the basic methods implemented
in the underlying modules.

InSilicoSpectro::Spectra::ExpSpectrum is the basic class for representing spectra, i.e. a list
of peaks (namely a list of pointers to peaks). Peaks are represented as list of attributes such
as mass, intensity, SN, etc. The order of the attributes in these lists is given by an object
of class InSilicoSpectro::Spectra::PeakDescriptor. See t/Spectra/testExpSpectrum.pl and
t/Spectra/testPeakDescriptor.pl.

By means of classes InSilicoSpectro::Spectra::MSSpectra, InSilicoSpectro::Spectra::MSMSSpectra,
InSilicoSpectro::Spectra::MSMSCmpd, and InSilicoSpectro::Spectra::MSRun we represent PMF (MS) and
MS/MS spectra, and HPLC runs. See t/Spectra/testSpectra.pl.

=head2 Utils

The module InSilicoSpectro::Utils::IO.pm contains miscellaneous utilities for accessing compressed files, defining a common
verbose variable, etc.

=head2 pI estimations

scripts/computePI.pl is a tool that exemplify the usage of the class InSilicoSpectro::InSilico::IsoelPoint.
Examples of how to use it can be found in t/InSilico/examples_rt_pi. See also the example in t/InSilico/testIsoelPoint.pl.
A CGI version of computePI.pl can be found in cgi folder.

=head2 Retention time prediction

scripts/computeRT.pl is a tool that exemplify the usage of the class InSilicoSpectro::InSilico::RetentionTimer.
Examples of how to use it can be found in t/InSilico/examples_rt_pi. See also the examples in t/InSilico/testPetritis.pl
and t/InSilico/testHodges.pl. A CGI version of computeRT.pl can be found in cgi folder.

=head2 Enzymes

Enzymes are modeled by class InSilicoSpectro::InSilico::CleavEnzyme. See t/InSilico/testCleavEnzyme.pl.

=head2 PTMs and other modifications

Modifications of residues are modeled by class InSilicoSpectro::InSilico::ModRes. See t/InSilico/testModRes.pl.

=head2 Protein and peptide sequences

The basic class for biological sequences is InSilicoSpectro::InSilico::Sequence. We then define
InSilicoSpectro::InSilico::AASequence to represent protein sequences with their modifications. A class
InSilicoSpectro::InSilico::Peptide is used for enzymatic digestion products as we need special data in
this case that are not part of a standard protein model.

Examples can be found in t/InSilico: testSequence.pl, testAASequence.pl, testPeptide.pl.

=head2 Protein digestion and mass computations

The main module for digestion and mass computations is InSilicoSpectro::InSilico::MassCalculator.
Examples of digestions and protein/peptide mass computations, including in the presence of fixed/variable
modifications, are found in t/InSilico: testCalcDigest.pl, testCalcDigestOOP.pl, and testCalcVarpept.pl. OOP
means an example with the OOP model as MassCalculator supports both an OOP and procedural interface.

=head2 PMF

The match between theoretical peptide masses and PMF experimental data is made by functions found
in InSilicoSpectro::InSilico::MassCalculator.
In the OOP model it is possible to represent PMF matches in objects of class InSilicoSpectro::InSilico::PMFMatch.
See t/InSilico/testCalcPMFMatch.pl and t/InSilico/testCalcPMFMatchOOP.pl.

=head2 Peptide fragmentation

Theoretical fragment masses are computed by functions found in InSilicoSpectro::InSilico::MassCalculator.
In the OOP model, theoretical MS/MS spectra can be represented as an object of class InSilicoSpectro::InSilico::MSMSTheoSpectrum,
which represents in turn the various ions as InSilicoSpectro::InSilico::InternIonSeries and InSilicoSpectro::InSilico::TermIonSeries.

The match between experimental and theoretical masses is also computed by InSilicoSpectro::InSilico::MassCalculator and in
the OOP model the class InSilicoSpectro::InSilico::MSMSTheoSpectrum can store the match in addition to the theoretical
spectrum.

See in t/InSilico: testCalcFrag.pl, testCalcFragOOP.pl, testCalcMatch.pl, testCalcMatchOOP.pl, getIonIntensities.pl, ionStat.R.

=head2 Graphical display of MS/MS spectra/matches

The class InSilicoSpectro::InSilico::MSMSOutput instanciates objects aimed at providing different formats
in order to represent MS/MS spectra and matches. See in t/InSilico: testMSMSOutText.pl,
testMSMSOutLatex.pl, testMSMSOutHtml.pl, testMSMSOutPlot.pl, testMSMSOutLegend.pl.

=head2 Mini web site

In folder miniweb we provide a perl script build-miniweb.pl that builds, from CGI scripts in folder cgi, a simple web site for protein digestion, mass computations,
and pI and retention time estimations.

=head1 MODULE DESCRIPTION

The module InSilicoSpectro.pm comprises generic functions that are useful for the whole project.

=head1 FUNCTIONS

=head3 saveInSilicoDef([$out])

Saves all registered definitions into the configuration file named $out, e.g. insilicodef.xml

=head3 getInSilicoDefFiles()

Returns the list of configuration files given by the operating system environment variable, whose name is stored in
$InSilicoSpectro::DEF_FILENAME_ENV (default "INSILICOSPECTRO_DEFFILE").

The environment variable can point more than one file (separated by ':'), or be a glob ('...*...' expression).

=head3 init([@files])

Loads a list of configuration files given as parameter or the default configuration files as returned by
getInSilicoDefFiles.

=head1 SEE ALSO

L<InSilicoSpectro::InSilico>, L<InSilicoSpectro::Spectra>, L<InSilicoSpectro::Utils>

=head1 COPYRIGHT

Copyright (C) 2004-2005  Geneva Bioinformatics (www.genebio.com) & Jacques Colinge (Upper Austria University of Applied Science at Hagenberg)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=head1 AUTHORS

Jacques Colinge, www.fhs-hagenberg.ac.at

Alexandre Masselot, www.genebio.com

=cut

use strict;

package InSilicoSpectro;
require Exporter;
use Carp;

use InSilicoSpectro::InSilico::CleavEnzyme;
use InSilicoSpectro::InSilico::ModRes;
use InSilicoSpectro::InSilico::MassCalculator;
our ( @ISA, @EXPORT, @EXPORT_OK, $VERSION );
@ISA = qw(Exporter);

@EXPORT =
  qw($VERSION &saveInSilicoDef &init &getInSilicoDefFiles $DEF_FILENAME_ENV);
@EXPORT_OK = ();
$VERSION   = "1.3.24";

our $DEF_FILENAME_ENV = 'INSILICOSPECTRO_DEFFILE';

sub saveInSilicoDef {
    my $out = shift;
    $out = ">$out" if ( ( defined $out ) and not $out =~ /^>/ );
    my $saver =
      ( defined $out )
      ? (
        new SelectSaver(
                 InSilicoSpectro::Utils::io->getFD($out)
              or CORE::die "cannot open [$out]: $!"
        )
      )
      : \*STDOUT;
    print <<EOT;
<inSilicoDefinitions>
  <elements/>
  <aminoAcids/>
  <codons/>
  <cleavEnzymes>
EOT
    foreach ( InSilicoSpectro::InSilico::CleavEnzyme::getList() ) {
        $_->getXMLTwigElt->print();
        print "\n";
    }
    print <<EOT;
  </cleavEnzymes>
  <fragTypeDescriptions/>
  <modRes>
EOT
    foreach ( InSilicoSpectro::InSilico::ModRes::getList() ) {
        $_->getXMLTwigElt->print();
        print "\n";
    }
    print <<EOT;
  </modRes>
</inSilicoDefinitions>
EOT
}

sub init {
    my @tmp = @_;
    push( @tmp, getInSilicoDefFiles() )
      if ( ( not @tmp ) and getInSilicoDefFiles() );

    unless (@tmp) {
        print STDERR
          "no default found, opening config file from Phenyx::Config::GlobalParam\n"
          if $InSilicoSpectro::Utils::io::VERBOSE;
        require Phenyx::Config::GlobalParam;
        Phenyx::Config::GlobalParam::readParam( undef, 1 );
        my $tmp;
        eval "
      use Phenyx::Manage::User;
      \$tmp=Phenyx::Manage::User->new(name=>'default')->getFile('insilicodef.xml');
    " or confess "error during eval statmeent: $!";
        push @tmp, $tmp;
    }
    @tmp or croak "must provide at least one  file argument";

    #warn  "init from files(@tmp)";
    InSilicoSpectro::InSilico::ModRes::init(@tmp);
    InSilicoSpectro::InSilico::CleavEnzyme::init(@tmp);
    InSilicoSpectro::InSilico::MassCalculator::init(@tmp);
}

use File::Basename;

sub getInSilicoDefFiles {
    my $env    = $ENV{$DEF_FILENAME_ENV};
    my $reldir = __FILE__;
    $reldir =~ s/\.pm//i;
    my @files;
    foreach ( split /(?<!\b[A-Za-z]):/, $env ) {
        if ( -f $_ ) {
            push @files, $_;
            next;
        }

        #add InSIlicoSpectro.pm dirname if file is not absolute
        $_ = "$reldir/$_" unless /^(\/|[A-Za-z]:\\)/;
        s/\\/\\\\/g;
        s/ /\\ /g;
        foreach ( glob "$_" ) {
            push @files, $_;
        }
    }
    return @files;
}

1;
