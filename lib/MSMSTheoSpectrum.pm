package MSMSTheoSpectrum;

require Exporter;

our ( @ISA, @EXPORT, @EXPORT_OK );

@ISA       = qw(Exporter);
@EXPORT    = qw();
@EXPORT_OK = ();

use strict;
use Carp;

our %visibleAttr = ( theoSpectrum => 1, massType => 1, tol => 1, minTol => 1 );
our %series;
my $massType;
our %immoAA;
my $immoDelta;
my $nTerm;
our %loss;
our %elMass;
my $cTerm;
my $invalidElementCall;

sub new {
    my $pkg = shift;

    my $tsp = {};
    my $class = ref($pkg) || $pkg;
    bless( $tsp, $class );
    my %h = @_;
    foreach ( keys(%h) ) {
        $tsp->$_( $h{$_} ) if ( $visibleAttr{$_} );
    }
    return $tsp;

}    # new

sub theoSpectrum {
    my $this = shift;

    if ( defined( $_[0] ) ) {
        if ( ref( $_[0] ) eq 'HASH' ) {
            $this->{theoSpectrum} = $_[0];
            my $peptide = $this->{theoSpectrum}{peptide};
            $this->{peptideLength} = length($peptide);
        }
        else {
            croak("Illegal data type for theoSpectrum");
        }
    }
    return $this->{theoSpectrum};

}    # theoSpectrum

sub massType {
    my ( $this, $mt ) = @_;

    if ( defined($mt) ) {
        $mt = int($mt);
        if ( ( $mt == 0 ) || ( $mt == 1 ) ) {
            $this->{massType} = $mt;
        }
        else {
            croak("Invalid massType value [$mt]");
        }
    }
    return $this->{massType};

}    # massType

=head2 tol([$mt])

Accessor/modifier of the attribute tol.

=cut

sub tol {
    my ( $this, $mt ) = @_;

    if ( defined($mt) ) {
        $this->{tol} = $mt;
    }
    return $this->{tol};

}    # tol

=head2 minTol([$mt])

Accessor/modifier of the attribute minTol.

=cut

sub minTol {
    my ( $this, $mt ) = @_;

    if ( defined($mt) ) {
        $this->{minTol} = $mt;
    }
    return $this->{minTol};

}    # minTol

=head2 getPeptide

Returns the peptide sequence or the Peptide object stored in the data structure produced
by getFragmentMasses.

=cut

sub getPeptide {
    my $this = shift;
    return $this->theoSpectrum()->{peptide};

}    # getPeptide

=head2 getPeptideLength

Returns the length of the fragmented peptide.

=cut

sub getPeptideLength {
    my $this = shift;
    return $this->{peptideLength};

}    # getPeptide Length

=head2 getModif

Returns a reference to the vector of modifications used for the theoretical masses
computation.

=cut

sub getModif {
    my $this = shift;
    return $this->theoSpectrum()->{modif};

}    # getModif

=head2 getPeptideMass

Returns the peptide mass.

=cut

sub getPeptideMass {
    my $this = shift;
    return $this->theoSpectrum()->{peptideMass};

}    # getPeptideMass

=head2 getFragType($name)

Returns a vector (series, charge, loss 1, repeat 1, loss 2, repeat 2, ...)
containing the parameters of a fragment type, where

=over 4

=item $name

is the name of the fragment type.

=item series

is the series on which this fragment type is based.

=item charge

is the charge state of a fragment of this type.

=item loss k

in case the fragment type includes losses, loss k is set to the name of
each loss possible for this fragment type.

=item repeat k

the maximum number of each loss, -1 means no maximum.

=back

=cut

sub getFragType {
    my $this     = shift;
    my $name     = shift;
    my $hashref  = shift;
    my %fragType = %{$hashref};

    my @tmp = ( $fragType{$name}{series}, $fragType{$name}{charge} );
    if ( defined( $fragType{$name}{loss} ) ) {
        for ( my $i = 0 ; $i < @{ $fragType{$name}{loss} } ; $i++ ) {
            push( @tmp,
                $fragType{$name}{loss}[$i],
                $fragType{$name}{repeat}[$i] );
        }
    }
    return @tmp;

}    # getFragType

=head2 getSeries($name)

Returns a vector (terminus, monoisotopic delta, average delta,
first fragment, last fragment) that contains the parameters
of series $name, where

=over 4

=item $name

is th name of the series, e.g. b or y.

=item terminus

is equal to 'N' or 'C'.

=item monoisotopic delta

is the mass delta to apply when computing the monoisotopic mass.
It is 0 for a b series and 1.007825 for y.

=item average delta

is the mass delta to apply when computing the average mass.
It is 0 for a b series and 1.007976 for y.

=item first fragment

is the number of the first fragment to compute. For instance,
fragment b1 is generally not detected hence first fragment
should be set to 2 for a b series. The rule is the same for
N- and C-term series.

=item last fragment

is the number of the last fragment counted from the end. For
instance, the last b fragment is normally not detected hence
last fragment should be set to 2. If a fragment containing
all the amino acids of the peptide is possible, then it
should be set to 1. The rule is the same for N- and C-term
series.

=back

=cut

sub getSeries {
    my $this    = shift;
    my $name    = shift;
    my $hashref = shift;
    my %series  = %{$hashref};

    return (
        $series{$name}{terminus}, $series{$name}{delta}[0],
        $series{$name}{delta}[1], $series{$name}{firstFrag},
        $series{$name}{lastFrag}
    );

}    # getSeries

=head2 cmpFragTypes

This function can be used in a sort of fragment type names. Fragment
type names are assumed to follow the rule:

=over 4

=item internal fragments

They are named after their generic name, only immonium ions
are supported so far and they are named 'immo'.

=item N-/C-terminal fragments

They must comply with the pattern

  ion&charge - loss1 -loss2 - ...

For instance, singly charged b ions are simply named 'b' and
their doubly and triply counterparts are names 'b++' and
'b+++'. This is the ion&charge part of the pattern above.

The losses may occur once or several times, multiple losses
are indicated in parentheses preceeded by multiplicity.
Examples are:

  b-H2O
  b-3(H2O)
  b++-H2O-NH3
  b++-3(H2O)-NH3
  y-H2O-2(H3PO4)-NH3

=back

The order on fragment type names is defined as follows: (1) immonium
ions always come after N-/C-terminal fragments; (2) N-/C-terminal
fragment types are compared by doing a sequence of comparisons
which continues as long as the compared values are equal. The first
comparison is on the ion type (a,b,y,...) followed by a comparison
on the charge. If ion types and charges are equal, comparisons
are made on the losses. The fragment that has less loss types is
considered smaller. If the two fragment types have the same number
of loss types then the losses are sorted lexicographically and the
first ones are compared on their name, if the names are the same then
the comparison is on the multiplicity, if the multiplicities are
the same then the second losses are compared, etc.

Asterisks that are used for signaling multiple possible losses are
ignored in the comparisons.

Since this function is defined in package MSMSOutput and it is used in
other packages with function sort (and predefined variables $a and $b),
we had to use prototypes ($$). Therefore it can no longer be exported
by the package MSMSOutput and you have to call it via MSMSOutput::cmpFragTypes.

Example:

foreach (sort MSMSOutput::cmpFragTypes ('y','b','y++','a','b-NH3','b-2(NH3)','b++-10(NH3)','b-H2O-NH3','immo(Y)', 'b++','y-NH3*','y-H2O*','z')){
  print $_,"\n";
}

=cut

=head2 getTermIons

Returns a vector containing all the terminal ion series sorted with respect to their name
by the function InSilicoSpectro::InSilico::MSMSOutput::cmpFragTypes.

=cut

sub ionType {
    my ( $this, $hashref, $it ) = @_;

    if ( defined($it) ) {
        $hashref->{ionType} = $it;
    }
    return $hashref->{ionType};
}

sub getTermIons {
    my $this        = shift;
    my %fragType    = %{ $_[0] };
    my %series_hash = %{ $_[1] };

    my $tol = $this->tol();
    my $minTol = $this->minTol() || 0.2;
    my @termIons;
    my $len       = length( $this->getPeptide() );
    my $names     = [ ( 1 .. $len ) ];
    my %spectrum  = %{ $this->{theoSpectrum} };
    my $massIndex = $spectrum{intensityIndex};
    foreach my $frag ( keys( %{ $spectrum{mass}{term} } ) ) {
        for ( my $i = 0 ; $i < @{ $spectrum{ionType}{$frag} } ; $i++ ) {
            my ( $series, $charge ) =
              ( $this->getFragType( $frag, \%fragType ) )[ 0, 1 ];
            my $terminus = ( $this->getSeries( $series, \%series_hash ) )[0];
            my $ionType = $spectrum{ionType}{$frag}[$i];
            my ( @masses, @matches );
            for ( my $j = $i * $len ; $j < ( $i + 1 ) * $len ; $j++ ) {
                push( @masses, $spectrum{mass}{term}{$frag}[$j] );
                if ( defined( $spectrum{match}{term}{$frag}[$j] ) ) {
                    if ( defined($tol) ) {
                        my $theoMass = $spectrum{mass}{term}{$frag}[$j];
                        my $expMass =
                          $spectrum{match}{term}{$frag}[$j][$massIndex];
                        my $error = $expMass - $theoMass;
                        if (
                            (
                                abs($error) / ( $theoMass + $expMass ) * 2e+6 <=
                                $tol
                            )
                            || ( abs($error) <= $minTol )
                          )
                        {
                            push( @matches, $spectrum{match}{term}{$frag}[$j] );
                        }
                        else {
                            push( @matches, undef );
                        }
                    }
                    else {
                        push( @matches, $spectrum{match}{term}{$frag}[$j] );
                    }
                }
                else {
                    push( @matches, undef );
                }
            }
            push(
                @termIons,
                {
                    ionType        => $ionType,
                    charge         => $charge,
                    series         => $series,
                    terminus       => $terminus,
                    names          => $names,
                    masses         => [@masses],
                    massType       => $this->massType(),
                    intensityIndex => $spectrum{intensityIndex},
                    massIndex      => $massIndex,
                    matches =>
                      ( defined( $spectrum{match} ) ? [@matches] : undef )
                }
            );
        }
    }
#### changed by xusheng;
    @termIons =
      sort { cmpFragTypes( $this->ionType($a), $this->ionType($b) ) } @termIons;
    return @termIons;
}    # getTermIons

sub toString {
    my $this    = shift;
    my $hashref = shift;
    my $string  = $this->ionType($hashref) . ": ";
    my @names   = @{ $this->names($hashref) };
    my @masses  = @{ $this->masses($hashref) };
    my @matches =
      defined( $this->matches($hashref) ) ? @{ $this->matches($hashref) } : ();
    for ( my $i = 0 ; $i < @masses ; $i++ ) {
        if ( defined( $names[$i] ) ) {
            $string .= $names[$i] . ':' . $masses[$i];
            if ( defined( $matches[$i] ) ) {
                $string .= "-match($matches[$i][0])";
            }
            $string .= "\t";
        }
    }
    $string .= '(z='
      . $this->charge($hashref) . ', '
      . ( $this->massType() == 0 ? 'mono' : 'avg' ) . ')';
    return $string;

}    # toString

sub matches {
    my ( $this, $hashref, $mp ) = @_;

    if ( defined($mp) ) {
        if ( ref($mp) eq 'ARRAY' ) {
            $hashref->{matches} = $mp;
        }
        else {
            croak("Illegal type for matches [$mp]");
        }
    }
    return $hashref->{matches};

}    # matches

sub charge {
    my ( $this, $hashref, $z ) = @_;

    if ( defined($z) ) {
        $z = int($z);
        if ( $z > 0 ) {
            $hashref->{charge} = $z;
        }
        else {
            croak("Invalid charge value [$z]");
        }
    }
    return $hashref->{charge};

}    # charge

sub masses {
    my ( $this, $hashref, $m ) = @_;

    if ( defined($m) ) {
        if ( ref($m) eq 'ARRAY' ) {
            $hashref->{masses} = $m;
        }
        else {
            croak("Invalid masses type [$m]");
        }
    }
    return $hashref->{masses};

}    # masses

sub names {
    my ( $this, $hashref, $n ) = @_;

    if ( defined($n) ) {
        if ( ref($n) eq 'ARRAY' ) {
            $hashref->{names} = $n;
        }
        else {
            croak("Invalid names type [$n]");
        }
    }
    return $hashref->{names};

}    # names

sub cmpFragTypes ($$) {
    my ( $fragA, $fragB ) = @_;

    # Presence of an immonium ion
    my $immoA = index( $fragA, 'immo' ) != -1;
    my $immoB = index( $fragB, 'immo' ) != -1;
    if ( $immoA && $immoB ) {
        return 0;
    }
    elsif ( $immoA && !$immoB ) {
        return 1;
    }
    elsif ( !$immoA && $immoB ) {
        return -1;
    }

    # Only N-/C-terminal fragments

    # Extracts ion types and charge, compare them
    my ( $ionA, @partA ) = split( /\-/, $fragA );
    my ( $ionB, @partB ) = split( /\-/, $fragB );
    $ionA =~ /(\++)/;
    my $chargeA = length($1) || 1;
    $ionA =~ s/\+//g;
    $ionB =~ /(\++)/;
    my $chargeB = length($1) || 1;
    $ionB =~ s/\+//g;
    my $comp = ( $ionA cmp $ionB ) || ( $chargeA <=> $chargeB );
    return $comp if ($comp);

    # Compares the number of losses
    $comp = @partA <=> @partB;
    return $comp if ($comp);

    # Prepares the losses for comparison
    my @lossA;
    foreach my $loss (@partA) {
        $loss =~ s/\*//g;
        if ( $loss =~ /(\d+)\((\w+)\)/ ) {

            # Multiple losses
            push( @lossA, [ $2, $1 ] );
        }
        else {

            # Single loss
            push( @lossA, [ $loss, 1 ] );
        }
    }
    @lossA = sort { $a->[0] cmp $b->[0] } @lossA;
    my @lossB;
    foreach my $loss (@partB) {
        $loss =~ s/\*//g;
        if ( $loss =~ /(\d+)\((\w+)\)/ ) {

            # Multiple losses
            push( @lossB, [ $2, $1 ] );
        }
        else {

            # Single loss
            push( @lossB, [ $loss, 1 ] );
        }
    }
    @lossB = sort { $a->[0] cmp $b->[0] } @lossB;

    # Compares the losses
    for ( my $i = 0 ; $i < @lossA ; $i++ ) {
        $comp = ( $lossA[$i][0] cmp $lossB[$i][0] )
          || ( $lossA[$i][1] <=> $lossB[$i][1] );
        return $comp if ($comp);
    }

    return 0;

}    # cmpFragTypes

sub getInternIons {
    my $this = shift;

    my $tol = $this->tol();
    my $minTol = $this->minTol() || 0.2;
    my @internIons;
    my %spectrum  = %{ $this->{theoSpectrum} };
    my $massIndex = $spectrum{intensityIndex};
    foreach my $frag ( keys( %{ $spectrum{mass}{intern} } ) ) {
        my $ionType = $spectrum{ionType}{$frag}[0];
        my ( @masses, @names, @matches );
        foreach my $aa ( sort keys( %{ $spectrum{mass}{intern}{$frag} } ) ) {
            push( @names,  $aa );
            push( @masses, $spectrum{mass}{intern}{$frag}{$aa} );
            if ( defined( $spectrum{match}{intern}{$frag}{$aa} ) ) {
                if ( defined($tol) ) {
                    my $theoMass = $spectrum{mass}{intern}{$frag}{$aa};
                    my $expMass =
                      $spectrum{match}{intern}{$frag}{$aa}[$massIndex];
                    my $error = $expMass - $theoMass;
                    if (
                        (
                            abs($error) / ( $theoMass + $expMass ) * 2e+6 <=
                            $tol
                        )
                        || ( abs($error) <= $minTol )
                      )
                    {
                        push( @matches, $spectrum{match}{intern}{$frag}{$aa} );
                    }
                    else {
                        push( @matches, undef );
                    }
                }
                else {
                    push( @matches, $spectrum{match}{intern}{$frag}{$aa} );
                }
            }
            else {
                push( @matches, undef );
            }
        }
        push(
            @internIons,
            {
                ionType        => $ionType,
                charge         => 1,
                names          => [@names],
                masses         => [@masses],
                massType       => $this->massType(),
                massIndex      => $massIndex,
                intensityIndex => $spectrum{intensityIndex},
                matches => ( defined( $spectrum{match} ) ? [@matches] : undef )
            }
        );
    }
#### changed by yiming;
    @internIons =
      sort { cmpFragTypes( $this->ionType($a), $this->ionType($b) ) }
      @internIons;
    return @internIons;

}    # getInternIons

sub getFragmentMasses {
    my ( %h, $hashref ) = @_;
    my %fragType = %{$hashref};
    my ( $pept, $modif, $frags, $spectrum ) =
      ( $h{pept}, $h{modif}, $h{fragTypes}, $h{spectrum} );

    # Cleans the spectrum hash just in case
    undef(%$spectrum);

    # Gets peptide sequence
    my $peptSeq;
    if ( ref($pept) ) {
        ## changed by xusheng
        #    if ($pept->isa('InSilicoSpectro::InSilico::Peptide') || $pept->isa('InSilicoSpectro::InSilico::AASequence')){
        $peptSeq = $pept->sequence();
        if ( !defined($modif) ) {
            $modif = $pept->modif();
        }

        #   }
        #  else{
        #    croak("Illegal peptide object [$pept]");
        # }
    }
    else {
        $peptSeq = $pept;
    }

    # Computes peptide mass
    my $len = length($peptSeq);
    my @modif;
    if ( defined($modif) ) {
        if ( ( my $ref = ref($modif) ) eq 'ARRAY' ) {

            # The modifs are given as a vector directly
            @modif = @$modif;
            croak( "Vector @$modif too long[" . join( ',', @modif ) . "]" )
              if ( @modif > $len + 2 );
        }
        elsif ( !$ref ) {

            # Scalar, assume string
            @modif = split( /:/, $modif );
        }
        else {
            croak("Unknown format for specifying modifications [$modif]");
        }
    }
    my $peptideMass = getPeptideMass( pept => $peptSeq, modif => \@modif );
    $spectrum->{peptideMass} = $peptideMass;
    $spectrum->{peptide}     = $pept;
    $spectrum->{modif}       = [@modif];

    # Computes the sums of amino acid masses
    my @pept = split( //, $peptSeq );
    my $mass = 0;
    $mass += getMass("mod_$modif[0]") if ( $modif[0] );    # N-term modif
    my @base;
    push( @base, 0 );    # for complete C-Term ions
    for ( my $i = 0 ; $i < @pept ; $i++ ) {
        $mass += getMass("aa_$pept[$i]");
        $mass += getMass("mod_$modif[$i+1]")
          if ( $modif[ $i + 1 ] );    # internal modif
        push( @base, $mass );
    }
    $base[-1] += getMass("mod_$modif[$len+1]")
      if ( $modif[ $len + 1 ] );      # C-term modif

    # Computes the fragments of each requested fragment type
    foreach my $frag (@$frags) {

        if ( $frag eq 'immo' ) {

            # Immonium ions

            my %already;
            for ( my $i = 1 ; $i < @pept - 1 ; $i++ ) {
                if ( defined( $immoAA{ $pept[$i] } ) ) {
                    my $actualAA = "$pept[$i]|$modif[$i+1]";
                    next if ( defined( $already{$actualAA} ) );

                    if (
                        !$modif[ $i + 1 ]
                        || (   ( $pept[$i] eq 'C' )
                            || ( $pept[$i] eq 'M' )
                            || ( $pept[$i] eq 'H' ) )
                      )
                    {
                        my $mass     = getMass("aa_$pept[$i]") + $immoDelta;
                        my $immoName = $pept[$i];
                        if ( $modif[ $i + 1 ] ) {
                            $immoName .= "+$modif[$i+1]";
                            $mass += getMass("mod_$modif[$i+1]");
                        }
                        $spectrum->{mass}{intern}{$frag}{$immoName} = $mass;
                        if ( $pept[$i] eq 'K' ) {

                            # Consider a possible extra mass with ammonia loss
                            $mass -= getMass('mol_NH3');
                            $spectrum->{mass}{intern}{$frag}{"$immoName-NH3"} =
                              $mass;
                        }
                        $already{$actualAA} = 1;
                        $spectrum->{ionType}{$frag}[0] = 'immo';
                    }
                }
            }
        }
        else {

            # Regular fragment types

            my $series    = $fragType{$frag}{series};
            my $charge    = $fragType{$frag}{charge};
            my $loss      = $fragType{$frag}{loss};
            my $firstFrag = $series{$series}{firstFrag};
            my $lastFrag  = $series{$series}{lastFrag};
            my $delta     = $series{$series}{delta}[$massType];

            if ( !defined($loss) ) {

                # no loss, straightforward computation
                $delta += ( $charge - 1 ) * getMass('el_H+');

                if ( $series{$series}{terminus} eq 'N' ) {

                    # N-term ions
                    $delta += $nTerm;
                    for (
                        my $i = $firstFrag ;
                        $i <= $len - $lastFrag + 1 ;
                        $i++
                      )
                    {
                        $spectrum->{mass}{term}{$frag}[ $i - 1 ] =
                          ( $base[$i] + $delta ) / $charge;
                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
                else {

                    # C-term ions, reverse and complement masses
                    for (
                        my $i = $firstFrag - 1 ;
                        $i < $len - $lastFrag + 1 ;
                        $i++
                      )
                    {
                        $spectrum->{mass}{term}{$frag}[$i] =
                          ( $peptideMass - $base[ $len - $i - 1 ] + $delta ) /
                          $charge;
                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
            }
            else {

                # Losses, possibly multiple and combined

                # Locates available positions for loss for each loss type (and checks all residues are different)
                my @loss = @$loss;
                my ( @avail, $nTotLoss, @distinctPos );
                my ( $startLoop, $endLoop ) =
                    ( $series{$series}{terminus} eq 'N' )
                  ? ( 0, $len - $lastFrag )
                  : ( $lastFrag - 1, $len - 1 );
                for ( my $i = $startLoop ; $i <= $endLoop ; $i++ ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
                        if ( $loss{ $loss[$j] }{residues}{ $pept[$i] } ) {
                            push( @{ $avail[$j] }, $i );
                            $nTotLoss++;
                            $distinctPos[$i]++;
                        }
                    }
                }
                next
                  if ( $nTotLoss == 0 )
                  ;    # this ion type is not possible for this peptide
                for ( my $i = 0 ; $i < @distinctPos ; $i++ ) {
                    if ( $distinctPos[$i] > 1 ) {
                        croak(
                            "Multiple loss at one single amino acid [$pept[$i]] in [$frag]"
                        );
                    }
                }

                # Computes maximum number of losses for each loss type
                my @maxLoss;
                for ( my $j = 0 ; $j < @loss ; $j++ ) {
                    my $repeat = $fragType{$frag}{repeat}[$j];
                    $repeat = $len if ( $repeat == -1 );
                    $maxLoss[$j] =
                      defined( $avail[$j] )
                      ? (
                        ( @{ $avail[$j] } < $repeat )
                        ? scalar( @{ $avail[$j] } )
                        : $repeat
                      )
                      : 0;
                }

                # Reverses the loss positions for C-term ions
                if ( $series{$series}{terminus} eq 'C' ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
                        @{ $avail[$j] } = reverse( @{ $avail[$j] } );
                    }
                }

                # Generates every combination of number of possible losses
                my $comb = 0;
                my @nLoss = split( //, '0' x @maxLoss );
                $nLoss[0] = 1;
                $delta += $nTerm if ( $series{$series}{terminus} eq 'N' );
                while (1) {

                    # Computes the fragments of the current combination

                    # Adapt delta to the number of losses
                    my $d = $delta;
                    for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                        $d += $nLoss[$i] * $loss{ $loss[$i] }{delta}[$massType];
                    }

                    if ( $series{$series}{terminus} eq 'N' ) {

                        # N-term ions

                        # First position that includes as many as required possible losses
                        my $rightMost = 0;
                        for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                            if (   ( $nLoss[$i] > 0 )
                                && ( $avail[$i][ $nLoss[$i] - 1 ] > $rightMost )
                              )
                            {
                                $rightMost = $avail[$i][ $nLoss[$i] - 1 ];
                            }
                        }

                        # Computes the fragments and check they are visible (firstFrag/lastFrag)
                        for (
                            my $i = $rightMost + 1 ;
                            $i <= $len - $lastFrag + 1 ;
                            $i++
                          )
                        {
                            if ( $i >= $firstFrag ) {
                                $spectrum->{mass}{term}{$frag}
                                  [ ( $comb * $len ) + $i - 1 ] =
                                  ( $base[$i] + $d ) / $charge;
                            }
                        }
                    }
                    else {

                        # C-term ions

                        # Last position that includes as many as required possible losses (from the right)
                        my $leftMost = $len - 1;
                        for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                            if (   ( $nLoss[$i] > 0 )
                                && ( $avail[$i][ $nLoss[$i] - 1 ] < $leftMost )
                              )
                            {
                                $leftMost = $avail[$i][ $nLoss[$i] - 1 ];
                            }
                        }

                        # Computes the fragments and check they are visible (firstFrag/lastFrag)
                        for (
                            my $i = $len - $leftMost - 1 ;
                            $i < $len - $lastFrag + 1 ;
                            $i++
                          )
                        {
                            if ( $i >= $firstFrag - 1 ) {
                                $spectrum->{mass}{term}{$frag}
                                  [ ( $comb * $len ) + $i ] =
                                  ( $peptideMass - $base[ $len - $i - 1 ] + $d )
                                  / $charge;
                            }
                        }
                    }

                    # Computes the exact ion type and saves its name
                    my $ionType = $series;
                    for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                        if ( $nLoss[$i] > 1 ) {
                            $ionType .= "-$nLoss[$i]($loss[$i])";
                        }
                        elsif ( $nLoss[$i] == 1 ) {
                            $ionType .= "-$loss[$i]";
                        }
                    }
                    $spectrum->{ionType}{$frag}[$comb] = $ionType;
                    $comb++;

                    # Gets the next combination
                    my $i;
                    for (
                        $i = 0 ;
                        ( $i < @maxLoss ) && ( $nLoss[$i] == $maxLoss[$i] ) ;
                        $i++
                      )
                    {
                    }
                    last if ( $i == @maxLoss );
                    $nLoss[$i]++;
                    for ( my $j = 0 ; $j < $i ; $j++ ) {
                        $nLoss[$j] = 0;
                    }
                }
            }
        }
    }

}    # getFragmentMasses

=pod
sub getPeptideMass {
    my (%h) = @_;
    my ( $pept, $modif, $termGain ) = ( $h{pept}, $h{modif}, $h{termGain} );
    croak("No peptide given in getPeptideMass") unless ( defined($pept) );

#  if ((index($pept, 'B') == -1) && (index($pept, 'Z') == -1) && (index($pept, 'X') == -1)){
    my %hEls;
    ## changed by xusheng;
    # unless ($pept=~/InSilicoSpectro::InSilico::AASequence::qrValidAASeq/){
    if ( defined($modif) ) {
        foreach (@$modif) {
            $hEls{"mod_$_"}++ if ($_);
        }
    }
    foreach ( split( //, $pept ) ) {
        $hEls{"aa_$_"}++;
    }
    my $m = $termGain || $nTerm + $cTerm;
    foreach ( keys %hEls ) {
        $m += $hEls{$_} * getMass($_);
    }
    return $m;
}

#else{
# No defined mass
# warn "no mass for peptide [$pept]";
# return -1.0;
#}
=cut

sub getMass {

    # Returns the mass of the given element/molecule, according to $massType or to the second argument if defined
    my ( $el, $mt ) = @_;

    $mt = $massType if ( !defined($mt) );
    defined( $elMass{$el} )
      || ( $invalidElementCall
        && $invalidElementCall->("Unknown element/molecule in getMass: [$el]")
      );
    $elMass{$el}[$mt];

}    # getMass

1;

__END__
