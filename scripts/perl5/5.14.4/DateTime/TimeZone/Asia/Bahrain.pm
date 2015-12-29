# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/_mzyzyR3wa/asia.  Olson data version 2014g
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Asia::Bahrain;
$DateTime::TimeZone::Asia::Bahrain::VERSION = '1.74';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Asia::Bahrain::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
60557747860, #      utc_end 1919-12-31 20:37:40 (Wed)
DateTime::TimeZone::NEG_INFINITY, #  local_start
60557760000, #    local_end 1920-01-01 00:00:00 (Thu)
12140,
0,
'LMT',
    ],
    [
60557747860, #    utc_start 1919-12-31 20:37:40 (Wed)
62211873600, #      utc_end 1972-05-31 20:00:00 (Wed)
60557762260, #  local_start 1920-01-01 00:37:40 (Thu)
62211888000, #    local_end 1972-06-01 00:00:00 (Thu)
14400,
0,
'GST',
    ],
    [
62211873600, #    utc_start 1972-05-31 20:00:00 (Wed)
DateTime::TimeZone::INFINITY, #      utc_end
62211884400, #  local_start 1972-05-31 23:00:00 (Wed)
DateTime::TimeZone::INFINITY, #    local_end
10800,
0,
'AST',
    ],
];

sub olson_version { '2014g' }

sub has_dst_changes { 0 }

sub _max_year { 2024 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}



1;

