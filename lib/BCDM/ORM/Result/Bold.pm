use utf8;
package BCDM::ORM::Result::Bold;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

BCDM::ORM::Result::Bold

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<bold>

=cut

__PACKAGE__->table("bold");

=head1 ACCESSORS

=head2 recordid

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 taxonid

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 processid

  data_type: 'text'
  is_nullable: 1

=head2 sampleid

  data_type: 'text'
  is_nullable: 1

=head2 fieldid

  data_type: 'text'
  is_nullable: 1

=head2 museumid

  data_type: 'text'
  is_nullable: 1

=head2 record_id

  data_type: 'text'
  is_nullable: 1

=head2 specimenid

  data_type: 'text'
  is_nullable: 1

=head2 processid_minted_date

  data_type: 'text'
  is_nullable: 1

=head2 bin_uri

  data_type: 'text'
  is_nullable: 1

=head2 bin_created_date

  data_type: 'text'
  is_nullable: 1

=head2 collection_code

  data_type: 'text'
  is_nullable: 1

=head2 inst

  data_type: 'text'
  is_nullable: 1

=head2 taxid

  data_type: 'text'
  is_nullable: 1

=head2 taxon_name

  data_type: 'text'
  is_nullable: 1

=head2 taxon_rank

  data_type: 'text'
  is_nullable: 1

=head2 kingdom

  data_type: 'text'
  is_nullable: 1

=head2 phylum

  data_type: 'text'
  is_nullable: 1

=head2 class

  data_type: 'text'
  is_nullable: 1

=head2 order

  data_type: 'text'
  is_nullable: 1

=head2 family

  data_type: 'text'
  is_nullable: 1

=head2 subfamily

  data_type: 'text'
  is_nullable: 1

=head2 tribe

  data_type: 'text'
  is_nullable: 1

=head2 genus

  data_type: 'text'
  is_nullable: 1

=head2 species

  data_type: 'text'
  is_nullable: 1

=head2 subspecies

  data_type: 'text'
  is_nullable: 1

=head2 species_reference

  data_type: 'text'
  is_nullable: 1

=head2 identification

  data_type: 'text'
  is_nullable: 1

=head2 identification_method

  data_type: 'text'
  is_nullable: 1

=head2 identification_rank

  data_type: 'text'
  is_nullable: 1

=head2 identified_by

  data_type: 'text'
  is_nullable: 1

=head2 identifier_email

  data_type: 'text'
  is_nullable: 1

=head2 taxonomy_notes

  data_type: 'text'
  is_nullable: 1

=head2 sex

  data_type: 'text'
  is_nullable: 1

=head2 reproduction

  data_type: 'text'
  is_nullable: 1

=head2 life_stage

  data_type: 'text'
  is_nullable: 1

=head2 short_note

  data_type: 'text'
  is_nullable: 1

=head2 notes

  data_type: 'text'
  is_nullable: 1

=head2 voucher_type

  data_type: 'text'
  is_nullable: 1

=head2 tissue_type

  data_type: 'text'
  is_nullable: 1

=head2 specimen_linkout

  data_type: 'text'
  is_nullable: 1

=head2 associated_specimens

  data_type: 'text'
  is_nullable: 1

=head2 associated_taxa

  data_type: 'text'
  is_nullable: 1

=head2 collectors

  data_type: 'text'
  is_nullable: 1

=head2 collection_date_start

  data_type: 'text'
  is_nullable: 1

=head2 collection_date_end

  data_type: 'text'
  is_nullable: 1

=head2 collection_event_id

  data_type: 'text'
  is_nullable: 1

=head2 collection_time

  data_type: 'text'
  is_nullable: 1

=head2 collection_notes

  data_type: 'text'
  is_nullable: 1

=head2 geoid

  data_type: 'text'
  is_nullable: 1

=head2 country/ocean

  accessor: 'country_ocean'
  data_type: 'text'
  is_nullable: 1

=head2 country_iso

  data_type: 'text'
  is_nullable: 1

=head2 province/state

  accessor: 'province_state'
  data_type: 'text'
  is_nullable: 1

=head2 region

  data_type: 'text'
  is_nullable: 1

=head2 sector

  data_type: 'text'
  is_nullable: 1

=head2 site

  data_type: 'text'
  is_nullable: 1

=head2 site_code

  data_type: 'text'
  is_nullable: 1

=head2 coord

  data_type: 'text'
  is_nullable: 1

=head2 coord_accuracy

  data_type: 'text'
  is_nullable: 1

=head2 coord_source

  data_type: 'text'
  is_nullable: 1

=head2 elev

  data_type: 'text'
  is_nullable: 1

=head2 elev_accuracy

  data_type: 'text'
  is_nullable: 1

=head2 depth

  data_type: 'text'
  is_nullable: 1

=head2 depth_accuracy

  data_type: 'text'
  is_nullable: 1

=head2 habitat

  data_type: 'text'
  is_nullable: 1

=head2 sampling_protocol

  data_type: 'text'
  is_nullable: 1

=head2 nuc

  data_type: 'text'
  is_nullable: 1

=head2 nuc_basecount

  data_type: 'text'
  is_nullable: 1

=head2 insdc_acs

  data_type: 'text'
  is_nullable: 1

=head2 funding_src

  data_type: 'text'
  is_nullable: 1

=head2 marker_code

  data_type: 'text'
  is_nullable: 1

=head2 primers_forward

  data_type: 'text'
  is_nullable: 1

=head2 primers_reverse

  data_type: 'text'
  is_nullable: 1

=head2 sequence_run_site

  data_type: 'text'
  is_nullable: 1

=head2 sequence_upload_date

  data_type: 'text'
  is_nullable: 1

=head2 bold_recordset_code_arr

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "recordid",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "taxonid",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "processid",
  { data_type => "text", is_nullable => 1 },
  "sampleid",
  { data_type => "text", is_nullable => 1 },
  "fieldid",
  { data_type => "text", is_nullable => 1 },
  "museumid",
  { data_type => "text", is_nullable => 1 },
  "record_id",
  { data_type => "text", is_nullable => 1 },
  "specimenid",
  { data_type => "text", is_nullable => 1 },
  "processid_minted_date",
  { data_type => "text", is_nullable => 1 },
  "bin_uri",
  { data_type => "text", is_nullable => 1 },
  "bin_created_date",
  { data_type => "text", is_nullable => 1 },
  "collection_code",
  { data_type => "text", is_nullable => 1 },
  "inst",
  { data_type => "text", is_nullable => 1 },
  "taxid",
  { data_type => "text", is_nullable => 1 },
  "taxon_name",
  { data_type => "text", is_nullable => 1 },
  "taxon_rank",
  { data_type => "text", is_nullable => 1 },
  "kingdom",
  { data_type => "text", is_nullable => 1 },
  "phylum",
  { data_type => "text", is_nullable => 1 },
  "class",
  { data_type => "text", is_nullable => 1 },
  "order",
  { data_type => "text", is_nullable => 1 },
  "family",
  { data_type => "text", is_nullable => 1 },
  "subfamily",
  { data_type => "text", is_nullable => 1 },
  "tribe",
  { data_type => "text", is_nullable => 1 },
  "genus",
  { data_type => "text", is_nullable => 1 },
  "species",
  { data_type => "text", is_nullable => 1 },
  "subspecies",
  { data_type => "text", is_nullable => 1 },
  "species_reference",
  { data_type => "text", is_nullable => 1 },
  "identification",
  { data_type => "text", is_nullable => 1 },
  "identification_method",
  { data_type => "text", is_nullable => 1 },
  "identification_rank",
  { data_type => "text", is_nullable => 1 },
  "identified_by",
  { data_type => "text", is_nullable => 1 },
  "identifier_email",
  { data_type => "text", is_nullable => 1 },
  "taxonomy_notes",
  { data_type => "text", is_nullable => 1 },
  "sex",
  { data_type => "text", is_nullable => 1 },
  "reproduction",
  { data_type => "text", is_nullable => 1 },
  "life_stage",
  { data_type => "text", is_nullable => 1 },
  "short_note",
  { data_type => "text", is_nullable => 1 },
  "notes",
  { data_type => "text", is_nullable => 1 },
  "voucher_type",
  { data_type => "text", is_nullable => 1 },
  "tissue_type",
  { data_type => "text", is_nullable => 1 },
  "specimen_linkout",
  { data_type => "text", is_nullable => 1 },
  "associated_specimens",
  { data_type => "text", is_nullable => 1 },
  "associated_taxa",
  { data_type => "text", is_nullable => 1 },
  "collectors",
  { data_type => "text", is_nullable => 1 },
  "collection_date_start",
  { data_type => "text", is_nullable => 1 },
  "collection_date_end",
  { data_type => "text", is_nullable => 1 },
  "collection_event_id",
  { data_type => "text", is_nullable => 1 },
  "collection_time",
  { data_type => "text", is_nullable => 1 },
  "collection_notes",
  { data_type => "text", is_nullable => 1 },
  "geoid",
  { data_type => "text", is_nullable => 1 },
  "country/ocean",
  { accessor => "country_ocean", data_type => "text", is_nullable => 1 },
  "country_iso",
  { data_type => "text", is_nullable => 1 },
  "province/state",
  { accessor => "province_state", data_type => "text", is_nullable => 1 },
  "region",
  { data_type => "text", is_nullable => 1 },
  "sector",
  { data_type => "text", is_nullable => 1 },
  "site",
  { data_type => "text", is_nullable => 1 },
  "site_code",
  { data_type => "text", is_nullable => 1 },
  "coord",
  { data_type => "text", is_nullable => 1 },
  "coord_accuracy",
  { data_type => "text", is_nullable => 1 },
  "coord_source",
  { data_type => "text", is_nullable => 1 },
  "elev",
  { data_type => "text", is_nullable => 1 },
  "elev_accuracy",
  { data_type => "text", is_nullable => 1 },
  "depth",
  { data_type => "text", is_nullable => 1 },
  "depth_accuracy",
  { data_type => "text", is_nullable => 1 },
  "habitat",
  { data_type => "text", is_nullable => 1 },
  "sampling_protocol",
  { data_type => "text", is_nullable => 1 },
  "nuc",
  { data_type => "text", is_nullable => 1 },
  "nuc_basecount",
  { data_type => "text", is_nullable => 1 },
  "insdc_acs",
  { data_type => "text", is_nullable => 1 },
  "funding_src",
  { data_type => "text", is_nullable => 1 },
  "marker_code",
  { data_type => "text", is_nullable => 1 },
  "primers_forward",
  { data_type => "text", is_nullable => 1 },
  "primers_reverse",
  { data_type => "text", is_nullable => 1 },
  "sequence_run_site",
  { data_type => "text", is_nullable => 1 },
  "sequence_upload_date",
  { data_type => "text", is_nullable => 1 },
  "bold_recordset_code_arr",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</recordid>

=back

=cut

__PACKAGE__->set_primary_key("recordid");

=head1 RELATIONS

=head2 bold_criterias

Type: has_many

Related object: L<BCDM::ORM::Result::BoldCriteria>

=cut

__PACKAGE__->has_many(
  "bold_criterias",
  "BCDM::ORM::Result::BoldCriteria",
  { "foreign.recordid" => "self.recordid" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 taxonid

Type: belongs_to

Related object: L<BCDM::ORM::Result::Taxa>

=cut

__PACKAGE__->belongs_to(
  "taxonid",
  "BCDM::ORM::Result::Taxa",
  { taxonid => "taxonid" },
  {
    is_deferrable => 0,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07049 @ 2024-07-12 17:27:42
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:EIuy1NZgdbvVBTQvqISLdA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
