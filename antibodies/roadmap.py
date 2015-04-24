from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from io import StringIO
import luigi

import requests
import pandas as pd
import logging
import numpy as np
import time
from antibodies.cleanup import cleanup_vendor_name, cleanup_catalog_id, cleanup_lot_number
from antibodies.validity.antibody_validation_database import AntibodyValidationDatabase
from antibodies.validity.citeab import parse_number_of_citations_from_citeab, CitationsTaskBase

LOGGER = logging.getLogger('roadmap_antibodies')

NCBI_DATA_URI = 'http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples&sort=acc&mode=csv'
GEO_QUERY_URL = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=brief'


def _fetch_data_from_ncbi():
    LOGGER.debug('Fetching {}'.format(NCBI_DATA_URI))
    response = requests.get(NCBI_DATA_URI)
    response.raise_for_status()

    text = StringIO(response.text)

    data = pd.read_csv(text)

    # Strip the hashes form headers
    data.columns = [x.strip(' #') for x in data.columns]

    return data


def _query_geo(geo_id):

    url = GEO_QUERY_URL.format(geo_id)

    LOGGER.debug('Fetching {}'.format(url))

    response = requests.get(url)
    response.raise_for_status()

    text = response.text

    row = {}
    row['GEO Accession'] = geo_id

    for line in text.split('\n'):
        line = line.strip()

        key, sep, description = line.partition(' = ')
        if key == '!Sample_characteristics_ch1':
            sub_key, __, sub_value = description.partition(':')
            sub_value = sub_value.strip()
            if sub_key in row:
                sub_value = row[sub_key] + '\n' + sub_value
            row[sub_key] = sub_value
    return row

def _process_roadmap_data():

    EXPERIMENTS_TO_SKIP = {'Bisulfite-Seq', 'Digital genomic footprinting',
                            'DNase hypersensitivity', 'Exon array', 'Expression array', 'Genotyping array',
                            'MeDIP-Seq', 'MRE-Seq', 'mRNA-Seq', 'RRBS', 'smRNA-Seq', 'TAB-Seq', 'Whole genome sequencing'}

    data_list = _fetch_data_from_ncbi()
    data_list = data_list[~data_list['Experiment'].isin(EXPERIMENTS_TO_SKIP)]

    data_list = data_list[['GEO Accession', 'Sample Name', 'Experiment', 'NA Accession', 'Center']]
    data_list = data_list.set_index('GEO Accession')

    antibody_data = []
    for i, geo_id in enumerate(data_list.index):
        LOGGER.info('GEO ID {}/{} ({:.2f}%)'.format(i+1, len(data_list), 100 * (i+1) / len(data_list)))
        antibody_data.append(_query_geo(geo_id))

    antibody_data = pd.DataFrame(antibody_data).set_index('GEO Accession')

    return data_list.join(antibody_data)

class RawRoadmapData(luigi.Task):


    def output(self):
        return luigi.File('roadmap_antibodies.raw.csv')

    def run(self):
        with self.output().open('w') as f:
            df = _process_roadmap_data()
            df.to_csv(f, encoding='utf-8')

class CleanRoadmapData(luigi.Task):

    def requires(self):
        return RawRoadmapData()

    def output(self):
        return luigi.File('roadmap_antibodies.csv')

    def _fix_sex(self, original_df):

        sex_df = original_df[['Sex', 'donor_sex', 'sex']]

        new_values = []
        index = []
        for ix, row in sex_df.iterrows():

            values = filter(lambda x: isinstance(x, basestring), row)
            values = map(lambda x: x.lower(), values)
            if 'male' in values:
                new_values.append('Male')
                index.append(ix)
            elif 'female' in values:
                new_values.append('Female')
                index.append(ix)

        return pd.Series(new_values, index=index)


    def _fixed_chip_protocol_df(self, original_df):

        protocol_df = original_df[['chip_protocol', 'chip_protocol_antibody_amount',
                                   'chip_protocol_bead_amount', 'chip_protocol_bead_type',
                                   'chip_protocol_chromatin_amount']].fillna('')

        index = []
        data = []

        protocol_substitutions = {'See http://bioinformatics-renlab.ucsd.edu/RenLabChipProtocolV1.pdf': 'Ren',
                                  'Farnham Lab Protocol': 'Farnham',
                                  'Farnham lab protocol': 'Farnham',
                                  'Input': 'Input',
                                  'Bernstein_BROAD_ENCODE_protocol': 'Bernstein',
                                  'BCCAGSC ChIP Standard Operating Procedure': 'BCCAGSC',
                                  'http://www.roadmapepigenomics.org/protocolstype/experimental': 'Roadmap Epigenomics'}

        for ix, row in protocol_df.iterrows():
            index.append(ix)
            protocol = row['chip_protocol']
            if protocol in {'standard'}:
                protocol = None
            if protocol:
                protocol = protocol_substitutions[protocol]


            antibody_amount_str = row['chip_protocol_antibody_amount']
            antibody_amount, antibody_units = None, None

            if antibody_amount_str in {'n/a', 'standard'}:
                antibody_amount_str = None

            if antibody_amount_str:
                antibody_amount_str = antibody_amount_str.strip('~').replace(' ', '')

                if 'micrograms' in antibody_amount_str:
                    antibody_amount = float(antibody_amount_str[:-len('micrograms')].strip())
                    antibody_units = 'ug'
                elif 'ug' in antibody_amount_str:
                    antibody_units = 'ug'
                    antibody_amount = antibody_amount_str[:-len('ug')]
                elif 'ul' in antibody_amount_str:
                    antibody_units = 'ul'
                    antibody_amount = antibody_amount_str[:-len('ul')]
                else:
                    raise Exception('Cannot parse antibody amount from {!r}'.format(antibody_amount_str))

            beads_amount_str = row['chip_protocol_bead_amount']
            if beads_amount_str in {'standard'}:
                beads_amount_str = None

            if beads_amount_str:
                beads_amount_str = beads_amount_str.lower().replace('(cst)', '')\
                                                           .replace('bed volume', '').replace('a/g bead slurry', '')
                beads_amount_str = beads_amount_str.replace(' ', '')

                if ',' in beads_amount_str:
                    beads_amount_units = 'count'
                    beads_amount = int(beads_amount_str.replace(',', ''))
                elif 'ul' in beads_amount_str:
                    beads_amount = float(beads_amount_str[:-len('ul')])
                    beads_amount_units = 'ul'
                else:
                    raise Exception('Cannot parse beads amount from {!r}'.format(beads_amount_str))

            chromatin_amount_str = row['chip_protocol_chromatin_amount']

            if chromatin_amount_str in {'standard'}:
                chromatin_amount_str = None

            if chromatin_amount_str:
                chromatin_amount_str = chromatin_amount_str.lower()
                chromatin_amount_str = chromatin_amount_str.replace(' ', '')
                chromatin_amount_str = chromatin_amount_str.replace('to', '-')

                if 'micrograms' in chromatin_amount_str:
                    chromatin_amount_units = 'ug'
                    chromatin_amount = float(chromatin_amount_str[:-len('micrograms')].strip())
                elif 'millioncells' in chromatin_amount_str:
                    chromatin_amount_units = 'million cells'
                    chromatin_amount = chromatin_amount_str[:-len('millioncells')].strip()
                elif 'ug' in chromatin_amount_str:
                    chromatin_amount_units = 'ug'
                    chromatin_amount = chromatin_amount_str[:-len('ug')]
                elif 'ng' in chromatin_amount_str:
                    chromatin_amount_units = 'ng'
                    chromatin_amount = chromatin_amount_str[:-len('ng')]
                else:
                    raise Exception('Cannot parse chromatin amount: {!r}'.format(chromatin_amount_str))

            data.append({'Chromatin Amount': chromatin_amount,
                         'Chromatin Amount Units': chromatin_amount_units,
                         'Beads Amount': beads_amount,
                         'Beads Amount Units': beads_amount_units,
                         'Antibody Amount': antibody_amount,
                         'Antibody Amount Units': antibody_units,
                         'Protocol': protocol})

        return pd.DataFrame(data, index=index)


    def _fixed_antibody_df(self, original_df):

        antibody_data = original_df[['chip_antibody', 'chip_antibody_catalog', 'chip_antibody_lot',
                                     'chip_antibody_provider']].fillna('')

        index = []
        data = []

        for ix, row in antibody_data.iterrows():
            index.append(ix)

            target = row['chip_antibody']
            if isinstance(target, basestring):
                target = target.replace('Me', 'me').replace('Ac', 'ac')
                if target.startswith('Histone'):
                    target = target[len('Histone'):].strip()

            vendor = row['chip_antibody_provider']
            vendor = cleanup_vendor_name(vendor)

            catalog_id = row['chip_antibody_catalog']
            catalog_id = cleanup_catalog_id(catalog_id, vendor)

            lot_number = row['chip_antibody_lot']
            lot_number = cleanup_lot_number(vendor, catalog_id, lot_number)

            data.append({'Antibody Target': target,
                         'Antibody Vendor': vendor,
                         'Antibody Lot Number': lot_number,
                         'Antibody Catalogue ID': catalog_id})

        return pd.DataFrame(data, index=index)


    def run(self):
        COLUMNS_TO_KEEP = ['Sample Name',
                           'Experiment', 'Center']

        COLUMNS_TO_RENAME = [('biomaterial_provider', 'Biomaterial Provider')]

        with self.input().open('r') as in_:
            original_df = pd.read_csv(in_).set_index('GEO Accession')

            df = original_df[COLUMNS_TO_KEEP].copy()

            df['Sex'] = self._fix_sex(original_df)

            for old_column, new_column in COLUMNS_TO_RENAME:
                df[new_column] = original_df[old_column]

            df = df.join(self._fixed_antibody_df(original_df))
            df = df.join(self._fixed_chip_protocol_df(original_df))

            column_order = ['Sample Name',
                            'Sex',
                            'Center',
                            'Experiment',
                            'Protocol',
                            'Antibody Target',
                            'Antibody Vendor',
                            'Antibody Catalogue ID',
                            'Antibody Lot Number',
                            'Antibody Amount',
                            'Antibody Amount Units',
                            'Beads Amount',
                            'Beads Amount Units',
                            'Chromatin Amount',
                            'Chromatin Amount Units',
                            ]
            df = df[column_order]

            with self.output().open('w') as out_:
                df.to_csv(out_)

class RoadmapCitations(CitationsTaskBase):

    def requires(self):
        return CleanRoadmapData()

    def output(self):
        return luigi.File('roadmap_citations.csv')

class RoadmapSummary(luigi.Task):

    def requires(self):
        return [CleanRoadmapData(), RoadmapCitations(), AntibodyValidationDatabase()]

    def output(self):
        return luigi.File('roadmap_summary.csv')

    def run(self):
        roadmap_data_output, roadmap_citations_output, antibody_validation_database_output = self.input()

        with roadmap_data_output.open('r') as f:
            roadmap_antibodies = pd.read_csv(f)

        with roadmap_citations_output.open('r') as f:
            roadmap_citations = pd.read_csv(f, index_col=['Antibody Vendor', 'Antibody Catalogue ID'])

        with antibody_validation_database_output.open('r') as f:
            antibody_validation_database = pd.read_csv(f)


        number_of_experiments_used = roadmap_antibodies.groupby(['Experiment', 'Antibody Target', 'Antibody Vendor',
                                                                 'Antibody Catalogue ID', 'Antibody Lot Number']).aggregate(len)
        number_of_experiments_used = number_of_experiments_used[number_of_experiments_used.columns[0]]
        number_of_experiments_used.name = 'Number of Roadmap Experiments'
        number_of_experiments_used = number_of_experiments_used

        h9_abs = roadmap_antibodies[roadmap_antibodies['Sample Name'] == 'H9 cell line'][['Antibody Vendor', 'Antibody Catalogue ID', 'Antibody Lot Number']].drop_duplicates()
        h9_abs = {tuple(x) for __, x in h9_abs.iterrows()}
        used_for_h9 = pd.Series([tuple(x) in h9_abs for __, x in number_of_experiments_used.reset_index()[['Antibody Vendor', 'Antibody Catalogue ID', 'Antibody Lot Number']].iterrows()], index=number_of_experiments_used.index, name='Used for H9')

        metadata = pd.concat((used_for_h9, number_of_experiments_used), axis=1)

        metadata = metadata.reset_index().join(roadmap_citations,
                                               on=['Antibody Vendor', 'Antibody Catalogue ID'],
                                               how='left')

        antibody_validation_database['Validation Result'] = antibody_validation_database['Validation Result'].fillna('Unknown')

        grouped_validations = antibody_validation_database.groupby(['Vendor', 'Catalogue ID', 'Lot Number',
                                           'Validation', 'Validation Result']).count().max(axis=1)

        grouped_validations.name = 'Count'
        validations_pivot_lot_number = pd.pivot_table(grouped_validations.reset_index(), values='Count',
                                                      index=['Vendor', 'Catalogue ID',
                                                             'Lot Number'
                                                            ],
                                                      columns=['Validation', 'Validation Result'],
                                                      aggfunc=np.sum, fill_value=0)
        validations_pivot_no_lot = pd.pivot_table(grouped_validations.reset_index(), values='Count',
                                                  index=['Vendor', 'Catalogue ID',
                                                        ],
                                                  columns=['Validation', 'Validation Result'],
                                                  aggfunc=np.sum, fill_value=0)

        full_data = metadata.join(validations_pivot_lot_number,
                                                          on=['Antibody Vendor',
                                                              'Antibody Catalogue ID', 'Antibody Lot Number'],
                                                          how='left')
        full_data.columns = [x if not isinstance(x, tuple) else ('Lot', x[0], x[1]) for x in full_data.columns]
        full_data = full_data.join(validations_pivot_no_lot, on=['Antibody Vendor', 'Antibody Catalogue ID'],
              how='left').set_index(['Experiment', 'Antibody Target', 'Antibody Vendor', 'Antibody Catalogue ID', 'Antibody Lot Number'])
        full_data.columns = [x if not isinstance(x, tuple) or x[0] == 'Lot' else ('Antibody', x[0], x[1])
                             for x in full_data.columns]


        full_data.columns = pd.MultiIndex.from_tuples([x if isinstance(x, tuple) else (None, None,x)
                                                       for x in full_data.columns])

        with self.output().open('w') as f:
            full_data.to_csv(f)

if __name__ == '__main__':
    LOGGER.setLevel(logging.DEBUG)
    logging.basicConfig()

    luigi.run()


