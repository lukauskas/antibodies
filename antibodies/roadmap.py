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
            if vendor in {'Cell Signaling', 'Cell Signalling', 'Cell Signaling Technology'}:
                vendor = 'Cell Signaling Technology'
            elif vendor in {'diagenode', 'Diagenode', 'Diagenonde'}:
                vendor = 'Diagenode'
            elif vendor in {'none', '0'}:
                vendor = None
            elif vendor in {'Millipore', 'Upstate', 'Upstate or Millipore'}:
                vendor = 'Millipore/Upstate'

            lot_number = row['chip_antibody_lot']

            if lot_number in {'none', 'Unknown', '0'}:
                lot_number = None
            elif lot_number.lower().startswith('lot'):
                lot_number = lot_number[3:]
                lot_number = lot_number.lstrip('#')

            catalog_id = row['chip_antibody_catalog']
            if catalog_id in {'0', 'none'}:
                catalog_id = None

            if catalog_id:
                if vendor == 'Abcam':
                    if not catalog_id.startswith('EP'):
                        catalog_id = catalog_id.lower()
                        if not catalog_id.startswith('ab'):
                            catalog_id = 'ab' + str(int(catalog_id))

                elif vendor == 'Cell Signaling Technology':
                    catalog_id = catalog_id.upper()
                    if not catalog_id.endswith('S'):
                        if catalog_id.endswith('B'):
                            catalog_id = catalog_id.rstrip('B')  # I think this is a typo

                        catalog_id = str(int(catalog_id)) + 'S'
                elif vendor == 'Active Motif':
                    catalog_id = catalog_id.upper()
                    if not catalog_id.startswith('AM'):
                        catalog_id = 'AM' + str(int(catalog_id))
                elif vendor == 'Millipore/Upstate':
                    if catalog_id.startswith('Upstate') or catalog_id.startswith('upstate'):
                        catalog_id = catalog_id[len('upstate'):].strip()
                    if catalog_id.startswith('Millipore'):
                        catalog_id = catalog_id[len('Millipore'):].strip()

            data.append({'Antibody Target': target,
                         'Antibody Vendor': vendor,
                         'Lot Number': lot_number,
                         'Catalog ID': catalog_id})

        return pd.DataFrame(data, index=index)


    def run(self):
        COLUMNS_TO_KEEP = ['Sample Name',
                           'Experiment', 'NA Accession', 'Center',
                           'chip_protocol',
                           'chip_protocol_antibody_amount',
                           'chip_protocol_bead_amount',
                           'chip_protocol_bead_type',
                           'chip_protocol_chromatin_amount',
                           'collection_method',
                           'differentiation_stage',
                           'disease',
                           'donor_age',
                           'donor_ethnicity',
                           'donor_id',
                           'extraction_protocol_sonication_cycles',
                           'extraction_protocol_type_of_sonicator',
                           'size_fraction',
                           ]

        COLUMNS_TO_RENAME = [('biomaterial_provider', 'Biomaterial Provider')]

        with self.input().open('r') as in_:
            original_df = pd.read_csv(in_).set_index('GEO Accession')

            df = original_df[COLUMNS_TO_KEEP].copy()

            df['Sex'] = self._fix_sex(original_df)

            for old_column, new_column in COLUMNS_TO_RENAME:
                df[new_column] = original_df[old_column]

            df = df.join(self._fixed_antibody_df(original_df))

            with self.output().open('w') as out_:
                df.to_csv(out_)



if __name__ == '__main__':
    LOGGER.setLevel(logging.DEBUG)
    logging.basicConfig()

    luigi.run(main_task_cls=CleanRoadmapData)


