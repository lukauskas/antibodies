from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from io import StringIO
import luigi

import requests
import pandas as pd
import logging

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

    data_list = _fetch_data_from_ncbi()
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
        return luigi.File('roadmap_data.raw.csv')

    def run(self):
        with self.output().open('w') as f:
            df = _process_roadmap_data()
            df.to_csv(f, encoding='utf-8')

if __name__ == '__main__':
    LOGGER.setLevel(logging.DEBUG)
    logging.basicConfig()

    luigi.run(main_task_cls=RawRoadmapData)


