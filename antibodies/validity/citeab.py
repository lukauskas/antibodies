from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import itertools
import luigi
import lxml.html
import requests
import time
import pandas as pd

SEARCH_URI = 'http://www.citeab.com/search?q={}'

LOGGER = logging.getLogger('antibodies.validity.citeab')

# FROM: http://whatsmyuseragent.com/CommonUserAgents
USER_AGENTS = itertools.cycle([
    'Mozilla/5.0 ;Windows NT 6.1; WOW64; Trident/7.0; rv:11.0; like Gecko',
    'Mozilla/5.0 ;Windows NT 6.2; WOW64; rv:27.0; Gecko/20100101 Firefox/27.0',
    'Mozilla/5.0 ;Windows NT 6.1; WOW64; AppleWebKit/537.36 ;KHTML, like Gecko; Chrome/36.0.1985.143 Safari/537.36',
    'Mozilla/5.0 ;iPhone; CPU iPhone OS 8_1_2 like Mac OS X; AppleWebKit/600.1.4 ;KHTML, like Gecko; Version/8.0 Mobile/12B440 Safari/600.1.4',
    'Mozilla/5.0 ;Macintosh; Intel Mac OS X 10_10_2; AppleWebKit/537.36 ;KHTML, like Gecko; Chrome/40.0.2214.111 Safari/537.36',
    'Mozilla/5.0 ;Windows NT 6.3; WOW64; Trident/7.0; rv:11.0; like Gecko',
    'Mozilla/5.0 ;iPhone; CPU iPhone OS 8_1_3 like Mac OS X; AppleWebKit/600.1.4 ;KHTML, like Gecko; Version/8.0 Mobile/12B466 Safari/600.1.4'
])

def parse_number_of_citations_from_citeab(antibody_vendor, antibody_code):
    logger = LOGGER
    query_vendor = antibody_vendor
    query_code = antibody_code

    vendor_is_correct = lambda x: antibody_vendor == x
    code_is_correct = lambda x: antibody_code == x

    if antibody_vendor == 'Active Motif' and antibody_code.startswith('AM'):
        query_code = antibody_code[2:]
        code_is_correct = lambda x: antibody_code[2:] == x

    elif antibody_vendor == 'Millipore/Upstate':
        query_vendor = 'Millipore'
        vendor_is_correct = lambda x: 'Millipore' in x

    elif antibody_vendor == 'Cell Signaling Technology':
        query_code = antibody_code.rstrip('S')
        code_is_correct = lambda x: x == antibody_code.rstrip('S')

    user_agent = USER_AGENTS.next()

    response = requests.get(SEARCH_URI.format('{} {}'.format(query_code, query_vendor)),
                            headers={'User-Agent': user_agent})

    response.raise_for_status()
    if 'Your search returned no matches.' in response.text:
        return None

    tree = lxml.html.fromstring(response.text)

    results_rows = tree.xpath('//div[@class="search-results-row"]')

    data = []
    for row in results_rows:
        row_antibody_code = row.xpath('./div[@class="ab"]/span[@class="ab-code"]/text()')[0]
        row_vendor = ' '.join(row.xpath('./div[@class="ab-company"]/text()')).strip()

        if not vendor_is_correct(row_vendor):
            logger.info('Skipping {}:{} as vendor is wrong'.format(row_vendor, row_antibody_code))
            continue
        if not code_is_correct(row_antibody_code):
            logger.info('Skipping {}:{} as vendor is wrong'.format(row_vendor, row_antibody_code))
            continue

        number_of_citations = int(row.xpath('./div[contains(@class, "ab-cite")]/div[@class="label radius"]/text()')[0])

        data.append(number_of_citations)

    assert data, 'No data parsed for {!r} {!r} but search returned results'.format(antibody_vendor, antibody_code)
    assert len(data) == 1, 'More than one number parsed for {!r} {!r}'.format(antibody_vendor, antibody_code)

    return data[0]

class CitationsTaskBase(luigi.Task):

    _DELAY_ON_TOO_MANY_REQUESTS_SECONDS = 60
    _N_ATTEMPTS = 5

    def requires(self):
        raise NotImplementedError

    def output(self):
        raise NotImplementedError

    def run(self):

        with self.input().open('r') as f:
            df = pd.read_csv(f)

        vendor_catalogues = df[['Antibody Vendor', 'Antibody Catalogue ID']].dropna().drop_duplicates()

        data = []
        for i, (ix, row) in enumerate(vendor_catalogues.iterrows(), start=1):
            vendor = row['Antibody Vendor']
            catalogue_id = row['Antibody Catalogue ID']

            LOGGER.info('Processing antibody {}/{} ({:.2f}%) -- {} {}'.format(i,
                                                                              len(vendor_catalogues),
                                                                              i*100.0 / len(vendor_catalogues),
                                                                              vendor,
                                                                              catalogue_id))
            for x in range(1, self._N_ATTEMPTS+1):
                if x > 1:
                    LOGGER.info('Attempt #{}/{}'.format(x, self._N_ATTEMPTS))

                try:
                    citations = parse_number_of_citations_from_citeab(vendor, catalogue_id)
                    break
                except requests.HTTPError as e:
                    if e.response.status_code == 429:
                        if x < self._N_ATTEMPTS:
                            LOGGER.info('Got 429 sleeping for {} seconds'.format(self._DELAY_ON_TOO_MANY_REQUESTS_SECONDS))
                            time.sleep(self._DELAY_ON_TOO_MANY_REQUESTS_SECONDS)
                            # TODO: change http agent
                            continue
                        else:
                            LOGGER.error('Exhausted number of requests retries')
                            raise
                    else:
                        raise

            data.append({'Antibody Vendor': vendor,
                         'Antibody Catalogue ID': catalogue_id,
                         'Number of Citations': citations})

        ans_df = pd.DataFrame(data)
        ans_df = ans_df[['Antibody Vendor', 'Antibody Catalogue ID', 'Number of Citations']]

        with self.output().open('w') as f:
            ans_df.to_csv(f, index=False)

