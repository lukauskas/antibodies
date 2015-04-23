from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import lxml.html
import requests
import pandas as pd

SEARCH_URI = 'http://www.citeab.com/search?q={}'

LOGGER = logging.getLogger('antibodies.validity.abcite')

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

    response = requests.get(SEARCH_URI.format('{} {}'.format(query_code, query_vendor)))
    response.raise_for_status()
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

        number_of_citations = int(row.xpath('./div[@class="ab-cite"]/div[@class="label radius"]/text()')[0])

        data.append({'Vendor': antibody_vendor,
                     'Antibody Catalogue ID': antibody_code,
                     'Number of Citations': number_of_citations})

    assert data, 'No data parsed for {!r} {!r}'.format(antibody_vendor, antibody_code)

    df = pd.DataFrame(data)
    df = df[['Vendor', 'Antibody Catalogue ID', 'Number of Citations']]

    return df


