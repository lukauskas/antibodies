from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import requests
import re
import lxml.html
import json
import pandas as pd

ABCAM_SEARCH_REQUEST = 'http://www.abcam.com/products?keywords={}'
ABCAM_PRICES_REQUEST = 'http://www.abcam.com/DatasheetProperties/Availability?abId={}'

CELL_SIGNALING_PRICES_REQUEST = 'http://www.neb.uk.com/product/{}'
ACTIVE_MOTIF_PRICES_REQUEST = 'http://www.activemotif.com/search?terms={}'

COUNTRY_CODE = 'gb'

def parse_abcam_antibody(antibody_id):

    session = requests.Session()

    abcam_response = requests.get(ABCAM_SEARCH_REQUEST.format(antibody_id))
    abcam_response.raise_for_status()

    tree = lxml.html.fromstring(abcam_response.text)
    product_id = tree.xpath('//input[@name="ProductCode"]/@value')[0]
    assert antibody_id == product_id


    # Parse the request ID for XHR request (otherwise it won't work)
    xpid = re.search('loader_config=\{xpid:"([^"]+)"\}', abcam_response.text)
    xpid = xpid.group(1)

    # Strip ab from product id
    integer_product_id = int(product_id.lstrip('ab'))

    # Fetch the price data
    prices_response = session.get(ABCAM_PRICES_REQUEST.format(integer_product_id),
                                  headers={'X-NewRelic-ID': xpid, 'X-Requested-With': 'XMLHttpRequest'},
                                  cookies={'C2LC': 'GB'})

    prices_response.raise_for_status()

    json_prices = json.loads(prices_response.text)

    assert json_prices['size-information']['CountryCode'] == COUNTRY_CODE.lower()

    data = []
    for size_info in json_prices['size-information']['Sizes']:
        size = size_info['Size']
        price = size_info['Price']

        if size.endswith('&micro;g'):
            weight_in_ug = float(size.partition('&micro;g')[0].strip())
        else:
            raise Exception('Don\'t know how to get weight in micrograms for {!r}'.format(size))

        currency = None
        if price.startswith('&pound;'):
            price = float(price.partition('&pound;')[2].strip())
            currency = 'GBP'
        else:
            raise Exception('Don\'t know how to parse price for {!r}'.format(size))

        data.append({'antibody_id': product_id,
                     'price': price, 'currency': currency,
                     'amount': weight_in_ug,
                     'unit': 'ug',
                     'vendor': 'Abcam'})

    df = pd.DataFrame(data)
    return df

def parse_cell_signalling_antibody(cell_signalling_antibody_id):

    integer_id = int(cell_signalling_antibody_id.strip('s'))
    session = requests.Session()

    response = session.get(CELL_SIGNALING_PRICES_REQUEST.format(integer_id))
    response.raise_for_status()

    tree = lxml.html.fromstring(response.text)

    price_header = tree.xpath('//table[@id="purchase"]//th[@id="th-price"]/text()')[0]
    assert u'\xa3' in price_header  # Assert pound symbol is in the header

    price_trs = tree.xpath('//table[@id="purchase"]/tbody/tr')

    data = []

    for price_tr in price_trs:
        # If we cannot parse the size-tests column we must be in one of the custom rows
        amount_str = price_tr.xpath('./td[@class="size-tests"]/text()')
        if len(amount_str) == 0:
            continue
        else:
            amount_str = amount_str[0]

        price_in_gbp = float(price_tr.xpath('./td[@class="price"]/text()')[0])
        currency = 'GBP'

        amount, __, __ = amount_str.partition(' ')
        amount = float(amount)

        if u'\xb5l' in amount_str:
            units = 'ul'
        else:
            raise Exception('Cannot parse units from {!r}'.format(amount_str))

        data.append({'antibody_id': cell_signalling_antibody_id.upper(),
                     'price': price_in_gbp, 'currency': currency,
                     'amount': amount,
                     'unit': units,
                     'vendor': 'Cell Signaling Technology'
                     })
    assert data

    return data

def parse_active_motif(antibody_id):
    if COUNTRY_CODE == 'gb':
        country_id = 'UK'
    else:
        raise Exception('Unsupported country code {!r}'.format(COUNTRY_CODE))
    integer_id = int(antibody_id.lstrip('am'))

    s = requests.Session()
    response = s.get(ACTIVE_MOTIF_PRICES_REQUEST.format(integer_id),
                     cookies=dict(country=country_id, language='en'))

    tree = lxml.html.fromstring(response.text)

    data = []
    rows = tree.xpath('//div[@id="title"]//table[@class="app-table"]//tr')
    for row in rows:

        amount_str = row.xpath('./td[2]/text()')[0].replace(u'\xa0', '')
        amount, __, units = amount_str.partition(' ')

        if u'\xb5g' in units:
            units = 'ug'
        elif u'\xb5l' in units:
            units = 'ul'
        else:
            raise Exception('Cannot parse units: {!r}'.format(units))

        amount = float(amount)
        price = row.xpath('./td[3]/strong/text()')[0]
        if u'\xa3' in price:
            currency = 'GBP'
            price = float(price.lstrip(u'\xa3'))
        else:
            raise Exception('Cannot parse currency from {!r}'.format(price))

        data.append({'antibody_id': antibody_id.upper(),
                     'price': price, 'currency': currency,
                     'amount': amount,
                     'unit': units,
                     'vendor': 'Active Motif'
                     })
    assert data
    return data



def parse_antibody(antibody_id):

    if antibody_id.startswith('ab'):
        return parse_abcam_antibody(antibody_id)
    elif antibody_id.endswith('s'):
        return parse_cell_signalling_antibody(antibody_id)
    elif antibody_id.startswith('am'):
        return parse_active_motif(antibody_id)
    else:
        raise Exception('Cannot determine antibody vendor from {!r}'.format(antibody_id))

if __name__ == '__main__':
    import sys
    product_ids = sys.argv[1:]
    product_ids = map(lambda x: x.lower(), product_ids)
    assert len(set(product_ids)) == len(product_ids)

    main_df = pd.DataFrame()

    for product_id in product_ids:
        df = parse_antibody(product_id)
        main_df = main_df.append(df, ignore_index=True)

     # Reorder columns
    main_df = main_df[['vendor', 'antibody_id', 'amount', 'unit', 'price', 'currency']]
    main_df.to_csv(sys.stdout, index=False)
