from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import luigi
import lxml.html
import requests
import pandas as pd
from antibodies.cleanup import cleanup_vendor_name, cleanup_catalog_id, cleanup_lot_number


def parse_antibody_id(antibody_id_in_database):
    antibody_id_in_database = int(antibody_id_in_database)
    url = 'http://compbio.med.harvard.edu/antibodies/antibodies/{}'.format(antibody_id_in_database)

    response = requests.get(url)
    response.raise_for_status()
    tree = lxml.html.fromstring(response.text)

    target = tree.xpath('.//div[@id="main"]/h1/a/text()')[0]
    if target.endswith('(multiple)'):
        target = target[:-len(' (multiple)')]

    info_row = tree.xpath('//table[@class="antibodies"]/tr[2]')[0]
    _extract = lambda result: result[0] if result else None
    source, catalogue_number, lot_number, host, clonality = [_extract(info_row.xpath('./td[{}]//text()'.format(i)))
                                                             for i in xrange(1, 6)]

    source = cleanup_vendor_name(source)
    catalogue_number = cleanup_catalog_id(catalogue_number, source)
    lot_number = cleanup_lot_number(source, catalogue_number, lot_number)

    data = []

    validation_rows = tree.xpath('//table[@class="validations"]/tr')[1:]

    for row in validation_rows:
        category, species, result, notes, conditions, validator = map(_extract, [row.xpath('./td[{}]//text()'.format(i)) for i in [1, 2, 3, 4, 5, 7]])
        image = _extract(row.xpath('./td[6]/a/@href'))

        d = {'Target': target,
             'Vendor': source,
             'Catalogue ID': catalogue_number,
             'Lot Number': lot_number,
             'Host': host,
             'Clonality': clonality,
             'Validation': category,
             'Validation Species': species,
             'Validation Result': result,
             'Validation Notes': notes,
             'Validation Conditions': conditions,
             'Validation Evidence Image': image,
             'Validation Lab': validator}
        data.append(d)

    df = pd.DataFrame(data)
    df = df[['Target', 'Vendor', 'Catalogue ID', 'Lot Number', 'Host', 'Clonality',
             'Validation', 'Validation Species', 'Validation Result', 'Validation Notes', 'Validation Conditions',
             'Validation Evidence Image', 'Validation Lab']]

    return df

def get_total_number_for_db():
    totals_url = 'http://compbio.med.harvard.edu/antibodies/antibodies/'
    totals_response = requests.get(totals_url)
    totals_response.raise_for_status()

    totals_tree = lxml.html.fromstring(totals_response.text)

    return int(totals_tree.xpath('//div[@class="pagination_info"]/b[2]/text()')[0])

class AntibodyValidationDatabase(luigi.Task):

    def output(self):
        return luigi.File('antibody_validation_database.csv')

    @classmethod
    def logger(cls):
        return logging.getLogger('antibodies.validity.antibody_validation_database.AntibodyValidationDatabase')

    def run(self):
        total = get_total_number_for_db()
        df = pd.DataFrame()
        logger = self.logger()
        for i in range(1, total+1):
            logger.info('Processing {}/{} ({:.2f}%)'.format(i, total, i * 100.0 / total))
            df = df.append(parse_antibody_id(i), ignore_index=True)

        logger.info('Done, outputting')

        with self.output().open('w') as f:
            df.to_csv(f, index=False)


if __name__ == '__main__':
    import logging
    AntibodyValidationDatabase.logger().setLevel(logging.INFO)
    logging.basicConfig()
    luigi.run(main_task_cls=AntibodyValidationDatabase)