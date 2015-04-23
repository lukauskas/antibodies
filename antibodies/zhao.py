from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import luigi
from antibodies.validity.citeab import CitationsTaskBase

class ZhaoAntibodies(luigi.Task):

    def output(self):
        return luigi.File('zhao_antibodies.csv')

    def run(self):
        raise Exception('Should not be run, zhao_antibodies.csv is in version control')

class ZhaoCitations(CitationsTaskBase):

    def requires(self):
        return ZhaoAntibodies()

    def output(self):
        return luigi.File('zhao_citations.csv')

if __name__ == '__main__':
    import antibodies.validity.citeab as citeab
    import logging
    citeab.LOGGER.setLevel(logging.INFO)
    logging.basicConfig()
    luigi.run(main_task_cls=ZhaoCitations)