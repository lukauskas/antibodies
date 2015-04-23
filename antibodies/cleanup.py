from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import re

KNOWN_VENDORS = {'Cell Signaling Technology', 'Diagenode', 'Millipore',
                 'Abcam', 'Active Motif',
                 #
                 'LP Bio',
                 # Don't know who these guys are, but they're here http://compbio.med.harvard.edu/antibodies/sources/6
                 'Hiroshi Kimura Lab',
                 'Thomas Jenuwein Lab',
                 'Wako',
                 'Abcam',
                 'Epitomics'
                 }


def cleanup_vendor_name(vendor_name):
    if vendor_name in {'Cell Signaling', 'Cell Signalling', 'Cell Signaling Technology', 'CST'}:
        vendor_name = 'Cell Signaling Technology'
    elif vendor_name in {'diagenode', 'Diagenode', 'Diagenonde'}:
        vendor_name = 'Diagenode'
    elif vendor_name in {'none', '0'} or not vendor_name:
        vendor_name = None
    elif vendor_name in {'Millipore', 'Upstate', 'Upstate or Millipore'}:
        vendor_name = 'Millipore'
    elif vendor_name.lower() == 'abcam':
        vendor_name = 'Abcam'

    if vendor_name:
        assert vendor_name in KNOWN_VENDORS, 'Unseen vendor name {!r}'.format(vendor_name)
    return vendor_name


def cleanup_catalog_id(catalog_id, vendor):
    if catalog_id in {'0', 'none'}:
        catalog_id = None
    if catalog_id:
        assert vendor in KNOWN_VENDORS, 'Unseen vendor name {!r}'.format(vendor)

        if catalog_id.endswith('(discontinued)'):
            catalog_id = catalog_id[:-len('(discontinued)')].strip()

        if vendor == 'Abcam':
            if not catalog_id.startswith('EP'):
                catalog_id = catalog_id.lower()
                if not catalog_id.startswith('ab'):
                    catalog_id = 'ab' + str(int(catalog_id))
                if '-' in catalog_id:
                    catalog_id = catalog_id.partition('-')[0]

        elif vendor == 'Cell Signaling Technology':
            # hack for GSM707004
            if catalog_id == 'pAb-037-050,9751S':
                catalog_id = '9751S'

            catalog_id = catalog_id.upper()
            if not catalog_id.endswith('S'):
                if catalog_id.endswith('B'):
                    catalog_id = catalog_id.rstrip('B')  # I think this is a typo

                catalog_id = str(int(catalog_id)) + 'S'

        elif vendor == 'Active Motif':
            if catalog_id.lower().startswith('am'):
                catalog_id = catalog_id[2:]
            catalog_id = int(catalog_id)
        elif vendor == 'Millipore':
            if catalog_id.startswith('Upstate') or catalog_id.startswith('upstate'):
                catalog_id = catalog_id[len('upstate'):].strip()
            elif catalog_id.startswith('Millipore'):
                catalog_id = catalog_id[len('Millipore'):].strip()
            elif catalog_id.startswith('up'):
                catalog_id = catalog_id[len('up'):].strip()
        elif vendor == 'Diagenode':
            if catalog_id.endswith(',9751S'):
                # Hack for GSM707003
                catalog_id = catalog_id[:-len(',9751S')]
            # Hack for GSM941732
            if catalog_id == 'pAb05605':
                catalog_id = 'pAb056050'

            match = re.match('(pab|cs|sn)-?(\d{3})-?(\d{3})$', catalog_id.lower())
            assert match, 'Cannot parse Diagenode catalog id {!r}'.format(catalog_id)
            group_one = match.group(1)
            if group_one == 'pab':
                group_one = 'pAb'
            elif group_one in {'cs', 'sn'}:
                group_one = group_one.upper()

            catalog_id = '{}-{}-{}'.format(group_one, match.group(2), match.group(3))

    return catalog_id

def cleanup_lot_number(vendor, catalog_id, lot_number):

    if lot_number in {'none', 'Unknown', '0', 'unknown'}:
        lot_number = None
    if lot_number:
        assert vendor in KNOWN_VENDORS, 'Unseen vendor name {!r}'.format(vendor)
        if lot_number.lower().startswith('lot'):
            lot_number = lot_number[3:]
            lot_number = lot_number.lstrip('#')

        if vendor == 'Abcam':
            lot_number = int(lot_number.lstrip('u'))
        elif vendor == 'Active Motif':
            lot_number = int(lot_number)
        elif vendor == 'Diagenode':
            # Strip all of the dashes as I cannot figure out how to re-add them from the ones that are stripped
            lot_number = lot_number.replace('-', '')
        elif vendor == 'Millipore':
            # Seems to be least structured
            if lot_number.startswith('dam'):
                # Some dams are lowercased
                lot_number = lot_number.upper()



    return lot_number