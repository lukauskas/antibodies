from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import re

KNOWN_VENDORS = {'Cell Signaling Technology', 'Diagenode', 'Millipore/Upstate',
                 'Abcam', 'Active Motif'}


def cleanup_vendor_name(vendor_name):
    if vendor_name in {'Cell Signaling', 'Cell Signalling', 'Cell Signaling Technology', 'CST'}:
        vendor_name = 'Cell Signaling Technology'
    elif vendor_name in {'diagenode', 'Diagenode', 'Diagenonde'}:
        vendor_name = 'Diagenode'
    elif vendor_name in {'none', '0'} or not vendor_name:
        vendor_name = None
    elif vendor_name in {'Millipore', 'Upstate', 'Upstate or Millipore'}:
        vendor_name = 'Millipore/Upstate'
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
            if catalog_id.lower().startswith('am'):
                catalog_id = catalog_id[2:]
            catalog_id = int(catalog_id)
        elif vendor == 'Millipore/Upstate':
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

            match = re.match('pab-?(\d{3})-?(\d{3})$', catalog_id)
            assert match, 'Canot parse Diagenode catalog id {!r}'.format(catalog_id)
            catalog_id = 'pAb-{}-{}'.format(match.group(1), match.group(2))

    return catalog_id