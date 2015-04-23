from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

KNOWN_VENDORS = {'Cell Signaling Technology', 'Diagenode', 'Millipore/Upstate',
                 'Abcam', 'Active Motif'}


def cleanup_vendor_name(vendor_name):
    if vendor_name in {'Cell Signaling', 'Cell Signalling', 'Cell Signaling Technology', 'CST'}:
        vendor_name = 'Cell Signaling Technology'
    elif vendor_name in {'diagenode', 'Diagenode', 'Diagenonde'}:
        vendor_name = 'Diagenode'
    elif vendor_name in {'none', '0'}:
        vendor_name = None
    elif vendor_name in {'Millipore', 'Upstate', 'Upstate or Millipore'}:
        vendor_name = 'Millipore/Upstate'

    assert vendor_name in KNOWN_VENDORS, 'Unseen vendor name {!r}'.format(vendor_name)
    return vendor_name


def cleanup_catalog_id(catalog_id, vendor):
    assert vendor in KNOWN_VENDORS

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
            if catalog_id.lower().startswith('am'):
                catalog_id = catalog_id[2:]
            catalog_id = int(catalog_id)
        elif vendor == 'Millipore/Upstate':
            if catalog_id.startswith('Upstate') or catalog_id.startswith('upstate'):
                catalog_id = catalog_id[len('upstate'):].strip()
            if catalog_id.startswith('Millipore'):
                catalog_id = catalog_id[len('Millipore'):].strip()
    return catalog_id