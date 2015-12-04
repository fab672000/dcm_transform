#!/usr/bin/python
#
# test_dataset_transforms.py
#
""" Shall Test the dcm tags manipulation
"""

import unittest
import dcm_transform

import os, os.path, time
#, sys, math, argparse
from datetime import datetime, timedelta
#import numpy as np
#from scipy import linalg

try:
    import dicom
except ImportError:
    import pydicom as dicom

class DcmTestCase(unittest.TestCase):
    """Base class for testing dcm operations, loads a dicom test dataset by default"""

    dcm_data_root = 'examples\\data\\'
    image1 = 'brain1.dcm'
    dataset = None
    
    input_ds_path = os.path.join(dcm_data_root, image1)
    output_ds_path = os.path.join(dcm_data_root, 'result.dcm')

    timestamp = str(int(time.time()))
    series_uid = '1.2.3.4.' + timestamp + '.0.0.0'
    sopiuid = series_uid + '.' + '42'
    frame_of_ref_uid = '2.3.4.0.' + timestamp + '.0.0.0'

    def setUp(self):
        self.dataset = dicom.read_file(self.input_ds_path)

class DcmTestTagChanges(DcmTestCase):
    """ test dcm_transform tag changing options"""  
    
    def test_A(self):
        self.assertNotEqual(self.series_uid, self.dataset.SeriesInstanceUID)
        self.assertNotEqual(self.sopiuid, self.dataset.SOPInstanceUID)
        self.assertNotEqual(self.frame_of_ref_uid, self.dataset.FrameOfReferenceUID)

        dcm_transform.generate_new_uids(self.dataset, self.series_uid, self.frame_of_ref_uid, self.sopiuid)

        self.assertEqual(self.series_uid, self.dataset.SeriesInstanceUID)
        self.assertEqual(self.sopiuid, self.dataset.SOPInstanceUID)
        self.assertEqual(self.frame_of_ref_uid, self.dataset.FrameOfReferenceUID)

if __name__ == '__main__':
    unittest.main(verbosity=2)