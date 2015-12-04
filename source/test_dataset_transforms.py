#!/usr/bin/python
#
# test_dataset_transforms.py
#
""" Shall Test the dcm tags manipulation
"""

import unittest
import source.dcm_transform

#import os, os.path, sys, math, argparse, time
#from datetime import datetime, timedelta
#import numpy as np
#from scipy import linalg

try:
    import dicom
except ImportError:
    import pydicom as dicom

class DcmTestCase(unittest.TestCase):
    dataset = None

    def load_dataset(self, input_filename):
        self.dataset = dicom.read_file(input_filename)

    def setUp(self):
        self.dataset = self.load_dataset(self, 'examples\\roi\\')

class Test_test_dataset_transforms(DcmTestCase):

    def test_A(self):
        self.fail("Not implemented")


if __name__ == '__main__':
    unittest.main()