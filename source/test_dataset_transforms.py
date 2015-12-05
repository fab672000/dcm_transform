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
#from datetime import datetime, timedelta
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
    image2 = 'brain2.dcm'
    dataset = None
    input_ds_path = os.path.join(dcm_data_root, image1)
    output_ds_path = os.path.join(dcm_data_root, 'result.dcm')

    timestamp = str(int(time.time()))
    series_uid = '1.2.3.4.' + timestamp + '.0.0.0'
    sopiuid = series_uid + '.' + '42'
    frame_of_ref_uid = '2.3.4.0.' + timestamp + '.0.0.0'
    in_args = None
    test_args = None
    file_count = 0

    def setUp(self):
        self.dataset = dicom.read_file(self.input_ds_path)
        self.in_args = [self.input_ds_path, self.output_ds_path]
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        
    def set_sample_images_io(self, in_filename, out_filename):
        self.input_ds_path = os.path.join(self.dcm_data_root, in_filename)
        self.output_ds_path = os.path.join(self.dcm_data_root, out_filename)
        self.in_args = [self.input_ds_path, self.output_ds_path]

    def instanciate_sut_transform(self, args, file_count=0, out_file=None, in_file=None):
        """Call the sut transform funtion with predef'd parameters for testing purpose """
        if in_file == None:
            in_file = self.input_ds_path
        if out_file == None:
            out_file = self.output_ds_path

        self.file_count, dataset = dcm_transform.transform(file_count, args, '', in_file, out_file)
        return self.file_count, dataset

    
class DcmTestPixelEditor(DcmTestCase):
    def test_draw_pixel(self):
        self.set_sample_images_io(self.image2,'result_pixel.dcm')
        self.in_args.extend(['-pixel', \
            '35', '12', '1023', '1.0', \
            '70', '92', '1023', '0.5'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

    def test_draw_roi(self):
        self.set_sample_images_io(self.image2,'result_roi.dcm')
        self.in_args.extend(['-roi','62','62','4','1023'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

    def test_draw_rect(self):
        self.set_sample_images_io(self.image2,'result_rect.dcm')
        self.in_args.extend(['-rect',\
            '82', '32', '30', '20', '3', '1023', '1.0', \
            '22', '22', '20', '40', '1  ', '1023', '1.0'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

    def test_draw_frect(self):
        self.set_sample_images_io(self.image2,'result_frect.dcm')
        self.in_args.extend(['-frect',\
            '50', '10', '20', '40', '500', '0.5', \
            '80', '80', '30', '30', '500', '.3'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

    def test_draw_elp(self):
        self.set_sample_images_io(self.image2,'result_elp.dcm')
        self.in_args.extend(['-elp','63.5', '63.5', '20' , '20', '1023', '1.0', '1'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

    def test_draw_crosshair(self):
        self.set_sample_images_io(self.image2,'result_crosshair.dcm')
        self.in_args.extend(['-crosshair', \
            '30', '82', '5', '1', '1024', '1.0', '60', '92', '10', '2', '1024', '1.0'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)
        file_count, dataset = self.instanciate_sut_transform(self.test_args)

class DcmTestTagChanges(DcmTestCase):
    """ Test dcm_transform tag changing options"""

    def test_anon1(self):
        """Test anonymization (with dpt) """
        str_anon = 'test_anon1'
        self.in_args.extend(['-an', str_anon, '-dpt'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)

        # needed because dataset.walk(PN_Callback) does not give an optional data param:
        dcm_transform.ARGS = self.test_args
        dcm_transform.anonymize_tags_if_anon(self.dataset, self.test_args, True, True)

    def test_patient_tags(self):
        """Test common patient tags settings"""
        self.in_args.extend(['-pid', '1234', '-pname', 'doe^john', '-dob', '19420402'])
        self.test_args = dcm_transform.parse_arguments(self.in_args)

        file_count, dataset = self.instanciate_sut_transform(self.test_args)

        self.assertEqual(file_count, 1)
        self.assertEqual(dataset.PatientName, 'doe^john')
        self.assertEqual(dataset.PatientID, '1234')
        self.assertEqual(dataset.PatientBirthDate, '19420402')

    def test_generate_uids(self):
        """Test common patient tags settings"""
        self.assertNotEqual(self.series_uid, self.dataset.SeriesInstanceUID)
        self.assertNotEqual(self.sopiuid, self.dataset.SOPInstanceUID)
        self.assertNotEqual(self.frame_of_ref_uid, self.dataset.FrameOfReferenceUID)

        dcm_transform.generate_new_uids(self.dataset, self.series_uid, \
            self.frame_of_ref_uid, self.sopiuid)

        self.assertEqual(self.series_uid, self.dataset.SeriesInstanceUID)
        self.assertEqual(self.sopiuid, self.dataset.SOPInstanceUID)
        self.assertEqual(self.frame_of_ref_uid, self.dataset.FrameOfReferenceUID)

if __name__ == '__main__':
    unittest.main(verbosity=2)
