#!/usr/bin/python

#
# dcm_transform.py
#
"""Read a dicom file (or directory of files), partially "transform" it (them),

  Allows transforming dicom files as follows:
    - Supports anonymisation, changing of various common tags as seriesuid, foruid,
        seriesnum generation and more.
    - Supports 3D transforms (ImagePositionPatient and ImageOrientationPatient)
        by allowing adding position offset
          or rotating the image (see -x,-y,-z,-ax,-ay &-az script parameters).
    - Supports changing any dicom tag value that can be expressed
        as a key, value pair (see -tags parameter).
    - Support recursive tree traversal in order to run in batch mode
        and replicate a complete dicom tree directory
          structure with the required modifications applied to it.

	- Requires Python 2.7
	- Requires Packages: pydicom (as well as scipy, numpy if you don't use anaconda)

    Created by:  fab672000
    Created on:  Sep. 23, 2015

	Original (Git) Repository: https://github.com/fab672000/dcmTransform

	Copyright (c) 2015 fab672000 under MIT License as described below:
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
"""

from __future__ import print_function
import os
import sys
#import platform
import os.path
import argparse
import time
from datetime import datetime, timedelta

import math

import numpy as np
#from scipy import linalg
try:
    import dicom
except ImportError:
    import pydicom as dicom

#------------------------------------------------------------------------------
def parse_arguments():
    """Parse all command line arguments"""
    version = '1.0.0'

    timestamp = str(int(time.time()))
    series_uid = '1.2.3.4.' + timestamp + '.0.0.0'
    frame_of_ref_uid = '3.2.1.0.' + timestamp + '.0.0.0'

    parser = argparse.ArgumentParser(description='dcm_transform, version ' \
                                     + version + ' (https://github.com/fab672000/dcmTransform). ')

    parser.add_argument('-x', nargs='?', type=float, default=0.0, help='X transform offset in mm')
    parser.add_argument('-y', nargs='?', type=float, default=0.0, help='Y transform offset in mm')
    parser.add_argument('-z', nargs='?', type=float, default=0.0, help='Z transform offset in mm')

    parser.add_argument('-ax', nargs='?', type=float, default=0.0, help='AX rotate angle in deg')
    parser.add_argument('-ay', nargs='?', type=float, default=0.0, help='AY rotate angle in deg')
    parser.add_argument('-az', nargs='?', type=float, default=0.0, help='AZ rotate angle in deg')

    parser.add_argument('-sn', nargs='?', type=int, default=66, help='Output Series  Number')
    parser.add_argument('-desc', nargs='?', type=str, default='', \
                        help='Set Custom (Series) Description')
    parser.add_argument('-sdesc', nargs='?', type=str, default='', \
                        help='Set Custom Study Description')
    parser.add_argument('-suid', nargs='?', type=str, default=series_uid, \
                        help='Set Custom Series Instance UID')
    parser.add_argument('-foruid', nargs='?', type=str, default=frame_of_ref_uid, \
                        help='Set Custom Frame Of Reference UID')
    parser.add_argument('-pid', nargs='?', type=str, default='id', \
                        help='Set Custom PatientID when anon is on (see -an)')

    parser.add_argument('-stime', nargs='?', type=str, default='', help='Change study time')
    parser.add_argument('-sdate', nargs='?', type=str, default='', help='Change study date')

    parser.add_argument('-time', nargs='?', type=str, default='', help='Change (series) time')
    parser.add_argument('-date', nargs='?', type=str, default='', help='Change (series) date')

    parser.add_argument('-atime', nargs='?', type=str, default='', help='Change acquisition time')
    parser.add_argument('-adelta', nargs='?', type=str, default='', \
                        help='Change acq. delta in seconds (i.e. 7.78), ' + \
                            'will beapplied sequentially to each file acq time', \
                        metavar=('DELTA_SECS'))
    parser.add_argument('-adate', nargs='?', type=str, default='', help='Change acquisition date')

    parser.add_argument('-tags', nargs='+', type=str, default='', \
                        help='Set your custom tags from  a <tag_name, tag_value> sequence !', \
                        metavar=('TAG_NAME', 'TAG_VALUE'))

    parser.add_argument('-an', nargs='?', type=str, default='', \
                        help='Anonymize series as well (default is none)', metavar=('ANON_NAME'))

    parser.add_argument('-r', '--recurse', action='store_true', \
                        help='traverse input dir tree and reproduce same subtree in output dir')

    parser.add_argument('input_series', help='Input Series folder (or file)  location')
    parser.add_argument('output_series', help='Output Series folder (or file) location')

    #parser.print_help()

    ret_args = parser.parse_args()
    return ret_args

#------------------------------------------------------------------------------
def to_radians(deg):
    """Convert (degrees) to radians"""
    return deg * math.pi / 180.0

#------------------------------------------------------------------------------
def matrix_set_rotation(matrix, xrot, yrot, zrot):
    """Set a rotation into the input matrix from 3 angles in degrees"""
    cos_x = 0.0
    sin_x = 0.0
    cos_y = 0.0
    sin_y = 0.0
    cos_z = 0.0
    sin_z = 0.0

    cos_x = math.cos(to_radians(xrot))
    sin_x = math.sin(to_radians(xrot))
    cos_y = math.cos(to_radians(yrot))
    sin_y = math.sin(to_radians(yrot))
    cos_z = math.cos(to_radians(zrot))
    sin_z = math.sin(to_radians(zrot))

    matrix[0][0] = cos_y * cos_z
    matrix[0][1] = (sin_x * sin_y * cos_z) - (cos_x * sin_z)
    matrix[0][2] = (cos_x * sin_y * cos_z) + (sin_x * sin_z)
    matrix[1][0] = cos_y * sin_z
    matrix[1][1] = (sin_x * sin_y * sin_z) + (cos_x * cos_z)
    matrix[1][2] = (cos_x * sin_y * sin_z) - (sin_x * cos_z)
    matrix[2][0] = -sin_y
    matrix[2][1] = sin_x * cos_y
    matrix[2][2] = cos_x * cos_y

    return matrix

#------------------------------------------------------------------------------
def matrix_set_translation(matrix, translate_x, translate_y, translate_z):
    """Set matrix to translate by a 3D vector."""
    matrix[0, 3] = translate_x
    matrix[1, 3] = translate_y
    matrix[2, 3] = translate_z
    return matrix

#------------------------------------------------------------------------------
def set_str_vec(sarray, val, dim):
    """#set string values for dataset array"""
    for i in range(0, dim):
        sarray[i] = str(val[i])

#------------------------------------------------------------------------------
def compute_3d_transforms(dataset):
    """Compute all 3d transforms (translation/rotation) and update"""

    pos = [dataset.ImagePositionPatient[0].real, \
           dataset.ImagePositionPatient[1].real, dataset.ImagePositionPatient[2].real, 1.]
    row = [dataset.ImageOrientationPatient[0].real, \
           dataset.ImageOrientationPatient[1].real, dataset.ImageOrientationPatient[2].real, 1.]
    col = [dataset.ImageOrientationPatient[3].real, \
           dataset.ImageOrientationPatient[4].real, dataset.ImageOrientationPatient[5].real, 1.]
    xform = np.identity(4)

    matrix_set_rotation(xform, ARGS.ax, ARGS.ay, ARGS.az)

    precision = 5
    new_row = np.around(xform.dot(row), precision)
    new_col = np.around(xform.dot(col), precision)
    new_orient = np.concatenate([new_row[:3], new_col[:3]])

    matrix_set_translation(xform, ARGS.x, ARGS.y, ARGS.z)
    new_pos = np.around([pos[0] + ARGS.x, pos[1] + ARGS.y, pos[2] + ARGS.z], precision)

    #update dataset with new values
    set_str_vec(dataset.ImagePositionPatient, new_pos[:3], 3)
    set_str_vec(dataset.ImageOrientationPatient, new_orient[:6], 6)

#------------------------------------------------------------------------------
def fmt_float(title, num, sep):
    """format all float numbers the same way but only if not zero"""
    try:
        if num != 0.0:
            val = title + "{:+06.1f}".format(num) + sep
        else:
            val = ''
    except Exception:
        val = '???' + sep

    return val

#------------------------------------------------------------------------------
def fmt_float3d(prefix, num1, num2, num3, sep):
    """Format a float 3d set for display"""
    try:
        val = fmt_float(prefix + 'X', num1, sep) + fmt_float(prefix + 'Y', num2, sep) \
            + fmt_float(prefix + 'Z', num3, sep)
    except Exception:
        val = "???"
    return val

#------------------------------------------------------------------------------
def modify_time(origdate, origtime, deltasecs):
    """Takes a dicom time stamp in str format, and a delta and then adds the delta
       to riginal time and returns that modified time
       returns a datetime result
    """
    try:
        torig = datetime.strptime(origdate + ' ' + origtime, "%Y%m%d %H%M%S.%f")
        dmicrosec = int((float(deltasecs) - int(float(deltasecs))) * 1000000)
        tdelta = timedelta(seconds=int(float(deltasecs)), microseconds=dmicrosec)

        tresult = torig + tdelta

        return tresult

    except Exception:
        print("modify_time: could not modify timestamp error during conversion ...")

#------------------------------------------------------------------------------
def get_dicom_date_from(dtime):
    """Format the current date inside the input datetime to a dicom str format"""
    try:
        res = dtime.strftime("%Y%m%d")
        return res
    except Exception:
        print("get_dicom_date_from: could not convert")
        return '19010101'

#------------------------------------------------------------------------------
def get_dicom_time_from(dtime):
    """Format the current time inside the input datetime to a dicom str format"""
    try:
        res = dtime.strftime("%H%M%S.%f")
        return res

    except Exception:
        print("get_dicom_time_from: could not convert")
        return '000000.000000'

#------------------------------------------------------------------------------
def transform(file_count, desc_prefix, input_filename, output_filename, \
              remove_curves=False, remove_private_tags=False):
    """Replace data element values to partly transform a DICOM file.
    Note: completely transforming a DICOM file is very complicated; there
    are many things this example code does not address. USE AT YOUR OWN RISK.
    """

    # Define call-back functions for the dataset.walk() function
    def PN_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements."""
        if ARGS.an != '':
            if data_element.VR == "PN":
                data_element.value = ARGS.an

            #print (data_element.value)

    def curves_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements."""
        if data_element.tag.group & 0xFF00 == 0x5000:
            del ds[data_element.tag]

    def transform_dates(file_count, dataset):
        """Transform acquisition date or time, or series date or time"""
        try:
            if ARGS.adate != '':
                dataset.AcquisitionDate = ARGS.adate
        except Exception as exc:
            print(exc)

        try:
            if ARGS.atime != '':
                dataset.AcquisitionTime = ARGS.atime
        except Exception as exc:
            print(exc)

        try:
            if ARGS.adelta != '':
                secs_offset = str(float(ARGS.adelta) * file_count)
                current_time = modify_time(dataset.AcquisitionDate, dataset.AcquisitionTime, \
                                           secs_offset)
                dataset.AcquisitionDate = get_dicom_date_from(current_time)
                dataset.AcquisitionTime = get_dicom_time_from(current_time)
        except Exception as exc:
            print(exc)

        # Series Date & Time
        try:
            if ARGS.date != '':
                dataset.SeriesDate = ARGS.date
        except Exception as exc:
            print(exc)

        try:
            if ARGS.time != '':
                dataset.SeriesTime = ARGS.time
        except Exception as exc:
            print(exc)

        # Study Date & Time
        try:
            if ARGS.stime != '':
                dataset.StudyTime = ARGS.stime
        except Exception as exc:
            print(exc)

        try:
            if ARGS.sdate != '':
                dataset.StudyDate = ARGS.sdate
        except Exception as exc:
            print(exc)

    try:
        file_count += 1
        # Load the current dicom file to 'transform'
        dataset = dicom.read_file(input_filename)

        compute_3d_transforms(dataset)

        # assign custom tags
        if ARGS.tags != '':
            try:
                tags_len = len(ARGS.tags)
                if tags_len > 1:
                    print()
                for i in xrange(0, tags_len / 2 * 2, 2):
                    print("Setting tag " + ARGS.tags[i] + " to " + ARGS.tags[i + 1] + ' ... ')
                    data_element = dataset.data_element(ARGS.tags[i])
                    if data_element == None:
                        print("  Could not find that tag, value won't be set ...")
                    else:
                        try:
                            data_element.value = ARGS.tags[i + 1]
                        except Exception as exc:
                            print('  Could not set that tag value  <' + ARGS.tags[i + 1] + \
                                  '>, value will not be set ...')

                if tags_len % 2 != 0:
                    print("  Warning: list of pair of <tags value> expected, "+ \
                          "but odd count was found instead, " + \
                          "found ending: <" + ARGS.tags[tags_len - 1] + '>')
            except Exception as exc:
                print(exc)

        # Deal with all sorts of dates and time if user asks for it:
        transform_dates(file_count, dataset)

        # Remove patient name and any other person names
        dataset.walk(PN_callback)

        dataset.SeriesInstanceUID = ARGS.suid
        dataset.SOPInstanceUID = ARGS.suid + '.' + str(dataset.InstanceNumber)
        dataset.FrameOfReferenceUID = ARGS.foruid

        dataset.SeriesNumber = ARGS.sn

        # Series Description gets automatically populated with useful transform
        # info by default:
        if ARGS.desc != '':
            try:
                series_description = desc_prefix + dataset.SeriesDescription
                if len(series_description) > 60:
                    series_description = (series_description[:60] + '..')
                else:
                    series_description = ARGS.desc
                dataset.SeriesDescription = series_description
            except Exception as exc:
                print(exc)

        #optionally change study desc
        if ARGS.sdesc != '':
            try:
                dataset.StudyDescription = ARGS.sdesc
            except Exception as exc:
                print(exc)


        if ARGS.an != '':
            dataset.PatientID = ARGS.pid
            # Remove data elements (should only do so if DICOM type 3 optional)
            # Use general loop so easy to add more later
            # Could also have done: del ds.OtherPatientIDs, etc.
            for name in ['OtherPatientIDs', 'OtherPatientIDsSequence']:
                if name in dataset:
                    delattr(dataset, name)

            # Same as above but for blanking data elements that are type 2.
            for name in ['PatientBirthDate', 'StudyDate', 'SeriesDate']:
                if name in dataset:
                    dataset.data_element(name).value = '19010101'

            for name in ['SeriesTime', 'StudyTime']:
                if name in dataset:
                    dataset.data_element(name).value = '000000.000000'

            ## Remove private tags if function argument says to do so.  Same
            ## for curves
            if remove_private_tags:
                dataset.remove_private_tags()
            if remove_curves:
                dataset.walk(curves_callback)

        # write the 'transformed' DICOM out under the new filename
        dataset.save_as(output_filename)

    except Exception as exc:
        print(exc)

    return file_count

#------------------------------------------------------------------------------
def iterate_once(in_args, input_dir, output_dir):
    """Execute the full script except the recursive option"""
    series_file_count = 0

    try:
        series_desc_prefix = 'T[' + fmt_float3d('', in_args.x, in_args.y, in_args.z, ' ') \
            + fmt_float3d('A', in_args.ax, in_args.ay, in_args.az, ' ') + '] '

    except Exception:
        print("Could not convert the x, y,z offsets")
        sys.exit()

    if os.path.isdir(input_dir):
        if os.path.exists(output_dir):
            if not os.path.isdir(output_dir):
                raise IOError("Input is directory; output name exists but is not a directory")
        else:  # out_dir does not exist; create it.
            os.makedirs(output_dir)

        fnames = os.listdir(input_dir)
        for filename in fnames:
            if not os.path.isdir(os.path.join(input_dir, filename)):
                print('Transforming ' + series_desc_prefix + filename + " ...", end='')
                series_file_count = transform(series_file_count, series_desc_prefix, \
                                       os.path.join(input_dir, filename), \
                                       os.path.join(output_dir, filename))
                print(" done\r")
    else:  # first arg not a directory, assume two files given
        in_filename = in_args.input_series
        out_filename = in_args.output_series
        transform(series_file_count, series_desc_prefix, in_filename, out_filename)
    print()

#------------------------------------------------------------------------------
# main program
#------------------------------------------------------------------------------
if __name__ == "__main__":
    ARGS = parse_arguments()

    #for timestamped offset computing
    if not ARGS.recurse:
        iterate_once(ARGS, ARGS.input_series, ARGS.output_series)
    else:
        IN_DIR = ARGS.input_series
        OUT_DIR = ARGS.output_series
        for dirpath, dirnames, filenames in os.walk(IN_DIR):
            cur_dir = os.path.join(OUT_DIR, dirpath[1 + len(IN_DIR):])
            print('Traversing: ', dirpath)
            iterate_once(ARGS, dirpath, cur_dir)
