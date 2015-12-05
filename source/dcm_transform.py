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
          or rotating the image (see -x, -y, -z, -ax, -ay &-az script parameters).
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
import os, sys, math, argparse, time
import os.path
from datetime import datetime, timedelta

import numpy as np
#from scipy import linalg
try:
    import dicom
except ImportError:
    import pydicom as dicom

#------------------------------------------------------------------------------
class PixelEditor():
    """ Pixel editing utility class for drawing simple geometries inside the 2D pixel array"""
    pixel_buffer = None

    #------------------------------------------------------------------------------
    def __init__(self, the_buf): 
        self.pixel_buffer = the_buf
    #------------------------------------------------------------------------------
    def buffer_length(self): 
        return len(self.pixel_buffer)
    #------------------------------------------------------------------------------
    def buffer_to_string(self): 
        return self.pixel_buffer.tostring()
    #------------------------------------------------------------------------------
    def draw_pixel(self, pos_x, width, xstep, pos_y, height, ystep, val, alpha=1.0):
        """Set a pixel buffer value val at x, y to x+w, y+h with an xstep and ystep increments """
        for col in range(int(pos_x), int(pos_x + width), xstep):
            for row in range(int(pos_y), int(pos_y + height), ystep):
    #            print(col, row)
                orig = self.pixel_buffer[row, col] * (1 - alpha)
                new = val * alpha
                blend = orig + new
                self.pixel_buffer[row, col] = blend
    #------------------------------------------------------------------------------
    def draw_hline(self, pos_x, pos_y, width, step, val, alpha=1.0):
        """Draw an horizontal line starting at x, y of width w with step step and value val """
        self.draw_pixel(pos_x, width, step, pos_y, 1, 1, val, alpha)
    #------------------------------------------------------------------------------
    def draw_vline(self, pos_x, pos_y, height, step, val, alpha=1.0):
        """Set a pixel buffer value val at x, y to x+w, y+h with an xstep and ystep increments """
        self.draw_pixel(pos_x, 1, 1, pos_y, height, step, val, alpha)
    #------------------------------------------------------------------------------
    def draw_elp(self, pos_x, pos_y, width, height, val, alpha=1.0, step=1):
        """ Draws an ellipse"""
        for theta in range(0, 360, step):
            point_x1 = pos_x + width / 2.0 * math.cos(theta / 180.0 * math.pi)
            point_y1 = pos_y - height / 2.0 * math.sin(theta / 180.0 * math.pi)
            self.draw_pixel(int(point_x1), 1, 1, int(point_y1), 1, 1, val, alpha)
    #------------------------------------------------------------------------------
    def draw_rect(self, pos_x, pos_y, width, height, step, val, alpha=1.0):
        """Draws a rect in the buffer at x, y to x+w, y+h with an xstep and ystep increments """
        self.draw_hline(pos_x, pos_y, width, step, val, alpha)
        self.draw_hline(pos_x, pos_y + height - 1, width, step, val, alpha)
        self.draw_vline(pos_x, pos_y, height, step, val, alpha)
        self.draw_vline(pos_x + width - 1, pos_y, height, step, val, alpha)
    #------------------------------------------------------------------------------
    def draw_frect(self, pos_x, pos_y, width, height, val, alpha=1.0):
        """Draws a rect in the buffer at x, y to x+w, y+h with an [0, ] transparency factor """
        for i in range(int(pos_y), int(pos_y + height)):
            self.draw_hline( pos_x, i, width, 1, val, alpha)
    #------------------------------------------------------------------------------
    def draw_xhair(self, pos_x, pos_y, pen_size, width, val, alpha=1.0):
        """ Draws a cross hair  in the buffer at x, y with size s and pen width w
            with intensity val and alpha blending alpha"""
        rounded_width = int(width / 2 * 2) + 1
        half_width = (rounded_width - 1) / 2
        line_len = int(pen_size)
        self.draw_frect(pos_x - line_len, pos_y - half_width, line_len - half_width, rounded_width, val, alpha)
        self.draw_frect(pos_x + half_width + 1, pos_y - half_width, line_len - half_width, rounded_width, val, alpha)
        self.draw_frect(pos_x - half_width, pos_y - line_len, rounded_width, line_len - half_width, val, alpha)
        self.draw_frect(pos_x - half_width, pos_y + 1 + half_width, rounded_width, line_len - half_width, val, alpha)

#------------------------------------------------------------------------------
def parse_arguments(the_args=None):
    """Parse all command line arguments"""
    version = '1.1.9'

    timestamp = str(int(time.time()))
    defaulf_series_uid = '1.2.3.4.' + timestamp + '.0.0.0'
    defaulf_frame_of_ref_uid = '2.3.4.0.' + timestamp + '.0.0.0'

    parser = argparse.ArgumentParser(description='dcm_transform, version ' \
                                     + version + ' (https://github.com/fab672000/dcmTransform). ')

    parser.add_argument('input_series', help='Input Series folder (or file)  location')
    parser.add_argument('output_series', help='Output Series folder (or file) location')

    parser.add_argument('-an', nargs='?', type=str, default='', \
                        help='Anonymize series as well (default is none)', metavar=('ANON_NAME'))

    parser.add_argument('-r', '--recurse', action='store_true', \
                        help='traverse input dir tree and reproduce same subtree in output dir')

    parser.add_argument('-dpt', '--delete_private_tags', action='store_true', \
                        help='Delete private tags. Can be useful when anonymizing.')


    parser.add_argument('-x', nargs='?', type=float, default=0.0, help='X transform offset in mm')
    parser.add_argument('-y', nargs='?', type=float, default=0.0, help='Y transform offset in mm')
    parser.add_argument('-z', nargs='?', type=float, default=0.0, help='Z transform offset in mm')

    parser.add_argument('-ax', nargs='?', type=float, default=0.0, help='AX rotate angle in deg')
    parser.add_argument('-ay', nargs='?', type=float, default=0.0, help='AY rotate angle in deg')
    parser.add_argument('-az', nargs='?', type=float, default=0.0, help='AZ rotate angle in deg')

    parser.add_argument('-sn', nargs='?', type=int, default=-1, help='Output Series  Number')
    parser.add_argument('-desc', nargs='?', type=str, default='', \
                        help='Set Custom (Series) Description')
    parser.add_argument('-sdesc', nargs='?', type=str, default='', \
                        help='Set Custom Study Description')
    parser.add_argument('-suid', nargs='?', type=str, default=defaulf_series_uid, \
                        help='Set Custom Series Instance UID')
    parser.add_argument('-foruid', nargs='?', type=str, default=defaulf_frame_of_ref_uid, \
                        help='Set Custom Frame Of Reference UID')
    parser.add_argument('-pid', nargs='?', type=str, default='', \
                        help='Set Custom PatientID')
    parser.add_argument('-pname', nargs='?', type=str, default='', \
                        help='Set Custom PatientName')
    parser.add_argument('-dob', nargs='?', type=str, default='', \
                        help='Set Custom (Patient) Date of Birth')
    parser.add_argument('-iname', nargs='?', type=str, default='', \
                        help='Set InstitutionName')
    parser.add_argument('-iaddr', nargs='?', type=str, default='', \
                        help='Set InstitutionAddress')
    parser.add_argument('-proto', nargs='?', type=str, default='', \
                        help='Set ProtocolName')
    parser.add_argument('-mname', nargs='?', type=str, default='', \
                        help='Set ManufacturerName')
    parser.add_argument('-mmname', nargs='?', type=str, default='', \
                        help='Set ManufacturerModelName')

    parser.add_argument('-time', nargs='?', type=str, default='', help='Change (series) time')
    parser.add_argument('-date', nargs='?', type=str, default='', help='Change (series) date')
    parser.add_argument('-stime', nargs='?', type=str, default='', help='Change study time')
    parser.add_argument('-sdate', nargs='?', type=str, default='', help='Change study date')
    parser.add_argument('-ctime', nargs='?', type=str, default='', help='Change content time')
    parser.add_argument('-cdate', nargs='?', type=str, default='', help='Change content date')


    parser.add_argument('-atime', nargs='?', type=str, default='', help='Change acquisition time')
    parser.add_argument('-adelta', nargs='?', type=str, default='', \
        help='Change acq. delta in seconds (i.e. 7.78), ' + \
            'will beapplied sequentially to each file acq time', \
        metavar=('DELTA_SECS'))
    parser.add_argument('-adate', nargs='?', type=str, default='', help='Change acquisition date')

    parser.add_argument('-tags', nargs='+', type=str, default='', \
        help='Set your custom tags from  a <tag_name, tag_value> sequence !', \
        metavar=('TAG_NAME', 'TAG_VALUE'))

    parser.add_argument('-pixel', nargs='+', type=str, default='', \
        help='Set one pixel value in each file ', \
        metavar=('ROW_COL_VAL_1', 'ROW_COL_VAL_2'))

    parser.add_argument('-roi', nargs='+', type=str, default='', \
        help='Set a squared ROI of width LEN starting at top-left position x, y ' + \
            'of pixel value val.', \
        metavar=('X, Y, LEN, VAL', '...'))

    parser.add_argument('-crosshair', nargs='+', type=str, default='', \
        help='Set a crosshair at top-left pos x, y of size S and line width W(odd number)' + \
            ' with an intensity I and alpha blending A.', \
        metavar=('X, Y, S, W, I, A', '...'))

    parser.add_argument('-elp', nargs='+', type=str, default='', \
        help='Set an ellipse at top-left pos x, y of width radius w and radius h' + \
            ' with intensity I and alpha blending A and cos/sin steps S(i.e. 120) .', \
        metavar=('X, Y, W, H, I, A, S', '...'))

    parser.add_argument('-rect', nargs='+', type=str, default='', \
        help='Set a rectangle at top-left pos x, y of width w and height h and step S ' + \
            'with intensity I and alpha blending A.', \
        metavar=('X, Y, W, H, S, I, A', '...'))

    parser.add_argument('-frect', nargs='+', type=str, default='', \
        help='Set a filled rectangle at top-left pos x, y of width w and height h, S ' + \
            'with intensity I and alpha blending A.', \
        metavar=('X, Y, W, H, I, A', '...'))

    #parser.print_help()

    ret_args = parser.parse_args(the_args)
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
def generate_soiud_from_seriesuid(defaulf_series_uid, instance_number):
    return defaulf_series_uid + '.' + str(instance_number)

#------------------------------------------------------------------------------
def generate_new_uids(dataset, suid, foruid, sopiuid):
    """Generate new series, FOR, SOP Instance uids"""
    #Generate new UIDs automatically when any transform changes the geometry
    dataset.SeriesInstanceUID = suid
    dataset.SOPInstanceUID = sopiuid
    dataset.file_meta.data_element("MediaStorageSOPInstanceUID").value = sopiuid
    # do not use: dataset.MediaStorageSOPInstanceUID = dataset.SOPInstanceUID!

    # Generate a new FORUID too
    dataset.FrameOfReferenceUID = foruid

#------------------------------------------------------------------------------
def compute_3d_transforms(dataset, args):
    """Compute all 3d transforms (translation/rotation) and update"""

    if not is_3d_tranformation(args):
        return

    #Generate new UIDs automatically when any transform changes the geometry
    generate_new_uids(dataset, args.suid, args.foruid, generate_soiud_from_seriesuid(defaulf_series_uid, dataset.InstanceNumber))

    pos = [dataset.ImagePositionPatient[0].real, \
           dataset.ImagePositionPatient[1].real, dataset.ImagePositionPatient[2].real, 1.]
    row = [dataset.ImageOrientationPatient[0].real, \
           dataset.ImageOrientationPatient[1].real, dataset.ImageOrientationPatient[2].real, 1.]
    col = [dataset.ImageOrientationPatient[3].real, \
           dataset.ImageOrientationPatient[4].real, dataset.ImageOrientationPatient[5].real, 1.]
    xform = np.identity(4)

    matrix_set_rotation(xform, args.ax, args.ay, args.az)

    precision = 5
    new_row = np.around(xform.dot(row), precision)
    new_col = np.around(xform.dot(col), precision)
    new_orient = np.concatenate([new_row[:3], new_col[:3]])

    matrix_set_translation(xform, args.x, args.y, args.z)
    new_pos = np.around([pos[0] + args.x, pos[1] + args.y, pos[2] + args.z], precision)

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
# Define call-back functions for the dataset.walk() function
def pn_callback(dataset, data_element):
    """Called from the dataset "walk" recursive function for all data elements."""
    if data_element.VR == "PN":
        data_element.value = ARGS.an
    #print (data_element.value)
#------------------------------------------------------------------------------
def curves_callback(dataset, data_element):
    """Called from the dataset "walk" recursive function for all data elements."""
    if data_element.tag.group & 0xFF00 == 0x5000:
        del dataset[data_element.tag]

#------------------------------------------------------------------------------
def assign_custom_tags(dataset, args):
    """ Given a dataset and the args key, val pair array  arguments set custom tags"""
    # assign custom tags
    if args == '':
        return
    try:
        tags_len = len(args)
        if tags_len > 1:
            print()
        for i in xrange(0, tags_len / 2 * 2, 2):
            print("Setting tag " + args[i] + " to " + args[i + 1] + ' ... ')
            try:
                data_element = dataset.data_element(args[i])
            except KeyError as exc: #OK, maybe user wants a tag that is in file_meta metadata structure
                                    #so give it a second chance:
                data_element = dataset.file_meta.data_element(args[i])

            if data_element == None:
                print("  Could not find that tag, value won't be set ...")
            else:
                try:
                    data_element.value = args[i + 1]
                except Exception as exc:
                    print('  Could not set that tag value  <' + args[i + 1] + \
                            '>, value will not be set ...')

        if tags_len % 2 != 0:
            print("  Warning: list of pair of <tags value> expected, " + \
                    "but odd count was found instead, " + \
                    "found ending: <" + args[tags_len - 1] + '>')
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def set_image_pixels(dataset, args):
    """ Given a dataset and the args key, val pair array  arguments set custom tags"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        pix_len = len(args)
        n_vals = 4
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = int(args[i + 0])
            pos_y = int(args[i + 1])
            val = int(args[i + 2])
            alpha = float(args[i + 3])
            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_pixel( pos_x, 1, 1, pos_y, 1, 1, val, alpha)
                except Exception as exc:
                    print('  Could not set that pixel value  <' + args[i + 2] + \
                            '>, value will not be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of quadruplets expected, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def draw_crosshair(dataset, args):
    """ Draws a rectangular region of interest"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        n_vals = 6
        pix_len = len(args)
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = float(args[i + 0]) # x pos
            pos_y = float(args[i + 1]) # y pos
            crosshair_size = int(args[i + 2]) # crosshair size (line length)
            pen_width = int(args[i + 3]) # pen width (number of pixels)
            intensity = int(args[i + 4]) # intensity (for dash lines purpose)
            alpha = float(args[i + 5]) # alpha blending (for dash lines purpose)

            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_xhair(pos_x, pos_y, crosshair_size, pen_width, intensity, alpha)
                except Exception as exc:
                    print('  Error while trying to draw the crosshair, pixel values wont  be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of 6 parameters sequences, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def draw_roi(dataset, args):
    """ Draws a rectangular region of interest"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        n_vals = 4
        pix_len = len(args)
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = float(args[i + 0])
            pos_y = float(args[i + 1])
            stepping = float(args[i + 2])
            val = int(args[i + 3])
            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_rect(pos_x, pos_y, stepping, stepping, 1, val)
                except Exception as exc:
                    print('  Could not set that pixel value  <' + args[i + 2] + \
                            '>, value will not be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of triplets expected, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def draw_ellipse(dataset, args):
    """ Draws a rectangular region of interest"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        n_vals = 7
        pix_len = len(args)
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = float(args[i + 0]) # x pos
            pos_y = float(args[i + 1]) # y pos
            width = int(args[i + 2]) # rect rad width
            height = int(args[i + 3]) # rect rad height
            pixel_intensity = int(args[i + 4]) # pixel intensity
            alpha = float(args[i + 5]) # pixel alpha transparency
            stepping = int(args[i + 6]) # stepping

            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_elp(pos_x, pos_y, width, height, pixel_intensity, alpha, stepping)
                except Exception as exc:
                    print('  Error while trying to draw the rect, pixel values wont  be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of triplets expected, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def draw_rectangle(dataset, args):
    """ Draws a rectangular region of interest"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        n_vals = 7
        pix_len = len(args)
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = float(args[i + 0]) # x pos
            pos_y = float(args[i + 1]) # y pos
            width = float(args[i + 2]) # rect width
            height = float(args[i + 3]) # rect height
            stepping = int(args[i + 4]) # stepping (for dash lines purpose)
            pixel_intensity = int(args[i + 5]) # pixel intensity
            alpha = float(args[i + 6]) # pixel alpha transparency

            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_rect(pos_x, pos_y, width, height, stepping, pixel_intensity, alpha)
                except Exception as exc:
                    print('  Error while trying to draw the rect, pixel values wont  be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of triplets expected, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def draw_frectangle(dataset, args):
    """ Draws a rectangular region of interest"""
    # assign custom tags
    if args == '':
        return
    try:
        pixel_editor = PixelEditor(dataset.pixel_array)
        n_vals = 6
        pix_len = len(args)
        #if pix_len > 1:
        #    print()
        for i in xrange(0, pix_len / n_vals * n_vals, n_vals):
            pos_x = float(args[i + 0]) # x pos
            pos_y = float(args[i + 1]) # y pos
            width = float(args[i + 2]) # rect width
            height = float(args[i + 3]) # rect height
            pixel_intensity = int(args[i + 4]) # pixel intensity
            alpha = float(args[i + 5]) # pixel alpha transparency

            if pixel_editor.buffer_length() == 0:
                print("  Could not find a pixel array, value won't be set ...")
            else:
                try:
                    pixel_editor.draw_frect(pos_x, pos_y, width, height, pixel_intensity, alpha)
                except Exception as exc:
                    print('  Error while trying to draw the rect, pixel values wont  be set ...')
                    print(exc)
        if pix_len % n_vals != 0:
            print("  Warning: list of triplets expected, but odd count was found instead, " + \
                "found ending: <" + args[pix_len - 1] + '>')

        dataset.PixelData = pixel_editor.buffer_to_string()
    except Exception as exc:
        print(exc)
#------------------------------------------------------------------------------
def change_tag_if_arg(dataset, tag, arg):
    """Change tag if arg is set"""
    if arg != '':
        try:
            dataset.data_element(tag).value = arg
        except Exception as exc:
            print(exc)
#------------------------------------------------------------------------------
def transform_dates(file_count, dataset, args):
    """Transform acquisition date or time, or series date or time"""
    change_tag_if_arg(dataset, "SeriesDate", args.date)
    change_tag_if_arg(dataset, "SeriesTime", args.time)
    change_tag_if_arg(dataset, "StudyDate", args.sdate)
    change_tag_if_arg(dataset, "StudyTime", args.stime)
    change_tag_if_arg(dataset, "ContentDate", args.cdate)
    change_tag_if_arg(dataset, "ContentTime", args.ctime)
    change_tag_if_arg(dataset, "AcquisitionDate", args.adate)
    change_tag_if_arg(dataset, "AcquisitionTime", args.atime)

    try:
        if args.adelta != '':
            secs_offset = str(float(args.adelta) * file_count)
            current_time = modify_time(dataset.AcquisitionDate, dataset.AcquisitionTime, \
                                        secs_offset)
            dataset.AcquisitionDate = get_dicom_date_from(current_time)
            dataset.AcquisitionTime = get_dicom_time_from(current_time)
    except Exception as exc:
        print(exc)

#------------------------------------------------------------------------------
def anonymize_tags_if_anon(dataset, args, remove_curves=False, remove_private_tags=False):
    """ Anonymize dataset tags"""
    # Series Description gets automatically populated with useful transform
    # info by default:
    # Remove patient name and any other person names

    if args.an != '':
        dataset.walk(pn_callback)

    if args.an != '':
        if args.pid != '':
            dataset.PatientID = args.pid
        else:
            dataset.PatientID = 'id'

        change_tag_if_arg(dataset, "InstitutionName", args.an)
        change_tag_if_arg(dataset, "InstitutionAddress", args.an)
        change_tag_if_arg(dataset, "StationName", args.an)
        change_tag_if_arg(dataset, "SequenceName", args.an)
        change_tag_if_arg(dataset, "ProtocolName", args.an)
        change_tag_if_arg(dataset, "ContentDate", '19010101')
        change_tag_if_arg(dataset, "ContentTime", '000000.000000')
        change_tag_if_arg(dataset, "PerformedProcedureStepStartDate", '19010101')
        change_tag_if_arg(dataset, "PerformedProcedureStepStartTime", '000000.000000')
        change_tag_if_arg(dataset, "PerformedProcedureStepID", "0")
        change_tag_if_arg(dataset, "PerformedProcedureStepDescription", args.an)

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

#------------------------------------------------------------------------------
def truncate_str(input_str, maxlen, ending='..'):
    """Truncate teh string if more than maxlen chars """
    if len(input_str) > maxlen:
        return input_str[:(maxlen - len(ending))] + ending
    return input_str
#------------------------------------------------------------------------------
def transform(file_count, args, desc_prefix, input_filename, output_filename):
    """Replace data element values to partly transform a DICOM file.
    Note: completely transforming a DICOM file is very complicated; there
    are many things this example code does not address. USE AT YOUR OWN RISK.
    """

    try:
        file_count += 1
        # Load the current dicom file to 'transform'
        dataset = dicom.read_file(input_filename)

        # 3d xforms user cmd options
        compute_3d_transforms(dataset, args)

        set_image_pixels(dataset, args.pixel)       # set pixels in image buffer
        draw_roi(dataset, args.roi)                 # set a ROI square in image buffer
        draw_ellipse(dataset, args.elp)             # set an ellipse
        draw_rectangle(dataset, args.rect)          # set a rectangle in image buffer
        draw_frectangle(dataset, args.frect)        # set a filled rectangle in image buffer
        draw_crosshair(dataset, args.crosshair)     # set a crosshair in image buffer

        if args.sn > 0:
            dataset.SeriesNumber = args.sn

        # Anonymize dataset tags
        anonymize_tags_if_anon(dataset, args, False, args.delete_private_tags)

        # Deal with all sorts of dates and time if user asks for it:
        transform_dates(file_count, dataset, args)

        #optional changes of useful tags if they have a non default / set value:
        change_tag_if_arg(dataset, "StudyDescription", args.sdesc)
        change_tag_if_arg(dataset, "InstitutionName", args.iname)
        change_tag_if_arg(dataset, "InstitutionAddress", args.iaddr)
        change_tag_if_arg(dataset, "ProtocolName", args.proto)
        change_tag_if_arg(dataset, "Manufacturer", args.mname)
        change_tag_if_arg(dataset, "ManufacturerModelName", args.mmname)
        change_tag_if_arg(dataset, "PatientID", args.pid)
        change_tag_if_arg(dataset, "PatientName", args.pname)
        change_tag_if_arg(dataset, "PatientBirthDate", args.dob)

        # custom DICOM tags settings alternative
        assign_custom_tags(dataset, args.tags)

        # do useful things with the series description
        if args.desc != '':
            change_tag_if_arg(dataset, "SeriesDescription", args.desc) #optionally change study desc
        else: # automatic tracking of transformations
            sdesc = dataset.dir('SeriesDescription');
            
            try:
                if desc_prefix != '':
                    if len(sdesc) != 0:
                        desc = ' ' + dataset.SeriesDescription
                    else:
                        desc = ''
                    dataset.SeriesDescription = desc_prefix + desc
                    sdesc = dataset.dir('SeriesDescription'); #refresh in case we just created it
            except Exception as exc:
                print(exc)
        if len(sdesc) != 0:
            dataset.SeriesDescription = truncate_str(dataset.SeriesDescription, 63)

        # write the 'transformed' DICOM out under the new filename
        dataset.save_as(output_filename)

    except Exception as exc:
        print(exc)

    return file_count, dataset
#------------------------------------------------------------------------------
def is_3d_tranformation(in_args):
    """Determine if any 3d transform on image position patient needs to be computed"""
    if in_args.x != 0.0 or in_args.y != 0.0 or in_args.z != 0.0 or \
        in_args.ax != 0.0 or in_args.ay != 0.0 or in_args.az != 0.0:
        return True
    return False
#------------------------------------------------------------------------------
def iterate_once(in_args, input_dir, output_dir):
    """Execute the full script except the recursive option"""
    series_file_count = 0

    try:
        if is_3d_tranformation(in_args):
            series_desc_prefix = 'T[' + fmt_float3d('', in_args.x, in_args.y, in_args.z, ' ') \
                + fmt_float3d('A', in_args.ax, in_args.ay, in_args.az, ' ') + ']'
        else:
            series_desc_prefix = ''

    except Exception:
        print("Could not convert the x, y, z offsets")
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
                series_file_count, dataset = transform(series_file_count, ARGS, series_desc_prefix, \
                                       os.path.join(input_dir, filename), \
                                       os.path.join(output_dir, filename))
                print(" done\r")
    else:  # first arg not a directory, assume two files given
        in_filename = in_args.input_series
        out_filename = in_args.output_series
        transform(series_file_count, ARGS, series_desc_prefix, in_filename, out_filename)
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
