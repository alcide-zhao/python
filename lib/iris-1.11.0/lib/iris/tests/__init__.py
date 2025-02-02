# (C) British Crown Copyright 2010 - 2016, Met Office
#
# This file is part of Iris.
#
# Iris is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Iris is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Iris.  If not, see <http://www.gnu.org/licenses/>.
"""
Provides testing capabilities and customisations specific to Iris.

.. note:: This module needs to control the matplotlib backend, so it
          **must** be imported before ``matplotlib.pyplot``.

The primary class for this module is :class:`IrisTest`.

By default, this module sets the matplotlib backend to "agg". But when
this module is imported it checks ``sys.argv`` for the flag "-d". If
found, it is removed from ``sys.argv`` and the matplotlib backend is
switched to "tkagg" to allow the interactive visual inspection of
graphical test results.

"""

from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa
import six

import codecs
import collections
import contextlib
import difflib
import filecmp
import functools
import gzip
import inspect
import json
import io
import logging
import math
import os
import os.path
import shutil
import subprocess
import sys
import unittest
import threading
import warnings
import xml.dom.minidom
import zlib

try:
    from unittest import mock
except ImportError:
    import mock

import filelock
import numpy as np
import numpy.ma as ma
import requests

import iris.cube
import iris.config
import iris.util

# Test for availability of matplotlib.
# (And remove matplotlib as an iris.tests dependency.)
try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.testing.compare as mcompare
    import matplotlib.pyplot as plt
except ImportError:
    MPL_AVAILABLE = False
else:
    MPL_AVAILABLE = True

try:
    from osgeo import gdal
except ImportError:
    GDAL_AVAILABLE = False
else:
    GDAL_AVAILABLE = True

try:
    import iris_grib
    GRIB_AVAILABLE = True
    from iris_grib.message import GribMessage
except ImportError:
    try:
        import gribapi
        GRIB_AVAILABLE = True
        from iris.fileformats.grib.message import GribMessage
    except ImportError:
        GRIB_AVAILABLE = False

try:
    import iris_sample_data
except ImportError:
    SAMPLE_DATA_AVAILABLE = False
else:
    SAMPLE_DATA_AVAILABLE = True

try:
    import nc_time_axis
    NC_TIME_AXIS_AVAILABLE = True
except ImportError:
    NC_TIME_AXIS_AVAILABLE = False

try:
    requests.get('https://github.com/SciTools/iris')
    INET_AVAILABLE = True
except requests.exceptions.ConnectionError:
    INET_AVAILABLE = False

#: Basepath for test results.
_RESULT_PATH = os.path.join(os.path.dirname(__file__), 'results')
#: Default perceptual hash size.
_HASH_SIZE = 16
#: Default maximum perceptual hash hamming distance.
_HAMMING_DISTANCE = 2

if '--data-files-used' in sys.argv:
    sys.argv.remove('--data-files-used')
    fname = '/var/tmp/all_iris_test_resource_paths.txt'
    print('saving list of files used by tests to %s' % fname)
    _EXPORT_DATAPATHS_FILE = open(fname, 'w')
else:
    _EXPORT_DATAPATHS_FILE = None


if '--create-missing' in sys.argv:
    sys.argv.remove('--create-missing')
    print('Allowing creation of missing test results.')
    os.environ['IRIS_TEST_CREATE_MISSING'] = 'true'


# A shared logger for use by unit tests
logger = logging.getLogger('tests')

# Whether to display matplotlib output to the screen.
_DISPLAY_FIGURES = False

if (MPL_AVAILABLE and '-d' in sys.argv):
    sys.argv.remove('-d')
    plt.switch_backend('tkagg')
    _DISPLAY_FIGURES = True

# Threading non re-entrant blocking lock to ensure thread-safe plotting.
_lock = threading.Lock()


def main():
    """A wrapper for unittest.main() which adds iris.test specific options to the help (-h) output."""
    if '-h' in sys.argv or '--help' in sys.argv:
        stdout = sys.stdout
        buff = io.StringIO()
        # NB. unittest.main() raises an exception after it's shown the help text
        try:
            sys.stdout = buff
            unittest.main()
        finally:
            sys.stdout = stdout
            lines = buff.getvalue().split('\n')
            lines.insert(9, 'Iris-specific options:')
            lines.insert(10, '  -d                   Display matplotlib figures (uses tkagg).')
            lines.insert(11, '                       NOTE: To compare results of failing tests, ')
            lines.insert(12, '                             use idiff.py instead')
            lines.insert(13, '  --data-files-used    Save a list of files used to a temporary file')
            lines.insert(
                14, '  -m                   Create missing test results')
            print('\n'.join(lines))
    else:
        unittest.main()


def get_data_path(relative_path):
    """
    Return the absolute path to a data file when given the relative path
    as a string, or sequence of strings.

    """
    if not isinstance(relative_path, six.string_types):
        relative_path = os.path.join(*relative_path)
    test_data_dir = iris.config.TEST_DATA_DIR
    if test_data_dir is None:
        test_data_dir = ''
    data_path = os.path.join(test_data_dir, relative_path)

    if _EXPORT_DATAPATHS_FILE is not None:
        _EXPORT_DATAPATHS_FILE.write(data_path + '\n')

    if isinstance(data_path, six.string_types) and not os.path.exists(data_path):
        # if the file is gzipped, ungzip it and return the path of the ungzipped
        # file.
        gzipped_fname = data_path + '.gz'
        if os.path.exists(gzipped_fname):
            with gzip.open(gzipped_fname, 'rb') as gz_fh:
                try:
                    with open(data_path, 'wb') as fh:
                        fh.writelines(gz_fh)
                except IOError:
                    # Put ungzipped data file in a temporary path, since we
                    # can't write to the original path (maybe it is owned by
                    # the system.)
                    _, ext = os.path.splitext(data_path)
                    data_path = iris.util.create_temp_filename(suffix=ext)
                    with open(data_path, 'wb') as fh:
                        fh.writelines(gz_fh)

    return data_path


class IrisTest(unittest.TestCase):
    """A subclass of unittest.TestCase which provides Iris specific testing functionality."""

    _assertion_counts = collections.defaultdict(int)

    @classmethod
    def setUpClass(cls):
        # Ensure that the CF profile if turned-off for testing.
        iris.site_configuration['cf_profile'] = None

    def _assert_str_same(self, reference_str, test_str, reference_filename, type_comparison_name='Strings'):
        if reference_str != test_str:
            diff = ''.join(difflib.unified_diff(reference_str.splitlines(1), test_str.splitlines(1),
                                                 'Reference', 'Test result', '', '', 0))
            self.fail("%s do not match: %s\n%s" % (type_comparison_name, reference_filename, diff))

    @staticmethod
    def get_result_path(relative_path):
        """
        Returns the absolute path to a result file when given the relative path
        as a string, or sequence of strings.

        """
        if not isinstance(relative_path, six.string_types):
            relative_path = os.path.join(*relative_path)
        return os.path.abspath(os.path.join(_RESULT_PATH, relative_path))

    def result_path(self, basename=None, ext=''):
        """
        Return the full path to a test result, generated from the \
        calling file, class and, optionally, method.

        Optional kwargs :

            * basename    - File basename. If omitted, this is \
                            generated from the calling method.
            * ext         - Appended file extension.

        """
        if ext and not ext.startswith('.'):
            ext = '.' + ext

        # Generate the folder name from the calling file name.
        path = os.path.abspath(inspect.getfile(self.__class__))
        path = os.path.splitext(path)[0]
        sub_path = path.rsplit('iris', 1)[1].split('tests', 1)[1][1:]

        # Generate the file name from the calling function name?
        if basename is None:
            stack = inspect.stack()
            for frame in stack[1:]:
                if 'test_' in frame[3]:
                    basename = frame[3].replace('test_', '')
                    break
        filename = basename + ext

        result = os.path.join(self.get_result_path(''),
                              sub_path.replace('test_', ''),
                              self.__class__.__name__.replace('Test_', ''),
                              filename)
        return result

    def assertCMLApproxData(self, cubes, reference_filename=None, **kwargs):
        # passes args and kwargs on to approx equal
        if isinstance(cubes, iris.cube.Cube):
            cubes = [cubes]
        if reference_filename is None:
            reference_filename = self.result_path(None, 'cml')
            reference_filename = [self.get_result_path(reference_filename)]
        for i, cube in enumerate(cubes):
            fname = list(reference_filename)
            # don't want the ".cml" for the json stats file
            if fname[-1].endswith(".cml"):
                fname[-1] = fname[-1][:-4]
            fname[-1] += '.data.%d.json' % i
            self.assertCubeDataAlmostEqual(cube, fname, **kwargs)
        self.assertCML(cubes, reference_filename, checksum=False)

    def assertCDL(self, netcdf_filename, reference_filename=None, flags='-h'):
        """
        Test that the CDL for the given netCDF file matches the contents
        of the reference file.

        If the environment variable IRIS_TEST_CREATE_MISSING is
        non-empty, the reference file is created if it doesn't exist.

        Args:

        * netcdf_filename:
            The path to the netCDF file.

        Kwargs:

        * reference_filename:
            The relative path (relative to the test results directory).
            If omitted, the result is generated from the calling
            method's name, class, and module using
            :meth:`iris.tests.IrisTest.result_path`.

        * flags:
            Command-line flags for `ncdump`, as either a whitespace
            separated string or an iterable. Defaults to '-h'.

        """
        if reference_filename is None:
            reference_path = self.result_path(None, 'cdl')
        else:
            reference_path = self.get_result_path(reference_filename)

        # Convert the netCDF file to CDL file format.
        cdl_filename = iris.util.create_temp_filename(suffix='.cdl')

        if flags is None:
            flags = []
        elif isinstance(flags, six.string_types):
            flags = flags.split()
        else:
            flags = list(map(str, flags))

        with open(cdl_filename, 'w') as cdl_file:
            subprocess.check_call(['ncdump'] + flags + [netcdf_filename],
                                  stderr=cdl_file, stdout=cdl_file)

        # Ingest the CDL for comparison, excluding first line.
        with open(cdl_filename, 'r') as cdl_file:
            lines = cdl_file.readlines()[1:]

        # Sort the dimensions (except for the first, which can be unlimited).
        # This gives consistent CDL across different platforms.
        sort_key = lambda line: ('UNLIMITED' not in line, line)
        dimension_lines = slice(lines.index('dimensions:\n') + 1,
                                lines.index('variables:\n'))
        lines[dimension_lines] = sorted(lines[dimension_lines], key=sort_key)
        cdl = ''.join(lines)

        os.remove(cdl_filename)
        self._check_same(cdl, reference_path, type_comparison_name='CDL')

    def assertCML(self, cubes, reference_filename=None, checksum=True):
        """
        Test that the CML for the given cubes matches the contents of
        the reference file.

        If the environment variable IRIS_TEST_CREATE_MISSING is
        non-empty, the reference file is created if it doesn't exist.

        Args:

        * cubes:
            Either a Cube or a sequence of Cubes.

        Kwargs:

        * reference_filename:
            The relative path (relative to the test results directory).
            If omitted, the result is generated from the calling
            method's name, class, and module using
            :meth:`iris.tests.IrisTest.result_path`.

        * checksum:
            When True, causes the CML to include a checksum for each
            Cube's data. Defaults to True.

        """
        if isinstance(cubes, iris.cube.Cube):
            cubes = [cubes]
        if reference_filename is None:
            reference_filename = self.result_path(None, 'cml')

        if isinstance(cubes, (list, tuple)):
            xml = iris.cube.CubeList(cubes).xml(checksum=checksum, order=False,
                                                byteorder=False)
        else:
            xml = cubes.xml(checksum=checksum, order=False, byteorder=False)
        reference_path = self.get_result_path(reference_filename)
        self._check_same(xml, reference_path)

    def assertTextFile(self, source_filename, reference_filename, desc="text file"):
        """Check if two text files are the same, printing any diffs."""
        with open(source_filename) as source_file:
            source_text = source_file.readlines()
        with open(reference_filename) as reference_file:
            reference_text = reference_file.readlines()
        if reference_text != source_text:
            diff = ''.join(difflib.unified_diff(reference_text, source_text, 'Reference', 'Test result', '', '', 0))
            self.fail("%s does not match reference file: %s\n%s" % (desc, reference_filename, diff))

    def assertCubeDataAlmostEqual(self, cube, reference_filename, **kwargs):
        reference_path = self.get_result_path(reference_filename)
        if self._check_reference_file(reference_path):
            kwargs.setdefault('err_msg', 'Reference file %s' % reference_path)
            with open(reference_path, 'r') as reference_file:
                stats = json.load(reference_file)
                self.assertEqual(stats.get('shape', []), list(cube.shape))
                self.assertEqual(stats.get('masked', False),
                                       isinstance(cube.data, ma.MaskedArray))
                nstats = np.array((stats.get('mean', 0.), stats.get('std', 0.),
                                   stats.get('max', 0.), stats.get('min', 0.)),
                                  dtype=np.float_)
                if math.isnan(stats.get('mean', 0.)):
                    self.assertTrue(math.isnan(cube.data.mean()))
                else:
                    cube_stats = np.array((cube.data.mean(), cube.data.std(),
                                           cube.data.max(), cube.data.min()),
                                          dtype=np.float_)
                    self.assertArrayAllClose(nstats, cube_stats, **kwargs)
        else:
            self._ensure_folder(reference_path)
            logger.warning('Creating result file: %s', reference_path)
            masked = False
            if isinstance(cube.data, ma.MaskedArray):
                masked = True
            stats = {'mean': np.float_(cube.data.mean()),
                     'std': np.float_(cube.data.std()),
                     'max': np.float_(cube.data.max()),
                     'min': np.float_(cube.data.min()),
                     'shape': cube.shape, 'masked': masked}
            with open(reference_path, 'w') as reference_file:
                reference_file.write(json.dumps(stats))

    def assertFilesEqual(self, test_filename, reference_filename):
        reference_path = self.get_result_path(reference_filename)
        if self._check_reference_file(reference_path):
            fmt = 'test file {!r} does not match reference {!r}.'
            self.assertTrue(filecmp.cmp(test_filename, reference_path),
                            fmt.format(test_filename, reference_path))
        else:
            self._ensure_folder(reference_path)
            logger.warning('Creating result file: %s', reference_path)
            shutil.copy(test_filename, reference_path)

    def assertString(self, string, reference_filename=None):
        """
        Test that `string` matches the contents of the reference file.

        If the environment variable IRIS_TEST_CREATE_MISSING is
        non-empty, the reference file is created if it doesn't exist.

        Args:

        * string:
            The string to check.

        Kwargs:

        * reference_filename:
            The relative path (relative to the test results directory).
            If omitted, the result is generated from the calling
            method's name, class, and module using
            :meth:`iris.tests.IrisTest.result_path`.

        """
        if reference_filename is None:
            reference_path = self.result_path(None, 'txt')
        else:
            reference_path = self.get_result_path(reference_filename)
        self._check_same(string, reference_path,
                         type_comparison_name='Strings')

    def assertRepr(self, obj, reference_filename):
        self.assertString(repr(obj), reference_filename)

    def _check_same(self, item, reference_path, type_comparison_name='CML'):
        if self._check_reference_file(reference_path):
            with open(reference_path, 'rb') as reference_fh:
                reference = ''.join(part.decode('utf-8')
                                    for part in reference_fh.readlines())
            self._assert_str_same(reference, item, reference_path,
                                  type_comparison_name)
        else:
            self._ensure_folder(reference_path)
            logger.warning('Creating result file: %s', reference_path)
            with open(reference_path, 'wb') as reference_fh:
                reference_fh.writelines(
                    part.encode('utf-8')
                    for part in item)

    def assertXMLElement(self, obj, reference_filename):
        """
        Calls the xml_element method given obj and asserts the result is the same as the test file.

        """
        doc = xml.dom.minidom.Document()
        doc.appendChild(obj.xml_element(doc))
        pretty_xml = doc.toprettyxml(indent="  ")
        reference_path = self.get_result_path(reference_filename)
        self._check_same(pretty_xml, reference_path,
                         type_comparison_name='XML')

    def assertArrayEqual(self, a, b, err_msg=''):
        np.testing.assert_array_equal(a, b, err_msg=err_msg)

    def _assertMaskedArray(self, assertion, a, b, strict, **kwargs):
        # Define helper function to extract unmasked values as a 1d
        # array.
        def unmasked_data_as_1d_array(array):
            if array.ndim == 0:
                if array.mask:
                    data = np.array([])
                else:
                    data = np.array([array.data])
            else:
                data = array.data[~ma.getmaskarray(array)]
            return data

        # Compare masks. This will also check that the array shapes
        # match, which is not tested when comparing unmasked values if
        # strict is False.
        a_mask, b_mask = ma.getmaskarray(a), ma.getmaskarray(b)
        np.testing.assert_array_equal(a_mask, b_mask)

        if strict:
            assertion(a.data, b.data, **kwargs)
        else:
            assertion(unmasked_data_as_1d_array(a),
                      unmasked_data_as_1d_array(b),
                      **kwargs)

    def assertMaskedArrayEqual(self, a, b, strict=False):
        """
        Check that masked arrays are equal. This requires the
        unmasked values and masks to be identical.

        Args:

        * a, b (array-like):
            Two arrays to compare.

        Kwargs:

        * strict (bool):
            If True, perform a complete mask and data array equality check.
            If False (default), the data array equality considers only unmasked
            elements.

        """
        self._assertMaskedArray(np.testing.assert_array_equal, a, b, strict)

    def assertArrayAlmostEqual(self, a, b, decimal=6):
        np.testing.assert_array_almost_equal(a, b, decimal=decimal)

    def assertMaskedArrayAlmostEqual(self, a, b, decimal=6, strict=False):
        """
        Check that masked arrays are almost equal. This requires the
        masks to be identical, and the unmasked values to be almost
        equal.

        Args:

        * a, b (array-like):
            Two arrays to compare.

        Kwargs:

        * strict (bool):
            If True, perform a complete mask and data array equality check.
            If False (default), the data array equality considers only unmasked
            elements.

        * decimal (int):
            Equality tolerance level for
            :meth:`numpy.testing.assert_array_almost_equal`, with the meaning
            'abs(desired-actual) < 0.5 * 10**(-decimal)'

        """
        self._assertMaskedArray(np.testing.assert_array_almost_equal, a, b,
                                strict, decimal=decimal)

    def assertArrayAllClose(self, a, b, rtol=1.0e-7, atol=0.0, **kwargs):
        """
        Check arrays are equal, within given relative + absolute tolerances.

        Args:

        * a, b (array-like):
            Two arrays to compare.

        Kwargs:

        * rtol, atol (float):
            Relative and absolute tolerances to apply.

        Any additional kwargs are passed to numpy.testing.assert_allclose.

        Performs pointwise toleranced comparison, and raises an assertion if
        the two are not equal 'near enough'.
        For full details see underlying routine numpy.testing.assert_allclose.

        """
        np.testing.assert_allclose(a, b, rtol=rtol, atol=atol, **kwargs)

    @contextlib.contextmanager
    def temp_filename(self, suffix=''):
        filename = iris.util.create_temp_filename(suffix)
        try:
            yield filename
        finally:
            os.remove(filename)

    def file_checksum(self, file_path):
        """
        Generate checksum from file.
        """
        with open(file_path, "rb") as in_file:
            return zlib.crc32(in_file.read())

    def _unique_id(self):
        """
        Returns the unique ID for the current assertion.

        The ID is composed of two parts: a unique ID for the current test
        (which is itself composed of the module, class, and test names), and
        a sequential counter (specific to the current test) that is incremented
        on each call.

        For example, calls from a "test_tx" routine followed by a "test_ty"
        routine might result in::
            test_plot.TestContourf.test_tx.0
            test_plot.TestContourf.test_tx.1
            test_plot.TestContourf.test_tx.2
            test_plot.TestContourf.test_ty.0

        """
        # Obtain a consistent ID for the current test.
        # NB. unittest.TestCase.id() returns different values depending on
        # whether the test has been run explicitly, or via test discovery.
        # For example:
        #   python tests/test_plot.py => '__main__.TestContourf.test_tx'
        #   ird -t => 'iris.tests.test_plot.TestContourf.test_tx'
        bits = self.id().split('.')
        if bits[0] == '__main__':
            floc = sys.modules['__main__'].__file__
            path, file_name = os.path.split(os.path.abspath(floc))
            bits[0] = os.path.splitext(file_name)[0]
            folder, location = os.path.split(path)
            bits = [location] + bits
            while location not in ['iris', 'example_tests']:
                folder, location = os.path.split(folder)
                bits = [location] + bits
        test_id = '.'.join(bits)

        # Derive the sequential assertion ID within the test
        assertion_id = self._assertion_counts[test_id]
        self._assertion_counts[test_id] += 1

        return test_id + '.' + str(assertion_id)

    def _check_reference_file(self, reference_path):
        reference_exists = os.path.isfile(reference_path)
        if not (reference_exists or
                os.environ.get('IRIS_TEST_CREATE_MISSING')):
            msg = 'Missing test result: {}'.format(reference_path)
            raise AssertionError(msg)
        return reference_exists

    def _ensure_folder(self, path):
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            logger.warning('Creating folder: %s', dir_path)
            os.makedirs(dir_path)

    def check_graphic(self):
        """
        Check the hash of the current matplotlib figure matches the expected
        image hash for the current graphic test.

        To create missing image test results, set the IRIS_TEST_CREATE_MISSING
        environment variable before running the tests. This will result in new
        and appropriately "<hash>.png" image files being generated in the image
        output directory, and the imagerepo.json file being updated.

        """
        import imagehash
        from PIL import Image

        dev_mode = os.environ.get('IRIS_TEST_CREATE_MISSING')
        unique_id = self._unique_id()
        repo_fname = os.path.join(_RESULT_PATH, 'imagerepo.json')
        with open(repo_fname, 'rb') as fi:
            repo = json.load(codecs.getreader('utf-8')(fi))

        try:
            #: The path where the images generated by the tests should go.
            image_output_directory = os.path.join(os.path.dirname(__file__),
                                                  'result_image_comparison')
            if not os.access(image_output_directory, os.W_OK):
                if not os.access(os.getcwd(), os.W_OK):
                    raise IOError('Write access to a local disk is required '
                                  'to run image tests.  Run the tests from a '
                                  'current working directory you have write '
                                  'access to to avoid this issue.')
                else:
                    image_output_directory = os.path.join(
                        os.getcwd(), 'iris_image_test_output')
            result_fname = os.path.join(image_output_directory,
                                        'result-' + unique_id + '.png')

            if not os.path.isdir(image_output_directory):
                # Handle race-condition where the directories are
                # created sometime between the check above and the
                # creation attempt below.
                try:
                    os.makedirs(image_output_directory)
                except OSError as err:
                    # Don't care about "File exists"
                    if err.errno != 17:
                        raise

            def _create_missing():
                fname = '{}.png'.format(phash)
                base_uri = ('https://scitools.github.io/test-iris-imagehash/'
                            'images/{}')
                uri = base_uri.format(fname)
                hash_fname = os.path.join(image_output_directory, fname)
                uris = repo.setdefault(unique_id, [])
                uris.append(uri)
                print('Creating image file: {}'.format(hash_fname))
                figure.savefig(hash_fname)
                msg = 'Creating imagerepo entry: {} -> {}'
                print(msg.format(unique_id, uri))
                lock = filelock.FileLock(os.path.join(_RESULT_PATH,
                                                      'imagerepo.lock'))
                # The imagerepo.json file is a critical resource, so ensure
                # thread safe read/write behaviour via platform independent
                # file locking.
                with lock.acquire(timeout=600):
                    with open(repo_fname, 'wb') as fo:
                        json.dump(repo, codecs.getwriter('utf-8')(fo),
                                  indent=4, sort_keys=True)

            # Calculate the test result perceptual image hash.
            buffer = io.BytesIO()
            figure = plt.gcf()
            figure.savefig(buffer, format='png')
            buffer.seek(0)
            phash = imagehash.phash(Image.open(buffer), hash_size=_HASH_SIZE)

            if unique_id not in repo:
                if dev_mode:
                    _create_missing()
                else:
                    figure.savefig(result_fname)
                    emsg = 'Missing image test result: {}.'
                    raise ValueError(emsg.format(unique_id))
            else:
                uris = repo[unique_id]
                # Create the expected perceptual image hashes from the uris.
                to_hash = imagehash.hex_to_hash
                expected = [to_hash(os.path.splitext(os.path.basename(uri))[0],
                                    hash_size=_HASH_SIZE)
                            for uri in uris]

                # Calculate the hamming distance vector for the result hash.
                distances = [e - phash for e in expected]

                if np.all([hd > _HAMMING_DISTANCE for hd in distances]):
                    if dev_mode:
                        _create_missing()
                    else:
                        figure.savefig(result_fname)
                        msg = ('Bad phash {} with hamming distance {} '
                               'for test {}.')
                        msg = msg.format(phash, distances, unique_id)
                        if _DISPLAY_FIGURES:
                            emsg = 'Image comparion would have failed: {}'
                            print(emsg.format(msg))
                        else:
                            emsg = 'Image comparison failed: {}'
                            raise ValueError(emsg.format(msg))

            if _DISPLAY_FIGURES:
                plt.show()

        finally:
            plt.close()

    def _remove_testcase_patches(self):
        """Helper to remove per-testcase patches installed by :meth:`patch`."""
        # Remove all patches made, ignoring errors.
        for p in self.testcase_patches:
            p.stop()
        # Reset per-test patch control variable.
        self.testcase_patches.clear()

    def patch(self, *args, **kwargs):
        """
        Install a mock.patch, to be removed after the current test.

        The patch is created with mock.patch(*args, **kwargs).

        Returns:
            The substitute object returned by patch.start().

        For example::

            mock_call = self.patch('module.Class.call', return_value=1)
            module_Class_instance.call(3, 4)
            self.assertEqual(mock_call.call_args_list, [mock.call(3, 4)])

        """
        # Make the new patch and start it.
        patch = mock.patch(*args, **kwargs)
        start_result = patch.start()

        # Create the per-testcases control variable if it does not exist.
        # NOTE: this mimics a setUp method, but continues to work when a
        # subclass defines its own setUp.
        if not hasattr(self, 'testcase_patches'):
            self.testcase_patches = {}

        # When installing the first patch, schedule remove-all at cleanup.
        if not self.testcase_patches:
            self.addCleanup(self._remove_testcase_patches)

        # Record the new patch and start object for reference.
        self.testcase_patches[patch] = start_result

        # Return patch replacement object.
        return start_result

    def assertArrayShapeStats(self, result, shape, mean, std_dev):
        """
        Assert that the result, a cube, has the provided shape and that the
        mean and standard deviation of the data array are also as provided.
        Thus build confidence that a cube processing operation, such as a
        cube.regrid, has maintained its behaviour.

        """
        self.assertEqual(result.shape, shape)
        self.assertAlmostEqual(result.data.mean(), mean, places=5)
        self.assertAlmostEqual(result.data.std(), std_dev, places=5)


get_result_path = IrisTest.get_result_path


class GraphicsTest(IrisTest):

    # nose directive: dispatch tests concurrently.
    _multiprocess_can_split_ = True

    def setUp(self):
        # Acquire threading non re-entrant blocking lock to ensure
        # thread-safe plotting.
        _lock.acquire()
        # Make sure we have no unclosed plots from previous tests before
        # generating this one.
        if MPL_AVAILABLE:
            plt.close('all')

    def tearDown(self):
        # If a plotting test bombs out it can leave the current figure
        # in an odd state, so we make sure it's been disposed of.
        if MPL_AVAILABLE:
            plt.close('all')
        # Release the non re-entrant blocking lock.
        _lock.release()


class TestGribMessage(IrisTest):
    def assertGribMessageContents(self, filename, contents):
        """
        Evaluate whether all messages in a GRIB2 file contain the provided
        contents.

        * filename (string)
            The path on disk of an existing GRIB file

        * contents
            An iterable of GRIB message keys and expected values.

        """
        messages = GribMessage.messages_from_filename(filename)
        for message in messages:
            for element in contents:
                section, key, val = element
                self.assertEqual(message.sections[section][key], val)

    def assertGribMessageDifference(self, filename1, filename2, diffs,
                                    skip_keys=(), skip_sections=()):
        """
        Evaluate that the two messages only differ in the ways specified.

        * filename[0|1] (string)
            The path on disk of existing GRIB files

        * diffs
            An dictionary of GRIB message keys and expected diff values:
            {key: (m1val, m2val),...} .

        * skip_keys
            An iterable of key names to ignore during comparison.

        * skip_sections
            An iterable of section numbers to ignore during comparison.

        """
        messages1 = list(GribMessage.messages_from_filename(filename1))
        messages2 = list(GribMessage.messages_from_filename(filename2))
        self.assertEqual(len(messages1), len(messages2))
        for m1, m2 in zip(messages1, messages2):
            m1_sect = set(m1.sections.keys())
            m2_sect = set(m2.sections.keys())

            for missing_section in (m1_sect ^ m2_sect):
                what = ('introduced'
                        if missing_section in m1_sect else 'removed')
                # Assert that an introduced section is in the diffs.
                self.assertIn(missing_section, skip_sections,
                              msg='Section {} {}'.format(missing_section,
                                                         what))

            for section in (m1_sect & m2_sect):
                # For each section, check that the differences are
                # known diffs.
                m1_keys = set(m1.sections[section]._keys)
                m2_keys = set(m2.sections[section]._keys)

                difference = m1_keys ^ m2_keys
                unexpected_differences = difference - set(skip_keys)
                if unexpected_differences:
                    self.fail("There were keys in section {} which \n"
                              "weren't in both messages and which weren't "
                              "skipped.\n{}"
                              "".format(section,
                                        ', '.join(unexpected_differences)))

                keys_to_compare = m1_keys & m2_keys - set(skip_keys)

                for key in keys_to_compare:
                    m1_value = m1.sections[section][key]
                    m2_value = m2.sections[section][key]
                    msg = '{} {} != {}'
                    if key not in diffs:
                        # We have a key which we expect to be the same for
                        # both messages.
                        if isinstance(m1_value, np.ndarray):
                            # A large tolerance appears to be required for
                            # gribapi 1.12, but not for 1.14.
                            self.assertArrayAlmostEqual(m1_value, m2_value,
                                                        decimal=2)
                        else:
                            self.assertEqual(m1_value, m2_value,
                                             msg=msg.format(key, m1_value,
                                                            m2_value))
                    else:
                        # We have a key which we expect to be different
                        # for each message.
                        self.assertEqual(m1_value, diffs[key][0],
                                         msg=msg.format(key, m1_value,
                                                         diffs[key][0]))

                        self.assertEqual(m2_value, diffs[key][1],
                                         msg=msg.format(key, m2_value,
                                                        diffs[key][1]))


def skip_data(fn):
    """
    Decorator to choose whether to run tests, based on the availability of
    external data.

    Example usage:
        @skip_data
        class MyDataTests(tests.IrisTest):
            ...

    """
    no_data = (not iris.config.TEST_DATA_DIR
               or not os.path.isdir(iris.config.TEST_DATA_DIR)
               or os.environ.get('IRIS_TEST_NO_DATA'))

    skip = unittest.skipIf(
        condition=no_data,
        reason='Test(s) require external data.')

    return skip(fn)


def skip_gdal(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    GDAL library.

    Example usage:
        @skip_gdal
        class MyGeoTiffTests(test.IrisTest):
            ...

    """
    skip = unittest.skipIf(
        condition=not GDAL_AVAILABLE,
        reason="Test requires 'gdal'.")
    return skip(fn)


def skip_plot(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    matplotlib library.

    Example usage:
        @skip_plot
        class MyPlotTests(test.GraphicsTest):
            ...

    """
    skip = unittest.skipIf(
        condition=not MPL_AVAILABLE,
        reason='Graphics tests require the matplotlib library.')

    return skip(fn)


skip_grib = unittest.skipIf(not GRIB_AVAILABLE, 'Test(s) require "gribapi", '
                                                'which is not available.')


skip_sample_data = unittest.skipIf(not SAMPLE_DATA_AVAILABLE,
                                   ('Test(s) require "iris_sample_data", '
                                    'which is not available.'))


skip_nc_time_axis = unittest.skipIf(
    not NC_TIME_AXIS_AVAILABLE,
    'Test(s) require "nc_time_axis", which is not available.')


skip_inet = unittest.skipIf(not INET_AVAILABLE,
                            ('Test(s) require an "internet connection", '
                             'which is not available.'))


def no_warnings(func):
    """
    Provides a decorator to ensure that there are no warnings raised
    within the test, otherwise the test will fail.

    """
    @functools.wraps(func)
    def wrapped(self, *args, **kwargs):
        with mock.patch('warnings.warn') as warn:
            result = func(self, *args, **kwargs)
        self.assertEqual(0, warn.call_count,
                         ('Got unexpected warnings.'
                          ' \n{}'.format(warn.call_args_list)))
        return result
    return wrapped
