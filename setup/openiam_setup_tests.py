# -*- coding: utf-8 -*-
"""
OpenIAM setup script to test Python
environment, compile components, and
test OpenIAM functionalility.
Created on Tue Dec 12 10:02:24 2017
Modified Feb 05 2020

@author: Seth King
AECOM supporting NETL
Seth.King@netl.doe.gov
"""
import os
import sys
from distutils.spawn import find_executable
import pkgutil
import unittest
import importlib
import logging
import pkg_resources

from make_libs import make_from_file


def setup_logging():
    """
    Set up logging for the console and for a file.
    """
    logging.basicConfig(level=logging.INFO,
                        format='from %(module)s %(funcName)s - %(levelname)s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='OpenIAM_setup.log',
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)


def is_tool(name):
    """Check whether `name` is on PATH."""
    return find_executable(name) is not None


def check_python():
    """
    Check whether the proper version of Python is being used.

    Currently python above version 3.6 is supported.
    """
    logging.info('Checking Python version')
    py_major, py_minor = sys.version_info[0:2]
    if py_major != 3 or py_minor < 6:
        error_msg = ''.join([
            'Python 3.6 or above is required.\nPython {maj}.{mnr} ',
            'has been detected']).format(maj=py_major, mnr=py_minor)
        logging.error(error_msg)
    else:
        debug_msg = 'Python version: {info}'.format(info=sys.version_info)
        logging.debug(debug_msg)


def get_package_name(package):
    """Get the package name from the distribution name."""
    return list(pkg_resources.get_distribution(package)._get_metadata('top_level.txt'))[0]


def test_dependencies():
    """
    Test whether all required Python libraries are installed.

    Test whether all needed libraries are installed and satisfy the
    minimum version requirements. Uses python_libraries.txt file.
    """
    logging.info('Checking NRAP-Open-IAM dependencies')
    pkgs = {}
    with open('python_libraries.txt', 'r') as req_file:
        for line in req_file:
            pkg, vers = line.strip().split('==')
            pkgs[pkg] = vers
    # Check PyYaml, because it does not link to package name correctly
    try:
        import yaml
    except ImportError:
        logging.error(''.join([
            'Pyyaml is not detected, install library using the command ',
            '"conda install -c conda-forge pyyaml" or "pip install pyyaml" and run setup again.']))
        sys.exit()
    # Check Pmw
    try:
        import Pmw
    except ImportError:
        logging.error(''.join([
            'Pmw is not detected, install library using the command ',
            '"conda install -c conda-forge pmw" or "pip install pmw" and run setup again.']))
        sys.exit()

    for package, version in pkgs.items():
        found = pkgutil.find_loader(package)
        if not found:
            try:
                pkg_name = get_package_name(pkg)
                found = pkgutil.find_loader(pkg_name)
            except ImportError:
                error_msg = ''.join([
                    'Package {pack} is not found. Try installing the package ',
                    'using "conda install {pack}" or "pip install {pack}" then rerun the setup']).format(
                        pack=package)
                logging.error(error_msg)
            except IndexError:
                warn_msg = ''.join([
                    'Index Error with package {}. Check whether package ',
                    'is installed']).format(package)
                logging.warning(warn_msg)
        else:
            ldpkg = importlib.import_module(package)
            info_msg = 'Testing package {} ...'.format(package)
            logging.info(info_msg)
            if hasattr(ldpkg, '__version__'):
                if ldpkg.__version__ < version:
                    warn_msg = ''.join([
                        'Version {} of package {} recommended, ',
                        'version {} installed']).format(
                            version, package, ldpkg.__version__)
                    logging.warning(warn_msg)


def check_compiler():
    """
    Check that gfortran and make are installed (required for macOS and Linux).
    """
    if is_tool('gfortran'):
        state = True
#        if is_tool('make'):
#            state = True
#        else:
#            state = False
#            raise AssertionError('The make utility program must be installed ' +
#                                 'for the OpenIAM to compile the needed libraries.  ' +
#                                 'Please install make and run the setup again.')
    else:
        state = False
        raise AssertionError(''.join([
            'The gfortran compiler must be installed for the NRAP-Open-IAM ',
            'to compile the needed libraries. Please install gfortran ',
            'and run the setup again.']))
    return state


def create_libraries():
    """
    Check the existence of internal NRAP-Open-IAM libraries.

    For Windows systems test whether dlls exist, swap for 32-bit version if needed.
    For Mac and Linux systems compile component model libraries.
    """
    logging.info('Checking internal NRAP-Open-IAM libraries')
    platform = sys.platform
    components_dir = os.path.join('..', 'source', 'components')
    carb_aqu_dir = os.path.join(components_dir, 'aquifer', 'carbonate')
    ca_make_filename = 'make_carbonate.yaml'
    atm_rom_dir = os.path.join(components_dir, 'atmosphere')
    atm_make_filename = 'make_atm_dis_rom.yaml'
    library_name = []
    if platform in ["linux", "linux2"]:
        check_compiler()
        # Linux
        make_from_file(ca_make_filename, carb_aqu_dir, 'linux')
        library_name.append("carbonate.so")
        make_from_file(atm_make_filename, atm_rom_dir, 'linux')
        library_name.append('atmdisrom.so')
    elif platform == "darwin":
        # macOS
        check_compiler()
        make_from_file(ca_make_filename, carb_aqu_dir, 'mac')
        library_name.append("carbonate.dylib")
        make_from_file(atm_make_filename, atm_rom_dir, 'mac')
        library_name.append('atmdisrom.dylib')
    elif platform == "win32":
        # Windows
        cdll = os.path.join(carb_aqu_dir, 'carbonate.dll')
        atmdll = os.path.join(atm_rom_dir, 'atmdisrom.dll')
        if sys.maxsize <= 2**32:  # 32 bit
            logging.debug('32 bit Python detected')
            os.rename(cdll, os.path.join(carb_aqu_dir, 'carbonate_x64.dll'))
            os.rename(os.path.join(carb_aqu_dir, 'carbonate_x32.dll'), cdll)
            os.rename(atmdll, os.path.join(atm_rom_dir, 'atmdisrom_x64.dll'))
            os.rename(os.path.join(atm_rom_dir, 'atmdisrom_x32.dll'), atmdll)
        else:
            cdll64 = os.path.join(carb_aqu_dir, 'carbonate_x64.dll')
            if os.path.exists(cdll64) and not os.path.exists(cdll):
                os.rename(cdll64, cdll)
                logging.debug('Renaming carbonate_x64.dll to carbonate.dll')
        library_name.append("carbonate.dll")
        library_name.append("atmdisrom.dll")
    else:
        library_name.append("carbonate.dll")
        library_name.append("atmdisrom.dll")
        logging.warning('OS platform is not recognized')

    for lib, lib_dir in zip(library_name, (carb_aqu_dir, atm_rom_dir)):
        lib_name = os.path.join(lib_dir, lib)
        if not os.path.exists(lib_name):
            error_msg = ''.join([
                'Cannot find dynamic library:\n {}\n ',
                'Platform: {}']).format(library_name, platform)
            logging.error(error_msg)
        else:
            debug_msg = 'Library {} is found'.format(library_name)
            logging.debug(debug_msg)


def run_tests():
    """
    Run the NRAP-Open-IAM test suite.
    """
    logging.info('Testing NRAP-Open-IAM Installation')
    main_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.join(main_dir, 'test'))
    import iam_test

    runner = unittest.TextTestRunner(verbosity=0)
    test_suite = iam_test.suite('base')
    test_results = runner.run(test_suite)
    if test_results.wasSuccessful():
        message = 'Test suite run successfully, NRAP-Open-IAM is ready to be used.'
        print(message)
        logging.debug(message)
    else:
        message = ''.join([
            'There were errors testing the NRAP-Open-IAM installations, ',
            'please resolve these issues.'])
        logging.warning(message)
        test_results.printErrors()
    return test_results

if __name__ == '__main__':
    __spec__ = None
    setup_logging()
    check_python()
    test_dependencies()
    create_libraries()
    test_results = run_tests()
    # Return non-zero if test suite failed for gitlab-ci
    if not test_results.wasSuccessful():
        sys.exit(1)
