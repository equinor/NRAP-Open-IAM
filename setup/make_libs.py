# -*- coding: utf-8 -*-
"""
The module is build to replace the need for the make utility to compile the
necessary OpenIAM component models.  This code will not duplicate the full
functionality of make, just the necessary functionality for the OpenIAM
component model compiling.

Currently the compiler is hard coded as gfortran.

Created on Thu Feb 01 09:37:03 2018

@author: Seth King
"""
import yaml
import os


def make_from_file(make_filename, directory, build):
    """
    Compile model using the specified setup.

    Compile model given yaml input filename, the directory of the model,
    and the build option.

    :param make_filename: yaml formatted file with build options for compiling
        the model. The file must include an entry named objects which lists the
        source filenames minus the .f90 extension (.f90 is a hard coded assumption).
        The file must also include a entry for each build option to be called with.
        The build option can specify a command to run with cmd or 'flags'
        for global flags, 'comp_flags' for compiling flags, 'link_flags'
        for linking flags and an entry for library for the name of the library to
        be compiled.
    :type make_filename: str

    :param directory: Directory where source files are located and build will occur.
        The code will be ran from this directory and change back to the run directory
        before execution is finished.
    :type directory: str

    :param build: The build option to execute; generally, these are: 'windows',
        'mac', 'linux', or 'clean'
    :type build: str
    """
    original_directory = os.getcwd()
    os.chdir(directory)
    with open(make_filename, 'r') as mkf:
        make_data = yaml.load(mkf, Loader=yaml.SafeLoader)
    make_from_data(make_data, build)

    os.chdir(original_directory)
    return


def make_from_data(make_data, build):
    """
    Compiles model given a dictionary of build options and the build option.

    :param make_data: dictionary with build options for compiling
        the model. The dictionary must include an entry named objects which lists the
        source filenames minus the .f90 extension (.f90 is a hard coded assumption).
        The dictionary must also include a entry for each build option to be called with.
        The build option can specify a command to run with cmd or 'flags'
        for global flags, 'comp_flags' for compiling flags, 'link_flags'
        for linking flags and an entry for library for the name of the library to
        be compiled.
    :type make_data: dict

    :param build: The build option to execute; generally, these are: 'windows',
        'mac', 'linux', or 'clean'
    :type build: str
    """
    if build.lower() not in make_data:
        raise AssertionError('Build option {build} not found in data'.format(build))
    build_data = make_data[build]
    objs = make_data['objects']
    if 'cmd' in build_data:
        cmd = build_data['cmd']
        if type(cmd) == str:
            os.system(cmd)
        elif type(cmd) == list:
            for c in cmd:
                os.system(c)
        else:
            raise AssertionError('cmd data type not understood for {}'.format(cmd))
    else:
        flags = build_flags(['comp_flags', 'flags'], build_data)
        for obj in objs:
            os.system('gfortran {flags} -c {obj}.f90'.format(flags=flags, obj=obj))
        flags = build_flags(['link_flags', 'flags'], build_data)
        obj_fns = ' '.join([obj + '.o' for obj in objs])
        os.system('gfortran {flags} -o {lib} {objs}'.format(flags=flags,
                  lib=build_data['library'],
                  objs=obj_fns))
    return


def build_flags(flg_lst, data):
    """
    Outputs a flag string from the data dictionary from a list of optional
    dictionary keywords in the input flg_lst.
    """
    flags = ''
    for fg in flg_lst:
        if (fg in data) and data[fg]:
            flags += ' ' + data[fg]
    return flags

if __name__ == '__main__':
    # directory = os.path.join('..', 'source', 'components', 'aquifer', 'carbonate')
    # make_filename = 'make_carbonate.yaml'
    # make_from_file(make_filename, directory, 'clean')
    # make_from_file(make_filename, directory, 'windows')
    directory = os.path.join('..', 'source', 'components', 'atmosphere')
    make_filename = 'make_atm_dis_rom.yaml'
    make_from_file(make_filename, directory, 'clean')
    make_from_file(make_filename, directory, 'windows')
