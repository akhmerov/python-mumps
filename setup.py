#!/usr/bin/env python3

# Copyright 2011-2016 Anton Akhmerov, Christoph Groth, and Michael Wimmer
# Copyright 2017 Joseph Weston and Bas Nijholt.
#
# This file is part of mumpy. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# mumpy authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at https://github.com/basnijholt/mumpy

import sys
if sys.version_info[:2] < (3, 4):
    sys.exit('Sorry, Python < 3.4 is not supported')

import os
import re
import subprocess
import configparser
import collections
from setuptools import setup, find_packages, Extension
from distutils.command.build import build
from setuptools.command.build_ext import build_ext

import jinja2 as j2


def configure_extensions(exts, aliases=(), build_summary=None):
    """Modify extension configuration according to the configuration file

    `exts` must be a dict of (name, kwargs) tuples that can be used like this:
    `Extension(name, **kwargs).  This function modifies the kwargs according to
    the configuration file.

    This function modifies `sys.argv`.
    """
    global config_file, config_file_present

    # Determine the name of the configuration file.
    config_file_option = '--configfile'
    # Handle command line option
    for i, opt in enumerate(sys.argv):
        if not opt.startswith(config_file_option):
            continue
        l, _, config_file = opt.partition('=')
        if l != config_file_option or not config_file:
            print('error: Expecting {}=PATH'.format(config_file_option),
                  file=sys.stderr)
            sys.exit(1)
        sys.argv.pop(i)
        break
    else:
        config_file = 'build.conf'

    # Read build configuration file.
    configs = configparser.ConfigParser()
    try:
        with open(config_file) as f:
            configs.read_file(f)
    except IOError:
        config_file_present = False
    else:
        config_file_present = True

    # Handle section aliases.
    for short, long in aliases:
        if short in configs:
            if long in configs:
                print('Error: both {} and {} sections present in {}.'.format(
                    short, long, config_file))
                sys.exit(1)
            configs[long] = configs[short]
            del configs[short]

    # Apply config from file.  Use [DEFAULT] section for missing sections.
    defaultconfig = configs.defaults()
    for name, kwargs in exts.items():
        config = configs[name] if name in configs else defaultconfig
        for key, value in config.items():

            # Most, but not all, keys are lists of strings
            if key == 'language':
                pass
            elif key == 'optional':
                value = bool(int(value))
            else:
                value = value.split()

            if key == 'define_macros':
                value = [tuple(entry.split('=', maxsplit=1))
                         for entry in value]
                value = [(entry[0], None) if len(entry) == 1 else entry
                         for entry in value]

            if key in kwargs:
                msg = 'Caution: user config in file {} shadows {}.{}.'
                if build_summary is not None:
                    build_summary.append(msg.format(config_file, name, key))
            kwargs[key] = value

        kwargs.setdefault('depends', []).append(config_file)
        if config is not defaultconfig:
            del configs[name]

    unknown_sections = configs.sections()
    if unknown_sections:
        print('Error: Unknown sections in file {}: {}'.format(
            config_file, ', '.join(unknown_sections)))
        sys.exit(1)

    return exts


def init_cython():
    global cython_help
    required_cython_version = (0, 22)

    import Cython
    from Cython.Build import cythonize

    # Get Cython version.
    match = re.match('([0-9.]*)(.*)', Cython.__version__)
    cython_version = [int(n) for n in match.group(1).split('.')]
    # Decrease version if the version string contains a suffix.
    if match.group(2):
        while cython_version[-1] == 0:
            cython_version.pop()
        cython_version[-1] -= 1
    cython_version = tuple(cython_version)

    if cython_version < required_cython_version:
        cythonize = None

    return cythonize


def search_libs(libs):
    cmd = ['gcc']
    cmd.extend(['-l' + lib for lib in libs])
    cmd.extend(['-o/dev/null', '-xc', '-'])
    try:
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        pass
    else:
        p.communicate(input=b'int main() {}\n')
        if p.wait() == 0:
            return libs


def search_mumps():
    """Return the configuration for MUMPS if it is available in a known way.

    This is known to work with the MUMPS provided by the Debian package
    libmumps-scotch-dev and the MUMPS binaries in the conda-forge channel."""
    lib_sets = [
        # Debian
        ['{}mumps_scotch'.format(p) for p in 'sdcz']
        + ['mumps_common_scotch', 'mpiseq_scotch'],
        # Conda (via conda-forge).
        # TODO: remove dependency libs (scotch, metis...) when conda-forge
        # packaged mumps/scotch are built as properly linked shared libs
        ['{}mumps'.format(p) for p in 'sdcz']
        + ['mumps_common', 'metis', 'esmumps', 'scotch',
           'scotcherr', 'mpiseq'],
    ]
    common_libs = ['pord', 'gfortran']

    for libs in lib_sets:
        found_libs = search_libs(libs + common_libs)
        if found_libs:
            return found_libs
    return []


def configure_special_extensions(exts, build_summary):
    # Special config for MUMPS.
    mumps = exts['mumpy.mumps']
    if 'libraries' in mumps:
        build_summary.append('User-configured MUMPS')
    else:
        mumps['libraries'] = search_mumps()
        build_summary.append('Auto-configured MUMPS')
    return exts


# Loads version.py module without importing the whole package.
def get_version_and_cmdclass(package_path):
    import os
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location('version',
                                   os.path.join(package_path, '_version.py'))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass


def main():
    mumps = {'mumpy.mumps':
             dict(sources=['mumpy/mumps.pyx'],
                  depends=['mumpy/_mumps.pxd'])}

    # Add NumPy header path to include_dirs of all the extensions.
    import numpy
    numpy_include = numpy.get_include()
    aliases = [('mumps', 'mumpy.mumps')]

    global build_summary

    build_summary = []
    mumps = configure_extensions(mumps, aliases, build_summary)
    mumps = configure_special_extensions(mumps, build_summary)
    mumps['mumpy.mumps'].setdefault('include_dirs', []).append(numpy_include)

    os.chdir('mumpy')
    for fn in (f for f in os.listdir() if f.endswith('.j2')):
        base_name = fn[:-3]
        with open(fn, 'r') as infile, open(base_name, 'w') as outfile:
            outfile.write(j2.Template(infile.read()).render())
    os.chdir('..')

    cythonize = init_cython()
    if cythonize:
        mumps = cythonize([Extension(name, **kwargs)
                          for name, kwargs in mumps.items()],
                          language_level=3,
                          compiler_directives={'linetrace': True})

    version, cmdclass = get_version_and_cmdclass('mumpy')

    classifiers = """\
        Development Status :: 3 - Alpha
        Intended Audience :: Science/Research
        Intended Audience :: Developers
        Programming Language :: Python :: 3 :: Only
        Topic :: Scientific/Engineering
        Operating System :: POSIX
        Operating System :: Unix
        Operating System :: MacOS :: MacOS X
        Operating System :: Microsoft :: Windows"""

    setup(name='mumpy',
          version=version,
          author='mumpy authors',
          author_email='authors@kwant-project.org',
          description=("Python bindings for MUMPS "),
          platforms=["Unix", "Linux", "Mac OS-X", "Windows"],
          url="https://gitlab.kwant-project.org/kwant/mumpy",
          license="BSD",
          packages=find_packages('.'),
          cmdclass={
            'build': build,
            'build_ext': build_ext,
            **cmdclass,
          },
          ext_modules=mumps,
          install_requires=['numpy'],
          classifiers=[c.strip() for c in classifiers.split('\n')])

if __name__ == '__main__':
    main()
    print(build_summary)
