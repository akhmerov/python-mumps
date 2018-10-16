#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2018 Mumpy Authors.
#
# This file is part of mumpy. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# mumpy authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at
# https://gitlab.kwant-project.org/kwant/mumpy.

from __future__ import print_function

import sys

import re
import os
import glob
import importlib
import subprocess
import configparser
import collections
import textwrap
from setuptools import setup, find_packages, Extension, Command
from distutils.errors import DistutilsError, CCompilerError
from setuptools.command.build_ext import build_ext as build_ext_orig


distr_root = os.path.dirname(os.path.abspath(__file__))


def check_python_version(min_version):
    installed_version = sys.version_info[:3]
    if installed_version < min_version:
        print('Error: Python {} required, but {} is installed'.format(
              '.'.join(map(str, min_version)),
              '.'.join(map(str, installed_version)))
        )
        sys.exit(1)


# Loads version.py module without importing the whole package.
def check_versions(package_path):
    global version, cmdclass, version_is_from_git

    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location('version',
                                   os.path.join(package_path, '_version.py'))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    # Set globals
    version = module.__version__
    cmdclass = module.cmdclass
    version_is_from_git = module.version_is_from_git()



def configure_from_file(exts, aliases=(), build_summary=None):
    """Modify extension configuration according to the configuration file

    `exts` must be a dict of (name, kwargs) tuples that can be used like this:
    `Extension(name, **kwargs).  This function modifies the kwargs according to
    the configuration file.

    This function modifies `sys.argv`.
    """
    global config_file, config_file_present

    #### Determine the name of the configuration file.
    config_file_option = '--configfile'
    # Handle command line option
    for i, opt in enumerate(sys.argv):
        if not opt.startswith(config_file_option):
            continue
        l, _, config_file = opt.partition('=')
        if l != config_file_option or not config_file:
            print('Error: Expecting {}=PATH'.format(config_file_option),
                  file=sys.stderr)
            sys.exit(1)
        sys.argv.pop(i)
        break
    else:
        config_file = 'build.conf'

    #### Read build configuration file.
    configs = configparser.ConfigParser()
    try:
        with open(config_file) as f:
            configs.read_file(f)
    except IOError:
        config_file_present = False
    else:
        config_file_present = True

    #### Handle section aliases.
    for short, long in aliases:
        if short in configs:
            if long in configs:
                print('Error: both {} and {} sections present in {}.'.format(
                    short, long, config_file))
                sys.exit(1)
            configs[long] = configs[short]
            del configs[short]

    #### Apply config from file.  Use [DEFAULT] section for missing sections.
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
    """Set the global variable `cythonize` (and other related globals).

    The variable `cythonize` can be in three states:

    * If Cython should be run and is ready, it contains the `cythonize()`
      function.

    * If Cython is not to be run, it contains `False`.

    * If Cython should, but cannot be run it contains `None`.  A help message
      on how to solve the problem is stored in `cython_help`.

    This function modifies `sys.argv`.
    """
    global cythonize, cython_help

    cython_option = '--cython'
    required_cython_version = (0, 24)
    try:
        sys.argv.remove(cython_option)
        cythonize = True
    except ValueError:
        cythonize = version_is_from_git

    if cythonize:
        try:
            import Cython
            from Cython.Build import cythonize
        except ImportError:
            cythonize = None
        else:
            #### Get Cython version.
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

        if cythonize is None:
            msg = ("Install Cython >= {0} or use"
                    " a source distribution (tarball) of Mumpy.")
            ver = '.'.join(str(e) for e in required_cython_version)
            cython_help = msg.format(ver)
    else:
        msg = "Run setup.py with the {} option to enable Cython."
        cython_help = msg.format(cython_option)


def init_jinja():
    """Set the global variable `jinjize`.

    The variable `jinjize` can be in two states:

    * If Jinja should be run and is ready, it contains a function that
      takes a template as a string and returns the rendered template as
      a string.

    * If Jinja is not to be run, it contains `False`.

    If Jinja should, but cannot be run it calls 'sys.exit()' with an
    error message.
    """

    global jinjize

    if not version_is_from_git:
        jinjize = False
    else:
        try:
            import jinja2
        except ImportError:
            msg = ("Install Jinja2 or use"
                    " a source distribution (tarball) of Mumpy.")
            sys.exit(msg)
        else:
            jinjize = lambda t: jinja2.Template(t).render()


def banner(title=''):
    starred = title.center(79, '*')
    return '\n' + starred if title else starred


class build_ext(build_ext_orig):
    def run(self):
        if not config_file_present:
            # Create an empty config file if none is present so that the
            # extensions will not be rebuilt each time.  Only depending on the
            # config file if it is present would make it impossible to detect a
            # necessary rebuild due to a deleted config file.
            with open(config_file, 'w') as f:
                f.write('# Build configuration created by setup.py '
                        '- feel free to modify.\n')

        try:
            super().run()
        except (DistutilsError, CCompilerError):
            error_msg = self.__error_msg.format(
                header=banner(' Error '), sep=banner())
            print(error_msg.format(file=config_file, summary=build_summary),
                  file=sys.stderr)
            raise
        print(banner(' Build summary '), *build_summary, sep='\n')
        print(banner())

    __error_msg = textwrap.dedent("""\
        {header}
        The compilation of Mumpy has failed.  Please examine the error message
        above and consult the installation instructions in README.rst.
        You might have to customize {{file}}.

        Build configuration was:

        {{summary}}
        {sep}""")


def long_description():
    text = []
    try:
        with open('README.rst', encoding='utf8') as f:
            for line in f:
                if line.startswith('See also in this directory:'):
                    break
                text.append(line.rstrip())
            while text[-1] == "":
                text.pop()
    except:
        return ''
    return '\n'.join(text)


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
        ['zmumps_scotch', 'mumps_common_scotch', 'mpiseq_scotch'],
        # Conda (via conda-forge).
        # TODO: remove dependency libs (scotch, metis...) when conda-forge
        # packaged mumps/scotch are built as properly linked shared libs
        ['zmumps', 'mumps_common', 'metis', 'esmumps', 'scotch',
         'scotcherr', 'mpiseq', 'blas'],
    ]
    common_libs = ['pord', 'gfortran']

    for libs in lib_sets:
        found_libs = search_libs(libs + common_libs)
        if found_libs:
            return found_libs
    return []


def configure_from_defaults(exts, build_summary):
    #### Special config for MUMPS.
    mumps = exts['mumpy.mumps']
    if 'libraries' in mumps:
        build_summary.append('User-configured MUMPS')
    else:
        mumps_libs = search_mumps()
        if mumps_libs:
            mumps['libraries'] = mumps_libs
            build_summary.append('Auto-configured MUMPS')
        else:
            sys.exit('MUMPS not found')

    return exts


def maybe_add_numpy_include(exts):
    # Add NumPy header path to include_dirs of all the extensions.
    try:
        import numpy
    except ImportError:
        print(banner(' Caution '), 'NumPy header directory cannot be determined'
              ' ("import numpy" failed).', banner(), sep='\n', file=sys.stderr)
    else:
        numpy_include = numpy.get_include()
        for ext in exts.values():
            ext.setdefault('include_dirs', []).append(numpy_include)
    return exts


def maybe_cythonize(exts):
    """Prepare a list of `Extension` instances, ready for `setup()`.

    The argument `exts` must be a mapping of names to kwargs to be passed
    on to `Extension`.

    If Cython is to be run, create the extensions and calls `cythonize()` on
    them.  If Cython is not to be run, replace .pyx file with .c or .cpp,
    check timestamps, and create the extensions.
    """
    if cythonize:
        return cythonize([Extension(name, **kwargs)
                          for name, kwargs in exts.items()],
                         language_level=3,
                         compiler_directives={'linetrace': True})

    # Cython is not going to be run: replace pyx extension by that of
    # the shipped translated file.

    result = []
    problematic_files = []
    for name, kwargs in exts.items():
        language = kwargs.get('language')
        if language is None:
            ext = '.c'
        elif language == 'c':
            ext = '.c'
        elif language == 'c++':
            ext = '.cpp'
        else:
            print('Unknown language: {}'.format(language), file=sys.stderr)
            sys.exit(1)

        pyx_files = []
        cythonized_files = []
        sources = []
        for f in kwargs['sources']:
            if f.endswith('.pyx'):
                pyx_files.append(f)
                f = f.rstrip('.pyx') + ext
                cythonized_files.append(f)
            sources.append(f)
        kwargs['sources'] = sources

        # Complain if cythonized files are older than Cython source files.
        try:
            cythonized_oldest = min(os.stat(f).st_mtime
                                    for f in cythonized_files)
        except OSError:
            msg = "Cython-generated file {} is missing."
            print(banner(" Error "), msg.format(f), "",
                  cython_help, banner(), sep="\n", file=sys.stderr)
            sys.exit(1)

        for f in pyx_files + kwargs.get('depends', []):
            if f == config_file:
                # The config file is only a dependency for the compilation
                # of the cythonized file, not for the cythonization.
                continue
            if os.stat(f).st_mtime > cythonized_oldest:
                problematic_files.append(f)

        result.append(Extension(name, **kwargs))

    if problematic_files:
        msg = ("Some Cython source files are newer than files that have "
               "been derived from them:\n{}")
        msg = msg.format(", ".join(problematic_files))

        # Cython should be run but won't.  Signal an error if this is because
        # Cython *cannot* be run, warn otherwise.
        error = cythonize is None
        if cythonize is False:
            dontworry = ('(Do not worry about this if you are building Mumpy '
                         'from unmodified sources,\n'
                         'e.g. with "pip install".)\n\n')
            msg = dontworry + msg

        print(banner(" Error " if error else " Caution "), msg, "",
              cython_help, banner(), sep="\n", file=sys.stderr)
        if error:
            sys.exit(1)

    return result


def maybe_jinjize(exts):
    """Run Jinja on any templated source files in 'exts'.

    The argument `exts` must be a mapping of names to kwargs to be passed
    on to `Extension`.

    For each source or dependency file with a '.j2' extension we: remove
    the '.j2' extension from the specification in 'exts', and run Jinja
    on the file, producing a file without the '.j2' extension (if Jinja
    is to be run).
    """

    to_jinjize = []

    def convert(source):
        target, ext = os.path.splitext(source)
        if ext == '.j2':
            to_jinjize.append((source, target))
            return target
        else:
            return source

    def should_rebuild(source_target):
        source, target = source_target
        return (not os.path.exists(target)
                or os.path.getmtime(source) > os.path.getmtime(target))

    # Remove '.j2' extension from files that have it
    for ext, kwargs in exts.items():
        for key in ('sources', 'depends'):
            if key not in kwargs:
                continue
            kwargs[key] = [convert(s) for s in kwargs[key]]

    # Jinjize if required
    if jinjize:
        for source, target in filter(should_rebuild, to_jinjize):
            try:
                with open(source, 'r') as src, open(target, 'w') as dst:
                    dst.write(jinjize(src.read()))
            except Exception:
                sys.exit('Failed to render template "{}"'.format(source))
            else:
                print('Rendered template {} → {} because it changed'
                      .format(source, target))

    return exts


def main():
    check_python_version((3, 5))
    check_versions('mumpy')
    cmdclass.update(build_ext=build_ext)

    exts = {
        'mumpy.mumps': dict(
             sources=['mumpy/mumps.pyx.j2'],
             depends=['mumpy/_mumps.pxd.j2'],
         ),
    }

    aliases = [('mumps', 'mumpy.mumps')]

    init_cython()
    init_jinja()

    global build_summary
    build_summary = []
    exts = configure_from_file(exts, aliases, build_summary)
    exts = configure_from_defaults(exts, build_summary)
    exts = maybe_add_numpy_include(exts)
    exts = maybe_jinjize(exts)
    exts = maybe_cythonize(exts)

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
          description="Bindings for the MUMPS sparse solver",
          long_description=long_description(),
          platforms=["Unix", "Linux", "Mac OS-X", "Windows"],
          url="https://gitlab.kwant-project.org/kwant/mumpy",
          license="BSD",
          packages=find_packages('.'),
          cmdclass=cmdclass,
          ext_modules=exts,
          install_requires=['numpy', 'scipy'],
          classifiers=[c.strip() for c in classifiers.split('\n')])

if __name__ == '__main__':
    main()
