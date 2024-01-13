project(
    'mumpy',
    'cython', 'c',
    license: 'BSD-2-Clause',
    license_files: 'LICENSE',
    meson_version: '>= 1.1.0',
    default_options: [
        'buildtype=release',
        'c_std=c99',
    ],
    version: '0.1.0',
)

py = import('python').find_installation(pure: false)
cc = meson.get_compiler('c')
tempita = files('tempita.py')
fs = import('fs')

incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
    incdir_numpy = run_command(py,
        [
            '-c',
            '''import os
import numpy as np
try:
    incdir = os.path.relpath(np.get_include())
except Exception:
    incdir = np.get_include()
print(incdir)
        '''
        ],
        check: true
    ).stdout().strip()
else
    _incdir_numpy_abs = incdir_numpy
endif
inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)


mumps_names_conda = [
    'mumps_common_seq',
    'dmumps_seq',
    'zmumps_seq',
    'cmumps_seq',
    'smumps_seq',
]

# Use python's ctypes.util.find_library to find the library
# Works for conda and Debian/Ubuntu
mumps_lib_path = fs.parent(run_command(py,
    [
    '-c',
    'import ctypes.util; print(ctypes.util.find_library("@0@"))'.format(mumps_names_conda[0])
    ],
    check: true
).stdout().strip())
mumps_include = include_directories(mumps_lib_path / '..' / 'include')
mumps_libs = []
foreach mumps_name : mumps_names_conda
    mumps_libs += cc.find_library(mumps_name, dirs: mumps_lib_path)
endforeach

subdir('mumpy')