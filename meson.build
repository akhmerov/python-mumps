project(
    'python-mumps',
    'cython', 'c',
    license: 'BSD-2-Clause',
    license_files: 'LICENSE',
    meson_version: '>= 1.1.0',
    default_options: [
        'buildtype=release',
        'c_std=c99',
    ],
    version: run_command('mumps/_version.py', check: true).stdout().strip(),
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
numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']
np_dep = declare_dependency(include_directories: inc_np, compile_args: numpy_nodepr_api)
mumps_names = get_option('mumps_names')

mumps_libs = []
foreach mumps_name : mumps_names
    if cc.find_library(mumps_name, required: false).found()
        mumps_libs += cc.find_library(mumps_name)
    else
        # Fall back to the name without '_seq'.
        mumps_libs += cc.find_library(mumps_name.substring(0, -4))
    endif
endforeach

subdir('mumps')
