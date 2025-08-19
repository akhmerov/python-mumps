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
    version: run_command('python3', '-m', 'setuptools_scm', check: true).stdout().strip(),
)

py = import('python').find_installation(pure: false)
cc = meson.get_compiler('c')
tempita = files('tempita.py')
fs = import('fs')

# Declare NumPy via meson's dependency() (no explicit method), and add
# the NPY_NO_DEPRECATED_API compile flag.
numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']
np_dep_base = dependency('numpy')
np_dep = declare_dependency(dependencies: [np_dep_base], compile_args: numpy_nodepr_api)

# Declare MUMPS dependencies using meson's dependency() resolution.
# Fail early if a named dependency is not available.
mumps_names = get_option('mumps_names')
mumps_deps = []
foreach mumps_name : mumps_names
    # 1) Try the package as given
    pkg = dependency(mumps_name, required: false)
    if pkg.found()
        mumps_deps += pkg
        continue
    endif

    # 2) If it ends with '_seq', try without that suffix (MPI variant)
    if mumps_name.endswith('_seq')
        short = mumps_name.substring(0, -4)
        pkg2 = dependency(short, required: false)
        if pkg2.found()
            mumps_deps += pkg2
            continue
        endif
    endif

    # Neither found -> error
    error('MUMPS dependency not found: ' + mumps_name + '. Also tried short name if applicable.')
endforeach

subdir('src/mumps')
