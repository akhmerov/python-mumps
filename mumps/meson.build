templates = files('_mumps.pxd.in', '_mumps.pyx.in')

pyx_files = []
foreach template : templates
    base_name = fs.name(template)
    rendered_name = base_name.replace('.in', '')

    pyx_file = custom_target(
        rendered_name,
        input: template,
        output: rendered_name,
        command: [py, tempita, '@INPUT@', '@OUTPUT@']
    )
    pyx_files += pyx_file
endforeach

py_files = files('__init__.py', 'mumps.py', 'fortran_helpers.py', '_version.py')
py.install_sources(py_files, subdir: 'mumps')

foreach py_file : py_files
    fs.copyfile(py_file)
endforeach

py.extension_module('_mumps',
    pyx_files,
    dependencies: [np_dep, mumps_libs],
    install: true,
    subdir: 'mumps'
)

subdir('tests')
