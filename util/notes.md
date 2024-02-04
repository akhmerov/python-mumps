# Dev notes

## Generation of Enum objects

- Download the latest user guide from [the
  website](https://mumps-solver.org/index.php?page=doc)
- As of the time of writing, it is only available as a
  [PDF](https://mumps-solver.org/doc/userguide_5.6.2.pdf)
- Use `pdftotext` to convert the PDF to text
- Copy the contents of the section 6.1
- Use github copilot api to generate the enum objects. To do this, clone
  [copilot-api](https://github.com/B00TK1D/copilot-api) and follow the
  instructions to launch the server.
- Run the code in `process_userguide.py` to generate the rough version of the
  enum objects
- Manually check the variable names and formatting of the resulting code.
- Add the generated code to `constants.py`.
