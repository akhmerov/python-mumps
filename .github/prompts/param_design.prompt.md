# Prompt: Convert a MUMPS snippet into a python-mumps parameter definition

Goal
- Edit the parameter class that immediately follows the snippet, matching the conventions in `src/mumps/enums.py`.

Output format (MUST)
- The `@param(index=N, page=M)` decorator is ALREADY present above the class and MUST NOT be modified or repeated.
- Produce ONLY the class definition that follows the decorator: `class <snake_case_name>:` with its body.
- You MAY rename the class to a better snake_case name if needed.
- Include a docstring. For discrete-meaning parameters, the docstring MUST NOT contain any integers or per-value meanings—refer only to the member names. For continuous/arbitrary parameters, the docstring is the only content.
- For discrete-meaning parameters, you MUST NOT mention in the docstring that the user should choose from the members; this is implied.
- For discrete integer parameters: define members as `name = int_value, "short description"`. The description MUST NOT mention the integer value.
- For continuous/float or arbitrary-valued parameters: define NO members; docstring only.
- For continuous/float or arbitrary-valued parameters: you MUST mention parameter values and describe their meaning in the docstring.
- Output exactly one class block (no decorator line, no imports, no extra prose).

Decision rules (MUST)
1) If the parameter has a finite set of distinct meanings → create members (enum-like).
2) Otherwise (floats, ranges, arbitrary inputs) → no members; docstring-only.
3) When multiple numeric forms map to the same meaning (e.g., `<1` and `0`) → choose the simplest explicit integer (e.g., `0`). If this normalization affects the default, mention the normalized default in the docstring.

Docstring requirements (MUST)
- Summarize role, default, constraints, typical ranges, and interactions with other parameters.
- Do NOT list per-value meanings or integers when members exist. Refer to member names only (e.g., "Default is `errors_and_warnings`").
- For continuous parameters, include units when known and any guidance on typical ranges and effects.
- Do NOT mention anything not present in the snippet.

Naming and style (MUST)
- Class name: short, descriptive, snake_case.
- Member names: lowercase with underscores; concise.
- Member descriptions: short fragments.

Information preservation (MUST)
- Preserve all user-relevant information from the snippet. Put general behavior, defaults, constraints, ranges, and cross-parameter interactions into the docstring. Put per-value meanings into member descriptions (when members exist).
- Keep references to other MUMPS parameters (ICNTL/CNTL/DKEEP/KEEP/INFO, etc.) in the docstring.

Hard constraints (DO NOT)
- Do NOT mention integer values or per-value meanings in the docstring when members exist.
- Do NOT emit more than one class unless explicitly requested.
- Do NOT add imports, module headers, comments, or explanatory prose outside the class.
- Do NOT change the numeric values from the snippet except for allowed normalization in Decision rule (3).
- Do NOT modify or re-emit the `@param(...)` decorator; keep it as-is in the file.

Example
- Input (excerpt):
  # ICNTL(4) is used by ZMUMPS to control printing of error,
  #    warning, and diagnostic messages. It has default value 2.
  #    Possible values are:
  #
  #   <1       __No messages output.
  #    1       __Only error messages printed.
  #    2       __Errors and warnings printed.
  #    3       __Errors and warnings and terse diagnostics
  #               (only first ten entries
  #              of arrays printed).
  #    4       __Errors and warnings and all information
  #              on input and output parameters printed.

- Expected output:

@param(index=4)
class message_level:
    """
    Printing/verbosity level

    Default: Errors and warnings are printed (`errors_and_warnings`)
    """

    no_messages = 0, "no messages"
    errors_only = 1, "only error messages printed"
    errors_and_warnings = 2, "errors and warnings printed"
    terse_diagnostics = 3, "errors, warnings, and terse diagnostics (first ten array entries)"
    verbose_all = 4, "all information on input and output parameters printed"

Success criteria
- The class can replace the existing placeholder class body in `src/mumps/enums.py` (decorator preserved) and follows the project's style.
- All necessary information is present: docstring for behavior/defaults/constraints/interactions; member descriptions for per-value meanings.
