# %%
from pathlib import Path
from subprocess import call
import re
from tempfile import NamedTemporaryFile

import requests
import jinja2
import clipboard

VERSION = "5.6.2"
COPILOT_API_PORT = 8000

# Download the user guide if the file does not exist
target = Path(f"userguide_{VERSION}.txt")
if not target.is_file():
    url = f"https://mumps-solver.org/doc/userguide_{VERSION}.pdf"
    data = requests.get(url).content
    # Save it to a temporary file
    with NamedTemporaryFile(delete_on_close=False) as source:
        source.write(data)
        source.close()
        # Convert to txt
        call(["pdftotext", "-layout", source.name, source.name + ".txt"])
        data = Path(source.name + ".txt").read_text()
        Path(source.name + ".txt").unlink()

    # Remove page numbers; they are marked by number followed by form feed
    # character
    data = re.sub(r"\s*\d+\s*\f", "\n", data)
    target.write_text(data)
else:
    data = target.read_text()

# %%

# Get the control parameters section; may need to be adjusted for different
# versions of the user guide
START_STR = "6.1     Integer control parameters"
END_STR = "6.2     Real/complex control parameters"
start = data.find(START_STR)
end = data.find(END_STR, start)
control_params = data[start + len(START_STR) : end]
param_def_regex = re.compile(r"• ICNTL\(([-\d]+)\) (.*)?\n")
param_defs = param_def_regex.findall(control_params)
exclude = (
    "reserved in current version",
    "not used in current version",
)
param_defs = [param for param in param_defs if param[1] not in exclude]
param_defs = {int(param[0]): param[1] for param in param_defs}
url = f"http://localhost:{COPILOT_API_PORT}/api"


# %%
query = """
# Generate a name of the control parameter given the description.
# Example:
#
# Description: "is the output stream for error messages" -> Name: "ERROR_STREAM"
#
"""
names = []
for param, desc in param_defs.items():
    snippet = query + f"# Description: {desc}\n# Name:"
    response = requests.post(url, json={"prompt": snippet, "language": "python"})
    names.append(re.sub(r"[\s\"]", "", response.text))

# %%
enum_template = jinja2.Template(
    """
class {{ enum_name }}(IntEnum):
    description = dict()
    {% for param, desc, name in data %}
    {{ name }} = {{ param }}
    description[{{ param }}] = "{{ desc }}"
    {% endfor %}
"""
)
# %%
clipboard.copy(
    enum_template.render(
        enum_name="ICNTL",
        data=zip(param_defs.keys(), param_defs.values(), names),
    )
)

# %%

# Read actual parameter descriptions.
full_specs = re.split(r"\nICNTL\(([-\d]+)\)", control_params)[1:]
# Zip together
full_specs = zip(full_specs[::2], full_specs[1::2])
# Skip unused parameters
full_specs = [
    spec for spec in full_specs if not any(desc in spec[1] for desc in exclude)
]
full_specs = {int(param): {"desc": spec.strip()} for param, spec in full_specs}
# Extract possible values. These come after "Possible values:" and typically
# end before "Other values" or "Default value".

# The regex below is tested to be good enough to capture all enumerated options.
possible_values_regex = re.compile(
    r"^\s*Possible values\s*?:\s?^(.*)?^\s*?"
    r"(Other values|Default value|Values different)",
    re.DOTALL | re.MULTILINE,
)

# Within each list of possible values, we want to capture the values and their
# descriptions. These follow a pattern <value>: <description>. Values are not
# always specified as integers, which makes parsing complicated. Sometimes it's
# an letter or an enumeration.
value_regex = re.compile(r"^\s*(.{1,7})\s?:", re.MULTILINE)

for param, spec in full_specs.items():
    match = possible_values_regex.search(spec["desc"])
    if match:
        values_and_descriptions = re.split(value_regex, match.group(1))[1:]
        values = {
            value: desc
            for value, desc in zip(
                values_and_descriptions[::2], values_and_descriptions[1::2]
            )
        }
        spec["values"] = values
    else:
        print("Possible values not found for", param)
        spec["values"] = None

# %%
