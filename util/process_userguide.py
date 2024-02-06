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
full_specs = re.split(r"\nICNTL\([-\d]+\)", control_params)[1:]
full_specs = [spec for spec in full_specs if not any(desc in spec for desc in exclude)]
# %%
