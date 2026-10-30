"""HTML report for dso_window.py.

The page itself lives in templates/report.html; this module only injects the
data (as JSON) in place of the /*DATA*/null placeholder.
"""

import json
import os

TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates", "report.html")


def write_report(path, payload):
    with open(TEMPLATE_PATH, encoding="utf-8") as f:
        template = f.read()
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    with open(path, "w", encoding="utf-8") as f:
        f.write(template.replace("/*DATA*/null", data))
