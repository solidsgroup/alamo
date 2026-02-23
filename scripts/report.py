import os
from pathlib import Path
from html import escape

HTML_REPORT = "report.html"

def init_html(title="Regression Test Report"):
    with open(HTML_REPORT, "w") as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
  <title>{escape(title)}</title>
  <style>
    body {{ font-family: sans-serif; margin: 20px; }}
    .test-entry {{ border-bottom: 1px solid #ccc; margin-bottom: 10px; padding-bottom: 10px; }}
    .status {{
      display: inline-block;
      font-weight: bold;
      padding: 2px 6px;
      border-radius: 4px;
      color: white;
    }}
    .pass {{ background-color: #4CAF50; }}
    .fail {{ background-color: #f44336; }}
    .timeout {{ background-color: #999999; }}
    .warn {{ background-color: #ff9800; }}
    .thumbnails {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 5px;
    }}
    .thumbnail {{
      width: 200px;
      border: 1px solid #ccc;
    }}
  </style>
</head>
<body>
<h1>{escape(title)}</h1>
""")

def append_html(record):
    path = Path(record['path'])
    run_status = record['runStatus'].lower()
    check_status = record.get('checkStatus', '').lower()
    testname = escape(record.get('test-section', path.name))
    stdout_path = path / "stdout"
    stderr_path = path / "stderr"

    # Determine status classes
    status_class = "fail"
    if run_status == "pass":
        status_class = "pass"
    elif run_status == "timeout":
        status_class = "timeout"

    check_class = ""
    if check_status == "pass":
        check_class = "pass"
    elif check_status == "fail":
        check_class = "fail"
    elif check_status == "warn":
        check_class = "warn"

    with open(HTML_REPORT, "a") as f:
        f.write(f"<div class='test-entry'>\n")
        f.write(f"<h3>{testname}</h3>\n")
        f.write(f"<p>Run: <span class='status {status_class}'>{run_status.upper()}</span>\n")
        if check_status:
            f.write(f" &nbsp; Check: <span class='status {check_class}'>{check_status.upper()}</span>")
        f.write("</p>\n")
        f.write(f"<p><a href='../{stdout_path}'>stdout</a> | <a href='../{stderr_path}'>stderr</a></p>\n")

        # Images
        image_extensions = {".svg", ".png", ".jpg", ".jpeg"}
        images = [f for f in path.iterdir() if f.suffix.lower() in image_extensions]
        if images:
            f.write("<div class='thumbnails'>\n")
            for img in images:
                f.write(f"<a href='{img}'><img src='{img}' class='thumbnail'></a>\n")
            f.write("</div>\n")

        f.write("</div>\n")

        f.write("</div>\n")

def finalize_html():
    with open(HTML_REPORT, "a") as f:
        f.write("</body></html>\n")
