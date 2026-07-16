import os
from pathlib import Path
from html import escape
import json

HTML_SUMMARY = "report.html"
HTML_REPORT = None
JSON_REPORT = None

def html_hdr(title):
    return f"""<!DOCTYPE html>
<html>
<head>
  <title>{escape(title)}</title>
  <meta charset="UTF-8">
  <style>
    * {{ box-sizing: border-box; }}
    body {{ 
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      margin: 0;
      padding: 0;
      background: #f8f9fa;
    }}
    .container {{
      max-width: 1400px;
      margin: 0 auto;
      padding: 20px;
    }}
    .header {{
      background: white;
      padding: 30px;
      margin-bottom: 20px;
      border: 1px solid #dee2e6;
      border-radius: 4px;
    }}
    .header h1 {{
      margin: 0;
      font-size: 28px;
      font-weight: 500;
      color: #212529;
    }}
    .content {{
      background: white;
      padding: 30px;
      border: 1px solid #dee2e6;
      border-radius: 4px;
    }}
    h2 {{
      color: #212529;
      font-size: 20px;
      margin: 0 0 20px 0;
      font-weight: 500;
      padding-bottom: 10px;
      border-bottom: 2px solid #dee2e6;
    }}
    .runs-list {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
      gap: 10px;
      margin-bottom: 30px;
    }}
    .run-link {{
      display: block;
      padding: 12px 16px;
      background: white;
      border: 1px solid #dee2e6;
      border-radius: 4px;
      text-decoration: none;
      color: #0d6efd;
      transition: all 0.15s;
      font-weight: 400;
    }}
    .run-link:hover {{
      background: #0d6efd;
      color: white;
      border-color: #0d6efd;
    }}
    .run-link.missing {{
      color: #6c757d;
      cursor: not-allowed;
    }}
    .run-link.missing:hover {{
      background: white;
      color: #6c757d;
      border-color: #dee2e6;
    }}
    .summary-table-wrapper {{
      overflow-x: auto;
      margin-top: 20px;
      border: 1px solid #dee2e6;
      border-radius: 4px;
    }}
    .summary-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
      background: white;
    }}
    .summary-table th {{
      background: #f8f9fa;
      color: #495057;
      padding: 12px 8px;
      text-align: left;
      font-weight: 600;
      border-bottom: 2px solid #dee2e6;
      white-space: nowrap;
      width: 120px;
    }}
    .summary-table th.test-header {{
      text-align: center;
      vertical-align: top;
      white-space: normal;
      line-height: 1.3;
    }}
    .test-family {{
      font-weight: 600;
      font-size: 13px;
      display: block;
      margin-bottom: 3px;
    }}
    .test-section {{
      font-weight: 400;
      font-size: 12px;
      color: #6c757d;
      display: block;
    }}
    .summary-table th:first-child {{
      position: sticky;
      left: 0;
      z-index: 11;
      min-width: 250px;
      background: #f8f9fa;
      border-right: 2px solid #dee2e6;
    }}
    .summary-table td {{
      padding: 10px 8px;
      border-bottom: 1px solid #dee2e6;
    }}
    .summary-table tr:hover td {{
      background: #f8f9fa;
    }}
    .summary-table td:first-child {{
      position: sticky;
      left: 0;
      background: white;
      border-right: 1px solid #dee2e6;
      z-index: 5;
    }}
    .summary-table tr:hover td:first-child {{
      background: #f8f9fa;
    }}
    .testid-cell {{
      display: flex;
      align-items: center;
      gap: 10px;
      justify-content: space-between;
    }}
    .testid-link {{
      color: #0d6efd;
      text-decoration: none;
      font-weight: 500;
      flex: 1;
    }}
    .testid-link:hover {{
      text-decoration: underline;
    }}
    .result-badges {{
      display: flex;
      gap: 5px;
      flex-wrap:wrap;
    }}
    .result-badge {{
      display: block;
      padding: 3px 8px;
      border-radius: 3px;
      font-size: 12px;
      font-weight: 600;
      line-height: 1;
    }}
    .result-badge.pass {{
      background: #FFFFFF00;
      color: #0f5132;
      border: 1px solid #badbcc;
    }}
    .result-badge.check.pass {{
      background: #d1e7dd;
    }}
    .result-badge.fail {{
      background: #FFFFFF00;
      color: #842029;
    }}
    .result-badge.check.fail {{
      background: #f8d7da;
    }}
    .result-badge.warn {{
      background: #FFFFFF00;
      color: #664d03;
    }}
    .result-badge.check.warn {{
      background: #fff3cd;
    }}
    .result-badge.timeout {{
      background: #FFFFFF00;
      color: #41464b;
    }}
    .result-badge.check.timeout {{
      background: #e2e3e5;
    }}
    .cell {{
      text-align: center;
      width: 120px;
      min-width: 120px;
      max-width: 120px !important;
    }}
    .cell-content {{
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 4px;
    }}
    .status-badge {{
      display: inline-block;
      padding: 4px 10px;
      border-radius: 3px;
      font-weight: 600;
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.3px;
      min-width: 65px;
    }}
    .status-PASS {{ 
      background: #FFFFFF00;
      color: #0f5132;
      border: 1px solid #badbcc;
    }}
    .check.status-PASS {{ 
      background: #d1e7dd;
    }}
    .status-FAIL {{ 
      background: #FFFFFF00;
      color: #842029;
      border: 1px solid #f5c2c7;
    }}
    .check.status-FAIL {{ 
      background: #f8d7da;
    }}
    .status-WARN {{ 
      background: #FFFFFF00;
      color: #664d03;
      border: 1px solid #ffecb5;
    }}
    .check.status-WARN {{ 
      background: #fff3cd;
    }}
    .status-TIMEOUT {{ 
      color: #41464b;
      border: 1px solid #d3d6d8;
    }}
    .check.status-TIMEOUT {{ 
      background: #e2e3e5;
    }}
    .runtime {{
      font-size: 11px;
      color: #6c757d;
      font-weight: 500;
    }}
    .runtime-bar {{
      width: 100%;
      height: 4px;
      background: #e9ecef;
      border-radius: 2px;
      overflow: hidden;
      border: 1px solid #dee2e6;
    }}
    .runtime-fill {{
      height: 100%;
      background: #0d6efd;
      transition: width 0.3s;
    }}
    .empty-cell {{
      background: #f8f9fa;
      color: #adb5bd;
    }}
    .legend {{
      display: flex;
      gap: 20px;
      margin-top: 20px;
      padding: 15px 20px;
      background: #f8f9fa;
      border: 1px solid #dee2e6;
      border-radius: 4px;
      flex-wrap: wrap;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      gap: 8px;
      font-size: 14px;
    }}
    .legend-color {{
      width: 18px;
      height: 18px;
      border-radius: 3px;
      border: 1px solid #dee2e6;
    }}
    .test-entry {{ 
      border-bottom: 1px solid #dee2e6;
      margin-bottom: 15px;
      padding-bottom: 15px;
    }}
    .status {{
      display: inline-block;
      font-weight: 600;
      padding: 3px 8px;
      border-radius: 3px;
      font-size: 12px;
    }}
    .pass {{ 
      background: #d1e7dd;
      color: #0f5132;
    }}
    .fail {{ 
      background: #f8d7da;
      color: #842029;
    }}
    .timeout {{ 
      background: #e2e3e5;
      color: #41464b;
    }}
    .warn {{ 
      background: #fff3cd;
      color: #664d03;
    }}
    .thumbnails {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 5px;
    }}
    .thumbnail {{
      width: 200px;
      border: 1px solid #dee2e6;
      border-radius: 4px;
    }}
  </style>
</head>
<body>
<div class="container">
"""

def html_ftr():
    return """</div>
</body>
</html>"""

def init_html(html_path, title="Regression Test Report"):
    global HTML_REPORT, JSON_REPORT
    HTML_REPORT = html_path
    JSON_REPORT = html_path.replace(".html", ".json")
    with open(HTML_REPORT, "w") as f:
        f.write(html_hdr(title))
        f.write(f"<div class='header'><h1>{escape(title)}</h1></div>\n")
        f.write("<div class='content'>\n")

def write_master_html():
    """
    Scan all report/*.json files and generate report/report.html
    with:
      - list of all runs
      - dense matrix of testid x test-section with PASS/FAIL/etc.
    """
    report_dir = Path("report")
    json_files = sorted(report_dir.glob("*.json"))
    
    # Load all JSON records keyed by testid
    db = {}
    for jf in json_files:
        try:
            rec = json.loads(jf.read_text())
        except Exception:
            continue
        for test_section, record in rec.items():
            testid = record.get("testid", "UNKNOWN")
            if testid not in db:
                db[testid] = {}
            db[testid][test_section] = record
    
    # Collect test names and testids
    all_tests = sorted({
        test_section
        for run in db.values()
        for test_section in run.keys()
    })
    testids = sorted(db.keys())
    
    # Calculate max runtime for each test for normalization
    max_runtimes = {}
    for test_section in all_tests:
        max_time = 0
        for testid in testids:
            rec = db[testid].get(test_section)
            if rec and "executionTime" in rec:
                try:
                    time = float(rec["executionTime"])
                    max_time = max(max_time, time)
                except:
                    pass
        max_runtimes[test_section] = max_time if max_time > 0 else 1
    
    # Build HTML
    html = []
    html.append(html_hdr("Regression Test Summary"))
    html.append("<div class='header'><h1>Regression Test Summary</h1></div>")
    html.append("<div class='content'>")
    
    # List of runs
    html.append("<h2>Test Runs</h2>")
    html.append("<div class='runs-list'>")
    for testid in testids:
        html_file = f"report/{testid}.html"
        if (report_dir / f"{testid}.html").exists():
            html.append(f"<a href='{html_file}' class='run-link'>{escape(testid)}</a>")
        else:
            html.append(f"<div class='run-link missing'>{escape(testid)} (missing)</div>")
    html.append("</div>")
    
    # Summary table
    html.append("<h2>Test Results Matrix</h2>")
    html.append("<div class='summary-table-wrapper'>")
    html.append("<table class='summary-table'>")
    html.append("<thead><tr><th>Test ID</th>")
    for t in all_tests:
        # Split test name into family and section
        if '/' in t:
            parts = t.rsplit('/', 1)
            family = parts[0]
            section = parts[1]
        else:
            family = t
            section = ""
        
        html.append(f"<th class='test-header' title='{escape(t)}'>")
        html.append(f"<span class='test-family'>{escape(family)}</span>")
        if section:
            html.append(f"<span class='test-section'>{escape(section)}</span>")
        html.append("</th>")
    html.append("</tr></thead><tbody>")
    
    for testid in reversed(testids):
        # Calculate summary statistics for this testid
        counts = {"PASS": 0, "FAIL": 0, "WARN": 0, "TIMEOUT": 0,
                  "CHECK-PASS": 0, "CHECK-FAIL": 0, "CHECK-WARN": 0, "CHECK-TIMEOUT": 0}
        for test_section in all_tests:
            rec = db[testid].get(test_section)
            if rec:
                runStatus = rec.get("runStatus", "?")
                checkStatus = rec.get("checkStatus","?")
                if checkStatus not in {"NONE","?"}:
                    counts[f"CHECK-{checkStatus}"] += 1
                else:
                    counts[runStatus] += 1
        
        html.append("<tr><td>")
        html.append("<div class='testid-cell'>")
        html.append(f"<a href='report/{testid}.html' class='testid-link'>{escape(testid)}</a>")
        html.append("<div class='result-badges'>")
        if counts["PASS"] > 0:
            html.append(f"<span class='result-badge pass'>{counts['PASS']} ✓</span>")
        if counts["CHECK-PASS"] > 0:
            html.append(f"<span class='result-badge check pass'>{counts['CHECK-PASS']} ✓</span>")
        if counts["FAIL"] > 0:
            html.append(f"<span class='result-badge fail'>{counts['FAIL']} ✗</span>")
        if counts["CHECK-FAIL"] > 0:
            html.append(f"<span class='result-badge check fail'>{counts['CHECK-FAIL']} ✗</span>")
        if counts["WARN"] > 0:
            html.append(f"<span class='result-badge warn'>{counts['WARN']} ⚠</span>")
        if counts["CHECK-WARN"] > 0:
            html.append(f"<span class='result-badge check warn'>{counts['CHECK-WARN']} ⚠</span>")
        if counts["TIMEOUT"] > 0:
            html.append(f"<span class='result-badge timeout'>{counts['TIMEOUT']} ⏱</span>")
        if counts["CHECK-TIMEOUT"] > 0:
            html.append(f"<span class='result-badge check timeout'>{counts['CHECK-TIMEOUT']} ⏱</span>")
        html.append("</div>")
        html.append("</div>")
        html.append("</td>")
        
        for test_section in all_tests:
            rec = db[testid].get(test_section)
            if rec:
                runStatus = rec.get("runStatus", "?")
                checkStatus = rec.get("checkStatus","?")
                runtime = rec.get("executionTime", "")
                
                html.append("<td class='cell'><div class='cell-content'>")
                if runStatus == "PASS":
                    if checkStatus != "NONE":
                        html.append(f"<span class='status-badge check status-{checkStatus}'>{checkStatus}</span>")
                    else:
                        html.append(f"<span class='status-badge status-{runStatus}'>{runStatus}</span>")
                        
                    
                
                if runtime:
                    try:
                        time_val = float(runtime)
                        time_str = f"{time_val:.2f}s"
                        # Calculate percentage for progress bar
                        max_time = max_runtimes.get(test_section, 1)
                        pct = min(100, (time_val / max_time) * 100) if max_time > 0 else 0
                        
                        html.append(f"<span class='runtime'>{time_str}</span>")
                        html.append(f"<div class='runtime-bar'><div class='runtime-fill' style='width:{pct}%'></div></div>")
                    except:
                        pass
                
                html.append("</div></td>")
            else:
                html.append("<td class='cell empty-cell'>—</td>")
        html.append("</tr>")
    
    html.append("</tbody></table>")
    html.append("</div>")
    
    # Legend
    html.append("<div class='legend'>")
    html.append("<div class='legend-item'><div class='legend-color' style='background:#d1e7dd;border-color:#badbcc'></div><span>Pass</span></div>")
    html.append("<div class='legend-item'><div class='legend-color' style='background:#f8d7da;border-color:#f5c2c7'></div><span>Fail</span></div>")
    html.append("<div class='legend-item'><div class='legend-color' style='background:#fff3cd;border-color:#ffecb5'></div><span>Warning</span></div>")
    html.append("<div class='legend-item'><div class='legend-color' style='background:#e2e3e5;border-color:#d3d6d8'></div><span>Timeout</span></div>")
    html.append("<div class='legend-item'><span>Runtime bars show relative execution time within each test</span></div>")
    html.append("</div>")
    
    html.append("</div>")
    html.append(html_ftr())
    
    # Write output
    Path(HTML_SUMMARY).write_text("\n".join(html))

def append_html(record):
    if Path(JSON_REPORT).exists():
        db = json.loads(Path(JSON_REPORT).read_text())
    else:
        db = {}
    testname = record['testdir'] + "/" + record['section']
    db[testname] = record
    Path(JSON_REPORT).write_text(json.dumps(db, indent=2))

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
                f.write(f"<a href='../{img}'><img src='../{img}' class='thumbnail'></a>\n")
            f.write("</div>\n")

        f.write("</div>\n")

    write_master_html()

def finalize_html():
    with open(HTML_REPORT, "a") as f:
        f.write("</div>\n")
        f.write(html_ftr())
