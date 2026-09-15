#!/usr/bin/env python3
"""Keep local pilot diagnostics current; never launch or approve production runs."""
import fcntl
import json
import os
import re
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from analyze_runs import analyze
from plot_pilots import main as plot_pilots

HERE=Path(__file__).resolve().parent
CASES=[HERE/'runs'/f'M03_p68.0000_seed101_pilot_{s}_v5_restart'
       for s in ('dx1_w4','dx0.5_w8','dx0.5_w4')]
CASES.append(HERE/'runs/M24_p67.3863_seed101_pilot_ignition')

def log_tail(path):
    if not path.exists():return ''
    with path.open('rb') as f:
        f.seek(max(0,path.stat().st_size-8192))
        return f.read().decode(errors='replace')

def main():
    lock=(HERE/'analysis/pilot_watch.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    lock.write(str(os.getpid()));lock.flush()
    signatures={}
    while True:
        rows=[];changed=False
        for case in CASES:
            receipt=case/'run.json'
            headers=sorted((case/'output').glob('*cell/Header'))
            signature=tuple((str(p),p.stat().st_mtime_ns) for p in headers)
            signature+=(('receipt',receipt.stat().st_mtime_ns if receipt.exists() else 0),)
            error=None
            if signature!=signatures.get(case):
                try:
                    analyze(case);signatures[case]=signature;changed=True
                except Exception:
                    error=traceback.format_exc()
                    print(error,flush=True)
            tail=log_tail(case/'stdout.log')
            steps=re.findall(r'STEP (\d+) ends\. TIME = ([\deE+.-]+) DT = ([\deE+.-]+)',tail)
            r=json.loads(receipt.read_text()) if receipt.exists() else {}
            summary=case/'analysis/summary.json'
            s=json.loads(summary.read_text()) if summary.exists() else {}
            rc=r.get('returncode',r.get('exit_code'))
            finalized='AMReX (26.06) finalized' in tail
            aborted='SIGABRT' in tail or 'MPI_ABORT' in tail
            state=('completed' if rc==0 else 'failed') if rc is not None else (
                'finalized; receipt pending' if finalized else 'aborted; receipt pending' if aborted else 'running')
            rows.append(dict(case=case.name,state=state,returncode=rc,
                step=int(steps[-1][0]) if steps else None,
                last_log_time_s=float(steps[-1][1]) if steps else s.get('last_time_s'),
                last_analyzed_time_s=s.get('last_time_s'),
                analysis_status=s.get('status'),analysis_error=error,
                terminal=rc is not None or finalized or aborted))
        if changed:plot_pilots()
        now=datetime.now(timezone.utc).isoformat()
        status=dict(updated_utc=now,watcher_pid=os.getpid(),cases=rows,
            production_launched=False,accepted_figure10_cases=0,
            note='Pilot diagnostics only. Root-model review of sustained burning, resolution, and packed recession is required before the production sweep.')
        p=HERE/'analysis/pilot_watch_status.json';tmp=p.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(status,indent=2)+'\n');tmp.replace(p)
        lines=['# Local pilot status','',f'Updated: {now}','',
               'Four pilots are running or being checked before the 84-case production sweep. No Figure 10 simulation rates have been accepted.','',
               '| Case | Process state | Simulated time (ms) | Analysis status |',
               '|---|---|---:|---|']
        for row in rows:
            tm=row['last_log_time_s']
            lines.append(f"| {row['case']} | {row['state']} | {tm*1e3:.4f} | {row['analysis_status']} |" if tm is not None else f"| {row['case']} | {row['state']} | — | {row['analysis_status']} |")
        lines+=['','[Ignition and resolution histories](analysis/ignition_resolution_pilots.png)',
                '','[Study assumptions and numerical checks](README.md)',
                '','The monitor refreshes root-authored postprocessing when new snapshots appear. It never launches production cases or treats a transient pilot slope as a validated burning rate.','']
        tmp=HERE/'STATUS.md.tmp';tmp.write_text('\n'.join(lines));tmp.replace(HERE/'STATUS.md')
        if all(row['terminal'] for row in rows):
            print('All pilots have terminal log or receipt status; root review remains required.',flush=True)
            break
        time.sleep(60)

if __name__=='__main__':main()
