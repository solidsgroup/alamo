#!/usr/bin/env python3
"""Launch-only serial worker; root's independent monitor makes rate decisions."""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def stamp():return datetime.now(timezone.utc).isoformat()


def write_json(path,data):
    temporary=path.with_suffix('.json.tmp');temporary.write_text(json.dumps(data,indent=2)+'\n');temporary.replace(path)


def monitor_available(max_age):
    try:
        d=json.loads((HERE/'analysis/transient_monitor_state.json').read_text())
        age=(datetime.now(timezone.utc)-datetime.fromisoformat(d['heartbeat_utc'])).total_seconds()
        return d['status']=='running' and 0<=age<=max_age
    except (OSError,ValueError,KeyError):return False


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('manifest',type=Path)
    p.add_argument('--check',action='store_true');args=p.parse_args()
    manifest=args.manifest.resolve()
    if manifest.parent!=HERE/'queues':p.error('Manifest must belong to this study/queues')
    config=json.loads(manifest.read_text());ranks=config['mpi_ranks']
    if not isinstance(ranks,int) or ranks<1:p.error('Invalid MPI rank count')
    cases=[ROOT/row['case'] for row in config['cases']]
    if not cases or len(cases)!=len(set(cases)):p.error('Cases must be distinct')
    def check(case,row):
        if case.resolve().parent!=HERE/'runs':raise RuntimeError('Case outside study')
        if digest(ROOT/'bin/lowmach-2d-clang++')!=config['binary_sha256']:raise RuntimeError('Executable changed')
        m=json.loads((case/'case.json').read_text())
        if digest(case/'input')!=row['input_sha256'] or m['input_sha256']!=row['input_sha256']:raise RuntimeError('Input changed')
        if not m.get('transient_steady_stop'):raise RuntimeError('Wrong integration/stop mode')
        if any((case/name).exists() for name in ('run.json','stdout.log','STOP')):raise RuntimeError('Case already launched or stopped')
    for case,row in zip(cases,config['cases']):check(case,row)
    if args.check:print(f'Validated {len(cases)} transient cases');return
    state_path=manifest.with_name(manifest.stem+'.state.json')
    with manifest.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if state_path.exists():raise RuntimeError('Queue already started')
        state=dict(status='running',pid=os.getpid(),started_utc=stamp(),binary_sha256=config['binary_sha256'],
                   manifest_sha256=digest(manifest),mpi_ranks=ranks,active_case=None,completed=[],
                   pending=[str(c.relative_to(ROOT)) for c in cases])
        write_json(state_path,state)
        try:
            for case,row in zip(cases,config['cases']):
                check(case,row)
                if (HERE/'PAUSED').exists():raise RuntimeError('Study paused')
                if not monitor_available(config['monitor_heartbeat_max_age_s']):raise RuntimeError('Statistical monitor unavailable; no launch')
                if shutil.disk_usage(HERE).free<8*1024**3:raise RuntimeError('Less than 8 GiB free disk')
                state['active_case']=str(case.relative_to(ROOT));write_json(state_path,state)
                child=subprocess.Popen([sys.executable,str(HERE/'launch_case.py'),str(case),'--ranks',str(ranks)],cwd=ROOT)
                lost_monitor=False
                while child.poll() is None:
                    if not monitor_available(config['monitor_heartbeat_max_age_s']):
                        (case/'STOP').write_text('Statistical monitor unavailable; stop without accepting a rate.\n')
                        lost_monitor=True
                    time.sleep(2)
                receipt=json.loads((case/'run.json').read_text())
                if receipt['status']=='running':raise RuntimeError('Child exited without final receipt')
                state['completed'].append(dict(case=str(case.relative_to(ROOT)),status=receipt['status'],returncode=receipt.get('returncode')))
                state['pending'].remove(str(case.relative_to(ROOT)));state['active_case']=None;write_json(state_path,state)
                if lost_monitor:raise RuntimeError('Stopped after loss of statistical monitor')
            state['status']='completed' if all(r['returncode']==0 for r in state['completed']) else 'completed_with_failures'
        except BaseException as exc:
            state.update(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'stopped',error=repr(exc));raise
        finally:
            state['finished_utc']=stamp();write_json(state_path,state);print(json.dumps(state),flush=True)


if __name__=='__main__':main()
