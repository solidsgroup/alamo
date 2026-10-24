#!/usr/bin/env python3
"""Root-authored live postprocessing and statistical stopping for transient runs."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import time
from assess_steady import assessment

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
STATE=HERE/'analysis/transient_monitor_state.json'


def now():return datetime.now(timezone.utc)


def write_json(path,data):
    path.parent.mkdir(exist_ok=True)
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n');temporary.replace(path)


def stop_case(case,reason):
    path=case/'STOP'
    if not path.exists():path.write_text(reason+'\n')


def completed_signature(case):
    visit=case/'output/celloutput.visit'
    if not visit.exists():return ()
    return tuple(line.strip() for line in visit.read_text().splitlines(keepends=True)
                 if line.endswith('\n') and line.strip().endswith('/Header'))


def record_acceptance(case,result,receipt,detected=None):
    """Record detection times before signaling the existing clean-stop hook."""
    if not result.get('accepted') or result.get('status')!='accepted_steady':return False
    path=case/'analysis/steady_stop.json'
    if path.exists():return False
    detected=detected or now();candidate=result['candidate']
    meta=json.loads((case/'case.json').read_text())
    restart_time=meta.get('continuation_start_s',0.0)
    write_json(path,dict(reason='confirmed_statistical_stationarity',
        detected_utc=detected.isoformat(),
        time_to_steady_simulation_s=result['last_time_s'],
        first_passing_simulation_s=result['first_passing_time_s'],
        restart_initial_time_s=restart_time,
        time_to_steady_since_restart_s=result['last_time_s']-restart_time,
        time_to_steady_wall_s=(detected-datetime.fromisoformat(receipt['started_utc'])).total_seconds(),
        rate_cm_s=candidate['block_mean_cm_s'],
        temporal_ci95_halfwidth_cm_s=candidate['temporal_ci95_halfwidth_cm_s'],
        acceptance=result,stop_requested_while_running=receipt['status']=='running',
        note='Detection is limited by output cadence and monitor polling; wall time includes this case initialization and restart loading.'))
    if receipt['status']=='running':stop_case(case,'Confirmed statistically steady regression; see analysis/steady_stop.json')
    return True


def watch_case(case,receipt):
    result=assessment(case,completed_only=True)
    record_acceptance(case,result,receipt)
    if result['status']=='cooling_or_extinguished' and receipt['status']=='running':
        stop_case(case,'Flame cooled/extinguished; no steady burning rate accepted.')
    record=case/'analysis/steady_stop.json'
    if receipt['status']!='running':
        final=dict(execution_status=receipt['status'],returncode=receipt.get('returncode'),
                   solver_finalized=receipt.get('amrex_finalized',False),
                   final_simulation_time_s=result.get('last_time_s'),
                   actual_wall_s=receipt.get('wall_seconds'),
                   accepted=bool(result.get('accepted') and receipt.get('returncode')==0 and receipt.get('amrex_finalized')),
                   assessment_status=result['status'])
        if record.exists():
            data=json.loads(record.read_text());data['final_run']=final;write_json(record,data)
        else:
            final['reason']='run_ended_without_confirmed_stationarity'
            write_json(case/'analysis/termination.json',final)
    return dict(execution_status=receipt['status'],assessment_status=result['status'],
                last_time_s=result.get('last_time_s'),accepted=bool(result.get('accepted')),
                stop_requested=(case/'STOP').exists())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--interval',type=float,default=15)
    args=p.parse_args()
    if args.interval<=0:p.error('interval must be positive')
    cases=[ROOT/path for path in (HERE/'production_cases.txt').read_text().splitlines() if path]
    for case in cases:
        if not json.loads((case/'case.json').read_text()).get('transient_steady_stop'):
            p.error('All listed cases must use transient statistical stopping')
    STATE.parent.mkdir(exist_ok=True)
    with STATE.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state=dict(status='running',pid=os.getpid(),started_utc=now().isoformat(),cases={})
        seen={};errors={}
        def heartbeat():
            state['heartbeat_utc']=now().isoformat();write_json(STATE,state)
        try:
            while True:
                heartbeat()
                if (HERE/'PAUSED').exists():
                    for case in cases:
                        rp=case/'run.json'
                        if rp.exists() and json.loads(rp.read_text())['status']=='running':
                            stop_case(case,'Study paused by user')
                    state['status']='paused';break
                changed=False
                for case in cases:
                    rp=case/'run.json'
                    if not rp.exists():continue
                    receipt=json.loads(rp.read_text());sig=(receipt['status'],completed_signature(case))
                    if seen.get(case)==sig:continue
                    try:
                        state['cases'][case.name]=watch_case(case,receipt)
                        seen[case]=sig;errors[case]=0;changed=True
                    except Exception as exc:
                        errors[case]=errors.get(case,0)+1
                        state['cases'][case.name]=dict(status='postprocessing_error',error=repr(exc),attempts=errors[case])
                        if errors[case]>=3:
                            stop_case(case,'Postprocessing repeatedly failed; no rate accepted. See monitor state.')
                            seen[case]=sig
                    heartbeat()
                if changed:
                    from summarize_sweep import main as summarize
                    summarize()
                queues=[HERE/'queues'/f'fig10_transient_{label}.state.json' for label in ('small','large')]
                if all(q.exists() and json.loads(q.read_text())['status']!='running' for q in queues):
                    pending_analysis=False
                    for case in cases:
                        rp=case/'run.json'
                        if rp.exists():
                            sig=(json.loads(rp.read_text())['status'],completed_signature(case))
                            pending_analysis |= seen.get(case)!=sig
                    if not pending_analysis:
                        state['status']='queues_finished';break
                time.sleep(args.interval)
        except BaseException as exc:
            state.update(status='stopped',error=repr(exc));raise
        finally:
            state['finished_utc']=now().isoformat();heartbeat()
            print(json.dumps(state),flush=True)


if __name__=='__main__':main()
