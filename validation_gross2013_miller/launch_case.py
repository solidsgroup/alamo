#!/usr/bin/env python3
"""Run one prepared case with an immutable binary and automatic provenance.

Simulation-launch agents invoke this script. It performs no scientific
postprocessing and never changes an input deck or reruns an existing receipt.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path, data):
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n');temporary.replace(path)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('--ranks',type=int,default=1)
    args=parser.parse_args()
    if (HERE/'PAUSED').exists():
        parser.error('Study is paused by the user; no simulations may launch until the user resumes it')
    if args.ranks<1:parser.error('--ranks must be positive')
    case=args.case.resolve()
    if case.parent!=HERE/'runs':parser.error('case must be a direct child of this study/runs')
    deck=case/'input';metadata=json.loads((case/'case.json').read_text())
    if digest(deck)!=metadata['input_sha256']:raise RuntimeError('Input differs from prepared case metadata')
    receipt_path=case/'run.json'
    if receipt_path.exists() or (case/'stdout.log').exists():raise RuntimeError('Case already launched; prepare a new case for another attempt')
    # Content-addressing prevents subsequent source rebuilds from invalidating
    # the launch hash. Other launcher instances may safely share this file.
    source=ROOT/'bin/lowmach-2d-clang++';binary_hash=digest(source)
    binary=HERE/'binaries'/f'lowmach-{binary_hash}'
    binary.parent.mkdir(exist_ok=True)
    if not binary.exists():
        temporary=binary.with_suffix(f'.{os.getpid()}.tmp')
        shutil.copy2(source,temporary)
        if digest(temporary)!=binary_hash:raise RuntimeError('Executable changed while snapshotting; wait for the build')
        temporary.chmod(0o555)
        temporary.replace(binary)
    if digest(binary)!=binary_hash:raise RuntimeError('Archived executable is corrupt')
    command=[str(binary),str(deck)]
    if args.ranks>1:command=['mpirun','-np',str(args.ranks)]+command
    receipt=dict(status='running',command=command,input_sha256=digest(deck),binary_sha256=binary_hash,
        binary_snapshot=str(binary.relative_to(ROOT)),hashes_captured_before_launch=True,
        started_utc=datetime.now(timezone.utc).isoformat(),returncode=None,mpi_ranks=args.ranks)
    write_json(receipt_path,receipt);started=time.monotonic()
    print(json.dumps(receipt),flush=True)
    env=os.environ.copy();env.setdefault('OMP_NUM_THREADS','1')
    try:
        with (case/'stdout.log').open('w') as log:
            result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=env,check=False)
        receipt.update(returncode=result.returncode,status='completed' if result.returncode==0 else 'failed',
            wall_seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat(),
            amrex_finalized='AMReX (26.06) finalized' in (case/'stdout.log').read_text(),
            input_unchanged=digest(deck)==receipt['input_sha256'])
    except BaseException as exc:
        receipt.update(status='launcher_error',error=repr(exc),wall_seconds=time.monotonic()-started,
                       finished_utc=datetime.now(timezone.utc).isoformat())
        write_json(receipt_path,receipt);raise
    write_json(receipt_path,receipt);print(json.dumps(receipt),flush=True)
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
