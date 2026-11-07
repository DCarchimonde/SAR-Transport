"""Observe existing authorized jobs; never launch a fit or alter its configuration."""
import datetime as dt,json,os,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/sar_transport_v1'
RUNS=ROOT/'runs/sar_transport_v1'

def snapshot():
    completed=[]
    for p in RUNS.rglob('COMPLETE.json'):
        try:r=json.loads(p.read_text())
        except json.JSONDecodeError:continue
        completed.append({'path':str(p.relative_to(ROOT)),'cpu_seconds':r.get('cpu_seconds',0)})
    completed_paths={x['path'] for x in completed}
    active=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            args=(p/'cmdline').read_bytes().split(b'\0')
            if not args[0].endswith(b'python') or not any(a.endswith(b'/run_sar_transport.py') for a in args):continue
            stat=(p/'stat').read_text().split()
            rss=int(stat[23])*os.sysconf('SC_PAGE_SIZE')
            decoded=[a.decode() for a in args if a]
            seed=decoded[decoded.index('--seed')+1] if '--seed' in decoded else None
            model=decoded[decoded.index('--model')+1] if '--model' in decoded else 'rf'
            already_counted=bool(seed and str((RUNS/f'{model}_seed{seed}'/'COMPLETE.json').relative_to(ROOT)) in completed_paths)
            offset_file=RUNS/f'{model}_seed{seed}'/'CPU_CONTINUATION.json';offset=0.
            if offset_file.exists():
                try:offset=json.loads(offset_file.read_text()).get('prior_saved_cpu_seconds',0.)
                except json.JSONDecodeError:pass
            active.append({'pid':int(p.name),'args':decoded,'completion_already_counted':already_counted,
                'cpu_seconds':(int(stat[13])+int(stat[14]))/os.sysconf('SC_CLK_TCK'),'prior_saved_cpu_seconds':offset,'rss_bytes':rss})
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    # Pipeline preparation is explicitly outside individual fit timers.
    prep=0
    for name in ['ELIGIBILITY_SUMMARY.json','FEATURES.json']:
        r=json.loads((OUT/name).read_text())
        prep+=r.get('cpu_seconds',0)
    recovery_file=OUT/'RECOVERY_CPU_CHARGES.json';recovery=json.loads(recovery_file.read_text()) if recovery_file.exists() else {}
    known=recovery.get('known_interrupted_cpu_seconds',0.);conservative=recovery.get('conservative_unobserved_gap_charge_seconds',0.)
    cpu=sum(x['cpu_seconds'] for x in completed)+sum(x['cpu_seconds']+x['prior_saved_cpu_seconds'] for x in active if not x['completion_already_counted'])+prep+known
    return {'utc':dt.datetime.now(dt.timezone.utc).isoformat(),'completed_cpu_seconds':sum(x['cpu_seconds'] for x in completed),
       'completed_records':len(completed),'active':active,'preparation_cpu_seconds':prep,
       'observed_total_cpu_seconds':cpu,'conservative_unobserved_gap_charge_seconds':conservative,'budget_charged_cpu_seconds':cpu+conservative,
       'known_interrupted_cpu_seconds':known,'aggregate_active_rss_bytes':sum(x['rss_bytes'] for x in active),
       'budget_cpu_seconds':48*3600,'old_discarded_engineering_fit_included':True}

def main():
    while True:
        r=snapshot()
        with open(OUT/'RESOURCE_MONITOR.jsonl','a') as f:f.write(json.dumps(r)+'\n')
        stop=r['budget_charged_cpu_seconds']>48*3600-600 or r['aggregate_active_rss_bytes']>15.5*1024**3
        if stop:
            (OUT/'RESOURCE_STOP.json').write_text(json.dumps(r,indent=2))
            # /proc exposes outer PIDs across command namespaces. They cannot
            # safely be signalled here. The supervising execution session must
            # interrupt its existing jobs; never claim that this monitor killed them.
            print('RESOURCE_BOUNDARY: supervising session must stop active jobs.',flush=True)
            return
        time.sleep(10)

if __name__=='__main__':main()
