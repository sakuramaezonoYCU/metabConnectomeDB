import psutil

print("Running Python processes:")
for p in psutil.process_iter(['pid', 'name', 'cmdline']):
    try:
        if 'python' in p.info['name'].lower() or (p.info['cmdline'] and any('python' in arg.lower() for arg in p.info['cmdline'])):
            cmd = " ".join(p.info['cmdline']) if p.info['cmdline'] else p.info['name']
            if 'run_all_cancers' in cmd or 'run_cancer_pipeline' in cmd:
                print(f"PID: {p.info['pid']} - CMD: {cmd}")
    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
        pass
