from subprocess import Popen, PIPE
import time
import sys
import os

binary_path = 'zxbarren/target/release/benchmark'
quimb_path = 'run_quimb.py'
ansatze = ['sim2', 'sim10', 'sim11', 'sim15']

num_cores = int(sys.argv[1]) if len(sys.argv) > 1 else 1

timeout = 15*60

QUIZX = 'quizx'
OURS = 'ours'
QUIMB = 'quimb'
QUIMB_JAX = 'quimbjax'


def jobs():
    for l in range(1, 4):
        for a in ansatze:
            for n in range(4, 21):
                for m in [OURS, QUIZX]:  # QUIMB, QUIMB_JAX
                    yield {'ansatz': a, 'qubits': n, 'layers': l, 'mode': m, 'parallel': 'false'}


def prefix(job):
    return f"{job['ansatz']}-{job['qubits']}-{job['layers']}-{job['mode']}-{job['parallel']}"


jobs = jobs()
running = []

f = open('result.txt', 'w')

jobs_remaining = True


def kill_all():
    for (_, p, _) in running:
        p.kill()

import atexit
atexit.register(kill_all)


while jobs_remaining or len(running) > 0:
    # Start new jobs
    while jobs_remaining and len(running) < num_cores:
        try:
            job = next(jobs)
        except:
            jobs_remaining = False
            break

        if job['mode'] == QUIZX or job['mode'] == OURS:
            p = Popen([binary_path, job['ansatz'], str(job['qubits']), str(job['layers']), job['mode'], job['parallel']], stdout=PIPE, stderr=PIPE)
            running.append((job, p, time.time()))
        else:
            env = os.environ.copy()
            if job["parallel"] == "false":
                env["QUIMB_NUM_THREAD_WORKERS"] = "1"
                env["QUIMB_NUM_PROCS"] = "1"
                env["OMP_NUM_THREADS"] = "1"
                env["MKL_NUM_THREADS"] = "1"
                env["OPENBLAS_NUM_THREADS"] = "1"
                env["XLA_FLAGS"] = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"
            p = Popen(["python3", quimb_path, job['ansatz'], str(job['qubits']), str(job['layers']), job['mode'], job['parallel']], stdout=PIPE, stderr=PIPE, env=env)
            running.append((job, p, time.time()))

    time.sleep(.01)
    
    # Check if anyone is done
    terminated = []
    for (job, proc, start_time) in running:
        retcode = proc.poll()
        running_time = time.time() - start_time
        if retcode is not None:
            data, err = proc.communicate()
            data = data.decode('ASCII').rstrip()
            terminated.append((job, proc, start_time))
            if retcode == 0:
                line = f"{prefix(job)}: {data}"
            else:
                line = f"{prefix(job)}: failed"
            f.write(line + "\n")
            print(line)
        elif running_time > timeout:
            proc.kill()
            terminated.append((job, proc, start_time))
            line = f"{prefix(job)}: timeout"
            f.write(line + "\n")
            print(line)

    for t in terminated:
        running.remove(t)

