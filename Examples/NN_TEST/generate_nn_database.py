import argparse
import re
import subprocess
from pathlib import Path


def temperature_count(filename):
    text = Path(filename).read_text()
    text = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    match = re.search(r"\bN_RTMP\s*=\s*(\d+)", text, re.IGNORECASE)
    return int(match.group(1)) if match else 44


def main():
    p = argparse.ArgumentParser()
    p.add_argument("executable")
    p.add_argument("-i", "--input", default="RADCAL_NN.IN")
    p.add_argument("-o", "--output-prefix", default="radcal_nn")
    p.add_argument("-n", "--processes", type=int, default=1)
    p.add_argument("--mpiexec", default="mpiexec")
    p.add_argument("--oversubscribe", action="store_true")
    a = p.parse_args()
    if a.processes < 1:
        p.error("--processes must be at least 1")
    n_temperature = temperature_count(a.input)
    if a.processes > n_temperature:
        print(f"Requested {a.processes} processes; using N_RTMP={n_temperature}.", flush=True)
        a.processes = n_temperature

    command = [a.executable]
    if a.processes > 1:
        command = [a.mpiexec]
        if a.oversubscribe:
            command += ["--map-by", ":OVERSUBSCRIBE"]
        command += ["-n", str(a.processes), a.executable]
    subprocess.run(
        command,
        input=f"{a.input}\n{a.output_prefix}\n",
        text=True,
        check=True
    )


if __name__ == "__main__":
    main()
