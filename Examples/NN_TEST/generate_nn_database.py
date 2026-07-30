import argparse
import subprocess


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
    if a.processes > 44:
        print(f"Requested {a.processes} processes; using the maximum of 44.", flush=True)
        a.processes = 44

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
