import argparse
import subprocess


def main():
    p = argparse.ArgumentParser()
    p.add_argument("executable")
    p.add_argument("-i", "--input", default="RADCAL_NN.IN")
    p.add_argument("-o", "--output-prefix", default="radcal_nn")
    a = p.parse_args()

    subprocess.run(
        [a.executable],
        input=f"{a.input}\n{a.output_prefix}\n",
        text=True,
        check=True
    )


if __name__ == "__main__":
    main()
