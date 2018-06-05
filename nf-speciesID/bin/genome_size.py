#!/usr/bin/env python3

from subprocess import Popen, PIPE, CalledProcessError
import shlex
import os
from argparse import ArgumentParser

def run_command(cmd):
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    out = stdout.decode("utf-8")
    err = stderr.decode("utf-8")
    if p.returncode != 0:
        print('STDERR from called program: {}'.format(stderr))
        print('STDOUT from called program: {}'.format(stdout))
        raise CalledProcessError(p.returncode, cmd)
    else:
        return [out, err]


def genome_size_bbtools(sample_id, forward, reverse):
    cmd = "kmercountexact.sh in={} in2={} peaks={}.peak".format(forward, reverse, sample_id)
    cmd_to_run = shlex.split(cmd)

    p = Popen(cmd_to_run, stderr=PIPE, stdout=PIPE)
    _, stderr = p.communicate()
    gs = cov = "undefined"
    if stderr != 0:
        with open("{}.peak".format(sample_id), "r") as fh:
            handle = fh.readlines()
            for line in handle:
                if "#haploid_genome_size" in line:
                    gs = line.split("\t")[1].strip()
                if "#haploid_fold_coverage" in line:
                    cov = line.split("\t")[1].strip()
                    break
        # with open("{}.tsv".format(sample_id), "w") as oh:
        #     oh.writelines("{},{}".format(gs, cov))
    return [gs, cov]

def genome_size(fastq, mincov=25, minkc=3):
    cmd = "mash sketch -o tmp -k 32 -p 4 -m {} -c {} -r {}".format(minkc, mincov, fastq)
    cmd_to_run = shlex.split(cmd)
    _, stderr = run_command(cmd_to_run)
    lines = stderr.split("\n")
    genome_size = lines[1].replace("Estimated genome size: ", "")
    gs = int(float(genome_size))
    os.unlink("tmp.msh")
    return gs


def estimate_coverage(fastq, genome_size):
    cmd = "seqkit stat {}".format(fastq)
    cmd_to_run = shlex.split(cmd)
    stdout, _ = run_command(cmd_to_run)
    lines = stdout.split("\n")
    bp = lines[1].split()[4].replace(",", "")
    est_cov = int(bp) * 2 / genome_size
    return est_cov

def find_refseq(forward, reverse, refseq, temp = "tmp"):
    if forward.endswith("gz"):
        extr = "zcat"
    else:
        extr = "cat"
    # Sketching read
    sketch_cmd = "{} \"{}\" \"{}\" | mash sketch -o {} -m 2 -r -".format(extr, forward, reverse, temp)
    sketch_cmd_to_run = shlex.split(sketch_cmd)
    run_command(sketch_cmd_to_run)
    # Estimate distance between refseq and read
    dist_cmd = "mash dist {} {}.msh | sort -gk3 | head -n1".format(refseq, temp)
    dist_cmd_to_run = shlex.split(dist_cmd)
    out, _ = run_command(dist_cmd_to_run)
    lines = out.split("\n")
    refseq = lines.split()[0]
    return refseq

def parse_args():
    parser = ArgumentParser(description="Estimate genome size and closely related refseq")
    parser.add_argument('--R', help="Read (Forward or Reverse or SE)", required=True)
#    parser.add_argument('--R2', help="Reverse read", required=False)
#    parser.add_argument('--refseq', help="RefSeq Sketch", required=False, default="/home/ubuntu/data2/mash/refseq.k21s1000.042018.msh")
    return parser.parse_args()

def main():
    args = parse_args()
    forward = args.R
#    reverse = args.R2
#    refseq = args.refseq

    gs = genome_size(forward)
    cv = estimate_coverage(forward, gs)
#    ref_seq = find_refseq(forward, reverse, refseq)
    print("{},{:.1f}".format(gs, cv))

if __name__ == "__main__":
    main()

