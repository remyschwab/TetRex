#! /usr/bin/python3

import argparse


def parse_truth(truth_handle):
    ids = set()
    lines = [line.split(",")[0] for line in truth_handle.readlines()]
    for line in lines:
        tmp_id = line.split("_")[1]
        id = int(tmp_id.split(".")[0])
        ids.add(id)
    return ids


def parse_output(t):
    return set([i for i, x in enumerate(t) if x])


def evaluate(truth, hits):
    tp_count = len(truth)
    fps = len(hits-truth)
    fns = len(truth-hits)
    output_string = "FP COUNT: {} | FN COUNT: {}".format(fps, fns)
    print(output_string)


def main(truth_filename, output):
    handle = open(truth_filename)
    true_bin_ids = parse_truth(handle)
    hit_ids = parse_output(output)
    evaluate(true_bin_ids, hit_ids)


if __name__ == "__main__":
    CLI = argparse.ArgumentParser()
    CLI.add_argument(
        "-p", "--path", required=True, help="The path to your bins.txt file"
    )
    CLI.add_argument(
        "-i", "--input",  # name on the CLI - drop the `--` for positional/required parameters
        nargs="*",  # 0 or more values expected => creates a list
        type=int
    )
    ARGS = CLI.parse_args()
    TRUTH = ARGS.path
    OUTPUT = ARGS.input
    main(TRUTH, OUTPUT)
