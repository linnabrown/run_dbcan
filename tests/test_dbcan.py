# run_dbcan testing file

import os
import shlex
from pathlib import Path

from dbcan.cli.run_dbcan import run_dbCAN, rundbCAN_parser

TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")


def run_dbcan(args):
    # wrapper for calling run_dbcan from parsed args

    run_dbCAN(
        inputFile=args.inputFile,
        inputType=args.inputType,
        cluster=args.cluster,
        dbCANFile=args.dbCANFile,
        dia_eval=args.dia_eval,
        dia_cpu=args.dia_cpu,
        hmm_eval=args.hmm_eval,
        hmm_cov=args.hmm_cov,
        hmm_cpu=args.hmm_cpu,
        dbcan_thread=args.dbcan_thread,
        tf_eval=args.tf_eval,
        tf_cov=args.tf_cov,
        tf_cpu=args.tf_cpu,
        stp_eval=args.stp_eval,
        stp_cov=args.stp_cov,
        stp_cpu=args.stp_cpu,
        prefix=args.out_pre,
        outDir=args.out_dir,
        dbDir=args.db_dir,
        cgc_dis=args.cgc_dis,
        cgc_sig_genes=args.cgc_sig_genes,
        tool_arg=args.tools,
        use_signalP=args.use_signalP,
        signalP_path=args.signalP_path,
        gram=args.gram,
        verbose=args.verbose,
    )


class Test_dbCAN:
    """Test Class for run_dbcan"""

    def test_dbcan_with_isolated_genome(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.fna prok"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)

        args.out_dir = "tests/test_isolated_genome_as_input"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_from_protein(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.faa protein"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.out_dir = "tests/test_protein_as_input"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_with_meta_genome(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.fna meta"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.out_dir = "tests/test_meta_genome_as_input"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_cazyme_and_cgc_with_meta_genome(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.cluster = os.path.join(DATA_ROOT, args.cluster)
        args.out_dir = "tests/test_cazyme_and_cgc_meta_genome_as_input"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_cgc_substrate_prediction(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff --cgc_substrate --only_sub"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.cluster = os.path.join(DATA_ROOT, args.cluster)
        args.out_dir = "tests/test_cgc_substrate_prediction"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_cazyme_cgc_substrate_prediction(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff --cgc_substrate"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.cluster = os.path.join(DATA_ROOT, args.cluster)
        args.out_dir = "tests/test_cazyme_cgc_substrate_prediction"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_cazyme_and_cgc_substrate_with_protein_and_gff(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --cgc_substrate"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.cluster = os.path.join(DATA_ROOT, args.cluster)
        args.out_dir = "tests/test_cazyme_cgc_substrate_with_protein_and_gff_as_input"

        # run dbcan pipeline
        run_dbcan(args)

    def test_dbcan_cazyme_and_cgc_substrate_with_protein_and_gff_hmmer(self):
        parser = rundbCAN_parser()

        # parse the command from string directly
        cmd = "run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --cgc_substrate --tools hmmer"
        arguments = shlex.split(cmd)[1:]

        args = parser.parse_args(arguments)
        args.inputFile = os.path.join(DATA_ROOT, args.inputFile)
        args.cluster = os.path.join(DATA_ROOT, args.cluster)
        args.out_dir = "tests/test_cazyme_cgc_substrate_with_protein_and_gff_as_input_using_hmmer"

        # run dbcan pipeline
        run_dbcan(args)
