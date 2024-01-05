SignalP Peptide Prediction Integration
======================================

Our program integrates peptide prediction functionality using SignalP. To enable this feature, please follow these steps:

1. Activate SignalP in your program by setting the parameter ``use_signalP=True``.
2. Acquire an academic license for SignalP and download it `from the official site <https://services.healthtech.dtu.dk/services/SignalP-4.1/>`_.
3. Extract the perl file from the downloaded tarball (signalp-4.1g.Linux.tar.gz) and move it to ``/usr/bin/signalp``.

To run the program with SignalP, use the following command:

.. code-block:: bash

    run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655 --use_signalP=TRUE

.. warning::
    If you lack permission to access `/usr/bin`, specify the path of the SignalP executable file using the `-sp` or `--signalP_path` parameter. Here's an example command:

    .. code-block:: bash

        run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655 --use_signalP=TRUE -sp /home/lehuang/Downloads/signalp-4.1/signalp

SignalP-4.1 Installation Instructions
-------------------------------------

Begin by decompressing the SignalP tarball and navigating to its directory:

.. code-block:: bash

   tar -xvf signalp-4.1g.Linux.tar.gz && cd signalp-4.1

Inside the `signalp-4.1` directory, you'll find the following files and directories:

.. code-block:: bash

   (base) lehuang@lehuang:~/Downloads/signalp-4.1$ ls
   bin  lib  signalp  signalp.1  signalp-4.1.readme  syn  test

The `signalp` file is the perl script that will be utilized in your program.

Customizing the SignalP Script
------------------------------

Modify the "GENERAL SETTINGS, CUSTOMIZE ..." section at the start of the `signalp` file. Ensure these mandatory variables are correctly set:

- **SIGNALP**: Specify the full path to the signalp-4.1 directory on your system.
- **outputDir**: Choose a directory for storing temporary files (must be writable by all users).
- **MAX_ALLOWED_ENTRIES**: Define the maximum number of input sequences allowed per run.

Here's an example of how to configure these settings in the `signalp` file:

.. code-block:: bash

   ##############################################################################
   #               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
   ##############################################################################

   # Full path to the signalp-4.1 directory (mandatory)
   BEGIN {
       $ENV{SIGNALP} = '/home/lehuang/Downloads/signalp-4.1';
   }

   # Directory for temporary files (writable by all users)
   my $outputDir = "/home/lehuang/Downloads/signalp-4.1/output";

   # Max number of sequences per run (flexible)
   my $MAX_ALLOWED_ENTRIES=100000;

Copying the SignalP Script to /usr/bin (if accessible)
------------------------------------------------------

If you have the necessary permissions, use these commands to copy the `signalp` script:

.. code-block:: bash

   sudo cp signalp /usr/bin/signalp
   sudo chmod 755 /usr/bin/signalp
