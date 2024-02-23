Installation
============

We support different ways to install `dbcan`_, including:

- Using `Anaconda`_ or `Miniconda`_ (**Recommended**)
- Using `PyPI`_
- Using `Docker`_

.. note::

   If you prefer not installing `dbcan`_ locally, you can also use it via our online `server <https://bcb.unl.edu/dbCAN2/index.php>`_.

Requirements
------------

- A Posix-compliant operating system, e.g. ``Linux`` or ``MacOS``.
- A ``Python`` 3.6 or later environment (you can use ``conda`` to create it).
- When using ``Conda`` or ``PyPI`` to install `dbcan`_, you also need to prepare the ``databases`` used by `dbcan`_ seperately (See :doc:`user_guide/database_preparation`).


Installing with Conda
---------------------

If you haven't already installed ``conda``, you need to install a ``conda`` environment. ``Conda`` is available through the `Anaconda <https://docs.anaconda.com/free/anaconda/>`_
or `Miniconda <https://docs.conda.io/projects/miniconda/en/latest/>`_. Then, you can create a new ``conda`` environment (optional but recommended) using the command:

.. code-block:: shell

    conda create --name dbcan python=3.8

If you already have a ``conda`` environment, you can skip the step above.

To install the `dbcan`_ package, use the ``conda install`` command:

.. code-block:: shell

    conda install dbcan -c conda-forge -c bioconda



Build database
--------------

You can build database via this command:

.. code-block:: shell

    dbcan_build --cpus 8 --db-dir db --clean


1. `--cpu` indicates count of cpu you can use. Try as many as possible for fast building.

2. `--db-dir` indicates database folder path. Default is `db` on your current database

3. `--clean` indicates clean the folder indicated by `--db-dir`. 
You can remove this parameter if you don't want to clean, but we recommend you add this to keep
away from index contamination.


Installing SignalP (Optional)
--------------------------------


- `SignalP` is optional and not essential for the core functionality of our software. Users requiring its specific features can integrate it as follows:
   1. Visit the `SignalP website <https://services.healthtech.dtu.dk/services/SignalP-4.1/>`_.
   2. Submit a download `request <https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=4.1&packageversion=4.1g&platform=Linux>`_.
   3. Post-download, add `SignalP` to your system's environmental variables to make it executable.

- For installation assistance, refer to the :doc:`faq/signalp_installation`.



.. Installing with PyPI
.. --------------------

.. To install the `dbcan`_ package via ``pip``, you first need to install a few executable
.. dependencies:

.. - `NCBI-BLAST+ <https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html>`_;
.. - `HMMER <http://hmmer.org/>`_ (:cite:`2011:hmmer`);
.. - `DIAMOND <https://github.com/bbuchfink/diamond>`_ (:cite:`2021:diamond`);
.. - `SignalP <https://services.healthtech.dtu.dk/services/SignalP-4.1/>`_ (:cite:`2017:nielsen`) (Optional).

.. .. warning::

..    **SignalP Integration Notice**

..    Due to the specific licensing terms of `SignalP`, it is not included directly as a dependency in our package. This requires users to undertake a separate installation process.

..    **Installing SignalP (Optional)**:

..       - `SignalP` is optional and not essential for the core functionality of our software. Users requiring its specific features can integrate it as follows:
..          1. Visit the `SignalP website <https://services.healthtech.dtu.dk/services/SignalP-4.1/>`_.
..          2. Submit a download `request <https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=4.1&packageversion=4.1g&platform=Linux>`_.
..          3. Post-download, add `SignalP` to your system's environmental variables to make it executable.
      
..       - For installation assistance, refer to the :doc:`faq/signalp_installation`.

..    This approach ensures compliance with `SignalP`'s licensing while offering the tool's functionality to those who need it.



.. After the dependencies are installed, `dbcan`_ can be installed via `PyPI <https://pypi.org/>`_:

.. .. code-block:: shell

..     pip install dbcan

.. .. note::

..    Since ``PyPI`` doesn't have an independent build system, the dependencies of dbcan need to be installed seperatedly.
..    Therefore, we recommended users to install ``dbcan`` via ``Conda`` which can resolve all dependencies automatically.

Installing with Docker
----------------------

To use `dbcan`_ via `Docker <https://www.docker.com/>`_, please follow these
steps:

1. Install ``Docker`` on your system (e.g. Linux, MacOS);
2. Pull the image `haidyi/run_dbcan <https://hub.docker.com/r/haidyi/run_dbcan>`_ from `Docker Hub <https://hub.docker.com/>`_;
3. Run the ``run_dbcan`` tool via Docker:

   .. code-block:: shell

      docker run -it haidyi/run_dbcan:latest <input_file> [args] --out_dir <output_dir>

   .. note::

      To use your own local files as input when using Docker, make sure the local files are ``mounted`` and visible to your container.

Check Installation
------------------

After installation, you can check if `dbcan`_ is successfully installed by running:

.. code-block:: shell

   run_dbcan -h

If it shows all the help information, congratulations! You are ready to annotate your own proteins right now.

.. _dbcan: https://github.com/linnabrown/run_dbcan/
.. _Anaconda: https://docs.anaconda.com/free/anaconda/
.. _Miniconda: https://docs.conda.io/projects/miniconda/en/latest/
.. _PyPI: https://pypi.org/
.. _Docker: https://www.docker.com/
