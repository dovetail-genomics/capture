.. _BYB:

Before you begin
================

Have a copy of the Pan Promoter Enrichment Panel capture scripts on your machine:
---------------------------------------------------------------------------------


Clone this repository:

.. code-block:: console

   git clone https://github.com/dovetail-genomics/capture.git


Dependencies
-------------

Make sure that the following dependencies are installed:

- `pysam <https://pysam.readthedocs.io/en/latest/>`_
- `tabulate <https://pypi.org/project/tabulate/>`_
- `bedtools <https://bedtools.readthedocs.io/en/latest/index.html>`_
- `matplotlib <https://matplotlib.org/>`_
- `pandas <https://pandas.pydata.org/pandas-docs/stable/dsintro.html>`_
- `bwa <https://github.com/lh3/bwa>`_
- `pairtools <https://github.com/open2c/pairtools>`_
- `samtools <https://github.com/samtools/samtools>`_
- `mosdepth <https://github.com/brentp/mosdepth>`_

If you are facing any issues with the installation of any of the dependencies, please contact the supporter of the relevant package. 

python3 and pip3 are required, if you don't already have them installed, you will need sudo privileges.

* Update and install python3 and pip3:

.. code-block:: console 

   sudo apt-get update
   sudo apt-get install python3 python3-pip


* To set python3 and pip3 as primary alternative:

.. code-block:: console

   sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
   sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1


If you are working on a new machine and installing all dependencies from scratch, you can use the ``installDep.sh`` script in this repository for updating your instance including python3. This process will take approximately 10 minutes and again requires sudo privileges to complete. The script was tested on Ubuntu 18.04 with the latest version as of 04/11/2020.

If you choose to run the provided installation script you will first need to set the permission to the file:

.. code-block:: console

   chmod +x ./capture/installDep.sh


And then run the installation script:

.. code-block:: console

   ./capture/installDep.sh


.. admonition:: Remember!

   Once the installation is completed, sign off and then sign back on to your instance to refresh the application database.


Input files
-----------

For this tutorial you will need: 

* **Fastq files** R1 and R2 (either fastq or fastq.gz are acceptable)
* **Your genome reference in a fasta file format**, e.g. hg38
* **Coordinates of targeted sequences** in bed format (provided in the :ref:`capture Data Sets section<DATASETS>`). 

.. Tip:: 

   If you don't already have your own input files or want to run a test, you can download sample fastq files from the :ref:`Pan Promoter Enrichment Panel Data Sets section<DATASETS>`. Links to list of probes, baits and reference genomes are also provided, enabling you to reproduce the results presented in this tutorial. 

   The datasets include: two replicas of human induced pleuipotent stem cells (iPSC rep1 and iPSC rep2) and  two replicas of neuronal stem cells (NSC), derived from the same iPSC cells (NSC rep1 and NSC rep 2). 

   The NSC rep1 dataset is used as the main example dataset throughout this tutorial (:ref:`From fastq to bam files <FTB>`, :ref:`QC section <LQ>` etc.), NSC rep1 and NSC rep2 are used in the :ref:`reproducibility section<RR>`, and results from both NSC replicas and iPSC replicas are used for demonstrating :ref:`interaction calling<INT>` and :ref:`identifying differential interactions <DIFF>`.   

