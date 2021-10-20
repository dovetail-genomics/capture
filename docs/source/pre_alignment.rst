.. _PA:

Pre-Alignment
=============

For downstream steps you will need a genome file. A genome file is a tab delimited file with chromosome names and their respective sizes. If you do not already have a genome file, follow these steps:

1. Generate an index file for your reference. A reference file with only the main chromosomes should be used (e.g. without alternative or unplaced chromosomes).

**Command:**

.. code-block:: console

   samtools faidx <ref.fasta>


**Example:**

.. code-block:: console

   samtools faidx hg38.fasta

faidx will index the ref file and create <ref.fasta>.fai in the reference directory.

.. _GENOME:

2. Use the index file to generate the genome file by printing the first two columns into a new file.

**Command:**

.. code-block:: console

   cut -f1,2 <ref.fasta.fai> > <ref.genome>


**Example:**

.. code-block:: console

   cut -f1,2 hg38.fasta.fai > hg38.genome


In line with 4DN project guidelines and from our own experience, optimal alignment results are obtained with the Burrows-Wheeler Aligner (bwa).
Prior to alignment, generate a bwa index file for the chosen reference.


.. code-block:: console

   bwa index <ref.fasta>


**Example:**

.. code-block:: console

   bwa index hg38.fasta



There is no need to specify an output path as the bwa index files are automatically generated in the reference directory. Please note that this step is time consuming, however, you need only run it once for your reference. 

To avoid memory issues, some of the steps require writing temporary files into a temp folder. Create a temp folder and remember its full path. Temp files may take up to 3x the space your final fastq.gz files. That is, if the total volume of the fastq files are 5Gb, make sure that the temp folder can store at least 15 Gb.

**Command:**

.. code-block:: console

   mkdir <full_path/to/tmpdir>


**Example:**

.. code-block:: console

   mkdir /home/ubuntu/ebs/temp


In this example the folder `temp` will be generated on a mounted volume called `ebs` on a user account `ubuntu`.
