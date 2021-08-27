.. _FTB:

From fastq to bam files
=======================

.. _Impatient:

.. admonition:: fastq to final valid pairs bam file - for the impatient!

   If you just want to give it a shot and run all the alignment and filtering steps without going over all the details, we made a shorter version for you, with all the steps piped. Otherwise move to the next section :ref:`fastq to final valid pairs bam file - step by step<step-by-step>`. 
   The piped commands outputs two different formats of final bam files, bam index file and a dup stats file. The example below is based on human NSC dataset, replica 1, you can find the fastq files and a link to the reference in the  :ref:`capture Data Sets section<DATASETS>`. Otherwise move to the next section :ref:`fastq to final valid pairs bam file - step by step<step-by-step>`

   **Command:**

   .. code-block:: console

      bwa mem -5SP -T0 -t<cores> <ref.fa> <capture.R1.fastq.gz> <capture.R2.fastq.gz>| \
      pairtools parse --min-mapq 40 --walks-policy 5unique \
      --max-inter-align-gap 30 --nproc-in <cores> --nproc-out <cores> --chroms-path <ref.genome> | \
      pairtools sort --tmpdir=<full_path/to/tmpdir> --nproc <cores>|pairtools dedup --nproc-in <cores> \
      --nproc-out <cores> --mark-dups --output-stats <stats.txt>|pairtools split --nproc-in <cores> \
      --nproc-out <cores> --output-pairs <mapped.pairs> --output-sam -|samtools view -bS -@<cores> | \
      samtools sort -@<cores> -o <mapped.PT.bam>;samtools index <mapped.PT.bam>;samtools view -@ <threads> -Shu -F 2048 <mapped.PT.bam>|samtools sort -n -T <path_to_temp_dir> --threads <threads> -o <chicago.bam> -

   **Example:**

   .. code-block:: console

      bwa mem -5SP -T0 -t16 hg38.fa NSC_rep1_R1.fastq.gz NSC_rep1_R2.fastq.gz| pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hg38.genome | pairtools sort --tmpdir=/home/ubuntu/ebs/temp/ --nproc 16|pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam -|samtools view -bS -@16 | samtools sort -@16 -o NSC_rep1_PT.bam;samtools index NSC_rep1_PT.bam;samtools view -@ 16 -Shu -F 2048 NSC_rep1_PT.bam|samtools sort -n -T /home/ubuntu/ebs/temp/ --threads 16 -o NSC_rep1_chicago.bam -


|clock| The full pipeline, with 250M read pairs on an Ubuntu 18.04 machine with 16 CPUs, 1TB storage and 64GiB memory takes about 8 hours to complete.


.. |clock| image:: /images/clock.jpg
           :scale: 5 %


.. _step-by-step:

fastq to final valid pairs bam file - step by step
--------------------------------------------------


Alignment 
+++++++++

Now that you have a genome file, index file and a reference fasta file you are all set to align your captured Micro-C\ :sup:`®` \ or Omni-C\ :sup:`®` \ library to the reference. Please note the specific settings that are needed to map mates independently and for optimal results with our proximity library reads.

To replicate the output files generated in this workflow, please use the NSC replica 1 fastq files from our :ref:`capture Data Sets section<DATASETS>`, for your convenience reference file is also included. if you're using your own fastq files, the results will be different from the example outputs displayed here. 


.. csv-table::
   :file: tables/alignment.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table


Bwa mem will output a sam file that you can either pipe or save to a path using -o option, as in the example below (please note that version 0.7.17 or higher should be used, older versions do not support the `-5` flag)

**Command:**

.. code-block:: console

   bwa mem -5SP -T0 -t<threads> <ref.fasta> <capture_R1.fastq> <capture_R2.fastq> -o <aligned.sam> 


**Example (one pair of fastq files):**

.. code-block:: console

   bwa mem -5SP -T0 -t16 hg38.fasta NSC_rep1_R1.fastq.gz NSC_rep1_R2.fastq.gz -o  NSC_rep1_aligned.sam


**Example (multiple pairs of fastq files):**

.. code-block:: console

   bwa mem -5SP -T0 -t16 hg38.fasta <(zcat file1.R1.fastq.gz file2.R1.fastq.gz file3.R1.fastq.gz) <(zcat file1.R2.fastq.gz file2.R2.fastq.gz file3.R2.fastq.gz) -o aligned.sam

.. note::

   The bwa command will work on either fastq files or fastq.gz files


Recording valid ligation events
+++++++++++++++++++++++++++++++

We use the ``parse`` module of the ``pairtools`` pipeline to find ligation junctions. When a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair and the strand of each one of the paired reads into ``.pairsam`` file (pairsam format captures SAM entries together with the Hi-C pair information). In addition, it will also assign a pair type for each event. e.g. if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. The following steps are necessary to identify the high quality valid pairs over low quality events (e.g. due to low mapping quality):


``pairtools parse`` options:


.. csv-table::
   :file: tables/parse.csv
   :header-rows: 1
   :widths: 20 20 60
   :class: tight-table


``pairtools parse`` command example for finding ligation events:

**Command:**

.. code-block:: console

   pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in <cores>\
  --nproc-out <cores> --chroms-path <ref.genome> <aligned.sam> > <parsed.pairsam>


**Example:**

.. code-block:: console

   pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hg38.genome NSC_rep1_aligned.sam >   NSC_rep1_parsed.pairsam


At the parsing step, pairs will be flipped such that regardless of read1 and read2, pairs are always recorded with first side of the pair having the lower genomic coordinates. 


Sorting the pairsam file
++++++++++++++++++++++++


The parsed pairs are then sorted using `pairtools sort`

``pairtools sort`` options:

.. csv-table::
   :file: tables/sort.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table

**Command:**

.. code-block:: console

   pairtools sort --nproc <cores> --tmpdir=<path/to/tmpdir> <parsed.pairsam> > <sorted.pairsam>


**Example:**

.. code-block:: console

   pairtools sort --nproc 16 --tmpdir=/home/ubuntu/ebs/temp/  NSC_rep1_parsed.pairsam > NSC_rep1_sorted.pairsam

.. admonition:: Important!

   Please note that an absolute path for the temp directory is required for ``pairtools sort``, e.g. path of the structure ~/ebs/temp/ or ./temp/ will not work, instead, something of this sort is needed /home/user/ebs/temp/

.. _DUPs:

Removing PCR duplicates
+++++++++++++++++++++++

``pairtools dedup`` detects molecules that could be formed via PCR duplication and tags them as “DD” pair type. These pairs should be excluded from downstream analysis. Use the pairtools dedup command with the `--output-stats` option to save the dup stats into a text file.

``pairtools dedup`` options:

.. csv-table::
   :file: tables/dedup.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table

**Command:**

.. code-block:: console

   pairtools dedup --nproc-in <cores> --nproc-out <cores> --mark-dups --output-stats <stats.txt> \
   --output <dedup.pairsam> <sorted.pairsam>


**Example:**

.. code-block:: console

   pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output NSC_rep1_dedup.pairsam NSC_rep1_sorted.pairsam

.. _GPB:

Generating .pairs and bam files
+++++++++++++++++++++++++++++++

The ``pairtools split`` command is used to split the final ``.pairsam`` into two files: ``.sam`` (or ``.bam``) and ``.pairs`` (``.pairsam`` has two extra columns containing the alignments from which the Omni-C or Micro-C pair was extracted, these two columns are not included in ``.pairs`` files)

``pairtools split`` options:

.. csv-table::
   :file: tables/split.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table


**Command:**

.. code-block:: console

   pairtools split --nproc-in <cores> --nproc-out <cores> --output-pairs <mapped.pairs> \
   --output-sam <unsorted.bam> <dedup.pairsam>


**Example:**

.. code-block:: console

   pairtools split --nproc-in 8 --nproc-out 8 --output-pairs NSC_rep1_mapped.pairs --output-sam NSC_rep1_unsorted.bam NSC_rep1_dedup.pairsam

The ``.pairs`` file can be used for generating :ref:`contact matrix <GCM>`

.. _FBAM:

Generating the dedup, sorted bam file
+++++++++++++++++++++++++++++++++++++

For downstream steps, the bam file should be sorted, using the command `samtools sort`

``samtools sort`` options:

.. csv-table::
   :file: tables/bam_sort.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table
 

**Command:**

.. code-block:: console

   samtools sort -@<threads> -T <path/to/tmpdir/>-o <mapped.PT.bam> <unsorted.bam>


**Example:**

.. code-block:: console

   samtools sort -@16 -T /home/ubuntu/ebs/temp/ -o NSC_rep1_PT.bam NSC_rep1_unsorted.bam


For future steps an index (.bai) of the bam file is also needed.
Index the bam file:

**Command:**

.. code-block:: console

   samtools index <mapped.PT.bam>


**Example:**

.. code-block:: console

   samtools index NSC_rep1_PT.bam


The above steps resulted in multiple intermediate files, to simplify the process and avoid intermediate files, you can pipe the steps.

The `*PT.bam` (PT stands for pair tools) is a key bam file that will be used for :ref:`library QC <LQ>`, generating :ref:`contact maps <GCM>` and more. Additional processing of the bam file will be required for :ref:`interaction calling <INT>`.

.. _CHIBAM:

CHiCAGO compatible bam file
---------------------------

As will be discussed in the :ref:`interaction calling <INT>` section, we will use the `CHiCAGO tool <http://functionalgenecontrol.group/chicago>`_ for calling P-E interactions. CHiCAGO is designed to work with bam files produced with `HiCUP <http://www.bioinformatics.babraham.ac.uk/projects/hicup/>`_ pipeline. To match the format of our bam file to the format expected by CHiCAGO, we will clean the bam file from alignments that won't be used by CHiCAGO (e.g. supplementary alignment) and modify the sorting from position based to read name based sorting. 

Samtools parameter for generating a CHiCAGO compatible bam format:

.. csv-table::
   :file: tables/samtools_chic.csv
   :header-rows: 1
   :widths: 20 20 60
   :class: tight-table


**Command:**

.. code-block:: console

   samtools view -@ <threads> -Shu -F 2048 <input bam file>|samtools sort -n -T <path to temp dir> --threads <threads> -o <output bam file> -


**Example:**

.. code-block:: console

   samtools view -@ 16 -Shu -F 2048 NSC_rep1_PT.bam|samtools sort -n -T /home/ubuntu/ebs/temp/temp.bam --threads 16 -o NSC_rep1_chicago.bam -
   

