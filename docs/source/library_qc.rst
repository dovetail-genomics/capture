.. _LQ:

Library QC
==========

Proximity ligation properties
-----------------------------

At step :ref:`Removing PCR duplicates<DUPs>` you used the flag ``--output-stats``, generating a stats file in addition to the pairsam output (e.g. --output-stats stats.txt). The stats file is an extensive output of pairs statistics as calculated by pairtools, including total reads, total mapped, total dups, total pairs for each pair of chromosomes etc'. Although you can use directly the pairtools stats file as is to get informed on the quality of the Omni-C\ :sup:`®` \ or Micro-C\ :sup:`®` \ library, we find it easier to focus on a few key metrics. We include in this repository the script ``get_qc.py`` that summarize the paired-tools stats file and present them in percentage values in addition to absolute values.

The images below explains how the values on the QC report are calculated:

.. image:: /images/QC_align.png

.. image:: /images/QC_cis_trans_valids.png

**Command:**

.. code-block:: console

   python3 ./capture/get_qc.py -p <stats.txt>


**Example:**

.. code-block:: console

   python3 ./capture/get_qc.py -p NSC_rep1_stats.txt 


After the script completes, it will print (These are the same results you are expected to see if you ran the NSC rep1 sample):

.. code-block:: console

   Total Read Pairs                              253,341,035  100%
   Unmapped Read Pairs                           14,393,599   5.68%
   Mapped Read Pairs                             195,566,613  77.2%
   PCR Dup Read Pairs                            51,471,702   20.32%
   No-Dup Read Pairs                             144,094,911  56.88%
   No-Dup Cis Read Pairs                         94,499,595   65.58%
   No-Dup Trans Read Pairs                       49,595,316   34.42%
   No-Dup Valid Read Pairs (cis >= 1kb + trans)  125,160,568  86.86%
   No-Dup Cis Read Pairs < 1kb                   18,934,343   13.14%
   No-Dup Cis Read Pairs >= 1kb                  75,565,252   52.44%
   No-Dup Cis Read Pairs >= 10kb                 46,482,376   32.26%

- We recommend a minimum of 150M **Total Read Pairs**
- **PCR Dup Read Pairs** may range from 10% up to 35%
- **No-Dup Trans Read Pairs** are typically around 20% to 30% but lowere or higher values are often observed in good quality libraries (12%-35%)
- **No-Dup Cis Read Pairs >= 1kb** is typically between 50% to 65%



   

Target enrichment QC
--------------------

.. image:: /images/capture_qc.png
   :scale: 20%

To evaluate the level of target enrichment we will inspect the coverage over probes (targeted regions) and the overall fraction of the captured reads in our library. 

On-target rate 
++++++++++++++

We define the on-target rate as the percentage of read pairs that maps to targeted regions. The on-target rate use read pairs and not reads since most of the read pairs in the library are expected to have a large insert size, such that only one read from the pair maps to the targeted region. 

Since DNA fragments can extend beyond the sequnced region and even 30bp match are sufficient for capturing, we treat reads that fall in close proximity to the targeted regions (up to 200bp upstream or downstream from a probe) as on-target reads.

To count the number of on target pairs we will use:

 - :ref:`\*PT.bam <FBAM>` file that you generated in the :ref:`previous section <FTB>` (not the chicago compatible bam) 

 - Bed file with the probes positions padded with 200bp at both side (the files h_probes_200bp_v1.0.bed for human and m_probes_200bp_v1.0.bed for mouse can be found in the :ref:`Data sets section <DATASETS>`).

Count on-target read pairs:

**Command:**

.. code-block:: console

   samtools view <*PT.bam file> -L <padded bed file> -@ <threads> \
   |awk -F "\t" '{print "@"$1}'|sort -u|wc -l 

**Example:**

.. code-block:: console

   samtools view NSC_rep1.PT.bam -L h_probes_200bp_v1.0.bed -@ 16|awk -F "\t" '{print "@"$1}'|sort -u|wc -l 
    
samtools view with the ``-L`` argument allows to extract only reads that mapped to the region of interest, awk command help us parse the file and extract the read ID information, sort command with ``-u`` (unique) argument will remove any multiple occurrences of the same read ID (to avoid counting read1 and read2 of the same pair if both mapped to the target region) and finally ``wc -l`` counts the read IDs in this list.

The example above will output the value: 93171111 (**On Target Read Pairs**)


There is no need to count the total read pairs in the bam file (which represents the total number of pairs, or 100%) as it was already reported by the QC script above, labeled as **No-Dup Read Pairs**, in our example: 144094911

Now you can calculate the on target rate:



.. math::

  \frac{On Target Read Pairs}{No Dup Read Pairs}*100



And in the example above:


.. math::

  \frac{93171111}{144094911}*100=64.7\%


The **on target rate** of the NSC replica1 example library is 64.7%. This is a typical on target rate, although occasionally lower values may be observed, as low as 40% 

Coverage depth
++++++++++++++


There are multiple methods and tools that allow extracting coverage depth from a bam file at different regions. We chose to use the tool `mosdepth <https://github.com/brentp/mosdepth>`_ as we found it to be easy to use and relatively fast.

Use the probe bed file (and bait bed file if desired) to calculate coverage using the position sorted bam file (e.g. mapped.PT.bam, do not use the chicago compatible bam file):

.. _MOS:

**Command:**

.. code-block:: console

   mosdepth -t <threads> -b <bed file> -x <output prefix> -n <bam file>


**Example:**

.. code-block:: console

   mosdepth -t 16 -b h_probes_v1.0.bed -x NSC_rep1_probes -n NSC_rep1.PT.bam

This command will yield multiple output files, specifically, two files that will be useful for QC-ing your libraries are a bed file with mean coverage in each region from the bed file (e.g. NSC_rep1_probes.regions.bed.gz) and a summary output file (e.g. NSC_rep1_probes.mosdepth.summary.txt). The summary file inform us on the mean coverage of the total genome (second to last row) and mean coverage of the total_region (targeted region of interest - the last row in the summary).
To print the header and two last summarizing rows, follow this example: 

.. code-block:: console

   head -1 NSC_rep1_probes.mosdepth.summary.txt;tail -2 NSC_rep1_probes.mosdepth.summary.txt

This will output the following:

.. code-block:: console

   chrom          length      bases       mean     min   max
   total          3088269832  39020721947 12.64    0     482767
   total_region   19337280    7835787504  405.22   0     8129


In this example (NSC rep1) the mean coverage over targeted regions is 405.22, while non-targeted regions have a mean coverage depth of only 12.64. Overall the coverae depth is 32 time higher at targeted regions vs non targeted regions: :math:`405.22/12.64 = 32`. The fold difference between the mean coverage depth of targeted regions and non-targeted regions is typically around 30, just as seen in this example. 

The bed files with mean coverage values at on-target regions (e.g. NSC_rep1_probes.regions.bed.gz and NSC_rep2_probes.regions.bed.gz) will be used to assess :ref:`replica reproducibility <RR>`

|clock| Running the QC steps can be completed in less than 2 hours on an Ubuntu 18.04 machine with 16 CPUs, 1TB storage and 64GiB memory.


.. |clock| image:: /images/clock.jpg
           :scale: 5 %

