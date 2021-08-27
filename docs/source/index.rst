.. Capture documentation master file, created by
   sphinx-quickstart on Tue Apr 13 13:11:19 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: /images/DOV_FINAL_LOGO_2017_RGB.svg
   :width: 100pt

Welcome to Dovetail Genomics\ :sup:`®` \ Promoter Panel data analysis documentation!
====================================================================================

.. image:: /images/main_schematic.png

Overview
========
- The Dovetail\ :sup:`®` \ Pan-genome promoter panels were designed to include all known coding genes and lncRNA promoters in the human and mouse genomes. The panels can be used in conjugation with either `Dovetail Omni-C <https://omni-c.readthedocs.io/en/latest/>`_\ :sup:`®` \ or `Dovetail Micro-C <https://micro-c.readthedocs.io/en/latest/>`_\ :sup:`®` \ libraries. 
  
- Regulatory Elements: 

  - Transcriptional regulation relies on interactions between protein complexes bound at the transcriptional start site (TSS) and at regulatory sites, including enhancers. 
  - Promoter and enhancer sites can be separated by kbp (or even Mbp) of intervening sequence. Understanding these interactions is key to unravelling gene regulation and related human disease. 
  - Promoter (P) enhancer (E) interactions cannot be detected by standard NGS approaches. For example, while ChIP-seq maps protein binding to promoter or enhancer it is blind to P-E interactions. Another established method, ATAC-Seq maps nucleosome-free regions but does not offer information about P, E or P-E interactions. Hi-C Sequencing methods can be used for chromatin loop calling but at a high depth and limited resolution, requiring about 1B reads for the human genome. 

- Key benefits of Dovetail™ Pan-genome promoter panel:

  - The panel was optimized to work with Dovetail's unique restriction enzyme-free  Hi-C approaches: Omni-C and Micro-C. This means that no promoter was ommited from the panel due to lack of restriction sites in its vicinty as needed in other Capture Hi-C methods. The Pan-genome promoter panel offers highest resolution view of P-E interactions with lowest sequencing burden.
  - Captured Libraries maintain Hi-C characteristics in terms of valid pairs legnth distribution while achieving uniform coverage over baited regions. 
  - The reproducibility of the captured libraries is very high, with high corellation of coverage over baits between replicas :math:`R^2 > 0.97`.
  - The robust coverage of the human and mouse panels enabls genome wide gene regulation and epigenetic research, developmental research, complex trait mapping and more:
  .. image:: /images/panel_design.png
   :width: 500pt
   :align: center


- This guide will take you step by step on how to QC your captured library, how to interpret the :ref:`QC <LQ>` results and evaluate :ref:`reproducibility <RR>` between replicas.  We will show you how to generate :ref:`contact maps <GCM>`, and most importantly how to find :ref:`significant interactions <INT>`, prioritize them and how to combine this information with additional data sets that you might have. If you don't yet have a sequenced Dovetail captured library and you want to get familiar with the data, you can download iPSC and NSC captured Micro-C libraries from our publicaly available :ref:`data sets<DATASETS>`. In the :ref:`data sets<DATASETS>` you will also find bed files with the probes and bait locations.

.. _PROBAIT:

.. image:: /images/QC_enrich.png
   :width: 300pt
   :align: center

- If this is your first time following this tutorial, please check the :ref:`Before you begin page <BYB>` first.


**Main steps in pan-genome captured library preparation and data processing:**

.. figure:: /images/panel_data_viz.png
   :width: 2600pt

.. raw:: html
  
   <iframe title="vimeo-player" src="https://player.vimeo.com/video/563224681?h=625a491c61" width="640" height="327" frameborder="0" allowfullscreen></iframe>


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   before_you_begin
   
   pre_alignment
   
   fastq_to_bam
   
   library_qc

   reproducibility

   interactions
   
   contact_map

   data_sets

   support


.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
