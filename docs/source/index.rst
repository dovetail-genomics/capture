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
The Dovetail\ :sup:`®` \ Pan Promoter Enrichment Panels were designed to include all known coding genes and lncRNA promoters in the human and mouse genomes, respectively. The panels are to be used in conjugation with either `Dovetail Omni-C <https://omni-c.readthedocs.io/en/latest/>`_\ :sup:`®` \ or `Dovetail Micro-C <https://micro-c.readthedocs.io/en/latest/>`_\ :sup:`®` \ libraries for the detection of long-range interactions between gene promoters and their various distal-regulatory elements. 
  
Spatial Control of Transcriptional Regulation: 

  - Transcriptional regulation relies on interactions between protein complexes bound at the transcriptional start site (TSS) and at regulatory sites, including enhancers. 
  - Promoter and enhancer sites can be separated by kilobases (or even megabases) of intervening sequence. Understanding these interactions is key to unravelling gene regulation and related human disease. 
  - Promoter-enhancer (P-E) interactions cannot be detected by standard NGS approaches. For example, transcription factor binding to promoters or enhancers sequences is detectable using ChIP-seq maps protein binding to promoter or enhancer this approach is blind to any higher order P-E interactions that may be mediated by the bound transcription factors. Another established method, ATAC-Seq, provides information about nucleosome-free regions but does not offer information about P, E binding or P-E interactions. Lastly, Hi-C sequencing approaches can be used for chromatin loop calling but require a high sequencing depth and offer limited resolution, requiring about 1B reads for the human genome. 

Key benefits of Dovetail\ :sup:`®` \ Pan Promoter Enrichment Panel:

  - The panel was optimized to work using Dovetail's unique restriction enzyme-free  Hi-C, Omni-C and Micro-C technologies. Extemely high promoter coverage was achieved using this approach as no promoters were ommited from the panel due to lack of restriction sites in their vicinity as needed in other Capture Hi-C methods. The Pan Promoter Enrichment Panel offers the highest resolution view of P-E interactions with lowest sequencing burden.
  - Captured Libraries maintain Hi-C characteristics in terms of valid read-pairs length distribution while achieving uniform coverage over baited regions. 
  - The reproducibility of the captured libraries is very high, with high correlation of coverage over baits between replicas (typical :math:`R^2` are > 0.97).
  - The robust coverage of the human and mouse panels enables genome-wide research in the field of gene regulation and epigenetics, developmental biology, complex trait mapping and more.
  .. image:: /images/panel_design.png
   :width: 500pt
   :align: center


- This guide will take you step-by-step through the :ref:`QC <LQ>` and interpretation of your captured library and evaluation of  :ref:`reproducibility <RR>` between replicas.  We will show you how to:
  
   - Generate :ref:`contact maps <GCM>`
   - :ref:`Find <INT>` and prioritize significant interactions
   - Combine this information with additional data sets that you might have available. 

   Main steps in Pan Promoter captured library preparation and data processing as will be discussed in this document:

   .. figure:: /images/panel_data_viz.png
      :width: 2600pt

  If you have not yet sequenced a Dovetail captured library and simply want to get familiarize yourself with the data, an example dataset from iPSC and NSC captured Micro-C libraries can be downloaded from our publicly-available :ref:`repository<DATASETS>`. In these :ref:`data sets<DATASETS>` you will find bed files with the probes and bait locations.

.. _PROBAIT:

.. image:: /images/QC_enrich.png
   :width: 300pt
   :align: center

- If this is your first time following this tutorial, please check the :ref:`Before you begin page <BYB>` first.

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
