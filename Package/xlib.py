#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''
This software has been developed by:

    GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiologia e Historia Forestal)
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politecnica de Madrid
    https://github.com/ggfhf/

Licence: GNU General Public Licence Version 3.
'''

#-------------------------------------------------------------------------------

'''
This source contains general functions and classes used in NGScloud2
software package used in both console mode and gui mode.
'''

#-------------------------------------------------------------------------------

import collections
import configparser
import datetime
import os
import re
import requests
import subprocess
import sys
import tkinter

import xconfiguration

#-------------------------------------------------------------------------------
    
def get_project_code():
    '''
    Get the project name.
    '''

    return 'ngscloud2'

#-------------------------------------------------------------------------------
    
def get_project_name():
    '''
    Get the project name.
    '''

    return 'NGScloud2'

#-------------------------------------------------------------------------------

def get_project_version():
    '''
    Get the project name.
    '''

    return '2.12'

#-------------------------------------------------------------------------------
    
def get_project_manual_file():
    '''
    Get the project name.
    '''

    return './NGScloud2-manual.pdf'

#-------------------------------------------------------------------------------
    
def get_project_image_file():
    '''
    Get the project name.
    '''

    return './image_NGScloud2.png'

#-------------------------------------------------------------------------------

def get_starcluster():
    '''
    Get the script to run StarCluster corresponding to the Operating System.
    '''

    # assign the StarCluster script
    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        starcluster = './starcluster.sh'
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        starcluster = r'.\starcluster.bat'

    # return the StarCluster script
    return starcluster

#-------------------------------------------------------------------------------

def get_sci():
    '''
    Get the script to run sci.py corresponding to the Operating System.
    '''

    # assign the sci.py script
    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        sci = './sci.py'
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        sci = r'.\sci.bat'

    # return the sci.py script
    return sci

#-------------------------------------------------------------------------------

def get_editor():
    '''
    Get the editor depending on the Operating System.
    '''

    # assign the editor
    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        editor = 'nano'
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        editor = 'notepad'

    # return the editor
    return editor

#-------------------------------------------------------------------------------

def get_volume_creator_name():
    '''
    Get the name of volume creator.
    '''

    # set the name of volume creator
    volume_creator_name = '{0}-volume-creator'.format(xconfiguration.environment)

    # return the name of volume creator
    return volume_creator_name

#-------------------------------------------------------------------------------

def get_all_applications_selected_code():
    '''
    Get the code that means all applications.
    '''

    return 'all_applications_selected'

#-------------------------------------------------------------------------------

def get_awscli_name():
    '''
    Get the AWS CLI 2 name used to title.
    '''

    return 'AWSCLI2'

#-------------------------------------------------------------------------------

def get_awscli_url():
    '''
    Get the AWS CLI 2 URL.
    '''

    return 'https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip'

#-------------------------------------------------------------------------------

def get_anaconda_name():
    '''
    Get the Anaconda name used to title.
    '''

    return 'Anaconda'

#-------------------------------------------------------------------------------

def get_anaconda_url():
    '''
    Get the Anaconda URL.
    '''

    return 'https://anaconda.org/'

#-------------------------------------------------------------------------------

def get_github_name():
    '''
    Get the GitHub name used to title.
    '''

    return 'GitHub'

#-------------------------------------------------------------------------------

def get_github_url():
    '''
    Get the GitHub URL.
    '''

    return 'https://github.com/'

#-------------------------------------------------------------------------------

def get_bcftools_code():
    '''
    Get the BCFTools code used to identify its processes.
    '''

    return 'bcftools'

#-------------------------------------------------------------------------------

def get_bcftools_name():
    '''
    Get the BCFTools name used to title.
    '''

    return 'BCFtools'

#-------------------------------------------------------------------------------

def get_bcftools_anaconda_code():
    '''
    Get the BCFTools code used to identify the Anaconda package.
    '''

    return 'bcftools'

#-------------------------------------------------------------------------------

def get_bedtools_code():
    '''
    Get the BEDTools code used to identify its processes.
    '''

    return 'bedtools'

#-------------------------------------------------------------------------------

def get_bedtools_name():
    '''
    Get the BEDTools name used to title.
    '''

    return 'BEDtools'

#-------------------------------------------------------------------------------

def get_bedtools_anaconda_code():
    '''
    Get the BEDTools code used to identify the Anaconda package.
    '''

    return 'bedtools'

#-------------------------------------------------------------------------------

def get_blastplus_code():
    '''
    Get the BLAST+ code used to identify its processes.
    '''

    return 'blast'

#-------------------------------------------------------------------------------

def get_blastplus_name():
    '''
    Get the BLAST+ name used to title.
    '''

    return 'BLAST+'

#-------------------------------------------------------------------------------

def get_blastplus_anaconda_code():
    '''
    Get the BLAST+ code used to identify the Anaconda package.
    '''

    return 'blast'

#-------------------------------------------------------------------------------

def get_bowtie2_code():
    '''
    Get the Bowtie2 code used to identify its processes.
    '''

    return 'bowtie2'

#-------------------------------------------------------------------------------

def get_bowtie2_name():
    '''
    Get the Bowtie2 name used to title.
    '''

    return 'Bowtie2'

#-------------------------------------------------------------------------------

def get_bowtie2_anaconda_code():
    '''
    Get the Bowtie2 code used to identify the Anaconda package.
    '''

    return 'bowtie2'

#-------------------------------------------------------------------------------

def get_busco_code():
    '''
    Get the BUSCO code used to identify its processes.
    '''

    return 'busco'

#-------------------------------------------------------------------------------

def get_busco_name():
    '''
    Get the BUSCO name used to title.
    '''

    return 'BUSCO'

#-------------------------------------------------------------------------------

def get_busco_anaconda_code():
    '''
    Get the BUSCO code used to identify the Anaconda package.
    '''

    return 'busco'

#-------------------------------------------------------------------------------

def get_cd_hit_code():
    '''
    Get the CD-HIT code used to identify its processes.
    '''

    return 'cdhit'

#-------------------------------------------------------------------------------

def get_cd_hit_name():
    '''
    Get the CD-HIT name used to title.
    '''

    return 'CD-HIT'

#-------------------------------------------------------------------------------

def get_cd_hit_anaconda_code():
    '''
    Get the CD-HIT code used to identify the Anaconda package.
    '''

    return 'cd-hit'

#-------------------------------------------------------------------------------

def get_cd_hit_est_code():
    '''
    Get the CD-HIT-EST code used to identify its processes.
    '''

    return 'cdhitest'

#-------------------------------------------------------------------------------

def get_cd_hit_est_name():
    '''
    Get the CD-HIT-EST name used to title.
    '''

    return 'CD-HIT-EST'

#-------------------------------------------------------------------------------

def get_conda_code():
    '''
    Get the Conda code used to identify its processes.
    '''

    return 'conda'

#-------------------------------------------------------------------------------

def get_conda_name():
    '''
    Get the Conda name used to title.
    '''

    return 'Conda'

#-------------------------------------------------------------------------------

def get_cuffdiff_code():
    '''
    Get the Cuffdiff (Cufflinks package) code used to identify its processes.
    '''

    return 'cuffdiff'

#-------------------------------------------------------------------------------

def get_cuffdiff_name():
    '''
    Get the Cuffdiff (Cufflinks package) name used to title.
    '''

    return 'Cuffdiff'

#-------------------------------------------------------------------------------

def get_cufflinks_code():
    '''
    Get the Cufflinks code used to identify its processes.
    '''

    return 'cufflinks'

#-------------------------------------------------------------------------------

def get_cufflinks_name():
    '''
    Get the Cufflinks name used to title.
    '''

    return 'Cufflinks'

#-------------------------------------------------------------------------------

def get_cufflinks_anaconda_code():
    '''
    Get the Cufflinks code used to identify the Anaconda package.
    '''

    return 'cufflinks'

#-------------------------------------------------------------------------------

def get_cufflinks_cuffmerge_code():
    '''
    Get the Cufflinks-Cuffmerge (Cufflinks package) code used to identify its processes.
    '''

    return 'cufflnkmrg'

#-------------------------------------------------------------------------------

def get_cufflinks_cuffmerge_name():
    '''
    Get the Cufflinks-Cuffmerge (Cufflinks package) name used to title.
    '''

    return 'Cufflinks-Cuffmerge'

#-------------------------------------------------------------------------------

def get_cuffnorm_code():
    '''
    Get the Cuffnorm (Cufflinks package) code used to identify its processes.
    '''

    return 'cuffnorm'

#-------------------------------------------------------------------------------

def get_cuffnorm_name():
    '''
    Get the Cuffnorm (Cufflinks package) name used to title.
    '''

    return 'Cuffnorm'

#-------------------------------------------------------------------------------

def get_cuffquant_code():
    '''
    Get the Cuffquant (Cufflinks package) code used to identify its processes.
    '''

    return 'cuffquant'

#-------------------------------------------------------------------------------

def get_cuffquant_name():
    '''
    Get the Cuffquant (Cufflinks package) name used to title.
    '''

    return 'Cuffquant'

#-------------------------------------------------------------------------------

def get_cutadapt_code():
    '''
    Get the cutadapt code used to identify its processes.
    '''

    return 'cutadapt'

#-------------------------------------------------------------------------------

def get_cutadapt_name():
    '''
    Get the cutadapt name used to title.
    '''

    return 'cutadapt'

#-------------------------------------------------------------------------------

def get_cutadapt_anaconda_code():
    '''
    Get the cutadapt code used to identify the Anaconda package.
    '''

    return 'cutadapt'

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_code():
    '''
    Get the ddRADseq simulation (ddRADseqTools package) code used to identify its
    processes.
    '''

    return 'ddradseqsm'

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_name():
    '''
    Get the ddRADseq simulation (ddRADseqTools package) name used to title.
    '''

    return 'ddRADseq simulation'

#-------------------------------------------------------------------------------

def get_ddradseqtools_code():
    '''
    Get the ddRADseqTools code used to identify its processes.
    '''

    return 'ddradseqtl'

#-------------------------------------------------------------------------------

def get_ddradseqtools_name():
    '''
    Get the ddRADseqTools name used to title.
    '''

    return 'ddRADseqTools'

#-------------------------------------------------------------------------------

def get_detonate_code():
    '''
    Get the DETONATE code used to identify its processes.
    '''

    return 'detonate'

#-------------------------------------------------------------------------------

def get_detonate_name():
    '''
    Get the DETONATE name used to title.
    '''

    return 'DETONATE'

#-------------------------------------------------------------------------------

def get_detonate_anaconda_code():
    '''
    Get the DETONATE code used to identify the Anaconda package.
    '''

    return 'detonate'

#-------------------------------------------------------------------------------

def get_diamond_code():
    '''
    Get the DIAMOND code used to identify its processes.
    '''

    return 'diamond'

#-------------------------------------------------------------------------------

def get_diamond_name():
    '''
    Get the DIAMOND name used to title.
    '''

    return 'DIAMOND'

#-------------------------------------------------------------------------------

def get_diamond_anaconda_code():
    '''
    Get the DIAMOND code used to identify the Anaconda package.
    '''

    return 'diamond'

#-------------------------------------------------------------------------------

def get_emboss_code():
    '''
    Get the EMBOSS code used to identify its processes.
    '''

    return 'emboss'

#-------------------------------------------------------------------------------

def get_emboss_name():
    '''
    Get the EMBOSS name used to title.
    '''

    return 'EMBOSS'

#-------------------------------------------------------------------------------

def get_emboss_anaconda_code():
    '''
    Get the EMBOSS code used to identify the Anaconda package
    '''

    return 'emboss'

#-------------------------------------------------------------------------------

def get_entrez_direct_code():
    '''
    Get the Entrez Direct code used to identify its processes.
    '''

    return 'edirect'

#-------------------------------------------------------------------------------

def get_entrez_direct_name():
    '''
    Get the Entrez Direct name used to title.
    '''

    return 'Entrez Direct'

#-------------------------------------------------------------------------------

def get_entrez_direct_anaconda_code():
    '''
    Get the Entrez Direct code used to the Anaconda package.
    '''

    return 'entrez-direct'

#-------------------------------------------------------------------------------

def get_express_code():
    '''
    Get the eXpress code used to identify its processes.
    '''

    return 'express'

#-------------------------------------------------------------------------------

def get_express_name():
    '''
    Get the eXpress name used to title.
    '''

    return 'eXpress'

#-------------------------------------------------------------------------------

def get_express_anaconda_code():
    '''
    Get the eXpress code used to identify the Anaconda package.
    '''

    return 'express'

#-------------------------------------------------------------------------------

def get_fastqc_code():
    '''
    Get the FastQC code used to identify its processes.
    '''

    return 'fastqc'

#-------------------------------------------------------------------------------

def get_fastqc_name():
    '''
    Get the FastQC name used to title.
    '''

    return 'FastQC'

#-------------------------------------------------------------------------------

def get_fastqc_anaconda_code():
    '''
    Get the FastQC code used to identify the Anaconda package.
    '''

    return 'fastqc'

#-------------------------------------------------------------------------------

def get_ggtrinity_code():
    '''
    Get the Genome-guided Trinity (Trinity package) code used to identify its
    processes.
    '''

    return 'ggtrinity'

#-------------------------------------------------------------------------------

def get_ggtrinity_name():
    '''
    Get the Genome-guided Trinity (Trinity package) name used to title.
    '''

    return 'Genome-guided Trinity'

#-------------------------------------------------------------------------------

def get_gmap_gsnap_code():
    '''
    Get the GMAP-GSNAP code used to identify its processes.
    '''

    return 'gmap_gsnap'

#-------------------------------------------------------------------------------

def get_gmap_gsnap_name():
    '''
    Get the GMAP-GSNAP name used to title.
    '''

    return 'GMAP-GSNAP'

#-------------------------------------------------------------------------------

def get_gmap_gsnap_anaconda_code():
    '''
    Get the GMAP-GSNAP code used to identify the Anaconda package.
    '''

    return 'gmap'

#-------------------------------------------------------------------------------

def get_gmap_code():
    '''
    Get the GMAP code used to identify its processes.
    '''

    return 'gmap'

#-------------------------------------------------------------------------------

def get_gmap_name():
    '''
    Get the GMAP name used to title.
    '''

    return 'GMAP'

#-------------------------------------------------------------------------------

def get_gsnap_code():
    '''
    Get the GSNAP code used to identify its processes.
    '''

    return 'gsnap'

#-------------------------------------------------------------------------------

def get_gsnap_name():
    '''
    Get the GSNAP name used to title.
    '''

    return 'GSNAP'

#-------------------------------------------------------------------------------

def get_gzip_code():
    '''
    Get the gzip code used to identify its processes.
    '''

    return 'gzip'

#-------------------------------------------------------------------------------

def get_gzip_name():
    '''
    Get the gzip name used to title.
    '''

    return 'gzip'

#-------------------------------------------------------------------------------

def get_hisat2_code():
    '''
    Get the HISAT2 code used to identify its processes.
    '''

    return 'hisat2'

#-------------------------------------------------------------------------------

def get_hisat2_name():
    '''
    Get the HISAT2 name used to title.
    '''

    return 'HISAT2'

#-------------------------------------------------------------------------------

def get_hisat2_anaconda_code():
    '''
    Get the HISAT2 code used to identify the Anaconda package.
    '''

    return 'hisat2'

#-------------------------------------------------------------------------------

def get_htseq_code():
    '''
    Get the HTSeq code used to identify its processes.
    '''

    return 'htseq'

#-------------------------------------------------------------------------------

def get_htseq_name():
    '''
    Get the HTSeq name used to title.
    '''

    return 'HTSeq'

#-------------------------------------------------------------------------------

def get_htseq_anaconda_code():
    '''
    Get the HTSeq code used to identify the Anaconda package.
    '''

    return 'htseq'

#-------------------------------------------------------------------------------

def get_htseq_count_code():
    '''
    Get the htseq-count code used to identify its processes.
    '''

    return 'htseqcount'

#-------------------------------------------------------------------------------

def get_htseq_count_name():
    '''
    Get the htseq-count name used to title.
    '''

    return 'htseq-count'

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_code():
    '''
    Get the insilico_read_normalization (Trinity package) code used to identify its
    processes.
    '''

    return 'insreadnor'

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_name():
    '''
    Get the insilico_read_normalization (Trinity package) name used to title.
    '''

    return 'insilico_read_normalization'

#-------------------------------------------------------------------------------

def get_ipyrad_code():
    '''
    Get the ipyrad code used to identify its processes.
    '''

    return 'ipyrad'

#-------------------------------------------------------------------------------

def get_ipyrad_name():
    '''
    Get the ipyrad name used to title.
    '''

    return 'ipyrad'

#-------------------------------------------------------------------------------

def get_ipyrad_conda_code():
    '''
    Get the ipyrad code used to identify the Anaconda package.
    '''

    return 'ipyrad'

#-------------------------------------------------------------------------------

def get_kallisto_code():
    '''
    Get the kallisto code used to identify its processes.
    '''

    return 'kallisto'

#-------------------------------------------------------------------------------

def get_kallisto_name():
    '''
    Get the kallisto name used to title.
    '''

    return 'kallisto'

#-------------------------------------------------------------------------------

def get_kallisto_anaconda_code():
    '''
    Get the kallisto code used to identify the Anaconda package.
    '''

    return 'kallisto'

#-------------------------------------------------------------------------------

def get_miniconda3_code():
    '''
    Get the Miniconda3 code used to identify its processes.
    '''

    return 'miniconda3'

#-------------------------------------------------------------------------------

def get_miniconda3_name():
    '''
    Get the Miniconda3 name used to title.
    '''

    return 'Miniconda3'

#-------------------------------------------------------------------------------

def get_ngshelper_code():
    '''
    Get the NGShelper code used to identify its processes.
    '''

    return 'ngshelper'

#-------------------------------------------------------------------------------

def get_ngshelper_name():
    '''
    Get the NGShelper name used to title.
    '''

    return 'NGShelper'

#-------------------------------------------------------------------------------

def get_quast_code():
    '''
    Get the QUAST code used to identify process.
    '''

    return 'quast'

#-------------------------------------------------------------------------------

def get_quast_name():
    '''
    Get the QUAST name used to title.
    '''

    return 'QUAST'

#-------------------------------------------------------------------------------

def get_quast_anaconda_code():
    '''
    Get the QUAST code used to identify the Anaconda package.
    '''

    return 'quast'

#-------------------------------------------------------------------------------

def get_r_code():
    '''
    Get the R code used to identify its processes.
    '''

    return 'r'

#-------------------------------------------------------------------------------

def get_r_name():
    '''
    Get the R name used to title.
    '''

    return 'R'

#-------------------------------------------------------------------------------

def get_r_anaconda_code():
    '''
    Get the R code used to identify the Anaconda package.
    '''

    return 'R'

#-------------------------------------------------------------------------------

def get_raddesigner_code():
    '''
    Get the RADdesigner code used to identify its processes.
    '''

    return 'raddesigner'

#-------------------------------------------------------------------------------

def get_raddesigner_name():
    '''
    Get the RADdesigner name used to title.
    '''

    return 'RADdesigner'

#-------------------------------------------------------------------------------

def get_ref_eval_code():
    '''
    Get the REF-EVAL (DETONATE package) code used to identify its processes.
    '''

    return 'refeval'

#-------------------------------------------------------------------------------

def get_ref_eval_name():
    '''
    Get the REF-EVAL (DETONATE package) name used to title.
    '''

    return 'REF-EVAL'

#-------------------------------------------------------------------------------

def get_rnaquast_code():
    '''
    Get the rnaQUAST code used to identify its processes.
    '''

    return 'rnaquast'

#-------------------------------------------------------------------------------

def get_rnaquast_name():
    '''
    Get the rnaQUAST name used to title.
    '''

    return 'rnaQUAST'

#-------------------------------------------------------------------------------

def get_rnaquast_anaconda_code():
    '''
    Get the rnaQUAST code used to the Anaconda package.
    '''

    return 'rnaquast'

#-------------------------------------------------------------------------------

def get_rsem_code():
    '''
    Get the RSEM code used to identify its processes.
    '''

    return 'rsem'

#-------------------------------------------------------------------------------

def get_rsem_name():
    '''
    Get the RSEM name used to title.
    '''

    return 'RSEM'

#-------------------------------------------------------------------------------

def get_rsem_anaconda_code():
    '''
    Get the RSEM code used to identify the Anaconda package.
    '''

    return 'rsem'

#-------------------------------------------------------------------------------

def get_rsem_eval_code():
    '''
    Get the RSEM-EVAL (DETONATE package) code used to identify its processes.
    '''

    return 'rsemeval'

#-------------------------------------------------------------------------------

def get_rsem_eval_name():
    '''
    Get the RSEM-EVAL (DETONATE package) name used to title.
    '''

    return 'RSEM-EVAL'

#-------------------------------------------------------------------------------

def get_rsitesearch_code():
    '''
    Get the rsitesearch (ddRADseqTools package) code used to identify its
    processes.
    '''

    return 'rsitesearc'

#-------------------------------------------------------------------------------

def get_rsitesearch_name():
    '''
    Get the rsitesearch (ddRADseqTools package) name used to title.
    '''

    return 'rsitesearch'

#-------------------------------------------------------------------------------

def get_samtools_code():
    '''
    Get the SAMtools code used to identify its processes.
    '''

    return 'samtools'

#-------------------------------------------------------------------------------

def get_samtools_name():
    '''
    Get the SAMtools name used to title.
    '''

    return 'SAMtools'

#-------------------------------------------------------------------------------

def get_samtools_anaconda_code():
    '''
    Get the SAMtools code used to identify the Anaconda package.
    '''

    return 'samtools'

#-------------------------------------------------------------------------------

def get_soapdenovo2_code():
    '''
    Get the SOAPdenovo2 code used to identify its processes.
    '''

    return 'sdn2'

#-------------------------------------------------------------------------------

def get_soapdenovo2_name():
    '''
    Get the SOAPdenovo2 name used to title.
    '''

    return 'SOAPdenovo2'

#-------------------------------------------------------------------------------

def get_soapdenovo2_anaconda_code():
    '''
    Get the SOAPdenovo2 code used to identify the Anaconda package.
    '''

    return 'soapdenovo2'

#-------------------------------------------------------------------------------

def get_soapdenovotrans_code():
    '''
    Get the SOAPdenovo-Trans code used to identify its processes.
    '''

    return 'sdnt'

#-------------------------------------------------------------------------------

def get_soapdenovotrans_name():
    '''
    Get the SOAPdenovo-Trans name used to title.
    '''

    return 'SOAPdenovo-Trans'

#-------------------------------------------------------------------------------

def get_soapdenovotrans_anaconda_code():
    '''
    Get the SOAPdenovo-Trans code used to identify the Anaconda package.
    '''

    return 'soapdenovo-trans'

#-------------------------------------------------------------------------------

def get_star_code():
    '''
    Get the STAR code used to identify its processes.
    '''

    return 'star'

#-------------------------------------------------------------------------------

def get_star_name():
    '''
    Get the STAR name used to title.
    '''

    return 'STAR'

#-------------------------------------------------------------------------------

def get_star_anaconda_code():
    '''
    Get the STAR code used to identify the Anaconda package.
    '''

    return 'star'

#-------------------------------------------------------------------------------

def get_starcode_code():
    '''
    Get the starcode code used to identify its processes.
    '''

    return 'starcode'

#-------------------------------------------------------------------------------

def get_starcode_name():
    '''
    Get the starcode name used to title.
    '''

    return 'starcode'

#-------------------------------------------------------------------------------

def get_starcode_anaconda_code():
    '''
    Get the starcode code used to the Anaconda package.
    '''

    return 'starcode'

#-------------------------------------------------------------------------------

def get_tabix_code():
    '''
    Get the tabix code used to identify its processes.
    '''

    return 'tabix'

#-------------------------------------------------------------------------------

def get_tabix_name():
    '''
    Get the tabix name used to title.
    '''

    return 'Tabix'

#-------------------------------------------------------------------------------

def get_tabix_anaconda_code():
    '''
    Get the tabix code used to identify the Anaconda package.
    '''

    return 'tabix'

#-------------------------------------------------------------------------------

def get_tophat_code():
    '''
    Get the TopHat code used to identify its processes.
    '''

    return 'tophat'

#-------------------------------------------------------------------------------

def get_tophat_name():
    '''
    Get the TopHat name used to title.
    '''

    return 'TopHat'

#-------------------------------------------------------------------------------

def get_tophat_anaconda_code():
    '''
    Get the TopHat code used to identify the Anaconda package.
    '''

    return 'tophat'

#-------------------------------------------------------------------------------

def get_transabyss_code():
    '''
    Get the Trans-ABySS code used to identify its processes.
    '''

    return 'transabyss'

#-------------------------------------------------------------------------------

def get_transabyss_name():
    '''
    Get the Trans-ABySS name used to title.
    '''

    return 'Trans-ABySS'

#-------------------------------------------------------------------------------

def get_transabyss_anaconda_code():
    '''
    Get the Trans-ABySS code used to the Anaconda package.
    '''

    return 'transabyss'

#-------------------------------------------------------------------------------

def get_transcript_filter_code():
    '''
    Get the transcripts-filter (NGShelper package) code used to identify its
    processes.
    '''

    return 'transfil'

#-------------------------------------------------------------------------------

def get_transcript_filter_name():
    '''
    Get the transcripts-filter (NGShelper package) name used to title.
    '''

    return 'transcript-filter'

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_code():
    '''
    Get the transcriptome-blastx (NGShelper package) code used to identify its
    processes.
    '''

    return 'transbastx'

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_name():
    '''
    Get the transcriptome-blastx (NGShelper package) name used to title.
    '''

    return 'transcriptome-blastx'

#-------------------------------------------------------------------------------

def get_transdecoder_code():
    '''
    Get the TransDecoder code used to identify its processes.
    '''

    return 'transdecod'

#-------------------------------------------------------------------------------

def get_transdecoder_name():
    '''
    Get the TransDecoder name used to title.
    '''

    return 'TransDecoder'

#-------------------------------------------------------------------------------

def get_transdecoder_anaconda_code():
    '''
    Get the TransDecoder code used to the Anaconda package.
    '''

    return 'transdecoder'

#-------------------------------------------------------------------------------

def get_transrate_code():
    '''
    Get the Transrate code used to identify its processes.
    '''

    return 'transrate'

#-------------------------------------------------------------------------------

def get_transrate_name():
    '''
    Get the Transrate name used to title.
    '''

    return 'Transrate'

#-------------------------------------------------------------------------------

def get_transrate_anaconda_code():
    '''
    Get the Transrate code used to the Anaconda package.
    '''

    return 'transrate'

#-------------------------------------------------------------------------------

def get_transrate_tools_code():
    '''
    Get the Transrate tools code used to identify its processes.
    '''

    return 'transratools'

#-------------------------------------------------------------------------------

def get_transrate_tools_name():
    '''
    Get the Transrate Tools name used to title.
    '''

    return 'Transrate tools'

#-------------------------------------------------------------------------------

def get_transrate_tools_anaconda_code():
    '''
    Get the Transrate Tools code used to the Anaconda package.
    '''

    return 'transrate-tools'

#-------------------------------------------------------------------------------

def get_trimmomatic_code():
    '''
    Get the Trimmomatic code used to identify its processes.
    '''

    return 'trimmo'

#-------------------------------------------------------------------------------

def get_trimmomatic_name():
    '''
    Get the Trimmomatic name used to title.
    '''

    return 'Trimmomatic'

#-------------------------------------------------------------------------------

def get_trimmomatic_anaconda_code():
    '''
    Get the Trimmomatic code used to the Anaconda package.
    '''

    return 'trimmomatic'

#-------------------------------------------------------------------------------

def get_trinity_code():
    '''
    Get the Trinity code used to identify its processes.
    '''

    return 'trinity'

#-------------------------------------------------------------------------------

def get_trinity_name():
    '''
    Get the Trinity name used to title.
    '''

    return 'Trinity'

#-------------------------------------------------------------------------------

def get_trinity_anaconda_code():
    '''
    Get the Trinity code used to the Anaconda package.
    '''

    return 'trinity'

#-------------------------------------------------------------------------------

def get_variant_calling_code():
    '''
    Get the Variant calling (ddRADseqTools package) code used to identify its
    processes.
    '''

    return 'varcalling'

#-------------------------------------------------------------------------------

def get_variant_calling_name():
    '''
    Get the Variant calling (ddRADseqTools package) name used to title.
    '''

    return 'Variant calling'

#-------------------------------------------------------------------------------

def get_vcftools_code():
    '''
    Get the VCFtools code used to identify its processes.
    '''

    return 'vcftools'

#-------------------------------------------------------------------------------

def get_vcftools_name():
    '''
    Get the VCFtools name used to title.
    '''

    return 'VCFtools'

#-------------------------------------------------------------------------------

def get_vcftools_anaconda_code():
    '''
    Get the VCFtools code used to identify the Anaconda package.
    '''

    return 'vcftools'

#-------------------------------------------------------------------------------

def get_vcftools_perl_libraries_code():
    '''
    Get the VCFtools Perl libraries code used to identify its processes.
    '''

    return 'vcfperl'

#-------------------------------------------------------------------------------

def get_vcftools_perl_libraries_name():
    '''
    Get the VCFtools Perl libraries name used to title.
    '''

    return 'VCFtools Perl libraries'

#-------------------------------------------------------------------------------

def get_vcftools_perl_libraries_anaconda_code():
    '''
    Get the VCFtools Perl libraries code used to identify the Anaconda package.
    '''

    return 'perl-vcftools-vcf'

#-------------------------------------------------------------------------------

def get_vsearch_code():
    '''
    Get the VSEARCH code used to identify process.
    '''

    return 'vsearch'

#-------------------------------------------------------------------------------

def get_vsearch_name():
    '''
    Get the VSEARCH name used to title.
    '''

    return 'VSEARCH'

#-------------------------------------------------------------------------------

def get_vsearch_anaconda_code():
    '''
    Get the VSEARCH code used to identify the Anaconda package.
    '''

    return 'vsearch'

#-------------------------------------------------------------------------------

def get_toa_code():
    '''
    Get the TOA code used to identify its processes.
    '''

    return 'toa'

#-------------------------------------------------------------------------------

def get_toa_name():
    '''
    Get the TOA name used to title.
    '''

    return 'TOA'

#-------------------------------------------------------------------------------

def get_toa_database_dir():
    '''
    Get the directory where database data are saved.
    '''

    return 'TOA-databases'

#-------------------------------------------------------------------------------

def get_toa_data_basic_data_code():
    '''
    Get the code used to identify basic data in TOA processes.
    '''

    return 'basic data'

#-------------------------------------------------------------------------------

def get_toa_data_basic_data_name():
    '''
    Get the code used to title basic data in TOA processes.
    '''

    return 'basic data'

#-------------------------------------------------------------------------------

def get_toa_data_dicots_04_code():
    '''
    Get the code used to identify Dicots PLAZA 4.0 in TOA processes.
    '''

    return 'dicots_04'

#-------------------------------------------------------------------------------

def get_toa_data_dicots_04_name():
    '''
    Get the code used to identify Dicots PLAZA 4.0 in TOA processes.
    '''

    return 'Dicots PLAZA 4.0'

#-------------------------------------------------------------------------------

def get_toa_data_gene_code():
    '''
    Get the code used to title NCBI Gene in TOA processes.
    '''

    return 'gene'

#-------------------------------------------------------------------------------

def get_toa_data_gene_name():
    '''
    Get the code used to identify NCBI Gene in TOA processes.
    '''

    return 'NCBI Gene'

#-------------------------------------------------------------------------------

def get_toa_data_go_code():
    '''
    Get the code used to identify Gene Ontology in TOA processes.
    '''

    return 'go'

#-------------------------------------------------------------------------------

def get_toa_data_go_name():
    '''
    Get the code used to title Gene Ontology in TOA processes.
    '''

    return 'Gene Ontology'

#-------------------------------------------------------------------------------

def get_toa_data_gymno_01_code():
    '''
    Get the code used to identify Gymno PLAZA 1.0 in TOA processes.
    '''

    return 'gymno_01'

#-------------------------------------------------------------------------------

def get_toa_data_gymno_01_name():
    '''
    Get the code used to title Gymno PLAZA 1.0 in TOA processes.
    '''

    return 'Gymno PLAZA 1.0'

#-------------------------------------------------------------------------------

def get_toa_data_interpro_code():
    '''
    Get the code used to title InterPro in TOA processes.
    '''

    return 'interpro'

#-------------------------------------------------------------------------------

def get_toa_data_interpro_name():
    '''
    Get the code used to title InterPro in TOA processes.
    '''

    return 'InterPro'

#-------------------------------------------------------------------------------

def get_toa_data_monocots_04_code():
    '''
    Get the code used to identify Monocots PLAZA 4.0 in TOA processes.
    '''

    return 'monocots_04'

#-------------------------------------------------------------------------------

def get_toa_data_monocots_04_name():
    '''
    Get the code used to title Monocots PLAZA 4.0 in TOA processes.
    '''

    return 'Monocots PLAZA 4.0'

#-------------------------------------------------------------------------------

def get_toa_data_nr_code():
    '''
    Get the code used to identify NCBI BLAST database NR in TOA processes.
    '''

    return 'nr'

#-------------------------------------------------------------------------------

def get_toa_data_nr_name():
    '''
    Get the code used to title NCBI BLAST database NR in TOA processes.
    '''

    return 'NCBI BLAST database NR'

#-------------------------------------------------------------------------------

def get_toa_data_nt_code():
    '''
    Get the code used to identify NCBI BLAST database NT in TOA processes.
    '''

    return 'nt'

#-------------------------------------------------------------------------------

def get_toa_data_nt_name():
    '''
    Get the code used to title NCBI BLAST database NT in TOA processes.
    '''

    return 'NCBI BLAST database NT'

#-------------------------------------------------------------------------------

def get_toa_data_refseq_plant_code():
    '''
    Get the code used to identify NCBI RefSeq Plant in TOA processes.
    '''

    return 'refseq_plant'

#-------------------------------------------------------------------------------

def get_toa_data_refseq_plant_name():
    '''
    Get the name used to title NCBI RefSeq Plant in TOA processes.
    '''

    return 'NCBI RefSeq Plant'

#-------------------------------------------------------------------------------

def get_toa_data_taxonomy_code():
    '''
    Get the code used to title NCBI Taxonomy in TOA processes.
    '''

    return 'taxonomy'

#-------------------------------------------------------------------------------

def get_toa_data_taxonomy_name():
    '''
    Get the code used to identify NCBI Taxonomy in TOA processes.
    '''

    return 'NCBI Taxonomy'

#-------------------------------------------------------------------------------

def get_toa_data_viridiplantae_nucleotide_gi_code():
    '''
    Get the code used to identify NCBI Nucleotide GenInfo viridiplantae identifier list in TOA processes.
    '''

    return 'viridiplantae_nucleotide_gi'

#-------------------------------------------------------------------------------

def get_toa_data_viridiplantae_nucleotide_gi_name():
    '''
    Get the code used to title NCBI Nucleotide GenInfo viridiplantae identifier list in TOA processes.
    '''

    return 'NCBI Nucleotide GenInfo viridiplantae identifier list'

#-------------------------------------------------------------------------------

def get_toa_data_viridiplantae_protein_gi_code():
    '''
    Get the code used to identify NCBI Protein GenInfo viridiplantae identifier list in TOA processes.
    '''

    return 'viridiplantae_protein_gi'

#-------------------------------------------------------------------------------

def get_toa_data_viridiplantae_protein_gi_name():
    '''
    Get the code used to title NCBI Protein GenInfo viridiplantae identifier list in TOA processes.
    '''

    return 'NCBI Protein GenInfo viridiplantae identifier list'

#-------------------------------------------------------------------------------

def get_toa_process_download_basic_data_code():
    '''
    Get the code used to identify processes to download Gene Ontology functional annotations.
    '''

    return 'toaddbasic'

#-------------------------------------------------------------------------------

def get_toa_process_download_basic_data_name():
    '''
    Get the name used to title processes to download Gene Ontology functional annotations.
    '''

    return 'Download basic data'

#-------------------------------------------------------------------------------

def get_toa_process_download_dicots_04_code():
    '''
    Get the code used to identify processes to download Dicots PLAZA 4.0 functional annotations.
    '''

    return 'toadddicots04'

#-------------------------------------------------------------------------------

def get_toa_process_download_dicots_04_name():
    '''
    Get the name used to title processes to download Dicots PLAZA 4.0 functional annotations.
    '''

    return 'Download Dicots PLAZA 4.0 funcional annotations'

#-------------------------------------------------------------------------------

def get_toa_process_download_gene_code():
    '''
    Get the code used to identify processes to download NCBI Gene functional annotations.
    '''

    return 'toaddgene'

#-------------------------------------------------------------------------------

def get_toa_process_download_gene_name():
    '''
    Get the name used to title processes to download NCBI Gene functional annotations.
    '''

    return 'Download NCBI Gene funcional annotations'

#-------------------------------------------------------------------------------

def get_toa_process_download_go_code():
    '''
    Get the code used to identify processes to download Gene Ontology functional annotations.
    '''

    return 'toaddgo'

#-------------------------------------------------------------------------------

def get_toa_process_download_go_name():
    '''
    Get the name used to title processes to download Gene Ontology functional annotations.
    '''

    return 'Download Gene Ontology funcional annotations'

#-------------------------------------------------------------------------------

def get_toa_process_download_gymno_01_code():
    '''
    Get the code used to identify processes to download Gymno PLAZA 1.0 functional annotations.
    '''

    return 'toaddgymno01'

#-------------------------------------------------------------------------------

def get_toa_process_download_gymno_01_name():
    '''
    Get the name used to title processes to download Gymno PLAZA 1.0 functional annotations.
    '''

    return 'Download Gymno PLAZA 1.0 funcional annotations'

#-------------------------------------------------------------------------------

def get_toa_process_download_interpro_code():
    '''
    Get the code used to identify process to download InterPro functional annotations.
    '''

    return 'toaddinterpro'

#-------------------------------------------------------------------------------

def get_toa_process_download_interpro_name():
    '''
    Get the name used to title process to download InterPro functional annotations.
    '''

    return 'Download InterPro funcional annotations'

#-------------------------------------------------------------------------------

def get_toa_process_download_monocots_04_code():
    '''
    Get the code used to identify process to download Monocots PLAZA 4.0 functional annotations.
    '''

    return 'toaddmonocots04'

#-------------------------------------------------------------------------------

def get_toa_process_download_monocots_04_name():
    '''
    Get the name used to title process to download Monocots PLAZA 4.0 functional annotations.
    '''

    return 'Download Monocots PLAZA 4.0 funcional annotations'
#-------------------------------------------------------------------------------

def get_toa_process_download_taxonomy_code():
    '''
    Get the code used to identify processes to download NCBI Taxonomy data.
    '''

    return 'toaddtaxo'

#-------------------------------------------------------------------------------

def get_toa_process_download_taxonomy_name():
    '''
    Get the name used to title processes to download NCBI Taxonomy data.
    '''

    return 'Download NCBI Taxonomy data'

#-------------------------------------------------------------------------------

def get_toa_process_gilist_viridiplantae_nucleotide_gi_code():
    '''
    Get the code used to identifiy processes to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    return 'toablvpntgi'

#-------------------------------------------------------------------------------

def get_toa_process_gilist_viridiplantae_nucleotide_gi_name():
    '''
    Get the name used to title processes to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    return 'Build NCBI Nucleotide GenInfo viridiplantae identifier list'

#-------------------------------------------------------------------------------

def get_toa_process_gilist_viridiplantae_protein_gi_code():
    '''
    Get the code used to identify processes to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    return 'toablvpprgi'

#-------------------------------------------------------------------------------

def get_toa_process_gilist_viridiplantae_protein_gi_name():
    '''
    Get the name used to title processes to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    return 'Build NCBI Protein GenInfo viridiplantae identifier list'

#-------------------------------------------------------------------------------

def get_toa_process_load_basic_data_code():
    '''
    Get the code used to identify processes to load basic data.
    '''

    return 'toaldbasic'

#-------------------------------------------------------------------------------

def get_toa_process_load_basic_data_name():
    '''
    Get the name to title processes to load basic data.
    '''

    return 'Load basic data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_dicots_04_code():
    '''
    Get the code used to identify processes to load the Dicots PLAZA 4.0  data.
    '''

    return 'toalddicots04'

#-------------------------------------------------------------------------------

def get_toa_process_load_dicots_04_name():
    '''
    Get the name used to title processes to load the Dicots PLAZA 4.0  data.
    '''

    return 'Load Dicots PLAZA 4.0  data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_gene_code():
    '''
    Get the code used to identify processes to load the NCBI Gene data load.
    '''

    return 'toaldgene'

#-------------------------------------------------------------------------------

def get_toa_process_load_gene_name():
    '''
    Get the name used to title processes to load the NCBI Gene data load.
    '''

    return 'Load NCBI Gene data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_go_code():
    '''
    Get the code used to identify processes to load the Gene Ontology data.
    '''

    return 'toaldgo'

#-------------------------------------------------------------------------------

def get_toa_process_load_go_name():
    '''
    Get the name used to title processes to load the Gene Ontology data.
    '''

    return 'Load Gene Ontology data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_gymno_01_code():
    '''
    Get the code used to identify processes to load the Gymno PLAZA 1.0 data.
    '''

    return 'toaldgymno01'

#-------------------------------------------------------------------------------

def get_toa_process_load_gymno_01_name():
    '''
    Get the name used to title processes to load the Gymno PLAZA 1.0 data.
    '''

    return 'Load Gymno PLAZA 1.0 data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_interpro_code():
    '''
    Get the code used to identify processes to load the Interpro data.
    '''

    return 'toaldinterpro'

#-------------------------------------------------------------------------------

def get_toa_process_load_interpro_name():
    '''
    Get the name used to title processes to load the Interpro data.
    '''

    return 'Load Interpro data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_load_monocots_04_code():
    '''
    Get the code used to identify processes to load the Monocots PLAZA 4.0 data.
    '''

    return 'toaldmonocots04'

#-------------------------------------------------------------------------------

def get_toa_process_load_monocots_04_name():
    '''
    Get the name used to title processes to load the Monocots PLAZA 4.0 data.
    '''

    return 'Load Monocots PLAZA 4.0 data into TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_merge_annotations_code():
    '''
    Get the code used to identify processes to merge pipeline annotations.
    '''

    return 'toamergeann'

#-------------------------------------------------------------------------------

def get_toa_process_merge_annotations_name():
    '''
    Get the name used to title processes to merge pipeline annotations.
    '''

    return 'Merge pipeline annotations'

#-------------------------------------------------------------------------------

def get_toa_process_nr_blastplus_db_code():
    '''
    Get the code of the BLAST database NR build process with BLAST+ used to identify its processes.
    '''

    return 'toabbnrbp'

#-------------------------------------------------------------------------------

def get_toa_process_nr_blastplus_db_name():
    '''
    Get the name of the BLAST database NR build process with BLAST+ used to title.
    '''

    return 'Build BLAST database NR for BLAST+'

#-------------------------------------------------------------------------------

def get_toa_process_nr_diamond_db_code():
    '''
    Get the code of the BLAST database NR build process with DIAMOND used to identify its processes.
    '''

    return 'toabbnrdn'

#-------------------------------------------------------------------------------

def get_toa_process_nr_diamond_db_name():
    '''
    Get the name of the BLAST database NR build process with DIAMOND used to title.
    '''

    return 'Build BLAST database NR for DIAMOND'

#-------------------------------------------------------------------------------

def get_toa_process_nt_blastplus_db_code():
    '''
    Get the code of the BLAST database NT build process with BLAST+ used to identify its processes.
    '''

    return 'toabbntbp'

#-------------------------------------------------------------------------------

def get_toa_process_nt_blastplus_db_name():
    '''
    Get the name of the BLAST database NT build process with BLAST+ used to title.
    '''

    return 'Build BLAST database NT for BLAST+'

#-------------------------------------------------------------------------------

def get_toa_process_pipeline_aminoacid_code():
    '''
    Get the code used to identify amino acid pipelines.
    '''

    return 'toapipelineaa'

#-------------------------------------------------------------------------------

def get_toa_process_pipeline_aminoacid_name():
    '''
    Get the name used to title amino acid pipelines.
    '''

    return 'amino acid pipeline'

#-------------------------------------------------------------------------------

def get_toa_process_pipeline_nucleotide_code():
    '''
    Get the code used to identify nucleotide pipelines.
    '''

    return 'toapipelinent'

#-------------------------------------------------------------------------------

def get_toa_process_pipeline_nucleotide_name():
    '''
    Get the name used to title nucleotide pipelines.
    '''

    return 'nucleotide pipeline'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_dicots_04_code():
    '''
    Get the code used to identify processes to build the Dicots PLAZA 4.0 proteome.
    '''

    return 'toabpdicots04'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_dicots_04_name():
    '''
    Get the name used to title processes to build the Dicots PLAZA 4.0 proteome.
    '''

    return 'Build Dicots PLAZA 4.0 proteome'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_gymno_01_code():
    '''
    Get the code used to identify processes to build the Gymno PLAZA 1.0 proteome.
    '''

    return 'toabpgymno01'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_gymno_01_name():
    '''
    Get the name to title processes to build the Gymno PLAZA 1.0 proteome.
    '''

    return 'Build Gymno PLAZA 1.0 proteome'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_monocots_04_code():
    '''
    Get the code used to identify processes to build the Monocots PLAZA 4.0 proteome.
    '''

    return 'toabpmonocots04'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_monocots_04_name():
    '''
    Get the name used to title processes to build the Monocots PLAZA 4.0 proteome.
    '''

    return 'Build Monocots PLAZA 4.0 proteome'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_refseq_plant_code():
    '''
    Get the code used to identify processes to build the NCBI RefSeq Plant proteome.
    '''

    return 'toabprefseqplt'

#-------------------------------------------------------------------------------

def get_toa_process_proteome_refseq_plant_name():
    '''
    Get the name used to title processes to build the NCBI RefSeq Plant proteome.
    '''

    return 'Build NCBI RefSeq Plant proteome'

#-------------------------------------------------------------------------------

def get_toa_process_recreate_toa_database_code():
    '''
    Get the code used to identify processes to recreate the TOA database.
    '''

    return 'toarecreatedb'

#-------------------------------------------------------------------------------

def get_get_toa_process_recreate_toa_database_name():
    '''
    Get the name used to title processes to recreate the TOA database.
    '''

    return 'Recreate TOA database'

#-------------------------------------------------------------------------------

def get_toa_process_rebuild_toa_database_code():
    '''
    Get the code used to identify processes to rebuild the TOA database.
    '''

    return 'toarebuilddb'

#-------------------------------------------------------------------------------

def get_get_toa_process_rebuild_toa_database_name():
    '''
    Get the name used to title processes to rebuild the TOA database.
    '''

    return 'Rebuild TOA database'

#-------------------------------------------------------------------------------

def get_toa_pipeline_name():
    '''
    Get the name used to title nucleotide or amino acid pipelines.
    '''

    return 'TOA nucleotide or amino acid pipeline'

#-------------------------------------------------------------------------------

def get_toa_result_dir():
    '''
    Get the result directory where results datasets are saved.
    '''

    return 'TOA-results'

#-------------------------------------------------------------------------------

def get_toa_result_database_dir():
    '''
    Get the result subdirectory where TOA process results related to the genomic database managment are saved.
    '''

    return 'TOA-databases'

#-------------------------------------------------------------------------------

def get_toa_result_installation_dir():
    '''
    Get the result subdirectory where installation process results are saved.
    '''

    return 'installations'

#-------------------------------------------------------------------------------

def get_toa_result_pipeline_dir():
    '''
    Get the result subdirectory where TOA process results related to pipelines are saved.
    '''

    return 'TOA-pipelines'

#-------------------------------------------------------------------------------

def get_toa_type_build_blastplus_db():
    '''
    Get the code used to identify processes to build BLAST databases.
    '''

    return 'build_blastplus_db'

#-------------------------------------------------------------------------------

def get_toa_type_build_diamond_db():
    '''
    Get the code used to identify processes to build DIAMOND databases.
    '''

    return 'build_diamond_db'

#-------------------------------------------------------------------------------

def get_toa_type_build_gilist():
    '''
    Get the code used to identify processes to build GeneId identifier list.
    '''

    return 'build_gilist'

#-------------------------------------------------------------------------------

def get_toa_type_build_proteome():
    '''
    Get the code used to identify processes to build proteomes.
    '''

    return 'build_proteome'

#-------------------------------------------------------------------------------

def get_toa_type_download_data():
    '''
    Get the code used to identify processes to download functional annotations from a genomic database server.
    '''

    return 'download_data'

#-------------------------------------------------------------------------------

def get_toa_type_load_data():
    '''
    Get the code used to identify processes to load data of a genomic database into TOA database.
    '''

    return 'load_data'

#-------------------------------------------------------------------------------

def get_toa_type_rebuild():
    '''
    Get the code used to identify processes to rebuild the TOA database.
    '''

    return 'rebuild'

#-------------------------------------------------------------------------------

def get_toa_type_recreate():
    '''
    Get the code used to identify processes to recreate the TOA database.
    '''

    return 'recreate'


#-------------------------------------------------------------------------------

def get_taxonomy_server():
    '''
    Get the taxonomy server URL.
    '''
    return 'https://taxonomy.jgi-psf.org/'

#-------------------------------------------------------------------------------

def get_taxonomy_dict(type, value):
    '''
    Get a taxonomy dictionary with the a species data downloaded from the taxonomy server.
    '''

    # initialize the taxonomy dictionary
    taxonomy_dict = {}

    # set the taxonomy server
    taxonomy_server = get_taxonomy_server()

    # replace spaces by underscores in value
    value = value.strip().replace(' ', '_')

    # inquire the taxonomy data to the server
    try:
        r = requests.get(f'{taxonomy_server}/{type}/{value}')
    except requests.exceptions.ConnectionError:
        raise ProgramException('W002', taxonomy_server)
    except:
        raise ProgramException('W001', taxonomy_server)

    # build the taxonomy dictionary
    if r.status_code == requests.codes.ok: #pylint: disable=no-member
        try:
            if r.json()[value].get('error','OK') == 'OK' :
                taxonomy_dict = r.json()[value]
        except:
            pass
    else:
        raise ProgramException('W003', taxonomy_server, r.status_code)

    # return taxonomy dictionary
    return taxonomy_dict

#-------------------------------------------------------------------------------

def get_alignment_tool_code_list():
    '''
    Get the code list of "alignment_tool".
    '''

    return [get_blastplus_name(), get_diamond_name()]

#-------------------------------------------------------------------------------

def get_alignment_tool_code_list_text():
    '''
    Get the code list of "alignment_tool" as text.
    '''

    return str(get_alignment_tool_code_list()).strip('[]').replace('\'', '').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_config_dir():
    '''
    Get the configuration directory in the local computer.
    '''

    return './config'

#-------------------------------------------------------------------------------

def get_keypairs_dir():
    '''
    Get the key pairs directory in the local computer.
    '''

    return './keypairs'

#-------------------------------------------------------------------------------

def get_temp_dir():
    '''
    Get the temporal directory in the local computer.
    '''

    return './temp'

#-------------------------------------------------------------------------------

def get_log_dir():
    '''
    Get the log directory in the local computer.
    '''

    return './logs'

#-------------------------------------------------------------------------------

def get_log_file(function_name=None):
    '''
    Get the log file name of in the local computer.
    '''
    # set the log file name
    now = datetime.datetime.now()
    date = datetime.datetime.strftime(now, '%y%m%d')
    time = datetime.datetime.strftime(now, '%H%M%S')
    if function_name is not None:
        log_file_name = '{0}/{1}-{2}-{3}-{4}.txt'.format(get_log_dir(), xconfiguration.environment, function_name, date, time)
    else:
        log_file_name = '{0}/{1}-x-{2}-{3}.txt'.format(get_log_dir(), xconfiguration.environment, date, time)

    # return the log file name
    return log_file_name

#-------------------------------------------------------------------------------

def list_log_files_command(local_process_id):
    '''
    Get the command to list log files in the local computer depending on the Operating System.
    '''
    # get log dir
    log_dir = get_log_dir()

    # assign the command
    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        if local_process_id == 'all':
            command = 'ls {0}/{1}-*.txt'.format(log_dir, xconfiguration.environment)
        else:
            command = 'ls {0}/{1}-{2}-*.txt'.format(log_dir, xconfiguration.environment, local_process_id)
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        log_dir = log_dir.replace('/','\\')
        if local_process_id == 'all':
            command = r'dir /B {0}\{1}-*.txt'.format(log_dir, xconfiguration.environment)
        else:
            command = r'dir /B {0}\{1}-{2}-*.txt'.format(log_dir, xconfiguration.environment, local_process_id)

    # return the command
    return command

#-------------------------------------------------------------------------------

def get_submission_process_dict():
    '''
    Get the submission process dictionary.
    '''

    # build the submission process dictionary
    submission_process_dict = {}
    submission_process_dict['add_node']= {'text': 'Add node in a cluster'}
    submission_process_dict['create_cluster']= {'text': 'Create cluster'}
    submission_process_dict['create_instance']= {'text': 'Create instance'}
    submission_process_dict['create_volume']= {'text': 'Create volume'}
    submission_process_dict['download_result_dataset']= {'text': 'Download result dataset from a cluster'}
    submission_process_dict['install_anaconda_package_list']= {'text': 'Install Anaconda package list'}
    submission_process_dict['install_ddradseqtools']= {'text': 'Install {0}'.format(get_ddradseqtools_name())}
    submission_process_dict['install_miniconda3']= {'text': 'Install {0}'.format(get_miniconda3_name())}
    submission_process_dict['install_ngshelper']= {'text': 'Install {0}'.format(get_ngshelper_name())}
    submission_process_dict['install_r']= {'text': 'Install {0}'.format(get_r_name())}
    submission_process_dict['install_raddesigner']= {'text': 'Install {0}'.format(get_raddesigner_name())}
    submission_process_dict['install_toa']= {'text': 'Install {0}'.format(get_toa_name())}
    submission_process_dict['install_transrate']= {'text': 'Install {0}'.format(get_transrate_name())}
    submission_process_dict['kill_batch_job']= {'text': 'Kill batch job'}
    submission_process_dict['link_volumes']= {'text': 'Link volumes'}
    submission_process_dict['list_clusters']= {'text': 'List clusters'}
    submission_process_dict['manage_genomic_database']= {'text': 'Manage genomic database processes'}
    submission_process_dict['manage_toa_database']= {'text': 'Manage {0} database'.format(get_toa_name())}
    submission_process_dict['manage_toa_pipeline']= {'text': 'Manage {0} pipelines'.format(get_toa_name())}
    submission_process_dict['mount_volume']= {'text': 'Mount volume in a node'}
    submission_process_dict['remove_node']= {'text': 'Remove node in a cluster'}
    submission_process_dict['remove_volume']= {'text': 'Remove volume'}
    submission_process_dict['replicate_volume']= {'text': 'Replicate volume to another zone'}
    submission_process_dict['resize_volume']= {'text': 'Resize volume'}
    submission_process_dict['restart_ddradseq_simulation_process']= {'text': 'Restart {0} process'.format(get_ddradseq_simulation_name())}
    submission_process_dict['restart_ggtrinity_process']= {'text': 'Restart {0} process'.format(get_ggtrinity_name())}
    submission_process_dict['restart_insilico_read_normalization_process']= {'text': 'Restart {0} process'.format(get_insilico_read_normalization_name())}
    submission_process_dict['restart_pipeline_process']= {'text': 'Restart {0} process'.format(get_toa_pipeline_name())}
    submission_process_dict['restart_raddesigner_process']= {'text': 'Restart {0} process'.format(get_raddesigner_name())}
    submission_process_dict['restart_soapdenovo2_process']= {'text': 'Restart {0} process'.format(get_soapdenovo2_name())}
    submission_process_dict['restart_soapdenovotrans_process']= {'text': 'Restart {0} process'.format(get_soapdenovotrans_name())}
    submission_process_dict['restart_trinity_process']= {'text': 'Restart {0} process'.format(get_trinity_name())}
    submission_process_dict['restart_variant_calling_process']= {'text': 'Restart {0} process'.format(get_variant_calling_name())}
    submission_process_dict['run_bowtie2_process']= {'text': 'Run {0} process'.format(get_bowtie2_name())}
    submission_process_dict['run_busco_process']= {'text': 'Run {0} process'.format(get_busco_name())}
    submission_process_dict['run_cd_hit_est_process']= {'text': 'Run {0} process'.format(get_cd_hit_est_name())}
    submission_process_dict['run_cuffdiff_process']= {'text': 'Run {0} process'.format(get_cuffdiff_name())}
    submission_process_dict['run_cufflinks_cuffmerge_process']= {'text': 'Run {0} process'.format(get_cufflinks_cuffmerge_name())}
    submission_process_dict['run_cuffnorm_process']= {'text': 'Run {0} process'.format(get_cuffnorm_name())}
    submission_process_dict['run_cuffquant_process']= {'text': 'Run {0} process'.format(get_cuffquant_name())}
    submission_process_dict['run_ddradseq_simulation_process']= {'text': 'Run ddRADseq simulation process'}
    submission_process_dict['run_entrez_direct_process']= {'text': 'Run {0} process'.format(get_entrez_direct_name())}
    submission_process_dict['run_express_process']= {'text': 'Run {0} process'.format(get_express_name())}
    submission_process_dict['run_fastqc_process']= {'text': 'Run {0} process'.format(get_fastqc_name())}
    submission_process_dict['run_ggtrinity_process']= {'text': 'Run {0} process'.format(get_ggtrinity_name())}
    submission_process_dict['run_gmap_process']= {'text': 'Run {0} process'.format(get_gmap_name())}
    submission_process_dict['run_gsnap_process']= {'text': 'Run {0} process'.format(get_gsnap_name())}
    submission_process_dict['run_gzip_process']= {'text': 'Run compression/decompression process'}
    submission_process_dict['run_hisat2_process']= {'text': 'Run {0} process'.format(get_hisat2_name())}
    submission_process_dict['run_htseq_count_process']= {'text': 'Run {0} process'.format(get_htseq_count_name())}
    submission_process_dict['run_insilico_read_normalization_process']= {'text': 'Run {0} process'.format(get_insilico_read_normalization_name())}
    submission_process_dict['run_ipyrad_process']= {'text': 'Run {0} pipeline process'.format(get_ipyrad_name())}
    submission_process_dict['run_kallisto_process']= {'text': 'Run {0} process'.format(get_kallisto_name())}
    submission_process_dict['run_pipeline_process']= {'text': 'Run {0} process'.format(get_toa_pipeline_name())}
    submission_process_dict['run_quast_process']= {'text': 'Run {0} process'.format(get_quast_name())}
    submission_process_dict['run_raddesigner_process']= {'text': 'Run {0} process'.format(get_raddesigner_name())}
    submission_process_dict['run_ref_eval_process']= {'text': 'Run {0} process'.format(get_ref_eval_name())}
    submission_process_dict['run_rnaquast_process']= {'text': 'Run {0} process'.format(get_rnaquast_name())}
    submission_process_dict['run_rsem_eval_process']= {'text': 'Run {0} process'.format(get_rsem_eval_name())}
    submission_process_dict['run_rsitesearch_process']= {'text': 'Run {0} process'.format(get_rsitesearch_name())}
    submission_process_dict['run_soapdenovo2_process']= {'text': 'Run {0} process'.format(get_soapdenovo2_name())}
    submission_process_dict['run_soapdenovotrans_process']= {'text': 'Run {0} process'.format(get_soapdenovotrans_name())}
    submission_process_dict['run_star_process']= {'text': 'Run {0} process'.format(get_star_name())}
    submission_process_dict['run_starcode_process']= {'text': 'Run {0} process'.format(get_starcode_name())}
    submission_process_dict['run_tophat_process']= {'text': 'Run {0} process'.format(get_tophat_name())}
    submission_process_dict['run_transabyss_process']= {'text': 'Run {0} process'.format(get_transabyss_name())}
    submission_process_dict['run_transcript_filter_process']= {'text': 'Run {0} process'.format(get_transcript_filter_name())}
    submission_process_dict['run_transcriptome_blastx_process']= {'text': 'Run {0} process'.format(get_transcriptome_blastx_name())}
    submission_process_dict['run_transdecoder_process']= {'text': 'Run {0} process'.format(get_transdecoder_name())}
    submission_process_dict['run_transrate_process']= {'text': 'Run {0} process'.format(get_transrate_name())}
    submission_process_dict['run_trimmomatic_process']= {'text': 'Run {0} process'.format(get_trimmomatic_name())}
    submission_process_dict['run_trinity_process']= {'text': 'Run {0} process'.format(get_trinity_name())}
    submission_process_dict['run_variant_calling_process']= {'text': 'Run {0} process'.format(get_variant_calling_name())}
    submission_process_dict['show_cluster_composition']= {'text': 'Show cluster composition'}
    submission_process_dict['show_status_batch_jobs']= {'text': 'Show status of batch jobs'}
    submission_process_dict['terminate_cluster']= {'text': 'Terminate cluster'}
    submission_process_dict['terminate_instance']= {'text': 'Terminate instance'}
    submission_process_dict['terminate_volume_creator']= {'text': 'Terminate volume creator'}
    submission_process_dict['unmount_volume']= {'text': 'Unmount volume in a node'}
    submission_process_dict['upload_database_dataset']= {'text': 'Upload database dataset to a cluster'}
    submission_process_dict['upload_read_dataset']= {'text': 'Upload read dataset to a cluster'}
    submission_process_dict['upload_reference_dataset']= {'text': 'Upload reference dataset to a cluster'}

    # return the submission process dictionary
    return submission_process_dict

#-------------------------------------------------------------------------------

def get_submission_process_id(submission_process_text):
    '''
    Get the submission process identification from the submission process text.
    '''

    # initialize the control variable
    submission_process_id_found = None

    # get the dictionary of the submission processes
    submission_process_dict = get_submission_process_dict()

    # search the submission process identification
    for submission_process_id in submission_process_dict.keys():
        if submission_process_dict[submission_process_id]['text'] == submission_process_text:
            submission_process_id_found = submission_process_id
            break

    # return the submission process identification
    return submission_process_id_found
    # return the submission process identification

#-------------------------------------------------------------------------------

def get_cluster_ngscloud_dir():
    '''
    Get the NGScloud dataset directory in the cluster (single-volume).
    '''

    return '/ngscloud2'

#-------------------------------------------------------------------------------

def get_cluster_ngscloud_device_file():
    '''
    Get the device file used by the NGScloud dataset volume in AWS console (single-volume).
    '''

    return '/dev/sdf'

#-------------------------------------------------------------------------------

def get_cluster_app_dir():
    '''
    Get the aplication directory in the cluster.
    '''

    return f'{get_cluster_ngscloud_dir()}/apps'

#-------------------------------------------------------------------------------

def get_cluster_app_device_file():
    '''
    Get the device file used by the application volume in AWS console (multi-volume).
    '''

    return '/dev/sdg'

#-------------------------------------------------------------------------------

def get_cluster_database_dir():
    '''
    Get the database directory in the cluster.
    '''

    return f'{get_cluster_ngscloud_dir()}/databases'

#-------------------------------------------------------------------------------

def get_cluster_database_device_file():
    '''
    Get the device file used by the database volume in AWS console (multi-volume).
    '''

    return '/dev/sdh'

#-------------------------------------------------------------------------------

def get_cluster_database_dataset_dir(database_dataset_id):
    '''
    Get the directory of a database dataset in the cluster.
    '''

    # set the database directory in the cluster
    cluster_database_dataset_dir = '{0}/{1}'.format(get_cluster_database_dir(), database_dataset_id)

    # return the database directory in the cluster
    return cluster_database_dataset_dir

#-------------------------------------------------------------------------------

def get_cluster_database_file(database_dataset_id, file_name):
    '''
    Get the database file path of a database dataset in the cluster.
    '''

    # set the path of the database file
    cluster_database_file = '{0}/{1}'.format(get_cluster_database_dataset_dir(database_dataset_id), os.path.basename(file_name))

    # return the path of the database file
    return cluster_database_file

#-------------------------------------------------------------------------------

def get_cluster_read_dir():
    '''
    Get the read directory in the cluster.
    '''

    return f'{get_cluster_ngscloud_dir()}/reads'

#-------------------------------------------------------------------------------

def get_cluster_read_device_file():
    '''
    Get the device file used by the read volume in AWS console (multi-volume).
    '''

    return '/dev/sdi'

#-------------------------------------------------------------------------------

def get_uploaded_read_dataset_name():
    '''
    Get the name of the raw read dataset in the cluster.
    '''

    return 'uploaded-reads'

#-------------------------------------------------------------------------------

def get_cluster_experiment_read_dataset_dir(experiment_id, read_dataset_id):
    '''
    Get the directory of a experiment read dataset in the cluster.
    '''

    # set the experiment read directory in the cluster
    cluster_experiment_read_dataset_dir = '{0}/{1}/{2}'.format(get_cluster_read_dir(), experiment_id, read_dataset_id)

    # return the experiment read directory in the cluster
    return cluster_experiment_read_dataset_dir

#-------------------------------------------------------------------------------

def get_cluster_read_file(experiment_id, read_dataset_id, file_name):
    '''
    Get the read file path of an experiment read dataset in the cluster.
    '''

    # set the path of the read file
    cluster_read_file = '{0}/{1}'.format(get_cluster_experiment_read_dataset_dir(experiment_id, read_dataset_id), os.path.basename(file_name))

    # return the path of the read file
    return cluster_read_file

#-------------------------------------------------------------------------------

def get_cluster_reference_dir():
    '''
    Get the reference directory in the cluster.
    '''

    return f'{get_cluster_ngscloud_dir()}/references'

#-------------------------------------------------------------------------------

def get_cluster_reference_device_file():
    '''
    Get the device file used by the reference volume in AWS console (multi-volume).
    '''

    return '/dev/sdj'

#-------------------------------------------------------------------------------

def get_cluster_reference_dataset_dir(reference_dataset_id):
    '''
    Get the directory of a reference dataset in the cluster.
    '''

    # set the reference directory in the cluster
    cluster_reference_dataset_dir = '{0}/{1}'.format(get_cluster_reference_dir(), reference_dataset_id)

    # return the reference directory in the cluster
    return cluster_reference_dataset_dir

#-------------------------------------------------------------------------------

def get_cluster_reference_file(reference_dataset_id, file_name):
    '''
    Get the reference file path of a reference dataset in the cluster.
    '''

    # set the path of the reference file
    cluster_reference_file = '{0}/{1}'.format(get_cluster_reference_dataset_dir(reference_dataset_id), os.path.basename(file_name))

    # return the path of the reference file
    return cluster_reference_file

#-------------------------------------------------------------------------------

def get_cluster_result_dir():
    '''
    Get the result directory in the cluster.
    '''

    return f'{get_cluster_ngscloud_dir()}/results'

#-------------------------------------------------------------------------------

def get_cluster_result_device_file():
    '''
    Get the device file used by the result volume in AWS console (multi-volume).
    '''

    return '/dev/sdk'

#-------------------------------------------------------------------------------

def get_design_dataset_name():
    '''
    Get the name of design dataset in the cluster.
    '''

    return 'designs'

#-------------------------------------------------------------------------------

def get_cluster_experiment_result_dir(experiment_id):
    '''
    Get the directory of run result datasets in the cluster.
    '''

    # set the run result directory in the cluster
    cluster_experiment_results_dir = '{0}/{1}'.format(get_cluster_result_dir(), experiment_id)

    # return the run result directory in the cluster
    return cluster_experiment_results_dir

#-------------------------------------------------------------------------------

def get_cluster_experiment_result_dataset_dir(experiment_id, result_dataset_id):
    '''
    Get the directory of an experiment result dataset in the cluster.
    '''

    # set the experiment result dataset directory in the cluster
    cluster_experiment_result_dataset_dir = '{0}/{1}/{2}'.format(get_cluster_result_dir(), experiment_id, result_dataset_id)

    # return the experiment result dataset directory in the cluster
    return cluster_experiment_result_dataset_dir

#-------------------------------------------------------------------------------

def get_cluster_result_file(experiment_id, result_dataset_id, file_name):
    '''
    Get the result file path of an experiment read dataset in the cluster.
    '''

    # set the path of the read file
    cluster_result_file = '{0}/{1}'.format(get_cluster_experiment_result_dataset_dir(experiment_id, result_dataset_id), os.path.basename(file_name))

    # return the path of the read file
    return cluster_result_file

#-------------------------------------------------------------------------------

def get_cluster_current_run_dir(experiment_id, process):
    '''
    Get the run directory of a process in the cluster.
    '''

    # set the run identificacion
    now = datetime.datetime.now()
    date = datetime.datetime.strftime(now, '%y%m%d')
    time = datetime.datetime.strftime(now, '%H%M%S')
    run_id = '{0}-{1}-{2}'.format(process, date, time)

    # set the run directory in the cluster
    cluster_current_run_dir = get_cluster_experiment_result_dir(experiment_id) + '/' + run_id

    # return the run directory in the cluster
    return cluster_current_run_dir

#-------------------------------------------------------------------------------

def get_status_dir(current_run_dir):
    '''
    Get the status directory of a process in the cluster.
    '''

    # set the status directory
    status_dir = '{0}/status'.format(current_run_dir)

    # return the status directory
    return status_dir

#-------------------------------------------------------------------------------

def get_status_ok(current_run_dir):
    '''
    Get the OK status file in the cluster.
    '''

    # set the OK status file
    ok_status = '{0}/status/script.ok'.format(current_run_dir)

    # return the OK status file
    return ok_status

#-------------------------------------------------------------------------------

def get_status_wrong(current_run_dir):
    '''
    Get the WRONG status file in the cluster.
    '''

    # set the WRONG status file
    wrong_status = '{0}/status/script.wrong'.format(current_run_dir)

    # return the WRONG status file
    return wrong_status

#-------------------------------------------------------------------------------

def get_mounting_point_list():
    '''
    Get the available mounting point list
    '''

    return [get_cluster_app_dir(), get_cluster_database_dir(), get_cluster_read_dir(), get_cluster_reference_dir(), get_cluster_result_dir()]

#-------------------------------------------------------------------------------

def get_cluster_log_file():
    '''
    Get the log file name of an experiment run in the cluster.
    '''

    return 'log.txt'

#-------------------------------------------------------------------------------

def change_extension(path, new_extension):
    '''Change the file extension.'''

    # get the path with the new extension
    i = path.rfind('.')
    if i >= 0:
        new_path = path[:i + 1] + new_extension
    else:
        new_path = path + new_extension

    # return the path with new extension
    return new_path

#-------------------------------------------------------------------------------

def is_valid_path(path, operating_system=sys.platform):
    '''
    Check if a path is a valid path.
    '''

    # initialize control variable
    valid = False

    # check if the path is valid
    if operating_system.startswith('linux') or operating_system.startswith('darwin'):
        # -- valid = re.match('^(/.+)(/.+)*/?$', path)
        valid = True
    elif operating_system.startswith('win32') or operating_system.startswith('cygwin'):
        valid = True

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def is_absolute_path(path, operating_system=sys.platform):
    '''
    Check if a path is a absolute path.
    '''

    # initialize control variable
    valid = False

    # check if the path is absolute
    if operating_system.startswith('linux') or operating_system.startswith('darwin'):
        if path != '':
            # -- valid = is_path_valid(path) and path[0] == '/'
            valid = True
    elif operating_system.startswith('win32') or operating_system.startswith('cygwin'):
        valid = True

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def is_relative_path(path, operating_system=sys.platform):
    '''
    Check if a path is a relative path.
    '''

    # initialize control variable
    valid = False

    # check if the path is valid
    if operating_system.startswith('linux') or operating_system.startswith('darwin'):
        valid = True
    elif operating_system.startswith('win32') or operating_system.startswith('cygwin'):
        valid = True

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def is_device_file(path, device_pattern):
    '''
    Check if a path is a valid device file, e.g. /dev/sdk.
    '''

    # initialize control variable
    valid = False

    # build the complete pattern
    pattern = '^{0}$'.format(device_pattern)

    # check if path is a valid device file
    valid = re.match(pattern, path)

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def get_machine_device_file(aws_device_file):
    '''
    Get de machine device file from AWS device file, e.g. /dev/sdb1 -> /dev/xvdb1.
    sudo lsblk --output NAME,TYPE,SIZE,FSTYPE,MOUNTPOINT,LABEL
    sudo blkid
    '''

    # determine the machine device file
    machine_device_file = aws_device_file[0:5] + 'xv' + aws_device_file[6:]

    # return the machine device file
    return machine_device_file

#-------------------------------------------------------------------------------

def is_name_valid(name):
    '''
    Check if a name is valid.
    '''

    # initialize control variable
    valid = True

    # check if name is valid
    for i in range(len(name)):
        if not name[i].isalnum() and name[i] not in ['_']:
            valid = False
            break

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def is_email_address_valid(email):
    '''
    Check if an e-mail address is valid.
    '''

    # initialize control variable
    valid = False

    # build the complete pattern
    pattern = r'^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,4})$'

    # check if the e-mail address is valid
    valid = re.match(pattern, email)

    # return control variable
    return valid

#-------------------------------------------------------------------------------

def get_option_dict(config_file):
    '''
    Get a dictionary with the options retrieved from a configuration file.
    '''

    # initialize the option dictionary
    option_dict = {}

    # create class to parse the configuration files
    config = configparser.ConfigParser()

    # read the configuration file
    config.read(config_file)

    # build the dictionary
    for section in config.sections():
        # get the keys dictionary
        keys_dict = option_dict.get(section, {})
        # for each key in the section
        for key in config[section]:
            # get the value of the key
            value = config.get(section, key, fallback='')
            # add a new enter in the keys dictionary
            keys_dict[key] = get_option_value(value)
        # update the section with its keys dictionary
        option_dict[section] = keys_dict

    # return the option dictionary
    return option_dict

#-------------------------------------------------------------------------------

def get_option_value(option):
    '''
    Remove comments ans spaces from an option retrieve from a configuration file.
    '''

    # Remove comments
    position = option.find('#')
    if position == -1:
        value = option
    else:
        value = option[:position]

    # Remove comments
    value = value.strip()

    # return the value without comments and spaces
    return value

#-------------------------------------------------------------------------------

def split_literal_to_integer_list(literal):
    '''
    Split a string literal in a integer value list which are separated by comma.
    '''
  
    # initialize the string values list and the interger values list
    strings_list = []
    integers_list = []
    
    # split the string literal in a string values list
    strings_list = split_literal_to_string_list(literal)

    # convert each value from string to integer
    for i in range(len(strings_list)):
        try:
            integers_list.append(int(strings_list[i]))
        except:
            integers_list = []
            break

    # return the integer values list
    return integers_list

#-------------------------------------------------------------------------------

def split_literal_to_float_list(literal):
    '''
    Split a string literal in a float value list which are separated by comma.
    '''
  
    # initialize the string values list and the float values list
    strings_list = []
    float_list = []
    
    # split the string literal in a string values list
    strings_list = split_literal_to_string_list(literal)

    # convert each value from string to float
    for i in range(len(strings_list)):
        try:
            float_list.append(float(strings_list[i]))
        except:
            float_list = []
            break

    # return the float values list
    return float_list

#-------------------------------------------------------------------------------

def split_literal_to_string_list(literal):
    '''
    Split a string literal in a string value list which are separated by comma.
    '''
  
    # initialize the string values list
    string_list = []

    # split the string literal in a string values list
    string_list = literal.split(',')

    # remove the leading and trailing whitespaces in each value
    for i in range(len(string_list)):
        string_list[i] = string_list[i].strip()

    # return the string values list
    return string_list

#-------------------------------------------------------------------------------

def pair_files(file_name_list, specific_chars_1, specific_chars_2):
    '''
    ...
    '''

    # initialize the file lists
    file_name_1_list = []
    file_name_2_list = []
    unpaired_file_name_list = []

    # for each file name, append it to the corresponding list
    for file_name in file_name_list:
        if file_name.find(specific_chars_1) >= 0:
            file_name_1_list.append(file_name)
        elif  file_name.find(specific_chars_2) >= 0:
            file_name_2_list.append(file_name)
        else:
            unpaired_file_name_list.append(file_name)
    file_name_1_list.sort()
    file_name_2_list.sort()

    # check the file pairing
    short_file_name_1 = ''
    short_file_name_2 = ''
    review_file_name_1_list = []
    review_file_name_2_list = []
    index_1 = 0
    index_2 = 0
    while index_1 < len(file_name_1_list) or index_2 < len(file_name_2_list):
        if index_1 < len(file_name_1_list):
            file_name_1 = file_name_1_list[index_1]
            short_file_name_1 = file_name_1.replace(specific_chars_1, '')
        if index_2 < len(file_name_2_list):
            file_name_2 = file_name_2_list[index_2]
            short_file_name_2 = file_name_2.replace(specific_chars_2, '')
        if short_file_name_1 == short_file_name_2:
            review_file_name_1_list.append(file_name_1)
            index_1 += 1
            review_file_name_2_list.append(file_name_2)
            index_2 += 1
        elif short_file_name_1 < short_file_name_2:
            unpaired_file_name_list.append(file_name_1)
            index_1 += 1
        elif short_file_name_1 > short_file_name_2:
            unpaired_file_name_list.append(file_name_2)
            index_2 += 1

    # return the file lists
    return (review_file_name_1_list, review_file_name_2_list, unpaired_file_name_list)

#-------------------------------------------------------------------------------
    
def get_fasta_merger_operation_code_list():
    '''
    Get the code list of "merger_operation" with FASTA.
    '''

    return ['1AND2', '1LESS2']

#-------------------------------------------------------------------------------
    
def get_fasta_merger_operation_code_list_text():
    '''
    Get the code list of "merger_operation" with FASTA files as text.
    '''

    return '1AND2 (sequences included in both files) or 1LESS2 (sequences in 1 and not in 2)'

#-------------------------------------------------------------------------------
    
def get_annotation_merger_operation_code_list():
    '''
    Get the code list of "merger_operation" with annotation files.
    '''

    return ['1AND2', '1BEST']

#-------------------------------------------------------------------------------
    
def get_annotation_merger_operation_code_list_text():
    '''
    Get the code list of "merger_operation" with annotation files as text.
    '''

    return '1AND2 (annotations included in both files) or 1BEST (all annotations of the first file and annotations of the second file if their seq id is not in the first)'

#-------------------------------------------------------------------------------

def check_startswith(literal, text_list, case_sensitive=False):
    '''
    Check if a literal starts with a text in a list.
    '''

    # initialize the control variable
    OK = False
  
    # initialize the working list
    list = []

    # if the codification is not case sensitive, convert the code and code list to uppercase
    if not case_sensitive:
        try:
            literal = literal.upper()
        except:
            pass
        try:
            list = [x.upper() for x in text_list]
        except:
            pass
    else:
        list = text_list

    # check if the literal starts with a text in the list
    for text in list:
        if literal.startswith(text):
            OK = True
            break

    # return control variable
    return OK

#-------------------------------------------------------------------------------

def check_code(literal, code_list, case_sensitive=False):
    '''
    Check if a literal is in a code list.
    '''
    
    # initialize the working list
    list = []

    # if the codification is not case sensitive, convert the code and code list to uppercase
    if not case_sensitive:
        try:
            literal = literal.upper()
        except:
            pass
        try:
            list = [x.upper() for x in code_list]
        except:
            pass
    else:
        list = code_list

    # check if the code is in the code list
    OK = literal in list

    # return control variable
    return OK

#-------------------------------------------------------------------------------

def check_int(literal, minimum=(-sys.maxsize - 1), maximum=sys.maxsize):
    '''
    Check if a numeric or string literal is an integer number.
    '''

    # initialize the control variable
    OK = True
  
    # check the number
    try:
        int(literal)
        int(minimum)
        int(maximum)
    except:
        OK = False
    else:
        if int(literal) < int(minimum) or int(literal) > int(maximum):
            OK = False

    # return control variable
    return OK

#-------------------------------------------------------------------------------

def check_float(literal, minimum=float(-sys.maxsize - 1), maximum=float(sys.maxsize), mne=0.0, mxe=0.0):
    '''
    Check if a numeric or string literal is a float number.
    '''

    # initialize the control variable
    OK = True
  
    # check the number
    try:
        float(literal)
        float(minimum)
        float(maximum)
        float(mne)
        float(mxe)
    except:
        OK = False
    else:
        if float(literal) < (float(minimum) + float(mne)) or float(literal) > (float(maximum) - float(mxe)):
            OK = False

    # return control variable
    return OK

#-------------------------------------------------------------------------------

def check_parameter_list(parameters, key, not_allowed_parameters_list):
    '''
    Check if a string contains a parameter list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the parameter list
    parameter_list = [x.strip() for x in parameters.split(';')]

    # check the parameter list
    for parameter in parameter_list:
        try:
            if parameter.find('=') > 0:
                pattern = r'^--(.+)=(.+)$'
                mo = re.search(pattern, parameter)
                parameter_name = mo.group(1).strip()
                # -- parameter_value = mo.group(2).strip()
            else:
                pattern = r'^--(.+)$'
                mo = re.search(pattern, parameter)
                parameter_name = mo.group(1).strip()
        except:
            error_list.append('*** ERROR: the value of the key "{0}" has to NONE or a valid parameter list.'.format(key))
            OK = False
            break
        if parameter_name in not_allowed_parameters_list:
            error_list.append('*** ERROR: the parameter {0} is not allowed in the key "{1}" because it is controled by {2}.'.format(parameter_name, key, get_project_name()))
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_command(command, log):
    '''
    Run a Bash shell command and redirect stdout and stderr to log.
    '''

    # run the command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    for line in iter(process.stdout.readline, b''):
        # replace non-ASCII caracters by one blank space
        line = re.sub(b'[^\x00-\x7F]+', b' ', line)
        # control return code and new line characters
        if not isinstance(log, DevStdOut):
            line = re.sub(b'\r\n', b'\r', line)
            line = re.sub(b'\r', b'\r\n', line)
        elif sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            pass
        elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
            line = re.sub(b'\r\n', b'\r', line)
            line = re.sub(b'\r', b'\r\n', line)
        # create a string from the bytes literal
        line = line.decode('utf-8')
        # write the line in log
        log.write('{0}'.format(line))
    rc = process.wait()

    # return the return code of the command run
    return rc

#-------------------------------------------------------------------------------

def get_nucleotide_dict():
    '''
    Get a dictionary with nucleotide data.
    '''

    # +----+------------------+------------------+-------------+
    # |Code|   Description    |   Translation    |Complementary|
    # +----+------------------+------------------+-------------+
    # | A  |Adenine           |A                 |     T/U     |
    # +----+------------------+------------------+-------------+
    # | C  |Cytosine          |C                 |      G      |
    # +----+------------------+------------------+-------------+
    # | G  |Guanine           |G                 |      C      |
    # +----+------------------+------------------+-------------+
    # | T  |Thymine           |T                 |      A      |
    # +----+------------------+------------------+-------------+
    # | U  |Uracil            |U                 |      A      |
    # +----+------------------+------------------+-------------+
    # | R  |puRine            |A or G            |      Y      |
    # +----+------------------+------------------+-------------+
    # | Y  |pYrimidine        |C or T/U          |      R      |
    # +----+------------------+------------------+-------------+
    # | S  |Strong interaction|C or G            |      S      |
    # +----+------------------+------------------+-------------+
    # | W  |Weak interaction  |A or T/U          |      W      |
    # +----+------------------+------------------+-------------+
    # | K  |Keto group        |G or T/U          |      M      |
    # +----+------------------+------------------+-------------+
    # | M  |aMino group       |A or C            |      K      |
    # +----+------------------+------------------+-------------+
    # | B  |not A             |C or G or T/U     |      V      |
    # +----+------------------+------------------+-------------+
    # | V  |not T             |A or C or G       |      B      |
    # +----+------------------+------------------+-------------+
    # | D  |not C             |A or G or T/U     |      H      |
    # +----+------------------+------------------+-------------+
    # | H  |not G             |A or C or T/U     |      D      |
    # +----+------------------+------------------+-------------+
    # | N  |aNy               |A or C or G or T/U|      N      |
    # +----+------------------+------------------+-------------+

    # build the nucleotide dictonary
    nucleotide_dict = {
        'A':{'code': 'A', 'nuclotide_list':['A'], 'complementary_code':'T', 'complementary_nuclotide_list':['T']},
        'a':{'code': 'a', 'nuclotide_list':['a'], 'complementary_code':'t', 'complementary_nuclotide_list':['t']},
        'C':{'code': 'C', 'nuclotide_list':['C'], 'complementary_code':'G', 'complementary_nuclotide_list':['G']},
        'c':{'code': 'c', 'nuclotide_list':['c'], 'complementary_code':'g', 'complementary_nuclotide_list':['g']},
        'G':{'code': 'G', 'nuclotide_list':['G'], 'complementary_code':'C', 'complementary_nuclotide_list':['C']},
        'g':{'code': 'g', 'nuclotide_list':['g'], 'complementary_code':'c', 'complementary_nuclotide_list':['c']},
        'T':{'code': 'T', 'nuclotide_list':['T'], 'complementary_code':'A', 'complementary_nuclotide_list':['A']},
        't':{'code': 't', 'nuclotide_list':['t'], 'complementary_code':'a', 'complementary_nuclotide_list':['a']},
        'R':{'code': 'R', 'nuclotide_list':['A','G'], 'complementary_code':'Y', 'complementary_nuclotide_list':['C','T']},
        'r':{'code': 'r', 'nuclotide_list':['a','g'], 'complementary_code':'y', 'complementary_nuclotide_list':['c','t']},
        'Y':{'code': 'Y', 'nuclotide_list':['C','T'], 'complementary_code':'R', 'complementary_nuclotide_list':['A','G']},
        'y':{'code': 'y', 'nuclotide_list':['c','t'], 'complementary_code':'r', 'complementary_nuclotide_list':['a','g']},
        'S':{'code': 'S', 'nuclotide_list':['C','G'], 'complementary_code':'S', 'complementary_nuclotide_list':['C','G']},
        's':{'code': 's', 'nuclotide_list':['c','G'], 'complementary_code':'s', 'complementary_nuclotide_list':['c','g']},
        'W':{'code': 'W', 'nuclotide_list':['A','T'], 'complementary_code':'W', 'complementary_nuclotide_list':['A','T']},
        'w':{'code': 'w', 'nuclotide_list':['a','t'], 'complementary_code':'w', 'complementary_nuclotide_list':['a','t']},
        'K':{'code': 'K', 'nuclotide_list':['G','T'], 'complementary_code':'M', 'complementary_nuclotide_list':['A','C']},
        'k':{'code': 'k', 'nuclotide_list':['g','t'], 'complementary_code':'m', 'complementary_nuclotide_list':['a','c']},
        'M':{'code': 'M', 'nuclotide_list':['A','C'], 'complementary_code':'K', 'complementary_nuclotide_list':['G','T']},
        'm':{'code': 'm', 'nuclotide_list':['a','c'], 'complementary_code':'k', 'complementary_nuclotide_list':['g','t']},
        'B':{'code': 'B', 'nuclotide_list':['C','G','T'], 'complementary_code':'V', 'complementary_nuclotide_list':['A','C','G']},
        'b':{'code': 'b', 'nuclotide_list':['c','G','T'], 'complementary_code':'v', 'complementary_nuclotide_list':['a','c','g']},
        'V':{'code': 'V', 'nuclotide_list':['A','C','G'], 'complementary_code':'B', 'complementary_nuclotide_list':['C','G','T']},
        'v':{'code': 'v', 'nuclotide_list':['a','c','g'], 'complementary_code':'b', 'complementary_nuclotide_list':['c','g','t']},
        'D':{'code': 'D', 'nuclotide_list':['A','G','T'], 'complementary_code':'H', 'complementary_nuclotide_list':['A','C','T']},
        'd':{'code': 'd', 'nuclotide_list':['a','g','t'], 'complementary_code':'h', 'complementary_nuclotide_list':['a','c','t']},
        'H':{'code': 'H', 'nuclotide_list':['A','C','T'], 'complementary_code':'D', 'complementary_nuclotide_list':['A','G','T']},
        'h':{'code': 'h', 'nuclotide_list':['a','C','t'], 'complementary_code':'d', 'complementary_nuclotide_list':['a','g','t']},
        'N':{'code': 'N', 'nuclotide_list':['A','C','G','T'], 'complementary_code':'N', 'complementary_nuclotide_list':['A','C','G','T']},
        'n':{'code': 'n', 'nuclotide_list':['a','c','g','t'], 'complementary_code':'n', 'complementary_nuclotide_list':['a','c','g','t']}
        }
    
    # return the nucleotide dictionary
    return nucleotide_dict

#-------------------------------------------------------------------------------

def get_nucleotide_list(allowed_ambiguity_codes, allowed_lowercase_code):
    '''
    Get a list with the nucleotide codes.
    '''

    # initialize the nucleotide list
    nucleotide_list = []

    # get the nucleotide dictonary
    nucleotide_dict = get_nucleotide_dict()

    # build the nucleotide list
    for code in nucleotide_dict.keys():
        lenght = len(nucleotide_dict[code]['nuclotide_list'])
        if (not allowed_ambiguity_codes and lenght == 1 or allowed_ambiguity_codes) and (code.isupper() or code.islower() and allowed_lowercase_code):
            nucleotide_list.append(code)

    # sort the nucleotide_list
    if nucleotide_list != []:
        nucleotide_list.sort()

    # return the nucleotide list
    return nucleotide_list

#-------------------------------------------------------------------------------

def is_valid_sequence(seq, allowed_ambiguity_codes, other_allowed_characters_list, cut_tag_check):
    '''
    Check if seq have a valid nucleotides sequence. In addition to standard codes,
    others allowed characters can be passed.
    '''

    # initialize the control variable
    OK = True

    # get nucleotide list
    nucleotide_list = get_nucleotide_list(allowed_ambiguity_codes, allowed_lowercase_code=True)

    # set cut tag and cut tag counter
    cut_tag = '*'
    cut_tag_counter = 0

    # check each nucleotide of the sequence
    for i in range(len(seq)):
        if cut_tag_check:
            if seq[i] not in nucleotide_list and seq[i] != cut_tag and seq[i] not in other_allowed_characters_list:
                OK = False
                break
            if seq[i] == cut_tag:
                cut_tag_counter += 1
        else:
            if seq[i] not in nucleotide_list and seq[i] not in other_allowed_characters_list:
                OK = False
                break

    # check the cut tag counter
    if cut_tag_check:
        if cut_tag_counter != 1:
            OK = False

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def get_complementary_sequence(seq):
    '''
    Get the complementary sequence of seq.
    '''

    # get the nucleotide dictionary
    nucleotide_dict =  get_nucleotide_dict()

    # convert the sequence to a list
    seq_list = list(seq)

    # get the list changing each nucleotide by its complementary nucleotide
    complementary_seq_list = [nucleotide_dict[nucleotide]['complementary_code'] for nucleotide in seq_list]

    # get a string from the complementary list 
    complementary_seq = ''.join(complementary_seq_list)

    # return the complementary sequence
    return complementary_seq

#-------------------------------------------------------------------------------

def get_reverse_sequence(seq):
    '''
    Get the reverse sequence of seq.
    '''

    # convert the sequence to a list and reverse the elements of the list
    seq_list = list(seq)
    seq_list.reverse()

    # get a string from the reverse list 
    reverse_seq = ''.join(seq_list)

    # return the reverse complementary sequence
    return reverse_seq

#-------------------------------------------------------------------------------

def get_reverse_complementary_sequence(seq):
    '''
    Get the reverse complementary sequence of seq.
    '''

    # get the nucleotide dictionary
    nucleotide_dict =  get_nucleotide_dict()

    # convert the sequence to a list and reverse the elements of the list
    seq_list = list(seq)
    seq_list.reverse()

    # get the reverse list changing each nucleotide by its complementary nucleotide
    revcompl_seq_list = [nucleotide_dict[nucleotide]['complementary_code'] for nucleotide in seq_list]

    # get a string from the reverse complementary list 
    revcompl_seq = ''.join(revcompl_seq_list)

    # return the reverse complementary sequence
    return revcompl_seq

#-------------------------------------------------------------------------------

def get_na():
    '''
    Get the characters to represent not available.
    '''

    return 'N/A'

#-------------------------------------------------------------------------------

def get_separator():
    '''
    Get the separation line between process steps.
    '''

    return '**************************************************'

#-------------------------------------------------------------------------------

def get_time_output_format(separator=True):
    '''
    Get the format of the command time.
    '''

    # set the format
    format = 'Elapsed real time (s): %e\\n' + \
             'CPU time in kernel mode (s): %S\\n' + \
             'CPU time in user mode (s): %U\\n' + \
             'Percentage of CPU: %P\\n' + \
             'Maximum resident set size(Kb): %M\\n' + \
             'Average total memory use (Kb):%K'
    if separator:
       format = '$SEP\\n' + format

    # return the format
    return format

#-------------------------------------------------------------------------------

def get_mail_message_ok(process_name, cluster_name):
    '''
    Get the message text of the mail sent when a process ends with errors.
    '''

    # build the message
    message = f'The {process_name} process ended OK in node $HOST_ADDRESS of cluster {cluster_name} ' + \
               'at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION). ' + \
               'Please review its log.<br/>' + \
               '<br/>' + \
               'Regards,<br/>' + \
               'GI Sistemas Naturales e Historia Forestal<br/>'+ \
               '(formerly known as GI Genetica, Fisiologia e Historia Forestal)<br/>' + \
               'Dpto. Sistemas y Recursos Naturales<br/>' + \
               'ETSI Montes, Forestal y del Medio Natural<br/>' + \
               'Universidad Politecnica de Madrid<br/>' + \
               'https://github.com/ggfhf/'


    # return the message
    return message

#-------------------------------------------------------------------------------

def get_mail_message_wrong(process_name, cluster_name):
    '''
    Get the message text of the mail sent when a process ends with errors.
    '''

    # build the message
    message = f'The {process_name} process ended *** WRONG *** in node $HOST_ADDRESS of cluster {cluster_name} ' + \
               'at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION). ' + \
               'Please review its log.<br/>' + \
               '<br/>' + \
               'Regards,<br/>' + \
               'GI Sistemas Naturales e Historia Forestal<br/>'+ \
               '(formerly known as GI Genetica, Fisiologia e Historia Forestal)<br/>' + \
               'Dpto. Sistemas y Recursos Naturales<br/>' + \
               'ETSI Montes, Forestal y del Medio Natural<br/>' + \
               'Universidad Politecnica de Madrid<br/>' + \
               'https://github.com/ggfhf/'


    # build the message
    return message

#-------------------------------------------------------------------------------

class DevStdOut(object):
    '''
    This class is used when it is necessary write in sys.stdout and in a log file
    '''

    #---------------

    def __init__(self, calling_function=None, print_stdout=True):
        '''
        Execute actions correspending to the creation of a "DevStdOut" instance.
        '''

        # save initial parameters in instance variables
        self.calling_function = calling_function
        self.print_stdout = print_stdout

        # get the local log file
        self.log_file = get_log_file(self.calling_function)

        # open the local log file
        try:
            if not os.path.exists(os.path.dirname(self.log_file)):
                os.makedirs(os.path.dirname(self.log_file))
            self.log_file_id = open(self.log_file, mode='w', encoding='iso-8859-1', newline='\n')
        except:
            print('*** ERROR: The file {0} can not be created'.format(self.log_file))

    #---------------

    def write(self, message):
        '''
        Write the message in sys.stadout and in the log file
        '''

        # write in sys.stdout
        if self.print_stdout:
            sys.stdout.write(message)

        # write in the log file
        self.log_file_id.write(message)
        self.log_file_id.flush()
        os.fsync(self.log_file_id.fileno())

    #---------------

    def get_log_file(self):
        '''
        Get the current log file name
        '''

        return self.log_file

    #---------------

    def __del__(self):
        '''
        Execute actions correspending to the object removal.
        '''

        # close the local log file
        self.log_file_id.close()

    #---------------

#-------------------------------------------------------------------------------

class DevNull(object):
    '''
    This class is used when it is necessary do not write a output
    '''

    #---------------

    def write(self, *_):
        '''
        Do not write anything.
        '''

        pass

    #---------------

#-------------------------------------------------------------------------------
 
class ProgramException(Exception):
    '''
    This class controls various exceptions that can occur in the execution of the application.
    '''

   #---------------

    def __init__(self, code_exception, param1='', param2='', param3=''):
        '''
        Execute actions correspending to the creation of an instance to manage a passed exception.
        ''' 

        if code_exception == 'C001':
            print('*** ERROR {0}: The application do not work if config files are not OK.'.format(code_exception), file=sys.stderr)
            sys.exit(1)
        elif code_exception == 'C002':
            print('*** ERROR {0}: The application do not work if the environment file is not OK.'.format(code_exception), file=sys.stderr)
            sys.exit(1)
        elif code_exception == 'EXIT':
            sys.exit(0)
        elif code_exception == 'P001':
            print('*** ERROR {0}: This program has parameters with invalid values.'.format(code_exception), file=sys.stderr)
            sys.exit(1)
        elif code_exception == 'S001':
            print('*** ERROR {0}: There are libraries are not installed.'.format(code_exception), file=sys.stderr)
            sys.exit(1)
        elif code_exception == 'S002':
            print('*** ERROR {0}: There is infrastructure software not installed.'.format(code_exception), file=sys.stderr)
            sys.exit(1)
        elif code_exception == 'S003':
            print('*** ERROR {0}: Boto3 - {1}.'.format(code_exception, param1), file=sys.stderr)
            sys.exit(1)
        else:
            print('*** ERROR {0}: This exception is not managed.'.format(code_exception), file=sys.stderr)
            sys.exit(1)

   #---------------

#-------------------------------------------------------------------------------

class NestedDefaultDict(collections.defaultdict):
    '''
    This class is used to create nested dictionaries.
    '''

    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))

#-------------------------------------------------------------------------------

class BreakAllLoops(Exception):
    '''
    This class is used to break out of nested loops
    '''

    pass

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This source contains general functions and classes used in {0} software package used in both console mode and gui mode.'.format(get_project_name()))
     sys.exit(0)

#-------------------------------------------------------------------------------
