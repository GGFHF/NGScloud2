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
This file contains the functions related to forms corresponding BioInfo application
menu items in console mode.
'''

#-------------------------------------------------------------------------------

import subprocess
import sys

import cinputs
import clib
import xbioinfoapp
import xbusco
import xcdhit
import xconfiguration
import xcufflinks
import xcutadapt
import xddradseqtools
import xdetonate
import xec2
import xfastqc
import xgmap
import xhisat2
import xhtseq
import xipyrad
import xkallisto
import xlib
import xngshelper
import xquast
import xrnaquast
import xsoapdenovo2
import xsoapdenovotrans
import xstar
import xstarcode
import xtoa
import xtophat
import xtransabyss
import xtransrate
import xtrimmomatic
import xtrinity
import xssh

#-------------------------------------------------------------------------------

def form_installation_bioinfo_app(app_code):
    '''
    Install the bioinfo application software in the cluster.
    '''

    # initialize the control variable
    OK = True

    # set the bioinfo application name
    if app_code == xlib.get_bedtools_code():
        app_name = xlib.get_bedtools_name()
    elif app_code == xlib.get_blastplus_code():
        app_name = xlib.get_blastplus_name()
    elif app_code == xlib.get_bowtie2_code():
        app_name = xlib.get_bowtie2_name()
    elif app_code == xlib.get_busco_code():
        app_name = xlib.get_busco_name()
    elif app_code == xlib.get_cd_hit_code():
        app_name = xlib.get_cd_hit_name()
    elif app_code == xlib.get_cufflinks_code():
        app_name = xlib.get_cufflinks_name()
    elif app_code == xlib.get_cutadapt_code():
        app_name = xlib.get_cutadapt_name()
    elif app_code == xlib.get_ddradseqtools_code():
        app_name = xlib.get_ddradseqtools_name()
    elif app_code == xlib.get_detonate_code():
        app_name = xlib.get_detonate_name()
    elif app_code == xlib.get_emboss_code():
        app_name = xlib.get_emboss_name()
    elif app_code == xlib.get_entrez_direct_code():
        app_name = xlib.get_entrez_direct_name()
    elif app_code == xlib.get_fastqc_code():
        app_name = xlib.get_fastqc_name()
    elif app_code == xlib.get_gmap_gsnap_code():
        app_name = xlib.get_gmap_gsnap_name()
    elif app_code == xlib.get_hisat2_code():
        app_name = xlib.get_hisat2_name()
    elif app_code == xlib.get_htseq_code():
        app_name = xlib.get_htseq_name()
    elif app_code == xlib.get_ipyrad_code():
        app_name = xlib.get_ipyrad_name()
    elif app_code == xlib.get_kallisto_code():
        app_name = xlib.get_kallisto_name()
    elif app_code == xlib.get_miniconda3_code():
        app_name = xlib.get_miniconda3_name()
    elif app_code == xlib.get_ngshelper_code():
        app_name = xlib.get_ngshelper_name()
    elif app_code == xlib.get_quast_code():
        app_name = xlib.get_quast_name()
    elif app_code == xlib.get_r_code():
        app_name = xlib.get_r_name()
    elif app_code == xlib.get_rnaquast_code():
        app_name = xlib.get_rnaquast_name()
    elif app_code == xlib.get_rsem_code():
        app_name = xlib.get_rsem_name()
    elif app_code == xlib.get_samtools_code():
        app_name = xlib.get_samtools_name()
    elif app_code == xlib.get_soapdenovo2_code():
        app_name = xlib.get_soapdenovo2_name()
    elif app_code == xlib.get_soapdenovotrans_code():
        app_name = xlib.get_soapdenovotrans_name()
    elif app_code == xlib.get_star_code():
        app_name = xlib.get_star_name()
    elif app_code == xlib.get_starcode_code():
        app_name = xlib.get_starcode_name()
    elif app_code == xlib.get_toa_code():
        app_name = xlib.get_toa_name()
    elif app_code == xlib.get_tophat_code():
        app_name = xlib.get_tophat_name()
    elif app_code == xlib.get_transabyss_code():
        app_name = xlib.get_transabyss_name()
    elif app_code == xlib.get_transdecoder_code():
        app_name = xlib.get_transdecoder_name()
    elif app_code == xlib.get_transrate_code():
        app_name = xlib.get_transrate_name()
    elif app_code == xlib.get_trimmomatic_code():
        app_name = xlib.get_trimmomatic_name()
    elif app_code == xlib.get_trinity_code():
        app_name = xlib.get_trinity_name()
    elif app_code == xlib.get_vsearch_code():
        app_name = xlib.get_vsearch_name()

    # get the version and download URL of the BioInfo application
    (bioinfoapp_version, bioinfoapp__url) = xconfiguration.get_bioinfo_app_data(app_name)

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Install software'.format(app_name))

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # get the version
    if bioinfoapp_version == 'any':
        version = cinputs.input_bioinfoapp_version(app_name)
    else:
        version = bioinfoapp_version

    # confirm the software install
    if OK:
        print(xlib.get_separator())
        if app_code == xlib.get_miniconda3_code():
            OK = clib.confirm_action('{0} (Bioconda infrastructure) is going to be installed in the cluster {1}. All Bioconda packages previously installed will be lost and they have to be reinstalled.'.format(app_name, cluster_name))
        elif app_code == xlib.get_r_code():
            OK = clib.confirm_action('{0} and analysis packages are going to be installed in the cluster {1}. The previous version will be lost, if it exists.'.format(app_name, cluster_name))
        elif app_code in [xlib.get_ddradseqtools_code(), xlib.get_ngshelper_code(), xlib.get_rnaquast_code(), xlib.get_transrate_code()]:
            OK = clib.confirm_action('{0} is going to be installed in the cluster {1}. The previous version will be lost, if it exists.'.format(app_name, cluster_name))
        else:
            OK = clib.confirm_action('The {0} Bioconda/Conda package is going to be installed in the cluster {1}. The previous version will be lost, if it exists.'.format(app_name, cluster_name))

    # install the software
    if OK:

        # install the BEDtools software
        if app_code == xlib.get_bedtools_code():
            package_code_list = [(xlib.get_bedtools_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the BLAST+ software
        elif app_code == xlib.get_blastplus_code():
            package_code_list = [(xlib.get_blastplus_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Bowtie2 software
        elif app_code == xlib.get_bowtie2_code():
            package_code_list = [(xlib.get_bowtie2_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the BUSCO software
        elif app_code == xlib.get_busco_code():
            package_code_list = [(xlib.get_busco_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the CD-HIT software
        elif app_code == xlib.get_cd_hit_code():
            package_code_list = [(xlib.get_cd_hit_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Cufflinks software
        elif app_code == xlib.get_cufflinks_code():
            package_code_list = [(xlib.get_cufflinks_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the cutadapt software
        elif app_code == xlib.get_cutadapt_code():
            package_code_list = [(xlib.get_cutadapt_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the ddRADseqTools software
        elif app_code == xlib.get_ddradseqtools_code():
            devstdout = xlib.DevStdOut(xddradseqtools.install_ddradseqtools.__name__)
            OK = xddradseqtools.install_ddradseqtools(cluster_name, devstdout, function=None)

        # install the DETONATE software
        elif app_code == xlib.get_detonate_code():
            package_code_list = [(xlib.get_detonate_bioconda_code(), version), (xlib.get_bowtie2_bioconda_code(), 'last')]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the EMBOSS software
        elif app_code == xlib.get_emboss_code():
            package_code_list = [(xlib.get_emboss_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Entrez Direct software
        elif app_code == xlib.get_entrez_direct_code():
            package_code_list = [(xlib.get_entrez_direct_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the FastQC software
        elif app_code == xlib.get_fastqc_code():
            package_code_list = [(xlib.get_fastqc_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the GMAP-GSNAP software
        elif app_code == xlib.get_gmap_gsnap_code():
            package_code_list = [(xlib.get_gmap_gsnap_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the HISAT2 software
        elif app_code == xlib.get_hisat2_code():
            package_code_list = [(xlib.get_hisat2_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the HTSeq software
        elif app_code == xlib.get_htseq_code():
            package_code_list = [(xlib.get_htseq_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the ipyrad software
        elif app_code == xlib.get_ipyrad_code():
            package_code_list = [xlib.get_ipyrad_conda_code()]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_conda_package_list.__name__)
            OK = xbioinfoapp.install_conda_package_list(app_code, app_name, 2, xlib.get_ipyrad_conda_code(), package_code_list, cluster_name, devstdout, function=None)

        # install the kallisto software
        elif app_code == xlib.get_kallisto_code():
            package_code_list = [(xlib.get_kallisto_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Miniconda3 software
        elif app_code == xlib.get_miniconda3_code():
            devstdout = xlib.DevStdOut(xbioinfoapp.install_miniconda3.__name__)
            OK = xbioinfoapp.install_miniconda3(cluster_name, devstdout, function=None)

        # install the NGShelper software
        elif app_code == xlib.get_ngshelper_code():
            devstdout = xlib.DevStdOut(xngshelper.install_ngshelper.__name__)
            OK = xngshelper.install_ngshelper(cluster_name, devstdout, function=None)

        # install the QUAST software
        elif app_code == xlib.get_quast_code():
            package_code_list = [(xlib.get_quast_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install R and analysis packages
        elif app_code == xlib.get_r_code():
            devstdout = xlib.DevStdOut(xbioinfoapp.install_r.__name__)
            OK = xbioinfoapp.install_r(cluster_name, devstdout, function=None)

        # install the rnaQUAST software
        elif app_code == xlib.get_rnaquast_code():
            package_code_list = [(xlib.get_rnaquast_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the RSEM software
        elif app_code == xlib.get_rsem_code():
            package_code_list = [(xlib.get_rsem_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the SAMtools software
        if app_code == xlib.get_samtools_code():
            package_code_list = [(xlib.get_samtools_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the SOAPdenovo2 software
        elif app_code == xlib.get_soapdenovo2_code():
            package_code_list = [(xlib.get_soapdenovo2_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the SOAPdenovo-Trans software
        elif app_code == xlib.get_soapdenovotrans_code():
            package_code_list = [(xlib.get_soapdenovotrans_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the STAR software
        elif app_code == xlib.get_star_code():
            package_code_list = [(xlib.get_star_bioconda_code(), version), (xlib.get_bowtie2_bioconda_code(), 'last')]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the starcode software
        elif app_code == xlib.get_starcode_code():
            package_code_list = [(xlib.get_starcode_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the TOA software
        elif app_code == xlib.get_toa_code():
            devstdout = xlib.DevStdOut(xtoa.install_toa.__name__)
            OK = xtoa.install_toa(cluster_name, devstdout, function=None)

        # install the TopHat software
        elif app_code == xlib.get_tophat_code():
            package_code_list = [(xlib.get_tophat_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Trans-ABySS software
        elif app_code == xlib.get_transabyss_code():
            package_code_list = [(xlib.get_transabyss_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the TransDecoder software
        elif app_code == xlib.get_transdecoder_code():
            package_code_list = [(xlib.get_transdecoder_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Transrate software
        elif app_code == xlib.get_transrate_code():
            devstdout = xlib.DevStdOut(xtransrate.install_transrate.__name__)
            OK = xtransrate.install_transrate(cluster_name, devstdout, function=None)
            # -- package_code_list = [(xlib.get_transrate_bioconda_code(), version)]
            # -- devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            # -- OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Trimmomatic software
        elif app_code == xlib.get_trimmomatic_code():
            package_code_list = [(xlib.get_trimmomatic_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the Trinity software
        elif app_code == xlib.get_trinity_code():
            package_code_list = [(xlib.get_trinity_bioconda_code(), version), (xlib.get_bowtie2_bioconda_code(), 'last')]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

        # install the vsearch software
        elif app_code == xlib.get_vsearch_code():
            package_code_list = [(xlib.get_vsearch_bioconda_code(), version)]
            devstdout = xlib.DevStdOut(xbioinfoapp.install_bioconda_package_list.__name__)
            OK = xbioinfoapp.install_bioconda_package_list(app_code, app_name, package_code_list, cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_busco_config_file():
    '''
    Recreate the BUSCO config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_busco_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the BUSCO config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xbusco.get_busco_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xbusco.create_busco_config_file(experiment_id, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_cd_hit_est_config_file():
    '''
    Recreate the CD-HIT-EST config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_cd_hit_est_name()))

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the CD-HIT-EST config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xcdhit.get_cd_hit_est_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xcdhit.create_cd_hit_est_config_file(experiment_id, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_cuffdiff_config_file():
    '''
    Recreate the Cuffdiff config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_cuffdiff_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_cufflinks_cuffmerge_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list,'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the quantitation dataset identification
    if OK:
        app_list = [xlib.get_cuffquant_code()]
        quantitation_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'quantitation', app_list,'uncompressed', help=True)
        if quantitation_dataset_id == '':
            print('WARNING: The cluster {0} does not have quantitation datasets.'.format(cluster_name))
            OK = False

    # recreate the Cuffdiff config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xcufflinks.get_cuffdiff_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xcufflinks.create_cuffdiff_config_file(experiment_id, assembly_dataset_id, quantitation_dataset_id)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_cufflinks_cuffmerge_config_file():
    '''
    Recreate the Cufflinks-Cuffmerge config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_cufflinks_cuffmerge_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the annotation file
    if OK:
        annotation_file = cinputs.input_gtf_file(ssh_client, reference_dataset_id, help=True)
        if annotation_file == '':
            print('ERROR: The reference dataset {0} does not have annotation files.'.format(reference_dataset_id))
            OK = False

    # get the mask file
    if OK:
        mask_file = cinputs.input_reference_file2(ssh_client, reference_dataset_id, type='mask', allowed_none=True, help=True)
        if mask_file == '':
            print('ERROR: The reference dataset {0} does not have files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the alignment dataset identification list
    if OK:
        app_list = [xlib.get_star_code(), xlib.get_tophat_code()]
        alignment_dataset_id_list = cinputs.input_result_dataset_id_list(ssh_client, experiment_id, 'alignment', app_list, 'uncompressed', help=True)
        if alignment_dataset_id_list == []:
            print('ERROR: The cluster {0} does not have alignment datasets or you did not select them.'.format(cluster_name))
            OK = False

    # recreate the Cufflinks-Cuffmerge config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xcufflinks.get_cufflinks_cuffmerge_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xcufflinks.create_cufflinks_cuffmerge_config_file(experiment_id, reference_dataset_id, reference_file, annotation_file, mask_file, alignment_dataset_id_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_cuffquant_config_file():
    '''
    Recreate the Cuffquant config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_cuffquant_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the mask file
    if OK:
        mask_file = cinputs.input_reference_file2(ssh_client, reference_dataset_id, type='mask', allowed_none=True, help=True)
        if mask_file == '':
            print('ERROR: The reference dataset {0} does not have files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the alignment dataset identification list
    if OK:
        app_list = [xlib.get_star_code(), xlib.get_tophat_code()]
        alignment_dataset_id_list = cinputs.input_result_dataset_id_list(ssh_client, experiment_id, 'alignment', app_list, 'uncompressed', help=True)
        if alignment_dataset_id_list == []:
            print('ERROR: The cluster {0} does not have alignment datasets or you did not select them.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_cufflinks_cuffmerge_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list,'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # recreate the Cuffquant config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xcufflinks.get_cuffquant_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xcufflinks.create_cuffquant_config_file(experiment_id, reference_dataset_id, mask_file, alignment_dataset_id_list, assembly_dataset_id)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_cutadapt_config_file():
    '''
    Recreate the cutadapt config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_cutadapt_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the cutadapt config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xcutadapt.get_cutadapt_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xcutadapt.create_cutadapt_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xcutadapt.create_cutadapt_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_ddradseq_simulation_config_file():
    '''
    Recreate the ddRADseq simulation config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_ddradseq_simulation_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the dictionary of restriction enzymes
    if OK:
        (OK, error_list, restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        for error in error_list:
            print(error)

    # get the enzime identification list
    if OK:
        enzyme_id_list = list(restriction_enzyme_dict.keys())
        enzyme_id_list.sort()

    # get the enzyme 1
    if OK:
        enzyme1 = cinputs.input_enzyme('1', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)

    # get the enzyme 2
    if OK:
        enzyme2 = cinputs.input_enzyme('2', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)

    # check that enzyme 1 has to be different from enzyme 2
    if OK:
        if enzyme1 in enzyme_id_list:
            enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
        else:
            enzyme1_seq = enzyme1
        if enzyme2 in enzyme_id_list:
            enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
        else:
            enzyme2_seq = enzyme2
        if enzyme1_seq.upper() == enzyme2_seq.upper():
            print('ERROR: Both enzymes have the same sequence. A ddRADseq experiment has to be performed with two different enzymes.')
            OK = False

    # recreate the ddRADseq simulation config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xddradseqtools.get_ddradseq_simulation_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xddradseqtools.create_ddradseq_simulation_config_file(reference_dataset_id, reference_file, enzyme1, enzyme2)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_fastqc_config_file():
    '''
    Recreate the FastQC config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_fastqc_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # create the FastQC config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xfastqc.get_fastqc_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xfastqc.create_fastqc_config_file(experiment_id, read_dataset_id, selected_file_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_ggtrinity_config_file():
    '''
    Recreate the Genome-guided Trinity config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_ggtrinity_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the alignment dataset identification
    if OK:
        app_list = [xlib.get_star_code(), xlib.get_tophat_code()]
        alignment_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'alignment', app_list, 'uncompressed', help=True)
        if alignment_dataset_id == '':
            print('ERROR: The cluster {0} does not have alignment datasets.'.format(cluster_name))
            OK = False

    # recreate the Genome-guided Trinity config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtrinity.get_ggtrinity_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xtrinity.create_ggtrinity_config_file(experiment_id, alignment_dataset_id)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_gmap_config_file():
    '''
    Recreate the GMAP config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_gmap_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('ERROR: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the GAMP config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xgmap.get_gmap_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xgmap.create_gmap_config_file(experiment_id, reference_dataset_id, reference_file, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_gsnap_config_file():
    '''
    Recreate the GSNAP config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_gsnap_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=True, help=True)
        if reference_dataset_id == '':
            print('WARNING: The cluster {0} does not have reference datasets. NONE is assumed as value You have to select an assembly dataset.'.format(cluster_name))
            reference_dataset_id = 'NONE'

    # get the reference file
    if OK:
        if reference_dataset_id != 'NONE':
            reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
            if reference_file == '':
                print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
                OK = False

        else:
            reference_file = 'NONE'

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        if reference_dataset_id == 'NONE':
            app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
            assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
            if assembly_dataset_id == '':
                print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
                OK = False
        else:
            assembly_dataset_id = 'NONE'

    # get the assembly type
    if OK:
        if reference_dataset_id == 'NONE':
            if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
                assembly_type = cinputs.input_assembly_type(help=True)
            elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
                assembly_type = 'NONE'
        else:
            assembly_type = 'NONE'

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the GSNAP config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xgmap.get_gsnap_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xgmap.create_gsnap_config_file(experiment_id, reference_dataset_id, reference_file, assembly_dataset_id, assembly_type, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xgmap.create_gsnap_config_file(experiment_id, reference_dataset_id, reference_file, assembly_dataset_id, assembly_type, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_hisat2_config_file():
    '''
    Recreate the HISAT2 config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_hisat2_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the splice site file
    if OK:
        splice_site_file = cinputs.input_reference_file2(ssh_client, reference_dataset_id, type='splice site', allowed_none=True, help=True)
        if splice_site_file == '':
            print('ERROR: The reference dataset {0} does not have files.'.format(reference_dataset_id))
            OK = False

    # get the exon file
    if OK:
        exon_file = cinputs.input_reference_file2(ssh_client, reference_dataset_id, type='exon', allowed_none=True, help=True)
        if exon_file == '':
            print('ERROR: The reference dataset {0} does not have files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the HISAT2 config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xhisat2.get_hisat2_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xhisat2.create_hisat2_config_file(experiment_id, reference_dataset_id, reference_file, splice_site_file, exon_file, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xhisat2.create_hisat2_config_file(experiment_id, reference_dataset_id, reference_file, splice_site_file, exon_file, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_htseq_count_config_file():
    '''
    Recreate the htseq-count config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_htseq_count_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the annotation file
    if OK:
        annotation_file = cinputs.input_gtf_file(ssh_client, reference_dataset_id, help=True)
        if annotation_file == '':
            print('ERROR: The reference dataset {0} does not have annotation files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the alignment dataset identification list
    if OK:
        app_list = [xlib.get_star_code(), xlib.get_tophat_code()]
        alignment_dataset_id_list = cinputs.input_result_dataset_id_list(ssh_client, experiment_id, 'alignment', app_list, 'uncompressed', help=True)
        if alignment_dataset_id_list == []:
            print('ERROR: The cluster {0} does not have alignment datasets or you did not select them.'.format(cluster_name))
            OK = False

    # recreate the htseq-count config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xhtseq.get_htseq_count_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xhtseq.create_htseq_count_config_file(experiment_id, reference_dataset_id, annotation_file, alignment_dataset_id_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_insilico_read_normalization_config_file():
    '''
    Recreate the insilico_read_normalization config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_insilico_read_normalization_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the insilico_read_normalization config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtrinity.get_insilico_read_normalization_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtrinity.create_insilico_read_normalization_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xtrinity.create_insilico_read_normalization_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_ipyrad_config_file():
    '''
    Recreate the ipyrad config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_ipyrad_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the assembly method
    if OK:
        assembly_method = cinputs.input_code(text='Assembly method', code_list=['DENOVO', 'REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE'], default_code=None).upper()

    # get the reference dataset identification
    if OK:
        if assembly_method in ['REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE']:
            reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
            if reference_dataset_id == '':
                print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
                OK = False
        elif ['DENOVO']:
            reference_dataset_id = 'NONE'

    # get the reference file
    if OK:
        if assembly_method in ['REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE']:
            reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
            if reference_file == '':
                print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
                OK = False
        else:
            reference_file = 'NONE'

    # get the datatype
    if OK:
        datatype = cinputs.input_code(text='RAD datatype', code_list=['RAD', 'DDRAD', 'PAIRDDRAD', 'GBS', 'PAIRGBS'], default_code=None).upper()

    # get the dictionary of restriction enzymes
    if OK:
        (OK, error_list, restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        for error in error_list:
            print(error)

    # get the enzime identification list
    if OK:
        enzyme_id_list = list(restriction_enzyme_dict.keys())
        enzyme_id_list.sort()

    # get the enzyme 1
    if OK:
        enzyme1 = cinputs.input_enzyme('1', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)

    # get the enzyme 2
    if OK:
        if datatype in ['DDRAD', 'PAIRDDRAD']:
            enzyme2 = cinputs.input_enzyme('2', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)
        else:
            enzyme2 = 'NONE'

    # check that enzyme 1 has to be different from enzyme 2
    if OK:
        if datatype in ['DDRAD', 'PAIRDDRAD']:
            if enzyme1 in enzyme_id_list:
                enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
            else:
                enzyme1_seq = enzyme1
            if enzyme2 in enzyme_id_list:
                enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
            else:
                enzyme2_seq = enzyme2
            if enzyme1_seq.upper() == enzyme2_seq.upper():
                print('ERROR: Both enzymes have the same sequence. A ddRADseq experiment has to be performed with two different enzymes.')
                OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('ERROR: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('ERROR: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get if data are demultiplexed
    if OK:
        are_data_demultiplexed = cinputs.input_code(text='Are data demultiplexed?', code_list=['y', 'n'], default_code=None).lower()
        are_data_demultiplexed = 'YES' if are_data_demultiplexed == 'y' else 'NO'

    # recreate the ipyrad config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xipyrad.get_ipyrad_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xipyrad.create_ipyrad_config_file(experiment_id, assembly_method, reference_dataset_id, reference_file, datatype, enzyme1, enzyme2,  read_dataset_id, file_pattern, are_data_demultiplexed)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_kallisto_config_file():
    '''
    Recreate the kallisto config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_kallisto_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=True, help=True)
        if reference_dataset_id == '':
            reference_dataset_id = 'NONE'
            print('WARNING: The cluster {0} does not have reference datasets. NONE is assumed as value.'.format(cluster_name))

    # get the annotation file
    if OK:
        if reference_dataset_id.upper() != 'NONE':
            annotation_file = cinputs.input_gtf_file(ssh_client, reference_dataset_id, help=True)
            if annotation_file == '':
                print('ERROR: The reference dataset {0} does not have annotation files.'.format(reference_dataset_id))
                OK = False
        else:
            annotation_file = 'NONE'

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the kallisto config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xkallisto.get_kallisto_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xkallisto.create_kallisto_config_file(experiment_id, reference_dataset_id, annotation_file, cluster_read_dir, read_type, selected_file_list, None, assembly_dataset_id, assembly_type)
            elif read_type == 'PE':
                (OK, error_list) = xkallisto.create_kallisto_config_file(experiment_id, reference_dataset_id, annotation_file, cluster_read_dir, read_type, file_1_list, file_2_list, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_quast_config_file():
    '''
    Recreate the QUAST config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_quast_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=True, help=True)
        if reference_dataset_id == '':
            reference_dataset_id = 'NONE'
            print('WARNING: The cluster {0} does not have reference datasets. NONE is assumed as value.'.format(cluster_name))

    # get the reference file
    if OK:
        if reference_dataset_id.upper() != 'NONE':
            reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
            if reference_file == '':
                print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
                OK = False
        else:
            reference_file = 'NONE'

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the QUAST config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xquast.get_quast_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xquast.create_quast_config_file(experiment_id, reference_dataset_id, reference_file, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_ref_eval_config_file():
    '''
    Recreate the REF-EVAL config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_ref_eval_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('WARNING: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the REF-EVAL config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xdetonate.get_ref_eval_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xdetonate.create_ref_eval_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None, assembly_dataset_id, assembly_type)
            elif read_type == 'PE':
                (OK, error_list) = xdetonate.create_ref_eval_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_rnaquast_config_file():
    '''
    Recreate the rnaQUAST config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_rnaquast_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=True, help=True)
        if reference_dataset_id == '':
            reference_dataset_id = 'NONE'
            print('WARNING: The cluster {0} does not have reference datasets. NONE is assumed as value.'.format(cluster_name))

    # get the reference file
    if OK:
        if reference_dataset_id.upper() != 'NONE':
            reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
            if reference_file == '':
                print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
                OK = False
        else:
            reference_file = 'NONE'

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the rnaQUAST config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xrnaquast.get_rnaquast_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xrnaquast.create_rnaquast_config_file(experiment_id, reference_dataset_id, reference_file, cluster_read_dir, read_type, selected_file_list, None, assembly_dataset_id, assembly_type)
            elif read_type == 'PE':
                (OK, error_list) = xrnaquast.create_rnaquast_config_file(experiment_id, reference_dataset_id, reference_file, cluster_read_dir, read_type, file_1_list, file_2_list, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_rsem_eval_config_file():
    '''
    Recreate the RSEM-EVAL config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_ref_eval_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the RSEM-EVAL config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xdetonate.get_rsem_eval_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xdetonate.create_rsem_eval_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None, assembly_dataset_id, assembly_type)
            elif read_type == 'PE':
                (OK, error_list) = xdetonate.create_rsem_eval_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_rsitesearch_config_file():
    '''
    Recreate the rsitesearch config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_rsitesearch_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the dictionary of restriction enzymes
    if OK:
        (OK, error_list, restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        for error in error_list:
            print(error)

    # get the enzime identification list
    if OK:
        enzyme_id_list = list(restriction_enzyme_dict.keys())
        enzyme_id_list.sort()

    # get the enzyme 1
    if OK:
        enzyme1 = cinputs.input_enzyme('1', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)

    # get the enzyme 2
    if OK:
        enzyme2 = cinputs.input_enzyme('2', restriction_enzyme_dict, allowed_ambiguity_codes=True, help=True)

    # check if enzyme 1 is different or not to enzyme 2
    if OK:
        if enzyme1 in enzyme_id_list:
            enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
        else:
            enzyme1_seq = enzyme1
        if enzyme2 in enzyme_id_list:
            enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
        else:
            enzyme2_seq = enzyme2
        if enzyme1_seq.upper() == enzyme2_seq.upper():
            print('Both enzymes have the same sequence. An enzyme analysis of a RAD-seq experiment will be performed.')
        else:
            print('Both enzymes have different sequences. An enzyme analysis of a ddRADseq experiment will be performed.')

    # recreate the rsitesearch config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xddradseqtools.get_rsitesearch_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xddradseqtools.create_rsitesearch_config_file(reference_dataset_id, reference_file, enzyme1, enzyme2)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_soapdenovo2_config_file():
    '''
    Recreate the SOAPdenovo2 application config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_soapdenovo2_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the SOAPdenovo2 config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xsoapdenovo2.get_soapdenovo2_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xsoapdenovo2.create_soapdenovo2_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xsoapdenovo2.create_soapdenovo2_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_soapdenovotrans_config_file():
    '''
    Recreate the SOAPdenovo-Trans application config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_soapdenovotrans_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the SOAPdenovo-Trans config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xsoapdenovotrans.get_soapdenovotrans_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xsoapdenovotrans.create_soapdenovotrans_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xsoapdenovotrans.create_soapdenovotrans_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_star_config_file():
    '''
    Recreate the STAR config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_star_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the annotation file
    if OK:
        annotation_file = cinputs.input_gtf_file(ssh_client, reference_dataset_id, help=True)
        if annotation_file == '':
            print('ERROR: The reference dataset {0} does not have annotation files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('ERROR: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('ERROR: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the STAR config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xstar.get_star_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xstar.create_star_config_file(experiment_id, reference_dataset_id, reference_file, annotation_file, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xstar.create_star_config_file(experiment_id, reference_dataset_id, reference_file, annotation_file, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_starcode_config_file():
    '''
    Recreate the starcode config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_starcode_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the starcode config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xstarcode.get_starcode_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xstarcode.create_starcode_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xstarcode.create_starcode_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_tophat_config_file():
    '''
    Recreate the TopHat config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_tophat_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('ERROR: The cluster {0} does not have reference datasets.'.format(cluster_name))
            OK = False

    # get the reference file
    if OK:
        reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
        if reference_file == '':
            print('ERROR: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
            OK = False

    # get the annotation file
    if OK:
        annotation_file = cinputs.input_gtf_file(ssh_client, reference_dataset_id, help=True)
        if annotation_file == '':
            print('ERROR: The reference dataset {0} does not have annotation files.'.format(reference_dataset_id))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('ERROR: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('ERROR: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('ERROR: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the TopHat config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtophat.get_tophat_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtophat.create_tophat_config_file(experiment_id, reference_dataset_id, reference_file, annotation_file, read_dataset_id, read_type, selected_file_list[0], None)
            elif read_type == 'PE':
                (OK, error_list) = xtophat.create_tophat_config_file(experiment_id, reference_dataset_id, reference_file, annotation_file, read_dataset_id, read_type, file_1_list[0], file_2_list[0])
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_transabyss_config_file():
    '''
    Recreate the Trans-ABySS config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_transabyss_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the Trans-ABySS config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtransabyss.get_transabyss_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtransabyss.create_transabyss_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xtransabyss.create_transabyss_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_transcript_filter_config_file():
    '''
    Recreate the transcript-filter config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_transcript_filter_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the RSEM-EVAL dataset identification
    if OK:
        app_list = [xlib.get_rsem_eval_code()]
        rsem_eval_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'RSEM-EVAL', app_list, 'uncompressed', help=True)
        if rsem_eval_dataset_id == '':
            print('WARNING: The cluster {0} does not have RSEM-EVAL datasets.'.format(cluster_name))
            OK = False

    # recreate the transcripts-filter config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xngshelper.get_transcript_filter_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xngshelper.create_transcript_filter_config_file(experiment_id, rsem_eval_dataset_id)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_transcriptome_blastx_config_file():
    '''
    Recreate the transcriptome-blastx config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_transcriptome_blastx_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the database dataset identification
    if OK:
        database_dataset_id = cinputs.input_database_dataset_id(ssh_client, help=True)
        if database_dataset_id == '':
            print('WARNING: The cluster {0} does not have any database dataset.'.format(cluster_name))
            OK = False

    # get the protein database name
    if OK:
        protein_database_name = cinputs.input_protein_database_name(ssh_client, database_dataset_id, help=True)
        if protein_database_name == '':
            print('WARNING: The dataset {0} in the cluster {1} does not have any protein database.'.format(database_dataset_id, cluster_name))
            OK = False

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the transcriptome-blastx config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xngshelper.get_transcriptome_blastx_config_file()))

        # recreate the config file
        if OK:
            (OK, error_list) = xngshelper.create_transcriptome_blastx_config_file(database_dataset_id, protein_database_name, experiment_id, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_transrate_config_file():
    '''
    Recreate the Transrate config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_transrate_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=True, help=True)
        if reference_dataset_id == '':
            reference_dataset_id = 'NONE'
            print('WARNING: The cluster {0} does not have reference datasets. NONE is assumed as value.'.format(cluster_name))

    # get the reference file
    if OK:
        if reference_dataset_id.upper() != 'NONE':
            reference_file = cinputs.input_reference_file(ssh_client, reference_dataset_id, help=True)
            if reference_file == '':
                print('WARNING: The reference dataset {0} does not have reference files.'.format(reference_dataset_id))
                OK = False
        else:
            reference_file = 'NONE'

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # get the assembly dataset identification
    if OK:
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
        if assembly_dataset_id == '':
            print('WARNING: The cluster {0} does not have assembly datasets.'.format(cluster_name))
            OK = False

    # get the assembly type
    if OK:
        if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            assembly_type = cinputs.input_assembly_type(help=True)
        elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            assembly_type = 'NONE'

    # recreate the Transrate config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtransrate.get_transrate_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtransrate.create_transrate_config_file(experiment_id, reference_dataset_id, reference_file, cluster_read_dir, read_type, selected_file_list, None, assembly_dataset_id, assembly_type)
            elif read_type == 'PE':
                (OK, error_list) = xtransrate.create_transrate_config_file(experiment_id, reference_dataset_id, reference_file, cluster_read_dir, read_type, file_1_list, file_2_list, assembly_dataset_id, assembly_type)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_trimmomatic_config_file():
    '''
    Recreate the Trimmomatic config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_trimmomatic_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the Trimmomatic config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtrimmomatic.get_trimmomatic_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtrimmomatic.create_trimmomatic_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xtrimmomatic.create_trimmomatic_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_trinity_config_file():
    '''
    Recreate the Trinity config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Recreate config file'.format(xlib.get_trinity_name()))

    # get the cluster name, experiment identification, read dataset identification and the file pattern
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print('WARNING: The cluster {0} does not have read datasets.'.format(cluster_name))
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*fastq')

    # build the cluster read directory path
    if OK:
        cluster_read_dir = '{0}/{1}/{2}'.format(xlib.get_cluster_read_dir(), experiment_id, read_dataset_id)

    # get the selected file list
    if OK:
        selected_file_list = []
        command = 'cd {0}; find . -type f -regex "./{1}"'.format(cluster_read_dir, file_pattern)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
        else:
            print('*** ERROR: Wrong command ---> {0}'.format(command))
        if selected_file_list == []:
            print('WARNING: There are not files in the cluster directory {0} with the pattern {1}'.format(cluster_read_dir, file_pattern))
            OK = False

    # get the read type
    if OK:
        read_type = cinputs.input_read_type()

    # get the specific_chars to identify files when the read type is paired 
    if OK:
        if read_type == 'SE':
            specific_chars_1 = None
            specific_chars_2 = None
        elif read_type == 'PE':
            specific_chars_1 = cinputs.input_file_pairing_specific_chars(1, '1.fastq')
            specific_chars_2 = cinputs.input_file_pairing_specific_chars(2, '2.fastq')

    # get the paired file list when the read type is paired
    if OK:
        if read_type == 'PE':
            (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, specific_chars_1, specific_chars_2)
            if unpaired_file_list != []:
                print('ERROR: There are unpaired files: {0}'.format(unpaired_file_list))
                OK = False

    # recreate the Trinity config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(xtrinity.get_trinity_config_file()))

        # recreate the config file
        if OK:
            if read_type == 'SE':
                (OK, error_list) = xtrinity.create_trinity_config_file(experiment_id, read_dataset_id, read_type, selected_file_list, None)
            elif read_type == 'PE':
                (OK, error_list) = xtrinity.create_trinity_config_file(experiment_id, read_dataset_id, read_type, file_1_list, file_2_list)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_bioinfo_config_file(app):
    '''
    Edit a bioinfo appliation config file to change the parameters of each process.
    '''

    # initialize the control variable
    OK = True

    # set the bioinfo application name
    if app == xlib.get_busco_code():
        name = xlib.get_busco_name()
    elif app == xlib.get_cd_hit_est_code():
        name = xlib.get_cd_hit_est_name()
    elif app == xlib.get_cuffdiff_code():
        name = xlib.get_cuffdiff_name()
    elif app == xlib.get_cufflinks_cuffmerge_code():
        name = xlib.get_cufflinks_cuffmerge_name()
    elif app == xlib.get_cuffquant_code():
        name = xlib.get_cuffquant_name()
    elif app == xlib.get_cutadapt_code():
        name = xlib.get_cutadapt_name()
    elif app == xlib.get_ddradseq_simulation_code():
        name = xlib.get_ddradseq_simulation_name()
    elif app == xlib.get_fastqc_code():
        name = xlib.get_fastqc_name()
    elif app == xlib.get_ggtrinity_code():
        name = xlib.get_ggtrinity_name()
    elif app == xlib.get_gmap_code():
        name = xlib.get_gmap_name()
    elif app == xlib.get_gsnap_code():
        name = xlib.get_gsnap_name()
    elif app == xlib.get_hisat2_code():
        name = xlib.get_hisat2_name()
    elif app == xlib.get_htseq_count_code():
        name = xlib.get_htseq_count_name()
    elif app == xlib.get_insilico_read_normalization_code():
        name = xlib.get_insilico_read_normalization_name()
    elif app == xlib.get_ipyrad_code():
        name = xlib.get_ipyrad_name()
    elif app == xlib.get_kallisto_code():
        name = xlib.get_kallisto_name()
    elif app == xlib.get_quast_code():
        name = xlib.get_quast_name()
    elif app == xlib.get_ref_eval_code():
        name = xlib.get_ref_eval_name()
    elif app == xlib.get_rnaquast_code():
        name = xlib.get_rnaquast_name()
    elif app == xlib.get_rsem_eval_code():
        name = xlib.get_rsem_eval_name()
    elif app == xlib.get_rsitesearch_code():
        name = xlib.get_rsitesearch_name()
    elif app == xlib.get_soapdenovo2_code():
        name = xlib.get_soapdenovo2_name()
    elif app == xlib.get_soapdenovotrans_code():
        name = xlib.get_soapdenovotrans_name()
    elif app == xlib.get_star_code():
        name = xlib.get_star_name()
    elif app == xlib.get_starcode_code():
        name = xlib.get_starcode_name()
    elif app == xlib.get_tophat_code():
        name = xlib.get_tophat_name()
    elif app == xlib.get_transabyss_code():
        name = xlib.get_transabyss_name()
    elif app == xlib.get_transcript_filter_code():
        name = xlib.get_transcript_filter_name()
    elif app == xlib.get_transcriptome_blastx_code():
        name = xlib.get_transcriptome_blastx_name()
    elif app == xlib.get_transrate_code():
        name = xlib.get_transrate_name()
    elif app == xlib.get_trimmomatic_code():
        name = xlib.get_trimmomatic_name()
    elif app == xlib.get_trinity_code():
        name = xlib.get_trinity_name()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Edit config file'.format(name))

    # get the config file
    if app == xlib.get_busco_code():
        config_file = xbusco.get_busco_config_file()
    elif app == xlib.get_cd_hit_est_code():
        config_file = xcdhit.get_cd_hit_est_config_file()
    elif app == xlib.get_cuffdiff_code():
        config_file = xcufflinks.get_cuffdiff_config_file()
    elif app == xlib.get_cufflinks_cuffmerge_code():
        config_file = xcufflinks.get_cufflinks_cuffmerge_config_file()
    elif app == xlib.get_cuffquant_code():
        config_file = xcufflinks.get_cuffquant_config_file()
    elif app == xlib.get_cutadapt_code():
        config_file = xcutadapt.get_cutadapt_config_file()
    elif app == xlib.get_ddradseq_simulation_code():
        config_file = xddradseqtools.get_ddradseq_simulation_config_file()
    elif app == xlib.get_fastqc_code():
        config_file = xfastqc.get_fastqc_config_file()
    elif app == xlib.get_ggtrinity_code():
        config_file = xtrinity.get_ggtrinity_config_file()
    elif app == xlib.get_gmap_code():
        config_file = xgmap.get_gmap_config_file()
    elif app == xlib.get_gsnap_code():
        config_file = xgmap.get_gsnap_config_file()
    elif app == xlib.get_hisat2_code():
        config_file = xhisat2.get_hisat2_config_file()
    elif app == xlib.get_htseq_count_code():
        config_file = xhtseq.get_htseq_count_config_file()
    elif app == xlib.get_insilico_read_normalization_code():
        config_file = xtrinity.get_insilico_read_normalization_config_file()
    elif app == xlib.get_ipyrad_code():
        config_file = xipyrad.get_ipyrad_config_file()
    elif app == xlib.get_kallisto_code():
        config_file = xkallisto.get_kallisto_config_file()
    elif app == xlib.get_quast_code():
        config_file = xquast.get_quast_config_file()
    elif app == xlib.get_ref_eval_code():
        config_file = xdetonate.get_ref_eval_config_file()
    elif app == xlib.get_rnaquast_code():
        config_file = xrnaquast.get_rnaquast_config_file()
    elif app == xlib.get_rsem_eval_code():
        config_file = xdetonate.get_rsem_eval_config_file()
    elif app == xlib.get_rsitesearch_code():
        config_file = xddradseqtools.get_rsitesearch_config_file()
    elif app == xlib.get_soapdenovo2_code():
        config_file = xsoapdenovo2.get_soapdenovo2_config_file()
    elif app == xlib.get_soapdenovotrans_code():
        config_file = xsoapdenovotrans.get_soapdenovotrans_config_file()
    elif app == xlib.get_star_code():
        config_file = xstar.get_star_config_file()
    elif app == xlib.get_starcode_code():
        config_file = xstarcode.get_starcode_config_file()
    elif app == xlib.get_tophat_code():
        config_file = xtophat.get_tophat_config_file()
    elif app == xlib.get_transabyss_code():
        config_file = xtransabyss.get_transabyss_config_file()
    elif app == xlib.get_transcript_filter_code():
        config_file = xngshelper.get_transcript_filter_config_file()
    elif app == xlib.get_transcriptome_blastx_code():
        config_file = xngshelper.get_transcriptome_blastx_config_file()
    elif app == xlib.get_transrate_code():
        config_file = xtransrate.get_transrate_config_file()
    elif app == xlib.get_trimmomatic_code():
        config_file = xtrimmomatic.get_trimmomatic_config_file()
    elif app == xlib.get_trinity_code():
        config_file = xtrinity.get_trinity_config_file()

    # edit the read transfer config file
    print(xlib.get_separator())
    print('Editing the {0} config file ...'.format(name))
    command = '{0} {1}'.format(xlib.get_editor(), config_file)
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print('*** ERROR: Return code {0} in command -> {1}'.format(rc, command))
        OK = False

    # check the config file
    if OK:
        print(xlib.get_separator())
        print('Checking the {0} config file ...'.format(name))
        if app == xlib.get_busco_code():
            (OK, error_list) = xbusco.check_busco_config_file(strict=False)
        elif app == xlib.get_cd_hit_est_code():
            (OK, error_list) = xcdhit.check_cd_hit_est_config_file(strict=False)
        elif app == xlib.get_cuffdiff_code():
            (OK, error_list) = xcufflinks.check_cuffdiff_config_file(strict=False)
        elif app == xlib.get_cufflinks_cuffmerge_code():
            (OK, error_list) = xcufflinks.check_cufflinks_cuffmerge_config_file(strict=False)
        elif app == xlib.get_cuffquant_code():
            (OK, error_list) = xcufflinks.check_cuffquant_config_file(strict=False)
        elif app == xlib.get_cutadapt_code():
            (OK, error_list) = xcutadapt.check_cutadapt_config_file(strict=False)
        elif app == xlib.get_ddradseq_simulation_code():
            (OK, error_list) = xddradseqtools.check_ddradseq_simulation_config_file(strict=False)
        elif app == xlib.get_fastqc_code():
            (OK, error_list) = xfastqc.check_fastqc_config_file(strict=False)
        elif app == xlib.get_ggtrinity_code():
            (OK, error_list) = xtrinity.check_ggtrinity_config_file(strict=False)
        elif app == xlib.get_gmap_code():
            (OK, error_list) = xgmap.check_gmap_config_file(strict=False)
        elif app == xlib.get_gsnap_code():
            (OK, error_list) = xgmap.check_gsnap_config_file(strict=False)
        elif app == xlib.get_hisat2_code():
            (OK, error_list) = xhisat2.check_hisat2_config_file(strict=False)
        elif app == xlib.get_htseq_count_code():
            (OK, error_list) = xhtseq.check_htseq_count_config_file(strict=False)
        elif app == xlib.get_insilico_read_normalization_code():
            (OK, error_list) = xtrinity.check_insilico_read_normalization_config_file(strict=False)
        elif app == xlib.get_ipyrad_code():
            (OK, error_list) = xipyrad.check_ipyrad_config_file(strict=False)
        elif app == xlib.get_kallisto_code():
            (OK, error_list) = xkallisto.check_kallisto_config_file(strict=False)
        elif app == xlib.get_quast_code():
            (OK, error_list) = xquast.check_quast_config_file(strict=False)
        elif app == xlib.get_ref_eval_code():
            (OK, error_list) = xdetonate.check_ref_eval_config_file(strict=False)
        elif app == xlib.get_rnaquast_code():
            (OK, error_list) = xrnaquast.check_rnaquast_config_file(strict=False)
        elif app == xlib.get_rsem_eval_code():
            (OK, error_list) = xdetonate.check_rsem_eval_config_file(strict=False)
        elif app == xlib.get_rsitesearch_code():
            (OK, error_list) = xddradseqtools.check_rsitesearch_config_file(strict=False)
        elif app == xlib.get_soapdenovo2_code():
            (OK, error_list) = xsoapdenovo2.check_soapdenovo2_config_file(strict=False)
        elif app == xlib.get_soapdenovotrans_code():
            (OK, error_list) = xsoapdenovotrans.check_soapdenovotrans_config_file(strict=False)
        elif app == xlib.get_star_code():
            (OK, error_list) = xstar.check_star_config_file(strict=False)
        elif app == xlib.get_starcode_code():
            (OK, error_list) = xstarcode.check_starcode_config_file(strict=False)
        elif app == xlib.get_tophat_code():
            (OK, error_list) = xtophat.check_tophat_config_file(strict=False)
        elif app == xlib.get_transabyss_code():
            (OK, error_list) = xtransabyss.check_transabyss_config_file(strict=False)
        elif app == xlib.get_transcript_filter_code():
            (OK, error_list) = xngshelper.check_transcript_filter_config_file(strict=False)
        elif app == xlib.get_transcriptome_blastx_code():
            (OK, error_list) = xngshelper.check_transcriptome_blastx_config_file(strict=False)
        elif app == xlib.get_transrate_code():
            (OK, error_list) = xtransrate.check_transrate_config_file(strict=False)
        elif app == xlib.get_trimmomatic_code():
            (OK, error_list) = xtrimmomatic.check_trimmomatic_config_file(strict=False)
        elif app == xlib.get_trinity_code():
            (OK, error_list) = xtrinity.check_trinity_config_file(strict=False)
        if OK:
            print('The file is OK.')
        else:
            print()
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_run_bioinfo_process(app):
    '''
    Run a bioinfo application process with the parameters in the corresponding config file.
    '''

    # initialize the control variable
    OK = True

    # set the bioinfo application name
    if app == xlib.get_busco_code():
        name = xlib.get_busco_name()
    elif app == xlib.get_cd_hit_est_code():
        name = xlib.get_cd_hit_est_name()
    elif app == xlib.get_cuffdiff_code():
        name = xlib.get_cuffdiff_name()
    elif app == xlib.get_cufflinks_cuffmerge_code():
        name = xlib.get_cufflinks_cuffmerge_name()
    elif app == xlib.get_cuffquant_code():
        name = xlib.get_cuffquant_name()
    elif app == xlib.get_cutadapt_code():
        name = xlib.get_cutadapt_name()
    elif app == xlib.get_ddradseq_simulation_code():
        name = xlib.get_ddradseq_simulation_name()
    elif app == xlib.get_fastqc_code():
        name = xlib.get_fastqc_name()
    elif app == xlib.get_ggtrinity_code():
        name = xlib.get_ggtrinity_name()
    elif app == xlib.get_gmap_code():
        name = xlib.get_gmap_name()
    elif app == xlib.get_gsnap_code():
        name = xlib.get_gsnap_name()
    elif app == xlib.get_hisat2_code():
        name = xlib.get_hisat2_name()
    elif app == xlib.get_htseq_count_code():
        name = xlib.get_htseq_count_name()
    elif app == xlib.get_insilico_read_normalization_code():
        name = xlib.get_insilico_read_normalization_name()
    elif app == xlib.get_ipyrad_code():
        name = xlib.get_ipyrad_name()
    elif app == xlib.get_kallisto_code():
        name = xlib.get_kallisto_name()
    elif app == xlib.get_quast_code():
        name = xlib.get_quast_name()
    elif app == xlib.get_ref_eval_code():
        name = xlib.get_ref_eval_name()
    elif app == xlib.get_rnaquast_code():
        name = xlib.get_rnaquast_name()
    elif app == xlib.get_rsem_eval_code():
        name = xlib.get_rsem_eval_name()
    elif app == xlib.get_rsitesearch_code():
        name = xlib.get_rsitesearch_name()
    elif app == xlib.get_soapdenovo2_code():
        name = xlib.get_soapdenovo2_name()
    elif app == xlib.get_soapdenovotrans_code():
        name = xlib.get_soapdenovotrans_name()
    elif app == xlib.get_star_code():
        name = xlib.get_star_name()
    elif app == xlib.get_starcode_code():
        name = xlib.get_starcode_name()
    elif app == xlib.get_tophat_code():
        name = xlib.get_tophat_name()
    elif app == xlib.get_transabyss_code():
        name = xlib.get_transabyss_name()
    elif app == xlib.get_transcript_filter_code():
        name = xlib.get_transcript_filter_name()
    elif app == xlib.get_transcriptome_blastx_code():
        name = xlib.get_transcriptome_blastx_name()
    elif app == xlib.get_transrate_code():
        name = xlib.get_transrate_name()
    elif app == xlib.get_trimmomatic_code():
        name = xlib.get_trimmomatic_name()
    elif app == xlib.get_trinity_code():
        name = xlib.get_trinity_name()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('{0} - Run process'.format(name))

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # confirm the process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The {0} process is going to be run.'.format(name))

    # run the process
    if OK:

        # execute the process when it is a BUSCO process
        if app == xlib.get_busco_code():
            devstdout = xlib.DevStdOut(xbusco.run_busco_process.__name__)
            OK = xbusco.run_busco_process(cluster_name, devstdout, function=None)

        # execute the process when it is a CD-HIT-EST process
        elif app == xlib.get_cd_hit_est_code():
            devstdout = xlib.DevStdOut(xcdhit.run_cd_hit_est_process.__name__)
            OK = xcdhit.run_cd_hit_est_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Cuffdiff process
        elif app == xlib.get_cuffdiff_code():
            devstdout = xlib.DevStdOut(xcufflinks.run_cuffdiff_process.__name__)
            OK = xcufflinks.run_cuffdiff_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Cufflinks-Cuffmerge process
        elif app == xlib.get_cufflinks_cuffmerge_code():
            devstdout = xlib.DevStdOut(xcufflinks.run_cufflinks_cuffmerge_process.__name__)
            OK = xcufflinks.run_cufflinks_cuffmerge_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Cuffquant process
        elif app == xlib.get_cuffquant_code():
            devstdout = xlib.DevStdOut(xcufflinks.run_cuffquant_process.__name__)
            OK = xcufflinks.run_cuffquant_process(cluster_name, devstdout, function=None)

        # execute the process when it is a cutadapt process
        elif app == xlib.get_cutadapt_code():
            devstdout = xlib.DevStdOut(xcutadapt.run_cutadapt_process.__name__)
            OK = xcutadapt.run_cutadapt_process(cluster_name, devstdout, function=None)

        # execute the process when it is a ddRADseq simulation process
        elif app == xlib.get_ddradseq_simulation_code():
            devstdout = xlib.DevStdOut(xddradseqtools.run_ddradseq_simulation_process.__name__)
            OK = xddradseqtools.run_ddradseq_simulation_process(cluster_name, devstdout, function=None)

        # execute the process when it is a FastQC process
        elif app == xlib.get_fastqc_code():
            devstdout = xlib.DevStdOut(xfastqc.run_fastqc_process.__name__)
            OK = xfastqc.run_fastqc_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Genome-guided Trinity process
        elif app == xlib.get_ggtrinity_code():
            devstdout = xlib.DevStdOut(xtrinity.run_ggtrinity_process.__name__)
            OK = xtrinity.run_ggtrinity_process(cluster_name, devstdout, function=None)

        # execute the process when it is a GMAP process
        elif app == xlib.get_gmap_code():
            devstdout = xlib.DevStdOut(xgmap.run_gmap_process.__name__)
            OK = xgmap.run_gmap_process(cluster_name, devstdout, function=None)

        # execute the process when it is a GSNAP process
        elif app == xlib.get_gsnap_code():
            devstdout = xlib.DevStdOut(xgmap.run_gsnap_process.__name__)
            OK = xgmap.run_gsnap_process(cluster_name, devstdout, function=None)

        # execute the process when it is a HISAT2 process
        elif app == xlib.get_hisat2_code():
            devstdout = xlib.DevStdOut(xhisat2.run_hisat2_process.__name__)
            OK = xhisat2.run_hisat2_process(cluster_name, devstdout, function=None)

        # execute the process when it is a htseq-count process
        elif app == xlib.get_htseq_count_code():
            devstdout = xlib.DevStdOut(xhtseq.run_htseq_count_process.__name__)
            OK = xhtseq.run_htseq_count_process(cluster_name, devstdout, function=None)

        # execute the process when it is a insilico_read_normalization process
        elif app == xlib.get_insilico_read_normalization_code():
            devstdout = xlib.DevStdOut(xtrinity.run_insilico_read_normalization_process.__name__)
            OK = xtrinity.run_insilico_read_normalization_process(cluster_name, devstdout, function=None)

        # execute the process when it is a ipyrad process
        elif app == xlib.get_ipyrad_code():
            devstdout = xlib.DevStdOut(xipyrad.run_ipyrad_process.__name__)
            OK = xipyrad.run_ipyrad_process(cluster_name, devstdout, function=None)

        # execute the process when it is a kallisto process
        elif app == xlib.get_kallisto_code():
            devstdout = xlib.DevStdOut(xkallisto.run_kallisto_process.__name__)
            OK = xkallisto.run_kallisto_process(cluster_name, devstdout, function=None)

        # execute the process when it is a QUAST process
        elif app == xlib.get_quast_code():
            devstdout = xlib.DevStdOut(xquast.run_quast_process.__name__)
            OK = xquast.run_quast_process(cluster_name, devstdout, function=None)

        # execute the process when it is a REF-EVAL process
        elif app == xlib.get_ref_eval_code():
            devstdout = xlib.DevStdOut(xdetonate.run_ref_eval_process.__name__)
            OK = xdetonate.run_ref_eval_process(cluster_name, devstdout, function=None)

        # execute the process when it is a rnaQUAST process
        elif app == xlib.get_rnaquast_code():
            devstdout = xlib.DevStdOut(xrnaquast.run_rnaquast_process.__name__)
            OK = xrnaquast.run_rnaquast_process(cluster_name, devstdout, function=None)

        # execute the process when it is a RSEM-EVAL process
        elif app == xlib.get_rsem_eval_code():
            devstdout = xlib.DevStdOut(xdetonate.run_rsem_eval_process.__name__)
            OK = xdetonate.run_rsem_eval_process(cluster_name, devstdout, function=None)

        # execute the process when it is a rsitesearch process
        elif app == xlib.get_rsitesearch_code():
            devstdout = xlib.DevStdOut(xddradseqtools.run_rsitesearch_process.__name__)
            OK = xddradseqtools.run_rsitesearch_process(cluster_name, devstdout, function=None)

        # execute the process when it is a SOAPdenovo2 process
        elif app == xlib.get_soapdenovo2_code():
            devstdout = xlib.DevStdOut(xsoapdenovo2.run_soapdenovo2_process.__name__)
            OK = xsoapdenovo2.run_soapdenovo2_process(cluster_name, devstdout, function=None)

        # execute the process when it is a SOAPdenovo-Trans process
        elif app == xlib.get_soapdenovotrans_code():
            devstdout = xlib.DevStdOut(xsoapdenovotrans.run_soapdenovotrans_process.__name__)
            OK = xsoapdenovotrans.run_soapdenovotrans_process(cluster_name, devstdout, function=None)

        # execute the process when it is a STAR process
        elif app == xlib.get_star_code():
            devstdout = xlib.DevStdOut(xstar.run_star_process.__name__)
            OK = xstar.run_star_process(cluster_name, devstdout, function=None)

        # execute the process when it is a starcode process
        elif app == xlib.get_starcode_code():
            devstdout = xlib.DevStdOut(xstarcode.run_starcode_process.__name__)
            OK = xstarcode.run_starcode_process(cluster_name, devstdout, function=None)

        # execute the process when it is a TopHat process
        elif app == xlib.get_tophat_code():
            devstdout = xlib.DevStdOut(xtophat.run_tophat_process.__name__)
            OK = xtophat.run_tophat_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Trans-ABySS process
        elif app == xlib.get_transabyss_code():
            devstdout = xlib.DevStdOut(xtransabyss.run_transabyss_process.__name__)
            OK = xtransabyss.run_transabyss_process(cluster_name, devstdout, function=None)

        # execute the process when it is a transcripts-filter process
        elif app == xlib.get_transcript_filter_code():
            devstdout = xlib.DevStdOut(xngshelper.run_transcript_filter_process.__name__)
            OK = xngshelper.run_transcript_filter_process(cluster_name, devstdout, function=None)

        # execute the process when it is a transcriptome-blastx process
        elif app == xlib.get_transcriptome_blastx_code():
            devstdout = xlib.DevStdOut(xngshelper.run_transcriptome_blastx_process.__name__)
            OK = xngshelper.run_transcriptome_blastx_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Transrate process
        elif app == xlib.get_transrate_code():
            devstdout = xlib.DevStdOut(xtransrate.run_transrate_process.__name__)
            OK = xtransrate.run_transrate_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Trimmomatic process
        elif app == xlib.get_trimmomatic_code():
            devstdout = xlib.DevStdOut(xtrimmomatic.run_trimmomatic_process.__name__)
            OK = xtrimmomatic.run_trimmomatic_process(cluster_name, devstdout, function=None)

        # execute the process when it is a Trinity process
        elif app == xlib.get_trinity_code():
            devstdout = xlib.DevStdOut(xtrinity.run_trinity_process.__name__)
            OK = xtrinity.run_trinity_process(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_data_file(data_file):
    '''
    Recreate a data file.
    '''

    # get the head
    if data_file == xddradseqtools.get_restriction_site_file():
        head = '{0} - Recreate the file of restriction sites'.format(xlib.get_ddradseqtools_name())
    elif data_file == xddradseqtools.get_end_file():
        head = '{0} - Recreate the file of ends'.format(xlib.get_ddradseqtools_name())
    elif data_file == xddradseqtools.get_individual_file():
        head = '{0} - Recreate the file of individuals'.format(xlib.get_ddradseqtools_name())
    elif data_file == xtoa.get_dataset_file():
        head = '{0} - Recreate the file of genomic dataset'.format(xlib.get_toa_name())
    elif data_file == xtoa.get_species_file():
        head = '{0} - Recreate the file of species'.format(xlib.get_toa_name())

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(head)

    # confirm the creation of the data file
    print(xlib.get_separator())
    OK = clib.confirm_action('The file {0} is going to be recreated. The previous files will be lost.'.format(data_file))

    # recreate the config file
    if OK:
        if data_file == xddradseqtools.get_restriction_site_file():
            (OK, error_list) = xddradseqtools.create_restriction_site_file()
        elif data_file == xddradseqtools.get_end_file():
            (OK, error_list) = xddradseqtools.create_end_file()
        elif data_file == xddradseqtools.get_individual_file():
            (OK, error_list) = xddradseqtools.create_individual_file()
        elif data_file == xtoa.get_dataset_file():
            (OK, error_list) = xtoa.create_dataset_file()
        elif data_file == xtoa.get_species_file():
            (OK, error_list) = xtoa.create_species_file()
        if OK:
            print('The file is recreated.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_data_file(data_file):
    '''
    Edit a data file.
    '''

    # initialize the control variable
    OK = True

    # get the head
    if data_file == xddradseqtools.get_restriction_site_file():
        head = '{0} - Edit the file of restriction sites'.format(xlib.get_ddradseqtools_name())
    elif data_file == xddradseqtools.get_end_file():
        head = '{0} - Edit the file of ends'.format(xlib.get_ddradseqtools_name())
    elif data_file == xddradseqtools.get_individual_file():
        head = '{0} - Edit the file of individuals'.format(xlib.get_ddradseqtools_name())
    elif data_file == xtoa.get_dataset_file():
        head = '{0} - Edit the file of genomic dataset'.format(xlib.get_toa_name())
    elif data_file == xtoa.get_species_file():
        head = '{0} - Edit the file of species'.format(xlib.get_toa_name())

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(head)

    # edit the read transfer config file
    print(xlib.get_separator())
    print('Editing the file {0} ...'.format(data_file))
    command = '{0} {1}'.format(xlib.get_editor(), data_file)
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print('*** ERROR: Return code {0} in command -> {1}'.format(rc, command))
        OK = False

    # check the data file
    if OK:
        print(xlib.get_separator())
        print('Checking the file {0} ...'.format(data_file))
        if data_file == xddradseqtools.get_restriction_site_file():
            (OK, error_list) = xddradseqtools.check_restriction_site_file(strict=False)
        elif data_file == xddradseqtools.get_end_file():
            (OK, error_list) = xddradseqtools.check_end_file(strict=False)
        elif data_file == xddradseqtools.get_individual_file():
            (OK, error_list) = xddradseqtools.check_individual_file(strict=False)
        elif data_file == xtoa.get_dataset_file():
            (OK, error_list) = xtoa.check_dataset_file(strict=False)
        elif data_file == xtoa.get_species_file():
            (OK, error_list) = xtoa.check_species_file(strict=False)
        if OK:
            print('The file is OK.')
        else:
            print()
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the functions related to forms corresponding BioInfo application menu items in console mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
