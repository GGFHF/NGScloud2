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
This file contains the functions related to forms corresponding to dataset
menu items in mode console.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys

import cinputs
import clib
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def form_list_submission_logs():
    '''
    List the submission logs.
    '''

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Logs - List submission logs')

    # get the submission process dictionary
    submission_process_dict = xlib.get_submission_process_dict()

    # build the log dictionary
    log_dict = {}
    log_file_list = [file for file in os.listdir(xlib.get_log_dir()) if os.path.isfile(os.path.join(xlib.get_log_dir(), file)) and file.startswith(xconfiguration.environment)]
    for log_file in log_file_list:
        try:
            pattern = r'^(.+)\-(.+)\-(.+)\-(.+).txt$'
            mo = re.search(pattern, log_file)
            environment = mo.group(1).strip()
            submission_process_id = mo.group(2).strip()
            yymmdd = mo.group(3)
            hhmmss = mo.group(4)
            submission_process_text = submission_process_dict[submission_process_id]['text']
            date = '20{0}-{1}-{2}'.format(yymmdd[:2], yymmdd[2:4], yymmdd[4:])
            time = '{0}:{1}:{2}'.format(hhmmss[:2], hhmmss[2:4], hhmmss[4:])
        except Exception as e:
            submission_process_text = 'unknown process'
            date = '0000-00-00'
            time = '00:00:00'
        key = '{0}-{1}'.format(submission_process_text, log_file)
        log_dict[key] = {'submission_process_text': submission_process_text, 'log_file': log_file, 'date': date, 'time': time}

    # print the submission log list
    print(xlib.get_separator())
    if log_dict == {}:
        print('*** WARNING: There is not any submission log.')
    else:
        # set data width
        submission_process_text_width = 55
        log_file_width = 50
        date_width = 10
        time_width = 8
        # set line template
        line_template = '{0:' + str(submission_process_text_width) + '}   {1:' + str(log_file_width) + '}   {2:' + str(date_width) + '}   {3:' + str(time_width) + '}'
        # print header
        print(line_template.format('Process', 'Log file', 'Date', 'Time'))
        print(line_template.format('=' * submission_process_text_width, '=' * log_file_width, '=' * date_width, '=' * time_width))
        # print detail lines
        for key in sorted(log_dict.keys()):
            print(line_template.format(log_dict[key]['submission_process_text'], log_dict[key]['log_file'], log_dict[key]['date'], log_dict[key]['time']))

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_submission_log():
    '''
    View the log of a submission.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Logs - View a submission log')

    # get the submission log
    submission_log_file = cinputs.input_submission_log_file()
    if submission_log_file == '':
        print('WARNING: There is not any submission log.')
        OK = False
    
    # view the log file
    if OK:
        text = 'Logs - View a submission log'
        OK = clib.view_file(os.path.join(xlib.get_log_dir(), submission_log_file), text)

    # show continuation message 
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_cluster_experiment_processes():
    '''
    List the processes of an experiment in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Logs - List experiment processes in the cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster {0} does not have experiment data.'.format(cluster_name))
            OK = False

    # get the result dataset list of the experiment
    if OK:
        command = 'cd  {0}/{1}; for list in `ls`; do ls -ld $list | grep -v ^- > /dev/null && echo $list; done;'.format(xlib.get_cluster_result_dir(), experiment_id)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            result_dataset_id_list = []
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    result_dataset_id_list.append(line)

    # print the result dataset identification list of the experiment
    if OK:
        print(xlib.get_separator())
        if result_dataset_id_list == []:
            print('*** WARNING: There is not any result dataset of the experiment {0}.'.format(experiment_id))
        else:
            result_dataset_id_list.sort()
            # set data width
            result_dataset_width = 30
            bioinfo_app_width = 25
            # set line template
            line_template = '{0:' + str(result_dataset_width) + '}   {1:' + str(bioinfo_app_width) + '}'
            # print header
            print(line_template.format('Result dataset', 'Bioinfo app / Utility'))
            print(line_template.format('=' * result_dataset_width, '=' * bioinfo_app_width))
            # print detail lines
            for result_dataset_id in result_dataset_id_list:
                if result_dataset_id.startswith(xlib.get_bedtools_code()+'-'):
                    bioinfo_app_name = xlib.get_bedtools_name()
                elif result_dataset_id.startswith(xlib.get_blastplus_code()+'-'):
                    bioinfo_app_name = xlib.get_blastplus_name()
                elif result_dataset_id.startswith(xlib.get_bowtie2_code()+'-'):
                    bioinfo_app_name = xlib.get_bowtie2_name()
                elif result_dataset_id.startswith(xlib.get_busco_code()+'-'):
                    bioinfo_app_name = xlib.get_busco_name()
                elif result_dataset_id.startswith(xlib.get_cd_hit_code()+'-'):
                    bioinfo_app_name = xlib.get_cd_hit_name()
                elif result_dataset_id.startswith(xlib.get_cd_hit_est_code()+'-'):
                    bioinfo_app_name = xlib.get_cd_hit_est_name()
                elif result_dataset_id.startswith(xlib.get_cuffdiff_code()+'-'):
                    bioinfo_app_name = xlib.get_cuffdiff_name()
                elif result_dataset_id.startswith(xlib.get_cufflinks_code()+'-'):
                    bioinfo_app_name = xlib.get_cufflinks_name()
                elif result_dataset_id.startswith(xlib.get_cufflinks_cuffmerge_code()+'-'):
                    bioinfo_app_name = xlib.get_cufflinks_cuffmerge_name()
                elif result_dataset_id.startswith(xlib.get_cuffquant_code()+'-'):
                    bioinfo_app_name = xlib.get_cuffquant_name()
                elif result_dataset_id.startswith(xlib.get_cutadapt_code()+'-'):
                    bioinfo_app_name = xlib.get_cutadapt_name()
                elif result_dataset_id.startswith(xlib.get_ddradseq_simulation_code()+'-'):
                    bioinfo_app_name = xlib.get_ddradseq_simulation_name()
                elif result_dataset_id.startswith(xlib.get_ddradseqtools_code()+'-'):
                    bioinfo_app_name = xlib.get_ddradseqtools_name()
                elif result_dataset_id.startswith(xlib.get_detonate_code()+'-'):
                    bioinfo_app_name = xlib.get_detonate_name()
                elif result_dataset_id.startswith(xlib.get_emboss_code()+'-'):
                    bioinfo_app_name = xlib.get_emboss_name()
                elif result_dataset_id.startswith(xlib.get_entrez_direct_code()+'-'):
                    bioinfo_app_name = xlib.get_entrez_direct_name()
                elif result_dataset_id.startswith(xlib.get_fastqc_code()+'-'):
                    bioinfo_app_name = xlib.get_fastqc_name()
                elif result_dataset_id.startswith(xlib.get_ggtrinity_code()+'-'):
                    bioinfo_app_name = xlib.get_ggtrinity_name()
                elif result_dataset_id.startswith(xlib.get_gmap_gsnap_code()+'-'):
                    bioinfo_app_name = xlib.get_gmap_gsnap_name()
                elif result_dataset_id.startswith(xlib.get_gmap_code()+'-'):
                    bioinfo_app_name = xlib.get_gmap_name()
                elif result_dataset_id.startswith(xlib.get_gsnap_code()+'-'):
                    bioinfo_app_name = xlib.get_gsnap_name()
                elif result_dataset_id.startswith(xlib.get_gzip_code()+'-'):
                    bioinfo_app_name = xlib.get_gzip_name()
                elif result_dataset_id.startswith(xlib.get_hisat2_code()+'-'):
                    bioinfo_app_name = xlib.get_hisat2_name()
                elif result_dataset_id.startswith(xlib.get_htseq_code()+'-'):
                    bioinfo_app_name = xlib.get_htseq_name()
                elif result_dataset_id.startswith(xlib.get_htseq_count_code()+'-'):
                    bioinfo_app_name = xlib.get_htseq_count_name()
                elif result_dataset_id.startswith(xlib.get_insilico_read_normalization_code()+'-'):
                    bioinfo_app_name = xlib.get_insilico_read_normalization_name()
                elif result_dataset_id.startswith(xlib.get_ipyrad_code()+'-'):
                    bioinfo_app_name = xlib.get_ipyrad_name()
                elif result_dataset_id.startswith(xlib.get_kallisto_code()+'-'):
                    bioinfo_app_name = xlib.get_kallisto_name()
                elif result_dataset_id.startswith(xlib.get_miniconda3_code()+'-'):
                    bioinfo_app_name = xlib.get_miniconda3_name()
                elif result_dataset_id.startswith(xlib.get_ngshelper_code()+'-'):
                    bioinfo_app_name = xlib.get_ngshelper_name()
                elif result_dataset_id.startswith(xlib.get_quast_code()+'-'):
                    bioinfo_app_name = xlib.get_quast_name()
                elif result_dataset_id.startswith(xlib.get_r_code()+'-'):
                    bioinfo_app_name = xlib.get_r_name()
                elif result_dataset_id.startswith(xlib.get_ref_eval_code()+'-'):
                    bioinfo_app_name = xlib.get_ref_eval_name()
                elif result_dataset_id.startswith(xlib.get_rnaquast_code()+'-'):
                    bioinfo_app_name = xlib.get_rnaquast_name()
                elif result_dataset_id.startswith(xlib.get_rsem_code()+'-'):
                    bioinfo_app_name = xlib.get_rsem_name()
                elif result_dataset_id.startswith(xlib.get_rsem_eval_code()+'-'):
                    bioinfo_app_name = xlib.get_rsem_eval_name()
                elif result_dataset_id.startswith(xlib.get_rsitesearch_code()+'-'):
                    bioinfo_app_name = xlib.get_rsitesearch_name()
                elif result_dataset_id.startswith(xlib.get_samtools_code()+'-'):
                    bioinfo_app_name = xlib.get_samtools_name()
                elif result_dataset_id.startswith(xlib.get_soapdenovo2_code()+'-'):
                    bioinfo_app_name = xlib.get_soapdenovo2_name()
                elif result_dataset_id.startswith(xlib.get_soapdenovotrans_code()+'-'):
                    bioinfo_app_name = xlib.get_soapdenovotrans_name()
                elif result_dataset_id.startswith(xlib.get_star_code()+'-'):
                    bioinfo_app_name = xlib.get_star_name()
                elif result_dataset_id.startswith(xlib.get_starcode_code()+'-'):
                    bioinfo_app_name = xlib.get_starcode_name()
                elif result_dataset_id.startswith(xlib.get_toa_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_blastdb_nr_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_blastdb_nr_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_blastdb_nt_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_blastdb_nt_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_basic_data_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_basic_data_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_dicots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_dicots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_gene_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_gene_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_go_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_go_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_gymno_01_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_gymno_01_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_interpro_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_interpro_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_download_monocots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_download_monocots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_gilist_viridiplantae_protein_gi_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_gilist_viridiplantae_protein_gi_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_basic_data_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_basic_data_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_dicots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_dicots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_gene_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_gene_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_go_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_go_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_gymno_01_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_gymno_01_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_interpro_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_interpro_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_load_monocots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_load_monocots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_pipeline_aminoacid_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_pipeline_aminoacid_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_pipeline_nucleotide_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_pipeline_nucleotide_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_proteome_dicots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_proteome_dicots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_proteome_gymno_01_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_proteome_gymno_01_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_proteome_monocots_04_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_proteome_monocots_04_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_proteome_refseq_plant_code()+'-'):
                    bioinfo_app_name = xlib.get_toa_process_proteome_refseq_plant_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_rebuild_toa_database_code()+'-'):
                    bioinfo_app_name = xlib.get_get_toa_process_rebuild_toa_database_name()
                elif result_dataset_id.startswith(xlib.get_toa_process_recreate_toa_database_code()+'-'):
                    bioinfo_app_name = xlib.get_get_toa_process_recreate_toa_database_name()
                elif result_dataset_id.startswith(xlib.get_tophat_code()+'-'):
                    bioinfo_app_name = xlib.get_tophat_name()
                elif result_dataset_id.startswith(xlib.get_transabyss_code()+'-'):
                    bioinfo_app_name = xlib.get_transabyss_name()
                elif result_dataset_id.startswith(xlib.get_transcript_filter_code()+'-'):
                    bioinfo_app_name = xlib.get_transcript_filter_name()
                elif result_dataset_id.startswith(xlib.get_transcriptome_blastx_code()+'-'):
                    bioinfo_app_name = xlib.get_transcriptome_blastx_name()
                elif result_dataset_id.startswith(xlib.get_transdecoder_code()+'-'):
                    bioinfo_app_name = xlib.get_transdecoder_name()
                elif result_dataset_id.startswith(xlib.get_transrate_code()+'-'):
                    bioinfo_app_name = xlib.get_transrate_name()
                elif result_dataset_id.startswith(xlib.get_trimmomatic_code()+'-'):
                    bioinfo_app_name = xlib.get_trimmomatic_name()
                elif result_dataset_id.startswith(xlib.get_trinity_code()+'-'):
                    bioinfo_app_name = xlib.get_trinity_name()
                elif result_dataset_id.startswith(xlib.get_vsearch_code()+'-'):
                    bioinfo_app_name = xlib.get_vsearch_name()
                else:
                    bioinfo_app_name = 'xxx'
                print(line_template.format(result_dataset_id, bioinfo_app_name))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_cluster_experiment_process_log():
    '''
    View the log of an experiment process in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Logs - View an experiment process log in the cluster')

    # get the clustner name
    if OK:
        print(xlib.get_separator())
        if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
            cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
        else:
            print('WARNING: There is not any running cluster.')
            OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # create the SSH transport connection
    if OK:
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        for error in error_list:
            print(error)

    # create the SFTP client 
    if OK:
        sftp_client = xssh.create_sftp_client(ssh_transport)

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster does not have experiment data.')
            OK = False

    # get the result_dataset identification
    if OK:
        result_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'result', xlib.get_all_applications_selected_code(), 'uncompressed', help=True)
        if result_dataset_id == '':
            print('WARNING: The experiment {0} does not have result datasets.'.format(experiment_id))
            OK = False

    # create the local path
    if not os.path.exists(xlib.get_temp_dir()):
        os.makedirs(xlib.get_temp_dir())

    # get the log file name and build local and cluster paths
    if OK:
        log_file = xlib.get_cluster_log_file()
        local_path = '{0}/{1}'.format(xlib.get_temp_dir(), log_file)
        cluster_path = '{0}/{1}/{2}'.format(xlib.get_cluster_experiment_result_dir(experiment_id), result_dataset_id, log_file)

    # download the log file from the cluster
    if OK:
        print(xlib.get_separator())
        print('The file {0} is being downloaded from {1} ...'.format(log_file, cluster_path))
        OK = xssh.get_file(sftp_client, cluster_path, local_path)
        if OK:
            print('The file has been uploaded.')

    # close the SSH transport connection
    if OK:
        xssh.close_ssh_transport_connection(ssh_transport)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)
    
    # view the log file
    if OK:
        text = 'Logs - View an experiment process log in the cluster'
        OK = clib.view_file(local_path, text)

    # show continuation message 
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the functions related to forms corresponding to log menu items in mode console.')
    sys.exit(0)

#-------------------------------------------------------------------------------
