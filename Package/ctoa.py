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
This file contains the functions related to forms corresponding to TOA (Tree-oriented Annotation)
menu items in console mode.
'''

#-------------------------------------------------------------------------------

import gzip
import os
import subprocess
import sys

import cinputs
import clib
import xec2
import xlib
import xtoa
import xssh

#-------------------------------------------------------------------------------

def form_create_toa_config_file():
    '''
    Create the TOA config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Recreate config file')

    # set the TOA directory
    print(xlib.get_separator())
    toa_dir = f'{xlib.get_cluster_app_dir()}/TOA/Package'
    print(f'TOA directory: {toa_dir}')

    # set the Miniconda3 directory
    miniconda3_dir = f'{xlib.get_cluster_app_dir()}/Miniconda3'
    print(f'Miniconda bin directory: {miniconda3_dir}')

    # set the database directory
    db_dir = xlib.get_cluster_database_dir()
    print(f'Databases directory: {db_dir}')

    # set the result directory
    result_dir = '/result'
    print(f'Result directory: {result_dir}')

    # create the TOA config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The file {xtoa.get_toa_config_file()} is going to be recreated. The previous files will be lost.')

        # recreate the config file
        if OK:
            (OK, error_list) = xtoa.create_toa_config_file(toa_dir, miniconda3_dir, db_dir, result_dir)
            if OK:
                print('The file is recreated.')
            else:
                for error in error_list:
                    print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_manage_toa_database(process_type):
    '''
    Manage the TOA database.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    if process_type == xlib.get_toa_type_recreate():
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Recreate database')
    elif process_type == xlib.get_toa_type_rebuild():
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Rebuild database')

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
        OK = clib.confirm_action(f'The {xlib.get_toa_name()} database is going to be {process_type}.')

    # run the process
    if OK:
        devstdout = xlib.DevStdOut(xtoa.manage_toa_database.__name__)
        OK = xtoa.manage_toa_database(cluster_name, process_type, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_toa_config_file():
    '''
    List the TOA config file corresponding to the environment.
    '''

    # get the TOA config file
    toa_config_file = xtoa.get_toa_config_file()

    # view the file
    text = f'{xlib.get_toa_name()} - View config file'
    clib.view_file(toa_config_file, text)

    # show continuation message 
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_manage_genomic_database(process_type, genomic_database):
    '''
    Manage processes of genomic database.
    '''

    # initialize the control variable
    OK = True

    # set the genomica database name
    if genomic_database == xlib.get_toa_data_basic_data_code():
        name = xlib.get_toa_data_basic_data_name()
    elif genomic_database == xlib.get_toa_data_gymno_01_code():
        name = xlib.get_toa_data_gymno_01_name()
    elif genomic_database == xlib.get_toa_data_dicots_04_code():
        name = xlib.get_toa_data_dicots_04_name()
    elif genomic_database == xlib.get_toa_data_monocots_04_code():
        name = xlib.get_toa_data_monocots_04_name()
    elif genomic_database == xlib.get_toa_data_refseq_plant_code():
        name = xlib.get_toa_data_refseq_plant_name()
    elif genomic_database == xlib.get_toa_data_taxonomy_code():
        name = xlib.get_toa_data_taxonomy_name()
    elif genomic_database == xlib.get_toa_data_nt_code():
        name = xlib.get_toa_data_nt_name()
    elif genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
        name = xlib.get_toa_data_viridiplantae_nucleotide_gi_name()
    elif genomic_database == xlib.get_toa_data_nr_code():
        name = xlib.get_toa_data_nr_name()
    elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
        name = xlib.get_toa_data_viridiplantae_protein_gi_name()
    elif genomic_database == xlib.get_toa_data_gene_code():
        name = xlib.get_toa_data_gene_name()
    elif genomic_database == xlib.get_toa_data_interpro_code():
        name = xlib.get_toa_data_interpro_name()
    elif genomic_database == xlib.get_toa_data_go_code():
        name = xlib.get_toa_data_go_name()

    # print the header
    clib.clear_screen()
    if process_type == xlib.get_toa_type_build_blastplus_db():
        clib.print_headers_with_environment(f'Build {name} for BLAST+')
    elif process_type == xlib.get_toa_type_build_diamond_db():
        clib.print_headers_with_environment(f'Build {name} for DIAMOND')
    elif process_type == xlib.get_toa_type_build_gilist():
        clib.print_headers_with_environment(f'Build {name}')
    elif process_type == xlib.get_toa_type_build_proteome():
        clib.print_headers_with_environment(f'Build {name} proteome')
    elif process_type == xlib.get_toa_type_download_data():
        clib.print_headers_with_environment(f'Download {name} functional annotations')
    elif process_type == xlib.get_toa_type_load_data():
        clib.print_headers_with_environment(f'Load {name} data in {xlib.get_toa_name()} database')

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
        OK = clib.confirm_action(f'The {name} process is going to be run.')

    # run the process
    if OK:
        devstdout = xlib.DevStdOut(xtoa.manage_genomic_database.__name__)
        OK = xtoa.manage_genomic_database(cluster_name, process_type, genomic_database, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_pipeline_config_file(pipeline_type):
    '''
    Recreate a pipeline config file.
    '''

    # initialize the control variable
    OK = True

    # set the pipeline name
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        name = xlib.get_toa_process_pipeline_nucleotide_name()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        name = xlib.get_toa_process_pipeline_aminoacid_name()

    # set the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = xtoa.get_nucleotide_pipeline_config_file()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = xtoa.get_aminoacid_pipeline_config_file()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{name} - Recreate config file')

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

    # get the transcriptome origin
    assembly_origin = cinputs.input_assembly_origin(help=True)

    # get the experiment identification
    if OK:
        if assembly_origin == 'NGSCLOUD':
            experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
            if experiment_id == '':
                print(f'WARNING: The cluster {cluster_name} does not have experiment data.')
                OK = False
        elif assembly_origin == 'EXTERNAL':
            experiment_id = 'NONE'

    # get the assembly dataset identification
    if OK:
        if assembly_origin == 'NGSCLOUD':
            app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
            assembly_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'assembly', app_list, 'uncompressed', help=True)
            if assembly_dataset_id == '':
                print(f'WARNING: The cluster {cluster_name} does not have assembly datasets.')
                OK = False
        elif assembly_origin == 'EXTERNAL':
            assembly_dataset_id = 'NONE'

    # get the assembly type
    if OK:
        if assembly_origin == 'NGSCLOUD':
            if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
                assembly_type = cinputs.input_assembly_type(help=True)
            elif assembly_dataset_id.startswith(xlib.get_transabyss_code()) or assembly_dataset_id.startswith(xlib.get_trinity_code()) or assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
                assembly_type = 'NONE'
        elif assembly_origin == 'EXTERNAL':
            assembly_type = 'NONE'

    # get the reference dataset identification
    if OK:
        if assembly_origin == 'NGSCLOUD':
            reference_dataset_id = 'NONE'
        elif assembly_origin == 'EXTERNAL':
            reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
            if reference_dataset_id == '':
                print(f'*** ERROR: The cluster {cluster_name} does not have reference datasets.')
                OK = False

    # get the transcriptome file
    if OK:
        if assembly_origin == 'NGSCLOUD':
            transcriptome_file = 'NONE'
        elif assembly_origin == 'EXTERNAL':
            transcriptome_file = cinputs.input_transcriptome_file(ssh_client, reference_dataset_id, help=True)
            if transcriptome_file == '':
                print(f'ERROR: The reference dataset {reference_dataset_id} does not have any file.'.format())
                OK = False

    # get the database list
    if OK:

        # nucleotide pipelines
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            database_list = cinputs.input_database_list(xtoa.get_nucleotide_annotation_database_code_list(), 'nt')

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            database_list = cinputs.input_database_list(xtoa.get_aminoacid_annotation_database_code_list(), 'nr')

    # recreate the pipeline config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The file {config_file} is going to be recreated. The previous files will be lost.')

        # recreate the config file
        if OK:
            (OK, error_list) = xtoa.create_pipeline_config_file(pipeline_type, assembly_origin, experiment_id, assembly_dataset_id, assembly_type, reference_dataset_id, transcriptome_file, database_list)
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

def form_recreate_annotation_merger_config_file():
    '''
    Recreate the annotation merger config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{xlib.get_toa_process_merge_annotations_name()} - Recreate config file')

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
        experiment_id = xlib.get_toa_database_dir()

    # get the identification of the first pipeline dataset
    app_list = [xlib.get_all_applications_selected_code()]
    print('First pipeline ...')
    pipeline_dataset_id_1 = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(),  'pipeline', app_list, 'uncompressed', help=True)
    if pipeline_dataset_id_1 == '':
        print( 'WARNING: There are not any pipeline datasets.')
        OK = False

    # get the identification of the second pipeline dataset
    app_list = [xlib.get_all_applications_selected_code()]
    print('Second pipeline ...')
    pipeline_dataset_id_2 = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), 'pipeline', app_list, 'uncompressed', help=True)
    if pipeline_dataset_id_2 == '':
        print( 'WARNING: There are not any pipeline datasets.')
        OK = False
    elif pipeline_dataset_id_1 == pipeline_dataset_id_2:
        print( 'ERROR: The first pipeline dataset is equal to the second one.')
        OK = False

    # get the merger operation
    if OK:
        merger_operation = cinputs.input_code(text='Merger operation', code_list=xlib.get_annotation_merger_operation_code_list(), default_code=None).upper()

    # recreate the pipeline config file
    if OK:

        # confirm the creation of the config file
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The file {xtoa.get_annotation_merger_config_file()} is going to be recreated. The previous files will be lost.')

        # recreate the config file
        if OK:
            (OK, error_list) = xtoa.create_annotation_merger_config_file(experiment_id, pipeline_dataset_id_1, pipeline_dataset_id_2, merger_operation)
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

def form_edit_pipeline_config_file(pipeline_type):
    '''
    Edit a pipeline config file to change the parameters of each process.
    '''

    # initialize the control variable
    OK = True

    # set the pipeline name
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        name = xlib.get_toa_process_pipeline_nucleotide_name()

    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        name = xlib.get_toa_process_pipeline_aminoacid_name()

    elif pipeline_type == xlib.get_toa_process_merge_annotations_code():
        name = xlib.get_toa_process_merge_annotations_name()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{name} - Edit config file')

    # get the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = xtoa.get_nucleotide_pipeline_config_file()

    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = xtoa.get_aminoacid_pipeline_config_file()

    elif pipeline_type == xlib.get_toa_process_merge_annotations_code():
        config_file = xtoa.get_annotation_merger_config_file()

    # edit the read transfer config file
    print(xlib.get_separator())
    print(f'Editing the {name} config file ...')
    command = f'{xlib.get_editor()} {config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the config file
    if OK:
        print(xlib.get_separator())
        print(f'Checking the {name} config file ...')

        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            (OK, error_list) = xtoa.check_pipeline_config_file(pipeline_type, strict=False)

        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            (OK, error_list) = xtoa.check_pipeline_config_file(pipeline_type, strict=False)

        elif pipeline_type == xlib.get_toa_process_merge_annotations_code():
            (OK, error_list) = xtoa.check_annotation_merger_config_file(strict=False)

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

def form_run_pipeline_process(pipeline_type):
    '''
    Run a pipeline process with the parameters in the corresponding config file.
    '''

    # initialize the control variable
    OK = True


    # set the pipeline name
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        name = xlib.get_toa_process_pipeline_nucleotide_name()

    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        name = xlib.get_toa_process_pipeline_aminoacid_name()

    elif pipeline_type == xlib.get_toa_process_merge_annotations_code():
        name = xlib.get_toa_process_merge_annotations_name()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{name} - Run process')

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
        OK = clib.confirm_action(f'The {name} process is going to be run.')

    # run the process
    if OK:

        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            devstdout = xlib.DevStdOut(xtoa.run_pipeline_process.__name__)
            OK = xtoa.run_pipeline_process(cluster_name, pipeline_type, devstdout, function=None)

        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            devstdout = xlib.DevStdOut(xtoa.run_pipeline_process.__name__)
            OK = xtoa.run_pipeline_process(cluster_name, pipeline_type, devstdout, function=None)

        elif pipeline_type == xlib.get_toa_process_merge_annotations_code():
            devstdout = xlib.DevStdOut(xtoa.run_annotation_merger_process.__name__)
            OK = xtoa.run_annotation_merger_process(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_restart_pipeline_process(pipeline_type):
    '''
    Restart a pipeline process from the last step ended OK.
    '''

    # initialize the control variable
    OK = True

    # set the pipeline name
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        name = xlib.get_toa_process_pipeline_nucleotide_name()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        name = xlib.get_toa_process_pipeline_aminoacid_name()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'{name} - Restart process')

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
        experiment_id = xlib.get_toa_result_pipeline_dir()
        print(f'... Enter the experiment id: {experiment_id} ...')

    # get the pipeline dataset identification
    if OK:
        app_list = [pipeline_type]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'pipeline', app_list, 'uncompressed', help=True)
        if pipeline_dataset_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have result datasets.')
            OK = False

    # confirm the process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The {name} process is going to be run.')

    # run the process
    if OK:

        devstdout = xlib.DevStdOut(xtoa.restart_pipeline_process.__name__)
        OK = xtoa.restart_pipeline_process(cluster_name, experiment_id, pipeline_type, pipeline_dataset_id, devstdout, function=None)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_x_per_y_data(stats_code):
    '''
    View the x per y data.
    '''

    # initialize the control variable
    OK = True

    # assign the text of the "name"
    if stats_code == 'hit_per_hsp':
        name = '# HITs per # HSPs'
    elif stats_code == 'seq_per_go':
        name = '# sequences per # GO terms'
    elif stats_code == 'seq_per_ec':
        name = '# sequences per # EC ids'
    elif stats_code == 'seq_per_interpro':
        name = '# sequences per # InterPro ids'
    elif stats_code == 'seq_per_kegg':
        name = '# sequences per # KEGG ids'
    elif stats_code == 'seq_per_mapman':
        name = '# sequences per # MapMan ids'
    elif stats_code == 'seq_per_metacyc':
        name = '# sequences per # MetaCyc ids'

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'Statistics - {name} data')

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

    # get the pipeline dataset identification
    if OK:
        app_list = [xlib.get_all_applications_selected_code()]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), '', app_list, status='uncompressed', help=True)
        if pipeline_dataset_id == '':
            print('WARNING: There are not any annotation pipeline result datasets.')
            OK = False

    # build distribution dictionary
    if OK:

        # initialize the distribution dictionary
        distribution_dict = {}

        # get the dictionary of TOA configuration
        toa_config_dict = xtoa.get_toa_config_dict()

        # get the statistics file path
        stats_file = f'{toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id}/{toa_config_dict["STATS_SUBDIR_NAME"]}/{stats_code}-{toa_config_dict["STATS_BASE_NAME"]}.csv'

        # open the statistics file
        if stats_file.endswith('.gz'):
            try:
                stats_file_id = gzip.open(stats_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F002', stats_file)
        else:
            try:
                stats_file_id = open(stats_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F001', stats_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = stats_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: "x_count";"y_count"
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end].strip('"'))
                    begin = end + 1
                data_list.append(record[begin:].strip('\n').strip('"'))
                try:
                    x_count = data_list[0]
                    y_count = data_list[1]
                except Exception as e:
                    print(f'*** EXCEPTION: "{e}".')
                    raise xlib.ProgramException('F006', os.path.basename(stats_file), record_counter)

                # add dato to the dictionary
                distribution_dict[record_counter] = {'x_count': x_count, 'y_count': y_count}

            # read the next record
            record = stats_file_id.readline()

    # print the distribution
    if OK:
        print(xlib.get_separator())
        if distribution_dict == {}:
            print('*** WARNING: There is not any stats data.')
        else:
            # set data width
            x_count_width = 15
            y_count_width = 15
            # set line
            line = '{0:' + str(x_count_width) + '}   {1:' + str(y_count_width) + '}'
            # print header
            if stats_code == 'hit_per_hsp':
                print(line.format('# HSPs', '# HITs'))
            elif stats_code == 'seq_per_go':
                print(line.format('GO terms', '# sequences'))
            elif stats_code == 'seq_per_ec':
                print(line.format('# EC ids', '# sequences'))
            elif stats_code == 'seq_per_interpro':
                print(line.format('# InterPro ids', '# sequences'))
            elif stats_code == 'seq_per_kegg':
                print(line.format('# KEGG ids', '# sequences'))
            elif stats_code == 'seq_per_mapman':
                print(line.format('# MapMan ids', '# sequences'))
            elif stats_code == 'seq_per_metacyc':
                print(line.format('# MetaCyc ids', '# sequences'))
            print(line.format('=' * x_count_width, '=' * y_count_width))
            # print detail lines
            for key in sorted(distribution_dict.keys()):
                print(line.format(distribution_dict[key]['x_count'], distribution_dict[key]['y_count']))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_dataset_data_frecuency():
    '''
    View the frecuency distribution of annotation dataset data.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Statistics - Annotation datasets - Frequency distribution data')

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

    # get the pipeline dataset identification
    if OK:
        app_list = [xlib.get_all_applications_selected_code()]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), '', app_list, status='uncompressed', help=True)
        if pipeline_dataset_id == '':
            print('WARNING: There are not any annotation pipeline result datasets.')
            OK = False

    # build distribution dictionary
    if OK:

        # initialize the distribution dictionary
        distribution_dict = {}

        # get the dictionary of TOA configuration
        toa_config_dict = xtoa.get_toa_config_dict()

        # get the statistics file path
        stats_file = f'{toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id}/{toa_config_dict["STATS_SUBDIR_NAME"]}/dataset-{toa_config_dict["STATS_BASE_NAME"]}.csv'

        # open the statistics file
        if stats_file.endswith('.gz'):
            try:
                stats_file_id = gzip.open(stats_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F002', stats_file)
        else:
            try:
                stats_file_id = open(stats_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                raise xlib.ProgramException('F001', stats_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = stats_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: "dataset_name";"annotated_seq_count";"remained_seq_count"
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end].strip('"'))
                    begin = end + 1
                data_list.append(record[begin:].strip('\n').strip('"'))
                try:
                    dataset_name = data_list[0]
                    annotated_seq_count = data_list[1]
                    remained_seq_count = data_list[2]
                except Exception as e:
                    raise xlib.ProgramException('F006', os.path.basename(stats_file), record_counter)

                # add dato to the dictionary
                distribution_dict[record_counter] = {'dataset_name': dataset_name, 'annotated_seq_count': annotated_seq_count, 'remained_seq_count': remained_seq_count}

            # read the next record
            record = stats_file_id.readline()

    # print the distribution
    if OK:
        print(xlib.get_separator())
        if distribution_dict == {}:
            print('*** WARNING: There is not any distribution.')
        else:
            # set data width
            dataset_name_width = 19
            annotated_seq_count_width = 14
            remained_seq_count_width = 14
            # set line
            line = '{0:' + str(dataset_name_width) + '}   {1:' + str(annotated_seq_count_width) + '}   {2:' + str(remained_seq_count_width) + '}'
            # print header
            print(line.format('Dataset', 'Annotated seqs', 'Remained seqs'))
            print(line.format('=' * dataset_name_width, '=' * annotated_seq_count_width, '=' * remained_seq_count_width))
            # print detail lines
            for key in sorted(distribution_dict.keys()):
                print(line.format(distribution_dict[key]['dataset_name'], distribution_dict[key]['annotated_seq_count'], distribution_dict[key]['remained_seq_count']))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_phylogenic_data_frecuency(stats_code):
    '''
    View the frecuency distribution of phylogenic data.
    '''

    # initialize the control variable
    OK = True

    # assign the text of the "name"
    if stats_code == 'species':
        name = 'Species - Frequency distribution'
    elif stats_code == 'family':
        name = 'Family - Frequency distribution'
    elif stats_code == 'phylum':
        name = 'Phylum - Frequency distribution'
    elif stats_code == 'namespace':
        name = 'GO - Frequency distribution per namespace'

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'Statistics - {name} data')

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

    # get the pipeline dataset identification
    if OK:
        app_list = [xlib.get_all_applications_selected_code()]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), '', app_list, status='uncompressed', help=True)
        if pipeline_dataset_id == '':
            print('WARNING: There are not any annotation pipeline result datasets.')
            OK = False

    # build distribution dictionary
    if OK:

        # initialize the distribution dictionary
        distribution_dict = {}

        # get the dictionary of TOA configuration
        toa_config_dict = xtoa.get_toa_config_dict()

        # get the statistics file path
        stats_file = f'{toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id}/{toa_config_dict["STATS_SUBDIR_NAME"]}/{stats_code}-{toa_config_dict["STATS_BASE_NAME"]}.csv'

        # open the statistics file
        if stats_file.endswith('.gz'):
            try:
                stats_file_id = gzip.open(stats_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F002', stats_file)
        else:
            try:
                stats_file_id = open(stats_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F001', stats_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = stats_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: "stats_code_id";"all_count";"first_hsp_count";"min_evalue_count"
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end].strip('"'))
                    begin = end + 1
                data_list.append(record[begin:].strip('\n').strip('"'))
                try:
                    id = data_list[0]
                    all_count = data_list[1]
                    first_hsp_count = data_list[2]
                    min_evalue_count = data_list[3]
                except Exception as e:
                    raise xlib.ProgramException('F006', os.path.basename(stats_file), record_counter)

                # add dato to the dictionary
                distribution_dict[id] = {'id': id, 'all_count': all_count, 'first_hsp_count': first_hsp_count, 'min_evalue_count': min_evalue_count}

            # read the next record
            record = stats_file_id.readline()

    # print the distribution
    if OK:
        print(xlib.get_separator())
        if distribution_dict == {}:
            print('*** WARNING: There is not any distribution.')
        else:
            # set data width
            id_width = 50
            all_count_width = 11
            first_hsp_count_width = 11
            min_evalue_count_width = 11
            # set line
            line = '{0:' + str(id_width) + '}   {1:' + str(all_count_width) + '}   {2:' + str(first_hsp_count_width) + '}   {3:' + str(min_evalue_count_width) + '}'
            # print header
            print(line.format(stats_code.capitalize(), 'All', 'First HSP', 'Min e-value'))
            print(line.format('=' * id_width, '=' * all_count_width, '=' * first_hsp_count_width, '=' * min_evalue_count_width))
            # print detail lines
            for key in sorted(distribution_dict.keys()):
                print(line.format(distribution_dict[key]['id'], distribution_dict[key]['all_count'], distribution_dict[key]['first_hsp_count'], distribution_dict[key]['min_evalue_count']))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_ontologic_data_frecuency(stats_code):
    '''
    View the frecuency distribution of ontologic data.
    '''

    # initialize the control variable
    OK = True

    # assign the text of the "name"
    if stats_code == 'go':
        name = 'Gene Ontology - Frequency distribution'
    elif stats_code == 'ec':
        name = 'EC - Frequency distribution'
    elif stats_code == 'interpro':
        name = 'InterPro - Frequency distribution'
    elif stats_code == 'kegg':
        name = 'KEGG - Frequency distribution'
    elif stats_code == 'mapman':
        name = 'MapMan - Frequency distribution'
    elif stats_code == 'metacyc':
        name = 'MetaCyc - Frequency distribution'

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment(f'Statistics - {name} data')

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

    # get the pipeline dataset identification
    if OK:
        app_list = [xlib.get_all_applications_selected_code()]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), '', app_list, status='uncompressed', help=True)
        if pipeline_dataset_id == '':
            print('WARNING: There are not any annotation pipeline result datasets.')
            OK = False

    # build distribution dictionary
    if OK:

        # initialize the distribution dictionary
        distribution_dict = {}

        # get the dictionary of TOA configuration
        toa_config_dict = xtoa.get_toa_config_dict()

        # get the statistics file path
        stats_file = f'{toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id}/{toa_config_dict["STATS_SUBDIR_NAME"]}/{stats_code}-{toa_config_dict["STATS_BASE_NAME"]}.csv'

        # open the statistics file
        if stats_file.endswith('.gz'):
            try:
                stats_file_id = gzip.open(stats_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F002', stats_file)
        else:
            try:
                stats_file_id = open(stats_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F001', stats_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = stats_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: "stats_code_id";"description";"all_count";"first_hsp_count";"min_evalue_count"
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end].strip('"'))
                    begin = end + 1
                data_list.append(record[begin:].strip('\n').strip('"'))
                try:
                    id = data_list[0]
                    desc = data_list[1]
                    all_count = data_list[2]
                    first_hsp_count = data_list[3]
                    min_evalue_count = data_list[4]
                except Exception as e:
                    raise xlib.ProgramException('F006', os.path.basename(stats_file), record_counter)

                # add dato to the dictionary
                distribution_dict[id] = {'id': id, 'desc': desc, 'all_count': all_count, 'first_hsp_count': first_hsp_count, 'min_evalue_count': min_evalue_count}

            # read the next record
            record = stats_file_id.readline()

    # print the distribution
    if OK:
        print(xlib.get_separator())
        if distribution_dict == {}:
            print('*** WARNING: There is not any distribution.')
        else:
            # set data width
            id_width = 30
            desc_width = 70
            all_count_width = 11
            first_hsp_count_width = 11
            min_evalue_count_width = 11
            # set line
            line = '{0:' + str(id_width) + '}   {1:' + str(desc_width) + '}   {2:' + str(all_count_width) + '}   {3:' + str(first_hsp_count_width) + '}   {4:' + str(min_evalue_count_width) + '}'
            # print header
            print(line.format(f'{stats_code.capitalize()} id', 'Description', 'All', 'First HSP', 'Min e-value'))
            print(line.format('=' * id_width, '=' * desc_width, '=' * all_count_width, '=' * first_hsp_count_width, '=' * min_evalue_count_width))
            # print detail lines
            for key in sorted(distribution_dict.keys()):
                print(line.format(distribution_dict[key]['id'], distribution_dict[key]['desc'], distribution_dict[key]['all_count'], distribution_dict[key]['first_hsp_count'], distribution_dict[key]['min_evalue_count']))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_go_data_frecuency():
    '''
    View the frecuency distribution of Gene Ontology data.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Statistics - Gene Ontology - Frequency distribution data')

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

    # get the pipeline dataset identification
    if OK:
        app_list = [xlib.get_all_applications_selected_code()]
        pipeline_dataset_id = cinputs.input_result_dataset_id(ssh_client, xlib.get_toa_result_pipeline_dir(), '', app_list, status='uncompressed', help=True)
        if pipeline_dataset_id == '':
            print('WARNING: There are not any annotation pipeline result datasets.')
            OK = False

    # build distribution dictionary
    if OK:

        # initialize the distribution dictionary
        distribution_dict = {}

        # get the dictionary of TOA configuration
        toa_config_dict = xtoa.get_toa_config_dict()

        # get the statistics file path
        stats_file = f'{toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id}/{toa_config_dict["STATS_SUBDIR_NAME"]}/go-{toa_config_dict["STATS_BASE_NAME"]}.csv'

        # open the statistics file
        if stats_file.endswith('.gz'):
            try:
                stats_file_id = gzip.open(stats_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F002', stats_file)
        else:
            try:
                stats_file_id = open(stats_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                print(f'*** EXCEPTION: "{e}".')
                raise xlib.ProgramException('F001', stats_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = stats_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: "go_id";"description";"namespace";"all_count";"first_hsp_count";"min_evalue_count"
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end].strip('"'))
                    begin = end + 1
                data_list.append(record[begin:].strip('\n').strip('"'))
                try:
                    id = data_list[0]
                    desc = data_list[1]
                    namespace = data_list[2]
                    all_count = data_list[3]
                    first_hsp_count = data_list[4]
                    min_evalue_count = data_list[5]
                except Exception as e:
                    raise xlib.ProgramException('F006', os.path.basename(stats_file), record_counter)

                # add dato to the dictionary
                distribution_dict[id] = {'id': id, 'desc': desc, 'namespace': namespace, 'all_count': all_count, 'first_hsp_count': first_hsp_count, 'min_evalue_count': min_evalue_count}

            # read the next record
            record = stats_file_id.readline()

    # print the distribution
    if OK:
        print(xlib.get_separator())
        if distribution_dict == {}:
            print('*** WARNING: There is not any distribution.')
        else:
            # set data width
            id_width = 10
            desc_width = 50
            namespace_width = 18
            all_count_width = 11
            first_hsp_count_width = 11
            min_evalue_count_width = 11
            # set line
            line = '{0:' + str(id_width) + '}   {1:' + str(desc_width) + '}   {2:' + str(namespace_width) + '}   {3:' + str(all_count_width) + '}   {4:' + str(first_hsp_count_width) + '}   {5:' + str(min_evalue_count_width) + '}'
            # print header
            print(line.format('GO id', 'Description', 'Namespace', 'All', 'First HSP', 'Min e-value'))
            print(line.format('=' * id_width, '=' * desc_width, '=' * namespace_width, '=' * all_count_width, '=' * first_hsp_count_width, '=' * min_evalue_count_width))
            # print detail lines
            for key in sorted(distribution_dict.keys()):
                print(line.format(distribution_dict[key]['id'], distribution_dict[key]['desc'], distribution_dict[key]['namespace'], distribution_dict[key]['all_count'], distribution_dict[key]['first_hsp_count'], distribution_dict[key]['min_evalue_count']))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the functions related to forms corresponding to TOA (Tree-oriented Annotation) menu items in console mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
