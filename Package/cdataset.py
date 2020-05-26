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
import subprocess
import sys

import cinputs
import clib
import xconfiguration
import xdatabase
import xec2
import xgzip
import xlib
import xread
import xreference
import xresult
import xssh

#-------------------------------------------------------------------------------

def form_recreate_reference_transfer_config_file():
    '''
    Recreate the reference transfer config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset transfer - Recreate config file')

    # get the local directory, cluster directory and file pattern
    print(xlib.get_separator())
    local_dir = cinputs.input_local_dir('reference')
    cluster_directoy = cinputs.input_cluster_directory('reference')
    file_pattern = cinputs.input_files_pattern('.*')

    # build the selected file list
    print(xlib.get_separator())
    print('Building the selected file list ...')
    selected_file_list = []
    for file in os.listdir(local_dir):
        if os.path.isfile(os.path.join(local_dir, file)) and re.match(file_pattern, file):
            selected_file_list.append(file)
    if selected_file_list != []:
        selected_file_list = str(selected_file_list).strip('[]').replace('\'','')
        print(f'The selected file list is built: {selected_file_list}')
    else:
        print(f'*** ERROR: There are not files in the directory {local_dir} with the pattern {file_pattern}')
        OK = False

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The reference transfer config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {xreference.get_reference_transfer_config_file()} is been created ...')
        (OK, error_list) = xreference.create_reference_transfer_config_file(local_dir, selected_file_list, cluster_directoy)
        if OK:
            print('The config file is created.')
        else:
            print()
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_reference_transfer_config_file():
    '''
    Edit reference transfer config file to change the parameters of each transfer.
    '''

    # initialize the control variable
    OK = True

    # print the head
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset transfer - Edit config file')

    # get the reference transfer config file
    reference_transfer_config_file = xreference.get_reference_transfer_config_file()

    # edit the reference transfer config file
    print(xlib.get_separator())
    print('The reference transfer config file is being edited ...')
    command = f'{xlib.get_editor()} {reference_transfer_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the reference transfer config file
    if OK:
        print(xlib.get_separator())
        print('Checking the reference transfer config file ...')
        (OK, error_list) = xreference.check_reference_transfer_config_file(strict=False)
        if OK:
            print('The file is OK.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_upload_reference_dataset():
    '''
    Upload reference dataset to a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset transfer - Upload dataset to a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the upload of the reference dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The reference dataset is going to be uploaded.')

    # upload the reference dataset to the cluster
    if OK:
        devstdout = xlib.DevStdOut(xreference.upload_reference_dataset.__name__)
        OK = xreference.upload_reference_dataset(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_reference_gzip_config_file():
    '''
    Recreate reference file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset file compression/decompression - Recreate config file')

    # get the compress/decompress config file
    gzip_config_file = xgzip.get_gzip_config_file('reference')

    # get the cluster name
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

    # get the action
    if OK:
        action = cinputs.input_code(text='Action', code_list=['compress', 'decompress'], default_code=None).lower()

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id == '':
            print('WARNING: There are not reference datasets.')
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*' if action == 'compress' else '.*gz')

    # get the selected file list
    if OK:
        print(xlib.get_separator())
        print('Building the selected file list ...')
        reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)
        selected_file_list = []
        command = f'cd {reference_dataset_dir}; find . -type f -regex "./{file_pattern}"'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
            if selected_file_list != []:
                selected_file_list = str(selected_file_list).strip('[]').replace('\'','')
                print(f'The selected file list is built: {selected_file_list}')
            else:
                print(f'*** ERROR: There are not files in the directory {reference_dataset_dir} with the pattern {file_pattern}')
                OK = False
        else:
            print(f'*** ERROR: Wrong command ---> {command}')

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The reference compression/decompression config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {gzip_config_file} is been created ...')
        (OK, error_list) = xgzip.create_gzip_config_file(action, 'reference', None, reference_dataset_id, selected_file_list)
        if OK:
            print('The config file is created.')
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

def form_edit_reference_gzip_config_file():
    '''
    Edit reference file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset file compression/decompression - Edit config file')

    # get the reference gzip config file
    reference_gzip_config_file = xgzip.get_gzip_config_file('reference')

    # edit the reference gzip config file
    print(xlib.get_separator())
    print('The reference gzip config file is being edited ...')
    command = f'{xlib.get_editor()} {reference_gzip_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the reference gzip config file
    if OK:
        print(xlib.get_separator())
        print('Checking the reference gzip config file ...')
        (OK, error_list) = xgzip.check_gzip_config_file('reference', strict=False)
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

def form_run_reference_gzip_process():
    '''
    Compress/decompress reference dataset files in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Reference dataset file compression/decompression - Run process')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the compression/decompression process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The reference files are going to be compressed/decompressed.')

    # compress/decompress reference dataset files in the cluster
    if OK:
        devstdout = xlib.DevStdOut(xgzip.run_gzip_process.__name__)
        OK = xgzip.run_gzip_process(cluster_name, 'reference', devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_reference_dataset():
    '''
    Remove reference dataset in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Datasets - Remove reference dataset')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the reference dataset identification
    if OK:
        reference_dataset_id = cinputs.input_reference_dataset_id(ssh_client, allowed_none=False, help=True)
        if reference_dataset_id  == '':
            print('WARNING: There are not any reference datasets.')
            OK = False

    # confirm the removal of the reference dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The reference dataset is going to be removed.')

    # removal the reference dataset in the cluster
    if OK:
        command = f'rm -fr {xlib.get_cluster_reference_dir()}/{reference_dataset_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            print(f'The dataset {reference_dataset_id} is removed.')
        else:
            print(f'*** ERROR: Wrong command ---> {command}')
            OK = False

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_database_transfer_config_file():
    '''
    Recreate the database transfer config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Database file transfer - Recreate config file')

    # get the local directory, cluster directory and file pattern
    print(xlib.get_separator())
    local_dir = cinputs.input_local_dir('database')
    cluster_directoy = cinputs.input_cluster_directory('database')
    file_pattern = cinputs.input_files_pattern('.*')

    # build the selected file list
    print(xlib.get_separator())
    print('Building the selected file list ...')
    selected_file_list = []
    for file in os.listdir(local_dir):
        if os.path.isfile(os.path.join(local_dir, file)) and re.match(file_pattern, file):
            selected_file_list.append(file)
    if selected_file_list != []:
        selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
        print(f'The selected file list is built: {selected_file_list_text}')
    else:
        print(f'*** ERROR: There are not files in the directory {local_dir} with the pattern {file_pattern}')
        OK = False

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The database transfer config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {xdatabase.get_database_transfer_config_file()} is been created ...')
        (OK, error_list) = xdatabase.create_database_transfer_config_file(local_dir, selected_file_list, cluster_directoy)
        if OK:
            print('The config file is created.')
        else:
            print()
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_database_transfer_config_file():
    '''
    Edit database transfer config file to change the parameters of each transfer.
    '''

    # initialize the control variable
    OK = True

    # print the head
    clib.clear_screen()
    clib.print_headers_with_environment('Database file transfer - Edit config file')

    # edit the database transfer config file
    print(xlib.get_separator())
    print('The database transfer config file is being edited ...')
    command = f'{xlib.get_editor()} {xdatabase.get_database_transfer_config_file()}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the database transfer config file
    if OK:
        print(xlib.get_separator())
        print('Checking the database transfer config file ...')
        (OK, error_list) = xdatabase.check_database_transfer_config_file(strict=False)
        if OK:
            print('The file is OK.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_upload_database_dataset():
    '''
    Upload database to a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Database file transfer - Upload dataset to a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the upload of the database dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The database files are going to be uploaded.')

    # upload the database dataset to the cluster
    if OK:
        devstdout = xlib.DevStdOut(xdatabase.upload_database_dataset.__name__)
        OK = xdatabase.upload_database_dataset(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_database_gzip_config_file():
    '''
    Recreate database file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Database file compression/decompression - Recreate config file')

    # get the cluster name
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

    # get the action
    if OK:
        action = cinputs.input_code(text='Action', code_list=['compress', 'decompress'], default_code=None).lower()

    # get the database dataset identification
    if OK:
        database_dataset_id = cinputs.input_database_dataset_id(ssh_client, help=True)
        if database_dataset_id == '':
            print('WARNING: There are not databases.')
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*' if action == 'compress' else '.*gz')

    # get the selected file list
    if OK:
        print(xlib.get_separator())
        print('Building the selected file list ...')
        database_dataset_dir = xlib.get_cluster_database_dataset_dir(database_dataset_id)
        selected_file_list = []
        command = f'cd {database_dataset_dir}; find . -type f -regex "./{file_pattern}"'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
            if selected_file_list != []:
                selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
                print(f'The selected file list is built: {selected_file_list_text}')
            else:
                print(f'*** ERROR: There are not files in the directory {database_dataset_dir} with the pattern {file_pattern}')
                OK = False
        else:
            print(f'*** ERROR: Wrong command ---> {command}')

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The database compression/decompression config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {xgzip.get_gzip_config_file("database")} is been created ...')
        (OK, error_list) = xgzip.create_gzip_config_file(action, 'database', None, database_dataset_id, selected_file_list)
        if OK:
            print('The config file is created.')
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

def form_edit_database_gzip_config_file():
    '''
    Edit database file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Database file compression/decompression - Edit config file')

    # edit the database gzip config file
    print(xlib.get_separator())
    print('The database gzip config file is being edited ...')
    command = f'{xlib.get_editor()} {xgzip.get_gzip_config_file("database")}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the database gzip config file
    if OK:
        print(xlib.get_separator())
        print('Checking the database gzip config file ...')
        (OK, error_list) = xgzip.check_gzip_config_file('database', strict=False)
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

def form_run_database_gzip_process():
    '''
    Compress/decompress database files in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Database file compression/decompression - Run process')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the compression/decompression process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The database files are going to be compressed/decompressed.')

    # compress/decompress database dataset files in the cluster
    if OK:
        devstdout = xlib.DevStdOut(xgzip.run_gzip_process.__name__)
        OK = xgzip.run_gzip_process(cluster_name, 'database', devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_database_dataset():
    '''
    Remove database in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Datasets - Remove database')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the database dataset identification
    if OK:
        database_dataset_id = cinputs.input_database_dataset_id(ssh_client, help=True)
        if database_dataset_id  == '':
            print('WARNING: There are not any databases.')
            OK = False

    # confirm the removal of the database dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The database is going to be removed.')

    # removal the database dataset in the cluster
    if OK:
        command = f'rm -fr {xlib.get_cluster_database_dir()}/{database_dataset_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            print(f'The dataset {database_dataset_id} is removed.')
        else:
            print(f'*** ERROR: Wrong command ---> {command}')
            OK = False

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_read_transfer_config_file():
    '''
    Recreate the read transfer config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset transfer - Recreate config file')

    # get the experiment identification, local directory and file pattern
    print(xlib.get_separator())
    experiment_id = cinputs.input_experiment_id(ssh_client=None, help=False)
    local_dir = cinputs.input_local_dir('read')
    file_pattern = cinputs.input_files_pattern('.*')

    # build the selected file list
    print(xlib.get_separator())
    print('Building the selected file list ...')
    selected_file_list = []
    for file in os.listdir(local_dir):
        if os.path.isfile(os.path.join(local_dir, file)) and re.match(file_pattern, file):
            selected_file_list.append(file)
    if selected_file_list != []:
        selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
        print(f'The selected file list is built: {selected_file_list_text}')
    else:
        print(f'*** ERROR: There are not files in the directory {local_dir} with the pattern {file_pattern}')
        OK = False

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The read transfer config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {xread.get_read_transfer_config_file()} is been created ...')
        (OK, error_list) = xread.create_read_transfer_config_file(experiment_id, local_dir, selected_file_list)
        if OK:
            print('The config file is created.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_read_transfer_config_file():
    '''
    Edit read transfer config file to change the parameters of each transfer.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset transfer - Edit config file')

    # get the read transfer config file
    read_transfer_config_file = xread.get_read_transfer_config_file()

    # edit the read transfer config file
    print(xlib.get_separator())
    print('The read transfer config file is being edited ...')
    command = f'{xlib.get_editor()} {read_transfer_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the config file
    if OK:
        print(xlib.get_separator())
        print('Checking the read transfer config file ...')
        (OK, error_list) = xread.check_read_transfer_config_file(strict=False)
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

def form_upload_read_dataset():
    '''
    Upload read dataset to a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset transfer - Upload dataset to a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the upload of the read dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The read dataset is going to be uploaded.')

    # upload the read dataset to the cluster
    if OK:
        devstdout = xlib.DevStdOut(xread.upload_read_dataset.__name__)
        OK = xread.upload_read_dataset(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_read_gzip_config_file():
    '''
    Recreate read file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset file compression/decompression - Recreate config file')

    # get the compress/decompress config file
    gzip_config_file = xgzip.get_gzip_config_file('read')

    # get the cluster name
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

    # get the action
    if OK:
        action = cinputs.input_code(text='Action', code_list=['compress', 'decompress'], default_code=None).lower()

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster does not have experiment data.')
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if read_dataset_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have result datasets.')
            OK = False

    # get the file pattern
    if OK:
        file_pattern = cinputs.input_files_pattern('.*' if action == 'compress' else '.*gz')

    # get the selected file list
    if OK:
        print(xlib.get_separator())
        print('Building the selected file list ...')
        read_dataset_dir = xlib.get_cluster_experiment_read_dataset_dir(experiment_id, read_dataset_id)
        selected_file_list = []
        command = f'cd {read_dataset_dir}; find . -type f -regex "./{file_pattern}"'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                selected_file_list.append(line.rstrip('\n'))
            if selected_file_list != []:
                selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
                print(f'The selected file list is built: {selected_file_list_text}')
            else:
                print(f'*** ERROR: There are not files in the directory {read_dataset_dir} with the pattern {file_pattern}')
                OK = False
        else:
            print(f'*** ERROR: Wrong command ---> {command}')

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The read compression/decompression config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {gzip_config_file} is been created ...')
        (OK, error_list) = xgzip.create_gzip_config_file(action, 'read', experiment_id, read_dataset_id, selected_file_list)
        if OK:
            print('The config file is created.')
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

def form_edit_read_gzip_config_file():
    '''
    Edit read file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset file compression/decompression - Edit config file')

    # get the read gzip config file
    read_gzip_config_file = xgzip.get_gzip_config_file('read')

    # edit the read gzip config file
    print(xlib.get_separator())
    print('The read gzip config file is being edited ...')
    command = f'{xlib.get_editor()} {read_gzip_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the read gzip config file
    if OK:
        print(xlib.get_separator())
        print('Checking the read gzip config file ...')
        (OK, error_list) = xgzip.check_gzip_config_file('read', strict=False)
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

def form_run_read_gzip_process():
    '''
    Compress/decompress read dataset files in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Read dataset file compression/decompression - Run process')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the compression/decompression process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The read files are going to be compressed/decompressed.')

    # compress/decompress read dataset files in the cluster
    if OK:
        devstdout = xlib.DevStdOut(xgzip.run_gzip_process.__name__)
        OK = xgzip.run_gzip_process(cluster_name, 'read', devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_read_dataset():
    '''
    Remove read dataset in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Datasets - Remove read dataset')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
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
            print(f'WARNING: The cluster {cluster_name} does not have experiment data.')
            OK = False

    # get the read dataset identification
    if OK:
        read_dataset_id = cinputs.input_read_dataset_id(ssh_client, experiment_id, help=True)
        if experiment_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have read datasets.')
            OK = False

    # confirm the removal of the read dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The read dataset is going to be removed.')

    # removal the read dataset in the cluster
    if OK:
        command = f'rm -fr {xlib.get_cluster_read_dir()}/{experiment_id}/{read_dataset_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            print(f'The dataset {read_dataset_id} is removed.')
        else:
            print(f'*** ERROR: Wrong command ---> {command}')
            OK = False

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_result_transfer_config_file():
    '''
    Recreate the result transfer config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Run result dataset transfer - Recreate config file')

    # get the cluster name
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

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster does not have experiment data.')
            OK = False

    # get the action and file pattern
    if OK:
        status = cinputs.input_code(text='Dataset status', code_list=['uncompressed', 'compressed'], default_code=None).lower()

    # get the run identification
    if OK:
        result_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'result', xlib.get_all_applications_selected_code(), status, help=True)
        if result_dataset_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have result datasets.')
            OK = False

    # get the file pattern and local directory
    if OK:
        if status == 'uncompressed':
            file_pattern = cinputs.input_files_pattern('.*')

    # get the file pattern and local directory
    if OK:
        local_dir = cinputs.input_local_dir('result')

    # build the selected file list
    if OK:
        if status == 'uncompressed':
            print(xlib.get_separator())
            print('Building the selected file list ...')
            selected_file_list = []
            cluster_result_dir = f'{xlib.get_cluster_result_dir()}/{experiment_id}/{result_dataset_id}'
            command = f'cd {cluster_result_dir}; find . -type f -regex "./{file_pattern}"'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                print(f'*** ERROR: Wrong command ---> {command}')
                OK = False
            if OK:
                if selected_file_list != []:
                    selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
                    print(f'The selected file list is built: {selected_file_list_text}')
                else:
                    print(f'*** WARNING: There are not files in the cluster directory {cluster_result_dir} with the pattern {file_pattern}')
                    OK = False
        elif status == 'compressed':
            selected_file_list = [None]

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The run result transfer config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {xresult.get_result_transfer_config_file()} is been created ...')
        (OK, error_list) = xresult.create_result_transfer_config_file(experiment_id, result_dataset_id, status, selected_file_list, local_dir)
        if OK:
            print('The config file is created.')
        else:
            print()
            for error in error_list:
                print(error)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_edit_result_transfer_config_file():
    '''
    Edit result transfer config file of a run to change the parameters of each transfer.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Run result dataset transfer - Edit config file')

    # get the result transfer config file
    result_transfer_config_file = xresult.get_result_transfer_config_file()

    # edit the result transfer config file
    print(xlib.get_separator())
    print('The result transfer config file is being edited ...')
    command = f'{xlib.get_editor()} {result_transfer_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the result transfer config file
    if OK:
        print(xlib.get_separator())
        print('Checking the result transfer config file ...')
        (OK, error_list) = xresult.check_result_transfer_config_file(strict=False)
        if OK:
            print('The file is OK.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_download_result_dataset():
    '''
    Download run result dataset from a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Run result dataset transfer - Download dataset from a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the download of the run result dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The run result dataset is going to be downloaded.')

    # download the run result dataset from the cluster
    if OK:
        devstdout = xlib.DevStdOut(xresult.download_result_dataset.__name__)
        OK = xresult.download_result_dataset(cluster_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_recreate_result_gzip_config_file():
    '''
    Recreate result file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Result dataset file compression/decompression - Recreate config file')

    # get the compress/decompress config file
    gzip_config_file = xgzip.get_gzip_config_file('result')

    # get the cluster name
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

    # get if whole dataset is going to be compress or decompress
    if OK:
        whole = cinputs.input_code(text='Is whole dataset going to be compressed/decompressed?', code_list=['y', 'n'], default_code='y').lower()

    # get the action and file pattern
    if OK:
        action = cinputs.input_code(text='Action', code_list=['compress', 'decompress'], default_code=None).lower()

    # get the experiment identification
    if OK:
        experiment_id = cinputs.input_experiment_id(ssh_client, help=True)
        if experiment_id == '':
            print('WARNING: The cluster does not have experiment data.')
            OK = False

    # get the result dataset identification
    if OK:
        if whole == 'y' and action == 'decompress':
            status = 'compressed'
        else:
            status = 'uncompressed'
        result_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'result', xlib.get_all_applications_selected_code(), status, help=True)
        if result_dataset_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have result datasets.')
            OK = False

    # get the file pattern
    if OK:
        if whole == 'n':
            file_pattern = cinputs.input_files_pattern('.*' if action == 'compress' else '.*gz')

    # get the selected file list
    if OK:
        if whole == 'n':
            print(xlib.get_separator())
            print('Building the selected file list ...')
            selected_file_list = []
            result_dataset_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, result_dataset_id)
            command = f'cd {result_dataset_dir}; find . -type f -regex "./{file_pattern}"'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
                if selected_file_list != []:
                    selected_file_list_text = str(selected_file_list).strip('[]').replace('\'','')
                    print(f'The selected file list is built: {selected_file_list_text}')
                else:
                    print(f'*** ERROR: There are not files in the directory {result_dataset_dir} with the pattern {file_pattern}')
                    OK = False
            else:
                print(f'*** ERROR: Wrong command ---> {command}')
        elif whole == 'y':
            selected_file_list = [None]

    # confirm the creation of the config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The result compression/decompression config file is going to be created. The previous files will be lost.')

    # create the config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {gzip_config_file} is been created ...')
        dataset_type = 'result' if whole == 'n' else 'whole-result'
        (OK, error_list) = xgzip.create_gzip_config_file(action, dataset_type, experiment_id, result_dataset_id, selected_file_list)
        if OK:
            print('The config file is created.')
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

def form_edit_result_gzip_config_file():
    '''
    Edit result file compression/decompression config file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Result dataset file compression/decompression - Edit config file')

    # get the result gzip config file
    result_gzip_config_file = xgzip.get_gzip_config_file('result')

    # edit the result gzip config file
    print(xlib.get_separator())
    print('The result gzip config file is being edited ...')
    command = f'{xlib.get_editor()} {result_gzip_config_file}'
    rc = subprocess.call(command, shell=True)
    if rc != 0:
        print(f'*** ERROR: Return code {rc} in command -> {command}')
        OK = False

    # check the result gzip config file
    if OK:
        print(xlib.get_separator())
        print('Checking the result gzip config file ...')
        (OK, error_list) = xgzip.check_gzip_config_file('result', strict=False)
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

def form_run_result_gzip_process():
    '''
    Compress/decompress result dataset files in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Result dataset file compression/decompression - Run process')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
        OK = False

    # confirm the compression/decompression process run
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The result files are going to be compressed/decompressed.')

    # compress/decompress result dataset files in the cluster
    if OK:
        devstdout = xlib.DevStdOut(xgzip.run_gzip_process.__name__)
        OK = xgzip.run_gzip_process(cluster_name, 'result', devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_result_dataset():
    '''
    Remove result dataset in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Datasets - Remove result dataset')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
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
            print(f'WARNING: The cluster {cluster_name} does not have experiment data.')
            OK = False

    # get the result dataset identification
    if OK:
        result_dataset_id = cinputs.input_result_dataset_id(ssh_client, experiment_id, 'result', xlib.get_all_applications_selected_code(), 'uncompressed', help=True)
        if result_dataset_id == '':
            print(f'WARNING: The experiment {experiment_id} does not have result datasets.')
            OK = False

    # confirm the removal of the result dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The result dataset is going to be removed.')

    # removal the result dataset in the cluster
    if OK:
        command = f'rm -fr {xlib.get_cluster_result_dir()}/{experiment_id}/{result_dataset_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            print(f'The dataset {result_dataset_id} is removed.')
        else:
            print(f'*** ERROR: Wrong command ---> {command}')
            OK = False

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_experiment():
    '''
    Remove result dataset in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Datasets - Remove result dataset')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('ERROR: There is not any running cluster.')
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
            print(f'WARNING: The cluster {cluster_name} does not have experiment data.')
            OK = False

    # confirm the removal of the result dataset
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action(f'Every read and run result dataset of the experiment {experiment_id} is going to be irreversibly removed.')

    # removal every read and result dataset in the cluster
    if OK:
        command = f'rm -fr {xlib.get_cluster_result_dir()}/{experiment_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            command = f'rm -fr {xlib.get_cluster_read_dir()}/{experiment_id}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            print(f'The datasets of experiment {experiment_id} are removed.')
        else:
            print(f'*** ERROR: Wrong command ---> {command}')
            OK = False

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the functions related to forms corresponding to dataset menu items in mode console.')
    sys.exit(0)

#-------------------------------------------------------------------------------
