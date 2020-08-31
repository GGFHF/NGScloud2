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
This file contains the general functions to data inputs in mode console.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys

import xconfiguration
import xcluster
import xddradseqtools
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def input_code(text, code_list, default_code):
    '''
    Input a code selected from a code list.
    '''

    # initialize the code
    code = ''

    # get the code list text
    code_list_text = str(code_list).strip('[]').replace('\'','').replace(',', ' or')

    # input and check the code
    while code == '':
        if default_code is None:
            code = input(f'{text} ({code_list_text}): ').lower()
        else:
            code = input(f'{text} ({code_list_text}) [{default_code}]: ').lower()
            if code == '': code = default_code
        found = False
        for i in range(len(code_list)):
            if code.lower() == code_list[i].lower():
                code = code_list[i]
                found = True
                break
        if not found:
            print(f'*** ERROR: {code} is not in {code_list_text}.')
            code = ''

    # return the code
    return code

#-------------------------------------------------------------------------------

def input_int(text, default=None, minimum=(-sys.maxsize - 1), maximum=sys.maxsize):
    '''
    Input a integer number.
    '''

    # initialize the number
    literal = None

    # input and check the integer number
    while literal is None:
        if default is None:
            literal = input(f'{text}: ')
        else:
            literal = input(f'{text} [{default}]: ')
            if literal == '': literal = default
        if not xlib.check_int(literal, minimum, maximum):
            print(f'*** ERROR: {literal} is not a valid value.')
            literal = None

    # return the integer value
    return int(literal)

#-------------------------------------------------------------------------------

def input_float(text, default=None, minimum=float(-sys.maxsize - 1), maximum=float(sys.maxsize), mne=0.0, mxe=0.0):
    '''
    Input a float number.
    '''

    # initialize the number
    literal = None

    # input and check the float number
    while literal is None:
        if default is None:
            literal = input(f'{text}: ')
        else:
            literal = input(f'{text} [{default}]: ')
            if literal == '': literal = default
        if not xlib.check_float(literal, minimum, maximum, mne, mxe):
            print(f'*** ERROR: {literal} is not a valid value.')
            literal = None

    # return the float value
    return float(literal)

#-------------------------------------------------------------------------------

def input_access_key_id(default_access_key_id):
    '''
    Input an access key identification.
    '''

    # initialize the access key identification
    access_key_id = ''

    # input and check the access key identification
    while access_key_id == '':
        if default_access_key_id != '':
            access_key_id = input(f'Enter the AWS access key id [{default_access_key_id}]: ')
            if access_key_id == '':
                access_key_id = default_access_key_id
        else:
            access_key_id = input('Enter the AWS access key id: ')

    # return the access key identification
    return access_key_id

#-------------------------------------------------------------------------------

def input_secret_access_key(default_secret_access_key):
    '''
    Input a secret access key.
    '''

    # initialize the secret access key
    secret_access_key = ''

    # input and check the secret access key
    while secret_access_key == '':
        if default_secret_access_key != '':
            secret_access_key = input(f'Enter the AWS secret access key [{default_secret_access_key[:5] + "*****" + default_secret_access_key[len(default_secret_access_key)-5:]}]: ')
            if secret_access_key == '':
                secret_access_key = default_secret_access_key
        else:
            secret_access_key = input('Enter the AWS secret access key: ')

    # return the secret access key
    return secret_access_key

#-------------------------------------------------------------------------------

def input_user_id(default_user_id):
    '''
    Input an user identification.
    '''

    # initialize the user identification
    user_id = ''

    # input and check the user identification
    while user_id == '':
        if default_user_id != '':
            user_id = input(f'Enter the AWS user id [{default_user_id}]: ')
            if user_id == '':
                user_id = default_user_id
        else:
            user_id = input('Enter the AWS user id: ')

    # return the user identification
    return user_id

#-------------------------------------------------------------------------------

def input_email(default_email):
    '''
    Input a contact e-mail address.
    '''

    # initialize the contact e-mail address
    email = ''

    # input and check the contact e-mail address
    while email == '':
        if default_email != '':
            email = input(f'Enter the contact e-mail [{default_email}]: ')
            if email == '':
                email = default_email
            elif not xlib.is_email_address_valid(email):
                print(f'*** ERROR: {email} is not a valid e-mail address.')
                email = ''
        else:
            email = input('Enter the contact e-mail: ')

    # return the contact e-mail address
    return email

#-------------------------------------------------------------------------------

def input_region_name(default_region_name, help):
    '''
    Input a region name.
    '''

    # initialize the control variable
    OK = True

    # initialize the region name
    region_name = ''

    # get the available region name list
    if help:
        available_region_names_list = xec2.get_available_region_list()
        if available_region_names_list == []:
            OK = False

    # check that the default region name is available
    if OK and help:
        if default_region_name not in available_region_names_list:
            default_region_name = ''

    # print the available region names
    if OK and help:
        available_region_names_list_test = str(available_region_names_list).strip('[]').replace('\'','')
        print(f'Available region names: {available_region_names_list_test} ...')

    # input and check the region name
    if OK:
        while region_name == '':
            if help:
                if default_region_name != '':
                    region_name = input(f'... Enter the region name [{default_region_name}]: ')
                else:
                    region_name = input('... Enter the region name: ')
                if region_name == '':
                    region_name = default_region_name
                elif region_name not in available_region_names_list:
                    print(f'*** ERROR: {region_name} is not available.')
                    region_name = ''
            else:
                if default_region_name != '':
                    region_name = input(f'Enter the region name [{default_region_name}]: ')
                else:
                    region_name = input('Enter the region name: ')
                if region_name == '':
                    region_name = default_region_name

    # return the region name
    return region_name

#-------------------------------------------------------------------------------

def input_zone_name(region_name, default_zone_name, help):
    '''
    Input a zone name.
    '''

    # initialize the control variable
    OK = True

    # initialize the zone name
    zone_name = ''

    # get the available zone name list of the region
    if help:
        available_zone_names_list = xec2.get_available_zone_list(region_name)
        if available_zone_names_list == []:
            OK = False

    # check that the default zone name is available
    if OK and help:
        if default_zone_name not in available_zone_names_list:
            default_zone_name = ''

    # print the available zone names
    if OK and help:
        available_zone_names_list_test = str(available_zone_names_list).strip('[]').replace('\'','')
        print(f'Available region names: {available_zone_names_list_test} ...')

    # input and check the zone name
    if OK:
        while zone_name == '':
            if help:
                if default_zone_name != '':
                    zone_name = input(f'... Enter the zone name [{default_zone_name}]: ')
                else:
                    zone_name = input('... Enter the zone name: ')
                if zone_name == '':
                    zone_name = default_zone_name
                elif zone_name not in available_zone_names_list:
                    print(f'*** ERROR: {zone_name} is not available.')
                    zone_name = ''
            else:
                if default_zone_name != '':
                    zone_name = input(f'Enter the zone name [{default_zone_name}]: ')
                else:
                    zone_name = input('Enter the zone name: ')
                if zone_name == '':
                    zone_name = default_zone_name

    # return the zone name
    return zone_name

#-------------------------------------------------------------------------------

def input_instance_type(cluster_mode, help):
    '''
    Input an instance type.
    '''

    # initialize the instance type
    instance_type = ''

    # initialize the instance type list
    instance_type_list = []

    # get the instance type dictionary
    if help:
        instance_type_dict = xconfiguration.get_instance_type_dict(cluster_mode)
        for key in sorted(instance_type_dict.keys()):
            instance_type_list.append(instance_type_dict[key]['id'])

    # print instance types
    if help:
        # set data width
        id_width = 12
        vcpu_width = 4
        memory_width = 12
        use_width = 17
        # set line
        line = '{0:' + str(id_width) + '} {1:' + str(vcpu_width) + '} {2:' + str(memory_width) + '} {3:' + str(use_width) + '}'
        # print header
        print('Available instance types:')
        print()
        print(line.format('Instance Type', 'vCPU', 'Memory (GiB)', 'Use'))
        print(line.format('=' * id_width, '=' * vcpu_width, '=' * memory_width, '=' * use_width))
        # print detail lines
        for key in sorted(instance_type_dict.keys()):
            print(line.format(instance_type_dict[key]['id'], instance_type_dict[key]['vcpu'], instance_type_dict[key]['memory'], instance_type_dict[key]['use']))
        print()

    # input and check the instance type
    while instance_type == '':
        if help:
            instance_type = input('... Enter the instance type: ')
            if instance_type not in instance_type_list:
                print(f'*** ERROR: {instance_type} is not available.')
                instance_type = ''
        else:
            instance_type = input('Enter the instance type: ')

    # return the instance type
    return instance_type

#-------------------------------------------------------------------------------

def input_cluster_name(volume_creator_included, help):
    '''
    Input a cluster name.
    '''

    # initialize the control variable
    OK = True

    # initialize the default cluster name
    default_cluster_name = ''

    # initialize the cluster name
    cluster_name = ''

    # get the running cluster list
    if help:
        clusters_running_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=volume_creator_included)
        if clusters_running_list == []:
            OK = False
        elif len(clusters_running_list) == 1:
            default_cluster_name = clusters_running_list[0]

    # print the running clusters
    if OK and help:
        clusters_running_list_text = str(clusters_running_list).strip('[]').replace('\'','')
        print(f'Running clusters: {clusters_running_list_text} ...')

    # input and check the cluster name
    if OK:
        while cluster_name == '':
            if help:
                if default_cluster_name != '':
                    cluster_name = input(f'... Enter the cluster name [{default_cluster_name}]: ')
                else:
                    cluster_name = input('... Enter the cluster name: ')
                if cluster_name == '':
                    cluster_name = default_cluster_name
                elif cluster_name not in clusters_running_list:
                    print(f'*** ERROR: {cluster_name} is not running.')
                    cluster_name = ''
            else:
                cluster_name = input('Enter the cluster name: ')

    # return the cluster name
    return cluster_name

#-------------------------------------------------------------------------------

def input_node_name(cluster_name, new, is_master_valid, help):
    '''
    Input a node name.
    '''

    # initialize the control variable
    OK = True

    # initialize the node name
    node_name = ''

    # get the running node list
    running_node_name_list = xec2.get_cluster_node_list(cluster_name)

    # get the valid node list
    valid_node_name_list = running_node_name_list
    if not is_master_valid:
        if valid_node_name_list.remove('master') == []:
            OK = False

    # print the running nodes
    if OK:
        if help:
            valid_node_name_list_text = str(valid_node_name_list).strip('[]').replace('\'','')
            print(f'Running nodes: {valid_node_name_list_text} ...')

    # input and check the node name
    if OK:
        while node_name == '':
            if help:
                node_name = input('... Enter the node name: ')
            else:
                node_name = input('Enter the node name: ')

            if not is_master_valid and node_name in ['master', 'instance']:
                print(f'*** ERROR: {node_name} is not a node name valid.')
                node_name = ''
            elif not node_name.isalnum() or not node_name[0].isalpha():
                print(f'*** ERROR: {node_name} is not an alphanumeric string or the first character is not alphabetic')
                node_name = ''
            elif new and node_name in running_node_name_list:
                print(f'*** ERROR: {node_name} is already running..')
                node_name = ''
            elif not new and node_name not in running_node_name_list:
                print(f'*** ERROR: {node_name} is not running.')
                node_name = ''

    # return the node name
    return node_name

#-------------------------------------------------------------------------------

def input_volume_name(text, zone_name, type, allowed_none, help):
    '''
    Input a volume name.
    '''

    # initialize the control variable
    OK = True

    # initialize the volume name
    volume_name = ''

    # get the available volume list
    if help:
        if type == 'created':
            volume_name_list = xec2.get_created_volume_name_list(zone_name)
        elif type == 'linked':
            volume_name_list = xconfiguration.get_volume_names_list()
        if volume_name_list == [] or volume_name_list == ['']:
            OK = False

    # print the available volumes
    if OK and help:
        volume_name_list_text = str(volume_name_list).strip('[]').replace('\'','')
        print(f'Available volume names: {volume_name_list_text} ...')

    # input and check the volume name
    if OK:
        if allowed_none:
            volume_name_list.append('NONE')
        while volume_name == '':
            if help:
                if allowed_none:
                    volume_name = input(f'... Enter the {text} or NONE: ')
                else:
                    volume_name = input(f'... Enter the {text}: ')
                found = False
                for i in range(len(volume_name_list)):
                    if volume_name.lower() == volume_name_list[i].lower():
                        volume_name = volume_name_list[i]
                        found = True
                        break
                if not found:
                    print(f'*** ERROR: {volume_name} is not available.')
                    volume_name = ''
            else:
                volume_name = input('Enter the volume name: ')

    # return the volume name
    return volume_name

#-------------------------------------------------------------------------------

def input_device_file(node_name, volume_name):
    '''
    Input a node device file where the volume is going to be attached.
    '''

    # initialize the device file pattern, the default device file and the device file
    device_file_pattern = '/dev/sd[m-p]'
    default_device_file = '/dev/sdm'
    device_file = ''

    # input and check the device file
    while device_file == '':
        device_file = input(f'Enter the device file of {node_name} where {volume_name} is going to be attached [{default_device_file}]: ')
        if device_file == '':
            device_file = default_device_file
        elif not xlib.is_device_file(device_file, device_file_pattern):
            print(f'***ERROR: The device file has to have a pattern {device_file_pattern}')
            device_file = ''

   # return the defice file
    return device_file

#-------------------------------------------------------------------------------

def input_local_dir(files_type):
    '''
    Input a local directory path.
    '''

    # initialize the local directory path
    local_dir_path = ''

    # set the input message text
    if files_type == 'reference':
        message_text = 'Enter the local directory path where the reference files are: '
    if files_type == 'database':
        message_text = 'Enter the local directory path where the built database files are: '
    elif files_type == 'read':
        message_text = 'Enter the local directory path where the read files are: '
    elif files_type == 'result':
        message_text = 'Enter the local directory path where the result files will be download: '

    # input and check the local directory path
    while local_dir_path == '':
        local_dir_path = input(message_text)
        if not os.path.isdir(local_dir_path):
            print(f'***ERROR: The directory {local_dir_path} does not exist.')
            local_dir_path = ''

   # return the local directory path
    return local_dir_path

#-------------------------------------------------------------------------------

def input_local_file_path(file_type):
    '''
    Input a local file path.
    '''

    # initialize the local file path
    local_file_path = ''

    # set the input message text
    if file_type == 'RADdesigner-conditions':
        message_text = 'Enter the local file path with RADdesigner conditions: '
    elif file_type == 'RADdesigner-samples':
        message_text = 'Enter the local file path with RADdesigner samples: '

    # input and check the local file path
    while local_file_path == '':
        local_file_path = input(message_text)
        if not os.path.isfile(local_file_path):
            print(f'***ERROR: The file {local_file_path} does not exist.')
            local_file_path = ''

   # return the local file path
    return local_file_path

#-------------------------------------------------------------------------------

def input_cluster_directory(files_type):
    '''
    Input a local directory path.
    '''

    # initialize the local directory path
    directoy_path = ''

    # set the input message text
    if files_type == 'reference':
        message_text = f'Enter the {xlib.get_cluster_reference_dir()} subdirectory where the references will be upload: '
    elif files_type == 'database':
        message_text = f'Enter the {xlib.get_cluster_database_dir()} subdirectory where the references will be upload: '

    # input and check the local directory path
    while directoy_path == '':
        directoy_path = input(message_text)
        if not xlib.is_valid_path(directoy_path):
            print(f'***ERROR: The directory {directoy_path} is not valid.')
            directoy_path = ''

   # return the directory path
    return directoy_path

#-------------------------------------------------------------------------------

def input_mount_path(node_name):
    '''
    Input a node directory path where the device file is going to be mounted.
    '''

    # initialize the mount path
    mount_path = ''

    # input and check the mount path
    while mount_path == '':
        mount_path = input(f'Enter the mount path in {node_name}: ')
        if not xlib.is_absolute_path(mount_path, 'linux'):
            print(f'***ERROR: {mount_path} is not an absolute path.')
            mount_path = ''

   # return the mount path
    return mount_path

#-------------------------------------------------------------------------------

def input_mounting_point():
    '''
    Input a mounting point of a volume in the cluster.
    '''

    # get mounting points list
    mounting_points_list = xlib.get_mounting_point_list()

    # print the available volumes
    if help:
        mounting_points_list_text = str(mounting_points_list).strip('[]').replace('\'','')
        print(f'Available mounting points: {mounting_points_list_text} ...')

    # initialize the mounting poing
    mounting_point = ''

    # input and check the mounting poing
    while mounting_point == '':
        mounting_point = input('... Enter the mounting point: ')
        if mounting_point not in mounting_points_list:
            print(f'*** ERROR: {mounting_point} is not valid.')
            mounting_point = ''

   # return the mounting point
    return mounting_point

#-------------------------------------------------------------------------------

def input_reference_dataset_id(ssh_client, allowed_none, help):
    '''
    Input a reference dataset identification in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the reference dataset identification
    reference_dataset_id = ''

    # initialize the reference dataset identification list
    reference_dataset_id_list = []

    # get the reference dataset identifications
    if help:
        command = f'ls {xlib.get_cluster_reference_dir()}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    reference_dataset_id_list.append(line)

    # print the reference identifications in the clusters
    if OK and help:
        if reference_dataset_id_list != []:
            reference_dataset_id_list_text = str(reference_dataset_id_list).strip('[]').replace('\'','')
            print(f'Reference dataset ids existing in the cluster: {reference_dataset_id_list_text} ...')
        else:
            OK = False

    # input and check the reference dataset identification
    if OK:
        if allowed_none:
            reference_dataset_id_list.append('NONE')
        while reference_dataset_id == '':
            if help:
                if allowed_none:
                    reference_dataset_id = input('... Enter the reference dataset id or NONE: ')
                else:
                    reference_dataset_id = input('... Enter the reference dataset id: ')
                if reference_dataset_id not in reference_dataset_id_list:
                    print(f'*** ERROR: {reference_dataset_id} does not exist.')
                    reference_dataset_id = ''
            else:
                reference_dataset_id = input('Enter the reference dataset id: ')
                if not reference_dataset_id.isidentifier():
                    print('*** ERROR: The reference id has some non-alphanumeric characters.')
                    reference_dataset_id = ''

    # return the reference dataset identification
    return reference_dataset_id

#-------------------------------------------------------------------------------

def input_reference_file(ssh_client, reference_dataset_id, help):
    '''
    Input a reference file of a reference dataset in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the reference file
    reference_file = ''

    # initialize the reference file list
    file_list = []

    # get the reference files of the reference dataset identification
    if help:
        command = f'find {xlib.get_cluster_reference_dataset_dir(reference_dataset_id)} -maxdepth 1 -type f'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found' and line.lower().find('gtf') == -1:
                    file_list.append(os.path.basename(line))

    # print the files in the clusters
    if OK and help:
        if file_list != []:
            file_list_text = str(file_list).strip('[]').replace('\'','')
            print(f'Files existing in the reference dataset: {file_list_text} ...')
        else:
            OK = False

    # input and check the reference file
    if OK:
        while reference_file == '':
            if help:
                reference_file = input('... Enter the reference file: ')
                if reference_file not in file_list:
                    print(f'*** ERROR: {reference_file} does not exist.')
                    reference_file = ''
            else:
                reference_file = input('Enter the reference file: ')
                if not reference_file.isidentifier():
                    print('*** ERROR: The reference file has some non-alphanumeric characters.')
                    reference_file = ''

    # return the reference file
    return reference_file

#-------------------------------------------------------------------------------

def input_transcriptome_file(ssh_client, reference_dataset_id, help):
    '''
    Input a transcriptome file of a reference dataset in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the transcriptome file
    transcriptome_file = ''

    # initialize the file list
    file_list = []

    # get the files of the reference dataset identification
    if help:
        command = f'find {xlib.get_cluster_reference_dataset_dir(reference_dataset_id)} -maxdepth 1 -type f'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found' and line.lower().find('gtf') == -1:
                    file_list.append(os.path.basename(line))

    # print the files in the clusters
    if OK and help:
        if file_list != []:
            file_list_text = str(file_list).strip('[]').replace('\'','')
            print(f'Files existing in the reference dataset: {file_list_text} ...')
        else:
            OK = False

    # input and check the transcriptome file
    if OK:
        while transcriptome_file == '':
            if help:
                transcriptome_file = input('... Enter the transcriptome file: ')
                if transcriptome_file not in file_list:
                    print(f'*** ERROR: {transcriptome_file} does not exist.')
                    transcriptome_file = ''
            else:
                transcriptome_file = input('Enter the transcriptome file: ')
                if not transcriptome_file.isidentifier():
                    print('*** ERROR: The transcriptome file has some non-alphanumeric characters.')
                    transcriptome_file = ''

    # return the transcriptome file
    return transcriptome_file

#-------------------------------------------------------------------------------

def input_gtf_file(ssh_client, reference_dataset_id, help):
    '''
    Input a GTF file of a reference dataset in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the GTF file
    gtf_file = ''

    # initialize the GTF file list
    gtf_file_list = []

    # get the GTF files of the reference dataset identification
    if help:
        command = f'find {xlib.get_cluster_reference_dataset_dir(reference_dataset_id)} -maxdepth 1 -type f'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found' and line.lower().find('gtf') > -1:
                    gtf_file_list.append(os.path.basename(line))

    # print the GTF files in the clusters
    if OK and help:
        if gtf_file_list != []:
            gtf_file_list_text = str(gtf_file_list).strip('[]').replace('\'','')
            print(f'GTF files existing in the reference dataset: {gtf_file_list_text} ...')
        else:
            OK = False

    # input and check the GTF file
    if OK:
        while gtf_file == '':
            if help:
                gtf_file = input('... Enter the GTF file: ')
                if gtf_file not in gtf_file_list:
                    print(f'*** ERROR: {gtf_file} does not exist.')
                    gtf_file = ''
            else:
                gtf_file = input('Enter the GTF file: ')
                if not gtf_file.isidentifier():
                    print('*** ERROR: The GTF file has some non-alphanumeric characters.')
                    gtf_file = ''

    # return the GTF dataset identification
    return gtf_file

#-------------------------------------------------------------------------------

def input_reference_file2(ssh_client, reference_dataset_id, type, allowed_none, help):
    '''
    Input a file (reference, transcriptome, mask, splice sites, exons) of a reference dataset in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the splice site file
    additional_reference_file = ''

    # initialize the file list
    file_list = []

    # get the splice site files of the reference dataset identification
    if help:
        command = f'find {xlib.get_cluster_reference_dataset_dir(reference_dataset_id)} -maxdepth 1 -type f'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    file_list.append(os.path.basename(line))

    # print the GTF files in the clusters
    if OK and help:
        if file_list != []:
            file_list_text = str(file_list).strip('[]').replace('\'','')
            print(f'Files existing in the reference dataset: {file_list_text} ...')
        else:
            OK = False

    # input and check the additional reference file
    if OK:
        if allowed_none:
            file_list.append('NONE')
        while additional_reference_file == '':
            if help:
                if allowed_none:
                    additional_reference_file = input(f'... Enter the {type} file or NONE: ')
                else:
                    additional_reference_file = input(f'... Enter the {type} file: ')
                if additional_reference_file not in file_list:
                    print(f'*** ERROR: {additional_reference_file} does not exist.')
                    additional_reference_file = ''
            else:
                if allowed_none:
                    additional_reference_file = input(f'Enter the {type} file or NONE: ')
                else:
                    additional_reference_file = input(f'Enter the {type} file: ')
                if not additional_reference_file.isidentifier():
                    print(f'*** ERROR: The {type} file has some non-alphanumeric characters.')
                    additional_reference_file = ''

    # return the GTF dataset identification
    return additional_reference_file

#-------------------------------------------------------------------------------

def input_database_dataset_id(ssh_client, help):
    '''
    Input a database dataset identification in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the database dataset identification
    database_dataset_id = ''

    # initialize the database dataset identification list
    database_dataset_id_list = []

    # get the database dataset identifications
    if help:
        command = f'ls {xlib.get_cluster_database_dir()}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    database_dataset_id_list.append(line)

    # print the database identifications in the clusters
    if OK and help:
        if database_dataset_id_list != []:
            database_dataset_id_list_text = str(database_dataset_id_list).strip('[]').replace('\'','')
            print(f'Database ids existing in the cluster: {database_dataset_id_list_text} ...')
        else:
            OK = False

    # input and check the database dataset identification
    if OK:
        while database_dataset_id == '':
            if help:
                database_dataset_id = input('... Enter the database id: ')
                if database_dataset_id not in database_dataset_id_list:
                    print(f'*** ERROR: {database_dataset_id} does not exist.')
                    database_dataset_id = ''
            else:
                database_dataset_id = input('Enter the database id: ')
                if not database_dataset_id.isidentifier():
                    print('*** ERROR: The database id has some non-alphanumeric characters.')
                    database_dataset_id = ''

    # return the database dataset identification
    return database_dataset_id

#-------------------------------------------------------------------------------

def input_protein_database_name(ssh_client, database_dataset_id, help):
    '''
    Input a protein database name of a database dataset in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the database file names
    protein_database_name = ''

    # initialize the database file list
    database_file_name_list = []

    # get the list of the database file names
    if help:
        command = f'ls {xlib.get_cluster_database_dataset_dir(database_dataset_id)}/*phr'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    database_file_name_list.append(os.path.basename(line))
            if database_file_name_list == []:
                OK = False

    # get the list of the database dataset names
    if OK and help:
        protein_database_name_list = []
        pattern = re.compile('^.*[0-9][0-9]$')
        for database_file_name in database_file_name_list:
            file_name, file_extension = os.path.splitext(database_file_name)
            if pattern.match(file_name):
                database_name = file_name[:-3]
            else:
                database_name = file_name
            if database_name not in protein_database_name_list:
                protein_database_name_list.append(database_name)

    # print the database name in the dataset
    if OK and help:
        protein_database_name_list_text = str(protein_database_name_list).strip('[]').replace('\'','')
        print(f'Protein databases existing in the dataset: {protein_database_name_list_text} ...')

    # input and check the protein database name
    if OK:
        while protein_database_name == '':
            if help:
                protein_database_name = input('... Enter the protein database: ')
                if protein_database_name not in protein_database_name_list:
                    print(f'*** ERROR: {protein_database_name} does not exist.')
                    protein_database_name = ''
            else:
                protein_database_name = input('Enter the protein database: ')

    # return the protein database name
    return protein_database_name

#-------------------------------------------------------------------------------

def input_experiment_id(ssh_client, help):
    '''
    Input an experiment/process identification of the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the experiment/process identification
    experiment_id = ''

    # initialize the experiment/process identification list
    experiment_id_list = []

    # get the experiment/process identifications
    if help:
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

    # print the experiment/process identifications in the clusters
    if OK and help:
        if experiment_id_list != []:
            experiment_id_list_text = str(experiment_id_list).strip('[]').replace('\'','')
            print(f'Experiment/process ids existing in the cluster: {experiment_id_list_text} ...')
        else:
            OK = False

    # input and check the experiment/process identification
    if OK:
        while experiment_id == '':
            if help:
                experiment_id = input('... Enter the experiment/process id: ')
                if experiment_id not in experiment_id_list:
                    print(f'*** ERROR: {experiment_id} does not exist.')
                    experiment_id = ''
            else:
                experiment_id = input('Enter the experiment/process id: ')
                if not experiment_id.isidentifier():
                    print('*** ERROR: The experiment/process id has some non-alphanumeric characters.')
                    experiment_id = ''

    # return the experiment identification
    return experiment_id

#-------------------------------------------------------------------------------

def input_read_dataset_id(ssh_client, experiment_id, help):
    '''
    Input a read dataset identification of an experimient in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the read dataset identification
    read_dataset_id = ''

    # initialize the read dataset identification list
    read_dataset_id_list = []

    # get the read dataset identifications of the experiment
    if help:
        command = f'ls {xlib.get_cluster_read_dir()}/{experiment_id}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                read_dataset_id_list.append(line.rstrip('\n'))

    # print the read dataset identifications in the experiment
    if OK and help:
        if read_dataset_id_list != []:
            read_dataset_id_list_text = str(read_dataset_id_list).strip('[]').replace('\'','')
            print(f'Read dataset ids existing in the experiment: {read_dataset_id_list_text} ...')
        else:
            OK = False

    # input and check the read dataset identification
    if OK:
        while read_dataset_id == '':
            if help:
                read_dataset_id = input('... Enter the read dataset id: ')
                if read_dataset_id not in read_dataset_id_list:
                    print(f'*** ERROR: {read_dataset_id} does not exist.')
                    read_dataset_id = ''
            else:
                read_dataset_id = input('Enter the read dataset id: ')

    # return the read dataset identification
    return read_dataset_id

#-------------------------------------------------------------------------------

def input_read_type():
    '''
    Input the read type of the reads files.
    '''

    # initialize the read type
    read_type = ''

    # input and check the read type
    while read_type == '':
        read_type = input('Read type (SE/PE)?: ').upper()
        if read_type not in ['SE', 'PE']:
            print(f'*** ERROR: {read_type} is not SE or PE.')
            read_type = ''

    # return the read type
    return read_type

#-------------------------------------------------------------------------------

def input_result_dataset_id(ssh_client, experiment_id, type, app_list, status, help):
    '''
    Input a result dataset identification of an experimient in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the result dataset identification
    result_dataset_id = ''

    # initialize the result dataset list
    result_dataset_id_list = []

    # get the result dataset identifications of the experiment
    if help:
        if status == 'uncompressed':
            command = f'cd {xlib.get_cluster_result_dir()}/{experiment_id}; for list in `ls`; do ls -ld $list | grep -v ^- > /dev/null && echo $list; done;'
        elif status == 'compressed':
            command = f'cd {xlib.get_cluster_result_dir()}/{experiment_id}; for list in `ls`; do ls -ld $list | grep -v ^d > /dev/null && echo $list; done;'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                for app in app_list:
                    if app == xlib.get_all_applications_selected_code() or line.startswith(app):
                        result_dataset_id_list.append(line.rstrip('\n'))
                      
                        break
            if result_dataset_id_list != []:
                result_dataset_id_list.sort()

    # print the result dataset identifications in the clusters
    if OK and help:
        if result_dataset_id_list != []:
            result_dataset_id_list_text = str(result_dataset_id_list).strip('[]').replace('\'','')
            print(f'{type} dataset ids existing in the experiment: {result_dataset_id_list_text} ...')
        else:
            OK = False

    # input and check the result dataset identification
    if OK:
        while result_dataset_id == '':
            if help:
                result_dataset_id = input(f'... Enter the {type} dataset id: ')
                if result_dataset_id not in result_dataset_id_list:
                    print(f'*** ERROR: {result_dataset_id} does not exist.')
                    result_dataset_id = ''
            else:
                result_dataset_id = input(f'Enter the {type} dataset id: ')

    # return the result dataset identification
    return result_dataset_id

#-------------------------------------------------------------------------------

def input_result_dataset_id_list(ssh_client, experiment_id, type, app_list, status, help):
    '''
    Input a list of result dataset identifications of an experimient in the cluster connected by a ssh client.
    '''

    # initialize the control variable
    OK = True

    # initialize the list with all result dataset identifications
    all_result_dataset_id_list = []

    # initialize the list with selected result dataset identifications
    selected_result_dataset_id_list = []

    # get the result dataset identifications of the experiment
    if help:
        if status == 'uncompressed':
            command = f'cd {xlib.get_cluster_result_dir()}/{experiment_id}; for list in `ls`; do ls -ld $list | grep -v ^- > /dev/null && echo $list; done;'
        elif status == 'compressed':
            command = f'cd {xlib.get_cluster_result_dir()}/{experiment_id}; for list in `ls`; do ls -ld $list | grep -v ^d > /dev/null && echo $list; done;'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                for app in app_list:
                    if app == xlib.get_all_applications_selected_code() or line.startswith(app):
                        all_result_dataset_id_list.append(line.rstrip('\n'))
                        break
            if result_dataset_id_list != []:
                result_dataset_id_list.sort()

    # print the result dataset identifications in the clusters
    if OK and help:
        if all_result_dataset_id_list != []:
            all_result_dataset_id_list_text = str(all_result_dataset_id_list).strip('[]').replace('\'','')
            print(f'{type} dataset ids existing in the experiment: {all_result_dataset_id_list_text} ...')
        else:
            OK = False

    # input result dataset identifications and check them
    if OK:
        result_dataset_id = ''
        while result_dataset_id.upper() != 'END':
            if help:
                result_dataset_id = input(f'... Enter the {type} dataset id or END to finish: ')
                if result_dataset_id.upper() != 'END' and result_dataset_id not in all_result_dataset_id_list:
                    print(f'*** ERROR: {result_dataset_id} does not exist.')
                    result_dataset_id = ''
            else:
                result_dataset_id = input(f'Enter the {type} dataset id or END to finish: ')
            if result_dataset_id != '' and result_dataset_id.upper() != 'END':
                if result_dataset_id not in selected_result_dataset_id_list:
                    selected_result_dataset_id_list.append(result_dataset_id)
                else:
                    print(f'*** ERROR: {result_dataset_id} is already entered.')

    # sort the list with selected result dataset identifications
    if OK:
        if selected_result_dataset_id_list != []:
            selected_result_dataset_id_list.sort()

    # return the list with selected result dataset identifications
    return selected_result_dataset_id_list

#-------------------------------------------------------------------------------

def input_assembly_type(help):
    '''
    Input an assembly type
    '''

    # get assembly type list
    assembly_type_list = ['CONTIGS', 'SCAFFOLDS']

    # print the assembly types
    if help:
        assembly_type_list_text = str(assembly_type_list).strip('[]').replace('\'','')
        print(f'Available assembly types: {assembly_type_list_text} ...')

    # initialize the assembly type
    assembly_type = ''

    # input and check the mounting poing
    while assembly_type == '':
        if help:
            assembly_type = input('... Enter the assembly type: ')
            if assembly_type not in assembly_type_list:
                print(f'*** ERROR: {assembly_type} is not valid.')
                assembly_type = ''
        else:
            assembly_type = input('Enter the assembly type: ')

   # return the assembly type
    return assembly_type

#-------------------------------------------------------------------------------

def input_enzyme(enzyme_number, restriction_enzyme_dict, allowed_ambiguity_codes, help):
    '''
    Input the identification of a restriction enzyme or its restriction site sequence
    '''

    # set the enzyme ordinal
    if enzyme_number == '1':
        enzyme_ordinal = '1st'
    elif enzyme_number == '2':
        enzyme_ordinal = '2nd'

    # get the enzime identification list
    enzyme_id_list = list(restriction_enzyme_dict.keys())
    enzyme_id_list.sort()

    # print the enzime identifications
    if help:
        enzyme_id_list_text = str(enzyme_id_list).strip('[]').replace('\'','')
        print(f'Available restriction enzyme ids: {enzyme_id_list_text} ...')

    # initialize the enzyme value
    enzyme = ''

    # input and check the mounting poing
    while enzyme == '':
        if help:
            enzyme = input(f'... Enter the id of {enzyme_ordinal} restriction enzyme or a restriction site sequence ("*" is the cut point): ')
            if enzyme not in enzyme_id_list and not xlib.is_valid_sequence(seq=enzyme, allowed_ambiguity_codes=allowed_ambiguity_codes, other_allowed_characters_list=[], cut_tag_check=True):
                print(f'*** ERROR: {enzyme} is not valid.')
                enzyme = ''
        else:
            enzyme = input('Enter the enzyme: ')

   # return the enzyme
    return enzyme

#-------------------------------------------------------------------------------

def input_files_pattern(default_file_pattern):
    '''
    Input a file pattern.
    '''

    # initialize the file pattern
    file_pattern = ''

    # input and check the file pattern
    while file_pattern == '':
        if default_file_pattern != '':
            file_pattern = input(f'Enter a file pattern [{default_file_pattern}]: ')
            if file_pattern == '':
                file_pattern = default_file_pattern
        else:
            file_pattern = input('Enter a file pattern: ')
        try:
            re.compile(file_pattern)
        except:
            print('Invalid pattern. It has to be a valid regular expression.')
            file_pattern = ''

   # return the file pattern
    return file_pattern

#-------------------------------------------------------------------------------

def input_file_pairing_specific_chars(file_order, default_specific_chars):
    '''
    Input specific characteres used to identify file 1 or file 2 when the read type is paired.
    '''

    # initialize the specific characters
    specific_chars = ''

    # input and check the file specific characters
    while specific_chars == '':
        if default_specific_chars != '':
            specific_chars = input(f'Enter file #{file_order} specific chars [{default_specific_chars}]: ')
            if specific_chars == '':
                specific_chars = default_specific_chars
        else:
            specific_chars = input(f'Enter file #{file_order} specific chars: ')

   # return the specific characters
    return specific_chars

#-------------------------------------------------------------------------------

def input_assembly_origin(help):
    '''
    Input an assembly origin
    '''

    # get assembly origin list
    assembly_origin_list = ['NGSCLOUD', 'EXTERNAL']

    # print the assembly origins
    if help:
        assembly_origin_list_text = str(assembly_origin_list).strip('[]').replace('\'','')
        print(f'Available assembly origins: {assembly_origin_list_text} ...')

    # initialize the assembly origin
    assembly_origin = ''

    # input and check the mounting poing
    while assembly_origin == '':
        if help:
            assembly_origin = input('... Enter the assembly origin: ')
            if assembly_origin not in assembly_origin_list:
                print(f'*** ERROR: {assembly_origin} is not valid.')
                assembly_origin = ''
        else:
            assembly_origin = input('Enter the assembly origin: ')

   # return the assembly origin
    return assembly_origin

#-------------------------------------------------------------------------------

def input_database_list(candidate_database_list, last_database):
    '''
    Input an ordered list of databases to annotate.
    '''

    # set the candidate database number 
    candidate_database_num = len(candidate_database_list)

    # initialize the list of selected databases 
    selected_database_list = []

    # initialize the order number
    order_number = 1

    # input database identifications and check them
    database = ''
    while database.upper() != 'END' and order_number <= candidate_database_num:

        # input a candidate database
        if order_number == 1:
            candidate_database_list_text = str(candidate_database_list).strip('[]').replace('\'','')
            print(f'All candidate databases: {candidate_database_list_text} ...')
            print(f'(if {last_database} is selected, it has to be the last)')
            database = input(f'Enter the database {order_number}: ')
        else:
            candidate_database_list_text = str(candidate_database_list).strip('[]').replace('\'','')
            print(f'Remaining candidate databases: {candidate_database_list_text} ...')
            print(f'(if {last_database} is selected, it has to be the last)')
            database = input(f'Enter the database {order_number} or END to finish: ')
        if database != '' and database.upper() != 'END':
            if database in candidate_database_list:
                selected_database_list.append(database)
                candidate_database_list.remove(database)
                order_number += 1
                if database == last_database:
                    database = 'END'
            else:
                print(f'*** ERROR: {database} is not in candidate database list.')
        elif database.upper() == 'END' and order_number == 1:
            database = ''
            print('*** ERROR: You have to input at least one database.')

    # return the selected database list
    return selected_database_list

#-------------------------------------------------------------------------------

def input_batch_job_id(ssh_client, help):
    '''
    Input a batch job identification.
    '''

    # initialize the control variable
    OK = True

    # initialize the list of identification of the batch jobs
    batch_job_id_list = []

    # initialize the batch job identification
    batch_job_id = ''

    # get the batch job dictionary and the batch job identification list
    if help:
        (OK, error_list, batch_job_dict) = xcluster.get_batch_job_dict(ssh_client)

    # build the list of identifications of the batch jobs
    if OK and help:
        for job_id in batch_job_dict.keys():
            batch_job_id_list.append(job_id)
        if batch_job_id_list != []:
            batch_job_id_list.sort()
        else:
            OK = False

    # print the batch jobs
    if OK and help:
        # set data width
        job_id_width = 6
        job_name_width = 10
        state_width = 15
        start_date_width = 10
        start_time_width = 10
        # set line
        line = '{0:' + str(job_id_width) + '} {1:' + str(job_name_width) + '} {2:' + str(state_width) + '} {3:' + str(start_date_width) + '} {4:' + str(start_time_width) + '}'
        # print header
        print('Batch jobs:')
        print()
        print(line.format('Job id', 'Job name', 'State', 'Start date', 'Start time'))
        print(line.format('=' * job_id_width, '=' * job_name_width, '=' * state_width, '=' * start_date_width, '=' * start_time_width))
        # print detail lines
        for job_id in batch_job_id_list:
            job_name = batch_job_dict[job_id]['job_name']
            state = batch_job_dict[job_id]['state']
            start_date = batch_job_dict[job_id]['start_date']
            start_time = batch_job_dict[job_id]['start_time']
            print(line.format(job_id, job_name, state, start_date, start_time))
        print()

    # input and check the batch job identification
    if OK:
        while batch_job_id == '':
            if help:
                batch_job_id = input('... Enter the batch job id: ')
                if batch_job_id not in batch_job_id_list:
                    print(f'*** ERROR: {batch_job_id} is not an id of a batch job.')
                    batch_job_id = ''
            else:
                batch_job_id = input('Enter the batch job id: ')

    # return the batch job identification
    return batch_job_id

#-------------------------------------------------------------------------------

def input_bioinfoapp_version(name):
    '''
    Input a version of a BioInfo application.
    '''

    # warm about the use of a previous version
    print (f'{xlib.get_project_name()} is designed to use the last version of {name} and its parameters. The use of some parameter in some previous versions could cause a malfunction.')

    # input and check the version
    version = input('Version (last or especific identification) [last]?: ')
    if version == '':
        version = 'last'

    # return the version
    return version

#-------------------------------------------------------------------------------

def input_submission_log_file():
    '''
    Input a submission process log.
    '''

    # initialize the control variable
    OK = True

    # initialize the submission log file
    submission_log_file = ''

    # get the submission log file list
    submission_log_file_list = [file for file in os.listdir(xlib.get_log_dir()) if os.path.isfile(os.path.join(xlib.get_log_dir(), file)) and file.startswith(xconfiguration.environment)]

    # print the submission log file list in the local computer
    if submission_log_file_list != []:
        submission_log_file_list_text = str(submission_log_file_list).strip('[]').replace('\'','')
        print(f'Submission log files: {submission_log_file_list_text} ...')
    else:
        OK = False

    # input and check the submission log file
    if OK:
        while submission_log_file == '':
            submission_log_file = input('... Enter the submission log file: ')
            if submission_log_file not in submission_log_file_list:
                print(f'*** ERROR: {submission_log_file} does not exist.')
                submission_log_file = ''

    # return the submission log file
    return submission_log_file

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the general functions to data inputs in mode console.')
    sys.exit(0)

#-------------------------------------------------------------------------------
