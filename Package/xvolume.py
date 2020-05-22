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
This file contains the functions related to the volume operation used in both
console mode and gui mode.
'''
#-------------------------------------------------------------------------------

import os
import re
import sys
import time

import xconfiguration
import xec2
import xinstance
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_volume(volume_name, volume_type, volume_size, terminate_indicator, log, function=None):
    '''
    Create a volume in the current zone.
    '''

    # initialize the control variable
    OK = True

    # get the volume creator name
    volume_creator_name = xlib.get_volume_creator_name()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # check the volume creator is running
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking if volume creator is running ...\n')
    (master_state_code, master_state_name) = xec2.get_node_state(volume_creator_name)
    if master_state_code == 16:
        log.write('The volume creator is running.\n')
    else:
        log.write('*** WARNING: The volume creator is not running. It will be created.\n')
        OK = xinstance.create_volume_creator(log, function=None, is_menu_call=False)

    # get the master node identification
    if OK:
        node_id = xec2.get_node_id(volume_creator_name)
        if node_id == '':
            log.write('*** ERROR: The master identification of the volume creator not has been got.\n')
            OK = False

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(volume_creator_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # warn that the requirements are being verified 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking process requirements ...\n')

    # get current region and zone names
    if OK:
        region_name = xconfiguration.get_current_region_name()
        zone_name = xconfiguration.get_current_zone_name()

    # check that the zone is available
    if OK:
        if not xec2.is_zone_available(region_name, zone_name):
            log.write('*** ERROR: The zone {0} is not available in the region {1}.\n'.format(zone_name, region_name))
            OK = False

    # check that the volume is not created
    if OK:
        if xec2.is_volume_created(volume_name, zone_name):
            log.write('*** WARNING: The volume {0} is already created.\n'.format(volume_name))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # create the volume
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Creating volume {0} ...\n'.format(volume_name))
        (OK, volume_id) = xec2.create_volume(volume_name, volume_type, volume_size)
        if OK:
            log.write('The volume is created with the identification {0}.\n'.format(volume_id))
        else:
            log.write('*** ERROR: The volume is not created.\n')

    # wait for the volume status to be available
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Waiting for the volume state to be available ...\n')
        i = 0
        available_state_timeout = 120
        while True:
            time.sleep(5)
            i += 5
            volume_status = xec2.get_volume_state(volume_name, zone_name)
            if volume_status == 'available':
                log.write('The volume is now available.\n')
                break
            elif i > available_state_timeout:
                log.write('*** The volume is not available after {0} s.\n'.format(available_state_timeout))
                Ok = False
                break

    # set the aws device file and get de machine device file
    if OK:
        aws_device_file = '/dev/sdp'
        machine_device_file = xlib.get_machine_device_file(aws_device_file)

    # attach the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Attaching volume {0} to node master of volume creator ...\n'.format(volume_name))
        OK = xec2.attach_volume(node_id, volume_id, aws_device_file)
        if OK:
            log.write('The volume is attached.\n')
        else:
            log.write('*** ERROR: The volume is not attached.\n')

    # wait for the volume attachment to be available
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Waiting for the volume attachment to be available ...\n')
        i = 0
        attachment_timeout = 120
        while True:
            time.sleep(5)
            i += 5
            volume_attachments = xec2.get_volume_attachments(volume_name, zone_name)
            if volume_attachments != []:
                log.write('The volume attachment is now available.\n')
                break
            elif i > attachment_timeout:
                log.write('*** ERROR: The volume attachment is not available after {0} s.\n'.format(attachment_timeout))
                Ok = False
                break

    # wait for the device availability
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Waiting for the availabity of the device {0} ...\n'.format(machine_device_file))
        i = 0
        availability_timeout = 120
        try:
            while True:
                time.sleep(5)
                i += 5
                command = 'hdparm -z {0}'.format(machine_device_file)
                (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
                command = 'lsblk --list --noheadings --output NAME'
                (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
                for line in stdout:
                    if line == os.path.basename(machine_device_file):
                        log.write('The device is available.\n')
                        raise xlib.BreakAllLoops
                if i > availability_timeout:
                    log.write('*** ERROR: The device is not available after {0} s.\n'.format(availability_timeout))
                    OK = False
                    break
        except xlib.BreakAllLoops:
            pass

    # format the volume
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Formating volume {0} to ext4 file system type ...\n'.format(volume_name))
        command = 'mkfs -t ext4 {0} 2>&1; RC=$?; echo "RC=$RC"'.format(machine_device_file)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            log.write('The volume is formatted.\n')
            OK = True
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')
            OK = False

    # detach the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Detaching volume {0} from node master of volume creator ...\n'.format(volume_name))
        OK = xec2.detach_volume(node_id, volume_id, aws_device_file)
        if OK:
            log.write('The volume is detached.\n')
        else:
            log.write('*** ERROR: The volume is not detached.\n')

    # terminate volume creator
    if OK:
        if terminate_indicator:
            OK = xinstance.terminate_volume_creator(log, function=None, is_menu_call=False)
        else:
            log.write(f'{xlib.get_separator()}\n')
            log.write('You do not indicate to terminate the volume creator.\n')
            log.write('Remember to terminate it when you finish to create the volumes!!!\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def remove_volume(volume_name, log, function=None):
    '''
    Delete a volume in the current zone.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # warn that the requirements are being verified 
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking process requirements ...\n')

    # get current region and zone names
    region_name = xconfiguration.get_current_region_name()
    zone_name = xconfiguration.get_current_zone_name()

    # check that the zone is available
    if not xec2.is_zone_available(region_name, zone_name):
        log.write('*** ERROR: The zone {0} is not available in the region {1}.\n'.format(zone_name, region_name))
        OK = False

    # check that the volume is available
    if OK:
        volume_status = xec2.get_volume_state(volume_name, zone_name)
        if volume_status != 'available':
            log.write('*** ERROR: The volume {0} is not available in the zone {1}.\n'.format(volume_name, zone_name))
            OK = False

    # check that the volume is not linked to any cluster templates
    if OK:
        if volume_name in xconfiguration.get_volume_names_list():
            log.write('*** ERROR: The volume is linked to some cluster templates.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # delete the volume
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Removing volume {0} ...\n'.format(volume_name))
        OK = xec2.delete_volume(volume_name)
        if OK:
            log.write('The volume is been deleted. It may remain in the deleting state for several minutes.\n')
        else:
            log.write('*** ERROR: The volume is not deleted.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def mount_volume(cluster_name, node_name, volume_name, aws_device_file, mount_path, log, function=None, is_menu_call=True):
    '''
    Mount a volume in a node.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    log.write(f'{xlib.get_separator()}\n')
    log.write('Connecting SSH client ...\n')
    (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
    if OK:
        log.write('The SSH client is connected.\n')
    else:
        for error in error_list:
            log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # warn that the requirements are being verified
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking process requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # get the zone name of the node
    if OK:
        zone_name = xec2.get_node_zone_name(cluster_name, node_name)

    # get the node identification
    if OK:
        node_id = xec2.get_node_id(cluster_name, node_name)
        if node_id == '':
            log.write('*** ERROR: The {0} identification of the cluster {1} not has been got.\n'.format(node_name, cluster_name))
            OK = False

    # check the volume is created
    if OK:
        if not xec2.is_volume_created(volume_name, zone_name):
            log.write('*** ERROR: The volume {0} is not created.\n'.format(volume_name))
            OK = False

    # check the state volume
    if OK:
        if xec2.get_volume_state(volume_name, zone_name) != 'available':
            log.write('*** ERROR: The volume {0} is not available.\n'.format(volume_name))
            OK = False

    # get the volume identification
    if OK:
        volume_id = xec2.get_volume_id(volume_name, zone_name)
        if volume_id == '':
            log.write('*** ERROR: The volume identification of {0} not has been got.\n'.format(volume_name))
            OK = False

    # get instance type
    if OK:
        node_data = xec2.get_node_data_dict(cluster_name, node_name)
        instance_type = node_data['instance_type']

    # get Nitro-based instance type indicator
    if OK:
        instance_type_data_dict = xconfiguration.get_instance_type_data_dict(instance_type)
        nitro = instance_type_data_dict['nitro']

    # get the machine device file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Getting the machine device file ...\n')

        # when instance type is Nitro-based
        if nitro == 'yes':
            node_mountpoint_file = './mountpoints.txt'
            command = f'lsblk --nodeps  --output NAME,MOUNTPOINT > {node_mountpoint_file}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                local_mountpoint_file = f'{xlib.get_temp_dir()}/mountpoints.txt'
                (OK, error_list) = xssh.get_file(sftp_client, node_mountpoint_file, local_mountpoint_file)
                if OK:
                    (OK, error_list, device_order) = get_next_device_order(local_mountpoint_file)
                    if OK:
                        machine_device_file = '/dev/nvme{0}n1'.format(device_order)
                    else:
                        for error in error_list:
                            log.write(f'{error}\n')
                else:
                    for error in error_list:
                        log.write(f'{error}\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

        # when instance type is not Nitro-based
        else:
            machine_device_file = xlib.get_machine_device_file(aws_device_file)

        if OK:
            log.write('The file is got.\n')

    # create the directory in the instance
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Creating directory {0} in node {1} of cluster {2} ...\n'.format(mount_path, node_name, cluster_name))
        command = 'mkdir --parents {0}'.format(mount_path)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory is created.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # attach the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Attaching volume {0} to node {1} of cluster {2} ...\n'.format(volume_name, node_name, cluster_name))
        OK = xec2.attach_volume(node_id, volume_id, aws_device_file)
        if OK:
            log.write('The volume is attached.\n')
            time.sleep(10)
        else:
            log.write('*** ERROR: The volume is not attached.\n')

    # mount the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Mounting volume {0} in directory {1} ...\n'.format(volume_name, mount_path))
        time.sleep(10)
        command = 'mount {0} {1}'.format(machine_device_file, mount_path)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The volume is mounted.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # close the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH transport connection ...\n')
        xssh.close_ssh_transport_connection(ssh_transport)
        log.write('The connection is closed.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def unmount_volume(cluster_name, node_name, volume_name, mount_path, log, function=None, is_menu_call=True):
    '''
    Unmount a volume in a node.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    log.write(f'{xlib.get_separator()}\n')
    log.write('Connecting SSH client ...\n')
    (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
    if OK:
        log.write('The SSH client is connected.\n')
    else:
        for error in error_list:
            log.write(f'{error}\n')

    # warn that the requirements are being verified 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking process requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name,master_state_code, master_state_name))
            OK = False

    # get the zone name of the node
    if OK:
        zone_name = xec2.get_node_zone_name(cluster_name, node_name)

    # get the node identification
    if OK:
        node_id = xec2.get_node_id(cluster_name, node_name)
        if node_id == '':
            log.write('*** ERROR: The {0} identification of the cluster {1} not has been got.\n'.format(node_name, cluster_name))
            OK = False

    # check the volume is created
    if OK:
        if not xec2.is_volume_created(volume_name, zone_name):
            log.write('*** ERROR: The volume {0} is not created.\n'.format(volume_name))
            OK = False

    # check the state volume
    if OK:
        if xec2.get_volume_state(volume_name, zone_name) != 'in-use':
            log.write('*** ERROR: The volume {0} is not in-use.\n'.format(volume_name))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # get the AWS device file
    if OK:
        aws_device_file = xec2.get_volume_device_file(cluster_name, node_name, volume_name)
        if aws_device_file == '':
            log.write('*** ERROR: the device file of the volume {volume_name} is not found.\n')
            OK = False

    # get the volume identification
    if OK:
        volume_id = xec2.get_volume_id(volume_name, zone_name)
        if volume_id == '':
            log.write(f'*** ERROR: The identificaction of volume {volume_name} not has been got.\n')
            OK = False

    # unmount the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Unmounting volume {volume_name} from the mount path {mount_path} ...\n')
        command = 'umount {0}'.format(mount_path)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The volume is unmounted.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # detach the volume to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Detaching volume {volume_name} on the device {aws_device_file} ...\n')
        OK = xec2.detach_volume(node_id, volume_id, aws_device_file)
        if OK:
            while True:
                time.sleep(10)
                volume_state = xec2.get_volume_state(volume_name, zone_name)
                if volume_state == 'available':
                    log.write('The volume is detached.\n')
                    break
                else:
                    log.write(f'... volume state: {volume_state} ...\n')
        else:
            log.write('*** ERROR: The volume is not detached.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_mount_point_dict(mountpoint_file):
    '''
    Build the mount point dictionary from a file with the output of the command:
        lsblk --nodeps  --output NAME,MOUNTPOINT
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the mount point dictionary
    mount_point_dict = {}

    # open the mount point file
    try:
        mountpoint_file_id = open(mountpoint_file, mode='r', encoding='iso-8859-1')
    except Exception as e:
        error_list.append(f'*** ERROR: {e}')
        OK = False

    # load mount point data into the mount point dictionary
    if OK:

        # read the first record
        record = mountpoint_file_id.readline()

        # while there are records
        while record != '':

            # for each mounting point
            for  mounting_point in xlib.get_mounting_point_list():
                if record.find(mounting_point) > -1:
                    end = record.find(' ')
                    mount_point_dict[mounting_point] = record[:end]
                    break

            # read the next record
            record = mountpoint_file_id.readline()

    # close the mount point file
    mountpoint_file_id.close()

    # return the control variable, error list and mount point dictionary
    return (OK, error_list, mount_point_dict)

#-------------------------------------------------------------------------------

def get_machine_device_file(mountpoint_file, mount_point):
    '''
    Get the machine device file of a mount point in a node with Nitro-based instance
    type from a file with the output of the command:
        lsblk --nodeps  --output NAME,MOUNTPOINT
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialice the machine device file
    machine_device_file = None

    # open the mount point file
    try:
        mountpoint_file_id = open(mountpoint_file, mode='r', encoding='iso-8859-1')
    except Exception as e:
        error_list.append(f'*** ERROR: {e}')
        OK = False

    # get the machine device file
    if OK:

        # read the first record
        record = mountpoint_file_id.readline()

        # while there are records
        while record != '':

            if record.find(mount_point) > -1:
                end = record.find(' ')
                machine_device_file = record[:end]
                break

            # read the next record
            record = mountpoint_file_id.readline()

    # close the mount point file
    mountpoint_file_id.close()

    # return the control variable, error list and machine device file
    return (OK, error_list, machine_device_file)

#-------------------------------------------------------------------------------

def get_next_device_order(mountpoint_file):
    '''
    Get the next device order in a node with Nitro-based instance type from a file
    with the output of the command:
        lsblk --nodeps  --output NAME,MOUNTPOINT
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialice the next device order
    next_device_order = -1

    # open the mount point file
    try:
        mountpoint_file_id = open(mountpoint_file, mode='r', encoding='iso-8859-1')
    except Exception as e:
        error_list.append(f'*** ERROR: {e}')
        OK = False

    # get the next device order
    if OK:

        # read the first record
        record = mountpoint_file_id.readline()

        # initializa the last order found
        last_order = -1

        # while there are records
        while record != '':

            if record.startswith('nvme'):
                end = record.find('n1')
                try:
                    order = int(record[4:end])
                    if order > last_order:
                        last_order = order
                except Exception as e:
                    error_list.append(f'*** ERROR: {e}')
                    OK = False
                    break

            # read the next record
            record = mountpoint_file_id.readline()

    # set the next device order
    if OK:
        next_device_order = last_order + 1

    # close the mount point file
    mountpoint_file_id.close()

    # return the control variable, error list and next device order
    return (OK, error_list, next_device_order)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the functions related to the volume operation used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
