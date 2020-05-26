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
This source contains the functions related to the instance operation used in
both console mode and gui mode.
'''
#-------------------------------------------------------------------------------

import configparser
import os
import re
import subprocess
import sys
import time

import xconfiguration
import xec2
import xlib
import xssh
import xvolume

#-------------------------------------------------------------------------------

def create_instance(instance_type, log, root_volume_type, root_volumen_size, purchasing_option, max_spot_price, interruption_behavior, function=None, is_menu_call=True):
    '''
    Create a instance corresponding to a instance type.
    '''

    # initialize the control variable
    OK = True

    # set the cluster name and node name
    cluster_name = f'{xconfiguration.environment}-{instance_type}'
    node_name = 'instance'

    # get current region and zone names
    region_name = xconfiguration.get_current_region_name()
    zone_name = xconfiguration.get_current_zone_name()

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # get the instance type dictionary
    instance_type_data_dict = xconfiguration.get_instance_type_data_dict(instance_type)

    # get the dictionary of volumes
    volumes_dict = xconfiguration.get_volumes_dict()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # warn that the requirements are being verified 
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking process requirements ...\n')

    # check that the zone is available
    if OK:
        if not xec2.is_zone_available(region_name, zone_name):
            log.write(f'*** ERROR: The zone name {zone_name} is not available.\n')
            OK = False

    # check that the cluster mode is None
    if OK:
        if xec2.get_cluster_mode(cluster_name) is not None:
            log.write('*** ERROR: There is a cluster or a instance running.\n')
            OK = False

    # check that the secuerity group is not created
    if OK:
        pass

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # create the security group
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Creating the security group {xec2.build_security_group_name(instance_type)} ... \n')
        (OK, error_list, security_group_id) = xec2.create_security_group(xec2.build_security_group_name(instance_type))
        if OK:
            log.write(f'The security group {security_group_id} is created.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the instance
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Creating the instance {xconfiguration.environment}-{instance_type} ... \n')
        (OK, error_list, instance_id) = xec2.create_instance(instance_type, node_name, security_group_id, root_volume_type, root_volumen_size, purchasing_option, max_spot_price, interruption_behavior)
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('Deleting the security group of the instance ...\n')
            (OK2, error_list) = xec2.delete_security_group(security_group_id)
            if OK2:
                log.write('The security group is deleted.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
        else:
            while True:
                time.sleep(10)
                instance_state = xec2.get_node_state(cluster_name, node_name)
                instance_state_code = instance_state[0]
                instance_state_name = instance_state[1]
                # state 0 (pending)
                if instance_state_code == 0:
                    pass
                # state 16 (running)
                elif instance_state_code == 16:
                    time.sleep(60)
                    log.write(f'The instance {instance_id} is created.\n')
                    break
                # state 32 (shutting-down) or 48 (terminated) or 64 (stopping) or 80 (stopped)
                elif instance_state_code in [32, 48, 64, 80]:
                    log.write(f'*** ERROR: The instance has status {instance_state_code} ({instance_state_name}).\n')
                    OK = False
                    break
                else:
                    log.write(f'*** ERROR: The instance status {instance_state_code} ({instance_state_name}) is not controled.\n')
                    OK = False
                    break

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name, user='ubuntu')
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name, user='ubuntu')
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

    # enable remote root login safely
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Enabling remote root login safely ...\n')
        command = 'sudo sed -i "s/#PermitRootLogin prohibit-password/PermitRootLogin yes/g" /etc/ssh/sshd_config'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            # -- command = r'''sudo sed -i $'s|no-port-forwarding,no-agent-forwarding,no-X11-forwarding,command="echo \'Please login as the user \\\\"ubuntu\\\\" rather than the user \\\\"root\\\\".\';echo;sleep 10\" ||g' /root/.ssh/authorized_keys'''
            command = 'sudo cp /home/ubuntu/.ssh/authorized_keys /root/.ssh/authorized_keys'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                command = 'sudo service ssh restart'
                (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
                if OK:
                    log.write('The login is enabled.\n')
                else:
                    log.write(f'*** ERROR: Wrong command ---> {command}\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # mount volumes
    if OK:

        # get Nitro-based instance type indicator
        nitro = instance_type_data_dict['nitro']

        # initializa the device order used when the instance type is Nitro-based
        device_order = 0

        # for each volume defined in the NGScloud config file
        for volume_name in volumes_dict.keys():
            volume_id = volumes_dict[volume_name]['volume_id']
            aws_device_file = volumes_dict[volume_name]['aws_device_file']
            mount_path = volumes_dict[volume_name]['mount_path']
            if nitro == 'yes':
                device_order += 1
                machine_device_file = f'/dev/nvme{device_order}n1'
            else:
                machine_device_file = xlib.get_machine_device_file(aws_device_file)
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Volume: {volume_name} - Id: {volume_id}\n')
            # create directory where the device is going to be mounted
            log.write('Creating directory {0} ...\n'.format(mount_path))
            command = 'sudo mkdir --parents {0}'.format(mount_path)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The directory is created.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            # attach the device
            if OK:
                log.write(f'Attaching the volume in the device {aws_device_file} ({machine_device_file}) ...\n')
                OK = xec2.attach_volume(instance_id, volume_id, aws_device_file)
                if OK:
                    log.write('The volume is attached.\n')
                    time.sleep(10)
                else:
                    log.write('*** ERROR: The volume is not attached.\n')
            # mount the device
            if OK:
                log.write(f'Mounting the device in {mount_path} ...\n')
                command = f'sudo mount {machine_device_file} {mount_path}'
                (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
                if OK:
                    log.write('The volume is mounted.\n')
                else:
                    log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the infrastructure software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the infrastructure software installation script {get_infrastructure_software_installation_script()} ...\n')
        (OK, error_list) = build_infrastructure_software_installation_script(cluster_name)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the infraestructe software installation script to the node
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the infraestructe software installation script to the instance ...\n')
        cluster_path = f'./{os.path.basename(get_infrastructure_software_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_infrastructure_software_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the infraestructe software installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision ...\n')
        command = f'chmod u+x ./{os.path.basename(get_infrastructure_software_installation_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the infraestructe software installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the infraestructe software installation starter {0} ...\n'.format(get_infrastructure_software_installation_starter()))
        (OK, error_list) = build_infrastructure_software_installation_starter()
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the infraestructe software installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the infraestructe software installation starter {get_infrastructure_software_installation_starter()} ...\n')
        cluster_path = f'./{os.path.basename(get_infrastructure_software_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_infrastructure_software_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the infraestructe software installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision ...\n')
        command = f'chmod u+x ./{os.path.basename(get_infrastructure_software_installation_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the infraestructe software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the infraestructe software installation starter ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, '.', os.path.basename(get_infrastructure_software_installation_starter()), log)

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

def build_infrastructure_software_installation_script(cluster_name):
    '''
    Build the infrastructure software installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the connetion data
    (user_id, access_key_id, secret_access_key)  = xconfiguration.get_basic_aws_data()

    # get the old region and user identification
    current_region_name = xconfiguration.get_current_region_name()

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the dataset structure and NGScloud_volume
    dataset_structure = ngscloud_options_dict['dataset info']['dataset_structure'].lower()

    # write the infrastructure software installation script
    try:
        if not os.path.exists(os.path.dirname(get_infrastructure_software_installation_script())):
            os.makedirs(os.path.dirname(get_infrastructure_software_installation_script()))
        with open(get_infrastructure_software_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export DEBIAN_FRONTEND=noninteractive\n')
            script_file_id.write( 'export DEBIAN_PRIORITY=critical\n')
            # -- script_file_id.write( 'export HOST_NAME=NGScloud2_instance\n')
            script_file_id.write( 'export HOST_NAME=master\n')
            # -- script_file_id.write( 'export HOST_NAME=node\n')
            script_file_id.write( 'export HOST_IP=`curl checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function modify_hostname\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Modifying the hostname ..."\n')
            script_file_id.write( '    sudo hostnamectl set-hostname $HOST_NAME\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error hostnamectl $RC; fi\n')
            script_file_id.write( '    echo "The file is modified."\n')
            script_file_id.write( '    sudo hostnamectl\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write( '    echo "HOST NAME: $HOSTNAME"\n')
            script_file_id.write( '    echo "HOST IP: $HOST_IP"\n')
            script_file_id.write( '    echo "HOST ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            if dataset_structure in [xconfiguration.get_dataset_structure_singlevolume(), xconfiguration.get_dataset_structure_none()]:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function create_dataset_structure\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Creating the dataset structure ..."\n')
                script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_app_dir()}\n')
                script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_database_dir()}\n')
                script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_read_dir()}\n')
                script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_reference_dir()}\n')
                script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_result_dir()}\n')
                script_file_id.write( '    echo "The dataset structure is created."\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_unzip\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package unzip ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install unzip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_awscli\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the AWS CLI ..."\n')
            script_file_id.write(f'    curl "{xlib.get_awscli_url()}" -o "{xlib.get_awscli_name()}.zip"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error curl $RC; fi\n')
            script_file_id.write(f'    unzip {xlib.get_awscli_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error unzip $RC; fi\n')
            script_file_id.write( '    sudo ./aws/install\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error install $RC; fi\n')
            script_file_id.write( '    rm -rf aws\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write(f'    rm {xlib.get_awscli_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function setup_aws\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Setting up AWS ..."\n')
            script_file_id.write( '    UBUNTU_AWS_DIR=/home/ubuntu/.aws\n')
            script_file_id.write( '    mkdir --parents $UBUNTU_AWS_DIR\n')
            script_file_id.write(f'    CONFIG_FILE=$UBUNTU_AWS_DIR/config\n')
            script_file_id.write( '    echo "[default]" > $CONFIG_FILE\n')
            script_file_id.write(f'    echo "region = {current_region_name}" >> $CONFIG_FILE\n')
            script_file_id.write( '    CREDENTIALS_FILE=$UBUNTU_AWS_DIR/credentials\n')
            script_file_id.write( '    echo "[default]" > $CREDENTIALS_FILE\n')
            script_file_id.write(f'    echo "aws_access_key_id = {access_key_id}" >> $CREDENTIALS_FILE\n')
            script_file_id.write(f'    echo "aws_secret_access_key = {secret_access_key}" >> $CREDENTIALS_FILE\n')
            script_file_id.write( '    echo "AWS is set up."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function upgrade_packages\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Updating the list of available packages ..."\n')
            script_file_id.write( '    sudo --preserve-env apt-get --assume-yes --quiet update\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The list is updated."\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "Upgrading packages ..."\n')
            script_file_id.write( '    sudo --preserve-env apt-get --assume-yes --quiet --option "Dpkg::Options::=--force-confdef" --option "Dpkg::Options::=--force-confold" upgrade\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "Packages are upgraded."\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "Deleting old downloaded archive files ..."\n')
            script_file_id.write( '    sudo --preserve-env apt-get --assume-yes --quiet autoclean\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "Files are deleted."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_libtbb2\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package libtbb2 ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install libtbb2\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_libxt6\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package libxt6 ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install libxt6\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_parallel\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package parallel ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install parallel\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_sge_in_master\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Preparing the debconf database with answers to be used in the install of the master node ..."\n')
            script_file_id.write( '    echo "gridengine-master       shared/gridenginemaster string  $HOSTNAME" | sudo debconf-set-selections\n')
            # -- script_file_id.write( '    sudo debconf-set-selections <<< "gridengine-master shared/gridenginemaster string $HOSTNAME"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-master       shared/gridenginecell   string  default" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-master       shared/gridengineconfig boolean false" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridenginemaster string  $HOSTNAME" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridenginecell   string  default" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridengineconfig boolean false" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-client       shared/gridenginemaster string  $HOSTNAME" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-client       shared/gridenginecell   string  default" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-client       shared/gridengineconfig boolean false" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            # -- script_file_id.write( '    echo "postfix postfix/main_mailer_type        select  No configuration" | sudo debconf-set-selections\n')
            # -- script_file_id.write( '    RC=$?\n')
            # -- script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    sudo debconf-set-selections <<< "postfix postfix/mailname string $HOST_ADDRESS"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    sudo debconf-set-selections <<< "postfix postfix/main_mailer_type string \'Internet Site\'"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "The database is prepared."\n')
            script_file_id.write( '    echo "Installing SGE packages ..."\n')
            # -- script_file_id.write( '    wget http://ftp.debian.org/debian/pool/main/g/gridengine/gridengine-exec_8.1.9+dfsg-9+b1_amd64.deb\n')
            # -- script_file_id.write( '    wget http://ftp.debian.org/debian/pool/main/g/gridengine/gridengine-common_8.1.9+dfsg-9_all.deb\n')
            # -- script_file_id.write( '    wget http://ftp.debian.org/debian/pool/main/g/gridengine/gridengine-master_8.1.9+dfsg-9+b1_amd64.deb\n')
            # -- script_file_id.write( '    wget http://ftp.debian.org/debian/pool/main/g/gridengine/gridengine-client_8.1.9+dfsg-9+b1_amd64.deb\n')
            script_file_id.write( '    sudo apt --assume-yes install gridengine-common gridengine-exec gridengine-master gridengine-client\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt $RC; fi\n')
            script_file_id.write( '    echo "Package are installed."\n')
            # -- script_file_id.write( '    adduser gsadmin --uid 500\n') ??????
            script_file_id.write( '    echo "Starting SGE ..."\n')
            script_file_id.write( '    sudo -u sgeadmin /usr/share/gridengine/scripts/init_cluster /var/lib/gridengine default /var/spool/gridengine/spooldb sgeadmin\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error init_cluster $RC; fi\n')
            script_file_id.write( '    sudo service gridengine-master restart\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error service $RC; fi\n')
            script_file_id.write( '    echo "SGE is started."\n')
            script_file_id.write( '    echo "Disabling Postfix ..."\n')
            script_file_id.write( '    sudo service postfix stop\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error service $RC; fi\n')
            script_file_id.write( '    sudo update-rc.d postfix disable\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error update-rc.d $RC; fi\n')
            script_file_id.write( '    echo "Postfix is disabled."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_sge_in_worker\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    export MASTER_HOSTNAME=master\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Preparing the debconf database with answers to be used in the install of the worker node ..."\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridenginemaster string  $MASTER_HOSTNAME" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridenginecell   string  default" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-common       shared/gridengineconfig boolean false" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-client       shared/gridenginemaster string  $MASTER_HOSTNAME" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "gridengine-client       shared/gridenginecell   string  default" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    eecho "gridengine-client       shared/gridengineconfig boolean false" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "postfix postfix/main_mailer_type        select  No configuration" | sudo debconf-set-selections\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "The database is prepared."\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    echo "Installing SGE packages ..."\n')
            script_file_id.write( '    sudo apt --assume-yes install gridengine-exec gridengine-client\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt $RC; fi\n')
            script_file_id.write( '    echo "Package are installed."\n')
            script_file_id.write( '    echo "Disabling Postfix ..."\n')
            script_file_id.write( '    sudo service postfix stop\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error service $RC; fi\n')
            script_file_id.write( '    sudo update-rc.d postfix disable\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error update-rc.d $RC; fi\n')
            script_file_id.write( '    echo "Postfix is disabled."\n')
            script_file_id.write( '    echo "Starting SGE ..."\n')
            script_file_id.write( '    echo $MASTER_HOSTNAME | sudo tee /var/lib/gridengine/default/common/act_qmaster\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tee $RC; fi\n')
            script_file_id.write( '    sudo service gridengine-exec restart\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error service $RC; fi\n')
            script_file_id.write( '    echo "SGE is started."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_xorg\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package xorg ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install xorg\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_texlive\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package texlive ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-latex-base\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-fonts-recommended\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-fonts-extra\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-latex-extra\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function end\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail ok\n')
            script_file_id.write( '    exit 0\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function manage_error\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "ERROR: $1 returned error $2"\n')
            script_file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail wrong\n')
            script_file_id.write( '    exit 3\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            process_name = 'Infrastructure software installation'
            mail_message_ok = xlib.get_mail_message_ok(process_name, cluster_name)
            mail_message_wrong = xlib.get_mail_message_wrong(process_name, cluster_name)
            script_file_id.write( 'function send_mail\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {process_name}"\n')
            script_file_id.write( '    if [ "$1" == "ok" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_ok}"\n')
            script_file_id.write( '    elif [ "$1" == "wrong" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_wrong}"\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '         MESSAGE=""\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '    DESTINATION_FILE=mail-destination.json\n')
            script_file_id.write( '    echo "{" > $DESTINATION_FILE\n')
            script_file_id.write(f'    echo "    \\\"ToAddresses\\\":  [\\\"{xconfiguration.get_contact_data()}\\\"]," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"CcAddresses\\\":  []," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"BccAddresses\\\":  []" >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "}" >> $DESTINATION_FILE\n')
            script_file_id.write( '    MESSAGE_FILE=mail-message.json\n')
            script_file_id.write( '    echo "{" > $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Subject\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Data\\\":  \\\"$SUBJECT\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Body\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Html\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Data\\\":  \\\"$MESSAGE\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "}" >> $MESSAGE_FILE\n')
            script_file_id.write(f'    aws ses send-email --from {xconfiguration.get_contact_data()} --destination file://$DESTINATION_FILE --message file://$MESSAGE_FILE\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_duration\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
            script_file_id.write( '    HH=`expr $DURATION / 3600`\n')
            script_file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
            script_file_id.write( '    SS=`expr $DURATION % 60`\n')
            script_file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'modify_hostname\n')
            script_file_id.write( 'init\n')
            if dataset_structure in [xconfiguration.get_dataset_structure_singlevolume(), xconfiguration.get_dataset_structure_none()]:
                script_file_id.write( 'create_dataset_structure\n')
            script_file_id.write( 'install_unzip\n')
            script_file_id.write( 'install_awscli\n')
            script_file_id.write( 'setup_aws\n')
            script_file_id.write( 'upgrade_packages\n')
            script_file_id.write( 'install_libtbb2\n')
            script_file_id.write( 'install_libxt6\n')
            script_file_id.write( 'install_parallel\n')
            script_file_id.write( '# -- install_sge_in_master\n')
            script_file_id.write( '# -- install_sge_in_worker\n')
            script_file_id.write( 'install_xorg\n')
            script_file_id.write( 'install_texlive\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_infrastructure_software_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_infrastructure_software_installation_starter():
    '''
    Build the starter of the infrastructure software installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the infrastructure software installation starter
    try:
        if not os.path.exists(os.path.dirname(get_infrastructure_software_installation_starter())):
            os.makedirs(os.path.dirname(get_infrastructure_software_installation_starter()))
        with open(get_infrastructure_software_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format(f'./{os.path.basename(get_infrastructure_software_installation_script())} &>./{get_infrastructure_software_installation_log()}'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_infrastructure_software_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_infrastructure_software_installation_script():
    '''
    Get the infrastructure software installation script path in the local computer.
    '''

    # assign infrastructure software installation script path
    infrastructure_software_installation_script = f'{xlib.get_temp_dir()}/infrastructure-software-installation.sh'

    # return the infrastructure software installation script path
    return infrastructure_software_installation_script

#-------------------------------------------------------------------------------

def get_infrastructure_software_installation_starter():
    '''
    Get the infrastructure software installation starter path in the local computer.
    '''

    # assign the infrastructure software installation starter path
    infrastructure_software_installation_starter = f'{xlib.get_temp_dir()}/infrastructure-software-installation-starter.sh'

    # return the infrastructure software installation starter path
    return infrastructure_software_installation_starter

#-------------------------------------------------------------------------------

def get_infrastructure_software_installation_log():
    '''
    Get the infrastructure software installation log path in the instance.
    '''

    # assign the infrastructure software installation log path
    infrastructure_software_installation_log = './infrastructure-software-installation.log'

    # return the infrastructure software installation log path
    return infrastructure_software_installation_log

#-------------------------------------------------------------------------------

def terminate_instance(cluster_name, log, function=None, is_menu_call=True):
    '''
    Terminate a instance.
    '''

    # initialize the control variable
    OK = True

    # get the instance type
    instance_type = xec2.get_instance_type(cluster_name)

    # get the instance identification
    instance_id = xec2.get_node_id(cluster_name)

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # check the instance is running
    (instance_state_code, instance_state_name) = xec2.get_node_state(cluster_name)
    if instance_state_code != 16:
        log.write(f'*** WARNING: The cluster {cluster_name} is not running. Its state is {instance_state_code} ({instance_state_name}).\n')

    # terminate the instance
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Terminating the instance {instance_id} ...\n')
    if instance_id is not None:
        (OK, error_list) = xec2.terminate_instance(instance_id)
        if OK:
            while True:
                time.sleep(10)
                instance_state = xec2.get_node_state(cluster_name)
                instance_state_code = instance_state[0]
                instance_state_name = instance_state[1]
                if instance_state_code in [48, -1]:
                    time.sleep(60)
                    log.write('The instance is terminated.\n')
                    break
                else:
                    log.write(f'... status instance: {instance_state_code} ({instance_state_name}) ....\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
    else:
        log.write('WARNING: The instance is not found.\n')

    # delete the security group
    log.write(f'{xlib.get_separator()}\n')
    log.write('Deleting the security group of the instance ...\n')
    (OK, error_list, security_group_id) = xec2.get_security_group_id(xec2.build_security_group_name(instance_type))
    if security_group_id is not None:
        (OK, error_list) = xec2.delete_security_group(security_group_id)
        if OK:
            log.write('The security group is deleted.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
    else:
        log.write('WARNING: The security group is not found.\n')

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

def create_volume_creator(log, function=None, is_menu_call=True):
    '''
    Create a instance of the volume creator.
    '''

    # initialize the control variable
    OK = True

    # set the instance type of the volume creator, cluster name and node name
    instance_type='t2.micro'
    cluster_name = f'{xconfiguration.environment}-volume-creator'
    node_name = 'instance'

    # get the volume creator name
    volume_creator_name = xlib.get_volume_creator_name()

    # get current region and zone names
    region_name = xconfiguration.get_current_region_name()
    zone_name = xconfiguration.get_current_zone_name()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # warn that the requirements are being verified 
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking process requirements ...\n')

    # check that the zone is available
    if OK:
        if not xec2.is_zone_available(region_name, zone_name):
            log.write(f'*** ERROR: The zone name {zone_name} is not available.\n')
            OK = False

    # check that the cluster mode is None
    if OK:
        if xec2.get_cluster_mode(cluster_name) is not None:
            log.write('*** ERROR: There is other volume creator running.\n')
            OK = False

    # check that the secuerity group is not created
    if OK:
        pass

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # create the security group
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Creating the security group {xec2.build_security_group_name("volume-creator")} ... \n')
        (OK, error_list, security_group_id) = xec2.create_security_group(xec2.build_security_group_name('volume-creator'))
        if OK:
            log.write(f'The security group {security_group_id} is created.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the instance
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Creating the volume creator ... \n')
        (OK, error_list, instance_id) = xec2.create_instance(instance_type, node_name, security_group_id, root_volume_type='standard', root_volumen_size=8, purchasing_option=xec2.get_purchasing_option_ondemand(), max_spot_price=0.0, interruption_behavior=xec2.get_interruption_behavior_terminate())
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('Deleting the security group of the instance ...\n')
            (OK2, error_list) = xec2.delete_security_group(security_group_id)
            if OK2:
                log.write('The security group is deleted.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
        else:
            while True:
                time.sleep(10)
                instance_state = xec2.get_node_state(cluster_name, node_name)
                instance_state_code = instance_state[0]
                instance_state_name = instance_state[1]
                # state 0 (pending)
                if instance_state_code == 0:
                    pass
                # state 16 (running)
                elif instance_state_code == 16:
                    time.sleep(60)
                    log.write(f'The volume creator is created.\n')
                    break
                # state 32 (shutting-down) or 48 (terminated) or 64 (stopping) or 80 (stopped)
                elif instance_state_code in [32, 48, 64, 80]:
                    log.write(f'*** ERROR: The volume creator has status {instance_state_code} ({instance_state_name}).\n')
                    OK = False
                    break
                else:
                    log.write(f'*** ERROR: The volume creator status {instance_state_code} ({instance_state_name}) is not controled.\n')
                    OK = False
                    break

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name, user='ubuntu')
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # enable remote root login safely
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Enabling remote root login safely ...\n')
        command = 'sudo sed -i "s/#PermitRootLogin prohibit-password/PermitRootLogin yes/g" /etc/ssh/sshd_config'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            command = 'sudo cp /home/ubuntu/.ssh/authorized_keys /root/.ssh/authorized_keys'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                command = 'sudo service ssh restart'
                (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
                if OK:
                    log.write('The login is enabled.\n')
                else:
                    log.write(f'*** ERROR: Wrong command ---> {command}\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

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

def terminate_volume_creator(log, function=None, is_menu_call=True):
    '''
    Terminate a instance.
    '''

    # initialize the control variable
    OK = True

    # set the instance type of the volume creator, cluster name and node name
    instance_type='t3a.micro'
    cluster_name = f'{xconfiguration.environment}-volume-creator'
    node_name = 'instance'

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # get the instance identification
    instance_id = xec2.get_node_id(cluster_name)

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut) and is_menu_call:
        log.write('This process might take a few minutes. Do not close this window, please wait!\n')

    # check the instance is running
    (instance_state_code, instance_state_name) = xec2.get_node_state(cluster_name)
    if instance_state_code != 16:
        log.write(f'*** WARNING: The volume_creator is not running. Its state is {instance_state_code} ({instance_state_name}).\n')

    # terminate the instance
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Terminating the volume creator {instance_id} ...\n')
    if instance_id is not None:
        (OK, error_list) = xec2.terminate_instance(instance_id)
        if OK:
            while True:
                time.sleep(10)
                instance_state = xec2.get_node_state(cluster_name)
                instance_state_code = instance_state[0]
                instance_state_name = instance_state[1]
                if instance_state_code in [48, -1]:
                    time.sleep(60)
                    log.write('The volume creator is terminated.\n')
                    break
                else:
                    log.write(f'... status volume creator: {instance_state_code} ({instance_state_name}) ....\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
    else:
        log.write('WARNING: The instance is not found.\n')

    # delete the security group
    log.write(f'{xlib.get_separator()}\n')
    log.write('Deleting the security group of the instance ...\n')
    (OK, error_list, security_group_id) = xec2.get_security_group_id(xec2.build_security_group_name('volume-creator'))
    if security_group_id is not None:
        (OK, error_list) = xec2.delete_security_group(security_group_id)
        if OK:
            log.write('The security group is deleted.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

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

if __name__ == '__main__':
     print('This file contains the functions related to the instane operation used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
