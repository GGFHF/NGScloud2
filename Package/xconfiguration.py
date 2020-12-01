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
This file contains the functions related to the NGScloud2 configuration used in
both console mode and gui mode.
'''

#-------------------------------------------------------------------------------

import configparser
import os
import sys

import xec2
import xlib
import xstarcluster

#-------------------------------------------------------------------------------

# Global variables

environment = ''    # the AWS environment where the processes are run

#-------------------------------------------------------------------------------

def get_environments_file():
    '''
    Get the defined environments file.
    '''

    # assign the defined environments file
    environments_file = f'{xlib.get_config_dir()}/environments.txt'

    # return the defined environments file
    return environments_file

#-------------------------------------------------------------------------------

def get_environments_list():
    '''
    Get the defined environments list
    '''

    # initialize the environments list
    environments_list = []

    # read the environments file
    try:
        environments_file = get_environments_file()
        with open(environments_file, mode='r', encoding='iso-8859-1', newline='\n') as file_id:
            records = file_id.readlines()
            for record in records:
                environments_list.append(record.strip())
    except Exception:
        pass

    # sort the environments list
    if environments_list != []:
        environments_list.sort()

    # return the environments list
    return environments_list

#-------------------------------------------------------------------------------

def add_environment(environment):
    '''
    Add a new environment in the environments file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # add environment in the environments file
    try:
        environments_file = get_environments_file()
        if not os.path.exists(os.path.dirname(environments_file)):
            os.makedirs(os.path.dirname(environments_file))
        with open(environments_file, mode='a', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write(f'{environment.strip()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {environments_file} can not be written.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_cluster_name(instance_type):
    '''
    Build the cluster name from the environment and the instace type.
    '''

    return f'{environment}-{instance_type}'

#-------------------------------------------------------------------------------

def get_cluster_mode_native():
    '''
    Get the mode value when the cluster is started with the managing of the native
    instance use developed in this project.
    '''

    return 'native'

#-------------------------------------------------------------------------------

def get_cluster_mode_starcluster():
    '''
    Get the mode value when the cluster is started using Starcluster.
    '''

    return 'StarCluster'

#-------------------------------------------------------------------------------

def get_cluster_mode_list():
    '''
    Get the cluster mode list.
    '''

    return [get_cluster_mode_native(), get_cluster_mode_starcluster()]

#-------------------------------------------------------------------------------

def get_dataset_structure_singlevolume():
    '''
    Get the value when the dataset structure is on single-volume.
    '''

    return 'single-volume'

#-------------------------------------------------------------------------------

def get_dataset_structure_multivolume():
    '''
    Get the value when the dataset structure is on multi-volume.
    '''

    return 'multi-volume'

#-------------------------------------------------------------------------------

def get_dataset_structure_none():
    '''
    Get the value when there is not any volume linked.
    '''

    return 'none'

#-------------------------------------------------------------------------------
    
def get_dataset_structure_list():
    '''
    Get the code list of "dataset_structure".
    '''

    return [get_dataset_structure_singlevolume(), get_dataset_structure_multivolume(), get_dataset_structure_none()]

#-------------------------------------------------------------------------------

def get_instance_type_dict(cluster_mode):
    '''
    Get the instance type dictionary used by the NGScloud.
    '''

    # initialize the instance type dictionary
    instance_type_dict = {}

    # add data of instance type defined in StarCluster by default to the instance type dictionary
    if cluster_mode in get_cluster_mode_list():

        # T2 instance types
        instance_type_dict['1_1_01'] = {'id': 't2.micro', 'vcpu': 1, 'memory': 1, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.3 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'current'}
        instance_type_dict['1_1_02'] = {'id': 't2.small', 'vcpu': 1, 'memory': 2, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.3 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'current'}
        instance_type_dict['1_1_03'] = {'id': 't2.medium', 'vcpu': 2, 'memory': 4, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.3 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'current'}

        # M3 instance types
        instance_type_dict['1_2_01'] = {'id': 'm3.medium', 'vcpu': 1, 'memory': 3.75, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['1_2_02'] = {'id': 'm3.large', 'vcpu': 2, 'memory': 7.5, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['1_2_03'] = {'id': 'm3.xlarge', 'vcpu': 4, 'memory': 15, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['1_2_04'] = {'id': 'm3.2xlarge', 'vcpu': 8, 'memory': 30, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}

        # C3 instance types
        instance_type_dict['2_1_01'] = {'id': 'c3.large', 'vcpu': 2, 'memory': 3.75, 'processor': 'Intel Xeon E5-2680 v2 (Ivy Bridge)', 'speed': '2.8 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['2_1_02'] = {'id': 'c3.xlarge', 'vcpu': 4, 'memory': 7.5, 'processor': 'Intel Xeon E5-2680 v2 (Ivy Bridge)', 'speed': '2.8 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['2_1_03'] = {'id': 'c3.2xlarge', 'vcpu': 8, 'memory': 15, 'processor': 'Intel Xeon E5-2680 v2 (Ivy Bridge)', 'speed': '2.8 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['2_1_04'] = {'id': 'c3.4xlarge', 'vcpu': 16, 'memory': 30, 'processor': 'Intel Xeon E5-2680 v2 (Ivy Bridge)', 'speed': '2.8 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['2_1_05'] = {'id': 'c3.8xlarge', 'vcpu': 32, 'memory': 30, 'processor': 'Intel Xeon E5-2680 v2 (Ivy Bridge)', 'speed': '2.8 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}

        # R3 instance types
        instance_type_dict['3_1_01'] = {'id': 'r3.large', 'vcpu': 2, 'memory': 15, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['3_1_02'] = {'id': 'r3.xlarge', 'vcpu': 4, 'memory': 30.5, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['3_1_03'] = {'id': 'r3.2xlarge', 'vcpu': 8, 'memory': 61, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['3_1_04'] = {'id': 'r3.4xlarge', 'vcpu': 16, 'memory': 122, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}
        instance_type_dict['3_1_05'] = {'id': 'r3.8xlarge', 'vcpu': 32, 'memory': 244, 'processor': 'Intel Xeon E5-2670 v2 (Ivy Bridge)', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'native', 'generation': 'previous'}

    # add data of additional instance types updated in StarCluster with the workaround
    if cluster_mode in get_cluster_mode_list():

        # t2 instance types
        instance_type_dict['1_1_04'] = {'id': 't2.large', 'vcpu': 2, 'memory': 8, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.0 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_1_05'] = {'id': 't2.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.0 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_1_06'] = {'id': 't2.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'Intel Xeon Scalable', 'speed': 'up to 3.0 GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}

        # m4 instance types
        instance_type_dict['1_2_11'] = {'id': 'm4.large', 'vcpu': 2, 'memory': 8, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_2_12'] = {'id': 'm4.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_2_13'] = {'id': 'm4.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_2_14'] = {'id': 'm4.4xlarge', 'vcpu': 16, 'memory': 64, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_2_15'] = {'id': 'm4.10xlarge', 'vcpu': 40, 'memory': 160, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['1_2_16'] = {'id': 'm4.16xlarge', 'vcpu': 64, 'memory': 256, 'processor': 'Intel Xeon E5-2676 v3 (Haswell) or Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.4  or 2.3  GHz', 'use': 'general purpose', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}

        # c4 instance types
        instance_type_dict['2_1_11'] = {'id': 'c4.large', 'vcpu': 2, 'memory': 3.75, 'processor': 'Intel Xeon E5-2666 v3 (Haswell)', 'speed': '2.9 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['2_1_12'] = {'id': 'c4.xlarge', 'vcpu': 4, 'memory': 7.5, 'processor': 'Intel Xeon E5-2666 v3 (Haswell)', 'speed': '2.9 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['2_1_13'] = {'id': 'c4.2xlarge', 'vcpu': 8, 'memory': 15, 'processor': 'Intel Xeon E5-2666 v3 (Haswell)', 'speed': '2.9 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['2_1_14'] = {'id': 'c4.4xlarge', 'vcpu': 16, 'memory': 30, 'processor': 'Intel Xeon E5-2666 v3 (Haswell)', 'speed': '2.9 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['2_1_15'] = {'id': 'c4.8xlarge', 'vcpu': 36, 'memory': 60, 'processor': 'Intel Xeon E5-2666 v3 (Haswell)', 'speed': '2.9 GHz', 'use': 'compute optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}

        # r4 instance types
        instance_type_dict['3_1_11'] = {'id': 'r4.large', 'vcpu': 2, 'memory': 15.250, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['3_1_12'] = {'id': 'r4.xlarge', 'vcpu': 4, 'memory': 30.5, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['3_1_13'] = {'id': 'r4.2xlarge', 'vcpu': 8, 'memory': 61, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['3_1_14'] = {'id': 'r4.4xlarge', 'vcpu': 16, 'memory': 122, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['3_1_15'] = {'id': 'r4.8xlarge', 'vcpu': 32, 'memory': 244, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}
        instance_type_dict['3_1_16'] = {'id': 'r4.16xlarge', 'vcpu': 64, 'memory': 488, 'processor': 'Intel Xeon E5-2686 v4 (Broadwell)', 'speed': '2.3 GHz', 'use': 'memory optimized', 'nitro': 'no', 'starcluster': 'workaround', 'generation': 'current'}

    # add data of other instance types only used when the cluster mode is "native"
    if cluster_mode == get_cluster_mode_native():

        # t3 instance types
        instance_type_dict['1_1_11'] = {'id': 't3.micro', 'vcpu': 2, 'memory': 1, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_12'] = {'id': 't3.small', 'vcpu': 2, 'memory': 2, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_13'] = {'id': 't3.medium', 'vcpu': 2, 'memory': 4, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_14'] = {'id': 't3.large', 'vcpu': 2, 'memory': 8, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_15'] = {'id': 't3.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_16'] = {'id': 't3.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'Intel Xeon Scalable', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # t3a instance types
        instance_type_dict['1_1_21'] = {'id': 't3a.micro', 'vcpu': 2, 'memory': 1, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_22'] = {'id': 't3a.small', 'vcpu': 2, 'memory': 2,  'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_23'] = {'id': 't3a.medium', 'vcpu': 2, 'memory': 4,  'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_24'] = {'id': 't3a.large', 'vcpu': 2, 'memory': 8,  'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_25'] = {'id': 't3a.xlarge', 'vcpu': 4, 'memory': 16,  'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_1_26'] = {'id': 't3a.2xlarge', 'vcpu': 8, 'memory': 32,  'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # t4g instance types
        # -- instance_type_dict['1_1_31'] = {'id': 't4g.micro', 'vcpu': 2, 'memory': 1, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_1_32'] = {'id': 't4g.small', 'vcpu': 2, 'memory': 2,  'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_1_33'] = {'id': 't4g.medium', 'vcpu': 2, 'memory': 4,  'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_1_34'] = {'id': 't4g.large', 'vcpu': 2, 'memory': 8,  'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_1_35'] = {'id': 't4g.xlarge', 'vcpu': 4, 'memory': 16,  'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_1_36'] = {'id': 't4g.2xlarge', 'vcpu': 8, 'memory': 32,  'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # m5 instance types
        instance_type_dict['1_2_21'] = {'id': 'm5.large', 'vcpu': 2, 'memory': 8, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_22'] = {'id': 'm5.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_23'] = {'id': 'm5.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_24'] = {'id': 'm5.4xlarge', 'vcpu': 16, 'memory': 64, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_25'] = {'id': 'm5.8xlarge', 'vcpu': 32, 'memory': 128, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_26'] = {'id': 'm5.12xlarge', 'vcpu': 48, 'memory': 192, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_27'] = {'id': 'm5.16xlarge', 'vcpu': 64, 'memory': 256, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_28'] = {'id': 'm5.24xlarge', 'vcpu': 96, 'memory': 384, 'processor': 'Intel Xeon Platinum 8175 with new Intel Advanced Vector Extension (AXV-512) instruction set', 'speed': 'up to 3.1 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # m5a instance types
        instance_type_dict['1_2_31'] = {'id': 'm5a.large', 'vcpu': 2, 'memory': 8, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_32'] = {'id': 'm5a.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_33'] = {'id': 'm5a.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_34'] = {'id': 'm5a.4xlarge', 'vcpu': 16, 'memory': 64, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_35'] = {'id': 'm5a.8xlarge', 'vcpu': 32, 'memory': 128, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_36'] = {'id': 'm5a.12xlarge', 'vcpu': 48, 'memory': 192, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_37'] = {'id': 'm5a.16xlarge', 'vcpu': 64, 'memory': 256, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['1_2_38'] = {'id': 'm5a.24xlarge', 'vcpu': 96, 'memory': 384, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # m6g instance types
        # -- instance_type_dict['1_2_41'] = {'id': 'm6g.medium', 'vcpu': 1, 'memory': 4, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_42'] = {'id': 'm6g.large', 'vcpu': 2, 'memory': 8, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_43'] = {'id': 'm6g.xlarge', 'vcpu': 4, 'memory': 16, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_44'] = {'id': 'm6g.2xlarge', 'vcpu': 8, 'memory': 32, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_45'] = {'id': 'm6g.4xlarge', 'vcpu': 16, 'memory': 64, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_46'] = {'id': 'm6g.8xlarge', 'vcpu': 32, 'memory': 128, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_47'] = {'id': 'm6g.12xlarge', 'vcpu': 48, 'memory': 192, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['1_2_48'] = {'id': 'm6g.16xlarge', 'vcpu': 64, 'memory': 256, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'general purpose', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # c5 instance types
        instance_type_dict['2_1_21'] = {'id': 'c5.large', 'vcpu': 2, 'memory': 4, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_22'] = {'id': 'c5.xlarge', 'vcpu': 4, 'memory': 8, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_23'] = {'id': 'c5.2xlarge', 'vcpu': 8, 'memory': 16, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_24'] = {'id': 'c5.4xlarge', 'vcpu': 16, 'memory': 32, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_25'] = {'id': 'c5.9xlarge', 'vcpu': 36, 'memory': 72, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_26'] = {'id': 'c5.12xlarge', 'vcpu': 48, 'memory': 96, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake)', 'speed': 'up to 3.9 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_27'] = {'id': 'c5.18xlarge', 'vcpu': 72, 'memory': 144, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake) or 1st generation Intel Xeon Platinum 8000 serie (Skylake-SP)', 'speed': 'up to 3.5 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_28'] = {'id': 'c5.24xlarge', 'vcpu': 96, 'memory': 192, 'processor': '2nd generation Intel Xeon Scalable (Cascade Lake)', 'speed': 'up to 3.9 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # c5a instance types
        instance_type_dict['2_1_31'] = {'id': 'c5a.large', 'vcpu': 2, 'memory': 4, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_32'] = {'id': 'c5a.xlarge', 'vcpu': 4, 'memory': 8, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_33'] = {'id': 'c5a.2xlarge', 'vcpu': 8, 'memory': 16, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_34'] = {'id': 'c5a.4xlarge', 'vcpu': 16, 'memory': 32, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_35'] = {'id': 'c5a.8xlarge', 'vcpu': 32, 'memory': 64, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_36'] = {'id': 'c5a.12xlarge', 'vcpu': 48, 'memory': 96, 'processor': 'AMD EPYC 7002 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_37'] = {'id': 'c5a.16xlarge', 'vcpu': 64, 'memory': 128, 'processor': 'AMD EPYC 7000 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['2_1_38'] = {'id': 'c5a.24xlarge', 'vcpu': 96, 'memory': 192, 'processor': 'AMD EPYC 7000 serie', 'speed': 'up to 3.3 GHz', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
    
        # c6g instance types
        # -- instance_type_dict['2_1_41'] = {'id': 'c6g.medium', 'vcpu': 1, 'memory': 2, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_42'] = {'id': 'c6g.large', 'vcpu': 2, 'memory': 4, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_43'] = {'id': 'c6g.xlarge', 'vcpu': 4, 'memory': 8, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_44'] = {'id': 'c6g.2xlarge', 'vcpu': 8, 'memory': 16, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_45'] = {'id': 'c6g.4xlarge', 'vcpu': 16, 'memory': 32, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_46'] = {'id': 'c6g.8xlarge', 'vcpu': 32, 'memory': 64, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_47'] = {'id': 'c6g.12xlarge', 'vcpu': 48, 'memory': 96, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['2_1_48'] = {'id': 'c6g.16xlarge', 'vcpu': 64, 'memory': 128, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'compute optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # r5 instance types
        instance_type_dict['3_1_21'] = {'id': 'r5.large', 'vcpu': 2, 'memory': 16, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_22'] = {'id': 'r5.xlarge', 'vcpu': 4, 'memory': 32, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_23'] = {'id': 'r5.2xlarge', 'vcpu': 8, 'memory': 64, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_24'] = {'id': 'r5.4xlarge', 'vcpu': 16, 'memory': 128, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_25'] = {'id': 'r5.8xlarge', 'vcpu': 32, 'memory': 256, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_26'] = {'id': 'r5.12xlarge', 'vcpu': 48, 'memory': 384, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_27'] = {'id': 'r5.16xlarge', 'vcpu': 64, 'memory': 512, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_28'] = {'id': 'r5.24xlarge', 'vcpu': 96, 'memory': 768, 'processor': 'Intel Xeon Platinum 8175', 'speed': 'up to 3.1 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # r5a instance types
        instance_type_dict['3_1_31'] = {'id': 'r5a.large', 'vcpu': 2, 'memory': 16, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_32'] = {'id': 'r5a.xlarge', 'vcpu': 4, 'memory': 32, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_33'] = {'id': 'r5a.2xlarge', 'vcpu': 8, 'memory': 64, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_34'] = {'id': 'r5a.4xlarge', 'vcpu': 16, 'memory': 128, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_35'] = {'id': 'r5a.8xlarge', 'vcpu': 32, 'memory': 256, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_36'] = {'id': 'r5a.12xlarge', 'vcpu': 48, 'memory': 384, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_37'] = {'id': 'r5a.16xlarge', 'vcpu': 64, 'memory': 512, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        instance_type_dict['3_1_38'] = {'id': 'r5a.24xlarge', 'vcpu': 96, 'memory': 768, 'processor': 'AMD EPYC 7000 serie', 'speed': '2.5 GHz', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

        # r6g instance types
        # -- instance_type_dict['3_1_41'] = {'id': 'r6g.medium', 'vcpu': 1, 'memory': 8, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_42'] = {'id': 'r6g.large', 'vcpu': 2, 'memory': 16, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_43'] = {'id': 'r6g.xlarge', 'vcpu': 4, 'memory': 32, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_44'] = {'id': 'r6g.2xlarge', 'vcpu': 8, 'memory': 64, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_45'] = {'id': 'r6g.4xlarge', 'vcpu': 16, 'memory': 128, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_46'] = {'id': 'r6g.8xlarge', 'vcpu': 32, 'memory': 256, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_47'] = {'id': 'r6g.12xlarge', 'vcpu': 48, 'memory': 384, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}
        # -- instance_type_dict['3_1_48'] = {'id': 'r6g.16xlarge', 'vcpu': 64, 'memory': 512, 'processor': 'Arm-based AWS Graviton2', 'speed': 'N/A', 'use': 'memory optimized', 'nitro': 'yes', 'starcluster': 'non-supported', 'generation': 'current'}

    # return the volumes data dictionary
    return instance_type_dict

#-------------------------------------------------------------------------------

def get_instance_type_data_dict(instance_type):
    '''
    Get the data dictionary of a instance type.
    '''

    # initialize the data dictionary
    data_dict = {}

    # get the instance type dictionary
    instance_type_dict = get_instance_type_dict(get_cluster_mode_native())

    for key in instance_type_dict.keys():
        if instance_type_dict[key]['id'] == instance_type:
            data_dict = instance_type_dict[key]
            break

    # return the data dictionary
    return data_dict

#-------------------------------------------------------------------------------

def create_ngscloud_config_file(user_id, access_key_id, secret_access_key, email):
    '''
    Create the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the default current region and zone names
    region_name = get_default_region_name()
    zone_name = get_default_zone_name()

    # get NGScloud Key
    ngscloud_key = get_ngscloud_key()

    # get the instance type dictionary
    instance_type_dict = get_instance_type_dict(get_cluster_mode_starcluster())

    # get the instance type key list
    instance_type_key_list = list(instance_type_dict.keys())
    instance_type_key_list.sort()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # initialize the instance type dictionary
    instance_type_list_in_static_py = xstarcluster.get_instance_type_list_in_static_py()

    # create the NGScloud config file corresponding to the current environment
    if OK:
        try:
            if not os.path.exists(os.path.dirname(ngscloud_config_file)):
                os.makedirs(os.path.dirname(ngscloud_config_file))
            with open(ngscloud_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
                file_id.write( '[global]\n')
                file_id.write(f'default_template = {environment}-t2.micro\n')
                file_id.write(f'environment = {environment}\n')
                file_id.write(f'current_region = {region_name}\n')
                file_id.write(f'current_zone = {zone_name}\n')
                file_id.write( '\n')
                file_id.write( '[aws info]\n')
                file_id.write(f'aws_user_id = {user_id}\n')
                file_id.write(f'aws_access_key_id = {access_key_id}\n')
                file_id.write(f'aws_secret_access_key = {secret_access_key}\n')
                file_id.write( '\n')
                file_id.write( '[contact info]\n')
                file_id.write(f'email = {email}\n')
                file_id.write( '\n')
                file_id.write(f'[key {ngscloud_key}]\n')
                file_id.write(f'key_location = {xlib.get_keypairs_dir()}/{ngscloud_key}-{user_id}-{region_name}.rsa\n')
                file_id.write( '\n')
                file_id.write( '[dataset info]\n')
                file_id.write(f'dataset_structure = {get_dataset_structure_none()}\n')
                file_id.write( 'ngscloud_volume =\n')
                file_id.write( 'app_volume =\n')
                file_id.write( 'database_volume =\n')
                file_id.write( 'read_volume =\n')
                file_id.write( 'reference_volume =\n')
                file_id.write( 'result_volume =\n')
                file_id.write( '\n')
                for instance_type_key in instance_type_key_list:
                    id = instance_type_dict[instance_type_key]['id']
                    vcpu = instance_type_dict[instance_type_key]['vcpu']
                    memory = instance_type_dict[instance_type_key]['memory']
                    processor = instance_type_dict[instance_type_key]['processor']
                    use = instance_type_dict[instance_type_key]['use']
                    starcluster = instance_type_dict[instance_type_key]['starcluster']
                    if id in instance_type_list_in_static_py:
                        file_id.write(f'[cluster {environment}-{id}]\n')
                        file_id.write(f'description = {id} - {vcpu} vCPU - {memory} GiB - {use}\n')
                        file_id.write(f'vcpu = {vcpu}\n')
                        file_id.write(f'memory = {memory}\n')
                        file_id.write(f'use = {use}\n')
                        file_id.write(f'processor = {processor}\n')
                        file_id.write(f'starcluster = {starcluster}\n')
                        file_id.write(f'keyname = {ngscloud_key}\n')
                        file_id.write( 'cluster_size = 1\n')
                        file_id.write( 'cluster_user = sgeadmin\n')
                        file_id.write( 'cluster_shell = bash\n')
                        file_id.write( 'master_image_id = ami-6b211202\n')
                        file_id.write(f'master_instance_type = {id}\n')
                        file_id.write( 'node_image_id = ami-6b211202\n')
                        file_id.write(f'node_instance_type = {id}\n')
                        file_id.write( 'volumes = \n')
                        file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_miniconda3_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write( 'url = https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh\n')
                file_id.write( 'channels = N/A\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_bcftools_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_bedtools_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_blastplus_name()}]\n')
                file_id.write( 'version = 2.9.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_bowtie2_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_busco_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda;conda-forge\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_cd_hit_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_cufflinks_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_cutadapt_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_ddradseqtools_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write( 'url = https://github.com/GGFHF/ddRADseqTools/archive/master.zip\n')
                file_id.write( 'channels = N/A\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_detonate_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda;conda-forge\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_diamond_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_emboss_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_entrez_direct_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_express_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_fastqc_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_gmap_gsnap_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_hisat2_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_htseq_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_ipyrad_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_kallisto_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_ngshelper_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write( 'url = https://github.com/GGFHF/NGShelper/archive/master.zip\n')
                file_id.write( 'channels = N/A\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_quast_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_r_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_raddesigner_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write( 'url = https://github.com/GGFHF/RADdesigner/archive/master.zip\n')
                file_id.write( 'channels = N/A\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_rnaquast_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_rsem_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_samtools_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_soapdenovo2_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_soapdenovotrans_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_star_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_starcode_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_tabix_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_toa_name()}]\n')
                file_id.write( 'version = last\n')
                file_id.write(f'url = https://github.com/GGFHF/TOA/archive/master.zip\n')
                file_id.write( 'channels = N/A\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_tophat_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_transabyss_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda;conda-forge\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_transdecoder_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_transrate_name()}]\n')
                file_id.write( 'version = 1.0.3\n')
                file_id.write(f'url = https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz\n')
                file_id.write( 'channels = N/A\n')
                # -- file_id.write( 'version = any\n')
                # -- file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                # -- file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                # -- file_id.write(f'[bioinfoapp {xlib.get_transrate_tools_name()}]\n')
                # -- file_id.write( 'version = any\n')
                # -- file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                # -- file_id.write( 'channels = bioconda\n')
                # -- file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_trimmomatic_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_trinity_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_vcftools_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_vcftools_perl_libraries_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write(f'[bioinfoapp {xlib.get_vsearch_name()}]\n')
                file_id.write( 'version = any\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-essentials]\n')
                file_id.write( 'version = 3.6\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware bioconductor-phyloseq]\n')
                file_id.write( 'version = 1.30.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda;conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-ade4]\n')
                file_id.write( 'version = 1.7_15\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-adegenet]\n')
                file_id.write( 'version = 2.1.2\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-adegraphics]\n')
                file_id.write( 'version = 1.0_15\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = r\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-devtools]\n')
                file_id.write( 'version = 2.3.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-ggmap]\n')
                file_id.write( 'version = 3.0.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-ggplot2]\n')
                file_id.write( 'version = 3.3.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-gplots]\n')
                file_id.write( 'version = 3.0.3\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-gridextra]\n')
                file_id.write( 'version = 2.3\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-rlist]\n')
                file_id.write( 'version = 0.4.6.1\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-sparkr]\n')
                file_id.write( 'version = 2.4.4\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = r\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-sqldf]\n')
                file_id.write( 'version = 0.4_11\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-stampp]\n')
                file_id.write( 'version = 1.5.1\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = bioconda;conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-tidyverse]\n')
                file_id.write( 'version = 1.3.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware r-vcfr]\n')
                file_id.write( 'version = 1.10.0\n')
                file_id.write(f'url = {xlib.get_anaconda_url()}\n')
                file_id.write( 'channels = conda-forge\n')
                file_id.write( '\n')
                file_id.write( '[rsoftware easyGgplot2]\n')
                file_id.write( 'version = last\n')
                file_id.write(f'url = {xlib.get_github_url()}\n')
                file_id.write( 'channels = kassambara/easyGgplot2\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {ngscloud_config_file} can not be created')
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def update_connection_data(user_id, access_key_id, secret_access_key):
    '''
    Update the user id, access key id and secret access key in the NGScloud
    config file corresponding to the environment.
    '''

    # initialize the control variable
    OK = True

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # update the connection data in the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict['aws info']['aws_user_id'] = user_id
    ngscloud_options_dict['aws info']['aws_access_key_id'] = access_key_id
    ngscloud_options_dict['aws info']['aws_secret_access_key'] = secret_access_key

    # save the option dictionary in the NGScloud config file corresponding to the environment
    (OK, error_list) = save_ngscloud_config_file(ngscloud_options_dict)

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def update_contact_data(email):
    '''
    Update the contact e-mail in the NGScloud config file corresponding to the
    environment.
    '''

    # initialize the control variable
    OK = True

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # update the contact e-mail in the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict['contact info']['email'] = email

    # save the option dictionary in the NGScloud config file corresponding to the environment
    (OK, error_list) = save_ngscloud_config_file(ngscloud_options_dict)

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def update_region_zone_data(region_name, zone_name):
    '''
    Update the current region and zone names in the NGScloud config file
    corresponding to the environment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get NGScloud Key
    ngscloud_key = get_ngscloud_key()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the old region and user identification
    old_region_name = ngscloud_options_dict['global']['current_region']
    user_id = ngscloud_options_dict['aws info']['aws_user_id']

    # chek the Ubuntu AMI identification of the new region
    ubuntu_ami_id = xec2.get_ubuntu_ami_id(region_name)
    if ubuntu_ami_id == xec2.get_unknown_ami_id():
        error_list.append(f'*** ERROR: The AMI {xec2.get_ubuntu_ami_name()} is not found in the region {region_name}. The region and zone are not modified.')
        OK = False

    # get the StarCluster AMI identification of the new region
    if OK:
        starcluster_ami_id = xec2.get_starcluster_ami_id(region_name)
        if starcluster_ami_id == xec2.get_unknown_ami_id():
            error_list.append(f'*** WARNING: The AMI {xec2.get_starcluster_ami_name()} is not found in the region {region_name}. The cluster mode StarCluster cannot be used.')

    # update the current region and zone in the option dictionary corresponding to the NGScloud config file
    if OK:
        ngscloud_options_dict['global']['current_region'] = region_name
        ngscloud_options_dict['global']['current_zone'] = zone_name
        ngscloud_options_dict[f'key {xlib.get_project_name()}Key']['key_location'] = f'{xlib.get_keypairs_dir()}/{ngscloud_key}-{user_id}-{region_name}.rsa'

    # update the master and node image identification and inicialize the asigned volumes in each cluster template
    if OK:
        if region_name != old_region_name:
            template_name_list = get_template_name_list()
            for template_name in template_name_list:
                ngscloud_options_dict['cluster {}'.format(template_name)]['master_image_id'] = starcluster_ami_id
                ngscloud_options_dict['cluster {}'.format(template_name)]['node_image_id'] = starcluster_ami_id
                ngscloud_options_dict['cluster {}'.format(template_name)]['volumes'] = ''

    # save the option dictionary in the NGScloud config file
    if OK:
        (OK, error_list2) = save_ngscloud_config_file(ngscloud_options_dict)
        error_list = error_list + error_list2

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def link_volumes(dataset_structure, ngscloud_volume, app_volume, database_volume, read_volume, reference_volume, result_volume, log, function=None):
    '''
    Link volumes in the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('Do not close this window, please wait!\n')

    # get current zone name
    zone_name = get_current_zone_name()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the defined cluster templates list
    template_names_list = get_template_name_list()

    # initialize the linked volume name list
    linked_volume_list = []

    # set dataset info data in the option dictionary of NGScloud config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'The dataset info data are going to be updated ...\n')
    ngscloud_options_dict['dataset info']['dataset_structure'] = dataset_structure
    ngscloud_options_dict['dataset info']['ngscloud_volume'] = ngscloud_volume
    ngscloud_options_dict['dataset info']['app_volume'] = app_volume
    ngscloud_options_dict['dataset info']['database_volume'] = database_volume
    ngscloud_options_dict['dataset info']['read_volume'] = read_volume
    ngscloud_options_dict['dataset info']['reference_volume'] = reference_volume
    ngscloud_options_dict['dataset info']['result_volume'] = result_volume
    log.write(f'Data are updated.\n')

    # remove previous volume sections
    previous_volume_section_list = []
    for section in ngscloud_options_dict.keys():
        if section.startswith('volume'):
            previous_volume_section_list.append(section)
    for section in previous_volume_section_list:
        del(ngscloud_options_dict[section])

    # set the NGScloud volume section in the option dictionary of NGScloud config file
    if ngscloud_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {ngscloud_volume} data are going to be added ...\n')
        linked_volume_list.append(ngscloud_volume)
        mount_path = xlib.get_cluster_ngscloud_dir()
        aws_device_file = xlib.get_cluster_ngscloud_device_file()
        ngscloud_options_dict[f'volume {ngscloud_volume}'] = {'volume_id': xec2.get_volume_id(ngscloud_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # set the application volume data in the option dictionary of NGScloud config file
    if app_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {app_volume} data are going to be added ...\n')
        linked_volume_list.append(app_volume)
        mount_path = xlib.get_cluster_app_dir()
        aws_device_file = xlib.get_cluster_app_device_file()
        ngscloud_options_dict[f'volume {app_volume}'] = {'volume_id': xec2.get_volume_id(app_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # set the database volume data in the option dictionary of NGScloud config file
    if database_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {database_volume} data are going to be added ...\n')
        linked_volume_list.append(database_volume)
        mount_path = xlib.get_cluster_database_dir()
        aws_device_file = xlib.get_cluster_database_device_file()
        ngscloud_options_dict[f'volume {database_volume}'] = {'volume_id': xec2.get_volume_id(database_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # set the read volume data in the option dictionary of NGScloud config file
    if read_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {read_volume} data are going to be added ...\n')
        linked_volume_list.append(read_volume)
        mount_path = xlib.get_cluster_read_dir()
        aws_device_file = xlib.get_cluster_read_device_file()
        ngscloud_options_dict[f'volume {read_volume}'] = {'volume_id': xec2.get_volume_id(read_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # set the reference volume data in the option dictionary of NGScloud config file
    if reference_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {reference_volume} data are going to be added ...\n')
        linked_volume_list.append(reference_volume)
        mount_path = xlib.get_cluster_reference_dir()
        aws_device_file = xlib.get_cluster_reference_device_file()
        ngscloud_options_dict[f'volume {reference_volume}'] = {'volume_id': xec2.get_volume_id(reference_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # set the result volume data in the option dictionary of NGScloud config file
    if result_volume != '':
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The volume {result_volume} data are going to be added ...\n')
        linked_volume_list.append(result_volume)
        mount_path = xlib.get_cluster_result_dir()
        aws_device_file = xlib.get_cluster_result_device_file()
        ngscloud_options_dict[f'volume {result_volume}'] = {'volume_id': xec2.get_volume_id(result_volume, zone_name), 'mount_path': mount_path, 'aws_device_file': aws_device_file}
        log.write('Data are added.\n')

    # add linked volume list to every template in the tamplate list
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Volumes are linked to cluster templates used by StarCluster ...\n')
        for template_name in template_names_list:
            ngscloud_options_dict[f'cluster {template_name}']['volumes'] = ', '.join(linked_volume_list)
        log.write('Volume are linked.\n')

    # save the option dictionary in the NGScloud config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'The file {ngscloud_config_file} is going to be saved ...\n')
        (OK, error_list) = save_ngscloud_config_file(ngscloud_options_dict)
        if OK:
            log.write('The config file has been saved.\n')
        else:
            for error in error_list:
                log.write(error)

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return (OK, error_list)

#-------------------------------------------------------------------------------

def save_ngscloud_config_file(ngscloud_options_dict):
    '''
    Save the NGScloud config options in the config file corresponding
    to the environment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the sections list with the new values
    sections_list = []
    for section in ngscloud_options_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build the section "global"
    section = 'global'
    config[section] = {}
    for k, v in ngscloud_options_dict[section].items():
        config[section][k] = v

    # build the section "aws info"
    section = 'aws info'
    config[section] = {}
    for k, v in ngscloud_options_dict[section].items():
        config[section][k] = v

    # build the section "contact info"
    section = 'contact info'
    config[section] = {}
    for k, v in ngscloud_options_dict[section].items():
        config[section][k] = v

    # build the sections "key *"
    for section in sections_list:
        if section.startswith('key'):
            config[section] = {}
            for k, v in ngscloud_options_dict[section].items():
                config[section][k] = v

    # build the section "dataset info"
    section = 'dataset info'
    config[section] = {}
    for k, v in ngscloud_options_dict[section].items():
        config[section][k] = v

    # build the sections "volume *"
    for section in sections_list:
        if section.startswith('volume'):
            config[section] = {}
            for k, v in ngscloud_options_dict[section].items():
                config[section][k] = v

    # build the sections "cluster *"
    for section in sections_list:
        if section.startswith('cluster'):
            config[section] = {}
            for k, v in ngscloud_options_dict[section].items():
                config[section][k] = v

    # build the sections "bioinfoapp *"
    for section in sections_list:
        if section.startswith('bioinfoapp'):
            config[section] = {}
            for k, v in ngscloud_options_dict[section].items():
                config[section][k] = v

    # build the sections "rsoftware *"
    for section in sections_list:
        if section.startswith('rsoftware'):
            config[section] = {}
            for k, v in ngscloud_options_dict[section].items():
                config[section][k] = v

    # write the NGScloud config file
    try:
        with open(ngscloud_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            config.write(file_id)
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {ngscloud_config_file} can not be written')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------
    
def get_default_region_name():
    '''
    Get the region name by default.
    '''

    return 'us-east-1'

#-------------------------------------------------------------------------------
    
def get_default_zone_name():
    '''
    Get the region name by default.
    '''

    return 'us-east-1c'

#-------------------------------------------------------------------------------

def get_basic_aws_data():
    '''
    Get the connection data to AWS: access the user identification, the key
    identification and the secret access key from NGScloud config file.
    '''

    # create class to parse the config file
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the connection data from NGScloud config file
    user_id= config.get('aws info', 'aws_user_id', fallback='')
    access_key_id = config.get('aws info', 'aws_access_key_id', fallback='')
    secret_access_key = config.get('aws info', 'aws_secret_access_key', fallback='')

    # return the connection data
    return (user_id, access_key_id, secret_access_key)

#-------------------------------------------------------------------------------

def get_contact_data():
    '''
    Get the access the contact e-mail from NGScloud config file.
    '''

    # create class to parse the config file
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the contact e-mail data from the NGScloud config file
    email = config.get('contact info', 'email', fallback='')

    # return the contact e-mail
    return email

#-------------------------------------------------------------------------------

def get_key_sections_dict():
    '''
    Get the key sections data dictionary from the NGScloud config file
    corresponding to the environment.
    '''

    # initialize the key sections data dictionary
    key_sections_dict = {}
    
    # create class to parse the config files
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the sections list
    sections_list = []
    for section in config.sections():
        sections_list.append(section)
    sections_list.sort()

    # build the key sections data dictionary
    for section in sections_list:
        section_type = 'key'
        if section.startswith(section_type):
            key_section_name = section[len(section_type) + 1:]
            key_sections_dict[key_section_name] = {}
            for k, v in config[section].items():
                key_sections_dict[key_section_name][k] = v

    # return the key sections data dictionary
    return key_sections_dict

#-------------------------------------------------------------------------------

def get_template_dict():
    '''
    Get the cluster templates data dictionary from the NGScloud config file
    corresponding to the environment.
    '''

    # initialize the key sections data dictionary
    template_dict = {}
    
    # create class to parse the config files
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the sections list
    sections_list = []
    for section in config.sections():
        sections_list.append(section)
    sections_list.sort()

    # build the cluster template data dictionary
    for section in sections_list:
        section_type = 'cluster'
        if section.startswith(section_type):
            key_section_name = section[len(section_type) + 1:]
            template_dict[key_section_name] = {}
            template_dict[key_section_name]['template_name'] = key_section_name
            for k, v in config[section].items():
                template_dict[key_section_name][k] = v

    # return the cluster templates  data dictionary
    return template_dict

#-------------------------------------------------------------------------------

def get_template_name_list():
    '''
    Get the template name list from the NGScloud config file corresponding
    to the environment.
    '''

    # initialize the template name list
    template_name_list = []

    # get the instance type dictionary
    instance_type_dict = get_instance_type_dict(get_cluster_mode_starcluster())

    # get the instance type key list
    instance_type_key_list = list(instance_type_dict.keys())
    instance_type_key_list.sort()
    
    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # build the template name list ordered
    for instance_type_key in instance_type_key_list:
        id = instance_type_dict[instance_type_key]['id']
        template_name = f'{environment}-{id}'
        if f'cluster {template_name}' in config.sections():
            template_name_list.append(template_name)

    # return the template name list
    return template_name_list

#-------------------------------------------------------------------------------

def get_volumes_dict():
    '''
    Get the volumes data dictionary from the NGScloud config file
    corresponding to the environment.
    '''

    # initialize the key sections data dictionary
    volumes_dict = {}
    
    # create class to parse the config files
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the sections list
    sections_list = []
    for section in config.sections():
        sections_list.append(section)
    sections_list.sort()

    # build the volumes data dictionary
    for section in sections_list:
        section_type = 'volume'
        if section.startswith(section_type):
            key_section_name = section[len(section_type) + 1:]
            volumes_dict[key_section_name] = {}
            for k, v in config[section].items():
                volumes_dict[key_section_name][k] = v

    # return the volumes data dictionary
    return volumes_dict

#-------------------------------------------------------------------------------

def get_volume_names_list():
    '''
    Get the list of the volume names from the NGScloud config file
    corresponding to the environment.
    '''

    # get the volume names list
    volume_names_list = get_volumes_dict().keys()

    # return the volume names list
    return sorted(volume_names_list)

#-------------------------------------------------------------------------------

def get_current_region_name():
    '''
    Get the current region name from the NGScloud config file corresponding
    to the environment.
    '''

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # get the current region name
    current_region_name = config.get('global', 'current_region', fallback='')

    # return  the current region name
    return current_region_name

#-------------------------------------------------------------------------------

def get_current_zone_name():
    '''
    Get the current zone name from the NGScloud config file corresponding to
    the environment.
    '''

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # get the current zone name
    current_zone_name = config.get('global', 'current_zone', fallback='')

    # return  the current zone name
    return current_zone_name

#-------------------------------------------------------------------------------

def get_keypair_file():
    '''
    Get the path of the key pair file to the current region from the NGScloud
    config file corresponding to the environment.
    '''

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # get the option dictionary
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the key pair file
    keypair_file = ngscloud_options_dict[f'key {xlib.get_project_name()}Key']['key_location']

    # return the key pair file
    return keypair_file

#-------------------------------------------------------------------------------

def get_bioinfo_app_data(bioinfo_app_name):
    '''
    Get the data of a bioinfo application from the NGScloud
    config file corresponding to the environment.
    '''

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # get the data
    bioinfo_app_version = config.get(f'bioinfoapp {bioinfo_app_name}', 'version', fallback='')
    bioinfo_app_url = config.get(f'bioinfoapp {bioinfo_app_name}', 'url', fallback='')
    bioinfo_app_channels = config.get(f'bioinfoapp {bioinfo_app_name}', 'channels', fallback='')

    # return the data
    return (bioinfo_app_version, bioinfo_app_url, bioinfo_app_channels)

#-------------------------------------------------------------------------------

def get_r_software_dict():
    '''
    Get the R software data dictionary from the NGScloud config file
    corresponding to the environment.
    '''

    # initialize the key sections data dictionary
    r_software_dict = {}
    
    # create class to parse the config files
    config = configparser.ConfigParser()

    # read the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()
    config.read(ngscloud_config_file)

    # get the sections list
    sections_list = []
    for section in config.sections():
        sections_list.append(section)
    sections_list.sort()

    # build the R software data dictionary
    for section in sections_list:
        section_type = 'rsoftware'
        if section.startswith(section_type):
            key_section_name = section[len(section_type) + 1:]
            r_software_dict[key_section_name] = {}
            for k, v in config[section].items():
                r_software_dict[key_section_name][k] = v

    # return the R software data dictionary
    return r_software_dict

#-------------------------------------------------------------------------------

def get_r_software_data(r_software_name):
    '''
    Get the data of a R software from the NGScloud
    config file corresponding to the environment.
    '''

    # create class to parse the config files
    config = configparser.ConfigParser()

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # read the NGScloud config file
    config.read(ngscloud_config_file)

    # get the data
    r_software_version = config.get(f'rsoftware {r_software_name}', 'version', fallback='')
    r_software_url = config.get(f'rsoftware {r_software_name}', 'url', fallback='')
    r_software_channels = config.get(f'rsoftware {r_software_name}', 'channels', fallback='')

    # return the data
    return (r_software_version, r_software_url, r_software_channels)

#-------------------------------------------------------------------------------

def is_template_defined(template_name):
    '''
    Check if a template name is defined in the NGScloud config file corresponding
    to the environment.
    '''

    # initialize of the control variable
    found = False

    # find the template definition
    if template_name in get_template_name_list():
        found = True

    # return the control variable
    return found

#-------------------------------------------------------------------------------

def is_ngscloud_config_file_created():
    '''
    Check if the NGScloud config file corresponding to the environment is created.
    '''

    # initialize the found variable
    found = True

    # check if NGScloud config file is created
    if not os.path.isfile(get_ngscloud_config_file()):
        found = False

    # return the found variable
    return found

#-------------------------------------------------------------------------------

def set_environment_variables():
    '''
    Set the environment variables corresponding to the NGScloud config file
    corresponding to the environment: the AWS access key identification, AWS
    secret access key and the current region name
    '''

    # get the NGScloud config file
    ngscloud_config_file = get_ngscloud_config_file()

    # set environment variable TRANCRIPTOMECLOUD_CONFIG_FILE
    os.environ[f'{xlib.get_project_code().upper()}_CONFIG_FILE'] = ngscloud_config_file

#-------------------------------------------------------------------------------

def get_ngscloud_config_file():
    '''
    Get the path of the NGScloud config file corresponding to the environment.
    '''

    # assign the NGScloud config file
    ngscloud_config_file = f'{xlib.get_config_dir()}/{environment}-{xlib.get_project_code()}-config.txt'

    # return the NGScloud config file
    return ngscloud_config_file

#-------------------------------------------------------------------------------

def get_ngscloud_key():
    '''
    Get the code of the NGScloud key.
    '''

    # assign the NGScloud config file
    ngscloud_key = f'{xlib.get_project_name()}Key'

    return ngscloud_key

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print(f'This file contains the functions related to the {xlib.get_project_name()} configuration used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
