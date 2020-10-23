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
This source contains the class Main corresponding to the graphical user interface of
the NGScloud software package.
'''

#-------------------------------------------------------------------------------

import os
import PIL.Image
import PIL.ImageTk
import sys
import threading
import time
import tkinter
import tkinter.ttk
import tkinter.messagebox
import webbrowser

import gbioinfoapp
import gcloud
import gdataset
import gdialogs
import glog
import gtoa
import xbowtie2
import xbusco
import xcdhit
import xcluster
import xconfiguration
import xcufflinks
import xcutadapt
import xdatabase
import xddradseqtools
import xdetonate
import xec2
import xexpress
import xfastqc
import xgmap
import xgzip
import xhisat2
import xhtseq
import xinstance
import xipyrad
import xkallisto
import xlib
import xngshelper
import xquast
import xraddesigner
import xsoapdenovo2
import xsoapdenovotrans
import xstar
import xstarcode
import xread
import xreference
import xresult
import xrnaquast
import xtoa
import xtophat
import xtransabyss
import xtransrate
import xtrimmomatic
import xtrinity
import xvolume

#-------------------------------------------------------------------------------

class Main():

    #---------------

    if sys.platform.startswith('linux'):
        WINDOW_HEIGHT = 620
        WINDOW_WIDTH = 925
    elif sys.platform.startswith('darwin'):
        WINDOW_HEIGHT = 650
        WINDOW_WIDTH = 1080
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        WINDOW_HEIGHT = 590
        WINDOW_WIDTH = 870

    #---------------

    def __init__(self):
        '''
        Execute actions correspending to the creation of a "App" instance.
        '''

        # call the init method of the parent class
        self.root = tkinter.Tk()

        # create the window
        self.create_window()

        # build the graphical user interface
        self.build_gui()
        # self.root.grid()

        # initialize the forms dictionary
        self.forms_dict = {}

        # create "form_welcome" and register it in "container" with the grid geometry manager
        self.form_welcome = FormWelcome(self)
        self.form_welcome.grid(row=0, column=0, sticky='nsew')

        # add "form_welcome" in the forms dictionary
        self.forms_dict['form_welcome'] = self.form_welcome

        # create and register "form_set_environment" in "container" with the grid geometry manager
        self.form_set_environment = gcloud.FormSetEnvironment(self)
        self.form_set_environment.grid(row=0, column=0, sticky='nsew')

        # set "form_set_environment" as current form and add it in the forms dictionary
        self.current_form = 'form_set_environment'
        self.forms_dict[self.current_form] = self.form_set_environment

        # raise "form_set_environment" to front
        self.form_set_environment.tkraise()

    #---------------

    def create_window(self):
        '''
        Create the window of "Main".
        '''

        # define the dimensions
        x = round((self.root.winfo_screenwidth() - self.WINDOW_WIDTH) / 2)
        y = round((self.root.winfo_screenheight() - self.WINDOW_HEIGHT) / 2)
        self.root.geometry('{}x{}+{}+{}'.format(self.WINDOW_WIDTH, self.WINDOW_HEIGHT, x, y))
        self.root.minsize(height=self.WINDOW_HEIGHT, width=self.WINDOW_WIDTH)
        self.root.maxsize(height=self.WINDOW_HEIGHT, width=self.WINDOW_WIDTH)

        # set default fondt
        if sys.platform.startswith('linux'):
            self.root.option_add('*Font', 'Verdana 10')
        elif sys.platform.startswith('darwin'):
            self.root.option_add('*Font', 'Verdana 10')
        elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
            self.root.option_add('*Font', 'Segoe 10')

        # set default language of MessageBox
        self.root.tk.eval('::msgcat::mclocale en')

        # set the title
        self.root.title(xlib.get_project_name())

        # set the icon
        image_app = PIL.Image.open(xlib.get_project_image_file())
        self.photoimage_app = PIL.ImageTk.PhotoImage(image_app)
        self.root.tk.call('wm', 'iconphoto', self.root._w, self.photoimage_app)

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "Main".
        '''

        # create "imagetk_exit"
        image_exit = PIL.Image.open('./image_exit.png')
        imagetk_exit = PIL.ImageTk.PhotoImage(image_exit)  

        # maximize the width of column 0
        self.root.grid_columnconfigure(0, weight=1)

        # create "menu_bar"
        self.menu_bar = tkinter.Menu(self.root)

        # create "menu_system" and add its menu items and links with "menu_configuration" and "menu_security"
        self.menu_system = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_system.add_command(label='Exit', command=self.exit, accelerator='Alt+F4', compound='left', image=imagetk_exit)

        # link "menu_system" to "menu_bar"
        self.menu_bar.add_cascade(label='System', menu=self.menu_system)

        # create "menu_configuration" and add its menu items
        self.menu_configuration = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_configuration.add_command(label=f'Recreate {xlib.get_project_name()} config file', command=self.recreate_ngscloud_config_file)
        self.menu_configuration.add_command(label=f'View {xlib.get_project_name()} config file', command=self.view_ngscloud_config_file)
        self.menu_configuration.add_separator()
        self.menu_configuration.add_command(label='List instance types', command=self.list_instance_types)
        self.menu_configuration.add_separator()
        self.menu_configuration.add_command(label='Update connection data and contact e-mail', command=self.update_connection_data)
        self.menu_configuration.add_command(label='Update region and zone', command=self.update_region_zone)
        self.menu_configuration.add_separator()
        self.menu_configuration.add_command(label='Link volumes', command=self.link_volumes)

        # create "menu_security" and add its menu items
        self.menu_security = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_security.add_command(label='List key pairs', command=self.list_keypairs)
        self.menu_security.add_separator()
        self.menu_security.add_command(label='Create key pairs', command=self.create_keypairs)
        self.menu_security.add_separator()
        self.menu_security.add_command(label='List cluster security groups', command=self.warn_unavailable_process)
        self.menu_security.add_separator()
        self.menu_security.add_command(label='Force removal of a cluster security group', command=self.warn_unavailable_process)

        # create "menu_cluster_operation" add add its menu items
        self.menu_cluster_operation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cluster_operation.add_command(label='List clusters', command=self.list_clusters)
        self.menu_cluster_operation.add_separator()
        self.menu_cluster_operation.add_command(label='Create cluster', command=self.create_cluster)
        self.menu_cluster_operation.add_command(label='Terminate cluster', command=self.terminate_cluster)
        self.menu_cluster_operation.add_separator()
        self.menu_cluster_operation.add_command(label='Force termination of a cluster', command=self.force_cluster_termination)
        self.menu_cluster_operation.add_separator()
        self.menu_cluster_operation.add_command(label='Show cluster composition', command=self.show_cluster_composition)
        self.menu_cluster_operation.add_separator()
        self.menu_cluster_operation.add_command(label='Show status of batch jobs', command=self.show_status_batch_jobs)
        self.menu_cluster_operation.add_command(label='Kill batch job', command=self.kill_batch_job)

        # create "menu_node_operation" add add its menu items
        self.menu_node_operation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_node_operation.add_command(label='List nodes', command=self.list_nodes)
        self.menu_node_operation.add_separator()
        self.menu_node_operation.add_command(label='Add node in a cluster', command=self.add_node)
        self.menu_node_operation.add_command(label='Remove node in acluster', command=self.remove_node)

        # create "menu_volume_operation" add add its menu items
        self.menu_volume_operation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_volume_operation.add_command(label='List volumes', command=self.list_volumes)
        self.menu_volume_operation.add_separator()
        self.menu_volume_operation.add_command(label='Create volume', command=self.create_volume)
        self.menu_volume_operation.add_command(label='Remove volume', command=self.remove_volume)
        self.menu_volume_operation.add_separator()
        self.menu_volume_operation.add_command(label='Terminate volume creator', command=self.terminate_volume_creator)
        self.menu_volume_operation.add_separator()
        self.menu_volume_operation.add_command(label='Mount volume in a node', command=self.mount_volume)
        self.menu_volume_operation.add_command(label='Unmount volume in a node', command=self.unmount_volume)

        # create "menu_bioinfo_software_installation_a_m" add add its menu items
        self.menu_bioinfo_software_installation_a_m = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_bcftools_name(), command=self.install_bcftools)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_bedtools_name(), command=self.install_bedtools)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_blastplus_name(), command=self.install_blastplus)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_bowtie2_name(), command=self.install_bowtie2)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_busco_name(), command=self.install_busco)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_cd_hit_name(), command=self.install_cd_hit)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_cufflinks_name(), command=self.install_cufflinks)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_cutadapt_name(), command=self.install_cutadapt)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_ddradseqtools_name(), command=self.install_ddradseqtools)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_detonate_name(), command=self.install_detonate)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_diamond_name(), command=self.install_diamond)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_entrez_direct_name(), command=self.install_entrez_direct)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_express_name(), command=self.install_express)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_fastqc_name(), command=self.install_fastqc)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_gmap_gsnap_name(), command=self.install_gmap_gsnap)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_hisat2_name(), command=self.install_hisat2)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_htseq_name(), command=self.install_htseq)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_ipyrad_name(), command=self.install_ipyrad)
        self.menu_bioinfo_software_installation_a_m.add_command(label=xlib.get_kallisto_name(), command=self.install_kallisto)

        # create "menu_bioinfo_software_installation_n_z" add add its menu items
        self.menu_bioinfo_software_installation_n_z = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_ngshelper_name(), command=self.install_ngshelper)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_quast_name(), command=self.install_quast)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_raddesigner_name(), command=self.install_raddesigner)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_rnaquast_name(), command=self.install_rnaquast)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_rsem_name(), command=self.install_rsem)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_samtools_name(), command=self.install_samtools)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_soapdenovo2_name(), command=self.install_soapdenovo2)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_soapdenovotrans_name(), command=self.install_soapdenovotrans)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_star_name(), command=self.install_star)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_starcode_name(), command=self.install_starcode)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_toa_name(), command=self.install_toa)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_tophat_name(), command=self.install_tophat)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_transabyss_name(), command=self.install_transabyss)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_transdecoder_name(), command=self.install_transdecoder)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_transrate_name(), command=self.install_transrate)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_trimmomatic_name(), command=self.install_trimmomatic)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_trinity_name(), command=self.install_trinity)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_vcftools_name(), command=self.install_vcftools)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_vcftools_perl_libraries_name(), command=self.install_vcftools_perl_libraries)
        self.menu_bioinfo_software_installation_n_z.add_command(label=xlib.get_vsearch_name(), command=self.install_vsearch)

        # create "menu_bioinfo_software_installation" add add its menu items
        self.menu_bioinfo_software_installation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_bioinfo_software_installation.add_command(label=f'{xlib.get_miniconda3_name()} (Bioconda infrastructure)', command=self.install_miniconda3)
        self.menu_bioinfo_software_installation.add_separator()
        self.menu_bioinfo_software_installation.add_cascade(label='A-M bioinfo software', menu=self.menu_bioinfo_software_installation_a_m)
        self.menu_bioinfo_software_installation.add_cascade(label='N-Z bioinfo software', menu=self.menu_bioinfo_software_installation_n_z)
        self.menu_bioinfo_software_installation.add_separator()
        self.menu_bioinfo_software_installation.add_command(label=f'{xlib.get_r_name()} & analysis packages', command=self.install_r)

        # create "menu_cloud_control" and add its menu items
        self.menu_cloud_control = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cloud_control.add_command(label='Set environment', command=self.set_environment)
        self.menu_cloud_control.add_separator()
        self.menu_cloud_control.add_cascade(label='Configuration', menu=self.menu_configuration)
        self.menu_cloud_control.add_cascade(label='Security', menu=self.menu_security)
        self.menu_cloud_control.add_separator()
        self.menu_cloud_control.add_cascade(label='Cluster operation', menu=self.menu_cluster_operation)
        self.menu_cloud_control.add_cascade(label='Node operation', menu=self.menu_node_operation)
        self.menu_cloud_control.add_cascade(label='Volume operation', menu=self.menu_volume_operation)
        self.menu_cloud_control.add_separator()
        self.menu_cloud_control.add_cascade(label='Bioinfo software installation', menu=self.menu_bioinfo_software_installation)
        self.menu_cloud_control.add_separator()
        self.menu_cloud_control.add_command(label='Open a terminal', command=self.open_terminal)

        # link "menu_cloud_control" with "menu_bar"
        self.menu_bar.add_cascade(label='Cloud control', menu=self.menu_cloud_control)

        # create "menu_restriction_site_file" and add its menu items
        self.menu_restriction_site_file = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_restriction_site_file.add_command(label='Recreate data file', command=self.recreate_restriction_site_file)
        self.menu_restriction_site_file.add_command(label='Edit data file', command=self.edit_restriction_site_file)

        # create "menu_end_file" and add its menu items
        self.menu_end_file = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_end_file.add_command(label='Recreate data file', command=self.recreate_end_file)
        self.menu_end_file.add_command(label='Edit data file', command=self.edit_end_file)

        # create "menu_individual_file" and add its menu items
        self.menu_individual_file = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_individual_file.add_command(label='Recreate data file', command=self.recreate_individual_file)
        self.menu_individual_file.add_command(label='Edit data file', command=self.edit_individual_file)

        # create "menu_vcf_sample_file" and add its menu items
        self.menu_vcf_sample_file = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_vcf_sample_file.add_command(label='Recreate data file', command=self.recreate_vcf_sample_file)
        self.menu_vcf_sample_file.add_command(label='Edit data file', command=self.edit_vcf_sample_file)

        # create "menu_raddesinger_condition_file" and add its menu items
        self.menu_raddesinger_condition_file = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_raddesinger_condition_file.add_command(label='Recreate data file', command=self.recreate_raddesinger_condition_file)
        self.menu_raddesinger_condition_file.add_command(label='Edit data file', command=self.edit_raddesinger_condition_file)

        # create "menu_bowtie2" and add its menu items
        self.menu_bowtie2 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_bowtie2.add_command(label='Recreate config file', command=self.recreate_bowtie2_config_file)
        self.menu_bowtie2.add_command(label='Edit config file', command=self.edit_bowtie2_config_file)
        self.menu_bowtie2.add_separator()
        self.menu_bowtie2.add_command(label='Run read alignment process', command=self.run_bowtie2_process)

        # create "menu_busco" and add its menu items
        self.menu_busco = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_busco.add_command(label='Recreate config file', command=self.recreate_busco_config_file)
        self.menu_busco.add_command(label='Edit config file', command=self.edit_busco_config_file)
        self.menu_busco.add_separator()
        self.menu_busco.add_command(label='Run assembly quality assessment process', command=self.run_busco_process)

        # create "menu_cd_hit_est" and add its menu items
        self.menu_cd_hit_est = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cd_hit_est.add_command(label='Recreate config file', command=self.recreate_cd_hit_est_config_file)
        self.menu_cd_hit_est.add_command(label='Edit config file', command=self.edit_cd_hit_est_config_file)
        self.menu_cd_hit_est.add_separator()
        self.menu_cd_hit_est.add_command(label='Run transcriptome filtering process', command=self.run_cd_hit_est_process)

        # create "menu_cuffdiff" and add its menu items
        self.menu_cuffdiff = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cuffdiff.add_command(label='Recreate config file', command=self.recreate_cuffdiff_config_file)
        self.menu_cuffdiff.add_command(label='Edit config file', command=self.edit_cuffdiff_config_file)
        self.menu_cuffdiff.add_separator()
        self.menu_cuffdiff.add_command(label='Run differential expression process', command=self.run_cuffdiff_process)

        # create "menu_cufflinks_cuffmerge" and add its menu items
        self.menu_cufflinks_cuffmerge = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cufflinks_cuffmerge.add_command(label='Recreate config file', command=self.recreate_cufflinks_cuffmerge_config_file)
        self.menu_cufflinks_cuffmerge.add_command(label='Edit config file', command=self.edit_cufflinks_cuffmerge_config_file)
        self.menu_cufflinks_cuffmerge.add_separator()
        self.menu_cufflinks_cuffmerge.add_command(label='Run assembly process', command=self.run_cufflinks_cuffmerge_process)

        # create "menu_cuffquant" and add its menu items
        self.menu_cuffquant = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cuffquant.add_command(label='Recreate config file', command=self.recreate_cuffquant_config_file)
        self.menu_cuffquant.add_command(label='Edit config file', command=self.edit_cuffquant_config_file)
        self.menu_cuffquant.add_separator()
        self.menu_cuffquant.add_command(label='Run quantitation and differential expression process', command=self.run_cuffquant_process)

        # create "menu_cutadapt" and add its menu items
        self.menu_cutadapt = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_cutadapt.add_command(label='Recreate config file', command=self.recreate_cutadapt_config_file)
        self.menu_cutadapt.add_command(label='Edit config file', command=self.edit_cutadapt_config_file)
        self.menu_cutadapt.add_separator()
        self.menu_cutadapt.add_command(label='Run trimming process', command=self.run_cutadapt_process)

        # create "menu_ddradseq_simulation" and add its menu items
        self.menu_ddradseq_simulation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_ddradseq_simulation.add_command(label='Recreate config file', command=self.recreate_ddradseq_simulation_config_file)
        self.menu_ddradseq_simulation.add_command(label='Edit config file', command=self.edit_ddradseq_simulation_config_file)
        self.menu_ddradseq_simulation.add_separator()
        self.menu_ddradseq_simulation.add_command(label='Run simulation process', command=self.run_ddradseq_simulation_process)
        self.menu_ddradseq_simulation.add_command(label='Restart simulation process', command=self.restart_ddradseq_simulation_process)

        # create "menu_express" and add its menu items
        self.menu_express = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_express.add_command(label='Recreate config file', command=self.recreate_express_config_file)
        self.menu_express.add_command(label='Edit config file', command=self.edit_express_config_file)
        self.menu_express.add_separator()
        self.menu_express.add_command(label='Run quantitation process', command=self.run_express_process)

        # create "menu_fastqc" and add its menu items
        self.menu_fastqc = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_fastqc.add_command(label='Recreate config file', command=self.recreate_fastqc_config_file)
        self.menu_fastqc.add_command(label='Edit config file', command=self.edit_fastqc_config_file)
        self.menu_fastqc.add_separator()
        self.menu_fastqc.add_command(label='Run read quality process', command=self.run_fastqc_process)

        # create "menu_ggtrinity" and add its menu items
        self.menu_ggtrinity = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_ggtrinity.add_command(label='Recreate config file', command=self.recreate_ggtrinity_config_file)
        self.menu_ggtrinity.add_command(label='Edit config file', command=self.edit_ggtrinity_config_file)
        self.menu_ggtrinity.add_separator()
        self.menu_ggtrinity.add_command(label='Run assembly process', command=self.run_ggtrinity_process)
        self.menu_ggtrinity.add_command(label='Restart assembly process', command=self.restart_ggtrinity_process)

        # create "menu_gmap" and add its menu items
        self.menu_gmap = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_gmap.add_command(label='Recreate config file', command=self.recreate_gmap_config_file)
        self.menu_gmap.add_command(label='Edit config file', command=self.edit_gmap_config_file)
        self.menu_gmap.add_separator()
        self.menu_gmap.add_command(label='Run transcriptome alignment process', command=self.run_gmap_process)

        # create "menu_gsnap" and add its menu items
        self.menu_gsnap = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_gsnap.add_command(label='Recreate config file', command=self.recreate_gsnap_config_file)
        self.menu_gsnap.add_command(label='Edit config file', command=self.edit_gsnap_config_file)
        self.menu_gsnap.add_separator()
        self.menu_gsnap.add_command(label='Run read alignment process', command=self.run_gsnap_process)

        # create "menu_hisat2" and add its menu items
        self.menu_hisat2 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_hisat2.add_command(label='Recreate config file', command=self.recreate_hisat2_config_file)
        self.menu_hisat2.add_command(label='Edit config file', command=self.edit_hisat2_config_file)
        self.menu_hisat2.add_separator()
        self.menu_hisat2.add_command(label='Run read alignment process', command=self.run_hisat2_process)

        # create "menu_htseq_count" and add its menu items
        self.menu_htseq_count = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_htseq_count.add_command(label='Recreate config file', command=self.recreate_htseq_count_config_file)
        self.menu_htseq_count.add_command(label='Edit config file', command=self.edit_htseq_count_config_file)
        self.menu_htseq_count.add_separator()
        self.menu_htseq_count.add_command(label='Run quantitation process', command=self.run_htseq_count_process)

        # create "menu_insilico_read_normalization" and add its menu items
        self.menu_insilico_read_normalization = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_insilico_read_normalization.add_command(label='Recreate config file', command=self.recreate_insilico_read_normalization_config_file)
        self.menu_insilico_read_normalization.add_command(label='Edit config file', command=self.edit_insilico_read_normalization_config_file)
        self.menu_insilico_read_normalization.add_separator()
        self.menu_insilico_read_normalization.add_command(label='Run read normalization process', command=self.run_insilico_read_normalization_process)
        self.menu_insilico_read_normalization.add_command(label='Restart read normalization process', command=self.restart_insilico_read_normalization_process)

        # create "menu_ipyrad" and add its menu items
        self.menu_ipyrad = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_ipyrad.add_command(label='Recreate config file', command=self.recreate_ipyrad_config_file)
        self.menu_ipyrad.add_command(label='Edit config file', command=self.edit_ipyrad_config_file)
        self.menu_ipyrad.add_separator()
        self.menu_ipyrad.add_command(label='Run pipeline', command=self.run_ipyrad_process)

        # create "menu_kallisto" and add its menu items
        self.menu_kallisto = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_kallisto.add_command(label='Recreate config file', command=self.recreate_kallisto_config_file)
        self.menu_kallisto.add_command(label='Edit config file', command=self.edit_kallisto_config_file)
        self.menu_kallisto.add_separator()
        self.menu_kallisto.add_command(label='Run quantitation process', command=self.run_kallisto_process)

        # create "menu_quast" and add its menu items
        self.menu_quast = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_quast.add_command(label='Recreate config file', command=self.recreate_quast_config_file)
        self.menu_quast.add_command(label='Edit config file', command=self.edit_quast_config_file)
        self.menu_quast.add_separator()
        self.menu_quast.add_command(label='Run assembly quality assessment process', command=self.run_quast_process)

        # create "menu_rad_designer" and add its menu items
        self.menu_raddesigner = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_raddesigner.add_command(label='Recreate config file', command=self.recreate_raddesigner_config_file)
        self.menu_raddesigner.add_command(label='Edit config file', command=self.edit_raddesigner_config_file)
        self.menu_raddesigner.add_separator()
        self.menu_raddesigner.add_command(label='Run RAD designer process', command=self.run_raddesigner_process)
        self.menu_raddesigner.add_command(label='Restart RAD designer process', command=self.restart_raddesigner_process)

        # create "menu_ref_eval" and add its menu items
        # -- self.menu_ref_eval = tkinter.Menu(self.menu_bar, tearoff=0)
        # -- self.menu_ref_eval.add_command(label='Recreate config file', command=self.recreate_ref_eval_config_file)
        # -- self.menu_ref_eval.add_command(label='Edit config file', command=self.edit_ref_eval_config_file)
        # -- self.menu_ref_eval.add_separator()
        # -- self.menu_ref_eval.add_command(label='Run assembly quality assessment process', command=self.run_ref_eval_process)

        # create "menu_rnaquast" and add its menu items
        self.menu_rnaquast = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_rnaquast.add_command(label='Recreate config file', command=self.recreate_rnaquast_config_file)
        self.menu_rnaquast.add_command(label='Edit config file', command=self.edit_rnaquast_config_file)
        self.menu_rnaquast.add_separator()
        self.menu_rnaquast.add_command(label='Run assembly quality assessment process', command=self.run_rnaquast_process)

        # create "menu_rsem_eval" and add its menu items
        self.menu_rsem_eval = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_rsem_eval.add_command(label='Recreate config file', command=self.recreate_rsem_eval_config_file)
        self.menu_rsem_eval.add_command(label='Edit config file', command=self.edit_rsem_eval_config_file)
        self.menu_rsem_eval.add_separator()
        self.menu_rsem_eval.add_command(label='Run assembly quality assessment process', command=self.run_rsem_eval_process)

        # create "menu_rsitesearch" and add its menu items
        self.menu_rsitesearch = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_rsitesearch.add_command(label='Recreate config file', command=self.recreate_rsitesearch_config_file)
        self.menu_rsitesearch.add_command(label='Edit config file', command=self.edit_rsitesearch_config_file)
        self.menu_rsitesearch.add_separator()
        self.menu_rsitesearch.add_command(label='Run enzyme analysis process', command=self.run_rsitesearch_process)

        # create "menu_star" and add its menu items
        self.menu_star = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_star.add_command(label='Recreate config file', command=self.recreate_star_config_file)
        self.menu_star.add_command(label='Edit config file', command=self.edit_star_config_file)
        self.menu_star.add_separator()
        self.menu_star.add_command(label='Run read alignment process', command=self.run_star_process)

        # create "menu_starcode" and add its menu items
        self.menu_starcode = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_starcode.add_command(label='Recreate config file', command=self.recreate_starcode_config_file)
        self.menu_starcode.add_command(label='Edit config file', command=self.edit_starcode_config_file)
        self.menu_starcode.add_separator()
        self.menu_starcode.add_command(label='Run pseudo assembly process', command=self.run_starcode_process)

        # create "menu_soapdenovo2" and add its menu items
        self.menu_soapdenovo2 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_soapdenovo2.add_command(label='Recreate config file', command=self.recreate_soapdenovo2_config_file)
        self.menu_soapdenovo2.add_command(label='Edit config file', command=self.edit_soapdenovo2_config_file)
        self.menu_soapdenovo2.add_separator()
        self.menu_soapdenovo2.add_command(label='Run assembly process', command=self.run_soapdenovo2_process)
        self.menu_soapdenovo2.add_command(label='Restart pseudo assembly process', command=self.restart_soapdenovo2_process)

        # create "menu_soapdenovotrans" and add its menu items
        self.menu_soapdenovotrans = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_soapdenovotrans.add_command(label='Recreate config file', command=self.recreate_soapdenovotrans_config_file)
        self.menu_soapdenovotrans.add_command(label='Edit config file', command=self.edit_soapdenovotrans_config_file)
        self.menu_soapdenovotrans.add_separator()
        self.menu_soapdenovotrans.add_command(label='Run assembly process', command=self.run_soapdenovotrans_process)
        self.menu_soapdenovotrans.add_command(label='Restart assembly process', command=self.restart_soapdenovotrans_process)

        # create "menu_tophat" and add its menu items
        self.menu_tophat = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_tophat.add_command(label='Recreate config file', command=self.recreate_tophat_config_file)
        self.menu_tophat.add_command(label='Edit config file', command=self.edit_tophat_config_file)
        self.menu_tophat.add_separator()
        self.menu_tophat.add_command(label='Run read alignment process', command=self.run_tophat_process)

        # create "menu_transabyss" and add its menu items
        self.menu_transabyss = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_transabyss.add_command(label='Recreate config file', command=self.recreate_transabyss_config_file)
        self.menu_transabyss.add_command(label='Edit config file', command=self.edit_transabyss_config_file)
        self.menu_transabyss.add_separator()
        self.menu_transabyss.add_command(label='Run assembly process', command=self.run_transabyss_process)

        # create "menu_transcript_filter" and add its menu items
        self.menu_transcript_filter = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_transcript_filter.add_command(label='Recreate config file', command=self.recreate_transcript_filter_config_file)
        self.menu_transcript_filter.add_command(label='Edit config file', command=self.edit_transcript_filter_config_file)
        self.menu_transcript_filter.add_separator()
        self.menu_transcript_filter.add_command(label='Run transcriptome filtering process', command=self.run_transcript_filter_process)

        # create "menu_transcriptome_blastx" and add its menu items
        self.menu_transcriptome_blastx = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_transcriptome_blastx.add_command(label='Recreate config file', command=self.recreate_transcriptome_blastx_config_file)
        self.menu_transcriptome_blastx.add_command(label='Edit config file', command=self.edit_transcriptome_blastx_config_file)
        self.menu_transcriptome_blastx.add_separator()
        self.menu_transcriptome_blastx.add_command(label='Run annotation process', command=self.run_transcriptome_blastx_process)

        # create "menu_transrate" and add its menu items
        self.menu_transrate = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_transrate.add_command(label='Recreate config file', command=self.recreate_transrate_config_file)
        self.menu_transrate.add_command(label='Edit config file', command=self.edit_transrate_config_file)
        self.menu_transrate.add_separator()
        self.menu_transrate.add_command(label='Run assembly quality assessment process', command=self.run_transrate_process)

        # create "menu_trimmomatic" and add its menu items
        self.menu_trimmomatic = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_trimmomatic.add_command(label='Recreate config file', command=self.recreate_trimmomatic_config_file)
        self.menu_trimmomatic.add_command(label='Edit config file', command=self.edit_trimmomatic_config_file)
        self.menu_trimmomatic.add_separator()
        self.menu_trimmomatic.add_command(label='Run trimming process', command=self.run_trimmomatic_process)

        # create "menu_trinity" and add its menu items
        self.menu_trinity = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_trinity.add_command(label='Recreate config file', command=self.recreate_trinity_config_file)
        self.menu_trinity.add_command(label='Edit config file', command=self.edit_trinity_config_file)
        self.menu_trinity.add_separator()
        self.menu_trinity.add_command(label='Run assembly process', command=self.run_trinity_process)
        self.menu_trinity.add_command(label='Restart assembly process', command=self.restart_trinity_process)

        # create "menu_variant_calling" and add its menu items
        self.menu_variant_calling = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_variant_calling.add_command(label='Recreate config file', command=self.recreate_variant_calling_config_file)
        self.menu_variant_calling.add_command(label='Edit config file', command=self.edit_variant_calling_config_file)
        self.menu_variant_calling.add_separator()
        self.menu_variant_calling.add_command(label='Run variant calling process', command=self.run_variant_calling_process)
        self.menu_variant_calling.add_command(label='Restart variant calling process', command=self.restart_variant_calling_process)

        # create "menu_denovo_rnaseq_read_quality" and add its menu items
        self.menu_denovo_rnaseq_read_quality = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_read_quality.add_cascade(label=xlib.get_fastqc_name(), menu=self.menu_fastqc)

        # create "menu_denovo_rnaseq_trimming" and add its menu items
        self.menu_denovo_rnaseq_trimming = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_trimming.add_cascade(label=xlib.get_cutadapt_name(), menu=self.menu_cutadapt)
        self.menu_denovo_rnaseq_trimming.add_cascade(label=xlib.get_trimmomatic_name(), menu=self.menu_trimmomatic)

        # create "menu_denovo_rnaseq_digital_normalization" and add its menu items
        self.menu_denovo_rnaseq_digital_normalization = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_digital_normalization.add_cascade(label=f'{xlib.get_insilico_read_normalization_name()} ({xlib.get_trinity_name()} package)', menu=self.menu_insilico_read_normalization)

        # create "menu_denovo_rnaseq_assembly" and add its menu items
        self.menu_denovo_rnaseq_assembly = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_assembly.add_cascade(label=xlib.get_soapdenovotrans_name(), menu=self.menu_soapdenovotrans)
        self.menu_denovo_rnaseq_assembly.add_cascade(label=xlib.get_transabyss_name(), menu=self.menu_transabyss)
        self.menu_denovo_rnaseq_assembly.add_cascade(label=xlib.get_trinity_name(), menu=self.menu_trinity)

        # create "menu_denovo_rnaseq_read_alignment" add add its menu items
        self.menu_denovo_rnaseq_read_alignment = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_read_alignment.add_cascade(label=xlib.get_bowtie2_name(), menu=self.menu_bowtie2)

        # create "menu_denovo_rnaseq_transcriptome_quality_assessment" and add its menu items
        self.menu_denovo_rnaseq_transcriptome_quality_assessment = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=xlib.get_busco_name(), menu=self.menu_busco)
        self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=xlib.get_quast_name(), menu=self.menu_quast)
        self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=xlib.get_rnaquast_name(), menu=self.menu_rnaquast)
        self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=f'{xlib.get_rsem_eval_name()} ({xlib.get_detonate_name()} package)', menu=self.menu_rsem_eval)
        # -- self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=f'{xlib.get_ref_eval_name()} ({xlib.get_detonate_name()} package)', menu=self.menu_ref_eval)
        self.menu_denovo_rnaseq_transcriptome_quality_assessment.add_cascade(label=xlib.get_transrate_name(), menu=self.menu_transrate)

        # create "menu_denovo_rnaseq_transcriptome_filtering" and add its menu items
        self.menu_denovo_rnaseq_transcriptome_filtering = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_transcriptome_filtering.add_cascade(label=f'{xlib.get_cd_hit_est_name()} ({xlib.get_cd_hit_name()} package)', menu=self.menu_cd_hit_est)
        self.menu_denovo_rnaseq_transcriptome_filtering.add_cascade(label=f'{xlib.get_transcript_filter_name()} ({xlib.get_ngshelper_name()} package)', menu=self.menu_transcript_filter)

        # create "menu_denovo_rnaseq_quantitation" and add its menu items
        self.menu_denovo_rnaseq_quantitation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_quantitation.add_cascade(label=xlib.get_express_name(), menu=self.menu_express)
        self.menu_denovo_rnaseq_quantitation.add_cascade(label=xlib.get_kallisto_name(), menu=self.menu_kallisto)

        # create "menu_denovo_rnaseq_variant_calling" and add its menu items
        self.menu_denovo_rnaseq_variant_calling = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_variant_calling.add_cascade(label=xlib.get_variant_calling_name(), menu=self.menu_variant_calling)

        # create "menu_denovo_rnaseq_annotation" and add its menu items
        self.menu_denovo_rnaseq_annotation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq_annotation.add_cascade(label=f'{xlib.get_transcriptome_blastx_name()} ({xlib.get_ngshelper_name()} package)', menu=self.menu_transcriptome_blastx)

        # create "menu_denovo_rnaseq" add add its menu items
        self.menu_denovo_rnaseq = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_denovo_rnaseq.add_cascade(label='Read quality', menu=self.menu_denovo_rnaseq_read_quality)
        self.menu_denovo_rnaseq.add_cascade(label='Trimming', menu=self.menu_denovo_rnaseq_trimming)
        self.menu_denovo_rnaseq.add_cascade(label='Digital normalization', menu=self.menu_denovo_rnaseq_digital_normalization)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Assembly', menu=self.menu_denovo_rnaseq_assembly)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Read alignment', menu=self.menu_denovo_rnaseq_read_alignment)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Transcriptome quality assessment', menu=self.menu_denovo_rnaseq_transcriptome_quality_assessment)
        self.menu_denovo_rnaseq.add_cascade(label='Transcriptome filtering', menu=self.menu_denovo_rnaseq_transcriptome_filtering)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Quantitation', menu=self.menu_denovo_rnaseq_quantitation)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Variant calling', menu=self.menu_denovo_rnaseq_variant_calling)
        self.menu_denovo_rnaseq.add_separator()
        self.menu_denovo_rnaseq.add_cascade(label='Annotation', menu=self.menu_denovo_rnaseq_annotation)

        # link "menu_denovo_rnaseq" with "menu_bar"
        self.menu_bar.add_cascade(label='De novo RNA-seq', menu=self.menu_denovo_rnaseq)

        # create "menu_denovo_rnaseq_read_quality" and add its menu items
        self.menu_reference_based_rnaseq_read_quality = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_read_quality.add_cascade(label=xlib.get_fastqc_name(), menu=self.menu_fastqc)

        # create "menu_reference_based_rnaseq_trimming" and add its menu items
        self.menu_reference_based_rnaseq_trimming = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_trimming.add_cascade(label=xlib.get_cutadapt_name(), menu=self.menu_cutadapt)
        self.menu_reference_based_rnaseq_trimming.add_cascade(label=xlib.get_trimmomatic_name(), menu=self.menu_trimmomatic)

        # create "menu_reference_based_rnaseq_read_alignment" add add its menu items
        self.menu_reference_based_rnaseq_read_alignment = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_read_alignment.add_cascade(label=xlib.get_bowtie2_name(), menu=self.menu_bowtie2)
        self.menu_reference_based_rnaseq_read_alignment.add_cascade(label=f'{xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)', menu=self.menu_gsnap)
        self.menu_reference_based_rnaseq_read_alignment.add_cascade(label=xlib.get_hisat2_name(), menu=self.menu_hisat2)
        self.menu_reference_based_rnaseq_read_alignment.add_cascade(label=xlib.get_star_name(), menu=self.menu_star)
        self.menu_reference_based_rnaseq_read_alignment.add_cascade(label=xlib.get_tophat_name(), menu=self.menu_tophat)

        # create "menu_reference_based_rnaseq_assembly" add add its menu items
        self.menu_reference_based_rnaseq_assembly = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_assembly.add_cascade(label=f'{xlib.get_cufflinks_cuffmerge_name()} ({xlib.get_cufflinks_name()} package)', menu=self.menu_cufflinks_cuffmerge)
        self.menu_reference_based_rnaseq_assembly.add_cascade(label=f'{xlib.get_ggtrinity_name()} ({xlib.get_trinity_name()} package)', menu=self.menu_ggtrinity)

        # create "menu_reference_based_rnaseq_transcriptome_alignment" add add its menu items
        self.menu_reference_based_rnaseq_transcriptome_alignment = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_transcriptome_alignment.add_cascade(label=f'{xlib.get_gmap_name()} ({xlib.get_gmap_gsnap_name()} package)', menu=self.menu_gmap)

        # create "menu_reference_based_rnaseq_quantitation" and add its menu items
        self.menu_reference_based_rnaseq_quantitation = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_quantitation.add_cascade(label=f'{xlib.get_cuffquant_name()} ({xlib.get_cufflinks_name()} package)', menu=self.menu_cuffquant)
        self.menu_reference_based_rnaseq_quantitation.add_cascade(label=f'{xlib.get_htseq_count_name()} ({xlib.get_htseq_name()} package)', menu=self.menu_htseq_count)

        # create "menu_reference_based_rnaseq_differential_expression" and add its menu items
        self.menu_reference_based_rnaseq_differential_expression = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_differential_expression.add_cascade(label=f'{xlib.get_cuffdiff_name()} ({xlib.get_cufflinks_name()} package)', menu=self.menu_cuffdiff)

        # create "menu_reference_based_rnaseq_variant_calling" and add its menu items
        self.menu_reference_based_rnaseq_variant_calling = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq_variant_calling.add_cascade(label=xlib.get_variant_calling_name(), menu=self.menu_variant_calling)

        # create "menu_reference_based_rnaseq" add add its menu items
        self.menu_reference_based_rnaseq = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_based_rnaseq.add_cascade(label='Read quality', menu=self.menu_reference_based_rnaseq_read_quality)
        self.menu_reference_based_rnaseq.add_cascade(label='Trimming', menu=self.menu_reference_based_rnaseq_trimming)
        self.menu_reference_based_rnaseq.add_separator()
        self.menu_reference_based_rnaseq.add_cascade(label='Read alignment', menu=self.menu_reference_based_rnaseq_read_alignment)
        self.menu_reference_based_rnaseq.add_separator()
        self.menu_reference_based_rnaseq.add_cascade(label='Assembly', menu=self.menu_reference_based_rnaseq_assembly)
        self.menu_reference_based_rnaseq.add_separator()
        self.menu_reference_based_rnaseq.add_cascade(label='Transcriptome alignment', menu=self.menu_reference_based_rnaseq_transcriptome_alignment)
        self.menu_reference_based_rnaseq.add_separator()
        self.menu_reference_based_rnaseq.add_cascade(label='Quantitation', menu=self.menu_reference_based_rnaseq_quantitation)
        self.menu_reference_based_rnaseq.add_cascade(label='Differential expression', menu=self.menu_reference_based_rnaseq_differential_expression)
        self.menu_reference_based_rnaseq.add_separator()
        self.menu_reference_based_rnaseq.add_cascade(label='Variant calling', menu=self.menu_reference_based_rnaseq_variant_calling)

        # link "menu_reference_based_rnaseq" with "menu_bar"
        self.menu_bar.add_cascade(label='Reference-based RNA-seq', menu=self.menu_reference_based_rnaseq)

        # create "menu_radseq_data_files_maintenance" and add its menu items
        self.menu_radseq_data_files_maintenance = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_data_files_maintenance.add_cascade(label='Restriction site file', menu=self.menu_restriction_site_file)
        self.menu_radseq_data_files_maintenance.add_cascade(label='End file', menu=self.menu_end_file)
        self.menu_radseq_data_files_maintenance.add_cascade(label='Individual file', menu=self.menu_individual_file)
        self.menu_radseq_data_files_maintenance.add_cascade(label='VCF sample file', menu=self.menu_vcf_sample_file)
        self.menu_radseq_data_files_maintenance.add_cascade(label='RADdesigner condition file', menu=self.menu_raddesinger_condition_file)

        # create "menu_radseq_enzyme_analysis" and add its menu items
        self.menu_radseq_enzyme_analysis = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_enzyme_analysis.add_cascade(label=f'{xlib.get_rsitesearch_name()} ({xlib.get_ddradseqtools_name()} package)', menu=self.menu_rsitesearch)

        # create "menu_radseq_in_silico_simulations" and add its menu items
        self.menu_radseq_in_silico_simulations = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_in_silico_simulations.add_cascade(label=f'{xlib.get_ddradseq_simulation_name()} ({xlib.get_ddradseqtools_name()} package)', menu=self.menu_ddradseq_simulation)

        # create "menu_radseq_rad_design" and add its menu items
        self.menu_radseq_rad_design = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_rad_design.add_cascade(label='RADdesigner', menu=self.menu_raddesigner)

        # create "menu_radseq_read_quality" and add its menu items
        self.menu_radseq_read_quality = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_read_quality.add_cascade(label=xlib.get_fastqc_name(), menu=self.menu_fastqc)

        # create "menu_radseq_trimming" and add its menu items
        self.menu_radseq_trimming = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_trimming.add_cascade(label=xlib.get_cutadapt_name(), menu=self.menu_cutadapt)
        self.menu_radseq_trimming.add_cascade(label=xlib.get_trimmomatic_name(), menu=self.menu_trimmomatic)

        # create "menu_pseudo_assembly" and add its menu items
        self.menu_radseq_pseudo_assembly = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_pseudo_assembly.add_cascade(label=xlib.get_soapdenovo2_name(), menu=self.menu_soapdenovo2)
        self.menu_radseq_pseudo_assembly.add_cascade(label=xlib.get_starcode_name(), menu=self.menu_starcode)

        # create "menu_radseq_read_alignment" add add its menu items
        self.menu_radseq_read_alignment = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_read_alignment.add_cascade(label=xlib.get_bowtie2_name(), menu=self.menu_bowtie2)
        self.menu_radseq_read_alignment.add_cascade(label=f'{xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)', menu=self.menu_gsnap)
        self.menu_radseq_read_alignment.add_cascade(label=xlib.get_hisat2_name(), menu=self.menu_hisat2)

        # create "menu_radseq_variant_calling" and add its menu items
        self.menu_radseq_variant_calling = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_variant_calling.add_cascade(label=xlib.get_variant_calling_name(), menu=self.menu_variant_calling)

        # create "menu_radseq_pipelines" and add its menu items
        self.menu_radseq_pipelines = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq_pipelines.add_cascade(label=xlib.get_ipyrad_name(), menu=self.menu_ipyrad)

        # create "menu_radseq" add add its menu items
        self.menu_radseq = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_radseq.add_cascade(label='Maintenance of data files', menu=self.menu_radseq_data_files_maintenance)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Enzyme analysis', menu=self.menu_radseq_enzyme_analysis)
        self.menu_radseq.add_cascade(label='In silico simulations', menu=self.menu_radseq_in_silico_simulations)
        self.menu_radseq.add_cascade(label='RAD design', menu=self.menu_radseq_rad_design)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Read quality', menu=self.menu_radseq_read_quality)
        self.menu_radseq.add_cascade(label='Trimming', menu=self.menu_radseq_trimming)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Pseudo assembly', menu=self.menu_radseq_pseudo_assembly)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Read alignment', menu=self.menu_radseq_read_alignment)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Variant calling', menu=self.menu_radseq_variant_calling)
        self.menu_radseq.add_separator()
        self.menu_radseq.add_cascade(label='Pipelines', menu=self.menu_radseq_pipelines)

        # link "menu_radseq" with "menu_bar"
        self.menu_bar.add_cascade(label='RAD-seq', menu=self.menu_radseq)

        # create "menu_toa_configuration" and add its menu items
        self.menu_toa_configuration = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_configuration.add_command(label=f'Recreate {xlib.get_toa_name()} config file', command=self.recreate_toa_config_file)
        self.menu_toa_configuration.add_command(label=f'View {xlib.get_toa_name()} config file', command=self.view_toa_config_file)
        self.menu_toa_configuration.add_separator()
        self.menu_toa_configuration.add_command(label=f'Recreate {xlib.get_toa_name()} database', command=self.recreate_toa_database)
        # -- self.menu_toa_configuration.add_command(label=f'Rebuild {xlib.get_toa_name()} database', command=self.rebuild_toa_database)

        # create "menu_toa_basic_data" and add its menu items
        self.menu_toa_basic_data = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_basic_data.add_command(label='Recreate genomic dataset file', command=self.recreate_dataset_file)
        self.menu_toa_basic_data.add_command(label='Edit genomic dataset file', command=self.edit_dataset_file)
        self.menu_toa_basic_data.add_separator()
        self.menu_toa_basic_data.add_command(label='Recreate species file', command=self.recreate_species_file)
        self.menu_toa_basic_data.add_command(label='Edit species file', command=self.edit_species_file)
        self.menu_toa_basic_data.add_separator()
        self.menu_toa_basic_data.add_command(label='Download other basic data', command=self.download_basic_data)
        self.menu_toa_basic_data.add_separator()
        self.menu_toa_basic_data.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_basic_data)

        # create "menu_toa_gymno_01" and add its menu items
        self.menu_toa_gymno_01 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_gymno_01.add_command(label='Build proteome', command=self.build_gymno_01_proteome)
        self.menu_toa_gymno_01.add_separator()
        self.menu_toa_gymno_01.add_command(label='Download functional annotations from PLAZA server', command=self.download_gymno_01_data)
        self.menu_toa_gymno_01.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_gymno_01_data)

        # create "menu_toa_dicots_04" and add its menu items
        self.menu_toa_dicots_04 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_dicots_04.add_command(label='Build proteome', command=self.build_dicots_04_proteome)
        self.menu_toa_dicots_04.add_separator()
        self.menu_toa_dicots_04.add_command(label='Download functional annotations from PLAZA server', command=self.download_dicots_04_data)
        self.menu_toa_dicots_04.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_dicots_04_data)

        # create "menu_toa_monocots_04" and add its menu items
        self.menu_toa_monocots_04 = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_monocots_04.add_command(label='Build proteome', command=self.build_monocots_04_proteome)
        self.menu_toa_monocots_04.add_separator()
        self.menu_toa_monocots_04.add_command(label='Download functional annotations from PLAZA server', command=self.download_monocots_04_data)
        self.menu_toa_monocots_04.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_monocots_04_data)

        # create "menu_toa_refseq_plant" and add its menu items
        self.menu_toa_refseq_plant = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_refseq_plant.add_command(label='Build proteome', command=self.build_refseq_plant_proteome)

        # create "menu_toa_taxonomy" and add its menu items
        self.menu_toa_taxonomy = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_taxonomy.add_command(label='Download taxonomy data from NCBI server', command=self.download_taxonomy_data)

        # create "menu_toa_nt" and add its menu items
        self.menu_toa_nt = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_nt.add_command(label='Build database for BLAST+', command=self.build_blastplus_nt_db)

        # create "menu_toa_nucleotide_gi" and add its menu items
        self.menu_toa_nucleotide_gi = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_nucleotide_gi.add_command(label='Build identifier list using NCBI server', command=self.build_viridiplantae_nucleotide_gi_gilist)

        # create "menu_toa_nr" and add its menu items
        self.menu_toa_nr = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_nr.add_command(label='Build database for BLAST+', command=self.build_blastplus_nr_db)
        self.menu_toa_nr.add_command(label='Build database for DIAMOND', command=self.build_diamond_nr_db)

        # create "menu_toa_protein_gi" and add its menu items
        self.menu_toa_protein_gi = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_protein_gi.add_command(label='Build identifier list using NCBI server', command=self.build_viridiplantae_protein_gi_gilist)

        # create "menu_toa_gene" and add its menu items
        self.menu_toa_gene = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_gene.add_command(label='Download functional annotations from NCBI server', command=self.download_gene_data)
        self.menu_toa_gene.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_gene_data)

        # create "menu_toa_interpro" and add its menu items
        self.menu_toa_interpro = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_interpro.add_command(label='Download functional annotations from InterPro server', command=self.download_interpro_data)
        self.menu_toa_interpro.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_interpro_data)

        # create "menu_toa_go" and add its menu items
        self.menu_toa_go = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_go.add_command(label='Download functional annotations from Gene Ontology server', command=self.download_go_data)
        self.menu_toa_go.add_command(label=f'Load data into {xlib.get_toa_name()} database', command=self.load_go_data)

        # create "menu_toa_databases" and add its menu items
        self.menu_toa_databases = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_basic_data_name(), menu=self.menu_toa_basic_data)
        self.menu_toa_databases.add_separator()
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_gymno_01_name(), menu=self.menu_toa_gymno_01)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_dicots_04_name(), menu=self.menu_toa_dicots_04)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_monocots_04_name(), menu=self.menu_toa_monocots_04)
        self.menu_toa_databases.add_separator()
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_refseq_plant_name(), menu=self.menu_toa_refseq_plant)
        # -- self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_taxonomy_name(), menu=self.menu_toa_taxonomy)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_nt_name(), menu=self.menu_toa_nt)
        # -- self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_viridiplantae_nucleotide_gi_name(), menu=self.menu_toa_nucleotide_gi)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_nr_name(), menu=self.menu_toa_nr)
        # -- self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_viridiplantae_protein_gi_name(), menu=self.menu_toa_protein_gi)
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_gene_name(), menu=self.menu_toa_gene)
        self.menu_toa_databases.add_separator()
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_interpro_name(), menu=self.menu_toa_interpro)
        self.menu_toa_databases.add_separator()
        self.menu_toa_databases.add_cascade(label=xlib.get_toa_data_go_name(), menu=self.menu_toa_go)

        # create "menu_toa_nucleotide_pipeline" and add its menu items
        self.menu_toa_nucleotide_pipeline = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_nucleotide_pipeline.add_command(label='Recreate config file', command=self.recreate_nucleotide_pipeline_config_file)
        self.menu_toa_nucleotide_pipeline.add_command(label='Edit config file', command=self.edit_nucleotide_pipeline_config_file)
        self.menu_toa_nucleotide_pipeline.add_separator()
        self.menu_toa_nucleotide_pipeline.add_command(label='Run pipeline', command=self.run_nucleotide_pipeline_process)
        self.menu_toa_nucleotide_pipeline.add_command(label='Restart pipeline', command=self.restart_nucleotide_pipeline_process)

        # create "menu_toa_aminoacid_pipeline" and add its menu items
        self.menu_toa_aminoacid_pipeline = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_aminoacid_pipeline.add_command(label='Recreate config file', command=self.recreate_aminoacid_pipeline_config_file)
        self.menu_toa_aminoacid_pipeline.add_command(label='Edit config file', command=self.edit_aminoacid_pipeline_config_file)
        self.menu_toa_aminoacid_pipeline.add_separator()
        self.menu_toa_aminoacid_pipeline.add_command(label='Run pipeline', command=self.run_aminoacid_pipeline_process)
        self.menu_toa_aminoacid_pipeline.add_command(label='Restart pipeline', command=self.restart_aminoacid_pipeline_process)

        # create "menu_toa_annotation_merger" and add its menu items
        self.menu_toa_annotation_merger = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_annotation_merger.add_command(label='Recreate config file', command=self.recreate_annotation_merger_config_file)
        self.menu_toa_annotation_merger.add_command(label='Edit config file', command=self.edit_annotation_merger_config_file)
        self.menu_toa_annotation_merger.add_separator()
        self.menu_toa_annotation_merger.add_command(label='Run process', command=self.run_annotation_merger_process)

        # create "menu_toa_pipelines" and add its menu items
        self.menu_toa_pipelines = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_pipelines.add_cascade(label=f'{xlib.get_toa_name()} {xlib.get_toa_process_pipeline_nucleotide_name()}', menu=self.menu_toa_nucleotide_pipeline)
        self.menu_toa_pipelines.add_cascade(label=f'{xlib.get_toa_name()} {xlib.get_toa_process_pipeline_aminoacid_name()}', menu=self.menu_toa_aminoacid_pipeline)
        self.menu_toa_pipelines.add_separator()
        self.menu_toa_pipelines.add_cascade(label=f'Annotation merger of {xlib.get_toa_name()} pipelines', menu=self.menu_toa_annotation_merger)

        # create "menu_toa_alignment_stats" and add its menu items
        self.menu_toa_alignment_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_alignment_stats.add_command(label='# HITs per # HSPs data', command=self.view_hit_per_hsp_data)
        self.menu_toa_alignment_stats.add_command(label='# HITs per # HSPs plot', command=self.plot_hit_per_hsp_data)

        # create "menu_toa_annotation_dataset_stats" and add its menu items
        self.menu_toa_annotation_dataset_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_annotation_dataset_stats.add_command(label='Frequency distribution data', command=self.view_annotation_dataset_frequency)
        self.menu_toa_annotation_dataset_stats.add_command(label='Frequency distribution plot', command=self.plot_annotation_dataset_frequency)

        # create "menu_toa_species_stats" and add its menu items
        self.menu_toa_species_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_species_stats.add_command(label='Frequency distribution data', command=self.view_species_frequency)
        self.menu_toa_species_stats.add_command(label='Frequency distribution plot', command=self.plot_species_frequency)

        # create "menu_toa_family_stats" and add its menu items
        self.menu_toa_family_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_family_stats.add_command(label='Frequency distribution data', command=self.view_family_frequency)
        self.menu_toa_family_stats.add_command(label='Frequency distribution plot', command=self.plot_family_frequency)

        # create "menu_toa_phylum_stats" and add its menu items
        self.menu_toa_phylum_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_phylum_stats.add_command(label='Frequency distribution data', command=self.view_phylum_frequency)
        self.menu_toa_phylum_stats.add_command(label='Frequency distribution plot', command=self.plot_phylum_frequency)

        # create "menu_toa_ec_stats" and add its menu items
        self.menu_toa_ec_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_ec_stats.add_command(label='Frequency distribution data', command=self.view_ec_frequency)
        self.menu_toa_ec_stats.add_command(label='Frequency distribution plot', command=self.plot_ec_frequency)
        self.menu_toa_ec_stats.add_separator()
        self.menu_toa_ec_stats.add_command(label='# sequences per # ids data', command=self.view_seq_per_ec_data)
        self.menu_toa_ec_stats.add_command(label='# sequences per # ids plot', command=self.plot_seq_per_ec_data)

        # create "menu_toa_go_stats" and add its menu items
        self.menu_toa_go_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_go_stats.add_command(label='Frequency distribution per term data', command=self.view_go_frequency)
        self.menu_toa_go_stats.add_command(label='Frequency distribution per term plot', command=self.plot_go_frequency)
        self.menu_toa_go_stats.add_separator()
        self.menu_toa_go_stats.add_command(label='Frequency distribution per namespace data', command=self.view_namespace_frequency)
        self.menu_toa_go_stats.add_command(label='Frequency distribution per namespace plot', command=self.plot_namespace_frequency)
        self.menu_toa_go_stats.add_separator()
        self.menu_toa_go_stats.add_command(label='# sequences per # terms data', command=self.view_seq_per_go_data)
        self.menu_toa_go_stats.add_command(label='# sequences per # terms plot', command=self.plot_seq_per_go_data)

        # create "menu_toa_interpro_stats" and add its menu items
        self.menu_toa_interpro_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_interpro_stats.add_command(label='Frequency distribution data', command=self.view_interpro_frequency)
        self.menu_toa_interpro_stats.add_command(label='Frequency distribution plot', command=self.plot_interpro_frequency)
        self.menu_toa_interpro_stats.add_separator()
        self.menu_toa_interpro_stats.add_command(label='# sequences per # ids data', command=self.view_seq_per_interpro_data)
        self.menu_toa_interpro_stats.add_command(label='# sequences per # ids plot', command=self.plot_seq_per_interpro_data)

        # create "menu_toa_kegg_stats" and add its menu items
        self.menu_toa_kegg_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_kegg_stats.add_command(label='Frequency distribution data', command=self.view_kegg_frequency)
        self.menu_toa_kegg_stats.add_command(label='Frequency distribution plot', command=self.plot_kegg_frequency)
        self.menu_toa_kegg_stats.add_separator()
        self.menu_toa_kegg_stats.add_command(label='# sequences per # ids data', command=self.view_seq_per_kegg_data)
        self.menu_toa_kegg_stats.add_command(label='# sequences per # ids plot', command=self.plot_seq_per_kegg_data)

        # create "menu_toa_mapman_stats" and add its menu items
        self.menu_toa_mapman_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_mapman_stats.add_command(label='Frequency distribution data', command=self.view_mapman_frequency)
        self.menu_toa_mapman_stats.add_command(label='Frequency distribution plot', command=self.plot_mapman_frequency)
        self.menu_toa_mapman_stats.add_separator()
        self.menu_toa_mapman_stats.add_command(label='# sequences per # ids data', command=self.view_seq_per_mapman_data)
        self.menu_toa_mapman_stats.add_command(label='# sequences per # ids plot', command=self.plot_seq_per_mapman_data)

        # create "menu_toa_metacyc_stats" and add its menu items
        self.menu_toa_metacyc_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_metacyc_stats.add_command(label='Distribution data', command=self.view_metacyc_frequency)
        self.menu_toa_metacyc_stats.add_command(label='Distribution plot', command=self.plot_metacyc_frequency)
        self.menu_toa_metacyc_stats.add_separator()
        self.menu_toa_metacyc_stats.add_command(label='# sequences per # ids data', command=self.view_seq_per_metacyc_data)
        self.menu_toa_metacyc_stats.add_command(label='# sequences per # ids plot', command=self.plot_seq_per_metacyc_data)

        # create "menu_toa_stats" and add its menu items
        self.menu_toa_stats = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa_stats.add_cascade(label='Alignment', menu=self.menu_toa_alignment_stats)
        self.menu_toa_stats.add_separator()
        self.menu_toa_stats.add_cascade(label='Annotation datasets', menu=self.menu_toa_annotation_dataset_stats)
        self.menu_toa_stats.add_separator()
        self.menu_toa_stats.add_cascade(label='Species', menu=self.menu_toa_species_stats)
        self.menu_toa_stats.add_cascade(label='Family', menu=self.menu_toa_family_stats)
        self.menu_toa_stats.add_cascade(label='Phylum', menu=self.menu_toa_phylum_stats)
        self.menu_toa_stats.add_separator()
        self.menu_toa_stats.add_cascade(label='EC', menu=self.menu_toa_ec_stats)
        self.menu_toa_stats.add_cascade(label='Gene Ontology', menu=self.menu_toa_go_stats)
        self.menu_toa_stats.add_cascade(label='InterPro', menu=self.menu_toa_interpro_stats)
        self.menu_toa_stats.add_cascade(label='KEGG', menu=self.menu_toa_kegg_stats)
        self.menu_toa_stats.add_cascade(label='MapMan', menu=self.menu_toa_mapman_stats)
        self.menu_toa_stats.add_cascade(label='MetaCyc', menu=self.menu_toa_metacyc_stats)

        # create "menu_toa" add add its menu items
        self.menu_toa = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_toa.add_cascade(label='TOA configuration', menu=self.menu_toa_configuration)
        self.menu_toa.add_separator()
        self.menu_toa.add_cascade(label='Genomic databases', menu=self.menu_toa_databases)
        self.menu_toa.add_separator()
        self.menu_toa.add_cascade(label='Annotation pipelines', menu=self.menu_toa_pipelines)
        self.menu_toa.add_cascade(label='Statistics', menu=self.menu_toa_stats)

        # link "menu_toa" with "menu_bar"
        self.menu_bar.add_cascade(label='Taxonomy-oriented annotation', menu=self.menu_toa)

        # create "menu_reference_file_transfer" add add its menu items
        self.menu_reference_file_transfer = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_file_transfer.add_command(label='Recreate config file', command=self.recreate_reference_transfer_config_file)
        self.menu_reference_file_transfer.add_command(label='Edit config file', command=self.edit_reference_transfer_config_file)
        self.menu_reference_file_transfer.add_separator()
        self.menu_reference_file_transfer.add_command(label='Upload dataset to a cluster', command=self.upload_reference_dataset)

        # create "menu_reference_file_compression_decompression" add add its menu items
        self.menu_reference_file_compression_decompression = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_reference_file_compression_decompression.add_command(label='Recreate config file', command=self.recreate_reference_gzip_config_file)
        self.menu_reference_file_compression_decompression.add_command(label='Edit config file', command=self.edit_reference_gzip_config_file)
        self.menu_reference_file_compression_decompression.add_separator()
        self.menu_reference_file_compression_decompression.add_command(label='Run compression/decompression process', command=self.run_reference_gzip_process)

        # create "menu_database_file_transfer" add add its menu items
        self.menu_database_file_transfer = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_database_file_transfer.add_command(label='Recreate config file', command=self.recreate_database_transfer_config_file)
        self.menu_database_file_transfer.add_command(label='Edit config file', command=self.edit_database_transfer_config_file)
        self.menu_database_file_transfer.add_separator()
        self.menu_database_file_transfer.add_command(label='Upload dataset to a cluster', command=self.upload_database_dataset)

        # create "menu_database_file_compression_decompression" add add its menu items
        self.menu_database_file_compression_decompression = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_database_file_compression_decompression.add_command(label='Recreate config file', command=self.recreate_database_gzip_config_file)
        self.menu_database_file_compression_decompression.add_command(label='Edit config file', command=self.edit_database_gzip_config_file)
        self.menu_database_file_compression_decompression.add_separator()
        self.menu_database_file_compression_decompression.add_command(label='Run compression/decompression process', command=self.run_database_gzip_process)

        # create "menu_read_file_transfer" add add its menu items
        self.menu_read_file_transfer = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_read_file_transfer.add_command(label='Recreate config file', command=self.recreate_read_transfer_config_file)
        self.menu_read_file_transfer.add_command(label='Edit config file', command=self.edit_read_transfer_config_file)
        self.menu_read_file_transfer.add_separator()
        self.menu_read_file_transfer.add_command(label='Upload dataset to a cluster', command=self.upload_read_dataset)

        # create "menu_read_file_compression_decompression" add add its menu items
        self.menu_read_file_compression_decompression = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_read_file_compression_decompression.add_command(label='Recreate config file', command=self.recreate_read_gzip_config_file)
        self.menu_read_file_compression_decompression.add_command(label='Edit config file', command=self.edit_read_gzip_config_file)
        self.menu_read_file_compression_decompression.add_separator()
        self.menu_read_file_compression_decompression.add_command(label='Run compression/decompression process', command=self.run_read_gzip_process)

        # create "menu_result_file_transfer" add add its menu items
        self.menu_result_file_transfer = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_result_file_transfer.add_command(label='Recreate config file', command=self.recreate_result_transfer_config_file)
        self.menu_result_file_transfer.add_command(label='Edit config file', command=self.edit_result_transfer_config_file)
        self.menu_result_file_transfer.add_separator()
        self.menu_result_file_transfer.add_command(label='Download dataset from a cluster', command=self.download_result_dataset)

        # create "menu_result_file_compression_decompression" add add its menu items
        self.menu_result_file_compression_decompression = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_result_file_compression_decompression.add_command(label='Recreate config file', command=self.recreate_result_gzip_config_file)
        self.menu_result_file_compression_decompression.add_command(label='Edit config file', command=self.edit_result_gzip_config_file)
        self.menu_result_file_compression_decompression.add_separator()
        self.menu_result_file_compression_decompression.add_command(label='Run compression/decompression process', command=self.run_result_gzip_process)

        # create "menu_datasets" add add its menu items
        self.menu_datasets = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_datasets.add_command(label='List dataset', command=self.list_dataset)
        self.menu_datasets.add_separator()
        self.menu_datasets.add_cascade(label='Reference dataset file transfer', menu=self.menu_reference_file_transfer)
        self.menu_datasets.add_cascade(label='Reference dataset file compression/decompression', menu=self.menu_reference_file_compression_decompression)
        self.menu_datasets.add_command(label='Remove reference dataset', command=self.remove_reference_dataset)
        self.menu_datasets.add_separator()
        self.menu_datasets.add_cascade(label='Database file transfer', menu=self.menu_database_file_transfer)
        self.menu_datasets.add_cascade(label='Database file compression/decompression', menu=self.menu_database_file_compression_decompression)
        self.menu_datasets.add_command(label='Remove database', command=self.remove_database_dataset)
        self.menu_datasets.add_separator()
        self.menu_datasets.add_cascade(label='Read dataset file transfer', menu=self.menu_read_file_transfer)
        self.menu_datasets.add_cascade(label='Read dataset file compression/decompression', menu=self.menu_read_file_compression_decompression)
        self.menu_datasets.add_command(label='Remove read dataset', command=self.remove_read_dataset)
        self.menu_datasets.add_separator()
        self.menu_datasets.add_cascade(label='Result dataset file transfer', menu=self.menu_result_file_transfer)
        self.menu_datasets.add_cascade(label='Result dataset file compression/decompression', menu=self.menu_result_file_compression_decompression)
        self.menu_datasets.add_command(label='Remove result dataset', command=self.remove_result_dataset)
        self.menu_datasets.add_separator()
        self.menu_datasets.add_command(label='Remove experiment', command=self.remove_experiment)

        # link "menu_datasets" with "menu_bar"
        self.menu_bar.add_cascade(label='Datasets', menu=self.menu_datasets)

        # create "menu_logs" add add its menu items
        self.menu_logs = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_logs.add_command(label='View the cluster start log', command=self.view_cluster_start_log)
        self.menu_logs.add_separator()
        self.menu_logs.add_command(label='View submission logs in the local computer', command=self.view_submission_logs)
        self.menu_logs.add_separator()
        self.menu_logs.add_command(label='View result logs in the cluster', command=self.view_result_logs)

        # link "menu_logs" with "menu_bar"
        self.menu_bar.add_cascade(label='Logs', menu=self.menu_logs)

        # create "menu_help" add add its menu items
        self.menu_help = tkinter.Menu(self.menu_bar, tearoff=0)
        self.menu_help.add_command(label='View help', command=self.open_help, accelerator='F1')
        self.menu_help.add_separator()
        self.menu_help.add_command(label='About...', command=self.show_dialog_about)

        # link "menu_help" with "menu_bar"
        self.menu_bar.add_cascade(label='Help', menu=self.menu_help)

        #  assign "menu_bar" as the window menu
        self.root.config(menu=self.menu_bar)

        # create "frame_toolbar" and register it in "Main" with the grid geometry manager
        self.frame_toolbar = tkinter.Frame(self.root, borderwidth=1, relief='raised')
        self.frame_toolbar.grid(row=0, column=0, sticky='ew')

        # create and register "button_exit" in "frame_toolbar" with the pack geometry manager
        self.button_exit = tkinter.Button(self.frame_toolbar, command=self.exit, relief='flat', image=imagetk_exit)
        self.button_exit.image = imagetk_exit
        self.button_exit.pack(side='left', padx=2, pady=5)

        # create "frame_information" and register it in "Main" with the grid geometry manager
        self.frame_information = tkinter.Frame(self.root, borderwidth=1, relief='raised')
        self.frame_information.grid(row=1, column=0, sticky='ew')

        # create "label_environment_text" and register it in "frame_information" with the pack geometry manager
        self.label_environment_text = tkinter.Label(self.frame_information, text='Environment:')
        self.label_environment_text.pack(side='left', padx=(5,0))

        # create "label_environment_value" and register it in "frame_information" with the pack geometry manager
        self.label_environment_value = tkinter.Label(self.frame_information, text='', foreground='dark olive green')
        self.label_environment_value.pack(side='left', padx=(0,5))

        # create "label_region_text" and register it in "frame_information" with the pack geometry manager
        self.label_region_text = tkinter.Label(self.frame_information, text='Region:')
        self.label_region_text.pack(side='left', padx=(5,0))

        # create "label_region_value" and register it in "frame_information" with the pack geometry manager
        self.label_region_value = tkinter.Label(self.frame_information, text='', foreground='dark olive green')
        self.label_region_value.pack(side='left', padx=(0,0))

        # create "label_zone_text" and register it in "frame_information" with the pack geometry manager
        self.label_zone_text = tkinter.Label(self.frame_information, text='Zone:')
        self.label_zone_text.pack(side='left', padx=(5,0))

        # create "label_zone_value" and register it in "frame_information" with the pack geometry manager
        self.label_zone_value = tkinter.Label(self.frame_information, text='', foreground='dark olive green')
        self.label_zone_value.pack(side='left', padx=(0,0))

        # create "label_process" and register it in "frame_information" with the pack geometry manager
        self.label_process = tkinter.Label(self.frame_information, text='')
        self.label_process.pack(side='right', padx=(0,10))

        # create "container" and register it in "Main" with the grid geometry manager
        self.container = tkinter.Frame(self.root)
        self.container.grid(row=2, column=0, sticky='nsew')

        # disable certain menus on startup
        self.menu_bar.entryconfig('Cloud control', state='disabled')
        self.menu_bar.entryconfig('De novo RNA-seq', state='disabled')
        self.menu_bar.entryconfig('Reference-based RNA-seq', state='disabled')
        self.menu_bar.entryconfig('RAD-seq', state='disabled')
        self.menu_bar.entryconfig('Taxonomy-oriented annotation', state='disabled')
        self.menu_bar.entryconfig('Datasets', state='disabled')
        self.menu_bar.entryconfig('Logs', state='disabled')

        # link a handler to events
        self.root.bind('<F1>', self.open_help)
        self.root.bind('<Alt-F4>', self.exit)

        # link a handler to interactions between the application and the window manager
        self.root.protocol('WM_DELETE_WINDOW', self.exit)

    #---------------

    def set_environment(self):
        '''
        Set the configuration corresponding to an environment.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_set_environment" in "container" with the grid geometry manager
        form_set_environment = gcloud.FormSetEnvironment(self)
        form_set_environment.grid(row=0, column=0, sticky='nsew')

        # set "form_set_environment" as current form and add it in the forms dictionary
        self.current_form = 'form_set_environment'
        self.forms_dict[self.current_form] = form_set_environment

        # raise "form_set_environment" to front
        form_set_environment.tkraise()

    #---------------

    def recreate_ngscloud_config_file(self):
        '''
        Recreate the NGScloud config file corresponding to the environment.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_ngscloud_config_file" in "container" with the grid geometry manager
        form_recreate_ngscloud_config_file = gcloud.FormRecreateNGScloudConfigFile(self)
        form_recreate_ngscloud_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_ngscloud_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_create_ngscloud_config_file'
        self.forms_dict[self.current_form] = form_recreate_ngscloud_config_file

        # raise "form_recreate_ngscloud_config_file" to front
        form_recreate_ngscloud_config_file.tkraise()

    #---------------

    def view_ngscloud_config_file(self):
        '''
        List the NGScloud config file corresponding to the environment.
        '''

        # close the current form
        self.close_current_form()

        # get the NGScloud config file
        ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

        # create and show a instance DialogViewer to view the NGScloud config file
        dialog_viewer = gdialogs.DialogViewer(self.root, ngscloud_config_file)
        self.root.wait_window(dialog_viewer)

    #---------------

    def list_instance_types(self):
        '''
        List the characteristics of the instance types.
        '''

        # close the current form
        self.close_current_form()

        # set the head
        head = 'List instance types'

        # get the instance type dictionary
        instance_type_dict = xconfiguration.get_instance_type_dict(xconfiguration.get_cluster_mode_native())

        # build the data list
        data_list = ['use', 'id', 'vcpu', 'memory', 'processor', 'speed', 'nitro', 'starcluster', 'generation']

        # build the data dictionary
        data_dict = {}
        data_dict['use'] = {'text': 'Use', 'width': 135, 'alignment': 'left'}
        data_dict['id'] = {'text': 'Instance type', 'width': 120, 'alignment': 'left'}
        data_dict['vcpu'] = {'text': 'vCPUs', 'width': 50, 'alignment': 'right'}
        data_dict['memory'] = {'text': 'Memory (GiB)', 'width': 100, 'alignment': 'right'}
        data_dict['processor'] = {'text': 'Processor', 'width': 290, 'alignment': 'left'}
        data_dict['speed'] = {'text': 'Clock speed', 'width': 100, 'alignment': 'right'}
        data_dict['nitro'] = {'text': 'Nitro System', 'width': 90, 'alignment': 'left'}
        data_dict['starcluster'] = {'text': 'StarCluster', 'width': 110, 'alignment': 'left'}
        data_dict['generation'] = {'text': 'Generation', 'width': 90, 'alignment': 'left'}

        # create the dialog Table to show the instance types
        dialog_table = gdialogs.DialogTable(self.root, head, 400, 1100, data_list, data_dict, instance_type_dict, sorted(instance_type_dict.keys()))
        self.root.wait_window(dialog_table)

        # show warnings about characteristics and pricing
        message = 'You can consult the characteristics of the EC2 intance types in:\n\n'
        message += '   aws.amazon.com/ec2/instance-types/\n\n'
        message += 'and the EC2 pricing is detailed in:\n\n'
        message += '   aws.amazon.com/ec2/pricing/'
        tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def update_connection_data(self):
        '''
        Update the user id, access key id,  secret access key and contact e-mail address
        in the NGScloud config file corresponding to the environment.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_update_connection_data" in "container" with the grid geometry manager
        form_update_connection_data = gcloud.FormUpdateConnectionData(self)
        form_update_connection_data.grid(row=0, column=0, sticky='nsew')

        # set "form_create_update_connection_data" as current form and add it in the forms dictionary
        self.current_form = 'form_create_update_connection_data'
        self.forms_dict[self.current_form] = form_update_connection_data

        # raise "form_update_connection_data" to front
        form_update_connection_data.tkraise()

    #---------------

    def update_region_zone(self):
        '''
        Update the current region and zone names in the NGScloud config file
        corresponding to the environment.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_update_region_zone" in "container" with the grid geometry manager
        form_update_region_zone = gcloud.FormUpdateRegionZone(self)
        form_update_region_zone.grid(row=0, column=0, sticky='nsew')

        # set "form_update_region_zone" as current form and add it in the forms dictionary
        self.current_form = 'form_update_region_zone'
        self.forms_dict[self.current_form] = form_update_region_zone

        # raise "form_update_region_zone" to front
        form_update_region_zone.tkraise()

    #---------------

    def link_volumes(self):
        '''
        Link volumes in the NGScloud config file corresponding to the environment.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_link_volumes" in "container" with the grid geometry manager
        form_link_volumes = gcloud.FormLinkVolumes(self)
        form_link_volumes.grid(row=0, column=0, sticky='nsew')

        # set "form_link_volumes" as current form and add it in the forms dictionary
        self.current_form = 'form_link_volumes'
        self.forms_dict[self.current_form] = form_link_volumes

        # raise "form_link_volumes" to front
        form_link_volumes.tkraise()

    #---------------

    def list_keypairs(self):
        '''
        List the key pairs of a region.
        '''

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Security -  List key pairs'

        # get key pair dictionary
        keypairs_dict = xec2.get_keypair_dict(xconfiguration.get_current_region_name())

        # check if there are any key pairs createdview_ngscloud_config_file
        if keypairs_dict == {}:
            message = f'There is not any key pair created in region {xconfiguration.get_current_region_name()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {head}', message)
            return

        # build the data list
        data_list = ['keypair_name', 'fingerprint']

        # build the data dictionary
        data_dict = {}
        data_dict['keypair_name'] = {'text': 'Key Pair Name', 'width': 250, 'alignment': 'left'}
        data_dict['fingerprint'] = {'text': 'Fingerprint', 'width': 590, 'alignment': 'left'}

        # create the dialog Table to show the volumes created
        dialog_table = gdialogs.DialogTable(self.root, f'{head} in region {xconfiguration.get_current_region_name()}', 400, 900, data_list, data_dict, keypairs_dict, sorted(keypairs_dict.keys()))
        self.root.wait_window(dialog_table)

    #---------------

    def create_keypairs(self):
        '''
        Create all the key pairs of a region.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Security - Create key pairs'

        # confirm the creation of the key pairs
        message = f'The key pairs of the region {xconfiguration.get_current_region_name()} are going to be created.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # create key pairs
        if OK:
            (OK, error_list) = xec2.create_keypairs(xconfiguration.get_current_region_name())
            if OK:
                message = f'The key pairs and their corresponding local files of the region {xconfiguration.get_current_region_name()} file are created.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def list_clusters(self):
        '''
        List clusters.
        '''

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Cluster operation - List clusters'

        # list clusters
        dialog_log = gdialogs.DialogLog(self.root, head, xcluster.list_clusters.__name__)
        threading.Thread(target=self.root.wait_window, args=(dialog_log,)).start()
        threading.Thread(target=xcluster.list_clusters, args=(dialog_log, lambda: dialog_log.enable_button_close())).start()

    #---------------

    def create_cluster(self):
        '''
        Create a cluster.
        '''

        # close the current form
        self.close_current_form()

        # check if there is a cluster running
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if len(running_cluster_list) > 0:
            message = 'There is already a cluster running in this environment.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - Cluster operation - Create cluster', message)
            return

        # create and register "form_create_cluster" in "container" with the grid geometry manager
        form_create_cluster = gcloud.FormCreateCluster(self)
        form_create_cluster.grid(row=0, column=0, sticky='nsew')

        # set "form_create_cluster" as current form and add it in the forms dictionary
        self.current_form = 'form_create_cluster'
        self.forms_dict[self.current_form] = form_create_cluster

        # raise "form_create_cluster" to front
        form_create_cluster.tkraise()

    #---------------

    def terminate_cluster(self):
        '''
        Terminate a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_terminate_cluster" in "container" with the grid geometry manager
        form_terminate_cluster = gcloud.FormTerminateCluster(self, force=False)
        form_terminate_cluster.grid(row=0, column=0, sticky='nsew')

        # set "form_terminate_cluster" as current form and add it in the forms dictionary
        self.current_form = 'form_terminate_cluster'
        self.forms_dict[self.current_form] = form_terminate_cluster

        # raise "form_terminate_cluster" to front
        form_terminate_cluster.tkraise()

    #---------------

    def force_cluster_termination(self):
        '''
        Force the termination of a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_force_cluster_termination" in "container" with the grid geometry manager
        form_force_cluster_termination = gcloud.FormTerminateCluster(self, force=True)
        form_force_cluster_termination.grid(row=0, column=0, sticky='nsew')

        # set "form_force_cluster_termination" as current form and add it in the forms dictionary
        self.current_form = 'form_force_cluster_termination'
        self.forms_dict[self.current_form] = form_force_cluster_termination

        # raise "form_force_cluster_termination" to front
        form_force_cluster_termination.tkraise()

    #---------------

    def show_cluster_composition(self):
        '''
        Show cluster information of every node: OS, CPU number and memory.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_show_cluster_composition" in "container" with the grid geometry manager
        form_show_cluster_composition = gcloud.FormShowClusterComposition(self)
        form_show_cluster_composition.grid(row=0, column=0, sticky='nsew')

        # set "form_show_cluster_composition" as current form and add it in the forms dictionary
        self.current_form = 'form_show_cluster_composition'
        self.forms_dict[self.current_form] = form_show_cluster_composition

        # raise "form_show_cluster_composition" to front
        form_show_cluster_composition.tkraise()

    #---------------

    def show_status_batch_jobs(self):
        '''
        Show the status of batch jobs in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_show_status_batch_jobs" in "container" with the grid geometry manager
        form_show_status_batch_jobs = gcloud.FormShowStatusBatchJobs(self)
        form_show_status_batch_jobs.grid(row=0, column=0, sticky='nsew')

        # set "form_show_status_batch_jobs" as current form and add it in the forms dictionary
        self.current_form = 'form_show_status_batch_jobs'
        self.forms_dict[self.current_form] = form_show_status_batch_jobs

        # raise "form_show_status_batch_jobs" to front
        form_show_status_batch_jobs.tkraise()

    #---------------

    def kill_batch_job(self):
        '''
        Kill a batch job in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_kill_batch_job" in "container" with the grid geometry manager
        form_kill_batch_job = gcloud.FormKillBatchJob(self)
        form_kill_batch_job.grid(row=0, column=0, sticky='nsew')

        # set "form_kill_batch_job" as current form and add it in the forms dictionary
        self.current_form = 'form_kill_batch_job'
        self.forms_dict[self.current_form] = form_kill_batch_job

        # raise "form_kill_batch_job" to front
        form_kill_batch_job.tkraise()

    #---------------

    def list_nodes(self):
        '''
        List nodes running.
        '''

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Node operation - List nodes'

        # get the node dictionary
        node_dict = xec2.get_node_dict()

        # check if there are any nodes running
        if node_dict == {}:
            message = 'There is not any node running.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {head}', message)
            return

        # build the data list
        data_list = ['security_group_name', 'zone_name', 'node_name', 'node_id', 'state']

        # build the data dictionary
        data_dict = {}
        data_dict['security_group_name'] = {'text': 'Security Group Name', 'width': 200, 'alignment': 'left'}
        data_dict['zone_name'] = {'text': 'Zone', 'width': 100, 'alignment': 'left'}
        data_dict['node_name'] = {'text': 'Node Name', 'width': 200, 'alignment': 'left'}
        data_dict['node_id'] = {'text': 'Node Id', 'width': 190, 'alignment': 'left'}
        data_dict['state'] = {'text': 'State', 'width': 150, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        dialog_table = gdialogs.DialogTable(self.root, head, 400, 900, data_list, data_dict, node_dict, sorted(node_dict.keys()))
        self.root.wait_window(dialog_table)

    #---------------

    def add_node(self):
        '''
        Add a node in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_add_node" in "container" with the grid geometry manager
        form_add_node = gcloud.FormAddNode(self)
        form_add_node.grid(row=0, column=0, sticky='nsew')

        # set "form_add_node" as current form and add it in the forms dictionary
        self.current_form = 'form_add_node'
        self.forms_dict[self.current_form] = form_add_node

        # raise "form_add_node" to front
        form_add_node.tkraise()

    #---------------

    def remove_node(self):
        '''
        Remove a node in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_node" in "container" with the grid geometry manager
        form_remove_node = gcloud.FormRemoveNode(self)
        form_remove_node.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_node" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_node'
        self.forms_dict[self.current_form] = form_remove_node

        # raise "form_remove_node" to front
        form_remove_node.tkraise()

    #---------------

    def list_volumes(self):
        '''
        List volumes created.
        '''

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Volume operation - List volumes'

        # get the volume dictionary
        volume_dict = xec2.get_volume_dict()

        # check if there are any volumes created
        if volume_dict == {}:
            message = 'There is not any volume created.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {head}', message)
            return

        # build the data list
        data_list = ['zone_name', 'volume_name', 'volume_id', 'volume_type', 'size', 'state', 'attachments_number']

        # build the data dictionary
        data_dict = {}
        data_dict['zone_name'] = {'text': 'Zone', 'width': 100, 'alignment': 'left'}
        data_dict['volume_name'] = {'text': 'Volume Name', 'width': 210, 'alignment': 'left'}
        data_dict['volume_id'] = {'text': 'Volume Id', 'width': 210, 'alignment': 'left'}
        data_dict['volume_type'] = {'text': 'Volume Type', 'width': 105, 'alignment': 'left'}
        data_dict['size'] = {'text': 'Size (GiB)', 'width': 95, 'alignment': 'right'}
        data_dict['state'] = {'text': 'State', 'width': 95, 'alignment': 'left'}
        data_dict['attachments_number'] = {'text': 'Attachments', 'width': 95, 'alignment': 'right'}

        # create the dialog Table to show the volumes created
        dialog_table = gdialogs.DialogTable(self.root, head, 400, 930, data_list, data_dict, volume_dict, sorted(volume_dict.keys()))
        self.root.wait_window(dialog_table)

    #---------------

    def create_volume(self):
        '''
        Create a volume in the current zone.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_create_volume" in "container" with the grid geometry manager
        form_create_volume = gcloud.FormCreateVolume(self)
        form_create_volume.grid(row=0, column=0, sticky='nsew')

        # set "form_create_volume" as current form and add it in the forms dictionary
        self.current_form = 'form_create_volume'
        self.forms_dict[self.current_form] = form_create_volume

        # raise "form_create_volume" to front
        form_create_volume.tkraise()

    #---------------

    def remove_volume(self):
        '''
        Remove a volume in the current zone.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_volume" in "container" with the grid geometry manager
        form_remove_volume = gcloud.FormRemoveVolume(self)
        form_remove_volume.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_volume" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_volume'
        self.forms_dict[self.current_form] = form_remove_volume

        # raise "form_remove_volume" to front
        form_remove_volume.tkraise()

    #---------------

    def terminate_volume_creator(self):
        '''
        Terminate de volume creator of the current zone.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Volume operation - Terminate volume creator'

        # confirm the review of volumes links
        message = 'The volume creator is going to be terminated.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # terminate the volume creator
        if OK:
            dialog_log = gdialogs.DialogLog(self.root, head, xinstance.terminate_volume_creator.__name__)
            threading.Thread(target=self.root.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xinstance.terminate_volume_creator, args=(dialog_log, lambda: dialog_log.enable_button_close())).start()

    #---------------

    def mount_volume(self):
        '''
        Mount a volume in a node.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_mount_volume" in "container" with the grid geometry manager
        form_mount_volume = gcloud.FormMountVolume(self)
        form_mount_volume.grid(row=0, column=0, sticky='nsew')

        # set "form_mount_volume" as current form and add it in the forms dictionary
        self.current_form = 'form_mount_volume'
        self.forms_dict[self.current_form] = form_mount_volume

        # raise "form_mount_volume" to front
        form_mount_volume.tkraise()

    #---------------

    def unmount_volume(self):
        '''
        Unmount a volume in a node.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_unmount_volume" in "container" with the grid geometry manager
        form_unmount_volume = gcloud.FormUnmountVolume(self)
        form_unmount_volume.grid(row=0, column=0, sticky='nsew')

        # set "form_unmount_volume" as current form and add it in the forms dictionary
        self.current_form = 'form_unmount_volume'
        self.forms_dict[self.current_form] = form_unmount_volume

        # raise "form_unmount_volume" to front
        form_unmount_volume.tkraise()

    #---------------

    def install_miniconda3(self):
        '''
        Install the Miniconda3 in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_miniconda3" in "container" with the grid geometry manager
        form_install_miniconda3 = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_miniconda3_code())
        form_install_miniconda3.grid(row=0, column=0, sticky='nsew')

        # set "form_install_miniconda3" as current form and add it in the forms dictionary
        self.current_form = 'form_install_miniconda3'
        self.forms_dict[self.current_form] = form_install_miniconda3

        # raise "form_install_miniconda3" to front
        form_install_miniconda3.tkraise()

    #---------------

    def install_bcftools(self):
        '''
        Install the BCFtools in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_bcftools" in "container" with the grid geometry manager
        form_install_bcftools = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_bcftools_code())
        form_install_bcftools.grid(row=0, column=0, sticky='nsew')

        # set "form_install_bcftools" as current form and add it in the forms dictionary
        self.current_form = 'form_install_bcftools'
        self.forms_dict[self.current_form] = form_install_bcftools

        # raise "form_install_bcftools" to front
        form_install_bcftools.tkraise()

    #---------------

    def install_bedtools(self):
        '''
        Install the BEDtools in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_bedtools" in "container" with the grid geometry manager
        form_install_bedtools = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_bedtools_code())
        form_install_bedtools.grid(row=0, column=0, sticky='nsew')

        # set "form_install_bedtools" as current form and add it in the forms dictionary
        self.current_form = 'form_install_bedtools'
        self.forms_dict[self.current_form] = form_install_bedtools

        # raise "form_install_bedtools" to front
        form_install_bedtools.tkraise()

    #---------------

    def install_blastplus(self):
        '''
        Install the BLAST+ in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_blastplus" in "container" with the grid geometry manager
        form_install_blastplus = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_blastplus_code())
        form_install_blastplus.grid(row=0, column=0, sticky='nsew')

        # set "form_install_blastplus" as current form and add it in the forms dictionary
        self.current_form = 'form_install_blastplus'
        self.forms_dict[self.current_form] = form_install_blastplus

        # raise "form_install_blastplus" to front
        form_install_blastplus.tkraise()

    #---------------

    def install_bowtie2(self):
        '''
        Install the Bowtie 2 software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_bowtie2" in "container" with the grid geometry manager
        form_install_bowtie2 = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_bowtie2_code())
        form_install_bowtie2.grid(row=0, column=0, sticky='nsew')

        # set "form_install_bowtie2" as current form and add it in the forms dictionary
        self.current_form = 'form_install_bowtie2'
        self.forms_dict[self.current_form] = form_install_bowtie2

        # raise "form_install_bowtie2" to front
        form_install_bowtie2.tkraise()

    #---------------

    def install_busco(self):
        '''
        Install the BUSCO in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_busco" in "container" with the grid geometry manager
        form_install_busco = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_busco_code())
        form_install_busco.grid(row=0, column=0, sticky='nsew')

        # set "form_install_busco" as current form and add it in the forms dictionary
        self.current_form = 'form_install_busco'
        self.forms_dict[self.current_form] = form_install_busco

        # raise "form_install_busco" to front
        form_install_busco.tkraise()

    #---------------

    def install_cd_hit(self):
        '''
        Install the CD-HIT in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_cd_hit" in "container" with the grid geometry manager
        form_install_cd_hit = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_cd_hit_code())
        form_install_cd_hit.grid(row=0, column=0, sticky='nsew')

        # set "form_install_cd_hit" as current form and add it in the forms dictionary
        self.current_form = 'form_install_cd_hit'
        self.forms_dict[self.current_form] = form_install_cd_hit

        # raise "form_install_cd_hit" to front
        form_install_cd_hit.tkraise()

    #---------------

    def install_cufflinks(self):
        '''
        Install the Cufflinks in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_cufflinks" in "container" with the grid geometry manager
        form_install_cufflinks = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_cufflinks_code())
        form_install_cufflinks.grid(row=0, column=0, sticky='nsew')

        # set "form_install_cufflinks" as current form and add it in the forms dictionary
        self.current_form = 'form_install_cufflinks'
        self.forms_dict[self.current_form] = form_install_cufflinks

        # raise "form_install_cufflinks" to front
        form_install_cufflinks.tkraise()

    #---------------

    def install_cutadapt(self):
        '''
        Install the cutadapt in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_cutadapt" in "container" with the grid geometry manager
        form_install_cutadapt = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_cutadapt_code())
        form_install_cutadapt.grid(row=0, column=0, sticky='nsew')

        # set "form_install_cutadapt" as current form and add it in the forms dictionary
        self.current_form = 'form_install_cutadapt'
        self.forms_dict[self.current_form] = form_install_cutadapt

        # raise "form_install_cutadapt" to front
        form_install_cutadapt.tkraise()

    #---------------

    def install_ddradseqtools(self):
        '''
        Install the ddRADseqTools in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_ddradseqtools" in "container" with the grid geometry manager
        form_install_ddradseqtools = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_ddradseqtools_code())
        form_install_ddradseqtools.grid(row=0, column=0, sticky='nsew')

        # set "form_install_ddradseqtools" as current form and add it in the forms dictionary
        self.current_form = 'form_install_ddradseqtools'
        self.forms_dict[self.current_form] = form_install_ddradseqtools

        # raise "form_install_ddradseqtools" to front
        form_install_ddradseqtools.tkraise()

    #---------------

    def install_detonate(self):
        '''
        Install the DETONATE in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_detonate" in "container" with the grid geometry manager
        form_install_detonate = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_detonate_code())
        form_install_detonate.grid(row=0, column=0, sticky='nsew')

        # set "form_install_detonate" as current form and add it in the forms dictionary
        self.current_form = 'form_install_detonate'
        self.forms_dict[self.current_form] = form_install_detonate

        # raise "form_install_detonate" to front
        form_install_detonate.tkraise()

    #---------------

    def install_diamond(self):
        '''
        Install the DIAMOND software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_diamond" in "container" with the grid geometry manager
        form_install_diamond = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_diamond_code())
        form_install_diamond.grid(row=0, column=0, sticky='nsew')

        # set "form_install_diamond" as current form and add it in the forms dictionary
        self.current_form = 'form_install_diamond'
        self.forms_dict[self.current_form] = form_install_diamond

        # raise "form_install_diamond" to front
        form_install_diamond.tkraise()

    #---------------

    def install_entrez_direct(self):
        '''
        Install the Entrez Direct software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_entrez_direct" in "container" with the grid geometry manager
        form_install_entrez_direct = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_entrez_direct_code())
        form_install_entrez_direct.grid(row=0, column=0, sticky='nsew')

        # set "form_install_entrez_direct" as current form and add it in the forms dictionary
        self.current_form = 'form_install_entrez_direct'
        self.forms_dict[self.current_form] = form_install_entrez_direct

        # raise "form_install_entrez_direct" to front
        form_install_entrez_direct.tkraise()

    #---------------

    def install_express(self):
        '''
        Install the eXpress software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_express" in "container" with the grid geometry manager
        form_install_express = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_express_code())
        form_install_express.grid(row=0, column=0, sticky='nsew')

        # set "form_install_express" as current form and add it in the forms dictionary
        self.current_form = 'form_install_express'
        self.forms_dict[self.current_form] = form_install_express

        # raise "form_install_express" to front
        form_install_express.tkraise()

    #---------------

    def install_fastqc(self):
        '''
        Install the FastQC software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_fastqc" in "container" with the grid geometry manager
        form_install_fastqc = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_fastqc_code())
        form_install_fastqc.grid(row=0, column=0, sticky='nsew')

        # set "form_install_fastqc" as current form and add it in the forms dictionary
        self.current_form = 'form_install_fastqc'
        self.forms_dict[self.current_form] = form_install_fastqc

        # raise "form_install_fastqc" to front
        form_install_fastqc.tkraise()

    #---------------

    def install_gmap_gsnap(self):
        '''
        Install the GMAP-GSNAP in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_gmap_gsnap" in "container" with the grid geometry manager
        form_install_gmap_gsnap = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_gmap_gsnap_code())
        form_install_gmap_gsnap.grid(row=0, column=0, sticky='nsew')

        # set "form_install_gmap_gsnap" as current form and add it in the forms dictionary
        self.current_form = 'form_install_gmap_gsnap'
        self.forms_dict[self.current_form] = form_install_gmap_gsnap

        # raise "form_install_gmap_gsnap" to front
        form_install_gmap_gsnap.tkraise()

    #---------------

    def install_hisat2(self):
        '''
        Install the HISAT2 software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_hisat2" in "container" with the grid geometry manager
        form_install_hisat2 = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_hisat2_code())
        form_install_hisat2.grid(row=0, column=0, sticky='nsew')

        # set "form_install_hisat2" as current form and add it in the forms dictionary
        self.current_form = 'form_install_hisat2'
        self.forms_dict[self.current_form] = form_install_hisat2

        # raise "form_install_hisat2" to front
        form_install_hisat2.tkraise()

    #---------------

    def install_htseq(self):
        '''
        Install the HTSeq software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_htseq" in "container" with the grid geometry manager
        form_install_htseq = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_htseq_code())
        form_install_htseq.grid(row=0, column=0, sticky='nsew')

        # set "form_install_htseq" as current form and add it in the forms dictionary
        self.current_form = 'form_install_htseq'
        self.forms_dict[self.current_form] = form_install_htseq

        # raise "form_install_htseq" to front
        form_install_htseq.tkraise()

    #---------------

    def install_ipyrad(self):
        '''
        Install the ipyrad in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_ipyrad" in "container" with the grid geometry manager
        form_install_ipyrad = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_ipyrad_code())
        form_install_ipyrad.grid(row=0, column=0, sticky='nsew')

        # set "form_install_ipyrad" as current form and add it in the forms dictionary
        self.current_form = 'form_install_ipyrad'
        self.forms_dict[self.current_form] = form_install_ipyrad

        # raise "form_install_ipyrad" to front
        form_install_ipyrad.tkraise()

    #---------------

    def install_kallisto(self):
        '''
        Install the kallisto in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_kallisto" in "container" with the grid geometry manager
        form_install_kallisto = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_kallisto_code())
        form_install_kallisto.grid(row=0, column=0, sticky='nsew')

        # set "form_install_kallisto" as current form and add it in the forms dictionary
        self.current_form = 'form_install_kallisto'
        self.forms_dict[self.current_form] = form_install_kallisto

        # raise "form_install_kallisto" to front
        form_install_kallisto.tkraise()

    #---------------

    def install_ngshelper(self):
        '''
        Install the NGShelper in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_ngshelper" in "container" with the grid geometry manager
        form_install_ngshelper = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_ngshelper_code())
        form_install_ngshelper.grid(row=0, column=0, sticky='nsew')

        # set "form_install_ngshelper" as current form and add it in the forms dictionary
        self.current_form = 'form_install_ngshelper'
        self.forms_dict[self.current_form] = form_install_ngshelper

        # raise "form_install_ngshelper" to front
        form_install_ngshelper.tkraise()

    #---------------

    def install_quast(self):
        '''
        Install the QUAST in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_quast" in "container" with the grid geometry manager
        form_install_quast = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_quast_code())
        form_install_quast.grid(row=0, column=0, sticky='nsew')

        # set "form_install_quast" as current form and add it in the forms dictionary
        self.current_form = 'form_install_quast'
        self.forms_dict[self.current_form] = form_install_quast

        # raise "form_install_quast" to front
        form_install_quast.tkraise()

    #---------------

    def install_raddesigner(self):
        '''
        Install the RADdesigner in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_raddesigner" in "container" with the grid geometry manager
        form_install_raddesigner = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_raddesigner_code())
        form_install_raddesigner.grid(row=0, column=0, sticky='nsew')

        # set "form_install_raddesigner" as current form and add it in the forms dictionary
        self.current_form = 'form_install_raddesigner'
        self.forms_dict[self.current_form] = form_install_raddesigner

        # raise "form_install_raddesigner" to front
        form_install_raddesigner.tkraise()

    #---------------

    def install_rnaquast(self):
        '''
        Install the rnaQUAST in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_rnaquast" in "container" with the grid geometry manager
        form_install_rnaquast = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_rnaquast_code())
        form_install_rnaquast.grid(row=0, column=0, sticky='nsew')

        # set "form_install_rnaquast" as current form and add it in the forms dictionary
        self.current_form = 'form_install_rnaquast'
        self.forms_dict[self.current_form] = form_install_rnaquast

        # raise "form_install_rnaquast" to front
        form_install_rnaquast.tkraise()

    #---------------

    def install_rsem(self):
        '''
        Install the RSEM in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_rsem" in "container" with the grid geometry manager
        form_install_rsem = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_rsem_code())
        form_install_rsem.grid(row=0, column=0, sticky='nsew')

        # set "form_install_rsem" as current form and add it in the forms dictionary
        self.current_form = 'form_install_rsem'
        self.forms_dict[self.current_form] = form_install_rsem

        # raise "form_install_rsem" to front
        form_install_rsem.tkraise()

    #---------------

    def install_samtools(self):
        '''
        Install the SAMtools in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_samtools" in "container" with the grid geometry manager
        form_install_samtools = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_samtools_code())
        form_install_samtools.grid(row=0, column=0, sticky='nsew')

        # set "form_install_samtools" as current form and add it in the forms dictionary
        self.current_form = 'form_install_samtools'
        self.forms_dict[self.current_form] = form_install_samtools

        # raise "form_install_samtools" to front
        form_install_samtools.tkraise()

    #---------------

    def install_soapdenovo2(self):
        '''
        Install the SOAPdenovo2 software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_soapdenovo2" in "container" with the grid geometry manager
        form_install_soapdenovo2 = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_soapdenovo2_code())
        form_install_soapdenovo2.grid(row=0, column=0, sticky='nsew')

        # set "form_install_soapdenovo2" as current form and add it in the forms dictionary
        self.current_form = 'form_install_soapdenovo2'
        self.forms_dict[self.current_form] = form_install_soapdenovo2

        # raise "form_install_soapdenovo2" to front
        form_install_soapdenovo2.tkraise()

    #---------------

    def install_soapdenovotrans(self):
        '''
        Install the SOAPdenovo-Trans software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_soapdenovotrans" in "container" with the grid geometry manager
        form_install_soapdenovotrans = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_soapdenovotrans_code())
        form_install_soapdenovotrans.grid(row=0, column=0, sticky='nsew')

        # set "form_install_soapdenovotrans" as current form and add it in the forms dictionary
        self.current_form = 'form_install_soapdenovotrans'
        self.forms_dict[self.current_form] = form_install_soapdenovotrans

        # raise "form_install_soapdenovotrans" to front
        form_install_soapdenovotrans.tkraise()

    #---------------

    def install_star(self):
        '''
        Install the STAR in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_star" in "container" with the grid geometry manager
        form_install_star = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_star_code())
        form_install_star.grid(row=0, column=0, sticky='nsew')

        # set "form_install_star" as current form and add it in the forms dictionary
        self.current_form = 'form_install_star'
        self.forms_dict[self.current_form] = form_install_star

        # raise "form_install_star" to front
        form_install_star.tkraise()

    #---------------

    def install_starcode(self):
        '''
        Install the starcode in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_starcode" in "container" with the grid geometry manager
        form_install_starcode = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_starcode_code())
        form_install_starcode.grid(row=0, column=0, sticky='nsew')

        # set "form_install_starcode" as current form and add it in the forms dictionary
        self.current_form = 'form_install_starcode'
        self.forms_dict[self.current_form] = form_install_starcode

        # raise "form_install_starcode" to front
        form_install_starcode.tkraise()

    #---------------

    def install_toa(self):
        '''
        Install the TOA in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_toa" in "container" with the grid geometry manager
        form_install_toa = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_toa_code())
        form_install_toa.grid(row=0, column=0, sticky='nsew')

        # set "form_install_toa" as current form and add it in the forms dictionary
        self.current_form = 'form_install_toa'
        self.forms_dict[self.current_form] = form_install_toa

        # raise "form_install_toa" to front
        form_install_toa.tkraise()

    #---------------

    def install_tophat(self):
        '''
        Install the TopHat in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_tophat" in "container" with the grid geometry manager
        form_install_tophat = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_tophat_code())
        form_install_tophat.grid(row=0, column=0, sticky='nsew')

        # set "form_install_tophat" as current form and add it in the forms dictionary
        self.current_form = 'form_install_tophat'
        self.forms_dict[self.current_form] = form_install_tophat

        # raise "form_install_tophat" to front
        form_install_tophat.tkraise()

    #---------------

    def install_transabyss(self):
        '''
        Install the Trans-ABySS in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_transabyss" in "container" with the grid geometry manager
        form_install_transabyss = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_transabyss_code())
        form_install_transabyss.grid(row=0, column=0, sticky='nsew')

        # set "form_install_transabyss" as current form and add it in the forms dictionary
        self.current_form = 'form_install_transabyss'
        self.forms_dict[self.current_form] = form_install_transabyss

        # raise "form_install_transabyss" to front
        form_install_transabyss.tkraise()

    #---------------

    def install_transdecoder(self):
        '''
        Install the TransDecoder software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_transdecoder" in "container" with the grid geometry manager
        form_install_transdecoder = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_transdecoder_code())
        form_install_transdecoder.grid(row=0, column=0, sticky='nsew')

        # set "form_install_transdecoder" as current form and add it in the forms dictionary
        self.current_form = 'form_install_transdecoder'
        self.forms_dict[self.current_form] = form_install_transdecoder

        # raise "form_install_transdecoder" to front
        form_install_transdecoder.tkraise()

    #---------------

    def install_transrate(self):
        '''
        Install the Transrate in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_transrate" in "container" with the grid geometry manager
        form_install_transrate = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_transrate_code())
        form_install_transrate.grid(row=0, column=0, sticky='nsew')

        # set "form_install_transrate" as current form and add it in the forms dictionary
        self.current_form = 'form_install_transrate'
        self.forms_dict[self.current_form] = form_install_transrate

        # raise "form_install_transrate" to front
        form_install_transrate.tkraise()

    #---------------

    def install_trimmomatic(self):
        '''
        Install the Trimmomatic software in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_trimmomatic" in "container" with the grid geometry manager
        form_install_trimmomatic = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_trimmomatic_code())
        form_install_trimmomatic.grid(row=0, column=0, sticky='nsew')

        # set "form_install_trimmomatic" as current form and add it in the forms dictionary
        self.current_form = 'form_install_trimmomatic'
        self.forms_dict[self.current_form] = form_install_trimmomatic

        # raise "form_install_trimmomatic" to front
        form_install_trimmomatic.tkraise()

    #---------------

    def install_trinity(self):
        '''
        Install the Trinity in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_trinity" in "container" with the grid geometry manager
        form_install_trinity = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_trinity_code())
        form_install_trinity.grid(row=0, column=0, sticky='nsew')

        # set "form_install_trinity" as current form and add it in the forms dictionary
        self.current_form = 'form_install_trinity'
        self.forms_dict[self.current_form] = form_install_trinity

        # raise "form_install_trinity" to front
        form_install_trinity.tkraise()

    #---------------

    def install_vcftools(self):
        '''
        Install the VCFtools in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_vcftools" in "container" with the grid geometry manager
        form_install_vcftools = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_vcftools_code())
        form_install_vcftools.grid(row=0, column=0, sticky='nsew')

        # set "form_install_vcftools" as current form and add it in the forms dictionary
        self.current_form = 'form_install_vcftools'
        self.forms_dict[self.current_form] = form_install_vcftools

        # raise "form_install_vcftools" to front
        form_install_vcftools.tkraise()

    #---------------

    def install_vcftools_perl_libraries(self):
        '''
        Install the VCFtools Perl libraries in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_vcftools_perl_libraries" in "container" with the grid geometry manager
        form_install_vcftools_perl_libraries = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_vcftools_perl_libraries_code())
        form_install_vcftools_perl_libraries.grid(row=0, column=0, sticky='nsew')

        # set "form_install_vcftools_perl_libraries" as current form and add it in the forms dictionary
        self.current_form = 'form_install_vcftools_perl_libraries'
        self.forms_dict[self.current_form] = form_install_vcftools_perl_libraries

        # raise "form_install_vcftools_perl_libraries" to front
        form_install_vcftools_perl_libraries.tkraise()

    #---------------

    def install_vsearch(self):
        '''
        Install the VSEARCH in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_vsearch" in "container" with the grid geometry manager
        form_install_vsearch = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_vsearch_code())
        form_install_vsearch.grid(row=0, column=0, sticky='nsew')

        # set "form_install_vsearch" as current form and add it in the forms dictionary
        self.current_form = 'form_install_vsearch'
        self.forms_dict[self.current_form] = form_install_vsearch

        # raise "form_install_vsearch" to front
        form_install_vsearch.tkraise()

    #---------------

    def install_r(self):
        '''
        Install the R in the cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_install_r" in "container" with the grid geometry manager
        form_install_r = gbioinfoapp.FormInstallBioinfoApp(self, app=xlib.get_r_code())
        form_install_r.grid(row=0, column=0, sticky='nsew')

        # set "form_install_r" as current form and add it in the forms dictionary
        self.current_form = 'form_install_r'
        self.forms_dict[self.current_form] = form_install_r

        # raise "form_install_r" to front
        form_install_r.tkraise()

    #---------------

    def open_terminal(self):
        '''
        Open a terminal windows of a node cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_open_terminal" in "container" with the grid geometry manager
        form_open_terminal = gcloud.FormOpenTerminal(self)
        form_open_terminal.grid(row=0, column=0, sticky='nsew')

        # set "form_open_terminal" as current form and add it in the forms dictionary
        self.current_form = 'form_open_terminal'
        self.forms_dict[self.current_form] = form_open_terminal

        # raise "form_open_terminal" to front
        form_open_terminal.tkraise()

    #---------------
    # Bowtie2
    #---------------

    def recreate_bowtie2_config_file(self):
        '''
        Recreate the Bowtie2 config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_bowtie2_config_file" in "container" with the grid geometry manager
        form_recreate_bowtie2_config_file = gbioinfoapp.FormRecreateBowtie2ConfigFile(self)
        form_recreate_bowtie2_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_bowtie2_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_bowtie2_config_file'
        self.forms_dict[self.current_form] = form_recreate_bowtie2_config_file

        # raise "form_recreate_bowtie2_config_file" to front
        form_recreate_bowtie2_config_file.tkraise()

    #---------------

    def edit_bowtie2_config_file(self):
        '''
        Edit the Bowtie2 config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_bowtie2_name()} - Edit config file'

        # edit the Bowtie2 config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xbowtie2.get_bowtie2_config_file())
        self.root.wait_window(dialog_editor)

        # check the Bowtie2 config file
        (OK, error_list) = xbowtie2.check_bowtie2_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_bowtie2_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_bowtie2_process(self):
        '''
        Run a Bowtie2 process corresponding to the options in Bowtie2 config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_bowtie2_process" in "container" with the grid geometry manager
        form_run_bowtie2_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_bowtie2_code())
        form_run_bowtie2_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_bowtie2_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_bowtie2_process'
        self.forms_dict[self.current_form] = form_run_bowtie2_process

        # raise "form_run_bowtie2_process" to front
        form_run_bowtie2_process.tkraise()

    #---------------
    # BUSCO
    #---------------

    def recreate_busco_config_file(self):
        '''
        Recreate the BUSCO config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_busco_config_file" in "container" with the grid geometry manager
        form_recreate_busco_config_file = gbioinfoapp.FormRecreateBuscoConfigFile(self)
        form_recreate_busco_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_busco_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_busco_config_file'
        self.forms_dict[self.current_form] = form_recreate_busco_config_file

        # raise "form_recreate_busco_config_file" to front
        form_recreate_busco_config_file.tkraise()

    #---------------

    def edit_busco_config_file(self):
        '''
        Edit the BUSCO config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_busco_name()} - Edit config file'

        # edit the BUSCO config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xbusco.get_busco_config_file())
        self.root.wait_window(dialog_editor)

        # check the BUSCO config file
        (OK, error_list) = xbusco.check_busco_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_busco_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_busco_process(self):
        '''
        Run a BUSCO process corresponding to the options in BUSCO config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_busco_process" in "container" with the grid geometry manager
        form_run_busco_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_busco_code())
        form_run_busco_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_busco_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_busco_process'
        self.forms_dict[self.current_form] = form_run_busco_process

        # raise "form_run_busco_process" to front
        form_run_busco_process.tkraise()

    #---------------
    # CD-HIT-EST
    #---------------

    def recreate_cd_hit_est_config_file(self):
        '''
        Create the CD-HIT-EST config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_cd_hit_est_config_file" in "container" with the grid geometry manager
        form_cd_hit_est_eval_config_file = gbioinfoapp.FormRecreateCdHitEstConfigFile(self)
        form_cd_hit_est_eval_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_cd_hit_est_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_cd_hit_est_config_file'
        self.forms_dict[self.current_form] = form_cd_hit_est_eval_config_file

        # raise "form_recreate_cd_hit_est_config_file" to front
        form_cd_hit_est_eval_config_file.tkraise()

    #---------------

    def edit_cd_hit_est_config_file(self):
        '''
        Edit the CD-HIT-EST config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_cd_hit_est_name()} - Edit config file'

        # edit the CD-HIT-EST config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xcdhit.get_cd_hit_est_config_file())
        self.root.wait_window(dialog_editor)

        # check the CD-HIT-EST config file
        (OK, error_list) = xcdhit.check_cd_hit_est_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_cd_hit_est_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_cd_hit_est_process(self):
        '''
        Run a CD-HIT-EST process corresponding to the options in CD-HIT-EST config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_cd_hit_est_process" in "container" with the grid geometry manager
        form_run_cd_hit_est_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_cd_hit_est_code())
        form_run_cd_hit_est_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_cd_hit_est_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_cd_hit_est_process'
        self.forms_dict[self.current_form] = form_run_cd_hit_est_process

        # raise "form_run_cd_hit_est_process" to front
        form_run_cd_hit_est_process.tkraise()

    #---------------
    # Cuffdiff
    #---------------

    def recreate_cuffdiff_config_file(self):
        '''
        Recreate the Cuffdiff config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_cuffdiff_config_file" in "container" with the grid geometry manager
        form_recreate_cuffdiff_config_file = gbioinfoapp.FormRecreateCuffdiffConfigFile(self)
        form_recreate_cuffdiff_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_cuffdiff_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_cuffdiff_config_file'
        self.forms_dict[self.current_form] = form_recreate_cuffdiff_config_file

        # raise "form_recreate_cuffdiff_config_file" to front
        form_recreate_cuffdiff_config_file.tkraise()

    #---------------

    def edit_cuffdiff_config_file(self):
        '''
        Edit the Cuffdiff config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_cuffdiff_name()} - Edit config file'

        # edit the Cuffdiff config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xcufflinks.get_cuffdiff_config_file())
        self.root.wait_window(dialog_editor)

        # check the Cuffdiff config file
        (OK, error_list) = xcufflinks.check_cuffdiff_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_cuffdiff_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_cuffdiff_process(self):
        '''
        Run a Cuffdiff process corresponding to the options in Cuffdiff config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_cuffdiff_process" in "container" with the grid geometry manager
        form_run_cuffdiff_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_cuffdiff_code())
        form_run_cuffdiff_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_cuffdiff_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_cuffdiff_process'
        self.forms_dict[self.current_form] = form_run_cuffdiff_process

        # raise "form_run_cuffdiff_process" to front
        form_run_cuffdiff_process.tkraise()

    #---------------
    # Cufflinks-Cuffmerge
    #---------------

    def recreate_cufflinks_cuffmerge_config_file(self):
        '''
        Recreate the Cufflinks-Cuffmerge config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_cufflinks_cuffmerge_config_file" in "container" with the grid geometry manager
        form_recreate_cufflinks_cuffmerge_config_file = gbioinfoapp.FormRecreateCufflinksCuffmergeConfigFile(self)
        form_recreate_cufflinks_cuffmerge_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_cufflinks_cuffmerge_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_cufflinks_cuffmerge_config_file'
        self.forms_dict[self.current_form] = form_recreate_cufflinks_cuffmerge_config_file

        # raise "form_recreate_cufflinks_cuffmerge_config_file" to front
        form_recreate_cufflinks_cuffmerge_config_file.tkraise()

    #---------------

    def edit_cufflinks_cuffmerge_config_file(self):
        '''
        Edit the Cufflinks-Cuffmerge config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_cufflinks_cuffmerge_name()} - Edit config file'

        # edit the Cufflinks-Cuffmerge config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xcufflinks.get_cufflinks_cuffmerge_config_file())
        self.root.wait_window(dialog_editor)

        # check the Cufflinks-Cuffmerge config file
        (OK, error_list) = xcufflinks.check_cufflinks_cuffmerge_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_cufflinks_cuffmerge_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_cufflinks_cuffmerge_process(self):
        '''
        Run a Cufflinks-Cuffmerge process corresponding to the options in Cufflinks-Cuffmerge config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_cufflinks_cuffmerge_process" in "container" with the grid geometry manager
        form_run_cufflinks_cuffmerge_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_cufflinks_cuffmerge_code())
        form_run_cufflinks_cuffmerge_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_cufflinks_cuffmerge_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_cufflinks_cuffmerge_process'
        self.forms_dict[self.current_form] = form_run_cufflinks_cuffmerge_process

        # raise "form_run_cufflinks_cuffmerge_process" to front
        form_run_cufflinks_cuffmerge_process.tkraise()

    #---------------
    # Cuffquant
    #---------------

    def recreate_cuffquant_config_file(self):
        '''
        Recreate the Cuffquant config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_cuffquant_config_file" in "container" with the grid geometry manager
        form_recreate_cuffquant_config_file = gbioinfoapp.FormRecreateCuffquantConfigFile(self)
        form_recreate_cuffquant_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_cuffquant_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_cuffquant_config_file'
        self.forms_dict[self.current_form] = form_recreate_cuffquant_config_file

        # raise "form_recreate_cuffquant_config_file" to front
        form_recreate_cuffquant_config_file.tkraise()

    #---------------

    def edit_cuffquant_config_file(self):
        '''
        Edit the Cuffquant config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_cuffquant_name()} - Edit config file'

        # edit the Cuffquant config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xcufflinks.get_cuffquant_config_file())
        self.root.wait_window(dialog_editor)

        # check the Cuffquant config file
        (OK, error_list) = xcufflinks.check_cuffquant_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_cuffquant_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_cuffquant_process(self):
        '''
        Run a Cuffquant process corresponding to the options in CuffQuant config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_cuffquant_process" in "container" with the grid geometry manager
        form_run_cuffquant_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_cuffquant_code())
        form_run_cuffquant_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_cuffquant_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_cuffquant_process'
        self.forms_dict[self.current_form] = form_run_cuffquant_process

        # raise "form_run_cuffquant_process" to front
        form_run_cuffquant_process.tkraise()

    #---------------
    # cutadapt
    #---------------

    def recreate_cutadapt_config_file(self):
        '''
        Recreate the cutadapt config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_cutadapt_config_file" in "container" with the grid geometry manager
        form_recreate_cutadapt_config_file = gbioinfoapp.FormRecreateCutadaptConfigFile(self)
        form_recreate_cutadapt_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_cutadapt_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_cutadapt_config_file'
        self.forms_dict[self.current_form] = form_recreate_cutadapt_config_file

        # raise "form_recreate_cutadapt_config_file" to front
        form_recreate_cutadapt_config_file.tkraise()

    #---------------

    def edit_cutadapt_config_file(self):
        '''
        Edit the cutadapt config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_cutadapt_name()} - Edit config file'

        # edit the cutadapt config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xcutadapt.get_cutadapt_config_file())
        self.root.wait_window(dialog_editor)

        # check the cutadapt config file
        (OK, error_list) = xcutadapt.check_cutadapt_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_cutadapt_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_cutadapt_process(self):
        '''
        Run a cutadapt process corresponding to the options in cutadapt config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_cutadapt_process" in "container" with the grid geometry manager
        form_run_cutadapt_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_cutadapt_code())
        form_run_cutadapt_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_cutadapt_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_cutadapt_process'
        self.forms_dict[self.current_form] = form_run_cutadapt_process

        # raise "form_run_cutadapt_process" to front
        form_run_cutadapt_process.tkraise()

    #---------------
    # ddRADseq simulation
    #---------------

    def recreate_ddradseq_simulation_config_file(self):
        '''
        Recreate the ddRADseq simulation config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_ddradseq_simulation_config_file" in "container" with the grid geometry manager
        form_recreate_ddradseq_simulation_config_file = gbioinfoapp.FormRecreateDdRadSeqSimulationConfigFile(self)
        form_recreate_ddradseq_simulation_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_ddradseq_simulation_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_ddradseq_simulation_config_file'
        self.forms_dict[self.current_form] = form_recreate_ddradseq_simulation_config_file

        # raise "form_recreate_ddradseq_simulation_config_file" to front
        form_recreate_ddradseq_simulation_config_file.tkraise()

    #---------------

    def edit_ddradseq_simulation_config_file(self):
        '''
        Edit the ddRADseq simulation config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseq_simulation_name()} - Edit config file'

        # edit the ddRADseq simulation config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_ddradseq_simulation_config_file())
        self.root.wait_window(dialog_editor)

        # check the ddRADseq simulation config file
        (OK, error_list) = xddradseqtools.check_ddradseq_simulation_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_ddradseq_simulation_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_ddradseq_simulation_process(self):
        '''
        Run a ddRADseq simulation process corresponding to the options in ddRADseq simulation config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_ddradseq_simulation_process" in "container" with the grid geometry manager
        form_run_ddradseq_simulation_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_ddradseq_simulation_code())
        form_run_ddradseq_simulation_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_ddradseq_simulation_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_ddradseq_simulation_process'
        self.forms_dict[self.current_form] = form_run_ddradseq_simulation_process

        # raise "form_run_ddradseq_simulation_process" to front
        form_run_ddradseq_simulation_process.tkraise()

    #---------------

    def restart_ddradseq_simulation_process(self):
        '''
        Restart a ddRADseq simulation process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_ddradseq_simulation_process" in "container" with the grid geometry manager
        form_restart_ddradseq_simulation_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_ddradseq_simulation_code())
        form_restart_ddradseq_simulation_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_ddradseq_simulation_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_ddradseq_simulation_process'
        self.forms_dict[self.current_form] = form_restart_ddradseq_simulation_process

        # raise "form_restart_ddradseq_simulation_process" to front
        form_restart_ddradseq_simulation_process.tkraise()

    #---------------
    # eXpress
    #---------------

    def recreate_express_config_file(self):
        '''
        Recreate the eXpress config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_express_config_file" in "container" with the grid geometry manager
        form_recreate_express_config_file = gbioinfoapp.FormRecreateExpressConfigFile(self)
        form_recreate_express_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_express_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_express_config_file'
        self.forms_dict[self.current_form] = form_recreate_express_config_file

        # raise "form_recreate_express_config_file" to front
        form_recreate_express_config_file.tkraise()

    #---------------

    def edit_express_config_file(self):
        '''
        Edit the eXpress config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_express_name()} - Edit config file'

        # edit the eXpress config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xexpress.get_express_config_file())
        self.root.wait_window(dialog_editor)

        # check the eXpress config file
        (OK, error_list) = xexpress.check_express_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_express_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_express_process(self):
        '''
        Run a eXpress process corresponding to the options in eXpress config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_express_process" in "container" with the grid geometry manager
        form_run_express_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_express_code())
        form_run_express_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_express_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_express_process'
        self.forms_dict[self.current_form] = form_run_express_process

        # raise "form_run_express_process" to front
        form_run_express_process.tkraise()

    #---------------
    # FastQC
    #---------------

    def recreate_fastqc_config_file(self):
        '''
        Recreate the FastQC config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_fastqc_config_file" in "container" with the grid geometry manager
        form_recreate_fastqc_config_file = gbioinfoapp.FormRecreateFastQCConfigFile(self)
        form_recreate_fastqc_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_fastqc_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_fastqc_config_file'
        self.forms_dict[self.current_form] = form_recreate_fastqc_config_file

        # raise "form_recreate_fastqc_config_file" to front
        form_recreate_fastqc_config_file.tkraise()

    #---------------

    def edit_fastqc_config_file(self):
        '''
        Edit the FastQC config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_fastqc_name()} - Edit config file'

        # edit the FastQC config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xfastqc.get_fastqc_config_file())
        self.root.wait_window(dialog_editor)

        # check the FastQC config file
        (OK, error_list) = xfastqc.check_fastqc_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_fastqc_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_fastqc_process(self):
        '''
        Run a FastQC process corresponding to the options in FastQC config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_fastqc_process" in "container" with the grid geometry manager
        form_run_fastqc_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_fastqc_code())
        form_run_fastqc_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_fastqc_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_fastqc_process'
        self.forms_dict[self.current_form] = form_run_fastqc_process

        # raise "form_run_fastqc_process" to front
        form_run_fastqc_process.tkraise()

    #---------------
    # Genome-guided Trinity
    #---------------

    def recreate_ggtrinity_config_file(self):
        '''
        Create the Genome-guided Trinity config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_ggtrinity_config_file" in "container" with the grid geometry manager
        form_recreate_ggtrinity_config_file = gbioinfoapp.FormRecreateGgTrinityConfigFile(self)
        form_recreate_ggtrinity_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_ggtrinity_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_ggtrinity_config_file'
        self.forms_dict[self.current_form] = form_recreate_ggtrinity_config_file

        # raise "form_recreate_ggtrinity_config_file" to front
        form_recreate_ggtrinity_config_file.tkraise()

    #---------------

    def edit_ggtrinity_config_file(self):
        '''
        Edit the Genome-guided Trinity config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ggtrinity_name()} - Edit config file'

        # edit the Genome-guided Trinity config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtrinity.get_ggtrinity_config_file())
        self.root.wait_window(dialog_editor)

        # check the Genome-guided Trinity config file
        (OK, error_list) = xtrinity.check_ggtrinity_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_ggtrinity_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_ggtrinity_process(self):
        '''
        Run a Genome-guided Trinity process corresponding to the options in Genome-guided Trinity config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_ggtrinity_process" in "container" with the grid geometry manager
        form_run_ggtrinity_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_ggtrinity_code())
        form_run_ggtrinity_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_ggtrinity_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_ggtrinity_process'
        self.forms_dict[self.current_form] = form_run_ggtrinity_process

        # raise "form_run_ggtrinity_process" to front
        form_run_ggtrinity_process.tkraise()

    #---------------

    def restart_ggtrinity_process(self):
        '''
        Restart a Genome-guided Trinity process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_ggtrinity_process" in "container" with the grid geometry manager
        form_restart_ggtrinity_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_ggtrinity_code())
        form_restart_ggtrinity_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_ggtrinity_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_ggtrinity_process'
        self.forms_dict[self.current_form] = form_restart_ggtrinity_process

        # raise "form_restart_ggtrinity_process" to front
        form_restart_ggtrinity_process.tkraise()

    #---------------
    # GMAP
    #---------------

    def recreate_gmap_config_file(self):
        '''
        Recreate the GMAP config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_gmap_config_file" in "container" with the grid geometry manager
        form_recreate_gmap_config_file = gbioinfoapp.FormRecreateGmapConfigFile(self)
        form_recreate_gmap_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_gmap_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_gmap_config_file'
        self.forms_dict[self.current_form] = form_recreate_gmap_config_file

        # raise "form_recreate_gmap_config_file" to front
        form_recreate_gmap_config_file.tkraise()

    #---------------

    def edit_gmap_config_file(self):
        '''
        Edit the GMAP config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_gmap_name()} - Edit config file'

        # edit the GMAP config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgmap.get_gmap_config_file())
        self.root.wait_window(dialog_editor)

        # check the GMAP config file
        (OK, error_list) = xgmap.check_gmap_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_gmap_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_gmap_process(self):
        '''
        Run a GMAP process corresponding to the options in GMAP config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_gmap_process" in "container" with the grid geometry manager
        form_run_gmap_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_gmap_code())
        form_run_gmap_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_gmap_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_gmap_process'
        self.forms_dict[self.current_form] = form_run_gmap_process

        # raise "form_run_gmap_process" to front
        form_run_gmap_process.tkraise()

    #---------------
    # GSNAP
    #---------------

    def recreate_gsnap_config_file(self):
        '''
        Recreate the GSNAP config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_gsnap_config_file" in "container" with the grid geometry manager
        form_recreate_gsnap_config_file = gbioinfoapp.FormRecreateGsnapConfigFile(self)
        form_recreate_gsnap_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_gsnap_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_gsnap_config_file'
        self.forms_dict[self.current_form] = form_recreate_gsnap_config_file

        # raise "form_recreate_gsnap_config_file" to front
        form_recreate_gsnap_config_file.tkraise()

    #---------------

    def edit_gsnap_config_file(self):
        '''
        Edit the GSNAP config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_gsnap_name()} - Edit config file'

        # edit the GSNAP config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgmap.get_gsnap_config_file())
        self.root.wait_window(dialog_editor)

        # check the GSNAP config file
        (OK, error_list) = xgmap.check_gsnap_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_gsnap_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_gsnap_process(self):
        '''
        Run a GSNAP process corresponding to the options in GSNAP config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_gsnap_process" in "container" with the grid geometry manager
        form_run_gsnap_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_gsnap_code())
        form_run_gsnap_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_gsnap_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_gsnap_process'
        self.forms_dict[self.current_form] = form_run_gsnap_process

        # raise "form_run_gsnap_process" to front
        form_run_gsnap_process.tkraise()

    #---------------
    # HISAT2
    #---------------

    def recreate_hisat2_config_file(self):
        '''
        Recreate the HISAT2 config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_hisat2_config_file" in "container" with the grid geometry manager
        form_recreate_hisat2_config_file = gbioinfoapp.FormRecreateHisat2ConfigFile(self)
        form_recreate_hisat2_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_hisat2_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_hisat2_config_file'
        self.forms_dict[self.current_form] = form_recreate_hisat2_config_file

        # raise "form_recreate_hisat2_config_file" to front
        form_recreate_hisat2_config_file.tkraise()

    #---------------

    def edit_hisat2_config_file(self):
        '''
        Edit the HISAT2 config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_hisat2_name()} - Edit config file'

        # edit the HISAT2 config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xhisat2.get_hisat2_config_file())
        self.root.wait_window(dialog_editor)

        # check the HISAT2 config file
        (OK, error_list) = xhisat2.check_hisat2_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_hisat2_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_hisat2_process(self):
        '''
        Run a HISAT2 process corresponding to the options in HISAT2 config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_hisat2_process" in "container" with the grid geometry manager
        form_run_hisat2_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_hisat2_code())
        form_run_hisat2_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_hisat2_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_hisat2_process'
        self.forms_dict[self.current_form] = form_run_hisat2_process

        # raise "form_run_hisat2_process" to front
        form_run_hisat2_process.tkraise()

    #---------------
    # htseq-count
    #---------------

    def recreate_htseq_count_config_file(self):
        '''
        Recreate the htseq-count config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_htseq_count_config_file" in "container" with the grid geometry manager
        form_recreate_htseq_count_config_file = gbioinfoapp.FormRecreateHtseqCountConfigFile(self)
        form_recreate_htseq_count_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_htseq_count_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_htseq_count_config_file'
        self.forms_dict[self.current_form] = form_recreate_htseq_count_config_file

        # raise "form_recreate_htseq_count_config_file" to front
        form_recreate_htseq_count_config_file.tkraise()

    #---------------

    def edit_htseq_count_config_file(self):
        '''
        Edit the htseq-count config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_htseq_count_name()} - Edit config file'

        # edit the htseq-count config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xhtseq.get_htseq_count_config_file())
        self.root.wait_window(dialog_editor)

        # check the htseq-count config file
        (OK, error_list) = xhtseq.check_htseq_count_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_htseq_count_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_htseq_count_process(self):
        '''
        Run a htseq-count process corresponding to the options in htseq-count config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_htseq_count_process" in "container" with the grid geometry manager
        form_run_htseq_count_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_htseq_count_code())
        form_run_htseq_count_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_htseq_count_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_htseq_count_process'
        self.forms_dict[self.current_form] = form_run_htseq_count_process

        # raise "form_run_htseq_count_process" to front
        form_run_htseq_count_process.tkraise()

    #---------------
    # insilico_read_normalization
    #---------------

    def recreate_insilico_read_normalization_config_file(self):
        '''
        Create the insilico_read_normalization config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_insilico_read_normalization_config_file" in "container" with the grid geometry manager
        form_recreate_insilico_read_normalization_config_file = gbioinfoapp.FormRecreateInsilicoReadNormalizationConfigFile(self)
        form_recreate_insilico_read_normalization_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_insilico_read_normalization_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_insilico_read_normalization_config_file'
        self.forms_dict[self.current_form] = form_recreate_insilico_read_normalization_config_file

        # raise "form_recreate_insilico_read_normalization_config_file" to front
        form_recreate_insilico_read_normalization_config_file.tkraise()

    #---------------

    def edit_insilico_read_normalization_config_file(self):
        '''
        Edit the insilico_read_normalization config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_insilico_read_normalization_name()} - Edit config file'

        # edit the insilico_read_normalization config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtrinity.get_insilico_read_normalization_config_file())
        self.root.wait_window(dialog_editor)

        # check the insilico_read_normalization config file
        (OK, error_list) = xtrinity.check_insilico_read_normalization_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_insilico_read_normalization_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_insilico_read_normalization_process(self):
        '''
        Run a insilico_read_normalization process corresponding to the options in insilico_read_normalization config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_insilico_read_normalization_process" in "container" with the grid geometry manager
        form_run_insilico_read_normalization_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_insilico_read_normalization_code())
        form_run_insilico_read_normalization_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_insilico_read_normalization_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_insilico_read_normalization_process'
        self.forms_dict[self.current_form] = form_run_insilico_read_normalization_process

        # raise "form_run_insilico_read_normalization_process" to front
        form_run_insilico_read_normalization_process.tkraise()

    #---------------

    def restart_insilico_read_normalization_process(self):
        '''
        Restart a insilico_read_normalization process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_insilico_read_normalization_process" in "container" with the grid geometry manager
        form_restart_insilico_read_normalization_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_insilico_read_normalization_code())
        form_restart_insilico_read_normalization_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_insilico_read_normalization_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_insilico_read_normalization_process'
        self.forms_dict[self.current_form] = form_restart_insilico_read_normalization_process

        # raise "form_restart_insilico_read_normalization_process" to front
        form_restart_insilico_read_normalization_process.tkraise()

    #---------------
    # ipyrad
    #---------------

    def recreate_ipyrad_config_file(self):
        '''
        Recreate the ipyrad config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_ipyrad_config_file" in "container" with the grid geometry manager
        form_recreate_ipyrad_config_file = gbioinfoapp.FormRecreateIpyradConfigFile(self)
        form_recreate_ipyrad_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_ipyrad_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_ipyrad_config_file'
        self.forms_dict[self.current_form] = form_recreate_ipyrad_config_file

        # raise "form_recreate_ipyrad_config_file" to front
        form_recreate_ipyrad_config_file.tkraise()

    #---------------

    def edit_ipyrad_config_file(self):
        '''
        Edit the ipyrad config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ipyrad_name()} - Edit config file'

        # edit the ipyrad config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xipyrad.get_ipyrad_config_file())
        self.root.wait_window(dialog_editor)

        # check the ipyrad config file
        (OK, error_list) = xipyrad.check_ipyrad_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_ipyrad_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_ipyrad_process(self):
        '''
        Run a ipyrad process corresponding to the options in ipyrad config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_ipyrad_process" in "container" with the grid geometry manager
        form_run_ipyrad_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_ipyrad_code())
        form_run_ipyrad_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_ipyrad_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_ipyrad_process'
        self.forms_dict[self.current_form] = form_run_ipyrad_process

        # raise "form_run_ipyrad_process" to front
        form_run_ipyrad_process.tkraise()

    #---------------
    # kallisto
    #---------------

    def recreate_kallisto_config_file(self):
        '''
        Recreate the kallisto config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_kallisto_config_file" in "container" with the grid geometry manager
        form_recreate_kallisto_config_file = gbioinfoapp.FormRecreateKallistoConfigFile(self)
        form_recreate_kallisto_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_kallisto_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_kallisto_config_file'
        self.forms_dict[self.current_form] = form_recreate_kallisto_config_file

        # raise "form_recreate_kallisto_config_file" to front
        form_recreate_kallisto_config_file.tkraise()

    #---------------

    def edit_kallisto_config_file(self):
        '''
        Edit the kallisto config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_kallisto_name()} - Edit config file'

        # edit the kallisto config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xkallisto.get_kallisto_config_file())
        self.root.wait_window(dialog_editor)

        # check the kallisto config file
        (OK, error_list) = xkallisto.check_kallisto_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_kallisto_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_kallisto_process(self):
        '''
        Run a kallisto process corresponding to the options in kallisto config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_kallisto_process" in "container" with the grid geometry manager
        form_run_kallisto_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_kallisto_code())
        form_run_kallisto_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_kallisto_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_kallisto_process'
        self.forms_dict[self.current_form] = form_run_kallisto_process

        # raise "form_run_kallisto_process" to front
        form_run_kallisto_process.tkraise()

    #---------------
    # QUAST
    #---------------

    def recreate_quast_config_file(self):
        '''
        Recreate the QUAST config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_quast_config_file" in "container" with the grid geometry manager
        form_recreate_quast_config_file = gbioinfoapp.FormRecreateQuastConfigFile(self)
        form_recreate_quast_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_quast_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_quast_config_file'
        self.forms_dict[self.current_form] = form_recreate_quast_config_file

        # raise "form_recreate_quast_config_file" to front
        form_recreate_quast_config_file.tkraise()

    #---------------

    def edit_quast_config_file(self):
        '''
        Edit the QUAST config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_quast_name()} - Edit config file'

        # edit the QUAST config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xquast.get_quast_config_file())
        self.root.wait_window(dialog_editor)

        # check the QUAST config file
        (OK, error_list) = xquast.check_quast_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_quast_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_quast_process(self):
        '''
        Run a QUAST process corresponding to the options in QUAST config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_quast_process" in "container" with the grid geometry manager
        form_run_quast_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_quast_code())
        form_run_quast_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_quast_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_quast_process'
        self.forms_dict[self.current_form] = form_run_quast_process

        # raise "form_run_quast_process" to front
        form_run_quast_process.tkraise()

    #---------------
    # RADdesigner
    #---------------

    def recreate_raddesigner_config_file(self):
        '''
        Recreate the raddesigner config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_raddesigner_config_file" in "container" with the grid geometry manager
        form_recreate_raddesigner_config_file = gbioinfoapp.FormRecreateRADdesignerConfigFile(self)
        form_recreate_raddesigner_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_raddesigner_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_raddesigner_config_file'
        self.forms_dict[self.current_form] = form_recreate_raddesigner_config_file

        # raise "form_recreate_raddesigner_config_file" to front
        form_recreate_raddesigner_config_file.tkraise()

    #---------------

    def edit_raddesigner_config_file(self):
        '''
        Edit the RADdesigner config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_raddesigner_name()} - Edit config file'

        # edit the RADdesigner config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xraddesigner.get_raddesigner_config_file())
        self.root.wait_window(dialog_editor)

        # check the RADdesigner config file
        (OK, error_list) = xraddesigner.check_raddesigner_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_raddesigner_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_raddesigner_process(self):
        '''
        Run a RADdesigner process corresponding to the options in RADdesigner config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_raddesigner_process" in "container" with the grid geometry manager
        form_run_raddesigner_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_raddesigner_code())
        form_run_raddesigner_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_raddesigner_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_raddesigner_process'
        self.forms_dict[self.current_form] = form_run_raddesigner_process

        # raise "form_run_raddesigner_process" to front
        form_run_raddesigner_process.tkraise()

    #---------------

    def restart_raddesigner_process(self):
        '''
        Restart a RADdesigner process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_raddesigner_process" in "container" with the grid geometry manager
        form_restart_raddesigner_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_raddesigner_code())
        form_restart_raddesigner_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_raddesigner_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_raddesigner_process'
        self.forms_dict[self.current_form] = form_restart_raddesigner_process

        # raise "form_restart_raddesigner_process" to front
        form_restart_raddesigner_process.tkraise()

    #---------------
    # REF-EVAL
    #---------------

    #def recreate_ref_eval_config_file(self):
    #    '''
    #    Create the REF-EVAL config file with the default options. It is necessary
    #    update the options in process run.
    #    '''

    #    # close the current form
    #    self.close_current_form()

    #    # create and register "form_recreate_ref_eval_config_file" in "container" with the grid geometry manager
    #    form_recreate_ref_eval_config_file = gbioinfoapp.FormRecreateRefEvalConfigFile(self)
    #    form_recreate_ref_eval_config_file.grid(row=0, column=0, sticky='nsew')

    #    # set "form_recreate_ref_eval_config_file" as current form and add it in the forms dictionary
    #    self.current_form = 'form_recreate_ref_eval_config_file'
    #    self.forms_dict[self.current_form] = form_recreate_ref_eval_config_file

    #    # raise "form_recreate_ref_eval_config_file" to front
    #    form_recreate_ref_eval_config_file.tkraise()

    #---------------

    #def edit_ref_eval_config_file(self):
    #    '''
    #    Edit the REF-EVAL config file to change the parameters of process run.
    #    '''

    #    # initialize the control variable
    #    OK = True

    #    # close the current form
    #    self.close_current_form()

    #    # set the head
    #    head = f'{xlib.get_ref_eval_name()} - Edit config file'

    #    # edit the RSEM-EVAL config file using "DialogEditor" 
    #    dialog_editor = gdialogs.DialogEditor(self.root, xdetonate.get_ref_eval_config_file())
    #    self.root.wait_window(dialog_editor)

    #    # check the RSEM-EVAL config file
    #    (OK, error_list) = xdetonate.check_ref_eval_config_file(strict=False)
    #    if OK:
    #        message = f'The {xlib.get_ref_eval_name()} config file is OK.'
    #        tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
    #    else:
    #        message = 'Detected errors:\n\n'
    #        for error in error_list:
    #            message = f'{message}{error}\n'
    #        tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    #def run_ref_eval_process(self):
    #    '''
    #    Run a REF-EVAL process corresponding to the options in REF-EVAL config file.
    #    '''

    #    # close the current form
    #    self.close_current_form()

    #    # create and register "form_run_ref_eval_process" in "container" with the grid geometry manager
    #    form_run_ref_eval_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_ref_eval_code())
    #    form_run_ref_eval_process.grid(row=0, column=0, sticky='nsew')

    #    # set "form_run_ref_eval_process" as current form and add it in the forms dictionary
    #    self.current_form = 'form_run_ref_eval_process'
    #    self.forms_dict[self.current_form] = form_run_ref_eval_process

    #    # raise "form_run_ref_eval_process" to front
    #    form_run_ref_eval_process.tkraise()

    #---------------
    # rnaQUAST
    #---------------

    def recreate_rnaquast_config_file(self):
        '''
        Recreate the rnaQUAST config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_rnaquast_config_file" in "container" with the grid geometry manager
        form_recreate_rnaquast_config_file = gbioinfoapp.FormRecreateRnaQuastConfigFile(self)
        form_recreate_rnaquast_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_rnaquast_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_rnaquast_config_file'
        self.forms_dict[self.current_form] = form_recreate_rnaquast_config_file

        # raise "form_recreate_rnaquast_config_file" to front
        form_recreate_rnaquast_config_file.tkraise()

    #---------------

    def edit_rnaquast_config_file(self):
        '''
        Edit the rnaQUAST config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_rnaquast_name()} - Edit config file'

        # edit the rnaQUAST config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xrnaquast.get_rnaquast_config_file())
        self.root.wait_window(dialog_editor)

        # check the rnaQUAST config file
        (OK, error_list) = xrnaquast.check_rnaquast_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_rnaquast_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_rnaquast_process(self):
        '''
        Run a rnaQUAST process corresponding to the options in rnaQUAST config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_rnaquast_process" in "container" with the grid geometry manager
        form_run_rnaquast_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_rnaquast_code())
        form_run_rnaquast_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_rnaquast_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_rnaquast_process'
        self.forms_dict[self.current_form] = form_run_rnaquast_process

        # raise "form_run_rnaquast_process" to front
        form_run_rnaquast_process.tkraise()

    #---------------
    # RSEM-EVAL
    #---------------

    def recreate_rsem_eval_config_file(self):
        '''
        Create the RSEM-EVAL config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_rsem_eval_config_file" in "container" with the grid geometry manager
        form_recreate_rsem_eval_config_file = gbioinfoapp.FormRecreateRsemEvalConfigFile(self)
        form_recreate_rsem_eval_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_rsem_eval_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_rsem_eval_config_file'
        self.forms_dict[self.current_form] = form_recreate_rsem_eval_config_file

        # raise "form_recreate_rsem_eval_config_file" to front
        form_recreate_rsem_eval_config_file.tkraise()

    #---------------

    def edit_rsem_eval_config_file(self):
        '''
        Edit the RSEM-EVAL config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_rsem_eval_name()} - Edit config file'

        # edit the RSEM-EVAL config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xdetonate.get_rsem_eval_config_file())
        self.root.wait_window(dialog_editor)

        # check the RSEM-EVAL config file
        (OK, error_list) = xdetonate.check_rsem_eval_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_rsem_eval_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_rsem_eval_process(self):
        '''
        Run a RSEM-EVAL process corresponding to the options in RSEM-EVAL config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_rsem_eval_process" in "container" with the grid geometry manager
        form_run_rsem_eval_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_rsem_eval_code())
        form_run_rsem_eval_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_rsem_eval_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_rsem_eval_process'
        self.forms_dict[self.current_form] = form_run_rsem_eval_process

        # raise "form_run_rsem_eval_process" to front
        form_run_rsem_eval_process.tkraise()

    #---------------
    # rsitesearch
    #---------------

    def recreate_rsitesearch_config_file(self):
        '''
        Recreate the rsitesearch config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_rsitesearch_config_file" in "container" with the grid geometry manager
        form_recreate_rsitesearch_config_file = gbioinfoapp.FormRecreateRSiteSearchConfigFile(self)
        form_recreate_rsitesearch_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_rsitesearch_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_rsitesearch_config_file'
        self.forms_dict[self.current_form] = form_recreate_rsitesearch_config_file

        # raise "form_recreate_rsitesearch_config_file" to front
        form_recreate_rsitesearch_config_file.tkraise()

    #---------------

    def edit_rsitesearch_config_file(self):
        '''
        Edit the rsitesearch config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_rsitesearch_name()} - Edit config file'

        # edit the rsitesearch config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_rsitesearch_config_file())
        self.root.wait_window(dialog_editor)

        # check the rsitesearch config file
        (OK, error_list) = xddradseqtools.check_rsitesearch_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_rsitesearch_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_rsitesearch_process(self):
        '''
        Run a rsitesearch process corresponding to the options in rsitesearch config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_rsitesearch_process" in "container" with the grid geometry manager
        form_run_rsitesearch_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_rsitesearch_code())
        form_run_rsitesearch_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_rsitesearch_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_rsitesearch_process'
        self.forms_dict[self.current_form] = form_run_rsitesearch_process

        # raise "form_run_rsitesearch_process" to front
        form_run_rsitesearch_process.tkraise()

    #---------------
    # SOAPdenovo2
    #---------------

    def recreate_soapdenovo2_config_file(self):
        '''
        Recreate the SOAPdenovo2 config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_soapdenovo2_config_file" in "container" with the grid geometry manager
        form_recreate_soapdenovo2_config_file = gbioinfoapp.FormRecreateSoapdenovo2ConfigFile(self)
        form_recreate_soapdenovo2_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_soapdenovo2_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_soapdenovo2_config_file'
        self.forms_dict[self.current_form] = form_recreate_soapdenovo2_config_file

        # raise "form_recreate_soapdenovo2_config_file" to front
        form_recreate_soapdenovo2_config_file.tkraise()

    #---------------

    def edit_soapdenovo2_config_file(self):
        '''
        Edit the SOAPdenovo2 config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_soapdenovo2_name()} - Edit config file'

        # edit the SOAPdenovo2 config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xsoapdenovo2.get_soapdenovo2_config_file())
        self.root.wait_window(dialog_editor)

        # check the SOAPdenovo2 config file
        (OK, error_list) = xsoapdenovo2.check_soapdenovo2_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_soapdenovo2_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_soapdenovo2_process(self):
        '''
        Run a SOAPdenovo2 process corresponding to the options in SOAPdenovo2 config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_soapdenovo2_process" in "container" with the grid geometry manager
        form_run_soapdenovo2_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_soapdenovo2_code())
        form_run_soapdenovo2_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_soapdenovo2_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_soapdenovo2_process'
        self.forms_dict[self.current_form] = form_run_soapdenovo2_process

        # raise "form_run_soapdenovo2_process" to front
        form_run_soapdenovo2_process.tkraise()

    #---------------

    def restart_soapdenovo2_process(self):
        '''
        Restart a SOAPdenovo2 process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_soapdenovo2_process" in "container" with the grid geometry manager
        form_restart_soapdenovo2_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_soapdenovo2_code())
        form_restart_soapdenovo2_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_soapdenovo2_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_soapdenovo2_process'
        self.forms_dict[self.current_form] = form_restart_soapdenovo2_process

        # raise "form_restart_soapdenovo2_process" to front
        form_restart_soapdenovo2_process.tkraise()

    #---------------
    # SOAPdenovo-Trans
    #---------------

    def recreate_soapdenovotrans_config_file(self):
        '''
        Recreate the SOAPdenovo-Trans config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_soapdenovotrans_config_file" in "container" with the grid geometry manager
        form_recreate_soapdenovotrans_config_file = gbioinfoapp.FormRecreateSoapdenovoTransConfigFile(self)
        form_recreate_soapdenovotrans_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_soapdenovotrans_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_soapdenovotrans_config_file'
        self.forms_dict[self.current_form] = form_recreate_soapdenovotrans_config_file

        # raise "form_recreate_soapdenovotrans_config_file" to front
        form_recreate_soapdenovotrans_config_file.tkraise()

    #---------------

    def edit_soapdenovotrans_config_file(self):
        '''
        Edit the SOAPdenovo-Trans config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_soapdenovotrans_name()} - Edit config file'

        # edit the SOAPdenovo-Trans config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xsoapdenovotrans.get_soapdenovotrans_config_file())
        self.root.wait_window(dialog_editor)

        # check the SOAPdenovo-Trans config file
        (OK, error_list) = xsoapdenovotrans.check_soapdenovotrans_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_soapdenovotrans_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_soapdenovotrans_process(self):
        '''
        Run a SOAPdenovo-Trans process corresponding to the options in SOAPdenovo-Trans config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_soapdenovotrans_process" in "container" with the grid geometry manager
        form_run_soapdenovotrans_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_soapdenovotrans_code())
        form_run_soapdenovotrans_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_soapdenovotrans_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_soapdenovotrans_process'
        self.forms_dict[self.current_form] = form_run_soapdenovotrans_process

        # raise "form_run_soapdenovotrans_process" to front
        form_run_soapdenovotrans_process.tkraise()

    #---------------

    def restart_soapdenovotrans_process(self):
        '''
        Restart a SOAPdenovo-Trans process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_soapdenovotrans_process" in "container" with the grid geometry manager
        form_restart_soapdenovotrans_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_soapdenovotrans_code())
        form_restart_soapdenovotrans_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_soapdenovotrans_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_soapdenovotrans_process'
        self.forms_dict[self.current_form] = form_restart_soapdenovotrans_process

        # raise "form_restart_soapdenovotrans_process" to front
        form_restart_soapdenovotrans_process.tkraise()

    #---------------
    # STAR
    #---------------

    def recreate_star_config_file(self):
        '''
        Recreate the STAR config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_star_config_file" in "container" with the grid geometry manager
        form_recreate_star_config_file = gbioinfoapp.FormRecreateSTARConfigFile(self)
        form_recreate_star_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_star_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_star_config_file'
        self.forms_dict[self.current_form] = form_recreate_star_config_file

        # raise "form_recreate_star_config_file" to front
        form_recreate_star_config_file.tkraise()

    #---------------

    def edit_star_config_file(self):
        '''
        Edit the STAR config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_star_name()} - Edit config file'

        # edit the STAR config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xstar.get_star_config_file())
        self.root.wait_window(dialog_editor)

        # check the STAR config file
        (OK, error_list) = xstar.check_star_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_star_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_star_process(self):
        '''
        Run a STAR process corresponding to the options in STAR config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_star_process" in "container" with the grid geometry manager
        form_run_star_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_star_code())
        form_run_star_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_star_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_star_process'
        self.forms_dict[self.current_form] = form_run_star_process

        # raise "form_run_star_process" to front
        form_run_star_process.tkraise()

    #---------------
    # starcode
    #---------------

    def recreate_starcode_config_file(self):
        '''
        Recreate the starcode config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_starcode_config_file" in "container" with the grid geometry manager
        form_recreate_starcode_config_file = gbioinfoapp.FormRecreateStarcodeConfigFile(self)
        form_recreate_starcode_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_starcode_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_starcode_config_file'
        self.forms_dict[self.current_form] = form_recreate_starcode_config_file

        # raise "form_recreate_starcode_config_file" to front
        form_recreate_starcode_config_file.tkraise()

    #---------------

    def edit_starcode_config_file(self):
        '''
        Edit the starcode config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_starcode_name()} - Edit config file'

        # edit the starcode config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xstarcode.get_starcode_config_file())
        self.root.wait_window(dialog_editor)

        # check the starcode config file
        (OK, error_list) = xstarcode.check_starcode_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_starcode_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_starcode_process(self):
        '''
        Run a starcode process corresponding to the options in starcode config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_starcode_process" in "container" with the grid geometry manager
        form_run_starcode_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_starcode_code())
        form_run_starcode_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_starcode_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_starcode_process'
        self.forms_dict[self.current_form] = form_run_starcode_process

        # raise "form_run_starcode_process" to front
        form_run_starcode_process.tkraise()

    #---------------
    # TopHat
    #---------------

    def recreate_tophat_config_file(self):
        '''
        Recreate the TopHat config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_tophat_config_file" in "container" with the grid geometry manager
        form_recreate_tophat_config_file = gbioinfoapp.FormRecreateTopHatConfigFile(self)
        form_recreate_tophat_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_tophat_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_tophat_config_file'
        self.forms_dict[self.current_form] = form_recreate_tophat_config_file

        # raise "form_recreate_tophat_config_file" to front
        form_recreate_tophat_config_file.tkraise()

    #---------------

    def edit_tophat_config_file(self):
        '''
        Edit the TopHat config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_tophat_name()} - Edit config file'

        # edit the TopHat config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtophat.get_tophat_config_file())
        self.root.wait_window(dialog_editor)

        # check the TopHat config file
        (OK, error_list) = xtophat.check_tophat_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_tophat_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_tophat_process(self):
        '''
        Run a TopHat process corresponding to the options in TopHat config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_tophat_process" in "container" with the grid geometry manager
        form_run_tophat_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_tophat_code())
        form_run_tophat_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_tophat_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_tophat_process'
        self.forms_dict[self.current_form] = form_run_tophat_process

        # raise "form_run_tophat_process" to front
        form_run_tophat_process.tkraise()

    #---------------
    # Transrate
    #---------------

    def recreate_transrate_config_file(self):
        '''
        Recreate the Transrate config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_transrate_config_file" in "container" with the grid geometry manager
        form_recreate_transrate_config_file = gbioinfoapp.FormRecreateTransrateConfigFile(self)
        form_recreate_transrate_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_transrate_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_transrate_config_file'
        self.forms_dict[self.current_form] = form_recreate_transrate_config_file

        # raise "form_recreate_transrate_config_file" to front
        form_recreate_transrate_config_file.tkraise()

    #---------------

    def edit_transrate_config_file(self):
        '''
        Edit the Transrate config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_transrate_name()} - Edit config file'

        # edit the Transrate config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtransrate.get_transrate_config_file())
        self.root.wait_window(dialog_editor)

        # check the Transrate config file
        (OK, error_list) = xtransrate.check_transrate_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_transrate_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_transrate_process(self):
        '''
        Run a Transrate process corresponding to the options in Transrate config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_transrate_process" in "container" with the grid geometry manager
        form_run_transrate_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_transrate_code())
        form_run_transrate_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_transrate_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_transrate_process'
        self.forms_dict[self.current_form] = form_run_transrate_process

        # raise "form_run_transrate_process" to front
        form_run_transrate_process.tkraise()

    #---------------
    # Trans-ABySS
    #---------------

    def recreate_transabyss_config_file(self):
        '''
        Recreate the Trans-ABySS config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_transabyss_config_file" in "container" with the grid geometry manager
        form_recreate_transabyss_config_file = gbioinfoapp.FormRecreateTransAbyssConfigFile(self)
        form_recreate_transabyss_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_transabyss_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_transabyss_config_file'
        self.forms_dict[self.current_form] = form_recreate_transabyss_config_file

        # raise "form_recreate_transabyss_config_file" to front
        form_recreate_transabyss_config_file.tkraise()

    #---------------

    def edit_transabyss_config_file(self):
        '''
        Edit the Trans-ABySS config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_transabyss_name()} - Edit config file'

        # edit the Trans-ABySS config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtransabyss.get_transabyss_config_file())
        self.root.wait_window(dialog_editor)

        # check the Trans-ABySS config file
        (OK, error_list) = xtransabyss.check_transabyss_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_transabyss_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_transabyss_process(self):
        '''
        Run a Trans-ABySS process corresponding to the options in Trans-ABySS config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_transabyss_process" in "container" with the grid geometry manager
        form_run_transabyss_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_transabyss_code())
        form_run_transabyss_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_transabyss_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_transabyss_process'
        self.forms_dict[self.current_form] = form_run_transabyss_process

        # raise "form_run_transabyss_process" to front
        form_run_transabyss_process.tkraise()

    #---------------
    # transcripts-filter
    #---------------

    def recreate_transcript_filter_config_file(self):
        '''
        Create the transcripts-filter config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_transcript_filter_config_file" in "container" with the grid geometry manager
        form_transcript_filter_eval_config_file = gbioinfoapp.FormRecreateTranscriptFilterConfigFile(self)
        form_transcript_filter_eval_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_transcript_filter_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_transcript_filter_config_file'
        self.forms_dict[self.current_form] = form_transcript_filter_eval_config_file

        # raise "form_recreate_transcript_filter_config_file" to front
        form_transcript_filter_eval_config_file.tkraise()

    #---------------

    def edit_transcript_filter_config_file(self):
        '''
        Edit the transcripts-filter config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_transcript_filter_name()} - Edit config file'

        # edit the transcripts-filter config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xngshelper.get_transcript_filter_config_file())
        self.root.wait_window(dialog_editor)

        # check the transcripts-filter config file
        (OK, error_list) = xngshelper.check_transcript_filter_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_transcript_filter_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_transcript_filter_process(self):
        '''
        Run a transcripts-filter process corresponding to the options in transcripts-filter config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_transcript_filter_process" in "container" with the grid geometry manager
        form_run_transcript_filter_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_transcript_filter_code())
        form_run_transcript_filter_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_transcript_filter_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_transcript_filter_process'
        self.forms_dict[self.current_form] = form_run_transcript_filter_process

        # raise "form_run_transcript_filter_process" to front
        form_run_transcript_filter_process.tkraise()

    #---------------
    # transcriptome-blastx
    #---------------

    def recreate_transcriptome_blastx_config_file(self):
        '''
        Create the transcriptome-blastx config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_transcriptome_blastx_config_file" in "container" with the grid geometry manager
        form_transcriptome_blastx_eval_config_file = gbioinfoapp.FormRecreateTranscriptomeBlastxConfigFile(self)
        form_transcriptome_blastx_eval_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_transcriptome_blastx_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_transcriptome_blastx_config_file'
        self.forms_dict[self.current_form] = form_transcriptome_blastx_eval_config_file

        # raise "form_recreate_transcriptome_blastx_config_file" to front
        form_transcriptome_blastx_eval_config_file.tkraise()

    #---------------

    def edit_transcriptome_blastx_config_file(self):
        '''
        Edit the transcriptome-blastx config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_transcriptome_blastx_name()} - Edit config file'

        # edit the transcripts-filter config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xngshelper.get_transcriptome_blastx_config_file())
        self.root.wait_window(dialog_editor)

        # check the transcripts-filter config file
        (OK, error_list) = xngshelper.check_transcriptome_blastx_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_transcriptome_blastx_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_transcriptome_blastx_process(self):
        '''
        Run a transcriptome-blastx process corresponding to the options in RSEM-EVAL config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_transcriptome_blastx_process" in "container" with the grid geometry manager
        form_run_transcriptome_blastx_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_transcriptome_blastx_code())
        form_run_transcriptome_blastx_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_transcriptome_blastx_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_transcriptome_blastx_process'
        self.forms_dict[self.current_form] = form_run_transcriptome_blastx_process

        # raise "form_run_transcriptome_blastx_process" to front
        form_run_transcriptome_blastx_process.tkraise()

    #---------------
    # Trimmomatic
    #---------------

    def recreate_trimmomatic_config_file(self):
        '''
        Recreate the Trimmomatic config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_trimmomatic_config_file" in "container" with the grid geometry manager
        form_recreate_trimmomatic_config_file = gbioinfoapp.FormRecreateTrimmomaticConfigFile(self)
        form_recreate_trimmomatic_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_trimmomatic_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_trimmomatic_config_file'
        self.forms_dict[self.current_form] = form_recreate_trimmomatic_config_file

        # raise "form_recreate_trimmomatic_config_file" to front
        form_recreate_trimmomatic_config_file.tkraise()

    #---------------

    def edit_trimmomatic_config_file(self):
        '''
        Edit the Trimmomatic config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_trimmomatic_name()} - Edit config file'

        # edit the Trimmomatic config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtrimmomatic.get_trimmomatic_config_file())
        self.root.wait_window(dialog_editor)

        # check the Trimmomatic config file
        (OK, error_list) = xtrimmomatic.check_trimmomatic_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_trimmomatic_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_trimmomatic_process(self):
        '''
        Run a Trimmomatic process corresponding to the options in Trimmomatic config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_trimmomatic_process" in "container" with the grid geometry manager
        form_run_trimmomatic_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_trimmomatic_code())
        form_run_trimmomatic_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_trimmomatic_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_trimmomatic_process'
        self.forms_dict[self.current_form] = form_run_trimmomatic_process

        # raise "form_run_trimmomatic_process" to front
        form_run_trimmomatic_process.tkraise()

    #---------------
    # Trinity
    #---------------

    def recreate_trinity_config_file(self):
        '''
        Create the Trinity config file with the default options. It is necessary
        update the options in process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_trinity_config_file" in "container" with the grid geometry manager
        form_recreate_trinity_config_file = gbioinfoapp.FormRecreateTrinityConfigFile(self)
        form_recreate_trinity_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_trinity_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_trinity_config_file'
        self.forms_dict[self.current_form] = form_recreate_trinity_config_file

        # raise "form_recreate_trinity_config_file" to front
        form_recreate_trinity_config_file.tkraise()

    #---------------

    def edit_trinity_config_file(self):
        '''
        Edit the Trinity config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_trinity_name()} - Edit config file'

        # edit the Trinity config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtrinity.get_trinity_config_file())
        self.root.wait_window(dialog_editor)

        # check the Trinity config file
        (OK, error_list) = xtrinity.check_trinity_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_trinity_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_trinity_process(self):
        '''
        Run a Trinity process corresponding to the options in Trinity config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_trinity_process" in "container" with the grid geometry manager
        form_run_trinity_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_trinity_code())
        form_run_trinity_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_trinity_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_trinity_process'
        self.forms_dict[self.current_form] = form_run_trinity_process

        # raise "form_run_trinity_process" to front
        form_run_trinity_process.tkraise()

    #---------------

    def restart_trinity_process(self):
        '''
        Restart a Trinity process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_trinity_process" in "container" with the grid geometry manager
        form_restart_trinity_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_trinity_code())
        form_restart_trinity_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_trinity_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_trinity_process'
        self.forms_dict[self.current_form] = form_restart_trinity_process

        # raise "form_restart_trinity_process" to front
        form_restart_trinity_process.tkraise()

    #---------------
    # variant_calling
    #---------------

    def recreate_variant_calling_config_file(self):
        '''
        Recreate the Variant calling config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_variant_calling_config_file" in "container" with the grid geometry manager
        form_recreate_variant_calling_config_file = gbioinfoapp.FormRecreateVariantCallingConfigFile(self)
        form_recreate_variant_calling_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_variant_calling_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_variant_calling_config_file'
        self.forms_dict[self.current_form] = form_recreate_variant_calling_config_file

        # raise "form_recreate_variant_calling_config_file" to front
        form_recreate_variant_calling_config_file.tkraise()

    #---------------

    def edit_variant_calling_config_file(self):
        '''
        Edit the Variant calling config file to change the parameters of each process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_variant_calling_name()} - Edit config file'

        # edit the Variant calling config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_variant_calling_config_file())
        self.root.wait_window(dialog_editor)

        # check the Variant calling config file
        (OK, error_list) = xddradseqtools.check_variant_calling_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_variant_calling_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_variant_calling_process(self):
        '''
        Run a Variant calling process corresponding to the options in Variant calling config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_variant_calling_process" in "container" with the grid geometry manager
        form_run_variant_calling_process = gbioinfoapp.FormRunBioinfoProcess(self, app=xlib.get_variant_calling_code())
        form_run_variant_calling_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_variant_calling_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_variant_calling_process'
        self.forms_dict[self.current_form] = form_run_variant_calling_process

        # raise "form_run_variant_calling_process" to front
        form_run_variant_calling_process.tkraise()

    #---------------

    def restart_variant_calling_process(self):
        '''
        Restart a Variant calling process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_variant_calling_process" in "container" with the grid geometry manager
        form_restart_variant_calling_process = gbioinfoapp.FormRestartBioinfoProcess(self, app=xlib.get_variant_calling_code())
        form_restart_variant_calling_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_variant_calling_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_variant_calling_process'
        self.forms_dict[self.current_form] = form_restart_variant_calling_process

        # raise "form_restart_variant_calling_process" to front
        form_restart_variant_calling_process.tkraise()

    #---------------
    # RAD-seq data files
    #---------------

    def recreate_end_file(self):
        '''
        Recreate the file of ends.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Recreate file of ends'

        # confirm the creation of the file of ends
        message = f'The file {xddradseqtools.get_end_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of ends
        if OK:
            (OK, error_list) = xddradseqtools.create_end_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of ends
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_end_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xddradseqtools.check_end_file(strict=False)
            if OK:
                message = f'The file {xddradseqtools.get_end_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_end_file(self):
        '''
        Edit the file of ends.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Edit file of ends'

        # edit the file of ends using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_end_file())
        self.root.wait_window(dialog_editor)

        # check the file of ends
        (OK, error_list) = xddradseqtools.check_end_file(strict=False)
        if OK:
            message = f'The file {xddradseqtools.get_end_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def recreate_individual_file(self):
        '''
        Recreate the file of individuals.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Recreate file of individuals'

        # confirm the creation of the file of individuals
        message = f'The file {xddradseqtools.get_individual_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of individuals
        if OK:
            (OK, error_list) = xddradseqtools.create_individual_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of individuals
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_individual_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xddradseqtools.check_individual_file(strict=False)
            if OK:
                message = f'The file {xddradseqtools.get_individual_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_individual_file(self):
        '''
        Edit the file of individuals.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Edit file of individuals'

        # edit the file of ends using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_individual_file())
        self.root.wait_window(dialog_editor)

        # check the file of individuals
        (OK, error_list) = xddradseqtools.check_individual_file(strict=False)
        if OK:
            message = f'The file {xddradseqtools.get_individual_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def recreate_raddesinger_condition_file(self):
        '''
        Recreate the file of RADdesigner conditions.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Recreate file of RADdesigner conditions'

        # confirm the creation of the file of RADdesigner conditions
        message = f'The file {xraddesigner.get_condition_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of RADdesigner conditions
        if OK:
            (OK, error_list) = xraddesigner.create_condition_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of RADdesigner conditions
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xraddesigner.get_condition_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xraddesigner.check_condition_file(strict=False)
            if OK:
                message = f'The file {xraddesigner.get_condition_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_raddesinger_condition_file(self):
        '''
        Edit the file of RADdesigner conditions.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Edit file of RADdesigner conditions'

        # edit the file of ends using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xraddesigner.get_condition_file())
        self.root.wait_window(dialog_editor)

        # check the file of RADdesigner conditions
        (OK, error_list) = xraddesigner.check_condition_file(strict=False)
        if OK:
            message = f'The file {xraddesigner.get_condition_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def recreate_restriction_site_file(self):
        '''
        Recreate the file of restriction sites.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Recreate file of restriction sites'

        # confirm the creation of the file of restriction sites
        message = f'The file {xddradseqtools.get_restriction_site_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of restriction sites
        if OK:
            (OK, error_list) = xddradseqtools.create_restriction_site_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of restriction sites
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_restriction_site_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xddradseqtools.check_restriction_site_file(strict=False)
            if OK:
                message = f'The file {xddradseqtools.get_restriction_site_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_restriction_site_file(self):
        '''
        Edit the file of restriction sites.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Edit file of restriction sites'

        # edit the file of restriction sites using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xddradseqtools.get_restriction_site_file())
        self.root.wait_window(dialog_editor)

        # check the file of restriction sites
        (OK, error_list) = xddradseqtools.check_restriction_site_file(strict=False)
        if OK:
            message = f'The file {xddradseqtools.get_restriction_site_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def recreate_vcf_sample_file(self):
        '''
        Recreate the file of VCF samples.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Recreate file of VCF samples'

        # confirm the creation of the file of VCF samples
        message = f'The file {xngshelper.get_vcf_sample_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of VCF samples
        if OK:
            (OK, error_list) = xngshelper.create_vcf_sample_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of VCF samples
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xngshelper.get_vcf_sample_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xngshelper.check_vcf_sample_file(strict=False)
            if OK:
                message = f'The file {xngshelper.get_vcf_sample_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_vcf_sample_file(self):
        '''
        Edit the file of VCF samples.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_ddradseqtools_name()} - Edit file of VCF samples'

        # edit the file of ends using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xngshelper.get_vcf_sample_file())
        self.root.wait_window(dialog_editor)

        # check the file of VCF samples
        (OK, error_list) = xngshelper.check_vcf_sample_file(strict=False)
        if OK:
            message = f'The file {xngshelper.get_vcf_sample_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------
    # TOA
    #---------------

    def recreate_toa_config_file(self):
        '''
        Recreate the TOA config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_toa_config_file" in "container" with the grid geometry manager
        form_recreate_toa_config_file = gtoa.FormRecreateToaConfigFile(self)
        form_recreate_toa_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_toa_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_create_toa_config_file'
        self.forms_dict[self.current_form] = form_recreate_toa_config_file

        # raise "form_recreate_toa_config_file" to front
        form_recreate_toa_config_file.tkraise()

    #---------------

    def view_toa_config_file(self):
        '''
        List the TOA config file.
        '''

        # close the current form
        self.close_current_form()

        # get the TOA config file
        toa_config_file = xtoa.get_toa_config_file()

        # create and show a instance DialogViewer to view the TOA config file
        dialog_viewer = gdialogs.DialogViewer(self.root, toa_config_file)
        self.root.wait_window(dialog_viewer)

    #---------------

    def recreate_toa_database(self):
        '''
        Recreate the TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_toa_database" in "container" with the grid geometry manager
        form_recreate_toa_database = gtoa.FormManageToaDatabase(self, process_type=xlib.get_toa_type_recreate())
        form_recreate_toa_database.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_toa_database" as current form and add it in the forms dictionary
        self.current_form = 'form_create_toa_config_file'
        self.forms_dict[self.current_form] = form_recreate_toa_database

        # raise "form_recreate_toa_database" to front
        form_recreate_toa_database.tkraise()

    #---------------

    def rebuild_toa_database(self):
        '''
        Rebuild the TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_rebuild_toa_database" in "container" with the grid geometry manager
        form_rebuild_toa_database = gtoa.FormManageToaDatabase(self, process_type=xlib.get_toa_type_rebuild())
        form_rebuild_toa_database.grid(row=0, column=0, sticky='nsew')

        # set "form_rebuild_toa_database" as current form and add it in the forms dictionary
        self.current_form = 'form_create_toa_config_file'
        self.forms_dict[self.current_form] = form_rebuild_toa_database

        # raise "form_rebuild_toa_database" to front
        form_rebuild_toa_database.tkraise()

    #---------------

    def recreate_dataset_file(self):
        '''
        Recreate the file of datasets.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_name()} - Recreate file of datasets'

        # confirm the creation of the file of datasets
        message = f'The file {xtoa.get_dataset_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of datasets
        if OK:
            (OK, error_list) = xtoa.create_dataset_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of datasets
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_dataset_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xtoa.check_dataset_file(strict=False)
            if OK:
                message = f'The file {xtoa.get_dataset_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_dataset_file(self):
        '''
        Edit the file of datasets.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_name()} - Edit file of datasets'

        # edit the file of datasets using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_dataset_file())
        self.root.wait_window(dialog_editor)

        # check the file of datasets
        (OK, error_list) = xtoa.check_dataset_file(strict=False)
        if OK:
            message = f'The file {xtoa.get_dataset_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def recreate_species_file(self):
        '''
        Recreate the file of species.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_name()} - Recreate file of species'

        # confirm the creation of the file of species
        message = f'The file {xtoa.get_species_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {head}', message)

        # recreate the file of species
        if OK:
            (OK, error_list) = xtoa.create_species_file()
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

        # edit the file of species
        if OK:

            # edit the data file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_species_file())
            self.root.wait_window(dialog_editor)

            # check the data file
            (OK, error_list) = xtoa.check_species_file(strict=False)
            if OK:
                message = f'The file {xtoa.get_species_file()} is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def edit_species_file(self):
        '''
        Edit the file of species.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_name()} - Edit file of species'

        # edit the file of species using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_species_file())
        self.root.wait_window(dialog_editor)

        # check the file of species
        (OK, error_list) = xtoa.check_species_file(strict=False)
        if OK:
            message = f'The file {xtoa.get_species_file()} is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def download_basic_data(self):
        '''
        Download other basic data.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_basic_data_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_basic_data(self):
        '''
        Load basic data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_basic_data_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def build_gymno_01_proteome(self):
        '''
        Build the Gymno PLAZA 1.0 proteome.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_proteome(), genomic_database=xlib.get_toa_data_gymno_01_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_gymno_01_data(self):
        '''
        Download Gymno PLAZA 1.0 functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_gymno_01_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_gymno_01_data(self):
        '''
        Load Gymno PLAZA 1.0 data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_gymno_01_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def build_dicots_04_proteome(self):
        '''
        Build the Dicots PLAZA 4.0 proteome.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_proteome(), genomic_database=xlib.get_toa_data_dicots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_dicots_04_data(self):
        '''
        Download Dicots PLAZA 4.0 functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_dicots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_dicots_04_data(self):
        '''
        Load Dicots PLAZA 4.0 data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_dicots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def build_monocots_04_proteome(self):
        '''
        Build the Monocots PLAZA 4.0 proteome.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_proteome(), genomic_database=xlib.get_toa_data_monocots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_monocots_04_data(self):
        '''
        Download Monocots PLAZA 4.0 functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_monocots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_monocots_04_data(self):
        '''
        Load Monocots PLAZA 4.0 data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_monocots_04_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def build_refseq_plant_proteome(self):
        '''
        Build the NCBI RefSeq Plant proteome.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_proteome(), genomic_database=xlib.get_toa_data_refseq_plant_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_taxonomy_data(self):
        '''
        Download NCBI Taxonomy data.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_download_taxonomy_data" in "container" with the grid geometry manager
        form_download_taxonomy_data = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_taxonomy_code())
        form_download_taxonomy_data.grid(row=0, column=0, sticky='nsew')

        # set "form_download_taxonomy_data" as current form and add it in the forms dictionary
        self.current_form = 'form_download_taxonomy_data'
        self.forms_dict[self.current_form] = form_download_taxonomy_data

        # raise "form_download_taxonomy_data" to front
        form_download_taxonomy_data.tkraise()

    #---------------

    def build_blastplus_nt_db(self):
        '''
        Build the NCBI BLAST database NT for BLAST+.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_build_blastplus_nt_db" in "container" with the grid geometry manager
        form_build_blastplus_nt_db = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_blastplus_db(), genomic_database=xlib.get_toa_data_nt_code())
        form_build_blastplus_nt_db.grid(row=0, column=0, sticky='nsew')

        # set "form_build_blastplus_nt_db" as current form and add it in the forms dictionary
        self.current_form = 'form_build_blastplus_nt_db'
        self.forms_dict[self.current_form] = form_build_blastplus_nt_db

        # raise "form_build_blastplus_nt_db" to front
        form_build_blastplus_nt_db.tkraise()

    #---------------

    def build_viridiplantae_nucleotide_gi_gilist(self):
        '''
        Build the NCBI Nucleotide GenInfo viridiplantae identifier list.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_gilist(), genomic_database=xlib.get_toa_data_viridiplantae_nucleotide_gi_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def build_blastplus_nr_db(self):
        '''
        Build the NCBI BLAST database NR for BLAST+.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_build_blastplus_nr_db" in "container" with the grid geometry manager
        form_build_blastplus_nr_db = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_blastplus_db(), genomic_database=xlib.get_toa_data_nr_code())
        form_build_blastplus_nr_db.grid(row=0, column=0, sticky='nsew')

        # set "form_build_blastplus_nr_db" as current form and add it in the forms dictionary
        self.current_form = 'form_build_blastplus_nr_db'
        self.forms_dict[self.current_form] = form_build_blastplus_nr_db

        # raise "form_build_blastplus_nr_db" to front
        form_build_blastplus_nr_db.tkraise()

    #---------------

    def build_diamond_nr_db(self):
        '''
        Build the NCBI BLAST database NR for DIAMOND.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_build_diamond_nr_db" in "container" with the grid geometry manager
        form_build_diamond_nr_db = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_diamond_db(), genomic_database=xlib.get_toa_data_nr_code())
        form_build_diamond_nr_db.grid(row=0, column=0, sticky='nsew')

        # set "form_build_diamond_nr_db" as current form and add it in the forms dictionary
        self.current_form = 'form_build_diamond_nr_db'
        self.forms_dict[self.current_form] = form_build_diamond_nr_db

        # raise "form_build_diamond_nr_db" to front
        form_build_diamond_nr_db.tkraise()

    #---------------

    def build_viridiplantae_protein_gi_gilist(self):
        '''
        Build the NCBI Protein GenInfo viridiplantae identifier list.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_build_gilist(), genomic_database=xlib.get_toa_data_viridiplantae_protein_gi_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_gene_data(self):
        '''
        Download NCBI Gene functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_gene_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_gene_data(self):
        '''
        Load NCBI Gene data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_gene_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_interpro_data(self):
        '''
        Download InterPro functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_interpro_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_interpro_data(self):
        '''
        Load InterPro data into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_interpro_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def download_go_data(self):
        '''
        Download Gene Ontology functional annotation.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_download_data(), genomic_database=xlib.get_toa_data_go_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def load_go_data(self):
        '''
        Load NCBI Gene Ontology into TOA database.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_manage_genomic_database" in "container" with the grid geometry manager
        form_manage_genomic_database = gtoa.FormManageGenomicDatabase(self, process_type=xlib.get_toa_type_load_data(), genomic_database=xlib.get_toa_data_go_code())
        form_manage_genomic_database.grid(row=0, column=0, sticky='nsew')

        # set "form_manage_genomic_database" as current form and add it in the forms dictionary
        self.current_form = 'form_manage_genomic_database'
        self.forms_dict[self.current_form] = form_manage_genomic_database

        # raise "form_manage_genomic_database" to front
        form_manage_genomic_database.tkraise()

    #---------------

    def recreate_nucleotide_pipeline_config_file(self):
        '''
        Recreate the nucleotide pipeline config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_nucleotide_pipeline_config_file" in "container" with the grid geometry manager
        form_recreate_nucleotide_pipeline_config_file = gtoa.FormRecreatePipelineConfigFile(self, pipeline_type=xlib.get_toa_process_pipeline_nucleotide_code())
        form_recreate_nucleotide_pipeline_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_nucleotide_pipeline_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_nucleotide_pipeline_config_file'
        self.forms_dict[self.current_form] = form_recreate_nucleotide_pipeline_config_file

        # raise "form_recreate_nucleotide_pipeline_config_file" to front
        form_recreate_nucleotide_pipeline_config_file.tkraise()

    #---------------

    def edit_nucleotide_pipeline_config_file(self):
        '''
        Edit the nucleotide pipeline config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_process_pipeline_nucleotide_name()} - Edit config file'

        # edit the nucleotide pipeline config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_nucleotide_pipeline_config_file())
        self.root.wait_window(dialog_editor)

        # check the nucleotide pipeline config file
        (OK, error_list) = xtoa.check_pipeline_config_file(pipeline_type=xlib.get_toa_process_pipeline_nucleotide_code(), strict=False)
        if OK:
            message = f'The {xlib.get_toa_process_pipeline_nucleotide_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_nucleotide_pipeline_process(self):
        '''
        Run a nucleotide pipeline process corresponding to the options in nucleotide pipeline config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_nucleotide_pipeline_process" in "container" with the grid geometry manager
        form_run_nucleotide_pipeline_process = gtoa.FormRunPipelineProcess(self, pipeline_type=xlib.get_toa_process_pipeline_nucleotide_code())
        form_run_nucleotide_pipeline_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_nucleotide_pipeline_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_nucleotide_pipeline_process'
        self.forms_dict[self.current_form] = form_run_nucleotide_pipeline_process

        # raise "form_run_nucleotide_pipeline_process" to front
        form_run_nucleotide_pipeline_process.tkraise()

    #---------------

    def restart_nucleotide_pipeline_process(self):
        '''
        Restart a nucleotide pipeline process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_nucleotide_pipeline_process" in "container" with the grid geometry manager
        form_restart_nucleotide_pipeline_process = gtoa.FormRestartPipelineProcess(self, pipeline_type=xlib.get_toa_process_pipeline_nucleotide_code())
        form_restart_nucleotide_pipeline_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_nucleotide_pipeline_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_nucleotide_pipeline_process'
        self.forms_dict[self.current_form] = form_restart_nucleotide_pipeline_process

        # raise "form_restart_nucleotide_pipeline_process" to front
        form_restart_nucleotide_pipeline_process.tkraise()

    #---------------

    def recreate_aminoacid_pipeline_config_file(self):
        '''
        Recreate the amino acid pipeline config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_aminoacid_pipeline_config_file" in "container" with the grid geometry manager
        form_recreate_aminoacid_pipeline_config_file = gtoa.FormRecreatePipelineConfigFile(self, pipeline_type=xlib.get_toa_process_pipeline_aminoacid_code())
        form_recreate_aminoacid_pipeline_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_aminoacid_pipeline_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_aminoacid_pipeline_config_file'
        self.forms_dict[self.current_form] = form_recreate_aminoacid_pipeline_config_file

        # raise "form_recreate_aminoacid_pipeline_config_file" to front
        form_recreate_aminoacid_pipeline_config_file.tkraise()

    #---------------

    def edit_aminoacid_pipeline_config_file(self):
        '''
        Edit the amino acid pipeline config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_process_pipeline_aminoacid_name()} - Edit config file'

        # edit the amino acid pipeline config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_aminoacid_pipeline_config_file())
        self.root.wait_window(dialog_editor)

        # check the amino acid pipeline config file
        (OK, error_list) = xtoa.check_pipeline_config_file(pipeline_type=xlib.get_toa_process_pipeline_aminoacid_code(), strict=False)
        if OK:
            message = f'The {xlib.get_toa_process_pipeline_aminoacid_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_aminoacid_pipeline_process(self):
        '''
        Run a amino acid pipeline process corresponding to the options in amino acid pipeline config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_aminoacid_pipeline_process" in "container" with the grid geometry manager
        form_run_aminoacid_pipeline_process = gtoa.FormRunPipelineProcess(self, pipeline_type=xlib.get_toa_process_pipeline_aminoacid_code())
        form_run_aminoacid_pipeline_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_aminoacid_pipeline_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_aminoacid_pipeline_process'
        self.forms_dict[self.current_form] = form_run_aminoacid_pipeline_process

        # raise "form_run_aminoacid_pipeline_process" to front
        form_run_aminoacid_pipeline_process.tkraise()

    #---------------

    def restart_aminoacid_pipeline_process(self):
        '''
        Restart a amino acid pipeline process from the last step ended OK.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_restart_aminoacid_pipeline_process" in "container" with the grid geometry manager
        form_restart_aminoacid_pipeline_process = gtoa.FormRestartPipelineProcess(self, pipeline_type=xlib.get_toa_process_pipeline_aminoacid_code())
        form_restart_aminoacid_pipeline_process.grid(row=0, column=0, sticky='nsew')

        # set "form_restart_aminoacid_pipeline_process" as current form and add it in the forms dictionary
        self.current_form = 'form_restart_aminoacid_pipeline_process'
        self.forms_dict[self.current_form] = form_restart_aminoacid_pipeline_process

        # raise "form_restart_aminoacid_pipeline_process" to front
        form_restart_aminoacid_pipeline_process.tkraise()

    #---------------

    def recreate_annotation_merger_config_file(self):
        '''
        Recreate the pipeline merger config file with the default options. It is necessary
        update the options in each process run.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_annotation_merger_config_file" in "container" with the grid geometry manager
        form_recreate_annotation_merger_config_file = gtoa.FormRecreateAnnotatioMergerConfigFile(self)
        form_recreate_annotation_merger_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_annotation_merger_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_annotation_merger_config_file'
        self.forms_dict[self.current_form] = form_recreate_annotation_merger_config_file

        # raise "form_recreate_annotation_merger_config_file" to front
        form_recreate_annotation_merger_config_file.tkraise()

    #---------------

    def edit_annotation_merger_config_file(self):
        '''
        Edit the pipeline merger config file to change the parameters of process run.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = f'{xlib.get_toa_process_merge_annotations_name()} - Edit config file'

        # edit the pipeline merger config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xtoa.get_annotation_merger_config_file())
        self.root.wait_window(dialog_editor)

        # check the nucleotide pipeline config file
        (OK, error_list) = xtoa.check_annotation_merger_config_file(strict=False)
        if OK:
            message = f'The {xlib.get_toa_process_merge_annotations_name()} config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n' 
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_annotation_merger_process(self):
        '''
        Run an annotation merger process corresponding to the options in nucleotide pipeline config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_annotation_merger_process" in "container" with the grid geometry manager
        form_run_annotation_merger_process = gtoa.FormRunPipelineProcess(self, pipeline_type=xlib.get_toa_process_merge_annotations_code())
        form_run_annotation_merger_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_annotation_merger_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_annotation_merger_process'
        self.forms_dict[self.current_form] = form_run_annotation_merger_process

        # raise "form_run_annotation_merger_process" to front
        form_run_annotation_merger_process.tkraise()

    #---------------

    def view_hit_per_hsp_data(self):
        '''
        View the # HITs per # HSPs data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_hit_per_hsp_data" in "container" with the grid geometry manager
        form_view_hit_per_hsp_data = gtoa.FormViewStats(self, stats_code='hit_per_hsp')
        form_view_hit_per_hsp_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_hit_per_hsp_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_hit_per_hsp_data'
        self.forms_dict[self.current_form] = form_view_hit_per_hsp_data

        # raise "form_view_hit_per_hsp_data" to front
        form_view_hit_per_hsp_data.tkraise()

    #---------------

    def plot_hit_per_hsp_data(self):
        '''
        Plot the # HITs per # HSPs data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_dataset_frequency" in "container" with the grid geometry manager
        form_plot_dataset_frequency = gtoa.FormPlotStats(self, stats_code='hit_per_hsp')
        form_plot_dataset_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_dataset_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_dataset_frequency'
        self.forms_dict[self.current_form] = form_plot_dataset_frequency

        # raise "form_plot_dataset_frequency" to front
        form_plot_dataset_frequency.tkraise()

    #---------------

    def view_annotation_dataset_frequency(self):
        '''
        View the dataset frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_dataset_frequency" in "container" with the grid geometry manager
        form_view_dataset_frequency = gtoa.FormViewStats(self, stats_code='dataset')
        form_view_dataset_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_dataset_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_dataset_frequency'
        self.forms_dict[self.current_form] = form_view_dataset_frequency

        # raise "form_view_dataset_frequency" to front
        form_view_dataset_frequency.tkraise()

    #---------------

    def plot_annotation_dataset_frequency(self):
        '''
        Plot the dataset frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_dataset_frequency" in "container" with the grid geometry manager
        form_plot_dataset_frequency = gtoa.FormPlotStats(self, stats_code='dataset')
        form_plot_dataset_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_dataset_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_dataset_frequency'
        self.forms_dict[self.current_form] = form_plot_dataset_frequency

        # raise "form_plot_dataset_frequency" to front
        form_plot_dataset_frequency.tkraise()

    #---------------

    def view_species_frequency(self):
        '''
        View the species frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_species_frequency" in "container" with the grid geometry manager
        form_view_species_frequency = gtoa.FormViewStats(self, stats_code='species')
        form_view_species_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_species_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_species_frequency'
        self.forms_dict[self.current_form] = form_view_species_frequency

        # raise "form_view_species_frequency" to front
        form_view_species_frequency.tkraise()

    #---------------

    def plot_species_frequency(self):
        '''
        Plot the species frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_species_frequency" in "container" with the grid geometry manager
        form_plot_species_frequency = gtoa.FormPlotStats(self, stats_code='species')
        form_plot_species_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_species_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_species_frequency'
        self.forms_dict[self.current_form] = form_plot_species_frequency

        # raise "form_plot_species_frequency" to front
        form_plot_species_frequency.tkraise()

    #---------------

    def view_family_frequency(self):
        '''
        View the family frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_family_frequency" in "container" with the grid geometry manager
        form_view_family_frequency = gtoa.FormViewStats(self, stats_code='family')
        form_view_family_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_family_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_family_frequency'
        self.forms_dict[self.current_form] = form_view_family_frequency

        # raise "form_view_family_frequency" to front
        form_view_family_frequency.tkraise()

    #---------------

    def plot_family_frequency(self):
        '''
        Plot the family frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_family_frequency" in "container" with the grid geometry manager
        form_plot_family_frequency = gtoa.FormPlotStats(self, stats_code='family')
        form_plot_family_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_family_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_family_frequency'
        self.forms_dict[self.current_form] = form_plot_family_frequency

        # raise "form_plot_family_frequency" to front
        form_plot_family_frequency.tkraise()

    #---------------

    def view_phylum_frequency(self):
        '''
        View the phylum frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_phylum_frequency" in "container" with the grid geometry manager
        form_view_phylum_frequency = gtoa.FormViewStats(self, stats_code='phylum')
        form_view_phylum_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_phylum_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_phylum_frequency'
        self.forms_dict[self.current_form] = form_view_phylum_frequency

        # raise "form_view_phylum_frequency" to front
        form_view_phylum_frequency.tkraise()

    #---------------

    def plot_phylum_frequency(self):
        '''
        Plot the phylum frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_phylum_frequency" in "container" with the grid geometry manager
        form_plot_phylum_frequency = gtoa.FormPlotStats(self, stats_code='phylum')
        form_plot_phylum_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_phylum_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_phylum_frequency'
        self.forms_dict[self.current_form] = form_plot_phylum_frequency

        # raise "form_plot_phylum_frequency" to front
        form_plot_phylum_frequency.tkraise()

    #---------------

    def view_ec_frequency(self):
        '''
        View the EC frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_ec_frequency" in "container" with the grid geometry manager
        form_view_ec_frequency = gtoa.FormViewStats(self, stats_code='ec')
        form_view_ec_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_ec_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_ec_frequency'
        self.forms_dict[self.current_form] = form_view_ec_frequency

        # raise "form_view_ec_frequency" to front
        form_view_ec_frequency.tkraise()

    #---------------

    def plot_ec_frequency(self):
        '''
        Plot the EC frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_ec_frequency" in "container" with the grid geometry manager
        form_plot_ec_frequency = gtoa.FormPlotStats(self, stats_code='ec')
        form_plot_ec_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_ec_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_ec_frequency'
        self.forms_dict[self.current_form] = form_plot_ec_frequency

        # raise "form_plot_ec_frequency" to front
        form_plot_ec_frequency.tkraise()

    #---------------

    def view_seq_per_ec_data(self):
        '''
        View the # sequences per # EC ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_ec_data" in "container" with the grid geometry manager
        form_view_seq_per_ec_data = gtoa.FormViewStats(self, stats_code='seq_per_ec')
        form_view_seq_per_ec_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_ec_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_ec_data'
        self.forms_dict[self.current_form] = form_view_seq_per_ec_data

        # raise "form_view_seq_per_ec_data" to front
        form_view_seq_per_ec_data.tkraise()

    #---------------

    def plot_seq_per_ec_data(self):
        '''
        Plot the # sequences per # EC ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_ec_data" in "container" with the grid geometry manager
        form_plot_seq_per_ec_data = gtoa.FormPlotStats(self, stats_code='seq_per_ec')
        form_plot_seq_per_ec_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_ec_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_ec_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_ec_data

        # raise "form_plot_seq_per_ec_data" to front
        form_plot_seq_per_ec_data.tkraise()

    #---------------

    def view_go_frequency(self):
        '''
        View the GO frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_go_frequency" in "container" with the grid geometry manager
        form_view_go_frequency = gtoa.FormViewStats(self, stats_code='go')
        form_view_go_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_go_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_go_frequency'
        self.forms_dict[self.current_form] = form_view_go_frequency

        # raise "form_view_go_frequency" to front
        form_view_go_frequency.tkraise()

    #---------------

    def plot_go_frequency(self):
        '''
        Plot the GO frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_go_frequency" in "container" with the grid geometry manager
        form_plot_go_frequency = gtoa.FormPlotStats(self, stats_code='go')
        form_plot_go_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_go_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_go_frequency'
        self.forms_dict[self.current_form] = form_plot_go_frequency

        # raise "form_plot_go_frequency" to front
        form_plot_go_frequency.tkraise()

    #---------------

    def view_namespace_frequency(self):
        '''
        View the namespace frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_namespace_frequency" in "container" with the grid geometry manager
        form_view_namespace_frequency = gtoa.FormViewStats(self, stats_code='namespace')
        form_view_namespace_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_namespace_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_namespace_frequency'
        self.forms_dict[self.current_form] = form_view_namespace_frequency

        # raise "form_view_namespace_frequency" to front
        form_view_namespace_frequency.tkraise()

    #---------------

    def plot_namespace_frequency(self):
        '''
        Plot the GO frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_namespace_frequency" in "container" with the grid geometry manager
        form_plot_namespace_frequency = gtoa.FormPlotStats(self, stats_code='namespace')
        form_plot_namespace_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_namespace_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_namespace_frequency'
        self.forms_dict[self.current_form] = form_plot_namespace_frequency

        # raise "form_plot_namespace_frequency" to front
        form_plot_namespace_frequency.tkraise()

    #---------------

    def view_seq_per_go_data(self):
        '''
        View the # sequences per GO terms data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_go_data" in "container" with the grid geometry manager
        form_view_seq_per_go_data = gtoa.FormViewStats(self, stats_code='seq_per_go')
        form_view_seq_per_go_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_go_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_go_data'
        self.forms_dict[self.current_form] = form_view_seq_per_go_data

        # raise "form_view_seq_per_go_data" to front
        form_view_seq_per_go_data.tkraise()

    #---------------

    def plot_seq_per_go_data(self):
        '''
        Plot the # sequences per GO terms data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_go_data" in "container" with the grid geometry manager
        form_plot_seq_per_go_data = gtoa.FormPlotStats(self, stats_code='seq_per_go')
        form_plot_seq_per_go_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_go_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_go_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_go_data

        # raise "form_plot_seq_per_go_data" to front
        form_plot_seq_per_go_data.tkraise()

    #---------------

    def view_interpro_frequency(self):
        '''
        View the InterPro frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_interpro_frequency" in "container" with the grid geometry manager
        form_view_interpro_frequency = gtoa.FormViewStats(self, stats_code='interpro')
        form_view_interpro_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_interpro_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_interpro_frequency'
        self.forms_dict[self.current_form] = form_view_interpro_frequency

        # raise "form_view_interpro_frequency" to front
        form_view_interpro_frequency.tkraise()

    #---------------

    def plot_interpro_frequency(self):
        '''
        Plot the InterPro frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_interpro_frequency" in "container" with the grid geometry manager
        form_plot_interpro_frequency = gtoa.FormPlotStats(self, stats_code='interpro')
        form_plot_interpro_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_interpro_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_interpro_frequency'
        self.forms_dict[self.current_form] = form_plot_interpro_frequency

        # raise "form_plot_interpro_frequency" to front
        form_plot_interpro_frequency.tkraise()

    #---------------

    def view_seq_per_interpro_data(self):
        '''
        View the # sequences per # InterPro ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_interpro_data" in "container" with the grid geometry manager
        form_view_seq_per_interpro_data = gtoa.FormViewStats(self, stats_code='seq_per_interpro')
        form_view_seq_per_interpro_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_interpro_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_interpro_data'
        self.forms_dict[self.current_form] = form_view_seq_per_interpro_data

        # raise "form_view_seq_per_interpro_data" to front
        form_view_seq_per_interpro_data.tkraise()

    #---------------

    def plot_seq_per_interpro_data(self):
        '''
        Plot the # sequences per # InterPro ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_interpro_data" in "container" with the grid geometry manager
        form_plot_seq_per_interpro_data = gtoa.FormPlotStats(self, stats_code='seq_per_interpro')
        form_plot_seq_per_interpro_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_interpro_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_interpro_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_interpro_data

        # raise "form_plot_seq_per_interpro_data" to front
        form_plot_seq_per_interpro_data.tkraise()

    #---------------

    def view_kegg_frequency(self):
        '''
        View the KEGG frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_kegg_frequency" in "container" with the grid geometry manager
        form_view_kegg_frequency = gtoa.FormViewStats(self, stats_code='kegg')
        form_view_kegg_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_kegg_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_kegg_frequency'
        self.forms_dict[self.current_form] = form_view_kegg_frequency

        # raise "form_view_kegg_frequency" to front
        form_view_kegg_frequency.tkraise()

    #---------------

    def plot_kegg_frequency(self):
        '''
        Plot the KEGG frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_kegg_frequency" in "container" with the grid geometry manager
        form_plot_kegg_frequency = gtoa.FormPlotStats(self, stats_code='kegg')
        form_plot_kegg_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_kegg_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_kegg_frequency'
        self.forms_dict[self.current_form] = form_plot_kegg_frequency

        # raise "form_plot_kegg_frequency" to front
        form_plot_kegg_frequency.tkraise()

    #---------------

    def view_seq_per_kegg_data(self):
        '''
        View the # sequences per # KEGG ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_kegg_data" in "container" with the grid geometry manager
        form_view_seq_per_kegg_data = gtoa.FormViewStats(self, stats_code='seq_per_kegg')
        form_view_seq_per_kegg_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_kegg_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_kegg_data'
        self.forms_dict[self.current_form] = form_view_seq_per_kegg_data

        # raise "form_view_seq_per_kegg_data" to front
        form_view_seq_per_kegg_data.tkraise()

    #---------------

    def plot_seq_per_kegg_data(self):
        '''
        Plot the # sequences per # KEGG ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_kegg_data" in "container" with the grid geometry manager
        form_plot_seq_per_kegg_data = gtoa.FormPlotStats(self, stats_code='seq_per_kegg')
        form_plot_seq_per_kegg_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_kegg_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_kegg_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_kegg_data

        # raise "form_plot_seq_per_kegg_data" to front
        form_plot_seq_per_kegg_data.tkraise()

    #---------------

    def view_mapman_frequency(self):
        '''
        View the MapMan frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_mapman_frequency" in "container" with the grid geometry manager
        form_view_mapman_frequency = gtoa.FormViewStats(self, stats_code='mapman')
        form_view_mapman_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_mapman_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_mapman_frequency'
        self.forms_dict[self.current_form] = form_view_mapman_frequency

        # raise "form_view_mapman_frequency" to front
        form_view_mapman_frequency.tkraise()

    #---------------

    def plot_mapman_frequency(self):
        '''
        Plot the MapMan frequency distribution of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_mapman_frequency" in "container" with the grid geometry manager
        form_plot_mapman_frequency = gtoa.FormPlotStats(self, stats_code='mapman')
        form_plot_mapman_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_mapman_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_mapman_frequency'
        self.forms_dict[self.current_form] = form_plot_mapman_frequency

        # raise "form_plot_mapman_frequency" to front
        form_plot_mapman_frequency.tkraise()

    #---------------

    def view_seq_per_mapman_data(self):
        '''
        View the # sequences per # MapMan ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_mapman_data" in "container" with the grid geometry manager
        form_view_seq_per_mapman_data = gtoa.FormViewStats(self, stats_code='seq_per_mapman')
        form_view_seq_per_mapman_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_mapman_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_mapman_data'
        self.forms_dict[self.current_form] = form_view_seq_per_mapman_data

        # raise "form_view_seq_per_mapman_data" to front
        form_view_seq_per_mapman_data.tkraise()

    #---------------

    def plot_seq_per_mapman_data(self):
        '''
        Plot the # sequences per # MapMan ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_mapman_data" in "container" with the grid geometry manager
        form_plot_seq_per_mapman_data = gtoa.FormPlotStats(self, stats_code='seq_per_mapman')
        form_plot_seq_per_mapman_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_mapman_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_mapman_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_mapman_data

        # raise "form_plot_seq_per_mapman_data" to front
        form_plot_seq_per_mapman_data.tkraise()

    #---------------

    def view_metacyc_frequency(self):
        '''
        View the MetaCyc frequency of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_metacyc_frequency" in "container" with the grid geometry manager
        form_view_metacyc_frequency = gtoa.FormViewStats(self, stats_code='metacyc')
        form_view_metacyc_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_view_metacyc_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_view_metacyc_frequency'
        self.forms_dict[self.current_form] = form_view_metacyc_frequency

        # raise "form_view_metacyc_frequency" to front
        form_view_metacyc_frequency.tkraise()

    #---------------

    def plot_metacyc_frequency(self):
        '''
        Plot the MetaCyc frequency of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_metacyc_frequency" in "container" with the grid geometry manager
        form_plot_metacyc_frequency = gtoa.FormPlotStats(self, stats_code='metacyc')
        form_plot_metacyc_frequency.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_metacyc_frequency" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_metacyc_frequency'
        self.forms_dict[self.current_form] = form_plot_metacyc_frequency

        # raise "form_plot_metacyc_frequency" to front
        form_plot_metacyc_frequency.tkraise()

    #---------------

    def view_seq_per_metacyc_data(self):
        '''
        View the # sequences per # MetaCyc ids data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_seq_per_metacyc_data" in "container" with the grid geometry manager
        form_view_seq_per_metacyc_data = gtoa.FormViewStats(self, stats_code='seq_per_metacyc')
        form_view_seq_per_metacyc_data.grid(row=0, column=0, sticky='nsew')

        # set "form_view_seq_per_metacyc_data" as current form and add it in the forms dictionary
        self.current_form = 'form_view_seq_per_metacyc_data'
        self.forms_dict[self.current_form] = form_view_seq_per_metacyc_data

        # raise "form_view_seq_per_metacyc_data" to front
        form_view_seq_per_metacyc_data.tkraise()

    #---------------

    def plot_seq_per_metacyc_data(self):
        '''
        Plot the # sequences per # MetaCyc data of an annotation pipeline.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_plot_seq_per_metacyc_data" in "container" with the grid geometry manager
        form_plot_seq_per_metacyc_data = gtoa.FormPlotStats(self, stats_code='seq_per_metacyc')
        form_plot_seq_per_metacyc_data.grid(row=0, column=0, sticky='nsew')

        # set "form_plot_seq_per_metacyc_data" as current form and add it in the forms dictionary
        self.current_form = 'form_plot_seq_per_metacyc_data'
        self.forms_dict[self.current_form] = form_plot_seq_per_metacyc_data

        # raise "form_plot_seq_per_metacyc_data" to front
        form_plot_seq_per_metacyc_data.tkraise()

    #---------------
    # Datasets
    #---------------

    def list_dataset(self):
        '''
        List datasets showing data of its directories and files.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_list_dataset" in "container" with the grid geometry manager
        form_list_dataset = gdataset.FormListDataset(self)
        form_list_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_list_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_list_dataset'
        self.forms_dict[self.current_form] = form_list_dataset

        # raise "form_list_dataset" to front
        form_list_dataset.tkraise()

    #---------------

    def recreate_reference_transfer_config_file(self):
        '''
        Recreate the reference transfer config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_reference_transfer_config_file" in "container" with the grid geometry manager
        form_recreate_reference_transfer_config_file = gdataset.FormRecreateReferenceTransferConfigFile(self)
        form_recreate_reference_transfer_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_reference_transfer_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_reference_transfer_config_file'
        self.forms_dict[self.current_form] = form_recreate_reference_transfer_config_file

        # raise "form_recreate_reference_transfer_config_file" to front
        form_recreate_reference_transfer_config_file.tkraise()

    #---------------

    def edit_reference_transfer_config_file(self):
        '''
        Edit the reference transfer config file to change the parameters of each transfer.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Reference dataset file transfer - Edit config file'

        # edit the reference transfer config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xreference.get_reference_transfer_config_file())
        self.root.wait_window(dialog_editor)

        # check the reference transfer config file
        (OK, error_list) = xreference.check_reference_transfer_config_file(strict=False)
        if OK:
            message = 'The reference transfer config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def upload_reference_dataset(self):
        '''
       Upload a reference dataset to a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_upload_reference_dataset" in "container" with the grid geometry manager
        form_upload_reference_dataset = gdataset.FormUploadReferenceDataSet(self)
        form_upload_reference_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_upload_reference_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_upload_reference_dataset'
        self.forms_dict[self.current_form] = form_upload_reference_dataset

        # raise "form_upload_reference_dataset" to front
        form_upload_reference_dataset.tkraise()

    #---------------

    def recreate_reference_gzip_config_file(self):
        '''
       Recreate the reference file compression/decompression config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_reference_gzip_config_file" in "container" with the grid geometry manager
        form_recreate_reference_gzip_config_file = gdataset.FormRecreateReferenceGzipConfigFile(self)
        form_recreate_reference_gzip_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_reference_gzip_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_reference_gzip_config_file'
        self.forms_dict[self.current_form] = form_recreate_reference_gzip_config_file

        # raise "form_recreate_reference_gzip_config_file" to front
        form_recreate_reference_gzip_config_file.tkraise()

    #---------------

    def edit_reference_gzip_config_file(self):
        '''
       Edit the reference file compression/decompression config file.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Reference dataset file compression/decompression - Edit config file'

        # edit the reference gzip config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgzip.get_gzip_config_file('reference'))
        self.root.wait_window(dialog_editor)

        # check the reference gzip config file
        (OK, error_list) = xgzip.check_gzip_config_file('reference', strict=False)
        if OK:
            message = 'The reference gzip config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_reference_gzip_process(self):
        '''
        Compress/decompress reference dataset files in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_reference_gzip_process" in "container" with the grid geometry manager
        form_run_reference_gzip_process = gdataset.FormRunGzipProcess(self, 'reference')
        form_run_reference_gzip_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_reference_gzip_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_reference_gzip_process'
        self.forms_dict[self.current_form] = form_run_reference_gzip_process

        # raise "form_run_reference_gzip_process" to front
        form_run_reference_gzip_process.tkraise()

    #---------------

    def remove_reference_dataset(self):
        '''
       Remove the reference dataset in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_reference_dataset" in "container" with the grid geometry manager
        form_remove_reference_dataset = gdataset.FormRemoveReferenceDataSet(self)
        form_remove_reference_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_reference_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_reference_dataset'
        self.forms_dict[self.current_form] = form_remove_reference_dataset

        # raise "form_remove_reference_dataset" to front
        form_remove_reference_dataset.tkraise()

    #---------------

    def recreate_database_transfer_config_file(self):
        '''
        Recreate the database transfer config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_database_transfer_config_file" in "container" with the grid geometry manager
        form_recreate_database_transfer_config_file = gdataset.FormRecreateDatabaseTransferConfigFile(self)
        form_recreate_database_transfer_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_database_transfer_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_database_transfer_config_file'
        self.forms_dict[self.current_form] = form_recreate_database_transfer_config_file

        # raise "form_recreate_database_transfer_config_file" to front
        form_recreate_database_transfer_config_file.tkraise()

    #---------------

    def edit_database_transfer_config_file(self):
        '''
        Edit the database transfer config file to change the parameters of each transfer.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Database file transfer - Edit config file'

        # edit the database transfer config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xdatabase.get_database_transfer_config_file())
        self.root.wait_window(dialog_editor)

        # check the database transfer config file
        (OK, error_list) = xdatabase.check_database_transfer_config_file(strict=False)
        if OK:
            message = 'The database transfer config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def upload_database_dataset(self):
        '''
       Upload a database dataset to a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_upload_database_dataset" in "container" with the grid geometry manager
        form_upload_database_dataset = gdataset.FormUploadDatabaseDataSet(self)
        form_upload_database_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_upload_database_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_upload_database_dataset'
        self.forms_dict[self.current_form] = form_upload_database_dataset

        # raise "form_upload_database_dataset" to front
        form_upload_database_dataset.tkraise()

    #---------------

    def recreate_database_gzip_config_file(self):
        '''
       Recreate the database file compression/decompression config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_database_gzip_config_file" in "container" with the grid geometry manager
        form_recreate_database_gzip_config_file = gdataset.FormRecreateDatabaseGzipConfigFile(self)
        form_recreate_database_gzip_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_database_gzip_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_database_gzip_config_file'
        self.forms_dict[self.current_form] = form_recreate_database_gzip_config_file

        # raise "form_recreate_database_gzip_config_file" to front
        form_recreate_database_gzip_config_file.tkraise()

    #---------------

    def edit_database_gzip_config_file(self):
        '''
       Edit the database file compression/decompression config file.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Database file compression/decompression - Edit config file'

        # edit the database gzip config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgzip.get_gzip_config_file('database'))
        self.root.wait_window(dialog_editor)

        # check the database gzip config file
        (OK, error_list) = xgzip.check_gzip_config_file('database', strict=False)
        if OK:
            message = 'The database gzip config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_database_gzip_process(self):
        '''
        Compress/decompress database dataset files in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_database_gzip_process" in "container" with the grid geometry manager
        form_run_database_gzip_process = gdataset.FormRunGzipProcess(self, 'database')
        form_run_database_gzip_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_database_gzip_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_database_gzip_process'
        self.forms_dict[self.current_form] = form_run_database_gzip_process

        # raise "form_run_database_gzip_process" to front
        form_run_database_gzip_process.tkraise()

    #---------------

    def remove_database_dataset(self):
        '''
       Remove the database dataset in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_database_dataset" in "container" with the grid geometry manager
        form_remove_database_dataset = gdataset.FormRemoveDatabaseDataSet(self)
        form_remove_database_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_database_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_database_dataset'
        self.forms_dict[self.current_form] = form_remove_database_dataset

        # raise "form_remove_database_dataset" to front
        form_remove_database_dataset.tkraise()

    #---------------

    def recreate_read_transfer_config_file(self):
        '''
        Recreate the read transfer config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_read_transfer_config_file" in "container" with the grid geometry manager
        form_recreate_read_transfer_config_file = gdataset.FormRecreateReadTransferConfigFile(self)
        form_recreate_read_transfer_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_read_transfer_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_read_transfer_config_file'
        self.forms_dict[self.current_form] = form_recreate_read_transfer_config_file

        # raise "form_recreate_read_transfer_config_file" to front
        form_recreate_read_transfer_config_file.tkraise()

    #---------------

    def edit_read_transfer_config_file(self):
        '''
        Edit the read transfer config file to change the parameters of each transfer.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Read dataset file transfer - Edit config file'

        # edit the read transfer config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xread.get_read_transfer_config_file())
        self.root.wait_window(dialog_editor)

        # check the read transfer config file
        (OK, error_list) = xread.check_read_transfer_config_file(strict=False)
        if OK:
            message = 'The read transfer config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def upload_read_dataset(self):
        '''
       Upload a read dataset to a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_upload_read_dataset" in "container" with the grid geometry manager
        form_upload_read_dataset = gdataset.FormUploadReadDataSet(self)
        form_upload_read_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_upload_read_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_upload_read_dataset'
        self.forms_dict[self.current_form] = form_upload_read_dataset

        # raise "form_upload_read_dataset" to front
        form_upload_read_dataset.tkraise()

    #---------------

    def recreate_read_gzip_config_file(self):
        '''
       Recreate the read file compression/decompression config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_read_gzip_config_file" in "container" with the grid geometry manager
        form_recreate_read_gzip_config_file = gdataset.FormRecreateReadGzipConfigFile(self)
        form_recreate_read_gzip_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_read_gzip_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_read_gzip_config_file'
        self.forms_dict[self.current_form] = form_recreate_read_gzip_config_file

        # raise "form_recreate_read_gzip_config_file" to front
        form_recreate_read_gzip_config_file.tkraise()

    #---------------

    def edit_read_gzip_config_file(self):
        '''
       Edit the read file compression/decompression config file.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Read dataset file compression/decompression - Edit config file'

        # edit the read gzip config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgzip.get_gzip_config_file('read'))
        self.root.wait_window(dialog_editor)

        # check the read gzip config file
        (OK, error_list) = xgzip.check_gzip_config_file('read', strict=False)
        if OK:
            message = 'The read gzip config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_read_gzip_process(self):
        '''
        Compress/decompress read dataset files in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_run_read_gzip_process" in "container" with the grid geometry manager
        form_run_read_gzip_process = gdataset.FormRunGzipProcess(self, 'read')
        form_run_read_gzip_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_read_gzip_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_read_gzip_process'
        self.forms_dict[self.current_form] = form_run_read_gzip_process

        # raise "form_run_read_gzip_process" to front
        form_run_read_gzip_process.tkraise()

    #---------------

    def remove_read_dataset(self):
        '''
       Remove a read dataset in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_read_dataset" in "container" with the grid geometry manager
        form_remove_read_dataset = gdataset.FormRemoveReadDataSet(self)
        form_remove_read_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_read_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_read_dataset'
        self.forms_dict[self.current_form] = form_remove_read_dataset

        # raise "form_remove_read_dataset" to front
        form_remove_read_dataset.tkraise()

    #---------------

    def recreate_result_transfer_config_file(self):
        '''
        Recreate the result transfer config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_result_transfer_config_file" in "container" with the grid geometry manager
        form_recreate_result_transfer_config_file = gdataset.FormRecreateResultTransferConfigFile(self)
        form_recreate_result_transfer_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_result_transfer_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_result_transfer_config_file'
        self.forms_dict[self.current_form] = form_recreate_result_transfer_config_file

        # raise "form_recreate_result_transfer_config_file" to front
        form_recreate_result_transfer_config_file.tkraise()

    #---------------

    def edit_result_transfer_config_file(self):
        '''
        Edit the result transfer config file of a run to change the parameters of each transfer.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Result dataset file transfer - Edit config file'

        # edit the result transfer config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xresult.get_result_transfer_config_file())
        self.root.wait_window(dialog_editor)

        # check the result transfer config file
        (OK, error_list) = xresult.check_result_transfer_config_file(strict=False)
        if OK:
            message = 'The result transfer config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def download_result_dataset(self):
        '''
        Download a result dataset from a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_download_result_dataset" in "container" with the grid geometry manager
        form_download_result_dataset = gdataset.FormDownloadResultDataSet(self)
        form_download_result_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_download_result_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_download_result_dataset'
        self.forms_dict[self.current_form] = form_download_result_dataset

        # raise "form_download_result_dataset" to front
        form_download_result_dataset.tkraise()

    #---------------

    def recreate_result_gzip_config_file(self):
        '''
       Recreate the result file compression/decompression config file.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_recreate_result_gzip_config_file" in "container" with the grid geometry manager
        form_recreate_result_gzip_config_file = gdataset.FormRecreateResultGzipConfigFile(self)
        form_recreate_result_gzip_config_file.grid(row=0, column=0, sticky='nsew')

        # set "form_recreate_result_file_compression_decompression_config_file" as current form and add it in the forms dictionary
        self.current_form = 'form_recreate_result_gzip_config_file'
        self.forms_dict[self.current_form] = form_recreate_result_gzip_config_file

        # raise "form_recreate_result_gzip_config_file" to front
        form_recreate_result_gzip_config_file.tkraise()

    #---------------

    def edit_result_gzip_config_file(self):
        '''
       Edit the result file compression/decompression config file.
        '''

        # initialize the control variable
        OK = True

        # close the current form
        self.close_current_form()

        # set the head
        head = 'Result dataset file compression/decompression - Edit config file'

        # edit the result gzip config file using "DialogEditor" 
        dialog_editor = gdialogs.DialogEditor(self.root, xgzip.get_gzip_config_file('result'))
        self.root.wait_window(dialog_editor)

        # check the result gzip config file
        (OK, error_list) = xgzip.check_gzip_config_file('result', strict=False)
        if OK:
            message = 'The result gzip config file is OK.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {head}', message)
        else:
            message = 'Detected errors:\n\n'
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {head}', message)

    #---------------

    def run_result_gzip_process(self):
        '''
        Compress/decompress result dataset files in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_result_result_gzip_process" in "container" with the grid geometry manager
        form_run_result_gzip_process = gdataset.FormRunGzipProcess(self, 'result')
        form_run_result_gzip_process.grid(row=0, column=0, sticky='nsew')

        # set "form_run_result_gzip_process" as current form and add it in the forms dictionary
        self.current_form = 'form_run_result_gzip_process'
        self.forms_dict[self.current_form] = form_run_result_gzip_process

        # raise "form_run_result_gzip_process" to front
        form_run_result_gzip_process.tkraise()

    #---------------

    def remove_result_dataset(self):
        '''
       Remove a run result dataset in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_result_dataset" in "container" with the grid geometry manager
        form_remove_result_dataset = gdataset.FormRemoveResultDataSet(self)
        form_remove_result_dataset.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_result_dataset" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_result_dataset'
        self.forms_dict[self.current_form] = form_remove_result_dataset

        # raise "form_remove_result_dataset" to front
        form_remove_result_dataset.tkraise()

    #---------------

    def remove_experiment(self):
        '''
       Remove all datasets of an experiment in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_remove_experiment_datasets" in "container" with the grid geometry manager
        form_remove_experiment_datasets = gdataset.FormRemoveExperiment(self)
        form_remove_experiment_datasets.grid(row=0, column=0, sticky='nsew')

        # set "form_remove_experiment_datasets" as current form and add it in the forms dictionary
        self.current_form = 'form_remove_experiment_datasets'
        self.forms_dict[self.current_form] = form_remove_experiment_datasets

        # raise "form_remove_experiment_datasets" to front
        form_remove_experiment_datasets.tkraise()

    #---------------

    def view_cluster_start_log(self):
        '''
        View the cluster start log.
        '''

        # close the current form
        self.close_current_form()

        # create and register "view_cluster_start_log" in container with the grid geometry manager
        form_view_cluster_start_log = glog.FormViewClusterStartLog(self)
        form_view_cluster_start_log.grid(row=0, column=0, sticky='nsew')

        # set "form_view_cluster_start_log" as current form and add it in the forms dictionary
        self.current_form = 'form_view_cluster_start_log'
        self.forms_dict[self.current_form] = form_view_cluster_start_log

        # raise "form_view_cluster_start_log" to front
        form_view_cluster_start_log.tkraise()

    #---------------

    def view_submission_logs(self):
        '''
        List logs of process submission in local computer.
        '''

        # close the current form
        self.close_current_form()

        # create and register "form_view_submission_logs" in container with the grid geometry manager
        form_view_submission_logs = glog.FormViewSubmissionLogs(self)
        form_view_submission_logs.grid(row=0, column=0, sticky='nsew')

        # set "form_view_submission_logs" as current form and add it in the forms dictionary
        self.current_form = 'form_view_submission_logs'
        self.forms_dict[self.current_form] = form_view_submission_logs

        # raise "form_view_submission_logs" to front
        form_view_submission_logs.tkraise()

    #---------------

    def view_result_logs(self):
        '''
        List logs of results in a cluster.
        '''

        # close the current form
        self.close_current_form()

        # create and register "view_result_logs" in container with the grid geometry manager
        form_view_result_logs = glog.FormViewResultLogs(self)
        form_view_result_logs.grid(row=0, column=0, sticky='nsew')

        # set "form_view_result_logs" as current form and add it in the forms dictionary
        self.current_form = 'form_view_result_logs'
        self.forms_dict[self.current_form] = form_view_result_logs

        # raise "form_view_result_logs" to front
        form_view_result_logs.tkraise()

    #---------------

    def open_help(self, event=None):
        '''
        Open the help file.
        '''

        try:
            manual_file = os.path.abspath(xlib.get_project_manual_file())
            webbrowser.open_new(f'file://{manual_file}')
        except:
            message = f'The document {manual_file}\n is not available.'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - Open help', message)

    #---------------

    def show_dialog_about(self):
        '''
        Show the application information.
        '''

        dialog_about = gdialogs.DialogAbout(self.root)
        self.root.wait_window(dialog_about)

    #---------------

    def warn_unavailable_process(self):

        message = 'This process is been built.\nIt is coming soon!'
        tkinter.messagebox.showwarning(xlib.get_project_name(), message)

    #---------------

    def close_current_form(self):
        '''
        Close the current form.
        '''

        # clear the label of the current process name
        self.label_process['text'] = ''

        # destroy the current form
        if self.current_form != 'form_welcome':
            self.forms_dict[self.current_form].destroy()
            self.forms_dict['form_welcome'].tkraise()

    #---------------

    def exit(self, event=None):
        '''
        Exit the application.
        '''

        message = f'Are you sure to exit {xlib.get_project_name()}?'
        if tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - Exit', message):
            self.close_current_form()
            self.root.quit()
            self.root.destroy()
            exit()

   #---------------

#-------------------------------------------------------------------------------

class FormWelcome(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormWelcome" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # build the graphical user interface
        self.build_gui()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormWelcome".
        '''

        # create "image_hierro"
        image_hierro = PIL.Image.open('./image_hierro.jpg')
        image_hierro.thumbnail((self.main.WINDOW_WIDTH,self.main.WINDOW_HEIGHT), PIL.Image.ANTIALIAS)

        # create "photoimage_hierro"
        self.photoimage_hierro = PIL.ImageTk.PhotoImage(image_hierro)  

        # create "canvas_photoimage_hierro" and register it with the grid geometry manager
        self.canvas_photoimage_hierro = tkinter.Canvas(self, width=self.main.WINDOW_WIDTH, height=self.main.WINDOW_HEIGHT)
        self.canvas_photoimage_hierro.create_image(round(self.main.WINDOW_WIDTH / 2), round(self.main.WINDOW_HEIGHT / 2), image=self.photoimage_hierro, anchor='center')
        if sys.platform.startswith('linux'):
            x_coordinate = 10
            y_coordinate = self.main.WINDOW_HEIGHT - 100
        elif sys.platform.startswith('darwin'):
            x_coordinate = 10
            y_coordinate = self.main.WINDOW_HEIGHT - 85
        elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
            x_coordinate = 10
            y_coordinate = self.main.WINDOW_HEIGHT - 70
        self.canvas_photoimage_hierro.create_text(x_coordinate, y_coordinate, anchor='w', text='Canary Island pines over the sea (El Hierro, Canary Islands, Spain)', fill='white') 
        self.canvas_photoimage_hierro.pack(side='left', fill='both', expand=True)

    #---------------

    def close(self):
        '''
        Close "FormWelcome".
        '''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print(f'This file contains the class Main corresponding to the graphical user interface of the {xlib.get_project_name()} software package.')
    sys.exit(0)

#-------------------------------------------------------------------------------
